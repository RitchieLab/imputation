#!/bin/bash

scriptname="impute2.sh"

usage(){
  echo ""
  echo "$scriptname utilizes the output from SHAPEIT to run impute2.  It submits an array job to the cluster."
  echo ""
  echo "	Usage:  $scriptname: -h hap_files_prefix [-h hap_files_directory] [-n cores] [-r reference_directory] [-c chromosomes] [-t walltime]"
  echo "		-p:  Prefix for hap files"
  echo "		-h:  Directory for hap files (default is current directory)"
  echo " 		-s:  SNPs per chunk (default=56000)" 
  echo "		-r:  Reference file location (default=/project/ritchie00/datasets/1KG_Phase3/Impute2_files/1000GP_Phase3)"
  echo "		-c:  Chromosome to process (optional -- default is to process all chromosomes)"
  echo "		-n:  Number of cores to use for each chromosome (default=8)"
  echo "		-t:  Wall time for bsub array job (default=\"12:00\")"
  echo 
  exit 0
}

happrefix=""
hap_dir="."
chromosome=""
ncores=8
ref_dir="/project/ritchie00/datasets/1KG_Phase3/Impute2_files/1000GP_Phase3"
walltime="12:00"
snps_per_chunk=56000

while getopts ":p:s:c:r:h:n:t:" opt; do
  case ${opt} in
	p ) happrefix=$OPTARG
     ;;
    s ) snps_per_chunk=$OPTARG
     ;;
    h ) hap_dir=$OPTARG
     ;;
    r ) ref_dir=$OPTARG
     ;;
    c ) chromosome=$OPTARG
     ;;
    n ) ncores=$OPTARG
     ;;
    t ) walltime=$OPTARG
     ;;
    \? ) usage
     ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
     ;;
  esac
done

if [ -z "$happrefix" ] 
then
	usage
fi

# when all chromosomes being processed use the matching hap files from the phasing run
# in the hap file directory
if [ -z $chromsome] ; then
# get all hap files in directory matching prefix
hapfiles=( $(find $hap_dir -type f -name "$happrefix*.haps"))
# get matching chromosome for each hapfile
	chroms=()
	regex="chr([0-9]+)"
	for h in "${hapfiles[@]}"
	do
		if [[ $h =~ $regex ]]
		then
			chroms+=("${BASH_REMATCH[1]}")
		fi
	done
else
	chroms=( $chromosome )
	hapfiles=( $(find $hap_dir -type -f -name "$happrefix*chr${chromosome}*.hap"))
fi

# create bsub array job for imputation
chromlength=${#chroms[@]}

bsubfile=${happrefix}_impute2.bsub

if (( chromlength > 1 )); then
	arraynumbers="1-$chromlength"
else
	arraynumbers="1"
fi

	cat > $bsubfile << EOL
#!/bin/bash
#BSUB -J "${happrefix}_impute2[$arraynumbers]"
#BSUB -o bsub.%J-%I.out
#BSUB -R "rusage[mem=16384]"
#BSUB -M 16384
#BSUB -W $walltime
#BSUB -n $ncores
#BSUB -N
#BSUB -B
#BSUB -R "select[ostype>=CENT7]"
	
if test "\${HOME}/ritchielab.bashrc" -nt "\${HOME}/group/ritchielab.bashrc" ; then
        . "\${HOME}/ritchielab.bashrc"
elif test -f "\${HOME}/group/ritchielab.bashrc" ; then
        . "\${HOME}/group/ritchielab.bashrc"
else
        echo "WARNING: Could not find Ritchie Lab bashrc group environment script."
fi
module load impute pigz

EOL


chromline="CHROMS=("
for i in "${chroms[@]}"
do
	chromline="$chromline $i"
done
chromline="$chromline )"
echo $chromline >> $bsubfile

hapsline="HAPS=("
for i in "${hapfiles[@]}"
do
	hapsline="$hapsline $i"
done
hapsline="$hapsline )"
echo $hapsline >> $bsubfile

cat >> $bsubfile << EOL

arrayindex=\$((\$LSB_JOBINDEX-1))
chromnum=\${CHROMS[\$arrayindex]}
sample_haps=\${HAPS[\$arrayindex]}

spc=$snps_per_chunk
sample_file=$ref_dir/1000GP_Phase3.sample
haplotype_file=$ref_dir/1000GP_Phase3_chr\${chromnum}.hap.gz
legend_file=$ref_dir/1000GP_Phase3_chr\${chromnum}.legend.gz
genetic_map=$ref_dir/genetic_map_chr\${chromnum}_combined_b37.txt

# First, how many chunks do I need to get ~56K reference variants per chunk
# (this is the average # of variants per 5Mb)
regex=".gz\$"
if [[ \$legend_file =~ \$regex ]] ; then
	N_VARS=\$(zcat \$legend_file | wc -l)
else
	N_VARS=\$(cat \$legend_file | wc -l)
fi
# number of chunks is number of vars / snps_per_chunk (+1 for remainder)
N_CHUNKS=\$((N_VARS / spc + 1))

# OK, now how many study variants per chunk and how many left over?
# VPC is variants per chunk
# LO are number left over 
N_SVARS=\$(cat \$sample_haps | wc -l)
VPC=\$((N_SVARS / N_CHUNKS))
LO=\$((N_SVARS % N_CHUNKS))

# Now, get the break points by splitting the study variants into the # per
# chunk above and get the midpoint between them
BREAK_FN=\$(mktemp)
cut -d' ' -f3 \$sample_haps | awk "NR>\$VPC-1 && NR < \$VPC*\$N_CHUNKS && NR % (\$VPC + (NR/\$VPC < (\$LO + 1))) <= 1" | awk "{C=\\\$1; if(NR % 2 == 1){X=C} else {print int((C+X)/2)} }" > \$BREAK_FN

# get the prefix of the haplotype file
#prefix="\$(echo \$sample_haps | sed 's/\..*//').\$chromnum"
prefix=$happrefix


output_prefix=\${prefix}_chr\$chromnum

CHUNK_FILE=\${output_prefix}_chunk_file

#Now, generate the first chunks by iterating over the break file
if [ -f \$CHUNK_FILE ]; then
	rm \$CHUNK_FILE
fi
touch \$CHUNK_FILE
LAST_BREAK=0
for b in \$(cat \$BREAK_FN); do
		echo "\$((LAST_BREAK + 1)) \$b" >> \$CHUNK_FILE
		LAST_BREAK=\$b
done

# get the last position by looking at the tail of the ref_legend file
regex=".gz\$"
if [[ \$legend_file =~ \$regex ]] ; then
LAST_POS=\$(zcat \$legend_file | tail -n 1 | cut -d' ' -f2)
else
LAST_POS=\$(tail -n 1 \$legend_file | cut -d' ' -f2)
fi
echo "\$((LAST_BREAK + 1)) \$((LAST_POS + 1))" >> \$CHUNK_FILE

echo "Contents of chunk file:"
cat \$CHUNK_FILE



CHUNK_DIR="\${output_prefix}_chunk"
if [ ! -d \$CHUNK_DIR ] ; then
	mkdir \$CHUNK_DIR
fi

CONCAT_DIR="\${output_prefix}_concat"
if [ ! -d \$CONCAT_DIR ] ; then
	mkdir \$CONCAT_DIR
fi



# run impute2
mpirun --mca mpi_cuda_support 0 -np $ncores mpi_impute \$CHUNK_FILE \$genetic_map \$haplotype_file \$legend_file \$sample_haps \$CHUNK_DIR/\$output_prefix

# merge all of the files in the CHUNK_DIR, but sort on the 3rd column, please
# NOTE: we are only sorting files based on the position of the 1st line
POS_F=\$(mktemp)

for f in \${CHUNK_DIR}/\$output_prefix.*.impute2; do
	# Read the 1st 1,000 bytes and print the 3rd field.  Should eliminate 
	# issues with too many fields for awk to handle
	echo -e "\$(head -1 \$f | head -c 1000 | awk '{print \$3}')\t\$f"
done | tee \$POS_F


POS_SORT_F=\$(mktemp)
sort -t\$'\t' -gk1,1 \$POS_F > \$POS_SORT_F

while read f; do
	cat "\$f"
done < <(cut -f2- \$POS_SORT_F) | pigz > \$CONCAT_DIR/\$output_prefix.impute2.gz

#OK, make sure this same logic is implemented for the impute2_info files
# we MUST match the lines exactly!
NLINE=1
while read f; do
	# This should print the header line in the first 
	tail -n+\${NLINE} "\$(echo "\$f" | sed 's/\$/_info/')"
	NLINE=2
done < <(cut -f2- \$POS_SORT_F) | pigz > \$CONCAT_DIR/\$output_prefix.impute2_info.gz

while read f; do
	# This should print the header line in the first 
	cat "\$(echo "\$f" | sed 's/\$/_summary/')"
done < <(cut -f2- \$POS_SORT_F) > \$CONCAT_DIR/\$output_prefix.summary

echo "Warnings:"
cat \${CHUNK_DIR}/\$output_prefix.*.impute2_warnings

rm -r \$CHUNK_DIR

EOL

echo "Created bsub file: $bsubfile"

# submit array job to queue
# bsub < $bsubfile
