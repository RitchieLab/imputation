#!/bin/bash

# combines functionality of the strand and phasing apps from DNAnexus

scriptname="strand_phase.sh"

usage(){
  echo ""
  echo "\t$scriptname uses SHAPEIT to run strand alignment followed by phasing.  The phasing step submits an array job to the cluster (one job per chromosome)"
  echo ""
  echo "	Usage:  $scriptname: -b binary_file [-s strand_file] [-p remove_palindromic] [-r reference_directory] [-c chromosome] [-n cores] [-t walltime]"
  echo "		-b:  Prefix for bed, bim, fam file (plink)"
  echo " 		-s:  Strand alignment file (optional)" 
  echo "		-p:  Remove palindromic SNPs"
  echo "		-r:  Reference file location (default=/project/ritchie00/datasets/1KG_Phase3/Impute2_files/1000GP_Phase3)"
  echo "		-c:  Chromosome to process (optional -- default is to process all chromosomes)"
  echo "		-n:  Cores to use for phasing (default=4)"
  echo "		-t:  Wall time for bsub array phasing job (default=\"12:00\")"
  echo 
  exit 0
}


prefix=""
chromosome=""
strand_file=""
remove_palindromic=false
ref_dir="/project/ritchie00/datasets/1KG_Phase3/Impute2_files/1000GP_Phase3"
ncores=4
walltime="12:00"

while getopts ":b:s:r:c:n:t:p" opt; do
  case ${opt} in
	b ) prefix=$OPTARG
     ;;
    s ) strand_file=$OPTARG
     ;;
    p ) remove_palindromic=true
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

if [ -z "$prefix" ]
then
	usage
fi



# create and submit an alignment/phasing job for a chromosome or chromosomes
run_alignment_phasing(){
	chrfile=$1
	fileprefix=$2
	thousandkdir=$3
	excl_pal=$4
	bsubfile=${fileprefix}_strand_phase.bsub

	chroms=()
	for i in $(cat $chr_file)
	do
		chroms+=($i)
	done
	
	chromlength=${#chroms[@]}
	if (( chromlength > 1 )); then
		arraynumbers="1-$chromlength"
	else
		arraynumbers="1"
	fi	
	
	# create and write bsub script for running alignment and phasing per chromosome
	
	cat > $bsubfile << EOL
#!/bin/bash
#BSUB -J "${fileprefix}_strand[$arraynumbers]"
#BSUB -o bsub.%J-%I.out
#BSUB -R "rusage[mem=8192]"
#BSUB -M 8192
#BSUB -W $walltime
#BSUB -R "span[ptile=$ncores]"
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
#module load shapeit
EOL
	
	chromline="CHROMS=("
	for i in $"${chroms[@]}"
	do
		chromline="$chromline $i"
	done
	chromline="$chromline )"
	echo $chromline >> $bsubfile

	
cat >> $bsubfile << EOL
arrayindex=\$((\$LSB_JOBINDEX-1))
chrnum=\${CHROMS[\$arrayindex]}
	
sample_file=$thousandkdir/1000GP_Phase3.sample
haplotype_file=$thousandkdir/1000GP_Phase3_chr\${chrnum}.hap.gz
legend_file=$thousandkdir/1000GP_Phase3_chr\${chrnum}.legend.gz
PLINK_ARGS=""
EOL

if test "$excl_pal" == "true"; then
cat >> $bsubfile << EOL	
## Split PLINK files
PLINK_ARGS="--exclude  ${fileprefix}_palindromic_SNPs.txt"		
EOL
fi

cat >> $bsubfile << EOL
# create plink files using chromosome selected and excluding palindromic (if indicated)
plink --allow-no-sex --bfile $fileprefix --chr \${chrnum} \$PLINK_ARGS --make-bed --out ${fileprefix}_chr\${chrnum} --noweb

## Run ShapeIt to check strands
shapeit.v2.r644.linux.x84_64 -check -B ${fileprefix}_chr\${chrnum} --input-ref \$haplotype_file \$legend_file \$sample_file --output-log ${fileprefix}_chr\${chrnum}.alignments -T $ncores || true

# check for error in log file and abort if found
#if test -f ${fileprefix}_chr\${chrnum}.alignments.log -a \$(grep "ERROR" ${fileprefix}_chr\${chrnum}.alignments.log | wc -l) -gt 0; then
#	cat ${fileprefix}_chr\${chrnum}.alignments
#	exit 1
#fi

# flip strands and then recheck
if test -f ${fileprefix}_chr\${chrnum}.alignments.snp.strand -a \$(grep "^strand" ${fileprefix}_chr\${chrnum}.alignments.snp.strand | wc -l) -gt 0; then
	## Store all Strand Alignment problems in a separate file
		
	echo "Log File:"
	cat ${fileprefix}_chr\${chrnum}.alignments.log
		
	echo "Strand problems"
	grep '^Strand' ${fileprefix}_chr\${chrnum}.alignments.snp.strand || true

	# add monomorphic alleles and flip simple strands
	fix_monomorphic.py ${fileprefix}_chr\${chrnum}.alignments.snp.strand > ${fileprefix}_allele_recode_chr\${chrnum}
		
	echo "Allele fixer:"
	cat ${fileprefix}_allele_recode_chr\${chrnum}
		
	## Use PLINK to flip strands for list of snps 
	plink --allow-no-sex --bfile ${fileprefix}_chr\${chrnum} --update-alleles ${fileprefix}_allele_recode_chr\${chrnum} --make-bed --out ${fileprefix}_chr\${chrnum}_flipped
		
	# recheck the strands using shapeit and report any errors!
	shapeit.v2.r644.linux.x84_64 -check -B ${fileprefix}_chr\${chrnum}_flipped --input-ref \$haplotype_file \$legend_file \$sample_file --output-log ${fileprefix}_chr\${chrnum}.alignments_recheck -T $ncores || true
		
	if test -f ${fileprefix}_chr\${chrnum}.alignments_recheck.snp.strand -a \$(grep "^strand" ${fileprefix}_chr\${chrnum}.alignments_recheck.snp.strand | wc -l) -gt 0; then
		echo "WARNING: SNPs that still have strand issues:"
		grep "^strand" ${fileprefix}_chr\${chrnum}.alignments_recheck.snp.strand
	fi
	
	# remove previous non-flipped plink files	
	rm ${fileprefix}_chr\${chrnum}.*
else
	# rename the files so they would be the same as if flipped
	mv ${fileprefix}_chr\${chrnum}.bed ${fileprefix}_chr\${chrnum}_flipped.bed
	mv ${fileprefix}_chr\${chrnum}.bim ${fileprefix}_chr\${chrnum}_flipped.bim
	mv ${fileprefix}_chr\${chrnum}.fam ${fileprefix}_chr\${chrnum}_flipped.fam
fi	
EOL


# continue on with the phasing app work on each bed/bim/fam chromosome file set
cat >> $bsubfile << EOL
genetic_map=$thousandkdir/genetic_map_chr\${chrnum}_combined_b37.txt
# set 5% missingness maximum so we don't get complaints from shapeit
plink --allow-no-sex --noweb --geno 0.05 --bfile ${fileprefix}_chr\${chrnum}_flipped --chr \${chrnum} --make-bed --out ${fileprefix}_chr\${chrnum}_5percent

# remove the flipped version of the plink files as not needed
rm ${fileprefix}_chr\${chrnum}_flipped.*

# run phase on the 5 percent threshold plink files
shapeit.v2.r644.linux.x84_64 -B  ${fileprefix}_chr\${chrnum}_5percent -M \$genetic_map -O ${fileprefix}_chr\${chrnum}.phased -T $ncores 

# remove the 5percent plink files
rm ${fileprefix}_chr\${chrnum}_5percent.*

EOL

echo "Created $bsubfile"
# submit array job to queue
# bsub < $bsubfile

}


# if given the strand file, process it and create replacement bed/bam/fam
if [ ! -z "$strand_file" ] ; then
	CHR_FILE=$(mktemp)
	POS_FILE=$(mktemp)
	FLIP_FILE=$(mktemp)
	
	cut -f 1,2 "$strand_file" > $CHR_FILE
	cut -f 1,3 "$strand_file" > $POS_FILE	
	awk '{if ($5=="-") print $0}' "$strand_file" | cut -f 1 > $FLIP_FILE

	#1. Apply the chr
	plink --noweb --allow-no-sex --bfile $prefix --update-map $CHR_FILE --update-chr --make-bed --out ${prefix}_data_chr
	#2. Apply the pos
	plink --noweb --allow-no-sex --bfile ${prefix}_data_chr --update-map $POS_FILE --make-bed --out ${prefix}_data_pos
	rm ${prefix}_data_chr.*
	#3. Apply the flip
	plink --noweb --allow-no-sex --bfile ${prefix}_data_pos --flip $FLIP_FILE --make-bed --out ${prefix}_data_flip
	rm ${prefix}_data_pos.*
	#4. Extract the SNPs in the pos file, we don't want SNPs that aren't in the strand file
	plink --noweb --allow-no-sex --bfile ${prefix}_data_flip --extract $POS_FILE --make-bed --out ${prefix}_final	
	rm  ${prefix}_data_flip.*
	
	rm $CHR_FILE
	rm $POS_FILE
	rm $FLIP_FILE	
	prefix=${prefix}_final
	
fi

if test "$remove_palindromic" == "true"; then
	# checking for palindromic SNPs and then adding option to exclude them
	cat ${prefix}.bim | awk '{if($5=="A" && $6=="T" || $5=="T" && $6=="A" || $5=="G" && $6=="C" || $5=="C" && $6=="G")print $2}' > ${prefix}_palindromic_SNPs.txt
fi

# create chromosome list from bim file or use one passed in as command-line argumebnt
chr_file=$(mktemp)
if [ -z $chromosome ] ; then	
	sed 's/ /\t/g' ${prefix}.bim | cut -f1 | sort -u | awk '$1>=1 && $1<=22' > $chr_file
else
	echo $chromosome > $chr_file
fi

# echo "21" > $chr_file
# echo "22" >> $chr_file

echo "Chromosomes included:"
cat $chr_file

run_alignment_phasing $chr_file $prefix $ref_dir $remove_palindromic









