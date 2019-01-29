// mpi_impute

#include <stdlib.h>  
#include <mpi.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

int main(int argc, char* argv[]){


	string program_name="mpi_impute";

	if(argc != 7){
		cout << "\nUsage:" << endl;
		cout << program_name << " <chunks file> <map file> <ref haplotypes file> <ref legend file> <phase samples haps file> <prefix for output files>\n\n";
		return 0;
	}	

	int myrank;
	int nproc;
	
	string chunkfile = argv[1];
	string refmapfile = argv[2];
	string refhapfile = argv[3];
	string reflegendfile = argv[4];
	string samplehapfile = argv[5];
	string outputprefix = argv[6]; 

	// set up MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

if(myrank==0){
  cout << "myrank=" << myrank << " nproc=" << nproc << endl;
}

	
	vector<string> chunks;
	// read in chunkfile
	ifstream infile;
	infile.open(chunkfile.c_str());
	if(!infile.is_open()){
		cerr  << "ERROR: Unable to open to open file " << chunkfile << endl << endl;
		MPI_Finalize();
		exit(1);
	}
	
	string line;
	while(getline(infile, line)){
		chunks.push_back(line);
	}
	
	infile.close();
	
	size_t chunk_index=myrank;
	size_t n_chunks = chunks.size();
	
	string NE="20000";
	
	
	string chunk_start, chunk_end;
	string base_command ="impute2 -m " + refmapfile + " -use_prephased_g -known_haps_g " +
		samplehapfile + " -h " + refhapfile + " -l " + reflegendfile + "  -Ne " + NE + " -int ";
	string end_command= " -buffer 250kb -call_thresh 0.9 -allow_large_regions -o ";
	string command, output_file;
	
	// location of impute 2
	// /appl/impute-2.3.2/impute2
	// each process works through its indicated chunks and uses the chunks assigned to it
	while(chunk_index < n_chunks){ 
		// run impute2
		stringstream ss(chunks[chunk_index]);
		ss >> chunk_start >> chunk_end;
		output_file = outputprefix + ".pos" + chunk_start + "-" + chunk_end + 
			".best_guess_haps_imputation.impute2";
		command = base_command + chunk_start + " " + chunk_end + end_command + output_file;
		
		system(command.c_str());
		chunk_index+=nproc;
	}
	
	
	MPI_Finalize();
	return 0;
}

