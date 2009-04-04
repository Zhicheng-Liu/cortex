#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <dB_graph.h>



int MAX_READ_LENGTH;


//these routines return the length of the read sequence, for the binary file is all the kmers conctenated

//for fasta
int load_fasta_data_from_filename_into_graph(char* filename, long long * bad_reads, int max_read_length, dBGraph* db_graph);

//for fastq
int load_fastq_data_from_filename_into_graph(char* filename, long long * bad_reads,  char quality_cut_off, int max_read_length, dBGraph* db_graph);

//for binary
int load_binary_data_from_filename_into_graph(char* filename, dBGraph* db_graph);

int load_binary_data_from_filename_into_graph_with_unique_kmers(char* filename, dBGraph* db_graph);

#endif /* FILE_READER_H_ */
