#include <stdio.h>
#include <time.h>
#include <CUnit.h>
#include <genome_complexity.h>
#include <limits.h>

void test_count_reads_where_snp_makes_clean_bubble1()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    printf("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  //  simple test. Use this refernece genome which I load into colour 0:

  /// file: data/test/genome_complexity/test_allele_clean_file1.fa

  //  >ref contains CAAGTTC CACGTTC and CAGGTTC
  //  CAAGTTCAGAGTTACTCACACCCGATCGATAAGCGGTACAGAGCACGTTCAGAAAAAAACAGGTTCAGA
  //  >ref + SNP
  //  TATCCATGTTCAGAGTTACTGACACCCGATCGATAAGCG


  //first set up the hash/graph
  int kmer_size = 7;
  int number_of_bits = 8;
  int bucket_size = 10;
  int max_retries = 10;

  dBGraph *hash_table = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);

  if (hash_table==NULL)
    {
      die("Unable to alloc the hash table.");
    }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/genome_complexity/pop_first_test.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0);

  //and use this file of reads: /data/test/genome_complexity/test_allele_clean_file2.fa
  //  >read lies entirely in graph defined by test_allele_clean_file1.fa, and is clean (forms supernode at k=7)
  //  GTTCAGAGTTACT
  //  >read lies in graph defined by test_allele_clean_file1.fa, but overlaps a junction so is not clean
  //  TTACTGACACCCGATCG
  //  >read lies in ref but any change to the 7th base results in an overlap with the ref, so
  //  GTTCAGAGTTACTCACA

  int col_genome=0;
  int reads_tested=0;
  int reads_where_snp_makes_clean_bubble = 0;
  dBNode* array_nodes[500];
  Orientation array_or[500];
  
    //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence\n");
  }
  int max_read_length=2000;
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      die("Failed to malloc kmer sliding window in test. Exit.\n");
    }
  

  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-hash_table->kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      die("Failed to malloc kmer_window->kmer in test. Exit.\n");
    }
  kmer_window->nkmers=0;
  
  
  //end of intialisation 
	  
	  
  //create file readers
  int file_reader_fasta(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = hash_table->kmer_size;
      //die("new_entry must be true in hsi test function");
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }

  count_reads_where_snp_makes_clean_bubble(hash_table, "../data/test/genome_complexity/test_allele_clean_file2.fa", true,
					   col_genome, &reads_tested, &reads_where_snp_makes_clean_bubble, 
					   file_reader_fasta,
					   array_nodes, array_or, seq, kmer_window, max_read_length);


  
  CU_ASSERT(reads_tested==3);
  CU_ASSERT(reads_where_snp_makes_clean_bubble==2);

    
}
