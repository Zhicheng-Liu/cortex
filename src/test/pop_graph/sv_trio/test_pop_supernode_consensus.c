#include <CUnit.h>
#include <Basic.h>
#include <file_reader.h>
#include <dB_graph_population.h>
#include <dB_graph_supernode.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <stdlib.h>
#include "supernode_cmp.h"


void test_find_first_node_in_supernode()
{

  printf("\nplease reimplement (not urgent as no current users of this code) conssensus code using db_graph_supernode NULL TEST\n. ");


  //*********************************************************************************
  // 1. Test two simple fasta, one each for a different person in the same pop
  //
  //       >person1_read1
  //        AAAAAAAAAAAAAAA
  //        >person1_read2
  //        ACGTT
  //
  //
  //       >person2_read1
  //       GGGGGGGGGGGGGGGGGG
  //       >person2_read2
  //       TTGACG


  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size    = 4;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  char tmp_seq[hash_table->kmer_size];

  
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }


  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/test_pop_load_and_print/two_individuals_simple.txt", &bad_reads, hash_table);
  //printf("Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 44);
  CU_ASSERT(bad_reads==0);

  // Now just see if it can correctly find the first node in a supernode


  // *******
  //first check all the kmers in person 2. all of them should only be found in person2, except ACG
  // ****


  // **** GGG first

  dBNode* query_node = hash_table_find(element_get_key(seq_to_binary_kmer("GGG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT((query_node!=NULL));
  dBNode* testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  
  //this kmer doesn't exist in person 1.
  CU_ASSERT(testnode==NULL);

  //try again for person 2:
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(testnode != NULL);  
  binary_kmer_to_seq(testnode ->kmer, hash_table->kmer_size, tmp_seq);
  CU_ASSERT( !strcmp(tmp_seq, "GGG") || !strcmp(tmp_seq, "CCC") );
  

  // *** then TTG
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("TTG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT((query_node!=NULL));
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  
  //this kmer doesn't exist in person 1.
  CU_ASSERT(testnode==NULL);

  //try again for person 2:
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(testnode != NULL);  
  binary_kmer_to_seq(testnode ->kmer, hash_table->kmer_size, tmp_seq);
  CU_ASSERT( !strcmp(tmp_seq, "TTG") || !strcmp(tmp_seq, "CAA") || !strcmp(tmp_seq, "CGT") || !strcmp(tmp_seq, "ACG") );
  
  

  // *****then TGA
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("TGA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT((query_node!=NULL));
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  
  //this kmer doesn't exist in person 1.
  CU_ASSERT(testnode==NULL);

  //try again for person 2:
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(testnode != NULL);  
  binary_kmer_to_seq(testnode ->kmer, hash_table->kmer_size, tmp_seq);
  CU_ASSERT( !strcmp(tmp_seq, "TTG") || !strcmp(tmp_seq, "CAA") || !strcmp(tmp_seq, "CGT") || !strcmp(tmp_seq, "ACG") );
  


  //*** then GAC
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("GAC",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT((query_node!=NULL));
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  
  //this kmer doesn't exist in person 1.
  CU_ASSERT(testnode==NULL);

  //try again for person 2:
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(testnode != NULL);  
  binary_kmer_to_seq(testnode ->kmer, hash_table->kmer_size, tmp_seq);
  CU_ASSERT( !strcmp(tmp_seq, "TTG") || !strcmp(tmp_seq, "CAA") || !strcmp(tmp_seq, "CGT") || !strcmp(tmp_seq, "ACG") );
  

  
  // then ACG - this is the one that I expect to see in both person1 and person2

  //check person1
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ACG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  CU_ASSERT(!(testnode==NULL));
  binary_kmer_to_seq(testnode->kmer, hash_table->kmer_size, tmp_seq);
  CU_ASSERT( !strcmp(tmp_seq, "ACG") || !strcmp(tmp_seq, "AAC")  || !strcmp(tmp_seq, "CGT") || !strcmp(tmp_seq, "GTT") );
  


  //then person2
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(!(testnode==NULL));
  binary_kmer_to_seq(testnode->kmer, hash_table->kmer_size, tmp_seq);
  CU_ASSERT( !strcmp(tmp_seq, "ACG") || !strcmp(tmp_seq, "CGT")  || !strcmp(tmp_seq, "TTG") || !strcmp(tmp_seq, "CAA") );
  


  //*****
  // then do the kmers you expect to find only in person1 (not person 2) (horrible having index from 0 to 1 and people labelled 1 and 2 - stupid of me)


  // *** AAA
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("AAA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  CU_ASSERT(!(testnode==NULL));
  binary_kmer_to_seq(testnode->kmer, hash_table->kmer_size, tmp_seq);
  CU_ASSERT( !strcmp(tmp_seq, "AAA") || !strcmp(tmp_seq, "TTT") );
  

  //confirm not seen in person1
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(testnode==NULL);



  // *** GTT

  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("GTT",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  CU_ASSERT(!(testnode==NULL));
  binary_kmer_to_seq(testnode->kmer, hash_table->kmer_size, tmp_seq);
  CU_ASSERT( !strcmp(tmp_seq, "ACG") || !strcmp(tmp_seq, "CGT")  ||  !strcmp(tmp_seq, "GTT") || !strcmp(tmp_seq, "AAC") );
  

  //confirm not seen in person1
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(testnode==NULL);
    
  
  hash_table_free(&hash_table);

}



void test_find_next_node_in_supernode()
{

  printf("\nplease reimplement (not urgent as no current users of this code) conssensus code using db_graph_supernode NULL TEST\n. ");

  //*********************************************************************************
  // 1. Test two simple fasta, one each for a different person in the same pop
  //
  //       >person1_read1
  //        AAAAAAAAAAAAAAA
  //        >person1_read2
  //        ACGTT
  //
  //
  //       >person2_read1
  //       GGGGGGGGGGGGGGGGGG
  //       >person2_read2
  //       TTGACG


  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size    = 4;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);


  char tmp_seq[hash_table->kmer_size];

  
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }


  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/test_pop_load_and_print/two_individuals_simple.txt", &bad_reads, hash_table);
  //printf("Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 44);

  // Now just see if it can correctly find the first node in a supernode, and then walk along to the end. Work with person 2 for the moment (=index 1 in the array)

  //first get the first node in the supernode
  dBNode* query_node = hash_table_find(element_get_key(seq_to_binary_kmer("GAC",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT((query_node!=NULL));
  dBNode* testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(testnode != NULL);  
  binary_kmer_to_seq(testnode ->kmer, hash_table->kmer_size, tmp_seq);
  CU_ASSERT( !strcmp(tmp_seq, "TTG") || !strcmp(tmp_seq, "CGT")  || !strcmp(tmp_seq, "CAA") || !strcmp(tmp_seq, "ACG") );//zax
  


  //now start walking 
  Orientation start_orientation, next_orientation, current_orientation;
  dBNode* current_node=testnode;

  if (db_node_is_supernode_end(testnode, forward, individual_edge_array, 1, hash_table))
    {
      start_orientation=reverse;
    }
  else if (db_node_is_supernode_end(testnode, reverse, individual_edge_array, 1, hash_table))
    {
      start_orientation=forward;
    }

  dBNode* next_node = db_graph_get_next_node_in_supernode_for_specific_person_or_pop(testnode, start_orientation, &next_orientation, individual_edge_array, 1, hash_table);
  char* next_kmer= binary_kmer_to_seq(next_node->kmer, hash_table->kmer_size, tmp_seq);
  CU_ASSERT( !strcmp(next_kmer,"TGA") || !strcmp(next_kmer,"GTC") || !strcmp(next_kmer,"TCA") || !strcmp(next_kmer,"GAC")); 
  


  current_node=next_node;
  current_orientation=next_orientation;

  next_node = db_graph_get_next_node_in_supernode_for_specific_person_or_pop(testnode, start_orientation, &next_orientation, individual_edge_array, 1, hash_table);
  next_kmer= binary_kmer_to_seq(next_node->kmer, hash_table->kmer_size, tmp_seq);
  CU_ASSERT( !strcmp(next_kmer,"GAC") || !strcmp(next_kmer,"TCA") || !strcmp(next_kmer,"GTC") || !strcmp(next_kmer,"TGA")); 


  //*********************
  //OK - that seems ok. Let's just check that it copes with a self-looping kmer
  //*********************
  
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("GGG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT((query_node!=NULL));
  
  //first - does it correctly get itself as the first node of the supernode. It should be end of supernode in both directions

  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(testnode != NULL);  
  binary_kmer_to_seq(testnode ->kmer, hash_table->kmer_size, tmp_seq);
  CU_ASSERT(  !strcmp(tmp_seq, "GGG") || !strcmp(tmp_seq, "CCC") );
  

  //now - can it get the next:

  next_node = db_graph_get_next_node_in_supernode_for_specific_person_or_pop(testnode, forward, &next_orientation, individual_edge_array, 1, hash_table);
  CU_ASSERT(next_node == NULL); //because next node is itself

  // OK - now happy with self-loops.
  // what about confusions with orientatio and kmers?



  hash_table_free(&hash_table);

}

void test_correctly_find_subsection_of_supernode()
{



  printf("\nplease reimplement (not urgent as no current users of this code) conssensus code using db_graph_supernode NULL TEST\n. ");


  //first set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 5;
  int bucket_size    = 7;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }

  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/two_people_test_consensus.txt", &bad_reads, hash_table);
  //printf("Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 23);

  
  //have just loaded the following
  // >cons_person1_read
  //  AAGACTGGAGAT
  // >cons_person2_read
  // ATCCTGGAGAA




  //notice the people share CTGGA
  //notice that each person basically has one big supernode, so the "breaking of supernodes" only occurs when you compare people.


  char subsection[24];

  dBNode* node = hash_table_find(element_get_key(seq_to_binary_kmer("GGAGA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(node != NULL);

  //test the subsections of person1's supernode, startng from the beginning

  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,1,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "AAGACT");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,2,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "AAGACTG");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,3,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "AAGACTGG");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,4,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "AAGACTGGA");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,5,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "AAGACTGGAG");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,6,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "AAGACTGGAGA");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,7,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "AAGACTGGAGAT");

  //what happens if you go too far?
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,8,individual_edge_array, 0, hash_table)==1);
  //what happens with negative numbers?
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,-1,individual_edge_array, 0, hash_table)==1);
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, -1,8,individual_edge_array, 0, hash_table)==1);

  //what happens if start=end - originally planned to return 1 in this case, but changed behaviour, and so had to fix test. function should work if start=end and return 0
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,0,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection,"AAGAC");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 1,1,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection,"AGACT");

  //what happens if start>end
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 5,4,individual_edge_array, 0, hash_table)==1);

  //now check a few others, not starting from 0
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 3,6,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "ACTGGAGA");

  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 2,3,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "GACTGG");


  //now do the same kind of thing for the other person

  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,1,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "ATCCTG");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,2,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "ATCCTGG");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,3,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "ATCCTGGA");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,4,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "ATCCTGGAG");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,5,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "ATCCTGGAGA");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,6,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "ATCCTGGAGAA");

  //what happens if you go too far?
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,18,individual_edge_array, 1, hash_table)==1);
  //what happens with negative numbers?
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,-1,individual_edge_array, 1, hash_table)==1);
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, -1,8,individual_edge_array, 1, hash_table)==1);

  //what happens if start=end - again - should not complain
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,0,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection,"ATCCT");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 1,1,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection,"TCCTG");

  //what happens if start>end
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 5,4,individual_edge_array, 1, hash_table)==1);
  
  //now check a few others, not starting from 0
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 1,3,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "TCCTGGA");

  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 3,5,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "CTGGAGA");
  

  hash_table_free(&hash_table);


}

void test_find_best_subsection_of_supernode_with_just_two_people()
{

  printf("\nplease reimplement (not urgent as no current users of this code) conssensus code using db_graph_supernode NULL TEST\n. ");


  //first set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 6;
  int bucket_size    = 4;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);


  char tmp_seq[hash_table->kmer_size];


  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }

  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/two_people_test_consensus.txt", &bad_reads, hash_table);
  //printf("Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 23);


  //have just loaded the following
  // >cons_person1_read
  //  AAGACTGGAGAT
  // >cons_person2_read
  // ATCCTGGAGAA

  //notice the people share CTGGA
  //notice that each person basically has one big supernode, so the "breaking of supernodes" only occurs when you compare people.

  //******** 1. *************************************************
  //Let's pick a node that happens to be in both people's graphs
  // ************************************************************
  dBNode* node = hash_table_find(element_get_key(seq_to_binary_kmer("TGGAG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(node != NULL);

  //OK, let's see if it can successfully find the best sub_supernode for each person.
  //There are only 2 people, so mke min people coverage 2, and don't have a min length.

  int index_of_start_of_best_sub_supernode;
  int length_of_best_sub_supernode;

  dBNode* first_node = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(node, individual_edge_array, 0, hash_table);
  char* first_kmer= binary_kmer_to_seq(first_node->kmer, hash_table->kmer_size, tmp_seq);
  //printf("Frst kmer is %s\n", first_kmer);
  CU_ASSERT( !strcmp(first_kmer,"GAGAT") || !strcmp(first_kmer,"ATCTC") || !strcmp(first_kmer,"AAGAC") || !strcmp(first_kmer,"GTCTT") );

  db_graph_get_best_sub_supernode_given_min_covg_and_length_for_specific_person_or_pop(first_node, &index_of_start_of_best_sub_supernode, &length_of_best_sub_supernode, 2, 0, individual_edge_array, 0, hash_table);

  //printf("\nindex of start is %d and length is %d", index_of_start_of_best_sub_supernode, length_of_best_sub_supernode);
  CU_ASSERT( ( index_of_start_of_best_sub_supernode==1)  || ( index_of_start_of_best_sub_supernode==4) ); 
  CU_ASSERT( length_of_best_sub_supernode==3);


  //OK - works for person1 (index 0) - what about person 2

  first_node = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(node, individual_edge_array, 1, hash_table);
  first_kmer= binary_kmer_to_seq(first_node->kmer, hash_table->kmer_size, tmp_seq);
  //printf("Frst kmer is %s\n", first_kmer);
  CU_ASSERT( !strcmp(first_kmer,"GAGAA") || !strcmp(first_kmer,"TTCTC") || !strcmp(first_kmer,"ATCCT") || !strcmp(first_kmer,"AGGAT") );

  db_graph_get_best_sub_supernode_given_min_covg_and_length_for_specific_person_or_pop(first_node, &index_of_start_of_best_sub_supernode, &length_of_best_sub_supernode, 2, 0, individual_edge_array, 1, hash_table);

  //printf("\nindex of start is %d and length is %d", index_of_start_of_best_sub_supernode, length_of_best_sub_supernode);
  CU_ASSERT( ( index_of_start_of_best_sub_supernode==1)  || ( index_of_start_of_best_sub_supernode==2) ); 
  CU_ASSERT( length_of_best_sub_supernode==3);



  hash_table_free(&hash_table);
  
  
}


// need to fix this so works for 3 people only
void test_get_population_consensus_supernode()
{


  printf("\nplease reimplement (not urgent as no current users of this code) conssensus code using db_graph_supernode NULL TEST\n. ");

   //first set up the hash/graph
   int kmer_size = 5;
   int number_of_bits = 6;
   int bucket_size    = 4;
   long long bad_reads = 0;
   int max_retries=10;
   
   dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
   
   
   char tmp_seq[hash_table->kmer_size];
   
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }

  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/five_people_test.txt", &bad_reads, hash_table);
  //printf("Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 155);


  // **** Have just loaded the following

  // person1: ATTGAGACTACCGTATTAACTAGGAGCACTG
  // person2:   TGAGA    CGTATTAACTAGGAGCACTG
  // person3:   TGAGA         TAACTAGGAGCACTG
  // person4:                                   identical to person 3
  // person5:   TGAGA                  GCACTG
  
  //spaces signify that they have other sequence in there that is different for each person.



  //All people share the node TGAGA. What is the population consensus supernode rooted at TGAGA?
  // That depends on what choice on minimum people coverage and minimum supernode length you demand.
  

  dBNode* node = hash_table_find(element_get_key(seq_to_binary_kmer("TGAGA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(node != NULL);


  Sequence* popseq_obj = malloc(sizeof(Sequence));
  if (popseq_obj==NULL)
    {
      printf("OOM alloc seq in test of pop consensus supernode");
      exit(1);
    }

  popseq_obj->seq= malloc(sizeof(char)*40);
  if (popseq_obj->seq==NULL)
    {
      printf("OOM alloc seq in test of pop consensus supernode");
      exit(1);
    }
  popseq_obj->name=NULL;
  popseq_obj->qual=NULL;

  //printf("*********\n\n***********\n\n Start of next test\n\n");
  //min covg 5, min length 6
  db_graph_find_population_consensus_supernode_based_on_given_node(popseq_obj,40,  node, 5, 6, hash_table);
  //printf("WE GET BACK %s\n", popseq_obj->seq);
  CU_ASSERT_STRING_EQUAL(popseq_obj->seq, "");

  //printf("*********\n\n***********\n\n Start of next test\n\n");

  //min covg 5, min length 5
  db_graph_find_population_consensus_supernode_based_on_given_node(popseq_obj, 40, node, 5, 5, hash_table);
  //printf("Answer is %s and expect TGAGA or TCTCA\n", popseq_obj->seq); 
  CU_ASSERT( !strcmp(popseq_obj->seq, "TGAGA") || !strcmp(popseq_obj->seq, "TCTCA") );


  //  printf("*********\n\n***********\n\n Start of next test\n\n");

  //min covg 4, min length 6
  db_graph_find_population_consensus_supernode_based_on_given_node(popseq_obj, 40, node, 4, 6, hash_table);
  //printf("Answer is %s and expect TAACTAGGA or TCCTAGTTA\n", popseq_obj->seq); 
  CU_ASSERT( !strcmp(popseq_obj->seq, "TAACTAGGA") || !strcmp(popseq_obj->seq, "TCCTAGTTA") );

  //  printf("*********\n\n***********\n\n Start of next test\n\n");

  //min covg 2, min length 14
  db_graph_find_population_consensus_supernode_based_on_given_node(popseq_obj, 40, node, 2, 14, hash_table);
  //printf("Answer is %s and expect CTGGCATCCTAGTTATCGTTAGAATCTCACC  or GGTGAGATTCTAACGATAACTAGGATGCCAG \n", popseq_obj->seq); 
  CU_ASSERT( !strcmp(popseq_obj->seq, "CTGGCATCCTAGTTATCGTTAGAATCTCACC") || !strcmp(popseq_obj->seq, "GGTGAGATTCTAACGATAACTAGGATGCCAG") );


  free_sequence(&popseq_obj);
  hash_table_free(&hash_table);

}

