
#include <CUnit.h>
#include <Basic.h>
#include <count_kmers.h>
#include "test_count_kmers.h"

#include <stdio.h>
#include <assert.h>

void test_count_kmers()
{

  //get output file 
  FILE* output = fopen("test_count_kmers_dir/test_log.txt","w");
  FILE* fp; 


   //1. take a fasta with one entry all A's and count 1,2,...length-of-sequence -mers  

  fp = fopen("test_count_kmers_dir/filelist_allAs", "r");
  CU_ASSERT_FALSE(count_kmers(fp,3,output));
  fclose(fp);
  

  //2. Do same for C,G,T

  fp = fopen("test_count_kmers_dir/filelist_allCs", "r");
  CU_ASSERT_FALSE(count_kmers(fp,3,output));
  fclose(fp);

  //note a file full of G's will lead to this program counting CCC..C kmers
  fp = fopen("test_count_kmers_dir/filelist_allGs", "r");
  CU_ASSERT_FALSE(count_kmers(fp,3,output));
  fclose(fp);

  fp = fopen("test_count_kmers_dir/filelist_allTs", "r");
  CU_ASSERT_FALSE(count_kmers(fp,3,output));
  fclose(fp);

  fclose(output);

  //Now checkt that these have all been read correctly.
  //TODO

  
  //3. Take a file that has a number of entries, each of which is 5 characters long, and count number of kmers totalled

  //4. Does it correctly count repeats that occur in different entries
  


  //5. What happens if one of the entries is smaller than the kmer length, but the others are all OK?

  //6. Is there any kind of odd-even difference

  //7. Is there a problem around k=32? <yes - not supported above 31>

  //8. Do you get the same answer if you cat together two fasta files as if you had them separately?

  //9. does it correctly count kmers when I insert an N into the middle - it should ignore kmers containing an N.

  //10. Take a known fastq where I have worked out kmer stats with perl, and compare results.



    fclose(output);
}


