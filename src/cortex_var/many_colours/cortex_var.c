/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo 
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */

#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <string.h>
#include <cmd_line.h>
#include <time.h>
#include <graph_info.h>
#include <db_differentiation.h>
#include <model_selection.h>
#include <experiment.h>
#include <genome_complexity.h>
#include <math.h>
#include <maths.h>

void timestamp();

/*
double  get_histogram_for_kmers_in_just_one_haplotype(long long** arr, int max, int in_colour, int not_in_colour, int data_colour, dBGraph* db_graph)
{
  long long total_kmers_possible=0;
  long long total_kmers_hit=0; 


  void analyse_kmer(dBGraph* db_graph, dBNode* node, int** array, int len, EdgeArrayType type, int colour)
  {
    if ( (node->coverage[in_colour]>0) && (node->coverage[not_in_colour]==0) )
      {
	total_kmers_possible++;
	//this node is informative about haplotype in_colour, so see how much covg it has in your sample data
	long long num = (long long) (node->coverage[data_colour]);
	if (num>0)
	  {
	    total_kmers_hit++;
	  }
	if (num<len)
	  {
	    *(array[num]) = *(array[num]) + 1;
	  }
	else
	  {
	    *(array[len]) = *(array[len]) + 1;
	  }
      }
  }
  db_graph_traverse_with_array_of_longlongs(&analyse_kmer, (HashTable*) db_graph, arr, max , individual_edge_array, -1);//the -1 is ignored
  return ((double) total_kmers_possible)/((double) total_kmers_hit);
}
*/

long long calculate_mean(long long* array, long long len)
{
  long long sum=0;
  long long num=0;
  long long i;
  for (i=0; i<len; i++)
    {
      sum += i*array[i];
      num += array[i];
    }
  return  (sum/num);
}

void run_novel_seq(CmdLine* cmd_line, dBGraph* db_graph, GraphAndModelInfo* model_info)
{
  
}


void run_pd_calls(CmdLine* cmd_line, dBGraph* db_graph, 
		  void (*print_some_extra_var_info)(VariantBranchesAndFlanks* var, FILE* fp),
		  GraphAndModelInfo* model_info)
{
  printf("Calling variants using Path Divergence Caller.\n");
  printf("Calls are made between the reference path as specified by the fasta in %s\n", cmd_line->ref_chrom_fasta_list);
  printf(" and the sample, which is the union of the following colour(s): ");
  int i;
  for (i=0; i< cmd_line->num_colours_in_pd_colour_list-1; i++)
    {
      printf("%d,", cmd_line->pd_colour_list[i]);
    }
  printf("%d\n", cmd_line->pd_colour_list[cmd_line->num_colours_in_pd_colour_list-1]);
  printf("The reference colour is %d\n", cmd_line->ref_colour);
  if (model_info->expt_type==Unspecified)
    {
      printf("Since you did not use --experiment_type, Cortex does not know how if each colour is a diploid/haploid sample or pool, so\nwill not calculate genotypes or likelihoods\n");
    }
  else if (cmd_line->genome_size==0)
    {
      printf("Since you did not specify the genome size/length, Cortex cannot calculate genotype likelihoods or call genotypes\n");
    }

  //this will also check that all the ref chrom fasta files exist
  int num_ref_chroms = get_number_of_files_and_check_existence_from_filelist(cmd_line->ref_chrom_fasta_list);

  char** ref_chroms = malloc( sizeof(char*) * num_ref_chroms);
  if (ref_chroms==NULL)
    {
      printf("OOM. Give up can't even allocate space for the names of the ref chromosome files\n");
      exit(1);
    }

  for (i=0; i< num_ref_chroms; i++)
    {
      ref_chroms[i] = malloc(sizeof(char)*500);
      if (ref_chroms[i]==NULL)
	{
	  printf("OOM. Giveup can't even allocate space for the names of the ref chromosome file i = %d\n",i);
	  exit(1);
	}
    }

  get_filenames_from_list(cmd_line->ref_chrom_fasta_list, ref_chroms, num_ref_chroms);

  //now set up output files

  char** output_files = malloc( sizeof(char*) * num_ref_chroms); //one for each chromosome
  if (output_files==NULL)
    {
      printf("OOM. Give up can't even allocate space for the names of the output  files \n");
      exit(1);
    }
	
  for (i=0; i< num_ref_chroms; i++)
    {
      output_files[i] = malloc(sizeof(char)*1000); 
      if (output_files[i]==NULL)
	{
	  printf("OOM. Giveup can't even allocate space for the names of the ref chromosome file i = %d\n",i);
	  exit(1);
	}
    }
  
  
  for (i=0; i<num_ref_chroms; i++)
    {
      sprintf(output_files[i], "%s_pd_chr_%i", cmd_line->path_divergence_caller_output_stub, i+1);
    }


  int min_fiveprime_flank_anchor = 2;
  int min_threeprime_flank_anchor= cmd_line->kmer_size;
  int max_anchor_span =  cmd_line -> max_var_len;
  int length_of_arrays = 2*max_anchor_span;
  int min_covg =1;
  int max_covg = 10000000;//this is ignored. will be changing API
  int max_expected_size_of_supernode=cmd_line -> max_var_len;
	
      
  for (i=0; i<num_ref_chroms; i++) 
    {
      printf("Call SV comparing individual with chromosome %s\n", ref_chroms[i]);
	    
      FILE* chrom_fptr = fopen(ref_chroms[i], "r");
      if (chrom_fptr==NULL)
	{
	  printf("Cannot open %s \n", ref_chroms[i]);
	  exit(1);
	}
      
      FILE* out_fptr = fopen(output_files[i], "w");
      if (out_fptr==NULL)
	{
	  printf("Cannot open %s for output\n", output_files[i]);
	  exit(1);
	}
      
      
      db_graph_make_reference_path_based_sv_calls_given_list_of_colours_for_indiv(cmd_line->pd_colour_list, cmd_line->num_colours_in_pd_colour_list,
										  chrom_fptr, cmd_line->ref_colour,
										  min_fiveprime_flank_anchor, min_threeprime_flank_anchor, 
										  max_anchor_span, min_covg, max_covg, 
										  max_expected_size_of_supernode, length_of_arrays, db_graph, out_fptr,
										  0, NULL, NULL, NULL, NULL, NULL, 
										  &make_reference_path_based_sv_calls_condition_always_true_in_subgraph_defined_by_func_of_colours, 
										  &db_variant_action_do_nothing,
										  print_some_extra_var_info, model_info);
      
      
      fclose(chrom_fptr);
      fclose(out_fptr);
    }
  //cleanup
  for(i=0; i<num_ref_chroms; i++)
    {
      free(ref_chroms[i]);
      free(output_files[i]);
    }
  free(ref_chroms);
  free(output_files);

  

}

void run_bubble_calls(CmdLine* cmd_line, int which, dBGraph* db_graph, 
		      void (*print_appropriate_extra_var_info)(VariantBranchesAndFlanks* var, FILE* fp),
		      Edges(*get_col_ref) (const dBNode* e),
		      int (*get_cov_ref)(const dBNode* e),
		      GraphInfo* db_graph_info, GraphAndModelInfo* model_info
		      )
{


  printf("Detecting bubbles between the union of this set of colours: ");
  int k;
  if (which==1)
	{
	  for (k=0; k<cmd_line->num_colours_in_detect_bubbles1_first_colour_list; k++)
	    {
	      printf("%d, ", cmd_line->detect_bubbles1_first_colour_list[k]);
	    }
	}
  else
    {
      for (k=0; k<cmd_line->num_colours_in_detect_bubbles2_first_colour_list; k++)
	{
	  printf("%d, ", cmd_line->detect_bubbles2_first_colour_list[k]);
	}
    }
  printf("\nand the union of this set of colours: ");
  if (which==1)
    {
      for (k=0; k<cmd_line->num_colours_in_detect_bubbles1_second_colour_list; k++)
	{
	  printf("%d, ", cmd_line->detect_bubbles1_second_colour_list[k]);
	}
    }
  else
    {
      for (k=0; k<cmd_line->num_colours_in_detect_bubbles2_second_colour_list; k++)
	{
	  printf("%d, ", cmd_line->detect_bubbles2_second_colour_list[k]);
	}
      printf("\n");
    }
  

  if (cmd_line->exclude_ref_bubbles==true)
    {
      printf("Will remove bubbles in the reference colour %d before doing any calling\n", cmd_line->ref_colour);
    }
  
  if (cmd_line->apply_model_selection_at_bubbles==true)
    {
      printf("Will compare likelihoods of two models (Repeat and  Variation (Hardy-Weinberg)) at all bubbles,\nand mark those more likely to be repeats for filtering\n");
    }
  
  if (model_info->expt_type==Unspecified)
    {
      printf("Since you did not use --experiment_type, Cortex does not know how if each colour is a diploid/haploid sample or pool, so\nwill not calculate genotypes or likelihoods\n");
    }
  else if (cmd_line->genome_size==0)
    {
      printf("Since you did not specify the genome size/length, Cortex cannot calculate genotype likelihoods\n");
    }
  FILE* fp;
  
  if (which==1)
    {
      fp = fopen(cmd_line->output_detect_bubbles1, "w");
    }
  else
    {
      fp = fopen(cmd_line->output_detect_bubbles2, "w");
    }
  
  if (fp==NULL)
    {
      if (which==1)
	{
	  printf("Cannot open %s. Exit.", cmd_line->output_detect_bubbles1);
	}
      else
	{
	  printf("Cannot open %s. Exit.", cmd_line->output_detect_bubbles2);
	}
      exit(1);
    }
  
  boolean (*mod_sel_criterion)(AnnotatedPutativeVariant* annovar,  GraphAndModelInfo* model_info)=NULL;
  
  if (cmd_line->apply_model_selection_at_bubbles==true)
    {
      mod_sel_criterion = &basic_model_selection_condition; 
    }
  
  
  if (cmd_line->exclude_ref_bubbles==true)
    {
      printf("(First exclude bubbles from ref colour %d) \n", cmd_line->ref_colour);
    }
  if (which==1)
    {
      db_graph_detect_vars_given_lists_of_colours(fp,cmd_line->max_var_len,db_graph, 
						  cmd_line->detect_bubbles1_first_colour_list, 
						  cmd_line->num_colours_in_detect_bubbles1_first_colour_list,
						  cmd_line->detect_bubbles1_second_colour_list, 
						  cmd_line->num_colours_in_detect_bubbles1_second_colour_list,
						  &detect_vars_condition_always_true, print_appropriate_extra_var_info,
						  cmd_line->exclude_ref_bubbles, get_col_ref, get_cov_ref, 
						  cmd_line->apply_model_selection_at_bubbles, mod_sel_criterion, 
						  model_info);
    }
  else
    {
      
      db_graph_detect_vars_given_lists_of_colours(fp,cmd_line->max_var_len,db_graph, 
						  cmd_line->detect_bubbles2_first_colour_list, 
						  cmd_line->num_colours_in_detect_bubbles2_first_colour_list,
						  cmd_line->detect_bubbles2_second_colour_list, 
						  cmd_line->num_colours_in_detect_bubbles2_second_colour_list,
						  &detect_vars_condition_always_true, print_appropriate_extra_var_info,
						  cmd_line->exclude_ref_bubbles, get_col_ref, get_cov_ref, 
						  cmd_line->apply_model_selection_at_bubbles, mod_sel_criterion,
						  model_info);
    }


}






int main(int argc, char **argv){

  timestamp();
  printf("Starting Cortex\n");

  CmdLine cmd_line = parse_cmdline(argc,argv,sizeof(Element));

  int hash_key_bits, bucket_size;
  dBGraph * db_graph = NULL;
  short kmer_size;

  //next two needed much later
  int num_kmers_dumped_after_alignment=0;
  char tmp_dump[300];


  //***************************************************************************
  //define local functions
  //***************************************************************************
  /*
  void apply_reset_to_all_edges(dBNode* node, Orientation or, Nucleotide nuc)
  {
    int j;
    for (j=0; j<NUMBER_OF_COLOURS; j++)
      {
	reset_one_edge(node, or, nuc, individual_edge_array, j);
      }
  }
  void apply_reset_to_all_edges_2(dBNode* node )
  {
    int j;
    for (j=0; j<NUMBER_OF_COLOURS; j++)
      {
	      db_node_reset_edges(node, individual_edge_array, j);
      }
  }
  */

  void print_appropriate_extra_variant_info(VariantBranchesAndFlanks* var, FILE* fp)
  {
    if (cmd_line.print_colour_coverages==true)
      {
	print_standard_extra_info(var, fp);
      }
    else
      {
	print_no_extra_info(var, fp);
      }
  }

  void print_appropriate_extra_supernode_info(dBNode** node_array, Orientation* or_array, int len, FILE* fout)
  {
    if (cmd_line.print_colour_coverages==true)
      {
	print_standard_extra_supernode_info(node_array, or_array, len, fout);
      }
    else
      {
	print_no_extra_supernode_info(node_array, or_array, len, fout);
      }
  }


  Edges get_colour_ref(const dBNode* e)
  {

    if (cmd_line.using_ref==false)
      {
	printf("Do not call get_colour_ref when --using_ref was not specified. Exiting - coding error\n");
	exit(1);
      }

    if ( (cmd_line.ref_colour<0) || (cmd_line.ref_colour>NUMBER_OF_COLOURS-1) )
      {
	printf("Calling get_colour_ref, but the reference colour %d has not been specified, or is > than the compile-time limit, of %d\n", 
	       cmd_line.ref_colour, NUMBER_OF_COLOURS-1);
	exit(1);
      }
    Edges ed = get_edge_copy(*e, individual_edge_array, cmd_line.ref_colour);
    return ed;
  }

  int get_covg_ref(const dBNode* e)
  {
    if (cmd_line.using_ref==false)
      {
	printf("Do not call get_coverage_ref when --using_ref was not specified. Exiting - coding error\n");
      }
    
    return e->coverage[cmd_line.ref_colour];
  }


  int get_covg_in_union_all_colours_except_ref(const dBNode* e)
  {
    int cov = element_get_covg_union_of_all_covgs(e);
    if (cmd_line.ref_colour!=-1)
      {
	cov -= e->coverage[cmd_line.ref_colour];
      }
    return cov;
  }

  //***************************************************************************
  // end local functions
  //***************************************************************************



  
  //set hash table variables:
  kmer_size        = cmd_line.kmer_size;
  hash_key_bits    = cmd_line.number_of_buckets_bits; //number of buckets: 2^hash_key_bits
  bucket_size      = cmd_line.bucket_size;

  int number_of_bitfields = ((kmer_size * 2) / (sizeof(bitfield_of_64bits)*8))+1;
  int max_kmer_size = (NUMBER_OF_BITFIELDS_IN_BINARY_KMER * sizeof(bitfield_of_64bits) * 4) -1;
  int min_kmer_size = ((NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1) * sizeof(bitfield_of_64bits) * 4) + 1;

  if (number_of_bitfields != NUMBER_OF_BITFIELDS_IN_BINARY_KMER){
    printf("K-mer %i  is not in current range of kmers [%i - %i] required for this executable.!\n",kmer_size,min_kmer_size,max_kmer_size);
    exit(0);
  }
  
  printf("Maximum k-mer size (compile-time setting): %ld\n", (NUMBER_OF_BITFIELDS_IN_BINARY_KMER * sizeof(bitfield_of_64bits) * 4) -1);
  
  if (cmd_line.kmer_size> (NUMBER_OF_BITFIELDS_IN_BINARY_KMER * sizeof(bitfield_of_64bits) * 4) -1){
    printf("k-mer size is too big [%i]!",cmd_line.kmer_size);
    exit(1);
  }
  printf("Actual K-mer size: %d\n", cmd_line.kmer_size);

  
  //Create the de Bruijn graph/hash table
  int max_retries=15;
  db_graph = hash_table_new(hash_key_bits,bucket_size, max_retries, kmer_size);
  if (db_graph==NULL)
    {
      printf("Giving up - unable to allocate memory for the hash table\n");
      exit(1);
    }
  printf("Hash table created, number of buckets: %d\n",1 << hash_key_bits);



  GraphInfo db_graph_info;
  graph_info_initialise(&db_graph_info);



  // input data:
  if (cmd_line.input_seq==true)
    {
      if (strcmp(cmd_line.se_list, "")!=0)
	{
	  fprintf(stdout,"Input file of single ended data filenames: %s\n",
		  cmd_line.se_list); 
	}
      else
	{
	  printf("No SE data\n");
	}
      if (strcmp(cmd_line.pe_list_lh_mates, "") !=0)
	{
	  fprintf(stdout,"Input file of paired end data: %s, and %s, \n",
		  cmd_line.pe_list_lh_mates, cmd_line.pe_list_rh_mates); 
	}
      else
	{
	  printf("No paired-end data\n");
	}
      if (cmd_line.quality_score_threshold>0)
	{
	  fprintf(stdout,"quality cut-off: %i\n",cmd_line.quality_score_threshold);
	}
      else
	{
	  fprintf(stdout, "No quality filtering\n");
	}
      if (cmd_line.remove_pcr_dups==true)
	{
	  printf("Removing duplicates from the paired end files when both mates start with the same kmer\n");
	}
      else
	{
	  printf("No PCR duplicate removal\n");
	}
      if (cmd_line.cut_homopolymers==true)
	{
	  printf("Breaking reads at homopolymer runs of length %d or greater. Read restarts at first base after the run\n", cmd_line.homopolymer_limit);
	}

      boolean there_is_se_data=false;
      boolean there_is_pe_data=false;
      if (strcmp(cmd_line.se_list, "")!=0)
	{
	  there_is_se_data=true;
	}
      if (strcmp(cmd_line.pe_list_lh_mates, "") != 0)
	{
	  there_is_pe_data=true;
	}
      
      //note, if we load fasta/fastq data it always goes into colour 0;
      long long bases_parsed=0;
      long long bases_pass_filters_and_loaded=0;
      //get read-len distribution (after filters):
      long long* readlen_distrib=(long long*) malloc(sizeof(long long) * (cmd_line.max_read_length+1));
      long long** readlen_distrib_ptrs = (long long**) malloc(sizeof(long long*) * (cmd_line.max_read_length+1));

      if (readlen_distrib==NULL)
	{
	  printf("Unable to malloc array to hold readlen distirbution!Exit.\n");
	  exit(1);
	}
      int i;
      for (i=0; i<=cmd_line.max_read_length; i++)
	{
	  readlen_distrib[i]=0;
	  readlen_distrib_ptrs[i]=&readlen_distrib[i];
	}

      timestamp();

      load_se_and_pe_filelists_into_graph_of_specific_person_or_pop(there_is_se_data, there_is_pe_data, 
								    cmd_line.se_list, cmd_line.pe_list_lh_mates, cmd_line.pe_list_rh_mates,
								    &bases_parsed, &bases_pass_filters_and_loaded,readlen_distrib_ptrs,
								    cmd_line.quality_score_threshold, false, cmd_line.remove_pcr_dups,
								    cmd_line.cut_homopolymers, cmd_line.homopolymer_limit, cmd_line.quality_score_offset,
								    cmd_line.format_of_input_seq,
								    cmd_line.max_read_length, 0, db_graph);

      //update the graph info object
      if (cmd_line.format_of_input_seq==FASTQ)
	{
	  graph_info_update_mean_readlen_and_total_seq(&db_graph_info, 0, calculate_mean(readlen_distrib, (long long) (cmd_line.max_read_length+1)), bases_pass_filters_and_loaded);
	}
      else//for FASTA we do not get read length distribution
	{
	  printf("When binaries are built, Cortex calculates the legnth distribution of reads\n");
	  printf("after quality filters, Ns, PCR duplicates, homopolymer filters have been taked into account\n");
	  printf("This data is saved in the binary header, if you dump a binary file.\n");
	  printf("However Cortex does not calculate mean read length of input data if it is FASTA\n");
	  printf("These stats are primarily used for \"real\" data from fastq files\n");
	  printf("Since you have used fasta, we just use the value you entered for --max_read_len and log that in the binary header as the read length\n");
	  printf("Set mean read len in colour 0 to %d\n", cmd_line.max_read_length);
	  graph_info_set_mean_readlen(&db_graph_info, 0, cmd_line.max_read_length);
	  graph_info_increment_seq(&db_graph_info, 0, bases_pass_filters_and_loaded);
	}



      //cleanup marks left on nodes by loading process (for PE reads)
      hash_table_traverse(&db_node_set_status_to_none, db_graph);


      timestamp();
      if (cmd_line.format_of_input_seq==FASTQ)
	{
	  printf("Fastq data loaded\nTotal bases parsed:%qd\nTotal bases passing filters and loaded into graph:%qd\nMean read length after filters applied:%d\n", 
		 bases_parsed, bases_pass_filters_and_loaded, db_graph_info.mean_read_length[0]);
	}
      else
	{
	  printf("Fasta data loaded\nTotal bases parsed:%qd\nTotal bases passing filters and loaded into graph:%qd\n", 
		 bases_parsed, bases_pass_filters_and_loaded);
	}

      if (cmd_line.dump_readlen_distrib==true)
	{
	  FILE* rd_distrib_fptr = fopen(cmd_line.readlen_distrib_outfile, "w");
	  if (rd_distrib_fptr==NULL)
	    {
	      printf("Cannot open %s, so will dumping distribution of filtered read-lengths to stdout\n", cmd_line.readlen_distrib_outfile);
	      for (i=db_graph->kmer_size; i<=cmd_line.max_read_length; i++)
		{
		  printf("%d\t%qd\n", i, readlen_distrib[i]);
		}
	    }
	  else
	    {
	      printf("Dumping distribution of effective read lengths (ie after quality, homopolymer and/or PCR duplicate filters to file %s.\n", cmd_line.readlen_distrib_outfile);
	      for (i=db_graph->kmer_size; i<=cmd_line.max_read_length; i++)
		{
		  fprintf(rd_distrib_fptr, "%d\t%qd\n", i, readlen_distrib[i]);
		}
	      fclose(rd_distrib_fptr);
	    }
	  
	}


      free(readlen_distrib_ptrs);
      free(readlen_distrib);
      
    }
  else
    {
      //if there is a multicolour binary, load that in first
      timestamp();      
      int first_colour_data_starts_going_into=0;
      boolean graph_has_had_no_other_binaries_loaded=true;

      if (cmd_line.input_multicol_bin==true)
	{
	  int mean_readlens[NUMBER_OF_COLOURS];
	  int* mean_readlens_ptrs[NUMBER_OF_COLOURS];
	  long long total_seq_in_that_colour[NUMBER_OF_COLOURS];
	  long long* total_seq_in_that_colour_ptrs[NUMBER_OF_COLOURS];
	  int j;
	  for (j=0; j<NUMBER_OF_COLOURS; j++)
	    {
	      mean_readlens[j]=0;
	      mean_readlens_ptrs[j]=&(mean_readlens[j]);
	      total_seq_in_that_colour[j]=0;
	      total_seq_in_that_colour_ptrs[j]=&(total_seq_in_that_colour[j]);
	    }
	  long long  bp_loaded = load_multicolour_binary_from_filename_into_graph(cmd_line.multicolour_bin,db_graph, &first_colour_data_starts_going_into,
										  mean_readlens_ptrs, total_seq_in_that_colour_ptrs);
	  //update graph_info object
	  for (j=0; j<first_colour_data_starts_going_into; j++)
            {
	      graph_info_update_mean_readlen_and_total_seq(&db_graph_info, j, mean_readlens[j], total_seq_in_that_colour[j]);
	    }
	  timestamp();
	  printf("Loaded the multicolour binary %s, and got %qd kmers\n", cmd_line.multicolour_bin, bp_loaded/db_graph->kmer_size);
	  graph_has_had_no_other_binaries_loaded=false;
	}

      if (cmd_line.input_colours==true)
	{
	  timestamp();

	  //normal use
	  if (cmd_line.successively_dump_cleaned_colours==false)
	    {
	      
	      printf("List of colours: %s (contains one filelist per colour). Load data into consecutive colours starting at %d\n", 
		     cmd_line.colour_list, first_colour_data_starts_going_into);
	      if (cmd_line.load_colours_only_where_overlap_clean_colour==true)
		{
		  printf("When loading the binaries specified in %s, we only load nodes that are already in colour %d\n", 
			 cmd_line.colour_list, cmd_line.clean_colour);
		}
	      load_population_as_binaries_from_graph(cmd_line.colour_list, first_colour_data_starts_going_into, 
						     graph_has_had_no_other_binaries_loaded, db_graph, &db_graph_info,
						     cmd_line.load_colours_only_where_overlap_clean_colour, cmd_line.clean_colour);
	      timestamp();
	      printf("Finished loading single_colour binaries\n");
	    }
	  else
	    {
	      //we have loaded a multicolour binary, and we have checked that the clean_colour is one of the colours in that binary
	      if (cmd_line.load_colours_only_where_overlap_clean_colour==false)
		{
		  printf("If you specify --successively_dump_cleaned_colours, you must also specify --load_colours_only_where_overlap_clean_colour\n");
		  printf("That should fix your problem, however, this should have been caught as soon as Cortex parsed your command-line. Please inform Zam Iqbal (zam@well.ox.ac.uk) so he can fix that UI bug\n");
		  exit(1);
		}
	      printf("For each colour in %s, load data into graph, cleaning by comparison with colour %d, then dump a single-colour binary\n",
		     cmd_line.colour_list,cmd_line.clean_colour);
	      dump_successive_cleaned_binaries(cmd_line.colour_list, first_colour_data_starts_going_into,cmd_line.clean_colour,
					       cmd_line.successively_dump_cleaned_colours_suffix, db_graph, &db_graph_info);
	      printf("Completed dumping of clean binaries\n");
	    }



	}
    }




  GraphAndModelInfo model_info;
  float repeat_geometric_param_mu = 0.8;
  float seq_err_rate_per_base;

  if (cmd_line.manually_override_error_rate==false)
    {
      seq_err_rate_per_base= 0.01;//default, but this may change below

      /*      if (cmd_line.genome_size !=0)
	{
	  long long num_cov1_kmers = db_graph_count_covg1_kmers_in_func_of_colours(db_graph, &get_covg_in_union_all_colours_except_ref);
	  long long num_cov2_kmers = db_graph_count_covg2_kmers_in_func_of_colours(db_graph, &get_covg_in_union_all_colours_except_ref);

	  double estim_err = num_cov1_kmers/(num_cov2_kmers*(db_graph->kmer_size) );
	  if ((estim_err<0.1) && (estim_err>0.001) )
	    {
	      printf("Estimating sequencing error rate of %f, using estimate num_kmers_covg_1/(kmer_size * num_kmers_covg_2), where num_kmers_covg1 = %qd, num_kmers_covg2 = %qd\n", estim_err, num_cov1_kmers, num_cov2_kmers);
	      printf("This is used in error-cleaning and likelihood calcs. Remember you can override this with --estimated_error_rate\n");
	      seq_err_rate_per_base=estim_err;
	    }
	  else
	    {
	      printf("Attempted to estimate the sequencing error rate, using estimate num_kmers_covg_1/(kmer_size * num_kmers_covg_2), where num_kmers_covg1 = %qd, num_kmers_covg2 = %qd. However the estimate we got, %f, was not believable so fell back on default value of 0.01. This is used in error-cleaning and likelihood calcs. Remember you can override this with --estimated_error_rate. One reason why this estimate may have failed is if the loaded binary was already cleaned\n",  num_cov1_kmers, num_cov2_kmers, estim_err);
	    }
	}
      */
    }
  else
    {
      seq_err_rate_per_base=cmd_line.manually_entered_seq_error_rate;
    }

  int num_chroms_in_expt=NUMBER_OF_COLOURS;
  if (cmd_line.expt_type==EachColourADiploidSample)
    {
      num_chroms_in_expt=2*NUMBER_OF_COLOURS;
    }
  else if (cmd_line.expt_type==EachColourADiploidSampleExceptTheRefColour)
    {
      num_chroms_in_expt=2*NUMBER_OF_COLOURS-2;
    }
  else if (cmd_line.expt_type==EachColourAHaploidSample)
    {
      num_chroms_in_expt=NUMBER_OF_COLOURS;
    }
  else if (cmd_line.expt_type==EachColourAHaploidSampleExceptTheRefColour)
    {
      num_chroms_in_expt=NUMBER_OF_COLOURS-1;
    }
  else if (cmd_line.expt_type==Unspecified)
    {
      num_chroms_in_expt=NUMBER_OF_COLOURS;
    }
  else
    {
      printf("Coding error - expt type not specified - not even as Unspecified");
      exit(1);
    }

  initialise_model_info(&model_info, &db_graph_info, cmd_line.genome_size, 
			repeat_geometric_param_mu, seq_err_rate_per_base, cmd_line.ref_colour, num_chroms_in_expt, cmd_line.expt_type);


  int j;
  /*
  printf("The following mod is for the sims for the paper only: i need the graphinfo to contain read lengths for the fasta based sim binaries\n");

  for (j=0; j<NUMBER_OF_COLOURS; j++)
    {
      if (db_graph_info.mean_read_length[j]==0)
	{
	  graph_info_set_mean_readlen(&db_graph_info, j, cmd_line.max_read_length);
	}
    }
  */

      
  printf("Total kmers in table: %qd\n", hash_table_get_unique_kmers(db_graph));	  
  printf("****************************************\n");
  printf("SUMMARY:\nColour:\tMeanReadLen\tTotalSeq\n");

  for (j=0; j<NUMBER_OF_COLOURS; j++)
    {
      printf("%d\t%d\t%qd\n", j, db_graph_info.mean_read_length[j], db_graph_info.total_sequence[j]);
    }
  printf("****************************************\n");

  if (cmd_line.health_check==true)
    {
      timestamp();
      printf("Run health check on loaded graph\n");
      db_graph_health_check(false, db_graph);
      printf("End of health check\n");
      timestamp();
    }


  
  // Error Correction actions
  if (cmd_line.remove_seq_errors==true)
    {
      timestamp();
      printf("Remove nodes that look like sequencing errors. Clip tips first\n");
      db_graph_clip_tips_in_union_of_all_colours(db_graph);
      
      /*
	printf("Then remove low coverage supernodes covg (<= %d) \n", cmd_line.remv_low_covg_sups_threshold);
      db_graph_remove_errors_considering_covg_and_topology(cmd_line.remv_low_covg_sups_threshold,db_graph, &element_get_covg_union_of_all_covgs, &element_get_colour_union_of_all_colours,
							   &apply_reset_to_specific_edge_in_union_of_all_colours, &apply_reset_to_all_edges_in_union_of_all_colours,      
						   cmd_line.max_var_len);
      */
      db_graph_remove_supernodes_more_likely_errors_than_sampling(db_graph, &db_graph_info, &model_info,
								  cmd_line.max_var_len, 
								  &element_get_covg_union_of_all_covgs, &element_get_colour_union_of_all_colours,
								  &apply_reset_to_specific_edge_in_union_of_all_colours, &apply_reset_to_all_edges_in_union_of_all_colours,      
								  cmd_line.max_var_len, cmd_line.remv_low_covg_sups_threshold);
      timestamp();
      printf("Error correction done\n");

    }
  else if (cmd_line.remv_low_covg_sups_threshold!=-1)
    {
      printf("Clip tips first\n");
      db_graph_clip_tips_in_union_of_all_colours(db_graph);

      printf("Remove low coverage supernodes covg (<= %d) \n", cmd_line.remv_low_covg_sups_threshold);
      db_graph_remove_errors_considering_covg_and_topology(cmd_line.remv_low_covg_sups_threshold,db_graph, &element_get_covg_union_of_all_covgs, &element_get_colour_union_of_all_colours,
							   &apply_reset_to_specific_edge_in_union_of_all_colours, &apply_reset_to_all_edges_in_union_of_all_colours,
							   cmd_line.max_var_len);
      timestamp();
      printf("Error correction done\n");

    }
  else if (cmd_line.remove_low_coverage_nodes==true)
    {
      timestamp();
      printf("Start to to remove nodes with covg (in union of all colours)  <= %d\n", cmd_line.node_coverage_threshold);
      db_graph_remove_low_coverage_nodes_ignoring_colours(cmd_line.node_coverage_threshold, db_graph);
      timestamp();
      printf("Error correction done\n");
      
    }

  if (cmd_line.health_check==true)
    {
      timestamp();
      printf("Run health check on cleaned graph\n");
      db_graph_health_check(false, db_graph);
      printf("Health check done\n");
      timestamp();
    }

  if (cmd_line.dump_binary==true)
    {
      if (cmd_line.input_seq==true)
	{
	  //dump single colour
	  timestamp();
	  printf("Input data was fasta/q, so dump single colour binary file: %s\n", cmd_line.output_binary_filename);
	  db_graph_dump_single_colour_binary_of_colour0(cmd_line.output_binary_filename, &db_node_check_status_not_pruned,db_graph, &db_graph_info);
	  timestamp();
	  printf("Binary dumped\n");


	}
      else
	{
	  timestamp();
	  printf("Dump multicolour binary with %d colours (compile-time setting)\n", NUMBER_OF_COLOURS);
	  db_graph_dump_binary(cmd_line.output_binary_filename, &db_node_check_status_not_pruned,db_graph, &db_graph_info);
	  timestamp();
	  printf("Binary dumped\n");
	}
    }

  if (cmd_line.print_supernode_fasta==true)
    {

      timestamp();


      /*
      printf("Get histograms of covg on supernodes specific to either PGF or COX. Calculate number of reads/length for each sup, and put in bins\n");
      
      int* bins_PGF = (int*) malloc( sizeof(int) * 10000);
      int* bins_COX = (int*) malloc( sizeof(int) * 10000);
      if ( (bins_PGF==NULL) || (bins_COX==NULL) )
	{
	  printf("Unable to alloc bins for PGF and COX\n");
	  exit(1);
	}
      int p;
      for (p=0; p<10000; p++)
	{
	  bins_PGF[p]=0;
	  bins_COX[p]=0;
	}

      db_graph_get_stats_of_supernodes_that_split_two_colour(cmd_line.max_var_len, 0,1, db_graph, &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs,
							     &condition_splits_COX_and_PGF, bins_PGF, bins_COX);


      FILE* fout = fopen(cmd_line.knight_output, "w");
      if (fout==NULL)
	{
	  printf("Cannot open output fule for knight\n");
	  exit(1);
	}
      fprintf(fout, "Covg\tPGF_not_COX\tCOX_not_PGF\n");

      for (p=0; p<10000; p++)
	{
	  fprintf(fout, "%d\t%d\t%d\n", p, bins_PGF[p], bins_COX[p]);
	}
      fclose(fout);
      */

      timestamp();
      printf("Finished COX PGF Julian Knight analysis\n");


    }


  if (cmd_line.dump_covg_distrib==true)
    {
      timestamp();
      printf("Dump kmer coverage distribution for colour 0 to file %s\n", cmd_line.covg_distrib_outfile);
      db_graph_get_covg_distribution(cmd_line.covg_distrib_outfile, db_graph, individual_edge_array, 0, &db_node_check_status_not_pruned);
      timestamp();
      printf("Covg distribution dumped\n");
    }

  // DETECT BUBBLES

  if (cmd_line.detect_bubbles1==true)
    {
      timestamp();
      printf("Start first set of bubble calls\n");
      run_bubble_calls(&cmd_line, 1, db_graph, &print_appropriate_extra_variant_info, &get_colour_ref, &get_covg_ref, &db_graph_info, &model_info);

      //unset the nodes marked as visited, but not those marked as to be ignored
      hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
      timestamp();
      printf("Detect Bubbles 1, completed\n");
    }

  //second detect bubbles
  if (cmd_line.detect_bubbles2==true)
    {
      timestamp();
      printf("Start second set of bubble calls\n");
      run_bubble_calls(&cmd_line, 2, db_graph, &print_appropriate_extra_variant_info, &get_colour_ref, &get_covg_ref, &db_graph_info, &model_info);
      //unset the nodes marked as visited, but not those marked as to be ignored
      hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
      timestamp();
      printf("Detect Bubbles 2, completed\n");
    }

  if (cmd_line.make_pd_calls==true)
    {
      timestamp();
      printf("Run Path-Divergence Calls\n");
      run_pd_calls(&cmd_line, db_graph, &print_appropriate_extra_variant_info, &model_info);
      hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
      timestamp();
      printf("Finished Path Divergence calls\n");
    }
  if (cmd_line.align_given_list==true)
    {
      timestamp();
      printf("Start aligning the fasta/q listed in this file: %s\n", cmd_line.list_fastaq_to_align);
      int array_of_colours[NUMBER_OF_COLOURS];
      int j;
      char* array_of_colournames[NUMBER_OF_COLOURS];
      for (j=0; j<NUMBER_OF_COLOURS; j++)
	{
	  array_of_colours[j]=j;
	  array_of_colournames[j]=(char*)malloc(sizeof(char) * 50);
	  if (array_of_colournames[j]==NULL)
	    {
	      printf("Severe lack of memory. Cannot even allocate 50 chars. Give up\n");
	      exit(1);
	    }
	  sprintf(array_of_colournames[j], "colour_%d", j);
	}
      align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours(cmd_line.format_of_files_to_align, cmd_line.list_fastaq_to_align,
								       cmd_line.max_read_length, array_of_colours, array_of_colournames,
								       NUMBER_OF_COLOURS,db_graph,cmd_line.quality_score_offset,
								       false, NULL, NULL, cmd_line.dump_aligned_overlap_binary);
      for (j=0; j<NUMBER_OF_COLOURS; j++)
	{
	  free(array_of_colournames[j]);
	}
      printf("Completed alignment of fasta/q to graph to print coverages in all colours\n");
      printf("Dumping a binary, %s,  of all the nodes which were hit by the alignment process\n", cmd_line.output_aligned_overlap_binname);
      sprintf(tmp_dump, "%s.temporary_delete_me", cmd_line.output_aligned_overlap_binname);
      printf("In the process we have to create a temporary file, %s, which you can/should delete when cortex has completed\n", tmp_dump);

      num_kmers_dumped_after_alignment = db_graph_dump_binary(tmp_dump, &db_node_check_status_to_be_dumped, db_graph, &db_graph_info);
      hash_table_traverse(&db_node_action_set_status_of_unpruned_to_none, db_graph);	

    }
  if (cmd_line.print_colour_overlap_matrix==true)
    {
      timestamp();
      printf("Do a direct graph comparison of each colour x,y,z.. listed in --colour_subgraph_overlap_matrix a,b,c../x,y,z... against the colours a,b,c.. (process is symmetrical)\n");
      
      db_graph_print_colour_overlap_matrix(cmd_line.colour_overlap_first_colour_list,
					   cmd_line.num_colours_in_colour_overlap_first_list,
					   cmd_line.colour_overlap_second_colour_list,
					   cmd_line.num_colours_in_colour_overlap_second_list,
					   db_graph);
      printf("\nCompleted graph overlap matrix printing\n");
    }

  if (cmd_line.knight_expt==true)
    {

      printf("Start COX PGF Julian Knight analysis\n");

      FILE* fout = fopen(cmd_line.knight_output, "w");
      if (fout==NULL)
	{
	  printf("Cannot open output fule for knight\n");
	  exit(1);
	}
      
      call_bubbles_distinguishing_cox_pgf(db_graph, 10000, fout, &model_info);
      fclose(fout);

      /*
      //this supernode is seen in pure COX RNA but not in pure PGF RNA, or vice-versa
      //PGF pure RNA is colour 0, and pure COX RNA is colour 1
      boolean condition_splits_COX_and_PGF_RNA(dBNode** path, int length, int* which_col)
      {
	int col_pure_PGF_RNA=0;
	int col_pure_COX_RNA=1;

	boolean all_nodes_in_COX=true;
	boolean no_nodes_in_COX = true;
	boolean all_nodes_in_PGF=true;
	boolean no_nodes_in_PGF=true;
	int p;
	for (p=1; p<length; p++)
	  {
	    if (db_node_get_coverage(path[p], individual_edge_array, col_pure_PGF_RNA)==0)
	      {
		all_nodes_in_PGF=false;
	      }
	    if (db_node_get_coverage(path[p], individual_edge_array, col_pure_COX_RNA)==0)
	      {
		all_nodes_in_COX=false;
	      }
	    if (db_node_get_coverage(path[p], individual_edge_array, col_pure_PGF_RNA)>0)
	      {
		no_nodes_in_PGF=false;
	      }
	    if (db_node_get_coverage(path[p], individual_edge_array, col_pure_COX_RNA)>0)
	      {
		no_nodes_in_COX=false;
	      }
	  }
	if ( (all_nodes_in_COX==true) && (no_nodes_in_PGF==true) )
	  {
	    *which_col = 1;
	    return true;
	  }
	else if ( (all_nodes_in_PGF==true) && (no_nodes_in_COX==true) )
	  {
	    *which_col=0;
	    return true;
	  }
	else
	  {
	    *which_col=-1;
	    return false;
	  }

      }

      //traverse supernodes in union of colour-pure-PGF_RNA and colour-pure-COX-RNA, and for each supernode,
      //calculate number of reads in pure COX and pure PGF, and add to two running totals, T_COX and T_PGF
      //while doing this, get running totals for my mixed-data. T1, T2,..T5 (T1 refers to 1:1 mixing, T2 to 1:1.125 mixing , etc)

      int tot_pgf=0;
      int tot_cox=0;
      int tot_data1=0;
      int tot_data2=0;
      int tot_data3=0;
      int tot_data4=0;
      int tot_data5=0;
      int max_length_sup=10000;
      db_graph_get_covgs_in_all_colours_of_col0union1_sups(max_length_sup, db_graph,
							   &tot_pgf, &tot_cox, 
							   &tot_data1,&tot_data2,&tot_data3,&tot_data4,&tot_data5);


      //then traverse the same supernodes one more time, and this time for each supernode,
      //IF it is distinct to one or other (COX or PGF), then print ratio of number of reads on it, to T_COX or T_PGF as appropriate.
      // and then for each data colour, print number of reads in that colour, to T1, T2 or whatever, as appropriate
      //The idea being for each supernode that splits the two haplotypes, we get an estimate of the mixing ratio.

      db_graph_get_proportion_of_cvg_on_each_sup(max_length_sup, db_graph, tot_pgf, tot_cox,
						 tot_data1, tot_data2, tot_data3, tot_data4, tot_data5,
						 &condition_splits_COX_and_PGF_RNA, fout);
      fclose(fout);

      */

      printf("Finished J Knight analysis\n");
    }
  

  if (cmd_line.genotype_complex_site==true)
    {
      printf("Genotyping a complex site with %d known alleles\n", cmd_line.num_alleles_of_site);
      printf("Of all possible genotypes, we will calculate likelihoods for those numbered %d to %d\n",
	     cmd_line.first_genotype_to_calc_likelihoods_for, cmd_line.last_genotype_to_calc_likelihoods_for);
      printf("We do this for these samples (which we assume are diploid): colours ");
      int k;
      for (k=0; k<cmd_line.num_colours_to_genotype; k++)
	{
	  printf("%d,", cmd_line.list_colours_to_genotype[k]);
	}
      printf("\n");

      //sanity check before we start
      //check that we have read lengths for the colours we want to genotype
      for (k=0; k<cmd_line.num_colours_to_genotype; k++)
	{
	  if (db_graph_info.mean_read_length[cmd_line.list_colours_to_genotype[k]] < db_graph->kmer_size )
	    {
	      printf("This will not work. If you scroll up to the summary of read-lengths and covgs in your colours, you will see that\n");
	      printf("at least one of the colours you want to genotype has mean read length < kmer_size. We only know of one way this can happen:\n");
	      printf("If you load fastA files, then cortex does not store read-length/covg data, essentially because the only real use-cases for\n");
	      printf("fasta files are a)testing of code and b) reference genomes. So I suggest that either you are\n");
	      printf("running a test on fasta files (in which case please use fastq) or you have accidentally specified the wrong colour to be genotyped\n");
	      exit(1);
	    }

	  if (db_graph_info.total_sequence[cmd_line.list_colours_to_genotype[k]]==0)
	    {
	      printf("This will not work. If you scroll up to the summary of read-lengths and covgs in your colours, you will see that\n");
	      printf("at least one of the colours you want to genotype has total sequence 0. We only know of one way this can happen:\n");
	      printf("If you load fastA files, then cortex does not store read-length/covg data, essentially because the only real use-cases for\n");
	      printf("fasta files are a)testing of code and b) reference genomes. So I suggest that either you are\n");
	      printf("running a test on fasta files (in which case please use fastq) or you have accidentally specified the wrong colour to be genotyped\n");
	      exit(1);
	    }
	}

      double* current_max_lik_array         = alloc_ML_results_array(cmd_line.num_colours_to_genotype);
      double* current_max_but_one_lik_array = alloc_ML_results_array(cmd_line.num_colours_to_genotype);
      char** name_current_max_lik_array             = alloc_ML_results_names_array(cmd_line.num_colours_to_genotype);
      char** name_current_max_but_one_lik_array     = alloc_ML_results_names_array(cmd_line.num_colours_to_genotype);

      calculate_max_and_max_but_one_llks_of_specified_set_of_genotypes_of_complex_site(cmd_line.list_colours_to_genotype, cmd_line.num_colours_to_genotype,
										       cmd_line.colour_of_reference_with_site_excised,
										       cmd_line.num_alleles_of_site,
										       cmd_line.first_genotype_to_calc_likelihoods_for,
										       cmd_line.last_genotype_to_calc_likelihoods_for,
										       cmd_line.max_var_len, cmd_line.fasta_alleles_for_complex_genotyping,
										       cmd_line.assump_for_genotyping,
										       current_max_lik_array, current_max_but_one_lik_array,
										       name_current_max_lik_array, name_current_max_but_one_lik_array,
										       true, &model_info, db_graph, cmd_line.working_colour1, cmd_line.working_colour2,
										       cmd_line.using_1net, cmd_line.using_2net,
										       cmd_line.min_acceptable_llk
										       );


    }
  if (cmd_line.estimate_genome_complexity==true)
    {
      int num_reads_used_in_estimate=0;
      double g = estimate_genome_complexity(db_graph, cmd_line.fastaq_for_estimating_genome_complexity,
					    true, 0, 1,cmd_line.max_read_length, cmd_line.format_of_input_seq,
					    cmd_line.quality_score_offset, &num_reads_used_in_estimate);
      printf("We estimate genome complexity at k=%d (for SNPs) as %f\n", db_graph->kmer_size, g);
      printf("This estimate used a sample of %d high-quality reads\n", num_reads_used_in_estimate);

    }

  if (cmd_line.print_novel_contigs==true)
    {
      timestamp();
      printf("Start to search for and print novel contigs\n");
      printf("Definition of novel: contig must lie in union of these colours: ");
      int j;
      for (j=0; j<cmd_line.numcols_novelseq_colours_search; j++)
	{
	  printf("%d,", cmd_line.novelseq_colours_search[j]);
	}
      printf("\n and at least %d percent of the kmers in this contig (excluding first and last) must have zero coverage in the union of these colours: ", cmd_line.novelseq_min_percentage_novel);
      for (j=0; j<cmd_line.numcols_novelseq_colours_avoid; j++)
	{
	  printf("%d,", cmd_line.novelseq_colours_avoid[j]);
	}
      printf("\nAlso contig must be at least %d bp long\n", cmd_line.novelseq_contig_min_len_bp);
      db_graph_print_novel_supernodes(cmd_line.novelseq_outfile, cmd_line.max_var_len, db_graph, 
				      cmd_line.novelseq_colours_search, cmd_line.numcols_novelseq_colours_search,
				      cmd_line.novelseq_colours_avoid, cmd_line.numcols_novelseq_colours_avoid,
				      cmd_line.novelseq_contig_min_len_bp, cmd_line.novelseq_min_percentage_novel,
				      &print_appropriate_extra_supernode_info);
      timestamp();
      printf("Finished printing novel contigs\n");
      
    }

  
  hash_table_free(&db_graph);
  timestamp();

  if (cmd_line.dump_aligned_overlap_binary==true)
    {

      //reload the binary you dumped, clean off the edges, and then dump again.
      //malloc a new hash table. Only needs to be as large as you need.
      float s = log(num_kmers_dumped_after_alignment/100)/log(2); //make width 100
      // now 2^s = num_kmers_dumped_after_alignment/100, so s is my height
      dBGraph* db_graph2;
      boolean try_smaller_hash=true;
      if (2*num_kmers_dumped_after_alignment < pow(2,hash_key_bits) * bucket_size)
	{
	  db_graph2 = hash_table_new(s+1,100, max_retries, kmer_size);
	  if (db_graph2==NULL)
	    {
	      try_smaller_hash=false;
	    }
	}
      if (try_smaller_hash==false)
	{
	  db_graph2 = hash_table_new(hash_key_bits,bucket_size, max_retries, kmer_size);
	  if (db_graph2==NULL)
	    {
	      printf("Cortex has nearly finished. It's done everything you asked it to do, and has dumped a binary of the overlap of your alignment with the graph. However, by \"ripping out\" nodes from the main graph, that dumped binary now has edges pointing out to nodes that are not in the binary. So the idea is that we have deallocated the main graph now, and we were going to load the dumped binary, clean it up and re-dump it. However that has failed, becaause we could not malloc the memory to do it - the most likely reason is that someone else is sharing your server and their mempory use has gone up.\n");
	      exit(1);
	    }
	}	  
      
      //this is all for the API - we wont use this info
      int mean_readlens2[NUMBER_OF_COLOURS];
      int* mean_readlens_ptrs2[NUMBER_OF_COLOURS];
      long long total_seq_in_that_colour2[NUMBER_OF_COLOURS];
      long long* total_seq_in_that_colour_ptrs2[NUMBER_OF_COLOURS];
      int j;
      for (j=0; j<NUMBER_OF_COLOURS; j++)
	{
	  mean_readlens2[j]=0;
	  mean_readlens_ptrs2[j]=&(mean_readlens2[j]);
	  total_seq_in_that_colour2[j]=0;
	  total_seq_in_that_colour_ptrs2[j]=&(total_seq_in_that_colour2[j]);
	}
      int num_c;//number of colours in binary
      
      long long  n  = load_multicolour_binary_from_filename_into_graph(tmp_dump,db_graph2, &num_c,
								       mean_readlens_ptrs2, total_seq_in_that_colour_ptrs2);
      printf("Loaded the temporary binary %s and found %qd kmers\n", tmp_dump, n);
      //this is why we are going to all this bother - cleaning edges
      db_graph_clean_orphan_edges(db_graph2);
      db_graph_dump_binary(cmd_line.output_aligned_overlap_binname, 
			   &db_node_condition_always_true,
			   db_graph2, &db_graph_info);//deliberately using original graph info - we want the same info
      hash_table_free(&db_graph2);
    }
  




  printf("Cortex completed - have a nice day!\n");
  return 0;
}



void timestamp(){
 time_t ltime;
 ltime = time(NULL);
 printf("\n-----\n%s",asctime(localtime(&ltime)));
 fflush(stdout);
}
