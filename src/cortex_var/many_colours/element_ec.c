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
/*
  element.c - defines the interface for the de Bruijn graph node. The
  implementation is complemented by a hash table that stores every node indexed
  by kmers (BinaryKmers). 

  The element routines, ie the one required by hash_table/priority queue, are
  prefixed with element_

  The de Bruijn based routines are prefixed with db_node_ec_
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <string.h>

// cortex_var headers
#include "element_ec.h"

//const Covg COVG_MAX = UINT_MAX;
//const Covg COVG_MAX = INT_MAX;

// Only print covg overflow warning once
//char overflow_warning_printed = 0;

//currently noone calls this in normal use
// In normal use, the priority queue allocates space to put the eloement directly within,
// and calls element_initialise
ElementEc* new_element_ec()
{
  ElementEc* e = malloc(sizeof(ElementEc));

  if (e==NULL)
  {
    die("Unable to allocate a new element");
  }

  binary_kmer_initialise_to_zero(&(e->kmer));

  int i;
  for (i=0; i< NUMBER_OF_COLOURS; i++)
  {
    e->individual_edges[i]=0;
  }

  e->status = (char)ec_none;

  return e;
}


void free_element_ec(ElementEc** element)
{
  free(*element);
  *element = NULL;
}

void element_ec_assign(ElementEc* e1, ElementEc* e2)
{
  binary_kmer_assignment_operator(e1->kmer, e2->kmer);

  int i;
  for (i=0; i< NUMBER_OF_COLOURS; i++)
  {
    e1->individual_edges[i] = e2->individual_edges[i];
  }

  e1->status = e2->status;
}


//return a copy of the edge you are referring to
Edges element_ec_get_edge_copy(const ElementEc e, int colour)
{
  assert(colour < NUMBER_OF_COLOURS);

  return e.individual_edges[colour];
}


Edges element_ec_get_union_of_edges(ElementEc e)
{
  int i;
  Edges edges=0;

  for (i=0; i< NUMBER_OF_COLOURS; i++)
  {
    edges |= e.individual_edges[i];
  }

  return edges;
}

Edges element_ec_get_colour_union_of_all_colours(const ElementEc* e)
{
  int i;
  Edges edges=0;
  
  for (i=0; i< NUMBER_OF_COLOURS; i++)
  {
    edges |= e->individual_edges[i];
  }

  return edges;
  
}



Edges element_ec_get_last_colour(const ElementEc* e)
{
  Edges edges =  element_ec_get_edge_copy(*e, NUMBER_OF_COLOURS-1);
  return edges;
}


Edges element_ec_get_colour0(const ElementEc* e)
{
  Edges edges=element_ec_get_edge_copy(*e,0);
  return edges;
}

Edges element_ec_get_colour1(const ElementEc* e)
{
  Edges edges=element_ec_get_edge_copy(*e,1);
  return edges;
}



Covg element_ec_get_covg_union_of_all_covgs(const dBNodeEc* e)
{
  Covg sum_covg = 0;
  
  return sum_covg;
}

Covg element_ec_get_covg_for_colourlist(const dBNodeEc* e, int* colour_list,
                                     int list_len)
{
  Covg sum_covg = 0;
  
  return sum_covg;
}


Covg element_ec_get_covg_colour0(const dBNodeEc* e)
{
  return 0;
}

Covg element_ec_get_covg_last_colour(const dBNodeEc* e)
{
  return 0;
}

#if NUMBER_OF_COLOURS > 1

Covg element_ec_get_covg_colour1(const dBNodeEc* e)
{
  return 0;
}

#endif /* NUMBER_OF_COLOURS > 1 */



//adds edges from edge_char to the appropriate person/population edgeset, without removing existing edges
void element_ec_add_edges(ElementEc* e, int colour, Edges edge_char)
{
  assert(colour < NUMBER_OF_COLOURS);

  e->individual_edges[colour] |= edge_char;
}


void element_ec_set_edges(ElementEc* e, int colour, Edges edge_char)
{
  assert(colour < NUMBER_OF_COLOURS);

  e->individual_edges[colour] = edge_char;
}


void db_node_ec_reset_all_edges_for_all_people_and_pops_to_zero(ElementEc* e)
{
  int i;

  for (i=0; i<NUMBER_OF_COLOURS; i++)
  {
    e->individual_edges[i]=0;
  }
}

void element_ec_reset_one_edge(ElementEc* e, Orientation orientation,
                    Nucleotide nucleotide, int colour)
{
  assert(colour < NUMBER_OF_COLOURS);

  char edge = 1 << nucleotide;

  if (orientation == reverse)
  {
    edge <<= 4;
  }

  // toggle 1->0 0->1
  // xor with all 1's, ie 00010000 -> 11101111
  edge ^= (unsigned char) 0xFF;

  // reset one edge
  e->individual_edges[colour] &= edge;
}

// DEV: why is this using edges and not coverage?
int element_ec_get_number_of_people_or_pops_containing_this_element(ElementEc* e)
{
  int i;
  int count = 0;
  
  for(i = 0; i < NUMBER_OF_COLOURS; i++)
  {
    if(e->individual_edges[i] != 0)
    {
      count++;
    }
  }

  return count;
}

boolean element_ec_smaller(ElementEc  e1, ElementEc e2)
{ 
  return element_ec_get_union_of_edges(e1)  <  element_ec_get_union_of_edges(e2); 
}



// WARNING: this gives you a pointer to a the binary kmer in the node.
// You could modify contents of the hash table
BinaryKmer* element_ec_get_kmer(ElementEc * e)
{
  return &(e->kmer);
}

boolean element_ec_is_key(Key key, ElementEc e)
{
  //  return key == e.kmer;
  return binary_kmer_comparison_operator(*key, e.kmer);
}

// For a given kmer, get the BinaryKmer 'key':
// the lower of the kmer vs reverse complement of itself
Key element_ec_get_key(BinaryKmer* kmer, short kmer_size, Key preallocated_key)
{
  // Get first and last nucleotides
  Nucleotide first = binary_kmer_get_first_nucleotide(kmer, kmer_size);
  Nucleotide last  = binary_kmer_get_last_nucleotide(kmer);
  Nucleotide rev_last = ~last & 0x3;

  if(first < rev_last)
  {
    binary_kmer_assignment_operator(*((BinaryKmer*)preallocated_key), *kmer);
  }
  else if(first > rev_last)
  {
    binary_kmer_reverse_complement(kmer, kmer_size, (BinaryKmer*)preallocated_key );
  }
  else
  {
    // Don't know which is going to be correct
    // This will happen 1 in 4 times
    binary_kmer_reverse_complement(kmer, kmer_size, (BinaryKmer*)preallocated_key );

    if(binary_kmer_less_than(*kmer, *((BinaryKmer*)preallocated_key), kmer_size))
    {
      binary_kmer_assignment_operator( *((BinaryKmer*)preallocated_key) , *kmer);
    }
  }

  return preallocated_key;
}


void element_ec_initialise(ElementEc * e, Key kmer, short kmer_size){

  if (e==NULL)
    {
      die("Called element_ec_initialise on NULL ptr");
    }

  BinaryKmer tmp_kmer;
  binary_kmer_initialise_to_zero(&tmp_kmer);
  binary_kmer_assignment_operator( e->kmer, *(element_ec_get_key(kmer, kmer_size, &tmp_kmer)));

  //hash table has calloc-ed all elements, so elements from the hash table are already initialised to zero.
  //however this function is used to reset to 0 ElementEcs that are reused,
  // - see below in the read_binary functions. Also in tests.

  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      e->individual_edges[i]=0;
    }

  db_node_ec_set_status(e, ec_none);
}



void element_ec_initialise_kmer_covgs_edges_and_status_to_zero(ElementEc * e)
{
  binary_kmer_initialise_to_zero(&(e->kmer));

  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
  {
    e->individual_edges[i] = 0;
  }

    printf("sizeof(Element)=%lu\n", sizeof(ElementEc));
    printf("sizeof(e->kmer)=%lu\n", sizeof(e->kmer));
    //printf("sizeof(e->coverage)=%lu\n", sizeof(e->coverage));
    printf("sizeof(e->individual_edges)=%lu\n", sizeof(e->individual_edges));
    //printf("sizeof(e->status)=%lu\n", sizeof(e->status));

  db_node_ec_set_status(e, ec_none);
}


void element_ec_set_kmer(ElementEc * e, Key kmer, short kmer_size)
{
  BinaryKmer tmp_kmer;
  binary_kmer_initialise_to_zero(&tmp_kmer);
  binary_kmer_assignment_operator( e->kmer, *(element_ec_get_key(kmer, kmer_size, &tmp_kmer)));
}



void db_node_ec_increment_coverage(dBNodeEc* e, int colour)
{
  db_node_ec_update_coverage(e, colour, 1);
}

void db_node_ec_update_coverage(dBNodeEc* e, int colour, long update)
{
  assert(colour < NUMBER_OF_COLOURS);
}


// Apparently some code calls to db_node_ec_get_coverage pass
// dBNodeEc e == NULL and so we need to test for NULL *sometimes*.
// For an example see: count_reads_on_allele_in_specific_colour
// (expects zero covg if it passes a NULL)
// Probably safer (catches bugs) to explicitly have a function that tolerates
// NULLs and one that throws an error
Covg db_node_ec_get_coverage_tolerate_null(const dBNodeEc* const e, int colour)
{
  assert(colour < NUMBER_OF_COLOURS);

  return 0;
}

Covg db_node_ec_get_coverage(const dBNodeEc* const e, int colour)
{
  assert(e != NULL);
  assert(colour < NUMBER_OF_COLOURS);

  return 0;
}


void db_node_ec_set_coverage(dBNodeEc* e, int colour, Covg covg)
{
  //e->coverage[colour] = covg;
}


// DEV: why is this needed?
// replace with (e == NULL ? 0 : get_covg(e))
Covg db_node_ec_get_coverage_in_subgraph_defined_by_func_of_colours(const dBNodeEc* const e, Covg (*get_covg)(const dBNodeEc*) )
{

  if (e==NULL)
    {
      return 0;
    }
  else
    {
      return get_covg(e);
    }
}



Orientation element_ec_opposite_orientation(Orientation o)
{
  return o ^ 1;
}

Orientation db_node_ec_get_orientation(BinaryKmer* k, dBNodeEc * e, short kmer_size){

  if (binary_kmer_comparison_operator(e->kmer,*k)==true)
    {
      return forward;
    }
  
  BinaryKmer tmp_kmer;

  if (binary_kmer_comparison_operator(e->kmer, *(binary_kmer_reverse_complement(k,kmer_size, &tmp_kmer)))==true)
    {
      return reverse;
    }

  char tmpseq1[kmer_size+1];
  char tmpseq2[kmer_size+1];

  die("programming error - you have called  db_node_ec_get_orientation with a kmer\n"
      "that is neither equal to the kmer in this node, nor its rev comp\n"
      "Arg 1 Kmer is %s and Arg 2 node kmer is %s\n",
      binary_kmer_to_seq(k, kmer_size, tmpseq1),
      binary_kmer_to_seq(&(e->kmer), kmer_size, tmpseq2));
}






// After specifying which individual or population you are talking about, this
// function adds one edge ("arrow") to the appropriate edge in the appropriate
// array in the element -- basically sets a bit in the correct edges char
void db_node_ec_add_labeled_edge(dBNodeEc * e, Orientation o, Nucleotide base, int edge_index)
{
  // set edge
  // A (0) -> 0001, C (1) -> 0010, G (2) -> 0100, T (3) -> 1000
  char edge = 1 << base;
  
  if (o == reverse){
    edge <<= 4; //move to next nibble 
  }

  //update node
  element_ec_add_edges(e, edge_index, edge);
}


//adding an edge between two nodes implies adding two labeled edges (one in each direction)
//be aware that in the case of self-loops in palindromes the two labeled edges collapse in one
// DEV: improve this!
boolean db_node_ec_add_edge(dBNodeEc * src_e, dBNodeEc * tgt_e,
                         Orientation src_o, Orientation tgt_o,
                         short kmer_size, int colour)
{
  BinaryKmer src_k, tgt_k, tmp_kmer; 
  char seq1[kmer_size+1];
  char seq2[kmer_size+1];

  binary_kmer_assignment_operator(src_k, src_e->kmer);
  binary_kmer_assignment_operator(tgt_k, tgt_e->kmer);

  if (src_o == reverse){
    binary_kmer_assignment_operator(src_k, *(binary_kmer_reverse_complement(&src_k,kmer_size, &tmp_kmer)));
  }
    
  if (tgt_o == reverse){
    binary_kmer_assignment_operator(tgt_k, *(binary_kmer_reverse_complement(&tgt_k,kmer_size, &tmp_kmer)));
  }
    
  
  if (DEBUG){
    printf("add edge %s -%c-> %s to edge type %d, and edge index %d\n",
	   binary_kmer_to_seq(&src_k,kmer_size, seq1),
	   binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(&tgt_k)),
	   binary_kmer_to_seq(&tgt_k,kmer_size, seq2), 0, colour);
  }

  db_node_ec_add_labeled_edge(src_e, src_o, binary_kmer_get_last_nucleotide(&tgt_k), colour);

  if (DEBUG){

    printf("add edge %s -%c-> %s to edge type %d, and edge index %d\n",
	   binary_kmer_to_seq(&tgt_k,kmer_size,seq1),
	   binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(binary_kmer_reverse_complement(&src_k,kmer_size, &tmp_kmer))),
	   binary_kmer_to_seq(&src_k,kmer_size, seq2),  0, colour);
  }

  Nucleotide base = binary_kmer_get_last_nucleotide(binary_kmer_reverse_complement(&src_k,kmer_size, &tmp_kmer));
  db_node_ec_add_labeled_edge(tgt_e, element_ec_opposite_orientation(tgt_o), base, colour);

  return true;
}



boolean db_node_ec_edge_exist(dBNodeEc * element, Nucleotide base,
                           Orientation orientation, int colour)
{
  //get the edge char for this specific person or pop:
  char edge = element_ec_get_edge_copy(*element, colour);

  edge >>= base;

  if (orientation == reverse)
  {
    edge >>= 4;
  }
  
  return (edge & 0x1);
}


//final argument f is a function that returns an Edge that is a function of the different colured edges in a node.
// eg we might be interested in the union of all the coloured edges in the graph, or just the colour/edge for the first person,
//    or the union of all edges except that corresponding to the reference.
boolean db_node_ec_edge_exist_within_specified_function_of_coloured_edges(
  dBNodeEc *element, Nucleotide base,
  Orientation orientation, Edges (*f)( const ElementEc* ))
{
  char edge = f(element);

  edge >>= base;

  if (orientation == reverse)
  {
    edge >>= 4;
  }

  return (edge & 1);
}


void db_node_ec_reset_edges(dBNodeEc * node, int colour)
{
  element_ec_set_edges(node, colour, 0);
}

void db_node_ec_reset_edge(dBNodeEc * node, Orientation orientation, Nucleotide nucleotide, int edge_index)
{
  element_ec_reset_one_edge(node, orientation, nucleotide, edge_index);
}

//pass in a function which will apply reset_one_edge to whichever of the edges it cares about
void db_node_ec_reset_specified_edges(dBNodeEc * node, Orientation orientation,
                                   Nucleotide nucleotide,  
				                           void (*f)(dBNodeEc*, Orientation, Nucleotide))
{
  f(node, orientation, nucleotide);
}


boolean db_node_ec_edges_reset(dBNodeEc *node, int colour)
{
  return (element_ec_get_edge_copy(*node, colour) == 0);
}


boolean db_node_ec_has_precisely_one_edge(dBNodeEc *node, Orientation orientation,
                                       Nucleotide *nucleotide, int colour)
{
  Nucleotide n;
  Edges edges = element_ec_get_edge_copy(*node, colour);
  short edges_count = 0;

  if (orientation == reverse){
    edges >>= 4;
  }

  for(n=0;n<4;n++)
  {
    if ((edges & 1) == 1)
    {
      *nucleotide = n;
      edges_count++;
    }
    
    edges >>= 1;    
  }
  
  return (edges_count == 1);
}



boolean db_node_ec_has_precisely_one_edge_in_subgraph_defined_by_func_of_colours(
  dBNodeEc * node, Orientation orientation, Nucleotide * nucleotide, 
  Edges (*get_colour)(const dBNodeEc*))
{
  
  Nucleotide n;
  Edges edges = get_colour(node);
  short edges_count = 0;

  if (orientation == reverse){
    edges >>= 4;
  }
 
  
  for(n=0;n<4;n++){
    
    if ((edges & 1) == 1){
      *nucleotide = n;
      edges_count++;
    }
    
    edges >>= 1;    
  }
  
  return (edges_count == 1);
  
}


boolean db_node_ec_has_precisely_one_edge_in_union_graph_over_all_people(
  dBNodeEc *node, Orientation orientation, Nucleotide *nucleotide)
{
  Nucleotide n;
  Edges edges = element_ec_get_union_of_edges(*node);
  short edges_count = 0;

  if (orientation == reverse){
    edges >>= 4;
  }
 
  
  for(n=0;n<4;n++){
    
    if ((edges & 1) == 1){
      *nucleotide = n;
      edges_count++;
    }
    
    edges >>= 1;    
  }
  
  return (edges_count == 1);
  
}

//a conflict - bifurcation
boolean db_node_ec_has_precisely_two_edges(dBNodeEc *node, Orientation orientation,
                                        Nucleotide *nucleotide1,
                                        Nucleotide *nucleotide2, int colour)
{
  Nucleotide n;

  Edges edges = element_ec_get_edge_copy(*node, colour);

  short edges_count = 0;

  if (orientation == reverse){
    edges >>= 4;
  }
  
  for(n=0;n<4;n++){
    
    if ((edges & 1) == 1){
      if (edges_count == 0){
	*nucleotide1 = n;
      }
      
      if (edges_count == 1){
	*nucleotide2 = n;
      }
      edges_count++;
    }
    
    edges >>= 1;    
  }
  
  return (edges_count == 2);
}



boolean db_node_ec_is_blunt_end(dBNodeEc * node, Orientation orientation, int colour)
{  
  Edges edges = element_ec_get_edge_copy(*node, colour);

  if (orientation == reverse)
  {
    edges >>= 4;
  }

  // AND with 00001111 so that we only look at the 4 least significant bits
  edges &= 15;

  return edges == 0;
}



boolean db_node_ec_is_blunt_end_in_subgraph_given_by_func_of_colours(
  dBNodeEc * node, Orientation orientation,  Edges (*get_colour)(const dBNodeEc*))
{
  Edges edges = get_colour(node);

  if (orientation == reverse)
  {
    edges >>= 4;
  }
  
  edges &= 15; // AND with 00001111 so that we only look at the 4 least significant bits
  
  return (edges == 0);
}



boolean db_node_ec_check_status(dBNodeEc * node, NodeStatusEc status)
{
  return (node->status == (char)status);
}

boolean db_node_ec_check_status_to_be_dumped(dBNodeEc * node)
{
  if ( db_node_ec_check_status(node, ec_to_be_dumped)==true )
    {
      return true;
    }
  else
    {
      return false;
    }
}

boolean db_node_ec_check_status_not_pruned(dBNodeEc * node)
{
  if ( db_node_ec_check_status(node, ec_none) || db_node_ec_check_status(node, ec_visited))
    {
      return true;
    }
  return false;
}


boolean db_node_ec_check_status_not_pruned_or_visited(dBNodeEc * node)
{
  if(db_node_ec_check_status(node,ec_visited) ||
     db_node_ec_check_status(node,ec_visited_and_exists_in_reference) ||
     db_node_ec_check_status(node, ec_pruned))
  {
    return false;
  }
  else
  {
    return true;
  }
}


void db_node_ec_set_status(dBNodeEc *node, NodeStatusEc status)
{
  node->status = (char)status;
}

void db_node_ec_set_status_to_none(dBNodeEc * node)
{
  node->status = (char)ec_none;
}




boolean db_node_ec_is_this_node_in_this_person_or_populations_graph(dBNodeEc* node, int colour)
{
  if(node == NULL)
  {
    return false;
  }

  if(db_node_ec_check_status(node, ec_pruned) == true)
  {
    return false;
  }

  Edges edge_for_this_person_or_pop = element_ec_get_edge_copy(*node, colour);

  return (edge_for_this_person_or_pop != 0);
}

boolean db_node_ec_is_this_node_in_subgraph_defined_by_func_of_colours(dBNodeEc* node, Edges (*get_colour)(const dBNodeEc*))
{

  if (node ==NULL)
  {
    return false;
  }

  if (db_node_ec_check_status(node, ec_pruned)==true)
  {
    return false;
  }

  // get_colour will return an edge which is some function of all the 
  // edges/colours at that node. eg, will AND all the edges, to consider the
  // union of all colours, etc
  Edges edge = get_colour(node);

  return (edge != 0);
}

// DEV: combine some of these functions, remove local copying
//prints all colours
void db_node_ec_print_multicolour_binary(FILE * fp, dBNodeEc * node)
{
  assert(node != NULL);

  BinaryKmer kmer;
  binary_kmer_assignment_operator(kmer, *element_ec_get_kmer(node) );
  Covg covg[NUMBER_OF_COLOURS];
  Edges individual_edges[NUMBER_OF_COLOURS]; 

  int i;
  for (i=0; i< NUMBER_OF_COLOURS; i++)                                                     
  {
    covg[i]             = (Covg) db_node_ec_get_coverage(node, i);
    individual_edges[i] = element_ec_get_edge_copy(*node, i);
  }
				  
  fwrite(kmer, NUMBER_OF_BITFIELDS_IN_BINARY_KMER*sizeof(bitfield_of_64bits), 1, fp);
  fwrite(covg, sizeof(uint32_t), NUMBER_OF_COLOURS, fp); 
  fwrite(individual_edges, sizeof(Edges), NUMBER_OF_COLOURS, fp);
}


void db_node_ec_print_single_colour_binary_of_colour0(FILE * fp, dBNodeEc * node)
{
  assert(node != NULL);

  BinaryKmer kmer;
  binary_kmer_assignment_operator(kmer, *element_ec_get_kmer(node) );
  Covg covg;
  Edges individual_edges; 

  covg             = (Covg) db_node_ec_get_coverage(node, 0);
  individual_edges = element_ec_get_edge_copy(*node, 0);
  
  fwrite(kmer, NUMBER_OF_BITFIELDS_IN_BINARY_KMER*sizeof(bitfield_of_64bits), 1, fp);
  fwrite(&covg, sizeof(uint32_t), 1, fp); 
  fwrite(&individual_edges, sizeof(Edges), 1, fp);

  //do NOT remove this. For some reason prevents an OS/disk
  //issue which results in writing zeroes. Basically have no idea why, but this fixes it.
  fflush(fp); 
              
}


void db_node_ec_print_single_colour_binary_of_specified_colour(FILE * fp, dBNodeEc * node, int colour)
{
  assert(node != NULL);

  BinaryKmer kmer;
  binary_kmer_assignment_operator(kmer, *element_ec_get_kmer(node) );
  Covg covg = db_node_ec_get_coverage(node, colour);
  Edges individual_edges = element_ec_get_edge_copy(*node, colour);

  fwrite(kmer, NUMBER_OF_BITFIELDS_IN_BINARY_KMER*sizeof(bitfield_of_64bits), 1, fp);
  fwrite(&covg, sizeof(uint32_t), 1, fp); 
  fwrite(&individual_edges, sizeof(Edges), 1, fp);
}

/*
boolean db_node_ec_read_multicolour_binary(FILE * fp, short kmer_size, dBNodeEc * node){

  BinaryKmer kmer;
  binary_kmer_assignment_operator(kmer, *(element_ec_get_kmer(node)) );
  int covg[NUMBER_OF_COLOURS];
  Edges individual_edges[NUMBER_OF_COLOURS]; 

  int read;

  read = fread(&kmer,sizeof(bitfield_of_64bits)*NUMBER_OF_BITFIELDS_IN_BINARY_KMER,1,fp);  

  if (read>0){


    read = fread(covg, sizeof(int), NUMBER_OF_COLOURS, fp);    
    if (read==0){
      die("error with input file - failed to read covg in db_node_ec_read_sv_trio_binary. You have tried to read an incompatible binary - this error message should not happen - you should have hit other checks first.. Please contact mario.caccamo@bbsrc.ac.uk and zam@well.ox.ac.uk\n");
    }

    read = fread(individual_edges, sizeof(Edges), NUMBER_OF_COLOURS, fp);
    if (read==0){
      die("error with input file - failed to read Edges in db_node_ec_read_sv_trio_binary. You have tried to read an incompatible binary - this error message should not happen - you should have hit other checks first.. Please contact mario.caccamo@bbsrc.ac.uk and zam@well.ox.ac.uk\n");
    }



  }
  else{
    return false;
  }

  element_ec_initialise(node,&kmer,kmer_size);

  int i;
  for (i=0; i< NUMBER_OF_COLOURS; i++)
    {
      node->coverage[i]         = covg[i];
      node->individual_edges[i] = individual_edges[i];
    }

  return true;
}
*/



// assumes we have the same NUMBER_OF_BITFIELDS_IN_BINARY_KMER in both this hash table and the binary we are loading
// this is checked when first opening the binary file
// assumes we ar loading a single multicolour binary into an empty hash table, so initialises each node to 0 before 
// adding info from the binary
boolean db_node_ec_read_multicolour_binary(FILE * fp, short kmer_size, dBNodeEc * node, int num_colours_in_binary, int binversion_in_binheader)
{

  if ( (num_colours_in_binary>NUMBER_OF_COLOURS)
    ||
       (num_colours_in_binary<=0)
       )
    {
      die("You should not call db_node_ec_read_multicolour_binary with %d as final argument.\n"
          "NUMBER_OF_COLOURS is %d", num_colours_in_binary, NUMBER_OF_COLOURS);
    }

  
  if (binversion_in_binheader==4)//legacy
    {
      BinaryKmer kmer;
      int covg_reading_from_binary[num_colours_in_binary];
      Edges    individual_edges_reading_from_binary[num_colours_in_binary];
      int read;
      read = fread(&kmer,sizeof(bitfield_of_64bits)*NUMBER_OF_BITFIELDS_IN_BINARY_KMER,1,fp);  
      
      if (read>0){
	
        read = fread(covg_reading_from_binary, sizeof(int), num_colours_in_binary, fp);    
        if (read==0){
          die("error with input file - failed to read covg in db_node_ec_read_multicolour_binary\n");
        }
        
        read = fread(individual_edges_reading_from_binary, sizeof(Edges), num_colours_in_binary, fp);
        if (read==0){
          die("error with input file - failed to read Edges in db_node_ec_read_multicolour_binary\n");
        }
      }
      else{
        return false;
      }

      //element_ec_initialise(node,&kmer,kmer_size);
      element_ec_set_kmer(node,&kmer,kmer_size);

      int i;
      for (i=0; i< num_colours_in_binary; i++)
        {
          node->individual_edges[i] = individual_edges_reading_from_binary[i];
        }      
    }
  else//higher than 4
    {
      BinaryKmer kmer;
      Covg covg_reading_from_binary[num_colours_in_binary];
      Edges    individual_edges_reading_from_binary[num_colours_in_binary];
      int read;      
      read = fread(&kmer,sizeof(bitfield_of_64bits)*NUMBER_OF_BITFIELDS_IN_BINARY_KMER,1,fp);  
      
      if (read>0){
	
        read = fread(covg_reading_from_binary, sizeof(uint32_t), num_colours_in_binary, fp);    
        if (read==0){
          die("Failed to read covg in db_node_ec_read_multicolour_binary\n");
        }
        
        read = fread(individual_edges_reading_from_binary, sizeof(Edges), num_colours_in_binary, fp);
        if (read==0){
          die("Failed to read Edges in db_node_ec_read_multicolour_binary\n");
        }
      }
      else{
        return false;
      }
      
      //element_ec_initialise(node,&kmer,kmer_size);
      element_ec_set_kmer(node,&kmer,kmer_size);

      int i;
      for (i=0; i< num_colours_in_binary; i++)
        {
          node->individual_edges[i] = individual_edges_reading_from_binary[i];
        }      
    }
  
  return true;
}



// the edge array type and index tell you which person you should load this data into
boolean db_node_ec_read_single_colour_binary(FILE *fp, short kmer_size, dBNodeEc *node,
                                          int index, int binversion_in_binheader)
{

  if ( (index<0) || (index>=NUMBER_OF_COLOURS))
    {
      die("Invalid index for which person to load binary into: %d. Exiting.", index);
    }

  if (binversion_in_binheader==4)//legacy
    {
      BinaryKmer kmer;
      Edges edges;
      int coverage;
      int read;
      
      read = fread(&kmer,sizeof(bitfield_of_64bits)*NUMBER_OF_BITFIELDS_IN_BINARY_KMER,1,fp);
      
      if (read>0){
	read = fread(&coverage,sizeof(int),1,fp);    
	if (read==0){
	  die("Failed to read cvg, incompatible binary format?");
	}
	read = fread(&edges,sizeof(Edges),1,fp);
	if (read==0){
	  die("Failed to read edges, incompatible binary format?");
	}   
      }
      else{
	return false;
      }
      
      element_ec_set_kmer(node,&kmer,kmer_size);
      //element_ec_initialise(node,&kmer,kmer_size);
      node->individual_edges[index]    = edges;
      db_node_ec_action_set_status_none(node);
      
    }
  else//proper BINARY VERSION
    {
      BinaryKmer kmer;
      Edges edges;
      Covg coverage;
      int read;
      
      read = fread(&kmer,sizeof(bitfield_of_64bits)*NUMBER_OF_BITFIELDS_IN_BINARY_KMER,1,fp);
      
      if (read>0){
	read = fread(&coverage,sizeof(uint32_t),1,fp);    
	if (read==0){
	  die("Failed to read cvg, incompatible binary format?");
	}
	read = fread(&edges,sizeof(Edges),1,fp);
	if (read==0){
	  die("Failed to read edges, incompatible binary format?");
	}   
      }
      else{
	return false;
      }
      
      element_ec_set_kmer(node,&kmer,kmer_size);
      //element_ec_initialise(node,&kmer,kmer_size);
      node->individual_edges[index]    = edges;  
      db_node_ec_action_set_status_none(node);    
    }
  
  return true;
}


void db_node_ec_action_set_status_none(dBNodeEc * node){
  db_node_ec_set_status(node,ec_none);
}

void db_node_ec_action_set_status_of_unpruned_to_none(dBNodeEc * node){

  if (db_node_ec_check_status(node, ec_pruned)==false)
    {
      db_node_ec_set_status(node,ec_none);
    }
}



void db_node_ec_action_set_status_pruned(dBNodeEc * node){
  db_node_ec_set_status(node,ec_pruned);
}


void db_node_ec_action_set_status_visited(dBNodeEc * node){
  db_node_ec_set_status(node,ec_visited);
}

void db_node_ec_action_set_status_visited_unless_marked_to_be_ignored(dBNodeEc * node)
{
  if (db_node_ec_check_status(node, ec_ignore_this_node)==false)
    {
      db_node_ec_set_status(node,ec_visited);
    }
}

void db_node_ec_action_set_status_special_visited(dBNodeEc * node){
  db_node_ec_set_status(node, ec_special_visited);
}

boolean db_node_ec_check_status_special(dBNodeEc* node)
{
  if (  (db_node_ec_check_status(node, ec_special_none)==true)
	||
	(db_node_ec_check_status(node, ec_special_visited)==true)
	||
	(db_node_ec_check_status(node, ec_special_pruned)==true)
	)
    {
      return true;
    }
  else
    {
      return false;
    }
}

void db_node_ec_action_specialise_status(dBNodeEc * node){
  if (db_node_ec_check_status(node, ec_visited)==true)
    {
      db_node_ec_set_status(node, ec_special_visited);      
    }
  else if (db_node_ec_check_status(node, ec_none)==true)
    {
      db_node_ec_set_status(node, ec_special_none);      
    }
  else if (db_node_ec_check_status(node, ec_pruned)==true)
    {
      db_node_ec_set_status(node, ec_special_pruned);      
    }
  else if ( (db_node_ec_check_status(node, ec_read_start_forward)==true)
	    ||
	    (db_node_ec_check_status(node, ec_read_start_reverse)==true)
	    ||
	    (db_node_ec_check_status(node, ec_read_start_forward_and_reverse)==true)
	    )
  {
    // happy to lose PCR dup info at this stage
    db_node_ec_set_status(node, ec_special_none);
    warn("Warn Zam (zam@well.ox.ac.uk) that you met a PCR dup status during "
         "genotyping. He knows how to fix it\n"); 
  }
  else if (db_node_ec_check_status_special(node)==false)
  {
    //NodeStatusEc ret = node->status;
    //warn("Warn Zam (zam@well.ox.ac.uk) that you met a status of %d. \n"
    //   "Could signal a subtle (but with your information, fixable) bug.\n",
    //   (int)ret); 
  }
}



void db_node_ec_action_unspecialise_status(dBNodeEc * node){
  if (db_node_ec_check_status(node, ec_special_visited)==true)
    {
      db_node_ec_set_status(node, ec_visited);      
    }
  else if (db_node_ec_check_status(node, ec_special_none)==true)
    {
      db_node_ec_set_status(node, ec_none);      
    }
  else if (db_node_ec_check_status(node, ec_special_pruned)==true)
    {
      db_node_ec_set_status(node, ec_pruned);      
    }
  
}


void db_node_ec_action_set_status_ignore_this_node(dBNodeEc * node)
{
  db_node_ec_set_status(node, ec_ignore_this_node);
}

void db_node_ec_action_set_status_visited_or_visited_and_exists_in_reference(dBNodeEc * node){

  if (db_node_ec_check_status(node, ec_exists_in_reference))
    {
      db_node_ec_set_status(node, ec_visited_and_exists_in_reference);      
    }
  else if (!db_node_ec_check_status(node, ec_unassigned)) 
    {
      db_node_ec_set_status(node, ec_visited);
    }

}


void db_node_ec_action_unset_status_visited_or_visited_and_exists_in_reference(dBNodeEc * node){
  if (db_node_ec_check_status_visited_and_exists_in_reference(node)==true)
  {
    db_node_ec_set_status(node, ec_exists_in_reference);
  }
  else if (db_node_ec_check_status(node, ec_visited)==true)
  {
    db_node_ec_set_status(node, ec_none);
  }
      
}

void db_node_ec_action_unset_status_visited_or_visited_and_exists_in_reference_or_ignore_this_node(dBNodeEc * node){
  if (db_node_ec_check_status_visited_and_exists_in_reference(node))
  {
    db_node_ec_set_status(node, ec_exists_in_reference);
  }
  else if (db_node_ec_check_status(node, ec_visited))
  {
    db_node_ec_set_status(node, ec_none);
  }
  else if (db_node_ec_check_status(node, ec_ignore_this_node))
  {
    db_node_ec_set_status(node, ec_none);
  }
      
}




void db_node_ec_action_do_nothing(dBNodeEc * node)
{
  // Let the compiler know that we're deliberately not using a paramater
  (void)node;
}


boolean db_node_ec_check_status_none(dBNodeEc * node){
  return db_node_ec_check_status(node,ec_none);
}


// DEV: replace calls to db_node_ec_check_for_flag_ALL_OFF
//      with db_node_ec_check_status(node, ec_unassigned)
//wrapper - in future this will be part of the flags introduced by Mario et al
boolean db_node_ec_check_for_flag_ALL_OFF(dBNodeEc * node) {

  if ( db_node_ec_check_status(node, ec_unassigned)==true)
    {
      return true;
    }
  else
    {
      return false;
    }
}


boolean db_node_ec_check_status_visited(dBNodeEc * node){
  return db_node_ec_check_status(node, ec_visited);
}

boolean db_node_ec_check_status_exists_in_reference(dBNodeEc * node){
  return db_node_ec_check_status(node, ec_exists_in_reference);
}

boolean db_node_ec_check_status_visited_and_exists_in_reference(dBNodeEc * node){
  return db_node_ec_check_status(node, ec_visited_and_exists_in_reference);
}

boolean db_node_ec_check_status_is_not_exists_in_reference(dBNodeEc * node){

  if (db_node_ec_check_status(node, ec_exists_in_reference)==true )
    {
      return false;
    }
  else
    {
      return true;
    }
}

boolean db_node_ec_check_status_is_not_visited(dBNodeEc* node)
{
   if (db_node_ec_check_status(node, ec_visited)  )
    {
      return false;
    }
  else
    {
      return true;
    }
  
}

boolean db_node_ec_check_status_is_not_visited_or_visited_and_exists_in_reference(dBNodeEc * node){
  if (db_node_ec_check_status(node, ec_visited_and_exists_in_reference) || db_node_ec_check_status(node, ec_visited)  )
    {
      return false;
    }
  else
    {
      return true;
    }
}


boolean db_node_ec_condition_is_not_marked_to_be_ignored(dBNodeEc* node)
{
  return !db_node_ec_check_status(node, ec_ignore_this_node);
}

boolean db_node_ec_condition_always_true(dBNodeEc* node)
{
  // Let the compiler know we're deliberately ignoring our parameters
  (void)node;
  return true;
}



// wrapper to allow comparison with read_start_forward OR read_start_forward_and_reverse
// ie. if you want to check if this node is the start of a read (forwards), then call this, with second argument forward
// and it will handle the case when it is a read start in both directions

boolean db_node_ec_check_read_start(dBNodeEc* node, Orientation ori)
{
  if ( (ori==forward) && 
       (db_node_ec_check_status(node, ec_read_start_forward) || db_node_ec_check_status(node, ec_read_start_forward_and_reverse) ) ) 
    {
      return true;
    }
  else if ( (ori==reverse) && 
       (db_node_ec_check_status(node, ec_read_start_reverse) || db_node_ec_check_status(node, ec_read_start_forward_and_reverse) ) ) 
    {
      return true;
    }
  else
    {
      return false;
    }
}

void db_node_ec_set_read_start_status(dBNodeEc* node, Orientation ori)
{
  if (db_node_ec_check_status(node, ec_visited) )
    {
      warn("setting status of a visisted node to read_start\n");
    }
  else if (db_node_ec_check_status(node, ec_unassigned) )
    {
      die("Setting status of an ec_unassigned node to read_start. Exit");
    }
  else if (db_node_ec_check_status(node, ec_read_start_forward) && (ori==reverse) )
    {
      db_node_ec_set_status(node, ec_read_start_forward_and_reverse);
    }
  else if (db_node_ec_check_status(node, ec_read_start_reverse) && (ori==forward) )
    {
      db_node_ec_set_status(node, ec_read_start_forward_and_reverse);
    }
  else if (ori==forward)
    {
      db_node_ec_set_status(node, ec_read_start_forward);
    }
  else if (ori==reverse)
    {
      db_node_ec_set_status(node, ec_read_start_reverse);
    }
  else
    {
      die("Programming error in read start setting of node");
    }
  
}


