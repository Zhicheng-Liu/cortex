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
  element_ec.h  - defines the interface for the de Bruijn graph node. The
  implementation is complemented by a hash table that stores every node indexed
  by kmers (BinaryKmers). 

  The element routines, ie the one required by hash_table/priority queue, are
  prefixed with element_

  The de Bruijn based routines are prefixed with db_node_ec_
*/

#ifndef ELEMENT_EC_H_
#define ELEMENT_EC_H_

#include <stdio.h>
#include <inttypes.h>

#include "global.h"
#include "binary_kmer.h"

// type definitions

typedef char Edges;

typedef uint32_t Covg;
//typedef int Covg
// COVG_MAX is defined as UINT_MAX in element.c
const Covg COVG_MAX;

typedef struct __attribute__((__packed__)) {
  BinaryKmer kmer;
  //  Covg       coverage[NUMBER_OF_COLOURS];
  Edges      individual_edges[NUMBER_OF_COLOURS];
  char       status; // will cast a NodeStatus to char
} ElementEc;


typedef enum
  {
    ec_unassigned   = 0,
    ec_none         = 1,
    ec_visited      = 2,
    ec_pruned       = 3,
    ec_exists_in_reference = 4,
    ec_visited_and_exists_in_reference = 5,
    ec_to_be_dumped = 6, //to be dumped as binary 
    ec_read_start_forward = 7,//used when removing duplicate reads
    ec_read_start_reverse = 8,//used when removing duplicate reads
    ec_read_start_forward_and_reverse = 9,//used when removing duplicate reads  
    ec_ignore_this_node = 10,
    ec_in_desired_genotype = 11,
    ec_special_visited = 12,
    ec_special_none = 13,
    ec_special_pruned = 14,
    ec_fw_strand = 15,
    ec_rv_strand = 16,
  } NodeStatusEc;

typedef ElementEc dBNodeEc;
typedef BinaryKmer* Key;


ElementEc* new_element_c();
void free_element_ec(ElementEc** element);
void element_ec_assign(ElementEc* e1, ElementEc* e2);


// utility function for getting the desired edge char,
// by specifying the appropriate colour

// gets copy of edge
Edges element_ec_get_edge_copy(const ElementEc e, int colour);
Edges element_ec_get_union_of_edges(ElementEc e);
Edges element_ec_get_colour_union_of_all_colours(const ElementEc*);

Edges element_ec_get_colour0(const ElementEc* e);
Edges element_ec_get_colour1(const ElementEc* e);
Edges element_ec_get_last_colour(const ElementEc* e);

Covg element_ec_get_covg_union_of_all_covgs(const dBNodeEc*);
Covg element_ec_get_covg_for_colourlist(const dBNodeEc* e, int* colour_list,
                                     int list_len);

Covg element_ec_get_covg_colour0(const dBNodeEc* e);

#if NUMBER_OF_COLOURS > 1
Covg element_ec_get_covg_colour1(const dBNodeEc* e);
#endif

Covg element_ec_get_covg_last_colour(const dBNodeEc* e);

void element_ec_add_edges(ElementEc* e, int colour, Edges edge_char);
void element_ec_set_edges(ElementEc* e, int colour, Edges edge_char);
void element_ec_reset_one_edge(ElementEc* e, Orientation orientation, Nucleotide nucleotide, int colour);

int element_ec_get_number_of_people_or_pops_containing_this_element(ElementEc* e);


boolean element_ec_smaller(ElementEc e1, ElementEc e2);
BinaryKmer* element_ec_get_kmer(ElementEc *e);
boolean element_ec_is_key(Key key, ElementEc e);
Key element_ec_get_key(BinaryKmer*,short kmer_size, Key preallocated_key);
void element_ec_initialise(ElementEc *,Key, short kmer_size);
void element_ec_initialise_kmer_covgs_edges_and_status_to_zero(ElementEc * e);

void element_ec_set_kmer(ElementEc *e, Key kmer, short kmer_size);


// reverse orientation
Orientation element_ec_opposite_orientation(Orientation);
Orientation db_node_ec_get_orientation(BinaryKmer*, dBNodeEc *, short kmer_size);

// add an edge between nodes -- NB: it adds both edges: forward and reverse
boolean db_node_ec_add_edge(dBNodeEc *, dBNodeEc *, Orientation, Orientation, short kmer_size, int colour); 


// returns yes if the label defined by the nucleotide coresponds to an 
// outgoing edge in the side defined by the orientation.   
boolean db_node_ec_edge_exist(dBNodeEc *e, Nucleotide n, Orientation o, int colour);

// final argument f is a function that returns an Edge that is a function of the
// different colured edges in a node. e.g. we might be interested in the union
// of all the coloured edges in the graph, or just the colour/edge for the first
// person, or the union of all edges except that corresponding to the reference.
boolean db_node_ec_edge_exist_within_specified_function_of_coloured_edges(
  dBNodeEc *element, Nucleotide base,
  Orientation orientation, Edges (*f)(const ElementEc* ));


//returns the label of the first outgoing edge -- leaving from the side 
//defined by orientation. 
boolean db_node_ec_has_precisely_one_edge(dBNodeEc *e, Orientation o, Nucleotide *n, int colour);


boolean db_node_ec_has_precisely_one_edge_in_subgraph_defined_by_func_of_colours(
  dBNodeEc * node, Orientation orientation, Nucleotide * nucleotide, 
	Edges (*get_colour)(const dBNodeEc*) );

boolean db_node_ec_has_precisely_one_edge_in_union_graph_over_all_people(
  dBNodeEc * node, Orientation orientation, Nucleotide * nucleotide);

boolean db_node_ec_has_precisely_two_edges(dBNodeEc * node, Orientation orientation,
                                        Nucleotide *n1, Nucleotide *n2, int colour);

void db_node_ec_reset_all_edges_for_all_people_and_pops_to_zero(ElementEc* e);

//forgets about the edges
void db_node_ec_reset_edges(dBNodeEc *e, int colour);

void db_node_ec_reset_edge(dBNodeEc *e, Orientation o, Nucleotide n, int colour);



//TODO - maybe do not need to export this:
void db_node_ec_reset_specified_edges(dBNodeEc * node, Orientation orientation,
                                   Nucleotide nucleotide,
                                   void (*f)(dBNodeEc*, Orientation, Nucleotide));


//check that the edges are 0's
boolean db_node_ec_edges_reset(dBNodeEc *node, int colour);

boolean db_node_ec_check_status(dBNodeEc *node, NodeStatusEc status);
boolean db_node_ec_check_status_not_pruned(dBNodeEc *node);
boolean db_node_ec_check_status_not_pruned_or_visited(dBNodeEc *node);
boolean db_node_ec_check_status_to_be_dumped(dBNodeEc *node);
boolean db_node_ec_check_status_is_not_visited(dBNodeEc *node);


void db_node_ec_set_status(dBNodeEc *node,NodeStatusEc status);
void db_node_ec_trio_aware_set_pruned_status(dBNodeEc *node, int colour);
void db_node_ec_set_status_to_none(dBNodeEc *node);



//actions and conditions 

void db_node_ec_action_set_status_none(dBNodeEc *node);
void db_node_ec_action_set_status_of_unpruned_to_none(dBNodeEc *node);

void db_node_ec_action_set_status_pruned(dBNodeEc *node);
void db_node_ec_action_set_status_visited(dBNodeEc *node);
void db_node_ec_action_set_status_visited_unless_marked_to_be_ignored(dBNodeEc * node);
void db_node_ec_action_set_status_special_visited(dBNodeEc *node);
boolean db_node_ec_check_status_special(dBNodeEc*node);
void db_node_ec_action_specialise_status(dBNodeEc *node);
void db_node_ec_action_unspecialise_status(dBNodeEc *node);

void db_node_ec_action_set_status_ignore_this_node(dBNodeEc *node);

void db_node_ec_action_set_status_visited_or_visited_and_exists_in_reference(dBNodeEc *node);

void db_node_ec_action_unset_status_visited_or_visited_and_exists_in_reference(dBNodeEc *node);

void db_node_ec_action_unset_status_visited_or_visited_and_exists_in_reference_or_ignore_this_node(dBNodeEc * node);


void db_node_ec_action_do_nothing(dBNodeEc * node);

boolean db_node_ec_check_status_none(dBNodeEc * node);
boolean db_node_ec_check_for_flag_ALL_OFF(dBNodeEc * node);


boolean db_node_ec_check_status_visited(dBNodeEc * node);

boolean db_node_ec_check_status_exists_in_reference(dBNodeEc * node);

boolean db_node_ec_check_status_visited_and_exists_in_reference(dBNodeEc * node);

boolean db_node_ec_check_status_is_not_exists_in_reference(dBNodeEc * node);

boolean db_node_ec_check_status_is_not_visited_or_visited_and_exists_in_reference(dBNodeEc * node);

boolean db_node_ec_condition_is_not_marked_to_be_ignored(dBNodeEc* node);

boolean db_node_ec_condition_always_true(dBNodeEc* node);





void db_node_ec_increment_coverage(dBNodeEc* e, int colour);
void db_node_ec_update_coverage(dBNodeEc* e, int colour, long update);
Covg db_node_ec_get_coverage_tolerate_null(const dBNodeEc* e, int colour);
Covg db_node_ec_get_coverage(const dBNodeEc* e, int colour);
void db_node_ec_set_coverage(dBNodeEc* e, int colour, Covg covg);

Covg db_node_ec_get_coverage_in_subgraph_defined_by_func_of_colours(
  const dBNodeEc* e, Covg (*get_covg)(const dBNodeEc*));



//check if node doesn't have any edges in a given orientation
boolean db_node_ec_is_blunt_end(dBNodeEc * node, Orientation orientation, int colour);

boolean db_node_ec_is_blunt_end_in_subgraph_given_by_func_of_colours(
  dBNodeEc * node, Orientation orientation,  Edges (*get_colour)(const dBNodeEc*));


boolean db_node_ec_is_this_node_in_this_person_or_populations_graph(dBNodeEc* node, int colour);


boolean db_node_ec_is_this_node_in_subgraph_defined_by_func_of_colours(
  dBNodeEc* node, Edges (*get_colour)(const dBNodeEc*));


//functions for binary format
void db_node_ec_print_multicolour_binary(FILE *fp, dBNodeEc *node);

void db_node_ec_print_single_colour_binary_of_colour0(FILE *fp, dBNodeEc *node);
void db_node_ec_print_single_colour_binary_of_specified_colour(FILE *fp, dBNodeEc *node, int colour);

//reading multicolour binaries
boolean db_node_ec_read_multicolour_binary(FILE *fp, short kmer_size, dBNodeEc *node,
                                        int num_colours_in_binary,
                                        int binversion_in_binheader);


// read a binary for an individual person, as dumped by the target "graph"
// load this data into given colour
boolean db_node_ec_read_single_colour_binary(FILE *fp, short kmer_size, dBNodeEc *node,
                                          int colour, int binversion_in_binheader);


boolean db_node_ec_check_read_start(dBNodeEc*node, Orientation ori);

void db_node_ec_set_read_start_status(dBNodeEc* node, Orientation ori);


#endif /* ELEMENT_EC_H_ */
