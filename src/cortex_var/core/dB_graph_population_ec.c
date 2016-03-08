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
  dB_graph_population.c - implementation
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <inttypes.h>

// cortex_var headers
#include "element_ec.h"
#include "open_hash/hash_table_ec.h"
#include "db_graph_ec.h"
#include "dB_graph_population_ec.h"
#include "seq.h"
#include "file_reader.h"
#include "model_selection.h"
#include "maths.h"
#include "db_variants.h"
#include "db_complex_genotyping.h"




//This function does not  check that it there is such an edge in the specified person/colour - but it does check if the target node is in the specific person.
//if you want to be sure the dge exists in that colour, then check it before calling this function
//The last argument allows you to apply the operation to some subgraph - eg you might take the unuiion of colours 2 and 3, or of all colours.
//If for example you wanted to get the next node in the graph irrespective of colour, ec_get_colour would return the union (bitwise AND) of all edges in a node.
dBNodeEc * db_graph_ec_get_next_node_in_subgraph_defined_by_func_of_colours(dBNodeEc * current_node, Orientation current_orientation, 
								       Orientation * next_orientation,
								       Nucleotide edge, Nucleotide * reverse_edge,dBGraphEc * db_graph, 
								       Edges (*ec_get_colour)(const dBNodeEc*)
								       )
{
  BinaryKmer local_copy_of_kmer;
  binary_kmer_assignment_operator(local_copy_of_kmer, current_node->kmer);
  
  BinaryKmer tmp_kmer;
  dBNodeEc * next_node=NULL;
  
  // after the following line tmp_kmer and rev_kmer are pointing to the same B Kmer
  BinaryKmer* rev_kmer = binary_kmer_reverse_complement(&local_copy_of_kmer,db_graph->kmer_size, &tmp_kmer);
  
  if (current_orientation == reverse){   
    *reverse_edge = binary_kmer_get_last_nucleotide(&local_copy_of_kmer);
    binary_kmer_assignment_operator(local_copy_of_kmer,*rev_kmer);
  }
  else{
    *reverse_edge = binary_kmer_get_last_nucleotide(rev_kmer);
  }

  binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&local_copy_of_kmer, edge, db_graph->kmer_size);

   //get node from table
  next_node = hash_table_ec_find(element_ec_get_key(&local_copy_of_kmer,db_graph->kmer_size, &tmp_kmer),db_graph);

  if (next_node != NULL)
    {
      *next_orientation = db_node_ec_get_orientation(&local_copy_of_kmer,next_node,db_graph->kmer_size);
    }
  else
    {
      //no else
    }

  //need to check the node is in the specified subgraph graph
  if (! (db_node_ec_is_this_node_in_subgraph_defined_by_func_of_colours(next_node, ec_get_colour)) )
    {
      return NULL;
    }

  return next_node;
  

}
