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
  dB_graph_ec.h 

  all the routines as prefixed with db_graph
*/

#ifndef DB_GRAPH_EC_H_
#define DB_GRAPH_EC_H_

#include <stdio.h>
  
#include "global.h"
#include "open_hash/hash_table_ec.h"

typedef HashTableEc dBGraphEc;

//pays no attention to whether there is an edge joining current_node to the node you would get by adding this nucleotide.
//just checksto see if such a node is in the graph
dBNodeEc * db_graph_ec_get_next_node(dBNodeEc * current_node, Orientation current_orientation, 
				Orientation * next_orientation,
				Nucleotide edge, Nucleotide * reverse_edge,dBGraphEc * db_graph);




//Functions applying to whole graph

int db_graph_ec_clip_tip(dBNodeEc * node, int limit,dBGraphEc * db_graph);

void db_graph_ec_set_all_visited_nodes_to_status_none(dBGraphEc* hash_table);


#endif /* DB_GRAPH_EC_H_ */
