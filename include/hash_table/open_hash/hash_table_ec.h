/*
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
  hash_table_ec.h 
  
  open hash table implementation - ie every bucket has a predifined size 
  overloads results in rehashing
  all the routines as prefixed with hash_table
*/

#ifndef OPEN_HASH_TABLE_EC_H_
#define OPEN_HASH_TABLE_EC_H_

#include "global.h"
#include "element_ec.h"

typedef struct
{
  short kmer_size;
  long long number_buckets;
  int bucket_size;
  ElementEc * table; 
  short * next_element; //keeps index of the next free element in bucket 
  long long * collisions;
  long long unique_kmers;
  int max_rehash_tries;
} HashTableEc;


HashTableEc * hash_table_ec_new(int number_bits, int bucket_size, int max_rehash_tries, short kmer_size);

void hash_table_ec_free(HashTableEc * * hash_table);

//if the key is present applies f otherwise adds a new element for kmer
boolean hash_table_ec_apply_or_insert(Key key, void (*f)(ElementEc*), HashTableEc *);

//applies f to every element of the table
void hash_table_ec_traverse(void (*f)(ElementEc *),HashTableEc *);
long long hash_table_ec_traverse_returning_sum(long long (*f)(ElementEc *),HashTableEc * hash_table);
void hash_table_ec_traverse_passing_int(void (*f)(ElementEc *, int*),HashTableEc * hash_table, int* num);
void hash_table_ec_traverse_passing_ints_and_path(
  void (*f)(ElementEc *, int*, int*, dBNodeEc**, Orientation*, Nucleotide*, char*, int),
  HashTableEc * hash_table, int* num1, int* num2, 
  dBNodeEc** p_n, Orientation* p_o, Nucleotide* p_lab, char* p_str, int len);

void hash_table_ec_traverse_passing_3ints_and_path(
  void (*f)(ElementEc *, int*, int*, int*, dBNodeEc**, Orientation*, Nucleotide*, char*, int),
  HashTableEc * hash_table, int* num1, int* num2, int* num3, 
  dBNodeEc** p_n, Orientation* p_o, Nucleotide* p_lab, char* p_str, int len);


//if the element is not in table create an element with key and adds it
ElementEc * hash_table_ec_find_or_insert(Key key, boolean * found, HashTableEc * hash_table);
ElementEc * hash_table_ec_insert(Key key, HashTableEc * hash_table);

void hash_table_ec_print_stats(HashTableEc *);

long long hash_table_ec_get_unique_kmers(HashTableEc *);

//return entry for kmer
ElementEc * hash_table_ec_find(Key key, HashTableEc * hash_table);

long long hash_table_ec_get_capacity(HashTableEc * hash_table);

#endif /* OPEN_HASH_TABLE_EC_H_ */
