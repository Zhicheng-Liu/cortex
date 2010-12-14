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
   binary_kmer.c - routines to manipulate binary kmers
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <event_encoding.h>
#include <global.h>
#include <string.h>

//returns Undefined if given non AGCT character
Nucleotide char_to_binary_nucleotide(char c)
{
	switch (c)
	{
	case '0':
	  return Zero;
	case '1':
	  return One;
	case '2':
	  return Two;
	case '3':
	  return Three;
	default:
	  return Undefined;
	}
}




//this method does basically nothing
Nucleotide reverse_binary_nucleotide(Nucleotide n)
{
  switch (n)
    {
    case Zero:
      return Zero;
    case One:
      return One;
    case Two:
      return Two;
    case Three:
      return Three;
    default:
      printf("Calling reverse_binary_nucleotide on non-existent nucleotide %i\n",n);
      exit(1);
    }
}

char binary_nucleotide_to_char(Nucleotide n)
{
	switch (n) {
	case Zero:
	  return '0';
	case One:
	  return '1';
	case Two:
	  return '2';
	case Three:
	  return '3';
	default:
	  printf("Non existent binary nucleotide %d\n",n);
	  assert(0); 
	  return 'N'; //Don't really return this, must fail before this point. But stops compiler warning.
	}
}


boolean good_symbol(char c){
  boolean ret;
  if (c  != '0' && c != '1' && 
      c != '2' && c != '3' && 
      c != 'N' && c != 'n' 
      ){
    ret = false;
  }	
  else{
    ret =  true;
  }
  
  return ret;
}
