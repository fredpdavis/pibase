/* altloc_filter.c

Purpose: filter PDB file to leave only highest occupancy, or first instance
         of atoms with ALTLOC flags set

C implementation of altloc_filter.pl - need to speed up to run on the fly


Fred P. Davis, HHMI-JFRC (davisf@janelia.hhmi.org)

Copyright 2005,2008 Fred P. Davis.
See the file COPYING for copying permission.

This file is part of PIBASE.

PIBASE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PIBASE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PIBASE.  If not, see <http://www.gnu.org/licenses/>.



*/


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>


//#define DEBUG 1
#define MAXLINELENGTH 80

#define INITNUMATOMS 100
#define ATOMBLOCKSIZE 100

#define INITSETSIZE 100
#define SETBLOCKSIZE 100


#define INITNUMLINES 100
#define LINEBLOCKSIZE 100



//STRUCTURES

enum boolean { NO, YES } ;
enum rectype_t { NON, ATOM, HETATM } ;

struct atom_Struct {
   char		rectype[7] ;
   int          atomno ;
   char         atomna[5] ;
   char         altloc[2] ;
   char         resna[4] ;
   char         chainid[2] ;
   signed int   resno ;
   char         inscode[2] ;

   float        coord[3] ;

   float        occup ;
} ;
typedef struct atom_Struct atom_t ;


struct atomsig_Struct {
   char		atomsig[50] ;
   int		linenum ;
} ;
typedef struct atomsig_Struct atomsig_t ;


struct pdbline_Struct {
   int		linenum ;
   char		line[81] ;
   boolean	print_fl ;
   rectype_t	rectype
} ;
typedef struct pdbline_Struct pdbline_t ;


struct readinlines_Struct {
   int		num_lines ;
   pdbline_t	*details ;
} ;
typedef struct readinlines_Struct readinlines_t ;


//FUNCTION DECLARATION
#define Error( Str )   fprintf( stderr, "%s\n", Str ), exit( 1 )

readinatoms_t *readinatoms( void ) ;



int main(int argc, char *argv[])
{
   readinatoms_t *atoms ;
   char *pdb_fn ;
   FILE *pdb_fp ;

   if (argc > 1 ) {
      pdb_fn = argv[1] ;
   } else {
      Error("usage: altloc_filter pdbfile");
   }


   pdb_fp = fopen(pdb_fn, "r") ;
   if (pdb_fp == NULL) {
      fprintf(stderr, "ERROR: PDB file %s does not exist\n", pdb_fn) ;
      exit(1) ;
   }

//   0. read the lines into an array
   pdblines = readinlines(pdb_fp) ;
   fclose(pdb_fp) ;


//   1. iterate over the lines: store sigs in a binary lookup tree together
//       with the line numbers that contain it
   atomsigtree = build_sigtree(pdblines) ;

//   2. then iterate over the leaf nodes of the tree to find duplicate entries;
//        determine which one to print; set the print flag
   highestoccup( atomsigtree, pdblines) ;


//   3. iterate over the lines and print those with the flags set
   displaylines( pdblines ) ;

#ifdef DEBUG
   fprintf(stderr, "read %d atoms\n", atoms->number) ;
#endif

   atoms_kdtree = call_build_kdtree(atoms) ;

   return 0;
}



readinlines_t *readinlines( FILE *pdb_fp )
{
   int i ;
   char line[MAXLINELENGTH] ;
   readinlines_t *results ;

   int linelistsize = INITNUMLINES ;

   results = malloc(sizeof(readlines_t)) ;
   if (results == NULL ) {
      Error("readinlines(): out of memory on malloc()") ;
   }

   results->details = malloc(linelistsize * sizeof(pdbline_t)) ;
   if (results->details == NULL) {
      Error("readinlnies(): out of memory on malloc()") ;
   }

   i = 0 ;
   while ((fgets(line, sizeof(line), pdb_fp)) && j ) {
      *(line+(strlen(line)-1)) = '\0';

      if ( i >= linelistsize ) {
         line_t *newp ;
	 linelistsize += LINEBLOCKSIZE ;
	 newp = realloc(results->details, linelistsize * sizeof(pdbline_t)) ;

	 if (newp == NULL ) {
	    Error("readinlines(): out of memory on realloc()\n") ;
	 }
      }

      strcpy(results->details[i].line, &line) ;
      results->details[i].print_fl = 1 ;

      if ( (line[0] == 'A') &&
           (line[1] == 'T') &&
	   (line[2] == 'O') &&
	   (line[3] == 'M') ) {
         results->details[i].rectype = ATOM ;
      } else if ( (line[0] == 'A') &&
                  (line[1] == 'T') &&
                  (line[2] == 'O') &&
                  (line[3] == 'M') ) {
         results->details[i].rectype = HETATM ;
      } else {
         results->details[i].rectype = NON ;
      }

      i++ ;
   }

   results->number = i ;

   return results ;
}


/* readinatoms: reads in ATOM records from STDIN and returns a pointer to a
   readinatoms)_t struct with a pointer to an array of atom_t structs,
   the number of atoms, and the coordinate minimum and maximum along each
   dimension*/
sigtree_t *build_sigtree (readinlines_t pdblines)
{
   char line[MAXLINELENGTH] ;
   char tempsubstr[MAXLINELENGTH] ;

   sigtree_t *t ;

   int i, j;
   int atomlistsize = INITNUMATOMS ;


   results = malloc(sizeof(sigtree_t) k) ;
   if (result == NULL) {
      Error("Out of memory on result malloc()\n") ;
   }

   results->details = malloc(atomlistsize * sizeof(atom_t)) ;
   if (results->details == NULL) {
      Error("Out of memory on details malloc()\n") ;
   }


   i = 0;
   j = 1;
   while ( i <= pdblines->num_lines) {

      if ( (pdblines->details[i].rectype == ATOM) ||
           (pdblines->details[i].rectype == HETATM) ) {
         
         strncpy(tempsubstr, (line + 6), 5) ;
         *(tempsubstr+(5)) = '\0';
         result->details[i].atomno = atoi( tempsubstr ) ;
      
         strncpy( result->details[i].atomna, (line + 12), 4 ) ;
         result->details[i].atomna[4] = '\0' ;
      
         strncpy(result->details[i].altloc, (line + 16), 1) ;
         result->details[i].altloc[1] = '\0' ;
      
         strncpy(result->details[i].resna, (line + 17), 3) ;
         result->details[i].resna[3] = '\0' ;
      
         strncpy(result->details[i].chainid, (line + 21), 1) ;
         result->details[i].chainid[1] = '\0' ;
      
         strncpy(tempsubstr, (line + 22), 4) ;
         *(tempsubstr+(4)) = '\0' ;
         result->details[i].resno = atoi( tempsubstr ) ;
      
         strncpy(result->details[i].inscode, (line + 26), 1) ;
         result->details[i].inscode[1] = '\0' ;
      
         strncpy(tempsubstr, (line+30), 8) ;
         *(tempsubstr+(8)) = '\0';
         result->details[i].coord[0] = atof( tempsubstr ) ;
      
         strncpy(tempsubstr, (line+38), 8) ;
         *(tempsubstr+(8)) = '\0';
         result->details[i].coord[1] = atof( tempsubstr ) ;
      
         strncpy(tempsubstr, (line+46), 8) ;
         *(tempsubstr+(8)) = '\0';
         result->details[i].coord[2] = atof( tempsubstr ) ;
      
         strncpy(tempsubstr, (line+54), 6) ;
         *(tempsubstr+(6)) = '\0';
         result->details[i].occup = atof( tempsubstr ) ;


         if (i == 0 ) {
            result->coord_min[0] = result->details[i].coord[0] ;
            result->coord_min[1] = result->details[i].coord[1] ;
            result->coord_min[2] = result->details[i].coord[2] ;
            result->coord_max[0] = result->details[i].coord[0] ;
            result->coord_max[1] = result->details[i].coord[1] ;
            result->coord_max[2] = result->details[i].coord[2] ;

         } else {

            if (result->coord_min[0] > result->details[i].coord[0]) {
               result->coord_min[0] = result->details[i].coord[0] ; }
            if (result->coord_min[1] > result->details[i].coord[1]) {
               result->coord_min[1] = result->details[i].coord[1] ; }
            if (result->coord_min[2] > result->details[i].coord[2]) {
               result->coord_min[2] = result->details[i].coord[2] ; }

            if (result->coord_max[0] < result->details[i].coord[0]) {
               result->coord_max[0] = result->details[i].coord[0] ; }
            if (result->coord_max[1] < result->details[i].coord[1]) {
               result->coord_max[1] = result->details[i].coord[1] ; }
            if (result->coord_max[2] < result->details[i].coord[2]) {
               result->coord_max[2] = result->details[i].coord[2] ; }
         }

      }

      i++ ;
   }

   result->number = i ;

   return result ;
}


/* eucliddist: determines the euclidean distance between two 3D points */
float eucliddist (float a[3], float b[3])
{
   float distance ;

   distance = sqrt( (a[0] - b[0]) * (a[0] - b[0]) +
                    (a[1] - b[1]) * (a[1] - b[1]) +
                    (a[2] - b[2]) * (a[2] - b[2]) ) ;

   return distance ;
}
