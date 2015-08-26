/* kdcontacts.c - Calculates interatomic distances of PDB ATOMs

Description: Lists all atom pairs where distance is below a given threshold.
Uses a three dimensional kd-tree to efficiently perform fixed radius queries
of PDB coordinates.

Usage: ./kdcontacts [sphere radius] < pdbfile
Sphere radius defaults to 5 Angstroms

NOTE: - only uses ^ATOM records



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



//STRUCTURES

struct atom_Struct {
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



struct readinatoms_Struct {
   int          number ;
   float        coord_max[3] ;
   float        coord_min[3] ;
   atom_t       *details ;
} ;
typedef struct readinatoms_Struct readinatoms_t ;


typedef enum {CHILD_LEFT, CHILD_RIGHT} lorr_t ;


struct kdtree_Struct {
   struct kdtree_Struct *parent;

   int          depth ; // will store point number if leaf

   int          *points ;
   int          numpoints ;

   float        splitval ;

   float        bound_min[3] ;  // stores bounds from ancestors
   float        bound_max[3] ;  // stores bounds from ancestors
   lorr_t       lorr ;  //is this node the LEFT or RIGHT child of the parent

   struct kdtree_Struct *left ;
   struct kdtree_Struct *right ;
} ;
typedef struct kdtree_Struct kdtree_t ;



struct splitpoints_Struct {
   int          *Lset ;
   int          *Rset ;
   int          Lsize ;
   int          Rsize ;
   float        splitval ;
} ;
typedef struct splitpoints_Struct splitpoints_t ;


typedef enum { SAME, LESS, MORE } point_point_t ;
typedef enum { EMPTY, INTERSECT, CONTAINED } range_range_t ;


//FUNCTION DECLARATION
#define Error( Str )   fprintf( stderr, "%s\n", Str ), exit( 1 )

readinatoms_t *readinatoms( void ) ;

kdtree_t *call_build_kdtree( readinatoms_t *atoms ) ;

kdtree_t *build_kdtree(kdtree_t *parent, int depth, int *points, int numpoints, readinatoms_t *atoms, lorr_t lorr ) ;

splitpoints_t split_points (int *points, int numpoints, int depth, readinatoms_t *atoms) ;

float median_quick_select(float arr[], int n ) ;

void display_contacts (kdtree_t *atoms_kdtree, readinatoms_t *atoms, float *radius) ;

void atom_contacts (int centerind, float *radius, kdtree_t *atoms_kdtree, readinatoms_t *atoms) ;

void search_kdtree (kdtree_t *cur_kdtree, float *rect_min, float *rect_max, int centerind, float *center, float *radius, readinatoms_t *atoms) ;

range_range_t query_vs_node( kdtree_t *cur_kdtree, float *rect_min, float *rect_max ) ;

void report_kdtree( kdtree_t *cur_kdtree, int centerind, float *center, float *radius, readinatoms_t *atoms) ;

float eucliddist( float a[3], float b[3]) ;




int main(int argc, char *argv[])
{
   readinatoms_t *atoms ;
   kdtree_t *atoms_kdtree ;
   float radius ;

   if (argc > 1 ) {
      radius = atof(argv[1]) ;
   } else {
      radius = 5.0 ;
   }

   atoms = readinatoms() ;
#ifdef DEBUG
   fprintf(stderr, "read %d atoms\n", atoms->number) ;
#endif

   atoms_kdtree = call_build_kdtree(atoms) ;

   display_contacts(atoms_kdtree, atoms, &radius) ;

   return 0;
}




/* call_build_kdtree: setups atom data for a call to build_kdtree */
kdtree_t *call_build_kdtree(readinatoms_t *atoms)
{
   int *points ;
   kdtree_t *atoms_kdtree ;
   kdtree_t *temp_parent = NULL ;
   int j ;

   points = malloc(atoms->number * sizeof(int)) ;
   if (points == NULL) {
      Error("Out of memory on points malloc()\n") ; }

   for (j = 0; j < atoms->number; j++ ) {
      points[j] = j ; }

   atoms_kdtree = build_kdtree(temp_parent, 0, points, atoms->number, atoms, CHILD_LEFT);

   return atoms_kdtree ;
}



/* build_kdtree: recursively builds a kd-tree */
kdtree_t *build_kdtree(kdtree_t *parent, int depth, int *points, int numpoints, readinatoms_t *atoms, lorr_t lorr )
{
   kdtree_t *t ;
   int parentdim ;
   splitpoints_t partition ;
   int j ;

   t = malloc(sizeof(kdtree_t)) ;
   if (t == NULL) {
      Error("Out of memory on kdtree malloc()\n") ; }


   if (parent != NULL) {
      t->parent = parent ;
      t->lorr = lorr ;
   }


   if (numpoints == 1 ) {

      t->numpoints = 1 ;

      t->points = malloc(sizeof(int)) ;
      if (t->points == NULL) {
         Error("Error in t->points malloc()\n") ; }

      t->points[0] = points[0] ;

      t->left = NULL ;
      t->right = NULL ;

      free(points) ;

#ifdef DEBUG
      fprintf(stderr, "%d", t->points[0]) ;
#endif

   } else {

      if (parent != NULL) {
         t->bound_min[0] = t->parent->bound_min[0];
         t->bound_min[1] = t->parent->bound_min[1];
         t->bound_min[2] = t->parent->bound_min[2];
         t->bound_max[0] = t->parent->bound_max[0];
         t->bound_max[1] = t->parent->bound_max[1];
         t->bound_max[2] = t->parent->bound_max[2];

         parentdim = t->parent->depth % 3 ;
         if (t->lorr == CHILD_LEFT) {
            t->bound_max[parentdim] = t->parent->splitval ; }
         if (t->lorr == CHILD_RIGHT) {
            t->bound_min[parentdim] = t->parent->splitval ; }
      } else {
         t->bound_min[0] = atoms->coord_min[0] ;
         t->bound_min[1] = atoms->coord_min[1] ;
         t->bound_min[2] = atoms->coord_min[2] ;

         t->bound_max[0] = atoms->coord_max[0] ;
         t->bound_max[1] = atoms->coord_max[1] ;
         t->bound_max[2] = atoms->coord_max[2] ;
      }

      partition = split_points( points, numpoints, depth, atoms ) ;
      free(points) ;

      if (partition.Rsize == 0 ) {
         t->numpoints = partition.Lsize ;
         t->points = malloc(t->numpoints * sizeof(int)) ;
         if (t->points == NULL) {
            Error("Error in t->points malloc()\n") ; }

         for (j = 0; j < t->numpoints; j++) {
            t->points[j] = partition.Lset[j] ; }
         
         t->left = NULL ;
         t->right = NULL ;
      } else {
         t->splitval = partition.splitval ;
         t->depth = depth ;

#ifdef DEBUG
         fprintf(stderr, "(") ;
#endif

         t->left = build_kdtree( t, (depth + 1), partition.Lset, partition.Lsize, atoms, CHILD_LEFT ) ;

#ifdef DEBUG
         fprintf(stderr, ",") ;
#endif

         t->right = build_kdtree( t, (depth + 1), partition.Rset, partition.Rsize, atoms, CHILD_RIGHT ) ;

#ifdef DEBUG
         fprintf(stderr, ")") ;
#endif
      }

   }

   return t ;
}



/* split_points: splits an array of point numbers on the median value along
   the dimension corresponding to the current depth */
splitpoints_t split_points (int *points, int numpoints, int depth, readinatoms_t *atoms)
{
   int Dimension = (depth % 3) ;
   float *curdim ;

   int L, R ;

   int curLsize = INITSETSIZE ;
   int curRsize = INITSETSIZE ;

   int j ;

   splitpoints_t result ;

   result.Lset = malloc(curLsize * sizeof(int)) ;
   if (result.Lset == NULL) {
      Error("Out of memory on Lset malloc()\n") ; }

   result.Rset = malloc(curRsize * sizeof(int)) ;
   if (result.Rset == NULL) {
      Error("Out of memory on Rset malloc()\n") ; }


   curdim = malloc(numpoints * sizeof(float)) ;
   if (curdim == NULL) {
      Error("Out of memory on curdim malloc()\n") ; }


   for (j = 0; j < numpoints; j++ ) {
      curdim[j] = atoms->details[points[j]].coord[Dimension] ; }

   result.splitval = median_quick_select(curdim, numpoints) ;

   free(curdim) ;


   L = 0;
   R = 0;

   for (j = 0; j < numpoints; j++ ) {

      if ( atoms->details[points[j]].coord[Dimension] <= result.splitval) {

         if (L >= curLsize) {
            int *newp;
            curLsize += SETBLOCKSIZE ;
            newp = realloc(result.Lset, curLsize * sizeof(atom_t)) ;
            if (newp == NULL) {
               Error("Out of Memmory on realloc()\n") ; }
            result.Lset = newp ;
         }

         result.Lset[L] = points[j] ;
         L++ ;

      } else {

         if (R >= curRsize) {
            int *newp;
            curRsize += SETBLOCKSIZE ;
            newp = realloc(result.Rset, curRsize * sizeof(atom_t)) ;
            if (newp == NULL) {
               Error("Out of Memmory on realloc()\n") ; }
            result.Rset = newp ;
         }

         result.Rset[R] = points[j] ;
         R++ ;
      }

   }

   result.Lsize = L ;
   result.Rsize = R ;

   return result ;
}


#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

/*
   median_quick_select:
   returns the median of an array of floats using the Quick Select algorithm 

 ``Fast median search: an ANSI C implementation'', Nicolas Devillard
 http://ndevilla.free.fr/median/median/src/quickselect.c

*/
float median_quick_select(float arr[], int n)  
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }

}

#undef ELEM_SWAP


/* display_contacts; for each atom displays a list of other atoms within a given radius */
void display_contacts (kdtree_t *atoms_kdtree, readinatoms_t *atoms, float *radius)
{
   int j ;

   printf("#resna1\tresno1\tinscode1\tchain_id1\tatomno1\tatomna1\tresna2\tresno2\tinscode2\tchain_id2\tatomno2\tatomna2\tdistance\n") ;

   for (j = 0; j < atoms->number ; j++) {

#ifdef DEBUG
      fprintf(stderr, "\n%s %d %s %d (%d):\n", atoms->details[j].resna,
         atoms->details[j].resno, atoms->details[j].atomna,
         atoms->details[j].atomno, j ) ;
#endif

      atom_contacts(j, radius, atoms_kdtree, atoms ) ;
   }

#ifdef DEBUG
   fprintf(stderr, "\n") ;
#endif
}


/* atom_contacts; displays all atoms within a radius of a given atom*/
void atom_contacts (int centerind, float *radius, kdtree_t *atoms_kdtree, readinatoms_t *atoms)
{
   float *center ;
   float rect_max[3] ;
   float rect_min[3] ;

   int i ;

   center = atoms->details[centerind].coord ;

#ifdef DEBUG
   fprintf(stderr, "searching for (%f, %f, %f) +/- %f\n", center[0], center[1], center[2], *radius) ;
#endif

   for (i = 0; i < 3; i++ ) {
      rect_min[i] = center[i] - *radius ;
      rect_max[i] = center[i] + *radius ;
   }
   
   search_kdtree(atoms_kdtree, rect_min, rect_max, centerind, center, radius, atoms) ;
}


/* search_kdtree: searches the kdtree with a hypercube query range*/
void search_kdtree (kdtree_t *cur_kdtree, float *rect_min, float *rect_max, int centerind, float *center, float *radius, readinatoms_t *atoms)
{
   range_range_t rangecomp ;

   if ((cur_kdtree->left == NULL) && (cur_kdtree->right == NULL) ) {
      report_kdtree(cur_kdtree, centerind, center, radius, atoms) ;
      return ;
   }

   rangecomp = query_vs_node( cur_kdtree, rect_min, rect_max) ;

   if (rangecomp == CONTAINED) {
      report_kdtree(cur_kdtree, centerind, center, radius, atoms) ;

   } else if (rangecomp == INTERSECT) {

      int d = cur_kdtree->depth % 3 ;
      
      if ( cur_kdtree->splitval < rect_min[d]) {

#ifdef DEBUG
         for (i = 0; i < cur_kdtree->depth ; i++) {
            fprintf(stderr, " ") ; }
         fprintf(stderr, "go right\n") ;
#endif

         search_kdtree(cur_kdtree->right, rect_min, rect_max,
                       centerind, center, radius, atoms) ;

      } else if ( cur_kdtree->splitval > rect_max[d]) {

#ifdef DEBUG
         for (i = 0; i < cur_kdtree->depth ; i++) {
            fprintf(stderr, " ") ; }
         fprintf(stderr, "go left\n") ;
#endif

         search_kdtree(cur_kdtree->left, rect_min, rect_max,
                       centerind, center, radius, atoms) ;

      } else {

#ifdef DEBUG
         for (i = 0; i < cur_kdtree->depth ; i++) {
            fprintf(stderr, " ") ; }
         fprintf(stderr, "go right\n") ;
#endif

         search_kdtree(cur_kdtree->right, rect_min, rect_max,
                       centerind, center, radius, atoms) ;

#ifdef DEBUG
         for (i = 0; i < cur_kdtree->depth ; i++) {
            fprintf(stderr, " ") ; }
         fprintf(stderr, "go left\n") ;
#endif

         search_kdtree(cur_kdtree->left, rect_min, rect_max,
                          centerind, center, radius, atoms) ;

      }

   }

#ifdef DEBUG
   else if (rangecomp == EMPTY) {
      for (i = 0; i < cur_kdtree->depth ; i++) {
         fprintf(stderr " ") ; }
      fprintf(stderr, "X\n") ;
   }
#endif

   return ;
}


/* report_kdtree: traverses a kd-tree and displays all terminal nodes */
void report_kdtree( kdtree_t *cur_kdtree, int queryind, float *center, float *radius, readinatoms_t *atoms)
{
   int j ;

   if ((cur_kdtree->left == NULL) && (cur_kdtree->right == NULL) ) {
      for (j = 0; j < cur_kdtree->numpoints; j++) {
         int targetind= cur_kdtree->points[j] ;
         float dist = eucliddist(atoms->details[targetind].coord, center ) ;

         if ((dist <= *radius) && (queryind != targetind)) {
	    printf("%s\t%d\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%f\n",
                   atoms->details[queryind].resna,
		   atoms->details[queryind].resno,
		   atoms->details[queryind].inscode,
		   atoms->details[queryind].chainid,
		   atoms->details[queryind].atomno,
		   atoms->details[queryind].atomna,

                   atoms->details[targetind].resna,
		   atoms->details[targetind].resno,
		   atoms->details[targetind].inscode,
		   atoms->details[targetind].chainid,
		   atoms->details[targetind].atomno,
		   atoms->details[targetind].atomna,

                   dist ) ;
         }

      }
   } else {
      report_kdtree( cur_kdtree->left, queryind, center, radius, atoms ) ;
      report_kdtree( cur_kdtree->right, queryind, center, radius, atoms ) ;
   }

   return ;
}


/* query_vs_node: determines the overlap between the query and the ndoe hypercubes ranges*/
range_range_t query_vs_node( kdtree_t *cur_kdtree, float *rect_min, float *rect_max )
{
   range_range_t overlap[3] ;
   range_range_t result ;

   point_point_t qmin_bmin ;
   point_point_t qmin_bmax ;
   point_point_t qmax_bmin ;
   point_point_t qmax_bmax ;

   int i ;


#ifdef DEBUG
   fprintf(stderr, "*** Node:\t%f - %f\t%f - %f\t%f - %f\n", cur_kdtree->bound_min[0], cur_kdtree->bound_max[0], cur_kdtree->bound_min[1], cur_kdtree->bound_max[1], cur_kdtree->bound_min[2], cur_kdtree->bound_max[2]) ;

   fprintf(stderr, "  - Query:\t%f - %f\t%f - %f\t%f - %f\n", rect_min[0], rect_max[0], rect_min[1], rect_max[1], rect_min[2], rect_max[2]) ;
#endif


   for (i = 0; i < 3; i++) {

      if (rect_min[i] < cur_kdtree->bound_min[i]) {
         qmin_bmin = LESS;
      } else if (rect_min[i] > cur_kdtree->bound_min[i]) {
         qmin_bmin = MORE;
      } else if (rect_min[i] == cur_kdtree->bound_min[i]) {
         qmin_bmin = SAME;
      }


      if (rect_min[i] < cur_kdtree->bound_max[i]) {
         qmin_bmax = LESS;
      } else if (rect_min[i] > cur_kdtree->bound_max[i]) {
         qmin_bmax = MORE;
      } else if (rect_min[i] == cur_kdtree->bound_max[i]) {
         qmin_bmax = SAME;
      }

      if (rect_max[i] < cur_kdtree->bound_min[i]) {
         qmax_bmin = LESS;
      } else if (rect_max[i] > cur_kdtree->bound_min[i]) {
         qmax_bmin = MORE;
      } else if (rect_max[i] == cur_kdtree->bound_min[i]) {
         qmax_bmin = SAME;
      }

      if (rect_max[i] < cur_kdtree->bound_max[i]) {
         qmax_bmax = LESS;
      } else if (rect_max[i] > cur_kdtree->bound_max[i]) {
         qmax_bmax = MORE;
      } else if (rect_max[i] == cur_kdtree->bound_max[i]) {
         qmax_bmax = SAME;
      }

      overlap[i] = EMPTY ;

      if ( ((qmin_bmin == LESS) || (qmin_bmin == SAME)) &&
           ((qmax_bmax == MORE) || (qmax_bmax == SAME))  ) {
         overlap[i] = CONTAINED ;

      } else if ( (((qmin_bmin == SAME) || (qmin_bmin == MORE)) &&
                   (qmax_bmax == LESS)) ||
                  (((qmax_bmax == SAME) || (qmax_bmax == LESS)) &&
                   (qmin_bmin == MORE)) ) {
         overlap[i] = INTERSECT ;

      } else if ( ( ((qmin_bmin == LESS) && (qmax_bmax == LESS)) &&
                    ((qmax_bmin == SAME) || (qmax_bmin == MORE))) ||

                  ( ((qmin_bmin == MORE) && (qmax_bmax == MORE)) &&
                    ((qmin_bmax == SAME) || (qmin_bmax == LESS)) ) ) {

         overlap[i] = INTERSECT ;
      }


#ifdef DEBUG
      switch (overlap[i]) {
         case EMPTY: fprintf(stderr, "%d is empty\n", i) ; break;
         case INTERSECT: fprintf(stderr, "%d is intersect\n", i) ; break;
         case CONTAINED: fprintf(stderr, "%d is contained\n", i) ; break;
      }
#endif

   }

   result = overlap[0] ;

   if (result != EMPTY) {
      for (i = 1; i < 3; i++) {

         if (overlap[i] == EMPTY) {
            result = EMPTY ;
            break;
         } else if (overlap[i] == INTERSECT) {
            result = INTERSECT;
         }

      }
   }

#ifdef DEBUG
   if (result == EMPTY) {
      fprintf(stderr, "   EMTPY\n") ;
   } else if (result == INTERSECT) {
      fprintf(stderr, "   INTERSECT\n") ;
   } else if (result == CONTAINED) {
      fprintf(stderr, "   CONTAINED\n") ;
   }
#endif

   return result ;
}



/* readinatoms: reads in ATOM records from STDIN and returns a pointer to a
   readinatoms)_t struct with a pointer to an array of atom_t structs,
   the number of atoms, and the coordinate minimum and maximum along each
   dimension*/
readinatoms_t *readinatoms (void)
{
   char line[MAXLINELENGTH] ;
   char tempsubstr[MAXLINELENGTH] ;

   int i, j;
   readinatoms_t *result ;
   int atomlistsize = INITNUMATOMS ;


   result = malloc(sizeof(readinatoms_t)) ;
   if (result == NULL) {
      Error("Out of memory on result malloc()\n") ;
   }

   result->details = malloc(atomlistsize * sizeof(atom_t)) ;
   if (result->details == NULL) {
      Error("Out of memory on details malloc()\n") ;
   }


   i = 0;
   j = 1;
   while ((fgets(line, sizeof(line), stdin)) && j ) {

      *(line+(strlen(line)-1)) = '\0';

      if ((line[0] == 'A') &&
          (line[1] == 'T') &&
          (line[2] == 'O') &&
          (line[3] == 'M')) {

         if( i >= atomlistsize) {
            atom_t *newp;
            atomlistsize += ATOMBLOCKSIZE ;
            newp = realloc(result->details, atomlistsize * sizeof(atom_t)) ;

            if (newp == NULL) {
               Error("Out of Memmory on realloc()\n") ;
            }

            result->details = newp ;
         }

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


         i++ ;

      } else if ( (line[0] == 'E') &&
                (line[1] == 'N') &&
                (line[2] == 'D') &&
                (line[3] == 'M') &&
                (line[4] == 'D') &&
                (line[5] == 'L') ) {

         j = 0;

      }

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
