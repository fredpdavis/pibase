/* subset_extractor.c - Extracts specified residues from a PDB file

Purpose: extracts specified residues from a PDB file
Usage: ./subset_extractor pdbfile < subset definitions
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
#include<string.h>


//#define DEBUG 1
#define MAXLINELENGTH 81
#define INITNUMSEGS 100
#define SEGBLOCKSIZE 100
#define Error( Str )   fprintf( stderr, "%s\n", Str ), exit( 1 )


//STRUCTURES
struct segment_Struct {
   int		read;
   char		chainid[2] ;
   char		start[6] ;
   char		end[6] ;
} ;
typedef struct segment_Struct segment_t ;

struct readinsegments_Struct {
   int		number ;
   int		num_unread ;
   segment_t	*details ;
} ;
typedef struct readinsegments_Struct readinsegments_t ;



//FUNCTIONS
readinsegments_t *readinsegments() ;
void extractsegments( FILE *pdb_fp, readinsegments_t *segments ) ;
void deblank( char *l ) ;
char *st_sep (char **stringp, const char *delim) ;



int main(int argc, char *argv[])
{
   readinsegments_t *segments ;
   char *pdb_fn ;
   FILE *pdb_fp ;


   if (argc > 1 ) {
      pdb_fn = argv[1] ;
   } else {
      Error("usage: subset_extractor pdbfile < subset_definition") ;
   }

   segments = readinsegments() ;

#ifdef DEBUG
   fprintf(stderr, "read %d segments\n", segments->number) ;
#endif

   pdb_fp = fopen(pdb_fn, "r") ;
   if (pdb_fp == NULL) {
      fprintf(stderr, "ERROR: PDB file %s does not exist\n", pdb_fn) ;
      exit(1) ;
   }
   extractsegments(pdb_fp, segments) ;
   fclose(pdb_fp) ;

   return 0;
}


readinsegments_t *readinsegments( FILE *pdb_fp )
{
   int i = 0 ;
   char line[MAXLINELENGTH] ;
   char *line_p ;
   char *tempsubstr ;

   readinsegments_t *results ;
   int seglistsize = INITNUMSEGS ;

   results = malloc(sizeof(readinsegments_t)) ;
   if (results == NULL ) {
      Error("Out of memory on result malloc()\n") ;
   }

   results->details = malloc(seglistsize * sizeof(segment_t)) ;
   if (results->details == NULL) {
      Error("Out of memory on details malloc()\n") ;
   }

   i = 0 ;
   while(fgets(line, sizeof(line), stdin)) {

      *(line+(strlen(line)-1)) = '\0' ;
      line_p = line ;

      if( i >= seglistsize) {
         segment_t *newp;
	 seglistsize += SEGBLOCKSIZE ;
	 newp = realloc(results->details, seglistsize * sizeof(segment_t)) ;

         if (newp == NULL) {
	    Error("Out of Memmory on realloc()\n") ;
         }

         results->details = newp ;
      }

      tempsubstr = st_sep(&line_p, "\t") ;
      strcpy(results->details[i].chainid, tempsubstr) ;

      tempsubstr = st_sep(&line_p, "\t") ;
      strcpy(results->details[i].start, tempsubstr) ;
      tempsubstr = st_sep(&line_p, "\n") ;
      strcpy(results->details[i].end, tempsubstr) ;

#ifdef DEBUG
      fprintf(stderr, "segment %d: chain %s, start %s, end %s\n", i,
		      results->details[i].chainid,
		      results->details[i].start,
		      results->details[i].end) ;
#endif

      results->details[i].read = 0 ;
      i++ ;
   }

   results->number = i ;
   results->num_unread = i ;

#ifdef DEBUG
   fprintf(stderr, "num_unread is %d\n", results->num_unread) ;
#endif

   return results ;
}


void extractsegments( FILE *pdb_fp, readinsegments_t *segments )
{
   char line[MAXLINELENGTH] ;

   int j  = 1;

   char chainid[2] ;
   char resno_full[6] ;

   int inseg_fl = 0 ;
   int curseg = -1 ;
   int lastres_fl = 0 ;

   char resno_last[6] ;
   char chainid_last[2] ;


   while ((fgets(line, sizeof(line), pdb_fp)) && j ) {

      *(line+(strlen(line))) = '\0';

// fpd080611_2014 added this line to fix extra blank lines when extracting
//   subsets from PQS files which have ATOM records less than full 81 char
      if (line[strlen(line) - 1] == '\n') {
         *(line+(strlen(line) - 1)) = '\0'; }

      if ((line[0] == 'A') &&
          (line[1] == 'T') &&
          (line[2] == 'O') &&
          (line[3] == 'M')) {


         strncpy(chainid, (line + 21), 1) ;
         chainid[1] = '\0' ;
      
         strncpy(resno_full, (line + 22), 5) ;
         *(resno_full+(5)) = '\0' ;

	 deblank(resno_full) ;
	 deblank(chainid) ;

	 if ( lastres_fl &&
	      ((strcmp(resno_full, resno_last) != 0) ||
               (strcmp(chainid, chainid_last) != 0))) {

	    lastres_fl = 0 ;
	    inseg_fl = 0 ;
	    segments->details[curseg].read = 1 ;
	    segments->num_unread-- ;

	    if (segments->num_unread == 0 ) {
	       exit(1) ;
	    }

	 }

	 if (inseg_fl) {
	    int chain_nullfl = 0 ;
	    int end_nullfl = 0 ;

	    int chain_match = 0 ;
	    int end_match = 0 ;

	    if ( strcmp(segments->details[curseg].end, "") == 0  ) {
	       end_nullfl = 1 ; 
	    }

	    if ( strcmp(segments->details[curseg].chainid, "") == 0 ) {
	       chain_nullfl = 1 ; 
	    }

	    if (strcmp(resno_full, segments->details[curseg].end) == 0) {
	       end_match = 1 ;
	    }

	    if (strcmp(chainid, segments->details[curseg].chainid) == 0) {
	       chain_match = 1 ;
	    }

	    if ((end_nullfl) && (! chain_nullfl) && (! chain_match)) {
//               fprintf(stderr, " set out of inseg_fl\n") ;
	       inseg_fl = 0 ;
	    } else if ((chain_nullfl) && (! end_nullfl) && (end_match)) {
               lastres_fl = 1 ;
	    } else if ((!chain_nullfl) && (!end_nullfl) &&
	               (chain_match) && (end_match) ) {
               lastres_fl = 1 ;
	    }

//if lastres_fl, need to delete the current segment out of the segment list
// or just set a flag that the segment has been read alreaxdy

	    strcpy(resno_last, resno_full) ;
	    strcpy(chainid_last, chainid) ;

#ifdef DEBUG
	       fprintf(stderr, "in segment\n") ;
#endif

	 }

	 if (! inseg_fl) {

	    int k ; 
	    for (k = 0; k < segments->number; k++) {

	       if (!segments->details[k].read) {

	          int chain_nullfl = 0 ;
		  int start_nullfl = 0 ;

	          int chain_match = 0 ;
		  int start_match = 0 ;

	          if ( strcmp(segments->details[k].start, "") == 0 ) {
	             start_nullfl = 1 ;
                  }

	          if ( strcmp(segments->details[k].chainid, "") == 0 ) {
	             chain_nullfl = 1 ; 
	          }


	          if ((chain_nullfl) && (start_nullfl)) {
	             inseg_fl = 1 ;
		     curseg = k ;
		     break ;
	          }

	          if ((strcmp(resno_full, segments->details[k].start) == 0) ||
	              (start_nullfl)) {
	             start_match = 1 ;
	          }

	          if ((strcmp(chainid, segments->details[k].chainid) == 0) ||
	              (chain_nullfl)) {
	             chain_match = 1 ;
	          }

#ifdef DEBUG
                  fprintf(stdout, "comparing %s on chain %s to segment res %s on chain %s (res match %d, chain match %d)\n", resno_full, chainid, segments->details[k].start, segments->details[k].chainid, start_match, chain_match) ;
#endif


	          if ( (start_match) && (chain_match) ) {
	             inseg_fl = 1 ;
		     curseg = k ;
		     break ;
	          }

	       }

	    }

	 }

	 if (inseg_fl) {
	    printf("%s\n",line) ;
	 }

	 j = 1 ;

      } else if ( (line[0] == 'E') &&
                  (line[1] == 'N') &&
                  (line[2] == 'D') &&
                  (line[3] == 'M') &&
                  (line[4] == 'D') &&
                  (line[5] == 'L') ) {

	 inseg_fl = 0 ;
         j = 0;

      }

   }

}


// FROM: http://www.harpercollege.edu/bus-ss/cis/166/mmckenzi/samples/lect10c.c
/* strip leading, trailing and double spaces */
void deblank( char *l )
{
   int i;
   char *p;

   i = strspn( l, " " ); /* find 1st non-blank */
   if ( i > 0 ) {
      for( p = l; ; p++ ) {
	 *p = *(p+i);
         if (*p == '\0') break;
      }
   }

   while( (p = strstr( l, "  ") ) != NULL ) {
      i = strspn( p, " "); /* find length of blank prefix */
      do {
         ++p;
         *p = *(p+i-1);
      } while( *p );
   }

   if ( ( i = strlen( l ) - 1 ) >= 0 ) {
      l[i] = ( l[i] == ' ' ) ? '\0' : l[i]; /* make trailing space a null */
   }

}


// edited version of GNU libc strsep() grokked from google gruops
char *st_sep (char **stringp, const char *delim)
{
   char *begin, *end;

   begin = *stringp;
   if (begin == NULL)
      return NULL;

  /* A frequent case is when the delimiter string contains only one
     character.  Here we don't need to call the expensive `strpbrk'
     function and instead work using `strchr'.  */
   if (delim[0] == '\0' || delim[1] == '\0') {
      char ch = delim[0];

      if (ch == '\0')
         end = NULL;
      else {
         if (*begin == ch)
	    end = begin;
         else
	    end = strchr (begin + 1, ch);
      }
   } else
    /* Find the end of the token.  */
      end = strpbrk (begin, delim);

   if (end) {
      /* Terminate the token and set *STRINGP past NUL character.  */
      *end++ = '\0';
      *stringp = end;
   } else
    /* No more delimiters; this is the last token.  */
       *stringp = NULL;

   return begin;
}
