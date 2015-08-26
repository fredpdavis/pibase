/* inscode_check.c - Parses a PDB file to determine if insertion code is used

Purpose: output 0 or 1 if a pdb file does or doesn't contain insertion codes
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


#define MAXLINELENGTH 80

int checkinscode( void ) ;


int main(int argc, char *argv[])
{

   int answer ;
   answer = checkinscode() ;
   printf("%d\n", answer) ;

   return 0;
}



/* readinatoms: reads in ATOM records from STDIN and returns a pointer to a
   readinatoms)_t struct with a pointer to an array of atom_t structs,
   the number of atoms, and the coordinate minimum and maximum along each
   dimension*/
int checkinscode(void)
{
   char line[MAXLINELENGTH] ;
   char tempsubstr[2] ;

   int i, j, uses_inscode ;


   i = 0;
   j = 1;
   uses_inscode = 0;
   while ((fgets(line, sizeof(line), stdin)) &&
          (! uses_inscode) &&
	  j) {

      *(line+(strlen(line)-1)) = '\0';

      if ((line[0] == 'A') &&
          (line[1] == 'T') &&
          (line[2] == 'O') &&
          (line[3] == 'M')) {

         strncpy(tempsubstr, (line + 26), 1) ;
         tempsubstr[1] = '\0' ;

	 if (tempsubstr[0] != ' ') {
//	    fprintf(stderr,"uses inscode on atom %d\n", i) ;
	    uses_inscode = 1 ; }
      
         i++ ;

      } else if ( (line[0] == 'E') &&
                (line[1] == 'N') &&
                (line[2] == 'D') &&
                (line[3] == 'M') &&
                (line[4] == 'D') &&
                (line[5] == 'L') ) {

	 j = 0 ;

      }

   }

   return uses_inscode ;
}
