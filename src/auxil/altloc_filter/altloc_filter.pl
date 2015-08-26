#!/usr/local/bin/perl
=head1 NAME

altloc.filter.pl - filters out duplicate ATOM and HETATOM records
(alternative loations)

=head1 DESCRIPTION

Filters a PDB file so only one instance of each ATOM and HETATOM exists. If occupancy information is available, keeps the highest occupancy instance; otherwise, keeps the first instance.

=head1 SYNOPSIS

B<altloc.filter.pl> pdb_file


=head1 AUTHOR

Fred P. Davis, HHMI-JFRC (davisf@janelia.hhmi.org)

=head1 LICENCE AND COPYRIGHT

Copyright 2008 Fred P. Davis (davisf@janelia.hhmi.org).
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

=cut


use warnings;
use strict ;

main() ;

=head2 SUB main()

=cut

sub main {

   my $usage = "usage: altloc.filter.pl pdbfile" ;
   if ($#ARGV < 0 ) { die $usage."\n" ;}
   my $pdb_fn = $ARGV[0] ;

   my %records_seen ;

   my $hidelines ;

   my $linenum = 0 ;
   open (PDBF, $pdb_fn) ;
   while (my $line = <PDBF>) {
      chomp $line;

      my %values ;

      if ( ($line =~ /^ATOM/) || ($line =~ /^HETATM/) ) {

	 my $sig = substr($line, 0, 6).
	           substr($line, 12, 4).
		   substr($line, 17, 3).
		   substr($line, 21, 1).
		   substr($line, 22, 4).
		   substr($line, 26, 1) ;

	 if (exists $records_seen{$sig}) {
	    my $curr_occup = substr($line, 55, 6) ;
	    if ($curr_occup > $records_seen{$sig}->[1]) {
	       $hidelines->[$records_seen{$sig}->[0]]++ ;
	       $records_seen{$sig}->[0] = $linenum ;
	       $records_seen{$sig}->[1] = $curr_occup ;
	    } else {
	       $hidelines->[$linenum]++ ;
	    }
	 } else {
            $records_seen{$sig} = [ $linenum, substr($line, 55, 6 ) ] ;
	 }

      }

      $linenum++ ;
   }
   close (PDBF) ;


   open(PDBF, $pdb_fn) ;
   my $j = 0 ;
   while (my $line = <PDBF>) {
      if (!exists $hidelines->[$j]) {
#set altloc code to blank
         if ( ($line =~ /^ATOM/) || ($line =~ /^HETATM/) ) {
            substr($line, 16, 1) = ' ' ; }
         print $line ;
      }
      $j++ ;
   }
   close(PDBF) ;


}
