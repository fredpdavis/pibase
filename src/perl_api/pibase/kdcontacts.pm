=head1 NAME

pibase::kdcontacts - package that parses kdcontacts contacts.

=head1 DESCRIPTION

Parses kdcontacts output and displays it in a tab-delimited format ready for
import into pibase.interatomic_contacts

=head1 AUTHOR

Fred P. Davis, HHMI-JFRC (davisf@janelia.hhmi.org)

=head1 LICENCE AND COPYRIGHT

Copyright 2005,2010 Fred P. Davis (davisf@janelia.hhmi.org).
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

package pibase::kdcontacts ;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT = qw/kdcontacts_2_pibase/ ;


=head2 kdcontacts_2_pibase()

   Title:       kdcontacts_2_pibase()
   Function:    parases kdcontacts output and reformats for PIBASE tables
   Args:        $_[0] = kdcontacts file name
                $_[1] = output file name
                $_[2] = bdp_id (optional)

   FILE OUT:    writes to $_[1] output file
      1. 'bdp_id'
      2. chain_id1
      3. resno1
      4. resno1_int
      5. resna1
      6. atomna1
      7. chain_id2
      8. resno2
      9. resno2_int
      10. resna2
      11. atomna2
      12. distance
      13. contact_type

   FILE IN:     $_[0] = kdcontacts file naem
      1. atom1.resna
      2. atom1.resno
      3. atom1.inscode
      4. atom1.chainid
      5. atom1.atomno
      6. atom1.atomna
      7. atom2.resna
      8. atom2.resno
      9. atom2.inscode
      10. atom2.chainid
      11. atom2.atomno
      12. atom2.atomna
      13. distance

=cut

sub kdcontacts_2_pibase {

   my ($kdcontactsf, $outfile, $bdp_id) = @_ ;

   if ($kdcontactsf eq $outfile) {
      croak("ERROR kdcontacts_cleaner(): the kdcontacts file and the output filehave the same name") ; }

   if (! defined $bdp_id) {
      $bdp_id = 'bdp_id' ; }

# Initialize inlist flag.

   my $inlist = 0;

   my $fields = [
      'resna1', 'resno1', 'inscode1', 'chain_id1', 'atomno1', 'atomna1',
      'resna2', 'resno2', 'inscode2', 'chain_id2', 'atomno2', 'atomna2',
      'dist' ] ;

# Read in a line of input from the kdcontacts output file.

   open (INFILE, $kdcontactsf) ;
   open (OUTFILE, ">$outfile") ;
   while (my $line = <INFILE>) {
      if ($line !~ /^#/) {

         chomp $line ;

# Set contact type to vdw.

         my $contact_type = 'vdw' ;

         my %values ;

         my @t = split(/\t/, $line) ;
         foreach my $j (0 .. $#t) {
            $values{$fields->[$j]} = $t[$j] ; }

         my $resno1 = $values{resno1}.$values{inscode1} ; $resno1 =~ s/ $// ;
	 my $resno1_int = $values{resno1} ;

         my $resno2 = $values{resno2}.$values{inscode2} ; $resno2 =~ s/ $// ;
	 my $resno2_int = $values{resno2} ;

# Generate signature for each atom

	 my $sig1 = $resno1."^".$values{chain_id1}."^".
                     $values{atomno1}."^".$values{atomna1} ;
	 my $sig2 = $resno2."^".$values{chain_id2}."^".
                     $values{atomno2}."^".$values{atomna2} ;

# Redundancy removal: if atom 1's signature is ascii-wise less than
#  atom 2's signature, designate and display output fields.

	 if ($sig1 lt $sig2) {
            my @values = ( $bdp_id,
                           $values{chain_id1}, $resno1, $resno1_int,
			   $values{resna1}, $values{atomna1},
			   $values{chain_id2}, $resno2, $resno2_int,
			   $values{resna2}, $values{atomna2},
			   $values{dist} ) ; 
            print OUTFILE join("\t", @values)."\n" ;

         }
      }
   }

   close(INFILE) ;

   return 1 ;

}

1 ;
