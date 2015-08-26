=head1 NAME

pibase::residue_math - package that handles residue number operations

=head1 DESCRIPTION

The pibase::resno module performs common operations on residue numbers.

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


=head1 SUBROUTINES

=cut

package pibase::residue_math ;
use strict;
use warnings;
require Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/residue_int residue_add residue_comp residue_inrange/ ;


=head2 residue_int(resno)

   Title:       residue_int()
   Function:    Seperates the integer and insertion code of a residue number
      NOTE: assumes insertion code is always alphanumeric
   Args:        residue number (5 character - number and insertion code)
   Returns:     $_[0] integer portion of the residue number
                $_[1] insertion code of the residue number

=cut

sub residue_int {

   my $res_no = shift ;

   if (!defined $res_no) {
      die "res_no was undefined\n"; }

# Extract the integer and insertion code from the residue number.

   my ($res_no_int, $ins_code) = ($res_no =~ /(-?[0-9]+)(.?)/) ;

# If the insertion code is undefined, set to ''.

   if (!defined $ins_code) {
      $ins_code = ''; }

   return ($res_no_int, $ins_code) ;

}


=head2 residue_add(residue number, increment)

   Title:       residue_add()
   Function:    Adds an integer increment to a residue number
   Args:        $_[0] = residue number (full)
                $_[1] = increment
   Returns:     $_ = new residue number

=cut

sub residue_add {

   my ($resno, $incr) = @_ ;

# Extract integer and insertion portion of the residue number

   my ($int, $ins) = residue_int($resno) ;

# Add the incrementer to the integer portion, and concatenate the insertion code

   my $resno_new = ($int + $incr).$ins ;

   return $resno_new ;

}


=head2 residue_comparison(resno1, resno2)

   Title:       residue_comparison(resno1, resno2)
   Function:    Compares 2 residues to determine which one is greater.
   Args:        $_[0] = residue number 1
                $_[2] = residue number 2
   Returns:     $_ = comparison result: 0 = equal,
                                        1 = first is greater,
                                        2 = second is greater.

=cut

sub residue_comp {

   my ($res1, $res2) = @_ ;
   my ($res1_int, $res1_ins) = residue_int($res1) ;
   my ($res2_int, $res2_ins) = residue_int($res2) ;


# First compare the integer portions.

   my $result ;

   if ($res1_int < $res2_int) {
      $result = 2; }
   elsif ($res1_int > $res2_int) {
      $result = 1; }
   elsif ($res1_int == $res2_int) {
      $result = 0; }

# If the integer portions are equal, compare the insertion codes via an
#  ascii comparator.

   if ($result == 0 ) {
      if ($res1_ins lt $res2_ins) {
         $result = 2; }
      elsif ($res1_ins gt $res2_ins) {
         $result = 1; }
      elsif ($res1_ins eq $res2_ins) {
         $result = 0; }
   }

   return $result ;

}


=head2 residue_inrange(resno, start, end)

   Title:       residue_inrange(resno, start, end)
   Function:    Checks if a residue number is within a range of residue numbers.
   Args:        $_[0] = residue number
                $_[1] = start residue number range
                $_[1] = end residue number range
   Returns:     Comparison result: 0 = out of range, 1 = in range.

=cut

sub residue_inrange {

   my ($resno, $start, $end) = @_ ;

   my $a = residue_comp($resno, $start) ;
   my $b = residue_comp($resno, $end) ;

   my $in_range = 0;
   if ( ($a == 0 || $a == 1) &&
        ($b == 0 || $b == 2) ) {
      $in_range = 1;
   }

   return $in_range ;

}

1 ;
