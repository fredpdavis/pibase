=head1 NAME

pibase::subsets - perl module that deals with PDB subsets

=head1 DESCRIPTION

Perl module to handle subset operations on PDB files

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

package pibase::PDB::subsets ;
use strict;
use warnings;
use Exporter;
use Carp ;
our @ISA = qw/Exporter/ ;
our @EXPORT = qw/subset_extract/ ;

use pibase qw/locate_binaries/ ;
use File::Temp qw/tempfile/ ;


=head2 SUB subset_extract()

   Function: extracts specified residues/chains from PDB file
   Args:        $_[0] = PDB file name
                $_[1] = output PDB file name
                $_[2] = chain identifier
                $_[3] = start residue number
                $_[4] = end residue number

   Returns:     $_->[i] arrayref of errors

   Files IN:    PDB file ($_[0])
   Files OUT:   subset PDB file ($_[1])

=cut

sub subset_extract {

   my $in = shift ;
   my $pdb_fn = $in->{in_fn} ;
   my $cutpdb_fn= $in->{out_fn} ;
   my $chain_id = $in->{chain} ;
   my $resno_start = $in->{start} ;
   my $resno_end = $in->{end} ;

   my $binaries = pibase::locate_binaries() ;

   my ($temp_fh, $localbdp) = tempfile(SUFFIX => ".pdb") ; close($temp_fh) ;
   if ($pdb_fn =~ /\.gz$/) {
      my $tcom = $binaries->{zcat}." $pdb_fn > ".$localbdp ;
      system($tcom) ;
   } else {
      pibase::safe_copy($pdb_fn, $localbdp) ;
   }

   if (!-e $localbdp) {
      return {error_fl => "ERROR: pdb file access error - couldnt copy locally"}; }


   my $errors ;

   if ($binaries->{'subset_extractor'} eq 'ERROR') {
      croak("ERROR subset_extract(): subset_extractor binary not found") ; }

   if ($pdb_fn eq $cutpdb_fn) {
      croak("ERROR subset_extract(): the original and cut pdb filenames are the same") ; }

   my $compress_out_fl = 0 ;
   if ($cutpdb_fn =~ /\.gz$/) {
      $compress_out_fl = 1 ;
      $cutpdb_fn =~ s/\.gz$// ;
   }

   my ($subsetdef_fh, $subsetdef_fn) = tempfile("subsetdef.XXXXXX") ;
   my $errfile_fn = $subsetdef_fn.".err" ;

   foreach my $j ( 0 .. $#{$chain_id} ) {
      my @outvals = ($chain_id->[$j], $resno_start->[$j], $resno_end->[$j]) ;
      print $subsetdef_fh join("\t", @outvals)."\n" ;
   }
   close ($subsetdef_fh) ;

   my $tcom = "$binaries->{subset_extractor} $localbdp < $subsetdef_fn 2> $errfile_fn >$cutpdb_fn" ;
   system($tcom) ;
   if ($compress_out_fl) {
      system("gzip ".$cutpdb_fn) ;
      $cutpdb_fn .= '.gz' ;
   }
   unlink $subsetdef_fn ;

   if (-s $errfile_fn) {
      push @{$errors}, "subset_extractor error, see: $errfile_fn" ;
   } else {
      unlink $errfile_fn ;
   }

   {
      my ($newdev, $newino) = stat($localbdp) ;
      my ($olddev, $oldino) = stat($pdb_fn) ;
      if ($newdev != $olddev || $newino != $oldino) {
         unlink $localbdp ; }
   }

   return $errors ;

}

1 ;
