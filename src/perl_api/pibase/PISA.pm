=head1 NAME

pibase::PISA - perl module of routines that operate on PISA (PDBe's Protein
Interfaces, Surfaces, and Assemblies) files.

=head1 DESCRIPTION

The PISA.pm module contains routines to process and format PISA release files for PIBASE import.

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

package pibase::PISA ;
use strict;
use warnings;
use Carp ;

require Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT = qw/altloc_check altloc_filter/ ;

use pibase qw/locate_binaries/;
use File::Temp qw/tempfile/ ;

=head2 pisa_clean_index()

   Title:       pisa_clean_index()
   Function:    cleans PISA index.txt file for pibase import
   Args:        none
   STDIN:       INDEX
   STDOUT:      pibase.pisa_index table

=cut

sub pisa_clean_index {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }

   my @headers = qw/pisa_id pdb_id split_no entry_type homo_hetero num_chain num_resid total_sasa buried_sasa num_hbonds num_salt_bridges num_disulfide diss_energy assembly_formula/ ;

   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }

   my %pisa_id_seen ;
   while (my $line = readline($fh->{in})) {
      if (($line =~ /^\#/) || ($line =~ /^$/)) {next;}
      chomp $line ;
      $line =~ s/TOO BIG/TOO_BIG/ ;
      my ($pisa_id, $entry_type, $homo_hetero, $num_chain, $num_resid, $total_sasa, $buried_sasa, $num_hbonds, $num_salt_bridges, $num_disulfide, $diss_energy, $assembly_formula) = split(' ', $line) ;

#If a line containing this pisa_id has been read before, skip to the next line. Otherwise, set the pisa_id_seen hash flag.

      if ($pisa_id_seen{$pisa_id}) {
         next;
      } else {
         $pisa_id_seen{$pisa_id} = 1 ;
      }

   
      my $calpha_only = 0 ;
   
#Initialize the PDB id to the PQS id, and split number to NULL string.

      my ($pdb_id, $split_no) = ($pisa_id, '\N') ;

#If the PQS id contains '_', Extract PDB id and Split number from the PQS id.

      if ($pisa_id =~ /\_/) {
         ($pdb_id, $split_no) = ($pisa_id =~ /([^\_]+)\_([0-9]*)/) ; }

#Remove leading and trialing brackets from the assembly formula.

      $assembly_formula =~ s/\[// ;
      $assembly_formula =~ s/\]// ;

#Store PISA id, PDB id, split number, entry type, homo/hetero, number of chains, number of residues, total SASA, buried SASA, number of h-bonds, number of salt bridges, number of disulfides, diss_energy, assembly formula as output fields.

      my @outfields = ($pisa_id, $pdb_id, $split_no, $entry_type, $homo_hetero,
                       $num_chain, $num_resid, $total_sasa, $buried_sasa,
                       $num_hbonds, $num_salt_bridges, $num_disulfide,
                       $diss_energy, $assembly_formula) ;

#Iterate through the output fields

      foreach my $j (0 .. $#outfields) {

         $outfields[$j] =~ s/^\s*// ;
         $outfields[$j] =~ s/\s*$// ;

         if ((!defined $outfields[$j]) ||
	     ($outfields[$j] eq '') ||
   	     ($outfields[$j] eq 'na')) {
               $outfields[$j] = '\N' ; }
      }
      print {$fh->{out}} join("\t", @outfields)."\n" ;
   }
   close($fh->{in}) ;
   close($fh->{out}) ;

}


=head2 get_pisa_filepath()

   Title:       get_pisa_filepath()
   Function:    returns file path to a PISA entry
   Args:	->{pisa_id} = PISA identifier
                [->{pibase_specs} = $pibase_specs] - optional
   Returns:     returns PISA entry filepath

=cut

sub get_pisa_filepath {

   my $in = shift ;

   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::complete_pibase_specs()
   } else {
      $pibase_specs = $in->{pibase_specs}; }

   my $filepath = $pibase_specs->{pisa_dir} ;
   if ($pibase_specs->{external_data}->{pisa}->{file_layout} eq 'wwpdb') {
      $filepath .= '/'.substr($in->{pisa_id},1,2) ;
   }
   $filepath .= '/'.lc($in->{pisa_id}).'_pisa.pdb' ;


   if (exists $pibase_specs->{external_data}->{pisa}->{compress_fl} &&
       $pibase_specs->{external_data}->{pisa}->{compress_fl} == 1) {
      $filepath .= '.gz' ; }

   return $filepath ;
}

