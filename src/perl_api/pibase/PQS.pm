=head1 NAME

pibase::PQS - perl module of routines that operate on PQS (EBI's Probable
Quaternary Structure server) release files.

=head1 DESCRIPTION

The PQS.pm module contains routines to process and format PQS release files for PIBASE import.

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

package pibase::PQS ;
use strict;
use warnings;
use Carp ;

require Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT = qw/altloc_check altloc_filter/ ;

use pibase qw/locate_binaries/;
use File::Temp qw/tempfile/ ;


=head2 pqs_clean_asalist()

   Title:       pqs_clean_asalist()
   Function:    cleans PQS ASALIST file for pibase import
   Args:        none
   STDIN:       ASALIST
   STDOUT:      pibase.pqs_asalist table

=cut

sub pqs_clean_asalist {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }

   my @headers = qw/pqs_id pdb_id split_no contact_type homo_hetero quaternary num_chain num_resid num_hetatm delta_SASA delta_Solvat_E num_buried num_salt_bridges num_disulfide RMSD assembly_formula calpha_only/ ;
   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }

   my %pqs_id_seen ;
   while (my $line = readline($fh->{in})) {
      if (($line =~ /^\#/) || ($line =~ /^$/)) {next;}
      chomp $line ;
      my ($pqs_id, $contact_type, $homo_hetero, $quaternary, $num_chain, $num_resid, $num_hetatm, $delta_SASA, $delta_Solvat_E, $num_buried, $num_salt_bridges, $num_disulfide, $RMSD, $assembly_formula) = split(' ', $line) ;

#If a line containing this pqs_id has been read before, skip to the next line. Otherwise, set the pqs_id_seen hash flag.

      if ($pqs_id_seen{$pqs_id}) {
         next;
      } else {
         $pqs_id_seen{$pqs_id} = 1 ;
      }

   
      my $calpha_only = 0 ;
   
#If the number of disulfide field is '**', set the calpha_only flag to tru.

      if ($num_disulfide eq '**') {
         $num_disulfide = '\N'; }

#Initialize the PDB id to the PQS id, and split number to NULL string.

      my ($pdb_id, $split_no) = ($pqs_id, '\N') ;

#If the PQS id contains '_', Extract PDB id and Split number from the PQS id.

      if ($pqs_id =~ /\_/) {
         ($pdb_id, $split_no) = ($pqs_id =~ /([^\_]+)\_([0-9]*)/) ; }

#Remove leading and trialing brackets from the assembly formula.

      $assembly_formula =~ s/\[// ;
      $assembly_formula =~ s/\]// ;

#If the assembly formula begins with 'na Ca-only', set the calpha only flag to 1 and the assesmbly formula to the NULL string.

      if ($assembly_formula =~ /^na Ca-only/) {
         $calpha_only = 1;
         $assembly_formula = '\N' ;
      }

#Store PQS id, PDB id, split number, contact type, homo/hetero, quaternary state, number of chains, number of residues, number of hetero atoms, change in ASA, change in Solvation energy, number of buried reisudes, number of salt bridges, number of disulfides, RMSD, assmebly formula, Calpha_only flag as output fields.

      my @outfields = ($pqs_id, $pdb_id, $split_no, $contact_type, $homo_hetero,
                       $quaternary, $num_chain, $num_resid, $num_hetatm,
                       $delta_SASA, $delta_Solvat_E, $num_buried,
                       $num_salt_bridges, $num_disulfide, $RMSD,
                       $assembly_formula, $calpha_only) ;

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


=head2 pqs_clean_biolist()

   Title:       pqs_clean_biolist()
   Function:    cleans PQS BIOLIST file for pibase import
   Args:        none
   STDIN:       BIOLIST (preferably pqs_clean_biolist_contract() fixed)
   STDOUT:      pibase.pqs_biolist table

=cut

sub pqs_clean_biolist {
   
   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }

   my @headers= qw/pdb_id action reason num_split quaternary_state memo/;
   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }

   my %pdb_id_seen ;
   while (my $line = readline($fh->{in})) {
      chomp $line ;
      if ($line =~ /^\#/) {next;}
      if ($line =~ /^IC$/) {next;}

      $line =~ s/\s+/ /g ;
      $line =~ s/[^[:ascii:]]//g ;
      
      my (@raw_line) = split(' ', $line) ;
      my $pdb_id = shift @raw_line ;

      if ($pdb_id_seen{$pdb_id}) {
         next;
      } else {
         $pdb_id_seen{$pdb_id} = 1 ;
      }

      my $action = shift @raw_line ;
      my $num_split = 0;
      my $quaternary_state = '\N';
      my $memo = '\N' ;
   
      my $reason ;

      if ($action eq 'SKIP') {

         if (shift @raw_line eq 'WATER') {
   	       $reason = 'WATER' ;

         } elsif (shift @raw_line eq 'coordinates') {
            $reason = 'no_coordinates' ;
   	 } else {
            $reason = 'no_REMARK_350' ;
         }

      } elsif ($action eq 'SKIP-NOT-XRAY') {

         $action = 'SKIP';
         $reason = join('_', @raw_line) ;

      } elsif (($action eq 'ASU-COMPLEX') ||
               ($action eq 'SYMMETRY-COMPLEX')) {

         my $raw_line = join(' ', @raw_line) ;
         $raw_line =~ s/Oligomeric file of type\s*// ;
         $quaternary_state = $raw_line ;

      } elsif (($action eq 'SPLIT-ASU') ||
	        ($action eq 'SPLIT-SYMMETRY')) {

         shift @raw_line ; # get rid of word into
         $num_split = shift @raw_line; 

         my $raw_line = join(' ', @raw_line) ;
   	 $raw_line =~ s/Oligomeric files of type\s*// ;
   	 $quaternary_state = $raw_line ;
      }

   # still possible: CONNECTIVE-FIBRE, HELICAL-VIRUS, ICOSAHEDRAL-VIRUS, MICROTUBULE, POLAR-VIRUS, VIRUS-COMPLEX
   
#If quaternary state contains 'memo', store everything before 'memo' as the quaternary state, and everything after 'memo' as thge memo field.

      if ($quaternary_state =~ /memo/) {
         ($quaternary_state, $memo) =
            ($quaternary_state =~ /(.+) memo (.+)$/);
      }

      my @outfields = ($pdb_id, $action, $reason, $num_split,
                       $quaternary_state, $memo) ;

      foreach my $j (0 .. $#outfields) {
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


=head2 pqs_clean_biolist_contraction()

   Title:       pqs_clean_biolist_contraction()
   Function:    fixes contraction errors in PQS BIOLIST
   Args:        none
   STDIN:       BIOLIST
   STDOUT:      contraction fixed BIOLIST

=cut

sub pqs_clean_biolist_contraction {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
   }

   my ($temp_fh, $temp_fn) = tempfile() ;
   while (my $line = readline($fh->{in})) {
      chomp $line ;
      if ($line =~ /type.*type/) {
         while ($line =~ /type.*type/) {
            $line =~ s/type (.*)type/type\n$1type/g ; }
      }
      print {$temp_fh} $line."\n" ;
   }
   close($fh->{in}) ;
   close($temp_fh) ;

   my $tcom = "cat $temp_fn | tr -d '".'\000-\011\013-\037\177-\377'."' " ;
   if (exists $in->{out_fn}) {
      $tcom .= " >".$in->{out_fn} }
   system($tcom) ;

}


=head2 pqs_clean_list()

   Title:       pqs_clean_list()
   Function:    reformats PQS LIST for pibase import
   Args:        none
   STDIN:       PQS LIST
   STDOUT:      pibase.pqs_list table

=cut

sub pqs_clean_list {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }

   my @headers = qw/pqs_file_base pdb_id/ ;
   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }

   while (my $line = readline($fh->{in})) {

#If the line is not a comment AND does not contain 'ceramide' or 'lipid', continue.

      if ( $line =~ /^\#/ ||
           $line =~ /^core$/ ||
	   $line =~ /^j$/  ||
           $line =~ /ceramide/ ||
	   $line =~ /lipid/ ) {next;}

      chomp $line ;
      my ($pqs_file_base) = ($line =~ /(\w+\.\w+)/) ;
      my $pdb_id = substr($pqs_file_base, 0, 4) ;
   
      my @outfields = ($pqs_file_base, $pdb_id) ;
      print {$fh->{out}} join("\t", @outfields)."\n" ;

   }
   close($fh->{in}) ;
   close($fh->{out}) ;

}


=head2 pqs_clean_ranking()

   Title:       pqs_clean_ranking()
   Function:    reformats PQS RANKING for pibase import
   Args:        none
   STDIN:       PQS RANKING
   STDOUT:      pibase.pqs_ranking table

=cut

sub pqs_clean_ranking {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }

   my @headers = qw/pqs_id pdb_id split_no ramach_most_fav PROCHECK_morris_g B_main B_side no_atoms no_hoh no_hetatm/ ;
   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }

   my %pqs_id_seen ;
   while (my $line = readline($fh->{in})) {

      if ($line =~ /^\#/) {next;}
      chomp $line ;
      my ($raw_pqs_id, $rama, $morris_g, $b_main, $b_side, $no_atom, $no_hoh,
             $no_hetatm) = split(' ', $line) ;

      my ($pqs_id) = ($raw_pqs_id =~ /(.+)\.mmol/ ) ;

      if ($pqs_id_seen{$pqs_id}) {
         next;
      } else {
         $pqs_id_seen{$pqs_id} = 1 ; 
      }

      my ($pdb_id, $split_no) = ($pqs_id, '\N') ;

      if ($pqs_id =~ /\_/) {
         ($pdb_id, $split_no) = ($pqs_id =~ /(.+)\_([0-9]*)/) ; } 

#Store PQS id, PDB id, split number, ramachandran most favored, Morris G factor, main chain mean B value, side chain mean B value, number of atoms, number of waters, and number of HETERO atoms as output fields.

      my @outfields = ($pqs_id, $pdb_id, $split_no, $rama, $morris_g,
                       $b_main, $b_side, $no_atom, $no_hoh, $no_hetatm) ;
 
      foreach my $j (0 .. $#outfields) {

         $outfields[$j] =~ s/^\s*// ;
         $outfields[$j] =~ s/\s*$// ;

         if ((!defined $outfields[$j]) ||
	        ($outfields[$j] eq '') ||
		($outfields[$j] eq 'na') ||
		($outfields[$j] eq 'nan') ) {
            $outfields[$j] = '\N' ;
         }       
      }
      print {$fh->{out}} join("\t", @outfields)."\n" ;
   }
   close($fh->{in}) ;
   close($fh->{out}) ;

}


=head2 get_pqs_filepath()

   Title:       get_pqs_filepath()
   Function:    returns file path to a PQS entry
   Args:	->{pqs_id} = PQS identifier
                [->{pibase_specs} = $pibase_specs] - optional
   Returns:     returns PQS entry filepath

=cut

sub get_pqs_filepath {
   my $in = shift ;

   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::complete_pibase_specs()
   } else {
      $pibase_specs = $in->{pibase_specs}; }

   my $filepath = $pibase_specs->{pqs_dir} ;
   if ($pibase_specs->{external_data}->{pqs}->{file_layout} eq 'wwpdb') {
      $filepath .= '/'.substr($in->{pqs_id},1,2) ;
   }
   $filepath .= '/'.lc($in->{pqs_id}) ;


   if (exists $pibase_specs->{external_data}->{pqs}->{compress_fl} &&
       $pibase_specs->{external_data}->{pqs}->{compress_fl} == 1) {
      $filepath .= '.gz' ; }

   return $filepath ;
}

1 ;
