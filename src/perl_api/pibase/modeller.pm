=head1 NAME

pibase::modeller - module containing routines to call MODELLER for pibase

=head1 DESCRIPTION

Performs MODELLER operations needed by pibase. (still old-school TOP format)

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

package pibase::modeller ;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT = qw/subsets_2_mod_pick subsetdef_2_mod_pick modeller_subset_sasa get_vol cutpdb parse_ali get_resequiv get_dihedrals get_salign calc_sasa/ ;

use pibase ;
use File::Temp qw/tempfile/ ;
use pibase::PDB::residues ;


=head2 subsets_2_modpick()

   Function:    Converts subsets_details to MODELLER SELECTION SEGMENTS
   Args:        $_[0] = subset_id
                $_[1] = DBI db handle to pibase
   Return:      $_->[] arrayref of MODELLER pick statements to select domain

=cut

sub subsets_2_mod_pick {

   my ($subset_id, $dbh) = @_ ;

   my ( $chain_id, $start_resno, $end_resno ) = pibase::mysql_fetchcols($dbh, "SELECT chain_id, start_resno, end_resno FROM subsets_details WHERE subset_id = \"$subset_id\"") ;

   my $picks ;

# Iterate through the subset segments

   foreach my $j ( 0.. $#{$chain_id}) {

# If the chain is undefined, set it to ''

      if (!defined $chain_id->[$j]) {
         $chain_id->[$j] = ''; }

# If the start or end residue is undefined or '', set it to 'FIRST' and
#  'LAST', respectively.

      if ((!defined $start_resno->[$j]) || ($start_resno->[$j] eq '')) {
         $start_resno->[$j] = 'FIRST'; }

      if ((!defined $end_resno->[$j]) || ($end_resno->[$j] eq '')) {
         $end_resno->[$j] = 'LAST'; }

# MODELLER SEGEMENT SELECTION format: 'start_resno:chain_id' 'end_resno:chain_id'

      $picks->[$j] = "\'".$start_resno->[$j].':'.$chain_id->[$j]."\' ".
                     "\'".$end_resno->[$j].':'.$chain_id->[$j]."\'" ;
   }

   return $picks ;

}


=head2 subsetdef_2_mod_pick()

   Function:    Converts domain definition to MODELLER pick statements
   Args:        $_[0] = arrayref of chain_id
                $_[1] = arrayref of start_resno
                $_[2] = arrayref of end_resno
   Return:      $_->[] arrayref of MODELLER pick statements to select domain

=cut

sub subsetdef_2_mod_pick {

   my ( $chain_id, $start_resno, $end_resno ) = @_ ;

   my $picks ;

# Iterate through the subset segments

   foreach my $j ( 0.. $#{$chain_id}) {

# If the chain is undefined, set it to ''

      if (!defined $chain_id->[$j]) {
         $chain_id->[$j] = ''; }

# If the start or end residue is undefined or '', set it to 'FIRST' and 'LAST',
#  respectively.

      if ((!defined $start_resno->[$j]) || ($start_resno->[$j] eq '')) {
         $start_resno->[$j] = 'FIRST'; }

      if ((!defined $end_resno->[$j]) || ($end_resno->[$j] eq '')) {
         $end_resno->[$j] = 'LAST'; }

# MODELLER SEGEMENT SELECTION format: 'start_resno:chain_id' 'end_resno:chain_id'

      $picks->[$j] = "\'".$start_resno->[$j].':'.$chain_id->[$j]."\' ".
                     "\'".$end_resno->[$j].':'.$chain_id->[$j]."\'" ;
   }

   return $picks ;

}


=head2 get_salign (modeller_bin, bdp_file)

   Title:       get_salign()
   Function:    Calls MODELLER.SALIGN to structurally align two pdb files
   Args:        $_->{pdb_fn_1} - name of pdb file 1
                $_->{pdb_fn_2} - name of pdb file 2
                $_->{modeller_bin} - name of MODELLER binary file

=cut

sub get_salign {

# dont takes picks directly - no facility in MODELLER to pick discontiguous
# segments in both structures

   my $params = shift ;

   my $modeller_bin = $params->{modeller_bin} ;
   my $fn ;
   $fn->{pdb}->[0] = $params->{pdb_fn_1} ;
   $fn->{pdb}->[1] = $params->{pdb_fn_2} ;

# Specify the temporary alignment TOP file, and the output alignment file.

   my ($ali_top_fn, $ali_top_fh) ;
   my $temp_fh ;

   my ($fh) ;
   ($fh->{top}, $fn->{top}) =
      tempfile( "resequiv.align.XXXXXXX", SUFFIX => ".top") ;

# Generate the actual TOP file.


   my @afd ;
   foreach my $j ( 0 .. 1) {
      $fn->{pdbn}->[$j] = $fn->{pdb}->[$j] ; $fn->{pdbn}->[$j] =~ s/^.*\///g ;
      if ($fn->{pdb}->[$j] =~ /\//) {
         $fn->{pdbb}->[$j] = $fn->{pdb}->[$j] ;
         $fn->{pdbb}->[$j] =~ s/\/[^\/]+$//g ;
         push @afd, $fn->{pdbb}->[$j] ;
      }
   }

   my $t = pibase::timestamp() ;
   print {$fh->{top}} "# pibase::modeller::salign()\n#".$t."\n\n";
   print {$fh->{top}} "SET OUTPUT_CONTROL = 1 1 1 1 1\n" ;
   if ($#afd >= 0) {
      my $atom_file_dir  = ''; $atom_file_dir = join(':', @afd);
      print {$fh->{top}} "SET ATOM_FILES_DIRECTORY = \'$atom_file_dir\'\n\n" ; }

   print {$fh->{top}} "SET RMS_CUTOFFS = 3.5 3.5 60 60 15 60 60 60 60 60 60\n";
   print {$fh->{top}} "SET GAP_PENALTIES_3D = 0 3\n";
   print {$fh->{top}} "SET GAP_GAP_SCORE = 0, GAP_RESIDUE_SCORE = 0\n";
   print {$fh->{top}} "SET RR_FILE = \'\$(LIB)/as1.sim.mat\'\n" ;


   print {$fh->{top}} "READ_MODEL FILE = \'$fn->{pdbn}->[0]\'\n";
   print {$fh->{top}} "SEQUENCE_TO_ALI ALIGN_CODES = \'1_$fn->{pdbn}->[0]\',";
   print {$fh->{top}} " ATOM_FILES = \'$fn->{pdbn}->[0]\'\n" ;

   print {$fh->{top}} "READ_MODEL FILE = \'$fn->{pdbn}->[1]\'\n";
   print {$fh->{top}} "SEQUENCE_TO_ALI ADD_SEQUENCE = on, ";
   print {$fh->{top}} "ALIGN_CODES = ALIGN_CODES \'2_$fn->{pdbn}->[1]\', ";
   print {$fh->{top}} "ATOM_FILES = ATOM_FILES \'$fn->{pdbn}->[1]\'\n" ;

   print {$fh->{top}} "SET FEATURE_WEIGHTS = 1. 0. 0. 0. 1. 0., ";
   print {$fh->{top}} "GAP_PENALTIES_1D = -450 -50\n";
   print {$fh->{top}} "SALIGN OUTPUT = \'ALIGNMENT QUALITY\', IMPROVE_ALIGNMENT = on, FIT = on\n" ;


   print {$fh->{top}} "SET FEATURE_WEIGHTS = 1. .5 1. 1. 1. 0., ";
   print {$fh->{top}} "GAP_PENALTIES_1D = -450 -50\n";
   print {$fh->{top}} "SALIGN OUTPUT = \'ALIGNMENT QUALITY\', IMPROVE_ALIGNMENT = on, FIT = on\n" ;


   print {$fh->{top}} "SET WRITE_FIT = off, WRITE_WHOLE_PDB = off\n";
   print {$fh->{top}} "SET FEATURE_WEIGHTS = 1. 1. 1. 1. 1. 0., GAP_PENALTIES_1D = -450. -50.\n";
   print {$fh->{top}} "SALIGN OUTPUT = \'ALIGNMENT QUALITY\', IMPROVE_ALIGNMENT = on, FIT = on\n";
   print {$fh->{top}} "COMPARE\n" ;

   close ($fh->{top}) ;

# Specify the location of the MODELLER LOG file.

   $fn->{'log'} = $fn->{top} ; $fn->{'log'} =~ s/top$/log/ ;

# Run the TOP file through MODELLER.

   system("$modeller_bin $fn->{top} >/dev/null 2>&1") ;

   if (-s $fn->{'log'}) {
      return {ali_log => $fn->{'log'}, ali_top => $fn->{'top'}} ;
   } else {
      return {error => 'error'} ;
   }

}


=head2 get_salign_seqseq (modeller_bin, bdp_file)

   Title:       get_salign_seqseq()
   Function:    Calls MODELLER.SALIGN to sequence align two pdb files
   Args:        $_->{pdb_fn_1} - name of pdb file 1
                $_->{pdb_fn_2} - name of pdb file 2
                $_->{modeller_bin} - name of MODELLER binary file

=cut

sub get_salign_seqseq {

# dont takes picks directly - no facility in MODELLER to pick discontiguous
# segments in both structures

   my $params = shift ;

   my $modeller_bin = $params->{modeller_bin} ;
   my $fn ;
   $fn->{pdb}->[0] = $params->{pdb_fn_1} ;
   $fn->{pdb}->[1] = $params->{pdb_fn_2} ;

# Specify the temporary alignment TOP file, and the output alignment file.

   my ($ali_top_fn, $ali_top_fh) ;
   my $temp_fh ;

   my ($fh) ;
   ($fh->{top}, $fn->{top}) =
      tempfile( "resequiv.align.XXXXXXX", SUFFIX => ".top") ;

# Generate the actual TOP file.


   foreach my $j ( 0 .. 1) {
      $fn->{pdbn}->[$j] = $fn->{pdb}->[$j] ; $fn->{pdbn}->[$j] =~ s/^.*\///g ;
      $fn->{pdbb}->[$j] = $fn->{pdb}->[$j] ; $fn->{pdbb}->[$j] =~ s/\/[^\/]+$//g ;
   }

   my $t = pibase::timestamp() ;
   print {$fh->{top}} "# pibase::modeller::salign()\n#".$t."\n\n";
   print {$fh->{top}} "SET OUTPUT_CONTROL = 1 1 1 1 1\n" ;
   print {$fh->{top}} "SET ATOM_FILES_DIRECTORY = \'$fn->{pdbb}->[0]:" ;
   print {$fh->{top}} "$fn->{pdbb}->[1]\'\n\n" ;

   print {$fh->{top}} "SET RR_FILE = \'\$(LIB)/as1.sim.mat\'\n" ;


   print {$fh->{top}} "READ_MODEL FILE = \'$fn->{pdbn}->[0]\'\n";
   print {$fh->{top}} "SEQUENCE_TO_ALI ALIGN_CODES = \'1_$fn->{pdbn}->[0]\',";
   print {$fh->{top}} " ATOM_FILES = \'$fn->{pdbn}->[0]\'\n" ;

   print {$fh->{top}} "READ_MODEL FILE = \'$fn->{pdbn}->[1]\'\n";
   print {$fh->{top}} "SEQUENCE_TO_ALI ADD_SEQUENCE = on, ";
   print {$fh->{top}} "ALIGN_CODES = ALIGN_CODES \'2_$fn->{pdbn}->[1]\', ";
   print {$fh->{top}} "ATOM_FILES = ATOM_FILES \'$fn->{pdbn}->[1]\'\n" ;

   print {$fh->{top}} "SET FEATURE_WEIGHTS = 1. 0. 0. 0. 0. 0., ";
   print {$fh->{top}}    "GAP_PENALTIES_1D = -450 -50\n";
   print {$fh->{top}} "SET GAP_GAP_SCORE = 0, GAP_RESIDUE_SCORE = 0\n";
   print {$fh->{top}} "SET SIMILARITY_FLAG = \'on\'\n" ;
   print {$fh->{top}} "SET OVERHANG = 0\n" ;
   print {$fh->{top}} "SALIGN OUTPUT = \'\', IMPROVE_ALIGNMENT = on\n" ;


   print {$fh->{top}} "SET FEATURE_WEIGHTS = 1. .5 1. 1. 1. 0., ";
   print {$fh->{top}} "GAP_PENALTIES_1D = -450 -50\n";
   print {$fh->{top}} "SALIGN OUTPUT = \'ALIGNMENT QUALITY\', IMPROVE_ALIGNMENT = on, FIT = on\n" ;


   print {$fh->{top}} "SET WRITE_FIT = off, WRITE_WHOLE_PDB = off\n";
   print {$fh->{top}} "SET FEATURE_WEIGHTS = 1. 1. 1. 1. 1. 0., GAP_PENALTIES_1D = -450. -50.\n";
   print {$fh->{top}} "SALIGN OUTPUT = \'ALIGNMENT QUALITY\', IMPROVE_ALIGNMENT = on, FIT = on\n";

   close ($fh->{top}) ;
   print STDERR "$fn->{top}\n" ;

# Specify the location of the MODELLER LOG file.

   $fn->{'log'} = $fn->{top} ; $fn->{'log'} =~ s/top$/log/ ;

# Run the TOP file through MODELLER.

   system("$modeller_bin $fn->{top} >/dev/null 2>&1") ;

   if (-s $fn->{'log'}) {
      return {ali_log => $fn->{'log'}} ;
   } else {
      return {error => 'error'} ;
   }

}


=head2 OLD_modeller_subset_sasa (modeller_bin, bdp_file, picks)

   Title:       OLD_modeller_subset_sasa()
   Function:    Calculate the solvent accessible surface area of a pdb file
   Args:        $_[0] = modeller_binary location
                $_[1] = pdb file location
                $_[2] = pick statements for domain
   Returns:     $_[0] = SASA of domain (get_sasa() data structure)
                $_[1] = error_fl - error flag

=cut


sub OLD_modeller_subset_sasa {

   my ($modeller_bin, $bdp_file, $picks)  = @_ ;

   my $pick_init_pre = "PICK_ATOMS PICK_ATOMS_SET = 1, ATOM_TYPES = \'ALL\', SELECTION_STATUS = \'INITIALIZE\', SELECTION_SEGMENT = " ;
   my $pick_add_pre = "PICK_ATOMS PICK_ATOMS_SET = 1, ATOM_TYPES = \'ALL\', SELECTION_STATUS = \'ADD\', SELECTION_SEGMENT = " ;
   my $pick_rem_pre = "PICK_ATOMS PICK_ATOMS_SET = 1, ATOM_TYPES = \'ALL\', SELECTION_STATUS = \'REMOVE\', SELECTION_SEGMENT = " ;

=pod

Specify the temporary alignment TOP file, and the output alignment file.

=cut

   my $filename = $bdp_file ; $filename =~ s/^.*\///g ;
   my $filebase = $bdp_file ; $filebase =~ s/\/[^\/]+$//g ;

   my ($cutout_top_fh, $cutout_top_fn) =
      tempfile("$filename.cutup.XXXXXX", SUFFIX => ".top");

   my ($cutpdb_fh, $cutpdb_fn) =
      tempfile("$filename.subset.XXXXXX", SUFFIX =>".pdb") ;
   close ($cutpdb_fh) ;

=pod

Generate the actual TOP file.

=cut

   print $cutout_top_fh "# modeller_subset_sasa() ".localtime(time())."\n\n" ;
   if ($filename ne $bdp_file) {
      print $cutout_top_fh "SET ATOM_FILES_DIRECTORY = \'$filebase\'\n" ; }

   print $cutout_top_fh 'SET OUTPUT_CONTROL = 1 1 1 1 1'."\n" ;
   print $cutout_top_fh "READ_MODEL FILE = \'$filename\'\n" ;

   print $cutout_top_fh "SET SELECTION_SEARCH = \'SEGMENT\'\n" ;
   print $cutout_top_fh "SET RES_TYPES = \'ALL\'\n" ;

   print $cutout_top_fh $pick_init_pre.$picks->[0]."\n" ;
   foreach my $j ( 1 .. $#{$picks}) {
      print $cutout_top_fh $pick_add_pre.$picks->[$j]."\n" ; }
   print $cutout_top_fh "SET WRITE_ALL_ATOMS = off\n" ;
   print $cutout_top_fh "WRITE_MODEL FILE = \'".$cutpdb_fn."\'\n" ;

   close ( $cutout_top_fh ) ;


=pod

Specify the location of the MODELLER LOG file.

=cut

   my $cutout_log_fn = $cutout_top_fn ;
   $cutout_log_fn =~ s/top$/log/ ;

=pod

Run the TOP file through MODELLER.

=cut

   my ($cut_sasa, $error_fl) ;

   system("$modeller_bin $cutout_top_fn >/dev/null 2>&1") ;

   if (!-s $cutpdb_fn) {
      push @{$error_fl}, "cutpdb not generated, see: $cutpdb_fn $cutout_top_fn $cutout_log_fn" ;
   } else {
      unlink($cutout_top_fn, $cutout_log_fn);
   }

   if ($#{$error_fl} < 0) {
      my ($t1, $t2, $sasa_error_fl) ;

      ($t1, $t2, $cut_sasa, $sasa_error_fl) =
         get_sasa($cutpdb_fn, $modeller_bin);

      foreach my $j ( 0 .. $#{$sasa_error_fl}) {
         push @{$error_fl}, "get_sasa() $sasa_error_fl->[$j]" ; }
      unlink($cutpdb_fn) ;
   }

   return ($cut_sasa, $error_fl) ;

}


=head2 get_sasa(modeller_bin, bdp_file, picks)

   Title:       get_sasa()
   Function:    parse modeller sasa (psa)
   Args:        $_->{surftyp} = type of MODELLER surface area
                $_->{pdb_fn} = pdb file location
                $_->{modeller_bin} = modeller binary file
   Returns:     $_[0] = results
                $_[1] = resno_rev
                $_[2] = sasa
                $_[3] = error_fl

=cut

sub OLD_get_sasa {

   my $in = shift ;

   my ($bdp_file, $modeller_bin) ;
   my $surftyp = '1' ; #default contact

   if (ref($in) ne '') {
      if (exists $in->{surftyp}) {
         $surftyp = $in->{surftyp} ; }

      $bdp_file = $in->{pdb_fn} ;
      $modeller_bin = $in->{modeller_bin} ;
   } else {
      $bdp_file = $in ;
      $modeller_bin = shift ;
   }

# Specify MODELLER temporary TOP file, and output files.

   my $filesmade ;

   my ($sasa_top_fh, $sasa_top_fn) = 
      tempfile("temp_sasa.XXXXXX", SUFFIX=>".top") ;
   $filesmade->{$sasa_top_fn}++ ;

   my $tempbase = $sasa_top_fn ; $tempbase =~ s/\.top$// ;
   my $sasa_log_fn = $tempbase.".log" ;
   my $sasa_atm_fn = $tempbase.".sol" ;
   my $sasa_res_fn = $tempbase.".psa" ;

   $filesmade->{$sasa_log_fn}++ ;
   $filesmade->{$sasa_atm_fn}++ ;
   $filesmade->{$sasa_res_fn}++ ;

# Write the MODELLER TOP file to calcaulte solvent accessibilty.

   my $filename = $bdp_file ; $filename =~ s/^.*\///g ;
   my $filebase = $bdp_file ; $filebase =~ s/\/[^\/]+$//g ;

   if ($filename ne $bdp_file) {
      print $sasa_top_fh "SET ATOM_FILES_DIRECTORY = \'$filebase\'\n" ; }

   print $sasa_top_fh 'READ_TOPOLOGY FILE = \'$(LIB)/top_heav.lib\''."\n" ;
   print $sasa_top_fh 'SET HETATM_IO = off, WATER_IO=off'."\n" ;
   print $sasa_top_fh 'READ_MODEL FILE = \''.$filename.'\''."\n" ;
   print $sasa_top_fh 'SET SURFTYP = '."$surftyp\n" ;
   print $sasa_top_fh 'SET RADII_FACTOR = 1.0'."\n" ;
   print $sasa_top_fh 'WRITE_DATA FILE = \''.$tempbase.'\', OUTPUT = \'PSA\''."\n" ;
   close($sasa_top_fh) ;

# Run the MODELLER TOP file.

   my $tcom = "$modeller_bin $sasa_top_fn >/dev/null 2>&1" ;
   system($tcom) ;

# Initialize solvent accessibility storage arrays

   my $results ;
   my $resno_rev ;
   my %fields = (
      'resno' => [7, 4] ,
      'resna' => [14, 3] ,
      'chain' => [18, 1],
      'all_sum' => [19, 7],
      'all_perc' => [27, 5],
      'nonp_sum' => [33, 7],
      'nonp_perc' => [41, 5],
      'p_sum' => [47, 7],
      'p_perc' => [55, 5],
      'sc_sum' => [61, 7],
      'sc_perc' => [69, 5],
      'mc_sum' => [75, 7],
      'mc_perc' => [83, 5]
   ) ;

   foreach my $key (keys %fields) {
      $results->{$key} = [] ; }

# If the MODELLER PSA file does not exist, display to STDERR and return.

   if (!-s $sasa_res_fn) {
      $results->{resno}->[0] = 'ERROR' ;
      $resno_rev->{ERROR} = 'ERROR' ;
      my $err = ['DANGER WILL ROBINSON'] ;
      foreach my $tfn (keys %{$filesmade}) {if (-e $tfn) {unlink $tfn;}}
      return ($results, $resno_rev, undef, $err);
   }

# Open the resulting MODELLER PSA file.

   open(PSA, $sasa_res_fn) ;

# Read in a line of the PSA file.

   my $sasa ;
   $sasa->{all} = 0 ;
   $sasa->{mc} = 0 ;
   $sasa->{sc} = 0 ;
   $sasa->{nonp} = 0 ;
   $sasa->{p} = 0 ;

   while (my $line = <PSA>) {
      chomp $line;

# If the line contains accessibility information,

      if ($line =~ /^ACCESS/) {

# Extract the residue information and solvent accessibility parameters

         foreach my $key (keys %fields) {

            my $t = substr($line, $fields{$key}->[0], $fields{$key}->[1]);
            unless ($key eq 'chain') {
               $t =~ s/ //g ; }

	    if ($key eq 'all_sum') {
	       $sasa->{all} += $t;
	    } elsif ($key eq 'mc_sum') {
	       $sasa->{mc} += $t;
	    } elsif ($key eq 'sc_sum') {
	       $sasa->{sc} += $t;
	    } elsif ($key eq 'nonp_sum') {
	       $sasa->{nonp} += $t;
	    } elsif ($key eq 'p_sum') {
	       $sasa->{p} += $t;
	    }

            push @{$results->{$key}}, $t ;
         }

# Create an entry in the reverse lookup hash pointing from the actual
#  PDB residue number and chain id to the index of the arrays where the
#  solvent accessibility values are stored.

         my $cur_recno = $#{$results->{resno}} ;
         my $res_sig = $results->{resno}->[$cur_recno]."_".
                       $results->{chain}->[$cur_recno] ;

         $resno_rev->{$res_sig} = $cur_recno ;

      }
   }

# Close the PSA file.

   close PSA ;

# If no accessibility lines have been read in, set error flag

   my $error_fl ;

   if ($#{$results->{resno}} < 0 ) {
      push @{$error_fl}, "no sasa entries calculated" ; }

# Remove the MODELLER TOP file and output files.

# Return pointers to teh arrays holding the values and the residue lookup hash.

# Note: accessibility values for ALL Cys read are retunred; However, the calling function will only ask for Cys involved in disulfide bridges.

   return ( $results, $resno_rev, $sasa, $error_fl ) ;

}


=head2 calc_sasa()

   Title:       calc_sasa()
   Function:    Run and parse modeller sasa (psa)
   Args:        $_->{surftyp} = type of MODELLER surface area
                $_->{pdb_fn} = pdb file location
                $_->{modeller_bin} = modeller binary file
   Returns:     ->{res_sasa}->{all_sum|mc_sum|sc_sum|p_sum|nonp_sum}->[i] =
                    SASA information for residue record i
                ->{sasa_resno_rev}->{"residuenumber_chain"} = record number
                ->{full_sasa}->{p|nonp|mc|sc|all} = totals of residue sasa records
                ->{error_fl} => $error_fl,
                ->{atm_sasa}->{p|nonp|mc|sc|all} = totals of ATOM sasa records

=cut

sub calc_sasa {

   my $in = shift ;

   my ($bdp_file, $modeller_bin) ;
   my $surftyp = '1' ; #default contact

   if (ref($in) ne '') {
      if (exists $in->{surftyp}) {
         $surftyp = $in->{surftyp} ; }

      $bdp_file = $in->{pdb_fn} ;
      $modeller_bin = $in->{modeller_bin} ;
   } else {
      $bdp_file = $in ;
      $modeller_bin = shift ;
   }

#Specify MODELLER temporary TOP file, and output files.

   my $filesmade ;

   my ($sasa_top_fh, $sasa_top_fn) = 
      tempfile("temp_sasa.XXXXXX", SUFFIX=>".top") ;
   $filesmade->{$sasa_top_fn}++ ;

   my $tempbase = $sasa_top_fn ; $tempbase =~ s/\.top$// ;
   my $sasa_log_fn = $tempbase.".log" ;
   my $sasa_atm_fn = $tempbase.".sol" ;
   my $sasa_res_fn = $tempbase.".psa" ;

   $filesmade->{$sasa_log_fn}++ ;
   $filesmade->{$sasa_atm_fn}++ ;
   $filesmade->{$sasa_res_fn}++ ;

#Write the MODELLER TOP file to calcaulte solvent accessibilty.

   my $filename = $bdp_file ; $filename =~ s/^.*\///g ;
   my $filebase = $bdp_file ; $filebase =~ s/\/[^\/]+$//g ;

   if ($filename ne $bdp_file) {
      print $sasa_top_fh "SET ATOM_FILES_DIRECTORY = \'$filebase\'\n" ; }

   print $sasa_top_fh 'READ_TOPOLOGY FILE = \'$(LIB)/top_heav.lib\''."\n" ;
   print $sasa_top_fh 'SET HETATM_IO = off, WATER_IO=off'."\n" ;
   print $sasa_top_fh 'READ_MODEL FILE = \''.$filename.'\''."\n" ;
   print $sasa_top_fh 'SET SURFTYP = '."$surftyp\n" ;
   print $sasa_top_fh 'SET RADII_FACTOR = 1.0'."\n" ;
   print $sasa_top_fh 'WRITE_DATA FILE = \''.$tempbase.'\', OUTPUT = \'PSA\''."\n" ;
   close($sasa_top_fh) ;

#Run the MODELLER TOP file.

   my $tcom = "$modeller_bin $sasa_top_fn >/dev/null 2>&1" ;
   system($tcom) ;

#Initialize solvent accessibility storage arrays

   my $results ;
   my $resno_rev ;
   my %fields = (
      'resno' => [7, 4] ,
      'resna' => [14, 3] ,
      'chain' => [18, 1],
      'all_sum' => [19, 7],
      'all_perc' => [27, 5],
      'nonp_sum' => [33, 7],
      'nonp_perc' => [41, 5],
      'p_sum' => [47, 7],
      'p_perc' => [55, 5],
      'sc_sum' => [61, 7],
      'sc_perc' => [69, 5],
      'mc_sum' => [75, 7],
      'mc_perc' => [83, 5]
   ) ;

   foreach my $key (keys %fields) {
      $results->{$key} = [] ; }

#If the MODELLER PSA file does not exist, display to STDERR and return.

   if (!-s $sasa_res_fn) {
      $results->{resno}->[0] = 'ERROR' ;
      $resno_rev->{ERROR} = 'ERROR' ;
      my $err = ['DANGER WILL ROBINSON'] ;
      foreach my $tfn (keys %{$filesmade}) {if (-e $tfn) {unlink $tfn;}}
      return {
         res_sasa => $results,
         sasa_resno_rev => $resno_rev,
         error_fl => $err,
      } ;
   }

#Open the resulting MODELLER PSA file.


#Read in a line of the PSA file.

   my $sasaatoms ;
   $sasaatoms->{all} = 0 ;
   $sasaatoms->{p} = 0 ;
   $sasaatoms->{nonp} = 0 ;

   open(SOL, $sasa_atm_fn) ;
   while (my $line = <SOL>) {
      chomp $line;
      if ($line !~ /^ATOM/) {next;}

      my $atomtype = substr($line,12,1) ;
      my $atomacc = substr($line,64,8) ;
      $atomacc =~ s/ //g ;
      $sasaatoms->{all} += $atomacc ;
      if ($atomtype =~ /[NO]/) {
         $sasaatoms->{p} += $atomacc ;
      } else {
         $sasaatoms->{nonp} += $atomacc ;
      }
   }
   close(SOL) ;


   my $sasa ;
   $sasa->{all} = 0 ;
   $sasa->{mc} = 0 ;
   $sasa->{sc} = 0 ;
   $sasa->{nonp} = 0 ;
   $sasa->{p} = 0 ;
   open(PSA, $sasa_res_fn) ;
   while (my $line = <PSA>) {

      chomp $line;

#If the line contains accessibility information,

      if ($line =~ /^ACCESS/) {

#Extract residue information and solvent accessibility parameters from the line.

         foreach my $key (keys %fields) {

            my $t = substr($line, $fields{$key}->[0], $fields{$key}->[1]);
            unless ($key eq 'chain') {
               $t =~ s/ //g ; }

	    if ($key eq 'all_sum') {
	       $sasa->{all} += $t;
	    } elsif ($key eq 'mc_sum') {
	       $sasa->{mc} += $t;
	    } elsif ($key eq 'sc_sum') {
	       $sasa->{sc} += $t;
	    } elsif ($key eq 'nonp_sum') {
	       $sasa->{nonp} += $t;
	    } elsif ($key eq 'p_sum') {
	       $sasa->{p} += $t;
	    }

            push @{$results->{$key}}, $t ;
         }

#Create an entry in the reverse lookup hash pointing from the actual PDB residue number and chain id to the index of the arrays where the solvent accessibility values are stored.

         my $cur_recno = $#{$results->{resno}} ;
         my $res_sig = $results->{resno}->[$cur_recno]."_".
                       $results->{chain}->[$cur_recno] ;

         $resno_rev->{$res_sig} = $cur_recno ;

      }

   }

#Close the PSA file.

   close PSA ;

#If no accessibility lines have been read in, set error flag

   my $error_fl = [];

   if ($#{$results->{resno}} < 0 ) {
      push @{$error_fl}, "no sasa entries calculated" ; }

#Remove the MODELLER TOP file and output files.

   unlink($sasa_top_fn, $sasa_log_fn, $sasa_atm_fn, $sasa_res_fn) ;

#Return pointers to teh arrays holding the values and the residue lookup hash.
#
#Note: accessibility values for ALL Cys read are retunred; However, the calling function will only ask for Cys involved in disulfide bridges.
#

   return {
      res_sasa => $results,
      sasa_resno_rev => $resno_rev,
      full_sasa => $sasa,
      error_fl => $error_fl,
      atm_sasa => $sasaatoms
   } ;
}



=head2 get_dihedrals()

   Title:       get_dihedrals()
   Function:    run and parse modeller dihedrals (dih)
   Args:        $_[0] = bdp_file
                $_[1] = modeller_bin
   Returns:     $_[0] = results
                $_[1] = resno_rev
                $_[2] = error_fl

=cut

sub get_dihedrals {

   my $bdp_file = shift ;
   my $modeller_bin = shift ;

# Specify MODELLER temporary TOP file, and output files.

   my ($dih_top_fh, $dih_top_fn) = 
      tempfile("temp_dih.XXXXXX", SUFFIX=>".top") ;

   my $tempbase = $dih_top_fn ; $tempbase =~ s/\.top$// ;
   my $dih_log_fn = $tempbase.".log" ;
   my $dih_res_fn = $tempbase.".dih" ;

# Write the MODELLER TOP file to calcaulte solvent accessibilty.

   my $filename = $bdp_file ; $filename =~ s/^.*\///g ;
   my $filebase = $bdp_file ; $filebase =~ s/\/[^\/]+$//g ;

   if ($filename ne $bdp_file) {
      print $dih_top_fh "SET ATOM_FILES_DIRECTORY = \'$filebase\'\n" ; }

   print $dih_top_fh 'READ_MODEL FILE = \''.$filename.'\''."\n" ;
   print $dih_top_fh 'READ_TOPOLOGY FILE = \'$(LIB)/top_heav.lib\''."\n" ;
   print $dih_top_fh 'SET HETATM_IO = off, WATER_IO=off'."\n" ;
   print $dih_top_fh 'IUPAC_MODEL'."\n" ;
   print $dih_top_fh 'WRITE_DATA FILE = \''.$tempbase.'\', OUTPUT = \'DIH\''."\n" ;
   close($dih_top_fh) ;

# Run the MODELLER TOP file.

   my $tcom = "$modeller_bin $dih_top_fn >/dev/null 2>&1" ;
   system($tcom) ;

# Initialize solvent accessibility storage arrays

   my $results ;
   my $resno_rev ;
   my %fields = (
      'resno' => [0, 5] , #does splitting on 4 cut off the resno in PSA? BUG?
      'resna' => [7, 3] , #does splitting on 4 cut off the resno in PSA? BUG?
      'chain' => [12, 1],
      'alpha' => [13, 7],
      'phi' => [20, 7],
      'psi' => [27, 7],
      'omega' => [34, 7],
      'chi_1' => [41, 7],
      'chi_2' => [48, 7],
      'chi_3' => [55, 7],
      'chi_4' => [62, 7],
      'chi_5' => [69, 7],
   ) ;

   foreach my $key (keys %fields) {
      $results->{$key} = [] ; }

# If the MODELLER PSA file does not exist, display to STDERR and return.

   if (!(-s $dih_res_fn)) {
      $results->{resno}->[0] = 'ERROR' ;
      $resno_rev->{ERROR} = 'ERROR' ;
      return ($results, $resno_rev, ['modeller dihedral error']) ;
   }

# Open the resulting MODELLER DIH file.

   open(DIH, $dih_res_fn) ;

# Read in a line of the DIH file.

   my $dih;

   my $imin = 0 ;
   while (my $line = <DIH>) {

      chomp $line;

# If the line contains accessibility information,

      if ($imin) {

# Extract the residue information and solvent accessibility parameters from the line.

	 my $curval ;
         foreach my $key (keys %fields) {

            my $t = substr($line, $fields{$key}->[0], $fields{$key}->[1]);
            unless ($key eq 'chain') {
               $t =~ s/ //g ; }

	    if ($key eq 'resno') {
	       if ($t eq 'HSD') { $t = 'HIS'; }
	       if ($t eq 'HSE') { $t = 'HIS'; }
	    }
            push @{$results->{$key}}, $t ;
	    $curval->{$key} = $t ;
         }

# Create an entry in the reverse lookup hash pointing from the actual PDB residue number and chain id to the index of the arrays where the solvent accessibility values are stored.

         my $cur_recno = $#{$results->{resno}} ;
         my $res_sig = $results->{resno}->[$cur_recno]."\n".
                       $results->{chain}->[$cur_recno] ;

         $resno_rev->{$res_sig}->{resna} = $curval->{resna} ;
         $resno_rev->{$res_sig}->{alpha} = $curval->{alpha} ;
         $resno_rev->{$res_sig}->{phi} = $curval->{phi} ;
         $resno_rev->{$res_sig}->{psi} = $curval->{psi} ;
         $resno_rev->{$res_sig}->{omega} = $curval->{omega} ;
         $resno_rev->{$res_sig}->{chi_1} = $curval->{chi_1} ;
         $resno_rev->{$res_sig}->{chi_2} = $curval->{chi_2} ;
         $resno_rev->{$res_sig}->{chi_3} = $curval->{chi_3} ;
         $resno_rev->{$res_sig}->{chi_4} = $curval->{chi_4} ;
         $resno_rev->{$res_sig}->{chi_5} = $curval->{chi_5} ;

      } elsif ($line =~ /ALPHA/) {
         $imin = 1 ;
      }

   }

# Close the PSA file.

   close DIH ;

# If no dihedrals have been read in, set error flag

   my $error_fl ;

   if ($#{$results->{resno}} < 0 ) {
      push @{$error_fl}, "no dihedral entries calculated" ; }

# Remove the MODELLER TOP file and output files.

   unlink($dih_top_fn, $dih_log_fn, $dih_res_fn) ;

# Return pointers to teh arrays holding the values and the residue lookup hash.
# Note: accessibility values for ALL Cys read are retunred; However, the calling function will only ask for Cys involved in disulfide bridges.

   return ( $results, $resno_rev, $error_fl ) ;

}


=head2 get_vol()

   Function:    Run and parse MODELLER volume (psa)
   Args:        $_[0] = bdp_file
                $_[1] = MODELLER binary file
   Returns:     $_[0] = results
                $_[1] = resno_rev
                $_[2] = sasa
                $_[3] = error_fl


=cut

sub get_vol {

   my $bdp_file = shift ;
   my $modeller_bin = shift ;

# Specify MODELLER temporary TOP file, and output files.

   my ($vol_top_fh, $vol_top_fn) = 
      tempfile("temp_vol.XXXXXX", SUFFIX=>".top") ;

   my $tempbase = $vol_top_fn ; $tempbase =~ s/\.top$// ;
   my $vol_log_fn = $tempbase.".log" ;
   my $vol_res_fn = $tempbase.".cav" ;

# Write the MODELLER TOP file to calcaulte solvent accessibilty.

   my $filename = $bdp_file ; $filename =~ s/^.*\///g ;
   my $filebase = $bdp_file ; $filebase =~ s/\/[^\/]+$//g ;

   if ($filename ne $bdp_file) {
      print $vol_top_fh "SET ATOM_FILES_DIRECTORY = \'$filebase\'\n" ; }

   print $vol_top_fh 'READ_MODEL FILE = \''.$filename.'\''."\n" ;
   print $vol_top_fh 'READ_TOPOLOGY FILE = \'$(LIB)/top_heav.lib\''."\n" ;
   print $vol_top_fh 'SET HETATM_IO = off, WATER_IO=off'."\n" ;
   print $vol_top_fh 'SET RADII_FACTOR = 1.0'."\n" ;
   print $vol_top_fh 'WRITE_DATA FILE = \''.$tempbase.'\', OUTPUT = \'CAV\''."\n" ;
   close($vol_top_fh) ;

# Run the MODELLER TOP file.

   my $tcom = "$modeller_bin $vol_top_fn >/dev/null 2>&1" ;
   system($tcom) ;
   die("bye bye; vol done\n") ;

# Initialize solvent accessibility storage arrays

   my $results ;
   my $resno_rev ;
   my %fields = (
      'resno' => [7, 4] ,
      'chain' => [18, 1],
      'all_sum' => [19, 7],
      'all_perc' => [27, 5],
      'nonp_sum' => [33, 7],
      'nonp_perc' => [41, 5],
      'p_sum' => [47, 7],
      'p_perc' => [55, 5],
      'sc_sum' => [61, 7],
      'sc_perc' => [69, 5],
      'mc_sum' => [75, 7],
      'mc_perc' => [83, 5]
   ) ;

   foreach my $key (keys %fields) {
      $results->{$key} = [] ; }

# If the MODELLER PSA file does not exist, display to STDERR and return.

   if (!(-s $vol_res_fn)) {
      $results->{resno}->[0] = 'ERROR' ;
      $resno_rev->{ERROR} = 'ERROR' ;
      return ($results, $resno_rev);
   }

# Open the resulting MODELLER PSA file.

   open(PSA, $vol_res_fn) ;

# Read in a line of the PSA file.

   my $vol ;
   $vol->{all} = 0 ;
   $vol->{nonp} = 0 ;
   $vol->{p} = 0 ;

   while (my $line = <PSA>) {

      chomp $line;

# If the line contains accessibility information,

      if ($line =~ /^ACCESS/) {

# Extract the residue information and solvent accessibility parameters from the line.

         foreach my $key (keys %fields) {

            my $t = substr($line, $fields{$key}->[0], $fields{$key}->[1]);
            unless ($key eq 'chain') {
               $t =~ s/ //g ; }

	    if ($key eq 'all_sum') {
	       $vol->{all} += $t;
	    } elsif ($key eq 'nonp_sum') {
	       $vol->{nonp} += $t;
	    } elsif ($key eq 'p_sum') {
	       $vol->{p} += $t;
	    }

            push @{$results->{$key}}, $t ;
         }

# Create an entry in the reverse lookup hash pointing from the actual PDB residue number and chain id to the index of the arrays where the solvent accessibility values are stored.

         my $cur_recno = $#{$results->{resno}} ;
         my $res_sig = $results->{resno}->[$cur_recno]."_".
                       $results->{chain}->[$cur_recno] ;

         $resno_rev->{$res_sig} = $cur_recno ;

      }

   }

# Close the PSA file.

   close PSA ;

# If no accessibility lines have been read in, set error flag

   my $error_fl ;

   if ($#{$results->{resno}} < 0 ) {
      push @{$error_fl}, "no vol entries calculated" ; }

# Remove the MODELLER TOP file and output files.

   unlink($vol_top_fn, $vol_log_fn, $vol_res_fn) ;

# Return pointers to teh arrays holding the values and the residue lookup hash.
# Note: accessibility values for ALL Cys read are retunred; However, the calling function will only ask for Cys involved in disulfide bridges.

   return ( $results, $resno_rev, $vol, $error_fl ) ;

}


=head2 cutpdb()

   Title:       cutpdb()
   Function:    Uses MODELLER to extrct domain from a pdb file
   Args:        $_[0] = MODELLER binary location
                $_[1] = pdb file location
                $_[2] = MODELLER pick statements
                $_[3] = output pdb file name

   Returns:     $_[0] = 
                $_[1] = resno_rev
                $_[2] = error_fl

   FILE in:     pdb file ($_[1])
   FILE out:    domain pdb file ($_[3])

=cut

sub cutpdb {

   my ($modeller_bin, $bdp_file, $picks, $cutpdb_fn)  = @_ ;

   my $pick_init_pre = "PICK_ATOMS PICK_ATOMS_SET = 1, ATOM_TYPES = \'ALL\', SELECTION_STATUS = \'INITIALIZE\', SELECTION_SEGMENT = " ;
   my $pick_add_pre = "PICK_ATOMS PICK_ATOMS_SET = 1, ATOM_TYPES = \'ALL\', SELECTION_STATUS = \'ADD\', SELECTION_SEGMENT = " ;
   my $pick_rem_pre = "PICK_ATOMS PICK_ATOMS_SET = 1, ATOM_TYPES = \'ALL\', SELECTION_STATUS = \'REMOVE\', SELECTION_SEGMENT = " ;

# Specify the temporary alignment TOP file, and the output alignment file.

   my $filename = $bdp_file ; $filename =~ s/^.*\///g ;
   my $filebase = $bdp_file ; $filebase =~ s/\/[^\/]+$//g ;

   my ($cutout_top_fh, $cutout_top_fn) =
      tempfile("$filename.cutup.XXXXXX", SUFFIX => ".top") ;

# Generate the actual TOP file.

   print $cutout_top_fh "# cutpdb() ".localtime(time())."\n\n" ;
   if ($filename ne $bdp_file) {
      print $cutout_top_fh "SET ATOM_FILES_DIRECTORY = \'$filebase\'\n" ; }

   print $cutout_top_fh 'SET OUTPUT_CONTROL = \'11111\''."\n" ;
   print $cutout_top_fh "READ_MODEL FILE = \'$filename\'\n" ;

   print $cutout_top_fh "SET SELECTION_SEARCH = \'SEGMENT\'\n" ;
   print $cutout_top_fh "SET RES_TYPES = \'ALL\'\n" ;

   print $cutout_top_fh $pick_init_pre.$picks->[0]."\n" ;
   foreach my $j ( 1 .. $#{$picks}) {
      print $cutout_top_fh $pick_add_pre.$picks->[$j]."\n" ; }
   print $cutout_top_fh "SET WRITE_ALL_ATOMS = off\n" ;
   print $cutout_top_fh "WRITE_MODEL FILE = \'".$cutpdb_fn."\'\n" ;

   close ( $cutout_top_fh ) ;

# Specify the location of the MODELLER LOG file.

   my $cutout_log_fn = $cutout_top_fn ; $cutout_log_fn =~ s/top$/log/g ;

# Run the TOP file through MODELLER.

   system("$modeller_bin $cutout_top_fn >/dev/null 2>&1") ;
   unlink ($cutout_top_fn, $cutout_log_fn) ;

}


=head2 parse_ali()

   Title:       parse_ali()
   Function:    reads in a modeller PIR format alignment and returns residue
                number equivalence hashes
   Args:        $_->{ali_fn} alignment file
                $_->{modpipe_newstyle_orderswitch}
                  - 1 (Default) if new style MODPIPE run
                  - reordered sequences in the alignment file
   Results:     ->{seq} = $seq ;
                ->{resno_start} = $resno_start ;
                ->{resno_end} = $resno_end ;
                ->{chain_start} = $chain_start ;
                ->{chain_end} = $chain_end ;
                ->{alipos_2_serial} = $alipos_2_serresno ;
                ->{alipos_2_chainno} = $alipos_2_chainno ;
                ->{alipos_2_resna} = $alipos_2_resna ;
                ->{maxlength} = $maxlength;

=cut

sub parse_ali {

   my $params = shift ;
   my $ali_fn = $params->{ali_fn} ;

   if ((!defined $ali_fn) || (!-e $ali_fn)) {
      return {error => ['PDB file not specified or does not exist']} ; }

   open (ALIF, $ali_fn) ;
   my $headers ;
   my $seq ;
   my $cur_seq = -1 ;

   while ( my $line = <ALIF> ) {

      chomp $line;

      if (($line =~ /^structure/) || ($line =~ /^sequence/)) {

         if (exists $params->{modpipe_newstyle_orderswitch} &&
            $params->{modpipe_newstyle_orderswitch} == 1) {

            if ($line =~ /^structure/) { $cur_seq = 0 ; }
            elsif ($line =~ /^sequence/) { $cur_seq = 1 ; }

         } else {
	    $cur_seq++ ;
         }

	 $seq->[$cur_seq] = '' ;
         $headers->[$cur_seq] = $line ;

      } elsif ( ($line !~ /^$/) && ($line !~ /^\>/) &&
                ($line !~ /^C;/) && ($line !~ /^R;/) ) {

         $seq->[$cur_seq] .= $line ;

      }
   }
   close(ALIF) ;

   foreach my $j ( 0 .. $#{$seq}) {
      $seq->[$j] =~ s/\*$// ; }

   my ( $resno_start, $resno_end, $chain_start, $chain_end ) ;
   foreach my $j ( 0 .. $#{$headers} ) {
      my @t = split(':', $headers->[$j]) ;

      $t[2] =~ s/ //g ; $t[4] =~ s/ //g ;
      $t[2] =~ s/\.//g ; $t[4] =~ s/\.//g ;
      $t[3] =~ s/\.//g ; $t[5] =~ s/\.//g ;

      $resno_start->[$j] = $t[2] ;
      $chain_start->[$j] = $t[3] ;

#      print STDERR " from $headers->[$j]:\n" ;
#      print STDERR " modeller.pm: entry $j: STARTER IS $t[2] on $t[3]\n" ;
#      print STDERR " modeller.pm: entry $j: END     IS $t[4] on $t[5]\n" ;

      $resno_end->[$j] = $t[4] ;
      $chain_end->[$j] = $t[5] ;
   }

   my $alipos_2_serresno ;
   my $alipos_2_chainno ;
   my $alipos_2_resna ;
   my $maxlength = length($seq->[0]) ;

   foreach my $j ( 0 .. $#{$seq} ) {
      my $curres = 1 ;
      my $curchain = 1 ;
      my $curlength = length($seq->[$j]) ;
      if ($curlength > $maxlength) {$maxlength = $curlength} ;
      foreach my $k ( 0 .. ($curlength - 1) ) {
         my $curchar = substr($seq->[$j], $k, 1) ;
	 if ( $curchar eq "\\") {
	    $curres = 1 ;
	    $curchain++ ;
	 } elsif ( $curchar ne '-' ) {
            $alipos_2_serresno->[$j]->[$k] = $curres;
	    $alipos_2_chainno->[$j]->[$k] = $curchain ;
	    $alipos_2_resna->[$j]->[$k] = $curchar;
            $curres++ ;
         }

      }
   }

   my $results ;
   $results->{seq} = $seq ;
   $results->{resno_start} = $resno_start ;
   $results->{resno_end} = $resno_end ;
   $results->{chain_start} = $chain_start ;
   $results->{chain_end} = $chain_end ;
   $results->{alipos_2_serial} = $alipos_2_serresno ;
   $results->{alipos_2_chainno} = $alipos_2_chainno ;
   $results->{alipos_2_resna} = $alipos_2_resna ;
   $results->{maxlength} = $maxlength;
   return $results ;

}


=head2 get_resequiv()

   Title:       get_resequiv()
   Function:    Determines a mapping between residues in two pdb files.
                  Returns a combine serial residue number/positioning from
                  get_resequiv_serial with residue_info()
   Args:        ->{modeller_bin} = MODELLER binary location
                ->{pdb_fn_1} = PDB file 1 location
                ->{pdb_fn_2} = PDB file 2 location
   Returns:     $->[0]->{resno1} = resno2. maps from resno1 in first pdb file
                  to the aligned residue in the second pdb file
                $->[1]->{resno2} = resno1. maps from resno2 in second pdb file
                  to the aligned residue in the first pdb file

=cut

sub get_resequiv {


   my $params = shift ;

   my $modeller_bin = $params->{modeller_bin} ;
   my $pdb_fn ;
   $pdb_fn->[0] = $params->{pdb_fn_1} ;
   $pdb_fn->[1] = $params->{pdb_fn_2} ;

   my $re_serial = get_resequiv_serial($params) ;

   my $resinfo ;
   ($resinfo->[0]) = residue_info({ pdb_fn => $pdb_fn->[0]}) ;
   ($resinfo->[1]) = residue_info({ pdb_fn => $pdb_fn->[1]}) ;

   my $serial_2_resno ;
   my $resno_2_serial ;
   my $chainno_2_id ;
   my $chainid_2_no ;
   foreach my $j ( 0 .. $#{$resinfo}) {
      foreach my $k ( 0 .. $#{$resinfo->[$j]->{resno_serial}} ) {
	 my $curserch = $resinfo->[$j]->{chain_no}->[$k] ;
	 my $curserres = $resinfo->[$j]->{resno_serial}->[$k];
	 my $curresno = $resinfo->[$j]->{resno}->[$k];
	 my $curch = $resinfo->[$j]->{chain_id}->[$k];

	 my $sig1 = $curserres."\n".$curserch ;
	 my $sig2 = $curresno."\n".$curch;

         $serial_2_resno->[$j]->{$sig1} = $sig2 ;
         $resno_2_serial->[$j]->{$sig2} = $sig1 ;

	 $chainid_2_no->[$j]->{$curch} = $curserch ;
	 $chainno_2_id->[$j]->{$curserch} = $curch ;
      }
   }


   my $resequiv ;
   foreach my $pos ( 0 .. ($re_serial->{maxlength} - 1)) {
      if ( (exists $re_serial->{alipos_2_serial}->[0]->[$pos]) &&
           (exists $re_serial->{alipos_2_serial}->[1]->[$pos]) ) {

	 my $real_chainid ;
	 my $real_resno ;
	 foreach my $j ( 0 .. $#{$resinfo}) {

	    my $cur_serialresno = $re_serial->{alipos_2_serial}->[$j]->[$pos] ;

	    my $cur_chainno = $re_serial->{alipos_2_chainno}->[$j]->[$pos] ;

	    if ($re_serial->{alipos_2_chainno}->[$j] == 1) {
	       my $t_startserres = $resno_2_serial->[$j]->{$re_serial->{resno_start}->[$j]."\n".$re_serial->{chain_start}} ;
	       my $t_startserch = $chainid_2_no->[$j]->{$re_serial->{chain_start}->[$j]} ;
	       $cur_serialresno += $t_startserres ;
	       $cur_chainno += $t_startserch ;
            }
            
            $real_chainid->[$j] = $chainno_2_id->[$j]->{$cur_chainno} ;
	    $real_resno->[$j] = $serial_2_resno->[$j]->{$cur_serialresno."\n".$cur_chainno} ;

	 }
         $resequiv->[0]->{$real_resno->[0]} = $real_resno->[1] ;
         $resequiv->[1]->{$real_resno->[1]} = $real_resno->[0] ;

      }
   }

   return $resequiv ;

}


=head2 get_resequiv_serial()

   Title:       get_resequiv_serial()
   Function:    Determines residue equivalencies between two pdb files by
                  aligning them structurally using MODELLER.SALIGN
   Args:        ->{modeller_bin} = MODELLER binary location
                ->{pdb_fn_1} = location of pdb file name 1
                ->{pdb_fn_2} = location of pdb file name 2
   Returns:     parse_ali() alignment structure

=cut

sub get_resequiv_serial {

   my $params = shift ;

   my $modeller_bin = $params->{modeller_bin} ;
   my $pdb_fn ;
   $pdb_fn->[0] = $params->{pdb_fn_1} ;
   $pdb_fn->[1] = $params->{pdb_fn_2} ;

# Specify the temporary alignment TOP file, and the output alignment file.

   my ($ali_top_fn, $ali_top_fh, $ali_fn) ;
   my $temp_fh ;

   ($temp_fh, $ali_fn) =
      tempfile("resequiv.align.XXXXXX", SUFFIX => ".ali") ;
   close($temp_fh) ;

   ($ali_top_fh, $ali_top_fn) =
      tempfile( "resequiv.align.XXXXXX", SUFFIX => ".top" ) ;

# Generate the actual TOP file.

   print $ali_top_fh "# modeller_subsets_resequiv()\n#".time()."\n\n" ;
   print $ali_top_fh "READ_MODEL FILE = \'$pdb_fn->[0]\'\n" ;
   print $ali_top_fh "SEQUENCE_TO_ALI ALIGN_CODES = \'1_$pdb_fn->[0]\',".                            " ATOM_FILES = \'$pdb_fn->[0]\'\n" ;
   print $ali_top_fh "READ_MODEL FILE = \'$pdb_fn->[1]\'\n" ;
   print $ali_top_fh "SEQUENCE_TO_ALI ADD_SEQUENCE = on, ".                                          "ALIGN_CODES = ALIGN_CODES \'2_$pdb_fn->[1]\', ".                               "ATOM_FILES = ATOM_FILES \'$pdb_fn->[1]\'\n" ;

   print $ali_top_fh "\n" ;
   print $ali_top_fh "ALIGN GAP_PENALTIES_1D = -600 -400\n" ;
   print $ali_top_fh "ALIGN3D GAP_PENALTIES_3D = 0 2.0\n" ;

   print $ali_top_fh "\n" ;
   print $ali_top_fh "READ_MODEL FILE = \'$pdb_fn->[0]\'\n" ;
   print $ali_top_fh "READ_MODEL2 FILE = \'$pdb_fn->[1]\'\n" ;

   print $ali_top_fh "\n" ;
   print $ali_top_fh "SUPERPOSE FIT = on\n" ;
   print $ali_top_fh "WRITE_ALIGNMENT FILE=\'$ali_fn\'\n" ;

   close ($ali_top_fh) ;

# Specify the location of the MODELLER LOG file.

   my $ali_log_fn = $ali_top_fn ;
   $ali_log_fn =~ s/top$/log/ ;

# Run the TOP file through MODELLER.

   system("$modeller_bin $ali_top_fn >/dev/null 2>&1") ;

# Open the MODELLER LOG file.

   my $ali_parsed = parse_ali({ali_fn => $ali_fn}) ;

   unlink ($ali_top_fn, $ali_log_fn, $ali_fn) ;

   my $results = $ali_parsed ;

   return $ali_parsed ;

}



#=head2 OLD_parse_ali()
#
#   Title:       OLD_parse_ali()
#   Function:    OLD code to parse alignment file and return
#   Args:        $_ = alignment file name
#   Returns:     $_[0]->[0] = sequence of first protein
#                $_[0]->[1] = sequence of second protein
#                $_[1]->[0]->{resno1} = alipos1 - alignment position of
#                  this residue in pdb file 1
#                $_[1]->[1]->{resno2} = alipos2 - alignment position of
#                  this residue in pdb file 2
#
#=cut

sub OLD_parse_ali {

   my $alif = shift ;

   open (ALIF, $alif) ;
   my @headers ;
   my @seq ;
   my $cur_seq = -1 ;
   while (my $line = <ALIF>) {
      chomp $line;
      if (($line =~ /^structure/) || ($line =~ /^sequence/)) {
         push @headers, $line ;
	 $cur_seq++ ;
	 $seq[$cur_seq] = '' ;
      } elsif (($line !~ /^$/) && ($line !~ /^\>/) && ($line !~ /^\#/)) {
         $seq[$cur_seq] .= $line ;
      }
   }

   $seq[0] =~ s/\*$// ; $seq[1] =~ s/\*$// ;


   my $maxlen = length($seq[0]) ;
   my $t = length($seq[1]) ; if ($t > $maxlen) { $maxlen = $t; }

   my ($startres, $endres) ;
   foreach my $j ( 0 .. $#headers) {
      my @t = split(':', $headers[$j]) ;
      $t[2] =~ s/ //g ;
      $startres->[$j] = $t[2] ;
      $endres->[$j] = $t[2] ;
   }

   my $resequiv ;
   my $curres_1 = $startres->[0] ;
   my $curres_2 = $startres->[1];
   my $resno_2_alipos ;

   foreach my $j ( 0 .. ($maxlen - 1)) {

      if (substr($seq[0], $j, 1) ne '-') {
         $resno_2_alipos->[0]->{$curres_1} = $j ;
         $curres_1++ ;
      }
      if (substr($seq[1], $j, 1) ne '-') {
         $resno_2_alipos->[1]->{$curres_2} = $j ;
         $curres_2++ ;
      }

   }

   my $a; $a->[0] = $seq[0]; $a->[1] = $seq[1] ;
   return ($a, $resno_2_alipos) ;

}


1 ;
