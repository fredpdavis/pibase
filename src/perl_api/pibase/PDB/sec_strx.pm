=head1 NAME

pibase::PDB::sec_strx - Obtain secondary structure assignments from dsspcmbi.

=head1 DESCRIPTION

Interface to DSSP to calculate secondary structure for PIBASE structures.

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

package pibase::PDB::sec_strx ;

require Exporter;
@ISA = qw/Exporter/ ;
@EXPORT = qw/run_dssp parse_dssp parse_dssp_basic/ ;

use strict;
use warnings;

use File::Temp qw/tempfile/ ;
use Sys::Hostname ;
use pibase qw/locate_binaries/ ;


=head2 SUB get_sec_strx()

   Function: calls DSSP to get a PDB file's residue secondary structure assignment
   Args: $_[0] - pdb file
         $_[1] - output file
         $_[2] - bdp identifier - e.g. bdp_id [optional]

   Returns:     nothing
   Output file:
      1 pdb identifier (e.g. bdp_id)
      2 Chain number
      3 Chain id
      4 residue number
      5 residue number (integer only)
      6 residue name
      7 Polymer type - p (protein) or n (nucleic acid)

=cut

sub run_dssp {

   my $params = shift ;

   my $binaries = pibase::locate_binaries()  ;
   my $altloc_filter = $binaries->{'altloc_filter'} ;
   my $altloc_check = $binaries->{'altloc_check'} ;
   my $dssp_bin = $binaries->{'dssp'} ;

   if (!exists ($params->{pdb_fn})) {
      return (0, 'PDB file not specified]') ; }

   my $pdb_file = $params->{pdb_fn} ;

#make local copy of the pdb file; zcat if necessary

   my ($temp_fh, $localbdp) = tempfile(SUFFIX => ".pdb") ; close($temp_fh) ;
   if ($pdb_file =~ /\.gz$/) {
      my $tcom = $binaries->{zcat}." $pdb_file > ".$localbdp ;
      system($tcom) ;
   } else {
      pibase::safe_copy($pdb_file, $localbdp) ;
   }

   if (!-e $localbdp) {
      return {error_fl => "ERROR: pdb file access error - couldnt copy locally"}; }

   my $outfile ;
   if (exists $params->{outfile}) {
      $outfile = $params->{outfile} ; }


# Create and open dssp output files.

   my $host = hostname() ;

   my ($fn, $fh) ;
   if ( (! defined $outfile) || (-e $outfile) ) {
      ($fh->{t}, $fn->{dssp_out}) =
         tempfile("dssp.$host.XXXX", SUFFIX => ".out");
      close ($fh->{t}) ;
   } else {
      $fn->{dssp_out} = $outfile ;
   }

   ($fh->{t}, $fn->{dssp_err}) =
      tempfile("dssp.$host.XXXX", SUFFIX => ".err");
   close($fh->{t}) ;

# Check if pdb file has altloc flags set


   my $altloc_fl = `$altloc_check < $localbdp` ; chomp $altloc_fl ;
   my $tcom ;

   if ($altloc_fl) {
      $tcom = "$altloc_filter $localbdp" ;
   } else {
      $tcom = "cat $localbdp" ;
   }

   $tcom .= " | $dssp_bin -- 2>$fn->{dssp_err} >$fn->{dssp_out}" ;
   system($tcom) ;

   {
      my ($newdev, $newino) = stat($localbdp) ;
      my ($olddev, $oldino) = stat($pdb_file) ;
      if ($newdev != $olddev || $newino != $oldino) {
         unlink $localbdp ; }
   }


   if (!-s $fn->{dssp_out}) {
      unlink $fn->{dssp_err}, $fn->{dssp_out} ;
      return (0, 'dssp execution error') ;
   }

   return (1, { out=>$fn->{dssp_out}, err=>$fn->{dssp_err} }) ;

}


=head2 SUB parse_dssp()

   Title:    parse_dssp()
   Function: parses DSSP output and returns a hash of assignment data
   Args: $_[0] - DSSP output file
         $_[1] - output file
         $_[2] - bdp identifier - e.g. bdp_id [optional]

   Returns:     $->{detail}->{resno."\n".chain_id} = H|G|I|B|E|T|S|' '
                $->{basic}->{resno."\n".chain_id} = H|B|T|' '
                $->{ordering} = [resno1."\n".chain_id1, resno2."\n".chain_id2...]
                $->{ssnum}->{resno."\n".chain_id} = secondary structure element #
                   counts contiguous stretches of same particular detailed sec
                   strx assignment)
                $->{ssnum_basic}->{resno."\n".chain_id} = secondary structure
                   element number - counts contiguous stretches of basic secondary
                   structurea ssignment

=cut

sub parse_dssp {

   my $dssp_fn = shift ;
#cols from Bio::Structure::SecStrs::DSSP::Res_raw
#http://www.ensembl.org/Docs/Pdoc/bioperl-live/Bio/Structure/SecStr/DSSP/Res_raw.html
#   my $dssp_cols = {
#      'pdb_resnum' => [ 5, 5 ],
#      'insertionco' => [ 10, 1 ],
#      'pdb_chain' => [ 11, 1 ],
#      'amino_acid' => [ 13, 1 ],
#      'term_sig' => [ 14, 1 ],
#      'ss_summary' => [ 16, 1 ],
#      '3tph' => [ 18, 1 ],
#      '4tph' => [ 19, 1 ],
#      '5tph' => [ 20, 1 ],
#      'geo_bend' => [ 21, 1 ],
#      'chirality' => [ 22, 1 ],
#      'beta_br1la' => [ 23, 1 ],
#      'beta_br2la' => [ 24, 1 ],
#      'bb_part1nu' => [ 25, 4 ],
#      'bb_part2nu' => [ 29, 4 ],
#      'betash_lab' => [ 33, 1 ],
#      'solv_acces' => [ 34, 4 ],
#      'hb1_nh_o_p' => [ 39, 6 ],
#      'hb1_nh_o_e' => [ 46, 4 ],
#      'hb1_o_hn_p' => [ 50, 6 ],
#      'hb1_o_hn_e' => [ 57, 4 ],
#      'hb2_nh_o_p' => [ 61, 6 ],
#      'hb2_nh_o_e' => [ 68, 4 ],
#      'hb2_o_hn_p' => [ 72, 6 ],
#      'hb2_o_hn_e' => [ 79, 4 ],
#      'tco' => [ 85, 6 ],
#      'kappa' => [ 91, 6 ],
#      'alpha' => [ 97, 6 ],
#      'phi' => [ 103, 6 ],
#      'psi' => [ 109, 6 ],
#      'x_ca' => [ 115, 7 ],
#      'y_ca' => [ 122, 7 ],
#      'z_ca' => [ 129, 7 ]
#   } ;

   my $basic_secstrx = {
      'H' => 'H',
      'G' => 'H',
      'I' => 'H',
      'B' => 'B',
      'E' => 'B',
      'T' => 'T',
      'S' => 'T',
      ' ' => ' ',
   } ;


   open (DSSP_FH, $dssp_fn) ;
   my $pastheaders_fl = 0 ;
   my $sec_strx ;
   my $lastassign = 0; my $ssnum = 0 ;
   my $lastassign_basic = 0; my $ssnum_basic = 0 ;
   while (my $line =  <DSSP_FH>) {
      chomp $line;
      if ($pastheaders_fl) {
# ressig is resno."\n".chainid
	 my $resno = substr($line, 5, 6) ; $resno =~ s/ //g ;

	 if ((!defined $resno) || ($resno eq '')) {
	    next ; }

	 my $chain = substr($line, 11, 1) ;
         my $t_ressig = $resno."\n".$chain ;

	 my $t_assign = substr($line, 16 , 1) ;
         $sec_strx->{detail}->{$t_ressig} = $t_assign ;
	 $sec_strx->{basic}->{$t_ressig} = $basic_secstrx->{$t_assign} ;
	 push @{$sec_strx->{ordering}}, $t_ressig ;

	 if ($lastassign ne $t_assign) {
	    $ssnum++ ; }
	 $sec_strx->{ssnum}->{$t_ressig} = $ssnum ;

	 if ($lastassign_basic ne $basic_secstrx->{$t_assign}) {
	    $ssnum_basic++ ; }
	 $sec_strx->{ssnum_basic}->{$t_ressig} = $ssnum_basic ;

         $lastassign = $t_assign;
         $lastassign_basic = $basic_secstrx->{$t_assign};
      } elsif ($line =~ /^  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O/) {
         $pastheaders_fl = 1 ;
      }
   }
   close(DSSP_FH) ;

   return $sec_strx ;
}


sub OLD_parse_dssp_basic { #rolled into parse_dssp

   my $dssp_fn = shift ;
#cols from Bio::Structure::SecStrs::DSSP::Res_raw
#http://www.ensembl.org/Docs/Pdoc/bioperl-live/Bio/Structure/SecStr/DSSP/Res_raw.html

   my $basic_secstrx = {
      'H' => 'H',
      'G' => 'H',
      'I' => 'H',
      'B' => 'B',
      'E' => 'B',
      'T' => 'T',
      'S' => 'T',
      ' ' => ' ',
   } ;


   open (DSSP_FH, $dssp_fn) ;
   my $pastheaders_fl = 0 ;
   my $sec_strx ;
   my $lastassign = 0; my $ssnum = 0 ;
   while (my $line =  <DSSP_FH>) {
      chomp $line;
      if ($pastheaders_fl) {
# ressig is resno."\n".chainid
	 my $resno = substr($line, 5, 6) ; $resno =~ s/ //g ;

	 if ((!defined $resno) || ($resno eq '')) {
	    next ; }

	 my $chain = substr($line, 11, 1) ;
         my $t_ressig = $resno."\n".$chain ;

	 my $t_assign = substr($line, 16 , 1) ;
         $sec_strx->{detail}->{$t_ressig} = $t_assign ;
	 $sec_strx->{basic}->{$t_ressig} = $basic_secstrx->{$t_assign} ;
	 push @{$sec_strx->{ordering}}, $t_ressig ;

	 if ($lastassign ne $basic_secstrx->{$t_assign}) {
	    $ssnum++ ; }
	 $sec_strx->{ssnum}->{$t_ressig} = $ssnum ;

         $lastassign = $basic_secstrx->{$t_assign};
      } elsif ($line =~ /^  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O/) {
         $pastheaders_fl = 1 ;
      }
   }
   close(DSSP_FH) ;

   return $sec_strx ;
}

1 ;
