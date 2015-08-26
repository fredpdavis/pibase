=head1 NAME

pibase::PDB::residues - Perl package to extract residue listing from a PDB file.

=head1 DESCRIPTION

Parses a pdb file and lists the residues.

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

package pibase::PDB::residues;

require Exporter;
@ISA = qw/Exporter/ ;
@EXPORT = qw/residue_info get_weirdo_resser/ ;

use strict;
use warnings;
use pibase::residue_math qw/residue_int/;


=head2 residue_info()

   Title:       residue_info()
   Function:    calculates a residue listing for a pdb file
   Args:        $_[0] - pdb file
                $_[1] - output file
                $_[2] - identifier - e.g. bdp_id [optional]
   Returns:     nothing

   Input file:  PDB file ($_[0]) 
      http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html
      http://msdlocal.ebi.ac.uk/docs/pdb_format/y_index.html

   Output file: Residue listing ($_[1])
      1 pdb identifier (e.g. bdp_id)
      2 Chain number
      3 Chain id
      4 residue number
      5 residue number (integer only)
      6 residue name
      7 Polymer type - p (protein) or n (nucleic acid)

=cut

sub residue_info {

   my $params = shift ;

   my $pibase_specs ;
   if (!exists $params->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $params->{pibase_specs} ;
   }

   if (!exists ($params->{pdb_fn})) {
      return (0, 'PDB file not specified]') ; }

   my $compress_fl = 0 ;
   if (exists $params->{gzip} && $params->{gzip} == 1) {
      $compress_fl = 1 ; }

   my $pdb_file = $params->{pdb_fn} ;

   my $outfile ;
   if (exists $params->{outfile}) {
      $outfile = $params->{outfile} ; }

   my $identifier ;
   if (exists $params->{identifier}) {
      $identifier = $params->{identifier} ;
   } else {
      $identifier = 'bdp_id' ;
   }

# Build hash with single letter abbreviations for standard residue names from PDB ATOM records.
   my %rescode = (
      'ALA' => 'A' ,
      'ARG' => 'R' ,
      'ASN' => 'N' ,
      'ASP' => 'D' ,
      'CYS' => 'C' ,
      'GLN' => 'Q' ,
      'GLU' => 'E' ,
      'GLY' => 'G' ,
      'HIS' => 'H' ,
      'HSD' => 'H' ,
      'ILE' => 'I' ,
      'LEU' => 'L' ,
      'LYS' => 'K' ,
      'MET' => 'M' ,
      'MSE' => 'M' ,
      'MEX' => 'M' ,
      'PHE' => 'F' ,
      'PRO' => 'P' ,
      'SER' => 'S' ,
      'THR' => 'T' ,
      'TRP' => 'W' ,
      'TYR' => 'Y' ,
      'VAL' => 'V',
      'UNK' => 'X',
      '  C' => 'c',
      '  G' => 'g',
      '  A' => 'a',
      '  T' => 't',
      '  U' => 'u',
      '  I' => 'i',
      ' +C' => 'c',
      ' +G' => 'g',
      ' +A' => 'a',
      ' +T' => 't',
      ' +U' => 'u',
      ' +I' => 'i'
   ) ;

   my @chain_id ;
   my @start_res ;
   my @end_res ;
   my @num_atm ;
   my @num_res ;
   my @num_het ;
   my @chain_seq ;
   my @chain_type ;

   my @chain_resno ; # array of residue list
   my @chain_resno_int ;
   my @chain_resna ;

   my $lastchain = '' ;
   my $lastresno = '' ;


   my $inTER = 0;

# Read a line from the pdb file.
   my $read_tcom = "cat" ;
   if ($pdb_file =~ /\.gz$/) {
      $read_tcom = $pibase_specs->{binaries}->{zcat} ;
   }

   my $tcom = "$read_tcom $pdb_file | " ;
   open (PDBFILE, $tcom) ;
   while (my $line = <PDBFILE>) {
      chomp $line;

      my $rectype ;

# If the line is an ATOM or HETATM record, set rectype to the record type.
      if ($line =~ /^ATOM/) {
         $rectype = 'ATOM';
      } elsif ($line =~ /^HETATM/) {
         $rectype = 'HETATM';
      } elsif ($line =~ /^TER/) {
         $inTER = 1;
         next;
      }

# If record type has been set,
      if (defined $rectype) {

# Extract:
# * chain id (22)
# * residue number (23 - 27) - note: includes insertion code
# * residue name (18 - 20)

         my $curchain = substr($line, 21, 1) ;

         if ($inTER) {
            if ($curchain eq $lastchain) {
               if ($rectype eq 'ATOM') { next; }
            } else {
               $inTER = 0 ;
            }
         }

         my $curresno = substr($line, 22, 5) ; $curresno =~ s/ //g ;
         my $curresna = substr($line, 17, 3) ;

         my $rescode ;
         if (exists ($rescode{$curresna})) {
   	    $rescode = $rescode{$curresna} ; }

# Initialize the new chain flag.
# If the current chain id is different from the last chain id, set the new chain flag.

         my $newchain = 0;
         if ($curchain ne $lastchain) {
   	    $newchain = 1; }

# Initialize the new residue flag.
# If the current residue number is different from the last residue number OR the new chain flag has been set, set the new residue flag.

         my $newres = 0;
         if (($curresno ne $lastresno) ||
	     $newchain) {
            $newres = 1; }

# If the new chain flag has been set,


         if ($newchain) {

# Create an entry in the chain arrays for the current chain.
# Store the current chain id.
   	    push @chain_id, $curchain ;

# Initialize the start residue, end residue, chain seqeuncee, number of residues, atoms, and HETATMS (defaults 0), and chaintype (default 'p')
   	    push @start_res, $curresno ;
   	    push @end_res, $curresno ;
   	    push @chain_seq, "" ;
   	    push @num_res, 0 ;
   	    push @num_atm, 0 ;
   	    push @num_het, 0;
   	    push @chain_type, 'p' ;

            push @chain_resno, [ ] ;
            push @chain_resno_int, [ ] ;
            push @chain_resna, [ ] ;

	 }

#If the new reisdue flag has been set and the current line is an ATOM record,
         if ($newres && ($rectype eq 'ATOM')) {

# Set the end residue number of the current chain to the current residue number.
   	    $end_res[$#end_res] = $curresno ;


# If the residue name is known (single letter abbreviation known),
   	    if (defined $rescode) {

# If the residue code contains any lower case letters (nucleotides), set the chain type to 'n'.
               if ($rescode =~ /[a-z]/) {
   	          $chain_type[$#chain_type] = 'n' ; }

# Concatenate the current residue letter abbreviation to the stored sequence for this chain.
	       push @{$chain_resno[$#chain_type]}, $curresno ;
	       my ($curresno_int, undef) = residue_int($curresno) ;
	       push @{$chain_resno_int[$#chain_type]}, $curresno_int ;
	       push @{$chain_resna[$#chain_type]}, $curresna ;

   	       $chain_seq[$#chain_id] .= $rescode{$curresna} ;

# Increment the number of residues in this chain.
   	       $num_res[$#num_res]++ ;

	    }

	 }

# If the record type is 'ATOM' and the single letter abbreviation is known, increment the number of atoms.
         if (($rectype eq 'ATOM') &&
	     (defined $rescode)) {
   	    $num_atm[$#num_atm]++;

# Elsif the record type is 'HETATM', increment the number of HETATMs.
         } elsif ($rectype eq 'HETATM') {
   	    $num_het[$#num_het]++; }

# Store the current residue number as the last residue number.
# Store the current chain id as the last chain id.

         $lastresno = $curresno ;
         $lastchain = $curchain ;

# Elsif the line is an ENDMDL entry, stop reading from the pdb file.

      } elsif ($line =~ /^ENDMDL/) {
         last ;
      } 
   }
   close(PDBFILE) ;

# Initialize the chain and residue count.
   my $ch_count = 1;

# Open the output file.
   my $outfile_fh ;
   if (defined $outfile) {
      $outfile =~ s/^\>// ;
      open($outfile_fh, "> $outfile") ; }

   my $results ;

# Iterate through the chains.
   foreach my $j (0 .. $#chain_id) {
# If the number of atoms is greater than 0,
      if (($num_atm[$j] > 0 ) && (defined $chain_resno_int[$j]->[0])) {

         my $res_count = 1 ;

# If the chain sequence is '', goto the next chain.
         if ($chain_seq[$j] eq '') {
            next; }

# Designate the fields to be displayed based on the current display format.
         foreach my $k ( 0 .. $#{$chain_resno[$j]} ) {
	    if (defined $outfile_fh) {
               my @outvals = ( $identifier, $ch_count, $chain_id[$j],$res_count,
                            $chain_resno[$j]->[$k], $chain_resno_int[$j]->[$k],
                            $chain_resna[$j]->[$k], $chain_type[$j] ) ;
               print $outfile_fh join("\t", @outvals)."\n" ;
            } else {
               push @{$results->{chain_no}}, $ch_count ;
               push @{$results->{chain_id}}, $chain_id[$j] ;
               push @{$results->{resno_serial}}, $res_count;
               push @{$results->{resno}}, $chain_resno[$j]->[$k];
               push @{$results->{resno_int}}, $chain_resno_int[$j]->[$k] ;
               push @{$results->{resna}}, $chain_resna[$j]->[$k];
               push @{$results->{chain_type}}, $chain_type[$j]  ;
            }
            $res_count++ ;
         }

# Increment the chain count.

         $ch_count++;

      }

   }
   close(OUTFILE) ;

   if (defined $compress_fl && defined $outfile) {
      system("gzip $outfile") ; }

   if (!defined $outfile) {
      return ($results, []); }

   if (!defined $outfile) {
      return (1, []); }

}


=head2 get_weirdo_resser()

   Title:       get_weirdo_resser()
   Function:    residue_info() altered to read in non-standard amino acids
                  (HETATOM OR ATOM) in as regular amino acids according to
                  MODELLER's restyp.lib mapping
   Args:        $_->{pdb_fn} - pdb file name
                $_->{outfile} - output file name [optional]
                $_->{identifier} - bdp_id identifier [optional]

   Returns: (if outfile is not specified)
      $_->{chain_no}->[i] - serial chain number
      $_->{chain_id}->[i] - chain identifier
      $_->{resno_serial}->[i] - serial residue number
      $_->{resno}->[i] - residue number
      $_->{resno_int}->[i] - integer portion of residue number
      $_->{resna}->[i] - residue name
      $_->{chain_type}->[i] - chain_type ('p' or 'n')
      $_->{weird_res_ser}->{resno_serial."\n".chain_no} = 1 if resna = MSE or MEX
      $_->{weird_res_raw}->{resno."\n".chain_id} = 1 if resna = MSE or MEX

   Output file: Residue listing ($_->{outfile}
      1. bdp identifier
      2. serial chain number
      3. chain identifier
      4. serial residue number
      5. residue number
      6. integer portion of residue number
      7. residue name
      8. chain type ('p' or 'n')

=cut

sub get_weirdo_resser {

   my $params = shift ;

   if (!exists ($params->{pdb_fn})) {
      return (0, 'PDB file not specified]') ; }

   my $pdb_file = $params->{pdb_fn} ;

   my $outfile ;
   if (exists $params->{outfile}) {
      $outfile = $params->{outfile} ; }

   my $identifier ;
   if (exists $params->{identifier}) {
      $identifier = $params->{identifier} ;
   } else {
      $identifier = 'bdp_id' ;
   }

# Build hash with single letter abbreviations for standard residue names from PDB ATOM records.

   my %rescode = (
      'ALA' => 'A' ,
      'ARG' => 'R' ,
      'ASN' => 'N' ,
      'ASP' => 'D' ,
      'CYS' => 'C' ,
      'GLN' => 'Q' ,
      'GLU' => 'E' ,
      'GLY' => 'G' ,
      'HIS' => 'H' ,
      'HSD' => 'H' ,
      'ILE' => 'I' ,
      'LEU' => 'L' ,
      'LYS' => 'K' ,
      'MET' => 'M' ,
      'PHE' => 'F' ,
      'PRO' => 'P' ,
      'SER' => 'S' ,
      'THR' => 'T' ,
      'TRP' => 'W' ,
      'TYR' => 'Y' ,
      'VAL' => 'V',
      'UNK' => 'X',
      '  C' => 'c',
      '  G' => 'g',
      '  A' => 'a',
      '  T' => 't',
      '  U' => 'u',
      '  I' => 'i',
      ' +C' => 'c',
      ' +G' => 'g',
      ' +A' => 'a',
      ' +T' => 't',
      ' +U' => 'u',
      ' +I' => 'i',
      'MSE' => 'M',
      'MEX' => 'C',
   ) ;

# got last 2 mappings from src/commands/read_model.F90
# special mappings from modlib/restyp.lib
#      'CSH' => 'C', 
#      'PR0' => 'P', 
#      'PRZ' => 'P',  #heterocyclic aromatic compound - not really a pro, but never used anyways - 
#      'ASX' => 'B', 
#      'GLX' => 'Z', 
##      'CSS' => 'C',  - not used see 1h32 models
#      'CYX' => 'C', 
#      'MSE' => 'M', 
#      'MEX' => 'C', 

   my @chain_id ;
   my @start_res ;
   my @end_res ;
   my @num_atm ;
   my @num_res ;
   my @num_het ;
   my @chain_seq ;
   my @chain_type ;

   my @chain_resno ; # array of residue list
   my @chain_resno_int ;
   my @chain_resna ;

   my $lastchain = '' ;
   my $lastresno = '' ;

# Read a line from the pdb file.

   open (PDBFILE, $pdb_file) ;
   while (my $line = <PDBFILE>) {

      chomp $line;

      my $rectype ;

# If the line is an ATOM or HETATM record, set rectype to the record type.
# - unless an MSE or MEX entry which is set as an ATOM record.

      if (($line =~ /^ATOM/) ||
          ( ($line =~ /^HETATM/) &&
	    length($line) > 20 &&
	    ((substr($line, 17, 3) eq 'MSE')||
	     (substr($line, 17, 3) eq 'MEX')) )) {
         $rectype = 'ATOM';
      } elsif ($line =~ /^HETATM/) {
         $rectype = 'HETATM'; }

# If record type has been set,

      if (defined $rectype) {

# Extract:
# * chain id (22)
# * residue number (23 - 27) - note: includes insertion code
# * residue name (18 - 20)

         my $curchain = substr($line, 21, 1) ;
         my $curresno = substr($line, 22, 5) ; $curresno =~ s/ //g ;
         my $curresna = substr($line, 17, 3) ;

         my $rescode ;
         if (exists ($rescode{$curresna})) {
   	    $rescode = $rescode{$curresna} ; }

# Initialize the new chain flag.
# If the current chain id is different from the last chain id, set the new chain flag.

         my $newchain = 0;
         if ($curchain ne $lastchain) {
   	    $newchain = 1; }

# Initialize the new residue flag.
# If the current residue number is different from the last residue number OR the new chain flag has been set, set the new residue flag.

         my $newres = 0;
         if (($curresno ne $lastresno) ||
	     $newchain) {
            $newres = 1; }

# If the new chain flag has been set,


         if ($newchain) {

# Create an entry in the chain arrays for the current chain.
# Store the current chain id.

   	    push @chain_id, $curchain ;

# Initialize the start residue, end residue, chain seqeuncee, number of residues, atoms, and HETATMS (defaults 0), and chaintype (default 'p')

   	    push @start_res, $curresno ;
   	    push @end_res, $curresno ;
   	    push @chain_seq, "" ;
   	    push @num_res, 0 ;
   	    push @num_atm, 0 ;
   	    push @num_het, 0;
   	    push @chain_type, 'p' ;

            push @chain_resno, [ ] ;
            push @chain_resno_int, [ ] ;
            push @chain_resna, [ ] ;

	 }


# If the new reisdue flag has been set and the current line is an ATOM record,
         if ($newres && ($rectype eq 'ATOM')) {

# Set the end residue number of the current chain to the current residue number.
   	    $end_res[$#end_res] = $curresno ;

# If the residue name is known (single letter abbreviation known),
   	    if (defined $rescode) {

# If the residue code contains any lower case letters (nucleotides), set the chain type to 'n'.
               if ($rescode =~ /[a-z]/) {
   	          $chain_type[$#chain_type] = 'n' ; }

# Concatenate the current residue letter abbreviation to the stored sequence for this chain.

	       push @{$chain_resno[$#chain_type]}, $curresno ;
	       my ($curresno_int, undef) = residue_int($curresno) ;
	       push @{$chain_resno_int[$#chain_type]}, $curresno_int ;
	       push @{$chain_resna[$#chain_type]}, $curresna ;

   	       $chain_seq[$#chain_id] .= $rescode{$curresna} ;

# Increment the number of residues in this chain.

   	       $num_res[$#num_res]++ ;

	    }

	 }

# If the record type is 'ATOM' and the single letter abbreviation is known, increment the number of atoms.


         if (($rectype eq 'ATOM') &&
	     (defined $rescode)) {
   	    $num_atm[$#num_atm]++; }

# Elsif the record type is 'HETATM', increment the number of HETATMs.

         elsif ($rectype eq 'HETATM') {
   	    $num_het[$#num_het]++; }

# Store the current residue number as the last residue number.
# Store the current chain id as the last chain id.

         $lastresno = $curresno ;
         $lastchain = $curchain ;

      }

# Elsif the line is an ENDMDL entry, stop reading from the pdb file.

      elsif ($line =~ /^ENDMDL/) {
         last ; }

   }
   close(PDBFILE) ;

# Initialize the chain and residue count.
   my $ch_count = 1;

# Open the output file.
   my $outfile_fh ;
   if (defined $outfile) {
      if ($outfile !~ /^\>/) {
         $outfile = '>'.$outfile ; }
      open($outfile_fh, $outfile) ; }

# Iterate through the chains.
   my $results ;


   my $weird_res_ser ;
   my $weird_res_raw ;
   foreach my $j (0 .. $#chain_id) {

# If the number of atoms is greater than 0,

      if (($num_atm[$j] > 0 ) && (defined $chain_resno_int[$j]->[0])) {

	 my $res_count = 1 ;

# If the chain sequence is '', goto the next chain.

         if ($chain_seq[$j] eq '') {
            next; }

# Designate the fields to be displayed based on the current display format.

         foreach my $k ( 0 .. $#{$chain_resno[$j]} ) {
	    if (defined $outfile_fh) {
	       my @outvals = ( $identifier, $ch_count, $chain_id[$j],$res_count,
	                    $chain_resno[$j]->[$k], $chain_resno_int[$j]->[$k],
			    $chain_resna[$j]->[$k], $chain_type[$j] ) ;
               print $outfile_fh join("\t", @outvals)."\n" ;
	    } else {
	       push @{$results->{chain_no}}, $ch_count ;
	       push @{$results->{chain_id}}, $chain_id[$j] ;
	       push @{$results->{resno_serial}}, $res_count;
	       push @{$results->{resno}}, $chain_resno[$j]->[$k];
	       push @{$results->{resno_int}}, $chain_resno_int[$j]->[$k];
	       push @{$results->{resna}}, $chain_resna[$j]->[$k];
	       push @{$results->{chain_type}}, $chain_type[$j]  ;

	       if ($chain_resna[$j]->[$k] eq 'MSE' ||
	           $chain_resna[$j]->[$k] eq 'MEX') {
	          $weird_res_ser->{$res_count."\n".$ch_count}++ ;
	          $weird_res_raw->{$chain_resno[$j]->[$k]."\n".$chain_id[$j]}++ ;
	       }
	    }
	    $res_count++ ;
	 }

# Increment the chain count.

         $ch_count++;

      }

   }
   close(OUTFILE) ;

   $results->{weird_res_ser} = $weird_res_ser ;
   $results->{weird_res_raw} = $weird_res_raw ;
   return ($results) ;

}


=head2 resno_2_serial()

   Function: calls residue_info() to return serial to actual residue number mapping
   Args:        $_->{pdb_fn} - pdb file name
   Returns:     $_->[0]->{serial_2_resno}->{serial_resno} = resno."\n".chain_id
                $_->[0]->{resno_2_serial}->{resno."\n".chain_id} =
                  serial residue number}

=cut

sub resno_2_serial {

   my $params = shift ;
   if (! exists $params->{pdb_fn}) {
      return ({}, ['PDB file not specified]']); }

   my $pdbfile = $params->{pdb_fn} ;
   if (! -e $pdbfile) {
      return ({}, ['PDB file does not exist']); }

   my ($res_info, $errors) = residue_info({pdb_fn => $pdbfile}) ;

   my $serial_2_resno ;
   my $resno_2_serial ;
   foreach my $j ( 0 .. $#{$res_info->{resno_serial}}) {
      my $sig = $res_info->{resno}->[$j]."\n".$res_info->{chain_id}->[$j] ;
      $serial_2_resno->{$res_info->{resno_serial}->[$j]} = $sig ;
      $resno_2_serial->{$sig} = $res_info->{resno_serial}->[$j] ;
   }

   my $results ;
   $results->{serial_2_resno} = $serial_2_resno ;
   $results->{resno_2_serial} = $resno_2_serial ;
  
   return ($results, []) ;

}

1 ;
