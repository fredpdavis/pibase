=head1 NAME

pibase::PDB::chains - Perl module to extract chain listing from a pdb file.

=head1 DESCRIPTION

Perl module to parse a pdb file and lists the chains.

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

package pibase::PDB::chains;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT = qw/chain_info/ ;

=head2 SUB chain_info()

   Title: chain_info()
   Function: calculates chain listing for a pdb file

   Args:        $_[0] = pdb_file
                $_[1] = outfile
                $_[2] = bdp_id identifier [optional]
   Returns:     nothing

   Input file:  PDB file ($_[0]) 
      http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html
      http://msdlocal.ebi.ac.uk/docs/pdb_format/y_index.html

   Output file: Chain listing ($_[1])
      1 pdb identifier (e.g. bdp_id)
      2 Chain number
      3 Chain id
      4 Chain type - p (protein) or n (nucleic acid)
      5 null
      6 null
      7 start residue number
      8 start residue number (integer only)
      9 end residue number
     10 end residue number (integer only)
     11 number of residues
     12 number of atoms
     13 number of het atoms
     14 chain sequence

=cut

sub chain_info {

   my $params = shift ;

   my $pibase_specs ;
   if (!exists $params->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $params->{pibase_specs} ;
   }

   if (!exists ($params->{pdb_fn})) {
      return {error_fl => 'PDB file not specified'} ; }

   my $pdb_file = $params->{pdb_fn}; #old 0

   my $outfile ;
   if (exists $params->{outfile}) {
      $outfile = $params->{outfile} ;
   } else {
      $outfile = '-' ;
   }

   my $identifier ;
   if (exists $params->{identifier}) {
      $identifier = $params->{identifier} ;
   } else {
      $identifier = 'bdp_id' ;
   }

# Build hash with single letter abbreviations for standard residue names form PDB ATOM records.
   my %rescode = (
      ALA => 'A' ,
      ARG => 'R' ,
      ASN => 'N' ,
      ASP => 'D' ,
      CYS => 'C' ,
      GLN => 'Q' ,
      GLU => 'E' ,
      GLY => 'G' ,
      HIS => 'H' ,
      HSD => 'H' ,
      ILE => 'I' ,
      LEU => 'L' ,
      LYS => 'K' ,
      MET => 'M' ,
      PHE => 'F' ,
      PRO => 'P' ,
      SER => 'S' ,
      THR => 'T' ,
      TRP => 'W' ,
      TYR => 'Y' ,
      VAL => 'V',
      UNK => 'X',
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


   my (@start_res, @end_res) ;
   my (@num_atm, @num_res, @num_het) ;
   my (@chain_id, @chain_seq, @chain_type) ;

   my $last_chain = '' ;
   my $last_resno = '' ;

   my $start_res_fl = 0 ;

   my $inTER = 0;

# Read a PDB line of input.
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

         my $cur_chain = substr($line, 21, 1) ;
         if (!defined $cur_chain) {
            print STDERR "ERROR (line $line): no chain_id\n" ;}

         if ($inTER) {
            if ($cur_chain eq $last_chain) {
               if ($rectype eq 'ATOM') { next; }
            } else {
               $inTER = 0 ;
            }
         }

         my $cur_resno = substr($line, 22, 5) ; #includes insertion code
         $cur_resno =~ s/ //g ;
         my $cur_resna = substr($line, 17, 3) ;

         my $rescode ;
         if (exists ($rescode{$cur_resna})) {
   	    $rescode = $rescode{$cur_resna} ; }

# Initialize the new chain flag.
# If the current chain id is different from the last chain id, set the new chain flag.
         my $newchain = 0;
         if ($cur_chain ne $last_chain) {
   	    $newchain = 1; }

# Initialize the new residue flag.
# If the current residue number is different from the last residue number OR the new chain flag has been set, set the new residue flag.
         my $newres = 0;
         if ( ($cur_resno ne $last_resno) ||
	       $newchain) {
            $newres = 1; }

# If the new chain flag has been set,
         if ($newchain) {

# Create an entry in the chain arrays for the current chain.
# Store the current chain id.
   	    push @chain_id, $cur_chain ;

# Initialize the start residue, end residue, chain seqeuncee, number of residues, atoms, and HETATMS (defaults 0), and chaintype (default 'p')
   	    push @start_res, $cur_resno ;
   	    push @end_res, $cur_resno ;

   	    push @chain_seq, "" ;
   	    push @num_res, 0 ;
   	    push @num_atm, 0 ;
   	    push @num_het, 0;
   	    push @chain_type, 'p' ;

# Initialize the start residue flag to 0
	    $start_res_fl = 0 ;

	 }

# If the start_residue flag has not been set AND the current entry is an ATOM entry AND is a residue of known type, set the start residue flag
	 if ( !$start_res_fl && $rectype eq 'ATOM' && defined $rescode) {

# Reset the start and end residue number for this chain to the current residue number
	    $start_res[$#start_res] = $cur_resno ;
	    $end_res[$#end_res] = $cur_resno ;

# Set the start residue flag
	    $start_res_fl = 1 ;

	 }

# If the new reisdue flag has been set and the current line is an ATOM record,
         if ($newres && ($rectype eq 'ATOM')) {

# Set the end residue number of the current chain to the current residue number.
   	    $end_res[$#end_res] = $cur_resno ;


# If the residue name is known (single letter abbreviation known),
   	    if (defined $rescode) {

# If the residue code contains any lower case letters (nucleotides), set the chain type to 'n'.
               if ($rescode =~ /[a-z]/) {
   	          $chain_type[$#chain_type] = 'n' ; }

# Concatenate the current residue letter abbreviation to the stored sequence for this chain.
   	       $chain_seq[$#chain_id] .= $rescode{$cur_resna} ;

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
   	    $num_het[$#num_het]++;
         }

# Store the current residue number as the last residue number.
# Store the current chain id as the last chain id.

         $last_resno = $cur_resno ;
         $last_chain = $cur_chain ;

      }

# Elsif the line is an ENDMDL entry, stop reading from the pdb file.
      elsif ($line =~ /^ENDMDL/) {
         last ; }

   }
   close(PDBFILE) ;

# Initialize the chain count.
   my $ch_count = 1;

# Open the output file.
   my $fh ;
   open($fh->{out}, ">$outfile") ;
# Iterate through the chains.
   foreach my $j (0 .. $#chain_id) {

# If the number of atoms is greater than 0,
      if ($num_atm[$j] > 0 ) {

# Remove spaces from the start and end residue numbers.
         $start_res[$j] =~ s/ //g ;
         $end_res[$j] =~ s/ //g ;

# Extract the integer portion of the start and end residue numbers.
         my ($start_res_int) = ($start_res[$j] =~ /(-?[0-9]+)/) ;
         my ($end_res_int) = ($end_res[$j] =~ /(-?[0-9]+)/) ;
	 if ((! defined $start_res_int) || (! defined $end_res_int)) {
	    next ; }

# Remove spaces from the chain id.
         $chain_id[$j] =~ s/ //g ;

# If the chain id is '', set it to ' '.
         if ($chain_id[$j] eq '') {
            $chain_id[$j] = ' ';}

# If the chain sequence is '', goto the next chain.
         if ($chain_seq[$j] eq '') {
            next; }

# Designate the fields to be displayed based on the current display format.
         my @outvals = ( $identifier, $ch_count, $chain_id[$j], $chain_type[$j],
	                 '\N', '\N',
			 $start_res[$j], $start_res_int, $end_res[$j],
			 $end_res_int, $num_res[$j], $num_atm[$j], $num_het[$j],
			 $chain_seq[$j] ) ;

         print {$fh->{out}} join("\t", @outvals)."\n" ;

# Increment the chain count.
         $ch_count++;

      }

   }
   close($fh->{out}) ;

   return {success => 1} ;
}

1 ;
