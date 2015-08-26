=head1 NAME

pibase::CATH - interface for CATH domain database data processing 

=head1 DESCRIPTION

This module contains routines to process and reformat CATH release files
for PIBASE import.  Files are from: http://www.cathdb.info

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

package pibase::CATH;
use Exporter;
use strict;
use warnings;

our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/cath_clean_cddf cath_parse_cdf_domainlist cath_clean_clf cath_clean_cnf_contraction cath_clean_cnf cath_pibase_import_domains/ ;

use pibase ;
use File::Temp qw/tempfile/ ;
use pibase::residue_math qw/residue_int/;


=head2 cath_clean_cddf()

   Title:       cath_clean_cddf()
   Function:    cleans the CATH CDDF file
   STDIN:       CATH CDDF file
   STDOUT:      cleaned up CATH CDDF file
   Args:        nothing
   Returns:     nothing

=cut

sub cath_clean_cddf {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }

   my @headers = qw/DOMAIN pdb_id pdb_chain_id domain_no class_id arch_id topol_id homol_id CLASS ARCH TOPOL HOMOL DLENGTH DSEQH DSEQS segment_no start_resno end_resno SLENGTH SSEQH SSEQS/ ;
   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }

   my @segment_entries = qw/SEGMENT SRANGE start_resno end_resno SLENGTH SSEQH SSEQS/ ;

   my %record_val ;

   while (my $line = readline($fh->{in})) {
      if ($line =~ /^\#/) {next;}
      chomp $line;

#If we are at the end of a domain entry ('//'), purge all variable contents.

      if ($line eq '//') {

         %record_val = ();
         undef %record_val ;

      } elsif ($line eq 'ENDSEG') {

         foreach my $key (@headers) {
         
            if ( !(exists $record_val{$key}) ||
                 !(defined $record_val{$key}) ||
                 ($record_val{$key} eq 'void') ||
                 ($record_val{$key} eq '') ) {

#If an undefined field is start_resno or end_resno,set to ' '. Otherwise set to the NULL string

               if (($key eq 'start_resno') || ($key eq 'end_resno') ||
		   ($key eq 'DSEQH') || ($key eq 'DSEQS') ||
		   ($key eq 'SSEQH') || ($key eq 'SSEQS')) {
		   $record_val{$key} = ' ';
               } elsif ( ($key eq 'DLENGTH') || ($key eq 'SLENGTH')) {
		     $record_val{$key} = '0' ;
               } else {
   	             $record_val{$key} = '\N' ;
               }
            } else {
               $record_val{$key} =~ s/^\s*// ;
               $record_val{$key} =~ s/\s*$// ;
   	    }
         }

   	 my @outfields = @record_val{@headers} ;
   	 print {$fh->{out}} join("\t", @outfields)."\n" ;

         delete @record_val{@segment_entries} ;

      } else {

         my $record_type = substr($line,0,10) || $line ;
         $record_type =~ s/\s*$//g ;
         my $record_val = substr($line,10) || '' ;

         if (exists $record_val{$record_type}) {
            $record_val{$record_type} .= $record_val ;
         } else {
            $record_val{$record_type} = $record_val ;
         }

   	 if ($record_type eq 'DOMAIN') {
            $record_val{pdb_id} = substr($record_val, 0, 4) ;
            $record_val{pdb_chain_id} = substr($record_val, 4, 1) ;
            $record_val{domain_no} = substr($record_val, 5, 2) ;

   	 } elsif ($record_type eq 'SEGMENT') {
            $record_val{segment_no} = substr($record_val, 8, 1) ;

         } elsif ($record_type eq 'SRANGE') {
   	    ($record_val{start_resno}) =
               ($record_val =~ /START\=(-?[0-9]+\S?)/);
   	    ($record_val{end_resno}) =
               ($record_val =~ /STOP\=(-?[0-9]+\S?)/);

   	 } elsif ($record_type eq 'CATHCODE') {

   	    ($record_val{class_id}, $record_val{arch_id},
             $record_val{topol_id}, $record_val{homol_id}) =
               split(/\./, $record_val) ;

         }
      }
   }
   close($fh->{in}) ;
   close($fh->{out}) ;

}


=head2 cath_parse_cdf_domainlist()

   Title:       cath_parse_cdf_domainlist()
   Function:    cleans the CATH CDDF file
   STDIN:       CATH CDDF file
   STDOUT:      cleaned up CATH CDDF file
   Args:        nothing
   Returns:     nothing

=cut


sub cath_parse_cdf_domainlist {

   my $in = shift ;
   my $fn = $in->{fn} ;

   my $fh;
   $fh->{in} = open($in->{fn}) ;
   if (!exists $in->{out_fn}) {
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{out}, ">$in->{out_fn}") ;
   }

   if (!exists $fn->{domall} || !exists $fn->{domainlist}) {
      return {error_fl => "ERROR: must specify domall and domainlist files"};}

   my @headers = qw/DOMAIN pdb_id pdb_chain_id domain_no segment_no start_resno end_resno/ ;
   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }

   my %record_val ;

   open (DOMALL, $fn->{domall}) ;
   while (my $line = <DOMALL>) {

      if ($line !~ /^\#/) {

         chomp $line;
         my @fields = split(' ', $line) ;

         my $domain_base = shift @fields ;
         my $pdb_id = substr($domain_base, 0, 4) ;
         my $chain_id = substr($domain_base, 4, 1) ;

         my $numdomains = shift @fields ; $numdomains =~ s/^D// ;
         my $numfrags = shift @fields ; $numfrags =~ s/^F// ;

         if ($numdomains == 0) {next;} #no domain entries

         my $defs;
         my $domnum = 1;

         while ($domnum <= $numdomains) {

            my $numsegs = shift @fields;
            my $segnum = 1 ;

            while ($segnum <= $numsegs) {

               my $start_ch = shift @fields ;
               my $start_res = shift @fields ;
               my $start_ins = shift @fields ;

               my $end_ch = shift @fields ;
               my $end_res = shift @fields ;
               my $end_ins = shift @fields ;

               if ($start_ins ne '-') {$start_res .= $start_ins}
               if ($end_ins ne '-') {$end_res .= $end_ins}

               my $domain_name = $pdb_id.$chain_id.$domnum ;

               my @outvals = ($domain_name, $pdb_id, $chain_id, $domnum,
                              $segnum, $start_res, $end_res) ;

               print join("\t", @outvals)."\n" ;
               $segnum++ ;
            }

            $domnum++ ;
         }

      }

   }
   close(DOMALL) ;

   open (DOMLIST, $fn->{domainlist}) ;
   while (my $line = <DOMLIST>) {
      if ($line =~ /^\#/) {next;}
      chomp $line;
      my @fields = split(' ', $line) ;

      my $domain_name = $fields[0] ;

#if domain_name ends in '0', extract (and display) pdb_id, chain_id, and domain number

      if ($domain_name !~ /0$/) {next;}


      my $pdb_id = substr($domain_name, 0, 4) ;
      my $chain_id = substr($domain_name, 4, 1) ;
      my $domnum = substr($domain_name, 5, 2) ;

      my @outvals = ($domain_name, $pdb_id, $chain_id, $domnum, 1, '', '') ;
      print join("\t", @outvals)."\n" ;
   }
   close(DOMLIST) ;

}


=head2 cath_clean_clf()

   Title:       cath_clean_clf()
   Function:    cleans the CATH CLF file
   STDIN:       CATH CLF file
   STDOUT:      cleaned up CATH CLF file
   Args:        nothing
   Returns:     nothing

=cut

sub cath_clean_clf {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }

   my @headers = qw/domain_name pdbcode chain domain_no class arch topol homol s35no s95no s100no domain_length/ ;
   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }
 
   while (my $line = readline($fh->{in})) {
      if ($line =~ /^\#/) {next;}

      my ( $domain_name, $classno, $archno, $topno, $homolno,
	   $s35no, $s95no, $s100no, $domain_length,
	   $resolution ) = split(' ', $line) ;

      my $pdb_id = substr($domain_name,0,4);
      my $chain_id = substr($domain_name,4,1) ;
#fpd080229_1825 - CLF format 2.0 has 2 characters for
#  the domain number - changed from original substr($domain_name,5,1) ;
# chaned on two more lines as well: 128,247
      my $domain_no = substr($domain_name,5,2) ;
      my $obsolete_fl = 0 ;

#If the resolution is 1000, set the obsolete flag.

      if ($resolution == 1000) {
         $obsolete_fl = 1 ; }

#If the resolution is equal to or greater than 999, set resolution to the NULL string.

      if ($resolution >= 999) {
         $resolution = '\N' ; }

      my @outfields = ( $domain_name, $pdb_id, $chain_id, $domain_no,
	                $classno, $archno, $topno, $homolno, $s35no,
			$s95no, $s100no, $domain_length ) ;

      foreach my $j (0 .. $#outfields) {

#If the field is (undefined) or ('-' && not the chain_id), set it to NULL.

   	 if ( (! defined $outfields[$j]) ||
                 (($outfields[$j] eq '-') && ($j != 2))) {
               $outfields[$j] = '\N' ;

#Otherwise, remove leading and trailing spaces.

         } else {
               $outfields[$j] =~ s/^\s*// ;
               $outfields[$j] =~ s/\s*$// ;
         }
      }

      print {$fh->{out}} join("\t", @outfields)."\n" ;
   }
   close($fh->{in}) ;
   close($fh->{out}) ;

}



=head2 cath_clean_cnf_contraction()

   Title:       cath_clean_cnf_contraction()
   Function:    fixes the concatenation problem with the CATH CNF file
   STDIN:       CATH CNF file
   STDOUT:      concat-fixed CATH CNF file
   Args:        nothing
   Returns:     nothing

=cut

sub cath_clean_cnf_contraction {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }

   while (my $line = readline($fh->{in})) {
      chomp $line;

      if ($line !~ /^#/) {

#While the line contains a CATH code in the middle of the line,

         while ($line =~ / .+[0-9]\.[0-9]+\.[0-9].+ /) {
#Split the line before the erroneous CATH code.

            $line =~ s/( .+)([0-9]\.[0-9]+\.[0-9].+ )/$1\n$2/g ; }

      }

      print {$fh->{out}} $line."\n" ;
   }
   close($fh->{in}) ;
   close($fh->{out}) ;

}


=head2 cath_clean_cnf_contraction()

   Title:       cath_clean_cnf()
   Function:    fixes the concatenation problem with the CATH CNF file
   STDIN:       CATH CNF file (preferably concat-fixed)
   STDOUT:      cleaned CATH CNF file
   Args:        nothing
   Returns:     nothing

=cut

sub cath_clean_cnf {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }

   my @headers = qw/class arch topol homol representative_dom description/ ;
   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }

   while (my $line = readline($fh->{in})) {

      if ($line =~ /^\#/) {next;}

      chomp $line;

      my ($raw_domain_id, $description) = ($line =~ /(.*?)\:(.*)/) ;
      my ($CATH_code, $represent_dom) = split(' ', $raw_domain_id) ;
      my @cathfields = split(/\./, $CATH_code) ;


      while ($#cathfields < 3) { push @cathfields, '0' ; }
      my @outfields = @cathfields;


      push @outfields, ($represent_dom, $description) ;

      foreach my $j (0 .. $#outfields) {

         if ((! (defined $outfields[$j])) || $outfields[$j] eq '-') {
            $outfields[$j] = '\N' ;
         } else {
            $outfields[$j] =~ s/^\s*// ;
            $outfields[$j] =~ s/\s*$// ;
         }

      }

      print {$fh->{out}} join("\t", @outfields)."\n" ;
   }
   close($fh->{in}) ;
   close($fh->{out}) ;

}


=head2 pibase_import_cath_domains()

   Title:       pibase_import_cath_domains()
   Function:    reformats raw CATH tables imported into pibase as generic
                  subsets tables
   In tables:   1. cath_domain_list
                2. cath_domall_boundaries
                3. cath_names

   Out tables:  1. subsets
                2. subsets_class
                3. subsets_details

   Args:        nothing
   Returns:     nothing

=cut

sub pibase_import_cath_domains {

   require DBI ;

   my $in = shift ;
   my $overwrite_fl = 0;
   if (exists $in->{overwrite_fl}) {$overwrite_fl = $in->{overwrite_fl};}

#Connect to the pibase database.
   my ($dbh, $pibase) = pibase::connect_pibase() ;
   my $pibase_specs = $in->{pibase_specs} ;

# remove CATH entry from subsets_source table
   pibase::mysql_runcom($dbh,
      'DELETE FROM subsets_source WHERE subset_source = "cath"') ;

   pibase::mysql_runcom($dbh,
      'INSERT INTO subsets_source'.
      '(subset_source_id, subset_source, version) '.
      'VALUES("", "cath","'.
      $in->{pibase_specs}->{external_data}->{cath}->{ver}.'") ;');

   my ($subset_source_id) = pibase::mysql_singleval($dbh,
      'SELECT subset_source_id FROM subsets_source '.
      'WHERE subset_source = "cath"') ;

#Set subsets, subsets_class, and subsets_details output files.

   if (!-s $pibase_specs->{external_data}->{cath}->{processed_files}->{"subsets.cath"} || $overwrite_fl == 1) {
      my $fh ;
      open($fh->{subsets}, ">".
   $pibase_specs->{external_data}->{cath}->{processed_files}->{"subsets.cath"});

      open($fh->{subsets_class}, ">".
   $pibase_specs->{external_data}->{cath}->{processed_files}->{"subsets_class.cath"});

      open($fh->{subsets_details}, ">".
   $pibase_specs->{external_data}->{cath}->{processed_files}->{"subsets_details.cath"});
   

#mysql(SELECT domain_name, pdb_id, class, arch, topol, homol FROM cath_domain_list)

      my $cath_dl ;
      ( $cath_dl->{domain_name},
        $cath_dl->{pdb_id},
        $cath_dl->{class},
        $cath_dl->{arch},
        $cath_dl->{topol},
        $cath_dl->{homol} ) = 
         pibase::mysql_fetchcols($dbh,"SELECT domain_name, pdb_id, class, arch, topol, homol FROM cath_domain_list") ;
   
#Iterate through the cath_domain_list SELECT results.
   
      foreach my $j (0 .. $#{$cath_dl->{domain_name}}) {
   
#Create a class string by concatenating the CATH class, arch, topol, and homol ids with a '.' delimiter.
   
         my $class = join('.', ( $cath_dl->{class}->[$j],
                                 $cath_dl->{arch}->[$j],
                                 $cath_dl->{topol}->[$j],
                                 $cath_dl->{homol}->[$j] ) ) ;
   
#Designate and display output fields. The subset name is the CATH domain name prefixed by 'CATH.', and the pdb_file_id is the NULL string.
   
         my @subsets_insert_values = ( "CATH.".$cath_dl->{domain_name}->[$j],
                                       '\N', $cath_dl->{pdb_id}->[$j], '\N',
   				    $subset_source_id, $class ) ;
   
         print {$fh->{subsets}} join("\t", @subsets_insert_values)."\n" ;
   
      }
      close($fh->{subsets}) ; #Close the subsets output file.
   
#mysql(SELECT domain_name, chain_id, segment_id, start_resno, end_resno FROM cath_domain)
   
      my $cath_dd ;
      ( $cath_dd->{domain_name},
        $cath_dd->{chain_id},
        $cath_dd->{segment_id},
        $cath_dd->{start_resno},
        $cath_dd->{end_resno} ) = 
         pibase::mysql_fetchcols($dbh, "SELECT domain_name, chain_id, segment_id, start_resno, end_resno FROM cath_domain_description") ;
#      pibase::mysql_fetchcols($dbh, "SELECT domain_name, chain_id, segment_id, start_resno, end_resno FROM cath_domall_boundaries") ;
   
#Iterate through the cath_domain_description SELECT results.
   
      foreach my $j ( 0 .. $#{$cath_dd->{domain_name}}) {
   
#Change undefined start and end residue numbers to '\N'
   
         if (! defined $cath_dd->{start_resno}->[$j]) {
            $cath_dd->{start_resno}->[$j] = '\N'; }
   
         if (! defined $cath_dd->{end_resno}->[$j]) {
            $cath_dd->{end_resno}->[$j] = '\N'; }
   
#Determine the unsigned integert parts of the start and end residue numbers.
   
         my $start_resno_int ;
         if ( $cath_dd->{start_resno}->[$j] ne ' ' ) {
            ($start_resno_int, undef) = residue_int($cath_dd->{start_resno}->[$j]) ; }
         if (! defined $start_resno_int) {
            $start_resno_int = '\N' ; }
   
         my $end_resno_int ;
         if ( $cath_dd->{end_resno}->[$j] ne ' ' ) {
            ($end_resno_int, undef) = residue_int($cath_dd->{end_resno}->[$j]) ; }
         if (! defined $end_resno_int) {
            $end_resno_int = '\N' ; }
   
#If the CATH chain_id is '0' (CATH's null chain_id convention), change to blank character (chain NULL conventions)
   
         if ( $cath_dd->{chain_id}->[$j] eq '0' ) {
            $cath_dd->{chain_id}->[$j] = ' ' ; }
            
#Designate the subset_details INSERT values. The subset name is the CATH domain name prefixed by 'CATH.'.
   
         my @subset_details_insert = ( "CATH.".$cath_dd->{domain_name}->[$j],
                                       $cath_dd->{segment_id}->[$j], '\N',
   				    $cath_dd->{chain_id}->[$j],
   				    $cath_dd->{start_resno}->[$j],
   				    $start_resno_int,
   				    $cath_dd->{end_resno}->[$j],
   				    $end_resno_int ) ;
   
#Display fields
   
         print {$fh->{subsets_details}} join("\t", @subset_details_insert)."\n" ;
   
      }
      close($fh->{subsets_details}) ; #Close the subsets details output file.
   
#mysql(SELECT class, arch, topol, homol, description FROM cath_names)
   
      my $cath_n ;
      ( $cath_n->{class},
        $cath_n->{arch},
        $cath_n->{topol},
        $cath_n->{homol},
        $cath_n->{description} ) =
         pibase::mysql_fetchcols($dbh, "SELECT class, arch, topol, homol, description FROM cath_names") ;
   
#Iterate through the cath_names SELECT results.
   
      foreach my $j ( 0 .. $#{$cath_n->{class}}) {
   
#Concatenate class, arch, topol, and homol with a '.' delimiter to form the class string.
   
         my $class = join('.', ( $cath_n->{class}->[$j],
                                 $cath_n->{arch}->[$j],
   			      $cath_n->{topol}->[$j],
   			      $cath_n->{homol}->[$j] ) ) ;
   
#Designate subsets_class insertion fields.
   
         my @insert_values = ( $subset_source_id, $class,
                               $cath_n->{description}->[$j] ) ;
   
#Display the insert values.
   
         print {$fh->{subsets_class}} join("\t", @insert_values)."\n" ;
      }
      close($fh->{subsets_class}) ; #Close the subsets_class output file.
   }

}

1;
