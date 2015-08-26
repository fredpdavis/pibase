=head1 NAME

pibase::SCOP - package that handles SCOP release files and pibase.

=head1 DESCRIPTION

Processes SCOP release files and 

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

package pibase::SCOP ;
use strict;
use warnings;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/altloc_check altloc_filter/ ;
use Exporter;
use Carp ;

use pibase qw/locate_binaries/;
use pibase::residue_math qw/residue_int/;
use File::Temp qw/tempfile/ ;


=head2 scop_clean_cla()

   Title:       scop_clean_cla()
   Function:    Processes the SCOP .cla file for pibase import
   STDIN:       SCOP CLA file
   STDOUT:      pibase.scop_cla table
   Args:        none
   Returns:     none

=cut

sub scop_clean_cla {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }

   my @headers = qw/sid pdb_id chain start_resno_pdb end_resno_pdb class fold superfam fam cl cf sf fa dm sp px/ ;
   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }

   while (my $line = readline($fh->{in})) {

      if ($line =~ /^\#/) {next;}
      chomp $line;

      my ($sid, $pdb_id, $raw_domaindef, $sccs, $raw_px_id, $raw_sunid ) =
         split(/\t/, $line) ;

      my ($der_class, $der_fold, $der_superfam, $der_fam) = split(/\./, $sccs) ;

      my @raw_domaindef = split(',', $raw_domaindef) ;
      my (@cur_chain, @cur_start, @cur_end) ;

      foreach my $j (0 .. $#raw_domaindef) {

#If the segment domain definition contains ':',

         if ($raw_domaindef[$j] =~ /:/) {

#Store everything before ':' as the current chain, and remove it, as well as the ':' delimiter, from the segment domain definition.

            ($cur_chain[$j]) = ($raw_domaindef[$j] =~ /^(.+)\:/ ) ;
            $raw_domaindef[$j] =~ s/^.+\:// ;
         } else {
   	       $cur_chain[$j] = ' ' ;
         }
   
#If the domain definition is '' or '-', set the start and end residue to ' '.

         if (($raw_domaindef[$j] eq '') ||
	     ($raw_domaindef[$j] eq '-')) {

   	    $cur_start[$j] = $cur_end[$j] = ' ' ;
   	 } else {

            $raw_domaindef[$j] =~ s/^\s*// ;
            $raw_domaindef[$j] =~ s/\s*$// ;

#If the domain definition contains '-', extract the start and stop residues.
#Otherwise, set the start and end residues to the contents of the domain definition. (only 1 residue in the segment)

            if ($raw_domaindef[$j] =~ /-/) {
               ($cur_start[$j], $cur_end[$j]) =
                  ($raw_domaindef[$j] =~ /(-?[0-9]+.?)-(-?[0-9]+.?)/) ;
            } else {
	          $cur_start[$j] = $cur_end[$j] = $raw_domaindef[$j] ;
            }

#If start or stop residues are undefined, set them to the NULL string.

            if (!(defined $cur_start[$j])) {
               $cur_start[$j] = '\N' ; }

            if (!(defined $cur_end[$j])) {
               $cur_end[$j] = '\N' ; }

         }
      }

#Extract cl, cf, sf, fa, dm, sp, and px sun ids from the sun id string.

      my ($cl_id, $cf_id, $sf_id, $fa_id, $dm_id, $sp_id, $px_id) =
( $raw_sunid =~ /cl=(.+),cf=(.+),sf=(.+),fa=(.+),dm=(.+),sp=(.+),px=(.+)/ ) ;

#Designate the segment-independent output fields: sun id, PDB id, class id, fold id, superfamily id, family id, px, cl, cf, sf, fa, dm and sp sun id's

      my @outfields_pre = ($sid, $pdb_id, $der_class, $der_fold,
	                   $der_superfam, $der_fam, $cl_id, $cf_id,
			   $sf_id, $fa_id, $dm_id, $sp_id, $px_id) ;

      foreach my $j (0 .. $#outfields_pre) {
         if ((! (defined $outfields_pre[$j])) ||
	        ($outfields_pre[$j] eq '-')) {
               $outfields_pre[$j] = '\N' ;
         } else {
               $outfields_pre[$j] =~ s/^\s*// ;
               $outfields_pre[$j] =~ s/\s*$// ;
         }
      }

      foreach my $j (0 .. $#cur_start) {
         print {$fh->{out}} join("\t", @outfields_pre[0..1])."\t" ;
         print {$fh->{out}} join("\t", ($cur_chain[$j],$cur_start[$j],
                                        $cur_end[$j]))."\t";
         print {$fh->{out}} join("\t", @outfields_pre[2..12])."\n" ;
      }

   }
   close($fh->{in}) ;
   close($fh->{out}) ;

}


=head2 scop_clean_des()

   Title:       scop_clean_des()
   Function:    Processes the SCOP .des file for pibase import
   STDIN:       SCOP DES file
   STDOUT:      pibase.scop_des table
   Args:        none
   Returns:     none

=cut

sub scop_clean_des {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }


   my @headers = qw/sunid entry_type class_id fold_id superfam_id fam_id scop_id description/ ;
   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }

   while (my $line = readline($fh->{in})) {

      if ($line =~ /^\#/) {next;}
   
      chomp $line;
      
      my ($sunid, $entry_type, $sccs, $scop_id, $description) =
         split(/\t/, $line) ;
 
      my ($der_class, $der_fold, $der_superfam, $der_fam) = split(/\./, $sccs) ;
      my @outfields =  ($sunid, $entry_type, $der_class, $der_fold,
                        $der_superfam, $der_fam, $scop_id, $description) ;
  
      foreach my $j (0 .. $#outfields) {
         if ((! (defined $outfields[$j])) ||
	     ($outfields[$j] eq '-')) {

#If the field is fold, superfam or fam_id, set to '0'.

            if ( ($j >= 3) && ($j < 6) ) {
               $outfields[$j] = '0' ;
            } else {
               $outfields[$j] = '\N' ;
            }

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


=head2 scop_clean_hie()

   Title:       scop_clean_hie()
   Function:    Processes the SCOP .hie file for pibase import
   STDIN:       SCOP HIE file
   STDOUT:      pibase.scop_hie table
   Args:        $_->{in_fn} = input file
                $_->{out_fn} = output file
                $_->{header_fl} = flag to generate header in output (default 1)
   Returns:     none

=cut

sub scop_clean_hie {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }

   my @headers = qw/self_sunid parent_sunid kids_sunid/ ;
   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }

   while (my $line = readline($fh->{in})) {
      if ($line =~ /^\#/) {next;}
      chomp $line;
      my ($self_sunid, $parent_sunid, $kids_sunid) = split(/\t/, $line) ;

      my @outfields = ($self_sunid, $parent_sunid, $kids_sunid) ;

      foreach my $j (0 .. $#outfields) {

         if ((! (defined $outfields[$j])) ||
	     ($outfields[$j] eq '-')) {
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


=head2 pibase_import_scop_domains()

   Title:       pibase_import_scop_domains()
   Function:    Processes the SCOP .hie file for pibase import
   Tables in:   pibase.scop_cla
                pibase.scop_des

   Tables out:  pibase.subsets
                pibase.subsets_class
                pibase.subsets_details
   Args:        none
   Returns:     none

=cut

sub pibase_import_scop_domains {
   require DBI ;

   my $in = shift ;
   my $overwrite_fl = 0;
   if (exists $in->{overwrite_fl}) {$overwrite_fl = $in->{overwrite_fl};}

# only allow one version of each particular domain def
   my ($dbh, $pibase) = pibase::connect_pibase() ;
   my $pibase_specs = $in->{pibase_specs} ;

   pibase::mysql_runcom($dbh,
      'DELETE FROM subsets_source WHERE subset_source = "scop"') ;

   pibase::mysql_runcom($dbh,
      'INSERT INTO subsets_source'.
      '(subset_source_id, subset_source, version) '.
      'VALUES("", "scop","'.
      $in->{pibase_specs}->{external_data}->{scop}->{ver}.'") ;');

   my ($subset_source_id) = pibase::mysql_singleval($dbh,
      'SELECT subset_source_id FROM subsets_source '.
      'WHERE subset_source = "scop"') ;

#Set subsets, subsets_class, and subsets_details output files.

   if (!-s $pibase_specs->{external_data}->{scop}->{processed_files}->{"subsets.scop"} || $overwrite_fl == 1) {
      my $fh ;
      open($fh->{subsets}, ">".
   $pibase_specs->{external_data}->{scop}->{processed_files}->{"subsets.scop"});

      open($fh->{subsets_class}, ">".
   $pibase_specs->{external_data}->{scop}->{processed_files}->{"subsets_class.scop"});

      open($fh->{subsets_details}, ">".
   $pibase_specs->{external_data}->{scop}->{processed_files}->{"subsets_details.scop"});
   
#mysql(SELECT scop_id, pdb_id, class_id, fold_id, superfam_id, fam_id FROM scop_cla)
   
      my $scop_cla ;
      ( $scop_cla->{scop_id},
        $scop_cla->{pdb_id},
        $scop_cla->{class_id},
        $scop_cla->{fold_id},
        $scop_cla->{superfam_id},
        $scop_cla->{fam_id} ) =
         pibase::mysql_fetchcols($dbh,"SELECT distinct scop_id, pdb_id, ".
            "class_id, fold_id, superfam_id, fam_id FROM scop_cla") ;
   
      foreach my $j (0 .. $#{$scop_cla->{scop_id}}) {
   
#Create a class string by concatenating the SCOP class, fold_id, superfam_id, and fam_id ids with a '.' delimiter.
   
         my $class = join('.', ( $scop_cla->{class_id}->[$j],
                                 $scop_cla->{fold_id}->[$j],
                                 $scop_cla->{superfam_id}->[$j],
                                 $scop_cla->{fam_id}->[$j] ) ) ;
   
#Designate and display output fields. The subset name is the SCOP domain name prefixed by 'SCOP.', and the pdb_file_id is the NULL string.
   
         my @subsets_insert_values = ( "SCOP.".$scop_cla->{scop_id}->[$j],
                                       '\N', $scop_cla->{pdb_id}->[$j], '\N',
   				    $subset_source_id, $class ) ;
   
         print {$fh->{subsets}} join("\t", @subsets_insert_values)."\n" ;
   
      }
      close($fh->{subsets});
   
#mysql(SELECT scop_id, chain_id, start_resno, end_resno FROM scop_cla)
   
      $scop_cla = {} ;
   
      ( $scop_cla->{scop_id},
        $scop_cla->{chain_id},
        $scop_cla->{start_resno},
        $scop_cla->{end_resno} ) =
         pibase::mysql_fetchcols($dbh, "SELECT scop_id, chain_id, start_resno, ".
         "end_resno FROM scop_cla") ;
   
      my %domain_seen ;
      foreach my $j (0 .. $#{$scop_cla->{scop_id}}) {
         if (!(defined $scop_cla->{chain_id}->[$j])) {
            $scop_cla->{chain_id}->[$j] = '\N'; }
         if (!(defined $scop_cla->{start_resno}->[$j])) {
            $scop_cla->{start_resno}->[$j] = '\N'; }
         if (!(defined $scop_cla->{end_resno}->[$j])) {
            $scop_cla->{end_resno}->[$j] = '\N'; }
   
         $domain_seen{$scop_cla->{scop_id}->[$j]}++ ;
   
         my $start_resno_int ;
         if ( $scop_cla->{start_resno}->[$j] ne ' ' ) {
            ($start_resno_int, undef) =
               residue_int($scop_cla->{start_resno}->[$j]) ;
         }
   
         if (! defined $start_resno_int) {
            $start_resno_int = '\N' ; }
   
         my $end_resno_int ;
         if ( $scop_cla->{end_resno}->[$j] ne ' ' ) {
            ($end_resno_int, undef) =
               residue_int($scop_cla->{end_resno}->[$j]) ;}

         if (! defined $end_resno_int) {
            $end_resno_int = '\N' ; }
   
#Designate the subset_details INSERT values. The subset name is the SCOP domain name prefixed by 'SCOP.'.
   
         my @subset_details_insert = ( "SCOP.".$scop_cla->{scop_id}->[$j],
      	                            $domain_seen{$scop_cla->{scop_id}->[$j]},
   				    '\N', $scop_cla->{chain_id}->[$j],
   				    $scop_cla->{start_resno}->[$j],
   				    $start_resno_int,
   				    $scop_cla->{end_resno}->[$j],
   				    $end_resno_int ) ;
         print {$fh->{subsets_details}} join("\t", @subset_details_insert)."\n" ;
   
      }
      close($fh->{subsets_details});
   
#mysql(SELECT class_id, fold_id, superfam_id, fam_id, description FROM scop_des)
   
      my $scop_des ;
      ( $scop_des->{class_id},
        $scop_des->{fold_id},
        $scop_des->{superfam_id},
        $scop_des->{fam_id},
        $scop_des->{description} ) =
         pibase::mysql_fetchcols($dbh, 'SELECT class_id, fold_id, superfam_id,'.
                                       ' fam_id, description FROM scop_des'.
   				    ' WHERE entry_type = "cl" OR '.
   				    'entry_type = "cf" OR entry_type = "sf"'.
   				    ' OR entry_type  = "fa"') ;
   
      foreach my $j (0 .. $#{$scop_des->{class_id}}) {
         my $class = join('.', ( $scop_des->{class_id}->[$j],
                                 $scop_des->{fold_id}->[$j],
                                 $scop_des->{superfam_id}->[$j],
                                 $scop_des->{fam_id}->[$j] ) ) ;
   
         my @insert_values = ( $subset_source_id, $class,
                               $scop_des->{description}->[$j] ) ;
   
         print {$fh->{subsets_class}} join("\t", @insert_values)."\n" ;
      }
      close($fh->{subsets_class});
   }

}

1;
