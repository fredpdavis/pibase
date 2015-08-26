=head1 NAME

pibase::PDB - module to handle PDB functions

=head1 DESCRIPTION

Handles general PDB functions: altloc check, filter

=head1 FILES

Operates on PDB file

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


package pibase::PDB ;
use strict;
use warnings;
use Carp ;

require Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT = qw/altloc_check altloc_filter pdb_clean_entries_idx/ ;

use pibase qw/locate_binaries/;
use File::Temp qw/tempdir tempfile/ ;


=head2 altloc_check()

   Title:       altloc_check()
   Function:    checks whether a PDB file for atoms with multiple locations
   Args:        $_ = pdb filename
   Return:      1 if containts multiple-occurrence atoms, 0 if not

=cut

sub altloc_check {

   my $pdb_fn = shift ;

   my $binaries = pibase::locate_binaries() ;

   if ($binaries->{'altloc_check'} eq 'ERROR') {
      croak("Error: altloc_check binary not found") ;
   }

   if (! -e $pdb_fn) {
      croak("Error: $pdb_fn not found") ; }

   my $altloc_fl = `$binaries->{altloc_check} < $pdb_fn` ;
   chomp $altloc_fl ;

   return $altloc_fl ;

}


=head2 altloc_filter()

   Title:       altloc_filter()
   Function:    calls the altloc_filter binary so that each atom occurs once
      Leaves the highest occupied (if occupancy defined) or the first location
      listed in the file.

   Args:        $_[0] = source pdb_filename
                $_[1] = output pdb_filename

   Returns:     nothing

=cut

sub altloc_filter {

   my ($pdb_fn, $outpdb_fn) = @_ ;

   my $binaries = pibase::locate_binaries() ;

   if ($binaries->{'altloc_filter'} eq 'ERROR') {
      croak("Error: altloc_filter binary not found") ;
   }

   if (! -e $pdb_fn) {
      croak("Error: $pdb_fn not found") ; }

   if ($pdb_fn eq $outpdb_fn) {
      croak("Error: Input PDB file is same as output: $pdb_fn") ; }

   system("$binaries->{altloc_filter} $pdb_fn > $outpdb_fn") ;

   return 1 ;
}


=head2 pdb_copy_entry_type()

   Title:       pdb_copy_entry_type()
   Function:    copies PDB pdb_entry_type for pibase import
   STDIN:       PDB pdb_entry_type
   STDOUT:      pdb_entry_type.pibase_id pibase table
   Returns:     nothing

=cut

sub pdb_copy_entry_type {

   my $in = shift ;
   my $tcom = "cp $in->{in_fn} $in->{out_fn}" ;
   system($tcom) ;

}


=head2 pdb_clean_entries_idx()

   Title:       pdb_clean_entries_idx()
   Function:    reformats the PDB entries.idx for pibase import
   STDIN:       PDB entries.idx
   STDOUT:      pibase.pdb_entries table
   Returns:     nothing

=cut

sub pdb_clean_entries_idx {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }


   my @headers = qw/idcode header date compound source author resolution experiment_type/ ;
   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }

   while (my $line = readline($fh->{in})) {

      if (($line =~ /^\#/) ||
          ($line =~ /^IDCODE/) ||
          ($line =~ /\-\-\-\-\-/)) {next;}
   
      chomp $line;

      my @tempfields = split(/\t/, $line) ;
      while ($#tempfields < 7) {
         my $newline = <STDIN> ; chomp $newline ;
         $line .= $newline ;
         @tempfields = split(/\t/, $line) ;
      }

      my ($pdb_id, $header, $raw_date, $compound, $source, $author,
          $resol, $exper_type) = @tempfields ;

      $resol =~ s/^\s*//;
      $resol =~ s/\s*$// ;

      if ($resol eq 'NOT' || $resol eq '') {
         $resol = '\N' ;
      } else {

         if ($resol =~ /\,/) { #eg PDB 2R24, 3HGN, 3KCO, 3KYY, 3L45 both X-ray and NeutD
            my @resols = split(/\,/, $resol) ;
            $resol = shift @resols ;
            map {if ($_ < $resol) {$resol = $_;}} @resols ;
         }

         $resol = sprintf("%.2f", $resol) ;
      } 

      $pdb_id = lc($pdb_id) ;

      my $date ;
      if ((defined $raw_date) && ($raw_date =~ /\//)) {
         $raw_date =~ s/^\s*// ;
         $raw_date =~ s/\s*$// ;

         my ($mon, $day, $year) = split(/\//, $raw_date) ;

         if ($year >= 71) {
            $year = '19'.$year ;
         } else {
            $year = '20'.$year ;
         }

         $date = $year.'-'.$mon.'-'.$day ;

      } else {

         $date = '\N' ;
      }

      my @outfields = ($pdb_id, $header, $date, $compound,
                       $source, $author, $resol, $exper_type) ;
  
      foreach my $j (0 .. $#outfields) {
         $outfields[$j] =~ s/^\s*// ;
         $outfields[$j] =~ s/\s*$// ;

         if ((!defined $outfields[$j]) ||
             ($outfields[$j] eq '')) {
             $outfields[$j] = '\N' ; }
      }
 
      print {$fh->{out}} join("\t", @outfields)."\n" ;
   }
   close($fh->{in}) ;
   close($fh->{out}) ;

}


=head2 pdb_clean_obsolete_dat()

   Title:       pdb_clean_obsolete_dat()
   Function:    reformats the PDB obsolete.dat for pibase import
   STDIN:       PDB obsolete.dat
   STDOUT:      pibase.pdb_obsolete table
   Returns:     nothing

=cut

sub pdb_clean_obsolete_dat {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }

   my @headers = qw/old_id new_id date/ ;
   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }

   my %months = (
      JAN => '01',
      FEB => '02',
      MAR => '03',
      APR => '04',
      MAY => '05',
      JUN => '06',
      JUL => '07',
      AUG => '08',
      SEP => '09',
      OCT => '10',
      NOV => '11',
      DEC => '12'
   ) ;
   

   while (my $line = readline($fh->{in})) {
      if ($line !~ /^OBSLTE/) {next;}

      chomp $line; chomp $line;
      my @f = split(' ', $line) ;

      my $old_pdb = lc($f[2]) ;
      my $new_pdb = '\N' ;

      if ($#f == 3) {
         $new_pdb = lc($f[3]) ;
         $new_pdb =~ s/ //g ;
      }

      if ($new_pdb eq '') {
         $new_pdb = '\N' ; }

      my $datestring = $f[1] ;
      my $day = substr($datestring, 0, 2 ) ;

      my $raw_mon = substr($datestring, 3, 3 ) ;
      my $month = $months{$raw_mon} ;

      my $year = substr($datestring, 7, 2 ) ;
      if ($year >= 71) {
         $year = '19'.$year ; }
      else {
         $year = '20'.$year ; }

      if (!defined $month) {
         print STDERR "$old_pdb $new_pdb MONTH : $month ($raw_mon)\n" ;}

      my $date = $year.'-'.$month.'-'.$day ;
  
      my @outvals = ($old_pdb, $new_pdb, $date) ;
      print {$fh->{out}} join("\t", @outvals)."\n" ;
   }
   close($fh->{in}) ;
   close($fh->{out}) ;

}


=head2 pdb_clean_release_date()

   Title:       pdb_clean_release_date()
   Function:    reformats the PDB release file for pibase import
   STDIN:       NOTDONE PDB obsolete.dat
   STDOUT:      pibase.pdb_release table
   Returns:     nothing

=cut

sub pdb_clean_release_date {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{in_fn}) {
      open($fh->{in}, "-") ;
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{in}, $in->{in_fn}) ;
      open($fh->{out}, ">".$in->{out_fn}) ;
   }

   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      my $fetch_external_data = fetch_external_data({quiet_fl => 1}) ;
      $pibase_specs = $fetch_external_data->{pibase_specs} ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }


   my @headers = qw/idcode release_date/ ;
   if (!exists $in->{header_fl} || $in->{header_fl} == 1) {
      print {$fh->{out}} "#".join("\t", @headers)."\n" ; }

   my $months = {
      "JAN" => '01' ,
      "FEB" => '02' ,
      "MAR" => '03' ,
      "APR" => '04' ,
      "MAY" => '05' ,
      "JUN" => '06' ,
      "JUL" => '07' ,
      "AUG" => '08' ,
      "SEP" => '09' ,
      "OCT" => '10' ,
      "NOV" => '11' ,
      "DEC" => '12'
   } ;


   my ($temp_fh, $temp_fn) = tempfile() ;
   while (my $line = readline($fh->{in})) {
      if ( ($line =~ /^\#/) ||
           ($line =~ /^IDCODE/) ||
           ($line =~ /^\-\-\-\-/) ) { next; }
   
      chomp $line;
      my ($pdb_id) = split(' ', $line) ;
      $pdb_id = lc($pdb_id) ;

#Grep the 'REVDAT   1' line from the pdb file and extract the PDB id and release date

      my $fileloc = get_pdb_filepath({pdb_id => $pdb_id,
                                      pibase_specs => $pibase_specs}) ;
      my $tcom ;
      if ($fileloc =~ /\.gz$/) {
         $tcom = $pibase_specs->{binaries}->{zcat} ;
      } else {
         $tcom = "cat" ;
      }
      if (!-s $fileloc) {
         print "#ERROR: ($pdb_id) can't find $fileloc\n";
         print STDERR "ERROR: ($pdb_id) can't find $fileloc\n";
      }
      $tcom .= " $fileloc 2>/dev/null | grep -m1 '^REVDAT   1' >$temp_fn" ;
      print "#$tcom\n" ;
      system($tcom) ;

      open(TREVDAT, $temp_fn) ; my $t_line = <TREVDAT> ; chomp $t_line;
      close(TREVDAT) ; unlink $temp_fn ;
      my @t = split(' ', $t_line) ;

# the pdb_id in the REVDAT entry is not always correct
# bugs found & reported 050516 (1b5x, 1bx0, 1bzj, 1cf8, 1ehy, 16bwh)
#	 my $pdb_id = $t[3] ; $pdb_id = lc($pdb_id) ;
# just use the pdb given in STDIN

      my $raw_date = $t[2] ;
      if (!defined $raw_date) {
         print {$fh->{out}} "$pdb_id\t\\N\n" ; next; }

      my $date ;

#If the date is defined and contains a '-',

      if ((defined $raw_date) && ($raw_date =~ /\-/)) {
#Remove leading and trailing spaces.
#raw_date: DD-mm-YY where mm is month in alpha

         $raw_date =~ s/^\s*// ;
         $raw_date =~ s/\s*$// ;
         my ($day, $mon, $year) = split(/\-/, $raw_date) ;

         if ($year >= 71) {
            $year = '19'.$year ;
         } else {
            $year = '20'.$year ;
         }

	 my $month = $months->{$mon} ;
	 if (! defined $month) {
	    print STDERR "ERROR: pdb $pdb_id unknown month $mon\n" ; next;}

         $date = $year.'-'.$month.'-'.$day ;
      } else {
         $date = '\N' ;
      }

      my @outvals = ($pdb_id, $date) ;
      foreach my $j (0 .. $#outvals) {
         $outvals[$j] =~ s/^\s*// ;
         $outvals[$j] =~ s/\s*$// ;

         if ((!defined $outvals[$j]) ||
	        ($outvals[$j] eq '')) {
            $outvals[$j] = '\N' ;
         }
      }
      print {$fh->{out}} join("\t", @outvals)."\n" ;
   }
   close($fh->{in}) ;
   close($fh->{out}) ;

}


=head2 pdb_clean_symop()

   Title:       pdb_clean_symop()
   Function:    Extracts and displays symmetry operators:
   STDIN:       PDB symop lines
   STDOUT:      pibase.pdb_release table
   Returns:     nothing

=cut

sub pdb_clean_symmop {

   while (my $line = <STDIN>) {
      chomp $line ;

      if ($line =~ /^REMARK 290     NNNMMM/) {

         my $keepread = 1;
         while ($keepread) {

            my $symline = <STDIN> ;

            if ($symline =~ /X/) {
   	       $symline =~ s/REMARK 290\s*// ;
   	       $symline =~ s/\s+//g ;
   	       my ($symop_name, $symdescr) = split(/ /, $symline) ;
	       print "$symop_name\t$symdescr\n" ;
	    } else {
   	       $keepread = 0;
            }
         }
         
      } elsif ($line =~ /^CRYST1/) {

#Extract and display symmetry operators:
#* a	7-16
#* b	16-25
#* c	25-34
#* alpha	34-40
#* beta	41-47
#* gamma	48-54
#* space group	56-66
#* z value	67-70

         my ($a, $b, $c, $alpha, $beta, $gamma, $sgroup, $zval) =
            ( substr($line, 6, 9), substr($line, 15, 9), substr($line, 24, 9), 
               substr($line, 33, 7), substr($line, 40, 7), substr($line, 47, 7),
	       substr($line, 55, 11), substr($line, 66,4)) ;

         my @outfields = ($a, $b, $c, $alpha, $beta, $gamma, $sgroup) ;
         print join("\t", @outfields)."\n" ;
      }
   }
}


=head2 get_pdb_filepath()

   Title:       get_pdb_filepath()
   Function:    returns file path to a PDB entry
   Args:	->{pdb_id} = pdb identifier
                [->{pibase_specs} = $pibase_specs] - optional
   Returns:     returns PDB entry filepath

=cut

sub get_pdb_filepath {
   my $in = shift ;

   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::complete_pibase_specs()
   } else {
      $pibase_specs = $in->{pibase_specs}; }

   my $pdbpath = $pibase_specs->{pdb_dir} ;
   if ($pibase_specs->{external_data}->{pdb}->{file_layout} eq 'wwpdb') {
      $pdbpath .= '/'.substr($in->{pdb_id},1,2) ;
   }
   $pdbpath .= '/pdb'.lc($in->{pdb_id}).'.ent';


   if (exists $pibase_specs->{external_data}->{pdb}->{compress_fl} &&
       $pibase_specs->{external_data}->{pdb}->{compress_fl} == 1) {
      $pdbpath .= '.gz' ; }

   return $pdbpath ;
}


=head2 nmr_model1_extractor

   Title:       nmr_model1_extractor()
   Function:    extracts first model from NMR PDB files and moves to 
                $specs->{pdbnmr_dir}
   Args:        none
   Returns:     nothing
   Tables in:   pibase.pdb_entries
   Files in:    foreach PDB NMR entry: <$specs->{pdb_dir}>/pdb<$pdb_id>.ent
   Files out:   foreach PDB NMR entry: <$specs->{pdbnmr_dir}>/<$pdb_id>_1.ent
   NOTE: expects PDBs to be stored as pdb<$pdb_id>.ent in one raw directory

=cut

sub nmr_model1_extractor {

#pibase2010 - moved to a wwPDB directory layout, and compressed files.

   require DBI ;

   my $movethreshold = 500 ;

   my $in = shift ;

   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = $in->{pibase_specs} ;
   } else {
      $pibase_specs = pibase::get_specs() ;
   }
   my ($dbh) = pibase::connect_pibase() ;
   my ($nmr_pdbs) = pibase::mysql_fetchcols($dbh,
      "SELECT pdb_id FROM pdb_entry_type where experiment_type = \"NMR\"") ;

   my $temp_dir = tempdir(CLEANUP => 1) ;
   chdir $temp_dir ;
   my (@movethese, @movedest) ;

   my $j = 1 ;
   print STDERR "now on: 0" ;
   foreach my $pdb_id ( @{$nmr_pdbs}) {
      print STDERR "\b"x(length($j - 1)) ; print STDERR $j ;
      chomp $pdb_id;
      $pdb_id = lc($pdb_id) ;

      my $fn ;
      $fn->{pdb_o} = get_pdb_filepath({pdb_id => $pdb_id,
                                       pibase_specs => $pibase_specs}) ;
      $fn->{pdb_new} = $pdb_id.'_1.ent' ;

      if (! -s $fn->{pdb_o}) {
         print STDERR "ERROR ($pdb_id): $fn->{pdb_o} not found\n" ; }

      my $compress_fl = 0 ;
      my $tcom_catold  ;
      if ($fn->{pdb_o} =~ /\.gz$/) {
         $compress_fl = 1 ;
         $tcom_catold = $pibase_specs->{binaries}->{zcat}.' '.$fn->{pdb_o}.' '.
            '2>/dev/null ' ;
      } else {
         $tcom_catold = "cat ".$fn->{pdb_o}." " ;
      }
      push @movethese, $fn->{pdb_new}.'.gz' ;
      push @movedest, $pibase_specs->{pdbnmr_dir}."/".substr($pdb_id,1,2) ;

      if (! -s $fn->{pdb_new}.'.gz') {
         my $tcom = $tcom_catold." | sed \'/^ENDMDL/q\' > $fn->{pdb_new}" ;
         print "# $tcom\n" ;
         system($tcom) ;
         system("gzip ".$fn->{pdb_new}) ;
      } else {
         print "ERROR ($pdb_id): $fn->{pdb_new} already exists\n";
      }

      if ($#movethese == $movethreshold) {
         foreach my $k (0 .. $#movethese) {
            my $t_fn = $movethese[$k] ;
            print "#mv $t_fn ".$movedest[$k]."\n" ;
            pibase::safe_move($t_fn, $movedest[$k]);
         }
         @movethese = () ;
         @movedest= () ;
      }
      $j++ ;
   }

   foreach my $k (0 .. $#movethese) {
      my $t_fn = $movethese[$k] ;
      print "#mv $t_fn ".$movedest[$k]."\n" ;
      pibase::safe_move($t_fn, $movedest[$k]);
   }

   print STDERR "\n" ;

}


1 ;
