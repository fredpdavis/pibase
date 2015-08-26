=head1 NAME

pibase::data::external::PISA - perl interface to PISA routines

=head1 DESCRIPTION

Perl package that provides interface to PISA data routines

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

package pibase::data::external::PISA ;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT = qw/bdp_path_2_id get_pisa_files/ ;

use pibase qw/connect_pibase mysql_hashload mysql_fetchcols mysql_hasharrload safe_move/;
use File::Temp qw/tempfile tempdir/ ;
use File::Path qw/mkpath/ ;
use File::Basename;
use POSIX qw/ceil/ ;

=head2 get_pisa_files

   Title:       get_pisa_files
   Function:    downloads entries in PISA index.txt from the EBI ftp site
      ftp://ftp.ebi.ac.uk/pub/databases/msd/pisa/data/iv/1ivo/1ivo.pdb.gz
   Args:        none
   Returns:     nothing
   Comand line: $ARGV[0] = PISA index.txt file
   Files in:    PISA index.txt file ($ARGV[0])
   STDOUT:      "File ".<PISA file path>." ".(exists|failed|downloaded)

=cut

sub get_pisa_files {

   my $pisa_url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/pisa/data/" ;

   my $pibase_specs = pibase::get_specs() ;

   open(INDEX, $pibase_specs->{external_data}->{pisa}->{files}->{"index"}) ;
   while (my $line = <INDEX>) {
      chomp $line;
      if ($line =~ /^\#/) {next;}
      my ($name, $entry_type, $homohet, $no_chains, $no_res, undef) =
         split(' ', $line) ;
      if ($entry_type eq 'ERROR' || $no_chains > 
          $pibase_specs->{external_data}->{pisa}->{max_num_chains}){ next;}

      print $name;
      chomp $name;

      my $pdb_id = substr($name,0,4) ;
      my $orig_name = $name.".pdb.gz" ;
      my $local_name = $name."_pisa.pdb.gz" ;
      my $full_dir = $pibase_specs->{pisa_dir}.'/'.substr($pdb_id,1,2) ;
      if (! -s $full_dir) {
         mkpath($full_dir);
         print STDERR "making directory $full_dir\n" ;
      }

      my $cur_url = $pisa_url.substr($pdb_id,1,2)."/$pdb_id/".$orig_name ;

#-- Check if file exists
      if ( -e "$full_dir/$local_name" ){
         print "File $full_dir/$local_name exists\n";
      } else {
         my $wget_com = "wget -q $cur_url -O $full_dir/$local_name";
         print "#$wget_com\n";
         my $result = system($wget_com) ;
         if ( $result ){
            print "File $full_dir/$local_name failed\n";
         } else {
            print "File $full_dir/$local_name downloaded\n";
         }
      }
   }
   close(INDEX);

}



=head2 check_pisa_files

   Title:       check_pisa_files()
   Function:    check integrity of PISA files (looks for ^END|^TER on last line)
   Args:        none
   Returns:     nothing
   Input file:  PISA index.txt file (basename only, not fullpath)
   STDOUT:      (incomplete|complete).": ".<PISA file basename>

=cut

sub check_pisa_files {

   my $in = shift ;
   my $fh ;
   if (!exists $in->{out_fn}) {
      open($fh->{out}, ">-") ;
   } else {
      open($fh->{out}, ">".$in->{out_fn}) ;
   }

   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      my $fetch_external_data = fetch_external_data({quiet_fl => 1}) ;
      $pibase_specs = $fetch_external_data->{pibase_specs} ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }

   open(INDEX, $pibase_specs->{external_data}->{pisa}->{files}->{"index"}) ;
   while (my $line = <INDEX>) {
      chomp $line;
      if ($line =~ /^\#/) {next;}
      my ($name, $entry_type, $homohet, $no_chains, $no_res, undef) =
         split(' ', $line) ;
      if ($entry_type eq 'ERROR' || $no_chains > 
          $pibase_specs->{external_data}->{pisa}->{max_num_chains}){ next;}

      my $pdb_id = substr($name,0,4) ;
      my $orig_name = $name.".pdb.gz" ;
      my $local_name = $name."_pisa.pdb.gz" ;
      my $pisa_dir = $pibase_specs->{pisa_dir}.'/'.substr($pdb_id,1,2) ;

      if (!-s "$pisa_dir/$local_name") {
         print {$fh->{out}} "$local_name\tnot found\n" ;
         next;
      }

      my $complete_fl =
         `zcat $pisa_dir/$local_name | tail -1 | grep -c '^[END\|TER]'` ;

#complete if :
# * file ends in END or TER (majority)
# note though, falsely calls some really complete ones
#  * file with an OXT ATOM record, eg 1htz_1.mmol
#  * file has HETATM records after TER-terminated ATOM records - eg 1nmc_1.mmol

      if ($complete_fl == 0 ) { print {$fh->{out}} "$local_name\tincomplete\n";}
      else { print {$fh->{out}} "$local_name\tcomplete\n" ; }
   }
   close(INDEX) ;
   close($fh->{out}) ;

}

1 ;
