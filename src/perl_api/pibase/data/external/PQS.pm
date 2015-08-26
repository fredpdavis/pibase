=head1 NAME

pibase::data::external::PQS - perl interface to PQS routines

=head1 DESCRIPTION

Perl package that provides interface to PQS data routines

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

package pibase::data::external::PQS ;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT = qw/bdp_path_2_id/ ;

use pibase qw/connect_pibase mysql_hashload mysql_fetchcols mysql_hasharrload safe_move/;
use File::Temp qw/tempfile tempdir/ ;
use File::Path qw/mkpath/ ;
use File::Basename;
use POSIX qw/ceil/ ;

=head2 get_pqs_files

   Title:       get_pqs_files
   Function:    downloads entries in the PQS LIST files from the EBI website
      http://pqs.ebi.ac.uk/pqs-doc/macmol
   Args:        none
   Returns:     nothing
   Comand line: $ARGV[0] = PQS LIST file
   Files in:    PQS LIST file ($ARGV[0])
   STDOUT:      "File ".<PQS file path>." ".(exists|failed|downloaded)

=cut

sub get_pqs_files {

   my $pqsurl = "http://pqs.ebi.ac.uk/pqs-doc/macmol/";

   my $specs = pibase::get_specs() ;
   my $pqs_dir = $specs->{pqs_dir} ;

   if ($#ARGV < 0) { die "usage: ".__FILE__." LIST;"}

   open(LIST, $ARGV[0]) ;
   while ( my $list = <LIST> ){
      if ($list =~ /^\#/) {next;}
      my ($name) = ($list =~ /(\w+\.\w+)/) ;
      $list = $name ;
      if ($list =~ /water$/) {next;}
      if ($list =~ /23fold$/) {next;}
      if ($list =~ /5fold$/) {next;}
      print $list;
      chomp $list;

      my ($prot, $dir, $ext) = fileparse($list, '\..*');

   #-- Check if file exists
      if ( -e "$pqs_dir/$prot$ext" ){
         print "File $pqs_dir/$prot$ext exists\n";
      } else {
         my $result = system("lynx -source $pqsurl/$prot$ext > $pqs_dir/$prot$ext");
         if ( $result ){
            print "File $pqs_dir/$prot$ext failed\n";
         } else {
            print "File $pqs_dir/$prot$ext downloaded\n";
         }
      }
   }
   close(LIST);

}


=head2 list_pqs_file_urls

   Title:       list_pqs_file_urls
   Function:    Makes list of URLS and mv commands to download PQS files
      in the rcsb-esque substr(pdb_id,1,2)/ from 
      http://pqs.ebi.ac.uk/pqs-doc/macmol

      Makes directory for final resting place as needed (mkpath)

   Args:        none
   Returns:     nothing
   Comand line: $ARGV[0] = PQS LIST file
   Files in:    uses PQS LIST file
      $pibase_specs->{external_data}->{pqs}->{files}->{list}

   STDOUT:      PQS file URL
                mv PQS_file substr(pdb_id,1,2)

=cut


sub list_pqs_file_urls {
   print STDERR "GOT HERE\n" ;

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

   my $pqsurl = "http://pqs.ebi.ac.uk/pqs-doc/macmol/";

   if ($#ARGV < 0) { die "usage: ".__FILE__." LIST;"}

   open(LIST, $pibase_specs->{external_data}->{pqs}->{files}->{list}) ;
   while (my $list = <LIST>) {
      if ($list =~ /^\#/) {next;}
      my ($name) = ($list =~ /(\w+\.\w+)/) ;
      $list = $name ;
      if ($list =~ /water$/) {next;}
      if ($list =~ /23fold$/) {next;}
      if ($list =~ /5fold$/) {next;}
      chomp $list;

      my ($prot, $dir, $ext) = fileparse($list, '\..*');
      my $pdb_id = substr($prot,0,4) ;

#-- Check if deposit subdirectory exists - if not, create
      my $pqs_dir = $pibase_specs->{pqs_dir}.'/'.substr($pdb_id,1,2) ;
      if (! -s $pqs_dir) {
         mkpath($pqs_dir);
         print STDERR "making directory $pqs_dir\n" ;
      }

#-- If file already exists, skip printing the URL
      if ( -e "$pqs_dir/$prot$ext" ) {
         print STDERR "File $pqs_dir/$prot$ext exists\n";
      } else {
         my $cur_url = $pqsurl.$prot.$ext ;
#         print {$fh->{out}} "-P $pqs_dir $cur_url\n";
         print {$fh->{out}} $cur_url."\n";
         print {$fh->{out}} "mv $prot$ext ".substr($pdb_id,1,2)."\n" ;
      }
   }
   close(LIST) ;

}

=head2 check_pqs_files

   Title:       check_pqs_files()
   Function:    check integrity of PQS files (looks for ^END|^TER on last line)
   Args:        none
   Returns:     nothing
   STDIN:       PQS FILE list (basename only, not fullpath)
   STDOUT:      (incomplete|complete).": ".<PQS file basename>

=cut

sub check_pqs_files {

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

   open(LIST, $pibase_specs->{external_data}->{pqs}->{files}->{list}) ;
   while (my $file = <LIST>) {
      chomp $file;
      if ($file =~ /water$/) {next;}
      if ($file =~ /23fold$/) {next;}
      if ($file =~ /5fold$/) {next;}

      my $pdb_id = substr($file,0,4) ;
      my $pqs_dir = $pibase_specs->{pqs_dir}.'/'.substr($pdb_id,1,2) ;

      if (!-s "$pqs_dir/$file") {
         print {$fh->{out}} "$file\tnot found\n" ;
         next;
      }

      my $complete_fl = `tail -1 $pqs_dir/$file | grep -c '^[END\|TER]'` ;

#complete if :
# * file ends in END or TER (majority)
# note though, falsely calls some really complete ones
#  * file with an OXT ATOM record, eg 1htz_1.mmol
#  * file has HETATM records after TER-terminated ATOM records - eg 1nmc_1.mmol

      if ($complete_fl == 0 ) { print {$fh->{out}} "$file\tincomplete\n" ;}
      else { print {$fh->{out}} "$file\tcomplete\n" ; }
   }
   close(LIST) ;
   close($fh->{out}) ;

}

1 ;
