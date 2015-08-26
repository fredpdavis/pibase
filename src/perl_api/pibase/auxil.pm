=head1 NAME

pibase::aux - perl interface to auxiliary pibase routines

=head1 DESCRIPTION

Perl package that has miscellaneous pibase routines for non-core functions

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

package pibase::aux ;
require Exporter;
@ISA = qw/Exporter/ ;
@EXPORT = qw/bdp_path_2_id/ ;

use strict;
use warnings;
use pibase qw/connect_pibase mysql_hashload mysql_fetchcols mysql_hasharrload safe_move sid_2_domdir/;
use File::Temp qw/tempfile tempdir/ ;
use File::Path qw/mkpath/ ;
use POSIX qw/ceil/ ;


=head2 setup_clusterrun()

   Title:       setup_clusterrun()
   Function:    Splits input file and sets up tasklist for SGE script.
   Args:        $_[0] = name of inputfile to split into chunks.
                $_[1] = name of outputfile for tasklist
   Returns:     nothing
   STDIN:       is interactive

=cut

sub setup_clusterrun {

   my $usage = "./setuprun.pl inputfile outputfile" ;
   my $inputlist = $ARGV[0] || die $usage;
   my $outputfile = $ARGV[1] || die $usage ;
   my $numlines = `wc -l $inputlist` ; chomp $numlines;
   $numlines =~ s/^ *// ; $numlines =~ s/ .*$//g ;
   print "Number of lines: $numlines\n" ;
   print "How many tasks? " ; my $tasks = <STDIN> ; chomp $tasks ;

   my $splitlines ;
   if (isnumeric($tasks) && ($tasks >= 2)) {
      $splitlines = ceil($numlines / $tasks);
      print "split on $splitlines\n" ;
   } else {
      die "number of tasks is invalid: $tasks\n" ; }

   my $splitpre = 'split.'.$inputlist.'.' ;
   print "split prefix of: $splitpre\n" ;
   print "is this ok? (y|n)"; my $split_ok = <STDIN>; chomp $split_ok ;

   if ($split_ok eq 'n') {
      print "   so what do you want? " ;
      $splitpre = <STDIN> ; chomp $splitpre ; }

   my $tcom = "split -l $splitlines $inputlist $splitpre" ;
   system($tcom) ;

   my @splitfiles ;
   open (SPLITLIST, "ls $splitpre* |") ;
   open (OUTF, ">$outputfile") ;
   print OUTF "set tasks1=(" ;
   while (my $line = <SPLITLIST>) {
      chomp $line;
      print OUTF " '".$line."'" ;
   }
   print OUTF " )" ;
   close (OUTF) ;
   close(SPLITLIST) ;

}


=head2 isnumeric()

   Title:       isnumeric()
   Function:    Checks if a scalar value is numeric or not
   Args:        $_ - scalar value
   Returns:     1 if number, 0 if not

=cut

sub isnumeric {
   my $a = shift ;
   if ($a =~ /^[0-9]+$/) {
      return 1;
   } else {
      return 0 ; }
}


=head2 denuller()

   Title:       denuller()
   Function:    Replaces null values (NULL) in a tab-delimited file with tab.
            NOT SURE
   Args:        $_ - scalar value
   Returns:     1 if number, 0 if not

=cut

sub denuller {

   while (my $line = <STDIN>) {
      $line =~ s/^\t/ \t/g ;
      $line =~ s/\t$/\t /g ;
      $line =~ s/(?<=\t)(?=\t)/ /g ;

      $line =~ s/^NULL\t/\t/ ;
      $line =~ s/(?<=\t)NULL(?=\t)//g ;
      $line =~ s/NULL$// ;
      print $line ;
   }

}

1;
