=head1 NAME

pibase::data::external::ASTRAL - perl module of ASTRAL data routines

=head1 DESCRIPTION

Perl package that provides interface to ASTRAL data files. Includes routine
to run MAFFT on domain sequences in the case that the ASTRAL alignments
are yet to be released for a new SCOP version.

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


package pibase::data::external::ASTRAL ;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT = qw/bdp_path_2_id/ ;

use pibase qw/connect_pibase mysql_hashload mysql_fetchcols mysql_hasharrload safe_move sid_2_domdir get_specs/;
use File::Temp qw/tempfile tempdir/ ;
use File::Path qw/mkpath/ ;
use File::Basename;
use POSIX qw/ceil/ ;

=head2 get_sf_fam_seq_astral()

   Title:       get_sf_fam_seq_astral()
   Function:    Create files of amino acid sequences for each SCOP family and
                  superfamily
   Args:        none
   Returns:     nothing
   Files out:   foreach SCOP family: <ASTRAL fam aln directory>/<family>.fa (fasta)
                foreach SCOP superfamily: <ASTRAL fam aln directory>/<family>.fa

=cut

sub get_sf_fam_seq_astral {

   my ($fn,$fh);
   my $pibase_specs = pibase::get_specs() ;
   $fn->{fam_seq} =  $pibase_specs->{astral}->{fam_seq} ;
   $fn->{sf_seq} =  $pibase_specs->{astral}->{sf_seq} ;
   $fn->{gd_seq} =  $pibase_specs->{astral}->{gd_seq} ;

   open(GDSEQF, $fn->{gd_seq}) ;
   my ($lastfam, $lastsf)  = ('', '');
   while (my $line = <GDSEQF>) {
      if ($line =~ /^>/) {
         my (undef, $fam, undef, undef) = split(/\s/, $line) ;
         my ($sf) = ($fam=~ /([a-z]\.[0-9]+\.[0-9]+)\./) ;

         if ($fam ne $lastfam) {
            if (exists $fh->{fam_seq}) {close $fh->{fam_seq};}
            open ($fh->{fam_seq}, ">>$fn->{fam_seq}/$fam.fa") ;
            $lastfam = $fam ;
         }

         if ($sf ne $lastsf) {
            if (exists $fh->{sf_seq}) {close $fh->{sf_seq};}
            open ($fh->{sf_seq}, ">>$fn->{sf_seq}/$sf.fa") ;
            $lastsf = $sf ;
         }
      }

      print {$fh->{fam_seq}} $line ;
      print {$fh->{sf_seq}} $line ;
   }
   close($fh->{fam_seq}) ;
   close($fh->{sf_seq}) ;

}


=head2 run_mafft_fillemtpy()

   Title:       run_mafft_fillempty()
   Function:    Fills up empty mafft output alignment files with the original
                  source fasta sequence file
   Args:        none
   Returns:     nothing

   Files in:    foreach SCOP family with empty MAFFT alignment file:
                  <ASTRAL sf aln directory>/<family>.fa (fasta)
                foreach SCOP superfamily:
                  <ASTRAL sf aln directory>/<superfamily>.fa (fasta)

   Files out:   foreach SCOP family with empty MAFFT alignment file:
                  <ASTRAL sf aln directory>/<family>.aln.fa (fasta)
                foreach SCOP superfamily with empty MAFFT alignment file:
                  <ASTRAL fam aln directory>/<wuperfamily>.aln.fa (fasta)

=cut

sub run_mafft_fillempty {
   require DBI ;

   my ($fn,$fh);
   my $pibase_specs = pibase::get_specs() ;
   $fn->{fam_seq} =  $pibase_specs->{astral}->{fam_seq} ;
   $fn->{sf_seq} =  $pibase_specs->{astral}->{sf_seq} ;
   $fn->{gd_seq} =  $pibase_specs->{astral}->{gd_seq} ;

   opendir(FAM, $fn->{fam_seq}) ;
   my @fams ;
   while (my $fam = readdir(FAM)) {
      if ($fam =~ /\.aln\.fa$/) {
         push @fams, $fam;}}
   closedir(FAM) ;

   chdir $fn->{fam_seq} ;
   foreach my $fam (@fams) {
      if (-z $fam) {
         print "FAM aln seqcopy: $fam\n" ;
         my $seqfn = $fam ;$seqfn =~ s/aln.fa$/fa/ ;
         open (FILL, ">$fam") ;
         open (SOURCE, $seqfn) ;
         while (my $line = <SOURCE>) {
            if ($line !~ /^\>/) {$line = uc($line);}
            print FILL $line;
         }
         close SOURCE; close FILL ;
      }
   }


   opendir(SF, $fn->{sf_seq}) ;
   my @sfs;
   while (my $sf = readdir(SF)) { if ($sf =~ /\.aln\.fa$/) {push @sfs, $sf;}}
   closedir(SF) ;

   chdir $fn->{sf_seq} ;
   foreach my $sf (@sfs) {
      if (-z $sf) {
         print "SF aln seqcopy: $sf\n" ;
         my $seqfn = $sf ;$seqfn =~ s/aln.fa$/fa/ ;
         open (FILL, ">$sf") ;
         open (SOURCE, $seqfn) ;
         while (my $line = <SOURCE>) {
            if ($line !~ /^\>/) {$line = uc($line);}
            print FILL $line;
         }
         close SOURCE; close FILL ;
      }
   }

}


=head2 run_mafft_local()

   Title:       run_mafft_local()
   Function:    Runs MAFFT locally to create sf and fam-wide ASTRAL alignments
     useful before ASTRAL releases the alignments for a new SCOP release
     (used for 1.69/modtie paper)
   Args:        none
   Returns:     nothing

   Files in:    foreach SCOP family:
                  <ASTRAL sf aln directory>/<family>.fa (fasta)
                foreach SCOP superfamily:
                  <ASTRAL sf aln directory>/<superfamily>.fa (fasta)

   Files out:   foreach SCOP family with empty MAFFT alignment file:
                  <ASTRAL sf aln directory>/<family>.aln.fa (fasta)
                foreach SCOP superfamily with empty MAFFT alignment file:
                  <ASTRAL fam aln directory>/<wuperfamily>.aln.fa (fasta)

=cut

sub run_mafft_local {

   require DBI ;

   my $mafft_com = "mafft --retree 2 --maxiterate 0 " ;

   my ($fn,$fh);
   my $pibase_specs = pibase::get_specs() ;
   $fn->{fam_seq} =  $pibase_specs->{astral}->{fam_seq} ;
   $fn->{sf_seq} =  $pibase_specs->{astral}->{sf_seq} ;
   $fn->{gd_seq} =  $pibase_specs->{astral}->{gd_seq} ;

   opendir(FAM, $fn->{fam_seq}) ;
   my @fams ;
   while (my $fam = readdir(FAM)) {
      if ($fam =~ /.fa/) {
         push @fams, $fam; }}
   closedir(FAM) ;

   my $count = 1 ;
   chdir $fn->{fam_seq} ;
   foreach my $fam (@fams) {
      print STDERR "FAM $count: $fam\n" ;
      my $out_fn = $fam ; $out_fn =~ s/fa$/aln.fa/ ;
      system("$mafft_com $fam > $out_fn") ;
      $count++ ;
   }


   opendir(SF, $fn->{sf_seq}) ;
   my @sfs;
   while (my $sf = readdir(SF)) { if ($sf =~ /.fa/) {push @sfs, $sf;}}
   closedir(SF) ;

   chdir $fn->{sf_seq} ;
   my $cuont = 1 ;
   foreach my $sf (@sfs) {
      print STDERR "SF $count: $sf\n" ;
      my $out_fn = $sf ; $out_fn =~ s/fa$/aln.fa/ ;
      system("$mafft_com $sf > $out_fn") ;
      $count++ ;
   }

}

1;
