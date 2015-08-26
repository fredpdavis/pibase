=head1 NAME

pibase::ASTRAL - perl module to access ASTRAL data

=head1 DESCRIPTION

Perl module that contains routines for accessing ASTRAL data, for
use in clustering PIBASE interfaces

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

package pibase::ASTRAL ;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/interface_detect_calc bdp_path_2_id _interface_detect_calc__calc_res_pairs/ ;

use pibase ;
use File::Temp qw/tempfile tempdir/ ;
use File::Path qw/mkpath/ ;
use POSIX qw/ceil/ ;

=head2 load_asteroids_aln()

   Title:       load_asteroids_aln()
   Function:    Loads an ASTRAL ASTEROIDS alignment file
   Args:        $_->{aln_fn}	name of ASTERIODS alignment file
                $_->{seq_fn}	name of corresponding ASTEROIDS sequence file
                $_->{allchains}	
                $_->{gdseqh}	data structure holding contents of gdseqh file
                $_->{seqclcont100}
                $_->{seqcl100}
                $_->{doms}	optional hash list of domains to load.
   Returns:     parse_aln_raf() alignment structure

=cut

sub load_asteroids_aln {

   my $in = shift ;

   my $gdseqh = $in->{gdseqh} ;
   my $seqclcont100 = $in->{seqclcont100} ;
   my $seqcl100 = $in->{seqcl100} ;

   my $doms_enriched ; #supplement with the domain ids of seqcl100 representatives if any of the domains are not themselves representatives
   if (exists $in->{doms}) {
      foreach my $dom (sort keys %{$in->{doms}}) {
         $doms_enriched->{$dom}++ ;
         if (!exists $seqcl100->{$dom}) {next;}
         my $seqcl_no = $seqcl100->{$dom} ;
         my $cluster_rep= $seqclcont100->{$seqcl_no}->[0] ;
         $doms_enriched->{$cluster_rep}++ ; #reprsentative
#      print STDERR "dom $dom, enriched with $cluster_rep\n" ;
      }
   }

   my $read_asteroids_aln_options = {
      aln_fn => $in->{aln_fn},
      seq_fn => $in->{seq_fn},
      allchains => $in->{allchains}
   } ;

   if (exists $in->{doms}) {
      $read_asteroids_aln_options->{doms} = $doms_enriched ; }

   my $data = read_asteroids_aln($read_asteroids_aln_options) ;

   my $added ;

   my @origdoms ;
   if (exists $in->{doms}) {
      @origdoms = keys %{$doms_enriched} ;
   } else { #used to always expect $in->{doms} - changed for pilig.pm code
      @origdoms = keys %{$data->{aln}} ;
   }

   foreach my $rep_scopid (@origdoms) {
      if (!exists $seqcl100->{$rep_scopid}) {next;}
      my $clno = $seqcl100->{$rep_scopid} ;
      if ($#{$seqclcont100->{$clno}} > 0) {
         foreach my $j ( 1 .. $#{$seqclcont100->{$clno}}) {
            my $scopid = $seqclcont100->{$clno}->[$j] ;
#            print STDERR "scopid is $scopid\n" ;
#            print STDERR "now on scopid $scopid (from $rep_scopid) ".__LINE__."\n ";
#            if ($rep_scopid eq 'd1j4za1') {print STDERR " post-aln supplement for d1j4za1 - copying over to $scopid\n";}

            if (exists $added->{$scopid}) {next;}

            $data->{seq}->{$scopid} = $data->{seq}->{$rep_scopid} ;
            $data->{aln}->{$scopid} = $data->{aln}->{$rep_scopid} ;
#            print STDERR "100SUPPed $scopid ($rep_scopid): $data->{aln}->{$scopid}\n" ;

            $data->{defstring}->{$scopid} = $gdseqh->{defstring}->{$scopid} ;
            $data->{class}->{$scopid} = $gdseqh->{class}->{$scopid} ;
            $data->{pdb}->{$scopid} = $gdseqh->{pdb}->{$scopid} ;
            $data->{frags}->{$scopid} = $gdseqh->{frags}->{$scopid} ;
         }
      }
   }

   my $aln = parse_aln_raf({
      alndata => $data,
      raf => $in->{raf}
   }) ;
   $aln->{alnlength} = $data->{alnlength} ;
   $aln->{aln} = $data->{aln} ; #added 090303_2024 for aln strings to return

# added this081229_0008 for entropy score to pass through
   $aln->{meta} = $data->{meta} ;

   return $aln ;
}


sub LASTPRODUCTIONVER_20080619_load_asteroids_aln {

   my $in = shift ;

   my $gdseqh = $in->{gdseqh} ;
   my $seqclcont100 = $in->{seqclcont100} ;
   my $seqcl100 = $in->{seqcl100} ;

   my $doms_enriched ; #supplement with the domain ids of seqcl100 representatives if any of the domains are not themselves representatives
   foreach my $dom (sort keys %{$in->{doms}}) {
      $doms_enriched->{$dom}++ ;
      if (!exists $seqcl100->{$dom}) {next;}
      my $seqcl_no = $seqcl100->{$dom} ;
      my $cluster_rep= $seqclcont100->{$seqcl_no}->[0] ;
      $doms_enriched->{$cluster_rep}++ ; #reprsentative
#      print STDERR "dom $dom, enriched with $cluster_rep\n" ;
   }

# changed doms arg to a list enriched with representatives
#      doms => $in->{doms},

   my $data = read_asteroids_aln({
      doms => $doms_enriched,
      aln_fn => $in->{aln_fn},
      seq_fn => $in->{seq_fn},
      allchains => $in->{allchains}
   });

   my $added ;
# changed from origdoms = data->{aln} to orgidoms = $in->{doms_enriched} ;
# was choking here, trying to load domain data that wasnt requested
#   my @origdoms = keys %{$data->{aln}} ;
   my @origdoms = keys %{$doms_enriched} ;
   foreach my $rep_scopid (@origdoms) {
      if (!exists $seqcl100->{$rep_scopid}) {next;}
      my $clno = $seqcl100->{$rep_scopid} ;
      if ($#{$seqclcont100->{$clno}} > 0) {
         foreach my $j ( 1 .. $#{$seqclcont100->{$clno}}) {
            my $scopid = $seqclcont100->{$clno}->[$j] ;
#            print STDERR "scopid is $scopid\n" ;
#            print STDERR "now on scopid $scopid (from $rep_scopid) ".__LINE__."\n ";
#            if ($rep_scopid eq 'd1j4za1') {print STDERR " post-aln supplement for d1j4za1 - copying over to $scopid\n";}

            if (exists $added->{$scopid}) {next;}

            $data->{seq}->{$scopid} = $data->{seq}->{$rep_scopid} ;
            $data->{aln}->{$scopid} = $data->{aln}->{$rep_scopid} ;
#            print STDERR "100SUPPed $scopid ($rep_scopid): $data->{aln}->{$scopid}\n" ;

            $data->{defstring}->{$scopid} = $gdseqh->{defstring}->{$scopid} ;
            $data->{class}->{$scopid} = $gdseqh->{class}->{$scopid} ;
            $data->{pdb}->{$scopid} = $gdseqh->{pdb}->{$scopid} ;
            $data->{frags}->{$scopid} = $gdseqh->{frags}->{$scopid} ;
         }
      }
   }

   my $aln = parse_aln_raf({
      alndata => $data,
      raf => $in->{raf}
   }) ;

   return $aln ;
}


sub LASTPRODUCTIONVER_20080619_read_asteroids_aln {

   my $in = shift;
   my $fn;
   $fn->{aln} = $in->{aln_fn} ;
   $fn->{seq} = $in->{seq_fn} ;
   my $doms = $in->{doms} ;
   my $allchains = $in->{allchains} ;

   my $data;
   my $cur_dom = '' ;
   my $seq ;
   open(SEQF, $fn->{seq}) ;
   while (my $line = <SEQF>) {
      chomp $line;
      if ($line =~ /^>/) {
         $line =~ s/^>// ;
         my ($t_dom, undef, undef, undef) = split(/\s/, $line) ;
         if (exists $doms->{$t_dom}) {
            $data->{seq}->{$t_dom} = ''; }
         $cur_dom = $t_dom ;
      } else {
         if (exists $doms->{$cur_dom}) {
            my $newseq = $line ;
            # DONT CHANGE CASE OF x; x = unknown res, X = fragment break.
            $newseq =~ tr/[a-w]/[A-W]/ ;
            $newseq =~ tr/[y-z]/[Y-Z]/ ;
            $data->{seq}->{$cur_dom} .= $newseq ;
         } #uc to handle 1-seq aln's
      }
   }
   close(SEQF) ;

   $cur_dom = '' ;
   open(ALNF, $fn->{aln}) ;
#   print STDERR "parsing $fn->{aln}\n" ;
   while (my $line = <ALNF>) {
      chomp $line;
      if ($line =~ /^>/) {
         $line =~ s/^>// ;
         my ($t_dom, $t_class, $t_def, undef) = split(/\s/, $line) ;
         $cur_dom = $t_dom ;
         if (!exists $in->{doms}->{$cur_dom}) {next;}

         $t_def =~ s/^\(// ; $t_def =~ s/\)$// ;
         $data->{defstring}->{$t_dom} = $t_def ;
         $data->{class}->{$t_dom} = $t_class;
         $data->{aln}->{$t_dom} = '';
         $data->{pdb}->{$t_dom} = substr($t_dom, 1, 4) ;
         {
            my @t = split(/\,/, $t_def) ;
            my (@ch, @b, @e) ;
            foreach my $j ( 0 .. $#t) {
               my $t_frag ;
               $t_frag->{chain} = '_' ;
               if ($t[$j] =~ /:/) {
                  $t_frag->{chain} = substr($t[$j],0,1) ; }
               $t[$j] =~ s/^.\:// ;
               my ($b, $e) = (' ', ' ');
               if ($t[$j] =~ /.+\-.+/) {
                  ($t_frag->{b}, $t_frag->{e}) =
                  ($t[$j] =~ /(.+)\-(.+)/) ; }
               push @{$data->{frags}->{$t_dom}}, $t_frag ;

               if ($t_frag->{chain} eq '-' || $t_frag->{chain} eq '_') {
                  $t_frag->{chain} = ' ' ;
               } elsif (!exists $allchains->{$data->{pdb}->{$t_dom}}->{$t_frag->{chain}}) {
                  if ($t_frag->{chain} eq uc($t_frag->{chain})) {
                     my $lc = lc($t_frag->{chain}) ;
                     if (exists  $allchains->{$data->{pdb}->{$t_dom}}->{$lc}) {
                        print STDERR "WARNING defstring change: dom $t_dom old chain $t_frag->{chain} changed to $lc\n" ;
                        $t_frag->{chain} = $lc ; }
                  } elsif ($t_frag->{chain} eq lc($t_frag->{chain})) {
                     my $uc = uc($t_frag->{chain}) ;
                     if (exists $allchains->{$data->{pdb}->{$t_dom}}->{$uc}) {
                        print STDERR "WARNING defstring change: dom $t_dom old chain $t_frag->{chain} changed to $uc\n" ;
                        $t_frag->{chain} = $uc ; }
                  }
               }
            }
         }
      } elsif (exists $in->{doms}->{$cur_dom}) {
         my $newseq = $line ;
         $newseq =~ tr/[a-w]/[A-W]/ ; #for 1-seq alns eg a.205.1.1 d1us7b_
         $newseq =~ tr/[y-z]/[Y-Z]/ ;
         $data->{aln}->{$cur_dom} .= $newseq ;
      }
   }
   close(ALNF) ;

   return $data ;

}


=head2 raf_preload()

   Title:       raf_preload()
   Function:    Loads the ASTRAL raf file
   Args:        $_->{fn}	name of the ASTRAL raf file
   Returns:     ->{pdb_chain} = "RAF line contents"

=cut

sub raf_preload {
   my $in = shift ;

   my $raf ;
   open (RAF, $in->{fn}) ;
   while (my $line = <RAF>) {
      chomp $line;

      if ($line =~ /^# Header:/) {
         ($raf->{header_length}) = ($line =~ /: ([0-9]+) Bytes/) ;
      } elsif ($line !~ /^#/) {
         my $pdbch = substr($line, 0, 5) ;
         $raf->{$pdbch} = $line ;
      }
   }
   close(RAF) ;
   return $raf ;
}


=head2 parse_raf_line()

   Title:       parse_raf_line()
   Function:    Parses a line from ASTRAL raf file
   Args:        $_->{line}	RAF line
                $_->{headlength}	length of RAF file header

   Returns:     ->{atomresno_first}	ATOM residue number of first residue
                ->{atomresno_last}	ATOM residue number of last residue
                ->{atomresna}	= [ ATOM residue name 1,2, ... ]
                ->{atomresno}	= [ ATOM residue number 1,2, ... ]
                ->{seqresna}	= [ sequence residue name 1,2, ... ]
                ->{seqresno}	= [ sequence residue number 1,2, ... ]
                ->{seqresno2ind}->{seqresno} = index in @{->{seqresno}}
                ->{ind2seqresno}->{index in @{->{seqresno}}} = seqresno
                ->{atom2seqresno_back}->{$atomresno} = seqresno
                ->{atomresno2ind_back}->{$atomresno} = seqresno index
                ->{seq2atomresno}->{$seqresno} = $atomresno ;
                ->{atom2seqresno}->{$atomresno} = $seqresno ;
                ->{atomresno2ind}->{$atomresno} = $#{$res->{seqresno}} ;
                ->{ind2atomresno}->{$#{$res->{seqresno}}} = $atomresno  ;

=cut

sub parse_raf_line {

   my $in = shift ;
   my $line = $in->{line}  ;
   my $headlength = $in->{headlength} ;

   my $res ;

   my $head = substr($line, 0, $headlength) ;
   $res->{atomresno_first} = substr($head, -10, 5) ;
   $res->{atomresno_first} =~ s/ //g ;
   $res->{atomresno_last} = substr($head, -5, 5) ;
   $res->{atomresno_last} =~ s/ //g ;
   my $pos = $headlength ;
   my $seqresno = 0 ;
   my $firstseqresno ;

   my @missing_atomno = ();

   while ($pos < length($line)) {
      my $t_res = substr($line,$pos,7) ;
      my $atomresno = substr($t_res,0,4) ; $atomresno =~ s/ //g ;
      my $inscode = substr($t_res,4,1) ;
      my $atomresna = substr($t_res,5,1) ;
      my $seqresna = substr($t_res,6,1) ;

      if ($inscode ne ' ') {$atomresno .= $inscode ; }
      push @{$res->{atomresna}}, $atomresna ;
      push @{$res->{seqresna}}, $seqresna ;
      push @{$res->{atomresno}}, $atomresno;

      if ($seqresna ne '.') {
         $seqresno++ ;
         push @{$res->{seqresno}}, $seqresno;
         $res->{seqresno2ind}->{$seqresno}= $#{$res->{seqresno}};
         $res->{ind2seqresno}->{$#{$res->{seqresno}}} = $seqresno  ;
      } else {
         if ($#{$res->{seqresno}} >= 0 ) {
            $res->{atom2seqresno_back}->{$atomresno} = ${$res->{seqresno}}[-1] ;
            $res->{atomresno2ind_back}->{$atomresno} = $#{$res->{seqresno}} ;
         }

         push @{$res->{seqresno}}, '.' ;
         push @missing_atomno, $atomresno ;
      }

      if ($seqresna ne '.' && ($#missing_atomno >= 0)) {
         foreach my $t_atomresno (@missing_atomno) {
            $res->{atom2seqresno_ahead}->{$t_atomresno} = $seqresno ;
            $res->{atomresno2ind_ahead}->{$t_atomresno} = $#{$res->{seqresno}} ;
         }
         @missing_atomno = () ;
      }

      if ($seqresna ne '.' && $atomresna ne '.') {
         $res->{seq2atomresno}->{$seqresno} = $atomresno ;
         $res->{atom2seqresno}->{$atomresno} = $seqresno ;
         $res->{atomresno2ind}->{$atomresno} = $#{$res->{seqresno}} ;
         $res->{ind2atomresno}->{$#{$res->{seqresno}}} = $atomresno  ;
      }
      $pos += 7 ;
   }

# debug check080616_1102 for debugging 1wf4j1 problems
#   foreach my $t_res ( sort {$a <=> $b} keys %{$res->{atomresno2ind}}) {
#      print STDERR "raf $t_res   = ".$res->{atomresno2ind}->{$t_res}."\n" ;
#   }

   return $res ;
}


=head2 read_asteroids_aln()

   Title:       read_asteroids_aln()
   Function:    Reads an ASTEROIDS alignment file
   Args:        $_->{aln_fn} ASTEROIDS alignment file name
                $_->{seq_fn} corresponding ASTEROIDS sequence file name
                $_->{allchains}->{domain} = pdb chain

   Returns:     ->{seq}->{domain} = 'DOMAINSEQVENCE';
                ->{defstring}->{domain} = definition line from alignment
                ->{class}->{domain} = SCOP class
                ->{aln}->{domain} = domain sequence froma alignment
                ->{pdb}->{domain} = PDB code for the domain
                ->{frags}->{domain} = [{b => startresidue, e => endresidue},...]
                ->{alnlength} = alignment length

=cut

sub read_asteroids_aln {

   my $in = shift;
   my $fn;
   $fn->{aln} = $in->{aln_fn} ;
   $fn->{seq} = $in->{seq_fn} ;
   my $allchains = $in->{allchains} ;


   my $standard_20aa = {
      'A' => 1,
      'R' => 1,
      'N' => 1,
      'D' => 1,
      'C' => 1,
      'Q' => 1,
      'E' => 1,
      'G' => 1,
      'H' => 1,
      'I' => 1,
      'L' => 1,
      'K' => 1,
      'M' => 1,
      'F' => 1,
      'P' => 1,
      'S' => 1,
      'T' => 1,
      'W' => 1,
      'Y' => 1,
      'V'=> 1,
   } ;

   my $data;
   my $cur_dom = '' ;
   my $seq ;
   open(SEQF, $fn->{seq}) ;
   while (my $line = <SEQF>) {
      chomp $line;
      if ($line =~ /^>/) {
         $line =~ s/^>// ;
         my ($t_dom, undef, undef, undef) = split(/\s/, $line) ;
         if (!exists $in->{doms} ||
             (exists $in->{doms} && exists $in->{doms}->{$t_dom})) {
            $data->{seq}->{$t_dom} = ''; }
         $cur_dom = $t_dom ;
      } elsif (!exists $in->{doms} ||
               (exists $in->{doms} && exists $in->{doms}->{$cur_dom})) {
            my $newseq = $line ;
            # DONT CHANGE CASE OF x; x = unknown res, X = fragment break.
            $newseq =~ tr/[a-w]/[A-W]/ ;
            $newseq =~ tr/[y-z]/[Y-Z]/ ;
            $data->{seq}->{$cur_dom} .= $newseq ;
 #uc to handle 1-seq aln's
      }
   }
   close(SEQF) ;

   $cur_dom = '' ;
   open(ALNF, $fn->{aln}) ;
#   print STDERR "parsing $fn->{aln}\n" ;
   while (my $line = <ALNF>) {
      chomp $line;
      if ($line =~ /^>/) {
         $line =~ s/^>// ;
         my ($t_dom, $t_class, $t_def, undef) = split(/\s/, $line) ;
         $cur_dom = $t_dom ;
         if (exists $in->{doms} && !exists $in->{doms}->{$cur_dom}) {
#            print STDERR "NEXTED: $cur_dom\n";
            next;
         }

         $t_def =~ s/^\(// ; $t_def =~ s/\)$// ;
         $data->{defstring}->{$t_dom} = $t_def ;
         $data->{class}->{$t_dom} = $t_class;
         $data->{aln}->{$t_dom} = '';
         $data->{pdb}->{$t_dom} = substr($t_dom, 1, 4) ;
         {
            my @t = split(/\,/, $t_def) ;
            my (@ch, @b, @e) ;
            foreach my $j ( 0 .. $#t) {
               my $t_frag ;
               $t_frag->{chain} = '_' ;
               if ($t[$j] =~ /:/) {
                  $t_frag->{chain} = substr($t[$j],0,1) ; }
               $t[$j] =~ s/^.\:// ;
               my ($b, $e) = (' ', ' ');
               if ($t[$j] =~ /.+\-.+/) {
                  ($t_frag->{b}, $t_frag->{e}) =
                  ($t[$j] =~ /(.+)\-(.+)/) ; }
               push @{$data->{frags}->{$t_dom}}, $t_frag ;

               if ($t_frag->{chain} eq '-' || $t_frag->{chain} eq '_') {
                  $t_frag->{chain} = ' ' ;
               } elsif (!exists $allchains->{$data->{pdb}->{$t_dom}}->{$t_frag->{chain}}) {
                  if ($t_frag->{chain} eq uc($t_frag->{chain})) {
                     my $lc = lc($t_frag->{chain}) ;
                     if (exists  $allchains->{$data->{pdb}->{$t_dom}}->{$lc}) {
                        print STDERR "WARNING defstring change: dom $t_dom old chain $t_frag->{chain} changed to $lc\n" ;
                        $t_frag->{chain} = $lc ; }
                  } elsif ($t_frag->{chain} eq lc($t_frag->{chain})) {
                     my $uc = uc($t_frag->{chain}) ;
                     if (exists $allchains->{$data->{pdb}->{$t_dom}}->{$uc}) {
                        print STDERR "WARNING defstring change: dom $t_dom old chain $t_frag->{chain} changed to $uc\n" ;
                        $t_frag->{chain} = $uc ; }
                  }
               }
            }
         }
      } elsif ((!exists $in->{doms}) ||
                (exists $in->{doms} && exists $in->{doms}->{$cur_dom})) {
         my $newseq = $line ;
         $newseq =~ tr/[a-w]/[A-W]/ ; #for 1-seq alns eg a.205.1.1 d1us7b_
         $newseq =~ tr/[y-z]/[Y-Z]/ ;
         $data->{aln}->{$cur_dom} .= $newseq ;
      }
#      print STDERR "CUR_DOM: $cur_dom; aln_string = ".$data->{aln}->{$cur_dom}."; line = $line\n" ;
   }
   {
      my $numdom = keys %{$data->{aln}} ;
      if ($numdom  == 0) {return $data;} #Eg, d1ncpc_ (g.40.1.1) - no data available in ASTRAL

      my ($tdom) = keys %{$data->{aln}};
      $data->{alnlength} = length($data->{aln}->{$tdom}) ;
   }
#   $data->{alnlength} = length($data->{aln}->{$cur_dom}) ;
   close(ALNF) ;

# fpd081229_0003 - added for pilig
   foreach my $pos ( 0 .. ($data->{alnlength} - 1) ) {
      my $freq ;
      my $freq_std20 ;
      my $total_std20 = 0;

      my $numdoms = keys %{$data->{aln}} ;
#      my @tp ;
      foreach my $dom (keys %{$data->{aln}}) {
         my $char = substr($data->{aln}->{$dom}, $pos, 1) ;
         if (exists $standard_20aa->{$char}) {
            $freq_std20->{$char}++ ;
            $total_std20++ ; }
#         push @tp, substr($data->{aln}->{$dom}, $pos, 1) ;
         $freq->{substr($data->{aln}->{$dom}, $pos, 1)}++ ;
      }

      %{$data->{meta}->{res_freq}->{$pos}} = %{$freq} ;
#      %{$data->{meta}->{res_freq_std20}->{$pos}} = %{$freq_std20} ;

      $data->{meta}->{numrestypes}->{$pos} = keys %{$freq} ;
      $data->{meta}->{numrestypes_std20}->{$pos} = keys %{$freq_std20} ;

      $data->{meta}->{shannone}->{$pos} = 0 ;
      $data->{meta}->{shannone_std20}->{$pos} = 0 ;
      foreach my $res (keys %{$freq}) {
         $data->{meta}->{shannone}->{$pos} +=
            ($freq->{$res} / $numdoms) *
            (log($freq->{$res} / $numdoms) / log(21)) ;

         if (exists $standard_20aa->{$res}) {
            $data->{meta}->{shannone_std20}->{$pos} +=
               ($freq->{$res} / $total_std20) *
               (log($freq->{$res} / $total_std20) / log(20)) ;
         }
      }
      $data->{meta}->{shannone}->{$pos} = -1 *
         $data->{meta}->{shannone}->{$pos}  ;

      $data->{meta}->{shannone_std20}->{$pos} = -1 *
         $data->{meta}->{shannone_std20}->{$pos}  ;

#      if ($data->{shannone}->{$pos} == 0) {
#         print STDERR " shannone pos $pos $in->{aln_fn} = 0 charc:".join(",", @tp)."\n" ; }

   }


   return $data ;

}

sub OLDpre081229_0000_read_asteroids_aln_WITHOUTENTROPY {

   my $in = shift;
   my $fn;
   $fn->{aln} = $in->{aln_fn} ;
   $fn->{seq} = $in->{seq_fn} ;
   my $allchains = $in->{allchains} ;

   my $data;
   my $cur_dom = '' ;
   my $seq ;
   open(SEQF, $fn->{seq}) ;
   while (my $line = <SEQF>) {
      chomp $line;
      if ($line =~ /^>/) {
         $line =~ s/^>// ;
         my ($t_dom, undef, undef, undef) = split(/\s/, $line) ;
         if (!exists $in->{doms} ||
             (exists $in->{doms} && exists $in->{doms}->{$t_dom})) {
            $data->{seq}->{$t_dom} = ''; }
         $cur_dom = $t_dom ;
      } elsif (!exists $in->{doms} ||
               (exists $in->{doms} && exists $in->{doms}->{$cur_dom})) {
            my $newseq = $line ;
            # DONT CHANGE CASE OF x; x = unknown res, X = fragment break.
            $newseq =~ tr/[a-w]/[A-W]/ ;
            $newseq =~ tr/[y-z]/[Y-Z]/ ;
            $data->{seq}->{$cur_dom} .= $newseq ;
 #uc to handle 1-seq aln's
      }
   }
   close(SEQF) ;

   $cur_dom = '' ;
   open(ALNF, $fn->{aln}) ;
#   print STDERR "parsing $fn->{aln}\n" ;
   while (my $line = <ALNF>) {
      chomp $line;
      if ($line =~ /^>/) {
         $line =~ s/^>// ;
         my ($t_dom, $t_class, $t_def, undef) = split(/\s/, $line) ;
         $cur_dom = $t_dom ;
         if (exists $in->{doms} && !exists $in->{doms}->{$cur_dom}) {next;}

         $t_def =~ s/^\(// ; $t_def =~ s/\)$// ;
         $data->{defstring}->{$t_dom} = $t_def ;
         $data->{class}->{$t_dom} = $t_class;
         $data->{aln}->{$t_dom} = '';
         $data->{pdb}->{$t_dom} = substr($t_dom, 1, 4) ;
         {
            my @t = split(/\,/, $t_def) ;
            my (@ch, @b, @e) ;
            foreach my $j ( 0 .. $#t) {
               my $t_frag ;
               $t_frag->{chain} = '_' ;
               if ($t[$j] =~ /:/) {
                  $t_frag->{chain} = substr($t[$j],0,1) ; }
               $t[$j] =~ s/^.\:// ;
               my ($b, $e) = (' ', ' ');
               if ($t[$j] =~ /.+\-.+/) {
                  ($t_frag->{b}, $t_frag->{e}) =
                  ($t[$j] =~ /(.+)\-(.+)/) ; }
               push @{$data->{frags}->{$t_dom}}, $t_frag ;

               if ($t_frag->{chain} eq '-' || $t_frag->{chain} eq '_') {
                  $t_frag->{chain} = ' ' ;
               } elsif (!exists $allchains->{$data->{pdb}->{$t_dom}}->{$t_frag->{chain}}) {
                  if ($t_frag->{chain} eq uc($t_frag->{chain})) {
                     my $lc = lc($t_frag->{chain}) ;
                     if (exists  $allchains->{$data->{pdb}->{$t_dom}}->{$lc}) {
                        print STDERR "WARNING defstring change: dom $t_dom old chain $t_frag->{chain} changed to $lc\n" ;
                        $t_frag->{chain} = $lc ; }
                  } elsif ($t_frag->{chain} eq lc($t_frag->{chain})) {
                     my $uc = uc($t_frag->{chain}) ;
                     if (exists $allchains->{$data->{pdb}->{$t_dom}}->{$uc}) {
                        print STDERR "WARNING defstring change: dom $t_dom old chain $t_frag->{chain} changed to $uc\n" ;
                        $t_frag->{chain} = $uc ; }
                  }
               }
            }
         }
      } elsif ((!exists $in->{doms}) ||
                (exists $in->{doms} && exists $in->{doms}->{$cur_dom})) {
         my $newseq = $line ;
         $newseq =~ tr/[a-w]/[A-W]/ ; #for 1-seq alns eg a.205.1.1 d1us7b_
         $newseq =~ tr/[y-z]/[Y-Z]/ ;
         $data->{aln}->{$cur_dom} .= $newseq ;
      }
   }
   $data->{alnlength} = length($data->{aln}->{$cur_dom}) ;
   close(ALNF) ;

   return $data ;

}

=head2 parse_aln_raf()

   Title:       parse_aln_raf()
   Function:    Reads an ASTEROIDS alignment file
   Args:        $_->{alndata}  - alignment data from read_asteroids_aln()
                $_->{raf} - RAF data from raf_preload()

   Returns:     ->{pos2resno}->{domain}->{alignment position} = ATOM resno
                ->{resno2pos}->{domain}->{ATOM resno} = alignment position

=cut

sub parse_aln_raf {

   my $in = shift;
   my $data = $in->{alndata} ;
   my $raf = $in->{raf} ;

#data->{seq} has raw FASTA sequence strings - check substr(curpos) to see
# if x is capital = frag or lower case = residue unk

   my ($alnpos2atomresno, $atomresno2alnpos) ;
   foreach my $dom (keys %{$data->{defstring}}) {
      if (!exists $data->{aln}->{$dom}) {
         print STDERR "NOTE: skipping $dom in parse_aln_raf(): no aln info\n" ;
         next;
      }

#      print STDERR "now on dom $dom\n" ;
#      if ($dom ne 'd1pcqf1') {next;}
#      if (DEBUGALN) {
#         print STDERR "now on $dom\n" ;
#         print STDERR "   defstring:\t$data->{defstring}->{$dom}\n" ;
#         print STDERR "   alnstring:\t$data->{aln}->{$dom}\n" ;
#         print STDERR "   pdb:\t$data->{pdb}->{$dom}\n" ;
#         print STDERR "   class:\t$data->{class}->{$dom}\n" ;
#      }

      my $gs2atomresno ;
      my $gs2resna ; #to verify correct positioning when parsing alignment

#get gsresno -> atomresno mapping from raf file
      my $gsresno = -1 ;
      my $abort_thisdom = 0 ;
      foreach my $fragdef ( @{$data->{frags}->{$dom}}) {
         $gsresno++ ; #to take care of X in between fragments
         my $tch = $fragdef->{chain} ; if ($tch eq ' ') {$tch = '_';}
         if (!exists $raf->{$data->{pdb}->{$dom}.$tch} ) {
            my $t_error = "WARNING: $dom chain (".$data->{pdb}->{$dom}." $tch) not found in raf file" ;
            print STDERR $t_error ;
            my $tch_lc = lc($tch) ;
            if (exists $raf->{$data->{pdb}->{$dom}.$tch_lc}) {
               print STDERR "; using chain $tch_lc instead\n" ;
               $tch = $tch_lc ;
            } else {
               print STDERR "; skipping the domain\n" ;
               $abort_thisdom = 1 ;
               last;
            }
         }
         my $line = $raf->{$data->{pdb}->{$dom}.$tch} ;
         my $rafinfo = parse_raf_line({
            line => $line,
            headlength => $raf->{header_length}
         }) ;

         my ($ind_b, $ind_e) ;
         if (exists $fragdef->{b}) {
            if (defined $rafinfo->{atomresno2ind}->{$fragdef->{b}}) {
               $ind_b = $rafinfo->{atomresno2ind}->{$fragdef->{b}} ;
            } else {
#               if (DEBUGALN) {print STDERR " on dom $dom: lookahead b find for $fragdef->{b}\n" ;}
               $ind_b = $rafinfo->{atomresno2ind_ahead}->{$fragdef->{b}} ;
            }

            if (defined $rafinfo->{atomresno2ind}->{$fragdef->{e}}) {
               $ind_e = $rafinfo->{atomresno2ind}->{$fragdef->{e}} ;
            } else {
#               if (DEBUGALN) {print STDERR " on dom $dom: lookback e find for $fragdef->{e}\n" ;}
               $ind_e = $rafinfo->{atomresno2ind_back}->{$fragdef->{e}} ;
            }

         } else {


            if (defined $rafinfo->{atomresno2ind}->{$rafinfo->{atomresno_first}}) {
               $ind_b= $rafinfo->{atomresno2ind}->{$rafinfo->{atomresno_first}};
            } else {
#               if (DEBUGALN) {print STDERR " on dom $dom: lookahead b find for $rafinfo->{atomresno_first}\n" ;}
               $ind_b= $rafinfo->{atomresno2ind_ahead}->{$rafinfo->{atomresno_first}} ;
            }

            if (defined $rafinfo->{atomresno2ind}->{$rafinfo->{atomresno_last}}) {
               $ind_e = $rafinfo->{atomresno2ind}->{$rafinfo->{atomresno_last}};
            } else {
#               if (DEBUGALN) {print STDERR " on dom $dom: lookback e find for $rafinfo->{atomresno_last}\n" ;}
               $ind_e = $rafinfo->{atomresno2ind_back}->{$rafinfo->{atomresno_last}} ;
            }
         }

         my $keepon = 1 ;
         my $pos = $ind_b ;
         if (!defined $pos) { print STDERR " on dom $dom (chain $fragdef->{chain}), pos is unefined\n" ;}
         if (!defined $ind_e) { print STDERR " on dom $dom (chain $fragdef->{chain}), ind_e is udnefinedpos\n";}

         while ($pos <= $ind_e) {
            if (exists $rafinfo->{ind2seqresno}->{$pos}) {
               if (exists $rafinfo->{ind2atomresno}->{$pos}) {
                  $gs2resna->{$gsresno} = uc($rafinfo->{seqresna}->[$pos]) ;
                  $gs2atomresno->{$gsresno} = $rafinfo->{ind2atomresno}->{$pos}."\n".$fragdef->{chain};
               }
               $gsresno++ ;
            }

            $pos++ ;
         }
      }
      if ($abort_thisdom) {next;}

#go through aln entry and assign aln pos -> atomresno mapping

#      print STDERR "\n\naln: ".$data->{aln}->{$dom}."\n" ;
#      print STDERR "seq: ".$data->{seq}->{$dom}."\n" ;
#      if (!exists $data->{seq}->{$dom} ||
#          !defined $data->{seq}->{$dom}) {
#         die "oh shit, seq doesn't exist for $dom\n" ;
#      }

      my $curseqpos = 0 ;
      foreach my $pos (0 .. (length($data->{aln}->{$dom}) - 1)) {

#         print STDERR "$dom, seqpos $curseqpos (of ".length($data->{seq}->{$dom})."), alnpos $pos (of ".length($data->{aln}->{$dom})."); " ;
#         print STDERR " seqchar: ".substr($data->{seq}->{$dom},$curseqpos,1).", alnchar: ".substr($data->{aln}->{$dom},$pos,1)."\n ";
#
#         if ($curseqpos >= length($data->{seq}->{$dom})) {
#            print STDERR "WE HAVE PROBLEMS: curseqpos is >= seq len: ".
#               length($data->{seq}->{$dom})."\n";
#         }

         my $alnchar = substr($data->{aln}->{$dom}, $pos, 1) ;

#         if (!defined $alnchar || !defined (substr($data->{seq}->{$dom}, $curseqpos, 1))) {
#            die "oh shit, undefined alnchar $alnchar of seqchar: ".
#substr($data->{seq}->{$dom}, $curseqpos, 1)."\n" ;
#         }

#         if (DEBUGALN && $dom eq DEBUGDOM) {
#            print STDERR "dom $dom (alnpos $pos, curseqpos $curseqpos) $alnchar\n" ;
#         }

         if ($alnchar eq 'X' &&
            (substr($data->{seq}->{$dom}, ($curseqpos), 1) eq 'X')) {

#            if (DEBUGALN) {
#               print STDERR "$dom: real frag $alnchar at seqpos $curseqpos (alnpos $pos)\n" ;}

            $curseqpos++ ;

         } elsif ($alnchar ne '-') {
#            if (DEBUGALN) {
#               if ($alnchar eq 'X') {print STDERR "$dom caught unk res: alnchar $alnchar at seq pos $curseqpos (alnpos $pos)\n";}
#            }

            if (exists $gs2atomresno->{$curseqpos}) {

               if ($alnchar ne $gs2resna->{$curseqpos}) {
                  print "ERROR $dom: alignment position $pos ($alnchar) mismatched with gs sequence position $curseqpos ($gs2resna->{$curseqpos}\n" ;
                  print STDERR "ERROR $dom: alignment position $pos ($alnchar) mismatched with gs sequence position $curseqpos ($gs2resna->{$curseqpos}\n" ;
               }
#               if (DEBUGALN) {
#                  if ($dom eq DEBUGDOM) {
#                     print STDERR "$dom\t$pos\t$alnchar\t$gs2atomresno->{$curseqpos}\n" ;}}

               $alnpos2atomresno->{$dom}->{$pos} = $gs2atomresno->{$curseqpos} ;
               $atomresno2alnpos->{$dom}->{$gs2atomresno->{$curseqpos}} = $pos ;
            }
            $curseqpos++ ;
         }
      }
   }


   return {
     pos2resno => $alnpos2atomresno,
     resno2pos => $atomresno2alnpos
   } ;

}


=head2 load_astral_headers()

   Title:       load_astral_headers()
   Function:    Loads ASTRAL headers
   Args:        $_->{fn} - alignment file name
                $_->{raf} - RAF data from raf_preload()

   Returns:     ->{gdseqh}->{defstring}->{domain} = domain definition string
                ->{gdseqh}->{class}->{domain} = domain class
                ->{gdseqh}->{pdb}->{domain} = domain PDB code
                ->{gdseqh}->{frags}->{domain} = [{b => startres, e => endres}...]
                ->{gdom}->{domain (w d prefix)} = domain (w g prefix)

=cut

sub load_astral_headers {

   my $in = shift;
   my $gdseq_fn = $in->{fn} ;
   open(GDSEQF, $gdseq_fn) ;
   my $gdseqh ;
   my $gdom_fl ;
   while (my $line = <GDSEQF>) {
      if ($line =~ /^>/) {
         $line =~ s/^>// ;
         my ($t_dom, $t_class, $t_def, undef) = split(/\s/, $line) ;
         $t_def =~ s/^\(// ; $t_def =~ s/\)$// ;
         $gdseqh->{defstring}->{$t_dom} = $t_def ;
         $gdseqh->{class}->{$t_dom} = $t_class;
         $gdseqh->{pdb}->{$t_dom} = substr($t_dom, 1, 4) ;
         if ($t_dom =~ /^g/) {
            my $t_ddom = $t_dom ; $t_ddom =~ s/^g/d/ ;
            $gdom_fl->{$t_ddom} = $t_dom ;
         }

         {
            my @t = split(/\,/, $t_def) ;
            my (@ch, @b, @e) ;
            foreach my $j ( 0 .. $#t) {
               my $t_frag ;
               $t_frag->{chain} = '_' ;
               if ($t[$j] =~ /:/) {
                  $t_frag->{chain} = substr($t[$j],0,1) ; }
               $t[$j] =~ s/^.\:// ;
               my ($b, $e) = (' ', ' ');
               if ($t[$j] =~ /.+\-.+/) {
                  ($t_frag->{b}, $t_frag->{e}) =
                  ($t[$j] =~ /(.+)\-(.+)/) ; }
               push @{$gdseqh->{frags}->{$t_dom}}, $t_frag ;
            }
         }
      }
   }
   close(GDSEQF) ;

   return {
      gdseqh => $gdseqh,
      gdom => $gdom_fl
   } ;

}


=head2 load_astral_clusters()

   Title:       load_astral_clusters()
   Function:    Loads ASTRAL sequence cluster definitions
   Args:        $_->{out} - pointer to hash to hold output
                $_->{pibase_specs} -  pibase_specs structure

   Returns:     Nothing - populates the specified $_->{out}
                {out}->{seqcl}->{seq identity}->{scop identifier} = cluster num
                {out}->{seqcl2cont}->{seq identity}->{cluster num} = [scop id,...]

=cut

sub load_astral_clusters {

   my $in = shift ;
   my $out = $in->{out} ;
   my $pibase_specs = $in->{pibase_specs} ;

   my $seqcl ;
   my $seqcl2cont ;
   {
      foreach my $seqid (keys %{$pibase_specs->{astral}->{seqcl}}) {
         open (ASTRALCL, $pibase_specs->{astral}->{seqcl}->{$seqid}) ;
         my $clusnum = 0 ;
         while (my $line = <ASTRALCL>) {
            chomp $line;
            if ($line =~ /^Rep/) {
               $line =~ s/^Rep: //g;
               $clusnum++ ; }

            if ($line =~ /score/) {
               $line =~ s/^\s+// ;
               my $scopid = substr($line, 0, 7) ;
               $seqcl->{$seqid}->{$scopid} = $clusnum ;
               push @{$seqcl2cont->{$seqid}->{$clusnum}}, $scopid ;
            }
         }
         close(ASTRALCL) ;
      }
   }

   $out->{seqcl} = $seqcl ;
   $out->{seqcl2cont} = $seqcl2cont ;

}


=head2 get_astral_classlist()

   Title:       get_astral_classlist()
   Function:    Get list of SCOP classes in the ASTRAL compendium
   Args:        $_->{pibase_specs} -  pibase_specs structure
   Returns:     Nothing - populates the specified $_->{out}
                ->{fam}->{scop_family} = number of domains in the family
                ->{sf}->{scop_superfamily} = number of domains in the superfamily

=cut

sub get_astral_classlist {

   my $in = shift;
   my $pibase_specs = $in->{pibase_specs} ;

   my $classlist ;
   foreach my $type (qw/fam sf/) {
      my @files=glob($pibase_specs->{asteroids}->{$type.'_aln'}."/*.fasta_aln");
      for (@files) {
         my $a = $_ ;
         $a =~ s/\.fasta_aln$// ;
         $a =~ s/.*\/// ;
         $classlist->{$type}->{$a}++ ;
       }
   }

   return $classlist ;
  
}


1 ;
