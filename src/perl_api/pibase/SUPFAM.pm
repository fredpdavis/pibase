=head1 NAME

pibase::SUPFAM.pm

=head1 DESCRIPTION

This module contains routines to interface with SUPFAM annotation files.
Goal is to map SCOP residue numbers onto target sequences using
SUPFAM alignments and ASTRAL mapping.

=head1 VERSION

fpd091013_0708

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

package pibase::SUPFAM ;
use strict;
use warnings;

use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/get_aln_string/ ;

use pibase ;
use pibase::pilig ;
use pibase::ASTRAL ;
use File::Temp qw/tempfile tempdir/ ;
use Bit::Vector ;
use fpdgenlib;
use fpdgenlib::SGE ;


=head2 set_supfam_specs

   Title:       set_supfam_specs()
   Function:    Sets configuration parameters
   Args:        None
   Returns:     $_->{option} = value; hash of parameters

=cut

sub set_supfam_specs {

   my $params ;

   my $pibase_specs = pibase::get_specs() ;
   $params->{pibase_specs} = $pibase_specs ;

   my $pilig_specs = pibase::pilig::set_pilig_specs() ;
   $params->{pilig_specs} = $pilig_specs ;

   $params->{standardres} = pibase::pilig::_list_standardres() ;

#class_levels => ['sf', 'fam'],

   $params->{thresholds} = {
      class_levels => ['fam'],
      'L' => { min_mw => 250, max_mw => 1000, bs_seqid => 0.3 },
      'P' => { bs_seqid => 0.4 },
      'p' => { bs_seqid => 0.4 },
   } ;

   $params->{scop_cla_fn} =
      $pibase_specs->{external_data}->{scop}->{files}->{cla} ;

   $params->{supfam_rootdir} =
      "/groups/eddy/home/davisf/work/databases/SUPFAM";
   $params->{selfhits_dir} = $params->{supfam_rootdir}.'/self_hits';

   $params->{ass_file_headers} = [
      'sp',
      'seq_id',
      'model_id',
      'res_range',
      'sf_eval',
      'alnstring',
      'fa_eval',
      'px_id',
      'fa_id'
   ] ;
   map {$params->{ass_file_f2i}->{$params->{ass_file_headers}->[$_]} = $_;}
      (0 .. $#{$params->{ass_file_headers}}) ;

   $params->{selfhits_headers} = [
      'model_id',
      'px_id',
      'alnstring',
   ] ;
   map {$params->{selfhits_f2i}->{$params->{selfhits_headers}->[$_]} = $_;}
      (0 .. $#{$params->{selfhits_headers}}) ;

# BLOSUM62 matrix: http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
   my $blosum62_matrix = "#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
" ;

   $params->{substitution_matrix} = readin_substitution_matrix({
      string => $blosum62_matrix}) ;

   return $params ;

}


=head2 readin_substitution_matrix

   Title:       readin_substitution_matrix()
   Function:    Parses a substitution matrix in matblas format
   Args:        ->{matrix_fn} = filename of substitution matrix
                ->{string} = string containing contents of matrix file

   Returns:     ->{raw}->{aa1}->{aa2} = substitution matrix score for aa1,aa2
                ->{nl}->{aa1}->{aa2} = normalized aa similarity score
                  (see karling_normalize_matrix() for normalization scheme)

=cut

sub readin_substitution_matrix {
   my $in = shift ;

   my @matrix_string ;
   if (exists $in->{string}) {
      @matrix_string = split(/\n/, $in->{string});
   } elsif (exists $in->{matrix_fn}) {
      open(MATF, $in->{matrix_fn});
      @matrix_string = <MATF> ;
      close(MATF);
   } else {
      die "substitution matrix not specified" ;
   }

   my @aa_order ;
   my $raw_matrix ;
   foreach my $mat_line (@matrix_string) {
      $mat_line =~ s/\s$// ;
      if ($mat_line =~ /^#/) {
         next; }
      my @t = split(' ', $mat_line) ;
      if ($mat_line !~ /[0-9]/) {
         @aa_order = @t ;
      } else {
         foreach my $j ( 1 .. $#t) {
            $raw_matrix->{$t[0]}->{$aa_order[($j -1)]} =
               $t[$j] ;
         }
      }
   }

   my $nl_matrix = karlin_normalize_matrix({
      matrix => $raw_matrix}); 

   return {
      raw => $raw_matrix,
      nl => $nl_matrix
   } ;
}


=head2 karlin_normalize_matrix

   Title:       karlin_normalize_matrix()
   Function:    Normalizes a substitution matrix per Karlin and Brocchieri,
                  J Bacteriol 1996; and rescaled to range from 0-1:

               0.5 * (1 + ( mat(i,j) / sqrt(|mat(i,i) * mat(j,j)|)))

   Args:        ->{matrix}->{aa1}->{aa2} = raw substitution matrix
   Returns:     ->{aa1}->{aa2} = normalized similarity score

=cut

sub karlin_normalize_matrix {

   my $in = shift ;
   my $raw_matrix = $in->{matrix};

   my $nl_matrix ;
   foreach my $aa1 ( keys %{$raw_matrix} ) {
      foreach my $aa2 ( keys %{$raw_matrix->{$aa1}} ) {
         $nl_matrix->{$aa1}->{$aa2} = ($raw_matrix->{$aa1}->{$aa2}) /
            sqrt( abs($raw_matrix->{$aa1}->{$aa1} *
                      $raw_matrix->{$aa2}->{$aa2})) ;

         $nl_matrix->{$aa1}->{$aa2} = ($nl_matrix->{$aa1}->{$aa2} + 1) / 2 ;
      }
   }

   return $nl_matrix ;
}


=head2 run_pilig_supfam_annotate()

   Title:       run_pilig_supfam_annotate()
   Function:    Maps PIBASE/LIGBASE binding sites onto target sequences
                  annotated with SUPERFAMILY domain assignments 

   Args:        ->{ARGV} = ARGV array reference; parsed to provide:
                ->{ass_fn} = name of SUPERFAMILY domain assignment file
                ->{out_fn} = name of output file
                ->{err_fn} = name of error file
                ->{matrix_fn} = optional substitution matrix, default BLOSUM62
                ->{cluster_fl} = run on an SGE cluster (options in SGE.pm)

   Returns:     NOTHING
   Displays:     1. seq_id
                 2. res_range
                 3. classtype
                 4. class
                 5. bs_type
                 6. bs_template
                 7. partner_descr
                 8. residues
                 9. bs_percseqident
                10. bs_percseqsim
                11. bs_numident
                12. bs_numgap
                13. bs_tmpl_numres
                14. bs_fracaln
                15. wholedom_numident
                16. wholedom_aln_length
                17. wholedom_percseqident
                18. wholedom_percseqsim

=cut

sub run_pilig_supfam_annotate {

   my $in = shift;

# Set parameters
   my $supfam_specs = set_supfam_specs() ;

# Set usage
   my $usage = __FILE__.
      " -ass_fn SUPERFAMILY_assignment_file\n".
      " [-out_fn output_file][-err_fn error_file]\n".
      " [-matrix_fn substitution_matrix_file] matblas format substitution matrix\n".
      " [-thresh_min_L_mw minimum_ligand_MW]\tdefault: 250\n".
      " [-thresh_max_L_mw maximum_ligand_MW]\tdefault: 1000\n".
      " [-thresh_L_bs_seqid N] ligand binding site sequence identity threshold, default: 0.3\n".
      " [-thresh_P_bs_seqid N] domain binding site sequence identity threshold, default: 0.4\n".
      " [-thresh_p_bs_seqid N] peptide binding site sequence identity threshold, default: 0.4";

   my @result_headers = qw/seq_id res_range classtype class bs_type bs_template partner_descr residues bs_percseqident bs_percseqsim bs_numident bs_numgap bs_tmpl_numres bs_fracaln wholedom_numident wholedom_aln_length wholedom_percseqident wholedom_percseqsim/ ;


# Read in command line options
   my $j = 0; my $user_thresh = {};
   while ($j < $#{$in->{ARGV}}) {
      $ARGV[$j] =~ s/^\-// ;
      if ($ARGV[$j] =~ /^thresh/) {
         $user_thresh->{$in->{ARGV}->[$j]} = $in->{ARGV}->[($j+1)] ;
      } else {
         $in->{$in->{ARGV}->[$j]} = $in->{ARGV}->[($j+1)] ;
      }
      $j += 2 ;
   }

   if (exists $in->{summarize_results}) {
      summarize_results({ results_fn => $in->{summarize_results} }) ;
      exit;
   }

   my $fpdgenlib_specs ;
   if (exists $in->{fpdgenlib_specs}) {
      $fpdgenlib_specs = $in->{fpdgenlib_specs} ;
   } else {
      $fpdgenlib_specs = fpdgenlib::get_specs() ;
   }


# Read in matrix file if user specified
   if (exists $in->{matrix_fn}) {
      $supfam_specs->{substitution_matrix} =
         readin_substitution_matrix({fn => $in->{matrix_fn}}) ;
   }

# Over-write default thresholds if user specified
   foreach my $type (keys %{$user_thresh}) {
      if ($type eq 'thresh_min_L_mw') {

         $supfam_specs->{thresholds}->{L}->{min_mw} = $user_thresh->{$type} ;

      } elsif ($type eq 'thresh_max_L_mw') {

         $supfam_specs->{thresholds}->{L}->{max_mw} = $user_thresh->{$type} ;

      } elsif ($type eq 'thresh_L_bs_seqid') {

         if ($user_thresh->{type} < 0 || $user_thresh->{type} > 1) {
            die "seqid threshold has to range from 0 to 1\nUSAGE: $usage" ; }
         $supfam_specs->{thresholds}->{L}->{bs_seqid} = $user_thresh->{$type} ;

      } elsif ($type eq 'thresh_P_bs_seqid') {

         if ($user_thresh->{type} < 0 || $user_thresh->{type} > 1) {
            die "seqid threshold has to range from 0 to 1\nUSAGE: $usage" ; }
         $supfam_specs->{thresholds}->{P}->{bs_seqid} = $user_thresh->{$type} ;

      } elsif ($type eq 'thresh_p_bs_seqid') {

         if ($user_thresh->{type} < 0 || $user_thresh->{type} > 1) {
            die "seqid threshold has to range from 0 to 1\nUSAGE: $usage" ; }
         $supfam_specs->{thresholds}->{p}->{bs_seqid} = $user_thresh->{$type} ;

      }
   }

   if (!exists $in->{ass_fn} || !-s $in->{ass_fn}) {
      die "SUPFAM annotation file not found (-ass_fn)\nUSAGE: $usage" ; }

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1 ) { # master node:

      print "* run_pilig_supfam_annotate() ".localtime()
         if(!exists $in->{quiet_fl});
      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{run_pilig_supfam_annotate_in},
       $temp_fn->{run_pilig_supfam_annotate_in}) =
       tempfile("splits_run_pilig_supfam_annotate_input.XXXXX") ;

      ($temp_fh->{run_pilig_supfam_annotate_out},
       $temp_fn->{run_pilig_supfam_annotate_out}) =
       tempfile("splits_run_pilig_supfam_annotate_SGEout.XXXXX") ;
      ($temp_fh->{run_pilig_supfam_annotate_err},
       $temp_fn->{run_pilig_supfam_annotate_err}) =
       tempfile("splits_run_pilig_supfam_annotate_SGEerr.XXXXX") ;

# Sort ASS file by family identifier - to reduce ASTRAL aln read in time
      my $tcom_sort_input = "cat ".$in->{ass_fn}." | sort -t'	' -k9,9nr >".
                            $temp_fn->{run_pilig_supfam_annotate_in} ;
      system($tcom_sort_input) ;

      my $split_dir = tempdir("splits_run_pilig_supfam_annotate.XXXXX") ;
      my $numjobs ;
      if (exists $in->{numjobs}) {
         $numjobs = $in->{numjobs} ;
      } else {
         $numjobs = $fpdgenlib_specs->{SGE}->{numjobs} ;
      }

      my $splits = fpdgenlib::SGE::_clust_split_ins({
         fn => $temp_fn->{run_pilig_supfam_annotate_in},
         dir => $split_dir,
         numjobs => $numjobs,
      }) ;

      my ($perlscript_fh, $perlscript_fn) =
         tempfile("pb.run_pilig_supfam_annotate.XXXXX", SUFFIX =>'.ami.pl') ;


      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase;
use fpdgenlib;
use pibase::SUPFAM ;
use Bit::Vector ;

main() ;

sub main {

   pibase::SUPFAM::run_pilig_supfam_annotate({
      cluster_fl => 0,
      ARGV => \\\@ARGV,
   }) ;

}

" ;


      my ($sgescript_fh, $sgescript_fn)  =
         tempfile("pb.run_pilig_supfam_annotate.XXXXX", SUFFIX => ".SGE.sh") ;
      my $sge_outdir = tempdir("SGEOUT.run_pilig_supfam_annotate.XXXXX") ;

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $fpdgenlib_specs->{SGE}->{nodespecs}) {
         print {$sgescript_fh} $fpdgenlib_specs->{SGE}->{nodespecs}."\n" ;}

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1 = \$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn -ass_fn \$input1 

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print STDERR "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = fpdgenlib::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;
      print STDERR "      job $qsub_job_id\n" ;

      while (1) {
         sleep $fpdgenlib_specs->{SGE}->{qstat_sleep} ;
         my $job_status= fpdgenlib::SGE::_clust_qstat({job_id => $qsub_job_id});
         if ($job_status) {last;}
      }

      fpdgenlib::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{run_pilig_supfam_annotate_out},
         err_fn => $temp_fn->{run_pilig_supfam_annotate_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;


      if (exists $in->{out_fn}) {
# remove internal header lines.
         open(REALOUTF, ">".$in->{out_fn}) ;
         open(FULLOUTF, $temp_fn->{run_pilig_supfam_annotate_out}) ;
         my $header_line = '#'.join("\t", @result_headers) ;
         print REALOUTF $header_line."\n" ;
         while (my $line = <FULLOUTF>) {
            chomp $line;
            if ($line eq $header_line) {next;}
            print REALOUTF $line."\n" ;
         }
         close(FULLOUTF) ;
         close(REALOUTF) ;
#         system("mv ".$temp_fn->{run_pilig_supfam_annotate_out}." ".
#                $in->{out_fn}) ;
         unlink $temp_fn->{run_pilig_supfam_annotate_out} ;
      }

      if (exists $in->{err_fn}) {
         system("mv ".$temp_fn->{run_pilig_supfam_annotate_err}." ".
                $in->{err_fn}) ; }

   } else { #processing node:

# if cluster run specified, split the input file and submit them individually
#  - nothing complicated...

      my $out_fh ;
      if (exists $in->{out_fn}) {
         open($out_fh, ">".$in->{out_fn}) ;
      } else {
         open($out_fh, ">-") ;
      }

      if (exists $in->{err_fn}) {
         open(STDERR, ">".$in->{err_fn}) ; }

      my $pibase_specs = $supfam_specs->{pibase_specs} ;
      my $pilig_specs = $supfam_specs->{pilig_specs} ;
      my $standardres = $supfam_specs->{standardres};

#Preload some ASTRAL data
      my $astral = pibase::pilig::_pilig_astral_preload({
         pibase_specs => $pibase_specs }) ;

# Load PIBASE and LIGBASE binding site info
      my $pb = pibase::pilig::_pilig_tod_pibase_preload() ;
      my $liginfo = pibase::pilig::_pilig_load_liginfo() ;
      my $class2alnlength = {};

      my $expbits = pibase::pilig::readin_expassignments({
         fn => $pilig_specs->{outfiles}->{assign_exp},
         standardres => $standardres,
         class2alnlength => $class2alnlength,
      }) ;

      my $pepnucibits = pibase::pilig::readin_pepnuciassignments({
         fn => $pilig_specs->{outfiles}->{assign_pepnuci_clusters},
         clustrep_fl => 1,
         pb => $pb,
         pilig_specs => $pilig_specs,
         standardres => $standardres,
         class2alnlength => $class2alnlength,
      }) ;
      my $ligbits = pibase::pilig::readin_ligassignments({
         fn => $pilig_specs->{outfiles}->{assign_lig_clusters},
         clustrep_fl => 1,
         liginfo => $liginfo,
         standardres => $standardres,
         class2alnlength => $class2alnlength,
      }) ;

      my $pibits_both = pibase::pilig::readin_piassignments({
         fn => $pilig_specs->{outfiles}->{assign_pi_clusters},
         clustrep_fl => 1,
         pb => $pb,
         dont_read_jdomains_fl => 1,
         standardres => $standardres,
         class2alnlength => $class2alnlength,
      }) ;
      my $pibits = $pibits_both->{pibits} ;
      my $interfaces = $pibits_both->{interfaces} ;

      my $alltogether = pibase::pilig::combine_piligpepexpbits({
         pb => $pb,
         liginfo => $liginfo,
         ligbits => $ligbits,
         pibits => $pibits,
         pepbits => $pepnucibits->{pep},
         expbits => $expbits,
      });
      my $class2bits = $alltogether->{class2bits} ;
      my $class2ligs = $alltogether->{class2ligs} ;
      my $class2allligbits = $alltogether->{class2allligbits} ;

# create a family indexed list of binding sites
      my $class2sid12_side ;
      my $class2pint;
      foreach my $classtype (@{$supfam_specs->{thresholds}->{class_levels}}) {
         foreach my $sid12 (keys %{$interfaces}) {
            my $sid ;
            ($sid->{1}, $sid->{2}) = split(/\t/, $sid12) ;
            foreach my $side (1,2) {
               $class2sid12_side->{$pb->{sid2class}->{$classtype}->{$sid->{$side}}}->{$sid12."\t".$side}++ ; }
         }

         foreach my $sid (keys %{$pepnucibits->{pep}->{$classtype}}) {
            $class2pint->{$pb->{sid2class}->{$classtype}->{$sid}}->{$sid}++ ; }
      }


#Read in domain information from scop_cla: fam2sf, px2scopid, scopid2class
      my $scopinfo = readin_scop_cla({fn => $supfam_specs->{scop_cla_fn}}) ;


#NOTE: sort ASSF by domain fam/sf so less time spent reading in asteroids alns.
      print STDERR "NOW performing homology transfer\n" ;
      my $cur_aln = {};
      open(ASSF, $in->{ass_fn}) ;

# Print headers
      print {$out_fh} '#'.join("\t", @result_headers)."\n";

# Read in a target-SUPFAM assignment entry
      while (my $line = <ASSF>) {
         chomp $line;
         if ($line =~ /^#/) {next;}
         my @t = split(/\t/, $line) ;
         my $curass ;
         map {$curass->{$supfam_specs->{ass_file_headers}->[$_]} = $t[$_]}
            (0 .. $#t) ;

# 0. Check if this domain class is covered by ASTRAL
         my $tmpl_scopid = $scopinfo->{px2scopid}->{$curass->{px_id}} ;
         my $cur_class;
         $cur_class->{fam} = $scopinfo->{scopid2class}->{$tmpl_scopid} ;
         $cur_class->{sf} = $cur_class->{fam} ;
         $cur_class->{sf} =~ s/\.[0-9]+$// ;

# (ASTRAL alignments only cover classes a-g)
         if ($cur_class->{fam} !~ /[a-g]/) {
            my @outvals = ($curass->{seq_id}, $curass->{res_range},
                           'fam', $cur_class->{fam}, 'all',
                           'domain family not covered by ASTRAL') ;
            print {$out_fh} '#'.join("\t", @outvals)."\n" ;
            next;
         }

# 1. Get seq-SUPFAM template alignment: read line from self_hits/sf_id.hits with ^submodel_id px_id (map from target_seq to tmpl_seq)
         $curass->{sf_id} = $scopinfo->{faid2sfid}->{$curass->{fa_id}} ;

         my $supfam_aln ;
         $supfam_aln->{seq} = $curass->{alnstring} ;
         my $target_seq = $supfam_aln->{seq}; $target_seq =~ s/\-// ;

         $supfam_aln->{strx} = get_SUPFAM_selfhit({
            supfam_specs => $supfam_specs,
            sf_id => $curass->{sf_id},
            model_id => $curass->{model_id},
            px_id => $curass->{px_id},
         }) ;

# Exit this domain if the template self-hit string not found
         if (!defined $supfam_aln->{strx}) {
#            print STDERR "ERROR (seqid ".$curass->{seq_id}.", ".
#                  "res_range ".$curass->{res_range}."): ".
#                  "self-hit alignment string not found for ".
#                  "px ".$curass->{px_id}.", ".
#                  "sf ".$curass->{sf_id}.", ".
#                  "model ".$curass->{model_id}."\n" ;
            my @outvals = ($curass->{seq_id}, $curass->{res_range},
                           'fam', $cur_class->{fam}, 'all',
                           'ERROR: self-hit alignment not found for '.
                           "px ".$curass->{px_id}.", ".
                           "sf ".$curass->{sf_id}.", ".
                           "model ".$curass->{model_id}) ;
            print {$out_fh} '#'.join("\t", @outvals)."\n" ;
            next;
         }

# 2. Get SUPFAM template's ASTRAL alignment string from ASTRAL
         my $astral_aln;
         my $target_astral_aln = {} ;
         my $abort_target = 0 ;
         foreach my $classtype (@{$supfam_specs->{thresholds}->{class_levels}}){
            my $class = $cur_class->{$classtype} ;

# Load ASTRAL alignment if necessary (ie, new fam or sf)
            if (!exists $cur_aln->{"cur_".$classtype}  ||
                $cur_aln->{"cur_".$classtype} ne $cur_class->{$classtype}) {
               $cur_aln->{$classtype."_aln"} = 
                pibase::ASTRAL::load_asteroids_aln({
                  aln_fn => $pibase_specs->{asteroids}->{$classtype.'_aln'}.
                     '/'.$cur_class->{$classtype}.'.fasta_aln' ,
                  seq_fn => $pibase_specs->{asteroids}->{$classtype.'_seq'}.
                     '/'.$cur_class->{$classtype}.'.fa' ,
                  raf => $astral->{raf},
                  gdseqh => $astral->{gdseqh},
                  seqclcont100 => $astral->{seqcl2cont}->{100},
                  seqcl100 => $astral->{seqcl}->{100},
                  allchains => $pb->{pdbchains}
               }) ;
               $cur_aln->{"cur_".$classtype} = $cur_class->{$classtype} ;
            }

            $astral_aln->{$classtype} = {
               strx => $cur_aln->{$classtype."_aln"}->{aln}->{$tmpl_scopid}} ;

            if (!defined $astral_aln->{$classtype}->{strx}) {
               my @outvals = ($curass->{seq_id}, $curass->{res_range},
                           $classtype, $cur_class->{$classtype}, 'all',
                           "ERROR: no $classtype ASTRAL aln for ".
                           "$tmpl_scopid SUPERFAMILY template domain");
               print {$out_fh} '#'.join("\t", @outvals)."\n";
               next;
            }

# Merge the two alignments.; NOTE: SUPFAM aln is not a true aln - includes
#  lower case characters that are not actually part of the alignment.
            $target_astral_aln->{$classtype} = merge_SUPFAM_ASTRAL_alignments({
               superfam_aln => $supfam_aln,
               astral_aln => $astral_aln->{$classtype},
               common => 'strx',
               target => 'seq',
            }) ;

            if (exists $target_astral_aln->{$classtype}->{error_fl}) {

               my @outvals = ($curass->{seq_id}, $curass->{res_range},
                           $classtype, $cur_class->{$classtype}, 'all',
                           "ERROR in merging SUPFAM/ASTRAL alignments, ".
                           "px ".$curass->{px_id}.", ".
                           "sf ".$curass->{sf_id}.", ".
                           "model ".$curass->{model_id}.": ".
                           $target_astral_aln->{$classtype}->{error_fl});
               print {$out_fh} '#'.join("\t", @outvals)."\n";

               $abort_target = 1; 
               last;
            }
         }

         if ($abort_target) {next;}

# Iterate over the class level for homology transfer, from (fam, sf)
         foreach my $classtype (@{$supfam_specs->{thresholds}->{class_levels}}){
            my $class = $cur_class->{$classtype} ;

# Initialize warnings to display if no annotations are made
         my $fail_type = {
            L => 'no templates',
            p => 'no templates',
            P => 'no templates',
         } ;

# Iterate over all ligand binding sites in this class;
#   homology transfer to target seq and calc bs/dom-seqid
#  (code adapated from pibase::pilig::collate_per_instance)
         my $target_bs ;
         if (exists $class2allligbits->{$classtype}->{$class}) {
            my $curligs = $class2allligbits->{$classtype}->{$class};

# Iterate over all ligand binding sites in this class;
            foreach my $j (0 .. $#{$curligs}) {
#               print STDERR "   NOW ON LIGAND SITE $j\n" ;
               my ($sid_origdom_l) =
                  ($curligs->[$j]->[2] =~ /SCOP\.(.+)/) ;
               my $curligsig = $curligs->[$j]->[0] ;
               my (undef, $ligcod, undef) = split(/\t/, $curligsig) ;
               my $outcurligsig = $curligsig; $outcurligsig =~ s/\t/:/g ;
               my $curligbits = $curligs->[$j]->[1] ;
               my $curligsid = $curligs->[$j]->[2] ; #fpd090505_1049 

#Check if ligand is within molecular weight range, if specified
               if (exists $supfam_specs->{thresholds}->{L}->{min_mw} &&
                    $liginfo->{mw}->{$ligcod} <
                       $supfam_specs->{thresholds}->{L}->{min_mw}) {next;}

               if (exists $supfam_specs->{thresholds}->{L}->{max_mw} &&
                    $liginfo->{mw}->{$ligcod} >
                       $supfam_specs->{thresholds}->{L}->{max_mw}) {next;}

               $fail_type->{L} = 'sub-threshold templates';

               my $bs_sig = $sid_origdom_l.":".$outcurligsig ;
               $target_bs->{L}->{$bs_sig}->{num_aln} = 0 ;
               if (exists $liginfo->{mw}->{$ligcod} &&
                   $liginfo->{mw}->{$ligcod} ne '') {
                  $target_bs->{L}->{$bs_sig}->{lig_mw} =
                     sprintf("%.3f",$liginfo->{mw}->{$ligcod});
               } else {
                  $target_bs->{L}->{$bs_sig}->{lig_mw} = '' ;
               }
               $target_bs->{L}->{$bs_sig}->{descr_field} =
                  'MW='.$target_bs->{L}->{$bs_sig}->{lig_mw} ;

#               print STDERR "CHECKING TRANSFER of $outcurligsig to ".
#                  $curass->{seq_id}."\n" ;

# iterate over all alignment positions to calc whole dom seqid
               $target_bs->{L}->{$bs_sig}->{wholedom_aln_length} = 
                 length($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_l});
               $target_bs->{L}->{$bs_sig}->{wholedom_sum_simscore} = 0 ;
               $target_bs->{L}->{$bs_sig}->{wholedom_num_ident} = 0 ;

               foreach my $cur_alnpos (0 ..
                  length($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_l})
                         - 1) {

#                  print STDERR " ASTRAL alnpos $cur_alnpos\n" ;
#                  print STDERR "   = template res ".
#                  substr($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_l},
#                         $cur_alnpos,1)."\n";

                  if (!exists 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos}) {
#                     print STDERR "   target seq has a gap at this position\n";
                     next;
                  }
                  my $target_res = 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos} ;

#                  print STDERR "   = target seq res ".
#   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos}.
#                       " ".substr($target_seq, ($target_res - 1),1)."\n";

                  my $cur_tmpl_char = substr(
                     $cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_l},
                     $cur_alnpos,1) ;

                  my $cur_targ_char = substr($target_seq, ($target_res - 1),1);
                  if ($cur_tmpl_char eq $cur_targ_char) {
                     $target_bs->{L}->{$bs_sig}->{wholedom_num_ident}++ ; }

                  if (exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char} &&
                      exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char}) {
                     $target_bs->{L}->{$bs_sig}->{wholedom_sum_simscore} +=
$supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char};
                  }
               }

# iterate over the binding alingment positions and transfer them onto target
               my @bs_alnpos = $curligbits->Index_List_Read() ;
               $target_bs->{L}->{$bs_sig}->{residues} = [] ;
               foreach my $bs_alnpos (@bs_alnpos) {
                  if (!exists
                   $target_astral_aln->{$classtype}->{alnpos2resno}->{$bs_alnpos}){
                     next;}

                  my $target_res = 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$bs_alnpos} ;
                  push @{$target_bs->{L}->{$bs_sig}->{residues}}, $target_res ;
                  $target_bs->{L}->{$bs_sig}->{num_aln}++ ;

                  my $cur_tmpl_char = substr($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_l}, $bs_alnpos,1) ;
                  my $cur_targ_char = substr($target_seq, ($target_res - 1),1);
#check if identical
                  if ($cur_tmpl_char eq $cur_targ_char) {
                     $target_bs->{L}->{$bs_sig}->{num_ident}++ ; }

                  if (exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char} &&
                      exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char}) {
                     $target_bs->{L}->{$bs_sig}->{sum_simscore} +=
$supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char};
                  }
               }

               if ($#{$target_bs->{L}->{$bs_sig}->{residues}} >= 0) {
                  if (!exists $target_bs->{L}->{$bs_sig}->{num_ident} ) {
                     $target_bs->{L}->{$bs_sig}->{num_ident} = 0 ;}
                  if (!exists $target_bs->{L}->{$bs_sig}->{sum_simscore} ) {
                     $target_bs->{L}->{$bs_sig}->{sum_simscore} = 0 ;}
                  $target_bs->{L}->{$bs_sig}->{tmpl_bs_size} = $#bs_alnpos + 1 ;
                  $target_bs->{L}->{$bs_sig}->{reslist} = join(", ",
                     @{$target_bs->{L}->{$bs_sig}->{residues}}) ;
               } else {
                  delete $target_bs->{L}->{$bs_sig} ;
               }
            }
         }

# Iterate over all protein binding sites in this class;
#   homology transfer to target seq and calc bs/dom-seqid
         if (exists $class2sid12_side->{$class}) {
            foreach my $sid12_side (keys %{$class2sid12_side->{$class}}) {
               my ($sid, $side) ;
               ($sid->{1}, $sid->{2}, $side) = split(/\t/, $sid12_side);
               my $sid12 = $sid->{1}."\t".$sid->{2} ;
               my $classes ;
               $classes->{1} = $pb->{sid2class}->{'fam'}->{$sid->{1}} ;
               $classes->{2} = $pb->{sid2class}->{'fam'}->{$sid->{2}} ;
               my ($sid_origdom_p) = ($sid->{$side} =~ /SCOP\.(.+)/) ;
               my $bs_sig = $sid12_side; $bs_sig =~ s/\t/:/g ;

               $target_bs->{P}->{$bs_sig}->{partner_class} = $classes->{2} ;
               if ($side == 2) {
                  $target_bs->{P}->{$bs_sig}->{partner_class} = $classes->{1}; }

               $target_bs->{P}->{$bs_sig}->{num_aln} = 0 ;
               $target_bs->{P}->{$bs_sig}->{residues} = [] ;

#               print STDERR "   NOW ON PROTEIN SITE $bs_sig\n" ;

# get alignment positions for this binding site.
               if (!exists $interfaces->{$sid12}->{$side} ||
                !exists $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}||
                ($interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Norm()
                    == 0)) {next;}

               $fail_type->{P} = 'sub-threshold templates';

               $target_bs->{P}->{$bs_sig}->{chains} =
                  $interfaces->{$sid12}->{chains} ;
               $target_bs->{P}->{$bs_sig}->{descr_field} =
                  'chains='.$target_bs->{P}->{$bs_sig}->{chains}.';'.
                  'class='.$target_bs->{P}->{$bs_sig}->{partner_class} ;

               my @bs_alnpos =
   $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Index_List_Read() ;

# iterate over all alignment positions to calc whole dom seqid
               $target_bs->{P}->{$bs_sig}->{wholedom_aln_length} =
                 length($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p});
               $target_bs->{P}->{$bs_sig}->{wholedom_num_ident} = 0 ;
               $target_bs->{P}->{$bs_sig}->{wholedom_sum_simscore} = 0 ;

               foreach my $cur_alnpos (0 ..
                  length($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p})
                         - 1) {

#                  print STDERR " ASTRAL alnpos $cur_alnpos\n" ;
#                  print STDERR "   = template res ".
#                  substr($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p},
#                         $cur_alnpos,1)."\n";

                  if (!exists 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos}) {
#                     print STDERR "   target seq has a gap at this position\n";
                     next;
                  }
                  my $target_res = 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos} ;

                  my $cur_tmpl_char = substr(
                     $cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p},
                     $cur_alnpos,1) ;
                  my $cur_targ_char = substr($target_seq, ($target_res - 1),1);

#                  print STDERR "   = target seq res ".
#   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos}.
#                       " ".substr($target_seq, ($target_res - 1),1)."\n";
                  if ($cur_tmpl_char eq $cur_targ_char) {
                     $target_bs->{P}->{$bs_sig}->{wholedom_num_ident}++ ; }

                  if (exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char} &&
                      exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char}) {
                     $target_bs->{P}->{$bs_sig}->{wholedom_sum_simscore} +=
$supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char};
                  }
               }


# iterate over alnpos2targetseqresno to homology transfer the bs
               foreach my $bs_alnpos (@bs_alnpos) {
                  if (!exists
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$bs_alnpos}){ next;}

                  my $target_res = 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$bs_alnpos} ;
                  push @{$target_bs->{P}->{$bs_sig}->{residues}}, $target_res ;
                  $target_bs->{P}->{$bs_sig}->{num_aln}++ ;

                  my $cur_tmpl_char = substr(
                     $cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p},
                     $bs_alnpos,1) ;
                  my $cur_targ_char = substr($target_seq, ($target_res - 1),1);

#check if identical
                  if ($cur_tmpl_char eq $cur_targ_char) {
                     $target_bs->{P}->{$bs_sig}->{num_ident}++ ; }

                  if (exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char} &&
                      exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char}) {
                     $target_bs->{P}->{$bs_sig}->{sum_simscore} +=
$supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char};
                  }
               }

               if ($#{$target_bs->{P}->{$bs_sig}->{residues}} >= 0) {
                  if (!exists $target_bs->{P}->{$bs_sig}->{num_ident} ) {
                     $target_bs->{P}->{$bs_sig}->{num_ident} = 0 ;}
                  if (!exists $target_bs->{P}->{$bs_sig}->{sum_simscore} ) {
                     $target_bs->{P}->{$bs_sig}->{sum_simscore} = 0 ;}
                  $target_bs->{P}->{$bs_sig}->{tmpl_bs_size} = $#bs_alnpos + 1 ;
                  $target_bs->{P}->{$bs_sig}->{reslist} = join(", ",
                     @{$target_bs->{P}->{$bs_sig}->{residues}}) ;
               } else {
                  delete $target_bs->{P}->{$bs_sig} ;
               }
            }
         }


# Iterate over all peptide binding sites in this class;
#   homology transfer to target seq and calc bs/dom-seqid
         if (exists $class2pint->{$class}) {
            foreach my $sid (keys %{$class2pint->{$class}}) {
               my ($sid_origdom_p) = ($sid =~ /SCOP\.(.+)/) ;
            foreach my $targetch (
                keys %{$pepnucibits->{pep}->{$classtype}->{$sid}} ){
                if ($targetch eq 'cumulative') {next;}

               my $bs_sig = $sid.":".$targetch ;
               $target_bs->{p}->{$bs_sig}->{num_aln} = 0 ;
               $target_bs->{p}->{$bs_sig}->{residues} = [] ;
#            print STDERR "   NOW ON PEPTIDE SITE $bs_sig\n" ;

# get alignment positions for this binding site.
               my @bs_alnpos = 
   $pepnucibits->{pep}->{$classtype}->{$sid}->{$targetch}->Index_List_Read();

               $fail_type->{p} = 'sub-threshold templates';

# iterate over all alignment positions to calc whole dom seqid
               $target_bs->{p}->{$bs_sig}->{wholedom_aln_length} = 
                 length($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p});
               $target_bs->{p}->{$bs_sig}->{wholedom_num_ident} = 0 ;
               $target_bs->{p}->{$bs_sig}->{wholedom_sum_simscore} = 0 ;

               $target_bs->{p}->{$bs_sig}->{chain_length} =
                  $pepnucibits->{chain_info}->{$targetch}->{chain_length};
               $target_bs->{p}->{$bs_sig}->{descr_field} =
                  'peptide_length='.$target_bs->{p}->{$bs_sig}->{chain_length};

               foreach my $cur_alnpos (0 ..
                  length($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p})
                         - 1) {

#                  print STDERR " ASTRAL alnpos $cur_alnpos\n" ;
#                  print STDERR "   = template res ".
#                  substr($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p},
#                         $cur_alnpos,1)."\n";

                  if (!exists 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos}) {
#                     print STDERR "   target seq has a gap at this position\n";
                     next;
                  }
                  my $target_res = 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos} ;

#                  print STDERR "   = target seq res ".
#   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos}.
#                       " ".substr($target_seq, ($target_res - 1),1)."\n";
                  my $cur_tmpl_char = substr(
                     $cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p},
                     $cur_alnpos,1) ;
                  my $cur_targ_char = substr($target_seq, ($target_res - 1),1);

                  if ($cur_tmpl_char eq $cur_targ_char) {
                     $target_bs->{p}->{$bs_sig}->{wholedom_num_ident}++ ; }

                  if (exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char} &&
                      exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char}) {
                     $target_bs->{p}->{$bs_sig}->{wholedom_sum_simscore} +=
$supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char};
                  }
               }


# iterate over alnpos2targetseqresno to homology transfer the bs
               foreach my $bs_alnpos (@bs_alnpos) {
#                  print STDERR "      alnpos: $bs_alnpos\n" ;
                  if (!exists
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$bs_alnpos}){ next;}

                  my $target_res = 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$bs_alnpos} ;

#                  print STDERR "         == target_res $target_res\n";
                  push @{$target_bs->{p}->{$bs_sig}->{residues}}, $target_res ;
                  $target_bs->{p}->{$bs_sig}->{num_aln}++ ;

                  my $cur_tmpl_char = substr(
                     $cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p},
                     $bs_alnpos,1) ;
                  my $cur_targ_char = substr($target_seq, ($target_res - 1),1);

#check if identical
                  if ($cur_tmpl_char eq $cur_targ_char) {
                     $target_bs->{p}->{$bs_sig}->{num_ident}++ ; }

                  if (exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char} &&
                      exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char}) {
                     $target_bs->{p}->{$bs_sig}->{sum_simscore} +=
$supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char};
                  }
               }

               if ($#{$target_bs->{p}->{$bs_sig}->{residues}} >= 0) {
                  if (!exists $target_bs->{p}->{$bs_sig}->{num_ident} ) {
                     $target_bs->{p}->{$bs_sig}->{num_ident} = 0 ;}
                  if (!exists $target_bs->{p}->{$bs_sig}->{sum_simscore} ) {
                     $target_bs->{p}->{$bs_sig}->{sum_simscore} = 0 ;}
                  $target_bs->{p}->{$bs_sig}->{tmpl_bs_size} = $#bs_alnpos + 1 ;
                  $target_bs->{p}->{$bs_sig}->{reslist} = join(", ",
                     @{$target_bs->{p}->{$bs_sig}->{residues}}) ;
               } else {
                  delete $target_bs->{p}->{$bs_sig} ;
               }
            }
            }
         }

# Display transferred binding sites
         foreach my $bs_type (qw/L p P/) {
            if (! exists $target_bs->{$bs_type}) {
                  my @outvals = ($curass->{seq_id}, $curass->{res_range},
                                 $classtype, $cur_class->{$classtype},
                                 $bs_type, $fail_type->{$bs_type}) ;
               print {$out_fh} '#'.join("\t", @outvals)."\n" ;
               next;
            }

            my $num_transf = 0 ;
            foreach my $bs_sig (keys %{$target_bs->{$bs_type}}) {
               if (! exists $target_bs->{$bs_type}->{$bs_sig}->{num_aln} ||
                   $target_bs->{$bs_type}->{$bs_sig}->{num_aln} == 0) {
                  next; }
               my $wholedom_percident = sprintf("%.3f",
                  ($target_bs->{$bs_type}->{$bs_sig}->{wholedom_num_ident} / 
                   $target_bs->{$bs_type}->{$bs_sig}->{wholedom_aln_length})) ;
               my $bs_percident = sprintf("%.3f",
                  ($target_bs->{$bs_type}->{$bs_sig}->{num_ident} / 
                   $target_bs->{$bs_type}->{$bs_sig}->{tmpl_bs_size})) ;

               my $wholedom_percsim = sprintf("%.3f",
                  ($target_bs->{$bs_type}->{$bs_sig}->{wholedom_sum_simscore} / 
                   $target_bs->{$bs_type}->{$bs_sig}->{wholedom_aln_length})) ;
               my $bs_percsim = sprintf("%.3f",
                  ($target_bs->{$bs_type}->{$bs_sig}->{sum_simscore} / 
                   $target_bs->{$bs_type}->{$bs_sig}->{tmpl_bs_size})) ;

               if (exists $supfam_specs->{thresholds}->{$bs_type}->{bs_seqid} &&
                   $bs_percident < 
                     $supfam_specs->{thresholds}->{$bs_type}->{bs_seqid}){next;}

               $num_transf++ ;
               my $bs_fracaln = sprintf("%.3f",
                  ($target_bs->{$bs_type}->{$bs_sig}->{num_aln} / 
                   $target_bs->{$bs_type}->{$bs_sig}->{tmpl_bs_size})) ;
               my @outvals = ($curass->{seq_id}, $curass->{res_range},
                              $classtype, $cur_class->{$classtype},
                              $bs_type, $bs_sig,
                              $target_bs->{$bs_type}->{$bs_sig}->{descr_field},
                              $target_bs->{$bs_type}->{$bs_sig}->{reslist},
                              $bs_percident,
                              $bs_percsim,
                              $target_bs->{$bs_type}->{$bs_sig}->{num_ident},
                              $target_bs->{$bs_type}->{$bs_sig}->{num_aln},
                              $target_bs->{$bs_type}->{$bs_sig}->{tmpl_bs_size},
                              $bs_fracaln,
   $target_bs->{$bs_type}->{$bs_sig}->{wholedom_num_ident},
   $target_bs->{$bs_type}->{$bs_sig}->{wholedom_aln_length},
                              $wholedom_percident,
                              $wholedom_percsim
                             );
               print {$out_fh} join("\t", @outvals)."\n" ;
            }
            if ($num_transf == 0) {
                  my @outvals = ($curass->{seq_id}, $curass->{res_range},
                                 $classtype, $cur_class->{$classtype},
                                 $bs_type, $fail_type->{$bs_type}) ;
               print {$out_fh} '#'.join("\t", @outvals)."\n" ;
            }
         }
         }
      }
   }

}


=head2 merge_SUPFAM_ASTRAL_alignments

   Title:       merge_SUPFAM_ASTRAL_alignments()
   Function:    Merges SUPERFAMILY alignment string with ASTRAL alignment to
                 get SUPERFAMILY annotated target sequence in the ASTRAL
                 alignment frame (where it can receive binding site annotations)

   Args:        ->{superfam_aln} = SUPERFAMILY alignment string
                ->{astral_aln} = ASTRAL alignment string
                ->{target} = target sequence name
                ->{common} = template sequence name present in both alignments

   Returns:     ->{resno2alnpos}->{target resno} = alignment position
                ->{alnpos2resno}->{alignment position} = target resno

=cut

sub merge_SUPFAM_ASTRAL_alignments {

   my $in = shift ;
   my $aln1= $in->{superfam_aln} ;
   my $aln2 = $in->{astral_aln} ;
   my $target = $in->{target} ;
   my $common = $in->{common} ;

# CHECK 1. MAKE SURE THE COMMON SEQUENCES ARE ACTUALLY IDENTICAL!
   my $seq1 = $aln1->{$common} ; $seq1 =~ s/\-//g ; $seq1 = uc($seq1) ;
   my $seq2 = $aln2->{$common} ; $seq2 =~ s/\-//g ; $seq2 = uc($seq2) ;
   if ($seq1 ne $seq2) {
      return {error_fl => "the \"common\" sequence differs between the two ".
              "alignments; CANNOT MERGE!\tseq1: $seq1\tseq2: $seq2"} ;
   }

#   print STDERR "trying to merge:\nTMPL1: ".$aln1->{$common}."\n" ;
#   print STDERR "TARG1: ".$aln1->{$target}."\n" ;
#   print STDERR "           with:\nALN2: ".$aln2->{$common}."\n" ;

# ALIGNMENT 1: SUPFAM alignment string
# NOTE NOT A TRUE ALIGNMENT: skip residues in lower case
   my $aln1_resno2alnpos = {}; my $aln1_alnpos2resno = {};
   foreach my $seq ($target, $common) {
      my $cur_resno = 0 ; my $cur_alnpos = 0 ;
      foreach my $j (0 .. (length($aln1->{$seq}) - 1)) {
         my $alnchar = substr($aln1->{$seq}, $j, 1) ;
         if ($alnchar ne '-') {$cur_resno++;}
         if ($alnchar =~ /[a-z]/) {next;}
         $cur_alnpos++ ;

         if ($alnchar eq '-') {next;}

         $aln1_resno2alnpos->{$seq}->{$cur_resno} = $cur_alnpos ;
         $aln1_alnpos2resno->{$seq}->{$cur_alnpos} = $cur_resno ;
      }
   }

#   print STDERR " ALN1 parsing done for $target residues: ".
#      join(", ", sort {$a <=> $b} keys %{$aln1_resno2alnpos->{$target}})."\n";

#   print STDERR " ALN1 parsing done for $common residues: ".
#      join(", ", sort {$a <=> $b} keys %{$aln1_resno2alnpos->{$common}})."\n";

# Parse aln2 for seqresno_2_alnpos mapping for common
# ALIGNMENT 2: ASTRAL alignment string
   my $aln2_resno2alnpos = {}; my $aln2_alnpos2resno = {};
   my $cur_resno = 0 ;
   foreach my $alnpos (0 .. (length($aln2->{$common}) - 1)) {
      my $alnchar = substr($aln2->{$common}, $alnpos, 1) ;
      if ($alnchar eq '-') {next; }
      $cur_resno++ ;
      $aln2_resno2alnpos->{$common}->{$cur_resno} = $alnpos ;
      $aln2_alnpos2resno->{$common}->{$alnpos} = $cur_resno ;
   }
#   print STDERR " ALN2 parsing done for $common residues: ".
#      join(", ", sort {$a <=> $b} keys %{$aln2_resno2alnpos->{$common}})."\n";

# Combine seqresno_2_alnpos mappings to get target_2_alnpos2 mapping
   foreach my $target_resno (sort {$a <=> $b}
                             keys %{$aln1_resno2alnpos->{$target}}) {
#      print STDERR "\ntarget_resno $target_resno " ;

      my $aln1_pos = $aln1_resno2alnpos->{$target}->{$target_resno} ;
      if (!exists $aln1_alnpos2resno->{$common}->{$aln1_pos}) {next;}
#      print STDERR "   = alnpos_1 $aln1_pos" ;
      my $common_resno = $aln1_alnpos2resno->{$common}->{$aln1_pos} ;
#      print STDERR "   = common_resno $common_resno" ;

      if (!exists $aln2_resno2alnpos->{$common}->{$common_resno}) {next;}
      my $aln2_pos = $aln2_resno2alnpos->{$common}->{$common_resno} ;

#      print STDERR "   = alnpos_2 $aln2_pos" ;
      $aln2_resno2alnpos->{$target}->{$target_resno} = $aln2_pos ;
      $aln2_alnpos2resno->{$target}->{$aln2_pos} = $target_resno ;
   }
#   print STDERR "\n"; 

   return {
      resno2alnpos => $aln2_resno2alnpos->{$target},
      alnpos2resno => $aln2_alnpos2resno->{$target},
   } ;

}


=head2 readin_scop_cla

   Title:       readin_scop_cla()
   Function:    Reads in SCOP cla file (parsing logic from pibase::SCOP)
   Args:        ->{fn} = SCOP cla filename
   Returns:     ->{px2scopid}->{px_id} = scopid
                ->{scopid2class}->{scopid} = class
                ->{faid2sfid}->{fa_id} = sf_id0

=cut

sub readin_scop_cla {

   my $in = shift ;
   my $scopinfo ;

   open(SCOPCLA, $in->{fn}) ;
   while (my $line = <SCOPCLA>) { #parsing logic from pibase::SCOP
      if ($line =~ /^\#/) {next;}
      chomp $line;
      my ($sid, $pdb_id, $raw_domaindef, $sccs, $raw_px_id, $raw_sunid ) =
         split(/\t/, $line) ;

      my ($cl_id, $cf_id, $sf_id, $fa_id, $dm_id, $sp_id, $px_id) =
( $raw_sunid =~ /cl=(.+),cf=(.+),sf=(.+),fa=(.+),dm=(.+),sp=(.+),px=(.+)/ ) ;

      $scopinfo->{px2scopid}->{$px_id} = $sid ;
      $scopinfo->{scopid2class}->{$sid} = $sccs ;
      $scopinfo->{faid2sfid}->{$fa_id} = $sf_id ;
   }
   close(SCOPCLA) ;

   return $scopinfo ;
}


=head2 get_SUPFAM_selfhit

   Title:       get_SUPFAM_selfhit()
   Function:    Retrieves alignment string from SUPFAM self-hit files for
                  a particular SCOP template domain

   Args:        ->{supfam_specs} = configuration parameters
                ->{px_id} = SCOP px id
                ->{model_id} = SUPERFAMILY model id

   Returns:     self-hit alignment string

=cut

sub get_SUPFAM_selfhit {

   my $in = shift ;
   my $supfam_specs = $in->{supfam_specs} ;

   open(SELFHITF, $supfam_specs->{selfhits_dir}.'/'.$in->{sf_id}.'.tab') ;
   my @selfhit_headers ;
   my $alnstring ;
   while (my $line = <SELFHITF>) {
      chomp $line;
      my @t = split(/\t/, $line) ;
      if ($t[$supfam_specs->{selfhits_f2i}->{px_id}] eq $in->{px_id} &&
          $t[$supfam_specs->{selfhits_f2i}->{model_id}] eq $in->{model_id}) {
         $alnstring = $t[$supfam_specs->{selfhits_f2i}->{alnstring}] ;
         last;
      }
   }

   return $alnstring ;

}


=head2 summarize_results()

   Title:       summarize_results()
   Function:    Parses run_pilig_supfam_annotate() output and reports
                  summary of annotations.

   Args:        ->{results_fn} = run_pilig_supfam_annotate() output file
   Returns:     nothing
   Displays:    table describing numbers of proteins,domains,families,residues
                 with each kind of annotation

=cut

sub summarize_results {

   my $in = shift ;
   open(RESF, $in->{results_fn}) ;
   my $f2i ;
   my $header_line ;
   {
      my $headers = <RESF> ; $header_line = $headers ;
      $headers =~ s/^\#// ;
      my @headers = split(/\t/, $headers) ;
      map {$f2i->{$headers[$_]} = $_} (0 .. $#headers) ;
   }

   my $protein2dom ;
   my $protein2dom2bs_res ;
   while (my $line = <RESF>) {
      chomp $line;
      if ($line eq $header_line) { #skip if an (internal) header line
         next;
      } elsif ($line =~ /^#/) { #domain was skipped

         $line =~ s/^\#// ;
         my @t = split(/\t/, $line) ;
         $protein2dom->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$t[$f2i->{'classtype'}]} = $t[$f2i->{'class'}] ;

      } else { #annotation line

         my @t = split(/\t/, $line) ;
         $protein2dom->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$t[$f2i->{'classtype'}]} = $t[$f2i->{'class'}] ;

         my @res = split(/\, /, $t[$f2i->{'residues'}]) ;

         map {$protein2dom2bs_res->{$t[$f2i->{'classtype'}]}->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$t[$f2i->{'bs_type'}]}->{$_}++ } @res ;
      }
   }
   close(RESF) ;

   my $results ;
   my @classtypes = keys %{$protein2dom2bs_res} ;
   $results->{all}->{num_proteins} = keys %{$protein2dom} ;
   $results->{all}->{num_domains} = 0 ;
   $results->{all}->{num_residues} = 0 ;

# Iterate over all proteins
   foreach my $seq_id (keys %{$protein2dom}) {

      $results->{all}->{num_domains} += keys %{$protein2dom->{$seq_id}} ;
      my $res_list ; #tally of binding site annotated residues

# Iterate over all domains
      foreach my $res_range (keys %{$protein2dom->{$seq_id}}) {

# Count all residues in domain
         foreach my $subrange (split(',', $res_range)) {
            if ($subrange =~ /\-/) {
               my ($start, $stop) = ($subrange =~ /^([0-9]+)\-([0-9]+)$/);
               map {$res_list->{all}->{$_}++} ($start .. $stop) ;
            } else {
               $res_list->{all}->{$subrange}++ ;
            }
         }

# Count family/superfamily type
         foreach my $classtype (keys %{$protein2dom->{$seq_id}->{$res_range}}){
            $results->{all}->{'list_'.$classtype}->{$protein2dom->{$seq_id}->{$res_range}->{$classtype}}++ ;
         }

# Count up annotated binding sites: protein, domains, fam/sf, residues
         foreach my $classtype (keys %{$protein2dom2bs_res}) {
            if (!exists $protein2dom2bs_res->{$classtype}->{$seq_id} ||
             !exists $protein2dom2bs_res->{$classtype}->{$seq_id}->{$res_range}){
               next;}
            foreach my $bs_type (keys %{$protein2dom2bs_res->{$classtype}->{$seq_id}->{$res_range}}) {

               map {$res_list->{bs}->{$classtype}->{annotated}->{$_}++ } (keys
   %{$protein2dom2bs_res->{$classtype}->{$seq_id}->{$res_range}->{$bs_type}}) ;

               map {$res_list->{bs}->{$classtype}->{$bs_type}->{$_}++ } (keys
   %{$protein2dom2bs_res->{$classtype}->{$seq_id}->{$res_range}->{$bs_type}}) ;

               $results->{bs}->{$classtype}->{$bs_type}->{'list_proteins'}->{$seq_id}++;
               $results->{bs}->{$classtype}->{$bs_type}->{'list_domains'}->{$seq_id."\t".$res_range}++ ;
               $results->{bs}->{$classtype}->{$bs_type}->{"list_$classtype"}->{$protein2dom->{$seq_id}->{$res_range}->{$classtype}}++ ;

               $results->{bs}->{$classtype}->{annotated}->{'list_proteins'}->{$seq_id}++;
               $results->{bs}->{$classtype}->{annotated}->{'list_domains'}->{$seq_id."\t".$res_range}++ ;
               $results->{bs}->{$classtype}->{annotated}->{"list_$classtype"}->{$protein2dom->{$seq_id}->{$res_range}->{$classtype}}++ ;
            }
         }

      }

      foreach my $classtype (keys %{$res_list->{bs}}) {
# count up L and (P or p) residues (can do protein/domain/family check post-hoc)
         if (exists $res_list->{bs}->{$classtype}->{L} &&
             (exists $res_list->{bs}->{$classtype}->{P} ||
               exists $res_list->{bs}->{$classtype}->{p})) {
         foreach my $res (keys %{$res_list->{bs}->{$classtype}->{L}}) {
            if ((exists $res_list->{bs}->{$classtype}->{P} &&
                 exists $res_list->{bs}->{$classtype}->{P}->{$res}) ||
                (exists $res_list->{bs}->{$classtype}->{p} &&
                 exists $res_list->{bs}->{$classtype}->{P}->{$res})){

               $res_list->{bs}->{$classtype}->{"L_and_P"}->{$res}++ ;
            }
         }
         }

# count up residues per binding site type
         foreach my $bs_type (keys %{$res_list->{bs}->{$classtype}}) {
            $results->{bs}->{$classtype}->{$bs_type}->{num_residues} +=
               keys %{$res_list->{bs}->{$classtype}->{$bs_type}} ;
         }
      }

      $results->{all}->{num_residues} += keys %{$res_list->{all}} ;
   }

# Count final summary statistics
   foreach my $classtype (keys %{$results->{bs}}) {

#Count up L_and_P domains
      foreach my $score ('domains', $classtype, 'proteins') {
      foreach my $val (keys
         %{$results->{bs}->{$classtype}->{L}->{"list_$score"}}) {
        if (!exists $results->{bs}->{$classtype}->{P}->{"list_$score"}->{$val}||
            !exists $results->{bs}->{$classtype}->{p}->{"list_$score"}->{$val}){
            next; }
         $results->{bs}->{$classtype}->{"L_and_P"}->{"list_$score"}->{$val}++ ;
      }
      }

# count up number of families/superfamilies/domains/proteins
      foreach my $bstype (keys %{$results->{bs}->{$classtype}}) {
         $results->{bs}->{$classtype}->{$bstype}->{num_proteins} =
            keys %{$results->{bs}->{$classtype}->{$bstype}->{list_proteins}} ;

         $results->{bs}->{$classtype}->{$bstype}->{num_domains} =
            keys %{$results->{bs}->{$classtype}->{$bstype}->{list_domains}} ;

         $results->{bs}->{$classtype}->{$bstype}->{'num_'.$classtype} =
           keys %{$results->{bs}->{$classtype}->{$bstype}->{'list_'.$classtype}};
      }
   }

# Display overall counts
   foreach my $field  (keys %{$results->{all}}) {
      if ($field =~ /list/) {
         my $newfield = $field ;$newfield =~ s/^list_/num_/ ;
         my $val = keys %{$results->{all}->{$field}} ;
         print "ALL $newfield: $val\n" ;
      } else {
         print "ALL $field: ".$results->{all}->{$field}."\n";
      }
   }

# Display counts for each kind of binding site
   foreach my $class_type (keys %{$results->{bs}}) {
      foreach my $bs_type (keys %{$results->{bs}->{$class_type}}) {
         foreach my $field (keys %{$results->{bs}->{$class_type}->{$bs_type}}){
            if ($field !~ /^num/) {next;}
         print "BS: $class_type $bs_type $field ".
            $results->{bs}->{$class_type}->{$bs_type}->{$field}."\n"; }}}

# Format summary numbers as a LaTeX table
   print "\n\n" ;
   print join(" & ", 'CLASSTYPE', 'group', 'num_proteins',
               'num_domains',
               'num_fam',
               'num_residues')."\\\\\n" ;
   print join(" & ", '', 'annotated_domains',
               $results->{all}->{num_proteins},
               $results->{all}->{num_domains},
               $results->{all}->{num_fam},
               $results->{all}->{num_residues})."\\\\\n";
   foreach my $class_type (keys %{$results->{bs}}) {
      foreach my $bs_type (keys %{$results->{bs}->{$class_type}}) {
         print join(" & ", $class_type, $bs_type,
         $results->{bs}->{$class_type}->{$bs_type}->{num_proteins},
         $results->{bs}->{$class_type}->{$bs_type}->{num_domains},
         $results->{bs}->{$class_type}->{$bs_type}->{num_fam},
         $results->{bs}->{$class_type}->{$bs_type}->{num_residues},
         )."\\\\\n" ;
      }
   }

}
