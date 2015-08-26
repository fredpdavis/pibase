=head1 NAME

pibase::SGE - perl module for SGE cluster interaction

=head1 DESCRIPTION

Perl module with routines to interact with an SGE cluster,
adapted from routines in modtie.pm

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

package pibase::SGE ;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/_clust_qsub _clust_qstat _clust_split_ins _clust_merge_outs/ ;

use Cwd qw/getcwd/ ;
use File::Temp qw/tempfile tempdir/ ;
use File::Path qw/mkpath/ ;


=head2 _clust_split_ins()

   Title:       _clust_split_ins()
   Function:    Splits an input file for a parallel cluster run.
   Usage:
         my $split_dir = tempdir("splits_assign_model_domains.XXXXX") ;
         my $splits = _clust_split_ins({
            fn => $temp_fn->{model_domains_in},
            dir => $split_dir,
            numjobs => $in->{cluster}->{numjobs}
         });

   gives back:  ->{numjobs} = number of tasks
                ->{tasklist} = space delimited list of split inputs

=cut

sub _clust_split_ins {

   use File::Basename qw/basename/ ;

   my $in = shift ;

   if (!exists $in->{fn}) {
      die "_clust_split_ins: input file not specified\n" ; }

   if (!exists $in->{dir}) {
      $in->{dir} = "./" ; }

   my $inbase = basename($in->{fn}) ;

   my $header_lines = '';
   my $num_lines = 0 ;
   open(INF, $in->{fn}) ;
   while (my $line = <INF>) {
      if ($line =~ /^#SET /) {
         $header_lines .= $line ;
      } else {
         $num_lines++ ;
      }
   }
   close(INF);

   my $splitlines = POSIX::ceil($num_lines / $in->{numjobs});

   if (exists $in->{minlines} && ($splitlines < $in->{minlines})) {
      $splitlines = $in->{minlines} ;
   }

   my @splitfiles ;

   $num_lines = 0 ;

   my $tasklist = '';
   my $cur_splitnum = 1 ;

   my $cur_fn = "split.$inbase.$cur_splitnum" ;
   push @splitfiles, $cur_fn ;
   $tasklist .= "'$cur_fn' " ;

   open(INF, $in->{fn}) ;
   open(OUTF, ">$in->{dir}/$cur_fn") ;
   print OUTF $header_lines ;
   while (my $line = <INF>) {
      if ($line =~ /^#SET/) {next;}
      if ($num_lines > $splitlines) {
         close(OUTF) ;
         $cur_splitnum++ ;
         $cur_fn = "split.$inbase.$cur_splitnum" ;
         push @splitfiles, $cur_fn ;
         $tasklist .= "'$cur_fn' " ;
         open(OUTF, ">$in->{dir}/$cur_fn") ;
         print OUTF $header_lines ;
         $num_lines = 0 ;
      }
      print OUTF $line ;
      $num_lines++ ;
   }
   close(OUTF) ;
   close(INF);
   $tasklist =~ s/ $// ;
#   print STDERR "OUT TASKLIST: $tasklist\n" ;


   return {
      numjobs => $cur_splitnum,
      splitfiles => \@splitfiles,
      tasklist => $tasklist
   } ;
}


=head2 _clust_merge_outs()

   Title:       _clust_merge_outs()
   Function:    Merge the output files from a cluster run
   Usage:
         _clust_merge_outs({
            script_fn => $sgescript_fn,
            out_fn => $temp_fn->{cutpaths_out},
            err_fn => $temp_fn->{cutpaths_err},
            job_id => $qsub_job_id,
            outdir => $sge_outdir,
            numjobs => $splits->{numjobs}
         }) ;

=cut

sub _clust_merge_outs {

   my $in = shift ;
   open(OUTFN, ">$in->{out_fn}") ;
   open(ERRFN, ">$in->{err_fn}") ;
   foreach my $j ( 1 .. $in->{numjobs}) {
      if (-s "$in->{outdir}/$in->{script_fn}.o$in->{job_id}.$j") {
         open(SGEOUT, "$in->{outdir}/$in->{script_fn}.o$in->{job_id}.$j") ;
         while (my $line = <SGEOUT>) {
            if ($line !~ /^Warning: no access to tty/ &&
                $line !~ /^Thus no job control/ &&
                $line !~ /^#sgejob run/) {
               print OUTFN $line ;
            }
         }
         close(SGEOUT) ;
      }

      if (-s "$in->{outdir}/$in->{script_fn}.e$in->{job_id}.$j") {
         open(SGEERR, "$in->{outdir}/$in->{script_fn}.e$in->{job_id}.$j") ;
         while (my $line = <SGEERR>) {
            print ERRFN $line ; }
         close(SGEERR) ;
      }
   }
   close(OUTFN) ;
   close(ERRFN) ;

   return ;
}


=head2 _clust_qsub()

   Title:       _clust_qsub()
   Function:    Submit job to the SGE queue
   Returns:     $_ = SGE job id
   Usage:
         _clust_qsub({
            sgescript_fn => $sgescript_fn,
         }) ;


         while (1) {
            sleep $modtie_specs->{cluster}->{qstat_sleep} ;
            my $job_status = _clust_qstat({job_id => $qsub_job_id}) ;
            if ($job_status) {last;}
         }

=cut

sub _clust_qsub {

   my $in = shift;

   my $cwd = getcwd() ;
   my $status = `ssh login-eddy \"cd $cwd; qsub $in->{sgescript_fn}\"` ;
   chomp $status;
   if ($status !~ /has been submitted/) {
      die "failed to qsub $in->{sgescript_fn}";}
#   my ($job_id) = ($status =~ /Your job ([0-9]+)\./) ;
# sometimes its Your job 12345 ("job.SGE.sh") has been submitted.
#others it is Your job 30614.1-31:1 ("job.SGE.sh") has been submitted.
# don't know why one or the other. keep it simple

# Janelia SGE output:
# Your job-array 9967277.1-100:1 ("pb.call_residue_info.1RrIi.SGE.sh") has been submitted
   my ($job_id) = ($status =~ /Your job-array ([0-9]+)/) ;
   print STDERR "submitted $in->{sgescript_fn} job $job_id\n" ;

   return $job_id ;
}


=head2 _clust_qstat()

   Title:       _clust_qstat()
   Function:    Check if a job_id is still in SGE queue
   Returns:     1 = finished, 0 = still going.
   Usage:
            my $job_status = _clust_qstat({job_id => $qsub_job_id}) ;

=cut

sub _clust_qstat {

   my $in = shift;
   system("ssh login-eddy qstat -j $in->{job_id} > /tmp/qstat.$$.out 2>&1") ;
   my $status = `cat /tmp/qstat.$$.out` ;
   system("rm /tmp/qstat.$$.out") ;
   if ($status =~ /Following jobs do not exist/) {
      return 1;
   } else {
      return 0;
   }

}

1 ;
