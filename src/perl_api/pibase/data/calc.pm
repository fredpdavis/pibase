=head1 NAME

pibase::data::calc - perl module for pibase data calculation routines

=head1 DESCRIPTION

Perl package that provides pibase data calculation routines

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

package pibase::data::calc ;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT = qw/bdp_path_2_id call_residue_info call_chain_info bdp_subset_translator calc_bdp_secstrx calc_subsets_sasa calc_interface_dsasa calc_interface_secstrx calc_interface_bs_secstrx_basic calc_subsets_sequence calc_interface_secstrx_contacts calc_interface_secstrx_profile calc_interface_resvector calc_interface_sse_topology calc_bdp_interaction_topology calc_interface_size calc_bdp_interaction_topology_graph/ ;

use pibase qw/connect_pibase get_specs mysql_hashload mysql_fetchcols mysql_hasharrload safe_move sid_2_domdir/;
use pibase::interatomic_contacts qw/raw_contacts_select/ ;
use pibase::modeller qw/calc_sasa/ ;
use pibase::PDB::chains qw/chain_info/;
use pibase::PDB::residues qw/residue_info/;
use pibase::PDB::sec_strx qw/run_dssp parse_dssp/;
use pibase::residue_math qw/residue_int residue_add residue_inrange/;
use pibase::calc::interfaces qw/_interface_detect_calc__calc_res_pairs/ ;
use pibase::specs ;
use pibase::benchmark qw/memusage/;
use File::Temp qw/tempfile tempdir/ ;
use File::Path qw/mkpath/ ;
use POSIX qw/ceil floor/ ;
use Sys::Hostname ;
use pibase::SGE ;

use constant FLAG_NOSUBSRES => 0 ;

=head2 call_residue_info()

   Name:        call_residue_info()
   OLD NAME:    bdp_resinfo_caller_calc()
   Function:    Routine to calculate residue_info data for list of bdp_ids
                  (calls residue_info())
   Args:        none
   Returns:     nothing
   STDIN:       tabbed: bdp_id, bdp_path
   STDOUT:      tabbed: bdp_id, "\N", residue info file
   Files out:   foreach bdp_id:
                o <$pibase_specs->{metatod_dir}->{bdp_residues}>/
                  bdp_residues_<bdp_id>.<hostname>.<XXXXX>.resinfo.out.gz
                     (format per pibase::PDB::residue_info())

=cut

sub call_residue_info {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }


# read in bdp_id, bdp_path info from input if specified, otherwise
#  get list from bdp_files, and split locally.

   my $bdpid2path = {};
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line =~ /^#/) {next;}
         chomp $line;
         my ($bdp_id, $bdp_path) = split(/\t/, $line) ;
         $bdpid2path->{$bdp_id} = $bdp_path ;
      }
      close(INF) ;
   } else {
      my ($dbh) = pibase::connect_pibase() ;
      $bdpid2path = pibase::mysql_hashload($dbh,
         "SELECT bdp_id, file_path FROM bdp_files") ;
   }

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
# if this is master script, split input, recall self with
# in_fn specified, cluster_fl = 0, import_fl = 0
# send out, cluster run, and return merged

      print "* call_residue_info() ".localtime() if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{call_residue_info_in}, $temp_fn->{call_residue_info_in}) =
         tempfile("splits_call_residue_info_input.XXXXX");
      ($temp_fh->{call_residue_info_out}, $temp_fn->{call_residue_info_out}) =
         tempfile("splits_call_residue_info_SGEout_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{call_residue_info_out}) ;
      ($temp_fh->{call_residue_info_err}, $temp_fn->{call_residue_info_err}) =
         tempfile("splits_call_residue_info_SGEerr_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{call_residue_info_err}) ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdpid2path}) {
         print {$temp_fh->{call_residue_info_in}}
            join("\t", $bdp_id, $bdpid2path->{$bdp_id})."\n" ; }
      close($temp_fh->{call_residue_info_in}) ;

      my $split_dir = tempdir("splits_call_residue_info.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{call_residue_info_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.call_residue_info.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/call_residue_info/ ;

main() ;

sub main {

         pibase::data::calc::call_residue_info({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.call_residue_info.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.call_residue_info.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id}) ;
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{call_residue_info_out},
         err_fn => $temp_fn->{call_residue_info_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{call_residue_info_out},
           $temp_fn->{call_residue_info_out}) ;
      open(REALOUTF,">".$pibase_specs->{buildfiles}->{bdp_residues_tables}) ;
      while (my $line = readline($temp_fh->{call_residue_info_out})) {
         if ($line =~ /^\#/) {next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{call_residue_info_out}) ;
      close(REALOUTF) ;

   } else {

      my $movethreshold = 100 ;
      my $temp_dir = tempdir(CLEANUP => 1) ;
      chdir $temp_dir ;
      my (@movethese, @moveto) ;

# if actual compute node, then just run it and print to stdout

      foreach my $bdp_id (keys %{$bdpid2path}) {
# Set output file for STDOUT and STDERR from pdb_resinfo.pl.
         my $resf = "bdp_residues_$bdp_id.resinfo.out" ;

         my $res_depositdir = $pibase_specs->{metatod_dir}->{bdp_residues}.'/'.
            POSIX::floor($bdp_id / 1000) ;
         if (!-s $res_depositdir) { mkpath $res_depositdir; }

# get residue information
         my $bdp_path = $bdpid2path->{$bdp_id} ;
#         print STDERR "calling residue_info on $bdp_path\n";
         pibase::PDB::residues::residue_info({
            pdb_fn => $bdp_path,
            outfile => $resf,
            identifier => $bdp_id,
            gzip => 1
         }) ;
         $resf .= '.gz' ;
         push @movethese, $resf ;
         push @moveto, $res_depositdir ;

# Initialize the error flag.
         my $err = 0 ;

# If the resinfo file is empty, set the error flag.
         if (!-s $resf) {
            $err = 1;
            print STDERR "ERROR $bdp_path: residue_info() error: empty $resf\n";
            unlink $resf ;
         }

# If the error flag is not set, display resinfo file location to STDOUT.
         if (!($err)) {
            my @outvals = ($bdp_id, '\N', $res_depositdir.'/'.$resf);
            print join("\t", @outvals)."\n" ;
         }

         if ($#movethese == $movethreshold) {
            foreach my $j (0 .. $#movethese) {
               my $t_fn = $movethese[$j] ;
               my $t_dir = $moveto[$j] ;
               pibase::safe_move($t_fn, $t_dir);}
            @movethese = () ; @moveto = () ;
         }
      }

      foreach my $j (0 .. $#movethese) {
         my $t_fn = $movethese[$j] ;
         my $t_dir = $moveto[$j] ;
         pibase::safe_move($t_fn, $t_dir);
      }
   }

   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_residues_tables file (specified in specs) into pibase
      $import_status = pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{bdp_residues_tables}}) ;
   }

   return $import_status ;

}


sub OLD_bdp_resinfo_caller_calc {

   my $specs = pibase::get_specs() ;
   my $res_depositdir = $specs->{metatod_dir}->{bdp_residues} ;
   my $movethreshold = 100 ;

   my $temp_dir = tempdir(CLEANUP => 1) ;
   chdir $temp_dir ;
   my @movethese ;

   while (my $line = <STDIN>) {
      if ($line =~ /^#/) {next;}

      chomp $line ;
      my ($bdp_id, $bdp_path) = split(/\t/, $line) ;

# Create temporary files for STDOUT and STDERR from pdb_resinfo.pl.
      my $hostname = hostname() ;
      my ($fh, $resf) = tempfile("bdp_residues_$bdp_id.$hostname.XXXX",
                            SUFFIX => '.resinfo.out') ; close($fh) ;

# get residue information
      residue_info({
         pdb_fn => $bdp_path,
         outfile => $resf,
         identifier => $bdp_id,
         gzip => 1
      }) ;
      $resf .= '.gz' ;

      push @movethese, $resf ;

# Initialize the error flag.
      my $err = 0 ;

# If the resinfo file is empty, set the error flag.
      if (!-s $resf) {
         $err = 1;
         print STDERR "ERROR $bdp_path: resinfo calculate error: empty $resf\n" ;
         unlink $resf ;
      }

# If the error flag is not set, display resinfo file location to STDOUT.
      if (!($err)) {
         my @outvals = ($bdp_id, '\N', $resf) ;
         print join("\t", @outvals)."\n" ;
      }

      if ($#movethese == $movethreshold) {
         foreach my $t_fn (@movethese) {
            pibase::safe_move($t_fn,$res_depositdir);}
         @movethese = () ;
      }
   }

   foreach my $t_fn (@movethese) {
      pibase::safe_move($t_fn,$res_depositdir);
   }

}

=head2 calc_bdp_secstrx()

   Function:    Routine to calculate secondary structure information for a list
                  of bdp_ids (calls pibase::PDB::sec_strx::run_dssp() and
                  parse_dssp())
   Args:        none
   Returns:     nothing
   STDIN:       tabbed: bdp_id, bdp_path
   STDOUT:      tabbed: bdp_id, sec_strx_fn
   Files out:   foreach bdp_id:
                o <$pibase_specs->{metatod_dir}->{bdp_secstrx}>/
                  secstrx.<bdp_id>.<hostname>.<XXXXX>.txt
                     (format per pibase::PDB::residue_info())
                     1. bdp_id
                     2. chain
                     3. resno
                     4. detailed sec strx (H, G, I, B, E, T, S, ' ')
                     5. basic sec strx (H, B, T, ' ')
                     6. sec strx element number

=cut

sub calc_bdp_secstrx {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }

# read in bdp_id, bdp_path info from input if specified, otherwise
#  get list from bdp_files, and split locally.

   my $bdpid2path = {};
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line =~ /^#/) {next;}
         chomp $line;
         my ($bdp_id, $bdp_path) = split(/\t/, $line) ;
         $bdpid2path->{$bdp_id} = $bdp_path ;
      }
      close(INF) ;
   } else {
      my ($dbh) = pibase::connect_pibase() ;
      $bdpid2path = pibase::mysql_hashload($dbh,
         "SELECT bdp_id, file_path FROM bdp_files") ;
   }

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
# if this is master script, split input, recall self with
# in_fn specified, cluster_fl = 0, import_fl = 0
# send out, cluster run, and return merged

      print "* calc_bdp_secstrx() ".localtime() if (!exists $in->{quiet_fl});
      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{calc_bdp_secstrx_in}, $temp_fn->{calc_bdp_secstrx_in}) =
         tempfile("splits_calc_bdp_secstrx_input.XXXXX");
      ($temp_fh->{calc_bdp_secstrx_out}, $temp_fn->{calc_bdp_secstrx_out}) =
         tempfile("splits_calc_bdp_secstrx_SGEout_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{calc_bdp_secstrx_out}) ;
      ($temp_fh->{calc_bdp_secstrx_err}, $temp_fn->{calc_bdp_secstrx_err}) =
         tempfile("splits_calc_bdp_secstrx_SGEerr_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{calc_bdp_secstrx_err}) ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdpid2path}) {
         print {$temp_fh->{calc_bdp_secstrx_in}}
            join("\t", $bdp_id, $bdpid2path->{$bdp_id})."\n" ; }
      close($temp_fh->{calc_bdp_secstrx_in}) ;

      my $split_dir = tempdir("splits_calc_bdp_secstrx.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{calc_bdp_secstrx_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.calc_bdp_secstrx.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/calc_bdp_secstrx/ ;

main() ;

sub main {

         pibase::data::calc::calc_bdp_secstrx({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.calc_bdp_secstrx.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.calc_bdp_secstrx.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id}) ;
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{calc_bdp_secstrx_out},
         err_fn => $temp_fn->{calc_bdp_secstrx_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{calc_bdp_secstrx_out},
           $temp_fn->{calc_bdp_secstrx_out}) ;
      open(REALOUTF,">".$pibase_specs->{buildfiles}->{bdp_secstrx_tables}) ;
      while (my $line = readline($temp_fh->{calc_bdp_secstrx_out})) {
         if ($line =~ /^\#/) {next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{calc_bdp_secstrx_out}) ;
      close(REALOUTF) ;

   } else {

      my $movethreshold = 100 ;
      my $temp_dir = tempdir(CLEANUP => 1) ;
      chdir $temp_dir ;
      my (@movethese, @moveto) ;

      foreach my $bdp_id (keys %{$bdpid2path}) {
         my $bdp_path = $bdpid2path->{$bdp_id} ;
   
# Run DSSPget residue information
         my $dssp_run ;
         ($dssp_run->{status}, $dssp_run->{res}) =
            pibase::PDB::sec_strx::run_dssp({pdb_fn=>$bdp_path }) ;
   
         if ($dssp_run->{status} == 0) {
            print STDERR "ERROR $bdp_path: DSSP run error: empty ".
                         " $dssp_run->{res}->{out}\n" ;
            unlink $dssp_run->{res}->{out}, $dssp_run->{res}->{err} ;
            next ;
         }
   
# Parse dssp
         my $secstrx = pibase::PDB::sec_strx::parse_dssp(
            $dssp_run->{res}->{out}) ;
   
         if ($#{$secstrx->{ordering}} < 0) {
            print STDERR "ERROR $bdp_path: DSSP output parse error\n" ;
            print STDERR "   $dssp_run->{res}->{out}, $dssp_run->{res}->{err}\n" ;
            unlink $dssp_run->{res}->{out}, $dssp_run->{res}->{err} ;
            next ;
         }
   
         unlink $dssp_run->{res}->{out}, $dssp_run->{res}->{err} ;
   
# Create temporary files for STDOUT and STDERR from pdb_resinfo.pl.
         my $secstrxf = "bdp_secstrx_$bdp_id.out" ;
         my $secstrx_depositdir = $pibase_specs->{metatod_dir}->{bdp_secstrx}.
            '/'.POSIX::floor($bdp_id / 1000) ;
         if (!-s $secstrx_depositdir) { mkpath $secstrx_depositdir; }

         open(SECSTRXF, ">$secstrxf") ;
         foreach my $j ( 0 .. $#{$secstrx->{ordering}}) {
            my $ressig = $secstrx->{ordering}->[$j] ;
            my ($resno, $chain) = split(/\n/, $ressig) ;
            my @outvals = ( $bdp_id, $chain, $resno,
                              $secstrx->{detail}->{$ressig},
                              $secstrx->{basic}->{$ressig},
                              $secstrx->{ssnum}->{$ressig},
                              $secstrx->{ssnum_basic}->{$ressig});
            print SECSTRXF join("\t", @outvals)."\n" ;
         }
         close(SECSTRXF) ;
         system("gzip $secstrxf") ; $secstrxf .= '.gz' ;
         push @movethese, $secstrxf ;
         push @moveto, $secstrx_depositdir ;
   
         my @outvals = ($bdp_id, '\N', $secstrx_depositdir.'/'.$secstrxf) ;
         print join("\t", @outvals)."\n" ;

         if ($#movethese == $movethreshold) {
            foreach my $j (0 .. $#movethese) {
               my $t_fn = $movethese[$j] ;
               my $t_dir = $moveto[$j] ;
               if (! -d $t_dir) { mkpath($t_dir) ;}
               pibase::safe_move($t_fn, $t_dir);}
            @movethese = () ; @moveto = () ;
         }
      }

      foreach my $j (0 .. $#movethese) {
         my $t_fn = $movethese[$j] ;
         my $t_dir = $moveto[$j] ;
         if (! -d $t_dir) { mkpath($t_dir) ;}

         pibase::safe_move($t_fn, $t_dir);
      }
      @movethese = () ; @moveto = () ;

   }

   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_secstrx_tables file (specified in specs) into pibase
      $import_status = pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{bdp_secstrx_tables}}) ;
   }

   return $import_status ;

}


=head2 bdp_subset_translator()

   Function:    Translates the raw domain definitions to the same framework for
                  PDB and PISA entries
   Args:        none
   Returns:     nothing
   STDIN:       bdp_path
   STDOUT:      tabbed: subsets_filename, subsets_details_filename
                tabbed: bdp_id, '\N', subset_residues file path
                  (if FLAG_NOSUBSRES is set to 0 - currently set to 1)

   Files out:   subsets.<hostname>.<XXX>.<pibase db name>
                  1. subset_id
                  2. bdp_id
                  3. '\N'
                  4. subset_source_id
                  5. domain classification

                subsets_details.<hostname>.<XXX>.<pibase db name>
                  1. subset_id
                  2. domain segment number
                  3. serial chain number
                  4. chain identifer
                  5. start residue number
                  6. integer part of start residue number
                  7. end residue number
                  8. integer part of end residue number

                foreach bdp_id (if FLAG_NOSUBSRES is 0; currently set to 1):
                  <$pibase_specs->{metatod_dir}->{subsets_residues}>/
                  subsets_residues.<bdp_id>.<XXXX>.<pibase db name>
                     (format per pibase::PDB::residue_info())
                     1. bdp_id
                     2. serial chain number
                     3. chain id
                     4. serial residue number
                     5. residue number
                     6. integer part of residue number
                     7. subset_id

=cut

sub bdp_subset_translator {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }
   my $pibase = pibase::dbname() ;
   my $movethreshold= 500 ;

   my ($bdp_2_pdb_id, $path_2_bdp_id, $pdb_id_2_bdp_id) =
      pibase::todload_bdp_ids(
      'bdp_id_2_pdb_id', 'path_2_bdp_id', 'pdb_id_2_bdp_id') ;

# Read in bdp_paths to be processed
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         chomp $line;
         my ($bdp_id, $bdp_path) = split(/\t/, $line) ;
         my $pdb_id = $bdp_2_pdb_id->{$bdp_id} ;
         $in->{bdps}->{$pdb_id}->{$bdp_id} = $bdp_path ;
      }
      close(INF) ;
   } else {
      foreach my $bdp_path (keys %{$path_2_bdp_id}) {
         my $bdp_id = $path_2_bdp_id->{$bdp_path} ;
         $in->{bdps}->{$bdp_2_pdb_id->{$bdp_id}}->{$bdp_id} = $bdp_path ;
      }
   }

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) { #master
# split bdp_ids and send off to compute - everything else is computing.

      print "* bdp_subset_translator() ".localtime() if (!exists
         $in->{quiet_fl});

# make sure subset_source_id exists for PDB_CHAINS entry
      {
         my ($dbh, undef) = pibase::connect_pibase() ;
         my ($t) = pibase::mysql_fetchcols($dbh,
            'SELECT * FROM subsets_source WHERE subset_source = "pdb_chains"') ;

         if ($#{$t} < 0) {
            my $cur_ver = pibase::timestamp() ;
            pibase::mysql_runcom($dbh,
               'INSERT INTO subsets_source '.
               '(subset_source_id, subset_source, version) '.
               'VALUES("0", "pdb_chains", "'.$cur_ver.'")'
            );

            pibase::data::access::copy_table_to_disk({
               tables => ['subsets_source'] }) ;
         }
      }

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{bdp_subset_translator_in},
       $temp_fn->{bdp_subset_translator_in}) =
         tempfile("splits_bdp_subset_translator_input.XXXXX");

      ($temp_fh->{bdp_subset_translator_out},
       $temp_fn->{bdp_subset_translator_out}) =
         tempfile("splits_bdp_subset_translator_SGEout_XXXXX",
         SUFFIX => '.pibase');
         close($temp_fh->{bdp_subset_translator_out}) ;

      ($temp_fh->{bdp_subset_translator_err},
       $temp_fn->{bdp_subset_translator_err}) =
         tempfile("splits_bdp_subset_translator_SGEerr_XXXXX",
         SUFFIX => '.pibase');
         close($temp_fh->{bdp_subset_translator_err}) ;

      foreach my $pdb_id (sort keys %{$in->{bdps}}) {
         foreach my $bdp_id (sort {$a <=> $b} keys %{$in->{bdps}->{$pdb_id}}) {
            print {$temp_fh->{bdp_subset_translator_in}}
               join("\t", $bdp_id, $in->{bdps}->{$pdb_id}->{$bdp_id})."\n"; } }
      close($temp_fh->{bdp_subset_translator_in}) ;

      my $cur_numjobs ;
      if (exists $in->{numtasks_cluster}) {
         $cur_numjobs = $in->{numtasks_cluster} ;
      } else {
         $cur_numjobs = $pibase_specs->{SGE}->{numjobs} ;
      }

      my $split_dir = tempdir("splits_bdp_subset_translator.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{bdp_subset_translator_in},
         dir => $split_dir,
         numjobs => $cur_numjobs,
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.bdp_subset_translator.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/bdp_subset_translator/ ;

main() ;

sub main {

         pibase::data::calc::bdp_subset_translator({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.bdp_subset_translator.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.bdp_subset_translator.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id});
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{bdp_subset_translator_out},
         err_fn => $temp_fn->{bdp_subset_translator_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{bdp_subset_translator_out},
           $temp_fn->{bdp_subset_translator_out}) ;
      open(REALOUTF_SUB,">".
         $pibase_specs->{buildfiles}->{subsets}) ;
      open(REALOUTF_SUBDET,">".
         $pibase_specs->{buildfiles}->{subsets_details}) ;
      open(REALOUTF_SUBREST,">".
         $pibase_specs->{buildfiles}->{subsets_residues_tables});
      while (my $line = readline($temp_fh->{bdp_subset_translator_out})) {
# need to grep out ^SUBSETS, ^SUBSETS_DETAILS, and ^SUBSETS_RESIDUES entries
         if ($line =~ /^\#/) {
            next;
         } elsif ($line =~ /^SUBSETS\t/) {
            $line =~ s/^SUBSETS\t// ;
            print REALOUTF_SUB $line ;
         } elsif ($line =~ /^SUBSETS_DETAILS\t/) {
            $line =~ s/^SUBSETS_DETAILS\t// ;
            print REALOUTF_SUBDET $line ;
         } elsif ($line =~ /^SUBSETS_RESIDUES_TABLES\t/) {
            $line =~ s/^SUBSETS_RESIDUES_TABLES\t// ;
            print REALOUTF_SUBREST $line ;
         }
      }
      close($temp_fh->{bdp_subset_translator_out}) ;
      close(REALOUTF_SUB) ;
      close(REALOUTF_SUBDET) ;
      close(REALOUTF_SUBREST) ;

      my $import_status ;
      if (exists $in->{import_fl} && $in->{import_fl} == 1) {
         foreach my $import_fn (qw/subsets subsets_details subsets_residues_tables/) {
            $import_status->{$import_fn} = pibase::mysqlimport({
               pibase_specs => $pibase_specs,
               fn => $pibase_specs->{buildfiles}->{$import_fn}}) ;
         }
      }
      return $import_status ;

   } else { #compute node

# Set names of tables that need to be generated.
      my $tables ;
      $tables->{subsets_residues}->{namepre} = 'subsets_residues' ;

# Load the entire bdp_chains array into memory and index with bdp_id.
      print STDERR "reading in bdp_chains: " ;
      my $allchains ;
      ( $allchains->{bdp_id},
        $allchains->{real_chain_no},
        $allchains->{real_chain_id},
        $allchains->{pdb_chain_no},
        $allchains->{pdb_chain_id},
        $allchains->{start_resno},
        $allchains->{end_resno} ) =
         pibase::rawselect_tod( "SELECT bdp_id, real_chain_no, ".
         "real_chain_id, pdb_chain_no, pdb_chain_id, start_resno, end_resno ".
         "FROM bdp_chains") ;
      print STDERR $#{$allchains->{bdp_id}}." chains\n" ;

      $allchains->{real_chain_id} =
         pibase::replace_undefs($allchains->{real_chain_id}, ' ') ;

      $allchains->{pdb_chain_id} =
         pibase::replace_undefs($allchains->{pdb_chain_id}, ' ') ;

      my $allchains_point ;
      my $allchains_pdb_chid_2_bdp_chid ;
      my $allchains_pdb_chid_2_bdp_chno ;
      foreach my $j (0 .. $#{$allchains->{bdp_id}}) {
         push @{$allchains_point->{$allchains->{bdp_id}->[$j]}}, $j ;
         push @{$allchains_pdb_chid_2_bdp_chid->{$allchains->{bdp_id}->[$j]}->{$allchains->{pdb_chain_id}->[$j]}},
            $allchains->{real_chain_id}->[$j] ;
         push @{$allchains_pdb_chid_2_bdp_chno->{$allchains->{bdp_id}->[$j]}->{$allchains->{pdb_chain_id}->[$j]}},
            $allchains->{real_chain_no}->[$j] ;
      }

# Load the entire subsets table into memory.

      print STDERR "reading in subsets: \n" ;
      my $allsubs ;
      ( $allsubs->{pdb_id},
        $allsubs->{subset_id},
        $allsubs->{description},
        $allsubs->{source_id},
        $allsubs->{class} ) =
         pibase::rawselect_tod("SELECT pdb_id, subset_id, description, ".
                               "subset_source_id, class FROM subsets") ;
      print STDERR "x\n" ;

      my $subset_source2id ;
      {
         my $a ;
         ( $a->{subset_source},
           $a->{subset_source_id}) =
         pibase::rawselect_tod("SELECT subset_source, subset_source_id ".
                               "FROM subsets_source") ;

         foreach my $j ( 0 .. $#{$a->{subset_source}}) {
            $subset_source2id->{$a->{subset_source}->[$j]} =
               $a->{subset_source_id}->[$j] ; }
      } 

      print STDERR $#{$allsubs->{pdb_id}}." subsets\n" ;

      my $allsubs_point ;
      foreach my $j (0 .. $#{$allsubs->{pdb_id}}) {
         push @{$allsubs_point->{$allsubs->{pdb_id}->[$j]}}, $j ; }

# Load the entire bdp_residues_tables into a hash.
      print STDERR "reading in bdp_residue_tables: ." ;
      my $all_bdp_residues_tables ;
      {
         my ($t_bdpid, $t_sourcefile) = pibase::rawselect_tod(
            "SELECT bdp_id, source_file FROM bdp_residues_tables") ;
         foreach my $j ( 0 .. $#{$t_bdpid}) {
            $all_bdp_residues_tables->{$t_bdpid->[$j]} = $t_sourcefile->[$j] ; }
      }
      print STDERR "x\n" ;

# Load the entire subsets_details table into memory.
      print STDERR "reading in subsets_details: " ;
      my $allsubsdet ;
      ( $allsubsdet->{subset_id},
        $allsubsdet->{chain_id},
        $allsubsdet->{start},
        $allsubsdet->{end} ) =
         pibase::rawselect_tod( "SELECT subset_id, chain_id, ".
         "start_resno, end_resno FROM subsets_details") ;
      print STDERR $#{$allsubsdet->{subset_id}}." subsets_details\n" ;

      $allsubsdet->{chain_id} =
         pibase::replace_undefs_blanks($allsubsdet->{chain_id}, ' ') ;

      my $allsubsdet_point ;
      foreach my $j (0 .. $#{$allsubsdet->{subset_id}}) {
         push @{$allsubsdet_point->{$allsubsdet->{subset_id}->[$j]}}, $j ; }

# Open output files for subsets and subsets_details entries.
      my $temp_dir = tempdir(CLEANUP=>1) ;
      chdir $temp_dir ;
      my (@movethese, @movehere) ;

# OUTPUT: Display the file names.
#   print "$subsets_f\t$subsdet_f\n" ;

#Iterate through the bdp_ids to be processed
      my $lastpdb = 'a';
#      print STDERR "now on:  " ;
      foreach my $pdb_id (keys %{$in->{bdps}}) {
#         print STDERR "\b"x(length($lastpdb)).$pdb_id ;

#ERRCATCH: If no subsets defined for this pdb_id, STDERR and abort BDP.
         if (! exists $allsubs_point->{$pdb_id}) {
            print STDERR "ERROR: pdb $pdb_id: no subsets found\n" ;
            next;
         }

#Load pdb chain number, id, start and end residue numbers for the current pdb_id
         my $raw_bdpid = $pdb_id_2_bdp_id->{$pdb_id} ;

         if (! defined $raw_bdpid) {
            print STDERR "ERROR: pdb $pdb_id: did not find bdp_id entry\n" ;
            next ;
         } elsif ( !exists $allchains_point->{$raw_bdpid} ) {
            print STDERR "ERROR: pdb $pdb_id: did not find chains entries\n";
            next ;
         }


         my ($pdb_chain_no, $pdb_chain_id, $pdb_start, $pdb_end) ;
         {
            my @tind = @{$allchains_point->{$raw_bdpid}} ;
            push @{$pdb_chain_no}, @{$allchains->{real_chain_no}}[@tind] ;
            push @{$pdb_chain_id}, @{$allchains->{real_chain_id}}[@tind] ;
            push @{$pdb_start}, @{$allchains->{start_resno}}[@tind] ;
            push @{$pdb_end}, @{$allchains->{end_resno}}[@tind] ;
         }

#ERRCATCH: If no PDB chains loaded, STDERR and abort current PDB.
         if ($#{$pdb_chain_no} < 0) {
            print STDERR "ERROR: PDB $pdb_id: no pdb_chains loaded\n" ;
            next ;
         }

#Declare hashes for pdb_chain_no -> array index and pdb_chain_id -> array index
         my $pdb_chain_no_rev = pibase::array2hash($pdb_chain_no) ;
         my $pdb_chain_id_rev = pibase::array2hash($pdb_chain_id) ;

#Load resno_serial -> resno mapping for current pdb_id
         if (!exists $all_bdp_residues_tables->{$raw_bdpid}) {
            print STDERR "\nERROR: PDB $pdb_id ($raw_bdpid): bdp_residues ".
                         "table ($raw_bdpid) not found\n" ;
            next ;
         } else {
            $tables->{rawpdb_residues}->{name} =
               $all_bdp_residues_tables->{$raw_bdpid} ;
         }

#Read in bdp_residues entries for the raw pdb entry corresponding to
# current bdp entry.

         my $pdb_allres ;
         ( $pdb_allres->{chain_no},
           $pdb_allres->{chain_id},
           $pdb_allres->{resno_serial},
           $pdb_allres->{resno}) =
            pibase::rawselect_metatod($tables->{rawpdb_residues}->{name},
            "SELECT chain_no, chain_id, resno_serial, resno ".
            "FROM bdp_residues WHERE bdp_id = $raw_bdpid" ) ;

# Iterate through the rawpdb bdp_residues entries.

         my $rawpdb_resno_2_serial  ;
         my $rawpdb_ressig_seen ;
         my $rawpdb_duplicate_fl = 0 ;

         foreach my $j ( 0 .. $#{$pdb_allres->{resno}} ) {

#If this resno:chain_id has been seen before, increment the duplicate_residue flag.

            my $sig = $pdb_allres->{resno}->[$j]."\n".
                      $pdb_allres->{chain_id}->[$j] ;

            if (exists $rawpdb_ressig_seen->{$sig}) {
               $rawpdb_duplicate_fl++ ;
            } else {
#Otherwise, create a hash entry that points from resno:chain_id to the serial residue no.
               $rawpdb_resno_2_serial->{$sig}=$pdb_allres->{resno_serial}->[$j];
            }

#Increment the seen counter for this resno:chain_id

            $rawpdb_ressig_seen->{$sig}++ ;
         }

#ERRCATCH: If the duplicate residue flag has been set, display to STDERR and abort the current BDP entry.

         if ($rawpdb_duplicate_fl) {
            print STDERR "\nERROR: PDB $pdb_id ($raw_bdpid): ".
            "$rawpdb_duplicate_fl duplicate residue number:chain ids\n" ;
            next ;
         }

# Load subset_id, description, and source_id for the current pdb_id

         my ($subset_id, $subset_description, $subset_source_id, $subset_class);
         {
            my @tind = @{$allsubs_point->{$pdb_id}} ;
            push @{$subset_id}, @{$allsubs->{subset_id}}[@tind] ;
            push @{$subset_description}, @{$allsubs->{description}}[@tind] ;
            push @{$subset_source_id}, @{$allsubs->{source_id}}[@tind] ;
            push @{$subset_class}, @{$allsubs->{class}}[@tind] ;
         }


         my @subsets_chain_id ;
         my (@subsets_start, @subsets_end) ;
         my (@subsets_start_ser, @subsets_end_ser) ;

#Iterate through the PDB's subsets
         foreach my $k (0 .. $#{$subset_id}) {
#Read in fragment information (chain id, start resno, end resno) for the current subset_id.
            {
               my @tind = @{$allsubsdet_point->{$subset_id->[$k]}} ;
               push @{$subsets_chain_id[$k]}, @{$allsubsdet->{chain_id}}[@tind] ;
               push @{$subsets_start[$k]}, @{$allsubsdet->{start}}[@tind] ;
               push @{$subsets_end[$k]}, @{$allsubsdet->{end}}[@tind] ;
            }
               
#Iterate through each subset segment.
            foreach my $l ( 0 .. $#{$subsets_start[$k]} ) {
#If the chain_id for this segment is not defined, set to the NULL string.
               if ((!defined $subsets_chain_id[$k]->[$l]) ||
                   ($subsets_chain_id[$k]->[$l] eq 'NULL')) {
                  $subsets_chain_id[$k]->[$l] = ' '; }

#ERRCATCH: If the PDB chain referred to by the current segment does not exist in the pdb_chains loaded from pibase.bdp_chains, display to STDERR and abort the current segment.

               my $curseg_pdb_chain_id =
                  $pdb_chain_id_rev->{$subsets_chain_id[$k]->[$l]} ;

               if ((!defined $curseg_pdb_chain_id) ||
                   ($curseg_pdb_chain_id eq 'NULL')) {
                  print STDERR "\nWARNING: PDB $pdb_id ($raw_bdpid): ".
                     "chain $subsets_chain_id[$k]->[$l] of subset ".
                     "$subset_id->[$k] does not have an entry in the".
                     " PDB chains; this chain skipped\n" ;
                  next ;
               }


#If the start residue is undefined, set to the start residue number stored for the current PDB chain.

               if ( (! defined $subsets_start[$k]->[$l] ) ||
                           ( $subsets_start[$k]->[$l] eq '' ) ||
                           ( $subsets_start[$k]->[$l] eq ' ') ) {
                  $subsets_start[$k]->[$l] =
                  $pdb_start->[$curseg_pdb_chain_id] ;
               }

#If the end residue is undefined, set to the end residue number for the current PDB chain.

               if ( (!defined $subsets_end[$k]->[$l]) ||
                    ($subsets_end[$k]->[$l] eq 'NULL') ||
                    ($subsets_end[$k]->[$l] eq '' ) ||
                    ($subsets_end[$k]->[$l] eq ' ' ) ) {
                    $subsets_end[$k]->[$l] = $pdb_end->[$curseg_pdb_chain_id] ;
               }


               my $start_serial = $rawpdb_resno_2_serial->{$subsets_start[$k]->[$l]."\n".$subsets_chain_id[$k]->[$l]} ;
               if (!defined $start_serial) {
                  print STDERR "\nERROR: PDB $pdb_id ($raw_bdpid): ".
                               "$subset_id->[$k] start residue serial number ".
                               "not found $subsets_start[$k]->[$l] on chain ".
                               "$subsets_chain_id[$k]->[$l]\n" ;
                  next;
               }
               $subsets_start_ser[$k]->[$l] = $start_serial ;

               my $end_serial = $rawpdb_resno_2_serial->{$subsets_end[$k]->[$l].
                                "\n".$subsets_chain_id[$k]->[$l]} ;
               if (!defined $end_serial) {
                  print STDERR "\nERROR: PDB $pdb_id ($raw_bdpid): ".
                               "$subset_id->[$k] start residue serial number ".
                               "not found $subsets_end[$k]->[$l] on chain ".
                               "$subsets_chain_id[$k]->[$l]\n";
                  next;
               }
               $subsets_end_ser[$k]->[$l] = $end_serial ;
            }
         }


         foreach my $bdp_id (keys %{$in->{bdps}->{$pdb_id}}) {


            my @chain_subsets_lines ;
            my @chain_subsdet_lines ;
            my @chain_subsres_lines ;
   
#ERRCATCH: If bdp_chains entry not loaded, STDERR and abort current BDP

            if (!exists $allchains_point->{$bdp_id}) {
               print STDERR "\nERROR: pdb $pdb_id (bdps ".
               join(",", (keys %{$in->{bdps}->{$pdb_id}})).
               ") could not find chain entries\n" ;
               next;
            }
   
   
#Load chain information for the current bdp_id entry into memory.
   
            my ($bdp_chain_no, $bdp_chain_id, $bdp_pdb_chain_no, $bdp_pdb_chain_id,
                $bdp_start, $bdp_end ) ;
            {
               my @tind = @{$allchains_point->{$bdp_id}} ;
               push @{$bdp_chain_no}, @{$allchains->{real_chain_no}}[@tind] ;
               push @{$bdp_chain_id}, @{$allchains->{real_chain_id}}[@tind] ;
               push @{$bdp_pdb_chain_no}, @{$allchains->{pdb_chain_no}}[@tind] ;
               push @{$bdp_pdb_chain_id}, @{$allchains->{pdb_chain_id}}[@tind] ;
               push @{$bdp_start}, @{$allchains->{start_resno}}[@tind] ;
               push @{$bdp_end}, @{$allchains->{end_resno}}[@tind] ;
            }
   
#ERRCATCH: If no BDP chains loaded, write to STDERR and abort current BDP entry.
   
            if ($#{$bdp_chain_no} < 0) {
               print STDERR "\nERROR: bdp $bdp_id ($pdb_id): ".
                  "no bdp_chains loaded\n" ;
               next ;
            }

#iterate over chains in the bdp file.
#  SUBSET: class undef - call it "pdbchain"
#  SUBSET_DETAILS: fields from bdp_chains table; 
#  SUBSET_RESIDUES: BDP_RESIDUES with an additional column at end with the
#   "subset_id" for the chain-domain
#
# have to roll in to translation routine - files are opened there...
#            iterate over chains = each chain = domain.
#            for each domain print to the files that have been opened

            foreach my $k (0 .. $#{$bdp_chain_no}) {
               my $cur_chain_id ;
               if ((!defined $bdp_chain_id->[$k]) ||
                   ($bdp_chain_id->[$k] eq ' ') ||
                   ($bdp_chain_id->[$k] eq '') ||
                   ($bdp_chain_id->[$k] eq 'NULL')) {
                  $cur_chain_id = '_' ;
               } else {
                  $cur_chain_id = $bdp_chain_id->[$k] ;
               }

               my $cur_chain_subset_id = "BDP".$bdp_id."-".$bdp_chain_no->[$k].
                  "_CHAIN-".$cur_chain_id ;

               my $subsets_vals = [ $cur_chain_subset_id, $bdp_id, '\N',
                                    '\N', $subset_source2id->{"pdb_chains"},
                                    'noclass' ] ;
               $subsets_vals = pibase::replace_undefs($subsets_vals, '\N') ;
               push @chain_subsets_lines, join("\t", @{$subsets_vals}) ;

               my ($start_int, undef) = residue_int($bdp_start->[$k])  ;
               my ($end_int, undef) = residue_int($bdp_end->[$k])  ;

               my $details_vals = [
                  $cur_chain_subset_id, 1,
                  $bdp_chain_no->[$k],
                  $bdp_chain_id->[$k],
                  $bdp_start->[$k],
                  $start_int,
                  $bdp_end->[$k],
                  $end_int,
               ] ;
               $details_vals = pibase::replace_undefs($details_vals, '\N') ;
               push @chain_subsdet_lines, join("\t", @{$details_vals}) ;
            }

# iterate over all residues in the bdp; essentially print directly to
#  subsets residues with an additional column describing subset_id


#Declare hashes for BDP chain_no -> array index and BDP chain_id -> array index.
   
            my $bdp_chain_no_rev = pibase::array2hash($bdp_chain_no);
            my $bdp_chain_id_rev = pibase::array2hash($bdp_chain_id, ' ');
   
#NOTE:
#* Initialize null chain_no (numnullchains) counter. The procedure ignores any BDP entries with at least one chain with an undefined pdb_chain_no ( bdp_ids with at least one chain with a null pdb_chain_no: 307, bdp_ids with at least one null and one non-null pdb_chain_no chain: 59)
#
#* Is it wise to ignore even BDP entries with at least one non-null PDB chain_no? Would it be better if we processed all non-null chains, regardless of unassigned sibling chains?

# Count how many of the real chains (bdp_chain_no) have undefined bdp_pdb_chain_no. if there is even one chain with an undefined pdb_chain_no, the current bdp_id is aborted.
   
            my $numnullchains = 0;
            foreach my $k (0 .. $#{$bdp_chain_no}) {
               if ((!defined $bdp_pdb_chain_no->[$k]) ||
                   ($bdp_pdb_chain_no->[$k] eq '') ||
                   ($bdp_pdb_chain_no->[$k] eq 'NULL')) {
                  $numnullchains++; } }
   
#ERRCATCH: If there are any chains with null chain_no (numnullchains > 0), display to STDERR the number of non-null and null chains read for this BDP entry, and abort the current BDP entry.
   
            if ($numnullchains) {
   
               my $numnotnull = $#{$bdp_chain_no} + 1 - $numnullchains ;
   
               print STDERR "\nERROR: bdp $bdp_id ($pdb_id): ".
                      "$numnullchains chains with NULL pdb_chain_no, ".
                      "$numnotnull chains with non-NULL pdb_chain_no\n";
               next ;
   
            }
   
#FPD051019_1636 - removed residue number diff code - bullshit
# new strategy: move to resno_pdb -> resno_serial_(pdb = real) -> resno_real
#
#=pod
#
#Declare the array to hold residue number differences between BDP and corresponding PDB chain.
#
#=cut
#
#      my $bdp_res_add ;
#
#=pod
#
#Iterate through the real chains (bdp_chain_no), Calculate the residue number difference.
#
#=cut
#
#      foreach my $k (0 .. $#{$bdp_chain_no}) {
#         my ($bdp_int, $bdp_inscode) = residue_int($bdp_start->[$k]);
#
#         my ($pdb_int, $pdb_inscode) = residue_int($pdb_start->[$pdb_chain_no_rev->{$bdp_pdb_chain_no->[$k]}]) ;
#         $bdp_res_add->[$k] = $bdp_int - $pdb_int ;
#      }

#Build reverse lookup hash from PDB chain id to real chain no
            my $pdb_chid_2_bdp_chno = $allchains_pdb_chid_2_bdp_chno->{$bdp_id} ;
   
   
#Determine bdp_residues table name for current bdp_id.
#ERRCATCH: if a bdp_residues table does note exist for the current bdp_id, display to STDERR, and abort current BDP entry.
            if (!exists $all_bdp_residues_tables->{$bdp_id}) {
               print STDERR "\nERROR: bdp $bdp_id ($pdb_id): ".
                            "bdp_residues table not found\n" ;
               next ;
            } else {
               $tables->{bdp_residues}->{name} = $all_bdp_residues_tables->{$bdp_id} ;
            }
   
#Read in bdp_residues entries for the current bdp entry.
   
            my $bdp_allres ;
            ( $bdp_allres->{chain_no},
              $bdp_allres->{chain_id},
              $bdp_allres->{resno_serial},
              $bdp_allres->{resno},
              $bdp_allres->{resno_int}) =
               pibase::rawselect_metatod($tables->{bdp_residues}->{name},
               "SELECT chain_no, chain_id, resno_serial, resno, resno_int ".
               "FROM bdp_residues WHERE bdp_id = $bdp_id" ) ;
   
#Iterate through the bdp_residues entries.
   
            my $resno_2_serial  ;
            my $resserial_2_no  ;
            my $resserial_chno_2_no  ;
            my $ressig_seen ;
            my $duplicate_res_fl = 0 ;
   
            foreach my $j ( 0 .. $#{$bdp_allres->{resno}} ) {

#print this to the chain domain subsets residues file
# determine chain subset id from the chain id, 

               my $cur_chain_id ;
               if ((!defined $bdp_allres->{chain_id}->[$j]) ||
                   ($bdp_allres->{chain_id}->[$j] eq ' ') ||
                   ($bdp_allres->{chain_id}->[$j] eq '') ||
                   ($bdp_allres->{chain_id}->[$j] eq 'NULL')) {
                  $cur_chain_id = '_' ;
               } else {
                  $cur_chain_id = $bdp_allres->{chain_id}->[$j] ;
               }
               my $cur_chain_subset_id = "BDP".$bdp_id."-".
                  $bdp_allres->{chain_no}->[$j]."_CHAIN-".$cur_chain_id ;

               my $subsres_outvals = [ $bdp_id,
                               $bdp_allres->{chain_no}->[$j],
                               $bdp_allres->{chain_id}->[$j],
                               $bdp_allres->{resno_serial}->[$j],
                               $bdp_allres->{resno}->[$j],
                               $bdp_allres->{resno_int}->[$j],
                               $cur_chain_subset_id ] ;
               push @chain_subsres_lines, join("\t", @{$subsres_outvals}) ;

#If this resno:chain_id has been seen before, increment the duplicate_residue flag.
   
               my $sig = $bdp_allres->{resno}->[$j]."\n".
                         $bdp_allres->{chain_id}->[$j] ;
   
               if (exists $ressig_seen->{$sig}) {
                  $duplicate_res_fl++ ;
               } else {
#Otherwise, create a hash entry that points from resno:chain_id to the serial residue no.
                  $resno_2_serial->{$sig} = $bdp_allres->{resno_serial}->[$j] ;
                  $resserial_2_no->{$bdp_allres->{chain_id}->[$j]}->{$bdp_allres->{resno_serial}->[$j]} = $sig ;
                  $resserial_chno_2_no->{$bdp_allres->{chain_no}->[$j]}->{$bdp_allres->{resno_serial}->[$j]} = $sig ;
               }
   
#Increment the seen counter for this resno:chain_id
               $ressig_seen->{$sig}++ ;
   
            }
   
#ERRCATCH: If the duplicate residue flag has been set, display to STDERR and abort the current BDP entry.
   
            if ($duplicate_res_fl) {
               print STDERR "\nERROR: bdp $bdp_id ($pdb_id): $duplicate_res_fl".
                            " duplicate residue number:chain ids\n" ;
               next ;
            }
   
#Summary: convert subsets_details residue numbering from pdb to bdp and subset_id numbering to allow for multiple BDP instances of a single PDB domain
   
#Create arrays to hold translated BDP subset definitions.
   
            my @bdp_subsets_id ;
            my %bdp_subsets_id ;
            my @bdp_subsets_descr ;
            my @bdp_subsets_class ;
            my @bdp_subsets_source_id ;
            
            my @bdp_subsets_chain_no ;
            my @bdp_subsets_chain_id ;
            my @bdp_subsets_start ;
            my @bdp_subsets_end ;
   
#Iterate through the PDB subsets (subset_id)
            foreach my $k ( 0 .. $#{$subset_id} ) {
   
#Iterate through the segments of the current PDB subset
               foreach my $l ( 0 .. $#{$subsets_chain_id[$k]}) {
                  my $pdb_ch_id = $subsets_chain_id[$k]->[$l] ;
                  my $m = 0;
   
# Iterate through the BDP chains equivalent to the current PDB chain.
#
#BUG Assumes that domain assignment are consecutive; that is
#
#PDB domain Chain_id
#CATH.1 A,B
#CATH.2 C
#
#
#PDBch BDPeq
#A 1,4
#B 2,5
#C 3,6
#
#Assumes BDP0.CATH.1 is 1,2 and BDP1.CATH.1 is 4,5.
   
                  foreach my $bdp_chno (@{$pdb_chid_2_bdp_chno->{$pdb_ch_id}}) {
   
                     if (!defined $subsets_start_ser[$k]->[$l] ||
                         !defined $subsets_end_ser[$k]->[$l]) {
                        print STDERR "\nERROR: bdp $bdp_id ($pdb_id): ".
                           " skipping $subset_id->[$k] derivative assignment ".
                           " since serial number of start or end is undefined\n";
                        next;
                     }
                     
                     my $t_bdp_charrno = $bdp_chain_no_rev->{$bdp_chno} ;
   
#Create a BDP subset entry with a subset_id of BDP(current ordinal BDP chain counter)_(PDB subset_id), subsets_class = PDB subset_class, subsets_chain_no = current BDP chain no, subsets_chain_id = current BDP chain id.
   
                     my $bdp_subset_id = "BDP".$bdp_id."-".$m."_".$subset_id->[$k] ;
   
                     if (!exists $bdp_subsets_id{$bdp_subset_id} ) {
                        push @bdp_subsets_id, $bdp_subset_id ;
                        $bdp_subsets_id{$bdp_subset_id} = $#bdp_subsets_id ;
   
                        push @bdp_subsets_descr, $subset_description->[$k] ;
                        push @bdp_subsets_source_id, $subset_source_id->[$k] ;
                        push @bdp_subsets_class, $subset_class->[$k] ;
   
                        my $t_bdp_subset_no = $bdp_subsets_id{$bdp_subset_id} ;
   
                        $bdp_subsets_chain_no[$t_bdp_subset_no] = [ ] ;
                        $bdp_subsets_chain_id[$t_bdp_subset_no] = [ ] ;
                        $bdp_subsets_start[$t_bdp_subset_no] = [ ] ;
                        $bdp_subsets_end[$t_bdp_subset_no] = [ ] ;
   
                     }
                     
                     my $bdp_subset_no = $bdp_subsets_id{$bdp_subset_id} ;

                     push @{$bdp_subsets_chain_no[$bdp_subset_no]}, $bdp_chno ;
                     push @{$bdp_subsets_chain_id[$bdp_subset_no]},
                        $bdp_chain_id->[$t_bdp_charrno] ;
   
#Generate the BDP subset start and end residue by adding the residue difference calculated earlier to the PDB subset start and end residue.
   
                     my ($t_startres, $t_endres, $t_startch, $t_endch) ;
                     if (!exists $resserial_chno_2_no->{$bdp_chain_no->[$t_bdp_charrno]}->{$subsets_start_ser[$k]->[$l]}) {
                        $t_startres = 'MISSING' ;
                        $t_endres = 'MISSING' ;
                     } elsif (!exists $resserial_chno_2_no->{$bdp_chain_no->[$t_bdp_charrno]}->{$subsets_end_ser[$k]->[$l]}) {
                        $t_startres = 'MISSING' ;
                        $t_endres = 'MISSING' ;
                     } else {
                        my $t_startsig = $resserial_chno_2_no->{$bdp_chain_no->[$t_bdp_charrno]}->{$subsets_start_ser[$k]->[$l]};
                        my $t_endsig = $resserial_chno_2_no->{$bdp_chain_no->[$t_bdp_charrno]}->{$subsets_end_ser[$k]->[$l]} ;
   
                        ($t_startres, $t_startch) = split(/\n/, $t_startsig) ;
                        ($t_endres, $t_endch) = split(/\n/, $t_endsig) ;
                     }
                     
                     push @{$bdp_subsets_start[$bdp_subset_no]}, $t_startres ;
                     push @{$bdp_subsets_end[$bdp_subset_no]}, $t_endres ;
   
                     $m++ ;
   
                  }
               }
            }
   
#Set and open output file for subsets_residues. Display filename to STDOUT
#Determine table name for subsets residue listing.
   
            $tables->{subsets_residues}->{name} =
               $tables->{subsets_residues}->{namepre}."_$bdp_id" ;
   
            my ($subsres_fh, $subsres_f, $subsres_depositdir) ;
            if (!FLAG_NOSUBSRES) {
               $subsres_f = $tables->{subsets_residues}->{name}.
                            ".$bdp_id.$pibase" ;
   
               open($subsres_fh, ">".$subsres_f) ;
               $subsres_depositdir =
                  $pibase_specs->{metatod_dir}->{subsets_residues}.'/'.
                  POSIX::floor($bdp_id / 1000) ;
               my @outvals = ($bdp_id, '\N',
                  $subsres_depositdir.'/'.$subsres_f.'.gz') ;
               print "SUBSETS_RESIDUES_TABLES\t".join("\t", @outvals)."\n" ;
            }
   
# Iterate through BDP subsets.
            foreach my $k ( 0 .. $#bdp_subsets_id ) {
   
# Initialize the problems flag for this subset;
               my $subset_problems_fl = 0 ;
   
# Setup temporary buffers to print out to proper files (subsets_fh, subsdets_fh, subres_fh) if no problems. To prevent partial subset succesfully printing (e.g. 1 segment of 5) ;
   
               my @subsets_lines ;
               my @subsdet_lines ;
               my @subsres_lines ;
   
   
# Set the display fields for pibase.subsets import: the BDP subset_id, bdp_id, null, description, subset_source_id, subset class. change undefined values to '\N'
               my $subsets_vals = [ $bdp_subsets_id[$k], $bdp_id, '\N',
                                 $bdp_subsets_descr[$k],
                                 $bdp_subsets_source_id[$k],
                                 $bdp_subsets_class[$k] ] ;
   
   
               $subsets_vals = pibase::replace_undefs($subsets_vals, '\N') ;
   
# Store subsets information into the pibase.subsets output buffer.
               push @subsets_lines, join("\t", @{$subsets_vals}) ;
   
# Iterate through the segments of this BDP subset id.
               foreach my $l ( 0 .. $#{$bdp_subsets_chain_id[$k]}) {
   
# Store the BDP subset_id, bdp_id, null, description, subset_source_id, subset class into the pibase.subsets_details buffer .
                  if ($bdp_subsets_start[$k]->[$l] eq 'MISSING') {
                     print STDERR "ERROR $bdp_id : the start residue is undefined on chain $bdp_subsets_chain_id[$k]->[$l] ($bdp_subsets_id[$k])\n" ;
                     $subset_problems_fl = 1;  last;
                  }
   
                  if ($bdp_subsets_end[$k]->[$l] eq 'MISSING') {
                     print STDERR "ERROR $bdp_id : the end residue is undefined on chain $bdp_subsets_chain_id[$k]->[$l] ($bdp_subsets_id[$k])\n" ;
                     $subset_problems_fl = 1;  last;
                  }
   
                  my ($start_int, undef) = residue_int($bdp_subsets_start[$k]->[$l]) ;
                  my ($end_int, undef) = residue_int($bdp_subsets_end[$k]->[$l]) ;
   
                  my $details_vals = [ $bdp_subsets_id[$k], $l,
                                       $bdp_subsets_chain_no[$k]->[$l],
                                       $bdp_subsets_chain_id[$k]->[$l],
                                       $bdp_subsets_start[$k]->[$l],
                                       $start_int, $bdp_subsets_end[$k]->[$l],
                                       $end_int ] ;
   
                  $details_vals = pibase::replace_undefs($details_vals, '\N') ;
   
                  push @subsdet_lines, join("\t", @{$details_vals}) ;
   
# Check whether the start and end residues for this subset segment are defined in bdp_residues. If not, set the problems flag.
                  my $start_ressig = $bdp_subsets_start[$k]->[$l]."\n".
                                     $bdp_subsets_chain_id[$k]->[$l] ;
   
                  my $end_ressig = $bdp_subsets_end[$k]->[$l]."\n".
                                    $bdp_subsets_chain_id[$k]->[$l] ;
   
                  if ( (! exists $resno_2_serial->{$start_ressig}) ||
                      (! exists $resno_2_serial->{$end_ressig}) ) {
                     print STDERR "ERROR $bdp_id :couldnt find resser entry for ".
                        "$bdp_subsets_start[$k]->[$l] on chain ".
                        "$bdp_subsets_chain_id[$k]->[$l] ($bdp_subsets_id[$k])\n";
                     $subset_problems_fl = 1 ; }
   
# If the problems flag has been set, abort this subset
                  if ($subset_problems_fl) {
                     last ; }
   
# Iterate over the residue listing for this bdp entry
                  foreach my $m ( 0 .. $#{$bdp_allres->{resno}} ) {
   
# If the current residue is covered by the current subset segment, store output in the pibase.subsets_residues buffer.
                     if ( ($bdp_allres->{chain_no}->[$m] ==
                           $bdp_subsets_chain_no[$k]->[$l]) &&
   
                           ($bdp_allres->{resno_serial}->[$m] >=
                           $resno_2_serial->{$start_ressig}) &&
                        
                           ($bdp_allres->{resno_serial}->[$m] <=
                           $resno_2_serial->{$end_ressig}) ) {
   
                        my $outvals = [ $bdp_id,
                                    $bdp_allres->{chain_no}->[$m],
                                    $bdp_allres->{chain_id}->[$m],
                                    $bdp_allres->{resno_serial}->[$m],
                                    $bdp_allres->{resno}->[$m],
                                    $bdp_allres->{resno_int}->[$m],
                                    $bdp_subsets_id[$k] ] ;
   
                        push @subsres_lines, join("\t", @{$outvals}) ;
   
                     }
                  }
               }
   
# If no problems were encountered, print out the information for this subset.
               if (! $subset_problems_fl) {
                  foreach my $l (0 .. $#subsets_lines) {
                     print "SUBSETS\t".$subsets_lines[$l]."\n"; }
   
                  foreach my $l (0 .. $#subsdet_lines) {
                     print "SUBSETS_DETAILS\t".$subsdet_lines[$l]."\n"; }
   
                  if (!FLAG_NOSUBSRES) {
                     foreach my $l (0 .. $#subsres_lines) {
                        print $subsres_fh $subsres_lines[$l]."\n"; } }
               }
            }

            foreach my $l (0 .. $#chain_subsets_lines) {
               print "SUBSETS\t".$chain_subsets_lines[$l]."\n"; }

            foreach my $l (0 .. $#chain_subsdet_lines) {
               print "SUBSETS_DETAILS\t".$chain_subsdet_lines[$l]."\n"; }

            if (!FLAG_NOSUBSRES) {
               foreach my $l (0 .. $#chain_subsres_lines) {
                  print $subsres_fh $chain_subsres_lines[$l]."\n"; } }

# Close subsets_residues file.
            if (!FLAG_NOSUBSRES) {
               close ($subsres_fh) ;
               system("gzip $subsres_f") ;
               my $t_dir = $pibase_specs->{metatod_dir}->{subsets_residues}.'/'.
                  POSIX::floor($bdp_id / 1000) ;
               if (!-s $t_dir) { mkpath $t_dir; }
               push @movethese, $subsres_f.'.gz' ;
               push @movehere, $t_dir ;
            }
   
            if ($#movethese == $movethreshold) {
               foreach my $j (0 .. $#movethese) {
                  pibase::safe_move($movethese[$j], $movehere[$j]) ; }
               @movethese = () ;
               @movehere = () ;
            }
         }

         $lastpdb = $pdb_id ;
      }

      foreach my $j (0 .. $#movethese) {
         pibase::safe_move($movethese[$j], $movehere[$j]) ; }

# Close subsets and subsets_details output files. - old style - now goes to STDOUT
#      close ($subsets_fh) ;
#      close ($subsdet_fh) ;

      print STDERR "\n" ;

   }

}


sub OLD_bdp_subset_translator {

#   use constant FLAG_NOSUBSRES => 1 ;

   my $movethreshold=500 ;
   my $specs = pibase::get_specs() ;
   my $tables_dir = $specs->{tod_dir} ;

   my $subsres_depositdir = $specs->{metatod_dir}->{subsets_residues} ;


   my $pibase = pibase::dbname() ;
   my $hostname = hostname() ;


# Set names of tables that need to be generated.
   my $tables ;

   $tables->{subsets}->{name} = 'subsets' ;
   $tables->{subsets_details}->{name} = 'subsets_details' ;
   $tables->{subsets_residues}->{namepre} = 'subsets_residues' ;


# Load the pdb_ids and bdp_ids for the BDP entries into memory.
   print STDERR "reading in bdp_ids: ." ;
   my ($bdp_2_pdb_id, $path_2_bdp_id, $pdb_id_2_bdp_id) =
      pibase::todload_bdp_ids(
      'bdp_id_2_pdb_id', 'path_2_bdp_id', 'pdb_id_2_bdp_id') ;
   print STDERR "x\n" ;

# Load the entire bdp_chains array into memory and index with bdp_id.
   print STDERR "reading in bdp_chains: " ;
   my $allchains ;
   ( $allchains->{bdp_id},
     $allchains->{real_chain_no},
     $allchains->{real_chain_id},
     $allchains->{pdb_chain_no},
     $allchains->{pdb_chain_id},
     $allchains->{start_resno},
     $allchains->{end_resno} ) =
      pibase::rawselect_tod( "SELECT bdp_id, real_chain_no, ".
      "real_chain_id, pdb_chain_no, pdb_chain_id, start_resno, end_resno ".
      "FROM bdp_chains") ;
   print STDERR $#{$allchains->{bdp_id}}." chains\n" ;

   $allchains->{real_chain_id} = pibase::replace_undefs($allchains->{real_chain_id}, ' ') ;
   $allchains->{pdb_chain_id} = pibase::replace_undefs($allchains->{pdb_chain_id}, ' ') ;


   my $allchains_point ;
   my $allchains_pdb_chid_2_bdp_chid ;
   my $allchains_pdb_chid_2_bdp_chno ;
   foreach my $j (0 .. $#{$allchains->{bdp_id}}) {
      push @{$allchains_point->{$allchains->{bdp_id}->[$j]}}, $j ;
      push @{$allchains_pdb_chid_2_bdp_chid->{$allchains->{bdp_id}->[$j]}->{$allchains->{pdb_chain_id}->[$j]}},
         $allchains->{real_chain_id}->[$j] ;
      push @{$allchains_pdb_chid_2_bdp_chno->{$allchains->{bdp_id}->[$j]}->{$allchains->{pdb_chain_id}->[$j]}},
         $allchains->{real_chain_no}->[$j] ;
   }

# Load the entire subsets table into memory.

   print STDERR "reading in subsets: \n" ;
   my $allsubs ;
   ( $allsubs->{pdb_id},
     $allsubs->{subset_id},
     $allsubs->{description},
     $allsubs->{source_id},
     $allsubs->{class} ) =
      pibase::rawselect_tod("SELECT pdb_id, subset_id, description, ".
      "subset_source_id, class FROM subsets") ;
   print STDERR "x\n" ;

   print STDERR $#{$allsubs->{pdb_id}}." subsets\n" ;

   my $allsubs_point ;
   foreach my $j (0 .. $#{$allsubs->{pdb_id}}) {
      push @{$allsubs_point->{$allsubs->{pdb_id}->[$j]}}, $j ; }

# Load the entire bdp_residues_tables into a hash.
   print STDERR "reading in bdp_residue_tables: ." ;
   my $all_bdp_residues_tables ;
   {
      my ($t_bdpid, $t_sourcefile) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM bdp_residues_tables") ;
      foreach my $j ( 0 .. $#{$t_bdpid}) {
         $all_bdp_residues_tables->{$t_bdpid->[$j]} = $t_sourcefile->[$j] ; }
   }
   print STDERR "x\n" ;

# Load the entire subsets_details table into memory.
   print STDERR "reading in subsets_details: " ;
   my $allsubsdet ;
   ( $allsubsdet->{subset_id},
     $allsubsdet->{chain_id},
     $allsubsdet->{start},
     $allsubsdet->{end} ) =
      pibase::rawselect_tod( "SELECT subset_id, chain_id, ".
      "start_resno, end_resno FROM subsets_details") ;
   print STDERR $#{$allsubsdet->{subset_id}}." subsets_details\n" ;

   $allsubsdet->{chain_id} = pibase::replace_undefs_blanks($allsubsdet->{chain_id}, ' ') ;

   my $allsubsdet_point ;
   foreach my $j (0 .. $#{$allsubsdet->{subset_id}}) {
      push @{$allsubsdet_point->{$allsubsdet->{subset_id}->[$j]}}, $j ; }

# Open output files for subsets and subsets_details entries.
   my ($subsets_fh, $subsets_f) =
      tempfile($tables->{subsets}->{name}.".$hostname.XXXX",
               SUFFIX => ".$pibase");

   my ($subsdet_fh, $subsdet_f) =
      tempfile($tables->{subsets_details}->{name}.".$hostname.XXXX",
               SUFFIX => ".$pibase");

   my $temp_dir = tempdir(CLEANUP=>1) ;
   chdir $temp_dir ;
   my @movethese ;

# OUTPUT: Display the file names.
   print "$subsets_f\t$subsdet_f\n" ;

# Read in bdp_paths to be processed
   my $in;
   while (my $bdp_path= <STDIN>) {
      chomp $bdp_path ;
      my $bdp_id = $path_2_bdp_id->{$bdp_path} ;
      my $pdb_id = $bdp_2_pdb_id->{$bdp_id} ;
      $in->{bdps}->{$pdb_id}->{$bdp_id} = $bdp_path ;
   }


#Iterate through the bdp_ids to be processed
   my $lastpdb = 'a';
   print STDERR "now on:  " ;
   foreach my $pdb_id (keys %{$in->{bdps}}) {

      print STDERR "\b"x(length($lastpdb)).$pdb_id ;

#ERRCATCH: If no subsets defined for this pdb_id, STDERR and abort BDP.
      if (! exists $allsubs_point->{$pdb_id}) {
         print STDERR "\nERROR: pdb $pdb_id: no subsets found\n" ;
         next;
      }

#Load pdb chain number, id, start and end residue numbers for the current pdb_id
      my $raw_bdpid = $pdb_id_2_bdp_id->{$pdb_id} ;

      if (! defined $raw_bdpid) {
         print STDERR "\nERROR: pdb $pdb_id: could not find bdp_id entry" ;
         next ;
      } elsif ( !exists $allchains_point->{$raw_bdpid} ) {
         print STDERR "\nERROR: pdb $pdb_id: could not find chains entries\n" ;
         next ;
      }


      my ($pdb_chain_no, $pdb_chain_id, $pdb_start, $pdb_end) ;
      {
         my @tind = @{$allchains_point->{$raw_bdpid}} ;
         push @{$pdb_chain_no}, @{$allchains->{real_chain_no}}[@tind] ;
         push @{$pdb_chain_id}, @{$allchains->{real_chain_id}}[@tind] ;
         push @{$pdb_start}, @{$allchains->{start_resno}}[@tind] ;
         push @{$pdb_end}, @{$allchains->{end_resno}}[@tind] ;
      }

#ERRCATCH: If no PDB chains loaded, STDERR and abort current PDB.
      if ($#{$pdb_chain_no} < 0) {
         print STDERR "\nERROR: PDB $pdb_id: no pdb_chains loaded\n" ;
         next ;
      }

#Declare hashes for pdb_chain_no -> array index and pdb_chain_id -> array index
      my $pdb_chain_no_rev = pibase::array2hash($pdb_chain_no) ;
      my $pdb_chain_id_rev = pibase::array2hash($pdb_chain_id) ;

#Load resno_serial -> resno mapping for current pdb_id
      if (!exists $all_bdp_residues_tables->{$raw_bdpid}) {
         print STDERR "\nERROR: PDB $pdb_id ($raw_bdpid): bdp_residues ".
                      "table ($raw_bdpid) not found\n" ;
         next ;
      } else {
         $tables->{rawpdb_residues}->{name} = $all_bdp_residues_tables->{$raw_bdpid} ;
      }

#Read in bdp_residues entries for the raw pdb entry corresponding to current bdp entry.

      my $pdb_allres ;
      ( $pdb_allres->{chain_no},
        $pdb_allres->{chain_id},
        $pdb_allres->{resno_serial},
        $pdb_allres->{resno}) =
         pibase::rawselect_metatod($tables->{rawpdb_residues}->{name},
         "SELECT chain_no, chain_id, resno_serial, resno ".
         "FROM bdp_residues WHERE bdp_id = $raw_bdpid" ) ;

# Iterate through the rawpdb bdp_residues entries.

      my $rawpdb_resno_2_serial  ;
      my $rawpdb_ressig_seen ;
      my $rawpdb_duplicate_fl = 0 ;

      foreach my $j ( 0 .. $#{$pdb_allres->{resno}} ) {

#If this resno:chain_id has been seen before, increment the duplicate_residue flag.

         my $sig = $pdb_allres->{resno}->[$j]."\n".
                   $pdb_allres->{chain_id}->[$j] ;

         if (exists $rawpdb_ressig_seen->{$sig}) {
            $rawpdb_duplicate_fl++ ;

         } else {

#Otherwise, create a hash entry that points from resno:chain_id to the serial residue no.

            $rawpdb_resno_2_serial->{$sig} = $pdb_allres->{resno_serial}->[$j] ;
         }

#Increment the seen counter for this resno:chain_id

         $rawpdb_ressig_seen->{$sig}++ ;
      }

#ERRCATCH: If the duplicate residue flag has been set, display to STDERR and abort the current BDP entry.

      if ($rawpdb_duplicate_fl) {
         print STDERR "\nERROR: PDB $pdb_id ($raw_bdpid): $rawpdb_duplicate_fl".
                      " duplicate residue number:chain ids\n" ;
         next ;
      }

# Load subset_id, description, and source_id for the current pdb_id

      my ( $subset_id, $subset_description, $subset_source_id, $subset_class );
      {
         my @tind = @{$allsubs_point->{$pdb_id}} ;
         push @{$subset_id}, @{$allsubs->{subset_id}}[@tind] ;
         push @{$subset_description}, @{$allsubs->{description}}[@tind] ;
         push @{$subset_source_id}, @{$allsubs->{source_id}}[@tind] ;
         push @{$subset_class}, @{$allsubs->{class}}[@tind] ;
      }


      my @subsets_chain_id ;
      my (@subsets_start, @subsets_end) ;
      my (@subsets_start_ser, @subsets_end_ser) ;

#Iterate through the PDB's subsets
      foreach my $k (0 .. $#{$subset_id}) {

#Read in fragment information (chain id, start resno, end resno) for the current subset_id.

         {
            my @tind = @{$allsubsdet_point->{$subset_id->[$k]}} ;
            push @{$subsets_chain_id[$k]}, @{$allsubsdet->{chain_id}}[@tind] ;
            push @{$subsets_start[$k]}, @{$allsubsdet->{start}}[@tind] ;
            push @{$subsets_end[$k]}, @{$allsubsdet->{end}}[@tind] ;
         }
            
#Iterate through each subset segment.
         foreach my $l ( 0 .. $#{$subsets_start[$k]} ) {

#If the chain_id for this segment is not defined, set to the NULL string.
            if ((!defined $subsets_chain_id[$k]->[$l]) ||
                ($subsets_chain_id[$k]->[$l] eq 'NULL')) {
               $subsets_chain_id[$k]->[$l] = ' '; }

#ERRCATCH: If the PDB chain referred to by the current segment does not exist in the pdb_chains loaded from pibase.bdp_chains, display to STDERR and abort the current segment.

            my $curseg_pdb_chain_id =
               $pdb_chain_id_rev->{$subsets_chain_id[$k]->[$l]} ;

            if ((!defined $curseg_pdb_chain_id) ||
                ($curseg_pdb_chain_id eq 'NULL')) {
               print STDERR "\nWARNING: PDB $pdb_id ($raw_bdpid): ".
                            "chain $subsets_chain_id[$k]->[$l] of subset ".
                            "$subset_id->[$k] does not have an entry in the".
                            " PDB chains; this chain skipped\n" ;
               next ;
            }


#If the start residue is undefined, set to the start residue number stored for the current PDB chain.

            if ( (! defined $subsets_start[$k]->[$l] ) ||
                  ( $subsets_start[$k]->[$l] eq '' ) ||
                  ( $subsets_start[$k]->[$l] eq ' ') ) {
               $subsets_start[$k]->[$l] = $pdb_start->[$curseg_pdb_chain_id] ;
            }

#If the end residue is undefined, set to the end residue number for the current PDB chain.

            if ( (!defined $subsets_end[$k]->[$l]) ||
                 ($subsets_end[$k]->[$l] eq 'NULL') ||
                 ($subsets_end[$k]->[$l] eq '' ) ||
                 ($subsets_end[$k]->[$l] eq ' ' ) ) {
               $subsets_end[$k]->[$l] = $pdb_end->[$curseg_pdb_chain_id] ;
            }


            my $start_serial = $rawpdb_resno_2_serial->{$subsets_start[$k]->[$l]."\n".$subsets_chain_id[$k]->[$l]} ;
            if (!defined $start_serial) {
               print STDERR "\nERROR: PDB $pdb_id ($raw_bdpid): ".
                            "$subset_id->[$k] start residue serial number ".
                            "not found $subsets_start[$k]->[$l] on chain ".
                            "$subsets_chain_id[$k]->[$l]\n" ;
               next;
            }
            $subsets_start_ser[$k]->[$l] = $start_serial ;

            my $end_serial = $rawpdb_resno_2_serial->{$subsets_end[$k]->[$l]."\n".$subsets_chain_id[$k]->[$l]} ;
            if (!defined $end_serial) {
               print STDERR "\nERROR: PDB $pdb_id ($raw_bdpid): ".
                            "$subset_id->[$k] start residue serial number ".
                            "not found $subsets_end[$k]->[$l] on chain ".
                            "$subsets_chain_id[$k]->[$l]\n";
               next;
            }
            $subsets_end_ser[$k]->[$l] = $end_serial ;

         }
      }


    foreach my $bdp_id (keys %{$in->{bdps}->{$pdb_id}}) {

#ERRCATCH: If bdp_chains entry not loaded, STDERR and abort current BDP

      if (!exists $allchains_point->{$bdp_id}) {
         print STDERR "\nERROR: pdb $pdb_id (bdps ".
         join(",", (keys %{$in->{bdps}->{$pdb_id}})).
         ") could not find chain entries\n" ;
         next;
      }


#Load chain information for the current bdp_id entry into memory.

      my ( $bdp_chain_no, $bdp_chain_id, $bdp_pdb_chain_no, $bdp_pdb_chain_id,
           $bdp_start, $bdp_end ) ;

      {
         my @tind = @{$allchains_point->{$bdp_id}} ;
         push @{$bdp_chain_no}, @{$allchains->{real_chain_no}}[@tind] ;
         push @{$bdp_chain_id}, @{$allchains->{real_chain_id}}[@tind] ;
         push @{$bdp_pdb_chain_no}, @{$allchains->{pdb_chain_no}}[@tind] ;
         push @{$bdp_pdb_chain_id}, @{$allchains->{pdb_chain_id}}[@tind] ;
         push @{$bdp_start}, @{$allchains->{start_resno}}[@tind] ;
         push @{$bdp_end}, @{$allchains->{end_resno}}[@tind] ;
      }

#ERRCATCH: If no BDP chains loaded, write to STDERR and abort current BDP entry.

      if ($#{$bdp_chain_no} < 0) {
         print STDERR "\nERROR: bdp $bdp_id ($pdb_id): ".
                      "no bdp_chains loaded\n" ;
         next ;
      }


#Declare hashes for BDP chain_no -> array index and BDP chain_id -> array index.

      my $bdp_chain_no_rev = pibase::array2hash($bdp_chain_no);
      my $bdp_chain_id_rev = pibase::array2hash($bdp_chain_id, ' ');

#NOTE:
#* Initialize null chain_no (numnullchains) counter. The procedure ignores any BDP entries with at least one chain with an undefined pdb_chain_no ( bdp_ids with at least one chain with a null pdb_chain_no: 307, bdp_ids with at least one null and one non-null pdb_chain_no chain: 59)
#
#* Is it wise to ignore even BDP entries with at least one non-null PDB chain_no? Would it be better if we processed all non-null chains, regardless of unassigned sibling chains?

# Count how many of the real chains (bdp_chain_no) have undefined bdp_pdb_chain_no. if there is even one chain with an undefined pdb_chain_no, the current bdp_id is aborted.

      my $numnullchains = 0;
      foreach my $k (0 .. $#{$bdp_chain_no}) {
         if ((!defined $bdp_pdb_chain_no->[$k]) ||
             ($bdp_pdb_chain_no->[$k] eq '') ||
             ($bdp_pdb_chain_no->[$k] eq 'NULL')) {
            $numnullchains++; } }

#ERRCATCH: If there are any chains with null chain_no (numnullchains > 0), display to STDERR the number of non-null and null chains read for this BDP entry, and abort the current BDP entry.

      if ($numnullchains) {

         my $numnotnull = $#{$bdp_chain_no} + 1 - $numnullchains ;

         print STDERR "\nERROR: bdp $bdp_id ($pdb_id): ".
                      "$numnullchains chains with NULL pdb_chain_no, ".
                      "$numnotnull chains with non-NULL pdb_chain_no\n";
         next ;

      }

#FPD051019_1636 - removed residue number diff code - bullshit
# new strategy: move to resno_pdb -> resno_serial_(pdb = real) -> resno_real
#
#=pod
#
#Declare the array to hold residue number differences between BDP and corresponding PDB chain.
#
#=cut
#
#      my $bdp_res_add ;
#
#=pod
#
#Iterate through the real chains (bdp_chain_no), Calculate the residue number difference.
#
#=cut
#
#      foreach my $k (0 .. $#{$bdp_chain_no}) {
#         my ($bdp_int, $bdp_inscode) = residue_int($bdp_start->[$k]);
#
#         my ($pdb_int, $pdb_inscode) = residue_int($pdb_start->[$pdb_chain_no_rev->{$bdp_pdb_chain_no->[$k]}]) ;
#         $bdp_res_add->[$k] = $bdp_int - $pdb_int ;
#      }

#Build reverse lookup hash from PDB chain id to real chain no
      my $pdb_chid_2_bdp_chno = $allchains_pdb_chid_2_bdp_chno->{$bdp_id} ;


#Determine bdp_residues table name for current bdp_id.
#ERRCATCH: if a bdp_residues table does note exist for the current bdp_id, display to STDERR, and abort current BDP entry.
      if (!exists $all_bdp_residues_tables->{$bdp_id}) {
         print STDERR "\nERROR: bdp $bdp_id ($pdb_id): ".
                      "bdp_residues table not found\n" ;
         next ;
      } else {
         $tables->{bdp_residues}->{name} = $all_bdp_residues_tables->{$bdp_id} ;
      }

#Read in bdp_residues entries for the current bdp entry.

      my $bdp_allres ;
      ( $bdp_allres->{chain_no},
        $bdp_allres->{chain_id},
        $bdp_allres->{resno_serial},
        $bdp_allres->{resno},
        $bdp_allres->{resno_int}) =
         pibase::rawselect_metatod($tables->{bdp_residues}->{name},
         "SELECT chain_no, chain_id, resno_serial, resno, resno_int ".
         "FROM bdp_residues WHERE bdp_id = $bdp_id" ) ;

#Iterate through the bdp_residues entries.

      my $resno_2_serial  ;
      my $resserial_2_no  ;
      my $resserial_chno_2_no  ;
      my $ressig_seen ;
      my $duplicate_res_fl = 0 ;

      foreach my $j ( 0 .. $#{$bdp_allres->{resno}} ) {

#If this resno:chain_id has been seen before, increment the duplicate_residue flag.

         my $sig = $bdp_allres->{resno}->[$j]."\n".
                   $bdp_allres->{chain_id}->[$j] ;

         if (exists $ressig_seen->{$sig}) {
            $duplicate_res_fl++ ;

         } else {
#Otherwise, create a hash entry that points from resno:chain_id to the serial residue no.
            $resno_2_serial->{$sig} = $bdp_allres->{resno_serial}->[$j] ;
            $resserial_2_no->{$bdp_allres->{chain_id}->[$j]}->{$bdp_allres->{resno_serial}->[$j]} = $sig ;
            $resserial_chno_2_no->{$bdp_allres->{chain_no}->[$j]}->{$bdp_allres->{resno_serial}->[$j]} = $sig ;
         }

#Increment the seen counter for this resno:chain_id
         $ressig_seen->{$sig}++ ;

      }

#ERRCATCH: If the duplicate residue flag has been set, display to STDERR and abort the current BDP entry.

      if ($duplicate_res_fl) {
         print STDERR "\nERROR: bdp $bdp_id ($pdb_id): $duplicate_res_fl".
                      " duplicate residue number:chain ids\n" ;
         next ;
      }

#Summary: convert subsets_details residue numbering from pdb to bdp and subset_id numbering to allow for multiple BDP instances of a single PDB domain

#Create arrays to hold translated BDP subset definitions.

      my @bdp_subsets_id ;
      my %bdp_subsets_id ;
      my @bdp_subsets_descr ;
      my @bdp_subsets_class ;
      my @bdp_subsets_source_id ;
      
      my @bdp_subsets_chain_no ;
      my @bdp_subsets_chain_id ;
      my @bdp_subsets_start ;
      my @bdp_subsets_end ;

#Iterate through the PDB subsets (subset_id)
      foreach my $k ( 0 .. $#{$subset_id} ) {

#Iterate through the segments of the current PDB subset
         foreach my $l ( 0 .. $#{$subsets_chain_id[$k]}) {

            my $pdb_ch_id = $subsets_chain_id[$k]->[$l] ;
            my $m = 0;

# Iterate through the BDP chains equivalent to the current PDB chain.
#
#BUG Assumes that domain assignment are consecutive; that is
#
#PDB domain Chain_id
#CATH.1 A,B
#CATH.2 C
#
#
#PDBch BDPeq
#A 1,4
#B 2,5
#C 3,6
#
#Assumes BDP0.CATH.1 is 1,2 and BDP1.CATH.1 is 4,5.

            foreach my $bdp_chno (@{$pdb_chid_2_bdp_chno->{$pdb_ch_id}}) {

               if (!defined $subsets_start_ser[$k]->[$l] ||
                   !defined $subsets_end_ser[$k]->[$l]) {
                  print STDERR "\nERROR: bdp $bdp_id ($pdb_id): ".
                     " skipping $subset_id->[$k] derivative assignment ".
                     " since serial number of start or end is undefined\n";
                  next;
               }

               my $t_bdp_charrno = $bdp_chain_no_rev->{$bdp_chno} ;

#Create a BDP subset entry with a subset_id of BDP(current ordinal BDP chain counter)_(PDB subset_id), subsets_class = PDB subset_class, subsets_chain_no = current BDP chain no, subsets_chain_id = current BDP chain id.

               my $bdp_subset_id = "BDP".$bdp_id."-".$m."_".$subset_id->[$k] ;

               if (!exists $bdp_subsets_id{$bdp_subset_id} ) {
                  push @bdp_subsets_id, $bdp_subset_id ;
                  $bdp_subsets_id{$bdp_subset_id} = $#bdp_subsets_id ;
                  
                  push @bdp_subsets_descr, $subset_description->[$k] ;
                  push @bdp_subsets_source_id, $subset_source_id->[$k] ;
                  push @bdp_subsets_class, $subset_class->[$k] ;

                  my $t_bdp_subset_no = $bdp_subsets_id{$bdp_subset_id} ;
                  $bdp_subsets_chain_no[$t_bdp_subset_no] = [ ] ;
                  $bdp_subsets_chain_id[$t_bdp_subset_no] = [ ] ;
                  $bdp_subsets_start[$t_bdp_subset_no] = [ ] ;
                  $bdp_subsets_end[$t_bdp_subset_no] = [ ] ;
               }
               
               my $bdp_subset_no = $bdp_subsets_id{$bdp_subset_id} ;
               
               push @{$bdp_subsets_chain_no[$bdp_subset_no]}, $bdp_chno ;
               push @{$bdp_subsets_chain_id[$bdp_subset_no]},
                  $bdp_chain_id->[$t_bdp_charrno] ;

#Generate the BDP subset start and end residue by adding the residue difference calculated earlier to the PDB subset start and end residue.

               my ($t_startres, $t_endres, $t_startch, $t_endch) ;
               if (!exists $resserial_chno_2_no->{$bdp_chain_no->[$t_bdp_charrno]}->{$subsets_start_ser[$k]->[$l]}) {
                  $t_startres = 'MISSING' ;
                  $t_endres = 'MISSING' ;
               } elsif (!exists $resserial_chno_2_no->{$bdp_chain_no->[$t_bdp_charrno]}->{$subsets_end_ser[$k]->[$l]}) {
                  $t_startres = 'MISSING' ;
                  $t_endres = 'MISSING' ;
               } else {
                  my $t_startsig = $resserial_chno_2_no->{$bdp_chain_no->[$t_bdp_charrno]}->{$subsets_start_ser[$k]->[$l]};
                  my $t_endsig = $resserial_chno_2_no->{$bdp_chain_no->[$t_bdp_charrno]}->{$subsets_end_ser[$k]->[$l]} ;

                  ($t_startres, $t_startch) = split(/\n/, $t_startsig) ;
                  ($t_endres, $t_endch) = split(/\n/, $t_endsig) ;
               }

               push @{$bdp_subsets_start[$bdp_subset_no]}, $t_startres ;
               push @{$bdp_subsets_end[$bdp_subset_no]}, $t_endres ;

               $m++ ;

            }
         }
      }

#Set and open output file for subsets_residues. Display filename to STDOUT
#Determine table name for subsets residue listing.

      $tables->{subsets_residues}->{name} =
         $tables->{subsets_residues}->{namepre}."_$bdp_id" ;

      my ($subsres_fh, $subsres_f) ;
      if (!FLAG_NOSUBSRES) {
         ($subsres_fh, $subsres_f) =
         tempfile( $tables->{subsets_residues}->{name}.".$bdp_id.XXXX",
                   SUFFIX => ".$pibase" );

         my @outvals = ($bdp_id, '\N', $subsres_depositdir.'/'.$subsres_f.'.gz');
         print join("\t", @outvals)."\n" ;
      }

# Iterate through BDP subsets.
      foreach my $k ( 0 .. $#bdp_subsets_id ) {

# Initialize the problems flag for this subset;
         my $subset_problems_fl = 0 ;

# Setup temporary buffers to print out to proper files (subsets_fh, subsdets_fh, subres_fh) if no problems. To prevent partial subset succesfully printing (e.g. 1 segment of 5) ;

         my @subsets_lines ;
         my @subsdet_lines ;
         my @subsres_lines ;


# Set the display fields for pibase.subsets import: the BDP subset_id, bdp_id, null, description, subset_source_id, subset class. change undefined values to '\N'
         my $subsets_vals = [ $bdp_subsets_id[$k], $bdp_id, '\N',
                              $bdp_subsets_descr[$k],
                              $bdp_subsets_source_id[$k],
                              $bdp_subsets_class[$k] ] ;


         $subsets_vals = pibase::replace_undefs($subsets_vals, '\N') ;

# Store subsets information into the pibase.subsets output buffer.
         push @subsets_lines, join("\t", @{$subsets_vals}) ;

# Iterate through the segments of this BDP subset id.
         foreach my $l ( 0 .. $#{$bdp_subsets_chain_id[$k]}) {

# Store the BDP subset_id, bdp_id, null, description, subset_source_id, subset class into the pibase.subsets_details buffer .
            if ($bdp_subsets_start[$k]->[$l] eq 'MISSING') {
               print STDERR "ERROR $bdp_id : the start residue is undefined on chain $bdp_subsets_chain_id[$k]->[$l] ($bdp_subsets_id[$k])\n" ;
               $subset_problems_fl = 1;  last;
            }

            if ($bdp_subsets_end[$k]->[$l] eq 'MISSING') {
               print STDERR "ERROR $bdp_id : the end residue is undefined on chain $bdp_subsets_chain_id[$k]->[$l] ($bdp_subsets_id[$k])\n" ;
               $subset_problems_fl = 1;  last;
            }

            my ($start_int, undef) = residue_int($bdp_subsets_start[$k]->[$l]) ;
            my ($end_int, undef) = residue_int($bdp_subsets_end[$k]->[$l]) ;

            my $details_vals = [ $bdp_subsets_id[$k], $l,
                                 $bdp_subsets_chain_no[$k]->[$l],
                                 $bdp_subsets_chain_id[$k]->[$l],
                                 $bdp_subsets_start[$k]->[$l],
                                 $start_int, $bdp_subsets_end[$k]->[$l],
                                 $end_int ] ;

            $details_vals = pibase::replace_undefs($details_vals, '\N') ;

            push @subsdet_lines, join("\t", @{$details_vals}) ;

# Check whether the start and end residues for this subset segment are defined in bdp_residues. If not, set the problems flag.
            my $start_ressig = $bdp_subsets_start[$k]->[$l]."\n".
                               $bdp_subsets_chain_id[$k]->[$l] ;

            my $end_ressig = $bdp_subsets_end[$k]->[$l]."\n".
                              $bdp_subsets_chain_id[$k]->[$l] ;

            if ( (! exists $resno_2_serial->{$start_ressig}) ||
                 (! exists $resno_2_serial->{$end_ressig}) ) {
               print STDERR "ERROR $bdp_id :couldnt find resser entry for $bdp_subsets_start[$k]->[$l] on chain $bdp_subsets_chain_id[$k]->[$l] ($bdp_subsets_id[$k])\n" ;
               $subset_problems_fl = 1 ; }

# If the problems flag has been set, abort this subset
            if ($subset_problems_fl) {
               last ; }

# Iterate over the residue listing for this bdp entry
            foreach my $m ( 0 .. $#{$bdp_allres->{resno}} ) {

# If the current residue is covered by the current subset segment, store output in the pibase.subsets_residues buffer.
               if ( ($bdp_allres->{chain_no}->[$m] ==
                     $bdp_subsets_chain_no[$k]->[$l]) &&
                    ($bdp_allres->{resno_serial}->[$m] >=
                     $resno_2_serial->{$start_ressig}) &&

                    ($bdp_allres->{resno_serial}->[$m] <=
                     $resno_2_serial->{$end_ressig}) ) {

                  my $outvals = [ $bdp_id,
                                 $bdp_allres->{chain_no}->[$m],
                                 $bdp_allres->{chain_id}->[$m],
                                 $bdp_allres->{resno_serial}->[$m],
                                 $bdp_allres->{resno}->[$m],
                                 $bdp_allres->{resno_int}->[$m],
                                 $bdp_subsets_id[$k] ] ;

                  push @subsres_lines, join("\t", @{$outvals}) ;

               }
            }
         }

# If no problems were encountered, print out the information for this subset.
         if (! $subset_problems_fl) {
            foreach my $l (0 .. $#subsets_lines) {
               print $subsets_fh $subsets_lines[$l]."\n"; }

            foreach my $l (0 .. $#subsdet_lines) {
               print $subsdet_fh $subsdet_lines[$l]."\n"; }

            if (!FLAG_NOSUBSRES) {
               foreach my $l (0 .. $#subsres_lines) {
                  print $subsres_fh $subsres_lines[$l]."\n"; } }
         }
      }

# Close subsets_residues file.
      if (!FLAG_NOSUBSRES) {
         close ($subsres_fh) ;
         system("gzip $subsres_f") ;
         push @movethese, $subsres_f.'.gz' ;
      }

      if ($#movethese == $movethreshold) {
         foreach my $t_fn (@movethese) {
            pibase::safe_move($t_fn, $subsres_depositdir) ; }
         @movethese = () ;
      }

    }

    $lastpdb = $pdb_id ;
   }

   foreach my $t_fn (@movethese) {
      pibase::safe_move($t_fn, $subsres_depositdir) ; }

# Close subsets and subsets_details output files.

   close ($subsets_fh) ;
   close ($subsdet_fh) ;

   print STDERR "\n" ;

}


=head2 calc_bdp_chains_equiv()

   Title:       calc_bdp_chain_equiv()
   Function:    determine correspondence between PISA and parent PDB chains
   Args:        none
   Returns:     nothing
   STDIN:       bdp_file_path
   STDOUT:      1. bdp_id
                2. serial number of actual chain
                3. chain identifier of actual chain
                4. pdb_id
                5. serial number of corresponding PDB chain
                6. chain identifier of corresponding PDB chain

=cut

sub calc_bdp_chains_equiv {

   require DBI;
 
   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }

# Connect to the pibase database.
   my ($dbh, $pibase) = pibase::connect_pibase() ;

# Load bdp_ids
   my ($bdp_ids, $path_2_bdp_id, $bdp_id_2_pdb_id, $bdp_id_2_raw,
       $pdb_id_2_bdp_id, $bdp_id_2_path) = pibase::load_bdp_ids($dbh,
qw/bdp_id path_2_bdp_id bdp_id_2_pdb_id bdp_id_2_raw_pdb pdb_id_2_bdp_id bdp_id_2_path/) ;

# Load the entire bdp_chains array into memory and index with bdp_id
   my $allchains ;
   ( $allchains->{bdp_id} ,
     $allchains->{real_chain_no} ,
     $allchains->{real_chain_id} ,
     $allchains->{pdb_chain_no},
     $allchains->{pdb_chain_id},
     $allchains->{start_resno},
     $allchains->{end_resno},
     $allchains->{sequence} ) = pibase::mysql_fetchcols($dbh,
      "SELECT bdp_id, real_chain_no, real_chain_id, pdb_chain_no, ".
      "pdb_chain_id, start_resno, end_resno, sequence FROM bdp_chains") ;
 
   $allchains->{real_chain_id} =
      pibase::replace_undefs($allchains->{real_chain_id}, ' ') ;
   $allchains->{real_chain_id} =
      pibase::replace_char($allchains->{real_chain_id}, '', ' ') ;
   $allchains->{real_chain_id} =
      pibase::replace_char($allchains->{real_chain_id}, 'NULL', ' ') ;

   $allchains->{pdb_chain_id} =
      pibase::replace_undefs($allchains->{pdb_chain_id}, ' ') ;
   $allchains->{pdb_chain_id} =
      pibase::replace_char($allchains->{pdb_chain_id}, '', ' ') ;
   $allchains->{pdb_chain_id} =
      pibase::replace_char($allchains->{pdb_chain_id}, 'NULL', ' ') ;

   my $allchains_point ;
   foreach my $j (0 .. $#{$allchains->{bdp_id}}) {
      push @{$allchains_point->{$allchains->{bdp_id}->[$j]}}, $j ; }

   my @outlines ;
   foreach my $bdp_id (keys %{$bdp_id_2_path}) {
      my $file_path = $bdp_id_2_path->{$bdp_id} ;
# Determine bdp_id and pdb_id; Display to STDERR and abort entry if undefined
      if (!exists $allchains_point->{$bdp_id}) {
         print "ERROR: bdp $bdp_id: could not find a chains entry\n" ;
         next ;
      }
      my $pdb_id = $bdp_id_2_pdb_id->{$bdp_id} ;

# Load number, id, and sequence for all of the current bdp file's chains
      my ($real_chain_no, $real_chain_id, $real_chain_seq) ;
      {
         my @tind = @{$allchains_point->{$bdp_id}} ;
         push @{$real_chain_no}, @{$allchains->{real_chain_no}}[@tind] ;
         push @{$real_chain_id}, @{$allchains->{real_chain_id}}[@tind] ;
         push @{$real_chain_seq}, @{$allchains->{sequence}}[@tind] ;
      }


# Load number, id, and seq for all of the corresponding raw PDB file's chains.
      my ($pdb_chain_no, $pdb_chain_id, $pdb_chain_seq) ;
      my $raw_bdpid = $pdb_id_2_bdp_id->{$pdb_id} ;

      if (!defined $raw_bdpid) {
         print "ERROR: bdp $bdp_id: no $pdb_id raw pdb entry\n" ;
         next ;
      }

      if (!defined $allchains_point->{$raw_bdpid}) {
         print "ERROR: bdp $bdp_id: no bdp_chains entry for $pdb_id ".
         "raw pdb entry\n" ;
         next ;
      }


      {
         my @tind = @{$allchains_point->{$raw_bdpid}} ;
         push @{$pdb_chain_seq}, @{$allchains->{sequence}}[@tind] ;
         push @{$pdb_chain_no}, @{$allchains->{real_chain_no}}[@tind];
         push @{$pdb_chain_id}, @{$allchains->{real_chain_id}}[@tind] ;
      }

# Initialize the gotonext hash. This will store reasons to skip the current pdb file id.
      my %gotonext ;

# If number of pdb_chains is < 1, set the pdbnotfound flag in the gotonext hash.
      if ($#{$pdb_chain_no} < 0 ) {
         $gotonext{pdbnotfound}++ ;}

# If number of real chains is < 1, set the realnotfound flag in the gotonext hash.
      if ($#{$real_chain_no} < 0 ) {
         $gotonext{realnotfound}++ ;}

# Iterate through PDB chains. If the pdb chain_id is undefined or '', set it to ' '.
      foreach my $k (0 .. $#{$pdb_chain_no}) {
         if (!(defined $pdb_chain_id->[$k]) || ($pdb_chain_id->[$k] eq '')) {
            $pdb_chain_id->[$k] = ' ' ; }}

# Iterate through real chains. If the real chain_id is undefined or '', set it to ' '.
      foreach my $k (0 .. $#{$real_chain_no}) {
         if (!(defined $real_chain_id->[$k]) || ($real_chain_id->[$k] eq '')) {
            $real_chain_id->[$k] = ' ' ; } }

# If any gotonext hash flags have been set, print to STDERR the reason for aborting current bdp_id and skip to the next bdp_id.
      if ((keys %gotonext) > 0) {
         print "SKIPPED: bdp_id $bdp_id (PDB $pdb_id): ".
               join(' ,', (keys %gotonext))."\n" ;
         next;
      }

# Declare a sequence match matrix. r[i][j] denotes a sequence match between real chain i and PDB chain j.
      my $seq_match ;

# Iterate through real chains.
      foreach my $k (0 .. $#{$real_chain_no}) {

# Iterate through pdb chains.
         foreach my $l (0 .. $#{$pdb_chain_no}) {

# If the real chain sequence is equal to the PDB chain sequence, set the match matrix entry to 1 if the chain_ids are different, 2 if same.
            if ($real_chain_seq->[$k] eq $pdb_chain_seq->[$l]) {
               if ($real_chain_id->[$k] eq $pdb_chain_id->[$l]) {
                  $seq_match->[$k][$l] = 2 ;
               } else {
                  $seq_match->[$k][$l] = 1 ;
               }
# Else, set the match matrix entry to 0.
            } else {
               $seq_match->[$k][$l] = 0 ;
            }
         }
      }

      my ($equiv_ch, $equiv_ch_id) ;
      my $mismatch = 0;
      my $misoutput = "ERROR: bdp $bdp_id ($pdb_id)\n" ;

# Iterate through real chains.
      foreach my $k (0 .. $#{$real_chain_no}) {

# Find highest scoring entry in the sequence match matrix row (call max_ind(sequence match matrix row))
         my $t_bestmatch = max_ind($seq_match->[$k]) ;
         my $score = $seq_match->[$k][$t_bestmatch] ;

# Set the equivalent chain of this real chain to the highest scoring pdb chain.
         $equiv_ch->[$k] = $pdb_chain_no->[$t_bestmatch] ;
         $equiv_ch_id->[$k] = $pdb_chain_id->[$t_bestmatch] ;


# If the score is 0, set the mismatch flag to 1, as no equivalent PDB chain has been found.
         if ($score == 0) {
            $mismatch = 1;
            $misoutput .= "   REAL $real_chain_id->[$k]\t(".($k + 1).")\t$real_chain_seq->[$k]\n" ;
         } else { # Otherwise, display chain equivalence information.

            my @outvals = ($bdp_id, $real_chain_no->[$k], $real_chain_id->[$k],
                           $pdb_id, $equiv_ch->[$k], $equiv_ch_id->[$k]) ;
            push @outlines, \@outvals ;
         }
      }

# If mismatch flag is set, display all real and PDB chain sequences to STDERR.
      if ($mismatch) {
         print $misoutput."\n" ;
         foreach my $l (0 .. $#{$pdb_chain_no}) {
            print "   PDB $pdb_chain_id->[$l]\t(".($l +1).
                  ")\t$pdb_chain_seq->[$l]\n" ; }
      }
   }

   my $fh ;
   open($fh->{out}, '>'.$pibase_specs->{buildfiles}->{bdp_chains_equiv_update});
   foreach my $j (0 .. $#outlines) {
      my $outline = $outlines[$j] ;
      my ($bdp_id, $real_chain_no, $real_chain_id, $pdb_id,
          $pdb_chain_no, $pdb_chain_id) = @{$outline} ;

      my $whereclause  =  "bdp_id = $bdp_id" ;
      $whereclause .= " and real_chain_no = $real_chain_no" ;
 
      my $updateval = "pdb_chain_no = $pdb_chain_no" ;
      if ($pdb_chain_id ne '') {
         $updateval .= ", pdb_chain_id = \"$pdb_chain_id\"" ; }

# UPDATE bdp_chains SET pdb_chain_no = <pdb_chain_no> [, pdb_chain_id = <pdb_chain_id>] WHERE bdp_id = <bdp_id> and real_chain_no = <real_chain_no>

      my $chain_update_sql = "UPDATE bdp_chains SET $updateval ".
                              "WHERE $whereclause" ;
      print {$fh->{out}} $chain_update_sql."\n" ;

      if (exists $in->{import_fl} && $in->{import_fl} == 1) {
         $chain_update_sql = qq{$chain_update_sql} ;
      
         my $sth = $dbh->prepare($chain_update_sql) ;
         $sth->execute() ;
      }
   }
   close($fh->{out}) ;

}


sub OLD_bdp_chains_equiv_calc {

   require DBI;
 
# Connect to the pibase database.
   my ($dbh, $pibase) = pibase::connect_pibase() ;

# Load bdp_ids
   my ($bdp_ids, $path_2_bdp_id, $bdp_id_2_pdb_id, $bdp_id_2_raw, $pdb_id_2_bdp_id) =
      pibase::load_bdp_ids($dbh, qw/bdp_id path_2_bdp_id bdp_id_2_pdb_id bdp_id_2_raw_pdb pdb_id_2_bdp_id/) ;

# Load the entire bdp_chains array into memory and index with bdp_id
   my $allchains ;
   ( $allchains->{bdp_id} ,
     $allchains->{real_chain_no} ,
     $allchains->{real_chain_id} ,
     $allchains->{pdb_chain_no},
     $allchains->{pdb_chain_id},
     $allchains->{start_resno},
     $allchains->{end_resno},
     $allchains->{sequence} ) = pibase::mysql_fetchcols($dbh,
      "SELECT bdp_id, real_chain_no, real_chain_id, pdb_chain_no, ".
      "pdb_chain_id, start_resno, end_resno, sequence FROM bdp_chains") ;
 
   $allchains->{real_chain_id} =
      pibase::replace_undefs($allchains->{real_chain_id}, ' ') ;
   $allchains->{real_chain_id} =
      pibase::replace_char($allchains->{real_chain_id}, '', ' ') ;
   $allchains->{real_chain_id} =
      pibase::replace_char($allchains->{real_chain_id}, 'NULL', ' ') ;

   $allchains->{pdb_chain_id} =
      pibase::replace_undefs($allchains->{pdb_chain_id}, ' ') ;
   $allchains->{pdb_chain_id} =
      pibase::replace_char($allchains->{pdb_chain_id}, '', ' ') ;
   $allchains->{pdb_chain_id} =
      pibase::replace_char($allchains->{pdb_chain_id}, 'NULL', ' ') ;

   my $allchains_point ;
   foreach my $j (0 .. $#{$allchains->{bdp_id}}) {
      push @{$allchains_point->{$allchains->{bdp_id}->[$j]}}, $j ; }


   while (my $file_path = <STDIN>) {
      if ($file_path =~ /^#/) {next;}
      chomp $file_path ;

# Determine the bdp_id and pdb_id; Display to STDERR and abort entry if undefined
      my $bdp_id = $path_2_bdp_id->{$file_path} ;
      if (!defined $bdp_id) {
         print STDERR "ERROR: file $file_path: could not find bdp entry\n" ;
         next ;
      } elsif (!exists $allchains_point->{$bdp_id}) {
         print STDERR "ERROR: bdp $bdp_id: could not find a chains entry\n" ;
         next ;
      }
      my $pdb_id = $bdp_id_2_pdb_id->{$bdp_id} ;

      print STDERR "now on bdp $bdp_id\n" ;

# Load number, id, and sequence for all of the current bdp file's chains
      my ($real_chain_no, $real_chain_id, $real_chain_seq) ;
      {
         my @tind = @{$allchains_point->{$bdp_id}} ;
         push @{$real_chain_no}, @{$allchains->{real_chain_no}}[@tind] ;
         push @{$real_chain_id}, @{$allchains->{real_chain_id}}[@tind] ;
         push @{$real_chain_seq}, @{$allchains->{sequence}}[@tind] ;
      }


# Load number, id, and seq for all of the corresponding raw PDB file's chains.
      my ($pdb_chain_no, $pdb_chain_id, $pdb_chain_seq) ;
      my $raw_bdpid = $pdb_id_2_bdp_id->{$pdb_id} ;

      if (!defined $raw_bdpid) {
         print STDERR "ERROR: bdp $bdp_id: could not find $pdb_id raw pdb entry\n" ;
         next ;
      }


      {
         my @tind = @{$allchains_point->{$raw_bdpid}} ;
         push @{$pdb_chain_seq}, @{$allchains->{sequence}}[@tind] ;
         push @{$pdb_chain_no}, @{$allchains->{real_chain_no}}[@tind];
         push @{$pdb_chain_id}, @{$allchains->{real_chain_id}}[@tind] ;
      }

# Initialize the gotonext hash. This will store reasons to skip the current pdb file id.
      my %gotonext ;

# If number of pdb_chains is < 1, set the pdbnotfound flag in the gotonext hash.
      if ($#{$pdb_chain_no} < 0 ) {
         $gotonext{pdbnotfound}++ ;}

# If number of real chains is < 1, set the realnotfound flag in the gotonext hash.
      if ($#{$real_chain_no} < 0 ) {
         $gotonext{realnotfound}++ ;}

# Iterate through PDB chains. If the pdb chain_id is undefined or '', set it to ' '.
      foreach my $k (0 .. $#{$pdb_chain_no}) {
         if (!(defined $pdb_chain_id->[$k]) || ($pdb_chain_id->[$k] eq '')) {
            $pdb_chain_id->[$k] = ' ' ; }}

# Iterate through real chains. If the real chain_id is undefined or '', set it to ' '.
      foreach my $k (0 .. $#{$real_chain_no}) {
         if (!(defined $real_chain_id->[$k]) || ($real_chain_id->[$k] eq '')) {
            $real_chain_id->[$k] = ' ' ; } }

# If any gotonext hash flags have been set, print to STDERR the reason for aborting current bdp_id and skip to the next bdp_id.
      if ((keys %gotonext) > 0) {
         print STDERR "SKIPPED: bdp_id $bdp_id (PDB $pdb_id): ".join(' ,', (keys %gotonext))."\n" ;
         next;
      }

# Declare a sequence match matrix. r[i][j] denotes a sequence match between real chain i and PDB chain j.
      my $seq_match ;

# Iterate through real chains.
      foreach my $k (0 .. $#{$real_chain_no}) {

# Iterate through pdb chains.
         foreach my $l (0 .. $#{$pdb_chain_no}) {

# If the real chain sequence is equal to the PDB chain sequence, set the match matrix entry to 1 if the chain_ids are different, 2 if same.
            if ($real_chain_seq->[$k] eq $pdb_chain_seq->[$l]) {
               if ($real_chain_id->[$k] eq $pdb_chain_id->[$l]) {
                  $seq_match->[$k][$l] = 2 ;
               } else {
                  $seq_match->[$k][$l] = 1 ;
               }
# Else, set the match matrix entry to 0.
            } else {
               $seq_match->[$k][$l] = 0 ;
            }
         }
      }

      my ($equiv_ch, $equiv_ch_id) ;
      my $mismatch = 0;
      my $misoutput = "ERROR: bdp $bdp_id ($pdb_id)\n" ;

# Iterate through real chains.
      foreach my $k (0 .. $#{$real_chain_no}) {

# Find highest scoring entry in the sequence match matrix row (call max_ind(sequence match matrix row))
         my $t_bestmatch = max_ind($seq_match->[$k]) ;
         my $score = $seq_match->[$k][$t_bestmatch] ;

# Set the equivalent chain of this real chain to the highest scoring pdb chain.
         $equiv_ch->[$k] = $pdb_chain_no->[$t_bestmatch] ;
         $equiv_ch_id->[$k] = $pdb_chain_id->[$t_bestmatch] ;


# If the score is 0, set the mismatch flag to 1, as no equivalent PDB chain has been found.
         if ($score == 0) {
            $mismatch = 1;
            $misoutput .= "   REAL $real_chain_id->[$k]\t(".($k + 1).")\t$real_chain_seq->[$k]\n" ;
         }

# Otherwise, display chain equivalence information.
         else {
            my @outvals = ($bdp_id, $real_chain_no->[$k], $real_chain_id->[$k],
                           $pdb_id, $equiv_ch->[$k], $equiv_ch_id->[$k]) ;
            print join("\t", @outvals)."\n" ;
         }
      }

# If mismatch flag is set, display all real and PDB chain sequences to STDERR.
      if ($mismatch) {
         print STDERR $misoutput."\n" ;
         foreach my $l (0 .. $#{$pdb_chain_no}) {
            print STDERR  "   PDB $pdb_chain_id->[$l]\t(".($l +1).")\t$pdb_chain_seq->[$l]\n" ; }
      }
   }

   while (my $line = <STDIN>) {
      if ($line =~ /^#/) {next;}

      chomp $line;

      my ($bdp_id, $real_chain_no, $real_chain_id, $pdb_id,
          $pdb_chain_no, $pdb_chain_id) = split(/\t/, $line) ;

      my $whereclause  =  "bdp_id = $bdp_id" ;
      $whereclause .= " and real_chain_no = $real_chain_no" ;
 
      my $updateval = "pdb_chain_no = $pdb_chain_no" ;
      if ($pdb_chain_id ne '') {
         $updateval .= ", pdb_chain_id = \"$pdb_chain_id\"" ; }

# UPDATE bdp_chains SET pdb_chain_no = <pdb_chain_no> [, pdb_chain_id = <pdb_chain_id>] WHERE bdp_id = <bdp_id> and real_chain_no = <real_chain_no>

      my $chain_update_sql = "UPDATE bdp_chains SET $updateval ".
                              "WHERE $whereclause" ;
      $chain_update_sql = qq{$chain_update_sql} ;
      
      my $sth = $dbh->prepare($chain_update_sql) ;
      $sth->execute() ;
   }

}

=head2 max_ind (arrayref)

   Title:       max_ind()
   Function:    Return the index of the highest value in an array
   Args:        1 array reference
   Returns:     index of highest value in the array

=cut

sub max_ind {

   my $array = shift;
   my $max = $array->[0];
   my $maxind = 0;

   foreach my $j (1 .. $#{$array}) {
      if ($array->[$j] > $max) {
         $maxind = $j ;
         $max = $array->[$j];
      }
   }

   return($maxind) ;

}


#=head2 bdp_chains_equiv_import()
#
#   Title:       bdp_chains_equiv_import()
#   Function:    Import the chain equivalency table into pibase
#   Args:        none
#   Returns:     nothing
#   STDIN:       tabbed: bdp_id, real_chain_no, real_chain_id, pdb_id,
#                  pdb_chain_no, pdb_chain_id
#   Tables out:  bdp_chains
#
#=cut

sub OLDbdp_chains_equiv_import {

   require DBI;

# Connect to the pibase database.
   my ($dbh, $pibase) = pibase::connect_pibase() ;

   while (my $line = <STDIN>) {
      if ($line =~ /^#/) {next;}

      chomp $line;

      my ($bdp_id, $real_chain_no, $real_chain_id, $pdb_id,
          $pdb_chain_no, $pdb_chain_id) = split(/\t/, $line) ;

      my $whereclause  =  "bdp_id = $bdp_id" ;
      $whereclause .= " and real_chain_no = $real_chain_no" ;
 
      my $updateval = "pdb_chain_no = $pdb_chain_no" ;
      if ($pdb_chain_id ne '') {
         $updateval .= ", pdb_chain_id = \"$pdb_chain_id\"" ; }

# UPDATE bdp_chains SET pdb_chain_no = <pdb_chain_no> [, pdb_chain_id = <pdb_chain_id>] WHERE bdp_id = <bdp_id> and real_chain_no = <real_chain_no>

      my $chain_update_sql = "UPDATE bdp_chains SET $updateval ".
                              "WHERE $whereclause" ;
      $chain_update_sql = qq{$chain_update_sql} ;
      
      my $sth = $dbh->prepare($chain_update_sql) ;
      $sth->execute() ;
   }

}


=head2 call_chain_info()

   Title:       call_chain_info()
   Function:    Calculate chain information for list of pdb files.
   Args:        $_>{in_fn} - optional tabbed list of bdp_id, bdp_path to process
                $_>{pibase_specs} - optional pibase_specs structure
                $_>{cluster_fl} - if 1 will send to SGE cluster
                $_>{import_fl} - flag to import output into PIBASE
   Files out:   foreach bdp_id:
                o bdp_chains.<base pdb file name>.<hostname>.<XXXXX>.chaininfo.out
                     format per pibase::PDB::chains::chain_info()
   Returns:     import_status

=cut

sub call_chain_info {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }

# Create output file for chaininfo.

   my $bdpid2path = {};
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line =~ /^#/) {next;}
         chomp $line;
         my ($bdp_id, $bdp_path) = split(/\t/, $line) ;
         $bdpid2path->{$bdp_id} = $bdp_path ;
      }
      close(INF) ;
   } else {
      my ($dbh) = pibase::connect_pibase() ;
      $bdpid2path = pibase::mysql_hashload($dbh,
         "SELECT bdp_id, file_path FROM bdp_files") ;
   }


   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
# if this is master script, split input, recall self with
# in_fn specified, cluster_fl = 0, import_fl = 0
# send out, cluster run, and return merged

      print "* call_chain_info() ".localtime() if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{call_chain_info_in}, $temp_fn->{call_chain_info_in}) =
         tempfile("splits_call_chain_info_input.XXXXX");
      ($temp_fh->{call_chain_info_out}, $temp_fn->{call_chain_info_out}) =
         tempfile("splits_call_chain_info_SGEout_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{call_chain_info_out}) ;
      ($temp_fh->{call_chain_info_err}, $temp_fn->{call_chain_info_err}) =
         tempfile("splits_call_chain_info_SGEerr_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{call_chain_info_err}) ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdpid2path}) {
         print {$temp_fh->{call_chain_info_in}}
            join("\t", $bdp_id, $bdpid2path->{$bdp_id})."\n" ; }
      close($temp_fh->{call_chain_info_in}) ;

      my $split_dir = tempdir("splits_call_chain_info.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{call_chain_info_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.call_chain_info.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/call_chain_info/ ;

main() ;

sub main {

         pibase::data::calc::call_chain_info({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.call_chain_info.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.call_chain_info.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id}) ;
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{call_chain_info_out},
         err_fn => $temp_fn->{call_chain_info_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{call_chain_info_out},
           $temp_fn->{call_chain_info_out}) ;
      open(REALOUTF,">".$pibase_specs->{buildfiles}->{bdp_chains}) ;
      while (my $line = readline($temp_fh->{call_chain_info_out})) {
         if ($line =~ /^\#/) {next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{call_chain_info_out}) ;
      close(REALOUTF) ;

   } else { # if actual compute node, then just run it and print to stdout

      foreach my $bdp_id (keys %{$bdpid2path}) {
         print STDERR "now on: $bdp_id\n" ;
# get chain information (prints to stdout)
         my $bdp_path = $bdpid2path->{$bdp_id} ;
#         print STDERR "calling chain_info on $bdp_path\n";
         my $chaininfo = pibase::PDB::chains::chain_info({
            pdb_fn => $bdp_path,
            identifier => $bdp_id,
         }) ;

# If returned error, display to STDERR and go to next bdp file
         if (exists $chaininfo->{error_fl}) {
            print STDERR "ERROR: $bdp_id pibase::PDB::chains::chain_info: ".
               $chaininfo->{error_fl}."\n" ;
            next;
         }
      }

   }

   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_chains file (specified in specs) into pibase
      $import_status = pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{bdp_chains}}) ;
   }

   return $import_status ;

}


sub OLDcall_chain_info {

   my $in = shift ;

# Create output file for chaininfo.

#   my $hostname = hostname() ;
#   my $basename = $inputf ;
#   if ($inputf =~ /\//) {
#      ($basename) = ($inputf =~ /.*\/(.+)$/) ; }
#
#   my ($fh, $chainf) = tempfile("bdp_chains.$basename.$hostname.XXXXX",
#                                SUFFIX => '.chaininfo.out') ;

   my ($fh, $chainf) = tempfile("bdp_chains.XXXXX", SUFFIX => '.chaininfo.out');
   close($fh) ;

   open (INF, $in->{in_fn}) ;
   while (my $line= <INF>) {
      if ($line =~ /^#/) {next;}

      chomp $line;

      my ($bdp_id, $bdp_path) = split(/\t/, $line) ;
      print "NOW on $bdp_id\n" ;

# Get chain information.
      pibase::PDB::chains::chain_info($bdp_path, ">>$chainf", $bdp_id) ;

# Initialize the error flag.
      my $err = 0 ;

# If the chaininfo file is empty, set the error flag.
      if (-z $chainf) {
         $err = 1;
         print "ERROR $bdp_path: chaininfo calculate error: empty $chainf\n" ;
      }
   }
   close(INF) ;

   print "OUTFILE: ".$chainf."\n" ;

}

=head2 calc_subset_sequence()

   Title:       calc_subset_sequence()
   Function:    Determine amino acid sequences of all domains with entries
                  in subset_residues_tables
   Args:        none
   Returns:     nothing
   Tables in:   1. bdp_residues_tables
                2. subsets_residues_tables

   STDOUT:      1. bdp_id
                2. subset_id
                3. number of amino acids
                4. sequence - string of 1-letter aa abreviations

=cut

sub calc_subsets_sequence {

   my $in = shift;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }

   my $rescode = {
      'ALA' => 'A' ,
      'ARG' => 'R' ,
      'ASN' => 'N' ,
      'ASP' => 'D' ,
      'CYS' => 'C' ,
      'GLN' => 'Q' ,
      'GLU' => 'E' ,
      'GLY' => 'G' ,
      'HIS' => 'H' ,
      'HSD' => 'H' ,
      'ILE' => 'I' ,
      'LEU' => 'L' ,
      'LYS' => 'K' ,
      'MET' => 'M' ,
      'PHE' => 'F' ,
      'PRO' => 'P' ,
      'SER' => 'S' ,
      'THR' => 'T' ,
      'TRP' => 'W' ,
      'TYR' => 'Y' ,
      'VAL' => 'V',
      'UNK' => 'X',
      '  C' => 'c',
      '  G' => 'g',
      '  A' => 'a',
      '  T' => 't',
      '  U' => 'u',
      '  I' => 'i',
      ' +C' => 'c',
      ' +G' => 'g',
      ' +A' => 'a',
      ' +T' => 't',
      ' +U' => 'u',
      ' +I' => 'i',
      'MSE' => 'M',
      'MEX' => 'C',
   } ;

   my $bdpres_tables ;
   {
      my ($t_bdp, $t_restable_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM bdp_residues_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdpres_tables->{$t_bdp->[$j]} = $t_restable_fn->[$j] ; }
   }

   my $subsres_tables ;
   {
      my ($t_bdp, $t_restable_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM subsets_residues_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $subsres_tables->{$t_bdp->[$j]} = $t_restable_fn->[$j] ; }
   }

   my $bdp_list;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line =~ /^\#/) {next;}
         chomp $line; my $bdp_id = $line;
         $bdp_list->{$bdp_id}++ ;
      }
      close(INF) ;
   } else {
      my ($dbh) = pibase::connect_pibase() ;
      $bdp_list = pibase::mysql_hashload( $dbh,
         "SELECT distinct bdp_id, 1 FROM subsets WHERE bdp_id IS NOT NULL") ;
   }

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      print "* calc_subsets_sequence() ".localtime() if
         (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{calc_subsets_sequence_in},
       $temp_fn->{calc_subsets_sequence_in}) =
         tempfile("splits_calc_subsets_sequence_input.XXXXX");

      ($temp_fh->{calc_subsets_sequence_out},
       $temp_fn->{calc_subsets_sequence_out}) =
         tempfile("splits_calc_subsets_sequence_SGEout_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_subsets_sequence_out}) ;

      ($temp_fh->{calc_subsets_sequence_err},
       $temp_fn->{calc_subsets_sequence_err}) =
         tempfile("splits_calc_subsets_sequence_SGEerr_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_subsets_sequence_err}) ;

      foreach my $bdp (keys %{$bdp_list}) {
         print {$temp_fh->{calc_subsets_sequence_in}} $bdp."\n" ;}
      close($temp_fh->{calc_subsets_sequence_in}) ;

      my $split_dir = tempdir("splits_calc_subsets_sequence.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{calc_subsets_sequence_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.calc_subsets_sequence.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/calc_subsets_sequence/ ;

main() ;

sub main {

         pibase::data::calc::calc_subsets_sequence({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.calc_subsets_sequence.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.calc_subsets_sequence.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id}) ;
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{calc_subsets_sequence_out},
         err_fn => $temp_fn->{calc_subsets_sequence_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{calc_subsets_sequence_out},
           $temp_fn->{calc_subsets_sequence_out}) ;
      open(REALOUTF,">".$pibase_specs->{buildfiles}->{subsets_sequence}) ;

      while (my $line = readline($temp_fh->{calc_subsets_sequence_out})) {
         if ($line =~ /^\#/) { next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{calc_subsets_sequence_out}) ;
      close(REALOUTF) ;

   } else {

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         if (!-s $bdpres_tables->{$bdp_id}) {
            print STDERR "WARNING: (bdp $bdp_id) no entries in ".
                         "bdp_residues table $bdpres_tables->{$bdp_id}\n";
            next;
         }
   
         if (!-s $subsres_tables->{$bdp_id}) {
            print STDERR "WARNING: (bdp $bdp_id) no entries in ".
                         "subsets_residues table $subsres_tables->{$bdp_id}\n";
            next;
         }
   
         my ($b_resno, $b_chain_id, $b_resna) =
            pibase::rawselect_metatod( $bdpres_tables->{$bdp_id},
               "SELECT resno, chain_id, resna ".
               "FROM $bdpres_tables->{$bdp_id}") ;
   
         my $resno_2_na ;
         foreach my $j ( 0 .. $#{$b_resno}) {
            $resno_2_na->{$b_resno->[$j]."\n".$b_chain_id->[$j]} =
               $b_resna->[$j];}
   
         my ($sid, $res_ser, $resno, $chain_id) =
         pibase::rawselect_metatod( $subsres_tables->{$bdp_id},
            "SELECT subset_id, resno_serial, resno, chain_id ".
            "FROM $subsres_tables->{$bdp_id}" ) ;
   
         my $lastresser ;
         my $lastch ;
         my $seq ;
         my $length ;
         my $sid2chains ;
         foreach my $j ( 0 .. $#{$sid}) {
            my $concat = '';
            if (defined $lastresser->{$sid->[$j]}) {
               if ($chain_id->[$j] ne $lastch->{$sid->[$j]}) {
                  $concat = '/' ;
               } elsif ($lastresser->{$sid->[$j]} < ( $res_ser->[$j] - 1 )) {
                  $concat = '.' ;
               }
            }
            $seq->{$sid->[$j]} .= $concat.$rescode->{
                                   $resno_2_na->{
                                    $resno->[$j]."\n".$chain_id->[$j]}} ;
   
            $length->{$sid->[$j]}++ ;
            $lastresser->{$sid->[$j]} = $res_ser->[$j] ;
            $lastch->{$sid->[$j]} = $chain_id->[$j] ;
            $sid2chains->{$sid->[$j]}->{$chain_id->[$j]}++ ;
         }
   
         foreach my $sid ( sort keys %{$seq}) {
            my $numchains = keys %{$sid2chains->{$sid}} ;
            if (!defined $seq->{$sid} || $seq->{$sid} eq '') {
              print STDERR "ERROR: $sid (bdp_id $bdp_id ) undefined sequence\n";
               next;
            }
            my @outvals = ($bdp_id, $sid, $numchains,
                           $length->{$sid}, $seq->{$sid}) ;
            print join("\t", @outvals)."\n" ;
         }
      }
   }

   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_residues_tables file (specified in specs) into pibase
      $import_status->{subsets_sequence} = pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{subsets_sequence}}) ;
   }

   return $import_status ;

}


=head2 calc_subset_sasa()

   Title:       calc_subset_sasa()
   Function:    Computes SASA for all (or specified) PIBASE domains
   Args:        none
   Returns:     nothing

=cut

sub calc_subsets_sasa {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }
   my $modeller_bin = $pibase_specs->{binaries}->{modeller} ;
   my $move_thresh = 100 ;

   my $sid2fn ;
   my $sid2bdp_id ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line =~ /^\#/) {next;}
         chomp $line;
         my ($bdp_id, $subset_id, $subset_fn) = split(/\t/, $line) ;
         $sid2fn->{$subset_id} = $subset_fn;
         $sid2bdp_id->{$subset_id} = $bdp_id ;
      }
      close(INF) ;
   } else {
      my ($dbh) = pibase::connect_pibase() ;
      $sid2fn = pibase::mysql_hashload($dbh,
         "SELECT subset_id, file_path FROM subsets_files") ;
      $sid2bdp_id = pibase::mysql_hashload( $dbh,
         "SELECT subset_id, bdp_id FROM subsets WHERE bdp_id IS NOT NULL") ;
   }

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      print "* calc_subsets_sasa() ".localtime() if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{calc_subsets_sasa_in}, $temp_fn->{calc_subsets_sasa_in}) =
         tempfile("splits_calc_subsets_sasa_input.XXXXX");
      ($temp_fh->{calc_subsets_sasa_out}, $temp_fn->{calc_subsets_sasa_out}) =
         tempfile("splits_calc_subsets_sasa_SGEout_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{calc_subsets_sasa_out}) ;
      ($temp_fh->{calc_subsets_sasa_err}, $temp_fn->{calc_subsets_sasa_err}) =
         tempfile("splits_calc_subsets_sasa_SGEerr_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{calc_subsets_sasa_err}) ;

      foreach my $sid (sort keys %{$sid2fn}) {
         my @outvals = ($sid2bdp_id->{$sid}, $sid, $sid2fn->{$sid}) ;
         print {$temp_fh->{calc_subsets_sasa_in}} join("\t",@outvals)."\n" ;
      }

      close($temp_fh->{calc_subsets_sasa_in}) ;

      my $split_dir = tempdir("splits_calc_subsets_sasa.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{calc_subsets_sasa_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.calc_subsets_sasa.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/calc_subsets_sasa/ ;

main() ;

sub main {

         pibase::data::calc::calc_subsets_sasa({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.calc_subsets_sasa.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.calc_subsets_sasa.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id}) ;
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{calc_subsets_sasa_out},
         err_fn => $temp_fn->{calc_subsets_sasa_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{calc_subsets_sasa_out},
           $temp_fn->{calc_subsets_sasa_out}) ;
      open(REALOUTF,">".$pibase_specs->{buildfiles}->{subsets_sasa}) ;

      while (my $line = readline($temp_fh->{calc_subsets_sasa_out})) {
         if ($line =~ /^\#/) { next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{calc_subsets_sasa_out}) ;
      close(REALOUTF) ;

   } else {

      my $tempdir = tempdir(CLEANUP => 1) ;
      chdir $tempdir ;

      foreach my $subset_id (sort keys %{$sid2fn}) {
         my $bdp_id = $sid2bdp_id->{$subset_id} ;
   
         my $subset_sasa = pibase::modeller::calc_sasa({
               pdb_fn => $sid2fn->{$subset_id},
               surftyp => 2,
               modeller_bin => $modeller_bin
         }) ;
   
         if ($#{$subset_sasa->{error_fl}} >= 0 ) {
            foreach my $j ( 0 .. $#{$subset_sasa->{error_fl}}) {
               print STDERR "ERROR: subset ($bdp_id $subset_id): ".
                           "calc_sasa() $subset_sasa->{error_fl}->[$j]\n";
            }
            next;
         }

         my @sasa_outvals = ($bdp_id, $subset_id, 
            sprintf("%.3f", $subset_sasa->{full_sasa}->{all}),
            sprintf("%.3f", $subset_sasa->{full_sasa}->{sc}),
            sprintf("%.3f", $subset_sasa->{full_sasa}->{mc}),
            sprintf("%.3f", $subset_sasa->{atm_sasa}->{p}),
            sprintf("%.3f", $subset_sasa->{atm_sasa}->{nonp})
         ) ;
         print join("\t", @sasa_outvals)."\n" ;
      }
   }

   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_residues_tables file (specified in specs) into pibase
      $import_status->{subsets_sasa} = pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{subsets_sasa}}) ;
   }

   return $import_status ;

}


=head2 calc_interface_dsasa()

   Title:       calc_interface_dsasa()
   Function:    Computes dSASA for all (or specified) PIBASE interfaces

=cut

sub calc_interface_dsasa {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }
   my $modeller_bin = $pibase_specs->{binaries}->{modeller} ;
   my $move_thresh = 100 ;

   my $sid2sasa ;
   {
      my ($t_sid, $t_sasa_all, $t_sasa_sc, $t_sasa_mc,
          $t_sasa_p, $t_sasa_nonp) = pibase::rawselect_tod(
         "SELECT subset_id, sasa_all, sasa_sc, sasa_mc, ".
         "sasa_polar, sasa_nonpolar FROM subsets_sasa") ;
      foreach my $j ( 0 .. $#{$t_sid}) {
         $sid2sasa->{$t_sid->[$j]} = {
            all => $t_sasa_all->[$j],
            sc => $t_sasa_sc->[$j],
            mc => $t_sasa_mc->[$j],
            p => $t_sasa_p->[$j],
            nonp => $t_sasa_nonp->[$j],
         } ;
      }
   }

   my $sid2fn ;
   {
      my ($t_sid, $t_fn) = pibase::rawselect_tod(
         "SELECT subset_id, file_path FROM subsets_files") ;
      foreach my $j ( 0 .. $#{$t_sid}) {
         $sid2fn->{$t_sid->[$j]} = $t_fn->[$j] ; }
   }

   my $sid12_2_bdp ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line =~ /^\#/) {next;}
         chomp $line;
         my ($bdp_id, $sid1, $sid2) = split(/\t/, $line) ;
         $sid12_2_bdp->{$sid1."\t".$sid2} = $bdp_id ;
      }
      close(INF) ;
   } else {
      my ($bdp_id, $sid1, $sid2) = pibase::rawselect_tod(
         "SELECT bdp_id, subset_id_1, subset_id_2 FROM intersubset_contacts") ;
      foreach my $j ( 0 .. $#{$sid1}) {
         $sid12_2_bdp->{$sid1->[$j]."\t".$sid2->[$j]} = $bdp_id->[$j] ; }
   }

   
   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      print "* calc_interface_dsasa() ".localtime() if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{calc_interface_dsasa_in},
       $temp_fn->{calc_interface_dsasa_in}) =
         tempfile("splits_calc_interface_dsasa_input.XXXXX");

      ($temp_fh->{calc_interface_dsasa_out},
       $temp_fn->{calc_interface_dsasa_out}) =
         tempfile("splits_calc_interface_dsasa_SGEout_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_dsasa_out}) ;

      ($temp_fh->{calc_interface_dsasa_err},
       $temp_fn->{calc_interface_dsasa_err}) =
         tempfile("splits_calc_interface_dsasa_SGEerr_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_dsasa_err}) ;

      foreach my $sid12 (sort keys %{$sid12_2_bdp}) {
         my ($sid1, $sid2) = split(/\t/, $sid12) ;
         my @outvals = ($sid12_2_bdp->{$sid12}, $sid1, $sid2) ;
         print {$temp_fh->{calc_interface_dsasa_in}} join("\t",@outvals)."\n" ;
      }
      close($temp_fh->{calc_interface_dsasa_in}) ;

      my $split_dir = tempdir("splits_calc_interface_dsasa.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{calc_interface_dsasa_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.calc_interface_dsasa.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/calc_interface_dsasa/ ;

main() ;

sub main {

         pibase::data::calc::calc_interface_dsasa({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.calc_interface_dsasa.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.calc_interface_dsasa.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id}) ;
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{calc_interface_dsasa_out},
         err_fn => $temp_fn->{calc_interface_dsasa_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{calc_interface_dsasa_out},
           $temp_fn->{calc_interface_dsasa_out}) ;
      open(REALOUTF,">".$pibase_specs->{buildfiles}->{interface_sasa});

      while (my $line = readline($temp_fh->{calc_interface_dsasa_out})){
         if ($line =~ /^\#/) { next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{calc_interface_dsasa_out}) ;
      close(REALOUTF) ;

   } else {

      my $tempdir = tempdir(CLEANUP => 1) ;
      chdir $tempdir ;

      foreach my $sid12 (sort keys %{$sid12_2_bdp}) {
#         print STDERR "NOW ON: bdp $bdp interface $sid1 -- $sid2\n" ;
         my $bdp_id = $sid12_2_bdp->{$sid12} ;
         my ($sid1, $sid2) = split(/\t/, $sid12) ;
         my $sid1_fn = $sid2fn->{$sid1} ;
         my $sid2_fn = $sid2fn->{$sid2} ;

         if (!exists $sid2sasa->{$sid1}) {
            print STDERR "ERROR $sid1 - $sid2 interface_dsasa() : ".
                         "$sid1 SASA missing\n" ; next; }

         if (!exists $sid2sasa->{$sid2}) {
            print STDERR "ERROR $sid1 - $sid2 interface_dsasa() : ".
                         "$sid2 SASA missing\n" ; next; }

         if (!defined $sid1_fn || !-s $sid1_fn) {
            print STDERR "ERROR $sid1 - $sid2 interface_dsasa() : ".
                         "$sid1 subsets file missing\n" ; next; }

         if (!defined $sid2_fn || !-s $sid2_fn) {
            print STDERR "ERROR $sid1 - $sid2 interface_dsasa() : ".
                         "$sid2 subsets file missing\n" ; next; }

         my ($pairpdb_fh, $pairpdb_fn) =
            tempfile("pairpdb.$bdp_id.XXXXXX", SUFFIX=>"pdb") ;
         close($pairpdb_fh) ;
         foreach my $t_sid_fn ($sid1_fn, $sid2_fn) {
            if ($t_sid_fn =~ /\.gz$/) {
               system($pibase_specs->{binaries}->{zcat}.
                      " $t_sid_fn >> $pairpdb_fn") ;
            } else {
               system("cat $t_sid_fn >> $pairpdb_fn") ;
            }
         }
#         system("cat $sid1_fn $sid2_fn > $pairpdb_fn") ;

         if (!-s $pairpdb_fn) {
            print STDERR "ERROR: couldnt make complex pdb file: $sid1, $sid2\n";
            next; }
   
         my $sid12_sasa = pibase::modeller::calc_sasa({
            pdb_fn => $pairpdb_fn,
            surftyp => 2,
            modeller_bin => $modeller_bin
         }) ;
         unlink $pairpdb_fn ;

         if ($#{$sid12_sasa->{error_fl}} >= 0 ) {
            foreach my $j ( 0 .. $#{$sid12_sasa->{error_fl}}) {
               print STDERR "ERROR: interface ($bdp_id $sid1 $sid2): ".
                           "calc_sasa() $sid12_sasa->{error_fl}->[$j]\n";
            }
            next;
         }

         my $dsasa ;
         $dsasa->{all} = $sid2sasa->{$sid1}->{all} + $sid2sasa->{$sid2}->{all} -
                         $sid12_sasa->{full_sasa}->{all} ;

         $dsasa->{sc} = $sid2sasa->{$sid1}->{sc} + $sid2sasa->{$sid2}->{sc} -
                         $sid12_sasa->{full_sasa}->{sc} ;

         $dsasa->{mc} = $sid2sasa->{$sid1}->{mc} + $sid2sasa->{$sid2}->{mc} -
                         $sid12_sasa->{full_sasa}->{mc} ;

         $dsasa->{p} = $sid2sasa->{$sid1}->{p} +
                       $sid2sasa->{$sid2}->{p} -
                       $sid12_sasa->{atm_sasa}->{p} ;
   
         $dsasa->{nonp} = $sid2sasa->{$sid1}->{nonp} +
                          $sid2sasa->{$sid2}->{nonp} -
                          $sid12_sasa->{atm_sasa}->{nonp} ;

         my @dsasa_outvals = ($bdp_id, $sid1, $sid2,
            sprintf("%.3f", $dsasa->{all}),
            sprintf("%.3f", $dsasa->{sc}),
            sprintf("%.3f", $dsasa->{mc}),
            sprintf("%.3f", $dsasa->{p}),
            sprintf("%.3f", $dsasa->{nonp})
         ) ;
         print join("\t", @dsasa_outvals)."\n" ;
      }
   }

   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_residues_tables file (specified in specs) into pibase
      $import_status->{subsets_sasa} = pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{interface_sasa}}) ;
   }

   return $import_status ;

}


=head2 calc_interface_secstrx()

   Title:       calc_interface_secstrx()
   Function:    Computes secondary structure for all (or specified) PIBASE domains

=cut

sub calc_interface_secstrx {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }
   my $modeller_bin = $pibase_specs->{binaries}->{modeller} ;
   my $move_thresh = 100 ;

   my $bdp2sid12 ;
   {
      my ($t_bdp, $t_sid1, $t_sid2) = pibase::rawselect_tod(
         "SELECT bdp_id, subset_id_1, subset_id_2 FROM intersubset_contacts") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2sid12->{$t_bdp->[$j]} = $t_sid1->[$j]."\t".$t_sid2->[$j] ;
      }
   }

   my $bdp2secstrx_fn;
   {
      my ($t_bdp, $t_secstrx_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM bdp_secstrx_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2secstrx_fn->{$t_bdp->[$j]} = $t_secstrx_fn->[$j]  ; }
   }

   my $bdp2intcon_fn;
   {
      my ($t_bdp, $t_intcon_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM interface_contacts_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2intcon_fn->{$t_bdp->[$j]} = $t_intcon_fn->[$j]  ; }
   }

   my $bdp_list ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line =~ /^\#/) {next;}
         chomp $line;
         my $bdp_id = $line ;
         $bdp_list->{$line}++  ;
      }
      close(INF) ;
   } else {
      my ($t_bdp) = pibase::rawselect_tod(
         "SELECT bdp_id FROM intersubset_contacts") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp_list->{$t_bdp->[$j]}++ ;}
   }

   
   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      print "* calc_interface_secstrx() ".localtime()
         if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{calc_interface_secstrx_in},
       $temp_fn->{calc_interface_secstrx_in}) =
         tempfile("splits_calc_interface_secstrx_input.XXXXX");

      ($temp_fh->{calc_interface_secstrx_out},
       $temp_fn->{calc_interface_secstrx_out}) =
         tempfile("splits_calc_interface_secstrx_SGEout_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_secstrx_out}) ;

      ($temp_fh->{calc_interface_secstrx_err},
       $temp_fn->{calc_interface_secstrx_err}) =
         tempfile("splits_calc_interface_secstrx_SGEerr_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_secstrx_err}) ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         print {$temp_fh->{calc_interface_secstrx_in}} $bdp_id."\n" ; }
      close($temp_fh->{calc_interface_secstrx_in}) ;

      my $split_dir = tempdir("splits_calc_interface_secstrx.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{calc_interface_secstrx_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.calc_interface_secstrx.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/calc_interface_secstrx/ ;

main() ;

sub main {

         pibase::data::calc::calc_interface_secstrx({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.calc_interface_secstrx.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.calc_interface_secstrx.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id}) ;
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{calc_interface_secstrx_out},
         err_fn => $temp_fn->{calc_interface_secstrx_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{calc_interface_secstrx_out},
           $temp_fn->{calc_interface_secstrx_out}) ;
      open(REALOUTF,
         ">".$pibase_specs->{buildfiles}->{interface_secstrx_tables});

      while (my $line = readline($temp_fh->{calc_interface_secstrx_out})){
         if ($line =~ /^\#/) { next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{calc_interface_secstrx_out}) ;
      close(REALOUTF) ;

   } else {

      my (@movethese, @moveto) ;
      my $movethreshold = 100 ;

      my $tempdir = tempdir(CLEANUP => 1) ;
      chdir $tempdir ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         print STDERR "now on bdp $bdp_id\n" ;

         if (! exists $bdp2intcon_fn->{$bdp_id} ) {
            print STDERR "ERROR: interface_contacts_tables not found ".
                         "for bdp $bdp_id\n" ;
            next;
         }
   
         if (!exists $bdp2secstrx_fn->{$bdp_id}) {
            print STDERR "ERROR: bdp_secstrx_tables not found ".
                         "for bdp $bdp_id\n" ;
            next ;
         }
   
   
#   Read in full residue listing for this bdp entry
         my $res2sse;
         my $sse_defs;
         {
            my ($t_resno, $t_chain_id, $t_sseid, $t_sse) =
               pibase::rawselect_metatod( $bdp2secstrx_fn->{$bdp_id},
               "SELECT resno, chain_id, sse_id, sse ".
               "FROM $bdp2secstrx_fn->{$bdp_id}") ;

            foreach my $j ( 0 .. $#{$t_resno}) {
               my $ressig = $t_resno->[$j]."\n".$t_chain_id->[$j] ;
               $res2sse->{$ressig} = $t_sseid->[$j] ;
               $sse_defs->{$t_sseid->[$j]} = $t_sse->[$j] ;
            }
         }

         my $out_fn = "interface_secstrx_$bdp_id.out" ;
         my $out_fh  ;
         open($out_fh, ">$out_fn") ;
   
# Get all interface residue contacts
         my ( $subset_id_1, $subset_id_2,
              $chain_id_1, $resno_1, $chain_id_2, $resno_2 ) =
            pibase::rawselect_metatod($bdp2intcon_fn->{$bdp_id},
      "SELECT subset_id_1, subset_id_2, chain_id_1, resno_1, chain_id_2, resno_2 FROM $bdp2intcon_fn->{$bdp_id}") ;
   
#Iterate over intersubset contacts

   #      my $res_content ;
         my $sse_content ;
         foreach my $j ( 0 .. $#{$subset_id_1} ) {
   
            my $sub1 = $subset_id_1->[$j] ;
            my $sub2 = $subset_id_2->[$j] ;
            my $sub12 = $sub1."\n".$sub2 ;
   
# Translate residue contacts into sse contacts
   
            if ((!defined $chain_id_1->[$j]) || ($chain_id_1->[$j] eq '')) {
               $chain_id_1->[$j] = ' ' ; }
   
            if ((!defined $chain_id_2->[$j]) || ($chain_id_2->[$j] eq '')) {
               $chain_id_2->[$j] = ' ' ; }
   
            my $ressig1 = $resno_1->[$j]."\n".$chain_id_1->[$j] ;
            my $ressig2 = $resno_2->[$j]."\n".$chain_id_2->[$j] ;
            my $ressig12 = $ressig1."\n".$ressig2 ;

#$res_content->{$sub12}->{sub1}->{$ressig1}++ ;
#$res_content->{$sub12}->{sub2}->{$ressig2}++ ;
#$res_content->{$sub12}->{sub12}->{$ressig12}++ ;
   
            if ((! exists $res2sse->{$ressig1}) ||
                 ($res2sse->{$ressig1} eq ''))  {
               my $t = "WARNING ($bdp_id): no SSE assignment for $ressig1" ;
               $t =~ s/\n/:/g ;
               print STDERR $t."\n" ;
               next;
   
            } elsif ((! exists $res2sse->{$ressig2}) ||
                      ($res2sse->{$ressig2} eq ''))  {
   
               my $t = "WARNING ($bdp_id): no SSE assignment for $ressig2" ;
               $t =~ s/\n/:/g ;
               print STDERR $t."\n" ;
               next;
            }
   
            my $sse1 = $res2sse->{$ressig1} ;
            my $sse2 = $res2sse->{$ressig2} ;
            
            $sse_content->{$sub12}->{"sub1"}->{$sse1}->{$ressig1}++ ;
            $sse_content->{$sub12}->{"sub2"}->{$sse2}->{$ressig2}++ ;
   
         }
   
# Iterate over intersubset contacts
   
         foreach my $sub12 (keys %{$sse_content}) {
            my ($sub1, $sub2) = split("\n", $sub12) ;
   
            foreach my $sse1 (keys %{$sse_content->{$sub12}->{sub1}}) {
               my @outvals ;
               push @outvals, ($bdp_id, $sub1, $sub2, $sub1, $sse1) ;

               my $numres = keys %{$sse_content->{$sub12}->{sub1}->{$sse1}} ;
               push @outvals, $numres ;
               push @outvals, $sse_defs->{$sse1} ;
               print $out_fh join("\t", @outvals)."\n" ;
            }
   
            foreach my $sse2 (keys %{$sse_content->{$sub12}->{sub2}}) {
               my @outvals ;
               push @outvals, ($bdp_id, $sub1, $sub2, $sub2, $sse2) ;

               my $numres = keys %{$sse_content->{$sub12}->{sub2}->{$sse2}} ;
               
               push @outvals, $numres ;
               push @outvals, $sse_defs->{$sse2} ;
               print $out_fh join("\t", @outvals)."\n" ;
            }
   
         }
         close($out_fh) ;
   
         if (-s $out_fn) {
            my $t_dir = $pibase_specs->{metatod_dir}->{interface_secstrx}.'/'.
                           POSIX::floor($bdp_id / 1000) ;
            if (!-s $t_dir) { mkpath $t_dir; }
            push @movethese, $out_fn ;
            push @moveto, $t_dir ;

            my $t_filena = $t_dir.'/'.$out_fn ;
            print join("\t", $bdp_id, '\N', $t_filena)."\n" ;
         } else {
            print STDERR "ERROR: bdp $bdp_id output is empty\n" ;
            unlink $out_fn ;
         }
   
         if ($#movethese >= $movethreshold) {
            foreach my $j (0 .. $#movethese) {
               my $t_fn = $movethese[$j] ;
               my $t_dir = $moveto[$j] ;
               pibase::safe_move($t_fn, $t_dir);}
            @movethese = () ; @moveto = () ;
         }
      }

      foreach my $j (0 .. $#movethese) {
         my $t_fn = $movethese[$j] ;
         my $t_dir = $moveto[$j] ;
         pibase::safe_move($t_fn, $t_dir);
      }
   }


   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_residues_tables file (specified in specs) into pibase
      $import_status->{interface_secstrx_tables} = pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{interface_secstrx_tables}}) ;
   }

   return $import_status ;

}


=head2 calc_interface_secstrx_basic()

   Title:       calc_interface_secstrx_basic()
   Function:    Computes basic secondary structure for all (or specified)
                  PIBASE domains (H|B|- instead of full 6/7-letter DSSP code)

=cut

sub calc_interface_bs_secstrx_basic {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }
   my $modeller_bin = $pibase_specs->{binaries}->{modeller} ;
   my $move_thresh = 100 ;

   my $param_dist_cutoff = 6.05 ;

   my $bdp2sid12 ;
   {
      my ($t_bdp, $t_sid1, $t_sid2) = pibase::rawselect_tod(
         "SELECT bdp_id, subset_id_1, subset_id_2 FROM intersubset_contacts") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2sid12->{$t_bdp->[$j]} = $t_sid1->[$j]."\t".$t_sid2->[$j] ;
      }
   }


   my $bdp2secstrx_fn;
   {
      my ($t_bdp, $t_secstrx_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM bdp_secstrx_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2secstrx_fn->{$t_bdp->[$j]} = $t_secstrx_fn->[$j]  ; }
   }

   my $bdp2intcon_fn;
   {
      my ($t_bdp, $t_intcon_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM interface_contacts_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2intcon_fn->{$t_bdp->[$j]} = $t_intcon_fn->[$j]  ; }
   }

   my $bdp2fn;
   {
      my ($t_bdp, $t_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, file_path FROM bdp_files") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2fn->{$t_bdp->[$j]} = $t_fn->[$j]  ; }
   }

   my $bdp_list ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line =~ /^\#/) {next;}
         chomp $line;
         my $bdp_id = $line ;
         $bdp_list->{$line}++  ;
      }
      close(INF) ;
   } else {
      my ($t_bdp) = pibase::rawselect_tod(
         "SELECT bdp_id FROM intersubset_contacts") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp_list->{$t_bdp->[$j]}++ ;}
   }

   
   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      print "* calc_interface_bs_secstrx_basic() ".localtime()
         if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{calc_interface_bs_secstrx_basic_in},
       $temp_fn->{calc_interface_bs_secstrx_basic_in}) =
         tempfile("splits_calc_interface_bs_secstrx_basic_input.XXXXX");

      ($temp_fh->{calc_interface_bs_secstrx_basic_out},
       $temp_fn->{calc_interface_bs_secstrx_basic_out}) =
         tempfile("splits_calc_interface_bs_secstrx_basic_SGEout_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_bs_secstrx_basic_out}) ;

      ($temp_fh->{calc_interface_bs_secstrx_basic_err},
       $temp_fn->{calc_interface_bs_secstrx_basic_err}) =
         tempfile("splits_calc_interface_bs_secstrx_basic_SGEerr_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_bs_secstrx_basic_err}) ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         print {$temp_fh->{calc_interface_bs_secstrx_basic_in}} $bdp_id."\n" ; }
      close($temp_fh->{calc_interface_bs_secstrx_basic_in}) ;

      my $split_dir = tempdir("splits_calc_interface_bs_secstrx_basic.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{calc_interface_bs_secstrx_basic_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.calc_interface_bs_secstrx_basic.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/calc_interface_bs_secstrx_basic/ ;

main() ;

sub main {

         pibase::data::calc::calc_interface_bs_secstrx_basic({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.calc_interface_bs_secstrx_basic.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.calc_interface_bs_secstrx_basic.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id}) ;
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{calc_interface_bs_secstrx_basic_out},
         err_fn => $temp_fn->{calc_interface_bs_secstrx_basic_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{calc_interface_bs_secstrx_basic_out},
           $temp_fn->{calc_interface_bs_secstrx_basic_out}) ;

      open(REALOUTF_ISBC,
 ">".$pibase_specs->{buildfiles}->{interface_secstrx_basic_contacts_tables});

      open(REALOUTF_BSBC,
 ">".$pibase_specs->{buildfiles}->{bindingsite_secstrx_basic_contacts_tables});

      open(REALOUTF_BSC,
 ">".$pibase_specs->{buildfiles}->{bindingsite_contacts_tables});

      while (my $line =
       readline($temp_fh->{calc_interface_bs_secstrx_basic_out}) ) {
         if ($line =~ /^\#/) {
            next;
         } elsif ($line =~ /^ISBC/) {
            $line =~ s/^ISBC\t// ;
            print REALOUTF_ISBC $line ;
         } elsif ($line =~ /^BSBC/) {
            $line =~ s/^BSBC\t// ;
            print REALOUTF_BSBC $line ;
         } elsif ($line =~ /^BSC/) {
            $line =~ s/^BSC\t// ;
            print REALOUTF_BSC $line ;
         }
      }
      close($temp_fh->{calc_interface_bs_secstrx_basic_out}) ;
      close(REALOUTF_ISBC) ;
      close(REALOUTF_BSBC) ;
      close(REALOUTF_BSC) ;

   } else {

      my (@movethese, @moveto) ;
      my $movethreshold = 100 ;

      my $tempdir = tempdir(CLEANUP => 1) ;
      chdir $tempdir ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         print STDERR "now on bdp $bdp_id\n" ;
   
         if (! exists $bdp2intcon_fn->{$bdp_id} ) {
            print STDERR "ERROR: interface_contacts_tables not found ".
                         "for bdp $bdp_id\n" ;
            next;
         }
   
         if (!exists $bdp2secstrx_fn->{$bdp_id}) {
            print STDERR "ERROR: bdp_secstrx_tables not found ".
                         "for bdp $bdp_id\n" ;
            next ;
         }

         my $res2sse;
         my $sse_defs;
         {
            my ($t_resno, $t_chain_id, $t_sseid, $t_sse) =
               pibase::rawselect_metatod( $bdp2secstrx_fn->{$bdp_id},
               "SELECT resno, chain_id, sse_id, sse_basic ".
               "FROM $bdp2secstrx_fn->{$bdp_id}") ;

            foreach my $j ( 0 .. $#{$t_resno}) {
               my $ressig = $t_resno->[$j]."\n".$t_chain_id->[$j] ;
               $res2sse->{$ressig} = $t_sseid->[$j] ;
               $sse_defs->{$t_sseid->[$j]} = $t_sse->[$j] ;
            }
         }
   
         my ($out_fh, $out_fn) ;
   
         $out_fn->{int_sse_basic_contacts} = 
            "interface_secstrx_basic_contacts_$bdp_id.out" ;
         open ($out_fh->{int_sse_basic_contacts},
           ">".$out_fn->{int_sse_basic_contacts}) ;
   
#   Get all interface residue contacts
   
         my ( $subset_id_1, $subset_id_2,
              $chain_id_1, $resno_1, $chain_id_2, $resno_2 ) =
            pibase::rawselect_metatod($bdp2intcon_fn->{$bdp_id},
      "SELECT subset_id_1, subset_id_2, chain_id_1, resno_1, chain_id_2, resno_2 FROM $bdp2intcon_fn->{$bdp_id}") ;
   
   
         my ($sse_contacts, $sse_halves) ;
 
#   Iterate over intersubset contacts
   
         my $intersubset_contacts ;
         my $bs_res ;
         my $res2sub ;
         foreach my $j ( 0 .. $#{$subset_id_1} ) {

            my $sub1 = $subset_id_1->[$j] ;
            my $sub2 = $subset_id_2->[$j] ;
            my $sub12 = $sub1."\n".$sub2 ;
   
#   Translate residue contacts into sse contacts
   
            if ((!defined $chain_id_1->[$j]) || ($chain_id_1->[$j] eq '')) {
               $chain_id_1->[$j] = ' ' ; }
   
            if ((!defined $chain_id_2->[$j]) || ($chain_id_2->[$j] eq '')) {
               $chain_id_2->[$j] = ' ' ; }
   
            my $ressig1 = $resno_1->[$j]."\n".$chain_id_1->[$j] ;
            my $ressig2 = $resno_2->[$j]."\n".$chain_id_2->[$j] ;
   
            $bs_res->{$sub12}->[0]->{$ressig1}++;
            $bs_res->{$sub12}->[1]->{$ressig2}++;
            $res2sub->{$ressig1}->{$sub1}++ ;
            $res2sub->{$ressig2}->{$sub2}++ ;
   
            if ((! exists $res2sse->{$ressig1}) ||
                ($res2sse->{$ressig1} eq ''))  {
               my $t = "WARNING: (bdp $bdp_id) no SSE assignment for $ressig1" ;
               $t =~ s/\n/:/g ;
               print STDERR $t."\n" ;
               next;
            } elsif ((! exists $res2sse->{$ressig2}) ||
                     ($res2sse->{$ressig2} eq ''))  {
   
               my $t = "WARNING: (bdp $bdp_id) no SSE assignment for $ressig2" ;
               $t =~ s/\n/:/g ;
               print STDERR $t."\n" ;
               next;
            }
   
            my $ressig12 = $ressig1."\n".$ressig2 ;
            my $sse1 = $res2sse->{$ressig1} ;
            my $sse2 = $res2sse->{$ressig2} ;
            my $sse12 = $sse1."\n".$sse2 ;
   
   
            if (! exists $sse_contacts->{$sub12} ) {
               push @{$intersubset_contacts}, $sub12 ; }
   
            $sse_contacts->{$sub12}->{$sse12}->{$ressig12}++ ;
            $sse_halves->{$sub12}->{$sse12}->{sub1}->{$ressig1}++ ;
            $sse_halves->{$sub12}->{$sse12}->{sub2}->{$ressig2}++ ;
   
         }
   
   
         foreach my $sub12 (@{$intersubset_contacts}) {
            my ($sub1, $sub2) = split("\n", $sub12) ;
            foreach my $sse12 (keys %{$sse_contacts->{$sub12}}) {
               my ($sse1, $sse2) = split("\n", $sse12) ;
               my @outvals ;
               push @outvals, ($bdp_id, $sub1, $sub2, $sse1, $sse2) ;

               my $a = keys %{$sse_halves->{$sub12}->{$sse12}->{sub1}} ;
               my $b = keys %{$sse_halves->{$sub12}->{$sse12}->{sub2}} ;
               my $c = keys %{$sse_contacts->{$sub12}->{$sse12}} ;
               
               push @outvals, $a, $b, $c ;
               push @outvals, $sse_defs->{$sse1}, $sse_defs->{$sse2} ;
               print {$out_fh->{int_sse_basic_contacts}}
                  join("\t", @outvals)."\n" ;
            }
         }
         close($out_fh->{int_sse_basic_contacts}) ;
         system("gzip ".$out_fn->{int_sse_basic_contacts}) ;
         $out_fn->{int_sse_basic_contacts} .= '.gz' ;
   
         if (-s $out_fn->{int_sse_basic_contacts}) {
            push @movethese, $out_fn->{int_sse_basic_contacts} ;

            my $isbc_depositdir =
              $pibase_specs->{metatod_dir}->{interface_secstrx_basic_contacts}.
              '/'.POSIX::floor($bdp_id / 1000) ;

            push @moveto, $isbc_depositdir ;

            my $t_filena = $isbc_depositdir."/".
                        $out_fn->{int_sse_basic_contacts} ;
            my @outvals = ("ISBC", $bdp_id, '\N', $t_filena) ;
            print join("\t", @outvals)."\n" ;
         } else {
            print STDERR "ERROR: (bdp $bdp_id) ".
               "interface_secstrx_basic_contacts file is empty\n" ;
            unlink $out_fn->{int_sse_basic_contacts} ;
         }
   
# iterate over all conatcts and see if they are intra-binding site

# have to calculate this on the fly. - subset pdbs are already generated,
#  just split them in two and load it.

# fpd080617_1722 changed from interatomic_contacts_tables to 
#   on the fly calculation
# originally:
#         my $tables ;
#         $tables->{interatomic_contacts}->{sourcefile} =
#            $all_interatomic_contacts_files->{$bdp_id} ;
#         my ($contacts_fh) = raw_contacts_select(
#            $tables->{interatomic_contacts}->{sourcefile},
#               'SELECT chain_id_1, resno_1, resna_1, atomna_1, '.
#               'chain_id_2, resno_2, resna_2, atomna_2, distance' ,
#               {maxdist => $param->{dist_cutoff}} ) ;
#
#            my ($chain_id_1, $resno_1, $resna_1, $atomna_1,
#                $chain_id_2, $resno_2, $resna_2, $atomna_2,
#                $dist) = split(/\t/, $line) ;

#just call back _interface

#HERENOW
         my $kdcont_out =
            pibase::calc::interfaces::_interface_detect_calc__calc_res_pairs({
               radius => $param_dist_cutoff,
               compress => 1,
               bdp_path => $bdp2fn->{$bdp_id}
            }) ;
         my $kdfield2no = $kdcont_out->{field2no} ;
         my $contacts_fh ;
         open($contacts_fh, $kdcont_out->{contacts_fn}) ;

         my $respair_dup ;
         my $respairsig ;
         my $resnames ;
         while (my $line = <$contacts_fh>) {
   
            chomp $line ;
            if ($line =~ /^\#/) {next;}
            my @f = split(/\t/, $line) ;

            my $chain_id_1 = $f[$kdfield2no->{chain_id1}] ;
            my $resno_1 = $f[$kdfield2no->{resno1}] ;
            my $resna_1 = $f[$kdfield2no->{resna1}] ;
            my $atomna_1 = $f[$kdfield2no->{atomna1}] ;
            my $chain_id_2 = $f[$kdfield2no->{chain_id2}] ;
            my $resno_2 = $f[$kdfield2no->{resno2}] ;
            my $resna_2 = $f[$kdfield2no->{resna2}] ;
            my $atomna_2 = $f[$kdfield2no->{atomna2}] ;
            my $dist = $f[$kdfield2no->{dist}] ;

            if ((!defined $chain_id_1) || ($chain_id_1 eq '')) {
               $chain_id_1 = ' '; }

            if ((!defined $chain_id_2) || ($chain_id_2 eq '')) {
               $chain_id_2 = ' '; }
   
   
            my $sig1 = $resno_1."\n".$chain_id_1 ;
            my $sig2 = $resno_2."\n".$chain_id_2 ;
            my $ressig = $sig1."\n".$sig2 ;
   
            if ( (substr($atomna_1, 1, 1) ne 'H') &&
                 (substr($atomna_2, 1, 1) ne 'H') &&
                 (substr($atomna_1, 1, 1) ne 'Q') &&
                 (substr($atomna_2, 1, 1) ne 'Q') &&
                 !( $resno_1 eq $resno_2 &&
                 $chain_id_1 eq $chain_id_2) ) {
   
               if (!exists $resnames->{$sig1}) {
                  $resnames->{$sig1} = $resna_1 ; }
               if (!exists $resnames->{$sig2}) {
                  $resnames->{$sig2} = $resna_2 ; }
                  
               if (!exists $respairsig->{$ressig}) {
                  $respairsig->{$ressig}->{contacts} = 0 ;
                  $respairsig->{$ressig}->{counts_4} = 0 ;
                  $respairsig->{$ressig}->{counts_4p5} = 0 ;
                  $respairsig->{$ressig}->{counts_5} = 0 ;
                  $respairsig->{$ressig}->{counts_5p5} = 0 ;
               }
   
               $respairsig->{$ressig}->{contacts}++ ;
               
               if ((!exists $respairsig->{$ressig}->{min_dist}) ||
                   ($dist <= $respairsig->{$ressig}->{min_dist})) {
                  $respairsig->{$ressig}->{min_dist} = $dist ;
               }
   
               if ($dist <= 4) {
                  $respairsig->{$ressig}->{counts_4}++ ;
                  $respairsig->{$ressig}->{counts_4p5}++ ;
                  $respairsig->{$ressig}->{counts_5}++ ;
                  $respairsig->{$ressig}->{counts_5p5}++ ;
               } elsif ($dist <= 4.5) {
                  $respairsig->{$ressig}->{counts_4p5}++ ;
                  $respairsig->{$ressig}->{counts_5}++ ;
                  $respairsig->{$ressig}->{counts_5p5}++ ;
               } elsif ($dist <= 5) {
                  $respairsig->{$ressig}->{counts_5}++ ;
                  $respairsig->{$ressig}->{counts_5p5}++ ;
               } elsif ($dist <= 5.5) {
                  $respairsig->{$ressig}->{counts_5p5}++ ;
               }
   
               $respair_dup->{$sig1}->{$sig2}++ ;
#            $respair_dup->{$sig2}->{$sig1}->{contacts}++ ;
# dont need to dup, as we are iterating over all residues in the domain,
#  we are bound to run across all contacts
            }
         }
         close($contacts_fh) ;
         unlink $kdcont_out->{contacts_fn} ;
   
         my $bscontacts ;
         my $bssse_con ;
         my $bssse_con_halves ;
         foreach my $sid12sig (keys %{$bs_res}) {
            my @sids = split(/\n/, $sid12sig) ;
            foreach my $side (0 .. 1) {
               my $sid = $sids[$side] ;
               my $seen ;
               foreach my $res (keys %{$bs_res->{$sid12sig}->[$side]}) {
                  foreach my $res2 (keys %{$respair_dup->{$res}}) {
                     if ( exists $res2sub->{$res2}->{$sid} &&
                          !exists $seen->{$res}->{$res2}) {

                        my $res12sig = $res."\n".$res2 ;
                        $bscontacts->{$sid12sig}->{$sid}->{$res12sig}++ ;
                        $seen->{$res}->{$res2}++ ; $seen->{$res2}->{$res}++ ;
   # dups;
                        if ( exists $res2sse->{$res} &&
                             exists $res2sse->{$res2} ) {

                           my @t = ($res2sse->{$res}, $res2sse->{$res2}) ;
                           my @tsort = sort {$a <=> $b} @t ;
                           my $sse12 = join("\n", @tsort);
                           $bssse_con->{$sid12sig}->{$sid}->{$sse12}++ ;
                           $bssse_con_halves->{$sid12sig}->{$sid}->{$sse12}->{$t[0]}->{$res}++;
                           $bssse_con_halves->{$sid12sig}->{$sid}->{$sse12}->{$t[1]}->{$res2}++;
                        }
                     }
                  }
               }
            }
         }
   
         $out_fn->{bs_contacts} =
            "bindingsite_contacts_$bdp_id.out" ;
         open($out_fh->{bs_contacts}, ">".$out_fn->{bs_contacts}) ;

         $out_fn->{bs_sse_basic_contacts} =
            "bindingsite_secstrx_basic_contacts_$bdp_id.out" ;
         open($out_fh->{bs_sse_basic_contacts},
              ">".$out_fn->{bs_sse_basic_contacts}) ;
   
         foreach my $sid12sig (keys %{$bscontacts}) {
            my ($sid1, $sid2) = split(/\n/, $sid12sig) ;
            foreach my $sid (keys %{$bscontacts->{$sid12sig}}) {
               foreach my $res12sig (keys %{$bscontacts->{$sid12sig}->{$sid}}) {
                  my ($o_resno1, $o_chain1, $o_resno2, $o_chain2)=
                     split(/\n/, $res12sig) ;
                  my $res1 = $o_resno1."\n".$o_chain1 ;
                  my $res2 = $o_resno2."\n".$o_chain2 ;
                  my $o_resna1 = $resnames->{$res1} ;
                  my $o_resna2 = $resnames->{$res2} ;
                  
                  my @outvals = ($bdp_id,
                                 $sid1,
                                 $sid2,
                                 $sid,
                                 $o_chain1,
                                 $o_resno1,
                                 $o_resna1,
                                 $o_chain2,
                                 $o_resno2,
                                 $o_resna2,
                                 $respairsig->{$res12sig}->{min_dist},
                                 $respairsig->{$res12sig}->{contacts},
                                 $respairsig->{$res12sig}->{counts_4},
                                 $respairsig->{$res12sig}->{counts_4p5},
                                 $respairsig->{$res12sig}->{counts_5},
                                 $respairsig->{$res12sig}->{counts_5p5}) ;
   
                  print {$out_fh->{bs_contacts}} join("\t", @outvals)."\n";
               }
               
               foreach my $sse12 (keys %{$bssse_con->{$sid12sig}->{$sid}}) {
                  my ($sse1, $sse2) = split(/\n/, $sse12) ;
                  my $numres1 = keys %{$bssse_con_halves->{$sid12sig}->{$sid}->{$sse12}->{$sse1}};
                  my $numres2 = keys %{$bssse_con_halves->{$sid12sig}->{$sid}->{$sse12}->{$sse2}};
                  my @outvals = ( $bdp_id, $sid1, $sid2, $sid,
                                  $sse1, $sse2,
                                  $numres1, $numres2,
                                  $bssse_con->{$sid12sig}->{$sid}->{$sse12},
                                  $sse_defs->{$sse1},
                                  $sse_defs->{$sse2}) ;
                  print {$out_fh->{bs_sse_basic_contacts}} 
                     join("\t", @outvals)."\n";
               }
            }
         }
         close ($out_fh->{bs_sse_basic_contacts}) ;
         close ($out_fh->{bs_contacts}) ;

         if (-s $out_fn->{bs_contacts}) {
            system("gzip ".$out_fn->{bs_contacts}) ;
            $out_fn->{bs_contacts} .= '.gz' ;
            push @movethese, $out_fn->{bs_contacts} ;

            my $bsc_depositdir =
               $pibase_specs->{metatod_dir}->{bindingsite_contacts}.
              '/'.POSIX::floor($bdp_id / 1000) ;
            push @moveto, $bsc_depositdir ;

            my $t_filena = $bsc_depositdir.'/'.
                           $out_fn->{bs_contacts} ;
            my @outvals = ("BSC", $bdp_id, $param_dist_cutoff,'\N',
                           $t_filena);
            print join("\t", @outvals)."\n" ;
         } else {
            print STDERR "BSC\t$out_fn->{bs_contacts} is empty\n" ;
            unlink $out_fn->{bs_contacts} ;
         }
   
         if (-s $out_fn->{bs_sse_basic_contacts}) {
            system("gzip ".$out_fn->{bs_sse_basic_contacts}) ;
            $out_fn->{bs_sse_basic_contacts} .= '.gz' ;
            push @movethese, $out_fn->{bs_sse_basic_contacts} ;

            my $bsbc_depositdir =
   $pibase_specs->{metatod_dir}->{bindingsite_secstrx_basic_contacts}.
              '/'.POSIX::floor($bdp_id / 1000) ;
            push @moveto, $bsbc_depositdir ;

            my $t_filena = $bsbc_depositdir.'/'.
                           $out_fn->{bs_sse_basic_contacts} ;

            my @outvals = ("BSBC", $bdp_id, '\N', $t_filena) ;
            print join("\t", @outvals)."\n" ;
         } else {
            print STDERR "BSBC\t$out_fn->{bs_sse_basic_contacts} is empty\n" ;
            unlink $out_fn->{bs_sse_basic_contacts} ;
         }
   
   
         if ($#movethese >= $movethreshold) {
            foreach my $j (0 .. $#movethese) {
               my $t_fn = $movethese[$j] ;
               my $t_dir = $moveto[$j] ;
               pibase::safe_move($t_fn, $t_dir);}
            @movethese = () ; @moveto = () ;
         }
      }

      foreach my $j (0 .. $#movethese) {
         my $t_fn = $movethese[$j] ;
         my $t_dir = $moveto[$j] ;
         pibase::safe_move($t_fn, $t_dir);
      }
   }


   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_residues_tables file (specified in specs) into pibase
      $import_status->{interface_secstrx_basic_contacts_tables} =
         pibase::mysqlimport({
         pibase_specs => $pibase_specs,
 fn => $pibase_specs->{buildfiles}->{interface_secstrx_basic_contacts_tables}
      }) ;

      $import_status->{bindingsite_secstrx_basic_contacts_tables} =
         pibase::mysqlimport({
         pibase_specs => $pibase_specs,
 fn => $pibase_specs->{buildfiles}->{bindingsite_secstrx_basic_contacts_tables}
      }) ;

      $import_status->{bindingsite_contacts_tables} =
         pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{bindingsite_contacts_tables}}) ;
   }

   return $import_status ;

}


=head2 calc_interface_secstrx_contacts()

   Title:       calc_interface_secstrx_contacts()
   Function:    Computes contacts between secondary structure elements
                  across PIBASE interfaces

=cut

sub calc_interface_secstrx_contacts {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }
   my $modeller_bin = $pibase_specs->{binaries}->{modeller} ;
   my $move_thresh = 100 ;

   my $bdp2secstrx_fn;
   {
      my ($t_bdp, $t_secstrx_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM bdp_secstrx_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2secstrx_fn->{$t_bdp->[$j]} = $t_secstrx_fn->[$j]  ; }
   }

   my $bdp2intcon_fn;
   {
      my ($t_bdp, $t_intcon_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM interface_contacts_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2intcon_fn->{$t_bdp->[$j]} = $t_intcon_fn->[$j]  ; }
   }

   my $bdp_list ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line =~ /^\#/) {next;}
         chomp $line;
         my $bdp_id = $line ;
         $bdp_list->{$line}++  ;
      }
      close(INF) ;
   } else {
      my ($t_bdp) = pibase::rawselect_tod(
         "SELECT bdp_id FROM intersubset_contacts") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp_list->{$t_bdp->[$j]}++ ;}
   }

   
   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      print "* calc_interface_secstrx_contacts() ".localtime()
         if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{calc_interface_secstrx_contacts_in},
       $temp_fn->{calc_interface_secstrx_contacts_in}) =
         tempfile("splits_calc_interface_secstrx_contacts_input.XXXXX");

      ($temp_fh->{calc_interface_secstrx_contacts_out},
       $temp_fn->{calc_interface_secstrx_contacts_out}) =
         tempfile("splits_calc_interface_secstrx_contacts_SGEout_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_secstrx_contacts_out}) ;

      ($temp_fh->{calc_interface_secstrx_contacts_err},
       $temp_fn->{calc_interface_secstrx_contacts_err}) =
         tempfile("splits_calc_interface_secstrx_contacts_SGEerr_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_secstrx_contacts_err}) ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         print {$temp_fh->{calc_interface_secstrx_contacts_in}} $bdp_id."\n" ; }
      close($temp_fh->{calc_interface_secstrx_contacts_in}) ;

      my $split_dir = tempdir("splits_calc_interface_secstrx_contacts.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{calc_interface_secstrx_contacts_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.calc_interface_secstrx_contacts.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/calc_interface_secstrx_contacts/ ;

main() ;

sub main {

         pibase::data::calc::calc_interface_secstrx_contacts({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.calc_interface_secstrx_contacts.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.calc_interface_secstrx_contacts.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id}) ;
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{calc_interface_secstrx_contacts_out},
         err_fn => $temp_fn->{calc_interface_secstrx_contacts_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{calc_interface_secstrx_contacts_out},
           $temp_fn->{calc_interface_secstrx_contacts_out}) ;

      open(REALOUTF,
 ">".$pibase_specs->{buildfiles}->{interface_secstrx_contacts_tables});

      while (my $line = readline($temp_fh->{calc_interface_secstrx_contacts_out})){
         if ($line =~ /^\#/) { next; }
         print REALOUTF $line ;
      }
      close($temp_fh->{calc_interface_secstrx_contacts_out}) ;
      close(REALOUTF) ;

   } else {

      my (@movethese, @moveto) ;
      my $movethreshold = 100 ;

      my $tempdir = tempdir(CLEANUP => 1) ;
      chdir $tempdir ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
   
         print STDERR "now on bdp $bdp_id\n" ;

         if (! exists $bdp2intcon_fn->{$bdp_id} ) {
            print STDERR "ERROR: interface_contacts_tables not found ".
                         "for bdp $bdp_id\n" ;
            next;
         }
   
         if (!exists $bdp2secstrx_fn->{$bdp_id}) {
            print STDERR "ERROR: bdp_secstrx_tables not found ".
                         "for bdp $bdp_id\n" ;
            next ;
         }
 

         my $res2sse;
         my $sse_defs;
         {
            my ($t_resno, $t_chain_id, $t_sseid, $t_sse) =
               pibase::rawselect_metatod( $bdp2secstrx_fn->{$bdp_id},
               "SELECT resno, chain_id, sse_id, sse ".
               "FROM $bdp2secstrx_fn->{$bdp_id}") ;

            foreach my $j ( 0 .. $#{$t_resno}) {
               my $ressig = $t_resno->[$j]."\n".$t_chain_id->[$j] ;
               $res2sse->{$ressig} = $t_sseid->[$j] ;
               $sse_defs->{$t_sseid->[$j]} = $t_sse->[$j] ;
            }
         }
   

         my ($out_fn, $out_fh) ;
         $out_fn = "interface_secstrx_contacts_$bdp_id.out",
         open($out_fh, ">$out_fn") ;
   
# Get all interface residue contacts

         my ( $subset_id_1, $subset_id_2,
              $chain_id_1, $resno_1, $chain_id_2, $resno_2 ) =
            pibase::rawselect_metatod($bdp2intcon_fn->{$bdp_id},
      "SELECT subset_id_1, subset_id_2, chain_id_1, resno_1, chain_id_2, resno_2 FROM $bdp2intcon_fn->{$bdp_id}") ;
   
   
         my ($sse_contacts, $sse_halves) ;
   
#   Iterate over intersubset contacts
   
         my $intersubset_contacts ;
         foreach my $j ( 0 .. $#{$subset_id_1} ) {
   
#   Translate residue contacts into sse contacts
   
            if ((!defined $chain_id_1->[$j]) || ($chain_id_1->[$j] eq '')) {
               $chain_id_1->[$j] = ' ' ; }
   
            if ((!defined $chain_id_2->[$j]) || ($chain_id_2->[$j] eq '')) {
               $chain_id_2->[$j] = ' ' ; }
   
            my $ressig1 = $resno_1->[$j]."\n".$chain_id_1->[$j] ;
            my $ressig2 = $resno_2->[$j]."\n".$chain_id_2->[$j] ;
            if ((! exists $res2sse->{$ressig1}) ||
                ($res2sse->{$ressig1} eq ''))  {

               my $t = "   warning: no SSE assignment for $ressig1" ;
               $t =~ s/\n/:/g ;
               print STDERR $t."\n" ;
               next;
   
            } elsif ((! exists $res2sse->{$ressig2}) ||
                      ($res2sse->{$ressig2} eq ''))  {

               my $t = "   warning: no SSE assignment for $ressig2" ;
               $t =~ s/\n/:/g ;
               print STDERR $t."\n" ;
               next;
            }
            my $ressig12 = $ressig1."\n".$ressig2 ;
            my $sse1 = $res2sse->{$ressig1} ;
            my $sse2 = $res2sse->{$ressig2} ;
            my $sse12 = $sse1."\n".$sse2 ;
   
            my $sub1 = $subset_id_1->[$j] ;
            my $sub2 = $subset_id_2->[$j] ;
            my $sub12 = $sub1."\n".$sub2 ;
   
            if (! exists $sse_contacts->{$sub12} ) {
               push @{$intersubset_contacts}, $sub12 ; }
   
            $sse_contacts->{$sub12}->{$sse12}->{$ressig12}++ ;
            $sse_halves->{$sub12}->{$sse12}->{sub1}->{$ressig1}++ ;
            $sse_halves->{$sub12}->{$sse12}->{sub2}->{$ressig2}++ ;
   
         }
   
   
         foreach my $sub12 (@{$intersubset_contacts}) {
            my ($sub1, $sub2) = split("\n", $sub12) ;
            foreach my $sse12 (keys %{$sse_contacts->{$sub12}}) {
               my ($sse1, $sse2) = split("\n", $sse12) ;
               my @outvals ;
               push @outvals, ($bdp_id, $sub1, $sub2, $sse1, $sse2) ;
   
               my $a = keys %{$sse_halves->{$sub12}->{$sse12}->{sub1}} ;
               my $b = keys %{$sse_halves->{$sub12}->{$sse12}->{sub2}} ;
               my $c = keys %{$sse_contacts->{$sub12}->{$sse12}} ;
   
               push @outvals, $a, $b, $c ;
               push @outvals, $sse_defs->{$sse1}, $sse_defs->{$sse2} ;
               print $out_fh join("\t", @outvals)."\n" ;
            }
         }
         close($out_fh) ;
   
         if (-s $out_fn) {
            system("gzip ".$out_fn) ;
            $out_fn .= '.gz' ;
            push @movethese, $out_fn ;

            my $isc_depositdir = 
               $pibase_specs->{metatod_dir}->{interface_secstrx_contacts}.
              '/'.POSIX::floor($bdp_id / 1000) ;
            push @moveto, $isc_depositdir ;

            my $t_filena = $isc_depositdir."/".$out_fn ;
            print join("\t", $bdp_id, '\N', $t_filena)."\n";
         } else {
            print STDERR "ERROR: (bdp $bdp_id) $out_fn is empty\n" ;
            unlink $out_fn ;
         }

         if ($#movethese >= $movethreshold) {
            foreach my $j (0 .. $#movethese) {
               my $t_fn = $movethese[$j] ;
               my $t_dir = $moveto[$j] ;
               pibase::safe_move($t_fn, $t_dir);}
            @movethese = () ; @moveto = () ;
         }
   
      }

      foreach my $j (0 .. $#movethese) {
         my $t_fn = $movethese[$j] ;
         my $t_dir = $moveto[$j] ;
         pibase::safe_move($t_fn, $t_dir);
      }

   }


   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_residues_tables file (specified in specs) into pibase
      $import_status->{interface_secstrx_contacts_tables} =
         pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{interface_secstrx_contacts_tables}
      }) ;
   }

   return $import_status ;

}


=head2 calc_interface_secstrx_profile()

   Title:       calc_interface_secstrx_profile()
   Function:    Generates a summary profile to describe secondary structure
                  content and contacts across PIBASE interfaces

=cut

sub calc_interface_secstrx_profile {

   require Bit::Vector ;

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }
   my $modeller_bin = $pibase_specs->{binaries}->{modeller} ;
   my $move_thresh = 100 ;

   my $intsecstrx_fn;
   {
      my ($t_bdp, $t_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM interface_secstrx_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $intsecstrx_fn->{$t_bdp->[$j]} = $t_fn->[$j]  ; }
   }

   my $intsecstrx_cont_fn;
   {
      my ($t_bdp, $t_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM interface_secstrx_contacts_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $intsecstrx_cont_fn->{$t_bdp->[$j]} = $t_fn->[$j]  ; }
   }

   my $bdp_list ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line =~ /^\#/) {next;}
         chomp $line;
         my $bdp_id = $line ;
         $bdp_list->{$line}++  ;
      }
      close(INF) ;
   } else {
      my ($t_bdp) = pibase::rawselect_tod(
         "SELECT bdp_id FROM intersubset_contacts") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp_list->{$t_bdp->[$j]}++ ;}
   }

   
   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      print "* calc_interface_secstrx_profile() ".localtime()
         if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{calc_interface_secstrx_profile_in},
       $temp_fn->{calc_interface_secstrx_profile_in}) =
         tempfile("splits_calc_interface_secstrx_profile_input.XXXXX");

      ($temp_fh->{calc_interface_secstrx_profile_out},
       $temp_fn->{calc_interface_secstrx_profile_out}) =
         tempfile("splits_calc_interface_secstrx_profile_SGEout_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_secstrx_profile_out}) ;

      ($temp_fh->{calc_interface_secstrx_profile_err},
       $temp_fn->{calc_interface_secstrx_profile_err}) =
         tempfile("splits_calc_interface_secstrx_profile_SGEerr_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_secstrx_profile_err}) ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         print {$temp_fh->{calc_interface_secstrx_profile_in}} $bdp_id."\n" ; }
      close($temp_fh->{calc_interface_secstrx_profile_in}) ;

      my $split_dir = tempdir("splits_calc_interface_secstrx_profile.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{calc_interface_secstrx_profile_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.calc_interface_secstrx_profile.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/calc_interface_secstrx_profile/ ;

main() ;

sub main {

         pibase::data::calc::calc_interface_secstrx_profile({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.calc_interface_secstrx_profile.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.calc_interface_secstrx_profile.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id}) ;
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{calc_interface_secstrx_profile_out},
         err_fn => $temp_fn->{calc_interface_secstrx_profile_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{calc_interface_secstrx_profile_out},
           $temp_fn->{calc_interface_secstrx_profile_out}) ;

      open(REALOUTF,
         ">".$pibase_specs->{buildfiles}->{interface_secstrx_profile});
      while (my $line=readline($temp_fh->{calc_interface_secstrx_profile_out})){
         if ($line =~ /^\#/) { next; }
         print REALOUTF $line ;
      }
      close($temp_fh->{calc_interface_secstrx_profile_out}) ;
      close(REALOUTF) ;

   } else {

      my $sse2num = {
         'H' => 0 ,
         'B' => 1 ,
         'E' => 2 ,
         'G' => 3 ,
         'I' => 4 ,
         'T' => 5 ,
         'S' => 6 ,
         ' ' => 7
      } ;

      my $ssesse2num ;
      my $ind = 0 ;
      foreach my $j ( 0 .. 7) {
         foreach my $k ( $j .. 7) {
            $ssesse2num->{$j."\n".$k} = $ind ;
            $ssesse2num->{$k."\n".$j} = $ind ;
            $ind++ ;
         }
      }

      my $bounds ;
      $bounds->{'sse2num'} = 7 ;
      $bounds->{'ssesse2num'} = ($ind - 1) ; 


      my $tempdir = tempdir(CLEANUP => 1) ;
      chdir $tempdir ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         print STDERR "now on bdp_id $bdp_id\n" ;
   
         if (! exists $intsecstrx_cont_fn->{$bdp_id}) {
            print STDERR "ERROR: (bdp $bdp_id) interface_secstrx_contacts".
                         " not found\n";
            next ;
         }
   
         if (! exists $intsecstrx_fn->{$bdp_id}) {
            print STDERR "ERROR: (bdp $bdp_id) interface_secstrx not found\n";
            next ;
         }
   
         my ($fn, $fh, $compress_fl) ;
         my ($subset_id_1, $subset_id_2, $sse_id_1, $sse_id_2, $sse_1,
             $sse_2, $numres_1, $numres_2, $numres_12) =
             pibase::rawselect_metatod( $intsecstrx_cont_fn->{$bdp_id},
            "SELECT subset_id_1, subset_id_2, sse_id_1, sse_id_2, ".
            "sse_1, sse_2, num_res_1, num_res_2, num_res_12 ".
            "FROM $intsecstrx_cont_fn->{$bdp_id}") ;
   
#HERENOW, change the variable names of the retrieval to reflect the contenv s contacts nature
   
         my ($s_subset_id_1, $s_subset_id_2, $s_subset_id, $s_sse_id,
             $s_numres) = pibase::rawselect_metatod( $intsecstrx_fn->{$bdp_id},
            "SELECT subset_id_1, subset_id_2, subset_id, sse_id, num_res ".
            "FROM $intsecstrx_fn->{$bdp_id}") ;
   
         my $sse_content ;
         foreach my $j ( 0 .. $#{$s_subset_id_1}) {
            my $sid_pair  = $s_subset_id_1->[$j]."\n".$s_subset_id_2->[$j] ;
            $sse_content->{$sid_pair}->{$s_subset_id->[$j]}->{$s_sse_id->[$j]} =
               $s_numres->[$j] ;
         }
   
         my $sse_ids ;
         my ($sse_types, $sse_single_types_1, $sse_single_types_2) ;
         my $sse_defs ;
   
         my $sse_id_pairs ;
         my $sseres ;
   
         foreach my $j ( 0 .. $#{$subset_id_1}) {
         
            if ((!defined $sse_1->[$j]) || ($sse_1->[$j] eq '')) {
               $sse_1->[$j] = ' '; }
   
            if ((!defined $sse_2->[$j]) || ($sse_2->[$j] eq '')) {
               $sse_2->[$j] = ' '; }
   
            my $sid_pair  = $subset_id_1->[$j]."\n".$subset_id_2->[$j] ;
            my $sse_id_pair = $sse_id_1->[$j]."\n".$sse_id_2->[$j] ;
   
#   Record the sse seen on each side as well as the contact
   
            $sse_defs->{$sid_pair}->{$sse_id_1->[$j]}= $sse2num->{$sse_1->[$j]};
            $sse_defs->{$sid_pair}->{$sse_id_2->[$j]}= $sse2num->{$sse_2->[$j]};
   
#	 print STDERR "defined for $subset_id_1->[$j] -- $subset_id_2->[$j] : $sse_id_1->[$j] -- $sse_id_2->[$j]\n" ;
   
            $sse_id_pairs->{$sid_pair}->{sub1}->{$sse_id_1->[$j]}++ ;
            $sse_id_pairs->{$sid_pair}->{sub2}->{$sse_id_2->[$j]}++ ;
            $sse_id_pairs->{$sid_pair}->{sub12}->{$sse_id_pair}++ ;
   
            $sseres->{$sid_pair}->{sub1}->{$sse_id_1->[$j]} = $numres_1->[$j];
            $sseres->{$sid_pair}->{sub2}->{$sse_id_2->[$j]} = $numres_2->[$j];
            $sseres->{$sid_pair}->{sub12}->{$sse_id_pair} = $numres_12->[$j];
   
         }
   
   
# Count up:
# num of SSE elements on each side of the interface
# num of SSE pairs across the interface
# num of SSE types on each side of the interface
# num of SSE pair types across the interface
# num of residues of each SSE type on each side of the interface
# num of SSE pair types across the interface
   
         my $sseres_freq ;
         my $ssetype_freq ;
         my $ssetype_norms ;
         my $sse_freq ;
   
         foreach my $sid_pair ( keys (%{$sse_id_pairs}) ) {
   
            my ($sub1, $sub2) = split(/\n/, $sid_pair) ;
   
            foreach my $sse_id_1 (keys (%{$sse_id_pairs->{$sid_pair}->{sub1}})){
   
#               print STDERR "sse_id_1 : $sse_id_1\n" ;
#               print STDERR "  defed to ".$sse_defs->{$sid_pair}->{$sse_id_1}."\n" ;
 $ssetype_freq->{$sid_pair}->{sub1}->{$sse_defs->{$sid_pair}->{$sse_id_1}}++ ;
 $sseres_freq->{$sid_pair}->{sub1}->{$sse_defs->{$sid_pair}->{$sse_id_1}} +=
    $sse_content->{$sid_pair}->{$sub1}->{$sse_id_1} ;
            }
   
            foreach my $sse_id_2 (keys (%{$sse_id_pairs->{$sid_pair}->{sub2}})){
 $ssetype_freq->{$sid_pair}->{sub2}->{$sse_defs->{$sid_pair}->{$sse_id_2}}++ ; 
 $sseres_freq->{$sid_pair}->{sub2}->{$sse_defs->{$sid_pair}->{$sse_id_2}} +=
    $sse_content->{$sid_pair}->{$sub2}->{$sse_id_2} ;
            }
   
            foreach my $sse_id_12 (keys
               (%{$sse_id_pairs->{$sid_pair}->{sub12}})) {

               my ($sse_id_1, $sse_id_2)  = split (/\n/, $sse_id_12) ;
               my $ssetype_pair = $ssesse2num->{$sse_defs->{$sid_pair}->{$sse_id_1}."\n".$sse_defs->{$sid_pair}->{$sse_id_1}} ;
               $ssetype_freq->{$sid_pair}->{sub12}->{$ssetype_pair}++ ;
   
               $sseres_freq->{$sid_pair}->{sub12}->{$ssetype_pair} +=
                  $sseres->{$sid_pair}->{sub12}->{$sse_id_12} ;
            }
         }
   
   
         foreach my $sid_pair (keys %{$sse_id_pairs}) {
            my ($subset_id_1, $subset_id_2) = split(/\n/, $sid_pair) ;
   
            my $nums ;
            $nums->{num_sse1} = keys %{$sse_id_pairs->{$sid_pair}->{sub1}} ;
            $nums->{num_sse2} = keys %{$sse_id_pairs->{$sid_pair}->{sub2}} ;
            $nums->{num_sse12} = keys %{$sse_id_pairs->{$sid_pair}->{sub12}} ;

            my $ssetype_norms ;
            $ssetype_norms->{$sid_pair}->{sub1} =
               keys %{$ssetype_freq->{$sid_pair}->{sub1}} ;
   
            $ssetype_norms->{$sid_pair}->{sub2} =
               keys %{$ssetype_freq->{$sid_pair}->{sub2}} ;
   
            $ssetype_norms->{$sid_pair}->{sub12} =
               keys %{$ssetype_freq->{$sid_pair}->{sub12}} ;
   
            my $contact_vec = Bit::Vector->new(36) ; $contact_vec->Empty() ;
            my $sse_vec_1 = Bit::Vector->new(8) ; $sse_vec_1->Empty() ;
            my $sse_vec_2 = Bit::Vector->new(8) ; $sse_vec_2->Empty() ;
   
            foreach my $sse_id_pair (
               keys %{$sse_id_pairs->{$sid_pair}->{sub12}}) {
   
               my ($sse_id_1, $sse_id_2) = split(/\n/, $sse_id_pair) ;
   
               $sse_vec_1->Bit_On($sse_defs->{$sid_pair}->{$sse_id_1}) ;
               $sse_vec_2->Bit_On($sse_defs->{$sid_pair}->{$sse_id_2}) ;
               $contact_vec->Bit_On($ssesse2num->{$sse_defs->{$sid_pair}->{$sse_id_1}."\n".$sse_defs->{$sid_pair}->{$sse_id_2}}) ;
   
            }
   
            my @outvals = ( $bdp_id, $subset_id_1, $subset_id_2,
               $nums->{num_sse1}, $nums->{num_sse2},
               $nums->{num_sse12}, $contact_vec->to_Bin(),
               $sse_vec_1->to_Bin(), $sse_vec_2->to_Bin(),
               $contact_vec->Norm(), $sse_vec_1->Norm(), $sse_vec_2->Norm() ) ;
   
   
            my $o_sseres ;
            my $o_ssetype ;
            foreach my $sset (0 .. $bounds->{'sse2num'}) {
               if (exists $ssetype_freq->{$sid_pair}->{sub1}->{$sset}) {
                  push @{$o_ssetype->[0]},
                     $ssetype_freq->{$sid_pair}->{sub1}->{$sset} ;
                  push @{$o_sseres->[0]},
                     $sseres_freq->{$sid_pair}->{sub1}->{$sset} ;
               } else {
                  push @{$o_ssetype->[0]}, 0 ;
                  push @{$o_sseres->[0]}, 0 ;
               }
               
               if (exists $ssetype_freq->{$sid_pair}->{sub2}->{$sset}) {
                  push @{$o_ssetype->[1]},
                     $ssetype_freq->{$sid_pair}->{sub2}->{$sset} ;
                  push @{$o_sseres->[1]},
                     $sseres_freq->{$sid_pair}->{sub2}->{$sset} ;
               } else {
                  push @{$o_ssetype->[1]}, 0 ;
                  push @{$o_sseres->[1]}, 0 ;
               }
            }
   
   
            foreach my $ssett (0 .. $bounds->{'ssesse2num'}) {
   
               if (exists $ssetype_freq->{$sid_pair}->{sub12}->{$ssett}) {
                  push @{$o_ssetype->[2]},
                     $ssetype_freq->{$sid_pair}->{sub12}->{$ssett} ;
                  push @{$o_sseres->[2]},
                     $sseres_freq->{$sid_pair}->{sub12}->{$ssett} ;
               } else {
                  push @{$o_ssetype->[2]}, 0 ;
                  push @{$o_sseres->[2]}, 0 ;
               }
            }
   
            push @outvals, @{$o_ssetype->[0]} ;
            push @outvals, @{$o_ssetype->[1]} ;
            push @outvals, @{$o_ssetype->[2]} ;
   
   
            push @outvals, @{$o_sseres->[0]} ;
            push @outvals, @{$o_sseres->[1]} ;
            push @outvals, @{$o_sseres->[2]} ;
   
            print join("\t", @outvals)."\n" ;
         }
      }
   }


   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_residues_tables file (specified in specs) into pibase
      $import_status->{interface_secstrx_profile} =
         pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{interface_secstrx_profile}
      }) ;
   }

   return $import_status ;

}


=head2 calc_interface_resvector()

   Title:       calc_interface_resvector()
   Function:    Calculates a bit vector to describe amino acid type pairs
                  across PIBASE interfaces

=cut

sub calc_interface_resvector {

   require Bit::Vector ;

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }

   my $bdp2intcon_fn;
   {
      my ($t_bdp, $t_intcon_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM interface_contacts_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2intcon_fn->{$t_bdp->[$j]} = $t_intcon_fn->[$j]  ; }
   }

   my $bdp_list ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line =~ /^\#/) {next;}
         chomp $line;
         my $bdp_id = $line ;
         $bdp_list->{$line}++  ;
      }
      close(INF) ;
   } else {
      my ($t_bdp) = pibase::rawselect_tod(
         "SELECT bdp_id FROM intersubset_contacts") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp_list->{$t_bdp->[$j]}++ ;}
   }

   
   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      print "* calc_interface_resvector() ".localtime()
         if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{calc_interface_resvector_in},
       $temp_fn->{calc_interface_resvector_in}) =
         tempfile("splits_calc_interface_resvector_input.XXXXX");

      ($temp_fh->{calc_interface_resvector_out},
       $temp_fn->{calc_interface_resvector_out}) =
         tempfile("splits_calc_interface_resvector_SGEout_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_resvector_out}) ;

      ($temp_fh->{calc_interface_resvector_err},
       $temp_fn->{calc_interface_resvector_err}) =
         tempfile("splits_calc_interface_resvector_SGEerr_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_resvector_err}) ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         print {$temp_fh->{calc_interface_resvector_in}} $bdp_id."\n" ; }
      close($temp_fh->{calc_interface_resvector_in}) ;

      my $split_dir = tempdir("splits_calc_interface_resvector.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{calc_interface_resvector_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.calc_interface_resvector.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/calc_interface_resvector/ ;

main() ;

sub main {

         pibase::data::calc::calc_interface_resvector({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.calc_interface_resvector.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.calc_interface_resvector.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id}) ;
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{calc_interface_resvector_out},
         err_fn => $temp_fn->{calc_interface_resvector_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{calc_interface_resvector_out},
           $temp_fn->{calc_interface_resvector_out}) ;

      open(REALOUTF,
         ">".$pibase_specs->{buildfiles}->{interface_resvector});
      while (my $line=readline($temp_fh->{calc_interface_resvector_out})){
         if ($line =~ /^\#/) { next; }
         print REALOUTF $line ;
      }
      close($temp_fh->{calc_interface_resvector_out}) ;
      close(REALOUTF) ;

   } else {

#Create pointers from amino acids to indices and amino acid pairs to indices

      my $resnum = {
      'ALA' => 0 ,
      'ARG' => 1 ,
      'ASN' => 2 ,
      'ASP' => 3 ,
      'CYS' => 4 ,
      'GLN' => 5 ,
      'GLU' => 6 ,
      'GLY' => 7 ,
      'HIS' => 8 ,
      'HSD' => 8 ,
      'ILE' => 9 ,
      'LEU' => 10 ,
      'LYS' => 11 ,
      'MET' => 12 ,
      'PHE' => 13 ,
      'PRO' => 14 ,
      'SER' => 15 ,
      'THR' => 16 ,
      'TRP' => 17 ,
      'TYR' => 18 ,
      'VAL' => 19
      } ;

      my $resresnum ;
      my $ind = 0 ;
      foreach my $j ( 0 .. 19) {
         foreach my $k ( $j .. 19) {
            $resresnum->{$j."\n".$k} = $ind ;
            $resresnum->{$k."\n".$j} = $ind ;
            $ind++ ;
         }
      }

      my $tempdir = tempdir(CLEANUP => 1) ;
      chdir $tempdir ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {

         print STDERR "NOW on $bdp_id\n" ;
   
         if (! exists $bdp2intcon_fn->{$bdp_id} ) {
            print STDERR "ERROR: (bdp $bdp_id) interface_contacts_tables entry".
                         "not found\n";
            next;
         }

#   Get all residue-residue names
   
         my ($subset_id_1, $subset_id_2, $resna_1, $resna_2) =
            pibase::rawselect_metatod( $bdp2intcon_fn->{$bdp_id},
               "SELECT subset_id_1, subset_id_2, resna_1, resna_2 ".
               "FROM $bdp2intcon_fn->{$bdp_id}" ) ;
   
         my $composition ;
         foreach my $j ( 0 .. $#{$subset_id_1}) {
            $composition->{$subset_id_1->[$j]."\n".$subset_id_2->[$j]}->{$resna_1->[$j]."\n".$resna_2->[$j]}++ ;
         }
   
         my $interface_rev ;
         foreach my $interface (keys %{$composition}) {
            my $nonprot = 0 ;
            my ($sub1, $sub2) = split(/\n/, $interface) ;
   
            my $contact_vec = Bit::Vector->new(210) ; $contact_vec->Empty() ;
            my $res_vec_1 = Bit::Vector->new(20) ; $res_vec_1->Empty() ;
            my $res_vec_2 = Bit::Vector->new(20) ; $res_vec_2->Empty() ; 

            foreach my $contact (keys %{$composition->{$interface}}) {
               my ($res1, $res2) = split(/\n/, $contact) ;

               if (!exists $resnum->{$res1} || !exists $resnum->{$res2}) {
                  print STDERR "NOTE ($bdp_id $sub1 $sub2): ".
                     "non-amino acid contacts noted ($res1 -- $res2)\n" ;
                  $nonprot++ ; last;
               }

               $res_vec_1->Bit_On($resnum->{$res1}) ;
               $res_vec_2->Bit_On($resnum->{$res2}) ;
               $contact_vec->Bit_On($resresnum->{$resnum->{$res1}."\n".$resnum->{$res2}}) ;
            }

            if ($nonprot) {last;}
   
            my @outvals = ($bdp_id, $sub1, $sub2, $contact_vec->to_Bin(),
                           $res_vec_1->to_Bin(), $res_vec_2->to_Bin(),
                           $contact_vec->Norm(),
                           $res_vec_1->Norm(), $res_vec_2->Norm()) ;
   
            print join("\t", @outvals)."\n" ;
         }
      }
   }


   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_residues_tables file (specified in specs) into pibase
      $import_status->{interface_resvector} =
         pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{interface_resvector}
      }) ;
   }

   return $import_status ;

}


=head2 calc_interface_sse_topology()

   Title:       calc_interface_sse_topology()
   Function:    Calculates topology of secondary structure element contacts
                  across PIBASE interfaces

=cut

sub calc_interface_sse_topology {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }

   my $sse2shape = {
      "H" => "fcircle",
      "B" => "ftriangle",
      "T" => "fbox",
      " " => "box",
      "" => "box"
   } ;

   my @sub_colors = qw/black red blue green brown grey orange AliceBlue RosyBrown beige coral sienna orchid plum bisque Peach khaki OliveDrab SaddleBrown LimeGreen turquoise SteelBlue CornflowerBlue DeepSkyBlue aquamarine chocolate SpringGreen IndianRed DarkViolet maroon/ ;

   my $tables;
   {
   {
      my ($t_bdp, $t_fn) = pibase::rawselect_tod("SELECT bdp_id, source_file ".
         " FROM bindingsite_secstrx_basic_contacts_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $tables->{bssse_con}->{$t_bdp->[$j]} = $t_fn->[$j]  ; }
   }
   {
      my ($t_bdp, $t_fn) = pibase::rawselect_tod("SELECT bdp_id, source_file ".
         " FROM interface_secstrx_basic_contacts_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $tables->{isse_con}->{$t_bdp->[$j]} = $t_fn->[$j]  ; }
   }
   }

   my $bdp2sids ;
   my $sid2class ;
   my $sid2source ;
   my $domcolors ;
   {
      my @tod_res ;
      @tod_res = pibase::rawselect_tod(
         "SELECT bdp_id, subset_id, subset_source_id, class FROM subsets") ;
      foreach my $j ( 0 .. $#{$tod_res[0]}) {
         if (!defined $tod_res[0]->[$j] ||
             $tod_res[0]->[$j] eq "") { next ; }
         $sid2source->{$tod_res[1]->[$j]} = $tod_res[2]->[$j] ;
         $sid2class->{$tod_res[1]->[$j]} = $tod_res[3]->[$j] ;
         push @{$bdp2sids->{$tod_res[0]->[$j]}->{$sid2source->{$tod_res[1]->[$j]}}}, $tod_res[1]->[$j] ;
      }
      foreach my $bdp (keys %{$bdp2sids}) {
         foreach my $source (keys %{$bdp2sids->{$bdp}}) {
            my $t_classes ;
            foreach my $j ( 0 .. $#{$bdp2sids->{$bdp}->{$source}}) {
               my $t_sid = $bdp2sids->{$bdp}->{$source}->[$j] ;
               $t_classes->{$sid2class->{$t_sid}}++ ;
            }
            my @sort_classes = sort keys %{$t_classes} ;
            foreach my $j ( 0 .. $#sort_classes ) {
               $domcolors->{$bdp}->{$sort_classes[$j]} = $j ;
            }
         }
      }
   }


   my $bdp_list ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line =~ /^\#/) {next;}
         chomp $line;
         my $bdp_id = $line ;
         $bdp_list->{$line}++  ;
      }
      close(INF) ;
   } else {
      my ($t_bdp) = pibase::rawselect_tod(
         "SELECT bdp_id FROM intersubset_contacts") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp_list->{$t_bdp->[$j]}++ ;}
   }

   
   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      print "* calc_interface_sse_topology() ".localtime()
         if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{calc_interface_sse_topology_in},
       $temp_fn->{calc_interface_sse_topology_in}) =
         tempfile("splits_calc_interface_sse_topology_input.XXXXX");

      ($temp_fh->{calc_interface_sse_topology_out},
       $temp_fn->{calc_interface_sse_topology_out}) =
         tempfile("splits_calc_interface_sse_topology_SGEout_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_sse_topology_out}) ;

      ($temp_fh->{calc_interface_sse_topology_err},
       $temp_fn->{calc_interface_sse_topology_err}) =
         tempfile("splits_calc_interface_sse_topology_SGEerr_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_sse_topology_err}) ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         print {$temp_fh->{calc_interface_sse_topology_in}} $bdp_id."\n" ; }
      close($temp_fh->{calc_interface_sse_topology_in}) ;

      my $split_dir = tempdir("splits_calc_interface_sse_topology.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{calc_interface_sse_topology_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.calc_interface_sse_topology.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/calc_interface_sse_topology/ ;

main() ;

sub main {

         pibase::data::calc::calc_interface_sse_topology({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.calc_interface_sse_topology.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.calc_interface_sse_topology.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id}) ;
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{calc_interface_sse_topology_out},
         err_fn => $temp_fn->{calc_interface_sse_topology_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{calc_interface_sse_topology_out},
           $temp_fn->{calc_interface_sse_topology_out}) ;

      open(REALOUTF_B,
         ">".$pibase_specs->{buildfiles}->{bindingsite_sse_topology});
      open(REALOUTF_I,
         ">".$pibase_specs->{buildfiles}->{interface_sse_topology});
      while (my $line=readline($temp_fh->{calc_interface_sse_topology_out})){
         if ($line =~ /^\#/) {
            next;
         } elsif ($line =~ /^BS/) {
            $line =~ s/^BS\t// ;
            print REALOUTF_B $line ;
         } elsif ($line =~ /^INT/) {
            $line =~ s/^INT\t// ;
            print REALOUTF_I $line ;
         }
      }
      close($temp_fh->{calc_interface_sse_topology_out}) ;
      close(REALOUTF_I) ;
      close(REALOUTF_B) ;

   } else {

#Create pointers from amino acids to indices and amino acid pairs to indices

      my $tempdir = tempdir(CLEANUP => 1) ;
      chdir $tempdir ;
      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         print STDERR "now on $bdp_id\n" ;
         if (!exists $tables->{bssse_con}->{$bdp_id}) {
            print STDERR "ERROR ($bdp_id): no bindingsite_secstrx_basic_contacts entry for BDP $bdp_id\n" ; next;}
   
         if (!exists $tables->{isse_con}->{$bdp_id}) {
            print STDERR "ERROR ($bdp_id): no interface_secstrx_basic_contacts entry for BDP $bdp_id\n" ; next;}
   
         my ( $bs_sid1, $bs_sid2, $bs_sid,
              $bs_sse1, $bs_sse2,
              $bs_numres1, $bs_numres2, $bs_numres12,
              $bs_ssedef1, $bs_ssedef2 ) =
            pibase::rawselect_metatod( $tables->{bssse_con}->{$bdp_id},
            "SELECT subset_id_1, subset_id_2, subset_id, sse_id_1, sse_id_2, ".
            "num_res_1, num_res_2, num_res_12, sse_1, sse_2 ".
            "FROM ".$tables->{bssse_con}->{$bdp_id} ) ;
   
         my $all_edges ;
         my $all_ncols ;
         my $all_nshapes ;
         my $all_estyles ;
   
         my ( $i_sid1, $i_sid2,
              $i_sse1, $i_sse2,
              $i_numres1, $i_numres2, $i_numres12,
              $i_ssedef1, $i_ssedef2 ) =
            pibase::rawselect_metatod( $tables->{isse_con}->{$bdp_id},
               "SELECT subset_id_1, subset_id_2, sse_id_1, sse_id_2, ".
               "num_res_1, num_res_2, num_res_12, sse_1, sse_2 ".
               "FROM $tables->{isse_con}->{$bdp_id}" ) ;
   
         my $ssedef ;
         my $i_edges ;
         my $i_edges_dup ;
         my $intedgeseen ;
         foreach my $j ( 0 .. $#{$i_sid1}) {
            my $cursidsource = $sid2source->{$i_sid1->[$j]} ;
   #BUGFIX041118 - dont include sse's spanning the interface
            if ( $i_sse1->[$j] ne $i_sse2->[$j] ) {
   #         print STDERR "on $bdp_id: ( of $i_sid1->[$j] -- $i_sid2->[$j]) \n" ;
               $i_edges->{$cursidsource}->{$i_sid1->[$j]."\n".$i_sid2->[$j]}->{$i_sse1->[$j]}->{$i_sse2->[$j]}++ ;
               $intedgeseen->{$cursidsource}->{$i_sid1->[$j]."\n".$i_sid2->[$j]}->{$i_sse1->[$j]}++ ;
               $intedgeseen->{$cursidsource}->{$i_sid1->[$j]."\n".$i_sid2->[$j]}->{$i_sse2->[$j]}++ ;
   
               $i_edges_dup->{$cursidsource}->{$i_sid1->[$j]."\n".$i_sid2->[$j]}->{$i_sse1->[$j]}->{$i_sse2->[$j]}++ ;
               $i_edges_dup->{$cursidsource}->{$i_sid1->[$j]."\n".$i_sid2->[$j]}->{$i_sse2->[$j]}->{$i_sse1->[$j]}++ ;
               $ssedef->{$i_sse1->[$j]} = $i_ssedef1->[$j] ;
               $ssedef->{$i_sse2->[$j]} = $i_ssedef2->[$j] ;
   
#BUG:  BEWARE OF SSE that span the interface  - i.e. belong both to subset_id_1 and to subset_id_2
   
#if ( $i_sse1->[$j] ne $i_sse2->[$j] )
               $all_edges->{$cursidsource}->{$i_sid1->[$j]."\n".$i_sid2->[$j]}->{$i_sse1->[$j]}->{$i_sse2->[$j]}++ ;
               $all_estyles->{$cursidsource}->{$i_sid1->[$j]."\n".$i_sid2->[$j]}->{$i_sse1->[$j]}->{$i_sse2->[$j]} = "dashed" ;
               $all_nshapes->{$cursidsource}->{$i_sid1->[$j]."\n".$i_sid2->[$j]}->{$i_sse1->[$j]} = $sse2shape->{$ssedef->{$i_sse1->[$j]}};
               $all_nshapes->{$cursidsource}->{$i_sid1->[$j]."\n".$i_sid2->[$j]}->{$i_sse2->[$j]} = $sse2shape->{$ssedef->{$i_sse2->[$j]}};
            }
         }
   
#BUGfix041118_1620 - delete SSE nodes with only self-BS edges; they must have been barely in the cutoff for an interface edge; but got edged out (hehe) cause H contacts where taken out in constructing the bindingsite_contacts
# all nodes must have at least one cross-interface (dashed) edge
# move bindingsite graphing to after interface; so edges b/w nodes sin an already processed interface edge can be kicked out as the self-edges are built
   
   
         my $bs_edges ;
         my $bs_edges_dup ;
         foreach my $j ( 0 .. $#{$bs_sid1}) {
            my $cursidsource = $sid2source->{$bs_sid1->[$j]} ;
   #BUGFIX041118 there should be no self-sse edges within a given binding site
            if ( ($bs_sse1->[$j] ne $bs_sse2->[$j]) &&
                 (exists $intedgeseen->{$cursidsource}->{$bs_sid1->[$j]."\n".$bs_sid2->[$j]}->{$bs_sse1->[$j]}) &&
                 (exists $intedgeseen->{$cursidsource}->{$bs_sid1->[$j]."\n".$bs_sid2->[$j]}->{$bs_sse2->[$j]}) ) {
   
   #         print STDERR "on $bdp_id: $bs_sid->[$j] ( of $bs_sid1->[$j] -- $bs_sid2->[$j]) \n" ;
               $bs_edges->{$cursidsource}->{$bs_sid1->[$j]."\n".$bs_sid2->[$j]}->{$bs_sid->[$j]}->{$bs_sse1->[$j]}->{$bs_sse2->[$j]}++ ;

               $bs_edges_dup->{$cursidsource}->{$bs_sid1->[$j]."\n".$bs_sid2->[$j]}->{$bs_sid->[$j]}->{$bs_sse1->[$j]}->{$bs_sse2->[$j]}++ ;
               $bs_edges_dup->{$cursidsource}->{$bs_sid1->[$j]."\n".$bs_sid2->[$j]}->{$bs_sid->[$j]}->{$bs_sse2->[$j]}->{$bs_sse1->[$j]}++ ;
   
               $ssedef->{$bs_sse1->[$j]} = $bs_ssedef1->[$j] ;
               $ssedef->{$bs_sse2->[$j]} = $bs_ssedef2->[$j] ;
               
               $all_edges->{$cursidsource}->{$bs_sid1->[$j]."\n".$bs_sid2->[$j]}->{$bs_sse1->[$j]}->{$bs_sse2->[$j]}++ ;
               $all_estyles->{$cursidsource}->{$bs_sid1->[$j]."\n".$bs_sid2->[$j]}->{$bs_sse1->[$j]}->{$bs_sse2->[$j]} = "solid" ;
   
# COLOR by side of interface
               my $tcolor = $sub_colors[$domcolors->{$bdp_id}->{$sid2class->{$bs_sid->[$j]}}] ;
               $all_ncols->{$cursidsource}->{$bs_sid1->[$j]."\n".$bs_sid2->[$j]}->{$bs_sse1->[$j]}= $tcolor ;
               $all_ncols->{$cursidsource}->{$bs_sid1->[$j]."\n".$bs_sid2->[$j]}->{$bs_sse2->[$j]}= $tcolor ;
   
# SHAPE by type of SSE
#	    print STDERR "SHAPE $bs_sse1->[$j] ($ssedef->{$bs_sse1->[$j]}) gets a $sse2shape->{$ssedef->{$bs_sse1->[$j]}}\n" ;
               $all_nshapes->{$cursidsource}->{$bs_sid1->[$j]."\n".$bs_sid2->[$j]}->{$bs_sse1->[$j]} = $sse2shape->{$ssedef->{$bs_sse1->[$j]}};
               $all_nshapes->{$cursidsource}->{$bs_sid1->[$j]."\n".$bs_sid2->[$j]}->{$bs_sse2->[$j]} = $sse2shape->{$ssedef->{$bs_sse2->[$j]}};

            }
         }


#      print STDERR "sending $bdp_id for printing\n" ;
         foreach my $sidsource (keys %{$all_edges}) {
            foreach my $sid12sig (keys %{$all_edges->{$sidsource}}) {
   
               my @sids = split(/\n/, $sid12sig) ;
               my ($sid1, $sid2) = split(/\n/, $sid12sig) ;
# TO PRINT sse topo for paper
#	 if ($sid1 ne 'BDP10759-0_SCOP.d1ivoa1' || $sid2 ne 'BDP10759-0_SCOP.d1ivoc_') {next;}
#	  print STDERR "printed for $bdp_id ($sid1--$sid2) ($sidsource)\n" ;
#	  my $coords = LGL::edges2coords({ edges => $all_edges->{$sidsource}->{$sid12sig} }) ;
#          LGL::coords2ps({
#             edges => $all_edges->{$sidsource}->{$sid12sig},
#             coords => $coords,
#	     mag => 300,
#	     nrad => 15,
#	     ethick => 2,
#             estyles => $all_estyles->{$sidsource}->{$sid12sig},
#	     ncols => $all_ncols->{$sidsource}->{$sid12sig},
#	     nshapes => $all_nshapes->{$sidsource}->{$sid12sig}
#          }) ;
   
# determine and print binding site topologies
               foreach my $j ( 0 .. $#sids) {
                  my @degs ;
                  my @degs_na ;
                  my $degs;
                  my @degconns ;
                  my $h1 =$bs_edges_dup->{$sidsource}->{$sid12sig}->{$sids[$j]};
                  my $numnodes = keys %{$h1} ;

                  foreach my $sse1 (keys %{$h1}) {
                     my $deg = keys %{$h1->{$sse1}} ;
                     $degs->{$sse1} = $deg ;
                     push @degs, $deg ;
                     
                     my $actsse1 = $ssedef->{$sse1} ;
                     if ($actsse1 eq '' ||
                         $actsse1 eq ' ') { $actsse1 = "_"; }
   
                     push @degs_na, $deg.$actsse1 ;
                  }
                  my @sd = sort {$b <=> $a} @degs ;
                  my $node_fp = join('.', @sd) ;
                  
                  my @sd_na = sort {$b cmp $a} @degs_na ;
                  my $node_fp_na = join('.', @sd_na) ;
                  
                  my @edges;
                  my @edges_na;
                  my $h2 = $bs_edges->{$sidsource}->{$sid12sig}->{$sids[$j]};
                  my $numedges = 0 ;
                  foreach my $sse1 (keys %{$h2}) {
                     my $actsse1 = $ssedef->{$sse1} ;
                     if ($actsse1 eq '' || $actsse1 eq ' ') {
                        $actsse1 = "_"; }
   
                     $numedges += keys %{$h2->{$sse1}} ;
                     foreach my $sse2 (keys %{$h2->{$sse1}}) {
                        my $actsse2 = $ssedef->{$sse2} ;

                        if ($actsse2 eq '' || $actsse2 eq ' ') {
                           $actsse2 = "_"; }
   
                        my @tdegs = ($degs->{$sse1}, $degs->{$sse2}) ;
                        my @stdegs = sort {$b <=> $a} @tdegs ;
                        push @edges, join("-", @stdegs) ;
                        
                        my @tdegs_na = ($degs->{$sse1}.$actsse1,
                                        $degs->{$sse2}.$actsse2) ;
                        my @stdegs_na = sort {$b cmp $a} @tdegs_na ;
                        push @edges_na, join("-", @stdegs_na) ;
                     }
                  }
                  my @sedges = sort @edges ;
                  my $edge_fp = join('.', @sedges) ;
                  
                  my @sedges_na = sort @edges_na ;
                  my $edge_fp_na = join('.', @sedges_na) ;
                  
                  print join("\t",
                     ( "BS", $bdp_id,$sids[0],$sids[1],$sids[$j],
                       $numnodes, $numedges,
                       $node_fp, $node_fp_na, $edge_fp, $edge_fp_na))."\n" ;
               }
   
   # determine and print interface topologies
               {
               my @degs ;
               my @degs_na ;
               my $degs;
               my @degconns ;
               my $h1 = $i_edges_dup->{$sidsource}->{$sid12sig};
               my $numnodes = keys %{$h1} ;
               foreach my $sse1 (keys %{$h1}) {
                  my $deg = keys %{$h1->{$sse1}} ;
                  $degs->{$sse1} = $deg ;
                  push @degs, $deg ;
                  my $actsse1 = $ssedef->{$sse1} ;

                  if ($actsse1 eq '' || $actsse1 eq ' ') {
                     $actsse1 = "_"; }
                  push @degs_na, $deg.$actsse1 ;
               }
               my @sd = sort {$b <=> $a} @degs ;
               my $node_fp = join('.', @sd) ;
               
               my @sd_na = sort {$b cmp $a} @degs_na ;
               my $node_fp_na = join('.', @sd_na) ;
               
               my @edges;
               my @edges_na;
               my $h2 = $i_edges->{$sidsource}->{$sid12sig};
               my $numedges = 0 ;
               foreach my $sse1 (keys %{$h2}) {
                  my $actsse1 = $ssedef->{$sse1} ;
                  if ($actsse1 eq '' || $actsse1 eq ' ') {
                     $actsse1 = "_"; }
                  $numedges += keys %{$h2->{$sse1}} ;
                  foreach my $sse2 (keys %{$h2->{$sse1}}) {
                     my $actsse2 = $ssedef->{$sse2} ;
                     if ($actsse2 eq '' || $actsse2 eq ' ') {
                        $actsse2 = "_"; }
                     my @tdegs = ($degs->{$sse1}, $degs->{$sse2}) ;
                     my @stdegs = sort {$b <=> $a} @tdegs ;
                     push @edges, join("-", @stdegs) ;

                     my @tdegs_na = ($degs->{$sse1}.$actsse1,
                                     $degs->{$sse2}.$actsse2) ;
                     my @stdegs_na = sort {$b cmp $a} @tdegs_na ;
                     push @edges_na, join("-", @stdegs_na) ;
                  }
               }
               my @sedges = sort @edges ;
               my $edge_fp = join('.', @sedges) ;

               my @sedges_na = sort @edges_na ;
               my $edge_fp_na = join('.', @sedges_na) ;
               
               print join("\t",
                  ("INT",$bdp_id,$sids[0],$sids[1], $numnodes, $numedges,
                   $node_fp, $node_fp_na, $edge_fp, $edge_fp_na))."\n" ;
               }
#	 my $coords = LGL::edges2coords({ edges => $all_edges->{$sidsource}->{$sid12sig} }) ;
#         LGL::coords2ps({
#            edges => $all_edges->{$sidsource}->{$sid12sig},
#            coords => $coords,
#	    mag => 300,
#	    nrad => 15,
#	    ethick => 2,
#            estyles => $all_estyles->{$sidsource}->{$sid12sig},
#	    ncols => $all_ncols->{$sidsource}->{$sid12sig},
#	    nshapes => $all_nshapes->{$sidsource}->{$sid12sig}
#         }) ;

            }
         }
   
#      my @degs ;
#      my $degs ;
#      my @degconns ;
#      foreach my $dom (keys %{$inthash->{$bdp_id}}) {
#         my $deg = keys %{$inthash->{$bdp_id}->{$dom}} ;
#	 $degs->{$dom} = $deg ;
#         push @degs, $deg ;
#      }
#      my @sd = sort {$b <=> $a} @degs ;
#      my $node_fp = join('.',@sd) ;
#
#      my @edges ;
#      foreach my $sidsig (@{$interfaces->{$bdp_id}}) {
#         my ($dom1, $dom2) = split(/\n/, $sidsig) ;
#	 my @tdegs = ($degs->{$dom1}, $degs->{$dom2}) ;
#	 my @stdegs = sort {$b <=> $a} @tdegs ;
#	 push @edges, join("-",@stdegs) ;
#      }
#      my @sedges = sort @edges ;
#      my $edge_fp = join('.', @sedges) ;
#
#
#      print join(" ", ($bdp_id, $bdp2numdom->{$bdp_id},
#                 $bdp2numdomtype->{$bdp_id}, 
#		 ($#{$interfaces->{$bdp_id}} + 1 ), $node_fp, $edge_fp))."\n" ;
      }
   }

   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_residues_tables file (specified in specs) into pibase
      $import_status->{interface_sse_topology} =
         pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{interface_sse_topology}
      }) ;
      $import_status->{bindingsite_sse_topology} =
         pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{bindingsite_sse_topology}
      }) ;
   }

   return $import_status ;

}


=head2 calc_bdp_interaction_topology()

   Title:       calc_bdp_interaction_topology()
   Function:    Calculates topology of domain interactions in a PIBASE complex
                  for each type of domain classification

=cut

sub calc_bdp_interaction_topology {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }

   my ($dbh) = pibase::connect_pibase() ;

   my ($bdps, $sid1, $sid2) = pibase::mysql_fetchcols($dbh,
      "SELECT bdp_id, subset_id_1, subset_id_2 ".
      "FROM intersubset_contacts WHERE num_contacts >= 50") ;

   my ($bdp2source2numdom) = pibase::mysql_hash2load($dbh,
      "SELECT bdp_id, subset_source_id, count(subset_id) FROM subsets ".
      "WHERE bdp_id IS NOT NULL GROUP BY bdp_id, subset_source_id") ;

   my ($bdp2source2numdomtype) = pibase::mysql_hash2load($dbh,
      "SELECT bdp_id, subset_source_id, count(distinct class) FROM subsets ".
      "WHERE bdp_id IS NOT NULL GROUP BY bdp_id, subset_source_id") ;

   my ($sid2source) = pibase::mysql_hashload($dbh,
      "SELECT subset_id, subset_source_id FROM subsets ".
      "WHERE bdp_id IS NOT NULL") ;

   my $interfaces ;
   my $inthash ;
   foreach my $j ( 0 .. $#{$sid1}) {
      my $source_id = $sid2source->{$sid1->[$j]} ;

      $inthash->{$bdps->[$j]}->{$source_id}->{$sid1->[$j]}->{$sid2->[$j]}++ ;
      $inthash->{$bdps->[$j]}->{$source_id}->{$sid2->[$j]}->{$sid1->[$j]}++ ;
      push @{$interfaces->{$bdps->[$j]}->{$source_id}}, $sid1->[$j]."\n".$sid2->[$j] ;
   }

   open(OUTF, ">".$pibase_specs->{buildfiles}->{bdp_interaction_topology});
   foreach my $bdp_id (keys %{$inthash}) {
    foreach my $source_id (keys %{$inthash->{$bdp_id}}) {
      my @degs ;
      my $degs ;
      my @degconns ;
      foreach my $dom (keys %{$inthash->{$bdp_id}->{$source_id}}) {
         my $deg = keys %{$inthash->{$bdp_id}->{$source_id}->{$dom}} ;
         $degs->{$dom} = $deg ;
         push @degs, $deg ;
      }
      my @sd = sort {$b <=> $a} @degs ;
      my $node_fp = join('.',@sd) ;
#      print STDERR "nodestring is $node_fp\n ";

      my @edges ;
      foreach my $sidsig (@{$interfaces->{$bdp_id}->{$source_id}}) {
         my ($dom1, $dom2) = split(/\n/, $sidsig) ;
         my @tdegs = ($degs->{$dom1}, $degs->{$dom2}) ;
         my @stdegs = sort {$b <=> $a} @tdegs ;
         push @edges, join("-",@stdegs) ;
      }
      my @sedges = sort @edges ;
      my $edge_fp = join('.', @sedges) ;


#      print STDERR "  nodestring is now $node_fp\n ";
#      print STDERR join("\t", ($bdp_id,
      print OUTF join("\t", ($bdp_id, $source_id,
                 $bdp2source2numdom->{$bdp_id}->{$source_id},
                 $bdp2source2numdomtype->{$bdp_id}->{$source_id},
                 ($#{$interfaces->{$bdp_id}->{$source_id}} + 1 ),
                 $node_fp, $edge_fp))."\n" ;
    }
   }
   close(OUTF) ;
   

   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_residues_tables file (specified in specs) into pibase
      $import_status->{interface_resvector} =
         pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{bdp_interaction_topology}
      }) ;
   }

   return $import_status ;

}


=head2 calc_interface_size()

   Title:       calc_interface_size()
   Function:    Calculates the number of residues from each domain at
                  each PIBASE interface

=cut

sub calc_interface_size {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }

   my $bdp2sid12 ;
   {
      my ($t_bdp, $t_sid1, $t_sid2) = pibase::rawselect_tod(
         "SELECT bdp_id, subset_id_1, subset_id_2 FROM intersubset_contacts") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2sid12->{$t_bdp->[$j]}->{$t_sid1->[$j]."\t".$t_sid2->[$j]}++ ; }
   }

   my $bdp2intcon_fn;
   {
      my ($t_bdp, $t_intcon_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM interface_contacts_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2intcon_fn->{$t_bdp->[$j]} = $t_intcon_fn->[$j]  ; }
   }

   my $bdp2subsres_fn;
   {
      my ($t_bdp, $t_subsres_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM subsets_residues_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2subsres_fn->{$t_bdp->[$j]} = $t_subsres_fn->[$j]  ; }
   }

   my $subset_sizes;
   {
      my ($t_sid, $t_length) = pibase::rawselect_tod(
         "SELECT subset_id, num_res FROM subsets_sequence") ;
      foreach my $j ( 0 .. $#{$t_sid}) {
         $subset_sizes->{$t_sid->[$j]} = $t_length->[$j]  ; }
   }


   my $bdp_list ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line =~ /^\#/) {next;}
         chomp $line;
         my $bdp_id = $line ;
         $bdp_list->{$bdp_id}++ ;
      }
      close(INF) ;
   } else {
      my ($dbh) = pibase::connect_pibase() ;
      $bdp_list = pibase::mysql_hashload($dbh,
         "SELECT bdp_id, 1 FROM intersubset_contacts") ;
   }

   
   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      print "* calc_interface_size() ".localtime() if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{calc_interface_size_in},
       $temp_fn->{calc_interface_size_in}) =
         tempfile("splits_calc_interface_size_input.XXXXX");

      ($temp_fh->{calc_interface_size_out},
       $temp_fn->{calc_interface_size_out}) =
         tempfile("splits_calc_interface_size_SGEout_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_size_out}) ;

      ($temp_fh->{calc_interface_size_err},
       $temp_fn->{calc_interface_size_err}) =
         tempfile("splits_calc_interface_size_SGEerr_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_interface_size_err}) ;

      foreach my $bdp (sort {$a <=> $b} keys %{$bdp_list}) {
         print {$temp_fh->{calc_interface_size_in}} $bdp."\n" ;
      }
      close($temp_fh->{calc_interface_size_in}) ;

      my $split_dir = tempdir("splits_calc_interface_size.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{calc_interface_size_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.calc_interface_size.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/calc_interface_size/ ;

main() ;

sub main {

         pibase::data::calc::calc_interface_size({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.calc_interface_size.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.calc_interface_size.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id}) ;
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{calc_interface_size_out},
         err_fn => $temp_fn->{calc_interface_size_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{calc_interface_size_out},
           $temp_fn->{calc_interface_size_out}) ;
      open(REALOUTF,">".$pibase_specs->{buildfiles}->{interface_size});

      while (my $line = readline($temp_fh->{calc_interface_size_out})){
         if ($line =~ /^\#/) { next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{calc_interface_size_out}) ;
      close(REALOUTF) ;

   } else {

      my $tempdir = tempdir(CLEANUP => 1) ;
      chdir $tempdir ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         print STDERR "NOW ON: bdp $bdp_id\n" ;
         if (!exists $bdp2sid12->{$bdp_id}) {
            print STDERR "WARNING: no interfaces in $bdp_id\n ";
            next;
         }

         my $interface_size ;
         my $intres ;
         {
         my ($t_subset_id_1, $t_subset_id_2,
             $t_chain_id_1, $t_chain_no_1, $t_resno_1,
             $t_chain_id_2, $t_chain_no_2, $t_resno_2 ) =
             pibase::rawselect_metatod($bdp2intcon_fn->{$bdp_id},
           "SELECT subset_id_1, subset_id_2, chain_id_1, chain_no_1, resno_1,".
           "chain_id_2, chain_no_2, resno_2 FROM $bdp2intcon_fn->{$bdp_id}") ;

         my $t_intsize ;
         foreach my $j ( 0 .. $#{$t_subset_id_1}) {
            my $sid12 = $t_subset_id_1->[$j]."\t".$t_subset_id_2->[$j] ;
            my $res1 = $t_resno_1->[$j]."\n".$t_chain_no_1->[$j] ;
            my $res2 = $t_resno_2->[$j]."\n".$t_chain_no_2->[$j]  ;
            $t_intsize->{$sid12}->{sub1}->{$res1}++ ;
            $t_intsize->{$sid12}->{sub2}->{$res2}++ ;
            $intres->{$sid12}->{$t_subset_id_1->[$j]}->{$t_resno_1->[$j].":".$t_chain_id_1->[$j]}++ ;
            $intres->{$sid12}->{$t_subset_id_2->[$j]}->{$t_resno_2->[$j].":".$t_chain_id_2->[$j]}++ ;
         }
         
         foreach my $sid12 (keys %{$t_intsize}) {
            my $a = keys %{$t_intsize->{$sid12}->{sub1}} ;
            my $b = keys %{$t_intsize->{$sid12}->{sub2}} ;
            $interface_size->{$sid12} = [ $a, $b ] ;
         }
         }

         foreach my $sid12 (keys %{$bdp2sid12->{$bdp_id}}) {
            my ($sid1, $sid2) = split(/\t/, $sid12) ;
            if (!exists $subset_sizes->{$sid1}) {
               print STDERR "ERROR: no subset size for $sid1\n"; next;}

            if (!exists $subset_sizes->{$sid2}) {
               print STDERR "ERROR: no subset size for $sid2\n"; next;}

            if (!exists $interface_size->{$sid12}) {
               print STDERR "ERROR: no interface_size computed for bdp $bdp_id:".
                            " $sid1 - $sid2\n"; next;}

            my ($num_res_1, $num_res_2) = @{$interface_size->{$sid12}} ;

            my @outvals = ( $bdp_id, $sid1, $sid2, $num_res_1, $num_res_2,
                            $subset_sizes->{$sid1}, $subset_sizes->{$sid2},
                            join(",", keys %{$intres->{$sid12}->{$sid1}}),
                            join(",", keys %{$intres->{$sid12}->{$sid2}})) ;
            print join("\t", @outvals)."\n" ;
         }
      }
   }

   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_residues_tables file (specified in specs) into pibase
      $import_status->{subsets_size} = pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{interface_size}}) ;
   }

   return $import_status ;

}


sub OLD_int_dsasa_calc_main {

   my $modeller_bin = "modSVN" ;

   my $usage = "perl ".__FILE__." < bdpid_list" ;
   if ($#ARGV >= 0) {die $usage;}

   my $bdpids ;
   while (my $line = <STDIN>) {
      chomp $line ;
      if ($line =~ /^#/) {next;}
      $bdpids->{$line}++ ;
   }

   my $bdp2sid ;
   print STDERR "\tsubsets load: " ;
   {
      my ($tsid, $tclass, $tssid, $tbdp)= pibase::rawselect_tod(
         "SELECT subset_id, class, subset_source_id, bdp_id FROM subsets") ;
      foreach my $j ( 0 .. $#{$tsid}) {
         if (!defined $tbdp->[$j] || $tbdp->[$j] eq '' || $tssid->[$j] != 1) {
            next;}
         $bdp2sid->{$tbdp->[$j]}->{$tsid->[$j]}++ ;
      }
   }
   print STDERR "X\n" ;

   my $bdp2sid12 ;
   print STDERR "\tinterface load: " ;
   {
      my ($tbdp, $tsid1, $tsid2)= pibase::rawselect_tod(
         "SELECT bdp_id, subset_id_1, subset_id_2 FROM intersubset_contacts") ;
      foreach my $j ( 0 .. $#{$tbdp}) {
         if ($tsid1->[$j] !~ /SCOP/) {next;}
         $bdp2sid12->{$tbdp->[$j]}->{$tsid1->[$j]."\n".$tsid2->[$j]}++ ;
      }
   }
   print STDERR "X\n" ;

# get a list of bdp ids to take care of
# reads constituent subsets, interactions from tables on disk

#Read in a bdp file path from STDIN

   foreach my $bdp (sort {$a <=> $b} keys %{$bdpids}) {
      print STDERR "NOW ON: bdp $bdp\n" ;

      my $sid2sasa ;
      foreach my $sid (sort keys %{$bdp2sid->{$bdp}}) {
         print STDERR "NOW ON: bdp $bdp sid $sid\n" ;
         my $sid_fn = pibase::sid_2_domdir($sid)."/$sid.pdb" ;
         if (!-s $sid_fn) {
            print STDERR "ERROR: can't find $sid file $sid_fn\n ";
            next;
         }
         my ($res_sasa, $sasa_resno_rev, $fullsasa, $sasa_error_fl,
             $atmsasa_total) = pibase::modeller::calc_sasa({
            pdb_fn => $sid_fn,
            surftyp => 2,
            modeller_bin => $modeller_bin
         }) ;

         if ($#{$sasa_error_fl} >= 0 ) {
            foreach my $j ( 0 .. $#{$sasa_error_fl}) {
               print STDERR "ERROR: subset ($bdp $sid): calc_sasa() $sasa_error_fl->[$j]" ; } next; }
         $sid2sasa->{$sid}->{res_sasa} = $res_sasa ;
         $sid2sasa->{$sid}->{sasa_resno_rev} = $sasa_resno_rev;
         $sid2sasa->{$sid}->{fullsasa} = $fullsasa ;
         $sid2sasa->{$sid}->{atmsasa_total} = $atmsasa_total;

         my @outvals = ("SUBSET_SASA", $bdp, $sid, 
            sprintf("%.3f", $sid2sasa->{$sid}->{fullsasa}->{all}),
            sprintf("%.3f", $sid2sasa->{$sid}->{fullsasa}->{sc}),
            sprintf("%.3f", $sid2sasa->{$sid}->{fullsasa}->{mc}),
            sprintf("%.3f", $sid2sasa->{$sid}->{atmsasa_total}->{p}),
            sprintf("%.3f", $sid2sasa->{$sid}->{atmsasa_total}->{nonp})
         ) ;
         print join("\t", @outvals)."\n" ;
      }

      foreach my $sid12 (keys %{$bdp2sid12->{$bdp}}) {
         my ($sid1, $sid2) = split(/\n/, $sid12) ;
         print STDERR "NOW ON: bdp $bdp interface $sid1 -- $sid2\n" ;
         my ($pairpdb_fh, $pairpdb_fn) =
            tempfile("pairpdb.$bdp.XXXXXX", SUFFIX=>"pdb") ;
         close($pairpdb_fh) ;

         my $sid1_fn = pibase::sid_2_domdir($sid1)."/$sid1.pdb" ;
         my $sid2_fn = pibase::sid_2_domdir($sid2)."/$sid2.pdb" ;
         system("cat $sid1_fn $sid2_fn > $pairpdb_fn") ;
         if (!-s $pairpdb_fn) {
            print STDERR "ERROR: couldnt make complex pdb file: $sid1, $sid2\n";
            next; }

         my ($res_sasa, $sasa_resno_rev, $fullsasa, $sasa_error_fl,
             $atmsasa_total) = get_sasa({
            pdb_fn => $pairpdb_fn,
            surftyp => 2,
            modeller_bin => $modeller_bin
         }) ;
         unlink $pairpdb_fn ;
         
         if ($#{$sasa_error_fl} >= 0 ) {
            foreach my $j ( 0 .. $#{$sasa_error_fl}) {
               print STDERR "ERROR: interaction ($bdp $sid1 -- $sid2): get_sasa() $sasa_error_fl->[$j]" ; } next; }

         my $dsasa ;
         $dsasa->{all} = $sid2sasa->{$sid1}->{fullsasa}->{all} +
            $sid2sasa->{$sid2}->{fullsasa}->{all} - $fullsasa->{all} ;
         
         $dsasa->{sc} = $sid2sasa->{$sid1}->{fullsasa}->{sc} +
            $sid2sasa->{$sid2}->{fullsasa}->{sc} - $fullsasa->{sc} ;

         $dsasa->{mc} = $sid2sasa->{$sid1}->{fullsasa}->{mc} +
            $sid2sasa->{$sid2}->{fullsasa}->{mc} - $fullsasa->{mc} ;

         $dsasa->{nonp} = $sid2sasa->{$sid1}->{atmsasa_total}->{nonp} +
           $sid2sasa->{$sid2}->{atmsasa_total}->{nonp} - $atmsasa_total->{nonp};

         $dsasa->{p} = $sid2sasa->{$sid1}->{atmsasa_total}->{p} +
           $sid2sasa->{$sid2}->{atmsasa_total}->{p} - $atmsasa_total->{p};

         my @outvals = ("INTERFACE_dSASA", $bdp, $sid1, $sid2,
                        sprintf("%.3f", $dsasa->{all}),
                        sprintf("%.3f", $dsasa->{sc}),
                        sprintf("%.3f", $dsasa->{mc}),
                        sprintf("%.3f", $dsasa->{p}),
                        sprintf("%.3f", $dsasa->{nonp})
                     ) ;
         print join("\t", @outvals)."\n" ;
      }
   }

}



#OLD CODE not included in the overhaul of 2008
sub int_sc_calc_main {

   my $aares = {
      'ALA' => 'A' ,
      'ARG' => 'R' ,
      'ASN' => 'N' ,
      'ASP' => 'D' ,
      'CYS' => 'C' ,
      'GLN' => 'Q' ,
      'GLU' => 'E' ,
      'GLY' => 'G' ,
      'HIS' => 'H' ,
      'HSD' => 'H' ,
      'ILE' => 'I' ,
      'LEU' => 'L' ,
      'LYS' => 'K' ,
      'MET' => 'M' ,
      'PHE' => 'F' ,
      'PRO' => 'P' ,
      'SER' => 'S' ,
      'THR' => 'T' ,
      'TRP' => 'W' ,
      'TYR' => 'Y' ,
      'VAL' => 'V',
   } ;

#Set binary locations.

   my $modeller_bin = "modSVN" ;
   my $binaries = pibase::locate_binaries() ;
   my $altloc_check = $binaries->{'altloc_check'} ;
   my $altloc_filter = $binaries->{'altloc_filter'} ;

   my $usage = "perl ".__FILE__." < bdpid_list" ;
   if ($#ARGV >= 0) {die $usage;}


   my $bdp2sid12;
   while (my $line = <STDIN>) {
      chomp $line ;
      if ($line =~ /^#/) {next;}
      my ($tbdp, $tsid1, $tsid2) = split(/\t/, $line) ;
      my $tsid12 = $tsid1."\n".$tsid2 ;
      $bdp2sid12->{$tbdp}->{$tsid12}++ ;
   }

   my $bdp2sid ;
   print STDERR "\tsubsets load: " ;
   {
      my ($tsid, $tclass, $tssid, $tbdp)= pibase::rawselect_tod(
         "SELECT subset_id, class, subset_source_id, bdp_id FROM subsets") ;
      foreach my $j ( 0 .. $#{$tsid}) {
         if (!defined $tbdp->[$j] || $tbdp->[$j] eq '' || $tssid->[$j] != 1) {
            next;}
         $bdp2sid->{$tbdp->[$j]}->{$tsid->[$j]}++ ;
      }
   }
   print STDERR "X\n" ;

# get a list of bdp ids to take care of
# reads constituent subsets, interactions from tables on disk

#Read in a bdp file path from STDIN

   foreach my $bdp (sort {$a <=> $b} keys %{$bdp2sid12}) {
      print STDERR "NOW ON: bdp $bdp\n" ;

      foreach my $sid12 (sort keys %{$bdp2sid12->{$bdp}}) {
         my ($sid1, $sid2) = split(/\n/, $sid12) ;
         print STDERR "NOW ON: bdp $bdp interface $sid1 -- $sid2\n" ;

         my ($pairpdb_fh, $pairpdb_fn) =
            tempfile("pairpdb.$bdp.XXXXXX", SUFFIX=>"pdb") ;

         my $sid1_fn = pibase::sid_2_domdir($sid1)."/$sid1.pdb" ;
         my $sid2_fn = pibase::sid_2_domdir($sid2)."/$sid2.pdb" ;
         open(SID1, $sid1_fn) ;
         while (my $line = <SID1>) {
            chomp $line;
            if ($line !~ /^ATOM/ ||
                (substr($line, 13, 1) eq 'H') ||
                (substr($line, 13, 1) eq 'Q')) {next;}
            my $resna = substr($line,17,3) ;
            if (!exists $aares->{$resna}) {next;}
            substr($line,21,1) = 'A' ;
            print $pairpdb_fh $line."\n" ;
         }
         close(SID1) ;

         open(SID2, $sid2_fn) ;
         while (my $line = <SID2>) {
            chomp $line;
            if ($line !~ /^ATOM/ ||
                (substr($line, 13, 1) eq 'H') ||
                (substr($line, 13, 1) eq 'Q')) {next;}
            my $resna = substr($line,17,3) ;
            if (!exists $aares->{$resna}) {next;}
            substr($line,21,1) = 'B' ;
            print $pairpdb_fh $line."\n" ;
         }
         close(SID2) ;
         close($pairpdb_fh) ;

         if (!-s $pairpdb_fn) {
            print STDERR "ERROR: couldnt make complex pdb file: $sid1, $sid2\n";
            next;
         }

         my $altloc_fl = `$altloc_check < $pairpdb_fn` ; chomp $altloc_fl ;
         if ($altloc_fl) {
            print STDERR "NOTE: ($bdp) $sid1 - $sid2 altloc filtering\n" ;
            my ($tfh, $tfn) = tempfile("altloctemp.XXXXX",SUFFIX=>".pdb") ;
            close($tfh) ;
            my $filtcom = "$altloc_filter $pairpdb_fn > $tfn" ;
            system($filtcom) ;
            system("mv $tfn $pairpdb_fn") ;
         }

         my $scinfo = int_sc_calc_calc_sc({pdb_fn => $pairpdb_fn,
                               sc_bin => $binaries->{ccp4sc},
                               chain1 => 'A',
                               chain2 => 'B'}) ;
         unlink $pairpdb_fn ;


         if (exists $scinfo->{error_fl}) {
            print STDERR "ERROR: $bdp $sid1 $sid2 threw calc_sc() error:".
               $scinfo->{error_fl}."\n";
         }

         my @outvals = ($bdp, $sid1, $sid2,
                        $scinfo->{'Sc'}, $scinfo->{'mediandist'}) ;
         foreach my $j ( 3 .. 4) {
            if (!defined $outvals[$j]) {
               $outvals[$j] = "";
            } else {
               $outvals[$j] = sprintf("%.3f", $outvals[$j]) ;
            }
         }
         print join("\t", @outvals)."\n" ;
      }
   }

}

#OLD CODE not included in the overhaul of 2008
sub int_sc_calc_calc_sc {
   my $in = shift ;
   my ($tempin_fh, $tempin_fn) = tempfile("XXXXX", SUFFIX =>".scin") ;
   if (!exists $in->{sc_bin}) {
      die "FATAL ERROR: CCP4 Sc binary not specified\n" ;}

   print $tempin_fh "MOLECULE 1\nCHAIN $in->{chain1}\n" ;
   print $tempin_fh "MOLECULE 2\nCHAIN $in->{chain2}\n" ;
   close ($tempin_fh) ;

   my ($tempout_fh, $tempout_fn) = tempfile("XXXXX", SUFFIX =>".scout") ;

   system("$in->{sc_bin} $in->{pdb_fn} < $tempin_fn > $tempout_fn") ;
   unlink $tempin_fn ;

   my $scinfo ;
   if (!-s $tempout_fn) {
      if (-e $tempout_fn) {unlink $tempout_fn;}
      return {error_fl => "Sc didnt run properly" } ;
   } else {
      open(SCOUT, $tempout_fn) ;
      my ($sc, $mediandist) ;
      while (my $line = <SCOUT>) {
         chomp $line;
         if ($line =~ /Summary of results/) {
            my $sumline ;
            $sumline = <SCOUT> ;
            $sumline = <SCOUT> ;
            $sumline = <SCOUT> ;
            $sumline = <SCOUT> ;
            my @t = split(' ',$sumline)  ; $scinfo->{mediandist} = $t[3] ;
            $sumline = <SCOUT> ;
            $sumline = <SCOUT> ;
            $sumline = <SCOUT> ;
            $sumline = <SCOUT> ;
            @t = split(' ',$sumline)  ; $scinfo->{Sc} = $t[3] ;
            last;
         }
      }
      unlink $tempout_fn ;
   }

   return $scinfo ;
}


#OLD CODE not included in the overhaul of 2008
sub int_planarity_main {

   my $binaries = pibase::locate_binaries() ;

   my $bdp2contacts_fn;
   print STDERR "\ninterface_contacts_tables load: " ;
   {
      my ($tbdp, $tcontacts_fn)= pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM interface_contacts_tables") ;
      foreach my $j ( 0 .. $#{$tbdp}) {
         $bdp2contacts_fn->{$tbdp->[$j]} = $tcontacts_fn->[$j] ;
      }
   }
   print STDERR "X\n" ;

   my $count = 1 ;
   while (my $line = <STDIN>) {
      if ($line =~ /^#/) {next;}
      chomp $line;
      my $bdp_id = $line ;
      my $contacts_fn = $bdp2contacts_fn->{$bdp_id} ;
      print STDERR "NOW ON: $bdp_id\n ";
      if (!-s $contacts_fn) {
         print STDERR "ERROR: bdp $bdp_id could not find $contacts_fn\n" ;
         next;
      }

      my ($sid1, $chain1, $resno1, $resna1, $sid2, $chain2,
          $resno2, $resna2) = pibase::rawselect_metatod(
          $contacts_fn,
          "SELECT subset_id_1, chain_id_1, resno_1, resna_1,".
          "subset_id_2, chain_id_2, resno_2, resna_2 ".
          "FROM $contacts_fn") ;

      my $intres ;
      foreach my $j ( 0 .. $#{$sid1}) {
         if ($sid1->[$j] !~ /SCOP/) {next;}
         my $sid12 = join("\n", sort ($sid1->[$j], $sid2->[$j])) ;
         $intres->{$sid12}->{$resno1->[$j]."\n".$chain1->[$j]."\n".$resna1->[$j]}++ ;
         $intres->{$sid12}->{$resno2->[$j]."\n".$chain2->[$j]."\n".$resna2->[$j]}++ ;
      }
      my $numints = keys %{$intres} ;
      if ($numints <= 0) {
         print STDERR "NOTE: no interface in $bdp_id\n" ;
         next;
      }

      foreach my $sid12 (keys %{$intres}) {
         my ($cursid1, $cursid2) = split(/\n/, $sid12) ;
         print STDERR "NOW ON: $bdp_id $cursid1 $cursid2\n" ;

         my $sid1_fn = pibase::sid_2_domdir($cursid1)."/$cursid1.pdb" ;
         my $sid2_fn = pibase::sid_2_domdir($cursid2)."/$cursid2.pdb" ;
         if (!-s $sid1_fn) {
            print STDERR "ERROR: can't find $cursid1 file $sid1_fn\n ";
            next;
         }

         if (!-s $sid2_fn) {
            print STDERR "ERROR: can't find $cursid2 file $sid2_fn\n ";
            next;
         }

         my ($temppdb_fh, $temppdb_fn) = tempfile("interfaceonly.XXXXX", SUFFIX => ".pdb") ;
         foreach my $sidfn ($sid1_fn, $sid2_fn) {
            open(SIDF, $sidfn) ;
            while (my $pdbline = <SIDF>) {
               chomp $pdbline;
               if ($pdbline =~ /^ATOM/) {
                  my $curchain = substr($pdbline, 21, 1) ;
                  my $curresno = substr($pdbline, 22, 5) ;
                  $curresno =~ s/ //g ;
                  my $curresna = substr($pdbline, 17, 3) ;

                  my $ressig = $curresno."\n".$curchain."\n".$curresna ;

                  if (exists $intres->{$sid12}->{$ressig}) {
                     $intres->{$sid12}->{$ressig} = 'DONE';
                     print $temppdb_fh $pdbline."\n" ;
                  }
               } elsif ($pdbline =~ /^ENDMDL/) {
                  last;
               }
            }
            close(SIDF) ;
         }
         close($temppdb_fh) ;
         if (!-s $temppdb_fn) {
            print STDERR "ERROR: $bdp_id $cursid1 $cursid2 could not make complex pdb file()\n" ;
            next ;
         }

         my $princip = int_planarity_calc_princip({pdb_fn => $temppdb_fn,
                                     princip_bin => $binaries->{princip}});
         unlink $temppdb_fn ;

         if (exists $princip->{error_fl}) {
            print STDERR "ERROR: $bdp_id $cursid1 $cursid2 ".
                         "calc_princip(): $princip->{error_fl}\n" ;
            next ;
         }

         my @outvals = ($bdp_id, $cursid1, $cursid2, $princip->{rmsd}) ;
         print join("\t", @outvals)."\n" ;

         foreach my $res (keys %{$intres->{$sid12}}) {
            if ($intres->{$sid12}->{$res} ne 'DONE') {
               my $t = $res; $t =~ s/\n/:/g ;
               print STDERR "WARNING $bdp_id $cursid1 -- $cursid2: ".
                            "$t not found in PDB file\n" ;
            }
         }
      }
      $count++ ;
   }
}


#OLD CODE not included in the overhaul of 2008
sub int_planarity_calc_princip {

   my $in = shift ;
   my ($tempout_fh, $tempout_fn) = tempfile("principout.XXXXX") ;
   my ($temperr_fh, $temperr_fn) = tempfile("principerr.XXXXX") ;
   close($tempout_fh) ;
   my $tcom = "echo $in->{pdb_fn} | $in->{princip_bin} 2>$temperr_fn >$tempout_fn" ;
   system($tcom) ;

   if (-s $temperr_fn) {
      print STDERR "\nerror file $tcom\n";
      open(ERR, $temperr_fn) ;
      while (my $line = <ERR>) {
         print STDERR $line ;
      }
      close(ERR) ;
      print STDERR "\n------------\n" ;
   }
   unlink $temperr_fn ;

   my $results ;
   if (!-s $tempout_fn) {
      $results->{error_fl} = "PRINCIP run failed" ;
      print STDERR "princip run $tcom failed\n" ;
   } else {
      my $rmsd ;
      open (OUTF, $tempout_fn) ;
      while (my $outline =<OUTF>) {
         chomp $outline;
         if ($outline =~ /RMS difference from best-fit/) {
            $outline =~ s/.*plane:// ; $outline =~ s/ //g ;
            $rmsd = $outline ;
            last;
         }
      }
      close(OUTF) ;

      if (!defined $rmsd) {
         $results->{error_fl} = "PRINCIP did not report rmsd" ;
      } else {
         $results->{rmsd} = $rmsd ;
      }
   }

   unlink $tempout_fn ;

   return $results ;
}


=head2 calc_bdp_interaction_topology_graph()

   Title:       calc_bdp_interaction_topology_graph()
   Function:    Creates LGL graphs of domain interaction topology for each
                  PIBASE complex

=cut

sub calc_bdp_interaction_topology_graph {

   require LGL ;

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }
   my $binaries = pibase::locate_binaries() ;

# read in bdp_id, bdp_path info from input if specified, otherwise
#  get list from bdp_files, and split locally.

   my $bdp_list = {};
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line=~ /^#/) {next;}
         chomp $line;
         my ($bdp_id, $bdp_path) = split(/\t/, $line) ;
         $bdp_list->{$bdp_id}++ ;
      }
      close(INF) ;
   } else {
      my ($dbh) = pibase::connect_pibase() ;
      $bdp_list = pibase::mysql_hashload($dbh,
         "SELECT bdp_id, file_path FROM bdp_files") ;
   }

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
# if this is master script, split input, recall self with
# in_fn specified, cluster_fl = 0, import_fl = 0
# send out, cluster run, and return merged

      print "* calc_bdp_interaction_topology_graph() ".localtime()
         if (!exists $in->{quiet_fl});
      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{calc_bdp_interaction_topology_graph_in},
       $temp_fn->{calc_bdp_interaction_topology_graph_in}) =
         tempfile("splits_calc_bdp_interaction_topology_graph_input.XXXXX");
      ($temp_fh->{calc_bdp_interaction_topology_graph_out},
       $temp_fn->{calc_bdp_interaction_topology_graph_out}) =
         tempfile("splits_calc_bdp_interaction_topology_graph_SGEout_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_bdp_interaction_topology_graph_out}) ;
      ($temp_fh->{calc_bdp_interaction_topology_graph_err},
       $temp_fn->{calc_bdp_interaction_topology_graph_err}) =
         tempfile("splits_calc_bdp_interaction_topology_graph_SGEerr_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{calc_bdp_interaction_topology_graph_err}) ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         print {$temp_fh->{calc_bdp_interaction_topology_graph_in}}
            join("\t", $bdp_id, $bdp_list->{$bdp_id})."\n" ; }
      close($temp_fh->{calc_bdp_interaction_topology_graph_in}) ;

      my $split_dir = tempdir("splits_calc_bdp_interaction_topology_graph.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{calc_bdp_interaction_topology_graph_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.calc_bdp_interaction_topology_graph.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::calc qw/calc_bdp_interaction_topology_graph/ ;

main() ;

sub main {

         pibase::data::calc::calc_bdp_interaction_topology_graph({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.calc_bdp_interaction_topology_graph.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.calc_bdp_interaction_topology_graph.XXXXX");

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $pibase_specs->{SGE}->{priority}) {
print {$sgescript_fh} "#\$ -p $pibase_specs->{SGE}->{priority}\n"; }

      if (exists $pibase_specs->{SGE}->{nodespecs}) {
print {$sgescript_fh} $pibase_specs->{SGE}->{nodespecs}."\n"; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = pibase::SGE::_clust_qsub({
         sgescript_fn => $sgescript_fn,
      }) ;

      while (1) {
         sleep $pibase_specs->{SGE}->{qstat_sleep} ;
         my $job_status = pibase::SGE::_clust_qstat({job_id => $qsub_job_id}) ;
         if ($job_status) {last;}
      }

      pibase::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{calc_bdp_interaction_topology_graph_out},
         err_fn => $temp_fn->{calc_bdp_interaction_topology_graph_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{calc_bdp_interaction_topology_graph_out},
           $temp_fn->{calc_bdp_interaction_topology_graph_out}) ;
      open(REALOUTF,
         ">".$pibase_specs->{buildfiles}->{bdp_interaction_topology_graph}) ;

      while (my $line =
         readline($temp_fh->{calc_bdp_interaction_topology_graph_out})) {

         if ($line =~ /^\#/) {next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{calc_bdp_interaction_topology_graph_out}) ;
      close(REALOUTF) ;

   } else {

      my $movethreshold = 100 ;
      my $temp_dir = tempdir(CLEANUP => 1) ; chdir $temp_dir ;
      my (@movethese, @moveto) ;

      my ($sid2class, $sid2source);
      my $bdp2source2sid ;
      {
         my ($t_bdp, $t_sid, $t_class, $t_subset_source)=pibase::rawselect_tod(
            "SELECT bdp_id, subset_id, class, subset_source_id FROM subsets") ;
         foreach my $j ( 0 .. $#{$t_bdp}) {
            $sid2source->{$t_sid->[$j]} = $t_subset_source->[$j] ;
            $sid2class->{$t_sid->[$j]} = $t_class->[$j] ;
   $bdp2source2sid->{$t_bdp->[$j]}->{$t_subset_source->[$j]}->{$t_sid->[$j]}++ ;
         }
      }
   
      my $bdp2source2int ;
      {
         my ($t_bdp, $t_sid1, $t_sid2, $t_chains) = pibase::rawselect_tod(
            "SELECT bdp_id, subset_id_1, subset_id_2, chains ".
            "FROM intersubset_contacts");

         foreach my $j ( 0 .. $#{$t_bdp}) {
            my $t_source = $sid2source->{$t_sid1->[$j]} ;
$bdp2source2int->{$t_bdp->[$j]}->{$t_source}->{$t_sid1->[$j]}->{$t_sid2->[$j]} = $t_chains->[$j] ;
         }
      }
   
      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         my $depositdir = $pibase_specs->{otherdata_dir}->{bdp_topology_graphs}.
            '/'.POSIX::floor($bdp_id / 1000) ;
   
         foreach my $subset_source (keys %{$bdp2source2int->{$bdp_id}}) {
            print STDERR "NOW ON: $bdp_id (subset source $subset_source)\n" ;

            my $coords = LGL::edges2coords({
               edges => $bdp2source2int->{$bdp_id}->{$subset_source} }) ;

            my $class2col ;
            my $j = 0 ;
            foreach my $class (sort
@$sid2class{(keys %{$bdp2source2sid->{$bdp_id}->{$subset_source}})}) {

               if (exists $class2col->{$class}) {next;}
               $class2col->{$class}=$pibase_specs->{web}->{domain_colors}->[$j];
               $j++ ;
            }

# sort the classes, assign colors in order of class ascending;
# use same colors for domains listed on the html table for complexes.
   
            my $ncols ;
            foreach my $sid (
               keys %{$bdp2source2sid->{$bdp_id}->{$subset_source}}) {

               $ncols->{$sid} = $class2col->{$sid2class->{$sid}} ; }

            my $ecols ;
            foreach my $sid1 (keys
               %{$bdp2source2int->{$bdp_id}->{$subset_source}}) {
               foreach my $sid2 (keys
                  %{$bdp2source2int->{$bdp_id}->{$subset_source}->{$sid1}}) {
if ($bdp2source2int->{$bdp_id}->{$subset_source}->{$sid1}->{$sid2} eq 'same') {
                  $ecols->{$sid1}->{$sid2} = 'black' ;
} else {
                  $ecols->{$sid1}->{$sid2} = 'grey67' ;
}
               }
            }
   
            my $out_fn ;
            $out_fn->{eps} = "complex_topology_".$bdp_id."_".
               $subset_source.".eps" ;
            LGL::coords2ps({
               coords => $coords,
               edges => $bdp2source2int->{$bdp_id}->{$subset_source},
               ecols => $ecols,
               ncols => $ncols,
               nshape => "fcircle",
               nrad => 12,
               ethick => 4,
               mag => 200,
               out_fn => $out_fn->{eps}
            }) ;

            if (!-s $out_fn->{eps}) {
               print STDERR "ERROR: no graph created for bdp $bdp_id, ".
                            "error in LGL::coords2ps()\n";
               next;
            }

#            $out_fn->{png}="complex_topology_".$bdp_id."_".
#               $subset_source.".png" ;
#
#            system($binaries->{imagemagick_convert}." ".
#               $out_fn->{eps}." ".$out_fn->{png}) ;
#            unlink $out_fn->{eps} ;
#            if (!-s $out_fn->{png}) {
#               print STDERR "ERROR: no graph created for bdp $bdp_id, ".
#                  "error in convert to png\n";
#               next;
#            }
#   
#            my $im = LGL::coords2png({
#               coords => $coords,
#               edges => $bdp2source2int->{$bdp_id}->{$subset_source},
#               ecols => $ecols,
#               ncols => $ncols,
#               nshape => "fcircle",
#               nrad => 7,
#               ethick => 2.5,
#               mag => 400,
#               out_fn => $out_fn->{png}
#            }) ;
#
#            if (!-s $out_fn->{png}) {
#               print STDERR "ERROR: no png graph created for bdp $bdp_id\n";
#               next;
#            }
#
#            $out_fn->{png} = $depositdir.'/'.$out_fn->{png} ;
#            push @movethese, $out_fn->{png} ;
#            push @moveto, $depositdir ;

            push @movethese, $out_fn->{eps} ;
            push @moveto, $depositdir ;

            $out_fn->{eps} = $depositdir.'/'.$out_fn->{eps} ;
            my @outvals = ($bdp_id, $subset_source, $out_fn->{eps}) ;
            print join("\t", @outvals)."\n" ;
         }
   
         if ($#movethese == $movethreshold) {
            foreach my $j (0 .. $#movethese) {
               my $t_fn = $movethese[$j] ;
               my $t_dir = $moveto[$j] ;
               if (! -d $t_dir) { mkpath($t_dir) ;}
               pibase::safe_move($t_fn, $t_dir);
            }
            @movethese = () ; @moveto = () ;
         }
      }
   
      foreach my $j (0 .. $#movethese) {
         my $t_fn = $movethese[$j] ;
         my $t_dir = $moveto[$j] ;
         if (! -d $t_dir) { mkpath($t_dir) ; }
         pibase::safe_move($t_fn, $t_dir);
      }
      @movethese = () ; @moveto = () ;

   }

   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_secstrx_tables file (specified in specs) into pibase
      $import_status = pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{bdp_interaction_topology_graph}}) ;
   }

   return $import_status ;


}


1 ;
