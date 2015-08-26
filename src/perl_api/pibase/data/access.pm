=head1 NAME

pibase::data::access - perl module of pibase data access routines

=head1 DESCRIPTION

Perl package that provides pibase data access routines

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

package pibase::data::access ;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT = qw/bdp_path_2_id maketable_bdp_files create_subset_pdbs/ ;

use pibase qw/connect_pibase mysql_hashload mysql_fetchcols mysql_hasharrload safe_move sid_2_domdir/;
use pibase::PDB::subsets qw/subset_extract/;
use File::Temp qw/tempfile tempdir/ ;
use File::Path qw/mkpath/ ;
use File::Basename ;
use POSIX qw/ceil/ ;

=head2 bdp_path_2_id()

   Title:       bdp_path_2_id()
   Function:    determine bdp_ids of specified bdp file paths
   Args:        none
   Returns:     nothing
   STDIN:       bdp_path
   STDOUT:      bdp_id."\t".bdp_path

=cut

sub bdp_path_2_id {

   my ($dbh, $pibase) = pibase::connect_pibase() ;
   my ($path_2_bdp_id) = pibase::load_bdp_ids($dbh, qw/path_2_bdp_id/) ;

   while (my $bdp_path = <STDIN>) {
    if ($bdp_path !~ /^\#/) {
      chomp $bdp_path ;
      if (! exists $path_2_bdp_id->{$bdp_path}) {
         print STDERR "ERROR $bdp_path not found in bdp_files\n" ; }
      else {
         print $path_2_bdp_id->{$bdp_path}."\t".$bdp_path."\n" ; } } }

}


=head2 maketable_bdp_files()

   Title:       maketable_bdp_files()
   Function:    determine bdp_ids of specified bdp file paths
   Args:        none
   Returns:     nothing
   comand line: -list [optional] - expect STDIN list of pdb_ids
                  to process (defaults to everything in pdb_entries)
   STDIN:       pdb_id [optional - if -list option set.
   STDOUT:      bdp_id."\t".bdp_path
   Files out:   $in->{fn_out}
                  1. '0'
                  2. pisa file - full path
                  3. pisa file - base name
                  4. pdb_id
                  5. '0'

=cut

sub maketable_bdp_files {
   require DBI ;
   my $in = shift ;

# Connect to pibase.
   my ($dbh, $pibase) = pibase::connect_pibase() ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = $in->{pibase_specs} ;
   } else {
      $pibase_specs = pibase::get_specs() ;
   }

   my ($fh, $fn) ;
   if (!exists $in->{out_fn}) {
      $fn->{out} = $pibase_specs->{buildfiles}->{bdp_files} ;
      my (undef, $out_dir, undef)= File::Basename::fileparse(
         $pibase_specs->{buildfiles}->{bdp_files}) ;
      if (!-s $out_dir) {mkpath $out_dir;}
   } else {
      $fn->{out} = $in->{out_fn} ;
   }
   open($fh->{out}, ">".$fn->{out}) ;

# Load experiment type entries from pibase.pdb_entry_type.
   my $experiment = pibase::mysql_hashload($dbh,
      'SELECT pdb_id, experiment_type FROM pdb_entry_type');

   my @pdb_ids = sort keys %{$experiment};

# Load PISA entries from pibase.pisa_index.
   my $pdb2pisa = pibase::mysql_hasharrload($dbh,
      'SELECT pdb_id, pisa_id FROM pisa_index WHERE entry_type != "ERROR" AND '.
      "num_chain <= ".$pibase_specs->{external_data}->{pisa}->{max_num_chains});

   my @pdblist ;
   if (exists $in->{in_fn}) { # If input file specified, read in pdb_id list
      open(PDBLIST, $in->{in_fn}) ;
      while (my $pdb_id = <PDBLIST>) {
         if ($pdb_id =~ /^#/) {next;} # Next if a comment line
         chomp $pdb_id;

         if (!exists $experiment->{$pdb_id}) {
            print STDERR "ERROR: $pdb_id skipped: not found in pdb_entries\n";
         } else {
            push @pdblist, $pdb_id ;
         }
      }
      close(PDBLIST) ;
   } else { # Otherwise, load entire pdb_id list from pdb_entries in @pdblist
      @pdblist = @pdb_ids ;
   }

   print STDERR "now on 0" ;
   my $counter = 1 ;
   foreach my $pdb_id ( @pdblist ) { # Iterate through the pdb_id array
      print STDERR "\b"x(length($counter - 1))."$counter" ;

      if ($experiment->{$pdb_id} eq 'NMR') {

         my $modelpdb = $pdb_id."_1.ent.gz" ;
         my $full_path = $pibase_specs->{pdbnmr_dir}."/".substr($pdb_id,1,2).'/'.
                         $modelpdb ;
         if (!-s $full_path) {
            print STDERR "ERROR: PDB_NMR $pdb_id not found in $full_path\n"; }

         my @outvals = ('0', $full_path, $modelpdb, $pdb_id, '1') ;
         print {$fh->{out}} join("\t", @outvals)."\n" ;

      } else {

         my $full_path = pibase::PDB::get_pdb_filepath({pdb_id => $pdb_id,
                                            pibase_specs => $pibase_specs}) ;
         my $pdbfile = $full_path; $pdbfile =~ s/.*\/// ;
         if (-e $full_path) {
            my @outvals = ('0', $full_path, $pdbfile , $pdb_id, '1') ;
            print {$fh->{out}} join("\t", @outvals)."\n" ;
         } else {
            print STDERR "ERROR: PDB $pdb_id not found in $full_path\n" ;
         }

# If there are no PISA entries for this pdb id, go to next
         if (!exists $pdb2pisa->{$pdb_id}) {next;}

# Iterate through the PISA entries for this pdb id.
         foreach my $k ( 0 .. $#{$pdb2pisa->{$pdb_id}}) {

# Set the full PISA path and BDPFILES output values.
            my $pisa_path = pibase::PISA::get_pisa_filepath({
               pibase_specs => $pibase_specs,
               pisa_id => $pdb2pisa->{$pdb_id}->[$k]
            }) ;
            my $filebase = basename($pisa_path) ;
            my @outvals = ('0', $pisa_path, $filebase, $pdb_id, '0');

            if (-s $pisa_path) {
               print {$fh->{out}} join("\t", @outvals)."\n" ;
            } else {
               if (!-s $pisa_path && -e $pisa_path) {
                  print STDERR "ERROR: PISA $pdb2pisa->{$pdb_id}->[$k] is empty ".
                              "at $pisa_path\n" ;
               } elsif  (!-e $pisa_path) {
                  print STDERR "ERROR: PISA $pdb2pisa->{$pdb_id}->[$k] not found ".
                              "in $pisa_path\n" ;
               }
            }
         }
      }
      $counter++ ;
   }
   print STDERR "\n" ;
   close($fh->{out}) ;

   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
      $import_status = pibase::mysqlimport({fn => $fn->{out},
                           pibase_specs => $pibase_specs}) ;
   }

   return $import_status ;

}


=head2 create_subset_pdbs()

   Title:       create_subset_pdbs()
   OLD NAME:    cutdom_calc()
   Function:    extracts list of domains from their source PDB files,
                  calling subset_extract()
   Args:        none
   Returns:     nothing
   STDIN:       tabbed: bdp_file_path, domain_id (subset_id), domain definition
                  domain defintiion: tabbed: chain_1, start_resno_1, end_resno_1,
                     ..., chain_n, start_resno_n, end_resno_n
   STDOUT:      tabbed: subset_id, subset pdb file path
   Files out:   foreach domain:
                  o <directory per sid_2_domdir() logic>/<subset_id>.pdb

=cut

sub create_subset_pdbs {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }

   my $sid2bdp_path = {} ;
   my $sid2def = {};
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         chomp $line ;
         my ($bdp_path, $sid, @defs) = split(/\t/, $line) ;
         $sid2bdp_path->{$sid} = $bdp_path ;

         my $cur_fn = $bdp_path;
         my ($chain, $start, $end) ;
         my $j = 0 ;
         while ($j <= $#defs) {
            push @{$chain}, $defs[$j] ; $j++ ;

#            my $t_range = $defs[$j] ;
            push @{$start}, $defs[$j] ; $j++ ;

#            $t_range .= "-$defs[$j]" ;
            push @{$end}, $defs[$j] ; $j++ ;
#            print STDERR "sid $sid ".$chain->[$#{$chain}]." ".
#                         $start->[$#{$start}]." - ".$end->[$#{$end}]."\n" ;
         }
         $sid2def->{$sid} = {
            chain => $chain,
            start => $start,
            end => $end
         } ;
      }
      close(INF) ;

   } else {

      my ($dbh, $pibase) = pibase::connect_pibase() ;

      my $sid2bdp= pibase::mysql_hashload( $dbh, "SELECT subset_id, bdp_id ".
                                 "FROM subsets WHERE bdp_id IS NOT NULL") ;

      my $bdp2path= pibase::mysql_hashload( $dbh, "SELECT bdp_id, file_path ".
                                 "FROM bdp_files WHERE bdp_id IS NOT NULL") ;

      my ($subset_id, $chain_id, $start_resno, $end_resno) =
         pibase::mysql_fetchcols( $dbh, "SELECT subset_id, chain_id, ".
            "start_resno, end_resno FROM subsets_details ".
            "WHERE subset_id LIKE \"BDP\%\"");

      foreach my $j ( 0 .. $#{$subset_id}) {
         my $bdp_id = $sid2bdp->{$subset_id->[$j]} ;
         $sid2bdp_path->{$subset_id->[$j]} = $bdp2path->{$bdp_id} ;

         push @{$sid2def->{$subset_id->[$j]}->{chain}}, $chain_id->[$j] ;
         push @{$sid2def->{$subset_id->[$j]}->{start}}, $start_resno->[$j] ;
         push @{$sid2def->{$subset_id->[$j]}->{end}}, $end_resno->[$j] ;
      }
   }


   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      print "* create_subset_pdbs() ".localtime() if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{create_subset_pdbs_in}, $temp_fn->{create_subset_pdbs_in}) =
         tempfile("splits_create_subset_pdbs_input.XXXXX");
      ($temp_fh->{create_subset_pdbs_out}, $temp_fn->{create_subset_pdbs_out}) =
         tempfile("splits_create_subset_pdbs_SGEout_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{create_subset_pdbs_out}) ;
      ($temp_fh->{create_subset_pdbs_err}, $temp_fn->{create_subset_pdbs_err}) =
         tempfile("splits_create_subset_pdbs_SGEerr_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{create_subset_pdbs_err}) ;

      foreach my $sid (sort keys %{$sid2def}) {
         my @outvals = ($sid2bdp_path->{$sid}, $sid) ;
         foreach my $j ( 0 .. $#{$sid2def->{$sid}->{chain}}) {
            push @outvals, $sid2def->{$sid}->{chain}->[$j] ;
            push @outvals, $sid2def->{$sid}->{start}->[$j] ;
            push @outvals, $sid2def->{$sid}->{end}->[$j] ;
         }
         print {$temp_fh->{create_subset_pdbs_in}} join("\t",@outvals)."\n" ;
      }

      close($temp_fh->{create_subset_pdbs_in}) ;

      my $split_dir = tempdir("splits_create_subset_pdbs.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{create_subset_pdbs_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.create_subset_pdbs.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::data::access qw/create_subset_pdbs/ ;

main() ;

sub main {

         pibase::data::access::create_subset_pdbs({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.create_subset_pdbs.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.create_subset_pdbs.XXXXX");

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
         out_fn => $temp_fn->{create_subset_pdbs_out},
         err_fn => $temp_fn->{create_subset_pdbs_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{create_subset_pdbs_out},
           $temp_fn->{create_subset_pdbs_out}) ;
      open(REALOUTF,">".$pibase_specs->{buildfiles}->{subsets_files}) ;
      while (my $line = readline($temp_fh->{create_subset_pdbs_out})) {
         if ($line =~ /^\#/) {next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{create_subset_pdbs_out}) ;
      close(REALOUTF) ;

   } else {

      my $binaries = pibase::locate_binaries() ;

      my $movethreshold = 500 ;
      my $tempdir = tempdir(CLEANUP => 1) ;
      chdir $tempdir ;
      my (@movethese, @movedest) ;

      foreach my $sid (sort keys %{$sid2def}) {

         my $bdp_path = $sid2bdp_path->{$sid} ;

         my $chain = $sid2def->{$sid}->{chain} ;
         my $start= $sid2def->{$sid}->{start} ;
         my $end= $sid2def->{$sid}->{end} ;
         my $cur_fn = $bdp_path;
   
         my ($compress_fl, $unc_fh, $unc_fn) ;

         if ($bdp_path =~ /\.gz$/) {
            $compress_fl = 1;
            ($unc_fh, $unc_fn) = tempfile() ; close($unc_fh) ;
            system("$binaries->{zcat} $cur_fn > $unc_fn") ;
            $cur_fn = $unc_fn ;
         }
   
   
         my $altloc_fl = `$binaries->{altloc_check} < $cur_fn` ;
         chomp $altloc_fl ;
   
         if ($altloc_fl) {
            my ($t_fh, $t_fn) =
               tempfile("tempsource.XXXXXX", SUFFIX=>"pdb") ; close($t_fh) ;
   
            print STDERR "note: $bdp_path has an altloc set, now filtering: ." ;
            system("$binaries->{altloc_filter} $cur_fn > $t_fn") ;
            print STDERR "\bx\n" ;
            if ($compress_fl) {unlink $unc_fn;}
            $cur_fn = $t_fn ;
         }
   
#changed post-build 080825_0747         my $cutpdb_fn = "$sid.pdb" ;
         my $cutpdb_fn = "$sid.pdb.gz" ;

         my $extract_errors = pibase::PDB::subsets::subset_extract({
            in_fn => $cur_fn,
            out_fn => $cutpdb_fn,
            chain => $chain,
            start => $start,
            end => $end
         }) ;
   
         if ($#{$extract_errors} >= 0 ) {
            foreach my $j ( 0 .. $#{$extract_errors}) {
               print STDERR "ERROR: $bdp_path, extract $sid: ".
               "pibase::subsets::subset_extract(): $extract_errors->[$j]\n" ;
            }
         } elsif (-z $cutpdb_fn) {
            print STDERR "ERROR: $bdp_path, extract $sid: ".
                         "empty subset extract pdb file\n" ;
            unlink $cutpdb_fn ;
         } else {
            my $deposit_dir = pibase::sid_2_domdir($sid) ;
            push @movethese, $cutpdb_fn ;
            push @movedest, $deposit_dir;
            print "$sid\t$deposit_dir/$cutpdb_fn\n" ;
         }
   
         if ( ($altloc_fl) && ($bdp_path ne $cur_fn) ) {
            unlink $cur_fn;
         }

         if ($compress_fl) { unlink $unc_fn;}
   
         if ( $#movethese == $movethreshold ) {
            foreach my $j ( 0 .. $#movethese) {
               if (! -d $movedest[$j]) {
                  mkpath($movedest[$j]) ;}
               pibase::safe_move($movethese[$j], $movedest[$j]) ;}
            @movethese = () ;
            @movedest = () ;
         }
      }

      foreach my $j ( 0 .. $#movethese) {
         if (! -d $movedest[$j]) {
            mkpath($movedest[$j]) ;}
         pibase::safe_move($movethese[$j], $movedest[$j]) ;}
   }

   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_residues_tables file (specified in specs) into pibase
      $import_status = pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{subsets_files}}) ;
   }

   return $import_status ;

}


=head2 cutdom_prep()

   Title:       cutdom_prep()
   Function:    prepare input for cutdom_calc() extraction of domains from
                  their source PDB files
   Args:        none
   Returns:     nothing
   Tables:      pibase.subsets
                pibase.bdp_files
                pibase.subsets_details
   STDIN:       subset_id
   STDOUT:      tabbed: subset_id, subset pdb file path
   Files out:   tabbed: source bdp file path, subset_id, domain definition
                  domain defintiion: tabbed: chain_1, start_resno_1, end_resno_1,
                     ..., chain_n, start_resno_n, end_resno_n

=cut

sub cutdom_prep {

   require DBI ;

   my $usage = __FILE__." < alilist" ;
   if ($#ARGV > 0 ) {
      die "usage: $usage\n" ; }

   my ($dbh, $pibase) = pibase::connect_pibase() ;

   my $sid2bdp=
      pibase::mysql_hashload( $dbh, "SELECT subset_id, bdp_id ".
                                 "FROM subsets WHERE bdp_id IS NOT NULL") ;

   my $bdp2path=
      pibase::mysql_hashload( $dbh, "SELECT bdp_id, file_path ".
                                 "FROM bdp_files WHERE bdp_id IS NOT NULL") ;

   my ($subset_id, $chain_id, $start_resno, $end_resno) =
      pibase::mysql_fetchcols( $dbh, "SELECT subset_id, chain_id, ".
                               "start_resno, end_resno FROM subsets_details" );

   my $subdef ;
   foreach my $j ( 0 .. $#{$subset_id}) {
      push @{$subdef->{$subset_id->[$j]}},
           ($chain_id->[$j], $start_resno->[$j],$end_resno->[$j]) ;
   }

   while (my $sid = <STDIN>) {
      chomp $sid;
      my @outx = ($bdp2path->{$sid2bdp->{$sid}}, $sid, @{$subdef->{$sid}}) ;
      print join("\t", @outx)."\n" ;
   }

}


=head2 sqlcreate_parser()

   Title:       sqlcreate_parser()
   Function:    Creates perl module raw_table_specs.pm with pibase table
                  structure descriptions from the PIBASE SQL DDL
   Args:        none
   Returns:     nothing
   STDIN:       PIBASE SQL DDL (from dia2sql output)
   STDOUT:      raw_table_specs.pm perl module

=cut

sub sqlcreate_parser {

   my $tables ;
   my $pibase = 'pibase';

   my $curtable = 0 ;
   while (my $line = <STDIN>) {
      if ($line !~ /^#/ && $line !~ /^$/) {
         chomp $line;
         if ($line =~ /^CREATE TABLE/) {
            ($curtable) = ($line =~ /^CREATE TABLE (.+) \(/) ;
            ($tables->{$curtable}->{name}) = $curtable ;
         } elsif ($line =~ /\);/) {
            $curtable = '' ;
         } else {
            if ($line =~ /^\s*PRIMARY KEY/) {
               $line =~ s/\,$// ;
               ($tables->{$curtable}->{prikey}) = ($line =~ /^\s*PRIMARY KEY (.+)$/);
            } else {
               $line =~ s/\,$// ;
               my ($fieldname, $spec) = ($line =~ /(\S+)\s+(.+)$/) ;
               push @{$tables->{$curtable}->{field_name}}, $fieldname ;
               push @{$tables->{$curtable}->{field_spec}}, $spec;
            }
         }
      }
   }

   print 'package '.$pibase.'::raw_table_specs ;'."\n" ;
   print 'require Exporter;'."\n" ;
   print "\n" ;
   print 'use strict ;'."\n" ;
   print 'use warnings ;'."\n" ;
   print 'my @ISA = qw/Exporter/ ;'."\n" ;
   print 'my @EXPORT = qw/full_table_specs/ ;'."\n" ;
   print "\n" ;
   print "\n" ;

   print "sub full_table_specs {"."\n\n" ;
   print '   my $tables'." ;\n" ;
   foreach my $j ( keys %{$tables}) {
      if (exists $tables->{$j}->{prikey}) {
         print '   $tables->{'.$j.'}->{prikey} = "'.
            $tables->{$j}->{prikey}."\" ;\n" ;}

      foreach my $k ( 0 .. $#{$tables->{$j}->{field_name}} ) {
         print '   $tables->{'.$j.'}->{field_name}->['.$k.'] = "'.
            $tables->{$j}->{field_name}->[$k]."\" ;\n" ;
         print '   $tables->{'.$j.'}->{field_spec}->['.$k.'] = "'.
            $tables->{$j}->{field_spec}->[$k]."\" ;\n" ;
      }
   }
   print "\n" ;
   print '   return $tables ;'."\n" ;
   print "}\n" ;
   print "\n" ;
   print "1 ;\n" ;

}


=head2 copy_table_to_disk()

   Title:       copy_table_to_disk()
   Function:    Copies PIBASE table to disk
   Args:        $_->{pibase_specs} - optional; will use pibase::get_specs() if not
                $_->{tables} = [tablename1, tablename2, ...]
   Returns:     1

=cut

sub copy_table_to_disk {

   my $in = shift ;

   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs} ;
   }

   foreach my $table (@{$in->{tables}}) {
      if (!-e $pibase_specs->{tod_dir}) {
         mkpath $pibase_specs->{tod_dir} ; }
      my $out_fn = $pibase_specs->{tod_dir}."/$table" ;
      my $tcom = "echo 'SELECT * FROM ".$table.";' | mysql ".
         "-u ".$pibase_specs->{user}." -p'".$pibase_specs->{pass}.
         "' ".$pibase_specs->{db}." | sed 1d > $out_fn" ;
      system($tcom) ;
   }

   return 1 ;

}

=head2 copy_table_to_disk_custom()

   Title:       copy_table_to_disk_custom()
   Function:    Copies PIBASE table to disk
   Args:        $_->{pibase_specs} - optional; will use pibase::get_specs() if not
                $_->{tables} = [tablename1, tablename2, ...]
                $_->{destination} = destination directory for table on disk
   Returns:     1

=cut

sub copy_table_to_disk_custom {

   my $in = shift ;
   if (!exists $in->{destination}) {
      die "destination not specified";}

   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs} ;
   }


   foreach my $table (@{$in->{tables}}) {
      if (exists $in->{verbose}) {
         print STDERR "Dumping $table: " ; }

      if (!-e $in->{destination}) {
         mkpath $in->{destination} ; }
      my $out_fn = $in->{destination}."/$table" ;
      my $tcom = "echo 'SELECT * FROM ".$table.";' | mysql ".
         "-u ".$pibase_specs->{user}." -p'".$pibase_specs->{pass}.
         "' ".$pibase_specs->{db}." | sed 1d > $out_fn" ;
      system($tcom) ;
      if (exists $in->{verbose}) {
         print STDERR "X\n "; }
   }

   return 1 ;

}
