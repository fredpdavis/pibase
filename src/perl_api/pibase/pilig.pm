=head1 NAME

pibase::pilig - perl module for pibase-ligbase overlap calculations

=head1 DESCRIPTION

Perl module with routines to cross-query pibase and ligbase
to get small molecule - protein interaction site overlap
statistics

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

=head2 logic

0a. prep 
1a. calculate domain SASA;
1b. assign_exp
2. assign_pi() from pibase
3. assign_lig() from ligbase data (can we get the tables locally)
4. collate_perfam
5. collate_perinstance

=cut

package pibase::pilig ;
use strict;
use warnings;
#use Devel::Size qw/total_size/ ;


use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/pilig_main run_pilig assign_lig assign_pi calc_sid_sasa assign_exp assign_pepnuci_realres assign_pepnuci/ ;

use pibase ;
use pibase::ASTRAL ;
use pibase::SGE ;
use Cwd qw/getcwd/ ;
use File::Temp qw/tempfile tempdir/ ;
use File::Path qw/mkpath/ ;

use File::Basename qw/basename/ ;
use Sys::Hostname qw/hostname/ ;
use File::Temp qw/tempfile tempdir/ ;
use pibase::interatomic_contacts qw/contacts_select contacts_select_inter special_contact raw_contacts_select/;
use pibase::modeller qw/calc_sasa/ ;


sub set_pilig_specs {

   my $params ;
   $params->{DEBUG} = 0;
   $params->{DEBUGALN} = 0;


   $params->{PARAM_MIN_LIGMW} = 250;
   $params->{PARAM_MAX_LIGMW} = 1000;
   $params->{CALC_CONSSCORE_FLAG} = 1;
   $params->{CALC_PVAL_FLAG} = 1;

#LIGBASE is built at 5A, use same cutoff for PIBASE interface residues
   $params->{PARAM_INTERFACE_DIST_THRESH} = 5.0 ;

   $params->{PARAM_MIN_PEPTIDE_LENGTH} = 5 ;
   $params->{PARAM_MIN_NUMRES_INTERACTING_WITH_PEPTIDE} = 5 ;

   $params->{PARAM_SCPERCACC_THRESH} = 7; #more than this is an exposed residue
   $params->{PARAM_MIN_NUMCONTACTS} = 500 ; #change to minimum SASA cutoff

   $params->{PARAM_MIN_BS_SIMILARITY} = 0.9 ; #threshold to merge binding sites.


   $params->{pilig_rootdir} = "/groups/eddy/home/davisf/work/pilig" ;
   $params->{ligbase}->{activeres} = $params->{pilig_rootdir}.
      '/data/ligbase_new__activeres.gz' ;

   $params->{msdchem_xmldir} = "/groups/eddy/home/davisf/work/".
                               "databases/msdchem/xml" ;
   $params->{superfam_fxn} = "/groups/eddy/home/davisf/work/pilig/data/".
                             "scop.annotation.1.73.txt" ;
   $params->{superfam_fxn_categories} = "/groups/eddy/home/davisf/work/".
                                        "pilig/data/scop.larger.categories" ;

   $params->{exp_dir} = $params->{pilig_rootdir}."/run/calc_sid_sasa" ;
   $params->{assign_pi_dir} = $params->{pilig_rootdir}."/run/assign_pi" ;
   $params->{assign_lig_dir} = $params->{pilig_rootdir}."/run/assign_lig" ;
   $params->{assign_exp_dir} = $params->{pilig_rootdir}."/run/assign_exp" ;

   $params->{outfiles}->{liginfo} = $params->{pilig_rootdir}."/run/liginfo.out";
   $params->{outfiles}->{calc_sid_sasa} = $params->{pilig_rootdir}.
                                       '/run/calc_sid_sasa.out' ;
   $params->{outfiles}->{assign_exp} = $params->{pilig_rootdir}.
                                       '/run/assign_exp.out' ;
   $params->{outfiles}->{assign_pi} = $params->{pilig_rootdir}.
                                       '/run/assign_pi.out' ;
   $params->{outfiles}->{assign_lig} = $params->{pilig_rootdir}.
                                       '/run/assign_lig.out' ;
   $params->{outfiles}->{assign_pepnuci_realres} = $params->{pilig_rootdir}.
                                       '/run/assign_pepnuci_realres.out' ;
   $params->{outfiles}->{assign_pepnuci} = $params->{pilig_rootdir}.
                                       '/run/assign_pepnuci.out' ;

   $params->{outfiles}->{assign_pi_clusters} = $params->{pilig_rootdir}.
                                       '/run/assign_pi_clusters.out' ;
   $params->{outfiles}->{assign_lig_clusters} = $params->{pilig_rootdir}.
                                       '/run/assign_lig_clusters.out' ;
   $params->{outfiles}->{assign_pepnuci_clusters} = $params->{pilig_rootdir}.
                                       '/run/assign_pepnuci_clusters.out' ;

   return $params ;

}


sub pilig_main {

   my $usage = __FILE__ ;
   my $mode_usage ;
   $mode_usage->{'prep'}->{'numjobs'}++ ;
   $mode_usage->{'prep'}->{'priority'}++ ;

   $mode_usage->{'prep_sasacalc'} = {};
   $mode_usage->{'sasacalc'} = {};
   $mode_usage->{'assign_exp'} = {};

   $mode_usage->{'assign_lig'} = {};
   $mode_usage->{'assign_pi'} = {} ;

   $mode_usage->{'collate_perfam'} = {} ;
   $mode_usage->{'collate_instance'} = {} ;

   if ($#ARGV < 0) {
      die "$usage mode (".join(',', sort keys %{$mode_usage}).")\n" ; }

   my $mode = shift @ARGV ; $mode =~ s/^-// ;

   if (!exists $mode_usage->{$mode}) {
      die "$usage mode (".join(',', sort keys %{$mode_usage}).")\n" ; }


   my $options ;
   my $j = 0 ;
   while ($j <= $#ARGV) {
      my $param = $ARGV[$j] ;
      $param =~ s/^-// ;
      if (!exists $mode_usage->{$mode}->{$param}) {
         die "$param not recognized, possible options: ".
            join(", ", sort keys %{$mode_usage->{$mode}})."\n" ;
      }
      my $val = $ARGV[($j + 1)] ;
      $options->{$param} = $val ;
      $j+=2 ;
   }
   
   my $mode2sub= {
      prep => \&run_prep ,
      prep_sasacalc => \&run_prep_sasacalc ,
      sasacalc => \&run_sasacalc ,
      assign_exp => \&run_assign_exp ,
      assign_lig => \&run_assign_lig ,
      assign_pi => \&run_assign_pi,
      collate_perfam => \&run_collate_perfam,
      collate_instance_all => \&run_collate_perinstance_all,
      collate_instance => \&run_collate_perinstance,
   } ;

   $mode2sub->{$mode}->($options) ;

}


sub run_pilig {

#Not initially setting up directory structure- do it online prn
#Step -1. build mysql table structures.

   my $in = shift ;

# Step 0. fetch external data
   my $build_status ;
#   my $fetch_external_data = fetch_external_data() ;
#   my $pibase_specs = $fetch_external_data->{pibase_specs} ;
#   $build_status->{fetch_external_data} = $fetch_external_data->{status} ;

   my $pibase_specs = pibase::get_specs() ;

   my $cluster_fl = 0;
   if (exists $in->{steps}->{cluster_fl} &&
       $in->{steps}->{cluster_fl} == 1) {
      $cluster_fl = 1 ; }

   my $num_specified_steps = keys %{$in->{steps}} ;

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{assign_pi} &&
      $in->{steps}->{assign_pi} == 1)) {
      $build_status->{assign_pi} = pibase::pilig::assign_pi({
         cluster_fl => $cluster_fl,
         pibase_specs => $pibase_specs
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{cluster_pi} &&
      $in->{steps}->{cluster_pi} == 1)) {
      $build_status->{cluster_pi} = pibase::pilig::run_cluster_pi({
         cluster_fl => $cluster_fl,
         pibase_specs => $pibase_specs
      }) ;
   }




   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{assign_pepnuci_realres} &&
      $in->{steps}->{assign_pepnuci_realres} == 1)) {
      $build_status->{assign_pepnuci_realres} =
         pibase::pilig::assign_pepnuci_realres({
         cluster_fl => $cluster_fl,
         pibase_specs => $pibase_specs
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{assign_pepnuci} &&
      $in->{steps}->{assign_pepnuci} == 1)) {
      $build_status->{assign_pepnuci} = pibase::pilig::assign_pepnuci({
         cluster_fl => $cluster_fl,
         pibase_specs => $pibase_specs
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{cluster_pep} &&
      $in->{steps}->{cluster_pep} == 1)) {
      $build_status->{cluster_pep} = pibase::pilig::run_cluster_pep({
         cluster_fl => $cluster_fl,
         pibase_specs => $pibase_specs
      }) ;
   }





   if (exists $in->{steps} &&
       exists $in->{steps}->{extract_msdchem_liginfo} &&
       $in->{steps}->{extract_msdchem_liginfo} == 1) {
      $build_status->{extract_msdchem_liginfo} =
         pibase::pilig::extract_msdchem_liginfo({
            pibase_specs => $pibase_specs }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
       exists $in->{steps}->{assign_lig} &&
       $in->{steps}->{assign_lig} == 1)) {
      $build_status->{assign_lig} = pibase::pilig::assign_lig({
         cluster_fl => $cluster_fl,
         pibase_specs => $pibase_specs
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{cluster_lig} &&
      $in->{steps}->{cluster_lig} == 1)) {
      $build_status->{cluster_lig} = pibase::pilig::run_cluster_lig() ;
   }



   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{calc_sid_sasa} &&
      $in->{steps}->{calc_sid_sasa} == 1)) {
      $build_status->{calc_sid_sasa} = pibase::pilig::calc_sid_sasa({
         numtasks_cluster => 20,
         cluster_fl => 1,
         pibase_specs => $pibase_specs
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{assign_exp} &&
      $in->{steps}->{assign_exp} == 1)) {
      $build_status->{assign_exp} = pibase::pilig::assign_exp({
         cluster_fl => $cluster_fl,
         pibase_specs => $pibase_specs
      }) ;
   }



   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{collate_perfam} &&
      $in->{steps}->{collate_perfam} == 1)) {
      $build_status->{collate_perfam} = pibase::pilig::collate_perfam({
         pibase_specs => $pibase_specs
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{collate_perinstance_PINT} &&
      $in->{steps}->{collate_perinstance_PINT} == 1)) {
      $build_status->{collate_perinstance_PINT}=
         pibase::pilig::collate_perinstance({
         mode => "PINT",
         pibase_specs => $pibase_specs
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{collate_perinstance_pINT} &&
      $in->{steps}->{collate_perinstance_pINT} == 1)) {
      $build_status->{collate_perinstance_pINT} =
         pibase::pilig::collate_perinstance({
         mode => "pINT",
         pibase_specs => $pibase_specs
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{make_paper_figures} &&
      $in->{steps}->{make_paper_figures} == 1)) {
      $build_status->{make_paper_figures}= pibase::pilig::make_paper_figures({
         pibase_specs => $pibase_specs
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{create_mysql_tables} &&
      $in->{steps}->{create_mysql_tables} == 1)) {
      $build_status->{create_mysql_tables}= pibase::pilig::create_mysql_tables({
         pibase_specs => $pibase_specs
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{create_resno_tables} &&
      $in->{steps}->{create_resno_tables} == 1)) {
      $build_status->{create_resno_tables} =
         pibase::pilig::create_resno_tables({
         pibase_specs => $pibase_specs,
         bstype => $in->{steps}->{bstype}
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{create_overlapaln_tables} &&
      $in->{steps}->{create_overlapaln_tables} == 1)) {
      $build_status->{create_overlapaln_tables} =
         pibase::pilig::create_overlapaln_tables({
         pibase_specs => $pibase_specs,
         bstype => $in->{steps}->{bstype}
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{calc_pilig_db_stats} &&
      $in->{steps}->{calc_pilig_db_stats} == 1)) {
      $build_status->{calc_pilig_db_stats} =
         pibase::pilig::calc_pilig_db_stats({
         pibase_specs => $pibase_specs,
      }) ;
   }

}


sub _pilig_run_prep_sasacalc {

   my $dbh = connect_pibase_ligbase() ;
   my ($sid, $bdp, $class) = pibase::mysql_fetchcols($dbh->{pi},
     "SELECT subset_id, a.bdp_id, class FROM subsets as a, bdp_files as b ".
     "WHERE subset_source_id = 1 and a.bdp_id = b.bdp_id and b.raw_pdb = 1");

   my $class2sid ; my $sid2bdp ;
   foreach my $j ( 0 .. $#{$sid}) {
      push @{$class2sid->{$class->[$j]}}, $sid->[$j] ;
      $sid2bdp->{$sid->[$j]} = $bdp->[$j] ;
   }

   foreach my $tclass (sort keys %{$class2sid}) {
      if ($tclass !~ /^[a-g]/) {next;}
      foreach my $tsid (@{$class2sid->{$tclass}}) {
         print join("\t",$tsid, $sid2bdp->{$tsid}, $tclass)."\n" ;
      }
   }

}

sub _pilig__getinput_sasacalc {

   my ($sidlist, $sid2bdp, $sid2class) ;
   while (my $line = <STDIN>) {
      chomp $line;
      if ($line =~ /^#/) {next;}
      my ($sid, $bdp, $class) = split(/\t/, $line) ;
      push @{$sidlist}, $sid ;
      $sid2bdp->{$sid} = $bdp ;
      $sid2class->{$sid} = $class ;
   }

   return {
      sids => $sidlist,
      sid2bdp => $sid2bdp,
      sid2class => $sid2class,
   } ;

}


sub _pilig_run_sasacalc {

   my $in = _getinput_sasacalc() ;
   my $pilig_specs = set_pilig_specs() ;

   my $pibase_specs = pibase::get_specs() ;
   my $modeller_bin = $pibase_specs->{binaries}->{modeller} ;

   foreach my $sid (@{$in->{sids}}) {
      print STDERR "NOW ON: $sid\n" ;
      my $exposed ;

      my $bdpdir = pibase::sid_2_domdir($sid) ;
      my $sid_fn = $bdpdir."/$sid.pdb" ;

      my $subset_sasa = pibase::modeller::calc_sasa(
            {pdb_fn => $sid_fn,
             surftyp => 2,
             modeller_bin => $modeller_bin}) ;

      if ($#{$subset_sasa->{error_fl}} >= 0 ) {
         foreach my $j ( 0 .. $#{$subset_sasa->{error_fl}}) {
            print STDERR "ERROR: subset $sid: ".
                         "calc_sasa() $subset_sasa->{error_fl}->[$j]\n";
            }
         next;
      }

      foreach my $j ( 0 .. $#{$subset_sasa->{res_sasa}->{resno}}) {
         my $resno = $subset_sasa->{res_sasa}->{resno}->[$j] ;
         my $chain = $subset_sasa->{res_sasa}->{chain}->[$j] ;
         my $ressig = $resno."_".$chain ;

         if ($subset_sasa->{res_sasa}->{sc_perc}->[$j] >
             $pilig_specs->{PARAM_SCPERCACC_THRESH}) {
            push @{$exposed}, $ressig ; }
      }

      print join("\t", $sid, $in->{sid2bdp}->{$sid},
                  $in->{sid2class}->{$sid}, @{$exposed})."\n" ;

   }

}


sub _pilig_run_assign_exp {

   my $in = _getinput_assign_exp() ;
   my $fn = set_locations() ;
   my $astral = astral_preload({ fn => $fn }) ;
   my $pb = _pilig_tod_pibase_preload() ;

   foreach my $classtype (qw/fam sf/) {
      foreach my $class (sort keys %{$astral->{classes}->{$classtype}}) {

         if (!exists $in->{class2sid}->{$classtype}->{$class}) { next;}

         print STDERR "NOW ON $classtype $class (line ".__LINE__.")\n" ;
         my $pos2bs ;
         my $class_aln = pibase::ASTRAL::load_asteroids_aln({
            aln_fn => $fn->{aln}->{$classtype}.'/'.$class.'.aln.fa',
            seq_fn => $fn->{aln}->{$classtype}.'/'.$class.'.fa',
            raf => $astral->{raf},
            gdseqh => $astral->{gdseqh},
            seqclcont100 => $astral->{seqcl2cont}->{100},
            seqcl100 => $astral->{seqcl}->{100},
            allchains => $pb->{pdbchains}
         }) ;

         foreach my $sid (sort 
                    keys %{$in->{class2sid}->{$classtype}->{$class}} ){

            my $osid = $pb->{sid2osid}->{$sid} ;
            my $pdb = $pb->{sid2pdb}->{$sid} ;

            my @talnres = () ;
            my @undefres = () ;
            foreach my $res (@{$in->{sid2expres}->{$sid}}) {
               my $alnres = $class_aln->{resno2pos}->{$osid}->{$res} ;
               if (!defined $alnres) {
#                  print STDERR "WARNING: $osid ($sid) residue $res not ".
#                               "found in ASTRAL alignment\n" ;
                  push @undefres, 'undef';
               } else {
                  push @talnres, $alnres ;
               }
#               $pos2bs->{$alnres}->{e}->{$sid}++ ;
            }

            my @salnres = ();
            push @salnres, sort {$a <=> $b} @talnres ;
            push @salnres, @undefres ;
            my $alnposstring = join(',', @salnres) ;

            my @outvals = ( $pdb, $sid, $osid, $classtype,
                                 $class, "E", $alnposstring,
                                 $class_aln->{alnlength} ) ;
            print join("\t", @outvals)."\n" ;

         }
      }
   }

}


sub collate_perinstance {
# fpd090505_1043 - modified routine used to generate ms numbers and figures so
# that the subset_id of the ligand binding site is specified for each overlap

   require Bit::Vector ;

   my $standardres = _list_standardres() ;

# Per-PPI
# 0. load all of ligands first - lig info
# - iterate over assignments, and have individual
#   bit vectors set for the ligand binding sites in their
#   respective domain families
#
# per-PPI: iterate over each domain--domain interface
#  1. list any ligands actually bound to this interface
#  2. iterate over all ligands in that family/superfamily, find
#     - quantify per ligand coverage
#     - quantify cumulative ligand coverage
#
# per-ligand:


   my $in = shift ;
   my $mode = "PINT" ;
   if (exists $in->{mode}) { $mode = $in->{mode} ; }

   my $pibase_specs = pibase::get_specs() ;
   my $pilig_specs = set_pilig_specs() ;
   my $fn = $pilig_specs->{fn} ;

   my $astral = _pilig_astral_preload({
      pibase_specs => $pibase_specs,
      pilig_specs => $pilig_specs,
   }) ;
   my $pb = _pilig_tod_pibase_preload() ;
   my $liginfo = _pilig_load_liginfo() ;

   my $class2alnlength = {};

   my $expbits = readin_expassignments({
      fn => $pilig_specs->{outfiles}->{assign_exp},
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;

   my $pepnucibits = readin_pepnuciassignments({
      fn => $pilig_specs->{outfiles}->{assign_pepnuci_clusters},
      clustrep_fl => 1,
      pb => $pb,
      pilig_specs => $pilig_specs,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;
   my $ligbits = readin_ligassignments({
      fn => $pilig_specs->{outfiles}->{assign_lig_clusters},
      clustrep_fl => 1,
      liginfo => $liginfo,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;

   my $pibits_both = readin_piassignments({
      fn => $pilig_specs->{outfiles}->{assign_pi_clusters},
      clustrep_fl => 1,
      pb => $pb,
      dont_read_jdomains_fl => 1,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;
   my $pibits = $pibits_both->{pibits} ;
   my $interfaces = $pibits_both->{interfaces} ;

   my $alltogether = combine_piligpepexpbits({
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



# herenow090303_1128 - rejiggered so that ASTRAL alignment
# is read in - then ask $aln->{aln}->{domain_id} for alignment seq
# to comput pairwise seq identity for intersecting bits
# between ligand and protein binding sites.
# have to do this before the per-interaction, since need to load
# the full alignment on a per-family basis..

# 1. PINT-LINT - iterate over all domain interfaces involving this class
#   foreach my $classtype (qw/fam/)
   my ($hist2d_binsize, $hist2d_maxbin);
   $hist2d_binsize->{ligbs_seqid} = 5 ; #last bin includes right-end
   $hist2d_maxbin->{ligbs_seqid} = POSIX::floor(99 /
                                          $hist2d_binsize->{ligbs_seqid}) ;


   if ($mode eq 'PINT') {
      my $classtype = 'fam' ;

      my @headers_sum = ("PINT_LSUM", "SID1", "SID2", "CHAINS",
                          "CLASS1", "CLASS2") ;

      foreach my $s (1, 2) {
         push @headers_sum, "numres_p_$s", "cumulative_numres_l_and_p_$s",
                            "max_l_and_p_$s", "lig_max_l_and_p_$s";
         foreach my $j ( 0 .. $hist2d_maxbin->{ligbs_seqid}) {
            my $t_seqid_cutoff = $j * $hist2d_binsize->{ligbs_seqid} ;
            push @headers_sum, "cum_l_and_p_".$s."_perseqid_".$t_seqid_cutoff;
            push @headers_sum, "max_l_and_p_".$s."_perseqid_".$t_seqid_cutoff;
         }
      }
      print '#'.join("\t",@headers_sum)."\n" ;

      @headers_sum = ("SID1", "SID2", "CHAINS",
                          "CLASS1", "CLASS2",
                          "side", "SID", "LIG",
                           "LIG_SID", #fpd090505_1317
                          "numres_p_side", "numres_p_IDENT",
                          "numres_l", "numres_l_IDENT",
                          "numres_l_and_p_side", "numres_l_and_p_side_IDENT",
                          "numres_domain_nongap",
                          "numres_domain_nongap_IDENT",
                         ) ;
      print '#'.join("\t",@headers_sum)."\n" ;

# create a family indexed list of binding sites
      my $class2sid12_side ;
      foreach my $sid12 (keys %{$interfaces}) {
         my $sid ;
         ($sid->{1}, $sid->{2}) = split(/\t/, $sid12) ;
         foreach my $side (1,2) {
            $class2sid12_side->{$pb->{sid2class}->{$classtype}->{$sid->{$side}}}->{$sid12."\t".$side}++ ; }
      }

      my $skip_sid12_side ;
      my ($pistats, $ligstats) ; 
#      my $CUTOUT = 5 ; print STDERR "WARNING: CUTS OUT AFTER $CUTOUT FAMS\n";
#      my $t_class_count = 0 ;
      foreach my $class (keys %{$class2sid12_side}) {
#         if ($t_class_count > $CUTOUT) {last;}
         
# load the ASTRAL alignment - so have sequence string for each dom,
# compute pairwise seqid for each bs-ligand match

         my $curclasstype = 'fam' ; my $curfam = $class ;

         if ($curfam !~ /[a-g]/) {next;} #ASTRAL alignments only available
                                         #for class a-g
         my $class_aln = pibase::ASTRAL::load_asteroids_aln({
                  aln_fn => $pibase_specs->{asteroids}->{$curclasstype.'_aln'}.
                  '/'.$curfam.'.fasta_aln' ,
                  seq_fn => $pibase_specs->{asteroids}->{$curclasstype.'_seq'}.
                  '/'.$curfam.'.fa' ,
                  raf => $astral->{raf},
                  gdseqh => $astral->{gdseqh},
                  seqclcont100 => $astral->{seqcl2cont}->{100},
                  seqcl100 => $astral->{seqcl}->{100},
                  allchains => $pb->{pdbchains}
         }) ;
#         print STDERR "ASTRAL alignment for $curfam contains: ".join(", ", 
#            keys %{$class_aln->{aln}})."\n" ;

         foreach my $sid12_side (keys %{$class2sid12_side->{$class}}) {
            my ($sid, $side) ;
            ($sid->{1}, $sid->{2}, $side) = split(/\t/, $sid12_side);
            my $sid12 = $sid->{1}."\t".$sid->{2} ;
            my $classes ;
            $classes->{1} = $pb->{sid2class}->{'fam'}->{$sid->{1}} ;
            $classes->{2} = $pb->{sid2class}->{'fam'}->{$sid->{2}} ;

            my ($sid_origdom_p) = ($sid->{$side} =~ /SCOP\.(.+)/) ;
#            print STDERR "$sid->{$side} origdomP: $sid_origdom_p\n" ;
            my $pdb = $interfaces->{$sid12}->{pdb} ;
            my $chains = $interfaces->{$sid12}->{chains} ;
      
            if (!exists $class2alnlength->{$classtype}->{$class}) {
               next ; }
            my $alnlength = $class2alnlength->{$classtype}->{$class};
            my $bs_bits;
            my $curligs ;
      
            if (!exists $interfaces->{$sid12}->{$side} ||
               !exists $interfaces->{$sid12}->{$side}->{pibits}->{$classtype} ||
               ($interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Norm()
                 == 0)) {
      
               $skip_sid12_side->{$sid12_side}++ ;
               print STDERR "WARNING: skipping $sid->{1} -- $sid->{2} ".
                  " undefined binding site alignment positions for ".
                  " $sid->{$side}\n";
               next;
            }
      
            $bs_bits = $interfaces->{$sid12}->{$side}->{pibits}->{$classtype} ;
            $pistats->{$sid12_side}->{p} = $bs_bits->Norm() ;

            if (exists $class2allligbits->{$classtype}->{$class}) {
               $curligs= $class2allligbits->{$classtype}->{$class};
            } else {
               next;
            }
      
# init vals
# fpd090305_1356 - HERENOW stratify max and cum by seqid cutoffs

            foreach my $type (qw/p cumlig_l cumlig_l_and_p max_l_and_p/) {
               $ligstats->{$sid12_side}->{$type} = 0 ; }
            foreach my $type (qw/lig_max_l_and_p/) {
               $ligstats->{$sid12_side}->{$type} = 'U'; }

            foreach my $type (qw/cumlig_l_and_p_perseqid max_l_and_p_perseqid/){
               my $seqid_cutoff = 0 ;
               while ($seqid_cutoff < 100) {
                  $ligstats->{$sid12_side}->{$type."_perseqid_$seqid_cutoff"}=0;
                  $seqid_cutoff += $hist2d_binsize->{ligbs_seqid} ;
               }
            }
      
            $ligstats->{$sid12_side}->{"max_l_and_p"} = 0 ;
            $ligstats->{$sid12_side}->{"lig_max_l_and_p"} = 'undef' ;
            $ligstats->{$sid12_side}->{cumlig} = Bit::Vector->new($alnlength) ;

            my $seqid_cutoff = 0;
            while ($seqid_cutoff < 100) {
               $ligstats->{$sid12_side}->{"cumlig_perseqid_$seqid_cutoff"} =
                  Bit::Vector->new($alnlength) ;
               $seqid_cutoff += $hist2d_binsize->{ligbs_seqid} ;
            }
      
            my $temp_bits_and = Bit::Vector->new($alnlength) ;
            my $temp_bits_or = Bit::Vector->new($alnlength) ;
      
            foreach my $j (0 .. $#{$curligs}) {
# get this domain's SCOP id
               my ($sid_origdom_l) =
                  ($curligs->[$j]->[2] =~ /SCOP\.(.+)/) ;
#               print STDERR $curligs->[$j]->[2]." origdomL: $sid_origdom_l\n" ;
   
               my $curligsig = $curligs->[$j]->[0] ;
               my $outcurligsig = $curligsig; $outcurligsig =~ s/\t/:/g ;
               my $curligbits = $curligs->[$j]->[1] ;
               my $curligsid = $curligs->[$j]->[2] ; #fpd090505_1049 
                  
               $ligstats->{$sid12_side}->{cumlig}->Or(
                  $ligstats->{$sid12_side}->{cumlig}, $curligbits) ;
      
#intersect the ligands bit vector and the bit vector for the binding
# site positions of this particular interface;
               $temp_bits_and->And($curligbits, $bs_bits) ;
               $temp_bits_or->Or($curligbits, $bs_bits) ;
   
               my $curs ;
               $curs->{l_and_p} = $temp_bits_and->Norm() ;
               $curs->{l_or_p} = $temp_bits_or->Norm() ;
               $curs->{l} = $curligbits->Norm() ;
               $curs->{p} = $bs_bits->Norm() ;
   
      
               if ($curs->{l_and_p} == 0) {next;}
   
# compute pairwise sequence identity over intersected residues.
# is the ASTRAL aln loaded 0- or 1-indexed?
               $curs->{l_and_p_IDENT} = 0 ;
               $curs->{l_IDENT} = 0 ;
               $curs->{p_IDENT} = 0 ;

               foreach my $aln_pos ($temp_bits_and->Index_List_Read()) {
                  if (substr($class_aln->{aln}->{$sid_origdom_l}, $aln_pos,1) eq
                      substr($class_aln->{aln}->{$sid_origdom_p}, $aln_pos,1)) {
                      $curs->{l_and_p_IDENT}++ ; } }

               foreach my $aln_pos ($curligbits->Index_List_Read()) {
                  if (substr($class_aln->{aln}->{$sid_origdom_l}, $aln_pos,1) eq
                      substr($class_aln->{aln}->{$sid_origdom_p}, $aln_pos,1)) {
                      $curs->{l_IDENT}++ ; } }

               foreach my $aln_pos ($bs_bits->Index_List_Read()) {
                  if (substr($class_aln->{aln}->{$sid_origdom_l}, $aln_pos,1) eq
                      substr($class_aln->{aln}->{$sid_origdom_p}, $aln_pos,1)) {
                      $curs->{p_IDENT}++ ; } }

               $curs->{wholedom_numnongap} = 0 ;
               $curs->{wholedom_IDENT} = 0 ;
               foreach my $aln_pos (0 .. $alnlength) {
                  my $a= substr($class_aln->{aln}->{$sid_origdom_l},$aln_pos,1);
                  my $b= substr($class_aln->{aln}->{$sid_origdom_p},$aln_pos,1);
                  if ($a ne '-' || $b ne '-') {
                     $curs->{wholedom_numnongap}++ ;
                     if ($a eq $b) { $curs->{wholedom_IDENT}++ ; }
                  }
               }

               my $curoverlap ;
               $curoverlap->{obi} = $curs->{l_and_p} / $curs->{l_or_p} ;
               $curoverlap->{opi} = $curs->{l_and_p} / $curs->{p} ;
               $curoverlap->{olig} = $curs->{l_and_p} / $curs->{l} ;
      
               if (!exists $ligstats->{$sid12_side}->{"max_l_and_p"} ||
                  $curs->{l_and_p} > $ligstats->{$sid12_side}->{"max_l_and_p"}){
                  $ligstats->{$sid12_side}->{"max_l_and_p"} = $curs->{l_and_p};
                  $ligstats->{$sid12_side}->{"lig_max_l_and_p"} = $outcurligsig;
               }

               my $cur_ligbs_seqid = ($curs->{l_IDENT} / $curs->{l}) * 100 ;
               my $seqid_cutoff = 0 ;
               my $seqid_bin_ceil = POSIX::floor($cur_ligbs_seqid /
                                        $hist2d_binsize->{ligbs_seqid});

               if ($seqid_bin_ceil > $hist2d_maxbin->{ligbs_seqid}) {
                   $seqid_bin_ceil = $hist2d_maxbin->{ligbs_seqid}; }

               foreach my $t_seqid_bin (0 ..  $seqid_bin_ceil) {
                  my $t_seqid_cutoff = $t_seqid_bin *
                                       $hist2d_binsize->{ligbs_seqid} ;
                  if (
  !exists $ligstats->{$sid12_side}->{"max_l_and_p_perseqid_$t_seqid_cutoff"} ||
                     $curs->{l_and_p} >
   $ligstats->{$sid12_side}->{"max_l_and_p_perseqid_$t_seqid_cutoff"}
                     ){
   $ligstats->{$sid12_side}->{"max_l_and_p_perseqid_$t_seqid_cutoff"} =
      $curs->{l_and_p};
                  }

   $ligstats->{$sid12_side}->{"cumlig_perseqid_$t_seqid_cutoff"}->Or(
      $ligstats->{$sid12_side}->{"cumlig_perseqid_$t_seqid_cutoff"},
      $curligbits) ;
               }
   
#HERENOW fpd090303_2042 - seqid computed fine - make sure that the
#  per-side / summary separation is proper; compare line by line
# to old: sub BK_PRE090303_1855_collate_perinstance
   
               my @outvals = ($sid->{1}, $sid->{2}, $chains,
                                 $classes->{1}, $classes->{2},
                                 $side, $sid->{$side}, $outcurligsig,
                                 $curligsid, #fpd090505_1051 
                                 $curs->{p}, $curs->{p_IDENT},
                                 $curs->{l}, $curs->{l_IDENT},
                                 $curs->{l_and_p}, $curs->{l_and_p_IDENT},
                                 $curs->{wholedom_numnongap},
                                 $curs->{wholedom_IDENT},
                                ) ;
               print join("\t", @outvals)."\n";
            }
      
            $ligstats->{$sid12_side}->{"cumlig_l"} =
               $ligstats->{$sid12_side}->{cumlig}->Norm() ;
            $temp_bits_and->And($ligstats->{$sid12_side}->{cumlig}, $bs_bits) ;
            $ligstats->{$sid12_side}->{"cumlig_l_and_p"} =
               $temp_bits_and->Norm() ;

            foreach my $t_seqid_bin (0 ..  $hist2d_maxbin->{ligbs_seqid}) {
               my $t_seqid_cutoff = $t_seqid_bin *
                                       $hist2d_binsize->{ligbs_seqid} ;
               $temp_bits_and->And(
                  $ligstats->{$sid12_side}->{"cumlig_perseqid_$t_seqid_cutoff"},
                  $bs_bits) ;
               $ligstats->{$sid12_side}->{"cumlig_l_and_p_perseqid_$t_seqid_cutoff"} = $temp_bits_and->Norm() ;
            }
         }
#         $t_class_count++ ;
      }


# print PINT_LSUM interface coverage summaries
      foreach my $sid12 (keys %{$interfaces}) {
         my $sid ; ($sid->{1}, $sid->{2}) = split(/\t/, $sid12) ;
         my $class ;
         $class->{1} = $pb->{sid2class}->{$classtype}->{$sid->{1}} ;
         $class->{2} = $pb->{sid2class}->{$classtype}->{$sid->{2}} ;
         my $chains = $interfaces->{$sid12}->{chains} ;
         my @outvals =  ("PINT_LSUM", $sid->{1}, $sid->{2}, $chains,
                         $class->{1}, $class->{2}) ;
   
         foreach my $s (1, 2) {
            my $sid12_side = $sid12."\t".$s ;
            if (exists $skip_sid12_side->{$sid12_side} ||
               !exists $ligstats->{$sid12_side} ||
               $ligstats->{$sid12_side}->{"cumlig_l"} == 0 ) {
   
               if (exists $skip_sid12_side->{$sid12_side} ||
                   $class->{$s} !~ /[a-g]/) {
                  push @outvals, 0, 0 ;
               } else {
                  push @outvals, $pistats->{$sid12_side}->{p}, 0 ;
               }
   
               push @outvals, 0, "U" ;
               push @outvals,
                  split('','0'x(($hist2d_maxbin->{ligbs_seqid} + 1 )* 2));
               next;
            }
   
            push @outvals, $pistats->{$sid12_side}->{p} ;
            push @outvals, $ligstats->{$sid12_side}->{cumlig_l_and_p} ;
            push @outvals, $ligstats->{$sid12_side}->{"max_l_and_p"} ;
            push @outvals, $ligstats->{$sid12_side}->{"lig_max_l_and_p"} ;

# print per seqid cutoff coverage numbers.
            foreach my $t_seqid_bin (0 .. $hist2d_maxbin->{ligbs_seqid}) {
               my $t_seqid_cutoff = $t_seqid_bin *
                                    $hist2d_binsize->{ligbs_seqid} ;
               if (!exists $ligstats->{$sid12_side}->{"cumlig_l_and_p_perseqid_".$t_seqid_cutoff}) {$ligstats->{$sid12_side}->{"cumlig_l_and_p_perseqid_".$t_seqid_cutoff} = 0; }
               if (!exists $ligstats->{$sid12_side}->{"max_l_and_p_perseqid_".$t_seqid_cutoff}) {$ligstats->{$sid12_side}->{"max_l_and_p_perseqid_".$t_seqid_cutoff} = 0; }

               push @outvals, $ligstats->{$sid12_side}->{"cumlig_l_and_p_perseqid_".$t_seqid_cutoff} ;
               push @outvals, $ligstats->{$sid12_side}->{"max_l_and_p_perseqid_".$t_seqid_cutoff} ;
            }

         }
         print join("\t", @outvals)."\n"; 
      }
   }

# 2. pINT-LINT - iterate over all peptide involving this class
# new081231_1812 
   if ($mode eq 'pINT') {
      my $classtype = 'fam' ;
      my @headers_sum = ("pINT_LSUM", "SID", "CHAIN", "CHAIN_LENGTH", "CLASS",
                          "numres_p", "cumulative_numres_l_and_p",
                          "max_l_and_p", "lig_max_l_and_p",
                         ) ;
      foreach my $j ( 0 .. $hist2d_maxbin->{ligbs_seqid}) {
         my $t_seqid_cutoff = $j * $hist2d_binsize->{ligbs_seqid} ;
         push @headers_sum, "cum_l_and_p_perseqid_".$t_seqid_cutoff;
         push @headers_sum, "max_l_and_p_perseqid_".$t_seqid_cutoff;
      }
      print '#'.join("\t",@headers_sum)."\n" ;

      @headers_sum = ("SID", "CHAIN", "CHAIN_LENGTH", "CLASS", "LIG",
                      "LIG_SID", #fpd090505_1317 
                      "numres_p", "numres_p_IDENT",
                      "numres_l", "numres_l_IDENT",
                      "numres_l_and_p", "numres_l_and_p_side_IDENT",
                      "numres_domain_nongap", "numres_domain_nongap_IDENT",
                     ) ;
      print '#'.join("\t",@headers_sum)."\n" ;

# group them by class so that can calc seqids - rejigger like did to PINT-LINT
      my $class2pint;
      foreach my $sid (keys %{$pepnucibits->{pep}->{$classtype}}) {
         $class2pint->{$pb->{sid2class}->{$classtype}->{$sid}}->{$sid}++ ; }

      foreach my $class (keys %{$class2pint}) {
         my $curclasstype = 'fam' ; my $curfam = $class ;

         if ($curfam !~ /[a-g]/) {next;} #ASTRAL alignments only available
                                         #for class a-g

# load ASTRAL alignment for this family
         my $class_aln = pibase::ASTRAL::load_asteroids_aln({
                  aln_fn => $pibase_specs->{asteroids}->{$curclasstype.'_aln'}.
                  '/'.$curfam.'.fasta_aln' ,
                  seq_fn => $pibase_specs->{asteroids}->{$curclasstype.'_seq'}.
                  '/'.$curfam.'.fa' ,
                  raf => $astral->{raf},
                  gdseqh => $astral->{gdseqh},
                  seqclcont100 => $astral->{seqcl2cont}->{100},
                  seqcl100 => $astral->{seqcl}->{100},
                  allchains => $pb->{pdbchains}
         }) ;

         foreach my $sid (keys %{$class2pint->{$class}}){
          my ($sid_origdom_p) = ($sid =~ /SCOP\.(.+)/) ;

          foreach my $targetch (
             keys %{$pepnucibits->{pep}->{$classtype}->{$sid}} ){

             if ($targetch eq 'cumulative') {next;}

         my $pdb = $pb->{sid2pdb}->{$sid} ;
         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;
         my $alnlength = $class2alnlength->{$classtype}->{$class};
   
         my $bs_bits = $pepnucibits->{pep}->{$classtype}->{$sid}->{$targetch};
         my $curligs ;
   
         my $pistats ;
         my $ligstats ; 
         if (exists $class2allligbits->{$classtype}->{$class}) {
            $curligs = $class2allligbits->{$classtype}->{$class} ; }
   
# init vals
         foreach my $type (qw/p cumlig_l cumlig_l_and_p max_l_and_p/) {
            $ligstats->{$type} = 0 ; }
         foreach my $type (qw/lig_max_l_and_p/) {
            $ligstats->{$type} = 'U'; }

         foreach my $type (qw/p cumlig_l cumlig_l_and_p max_l_and_p/) {
            $ligstats->{$type} = 0 ; }
         foreach my $type (qw/lig_max_l_and_p/) {
            $ligstats->{$type} = 'U'; }

         foreach my $type (qw/cumlig_l_and_p_perseqid max_l_and_p_perseqid/) {
            my $seqid_cutoff = 0 ;
            while ($seqid_cutoff < 100) {
               $ligstats->{$type."_perseqid_$seqid_cutoff"}=0;
               $seqid_cutoff += $hist2d_binsize->{ligbs_seqid} ;
            }
         }
   
         $pistats->{p} = $bs_bits->Norm() ;
         $ligstats->{cumlig} = Bit::Vector->new($alnlength) ;
         $ligstats->{"max_l_and_p"} = 0 ;
         $ligstats->{"lig_max_l_and_p"} = 'undef' ;

         my $seqid_cutoff = 0;
         while ($seqid_cutoff < 100) {
            $ligstats->{"cumlig_perseqid_$seqid_cutoff"} =
               Bit::Vector->new($alnlength) ;
            $seqid_cutoff += $hist2d_binsize->{ligbs_seqid} ;
         }

   
            my $temp_bits_and = Bit::Vector->new($alnlength) ;
            my $temp_bits_or = Bit::Vector->new($alnlength) ;
   
            foreach my $j (0 .. $#{$curligs}) {
               my $curligsig = $curligs->[$j]->[0] ;
               my $outcurligsig = $curligsig; $outcurligsig =~ s/\t/:/g ;
               my $curligbits = $curligs->[$j]->[1] ;
               my $curligsid = $curligs->[$j]->[2] ; #fpd090505_1049 

               my ($sid_origdom_l) =
                  ($curligs->[$j]->[2] =~ /SCOP\.(.+)/) ;
               
               $ligstats->{cumlig}->Or($ligstats->{cumlig}, $curligbits) ;
   
#intersect the ligands bit vector and the bit vector for the bindin
# site positions of this particular interface;

               $temp_bits_and->And($curligbits, $bs_bits) ;
               $temp_bits_or->Or($curligbits, $bs_bits) ;
   
               my $curs ;
               $curs->{l_and_p} = $temp_bits_and->Norm() ;
               $curs->{l_or_p} = $temp_bits_or->Norm() ;
               $curs->{l} = $curligbits->Norm() ;
               $curs->{p} = $bs_bits->Norm() ;
   
               if ($curs->{l_and_p} == 0) {next;}

# compute pairwise sequence identity over intersected residues.
# is the ASTRAL aln loaded 0- or 1-indexed?
               $curs->{l_and_p_IDENT} = 0 ;
               $curs->{l_IDENT} = 0 ;
               $curs->{p_IDENT} = 0 ;

               foreach my $aln_pos ($temp_bits_and->Index_List_Read()) {
                  if (substr($class_aln->{aln}->{$sid_origdom_l}, $aln_pos,1) eq
                      substr($class_aln->{aln}->{$sid_origdom_p}, $aln_pos,1)) {
                      $curs->{l_and_p_IDENT}++ ; } }

               foreach my $aln_pos ($curligbits->Index_List_Read()) {
                  if (substr($class_aln->{aln}->{$sid_origdom_l}, $aln_pos,1) eq
                      substr($class_aln->{aln}->{$sid_origdom_p}, $aln_pos,1)) {
                      $curs->{l_IDENT}++ ; } }

               foreach my $aln_pos ($bs_bits->Index_List_Read()) {
                  if (substr($class_aln->{aln}->{$sid_origdom_l}, $aln_pos,1) eq
                      substr($class_aln->{aln}->{$sid_origdom_p}, $aln_pos,1)) {
                      $curs->{p_IDENT}++ ; } }

               $curs->{wholedom_numnongap} = 0 ;
               $curs->{wholedom_IDENT} = 0 ;
               foreach my $aln_pos (0 .. $alnlength) {
                  my $a= substr($class_aln->{aln}->{$sid_origdom_l},$aln_pos,1);
                  my $b= substr($class_aln->{aln}->{$sid_origdom_p},$aln_pos,1);
                  if ($a ne '-' || $b ne '-') {
                     $curs->{wholedom_numnongap}++ ;
                     if ($a eq $b) { $curs->{wholedom_IDENT}++ ; }
                  }
               }
   
               if (!exists $ligstats->{"max_l_and_p"} ||
                   $curs->{l_and_p} > $ligstats->{"max_l_and_p"}) {
                  $ligstats->{"max_l_and_p"} = $curs->{l_and_p} ;
                  $ligstats->{"lig_max_l_and_p"} = $outcurligsig ;
               }

               my $cur_ligbs_seqid = ($curs->{l_IDENT} / $curs->{l}) * 100 ;
               my $seqid_cutoff = 0 ;
               my $seqid_bin_ceil = POSIX::floor($cur_ligbs_seqid /
                                        $hist2d_binsize->{ligbs_seqid});

               if ($seqid_bin_ceil > $hist2d_maxbin->{ligbs_seqid}) {
                   $seqid_bin_ceil = $hist2d_maxbin->{ligbs_seqid}; }

               foreach my $t_seqid_bin (0 ..  $seqid_bin_ceil) {
                  my $t_seqid_cutoff = $t_seqid_bin *
                                       $hist2d_binsize->{ligbs_seqid} ;
                 if (!exists $ligstats->{"max_l_and_p_perseqid_$t_seqid_cutoff"}
                 || $curs->{l_and_p} >
                     $ligstats->{"max_l_and_p_perseqid_$t_seqid_cutoff"}){
                  $ligstats->{"max_l_and_p_perseqid_$t_seqid_cutoff"} =
                     $curs->{l_and_p};
                  }

                  $ligstats->{"cumlig_perseqid_$t_seqid_cutoff"}->Or(
                     $ligstats->{"cumlig_perseqid_$t_seqid_cutoff"},
                     $curligbits) ;
               }
   
               my @outvals = ($sid, $targetch,
                       $pepnucibits->{chain_info}->{$targetch}->{chain_length},
                       $class,
                       $outcurligsig,
                       $curligsid, #fpd090505_1058 
                       $curs->{p}, $curs->{p_IDENT},
                       $curs->{l}, $curs->{l_IDENT},
                       $curs->{l_and_p}, $curs->{l_and_p_IDENT},
                       $curs->{wholedom_numnongap}, $curs->{wholedom_IDENT},
               ) ;

               print join("\t", @outvals)."\n";
            }
   
            $ligstats->{"cumlig_l"} = $ligstats->{cumlig}->Norm() ;
            $temp_bits_and->And($ligstats->{cumlig}, $bs_bits) ;
            $ligstats->{"cumlig_l_and_p"} = $temp_bits_and->Norm() ;

            foreach my $t_seqid_bin (0 ..  $hist2d_maxbin->{ligbs_seqid}) {
               my $t_seqid_cutoff = $t_seqid_bin *
                                    $hist2d_binsize->{ligbs_seqid};
               $temp_bits_and->And(
                  $ligstats->{"cumlig_perseqid_$t_seqid_cutoff"},
                  $bs_bits) ;
               $ligstats->{"cumlig_l_and_p_perseqid_$t_seqid_cutoff"} =
                  $temp_bits_and->Norm() ;
            }
   
            my @outvals =  ("pINT_LSUM", $sid, $targetch,
                  $pepnucibits->{chain_info}->{$targetch}->{chain_length},
                  $class) ;
   
            if ($ligstats->{"cumlig_l"} == 0 ) {
               push @outvals, $pistats->{p}, 0 ;
               push @outvals, 0, "U" ;
            } else {
               push @outvals, $pistats->{p} ;
               push @outvals, $ligstats->{cumlig_l_and_p} ;

               push @outvals, $ligstats->{"max_l_and_p"} ;
               push @outvals, $ligstats->{"lig_max_l_and_p"} ;
            }

            foreach my $t_seqid_bin (0 .. $hist2d_maxbin->{ligbs_seqid}) {
               my $t_seqid_cutoff = $t_seqid_bin * $hist2d_binsize->{ligbs_seqid} ;
               if (!exists $ligstats->{"cumlig_l_and_p_perseqid_".$t_seqid_cutoff}) {
                  $ligstats->{"cumlig_l_and_p_perseqid_".$t_seqid_cutoff} = 0; }
               if (!exists $ligstats->{"max_l_and_p_perseqid_".$t_seqid_cutoff}) {
                  $ligstats->{"max_l_and_p_perseqid_".$t_seqid_cutoff} = 0; }

               push @outvals,
                  $ligstats->{"cumlig_l_and_p_perseqid_".$t_seqid_cutoff},
                  $ligstats->{"max_l_and_p_perseqid_".$t_seqid_cutoff};
            }

            print join("\t", @outvals)."\n"; 

          }
         }
      }


   }

}


sub PAPERNUMBERS_20090505_collate_perinstance {
# fpd090505_1043 - this is the routine that generate numbers for the manuscript
# it was modified, above, so that the subset_id of the ligand binding site
# is specified per each actual overlap

   require Bit::Vector ;

   my $standardres = _list_standardres() ;

# Per-PPI
# 0. load all of ligands first - lig info
# - iterate over assignments, and have individual
#   bit vectors set for the ligand binding sites in their
#   respective domain families
#
# per-PPI: iterate over each domain--domain interface
#  1. list any ligands actually bound to this interface
#  2. iterate over all ligands in that family/superfamily, find
#     - quantify per ligand coverage
#     - quantify cumulative ligand coverage
#
# per-ligand:


   my $in = shift ;
   my $mode = "PINT" ;
   if (exists $in->{mode}) { $mode = $in->{mode} ; }

   my $pibase_specs = pibase::get_specs() ;
   my $pilig_specs = set_pilig_specs() ;
   my $fn = $pilig_specs->{fn} ;

   my $astral = _pilig_astral_preload({
      pibase_specs => $pibase_specs,
      pilig_specs => $pilig_specs,
   }) ;
   my $pb = _pilig_tod_pibase_preload() ;
   my $liginfo = _pilig_load_liginfo() ;

   my $class2alnlength = {};

   my $expbits = readin_expassignments({
      fn => $pilig_specs->{outfiles}->{assign_exp},
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;

   my $pepnucibits = readin_pepnuciassignments({
      fn => $pilig_specs->{outfiles}->{assign_pepnuci_clusters},
      clustrep_fl => 1,
      pb => $pb,
      pilig_specs => $pilig_specs,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;
   my $ligbits = readin_ligassignments({
      fn => $pilig_specs->{outfiles}->{assign_lig_clusters},
      clustrep_fl => 1,
      liginfo => $liginfo,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;

   my $pibits_both = readin_piassignments({
      fn => $pilig_specs->{outfiles}->{assign_pi_clusters},
      clustrep_fl => 1,
      pb => $pb,
      dont_read_jdomains_fl => 1,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;
   my $pibits = $pibits_both->{pibits} ;
   my $interfaces = $pibits_both->{interfaces} ;

   my $alltogether = combine_piligpepexpbits({
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



# herenow090303_1128 - rejiggered so that ASTRAL alignment
# is read in - then ask $aln->{aln}->{domain_id} for alignment seq
# to comput pairwise seq identity for intersecting bits
# between ligand and protein binding sites.
# have to do this before the per-interaction, since need to load
# the full alignment on a per-family basis..

# 1. PINT-LINT - iterate over all domain interfaces involving this class
#   foreach my $classtype (qw/fam/)
   my ($hist2d_binsize, $hist2d_maxbin);
   $hist2d_binsize->{ligbs_seqid} = 5 ; #last bin includes right-end
   $hist2d_maxbin->{ligbs_seqid} = POSIX::floor(99 /
                                          $hist2d_binsize->{ligbs_seqid}) ;


   if ($mode eq 'PINT') {
      my $classtype = 'fam' ;

      my @headers_sum = ("PINT_LSUM", "SID1", "SID2", "CHAINS",
                          "CLASS1", "CLASS2") ;

      foreach my $s (1, 2) {
         push @headers_sum, "numres_p_$s", "cumulative_numres_l_and_p_$s",
                            "max_l_and_p_$s", "lig_max_l_and_p_$s";
         foreach my $j ( 0 .. $hist2d_maxbin->{ligbs_seqid}) {
            my $t_seqid_cutoff = $j * $hist2d_binsize->{ligbs_seqid} ;
            push @headers_sum, "cum_l_and_p_".$s."_perseqid_".$t_seqid_cutoff;
            push @headers_sum, "max_l_and_p_".$s."_perseqid_".$t_seqid_cutoff;
         }
      }
      print '#'.join("\t",@headers_sum)."\n" ;

      @headers_sum = ("SID1", "SID2", "CHAINS",
                          "CLASS1", "CLASS2",
                          "side", "SID", "LIG",
                          "numres_p_side", "numres_p_IDENT",
                          "numres_l", "numres_l_IDENT",
                          "numres_l_and_p_side", "numres_l_and_p_side_IDENT",
                          "numres_domain_nongap",
                          "numres_domain_nongap_IDENT",
                         ) ;
      print '#'.join("\t",@headers_sum)."\n" ;

# create a family indexed list of binding sites
      my $class2sid12_side ;
      foreach my $sid12 (keys %{$interfaces}) {
         my $sid ;
         ($sid->{1}, $sid->{2}) = split(/\t/, $sid12) ;
         foreach my $side (1,2) {
            $class2sid12_side->{$pb->{sid2class}->{$classtype}->{$sid->{$side}}}->{$sid12."\t".$side}++ ; }
      }

      my $skip_sid12_side ;
      my ($pistats, $ligstats) ; 
#      my $CUTOUT = 5 ; print STDERR "WARNING: CUTS OUT AFTER $CUTOUT FAMS\n";
#      my $t_class_count = 0 ;
      foreach my $class (keys %{$class2sid12_side}) {
#         if ($t_class_count > $CUTOUT) {last;}
         
# load the ASTRAL alignment - so have sequence string for each dom,
# compute pairwise seqid for each bs-ligand match

         my $curclasstype = 'fam' ; my $curfam = $class ;

         if ($curfam !~ /[a-g]/) {next;} #ASTRAL alignments only available
                                         #for class a-g
         my $class_aln = pibase::ASTRAL::load_asteroids_aln({
                  aln_fn => $pibase_specs->{asteroids}->{$curclasstype.'_aln'}.
                  '/'.$curfam.'.fasta_aln' ,
                  seq_fn => $pibase_specs->{asteroids}->{$curclasstype.'_seq'}.
                  '/'.$curfam.'.fa' ,
                  raf => $astral->{raf},
                  gdseqh => $astral->{gdseqh},
                  seqclcont100 => $astral->{seqcl2cont}->{100},
                  seqcl100 => $astral->{seqcl}->{100},
                  allchains => $pb->{pdbchains}
         }) ;
#         print STDERR "ASTRAL alignment for $curfam contains: ".join(", ", 
#            keys %{$class_aln->{aln}})."\n" ;

         foreach my $sid12_side (keys %{$class2sid12_side->{$class}}) {
            my ($sid, $side) ;
            ($sid->{1}, $sid->{2}, $side) = split(/\t/, $sid12_side);
            my $sid12 = $sid->{1}."\t".$sid->{2} ;
            my $classes ;
            $classes->{1} = $pb->{sid2class}->{'fam'}->{$sid->{1}} ;
            $classes->{2} = $pb->{sid2class}->{'fam'}->{$sid->{2}} ;

            my ($sid_origdom_p) = ($sid->{$side} =~ /SCOP\.(.+)/) ;
#            print STDERR "$sid->{$side} origdomP: $sid_origdom_p\n" ;
            my $pdb = $interfaces->{$sid12}->{pdb} ;
            my $chains = $interfaces->{$sid12}->{chains} ;
      
            if (!exists $class2alnlength->{$classtype}->{$class}) {
               next ; }
            my $alnlength = $class2alnlength->{$classtype}->{$class};
            my $bs_bits;
            my $curligs ;
      
            if (!exists $interfaces->{$sid12}->{$side} ||
               !exists $interfaces->{$sid12}->{$side}->{pibits}->{$classtype} ||
               ($interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Norm()
                 == 0)) {
      
               $skip_sid12_side->{$sid12_side}++ ;
               print STDERR "WARNING: skipping $sid->{1} -- $sid->{2} ".
                  " undefined binding site alignment positions for ".
                  " $sid->{$side}\n";
               next;
            }
      
            $bs_bits = $interfaces->{$sid12}->{$side}->{pibits}->{$classtype} ;
            $pistats->{$sid12_side}->{p} = $bs_bits->Norm() ;

            if (exists $class2allligbits->{$classtype}->{$class}) {
               $curligs= $class2allligbits->{$classtype}->{$class};
            } else {
               next;
            }
      
# init vals
# fpd090305_1356 - HERENOW stratify max and cum by seqid cutoffs

            foreach my $type (qw/p cumlig_l cumlig_l_and_p max_l_and_p/) {
               $ligstats->{$sid12_side}->{$type} = 0 ; }
            foreach my $type (qw/lig_max_l_and_p/) {
               $ligstats->{$sid12_side}->{$type} = 'U'; }

            foreach my $type (qw/cumlig_l_and_p_perseqid max_l_and_p_perseqid/){
               my $seqid_cutoff = 0 ;
               while ($seqid_cutoff < 100) {
                  $ligstats->{$sid12_side}->{$type."_perseqid_$seqid_cutoff"}=0;
                  $seqid_cutoff += $hist2d_binsize->{ligbs_seqid} ;
               }
            }
      
            $ligstats->{$sid12_side}->{"max_l_and_p"} = 0 ;
            $ligstats->{$sid12_side}->{"lig_max_l_and_p"} = 'undef' ;
            $ligstats->{$sid12_side}->{cumlig} = Bit::Vector->new($alnlength) ;

            my $seqid_cutoff = 0;
            while ($seqid_cutoff < 100) {
               $ligstats->{$sid12_side}->{"cumlig_perseqid_$seqid_cutoff"} =
                  Bit::Vector->new($alnlength) ;
               $seqid_cutoff += $hist2d_binsize->{ligbs_seqid} ;
            }
      
            my $temp_bits_and = Bit::Vector->new($alnlength) ;
            my $temp_bits_or = Bit::Vector->new($alnlength) ;
      
            foreach my $j (0 .. $#{$curligs}) {
# get this domain's SCOP id
               my ($sid_origdom_l) =
                  ($curligs->[$j]->[2] =~ /SCOP\.(.+)/) ;
#               print STDERR $curligs->[$j]->[2]." origdomL: $sid_origdom_l\n" ;
   
               my $curligsig = $curligs->[$j]->[0] ;
               my $outcurligsig = $curligsig; $outcurligsig =~ s/\t/:/g ;
               my $curligbits = $curligs->[$j]->[1] ;
                  
               $ligstats->{$sid12_side}->{cumlig}->Or(
                  $ligstats->{$sid12_side}->{cumlig}, $curligbits) ;
      
#intersect the ligands bit vector and the bit vector for the binding
# site positions of this particular interface;
               $temp_bits_and->And($curligbits, $bs_bits) ;
               $temp_bits_or->Or($curligbits, $bs_bits) ;
   
               my $curs ;
               $curs->{l_and_p} = $temp_bits_and->Norm() ;
               $curs->{l_or_p} = $temp_bits_or->Norm() ;
               $curs->{l} = $curligbits->Norm() ;
               $curs->{p} = $bs_bits->Norm() ;
   
      
               if ($curs->{l_and_p} == 0) {next;}
   
# compute pairwise sequence identity over intersected residues.
# is the ASTRAL aln loaded 0- or 1-indexed?
               $curs->{l_and_p_IDENT} = 0 ;
               $curs->{l_IDENT} = 0 ;
               $curs->{p_IDENT} = 0 ;

               foreach my $aln_pos ($temp_bits_and->Index_List_Read()) {
                  if (substr($class_aln->{aln}->{$sid_origdom_l}, $aln_pos,1) eq
                      substr($class_aln->{aln}->{$sid_origdom_p}, $aln_pos,1)) {
                      $curs->{l_and_p_IDENT}++ ; } }

               foreach my $aln_pos ($curligbits->Index_List_Read()) {
                  if (substr($class_aln->{aln}->{$sid_origdom_l}, $aln_pos,1) eq
                      substr($class_aln->{aln}->{$sid_origdom_p}, $aln_pos,1)) {
                      $curs->{l_IDENT}++ ; } }

               foreach my $aln_pos ($bs_bits->Index_List_Read()) {
                  if (substr($class_aln->{aln}->{$sid_origdom_l}, $aln_pos,1) eq
                      substr($class_aln->{aln}->{$sid_origdom_p}, $aln_pos,1)) {
                      $curs->{p_IDENT}++ ; } }

               $curs->{wholedom_numnongap} = 0 ;
               $curs->{wholedom_IDENT} = 0 ;
               foreach my $aln_pos (0 .. $alnlength) {
                  my $a= substr($class_aln->{aln}->{$sid_origdom_l},$aln_pos,1);
                  my $b= substr($class_aln->{aln}->{$sid_origdom_p},$aln_pos,1);
                  if ($a ne '-' || $b ne '-') {
                     $curs->{wholedom_numnongap}++ ;
                     if ($a eq $b) { $curs->{wholedom_IDENT}++ ; }
                  }
               }

               my $curoverlap ;
               $curoverlap->{obi} = $curs->{l_and_p} / $curs->{l_or_p} ;
               $curoverlap->{opi} = $curs->{l_and_p} / $curs->{p} ;
               $curoverlap->{olig} = $curs->{l_and_p} / $curs->{l} ;
      
               if (!exists $ligstats->{$sid12_side}->{"max_l_and_p"} ||
                  $curs->{l_and_p} > $ligstats->{$sid12_side}->{"max_l_and_p"}){
                  $ligstats->{$sid12_side}->{"max_l_and_p"} = $curs->{l_and_p};
                  $ligstats->{$sid12_side}->{"lig_max_l_and_p"} = $outcurligsig;
               }

               my $cur_ligbs_seqid = ($curs->{l_IDENT} / $curs->{l}) * 100 ;
               my $seqid_cutoff = 0 ;
               my $seqid_bin_ceil = POSIX::floor($cur_ligbs_seqid /
                                        $hist2d_binsize->{ligbs_seqid});

               if ($seqid_bin_ceil > $hist2d_maxbin->{ligbs_seqid}) {
                   $seqid_bin_ceil = $hist2d_maxbin->{ligbs_seqid}; }

               foreach my $t_seqid_bin (0 ..  $seqid_bin_ceil) {
                  my $t_seqid_cutoff = $t_seqid_bin *
                                       $hist2d_binsize->{ligbs_seqid} ;
                  if (
  !exists $ligstats->{$sid12_side}->{"max_l_and_p_perseqid_$t_seqid_cutoff"} ||
                     $curs->{l_and_p} >
   $ligstats->{$sid12_side}->{"max_l_and_p_perseqid_$t_seqid_cutoff"}
                     ){
   $ligstats->{$sid12_side}->{"max_l_and_p_perseqid_$t_seqid_cutoff"} =
      $curs->{l_and_p};
                  }

   $ligstats->{$sid12_side}->{"cumlig_perseqid_$t_seqid_cutoff"}->Or(
      $ligstats->{$sid12_side}->{"cumlig_perseqid_$t_seqid_cutoff"},
      $curligbits) ;
               }
   
#HERENOW fpd090303_2042 - seqid computed fine - make sure that the
#  per-side / summary separation is proper; compare line by line
# to old: sub BK_PRE090303_1855_collate_perinstance
   
               my @outvals = ($sid->{1}, $sid->{2}, $chains,
                                 $classes->{1}, $classes->{2},
                                 $side, $sid->{$side}, $outcurligsig,
                                 $curs->{p}, $curs->{p_IDENT},
                                 $curs->{l}, $curs->{l_IDENT},
                                 $curs->{l_and_p}, $curs->{l_and_p_IDENT},
                                 $curs->{wholedom_numnongap},
                                 $curs->{wholedom_IDENT},
                                ) ;
               print join("\t", @outvals)."\n";
            }
      
            $ligstats->{$sid12_side}->{"cumlig_l"} =
               $ligstats->{$sid12_side}->{cumlig}->Norm() ;
            $temp_bits_and->And($ligstats->{$sid12_side}->{cumlig}, $bs_bits) ;
            $ligstats->{$sid12_side}->{"cumlig_l_and_p"} =
               $temp_bits_and->Norm() ;

            foreach my $t_seqid_bin (0 ..  $hist2d_maxbin->{ligbs_seqid}) {
               my $t_seqid_cutoff = $t_seqid_bin *
                                       $hist2d_binsize->{ligbs_seqid} ;
               $temp_bits_and->And(
                  $ligstats->{$sid12_side}->{"cumlig_perseqid_$t_seqid_cutoff"},
                  $bs_bits) ;
               $ligstats->{$sid12_side}->{"cumlig_l_and_p_perseqid_$t_seqid_cutoff"} = $temp_bits_and->Norm() ;
            }
         }
#         $t_class_count++ ;
      }


# print PINT_LSUM interface coverage summaries
      foreach my $sid12 (keys %{$interfaces}) {
         my $sid ; ($sid->{1}, $sid->{2}) = split(/\t/, $sid12) ;
         my $class ;
         $class->{1} = $pb->{sid2class}->{$classtype}->{$sid->{1}} ;
         $class->{2} = $pb->{sid2class}->{$classtype}->{$sid->{2}} ;
         my $chains = $interfaces->{$sid12}->{chains} ;
         my @outvals =  ("PINT_LSUM", $sid->{1}, $sid->{2}, $chains,
                         $class->{1}, $class->{2}) ;
   
         foreach my $s (1, 2) {
            my $sid12_side = $sid12."\t".$s ;
            if (exists $skip_sid12_side->{$sid12_side} ||
               !exists $ligstats->{$sid12_side} ||
               $ligstats->{$sid12_side}->{"cumlig_l"} == 0 ) {
   
               if (exists $skip_sid12_side->{$sid12_side} ||
                   $class->{$s} !~ /[a-g]/) {
                  push @outvals, 0, 0 ;
               } else {
                  push @outvals, $pistats->{$sid12_side}->{p}, 0 ;
               }
   
               push @outvals, 0, "U" ;
               push @outvals,
                  split('','0'x(($hist2d_maxbin->{ligbs_seqid} + 1 )* 2));
               next;
            }
   
            push @outvals, $pistats->{$sid12_side}->{p} ;
            push @outvals, $ligstats->{$sid12_side}->{cumlig_l_and_p} ;
            push @outvals, $ligstats->{$sid12_side}->{"max_l_and_p"} ;
            push @outvals, $ligstats->{$sid12_side}->{"lig_max_l_and_p"} ;

# print per seqid cutoff coverage numbers.
            foreach my $t_seqid_bin (0 .. $hist2d_maxbin->{ligbs_seqid}) {
               my $t_seqid_cutoff = $t_seqid_bin *
                                    $hist2d_binsize->{ligbs_seqid} ;
               if (!exists $ligstats->{$sid12_side}->{"cumlig_l_and_p_perseqid_".$t_seqid_cutoff}) {$ligstats->{$sid12_side}->{"cumlig_l_and_p_perseqid_".$t_seqid_cutoff} = 0; }
               if (!exists $ligstats->{$sid12_side}->{"max_l_and_p_perseqid_".$t_seqid_cutoff}) {$ligstats->{$sid12_side}->{"max_l_and_p_perseqid_".$t_seqid_cutoff} = 0; }

               push @outvals, $ligstats->{$sid12_side}->{"cumlig_l_and_p_perseqid_".$t_seqid_cutoff} ;
               push @outvals, $ligstats->{$sid12_side}->{"max_l_and_p_perseqid_".$t_seqid_cutoff} ;
            }

         }
         print join("\t", @outvals)."\n"; 
      }
   }

# 2. pINT-LINT - iterate over all peptide involving this class
# new081231_1812 
   if ($mode eq 'pINT') {
      my $classtype = 'fam' ;
      my @headers_sum = ("pINT_LSUM", "SID", "CHAIN", "CHAIN_LENGTH", "CLASS",
                          "numres_p", "cumulative_numres_l_and_p",
                          "max_l_and_p", "lig_max_l_and_p",
                         ) ;
      foreach my $j ( 0 .. $hist2d_maxbin->{ligbs_seqid}) {
         my $t_seqid_cutoff = $j * $hist2d_binsize->{ligbs_seqid} ;
         push @headers_sum, "cum_l_and_p_perseqid_".$t_seqid_cutoff;
         push @headers_sum, "max_l_and_p_perseqid_".$t_seqid_cutoff;
      }
      print '#'.join("\t",@headers_sum)."\n" ;

      @headers_sum = ("SID", "CHAIN", "CHAIN_LENGTH", "CLASS", "LIG",
                      "numres_p", "numres_p_IDENT",
                      "numres_l", "numres_l_IDENT",
                      "numres_l_and_p", "numres_l_and_p_side_IDENT",
                      "numres_domain_nongap", "numres_domain_nongap_IDENT",
                     ) ;
      print '#'.join("\t",@headers_sum)."\n" ;

# group them by class so that can calc seqids - rejigger like did to PINT-LINT
      my $class2pint;
      foreach my $sid (keys %{$pepnucibits->{pep}->{$classtype}}) {
         $class2pint->{$pb->{sid2class}->{$classtype}->{$sid}}->{$sid}++ ; }

      foreach my $class (keys %{$class2pint}) {
         my $curclasstype = 'fam' ; my $curfam = $class ;

         if ($curfam !~ /[a-g]/) {next;} #ASTRAL alignments only available
                                         #for class a-g

# load ASTRAL alignment for this family
         my $class_aln = pibase::ASTRAL::load_asteroids_aln({
                  aln_fn => $pibase_specs->{asteroids}->{$curclasstype.'_aln'}.
                  '/'.$curfam.'.fasta_aln' ,
                  seq_fn => $pibase_specs->{asteroids}->{$curclasstype.'_seq'}.
                  '/'.$curfam.'.fa' ,
                  raf => $astral->{raf},
                  gdseqh => $astral->{gdseqh},
                  seqclcont100 => $astral->{seqcl2cont}->{100},
                  seqcl100 => $astral->{seqcl}->{100},
                  allchains => $pb->{pdbchains}
         }) ;

         foreach my $sid (keys %{$class2pint->{$class}}){
          my ($sid_origdom_p) = ($sid =~ /SCOP\.(.+)/) ;

          foreach my $targetch (
             keys %{$pepnucibits->{pep}->{$classtype}->{$sid}} ){

             if ($targetch eq 'cumulative') {next;}

         my $pdb = $pb->{sid2pdb}->{$sid} ;
         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;
         my $alnlength = $class2alnlength->{$classtype}->{$class};
   
         my $bs_bits = $pepnucibits->{pep}->{$classtype}->{$sid}->{$targetch};
         my $curligs ;
   
         my $pistats ;
         my $ligstats ; 
         if (exists $class2allligbits->{$classtype}->{$class}) {
            $curligs = $class2allligbits->{$classtype}->{$class} ; }
   
# init vals
         foreach my $type (qw/p cumlig_l cumlig_l_and_p max_l_and_p/) {
            $ligstats->{$type} = 0 ; }
         foreach my $type (qw/lig_max_l_and_p/) {
            $ligstats->{$type} = 'U'; }

         foreach my $type (qw/p cumlig_l cumlig_l_and_p max_l_and_p/) {
            $ligstats->{$type} = 0 ; }
         foreach my $type (qw/lig_max_l_and_p/) {
            $ligstats->{$type} = 'U'; }

         foreach my $type (qw/cumlig_l_and_p_perseqid max_l_and_p_perseqid/) {
            my $seqid_cutoff = 0 ;
            while ($seqid_cutoff < 100) {
               $ligstats->{$type."_perseqid_$seqid_cutoff"}=0;
               $seqid_cutoff += $hist2d_binsize->{ligbs_seqid} ;
            }
         }
   
         $pistats->{p} = $bs_bits->Norm() ;
         $ligstats->{cumlig} = Bit::Vector->new($alnlength) ;
         $ligstats->{"max_l_and_p"} = 0 ;
         $ligstats->{"lig_max_l_and_p"} = 'undef' ;

         my $seqid_cutoff = 0;
         while ($seqid_cutoff < 100) {
            $ligstats->{"cumlig_perseqid_$seqid_cutoff"} =
               Bit::Vector->new($alnlength) ;
            $seqid_cutoff += $hist2d_binsize->{ligbs_seqid} ;
         }

   
            my $temp_bits_and = Bit::Vector->new($alnlength) ;
            my $temp_bits_or = Bit::Vector->new($alnlength) ;
   
            foreach my $j (0 .. $#{$curligs}) {
               my $curligsig = $curligs->[$j]->[0] ;
               my $outcurligsig = $curligsig; $outcurligsig =~ s/\t/:/g ;
               my $curligbits = $curligs->[$j]->[1] ;

               my ($sid_origdom_l) =
                  ($curligs->[$j]->[2] =~ /SCOP\.(.+)/) ;
               
               $ligstats->{cumlig}->Or($ligstats->{cumlig}, $curligbits) ;
   
#intersect the ligands bit vector and the bit vector for the bindin
# site positions of this particular interface;

               $temp_bits_and->And($curligbits, $bs_bits) ;
               $temp_bits_or->Or($curligbits, $bs_bits) ;
   
               my $curs ;
               $curs->{l_and_p} = $temp_bits_and->Norm() ;
               $curs->{l_or_p} = $temp_bits_or->Norm() ;
               $curs->{l} = $curligbits->Norm() ;
               $curs->{p} = $bs_bits->Norm() ;
   
               if ($curs->{l_and_p} == 0) {next;}

# compute pairwise sequence identity over intersected residues.
# is the ASTRAL aln loaded 0- or 1-indexed?
               $curs->{l_and_p_IDENT} = 0 ;
               $curs->{l_IDENT} = 0 ;
               $curs->{p_IDENT} = 0 ;

               foreach my $aln_pos ($temp_bits_and->Index_List_Read()) {
                  if (substr($class_aln->{aln}->{$sid_origdom_l}, $aln_pos,1) eq
                      substr($class_aln->{aln}->{$sid_origdom_p}, $aln_pos,1)) {
                      $curs->{l_and_p_IDENT}++ ; } }

               foreach my $aln_pos ($curligbits->Index_List_Read()) {
                  if (substr($class_aln->{aln}->{$sid_origdom_l}, $aln_pos,1) eq
                      substr($class_aln->{aln}->{$sid_origdom_p}, $aln_pos,1)) {
                      $curs->{l_IDENT}++ ; } }

               foreach my $aln_pos ($bs_bits->Index_List_Read()) {
                  if (substr($class_aln->{aln}->{$sid_origdom_l}, $aln_pos,1) eq
                      substr($class_aln->{aln}->{$sid_origdom_p}, $aln_pos,1)) {
                      $curs->{p_IDENT}++ ; } }

               $curs->{wholedom_numnongap} = 0 ;
               $curs->{wholedom_IDENT} = 0 ;
               foreach my $aln_pos (0 .. $alnlength) {
                  my $a= substr($class_aln->{aln}->{$sid_origdom_l},$aln_pos,1);
                  my $b= substr($class_aln->{aln}->{$sid_origdom_p},$aln_pos,1);
                  if ($a ne '-' || $b ne '-') {
                     $curs->{wholedom_numnongap}++ ;
                     if ($a eq $b) { $curs->{wholedom_IDENT}++ ; }
                  }
               }
   
               if (!exists $ligstats->{"max_l_and_p"} ||
                   $curs->{l_and_p} > $ligstats->{"max_l_and_p"}) {
                  $ligstats->{"max_l_and_p"} = $curs->{l_and_p} ;
                  $ligstats->{"lig_max_l_and_p"} = $outcurligsig ;
               }

               my $cur_ligbs_seqid = ($curs->{l_IDENT} / $curs->{l}) * 100 ;
               my $seqid_cutoff = 0 ;
               my $seqid_bin_ceil = POSIX::floor($cur_ligbs_seqid /
                                        $hist2d_binsize->{ligbs_seqid});

               if ($seqid_bin_ceil > $hist2d_maxbin->{ligbs_seqid}) {
                   $seqid_bin_ceil = $hist2d_maxbin->{ligbs_seqid}; }

               foreach my $t_seqid_bin (0 ..  $seqid_bin_ceil) {
                  my $t_seqid_cutoff = $t_seqid_bin *
                                       $hist2d_binsize->{ligbs_seqid} ;
                 if (!exists $ligstats->{"max_l_and_p_perseqid_$t_seqid_cutoff"}
                 || $curs->{l_and_p} >
                     $ligstats->{"max_l_and_p_perseqid_$t_seqid_cutoff"}){
                  $ligstats->{"max_l_and_p_perseqid_$t_seqid_cutoff"} =
                     $curs->{l_and_p};
                  }

                  $ligstats->{"cumlig_perseqid_$t_seqid_cutoff"}->Or(
                     $ligstats->{"cumlig_perseqid_$t_seqid_cutoff"},
                     $curligbits) ;
               }
   
               my @outvals = ($sid, $targetch,
                       $pepnucibits->{chain_info}->{$targetch}->{chain_length},
                       $class,
                       $outcurligsig,
                       $curs->{p}, $curs->{p_IDENT},
                       $curs->{l}, $curs->{l_IDENT},
                       $curs->{l_and_p}, $curs->{l_and_p_IDENT},
                       $curs->{wholedom_numnongap}, $curs->{wholedom_IDENT},
               ) ;

               print join("\t", @outvals)."\n";
            }
   
            $ligstats->{"cumlig_l"} = $ligstats->{cumlig}->Norm() ;
            $temp_bits_and->And($ligstats->{cumlig}, $bs_bits) ;
            $ligstats->{"cumlig_l_and_p"} = $temp_bits_and->Norm() ;

            foreach my $t_seqid_bin (0 ..  $hist2d_maxbin->{ligbs_seqid}) {
               my $t_seqid_cutoff = $t_seqid_bin *
                                    $hist2d_binsize->{ligbs_seqid};
               $temp_bits_and->And(
                  $ligstats->{"cumlig_perseqid_$t_seqid_cutoff"},
                  $bs_bits) ;
               $ligstats->{"cumlig_l_and_p_perseqid_$t_seqid_cutoff"} =
                  $temp_bits_and->Norm() ;
            }
   
            my @outvals =  ("pINT_LSUM", $sid, $targetch,
                  $pepnucibits->{chain_info}->{$targetch}->{chain_length},
                  $class) ;
   
            if ($ligstats->{"cumlig_l"} == 0 ) {
               push @outvals, $pistats->{p}, 0 ;
               push @outvals, 0, "U" ;
            } else {
               push @outvals, $pistats->{p} ;
               push @outvals, $ligstats->{cumlig_l_and_p} ;

               push @outvals, $ligstats->{"max_l_and_p"} ;
               push @outvals, $ligstats->{"lig_max_l_and_p"} ;
            }

            foreach my $t_seqid_bin (0 .. $hist2d_maxbin->{ligbs_seqid}) {
               my $t_seqid_cutoff = $t_seqid_bin * $hist2d_binsize->{ligbs_seqid} ;
               if (!exists $ligstats->{"cumlig_l_and_p_perseqid_".$t_seqid_cutoff}) {
                  $ligstats->{"cumlig_l_and_p_perseqid_".$t_seqid_cutoff} = 0; }
               if (!exists $ligstats->{"max_l_and_p_perseqid_".$t_seqid_cutoff}) {
                  $ligstats->{"max_l_and_p_perseqid_".$t_seqid_cutoff} = 0; }

               push @outvals,
                  $ligstats->{"cumlig_l_and_p_perseqid_".$t_seqid_cutoff},
                  $ligstats->{"max_l_and_p_perseqid_".$t_seqid_cutoff};
            }

            print join("\t", @outvals)."\n"; 

          }
         }
      }


   }

}

sub BK_PRE090303_1855_collate_perinstance {

   require Bit::Vector ;

   my $standardres = _list_standardres() ;

# Per-PPI
# 0. load all of ligands first - lig info
# - iterate over assignments, and have individual
#   bit vectors set for the ligand binding sites in their
#   respective domain families
#
# per-PPI: iterate over each domain--domain interface
#  1. list any ligands actually bound to this interface
#  2. iterate over all ligands in that family/superfamily, find
#     - quantify per ligand coverage
#     - quantify cumulative ligand coverage
#
# per-ligand:


   my $pibase_specs = pibase::get_specs() ;
   my $pilig_specs = set_pilig_specs() ;
   my $fn = $pilig_specs->{fn} ;

   my $astral = _pilig_astral_preload({
      pibase_specs => $pibase_specs,
      pilig_specs => $pilig_specs,
   }) ;
   my $pb = _pilig_tod_pibase_preload() ;
   my $liginfo = _pilig_load_liginfo() ;

   my $class2alnlength = {};

   my $expbits = readin_expassignments({
      fn => $pilig_specs->{outfiles}->{assign_exp},
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;

   my $pepnucibits = readin_pepnuciassignments({
      fn => $pilig_specs->{outfiles}->{assign_pepnuci_clusters},
      clustrep_fl => 1,
      pb => $pb,
      pilig_specs => $pilig_specs,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;
   my $ligbits = readin_ligassignments({
      fn => $pilig_specs->{outfiles}->{assign_lig_clusters},
      clustrep_fl => 1,
      liginfo => $liginfo,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;

   my $pibits_both = readin_piassignments({
      fn => $pilig_specs->{outfiles}->{assign_pi_clusters},
      clustrep_fl => 1,
      pb => $pb,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;
   my $pibits = $pibits_both->{pibits} ;
   my $interfaces = $pibits_both->{interfaces} ;

   my $alltogether = combine_piligpepexpbits({
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


# 1. PINT-LINT - iterate over all domain interfaces involving this class
#   foreach my $classtype (qw/fam/)
   if (1) {
      my $classtype = 'fam' ;

      my @headers_sum = ("PINT_LSUM", "SID1", "SID2", "CHAINS",
                          "CLASS1", "CLASS2",
                          "numres_p_1", "cumulative_numres_l_and_p_1",
                          "max_l_and_p_1", "lig_max_l_and_p_1",
                          "numres_p_2", "cumulative_numres_l_and_p_2",
                          "max_l_and_p_2", "lig_max_l_and_p_2",
                         ) ;
      print '#'.join("\t",@headers_sum)."\n" ;

      @headers_sum = ("SID1", "SID2", "CHAINS",
                          "CLASS1", "CLASS2",
                          "side", "SID", "LIG",
                          "numres_p_side",
                          "numres_l",
                          "numres_l_and_p_side"
                         ) ;
      print '#'.join("\t",@headers_sum)."\n" ;

#   foreach my $classtype (qw/fam/)
      foreach my $sid12 (keys %{$interfaces}) {
         my $sid ;
         ($sid->{1}, $sid->{2}) = split(/\t/, $sid12) ;
         my $pdb = $interfaces->{$sid12}->{pdb} ;
         my $chains = $interfaces->{$sid12}->{chains} ;
   
         my ($alnlength, $class) ;
         foreach my $side (1, 2) {
            $class->{$side} = $pb->{sid2class}->{$classtype}->{$sid->{$side}} ;
            if (!exists $class2alnlength->{$classtype}->{$class->{$side}}) {
               next ; }
            $alnlength->{$side} = $class2alnlength->{$classtype}->{$class->{$side}};
#            $surf_pos->{$side} = $class2bits->{$classtype}->{$class->{$side}}->{'E'}->Norm() ;
         }
   
         my $bs_bits;
         my $curligs ;
   
         my $alnpos ;
         my $pistats ;
         my $ligstats ; 
         my $skipside ;
         foreach my $s (1, 2) {
            if( !exists $interfaces->{$sid12}->{$s} ||
               !exists $interfaces->{$sid12}->{$s}->{pibits}->{$classtype} ||
               ($interfaces->{$sid12}->{$s}->{pibits}->{$classtype}->Norm() == 0)){
   
               $skipside->{$s}++ ;
               print STDERR "WARNING: skipping $sid->{1} -- $sid->{2}".
                  " undefined binding site alignment positions for $sid->{$s}\n";
               next;
            }
   
#            if (exists $class2anligbits->{$classtype}->{$class->{$s}}) {
#               $curligs->{$s} = $class2anligbits->{$classtype}->{$class->{$s}};}

            if (exists $class2allligbits->{$classtype}->{$class->{$s}}) {
               $curligs->{$s}= $class2allligbits->{$classtype}->{$class->{$s}};}

            $bs_bits->{$s} = $interfaces->{$sid12}->{$s}->{pibits}->{$classtype} ;
         }
   
# init vals
         foreach my $s (1, 2) {
            foreach my $type (qw/p cumlig_l cumlig_l_and_p max_l_and_p/) {
               $ligstats->{$s}->{$type} = 0 ; }
            foreach my $type (qw/lig_max_l_and_p/) {
               $ligstats->{$s}->{$type} = 'U'; }
         }
   
         foreach my $s (1, 2) {
            if (exists $skipside->{$s}) {next;}
   
            $pistats->{$s}->{p} = $bs_bits->{$s}->Norm() ;
            $ligstats->{$s}->{cumlig} = Bit::Vector->new($alnlength->{$s}) ;
            foreach my $t (qw/l_and_p/) {
               $ligstats->{$s}->{"max_$t"} = 0 ;
               $ligstats->{$s}->{"lig_max_$t"} = 'undef' ;
            }
   
            my $temp_bits_and = Bit::Vector->new($alnlength->{$s}) ;
            my $temp_bits_or = Bit::Vector->new($alnlength->{$s}) ;
   
            if (!exists $curligs->{$s}) {next;}
            foreach my $j (0 .. $#{$curligs->{$s}}) {
               my $curligsig = $curligs->{$s}->[$j]->[0] ;
               my $outcurligsig = $curligsig; $outcurligsig =~ s/\t/:/g ;
               my $curligbits = $curligs->{$s}->[$j]->[1] ;
               
               $ligstats->{$s}->{cumlig}->Or($ligstats->{$s}->{cumlig},
                                             $curligbits) ;
   
               #intersect the ligands bit vector and the bit vector for the bindin
               # site positions of this particular interface;
               $temp_bits_and->And($curligbits, $bs_bits->{$s}) ;
               $temp_bits_or->Or($curligbits, $bs_bits->{$s}) ;
   
               my $curs ;
               $curs->{l_and_p} = $temp_bits_and->Norm() ;
               $curs->{l_or_p} = $temp_bits_or->Norm() ;
               $curs->{l} = $curligbits->Norm() ;
               $curs->{p} = $bs_bits->{$s}->Norm() ;
   
               if ($curs->{l_and_p} == 0) {next;}
   
               my $curoverlap ;
               $curoverlap->{obi} = $curs->{l_and_p} / $curs->{l_or_p} ;
               $curoverlap->{opi} = $curs->{l_and_p} / $curs->{p} ;
               $curoverlap->{olig} = $curs->{l_and_p} / $curs->{l} ;
   
   
               if (!exists $ligstats->{$s}->{"max_l_and_p"} ||
                   $curs->{l_and_p} > $ligstats->{$s}->{"max_l_and_p"}) {
                  $ligstats->{$s}->{"max_l_and_p"} = $curs->{l_and_p} ;
                  $ligstats->{$s}->{"lig_max_l_and_p"} = $outcurligsig ;
               }
   
               my (undef, $ligcod, undef) = split(/\t/, $curligsig) ;
               my @outvals = ($sid->{1}, $sid->{2}, $chains,
                              $class->{1}, $class->{2},
                              $s, $sid->{$s}, $outcurligsig,
                              $curs->{p}, $curs->{l},
                              $curs->{l_and_p}
                             ) ;
               print join("\t", @outvals)."\n";
            }
   
            $ligstats->{$s}->{"cumlig_l"} = $ligstats->{$s}->{cumlig}->Norm() ;
            $temp_bits_and->And($ligstats->{$s}->{cumlig}, $bs_bits->{$s}) ;
            $ligstats->{$s}->{"cumlig_l_and_p"} = $temp_bits_and->Norm() ;
         }
   
         my @outvals =  ("PINT_LSUM", $sid->{1}, $sid->{2}, $chains,
                         $class->{1}, $class->{2}) ;
   
         foreach my $s (1, 2) {
            if (exists $skipside->{$s} ||
               !exists $ligstats->{$s} ||
               $ligstats->{$s}->{"cumlig_l"} == 0 ) {
   
               if (exists $skipside->{$s}) {
                  push @outvals, 0, 0 ;
               } else {
                  push @outvals, $pistats->{$s}->{p}, 0 ;
               }
   
               push @outvals, 0, "U" ;
               next;
            }
   
            push @outvals, $pistats->{$s}->{p} ;
            push @outvals, $ligstats->{$s}->{cumlig_l_and_p} ;
            push @outvals, $ligstats->{$s}->{"max_l_and_p"} ;
            push @outvals, $ligstats->{$s}->{"lig_max_l_and_p"} ;
         }
         print join("\t", @outvals)."\n"; 
      }
   }

# 2. pINT-LINT - iterate over all peptide involving this class
# new081231_1812 
   if (0) {
      my @headers_sum = ("pINT_LSUM", "SID", "CHAIN", "CHAIN_LENGTH", "CLASS",
                          "numres_p", "cumulative_numres_l_and_p",
                          "max_l_and_p", "lig_max_l_and_p",
                         ) ;
      print '#'.join("\t",@headers_sum)."\n" ;

      @headers_sum = ("SID", "CHAIN", "CHAIN_LENGTH", "CLASS", "LIG",
                      "numres_p", "numres_l", "numres_l_and_p") ;
      print '#'.join("\t",@headers_sum)."\n" ;

      my $classtype = 'fam' ;
      foreach my $sid (keys %{$pepnucibits->{pep}->{$classtype}}) {
       foreach my $targetch (keys %{$pepnucibits->{pep}->{$classtype}->{$sid}}){
         if ($targetch eq 'cumulative') {next;}
         print "#IN ON $sid $targetch\n" ;

         my $pdb = $pb->{sid2pdb}->{$sid} ;
         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;
         my $alnlength = $class2alnlength->{$classtype}->{$class};
   
         my $bs_bits = $pepnucibits->{pep}->{$classtype}->{$sid}->{$targetch};
         my $curligs ;
   
         my $alnpos ;
         my $pistats ;
         my $ligstats ; 
         if (exists $class2allligbits->{$classtype}->{$class}) {
            $curligs = $class2allligbits->{$classtype}->{$class} ; }
   
# init vals
         foreach my $type (qw/p cumlig_l cumlig_l_and_p max_l_and_p/) {
            $ligstats->{$type} = 0 ; }
         foreach my $type (qw/lig_max_l_and_p/) {
            $ligstats->{$type} = 'U'; }
   
            $pistats->{p} = $bs_bits->Norm() ;
            $ligstats->{cumlig} = Bit::Vector->new($alnlength) ;
            foreach my $t (qw/l_and_p/) {
               $ligstats->{"max_$t"} = 0 ;
               $ligstats->{"lig_max_$t"} = 'undef' ;
            }
   
            my $temp_bits_and = Bit::Vector->new($alnlength) ;
            my $temp_bits_or = Bit::Vector->new($alnlength) ;
   
            foreach my $j (0 .. $#{$curligs}) {
               my $curligsig = $curligs->[$j]->[0] ;
               my $outcurligsig = $curligsig; $outcurligsig =~ s/\t/:/g ;
               my $curligbits = $curligs->[$j]->[1] ;
               
               $ligstats->{cumlig}->Or($ligstats->{cumlig},
                                             $curligbits) ;
   
               #intersect the ligands bit vector and the bit vector for the bindin
               # site positions of this particular interface;
               $temp_bits_and->And($curligbits, $bs_bits) ;
               $temp_bits_or->Or($curligbits, $bs_bits) ;
   
               my $curs ;
               $curs->{l_and_p} = $temp_bits_and->Norm() ;
               $curs->{l_or_p} = $temp_bits_or->Norm() ;
               $curs->{l} = $curligbits->Norm() ;
               $curs->{p} = $bs_bits->Norm() ;
   
               if ($curs->{l_and_p} == 0) {next;}
   
               if (!exists $ligstats->{"max_l_and_p"} ||
                   $curs->{l_and_p} > $ligstats->{"max_l_and_p"}) {
                  $ligstats->{"max_l_and_p"} = $curs->{l_and_p} ;
                  $ligstats->{"lig_max_l_and_p"} = $outcurligsig ;
               }
   
               my (undef, $ligcod, undef) = split(/\t/, $curligsig) ;
   
               my @outvals = ($sid, $targetch,
                       $pepnucibits->{chain_info}->{$targetch}->{chain_length},
                       $class,
                       $outcurligsig,
                       $curs->{p},
                       $curs->{l},
                       $curs->{l_and_p}
               ) ;

               print join("\t", @outvals)."\n";
            }
   
         $ligstats->{"cumlig_l"} = $ligstats->{cumlig}->Norm() ;
         $temp_bits_and->And($ligstats->{cumlig}, $bs_bits) ;
         $ligstats->{"cumlig_l_and_p"} = $temp_bits_and->Norm() ;
   
         my @outvals =  ("#pINT_LSUM", $sid, $targetch,
            $pepnucibits->{chain_info}->{$targetch}->{chain_length},
                                 $class) ;
   
         if ($ligstats->{"cumlig_l"} == 0 ) {
            push @outvals, $pistats->{p}, 0 ;
            push @outvals, 0, "U" ;
         } else {
            push @outvals, $pistats->{p} ;
            push @outvals, $ligstats->{cumlig_l_and_p} ;

            push @outvals, $ligstats->{"max_l_and_p"} ;
            push @outvals, $ligstats->{"lig_max_l_and_p"} ;
         }
         print join("\t", @outvals)."\n"; 
       }
      }
   }

#3. NINT-LINT 
   if (0) {
      my $classtype = 'fam' ;
      foreach my $sid (keys %{$pepnucibits->{nuc}->{$classtype}}) {
       foreach my $targetch (keys %{$pepnucibits->{nuc}->{$classtype}->{$sid}}){
         if ($targetch eq 'cumulative') {next;}

         my $pdb = $pb->{sid2pdb}->{$sid} ;
         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;
         my $alnlength = $class2alnlength->{$classtype}->{$class};
   
         my $bs_bits = $pepnucibits->{nuc}->{$classtype}->{$sid}->{$targetch};
         my $curligs ;
   
         my $alnpos ;
         my $pistats ;
         my $ligstats ; 
#         if (exists $class2anligbits->{$classtype}->{$class}) {
#            $curligs = $class2anligbits->{$classtype}->{$class} ; }
         if (exists $class2allligbits->{$classtype}->{$class}) {
            $curligs = $class2allligbits->{$classtype}->{$class} ; }
   
# init vals
         foreach my $type (qw/p cumlig_l cumlig_l_and_p max_liginrange max_l_and_p max_obi max_opi max_olig/) {
            $ligstats->{$type} = 0 ; }
         foreach my $type (qw/lig_max_l_and_p lig_max_obi lig_max_opi lig_maxolig/) {
            $ligstats->{$type} = 'U'; }
   
            $pistats->{p} = $bs_bits->Norm() ;
            $pistats->{max_liginrange} = 0 ;
            $ligstats->{cumlig} = Bit::Vector->new($alnlength) ;
            foreach my $t (qw/obi opi olig l_and_p/) {
               $ligstats->{"max_$t"} = 0 ;
               $ligstats->{"lig_max_$t"} = 'undef' ;
            }
   
            my $temp_bits_and = Bit::Vector->new($alnlength) ;
            my $temp_bits_or = Bit::Vector->new($alnlength) ;
   
            foreach my $j (0 .. $#{$curligs}) {
               my $curligsig = $curligs->[$j]->[0] ;
               my $outcurligsig = $curligsig; $outcurligsig =~ s/\t/:/g ;
               my $curligbits = $curligs->[$j]->[1] ;
               
               $ligstats->{cumlig}->Or($ligstats->{cumlig},
                                             $curligbits) ;
   
               #intersect the ligands bit vector and the bit vector for the bindin
               # site positions of this particular interface;
               $temp_bits_and->And($curligbits, $bs_bits) ;
               $temp_bits_or->Or($curligbits, $bs_bits) ;
   
               my $curs ;
               $curs->{l_and_p} = $temp_bits_and->Norm() ;
               $curs->{l_or_p} = $temp_bits_or->Norm() ;
               $curs->{l} = $curligbits->Norm() ;
               $curs->{p} = $bs_bits->Norm() ;
   
               if ($curs->{l_and_p} == 0) {next;}
   
               my $curoverlap ;
               $curoverlap->{obi} = $curs->{l_and_p} / $curs->{l_or_p} ;
               $curoverlap->{opi} = $curs->{l_and_p} / $curs->{p} ;
               $curoverlap->{olig} = $curs->{l_and_p} / $curs->{l} ;
   
   
               if (!exists $ligstats->{"max_l_and_p"} ||
                   $curs->{l_and_p} > $ligstats->{"max_l_and_p"}) {
                  $ligstats->{"max_l_and_p"} = $curs->{l_and_p} ;
                  $ligstats->{"lig_max_l_and_p"} = $outcurligsig ;
               }
   
               foreach my $t (qw/obi opi olig/) {
                  if (!exists $ligstats->{"max_$t"} ||
                         $curoverlap->{$t} > $ligstats->{"max_$t"}) {
                     $ligstats->{"max_$t"} = $curoverlap->{$t} ;
                     $ligstats->{"lig_max_$t"} = $outcurligsig ;
                  }
               }
   
               my (undef, $ligcod, undef) = split(/\t/, $curligsig) ;
               my $liginrange = 0 ;
               if (exists $liginfo->{mwinrange}->{$ligcod}) {
                   $ligstats->{max_liginrange} = 1 ;
                   $liginrange = 1 ; }
   
               my @outvals = ("nINT_LINT",
                                 $pdb, $sid, $targetch,
                       $pepnucibits->{chain_info}->{$targetch}->{chain_length},
                                 $classtype,
                                 $class,
                                 $outcurligsig,
                                 $liginrange,
                                 $curs->{p}, $curs->{l},
                                 $curs->{l_and_p}, $curs->{l_or_p},
                                 sprintf("%.3f", $curoverlap->{obi}),
                                 sprintf("%.3f", $curoverlap->{opi}),
                                 sprintf("%.3f", $curoverlap->{olig})
                              ) ;
               print join("\t", @outvals)."\n";
            }
   
         $ligstats->{"cumlig_l"} = $ligstats->{cumlig}->Norm() ;
         $temp_bits_and->And($ligstats->{cumlig}, $bs_bits) ;
         $ligstats->{"cumlig_l_and_p"} = $temp_bits_and->Norm() ;
   
         my @outvals =  ("nINT_LSUM", $pdb, $sid, $targetch,
            $pepnucibits->{chain_info}->{$targetch}->{chain_length},
                                 $classtype, $class) ;
   
         if ($ligstats->{"cumlig_l"} == 0 ) {
               push @outvals, $pistats->{p}, 0, 0 ;
               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
         } else {
   
            push @outvals, $pistats->{p} ;
            push @outvals, $ligstats->{cumlig_l_and_p} ;
            push @outvals, $ligstats->{max_liginrange} ;
            push @outvals, $ligstats->{"max_l_and_p"} ;
            push @outvals, $ligstats->{"lig_max_l_and_p"} ;
            foreach my $t (qw/obi opi olig/)  {
               push @outvals, sprintf("%.3f", $ligstats->{"max_$t"}) ;
               push @outvals, $ligstats->{"lig_max_$t"} ;
            }
         }
         print join("\t", @outvals)."\n"; 
       }
      }
   }

#4. NINT-pINT  - HERENOW HAVE TO EDIT STILL
   if (0) {
      my $classtype = 'fam' ;
#   foreach my $classtype (qw/fam/)
      foreach my $sid (keys %{$pepnucibits->{nuc}->{$classtype}}) {
       foreach my $targetch (keys %{$pepnucibits->{nuc}->{$classtype}->{$sid}}){
         if ($targetch eq 'cumulative') {next;}

         my $pdb = $pb->{sid2pdb}->{$sid} ;
         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;
         my $alnlength = $class2alnlength->{$classtype}->{$class};
   
         my $bs_bits = $pepnucibits->{nuc}->{$classtype}->{$sid}->{$targetch};
         my $curligs ;
   
         my $alnpos ;
         my $pistats ;
         my $ligstats ; 
         if (exists $class2allligbits->{$classtype}->{$class}) {
            $curligs = $class2allligbits->{$classtype}->{$class} ; }
   
# init vals
         foreach my $type (qw/p cumlig_l cumlig_l_and_p max_liginrange max_l_and_p max_obi max_opi max_olig/) {
            $ligstats->{$type} = 0 ; }
         foreach my $type (qw/lig_max_l_and_p lig_max_obi lig_max_opi lig_maxolig/) {
            $ligstats->{$type} = 'U'; }
   
            $pistats->{p} = $bs_bits->Norm() ;
            $pistats->{max_liginrange} = 0 ;
            $ligstats->{cumlig} = Bit::Vector->new($alnlength) ;
            foreach my $t (qw/obi opi olig l_and_p/) {
               $ligstats->{"max_$t"} = 0 ;
               $ligstats->{"lig_max_$t"} = 'undef' ;
            }
   
            my $temp_bits_and = Bit::Vector->new($alnlength) ;
            my $temp_bits_or = Bit::Vector->new($alnlength) ;
   
            foreach my $j (0 .. $#{$curligs}) {
               my $curligsig = $curligs->[$j]->[0] ;
               my $outcurligsig = $curligsig; $outcurligsig =~ s/\t/:/g ;
               my $curligbits = $curligs->[$j]->[1] ;
               
               $ligstats->{cumlig}->Or($ligstats->{cumlig},
                                             $curligbits) ;
   
               #intersect the ligands bit vector and the bit vector for the bindin
               # site positions of this particular interface;
               $temp_bits_and->And($curligbits, $bs_bits) ;
               $temp_bits_or->Or($curligbits, $bs_bits) ;
   
               my $curs ;
               $curs->{l_and_p} = $temp_bits_and->Norm() ;
               $curs->{l_or_p} = $temp_bits_or->Norm() ;
               $curs->{l} = $curligbits->Norm() ;
               $curs->{p} = $bs_bits->Norm() ;
   
               if ($curs->{l_and_p} == 0) {next;}
   
               my $curoverlap ;
               $curoverlap->{obi} = $curs->{l_and_p} / $curs->{l_or_p} ;
               $curoverlap->{opi} = $curs->{l_and_p} / $curs->{p} ;
               $curoverlap->{olig} = $curs->{l_and_p} / $curs->{l} ;
   
   
               if (!exists $ligstats->{"max_l_and_p"} ||
                   $curs->{l_and_p} > $ligstats->{"max_l_and_p"}) {
                  $ligstats->{"max_l_and_p"} = $curs->{l_and_p} ;
                  $ligstats->{"lig_max_l_and_p"} = $outcurligsig ;
               }
   
               foreach my $t (qw/obi opi olig/) {
                  if (!exists $ligstats->{"max_$t"} ||
                         $curoverlap->{$t} > $ligstats->{"max_$t"}) {
                     $ligstats->{"max_$t"} = $curoverlap->{$t} ;
                     $ligstats->{"lig_max_$t"} = $outcurligsig ;
                  }
               }
   
               my (undef, $ligcod, undef) = split(/\t/, $curligsig) ;
               my $liginrange = 0 ;
               if (exists $liginfo->{mwinrange}->{$ligcod}) {
                   $ligstats->{max_liginrange} = 1 ;
                   $liginrange = 1 ; }
   
               my @outvals = ("nINT_LINT",
                                 $pdb, $sid, $targetch,
                       $pepnucibits->{chain_info}->{$targetch}->{chain_length},
                                 $classtype,
                                 $class,
                                 $outcurligsig,
                                 $liginrange,
                                 $curs->{p}, $curs->{l},
                                 $curs->{l_and_p}, $curs->{l_or_p},
                                 sprintf("%.3f", $curoverlap->{obi}),
                                 sprintf("%.3f", $curoverlap->{opi}),
                                 sprintf("%.3f", $curoverlap->{olig})
                              ) ;
               print join("\t", @outvals)."\n";
            }
   
         $ligstats->{"cumlig_l"} = $ligstats->{cumlig}->Norm() ;
         $temp_bits_and->And($ligstats->{cumlig}, $bs_bits) ;
         $ligstats->{"cumlig_l_and_p"} = $temp_bits_and->Norm() ;
   
         my @outvals =  ("nINT_LSUM", $pdb, $sid, $targetch,
            $pepnucibits->{chain_info}->{$targetch}->{chain_length},
                                 $classtype, $class) ;
   
         if ($ligstats->{"cumlig_l"} == 0 ) {
               push @outvals, $pistats->{p}, 0, 0 ;
               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
         } else {
   
            push @outvals, $pistats->{p} ;
            push @outvals, $ligstats->{cumlig_l_and_p} ;
            push @outvals, $ligstats->{max_liginrange} ;
            push @outvals, $ligstats->{"max_l_and_p"} ;
            push @outvals, $ligstats->{"lig_max_l_and_p"} ;
            foreach my $t (qw/obi opi olig/)  {
               push @outvals, sprintf("%.3f", $ligstats->{"max_$t"}) ;
               push @outvals, $ligstats->{"lig_max_$t"} ;
            }
         }
         print join("\t", @outvals)."\n"; 
       }
      }
   }

}

sub PRE20081231_collate_perinstance {

   require Bit::Vector ;

   my $standardres = {
   ALA => 'A' ,
   ARG => 'R' ,
   ASN => 'N' ,
   ASP => 'D' ,
   CYS => 'C' ,
   GLN => 'Q' ,
   GLU => 'E' ,
   GLY => 'G' ,
   HIS => 'H' ,
   HSD => 'H' ,
   HSE => 'H' ,
   ILE => 'I' ,
   LEU => 'L' ,
   LYS => 'K' ,
   MET => 'M' ,
   PHE => 'F' ,
   PRO => 'P' ,
   SER => 'S' ,
   THR => 'T' ,
   TRP => 'W' ,
   TYR => 'Y' ,
   VAL => 'V',
   UNK => 'X',
   '  C' => 'c',
   '  G' => 'g',
   '  A' => 'a',
   '  T' => 't',
   '  U' => 'u',
   '  I' => 'i',
   'C' => 'c',
   'G' => 'g',
   'A' => 'a',
   'T' => 't',
   'U' => 'u',
   'I' => 'i',
   '+C' => 'c',
   '+G' => 'g',
   '+A' => 'a',
   '+T' => 't',
   '+U' => 'u',
   '+I' => 'i'
   } ;

# Per-PPI
# 0. load all of ligands first - lig info
# - iterate over assignments, and have individual
#   bit vectors set for the ligand binding sites in their
#   respective domain families
#
# per-PPI: iterate over each domain--domain interface
#  1. list any ligands actually bound to this interface
#  2. iterate over all ligands in that family/superfamily, find
#     - quantify per ligand coverage
#     - quantify cumulative ligand coverage
#
# per-ligand:


   my $pibase_specs = pibase::get_specs() ;
   my $pilig_specs = set_pilig_specs() ;
   my $fn = $pilig_specs->{fn} ;

   my $astral = _pilig_astral_preload({
      pibase_specs => $pibase_specs,
      pilig_specs => $pilig_specs,
   }) ;
   my $pb = _pilig_tod_pibase_preload() ;
   my $liginfo = _pilig_load_liginfo() ;

   my $class2alnlength = {};

   my $ligbits = readin_ligassignments({
      fn => $pilig_specs->{outfiles}->{assign_lig},
      liginfo => $liginfo,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;

   my $expbits = readin_expassignments({
      fn => $pilig_specs->{outfiles}->{assign_exp},
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;

   my $pibits_both = readin_piassignments({
      fn => $pilig_specs->{outfiles}->{assign_pi},
      pb => $pb,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;
   my $pibits = $pibits_both->{pibits} ;
   my $interfaces = $pibits_both->{interfaces} ;

   my $alltogether = combine_piligexpbits({
      pb => $pb,
      liginfo => $liginfo,
      ligbits => $ligbits,
      pibits => $pibits,
      expbits => $expbits,
   });
   my $class2bits = $alltogether->{class2bits} ;
   my $class2ligs = $alltogether->{class2ligs} ;
   my $class2anligbits = $alltogether->{class2anligbits} ;


   my @headers_sum=  ("PINT_LSUM", "PDB", "SID1", "SID2", "CHAINS",
                       "CLASSTYPE", "CLASS1", "CLASS2",
                       "numres_p_1", "cumulative_numres_l_and_p_1",
                       "max_liginrange_1",
                       "max_l_and_p_1", "lig_max_l_and_p_1",
                       "max_obi_1", "lig_max_obi_1",
                       "max_opi_1", "lig_max_opi_1",
                       "max_olig_1", "lig_max_olig_1",
                       "numres_p_2", "cumulative_numres_l_and_p_2",
                       "max_liginrange_2",
                       "max_l_and_p_2", "lig_max_l_and_p_2",
                       "max_obi_2", "lig_max_obi_2",
                       "max_opi_2", "lig_max_opi_2",
                       "max_olig_2", "lig_max_olig_2",
                       ) ;
   print '#'.join("\t",@headers_sum)."\n" ;

#   foreach my $classtype (qw/sf fam/)
   foreach my $classtype (qw/fam/) {
      foreach my $sid12 (keys %{$interfaces}) {
         my $sid ;
         ($sid->{1}, $sid->{2}) = split(/\t/, $sid12) ;
         my $pdb = $interfaces->{$sid12}->{pdb} ;
         my $chains = $interfaces->{$sid12}->{chains} ;
   
         my ($alnlength, $class) ;
         foreach my $side (1, 2) {
            $class->{$side} = $pb->{sid2class}->{$classtype}->{$sid->{$side}} ;
            if (!exists $class2alnlength->{$classtype}->{$class->{$side}}) {
               next ; }
            $alnlength->{$side} = $class2alnlength->{$classtype}->{$class->{$side}};
#            $surf_pos->{$side} = $class2bits->{$classtype}->{$class->{$side}}->{'E'}->Norm() ;
         }
   
         my $bs_bits;
         my $curligs ;
   
         my $alnpos ;
         my $pistats ;
         my $ligstats ; 
         my $skipside ;
         foreach my $s (1, 2) {
            if( !exists $interfaces->{$sid12}->{$s} ||
               !exists $interfaces->{$sid12}->{$s}->{pibits}->{$classtype} ||
               ($interfaces->{$sid12}->{$s}->{pibits}->{$classtype}->Norm() == 0)){
   
               $skipside->{$s}++ ;
               print STDERR "WARNING: skipping $sid->{1} -- $sid->{2}".
                  " undefined binding site alignment positions for $sid->{$s}\n";
               next;
            }
   
            if (exists $class2anligbits->{$classtype}->{$class->{$s}}) {
               $curligs->{$s} =
                  $class2anligbits->{$classtype}->{$class->{$s}} ; }
            $bs_bits->{$s} = $interfaces->{$sid12}->{$s}->{pibits}->{$classtype} ;
         }
   
   # init vals
         foreach my $s (1, 2) {
            foreach my $type (qw/p cumlig_l cumlig_l_and_p max_liginrange max_l_and_p max_obi max_opi max_olig/) {
               $ligstats->{$s}->{$type} = 0 ; }
            foreach my $type (qw/lig_max_l_and_p lig_max_obi lig_max_opi lig_maxolig/) {
               $ligstats->{$s}->{$type} = 'U'; }
         }
   
         foreach my $s (1, 2) {
            if (exists $skipside->{$s}) {next;}
   
            $pistats->{$s}->{p} = $bs_bits->{$s}->Norm() ;
            $pistats->{$s}->{max_liginrange} = 0 ;
            $ligstats->{$s}->{cumlig} = Bit::Vector->new($alnlength->{$s}) ;
            foreach my $t (qw/obi opi olig l_and_p/) {
               $ligstats->{$s}->{"max_$t"} = 0 ;
               $ligstats->{$s}->{"lig_max_$t"} = 'undef' ;
            }
   
            my $temp_bits_and = Bit::Vector->new($alnlength->{$s}) ;
            my $temp_bits_or = Bit::Vector->new($alnlength->{$s}) ;
   
            if (!exists $curligs->{$s}) {next;}
            foreach my $j (0 .. $#{$curligs->{$s}}) {
               my $curligsig = $curligs->{$s}->[$j]->[0] ;
               my $outcurligsig = $curligsig; $outcurligsig =~ s/\t/:/g ;
               my $curligbits = $curligs->{$s}->[$j]->[1] ;
               
               $ligstats->{$s}->{cumlig}->Or($ligstats->{$s}->{cumlig},
                                             $curligbits) ;
   
               #intersect the ligands bit vector and the bit vector for the bindin
               # site positions of this particular interface;
               $temp_bits_and->And($curligbits, $bs_bits->{$s}) ;
               $temp_bits_or->Or($curligbits, $bs_bits->{$s}) ;
   
               my $curs ;
               $curs->{l_and_p} = $temp_bits_and->Norm() ;
               $curs->{l_or_p} = $temp_bits_or->Norm() ;
               $curs->{l} = $curligbits->Norm() ;
               $curs->{p} = $bs_bits->{$s}->Norm() ;
   
               if ($curs->{l_and_p} == 0) {next;}
   
               my $curoverlap ;
               $curoverlap->{obi} = $curs->{l_and_p} / $curs->{l_or_p} ;
               $curoverlap->{opi} = $curs->{l_and_p} / $curs->{p} ;
               $curoverlap->{olig} = $curs->{l_and_p} / $curs->{l} ;
   
   
               if (!exists $ligstats->{$s}->{"max_l_and_p"} ||
                   $curs->{l_and_p} > $ligstats->{$s}->{"max_l_and_p"}) {
                  $ligstats->{$s}->{"max_l_and_p"} = $curs->{l_and_p} ;
                  $ligstats->{$s}->{"lig_max_l_and_p"} = $outcurligsig ;
               }
   
               foreach my $t (qw/obi opi olig/) {
                  if (!exists $ligstats->{$s}->{"max_$t"} ||
                         $curoverlap->{$t} > $ligstats->{$s}->{"max_$t"}) {
                     $ligstats->{$s}->{"max_$t"} = $curoverlap->{$t} ;
                     $ligstats->{$s}->{"lig_max_$t"} = $outcurligsig ;
                  }
               }
   
               my (undef, $ligcod, undef) = split(/\t/, $curligsig) ;
               my $liginrange = 0 ;
               if (exists $liginfo->{mwinrange}->{$ligcod}) {
                   $ligstats->{$s}->{max_liginrange} = 1 ;
                   $liginrange = 1 ; }
   
               my @outvals = ("PINT_LINT",
                                 $pdb, $sid->{1}, $sid->{2},
                                 $chains, $classtype,
                                 $class->{1}, $class->{2},
                                 $s, $sid->{$s}, $outcurligsig,
                                 $liginrange,
                                 $curs->{p}, $curs->{l},
                                 $curs->{l_and_p}, $curs->{l_or_p},
                                 sprintf("%.3f", $curoverlap->{obi}),
                                 sprintf("%.3f", $curoverlap->{opi}),
                                 sprintf("%.3f", $curoverlap->{olig})
                              ) ;
               print join("\t", @outvals)."\n";
            }
   
            $ligstats->{$s}->{"cumlig_l"} = $ligstats->{$s}->{cumlig}->Norm() ;
            $temp_bits_and->And($ligstats->{$s}->{cumlig}, $bs_bits->{$s}) ;
            $ligstats->{$s}->{"cumlig_l_and_p"} = $temp_bits_and->Norm() ;
         }
   
         my @outvals =  ("PINT_LSUM", $pdb, $sid->{1}, $sid->{2}, $chains,
                                 $classtype, $class->{1}, $class->{2}) ;
   
         foreach my $s ( 1, 2) {
            if (exists $skipside->{$s} ||
               !exists $ligstats->{$s} ||
               $ligstats->{$s}->{"cumlig_l"} == 0 ) {
   
               if (exists $skipside->{$s}) {
                  push @outvals, 0, 0, 0 ;
               } else {
                  push @outvals, $pistats->{$s}->{p}, 0, 0 ;
               }
   
               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
               next;
            }
   
            push @outvals, $pistats->{$s}->{p} ;
            push @outvals, $ligstats->{$s}->{cumlig_l_and_p} ;
            push @outvals, $ligstats->{$s}->{max_liginrange} ;
            push @outvals, $ligstats->{$s}->{"max_l_and_p"} ;
            push @outvals, $ligstats->{$s}->{"lig_max_l_and_p"} ;
            foreach my $t (qw/obi opi olig/)  {
               push @outvals, sprintf("%.3f", $ligstats->{$s}->{"max_$t"}) ;
               push @outvals, $ligstats->{$s}->{"lig_max_$t"} ;
            }
         }
         print join("\t", @outvals)."\n"; 
      }
   }

}


sub readin_ligassignments {

   my $in = shift ;

   my $fn = $in->{fn} ;

   my $class2alnlength = $in->{class2alnlength};
   my $standardres= $in->{standardres};
   my $liginfo = $in->{liginfo};

# NOTE: different frmo the in->{cluster_fl} flag that actually _does_ the
#       clusternig and prints out memberships to the outfiles specified file

# if this is the clustered assignment list, only read in cluster representative
# (arbitrarily chosen - first one in the list)
   my $clustrep_fl = 0;
   my $clusters_seen ;  # keep track of what clusters have been seen already
   if (exists $in->{clustrep_fl} && $in->{clustrep_fl} == 1) {
      $clustrep_fl = 1 ; }


   my $class2sid ;

   my $ligbits ;
   print STDERR "NOW reading LIG assignment\n" ;
   open(LIGFH, $fn) ;
   while (my $line = <LIGFH>) {
      if ($line =~ /^#/ ||
          $line =~ /^Warning: no access to/ ||
          $line =~ /^Thus no job control/ ) {next;}

      chomp $line;
      my @fields = split(/\t/, $line) ;

      my $cluster_num ;
      if ($clustrep_fl) {
         $cluster_num = shift @fields ; }

      my ($pdb, $sid, $osid, $classtype, $class, $btype, $alnposstring,
          $alnlength, $ligcod, $ligid ) = @fields ;
      $class2sid->{$classtype}->{$class}->{$sid}++ ;

      $ligcod =~ s/ //g ;
      my $ligsig = $pdb."\t".$ligcod."\t".$ligid ;

      if (exists $standardres->{$ligcod}) {next;}
      if ($ligcod eq 'DUM' || $ligcod eq 'UNX' ||
             $ligcod eq 'UNK' || $ligcod eq 'UNL') {next;}

      if (!exists $liginfo->{mw}->{$ligcod}) {
#         print STDERR "WARNING: MW not found for ligand $ligcod\n" ;
         next; }

      if (!exists $liginfo->{mwinrange}->{$ligcod}) {
#         print STDERR "SKIPPING: MW out of range for ligand $ligcod\n";
         next;}


      if ($clustrep_fl) {
         if (!exists $clusters_seen->{$classtype}->{$cluster_num}) {
            $clusters_seen->{$classtype}->{$cluster_num}++ ;
# don't have to keep track of what was chosen, since only 1 line is needed
# if a line has been read, that cluster is done...
         } else {
            next;
         }
      }

      $class2alnlength->{$classtype}->{$class} = $alnlength ;

      $ligbits->{$classtype}->{$sid}->{$ligsig} = Bit::Vector->new($alnlength) ;

      my @alnpos = split(/\,/, $alnposstring) ;
      foreach my $alnpos (@alnpos) {
         $ligbits->{$classtype}->{$sid}->{$ligsig}->Bit_On($alnpos); }


      if (!exists $ligbits->{$classtype}->{$sid}->{cumulative}) {
         $ligbits->{$classtype}->{$sid}->{cumulative} =
            $ligbits->{$classtype}->{$sid}->{$ligsig}->Clone();
      } else {
         $ligbits->{$classtype}->{$sid}->{cumulative}->Union(
            $ligbits->{$classtype}->{$sid}->{cumulative},
            $ligbits->{$classtype}->{$sid}->{$ligsig}) ;
      }
   }
   close(LIGFH) ;

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      my $newligbits = cluster_ligassignments({
         bs_type => "lig",
         ligbits => $ligbits,
         class2sid => $class2sid,
      }) ;
   }

   return $ligbits ;
}



sub run_cluster_lig {
   require Bit::Vector ;

   my $in = shift ;

   my $pibase_specs = pibase::get_specs() ;
   my $pilig_specs = set_pilig_specs() ;
   my $fn = $pilig_specs->{fn} ;
   my $standardres =  _list_standardres();
   my $liginfo = _pilig_load_liginfo() ;
   my $class2alnlength ;

   my $ligbits = readin_ligassignments({
      fn => $pilig_specs->{outfiles}->{assign_lig},
      liginfo => $liginfo,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
      cluster_fl => 1
   }) ;
}

sub cluster_ligassignments {

   my $in = shift ;
# cluster within family classes
   my $class2sid = $in->{class2sid} ;
   my $ligbits = $in->{ligbits} ;
   my $pilig_specs = set_pilig_specs() ;

   my $clusters;
#test case for implementing speedup: superfam a.132.1: 97 bs clustered to 16.
   foreach my $classtype (keys %{$class2sid}) {
      foreach my $class (sort keys %{$class2sid->{$classtype}}) {
         my $cur_ligs ;
         my ($bsclust, $bs2clust) ;
         foreach my $sid (keys %{$class2sid->{$classtype}->{$class}}) {
            foreach my $lig (keys %{$ligbits->{$classtype}->{$sid}}) {
               if ($lig eq 'cumulative') {next;}
               push @{$cur_ligs}, $sid."\n".$lig ;
               $bs2clust->{$sid."\n".$lig} = $#{$cur_ligs} ;
               $bsclust->{$#{$cur_ligs}} = [$sid."\n".$lig] ;
         } }

         if ($#{$cur_ligs} < 0) {next;}

         my $t_and;
         my $t_or;
         {
            my ($t_sid, $t_lig) = split(/\n/, $cur_ligs->[0]) ;
            $t_and = $ligbits->{$classtype}->{$t_sid}->{$t_lig}->Shadow() ;
            $t_or = $ligbits->{$classtype}->{$t_sid}->{$t_lig}->Shadow() ;
         }

         my $tj = 0 ;
         my @clust_list = sort {$a <=> $b} keys %{$bsclust} ;
         while ($tj <= ($#clust_list - 1)) {
            
            my $clust1 = $clust_list[$tj] ;

            my @merge_list ;
            foreach my $tk (($tj + 1) .. $#clust_list) {
               my $clust2 = $clust_list[$tk] ;
               
               my $mergethis = 0 ;
               foreach my $sidlig1 (@{$bsclust->{$clust1}}) {
                  my ($t_sid1, $t_lig1) = split(/\n/, $sidlig1) ;

                  foreach my $sidlig2 (@{$bsclust->{$clust2}}) {
                     my ($t_sid2, $t_lig2) = split(/\n/, $sidlig2) ;

                     $t_and->And(
                        $ligbits->{$classtype}->{$t_sid1}->{$t_lig1},
                        $ligbits->{$classtype}->{$t_sid2}->{$t_lig2},
                     ) ;
                     $t_or->Or(
                        $ligbits->{$classtype}->{$t_sid1}->{$t_lig1},
                        $ligbits->{$classtype}->{$t_sid2}->{$t_lig2},
                     ) ;
                     my $score = $t_and->Norm() / $t_or->Norm() ;
                     if ($score >= $pilig_specs->{PARAM_MIN_BS_SIMILARITY}) {
                        $mergethis = 1 ;
                        last ;
                     }
                  }
                  if ($mergethis) {last;}
               }
               if ($mergethis) { push @merge_list, $clust1."\n".$clust2; }
            }

            if ($#merge_list >= 0) {
               foreach my $j ( 0 .. $#merge_list) {
                  my ($clust1, $clust2) = split(/\n/, $merge_list[$j]) ;
                  foreach my $oldbs2 (@{$bsclust->{$clust2}}) {
                     $bs2clust->{$oldbs2} = $clust1 ; }

                  push @{$bsclust->{$clust1}}, @{$bsclust->{$clust2}} ;
                  delete $bsclust->{$clust2} ;
               }
               @clust_list = sort {$a <=> $b} keys %{$bsclust} ;
            } else {
               $tj++ ;
            }
         }

         foreach my $clust (keys %{$bsclust}) {
            foreach my $bs (@{$bsclust->{$clust}}) {
               my ($sid, $lig) = split(/\n/, $bs) ;
               $clusters->{$classtype}->{$class}->{$sid}->{$lig} = $clust ;
#               print join("\t", $classtype, $class, $clust, $sid, $lig)."\n" ;
            }
         }
      }
   }

# assume that anytime clusters are generated, should be displayed and
#  saved to the appropriate outfiles place...
   if (!exists $in->{bs_type} || $in->{bs_type} eq 'lig') {
      open(ASSIGNLIG, $pilig_specs->{outfiles}->{assign_lig}) ;
      open(OUTF_LIG, ">".$pilig_specs->{outfiles}->{assign_lig_clusters}) ;
      my $clust_counter = { sf => 1, fam => 1 } ;
      my $clustsseen ;
      while (my $line = <ASSIGNLIG>) {
         chomp $line;
         if ($line =~ /^#/) {print $line."\n" ;next;}
         my @t = split(/\t/, $line) ;

         my ( $pdb, $sid, $osid, $classtype, $class, 
              $btype, $alnposstring, $alnlength,
              $ligcod, $ligid ) = @t ;
         $ligcod =~ s/ //g ;
         my $ligsig = $pdb."\t".$ligcod."\t".$ligid ;
         my $clustnum ;
         if (exists $clusters->{$classtype}->{$class}->{$sid}->{$ligsig}) {
            if (exists
              $clustsseen->{$classtype}->{$class}->{$clusters->{$classtype}->{$class}->{$sid}->{$ligsig}}){
               $clustnum =
               $clustsseen->{$classtype}->{$class}->{$clusters->{$classtype}->{$class}->{$sid}->{$ligsig}};
            } else {
               $clustsseen->{$classtype}->{$class}->{$clusters->{$classtype}->{$class}->{$sid}->{$ligsig}}=
                  $clust_counter->{$classtype} ;
               $clustnum = $clust_counter->{$classtype} ;
               $clust_counter->{$classtype}++ ;
            }
         } else {
            $clustnum = $clust_counter->{$classtype} ;
            $clust_counter->{$classtype}++ ;
         }
         print OUTF_LIG join("\t", $clustnum, @t)."\n" ;
      }
      close(ASSIGNLIG);
      close(OUTF_LIG) ;

   } elsif ($in->{bs_type} eq 'pepnuci') {

      open(ASSIGNPEPNUCI, $pilig_specs->{outfiles}->{assign_pepnuci}) ;
      open(OUTF_PEPNUCI,
         ">".$pilig_specs->{outfiles}->{assign_pepnuci_clusters}) ;
      my $clust_counter = { sf => 1, fam => 1 } ;
      my $clustsseen ;
      while (my $line = <ASSIGNPEPNUCI>) {
         chomp $line;
         if ($line =~ /^#/) {print $line."\n" ;next;}
         my @t = split(/\t/, $line) ;

         my ($pdb, $sid, $osid, $classtype, $class, 
             $targetch_type, $targetch_sid,
             $targetch_len, $alnposstring, $alnlength) = @t ;

         my $clustnum ;
         if (exists $clusters->{$classtype}->{$class}->{$sid}->{$targetch_sid}){
            if (exists
   $clustsseen->{$classtype}->{$class}->{$clusters->{$classtype}->{$class}->{$sid}->{$targetch_sid}}){
               $clustnum =
   $clustsseen->{$classtype}->{$class}->{$clusters->{$classtype}->{$class}->{$sid}->{$targetch_sid}} ;
            } else {
   $clustsseen->{$classtype}->{$class}->{$clusters->{$classtype}->{$class}->{$sid}->{$targetch_sid}} =
                  $clust_counter->{$classtype} ;
               $clustnum = $clust_counter->{$classtype} ;
               $clust_counter->{$classtype}++ ;
            }
         } else {
            $clustnum = $clust_counter->{$classtype} ;
            $clust_counter->{$classtype}++ ;
         }
         print OUTF_PEPNUCI join("\t", $clustnum, $line)."\n" ;
      }
      close(ASSIGNPEPNUCI);
      close(OUTF_PEPNUCI) ;
      
   }

   return $clusters ;

}

sub run_cluster_pep {
   require Bit::Vector ;

   my $in = shift ;

   my $pibase_specs = pibase::get_specs() ;
   my $pilig_specs = set_pilig_specs() ;

   my $standardres =  _list_standardres();
   my $pb = _pilig_tod_pibase_preload() ;
   my $class2alnlength ;

   readin_pepnuciassignments({
      fn => $pilig_specs->{outfiles}->{assign_pepnuci},
      pb => $pb,
      pilig_specs => $pilig_specs,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
      cluster_fl => 1,
   }) ;
}


sub readin_pepnuciassignments {

   my $in = shift ;

   my $fn = $in->{fn} ;

   my $class2alnlength = $in->{class2alnlength};
   my $standardres= $in->{standardres};
   my $pilig_specs = $in->{pilig_specs} ;

# if this is the clustered assignment list, only read in cluster representative
# (arbitrarily chosen - first one in the list)
   my $clustrep_fl = 0;
   my $clusters_seen ;  # keep track of what clusters have been seen already
   if (exists $in->{clustrep_fl} && $in->{clustrep_fl} == 1) {
      $clustrep_fl = 1 ; }

   my $targetch_info ;
   my $pep_bits = {};
   my $nuc_bits = {};
   print STDERR "NOW reading PEPNUCI assignment\n" ;
   open(PEPNUCIF, $fn) ;
   my $class2sid ;
   my $pepnuci_info ;
   while (my $line = <PEPNUCIF>) {
      if ($line =~ /^#/ ||
          $line =~ /^Warning: no access to/ ||
          $line =~ /^Thus no job control/ ) {next;}

      chomp $line;
      my @fields = split(/\t/, $line) ;

      my $cluster_num ;
      if ($clustrep_fl) {
         $cluster_num = shift @fields ; }

      my ( $pdb, $sid, $osid, $classtype, $class, 
           $targetch_type, $targetch_sid,
           $targetch_len, $alnposstring, $alnlength ) = @fields ;
      $class2sid->{$classtype}->{$class}->{$sid}++ ;

      my $t_alnposstring = $alnposstring ;
      $t_alnposstring =~ s/\,//g ; $t_alnposstring =~ s/undef//g ;
      if ($t_alnposstring eq '') {next;}

      if ($clustrep_fl) {
         if (!exists $clusters_seen->{$classtype}->{$cluster_num}) {
            $clusters_seen->{$classtype}->{$cluster_num}++ ;
         } else {
            next ;
         }
      }


      if ($targetch_type eq 'p' &&
          $targetch_len < $pilig_specs->{PARAM_MIN_PEPTIDE_LENGTH}) {
         next;
      }
      my @t_alnpos = split(/\,/, $alnposstring) ;
      if (($#t_alnpos + 1) <
          $pilig_specs->{PARAM_MIN_NUMRES_INTERACTING_WITH_PEPTIDE}) {
         next; }

      $class2alnlength->{$classtype}->{$class} = $alnlength ;
      $targetch_info->{$targetch_sid}->{chain_type} = $targetch_type ;
      $targetch_info->{$targetch_sid}->{chain_length} = $targetch_len ;

      my $bitref = $pep_bits ;
      if ($targetch_type eq 'n') {
         $bitref = $nuc_bits ; }

      $bitref->{$classtype}->{$sid}->{$targetch_sid} =
         Bit::Vector->new($alnlength) ;

      my @alnpos = split(/\,/, $alnposstring) ;
      foreach my $alnpos (@alnpos) {
         if ($alnpos eq 'undef') {next;}
         $bitref->{$classtype}->{$sid}->{$targetch_sid}->Bit_On($alnpos); }


      if (!exists $bitref->{$classtype}->{$sid}->{cumulative}) {
         $bitref->{$classtype}->{$sid}->{cumulative} =
            $bitref->{$classtype}->{$sid}->{$targetch_sid}->Clone();
      } else {
         $bitref->{$classtype}->{$sid}->{cumulative}->Union(
            $bitref->{$classtype}->{$sid}->{cumulative},
            $bitref->{$classtype}->{$sid}->{$targetch_sid}) ;
      }
   }
   close(PEPNUCIF) ;

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      my $clustered_bits = cluster_ligassignments({
         bs_type => "pepnuci",
         ligbits => $pep_bits,
         class2sid => $class2sid,
      }) ;
   }

   return {
      pep => $pep_bits,
      nuc => $nuc_bits,
      chain_info => $targetch_info
   } ;

}


sub readin_piassignments {

   my $in = shift;
   my $fn = $in->{fn} ;
   my $class2alnlength = $in->{class2alnlength};
   my $pb = $in->{pb} ;
   my $pibits ;
   my $interfaces ;

# if this is the clustered assignment list, only read in cluster representative
# (arbitrarily chosen - first one in the list)
   my $clustrep_fl = 0;
   my $clusters_seen ;  # keep track of what clusters have been seen already
   if (exists $in->{clustrep_fl} && $in->{clustrep_fl} == 1) {
      $clustrep_fl = 1 ; }

   print STDERR "NOW reading PI assignment\n" ;
   my $classes2int ;
   my $class2chains2intside;
   open(PIFH, $fn) ;
   while (my $line = <PIFH>) {
      if ($line =~ /^#/ || $line =~ /^Warning: no access to/ ||
         $line =~ /^Thus no job control/ ) {next;}

      chomp $line;

      my @t = split(/\t/, $line) ;
      my $cluster_num ;
      if ($clustrep_fl) {
         $cluster_num = shift @t ; }

      my ($pdb, $sid, $osid, $classtype, $class, $obtype, $alnposstring,
          $alnlength, $sid1, $sid2, $fam1, $fam2, $chains) = @t ;

#if flag specified, don't read in data for j.* peptide domains
      if (exists $in->{dont_read_jdomains_fl} &&
          ($pb->{sid2class}->{fam}->{$sid1} =~ /^j/ ||
           $pb->{sid2class}->{fam}->{$sid2} =~ /^j/)) {
         next;
      }

      my @alnpos = split(/\,/, $alnposstring) ;
      {
         my @t = ();
         foreach my $p ( @alnpos) {
            if ($p ne 'undef') {push @t, $p;} }
         @alnpos = @t ;
      }
      if ($#alnpos < 0) {next;}

      my $sid12 = $sid1."\t".$sid2 ;
      my $side = 1; if ($sid eq $sid2) {$side = 2;}

      if ($clustrep_fl) {
         if ( !exists $clusters_seen->{$classtype}->{$cluster_num} ||
              $clusters_seen->{$classtype}->{$cluster_num} eq $sid12 ) {
            $clusters_seen->{$classtype}->{$cluster_num} = $sid12 ;
# have to keep track of what was chosen, since need to read in both sides
# of the interface...
         } else {
            next;
         }
      }

      $class2chains2intside->{$classtype}->{$class}->{$chains}->{$sid12."\n".$side}++ ;
      {
         my $temp_class1 = $pb->{sid2class}->{fam}->{$sid1} ;
         my $temp_class2 = $pb->{sid2class}->{fam}->{$sid2} ;
         my $sf1 = $pb->{sid2class}->{sf}->{$sid1} ;
         my $sf2 = $pb->{sid2class}->{sf}->{$sid2} ;
         my $temp_revfl = 0 ;
         if ($temp_class2 lt $temp_class1) {
            $temp_class1 = $pb->{sid2class}->{fam}->{$sid2} ;
            $temp_class2 = $pb->{sid2class}->{fam}->{$sid1} ;
            $sf1 = $pb->{sid2class}->{sf}->{$sid2} ;
            $sf2 = $pb->{sid2class}->{sf}->{$sid1} ;
            $temp_revfl = 1 ;
         }
         $classes2int->{sf}->{$sf1."\t".$sf2}->{$sid12} = $temp_revfl ;
         $classes2int->{fam}->{$temp_class1."\t".$temp_class2}->{$sid12} =
            $temp_revfl ;
      }

      $interfaces->{$sid12}->{$side}->{sid} = $sid ;
      $interfaces->{$sid12}->{chains} = $chains ;
      $interfaces->{$sid12}->{pdb} = $pdb;

      $class2alnlength->{$classtype}->{$class} = $alnlength ;

#      my @btypes = ($obtype) ;
# fpd090104_1840  - changed this so that when P bits are set in
#                   combine_piligpepexpbits() it includes intra,inter,and peptide
      my $btype = 'Pinter';
      if ($chains eq 'same') {
         $btype = "Pintra" ; }

      $interfaces->{$sid12}->{$side}->{pibits}->{$classtype} =
         Bit::Vector->new($alnlength) ;

      foreach my $alnpos (@alnpos)  {
         if ($alnpos eq 'undef') {next;}
         $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Bit_On($alnpos);}

      if (!exists $pibits->{$classtype}->{$sid}->{$btype}) {
         $pibits->{$classtype}->{$sid}->{$btype} =
            $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Clone();
      } else {
         $pibits->{$classtype}->{$sid}->{$btype}->Union(
            $pibits->{$classtype}->{$sid}->{$btype},
            $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}) ;
      }
   }
   close(PIFH) ;

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      my $newligbits = cluster_piassignments({
         classes2int => $classes2int,
         pibits => $pibits,
         interfaces => $interfaces,
         class2chains2intside => $class2chains2intside,
      }) ;
   }


   return {
      pibits => $pibits,
      interfaces => $interfaces
   } ;

}

sub run_cluster_pi {
   require Bit::Vector ;

   my $in = shift ;
   my $pibase_specs = pibase::get_specs() ;
   my $pilig_specs = set_pilig_specs() ;

   my $standardres =  _list_standardres();
   my $pb = _pilig_tod_pibase_preload() ;
   my $class2alnlength ;

   readin_piassignments({
      fn => $pilig_specs->{outfiles}->{assign_pi},
      pb => $pb,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
      cluster_fl => 1,
   }) ;

}

sub cluster_piassignments {

   my $in = shift ;
# cluster within family classes
   my $class2chains2intside = $in->{class2chains2intside} ;

   my $interfaces = $in->{interfaces} ;
   my $pilig_specs = set_pilig_specs() ;
   my $classes2int = $in->{classes2int} ;

   my $newligbits ;

#cluster binding sites
   my $glob_bsclust ;
   foreach my $classtype (sort keys %{$class2chains2intside}) {
      foreach my $class (sort keys %{$class2chains2intside->{$classtype}}) {
       foreach my $chains (sort
       keys %{$class2chains2intside->{$classtype}->{$class}}){
         my $cur_bs;
         my ($bsclust, $bs2clust) ;
         foreach my $sid12side (
          keys %{$class2chains2intside->{$classtype}->{$class}->{$chains}}) {
            push @{$cur_bs}, $sid12side ;
            $bs2clust->{$sid12side} = $#{$cur_bs} ;
            $bsclust->{$#{$cur_bs}} = [$sid12side] ;
         }

         if ($#{$cur_bs} < 0) {next;}

         my $t_and;
         my $t_or;
         {
            my ($t_sid12, $t_side) = split(/\n/, $cur_bs->[0]) ;
            $t_and = $interfaces->{$t_sid12}->{$t_side}->{pibits}->{$classtype}->Shadow();
            $t_or = $interfaces->{$t_sid12}->{$t_side}->{pibits}->{$classtype}->Shadow() ;
         }

         my $tj = 0 ;
         my @clust_list = sort {$a <=> $b} keys %{$bsclust} ;
         while ($tj <= ($#clust_list - 1)) {
            my $clust1 = $clust_list[$tj] ;

            my @merge_list ;
            foreach my $tk (($tj + 1) .. $#clust_list) {
               my $clust2 = $clust_list[$tk] ;
               
               my $mergethis = 0 ;
               foreach my $sid12side1 (@{$bsclust->{$clust1}}) {
                  my ($t_sidpair_1, $t_side_1) = split(/\n/, $sid12side1) ;

                  foreach my $sid12side2 (@{$bsclust->{$clust2}}) {
                     my ($t_sidpair_2, $t_side_2) = split(/\n/, $sid12side2) ;

                     $t_and->And(
         $interfaces->{$t_sidpair_1}->{$t_side_1}->{pibits}->{$classtype},
         $interfaces->{$t_sidpair_2}->{$t_side_2}->{pibits}->{$classtype}
                     ) ;
                     $t_or->Or(
         $interfaces->{$t_sidpair_1}->{$t_side_1}->{pibits}->{$classtype},
         $interfaces->{$t_sidpair_2}->{$t_side_2}->{pibits}->{$classtype}
                     ) ;
                     my $score = $t_and->Norm() / $t_or->Norm() ;
                     if ($score >= $pilig_specs->{PARAM_MIN_BS_SIMILARITY}) {
                        $mergethis = 1 ;
                        last ;
                     }
                  }
                  if ($mergethis) {last;}
               }
               if ($mergethis) { push @merge_list, $clust1."\n".$clust2; }
            }

            if ($#merge_list >= 0) {
               foreach my $j ( 0 .. $#merge_list) {
                  my ($clust1, $clust2) = split(/\n/, $merge_list[$j]) ;
                  foreach my $oldbs2 (@{$bsclust->{$clust2}}) {
                     $bs2clust->{$oldbs2} = $clust1 ; }

                  push @{$bsclust->{$clust1}}, @{$bsclust->{$clust2}} ;
                  delete $bsclust->{$clust2} ;
               }
               @clust_list = sort {$a <=> $b} keys %{$bsclust} ;
            } else {
               $tj++ ;
            }
         }

         foreach my $clust (keys %{$bsclust}) {
            foreach my $bs (@{$bsclust->{$clust}}) {
               my ($sid12, $side) = split(/\n/, $bs) ;
               my ($sid1, $sid2) = split(/\t/, $sid12) ;

               $glob_bsclust->{$classtype}->{$sid12}->{$side} =
                  $class."\t".$chains."\t".$clust ;

#               print join("\t", "BS", $classtype, $class,
#                          $clust, $sid1, $sid2, $side)."\n" ;
            }
         }

       }
      }
   }

# cluster interfaces
   my $intclust ;
   foreach my $classtype (keys %{$classes2int}) {
      foreach my $classpair (keys %{$classes2int->{$classtype}}) {
         foreach my $sid12 (keys %{$classes2int->{$classtype}->{$classpair}}) {
            my $t1 = $glob_bsclust->{$classtype}->{$sid12}->{1} ;
            my $t2 = $glob_bsclust->{$classtype}->{$sid12}->{2} ;

            if (!defined $t1) {$t1 = "NOCLUST\t$sid12\t1" ;}
            if (!defined $t2) {$t2 = "NOCLUST\t$sid12\t2" ;}

            my $intclustype = $t1."\t".$t2 ;
            my $revfl = 0 ;
            if ($t2 lt $t1) {
               $revfl = 1 ;
               $intclustype = $t2."\t".$t1 ;
            }
            $intclust->{$classtype}->{$intclustype}->{$sid12} = $revfl ;
         }
      }
   }

   my $cluster_ass ;
   foreach my $classtype (keys %{$intclust}) {
      my $intclust_num = 0 ;
      foreach my $intclust_type (keys %{$intclust->{$classtype}}) {
         foreach my $sid12 (keys %{$intclust->{$classtype}->{$intclust_type}}) {
            $cluster_ass->{$classtype}->{$sid12} = $intclust_type ;
               my ($sid1, $sid2) = split(/\t/, $sid12) ;
               my $revfl = $intclust->{$classtype}->{$intclust_type}->{$sid12} ;
#               print join("\t", "INT", $classtype, $intclust_num,
#                          $intclust_type, $revfl, $sid1, $sid2)."\n" ;
         }
         $intclust_num++ ;
      }
   }

   open(ASSIGNPI, $pilig_specs->{outfiles}->{assign_pi}) ;
   open(OUTF_PI, ">".$pilig_specs->{outfiles}->{assign_pi_clusters}) ;
   my $clustsseen ;
   my $clust_counter = { sf => 1, fam => 1} ;
   while (my $line = <ASSIGNPI>) {
      chomp $line;
      if ($line =~ /^#/) {print $line."\n" ;next;}
      my @t = split(/\t/, $line) ;

      my ($pdb, $sid, $osid, $classtype, $class,
          $obtype, $alnposstring, $alnlength,
          $sid1, $sid2, $fam1, $fam2, $chains) = @t ;
      
      my $sid12 = $sid1."\t".$sid2 ;
      my $clustnum ;
      if (exists $cluster_ass->{$classtype}->{$sid12}) {
         if (exists $clustsseen->{$classtype}->{$cluster_ass->{$classtype}->{$sid12}}) {
            $clustnum = $clustsseen->{$classtype}->{$cluster_ass->{$classtype}->{$sid12}} ;
         } else {
            $clustsseen->{$classtype}->{$cluster_ass->{$classtype}->{$sid12}} =
               $clust_counter->{$classtype} ;
            $clustnum = $clust_counter->{$classtype} ;
            $clust_counter->{$classtype}++ ;
         }
      } else {
         $clustnum = $clust_counter->{$classtype} ;
         $clust_counter->{$classtype}++ ;
      }

      print OUTF_PI join("\t", $clustnum, $line)."\n" ;
   }
   close(ASSIGNPI) ;
   close(OUTF_PI) ;

}


sub readin_expassignments {

   my $in = shift ;
   my $fn = $in->{fn};
   my $class2alnlength = $in->{class2alnlength};

   print STDERR "NOW reading EXP assignment\n" ;
   open(EXPFH, $fn) ;
   my $expbits ;
   while (my $line = <EXPFH>) {
      if ($line =~ /^#/ ||
            $line =~ /^Warning: no access to/ ||
            $line =~ /^Thus no job control/ ) {next;}
      chomp $line;

      my ( $pdb, $sid, $osid, $classtype,
           $class, $btype, $alnposstring, $alnlength) = split(/\t/, $line) ;

      $class2alnlength->{$classtype}->{$class} = $alnlength ;

      if (!exists $expbits->{$classtype}->{$class}) {
         $expbits->{$classtype}->{$class} =
            Bit::Vector->new($alnlength) ;
      }

      my @alnpos = split(/\,/, $alnposstring) ;
      foreach my $alnpos (@alnpos)  {
         if ($alnpos eq "undef") {next;}
         $expbits->{$classtype}->{$class}->Bit_On($alnpos);
      }
   }
   close(EXPFH) ;

   return $expbits ;

}

sub _list_standardres {

   return {
   ALA => 'A' ,
   ARG => 'R' ,
   ASN => 'N' ,
   ASP => 'D' ,
   CYS => 'C' ,
   GLN => 'Q' ,
   GLU => 'E' ,
   GLY => 'G' ,
   HIS => 'H' ,
   HSD => 'H' ,
   HSE => 'H' ,
   ILE => 'I' ,
   LEU => 'L' ,
   LYS => 'K' ,
   MET => 'M' ,
   PHE => 'F' ,
   PRO => 'P' ,
   SER => 'S' ,
   THR => 'T' ,
   TRP => 'W' ,
   TYR => 'Y' ,
   VAL => 'V',
   UNK => 'X',
   UNX => 'X',
   '  C' => 'c',
   '  G' => 'g',
   '  A' => 'a',
   '  T' => 't',
   '  U' => 'u',
   '  I' => 'i',
   'C' => 'c',
   'G' => 'g',
   'A' => 'a',
   'T' => 't',
   'U' => 'u',
   'I' => 'i',
   '+C' => 'c',
   '+G' => 'g',
   '+A' => 'a',
   '+T' => 't',
   '+U' => 'u',
   '+I' => 'i'
   } ;

}

sub _list_standardres_20aa {

   return {
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

}


sub _pilig_run_collate_perinstance {

   require Bit::Vector ;

   my $standardres = _list_standardres() ;

# Per-PPI
# 0. load all of ligands first - lig info
# - iterate over assignments, and have individual
#   bit vectors set for the ligand binding sites in their
#   respective domain families
#
# per-PPI: iterate over each domain--domain interface
#  1. list any ligands actually bound to this interface
#  2. iterate over all ligands in that family/superfamily, find
#     - quantify per ligand coverage
#     - quantify cumulative ligand coverage
#
# per-ligand:


   my $pilig_specs = set_pilig_specs() ;
#   my $fn = set_locations() ;
   my $fn = $pilig_specs->{fn} ;
   my $astral = astral_preload({ fn => $fn }) ;
   my $pb = _pilig_tod_pibase_preload() ;
   my $liginfo = _pilig_load_liginfo() ;

   my $class2alnlength = {};

   my $ligbits = readin_ligassignments({
      fn => $fn,
      liginfo => $liginfo,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;

   my $expbits = readin_expassignments({
      fn => $fn,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;

   my $pibits_both = readin_piassignments({
      fn => $fn,
      pb => $pb,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;
   my $pibits = $pibits_both->{pibits} ;
   my $interfaces = $pibits_both->{interfaces} ;

   my $alltogether = combine_piligexpbits({
      pb => $pb,
      liginfo => $liginfo,
      ligbits => $ligbits,
      pibits => $pibits,
      expbits => $expbits,
   });
   my $class2bits = $alltogether->{class2bits} ;
   my $class2ligs = $alltogether->{class2ligs} ;
   my $class2anligbits = $alltogether->{class2anligbits} ;


   my @headers_sum=  ("PINT_LSUM", "PDB", "SID1", "SID2", "CHAINS",
                       "CLASSTYPE", "CLASS1", "CLASS2",
                       "numres_p_1", "cumulative_numres_l_and_p_1",
                       "max_liginrange_1",
                       "max_l_and_p_1", "lig_max_l_and_p_1",
                       "max_obi_1", "lig_max_obi_1",
                       "max_opi_1", "lig_max_opi_1",
                       "max_olig_1", "lig_max_olig_1",
                       "numres_p_2", "cumulative_numres_l_and_p_2",
                       "max_liginrange_2",
                       "max_l_and_p_2", "lig_max_l_and_p_2",
                       "max_obi_2", "lig_max_obi_2",
                       "max_opi_2", "lig_max_opi_2",
                       "max_olig_2", "lig_max_olig_2",
                       ) ;
   print '#'.join("\t",@headers_sum)."\n" ;

#   foreach my $classtype (qw/sf fam/)
   foreach my $classtype (qw/fam/) {
      foreach my $sid12 (keys %{$interfaces}) {
         my $sid ;
         ($sid->{1}, $sid->{2}) = split(/\t/, $sid12) ;
         my $pdb = $interfaces->{$sid12}->{pdb} ;
         my $chains = $interfaces->{$sid12}->{chains} ;
   
         my ($alnlength, $class) ;
         foreach my $side (1, 2) {
            $class->{$side} = $pb->{sid2class}->{$classtype}->{$sid->{$side}} ;
            if (!exists $class2alnlength->{$classtype}->{$class->{$side}}) {
               next ; }
            $alnlength->{$side} = $class2alnlength->{$classtype}->{$class->{$side}};
#            $surf_pos->{$side} = $class2bits->{$classtype}->{$class->{$side}}->{'E'}->Norm() ;
         }
   
         my $bs_bits;
         my $curligs ;
   
         my $alnpos ;
         my $pistats ;
         my $ligstats ; 
         my $skipside ;
         foreach my $s (1, 2) {
            if( !exists $interfaces->{$sid12}->{$s} ||
               !exists $interfaces->{$sid12}->{$s}->{pibits}->{$classtype} ||
               ($interfaces->{$sid12}->{$s}->{pibits}->{$classtype}->Norm() == 0)){
   
               $skipside->{$s}++ ;
               print STDERR "WARNING: skipping $sid->{1} -- $sid->{2}".
                  " undefined binding site alignment positions for $sid->{$s}\n";
               next;
            }
   
            if (exists $class2anligbits->{$classtype}->{$class->{$s}}) {
               $curligs->{$s} =
                  $class2anligbits->{$classtype}->{$class->{$s}} ; }
            $bs_bits->{$s} = $interfaces->{$sid12}->{$s}->{pibits}->{$classtype} ;
         }
   
   # init vals
         foreach my $s (1, 2) {
            foreach my $type (qw/p cumlig_l cumlig_l_and_p max_liginrange max_l_and_p max_obi max_opi max_olig/) {
               $ligstats->{$s}->{$type} = 0 ; }
            foreach my $type (qw/lig_max_l_and_p lig_max_obi lig_max_opi lig_maxolig/) {
               $ligstats->{$s}->{$type} = 'U'; }
         }
   
         foreach my $s (1, 2) {
            if (exists $skipside->{$s}) {next;}
   
            $pistats->{$s}->{p} = $bs_bits->{$s}->Norm() ;
            $pistats->{$s}->{max_liginrange} = 0 ;
            $ligstats->{$s}->{cumlig} = Bit::Vector->new($alnlength->{$s}) ;
            foreach my $t (qw/obi opi olig l_and_p/) {
               $ligstats->{$s}->{"max_$t"} = 0 ;
               $ligstats->{$s}->{"lig_max_$t"} = 'undef' ;
            }
   
            my $temp_bits_and = Bit::Vector->new($alnlength->{$s}) ;
            my $temp_bits_or = Bit::Vector->new($alnlength->{$s}) ;
   
            if (!exists $curligs->{$s}) {next;}
            foreach my $j (0 .. $#{$curligs->{$s}}) {
               my $curligsig = $curligs->{$s}->[$j]->[0] ;
               my $outcurligsig = $curligsig; $outcurligsig =~ s/\t/:/g ;
               my $curligbits = $curligs->{$s}->[$j]->[1] ;
               
               $ligstats->{$s}->{cumlig}->Or($ligstats->{$s}->{cumlig},
                                             $curligbits) ;
   
               #intersect the ligands bit vector and the bit vector for the bindin
               # site positions of this particular interface;
               $temp_bits_and->And($curligbits, $bs_bits->{$s}) ;
               $temp_bits_or->Or($curligbits, $bs_bits->{$s}) ;
   
               my $curs ;
               $curs->{l_and_p} = $temp_bits_and->Norm() ;
               $curs->{l_or_p} = $temp_bits_or->Norm() ;
               $curs->{l} = $curligbits->Norm() ;
               $curs->{p} = $bs_bits->{$s}->Norm() ;
   
               if ($curs->{l_and_p} == 0) {next;}
   
               my $curoverlap ;
               $curoverlap->{obi} = $curs->{l_and_p} / $curs->{l_or_p} ;
               $curoverlap->{opi} = $curs->{l_and_p} / $curs->{p} ;
               $curoverlap->{olig} = $curs->{l_and_p} / $curs->{l} ;
   
   
               if (!exists $ligstats->{$s}->{"max_l_and_p"} ||
                   $curs->{l_and_p} > $ligstats->{$s}->{"max_l_and_p"}) {
                  $ligstats->{$s}->{"max_l_and_p"} = $curs->{l_and_p} ;
                  $ligstats->{$s}->{"lig_max_l_and_p"} = $outcurligsig ;
               }
   
               foreach my $t (qw/obi opi olig/) {
                  if (!exists $ligstats->{$s}->{"max_$t"} ||
                         $curoverlap->{$t} > $ligstats->{$s}->{"max_$t"}) {
                     $ligstats->{$s}->{"max_$t"} = $curoverlap->{$t} ;
                     $ligstats->{$s}->{"lig_max_$t"} = $outcurligsig ;
                  }
               }
   
               my (undef, $ligcod, undef) = split(/\t/, $curligsig) ;
               my $liginrange = 0 ;
               if (exists $liginfo->{mwinrange}->{$ligcod}) {
                   $ligstats->{$s}->{max_liginrange} = 1 ;
                   $liginrange = 1 ; }
   
               my @outvals = ("PINT_LINT",
                                 $pdb, $sid->{1}, $sid->{2},
                                 $chains, $classtype,
                                 $class->{1}, $class->{2},
                                 $s, $sid->{$s}, $outcurligsig,
                                 $liginrange,
                                 $curs->{p}, $curs->{l},
                                 $curs->{l_and_p}, $curs->{l_or_p},
                                 sprintf("%.3f", $curoverlap->{obi}),
                                 sprintf("%.3f", $curoverlap->{opi}),
                                 sprintf("%.3f", $curoverlap->{olig})
                              ) ;
               print join("\t", @outvals)."\n";
            }
   
            $ligstats->{$s}->{"cumlig_l"} = $ligstats->{$s}->{cumlig}->Norm() ;
            $temp_bits_and->And($ligstats->{$s}->{cumlig}, $bs_bits->{$s}) ;
            $ligstats->{$s}->{"cumlig_l_and_p"} = $temp_bits_and->Norm() ;
         }
   
         my @outvals =  ("PINT_LSUM", $pdb, $sid->{1}, $sid->{2}, $chains,
                                 $classtype, $class->{1}, $class->{2}) ;
   
         foreach my $s ( 1, 2) {
            if (exists $skipside->{$s} ||
               !exists $ligstats->{$s} ||
               $ligstats->{$s}->{"cumlig_l"} == 0 ) {
   
               if (exists $skipside->{$s}) {
                  push @outvals, 0, 0, 0 ;
               } else {
                  push @outvals, $pistats->{$s}->{p}, 0, 0 ;
               }
   
               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
               next;
            }
   
            push @outvals, $pistats->{$s}->{p} ;
            push @outvals, $ligstats->{$s}->{cumlig_l_and_p} ;
            push @outvals, $ligstats->{$s}->{max_liginrange} ;
            push @outvals, $ligstats->{$s}->{"max_l_and_p"} ;
            push @outvals, $ligstats->{$s}->{"lig_max_l_and_p"} ;
            foreach my $t (qw/obi opi olig/)  {
               push @outvals, sprintf("%.3f", $ligstats->{$s}->{"max_$t"}) ;
               push @outvals, $ligstats->{$s}->{"lig_max_$t"} ;
            }
         }
         print join("\t", @outvals)."\n"; 
      }
   }

}


sub collate_perfam {

   require Bit::Vector ;
   require R;
   require RReferences ;
   R::initR("--silent") ;

   my $standardres =  _list_standardres();

# given deposit dirs,
# load list of all ASTEROIDS fam/sf

   my $pibase_specs = pibase::get_specs() ;
   my $pilig_specs = set_pilig_specs() ;
   my $fn = $pilig_specs->{fn} ;

   my $astral = _pilig_astral_preload({
      pibase_specs => $pibase_specs,
      pilig_specs => $pilig_specs,
   }) ;
   my $pb = _pilig_tod_pibase_preload() ;
   my $liginfo = _pilig_load_liginfo() ;

   my $class2alnlength = {};

   my ($pibits_both, $pibits, $interfaces,
       $pepnucibits, $ligbits) ;

#      $pibits_both = readin_piassignments({
#         fn => $pilig_specs->{outfiles}->{assign_pi_clusters},
#         clustrep_fl => 1,
#         pb => $pb,
#         standardres => $standardres,
#         class2alnlength => $class2alnlength,
#      }) ;
#
#      $pibits = $pibits_both->{pibits} ;
#      $interfaces = $pibits_both->{interfaces} ;
#
#      $pepnucibits = readin_pepnuciassignments({
#         fn => $pilig_specs->{outfiles}->{assign_pepnuci_clusters},
#         clustrep_fl => 1,
#         pb => $pb,
#         pilig_specs => $pilig_specs,
#         standardres => $standardres,
#         class2alnlength => $class2alnlength,
#      }) ;
#
#      $ligbits = readin_ligassignments({
#         fn => $pilig_specs->{outfiles}->{assign_lig_clusters},
#         clustrep_fl => 1,
#         liginfo => $liginfo,
#         standardres => $standardres,
#         class2alnlength => $class2alnlength,
#      }) ;

      $pibits_both = readin_piassignments({
         fn => $pilig_specs->{outfiles}->{assign_pi_clusters},
         clustrep_fl => 1,
         pb => $pb,
         standardres => $standardres,
         class2alnlength => $class2alnlength,
         dont_read_jdomains_fl => 1,
      }) ;

      $pibits = $pibits_both->{pibits} ;
      $interfaces = $pibits_both->{interfaces} ;

      $pepnucibits = readin_pepnuciassignments({
         fn => $pilig_specs->{outfiles}->{assign_pepnuci_clusters},
         clustrep_fl => 1,
         pb => $pb,
         pilig_specs => $pilig_specs,
         standardres => $standardres,
         class2alnlength => $class2alnlength,
      }) ;

      $ligbits = readin_ligassignments({
         fn => $pilig_specs->{outfiles}->{assign_lig_clusters},
         clustrep_fl => 1,
         liginfo => $liginfo,
         standardres => $standardres,
         class2alnlength => $class2alnlength,
      }) ;

   my $expbits = readin_expassignments({
      fn => $pilig_specs->{outfiles}->{assign_exp},
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;


   my $alltogether = combine_piligpepexpbits({
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

   {
      my @headers ;
      push @headers, qw/classtype class alnlength/ ;
      my @tt = qw/L P Pinter Pintra p E/ ;
      foreach my $t (@tt) {
         push @headers, "cons_shanonne_std20_$t" ;
         push @headers, "cons_numtypes_std20_$t" ;
         push @headers, "bits_$t" ;
         push @headers, "norm_$t" ;
      }

# don't need E fields.
      foreach my $t (@tt[1..($#tt - 1)]) {
         push @headers, "cons_shannone_std20_".$tt[0]."_and_$t" ;
         push @headers, "cons_numtypes_std20_".$tt[0]."_and_$t" ;
         push @headers, "bits_".$tt[0]."_and_$t" ;
         push @headers, "norm_".$tt[0]."_and_$t" ;
         push @headers, "bits_".$tt[0]."_or_$t" ;
         push @headers, "norm_".$tt[0]."_or_$t" ;
         push @headers, "pval_".$tt[0]."_$t" ;
         push @headers, "pval_".$tt[0]."_$t".'_less' ;
      }
      push @headers, 'ligcodes' ;

      foreach my $j (0 .. $#headers) {
         print '#'.($j + 1)."\t".$headers[$j]."\n" ; }
      print '#'.join("\t", @headers)."\n" ;
   }

   my $astral_classes = $astral->{classes} ;
#   foreach my $classtype (qw/seqcl90 seqcl95 seqcl100 fam sf/)
   if ($pilig_specs->{CALC_CONSSCORE_FLAG}) {
      open(OUTCONSSCORE, ">>outconsscore.$$.out") ;}
   foreach my $classtype (qw/fam/) {
      my @classes ;
      my $curseqidlevel ;
      if ($classtype eq 'fam' || $classtype eq 'sf') {
         @classes = sort keys %{$astral_classes->{$classtype}} ;
      } elsif ($classtype =~ /^seqcl/) {
         my ($seqid) = ($classtype =~ /seqcl([0-9]+)/) ;
         $curseqidlevel = $seqid ;
         @classes = sort keys %{$astral->{seqcl2cont}->{$seqid}} ;
      }

      foreach my $class (@classes) {
         my $allprocalns ;
         my $outbits ;
         my $outnorm ;

         if (!exists $class2bits->{$classtype}->{$class}->{L} &&
             !exists $class2bits->{$classtype}->{$class}->{P}) {
            my @outvals = ( $classtype, $class, 'NOINTERACTIONS' );
            print join("\t", @outvals)."\n" ;
            next;
         }

         print STDERR "NOW ON $classtype $class (line ".__LINE__.")\n" ;
         my $class_aln ;
         if ($pilig_specs->{CALC_CONSSCORE_FLAG}) {
            my ($curfam, $curclasstype) ;
            if ($classtype eq 'fam' || $classtype eq 'sf') {
               $curfam = $class;
#BUG070503_1511                $curclasstype = $class ;
               $curclasstype = $classtype ;
            } elsif ($classtype =~ /^seqcl/) {
               $curfam = $pb->{osid2class}->{fam}->{$class} ;
               $curclasstype = 'fam' ;
            }

            if (! exists $allprocalns->{$curclasstype}->{$curfam}) {
#               print STDERR "OHFUCK $curclasstype $curfam aln file is ".
#                  $fn->{aln}->{$curclasstype}.'/'.$curfam.'.aln.fa'."\n" ;
               $allprocalns->{$curclasstype}->{$curfam} =
                  pibase::ASTRAL::load_asteroids_aln({
                  aln_fn => $pibase_specs->{asteroids}->{$curclasstype.'_aln'}.
                  '/'.$curfam.'.fasta_aln' ,
                  seq_fn => $pibase_specs->{asteroids}->{$curclasstype.'_seq'}.
                  '/'.$curfam.'.fa' ,
                  raf => $astral->{raf},
                  gdseqh => $astral->{gdseqh},
                  seqclcont100 => $astral->{seqcl2cont}->{100},
                  seqcl100 => $astral->{seqcl}->{100},
                  allchains => $pb->{pdbchains}
               }) ;
            }

            $class_aln = $allprocalns->{$curclasstype}->{$curfam} ;
         }

         my $consscore_se ;
         my $consscore_nt ;

         foreach my $btype (qw/L P Pinter Pintra p E/) {
            if (exists $class2bits->{$classtype}->{$class}->{$btype}) {
               $outnorm->{$btype} = 
               $class2bits->{$classtype}->{$class}->{$btype}->Norm();

               $outbits->{$btype} =
               $class2bits->{$classtype}->{$class}->{$btype}->to_Bin();

#               if ($classtype eq 'fam' || $classtype eq 'sf')
               if ($pilig_specs->{CALC_CONSSCORE_FLAG} && 
                   ($classtype eq 'fam' || $classtype eq 'sf')) {
                  my $tsum_se = 0;
                  my $tsum_nt = 0;
                  foreach my $tpos (
   $class2bits->{$classtype}->{$class}->{$btype}->Index_List_Read()) {


                     my $resfreq_string = '';
                     foreach my $aatype ( sort
                        keys %{$class_aln->{meta}->{res_freq}->{$tpos}}) {
                        $resfreq_string .= $aatype.
                           $class_aln->{meta}->{res_freq}->{$tpos}->{$aatype} ;
                     }
                     print OUTCONSSCORE join("\t",$class, $classtype,$btype,
                        $tpos,
                        $class_aln->{meta}->{shannone_std20}->{$tpos},
                        $class_aln->{meta}->{numrestypes_std20}->{$tpos},
                        $resfreq_string
                     )."\n" ;

                     $tsum_se += $class_aln->{meta}->{shannone_std20}->{$tpos};
                     $tsum_nt+=$class_aln->{meta}->{numrestypes_std20}->{$tpos};
                  }
                  if ($outnorm->{$btype} > 0) {
                     $consscore_se->{$btype} = $tsum_se / $outnorm->{$btype} ;
                     $consscore_nt->{$btype} = $tsum_nt / $outnorm->{$btype} ;
                  } else {
                     $consscore_se->{$btype} = 'undef' ;
                     $consscore_nt->{$btype} = 'undef' ;
                  }
               } else {
                  $consscore_se->{$btype} = "undefseqcl" ;
                  $consscore_nt->{$btype} = "undefseqcl" ;
               }

            } else {
               my $t = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
               $consscore_se->{$btype} = "undef" ;
               $consscore_nt->{$btype} = "undef" ;
               $outbits->{$btype} = $t->to_Bin() ;
               $outnorm->{$btype} = 0 ;
            }
         }

         my $overlap_pvals ;
         foreach my $ptype (qw/P Pinter Pintra p/) {
            $overlap_pvals->{"L_$ptype"} = "undef" ;
            $overlap_pvals->{"L_$ptype"."_less"} = "undef" ;
            if ($outnorm->{L} > 0 && $outnorm->{$ptype} > 0) {
               my $tor = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
               $tor->Or(
                  $class2bits->{$classtype}->{$class}->{$ptype},
                  $class2bits->{$classtype}->{$class}->{L}) ;

               $outbits->{"L_or_$ptype"} = $tor->to_Bin() ;
               $outnorm->{"L_or_$ptype"} = $tor->Norm() ;

               my $tand = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
               $tand->And(
                     $class2bits->{$classtype}->{$class}->{$ptype},
                     $class2bits->{$classtype}->{$class}->{L}) ;

               $outbits->{"L_and_$ptype"} = $tand->to_Bin() ;
               $outnorm->{"L_and_$ptype"} = $tand->Norm() ;
               $consscore_se->{"L_and_$ptype"} = "undef" ;
               $consscore_nt->{"L_and_$ptype"} = "undef" ;

               my $tsum_se = 0;
               my $tsum_nt = 0;
               $overlap_pvals->{"L_$ptype"} = 'pvalnotcalc' ;
               $overlap_pvals->{"L_$ptype"."_less"} = 'pvalnotcalc' ;

               if ($pilig_specs->{CALC_PVAL_FLAG}) {
                  my $tnorm ;
                  $tnorm->{lp} = $outnorm->{"E"} - $outnorm->{"L_or_$ptype"} ;
                  $tnorm->{Lp} = $outnorm->{"L"} - $outnorm->{"L_and_$ptype"} ;
                  $tnorm->{lP} = $outnorm->{$ptype} -
                                 $outnorm->{"L_and_$ptype"} ;
                  $tnorm->{LP} = $outnorm->{"L_and_$ptype"} ;

                  my @x = &R::eval("capture.output(fisher.test(matrix(c($tnorm->{lp}, $tnorm->{Lp}, $tnorm->{lP}, $tnorm->{LP}), 2, 2),workspace=2e7,alternative='g'))") ;
                  my $x = join("\n", @x) ;
                  my ($pval) = ($x =~ /p-value [\<=] (.*)\nalternative/) ;
                  $overlap_pvals->{"L_$ptype"} = $pval ;

                  my @y = &R::eval("capture.output(fisher.test(matrix(c($tnorm->{lp}, $tnorm->{Lp}, $tnorm->{lP}, $tnorm->{LP}), 2, 2),workspace=2e7,alternative='l'))") ;
                  my $y = join("\n", @y) ;
                  my ($pval_less) = ($y =~ /p-value [\<=] (.*)\nalternative/) ;
                  $overlap_pvals->{"L_$ptype"."_less"} = $pval_less ;
               }

               if ($outnorm->{"L_and_$ptype"} > 0 ) {
                  if ($pilig_specs->{CALC_CONSSCORE_FLAG}) {
                     foreach my $tpos ($tand->Index_List_Read()) {
                        $tsum_nt += $class_aln->{meta}->{numrestypes_std20}->{$tpos};
                        $tsum_se += $class_aln->{meta}->{shannone_std20}->{$tpos} ;

                        my $resfreq_string = '';
                        foreach my $aatype (sort
                           keys %{$class_aln->{meta}->{res_freq}->{$tpos}}) {
                           $resfreq_string .= $aatype.
                            $class_aln->{meta}->{res_freq}->{$tpos}->{$aatype};
                        }

                        print OUTCONSSCORE join("\t",$class, $classtype,
                           "L_and_$ptype",
                           $tpos,
                           $class_aln->{meta}->{shannone_std20}->{$tpos},
                           $class_aln->{meta}->{numrestypes_std20}->{$tpos},
                           $resfreq_string
                        )."\n" ;

                     }
                     $consscore_se->{"L_and_$ptype"} =
                        $tsum_se / $outnorm->{"L_and_$ptype"} ;
                     $consscore_nt->{"L_and_$ptype"} =
                        $tsum_nt / $outnorm->{"L_and_$ptype"} ;
                  }
               }
            } elsif ($outnorm->{L} > 0) {

                  $outbits->{"L_or_$ptype"} = $outbits->{L} ;
                  $outnorm->{"L_or_$ptype"} = $outnorm->{L} ;

                  my $t = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
                  $consscore_se->{"L_and_$ptype"} = "undef" ;
                  $consscore_nt->{"L_and_$ptype"} = "undef" ;
                  $outbits->{"L_and_$ptype"} = $t->to_Bin() ;
                  $outnorm->{"L_and_$ptype"} = 0 ;

            } elsif ($outnorm->{$ptype} > 0) {

                  $outbits->{"L_or_$ptype"} = $outbits->{$ptype} ;
                  $outnorm->{"L_or_$ptype"} = $outnorm->{$ptype} ;

                  my $t = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
                  $consscore_se->{"L_and_$ptype"} = "undef" ;
                  $consscore_nt->{"L_and_$ptype"} = "undef" ;
                  $outbits->{"L_and_$ptype"} = $t->to_Bin() ;
                  $outnorm->{"L_and_$ptype"} = 0 ;

            } else {

                  my $t = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
                  $outbits->{"L_or_$ptype"} = $t->to_Bin() ;
                  $outnorm->{"L_or_$ptype"} = 0 ;
                  $outbits->{"L_and_$ptype"} = $t->to_Bin() ;
                  $outnorm->{"L_and_$ptype"} = 0 ;
                  $consscore_se->{"L_and_$ptype"} = "undef" ;
                  $consscore_nt->{"L_and_$ptype"} = "undef" ;

            }
         }

         my @outvals = ( $classtype, $class,
                         $class2alnlength->{$classtype}->{$class});

         foreach my $btype (qw/L P Pinter Pintra p E/) {
            push @outvals, $consscore_se->{$btype} ;
            push @outvals, $consscore_nt->{$btype} ;
            push @outvals, $outbits->{$btype} ;
            push @outvals, $outnorm->{$btype} ;
         }

         foreach my $btype (qw/P Pinter Pintra p/) {
            push @outvals, $consscore_se->{"L_and_$btype"} ;
            push @outvals, $consscore_nt->{"L_and_$btype"} ;
            push @outvals, $outbits->{"L_and_$btype"} ;
            push @outvals, $outnorm->{"L_and_$btype"} ;
            push @outvals, $outbits->{"L_or_$btype"} ;
            push @outvals, $outnorm->{"L_or_$btype"} ;
            push @outvals, $overlap_pvals->{"L_$btype"} ;
            push @outvals, $overlap_pvals->{"L_$btype".'_less'} ;
         }


         my $ligstring = '';
         if (exists $class2ligs->{$classtype}->{$class}) {
            my @tl ;
            foreach my $l (keys %{$class2ligs->{$classtype}->{$class}}){
               my $t = $l ; $t =~ s/\t/:/g ;
               push @tl, $t ;
            }
            $ligstring = join(',', @tl) ;
         }
         push @outvals, $ligstring ;
         foreach my $j ( 0 .. $#outvals) {
            if (!defined $outvals[$j]) {
               print STDERR " WARNING: $classtype $class field $j is undefined\n "; } }

         print join("\t", @outvals)."\n" ;
      }
   }
   if ($pilig_specs->{CALC_CONSSCORE_FLAG}) { close(OUTCONSSCORE) ;}

}


sub OLDPRE20090104_collate_perfam {

   require Bit::Vector ;
   require R;
   require RReferences ;
   R::initR("--silent") ;

   my $standardres = {
   ALA => 'A' ,
   ARG => 'R' ,
   ASN => 'N' ,
   ASP => 'D' ,
   CYS => 'C' ,
   GLN => 'Q' ,
   GLU => 'E' ,
   GLY => 'G' ,
   HIS => 'H' ,
   HSD => 'H' ,
   HSE => 'H' ,
   ILE => 'I' ,
   LEU => 'L' ,
   LYS => 'K' ,
   MET => 'M' ,
   PHE => 'F' ,
   PRO => 'P' ,
   SER => 'S' ,
   THR => 'T' ,
   TRP => 'W' ,
   TYR => 'Y' ,
   VAL => 'V',
   UNK => 'X',
   UNX => 'X',
   '  C' => 'c',
   '  G' => 'g',
   '  A' => 'a',
   '  T' => 't',
   '  U' => 'u',
   '  I' => 'i',
   'C' => 'c',
   'G' => 'g',
   'A' => 'a',
   'T' => 't',
   'U' => 'u',
   'I' => 'i',
   '+C' => 'c',
   '+G' => 'g',
   '+A' => 'a',
   '+T' => 't',
   '+U' => 'u',
   '+I' => 'i'
   } ;



# given deposit dirs,
# load list of all ASTEROIDS fam/sf

   my $pibase_specs = pibase::get_specs() ;
   my $pilig_specs = set_pilig_specs() ;
   my $fn = $pilig_specs->{fn} ;

   my $astral = _pilig_astral_preload({
      pibase_specs => $pibase_specs,
      pilig_specs => $pilig_specs,
   }) ;
   my $pb = _pilig_tod_pibase_preload() ;
   my $liginfo = _pilig_load_liginfo() ;

   my $class2alnlength = {};

   my $ligbits = readin_ligassignments({
      fn => $pilig_specs->{outfiles}->{assign_lig},
      liginfo => $liginfo,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;

   my $expbits = readin_expassignments({
      fn => $pilig_specs->{outfiles}->{assign_exp},
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;

   my $pibits_both = readin_piassignments({
      fn => $pilig_specs->{outfiles}->{assign_pi},
      pb => $pb,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;
   my $pibits = $pibits_both->{pibits} ;
   my $interfaces = $pibits_both->{interfaces} ;

   my $alltogether = combine_piligexpbits({
      pb => $pb,
      liginfo => $liginfo,
      ligbits => $ligbits,
      pibits => $pibits,
      expbits => $expbits,
   });
   my $class2bits = $alltogether->{class2bits} ;
   my $class2ligs = $alltogether->{class2ligs} ;
   my $class2anligbits = $alltogether->{class2anligbits} ;

   {
      my @headers ;
      push @headers, qw/classtype class alnlength/ ;
      my @tt = qw/Lmwinrange P Pinter Pintra E/ ;
      foreach my $t (@tt) {
         push @headers, "cons_shanonne_$t" ;
         push @headers, "cons_numtypes_$t" ;
         push @headers, "bits_$t" ;
         push @headers, "norm_$t" ;
      }

      foreach my $t (@tt[1..$#tt]) {
         push @headers, "cons_shannone_".$tt[0]."_and_$t" ;
         push @headers, "cons_numtypes_".$tt[0]."_and_$t" ;
         push @headers, "bits_".$tt[0]."_and_$t" ;
         push @headers, "norm_".$tt[0]."_and_$t" ;
         push @headers, "bits_".$tt[0]."_or_$t" ;
         push @headers, "norm_".$tt[0]."_or_$t" ;
         push @headers, "pval_".$tt[0]."_$t" ;
      }
      push @headers, 'ligcodes' ;

      foreach my $j (0 .. $#headers) {
         print '#'.($j + 1)."\t".$headers[$j]."\n" ; }
      print '#'.join("\t", @headers)."\n" ;
   }

   my $astral_classes = $astral->{classes} ;
#   foreach my $classtype (qw/seqcl90 seqcl95 seqcl100 fam sf/)
   if ($pilig_specs->{CALC_CONSSCORE_FLAG}) {
      open(OUTCONSSCORE, ">>outconsscore.$$.out") ;}
#   foreach my $classtype (qw/fam sf/)
   foreach my $classtype (qw/fam/) {
      my @classes ;
      my $curseqidlevel ;
      if ($classtype eq 'fam' || $classtype eq 'sf') {
         @classes = sort keys %{$astral_classes->{$classtype}} ;
      } elsif ($classtype =~ /^seqcl/) {
         my ($seqid) = ($classtype =~ /seqcl([0-9]+)/) ;
         $curseqidlevel = $seqid ;
         @classes = sort keys %{$astral->{seqcl2cont}->{$seqid}} ;
      }

      foreach my $class (@classes) {
         my $allprocalns ;
         my $outbits ;
         my $outnorm ;

         if (!exists $class2bits->{$classtype}->{$class}->{Lmwinrange} &&
             !exists $class2bits->{$classtype}->{$class}->{P}) {
            my @outvals = ( $classtype, $class, 'NOINTERACTIONS' );
            print join("\t", @outvals)."\n" ;
            next;
         }

         print STDERR "NOW ON $classtype $class (line ".__LINE__.")\n" ;
         my $class_aln ;
         if ($pilig_specs->{CALC_CONSSCORE_FLAG}) {
            my ($curfam, $curclasstype) ;
            if ($classtype eq 'fam' || $classtype eq 'sf') {
               $curfam = $class;
#BUG070503_1511                $curclasstype = $class ;
               $curclasstype = $classtype ;
            } elsif ($classtype =~ /^seqcl/) {
               $curfam = $pb->{osid2class}->{fam}->{$class} ;
               $curclasstype = 'fam' ;
            }

            if (! exists $allprocalns->{$curclasstype}->{$curfam}) {
#               print STDERR "OHFUCK $curclasstype $curfam aln file is ".
#                  $fn->{aln}->{$curclasstype}.'/'.$curfam.'.aln.fa'."\n" ;
               $allprocalns->{$curclasstype}->{$curfam} =
                  pibase::ASTRAL::load_asteroids_aln({
                  aln_fn => $pibase_specs->{asteroids}->{$curclasstype.'_aln'}.
                  '/'.$curfam.'.fasta_aln' ,
                  seq_fn => $pibase_specs->{asteroids}->{$curclasstype.'_seq'}.
                  '/'.$curfam.'.fa' ,
                  raf => $astral->{raf},
                  gdseqh => $astral->{gdseqh},
                  seqclcont100 => $astral->{seqcl2cont}->{100},
                  seqcl100 => $astral->{seqcl}->{100},
                  allchains => $pb->{pdbchains}
               }) ;
            }

            $class_aln = $allprocalns->{$curclasstype}->{$curfam} ;
         }

         my $consscore_se ;
         my $consscore_nt ;

         foreach my $btype (qw/Lmwinrange Pinter Pintra P E/) {
            if (exists $class2bits->{$classtype}->{$class}->{$btype}) {
               $outnorm->{$btype} = 
               $class2bits->{$classtype}->{$class}->{$btype}->Norm();

               $outbits->{$btype} =
               $class2bits->{$classtype}->{$class}->{$btype}->to_Bin();

#               if ($classtype eq 'fam' || $classtype eq 'sf')
               if ($pilig_specs->{CALC_CONSSCORE_FLAG} && 
                   ($classtype eq 'fam' || $classtype eq 'sf')) {
                  my $tsum_se = 0;
                  my $tsum_nt = 0;
                  foreach my $tpos ($class2bits->{$classtype}->{$class}->{$btype}->Index_List_Read()) {
                     print OUTCONSSCORE join("\t",$class, $classtype,$btype,
                        $tpos,
                        $class_aln->{meta}->{shannone}->{$tpos},
                        $class_aln->{meta}->{numrestypes}->{$tpos})."\n" ;

                     $tsum_se += $class_aln->{meta}->{shannone}->{$tpos} ;
                     $tsum_nt += $class_aln->{meta}->{numrestypes}->{$tpos} ;
                  }
                  if ($outnorm->{$btype} > 0) {
                     $consscore_se->{$btype} = $tsum_se / $outnorm->{$btype} ;
                     $consscore_nt->{$btype} = $tsum_nt / $outnorm->{$btype} ;
                  } else {
                     $consscore_se->{$btype} = 'undef' ;
                     $consscore_nt->{$btype} = 'undef' ;
                  }
               } else {
                  $consscore_se->{$btype} = "undefseqcl" ;
                  $consscore_nt->{$btype} = "undefseqcl" ;
               }

            } else {
               my $t = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
               $consscore_se->{$btype} = "undef" ;
               $consscore_nt->{$btype} = "undef" ;
               $outbits->{$btype} = $t->to_Bin() ;
               $outnorm->{$btype} = 0 ;
            }
         }

         my $overlap_pvals ;
         foreach my $ptype (qw/P Pinter Pintra/) {
            $overlap_pvals->{"Lmwinrange_$ptype"} = "undef" ;
            if ($outnorm->{Lmwinrange} > 0 && $outnorm->{$ptype} > 0) {
               my $tor = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
               $tor->Or(
                  $class2bits->{$classtype}->{$class}->{$ptype},
                  $class2bits->{$classtype}->{$class}->{Lmwinrange}) ;

               $outbits->{"Lmwinrange_or_$ptype"} = $tor->to_Bin() ;
               $outnorm->{"Lmwinrange_or_$ptype"} = $tor->Norm() ;

               my $tand = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
               $tand->And(
                     $class2bits->{$classtype}->{$class}->{$ptype},
                     $class2bits->{$classtype}->{$class}->{Lmwinrange}) ;

               $outbits->{"Lmwinrange_and_$ptype"} = $tand->to_Bin() ;
               $outnorm->{"Lmwinrange_and_$ptype"} = $tand->Norm() ;
               $consscore_se->{"Lmwinrange_and_$ptype"} = "undef" ;
               $consscore_nt->{"Lmwinrange_and_$ptype"} = "undef" ;

               my $tsum_se = 0;
               my $tsum_nt = 0;
               $overlap_pvals->{"Lmwinrange_$ptype"} = 'pvalnotcalc' ;
               if ($outnorm->{"Lmwinrange_and_$ptype"} > 0 ) {
                  if ($pilig_specs->{CALC_CONSSCORE_FLAG}) {
                     foreach my $tpos ($tand->Index_List_Read()) {
                        $tsum_nt += $class_aln->{meta}->{numrestypes}->{$tpos} ;
                        $tsum_se += $class_aln->{meta}->{shannone}->{$tpos} ;

                     print OUTCONSSCORE join("\t",$class, $classtype,
                        "Lmwinrange_and_$ptype",
                        $tpos,
                        $class_aln->{meta}->{shannone}->{$tpos},
                        $class_aln->{meta}->{numrestypes}->{$tpos})."\n" ;

                     }
                     $consscore_se->{"Lmwinrange_and_$ptype"} =
                        $tsum_se / $outnorm->{"Lmwinrange_and_$ptype"} ;
                     $consscore_nt->{"Lmwinrange_and_$ptype"} =
                        $tsum_nt / $outnorm->{"Lmwinrange_and_$ptype"} ;
                  }

                  if ($pilig_specs->{CALC_PVAL_FLAG}) {
                     my $tnorm ;
                     $tnorm->{lp} = $outnorm->{"E"} - $outnorm->{"Lmwinrange_or_$ptype"} ;
                     $tnorm->{Lp} = $outnorm->{"Lmwinrange"} - $outnorm->{"Lmwinrange_and_$ptype"} ;
                     $tnorm->{lP} = $outnorm->{$ptype} - $outnorm->{"Lmwinrange_and_$ptype"} ;
                     $tnorm->{LP} = $outnorm->{"Lmwinrange_and_$ptype"} ;

                     my @x = &R::eval("capture.output(fisher.test(matrix(c($tnorm->{lp}, $tnorm->{Lp}, $tnorm->{lP}, $tnorm->{LP}), 2, 2),workspace=2e7))") ;
                     my $x = join("\n", @x) ;
                     my ($pval) = ($x =~ /p-value [\<=] (.*)\nalternative/) ;
                     $overlap_pvals->{"Lmwinrange_$ptype"} = $pval ;
                  }
               }
            } elsif ($outnorm->{Lmwinrange} > 0) {

                  $outbits->{"Lmwinrange_or_$ptype"} = $outbits->{Lmwinrange} ;
                  $outnorm->{"Lmwinrange_or_$ptype"} = $outnorm->{Lmwinrange} ;

                  my $t = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
                  $consscore_se->{"Lmwinrange_and_$ptype"} = "undef" ;
                  $consscore_nt->{"Lmwinrange_and_$ptype"} = "undef" ;
                  $outbits->{"Lmwinrange_and_$ptype"} = $t->to_Bin() ;
                  $outnorm->{"Lmwinrange_and_$ptype"} = 0 ;

            } elsif ($outnorm->{$ptype} > 0) {

                  $outbits->{"Lmwinrange_or_$ptype"} = $outbits->{$ptype} ;
                  $outnorm->{"Lmwinrange_or_$ptype"} = $outnorm->{$ptype} ;

                  my $t = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
                  $consscore_se->{"Lmwinrange_and_$ptype"} = "undef" ;
                  $consscore_nt->{"Lmwinrange_and_$ptype"} = "undef" ;
                  $outbits->{"Lmwinrange_and_$ptype"} = $t->to_Bin() ;
                  $outnorm->{"Lmwinrange_and_$ptype"} = 0 ;

            } else {

                  my $t = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
                  $outbits->{"Lmwinrange_or_$ptype"} = $t->to_Bin() ;
                  $outnorm->{"Lmwinrange_or_$ptype"} = 0 ;
                  $outbits->{"Lmwinrange_and_$ptype"} = $t->to_Bin() ;
                  $outnorm->{"Lmwinrange_and_$ptype"} = 0 ;
                  $consscore_se->{"Lmwinrange_and_$ptype"} = "undef" ;
                  $consscore_nt->{"Lmwinrange_and_$ptype"} = "undef" ;

            }
         }

         my @outvals = ( $classtype, $class,
                         $class2alnlength->{$classtype}->{$class});

         foreach my $btype (qw/Lmwinrange P Pinter Pintra E/) {
            push @outvals, $consscore_se->{$btype} ;
            push @outvals, $consscore_nt->{$btype} ;
            push @outvals, $outbits->{$btype} ;
            push @outvals, $outnorm->{$btype} ;
         }

         foreach my $btype (qw/P Pinter Pintra/) {
            push @outvals, $consscore_se->{"Lmwinrange_and_$btype"} ;
            push @outvals, $consscore_nt->{"Lmwinrange_and_$btype"} ;
            push @outvals, $outbits->{"Lmwinrange_and_$btype"} ;
            push @outvals, $outnorm->{"Lmwinrange_and_$btype"} ;
            push @outvals, $outbits->{"Lmwinrange_or_$btype"} ;
            push @outvals, $outnorm->{"Lmwinrange_or_$btype"} ;
            push @outvals, $overlap_pvals->{"Lmwinrange_$btype"} ;
         }


         my $ligstring = '';
         if (exists $class2ligs->{$classtype}->{$class}->{Lmwinrange}) {
            my @tl ;
            foreach my $l (keys %{$class2ligs->{$classtype}->{$class}->{Lmwinrange}}){
               my $t = $l ; $t =~ s/\t/:/g ;
               push @tl, $t ;
            }
            $ligstring = join(',', @tl) ;
         }
         push @outvals, $ligstring ;
         foreach my $j ( 0 .. $#outvals) {
            if (!defined $outvals[$j]) {
               print STDERR " WARNING: $classtype $class field $j is undefined\n "; } }

         print join("\t", @outvals)."\n" ;
      }
   }
   if ($pilig_specs->{CALC_CONSSCORE_FLAG}) { close(OUTCONSSCORE) ;}
}


sub _pilig__getinput_assign_pi {

   my $bdp2contacts_fn;
   my $bdp2sid;
   my $class2sid ;
   my $sid2class ;
   while (my $line = <STDIN>) {
      if ($line =~ /^\#/) {next;} ;
      chomp $line ;
      my ($bdp, $contacts_fn, $sid1, $class1, $sid2, $class2) =
         split(/\t/, $line);
      $bdp2contacts_fn->{$bdp} = $contacts_fn ;
      $class2sid->{$class1}->{$sid1}++ ;
      $class2sid->{$class2}->{$sid2}++ ;
      $sid2class->{$sid1}->{$class1}++ ;
      $sid2class->{$sid2}->{$class2}++ ;
   }

   return {
      sid2class => $sid2class,
      class2sid => $class2sid,
      bdp2contacts_fn=> $bdp2contacts_fn
   } ;

}

sub _pilig__getinput_assign_lig {

   my $pdb2res2ligid ;
   while (my $line = <STDIN>) {
      if ($line =~ /^\#/) {next;}
      chomp $line ;
      my ($pdb, $resno, $resch, $ligcod, $ligid) = split(/\t/, $line) ;
      $pdb2res2ligid->{$pdb}->{$resno."\n".$resch}->{$ligcod."\n".$ligid}++;
   }

   return {
      pdb2res2ligid => $pdb2res2ligid
   } ;
}


sub _pilig__getinput_assign_exp {

   my $sid2class ; my $class2sid ;
   my $sid2bdp ;
   my $sid2expres ;
   while (my $line = <STDIN>) {
      if ($line =~ /^\#/) {next;}
      chomp $line ;
      my @t = split(/\t/, $line) ;
      my $sid = shift @t ;
      my $bdp = shift @t ;
      my $class = shift @t ;

      foreach my $res ( @t) {
         substr($res,-2,1) = "\n" ;
         push @{$sid2expres->{$sid}}, $res ;
      }

      $sid2class->{$sid} = $class ;

      my ($sf) = ($class =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
      $sid2class->{fam}->{$sid} = $class;
      $sid2class->{sf}->{$sid} = $sf ;
      $class2sid->{fam}->{$class}->{$sid}++ ;
      $class2sid->{sf}->{$sf}->{$sid}++ ;
   }

   return {
      sid2class => $sid2class,
      class2sid => $class2sid,
      sid2expres=> $sid2expres,
   } ;

}

sub _pilig_load_sid_sasa {

   my $in = shift ;

   my $sid2class ; my $class2sid ;
   my $sid2bdp ;
   my $sid2expres ;
   open(SASASID_FN, $in->{fn}) ;
   while (my $line = <SASASID_FN>) {
      if ($line =~ /^\#/) {next;}
      chomp $line ;
      my @t = split(/\t/, $line) ;
      my $sid = shift @t ;
      my $bdp = shift @t ;
      my $class = shift @t ;

      if (exists $in->{sids} &&
          !exists $in->{sids}->{$sid}) {next;}

      foreach my $res ( @t) {
         substr($res,-2,1) = "\n" ;
         push @{$sid2expres->{$sid}}, $res ;
      }

      $sid2class->{$sid} = $class ;

      my ($sf) = ($class =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
      $sid2class->{fam}->{$sid} = $class;
      $sid2class->{sf}->{$sid} = $sf ;
      $class2sid->{fam}->{$class}->{$sid}++ ;
      $class2sid->{sf}->{$sf}->{$sid}++ ;
   }
   close(SASASID_FN) ;

   return {
      sid2class => $sid2class,
      class2sid => $class2sid,
      sid2expres=> $sid2expres,
   } ;

}


sub _pilig_load_pepnuci_realres {

   my $in = shift ;

   my $sid2class ; my $class2sid ;
   my $sid2bdp ;
   my $sid2pepnucres;
   my $target_chain ;
   open(PEPNUCI_REALRESF, $in->{fn}) ;
   while (my $line = <PEPNUCI_REALRESF>) {
      if ($line =~ /^\#/) {next;}
      chomp $line ;
      my @t = split(/\t/, $line) ;
      my $bdp = shift @t ;
      my $sid = shift @t ;
      my $osid = shift @t ;
      my $classtype = shift @t ;
      my $class = shift @t ;
      my $target_chain_type = shift @t ;
      my $target_chain_length = shift @t ;
      my $target_chain_no = shift @t ;
      my $target_chain_dom = shift @t ;

      $target_chain->{$target_chain_dom}->{chain_length} = $target_chain_length;
      $target_chain->{$target_chain_dom}->{chain_type} = $target_chain_type;

      if (exists $in->{sids} &&
          !exists $in->{sids}->{$sid}) {next;}

      foreach my $res ( @t) {
         substr($res,-2,1) = "\n" ;
         push @{$sid2pepnucres->{$sid}->{$target_chain_dom}}, $res ;
      }

      $sid2class->{$sid} = $class ;

      my ($sf) = ($class =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
      $sid2class->{fam}->{$sid} = $class;
      $sid2class->{sf}->{$sid} = $sf ;
      $class2sid->{fam}->{$class}->{$sid}++ ;
      $class2sid->{sf}->{$sf}->{$sid}++ ;
   }
   close(PEPNUCI_REALRESF) ;

   return {
      sid2class => $sid2class,
      target_chain => $target_chain,
      class2sid => $class2sid,
      sid2pepnucres => $sid2pepnucres,
   } ;

}


sub _pilig_run_assign_lig {

   my $in = _getinput_assign_lig() ;
   my $lb = { pdb2res2ligid => $in->{pdb2res2ligid} } ;
   my $fn = set_locations() ;
   my $astral = astral_preload({ fn => $fn }) ;
   my $pb = _pilig_tod_pibase_preload() ;
   my $lbdoms = get_ligbase_domains({lb => $lb, pb => $pb}) ;


   foreach my $classtype (qw/fam sf/) {
      foreach my $class (sort keys %{$astral->{classes}->{$classtype}}) {
         if (!exists $lbdoms->{class2ligsid}->{$classtype}->{$class}) {next;}

         print STDERR "NOW ON $classtype $class (line ".__LINE__.")\n" ;
         my $pos2bs ;
         my $class_aln = pibase::ASTRAL::load_asteroids_aln({
            aln_fn => $fn->{aln}->{$classtype}.'/'.$class.'.aln.fa',
            seq_fn => $fn->{aln}->{$classtype}.'/'.$class.'.fa',
            raf => $astral->{raf},
            gdseqh => $astral->{gdseqh},
            seqclcont100 => $astral->{seqcl2cont}->{100},
            seqcl100 => $astral->{seqcl}->{100},
            allchains => $pb->{pdbchains}
         }) ;

         foreach my $ligsid (sort keys %{$lbdoms->{class2ligsid}->{$classtype}->{$class}}) {
            my $osid = $pb->{sid2osid}->{$ligsid} ;
            my $pdb = $pb->{sid2pdb}->{$ligsid} ;
            my ($talnres, $undefres) ;
            foreach my $res (sort keys %{$lbdoms->{sid2ligres}->{$ligsid}}) {
               my $alnres = $class_aln->{resno2pos}->{$osid}->{$res} ;
               if (!defined $alnres) {
#                  print STDERR "WARNING: $osid ($ligsid) residue $res not ".
#                               "found in ASTRAL alignment\n" ;
                  $alnres = 'undef' ;
               }
#               $pos2bs->{$alnres}->{l}->{$ligsid}++ ;
               foreach my $lig (sort keys %{$lb->{pdb2res2ligid}->{$pdb}->{$res}}) {
                  if ($alnres ne 'undef' ) {
                     push @{$talnres->{$lig}}, $alnres ;
                  } else {
                     push @{$undefres->{$lig}}, $alnres ;
                  }
               }
            }

            foreach my $lig (keys %{$talnres}) {
               my ($ligcod, $ligid) = split(/\n/, $lig) ;

               my @salnres = ();
               if (exists $talnres->{$lig}) {
                  push @salnres, sort {$a <=> $b} @{$talnres->{$lig}} ; }
               if (exists $undefres->{$lig}) {
                  push @salnres, @{$undefres->{$lig}} ; }
               my $alnposstring = join(',', @salnres) ;

               my @outvals = ( $pdb, $ligsid, $osid, $classtype,
                                 $class, "L", $alnposstring,
                                 $class_aln->{alnlength},
                                 $ligcod, $ligid) ;
               print join("\t", @outvals)."\n" ;
            }
         }
      }
   }

}


sub assign_lig {

   my $in = shift ;

   my $pilig_specs = set_pilig_specs() ;
   my $pibase_specs = pibase::get_specs() ;

#old way:
#   my $astral = astral_preload({ fn => $fn }) ;

   my $astral = _pilig_astral_preload({
      pibase_specs => $pibase_specs,
      pilig_specs => $pilig_specs,
   }) ;

   my $class_list ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         chomp $line;
         my ($classtype, $class) = split(/\t/, $line);
         $class_list->{$classtype}->{$class}++ ;
      }
   } else {
      foreach my $classtype (keys %{$astral->{classes}}) {
         foreach my $class (keys %{$astral->{classes}->{$classtype}}) {
            $class_list->{$classtype}->{$class}++ ; } }
   }


   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
# if this is master script, split input, recall self with
# in_fn specified, cluster_fl = 0, import_fl = 0
# send out, cluster run, and return merged

      print "* assign_lig() ".localtime() if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{assign_lig_in}, $temp_fn->{assign_lig_in}) =
         tempfile("splits_assign_lig_input.XXXXX");
      ($temp_fh->{assign_lig_out}, $temp_fn->{assign_lig_out}) =
         tempfile("splits_assign_lig_SGEout_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{assign_lig_out}) ;
      ($temp_fh->{assign_lig_err}, $temp_fn->{assign_lig_err}) =
         tempfile("splits_assign_lig_SGEerr_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{assign_lig_err}) ;

      foreach my $class_type (sort keys %{$class_list}) {
         foreach my $class (sort keys %{$class_list->{$class_type}}) {
         print {$temp_fh->{assign_lig_in}}
            join("\t", $class_type, $class)."\n" ; } }
      close($temp_fh->{assign_lig_in}) ;

      my $split_dir = tempdir("splits_assign_lig.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{assign_lig_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.assign_lig.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::pilig qw/assign_lig/ ;

main() ;

sub main {

         pibase::pilig::assign_lig({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.assign_lig.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.assign_lig.XXXXX");

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
         out_fn => $temp_fn->{assign_lig_out},
         err_fn => $temp_fn->{assign_lig_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{assign_lig_out},
           $temp_fn->{assign_lig_out}) ;
      open(REALOUTF,">".$pilig_specs->{outfiles}->{assign_lig}) ;
      while (my $line = readline($temp_fh->{assign_lig_out})) {
         if ($line =~ /^\#/) {next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{assign_lig_out}) ;
      close(REALOUTF) ;

   } else {

#      print STDERR "loading ligbase: ".localtime() ;
#      my $lb = _pilig_tod_ligbase_preload() ;
#      print STDERR ", DONE: ".localtime()."\t" ;
#      print STDERR " cur size ".total_size($lb)." bytes\n";

# switched to per-class loading of ligbase data

      print STDERR "loading pibase: ".localtime() ;
      my $pb = _pilig_tod_pibase_preload() ;
      print STDERR ", DONE: ".localtime()."\n" ;
#      print STDERR " cur size ".total_size($pb)." bytes\n";


# TESTING FUCKER fpd080630_0210 
#   my ($t_subset_id, $t_chain_id, $t_resno_serial, $t_resno) =
#      pibase::rawselect_metatod("/groups/eddy/home/davisf/work/pibase/pibase200708/data/metatod/subsets_residues/32/subsets_residues_32533.32533.pibase.gz",
#      "SELECT subset_id, chain_id, resno_serial, resno FROM /groups/eddy/home/davisf/work/pibase/pibase200708/data/metatod/subsets_residues/32/subsets_residues_32533.32533.pibase.gz") ;
#   foreach my $j ( 0 .. $#{$t_subset_id}) {
#      print STDERR join("\t",$j,$t_subset_id->[$j], $t_chain_id->[$j], $t_resno->[$j])."\n" ;
#   }
#   die ;
# ---TESTING FUCKER---

# now loading per-family
#      print STDERR "loading ligbase domains : " ;
#      my $lbdoms = _pilig_get_ligbase_domains({ lb => $lb, pb => $pb }) ;
#      print STDERR "X\n ";
#
#      print STDERR "ALL PRELOADS DONE - HERE NOW ".__LINE__."\n" ;

      foreach my $classtype (sort keys %{$class_list}) {
         foreach my $class (sort keys %{$class_list->{$classtype}}) {

            my $lb = _pilig_tod_ligbase_preload_loadperclass({
               classtype => $classtype,
               class => $class,
               class2pdb => $pb->{class2pdb}
            }) ;
            my $lbdoms = _pilig_get_ligbase_domains({ lb => $lb,
                                                      pb => $pb });

            if (!exists $lbdoms->{class2ligsid}->{$classtype}->{$class}) {next;}
   
            print STDERR "NOW ON $classtype $class (line ".__LINE__.")\n" ;
            my $pos2bs ;
            my $class_aln = pibase::ASTRAL::load_asteroids_aln({
               aln_fn => $pibase_specs->{asteroids}->{$classtype.'_aln'}.
                  '/'.$class.'.fasta_aln' ,
               seq_fn => $pibase_specs->{asteroids}->{$classtype.'_seq'}.
                  '/'.$class.'.fa' ,
               raf => $astral->{raf},
               gdseqh => $astral->{gdseqh},
               seqclcont100 => $astral->{seqcl2cont}->{100},
               seqcl100 => $astral->{seqcl}->{100},
               allchains => $pb->{pdbchains}
            }) ;

            foreach my $ligsid (sort
                  keys %{$lbdoms->{class2ligsid}->{$classtype}->{$class}}) {

               my $osid = $pb->{sid2osid}->{$ligsid} ;
               my $pdb = $pb->{sid2pdb}->{$ligsid} ;
               my ($talnres, $undefres) ;
               foreach my $res (sort keys %{$lbdoms->{sid2ligres}->{$ligsid}}) {
                  my $alnres = $class_aln->{resno2pos}->{$osid}->{$res} ;
                  if (!defined $alnres) {
   #                  print STDERR "WARNING: $osid ($ligsid) residue $res not ".
   #                               "found in ASTRAL alignment\n" ;
                     $alnres = 'undef' ;
                  }
   #               $pos2bs->{$alnres}->{l}->{$ligsid}++ ;
                  foreach my $lig (sort keys %{$lb->{pdb2res2ligid}->{$pdb}->{$res}}) {
                     if ($alnres ne 'undef' ) {
                        push @{$talnres->{$lig}}, $alnres ;
                     } else {
                        push @{$undefres->{$lig}}, $alnres ;
                     }
                  }
               }
   
               foreach my $lig (keys %{$talnres}) {
                  my ($ligcod, $ligid) = split(/\n/, $lig) ;
   
                  my @salnres = ();
                  if (exists $talnres->{$lig}) {
                     push @salnres, sort {$a <=> $b} @{$talnres->{$lig}} ; }
                  if (exists $undefres->{$lig}) {
                     push @salnres, @{$undefres->{$lig}} ; }
                  my $alnposstring = join(',', @salnres) ;
   
                  my @outvals = ( $pdb, $ligsid, $osid, $classtype,
                                    $class, "L", $alnposstring,
                                    $class_aln->{alnlength},
                                    $ligcod, $ligid) ;
                  print join("\t", @outvals)."\n" ;
               }
            }
         }
      }
   }

}


sub assign_pi {

   my $in = shift ;

   my $pilig_specs = set_pilig_specs() ;
#   my $fn = $pilig_specs->{fn} ;
   my $pibase_specs = pibase::get_specs() ;

#old way:
#   my $astral = astral_preload({ fn => $fn }) ;

   my $astral = _pilig_astral_preload({
      pibase_specs => $pibase_specs,
      pilig_specs => $pilig_specs,
   }) ;

   my $class_list ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         chomp $line;
         my ($classtype, $class) = split(/\t/, $line);
         $class_list->{$classtype}->{$class}++ ;
      }
   } else {
      foreach my $classtype (keys %{$astral->{classes}}) {
         foreach my $class (keys %{$astral->{classes}->{$classtype}}) {
            $class_list->{$classtype}->{$class}++ ; } }
   }

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
# if this is master script, split input, recall self with
# in_fn specified, cluster_fl = 0
# send out, cluster run, and return merged

      print "* assign_pi() ".localtime() if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{assign_pi_in}, $temp_fn->{assign_pi_in}) =
         tempfile("splits_assign_pi_input.XXXXX");
      ($temp_fh->{assign_pi_out}, $temp_fn->{assign_pi_out}) =
         tempfile("splits_assign_pi_SGEout_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{assign_pi_out}) ;
      ($temp_fh->{assign_pi_err}, $temp_fn->{assign_pi_err}) =
         tempfile("splits_assign_pi_SGEerr_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{assign_pi_err}) ;

      foreach my $class_type (sort keys %{$class_list}) {
         foreach my $class (sort keys %{$class_list->{$class_type}}) {
         print {$temp_fh->{assign_pi_in}}
            join("\t", $class_type, $class)."\n" ; } }
      close($temp_fh->{assign_pi_in}) ;

      my $split_dir = tempdir("splits_assign_pi.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{assign_pi_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.assign_pi.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::pilig qw/assign_pi/ ;

main() ;

sub main {

         pibase::pilig::assign_pi({
            cluster_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.assign_pi.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.assign_pi.XXXXX");

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
         out_fn => $temp_fn->{assign_pi_out},
         err_fn => $temp_fn->{assign_pi_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{assign_pi_out},
           $temp_fn->{assign_pi_out}) ;
      open(REALOUTF,">".$pilig_specs->{outfiles}->{assign_pi}) ;
      while (my $line = readline($temp_fh->{assign_pi_out})) {
         if ($line =~ /^\#/) {next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{assign_pi_out}) ;
      close(REALOUTF) ;

   } else {
   
      my $pb = _pilig_tod_pibase_preload() ;

#      foreach my $classtype (qw/fam sf/)
#         foreach my $class (sort keys %{$astral->{classes}->{$classtype}})
      foreach my $classtype (sort keys %{$class_list}){
         foreach my $class (sort keys %{$class_list->{$classtype}}){
            if (!exists $pb->{class2sidpairs}->{$classtype}->{$class}) { next;}
   
            print STDERR "NOW ON $classtype $class (line ".__LINE__.")\n" ;
            my $pos2bs ;
            my $class_aln = pibase::ASTRAL::load_asteroids_aln({
               aln_fn => $pibase_specs->{asteroids}->{$classtype.'_aln'}.
                  '/'.$class.'.fasta_aln' ,
               seq_fn => $pibase_specs->{asteroids}->{$classtype.'_seq'}.
                  '/'.$class.'.fa' ,
               raf => $astral->{raf},
               gdseqh => $astral->{gdseqh},
               seqclcont100 => $astral->{seqcl2cont}->{100},
               seqcl100 => $astral->{seqcl}->{100},
               allchains => $pb->{pdbchains}
            }) ;
   
            foreach my $j (sort {$a <=> $b}
                       keys %{$pb->{class2sidpairs}->{$classtype}->{$class}} ){
   
               my $pdb = $pb->{sid2pdb}->{$pb->{sid1}->[$j]} ;
   
               my $interface = _pilig_load_interface_contacts({
                  sid1 => $pb->{sid1}->[$j],
                  sid2 => $pb->{sid2}->[$j],
                  dist_thresh => $pilig_specs->{PARAM_INTERFACE_DIST_THRESH},
                  fn => $pb->{bdp2contactsfn}->{$pb->{bdp_id}->[$j]}
               }) ;
               my $sid12 = $pb->{sid1}->[$j]."\t".$pb->{sid2}->[$j] ;
   
               my @dothese = () ;
               if ($pb->{sid2class}->{$classtype}->{$pb->{sid1}->[$j]} eq
                   $class) {push @dothese, $pb->{sid1}->[$j];}
   
               if ($pb->{sid2class}->{$classtype}->{$pb->{sid2}->[$j]} eq
                   $class) {push @dothese, $pb->{sid2}->[$j];}
   
               foreach my $sid (@dothese) {
                  my $osid = $pb->{sid2osid}->{$sid} ;
                  my @talnres = () ;
                  my @undefres = () ;
                  foreach my $res (sort keys %{$interface->{intres}->{$sid}}) {
                     my $alnres = $class_aln->{resno2pos}->{$osid}->{$res} ;
                     if (!defined $alnres) {
                        push @undefres, 'undef' ;
   #                     print STDERR "WARNING: $osid ($sid) residue $res not ".
   #                                  "found in ASTRAL alignment\n" ;
                     } else {
                        push @talnres, $alnres ;
                     }
   #                  $pos2bs->{$alnres}->{p}->{$sid}++ ;
                  }
                  my @salnres;
                  push @salnres, sort {$a <=> $b} @talnres ;
                  push @salnres, @undefres ;
                  my $alnposstring = join(',', @salnres) ;
                  my @outvals = ($pdb,  $sid, $osid, $classtype,
                                 $class, "P", $alnposstring,
                                 $class_aln->{alnlength},
                                 $pb->{sid1}->[$j], $pb->{sid2}->[$j],
                                 $pb->{sid2class}->{'fam'}->{$pb->{sid1}->[$j]},
                                 $pb->{sid2class}->{'fam'}->{$pb->{sid2}->[$j]},
                                 $pb->{sid12chain}->{$sid12},
                  ) ;
                  print join("\t", @outvals)."\n" ;
               }
            }
         }
      }
   }

}


#fpd081229_1639 - code to annotate domain - peptide chain interactions.
#  takes advantage of enumerated chain-chain contacts in new pibase.
#
# this calculates comprehensively for all BDP entries (raw and non-raw)
#  and in real residue numbers - consider as regular part of PIBASE dataset
#
# at least make available for download: all 
#
sub assign_pepnuci_realres {

   my $in = shift ;

   my $pilig_specs = set_pilig_specs() ;
   my $pibase_specs = pibase::get_specs() ;

   my $bdp_list ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         chomp $line;
         my ($bdp) = split(/\t/, $line);
         $bdp_list->{$bdp}++ ;
      }
   } else {
            #load all BDPs with at least 1 non-chain domain defined.
            # keep it scop only for now - consider generalizing to CATH also
            # for integration into PIBASE

      my ($t_sid, $t_bdp, $t_ssid) = pibase::rawselect_tod(
         "SELECT subset_id, bdp_id, subset_source_id FROM subsets") ;
      
      foreach my $j (0 .. $#{$t_sid}) {
         if ($t_bdp->[$j] eq 'NULL' ||
             $t_sid->[$j] =~ /CHAIN/) {next;}
         $bdp_list->{$t_bdp->[$j]}++ ;
      }
   }

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
# if this is master script, split input, recall self with
# in_fn specified, cluster_fl = 0, import_fl = 0
# send out, cluster run, and return merged

      print "* assign_pepnuci_realres() ".localtime()
         if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{assign_pepnuci_realres_in},
       $temp_fn->{assign_pepnuci_realres_in}) =
         tempfile("splits_assign_pepnuci_realres_input.XXXXX");

      ($temp_fh->{assign_pepnuci_realres_out},
       $temp_fn->{assign_pepnuci_realres_out}) =
         tempfile("splits_assign_pepnuci_realres_SGEout_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{assign_pepnuci_realres_out}) ;
      ($temp_fh->{assign_pepnuci_realres_err},
       $temp_fn->{assign_pepnuci_realres_err}) =
         tempfile("splits_assign_pepnuci_realres_SGEerr_XXXXX",
                  SUFFIX => '.pibase');
         close($temp_fh->{assign_pepnuci_realres_err}) ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}) {
         print {$temp_fh->{assign_pepnuci_realres_in}} $bdp_id."\n";}
      close($temp_fh->{assign_pepnuci_realres_in}) ;

      my $split_dir = tempdir("splits_assign_pepnuci_realres.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{assign_pepnuci_realres_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.assign_pepnuci_realres.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::pilig qw/assign_pepnuci_realres/ ;

main() ;

sub main {

         pibase::pilig::assign_pepnuci_realres({
            cluster_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.assign_pepnuci_realres.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.assign_pepnuci_realres.XXXXX");

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
         out_fn => $temp_fn->{assign_pepnuci_realres_out},
         err_fn => $temp_fn->{assign_pepnuci_realres_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{assign_pepnuci_realres_out},
           $temp_fn->{assign_pepnuci_realres_out}) ;
      open(REALOUTF,">".$pilig_specs->{outfiles}->{assign_pepnuci_realres}) ;
      while (my $line = readline($temp_fh->{assign_pepnuci_realres_out})) {
         if ($line =~ /^\#/) {next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{assign_pepnuci_realres_out}) ;
      close(REALOUTF) ;

   } else {
   
      my $pb = _pilig_tod_pibase_preload() ;

      my $chainchain_contacts ;
      {
         my ($t_bdp, $t_sid1, $t_sid2) = pibase::rawselect_tod(
         "SELECT bdp_id, subset_id_1, subset_id_2 FROM intersubset_contacts") ;
         foreach my $j ( 0.. $#{$t_bdp}) {
            if ($t_sid1->[$j] !~ /CHAIN/) {next;}
            my $sid12 = $t_sid1->[$j]."\t".$t_sid2->[$j] ;
            $chainchain_contacts->{$t_bdp->[$j]}->{$t_sid1->[$j]}->{$sid12} =
               $t_sid2->[$j] ;
            $chainchain_contacts->{$t_bdp->[$j]}->{$t_sid2->[$j]}->{$sid12} =
               $t_sid1->[$j] ;
         }
      }

      foreach my $bdp_id (sort {$a <=> $b} keys %{$bdp_list}){
# 1. identify peptide chains (PC)
# 2. identify PC-interacting chains (PCIC)
# 3. load all PC-PCIC residue contacts
# 4. identify residues in PCIC with domain definitions
# 5. list domain-PC pairs, and the residues involved in the interaction.
#

# 1. identify peptide chains (PC)
         print STDERR "NOW ON $bdp_id\n" ;
         my $cur_subsres = _identify_peptide_and_nucleic_chains({
            bdp_id => $bdp_id,
            pb => $pb,
            subsetsres_fn => $pb->{bdp2subsetsresfn}->{$bdp_id}
         }) ;

# 2. determine peptide chain-interacting chains (PCIC) with domain definitions
         foreach my $t_chain (keys %{$cur_subsres->{peptide_chains}}) {
#            if ($cur_subsres->{chain_type}->{$t_chain} ne 'p') {next;}
# NOTE
            my $t_chaindom = $cur_subsres->{chain2chaindom}->{$t_chain} ;
            my $pcic_chains;
            my $sid2pc ;
            my $dom2pcires ; #domain residues that interact with PC t_chain

            foreach my $t_chain12 (
             keys %{$chainchain_contacts->{$bdp_id}->{$t_chaindom}}) {
               my $t_chain2dom =
                  $chainchain_contacts->{$bdp_id}->{$t_chaindom}->{$t_chain12} ;

               my $t_chain2 = $cur_subsres->{chaindom2chain}->{$t_chain2dom} ;

               if (exists $cur_subsres->{chain2domains}->{$t_chain2}) {
                  $pcic_chains->{$t_chain}->{$t_chain2dom} = $t_chain12; }
            }

# 3. iterate over chain-chain contacts of interest, loading residue pairs
#    and translating to domain definition.
            foreach my $t_chain2dom (keys %{$pcic_chains->{$t_chain}}) {
               my $t_chain12 = $pcic_chains->{$t_chain}->{$t_chain2dom} ;
   
               my ($t_sid1, $t_sid2) = split(/\t/, $t_chain12) ;
               my $interface = _pilig_load_interface_contacts({
                  sid1 => $t_sid1,
                  sid2 => $t_sid2,
                  dist_thresh => $pilig_specs->{PARAM_INTERFACE_DIST_THRESH},
                  fn => $pb->{bdp2contactsfn}->{$bdp_id}
               }) ;

# does this res also have chain (doesnt matter, it is on chain2...)
               foreach my $res (sort
                keys %{$interface->{intres}->{$t_chain2dom}}) {
                  my $t_res = $res ; $t_res =~ s/\n/_/ ;
                  if (!exists $cur_subsres->{res2sid}->{$res}) {
                     next; }
                  my $c_dom = $cur_subsres->{res2sid}->{$res};
                  push @{$dom2pcires->{$c_dom}}, $t_res ;
               }
            }

            foreach my $sid (keys %{$dom2pcires}) {
               my $osid = $pb->{sid2osid}->{$sid} ;
               my $res_string = join("\t", @{$dom2pcires->{$sid}}) ;
               my @outvals = ($bdp_id, $sid, $osid, 'fam',
                              $pb->{sid2class}->{fam}->{$sid},
                              $pb->{chain_2_type}->{$bdp_id}->{$t_chain},
                              $pb->{chain_2_length}->{$bdp_id}->{$t_chain},
                              $t_chain, $t_chaindom,
                              $res_string,
               ) ;
               print join("\t", @outvals)."\n" ;
            }
         }
      }
   }

}


# translate real res to aln positions for pep and nuc bits - ala assign_exp
sub assign_pepnuci {
# both raw and PQS entries qre in the assign_pepnuci_realres file...
#  so make sure that the BDP of interest is actually a raw pdb file before
#  converting to pep and nuc bits

   my $in = shift ;

   my $pilig_specs = set_pilig_specs() ;
   my $pibase_specs = pibase::get_specs() ;

   my $astral = _pilig_astral_preload({
      pibase_specs => $pibase_specs,
      pilig_specs => $pilig_specs,
   }) ;

   if (!-s $pilig_specs->{outfiles}->{assign_pepnuci_realres}) {
      die "ERROR: can't find assign_pepnuci_realres file: ".
         $pilig_specs->{outfiles}->{assign_pepnuci_realres} ;
   }

   my $sid_list ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         chomp $line;
         my ($sid, $bdp, $class) = split(/\t/, $line);
         $sid_list->{sid2bdp}->{$sid} = $bdp ;
         $sid_list->{sid2class}->{$sid} = $class;

         $sid_list->{class2sid}->{fam}->{$class}->{$sid}++ ;
         my ($sf) = ($class =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
         $sid_list->{class2sid}->{sf}->{$sf}->{$sid}++ ;
      }
      close(INF) ;
   } else {
      my $bdp_raw ;
      {
         my ($t2_bdp, $t2_raw) = pibase::rawselect_tod(
            "SELECT bdp_id, raw_pdb FROM bdp_files") ;

         foreach my $j ( 0 .. $#{$t2_bdp}) {
            if ($t2_raw->[$j] == 1) {
               $bdp_raw->{$t2_bdp->[$j]}++ ; } }
      }

      my ($t_sid, $t_bdp, $t_class) = pibase::rawselect_tod(
         "SELECT subset_id, bdp_id, class FROM subsets") ;

      foreach my $j ( 0 .. $#{$t_sid}) {
         if (!exists $bdp_raw->{$t_bdp->[$j]} ||
             $t_sid->[$j] !~ /SCOP/ ||
             $t_class->[$j] !~ /^[a-g]/) {next;}
         $sid_list->{sid2bdp}->{$t_sid->[$j]} = $t_bdp->[$j] ;
         $sid_list->{sid2class}->{$t_sid->[$j]} = $t_class->[$j] ;

         $sid_list->{class2sid}->{fam}->{$t_class->[$j]}->{$t_sid->[$j]}++ ;
         my ($sf) = ($t_class->[$j] =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
         $sid_list->{class2sid}->{sf}->{$sf}->{$t_sid->[$j]}++ ;
      }
   }


   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
# if this is master script, split input, recall self with
# in_fn specified, cluster_fl = 0, import_fl = 0
# send out, cluster run, and return merged

      print "* assign_pepnuci() ".localtime() if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{assign_pepnuci_in}, $temp_fn->{assign_pepnuci_in}) =
         tempfile("splits_assign_pepnuci_input.XXXXX");
      ($temp_fh->{assign_pepnuci_out}, $temp_fn->{assign_pepnuci_out}) =
         tempfile("splits_assign_pepnuci_SGEout_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{assign_pepnuci_out}) ;
      ($temp_fh->{assign_pepnuci_err}, $temp_fn->{assign_pepnuci_err}) =
         tempfile("splits_assign_pepnuci_SGEerr_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{assign_pepnuci_err}) ;

      foreach my $class (sort keys %{$sid_list->{class2sid}->{fam}}) {
         foreach my $sid (sort keys %{$sid_list->{class2sid}->{fam}->{$class}}){
            print {$temp_fh->{assign_pepnuci_in}}
               join("\t", $sid, $sid_list->{sid2bdp}->{$sid}, $class)."\n" ; }}
      close($temp_fh->{assign_pepnuci_in}) ;

      my $split_dir = tempdir("splits_assign_pepnuci.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{assign_pepnuci_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.assign_pepnuci.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::pilig qw/assign_pepnuci/ ;

main() ;

sub main {

         pibase::pilig::assign_pepnuci({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.assign_pepnuci.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.assign_pepnuci.XXXXX");

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
         out_fn => $temp_fn->{assign_pepnuci_out},
         err_fn => $temp_fn->{assign_pepnuci_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{assign_pepnuci_out},
           $temp_fn->{assign_pepnuci_out}) ;
      open(REALOUTF,">".$pilig_specs->{outfiles}->{assign_pepnuci}) ;
      while (my $line = readline($temp_fh->{assign_pepnuci_out})) {
         if ($line =~ /^\#/) {next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{assign_pepnuci_out}) ;
      close(REALOUTF) ;

   } else {
   
      my $pb = _pilig_tod_pibase_preload() ;
      my $sid_pepnuci = _pilig_load_pepnuci_realres({
         fn => $pilig_specs->{outfiles}->{assign_pepnuci_realres},
         sids => $sid_list->{sid2class}
      }) ;

      foreach my $classtype (sort keys %{$sid_list->{class2sid}}) {
         foreach my $class (sort keys %{$sid_list->{class2sid}->{$classtype}}){
            print STDERR "NOW ON $classtype $class (line ".__LINE__.")\n" ;
            my $pos2bs ;
            my $class_aln = pibase::ASTRAL::load_asteroids_aln({
               aln_fn => $pibase_specs->{asteroids}->{$classtype.'_aln'}.
                  '/'.$class.'.fasta_aln' ,
               seq_fn => $pibase_specs->{asteroids}->{$classtype.'_seq'}.
                  '/'.$class.'.fa' ,
               raf => $astral->{raf},
               gdseqh => $astral->{gdseqh},
               seqclcont100 => $astral->{seqcl2cont}->{100},
               seqcl100 => $astral->{seqcl}->{100},
               allchains => $pb->{pdbchains}
            }) ;
   
            foreach my $sid (sort 
               keys %{$sid_list->{class2sid}->{$classtype}->{$class}}) {

               foreach my $target_chaindom (sort keys
                  %{$sid_pepnuci->{sid2pepnucres}->{$sid}}) {

               my $osid = $pb->{sid2osid}->{$sid} ;
               my $pdb = $pb->{sid2pdb}->{$sid} ;
   
               my @talnres = () ;
               my @undefres = () ;
               foreach my $res (
                @{$sid_pepnuci->{sid2pepnucres}->{$sid}->{$target_chaindom}}) {
                  my $alnres = $class_aln->{resno2pos}->{$osid}->{$res} ;
                  if (!defined $alnres) {
                     push @undefres, 'undef';
                  } else {
                     push @talnres, $alnres ;
                  }
               }
   
               my @salnres = ();
               push @salnres, sort {$a <=> $b} @talnres ;
               push @salnres, @undefres ;
               my $alnposstring = join(',', @salnres) ;
   
               my @outvals = ( $pdb, $sid, $osid, $classtype, $class,
               $sid_pepnuci->{target_chain}->{$target_chaindom}->{chain_type},
                               $target_chaindom,
               $sid_pepnuci->{target_chain}->{$target_chaindom}->{chain_length},
                               $alnposstring,
                               $class_aln->{alnlength} ) ;
               print join("\t", @outvals)."\n" ;
               }
            }
         }
      }
   }

}


sub assign_exp {

   my $in = shift ;

   my $pilig_specs = set_pilig_specs() ;
   my $pibase_specs = pibase::get_specs() ;

#old way:
#   my $astral = astral_preload({ fn => $fn }) ;

   my $astral = _pilig_astral_preload({
      pibase_specs => $pibase_specs,
      pilig_specs => $pilig_specs,
   }) ;

   if (!-s $pilig_specs->{outfiles}->{calc_sid_sasa}) {
      die "ERROR: can't find calc_sid_sasa file: ".
         $pilig_specs->{outfiles}->{calc_sid_sasa} ;
   }

   my $sid_list ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         chomp $line;
         my ($sid, $bdp, $class) = split(/\t/, $line);
         $sid_list->{sid2bdp}->{$sid} = $bdp ;
         $sid_list->{sid2class}->{$sid} = $class;

         $sid_list->{class2sid}->{fam}->{$class}->{$sid}++ ;
         my ($sf) = ($class =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
         $sid_list->{class2sid}->{sf}->{$sf}->{$sid}++ ;
      }
      close(INF) ;
   } else {
      my $bdp_raw ;
      {
         my ($t2_bdp, $t2_raw) = pibase::rawselect_tod(
            "SELECT bdp_id, raw_pdb FROM bdp_files") ;

         foreach my $j ( 0 .. $#{$t2_bdp}) {
            if ($t2_raw->[$j] == 1) {
               $bdp_raw->{$t2_bdp->[$j]}++ ; } }
      }

      my ($t_sid, $t_bdp, $t_class) = pibase::rawselect_tod(
         "SELECT subset_id, bdp_id, class FROM subsets") ;

      foreach my $j ( 0 .. $#{$t_sid}) {
         if (!exists $bdp_raw->{$t_bdp->[$j]} ||
             $t_sid->[$j] !~ /SCOP/ ||
             $t_class->[$j] !~ /^[a-g]/) {next;}
         $sid_list->{sid2bdp}->{$t_sid->[$j]} = $t_bdp->[$j] ;
         $sid_list->{sid2class}->{$t_sid->[$j]} = $t_class->[$j] ;

         $sid_list->{class2sid}->{fam}->{$t_class->[$j]}->{$t_sid->[$j]}++ ;
         my ($sf) = ($t_class->[$j] =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
         $sid_list->{class2sid}->{sf}->{$sf}->{$t_sid->[$j]}++ ;
      }
   }


   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
# if this is master script, split input, recall self with
# in_fn specified, cluster_fl = 0, import_fl = 0
# send out, cluster run, and return merged

      print "* assign_exp() ".localtime() if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{assign_exp_in}, $temp_fn->{assign_exp_in}) =
         tempfile("splits_assign_exp_input.XXXXX");
      ($temp_fh->{assign_exp_out}, $temp_fn->{assign_exp_out}) =
         tempfile("splits_assign_exp_SGEout_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{assign_exp_out}) ;
      ($temp_fh->{assign_exp_err}, $temp_fn->{assign_exp_err}) =
         tempfile("splits_assign_exp_SGEerr_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{assign_exp_err}) ;

      foreach my $class (sort keys %{$sid_list->{class2sid}->{fam}}) {
         foreach my $sid (sort keys %{$sid_list->{class2sid}->{fam}->{$class}}){
            print {$temp_fh->{assign_exp_in}}
               join("\t", $sid, $sid_list->{sid2bdp}->{$sid}, $class)."\n" ; }}
      close($temp_fh->{assign_exp_in}) ;

      my $split_dir = tempdir("splits_assign_exp.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{assign_exp_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.assign_exp.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::pilig qw/assign_exp/ ;

main() ;

sub main {

         pibase::pilig::assign_exp({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.assign_exp.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.assign_exp.XXXXX");

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
         out_fn => $temp_fn->{assign_exp_out},
         err_fn => $temp_fn->{assign_exp_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{assign_exp_out},
           $temp_fn->{assign_exp_out}) ;
      open(REALOUTF,">".$pilig_specs->{outfiles}->{assign_exp}) ;
      while (my $line = readline($temp_fh->{assign_exp_out})) {
         if ($line =~ /^\#/) {next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{assign_exp_out}) ;
      close(REALOUTF) ;

   } else {
   
      my $pb = _pilig_tod_pibase_preload() ;
      my $sid_sasa = _pilig_load_sid_sasa({
         fn => $pilig_specs->{outfiles}->{calc_sid_sasa},
         sids => $sid_list->{sid2class}
      }) ;

      foreach my $classtype (sort keys %{$sid_list->{class2sid}}) {
         foreach my $class (sort keys %{$sid_list->{class2sid}->{$classtype}}){
            print STDERR "NOW ON $classtype $class (line ".__LINE__.")\n" ;
            my $pos2bs ;
            my $class_aln = pibase::ASTRAL::load_asteroids_aln({
               aln_fn => $pibase_specs->{asteroids}->{$classtype.'_aln'}.
                  '/'.$class.'.fasta_aln' ,
               seq_fn => $pibase_specs->{asteroids}->{$classtype.'_seq'}.
                  '/'.$class.'.fa' ,
               raf => $astral->{raf},
               gdseqh => $astral->{gdseqh},
               seqclcont100 => $astral->{seqcl2cont}->{100},
               seqcl100 => $astral->{seqcl}->{100},
               allchains => $pb->{pdbchains}
            }) ;
   
            foreach my $sid (sort 
               keys %{$sid_list->{class2sid}->{$classtype}->{$class}}) {
   
               my $osid = $pb->{sid2osid}->{$sid} ;
               my $pdb = $pb->{sid2pdb}->{$sid} ;
   
               my @talnres = () ;
               my @undefres = () ;
               foreach my $res (@{$sid_sasa->{sid2expres}->{$sid}}) {
                  my $alnres = $class_aln->{resno2pos}->{$osid}->{$res} ;
                  if (!defined $alnres) {
   #                  print STDERR "WARNING: $osid ($sid) residue $res not ".
   #                               "found in ASTRAL alignment\n" ;
                     push @undefres, 'undef';
                  } else {
                     push @talnres, $alnres ;
                  }
   #               $pos2bs->{$alnres}->{e}->{$sid}++ ;
               }
   
               my @salnres = ();
               push @salnres, sort {$a <=> $b} @talnres ;
               push @salnres, @undefres ;
               my $alnposstring = join(',', @salnres) ;
   
               my @outvals = ( $pdb, $sid, $osid, $classtype,
                                    $class, "E", $alnposstring,
                                    $class_aln->{alnlength} ) ;
               print join("\t", @outvals)."\n" ;
            }
         }
      }
   }

}


sub calc_sid_sasa {

   my $in = shift ;

   my $pilig_specs = set_pilig_specs() ;
   my $pibase_specs = pibase::get_specs() ;
   my $modeller_bin = $pibase_specs->{binaries}->{modeller} ;

   my $astral = _pilig_astral_preload({
      pibase_specs => $pibase_specs,
      pilig_specs => $pilig_specs,
   }) ;

   my $sid2fn ;
   {
      my ($t_sid, $t_sidfn) = pibase::rawselect_tod(
         "SELECT subset_id, file_path FROM subsets_files") ;

      foreach my $j ( 0 .. $#{$t_sid}) {
         $sid2fn->{$t_sid->[$j]} = $t_sidfn->[$j];}
   }

   my $sid_list ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         chomp $line;
         my ($sid, $bdp, $class) = split(/\t/, $line);
         push @{$sid_list->{sids}}, $sid ;
         $sid_list->{sid2bdp}->{$sid} = $bdp ;
         $sid_list->{sid2class}->{$sid} = $class;
      }
      close(INF) ;
   } else {
      my $bdp_raw ;
      {
         my ($t2_bdp, $t2_raw) = pibase::rawselect_tod(
            "SELECT bdp_id, raw_pdb FROM bdp_files") ;

         foreach my $j ( 0 .. $#{$t2_bdp}) {
            if ($t2_raw->[$j] == 1) {
               $bdp_raw->{$t2_bdp->[$j]}++ ; } }
      }

      my ($t_sid, $t_bdp, $t_class) = pibase::rawselect_tod(
         "SELECT subset_id, bdp_id, class FROM subsets") ;

      foreach my $j ( 0 .. $#{$t_sid}) {
         if (!exists $bdp_raw->{$t_bdp->[$j]} ||
             $t_sid->[$j] !~ /SCOP/ ||
             $t_class->[$j] !~ /^[a-g]/) {next;}
         $sid_list->{sid2bdp}->{$t_sid->[$j]} = $t_bdp->[$j] ;
         $sid_list->{sid2class}->{$t_sid->[$j]} = $t_class->[$j] ;
         push @{$sid_list->{sids}}, $t_sid->[$j] ;
      }
   }

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
# if this is master script, split input, recall self with
# in_fn specified, cluster_fl = 0, import_fl = 0
# send out, cluster run, and return merged

      print "* calc_sid_sasa() ".localtime() if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{calc_sid_sasa_in}, $temp_fn->{calc_sid_sasa_in}) =
         tempfile("calc_sid_sasa_input.XXXXX");
      ($temp_fh->{calc_sid_sasa_out}, $temp_fn->{calc_sid_sasa_out}) =
         tempfile("splits_calc_sid_sasa_SGEout_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{calc_sid_sasa_out}) ;
      ($temp_fh->{calc_sid_sasa_err}, $temp_fn->{calc_sid_sasa_err}) =
         tempfile("splits_calc_sid_sasa_SGEerr_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{calc_sid_sasa_err}) ;

      foreach my $sid (sort keys %{$sid_list->{sid2bdp}}) {
         print {$temp_fh->{calc_sid_sasa_in}} join("\t", $sid,
            $sid_list->{sid2bdp}->{$sid},
            $sid_list->{sid2class}->{$sid})."\n";}
      close($temp_fh->{calc_sid_sasa_in}) ;

      my $split_dir = tempdir("splits_calc_sid_sasa.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{calc_sid_sasa_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.calc_sid_sasa.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::pilig qw/calc_sid_sasa/ ;

main() ;

sub main {

         pibase::pilig::calc_sid_sasa({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.calc_sid_sasa.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.calc_sid_sasa.XXXXX");

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
         out_fn => $temp_fn->{calc_sid_sasa_out},
         err_fn => $temp_fn->{calc_sid_sasa_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{calc_sid_sasa_out},
           $temp_fn->{calc_sid_sasa_out}) ;
      open(REALOUTF,">".$pilig_specs->{outfiles}->{calc_sid_sasa}) ;
      while (my $line = readline($temp_fh->{calc_sid_sasa_out})) {
         if ($line =~ /^\#/) {next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{calc_sid_sasa_out}) ;
      close(REALOUTF) ;

   } else {

      foreach my $sid (@{$sid_list->{sids}}) {
         print STDERR "NOW ON: $sid\n" ;
         my $exposed ;
   
         my $sid_fn = $sid2fn->{$sid} ;
   
         my $subset_sasa = pibase::modeller::calc_sasa({
            pdb_fn => $sid_fn,
            surftyp => 2,
            modeller_bin => $modeller_bin}) ; 

         if ($#{$subset_sasa->{error_fl}} >= 0 ) {
            foreach my $j ( 0 .. $#{$subset_sasa->{error_fl}}) {
               print STDERR "ERROR: subset $sid: ".
                            "calc_sasa() $subset_sasa->{error_fl}->[$j]\n"; }
            next;
         }
   
         foreach my $j ( 0 .. $#{$subset_sasa->{res_sasa}->{resno}}) {
            my $resno = $subset_sasa->{res_sasa}->{resno}->[$j] ;
            my $chain = $subset_sasa->{res_sasa}->{chain}->[$j] ;
            my $ressig = $resno."_".$chain ;
   
            if ($subset_sasa->{res_sasa}->{sc_perc}->[$j] >
                $pilig_specs->{PARAM_SCPERCACC_THRESH}) {
               push @{$exposed}, $ressig ; }
         }
   
         print join("\t", $sid, $sid_list->{sid2bdp}->{$sid},
                     $sid_list->{sid2class}->{$sid}, @{$exposed})."\n" ;
      }

   }

}

sub _pilig_run_prep {
   my $in = shift ;

   my $fn = set_locations() ;
   my $dbh = connect_pibase_ligbase() ;
   my $astral = astral_preload({fn => $fn});

   my $splits_dir = tempdir($fn->{pilig}->{prep_dir}."/assign_splits.XXXXX") ;
   my $sgeout_dir_lb = tempdir($fn->{pilig}->{assign_lig_dir}."/SGEassign_lig.XXXXX");
   my $sgeout_dir_pb = tempdir($fn->{pilig}->{assign_pi_dir}."/SGEassign_pi.XXXXX");


   print STDERR "preparing LIGBASE data: " ;
   my $lb = ligbase_preload({dbh => $dbh->{lig}}) ;
   my $lb_sge_fn = "assign_lig.$$.SGE.sh" ;
   split_ligs({lb => $lb,
               splits_dir => $splits_dir,
               out_fn_prefix => "ligbase_active_res",
               SGE_dir => $sgeout_dir_lb,
               SGE_fn => $lb_sge_fn }) ;
   print STDERR "X\n" ;

   print STDERR "preparing PIBASE data: " ;
   my $pb = pibase_preload({dbh => $dbh->{pi}, astral => $astral}) ;
   my $pb_sge_fn = "assign_pi.$$.SGE.sh" ;
   split_pi({  pb => $pb,
               splits_dir => $splits_dir,
               SGE_dir => $sgeout_dir_pb,
               out_fn_prefix => "bdp_ids",
               SGE_fn => $pb_sge_fn}) ;
   print STDERR "X\n" ;

   print STDERR "SGE scripts:\nqsub $pb_sge_fn\nqsub $lb_sge_fn\n ";
}



sub _pilig_astral_preload {

   my $in = shift;
   my $pibase_specs = $in->{pibase_specs} ;

#   print STDERR "Load ASTRAL fasta headers (gdseq): ".localtime() ;
   my $astral = pibase::ASTRAL::load_astral_headers({
      fn => $pibase_specs->{astral}->{gd_seq} }) ;
#   print STDERR ", DONE: ".localtime()."\t" ;
#   print STDERR " cur size ".total_size($astral)." bytes\n";
   

   my $astral_classes = pibase::ASTRAL::get_astral_classlist({
      pibase_specs => $pibase_specs}) ;
   $astral->{classes} = $astral_classes ;

#   print STDERR "Load ASTRAL raf: ".localtime() ;
   my $raf = pibase::ASTRAL::raf_preload({
      fn => $pibase_specs->{astral}->{raf}
   }) ;
   $astral->{raf} = $raf ;
#   print STDERR ", DONE: ".localtime()."\t" ;
#   print STDERR " cur size ".total_size($astral)." bytes\n";

#   print STDERR "Load ASTRAL sequence clusters: ".localtime() ;
   pibase::ASTRAL::load_astral_clusters({
      pibase_specs => $pibase_specs,
      out => $astral
   }) ;
#   print STDERR ", DONE: ".localtime()."\t" ;
#   print STDERR " cur size ".total_size($astral)." bytes\n";

   return $astral ;

}

sub _pilig_set_locations {

   my $fn ;
   $fn->{astraldir} = "/alto2/home/fred/sali/projects/domain_interfaces/work/PIBASE0510/data/preexisting/ASTRAL/" ;

   $fn->{astral_raf} = $fn->{astraldir}.
      "/astral-rapid-access-1.69.raf" ;

   $fn->{gdseq}= $fn->{astraldir}.'/scopseq-1.69/'.
      'astral-scopdom-seqres-gd-all-1.69.fa' ;

   foreach my $seqid ( qw/100 95 90 70 50 30 20 10/) {
      $fn->{seqcl}->{$seqid} = $fn->{astraldir}.'/scopseq-1.69/'.
      '/astral-scopdom-seqres-gd-sel-gs-bib-verbose-'.$seqid.'-1.69.txt' ;}

   $fn->{aln}->{fam} = $fn->{astraldir}."/aln/fam/" ;
   $fn->{aln}->{sf} = $fn->{astraldir}."/aln/sf/" ;

   $fn->{pilig}->{rootdir} = "/netapp/home/fred/pilig" ;
   $fn->{pilig}->{park_rootdir} = "/park1/fred/pilig" ;
   $fn->{pilig}->{prep_dir} = $fn->{pilig}->{rootdir}."/prepin";
   $fn->{pilig}->{assign_lig_dir} = $fn->{pilig}->{rootdir}."/assign_lig";
   $fn->{pilig}->{assign_pi_dir} = $fn->{pilig}->{rootdir}."/assign_pi";

   $fn->{pilig}->{assign_lig_fn} = $fn->{pilig}->{rootdir}."/assign_lig.out";
   $fn->{pilig}->{assign_pi_fn} = $fn->{pilig}->{rootdir}."/assign_pi.out";
   $fn->{pilig}->{assign_exp_fn} = $fn->{pilig}->{rootdir}."/assign_exp.out";

   $fn->{pilig}->{park_assign_lig_dir}=$fn->{pilig}->{park_rootdir}."/assign_lig";
   $fn->{pilig}->{park_assign_pi_dir}=$fn->{pilig}->{park_rootdir}."/assign_pi";

   $fn->{pilig}->{park_assign_lig_fn}=$fn->{pilig}->{park_rootdir}."/assign_lig.out";
   $fn->{pilig}->{park_assign_pi_fn}=$fn->{pilig}->{park_rootdir}."/assign_pi.out";
   $fn->{pilig}->{park_assign_exp_fn}=$fn->{pilig}->{park_rootdir}."/assign_exp.out";

   
   $fn->{msdchem}->{liginfo} =$fn->{pilig}->{park_rootdir}."/aux/ligandinfo.msdchem.txt" ;

   return $fn ;

}

sub _pilig_connect_pibase_ligbase {

   require DBI ;

   my $dbh ;

   ($dbh->{lig}) = pibase::connect_pibase({
      'db' => 'ligbase_new',
      'host' => 'modbase',
      'user' => 'modbase',
      'pass' => 'modbasesecret'
   }) ;

   ($dbh->{pi}) = pibase::connect_pibase() ;

   return $dbh ;
}


sub _pilig_compare_residue_sets {

   my $in = shift;;
   my $aln = $in->{aln} ;
   my $sid1 = $in->{sid1} ;
   my $sid2 = $in->{sid2} ;
   my $res1 = $in->{res1} ;
   my $res2 = $in->{res2} ;

   my $intpos ;
   foreach my $res1 (keys %{$res1}) {
      my $pos = $aln->{resno2pos}->{$sid1}->{$res1} ;
      if (!defined $pos) {
         print STDERR "ERROR: compare_residue_sets(): $res1 from $sid1 undefined pos; aborting $sid1 -- $sid2 comparison\n" ; return 0;
      } else {
         $intpos->{$pos}++ ;
      }}


   foreach my $res2 (keys %{$res2}) {
      my $pos = $aln->{resno2pos}->{$sid2}->{$res2} ;
      if (!defined $pos) {
         print STDERR "ERROR: compare_residue_sets(): $res2 from $sid2 undefined pos; aborting $sid1 -- $sid2 comparison\n" ; return 0;
      } else {
         $intpos->{$pos}++ ;
      }}

   my $union = 0 ;
   my $diff = 0 ;
   foreach my $pos (keys %{$intpos}) {
      if ($intpos->{$pos} == 2) {
         $union++;
      } else {
         $diff++ ; }
   }

   if (($union + $diff) > 0) {
      return ($union / ($union + $diff)) ;
   } else {
      return 0 ;
   }

}


sub _pilig_load_interface_contacts {

   my $in = shift ;
   my $sid1 = $in->{sid1} ;
   my $sid2 = $in->{sid2} ;
   my $fn = $in->{fn} ;
   my $pilig_specs = set_pilig_specs() ;

   if (!exists $in->{dist_thresh}) {
      $in->{dist_thresh} = $pilig_specs->{PARAM_INTERFACE_DIST_THRESH} ; }

   my ( $subset_id_1, $subset_id_2,
        $chain_id_1, $resno_1, $chain_id_2, $resno_2, $min_dist ) =
       pibase::rawselect_metatod($fn, "SELECT subset_id_1, subset_id_2, chain_id_1, resno_1, chain_id_2, resno_2, min_dist FROM $fn") ;

   my $contacts ;
   my $intres ;
   foreach my $j ( 0 .. $#{$subset_id_1}) {
      if ($subset_id_1->[$j] ne $sid1) {next;}
      if ($subset_id_2->[$j] ne $sid2) {next;}
      if ($min_dist->[$j] > $in->{dist_thresh}) {next;}

#      print STDERR "$subset_id_1->[$j] -- $subset_id_2->[$j]: $resno_1->[$j] $chain_id_1->[$j] to $resno_2->[$j] $chain_id_2->[$j]\n" ;

      $intres->{$sid1}->{$resno_1->[$j]."\n".$chain_id_1->[$j]}++ ;
      $intres->{$sid2}->{$resno_2->[$j]."\n".$chain_id_2->[$j]}++ ;
      $contacts->{$resno_1->[$j]."\n".$chain_id_1->[$j]}->{$resno_2->[$j]."\n".$chain_id_2->[$j]}++ ;
   }

   return {
      intres => $intres,
      contacts => $contacts
   } ;

}


sub _pilig_ligbase_preload {
   my $in = shift ;
   my $dbh = $in->{dbh} ;

   my $query="SELECT pdb_code, a.ligand_code, a.ligand_idnum, a.ligand_chain, ".
         "a.residue_num, a.ligand_chain FROM active_residues as a";

   my ($pdb, $ligcod, $ligid, $ligchain, $resno, $pdbchain) =
      pibase::mysql_fetchcols($dbh,$query) ;

   my ($pdb2ligid, $pdb2res2ligid) ;
   foreach my $j ( 0 .. $#{$ligcod}) {
      if ($pdbchain->[$j] eq '') {
         $pdbchain->[$j] = ' '; }

      $pdb2ligid->{$pdb->[$j]}->{$ligcod->[$j]."\n".$ligid->[$j]}++ ;
      $pdb2res2ligid->{$pdb->[$j]}->{$resno->[$j]."\n".$pdbchain->[$j]}->{$ligcod->[$j]."\n".$ligid->[$j]}++ ;
   }

   my $lb = {
      pdb2ligid => $pdb2ligid,
      pdb2res2ligid => $pdb2res2ligid
   } ;

   return $lb ;
}


sub _pilig_tod_ligbase_preload {

   my $in = shift ;
   my $pilig_specs = set_pilig_specs() ;
   my $binaries = pibase::locate_binaries()  ;
   my ($pdb2ligid, $pdb2res2ligid) ;

#   my $query="SELECT pdb_code, a.ligand_code, a.ligand_idnum, a.ligand_chain,".
#         " a.residue_num, a.ligand_chain FROM active_residues as a";
#
#   my ($pdb, $ligcod, $ligid, $ligchain, $resno, $pdbchain) =
#      pibase::mysql_fetchcols($dbh,$query) ;
#
#   foreach my $j ( 0 .. $#{$ligcod}) {
#      if ($pdbchain->[$j] eq '') {
#         $pdbchain->[$j] = ' '; }
#
#      $pdb2ligid->{$pdb->[$j]}->{$ligcod->[$j]."\n".$ligid->[$j]}++ ;
#      $pdb2res2ligid->{$pdb->[$j]}->{$resno->[$j]."\n".$pdbchain->[$j]}->{$ligcod->[$j]."\n".$ligid->[$j]}++ ;
#   }

   my $data_fn = $pilig_specs->{ligbase}->{activeres} ;
#   my ($t_fh, $t_fn) = tempfile("unc_activeres.XXXXX", SUFFIX=>".out") ;
   my ($t_fh, $t_fn) = tempfile() ; close($t_fh) ; 
#      $data_fn = $binaries->{zcat}." $data_fn |" ;
   my $activeres_fh ;
   if ($data_fn =~ /gz$/) {
      my $tcom = $binaries->{zcat}." $data_fn > $t_fn" ;
#      print STDERR "TCOM: $tcom\n" ;
      system($tcom) ;
      open($activeres_fh, $t_fn) ;
   } else {
      open($activeres_fh, $data_fn) ;
   }

   my $line = readline($activeres_fh);
   my $activeres_f2i ;
   if ($line !~ /^\#/) {
      die "ERROR: corrupted ligbase.active_res data: ".
         $pilig_specs->{ligbase}->{activeres};
   } else {
      $line =~ s/^\#// ;
      my @t = split(/\t/, $line) ;
      foreach my $j ( 0 .. $#t) {
         $activeres_f2i->{$t[$j]} = $j ; }
   }

   while (my $line = readline($activeres_fh)) {
      chomp $line;

      my @t = split(/\t/, $line) ;
      my $pdbchain = $t[$activeres_f2i->{ligand_chain}] ;
      if ($pdbchain eq '') { $pdbchain = ' '; }
      my $pdb = $t[$activeres_f2i->{pdb_code}] ;
      my $ligcod = $t[$activeres_f2i->{ligand_code}] ;
      my $ligid = $t[$activeres_f2i->{ligand_idnum}] ;
      my $resno = $t[$activeres_f2i->{residue_num}] ;

      $pdb2ligid->{$pdb}->{$ligcod."\n".$ligid}++ ;
      $pdb2res2ligid->{$pdb}->{$resno."\n".$pdbchain}->{$ligcod."\n".$ligid}++ ;
   }
   close($activeres_fh) ;


   my $lb = {
      pdb2ligid => $pdb2ligid,
      pdb2res2ligid => $pdb2res2ligid
   } ;


   return $lb ;
}


#fpd081226_2001  edit to only load a single SCOP class; ease memory reqs
sub _pilig_tod_ligbase_preload_loadperclass {

   my $in = shift ;
   my $pilig_specs = set_pilig_specs() ;
   my $binaries = pibase::locate_binaries()  ;
   my ($pdb2ligid, $pdb2res2ligid) ;

   my $classtype = $in->{classtype} ;
   my $class = $in->{class} ;
   my $class2pdb = $in->{class2pdb} ;

#   my $query="SELECT pdb_code, a.ligand_code, a.ligand_idnum, a.ligand_chain,".
#         " a.residue_num, a.ligand_chain FROM active_residues as a";
#
#   my ($pdb, $ligcod, $ligid, $ligchain, $resno, $pdbchain) =
#      pibase::mysql_fetchcols($dbh,$query) ;
#
#   foreach my $j ( 0 .. $#{$ligcod}) {
#      if ($pdbchain->[$j] eq '') {
#         $pdbchain->[$j] = ' '; }
#
#      $pdb2ligid->{$pdb->[$j]}->{$ligcod->[$j]."\n".$ligid->[$j]}++ ;
#      $pdb2res2ligid->{$pdb->[$j]}->{$resno->[$j]."\n".$pdbchain->[$j]}->{$ligcod->[$j]."\n".$ligid->[$j]}++ ;
#   }

   my $data_fn = $pilig_specs->{ligbase}->{activeres} ;
   my ($t_fh, $t_fn) = tempfile() ; close($t_fh) ; 
   my $activeres_fh ;
   if ($data_fn =~ /gz$/) {
      my $tcom = $binaries->{zcat}." $data_fn > $t_fn" ;
      system($tcom) ;
      open($activeres_fh, $t_fn) ;
   } else {
      open($activeres_fh, $data_fn) ;
   }

   my $line = readline($activeres_fh);
   my $activeres_f2i ;
   if ($line !~ /^\#/) {
      die "ERROR: corrupted ligbase.active_res data: ".
         $pilig_specs->{ligbase}->{activeres};
   } else {
      $line =~ s/^\#// ;
      my @t = split(/\t/, $line) ;
      foreach my $j ( 0 .. $#t) {
         $activeres_f2i->{$t[$j]} = $j ; }
   }

   while (my $line = readline($activeres_fh)) {
      chomp $line;

      my @t = split(/\t/, $line) ;
      my $pdbchain = $t[$activeres_f2i->{ligand_chain}] ;
      if ($pdbchain eq '') { $pdbchain = ' '; }
      my $pdb = $t[$activeres_f2i->{pdb_code}] ;

      if (!exists $class2pdb->{$classtype}->{$class}->{$pdb}) {next;}

      my $ligcod = $t[$activeres_f2i->{ligand_code}] ;
      my $ligid = $t[$activeres_f2i->{ligand_idnum}] ;
      my $resno = $t[$activeres_f2i->{residue_num}] ;

      $pdb2ligid->{$pdb}->{$ligcod."\n".$ligid}++ ;
      $pdb2res2ligid->{$pdb}->{$resno."\n".$pdbchain}->{$ligcod."\n".$ligid}++ ;
   }
   close($activeres_fh) ;


   my $lb = {
      pdb2ligid => $pdb2ligid,
      pdb2res2ligid => $pdb2res2ligid
   } ;

   return $lb ;
}

sub _pilig_tod_pibase_preload {

   my $in = shift ;
   my $astral = $in->{astral} ;
   my $pilig_specs = set_pilig_specs() ;

   my ($bdp_2_raw, $pdb2bdp, $bdp2pdb) = pibase::todload_bdp_ids(
      "bdp_id_2_raw_pdb", "pdb_id_2_bdp_id", "bdp_id_2_pdb_id") ;

   my $bdp2contactsfn ;
   {
      my ($t_bdp, $t_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM interface_contacts_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2contactsfn->{$t_bdp->[$j]} = $t_fn->[$j] ; }
   }

   my $subsetsource_name2id ;
   {
      my ($t_ssid, $t_name) = pibase::rawselect_tod(
         "SELECT subset_source_id, subset_source FROM subsets_source") ;
      foreach my $j ( 0 .. $#{$t_ssid}) {
         $subsetsource_name2id->{$t_name->[$j]} = $t_ssid->[$j]  ; }
   }

   my $bdp2subsetsresfn ;
   {
      my ($t_bdp, $t_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM subsets_residues_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2subsetsresfn->{$t_bdp->[$j]} = $t_fn->[$j] ; }
   }

   my $sid2class ;
   my $osid2class ;
   my $sid2osid ;
   my $sid2pdb ;
   my $class2pdb ;
   {
      my ($sid, $class, $bdp_id, $ssid) =pibase::rawselect_tod(
         "SELECT subset_id, class, bdp_id, subset_source_id FROM subsets") ;
      foreach my $j ( 0 .. $#{$sid}) {
         if ($ssid->[$j] != $subsetsource_name2id->{'scop'}) { next;}

         if (!defined $bdp_id->[$j] || $bdp_id->[$j] eq '') {next;}
         $sid2pdb->{$sid->[$j]} = $bdp2pdb->{$bdp_id->[$j]} ;

         $sid2class->{fam}->{$sid->[$j]} = $class->[$j] ;
         my $fam = $class->[$j] ;
         my ($sf) = ($fam =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
         $sid2class->{sf}->{$sid->[$j]} = $sf ;

         if ( $bdp_id->[$j] ne 'NULL') {
            $class2pdb->{sf}->{$sf}->{$bdp2pdb->{$bdp_id->[$j]}}++ ;
            $class2pdb->{fam}->{$class->[$j]}->{$bdp2pdb->{$bdp_id->[$j]}}++ ;
         }

         my $osid = substr($sid->[$j], -7, 7) ;
         if (exists $astral->{gdom}->{$osid}) {
            $osid = $astral->{gdom}->{$osid} ; }
         $sid2osid->{$sid->[$j]} = $osid ;
         $osid2class->{fam}->{$osid} = $class->[$j] ;
      }
   }


   my $classes2int ;
   my ($bdp_id, $sid1, $sid2, $class1, $class2, $abchains) ; 
   {
      my ($t_bdp, $t_sid1, $t_sid2, $t_class1, $t_class2, $t_numcon, $t_chains)=
      pibase::rawselect_tod(
         "SELECT bdp_id, subset_id_1, subset_id_2, class_1, class_2, ".
         " num_contacts, chains FROM intersubset_contacts") ;

## TESTING FUCKER fpd080630_0210 
#   my ($t_subset_id, $t_chain_id, $t_resno_serial, $t_resno) =
#      pibase::rawselect_metatod("/groups/eddy/home/davisf/work/pibase/pibase200708/data/metatod/subsets_residues/32/subsets_residues_32533.32533.pibase.gz",
#      "SELECT subset_id, chain_id, resno_serial, resno FROM /groups/eddy/home/davisf/work/pibase/pibase200708/data/metatod/subsets_residues/32/subsets_residues_32533.32533.pibase.gz") ;
#   foreach my $j ( 0 .. $#{$t_subset_id}) {
#      print STDERR join("\t",$j,$t_subset_id->[$j], $t_chain_id->[$j], $t_resno->[$j])."\n" ;
#   }
#   die ;
# ---TESTING FUCKER---


      foreach my $j ( 0 .. $#{$t_bdp}) {
         if ($t_sid1->[$j] !~ /SCOP/) {next;}
         if ($t_numcon->[$j] < $pilig_specs->{PARAM_MIN_NUMCONTACTS}) {next;}
         if ($bdp_2_raw->{$t_bdp->[$j]} != 1) {next;}

         push @{$bdp_id}, $t_bdp->[$j] ;
         push @{$sid1}, $t_sid1->[$j] ;
         push @{$sid2}, $t_sid2->[$j] ;
         push @{$class1}, $t_class1->[$j] ;
         push @{$class2}, $t_class2->[$j] ;
         push @{$abchains}, $t_chains->[$j] ;
      }
   }

   my $bdp2sid12 ;
   my ($interacting_classes, $class2sidpairs, $sid12chain) ;
   my ($osid1, $osid2, $revfl) ; #revfl if fam2 <fam1 and flipped to get fam21 isntead of fam12, let it be known so later it can be dealt with


   foreach my $j ( 0 .. $#{$bdp_id}) {
# ONLY ONE OF THEM HAS TO BE OF PROPER CLASS
#      if ($class1->[$j] !~ /^[a-g]/ || $class2->[$j] !~ /^[a-g]/) {next;}
      if ($class1->[$j] !~ /^[a-g]/ && $class2->[$j] !~ /^[a-g]/) {next;}
      $bdp2sid12->{$bdp_id->[$j]}->{$sid1->[$j]."\t".$sid2->[$j]}++ ;

      my $sid12 = $sid1->[$j]."\t".$sid2->[$j] ;
      $sid12chain->{$sid12} = $abchains->[$j] ;

      my $fam1 = $class1->[$j];
      my $fam2 = $class2->[$j] ;

      my ($sf1) = ($fam1 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
      my ($sf2) = ($fam2 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;

      $osid1->[$j] = substr($sid1->[$j], -7, 7) ;
      $osid2->[$j] = substr($sid2->[$j], -7, 7) ;

      if (exists $astral->{gdom}->{$osid1->[$j]}) {
         $osid1->[$j] = $astral->{gdom}->{$osid1->[$j]} ; }

      if (exists $astral->{gdom}->{$osid2->[$j]}) {
         $osid2->[$j] = $astral->{gdom}->{$osid2->[$j]} ; }

      my $fam12;
      if ($class2->[$j] lt $class1->[$j]) {
         $revfl->[$j] = 1 ;
         $fam12 = $class2->[$j]."\t".$class1->[$j];
      } else {
         $revfl->[$j] = 0 ;
         $fam12 = $class1->[$j]."\t".$class2->[$j];
      }
      $interacting_classes->{fam}->{$fam1}->{$sid1->[$j]}->{$sid2->[$j]} = $j;
      $interacting_classes->{fam}->{$fam2}->{$sid2->[$j]}->{$sid1->[$j]} = $j;
      $interacting_classes->{sf}->{$sf1}->{$sid1->[$j]}->{$sid2->[$j]} = $j;
      $interacting_classes->{sf}->{$sf2}->{$sid2->[$j]}->{$sid1->[$j]} = $j;

      $class2sidpairs->{fam}->{$fam1}->{$j}++ ;
      $class2sidpairs->{fam}->{$fam2}->{$j}++ ;
      $class2sidpairs->{sf}->{$sf1}->{$j}++ ;
      $class2sidpairs->{sf}->{$sf2}->{$j}++ ;
   }

   my $chain_2_pdbchain ;
   my ($chain_2_start, $chain_2_end) ;
   my $chain_2_length ; my $chain_2_type ;
   my $pdbchains ;
   {
      my ($tbdp_id, $real_chain_id, $pdb_chain_id, $resno_start, $resno_end) ;
      {
         my ($t_bdp, $t_real_cid, $t_pdb_cid, $t_start, $t_end,
             $t_type, $t_numres, $t_chain_no) = pibase::rawselect_tod(
            "SELECT bdp_id, real_chain_id, pdb_chain_id, ".
            "start_resno, end_resno, chain_type, num_res, ".
            "real_chain_no FROM bdp_chains") ;
    
         foreach my $j ( 0 .. $#{$t_bdp}) {
            $chain_2_length->{$t_bdp->[$j]}->{$t_chain_no->[$j]} =
               $t_numres->[$j];
            $chain_2_type->{$t_bdp->[$j]}->{$t_chain_no->[$j]} =
               $t_type->[$j] ;

            if ($bdp_2_raw->{$t_bdp->[$j]} != 1) {next;}
            if (!defined $t_pdb_cid->[$j] || $t_pdb_cid->[$j] eq '') {next;}

            if ($t_type->[$j] ne 'p') {next;}
            push @{$tbdp_id}, $t_bdp->[$j] ;
            push @{$real_chain_id}, $t_real_cid->[$j] ;
            push @{$pdb_chain_id}, $t_pdb_cid->[$j] ;
            push @{$resno_start}, $t_start->[$j] ;
            push @{$resno_end}, $t_end->[$j] ;
         }
      }

#      print STDERR ($#{$tbdp_id} + 1)." chains loaded entries\n " ;

      foreach my $j ( 0 .. $#{$tbdp_id}) {
         $pdbchains->{$bdp2pdb->{$tbdp_id->[$j]}}->{$real_chain_id->[$j]}++;

         $chain_2_pdbchain->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
            $pdb_chain_id->[$j] ;
         $chain_2_start->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
            $resno_start->[$j] ;
         $chain_2_end->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
            $resno_end->[$j] ;
      }
   }

   my $pb = {
      bdp2contactsfn => $bdp2contactsfn,
      bdp2subsetsresfn => $bdp2subsetsresfn,
      bdp2pdb => $bdp2pdb,
      pdb2bdp => $pdb2bdp,
      class2pdb => $class2pdb,
      bdp_id => $bdp_id,
      sid2pdb => $sid2pdb,
      sid2class => $sid2class,
      osid2class => $osid2class,
      sid2osid => $sid2osid,
      sid1 => $sid1,
      sid2 => $sid2,
      class1 => $class1,
      class2 => $class2,
      osid1 => $osid1,
      osid2 => $osid2,
      revfl => $revfl,
      interacting_classes => $interacting_classes,
      class2sidpairs => $class2sidpairs,
      classes2int => $classes2int,
      pdbchains => $pdbchains,
      chain_2_type => $chain_2_type,
      chain_2_length => $chain_2_length,
      chain_2_pdbchain => $chain_2_pdbchain,
      chain_2_start => $chain_2_start,
      chain_2_end => $chain_2_end,
      bdp2sid12 => $bdp2sid12,
      sid12chain => $sid12chain
   } ;

   if (exists $in->{onlythese}) {
      my $realpb ;
      foreach my $type (keys %{$in->{onlythese}}) {
         $realpb->{$type} = $pb->{$type}; }
      return $realpb ;
   } else {
      return $pb ;
   }

}


sub RAWCONVERT_PRE080630_pilig_tod_pibase_preload {

   my $in = shift ;
   my $astral = $in->{astral} ;
   my $pilig_specs = set_pilig_specs() ;

   my ($bdp_2_raw, $pdb2bdp, $bdp2pdb) = pibase::todload_bdp_ids(
      "bdp_id_2_raw_pdb", "pdb_id_2_bdp_id", "bdp_id_2_pdb_id") ;

   my $bdp2contactsfn ;
   {
      my ($t_bdp, $t_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM interface_contacts_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2contactsfn->{$t_bdp->[$j]} = $t_fn->[$j] ; }
   }

   my $subsetsource_name2id ;
   {
      my ($t_ssid, $t_name) = pibase::rawselect_tod(
         "SELECT subset_source_id, subset_source FROM subsets_source") ;
      foreach my $j ( 0 .. $#{$t_ssid}) {
         $subsetsource_name2id->{$t_name->[$j]} = $t_ssid->[$j]  ; }
   }

   my $bdp2subsetsresfn ;
   {
      my ($t_bdp, $t_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM subsets_residues_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2subsetsresfn->{$t_bdp->[$j]} = $t_fn->[$j] ; }
   }

   my $sid2class ;
   my $osid2class ;
   my $sid2osid ;
   my $sid2pdb ;
   {
      my ($sid, $class, $bdp_id, $ssid) =pibase::rawselect_tod(
         "SELECT subset_id, class, bdp_id, subset_source_id FROM subsets") ;
      foreach my $j ( 0 .. $#{$sid}) {
         if ($ssid->[$j] != $subsetsource_name2id->{'scop'}) { next;}

         if (!defined $bdp_id->[$j] || $bdp_id->[$j] eq '') {next;}
         $sid2pdb->{$sid->[$j]} = $bdp2pdb->{$bdp_id->[$j]} ;

         $sid2class->{fam}->{$sid->[$j]} = $class->[$j] ;
         my $fam = $class->[$j] ;
         my ($sf) = ($fam =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
         $sid2class->{sf}->{$sid->[$j]} = $sf ;

         my $osid = substr($sid->[$j], -7, 7) ;
         if (exists $astral->{gdom}->{$osid}) {
            $osid = $astral->{gdom}->{$osid} ; }
         $sid2osid->{$sid->[$j]} = $osid ;
         $osid2class->{fam}->{$osid} = $class->[$j] ;
      }
   }


   my ($bdp_id, $sid1, $sid2, $class1, $class2, $abchains) ; 
   {
      my ($t_bdp, $t_sid1, $t_sid2, $t_class1, $t_class2, $t_numcon, $t_chains)=
      pibase::rawselect_tod(
         "SELECT bdp_id, subset_id_1, subset_id_2, class_1, class_2, ".
         " num_contacts, chains FROM intersubset_contacts") ;

      foreach my $j ( 0 .. $#{$t_bdp}) {
         if ($t_sid1->[$j] !~ /SCOP/) {next;}
         if ($t_numcon->[$j] < $pilig_specs->{PARAM_MIN_NUMCONTACTS}) {next;}
         if ($bdp_2_raw->{$t_bdp->[$j]} != 1) {next;}
         push @{$bdp_id}, $t_bdp->[$j] ;
         push @{$sid1}, $t_sid1->[$j] ;
         push @{$sid2}, $t_sid2->[$j] ;
         push @{$class1}, $t_class1->[$j] ;
         push @{$class2}, $t_class2->[$j] ;
         push @{$abchains}, $t_chains->[$j] ;
      }
   }

   my $bdp2sid12 ;
   my ($interacting_classes, $class2sidpairs, $sid12chain) ;
   my ($osid1, $osid2, $revfl) ; #revfl if fam2 <fam1 and flipped to get fam21 isntead of fam12, let it be known so later it can be dealt with

   foreach my $j ( 0 .. $#{$bdp_id}) {
# ONLY ONE OF THEM HAS TO BE OF PROPER CLASS
#      if ($class1->[$j] !~ /^[a-g]/ || $class2->[$j] !~ /^[a-g]/) {next;}
      if ($class1->[$j] !~ /^[a-g]/ && $class2->[$j] !~ /^[a-g]/) {next;}
      $bdp2sid12->{$bdp_id->[$j]}->{$sid1->[$j]."\t".$sid2->[$j]}++ ;

      my $sid12 = $sid1->[$j]."\t".$sid2->[$j] ;
      $sid12chain->{$sid12} = $abchains->[$j] ;

      my $fam1 = $class1->[$j];
      my $fam2 = $class2->[$j] ;

      my ($sf1) = ($fam1 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
      my ($sf2) = ($fam2 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;

      $osid1->[$j] = substr($sid1->[$j], -7, 7) ;
      $osid2->[$j] = substr($sid2->[$j], -7, 7) ;

      if (exists $astral->{gdom}->{$osid1->[$j]}) {
         $osid1->[$j] = $astral->{gdom}->{$osid1->[$j]} ; }

      if (exists $astral->{gdom}->{$osid2->[$j]}) {
         $osid2->[$j] = $astral->{gdom}->{$osid2->[$j]} ; }

      my $fam12;
      if ($class2->[$j] lt $class1->[$j]) {
         $revfl->[$j] = 1 ;
         $fam12 = $class2->[$j]."\t".$class1->[$j];
      } else {
         $revfl->[$j] = 0 ;
         $fam12 = $class1->[$j]."\t".$class2->[$j];
      }
      $interacting_classes->{fam}->{$fam1}->{$sid1->[$j]}->{$sid2->[$j]} = $j;
      $interacting_classes->{fam}->{$fam2}->{$sid2->[$j]}->{$sid1->[$j]} = $j;
      $interacting_classes->{sf}->{$sf1}->{$sid1->[$j]}->{$sid2->[$j]} = $j;
      $interacting_classes->{sf}->{$sf2}->{$sid2->[$j]}->{$sid1->[$j]} = $j;

      $class2sidpairs->{fam}->{$fam1}->{$j}++ ;
      $class2sidpairs->{fam}->{$fam2}->{$j}++ ;
      $class2sidpairs->{sf}->{$sf1}->{$j}++ ;
      $class2sidpairs->{sf}->{$sf2}->{$j}++ ;
   }


   my $chain_2_pdbchain ;
   my ($chain_2_start, $chain_2_end) ;
   my $pdbchains ;
   {
      my ($tbdp_id, $real_chain_id, $pdb_chain_id, $resno_start, $resno_end) ;
      {
         my ($t_bdp, $t_real_cid, $t_pdb_cid, $t_start, $t_end, $t_type) =
         pibase::rawselect_tod(
            "SELECT bdp_id, real_chain_id, pdb_chain_id, ".
            "start_resno, end_resno, chain_type FROM bdp_chains") ;
    
         foreach my $j ( 0 .. $#{$t_bdp}) {
            if ($bdp_2_raw->{$t_bdp->[$j]} != 1) {next;}
            if (!defined $t_pdb_cid->[$j] || $t_pdb_cid->[$j] eq '') {next;}
            if ($t_type->[$j] ne 'p') {next;}
            push @{$tbdp_id}, $t_bdp->[$j] ;
            push @{$real_chain_id}, $t_real_cid->[$j] ;
            push @{$pdb_chain_id}, $t_pdb_cid->[$j] ;
            push @{$resno_start}, $t_start->[$j] ;
            push @{$resno_end}, $t_end->[$j] ;
         }
      }

#      print STDERR ($#{$tbdp_id} + 1)." chains loaded entries\n " ;

      foreach my $j ( 0 .. $#{$tbdp_id}) {
         $pdbchains->{$bdp2pdb->{$tbdp_id->[$j]}}->{$real_chain_id->[$j]}++;

         $chain_2_pdbchain->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
            $pdb_chain_id->[$j] ;
         $chain_2_start->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
            $resno_start->[$j] ;
         $chain_2_end->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
            $resno_end->[$j] ;
      }
   }

   my $pb = {
      bdp2contactsfn => $bdp2contactsfn,
      bdp2subsetsresfn => $bdp2subsetsresfn,
      bdp2pdb => $bdp2pdb,
      pdb2bdp => $pdb2bdp,
      bdp_id => $bdp_id,
      sid2pdb => $sid2pdb,
      sid2class => $sid2class,
      osid2class => $osid2class,
      sid2osid => $sid2osid,
      sid1 => $sid1,
      sid2 => $sid2,
      class1 => $class1,
      class2 => $class2,
      osid1 => $osid1,
      osid2 => $osid2,
      revfl => $revfl,
      interacting_classes => $interacting_classes,
      class2sidpairs => $class2sidpairs,
      pdbchains => $pdbchains,
      chain_2_pdbchain => $chain_2_pdbchain,
      chain_2_start => $chain_2_start,
      chain_2_end => $chain_2_end,
      bdp2sid12 => $bdp2sid12,
      sid12chain => $sid12chain
   } ;

   return $pb ;
}

sub _pilig_pibase_preload {

   my $in = shift ;
   my $astral = $in->{astral} ;
   my $pilig_specs = set_pilig_specs() ;

   my $dbh = $in->{dbh} ;

   my $bdp2contactsfn = pibase::mysql_hashload($dbh,
      "SELECT bdp_id, source_file FROM interface_contacts_tables") ;

   my $bdp2subsetsresfn = pibase::mysql_hashload($dbh,
      "SELECT bdp_id, source_file FROM subsets_residues_tables") ;

   my $sid2class ;
   my $osid2class ;
   $sid2class->{fam} = pibase::mysql_hashload($dbh,
      "SELECT subset_id, class FROM subsets WHERE bdp_id IS NOT NULL
      and subset_source_id = 1") ;

   my $sid2osid ;
   foreach my $sid (keys %{$sid2class->{fam}}) {
      my $fam = $sid2class->{fam}->{$sid} ;
      my $osid = substr($sid, -7, 7) ;
      if (exists $astral->{gdom}->{$sid}) {
         $osid = $astral->{gdom}->{$sid} ; }
      $sid2osid->{$sid} = $osid ;
      my ($sf) = ($fam =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
      $sid2class->{sf}->{$sid} = $sf ;
      $osid2class->{sf}->{$osid} = $sf ;
      $osid2class->{fam}->{$osid} = $fam ;
   }

   my ($bdp_id, $sid1, $sid2, $class1, $class2) = 
      pibase::mysql_fetchcols($dbh,
      "SELECT a.bdp_id, subset_id_1, subset_id_2, class_1, class_2 ".
      "FROM intersubset_contacts as a, bdp_files as b ".
      "WHERE subset_id_1 LIKE \"\%SCOP\%\" AND a.bdp_id = b.bdp_id ".
      "AND num_contacts >= ".$pilig_specs->{PARAM_MIN_NUMCONTACTS}." ".
      "AND b.raw_pdb = 1") ;

   my $bdp2sid12 ;
   my ($interacting_classes, $class2sidpairs) ;
   my ($osid1, $osid2, $revfl) ; #revfl if fam2 <fam1 and flipped to get fam21 isntead of fam12, let it be known so later it can be dealt with

   foreach my $j ( 0 .. $#{$bdp_id}) {
# ONLY ONE OF THEM HAS TO BE OF PROPER CLASS
#      if ($class1->[$j] !~ /^[a-g]/ || $class2->[$j] !~ /^[a-g]/) {next;}
      if ($class1->[$j] !~ /^[a-g]/ && $class2->[$j] !~ /^[a-g]/) {next;}
      $bdp2sid12->{$bdp_id->[$j]}->{$sid1->[$j]."\t".$sid2->[$j]}++ ;

      my $fam1 = $class1->[$j];
      my $fam2 = $class2->[$j] ;

      my ($sf1) = ($fam1 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
      my ($sf2) = ($fam2 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;

      $osid1->[$j] = substr($sid1->[$j], -7, 7) ;
      $osid2->[$j] = substr($sid2->[$j], -7, 7) ;

      if (exists $astral->{gdom}->{$osid1->[$j]}) {
         $osid1->[$j] = $astral->{gdom}->{$osid1->[$j]} ; }

      if (exists $astral->{gdom}->{$osid2->[$j]}) {
         $osid2->[$j] = $astral->{gdom}->{$osid2->[$j]} ; }

      my $fam12;
      if ($class2->[$j] lt $class1->[$j]) {
         $revfl->[$j] = 1 ;
         $fam12 = $class2->[$j]."\t".$class1->[$j];
      } else {
         $revfl->[$j] = 0 ;
         $fam12 = $class1->[$j]."\t".$class2->[$j];
      }
      $interacting_classes->{fam}->{$fam1}->{$sid1->[$j]}->{$sid2->[$j]} = $j;
      $interacting_classes->{fam}->{$fam2}->{$sid2->[$j]}->{$sid1->[$j]} = $j;
      $interacting_classes->{sf}->{$sf1}->{$sid1->[$j]}->{$sid2->[$j]} = $j;
      $interacting_classes->{sf}->{$sf2}->{$sid2->[$j]}->{$sid1->[$j]} = $j;

      $class2sidpairs->{fam}->{$fam1}->{$j}++ ;
      $class2sidpairs->{fam}->{$fam2}->{$j}++ ;
      $class2sidpairs->{sf}->{$sf1}->{$j}++ ;
      $class2sidpairs->{sf}->{$sf2}->{$j}++ ;
   }

   my $pdb2bdp = pibase::mysql_hashload($dbh,
      "SELECT pdb_id, bdp_id FROM bdp_files WHERE raw_pdb =1 ") ;

   my $bdp2pdb = pibase::mysql_hashload($dbh,
      "SELECT bdp_id, pdb_id FROM bdp_files") ;

   my $chain_2_pdbchain ;
   my ($chain_2_start, $chain_2_end) ;
   my ($chain_2_startser, $chain_2_endser) ;
   my $pdbchains ;
   {
      my ($tbdp_id, $real_chain_id, $pdb_chain_id, $resno_start,
         $resno_start_serial, $resno_end, $resno_end_serial) =
      pibase::mysql_fetchcols($dbh,
      "SELECT a.bdp_id, real_chain_id, pdb_chain_id, ".
      "start_resno, start_resno_int, ".
      "end_resno, end_resno_int FROM bdp_chains as a, bdp_files as b ".
      "WHERE pdb_chain_id IS NOT NULL AND chain_type = \"p\" AND ".
      "a.bdp_id = b.bdp_id AND raw_pdb = 1");

#      print STDERR ($#{$tbdp_id} + 1)." chains loaded entries\n " ;

      foreach my $j ( 0 .. $#{$tbdp_id}) {
         $pdbchains->{$bdp2pdb->{$tbdp_id->[$j]}}->{$real_chain_id->[$j]}++;

         $chain_2_pdbchain->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
            $pdb_chain_id->[$j] ;
         $chain_2_start->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
            $resno_start->[$j] ;
         $chain_2_end->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
            $resno_end->[$j] ;
         $chain_2_startser->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
            $resno_start_serial->[$j] ;
         $chain_2_endser->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
            $resno_end_serial->[$j] ;
      }
   }

   my $pb = {
      bdp2contactsfn => $bdp2contactsfn,
      bdp2subsetsresfn => $bdp2subsetsresfn,
      bdp2pdb => $bdp2pdb,
      pdb2bdp => $pdb2bdp,
      bdp_id => $bdp_id,
      sid2class => $sid2class,
      sid2osid => $sid2osid,
      osid2class => $osid2class,
      sid1 => $sid1,
      sid2 => $sid2,
      class1 => $class1,
      class2 => $class2,
      osid1 => $osid1,
      osid2 => $osid2,
      revfl => $revfl,
      interacting_classes => $interacting_classes,
      class2sidpairs => $class2sidpairs,
      pdbchains => $pdbchains,
      chain_2_pdbchain => $chain_2_pdbchain,
      chain_2_start => $chain_2_start,
      chain_2_end => $chain_2_end,
      chain_2_startser => $chain_2_startser,
      chain_2_endser => $chain_2_endser,
      bdp2sid12 => $bdp2sid12
   } ;

   return $pb ;
}


sub _pilig_get_astral_classlist {

   my $in = shift;

   my $classlist ;
   foreach my $type (qw/fam sf/) {
      my @files = glob($in->{fn}->{aln}->{$type}."/*.aln.fa") ;
      for (@files) {
         my $a = $_ ;
         $a =~ s/\.aln\.fa$// ; $a =~ s/.*\/// ;
         $classlist->{$type}->{$a}++ ;
       }
   }

   return $classlist ;
  
}




sub _pilig_get_subsres {
   
   my $in = shift ;
   my $fn = $in->{fn} ;
      
   my ($t_subset_id, $t_chain_id, $t_resno_serial, $t_resno) =
      pibase::rawselect_metatod($fn,
      "SELECT subset_id, chain_id, resno_serial, resno FROM $fn") ;
      
   my $res2sub ;
   foreach my $j ( 0 .. $#{$t_subset_id}) {
      if ($t_subset_id->[$j] =~ /SCOP/) {
         my $sig = $t_resno->[$j]."\n".$t_chain_id->[$j] ;
         $res2sub->{$sig} = $t_subset_id->[$j] ; }
   }

   return $res2sub ;
   
}


sub _pilig_get_ligbase_domains {

   my $in = shift ;
   my $pilig_specs = set_pilig_specs() ;

   my $pb = $in->{pb} ;
   my $lb = $in->{lb} ;

   my ($ligres2sid, $sid2ligres, $class2ligsid) ;

   my $count = 0 ;
   foreach my $pdb_id (sort keys %{$lb->{pdb2res2ligid}}) {
#      print STDERR "PROGRESS: get_ligbase_domains() now on $pdb_id\n" ;
      if (!exists $pb->{pdb2bdp}->{$pdb_id}) {
         print STDERR "ERROR get_ligbase_domains() ".__LINE__.
            ": $pdb_id no bdp entry\n"; next ; }
      my $bdp_id = $pb->{pdb2bdp}->{$pdb_id} ;

      if (!exists $pb->{bdp2subsetsresfn}->{$bdp_id}) {
         print STDERR "ERROR get_ligbase_domains() ".__LINE__.
            ": $pdb_id ($bdp_id) no subsets entry\n";
	 next ; }

      if (!-s $pb->{bdp2subsetsresfn}->{$bdp_id}) {
         print STDERR "ERROR get_ligbase_domains() ".__LINE__.": $pdb_id ".
            "empty subets_res file $pb->{bdp2subsetsresfn}->{$bdp_id}\n" ;
         next ;
      }

      my $curscopdef = _pilig_get_subsres({
         fn => $pb->{bdp2subsetsresfn}->{$bdp_id} }) ;

      foreach my $orig_ressig (keys %{$lb->{pdb2res2ligid}->{$pdb_id}}) {
         my $ressig = $orig_ressig ;
         if (!exists $curscopdef->{$ressig} && 
             $ressig =~ / $/) {$ressig =~ s/ $/A/;}

         if (exists $curscopdef->{$ressig}) {
            my $sid = $curscopdef->{$ressig} ;
	    $ligres2sid->{$pdb_id}->{$ressig}->{$sid}++ ;
	    $sid2ligres->{$sid}->{$ressig}++ ;
            $class2ligsid->{fam}->{$pb->{sid2class}->{fam}->{$sid}}->{$sid}++ ;
            $class2ligsid->{sf}->{$pb->{sid2class}->{sf}->{$sid}}->{$sid}++ ;
#            print STDERR "FOUND: $pdb_id $ressig in domain $sid is ligand bound\n" if DEBUG;
	 } else {
            my $tres = $ressig ; $tres =~ s/\n/:/ ;
            print STDERR "WARNING get_ligbase_domains() ".__LINE__.
               ": $pdb_id no domain definition for ligbase residue $tres\n" ;
         }
      }
      $count++ ;
      if ($pilig_specs->{DEBUG} && $count == 100) {last;}
   }

   return {
      ligres2sid => $ligres2sid,
      sid2ligres => $sid2ligres,
      class2ligsid => $class2ligsid
   } ;
}

sub _pilig_load_liginfo {

   my $in = shift ;
   my $pilig_specs = set_pilig_specs() ;

   my $liginfo ;
   open(LIGINFOF, $pilig_specs->{outfiles}->{liginfo}) ;
   while (my $line = <LIGINFOF>) {
      if ($line =~ /^\#/) {next;}
      chomp $line;

      my @t = split(/\t/, $line) ;
      my $lig = shift @t ;
      $liginfo->{ligs}->{$lig}++ ;

      ($liginfo->{name}->{$lig},
       $liginfo->{formula}->{$lig},
       $liginfo->{mw}->{$lig},
       $liginfo->{numatoms}->{$lig},
       $liginfo->{numatoms_nonh}->{$lig} ) = @t ;

      if ($liginfo->{mw}->{$lig} eq '') {
         delete $liginfo->{mw}->{$lig};
      } elsif ($liginfo->{mw}->{$lig} >= $pilig_specs->{PARAM_MIN_LIGMW} &&
               $liginfo->{mw}->{$lig} <= $pilig_specs->{PARAM_MAX_LIGMW}) {
         $liginfo->{mwinrange}->{$lig}++ ;
      }
   }
   close(LIGINFOF) ;

   return $liginfo ;

}


sub extract_msdchem_liginfo {

   require Chemistry::MolecularMass ;

   my $mm = new Chemistry::MolecularMass ;
   my $pilig_specs = set_pilig_specs() ;

   my $msdchemdir = $pilig_specs->{msdchem_xmldir} ;
   my @xml = <$msdchemdir/*.xml> ;
   my $data ;
   my $ligands ;
   foreach my $file (@xml) {
      open(CHEMF, $file) ;
      my $curcode = '';
      while (my $line = <CHEMF>) {
         chomp $line;
         $line =~ s/^\s+// ;
         if ($line =~ /<code3Letter>/) {
            $curcode = $line ;
            $curcode =~ s/<code3Letter>// ;
            $curcode =~ s/<\/code3Letter>// ;
            $ligands->{$curcode}++ ;

         } elsif ($line =~ /<formula>/ &&
            (!exists $data->{formula}->{$curcode} ||
             !exists $data->{mw}->{$curcode})) {

            my $t = $line;
            $t =~ s/<formula>// ;
            $t =~ s/<\/formula>// ;
            $t = lc($t) ;
            my @bits = split(' ', $t) ;
            foreach my $j ( 0 .. $#bits) {
               $bits[$j] =~ s/^([a-z])/\U$1\E/ ; }
            $data->{formula}->{$curcode} = join(' ', @bits) ;
            $t = $data->{formula}->{$curcode} ; $t =~ s/ //g ;
            $data->{mw}->{$curcode} = $mm->calc_mass($t) ;

         } elsif ($line =~ /<numAtomsAll>/ &&
             !exists $data->{numatoms}->{$curcode}) {

            my $t = $line;
            $t =~ s/<numAtomsAll>// ;
            $t =~ s/<\/numAtomsAll>// ;
            $data->{numatoms}->{$curcode} = $t ;

         } elsif ($line =~ /<numAtomsNonH>/ &&
             !exists $data->{numatoms_nonh}->{$curcode}) {

            my $t = $line;
            $t =~ s/<numAtomsNonH>// ;
            $t =~ s/<\/numAtomsNonH>// ;
            $data->{numatoms_nonh}->{$curcode}  = $t ;

         } elsif ($line =~ /<name>/ &&
            !exists $data->{name}->{$curcode}) {

            my $t = $line;
            $t =~ s/<name>// ;
            $t =~ s/<\/name>// ;
            $data->{name}->{$curcode}  = $t ;
         }
      }
      close(CHEMF) ;
   }

   open(LIGINFOF, ">".$pilig_specs->{outfiles}->{liginfo}) ;
   print LIGINFOF '#'.
      join("\t", qw/ligand name formula mw numatoms numatoms_nonh/)."\n";
   foreach my $lig (sort keys %{$ligands}) {
      my @outvals = ($lig,
         $data->{name}->{$lig},
         $data->{formula}->{$lig},
         $data->{mw}->{$lig},
         $data->{numatoms}->{$lig},
         $data->{numatoms_nonh}->{$lig},
      ) ;
      print LIGINFOF join("\t", @outvals)."\n" ;
   }
   close(LIGINFOF) ;

}


sub _pilig_get_sasa { # purpose: run and parse modeller sasa (psa)

   my $in = shift ;

   my ($bdp_file, $modeller_bin) ;
   my $surftyp = '1' ; #default contact

   if (ref($in) ne '') {
      if (exists $in->{surftyp}) {
         $surftyp = $in->{surftyp} ; }

      $bdp_file = $in->{pdb_fn} ;
      $modeller_bin = $in->{modeller_bin} ;
   } else {
      $bdp_file = $in ;
      $modeller_bin = shift ;
   }

#Specify MODELLER temporary TOP file, and output files.

   my $filesmade ;

   my ($sasa_top_fh, $sasa_top_fn) = 
      tempfile("temp_sasa.XXXXXX", SUFFIX=>".top") ;
   $filesmade->{$sasa_top_fn}++ ;

   my $tempbase = $sasa_top_fn ; $tempbase =~ s/\.top$// ;
   my $sasa_log_fn = $tempbase.".log" ;
   my $sasa_atm_fn = $tempbase.".sol" ;
   my $sasa_res_fn = $tempbase.".psa" ;

   $filesmade->{$sasa_log_fn}++ ;
   $filesmade->{$sasa_atm_fn}++ ;
   $filesmade->{$sasa_res_fn}++ ;

#Write the MODELLER TOP file to calcaulte solvent accessibilty.

   my $filename = $bdp_file ; $filename =~ s/^.*\///g ;
   my $filebase = $bdp_file ; $filebase =~ s/\/[^\/]+$//g ;

   if ($filename ne $bdp_file) {
      print $sasa_top_fh "SET ATOM_FILES_DIRECTORY = \'$filebase\'\n" ; }

   print $sasa_top_fh 'READ_TOPOLOGY FILE = \'$(LIB)/top_heav.lib\''."\n" ;
   print $sasa_top_fh 'SET HETATM_IO = off, WATER_IO=off'."\n" ;
   print $sasa_top_fh 'READ_MODEL FILE = \''.$filename.'\''."\n" ;
   print $sasa_top_fh 'SET SURFTYP = '."$surftyp\n" ;
   print $sasa_top_fh 'SET RADII_FACTOR = 1.0'."\n" ;
   print $sasa_top_fh 'WRITE_DATA FILE = \''.$tempbase.'\', OUTPUT = \'PSA\''."\n" ;
   close($sasa_top_fh) ;

#Run the MODELLER TOP file.

   my $tcom = "$modeller_bin $sasa_top_fn >/dev/null 2>&1" ;
   system($tcom) ;

#Initialize solvent accessibility storage arrays

   my $results ;
   my $resno_rev ;
   my %fields = (
      'resno' => [7, 4] ,
      'resna' => [14, 3] ,
      'chain' => [18, 1],
      'all_sum' => [19, 7],
      'all_perc' => [27, 5],
      'nonp_sum' => [33, 7],
      'nonp_perc' => [41, 5],
      'p_sum' => [47, 7],
      'p_perc' => [55, 5],
      'sc_sum' => [61, 7],
      'sc_perc' => [69, 5],
      'mc_sum' => [75, 7],
      'mc_perc' => [83, 5]
   ) ;

   foreach my $key (keys %fields) {
      $results->{$key} = [] ; }

#If the MODELLER PSA file does not exist, display to STDERR and return.

   if (!-s $sasa_res_fn) {
      $results->{resno}->[0] = 'ERROR' ;
      $resno_rev->{ERROR} = 'ERROR' ;
      my $err = ['DANGER WILL ROBINSON'] ;
      foreach my $tfn (keys %{$filesmade}) {if (-e $tfn) {unlink $tfn;}}
      return ($results, $resno_rev, undef, $err);
   }

#Open the resulting MODELLER PSA file.


#Read in a line of the PSA file.

   my $sasaatoms ;
   $sasaatoms->{all} = 0 ;
   $sasaatoms->{p} = 0 ;
   $sasaatoms->{nonp} = 0 ;

   open(SOL, $sasa_atm_fn) ;
   while (my $line = <SOL>) {
      chomp $line;
      if ($line !~ /^ATOM/) {next;}

      my $atomtype = substr($line,12,1) ;
      my $atomacc = substr($line,64,8) ;
      $atomacc =~ s/ //g ;
      $sasaatoms->{all} += $atomacc ;
      if ($atomtype =~ /[NO]/) {
         $sasaatoms->{p} += $atomacc ;
      } else {
         $sasaatoms->{nonp} += $atomacc ;
      }
   }
   close(SOL) ;


   my $sasa ;
   $sasa->{all} = 0 ;
   $sasa->{mc} = 0 ;
   $sasa->{sc} = 0 ;
   $sasa->{nonp} = 0 ;
   $sasa->{p} = 0 ;
   open(PSA, $sasa_res_fn) ;
   while (my $line = <PSA>) {

      chomp $line;

#If the line contains accessibility information,

      if ($line =~ /^ACCESS/) {

#Extract the residue information and solvent accessibility parameters from the line

         foreach my $key (keys %fields) {

            my $t = substr($line, $fields{$key}->[0], $fields{$key}->[1]);
            unless ($key eq 'chain') {
               $t =~ s/ //g ; }

	    if ($key eq 'all_sum') {
	       $sasa->{all} += $t;
	    } elsif ($key eq 'mc_sum') {
	       $sasa->{mc} += $t;
	    } elsif ($key eq 'sc_sum') {
	       $sasa->{sc} += $t;
	    } elsif ($key eq 'nonp_sum') {
	       $sasa->{nonp} += $t;
	    } elsif ($key eq 'p_sum') {
	       $sasa->{p} += $t;
	    }

            push @{$results->{$key}}, $t ;
         }

#Create an entry in the reverse lookup hash pointing from the actual PDB residue number and chain id to the index of the arrays where the solvent accessibility values are stored.

         my $cur_recno = $#{$results->{resno}} ;
         my $res_sig = $results->{resno}->[$cur_recno]."_".
                       $results->{chain}->[$cur_recno] ;

         $resno_rev->{$res_sig} = $cur_recno ;

      }

   }

#Close the PSA file.

   close PSA ;

#If no accessibility lines have been read in, set error flag

   my $error_fl ;

   if ($#{$results->{resno}} < 0 ) {
      push @{$error_fl}, "no sasa entries calculated" ; }

#Remove the MODELLER TOP file and output files.

   unlink($sasa_top_fn, $sasa_log_fn, $sasa_atm_fn, $sasa_res_fn) ;

#Return pointers to teh arrays holding the values and the residue lookup hash.
#
#Note: accessibility values for ALL Cys read are retunred; However, the calling function will only ask for Cys involved in disulfide bridges.

   return ( $results, $resno_rev, $sasa, $error_fl, $sasaatoms ) ;

}


sub _pilig_readin_ligassignments {

   my $in = shift ;

   my $fn = $in->{fn} ;

   my $class2alnlength = $in->{class2alnlength};
   my $standardres= $in->{standardres};
   my $liginfo = $in->{liginfo};
   my $ligbits ;
   print STDERR "NOW reading LIG assignment\n" ;
   open(LIGFH, $fn->{pilig}->{park_assign_lig_fn}) ;
   while (my $line = <LIGFH>) {
      if ($line =~ /^#/ ||
          $line =~ /^Warning: no access to/ ||
          $line =~ /^Thus no job control/ ) {next;}
      chomp $line;

      my ( $pdb, $sid, $osid, $classtype, $class, 
              $btype, $alnposstring, $alnlength,
              $ligcod, $ligid ) = split(/\t/, $line) ;
      $ligcod =~ s/ //g ;
      my $ligsig = $pdb."\t".$ligcod."\t".$ligid ;

      if (exists $standardres->{$ligcod}) {next;}
      if ($ligcod eq 'DUM' || $ligcod eq 'UNX' ||
             $ligcod eq 'UNK' || $ligcod eq 'UNL') {next;}

      if (!exists $liginfo->{mw}->{$ligcod}) {
            print STDERR "WARNING: MW not found for ligand $ligcod\n" ;}

      $class2alnlength->{$classtype}->{$class} = $alnlength ;

      $ligbits->{$classtype}->{$sid}->{$ligsig} = Bit::Vector->new($alnlength) ;

      my @alnpos = split(/\,/, $alnposstring) ;
      foreach my $alnpos (@alnpos) {
         $ligbits->{$classtype}->{$sid}->{$ligsig}->Bit_On($alnpos); }


      if (!exists $ligbits->{$classtype}->{$sid}->{cumulative}) {
         $ligbits->{$classtype}->{$sid}->{cumulative} =
            $ligbits->{$classtype}->{$sid}->{$ligsig}->Clone();
      } else {
         $ligbits->{$classtype}->{$sid}->{cumulative}->Union(
            $ligbits->{$classtype}->{$sid}->{cumulative},
            $ligbits->{$classtype}->{$sid}->{$ligsig}) ;
      }
#      if (exists $liginfo->{mwinrange}->{$ligcod})
   }
   close(LIGFH) ;

   return $ligbits ;
}


sub _pilig_readin_expassignments {

   my $in = shift ;
   my $fn = $in->{fn};
   my $class2alnlength = $in->{class2alnlength};

   print STDERR "NOW reading EXP assignment\n" ;
   open(EXPFH, $fn->{pilig}->{park_assign_exp_fn}) ;
   my $expbits ;
   while (my $line = <EXPFH>) {
      if ($line =~ /^#/ ||
            $line =~ /^Warning: no access to/ ||
            $line =~ /^Thus no job control/ ) {next;}
      chomp $line;

      my ( $pdb, $sid, $osid, $classtype,
           $class, $btype, $alnposstring, $alnlength) = split(/\t/, $line) ;

      $class2alnlength->{$classtype}->{$class} = $alnlength ;

      if (!exists $expbits->{$classtype}->{$class}) {
         $expbits->{$classtype}->{$class} =
            Bit::Vector->new($alnlength) ;
      }

      my @alnpos = split(/\,/, $alnposstring) ;
      foreach my $alnpos (@alnpos)  {
         if ($alnpos eq "undef") {next;}
         $expbits->{$classtype}->{$class}->Bit_On($alnpos);
      }
   }
   close(EXPFH) ;

   return $expbits ;

}


sub _pilig_readin_piassignments {

   my $in = shift;
   my $fn = $in->{fn} ;
   my $class2alnlength = $in->{class2alnlength};
   my $pb = $in->{pb} ;
   my $pibits ;
   my $interfaces ;

   print STDERR "NOW reading PI assignment\n" ;
   open(PIFH, $fn->{pilig}->{park_assign_pi_fn}) ;
   while (my $line = <PIFH>) {
      if ($line =~ /^#/ ||
         $line =~ /^Warning: no access to/ ||
         $line =~ /^Thus no job control/ ) {next;}
      chomp $line;

      my ($pdb, $sid, $osid, $classtype, $class,
          $obtype, $alnposstring, $alnlength,
          $sid1, $sid2, $fam1, $fam2, $chains) = split(/\t/, $line) ;

      my $sid12 = $sid1."\t".$sid2 ;
      my $side = 1; if ($sid eq $sid2) {$side = 2;}

      $interfaces->{$sid12}->{$side}->{sid} = $sid ;
      $interfaces->{$sid12}->{chains} = $chains ;
      $interfaces->{$sid12}->{pdb} = $pdb;

      my @alnpos = split(/\,/, $alnposstring) ;
      {
            my @t = ();
            foreach my $p ( @alnpos) {
               if ($p ne 'undef') {push @t, $p;} }
            @alnpos = @t ;
      }
      if ($#alnpos < 0) {next;}

      $class2alnlength->{$classtype}->{$class} = $alnlength ;

      my @btypes = ($obtype) ;
      if ($chains eq 'same') {
         push @btypes, "Pintra" ;
      } else {
         push @btypes, "Pinter" ;
      }

      $interfaces->{$sid12}->{$side}->{pibits}->{$classtype} =
         Bit::Vector->new($alnlength) ;

      foreach my $alnpos (@alnpos)  {
         if ($alnpos eq 'undef') {next;}
         $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Bit_On($alnpos); }

      foreach my $btype (@btypes) {
         if (!exists $pibits->{$classtype}->{$sid}->{$btype}) {
            $pibits->{$classtype}->{$sid}->{$btype} =
               $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Clone();
         } else {
            $pibits->{$classtype}->{$sid}->{$btype}->Union(
               $pibits->{$classtype}->{$sid}->{$btype},
               $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}) ;
         }
      }
   }
   close(PIFH) ;

   return {
      pibits => $pibits,
      interfaces => $interfaces
   } ;
}


#      if ($classtype eq 'fam') {
#         foreach my $seqid (qw/90 95 100/) {
#            push @curclasstypes, 'seqcl'.$seqid ;
#            my $tosid = $osid;
#            if (substr($tosid,5,1) eq '.') { substr($tosid,0,1) = 'g' ; }
#            my $thiscl = $tosid ;
#            if (defined $astral->{seqcl}->{$seqid}->{$tosid}) {
#                  $thiscl = $astral->{seqcl}->{$seqid}->{$tosid} ;
#            } else {
#                  $astral->{seqcl}->{$seqid}->{$tosid} = $thiscl;
#                  $astral->{seqcl2cont}->{$seqid}->{$thiscl} = [$tosid];
#            }
#            push @curclasses, $thiscl ;
#         }
#      }


sub _pilig_BK061111_run_collate_perinstance {

   require Bit::Vector ;


   my $standardres = {
   ALA => 'A' ,
   ARG => 'R' ,
   ASN => 'N' ,
   ASP => 'D' ,
   CYS => 'C' ,
   GLN => 'Q' ,
   GLU => 'E' ,
   GLY => 'G' ,
   HIS => 'H' ,
   HSD => 'H' ,
   HSE => 'H' ,
   ILE => 'I' ,
   LEU => 'L' ,
   LYS => 'K' ,
   MET => 'M' ,
   PHE => 'F' ,
   PRO => 'P' ,
   SER => 'S' ,
   THR => 'T' ,
   TRP => 'W' ,
   TYR => 'Y' ,
   VAL => 'V',
   UNK => 'X',
   '  C' => 'c',
   '  G' => 'g',
   '  A' => 'a',
   '  T' => 't',
   '  U' => 'u',
   '  I' => 'i',
   'C' => 'c',
   'G' => 'g',
   'A' => 'a',
   'T' => 't',
   'U' => 'u',
   'I' => 'i',
   '+C' => 'c',
   '+G' => 'g',
   '+A' => 'a',
   '+T' => 't',
   '+U' => 'u',
   '+I' => 'i'
   } ;

# Per-PPI
# 0. load all of ligands first - lig info
# - iterate over assignments, and have individual
#   bit vectors set for the ligand binding sites in their
#   respective domain families
#
# per-PPI: iterate over each domain--domain interface
#  1. list any ligands actually bound to this interface
#  2. iterate over all ligands in that family/superfamily, find
#     - quantify per ligand coverage
#     - quantify cumulative ligand coverage
#
# per-ligand:


   my $fn = set_locations() ;
   my $liginfo = load_liginfo({fn => $fn}) ;


   my $class2lig2bit ;
   my $lig2class ;

   my $class2bits ;
   my $class2alnlength;
   my $class2ligbits ;
      print STDERR "NOW reading LIG assignment\n" ;
      open(LIGFH, $fn->{pilig}->{park_assign_lig_fn}) ;
      while (my $line = <LIGFH>) {
         if ($line =~ /^#/ ||
            $line =~ /^Warning: no access to/ ||
            $line =~ /^Thus no job control/ ) {next;}
         chomp $line;

         my ( $pdb, $sid, $osid, $classtype, $class, 
              $btype, $alnposstring, $alnlength,
              $ligcod, $ligid ) = split(/\t/, $line) ;

         if (exists $standardres->{$ligcod}) {next;}
         if ($ligcod eq 'DUM' || $ligcod eq 'UNX' ||
             $ligcod eq 'UNK' || $ligcod eq 'UNL') {next;}

#         if (!exists $liginfo->{mwinrange}->{$ligcod}) {next;}

         $class2alnlength->{$classtype}->{$class} = $alnlength ;
         my $lig = $pdb."\t".$ligcod."\t".$ligid ;
         $class2ligbits->{all}->{$classtype}->{$class}->{$lig} = Bit::Vector->new($alnlength) ;

         $lig2class->{$ligcod}->{$classtype}->{$class}++ ;

         $ligcod =~ s/ //g ;
         if (!exists $liginfo->{mw}->{$ligcod}) {
            print STDERR "WARNING: $fn MW not found for ligand $ligcod\n" ;}

         if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
            $class2bits->{$classtype}->{$class}->{$btype} =
               Bit::Vector->new($alnlength) ;
         }

         my @alnpos = split(/\,/, $alnposstring) ;
         foreach my $alnpos (@alnpos)  {
            $class2ligbits->{all}->{$classtype}->{$class}->{$lig}->Bit_On($alnpos) ;}
      }
      close(LIGFH) ;


      print STDERR "NOW reading EXP assignment\n" ;
      open(EXPFH, $fn->{pilig}->{park_assign_exp_fn}) ;
      while (my $line = <EXPFH>) {
         if ($line =~ /^#/ ||
            $line =~ /^Warning: no access to/ ||
            $line =~ /^Thus no job control/ ) {next;}
         chomp $line;

         my ( $pdb, $sid, $osid, $classtype,
              $class, $btype, $alnposstring, $alnlength) = split(/\t/, $line) ;

         $class2alnlength->{$classtype}->{$class} = $alnlength ;

         if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
            $class2bits->{$classtype}->{$class}->{$btype} =
               Bit::Vector->new($alnlength) ;
         }

         my @alnpos = split(/\,/, $alnposstring) ;
         foreach my $alnpos (@alnpos)  {
            if ($alnpos eq "undef") {next;}
            $class2bits->{$classtype}->{$class}->{$btype}->Bit_On($alnpos); }
      }
      close(EXPFH) ;


 
   my @headers_sum=  ("PINT_LSUM", "PDB", "SID1", "SID2", "CHAINS",
                       "CLASSTYPE1", "CLASS1", "CLASSTYPE2", "CLASS2",
                       "numres_p_1", "cumulative_numres_l_and_p_1",
                       "max_liginrange_1",
                       "max_l_and_p_1", "lig_max_l_and_p_1",
                       "max_obi_1", "lig_max_obi_1",
                       "max_opi_1", "lig_max_opi_1",
                       "max_olig_1", "lig_max_olig_1",
                       "numres_p_2", "cumulative_numres_l_and_p_2",
                       "max_liginrange_2",
                       "max_l_and_p_2", "lig_max_l_and_p_2",
                       "max_obi_2", "lig_max_obi_2",
                       "max_opi_2", "lig_max_opi_2",
                       "max_olig_2", "lig_max_olig_2",
                       ) ;
   print '#'.join("\t",@headers_sum)."\n" ;

   print STDERR "NOW reading PI assignment\n" ;
   open(PIFH, $fn->{pilig}->{park_assign_pi_fn}) ;
   my $interfaces ;
   while (my $line = <PIFH>) {
      if ($line =~ /^#/ ||
            $line =~ /^Warning: no access to/ ||
            $line =~ /^Thus no job control/ ) {next;}
      chomp $line;

      my ($pdb, $sid, $osid, $classtype, $class, $bstype, $alnposstring,
             $alnlength, $sid1, $sid2, $fam1, $fam2, $chains)=
             split(/\t/, $line) ;
      my $sid12 = $sid1."\t".$sid2 ;

      my $side = 1; if ($sid eq $sid2) {$side = 2;}
      $interfaces->{$sid12}->{$side}->{osid} = $osid ;
      $interfaces->{$sid12}->{$side}->{$classtype}->{class} = $class ;
      $interfaces->{$sid12}->{chains} = $chains ;
      $interfaces->{$sid12}->{pdb} = $pdb ;

#      print STDERR "sid12 = $sid12\nsid =$side\nclasstype=$classtype\nalnposstring=$alnposstring\n" ;
      $interfaces->{$sid12}->{$side}->{$classtype}->{alnposstring} =
            $alnposstring;

      $interfaces->{$sid12}->{$side}->{$classtype}->{alnlength} =
            $alnlength;

      $interfaces->{$sid12}->{$side}->{sid} = $sid2 ;
   }
   close(PIFH) ;

   foreach my $sid12 (keys %{$interfaces}) {
      my $curi = $interfaces->{$sid12} ;
      my $pdb = $curi->{pdb} ;
      my $chains = $curi->{chains} ;

      my ($sid, $alnposstring, $alnlength, $class, $classtype) ;

#fam level analysis only
      $alnlength->{1} = $curi->{1}->{fam}->{alnlength} ;
      $alnposstring->{1} = $curi->{1}->{fam}->{alnposstring} ;

      $alnlength->{2} = $curi->{2}->{fam}->{alnlength} ;
      $alnposstring->{2} = $curi->{2}->{fam}->{alnposstring} ;

      ($sid->{1}, $sid->{2}) = split(/\t/, $sid12) ;

      $classtype->{1} = 'fam' ;
      $classtype->{2} = 'fam' ;
      $class->{1} = $curi->{1}->{fam}->{class} || 'undef';
      $class->{2} = $curi->{2}->{fam}->{class} || 'undef';

      my $surf_pos ;
      $surf_pos->{1} = $class2bits->{$classtype->{1}}->{$class->{1}}->{'E'}->Norm() ;
      $surf_pos->{2} = $class2bits->{$classtype->{2}}->{$class->{2}}->{'E'}->Norm() ;

#         ($pdb, $sid->{1}, $osid->{1}, $classtype->{1}, $class->{1},
#             undef, $alnposstring->{1}, $alnlength->{1},
#             undef, undef, undef, undef, $chains) = split(/\t/, $line1) ;

         my $bs_bits;
         my $curligs ;

         my $alnpos ;
         my $pistats ;
         my $ligstats ; 
         my $skipside ;
         foreach my $s (1, 2) {
            if (!defined $alnposstring->{$s}) {
               $skipside->{$s}++ ;
               print STDERR "WARNING: skipping $sid->{1} -- $sid->{2}".
               " undefined binding site alignment positions for $sid->{$s}\n";
               next;
            }

            $curligs->{$s} =
               $class2ligbits->{all}->{$classtype->{$s}}->{$class->{$s}} ;

            $bs_bits->{$s} = Bit::Vector->new($alnlength->{$s}) ;
            my @t = split(/\,/, $alnposstring->{$s}) ;
            {
               @{$alnpos->{$s}} = () ;
               foreach my $p ( @t) {
                  if ($p ne 'undef') {
                     $bs_bits->{$s}->Bit_On($p) ;
                     push @{$alnpos->{$s}}, $p } }
            }

            if ($#{$alnpos->{$s}} < 0 ) {
               $skipside->{$s}++ ;
               print STDERR "WARNING: skipping $sid->{1} -- $sid->{2}".
               " undefined binding site alignment positions for $sid->{$s}\n"; }
         }

# init vals
         foreach my $s (1, 2) {
            foreach my $type (qw/p cumlig_l_and_p max_liginrange max_l_and_p max_obi max_opi max_olig/) {
               $ligstats->{$s}->{$type} = 0 ; }
            foreach my $type (qw/lig_max_l_and_p lig_max_obi lig_max_opi lig_maxolig/) {
               $ligstats->{$s}->{$type} = 'U'; }
         }

         foreach my $s (1, 2) {
            if (exists $skipside->{$s}) {next;}
            $pistats->{$s}->{p} = $bs_bits->{$s}->Norm() ;
            $pistats->{$s}->{max_liginrange} = 0 ;
            $ligstats->{$s}->{cumlig} = Bit::Vector->new($alnlength->{$s}) ;
            foreach my $t (qw/obi opi olig l_and_p/) {
               $ligstats->{$s}->{"max_$t"} = 0 ;
               $ligstats->{$s}->{"lig_max_$t"} = 'undef' ;
            }

            my $temp_bits_and = Bit::Vector->new($alnlength->{$s}) ;
            my $temp_bits_or = Bit::Vector->new($alnlength->{$s}) ;

            if (!exists $curligs->{$s}) {next;}
            foreach my $lig (keys %{$curligs->{$s}}) {
               $ligstats->{$s}->{cumlig}->Or($ligstats->{$s}->{cumlig},
                                             $curligs->{$s}->{$lig}) ;

            #intersect the ligands bit vector and the bit vector for the bindin
            # site positions of this particular interface;
               $temp_bits_and->And($curligs->{$s}->{$lig}, $bs_bits->{$s}) ;
               $temp_bits_or->Or($curligs->{$s}->{$lig}, $bs_bits->{$s}) ;

               my $curs ;
               $curs->{l_and_p} = $temp_bits_and->Norm() ;
               $curs->{l_or_p} = $temp_bits_or->Norm() ;
               $curs->{l} = $curligs->{$s}->{$lig}->Norm() ;
               $curs->{p} = $bs_bits->{$s}->Norm() ;

               if ($curs->{l_and_p} == 0) {next;}

               my $curoverlap ;
               $curoverlap->{obi} = $curs->{l_and_p} / $curs->{l_or_p} ;
               $curoverlap->{opi} = $curs->{l_and_p} / $curs->{p} ;
               $curoverlap->{olig} = $curs->{l_and_p} / $curs->{l} ;


               if (!exists $ligstats->{$s}->{"max_l_and_p"} ||
                   $curs->{l_and_p} > $ligstats->{$s}->{"max_l_and_p"}) {
                  $ligstats->{$s}->{"max_l_and_p"} = $curs->{l_and_p} ;
                  $ligstats->{$s}->{"lig_max_l_and_p"} = $lig ;
               }

               foreach my $t (qw/obi opi olig/) {
                  if (!exists $ligstats->{$s}->{"max_$t"} ||
                      $curoverlap->{$t} > $ligstats->{$s}->{"max_$t"}) {
                     $ligstats->{$s}->{"max_$t"} = $curoverlap->{$t} ;
                     $ligstats->{$s}->{"lig_max_$t"} = $lig ;
                  }
               }

               my (undef, $ligcod, undef) = split(/\t/, $lig) ;
               my $liginrange = 0 ;
               if (exists $liginfo->{mwinrange}->{$ligcod}) {
                  $ligstats->{$s}->{max_liginrange} = 1 ;
                  $liginrange = 1 ; }

               my @outvals = ("PINT_LINT",
                              $pdb, $sid->{1}, $sid->{2},
                              $classtype->{1}, $class->{1},
                              $classtype->{2}, $class->{2},
                              $s, $sid->{$s}, $lig, $liginrange,
                              $curs->{p}, $curs->{l},
                              $curs->{l_and_p}, $curs->{l_or_p},
                              sprintf("%.3f", $curoverlap->{obi}),
                              sprintf("%.3f", $curoverlap->{opi}),
                              sprintf("%.3f", $curoverlap->{olig})
                              ) ;
               print join("\t", @outvals)."\n";
            }

            $ligstats->{$s}->{"cumlig_l"} = $ligstats->{$s}->{cumlig}->Norm() ;
            $temp_bits_and->And($ligstats->{$s}->{cumlig}, $bs_bits->{$s}) ;
            $ligstats->{$s}->{"cumlig_l_and_p"} = $temp_bits_and->Norm() ;
         }

         my @outvals =  ("PINT_LSUM", $pdb, $sid->{1}, $sid->{2}, $chains,
                              $classtype->{1}, $class->{1},
                              $classtype->{2}, $class->{2}) ;

         foreach my $s ( 1, 2) {
            if (exists $skipside->{$s} ||
                $ligstats->{$s}->{"cumlig_l"} == 0 ) {

               if (exists $skipside->{$s}) {
                  push @outvals, 0, 0, 0 ;
               } else {
                  push @outvals, $pistats->{$s}->{p}, 0, 0 ;
               }

               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
               push @outvals, 0, "U", "U", "U";
               next;
            }

            push @outvals, $pistats->{$s}->{p} ;
            push @outvals, $ligstats->{$s}->{cumlig_l_and_p} ;
            push @outvals, $ligstats->{$s}->{max_liginrange} ;
            push @outvals, $ligstats->{$s}->{"max_l_and_p"} ;
            push @outvals, $ligstats->{$s}->{"lig_max_l_and_p"} ;
            foreach my $t (qw/obi opi olig/)  {
               push @outvals, sprintf("%.3f", $ligstats->{$s}->{"max_$t"}) ;
               push @outvals, $ligstats->{$s}->{"lig_max_$t"} ;
            }
         }
      print join("\t", @outvals)."\n"; 
   }

}


sub OLDPREPEPTIDE_combine_piligexpbits {

   my $in = shift ;
   my $pb = $in->{pb} ;
   my $ligbits = $in->{ligbits} ;
   my $liginfo = $in->{liginfo} ;
   my $pibits = $in->{pibits} ;
   my $expbits = $in->{expbits} ;

   my $class2bits ;
   my $class2ligs ;
   my $class2anligbits ;
   my $class2allligbits ;

   foreach my $classtype (qw/fam sf/) {
      foreach my $sid (keys %{$ligbits->{$classtype}}) {
         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;

         foreach my $ligsig (keys %{$ligbits->{$classtype}->{$sid}}) {
            if ($ligsig eq 'cumulative') {next;}
            push @{$class2allligbits->{$classtype}->{$class}},
               [$ligsig, $ligbits->{$classtype}->{$sid}->{$ligsig}] ; }

         if (exists $pibits->{$classtype}->{$sid}) {
            foreach my $ligsig (keys %{$ligbits->{$classtype}->{$sid}}) {
               if ($ligsig eq 'cumulative') {next;}
               my ($pdb, $ligcod, $ligid) = split(/\t/, $ligsig) ;
               my $anligbits =
                  $ligbits->{$classtype}->{$sid}->{$ligsig}->Shadow() ;

               $anligbits->AndNot(
                  $ligbits->{$classtype}->{$sid}->{$ligsig},
                  $pibits->{$classtype}->{$sid}->{'P'}) ;

               if ($anligbits->Norm() > 0) {
                  push @{$class2anligbits->{$classtype}->{$class}},
                     [$ligsig, $anligbits] ;

                  my @btypes = ("L") ;
                  if (exists $liginfo->{mwinrange}->{$ligcod}) {
                     push @btypes, "Lmwinrange"; }

                  foreach my $btype (@btypes) {
                     $class2ligs->{$classtype}->{$class}->{$btype}->{$ligsig}++;
                     if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
                        $class2bits->{$classtype}->{$class}->{$btype} = 
                           $anligbits->Clone() ;
                     } else {
                        $class2bits->{$classtype}->{$class}->{$btype}->Union(
                          $class2bits->{$classtype}->{$class}->{$btype},
                          $anligbits);
                     }
                  }
               }
            }

            foreach my $btype (keys %{$pibits->{$classtype}->{$sid}}) {
               my $anpibits =
                  $pibits->{$classtype}->{$sid}->{$btype}->Shadow() ;

               $anpibits->AndNot(
                  $pibits->{$classtype}->{$sid}->{$btype},
                  $ligbits->{$classtype}->{$sid}->{cumulative}) ;

               if ($anpibits->Norm() > 0) {
                  if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
                     $class2bits->{$classtype}->{$class}->{$btype} = 
                        $anpibits->Clone() ;
                  } else {
                     $class2bits->{$classtype}->{$class}->{$btype}->Union(
                        $class2bits->{$classtype}->{$class}->{$btype}, $anpibits) ;
                  }
               }
            }

         } else {
            foreach my $ligsig (keys %{$ligbits->{$classtype}->{$sid}}) {
#               print STDERR "   now on ligand $ligsig\n" ;
               if ($ligsig eq 'cumulative') {next;}
               my ($pdb, $ligcod, $ligid) = split(/\t/, $ligsig) ;
               my @btypes = ("L") ;
               if (exists $liginfo->{mwinrange}->{$ligcod}) {
                  push @btypes, "Lmwinrange"; }

               foreach my $btype (@btypes) {
                  $class2ligs->{$classtype}->{$btype}->{$ligsig}++ ;
                  if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
                     $class2bits->{$classtype}->{$class}->{$btype} = 
                        $ligbits->{$classtype}->{$sid}->{$ligsig}->Clone() ;
                  } else {
                     $class2bits->{$classtype}->{$class}->{$btype}->Union(
                        $class2bits->{$classtype}->{$class}->{$btype},
                        $ligbits->{$classtype}->{$sid}->{$ligsig}) ;
                  }
               }
            }
         }
      }

      foreach my $sid (keys %{$pibits->{$classtype}}) {
         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;
         if (!exists $ligbits->{$classtype}->{$sid}) {
            foreach my $btype (keys %{$pibits->{$classtype}->{$sid}}) {
               if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
                  $class2bits->{$classtype}->{$class}->{$btype} = 
                     $pibits->{$classtype}->{$sid}->{$btype}->Clone() ;
               } else {
                  $class2bits->{$classtype}->{$class}->{$btype}->Union(
                     $class2bits->{$classtype}->{$class}->{$btype},
                     $pibits->{$classtype}->{$sid}->{$btype}) ;
               }
            }
         }
      }
   }

   foreach my $classtype (keys %{$expbits}) {
      foreach my $class (keys %{$expbits->{$classtype}}) {
         $class2bits->{$classtype}->{$class}->{E} =
            $expbits->{$classtype}->{$class}->Clone() ;
      }
   }

   return {
      class2bits => $class2bits,
      class2ligs => $class2ligs,
      class2anligbits => $class2anligbits,
      class2allligbits => $class2allligbits,
   } ;
   
}


sub combine_piligpepexpbits {

   my $in = shift ;
   my $pb = $in->{pb} ;
   my $ligbits = $in->{ligbits} ;
   my $liginfo = $in->{liginfo} ;
   my $pibits = $in->{pibits} ;
   my $pepbits = $in->{pepbits} ;
   my $expbits = $in->{expbits} ;

   my $class2bits ;
   my $class2ligs ;
   my $class2allligbits ;

   foreach my $classtype (qw/fam sf/) {
      foreach my $sid (keys %{$ligbits->{$classtype}}) {
         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;

         foreach my $ligsig (keys %{$ligbits->{$classtype}->{$sid}}) {
            if ($ligsig eq 'cumulative') {next;}
            $class2ligs->{$classtype}->{$class}->{$ligsig}++;

# fpd090303_2034: changed the elements of this array from (ligsig, ligbits) to
# (ligsig, ligbits, sid) - have to change collate_perfam and collate_perinst
# to respect this

            push @{$class2allligbits->{$classtype}->{$class}},
               [$ligsig,
                $ligbits->{$classtype}->{$sid}->{$ligsig}, $sid] ;

            my ($pdb, $ligcod, $ligid) = split(/\t/, $ligsig) ;

            $class2ligs->{$classtype}->{$class}->{$ligsig}++;
            if (!exists $class2bits->{$classtype}->{$class}->{'L'}) {
               $class2bits->{$classtype}->{$class}->{'L'} =
                  $ligbits->{$classtype}->{$sid}->{$ligsig}->Clone() ;
            } else {
               $class2bits->{$classtype}->{$class}->{'L'}->Union(
                  $class2bits->{$classtype}->{$class}->{'L'},
                  $ligbits->{$classtype}->{$sid}->{$ligsig}) ;
            }
         }
      }

      foreach my $sid (keys %{$pibits->{$classtype}}) {
         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;
         foreach my $btype (keys %{$pibits->{$classtype}->{$sid}}) {
            if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
               $class2bits->{$classtype}->{$class}->{$btype} =
                  $pibits->{$classtype}->{$sid}->{$btype}->Clone() ;
            } else {
               $class2bits->{$classtype}->{$class}->{$btype}->Union(
                  $class2bits->{$classtype}->{$class}->{$btype},
                  $pibits->{$classtype}->{$sid}->{$btype}) ;
            }

            if (!exists $class2bits->{$classtype}->{$class}->{'P'}) {
               $class2bits->{$classtype}->{$class}->{'P'} =
                  $pibits->{$classtype}->{$sid}->{$btype}->Clone() ;
            } else {
               $class2bits->{$classtype}->{$class}->{'P'}->Union(
                  $class2bits->{$classtype}->{$class}->{'P'},
                  $pibits->{$classtype}->{$sid}->{$btype}) ;
            }
         }
      }

      foreach my $sid (keys %{$pepbits->{$classtype}}) {
         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;
         foreach my $targetch (keys %{$pepbits->{$classtype}->{$sid}}) {
            foreach my $btype ('p', 'P') {
               if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
                  $class2bits->{$classtype}->{$class}->{$btype} = 
                     $pepbits->{$classtype}->{$sid}->{$targetch}->Clone() ;
               } else {
                  $class2bits->{$classtype}->{$class}->{$btype}->Union(
                     $class2bits->{$classtype}->{$class}->{$btype},
                     $pepbits->{$classtype}->{$sid}->{$targetch}) ;
               }
            }
         }
      }
   }

   foreach my $classtype (keys %{$expbits}) {
      foreach my $class (keys %{$expbits->{$classtype}}) {
         $class2bits->{$classtype}->{$class}->{E} =
            $expbits->{$classtype}->{$class}->Clone() ;
      }
   }

   return {
      class2bits => $class2bits,
      class2ligs => $class2ligs,
      class2allligbits => $class2allligbits,
   } ;
   
}

sub _identify_peptide_and_nucleic_chains {

   my $in = shift ;
   my $fn = $in->{subsetsres_fn} ;
   my $pb = $in->{pb} ;
   my $bdp_id = $in->{bdp_id} ;

   my ($t_subset_id, $t_chain_no, $t_chain_id, $t_resno_serial, $t_resno) =
      pibase::rawselect_metatod($fn,
      "SELECT subset_id, chain_no, chain_id, resno_serial, resno FROM $fn") ;

# ONLY CONSIDERING SCOP HERE: peptide chain is subset source specific
# - ie a SCOP peptide chain may not be a CATH peptide chain

   my $nucleic_chains ;
   my $res2sid ;
   my $chain2domains ;
   my $chain_list ;
   my $chain_type ;
   my $chain2chaindom ;
   my $chaindom2chain ;
   foreach my $j ( 0 .. $#{$t_subset_id}) {
      my $sig = $t_resno->[$j]."\n".$t_chain_id->[$j] ;
      $chain_list->{$t_chain_no->[$j]}++ ;

# without the type check, will let in DNA chain ``domains'' as well...

      if ($t_subset_id->[$j] =~ /CHAIN/) {

         if (!exists $pb->{chain_2_type}->{$bdp_id}->{$t_chain_no->[$j]}) {
            print STDERR "WARNING: can't find chain type for $bdp_id ".
                         $t_chain_no->[$j]."\n" ;
         }

         if ($pb->{chain_2_type}->{$bdp_id}->{$t_chain_no->[$j]} eq 'n') {
            $chain_type->{$t_chain_no->[$j]} = 'n' ;
         } else {
            $chain_type->{$t_chain_no->[$j]} = 'p' ;
         }

         $chain2chaindom->{$t_chain_no->[$j]} = $t_subset_id->[$j] ;
         $chaindom2chain->{$t_subset_id->[$j]} = $t_chain_no->[$j] ;

      } elsif ($t_subset_id->[$j] =~ /SCOP/ &&
               $pb->{sid2class}->{fam}->{$t_subset_id->[$j]} !~ /^j/) {

# assumes only 1 SCOP domain membership...
         $res2sid->{$sig} = $t_subset_id->[$j] ;
         $chain2domains->{$t_chain_no->[$j]}->{$t_subset_id->[$j]}++ ;
      }
   }

   my $numchains_withdoms = keys %{$chain2domains} ;
   if ($numchains_withdoms == 0) {return; }

   my $peptide_chains ;
   foreach my $chain (keys %{$chain_list}) {
      if (!exists $chain2domains->{$chain}) {
         $peptide_chains->{$chain}++ ; } }

   return {
      res2sid => $res2sid,
      chain_type => $chain_type,
      peptide_chains => $peptide_chains,
      chain2domains => $chain2domains,
      chain2chaindom => $chain2chaindom,
      chaindom2chain => $chaindom2chain,
   }

}


sub OLD_pilig_split_ligs {
   my $in = shift ;

   my $lb = $in->{lb} ;
   my $total_numpdb = keys %{$lb->{pdb2res2ligid}} ;

   my $splitlines = POSIX::ceil($total_numpdb / $in->{clustspecs}->{numjobs}) ;
   my @splitfiles ;

   my $cur_splitnum = 1 ;
   my $cur_fn = "split.".$in->{out_fn_prefix}.".$cur_splitnum" ;
   push @splitfiles, $cur_fn ;
   open(OUTF, ">$in->{splits_dir}/$cur_fn") ;

   my $tasklist = "'$cur_fn'" ;

   my $num_pdb = 0;
   foreach my $pdb (keys %{$lb->{pdb2res2ligid}}) {
      if ($num_pdb > $splitlines) {
         $cur_splitnum++ ;
         $cur_fn = "split.".$in->{out_fn_prefix}.".$cur_splitnum" ;
         close(OUTF) ;
         open(OUTF, ">$in->{splits_dir}/$cur_fn") ;
         push @splitfiles, $cur_fn ;
         $tasklist .= " '$cur_fn'" ;
         $num_pdb = 0 ;
      }

      foreach my $res ( keys %{$lb->{pdb2res2ligid}->{$pdb}}) {
         my ($resno, $resch) = split(/\n/, $res) ;
         foreach my $lig ( keys %{$lb->{pdb2res2ligid}->{$pdb}->{$res}}) {
            my ($ligcod, $ligid) = split(/\n/, $lig) ;
            my @outvals = ($pdb, $resno, $resch, $ligcod, $ligid) ;
            print OUTF join("\t", @outvals)."\n";
         }
      }
      $num_pdb++ ;
   }
   close(OUTF) ;
   my $numjobs = $cur_splitnum ;

   my $perlscript_fn = __FILE__;

   open (SGEFH, ">$in->{SGE_fn}") ;
   print SGEFH "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $in->{SGE_dir}
#\$ -e $in->{SGE_dir}
#\$ -r y
$in->{clustspecs}->{nodespecs}
#\$ -p $in->{clustspecs}->{priority}
#\$ -t 1-$numjobs

set tasks1=( $tasklist )

set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/scratch/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $in->{splits_dir}/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn -assign_lig < \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;

   close(SGEFH) ;



#   return {
#      numjobs => $cur_splitnum,
#      splitfiles => \@splitfiles,
#      tasklist => $tasklist
#   } ;
}

sub OLD_pilig_split_pi {

   my $in = shift ;
   my $pb = $in->{pb} ;
   my $numbdp = keys %{$pb} ;

   my ($tempfh, $tempfn) = tempfile("bdp_ids.XXXXX") ;
   foreach my $bdp (keys %{$in->{pb}->{bdp2sid12}}) {
      foreach my $sid12 (keys %{$in->{pb}->{bdp2sid12}->{$bdp}}) {
         my ($sid1, $sid2) = split(/\t/, $sid12) ;
         my @outvals = ($bdp, $in->{pb}->{bdp2contactsfn}->{$bdp},
                        $sid1, $in->{pb}->{sid2class}->{fam}->{$sid1},
                        $sid2, $in->{pb}->{sid2class}->{fam}->{$sid2}) ;
         print {$tempfh} join("\t", @outvals)."\n";
      }
   }
   close($tempfh) ;

   my $splits = _clust_split_ins({
      fn => $tempfn,
      dir => $in->{splits_dir},
      numjobs => $in->{clustspecs}->{numjobs}
   });
   unlink $tempfn ;

   my $perlscript_fn = __FILE__;

   open (SGEFH, ">$in->{SGE_fn}") ;
   print SGEFH "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $in->{SGE_dir}
#\$ -e $in->{SGE_dir}
#\$ -r y
$in->{clustspecs}->{nodespecs}
#\$ -p $in->{clustspecs}->{priority}
#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )

set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/scratch/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $in->{splits_dir}/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn -assign_pi < \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;

   close(SGEFH) ;


#   return {
#      splitfiles => $splits->{splitfiles},
#      numjobs => $splits->{numjobs},
#      tasklist => $splits->{tasklist},
#   } ;

}


sub OLD_pilig__clust_split_ins {

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


sub OLD_pilig_parse_raf_line {

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

   return $res ;
}


sub OLD_pilig_raf_preload {
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


sub OLD_pilig_get_raf_line {
   my $in = shift ;

   my $pdb = $in->{pdb_id} ;
   my $ch = $in->{chain} ;
   my $fn = $in->{fn} ;

#   if ($ch eq '-') {$ch = '_';}
   if ($ch eq ' ') {$ch = '_';}
   my $pdbch = $pdb.$ch ;

   my $tcom = "grep '^$pdbch' $fn" ;
   my $line = `$tcom` ;
   chomp $line;
#   if ($pdb eq '1a0w') {
#      print STDERR "**grepped $tcom to get $line\n" ; }

   return $line ;
}


sub OLD_pilig_read_asteroids_aln {

   my $in = shift;
   my $fn;
   $fn->{aln} = $in->{aln_fn} ;
   $fn->{seq} = $in->{seq_fn} ;
   my $allchains = $in->{allchains} ;

   my $data;
   my $cur_dom = '' ;
   my $seq ;
   open(SEQF, $fn->{seq}) ;
#070503_1511   print STDERR "tried to open $fn->{seq}\n" ;
   while (my $line = <SEQF>) {
      chomp $line;
      if ($line =~ /^>/) {
         $line =~ s/^>// ;
         my ($t_dom, undef, undef, undef) = split(/\s/, $line) ;
         $data->{seq}->{$t_dom} = '';
         $cur_dom = $t_dom ;
      } else {
         $data->{seq}->{$cur_dom} .= $line ;
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
         $cur_dom = $t_dom ;
      } else {
         $data->{aln}->{$cur_dom} .= $line ;
      }
   }
   $data->{alnlength} = length($data->{aln}->{$cur_dom}) ;
   close(ALNF) ;

   my $seqweights ;
   foreach my $pos ( 0 .. ($data->{alnlength} - 1) ) {
#      my $freq = {
#         'A' => 0,
#         'R' => 0,
#         'N' => 0,
#         'D' => 0,
#         'C' => 0,
#         'Q' => 0,
#         'E' => 0,
#         'G' => 0,
#         'H' => 0,
#         'I' => 0,
#         'L' => 0,
#         'K' => 0,
#         'M' => 0,
#         'F' => 0,
#         'P' => 0,
#         'S' => 0,
#         'T' => 0,
#         'W' => 0,
#         'Y' => 0,
#         'V' => 0,
#         '-' => 0,
#      } ;

      my $freq ;

      my $numdoms = keys %{$data->{aln}} ;
      my @tp ;
      foreach my $dom (keys %{$data->{aln}}) {
         push @tp, substr($data->{aln}->{$dom}, $pos, 1) ;
         $freq->{substr($data->{aln}->{$dom}, $pos, 1)}++ ; }

      $data->{meta}->{numrestypes}->{$pos} = keys %{$freq} ;
      $data->{meta}->{shannone}->{$pos} = 0 ;
      foreach my $res (keys %{$freq}) {
         $freq->{$res} /= $numdoms ;
         $data->{meta}->{shannone}->{$pos} += $freq->{$res} *
                                       (log($freq->{$res}) / log(21)) ;
      }
#      if ($data->{shannone}->{$pos} == 0) {
#         print STDERR " shannone pos $pos $in->{aln_fn} = 0 charc:".join(",", @tp)."\n" ; }

   }

   return $data ;

}


sub OLD_pilig_parse_aln_raf {

   my $in = shift;
   my $data = $in->{alndata} ;
   my $raf = $in->{raf} ;

#data->{seq} has raw FASTA sequence strings - check substr(curpos) to see
# if x is capital = frag or lower case = residue unk

   my $DEBUGALN = 0 ;
   my ($alnpos2atomresno, $atomresno2alnpos);
   foreach my $dom (keys %{$data->{defstring}}) {

      if ($DEBUGALN) {
         print STDERR "now on $dom\n" ;
         print STDERR "   defstring:\t$data->{defstring}->{$dom}\n" ;
         print STDERR "   alnstring:\t$data->{aln}->{$dom}\n" ;
         print STDERR "   pdb:\t$data->{pdb}->{$dom}\n" ;
         print STDERR "   class:\t$data->{class}->{$dom}\n" ;
      }

      my $gs2atomresno ;
      my $gs2resna ; #to verify correct positioning when parsing alignment

#get gsresno -> atomresno mapping from raf file
      my $gsresno = -1 ;
      foreach my $fragdef ( @{$data->{frags}->{$dom}}) {
         $gsresno++ ; #to take care of X in between fragments
         my $tch = $fragdef->{chain} ; if ($tch eq ' ') {$tch = '_';}
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
               if ($DEBUGALN) {print STDERR " on dom $dom: lookahead b find for $fragdef->{b}\n" ;}
               $ind_b = $rafinfo->{atomresno2ind_ahead}->{$fragdef->{b}} ;
            }

            if (defined $rafinfo->{atomresno2ind}->{$fragdef->{e}}) {
               $ind_e = $rafinfo->{atomresno2ind}->{$fragdef->{e}} ;
            } else {
               if ($DEBUGALN) {print STDERR " on dom $dom: lookback e find for $fragdef->{e}\n" ;}
               $ind_e = $rafinfo->{atomresno2ind_back}->{$fragdef->{e}} ;
            }

         } else {


            if (defined $rafinfo->{atomresno2ind}->{$rafinfo->{atomresno_first}}) {
               $ind_b= $rafinfo->{atomresno2ind}->{$rafinfo->{atomresno_first}};
            } else {
               if ($DEBUGALN) {print STDERR " on dom $dom: lookahead b find for $rafinfo->{atomresno_first}\n" ;}
               $ind_b= $rafinfo->{atomresno2ind_ahead}->{$rafinfo->{atomresno_first}} ;
            }

            if (defined $rafinfo->{atomresno2ind}->{$rafinfo->{atomresno_last}}) {
               $ind_e = $rafinfo->{atomresno2ind}->{$rafinfo->{atomresno_last}};
            } else {
               if ($DEBUGALN) {print STDERR " on dom $dom: lookback e find for $rafinfo->{atomresno_last}\n" ;}
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

#go through aln entry and assign aln pos -> atomresno mapping

      my $curseqpos = 0 ;
      foreach my $pos (0 .. (length($data->{aln}->{$dom}) - 1)) {
         my $alnchar = substr($data->{aln}->{$dom}, $pos, 1) ;

#         if (DEBUGALN && $dom eq DEBUGDOM) {
#            print STDERR "dom $dom (alnpos $pos, curseqpos $curseqpos) $alnchar\n" ; }

         if ($alnchar eq 'X' &&
            (substr($data->{seq}->{$dom}, ($curseqpos), 1) eq 'X')) {

            if ($DEBUGALN) {
               print STDERR "$dom: real frag $alnchar at seqpos $curseqpos (alnpos $pos)\n" ;}

            $curseqpos++ ;

         } elsif ($alnchar ne '-') {
            if ($DEBUGALN) {
               if ($alnchar eq 'X') {print STDERR "$dom caught unk res: alnchar $alnchar at seq pos $curseqpos (alnpos $pos)\n";}
            }

#            if ($alnchar eq 'X') {print STDERR "$dom caught unk res: alnchar $alnchar at seq pos $curseqpos (alnpos $pos) TRYING A FIX\n";next;}

            if (exists $gs2atomresno->{$curseqpos}) {

               if ($alnchar ne $gs2resna->{$curseqpos}) {
#                  print "ERROR $dom: alignment position $pos ($alnchar) mismatched with gs sequence position $curseqpos ($gs2resna->{$curseqpos}\n" ;
                  print STDERR "ERROR $dom: alignment position $pos ($alnchar) mismatched with gs sequence position $curseqpos ($gs2resna->{$curseqpos}\n" ;
               }

#               if (DEBUGALN && $dom eq DEBUGDOM) {
#                     print STDERR "$dom\t$pos\t$alnchar\t$gs2atomresno->{$curseqpos}\n" ;}

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


sub OLD_pilig_load_asteroids_aln {

   my $in = shift ;

   my $gdseqh = $in->{gdseqh} ;
   my $seqclcont100 = $in->{seqclcont100} ;
   my $seqcl100 = $in->{seqcl100} ;

   my $data = read_asteroids_aln({
      aln_fn => $in->{aln_fn},
      seq_fn => $in->{seq_fn},
      allchains => $in->{allchains}
   });

   my $aln = parse_aln_raf({
      alndata => $data,
      raf => $in->{raf}
   }) ;
   $aln->{alnlength} = $data->{alnlength} ;
   $aln->{meta} = $data->{meta} ;

   return $aln ;

}


sub OLD_pilig_load_astral_headers {

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

sub OLD_pilig_load_astral_clusters {

   my $in = shift ;
   my $out = $in->{out} ;
   my $fn = $in->{fn} ;

   my $seqcl ;
   my $seqcl2cont ;
   my $clusnum2rep;
   my $rep2clusnum;
   print STDERR "load ASTRAL sequence clusters: " ;
   {
      foreach my $seqid (keys %{$fn->{seqcl}}) {
         open (ASTRALCL, $fn->{seqcl}->{$seqid}) ;
         my $clusnum = 0 ;
         while (my $line = <ASTRALCL>) {
            chomp $line;
            if ($line =~ /^Rep/) {
               $line =~ s/^Rep: //g;
               $clusnum++ ;
               $line =~ s/^\s+// ;
               my $scopid = substr($line, 0, 7) ;
               $clusnum2rep->{$clusnum} = $scopid ;
               $rep2clusnum->{$scopid} = $clusnum;
            }

            if ($line =~ /score/) {
               $line =~ s/^\s+// ;
               my $scopid = substr($line, 0, 7) ;
               $seqcl->{$seqid}->{$scopid} = $clusnum2rep->{$clusnum} ;
               push @{$seqcl2cont->{$seqid}->{$clusnum2rep->{$clusnum}}},
                  $scopid ;
            }
         }
         close(ASTRALCL) ;
      }
   }
   print STDERR "X\n" ;

   $out->{seqcl} = $seqcl ;
   $out->{clusnum2rep} = $clusnum2rep ;
   $out->{seqcl2cont} = $seqcl2cont ;

}

sub pilig_htmlsummary_main {
#1pgj BDP23163-0_SCOP.d1pgja1 BDP23163-0_SCOP.d1pgjb1 a.100.1.1 a.100.1.1 1pgp:6PG:001 133 8 5 0.037 0.038 0.625

   print "<table border=\"1\">\n" ;
   my $count =1 ;
   while (my $line = <STDIN>) {
      chomp $line;
      my @t = split(' ', $line);
      my ($pdb, $sid1, $sid2, $class1, $class2, $ligsig, $nump, $numl, $numlandp, $obi, $opi, $olig) = @t ;

      my ($ligpdb,$ligcod,$lignum) = split(/\:/, $ligsig) ;

      my @outvals = (
         $count,
         '<a href="http://www.rcsb.org/pdb/explore/explore.do?structureId='.
            uc($pdb).'">'.$pdb.'</a>',
         $sid1, $sid2, $class1, $class2,
         '<a href="http://www.rcsb.org/pdb/explore/explore.do?structureId='.
            uc($ligpdb).'">'.$ligpdb.'</a>',
         $ligcod,$lignum, $nump, $numl, $numlandp, $obi, $opi, $olig
      ) ;

      if ($class1 eq $class2) {
         next;
         print "<tr>\n" ;
      } else {
         print "<tr bgcolor=\"lightblue\">\n" ;
      }

      foreach my $j ( 0 .. $#outvals) {
         print "<td>$outvals[$j]</td>" ;
      }
      print "</tr>\n" ;
      $count++ ;
   }
   print "</table>\n" ;
   
}


sub readin_superfam_fxn {

   my $in ;
   my $pilig_specs = set_pilig_specs() ;

   my $sf2fxn; my $fxn2broad ; my $sf2fxn_broad;
   open(SUPERFAM_FXN, $pilig_specs->{superfam_fxn}) ;
   while (my $line = <SUPERFAM_FXN>) {
      chomp $line ;
      my ($fxn, $sunid, $classtype, $class, undef,$name) = split(/\t/, $line) ;
      if ($fxn eq "S'") {$fxn = 'S';}
      $sf2fxn->{$class} = $fxn ;
   }
   close(SUPERFAM_FXN) ;

   open(SUPERFAM_FXN_CATEGORIES, $pilig_specs->{superfam_fxn_categories}) ;
   while (my $line = <SUPERFAM_FXN_CATEGORIES>) {
      chomp $line ;
      my ($broad, $fxn, $fxn_name, $fxn_descr) = split(/\t/, $line) ;
      if ($fxn eq "S'") {$fxn = 'S';}
      $fxn2broad->{$fxn} = $broad;
   }
   close(SUPERFAM_FXN_CATEGORIES) ;

   foreach my $sf (keys %{$sf2fxn}) {
      if (!exists $fxn2broad->{$sf2fxn->{$sf}}) {
         $sf2fxn_broad->{$sf} = 'N_A' ;
         print STDERR "WARNING: Fxn class for $sf ".
                      "(".$sf2fxn->{$sf}.") not found\n";
      } else {
         $sf2fxn_broad->{$sf} = $fxn2broad->{$sf2fxn->{$sf}} ;
      }
   }

   return {
      sf2fxn => $sf2fxn,
      fxn2broad => $fxn2broad,
      sf2fxn_broad => $sf2fxn_broad,
   }

}


sub make_paper_figures {

   require Bit::Vector;

#   require R;
#   require RReferences ;
#   R::initR("--silent") ;

   my $pilig_specs = set_pilig_specs() ;
   my $pibase_specs = pibase::get_specs() ;
   my $results ;

   my $superfam_fxn = readin_superfam_fxn() ;

   my $input_fn = {
      "PINT_LINT" => "collate_perinstance_PINT.20090319.out",
      "pINT_LINT" => "collate_perinstance_pINT.20090319.out",
      "perfam" => "collate_perfam.20090319.out",
      "consscore" => "outconsscore.31502.out"
   } ;

# 1. number of families/domains/interactions (in range (right?) - ie SCOP classes a-g; ie number of families with ASTEROIDS alignments)
# number of alignments in /groups/eddy/home/davisf/work/pibase/pibase200708/external_data/asteroids/1.73/alignments/fam

   my $dir_astral_fam_aln = $pibase_specs->{asteroids}->{'fam_aln'};
   my @all_fams = glob($dir_astral_fam_aln.'/*aln') ;
   my $num_fam = $#all_fams + 1 ;
   $results->{num_families}->{total} = $num_fam ;

   my $pb = _pilig_tod_pibase_preload() ;

# 2. number of families/domains/interactions that bind in range small molecules

   my $liginfo = _pilig_load_liginfo() ;
   my $standardres =  _list_standardres();
   my $standardres_20aa =  _list_standardres_20aa();

   my $class2alnlength ;
   my $ligbits = readin_ligassignments({
      fn => $pilig_specs->{outfiles}->{assign_lig},
      liginfo => $liginfo,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   });
   {
      foreach my $sid (keys %{$ligbits->{fam}}) {
         my $cur_class = $pb->{sid2class}->{fam}->{$sid} ;
         $results->{families}->{bind_ligands}->{$cur_class}++ ;
         foreach my $ligsig (keys %{$ligbits->{fam}->{$sid}}) {
            if ($ligsig eq 'cumulative') {next;}
            $results->{interactions}->{'all'}->{l}->{$sid."\t".$ligsig}++; }
      }
      $results->{num_families}->{bind_ligands} = 
         keys %{$results->{families}->{bind_ligands}} ;

      $results->{num_interactions}->{'all'}->{l} =
         keys %{$results->{interactions}->{'all'}->{l}} ;
   }

   my $lig_clus_assignments = readin_ligassignments({
      fn => $pilig_specs->{outfiles}->{assign_lig_clusters},
      liginfo => $liginfo,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
      clustrep_fl => 1,
   });
   {
      my $curbits = $lig_clus_assignments->{fam} ;
      foreach my $sid (keys %{$curbits}) {
         foreach my $ligsig (keys %{$curbits->{$sid}}) {
            if ($ligsig eq 'cumulative') {next;}
            $results->{interactions}->{'clustered'}->{l}->{$sid."\t".$ligsig}++ ;}}
      $results->{num_interactions}->{'clustered'}->{l} =
         keys %{$results->{interactions}->{'clustered'}->{l}} ;
   }


# 3. number of families/domains/interactions that bind peptides
   my $pepnucibits = readin_pepnuciassignments({
      fn => $pilig_specs->{outfiles}->{assign_pepnuci},
      pb => $pb,
      pilig_specs => $pilig_specs,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   });
   {
      foreach my $sid (keys %{$pepnucibits->{pep}->{fam}}) {
         my $cur_class = $pb->{sid2class}->{fam}->{$sid} ;
         $results->{families}->{bind_pep}->{$cur_class}++ ;
         $results->{families}->{bind_protein}->{$cur_class}++ ;
      }
      $results->{num_families}->{bind_pep} =
         keys %{$results->{families}->{bind_pep}} ;
   }


# 4. number of families/domains/interactions that bind domains (inter)
# 5. number of families/domains/interactions that bind domains (intra)


   my $piassignments= readin_piassignments({
      fn => $pilig_specs->{outfiles}->{assign_pi},
      pb => $pb,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
      dont_read_jdomains_fl => 1
   });
   my $pibits = $piassignments->{pibits} ;
   {
      foreach my $sid (keys %{$pibits->{fam}}) {
         my $cur_class = $pb->{sid2class}->{fam}->{$sid} ;
         foreach my $btype (keys %{$pibits->{fam}->{$sid}}) {
            $results->{families}->{"bind_$btype"}->{$cur_class}++ ;
            $results->{families}->{bind_protein}->{$cur_class}++ ;
         }
      }
      $results->{num_families}->{bind_Pinter} =
         keys %{$results->{families}->{bind_Pinter}} ;

      $results->{num_families}->{bind_Pintra} =
         keys %{$results->{families}->{bind_Pintra}} ;
   }

   my $pi_clus_assignments= readin_piassignments({
      fn => $pilig_specs->{outfiles}->{assign_pi_clusters},
      clustrep_fl => 1,
      pb => $pb,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
      dont_read_jdomains_fl => 1
   });

   foreach my $pi_type (qw/all clustered/) {
      my $cur_interfaces ;
      if ($pi_type eq 'all') {
         $cur_interfaces = $piassignments->{interfaces} ;
      } else {
         $cur_interfaces = $pi_clus_assignments->{interfaces} ;
      }

      foreach my $sid12 (keys %{$cur_interfaces}) {
         my $chains = $cur_interfaces->{$sid12}->{chains} ;
         my $cur_type = 'Pinter';
         if ($chains eq 'same') { $cur_type = 'Pintra'; }
         $results->{interactions}->{$pi_type}->{$cur_type}->{$sid12}++ ;
      }
   }

   foreach my $pi_type ( sort keys %{$results->{interactions}}) {
      foreach my $int_type (sort keys %{$results->{interactions}->{$pi_type}}) {
         $results->{num_interactions}->{$pi_type}->{$int_type} = 
            keys %{$results->{interactions}->{$pi_type}->{$int_type}} ;
         print "Number of $pi_type interactions that $int_type: ".
            $results->{num_interactions}->{$pi_type}->{$int_type}."\n" ; }}

   my $pepnuci_clus_bits = readin_pepnuciassignments({
      fn => $pilig_specs->{outfiles}->{assign_pepnuci_clusters},
      pb => $pb,
      clustrep_fl => 1,
      pilig_specs => $pilig_specs,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   });
   foreach my $pi_type (qw/all clustered/) {
      my $curbits ;
      if ($pi_type eq 'all') {
         $curbits = $pepnucibits->{pep}->{fam} ;
      } else {
         $curbits = $pepnuci_clus_bits->{pep}->{fam} ;
      }
      
      foreach my $sid (keys %{$curbits}) {
         foreach my $chain (keys %{$curbits->{$sid}}) {
            if ($chain eq 'cumulative') {next;}
            $results->{interactions}->{$pi_type}->{p}->{$sid."\t".$chain}++ ;}}

      $results->{num_interactions}->{$pi_type}->{p} =
         keys %{$results->{interactions}->{$pi_type}->{p}} ;
   }

   print "\n\n\n" ;
   foreach my $pi_type ( sort keys %{$results->{interactions}}) {
      foreach my $int_type (sort keys %{$results->{interactions}->{$pi_type}}) {
         $results->{num_interactions}->{$pi_type}->{$int_type} = 
            keys %{$results->{interactions}->{$pi_type}->{$int_type}} ;
         print "Number of $pi_type interactions that $int_type: ".
            $results->{num_interactions}->{$pi_type}->{$int_type}."\n" ; }}


# 3b. number of families/domains/interactions that bind peptides and small molecules
# 4b. number of families/domains/interactions that bind domains (inter) and small molecules
# 5b. number of families/domains/interactions that bind domains (intra) and small molecules
   $results->{num_families}->{bind_protein} =
      keys %{$results->{families}->{bind_protein}} ;

   foreach my $query_type (qw/bind_protein bind_pep bind_Pinter bind_Pintra/) {
      foreach my $ligand_fam (keys %{$results->{families}->{bind_ligands}}){
         if (exists $results->{families}->{$query_type}->{$ligand_fam}) {
            $results->{families_overlap}->{$query_type."_and_bind_ligands"}->{$ligand_fam}++ ; }
         $results->{num_families}->{$query_type."_and_bind_ligands"} =
keys %{$results->{families_overlap}->{$query_type."_and_bind_ligands"}}
      }
   }

   foreach my $fam_type ( sort keys %{$results->{num_families}}) {
      print "Number of families that ".$fam_type.": ".
         $results->{num_families}->{$fam_type}."\n" ; }

#   open(PERFAM, "collate_perfam_cluster.20090204.out") ;
   open(PERFAM, $input_fn->{perfam}) ;
   my $field2i ;
   my $histograms ;
   my $fam2pval;
   my $fam_fxncounts ;
   my $fam_fxnbroadcounts ;
   while (my $line = <PERFAM>) {
      chomp $line ;

      if ($line =~ /^\#/) {
         if ($line =~ /^#[0-9]/) { next; }
         $line =~ s/\#// ;
         my @t = split(/\t/, $line) ;
         foreach my $j ( 0 .. $#t) {
            $field2i->{$t[$j]} = $j ;
         }
         next;
      }

      my @t = split(/\t/, $line) ;

      my $curfam = $t[$field2i->{class}] ;

      $results->{collate_perfam}->{families}->{total}->{$curfam}++ ;

      if ($line =~ /NOINTERACTIONS/) {
         $results->{collate_perfam}->{families}->{"No interactions"}->{$curfam}++ ;
         next ;
      }

      if ($t[$field2i->{norm_L}] > 0) {
         $results->{collate_perfam}->{num_families}->{"Bind ligands 250-1000MW"}++ ;
         $results->{collate_perfam}->{families}->{"Bind ligands 250-1000MW"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_L}] >= 5) {
         $results->{collate_perfam}->{num_families}->{"Bind ligands 250-1000MW, >=5 positions"}++ ;
         $results->{collate_perfam}->{families}->{"Bind ligands 250-1000MW, >=5 positions"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_P}] > 0) {
         $results->{collate_perfam}->{num_families}->{"Bind proteins"}++ ;
         $results->{collate_perfam}->{families}->{"Bind proteins"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_Pintra}] > 0) {
         $results->{collate_perfam}->{num_families}->{"Bind domains(intra)"}++ ;
         $results->{collate_perfam}->{families}->{"Bind domains (intra)"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_Pinter}] > 0) {
         $results->{collate_perfam}->{num_families}->{"Bind domains(inter)"}++ ;
         $results->{collate_perfam}->{families}->{"Bind domains (inter)"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_p}] > 0) {
         $results->{collate_perfam}->{num_families}->{"Bind peptides"}++ ;
         $results->{collate_perfam}->{families}->{"Bind peptides"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_L}] > 0 &&
          $t[$field2i->{norm_P}] > 0) {

         my ($cursf) = ($curfam =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
         if (!exists $superfam_fxn->{sf2fxn_broad}->{$cursf}) {
            print STDERR "WARNING: no function found for $curfam ($cursf)\n" ;
         } else {
            $fam_fxnbroadcounts->{all_fam}->{$superfam_fxn->{sf2fxn_broad}->{$cursf}}++;
            $fam_fxnbroadcounts->{totals_all_fam}++ ;

            $fam_fxncounts->{all_fam}->{$superfam_fxn->{sf2fxn}->{$cursf}}++;
            $fam_fxncounts->{totals_all_fam}++ ;
         }


         $fam2pval->{over}->{$t[$field2i->{class}]} = $t[$field2i->{pval_L_P}];
         $fam2pval->{under}->{$t[$field2i->{class}]} =
            $t[$field2i->{pval_L_P_less}];
         $results->{collate_perfam}->{num_families}->{"Bind proteins and ligands"}++ ;
         $results->{collate_perfam}->{families}->{"Bind proteins and ligands"}++ ;
      }

      if ($t[$field2i->{norm_L_and_P}] >= 5) {
         $results->{collate_perfam}->{num_families}->{"Bind proteins and ligands, >=5 multifunc positions"}++ ;
         $results->{collate_perfam}->{families}->{"Bind proteins and ligands, >=5 multifunc positions"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_L_and_Pinter}] >= 5) {
         $results->{collate_perfam}->{num_families}->{"Bind domains(inter) and ligands, >=5 multifunc positions"}++ ;
         $results->{collate_perfam}->{families}->{"Bind domains(inter) and ligands, >=5 multifunc positions"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_L_and_Pintra}] >= 5) {
         $results->{collate_perfam}->{num_families}->{"Bind domains(intra) and ligands, >=5 multifunc positions"}++ ;
         $results->{collate_perfam}->{families}->{"Bind domains(intra) and ligands, >=5 multifunc positions"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_L_and_p}] >= 5) {
         $results->{collate_perfam}->{num_families}->{"Bind peptides and ligands, >=5 multifunc positions"}++ ;
         $results->{collate_perfam}->{families}->{"Bind peptides and ligands, >=5 multifunc positions"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_L}] > 0 &&
          $t[$field2i->{norm_Pintra}] > 0) {
         $results->{collate_perfam}->{families}->{"Bind domains(intra) and ligands"}->{$curfam}++;
         $results->{collate_perfam}->{num_families}->{"Bind domains(intra) and ligands"}++;
      }

      if ($t[$field2i->{norm_L}] > 0 &&
          $t[$field2i->{norm_Pinter}] > 0) {
         $results->{collate_perfam}->{families}->{"Bind proteins(inter) and ligands"}->{$curfam}++;
         $results->{collate_perfam}->{num_families}->{"Bind proteins(inter) and ligands"}++;
      }

      if ($t[$field2i->{norm_L}] > 0 &&
          $t[$field2i->{norm_p}] > 0) {
         $results->{collate_perfam}->{families}->{"Bind peptides and Ligands"}->{$curfam}++;
         $results->{collate_perfam}->{num_families}->{"Bind peptides and ligands"}++;
      }

      push @{$histograms->{num_alnpos_E}}, 
         $t[$field2i->{norm_E}] ;

#      push @{$histograms->{num_alnpos_E}}, 
#         $t[$field2i->{norm_E}] ;
#
#      push @{$histograms->{num_alnpos_Lmwinrange}},
#         $t[$field2i->{norm_Lmwinrange}];
#
#      push @{$histograms->{num_alnpos_Pinter}},
#         $t[$field2i->{norm_Pinter}];
#
#      push @{$histograms->{num_alnpos_Pintra}},
#         $t[$field2i->{norm_Pintra}];
#
#      push @{$histograms->{num_alnpos_P}},
#         $t[$field2i->{norm_P}];
#
#      push @{$histograms->{num_alnpos_Lmwinrange_Pinter}},
#         $t[$field2i->{norm_Lmwinrange_and_Pinter}];
#
#      push @{$histograms->{num_alnpos_Lmwinrange_Pintra}},
#         $t[$field2i->{norm_Lmwinrange_and_Pintra}];
#
#      push @{$histograms->{num_alnpos_Lmwinrange_P}},
#         $t[$field2i->{norm_Lmwinrange_and_P}];
#
#      push @{$histograms->{cons_shannone_L}},
#         $t[$field2i->{norm_Lmwinrange_and_P}];
   }
   $field2i = {} ;
   close(PERFAM) ;

   foreach my $fam (keys %{$fam2pval->{under}}) {
      my ($sf) = ($fam =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
      my $cur_fxn ; my $cur_fxn_broad ;
      if (exists $superfam_fxn->{sf2fxn_broad}->{$sf}) {
         $cur_fxn = $superfam_fxn->{sf2fxn}->{$sf} ;
         $cur_fxn_broad = $superfam_fxn->{sf2fxn_broad}->{$sf} ;
      }

      if ($fam2pval->{under}->{$fam} < 0.01) {
         $results->{collate_perfam}->{families}->{"pval0.01_under_overlap"}->{$fam}++ ;
         if (defined $cur_fxn) {
            $fam_fxnbroadcounts->{under}->{$cur_fxn_broad}++ ;
            $fam_fxnbroadcounts->{totals_under}++ ;

            $fam_fxncounts->{under}->{$cur_fxn}++ ;
            $fam_fxncounts->{totals_under}++ ;
         }
      }

      if ($fam2pval->{over}->{$fam} < 0.01) {
         $results->{collate_perfam}->{families}->{"pval0.01_over_overlap"}->{$fam}++ ;
         if (defined $cur_fxn) {
            $fam_fxnbroadcounts->{over}->{$cur_fxn_broad}++ ;
            $fam_fxnbroadcounts->{totals_over}++ ;

            $fam_fxncounts->{over}->{$cur_fxn}++ ;
            $fam_fxncounts->{totals_over}++ ;
         }
      }
   }

   my $fxn_propensities ;
   open(FXNPROP, ">fxn_propensities_detailed.sigunder_sigover.txt") ;
   foreach my $t_fxn ( sort keys %{$fam_fxncounts->{all_fam}}) {
#HERENOW090225_1953 
      foreach my $o_type (qw/under over/) {
         if (!exists $fam_fxncounts->{$o_type}->{$t_fxn}) {
            $fam_fxncounts->{$o_type}->{$t_fxn} = 0 ;
            $fxn_propensities->{$o_type}->{$t_fxn} = 0 ;
         } else {
            $fxn_propensities->{$o_type}->{$t_fxn} =
            ($fam_fxncounts->{$o_type}->{$t_fxn} / 
             $fam_fxncounts->{"totals_".$o_type}) /
            ($fam_fxncounts->{all_fam}->{$t_fxn} / 
             $fam_fxncounts->{totals_all_fam}) ;
         }
      }
      print FXNPROP join("\t", $t_fxn,
                           $fam_fxncounts->{under}->{$t_fxn},
                           $fam_fxncounts->{totals_under},
                           $fam_fxncounts->{over}->{$t_fxn},
                           $fam_fxncounts->{totals_over},
                           $fam_fxncounts->{all_fam}->{$t_fxn},
                           $fam_fxncounts->{totals_all_fam},
                           $fxn_propensities->{under}->{$t_fxn},
                           $fxn_propensities->{over}->{$t_fxn})."\n";
   }
   close(FXNPROP) ;

   my $fxnbroad_propensities ;
   open(FXNPROP, ">fxn_propensities_broad.sigunder_sigover.txt") ;
   foreach my $t_fxn ( sort keys %{$fam_fxnbroadcounts->{all_fam}}) {
#HERENOW090225_1953 
      foreach my $o_type (qw/under over/) {
         if (!exists $fam_fxnbroadcounts->{$o_type}->{$t_fxn}) {
            $fxnbroad_propensities->{$o_type}->{$t_fxn} = 0 ;
            $fam_fxnbroadcounts->{$o_type}->{$t_fxn} = 0 ;
         } else {
            $fxnbroad_propensities->{$o_type}->{$t_fxn} =
            ($fam_fxnbroadcounts->{$o_type}->{$t_fxn} / 
             $fam_fxnbroadcounts->{"totals_".$o_type}) /
            ($fam_fxnbroadcounts->{all_fam}->{$t_fxn} / 
             $fam_fxnbroadcounts->{totals_all_fam}) ;
         }
      }
      print FXNPROP join("\t", $t_fxn,
                           $fam_fxnbroadcounts->{under}->{$t_fxn},
                           $fam_fxnbroadcounts->{totals_under},
                           $fam_fxnbroadcounts->{over}->{$t_fxn},
                           $fam_fxnbroadcounts->{totals_over},
                           $fam_fxnbroadcounts->{all_fam}->{$t_fxn},
                           $fam_fxnbroadcounts->{totals_all_fam},
                           $fxnbroad_propensities->{under}->{$t_fxn},
                           $fxnbroad_propensities->{over}->{$t_fxn})."\n";
   }
   close(FXNPROP) ;

   foreach my $type (qw/pval0.01_under_overlap pval0.01_over_overlap/) {
      $results->{collate_perfam}->{num_families}->{$type} =
         keys %{$results->{collate_perfam}->{families}->{$type}} ; }


   foreach my $result_type (sort
      keys %{$results->{collate_perfam}->{num_families}}) {
      print "Families that $result_type: ".
         $results->{collate_perfam}->{num_families}->{$result_type}."\n" ; }


#


# Fig 3 - assign the alignment positions to non-overlapping sets L, P, L_and_P
# Fig 3a. bifunctional residue composition
# Fig 3b. bifunctional residue conservation

#   open(OUTCONSSCORE, "outconsscore.11773.out") ;
   open(OUTCONSSCORE, $input_fn->{consscore}) ;
   my $lastclass;
   my $local_posinfo;
   while (my $line = <OUTCONSSCORE>) {
      chomp $line;
      my ($class, $classtype, $btype, $aln_pos, $shannone, $numrestypes,
          $aa_freq) = split(/\t/, $line) ;
      $local_posinfo->{$aln_pos}->{type}->{$btype}++ ;
      $local_posinfo->{$aln_pos}->{shannone} = $shannone;
      $local_posinfo->{$aln_pos}->{numrestypes} = $numrestypes ;
      $local_posinfo->{$aln_pos}->{aa_freq} = $aa_freq;

      if (defined $lastclass &&
          $class ne $lastclass) {
         # do proper binning into non-overlapping classes
         foreach my $tpos (sort {$a <=> $b} keys %{$local_posinfo}) {
            my $cur_category;
            if ( exists $local_posinfo->{$tpos}->{type}->{E}) {
               push @{$cur_category}, "E" ;
            }

            if ( exists $local_posinfo->{$tpos}->{type}->{L_and_P}) {
               $histograms->{num_alnpos_type}->{LP}++ ;
               push @{$cur_category}, "LP" ;
            } elsif ( exists $local_posinfo->{$tpos}->{type}->{P} &&
                     !exists $local_posinfo->{$tpos}->{type}->{L}) {
               $histograms->{num_alnpos_type}->{P}++ ;
               push @{$cur_category}, "P" ;
            } elsif ( exists $local_posinfo->{$tpos}->{type}->{L} &&
                     !exists $local_posinfo->{$tpos}->{type}->{P}) {
               $histograms->{num_alnpos_type}->{L}++ ;
               push @{$cur_category}, "L" ;
            }

            my @aa_freq = ($local_posinfo->{$tpos}->{aa_freq} =~
                           /([A-Z\-])([0-9]+)/g) ;
            my $j = 0 ;
            while ($j <= ($#aa_freq - 1)) {
#               if ($aa_freq[$j] eq '-') {$j+=2; next;}
               if (!exists $standardres_20aa->{$aa_freq[$j]}) {$j += 2; next;}
               foreach my $t_category (@{$cur_category}) {
                  push @{$histograms->{shannone}->{$t_category}},
                     $local_posinfo->{$tpos}->{shannone} ;
                  push @{$histograms->{numrestypes}->{$t_category}},
                     $local_posinfo->{$tpos}->{numrestypes} ;
                  $histograms->{aa_freq}->{$t_category}->{$aa_freq[$j]} +=
                     $aa_freq[($j + 1)] ;
                  $histograms->{aa_freq}->{$t_category}->{total} +=
                     $aa_freq[($j + 1)] ;
               }
               $j += 2 ;
            }
         }
         $local_posinfo = {} ;
      }

      $lastclass = $class ;
   }
   close(OUTCONSSCORE);

   foreach my $btype (keys %{$histograms->{aa_freq}}) {
      open(AAHISTO, ">aahisto_nongap.$btype.txt") ;
      foreach my $aatype (sort keys %{$histograms->{aa_freq}->{$btype}}) {
         if ($aatype eq 'total') {next;}
         my $p_aa_given_btype = $histograms->{aa_freq}->{$btype}->{$aatype} /
             $histograms->{aa_freq}->{$btype}->{total} ;
         my $prop_aa_btype = $p_aa_given_btype / 
             ($histograms->{aa_freq}->{'E'}->{$aatype} /
             $histograms->{aa_freq}->{'E'}->{total}) ;
         my @outvals = ($aatype,
                        sprintf("%.4f",$p_aa_given_btype),
                        sprintf("%.4f",$prop_aa_btype)) ;
         print AAHISTO join("\t", @outvals)."\n" ;
      }
      close(AAHISTO) ;

# herenow090225_1215 - make sure that non-standard residues not counted
# or that the shannon entropy calculation in ASTRAL.pm is only over
# the 20 standard amino acids
      open(CONSTYPE, ">category_conservation.$btype.txt") ;
      foreach my $j ( 0 .. $#{$histograms->{shannone}->{$btype}}) {
         print CONSTYPE join("\t", $histograms->{shannone}->{$btype}->[$j],
                             $histograms->{numrestypes}->{$btype}->[$j])."\n";
      }
      close(CONSTYPE) ;
   }

# should also serve as a check on the clustered interaction numbers from
# the assign_pi_cluster and assign_pepnuci_cluster files.

# the runs are independent for peptides and the pi assignments...
   open(PERINST_PEP, $input_fn->{pINT_LINT}) ;
# only read in the SUMMARY LINES
   while (my $line = <PERINST_PEP>) {
      chomp $line ;
      if ($line =~ /^\#pINT_LSUM/) {
         $line =~ s/\#pINT_LSUM\t// ;
         my @t = split(/\t/, $line) ;
         foreach my $j ( 0 .. $#t) {
            $field2i->{$t[$j]} = $j ; }
         next;
      } elsif ($line !~ /^pINT_LSUM/) {
         next;
      }
      $line =~ s/^pINT_LSUM\t// ;

      my @t = split(/\t/, $line) ;

      my $int_sig = $t[$field2i->{SID}]."\t".$t[$field2i->{CHAIN}] ;
      $results->{collate_perinst}->{p}->{all_interactions}->{$int_sig}++ ;

      if ($t[$field2i->{cumulative_numres_l_and_p}] >=
          0.2 * $t[$field2i->{numres_p}]) {
         $results->{collate_perinst}->{p}->{min_20perc_cum_coverage}->{$int_sig}++ ;
      }

      if ($t[$field2i->{max_l_and_p}] >=
          0.2 * $t[$field2i->{numres_p}]) {
         $results->{collate_perinst}->{p}->{min_20perc_max_coverage}->{$int_sig}++ ;
      }
   }
   $field2i = {} ;
   close(PERINST_PEP) ;

   foreach my $cov_type (keys %{$results->{collate_perinst}->{p}}) {
      my $cur_num = keys %{$results->{collate_perinst}->{p}->{$cov_type}} ;
      print "peptide interactions of type $cov_type: $cur_num\n" ;
   }

# TODO (future...090204_1133)- calc seqid given family, two sids, alnpos set

   open(PERINST_P, $input_fn->{PINT_LINT}) ;
# only read in the SUMMARY LINES
   while (my $line = <PERINST_P>) {
      chomp $line ;
#pINT_LSUM      BDP91371-0_SCOP.d2oslb2 BDP91371-2_CHAIN-H      224     b.1.1.2 5       1       1       1nc4:DOF:006
      if ($line !~ /PINT_LSUM/) {next;}
      $line =~ s/PINT_LSUM\t// ;
      if ($line =~ /^#/) {
         $line =~ s/^\#// ;
         my @t = split(/\t/, $line) ;
         foreach my $j ( 0 .. $#t) {
            $field2i->{$t[$j]} = $j ; }
         next;
      }
      my @t = split(/\t/, $line) ;

#      push @{$results->{collate_perinst}->{"p"}->{histogram}->{max_l_and_p}},
#         $t[$field2i->{max_l_and_p}] ;
#
#      push @{$results->{collate_perinst}->{"p"}->{histogram}->{cum_l_and_p}},
#         $t[$field2i->{cumulative_numres_l_and_p}] ;

      my $int_sig = $t[$field2i->{SID1}]."\t".$t[$field2i->{SID2}] ;

      my $ptype = 'Pinter' ;
      if ($t[$field2i->{CHAINS}] eq 'same') {
         $ptype = 'Pintra' ; }

      $results->{collate_perinst}->{$ptype}->{all_interactions}->{$int_sig}++ ;

      if ( ($t[$field2i->{cumulative_numres_l_and_p_1}] >=
            0.2 * $t[$field2i->{numres_p_1}]) ||
           ($t[$field2i->{cumulative_numres_l_and_p_2}] >=
            0.2 * $t[$field2i->{numres_p_2}])) {
         $results->{collate_perinst}->{$ptype}->{min_20perc_cum_coverage}->{$int_sig}++;
      }

      if ( ($t[$field2i->{max_l_and_p_1}] >=
            0.2 * $t[$field2i->{numres_p_1}]) ||
           ($t[$field2i->{max_l_and_p_2}] >=
            0.2 * $t[$field2i->{numres_p_2}])) {
         $results->{collate_perinst}->{$ptype}->{min_20perc_max_coverage}->{$int_sig}++;
      }
   }
   $field2i = {} ;
   close(PERINST_P) ;

   foreach my $ptype (sort keys %{$results->{collate_perinst}}) {
      foreach my $cov_type (sort keys %{$results->{collate_perinst}->{$ptype}}){
         my $cur_num=keys %{$results->{collate_perinst}->{$ptype}->{$cov_type}};
         print "$ptype interactions of type $cov_type: $cur_num\n" ;
      }
   }


# PER LIGAND ANALYSIS - get maximal d-d and d-p coverages for each ligand.
# dump info with liginfo->{MW} and liginfo->{numatoms_nonh}
# per ligand analysis not very useful... fuck that and make the full
# sequence id cutoff scan through to get intrface coverage as a fxn of seqid
# - more convincing argument

# HERENOW - fpd090305_1354  editing for 2D histogram of coverage vs seqid c/o
#  can't do that calculation here for cumulative coverage - go back to
#  collation routine and keep the processing logic there.
#  max/lig_giving_max/cum coverage at 0-100 q 5 cutoffs
#  - this expanded {pP}INT_LSUM will serve as data for the interface
#    details page to display on PIBASE.

   my ($hist2d_binsize, $hist2d_maxbin);
   $hist2d_binsize->{ligbs_seqid} = 5 ;
   $hist2d_binsize->{coverage} = 5 ;
   foreach my $bin_type (keys %{$hist2d_binsize}) {
      $hist2d_maxbin->{$bin_type} = POSIX::floor(99 /
         $hist2d_binsize->{$bin_type}) ;}

   my $perlig_maxpcoverage ;
   if (0) {
   open(PERINST_PEP, $input_fn->{pINT_LINT}) ;
#   $field2i->{SUM} ; # summary lines
#   $field2i->{ENTRY} ; # individual entry lines
   while (my $line = <PERINST_PEP>) {
      chomp $line;
      if ($line =~ /^pINT_LSUM/ || $line =~ /^#pINT_LSUM/) {next;}
      if ($line =~ /^#SID/) {
         $line =~ s/^#// ;
         my @t = split(/\t/, $line) ;
         foreach my $j ( 0 .. $#t) {
            $field2i->{$t[$j]} = $j ; }
         next;
      }
      if ($line =~ /^#/) {next;}

      my @t = split(/\t/, $line) ;
      my $curlig = $t[$field2i->{LIG}] ;
      my (undef, $curligcode, undef) = split(/\:/, $curlig) ;
      my $cur_pcover = $t[$field2i->{numres_l_and_p}] /
                       $t[$field2i->{numres_p}] ;
      if (!exists $perlig_maxpcoverage->{p}->{$curligcode} ||
          $cur_pcover > $perlig_maxpcoverage->{p}->{$curligcode}) {
         $perlig_maxpcoverage->{p}->{$curligcode} = $cur_pcover ;
      }
   }
   $field2i = {} ;
   close(PERINST_PEP) ;
   }


# herenow 090227_1549 - change file opens so that actual name only
#  specified one place, and keep the final collate results in the data
# directory rather than in the mishmash of the run directory.
   if (0) {
   open(PERINST_P, $input_fn->{PINT_LINT}) ;
   while (my $line = <PERINST_P>) {
      chomp $line ;
      if ($line =~ /^#PINT_LSUM/ ||
          $line =~ /^PINT_LSUM/) {
         next;
      } elsif ($line =~ /^#SID1/) {
         my @t = split(/\t/, $line) ;
         foreach my $j ( 0 .. $#t) {
            $field2i->{$t[$j]} = $j ; }
         next;
      }
      my @t = split(/\t/, $line) ;
      my $curlig = $t[$field2i->{LIG}] ;
      my (undef, $curligcode, undef) = split(/\:/, $curlig) ;
      my $cur_pcover = $t[$field2i->{numres_l_and_p_side}] /
                       $t[$field2i->{numres_p_side}] ;
      my $btype ;
      if ($t[$field2i->{CHAINS}] eq 'same') {
         $btype = "Pintra" ;
      } else {
         $btype = "Pinter" ;
      }

      if (!exists $perlig_maxpcoverage->{$btype}->{$curligcode} ||
          $cur_pcover > $perlig_maxpcoverage->{$btype}->{$curligcode}) {
         $perlig_maxpcoverage->{$btype}->{$curligcode} = $cur_pcover ;
      }
   }
   $field2i = {} ;
   close(PERINST_P) ;

   foreach my $int_type (keys %{$perlig_maxpcoverage}) {
      open(PERLIGMAXPCOVERAGE, ">perlig_maxpcoverage.$int_type.txt") ;
      foreach my $lig (sort keys %{$perlig_maxpcoverage->{$int_type}}) {
         print PERLIGMAXPCOVERAGE join("\t", $lig, $int_type,
             $perlig_maxpcoverage->{$int_type}->{$lig},
             $liginfo->{mw}->{$lig}, $liginfo->{numatoms_nonh}->{$lig}
             )."\n" ;
      }
      close(PERLIGMAXPCOVERAGE) ;
   }
   }


# reparse - keep pint_LSUM lines only and bin up the per-seqid coverage #s
   my ($out_fh, $out_fn) ;
   my $seqid_cutoffs = {30 => 1,50 => 1,90 => 1} ;
   foreach my $int_type (qw/Pinter Pintra p/) {
      foreach my $cover_type (qw/cum max/) {
         $out_fn->{$int_type."_".$cover_type} = $cover_type.
               "_coverage_vs_seqid.$int_type".".points" ;
         foreach my $seqid (keys %{$seqid_cutoffs}) {
            $out_fn->{$int_type."_".$cover_type."_".$seqid} =
               $cover_type."_coverage_vs_seqid.$int_type.$seqid.points" ; }}}

   foreach my $fn_type (keys %{$out_fn}) {
      open($out_fh->{$fn_type}, ">".$out_fn->{$fn_type}) ; }

   my $hist2d_ligbsseqid_cover ; #->{p|Pinter|Pintra}->{seqid}

# per seqid analysis
   open(PERINST_PEP, $input_fn->{pINT_LINT}) ;
#   $field2i->{SUM} ; # summary lines
#   $field2i->{ENTRY} ; # individual entry lines
   while (my $line = <PERINST_PEP>) {
      chomp $line;
      if ($line !~ /pINT_LSUM/) {next;}
      if ($line =~ /^#pINT_LSUM/) {
          $line =~ s/^#pINT_LSUM\t// ;
         my @t = split(/\t/, $line) ;
         foreach my $j ( 0 .. $#t) {
            $field2i->{$t[$j]} = $j ; }
         next;
      }
      $line =~ s/^pINT_LSUM\t// ;
      my @t = split(/\t/, $line) ;

      foreach my $j ( 0 .. $hist2d_maxbin->{ligbs_seqid}) {
         my $t_seqid_cutoff = $j * $hist2d_binsize->{ligbs_seqid} ;
         foreach my $t_cover (qw/max cum/) {
            my $cur_cover = 0;
            if ($t[$field2i->{numres_p}] > 0) {
               $cur_cover = 
         ($t[$field2i->{$t_cover."_l_and_p_perseqid_".$t_seqid_cutoff}] /
          $t[$field2i->{numres_p}]) ;
            }
            print {$out_fh->{"p_".$t_cover}} join ("\t",
               $t_seqid_cutoff, $cur_cover)."\n";

            my $cur_cover_bin = POSIX::floor(100 * $cur_cover /
                                             $hist2d_binsize->{coverage}) ;

            $hist2d_ligbsseqid_cover->{$t_cover}->{p}->[$j]->[$cur_cover_bin]++;

            if (exists $seqid_cutoffs->{$t_seqid_cutoff}) {
               print {$out_fh->{"p_$t_cover"."_".$t_seqid_cutoff}}
                  $cur_cover."\n";}
         }

      }
   }
   close($out_fh->{"p_max"}) ; close($out_fh->{"p_cum"}) ;
   $field2i = {} ;
   close(PERINST_PEP) ;


   open(PERINST_P, $input_fn->{PINT_LINT}) ;
   while (my $line = <PERINST_P>) {
      chomp $line ;
      if ($line !~ /PINT_LSUM/) {next;}
      if ($line =~ /^#PINT_LSUM/) {
         $line =~ s/^#PINT_LSUM\t// ;
         my @t = split(/\t/, $line) ;
         foreach my $j ( 0 .. $#t) {
            $field2i->{$t[$j]} = $j ; }
         next;
      }
      $line =~ s/^PINT_LSUM\t// ;
      my @t = split(/\t/, $line) ;
      my $btype ;
      if ($t[$field2i->{CHAINS}] eq 'same') {
         $btype = "Pintra" ;
      } else {
         $btype = "Pinter" ;
      }

      foreach my $j ( 0 .. $hist2d_maxbin->{ligbs_seqid}) {
         my $t_seqid_cutoff = $j * $hist2d_binsize->{ligbs_seqid} ;
         foreach my $t_cover (qw/max cum/) {
            my $cur_cover_max = 0;
            foreach my $s (1, 2) {
               if ($t[$field2i->{"numres_p_$s"}] > 0) {
                  my $cur_cover_side = 
         ($t[$field2i->{$t_cover."_l_and_p_".$s."_perseqid_".$t_seqid_cutoff}]/
          $t[$field2i->{"numres_p_$s"}]) ;
                  if ($cur_cover_side > $cur_cover_max) {
                     $cur_cover_max = $cur_cover_side ; }
               }
            }
            print {$out_fh->{$btype."_".$t_cover}} join ("\t",
               $t_seqid_cutoff, $cur_cover_max)."\n";

            my $cur_cover_bin = POSIX::floor(100 * $cur_cover_max /
                                             $hist2d_binsize->{coverage}) ;

            $hist2d_ligbsseqid_cover->{$t_cover}->{$btype}->[$j]->[$cur_cover_bin]++ ;

            if (exists $seqid_cutoffs->{$t_seqid_cutoff}) {
               print {$out_fh->{$btype."_".$t_cover."_".$t_seqid_cutoff}}
                  $cur_cover_max."\n";}
         }
      }
   }
   close(OUTF) ;
   $field2i = {} ;
   close(PERINST_P) ;
   close($out_fh->{"Pinter_max"}) ; close($out_fh->{"Pinter_cum"}) ;
   close($out_fh->{"Pintra_max"}) ; close($out_fh->{"Pintra_cum"}) ;

   foreach my $cover_type (qw/max cum/) {
      foreach my $int_type (keys %{$hist2d_ligbsseqid_cover->{$cover_type}}) {
         open(OUTF, ">$cover_type"."_coverage_vs_seqid.$int_type".".out_gnuplot") ;
         foreach my $i (0 .. $hist2d_maxbin->{coverage}) {
            foreach my $j (0 .. $hist2d_maxbin->{ligbs_seqid}) {
               my @outvals = ($j * $hist2d_binsize->{ligbs_seqid},
                              $i * $hist2d_binsize->{coverage}) ;
               if (exists $hist2d_ligbsseqid_cover->{$cover_type}->{$int_type}->[$j] &&
     defined $hist2d_ligbsseqid_cover->{$cover_type}->{$int_type}->[$j]->[$i]) {
                  push @outvals,
     $hist2d_ligbsseqid_cover->{$cover_type}->{$int_type}->[$j]->[$i] ;
               } else {
                  push @outvals, 0 ;
               }
               print OUTF join("\t", @outvals)."\n" ;
            }
            print OUTF "\n" ;
         }
         close(OUTF) ;
      }
   }

}


sub BK_PRE20090305_make_paper_figures {

   require Bit::Vector;

#   require R;
#   require RReferences ;
#   R::initR("--silent") ;

   my $pilig_specs = set_pilig_specs() ;
   my $pibase_specs = pibase::get_specs() ;
   my $results ;

   my $superfam_fxn = readin_superfam_fxn() ;

# 1. number of families/domains/interactions (in range (right?) - ie SCOP classes a-g; ie number of families with ASTEROIDS alignments)
# number of alignments in /groups/eddy/home/davisf/work/pibase/pibase200708/external_data/asteroids/1.73/alignments/fam

   my $dir_astral_fam_aln = $pibase_specs->{asteroids}->{'fam_aln'};
   my @all_fams = glob($dir_astral_fam_aln.'/*aln') ;
   my $num_fam = $#all_fams + 1 ;
   $results->{num_families}->{total} = $num_fam ;

   my $pb = _pilig_tod_pibase_preload() ;

# 2. number of families/domains/interactions that bind in range small molecules

   my $liginfo = _pilig_load_liginfo() ;
   my $standardres =  _list_standardres();
   my $standardres_20aa =  _list_standardres_20aa();

   my $class2alnlength ;
   my $ligbits = readin_ligassignments({
      fn => $pilig_specs->{outfiles}->{assign_lig},
      liginfo => $liginfo,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   });
   {
      foreach my $sid (keys %{$ligbits->{fam}}) {
         my $cur_class = $pb->{sid2class}->{fam}->{$sid} ;
         $results->{families}->{bind_ligands}->{$cur_class}++ ;
      }
      $results->{num_families}->{bind_ligands} = 
         keys %{$results->{families}->{bind_ligands}} ;
   }


# 3. number of families/domains/interactions that bind peptides
   my $pepnucibits = readin_pepnuciassignments({
      fn => $pilig_specs->{outfiles}->{assign_pepnuci},
      pb => $pb,
      pilig_specs => $pilig_specs,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   });
   {
      foreach my $sid (keys %{$pepnucibits->{pep}->{fam}}) {
         my $cur_class = $pb->{sid2class}->{fam}->{$sid} ;
         $results->{families}->{bind_pep}->{$cur_class}++ ;
         $results->{families}->{bind_protein}->{$cur_class}++ ;
      }
      $results->{num_families}->{bind_pep} =
         keys %{$results->{families}->{bind_pep}} ;
   }


# 4. number of families/domains/interactions that bind domains (inter)
# 5. number of families/domains/interactions that bind domains (intra)


   my $piassignments= readin_piassignments({
      fn => $pilig_specs->{outfiles}->{assign_pi},
      pb => $pb,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   });
   my $pibits = $piassignments->{pibits} ;
   {
      foreach my $sid (keys %{$pibits->{fam}}) {
         my $cur_class = $pb->{sid2class}->{fam}->{$sid} ;
         foreach my $btype (keys %{$pibits->{fam}->{$sid}}) {
            $results->{families}->{"bind_$btype"}->{$cur_class}++ ;
            $results->{families}->{bind_protein}->{$cur_class}++ ;
         }
      }
      $results->{num_families}->{bind_Pinter} =
         keys %{$results->{families}->{bind_Pinter}} ;

      $results->{num_families}->{bind_Pintra} =
         keys %{$results->{families}->{bind_Pintra}} ;
   }

   my $pi_clus_assignments= readin_piassignments({
      fn => $pilig_specs->{outfiles}->{assign_pi_clusters},
      clustrep_fl => 1,
      pb => $pb,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   });
   foreach my $pi_type (qw/all clustered/) {
      my $cur_interfaces ;
      if ($pi_type eq 'all') {
         $cur_interfaces = $piassignments->{interfaces} ;
      } else {
         $cur_interfaces = $pi_clus_assignments->{interfaces} ;
      }

      foreach my $sid12 (keys %{$cur_interfaces}) {
         my $chains = $cur_interfaces->{$sid12}->{chains} ;
         my $cur_type = 'Pinter';
         if ($chains eq 'same') { $cur_type = 'Pintra'; }
         $results->{interactions}->{$pi_type}->{$cur_type}->{$sid12}++ ;
      }
   }

   foreach my $pi_type ( sort keys %{$results->{interactions}}) {
      foreach my $int_type (sort keys %{$results->{interactions}->{$pi_type}}) {
         $results->{num_interactions}->{$pi_type}->{$int_type} = 
            keys %{$results->{interactions}->{$pi_type}->{$int_type}} ;
         print "Number of $pi_type interactions that $int_type: ".
            $results->{num_interactions}->{$pi_type}->{$int_type}."\n" ; }}

   my $pepnuci_clus_bits = readin_pepnuciassignments({
      fn => $pilig_specs->{outfiles}->{assign_pepnuci_clusters},
      pb => $pb,
      clustrep_fl => 1,
      pilig_specs => $pilig_specs,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   });
   foreach my $pi_type (qw/all clustered/) {
      my $curbits ;
      if ($pi_type eq 'all') {
         $curbits = $pepnucibits->{pep}->{fam} ;
      } else {
         $curbits = $pepnuci_clus_bits->{pep}->{fam} ;
      }
      
      foreach my $sid (keys %{$curbits}) {
         foreach my $chain (keys %{$curbits->{$sid}}) {
            if ($chain eq 'cumulative') {next;}
            $results->{interactions}->{$pi_type}->{p}->{$sid."\t".$chain}++ ;}}

      $results->{num_interactions}->{$pi_type}->{p} =
         keys %{$results->{interactions}->{$pi_type}->{p}} ;
   }

   print "\n\n\n" ;
   foreach my $pi_type ( sort keys %{$results->{interactions}}) {
      foreach my $int_type (sort keys %{$results->{interactions}->{$pi_type}}) {
         $results->{num_interactions}->{$pi_type}->{$int_type} = 
            keys %{$results->{interactions}->{$pi_type}->{$int_type}} ;
         print "Number of $pi_type interactions that $int_type: ".
            $results->{num_interactions}->{$pi_type}->{$int_type}."\n" ; }}


# 3b. number of families/domains/interactions that bind peptides and small molecules
# 4b. number of families/domains/interactions that bind domains (inter) and small molecules
# 5b. number of families/domains/interactions that bind domains (intra) and small molecules
   $results->{num_families}->{bind_protein} =
      keys %{$results->{families}->{bind_protein}} ;

   foreach my $query_type (qw/bind_protein bind_pep bind_Pinter bind_Pintra/) {
      foreach my $ligand_fam (keys %{$results->{families}->{bind_ligands}}){
         if (exists $results->{families}->{$query_type}->{$ligand_fam}) {
            $results->{families_overlap}->{$query_type."_and_bind_ligands"}->{$ligand_fam}++ ; }
         $results->{num_families}->{$query_type."_and_bind_ligands"} =
keys %{$results->{families_overlap}->{$query_type."_and_bind_ligands"}}
      }
   }

   foreach my $fam_type ( sort keys %{$results->{num_families}}) {
      print "Number of families that ".$fam_type.": ".
         $results->{num_families}->{$fam_type}."\n" ; }

   open(PERFAM, "collate_perfam_cluster.20090204.out") ;
   my $field2i ;
   my $histograms ;
   my $fam2pval;
   my $fam_fxncounts ;
   my $fam_fxnbroadcounts ;
   while (my $line = <PERFAM>) {
      chomp $line ;

      if ($line =~ /^\#/) {
         if ($line =~ /^#[0-9]/) { next; }
         $line =~ s/\#// ;
         my @t = split(/\t/, $line) ;
         foreach my $j ( 0 .. $#t) {
            $field2i->{$t[$j]} = $j ;
         }
         next;
      }

      my @t = split(/\t/, $line) ;

      my $curfam = $t[$field2i->{class}] ;

      $results->{collate_perfam}->{families}->{total}->{$curfam}++ ;

      if ($line =~ /NOINTERACTIONS/) {
         $results->{collate_perfam}->{families}->{"No interactions"}->{$curfam}++ ;
         next ;
      }

      if ($t[$field2i->{norm_L}] > 0) {
         $results->{collate_perfam}->{num_families}->{"Bind ligands 250-1000MW"}++ ;
         $results->{collate_perfam}->{families}->{"Bind ligands 250-1000MW"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_L}] >= 5) {
         $results->{collate_perfam}->{num_families}->{"Bind ligands 250-1000MW, >=5 positions"}++ ;
         $results->{collate_perfam}->{families}->{"Bind ligands 250-1000MW, >=5 positions"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_P}] > 0) {
         $results->{collate_perfam}->{num_families}->{"Bind proteins"}++ ;
         $results->{collate_perfam}->{families}->{"Bind proteins"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_Pintra}] > 0) {
         $results->{collate_perfam}->{num_families}->{"Bind domains(intra)"}++ ;
         $results->{collate_perfam}->{families}->{"Bind domains (intra)"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_Pinter}] > 0) {
         $results->{collate_perfam}->{num_families}->{"Bind domains(inter)"}++ ;
         $results->{collate_perfam}->{families}->{"Bind domains (inter)"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_p}] > 0) {
         $results->{collate_perfam}->{num_families}->{"Bind peptides"}++ ;
         $results->{collate_perfam}->{families}->{"Bind peptides"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_L}] > 0 &&
          $t[$field2i->{norm_P}] > 0) {

         my ($cursf) = ($curfam =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
         if (!exists $superfam_fxn->{sf2fxn_broad}->{$cursf}) {
            print STDERR "WARNING: no function found for $curfam ($cursf)\n" ;
         } else {
            $fam_fxnbroadcounts->{all_fam}->{$superfam_fxn->{sf2fxn_broad}->{$cursf}}++;
            $fam_fxnbroadcounts->{totals_all_fam}++ ;

            $fam_fxncounts->{all_fam}->{$superfam_fxn->{sf2fxn}->{$cursf}}++;
            $fam_fxncounts->{totals_all_fam}++ ;
         }


         $fam2pval->{over}->{$t[$field2i->{class}]} = $t[$field2i->{pval_L_P}];
         $fam2pval->{under}->{$t[$field2i->{class}]} =
            $t[$field2i->{pval_L_P_less}];
         $results->{collate_perfam}->{num_families}->{"Bind proteins and ligands"}++ ;
         $results->{collate_perfam}->{families}->{"Bind proteins and ligands"}++ ;
      }

      if ($t[$field2i->{norm_L_and_P}] >= 5) {
         $results->{collate_perfam}->{num_families}->{"Bind proteins and ligands, >=5 multifunc positions"}++ ;
         $results->{collate_perfam}->{families}->{"Bind proteins and ligands, >=5 multifunc positions"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_L_and_Pinter}] >= 5) {
         $results->{collate_perfam}->{num_families}->{"Bind domains(inter) and ligands, >=5 multifunc positions"}++ ;
         $results->{collate_perfam}->{families}->{"Bind domains(inter) and ligands, >=5 multifunc positions"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_L_and_Pintra}] >= 5) {
         $results->{collate_perfam}->{num_families}->{"Bind domains(intra) and ligands, >=5 multifunc positions"}++ ;
         $results->{collate_perfam}->{families}->{"Bind domains(intra) and ligands, >=5 multifunc positions"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_L_and_p}] >= 5) {
         $results->{collate_perfam}->{num_families}->{"Bind peptides and ligands, >=5 multifunc positions"}++ ;
         $results->{collate_perfam}->{families}->{"Bind peptides and ligands, >=5 multifunc positions"}->{$curfam}++ ;
      }

      if ($t[$field2i->{norm_L}] > 0 &&
          $t[$field2i->{norm_Pintra}] > 0) {
         $results->{collate_perfam}->{families}->{"Bind domains(intra) and ligands"}->{$curfam}++;
         $results->{collate_perfam}->{num_families}->{"Bind domains(intra) and ligands"}++;
      }

      if ($t[$field2i->{norm_L}] > 0 &&
          $t[$field2i->{norm_Pinter}] > 0) {
         $results->{collate_perfam}->{families}->{"Bind proteins(inter) and ligands"}->{$curfam}++;
         $results->{collate_perfam}->{num_families}->{"Bind proteins(inter) and ligands"}++;
      }

      if ($t[$field2i->{norm_L}] > 0 &&
          $t[$field2i->{norm_p}] > 0) {
         $results->{collate_perfam}->{families}->{"Bind peptides and Ligands"}->{$curfam}++;
         $results->{collate_perfam}->{num_families}->{"Bind peptides and ligands"}++;
      }

      push @{$histograms->{num_alnpos_E}}, 
         $t[$field2i->{norm_E}] ;

#      push @{$histograms->{num_alnpos_E}}, 
#         $t[$field2i->{norm_E}] ;
#
#      push @{$histograms->{num_alnpos_Lmwinrange}},
#         $t[$field2i->{norm_Lmwinrange}];
#
#      push @{$histograms->{num_alnpos_Pinter}},
#         $t[$field2i->{norm_Pinter}];
#
#      push @{$histograms->{num_alnpos_Pintra}},
#         $t[$field2i->{norm_Pintra}];
#
#      push @{$histograms->{num_alnpos_P}},
#         $t[$field2i->{norm_P}];
#
#      push @{$histograms->{num_alnpos_Lmwinrange_Pinter}},
#         $t[$field2i->{norm_Lmwinrange_and_Pinter}];
#
#      push @{$histograms->{num_alnpos_Lmwinrange_Pintra}},
#         $t[$field2i->{norm_Lmwinrange_and_Pintra}];
#
#      push @{$histograms->{num_alnpos_Lmwinrange_P}},
#         $t[$field2i->{norm_Lmwinrange_and_P}];
#
#      push @{$histograms->{cons_shannone_L}},
#         $t[$field2i->{norm_Lmwinrange_and_P}];
   }
   $field2i = {} ;
   close(PERFAM) ;

   foreach my $fam (keys %{$fam2pval->{under}}) {
      my ($sf) = ($fam =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
      my $cur_fxn ; my $cur_fxn_broad ;
      if (exists $superfam_fxn->{sf2fxn_broad}->{$sf}) {
         $cur_fxn = $superfam_fxn->{sf2fxn}->{$sf} ;
         $cur_fxn_broad = $superfam_fxn->{sf2fxn_broad}->{$sf} ;
      }

      if ($fam2pval->{under}->{$fam} < 0.01) {
         $results->{collate_perfam}->{num_families}->{"pval0.01_under_overlap"}->{$fam}++ ;
         if (defined $cur_fxn) {
            $fam_fxnbroadcounts->{under}->{$cur_fxn_broad}++ ;
            $fam_fxnbroadcounts->{totals_under}++ ;

            $fam_fxncounts->{under}->{$cur_fxn}++ ;
            $fam_fxncounts->{totals_under}++ ;
         }
      }

      if ($fam2pval->{over}->{$fam} < 0.01) {
         $results->{collate_perfam}->{num_families}->{"pval0.01_over_overlap"}->{$fam}++ ;
         if (defined $cur_fxn) {
            $fam_fxnbroadcounts->{over}->{$cur_fxn_broad}++ ;
            $fam_fxnbroadcounts->{totals_over}++ ;

            $fam_fxncounts->{over}->{$cur_fxn}++ ;
            $fam_fxncounts->{totals_over}++ ;
         }
      }
   }

   my $fxn_propensities ;
   open(FXNPROP, ">fxn_propensities_detailed.sigunder_sigover.txt") ;
   foreach my $t_fxn ( sort keys %{$fam_fxncounts->{all_fam}}) {
#HERENOW090225_1953 
      foreach my $o_type (qw/under over/) {
         if (!exists $fam_fxncounts->{$o_type}->{$t_fxn}) {
            $fam_fxncounts->{$o_type}->{$t_fxn} = 0 ;
            $fxn_propensities->{$o_type}->{$t_fxn} = 0 ;
         } else {
            $fxn_propensities->{$o_type}->{$t_fxn} =
            ($fam_fxncounts->{$o_type}->{$t_fxn} / 
             $fam_fxncounts->{"totals_".$o_type}) /
            ($fam_fxncounts->{all_fam}->{$t_fxn} / 
             $fam_fxncounts->{totals_all_fam}) ;
         }
      }
      print FXNPROP join("\t", $t_fxn,
                           $fam_fxncounts->{under}->{$t_fxn},
                           $fam_fxncounts->{totals_under},
                           $fam_fxncounts->{over}->{$t_fxn},
                           $fam_fxncounts->{totals_over},
                           $fam_fxncounts->{all_fam}->{$t_fxn},
                           $fam_fxncounts->{totals_all_fam},
                           $fxn_propensities->{under}->{$t_fxn},
                           $fxn_propensities->{over}->{$t_fxn})."\n";
   }
   close(FXNPROP) ;

   my $fxnbroad_propensities ;
   open(FXNPROP, ">fxn_propensities_broad.sigunder_sigover.txt") ;
   foreach my $t_fxn ( sort keys %{$fam_fxnbroadcounts->{all_fam}}) {
#HERENOW090225_1953 
      foreach my $o_type (qw/under over/) {
         if (!exists $fam_fxnbroadcounts->{$o_type}->{$t_fxn}) {
            $fxnbroad_propensities->{$o_type}->{$t_fxn} = 0 ;
            $fam_fxnbroadcounts->{$o_type}->{$t_fxn} = 0 ;
         } else {
            $fxnbroad_propensities->{$o_type}->{$t_fxn} =
            ($fam_fxnbroadcounts->{$o_type}->{$t_fxn} / 
             $fam_fxnbroadcounts->{"totals_".$o_type}) /
            ($fam_fxnbroadcounts->{all_fam}->{$t_fxn} / 
             $fam_fxnbroadcounts->{totals_all_fam}) ;
         }
      }
      print FXNPROP join("\t", $t_fxn,
                           $fam_fxnbroadcounts->{under}->{$t_fxn},
                           $fam_fxnbroadcounts->{totals_under},
                           $fam_fxnbroadcounts->{over}->{$t_fxn},
                           $fam_fxnbroadcounts->{totals_over},
                           $fam_fxnbroadcounts->{all_fam}->{$t_fxn},
                           $fam_fxnbroadcounts->{totals_all_fam},
                           $fxnbroad_propensities->{under}->{$t_fxn},
                           $fxnbroad_propensities->{over}->{$t_fxn})."\n";
   }
   close(FXNPROP) ;


   foreach my $result_type (sort
      keys %{$results->{collate_perfam}->{num_families}}) {
      print "Families that $result_type: ".
         $results->{collate_perfam}->{num_families}->{$result_type}."\n" ;
   }


#


# Fig 3 - assign the alignment positions to non-overlapping sets L, P, L_and_P
# Fig 3a. bifunctional residue composition
# Fig 3b. bifunctional residue conservation

   open(OUTCONSSCORE, "outconsscore.11773.out") ;
   my $lastclass;
   my $local_posinfo;
   while (my $line = <OUTCONSSCORE>) {
      chomp $line;
      my ($class, $classtype, $btype, $aln_pos, $shannone, $numrestypes,
          $aa_freq) = split(/\t/, $line) ;
      $local_posinfo->{$aln_pos}->{type}->{$btype}++ ;
      $local_posinfo->{$aln_pos}->{shannone} = $shannone;
      $local_posinfo->{$aln_pos}->{numrestypes} = $numrestypes ;
      $local_posinfo->{$aln_pos}->{aa_freq} = $aa_freq;

      if (defined $lastclass &&
          $class ne $lastclass) {
         # do proper binning into non-overlapping classes
         foreach my $tpos (sort {$a <=> $b} keys %{$local_posinfo}) {
            my $cur_category;
            if ( exists $local_posinfo->{$tpos}->{type}->{E}) {
               push @{$cur_category}, "E" ;
            }

            if ( exists $local_posinfo->{$tpos}->{type}->{L_and_P}) {
               $histograms->{num_alnpos_type}->{LP}++ ;
               push @{$cur_category}, "LP" ;
            } elsif ( exists $local_posinfo->{$tpos}->{type}->{P} &&
                     !exists $local_posinfo->{$tpos}->{type}->{L}) {
               $histograms->{num_alnpos_type}->{P}++ ;
               push @{$cur_category}, "P" ;
            } elsif ( exists $local_posinfo->{$tpos}->{type}->{L} &&
                     !exists $local_posinfo->{$tpos}->{type}->{P}) {
               $histograms->{num_alnpos_type}->{L}++ ;
               push @{$cur_category}, "L" ;
            }

            my @aa_freq = ($local_posinfo->{$tpos}->{aa_freq} =~
                           /([A-Z\-])([0-9]+)/g) ;
            my $j = 0 ;
            while ($j <= ($#aa_freq - 1)) {
#               if ($aa_freq[$j] eq '-') {$j+=2; next;}
               if (!exists $standardres_20aa->{$aa_freq[$j]}) {$j += 2; next;}
               foreach my $t_category (@{$cur_category}) {
                  push @{$histograms->{shannone}->{$t_category}},
                     $local_posinfo->{$tpos}->{shannone} ;
                  push @{$histograms->{numrestypes}->{$t_category}},
                     $local_posinfo->{$tpos}->{numrestypes} ;
                  $histograms->{aa_freq}->{$t_category}->{$aa_freq[$j]} +=
                     $aa_freq[($j + 1)] ;
                  $histograms->{aa_freq}->{$t_category}->{total} +=
                     $aa_freq[($j + 1)] ;
               }
               $j += 2 ;
            }
         }
         $local_posinfo = {} ;
      }

      $lastclass = $class ;
   }
   close(OUTCONSSCORE);

   foreach my $btype (keys %{$histograms->{aa_freq}}) {
      open(AAHISTO, ">aahisto_nongap.$btype.txt") ;
      foreach my $aatype (sort keys %{$histograms->{aa_freq}->{$btype}}) {
         if ($aatype eq 'total') {next;}
         my $p_aa_given_btype = $histograms->{aa_freq}->{$btype}->{$aatype} /
             $histograms->{aa_freq}->{$btype}->{total} ;
         my $prop_aa_btype = $p_aa_given_btype / 
             ($histograms->{aa_freq}->{'E'}->{$aatype} /
             $histograms->{aa_freq}->{'E'}->{total}) ;
         my @outvals = ($aatype,
                        sprintf("%.4f",$p_aa_given_btype),
                        sprintf("%.4f",$prop_aa_btype)) ;
         print AAHISTO join("\t", @outvals)."\n" ;
      }
      close(AAHISTO) ;

# herenow090225_1215 - make sure that non-standard residues not counted
# or that the shannon entropy calculation in ASTRAL.pm is only over
# the 20 standard amino acids
      open(CONSTYPE, ">category_conservation.$btype.txt") ;
      foreach my $j ( 0 .. $#{$histograms->{shannone}->{$btype}}) {
         print CONSTYPE join("\t", $histograms->{shannone}->{$btype}->[$j],
                             $histograms->{numrestypes}->{$btype}->[$j])."\n";
      }
      close(CONSTYPE) ;
   }

# should also serve as a check on the clustered interaction numbers from
# the assign_pi_cluster and assign_pepnuci_cluster files.

# the runs are independent for peptides and the pi assignments...
   open(PERINST_PEP, "pINT_collate_perinstance.20090128.out") ;
# only read in the SUMMARY LINES
   while (my $line = <PERINST_PEP>) {
      chomp $line ;
#pINT_LSUM      BDP91371-0_SCOP.d2oslb2 BDP91371-2_CHAIN-H      224     b.1.1.2 5       1       1       1nc4:DOF:006
      if ($line !~ /^\#pINT_LSUM/) {next;}
      $line =~ s/\#pINT_LSUM\t// ;
      my @t = split(/\t/, $line) ;

      if ($line =~ /^SID/) {
         foreach my $j ( 0 .. $#t) {
            $field2i->{$t[$j]} = $j ; }
         next;
      }

#      push @{$results->{collate_perinst}->{"p"}->{histogram}->{max_l_and_p}},
#         $t[$field2i->{max_l_and_p}] ;
#
#      push @{$results->{collate_perinst}->{"p"}->{histogram}->{cum_l_and_p}},
#         $t[$field2i->{cumulative_numres_l_and_p}] ;

      my $int_sig = $t[$field2i->{SID}]."\t".$t[$field2i->{CHAIN}] ;
      $results->{collate_perinst}->{p}->{all_interactions}->{$int_sig}++ ;

      if ($t[$field2i->{cumulative_numres_l_and_p}] >=
          0.2 * $t[$field2i->{numres_p}]) {
         $results->{collate_perinst}->{p}->{min_20perc_cum_coverage}->{$int_sig}++ ;
      }

      if ($t[$field2i->{max_l_and_p}] >=
          0.2 * $t[$field2i->{numres_p}]) {
         $results->{collate_perinst}->{p}->{min_20perc_max_coverage}->{$int_sig}++ ;
      }
   }
   $field2i = {} ;
   close(PERINST_PEP) ;

   foreach my $cov_type (keys %{$results->{collate_perinst}->{p}}) {
      my $cur_num = keys %{$results->{collate_perinst}->{p}->{$cov_type}} ;
      print "peptide interactions of type $cov_type: $cur_num\n" ;
   }


# TODO (future...090204_1133)- calc seqid given family, two sids, alnpos set

   open(PERINST_P, "PINT_collate_perinstance.20090128.out") ;
# only read in the SUMMARY LINES
   while (my $line = <PERINST_P>) {
      chomp $line ;
#pINT_LSUM      BDP91371-0_SCOP.d2oslb2 BDP91371-2_CHAIN-H      224     b.1.1.2 5       1       1       1nc4:DOF:006
      if ($line !~ /PINT_LSUM/) {next;}
      $line =~ s/PINT_LSUM\t// ;
      if ($line =~ /^#/) {
         $line =~ s/^\#// ;
         my @t = split(/\t/, $line) ;
         foreach my $j ( 0 .. $#t) {
            $field2i->{$t[$j]} = $j ; }
         next;
      }
      my @t = split(/\t/, $line) ;

#      push @{$results->{collate_perinst}->{"p"}->{histogram}->{max_l_and_p}},
#         $t[$field2i->{max_l_and_p}] ;
#
#      push @{$results->{collate_perinst}->{"p"}->{histogram}->{cum_l_and_p}},
#         $t[$field2i->{cumulative_numres_l_and_p}] ;

      my $int_sig = $t[$field2i->{SID1}]."\t".$t[$field2i->{SID2}] ;

      my $ptype = 'Pinter' ;
      if ($t[$field2i->{CHAINS}] eq 'same') {
         $ptype = 'Pintra' ; }

      $results->{collate_perinst}->{$ptype}->{all_interactions}->{$int_sig}++ ;

      if ( ($t[$field2i->{cumulative_numres_l_and_p_1}] >=
            0.2 * $t[$field2i->{numres_p_1}]) ||
           ($t[$field2i->{cumulative_numres_l_and_p_2}] >=
            0.2 * $t[$field2i->{numres_p_2}])) {
         $results->{collate_perinst}->{$ptype}->{min_20perc_cum_coverage}->{$int_sig}++;
      }

      if ( ($t[$field2i->{max_l_and_p_1}] >=
            0.2 * $t[$field2i->{numres_p_1}]) ||
           ($t[$field2i->{max_l_and_p_2}] >=
            0.2 * $t[$field2i->{numres_p_2}])) {
         $results->{collate_perinst}->{$ptype}->{min_20perc_max_coverage}->{$int_sig}++;
      }
   }
   $field2i = {} ;
   close(PERINST_P) ;

   foreach my $ptype (sort keys %{$results->{collate_perinst}}) {
      foreach my $cov_type (sort keys %{$results->{collate_perinst}->{$ptype}}) {
         my $cur_num = keys %{$results->{collate_perinst}->{$ptype}->{$cov_type}} ;
         print "$ptype interactions of type $cov_type: $cur_num\n" ;
      }
   }


# PER LIGAND ANALYSIS - get maximal d-d and d-p coverages for each ligand.
# dump info with liginfo->{MW} and liginfo->{numatoms_nonh}

   open(PERINST_PEP, "pINT_collate_perinstance.20090128.out") ;
   my $perlig_maxpcoverage ;
   while (my $line = <PERINST_PEP>) {
      chomp $line;
      if ($line =~ /^#pINT_LSUM/) {next;}
      if ($line =~ /^#SID/) {
         $line =~ s/^#// ;
         my @t = split(/\t/, $line) ;
         foreach my $j ( 0 .. $#t) {
            $field2i->{$t[$j]} = $j ; }
         next;
      }
      if ($line =~ /^#/) {next;}

      my @t = split(/\t/, $line) ;
      my $curlig = $t[$field2i->{LIG}] ;
      my (undef, $curligcode, undef) = split(/\:/, $curlig) ;
      my $cur_pcover = $t[$field2i->{numres_l_and_p}] /
                       $t[$field2i->{numres_p}] ;
      if (!exists $perlig_maxpcoverage->{p}->{$curligcode} ||
          $cur_pcover > $perlig_maxpcoverage->{p}->{$curligcode}) {
         $perlig_maxpcoverage->{p}->{$curligcode} = $cur_pcover ;
      }
   }
   $field2i = {} ;
   close(PERINST_PEP) ;


# herenow 090227_1549 - change file opens so that actual name only
#  specified one place, and keep the final collate results in the data
# directory rather than in the mishmash of the run directory.
   open(PERINST_P, "PINT_collate_perinstance.20090128.out") ;
   while (my $line = <PERINST_P>) {
      chomp $line ;
      if ($line =~ /^#PINT_LSUM/ ||
          $line =~ /^PINT_LSUM/) {
         next;
      } elsif ($line =~ /^#SID1/) {
         my @t = split(/\t/, $line) ;
         foreach my $j ( 0 .. $#t) {
            $field2i->{$t[$j]} = $j ; }
         next;
      }
      my @t = split(/\t/, $line) ;
      my $curlig = $t[$field2i->{LIG}] ;
      my (undef, $curligcode, undef) = split(/\:/, $curlig) ;
      my $cur_pcover = $t[$field2i->{numres_l_and_p_side}] /
                       $t[$field2i->{numres_p_side}] ;
      my $btype ;
      if ($t[$field2i->{CHAINS}] eq 'same') {
         $btype = "Pintra" ;
      } else {
         $btype = "Pinter" ;
      }

      if (!exists $perlig_maxpcoverage->{$btype}->{$curligcode} ||
          $cur_pcover > $perlig_maxpcoverage->{$btype}->{$curligcode}) {
         $perlig_maxpcoverage->{$btype}->{$curligcode} = $cur_pcover ;
      }
   }
   $field2i = {} ;
   close(PERINST_P) ;

   foreach my $int_type (keys %{$perlig_maxpcoverage}) {
      open(PERLIGMAXPCOVERAGE, ">perlig_maxpcoverage.$int_type.txt") ;
      foreach my $lig (sort keys %{$perlig_maxpcoverage->{$int_type}}) {
         print PERLIGMAXPCOVERAGE join("\t", $lig, $int_type,
             $perlig_maxpcoverage->{$int_type}->{$lig},
             $liginfo->{mw}->{$lig}, $liginfo->{numatoms_nonh}->{$lig}
             )."\n" ;
      }
      close(PERLIGMAXPCOVERAGE) ;
   }

}

##ORIGINAL HEADERS - CALC_PILIG.PL ; last working version from ALTO days
##ORIGINAL HEADERS - CALC_PILIG.PL ; last working version from ALTO days
##ORIGINAL HEADERS - CALC_PILIG.PL ; last working version from ALTO days
##
##
##
###1. ligbase_new.active_res -> domain assn -> alnpos assignment
###2. pibase_itnerfac_contacts_tablrs -> alnpos assignment
###3. Ligand filters - size (msd chem), drug-like (pubchem)
###4. per-aln view; per- fam pair view...
###5. Viewer script for #4 ; map a family's aln postions on to
###all fam/sf members
##
##
###HERENOW: figure out why d1jsa_ fails (a.39.1.5)
##
##use strict;
##use warnings ;
##use lib '/alto2/home/fred/sali/projects/domain_interfaces/work/PIBASE0510/progs/perl_api' ;
##use pibase ;
##use File::Basename qw/basename/ ;
##use Sys::Hostname qw/hostname/ ;
##use File::Temp qw/tempfile tempdir/ ;
##use pibase::interatomic_contacts qw/contacts_select contacts_select_inter special_contact raw_contacts_select/;
##
##
##use constant PARAM_MIN_LIGMW => 250;
##use constant PARAM_MAX_LIGMW => 1000;
##use constant CALC_CONSSCORE_FLAG => 1;
##use constant CALC_PVAL_FLAG => 1;
##
##use constant PARAM_SCPERCACC_THRESH => 7; #more than this is an exposed residue
##
##use constant PARAM_MIN_NUMCONTACTS => 1000 ;
##
###LIGBASE is built at 5A, use same cutoff for PIBASE interface residues
##use constant PARAM_INTERFACE_DIST_THRESH => 5.0 ;
##
##use constant BIN_MOD => "modSVN" ;
##
##use constant DEBUG => 0;
##use constant DEBUGALN => 0;
##use constant PARAM_DEFAULT_NUMJOBS => 20;
##use constant PARAM_DEFAULT_PRIORITY=> -2;
##use constant PARAM_DEFAULT_NODESPECS => "#\$ -l cpu500=false,cpu600=false,cpu933=false,cpu1500=false
###\$ -l pansali=1G,alto1=1G,alto2=1G,alto3=1G,diva1=1G,diva3=1G,scratch=1G" ;
###use constant DEBUGDOM => 'd1wyka_';
###use constant DEBUGFAM => 'a.129.1.1' ;
###use constant DEBUGFAM => 'a.129.1.1' ;
##use R;
##
##
##main() ;
##
##sub main {
##
##   my $usage = __FILE__ ;
##   my $mode_usage ;
##   $mode_usage->{'prep'}->{'numjobs'}++ ;
##   $mode_usage->{'prep'}->{'priority'}++ ;
##
##   $mode_usage->{'prep_sasacalc'} = {};
##   $mode_usage->{'sasacalc'} = {};
##   $mode_usage->{'assign_exp'} = {};
##
##   $mode_usage->{'assign_lig'} = {};
##   $mode_usage->{'assign_pi'} = {} ;
##
##   $mode_usage->{'collate_perfam'} = {} ;
##   $mode_usage->{'collate_instance'} = {} ;
##
##   if ($#ARGV < 0) {
##      die "$usage mode (".join(',', sort keys %{$mode_usage}).")\n" ; }
##
##   my $mode = shift @ARGV ; $mode =~ s/^-// ;
##
##   if (!exists $mode_usage->{$mode}) {
##      die "$usage mode (".join(',', sort keys %{$mode_usage}).")\n" ; }
##
##
##   my $options ;
##   my $j = 0 ;
##   while ($j <= $#ARGV) {
##      my $param = $ARGV[$j] ;
##      $param =~ s/^-// ;
##      if (!exists $mode_usage->{$mode}->{$param}) {
##         die "$param not recognized, possible options: ".
##            join(", ", sort keys %{$mode_usage->{$mode}})."\n" ;
##      }
##      my $val = $ARGV[($j + 1)] ;
##      $options->{$param} = $val ;
##      $j+=2 ;
##   }
##   
##   my $mode2sub= {
##      prep => \&run_prep ,
##      prep_sasacalc => \&run_prep_sasacalc ,
##      sasacalc => \&run_sasacalc ,
##      assign_exp => \&run_assign_exp ,
##      assign_lig => \&run_assign_lig ,
##      assign_pi => \&run_assign_pi,
##      collate_perfam => \&run_collate_perfam,
##      collate_instance_all => \&run_collate_perinstance_all,
##      collate_instance => \&run_collate_perinstance,
##   } ;
##
##   $mode2sub->{$mode}->($options) ;
##
##}
##
##
##sub run_prep_sasacalc {
##
##   my $dbh = connect_pibase_ligbase() ;
##   my ($sid, $bdp, $class) = pibase::mysql_fetchcols($dbh->{pi},
##     "SELECT subset_id, a.bdp_id, class FROM subsets as a, bdp_files as b ".
##     "WHERE subset_source_id = 1 and a.bdp_id = b.bdp_id and b.raw_pdb = 1");
##
##   my $class2sid ; my $sid2bdp ;
##   foreach my $j ( 0 .. $#{$sid}) {
##      push @{$class2sid->{$class->[$j]}}, $sid->[$j] ;
##      $sid2bdp->{$sid->[$j]} = $bdp->[$j] ;
##   }
##
##   foreach my $tclass (sort keys %{$class2sid}) {
##      if ($tclass !~ /^[a-g]/) {next;}
##      foreach my $tsid (@{$class2sid->{$tclass}}) {
##         print join("\t",$tsid, $sid2bdp->{$tsid}, $tclass)."\n" ;
##      }
##   }
##
##}
##
##sub _getinput_sasacalc {
##
##   my ($sidlist, $sid2bdp, $sid2class) ;
##   while (my $line = <STDIN>) {
##      chomp $line;
##      if ($line =~ /^#/) {next;}
##      my ($sid, $bdp, $class) = split(/\t/, $line) ;
##      push @{$sidlist}, $sid ;
##      $sid2bdp->{$sid} = $bdp ;
##      $sid2class->{$sid} = $class ;
##   }
##
##   return {
##      sids => $sidlist,
##      sid2bdp => $sid2bdp,
##      sid2class => $sid2class,
##   } ;
##
##}
##
##
##sub run_sasacalc {
##
##   my $mod_bin = BIN_MOD ;
##
##   my $in = _getinput_sasacalc() ;
##   foreach my $sid (@{$in->{sids}}) {
##      print STDERR "NOW ON: $sid\n" ;
##      my $exposed ;
##
##      my $bdpdir = pibase::sid_2_domdir($sid) ;
##      my $sid_fn = $bdpdir."/$sid.pdb" ;
##
##      my ($res_sasa, $sasa_resno_rev, $sasa_total, $sasa_errfl,
##          $atmsasa_total) = get_sasa(
##            {pdb_fn => $sid_fn,
##             surftyp => 2,
##             modeller_bin => $mod_bin}) ;
##
##      if($#{$sasa_errfl} >= 0) {
##         print STDERR "ERRORS: get_sasa() threw error on $sid ($sid_fn): ".
##         join(',', @{$sasa_errfl})."\n" ;
##         next;
##      }
##
##
##      foreach my $j ( 0 .. $#{$res_sasa->{resno}}) {
##         my $resno = $res_sasa->{resno}->[$j] ;
##         my $chain = $res_sasa->{chain}->[$j] ;
##         my $ressig = $resno."_".$chain ;
##
##         if ($res_sasa->{sc_perc}->[$j] > PARAM_SCPERCACC_THRESH) {
##            push @{$exposed}, $ressig ; }
##      }
##
##      print join("\t", $sid, $in->{sid2bdp}->{$sid},
##                  $in->{sid2class}->{$sid}, @{$exposed})."\n" ;
##
##   }
##
##}
##
##
##sub run_assign_exp {
##
##   my $in = _getinput_assign_exp() ;
##   my $fn = set_locations() ;
##   my $astral = astral_preload({ fn => $fn }) ;
##   my $pb = tod_pibase_preload() ;
##
##   foreach my $classtype (qw/fam sf/) {
##      foreach my $class (sort keys %{$astral->{classes}->{$classtype}}) {
##
##         if (!exists $in->{class2sid}->{$classtype}->{$class}) { next;}
##
##         print STDERR "NOW ON $classtype $class (line ".__LINE__.")\n" ;
##         my $pos2bs ;
##         my $class_aln = load_asteroids_aln({
##            aln_fn => $fn->{aln}->{$classtype}.'/'.$class.'.aln.fa',
##            seq_fn => $fn->{aln}->{$classtype}.'/'.$class.'.fa',
##            raf => $astral->{raf},
##            gdseqh => $astral->{gdseqh},
##            seqclcont100 => $astral->{seqcl2cont}->{100},
##            seqcl100 => $astral->{seqcl}->{100},
##            allchains => $pb->{pdbchains}
##         }) ;
##
##         foreach my $sid (sort 
##                    keys %{$in->{class2sid}->{$classtype}->{$class}} ){
##
##            my $osid = $pb->{sid2osid}->{$sid} ;
##            my $pdb = $pb->{sid2pdb}->{$sid} ;
##
##            my @talnres = () ;
##            my @undefres = () ;
##            foreach my $res (@{$in->{sid2expres}->{$sid}}) {
##               my $alnres = $class_aln->{resno2pos}->{$osid}->{$res} ;
##               if (!defined $alnres) {
###                  print STDERR "WARNING: $osid ($sid) residue $res not ".
###                               "found in ASTRAL alignment\n" ;
##                  push @undefres, 'undef';
##               } else {
##                  push @talnres, $alnres ;
##               }
###               $pos2bs->{$alnres}->{e}->{$sid}++ ;
##            }
##
##            my @salnres = ();
##            push @salnres, sort {$a <=> $b} @talnres ;
##            push @salnres, @undefres ;
##            my $alnposstring = join(',', @salnres) ;
##
##            my @outvals = ( $pdb, $sid, $osid, $classtype,
##                                 $class, "E", $alnposstring,
##                                 $class_aln->{alnlength} ) ;
##            print join("\t", @outvals)."\n" ;
##
##         }
##      }
##   }
##
##}
##
##
##sub run_collate_perinstance {
##
##   require Bit::Vector ;
##
##   my $standardres = {
##   ALA => 'A' ,
##   ARG => 'R' ,
##   ASN => 'N' ,
##   ASP => 'D' ,
##   CYS => 'C' ,
##   GLN => 'Q' ,
##   GLU => 'E' ,
##   GLY => 'G' ,
##   HIS => 'H' ,
##   HSD => 'H' ,
##   HSE => 'H' ,
##   ILE => 'I' ,
##   LEU => 'L' ,
##   LYS => 'K' ,
##   MET => 'M' ,
##   PHE => 'F' ,
##   PRO => 'P' ,
##   SER => 'S' ,
##   THR => 'T' ,
##   TRP => 'W' ,
##   TYR => 'Y' ,
##   VAL => 'V',
##   UNK => 'X',
##   '  C' => 'c',
##   '  G' => 'g',
##   '  A' => 'a',
##   '  T' => 't',
##   '  U' => 'u',
##   '  I' => 'i',
##   'C' => 'c',
##   'G' => 'g',
##   'A' => 'a',
##   'T' => 't',
##   'U' => 'u',
##   'I' => 'i',
##   '+C' => 'c',
##   '+G' => 'g',
##   '+A' => 'a',
##   '+T' => 't',
##   '+U' => 'u',
##   '+I' => 'i'
##   } ;
##
### Per-PPI
### 0. load all of ligands first - lig info
### - iterate over assignments, and have individual
###   bit vectors set for the ligand binding sites in their
###   respective domain families
###
### per-PPI: iterate over each domain--domain interface
###  1. list any ligands actually bound to this interface
###  2. iterate over all ligands in that family/superfamily, find
###     - quantify per ligand coverage
###     - quantify cumulative ligand coverage
###
### per-ligand:
##
##
##   my $fn = set_locations() ;
##   my $astral = astral_preload({ fn => $fn }) ;
##   my $pb = tod_pibase_preload() ;
##   my $liginfo = load_liginfo({fn => $fn}) ;
##
##   my $class2alnlength = {};
##
##   my $ligbits = readin_ligassignments({
##      fn => $fn,
##      liginfo => $liginfo,
##      standardres => $standardres,
##      class2alnlength => $class2alnlength,
##   }) ;
##
##   my $expbits = readin_expassignments({
##      fn => $fn,
##      standardres => $standardres,
##      class2alnlength => $class2alnlength,
##   }) ;
##
##   my $pibits_both = readin_piassignments({
##      fn => $fn,
##      pb => $pb,
##      standardres => $standardres,
##      class2alnlength => $class2alnlength,
##   }) ;
##   my $pibits = $pibits_both->{pibits} ;
##   my $interfaces = $pibits_both->{interfaces} ;
##
##   my $alltogether = combine_piligexpbits({
##      pb => $pb,
##      liginfo => $liginfo,
##      ligbits => $ligbits,
##      pibits => $pibits,
##      expbits => $expbits,
##   });
##   my $class2bits = $alltogether->{class2bits} ;
##   my $class2ligs = $alltogether->{class2ligs} ;
##   my $class2anligbits = $alltogether->{class2anligbits} ;
##
##
##   my @headers_sum=  ("PINT_LSUM", "PDB", "SID1", "SID2", "CHAINS",
##                       "CLASSTYPE", "CLASS1", "CLASS2",
##                       "numres_p_1", "cumulative_numres_l_and_p_1",
##                       "max_liginrange_1",
##                       "max_l_and_p_1", "lig_max_l_and_p_1",
##                       "max_obi_1", "lig_max_obi_1",
##                       "max_opi_1", "lig_max_opi_1",
##                       "max_olig_1", "lig_max_olig_1",
##                       "numres_p_2", "cumulative_numres_l_and_p_2",
##                       "max_liginrange_2",
##                       "max_l_and_p_2", "lig_max_l_and_p_2",
##                       "max_obi_2", "lig_max_obi_2",
##                       "max_opi_2", "lig_max_opi_2",
##                       "max_olig_2", "lig_max_olig_2",
##                       ) ;
##   print '#'.join("\t",@headers_sum)."\n" ;
##
###   foreach my $classtype (qw/sf fam/)
##   foreach my $classtype (qw/fam/) {
##      foreach my $sid12 (keys %{$interfaces}) {
##         my $sid ;
##         ($sid->{1}, $sid->{2}) = split(/\t/, $sid12) ;
##         my $pdb = $interfaces->{$sid12}->{pdb} ;
##         my $chains = $interfaces->{$sid12}->{chains} ;
##   
##         my ($alnlength, $class) ;
##         foreach my $side (1, 2) {
##            $class->{$side} = $pb->{sid2class}->{$classtype}->{$sid->{$side}} ;
##            if (!exists $class2alnlength->{$classtype}->{$class->{$side}}) {
##               next ; }
##            $alnlength->{$side} = $class2alnlength->{$classtype}->{$class->{$side}};
###            $surf_pos->{$side} = $class2bits->{$classtype}->{$class->{$side}}->{'E'}->Norm() ;
##         }
##   
##         my $bs_bits;
##         my $curligs ;
##   
##         my $alnpos ;
##         my $pistats ;
##         my $ligstats ; 
##         my $skipside ;
##         foreach my $s (1, 2) {
##            if( !exists $interfaces->{$sid12}->{$s} ||
##               !exists $interfaces->{$sid12}->{$s}->{pibits}->{$classtype} ||
##               ($interfaces->{$sid12}->{$s}->{pibits}->{$classtype}->Norm() == 0)){
##   
##               $skipside->{$s}++ ;
##               print STDERR "WARNING: skipping $sid->{1} -- $sid->{2}".
##                  " undefined binding site alignment positions for $sid->{$s}\n";
##               next;
##            }
##   
##            if (exists $class2anligbits->{$classtype}->{$class->{$s}}) {
##               $curligs->{$s} =
##                  $class2anligbits->{$classtype}->{$class->{$s}} ; }
##            $bs_bits->{$s} = $interfaces->{$sid12}->{$s}->{pibits}->{$classtype} ;
##         }
##   
##   # init vals
##         foreach my $s (1, 2) {
##            foreach my $type (qw/p cumlig_l cumlig_l_and_p max_liginrange max_l_and_p max_obi max_opi max_olig/) {
##               $ligstats->{$s}->{$type} = 0 ; }
##            foreach my $type (qw/lig_max_l_and_p lig_max_obi lig_max_opi lig_maxolig/) {
##               $ligstats->{$s}->{$type} = 'U'; }
##         }
##   
##         foreach my $s (1, 2) {
##            if (exists $skipside->{$s}) {next;}
##   
##            $pistats->{$s}->{p} = $bs_bits->{$s}->Norm() ;
##            $pistats->{$s}->{max_liginrange} = 0 ;
##            $ligstats->{$s}->{cumlig} = Bit::Vector->new($alnlength->{$s}) ;
##            foreach my $t (qw/obi opi olig l_and_p/) {
##               $ligstats->{$s}->{"max_$t"} = 0 ;
##               $ligstats->{$s}->{"lig_max_$t"} = 'undef' ;
##            }
##   
##            my $temp_bits_and = Bit::Vector->new($alnlength->{$s}) ;
##            my $temp_bits_or = Bit::Vector->new($alnlength->{$s}) ;
##   
##            if (!exists $curligs->{$s}) {next;}
##            foreach my $j (0 .. $#{$curligs->{$s}}) {
##               my $curligsig = $curligs->{$s}->[$j]->[0] ;
##               my $outcurligsig = $curligsig; $outcurligsig =~ s/\t/:/g ;
##               my $curligbits = $curligs->{$s}->[$j]->[1] ;
##               
##               $ligstats->{$s}->{cumlig}->Or($ligstats->{$s}->{cumlig},
##                                             $curligbits) ;
##   
##               #intersect the ligands bit vector and the bit vector for the bindin
##               # site positions of this particular interface;
##               $temp_bits_and->And($curligbits, $bs_bits->{$s}) ;
##               $temp_bits_or->Or($curligbits, $bs_bits->{$s}) ;
##   
##               my $curs ;
##               $curs->{l_and_p} = $temp_bits_and->Norm() ;
##               $curs->{l_or_p} = $temp_bits_or->Norm() ;
##               $curs->{l} = $curligbits->Norm() ;
##               $curs->{p} = $bs_bits->{$s}->Norm() ;
##   
##               if ($curs->{l_and_p} == 0) {next;}
##   
##               my $curoverlap ;
##               $curoverlap->{obi} = $curs->{l_and_p} / $curs->{l_or_p} ;
##               $curoverlap->{opi} = $curs->{l_and_p} / $curs->{p} ;
##               $curoverlap->{olig} = $curs->{l_and_p} / $curs->{l} ;
##   
##   
##               if (!exists $ligstats->{$s}->{"max_l_and_p"} ||
##                   $curs->{l_and_p} > $ligstats->{$s}->{"max_l_and_p"}) {
##                  $ligstats->{$s}->{"max_l_and_p"} = $curs->{l_and_p} ;
##                  $ligstats->{$s}->{"lig_max_l_and_p"} = $outcurligsig ;
##               }
##   
##               foreach my $t (qw/obi opi olig/) {
##                  if (!exists $ligstats->{$s}->{"max_$t"} ||
##                         $curoverlap->{$t} > $ligstats->{$s}->{"max_$t"}) {
##                     $ligstats->{$s}->{"max_$t"} = $curoverlap->{$t} ;
##                     $ligstats->{$s}->{"lig_max_$t"} = $outcurligsig ;
##                  }
##               }
##   
##               my (undef, $ligcod, undef) = split(/\t/, $curligsig) ;
##               my $liginrange = 0 ;
##               if (exists $liginfo->{mwinrange}->{$ligcod}) {
##                   $ligstats->{$s}->{max_liginrange} = 1 ;
##                   $liginrange = 1 ; }
##   
##               my @outvals = ("PINT_LINT",
##                                 $pdb, $sid->{1}, $sid->{2},
##                                 $chains, $classtype,
##                                 $class->{1}, $class->{2},
##                                 $s, $sid->{$s}, $outcurligsig,
##                                 $liginrange,
##                                 $curs->{p}, $curs->{l},
##                                 $curs->{l_and_p}, $curs->{l_or_p},
##                                 sprintf("%.3f", $curoverlap->{obi}),
##                                 sprintf("%.3f", $curoverlap->{opi}),
##                                 sprintf("%.3f", $curoverlap->{olig})
##                              ) ;
##               print join("\t", @outvals)."\n";
##            }
##   
##            $ligstats->{$s}->{"cumlig_l"} = $ligstats->{$s}->{cumlig}->Norm() ;
##            $temp_bits_and->And($ligstats->{$s}->{cumlig}, $bs_bits->{$s}) ;
##            $ligstats->{$s}->{"cumlig_l_and_p"} = $temp_bits_and->Norm() ;
##         }
##   
##         my @outvals =  ("PINT_LSUM", $pdb, $sid->{1}, $sid->{2}, $chains,
##                                 $classtype, $class->{1}, $class->{2}) ;
##   
##         foreach my $s ( 1, 2) {
##            if (exists $skipside->{$s} ||
##               !exists $ligstats->{$s} ||
##               $ligstats->{$s}->{"cumlig_l"} == 0 ) {
##   
##               if (exists $skipside->{$s}) {
##                  push @outvals, 0, 0, 0 ;
##               } else {
##                  push @outvals, $pistats->{$s}->{p}, 0, 0 ;
##               }
##   
##               push @outvals, 0, "U", "U", "U";
##               push @outvals, 0, "U", "U", "U";
##               push @outvals, 0, "U", "U", "U";
##               push @outvals, 0, "U", "U", "U";
##               next;
##            }
##   
##            push @outvals, $pistats->{$s}->{p} ;
##            push @outvals, $ligstats->{$s}->{cumlig_l_and_p} ;
##            push @outvals, $ligstats->{$s}->{max_liginrange} ;
##            push @outvals, $ligstats->{$s}->{"max_l_and_p"} ;
##            push @outvals, $ligstats->{$s}->{"lig_max_l_and_p"} ;
##            foreach my $t (qw/obi opi olig/)  {
##               push @outvals, sprintf("%.3f", $ligstats->{$s}->{"max_$t"}) ;
##               push @outvals, $ligstats->{$s}->{"lig_max_$t"} ;
##            }
##         }
##         print join("\t", @outvals)."\n"; 
##      }
##   }
##
##}
##
##
##sub run_collate_perfam {
##
##   require Bit::Vector ;
##   require R;
##   require RReferences ;
##   R::initR("--silent") ;
##
##   my $standardres = {
##   ALA => 'A' ,
##   ARG => 'R' ,
##   ASN => 'N' ,
##   ASP => 'D' ,
##   CYS => 'C' ,
##   GLN => 'Q' ,
##   GLU => 'E' ,
##   GLY => 'G' ,
##   HIS => 'H' ,
##   HSD => 'H' ,
##   HSE => 'H' ,
##   ILE => 'I' ,
##   LEU => 'L' ,
##   LYS => 'K' ,
##   MET => 'M' ,
##   PHE => 'F' ,
##   PRO => 'P' ,
##   SER => 'S' ,
##   THR => 'T' ,
##   TRP => 'W' ,
##   TYR => 'Y' ,
##   VAL => 'V',
##   UNK => 'X',
##   UNX => 'X',
##   '  C' => 'c',
##   '  G' => 'g',
##   '  A' => 'a',
##   '  T' => 't',
##   '  U' => 'u',
##   '  I' => 'i',
##   'C' => 'c',
##   'G' => 'g',
##   'A' => 'a',
##   'T' => 't',
##   'U' => 'u',
##   'I' => 'i',
##   '+C' => 'c',
##   '+G' => 'g',
##   '+A' => 'a',
##   '+T' => 't',
##   '+U' => 'u',
##   '+I' => 'i'
##   } ;
##
##
### given deposit dirs,
### load list of all ASTEROIDS fam/sf
##
##   my $fn = set_locations() ;
##   my $astral = astral_preload({ fn => $fn }) ;
##   my $pb = tod_pibase_preload() ;
##   my $liginfo = load_liginfo({fn => $fn}) ;
##
##   my $class2alnlength = {};
##
##   my $ligbits = readin_ligassignments({
##      fn => $fn,
##      liginfo => $liginfo,
##      standardres => $standardres,
##      class2alnlength => $class2alnlength,
##   }) ;
##
##   my $expbits = readin_expassignments({
##      fn => $fn,
##      standardres => $standardres,
##      class2alnlength => $class2alnlength,
##   }) ;
##
##   my $pibits_both = readin_piassignments({
##      fn => $fn,
##      pb => $pb,
##      standardres => $standardres,
##      class2alnlength => $class2alnlength,
##   }) ;
##   my $pibits = $pibits_both->{pibits} ;
##   my $interfaces = $pibits_both->{interfaces} ;
##
##   my $alltogether = combine_piligexpbits({
##      pb => $pb,
##      liginfo => $liginfo,
##      ligbits => $ligbits,
##      pibits => $pibits,
##      expbits => $expbits,
##   });
##   my $class2bits = $alltogether->{class2bits} ;
##   my $class2ligs = $alltogether->{class2ligs} ;
##   my $class2anligbits = $alltogether->{class2anligbits} ;
##
##   {
##      my @headers ;
##      push @headers, qw/classtype class alnlength/ ;
##      my @tt = qw/Lmwinrange P Pinter Pintra E/ ;
##      foreach my $t (@tt) {
##         push @headers, "cons_shanonne_$t" ;
##         push @headers, "cons_numtypes_$t" ;
##         push @headers, "bits_$t" ;
##         push @headers, "norm_$t" ;
##      }
##
##      foreach my $t (@tt[1..$#tt]) {
##         push @headers, "cons_shannone_".$tt[0]."_and_$t" ;
##         push @headers, "cons_numtypes_".$tt[0]."_and_$t" ;
##         push @headers, "bits_".$tt[0]."_and_$t" ;
##         push @headers, "norm_".$tt[0]."_and_$t" ;
##         push @headers, "bits_".$tt[0]."_or_$t" ;
##         push @headers, "norm_".$tt[0]."_or_$t" ;
##         push @headers, "pval_".$tt[0]."_$t" ;
##      }
##      push @headers, 'ligcodes' ;
##
##      foreach my $j (0 .. $#headers) {
##         print '#'.($j + 1)."\t".$headers[$j]."\n" ; }
##      print '#'.join("\t", @headers)."\n" ;
##   }
##
##   my $astral_classes = get_astral_classlist({fn => $fn}) ;
##   my $allprocalns ;
###   foreach my $classtype (qw/seqcl90 seqcl95 seqcl100 fam sf/)
##   if (CALC_CONSSCORE_FLAG) { open(OUTCONSSCORE, ">>outconsscore.$$.out") ;}
###   foreach my $classtype (qw/fam sf/) {
##   foreach my $classtype (qw/fam/) {
##      my @classes ;
##      my $curseqidlevel ;
##      if ($classtype eq 'fam' || $classtype eq 'sf') {
##         @classes = sort keys %{$astral_classes->{$classtype}} ;
##      } elsif ($classtype =~ /^seqcl/) {
##         my ($seqid) = ($classtype =~ /seqcl([0-9]+)/) ;
##         $curseqidlevel = $seqid ;
##         @classes = sort keys %{$astral->{seqcl2cont}->{$seqid}} ;
##      }
##
##      foreach my $class (@classes) {
##         my $outbits ;
##         my $outnorm ;
##
##         if (!exists $class2bits->{$classtype}->{$class}->{Lmwinrange} &&
##             !exists $class2bits->{$classtype}->{$class}->{P}) {
##            my @outvals = ( $classtype, $class, 'NOINTERACTIONS' );
##            print join("\t", @outvals)."\n" ;
##            next;
##         }
##
##         print STDERR "NOW ON $classtype $class (line ".__LINE__.")\n" ;
##         my $class_aln ;
##         if (CALC_CONSSCORE_FLAG) {
##            my ($curfam, $curclasstype) ;
##            if ($classtype eq 'fam' || $classtype eq 'sf') {
##               $curfam = $class;
###BUG070503_1511                $curclasstype = $class ;
##               $curclasstype = $classtype ;
##            } elsif ($classtype =~ /^seqcl/) {
##               $curfam = $pb->{osid2class}->{fam}->{$class} ;
##               $curclasstype = 'fam' ;
##            }
##
##            if (! exists $allprocalns->{$curclasstype}->{$curfam}) {
###               print STDERR "OHFUCK $curclasstype $curfam aln file is ".
###                  $fn->{aln}->{$curclasstype}.'/'.$curfam.'.aln.fa'."\n" ;
##               $allprocalns->{$curclasstype}->{$curfam} = load_asteroids_aln({
##                  aln_fn => $fn->{aln}->{$curclasstype}.'/'.$curfam.'.aln.fa',
##                  seq_fn => $fn->{aln}->{$curclasstype}.'/'.$curfam.'.fa',
##                  raf => $astral->{raf},
##                  gdseqh => $astral->{gdseqh},
##                  seqclcont100 => $astral->{seqcl2cont}->{100},
##                  seqcl100 => $astral->{seqcl}->{100},
##                  allchains => $pb->{pdbchains}
##               }) ;
##            }
##
##            $class_aln = $allprocalns->{$curclasstype}->{$curfam} ;
##         }
##
##         my $consscore_se ;
##         my $consscore_nt ;
##
##         foreach my $btype (qw/Lmwinrange Pinter Pintra P E/) {
##            if (exists $class2bits->{$classtype}->{$class}->{$btype}) {
##               $outnorm->{$btype} = 
##               $class2bits->{$classtype}->{$class}->{$btype}->Norm();
##
##               $outbits->{$btype} =
##               $class2bits->{$classtype}->{$class}->{$btype}->to_Bin();
##
###               if ($classtype eq 'fam' || $classtype eq 'sf')
##               if (CALC_CONSSCORE_FLAG && 
##                   ($classtype eq 'fam' || $classtype eq 'sf')) {
##                  my $tsum_se = 0;
##                  my $tsum_nt = 0;
##                  foreach my $tpos ($class2bits->{$classtype}->{$class}->{$btype}->Index_List_Read()) {
##                     print OUTCONSSCORE join("\t",$class, $classtype,$btype,
##                        $tpos,
##                        $class_aln->{meta}->{shannone}->{$tpos},
##                        $class_aln->{meta}->{numrestypes}->{$tpos})."\n" ;
##
##                     $tsum_se += $class_aln->{meta}->{shannone}->{$tpos} ;
##                     $tsum_nt += $class_aln->{meta}->{numrestypes}->{$tpos} ;
##                  }
##                  $consscore_se->{$btype} = $tsum_se / $outnorm->{$btype} ;
##                  $consscore_nt->{$btype} = $tsum_nt / $outnorm->{$btype} ;
##               } else {
##                  $consscore_se->{$btype} = "undefseqcl" ;
##                  $consscore_nt->{$btype} = "undefseqcl" ;
##               }
##
##            } else {
##               my $t = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
##               $consscore_se->{$btype} = "undef" ;
##               $consscore_nt->{$btype} = "undef" ;
##               $outbits->{$btype} = $t->to_Bin() ;
##               $outnorm->{$btype} = 0 ;
##            }
##         }
##
##         my $overlap_pvals ;
##         foreach my $ptype (qw/P Pinter Pintra/) {
##            $overlap_pvals->{"Lmwinrange_$ptype"} = "undef" ;
##            if ($outnorm->{Lmwinrange} > 0 && $outnorm->{$ptype} > 0) {
##               my $tor = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
##               $tor->Or(
##                  $class2bits->{$classtype}->{$class}->{$ptype},
##                  $class2bits->{$classtype}->{$class}->{Lmwinrange}) ;
##
##               $outbits->{"Lmwinrange_or_$ptype"} = $tor->to_Bin() ;
##               $outnorm->{"Lmwinrange_or_$ptype"} = $tor->Norm() ;
##
##               my $tand = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
##               $tand->And(
##                     $class2bits->{$classtype}->{$class}->{$ptype},
##                     $class2bits->{$classtype}->{$class}->{Lmwinrange}) ;
##
##               $outbits->{"Lmwinrange_and_$ptype"} = $tand->to_Bin() ;
##               $outnorm->{"Lmwinrange_and_$ptype"} = $tand->Norm() ;
##               $consscore_se->{"Lmwinrange_and_$ptype"} = "undef" ;
##               $consscore_nt->{"Lmwinrange_and_$ptype"} = "undef" ;
##
##               my $tsum_se = 0;
##               my $tsum_nt = 0;
##               $overlap_pvals->{"Lmwinrange_$ptype"} = 'pvalnotcalc' ;
##               if ($outnorm->{"Lmwinrange_and_$ptype"} > 0 ) {
##                  if (CALC_CONSSCORE_FLAG) {
##                     foreach my $tpos ($tand->Index_List_Read()) {
##                        $tsum_nt += $class_aln->{meta}->{numrestypes}->{$tpos} ;
##                        $tsum_se += $class_aln->{meta}->{shannone}->{$tpos} ;
##
##                     print OUTCONSSCORE join("\t",$class, $classtype,
##                        "Lmwinrange_and_$ptype",
##                        $tpos,
##                        $class_aln->{meta}->{shannone}->{$tpos},
##                        $class_aln->{meta}->{numrestypes}->{$tpos})."\n" ;
##
##                     }
##                     $consscore_se->{"Lmwinrange_and_$ptype"} =
##                        $tsum_se / $outnorm->{"Lmwinrange_and_$ptype"} ;
##                     $consscore_nt->{"Lmwinrange_and_$ptype"} =
##                        $tsum_nt / $outnorm->{"Lmwinrange_and_$ptype"} ;
##                  }
##
##                  if (CALC_PVAL_FLAG) {
##                     my $tnorm ;
##                     $tnorm->{lp} = $outnorm->{"E"} - $outnorm->{"Lmwinrange_or_$ptype"} ;
##                     $tnorm->{Lp} = $outnorm->{"Lmwinrange"} - $outnorm->{"Lmwinrange_and_$ptype"} ;
##                     $tnorm->{lP} = $outnorm->{$ptype} - $outnorm->{"Lmwinrange_and_$ptype"} ;
##                     $tnorm->{LP} = $outnorm->{"Lmwinrange_and_$ptype"} ;
##
##                     my @x = &R::eval("capture.output(fisher.test(matrix(c($tnorm->{lp}, $tnorm->{Lp}, $tnorm->{lP}, $tnorm->{LP}), 2, 2),workspace=2e7))") ;
##                     my $x = join("\n", @x) ;
##                     my ($pval) = ($x =~ /p-value [\<=] (.*)\nalternative/) ;
##                     $overlap_pvals->{"Lmwinrange_$ptype"} = $pval ;
##                  }
##               }
##            } elsif ($outnorm->{Lmwinrange} > 0) {
##
##                  $outbits->{"Lmwinrange_or_$ptype"} = $outbits->{Lmwinrange} ;
##                  $outnorm->{"Lmwinrange_or_$ptype"} = $outnorm->{Lmwinrange} ;
##
##                  my $t = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
##                  $consscore_se->{"Lmwinrange_and_$ptype"} = "undef" ;
##                  $consscore_nt->{"Lmwinrange_and_$ptype"} = "undef" ;
##                  $outbits->{"Lmwinrange_and_$ptype"} = $t->to_Bin() ;
##                  $outnorm->{"Lmwinrange_and_$ptype"} = 0 ;
##
##            } elsif ($outnorm->{$ptype} > 0) {
##
##                  $outbits->{"Lmwinrange_or_$ptype"} = $outbits->{$ptype} ;
##                  $outnorm->{"Lmwinrange_or_$ptype"} = $outnorm->{$ptype} ;
##
##                  my $t = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
##                  $consscore_se->{"Lmwinrange_and_$ptype"} = "undef" ;
##                  $consscore_nt->{"Lmwinrange_and_$ptype"} = "undef" ;
##                  $outbits->{"Lmwinrange_and_$ptype"} = $t->to_Bin() ;
##                  $outnorm->{"Lmwinrange_and_$ptype"} = 0 ;
##
##            } else {
##
##                  my $t = Bit::Vector->new($class2alnlength->{$classtype}->{$class}) ;
##                  $outbits->{"Lmwinrange_or_$ptype"} = $t->to_Bin() ;
##                  $outnorm->{"Lmwinrange_or_$ptype"} = 0 ;
##                  $outbits->{"Lmwinrange_and_$ptype"} = $t->to_Bin() ;
##                  $outnorm->{"Lmwinrange_and_$ptype"} = 0 ;
##                  $consscore_se->{"Lmwinrange_and_$ptype"} = "undef" ;
##                  $consscore_nt->{"Lmwinrange_and_$ptype"} = "undef" ;
##
##            }
##         }
##
##         my @outvals = ( $classtype, $class,
##                         $class2alnlength->{$classtype}->{$class});
##
##         foreach my $btype (qw/Lmwinrange P Pinter Pintra E/) {
##            push @outvals, $consscore_se->{$btype} ;
##            push @outvals, $consscore_nt->{$btype} ;
##            push @outvals, $outbits->{$btype} ;
##            push @outvals, $outnorm->{$btype} ;
##         }
##
##         foreach my $btype (qw/P Pinter Pintra/) {
##            push @outvals, $consscore_se->{"Lmwinrange_and_$btype"} ;
##            push @outvals, $consscore_nt->{"Lmwinrange_and_$btype"} ;
##            push @outvals, $outbits->{"Lmwinrange_and_$btype"} ;
##            push @outvals, $outnorm->{"Lmwinrange_and_$btype"} ;
##            push @outvals, $outbits->{"Lmwinrange_or_$btype"} ;
##            push @outvals, $outnorm->{"Lmwinrange_or_$btype"} ;
##            push @outvals, $overlap_pvals->{"Lmwinrange_$btype"} ;
##         }
##
##
###061111_1310HERENOW
##         my $ligstring = '';
##         if (exists $class2ligs->{$classtype}->{$class}->{Lmwinrange}) {
##            my @tl ;
##            foreach my $l (keys %{$class2ligs->{$classtype}->{$class}->{Lmwinrange}}){
##               my $t = $l ; $t =~ s/\t/:/g ;
##               push @tl, $t ;
##            }
##            $ligstring = join(',', @tl) ;
##         }
##         push @outvals, $ligstring ;
##         foreach my $j ( 0 .. $#outvals) {
##            if (!defined $outvals[$j]) {
##               print STDERR " WARNING: $classtype $class field $j is undefined\n "; } }
##
##         print join("\t", @outvals)."\n" ;
##      }
##   }
##   if (CALC_CONSSCORE_FLAG) { close(OUTCONSSCORE) ;}
##}
##
##
##sub _getinput_assign_pi {
##
##   my $bdp2contacts_fn;
##   my $bdp2sid;
##   my $class2sid ;
##   my $sid2class ;
##   while (my $line = <STDIN>) {
##      if ($line =~ /^\#/) {next;} ;
##      chomp $line ;
##      my ($bdp, $contacts_fn, $sid1, $class1, $sid2, $class2) =
##         split(/\t/, $line);
##      $bdp2contacts_fn->{$bdp} = $contacts_fn ;
##      $class2sid->{$class1}->{$sid1}++ ;
##      $class2sid->{$class2}->{$sid2}++ ;
##      $sid2class->{$sid1}->{$class1}++ ;
##      $sid2class->{$sid2}->{$class2}++ ;
##   }
##
##   return {
##      sid2class => $sid2class,
##      class2sid => $class2sid,
##      bdp2contacts_fn=> $bdp2contacts_fn
##   } ;
##
##}
##
##sub _getinput_assign_lig {
##
##   my $pdb2res2ligid ;
##   while (my $line = <STDIN>) {
##      if ($line =~ /^\#/) {next;}
##      chomp $line ;
##      my ($pdb, $resno, $resch, $ligcod, $ligid) = split(/\t/, $line) ;
##      $pdb2res2ligid->{$pdb}->{$resno."\n".$resch}->{$ligcod."\n".$ligid}++;
##   }
##
##   return {
##      pdb2res2ligid => $pdb2res2ligid
##   } ;
##}
##
##
##sub _getinput_assign_exp {
##
##   my $sid2class ; my $class2sid ;
##   my $sid2bdp ;
##   my $sid2expres ;
##   while (my $line = <STDIN>) {
##      if ($line =~ /^\#/) {next;}
##      chomp $line ;
##      my @t = split(/\t/, $line) ;
##      my $sid = shift @t ;
##      my $bdp = shift @t ;
##      my $class = shift @t ;
##
##      foreach my $res ( @t) {
##         substr($res,-2,1) = "\n" ;
##         push @{$sid2expres->{$sid}}, $res ;
##      }
##
##      $sid2class->{$sid} = $class ;
##
##      my ($sf) = ($class =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
##      $sid2class->{fam}->{$sid} = $class;
##      $sid2class->{sf}->{$sid} = $sf ;
##      $class2sid->{fam}->{$class}->{$sid}++ ;
##      $class2sid->{sf}->{$sf}->{$sid}++ ;
##   }
##
##   return {
##      sid2class => $sid2class,
##      class2sid => $class2sid,
##      sid2expres=> $sid2expres,
##   } ;
##
##}
##
##sub run_assign_lig {
##
##   my $in = _getinput_assign_lig() ;
##   my $lb = { pdb2res2ligid => $in->{pdb2res2ligid} } ;
##   my $fn = set_locations() ;
##   my $astral = astral_preload({ fn => $fn }) ;
##   my $pb = tod_pibase_preload() ;
##   my $lbdoms = get_ligbase_domains({lb => $lb, pb => $pb}) ;
##
##
##   foreach my $classtype (qw/fam sf/) {
##      foreach my $class (sort keys %{$astral->{classes}->{$classtype}}) {
##         if (!exists $lbdoms->{class2ligsid}->{$classtype}->{$class}) {next;}
##
##         print STDERR "NOW ON $classtype $class (line ".__LINE__.")\n" ;
##         my $pos2bs ;
##         my $class_aln = load_asteroids_aln({
##            aln_fn => $fn->{aln}->{$classtype}.'/'.$class.'.aln.fa',
##            seq_fn => $fn->{aln}->{$classtype}.'/'.$class.'.fa',
##            raf => $astral->{raf},
##            gdseqh => $astral->{gdseqh},
##            seqclcont100 => $astral->{seqcl2cont}->{100},
##            seqcl100 => $astral->{seqcl}->{100},
##            allchains => $pb->{pdbchains}
##         }) ;
##
##         foreach my $ligsid (sort keys %{$lbdoms->{class2ligsid}->{$classtype}->{$class}}) {
##            my $osid = $pb->{sid2osid}->{$ligsid} ;
##            my $pdb = $pb->{sid2pdb}->{$ligsid} ;
##            my ($talnres, $undefres) ;
##            foreach my $res (sort keys %{$lbdoms->{sid2ligres}->{$ligsid}}) {
##               my $alnres = $class_aln->{resno2pos}->{$osid}->{$res} ;
##               if (!defined $alnres) {
###                  print STDERR "WARNING: $osid ($ligsid) residue $res not ".
###                               "found in ASTRAL alignment\n" ;
##                  $alnres = 'undef' ;
##               }
###               $pos2bs->{$alnres}->{l}->{$ligsid}++ ;
##               foreach my $lig (sort keys %{$lb->{pdb2res2ligid}->{$pdb}->{$res}}) {
##                  if ($alnres ne 'undef' ) {
##                     push @{$talnres->{$lig}}, $alnres ;
##                  } else {
##                     push @{$undefres->{$lig}}, $alnres ;
##                  }
##               }
##            }
##
##            foreach my $lig (keys %{$talnres}) {
##               my ($ligcod, $ligid) = split(/\n/, $lig) ;
##
##               my @salnres = ();
##               if (exists $talnres->{$lig}) {
##                  push @salnres, sort {$a <=> $b} @{$talnres->{$lig}} ; }
##               if (exists $undefres->{$lig}) {
##                  push @salnres, @{$undefres->{$lig}} ; }
##               my $alnposstring = join(',', @salnres) ;
##
##               my @outvals = ( $pdb, $ligsid, $osid, $classtype,
##                                 $class, "L", $alnposstring,
##                                 $class_aln->{alnlength},
##                                 $ligcod, $ligid) ;
##               print join("\t", @outvals)."\n" ;
##            }
##         }
##      }
##   }
##
##}
##
##
##sub run_assign_pi {
##
##   my $in = _getinput_assign_pi() ;
##   my $fn = set_locations() ;
##   my $astral = astral_preload({ fn => $fn }) ;
##   my $pb = tod_pibase_preload() ;
##
##   foreach my $classtype (qw/fam sf/) {
##      foreach my $class (sort keys %{$astral->{classes}->{$classtype}}) {
##
##         if (!exists $pb->{class2sidpairs}->{$classtype}->{$class}) { next;}
##
##         print STDERR "NOW ON $classtype $class (line ".__LINE__.")\n" ;
##         my $pos2bs ;
##         my $class_aln = load_asteroids_aln({
##            aln_fn => $fn->{aln}->{$classtype}.'/'.$class.'.aln.fa',
##            seq_fn => $fn->{aln}->{$classtype}.'/'.$class.'.fa',
##            raf => $astral->{raf},
##            gdseqh => $astral->{gdseqh},
##            seqclcont100 => $astral->{seqcl2cont}->{100},
##            seqcl100 => $astral->{seqcl}->{100},
##            allchains => $pb->{pdbchains}
##         }) ;
##
##         foreach my $j (sort {$a <=> $b}
##                    keys %{$pb->{class2sidpairs}->{$classtype}->{$class}} ){
##
##            my $pdb = $pb->{sid2pdb}->{$pb->{sid1}->[$j]} ;
##
##            my $interface = load_interface_contacts({
##               sid1 => $pb->{sid1}->[$j],
##               sid2 => $pb->{sid2}->[$j],
##               fn => $pb->{bdp2contactsfn}->{$pb->{bdp_id}->[$j]}
##            }) ;
##            my $sid12 = $pb->{sid1}->[$j]."\t".$pb->{sid2}->[$j] ;
##
##            my @dothese = () ;
##            if ($pb->{sid2class}->{$classtype}->{$pb->{sid1}->[$j]} eq
##                $class) {push @dothese, $pb->{sid1}->[$j];}
##
##            if ($pb->{sid2class}->{$classtype}->{$pb->{sid2}->[$j]} eq
##                $class) {push @dothese, $pb->{sid2}->[$j];}
##
##            foreach my $sid (@dothese) {
##               my $osid = $pb->{sid2osid}->{$sid} ;
##               my @talnres = () ;
##               my @undefres = () ;
##               foreach my $res (sort keys %{$interface->{intres}->{$sid}}) {
##                  my $alnres = $class_aln->{resno2pos}->{$osid}->{$res} ;
##                  if (!defined $alnres) {
##                     push @undefres, 'undef' ;
###                     print STDERR "WARNING: $osid ($sid) residue $res not ".
###                                  "found in ASTRAL alignment\n" ;
##                  } else {
##                     push @talnres, $alnres ;
##                  }
###                  $pos2bs->{$alnres}->{p}->{$sid}++ ;
##               }
##               my @salnres;
##               push @salnres, sort {$a <=> $b} @talnres ;
##               push @salnres, @undefres ;
##               my $alnposstring = join(',', @salnres) ;
##               my @outvals = ($pdb,  $sid, $osid, $classtype,
##                              $class, "P", $alnposstring,
##                              $class_aln->{alnlength},
##                              $pb->{sid1}->[$j], $pb->{sid2}->[$j],
##                              $pb->{sid2class}->{'fam'}->{$pb->{sid1}->[$j]},
##                              $pb->{sid2class}->{'fam'}->{$pb->{sid2}->[$j]},
##                              $pb->{sid12chain}->{$sid12},
##               ) ;
##               print join("\t", @outvals)."\n" ;
##            }
##
##         }
##      }
##   }
##
##}
##
##
##sub run_prep {
##   my $in = shift ;
##
##   my $clustspecs ;
##   $clustspecs->{nodespecs} = PARAM_DEFAULT_NODESPECS ;
##   if (!exists $in->{numjobs}) {
##      $clustspecs->{numjobs} = PARAM_DEFAULT_NUMJOBS ;
##   } else {
##      $clustspecs->{numjobs} = $in->{numjobs} ;
##   }
##
##   if (!exists $in->{priority}) {
##      $clustspecs->{priority} = PARAM_DEFAULT_PRIORITY;
##   } else {
##      $clustspecs->{priority} = $in->{priority} ;
##   }
##
##   my $fn = set_locations() ;
##   my $dbh = connect_pibase_ligbase() ;
##   my $astral = astral_preload({fn => $fn});
##
##
##
##   my $splits_dir = tempdir($fn->{pilig}->{prep_dir}."/assign_splits.XXXXX") ;
##   my $sgeout_dir_lb = tempdir($fn->{pilig}->{assign_lig_dir}."/SGEassign_lig.XXXXX");
##   my $sgeout_dir_pb = tempdir($fn->{pilig}->{assign_pi_dir}."/SGEassign_pi.XXXXX");
##
##
##   print STDERR "preparing LIGBASE data: " ;
##   my $lb = ligbase_preload({dbh => $dbh->{lig}}) ;
##   my $lb_sge_fn = "assign_lig.$$.SGE.sh" ;
##   split_ligs({lb => $lb,
##               clustspecs => $clustspecs,
##               splits_dir => $splits_dir,
##               out_fn_prefix => "ligbase_active_res",
##               SGE_dir => $sgeout_dir_lb,
##               SGE_fn => $lb_sge_fn }) ;
##   print STDERR "X\n" ;
##
##   print STDERR "preparing PIBASE data: " ;
##   my $pb = pibase_preload({dbh => $dbh->{pi}, astral => $astral}) ;
##   my $pb_sge_fn = "assign_pi.$$.SGE.sh" ;
##   split_pi({  pb => $pb,
##               clustspecs => $clustspecs,
##               splits_dir => $splits_dir,
##               SGE_dir => $sgeout_dir_pb,
##               out_fn_prefix => "bdp_ids",
##               SGE_fn => $pb_sge_fn}) ;
##   print STDERR "X\n" ;
##
##   print STDERR "SGE scripts:\nqsub $pb_sge_fn\nqsub $lb_sge_fn\n ";
##}
##
##
##sub split_ligs {
##   my $in = shift ;
##
##   my $lb = $in->{lb} ;
##   my $total_numpdb = keys %{$lb->{pdb2res2ligid}} ;
##
##   my $splitlines = POSIX::ceil($total_numpdb / $in->{clustspecs}->{numjobs}) ;
##   my @splitfiles ;
##
##   my $cur_splitnum = 1 ;
##   my $cur_fn = "split.".$in->{out_fn_prefix}.".$cur_splitnum" ;
##   push @splitfiles, $cur_fn ;
##   open(OUTF, ">$in->{splits_dir}/$cur_fn") ;
##
##   my $tasklist = "'$cur_fn'" ;
##
##   my $num_pdb = 0;
##   foreach my $pdb (keys %{$lb->{pdb2res2ligid}}) {
##      if ($num_pdb > $splitlines) {
##         $cur_splitnum++ ;
##         $cur_fn = "split.".$in->{out_fn_prefix}.".$cur_splitnum" ;
##         close(OUTF) ;
##         open(OUTF, ">$in->{splits_dir}/$cur_fn") ;
##         push @splitfiles, $cur_fn ;
##         $tasklist .= " '$cur_fn'" ;
##         $num_pdb = 0 ;
##      }
##
##      foreach my $res ( keys %{$lb->{pdb2res2ligid}->{$pdb}}) {
##         my ($resno, $resch) = split(/\n/, $res) ;
##         foreach my $lig ( keys %{$lb->{pdb2res2ligid}->{$pdb}->{$res}}) {
##            my ($ligcod, $ligid) = split(/\n/, $lig) ;
##            my @outvals = ($pdb, $resno, $resch, $ligcod, $ligid) ;
##            print OUTF join("\t", @outvals)."\n";
##         }
##      }
##      $num_pdb++ ;
##   }
##   close(OUTF) ;
##   my $numjobs = $cur_splitnum ;
##
##   my $perlscript_fn = __FILE__;
##
##   open (SGEFH, ">$in->{SGE_fn}") ;
##   print SGEFH "#!/bin/csh
###\$ -S /bin/csh
###\$ -cwd
###\$ -o $in->{SGE_dir}
###\$ -e $in->{SGE_dir}
###\$ -r y
##$in->{clustspecs}->{nodespecs}
###\$ -p $in->{clustspecs}->{priority}
###\$ -t 1-$numjobs
##
##set tasks1=( $tasklist )
##
##set input1=\$tasks1[\$SGE_TASK_ID\]
##
##set curdir=`pwd`
##set curhost=`hostname`
##set curtime=`date`
##set scratchdir=/scratch/fred/\$input1.\$\$\
##
##rm -rf \$scratchdir
##mkdir -p \$scratchdir
##
##cp $perlscript_fn \$scratchdir
##cp $in->{splits_dir}/\$input1 \$scratchdir
##
##cd \$scratchdir
##
##echo \"#sgejob run started on \$curhost at \$curtime\"
##perl $perlscript_fn -assign_lig < \$input1
##
##set curtime=`date`
##echo \"#sgejob run finished on \$curhost at \$curtime\"
##
##rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
##cd \$curdir
##rmdir \$scratchdir\n" ;
##
##   close(SGEFH) ;
##
##
##
###   return {
###      numjobs => $cur_splitnum,
###      splitfiles => \@splitfiles,
###      tasklist => $tasklist
###   } ;
##}
##
##
##sub split_pi {
##
##   my $in = shift ;
##   my $pb = $in->{pb} ;
##   my $numbdp = keys %{$pb} ;
##
##   my ($tempfh, $tempfn) = tempfile("bdp_ids.XXXXX") ;
##   foreach my $bdp (keys %{$in->{pb}->{bdp2sid12}}) {
##      foreach my $sid12 (keys %{$in->{pb}->{bdp2sid12}->{$bdp}}) {
##         my ($sid1, $sid2) = split(/\t/, $sid12) ;
##         my @outvals = ($bdp, $in->{pb}->{bdp2contactsfn}->{$bdp},
##                        $sid1, $in->{pb}->{sid2class}->{fam}->{$sid1},
##                        $sid2, $in->{pb}->{sid2class}->{fam}->{$sid2}) ;
##         print {$tempfh} join("\t", @outvals)."\n";
##      }
##   }
##   close($tempfh) ;
##
##   my $splits = _clust_split_ins({
##      fn => $tempfn,
##      dir => $in->{splits_dir},
##      numjobs => $in->{clustspecs}->{numjobs}
##   });
##   unlink $tempfn ;
##
##   my $perlscript_fn = __FILE__;
##
##   open (SGEFH, ">$in->{SGE_fn}") ;
##   print SGEFH "#!/bin/csh
###\$ -S /bin/csh
###\$ -cwd
###\$ -o $in->{SGE_dir}
###\$ -e $in->{SGE_dir}
###\$ -r y
##$in->{clustspecs}->{nodespecs}
###\$ -p $in->{clustspecs}->{priority}
###\$ -t 1-$splits->{numjobs}
##
##set tasks1=( $splits->{tasklist} )
##
##set input1=\$tasks1[\$SGE_TASK_ID\]
##
##set curdir=`pwd`
##set curhost=`hostname`
##set curtime=`date`
##set scratchdir=/scratch/fred/\$input1.\$\$\
##
##rm -rf \$scratchdir
##mkdir -p \$scratchdir
##
##cp $perlscript_fn \$scratchdir
##cp $in->{splits_dir}/\$input1 \$scratchdir
##
##cd \$scratchdir
##
##echo \"#sgejob run started on \$curhost at \$curtime\"
##perl $perlscript_fn -assign_pi < \$input1
##
##set curtime=`date`
##echo \"#sgejob run finished on \$curhost at \$curtime\"
##
##rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
##cd \$curdir
##rmdir \$scratchdir\n" ;
##
##   close(SGEFH) ;
##
##
###   return {
###      splitfiles => $splits->{splitfiles},
###      numjobs => $splits->{numjobs},
###      tasklist => $splits->{tasklist},
###   } ;
##
##}
##
##sub _clust_split_ins {
##
##   use File::Basename qw/basename/ ;
##
##   my $in = shift ;
##
##   if (!exists $in->{fn}) {
##      die "_clust_split_ins: input file not specified\n" ; }
##
##   if (!exists $in->{dir}) {
##      $in->{dir} = "./" ; }
##
##   my $inbase = basename($in->{fn}) ;
##
##   my $header_lines = '';
##   my $num_lines = 0 ;
##   open(INF, $in->{fn}) ;
##   while (my $line = <INF>) {
##      if ($line =~ /^#SET /) {
##         $header_lines .= $line ;
##      } else {
##         $num_lines++ ;
##      }
##   }
##   close(INF);
##
##   my $splitlines = POSIX::ceil($num_lines / $in->{numjobs});
##   my @splitfiles ;
##
##   $num_lines = 0 ;
##
##   my $tasklist = '';
##   my $cur_splitnum = 1 ;
##
##   my $cur_fn = "split.$inbase.$cur_splitnum" ;
##   push @splitfiles, $cur_fn ;
##   $tasklist .= "'$cur_fn' " ;
##
##   open(INF, $in->{fn}) ;
##   open(OUTF, ">$in->{dir}/$cur_fn") ;
##   print OUTF $header_lines ;
##   while (my $line = <INF>) {
##      if ($line =~ /^#SET/) {next;}
##      if ($num_lines > $splitlines) {
##         close(OUTF) ;
##         $cur_splitnum++ ;
##         $cur_fn = "split.$inbase.$cur_splitnum" ;
##         push @splitfiles, $cur_fn ;
##         $tasklist .= "'$cur_fn' " ;
##         open(OUTF, ">$in->{dir}/$cur_fn") ;
##         print OUTF $header_lines ;
##         $num_lines = 0 ;
##      }
##      print OUTF $line ;
##      $num_lines++ ;
##   }
##   close(OUTF) ;
##   close(INF);
##   $tasklist =~ s/ $// ;
###   print STDERR "OUT TASKLIST: $tasklist\n" ;
##
##
##   return {
##      numjobs => $cur_splitnum,
##      splitfiles => \@splitfiles,
##      tasklist => $tasklist
##   } ;
##}
##
##
##
##sub astral_preload {
##
##   my $in = shift;
##   my $fn = $in->{fn} ;
##
##   print STDERR "Load ASTRAL fasta headers (gdseq): " ;
##   my $astral = load_astral_headers({fn => $fn->{gdseq}}) ;
##   print STDERR "X\n" ;
##
##   my $astral_classes = get_astral_classlist({fn => $fn}) ;
##   $astral->{classes} = $astral_classes ;
##
##   print STDERR "Load ASTRAL raf: " ;
##   my $raf = raf_preload({fn => $fn->{astral_raf}}) ;
##   $astral->{raf} = $raf ;
##   print STDERR "X\n" ;
##
##   load_astral_clusters({fn => $fn, out => $astral}) ;
##
##   return $astral ;
##
##}
##
##sub set_locations {
##
##   my $fn ;
##   $fn->{astraldir} = "/alto2/home/fred/sali/projects/domain_interfaces/work/PIBASE0510/data/preexisting/ASTRAL/" ;
##
##   $fn->{astral_raf} = $fn->{astraldir}.
##      "/astral-rapid-access-1.69.raf" ;
##
##   $fn->{gdseq}= $fn->{astraldir}.'/scopseq-1.69/'.
##      'astral-scopdom-seqres-gd-all-1.69.fa' ;
##
##   foreach my $seqid ( qw/100 95 90 70 50 30 20 10/) {
##      $fn->{seqcl}->{$seqid} = $fn->{astraldir}.'/scopseq-1.69/'.
##      '/astral-scopdom-seqres-gd-sel-gs-bib-verbose-'.$seqid.'-1.69.txt' ;}
##
##   $fn->{aln}->{fam} = $fn->{astraldir}."/aln/fam/" ;
##   $fn->{aln}->{sf} = $fn->{astraldir}."/aln/sf/" ;
##
##   $fn->{pilig}->{rootdir} = "/netapp/home/fred/pilig" ;
##   $fn->{pilig}->{park_rootdir} = "/park1/fred/pilig" ;
##   $fn->{pilig}->{prep_dir} = $fn->{pilig}->{rootdir}."/prepin";
##   $fn->{pilig}->{assign_lig_dir} = $fn->{pilig}->{rootdir}."/assign_lig";
##   $fn->{pilig}->{assign_pi_dir} = $fn->{pilig}->{rootdir}."/assign_pi";
##
##   $fn->{pilig}->{assign_lig_fn} = $fn->{pilig}->{rootdir}."/assign_lig.out";
##   $fn->{pilig}->{assign_pi_fn} = $fn->{pilig}->{rootdir}."/assign_pi.out";
##   $fn->{pilig}->{assign_exp_fn} = $fn->{pilig}->{rootdir}."/assign_exp.out";
##
##   $fn->{pilig}->{park_assign_lig_dir}=$fn->{pilig}->{park_rootdir}."/assign_lig";
##   $fn->{pilig}->{park_assign_pi_dir}=$fn->{pilig}->{park_rootdir}."/assign_pi";
##
##   $fn->{pilig}->{park_assign_lig_fn}=$fn->{pilig}->{park_rootdir}."/assign_lig.out";
##   $fn->{pilig}->{park_assign_pi_fn}=$fn->{pilig}->{park_rootdir}."/assign_pi.out";
##   $fn->{pilig}->{park_assign_exp_fn}=$fn->{pilig}->{park_rootdir}."/assign_exp.out";
##
##   
##   $fn->{msdchem}->{liginfo} =$fn->{pilig}->{park_rootdir}."/aux/ligandinfo.msdchem.txt" ;
##
##   return $fn ;
##
##}
##
##sub connect_pibase_ligbase {
##
##   require DBI ;
##
##   my $dbh ;
##
##   ($dbh->{lig}) = pibase::connect_pibase({
##      'db' => 'ligbase_new',
##      'host' => 'modbase',
##      'user' => 'modbase',
##      'pass' => 'modbasesecret'
##   }) ;
##
##   ($dbh->{pi}) = pibase::connect_pibase() ;
##
##   return $dbh ;
##}
##
##
##sub compare_residue_sets {
##
##   my $in = shift;;
##   my $aln = $in->{aln} ;
##   my $sid1 = $in->{sid1} ;
##   my $sid2 = $in->{sid2} ;
##   my $res1 = $in->{res1} ;
##   my $res2 = $in->{res2} ;
##
##   my $intpos ;
##   foreach my $res1 (keys %{$res1}) {
##      my $pos = $aln->{resno2pos}->{$sid1}->{$res1} ;
##      if (!defined $pos) {
##         print STDERR "ERROR: compare_residue_sets(): $res1 from $sid1 undefined pos; aborting $sid1 -- $sid2 comparison\n" ; return 0;
##      } else {
##         $intpos->{$pos}++ ;
##      }}
##
##
##   foreach my $res2 (keys %{$res2}) {
##      my $pos = $aln->{resno2pos}->{$sid2}->{$res2} ;
##      if (!defined $pos) {
##         print STDERR "ERROR: compare_residue_sets(): $res2 from $sid2 undefined pos; aborting $sid1 -- $sid2 comparison\n" ; return 0;
##      } else {
##         $intpos->{$pos}++ ;
##      }}
##
##   my $union = 0 ;
##   my $diff = 0 ;
##   foreach my $pos (keys %{$intpos}) {
##      if ($intpos->{$pos} == 2) {
##         $union++;
##      } else {
##         $diff++ ; }
##   }
##
##   if (($union + $diff) > 0) {
##      return ($union / ($union + $diff)) ;
##   } else {
##      return 0 ;
##   }
##
##}
##
##
##sub load_asteroids_aln {
##
##   my $in = shift ;
##
##   my $gdseqh = $in->{gdseqh} ;
##   my $seqclcont100 = $in->{seqclcont100} ;
##   my $seqcl100 = $in->{seqcl100} ;
##
##   my $data = read_asteroids_aln({
##      aln_fn => $in->{aln_fn},
##      seq_fn => $in->{seq_fn},
##      allchains => $in->{allchains}
##   });
##
##   my $aln = parse_aln_raf({
##      alndata => $data,
##      raf => $in->{raf}
##   }) ;
##   $aln->{alnlength} = $data->{alnlength} ;
##   $aln->{meta} = $data->{meta} ;
##
##   return $aln ;
##
##}
##
##
##sub get_raf_line {
##   my $in = shift ;
##
##   my $pdb = $in->{pdb_id} ;
##   my $ch = $in->{chain} ;
##   my $fn = $in->{fn} ;
##
###   if ($ch eq '-') {$ch = '_';}
##   if ($ch eq ' ') {$ch = '_';}
##   my $pdbch = $pdb.$ch ;
##
##   my $tcom = "grep '^$pdbch' $fn" ;
##   my $line = `$tcom` ;
##   chomp $line;
###   if ($pdb eq '1a0w') {
###      print STDERR "**grepped $tcom to get $line\n" ; }
##
##   return $line ;
##}
##
##
##sub raf_preload {
##   my $in = shift ;
##
##   my $raf ;
##   open (RAF, $in->{fn}) ;
##   while (my $line = <RAF>) {
##      chomp $line;
##
##      if ($line =~ /^# Header:/) {
##         ($raf->{header_length}) = ($line =~ /: ([0-9]+) Bytes/) ;
##      } elsif ($line !~ /^#/) {
##         my $pdbch = substr($line, 0, 5) ;
##         $raf->{$pdbch} = $line ;
##      }
##   }
##   close(RAF) ;
##   return $raf ;
##}
##
##
##sub parse_raf_line {
##
##   my $in = shift ;
##   my $line = $in->{line}  ;
##   my $headlength = $in->{headlength} ;
##
##   my $res ;
##
##   my $head = substr($line, 0, $headlength) ;
##   $res->{atomresno_first} = substr($head, -10, 5) ;
##   $res->{atomresno_first} =~ s/ //g ;
##   $res->{atomresno_last} = substr($head, -5, 5) ;
##   $res->{atomresno_last} =~ s/ //g ;
##   my $pos = $headlength ;
##   my $seqresno = 0 ;
##   my $firstseqresno ;
##
##   my @missing_atomno = ();
##
##   while ($pos < length($line)) {
##      my $t_res = substr($line,$pos,7) ;
##      my $atomresno = substr($t_res,0,4) ; $atomresno =~ s/ //g ;
##      my $inscode = substr($t_res,4,1) ;
##      my $atomresna = substr($t_res,5,1) ;
##      my $seqresna = substr($t_res,6,1) ;
##
##      if ($inscode ne ' ') {$atomresno .= $inscode ; }
##      push @{$res->{atomresna}}, $atomresna ;
##      push @{$res->{seqresna}}, $seqresna ;
##      push @{$res->{atomresno}}, $atomresno;
##
##      if ($seqresna ne '.') {
##         $seqresno++ ;
##         push @{$res->{seqresno}}, $seqresno;
##         $res->{seqresno2ind}->{$seqresno}= $#{$res->{seqresno}};
##         $res->{ind2seqresno}->{$#{$res->{seqresno}}} = $seqresno  ;
##      } else {
##         if ($#{$res->{seqresno}} >= 0 ) {
##            $res->{atom2seqresno_back}->{$atomresno} = ${$res->{seqresno}}[-1] ;
##            $res->{atomresno2ind_back}->{$atomresno} = $#{$res->{seqresno}} ;
##         }
##
##         push @{$res->{seqresno}}, '.' ;
##         push @missing_atomno, $atomresno ;
##      }
##
##      if ($seqresna ne '.' && ($#missing_atomno >= 0)) {
##         foreach my $t_atomresno (@missing_atomno) {
##            $res->{atom2seqresno_ahead}->{$t_atomresno} = $seqresno ;
##            $res->{atomresno2ind_ahead}->{$t_atomresno} = $#{$res->{seqresno}} ;
##         }
##         @missing_atomno = () ;
##      }
##
##      if ($seqresna ne '.' && $atomresna ne '.') {
##         $res->{seq2atomresno}->{$seqresno} = $atomresno ;
##         $res->{atom2seqresno}->{$atomresno} = $seqresno ;
##         $res->{atomresno2ind}->{$atomresno} = $#{$res->{seqresno}} ;
##         $res->{ind2atomresno}->{$#{$res->{seqresno}}} = $atomresno  ;
##      }
##      $pos += 7 ;
##   }
##
##   return $res ;
##}
##
##
##sub read_asteroids_aln {
##
##   my $in = shift;
##   my $fn;
##   $fn->{aln} = $in->{aln_fn} ;
##   $fn->{seq} = $in->{seq_fn} ;
##   my $allchains = $in->{allchains} ;
##
##   my $data;
##   my $cur_dom = '' ;
##   my $seq ;
##   open(SEQF, $fn->{seq}) ;
###070503_1511   print STDERR "tried to open $fn->{seq}\n" ;
##   while (my $line = <SEQF>) {
##      chomp $line;
##      if ($line =~ /^>/) {
##         $line =~ s/^>// ;
##         my ($t_dom, undef, undef, undef) = split(/\s/, $line) ;
##         $data->{seq}->{$t_dom} = '';
##         $cur_dom = $t_dom ;
##      } else {
##         $data->{seq}->{$cur_dom} .= $line ;
##      }
##   }
##   close(SEQF) ;
##
##   $cur_dom = '' ;
##   open(ALNF, $fn->{aln}) ;
###   print STDERR "parsing $fn->{aln}\n" ;
##   while (my $line = <ALNF>) {
##      chomp $line;
##      if ($line =~ /^>/) {
##         $line =~ s/^>// ;
##         my ($t_dom, $t_class, $t_def, undef) = split(/\s/, $line) ;
##
##         $t_def =~ s/^\(// ; $t_def =~ s/\)$// ;
##         $data->{defstring}->{$t_dom} = $t_def ;
##         $data->{class}->{$t_dom} = $t_class;
##         $data->{aln}->{$t_dom} = '';
##         $data->{pdb}->{$t_dom} = substr($t_dom, 1, 4) ;
##         {
##            my @t = split(/\,/, $t_def) ;
##            my (@ch, @b, @e) ;
##            foreach my $j ( 0 .. $#t) {
##               my $t_frag ;
##               $t_frag->{chain} = '_' ;
##               if ($t[$j] =~ /:/) {
##                  $t_frag->{chain} = substr($t[$j],0,1) ; }
##               $t[$j] =~ s/^.\:// ;
##               my ($b, $e) = (' ', ' ');
##               if ($t[$j] =~ /.+\-.+/) {
##                  ($t_frag->{b}, $t_frag->{e}) =
##                  ($t[$j] =~ /(.+)\-(.+)/) ; }
##               push @{$data->{frags}->{$t_dom}}, $t_frag ;
##
##               if ($t_frag->{chain} eq '-' || $t_frag->{chain} eq '_') {
##                  $t_frag->{chain} = ' ' ;
##               } elsif (!exists $allchains->{$data->{pdb}->{$t_dom}}->{$t_frag->{chain}}) {
##                  if ($t_frag->{chain} eq uc($t_frag->{chain})) {
##                     my $lc = lc($t_frag->{chain}) ;
##                     if (exists  $allchains->{$data->{pdb}->{$t_dom}}->{$lc}) {
##                        print STDERR "WARNING defstring change: dom $t_dom old chain $t_frag->{chain} changed to $lc\n" ;
##                        $t_frag->{chain} = $lc ; }
##                  } elsif ($t_frag->{chain} eq lc($t_frag->{chain})) {
##                     my $uc = uc($t_frag->{chain}) ;
##                     if (exists $allchains->{$data->{pdb}->{$t_dom}}->{$uc}) {
##                        print STDERR "WARNING defstring change: dom $t_dom old chain $t_frag->{chain} changed to $uc\n" ;
##                        $t_frag->{chain} = $uc ; }
##                  }
##               }
##            }
##         }
##         $cur_dom = $t_dom ;
##      } else {
##         $data->{aln}->{$cur_dom} .= $line ;
##      }
##   }
##   $data->{alnlength} = length($data->{aln}->{$cur_dom}) ;
##   close(ALNF) ;
##
##   my $seqweights ;
##   foreach my $pos ( 0 .. ($data->{alnlength} - 1) ) {
###      my $freq = {
###         'A' => 0,
###         'R' => 0,
###         'N' => 0,
###         'D' => 0,
###         'C' => 0,
###         'Q' => 0,
###         'E' => 0,
###         'G' => 0,
###         'H' => 0,
###         'I' => 0,
###         'L' => 0,
###         'K' => 0,
###         'M' => 0,
###         'F' => 0,
###         'P' => 0,
###         'S' => 0,
###         'T' => 0,
###         'W' => 0,
###         'Y' => 0,
###         'V' => 0,
###         '-' => 0,
###      } ;
##
##      my $freq ;
##
##      my $numdoms = keys %{$data->{aln}} ;
##      my @tp ;
##      foreach my $dom (keys %{$data->{aln}}) {
##         push @tp, substr($data->{aln}->{$dom}, $pos, 1) ;
##         $freq->{substr($data->{aln}->{$dom}, $pos, 1)}++ ; }
##
##      $data->{meta}->{numrestypes}->{$pos} = keys %{$freq} ;
##      $data->{meta}->{shannone}->{$pos} = 0 ;
##      foreach my $res (keys %{$freq}) {
##         $freq->{$res} /= $numdoms ;
##         $data->{meta}->{shannone}->{$pos} += $freq->{$res} *
##                                       (log($freq->{$res}) / log(21)) ;
##      }
###      if ($data->{shannone}->{$pos} == 0) {
###         print STDERR " shannone pos $pos $in->{aln_fn} = 0 charc:".join(",", @tp)."\n" ; }
##
##   }
##
##   return $data ;
##
##}
##
##
##sub parse_aln_raf {
##
##   my $in = shift;
##   my $data = $in->{alndata} ;
##   my $raf = $in->{raf} ;
##
###data->{seq} has raw FASTA sequence strings - check substr(curpos) to see
### if x is capital = frag or lower case = residue unk
##
##   my ($alnpos2atomresno, $atomresno2alnpos);
##   foreach my $dom (keys %{$data->{defstring}}) {
##
##      if (DEBUGALN) {
##         print STDERR "now on $dom\n" ;
##         print STDERR "   defstring:\t$data->{defstring}->{$dom}\n" ;
##         print STDERR "   alnstring:\t$data->{aln}->{$dom}\n" ;
##         print STDERR "   pdb:\t$data->{pdb}->{$dom}\n" ;
##         print STDERR "   class:\t$data->{class}->{$dom}\n" ;
##      }
##
##      my $gs2atomresno ;
##      my $gs2resna ; #to verify correct positioning when parsing alignment
##
###get gsresno -> atomresno mapping from raf file
##      my $gsresno = -1 ;
##      foreach my $fragdef ( @{$data->{frags}->{$dom}}) {
##         $gsresno++ ; #to take care of X in between fragments
##         my $tch = $fragdef->{chain} ; if ($tch eq ' ') {$tch = '_';}
##         my $line = $raf->{$data->{pdb}->{$dom}.$tch} ;
##         my $rafinfo = parse_raf_line({
##            line => $line,
##            headlength => $raf->{header_length}
##         }) ;
##
##         my ($ind_b, $ind_e) ;
##         if (exists $fragdef->{b}) {
##            if (defined $rafinfo->{atomresno2ind}->{$fragdef->{b}}) {
##               $ind_b = $rafinfo->{atomresno2ind}->{$fragdef->{b}} ;
##            } else {
##               if (DEBUGALN) {print STDERR " on dom $dom: lookahead b find for $fragdef->{b}\n" ;}
##               $ind_b = $rafinfo->{atomresno2ind_ahead}->{$fragdef->{b}} ;
##            }
##
##            if (defined $rafinfo->{atomresno2ind}->{$fragdef->{e}}) {
##               $ind_e = $rafinfo->{atomresno2ind}->{$fragdef->{e}} ;
##            } else {
##               if (DEBUGALN) {print STDERR " on dom $dom: lookback e find for $fragdef->{e}\n" ;}
##               $ind_e = $rafinfo->{atomresno2ind_back}->{$fragdef->{e}} ;
##            }
##
##         } else {
##
##
##            if (defined $rafinfo->{atomresno2ind}->{$rafinfo->{atomresno_first}}) {
##               $ind_b= $rafinfo->{atomresno2ind}->{$rafinfo->{atomresno_first}};
##            } else {
##               if (DEBUGALN) {print STDERR " on dom $dom: lookahead b find for $rafinfo->{atomresno_first}\n" ;}
##               $ind_b= $rafinfo->{atomresno2ind_ahead}->{$rafinfo->{atomresno_first}} ;
##            }
##
##            if (defined $rafinfo->{atomresno2ind}->{$rafinfo->{atomresno_last}}) {
##               $ind_e = $rafinfo->{atomresno2ind}->{$rafinfo->{atomresno_last}};
##            } else {
##               if (DEBUGALN) {print STDERR " on dom $dom: lookback e find for $rafinfo->{atomresno_last}\n" ;}
##               $ind_e = $rafinfo->{atomresno2ind_back}->{$rafinfo->{atomresno_last}} ;
##            }
##         }
##
##         my $keepon = 1 ;
##         my $pos = $ind_b ;
##         if (!defined $pos) { print STDERR " on dom $dom (chain $fragdef->{chain}), pos is unefined\n" ;}
##         if (!defined $ind_e) { print STDERR " on dom $dom (chain $fragdef->{chain}), ind_e is udnefinedpos\n";}
##
##         while ($pos <= $ind_e) {
##            if (exists $rafinfo->{ind2seqresno}->{$pos}) {
##               if (exists $rafinfo->{ind2atomresno}->{$pos}) {
##                  $gs2resna->{$gsresno} = uc($rafinfo->{seqresna}->[$pos]) ;
##                  $gs2atomresno->{$gsresno} = $rafinfo->{ind2atomresno}->{$pos}."\n".$fragdef->{chain};
##               }
##               $gsresno++ ;
##            }
##
##            $pos++ ;
##         }
##      }
##
###go through aln entry and assign aln pos -> atomresno mapping
##
##      my $curseqpos = 0 ;
##      foreach my $pos (0 .. (length($data->{aln}->{$dom}) - 1)) {
##         my $alnchar = substr($data->{aln}->{$dom}, $pos, 1) ;
##
###         if (DEBUGALN && $dom eq DEBUGDOM) {
###            print STDERR "dom $dom (alnpos $pos, curseqpos $curseqpos) $alnchar\n" ; }
##
##         if ($alnchar eq 'X' &&
##            (substr($data->{seq}->{$dom}, ($curseqpos), 1) eq 'X')) {
##
##            if (DEBUGALN) {
##               print STDERR "$dom: real frag $alnchar at seqpos $curseqpos (alnpos $pos)\n" ;}
##
##            $curseqpos++ ;
##
##         } elsif ($alnchar ne '-') {
##            if (DEBUGALN) {
##               if ($alnchar eq 'X') {print STDERR "$dom caught unk res: alnchar $alnchar at seq pos $curseqpos (alnpos $pos)\n";}
##            }
##
###            if ($alnchar eq 'X') {print STDERR "$dom caught unk res: alnchar $alnchar at seq pos $curseqpos (alnpos $pos) TRYING A FIX\n";next;}
##
##            if (exists $gs2atomresno->{$curseqpos}) {
##
##               if ($alnchar ne $gs2resna->{$curseqpos}) {
###                  print "ERROR $dom: alignment position $pos ($alnchar) mismatched with gs sequence position $curseqpos ($gs2resna->{$curseqpos}\n" ;
##                  print STDERR "ERROR $dom: alignment position $pos ($alnchar) mismatched with gs sequence position $curseqpos ($gs2resna->{$curseqpos}\n" ;
##               }
##
###               if (DEBUGALN && $dom eq DEBUGDOM) {
###                     print STDERR "$dom\t$pos\t$alnchar\t$gs2atomresno->{$curseqpos}\n" ;}
##
##               $alnpos2atomresno->{$dom}->{$pos} = $gs2atomresno->{$curseqpos} ;
##               $atomresno2alnpos->{$dom}->{$gs2atomresno->{$curseqpos}} = $pos ;
##            }
##            $curseqpos++ ;
##         }
##      }
##   }
##
##
##   return {
##     pos2resno => $alnpos2atomresno,
##     resno2pos => $atomresno2alnpos
##   } ;
##
##}
##
##sub load_interface_contacts {
##
##   my $in = shift ;
##   my $sid1 = $in->{sid1} ;
##   my $sid2 = $in->{sid2} ;
##   my $fn = $in->{fn} ;
##
##   if (!exists $in->{dist_thresh}) {
##      $in->{dist_thresh} = PARAM_INTERFACE_DIST_THRESH ; }
##
##   my ( $subset_id_1, $subset_id_2,
##        $chain_id_1, $resno_1, $chain_id_2, $resno_2, $min_dist ) =
##       pibase::rawselect_metatod($fn, "SELECT subset_id_1, subset_id_2, chain_id_1, resno_1, chain_id_2, resno_2, min_dist FROM $fn") ;
##
##   my $contacts ;
##   my $intres ;
##   foreach my $j ( 0 .. $#{$subset_id_1}) {
##      if ($subset_id_1->[$j] ne $sid1) {next;}
##      if ($subset_id_2->[$j] ne $sid2) {next;}
##      if ($min_dist->[$j] > $in->{dist_thresh}) {next;}
##
###      print STDERR "$subset_id_1->[$j] -- $subset_id_2->[$j]: $resno_1->[$j] $chain_id_1->[$j] to $resno_2->[$j] $chain_id_2->[$j]\n" ;
##
##      $intres->{$sid1}->{$resno_1->[$j]."\n".$chain_id_1->[$j]}++ ;
##      $intres->{$sid2}->{$resno_2->[$j]."\n".$chain_id_2->[$j]}++ ;
##      $contacts->{$resno_1->[$j]."\n".$chain_id_1->[$j]}->{$resno_2->[$j]."\n".$chain_id_2->[$j]}++ ;
##   }
##
##   return {
##      intres => $intres,
##      contacts => $contacts
##   } ;
##
##}
##
##
##sub load_astral_headers {
##
##   my $in = shift;
##   my $gdseq_fn = $in->{fn} ;
##   open(GDSEQF, $gdseq_fn) ;
##   my $gdseqh ;
##   my $gdom_fl ;
##   while (my $line = <GDSEQF>) {
##      if ($line =~ /^>/) {
##         $line =~ s/^>// ;
##         my ($t_dom, $t_class, $t_def, undef) = split(/\s/, $line) ;
##         $t_def =~ s/^\(// ; $t_def =~ s/\)$// ;
##         $gdseqh->{defstring}->{$t_dom} = $t_def ;
##         $gdseqh->{class}->{$t_dom} = $t_class;
##         $gdseqh->{pdb}->{$t_dom} = substr($t_dom, 1, 4) ;
##         if ($t_dom =~ /^g/) {
##            my $t_ddom = $t_dom ; $t_ddom =~ s/^g/d/ ;
##            $gdom_fl->{$t_ddom} = $t_dom ;
##         }
##
##         {
##            my @t = split(/\,/, $t_def) ;
##            my (@ch, @b, @e) ;
##            foreach my $j ( 0 .. $#t) {
##               my $t_frag ;
##               $t_frag->{chain} = '_' ;
##               if ($t[$j] =~ /:/) {
##                  $t_frag->{chain} = substr($t[$j],0,1) ; }
##               $t[$j] =~ s/^.\:// ;
##               my ($b, $e) = (' ', ' ');
##               if ($t[$j] =~ /.+\-.+/) {
##                  ($t_frag->{b}, $t_frag->{e}) =
##                  ($t[$j] =~ /(.+)\-(.+)/) ; }
##               push @{$gdseqh->{frags}->{$t_dom}}, $t_frag ;
##            }
##         }
##      }
##   }
##   close(GDSEQF) ;
##
##   return {
##      gdseqh => $gdseqh,
##      gdom => $gdom_fl
##   } ;
##
##}
##
##
##sub ligbase_preload {
##   my $in = shift ;
##   my $dbh = $in->{dbh} ;
##
##   my $query="SELECT pdb_code, a.ligand_code, a.ligand_idnum, a.ligand_chain, ".
##         "a.residue_num, a.ligand_chain FROM active_residues as a";
##
##   my ($pdb, $ligcod, $ligid, $ligchain, $resno, $pdbchain) =
##      pibase::mysql_fetchcols($dbh,$query) ;
##
##   my ($pdb2ligid, $pdb2res2ligid) ;
##   foreach my $j ( 0 .. $#{$ligcod}) {
##      if ($pdbchain->[$j] eq '') {
##         $pdbchain->[$j] = ' '; }
##
##      $pdb2ligid->{$pdb->[$j]}->{$ligcod->[$j]."\n".$ligid->[$j]}++ ;
##      $pdb2res2ligid->{$pdb->[$j]}->{$resno->[$j]."\n".$pdbchain->[$j]}->{$ligcod->[$j]."\n".$ligid->[$j]}++ ;
##   }
##
##   my $lb = {
##      pdb2ligid => $pdb2ligid,
##      pdb2res2ligid => $pdb2res2ligid
##   } ;
##
##   return $lb ;
##}
##
##
##sub tod_pibase_preload {
##
##   my $in = shift ;
##   my $astral = $in->{astral} ;
##
##   my ($bdp_2_raw, $pdb2bdp, $bdp2pdb) = pibase::todload_bdp_ids(
##      "bdp_id_2_raw_pdb", "pdb_id_2_bdp_id", "bdp_id_2_pdb_id") ;
##
##   my $bdp2contactsfn ;
##   {
##      my ($t_bdp, $t_fn) = pibase::rawselect_tod(
##         "SELECT bdp_id, source_file FROM interface_contacts_tables") ;
##      foreach my $j ( 0 .. $#{$t_bdp}) {
##         $bdp2contactsfn->{$t_bdp->[$j]} = $t_fn->[$j] ; }
##   }
##
##   my $bdp2subsetsresfn ;
##   {
##      my ($t_bdp, $t_fn) = pibase::rawselect_tod(
##         "SELECT bdp_id, source_file FROM subsets_residues_tables") ;
##      foreach my $j ( 0 .. $#{$t_bdp}) {
##         $bdp2subsetsresfn->{$t_bdp->[$j]} = $t_fn->[$j] ; }
##   }
##
##   my $sid2class ;
##   my $osid2class ;
##   my $sid2osid ;
##   my $sid2pdb ;
##   {
##      my ($sid, $class, $bdp_id, $ssid) =pibase::rawselect_tod(
##         "SELECT subset_id, class, bdp_id, subset_source_id FROM subsets") ;
##      foreach my $j ( 0 .. $#{$sid}) {
##         if ($ssid->[$j] != 1) {next;}
##         if (!defined $bdp_id->[$j] || $bdp_id->[$j] eq '') {next;}
##         $sid2pdb->{$sid->[$j]} = $bdp2pdb->{$bdp_id->[$j]} ;
##
##         $sid2class->{fam}->{$sid->[$j]} = $class->[$j] ;
##         my $fam = $class->[$j] ;
##         my ($sf) = ($fam =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
##         $sid2class->{sf}->{$sid->[$j]} = $sf ;
##
##         my $osid = substr($sid->[$j], -7, 7) ;
##         if (exists $astral->{gdom}->{$osid}) {
##            $osid = $astral->{gdom}->{$osid} ; }
##         $sid2osid->{$sid->[$j]} = $osid ;
##         $osid2class->{fam}->{$osid} = $class->[$j] ;
##
##      }
##   }
##
##
##   my ($bdp_id, $sid1, $sid2, $class1, $class2, $abchains) ; 
##   {
##      my ($t_bdp, $t_sid1, $t_sid2, $t_class1, $t_class2, $t_numcon, $t_chains)=
##      pibase::rawselect_tod(
##         "SELECT bdp_id, subset_id_1, subset_id_2, class_1, class_2, ".
##         " num_contacts, chains FROM intersubset_contacts") ;
##
##      foreach my $j ( 0 .. $#{$t_bdp}) {
##         if ($t_sid1->[$j] !~ /SCOP/) {next;}
##         if ($t_numcon->[$j] < PARAM_MIN_NUMCONTACTS) {next;}
##         if ($bdp_2_raw->{$t_bdp->[$j]} != 1) {next;}
##         push @{$bdp_id}, $t_bdp->[$j] ;
##         push @{$sid1}, $t_sid1->[$j] ;
##         push @{$sid2}, $t_sid2->[$j] ;
##         push @{$class1}, $t_class1->[$j] ;
##         push @{$class2}, $t_class2->[$j] ;
##         push @{$abchains}, $t_chains->[$j] ;
##      }
##   }
##
##   my $bdp2sid12 ;
##   my ($interacting_classes, $class2sidpairs, $sid12chain) ;
##   my ($osid1, $osid2, $revfl) ; #revfl if fam2 <fam1 and flipped to get fam21 isntead of fam12, let it be known so later it can be dealt with
##
##   foreach my $j ( 0 .. $#{$bdp_id}) {
### ONLY ONE OF THEM HAS TO BE OF PROPER CLASS
###      if ($class1->[$j] !~ /^[a-g]/ || $class2->[$j] !~ /^[a-g]/) {next;}
##      if ($class1->[$j] !~ /^[a-g]/ && $class2->[$j] !~ /^[a-g]/) {next;}
##      $bdp2sid12->{$bdp_id->[$j]}->{$sid1->[$j]."\t".$sid2->[$j]}++ ;
##
##      my $sid12 = $sid1->[$j]."\t".$sid2->[$j] ;
##      $sid12chain->{$sid12} = $abchains->[$j] ;
##
##      my $fam1 = $class1->[$j];
##      my $fam2 = $class2->[$j] ;
##
##      my ($sf1) = ($fam1 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
##      my ($sf2) = ($fam2 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
##
##      $osid1->[$j] = substr($sid1->[$j], -7, 7) ;
##      $osid2->[$j] = substr($sid2->[$j], -7, 7) ;
##
##      if (exists $astral->{gdom}->{$osid1->[$j]}) {
##         $osid1->[$j] = $astral->{gdom}->{$osid1->[$j]} ; }
##
##      if (exists $astral->{gdom}->{$osid2->[$j]}) {
##         $osid2->[$j] = $astral->{gdom}->{$osid2->[$j]} ; }
##
##      my $fam12;
##      if ($class2->[$j] lt $class1->[$j]) {
##         $revfl->[$j] = 1 ;
##         $fam12 = $class2->[$j]."\t".$class1->[$j];
##      } else {
##         $revfl->[$j] = 0 ;
##         $fam12 = $class1->[$j]."\t".$class2->[$j];
##      }
##      $interacting_classes->{fam}->{$fam1}->{$sid1->[$j]}->{$sid2->[$j]} = $j;
##      $interacting_classes->{fam}->{$fam2}->{$sid2->[$j]}->{$sid1->[$j]} = $j;
##      $interacting_classes->{sf}->{$sf1}->{$sid1->[$j]}->{$sid2->[$j]} = $j;
##      $interacting_classes->{sf}->{$sf2}->{$sid2->[$j]}->{$sid1->[$j]} = $j;
##
##      $class2sidpairs->{fam}->{$fam1}->{$j}++ ;
##      $class2sidpairs->{fam}->{$fam2}->{$j}++ ;
##      $class2sidpairs->{sf}->{$sf1}->{$j}++ ;
##      $class2sidpairs->{sf}->{$sf2}->{$j}++ ;
##   }
##
##
##   my $chain_2_pdbchain ;
##   my ($chain_2_start, $chain_2_end) ;
##   my $pdbchains ;
##   {
##      my ($tbdp_id, $real_chain_id, $pdb_chain_id, $resno_start, $resno_end) ;
##      {
##         my ($t_bdp, $t_real_cid, $t_pdb_cid, $t_start, $t_end, $t_type) =
##         pibase::rawselect_tod(
##            "SELECT bdp_id, real_chain_id, pdb_chain_id, ".
##            "start_resno, end_resno, chain_type FROM bdp_chains") ;
##    
##         foreach my $j ( 0 .. $#{$t_bdp}) {
##            if ($bdp_2_raw->{$t_bdp->[$j]} != 1) {next;}
##            if (!defined $t_pdb_cid->[$j] || $t_pdb_cid->[$j] eq '') {next;}
##            if ($t_type->[$j] ne 'p') {next;}
##            push @{$tbdp_id}, $t_bdp->[$j] ;
##            push @{$real_chain_id}, $t_real_cid->[$j] ;
##            push @{$pdb_chain_id}, $t_pdb_cid->[$j] ;
##            push @{$resno_start}, $t_start->[$j] ;
##            push @{$resno_end}, $t_end->[$j] ;
##         }
##      }
##
###      print STDERR ($#{$tbdp_id} + 1)." chains loaded entries\n " ;
##
##      foreach my $j ( 0 .. $#{$tbdp_id}) {
##         $pdbchains->{$bdp2pdb->{$tbdp_id->[$j]}}->{$real_chain_id->[$j]}++;
##
##         $chain_2_pdbchain->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
##            $pdb_chain_id->[$j] ;
##         $chain_2_start->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
##            $resno_start->[$j] ;
##         $chain_2_end->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
##            $resno_end->[$j] ;
##      }
##   }
##
##   my $pb = {
##      bdp2contactsfn => $bdp2contactsfn,
##      bdp2subsetsresfn => $bdp2subsetsresfn,
##      bdp2pdb => $bdp2pdb,
##      pdb2bdp => $pdb2bdp,
##      bdp_id => $bdp_id,
##      sid2pdb => $sid2pdb,
##      sid2class => $sid2class,
##      osid2class => $osid2class,
##      sid2osid => $sid2osid,
##      sid1 => $sid1,
##      sid2 => $sid2,
##      class1 => $class1,
##      class2 => $class2,
##      osid1 => $osid1,
##      osid2 => $osid2,
##      revfl => $revfl,
##      interacting_classes => $interacting_classes,
##      class2sidpairs => $class2sidpairs,
##      pdbchains => $pdbchains,
##      chain_2_pdbchain => $chain_2_pdbchain,
##      chain_2_start => $chain_2_start,
##      chain_2_end => $chain_2_end,
##      bdp2sid12 => $bdp2sid12,
##      sid12chain => $sid12chain
##   } ;
##
##   return $pb ;
##}
##
##
##sub pibase_preload {
##
##   my $in = shift ;
##   my $astral = $in->{astral} ;
##
##   my $dbh = $in->{dbh} ;
##
##   my $bdp2contactsfn = pibase::mysql_hashload($dbh,
##      "SELECT bdp_id, source_file FROM interface_contacts_tables") ;
##
##   my $bdp2subsetsresfn = pibase::mysql_hashload($dbh,
##      "SELECT bdp_id, source_file FROM subsets_residues_tables") ;
##
##   my $sid2class ;
##   my $osid2class ;
##   $sid2class->{fam} = pibase::mysql_hashload($dbh,
##      "SELECT subset_id, class FROM subsets WHERE bdp_id IS NOT NULL
##      and subset_source_id = 1") ;
##
##   my $sid2osid ;
##   foreach my $sid (keys %{$sid2class->{fam}}) {
##      my $fam = $sid2class->{fam}->{$sid} ;
##      my $osid = substr($sid, -7, 7) ;
##      if (exists $astral->{gdom}->{$sid}) {
##         $osid = $astral->{gdom}->{$sid} ; }
##      $sid2osid->{$sid} = $osid ;
##      my ($sf) = ($fam =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
##      $sid2class->{sf}->{$sid} = $sf ;
##      $osid2class->{sf}->{$osid} = $sf ;
##      $osid2class->{fam}->{$osid} = $fam ;
##   }
##
##   my ($bdp_id, $sid1, $sid2, $class1, $class2) = 
##      pibase::mysql_fetchcols($dbh,
##      "SELECT a.bdp_id, subset_id_1, subset_id_2, class_1, class_2 ".
##      "FROM intersubset_contacts as a, bdp_files as b ".
##      "WHERE subset_id_1 LIKE \"\%SCOP\%\" AND a.bdp_id = b.bdp_id ".
##      "AND num_contacts >= ".PARAM_MIN_NUMCONTACTS." ".
##      "AND b.raw_pdb = 1") ;
##
##   my $bdp2sid12 ;
##   my ($interacting_classes, $class2sidpairs) ;
##   my ($osid1, $osid2, $revfl) ; #revfl if fam2 <fam1 and flipped to get fam21 isntead of fam12, let it be known so later it can be dealt with
##
##   foreach my $j ( 0 .. $#{$bdp_id}) {
### ONLY ONE OF THEM HAS TO BE OF PROPER CLASS
###      if ($class1->[$j] !~ /^[a-g]/ || $class2->[$j] !~ /^[a-g]/) {next;}
##      if ($class1->[$j] !~ /^[a-g]/ && $class2->[$j] !~ /^[a-g]/) {next;}
##      $bdp2sid12->{$bdp_id->[$j]}->{$sid1->[$j]."\t".$sid2->[$j]}++ ;
##
##      my $fam1 = $class1->[$j];
##      my $fam2 = $class2->[$j] ;
##
##      my ($sf1) = ($fam1 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
##      my ($sf2) = ($fam2 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
##
##      $osid1->[$j] = substr($sid1->[$j], -7, 7) ;
##      $osid2->[$j] = substr($sid2->[$j], -7, 7) ;
##
##      if (exists $astral->{gdom}->{$osid1->[$j]}) {
##         $osid1->[$j] = $astral->{gdom}->{$osid1->[$j]} ; }
##
##      if (exists $astral->{gdom}->{$osid2->[$j]}) {
##         $osid2->[$j] = $astral->{gdom}->{$osid2->[$j]} ; }
##
##      my $fam12;
##      if ($class2->[$j] lt $class1->[$j]) {
##         $revfl->[$j] = 1 ;
##         $fam12 = $class2->[$j]."\t".$class1->[$j];
##      } else {
##         $revfl->[$j] = 0 ;
##         $fam12 = $class1->[$j]."\t".$class2->[$j];
##      }
##      $interacting_classes->{fam}->{$fam1}->{$sid1->[$j]}->{$sid2->[$j]} = $j;
##      $interacting_classes->{fam}->{$fam2}->{$sid2->[$j]}->{$sid1->[$j]} = $j;
##      $interacting_classes->{sf}->{$sf1}->{$sid1->[$j]}->{$sid2->[$j]} = $j;
##      $interacting_classes->{sf}->{$sf2}->{$sid2->[$j]}->{$sid1->[$j]} = $j;
##
##      $class2sidpairs->{fam}->{$fam1}->{$j}++ ;
##      $class2sidpairs->{fam}->{$fam2}->{$j}++ ;
##      $class2sidpairs->{sf}->{$sf1}->{$j}++ ;
##      $class2sidpairs->{sf}->{$sf2}->{$j}++ ;
##   }
##
##   my $pdb2bdp = pibase::mysql_hashload($dbh,
##      "SELECT pdb_id, bdp_id FROM bdp_files WHERE raw_pdb =1 ") ;
##
##   my $bdp2pdb = pibase::mysql_hashload($dbh,
##      "SELECT bdp_id, pdb_id FROM bdp_files") ;
##
##   my $chain_2_pdbchain ;
##   my ($chain_2_start, $chain_2_end) ;
##   my ($chain_2_startser, $chain_2_endser) ;
##   my $pdbchains ;
##   {
##      my ($tbdp_id, $real_chain_id, $pdb_chain_id, $resno_start,
##         $resno_start_serial, $resno_end, $resno_end_serial) =
##      pibase::mysql_fetchcols($dbh,
##      "SELECT a.bdp_id, real_chain_id, pdb_chain_id, ".
##      "start_resno, start_resno_int, ".
##      "end_resno, end_resno_int FROM bdp_chains as a, bdp_files as b ".
##      "WHERE pdb_chain_id IS NOT NULL AND chain_type = \"p\" AND ".
##      "a.bdp_id = b.bdp_id AND raw_pdb = 1");
##
###      print STDERR ($#{$tbdp_id} + 1)." chains loaded entries\n " ;
##
##      foreach my $j ( 0 .. $#{$tbdp_id}) {
##         $pdbchains->{$bdp2pdb->{$tbdp_id->[$j]}}->{$real_chain_id->[$j]}++;
##
##         $chain_2_pdbchain->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
##            $pdb_chain_id->[$j] ;
##         $chain_2_start->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
##            $resno_start->[$j] ;
##         $chain_2_end->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
##            $resno_end->[$j] ;
##         $chain_2_startser->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
##            $resno_start_serial->[$j] ;
##         $chain_2_endser->{$tbdp_id->[$j]}->{$real_chain_id->[$j]} =
##            $resno_end_serial->[$j] ;
##      }
##   }
##
##   my $pb = {
##      bdp2contactsfn => $bdp2contactsfn,
##      bdp2subsetsresfn => $bdp2subsetsresfn,
##      bdp2pdb => $bdp2pdb,
##      pdb2bdp => $pdb2bdp,
##      bdp_id => $bdp_id,
##      sid2class => $sid2class,
##      sid2osid => $sid2osid,
##      osid2class => $osid2class,
##      sid1 => $sid1,
##      sid2 => $sid2,
##      class1 => $class1,
##      class2 => $class2,
##      osid1 => $osid1,
##      osid2 => $osid2,
##      revfl => $revfl,
##      interacting_classes => $interacting_classes,
##      class2sidpairs => $class2sidpairs,
##      pdbchains => $pdbchains,
##      chain_2_pdbchain => $chain_2_pdbchain,
##      chain_2_start => $chain_2_start,
##      chain_2_end => $chain_2_end,
##      chain_2_startser => $chain_2_startser,
##      chain_2_endser => $chain_2_endser,
##      bdp2sid12 => $bdp2sid12
##   } ;
##
##   return $pb ;
##}
##
##sub get_astral_classlist {
##
##   my $in = shift;
##
##   my $classlist ;
##   foreach my $type (qw/fam sf/) {
##      my @files = glob($in->{fn}->{aln}->{$type}."/*.aln.fa") ;
##      for (@files) {
##         my $a = $_ ;
##         $a =~ s/\.aln\.fa$// ; $a =~ s/.*\/// ;
##         $classlist->{$type}->{$a}++ ;
##       }
##   }
##
##   return $classlist ;
##  
##}
##
##
##sub load_astral_clusters {
##
##   my $in = shift ;
##   my $out = $in->{out} ;
##   my $fn = $in->{fn} ;
##
##   my $seqcl ;
##   my $seqcl2cont ;
##   my $clusnum2rep;
##   my $rep2clusnum;
##   print STDERR "load ASTRAL sequence clusters: " ;
##   {
##      foreach my $seqid (keys %{$fn->{seqcl}}) {
##         open (ASTRALCL, $fn->{seqcl}->{$seqid}) ;
##         my $clusnum = 0 ;
##         while (my $line = <ASTRALCL>) {
##            chomp $line;
##            if ($line =~ /^Rep/) {
##               $line =~ s/^Rep: //g;
##               $clusnum++ ;
##               $line =~ s/^\s+// ;
##               my $scopid = substr($line, 0, 7) ;
##               $clusnum2rep->{$clusnum} = $scopid ;
##               $rep2clusnum->{$scopid} = $clusnum;
##            }
##
##            if ($line =~ /score/) {
##               $line =~ s/^\s+// ;
##               my $scopid = substr($line, 0, 7) ;
##               $seqcl->{$seqid}->{$scopid} = $clusnum2rep->{$clusnum} ;
##               push @{$seqcl2cont->{$seqid}->{$clusnum2rep->{$clusnum}}},
##                  $scopid ;
##            }
##         }
##         close(ASTRALCL) ;
##      }
##   }
##   print STDERR "X\n" ;
##
##   $out->{seqcl} = $seqcl ;
##   $out->{clusnum2rep} = $clusnum2rep ;
##   $out->{seqcl2cont} = $seqcl2cont ;
##
##}
##
##
##
### old code snippets
###
### 1. supplementing the asteroids alignment with the sequences from 100%
###    identical domains that arent included.
### location: in load_asteroids_aln() right after read_asteroids_aln() returns.
###
###no longer necessary - aln are built locally with all domains
###   #supplement alignment with the 100% sequence matches
###   my $added ;
###   my @origdoms = keys %{$data->{aln}} ;
###   foreach my $rep_scopid (@origdoms) {
###      if (!exists $seqcl100->{$rep_scopid}) {next;}
###      my $clno = $seqcl100->{$rep_scopid} ;
###      if ($#{$seqclcont100->{$clno}} > 0) {
###         foreach my $j ( 1 .. $#{$seqclcont100->{$clno}}) {
###            my $scopid = $seqclcont100->{$clno}->[$j] ;
###
###            if (exists $added->{$scopid}) {next;}
###
###            $data->{aln}->{$scopid} = $data->{aln}->{$rep_scopid} ;
####            print STDERR "100SUPPed $scopid ($rep_scopid): $data->{aln}->{$scopid}\n" ;
###
###            $data->{defstring}->{$scopid} = $gdseqh->{defstring}->{$scopid} ;
###            $data->{class}->{$scopid} = $gdseqh->{class}->{$scopid} ;
###            $data->{pdb}->{$scopid} = $gdseqh->{pdb}->{$scopid} ;
###            $data->{frags}->{$scopid} = $gdseqh->{frags}->{$scopid} ;
###         }
###      }
###   }
##
##
###sub pibase_preload {
###
###   my $dbh = shift ;
###   my $data ;
###
###   print STDERR "subsres hashload\n" ;
###   $data->{tn}->{subsets_res} = pibase::mysql_hashload($dbh->{pi},
###      "SELECT bdp_id, source_file FROM subsets_residues_tables") ;
###
###   print STDERR "interface contacts\n" ;
###   ($data->{tn}->{interfacecon_res}) = pibase::mysql_hashload($dbh->{pi},
###      "SELECT bdp_id, source_file FROM interface_contacts_tables") ;
###
###   print STDERR "pdb->bdp\n" ;
###   ($data->{pdb_2_bdp}) = pibase::load_bdp_ids($dbh->{pi}, qw/pdb_id_2_bdp_id/) ;
###
###   ($data->{scop_fam}) = pibase::mysql_fetchcols($dbh->{pi},
###      "SELECT distinct class FROM subsets where subset_source_id = 1") ;
###
###   ($data->{sid2scop}) = pibase::mysql_hashload($dbh->{pi},
###         "SELECT subset_id, class FROM subsets WHERE subset_source_id = 1") ;
###
###   return $data ;
###}
##
##
##sub get_subsres {
##   
##   my $fn = shift ;
##      
##   my $assign ;
##   ($assign->{subset_id}, $assign->{chain_id}, $assign->{resno_serial},
##    $assign->{resno}) = pibase::rawselect_metatod($fn,
##      "SELECT subset_id, chain_id, resno_serial, resno FROM $fn") ;
##      
##   my $res2sub ;
##   foreach my $j ( 0 .. $#{$assign->{subset_id}}) {
##      if ($assign->{subset_id}->[$j] =~ /SCOP/) {
##         my $sig = $assign->{resno}->[$j]."\n".$assign->{chain_id}->[$j] ;
##         $res2sub->{$sig} = $assign->{subset_id}->[$j] ; }
##   }
##
##   return $res2sub ;
##   
##}
##
##
##sub get_ligbase_domains {
##   my $thissub = 'get_ligbase_domains' ;
##
##   my $in = shift ;
##
##   my $pb = $in->{pb} ;
##   my $lb = $in->{lb} ;
##
##   my ($ligres2sid, $sid2ligres, $class2ligsid) ;
##
##   my $count = 0 ;
##   foreach my $pdb_id (keys %{$lb->{pdb2res2ligid}}) {
##      print STDERR "PROGRESS: $thissub now on $pdb_id\n" ;
###      print STDERR "DEBUG $thissub ".__LINE__." NOW assigning to $pdb_id\n " if DEBUG;
##      if (!exists $pb->{pdb2bdp}->{$pdb_id}) {
##         print STDERR "ERROR $thissub ".__LINE__.": $pdb_id no bdp entry\n"; next ; }
##      my $bdp_id = $pb->{pdb2bdp}->{$pdb_id} ;
##
##      if (!exists $pb->{bdp2subsetsresfn}->{$bdp_id}) {
##         print STDERR "ERROR $thissub ".__LINE__.": $pdb_id ($bdp_id) no subsets entry\n";
##	 next ; }
##
##      if (!-s $pb->{bdp2subsetsresfn}->{$bdp_id}) {
##         print STDERR "ERROR $thissub ".__LINE__.": $pdb_id empty subets_res file $pb->{bdp2subsetsresfn}->{$bdp_id}\n" ;
##         next ;
##      }
##
##      my $curscopdef = get_subsres($pb->{bdp2subsetsresfn}->{$bdp_id}) ;
##
##      foreach my $ressig (keys %{$lb->{pdb2res2ligid}->{$pdb_id}}) {
##         if (exists $curscopdef->{$ressig}) {
##            my $sid = $curscopdef->{$ressig} ;
##	    $ligres2sid->{$pdb_id}->{$ressig}->{$sid}++ ;
##	    $sid2ligres->{$sid}->{$ressig}++ ;
##            $class2ligsid->{fam}->{$pb->{sid2class}->{fam}->{$sid}}->{$sid}++ ;
##            $class2ligsid->{sf}->{$pb->{sid2class}->{sf}->{$sid}}->{$sid}++ ;
###            print STDERR "FOUND: $pdb_id $ressig in domain $sid is ligand bound\n" if DEBUG;
##	 } else {
##            my $tres = $ressig ; $tres =~ s/\n/:/ ;
##            print STDERR "WARNING $thissub ".__LINE__.": $pdb_id no domain definition for ligbase residue $tres\n" ;
##         }
##      }
##      $count++ ;
##      if (DEBUG && $count == 100) {last;}
##   }
##
##   return {
##      ligres2sid => $ligres2sid,
##      sid2ligres => $sid2ligres,
##      class2ligsid => $class2ligsid
##   } ;
##}
##
##sub load_liginfo {
##
##   my $in = shift ;
##
##   my $liginfo ;
##   open(LIGINFOF, $in->{fn}->{msdchem}->{liginfo}) ;
##   while (my $line = <LIGINFOF>) {
##      if ($line =~ /^\#/) {next;}
##      chomp $line;
##
##      my @t = split(/\t/, $line) ;
##      my $lig = shift @t ;
##      $liginfo->{ligs}->{$lig}++ ;
##
##      ($liginfo->{name}->{$lig},
##       $liginfo->{formula}->{$lig},
##       $liginfo->{mw}->{$lig},
##       $liginfo->{numatoms}->{$lig},
##       $liginfo->{numatoms_nonh}->{$lig} ) = @t ;
##
##      if ($liginfo->{mw}->{$lig} eq '') {
##         delete $liginfo->{mw}->{$lig};
##      } elsif ($liginfo->{mw}->{$lig} >= PARAM_MIN_LIGMW &&
##               $liginfo->{mw}->{$lig} <= PARAM_MAX_LIGMW) {
##         $liginfo->{mwinrange}->{$lig}++ ;
##      }
##   }
##   close(LIGINFOF) ;
##
##   return $liginfo ;
##
##}
##
##sub get_sasa { # purpose: run and parse modeller sasa (psa)
##
##   my $in = shift ;
##
##   my ($bdp_file, $modeller_bin) ;
##   my $surftyp = '1' ; #default contact
##
##   if (ref($in) ne '') {
##      if (exists $in->{surftyp}) {
##         $surftyp = $in->{surftyp} ; }
##
##      $bdp_file = $in->{pdb_fn} ;
##      $modeller_bin = $in->{modeller_bin} ;
##   } else {
##      $bdp_file = $in ;
##      $modeller_bin = shift ;
##   }
##
##=pod
##
##Specify MODELLER temporary TOP file, and output files.
##
##=cut
##
##   my $filesmade ;
##
##   my ($sasa_top_fh, $sasa_top_fn) = 
##      tempfile("temp_sasa.XXXXXX", SUFFIX=>".top") ;
##   $filesmade->{$sasa_top_fn}++ ;
##
##   my $tempbase = $sasa_top_fn ; $tempbase =~ s/\.top$// ;
##   my $sasa_log_fn = $tempbase.".log" ;
##   my $sasa_atm_fn = $tempbase.".sol" ;
##   my $sasa_res_fn = $tempbase.".psa" ;
##
##   $filesmade->{$sasa_log_fn}++ ;
##   $filesmade->{$sasa_atm_fn}++ ;
##   $filesmade->{$sasa_res_fn}++ ;
##
##=pod
##
##Write the MODELLER TOP file to calcaulte solvent accessibilty.
##
##=cut
##
##   my $filename = $bdp_file ; $filename =~ s/^.*\///g ;
##   my $filebase = $bdp_file ; $filebase =~ s/\/[^\/]+$//g ;
##
##   if ($filename ne $bdp_file) {
##      print $sasa_top_fh "SET ATOM_FILES_DIRECTORY = \'$filebase\'\n" ; }
##
##   print $sasa_top_fh 'READ_TOPOLOGY FILE = \'$(LIB)/top_heav.lib\''."\n" ;
##   print $sasa_top_fh 'SET HETATM_IO = off, WATER_IO=off'."\n" ;
##   print $sasa_top_fh 'READ_MODEL FILE = \''.$filename.'\''."\n" ;
##   print $sasa_top_fh 'SET SURFTYP = '."$surftyp\n" ;
##   print $sasa_top_fh 'SET RADII_FACTOR = 1.0'."\n" ;
##   print $sasa_top_fh 'WRITE_DATA FILE = \''.$tempbase.'\', OUTPUT = \'PSA\''."\n" ;
##   close($sasa_top_fh) ;
##
##=pod
##
##Run the MODELLER TOP file.
##
##=cut
##
##   my $tcom = "$modeller_bin $sasa_top_fn >/dev/null 2>&1" ;
##   system($tcom) ;
##
##=pod
##
##Initialize solvent accessibility storage arrays
##
##=cut
##
##   my $results ;
##   my $resno_rev ;
##   my %fields = (
##      'resno' => [7, 4] ,
##      'resna' => [14, 3] ,
##      'chain' => [18, 1],
##      'all_sum' => [19, 7],
##      'all_perc' => [27, 5],
##      'nonp_sum' => [33, 7],
##      'nonp_perc' => [41, 5],
##      'p_sum' => [47, 7],
##      'p_perc' => [55, 5],
##      'sc_sum' => [61, 7],
##      'sc_perc' => [69, 5],
##      'mc_sum' => [75, 7],
##      'mc_perc' => [83, 5]
##   ) ;
##
##   foreach my $key (keys %fields) {
##      $results->{$key} = [] ; }
##
##=pod
##
##If the MODELLER PSA file does not exist, display to STDERR and return.
##
##=cut
##
##   if (!-s $sasa_res_fn) {
##      $results->{resno}->[0] = 'ERROR' ;
##      $resno_rev->{ERROR} = 'ERROR' ;
##      my $err = ['DANGER WILL ROBINSON'] ;
##      foreach my $tfn (keys %{$filesmade}) {if (-e $tfn) {unlink $tfn;}}
##      return ($results, $resno_rev, undef, $err);
##   }
##
##=pod
##
##Open the resulting MODELLER PSA file.
##
##=cut
##
##
##=pod
##
##Read in a line of the PSA file.
##
##=cut
##
##   my $sasaatoms ;
##   $sasaatoms->{all} = 0 ;
##   $sasaatoms->{p} = 0 ;
##   $sasaatoms->{nonp} = 0 ;
##
##   open(SOL, $sasa_atm_fn) ;
##   while (my $line = <SOL>) {
##      chomp $line;
##      if ($line !~ /^ATOM/) {next;}
##
##      my $atomtype = substr($line,12,1) ;
##      my $atomacc = substr($line,64,8) ;
##      $atomacc =~ s/ //g ;
##      $sasaatoms->{all} += $atomacc ;
##      if ($atomtype =~ /[NO]/) {
##         $sasaatoms->{p} += $atomacc ;
##      } else {
##         $sasaatoms->{nonp} += $atomacc ;
##      }
##   }
##   close(SOL) ;
##
##
##   my $sasa ;
##   $sasa->{all} = 0 ;
##   $sasa->{mc} = 0 ;
##   $sasa->{sc} = 0 ;
##   $sasa->{nonp} = 0 ;
##   $sasa->{p} = 0 ;
##   open(PSA, $sasa_res_fn) ;
##   while (my $line = <PSA>) {
##
##      chomp $line;
##
##=over
##
##If the line contains accessibility information,
##
##=cut
##
##      if ($line =~ /^ACCESS/) {
##
##=over
##
##Extract the residue information and solvent accessibility parameters from the line.
##
##=cut
##
##         foreach my $key (keys %fields) {
##
##            my $t = substr($line, $fields{$key}->[0], $fields{$key}->[1]);
##            unless ($key eq 'chain') {
##               $t =~ s/ //g ; }
##
##	    if ($key eq 'all_sum') {
##	       $sasa->{all} += $t;
##	    } elsif ($key eq 'mc_sum') {
##	       $sasa->{mc} += $t;
##	    } elsif ($key eq 'sc_sum') {
##	       $sasa->{sc} += $t;
##	    } elsif ($key eq 'nonp_sum') {
##	       $sasa->{nonp} += $t;
##	    } elsif ($key eq 'p_sum') {
##	       $sasa->{p} += $t;
##	    }
##
##            push @{$results->{$key}}, $t ;
##         }
##
##=pod
##
##Create an entry in the reverse lookup hash pointing from the actual PDB residue number and chain id to the index of the arrays where the solvent accessibility values are stored.
##
##=cut
##
##         my $cur_recno = $#{$results->{resno}} ;
##         my $res_sig = $results->{resno}->[$cur_recno]."_".
##                       $results->{chain}->[$cur_recno] ;
##
##         $resno_rev->{$res_sig} = $cur_recno ;
##
##      }
##
##=back
##
##=cut
##
##   }
##
##=back
##
##Close the PSA file.
##
##=cut
##
##   close PSA ;
##
##=pod
##
##If no accessibility lines have been read in, set error flag
##
##=cut
##
##   my $error_fl ;
##
##   if ($#{$results->{resno}} < 0 ) {
##      push @{$error_fl}, "no sasa entries calculated" ; }
##
##=pod
##
##Remove the MODELLER TOP file and output files.
##
##=cut
##
##   unlink($sasa_top_fn, $sasa_log_fn, $sasa_atm_fn, $sasa_res_fn) ;
##
##=pod
##
##Return pointers to teh arrays holding the values and the residue lookup hash.
##
##Note: accessibility values for ALL Cys read are retunred; However, the calling function will only ask for Cys involved in disulfide bridges.
##
##=cut
##
##   return ( $results, $resno_rev, $sasa, $error_fl, $sasaatoms ) ;
##
##}
##
##
##sub readin_ligassignments {
##
##   my $in = shift ;
##
##   my $fn = $in->{fn} ;
##
##   my $class2alnlength = $in->{class2alnlength};
##   my $standardres= $in->{standardres};
##   my $liginfo = $in->{liginfo};
##   my $ligbits ;
##   print STDERR "NOW reading LIG assignment\n" ;
##   open(LIGFH, $fn->{pilig}->{park_assign_lig_fn}) ;
##   while (my $line = <LIGFH>) {
##      if ($line =~ /^#/ ||
##          $line =~ /^Warning: no access to/ ||
##          $line =~ /^Thus no job control/ ) {next;}
##      chomp $line;
##
##      my ( $pdb, $sid, $osid, $classtype, $class, 
##              $btype, $alnposstring, $alnlength,
##              $ligcod, $ligid ) = split(/\t/, $line) ;
##      $ligcod =~ s/ //g ;
##      my $ligsig = $pdb."\t".$ligcod."\t".$ligid ;
##
##      if (exists $standardres->{$ligcod}) {next;}
##      if ($ligcod eq 'DUM' || $ligcod eq 'UNX' ||
##             $ligcod eq 'UNK' || $ligcod eq 'UNL') {next;}
##
##      if (!exists $liginfo->{mw}->{$ligcod}) {
##            print STDERR "WARNING: MW not found for ligand $ligcod\n" ;}
##
##      $class2alnlength->{$classtype}->{$class} = $alnlength ;
##
##      $ligbits->{$classtype}->{$sid}->{$ligsig} = Bit::Vector->new($alnlength) ;
##
##      my @alnpos = split(/\,/, $alnposstring) ;
##      foreach my $alnpos (@alnpos) {
##         $ligbits->{$classtype}->{$sid}->{$ligsig}->Bit_On($alnpos); }
##
##
##      if (!exists $ligbits->{$classtype}->{$sid}->{cumulative}) {
##         $ligbits->{$classtype}->{$sid}->{cumulative} =
##            $ligbits->{$classtype}->{$sid}->{$ligsig}->Clone();
##      } else {
##         $ligbits->{$classtype}->{$sid}->{cumulative}->Union(
##            $ligbits->{$classtype}->{$sid}->{cumulative},
##            $ligbits->{$classtype}->{$sid}->{$ligsig}) ;
##      }
###      if (exists $liginfo->{mwinrange}->{$ligcod})
##   }
##   close(LIGFH) ;
##
##   return $ligbits ;
##}
##
##
##sub readin_expassignments {
##
##   my $in = shift ;
##   my $fn = $in->{fn};
##   my $class2alnlength = $in->{class2alnlength};
##
##   print STDERR "NOW reading EXP assignment\n" ;
##   open(EXPFH, $fn->{pilig}->{park_assign_exp_fn}) ;
##   my $expbits ;
##   while (my $line = <EXPFH>) {
##      if ($line =~ /^#/ ||
##            $line =~ /^Warning: no access to/ ||
##            $line =~ /^Thus no job control/ ) {next;}
##      chomp $line;
##
##      my ( $pdb, $sid, $osid, $classtype,
##           $class, $btype, $alnposstring, $alnlength) = split(/\t/, $line) ;
##
##      $class2alnlength->{$classtype}->{$class} = $alnlength ;
##
##      if (!exists $expbits->{$classtype}->{$class}) {
##         $expbits->{$classtype}->{$class} =
##            Bit::Vector->new($alnlength) ;
##      }
##
##      my @alnpos = split(/\,/, $alnposstring) ;
##      foreach my $alnpos (@alnpos)  {
##         if ($alnpos eq "undef") {next;}
##         $expbits->{$classtype}->{$class}->Bit_On($alnpos);
##      }
##   }
##   close(EXPFH) ;
##
##   return $expbits ;
##
##}
##
##
##sub readin_piassignments {
##
##   my $in = shift;
##   my $fn = $in->{fn} ;
##   my $class2alnlength = $in->{class2alnlength};
##   my $pb = $in->{pb} ;
##   my $pibits ;
##   my $interfaces ;
##
##   print STDERR "NOW reading PI assignment\n" ;
##   open(PIFH, $fn->{pilig}->{park_assign_pi_fn}) ;
##   while (my $line = <PIFH>) {
##      if ($line =~ /^#/ ||
##         $line =~ /^Warning: no access to/ ||
##         $line =~ /^Thus no job control/ ) {next;}
##      chomp $line;
##
##      my ($pdb, $sid, $osid, $classtype, $class,
##          $obtype, $alnposstring, $alnlength,
##          $sid1, $sid2, $fam1, $fam2, $chains) = split(/\t/, $line) ;
##
##      my $sid12 = $sid1."\t".$sid2 ;
##      my $side = 1; if ($sid eq $sid2) {$side = 2;}
##
##      $interfaces->{$sid12}->{$side}->{sid} = $sid ;
##      $interfaces->{$sid12}->{chains} = $chains ;
##      $interfaces->{$sid12}->{pdb} = $pdb;
##
##      my @alnpos = split(/\,/, $alnposstring) ;
##      {
##            my @t = ();
##            foreach my $p ( @alnpos) {
##               if ($p ne 'undef') {push @t, $p;} }
##            @alnpos = @t ;
##      }
##      if ($#alnpos < 0) {next;}
##
##      $class2alnlength->{$classtype}->{$class} = $alnlength ;
##
##      my @btypes = ($obtype) ;
##      if ($chains eq 'same') {
##         push @btypes, "Pintra" ;
##      } else {
##         push @btypes, "Pinter" ;
##      }
##
##      $interfaces->{$sid12}->{$side}->{pibits}->{$classtype} =
##         Bit::Vector->new($alnlength) ;
##
##      foreach my $alnpos (@alnpos)  {
##         if ($alnpos eq 'undef') {next;}
##         $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Bit_On($alnpos); }
##
##      foreach my $btype (@btypes) {
##         if (!exists $pibits->{$classtype}->{$sid}->{$btype}) {
##            $pibits->{$classtype}->{$sid}->{$btype} =
##               $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Clone();
##         } else {
##            $pibits->{$classtype}->{$sid}->{$btype}->Union(
##               $pibits->{$classtype}->{$sid}->{$btype},
##               $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}) ;
##         }
##      }
##   }
##   close(PIFH) ;
##
##   return {
##      pibits => $pibits,
##      interfaces => $interfaces
##   } ;
##}
##
##
###      if ($classtype eq 'fam') {
###         foreach my $seqid (qw/90 95 100/) {
###            push @curclasstypes, 'seqcl'.$seqid ;
###            my $tosid = $osid;
###            if (substr($tosid,5,1) eq '.') { substr($tosid,0,1) = 'g' ; }
###            my $thiscl = $tosid ;
###            if (defined $astral->{seqcl}->{$seqid}->{$tosid}) {
###                  $thiscl = $astral->{seqcl}->{$seqid}->{$tosid} ;
###            } else {
###                  $astral->{seqcl}->{$seqid}->{$tosid} = $thiscl;
###                  $astral->{seqcl2cont}->{$seqid}->{$thiscl} = [$tosid];
###            }
###            push @curclasses, $thiscl ;
###         }
###      }
##
##
##sub BK061111_run_collate_perinstance {
##
##   require Bit::Vector ;
##
##
##   my $standardres = {
##   ALA => 'A' ,
##   ARG => 'R' ,
##   ASN => 'N' ,
##   ASP => 'D' ,
##   CYS => 'C' ,
##   GLN => 'Q' ,
##   GLU => 'E' ,
##   GLY => 'G' ,
##   HIS => 'H' ,
##   HSD => 'H' ,
##   HSE => 'H' ,
##   ILE => 'I' ,
##   LEU => 'L' ,
##   LYS => 'K' ,
##   MET => 'M' ,
##   PHE => 'F' ,
##   PRO => 'P' ,
##   SER => 'S' ,
##   THR => 'T' ,
##   TRP => 'W' ,
##   TYR => 'Y' ,
##   VAL => 'V',
##   UNK => 'X',
##   '  C' => 'c',
##   '  G' => 'g',
##   '  A' => 'a',
##   '  T' => 't',
##   '  U' => 'u',
##   '  I' => 'i',
##   'C' => 'c',
##   'G' => 'g',
##   'A' => 'a',
##   'T' => 't',
##   'U' => 'u',
##   'I' => 'i',
##   '+C' => 'c',
##   '+G' => 'g',
##   '+A' => 'a',
##   '+T' => 't',
##   '+U' => 'u',
##   '+I' => 'i'
##   } ;
##
### Per-PPI
### 0. load all of ligands first - lig info
### - iterate over assignments, and have individual
###   bit vectors set for the ligand binding sites in their
###   respective domain families
###
### per-PPI: iterate over each domain--domain interface
###  1. list any ligands actually bound to this interface
###  2. iterate over all ligands in that family/superfamily, find
###     - quantify per ligand coverage
###     - quantify cumulative ligand coverage
###
### per-ligand:
##
##
##   my $fn = set_locations() ;
##   my $liginfo = load_liginfo({fn => $fn}) ;
##
##
##   my $class2lig2bit ;
##   my $lig2class ;
##
##   my $class2bits ;
##   my $class2alnlength;
##   my $class2ligbits ;
##      print STDERR "NOW reading LIG assignment\n" ;
##      open(LIGFH, $fn->{pilig}->{park_assign_lig_fn}) ;
##      while (my $line = <LIGFH>) {
##         if ($line =~ /^#/ ||
##            $line =~ /^Warning: no access to/ ||
##            $line =~ /^Thus no job control/ ) {next;}
##         chomp $line;
##
##         my ( $pdb, $sid, $osid, $classtype, $class, 
##              $btype, $alnposstring, $alnlength,
##              $ligcod, $ligid ) = split(/\t/, $line) ;
##
##         if (exists $standardres->{$ligcod}) {next;}
##         if ($ligcod eq 'DUM' || $ligcod eq 'UNX' ||
##             $ligcod eq 'UNK' || $ligcod eq 'UNL') {next;}
##
###         if (!exists $liginfo->{mwinrange}->{$ligcod}) {next;}
##
##         $class2alnlength->{$classtype}->{$class} = $alnlength ;
##         my $lig = $pdb."\t".$ligcod."\t".$ligid ;
##         $class2ligbits->{all}->{$classtype}->{$class}->{$lig} = Bit::Vector->new($alnlength) ;
##
##         $lig2class->{$ligcod}->{$classtype}->{$class}++ ;
##
##         $ligcod =~ s/ //g ;
##         if (!exists $liginfo->{mw}->{$ligcod}) {
##            print STDERR "WARNING: $fn MW not found for ligand $ligcod\n" ;}
##
##         if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
##            $class2bits->{$classtype}->{$class}->{$btype} =
##               Bit::Vector->new($alnlength) ;
##         }
##
##         my @alnpos = split(/\,/, $alnposstring) ;
##         foreach my $alnpos (@alnpos)  {
##            $class2ligbits->{all}->{$classtype}->{$class}->{$lig}->Bit_On($alnpos) ;}
##      }
##      close(LIGFH) ;
##
##
##      print STDERR "NOW reading EXP assignment\n" ;
##      open(EXPFH, $fn->{pilig}->{park_assign_exp_fn}) ;
##      while (my $line = <EXPFH>) {
##         if ($line =~ /^#/ ||
##            $line =~ /^Warning: no access to/ ||
##            $line =~ /^Thus no job control/ ) {next;}
##         chomp $line;
##
##         my ( $pdb, $sid, $osid, $classtype,
##              $class, $btype, $alnposstring, $alnlength) = split(/\t/, $line) ;
##
##         $class2alnlength->{$classtype}->{$class} = $alnlength ;
##
##         if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
##            $class2bits->{$classtype}->{$class}->{$btype} =
##               Bit::Vector->new($alnlength) ;
##         }
##
##         my @alnpos = split(/\,/, $alnposstring) ;
##         foreach my $alnpos (@alnpos)  {
##            if ($alnpos eq "undef") {next;}
##            $class2bits->{$classtype}->{$class}->{$btype}->Bit_On($alnpos); }
##      }
##      close(EXPFH) ;
##
##
## 
##   my @headers_sum=  ("PINT_LSUM", "PDB", "SID1", "SID2", "CHAINS",
##                       "CLASSTYPE1", "CLASS1", "CLASSTYPE2", "CLASS2",
##                       "numres_p_1", "cumulative_numres_l_and_p_1",
##                       "max_liginrange_1",
##                       "max_l_and_p_1", "lig_max_l_and_p_1",
##                       "max_obi_1", "lig_max_obi_1",
##                       "max_opi_1", "lig_max_opi_1",
##                       "max_olig_1", "lig_max_olig_1",
##                       "numres_p_2", "cumulative_numres_l_and_p_2",
##                       "max_liginrange_2",
##                       "max_l_and_p_2", "lig_max_l_and_p_2",
##                       "max_obi_2", "lig_max_obi_2",
##                       "max_opi_2", "lig_max_opi_2",
##                       "max_olig_2", "lig_max_olig_2",
##                       ) ;
##   print '#'.join("\t",@headers_sum)."\n" ;
##
##   print STDERR "NOW reading PI assignment\n" ;
##   open(PIFH, $fn->{pilig}->{park_assign_pi_fn}) ;
##   my $interfaces ;
##   while (my $line = <PIFH>) {
##      if ($line =~ /^#/ ||
##            $line =~ /^Warning: no access to/ ||
##            $line =~ /^Thus no job control/ ) {next;}
##      chomp $line;
##
##      my ($pdb, $sid, $osid, $classtype, $class, $bstype, $alnposstring,
##             $alnlength, $sid1, $sid2, $fam1, $fam2, $chains)=
##             split(/\t/, $line) ;
##      my $sid12 = $sid1."\t".$sid2 ;
##
##      my $side = 1; if ($sid eq $sid2) {$side = 2;}
##      $interfaces->{$sid12}->{$side}->{osid} = $osid ;
##      $interfaces->{$sid12}->{$side}->{$classtype}->{class} = $class ;
##      $interfaces->{$sid12}->{chains} = $chains ;
##      $interfaces->{$sid12}->{pdb} = $pdb ;
##
###      print STDERR "sid12 = $sid12\nsid =$side\nclasstype=$classtype\nalnposstring=$alnposstring\n" ;
##      $interfaces->{$sid12}->{$side}->{$classtype}->{alnposstring} =
##            $alnposstring;
##
##      $interfaces->{$sid12}->{$side}->{$classtype}->{alnlength} =
##            $alnlength;
##
##      $interfaces->{$sid12}->{$side}->{sid} = $sid2 ;
##   }
##   close(PIFH) ;
##
##   foreach my $sid12 (keys %{$interfaces}) {
##      my $curi = $interfaces->{$sid12} ;
##      my $pdb = $curi->{pdb} ;
##      my $chains = $curi->{chains} ;
##
##      my ($sid, $alnposstring, $alnlength, $class, $classtype) ;
##
###fam level analysis only
##      $alnlength->{1} = $curi->{1}->{fam}->{alnlength} ;
##      $alnposstring->{1} = $curi->{1}->{fam}->{alnposstring} ;
##
##      $alnlength->{2} = $curi->{2}->{fam}->{alnlength} ;
##      $alnposstring->{2} = $curi->{2}->{fam}->{alnposstring} ;
##
##      ($sid->{1}, $sid->{2}) = split(/\t/, $sid12) ;
##
##      $classtype->{1} = 'fam' ;
##      $classtype->{2} = 'fam' ;
##      $class->{1} = $curi->{1}->{fam}->{class} || 'undef';
##      $class->{2} = $curi->{2}->{fam}->{class} || 'undef';
##
##      my $surf_pos ;
##      $surf_pos->{1} = $class2bits->{$classtype->{1}}->{$class->{1}}->{'E'}->Norm() ;
##      $surf_pos->{2} = $class2bits->{$classtype->{2}}->{$class->{2}}->{'E'}->Norm() ;
##
###         ($pdb, $sid->{1}, $osid->{1}, $classtype->{1}, $class->{1},
###             undef, $alnposstring->{1}, $alnlength->{1},
###             undef, undef, undef, undef, $chains) = split(/\t/, $line1) ;
##
##         my $bs_bits;
##         my $curligs ;
##
##         my $alnpos ;
##         my $pistats ;
##         my $ligstats ; 
##         my $skipside ;
##         foreach my $s (1, 2) {
##            if (!defined $alnposstring->{$s}) {
##               $skipside->{$s}++ ;
##               print STDERR "WARNING: skipping $sid->{1} -- $sid->{2}".
##               " undefined binding site alignment positions for $sid->{$s}\n";
##               next;
##            }
##
##            $curligs->{$s} =
##               $class2ligbits->{all}->{$classtype->{$s}}->{$class->{$s}} ;
##
##            $bs_bits->{$s} = Bit::Vector->new($alnlength->{$s}) ;
##            my @t = split(/\,/, $alnposstring->{$s}) ;
##            {
##               @{$alnpos->{$s}} = () ;
##               foreach my $p ( @t) {
##                  if ($p ne 'undef') {
##                     $bs_bits->{$s}->Bit_On($p) ;
##                     push @{$alnpos->{$s}}, $p } }
##            }
##
##            if ($#{$alnpos->{$s}} < 0 ) {
##               $skipside->{$s}++ ;
##               print STDERR "WARNING: skipping $sid->{1} -- $sid->{2}".
##               " undefined binding site alignment positions for $sid->{$s}\n"; }
##         }
##
### init vals
##         foreach my $s (1, 2) {
##            foreach my $type (qw/p cumlig_l_and_p max_liginrange max_l_and_p max_obi max_opi max_olig/) {
##               $ligstats->{$s}->{$type} = 0 ; }
##            foreach my $type (qw/lig_max_l_and_p lig_max_obi lig_max_opi lig_maxolig/) {
##               $ligstats->{$s}->{$type} = 'U'; }
##         }
##
##         foreach my $s (1, 2) {
##            if (exists $skipside->{$s}) {next;}
##            $pistats->{$s}->{p} = $bs_bits->{$s}->Norm() ;
##            $pistats->{$s}->{max_liginrange} = 0 ;
##            $ligstats->{$s}->{cumlig} = Bit::Vector->new($alnlength->{$s}) ;
##            foreach my $t (qw/obi opi olig l_and_p/) {
##               $ligstats->{$s}->{"max_$t"} = 0 ;
##               $ligstats->{$s}->{"lig_max_$t"} = 'undef' ;
##            }
##
##            my $temp_bits_and = Bit::Vector->new($alnlength->{$s}) ;
##            my $temp_bits_or = Bit::Vector->new($alnlength->{$s}) ;
##
##            if (!exists $curligs->{$s}) {next;}
##            foreach my $lig (keys %{$curligs->{$s}}) {
##               $ligstats->{$s}->{cumlig}->Or($ligstats->{$s}->{cumlig},
##                                             $curligs->{$s}->{$lig}) ;
##
##            #intersect the ligands bit vector and the bit vector for the bindin
##            # site positions of this particular interface;
##               $temp_bits_and->And($curligs->{$s}->{$lig}, $bs_bits->{$s}) ;
##               $temp_bits_or->Or($curligs->{$s}->{$lig}, $bs_bits->{$s}) ;
##
##               my $curs ;
##               $curs->{l_and_p} = $temp_bits_and->Norm() ;
##               $curs->{l_or_p} = $temp_bits_or->Norm() ;
##               $curs->{l} = $curligs->{$s}->{$lig}->Norm() ;
##               $curs->{p} = $bs_bits->{$s}->Norm() ;
##
##               if ($curs->{l_and_p} == 0) {next;}
##
##               my $curoverlap ;
##               $curoverlap->{obi} = $curs->{l_and_p} / $curs->{l_or_p} ;
##               $curoverlap->{opi} = $curs->{l_and_p} / $curs->{p} ;
##               $curoverlap->{olig} = $curs->{l_and_p} / $curs->{l} ;
##
##
##               if (!exists $ligstats->{$s}->{"max_l_and_p"} ||
##                   $curs->{l_and_p} > $ligstats->{$s}->{"max_l_and_p"}) {
##                  $ligstats->{$s}->{"max_l_and_p"} = $curs->{l_and_p} ;
##                  $ligstats->{$s}->{"lig_max_l_and_p"} = $lig ;
##               }
##
##               foreach my $t (qw/obi opi olig/) {
##                  if (!exists $ligstats->{$s}->{"max_$t"} ||
##                      $curoverlap->{$t} > $ligstats->{$s}->{"max_$t"}) {
##                     $ligstats->{$s}->{"max_$t"} = $curoverlap->{$t} ;
##                     $ligstats->{$s}->{"lig_max_$t"} = $lig ;
##                  }
##               }
##
##               my (undef, $ligcod, undef) = split(/\t/, $lig) ;
##               my $liginrange = 0 ;
##               if (exists $liginfo->{mwinrange}->{$ligcod}) {
##                  $ligstats->{$s}->{max_liginrange} = 1 ;
##                  $liginrange = 1 ; }
##
##               my @outvals = ("PINT_LINT",
##                              $pdb, $sid->{1}, $sid->{2},
##                              $classtype->{1}, $class->{1},
##                              $classtype->{2}, $class->{2},
##                              $s, $sid->{$s}, $lig, $liginrange,
##                              $curs->{p}, $curs->{l},
##                              $curs->{l_and_p}, $curs->{l_or_p},
##                              sprintf("%.3f", $curoverlap->{obi}),
##                              sprintf("%.3f", $curoverlap->{opi}),
##                              sprintf("%.3f", $curoverlap->{olig})
##                              ) ;
##               print join("\t", @outvals)."\n";
##            }
##
##            $ligstats->{$s}->{"cumlig_l"} = $ligstats->{$s}->{cumlig}->Norm() ;
##            $temp_bits_and->And($ligstats->{$s}->{cumlig}, $bs_bits->{$s}) ;
##            $ligstats->{$s}->{"cumlig_l_and_p"} = $temp_bits_and->Norm() ;
##         }
##
##         my @outvals =  ("PINT_LSUM", $pdb, $sid->{1}, $sid->{2}, $chains,
##                              $classtype->{1}, $class->{1},
##                              $classtype->{2}, $class->{2}) ;
##
##         foreach my $s ( 1, 2) {
##            if (exists $skipside->{$s} ||
##                $ligstats->{$s}->{"cumlig_l"} == 0 ) {
##
##               if (exists $skipside->{$s}) {
##                  push @outvals, 0, 0, 0 ;
##               } else {
##                  push @outvals, $pistats->{$s}->{p}, 0, 0 ;
##               }
##
##               push @outvals, 0, "U", "U", "U";
##               push @outvals, 0, "U", "U", "U";
##               push @outvals, 0, "U", "U", "U";
##               push @outvals, 0, "U", "U", "U";
##               next;
##            }
##
##            push @outvals, $pistats->{$s}->{p} ;
##            push @outvals, $ligstats->{$s}->{cumlig_l_and_p} ;
##            push @outvals, $ligstats->{$s}->{max_liginrange} ;
##            push @outvals, $ligstats->{$s}->{"max_l_and_p"} ;
##            push @outvals, $ligstats->{$s}->{"lig_max_l_and_p"} ;
##            foreach my $t (qw/obi opi olig/)  {
##               push @outvals, sprintf("%.3f", $ligstats->{$s}->{"max_$t"}) ;
##               push @outvals, $ligstats->{$s}->{"lig_max_$t"} ;
##            }
##         }
##      print join("\t", @outvals)."\n"; 
##   }
##
##}
##
##sub combine_piligexpbits {
##
##
##   my $in = shift ;
##   my $pb = $in->{pb} ;
##   my $ligbits = $in->{ligbits} ;
##   my $liginfo = $in->{liginfo} ;
##   my $pibits = $in->{pibits} ;
##   my $expbits = $in->{expbits} ;
##
##   my $class2bits ;
##   my $class2ligs ;
##   my $class2anligbits ;
##
##   foreach my $classtype (qw/fam sf/) {
##      foreach my $sid (keys %{$ligbits->{$classtype}}) {
##         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;
##         if (exists $pibits->{$classtype}->{$sid}) {
##            foreach my $ligsig (keys %{$ligbits->{$classtype}->{$sid}}) {
##               if ($ligsig eq 'cumulative') {next;}
##               my ($pdb, $ligcod, $ligid) = split(/\t/, $ligsig) ;
##               my $anligbits =
##                  $ligbits->{$classtype}->{$sid}->{$ligsig}->Shadow() ;
##
##               $anligbits->AndNot(
##                  $ligbits->{$classtype}->{$sid}->{$ligsig},
##                  $pibits->{$classtype}->{$sid}->{'P'}) ;
##
##               if ($anligbits->Norm() > 0) {
##                  push @{$class2anligbits->{$classtype}->{$class}},
##                     [$ligsig, $anligbits] ;
##
##                  my @btypes = ("L") ;
##                  if (exists $liginfo->{mwinrange}->{$ligcod}) {
##                     push @btypes, "Lmwinrange"; }
##
##                  foreach my $btype (@btypes) {
##                     $class2ligs->{$classtype}->{$class}->{$btype}->{$ligsig}++;
##                     if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
##                        $class2bits->{$classtype}->{$class}->{$btype} = 
##                           $anligbits->Clone() ;
##                     } else {
##                        $class2bits->{$classtype}->{$class}->{$btype}->Union(
##                          $class2bits->{$classtype}->{$class}->{$btype},
##                          $anligbits);
##                     }
##                  }
##               }
##            }
##
##            foreach my $btype (keys %{$pibits->{$classtype}->{$sid}}) {
##               my $anpibits =
##                  $pibits->{$classtype}->{$sid}->{$btype}->Shadow() ;
##
##               $anpibits->AndNot(
##                  $pibits->{$classtype}->{$sid}->{$btype},
##                  $ligbits->{$classtype}->{$sid}->{cumulative}) ;
##
##               if ($anpibits->Norm() > 0) {
##                  if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
##                     $class2bits->{$classtype}->{$class}->{$btype} = 
##                        $anpibits->Clone() ;
##                  } else {
##                     $class2bits->{$classtype}->{$class}->{$btype}->Union(
##                        $class2bits->{$classtype}->{$class}->{$btype}, $anpibits) ;
##                  }
##               }
##            }
##
##         } else {
##            foreach my $ligsig (keys %{$ligbits->{$classtype}->{$sid}}) {
###               print STDERR "   now on ligand $ligsig\n" ;
##               if ($ligsig eq 'cumulative') {next;}
##               my ($pdb, $ligcod, $ligid) = split(/\t/, $ligsig) ;
##               my @btypes = ("L") ;
##               if (exists $liginfo->{mwinrange}->{$ligcod}) {
##                  push @btypes, "Lmwinrange"; }
##
##               foreach my $btype (@btypes) {
##                  $class2ligs->{$classtype}->{$btype}->{$ligsig}++ ;
##                  if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
##                     $class2bits->{$classtype}->{$class}->{$btype} = 
##                        $ligbits->{$classtype}->{$sid}->{$ligsig}->Clone() ;
##                  } else {
##                     $class2bits->{$classtype}->{$class}->{$btype}->Union(
##                        $class2bits->{$classtype}->{$class}->{$btype},
##                        $ligbits->{$classtype}->{$sid}->{$ligsig}) ;
##                  }
##               }
##            }
##         }
##      }
##
##      foreach my $sid (keys %{$pibits->{$classtype}}) {
##         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;
##         if (!exists $ligbits->{$classtype}->{$sid}) {
##            foreach my $btype (keys %{$pibits->{$classtype}->{$sid}}) {
##               if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
##                  $class2bits->{$classtype}->{$class}->{$btype} = 
##                     $pibits->{$classtype}->{$sid}->{$btype}->Clone() ;
##               } else {
##                  $class2bits->{$classtype}->{$class}->{$btype}->Union(
##                     $class2bits->{$classtype}->{$class}->{$btype},
##                     $pibits->{$classtype}->{$sid}->{$btype}) ;
##               }
##            }
##         }
##      }
##   }
##
##   foreach my $classtype (keys %{$expbits}) {
##      foreach my $class (keys %{$expbits->{$classtype}}) {
##         $class2bits->{$classtype}->{$class}->{E} =
##            $expbits->{$classtype}->{$class}->Clone() ;
##      }
##   }
##
##   return {
##      class2bits => $class2bits,
##      class2ligs => $class2ligs,
##      class2anligbits => $class2anligbits,
##   } ;
##   
##}


sub create_mysql_tables {
# fpd090504_1630 
# purpose: create 8 pilig mysql tables from the 6 pilig output files
#          for the web and mysql interfaces.

   my $in = shift ;
   my $specs ;
   $specs->{results_dir} = "/groups/eddy/home/davisf/work/pilig".
                           "/run/finished_results" ;
   $specs->{results_fn} = {
      'liginfo' => 'liginfo.out',
      'assign_pi_clusters' => 'assign_pi_clusters.out',
      'assign_pepnuci_clusters' => 'assign_pepnuci_clusters.out',
      'assign_lig_clusters' => 'assign_lig_clusters.out',
      'collate_pINT' => 'collate_perinstance_pINT.20090505.out',
      'collate_PINT' => 'collate_perinstance_PINT.20090505.out',
      'collate_pINT_summary' => 'collate_perinstance_pINT.20090505.out',
      'collate_PINT_summary' => 'collate_perinstance_PINT.20090505.out',
   } ;

   $specs->{out_fn} = {
      'ligand_info' => 'pilig_ligand_info.pibase',
      'lig' => 'pilig_lig.pibase',
      'pi' => 'pilig_pi.pibase',
      'pep' => 'pilig_pep.pibase',
      'pi_lig_overlap' => 'pilig_pi_lig_overlap.pibase',
      'pi_lig_overlap_summary' => 'pilig_pi_lig_overlap_summary.pibase',
      'pep_lig_overlap' => 'pilig_pep_lig_overlap.pibase',
      'pep_lig_overlap_summary' => 'pilig_pep_lig_overlap_summary.pibase',
   } ;

   my $tables_to_make = {
      pilig_lig => 1,
   } ;

   my $results_fields ;
   $results_fields->{assign_lig_clusters} = [
      'cluster_num',
      'pdb_id',
      'subset_id',
      'scop_sid',
      'class_type',
      'class',
      'entry_type',
      'alnpos_string',
      'aln_length',
      'ligcode',
      'ligcode_id',
   ] ;

   $results_fields->{assign_pi_clusters} = [
      'cluster_num',
      'pdb_id',
      'subset_id',
      'scop_sid',
      'class_type',
      'class',
      'entry_type',
      'alnpos_string',
      'aln_length',
      'subset_id_1',
      'subset_id_2',
      'class_1',
      'class_2',
      'chains'
   ] ;

   $results_fields->{assign_pepnuci_clusters} = [
      'cluster_num',
      'pdb_id',
      'subset_id',
      'scop_sid',
      'class_type',
      'class',
      'entry_type',
      'chain_subset_id',
      'chain_length',
      'alnpos_string',
      'aln_length',
   ] ;


   $results_fields->{collate_PINT_summary} = [
      'PINT_LSUM',
      'subset_id_1',
      'subset_id_2',
      'chains',
      'class_1',
      'class_2',
   ] ;

   foreach my $side (1, 2) {
      push @{$results_fields->{collate_PINT_summary}},
         "numres_p_$side",
         "cum_l_and_p_$side",
         "max_l_and_p_$side",
         "lig_max_l_and_p_$side";

      foreach (my $j = 0; $j < 100; $j+= 5) {
         push @{$results_fields->{collate_PINT_summary}},
            "cum_l_and_p_".$side."_perseqid_".$j ;
         push @{$results_fields->{collate_PINT_summary}},
            "max_l_and_p_".$side."_perseqid_".$j ;
      }
   }

   $results_fields->{collate_PINT} = [
      'subset_id_1',
      'subset_id_2',
      'chains',
      'class_1',
      'class_2',
      'side',
      'subset_id',
      'ligbs_id',
      'ligbs_subset_id',
      'numres_p_side',
      'numres_p_ident',
      'numres_l',
      'numres_l_ident',
      'numres_l_and_p_side',
      'numres_l_and_p_ident',
      'numres_domain_nongap',
      'numres_domain_ident',
   ] ;


   $results_fields->{collate_pINT_summary} = [
      'pINT_LSUM',
      'subset_id',
      'chain_subset_id',
      'chain_length',
      'class',
      'numres_p',
      'cum_l_and_p',
      'max_l_and_p',
      'lig_max_l_and_p',
   ] ;

   foreach (my $j = 0; $j < 100; $j+= 5) {
      push @{$results_fields->{collate_pINT_summary}},
         "cum_l_and_p_perseqid_".$j ;
      push @{$results_fields->{collate_pINT_summary}},
         "max_l_and_p_perseqid_".$j ;
   }

   $results_fields->{collate_pINT} = [
      'subset_id',
      'chain_subset_id',
      'chain_length',
      'class',
      'ligbs_id',
      'ligbs_subset_id',
      'numres_p',
      'numres_p_ident',
      'numres_l',
      'numres_l_ident',
      'numres_l_and_p',
      'numres_l_and_p_ident',
      'numres_domain_nongap',
      'numres_domain_ident',
   ] ;

   $results_fields->{liginfo} = [
      'ligcode',
      'name',
      'formula',
      'mw',
      'numatoms',
      'numatoms_nonh',
   ] ;

   my $results_f2i = {};
   foreach my $file (keys %{$results_fields}) {
      foreach my $j ( 0 .. $#{$results_fields->{$file}}) {
         $results_f2i->{$file}->{$results_fields->{$file}->[$j]} = $j ; } }


   my $outtable_fields ;
   $outtable_fields->{pilig_ligand_info} = [
      'ligcode',
      'name',
      'formula',
      'mw',
      'numatoms',
      'numatoms_nonh',
   ] ;

   $outtable_fields->{pilig_lig} = [
      'ligbs_id',
      'subset_id',
      'cluster_num',
      'class_type',
      'class',
      'alnpos_string',
      'numres_bs',
      'aln_length',
      'ligcode',
   ] ;

   $outtable_fields->{pilig_pi} = [
      'subset_id_1',
      'subset_id_2',
      'subset_id',
      'cluster_num',
      'class_type',
      'class',
      'alnpos_string',
      'numres_bs',
      'aln_length',
   ] ;

   $outtable_fields->{pilig_pep} = [
      'subset_id',
      'chain_subset_id',
      'cluster_num',
      'class_type',
      'class',
      'chain_length',
      'alnpos_string',
      'numres_bs',
      'aln_length',
   ] ;

   $outtable_fields->{pilig_pi_lig_overlap_summary} = [
      'subset_id_1',
      'subset_id_2',
      'numres_p_1',
      'cum_l_and_p_1',
      'max_l_and_p_1',
      'lig_max_l_and_p_1',
      'cum_l_and_p_1_perseqid_20',
      'max_l_and_p_1_perseqid_20',
      'cum_l_and_p_1_perseqid_50',
      'max_l_and_p_1_perseqid_50',
      'cum_l_and_p_1_perseqid_90',
      'max_l_and_p_1_perseqid_90',
      'numres_p_2',
      'cum_l_and_p_2',
      'max_l_and_p_2',
      'lig_max_l_and_p_2',
      'cum_l_and_p_2_perseqid_20',
      'max_l_and_p_2_perseqid_20',
      'cum_l_and_p_2_perseqid_50',
      'max_l_and_p_2_perseqid_50',
      'cum_l_and_p_2_perseqid_90',
      'max_l_and_p_2_perseqid_90',
   ] ;

   $outtable_fields->{pilig_pep_lig_overlap_summary} = [
      'subset_id',
      'chain_subset_id',
      'chain_length',
      'numres_p',
      'cum_l_and_p',
      'max_l_and_p',
      'lig_max_l_and_p',
      'cum_l_and_p_perseqid_20',
      'max_l_and_p_perseqid_20',
      'cum_l_and_p_perseqid_50',
      'max_l_and_p_perseqid_50',
      'cum_l_and_p_perseqid_90',
      'max_l_and_p_perseqid_90',
   ] ;

   $outtable_fields->{pilig_pep_lig_overlap} = [
      'subset_id',
      'chain_subset_id',
      'chain_length',
      'ligbs_id',
      'ligbs_subset_id',
      'ligcode',
      'numres_p',
      'numres_p_ident',
      'numres_l',
      'numres_l_ident',
      'numres_l_and_p',
      'numres_l_and_p_ident',
      'numres_domain_nongap',
      'numres_domain_ident',
   ] ;

   $outtable_fields->{pilig_pi_lig_overlap} = [
      'subset_id_1',
      'subset_id_2',
      'subset_id',
      'ligbs_id',
      'ligbs_subset_id',
      'ligcode',
      'numres_p_side',
      'numres_p_ident',
      'numres_l',
      'numres_l_ident',
      'numres_l_and_p_side',
      'numres_l_and_p_ident',
      'numres_domain_nongap',
      'numres_domain_ident',
   ] ;

# table 1. liginfo.out;
   if (exists $tables_to_make->{pilig_ligand_info} &&
       $tables_to_make->{pilig_ligand_info} == 1) {
   print STDERR "NOW CREATING: pilig_ligand_info\n" ;
   open(INF, $specs->{results_dir}.'/'.$specs->{results_fn}->{liginfo}) ;
   open(OUTF, ">".$specs->{out_fn}->{ligand_info}) ;
   while (my $line = <INF>) {
      if ($line =~ /^#/) {next;}

      chomp $line;
      my @t = split(/\t/, $line) ;
      my @outvals ;
      foreach my $f ( @{$outtable_fields->{pilig_ligand_info}}) {
         if (!exists $results_f2i->{liginfo}->{$f}) {
            die "field $f not recognized" ; }

         if (!defined $t[$results_f2i->{liginfo}->{$f}] ||
             $t[$results_f2i->{liginfo}->{$f}] eq '') {
            push @outvals, '\N' ;
         } else {
            push @outvals, $t[$results_f2i->{liginfo}->{$f}] ;
         }
      }
      print OUTF join("\t", @outvals)."\n" ;
   }
   close(OUTF) ;
   close(INF) ;
   }


# table 2. assign_lig
   if (exists $tables_to_make->{pilig_lig} &&
       $tables_to_make->{pilig_lig} == 1) {
   print STDERR "NOW CREATING: pilig_lig\n" ;
   open(INF, $specs->{results_dir}.'/'.
             $specs->{results_fn}->{assign_lig_clusters}) ;
   open(OUTF, ">".$specs->{out_fn}->{lig}) ;
   while (my $line = <INF>) {
      if ($line =~ /^#/) {next;}

      chomp $line;
      my @t = split(/\t/, $line) ;
      my @outvals ;
      foreach my $f ( @{$outtable_fields->{pilig_lig}}) {
         if (exists $results_f2i->{assign_lig_clusters}->{$f}) {
            push @outvals, $t[$results_f2i->{assign_lig_clusters}->{$f}] ;
         } elsif ($f eq 'ligbs_id') {
            my $ligbs_id = join(':',
               $t[$results_f2i->{assign_lig_clusters}->{pdb_id}],
               $t[$results_f2i->{assign_lig_clusters}->{ligcode}],
               $t[$results_f2i->{assign_lig_clusters}->{ligcode_id}],
            ) ;
            push @outvals, $ligbs_id ;
         } elsif ($f eq 'numres_bs') {
            my $alnpos_string =
             $t[$results_f2i->{assign_lig_clusters}->{alnpos_string}] ;
            my $commas = 0; $commas += $alnpos_string =~ tr/,/,/; 
            my $numres_bs = $commas + 1 ;
            push @outvals, $numres_bs ;
         } else {
            die "output field $f not recognized" ;
         }
      }
      print OUTF join("\t", @outvals)."\n" ;
   }
   close(OUTF) ;
   close(INF) ;
   }

# table 3. assign_pep
   if (exists $tables_to_make->{pilig_pep} &&
       $tables_to_make->{pilig_pep} == 1) {
   print STDERR "NOW CREATING: pilig_pep\n" ;
   open(INF, $specs->{results_dir}.'/'.
             $specs->{results_fn}->{assign_pepnuci_clusters}) ;
   open(OUTF, ">".$specs->{out_fn}->{pep}) ;
   while (my $line = <INF>) {
      if ($line =~ /^#/) {next;}

      chomp $line;
      my @t = split(/\t/, $line) ;
      if ($t[$results_f2i->{assign_pepnuci_clusters}->{entry_type}] ne 'p') {
         next;}

      my @outvals ;
      foreach my $f ( @{$outtable_fields->{pilig_pep}}) {
         if (exists $results_f2i->{assign_pepnuci_clusters}->{$f}) {
            push @outvals, $t[$results_f2i->{assign_pepnuci_clusters}->{$f}] ;
         } elsif ($f eq 'numres_bs') {
            my $alnpos_string =
             $t[$results_f2i->{assign_pepnuci_clusters}->{alnpos_string}] ;
            my $commas = 0; $commas += $alnpos_string =~ tr/,/,/; 
            my $numres_bs = $commas + 1 ;
            push @outvals, $numres_bs ;
         } else {
            die "output field $f not recognized" ;
         }
      }
      print OUTF join("\t", @outvals)."\n" ;
   }
   close(OUTF) ;
   close(INF) ;
   }

# table 4. assign_pi
   if (exists $tables_to_make->{pilig_pi} &&
       $tables_to_make->{pilig_pi} == 1) {
   print STDERR "NOW CREATING: pilig_pi\n" ;
   open(INF, $specs->{results_dir}.'/'.
             $specs->{results_fn}->{assign_pi_clusters}) ;
   open(OUTF, ">".$specs->{out_fn}->{pi}) ;
   while (my $line = <INF>) {
      if ($line =~ /^#/) {next;}

      chomp $line;
      my @t = split(/\t/, $line) ;
      my @outvals ;
      foreach my $f ( @{$outtable_fields->{pilig_pi}}) {
         if (exists $results_f2i->{assign_pi_clusters}->{$f}) {
            push @outvals, $t[$results_f2i->{assign_pi_clusters}->{$f}] ;
         } elsif ($f eq 'numres_bs') {
            my $alnpos_string =
             $t[$results_f2i->{assign_pi_clusters}->{alnpos_string}] ;
            my $commas = 0; $commas += $alnpos_string =~ tr/,/,/; 
            my $numres_bs = $commas + 1 ;
            push @outvals, $numres_bs ;
         } else {
            die "output field $f not recognized" ;
         }
      }
      print OUTF join("\t", @outvals)."\n" ;
   }
   close(OUTF) ;
   close(INF) ;
   }

   if (exists $tables_to_make->{pilig_pi_lig_overlap_summary} &&
       $tables_to_make->{pilig_pi_lig_overlap_summary} == 1) {
   print STDERR "NOW CREATING: pilig_pi_lig_overlap_summary\n" ;
# table 5. pilig_pi_lig_overlap_summary
   open(INF, $specs->{results_dir}.'/'.
             $specs->{results_fn}->{collate_PINT_summary}) ;
   open(OUTF, ">".$specs->{out_fn}->{pi_lig_overlap_summary}) ;
   while (my $line = <INF>) {
      if ($line !~ /^PINT_LSUM/) {next;}

      chomp $line;
      my @t = split(/\t/, $line) ;
      my @outvals ;
      foreach my $f ( @{$outtable_fields->{pilig_pi_lig_overlap_summary}}) {
         if (exists $results_f2i->{collate_PINT_summary}->{$f}) {
            push @outvals, $t[$results_f2i->{collate_PINT_summary}->{$f}] ;
         } else {
            die "output field $f not recognized" ;
         }
      }
      print OUTF join("\t", @outvals)."\n" ;
   }
   close(OUTF) ;
   close(INF) ;
   }


   if (exists $tables_to_make->{pilig_pep_lig_overlap_summary} &&
       $tables_to_make->{pilig_pep_lig_overlap_summary} == 1) {
   print STDERR "NOW CREATING: pilig_pep_lig_overlap_summary\n" ;
# table 6. pilig_pep_lig_overlap_summary
   open(INF, $specs->{results_dir}.'/'.
             $specs->{results_fn}->{collate_pINT_summary}) ;
   open(OUTF, ">".$specs->{out_fn}->{pep_lig_overlap_summary}) ;
   while (my $line = <INF>) {
      if ($line !~ /^pINT_LSUM/) {next;}

      chomp $line;
      my @t = split(/\t/, $line) ;
      my @outvals ;
      foreach my $f ( @{$outtable_fields->{pilig_pep_lig_overlap_summary}}) {
         if (exists $results_f2i->{collate_pINT_summary}->{$f}) {
            push @outvals, $t[$results_f2i->{collate_pINT_summary}->{$f}] ;
         } else {
            die "output field $f not recognized" ;
         }
      }
      print OUTF join("\t", @outvals)."\n" ;
   }
   close(OUTF) ;
   close(INF) ;
   }


# table 7. pilig_pep_lig_overlap
   if (exists $tables_to_make->{pilig_pep_lig_overlap} &&
       $tables_to_make->{pilig_pep_lig_overlap} == 1) {
   print STDERR "NOW CREATING: pilig_pep_lig_overlap\n" ;
   open(INF, $specs->{results_dir}.'/'.$specs->{results_fn}->{collate_pINT}) ;
   open(OUTF, ">".$specs->{out_fn}->{pep_lig_overlap}) ;
   while (my $line = <INF>) {
      if ($line =~ /^\#/ || $line =~ /^pINT_LSUM/) {next;}

      chomp $line;
      my @t = split(/\t/, $line) ;

      my @outvals ;
      foreach my $f ( @{$outtable_fields->{pilig_pep_lig_overlap}}) {
         if (exists $results_f2i->{collate_pINT}->{$f}) {
            push @outvals, $t[$results_f2i->{collate_pINT}->{$f}] ;
         } elsif ($f eq 'ligcode') {
            my (undef, $ligcode,undef) =
               split(':',$t[$results_f2i->{collate_pINT}->{ligbs_id}]) ;
            push @outvals, $ligcode ;
         } else {
            die "output field $f not recognized" ;
         }
      }
      print OUTF join("\t", @outvals)."\n" ;
   }
   close(OUTF) ;
   close(INF) ;
   }


# table 8. pilig_pi_lig_overlap
   if (exists $tables_to_make->{pilig_pi_lig_overlap}) {
   print STDERR "NOW CREATING: pilig_pi_lig_overlap\n" ;
   open(INF, $specs->{results_dir}.'/'.$specs->{results_fn}->{collate_PINT}) ;
   open(OUTF, ">".$specs->{out_fn}->{pi_lig_overlap}) ;
   my $sid12_sid_ligbs_seen ;
   while (my $line = <INF>) {
      if ($line =~ /^\#/ || $line =~ /^PINT_LSUM/) {next;}

      chomp $line;
      my @t = split(/\t/, $line) ;
      my @outvals ;
      foreach my $f ( @{$outtable_fields->{pilig_pi_lig_overlap}}) {
         if (exists $results_f2i->{collate_PINT}->{$f}) {
            push @outvals, $t[$results_f2i->{collate_PINT}->{$f}] ;
         } elsif ($f eq 'lig_pose_num') {
            my $sid12_sid_ligbs_sig = join("\t",
               $t[$results_f2i->{collate_PINT}->{subset_id_1}],
               $t[$results_f2i->{collate_PINT}->{subset_id_2}],
               $t[$results_f2i->{collate_PINT}->{subset_id}],
               $t[$results_f2i->{collate_PINT}->{ligbs_id}],
            ) ;
            $sid12_sid_ligbs_seen->{$sid12_sid_ligbs_sig}++ ;
            push @outvals, $sid12_sid_ligbs_seen->{$sid12_sid_ligbs_sig} ;
         } elsif ($f eq 'ligcode') {
            my (undef, $ligcode,undef) =
               split(':',$t[$results_f2i->{collate_PINT}->{ligbs_id}]) ;
            push @outvals, $ligcode ;
         } else {
            die "output field $f not recognized" ;
         }
      }
      print OUTF join("\t", @outvals)."\n" ;
   }
   close(OUTF) ;
   close(INF) ;
   }

}

sub create_resno_tables {
# fpd090507_1048
# goal: make a string of realresno that correspond to alnpos_string
#   so that we can easily send binding site positions and overlaps
#   onto rasmol scripts
#Note: assumes everything is fam level

   print STDERR "creating resno tables ".__LINE__."\n";

   my $in = shift;
   my $pilig_specs = set_pilig_specs() ;
   my $pibase_specs = pibase::get_specs() ;

   my $bstype ;
   if (exists $in->{bstype}) {
      $bstype = $in->{bstype} ;
   } else {
      $bstype = 'pilig_lig' ;
   }

   my $dbh ;
   if (exists $in->{dbh}) {
      $dbh->{pi} = $in->{dbh} ;
   } else {
      ($dbh->{pi}) = pibase::connect_pibase() ;
   }

   my $sid2class = pibase::mysql_hashload($dbh->{pi},
      "SELECT subset_id, class FROM subsets") ;

   print STDERR "loading astral: ".localtime() ;
   my $astral = _pilig_astral_preload({
      pibase_specs => $pibase_specs,
      pilig_specs => $pilig_specs,
   }) ;
   print STDERR ", DONE: ".localtime()."\n" ;

   print STDERR "loading pibase: ".localtime() ;
   my $pb = _pilig_tod_pibase_preload() ;
   print STDERR ", DONE: ".localtime()."\n" ;

   if ($bstype eq 'pilig_lig') {
      my ($all_subset_id, $all_ligbs_id, $all_alnpos_string, $all_cluster_num)=
         pibase::mysql_fetchcols($dbh->{pi},
            "SELECT subset_id, ligbs_id, alnpos_string, cluster_num ".
            "FROM pilig_lig WHERE scop_level = \"fam\"") ;
      my $class2sid ;

      foreach my $j (0 .. $#{$all_subset_id}) {
         my $t_sid = $all_subset_id->[$j] ;
         push @{$class2sid->{$sid2class->{$t_sid}}->{$t_sid}}, $j ;
      }

      foreach my $class (sort keys %{$class2sid}) {
         my $class_aln = pibase::ASTRAL::load_asteroids_aln({
               aln_fn => $pibase_specs->{asteroids}->{'fam_aln'}.
                  '/'.$class.'.fasta_aln' ,
               seq_fn => $pibase_specs->{asteroids}->{'fam_seq'}.
                  '/'.$class.'.fa' ,
               raf => $astral->{raf},
               gdseqh => $astral->{gdseqh},
               seqclcont100 => $astral->{seqcl2cont}->{100},
               seqcl100 => $astral->{seqcl}->{100},
               allchains => $pb->{pdbchains}
         }) ;
         foreach my $sid (sort keys %{$class2sid->{$class}}) {
            my $osid = $pb->{sid2osid}->{$sid} ;

            foreach my $j (@{$class2sid->{$class}->{$sid}}) {
               my $alnpos_string = $all_alnpos_string->[$j] ;
               my @alnpos = split(',', $alnpos_string) ;
               my @resno ;
               foreach my $k ( 0 .. $#alnpos) {
                  if (exists $class_aln->{pos2resno}->{$osid}->{$alnpos[$k]}) {
                     push @resno,
                        $class_aln->{pos2resno}->{$osid}->{$alnpos[$k]} ;
                  } else {
                     print STDERR "WARNING: no res found for alignment ".
                     " position ".$alnpos[$k]." (fam $class) in $sid\n";
                  }
               }
               my $resno_string = join(',', @resno) ;
               $resno_string =~ s/\n//g ;
               my @outvals = ($all_ligbs_id->[$j], $all_subset_id->[$j],
                              $all_cluster_num->[$j], 'fam', $resno_string) ;
               print join("\t", @outvals)."\n" ;
            }
         }
      }

   } elsif ($bstype eq 'pilig_pep') {

      my ($all_subset_id, $all_chain_subset_id, $all_alnpos_string,
          $all_cluster_num) = pibase::mysql_fetchcols($dbh->{pi},
            "SELECT subset_id, chain_subset_id, alnpos_string, cluster_num ".
            "FROM pilig_pep WHERE scop_level = \"fam\" AND ".
            "alnpos_string NOT LIKE \"\%undef\%\"") ;

      my $class2sid ;
      foreach my $j (0 .. $#{$all_subset_id}) {
         my $t_sid = $all_subset_id->[$j] ;
         push @{$class2sid->{$sid2class->{$t_sid}}->{$t_sid}}, $j ;
      }

      foreach my $class (sort keys %{$class2sid}) {
         my $class_aln = pibase::ASTRAL::load_asteroids_aln({
               aln_fn => $pibase_specs->{asteroids}->{'fam_aln'}.
                  '/'.$class.'.fasta_aln' ,
               seq_fn => $pibase_specs->{asteroids}->{'fam_seq'}.
                  '/'.$class.'.fa' ,
               raf => $astral->{raf},
               gdseqh => $astral->{gdseqh},
               seqclcont100 => $astral->{seqcl2cont}->{100},
               seqcl100 => $astral->{seqcl}->{100},
               allchains => $pb->{pdbchains}
         }) ;

         foreach my $sid (sort keys %{$class2sid->{$class}}) {
            my $osid = $pb->{sid2osid}->{$sid} ;
            foreach my $j (@{$class2sid->{$class}->{$sid}}) {
               my $alnpos_string = $all_alnpos_string->[$j] ;
               my @alnpos = split(',', $alnpos_string) ;
               my @resno ;
               foreach my $k ( 0 .. $#alnpos) {
                  if (exists $class_aln->{pos2resno}->{$osid}->{$alnpos[$k]}) {
                     push @resno,
                        $class_aln->{pos2resno}->{$osid}->{$alnpos[$k]} ;
                  } else {
                     print STDERR "WARNING: no res found for alignment ".
                     " position ".$alnpos[$k]." (fam $class) in $sid\n";
                  }
               }
               my $resno_string = join(',', @resno) ;
               $resno_string =~ s/\n//g ;
               my @outvals = ($all_subset_id->[$j], $all_chain_subset_id->[$j],
                              $all_cluster_num->[$j], 'fam', $resno_string) ;
               print join("\t", @outvals)."\n" ;
            }
         }
      }

   } elsif ($bstype eq 'pilig_pi') {

      my ($all_subset_id_1, $all_subset_id_2, $all_subset_id,
          $all_alnpos_string, $all_cluster_num) =
          pibase::mysql_fetchcols($dbh->{pi},
            "SELECT subset_id_1, subset_id_2, subset_id, ".
            "alnpos_string, cluster_num ".
            "FROM pilig_pi WHERE scop_level = \"fam\" AND ".
            "alnpos_string NOT LIKE \"\%undef\%\"") ;
      my $class2sid ;

      foreach my $j (0 .. $#{$all_subset_id}) {
         my $t_sid = $all_subset_id->[$j] ;
         push @{$class2sid->{$sid2class->{$t_sid}}->{$t_sid}}, $j ;
      }

      foreach my $class (sort keys %{$class2sid}) {
         my $class_aln = pibase::ASTRAL::load_asteroids_aln({
               aln_fn => $pibase_specs->{asteroids}->{'fam_aln'}.
                  '/'.$class.'.fasta_aln' ,
               seq_fn => $pibase_specs->{asteroids}->{'fam_seq'}.
                  '/'.$class.'.fa' ,
               raf => $astral->{raf},
               gdseqh => $astral->{gdseqh},
               seqclcont100 => $astral->{seqcl2cont}->{100},
               seqcl100 => $astral->{seqcl}->{100},
               allchains => $pb->{pdbchains}
         }) ;
         foreach my $sid (sort keys %{$class2sid->{$class}}) {
            my $osid = $pb->{sid2osid}->{$sid} ;

            foreach my $j (@{$class2sid->{$class}->{$sid}}) {
               my $alnpos_string = $all_alnpos_string->[$j] ;
               my @alnpos = split(',', $alnpos_string) ;
               my @resno ;
               foreach my $k ( 0 .. $#alnpos) {
                  if (exists $class_aln->{pos2resno}->{$osid}->{$alnpos[$k]}) {
                     push @resno,
                        $class_aln->{pos2resno}->{$osid}->{$alnpos[$k]} ;
                  } else {
                     print STDERR "WARNING: no res found for alignment ".
                     " position ".$alnpos[$k]." (fam $class) in $sid\n";
                  }
               }
               my $resno_string = join(',', @resno) ;
               $resno_string =~ s/\n//g ;
               my @outvals = ($all_subset_id_1->[$j], $all_subset_id_2->[$j],
                  $all_subset_id->[$j], $all_cluster_num->[$j],
                  'fam', $resno_string) ;
               print join("\t", @outvals)."\n" ;
            }
         }
      }

   } elsif ($bstype eq 'pilig_pi_lig_overlap') {

# just print this query to a file directly command line mysql
# load up line by line - dont need it all in memory...
      my ($sorted_overlaps_fh, $sorted_overlaps_fn) =
         tempfile("pilig_pi_lig_overlap_classsorted.XXXXX", SUFFIX => ".out") ;
         close($sorted_overlaps_fh) ;
      pibase::mysql_commandline_query({
         db_name => "pibase",
         sql => "SELECT a.subset_id_1, a.subset_id_2, a.subset_id, ".
                "a.ligbs_id, a.ligbs_subset_id, b.class ".
                "FROM pilig_pi_lig_overlap as a, ".
                "subsets as b WHERE a.subset_id = b.subset_id",
         post_sort => '-k6,6',
         out_fn => $sorted_overlaps_fn,
      }) ;
#HERENOW090507_1915  - too large to do the class sort here;
# run the query and sort using command line before reading it in line
# by line...

      open(OVERLAPF, $sorted_overlaps_fn) ;
      my $class_aln ;
      my $last_class ;
      while (my $line = <OVERLAPF>) {
         chomp $line;
         my ($sid1, $sid2, $sid, $ligbs_id, $ligbs_sid, $class) =
            split(/\t/, $line) ;

         if (!defined $last_class || $class ne $last_class) {
            print STDERR "Loading $class alignment\n" ;
            $class_aln = pibase::ASTRAL::load_asteroids_aln({
               aln_fn => $pibase_specs->{asteroids}->{'fam_aln'}.
                  '/'.$class.'.fasta_aln' ,
               seq_fn => $pibase_specs->{asteroids}->{'fam_seq'}.
                  '/'.$class.'.fa' ,
               raf => $astral->{raf},
               gdseqh => $astral->{gdseqh},
               seqclcont100 => $astral->{seqcl2cont}->{100},
               seqcl100 => $astral->{seqcl}->{100},
               allchains => $pb->{pdbchains}
            }) ;
         }

         my $p_alnpos_string = pibase::mysql_singleval($dbh->{pi},
            "SELECT alnpos_string FROM pilig_pi WHERE ".
            "subset_id_1 = \"$sid1\" AND subset_id_2 = \"$sid2\" AND ".
            "subset_id = \"$sid\" AND scop_level = \"fam\"") ;

         my $l_alnpos_string = pibase::mysql_singleval($dbh->{pi},
            "SELECT alnpos_string FROM pilig_lig WHERE ".
            "ligbs_id = \"$ligbs_id\" AND subset_id = \"$ligbs_sid\" ".
            "AND scop_level = \"fam\"");

         my @p_alnpos = split(',', $p_alnpos_string) ;
         my $p_alnpos ; map{$p_alnpos->{$_}++;} @p_alnpos ;
         my @l_alnpos = split(',', $l_alnpos_string) ;

         my @lp_alnpos ;
         foreach my $t_pos (@l_alnpos) {
            if (exists $p_alnpos->{$t_pos}) {push @lp_alnpos, $t_pos;} }

# Need 4 resno strings:
# ligbs_subset_id's alnpos_string mapped onto the subset_id
# ligbs_alnpos_string onto the ligbs_subset_id
# overlapping alnpos mapped onto the subset_id
# overlapping alnpos mapped onto the ligbs_subset_id

         my $tomap_sid2alnpos;
         $tomap_sid2alnpos->{$sid}->{resno_l_onto_pbs} = \@l_alnpos ;
         $tomap_sid2alnpos->{$sid}->{resno_lp_onto_pbs} = \@lp_alnpos ;
         $tomap_sid2alnpos->{$ligbs_sid}->{resno_p_onto_lbs} = \@p_alnpos ;
         $tomap_sid2alnpos->{$ligbs_sid}->{resno_lp_onto_lbs} = \@lp_alnpos ;
         my $resno_map ;

# 1. determine overlapping ligbs:
         foreach my $sid (sort keys %{$tomap_sid2alnpos}) {
            my $osid = $pb->{sid2osid}->{$sid} ;
            foreach my $maptype (keys %{$tomap_sid2alnpos->{$sid}}) {
               my $alnpos = $tomap_sid2alnpos->{$sid}->{$maptype} ;
               my @resno ;
               foreach my $k ( 0 .. $#{$alnpos}) {
                  if (exists $class_aln->{pos2resno}->{$osid}->{$alnpos->[$k]}){
                     push @resno,
                        $class_aln->{pos2resno}->{$osid}->{$alnpos->[$k]} ;
                  } else {
#                     print STDERR "WARNING: no res found for alignment ".
#                     " position ".$alnpos->[$k]." (fam $class) in $sid\n";
                  }
               }
               $resno_map->{$maptype} = join(',', @resno) ;
               $resno_map->{$maptype} =~ s/\n//g ;
            }
         }
         my @outvals = ($sid1, $sid2, $sid, $ligbs_id, $ligbs_sid, 
            $resno_map->{resno_l_onto_pbs},
            $resno_map->{resno_lp_onto_pbs},
            $resno_map->{resno_p_onto_lbs},
            $resno_map->{resno_lp_onto_lbs},
         ) ;
         print join("\t", @outvals)."\n" ;
         $last_class = $class ;
      }
      unlink $sorted_overlaps_fn ;
      close(OVERLAPF) ;

   } elsif ($bstype eq 'pilig_pep_lig_overlap') {

# just print this query to a file directly command line mysql
# load up line by line - dont need it all in memory...
      my ($sorted_overlaps_fh, $sorted_overlaps_fn) =
         tempfile("pilig_pep_lig_overlap_classsorted.XXXXX", SUFFIX => ".out") ;
         close($sorted_overlaps_fh) ;
      pibase::mysql_commandline_query({
         db_name => "pibase",
         sql => "SELECT a.subset_id, a.chain_subset_id, ".
                "a.ligbs_id, a.ligbs_subset_id, b.class ".
                "FROM pilig_pep_lig_overlap as a, ".
                "subsets as b WHERE a.subset_id = b.subset_id",
         post_sort => '-k5,5',
         out_fn => $sorted_overlaps_fn,
      }) ;
#HERENOW090507_1915  - too large to do the class sort here;
# run the query and sort using command line before reading it in line
# by line...

      open(OVERLAPF, $sorted_overlaps_fn) ;
      my $class_aln ;
      my $last_class ;
      while (my $line = <OVERLAPF>) {
         chomp $line;
         my ($sid, $chain_sid, $ligbs_id, $ligbs_sid, $class) =
            split(/\t/, $line) ;

         if (!defined $last_class || $class ne $last_class) {
            print STDERR "Loading $class alignment\n" ;
            $class_aln = pibase::ASTRAL::load_asteroids_aln({
               aln_fn => $pibase_specs->{asteroids}->{'fam_aln'}.
                  '/'.$class.'.fasta_aln' ,
               seq_fn => $pibase_specs->{asteroids}->{'fam_seq'}.
                  '/'.$class.'.fa' ,
               raf => $astral->{raf},
               gdseqh => $astral->{gdseqh},
               seqclcont100 => $astral->{seqcl2cont}->{100},
               seqcl100 => $astral->{seqcl}->{100},
               allchains => $pb->{pdbchains}
            }) ;
         }

         my $p_alnpos_string = pibase::mysql_singleval($dbh->{pi},
            "SELECT alnpos_string FROM pilig_pep WHERE ".
            "subset_id = \"$sid\" AND chain_subset_id = \"$chain_sid\" AND ".
            "scop_level = \"fam\"") ;

         my $l_alnpos_string = pibase::mysql_singleval($dbh->{pi},
            "SELECT alnpos_string FROM pilig_lig WHERE ".
            "ligbs_id = \"$ligbs_id\" AND subset_id = \"$ligbs_sid\" ".
            "AND scop_level = \"fam\"");

         my @p_alnpos = split(',', $p_alnpos_string) ;
         my $p_alnpos ; map{$p_alnpos->{$_}++;} @p_alnpos ;
         my @l_alnpos = split(',', $l_alnpos_string) ;

         my @lp_alnpos ;
         foreach my $t_pos (@l_alnpos) {
            if (exists $p_alnpos->{$t_pos}) {push @lp_alnpos, $t_pos;} }

# Need 4 resno strings:
# ligbs_subset_id's alnpos_string mapped onto the subset_id
# ligbs_alnpos_string onto the ligbs_subset_id
# overlapping alnpos mapped onto the subset_id
# overlapping alnpos mapped onto the ligbs_subset_id

         my $tomap_sid2alnpos;
         $tomap_sid2alnpos->{$sid}->{resno_l_onto_pbs} = \@l_alnpos ;
         $tomap_sid2alnpos->{$sid}->{resno_lp_onto_pbs} = \@lp_alnpos ;
         $tomap_sid2alnpos->{$ligbs_sid}->{resno_p_onto_lbs} = \@p_alnpos ;
         $tomap_sid2alnpos->{$ligbs_sid}->{resno_lp_onto_lbs} = \@lp_alnpos ;
         my $resno_map ;

# 1. determine overlapping ligbs:
         foreach my $sid (sort keys %{$tomap_sid2alnpos}) {
            my $osid = $pb->{sid2osid}->{$sid} ;
            foreach my $maptype (keys %{$tomap_sid2alnpos->{$sid}}) {
               my $alnpos = $tomap_sid2alnpos->{$sid}->{$maptype} ;
               my @resno ;
               foreach my $k ( 0 .. $#{$alnpos}) {
                  if (exists $class_aln->{pos2resno}->{$osid}->{$alnpos->[$k]}){
                     push @resno,
                        $class_aln->{pos2resno}->{$osid}->{$alnpos->[$k]} ;
                  } else {
#                     print STDERR "WARNING: no res found for alignment ".
#                     " position ".$alnpos->[$k]." (fam $class) in $sid\n";
                  }
               }
               $resno_map->{$maptype} = join(',', @resno) ;
               $resno_map->{$maptype} =~ s/\n//g ;
            }
         }
         my @outvals = ($sid, $chain_sid, $ligbs_id, $ligbs_sid, 
            $resno_map->{resno_l_onto_pbs},
            $resno_map->{resno_lp_onto_pbs},
            $resno_map->{resno_p_onto_lbs},
            $resno_map->{resno_lp_onto_lbs},
         ) ;
         print join("\t", @outvals)."\n" ;
         $last_class = $class ;
      }
#      unlink $sorted_overlaps_fn ;
      close(OVERLAPF) ;
   }

}


sub create_overlapaln_tables {
# fpd090507_1048
# goal: make a string of realresno that correspond to alnpos_string
#   so that we can easily send binding site positions and overlaps
#   onto rasmol scripts
#Note: assumes everything is fam level

   print STDERR "creating overlap tables ".__LINE__."\n";

   my $in = shift;
   my $pilig_specs = set_pilig_specs() ;
   my $pibase_specs = pibase::get_specs() ;

   my $bstype ;
   if (exists $in->{bstype}) {
      $bstype = $in->{bstype} ;
   } else {
      $bstype = 'pilig_lig' ;
   }

   my $dbh ;
   if (exists $in->{dbh}) {
      $dbh->{pi} = $in->{dbh} ;
   } else {
      ($dbh->{pi}) = pibase::connect_pibase() ;
   }

   my $sid2class = pibase::mysql_hashload($dbh->{pi},
      "SELECT subset_id, class FROM subsets") ;

   print STDERR "loading astral: ".localtime() ;
   my $astral = _pilig_astral_preload({
      pibase_specs => $pibase_specs,
      pilig_specs => $pilig_specs,
   }) ;
   print STDERR ", DONE: ".localtime()."\n" ;

   print STDERR "loading pibase: ".localtime() ;
   my $pb = _pilig_tod_pibase_preload({
      onlythese => {
         pdbchains => 1,
         sid2osid => 1,
      }
   }) ;
   print STDERR ", DONE: ".localtime()."\n" ;

   if ($bstype eq 'pilig_pi_lig_overlap') {
      my ($sorted_overlaps_fh, $sorted_overlaps_fn) =
         tempfile("pilig_pi_lig_overlaalnp_classsorted.XXXXX",
                  SUFFIX => ".out");
         close($sorted_overlaps_fh) ;
      pibase::mysql_commandline_query({
         db_name => "pibase",
         sql => "SELECT a.subset_id_1, a.subset_id_2, a.subset_id, ".
                "a.ligbs_id, a.ligbs_subset_id, b.class ".
                "FROM pilig_pi_lig_overlap as a, ".
                "subsets as b WHERE a.subset_id = b.subset_id",
         post_sort => '-k6,6',
         out_fn => $sorted_overlaps_fn,
      }) ;

      open(OVERLAPF, $sorted_overlaps_fn) ;
      my $class_aln ;
      my $last_class ;
      while (my $line = <OVERLAPF>) {
         chomp $line;
         my ($sid1, $sid2, $sid, $ligbs_id, $ligbs_sid, $class) =
            split(/\t/, $line) ;

         if (!defined $last_class || $class ne $last_class) {
            print STDERR "Loading $class alignment\n" ;
            $class_aln = pibase::ASTRAL::load_asteroids_aln({
               aln_fn => $pibase_specs->{asteroids}->{'fam_aln'}.
                  '/'.$class.'.fasta_aln' ,
               seq_fn => $pibase_specs->{asteroids}->{'fam_seq'}.
                  '/'.$class.'.fa' ,
               raf => $astral->{raf},
               gdseqh => $astral->{gdseqh},
               seqclcont100 => $astral->{seqcl2cont}->{100},
               seqcl100 => $astral->{seqcl}->{100},
               allchains => $pb->{pdbchains}
            }) ;
         }

         my $p_alnpos_string = pibase::mysql_singleval($dbh->{pi},
            "SELECT alnpos_string FROM pilig_pi WHERE ".
            "subset_id_1 = \"$sid1\" AND subset_id_2 = \"$sid2\" AND ".
            "subset_id = \"$sid\" AND scop_level = \"fam\"") ;

         my $l_alnpos_string = pibase::mysql_singleval($dbh->{pi},
            "SELECT alnpos_string FROM pilig_lig WHERE ".
            "ligbs_id = \"$ligbs_id\" AND subset_id = \"$ligbs_sid\" ".
            "AND scop_level = \"fam\"");

# get alignment string from $class_aln and build the 

         my $alnpos ;
         @{$alnpos->{p}} = split(',', $p_alnpos_string) ;
         @{$alnpos->{l}} = split(',', $l_alnpos_string) ;

         my $osid_p = $pb->{sid2osid}->{$sid} ;
         my $osid_l = $pb->{sid2osid}->{$ligbs_sid} ;

         my $seqalnstring ;
         $seqalnstring->{p} = $class_aln->{aln}->{$osid_p} ;
         $seqalnstring->{l} = $class_aln->{aln}->{$osid_l} ;
         foreach my $t (qw/p l/) {
            $seqalnstring->{$t} =~ tr/[A-W]/[a-w]/ ;
            $seqalnstring->{$t} =~ tr/[Y-Z]/[y-z]/ ;
            foreach my $k ( @{$alnpos->{$t}}) {
               substr($seqalnstring->{$t}, $k , 1) = 
                  uc(substr($seqalnstring->{$t}, $k, 1)) ; } }

         my @outvals = ($sid1, $sid2, $sid, $ligbs_id, $ligbs_sid, 
            $seqalnstring->{p}, $seqalnstring->{l}) ;
         print join("\t", @outvals)."\n" ;
         $last_class = $class ;
      }
      unlink $sorted_overlaps_fn ;
      close(OVERLAPF) ;

   } elsif ($bstype eq 'pilig_pep_lig_overlap') {

# just print this query to a file directly command line mysql
# load up line by line - dont need it all in memory...
      my ($sorted_overlaps_fh, $sorted_overlaps_fn) =
         tempfile("pilig_pep_lig_overlapaln_classsorted.XXXXX",
                  SUFFIX => ".out");
         close($sorted_overlaps_fh) ;
      pibase::mysql_commandline_query({
         db_name => "pibase",
         sql => "SELECT a.subset_id, a.chain_subset_id, ".
                "a.ligbs_id, a.ligbs_subset_id, b.class ".
                "FROM pilig_pep_lig_overlap as a, ".
                "subsets as b WHERE a.subset_id = b.subset_id",
         post_sort => '-k5,5',
         out_fn => $sorted_overlaps_fn,
      }) ;
#HERENOW090507_1915  - too large to do the class sort here;
# run the query and sort using command line before reading it in line
# by line...

      open(OVERLAPF, $sorted_overlaps_fn) ;
      my $class_aln ;
      my $last_class ;
      while (my $line = <OVERLAPF>) {
         chomp $line;
         my ($sid, $chain_sid, $ligbs_id, $ligbs_sid, $class) =
            split(/\t/, $line) ;

         if (!defined $last_class || $class ne $last_class) {
            print STDERR "Loading $class alignment\n" ;
            $class_aln = pibase::ASTRAL::load_asteroids_aln({
               aln_fn => $pibase_specs->{asteroids}->{'fam_aln'}.
                  '/'.$class.'.fasta_aln' ,
               seq_fn => $pibase_specs->{asteroids}->{'fam_seq'}.
                  '/'.$class.'.fa' ,
               raf => $astral->{raf},
               gdseqh => $astral->{gdseqh},
               seqclcont100 => $astral->{seqcl2cont}->{100},
               seqcl100 => $astral->{seqcl}->{100},
               allchains => $pb->{pdbchains}
            }) ;
         }

         my $p_alnpos_string = pibase::mysql_singleval($dbh->{pi},
            "SELECT alnpos_string FROM pilig_pep WHERE ".
            "subset_id = \"$sid\" AND chain_subset_id = \"$chain_sid\" AND ".
            "scop_level = \"fam\"") ;

         my $l_alnpos_string = pibase::mysql_singleval($dbh->{pi},
            "SELECT alnpos_string FROM pilig_lig WHERE ".
            "ligbs_id = \"$ligbs_id\" AND subset_id = \"$ligbs_sid\" ".
            "AND scop_level = \"fam\"");

         my $alnpos ;
         @{$alnpos->{p}} = split(',', $p_alnpos_string) ;
         @{$alnpos->{l}} = split(',', $l_alnpos_string) ;

         my $osid_p = $pb->{sid2osid}->{$sid} ;
         my $osid_l = $pb->{sid2osid}->{$ligbs_sid} ;

         my $seqalnstring ;
         $seqalnstring->{p} = $class_aln->{aln}->{$osid_p} ;
         $seqalnstring->{l} = $class_aln->{aln}->{$osid_l} ;
         foreach my $t (qw/p l/) {
            $seqalnstring->{$t} =~ tr/[A-W]/[a-w]/ ;
            $seqalnstring->{$t} =~ tr/[Y-Z]/[y-z]/ ;
            foreach my $k ( @{$alnpos->{$t}}) {
               substr($seqalnstring->{$t}, $k , 1) = 
                  uc(substr($seqalnstring->{$t}, $k, 1)) ; } }

         my @outvals = ($sid, $chain_sid, $ligbs_id, $ligbs_sid, 
            $seqalnstring->{p}, $seqalnstring->{l}) ;
         print join("\t", @outvals)."\n" ;
         $last_class = $class ;
      }
      unlink $sorted_overlaps_fn ;
      close(OVERLAPF) ;
   }

}


sub calc_pilig_db_stats {
#fpd090514_1726
# Purpose: calculate statistics of pilig table entries for web page

   my $in = shift ;

   my $dbh ;
   if (exists $in->{dbh}) {
      $dbh->{pi} = $in->{dbh} ;
   } else {
      ($dbh->{pi}) = pibase::connect_pibase() ;
   }

#tables: select * from pilig_pep; select * from pili

   my $numentries ;
   ($numentries->{"domain-peptide interfaces"}) = pibase::mysql_singleval($dbh->{pi},
      "SELECT count(*) FROM pilig_pep WHERE scop_level = \"fam\"") ;

   ($numentries->{"ligand binding sites"}) = pibase::mysql_singleval($dbh->{pi},
      "SELECT count(*) FROM pilig_lig WHERE scop_level = \"fam\"") ;

   ($numentries->{"ligand binding site clusters"}) =
      pibase::mysql_singleval($dbh->{pi},
      "SELECT count(distinct cluster_num, class) FROM pilig_lig ".
      "WHERE scop_level = \"fam\"") ;

   ($numentries->{"overlapping pairs of ligand and peptide binding sites"}) =
      pibase::mysql_singleval($dbh->{pi},
      "SELECT count(*) FROM pilig_pep_lig_overlap");

   ($numentries->{"overlapping pairs of ligand and domain binding sites"}) =
      pibase::mysql_singleval($dbh->{pi},
      "SELECT count(*) FROM pilig_pi_lig_overlap");

   print '<TABLE frame="hsides">'."\n" ;
   print '<tr><th><DIV class=bodytext>entry type</DIV></th>'.
         '<th><DIV class=bodytext>number of entries</th><tr>'."\n" ;

   foreach my $type (sort keys %{$numentries}) {
      print '<tr><td><DIV class=bodytext>'.$type.'</DIV></td>'."\n" ;
      my $t_num = $numentries->{$type} ;
      $t_num =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g;
      print '<tr><td><DIV class=bodytext>'.$t_num.'</DIV></td>'."\n" ;
   }
   print '</TABLE>'."\n" ;

}

1 ;
