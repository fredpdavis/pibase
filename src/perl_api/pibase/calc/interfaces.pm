=head1 NAME

pibase::calc::interfaces - perl module to compute protein interfaces.

=head1 DESCRIPTION

Perl module that contains routines for computing structural interfaces

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

package pibase::calc::interfaces;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/cluster_scop_interfaces interface_detect_calc bdp_path_2_id _interface_detect_calc__calc_res_pairs/ ;

use pibase qw/connect_pibase mysql_hashload mysql_fetchcols mysql_hasharrload safe_move locate_binaries/;
use File::Temp qw/tempfile tempdir/ ;
use File::Path qw/mkpath/ ;
use POSIX qw/ceil/ ;
use Sys::Hostname qw/hostname/;

use pibase::ASTRAL ;
use pibase::kdcontacts ;
use pibase::interatomic_contacts qw/contacts_select contacts_select_inter special_contact raw_contacts_select/;


sub interface_detct_BUG_inscode_calc {

   use pibase::kdcontacts ;
   use pibase::interatomic_contacts qw/contacts_select contacts_select_inter special_contact raw_contacts_select/ ;

# Set parameters

   my $movethreshold = 200 ;
   my $param ;
   $param->{'dist_cutoff'} = 6.05 ; ## from Janin 1998 JMB captures H20-mediated

   my $dbspecs = pibase::get_specs() ;
   my $pibase = $dbspecs->{db} ;

# Get parameters for special contactsf rom interatomic_cotnacts,.pm

   my $t_special = pibase::interatomic_contacts::special_params() ;
   my $special_thresh= $t_special->{thresh} ;
   my $special_hbres = $t_special->{hbond} ;
   my $special_sbres = $t_special->{salt} ;


# Set table names and SQL ddl specs for tables that need to be generated. (call load_table_defs())

   my $host = hostname() ;
   my $time = pibase::timestamp() ;

   my $dbh_tod ;

# Load subsets_residues_table and interatomic_contacts_tables into memory indexed with bdp_id

   my $meta_tables;

   print STDERR "reading in subsets_residues_tables: " ;
   {
      my @tod_res = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM subsets_residues_tables") ;
      foreach my $j ( 0 .. $#{$tod_res[0]}) {
         $meta_tables->{subsets_residues}->{$tod_res[0]->[$j]} =
            $tod_res[1]->[$j];
      }
   }
   print STDERR "X\n" ;

# Load the entire bdp_chains table into memory and index with bdp_id

   print STDERR "reading in bdp_chains:" ;
   my $allchains ;
   ( $allchains->{bdp_id},
     $allchains->{real_chain_no},
     $allchains->{real_chain_id} ) =
      pibase::rawselect_tod(
      "SELECT bdp_id, real_chain_no, real_chain_id FROM bdp_chains") ;


   my $allchains_point ;
   foreach my $j (0 .. $#{$allchains->{bdp_id}}) {
      push @{$allchains_point->{$allchains->{bdp_id}->[$j]}}, $j ; }
   print STDERR "X\n" ;


# Load the entire subsets table into memory and index with bdp_id

   print STDERR "reading in subsets:" ;
   my $allsubs ;
   my $allsubs_point ;
   {
      my @tod_res = pibase::rawselect_tod(
         "SELECT bdp_id, subset_id, description, subset_source_id, class ".
         "FROM subsets") ;

      foreach my $j ( 0 .. $#{$tod_res[0]}) {
         if (defined $tod_res[0]->[$j]) { # if bdp_id IS NOT NULL
            push @{$allsubs->{bdp_id}}, $tod_res[0]->[$j] ;
            push @{$allsubs->{subset_id}}, $tod_res[1]->[$j] ;
            push @{$allsubs->{description}}, $tod_res[2]->[$j] ;
            push @{$allsubs->{source_id}}, $tod_res[3]->[$j] ;
            push @{$allsubs->{class}}, $tod_res[4]->[$j] ;
         }
      }

      foreach my $j (0 .. $#{$allsubs->{bdp_id}}) {
         push @{$allsubs_point->{$allsubs->{bdp_id}->[$j]}}, $j ; }
   }
   print STDERR "X\n" ;

# Set and open output files.

   my ($out_fh, $out_fn, $out_dir) ;
   my @outqueue ;

   ($out_fh->{intersub_cont}, $out_fn->{intersub_cont}) =
     tempfile( "intersubset_contacts.$host.$time.XXXXXX",
               SUFFIX => ".$pibase");
   print $out_fn->{intersub_cont}."\n" ;

   ($out_fh->{patch_residues_meta}, $out_fn->{patch_residues_meta}) =
     tempfile( "patch_residues_tables.meta.$host.$time.XXXXXX",
               SUFFIX => ".$pibase");
   print $out_fn->{patch_residues_meta}."\n" ;

   ($out_fh->{interface_contacts_meta}, $out_fn->{interface_contacts_meta}) =
     tempfile( "interface_contacts_tables.meta.$host.$time.XXXXXX",
               SUFFIX => ".$pibase");
   print $out_fn->{interface_contacts_meta}."\n" ;

   ($out_fh->{interface_contacts_special_meta},
    $out_fn->{interface_contacts_special_meta}) =
     tempfile( "interface_contacts_special_tables.meta.$host.$time.XXXXXX",
               SUFFIX => ".$pibase");
   print $out_fn->{interface_contacts_special_meta}."\n" ;

   my $tempdir = tempdir(CLEANUP=>1) ;
   chdir $tempdir ;

# Iterate through the bdp_ids.

   my $out_meta ;
   my $lastnowonsig = ' ';
   print STDERR "now on:  " ;

   while (my $line = <STDIN>) {
      chomp $line ;

      my ($bdp_id, $bdp_path) = split(/\t/, $line) ;
      my $tables ;

# Get the file name that holds subsets_residues

      $tables->{subsets_residues}->{sourcefile} =
         $meta_tables->{subsets_residues}->{$bdp_id} ;

# Create interface_contacts file

      $out_fn->{interface_contacts} = '' ;
      $out_fh->{interface_contacts} = '' ;
      ($out_fh->{interface_contacts}, $out_fn->{interface_contacts}) =
        tempfile( "interface_contacts_$bdp_id.XXXXXX", SUFFIX => ".$pibase") ;
      $out_dir->{interface_contacts} = $dbspecs->{metatod_dir}->{interface_contacts} ;
      $out_meta->{interface_contacts} = $out_fh->{interface_contacts_meta} ;


      $out_fn->{interface_contacts_special} = '' ;
      $out_fh->{interface_contacts_special} = '' ;
      ($out_fh->{interface_contacts_special}, $out_fn->{interface_contacts_special}) =
        tempfile("interface_contacts_special_$bdp_id.XXXXXX",SUFFIX => ".$pibase");
      $out_dir->{interface_contacts_special} = $dbspecs->{metatod_dir}->{interface_contacts_special} ;
      $out_meta->{interface_contacts_special} = $out_fh->{interface_contacts_special_meta} ;

# Create table for patch_residues and deposit name in patch_residues_tables.

      $out_fn->{patch_residues} = '' ; $out_fh->{patch_residues} = '' ;
      ($out_fh->{patch_residues}, $out_fn->{patch_residues}) =
         tempfile( "patch_residues_$bdp_id.XXXXXX", SUFFIX => ".$pibase") ;
      $out_dir->{patch_residues} = $dbspecs->{metatod_dir}->{patch_residues} ;
      $out_meta->{patch_residues} = $out_fh->{patch_residues_meta} ;

# OUTPUT (STDERR): Display "bdp_id pdb_id"

      my $nowonsig = $bdp_id ;

#      print STDERR "\b"x(length($lastnowonsig)) ;
      print STDERR $nowonsig ;

# If no subsets loaded for this bdp_id, display to STDERR and abort this bdp_id.

      if (!exists $allsubs_point->{$bdp_id} ) {
         print STDERR "\nSKIPPED bdp $bdp_id: no subsets entry found\n".$nowonsig ;
         $lastnowonsig = $nowonsig ;
         print STDERR "X\n" ;
         next;
      }

# mysql(SELECT subset_id, descsription, subset_source_id, class FROM subsets WHERE bdp_id = <bdp_id>)


      my ( $subset_id, $subset_description, $subset_source_id, $subset_class );
      {
      my @tind = @{$allsubs_point->{$bdp_id}} ;
      push @{$subset_id}, @{$allsubs->{subset_id}}[@tind] ;
      push @{$subset_description}, @{$allsubs->{description}}[@tind] ;
      push @{$subset_source_id}, @{$allsubs->{source_id}}[@tind] ;
      push @{$subset_class}, @{$allsubs->{class}}[@tind] ;
      }

# Iterate through subsets and creating hashes pointing from subset_id to subset_class and subset_id to subset_source_id

      my $rev_subset_class ;
      my $subset_id_2_index ;
      my @subsets_source_id ;
      my @subsets_id ;
      my @subsets_class ;
      foreach my $j ( 0 .. $#{$subset_id}) {
         $rev_subset_class->{$subset_id->[$j]} = $subset_class->[$j] ;
         $subsets_id[$j] = $subset_id->[$j] ;
         $subsets_class[$j] = $subset_class->[$j] ;
         $subsets_source_id[$j] = $subset_source_id->[$j] ;
         $subset_id_2_index->{$subset_id->[$j]} = $j ;
      }

# Read in residue subset assignments from subsets_residues

      my ( $subres_chain_no, $subres_chain_id, $subres_resno,
           $subres_resno_int, $subres_subset_id ) =
         pibase::rawselect_metatod( $tables->{subsets_residues}->{sourcefile},
            "SELECT chain_no, chain_id, ".
            "resno, resno_int, subset_id FROM subsets_residues" ) ;

# If the residue subset assignments is empty, display to STDERR and abort this bdp_id.

      if ($#{$subres_subset_id} < 0 ) {
         print STDERR "\nSKIPPED bdp $bdp_id: no residue subset assignments found\n".$nowonsig ;
         $lastnowonsig = $nowonsig ;
         print STDERR "X\n" ;
         next;
      }

# Load subsets residue assignments into a hash pointing from residue to subsets

      my $subset_assign ;
      foreach my $j (0 .. $#{$subres_chain_no}) {
         if (!defined $subres_chain_id->[$j]) {
            $subres_chain_id->[$j] = ' '; }

# BUG in tables holding subset residue number: blank space issues.

         $subres_resno->[$j] =~ s/ //g ;
         my $sig = $subres_resno->[$j]."\n".$subres_chain_id->[$j] ;
         push @{$subset_assign->{$sig}},
            $subset_id_2_index->{$subres_subset_id->[$j]} ;

      }


      my $subset_assign_flat ;
      foreach my $t_sig (keys %{$subset_assign}) {
         my @t = @{$subset_assign->{$t_sig}} ;
         my @sort_t = sort{$a <=> $b} @t ;
         $subset_assign_flat->{$t_sig} = join(',', @sort_t) ;
      }

# Laod bdp_chains info into a reverse hash  pointing from chain_id to chain_no ;


      my ($t_ch_no, $t_ch_id) ;

      {
      my @tind = @{$allchains_point->{$bdp_id}} ;
      push @{$t_ch_no}, @{$allchains->{real_chain_no}}[@tind] ;
      push @{$t_ch_id}, @{$allchains->{real_chain_id}}[@tind] ;
      }
      pibase::replace_undefs($t_ch_id, ' ') ;
      pibase::replace_char($t_ch_id, '', ' ') ;


      my $bdp_chain_no ;
      foreach my $j ( 0 .. $#{$t_ch_no}) {
         $bdp_chain_no->{$t_ch_id->[$j]} = $t_ch_no->[$j] ; }

# mysql(SELECT chain_id_1, resno_1, resno_1_int, resna_1, chain_id_2, resno_2, resno_2_int, resna_2, count(distance) from interatomic_contacts_<n> WHERE bdp_id = <bdp_id>) group by chain_id_1, resno_1, chain_id_2, resno_2


      my @cont_fields = ('resna1', 'resno1', 'inscode1', 'chain_id1', 'atomno1', 'atomna1', 'resna2', 'resno2', 'inscode2', 'chain_id2', 'atomno2', 'atomna2', 'dist' ) ;
      my $kdfield2no ;
      foreach my $j ( 0 .. $#cont_fields) { $kdfield2no->{$cont_fields[$j]} = $j;}

# Iterate through the contacts, build all res - res contact info vectors

      my $kdcont_out = _interface_detect_BUG_inscode_calc__calc_res_pairs({
         radius => 6.05,
         compress => 1,
         bdp_path => $bdp_path
      }) ;

      if (exists $kdcont_out->{error_fl})  {
         print STDERR "ERROR (bdp_id $bdp_id): $kdcont_out->{error_fl}\n" ;
         next;
      }
      my $kdcontacts_out = $kdcont_out->{contacts_fn} ;


      my $respairs ;
      my $resnames;
      my $specials ;
      my $contacts ;
      open (KDCONT, $kdcontacts_out) ;
      while (my $line = <KDCONT>) {
         chomp $line;
         if ($line =~ /^#/) {next;}
         my @f = split(/\t/, $line) ;
         my $resno_1 = $f[$kdfield2no->{resno1}] ;
         my $resna_1 = $f[$kdfield2no->{resna1}] ;
         my $chain_id_1 = $f[$kdfield2no->{chain_id1}];
         my $atomna_1 = $f[$kdfield2no->{atomna1}];

         my $resno_2 = $f[$kdfield2no->{resno2}] ;
         my $resna_2 = $f[$kdfield2no->{resna2}] ;
         my $chain_id_2 = $f[$kdfield2no->{chain_id2}];
         my $atomna_2 = $f[$kdfield2no->{atomna2}];

         if ((!defined $chain_id_1) || ($chain_id_1 eq '')) {
            $chain_id_1 = ' '; }
         if ((!defined $chain_id_2) || ($chain_id_2 eq '')) {
            $chain_id_2 = ' '; }

# take care of double counting

         my $sig1 = $resno_1."\n".$chain_id_1 ;
         my $sig2 = $resno_2."\n".$chain_id_2 ;
         if ($sig1 lt $sig2) { next; }

         my $dist = $f[$kdfield2no->{dist}] ;
         my $ressig = $sig1."\n".$sig2 ;

         my ($in1, $in2) = (0, 0);
         if (exists $subset_assign_flat->{$sig1}) {
            $in1 = 1; }
         if (exists $subset_assign_flat->{$sig2}) {
            $in2 = 1; }

         my $isdiff =1  ;
         if ((($in1 + $in2) == 2) &&
             ($subset_assign_flat->{$sig1} eq $subset_assign_flat->{$sig2})) {
            $isdiff = 0; }


         if ( (($in1 + $in2) > 0) && $isdiff &&
              (substr($atomna_1, 1, 1) ne 'H') &&
              (substr($atomna_2, 1, 1) ne 'H') &&
              (substr($atomna_1, 1, 1) ne 'Q') &&
              (substr($atomna_2, 1, 1) ne 'Q') &&
              !( $resno_1 eq $resno_2 &&
                 $chain_id_1 eq $chain_id_2) ) {

            if (!exists $resnames->{$sig1}) {
               $resnames->{$sig1} = $resna_1 ; }
            if (!exists $resnames->{$sig2}) {
               $resnames->{$sig2} = $resna_2 ; }


            if (!exists $respairs->{$ressig}) {
               $respairs->{$ressig}->{counts_4} = 0 ;
               $respairs->{$ressig}->{counts_4p5} = 0 ;
               $respairs->{$ressig}->{counts_5} = 0 ;
               $respairs->{$ressig}->{counts_5p5} = 0 ;
            }

            if ((!exists $respairs->{$ressig}->{min_dist}) ||
                ($dist <= $respairs->{$ressig}->{min_dist})) {
               $respairs->{$ressig}->{min_dist} = $dist ;
            }

            $respairs->{$ressig}->{contacts}++ ;
            if ($dist <= 4) {
               $respairs->{$ressig}->{counts_4}++ ;
               $respairs->{$ressig}->{counts_4p5}++ ;
               $respairs->{$ressig}->{counts_5}++ ;
               $respairs->{$ressig}->{counts_5p5}++ ;
            } elsif ($dist <= 4.5) {
               $respairs->{$ressig}->{counts_4p5}++ ;
               $respairs->{$ressig}->{counts_5}++ ;
               $respairs->{$ressig}->{counts_5p5}++ ;
            } elsif ($dist <= 5) {
               $respairs->{$ressig}->{counts_5}++ ;
               $respairs->{$ressig}->{counts_5p5}++ ;
            } elsif ($dist <= 5.5) {
               $respairs->{$ressig}->{counts_5p5}++ ;
            }

# internalized special_contact from itneratmoic_contacts.pm

            {
               my $poss ;
               my $hbflag1 ; my $hbflag2 ;

               if (exists $special_hbres->{$atomna_1}) {
                  $hbflag1 = $special_hbres->{$atomna_1} ;
               } else {
                  $hbflag1 = $special_hbres->{$resna_1}->{$atomna_1} ; }

               if (exists $special_hbres->{$atomna_2}) {
                  $hbflag2 = $special_hbres->{$atomna_2} ;
               } else {
                  $hbflag2 = $special_hbres->{$resna_2}->{$atomna_2} ; }

               my $sbflag1 ; my $sbflag2 ;
               if (exists $special_sbres->{$atomna_1}) {
                  $sbflag1 = $special_sbres->{$atomna_1} ;
               } else {
                  $sbflag1 = $special_sbres->{$resna_1}->{$atomna_1} ; }

               if (exists $special_sbres->{$atomna_2}) {
                  $sbflag2 = $special_sbres->{$atomna_2} ;
               } else {
                  $sbflag2 = $special_sbres->{$resna_2}->{$atomna_2} ; }


               if (($dist <= $special_thresh->{ssbond}) &&
                   ($atomna_1.$atomna_2 eq ' SG  SG ' )) {
                  $specials->{$ressig}->{'ssbond'}++ ;
               } elsif (($dist <= $special_thresh->{salt}) &&
                        (defined $sbflag1 &&
                         defined $sbflag2 &&
                         $sbflag1 + $sbflag2 == 0 ) ) {
                  $specials->{$ressig}->{'salt'}++ ;
               } else {

                  my $thresh = $special_thresh->{hbond} ;
                  if ($atomna_1 eq ' SG ' || $atomna_1 eq ' SD ' ||
                      $atomna_2 eq ' SG ' || $atomna_2 eq ' SD ' ) {
                      $thresh = $special_thresh->{hbond_s} ; }

                  if ($dist < $thresh) {
                     if ((defined $hbflag1 && defined $hbflag2) &&
                        (abs($hbflag1 + $hbflag2) < 2 )) {
                        $specials->{$ressig}->{'hbond'}++ ;
                     }
                  }
               }
            }
         }
      }

      close(KDCONT) ;
      unlink $kdcontacts_out ;

# Iterate through the residue pairs
# For each new interacting atom (residue? fpd030714_0531), iterate through BDP subsets, and find ALL correspding domains, then output for import into subsets_residue.

      my $intersubset_contacts ;

      my $contact_counts ;
      $contact_counts->{4} = 0 ;
      $contact_counts->{4.5} = 0 ;
      $contact_counts->{5} = 0 ;
      $contact_counts->{5.5} = 0 ;
      $contact_counts->{real} = 0 ;

      $contact_counts->{ssbond} = 0 ;
      $contact_counts->{hbond} = 0 ;
      $contact_counts->{saltbridge} = 0 ;

      my $patch_residues ; #->{subset_id_index}->{res\nchain_no} = # contacts

      foreach my $respair (keys %{$respairs}) {

         my ($resno_1, $chain_id_1, $resno_2, $chain_id_2) =
            split(/\n/, $respair) ;

         my $sig1 = $resno_1."\n".$chain_id_1 ;
         my $sig2 = $resno_2."\n".$chain_id_2 ;

# Go looking for subset membership of res 1 and 2.

         my @con_resno = ( $resno_1, $resno_2 ) ;
         my @con_resna = ( $resnames->{$sig1}, $resnames->{$sig2} ) ;
         my @con_chain_id = ( $chain_id_1, $chain_id_2 ) ;
         my @sigs = ( $sig1, $sig2) ;

#fixed BUG: make sure CATH v SCOP merging of domains doesn't preclude
# SCOP interface from being detected
# i.e. The CATH domain def may be the same, but the SCOPs are different

         my $diffsubset = 0;
         my $diffsubsets ;
         foreach my $l ( 0 .. $#{$subset_assign->{$sigs[0]}} ) {
            my $a = $subset_assign->{$sigs[0]}->[$l] ;
            my $a_source = $subsets_source_id[$a] ;
            foreach my $m ( 0 .. $#{$subset_assign->{$sigs[1]}} ) {
               my $b = $subset_assign->{$sigs[1]}->[$m] ;
               my $b_source = $subsets_source_id[$b] ;
               if ($a_source eq $b_source) {
                  if ($a != $b) {
                     $diffsubsets->{$a} = 1 ; $diffsubsets->{$b} = 1 ;
                     $diffsubset++ ;
                  } else {
                     $diffsubsets->{$a} = 0 ; $diffsubsets->{$b} = 0 ;
                  }
               }
            }
         }

         foreach my $l ( 0 .. $#{$subset_assign->{$sigs[0]}}) {
            my $a = $subset_assign->{$sigs[0]}->[$l] ;
            if ($diffsubsets->{$a}) {
               $patch_residues->{$a}->{$con_chain_id[0]}->{$con_resno[0]}->{$con_resna[0]} += $respairs->{$respair}->{contacts} ;
            }
         }

         foreach my $l ( 0 .. $#{$subset_assign->{$sigs[1]}}) {
            my $a = $subset_assign->{$sigs[1]}->[$l] ;
            if ($diffsubsets->{$a}) {
               $patch_residues->{$a}->{$con_chain_id[1]}->{$con_resno[1]}->{$con_resna[1]} += $respairs->{$respair}->{contacts} ;
            }
         }

         if ( ( $#{$subset_assign->{$sigs[0]}} >= 0 ) &&
              ( $#{$subset_assign->{$sigs[1]}} >= 0 ) ) {

            foreach my $l ( 0 .. $#{$subset_assign->{$sigs[0]}} ) {

               my $a = $subset_assign->{$sigs[0]}->[$l] ;
               foreach my $m ( 0 .. $#{$subset_assign->{$sigs[1]}} ) {

                  my $b = $subset_assign->{$sigs[1]}->[$m] ;

                  if ( ( $a != $b ) &&
                        ( $subsets_source_id[$a] ==
                          $subsets_source_id[$b]) ) {

		     my $one_n; my $two_n ;
		     my $one_a; my $two_a ;
		     if ($subsets_id[$a] lt $subsets_id[$b]) {
		        $one_n = 0 ; $two_n = 1 ;
			$one_a = $a ; $two_a = $b ;
                     } else {
		        $one_n = 1 ; $two_n = 0 ;
			$one_a = $b ; $two_a = $a ;
		     }
		     
		     
		     if (!(exists $intersubset_contacts->{$one_a}->{$two_a})) {
			$intersubset_contacts->{$one_a}->{$two_a}->{contacts}=0;
			$intersubset_contacts->{$one_a}->{$two_a}->{counts_4}=0;
			$intersubset_contacts->{$one_a}->{$two_a}->{counts_4p5}=0;
			$intersubset_contacts->{$one_a}->{$two_a}->{counts_5}=0;
			$intersubset_contacts->{$one_a}->{$two_a}->{counts_5p5}=0;
			$intersubset_contacts->{$one_a}->{$two_a}->{hbond}=0;
			$intersubset_contacts->{$one_a}->{$two_a}->{ssbond}=0;
			$intersubset_contacts->{$one_a}->{$two_a}->{salt}=0;
			$intersubset_contacts->{$one_a}->{$two_a}->{min_dist}=
			   $respairs->{$respair}->{min_dist} ;
		     }

		     if ($con_chain_id[$one_n] eq $con_chain_id[$two_n]) {
			$intersubset_contacts->{$one_a}->{$two_a}->{same}++ ;
		     } else {
			$intersubset_contacts->{$one_a}->{$two_a}->{diff}++ ;
		     }

		     $intersubset_contacts->{$one_a}->{$two_a}->{contacts} +=
		        $respairs->{$respair}->{contacts};

		     $intersubset_contacts->{$one_a}->{$two_a}->{counts_4} +=
		        $respairs->{$respair}->{counts_4};
		     $intersubset_contacts->{$one_a}->{$two_a}->{counts_4p5} +=
		        $respairs->{$respair}->{counts_4p5};
		     $intersubset_contacts->{$one_a}->{$two_a}->{counts_5} +=
		        $respairs->{$respair}->{counts_5};
		     $intersubset_contacts->{$one_a}->{$two_a}->{counts_5p5} +=
		        $respairs->{$respair}->{counts_5p5};

		     if (exists $specials->{$respair}) {
		        foreach my $type (keys %{$specials->{$respair}}) {
			   if (exists $specials->{$respair}->{$type}) {
			      $intersubset_contacts->{$one_a}->{$two_a}->{$type}++; }
		           my $special_fields = [
			                 $bdp_id,
					 $subsets_id[$one_a],
					 $subsets_id[$two_a],
					 $bdp_chain_no->{$con_chain_id[$one_n]},
					 $con_chain_id[$one_n],
					 $con_resno[$one_n],
					 $con_resna[$one_n],
					 $bdp_chain_no->{$con_chain_id[$two_n]},
					 $con_chain_id[$two_n],
					 $con_resno[$two_n],
					 $con_resna[$two_n] ,
					 $type
				] ;
		           print {$out_fh->{interface_contacts_special}}
			      join("\t", @{$special_fields})."\n" ;
		        }
		     }

		     my $interface_contacts_fields=[ $bdp_id,
					 $subsets_id[$one_a],
					 $subsets_id[$two_a],
					 $bdp_chain_no->{$con_chain_id[$one_n]},
					 $con_chain_id[$one_n],
					 $con_resno[$one_n],
					 $con_resna[$one_n],
					 $bdp_chain_no->{$con_chain_id[$two_n]},
					 $con_chain_id[$two_n],
					 $con_resno[$two_n],
					 $con_resna[$two_n] ,
	                                 $respairs->{$respair}->{min_dist},
	                                 $respairs->{$respair}->{contacts},
	                                 $respairs->{$respair}->{counts_4},
	                                 $respairs->{$respair}->{counts_4p5},
	                                 $respairs->{$respair}->{counts_5},
	                                 $respairs->{$respair}->{counts_5p5} ] ;

		     foreach my $ohshit (0 .. $#{$interface_contacts_fields}){
		        if (!defined $interface_contacts_fields->[$ohshit]) {
			   print STDERR " field $ohshit undef\n" ; } }
		     print {$out_fh->{interface_contacts}}
		        join("\t", @{$interface_contacts_fields})."\n";

                  }
	       }
	    }
	 }

      }
      close ($out_fh->{interface_contacts}) ;
      close ($out_fh->{interface_contacts_special}) ;

      system("gzip $out_fn->{interface_contacts}") ;
      $out_fn->{interface_contacts} .= '.gz' ;

      system("gzip $out_fn->{interface_contacts_special}") ;
      $out_fn->{interface_contacts_special} .= '.gz' ;

# Iterate through patch_residues hash and display to PATCH_RESIDUES file.

      foreach my $a ( keys %{$patch_residues} ) {
         foreach my $b ( keys %{$patch_residues->{$a}}) {
            foreach my $c ( keys %{$patch_residues->{$a}->{$b}}) {
               foreach my $d ( keys %{$patch_residues->{$a}->{$b}->{$c}}) {
		        my @fields = ( $bdp_id,
			               $subsets_id[$a],
			               $bdp_chain_no->{$b}, $b, $c, $d,
				       $patch_residues->{$a}->{$b}->{$c}->{$d}
			             ) ;
	                print {$out_fh->{patch_residues}} join("\t", @fields)."\n" ; } } } }
      close ($out_fh->{patch_residues}) ;
      system("gzip $out_fn->{patch_residues}") ;
      $out_fn->{patch_residues} .= '.gz' ;


# OUTPUT: Display output file names of non-empty output files.

      my @output = qw/patch_residues interface_contacts interface_contacts_special/ ;
      foreach my $outf (@output) {
         if (-z $out_fn->{$outf}) {
	    unlink $out_fn->{$outf} ;
	 } else {
            push @outqueue, {file => $out_fn->{$outf}, dir=>$out_dir->{$outf}};
            my @outvals = ($bdp_id, $param->{dist_cutoff},
	                   '\N',$out_dir->{$outf}."/".$out_fn->{$outf});
            print {$out_meta->{$outf}} join("\t", @outvals)."\n";
	 }
      }


# Iterate through the first subset partner of all contact

      foreach my $subset1 (keys %{$intersubset_contacts}) {

# Iterate through all interaction partners of subset1

         foreach my $subset2 (keys %{$intersubset_contacts->{$subset1}}) {

# Display contacting subsets for intersubset_contacts import.
# print(INTERSUBSET_CONTACTS	bdp_id	subset1	subset2)

            my $t_chainstat ;
            if (exists $intersubset_contacts->{$subset1}->{$subset2}->{same}) {
               $t_chainstat = 'same';
               if (exists $intersubset_contacts->{$subset1}->{$subset2}->{diff}) {
                  $t_chainstat = 'both' ;
               }
            } else {
               $t_chainstat = 'diff' ;
            }

            my @outfields = (
               $bdp_id,
               $subsets_id[$subset1],
               $subsets_class[$subset1] ,
               $subsets_id[$subset2] ,
               $subsets_class[$subset2],
               $intersubset_contacts->{$subset1}->{$subset2}->{contacts},
               $param->{dist_cutoff},
               $intersubset_contacts->{$subset1}->{$subset2}->{counts_4},
               $intersubset_contacts->{$subset1}->{$subset2}->{counts_4p5},
               $intersubset_contacts->{$subset1}->{$subset2}->{counts_5},
               $intersubset_contacts->{$subset1}->{$subset2}->{counts_5p5},
               $intersubset_contacts->{$subset1}->{$subset2}->{hbond},
               $intersubset_contacts->{$subset1}->{$subset2}->{salt},
               $intersubset_contacts->{$subset1}->{$subset2}->{ssbond},
               $t_chainstat
            ) ;

            print {$out_fh->{intersub_cont}} join("\t", @outfields)."\n" ;

         }

      }

      $lastnowonsig = $nowonsig ;

      if ($#outqueue >= $movethreshold) {
         foreach my $j ( 0 .. $#outqueue) {
	    pibase::safe_move($outqueue[$j]->{file}, $outqueue[$j]->{dir}) ; }
         @outqueue = () ;
      }

      print STDERR "X\n" ;
   }
   
   foreach my $j ( 0 .. $#outqueue) {
      pibase::safe_move($outqueue[$j]->{file}, $outqueue[$j]->{dir}) ; }

   close($out_fh->{intersub_cont}) ;
   close($out_fh->{patch_residues_meta}) ;
   close($out_fh->{interface_contacts_meta}) ;

   print STDERR "\n" ;

}

sub _interface_detect_BUG_inscode_calc__calc_res_pairs {
   my $params = shift ;

   my $kdcontacts_radius = $params->{radius} || "6.6" ;
   my $compress_fl = $params->{compress} || 1 ;

   my $binaries = pibase::locate_binaries() ;
   my $kdcontacts_bin = $binaries->{'kdcontacts'}." $kdcontacts_radius" ;
   my $altloc_check = $binaries->{'altloc_check'} ;
   my $altloc_filter = $binaries->{'altloc_filter'} ;

   my $bdp_file_path = $params->{bdp_path} ;
   my $localbdp = $bdp_file_path ; $localbdp =~ s/.*\/// ;

   pibase::safe_copy($bdp_file_path, $localbdp) ;
   
   if (!-e $localbdp) {
      return {error_fl => "ERROR: pdb file access error - couldnt copy locally"} ;
   }

   my $host = hostname() ;
   my ($fh, $kdcontacts_out) =
      tempfile("kdcontacts.$host.XXXXXX", SUFFIX => ".out"); close($fh) ;
   my ($fh2, $kdcontacts_err) =
      tempfile("kdcontacts.$host.XXXXXX", SUFFIX => ".err"); close($fh2) ;

# Check if the pdb file contains altloc identifiers. If so, first filter with altloc_filter.pl

   my $altloc_fl = `$altloc_check < $localbdp` ; chomp $altloc_fl ;
   my $tcom ;
   if ($altloc_fl) {
      $tcom = "$altloc_filter $localbdp" ;
   } else {
      $tcom = "cat $localbdp" ;
   }

# system(cat <bdp_file_path> | kdcontacts_bin 2><kdcontacts_err> ><kdcontacts_out>)

   $tcom .= " | $kdcontacts_bin 2>$kdcontacts_err >$kdcontacts_out" ;
   system($tcom) ;

# If the kdcontacts_out file is empty or the _err file exists, display to STDERR and remove the files.

   if (!(-e $kdcontacts_out) || (! -z $kdcontacts_err))  {
      unlink $kdcontacts_err, $kdcontacts_out;
      return {error_fl => "ERROR: $bdp_file_path kdcontacts execution error"} ;

# Otherwise remove the _err file

   } else {

      unlink $kdcontacts_err ;

   }

# Read in kdcontacts output
#  read in kdcontacts output into an array

   my @cont_fields = ('resna1', 'resno1', 'inscode1', 'chain_id1', 'atomno1', 'atomna1', 'resna2', 'resno2', 'inscode2', 'chain_id2', 'atomno2', 'atomna2', 'dist' ) ;
   my $kdfield2no ;
   foreach my $j ( 0 .. $#cont_fields) { $kdfield2no->{$cont_fields[$j]} = $j;}

   my $results = {
      contacts_fn => $kdcontacts_out,
      fields => \@cont_fields,
      field2no => $kdfield2no
   };
   unlink $localbdp ;
   return $results ;

}

=head2 interface_detect_calc()

   Title:       interface_detect_calc()
   Function:    Detect interfaces in bdp files.
   Args:        None
   Returns:     Nothing
   STDIN:       bdp_id."\t".bdp_path
   Files in:    PDB files (as specified in STDIN column 2 (bdp_path))
   Files out:
   o intersubset_contacts.<hostname>.<timestamp>.<XXXXXX>.<pibase db name>
   o patch_residues_tables.meta.<hostname>.<timestamp>.<XXXXXX>.<pibase db name>
   o interface_contacts_tables.meta.<hostname>.<timestamp>.<XXXXXX>.
      <pibase db name>
   o interface_contacts_special_tables.meta.<hostname>.<timestamp>.<XXXXXX>.
      <pibase db name>

   o foreach bdp_id:
      o patch_residues_<bdp_id>.<XXXXXX>.<pibase db name>
      o interface_contacts_<bdp_id>.<XXXXXX>.<pibase db name>
      o interface_contacts_special_<bdp_id>.<XXXXXX>.<pibase db name>

=cut

sub interface_detect_calc {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }
   my $pibase = $pibase_specs->{db} ;
   my $movethreshold = 200 ;

# Set parameters
   if (!exists $in->{dist_cutoff}) {
      $in->{dist_cutoff}= $pibase_specs->{interface_detect_calc}->{dist_cutoff};
   } else {
      $in->{dist_cutoff} = $in->{dist_cutoff} ;
   }

   my ($bdp_id_2_path) = pibase::todload_bdp_ids('bdp_id_2_path') ;

   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         chomp $line;
         my ($bdp_id, $bdp_path) = split(/\t/, $line) ;
         $in->{bdps}->{$bdp_id} = $bdp_path ;
      }
      close(INF) ;
   } else {
      foreach my $bdp_id (keys %{$bdp_id_2_path}) {
         $in->{bdps}->{$bdp_id} = $bdp_id_2_path->{$bdp_id} ;
      }
   }

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) { #master node
# split bdp_ids and send off to compute - everything else is computing.

      print "* interface_detect_calc() ".localtime() if (!exists
         $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{interface_detect_calc_in},
       $temp_fn->{interface_detect_calc_in}) =
         tempfile("splits_interface_detect_calc_input.XXXXX");

      ($temp_fh->{interface_detect_calc_out},
       $temp_fn->{interface_detect_calc_out}) =
         tempfile("splits_interface_detect_calc_SGEout_XXXXX",
         SUFFIX => '.pibase');
         close($temp_fh->{interface_detect_calc_out}) ;

      ($temp_fh->{interface_detect_calc_err},
       $temp_fn->{interface_detect_calc_err}) =
         tempfile("splits_interface_detect_calc_SGEerr_XXXXX",
         SUFFIX => '.pibase');
         close($temp_fh->{interface_detect_calc_err}) ;

      foreach my $bdp_id (sort {$a <=> $b} keys %{$in->{bdps}}) {
         print {$temp_fh->{interface_detect_calc_in}}
            join("\t", $bdp_id, $in->{bdps}->{$bdp_id})."\n" ; }
      close($temp_fh->{interface_detect_calc_in}) ;

      my $cur_numjobs ;
      if (exists $in->{numtasks_cluster}) {
         $cur_numjobs = $in->{numtasks_cluster} ;
      } else {
         $cur_numjobs = $pibase_specs->{SGE}->{numjobs} ;
      }

      my $split_dir = tempdir("splits_interface_detect_calc.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{interface_detect_calc_in},
         dir => $split_dir,
         numjobs => $cur_numjobs
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.interface_detect_calc.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::calc::interfaces qw/interface_detect_calc/ ;

main() ;

sub main {

         pibase::calc::interfaces::interface_detect_calc({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.interface_detect_calc.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.interface_detect_calc.XXXXX");

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
         out_fn => $temp_fn->{interface_detect_calc_out},
         err_fn => $temp_fn->{interface_detect_calc_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open(REALOUTF_INTERSUBCONT,">".
         $pibase_specs->{buildfiles}->{intersubset_contacts}) ;
      open(REALOUTF_PATCHRES_T,">".
         $pibase_specs->{buildfiles}->{patch_residues_tables}) ;
      open(REALOUTF_INTERFACECONT_T,">".
         $pibase_specs->{buildfiles}->{interface_contacts_tables}) ;
      open(REALOUTF_INTERFACECONT_SPECIAL_T,">".
         $pibase_specs->{buildfiles}->{interface_contacts_special_tables}) ;

      open($temp_fh->{interface_detect_calc_out},
           $temp_fn->{interface_detect_calc_out}) ;
      while (my $line = readline($temp_fh->{interface_detect_calc_out})) {
# need to grep out ^SUBSETS, ^SUBSETS_DETAILS, and ^SUBSETS_RESIDUES entries
         if ($line =~ /^\#/) {
            next;
         } elsif ($line =~ /^intersubset_contacts\t/) {
            $line =~ s/^intersubset_contacts\t// ;
            print REALOUTF_INTERSUBCONT $line ;
         } elsif ($line =~ /^patch_residues_tables\t/) {
            $line =~ s/^patch_residues_tables\t// ;
            print REALOUTF_PATCHRES_T $line ;
         } elsif ($line =~ /^interface_contacts_tables\t/) {
            $line =~ s/^interface_contacts_tables\t// ;
            print REALOUTF_INTERFACECONT_T $line ;
         } elsif ($line =~ /^interface_contacts_special_tables\t/) {
            $line =~ s/^interface_contacts_special_tables\t// ;
            print REALOUTF_INTERFACECONT_SPECIAL_T $line ;
         }
      }
      close($temp_fh->{interface_detect_calc_out}) ;
      close(REALOUTF_INTERSUBCONT) ;
      close(REALOUTF_PATCHRES_T) ;
      close(REALOUTF_INTERFACECONT_T) ;
      close(REALOUTF_INTERFACECONT_SPECIAL_T) ;


      my $import_status ;
      if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_chains file (specified in specs) into pibase
         foreach my $import_fn (qw/intersubset_contacts interface_contacts_tables interface_contacts_special_tables patch_residues_tables/) {
            $import_status->{$import_fn} = pibase::mysqlimport({
               pibase_specs => $pibase_specs,
               fn => $pibase_specs->{buildfiles}->{$import_fn}}) ;
         }
      }
      return $import_status ;

   } else { #compute node
   
# Get parameters for special contactsf rom interatomic_contacts.pm
      my $t_special = pibase::interatomic_contacts::special_params() ;
      my $special_thresh= $t_special->{thresh} ;
      my $special_hbres = $t_special->{hbond} ;
      my $special_sbres = $t_special->{salt} ;
   
# Set table names and SQL ddl specs for tables that need to be generated. (call load_table_defs())
      my $host = hostname() ;
      my $time = pibase::timestamp() ;
   
      my $dbh_tod ;
   
# Load subsets_residues_table and interatomic_contacts_tables into memory indexed with bdp_id
      my $meta_tables;
   
      print STDERR "reading in subsets_residues_tables: " ;
      {
         my @tod_res = pibase::rawselect_tod(
            "SELECT bdp_id, source_file FROM subsets_residues_tables") ;
         foreach my $j ( 0 .. $#{$tod_res[0]}) {
            $meta_tables->{subsets_residues}->{$tod_res[0]->[$j]} =
               $tod_res[1]->[$j];
         }
      }
      print STDERR "X\n" ;
   
# Load the entire bdp_chains table into memory and index with bdp_id
   
      print STDERR "reading in bdp_chains:" ;
      my $allchains ;
      ( $allchains->{bdp_id},
        $allchains->{real_chain_no},
        $allchains->{real_chain_id} ) =
         pibase::rawselect_tod(
         "SELECT bdp_id, real_chain_no, real_chain_id FROM bdp_chains") ;
   
      my $allchains_point ;
      foreach my $j (0 .. $#{$allchains->{bdp_id}}) {
         push @{$allchains_point->{$allchains->{bdp_id}->[$j]}}, $j ; }
      print STDERR "X\n" ;
   
   
# Load the entire subsets table into memory and index with bdp_id
      print STDERR "reading in subsets:" ;
      my $allsubs ;
      my $allsubs_point ;
      {
         my @tod_res = pibase::rawselect_tod(
            "SELECT bdp_id, subset_id, description, subset_source_id, class ".
            "FROM subsets") ;
   
         foreach my $j ( 0 .. $#{$tod_res[0]}) {
            if (defined $tod_res[0]->[$j]) { # if bdp_id IS NOT NULL
               push @{$allsubs->{bdp_id}}, $tod_res[0]->[$j] ;
               push @{$allsubs->{subset_id}}, $tod_res[1]->[$j] ;
               push @{$allsubs->{description}}, $tod_res[2]->[$j] ;
               push @{$allsubs->{source_id}}, $tod_res[3]->[$j] ;
               push @{$allsubs->{class}}, $tod_res[4]->[$j] ;
            }
         }
   
         foreach my $j (0 .. $#{$allsubs->{bdp_id}}) {
            push @{$allsubs_point->{$allsubs->{bdp_id}->[$j]}}, $j ; }
      }
      print STDERR "X\n" ;
   
# Set and open output files.
      my ($out_fh, $out_fn, $out_dir) ;
      my @outqueue ;
   
#Switched to direct STDOUT output
#   ($out_fh->{intersub_cont}, $out_fn->{intersub_cont}) =
#     tempfile( "intersubset_contacts.$host.$time.XXXXXX",
#               SUFFIX => ".$pibase");
#   print $out_fn->{intersub_cont}."\n" ;
#
#   ($out_fh->{patch_residues_meta}, $out_fn->{patch_residues_meta}) =
#     tempfile( "patch_residues_tables.meta.$host.$time.XXXXXX",
#               SUFFIX => ".$pibase");
#   print $out_fn->{patch_residues_meta}."\n" ;
#
#   ($out_fh->{interface_contacts_meta}, $out_fn->{interface_contacts_meta}) =
#     tempfile( "interface_contacts_tables.meta.$host.$time.XXXXXX",
#               SUFFIX => ".$pibase");
#   print $out_fn->{interface_contacts_meta}."\n" ;
#
#   ($out_fh->{interface_contacts_special_meta},
#    $out_fn->{interface_contacts_special_meta}) =
#     tempfile( "interface_contacts_special_tables.meta.$host.$time.XXXXXX",
#               SUFFIX => ".$pibase");
#   print $out_fn->{interface_contacts_special_meta}."\n" ;
   
      my $tempdir = tempdir(CLEANUP=>1) ;
      chdir $tempdir ;
   
      my $out_meta ;
   
# Iterate through the bdp_ids.
      foreach my $bdp_id (keys %{$in->{bdps}}) {
         my $bdp_path = $in->{bdps}->{$bdp_id} ;
         print STDERR "now on bdp $bdp_id ($bdp_path)\n" ;
         my $tables ;
   
# Get the file name that holds subsets_residues
         $tables->{subsets_residues}->{sourcefile} =
            $meta_tables->{subsets_residues}->{$bdp_id} ;
   
# Create interface_contacts file
         $out_fn->{interface_contacts} = '' ;
         $out_fh->{interface_contacts} = '' ;
         ($out_fh->{interface_contacts}, $out_fn->{interface_contacts}) =
           tempfile( "interface_contacts_$bdp_id.XXXXXX", SUFFIX => ".$pibase") ;
         $out_dir->{interface_contacts} =
            $pibase_specs->{metatod_dir}->{interface_contacts}.'/'.
                  POSIX::floor($bdp_id / 1000) ;
#      $out_meta->{interface_contacts} =
#         $out_fh->{interface_contacts_meta} ;
   
   
         $out_fn->{interface_contacts_special} = '' ;
         $out_fh->{interface_contacts_special} = '' ;
         ($out_fh->{interface_contacts_special},
          $out_fn->{interface_contacts_special}) =
           tempfile("interface_contacts_special_$bdp_id.XXXXXX",
                    SUFFIX => ".$pibase");
   
         $out_dir->{interface_contacts_special} =
            $pibase_specs->{metatod_dir}->{interface_contacts_special}.'/'.
                  POSIX::floor($bdp_id / 1000) ;
#      $out_meta->{interface_contacts_special} =
#         $out_fh->{interface_contacts_special_meta} ;
   
# Create table for patch_residues and deposit name in patch_residues_tables.
         $out_fn->{patch_residues} = '' ; $out_fh->{patch_residues} = '' ;
         ($out_fh->{patch_residues}, $out_fn->{patch_residues}) =
            tempfile( "patch_residues_$bdp_id.XXXXXX", SUFFIX => ".$pibase") ;
         $out_dir->{patch_residues} =
            $pibase_specs->{metatod_dir}->{patch_residues}.'/'.
                  POSIX::floor($bdp_id / 1000) ;
#      $out_meta->{patch_residues} = $out_fh->{patch_residues_meta} ;
   
# OUTPUT (STDERR): Display "bdp_id pdb_id"
# If no subsets loaded for this bdp_id, display to STDERR and abort this bdp_id.
         if (!exists $allsubs_point->{$bdp_id} ) {
            print STDERR "SKIPPED bdp $bdp_id: no subsets entry found\n";
            next;
         }
   
# mysql(SELECT subset_id, descsription, subset_source_id, class FROM subsets WHERE bdp_id = <bdp_id>)
         my ( $subset_id, $subset_description, $subset_source_id, $subset_class );
         {
         my @tind = @{$allsubs_point->{$bdp_id}} ;
         push @{$subset_id}, @{$allsubs->{subset_id}}[@tind] ;
         push @{$subset_description}, @{$allsubs->{description}}[@tind] ;
         push @{$subset_source_id}, @{$allsubs->{source_id}}[@tind] ;
         push @{$subset_class}, @{$allsubs->{class}}[@tind] ;
         }
   
# Iterate through subsets and creating hashes pointing from subset_id to subset_class and subset_id to subset_source_id
         my $rev_subset_class ;
         my $subset_id_2_index ;
         my @subsets_source_id ;
         my @subsets_id ;
         my @subsets_class ;
         foreach my $j ( 0 .. $#{$subset_id}) {
            $rev_subset_class->{$subset_id->[$j]} = $subset_class->[$j] ;
            $subsets_id[$j] = $subset_id->[$j] ;
            $subsets_class[$j] = $subset_class->[$j] ;
            $subsets_source_id[$j] = $subset_source_id->[$j] ;
            $subset_id_2_index->{$subset_id->[$j]} = $j ;
         }
   
# Read in residue subset assignments from subsets_residues
         my ( $subres_chain_no, $subres_chain_id, $subres_resno,
              $subres_resno_int, $subres_subset_id ) =
   	 pibase::rawselect_metatod( $tables->{subsets_residues}->{sourcefile},
   	 "SELECT chain_no, chain_id, ".
   	 "resno, resno_int, subset_id FROM subsets_residues" ) ;
   
# If the residue subset assignments is empty, display to STDERR and abort this bdp_id.
         if ($#{$subres_subset_id} < 0 ) {
            print STDERR "SKIPPED bdp $bdp_id: no residue ".
                         "subset assignments found\n";
            next;
         }
   
# Load subsets residue assignments into a hash pointing from residue to subsets
         my $subset_assign ;
         foreach my $j (0 .. $#{$subres_chain_no}) {
            if (!defined $subres_chain_id->[$j]) {
               $subres_chain_id->[$j] = ' '; }
   
# BUG in tables holding subset residue number: blank space issues.
               $subres_resno->[$j] =~ s/ //g ;
               my $sig = $subres_resno->[$j]."\n".$subres_chain_id->[$j] ;
               push @{$subset_assign->{$sig}},
                  $subset_id_2_index->{$subres_subset_id->[$j]} ;
         }
   
         my $subset_assign_flat ;
         foreach my $t_sig (keys %{$subset_assign}) {
            my @t = @{$subset_assign->{$t_sig}} ;
            my @sort_t = sort{$a <=> $b} @t ;
            $subset_assign_flat->{$t_sig} = join(',', @sort_t) ;
         }
   
# Load bdp_chains info into a reverse hash  pointing from chain_id to chain_no ;
         my ($t_ch_no, $t_ch_id) ;
   
         {
         my @tind = @{$allchains_point->{$bdp_id}} ;
         push @{$t_ch_no}, @{$allchains->{real_chain_no}}[@tind] ;
         push @{$t_ch_id}, @{$allchains->{real_chain_id}}[@tind] ;
         }
         pibase::replace_undefs($t_ch_id, ' ') ;
         pibase::replace_char($t_ch_id, '', ' ') ;
   
   
         my $bdp_chain_no ;
         foreach my $j ( 0 .. $#{$t_ch_no}) {
            $bdp_chain_no->{$t_ch_id->[$j]} = $t_ch_no->[$j] ; }
   
# mysql(SELECT chain_id_1, resno_1, resno_1_int, resna_1, chain_id_2, resno_2, resno_2_int, resna_2, count(distance) from interatomic_contacts_<n> WHERE bdp_id = <bdp_id>) group by chain_id_1, resno_1, chain_id_2, resno_2
         my @cont_fields = ('resna1', 'resno1', 'inscode1', 'chain_id1', 'atomno1', 'atomna1', 'resna2', 'resno2', 'inscode2', 'chain_id2', 'atomno2', 'atomna2', 'dist' ) ;
         my $kdfield2no ;
         foreach my $j ( 0 .. $#cont_fields) { $kdfield2no->{$cont_fields[$j]} = $j;}
   
# Iterate through the contacts, build all res - res contact info vectors
         my $kdcont_out  = _interface_detect_calc__calc_res_pairs({
            radius => 6.05,
            compress => 1,
            bdp_path => $bdp_path
         }) ;
   
         if (exists $kdcont_out->{error_fl})  {
            print STDERR "ERROR (bdp_id $bdp_id): $kdcont_out->{error_fl}\n" ;
            next;
         }
         my $kdcontacts_out = $kdcont_out->{contacts_fn} ;
   
   
         my $respairs ;
         my $resnames;
         my $specials ;
         my $contacts ;
         open (KDCONT, $kdcontacts_out) ;
         while (my $line = <KDCONT>) {
            chomp $line;
            if ($line =~ /^#/) {next;}
            my @f = split(/\t/, $line) ;
            my $resno_1 = $f[$kdfield2no->{resno1}].$f[$kdfield2no->{inscode1}] ; $resno_1 =~ s/ //g ;
            my $resna_1 = $f[$kdfield2no->{resna1}] ;
            my $chain_id_1 = $f[$kdfield2no->{chain_id1}];
            my $atomna_1 = $f[$kdfield2no->{atomna1}];
   
            my $resno_2 = $f[$kdfield2no->{resno2}].$f[$kdfield2no->{inscode2}] ; $resno_2 =~ s/ //g ;
            my $resna_2 = $f[$kdfield2no->{resna2}] ;
            my $chain_id_2 = $f[$kdfield2no->{chain_id2}];
            my $atomna_2 = $f[$kdfield2no->{atomna2}];
   
            if ((!defined $chain_id_1) || ($chain_id_1 eq '')) {
               $chain_id_1 = ' '; }
            if ((!defined $chain_id_2) || ($chain_id_2 eq '')) {
               $chain_id_2 = ' '; }
   
# take care of double counting
            my $sig1 = $resno_1."\n".$chain_id_1 ;
            my $sig2 = $resno_2."\n".$chain_id_2 ;
            if ($sig1 lt $sig2) { next; }
   
            my $dist = $f[$kdfield2no->{dist}] ;
            my $ressig = $sig1."\n".$sig2 ;
   
   	 my ($in1, $in2) = (0, 0);
   	 if (exists $subset_assign_flat->{$sig1}) {
   	    $in1 = 1; }
   	 if (exists $subset_assign_flat->{$sig2}) {
   	    $in2 = 1; }
   
   	 my $isdiff =1  ;
            if ((($in1 + $in2) == 2) &&
   	     ($subset_assign_flat->{$sig1} eq $subset_assign_flat->{$sig2})) {
   	    $isdiff = 0; }
   
   
   	 if ( (($in1 + $in2) > 0) && $isdiff &&
                 (substr($atomna_1, 1, 1) ne 'H') &&
                 (substr($atomna_2, 1, 1) ne 'H') &&
                 (substr($atomna_1, 1, 1) ne 'Q') &&
                 (substr($atomna_2, 1, 1) ne 'Q') &&
                 !( $resno_1 eq $resno_2 &&
   	         $chain_id_1 eq $chain_id_2) ) {
   
   	    if (!exists $resnames->{$sig1}) {
   	       $resnames->{$sig1} = $resna_1 ; }
   	    if (!exists $resnames->{$sig2}) {
   	       $resnames->{$sig2} = $resna_2 ; }
   
   
   	    if (!exists $respairs->{$ressig}) {
   	       $respairs->{$ressig}->{counts_4} = 0 ;
   	       $respairs->{$ressig}->{counts_4p5} = 0 ;
   	       $respairs->{$ressig}->{counts_5} = 0 ;
   	       $respairs->{$ressig}->{counts_5p5} = 0 ;
               }
   
   	    if ((!exists $respairs->{$ressig}->{min_dist}) ||
   	        ($dist <= $respairs->{$ressig}->{min_dist})) {
   	       $respairs->{$ressig}->{min_dist} = $dist ;
   	    }
   
               $respairs->{$ressig}->{contacts}++ ;
   	    if ($dist <= 4) {
   	       $respairs->{$ressig}->{counts_4}++ ;
   	       $respairs->{$ressig}->{counts_4p5}++ ;
   	       $respairs->{$ressig}->{counts_5}++ ;
   	       $respairs->{$ressig}->{counts_5p5}++ ;
   	    } elsif ($dist <= 4.5) {
   	       $respairs->{$ressig}->{counts_4p5}++ ;
   	       $respairs->{$ressig}->{counts_5}++ ;
   	       $respairs->{$ressig}->{counts_5p5}++ ;
   	    } elsif ($dist <= 5) {
   	       $respairs->{$ressig}->{counts_5}++ ;
   	       $respairs->{$ressig}->{counts_5p5}++ ;
   	    } elsif ($dist <= 5.5) {
   	       $respairs->{$ressig}->{counts_5p5}++ ;
   	    }
   
# internalized special_contact from itneratmoic_contacts.pm
   	    {
   	       my $poss ;
   
   	       my $hbflag1 ;
   	       my $hbflag2 ;
   	       if (exists $special_hbres->{$atomna_1}) {
   	          $hbflag1 = $special_hbres->{$atomna_1} ;
   	       } else {
   	          $hbflag1 = $special_hbres->{$resna_1}->{$atomna_1} ; }
   
   	       if (exists $special_hbres->{$atomna_2}) {
   	          $hbflag2 = $special_hbres->{$atomna_2} ;
   	       } else {
   	          $hbflag2 = $special_hbres->{$resna_2}->{$atomna_2} ; }
   
   	       my $sbflag1 ;
   	       my $sbflag2 ;
   	       if (exists $special_sbres->{$atomna_1}) {
   	          $sbflag1 = $special_sbres->{$atomna_1} ;
   	       } else {
   	          $sbflag1 = $special_sbres->{$resna_1}->{$atomna_1} ; }
   
   	       if (exists $special_sbres->{$atomna_2}) {
   	          $sbflag2 = $special_sbres->{$atomna_2} ;
   	       } else {
   	          $sbflag2 = $special_sbres->{$resna_2}->{$atomna_2} ; }
   
   
   
   	       if (($dist <= $special_thresh->{ssbond}) &&
                      ($atomna_1.$atomna_2 eq ' SG  SG ' )) {
   	          $specials->{$ressig}->{'ssbond'}++ ;
   	       } elsif (($dist <= $special_thresh->{salt}) &&
   	                (defined $sbflag1 &&
   	                 defined $sbflag2 &&
   		         $sbflag1 + $sbflag2 == 0 ) ) {
   	          $specials->{$ressig}->{'salt'}++ ;
   	       } else {
   
   		  my $thresh = $special_thresh->{hbond} ;
   		  if ($atomna_1 eq ' SG ' || $atomna_1 eq ' SD ' ||
   		      $atomna_2 eq ' SG ' || $atomna_2 eq ' SD ' ) {
   		      $thresh = $special_thresh->{hbond_s} ; }
   
                     if ($dist < $thresh) {
                        if ((defined $hbflag1 && defined $hbflag2) &&
                           (abs($hbflag1 + $hbflag2) < 2 )) {
                           $specials->{$ressig}->{'hbond'}++ ;
                        }
   	          }
   	       }
               }
   	 }
         }
   
         close(KDCONT) ;
         unlink $kdcontacts_out ;
   
# Iterate through the residue pairs
# For each new interacting atom (residue? fpd030714_0531), iterate through BDP subsets, and find ALL correspding domains, then output for import into subsets_residue.
         my $intersubset_contacts ;
   
         my $contact_counts ;
         $contact_counts->{4} = 0 ;
         $contact_counts->{4.5} = 0 ;
         $contact_counts->{5} = 0 ;
         $contact_counts->{5.5} = 0 ;
         $contact_counts->{real} = 0 ;
   
         $contact_counts->{ssbond} = 0 ;
         $contact_counts->{hbond} = 0 ;
         $contact_counts->{saltbridge} = 0 ;
   
         my $patch_residues ; #->{subset_id_index}->{res\nchain_no} = # contacts
   
         foreach my $respair (keys %{$respairs}) {
   
            my ($resno_1, $chain_id_1, $resno_2, $chain_id_2) = 
   	    split(/\n/, $respair) ;
   
            my $sig1 = $resno_1."\n".$chain_id_1 ;
            my $sig2 = $resno_2."\n".$chain_id_2 ;
   
# Go looking for subset membership of res 1 and 2.
   	 my @con_resno = ( $resno_1, $resno_2 ) ;
   	 my @con_resna = ( $resnames->{$sig1}, $resnames->{$sig2} ) ;
   	 my @con_chain_id = ( $chain_id_1, $chain_id_2 ) ;
   	 my @sigs = ( $sig1, $sig2) ;
   
#fixed BUG: make sure CATH v SCOP merging of domains doesn't preclude
# SCOP interface from being detected
# i.e. The CATH domain def may be the same, but the SCOPs are different
   	 my $diffsubset = 0;
   	 my $diffsubsets ;
   	 foreach my $l ( 0 .. $#{$subset_assign->{$sigs[0]}} ) {
   	    my $a = $subset_assign->{$sigs[0]}->[$l] ;
   	    my $a_source = $subsets_source_id[$a] ;
   	    foreach my $m ( 0 .. $#{$subset_assign->{$sigs[1]}} ) {
   	       my $b = $subset_assign->{$sigs[1]}->[$m] ;
   	       my $b_source = $subsets_source_id[$b] ;
   	       if ($a_source eq $b_source) {
   	          if ($a != $b) {
   	             $diffsubsets->{$a} = 1 ; $diffsubsets->{$b} = 1 ;
   		     $diffsubset++ ;
   	          } else {
   		     $diffsubsets->{$a} = 0 ; $diffsubsets->{$b} = 0 ;
   	          }
   	       }
   	    }
   	 }
   
   	 foreach my $l ( 0 .. $#{$subset_assign->{$sigs[0]}}) {
   	    my $a = $subset_assign->{$sigs[0]}->[$l] ;
   	    if ($diffsubsets->{$a}) {
   	       $patch_residues->{$a}->{$con_chain_id[0]}->{$con_resno[0]}->{$con_resna[0]} += $respairs->{$respair}->{contacts} ;
   	    }
   	 }
   
   	 foreach my $l ( 0 .. $#{$subset_assign->{$sigs[1]}}) {
   	    my $a = $subset_assign->{$sigs[1]}->[$l] ;
   	    if ($diffsubsets->{$a}) {
   	       $patch_residues->{$a}->{$con_chain_id[1]}->{$con_resno[1]}->{$con_resna[1]} += $respairs->{$respair}->{contacts} ;
   	    }
            }
   
            if ( ( $#{$subset_assign->{$sigs[0]}} >= 0 ) &&
   	      ( $#{$subset_assign->{$sigs[1]}} >= 0 ) ) {
   
   	    foreach my $l ( 0 .. $#{$subset_assign->{$sigs[0]}} ) {
   	    
   	       my $a = $subset_assign->{$sigs[0]}->[$l] ;
   
   	       foreach my $m ( 0 .. $#{$subset_assign->{$sigs[1]}} ) {
   
   		  my $b = $subset_assign->{$sigs[1]}->[$m] ;
   
   	          if ( ( $a != $b ) &&
   		       ( $subsets_source_id[$a] ==
   		         $subsets_source_id[$b]) ) {
   
   		     my $one_n; my $two_n ;
   		     my $one_a; my $two_a ;
   		     if ($subsets_id[$a] lt $subsets_id[$b]) {
   		        $one_n = 0 ; $two_n = 1 ;
   			$one_a = $a ; $two_a = $b ;
                        } else {
   		        $one_n = 1 ; $two_n = 0 ;
   			$one_a = $b ; $two_a = $a ;
   		     }
   		     
   		     
   		     if (!(exists $intersubset_contacts->{$one_a}->{$two_a})) {
   			$intersubset_contacts->{$one_a}->{$two_a}->{contacts}=0;
   			$intersubset_contacts->{$one_a}->{$two_a}->{counts_4}=0;
   			$intersubset_contacts->{$one_a}->{$two_a}->{counts_4p5}=0;
   			$intersubset_contacts->{$one_a}->{$two_a}->{counts_5}=0;
   			$intersubset_contacts->{$one_a}->{$two_a}->{counts_5p5}=0;
   			$intersubset_contacts->{$one_a}->{$two_a}->{hbond}=0;
   			$intersubset_contacts->{$one_a}->{$two_a}->{ssbond}=0;
   			$intersubset_contacts->{$one_a}->{$two_a}->{salt}=0;
   			$intersubset_contacts->{$one_a}->{$two_a}->{min_dist}=
   			   $respairs->{$respair}->{min_dist} ;
   		     }
   
   		     if ($con_chain_id[$one_n] eq $con_chain_id[$two_n]) {
   			$intersubset_contacts->{$one_a}->{$two_a}->{same}++ ;
   		     } else {
   			$intersubset_contacts->{$one_a}->{$two_a}->{diff}++ ;
   		     }
   
   		     $intersubset_contacts->{$one_a}->{$two_a}->{contacts} +=
   		        $respairs->{$respair}->{contacts};
   
   		     $intersubset_contacts->{$one_a}->{$two_a}->{counts_4} +=
   		        $respairs->{$respair}->{counts_4};
   		     $intersubset_contacts->{$one_a}->{$two_a}->{counts_4p5} +=
   		        $respairs->{$respair}->{counts_4p5};
   		     $intersubset_contacts->{$one_a}->{$two_a}->{counts_5} +=
   		        $respairs->{$respair}->{counts_5};
   		     $intersubset_contacts->{$one_a}->{$two_a}->{counts_5p5} +=
   		        $respairs->{$respair}->{counts_5p5};
   
   		     if (exists $specials->{$respair}) {
   		        foreach my $type (keys %{$specials->{$respair}}) {
   			   if (exists $specials->{$respair}->{$type}) {
   			      $intersubset_contacts->{$one_a}->{$two_a}->{$type}++; }
   		           my $special_fields = [
   			                 $bdp_id,
   					 $subsets_id[$one_a],
   					 $subsets_id[$two_a],
   					 $bdp_chain_no->{$con_chain_id[$one_n]},
   					 $con_chain_id[$one_n],
   					 $con_resno[$one_n],
   					 $con_resna[$one_n],
   					 $bdp_chain_no->{$con_chain_id[$two_n]},
   					 $con_chain_id[$two_n],
   					 $con_resno[$two_n],
   					 $con_resna[$two_n] ,
   					 $type
   				] ;
   		           print {$out_fh->{interface_contacts_special}}
   			      join("\t", @{$special_fields})."\n" ;
   		        }
   		     }
   
   		     my $interface_contacts_fields=[ $bdp_id,
   					 $subsets_id[$one_a],
   					 $subsets_id[$two_a],
   					 $bdp_chain_no->{$con_chain_id[$one_n]},
   					 $con_chain_id[$one_n],
   					 $con_resno[$one_n],
   					 $con_resna[$one_n],
   					 $bdp_chain_no->{$con_chain_id[$two_n]},
   					 $con_chain_id[$two_n],
   					 $con_resno[$two_n],
   					 $con_resna[$two_n] ,
   	                                 $respairs->{$respair}->{min_dist},
   	                                 $respairs->{$respair}->{contacts},
   	                                 $respairs->{$respair}->{counts_4},
   	                                 $respairs->{$respair}->{counts_4p5},
   	                                 $respairs->{$respair}->{counts_5},
   	                                 $respairs->{$respair}->{counts_5p5} ] ;
   
   		     foreach my $ohshit (0 .. $#{$interface_contacts_fields}){
   		        if (!defined $interface_contacts_fields->[$ohshit]) {
   			   print STDERR " field $ohshit undef\n" ; } }
   		     print {$out_fh->{interface_contacts}}
   		        join("\t", @{$interface_contacts_fields})."\n";
   
                     }
   	       }
   	    }
   	 }
   
         }
         close ($out_fh->{interface_contacts}) ;
         close ($out_fh->{interface_contacts_special}) ;
   
         system("gzip $out_fn->{interface_contacts}") ;
         $out_fn->{interface_contacts} .= '.gz' ;
   
         system("gzip $out_fn->{interface_contacts_special}") ;
         $out_fn->{interface_contacts_special} .= '.gz' ;
   
# Iterate through patch_residues hash and display to PATCH_RESIDUES file.
         foreach my $a ( keys %{$patch_residues} ) {
            foreach my $b ( keys %{$patch_residues->{$a}}) {
               foreach my $c ( keys %{$patch_residues->{$a}->{$b}}) {
                  foreach my $d ( keys %{$patch_residues->{$a}->{$b}->{$c}}) {
   		        my @fields = ( $bdp_id,
   			               $subsets_id[$a],
   			               $bdp_chain_no->{$b}, $b, $c, $d,
   				       $patch_residues->{$a}->{$b}->{$c}->{$d}
   			             ) ;
   	                print {$out_fh->{patch_residues}} join("\t", @fields)."\n" ; } } } }
         close ($out_fh->{patch_residues}) ;
         system("gzip $out_fn->{patch_residues}") ;
         $out_fn->{patch_residues} .= '.gz' ;
   
# OUTPUT: Display output file names of non-empty output files.
         my @output = qw/patch_residues interface_contacts interface_contacts_special/ ;
         foreach my $outf (@output) {
            if (-z $out_fn->{$outf}) {
               unlink $out_fn->{$outf} ;
            } else {
               push @outqueue, {file => $out_fn->{$outf},
                                dir  => $out_dir->{$outf}} ;
               my @outvals = ($bdp_id, $in->{dist_cutoff},
                              '\N',$out_dir->{$outf}."/".$out_fn->{$outf});
#            print {$out_meta->{$outf}} join("\t", @outvals)."\n";
               print $outf.'_tables'."\t".join("\t", @outvals)."\n";
            }
         }
   
# Iterate through the first subset partner of all contact
         foreach my $subset1 (keys %{$intersubset_contacts}) {
   
# Iterate through all interaction partners of subset1
            foreach my $subset2 (keys %{$intersubset_contacts->{$subset1}}) {
   
# Display contacting subsets for intersubset_contacts import.
# print(INTERSUBSET_CONTACTS	bdp_id	subset1	subset2)
               my $t_chainstat ;
               if (exists $intersubset_contacts->{$subset1}->{$subset2}->{same}) {
                  $t_chainstat = 'same';
                  if (exists $intersubset_contacts->{$subset1}->{$subset2}->{diff}) {
                     $t_chainstat = 'both' ;
                  }
               } else {
                  $t_chainstat = 'diff' ;
               }
   
               my @outfields = (
                  $bdp_id,
                  $subsets_id[$subset1],
                  $subsets_class[$subset1] ,
                  $subsets_id[$subset2] ,
                  $subsets_class[$subset2],
                  $intersubset_contacts->{$subset1}->{$subset2}->{contacts},
                  $in->{dist_cutoff},
                  $intersubset_contacts->{$subset1}->{$subset2}->{counts_4},
                  $intersubset_contacts->{$subset1}->{$subset2}->{counts_4p5},
                  $intersubset_contacts->{$subset1}->{$subset2}->{counts_5},
                  $intersubset_contacts->{$subset1}->{$subset2}->{counts_5p5},
                  $intersubset_contacts->{$subset1}->{$subset2}->{hbond},
                  $intersubset_contacts->{$subset1}->{$subset2}->{salt},
                  $intersubset_contacts->{$subset1}->{$subset2}->{ssbond},
                  $t_chainstat
               ) ;
   
#            print {$out_fh->{intersub_cont}} join("\t", @outfields)."\n" ;
               print "intersubset_contacts\t".join("\t", @outfields)."\n" ;
   
            }
   
         }
   
         if ($#outqueue >= $movethreshold) {
            foreach my $j ( 0 .. $#outqueue) {
   	    pibase::safe_move($outqueue[$j]->{file}, $outqueue[$j]->{dir}) ; }
            @outqueue = () ;
         }
   
         print STDERR " X\n" ;
      }
      
      foreach my $j ( 0 .. $#outqueue) {
         pibase::safe_move($outqueue[$j]->{file}, $outqueue[$j]->{dir}) ; }
   
# commented out fpd20080509_1202 - switch over to STDOUT display, post-grep
#   close($out_fh->{intersub_cont}) ;
#   close($out_fh->{patch_residues_meta}) ;
#   close($out_fh->{interface_contacts_meta}) ;
   
      print STDERR "\n" ;
   
   }

}


=head2 _interface_detect_calc__calc_res_pairs().

   Title:       _interface_detect_calc__calc_res_pairs()
   Function:    Calculates residue contacts in a given pdb file.
   Args:        $_->{radius} - upper distance limit on inter-atomic contacts calculation [default=6.6 Ang]
                $_->{compress} - compression flag
                $_->{bdp_path} - bdp file path

   Return:      $_->{contacts_fn} - kdcontacts output file
                $_->{fields} - kdcontacts file field names
                $_->{field2no}->{field} = i - hash mapping kdcontacts field names to field number

   Files IN:    PDB file ($_->{bdp_path})
   Files OUT:   kdcontacts output file ($_->{contacts_fn})

=cut

sub _interface_detect_calc__calc_res_pairs {
   my $params = shift ;

   my $kdcontacts_radius = $params->{radius} || "6.6" ;
   my $compress_fl = $params->{compress} || 1 ;

   my $binaries = pibase::locate_binaries() ;
   my $kdcontacts_bin = $binaries->{'kdcontacts'}." $kdcontacts_radius" ;
   my $altloc_check = $binaries->{'altloc_check'} ;
   my $altloc_filter = $binaries->{'altloc_filter'} ;

   my $bdp_file_path = $params->{bdp_path} ;
   my $localbdp = $bdp_file_path ; $localbdp =~ s/.*\/// ;

   if ($bdp_file_path =~ /gz$/) {
      my $tcom = "$binaries->{'zcat'} $bdp_file_path > $localbdp" ;
      system($tcom) ;
   } else {
      pibase::safe_copy($bdp_file_path, $localbdp) ;
   }
   
   if (!-e $localbdp) {
      return {error_fl => "ERROR: pdb file access error - couldnt copy locally"} ;
   }

   my $host = hostname() ;
   my ($fh, $kdcontacts_out) =
      tempfile("kdcontacts.$host.XXXXXX", SUFFIX => ".out"); close($fh) ;
   my ($fh2, $kdcontacts_err) =
      tempfile("kdcontacts.$host.XXXXXX", SUFFIX => ".err"); close($fh2) ;

# Check if the pdb file contains altloc identifiers. If so, first filter with altloc_filter.pl
   my $altloc_fl = `$altloc_check < $localbdp` ; chomp $altloc_fl ;
   my $tcom ;
   if ($altloc_fl) {
      $tcom = "$altloc_filter $localbdp" ;
   } else {
      $tcom = "cat $localbdp" ;
   }

# system(cat <bdp_file_path> | kdcontacts_bin 2><kdcontacts_err> ><kdcontacts_out>)
   $tcom .= " | $kdcontacts_bin 2>$kdcontacts_err >$kdcontacts_out" ;
   system($tcom) ;

# If the kdcontacts_out file is empty or the _err file exists, display to STDERR and remove the files.
   if (!(-e $kdcontacts_out) || (! -z $kdcontacts_err))  {
      unlink $kdcontacts_err, $kdcontacts_out;
      return {error_fl => "ERROR: $bdp_file_path kdcontacts execution error"} ;

# Otherwise remove the _err file
   } else {

      unlink $kdcontacts_err ;

   }

# Read in kdcontacts output
#  read in kdcontacts output into an array

   my @cont_fields = ('resna1', 'resno1', 'inscode1', 'chain_id1', 'atomno1', 'atomna1', 'resna2', 'resno2', 'inscode2', 'chain_id2', 'atomno2', 'atomna2', 'dist' ) ;
   my $kdfield2no ;
   foreach my $j ( 0 .. $#cont_fields) { $kdfield2no->{$cont_fields[$j]} = $j;}

   my $results = {
      contacts_fn => $kdcontacts_out,
      fields => \@cont_fields,
      field2no => $kdfield2no
   };

   {
      my ($newdev, $newino) = stat($localbdp) ;
      my ($olddev, $oldino) = stat($bdp_file_path) ;
      if ($newdev != $olddev || $newino != $oldino) {
         unlink $localbdp ;
      }
   }

   return $results ;

}



################### OLD INTERFACE CONTACT CODE
#BUG potential
#
#* fpd040224_0711 - ignoring residue contacts which form part of the same domain in ANY assignment. i.e. valid SCOP interfacial contact, but intra-domain CATH contact; yet also ignored for SCOP interface.
#
#most likely doesnt affect interface_contacts tables; JUST patch_residues
#
#see in CHANGE list
#
#* subset residues tables have a heading space; trace back to source
#

sub _interface_detect__calc_interatomic_contacts {

   use pibase::interatomic_contacts qw/contacts_select contacts_select_inter special_contact raw_contacts_select/ ;

   my $params = shift ;

   my $kdcontacts_radius = $params->{radius} || "6.6" ;
   my $compress_fl = $params->{compress} || 1 ;

   my $binaries = pibase::locate_binaries() ;
   my $kdcontacts_bin = $binaries->{'kdcontacts'}." $kdcontacts_radius" ;
   my $altloc_check = $binaries->{'altloc_check'} ;
   my $altloc_filter = $binaries->{'altloc_filter'} ;

   my $bdp_file_path = $params->{bdp_path} ;
   my $localbdp = $bdp_file_path ; $localbdp =~ s/.*\/// ;

   pibase::safe_copy($bdp_file_path, $localbdp) ;
   
   if (!-e $localbdp) {
      return {error_fl => "ERROR: pdb file access error - couldnt copy locally"} ;
   }

   my $host = hostname() ;
   my ($fh, $kdcontacts_out) =
      tempfile("kdcontacts.$host.XXXXXX", SUFFIX => ".out"); close($fh) ;
   my ($fh2, $kdcontacts_err) =
      tempfile("kdcontacts.$host.XXXXXX", SUFFIX => ".err"); close($fh2) ;

# Check if the pdb file contains altloc identifiers. If so, first filter with altloc_filter.pl

   my $altloc_fl = `$altloc_check < $localbdp` ; chomp $altloc_fl ;
   my $tcom ;
   if ($altloc_fl) {
      $tcom = "$altloc_filter $localbdp" ;
   } else {
      $tcom = "cat $localbdp" ;
   }

# system(cat <bdp_file_path> | kdcontacts_bin 2><kdcontacts_err> ><kdcontacts_out>)

   $tcom .= " | $kdcontacts_bin 2>$kdcontacts_err >$kdcontacts_out" ;
   system($tcom) ;

# If the kdcontacts_out file is empty or the _err file exists, display to STDERR and remove the files.

   if (!(-e $kdcontacts_out) || (! -z $kdcontacts_err))  {
      unlink $kdcontacts_err, $kdcontacts_out;
      return {error_fl => "ERROR: $bdp_file_path kdcontacts execution error"} ;

# Otherwise remove the _err file

   } else {

      unlink $kdcontacts_err ; }

# Read in kdcontacts output

# read in kdcontcts output into an array

   my $contacts ;
   open (KDCONT, $kdcontacts_out) ;
   while (my $line = <KDCONT>) {
      chomp $line;
      if ($line !~ /^#/) {
         my @fields = split(/\t/, $line) ;
         push @{$contacts}, \@fields ;
      }
   }
   close(KDCONT) ;

   unlink $kdcontacts_out ;
   unlink $localbdp ;

   return {contacts => $contacts };

}


=head2 cluster_scop_interfaces()

   Title:       cluster_scop_interfacesj()
   Function:    Clusters SCOP-SCOP interfaces using ASTRAL ASTEROIDS alignments.
                 and imports to PIBASE if specified
   Args:        ->{pibase_specs} = optional
                ->{import_fl} = 1 if to be imported into PIBASE
   Returns:     Prints out cluster membership to
                 $pibase_specs->{buildfiles}->{scop_interface_clusters}) ;

=cut

sub cluster_scop_interfaces {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }

   my $run_options ;
   $run_options->{DEBUGMERGER} = 0 ;
   $run_options->{DEBUG} = 0 ;
   $run_options->{DEBUGFAM1} = 0 ;
   $run_options->{DEBUGFAM2} = 0 ;
   $run_options->{DEBUGDOM} = 0 ;
   $run_options->{DEBUGDOMFAM} = 0 ;

   my $thresh ;
   $thresh->{aligned_intres} = 0.7 ;
   $thresh->{aligned_contacts} = 0.7 ;

   print STDERR "Load ASTRAL fasta headers (gdseq): " ;
   my $astral = pibase::ASTRAL::load_astral_headers({
      fn => $pibase_specs->{astral}->{gd_seq}
   }) ;
   print STDERR "X\n" ;

   print STDERR "Load ASTRAL raf: " ;
   my $raf = pibase::ASTRAL::raf_preload({
      fn => $pibase_specs->{astral}->{raf}
   }) ;
   print STDERR "X\n" ;

# 0. load sid1,sid2 pairs

   print STDERR "Load PIBASE data: " ;
   my $pb = _cluster_interfaces_pibase_preload({astral => $astral}) ;
   print STDERR " X\n" ;

   my $goal ;
   if (exists $in->{in_fn}) {
      open(INF, $in->{in_fn}) ;
      while (my $line = <INF>) {
         if ($line =~ /^#/) {next;}
         chomp $line;
         my ($fam1, $fam2) = split(/\t/, $line) ;
         my $fam12 ;
         if ($fam1 lt $fam2) {
            $fam12 = $fam1."\t".$fam2;
         } else {
            $fam12 = $fam2."\t".$fam1;
         }
         $goal->{fampairs}->{$fam12}++ ;
         $goal->{fam_2_pairs}->{$fam1}->{$fam12}++ ;
         $goal->{fam_2_pairs}->{$fam2}->{$fam12}++ ;
      }
      close(INF) ;
   } else {
      my ($dbh) = pibase::connect_pibase() ;
      my ($fam1, $fam2) = pibase::mysql_fetchcols($dbh,
         "SELECT distinct class_1, class_2 FROM intersubset_contacts as a,".
         "bdp_files as b WHERE subset_id_1 LIKE \"\%SCOP\%\" AND ".
         "a.bdp_id = b.bdp_id and b.raw_pdb = 1 AND ".
         "class_1 regexp \"(a|b|c|d|e|f|g)\" AND ".
         "class_2 regexp \"(a|b|c|d|e|f|g)\""
         ) ;
      foreach my $j ( 0 .. $#{$fam1}) {
         my $fam12 ;
         if ($fam1->[$j] lt $fam2->[$j]) {
            $fam12 = $fam1->[$j]."\t".$fam2->[$j];
         } else {
            $fam12 = $fam2->[$j]."\t".$fam1->[$j];
         }
         $goal->{fampairs}->{$fam12}++ ;
         $goal->{fam_2_pairs}->{$fam1->[$j]}->{$fam12}++ ;
         $goal->{fam_2_pairs}->{$fam2->[$j]}->{$fam12}++ ;
      }
   }


# 1. load ASTRAL seq clusters

   pibase::ASTRAL::load_astral_clusters({
      pibase_specs => $pibase_specs,
      out => $astral
   }) ;

# 1.5. read in class pairs to cluster from STDIN - for cluster move

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1) {
      print "* cluster_interfaces() ".localtime() if (!exists $in->{quiet_fl});

      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{cluster_interfaces_in}, $temp_fn->{cluster_interfaces_in}) =
         tempfile("splits_cluster_interfaces_input.XXXXX");
      ($temp_fh->{cluster_interfaces_out}, $temp_fn->{cluster_interfaces_out}) =
         tempfile("splits_cluster_interfaces_SGEout_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{cluster_interfaces_out}) ;
      ($temp_fh->{cluster_interfaces_err}, $temp_fn->{cluster_interfaces_err}) =
         tempfile("splits_cluster_interfaces_SGEerr_XXXXX", SUFFIX => '.pibase');
         close($temp_fh->{cluster_interfaces_err}) ;

      foreach my $fam12 (sort keys %{$goal->{fampairs}}) {
         my ($fam1, $fam2) = split(/\t/, $fam12) ;
         print {$temp_fh->{cluster_interfaces_in}}
            join("\t", $fam1, $fam2)."\n" ;
      }
      close($temp_fh->{cluster_interfaces_in}) ;

      my $split_dir = tempdir("splits_cluster_interfaces.XXXXX") ;
      my $splits = pibase::SGE::_clust_split_ins({
         fn => $temp_fn->{cluster_interfaces_in},
         dir => $split_dir,
         numjobs => $pibase_specs->{SGE}->{numjobs}
      });

      my ($perlscript_fh, $perlscript_fn) =
            tempfile("pb.cluster_interfaces.XXXXX",
                     SUFFIX => ".pbi.pl") ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase::calc::interfaces qw/cluster_scop_interfaces/ ;

main() ;

sub main {

         pibase::calc::interfaces::cluster_scop_interfaces({
            cluster_fl => 0,
            import_fl => 0,
            in_fn => \$ARGV[0],
         }) ;

}\n" ;
      close($perlscript_fh) ;

      my ($sgescript_fh, $sgescript_fn) =
         tempfile("pb.cluster_scop_interfaces.XXXXX", SUFFIX=>".SGE.sh");
      my $sge_outdir = tempdir("SGEOUT.cluster_scop_interfaces.XXXXX");

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
         out_fn => $temp_fn->{cluster_interfaces_out},
         err_fn => $temp_fn->{cluster_interfaces_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      open($temp_fh->{cluster_interfaces_out},
           $temp_fn->{cluster_interfaces_out}) ;
      open(REALOUTF,">".$pibase_specs->{buildfiles}->{scop_interface_clusters}) ;
      while (my $line = readline($temp_fh->{cluster_interfaces_out})) {
         if ($line =~ /^\#/) {next;}
         print REALOUTF $line ;
      }
      close($temp_fh->{cluster_interfaces_out}) ;
      close(REALOUTF) ;

   } else {

   # make list of interface instances to be studied for each family, family pair
   # make hash of binding sites and interfaces
   # interface: fam1.fam2->{pbnum}
   # bindingiste: fam1->{pbnum.{0|1}}
   
      my $dataset ;
      foreach my $fam1 (keys %{$goal->{fam_2_pairs}}) {
#         print STDERR "on fam1 $fam1\n" ;
         foreach my $fam12 (keys %{$goal->{fam_2_pairs}->{$fam1}}) {
#            print STDERR "   on fam12 $fam12\n" ;
            if (!exists $pb->{fampairs}->{$fam12}) {
               print STDERR "WARNING: no pibase entry for $fam12 (line ".
                            __LINE__.")" ;
               next;
            }

            if (!exists $pb->{fampairs}->{$fam12}) {
            }

            push @{$dataset->{perfam_ints}->{$fam1}},
               @{$pb->{fampairs}->{$fam12}} ;
   
            foreach my $j ( @{$pb->{fampairs}->{$fam12}}) {
               $dataset->{fam_osid}->{$pb->{class1}->[$j]}->{$pb->{osid1}->[$j]}++;
               $dataset->{fam_osid}->{$pb->{class2}->[$j]}->{$pb->{osid2}->[$j]}++;
   
               $dataset->{bs}->{$pb->{class1}->[$j]}->{$j.".".0}++ ;
               $dataset->{bs}->{$pb->{class2}->[$j]}->{$j.".".1}++ ;
   
               $dataset->{interf}->{$pb->{class1}->[$j]}->{$j}++ ;
               $dataset->{interf}->{$pb->{class2}->[$j]}->{$j}++ ;
            }
         }
      }
   
   # 1.75 iterate over families: cluster binding sites
   
      my $clust_bs ; # global bs cluster assignment
      print STDERR "Clustering binding sites:\n" ;
      foreach my $fam (sort keys %{$dataset->{bs}}) {
   
         if ($run_options->{DEBUG} &&
             ($fam ne $run_options->{DEBUGFAM1} &&
              $fam ne $run_options->{DEBUGFAM2})) {
            next;}
   
         print STDERR "on family $fam\n" ;

#         print STDERR "   load_asteroids_aln() on ".join(" ", sort keys %{$dataset->{fam_osid}->{$fam}})."\n" ;
   
         my $fam_aln = pibase::ASTRAL::load_asteroids_aln({
            doms => $dataset->{fam_osid}->{$fam},
            aln_fn => $pibase_specs->{asteroids}->{fam_aln}.'/'.$fam.
               '.fasta_aln' ,
            seq_fn => $pibase_specs->{asteroids}->{fam_seq}.'/'.$fam.'.fa' ,
            raf => $raf,
            gdseqh => $astral->{gdseqh},
            seqclcont100 => $astral->{seqcl2cont}->{100},
            seqcl100 => $astral->{seqcl}->{100},
            allchains => $pb->{pdbchains}
         }) ;
   
   # load all binding site residues for family $fam
   
         my $intres ;
         foreach my $j ( sort keys %{$dataset->{interf}->{$fam}} ) {
            $intres->{$j} = _cluster_interfaces_load_interface_contacts({
               sid1 => $pb->{sid1}->[$j],
               sid2 => $pb->{sid2}->[$j],
               fn => $pb->{bdp2contactsfn}->{$pb->{bdp_id}->[$j]}
            }) ;
         }
   
         my @fam_bs = keys %{$dataset->{bs}->{$fam}} ;
   
         my ($bsclust, $bs2clust) ;
         foreach my $j ( 0 .. $#fam_bs) {
            $bsclust->{$j} = [$fam_bs[$j]] ;
            $bs2clust->{$fam_bs[$j]} = $j ;
         }
   
         my @seqids = (sort {$b <=> $a} keys %{$pibase_specs->{astral}->{seqcl}});
         push @seqids, 0 ;
         my $clusround = 1 ;
         my $notthisone ;
   
         foreach my $seqid (@seqids) {
            my @mergelist ;
            my @bsclusts = sort keys %{$bsclust} ;
   
            print STDERR " round $clusround: fresh active clusters: ".
               join(",", @bsclusts)."\n" if $run_options->{DEBUGMERGER} ; 
   
            my ($csid, $cosid) ;
   
            foreach my $tj (0 .. ($#bsclusts - 1)) {
               print STDERR "round $clusround: clust1 = $bsclusts[$tj]\n"
                  if $run_options->{DEBUGMERGER};
   
               my ($j, $j_side) = split(/\./, $bsclust->{$bsclusts[$tj]}->[0]) ;
   
               my @tjs = ($pb->{sid1}->[$j], $pb->{sid2}->[$j]) ;
               $csid->[0] = $tjs[$j_side] ;
   
               my @tjos = ($pb->{osid1}->[$j], $pb->{osid2}->[$j]) ;
               $cosid->[0] = $tjos[$j_side] ;
   
               if (exists $notthisone->{$cosid->[0]}) {next;}
   
               foreach my $tk (($tj+1) .. $#bsclusts) {
                  my ($k, $k_side) = split(/\./, $bsclust->{$bsclusts[$tk]}->[0]) ;
   
                  my @tks = ($pb->{sid1}->[$k], $pb->{sid2}->[$k]) ;
                  $csid->[1] = $tks[$k_side] ;
   
                  my @tkos = ($pb->{osid1}->[$k], $pb->{osid2}->[$k]) ;
                  $cosid->[1] = $tkos[$k_side] ;
   
                  print STDERR "round $clusround: clust2 = $bsclusts[$tk]"
                     if $run_options->{DEBUGMERGER};
   
                  print STDERR "round $clusround: $j vs $k\n"
                     if $run_options->{DEBUGMERGER};
   
                  if (exists $notthisone->{$cosid->[1]}) {next;}
                  
                  if ($seqid > 0) {
                     my @cerrors = () ;
                     if (!exists $astral->{seqcl}->{$seqid}->{$cosid->[0]}) {
                        $notthisone->{$cosid->[0]}++ ;
                        push @cerrors, "$cosid->[0] not found in ASTRAL cluster";}
   
                     if (!exists $astral->{seqcl}->{$seqid}->{$cosid->[1]}) {
                        $notthisone->{$cosid->[1]}++ ;
                        push @cerrors, "$cosid->[1] not found in ASTRAL cluster";}
   
                     if ($#cerrors >= 0 ) {
                        print STDERR "ERROR: aborting comparison of ".
                           "$pb->{sid1}->[$j] -- $pb->{sid2}->[$j] ($j_side) to ".
                           "$pb->{sid1}->[$k] -- $pb->{sid2}->[$k] ($k_side) ".
                           "($seqid seqid) dt:".join(",", @cerrors)."\n" ; next; }
   
                     if ($astral->{seqcl}->{$seqid}->{$cosid->[0]} !=
                         $astral->{seqcl}->{$seqid}->{$cosid->[1]}) { next; }
   
                     @cerrors = ();
   
                     if (!exists $fam_aln->{resno2pos}->{$cosid->[0]}) {
                        $notthisone->{$cosid->[0]}++ ;
                        push @cerrors, "$cosid->[0] not in fam alignment";}
   
                     if (!exists $fam_aln->{resno2pos}->{$cosid->[1]}) {
                        $notthisone->{$cosid->[1]}++ ;
                        push @cerrors, "$cosid->[1] not in fam alignment";}
   
                     if ($#cerrors >= 0 ) {
                        print STDERR "ERROR: aborting comparison of ".
                           "$pb->{sid1}->[$j] -- $pb->{sid2}->[$j] ($j_side) to ".
                           "$pb->{sid1}->[$k] -- $pb->{sid2}->[$k] ($k_side) ".
                           "($seqid seqid) dt:".join(",", @cerrors)."\n" ; next; }
                  }
   
   
                  my $bscomp ;
                  $bscomp = _cluster_interfaces_compare_residue_sets({
                     aln => $fam_aln,
                     sid1 => $cosid->[0],
                     sid2 => $cosid->[1],
                     res1 => $intres->{$j}->{intres}->{$csid->[0]},
                     res2 => $intres->{$k}->{intres}->{$csid->[1]},
                  }) ;
                  
                  if ($bscomp < $thresh->{aligned_intres}) { next; }
   
                  print STDERR "BS COMP $j v $k: ($bscomp)\n"
                     if $run_options->{DEBUGMERGER} ;
   
                  print STDERR "merge $j.$j_side with $k.$k_side\n"
                     if $run_options->{DEBUGMERGER};
                  push @mergelist, [$j.'.'.$j_side, $k.'.'.$k_side] ;
   
               }
            }
            #do mergelist
   
   
            my $currloc ;
            foreach my $cllabel (keys %{$bsclust}) {
               $currloc->{$cllabel} = $cllabel ;  }
   
            foreach my $j ( 0 .. $#mergelist) {
               print STDERR "round $clusround MERGE #$j: $mergelist[$j]->[0] ".
                  "to $mergelist[$j]->[1]\n" if $run_options->{DEBUGMERGER};
   
               my $bs1 = $mergelist[$j]->[0] ;
               my $bs2 = $mergelist[$j]->[1] ;
   
               my $curr1 = $bs2clust->{$bs1} ;
               my $curr2 = $bs2clust->{$bs2} ;
               if ($curr1 == $curr2) {next;}
   
               print STDERR "  currents: $curr1 and $curr2\n"
                  if $run_options->{DEBUGMERGER};
   
               foreach my $oldbs2 ( @{$bsclust->{$curr2}}) {
                  $bs2clust->{$oldbs2} = $curr1 ;
               }
   
               push @{$bsclust->{$curr1}}, @{$bsclust->{$curr2}} ;
               delete $bsclust->{$curr2} ;
   
               if ($run_options->{DEBUGMERGER}) {
                  print STDERR " merging cluster $curr2 with $curr1\n" ;
                  print STDERR "set $curr1 now ".join(",", @{$bsclust->{$curr1}})."\n";
                  print STDERR "   Deleted set $curr2\n" ;
   
                  print STDERR "\tactive clusters: ".
                     join(",", sort {$a <=> $b} keys %{$bsclust})."\n" ;
               }
            }
            $clusround++ ;
         }
   
   # move final cluster memberships to global clust_bs hash
         my $nn = 1 ;
         foreach my $label (keys %{$bsclust}) {
            my $nm = 1 ;
            foreach my $bs (@{$bsclust->{$label}}) {
               my ($int, $side) = split(/\./, $bs) ;
               $clust_bs->{$bs} = "$fam.$nn" ;
               my $bs_cluster = $fam."_".$nn ;
#               my $cluster_no = $nn ;
#               my $member_no = $nm ;
# NOTE: not a valid global measure unless everything run together;
#   bs clustering is only done locally for those members of a family
#   that are involved in at least one of the goal family pairs.
#               my @outvals = ("BSCLUSTER", $sidinquestion, $fam
#                              $pb->{sid1}->[$int], $pb->{sid2}->[$int],
#                              $bs_cluster,
#                              $nn,
#                              $pb->{class1}->[$int], $pb->{class2}->[$int],
#                              $cluster_no,
#                              $member_no) ;
#               print join("\t", @outvals)."\n" ;
               $nm++ ;
            }
            $nn++ ;
         }
      }
   
   # 2. iterate over family pairs: cluster interfaces
   
      print STDERR "Clustering interfaces:\n" ;
      foreach my $fam12 (sort keys %{$goal->{fampairs}}) {
         my ($fam1, $fam2) = split(/\t/, $fam12) ;
         if ($run_options->{DEBUG} &&
             ($fam1 ne $run_options->{DEBUGFAM1} ||
              $fam2 ne $run_options->{DEBUGFAM2})) {
            next;}
   
         print STDERR "now on $fam1 -- $fam2\n";
   
# 2.1 load fam1 and fam2 alignments
   
         my $fam1_aln = pibase::ASTRAL::load_asteroids_aln({
            doms => $pb->{fampairs_single_osid}->{$fam12}->{$fam1},
            aln_fn => $pibase_specs->{asteroids}->{fam_aln}.
               '/'.$fam1.'.fasta_aln' ,
            seq_fn => $pibase_specs->{asteroids}->{fam_seq}.
               '/'.$fam1.'.fa' ,
            raf => $raf,
            gdseqh => $astral->{gdseqh},
            seqclcont100 => $astral->{seqcl2cont}->{100},
            seqcl100 => $astral->{seqcl}->{100},
            allchains => $pb->{pdbchains}
         }) ;
   
         my $fam2_aln = pibase::ASTRAL::load_asteroids_aln({
            doms => $pb->{fampairs_single_osid}->{$fam12}->{$fam2},
            aln_fn => $pibase_specs->{asteroids}->{fam_aln}.
               '/'.$fam2.'.fasta_aln' ,
            seq_fn => $pibase_specs->{asteroids}->{fam_seq}.
               '/'.$fam2.'.fa' ,
            raf => $raf,
            gdseqh => $astral->{gdseqh},
            seqclcont100 => $astral->{seqcl2cont}->{100},
            seqcl100 => $astral->{seqcl}->{100},
            allchains => $pb->{pdbchains}
         }) ;
   
         if ($run_options->{DEBUG} && $run_options->{DEBUGDOM}) {
            my $famnum = $run_options->{DEBUGDOMFAM} ;
            my $famx_aln = $fam1_aln ;
            if ($famnum == 2) {$famx_aln = $fam2_aln;}
   
            my @doms = sort keys %{$famx_aln->{resno2pos}} ;
            print STDERR "fam $famnum has ".join(",", @doms)."\n" ;
            my $cdom = $run_options->{DEBUGDOM} ;
            if (exists $famx_aln->{resno2pos}->{$cdom}) {
               foreach my $res (sort keys %{$famx_aln->{resno2pos}->{$cdom}}) {
                  print STDERR "$res has a pos of ".$famx_aln->{resno2pos}->{$cdom}->{$res}."\n" ;
               }
            } else {
               print STDERR " WARNING: $cdom found in family $famnum aln\n" ;
            }
         }
   
   
   # 2.2 load interface contacts for all fam12 interfaces
   
         my $intconts ;
         foreach my $j ( @{$pb->{fampairs}->{$fam12}} ) {
            $intconts->{$j} = _cluster_interfaces_load_interface_contacts({
               sid1 => $pb->{sid1}->[$j],
               sid2 => $pb->{sid2}->[$j],
               fn => $pb->{bdp2contactsfn}->{$pb->{bdp_id}->[$j]}
            }) ;
         }
   
   # 2.3 do comparisons
   
         my ($iclust, $i2clust, $i2clust_order) ;
         foreach my $j ( 0.. $#{$pb->{fampairs}->{$fam12}}) {
            $iclust->{$j} = [$pb->{fampairs}->{$fam12}->[$j]] ;
            $i2clust->{$pb->{fampairs}->{$fam12}->[$j]} = $j ;
            $i2clust_order->{$pb->{fampairs}->{$fam12}->[$j]} = 0 ;
   #0 = forward, 1 = backward
         }
   
         my @seqids = (sort {$b <=> $a}
            keys %{$pibase_specs->{astral}->{seqcl}});

         push @seqids, 0 ;
         my $clusround = 1 ;
   
         my $notthisone ;
         foreach my $seqid (@seqids) {
            my @mergelist ;
            my @iclusts = sort {$a <=> $b} keys %{$iclust} ;
   
            print STDERR " round $clusround: fresh active clusters: ".
               join(",",@iclusts)."\n" if $run_options->{DEBUGMERGER} ; 
   
            my ($csid1, $csid2) ;
            my ($cosid1, $cosid2) ;
            my ($curbs1, $curbs2) ;
   
            foreach my $tj (0 .. ($#iclusts - 1)) {
               print STDERR "round $clusround: clust1 = $iclusts[$tj]\n"
                  if $run_options->{DEBUGMERGER};
   
               my $j = $iclust->{$iclusts[$tj]}->[0] ;
   
               if ($pb->{revfl}->[$j] == 1) {
                  $curbs2->[0] = $j.'.0' ; $curbs1->[0] = $j.'.1' ;
   
                  ($csid2->[0], $csid1->[0]) = ($pb->{sid1}->[$j],
                                             $pb->{sid2}->[$j]) ;
   
                  ($cosid2->[0], $cosid1->[0]) = ($pb->{osid1}->[$j],
                                                $pb->{osid2}->[$j]) ;
               } else {
                  $curbs1->[0] = $j.'.0' ; $curbs2->[0] = $j.'.1' ;
   
                  ($csid1->[0], $csid2->[0]) = ($pb->{sid1}->[$j],
                                                $pb->{sid2}->[$j]) ;
   
                  ($cosid1->[0], $cosid2->[0]) = ($pb->{osid1}->[$j],
                                                   $pb->{osid2}->[$j]) ;
               }
   
               if (exists $notthisone->{$cosid1->[0]} ||
                   exists $notthisone->{$cosid2->[0]}) {next;}
   
               foreach my $tk (($tj+1) .. $#iclusts) {
                  my $k = $iclust->{$iclusts[$tk]}->[0] ;
   
                  print STDERR "round $clusround: clust2 = $iclusts[$tk]"
                     if $run_options->{DEBUGMERGER};
   
                  print STDERR "round $clusround: $j vs $k\n"
                     if $run_options->{DEBUGMERGER};
   
                  if ($pb->{revfl}->[$k] == 1) {
                     $curbs2->[1] = $k.'.0' ; $curbs1->[1] = $k.'.1' ;
                     ($csid2->[1], $csid1->[1]) = ($pb->{sid1}->[$k],
                                                   $pb->{sid2}->[$k]);
   
                     ($cosid2->[1], $cosid1->[1])=($pb->{osid1}->[$k],
                                                   $pb->{osid2}->[$k]) ;
                  } else {
                     $curbs1->[1] = $k.'.0' ; $curbs2->[1] = $k.'.1' ;
                     ($csid1->[1], $csid2->[1]) = ($pb->{sid1}->[$k],
                                                   $pb->{sid2}->[$k]);
   
                     ($cosid1->[1], $cosid2->[1])=($pb->{osid1}->[$k],
                                                   $pb->{osid2}->[$k]) ;
                  }
                  
                  if (exists $notthisone->{$cosid1->[1]} ||
                      exists $notthisone->{$cosid2->[1]}) {next;}
                  
   #               if (DEBUGMERGER) {
   #                  my @tt = ($astral->{seqcl}->{$seqid}->{$cosid1->[0]},
   #                            $astral->{seqcl}->{$seqid}->{$cosid1->[1]},
   #                            $astral->{seqcl}->{$seqid}->{$cosid2->[0]},
   #                            $astral->{seqcl}->{$seqid}->{$cosid2->[1]}) ;
   #
   #                  print STDERR "check seqid $seqid: compare ".
   #                     "$j ( $pb->{osid1}->[$j] - $pb->{osid2}->[$j] ) to ".
   #                     "$k ( $pb->{osid1}->[$k] - $pb->{osid2}->[$k] )" ;
   #
   #                  if ($seqid > 0 ) {
   #                     print STDERR "\tclusters: ".join(", ", @tt);}
   #                  print STDERR "\n" ;
   #                  print STDERR "comparing set2: ".join(keys %{$intconts->[$j]->{intres}->{$pb->{sid1}->[$j]}})."\n" ;
   #                  print STDERR "FAM1 is $fam1\n" ;
   #                  print STDERR "   $csid1->[0] ($cosid1->[0]) vs ".
   #                               "$csid1->[1] ($cosid1->[1])\n" ;
   #               }
   
                  if ($seqid > 0) {
                     my @cerrors = () ;
                     if (!exists $astral->{seqcl}->{$seqid}->{$cosid1->[0]}) {
                        $notthisone->{$cosid1->[0]}++ ;
                        push @cerrors, "$cosid1->[0] not found in ASTRAL cluster";}
   
                     if (!exists $astral->{seqcl}->{$seqid}->{$cosid1->[1]}) {
                        $notthisone->{$cosid1->[1]}++ ;
                        push @cerrors, "$cosid1->[1] not found in ASTRAL cluster";}
   
                     if (!exists $astral->{seqcl}->{$seqid}->{$cosid2->[0]}) {
                        $notthisone->{$cosid2->[0]}++ ;
                        push @cerrors, "$cosid2->[0] not found in ASTRAL cluster";}
   
                     if (!exists $astral->{seqcl}->{$seqid}->{$cosid2->[1]}) {
                        $notthisone->{$cosid2->[1]}++ ;
                        push @cerrors,
                           "$cosid2->[1] not found in ASTRAL cluster";}
   
                     if ($#cerrors >= 0 ) {
                        print STDERR "ERROR: aborting comparison of ".
                           "$cosid1->[0] - $cosid2->[0] to ".
                           "$cosid1->[1] - $cosid2->[1] ".
                           "($seqid seqid) dt:".join(",", @cerrors)."\n" ;
                           next; }
   
                     if (($astral->{seqcl}->{$seqid}->{$cosid1->[0]} !=
                          $astral->{seqcl}->{$seqid}->{$cosid1->[1]}) ||
                         ($astral->{seqcl}->{$seqid}->{$cosid2->[0]} !=
                          $astral->{seqcl}->{$seqid}->{$cosid2->[1]})) { next; }
   
                     @cerrors = ();
   
                     if (!exists $fam1_aln->{resno2pos}->{$cosid1->[0]}) {
                        $notthisone->{$cosid1->[0]}++ ;
                        push @cerrors, "$cosid1->[0] not in fam1 alignment";}
   
                     if (!exists $fam1_aln->{resno2pos}->{$cosid1->[1]}) {
                        $notthisone->{$cosid1->[1]}++ ;
                        push @cerrors, "$cosid1->[1] not in fam1 alignment";}
   
                     if (!exists $fam2_aln->{resno2pos}->{$cosid2->[0]}) {
                        $notthisone->{$cosid2->[0]}++ ;
                        push @cerrors, "$cosid2->[0] not in fam2 alignment";}
   
                     if (!exists $fam2_aln->{resno2pos}->{$cosid2->[1]}) {
                        $notthisone->{$cosid2->[1]}++ ;
                        push @cerrors, "$cosid2->[1] not in fam2 alignment";}
   
                     if ($#cerrors >= 0 ) {
                        print STDERR "ERROR: aborting comparison of ".
                           "$cosid1->[0] - $cosid2->[0] to ".
                           "$cosid1->[1] - $cosid2->[1] ".
                           "($seqid seqid) dt:".join(",", @cerrors)."\n" ;
                        next; }
                  }
   
   
                  if (!($clust_bs->{$curbs1->[0]} eq $clust_bs->{$curbs1->[1]} &&
                      $clust_bs->{$curbs2->[0]} eq $clust_bs->{$curbs2->[1]}) &&
                     !($clust_bs->{$curbs1->[0]} eq $clust_bs->{$curbs2->[1]} &&
                      $clust_bs->{$curbs2->[0]} eq $clust_bs->{$curbs1->[1]})) {
                     next;}
   
                  my $contcomp ;
                  $contcomp->{f} = _cluster_interfaces_compare_contact_sets({
                     aln1 => $fam1_aln,
                     aln2 => $fam2_aln,
                     revfl => [$pb->{revfl}->[$j], $pb->{revfl}->[$k]],
                     sid1 => [$cosid1->[0], $cosid1->[1]],
                     sid2 => [$cosid2->[0], $cosid2->[1]],
                     cont1 => $intconts->{$j}->{contacts},
                     cont2 => $intconts->{$k}->{contacts}
                  }) ;
                  
                  $contcomp->{b} = $contcomp->{f} ;
                  if ($fam1 eq $fam2) {
                     my $newrevfl2 = 1 ;
                     if ($pb->{revfl}->[$k]) { $newrevfl2 = 0 ; }
                     $contcomp->{b} = _cluster_interfaces_compare_contact_sets({
                        aln1 => $fam1_aln,
                        aln2 => $fam2_aln,
                        revfl => [$pb->{revfl}->[$j], $newrevfl2],
                        sid1 => [$cosid1->[0], $cosid2->[1]],
                        sid2 => [$cosid2->[0], $cosid1->[1]],
                        cont1 => $intconts->{$j}->{contacts},
                        cont2 => $intconts->{$k}->{contacts}
                     }) ;
                  }
   
                  if ($contcomp->{f} < $thresh->{aligned_contacts} &&
                      $contcomp->{b} < $thresh->{aligned_contacts}) { next; }
   
                  print STDERR "COMP $j v $k: f ($contcomp->{f})\n"
                     if $run_options->{DEBUGMERGER};
                  print STDERR "              b ($contcomp->{b})\n"
                     if $run_options->{DEBUGMERGER};
   
   
                  if ($contcomp->{b} > $contcomp->{f}) {
                     push @mergelist, [$j, $k, 1] ;
                  } else {
                     push @mergelist, [$j, $k, 0] ;
                  }
   
               }
            }
            #do mergelist
   
   
            my $currloc ;
            foreach my $cllabel (keys %{$iclust}) {
               $currloc->{$cllabel} = $cllabel ;  }
   
            foreach my $j ( 0 .. $#mergelist) {
               print STDERR "round $clusround MERGE #$j: $mergelist[$j]->[0] ".
                  "to $mergelist[$j]->[1]\n" if $run_options->{DEBUGMERGER};
   
               my $i1 = $mergelist[$j]->[0] ;
               my $i2 = $mergelist[$j]->[1] ;
               my $order = $mergelist[$j]->[2] ;
   
               my $curr1 = $i2clust->{$i1} ;
               my $curr2 = $i2clust->{$i2} ;
               if ($curr1 == $curr2) {next;}
   
               my $curr1_ord = $i2clust_order->{$i1} ;
               my $curr2_ord = $i2clust_order->{$i2} ;
   
               print STDERR "  currents: $curr1 and $curr2\n" if 
                  $run_options->{DEBUGMERGER};
   
               my $mergeorder = 0 ;
               if (($order == 1 && $curr1_ord == $curr2_ord) || 
                   ($order == 0 && $curr1_ord != $curr2_ord)) {
                  $mergeorder = 1 ; }
   
               foreach my $oldi2 ( @{$iclust->{$curr2}}) {
                  if ($mergeorder) {
                     if ($i2clust_order->{$oldi2} == 1) {
                        $i2clust_order->{$oldi2} = 0 ;
                     } else {
                        $i2clust_order->{$oldi2} = 1 ;
                     }
                  }
                  $i2clust->{$oldi2} = $curr1 ;
               }
   
               push @{$iclust->{$curr1}}, @{$iclust->{$curr2}} ;
               delete $iclust->{$curr2} ;
   
   
               if ($run_options->{DEBUGMERGER}) {
                  print STDERR " merging cluster $curr2 with $curr1\n" ;
                  print STDERR "set $curr1 now ".join(",", @{$iclust->{$curr1}})."\n";
                  print STDERR "   Deleted set $curr2\n" ;
   
                  print STDERR "\tactive clusters: ".
                     join(",", sort {$a <=> $b} keys %{$iclust})."\n" ;
               }
            }
            $clusround++ ;
         }
   
         #printout cluster memberships for fam12 interfaces

         my $cluster_no = 1 ;
         foreach my $label (sort keys %{$iclust}) {
            my $interface_class = $fam1."_".$fam2."_".$cluster_no ;
            my $scopclass_pair = $fam1."_".$fam2 ;
            my $member_no = 1 ;
            foreach my $int (@{$iclust->{$label}}) {
               my @outvals = ( $pb->{bdp_id}->[$int],
                               $pb->{sid1}->[$int],
                               $pb->{sid2}->[$int],
                               $scopclass_pair,
                               'fam',
                               $interface_class,
                               $cluster_no,
                               $member_no ) ;
               print join("\t", @outvals)."\n" ;
               $member_no++ ;
            }
            $cluster_no++ ;
         }

#ORIGINAL OUTPUT
#         my $nn = 1 ;
#         foreach my $label (keys %{$iclust}) {
#            my $nm = 1 ;
#            foreach my $int (@{$iclust->{$label}}) {
#               my @outvals = ($fam1, $fam2, $nn,
#                              $pb->{sid1}->[$int],
#                              $pb->{sid2}->[$int], $nm,
#                              $i2clust_order->{$int});
#               print join("\t", @outvals)."\n" ;
#               $nm++ ;
#            }
#            $nn++ ;
#         }
      }
   
   }

   my $import_status ;
   if (exists $in->{import_fl} && $in->{import_fl} == 1) {
# upload bdp_residues_tables file (specified in specs) into pibase
      $import_status = pibase::mysqlimport({
         pibase_specs => $pibase_specs,
         fn => $pibase_specs->{buildfiles}->{scop_interface_clusters}}) ;
   }

}


=head2 _cluster_interfaces_compare_residue_sets()

   Title:       _cluster_interfaces_compare_residue_sets()
   Function:    Venn comparison of interface residue sets 
   Args:        ->{aln} = alignment data
                ->{sid1} = domain 1 identifier
                ->{sid2} = domain 2 identifier
                ->{res1}->{ resno1 => 1, resno2 => 1} - residues in set 1
                ->{res2}->{ resno1 => 1, resno2 => 1} - residues in set 2
   Returns:     union / (union + diff) of the two sets

=cut

sub _cluster_interfaces_compare_residue_sets {

   my $in = shift;;
   my $aln = $in->{aln} ;
   my $sid1 = $in->{sid1} ;
   my $sid2 = $in->{sid2} ;
   my $res1 = $in->{res1} ;
   my $res2 = $in->{res2} ;

   my $intpos ;
   foreach my $t_res1 (keys %{$res1}) {
      my $pos = $aln->{resno2pos}->{$sid1}->{$t_res1} ;
      if (!defined $pos) {
         print STDERR "ERROR: compare_residue_sets(): $t_res1 from $sid1 undefined pos; aborting $sid1 -- $sid2 comparison\n" ; return 0;
      } else {
         $intpos->{$pos}++ ;
      }}


   foreach my $t_res2 (keys %{$res2}) {
      my $pos = $aln->{resno2pos}->{$sid2}->{$t_res2} ;
      if (!defined $pos) {
         print STDERR "ERROR: compare_residue_sets(): $t_res2 from $sid2 undefined pos; aborting $sid1 -- $sid2 comparison\n" ; return 0;
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


=head2 _cluster_interfaces_load_interface_contacts()

   Title:       _cluster_interfaces_load_interface_contacts()
   Function:    Venn comparison of interface residue sets 
   Args:        ->{sid1} = domain 1 identifier
                ->{sid2} = domain 2 identifier
   Returns:     ->{intres}->{sid1}->{resno1."\n".chain1}= domain 1 res in interface
                ->{intres}->{sid2}->{resno2."\n".chain2}= domain 2 res in interface
                ->{contacts}->{resno1."\n".chain1."\n".resno2."\n".chain2} -
                   inter-domain contact

=cut

sub _cluster_interfaces_load_interface_contacts {

   my $in = shift ;
   my $sid1 = $in->{sid1} ;
   my $sid2 = $in->{sid2} ;
   my $fn = $in->{fn} ;

   my ( $subset_id_1, $subset_id_2,
        $chain_id_1, $resno_1, $chain_id_2, $resno_2 ) =
       pibase::rawselect_metatod($fn, "SELECT subset_id_1, subset_id_2, chain_id_1, resno_1, chain_id_2, resno_2 FROM $fn") ;

   my $contacts ;
   my $intres ;
   foreach my $j ( 0 .. $#{$subset_id_1}) {
      if ($subset_id_1->[$j] ne $sid1) {next;}
      if ($subset_id_2->[$j] ne $sid2) {next;}

#      print STDERR "$subset_id_1->[$j] -- $subset_id_2->[$j]: $resno_1->[$j] $chain_id_1->[$j] to $resno_2->[$j] $chain_id_2->[$j]\n" ;

      $intres->{$sid1}->{$resno_1->[$j]."\n".$chain_id_1->[$j]}++ ;
      $intres->{$sid2}->{$resno_2->[$j]."\n".$chain_id_2->[$j]}++ ;
      $contacts->{$resno_1->[$j]."\n".$chain_id_1->[$j]."\n".$resno_2->[$j]."\n".$chain_id_2->[$j]}++ ;
   }

   return {
      intres => $intres,
      contacts => $contacts
   } ;

}


=head2 _cluster_interfaces_compare_contact_sets()

   Title:       _cluster_interfaces_compare_contact_ses()
   Function:    Venn comparison of interface contact sets 
   Args:        ->{aln1} = alignment data for domain type 1
                ->{aln2} = alignment data for domain type 2
                ->{sid1} = domain 1 identifier
                ->{sid2} = domain 2 identifier
                ->{cont1} = contacts for interface 1
                ->{cont2} = contacts for interface 2
                ->{revfl}->[i] = reversal flag;
                   if 1 switch sid1/2 in ith interface

   Returns:     union / (union + diff) of contacts

=cut

sub _cluster_interfaces_compare_contact_sets {

   my $in = shift;
   my $aln1 = $in->{aln1} ;
   my $aln2 = $in->{aln2} ;
   my $sid1 = $in->{sid1} ;
   my $sid2 = $in->{sid2} ;
   my $cont1 = $in->{cont1} ;
   my $cont2 = $in->{cont2} ;
   my $revfl = $in->{revfl} ;


   my $contacts;

   foreach my $res12 (keys %{$cont1}) {
      my ($pos1, $pos2) ;
      my @t = split(/\n/, $res12) ;
      my $res1 = $t[0]."\n".$t[1] ; my $res2 = $t[2]."\n".$t[3] ;

      if ($revfl->[0]) {
         $pos1 = $aln1->{resno2pos}->{$sid1->[0]}->{$res2} ;
         $pos2 = $aln2->{resno2pos}->{$sid2->[0]}->{$res1} ;
      } else {
         $pos1 = $aln1->{resno2pos}->{$sid1->[0]}->{$res1} ;
         $pos2 = $aln2->{resno2pos}->{$sid2->[0]}->{$res2} ;
      }
      
      if (!defined $pos1) {
         print STDERR "ERROR: compare_contacts(): undefined pos1 in $sid1->[0] (--$sid2->[0]) res $res1 (aborting $sid1->[0] -- $sid2->[0] -- $sid1->[1] -- $sid2->[1] comparison)\n" ;
         return 0 ;
      }

      if (!defined $pos2) {
         print STDERR "ERROR: compare_contacts(): undefined pos2 in $sid2->[0] (-- 1: $sid1->[0]) res $res2 (aborting $sid1->[0] -- $sid2->[0] -- $sid1->[1] -- $sid2->[1] comparison)\n" ; return 0 ;
      }

      $contacts->{$pos1."\n".$pos2}++;
   }


   foreach my $res12 (keys %{$cont2}) {
      my ($pos1, $pos2) ;
      my @t = split(/\n/, $res12) ;
      my $res1 = $t[0]."\n".$t[1] ; my $res2 = $t[2]."\n".$t[3] ;

      if ($revfl->[1]) {
         $pos1 = $aln1->{resno2pos}->{$sid1->[1]}->{$res2} ;
         $pos2 = $aln2->{resno2pos}->{$sid2->[1]}->{$res1} ;
      } else {
         $pos1 = $aln1->{resno2pos}->{$sid1->[1]}->{$res1} ;
         $pos2 = $aln2->{resno2pos}->{$sid2->[1]}->{$res2} ;
      }

      if (!defined $pos1) {
         print STDERR "compare_contacts(): undefined pos1 in $sid1->[1] (-- $sid2->[1]) res $res1 (aborting $sid1->[0] -- $sid2->[0] -- $sid1->[1] -- $sid2->[1] comparison)\n" ; return 0 ;
      }
      
      if (!defined $pos2) {
         print STDERR "compare_contacts(): undefined pos2 in $sid2->[1] (-- 1: $sid1->[1]) res $res2 (aborting $sid1->[0] -- $sid2->[0] -- $sid1->[1] -- $sid2->[1] comparison)\n" ; return 0;
      }

      $contacts->{$pos1."\n".$pos2}++;
   }


   my $union = 0;
   my $diff = 0 ;
   foreach my $contact (keys %{$contacts}) {
      if ($contacts->{$contact} == 2) {
         $union++;
      } else {
         $diff++ ;
      }
   }

   if (($union + $diff) > 0) {
      return ($union / ($union + $diff)) ;
   } else {
      return 0 ;
   }
}


=head2 _cluster_interface_pibase_preload()

   Title:       _cluster_interface_pibase_preload()
   Function:    Preloads PIBASE data necessary for SCOP interface clustering
   Args:        ->{astral} = astral data
                ->{aln2} = alignment data for domain type 2
                ->{sid1} = domain 1 identifier
                ->{sid2} = domain 2 identifier
                ->{cont1} = contacts for interface 1
                ->{cont2} = contacts for interface 2
                ->{revfl}->[i] = reversal flag; if 1 switch sid1/2 in ith interface

   Returns:     Huge hash of PIBASE data, see code for field explanations
      bdp2contactsfn => $bdp2contactsfn,
      bdp2pdb => $bdp2pdb,
      pdb2bdp => $pdb2bdp,
      bdp_id => $bdp_id,
      sid1 => $sid1,
      sid2 => $sid2,
      class1 => $class1,
      class2 => $class2,
      osid1 => $osid1,
      osid2 => $osid2,
      revfl => $revfl,
      fampairs => $fampairs,
      fampairs_single => $fampairs_single,
      fampairs_single_osid => $fampairs_single_osid,
      pdbchains => $pdbchains,
      chain_2_pdbchain => $chain_2_pdbchain,
      chain_2_start => $chain_2_start,
      chain_2_end => $chain_2_end,
      chain_2_startser => $chain_2_startser,
      chain_2_endser => $chain_2_endser

=cut

sub _cluster_interfaces_pibase_preload {

   my $in = shift ;
   my $astral = $in->{astral} ;

   my $bdp2contactsfn ;
   {
#      "SELECT bdp_id, source_file FROM interface_contacts_tables"
      my ($t_bdp, $t_fn) = pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM interface_contacts_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         $bdp2contactsfn->{$t_bdp->[$j]} = $t_fn->[$j] ;
      }
   }

   my $pdb2bdp ; my $bdp2pdb;
   my $rawbdp_list ;
   {
#   pdb2bdp: "SELECT pdb_id, bdp_id FROM bdp_files WHERE raw_pdb =1 " ;
#   bdp2pdb: "SELECT bdp_id, pdb_id FROM bdp_files") ;
      my ($t_pdb, $t_bdp, $t_raw) = pibase::rawselect_tod(
         "SELECT pdb_id, bdp_id, raw_pdb FROM bdp_files") ;
      foreach my $j ( 0.. $#{$t_pdb}) {
         $bdp2pdb->{$t_bdp->[$j]} = $t_pdb->[$j] ;

         if ($t_raw->[$j] != 1) {next;}
         $pdb2bdp->{$t_pdb->[$j]} = $t_bdp->[$j] ;
         $rawbdp_list->{$t_bdp->[$j]}++ ;
      }
   }


   my ($bdp_id, $sid1, $sid2, $class1, $class2) ; 
   {
#      "SELECT a.bdp_id, subset_id_1, subset_id_2, class_1, class_2 ".
#      "FROM intersubset_contacts as a, bdp_files as b ".
#      "WHERE subset_id_1 LIKE \"\%SCOP\%\" AND a.bdp_id = b.bdp_id ".
#      "AND b.raw_pdb = 1" ;
      my ($t_bdp, $t_sid1, $t_sid2, $t_class1, $t_class2) = 
         pibase::rawselect_tod(
         "SELECT bdp_id, subset_id_1, subset_id_2, class_1, class_2 ".
         "FROM intersubset_contacts") ;
      foreach my $j ( 0 .. $#{$t_bdp}) {
         if (!exists $rawbdp_list->{$t_bdp->[$j]}) {next;}
         push @{$bdp_id}, $t_bdp->[$j] ;
         push @{$sid1}, $t_sid1->[$j] ;
         push @{$sid2}, $t_sid2->[$j] ;
         push @{$class1}, $t_class1->[$j] ;
         push @{$class2}, $t_class2->[$j] ;
      }
   }

   my ($fampairs, $fampairs_single,$fampairs_single_osid) ;
   my ($osid1, $osid2, $revfl) ; #revfl if fam2 <fam1 and flipped to get fam21 isntead of fam12, let it be known so later it can be dealt with

   foreach my $j ( 0 .. $#{$bdp_id}) {
      if ($class1->[$j] !~ /^[a-g]/ || $class2->[$j] !~ /^[a-g]/) {next;}

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
      push @{$fampairs->{$fam12}}, $j ;
      $fampairs_single->{$fam12}->{$class1->[$j]}->{$sid1->[$j]}++ ;
      $fampairs_single->{$fam12}->{$class2->[$j]}->{$sid2->[$j]}++ ;

      $fampairs_single_osid->{$fam12}->{$class1->[$j]}->{$osid1->[$j]}++ ;
      $fampairs_single_osid->{$fam12}->{$class2->[$j]}->{$osid2->[$j]}++ ;
   }
   print STDERR "X\n" ;


   my $chain_2_pdbchain ;
   my ($chain_2_start, $chain_2_end) ;
   my ($chain_2_startser, $chain_2_endser) ;
#   my $chain_incr ;
   my $pdbchains ;
   {
#      my ($tbdp_id, $real_chain_id, $pdb_chain_id, $resno_start,
#         $resno_start_serial, $resno_end, $resno_end_serial) =
#      pibase::mysql_fetchcols($dbh,
#      "SELECT a.bdp_id, real_chain_id, pdb_chain_id, ".
#      "start_resno, start_resno_int, ".
#      "end_resno, end_resno_int FROM bdp_chains as a, bdp_files as b ".
#      "WHERE pdb_chain_id IS NOT NULL AND chain_type = \"p\" AND ".
#      "a.bdp_id = b.bdp_id AND raw_pdb = 1");

      my ($tbdp_id, $real_chain_id, $pdb_chain_id, $resno_start,
         $resno_start_serial, $resno_end, $resno_end_serial,
         $t_chain_type) =
      pibase::rawselect_tod(
         "SELECT bdp_id, real_chain_id, pdb_chain_id, ".
         "start_resno, start_resno_int, end_resno, end_resno_int, ".
         "chain_type FROM bdp_chains") ;
#         as a, bdp_files as b ".
#         "WHERE pdb_chain_id IS NOT NULL AND chain_type = \"p\" AND ".
#         "a.bdp_id = b.bdp_id AND raw_pdb = 1");

#      print STDERR ($#{$tbdp_id} + 1)." chains loaded entries\n " ;

      foreach my $j ( 0 .. $#{$tbdp_id}) {
         if ($pdb_chain_id->[$j] eq 'NULL' ||
             $t_chain_type->[$j] ne 'p' ||
             !exists $rawbdp_list->{$tbdp_id->[$j]}) { next;}

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

#      foreach my $j ( 0 .. $#{$tbdp_id}) {
#         print STDERR "NOW ON $tbdp_id->[$j] chain $real_chain_id->[$j]\n" ;
#         if ($real_chain_id->[$j] ne $pdb_chain_id->[$j]) {
#
#            my $rawbdp = $pdb2bdp->{$bdp2pdb->{$tbdp_id->[$j]}} ;
#
#            $chain_incr->{$tbdp_id->[$j]} = $resno_start_serial->[$j] - 
#               $chain_2_startser->{$rawbdp}->{$pdb_chain_id->[$j]} ;
#
#            if ($chain_incr->{$tbdp_id->[$j]} != ($resno_end_serial->[$j] - 
#               $chain_2_endser->{$rawbdp}->{$pdb_chain_id->[$j]})) {
#               print STDERR " OH SHIT (bdp $tbdp_id->[$j] $bdp2pdb->{$tbdp_id->[$j]}): chain $real_chain_id->[$j] increment not so clear - start and ends vary\n" ;
#            }
#         }
#      }
   }

   my $pb = {
      bdp2contactsfn => $bdp2contactsfn,
      bdp2pdb => $bdp2pdb,
      pdb2bdp => $pdb2bdp,
      bdp_id => $bdp_id,
      sid1 => $sid1,
      sid2 => $sid2,
      class1 => $class1,
      class2 => $class2,
      osid1 => $osid1,
      osid2 => $osid2,
      revfl => $revfl,
      fampairs => $fampairs,
      fampairs_single => $fampairs_single,
      fampairs_single_osid => $fampairs_single_osid,
      pdbchains => $pdbchains,
      chain_2_pdbchain => $chain_2_pdbchain,
      chain_2_start => $chain_2_start,
      chain_2_end => $chain_2_end,
      chain_2_startser => $chain_2_startser,
      chain_2_endser => $chain_2_endser
   } ;

#   my @loaded_fampairs = sort keys %{$fampairs} ;
#   foreach my $j ( 0 .. $#loaded_fampairs) {
#      print "$j $loaded_fampairs[$j]\n";
#   }

   return $pb ;

}


sub OLDDBI__cluster_interfaces_pibase_preload {

   my $in = shift ;
   my $astral = $in->{astral} ;

   my ($dbh) = pibase::connect_pibase() ;

   my $bdp2contactsfn = pibase::mysql_hashload($dbh,
      "SELECT bdp_id, source_file FROM interface_contacts_tables") ;

   my ($bdp_id, $sid1, $sid2, $class1, $class2) = 
      pibase::mysql_fetchcols($dbh,
      "SELECT a.bdp_id, subset_id_1, subset_id_2, class_1, class_2 ".
      "FROM intersubset_contacts as a, bdp_files as b ".
      "WHERE subset_id_1 LIKE \"\%SCOP\%\" AND a.bdp_id = b.bdp_id ".
      "AND b.raw_pdb = 1") ;

   my ($fampairs, $fampairs_single,$fampairs_single_osid) ;
   my ($osid1, $osid2, $revfl) ; #revfl if fam2 <fam1 and flipped to get fam21 isntead of fam12, let it be known so later it can be dealt with

   foreach my $j ( 0 .. $#{$bdp_id}) {
      if ($class1->[$j] !~ /^[a-g]/ || $class2->[$j] !~ /^[a-g]/) {next;}

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
      push @{$fampairs->{$fam12}}, $j ;
      $fampairs_single->{$fam12}->{$class1->[$j]}->{$sid1->[$j]}++ ;
      $fampairs_single->{$fam12}->{$class2->[$j]}->{$sid2->[$j]}++ ;

      $fampairs_single_osid->{$fam12}->{$class1->[$j]}->{$osid1->[$j]}++ ;
      $fampairs_single_osid->{$fam12}->{$class2->[$j]}->{$osid2->[$j]}++ ;
   }
   print STDERR "X\n" ;

   my $pdb2bdp = pibase::mysql_hashload($dbh,
      "SELECT pdb_id, bdp_id FROM bdp_files WHERE raw_pdb =1 ") ;

   my $bdp2pdb = pibase::mysql_hashload($dbh,
      "SELECT bdp_id, pdb_id FROM bdp_files") ;

   my $chain_2_pdbchain ;
   my ($chain_2_start, $chain_2_end) ;
   my ($chain_2_startser, $chain_2_endser) ;
#   my $chain_incr ;
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

#      foreach my $j ( 0 .. $#{$tbdp_id}) {
#         print STDERR "NOW ON $tbdp_id->[$j] chain $real_chain_id->[$j]\n" ;
#         if ($real_chain_id->[$j] ne $pdb_chain_id->[$j]) {
#
#            my $rawbdp = $pdb2bdp->{$bdp2pdb->{$tbdp_id->[$j]}} ;
#
#            $chain_incr->{$tbdp_id->[$j]} = $resno_start_serial->[$j] - 
#               $chain_2_startser->{$rawbdp}->{$pdb_chain_id->[$j]} ;
#
#            if ($chain_incr->{$tbdp_id->[$j]} != ($resno_end_serial->[$j] - 
#               $chain_2_endser->{$rawbdp}->{$pdb_chain_id->[$j]})) {
#               print STDERR " OH SHIT (bdp $tbdp_id->[$j] $bdp2pdb->{$tbdp_id->[$j]}): chain $real_chain_id->[$j] increment not so clear - start and ends vary\n" ;
#            }
#         }
#      }
   }


   my $pb = {
      bdp2contactsfn => $bdp2contactsfn,
      bdp2pdb => $bdp2pdb,
      pdb2bdp => $pdb2bdp,
      bdp_id => $bdp_id,
      sid1 => $sid1,
      sid2 => $sid2,
      class1 => $class1,
      class2 => $class2,
      osid1 => $osid1,
      osid2 => $osid2,
      revfl => $revfl,
      fampairs => $fampairs,
      fampairs_single => $fampairs_single,
      fampairs_single_osid => $fampairs_single_osid,
      pdbchains => $pdbchains,
      chain_2_pdbchain => $chain_2_pdbchain,
      chain_2_start => $chain_2_start,
      chain_2_end => $chain_2_end,
      chain_2_startser => $chain_2_startser,
      chain_2_endser => $chain_2_endser
   } ;

}



1;
