=head1 NAME

pibase::build - perl module to build the pibase database

=head1 DESCRIPTION

Perl module that executes the PIBASE build protocol.

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


package pibase::build ;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/build_pibase/ ;

use Cwd ;

use pibase qw/get_specs/;
use pibase::SCOP ;
use pibase::CATH ;
use pibase::PDB ;
#use pibase::PQS ;
use pibase::PISA ;
use pibase::calc::interfaces ;
#use pibase::data::external::PQS ;
use pibase::data::external::PISA ;
use pibase::data::calc ;
use pibase::data::access ;
use pibase::web ;

use File::Temp qw/tempfile tempdir/ ;
use File::Path qw/mkpath/ ;


=head2 build_pibase()

   Title:       build_pibase()
   Function:    Builds the PIBASE database.
   Args:        $_->{steps}->{pibasebuild_step} = 1 - optional, if specified,
                 only performs this particular step of the build process
   Returns:     Nothing

=cut

sub build_pibase {

#Not initially setting up directory structure- do it online prn
#Step -1. build mysql table structures.

   my $in = shift ;

# Step 0. fetch external data
   my $build_status ;
   my $fetch_external_data = fetch_external_data() ;
   my $pibase_specs = $fetch_external_data->{pibase_specs} ;
   $build_status->{fetch_external_data} = $fetch_external_data->{status} ;
   if (exists $build_status->{fetch_external_data}->{error_fl}) {
      die $build_status->{fetch_external_data}->{error_fl} ; }
   
   if (exists $in->{steps}->{download_asteroids} ||
       !-d $pibase_specs->{asteroids}->{sf_aln}) {
      download_asteroids({pibase_specs => $pibase_specs}) ;
   }

#commented out100906_1732  - PQS deprecated
## Step -1. get PQS mirror going - only if explicitly specified
#   if (exists $in->{steps}->{mirror_pqs}) {
#      pibase::data::external::PQS::list_pqs_file_urls({
#         pibase_specs => $pibase_specs}) ;
#      die ;
#   }
#
#   if (exists $in->{steps}->{check_pqs_files}) {
#      pibase::data::external::PQS::check_pqs_files({
#         pibase_specs => $pibase_specs}) ;
#      die ;
#   }

# Step 0. setup PISA mirror
   if (exists $in->{steps}->{mirror_pisa}) {
      pibase::data::external::PISA::get_pisa_files({
         pibase_specs => $pibase_specs}) ;
      die "mirror_pisa done";
   }

   if (exists $in->{steps}->{check_pisa_files}) {
      pibase::data::external::PISA::check_pisa_files({
         pibase_specs => $pibase_specs}) ;
      die "check_pisa_files done";
   }

# Step 1. process and import external data
   my $num_specified_steps = keys %{$in->{steps}} ;
   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{process_external_data} &&
      $in->{steps}->{process_external_data} == 1)) {
      $build_status->{process_external_data} = process_external_data({
         pibase_specs => $pibase_specs }) ;
   }

# Step 2. Extract NMR model 1 structures

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{extract_nmr_model1} &&
      $in->{steps}->{extract_nmr_model1} == 1)) {
      $build_status->{extract_nmr_model1} = pibase::PDB::nmr_model1_extractor({
         pibase_specs => $pibase_specs }) ;
   }

# Step 3. Make bdp_files table - internal structure list

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{maketable_bdp_files} &&
      $in->{steps}->{maketable_bdp_files} == 1)) {
      $build_status->{maketable_bdp_files} =
         pibase::data::access::maketable_bdp_files({
            import_fl => 1,
            pibase_specs => $pibase_specs }) ;
   }

# Step 4. data::calc::pdb_chaininfo_caller_calc

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{call_chain_info} &&
      $in->{steps}->{call_chain_info} == 1)) {
      $build_status->{call_chain_info} =
         pibase::data::calc::call_chain_info({
            numtasks_cluster => 20,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs }) ;
   }

# Step 5. data::calc::bdp_resinfo_caller_calc

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{call_residue_info} &&
      $in->{steps}->{call_residue_info} == 1)) {
      $build_status->{call_residue_info} =
         pibase::data::calc::call_residue_info({
            numtasks_cluster => 20,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs }) ;
   }


# Step 6. determine pisa-pdb chain equivalencies

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
      exists $in->{steps}->{calc_bdp_chains_equiv} &&
      $in->{steps}->{calc_bdp_chains_equiv} == 1)) {
      $build_status->{calc_bdp_chains_equiv} =
         pibase::data::calc::calc_bdp_chains_equiv({
            numtasks_cluster => 20,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs }) ;
   }

# Step 7. Determine subset definitions for PISA entries
#         and generate pdb chain-based "subsets"

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
       exists $in->{steps}->{bdp_subset_translator} &&
       $in->{steps}->{bdp_subset_translator} == 1)) {

# Step 7a. Dump tables in tod directory

# check that there is an entry for subset_source = 'pdb_chain' ;

      pibase::data::access::copy_table_to_disk({
         tables => ['bdp_chains', 'bdp_files',
            'subsets', 'subsets_details', 'subsets_class', 'subsets_source',
            'bdp_residues_tables',
            'pdb_entries', 'pdb_entry_type']
      }) ;

      $build_status->{bdp_subset_translator} =
         pibase::data::calc::bdp_subset_translator({
            numtasks_cluster => 20,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['subsets', 'subsets_details', 'subsets_source',
                    'subsets_residues_tables' ]
      }) ;
   }

# Step 8. Detect interfaces (on the cluster)
#         and generate pdb chain-based "subsets"

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         exists $in->{steps}->{interface_detect_calc} &&
         $in->{steps}->{interface_detect_calc} == 1)) {

      pibase::data::access::copy_table_to_disk({
         tables => ['subsets', 'subsets_details', 'subsets_source',
                    'subsets_residues_tables' ]
      }) ;

      $build_status->{interface_detect} =
         pibase::calc::interfaces::interface_detect_calc({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['intersubset_contacts',
                    'patch_residues_tables',
                    'interface_contacts_tables',
                    'interface_contacts_special_tables' ]
      }) ;
   }


# Step 9. Calculate BDP entry properties (on the cluster)

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         ((exists $in->{steps}->{calc_bdp_properties} &&
           $in->{steps}->{calc_bdp_properties} == 1) ||
          (exists $in->{steps}->{calc_bdp_secstrx} &&
           $in->{steps}->{calc_bdp_secstrx} == 1)))) {

      $build_status->{calc_bdp_secstrx} =
         pibase::data::calc::calc_bdp_secstrx({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['bdp_secstrx_tables']
      }) ;

   }

# Step 10a. Create subset pdb files (on the cluster)

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         exists $in->{steps}->{create_subset_pdbs} &&
         $in->{steps}->{create_subset_pdbs} == 1)) {

      $build_status->{create_subset_pdbs} =
         pibase::data::access::create_subset_pdbs({
            numtasks_cluster => 20,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['subsets_files']
      }) ;

   }

# Step 10b. Calculate subset properties (on the cluster)

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         ((exists $in->{steps}->{calc_subset_properties} &&
           $in->{steps}->{calc_subset_properties} == 1) ||
          (exists $in->{steps}->{calc_subsets_sasa} &&
           $in->{steps}->{calc_subsets_sasa} == 1)))) {

      $build_status->{calc_subsets_sasa} =
         pibase::data::calc::calc_subsets_sasa({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['subsets_sasa']
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         ((exists $in->{steps}->{calc_subset_properties} &&
           $in->{steps}->{calc_subset_properties} == 1) ||
          (exists $in->{steps}->{calc_subsets_sequence} &&
           $in->{steps}->{calc_subsets_sequence} == 1)))) {

      $build_status->{calc_subsets_sequence} =
         pibase::data::calc::calc_subsets_sequence({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['subsets_sequence']
      }) ;

   }

# Step 11. Calculate interface properties (on the cluster)
   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         ((exists $in->{steps}->{calc_interface_properties} &&
           $in->{steps}->{calc_interface_properties} == 1) ||
          (exists $in->{steps}->{calc_interface_dsasa} &&
           $in->{steps}->{calc_interface_dsasa} == 1)))) {

      $build_status->{calc_interface_dsasa} =
         pibase::data::calc::calc_interface_dsasa({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['interface_sasa']
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         ((exists $in->{steps}->{calc_interface_properties} &&
           $in->{steps}->{calc_interface_properties} == 1) ||
          (exists $in->{steps}->{calc_interface_secstrx} &&
           $in->{steps}->{calc_interface_secstrx} == 1)))) {

      $build_status->{calc_interface_secstrx} =
         pibase::data::calc::calc_interface_secstrx({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['interface_secstrx_tables']
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         ((exists $in->{steps}->{calc_interface_properties} &&
           $in->{steps}->{calc_interface_properties} == 1) ||
          (exists $in->{steps}->{calc_interface_bs_secstrx_basic} &&
           $in->{steps}->{calc_interface_bs_secstrx_basic} == 1)))) {

      $build_status->{calc_interface_bs_secstrx_basic} =
         pibase::data::calc::calc_interface_bs_secstrx_basic({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['interface_secstrx_basic_contacts_tables',
                    'bindingsite_secstrx_basic_contacts_tables',
                    'bindingsite_contacts_tables']
      }) ;
   }


   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         ((exists $in->{steps}->{calc_interface_properties} &&
           $in->{steps}->{calc_interface_properties} == 1) ||
          (exists $in->{steps}->{calc_interface_secstrx_contacts} &&
           $in->{steps}->{calc_interface_secstrx_contacts} == 1)))) {

      $build_status->{calc_interface_secstrx_contacts} =
         pibase::data::calc::calc_interface_secstrx_contacts({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['interface_secstrx_contacts_tables']
      }) ;
   }


   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         ((exists $in->{steps}->{calc_interface_properties} &&
           $in->{steps}->{calc_interface_properties} == 1) ||
          (exists $in->{steps}->{calc_interface_secstrx_profile} &&
           $in->{steps}->{calc_interface_secstrx_profile} == 1)))) {

      $build_status->{calc_interface_secstrx_profile} =
         pibase::data::calc::calc_interface_secstrx_profile({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['interface_secstrx_profile']
      }) ;
   }


# calc_interface_resvector
   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         ((exists $in->{steps}->{calc_interface_properties} &&
           $in->{steps}->{calc_interface_properties} == 1) ||
          (exists $in->{steps}->{calc_interface_resvector} &&
           $in->{steps}->{calc_interface_resvector} == 1)))) {

      $build_status->{calc_interface_secstrx_profile} =
         pibase::data::calc::calc_interface_resvector({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['interface_resvector']
      }) ;
   }

# calc_interface_sse_topology
   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         ((exists $in->{steps}->{calc_interface_properties} &&
           $in->{steps}->{calc_interface_properties} == 1) ||
          (exists $in->{steps}->{calc_interface_sse_topology} &&
           $in->{steps}->{calc_interface_sse_topology} == 1)))) {

      $build_status->{calc_interface_sse_topology} =
         pibase::data::calc::calc_interface_sse_topology({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['interface_sse_topology',
                    'bindingsite_sse_topology']
      }) ;
   }

# calc_interface_size
   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         ((exists $in->{steps}->{calc_interface_properties} &&
           $in->{steps}->{calc_interface_properties} == 1) ||
          (exists $in->{steps}->{calc_interface_size} &&
           $in->{steps}->{calc_interface_size} == 1)))) {

      $build_status->{calc_interface_size} =
         pibase::data::calc::calc_interface_size({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['interface_size']
      }) ;
   }


# calc_bdp_interaction_topology
   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         ((exists $in->{steps}->{calc_interface_properties} &&
           $in->{steps}->{calc_interface_properties} == 1) ||
          (exists $in->{steps}->{calc_bdp_interaction_topology} &&
           $in->{steps}->{calc_bdp_interaction_topology} == 1)))) {

      $build_status->{calc_bdp_interaction_topology} =
         pibase::data::calc::calc_bdp_interaction_topology({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['bdp_interaction_topology']
      }) ;
   }

   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         ((exists $in->{steps}->{calc_interface_properties} &&
           $in->{steps}->{calc_interface_properties} == 1) ||
          (exists $in->{steps}->{calc_bdp_interaction_topology_graph} &&
           $in->{steps}->{calc_bdp_interaction_topology_graph} == 1)))) {

      $build_status->{calc_bdp_interaction_topology_graph} =
         pibase::data::calc::calc_bdp_interaction_topology_graph({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['bdp_interaction_topology_graph']
      }) ;
   }

# cluster scop interfaces
   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         exists $in->{steps}->{cluster_scop_interfaces} &&
         $in->{steps}->{cluster_scop_interfaces} == 1)) {

      $build_status->{cluster_scop_interfaces} =
         pibase::calc::interfaces::cluster_scop_interfaces({
            numtasks_cluster => 100,
            cluster_fl => 1,
            import_fl => 1,
            pibase_specs => $pibase_specs
         }) ;

      pibase::data::access::copy_table_to_disk({
         tables => ['scop_interface_clusters']
      }) ;
   }

# Step 12. Calculate complex properties (on the cluster)
   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         exists $in->{steps}->{calc_interface_properties} &&
         $in->{steps}->{calc_interface_properties} == 1)) {
   }

# Step 13. Convert topology graphs into png for web page use
   if ($num_specified_steps == 0 || (exists $in->{steps} &&
         ((exists $in->{steps}->{calc_interface_properties} &&
           $in->{steps}->{calc_interface_properties} == 1) ||
          (exists $in->{steps}->{pngconvert_bdp_interaction_topology_graph} &&
           $in->{steps}->{pngconvert_bdp_interaction_topology_graph} == 1)))) {

      $build_status->{calc_bdp_interaction_topology_graph} =
         pibase::web::pngconvert_bdp_interaction_topology_graph({
            pibase_specs => $pibase_specs
         }) ;
   }

}


=head2 fetch_external_data()

   Title:       fetch_external_data()
   Function:    Gets external data necessary for a PIBASE build.
   In:          ->{pibase_specs} = pibase_specs - optional, uses get_specs() if not
   Returns:     ->{pibase_specs} = pibase_specs structure
                ->{status}->{resource}->{status} = 1 if resource is good
                ->{status}->{resource}->{missing} = [fileurl of missing,...]

=cut

sub fetch_external_data {

   my $in = shift ;
   my $cursubroutine = (caller(0))[3];

   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) { $pibase_specs = pibase::get_specs() ;}
   else { $pibase_specs = $in->{pibase_specs};}

   my $status ;
   print "CHECKING EXTERNAL DATA ($cursubroutine) ".localtime()."\n"
      if (!exists $in->{quiet_fl});
   foreach my $resource (sort keys %{$pibase_specs->{external_data}}) {
      print "* $resource files:\n" if (!exists $in->{quiet_fl});
      my $resource_info = $pibase_specs->{external_data}->{$resource};

      if (! -s $resource_info->{local_dir}) {
         mkpath($resource_info->{local_dir}) ; }

      $status->{$resource}->{status} = 1 ;
      foreach my $file_url (keys %{$resource_info->{urls}}) {
         my $base = $file_url; $base =~ s/.*\/// ;
         my $local_copy = $resource_info->{local_dir}.'/'.$base ;
         print "   * $base: " if (!exists $in->{quiet_fl});
         if (! -s $local_copy) {
            my $tcom = "wget -P $resource_info->{local_dir} $file_url" ;

            print " downloading from $file_url\n"  if (!exists $in->{quiet_fl});
            print "# $tcom\n"  if (!exists $in->{quiet_fl});
            system($tcom) ;

            if (!-s $local_copy) {
               $status->{$resource}->{status} = 0 ;
               print " error\n"  if (!exists $in->{quiet_fl});
               print STDERR "ERROR: didn't get $file_url\n";
               push @{$status->{$resource}->{missing}}, $file_url ;
               $status->{error_fl} .= "ERROR $resource: didn't get $file_url\n";
            } else {
               if (exists $resource_info->{postcommands}->{$file_url}) {
                  my $orig_dir = getcwd() ;
                  chdir $resource_info->{local_dir} ;
                  system($resource_info->{postcommands}->{$file_url}) ;
                  chdir $orig_dir ;
               }
            }
         } else {
            print " x\n" if (!exists $in->{quiet_fl});
         }
      }
   }

   return {
      pibase_specs => $pibase_specs,
      status => $status
   } ;

}


=head2 process_external_data()

   Title:       process_external_data()
   Function:    Processes external data necessary for a PIBASE build.
   Returns:     ->{status}->{resource}->{status} = 1 if resource is good
                ->{status}->{resource}->{missing} = [fileurl of missing,...]

=cut

sub process_external_data {

   my $in = shift ;
   my $cursubroutine = (caller(0))[3];
   my $pibase_specs = $in->{pibase_specs} ;
   my $overwrite_fl = 0;
   if (exists $in->{overwrite_fl}) {$overwrite_fl = $in->{overwrite_fl};}
#   $in->{resources}->{pdb}++ ;

   my $status ;
   my ($temp_fh, $temp_fn) ;

   ($temp_fh->{cath}->{cnf_fixcontraction},
    $temp_fn->{cath}->{cnf_fixcontraction}) = tempfile(CLEANUP=>1);
   close($temp_fh->{cath}->{cnf_fixcontraction}) ;


   my $procsteps;
   $procsteps->{pisa} = [
      {
         name => 'pisa_clean_index',
         processor => \&pibase::PISA::pisa_clean_index,
         in_fn => $pibase_specs->{external_data}->{pisa}->{files}->{"index"},
         out_fn =>
            $pibase_specs->{external_data}->{pisa}->{processed_files}->{pisa_index},
         header_fl => 0,
      },
   ] ;

#   $procsteps->{pqs} = [
#      {
#         name => 'pqs_clean_list',
#         processor => \&pibase::PQS::pqs_clean_list,
#         in_fn => $pibase_specs->{external_data}->{pqs}->{files}->{list},
#         out_fn => $pibase_specs->{external_data}->{pqs}->{processed_files}->{pqs_list},
#         header_fl => 0,
#      },
#      {
#         name => 'pqs_clean_asalist',
#         processor => \&pibase::PQS::pqs_clean_asalist,
#         in_fn => $pibase_specs->{external_data}->{pqs}->{files}->{asalist},
#         out_fn => $pibase_specs->{external_data}->{pqs}->{processed_files}->{pqs_asalist},
#         header_fl => 0,
#      },
#      {
#         name => 'pqs_clean_biolist_contraction',
#         processor => \&pibase::PQS::pqs_clean_biolist_contraction,
#         in_fn => $pibase_specs->{external_data}->{pqs}->{files}->{biolist},
#         out_fn => $temp_fn->{pqs}->{biolist_fixcontraction},
#         header_fl => 0,
#      },
#      {
#         name => 'pqs_clean_biolist',
#         processor => \&pibase::PQS::pqs_clean_biolist,
#         in_fn => $temp_fn->{pqs}->{biolist_fixcontraction},
#         out_fn => $pibase_specs->{external_data}->{pqs}->{processed_files}->{pqs_biolist},
#         header_fl => 0,
#      },
#      {
#         name => 'pqs_clean_ranking',
#         processor => \&pibase::PQS::pqs_clean_ranking,
#         in_fn => $pibase_specs->{external_data}->{pqs}->{files}->{ranking},
#         out_fn => $pibase_specs->{external_data}->{pqs}->{processed_files}->{pqs_ranking},
#         header_fl => 0,
#      },
#   ] ;

#      {
#         name => "pibase_import_pqs_domains",
#         processor => \&pibase::CATH::pibase_import_pqs_domains,
#         in_fn => '',
#         out_fn => '',
#         header_fl => 0,
#      }

   $procsteps->{pdb} = [
      {
         name => 'pdb_clean_entries_idx',
         processor => \&pibase::PDB::pdb_clean_entries_idx,
         in_fn => $pibase_specs->{external_data}->{pdb}->{files}->{entries},
         out_fn => $pibase_specs->{external_data}->{pdb}->{processed_files}->{pdb_entries},
         header_fl => 0,
      },
      {
         name => 'pdb_clean_obsolete_dat',
         processor => \&pibase::PDB::pdb_clean_obsolete_dat,
         in_fn => $pibase_specs->{external_data}->{pdb}->{files}->{obsolete},
         out_fn => $pibase_specs->{external_data}->{pdb}->{processed_files}->{pdb_obsolete},
         header_fl => 0,
      },
      {
         name => 'pdb_clean_release_date',
         processor => \&pibase::PDB::pdb_clean_release_date,
         in_fn => $pibase_specs->{external_data}->{pdb}->{files}->{entries},
         out_fn => $pibase_specs->{external_data}->{pdb}->{processed_files}->{pdb_release},
         header_fl => 0,
      },
      {
         name => 'pdb_copy_entry_type',
         processor => \&pibase::PDB::pdb_copy_entry_type,
         in_fn => $pibase_specs->{external_data}->{pdb}->{files}->{pdb_entry_type},
         out_fn => $pibase_specs->{external_data}->{pdb}->{processed_files}->{pdb_entry_type},
         header_fl => 0,
      },
   ] ;

#      {
#         name => "pibase_import_pdb_domains",
#         processor => \&pibase::CATH::pibase_import_pdb_domains,
#         in_fn => '',
#         out_fn => '',
#         header_fl => 0,
#      }

#      {
#         name => 'cath_clean_clf CathChainList',
#         processor => \&pibase::CATH::cath_clean_clf,
#         in_fn => $pibase_specs->{external_data}->{cath}->{files}->{CathChainList},
#         out_fn => $pibase_specs->{external_data}->{cath}->{processed_files}->{cath_chain_list},
#         header_fl => 0,
#      },

   $procsteps->{cath} = [
      {
         name => 'cath_clean_cddf',
         processor => \&pibase::CATH::cath_clean_cddf,
         in_fn => $pibase_specs->{external_data}->{cath}->{files}->{CathDomainDescription},
         out_fn => $pibase_specs->{external_data}->{cath}->{processed_files}->{cath_domain_description},
         header_fl => 0,
      },
      {
         name => 'cath_clean_clf CathDomainList',
         processor => \&pibase::CATH::cath_clean_clf,
         in_fn => $pibase_specs->{external_data}->{cath}->{files}->{CathDomainList},
         out_fn => $pibase_specs->{external_data}->{cath}->{processed_files}->{cath_domain_list},
         header_fl => 0,
      },
      { 
         name => "cath_clean_cnf_contraction CathNames",
         processor => \&pibase::CATH::cath_clean_cnf_contraction,
         in_fn => $pibase_specs->{external_data}->{cath}->{files}->{CathNames},
         out_fn => $temp_fn->{cath}->{cnf_fixcontraction},
         header_fl => 0,
      },
      { 
         name => "cath_clean_cnf CathNames",
         processor => \&pibase::CATH::cath_clean_cnf,
         in_fn => $temp_fn->{cath}->{cnf_fixcontraction},
         out_fn => $pibase_specs->{external_data}->{cath}->{processed_files}->{cath_names},
         header_fl => 0,
      },
      {
         name => "pibase_import_cath_domains",
         processor => \&pibase::CATH::pibase_import_cath_domains,
         in_fn => '',
         out_fn => '',
         header_fl => 0,
      }
   ] ;


   $procsteps->{scop} = [
      {
         name => 'scop_hie',
         processor => \&pibase::SCOP::scop_clean_hie,
         in_fn => $pibase_specs->{external_data}->{scop}->{files}->{hie},
         out_fn =>
         $pibase_specs->{external_data}->{scop}->{processed_files}->{scop_hie},
         header_fl => 0,
      },
      {
         name => 'scop_cla',
         processor => \&pibase::SCOP::scop_clean_cla,
         in_fn => $pibase_specs->{external_data}->{scop}->{files}->{cla},
         out_fn =>
          $pibase_specs->{external_data}->{scop}->{processed_files}->{scop_cla},
         header_fl => 0,
      },
      {
         name => 'scop_des',
         processor => \&pibase::SCOP::scop_clean_des,
         in_fn => $pibase_specs->{external_data}->{scop}->{files}->{des},
         out_fn =>
        $pibase_specs->{external_data}->{scop}->{processed_files}->{scop_des},
         header_fl => 0,
      },
      {
         name => "pibase_import_scop_domains",
         processor => \&pibase::SCOP::pibase_import_scop_domains,
         in_fn => '',
         out_fn => '',
         header_fl => 0,
      }
   ] ;

   print "#PROCESSING EXTERNAL DATA ".$cursubroutine." ".localtime()."\n"
      if (!exists $in->{quiet_fl});

   my @resources = sort keys %{$procsteps} ;
   if (exists $in->{resources}) {
      @resources = sort keys %{$in->{resources}} ; }

   foreach my $resource (@resources) {
      print "# NOW ON resource $resource\n" if (!exists $in->{quiet_fl});
      foreach my $step (@{$procsteps->{$resource}}) {
         print STDERR "processing $resource ($step->{name})\n" ;
         if (!-s $step->{out_fn} || $overwrite_fl == 1) {
            print "   * $step->{name} ".localtime()."\n" if
               (!exists $in->{quiet_fl});

            &{$step->{processor}}({
               pibase_specs => $pibase_specs,
               in_fn => $step->{in_fn},
               out_fn => $step->{out_fn},
               header_fl => $step->{header_fl}}) ;
         }
      }
      print STDERR "importing external data $resource\n" ;
      import_external_data({pibase_specs => $pibase_specs,
                            status => $status,
                            resource => $resource}) ;
   }

   return $status ;

}


=head2 import_external_data()

   Title:       import_external_data()
   Function:    Imports processed external data into PIBASE
   Args:        ->{pibase_specs} = pibase_specs data
                ->{resource} = resource information from fetch_external_data()
                ->{status} = resource file status info from fetch_external_data()
   Return:      Nothing

=cut

sub import_external_data {
   my $in = shift ;

   my $pibase_specs = $in->{pibase_specs} ;
   my $resource = $in->{resource} ;
   my $status = $in->{status} ;

   my ($temp_fh, $temp_fn) ;
   ($temp_fh->{external_data_sources}, $temp_fn->{external_data_sources}) =
       tempfile("external_data_sources.XXXXX",
                SUFFIX => $pibase_specs->{pibase_id}, CLEANUP => 1) ;

   foreach my $proc_file (
    keys %{$pibase_specs->{external_data}->{$resource}->{processed_files}}) {
      print localtime()."   * IMPORTING $proc_file\n" if
         (!exists $in->{quiet_fl});
      my $t_status = pibase::mysqlimport({fn =>
$pibase_specs->{external_data}->{$resource}->{processed_files}->{$proc_file},
         pibase_specs => $pibase_specs});

      print "#IMPORT ".join("\t", $proc_file, $t_status->{command_output})."\n";
      if ($t_status->{skipped} > 0 || $t_status->{warnings} > 0) {
         print STDERR "WARNING: mysqlimport of $proc_file: ".
            $t_status->{command_output}."\n" ;
      }
      $status->{$resource}->{mysqlimport}->{$proc_file} = $t_status ;

      print {$temp_fh->{external_data_sources}} join("\t",
         $resource, $proc_file,
         pibase::get_current_date_mysql(), pibase::get_current_date_mysql(),
         $pibase_specs->{external_data}->{$resource}->{ver},
         $pibase_specs->{external_data}->{$resource}->{base_url}, "")."\n" ;
   }
   close($temp_fh->{external_data_sources}) ;

   pibase::mysqlimport({fn => $temp_fn->{external_data_sources},
                        pibase_specs => $pibase_specs}) ;
   unlink($temp_fn->{external_data_sources}) ;
}


=head2 download_asteroids()

   Title:       download_asteroids()
   Function:    Downloads ASTRAL ASTEROIDS data
   Args:        ->{pibase_specs} = pibase_specs data
   Return:      Nothing

=cut

sub download_asteroids {

   my $in = shift;
   my $pibase_specs = $in->{pibase_specs} ;

   system("mkdir -p ".$pibase_specs->{asteroids}->{dir}) ;
   system("mkdir -p ".$pibase_specs->{asteroids}->{dir}.'/sequences') ;
   system("mkdir -p ".$pibase_specs->{asteroids}->{dir}.'/alignments') ;

   my $orig_dir = getcwd() ;
   chdir $pibase_specs->{asteroids}->{dir} ;

   system("wget -r -l1 --no-parent -A.fa http://astral.berkeley.edu/asteroids-1.75/sequences/fam/") ;
   system("mv astral.berkeley.edu/asteroids-1.75/sequences/fam sequences") ;

   system("wget -r -l1 --no-parent -A.fa http://astral.berkeley.edu/asteroids-1.75/sequences/sf/") ;
   system("mv astral.berkeley.edu/asteroids-1.75/sequences/sf sequences") ;

   system("wget -r -l1 --no-parent -A.fasta_aln http://astral.berkeley.edu/asteroids-1.75/alignments/fam/") ;
   system("mv astral.berkeley.edu/asteroids-1.75/alignments/fam alignments") ;

   system("wget -r -l1 --no-parent -A.fasta_aln http://astral.berkeley.edu/asteroids-1.75/alignments/sf/") ;
   system("mv astral.berkeley.edu/asteroids-1.75/alignments/sf alignments") ;

   system("rm -rf astral.berkeley.edu") ;
   chdir $orig_dir ;

}

1;
