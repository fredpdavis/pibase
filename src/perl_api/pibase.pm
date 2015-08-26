=head1 NAME

pibase.pm - perl interface to the pibase database

=head1 DESCRIPTION

The pibase.pm perl library contains subroutines used to build and access the
PIBASE database.

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

package pibase ;
use Exporter ;

use strict ;
use warnings ;
use Sys::Hostname qw/hostname/ ;
use pibase::specs qw/table_spec/ ;
use File::Copy qw/move/ ;
use POSIX qw/floor/ ;
use File::Path qw/mkpath/ ;
use File::Temp qw/tempfile/ ;

our $VERSION = "201009" ;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/connect_pibase connect_tod connect_metatod load_bdp_ids mysql_fetchcols mysql_hashload mysql_hashindload mysql_hasharrload timestamp timestampsec replace_char replace_undefs mysql_createtable locate_binaries get_specs safe_move safe_copy complete_pibase_specs $pibase_specs sid_2_domdir mysql_singleval/ ;

my $pibase_specs = {
   db => 'pibasemysqldatabasename',
   host => 'localhost',
   user => 'pibaseusername' ,
   pass => 'pibasepassword',
   root => '/groups/eddy/home/davisf/work/pibase/pibase2010',
} ;

{
   
## from Janin 1998 JMB captures H20-mediated
   $pibase_specs->{interface_detect_calc}->{dist_cutoff} = 6.05 ;

   $pibase_specs->{pibase_version} = $VERSION ;
   $pibase_specs->{pibase_id} = 'pibase'.$pibase_specs->{pibase_version} ;
   $pibase_specs->{dataroot} = $pibase_specs->{root}.'/data' ;
   $pibase_specs->{externaldata_dir} = $pibase_specs->{root}.'/external_data' ;

   $pibase_specs->{external_data}->{scop} = {
      ver => '1.75',
      base_url => 'http://scop.mrc-lmb.cam.ac.uk/scop',
      cgi_base_url => 'http://scop.mrc-lmb.cam.ac.uk/scop',
      processed_files => {
         scop_cla => 0,
         scop_des => 0,
         scop_hie => 0,
         "subsets.scop" => 0,
         "subsets_class.scop" => 0,
         "subsets_details.scop" => 0,
      },
      files => {
         cla => 'dir.cla.scop.txt_1.75',
         des => 'dir.des.scop.txt_1.75',
         hie => 'dir.hie.scop.txt_1.75',
      },
      urls => {
       'http://scop.mrc-lmb.cam.ac.uk/scop/parse/dir.hie.scop.txt_1.75' => 0,
       'http://scop.mrc-lmb.cam.ac.uk/scop/parse/dir.cla.scop.txt_1.75' => 0,
       'http://scop.mrc-lmb.cam.ac.uk/scop/parse/dir.des.scop.txt_1.75' => 0,
      }
   } ;

   $pibase_specs->{external_data}->{astral} = {
      ver => '1.75',
      urls => {
         'http://astral.berkeley.edu/scopseq-1.75.tgz' => 0,
         'http://astral.berkeley.edu/raf/astral-rapid-access-1.75.raf' => 0,
      },
      postcommands => {
         'http://astral.berkeley.edu/scopseq-1.75.tgz' =>
            "tar xvfz scopseq-1.75.tgz"
      },
   } ;

   $pibase_specs->{external_data}->{cath} = {
      ver => '3.3.0',
      base_url => "http://www.cathdb.info",
      files => {
         CathNames => 'CathNames.v3.3.0',
         CathDomainList => 'CathDomainList.v3.3.0',
         CathDomainDescription => 'CathDomainDescriptionFile.v3.3.0',
      },
      processed_files => {
         cath_names => 0,
         cath_domain_description => 0,
         cath_domain_list => 0,
         "subsets.cath" => 0,
         "subsets_class.cath" => 0,
         "subsets_details.cath" => 0,
      },
      urls => {
'http://release.cathdb.info/v3.3.0/CathChainList.v3.3.0' => 0,
'http://release.cathdb.info/v3.3.0/CathDomainList.v3.3.0' => 0,
'http://release.cathdb.info/v3.3.0/CathNames.v3.3.0' => 0,
'http://release.cathdb.info/v3.3.0/CathDomall.v3.3.0' => 0,
'http://release.cathdb.info/v3.3.0/CathDomainDescriptionFile.v3.3.0' => 0,
      }
   } ;

   $pibase_specs->{external_data}->{pdb} = {
      ver => '20100906',
      base_url => 'http://www.wwpdb.org',
      file_layout => 'wwpdb', #alternative: onedir
      compress_fl => 1, #0 = not compressed
      files => {
         entries => 'entries.idx',
         pdb_entry_type => 'pdb_entry_type.txt',
         obsolete => 'obsolete.dat',
      },
      processed_files => {
         pdb_entries => 0,
         pdb_entry_type => 0, #note, still requires putting a pibase_id suffix before upload
         pdb_obsolete => 0,
         pdb_release => 0,
      },
      urls => {
         'ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt' => 0,
         'ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx' => 0,
         'ftp://ftp.wwpdb.org/pub/pdb/data/status/obsolete.dat' => 0,
      }
   } ;

#   $pibase_specs->{external_data}->{pqs} = {
#      ver => '20080301',
#      base_url => "http://pqs.ebi.ac.uk",
#      file_layout => 'wwpdb', #alternative: onedir
#      compress_fl => 0,
#      files => {
#         list => 'LIST',
#         asalist => 'ASALIST',
#         biolist => 'BIOLIST',
#         ranking => 'RANKING',
#      },
#      processed_files => {
#         pqs_list => 0,
#         pqs_asalist => 0,
#         pqs_biolist => 0,
#         pqs_ranking => 0,
#      },
#      urls => {
#         'ftp://ftp.ebi.ac.uk/pub/databases/msd/pqs/LIST' => 0,
#         'ftp://ftp.ebi.ac.uk/pub/databases/msd/pqs/ASALIST' => 0,
#         'ftp://ftp.ebi.ac.uk/pub/databases/msd/pqs/BIOLIST' => 0,
#         'ftp://ftp.ebi.ac.uk/pub/databases/msd/pqs/RANKING' => 0,
#      }
#   } ;

   $pibase_specs->{external_data}->{pisa} = {
      ver => '20100906',
      base_url => "ftp://ftp.ebi.ac.uk/pub/databases/msd/pisa",
      file_layout => 'wwpdb', #alternative: onedir
      compress_fl => 1,
      max_num_chains => 100,
      files => {
         "index" => 'index.txt',
      },
      processed_files => {
         pisa_index => 0,
      },
      urls => {
         'ftp://ftp.ebi.ac.uk/pub/databases/msd/pisa/index.txt' => 0,
      }
   } ;

   foreach my $resource (keys %{$pibase_specs->{external_data}}) {
      my $dirname = $pibase_specs->{externaldata_dir}.'/'.$resource.'/'.
         $pibase_specs->{external_data}->{$resource}->{ver} ;

      $pibase_specs->{external_data}->{$resource}->{local_dir} = $dirname ;
      foreach my $file (keys
         %{$pibase_specs->{external_data}->{$resource}->{files}}) {
         $pibase_specs->{external_data}->{$resource}->{files}->{$file}=$dirname.
            '/'.$pibase_specs->{external_data}->{$resource}->{files}->{$file} ;
      }

      foreach my $proc_file (keys
         %{$pibase_specs->{external_data}->{$resource}->{processed_files}}) {
        $pibase_specs->{external_data}->{$resource}->{processed_files}->{$proc_file}=
          $pibase_specs->{external_data}->{$resource}->{local_dir}."/$proc_file.".
          $pibase_specs->{pibase_id} ;
      }
   }


   $pibase_specs->{web}->{html_dir} = "/var/www/websites/pibase/html/pibase2010/" ;
   $pibase_specs->{web}->{cgi_dir} = "/var/www/websites/pibase/cgi-bin/pibase2010/" ;

   $pibase_specs->{web}->{base_url} = "http://localhost/pibase/pibase2010";
   $pibase_specs->{web}->{basecgi_url} = "http://localhost/pibase-cgi/pibase2010";
   $pibase_specs->{web}->{domain_colors} = [qw/black cyan4 DarkOrange MediumPurple2 blue grey SeaGreen PowderBlue RosyBrown beige brown coral sienna orchid plum bisque PeachPuff khaki OliveDrab SaddleBrown LimeGreen turquoise SteelBlue CornflowerBlue DeepSkyBlue aquamarine chocolate SpringGreen IndianRed DarkViolet maroon/] ;
   $pibase_specs->{web}->{bdp_topology_graphs_baseurl} = 
      $pibase_specs->{web}->{base_url}.'/data_files/bdp_topology_graphs' ;
   $pibase_specs->{web}->{subsets_files_baseurl} = 
      $pibase_specs->{web}->{base_url}.'/data_files/subsets_files' ;
   $pibase_specs->{web}->{subsets_files_basedir} = 
      $pibase_specs->{web}->{html_dir}.'/data_files/subsets_files' ;

   $pibase_specs->{web_pilig}->{html_dir} = "/var/www/websites/pibase/html/pibase2010/" ;
   $pibase_specs->{web_pilig}->{cgi_dir} = "/var/www/websites/pibase/cgi-bin/pibase2010/" ;

   $pibase_specs->{web_pilig}->{base_url} = "http://localhost/pibase/pibase2010";
   $pibase_specs->{web_pilig}->{basecgi_url} = "http://localhost/pibase-cgi/pibase2010";
   $pibase_specs->{web_pilig}->{domain_colors} = [qw/black cyan4 DarkOrange MediumPurple2 blue grey SeaGreen PowderBlue RosyBrown beige brown coral sienna orchid plum bisque PeachPuff khaki OliveDrab SaddleBrown LimeGreen turquoise SteelBlue CornflowerBlue DeepSkyBlue aquamarine chocolate SpringGreen IndianRed DarkViolet maroon/] ;
   $pibase_specs->{web_pilig}->{bdp_topology_graphs_baseurl} = 
      $pibase_specs->{web_pilig}->{base_url}.'/data_files/bdp_topology_graphs' ;
   $pibase_specs->{web_pilig}->{subsets_files_baseurl} = 
      $pibase_specs->{web_pilig}->{base_url}.'/data_files/subsets_files' ;
   $pibase_specs->{web_pilig}->{subsets_files_basedir} = 
      $pibase_specs->{web_pilig}->{html_dir}.'/data_files/subsets_files' ;

   $pibase_specs->{pdb_dir}=
      '/groups/eddy/home/davisf/work/databases/pdb/data/structures' ;
   $pibase_specs->{pdbnmr_dir}= $pibase_specs->{externaldata_dir}.
      '/structures/pdb_nmr';
#   $pibase_specs->{pqs_dir}= '/groups/eddy/home/davisf/work/databases/pqs' ;
   $pibase_specs->{pisa_dir}= '/groups/eddy/home/davisf/work/databases/pisa' ;

   $pibase_specs->{astral}->{dir} =  $pibase_specs->{externaldata_dir}.
      '/astral/'.$pibase_specs->{external_data}->{astral}->{ver} ;

   $pibase_specs->{astral}->{gd_seq} = $pibase_specs->{astral}->{dir}.
      '/scopseq-'.$pibase_specs->{external_data}->{astral}->{ver}.
      '/astral-scopdom-seqres-gd-all-'.
      $pibase_specs->{external_data}->{astral}->{ver}.'.fa' ;

   $pibase_specs->{astral}->{raf} = $pibase_specs->{astral}->{dir}.
      '/astral-rapid-access-'.
      $pibase_specs->{external_data}->{astral}->{ver}.
      '.raf' ;

   foreach my $seqid ( qw/100 95 90 70 50 30 20 10/) {
      $pibase_specs->{astral}->{seqcl}->{$seqid} =
         $pibase_specs->{astral}->{dir}.'/scopseq-'.
         $pibase_specs->{external_data}->{astral}->{ver}.'/'.
         '/astral-scopdom-seqres-gd-sel-gs-bib-verbose-'.$seqid.'-'.
         $pibase_specs->{external_data}->{astral}->{ver}.'.txt';
   }

   $pibase_specs->{asteroids}->{dir} = $pibase_specs->{externaldata_dir}.
      '/asteroids/1.75' ;
   $pibase_specs->{asteroids}->{fam_aln} = $pibase_specs->{asteroids}->{dir}.
      '/alignments/fam';
   $pibase_specs->{asteroids}->{sf_aln}=$pibase_specs->{asteroids}->{dir}.
      '/alignments/sf';
   $pibase_specs->{asteroids}->{fam_seq} = $pibase_specs->{asteroids}->{dir}.
      '/sequences/fam' ;
   $pibase_specs->{asteroids}->{sf_seq} = $pibase_specs->{asteroids}->{dir}.
      '/sequences/sf' ;

   $pibase_specs->{subsets_dir} = $pibase_specs->{dataroot}.'/subsets_files';

   $pibase_specs->{buildfiles_dir} = $pibase_specs->{dataroot}.'/build_files';
   foreach my $t_buildfile (qw/bdp_files bdp_residues_tables bdp_chains/,
    qw/bdp_chains_equiv_update subsets subsets_details subsets_residues_tables/,
    qw/intersubset_contacts patch_residues_tables interface_contacts_tables/,
    qw/interface_contacts_special_tables bdp_secstrx_tables subsets_files/,
    qw/subsets_sasa_details subsets_sasa interface_sasa /,
    qw/scop_interface_clusters interface_secstrx_tables/,
    qw/bindingsite_secstrx_basic_contacts_tables bindingsite_contacts_tables/,
    qw/interface_secstrx_basic_contacts_tables/,
    qw/interface_secstrx_contacts_tables interface_secstrx_profile/,
    qw/interface_resvector bdp_interaction_topology/,
    qw/bindingsite_sse_topology interface_sse_topology subsets_sequence/,
    qw/interface_size bdp_interaction_topology_graph/
     ) {
      $pibase_specs->{buildfiles}->{$t_buildfile}=
         $pibase_specs->{buildfiles_dir}.'/'.
         $t_buildfile.'.'.$pibase_specs->{pibase_id} ;
   }

   $pibase_specs->{otherdata_root} = $pibase_specs->{dataroot}.'/other_data' ;
   $pibase_specs->{otherdata_dir} = {
bdp_topology_graphs => $pibase_specs->{otherdata_root}.'/bdp_topology_graphs',
   } ;

   $pibase_specs->{tod_dir} = $pibase_specs->{dataroot}.'/tod';
   $pibase_specs->{metatod_root} = $pibase_specs->{dataroot}.'/metatod' ;
   $pibase_specs->{metatod_dir} = {
      bdp_residues => $pibase_specs->{metatod_root}.'/bdp_residues',
      interatomic_contacts => $pibase_specs->{metatod_root}.
         '/interatomic_contacts',
      interface_contacts_special => $pibase_specs->{metatod_root}.
         '/interface_contacts_special',
      interface_contacts => $pibase_specs->{metatod_root}.'/interface_contacts',
      patch_residues => $pibase_specs->{metatod_root}.'/patch_residues',
      bdp_secstrx => $pibase_specs->{metatod_root}.'/bdp_sec_strx',
      subsets_residues => $pibase_specs->{metatod_root}.'/subsets_residues',
      interface_secstrx => $pibase_specs->{metatod_root}.'/interface_secstrx',
      interface_secstrx_contacts => $pibase_specs->{metatod_root}.
         '/interface_secstrx_contacts',
      interface_secstrx_basic_contacts => $pibase_specs->{metatod_root}.
         '/interface_secstrx_basic_contacts',
      bindingsite_secstrx_basic_contacts => $pibase_specs->{metatod_root}.
         '/bindingsite_secstrx_basic_contacts',
      bindingsite_contacts => $pibase_specs->{metatod_root}.
         '/bindingsite_contacts',
   } ;
   $pibase_specs->{binaries} = locate_binaries() ;

   $pibase_specs->{SGE}->{cluster_mode} = 1 ;
   $pibase_specs->{SGE}->{numjobs} = 100 ;
   $pibase_specs->{SGE}->{qstat_sleep} = 120 ;

}


=head2 dbname()

   Title:       dbname()
   Function:    gets the name of the pibase database
   Args:	pibase $specs data structure
   Returns:     returns the name of the pibase mysql database

=cut

sub dbname {

   my $specs = shift ;
   $specs = complete_pibase_specs() ;

   return $specs->{db} ;

}


=head2 get_specs()

   Title:       get_specs()
   Function:    gives pibase specifications
   Args:	pibase $specs data structure
   Returns:     returns completed pibase data specifications

=cut

sub get_specs {

   my $specs = $pibase_specs ;
   return $specs ;

}


=head2 connect_pibase($dbspecs)

   Title:       connect_pibase()
   Function:    Connects to the pibase database.
   Args:        $_->{db}	database name
                $_->{user}	user name
                $_->{pass}	password
   Returns:     DBI database handle to pibaes

=cut

sub connect_pibase {

   require DBI ;

   my $specs = shift ;
   $specs = complete_pibase_specs($specs) ;

   my $curhost = hostname()  ;
   my $dbname ;

   if ($specs->{host} =~ /^$curhost/) {
      $dbname = "DBI:mysql:database=$specs->{db}";
   } else {
      $dbname = "DBI:mysql:database=$specs->{db}".';'."host=$specs->{host}"; } 

   if (exists $specs->{mysql_socket}) {
      $dbname.=';mysql_socket='.$specs->{mysql_socket};}

   my $dbh = DBI->connect( $dbname, $specs->{user}, $specs->{pass},
                           {RaiseError => 1, AutoCommit => 1} ) ;

   return ($dbh, $specs->{db}, $specs->{root}) ;

}


=head2 connect_tod($dbspecs)

   Title:       connect_tod()
   Function:    Connects to pibase tables on disk (POD).
   Args:        $_[0] = tablename
                $_[1] = dbspecs
                $_[2] = tod_dir = location of tables on disk

   Returns:     DBI database handle to table on disk

=cut

sub connect_tod {

   my $table = shift ;
   my $specs = shift ;

   $specs = complete_pibase_specs() ;

   my $tablespecs = pibase::specs::table_spec($table) ;
   my $format= join( ',', @{$tablespecs->[0]->{specs}->{field_name}});

   my $dbh = DBI->connect('dbi:AnyData(RaiseError=>1):') ;
   $dbh->func( $table, 'Tab', $specs->{tod_dir}.'/'.$table,
               {col_names => $format} , 'ad_import') ;

   return $dbh ;

}


=head2 rawselect_tod()

   Title:       rawselect_tod()
   Function:    performs basic SELECT statement queries on a table on disk
   Args:        $_[0] = SQL SELECT-like statement
                $_[1] = fullfile

=cut

sub rawselect_tod {

   my $select_sql = shift ;
   my $fullfile ;
   if ($#_ >= 0) {
      $fullfile = shift ; }
   my $specs = complete_pibase_specs() ;
   my $bins = locate_binaries() ;

   my ($fields, $table) =
      ($select_sql =~ /(?:SELECT|select) (.+) (?:FROM|from) (\w+)$/) ;

   $fields =~ s/ //g ;
   my @fields = split(/\,/, $fields) ;

   my $tablespecs = pibase::specs::table_spec($table) ;
   my $field_rev ;
   foreach my $j (0 .. $#{$tablespecs->[0]->{specs}->{field_name}}) {
      $field_rev->{$tablespecs->[0]->{specs}->{field_name}->[$j]} = $j ; }

   my @fields_id ;
   foreach my $j ( 0 .. $#fields) {
      if (!exists $field_rev->{$fields[$j]}) {
         print STDERR "ERROR rawselect_tod(): field $fields[$j] doesn't exist\n";
	 return "ERROR" ;
      }
      $fields_id[$j] = $field_rev->{$fields[$j]} ;
   }

   if (!defined $fullfile) {
      $fullfile = $specs->{tod_dir}.'/'.$table ;
      if (!-s $fullfile) {
         print STDERR "ERROR rawselect_tod(): $fullfile not found\n" ;
         return "ERROR" ;
      }
   }

   if ($fullfile =~ /gz$/) { $fullfile = "$bins->{zcat} $fullfile |" ; }
   my @results ;
   open(TABLEF, $fullfile) ;
   while (my $line = <TABLEF>) {
      chomp $line;
      my @curf = split(/\t/, $line) ;
      foreach my $j ( 0 .. $#fields_id) {
         push @{$results[$j]}, $curf[$fields_id[$j]] ; }
   }
   close(TABLEF) ;

   return @results ;
}


=head2 connect_metatod()

   Title:       connect_metatod()
   Function:    performs basic SELECT statement queries on a meta-table on disk
   Args:        $_[0] = filename
                $_[1] = tablename
                $_[2] = pibase db specs
   Returns:     dbh DBI:AnyData database handle

=cut

sub connect_metatod {

   my $filename = shift ;
   my $tablename = shift ;

   my $specs = shift ;
   $specs = complete_pibase_specs() ;

   my $proto = $filename;
   $proto =~ s/.*\/// ;
   $proto =~ s/_[0-9]+.*/_prototype/ ;

   my $tablespecs = pibase::specs::table_spec($proto) ;
   my $format = join( ',', @{$tablespecs->[0]->{specs}->{field_name}} ) ;

   my $dbh = DBI->connect('dbi:AnyData(RaiseError=>1):') ;
   $dbh->func( $tablename, 'Tab', $filename,
               {col_names => $format} , 'ad_import') ;

   return $dbh ;

}


=head2 rawselect_metatod()

   Title:       rawselect_metatod()
   Function:    selects specified fields from a table-on-disk.
      Note: WHERE clause does not work, this command just recognizes
            the field names and returns the appropriate columns.
   Args:        $_[0] = filename
                $_[1] = SELECT sql command
   Returns:     array of query results

=cut


sub rawselect_metatod {

   my $filename = shift ;
   my $select_sql = shift ;
   my $specs = complete_pibase_specs() ;
   my $bins = locate_binaries() ;

   my $proto = $filename;
   $proto =~ s/.*\/// ;
   $proto =~ s/_[0-9]+.*// ;
   $proto =~ s/\.[0-9]+.*// ;
   $proto .= '_prototype' ;

   if ($proto eq 'secstrx_prototype') {
      $proto = 'bdp_secstrx_prototype' ; }

   if ($proto eq 'secstrx_basic_prototype') {
      $proto = 'bdp_secstrx_basic_prototype' ; }

   my $tablespecs = pibase::specs::table_spec($proto) ;

   my $fields = $select_sql ;
   $fields =~ s/(FROM|from).*$// ;
   $fields =~ s/^(SELECT|select)// ;
   $fields =~ s/ //g ;
   my @fields = split(/\,/, $fields) ;

   my $field_rev ;
   foreach my $j (0 .. $#{$tablespecs->[0]->{specs}->{field_name}}) {
      $field_rev->{$tablespecs->[0]->{specs}->{field_name}->[$j]} = $j ; }

   my @fields_id ;
   foreach my $j ( 0 .. $#fields) {
      if (!exists $field_rev->{$fields[$j]}) {
         print STDERR "rawselect_metatod() error: field $fields[$j] doesn't exist\n" ;
	 return "ERROR" ;
      }
      $fields_id[$j] = $field_rev->{$fields[$j]} ;
   }

   my @results ;
   if (!-s $filename) {
      print STDERR "ERROR rawselect_metatod(): $filename not found\n" ;
      return "ERROR" ;
   }

   if ($filename =~ /gz$/) {
      $filename = "$bins->{zcat} $filename |" ; }

   open(TABLEF, $filename) ;
   while (my $line = <TABLEF>) {
      chomp $line;
      my @curf = split(/\t/, $line) ;
      foreach my $j ( 0 .. $#fields_id) {
         push @{$results[$j]}, $curf[$fields_id[$j]] ; }
   }
   close(TABLEF) ;

   return @results ;

}



=head2 sid_2_domdir()

   Title:       sid_2_domdir()
   Function:    returns the directory name where the PDB file of the
      specified domain resides.
   Args:        $_ = subset_id
   Returns:     directory name

=cut

sub sid_2_domdir {

   my $sid = shift ;
   my ($bdp) = ($sid =~ /BDP([0-9]+)/) ;
   my $dirnum = POSIX::floor($bdp / 100) ;
   my $dir = $pibase_specs->{subsets_dir}."/$dirnum/$bdp" ;

   return $dir ;

}


=head2 complete_pibase_specs(specs)

   Title:       complete_pibase_specs
   Function:    Fills in blanks in specs with default values.
   Input:       $_ = specs hashref
                $_->{db} = database_name
                $_->{user} = user name
                $_->{pass} = password
   Return:      specs - hashref

=cut

sub complete_pibase_specs {

   my $specs = shift  ;

   my @keys = keys %{$pibase_specs} ;

   foreach my $key (@keys) {
      if (!exists $specs->{$key}) {
         $specs->{$key} = $pibase_specs->{$key} ; } }

   return $specs ;

}


=head2 load_bdp_ids($dbh, @results_type)

   Name: load_bdp_ids()
   Function:    Returns bdp_id and depending on results_type
                specified, its relation to bdp_path and pdb_id
                in a variety of forms.
   Return:      results
   Args:        $_[0] = DBI dbh handle
                $_[1] = results type

=over

=item * path_2_bdp_id (hash) [default]

=item * bdp_id_2_path (hash)

=item * bdp_id_2_pdb_id (hash)

=item * bdp_id_2_raw_pdb (hash)

=item * pdb_id_2_bdp_id (hash)

=item * bdp_id (array)

=back

=cut

sub load_bdp_ids {

   my $dbh = shift ;

   my @result_types = @_ ;

   if ($#result_types < 0) {
      push @result_types, 'path_2_bdp_id' ; }

   my @results ;

   foreach my $j ( 0 .. $#result_types) {
      my $ans ;
      if ($result_types[$j] eq 'path_2_bdp_id') {
	 $ans = mysql_hashload($dbh, 'SELECT file_path, bdp_id FROM bdp_files');
      } elsif ($result_types[$j] eq 'bdp_id_2_path') {
	 $ans = mysql_hashload($dbh, 'SELECT bdp_id, file_path FROM bdp_files');
      } elsif ($result_types[$j] eq 'bdp_id_2_pdb_id') {
	 $ans = mysql_hashload($dbh, 'SELECT bdp_id, pdb_id FROM bdp_files');
      } elsif ($result_types[$j] eq 'bdp_id_2_raw_pdb') {
	 $ans = mysql_hashload($dbh, 'SELECT bdp_id, raw_pdb FROM bdp_files');
      } elsif ($result_types[$j] eq 'pdb_id_2_bdp_id') {
	 $ans = mysql_hashload($dbh, 'SELECT pdb_id, bdp_id FROM bdp_files WHERE raw_pdb = 1');
      } elsif ($result_types[$j] eq 'bdp_id') {
	 ($ans) = mysql_fetchcols($dbh, 'SELECT bdp_id FROM bdp_files') ;
      }
      push @results, $ans ;
   }

   return @results ;

}


=head2 todload_bdp_ids(@results_type)

   Name: todload_bdp_ids()
   Function:    Returns bdp_id and depending on results_type
                specified, its relation to bdp_path and pdb_id
                in a variety of forms.
      Analogous to load_bdp_ids() with tables-on-disk instead of DBI
   Return:      query results
   Args:        $_[0] = DBI dbh handle
                $_[1] = results type

=over

=item * path_2_bdp_id (hash) [default]

=item * bdp_id_2_path (hash)

=item * bdp_id_2_pdb_id (hash)

=item * bdp_id_2_raw_pdb (hash)

=item * pdb_id_2_bdp_id (hash)

=item * bdp_id (array)

=back

=cut

sub todload_bdp_ids {

   my @result_types = @_ ;

   my $tablespecs = pibase::specs::table_spec("bdp_files") ;
   my ($bdp_id, $file_path, $file_base, $pdb_id, $raw_pdb)=
      rawselect_tod("SELECT bdp_id, file_path, file_base, pdb_id, raw_pdb ".
                    "FROM bdp_files") ;

   if ($#result_types < 0) {
      push @result_types, 'path_2_bdp_id' ; }

   my @results ;

   foreach my $j ( 0 .. $#result_types) {

      my $ans ;
      if ($result_types[$j] eq 'path_2_bdp_id') {
         foreach my $k ( 0 .. $#{$bdp_id}) {
            $ans->{$file_path->[$k]} = $bdp_id->[$k] ; }
      } elsif ($result_types[$j] eq 'bdp_id_2_path') {
         foreach my $k ( 0 .. $#{$bdp_id}) {
            $ans->{$bdp_id->[$k]} = $file_path->[$k] ; }
      } elsif ($result_types[$j] eq 'bdp_id_2_pdb_id') {
         foreach my $k ( 0 .. $#{$bdp_id}) {
            $ans->{$bdp_id->[$k]} = $pdb_id->[$k] ; }
      } elsif ($result_types[$j] eq 'bdp_id_2_raw_pdb') {
         foreach my $k ( 0 .. $#{$bdp_id}) {
            $ans->{$bdp_id->[$k]} = $raw_pdb->[$k] ; }
      } elsif ($result_types[$j] eq 'pdb_id_2_bdp_id') {
         foreach my $k ( 0 .. $#{$bdp_id}) {
            if ($raw_pdb->[$k] == 1) {
               $ans->{$pdb_id->[$k]} = $bdp_id->[$k] ; }}
      } elsif ($result_types[$j] eq 'bdp_id') {
         push @{$ans}, @{$bdp_id} ;
      }

      push @results, $ans ;
   }

   return @results ;

}


=head2 mysql_fetchcols(dbh, query)

   Function:    Processes an n column query and returns a list of array references, where each array reference holds all the valeus for a given column.
   Args:        $_[0] = dbh - DBI database handle
                $_[1] = SQL select query

   Returns:     @_ - list of arrayref
                $a[i]->[j] ith column, jth row

=cut

sub mysql_fetchcols {

   my $dbh = shift ;
   my $query = shift ;
   my $vals = shift ;

   my $sth ;
   if ($query->can("execute")) {
      $sth = $query ;
      $sth->execute(@{$vals}) ;
   } else {
      $query = qq{$query} ;
      $sth = $dbh->prepare($query) ;
      $sth->execute() ;
   }

   my @results ;
   while (my @currow = $sth->fetchrow()) {
      foreach my $j (0 .. $#currow) {
         push @{$results[$j]}, $currow[$j] ; } }

   return @results ;

}


=head2 mysql_hashindload(dbh, query)

   Name:        mysql_hashindload() ;
   Function:    Processes a 1 column query and returns a hash pointing from
                  column1 values to row number.
   Args:        $_[0] = dbh - DBI database handle
                $_[1] = SQL SELECT query
                $_[2] = substitutor for undefined value
   Returns:	$_ = hashref
                $a->{col1} = row number

=cut

sub mysql_hashindload {

   my $dbh = shift;
   my $query = shift;
   my $undef_sub = shift;

   $query = qq{$query} ;
   my $sth = $dbh->prepare($query) ;
   $sth->execute() ;

   my $ans ;
   my $j = 0;
   while ( my @currow = $sth->fetchrow() ) {
      if ((! defined $currow[0]) && (defined $undef_sub)) { $currow[0] = $undef_sub; }
      $ans->{$currow[0]} = $j ;
      $j++ ;
   }

   return $ans ;

}


=head2 array2hash(arrayref, undef_sub)

   Name:        array2hash() ;
   Function:    Takes an array reference and returns a hashref with value,
                index pairs.
   Args:        $_[0] = array references
                $_[1] = undefined substitution - if the cell contains an
                  undefined value, use this value as the hash key
   Returns:	$_ = hashref
                $a->{value} = row number

=cut

sub array2hash {

   my $array = shift ;
   my $undef_sub = shift;

   my $ans ;
   foreach my $j ( 0 .. $#{$array}) {
      my $t = $array->[$j] ;
      if ((! defined $array->[$j]) && (defined $undef_sub)) { $t = $undef_sub; }
      $ans->{$t} = $j ;
   }

   return $ans ;

}


=head2 replace_undefs(arrayref, undef_sub)

   Name:        array2hash() ;
   Function:    Takes an array reference and replaces (inplace) undefined
                  values to a specified substitution value.
   Args:        $_[0] = array references
                $_[1] = undefined substitution - if the cell contains an
                  undefined value, replace with this value
   Returns:	$_ = array reference

=cut

sub replace_undefs {

   my $array = shift ;
   my $undef_sub = shift;

   foreach my $j ( 0 .. $#{$array}) {
      if (! defined $array->[$j]) { $array->[$j] = $undef_sub; }}

   return $array ;

}


=head2 replace_undefs_blanks(arrayref, undef_sub)

   Name:        array2hash() ;
   Function:    Takes an array reference and replaces (inplace) undefined and
                  blank values to a specified substitution value.
   Args:        $_[0] = array references
                $_[1] = undefined/blank substitution - if the cell
                  contains an undefined or blank value, replace with this value
   Returns:	$_ = array reference

=cut

sub replace_undefs_blanks {

   my $array = shift ;
   my $undef_sub = shift;

   foreach my $j ( 0 .. $#{$array}) {
      if (! defined $array->[$j] || $array->[$j] eq '') {
         $array->[$j] = $undef_sub; }}

   return $array ;

}


=head2 replace_char(arrayref, target, replacement)

   Name:        array2hash() ;
   Function:    Takes an array reference and replaces (inplace) target values to a specified substitution value.
   Args:        $_[0] = array reference
                $_[1] = target value
                $_[2] = substitution value
   Returns:	$_ = array reference

=cut

sub replace_char {

   my $array = shift ;
   my $target = shift;
   my $replacement = shift;

   foreach my $j ( 0 .. $#{$array}) {
      if ($array->[$j] eq $target) { $array->[$j] = $replacement; } }

   return $array ;

}



=head2 mysql_hashload(dbh, query)

   Name:        mysql_hashload()
   Function:    Processes a 2 column query and returns a hash pointing from
                  column1 values to column2 values (use for 1:1 relationships)
   Args:        $_[0] = dbh - DBI database handle
                $_[1] = SQL SELECT command
   Returns:     $a - hashref
                $a->{col1} = col2

=cut

sub mysql_hashload {

   my ($dbh, $query) = @_ ;

   $query = qq{$query} ;
   my $sth = $dbh->prepare($query) ;
   $sth->execute() ;

   my $ans ;
   while ( my @currow = $sth->fetchrow() ) {
      $ans->{$currow[0]} = $currow[1] ; }

   return $ans ;

}


=head2 mysql_hash2load(dbh, query)

   Name:        mysql_hash2load()
   Function:    Processes a 3 column query and returns a hash pointing from
                  column1 values to column2 to column3 values
                  (use for 1:1:1 relationships)
   Args:        $_[0] = dbh - DBI database handle
                $_[1] = SQL SELECT command
   Returns:     $a - hashref
                $a->{col1}->{col2} = col3 ;

=cut


sub mysql_hash2load {

   my ($dbh, $query) = @_ ;

   $query = qq{$query} ;
   my $sth = $dbh->prepare($query) ;
   $sth->execute() ;

   my $ans ;
   while ( my @currow = $sth->fetchrow() ) {
      $ans->{$currow[0]}->{$currow[1]} = $currow[2] ; }

   return $ans ;

}

=head2 mysql_hasharrload(dbh, query, undef_subs)

   Name:        mysql_hashload()
   Function:    Processes a 2 column query and returns a hash pointing from
                  column1 values to an array of column2 values
                  (use for 1:n relationships)
   Args:        $_[0] dbh - DBI database handle
                $_[1] query - SQL format
                $_[2] undefined value substitutors - list
   Returns:     $a - hashref to arrayrefs
   $a->{col1} = [ col2_1, col2_2,... col2_n ]

=cut

sub mysql_hasharrload {

   my $dbh = shift;
   my $query = shift;

   my $undef_sub = shift;

   $query = qq{$query} ;
   my $sth = $dbh->prepare($query) ;
   $sth->execute() ;

   my $ans ;
   while ( my @currow = $sth->fetchrow() ) {
      if ((! defined $currow[0]) && (defined $undef_sub->[0])) {
         $currow[0] = $undef_sub->[0] ; }
      if ((! defined $currow[1]) && (defined $undef_sub->[1])) {
         $currow[1] = $undef_sub->[1] ; }
      push @{$ans->{$currow[0]}}, $currow[1] ;
   }

   return $ans ;

}


=head2 mysql_singleval(dbh, query)

   Function:    Processes a 1 column, 1 row query and returns a scalar
                  containing the value.
   Args:        $_[0] = dbh - DBI database handle
                $_[1] = query - SQL format
   Return:      $a - scalar
                $a = "result"

=cut

sub mysql_singleval {

   my ($dbh, $query) = @_ ;

   $query = qq{$query} ;
   my $sth = $dbh->prepare($query) ;
   $sth->execute() ;

   my $ans = $sth->fetchrow() ;

   return $ans ;

}


=head2 timestamp()

   Function:    Returns a timestamp
   Args:        none
   Return:      $_[0] = timestamp: <4-digit YEAR><2-digit MONTH><2-digit DAY>_
                  <2-digit HOUR><2-digit MINUTE>

=cut

sub timestamp {

   my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);

   my $f_day = '0'x(2 - length($mday)).$mday;
   my $f_mon = '0'x(2 - length(($mon + 1))).($mon + 1);
   if ($year > 100) {
      $year = $year - 100 ; }
   my $f_year = '0'x(2-length($year)).$year;

   my $f_hour = '0'x(2 - length($hour)).$hour;
   my $f_min = '0'x(2 - length($min)).$min;


   my $time_stamp = $f_year.$f_mon.$f_day.'_'.$f_hour.$f_min ;
   return $time_stamp ;

}


=head2 get_current_date_mysql()

   Function:    Returns current date in YYYY-MM-DD mysql format
   Args:        none
   Return:      $_[0] = date: <4-digit YEAR>-<2-digit MONTH>-<2-digit DAY>

=cut

sub get_current_date_mysql {
   my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);

   my $f_day = '0'x(2 - length($mday)).$mday;
   my $f_mon = '0'x(2 - length(($mon + 1))).($mon + 1);
   if ($year > 100) {
      $year = $year - 100 ; }
   my $f_year = '0'x(2-length($year)).$year;

   my $time_stamp = $f_year.'-'.$f_mon.'-'.$f_day ;
   return $time_stamp ;
}


=head2 timestampsec()

   Function:    Returns a second-resolution timestamp
   Args:        none
   Return:      timestamp: <4-digit YEAR><2-digit MONTH><2-digit DAY>_
                  <2-digit HOUR><2-digit MINUTE><2-digit SECOND>

=cut

sub timestampsec {

   my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);

   my $f_day = '0'x(2 - length($mday)).$mday;
   my $f_mon = '0'x(2 - length(($mon + 1))).($mon + 1);
   if ($year > 100) {
      $year = $year - 100 ; }
   my $f_year = '0'x(2-length($year)).$year;

   my $f_hour = '0'x(2 - length($hour)).$hour;
   my $f_min = '0'x(2 - length($min)).$min;


   my $time_stamp = $f_year.$f_mon.$f_day.'_'.$f_hour.$f_min."_$sec" ;
   return $time_stamp ;

}

=head2 mysqlimport(file,dbspecs)

   Function:    load a file into a mysql database using system(mysqlimport)
   Return:	$_[0] = records imported
                $_[1] = records deleted
                $_[2] = records skipped
                $_[3] = number of warnings

   Args:        $_[0] = filename
                $_[1] = database specs - hashref

=over

=item * db => database name

=item * user => user name

=item * pass => password

=back

=cut

sub mysqlimport {

# mysqlimport -L <database> -u <username> -p <wordcode> <datafile>

   my $in = shift ;
   my $file = $in->{fn} ;
   my $specs = $in->{pibase_specs} ;
   $specs = complete_pibase_specs($specs) ;

   my $import_comm = "mysqlimport -L $specs->{db} -u $specs->{user} -p$specs->{pass} $file" ;
   my $output = `$import_comm` ; chomp $output ;

#   print STDERR "Output: $output e\n";

# Parse mysqlimport statistics: number of records imported, deleted, skipped, and number of warnings.

   my ($records) = ($output =~ /Records:\s?([0-9]+)/ ) ;
   my ($deleted) = ($output =~ /Deleted:\s?([0-9]+)/ ) ;
   my ($skipped) = ($output =~ /Skipped:\s?([0-9]+)/ ) ;
   my ($warnings) = ($output =~ /Warnings:\s?([0-9]+)/ ) ;

   return {
      command_output => $output,
      records => $records,
      skipped => $skipped,
      warnings => $warnings,
      deleted => $deleted
   } ;

}


=head2 mysql_runcom(dbh, query, vaues)

   Function:    Runs an non-SELECT SQL statement
   Return:	none
   Args:        $_[0] = dbh - DBI database handle
                $_[1] = SQL query command or DBI statement handle
                $_[2] = values- arrayref that hold values for '?' holders
                        in DBI statement handle

=cut

sub mysql_runcom {

   my $dbh = shift ;
   my $sql = shift ;
   my $vals = shift ;

   my $sth ;
   if ($sql->can("execute")) {
      $sth = $sql;
      $sth->execute(@{$vals}) ;
   } else {
      $sql= qq{$sql} ;
      $sth = $dbh->prepare($sql) ;
      $sth->execute() ;
   }

   return 1 ;
}


=head2 mysql_createtable(dbh, tablename, spec)

   Function:    Creates a table (NOTE: if table already exists, drops it first).
   Return:	none
   Args:        $_[0] = dbh - DBI database handle
                $_[1] = table name
                $_[2] = spec - SQL DDL string

   Example uasge:

      pibase::mysql_createtable($dbh,
         $tables->{interface_contacts}->{name},
         $tables->{interface_contacts}->{spec}) ;

      pibase::mysql_runcom($dbh,
         "REPLACE INTO $tables->{interface_contacts}->{meta} ".
         "( bdp_id, table_name) values($bdp_id, ".
         "\"$tables->{interface_contacts}->{name}\")") ;

=cut

sub mysql_createtable {


   my $dbh  = shift ;
   my $name = shift ;
   my $spec = shift ;

   mysql_runcom($dbh, "DROP TABLE IF EXISTS $name") ;
   mysql_runcom($dbh, "CREATE TABLE $name $spec") ;

   return 1 ;
}


=head2 mysql_commandline_query(dbh, query, vaues)

   Function:    Runs an SQL SELECT statement on the command line,
                optionally performs a sort (on command line), and 
                displays the output to a specified file
   Return:	none
   Args:        $->{db_name} = database name
                $->{sql} = SQL query
                $->{out_fn} = file to display results to
                $->{post_sort} = optional specify field to sort
                $->{sort_order} = optional sort order (defaults ASC)

=cut

sub mysql_commandline_query {

   my $in = shift ;
   my $sql = $in->{sql} ;
   my $out_fn = $in->{out_fn} ;

   my $db_specs ;
   foreach my $t (qw/db host user pass/) {
      if (!exists $in->{$t}) {
         $db_specs->{$t} = $pibase_specs->{$t} ;
      } else {
         $db_specs->{$t} = $in->{$t} ;
      }
   }

   my ($t_sql_fh, $t_sql_fn) = tempfile() ;
   print $t_sql_fh $sql ;
   if ($in->{sql} !~ /\;$/) {
      print $t_sql_fh ';';}
   print $t_sql_fh "\n" ;
   close($t_sql_fh) ;

   my ($t_out_fh, $t_out_fn) = tempfile(); close($t_out_fh) ;
   my $tcom_mysql = "cat $t_sql_fn | mysql ".$db_specs->{db}.
      " -u ".$db_specs->{user}.
      " -p".$db_specs->{pass}.
      " -h ".$db_specs->{host}." | sed 1d > $t_out_fn" ;
   my $status1 = system($tcom_mysql) ;
   if ($status1 == -1) { die "system call failed: $tcom_mysql"; }

   if (exists $in->{post_sort}) {
      my $tcom_sort = "sort ".$in->{post_sort}." -t\"	\" $t_out_fn > $out_fn";
      my $status2 = system($tcom_sort) ;
      if ($status2 == -1) { die "system call failed: $tcom_sort"; }
   } else {
      my $tcom_move = "mv $t_out_fn $out_fn" ;
      my $status2 = system($tcom_move) ;
      if ($status2 == -1) { die "system call failed"; }
   }

   return 1 ;

}




=head2 locate_binaries()

   Function:    Returns location of binaries used in PIBASE associated activities
   Return:	$_->{program} = program location.
      perl, zcat, rigor, subset_extractor, altloc_check
   Args:        none

=cut

sub locate_binaries {

   use Sys::Hostname qw/hostname/ ;

   my $pibase_specs = complete_pibase_specs() ;
   my $hostname = Sys::Hostname::hostname ;

   my $rootdir = $pibase_specs->{'root'};

   my $binaries ;

   my $mach = 'i386' ;
   $binaries->{'perl'} = "perl" ;
   $binaries->{'zcat'} = "zcat" ;

   $binaries->{'imagemagick_convert'} = "convert" ;
   if ($hostname =~ /^c[0-9][0-9]u[0-9][0-9]/) { #janelia cluster
      $binaries->{'imagemagick_convert'} =
         '/usr/local/ImageMagick-6.2.9/bin/convert' ;
   }

   my $proc_type = `uname -p` ; chomp $proc_type;
   if (!defined $proc_type ) {$proc_type = 'o64' ;}
#fpd080630_0029  weird bug - sometimes uname -p returns nothing??

   if ( ($hostname eq 'alto') || ($hostname eq 'diva') ) {
      $binaries->{'perl'} = 'perl5.6.1' ;
      $binaries->{'zcat'} = 'gzcat' ;
      $mach = 'sun4u' ;
   }

   if ($proc_type eq 'x86_64'||
       $hostname =~ /^o64/ || $hostname =~ /marimba/ ||
       $hostname =~ /^opt/ || $hostname =~ /lyre/) {
      $mach = 'o64' ;
   }

   if ($hostname =~ /^intel/) {
      $mach = 'ia64' ;
   }

   $binaries->{'modeller'} = "/groups/eddy/home/davisf/bin/modeller9v4/bin/mod9v4" ;

   $binaries->{'rigor'} = "$rootdir/auxil/".
      "rigor/rigor.$mach" ;

   $binaries->{'subset_extractor'} = "$rootdir/auxil/".
      "subset_extractor/subset_extractor.$mach" ;

   if (! -e $binaries->{'subset_extractor'}) {
      $binaries->{'subset_extractor'} = "ERROR" ;
   }


   $binaries->{'altloc_check'} = "$rootdir/auxil/".
      "altloc_check/altloc_check.$mach" ;
   if (! -e $binaries->{'altloc_check'}) {
      $binaries->{'altloc_check'} = "ERROR" ;
   }


   $binaries->{'altloc_filter'} = "$rootdir/auxil/".
      "altloc_filter/altloc_filter.pl" ;
   if (! -e $binaries->{'altloc_filter'}) {
      $binaries->{'altloc_filter'} = "ERROR" ;
   } else {
      $binaries->{'altloc_filter'} = $binaries->{perl}." ".
         $binaries->{'altloc_filter'} ;
   }


   $binaries->{'inscode_check'} = "$rootdir/auxil/".
      "inscode_check/inscode_check.$mach" ;
   if (! -e $binaries->{'inscode_check'}) {
      $binaries->{'inscode_check'} = "ERROR" ;
   }


   $binaries->{'ccp4sc'} = "$rootdir/auxil/ccp4sc/ccp4sc.$mach" ;
   if (! -e $binaries->{'ccp4sc'}) {
      $binaries->{'ccp4sc'} = "ERROR" ;
   } else {
      $binaries->{'ccp4sc'} .= " -nosummary -nohtml -n SCRADII $rootdir/auxil/ccp4sc/sc_radii.lib XYZIN" ;
   }

   $binaries->{'princip'} = "$rootdir/auxil/princip/princip.$mach" ;
   if (! -e $binaries->{'princip'}) {
      $binaries->{'princip'} = "ERROR" ;
   }

   $binaries->{'kdcontacts'} = "$rootdir/auxil/".
      "kdcontacts/kdcontacts.$mach" ;

   if (! -e $binaries->{'kdcontacts'}) {
      $binaries->{'kdcontacts'} = "ERROR" ;
   }

   $binaries->{'dssp'} = "$rootdir/auxil/".
      "dssp/dsspcmbi.$mach" ;

   if (! -e $binaries->{'dssp'}) {
      $binaries->{'dssp'} = "ERROR" ;
   }

   return $binaries ;
}



=head2 safe_move()

   Function:    Safely move a file to a directory (using File::Copy::move),
                  retries 14 times, and prints an error if it didnt work
   Return:	nothing
   Args:        $_[0] = source filename
                $_[1] = target directory

=cut

sub safe_move {

   my $file = shift ;
   my $dir = shift ;
   my $tries = 15 ;

   my $res = 0 ;
   if (!-s $dir) { mkpath($dir) ; }

   while (($tries > 0) && ($res == 0 )) {
      $res = File::Copy::move($file, $dir) ;
      $tries-- ;
   }

   if (-s $file) {
      print STDERR "ERROR: couldnt move $file to $dir\t$!\n";
   }

   return ;

}


=head2 safe_copy()

   Function:    Safely copy a file to a directory (using File::Copy::copy),
                  retries 14 times, and prints an error if it didnt work
   Return:	nothing
   Args:        $_[0] = source filename
                $_[1] = target directory

=cut

sub safe_copy {

   my $file = shift ;
   my $dest = shift ;
   my $tries = 15 ;

   if (!-s $file) {
      print STDERR "ERROR: couldnt find $file\n"; return}

   my $res = 0 ;
   while (($tries > 0) && ($res == 0 )) {
      $res = File::Copy::copy($file, $dest) ;
      $tries-- ;
   }

   if (!-s $dest) {
      print STDERR "ERROR: couldnt copy $file to $dest\t$!\n"; }

   return ;

}

1 ;
