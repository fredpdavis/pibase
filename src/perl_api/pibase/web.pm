=head1 NAME

pibase::web - collection of routines for PIBASE web interface

=head1 DESCRIPTION

This module provides routines for the PIBASE web interface.

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

package pibase::web ;
use strict;
use warnings;
use Exporter;
use Carp ;

use pibase qw/locate_binaries/;
our @ISA = qw/Exporter/ ;
our @EXPORT = qw/greeting closing find_object get_object_details pngconvert_bdp_interaction_topology_graph view_object/ ;


=head2 greeting()

   Title:       greeting()
   Function:    Provides pibase web site greeting
   Args:        $_->{base} [optional] = website base url
   Return:      html code for greeting section of pibase webpage

=cut

sub greeting {

   my $in = shift;
   my $base ;
   if (exists $in->{base}) {
      $base = $in->{base} ;
   } else {
      my $pb_specs = pibase::get_specs() ;
      $base = $pb_specs->{web}->{base_url}
   }

my $greeting = '<DIV class=bodytext>
<table width=65% align="center" border="0" width=100%>
<tr>
<td valign="center" align="right">
<IMG height=50 width=50 SRC="'.$base.'/figs/pibase_blue_web.png" BORDER=0>
<DIV class=large><B>PIBASE v2010</font></B></DIV>
</td>
<td align="left">
<a href="'.$base.'/queries.html">Queries</a> <a href="'.$base.'/queries_detailed.html">(detailed)</a> |
<a href="'.$base.'/introduction.html">Introduction</a> |
<a href="'.$base.'/statistics.html">Statistics</a> |
<a href="'.$base.'/download.html">Download</a> |
<a href="http://research.janelia.org/davis/projects.html">related work</a>
</td>
</tr>
</table>
</DIV>
<br>
' ;

   return $greeting ;

}


=head2 closing()

   Title:       closing()
   Function:    Provides pibase web site closing
   Args:        $_->{base} [optional] = website base url
   Return:      html code for closing section of pibase webpage

=cut


sub closing {

   my $in = shift;
   my $base ;
   if (exists $in->{base}) {
      $base = $in->{base} ;
   } else {
      my $pb_specs = pibase::get_specs() ;
      $base = $pb_specs->{www}->{base_url}
   }


   my $closing ;
   $closing  = '<p>
<DIV class=bodytext><a href=\"http://research.janelia.org/davis\">davisf</a> at <a href =\"http://hhmi.org/janelia\">janelia dot hhmi</a> dot org</DIV>
' ;

   return $closing ;

}


=head2 find_object()

   Title:       find_object()
   Function:    Routine for finding complex, interface, or domain
   Args:        $_->{base} [optional] = website base url
   Return:      STDOUT html page of search results

=cut

sub find_object {

   my $in = shift ;
   my $q = new CGI ;

   my $base ;
   if (exists $in->{base}) {
      $base = $in->{base} ;
   } else {
      my $pb_specs = pibase::get_specs() ;
      $base = $pb_specs->{web}->{base_url}
   }

   my $basecgi ;
   if (exists $in->{basecgi}) {
      $basecgi = $in->{basecgi} ;
   } else {
      my $pb_specs = pibase::get_specs() ;
      $basecgi = $pb_specs->{web}->{basecgi_url}
   }

   # get all options into a hash

   my $input = $q->Vars ;
   my $object_type = $input->{object_type} ;
   print $q->header('text/html') ;
   print $q->start_html(-title => 'PIBASE Query results',
             -style=>{'src'=>'http://salilab.org/~fred/style/styles.css'}) ;

   my $greeting = pibase::web::greeting() ;
   print $greeting ;
   print "<DIV class=bodytext>\n" ;
   print "<table width=95% align=\"center\"><td>\n" ;

   my $objecttype2table = {
      domains => "subsets",
      interfaces => "intersubset_contacts",
      complexes => "bdp_files",
   } ;

   my $results_fields = {
      domains => [
         'subsets.subset_id',
         'subsets.bdp_id',
         'bdp_files.pdb_id',
      ],
      interfaces => [
         'intersubset_contacts.subset_id_1',
         'intersubset_contacts.class_1',
         'intersubset_contacts.subset_id_2',
         'intersubset_contacts.class_2',
         'intersubset_contacts.num_contacts',
         'intersubset_contacts.chains',
         'intersubset_contacts.bdp_id',
         'bdp_files.pdb_id'
      ],
      complexes => [
         'bdp_files.bdp_id',
         'bdp_files.pdb_id',
         'pdb_entries.compound',
      ]
   } ;

   my $details_link =  {
      complexes => "object_type=complexes&bdp_id=F0",
      interfaces => "object_type=interfaces&subset_id_1=F0&subset_id_2=F2",
      domains => "object_type=domains&subset_id=F0",
   } ;

   my $view_link =  {
      interfaces => "subset_id_1=F0&subset_id_2=F2",
      domains => "subset_id=F0",
   } ;


   # begin the SQL query.

   my $master_query ;
   push @{$master_query},
      {table => $objecttype2table->{$object_type},
       fields => join(", ", @{$results_fields->{$object_type}})} ;
#   print "fields is ".$master_query->[0]->{fields}."\n" ;

   my $where_connector = {
      domains => {
         pdb_entries => "pdb_entries.pdb_id = bdp_files.pdb_id",
         bdp_files => "(bdp_files.bdp_id = subsets.bdp_id)",
         interface_sasa => "(interface_sasa.subset_id_1 = subsets.subset_id ".
            " OR interface_sasa.subset_id_2 = subsets.subset_id)",
         intersubset_contacts =>
            "(intersubset_contacts.subset_id_1 = subsets.subset_id ".
            " OR intersubset_contacts.subset_id_2 = subsets.subset_id)",
         bdp_interaction_topology =>
"(bdp_interaction_topology.bdp_id = subsets.bdp_id AND ".
" bdp_interaction_topology.subste_source_id = subsets.subset_source_id)",
      },
      interfaces => {
         subsets => "subsets.subset_id = intersubset_contacts.subset_id_1",
         bdp_files => "bdp_files.bdp_id = intersubset_contacts.bdp_id",
         pdb_entries => "pdb_entries.pdb_id = bdp_files.pdb_id",
         interface_sasa =>
"(interface_sasa.subset_id_1 = intersubset_contacts.subset_id_1 AND".
" interface_sasa.subset_id_2 = "."intersubset_contacts.subset_id_2)",
         interface_size =>
"(interface_size.subset_id_1 = intersubset_contacts.subset_id_1 AND".
" interface_size.subset_id_2 = "."intersubset_contacts.subset_id_2)",
         bdp_interaction_topology =>
"bdp_interaction_topology.bdp_id = intersubset_contacts.bdp_id"
      },
      complexes => {
         subsets => "subsets.bdp_id = bdp_files.bdp_id",
         pdb_entries => "pdb_entries.pdb_id = bdp_files.pdb_id",
         interface_sasa => "interface_sasa.bdp_id = bdp_files.bdp_id",
         intersubset_contacts=>"intersubset_contacts.bdp_id = bdp_files.bdp_id",
         bdp_interaction_topology =>
"bdp_interaction_topology.bdp_id = bdp_files.bdp_id"
      }
   } ;

   my $criteria = {
      complex_bdp_id => {
         table => "bdp_files",
         where => "bdp_files.bdp_id = \"VAL\""
      },
      complex_pdb_id => {
         table => "bdp_files",
         where => "bdp_files.pdb_id = \"VAL\""
      },
      complex_resolution => {
         table => "pdb_entries",
         where => "pdb_entries.resolution <= VAL",
         requires => ["bdp_files"]
      },
      complex_domain_type => {
         table => "bdp_interaction_topology",
         where => "bdp_interaction_topology.subset_source_id = \"VAL\""
      },
      complex_pdb_keyword=> {
         table => "pdb_entries",
         where => "MATCH(pdb_entries.header,pdb_entries.compound) ".
                  "AGAINST(\'VAL\' IN BOOLEAN MODE)",
         requires => ["bdp_files"]
      },
      complex_topology => {
         table => "bdp_interaction_topology",
         where => "(bdp_interaction_topology.nodestring = \"VAL\" OR ".
                  "bdp_interaction_topology.edgestring = \"VAL\")"
      },
      complex_num_domains => {
         table => "bdp_interaction_topology",
         where => "bdp_interaction_topology.num_domains OPERATOR VAL"
      },
      complex_num_interfaces => {
         table => "bdp_interaction_topology",
         where => "bdp_interaction_topology.num_edges OPERATOR VAL"
      },
      interface_subset_id_1 => {
         table => "intersubset_contacts",
         where => "(intersubset_contacts.subset_id_1 = \"VAL\" OR ".
                  "intersubset_contacts.subset_id_2 = \"VAL\")"
      },
      interface_subset_id_2 => {
         table => "intersubset_contacts",
         where => "(intersubset_contacts.subset_id_1 = \"VAL\" OR ".
                  "intersubset_contacts.subset_id_2 = \"VAL\")"
      },
      interface_class_1 => {
         table => "intersubset_contacts",
         where => "(intersubset_contacts.class_1 LIKE \"VAL\" OR ".
                  "intersubset_contacts.class_2 LIKE \"VAL\")"
      },
      interface_class_2 => {
         table => "intersubset_contacts",
         where => "(intersubset_contacts.class_1 LIKE \"VAL\" OR ".
                  "intersubset_contacts.class_2 LIKE \"VAL\")"
      },
      interface_chains => {
         table => "intersubset_contacts",
         where => "chains = \"VAL\""
      },
      interface_domain_type => {
         table => "subsets",
         where => "subsets.subset_source_id = \"VAL\"",
         requires => ["subsets"] 
      },
      interface_dsasa => {
         table => "interface_sasa",
         where => "interface_sasa.dsasa_all OPERATOR VAL"
      },
      interface_perc_dsasa_polar => {
         table => "interface_sasa",
         where => "(100 * interface_sasa.dsasa_polar / interface_sasa.dsasa_all) OPERATOR VAL"
      },
      domain_subset_id => {
         table => "subsets",
         where => "subsets.subset_id = \"VAL\""
      },
      domain_type => {
         table => "subsets",
         where => "subsets.subset_source_id = \"VAL\""
      },
      domain_class => {
         table => "subsets",
         where => "subsets.class LIKE \"VAL\""
      },
   } ;

   my ($tables, $where_clauses) ;
   my @table_list ; my $table_connected ;

# any default connections - eg interfaces: 
   if ($object_type eq 'complexes') {
      push @table_list, 'pdb_entries' ;
      $tables->{'pdb_entries'}++ ;
      push @{$where_clauses},
         $where_connector->{$object_type}->{'pdb_entries'} ;
      $table_connected->{'pdb_entries'}++ ;
   }

   if ($object_type eq 'subsets' || $object_type eq 'interfaces') {
      push @table_list, 'bdp_files';
      $tables->{'bdp_files'}++ ;
      push @{$where_clauses},
         $where_connector->{$object_type}->{'bdp_files'} ;
      $table_connected->{'bdp_files'}++ ;
   }

   foreach my $field (keys %{$input}) {
      if (!defined $input->{$field} || $input->{$field} eq '') {next;}
      my $val = $input->{$field} ;
      if (exists $criteria->{$field}) {
         my $cur_table = $criteria->{$field}->{table} ;
         my $cur_field = $criteria->{$field}->{where} ;

         $cur_field =~ s/VAL/$val/g ;
         if ($cur_field =~ /OPERATOR/) {
            my $operator = $input->{$field.'_operator'} ;
            $cur_field =~ s/OPERATOR/$operator/g ;
         }
         push @{$where_clauses}, $cur_field ;
         if ($cur_table ne $master_query->[0]->{table}) {
            if (!exists $tables->{$cur_table}) {
               push @table_list, $cur_table ; }
            $tables->{$cur_table}++ ;
         }
         if (!exists $table_connected->{$cur_table} &&
             exists $where_connector->{$object_type}->{$cur_table}) {
            push @{$where_clauses},
               $where_connector->{$object_type}->{$cur_table} ;
            $table_connected->{$cur_table}++ ;
         }

         if (exists $criteria->{$field}->{requires}) {
            foreach my $req_table (@{$criteria->{$field}->{requires}}) {
               if (!exists $table_connected->{$req_table} &&
                   exists $where_connector->{$object_type}->{$req_table}) {
                  push @{$where_clauses},
                     $where_connector->{$object_type}->{$req_table} ;
                  if (!exists $tables->{$req_table}) {
                     push @table_list, $req_table ; }
                  $tables->{$req_table}++ ;
                  $table_connected->{$req_table}++ ;
               }
            }
         }

      }
#      else {
#         print STDERR "$field not recognized\n" ;
#      }
   }

   # based on what type of object, call list_interfaces, list_domains, or
   #  list_complexes

   my $select_clause = "SELECT distinct ".$master_query->[0]->{fields};

   my $from_clause = " FROM ".$master_query->[0]->{table}.', ' ;
   $from_clause .= join(", ", @table_list) ;
   $from_clause =~ s/, $// ;

   my $where_clause = "WHERE ".join(" AND ", @{$where_clauses}) ;

   my $full_query = $select_clause." ".$from_clause." ".$where_clause." ".
      "ORDER BY bdp_files.bdp_id" ;

   if (!exists $input->{view_query_only}) {
      my ($dbh) = pibase::connect_pibase ;
      my @results = pibase::mysql_fetchcols($dbh, $full_query) ;
      display_results({
         object_type => $object_type,
         results_fields => $results_fields,
         details_link => $details_link,
         view_link => $view_link,
         cgi_h => $q,
         input => $input,
         data => \@results,
         dbh => $dbh
      }) ;
   }
   print "<p>query: $full_query\n" ;
   print "</td></table>\n" ;

   my $closing = pibase::web::closing() ;
   print $closing ;
   $q->end_html ;
    
}


=head2 display_results()

   Title:       display_results()
   Function:    Formats a list of query results (find_object()) for html viewing
   Args:        $_->{input}->{object_type} = complexes|interfaces|domains
                $_->{details_link}->{object_type} = construct for detail URL link
                $_->{view_link}->{object_type} = construct for viewer URL link
                $_->{results_field}->[i] = ith result field header
                $_->{data}->[i]->[j] = jth field of ith record
                $_->{cgi_h} = CGI handle
                $_->{base} [optional] = website base url
                $_->{basecgi} [optional] = website base url
   Returns:     STDOUT html page of search results

=cut

sub display_results {

   my $in = shift ;
   my $input = $in->{input};
   my $data = $in->{data} ;
   my $q = $in->{cgi_h} ;

   my $base ;
   if (exists $in->{base}) {
      $base = $in->{base} ;
   } else {
      my $pb_specs = pibase::get_specs() ;
      $base = $pb_specs->{web}->{base_url}
   }
   my $basecgi ;
   if (exists $in->{basecgi}) {
      $basecgi = $in->{basecgi} ;
   } else {
      my $pb_specs = pibase::get_specs() ;
      $basecgi = $pb_specs->{web}->{basecgi_url}
   }

   my $object_type = $input->{object_type} ;
   print "<B>" ;
   if ($#{$data->[0]} < 0) {
      print CGI::div({-class=>"large"},
         $q->p("No matching $object_type found"));
   } else {
      my $numhits = $#{$data->[0]} + 1 ;
      print CGI::div({-class=>"large"},
         $q->p("Found $numhits matching $object_type"));
   }
   print "</B>" ;
 
   my $numfields = $#{$data} ;
   print "<TABLE bordercolor=black frame=\"hsides\">\n" ;

# print headers
   my @header_fields= ('','',@{$in->{results_fields}->{$object_type}});
   foreach my $j ( 0 .. $#header_fields) { $header_fields[$j] =~ s/.*\.// ; }
   print $q->Tr( $q->th( \@header_fields)) ;

# figure out if even or odd numbered results entry
   foreach my $j (0 .. $#{$data->[0]}) {
      my $prediv; my $postdiv ;
      if (($j % 2) == 0) {
         $prediv = "<DIV class=treven>" ;
      } else {
         $prediv = "<DIV class=trodd>" ;
      }
      $postdiv = "</DIV>\n" ;

# print results.
      my @outvals = (($j + 1));
      my $details_link = $in->{details_link}->{$object_type} ;
      $details_link =~ s/F([0-9]+)/$data->[$1]->[$j]/g ;
      $details_link = "<a href=$basecgi/get_details.pl?".
         "$details_link>details</a>" ;
      my $view_link = '';
      if (exists $in->{view_link}->{$object_type}) {
         $view_link = $in->{view_link}->{$object_type} ;
         $view_link =~ s/F([0-9]+)/$data->[$1]->[$j]/g ;
#         $view_link = "<a href=$basecgi/view.pl?$view_link>view</a>" ;
         $view_link = "<a href=$basecgi/view.pl/out.raspdb?".
            "$view_link>view</a>" ;
         $view_link = ' | '.$view_link ;
      }
      push @outvals, $details_link.$view_link ;

      foreach my $k ( 0 .. $#{$data}) {
         push @outvals, $prediv.$data->[$k]->[$j].$postdiv ; }

# make a details link - according to specified specs in $in, push onto outvals

      print $q->Tr( $q->td( \@outvals)) ;
   }
   print "</TABLE>\n" ;

}


=head2 get_object_details()

   Title:       get_object_details()
   Function:    Formats a list of query results (find_object()) for html viewing
   Args:        $_->{input}->{object_type} = complexes|interfaces|domains
                $_->{details_link}->{object_type} = construct for detail URL link
                $_->{view_link}->{object_type} = construct for viewer URL link
                $_->{results_field}->[i] = ith result field header
                $_->{data}->[i]->[j] = jth field of ith record
                $_->{cgi_h} = CGI handle
                $_->{base} [optional] = website base url
                $_->{basecgi} [optional] = website base url
   Returns:     STDOUT html page of search results

=cut

sub get_object_details { #routine to get details of domains/interfaces/complexes

   my $in = shift ;
   my $q = new CGI ;

   my $pb_specs = pibase::get_specs() ;
   my $base ;
   if (exists $in->{base}) {
      $base = $in->{base} ;
   } else {
      $base = $pb_specs->{web}->{base_url} ;
   }

   my $basecgi ;
   if (exists $in->{basecgi}) {
      $basecgi = $in->{basecgi} ;
   } else {
      $basecgi = $pb_specs->{web}->{basecgi_url}
   }


   # get all options into a hash

   my $input = $q->Vars ;
   my $object_type = $input->{object_type} ;
   print $q->header('text/html') ;
   print $q->start_html(-title => 'PIBASE '.$object_type.' details',
             -style=>{'src'=>'http://salilab.org/~fred/style/styles.css'}) ;

   my $greeting = pibase::web::greeting() ;
   print $greeting ;
   print "<DIV class=bodytext>\n" ;
   print "<table width=95% align=\"center\"><td>\n" ;

   my ($dbh) = pibase::connect_pibase ;

   if ($object_type eq 'complexes') {
      my $bdp_id = $input->{bdp_id} ;

      print $q->h3("Details of complex ".$bdp_id) ;

      my ($file_base, $pdb_id, $raw_pdb) ;
      {
         my ($t_file_base, $t_pdb_id, $t_raw_pdb)=pibase::mysql_fetchcols($dbh,
            "SELECT file_base, pdb_id, raw_pdb FROM bdp_files ".
            "where bdp_id = ".$bdp_id) ;
         $file_base = $t_file_base->[0] ;
         $pdb_id = $t_pdb_id->[0] ;
         $raw_pdb = $t_raw_pdb->[0] ;
      }
      my $uc_pdb_id = uc($pdb_id) ;

      my $compound = pibase::mysql_singleval($dbh,
         "SELECT compound FROM pdb_entries where pdb_id = \"$pdb_id\"") ;
      print "<b>Compound:</b> $compound<p>\n" ;


      if ($raw_pdb == 0) {
         print "<b>Source:</b> This complex is an entry in ".
"<a href=\"http://www.ebi.ac.uk/msd-srv/pisa/\">PISA</a> (file ".
"<a href =\"http://www.ebi.ac.uk/msd-srv/pisa/cgi-bin/piserver?qa=$pdb_id\">".
$file_base."</a>) that is derived from <a href=\"http://wwpdb.org\">PDB</a> ".
"entry <a href =\"http://www.pdb.org/pdb/explore/explore.do?structureId=".
$uc_pdb_id."\">$uc_pdb_id</a><br>\n" ;
      } else {
         print "<b>Source:</b> This complex is ".
"<a href=\"http://wwpdb.org\">PDB</a> entry ".
"<a href =\"http://www.pdb.org/pdb/explore/explore.do?structureId=".
$uc_pdb_id."\">$uc_pdb_id</a><br>\n" ;
      }

      my ($subsetsource_2_name) = pibase::mysql_hashload($dbh,
         "SELECT subset_source_id, concat(subset_source, \" ver \",".
         "version) FROM subsets_source") ;

      my ($subset_id, $class, $class_desc, $t_subset_source_id) =
         pibase::mysql_fetchcols($dbh,
         "SELECT a.subset_id, a.class, b.description, a.subset_source_id FROM ".
"subsets as a LEFT JOIN subsets_class as b ON (a.class = b.class) ".
"WHERE a.bdp_id = ".$bdp_id);
         
      my $domain_types ;
      foreach my $j ( 0 .. $#{$subset_id}) {
         push @{$domain_types->{$t_subset_source_id->[$j]}},$j ; }

      my ($subset_id_1, $class_1, $subset_id_2, $class_2, $num_contacts,
         $chains, $t_int_subset_source_id) = pibase::mysql_fetchcols($dbh,
"SELECT subset_id_1, class_1, subset_id_2, class_2, num_contacts, ".
"chains, subset_source_id FROM intersubset_contacts as a, subsets as b ".
"WHERE a.bdp_id = $bdp_id AND a.subset_id_1 = b.subset_id") ;

      my $interface_types ;
      foreach my $j ( 0 .. $#{$subset_id_1}) {
         push @{$interface_types->{$t_int_subset_source_id->[$j]}},$j ; }

      my @domainclass_string ;
      foreach my $subset_source_id (keys %{$domain_types}) {
         push @domainclass_string,
            "<a href =\"#subsetsource_$subset_source_id\">".
            $subsetsource_2_name->{$subset_source_id}."</a>" ;
      }
      if ($#domainclass_string >= 0) {
         print "<p>This complex has been parsed using domain definitions from ".
         join(", ", @domainclass_string)."<br>\n";
      } else {
         print "<p>No domain definitions are available for this complex<br>\n";
      }

      foreach my $subset_source_id (keys %{$domain_types}) {
      my $domsource_link ;
      if ($subsetsource_2_name->{$subset_source_id} =~ /scop/) {
         $domsource_link = $pb_specs->{external_data}->{scop}->{base_url} ;
      } elsif ($subsetsource_2_name->{$subset_source_id} =~ /cath/) {
         $domsource_link = $pb_specs->{external_data}->{cath}->{base_url} ;
      } elsif ($subsetsource_2_name->{$subset_source_id} =~ /chain/) {
         $domsource_link = $pb_specs->{external_data}->{pdb}->{base_url} ;
      }

         print $q->h3("Interactions based on ".
            "<a href=\"$domsource_link\">".
            $subsetsource_2_name->{$subset_source_id}."</a>") ;
         print "<a name=\"subsetsource_".$subset_source_id."\"</a>" ;

         my $view_link = $basecgi."/view.pl/out.raspdb?bdp_id=$bdp_id&".
            "subset_source_id=$subset_source_id" ;
         print "<a href =\"$view_link\">View complex</a><br>\n" ;

         my ($topo_graph_eps_fn)= pibase::mysql_singleval($dbh,
            "SELECT file_path FROM bdp_interaction_topology_graph ".
            "WHERE bdp_id = $bdp_id and subset_source_id = $subset_source_id") ;
         my $topo_graph_url = $topo_graph_eps_fn ;
         $topo_graph_url=~ s/eps$/png/ ;
         my $t_1 = $pb_specs->{otherdata_dir}->{bdp_topology_graphs} ;
         my $t_2 = $pb_specs->{web}->{bdp_topology_graphs_baseurl} ;
         $topo_graph_url =~ s/$t_1/$t_2/ ;

         print "<IMG src=\"$topo_graph_url\"/>\n" ;
         print "<TABLE bordercolor=black frame=\"hsides\">\n" ;
         print $q->Tr($q->th(["<i>Domains</i>","","subset_id",
            "class", "description"])) ;

# get domain colors;
         my $class2bgcol ; my $class2fgcol ;
         my $class_counter = 0;
         foreach my $t_class (sort
            @$class[@{$domain_types->{$subset_source_id}}]) {
            if (exists $class2bgcol->{$t_class}) {next;}
            my $t_col = get_color_codes({
               name => $pb_specs->{web}->{domain_colors}->[$class_counter]}) ;
            $class2bgcol->{$t_class}= $t_col->{hexcode} ;
            my @rgb = split(' ' , $t_col->{rgb}) ;
            if (($rgb[0] + $rgb[1] + $rgb[2]) < 1.2) {
               $class2fgcol->{$t_class} = "FFFFFF" ;
            } else {
               $class2fgcol->{$t_class} = "000000" ;
            }
            $class_counter++ ;
         }

         my $row_counter = 1 ;
         foreach my $entry_no (@{$domain_types->{$subset_source_id}}) {

            my $d_details_link = $basecgi."/get_details.pl?object_type=domains".
               "&subset_id=$subset_id->[$entry_no]" ;
            my $d_view_link=$basecgi."/view.pl/out.raspdb?".
               "subset_id=$subset_id->[$entry_no]" ;

            my $prediv; my $postdiv ;
            if (($row_counter % 2) == 0) {
               $prediv = "<DIV class=treven>" ;
            } else {
               $prediv = "<DIV class=trodd>" ;
            }
            $postdiv = "</DIV>\n" ;

            my @data = ($row_counter, "<a href=\"$d_details_link\">details</a>".
               " | <a href=\"$d_view_link\">view</a>", $subset_id->[$entry_no],
               $class->[$entry_no]) ;

            if (!defined $class_desc->[$entry_no]) {$class_desc->[$entry_no] = '';}
            my $col_desc = "<td bgcolor=#".$class2bgcol->{$class->[$entry_no]}.
               "><font color=#".$class2fgcol->{$class->[$entry_no]}.
               ">".$class_desc->[$entry_no]."</font></td>" ;

            my @outvals ;
            foreach my $k ( 0 .. $#data) {
               if (!defined $data[$k]) {$data[$k] = '';}
               push @outvals, $prediv.$data[$k].$postdiv ;
            }
            print $q->Tr($q->td(\@outvals),$col_desc) ;
            $row_counter++ ;
         }
         print "</TABLE>\n" ;

         print "<TABLE bordercolor=black frame=\"hsides\">\n" ;
         print $q->Tr($q->th(["<i>Interfaces</i>",'',"subset_id_1", "class_1",
                  "subset_id_2", "class_2", "num_contacts", "chains"])) ;
         $row_counter = 1 ;
         foreach my $entry_no (@{$interface_types->{$subset_source_id}}) {
            my $i_details_link = $basecgi.
               "/get_details.pl?object_type=interfaces".
               "&subset_id_1=$subset_id_1->[$entry_no]".
               "&subset_id_2=$subset_id_2->[$entry_no]" ;

            my $i_view_link=$basecgi."/view.pl/out.raspdb?".
               "subset_id_1=$subset_id_1->[$entry_no]".
               "&subset_id_2=$subset_id_2->[$entry_no]" ;

            my $prediv; my $postdiv ;
            if (($row_counter % 2) == 0) {
               $prediv = "<DIV class=treven>" ;
            } else {
               $prediv = "<DIV class=trodd>" ;
            }
            $postdiv = "</DIV>\n" ;

            my @data = ($row_counter,"<a href=\"$i_details_link\">details</a> ".
                        " | <a href=\"$i_view_link\">view</a>",
                        $subset_id_1->[$entry_no], $class_1->[$entry_no],
                        $subset_id_2->[$entry_no], $class_2->[$entry_no],
                        $num_contacts->[$entry_no],
                        $chains->[$entry_no]) ;
            my @outvals ;
            foreach my $k ( 0 .. $#data) {
               push @outvals, $prediv.$data[$k].$postdiv ; }
            print $q->Tr($q->td(\@outvals)) ;
            $row_counter++ ;
         }
         print "</TABLE>\n" ;

      }
      
   } elsif ($object_type eq 'interfaces') {


      my $subset_id_1= $input->{subset_id_1} ;
      my $subset_id_2= $input->{subset_id_2} ;

      print $q->h3("Details of interface $subset_id_1 -- $subset_id_2");

      my ($bdp_id)=pibase::mysql_singleval($dbh,
         "SELECT bdp_id FROM intersubset_contacts WHERE ".
         "subset_id_1 = \"$subset_id_1\" AND subset_id_2 = \"$subset_id_2\"") ;

      print "<b>Source:</b> This interface is from PIBASE complex ".
"<a href =\"$basecgi/get_details.pl?object_type=complexes&bdp_id=$bdp_id\">".
$bdp_id."</a><br>\n" ;

      my ($class, $class_desc) ;
      foreach my $subset_id ($subset_id_1, $subset_id_2) {

         ($class->{$subset_id}, $class_desc->{$subset_id}) =
            pibase::mysql_fetchcols($dbh,
"SELECT a.class, b.description, a.subset_source_id FROM ".
"subsets as a LEFT JOIN subsets_class as b ON (a.class = b.class) ".
"WHERE a.subset_id = \"$subset_id\"");
      } 
      print "<TABLE bordercolor=black frame=\"hsides\">\n" ;
      print $q->Tr($q->th(["<i>Domains</i>","","subset_id",
         "class", "description"])) ;
      my $j = 1 ;

      foreach my $subset_id ($subset_id_1, $subset_id_2) {

         my $details_link = $basecgi."/get_details.pl?object_type=domains".
            "&subset_id=$subset_id" ;

         my $prediv; my $postdiv ;
         if (($j % 2) == 0) {
            $prediv = "<DIV class=treven>" ;
         } else {
            $prediv = "<DIV class=trodd>" ;
         }
         $postdiv = "</DIV>\n" ;

         if (!defined $class_desc->{$subset_id}->[0]) {
            $class_desc->{$subset_id}->[0] = '';}
         my @data = ($j, "<a href=\"$details_link\">details</a>",
                     $subset_id, $class->{$subset_id}->[0],
                     $class_desc->{$subset_id}->[0]) ;
         my @outvals ;
         foreach my $k ( 0 .. $#data) {
            if (!defined $data[$k]) {$data[$k] = '';}
            push @outvals, $prediv.$data[$k].$postdiv ;
         }
         print $q->Tr($q->td(\@outvals)) ;
         $j++ ;
      }
      print "</TABLE>\n" ;

      {
      print $q->h4("Contacts across the interface:") ;
      my ($num_contacts, $cutoff, $num_contacts_4, $num_contacts_4p5,
          $num_contacts_5, $num_contacts_5p5, $hbond, $salt, $ssbond, $chains)=
         pibase::mysql_fetchcols($dbh,
"SELECT num_contacts, cutoff, num_contacts_4, num_contacts_4p5,".
"num_contacts_5, num_contacts_5p5, hbond, salt, ssbond, chains FROM ".
"intersubset_contacts WHERE subset_id_1 = \"$subset_id_1\" AND ".
"subset_id_2 = \"$subset_id_2\"") ;

      if ($chains->[0] eq 'both') {
         print "This interface is both intra- and inter-chain<br>\n" ;
      } elsif ($chains->[0] eq 'same') {
         print "This is an intra-chain interface<br>\n" ;
      } elsif ($chains->[0] eq 'diff') {
         print "This is an inter-chain interface<br>\n" ;
      }
      print $num_contacts->[0]." inter-atomic contacts are made across the interface ".
"within ".$cutoff->[0]." &#8491;<br>\n".
" ( ".$num_contacts_4->[0]." inter-atomic contacts within 4 &#8491; ".
", ".$num_contacts_4p5->[0]." within 4.5 &#8491; ".
", ".$num_contacts_5->[0]." within 5 &#8491; ".
", ".$num_contacts_5p5->[0]." within 5.5 &#8491;)<br>\n";


      print $hbond->[0]." possible hydrogen-bonds (".$salt->[0]." salt bridges) and ".$ssbond->[0]." possible disulfide bridges<br>\n" ;
      }

      {
      print $q->h4("Buried solvent accessible surface area:") ;
      my ($dsasa_all, $dsasa_sc, $dsasa_mc, $dsasa_polar, $dsasa_nonpolar)=
         pibase::mysql_fetchcols($dbh,
"SELECT dsasa_all, dsasa_sc, dsasa_mc, dsasa_polar, dsasa_nonpolar ".
"FROM interface_sasa WHERE subset_id_1 = \"$subset_id_1\" AND ".
"subset_id_2 = \"$subset_id_2\"") ;
      print $dsasa_all->[0]." &#8491;&#178; of solvent accessible surface area".
" is buried at the interface<br>\nof which ".
$dsasa_polar->[0]." &#8491;&#178; is polar, ".
$dsasa_nonpolar->[0]." &#8491;&#178; is nonpolar\n and ".
$dsasa_sc->[0]." &#8491;&#178; is from side chain".
" and ".$dsasa_mc->[0]." &#8491;&#178; is from main chain atoms\n" ;
      }

      {
      print $q->h4("Residues at the interface:") ;
      print "<b>NOTE:</b> The lists below ARE NOT pairwise contacts, but rather independent lists of residues<br>\n";
      my ($num_res_1, $num_res_2, $subset_size_1, $subset_size_2,
          $intres_1, $intres_2) = pibase::mysql_fetchcols($dbh,
"SELECT num_res_1, num_res_2, subset_size_1, subset_size_2, ".
"interface_res_1, interface_res_2 ".
"FROM interface_size WHERE subset_id_1 = \"$subset_id_1\" AND ".
"subset_id_2 = \"$subset_id_2\"") ;

      print "<TABLE bordercolor=black frame=\"hsides\">\n" ;
      print $q->Tr($q->th(["<i>feature</i>","Domain 1","Domain 2"])) ;
      my $resstring1= join("<br>", sort {$a <=> $b} split(',', $intres_1->[0]));
      my $resstring2= join("<br>", sort {$a <=> $b} split(',', $intres_2->[0]));
      print $q->Tr($q->td(["domain length", $subset_size_1->[0],
                                            $subset_size_2->[0]])) ;
      print $q->Tr($q->td(["number of interface residues",
                           $num_res_1->[0], $num_res_2->[0]]));
      print $q->Tr($q->td({valign => 'top'},
                   ["interface residues", $resstring1, $resstring2])) ;
      print "</TABLE>\n" ;

#      print $num_res_1->[0]." of the ".$subset_size_1->[0].
#         " residues in domain 1 participate in the interface:\n";
#      print $num_res_2->[0]." of the ".$subset_size_2->[0].
#         " residues in domain 2 participate in the interface:\n".
#         $intres_2->[0]."\n";
      }


   } elsif ($object_type eq 'domains') {

      my $subset_id = $input->{subset_id} ;

      print $q->h3("Details of domain $subset_id");

      my ($bdp_id)=pibase::mysql_singleval($dbh,
         "SELECT bdp_id FROM subsets WHERE subset_id = \"$subset_id\"") ;

      my ($class, $class_desc, $subset_source_id) =
            pibase::mysql_fetchcols($dbh,
"SELECT a.class, b.description, a.subset_source_id FROM ".
"subsets as a LEFT JOIN subsets_class as b ON (a.class = b.class) ".
"WHERE a.subset_id = \"$subset_id\"");

      my ($subsetsource_2_name) = pibase::mysql_hashload($dbh,
         "SELECT subset_source_id, concat(subset_source, \" ver \",".
         "version) FROM subsets_source") ;

      my $domsource_link ;
      if ($subsetsource_2_name->{$subset_source_id->[0]} =~ /scop/) {
         $domsource_link = $pb_specs->{external_data}->{scop}->{base_url} ;
      } elsif ($subsetsource_2_name->{$subset_source_id->[0]} =~ /cath/) {
         $domsource_link = $pb_specs->{external_data}->{cath}->{base_url} ;
      } elsif ($subsetsource_2_name->{$subset_source_id->[0]} =~ /chain/) {
         $domsource_link = $pb_specs->{external_data}->{pdb}->{base_url} ;
      }

      print "<b>Source:</b> This domain is from PIBASE complex ".
"<a href =\"$basecgi/get_details.pl?object_type=complexes&bdp_id=$bdp_id\">".
$bdp_id."</a>, the definition is based on <a href=\"$domsource_link\">".$subsetsource_2_name->{$subset_source_id->[0]}."</a><br>\n" ;



      print "<TABLE bordercolor=black frame=\"hsides\">\n" ;
      print $q->Tr($q->th(["class", "description"])) ;
      my $j = 1 ;

      my $prediv; my $postdiv ;
      if (($j % 2) == 0) {
         $prediv = "<DIV class=treven>" ;
      } else {
         $prediv = "<DIV class=trodd>" ;
      }
      $postdiv = "</DIV>\n" ;

      if (!defined $class_desc->[0] ) {
         $class_desc->[0] = '' ; }
      my @data = ($class->[0], $class_desc->[0]) ;
      my @outvals ;
      foreach my $k ( 0 .. $#data) {
            if (!defined $data[$k]) {$data[$k] = '';}
            push @outvals, $prediv.$data[$k].$postdiv ;
      }
      print $q->Tr($q->td(\@outvals)) ;
      print "</TABLE>\n" ;

      print "<TABLE bordercolor=black frame=\"hsides\">\n" ;
      print $q->Tr($q->th(["segment","chain_id",
         "start_resno", "end_resno"])) ;
      {
      print $q->h4("Domain definition:") ;
      my ($segment_id, $chain_id, $start_resno, $end_resno) =
         pibase::mysql_fetchcols($dbh,
"SELECT segment_id, chain_id, start_resno, end_resno ".
"FROM subsets_details WHERE subset_id = \"$subset_id\"") ;

      foreach my $j ( 0 .. $#{$segment_id}) {
         my $prediv; my $postdiv ;
         if (($j % 2) == 0) {
            $prediv = "<DIV class=treven>" ;
         } else {
            $prediv = "<DIV class=trodd>" ;
         }
         $postdiv = "</DIV>\n" ;

         my @data = ($segment_id->[$j], $chain_id->[$j],
            $start_resno->[$j], $end_resno->[$j]) ;
         my @outvals ;
         foreach my $k ( 0 .. $#data) {
            if (!defined $data[$k]) {$data[$k] = '';}
            push @outvals, $prediv.$data[$k].$postdiv ;
         }
         print $q->Tr($q->td(\@outvals)) ;
      }
      print "</TABLE>\n" ;


      my ($num_chains, $num_res, $sequence) = pibase::mysql_fetchcols($dbh,
"SELECT num_chains, num_res, sequence FROM subsets_sequence ".
"WHERE subset_id = \"$subset_id\"") ;

      print "Number of residues: ".$num_res->[0]."<br>\n" ;
      print "Number of chains: ".$num_chains->[0]."<br>\n" ;
      print $q->h4("Sequence:") ;
      print "<tt>".foldtext({string => $sequence->[0],
         foldchar=> "<br>\n", wrap=>60})."</tt><br>\n" ;
      }

      {
      print $q->h4("Solvent accessible surface area:") ;
      my ($sasa_all, $sasa_sc, $sasa_mc, $sasa_polar, $sasa_nonpolar)=
         pibase::mysql_fetchcols($dbh,
"SELECT sasa_all, sasa_sc, sasa_mc, sasa_polar, sasa_nonpolar ".
"FROM subsets_sasa WHERE subset_id = \"$subset_id\"") ;
      print "Total solvent accessible surface area: ".$sasa_all->[0].
         " &#8491;&#178;\n<br> of which ".
$sasa_polar->[0]." &#8491;&#178; is polar, ".
$sasa_nonpolar->[0]." &#8491;&#178; is nonpolar<br>\n and ".
$sasa_sc->[0]." &#8491;&#178; is from side chain".
" and ".$sasa_mc->[0]." &#8491;&#178; is from main chain atoms\n" ;
      }

   }
   print "</td></table>\n";
   

   # begin the SQL query.

#   if (exists $input->{run_query}) {
#      my ($dbh) = pibase::connect_pibase ;
#      my @results = pibase::mysql_fetchcols($dbh, $full_query) ;
#      display_results({
#         object_type => $object_type,
#         results_fields => $results_fields,
#         details_link => $details_link,
#         cgi_h => $q,
#         input => $input,
#         data => \@results,
#         dbh => $dbh
#      }) ;
#   }
#   print "<p>query: $full_query\n" ;

   my $closing = pibase::web::closing() ;
   print $closing ;
   $q->end_html ;
    
}


=head2 foldtext()

   Title:       foldtext()
   Function:    Wrapping code for a string with specified folding characters
   Args:        $_->{string} = 'ORIGINALSTRING'
                $_->{wrap} = wrapping length
                $_->{foldchar} = folding charcter [optional - defaults to newline]
   Returns:     $_ = "ORIG\nINAL\nSTRI\nNG" ;

=cut

sub foldtext {

   my $in = shift ;
   my $folded = '';
   my $j = 0 ;
   my $foldchar ;
   if (exists $in->{foldchar}) {
      $foldchar = $in->{foldchar} ;
   } else {
      $foldchar = "\n" ;
   }

   while ($j < length($in->{string})) {
      if (($j > 0) && ($j % $in->{wrap} == 0)) {
         $folded .= $foldchar ;
      }
      $folded .= substr($in->{string},$j,1) ;
      $j++ ;
   }

   return $folded ;
}


=head2 get_color_codes()

   Title:       get_color_codes()
   Function:    Returns rgb decimal and hex values for a specified color name
   Args:        $_->{name} = color name
   Returns:     $_->{rgb} = "R G B" (0-1 scale)
                $_->{rgb255} => "R,G,B" (0-255 scale)
                $_->{hexcode} => "RRGGBB" hex code

=cut

sub get_color_codes {

   my $in = shift ;
   my $color_name = $in->{name} ;

my $color2rgb = {
   black => "0.000 0.000 0.000",
   grey => "0.745 0.745 0.745",
   DimGrey => "0.412 0.412 0.412",
   LightGray => "0.827 0.827 0.827",
   LightSlateGrey => "0.467 0.533 0.600",
   SlateGray => "0.439 0.502 0.565",
   SlateGray1 => "0.776 0.886 1.000",
   SlateGray2 => "0.725 0.827 0.933",
   SlateGray3 => "0.624 0.714 0.804",
   SlateGray4 => "0.424 0.482 0.545",
   SlateGrey => "0.439 0.502 0.565",
   grey0 => "0.000 0.000 0.000",
   grey1 => "0.012 0.012 0.012",
   grey2 => "0.020 0.020 0.020",
   grey3 => "0.031 0.031 0.031",
   grey4 => "0.039 0.039 0.039",
   grey5 => "0.051 0.051 0.051",
   grey6 => "0.059 0.059 0.059",
   grey7 => "0.071 0.071 0.071",
   grey8 => "0.078 0.078 0.078",
   grey9 => "0.090 0.090 0.090",
   grey10 => "0.102 0.102 0.102",
   grey11 => "0.110 0.110 0.110",
   grey12 => "0.122 0.122 0.122",
   grey13 => "0.129 0.129 0.129",
   grey14 => "0.141 0.141 0.141",
   grey15 => "0.149 0.149 0.149",
   grey16 => "0.161 0.161 0.161",
   grey17 => "0.169 0.169 0.169",
   grey18 => "0.180 0.180 0.180",
   grey19 => "0.188 0.188 0.188",
   grey20 => "0.200 0.200 0.200",
   grey21 => "0.212 0.212 0.212",
   grey22 => "0.220 0.220 0.220",
   grey23 => "0.231 0.231 0.231",
   grey24 => "0.239 0.239 0.239",
   grey25 => "0.251 0.251 0.251",
   grey26 => "0.259 0.259 0.259",
   grey27 => "0.271 0.271 0.271",
   grey28 => "0.278 0.278 0.278",
   grey29 => "0.290 0.290 0.290",
   grey30 => "0.302 0.302 0.302",
   grey31 => "0.310 0.310 0.310",
   grey32 => "0.322 0.322 0.322",
   grey33 => "0.329 0.329 0.329",
   grey34 => "0.341 0.341 0.341",
   grey35 => "0.349 0.349 0.349",
   grey36 => "0.361 0.361 0.361",
   grey37 => "0.369 0.369 0.369",
   grey38 => "0.380 0.380 0.380",
   grey39 => "0.388 0.388 0.388",
   grey40 => "0.400 0.400 0.400",
   grey41 => "0.412 0.412 0.412",
   grey42 => "0.420 0.420 0.420",
   grey43 => "0.431 0.431 0.431",
   grey44 => "0.439 0.439 0.439",
   grey45 => "0.451 0.451 0.451",
   grey46 => "0.459 0.459 0.459",
   grey47 => "0.471 0.471 0.471",
   grey48 => "0.478 0.478 0.478",
   grey49 => "0.490 0.490 0.490",
   grey50 => "0.498 0.498 0.498",
   grey51 => "0.510 0.510 0.510",
   grey52 => "0.522 0.522 0.522",
   grey53 => "0.529 0.529 0.529",
   grey54 => "0.541 0.541 0.541",
   grey55 => "0.549 0.549 0.549",
   grey56 => "0.561 0.561 0.561",
   grey57 => "0.569 0.569 0.569",
   grey58 => "0.580 0.580 0.580",
   grey59 => "0.588 0.588 0.588",
   grey60 => "0.600 0.600 0.600",
   grey61 => "0.612 0.612 0.612",
   grey62 => "0.620 0.620 0.620",
   grey63 => "0.631 0.631 0.631",
   grey64 => "0.639 0.639 0.639",
   grey65 => "0.651 0.651 0.651",
   grey66 => "0.659 0.659 0.659",
   grey67 => "0.671 0.671 0.671",
   grey68 => "0.678 0.678 0.678",
   grey69 => "0.690 0.690 0.690",
   grey70 => "0.702 0.702 0.702",
   grey71 => "0.710 0.710 0.710",
   grey72 => "0.722 0.722 0.722",
   grey73 => "0.729 0.729 0.729",
   grey74 => "0.741 0.741 0.741",
   grey75 => "0.749 0.749 0.749",
   grey76 => "0.761 0.761 0.761",
   grey77 => "0.769 0.769 0.769",
   grey78 => "0.780 0.780 0.780",
   grey79 => "0.788 0.788 0.788",
   grey80 => "0.800 0.800 0.800",
   grey81 => "0.812 0.812 0.812",
   grey82 => "0.820 0.820 0.820",
   grey83 => "0.831 0.831 0.831",
   grey84 => "0.839 0.839 0.839",
   grey85 => "0.851 0.851 0.851",
   grey86 => "0.859 0.859 0.859",
   grey87 => "0.871 0.871 0.871",
   grey88 => "0.878 0.878 0.878",
   grey89 => "0.890 0.890 0.890",
   grey90 => "0.898 0.898 0.898",
   grey91 => "0.910 0.910 0.910",
   grey92 => "0.922 0.922 0.922",
   grey93 => "0.929 0.929 0.929",
   grey94 => "0.941 0.941 0.941",
   grey95 => "0.949 0.949 0.949",
   grey96 => "0.961 0.961 0.961",
   grey97 => "0.969 0.969 0.969",
   grey98 => "0.980 0.980 0.980",
   grey99 => "0.988 0.988 0.988",
   grey100 => "1.000 1.000 1.000",
   AliceBlue => "0.941 0.973 1.000",
   BlueViolet => "0.541 0.169 0.886",
   CadetBlue => "0.373 0.620 0.627",
   CadetBlue1 => "0.596 0.961 1.000",
   CadetBlue2 => "0.557 0.898 0.933",
   CadetBlue3 => "0.478 0.773 0.804",
   CadetBlue4 => "0.325 0.525 0.545",
   CornflowerBlue => "0.392 0.584 0.929",
   DarkSlateBlue => "0.282 0.239 0.545",
   DarkTurquoise => "0.000 0.808 0.820",
   DeepSkyBlue => "0.000 0.749 1.000",
   DeepSkyBlue1 => "0.000 0.749 1.000",
   DeepSkyBlue2 => "0.000 0.698 0.933",
   DeepSkyBlue3 => "0.000 0.604 0.804",
   DeepSkyBlue4 => "0.000 0.408 0.545",
   DodgerBlue => "0.118 0.565 1.000",
   DodgerBlue1 => "0.118 0.565 1.000",
   DodgerBlue2 => "0.110 0.525 0.933",
   DodgerBlue3 => "0.094 0.455 0.804",
   DodgerBlue4 => "0.063 0.306 0.545",
   LightBlue => "0.678 0.847 0.902",
   LightBlue1 => "0.749 0.937 1.000",
   LightBlue2 => "0.698 0.875 0.933",
   LightBlue3 => "0.604 0.753 0.804",
   LightBlue4 => "0.408 0.514 0.545",
   LightCyan => "0.878 1.000 1.000",
   LightCyan1 => "0.878 1.000 1.000",
   LightCyan2 => "0.820 0.933 0.933",
   LightCyan3 => "0.706 0.804 0.804",
   LightCyan4 => "0.478 0.545 0.545",
   LightSkyBlue => "0.529 0.808 0.980",
   LightSkyBlue1 => "0.690 0.886 1.000",
   LightSkyBlue2 => "0.643 0.827 0.933",
   LightSkyBlue3 => "0.553 0.714 0.804",
   LightSkyBlue4 => "0.376 0.482 0.545",
   LightSlateBlue => "0.518 0.439 1.000",
   LightSteelBlue => "0.690 0.769 0.871",
   LightSteelBlue1 => "0.792 0.882 1.000",
   LightSteelBlue2 => "0.737 0.824 0.933",
   LightSteelBlue3 => "0.635 0.710 0.804",
   LightSteelBlue4 => "0.431 0.482 0.545",
   MediumAquamarine => "0.400 0.804 0.667",
   MediumBlue => "0.000 0.000 0.804",
   MediumSlateBlue => "0.482 0.408 0.933",
   MediumTurquoise => "0.282 0.820 0.800",
   MidnightBlue => "0.098 0.098 0.439",
   NavyBlue => "0.000 0.000 0.502",
   PaleTurquoise => "0.686 0.933 0.933",
   PaleTurquoise1 => "0.733 1.000 1.000",
   PaleTurquoise2 => "0.682 0.933 0.933",
   PaleTurquoise3 => "0.588 0.804 0.804",
   PaleTurquoise4 => "0.400 0.545 0.545",
   PowderBlue => "0.690 0.878 0.902",
   RoyalBlue => "0.255 0.412 0.882",
   RoyalBlue1 => "0.282 0.463 1.000",
   RoyalBlue2 => "0.263 0.431 0.933",
   RoyalBlue3 => "0.227 0.373 0.804",
   RoyalBlue4 => "0.153 0.251 0.545",
   RoyalBlue5 => "0.000 0.133 0.400",
   SkyBlue => "0.529 0.808 0.922",
   SkyBlue1 => "0.529 0.808 1.000",
   SkyBlue2 => "0.494 0.753 0.933",
   SkyBlue3 => "0.424 0.651 0.804",
   SkyBlue4 => "0.290 0.439 0.545",
   SlateBlue => "0.416 0.353 0.804",
   SlateBlue1 => "0.514 0.435 1.000",
   SlateBlue2 => "0.478 0.404 0.933",
   SlateBlue3 => "0.412 0.349 0.804",
   SlateBlue4 => "0.278 0.235 0.545",
   SteelBlue => "0.275 0.510 0.706",
   SteelBlue1 => "0.388 0.722 1.000",
   SteelBlue2 => "0.361 0.675 0.933",
   SteelBlue3 => "0.310 0.580 0.804",
   SteelBlue4 => "0.212 0.392 0.545",
   aquamarine => "0.498 1.000 0.831",
   aquamarine1 => "0.498 1.000 0.831",
   aquamarine2 => "0.463 0.933 0.776",
   aquamarine3 => "0.400 0.804 0.667",
   aquamarine4 => "0.271 0.545 0.455",
   azure => "0.941 1.000 1.000",
   azure1 => "0.941 1.000 1.000",
   azure2 => "0.878 0.933 0.933",
   azure3 => "0.757 0.804 0.804",
   azure4 => "0.514 0.545 0.545",
   blue => "0.000 0.000 1.000",
   blue1 => "0.000 0.000 1.000",
   blue2 => "0.000 0.000 0.933",
   blue3 => "0.000 0.000 0.804",
   blue4 => "0.000 0.000 0.545",
   cyan => "0.000 1.000 1.000",
   cyan1 => "0.000 1.000 1.000",
   cyan2 => "0.000 0.933 0.933",
   cyan3 => "0.000 0.804 0.804",
   cyan4 => "0.000 0.545 0.545",
   navy => "0.000 0.000 0.502",
   turquoise => "0.251 0.878 0.816",
   turquoise1 => "0.000 0.961 1.000",
   turquoise2 => "0.000 0.898 0.933",
   turquoise3 => "0.000 0.773 0.804",
   turquoise4 => "0.000 0.525 0.545",
   DarkSlateGray => "0.184 0.310 0.310",
   DarkSlateGray1 => "0.592 1.000 1.000",
   DarkSlateGray2 => "0.553 0.933 0.933",
   DarkSlateGray3 => "0.475 0.804 0.804",
   DarkSlateGray4 => "0.322 0.545 0.545",
   RosyBrown => "0.737 0.561 0.561",
   RosyBrown1 => "1.000 0.757 0.757",
   RosyBrown2 => "0.933 0.706 0.706",
   RosyBrown3 => "0.804 0.608 0.608",
   RosyBrown4 => "0.545 0.412 0.412",
   SaddleBrown => "0.545 0.271 0.075",
   SandyBrown => "0.957 0.643 0.376",
   beige => "0.961 0.961 0.863",
   brown => "0.647 0.165 0.165",
   brown1 => "1.000 0.251 0.251",
   brown2 => "0.933 0.231 0.231",
   brown3 => "0.804 0.200 0.200",
   brown4 => "0.545 0.137 0.137",
   burlywood => "0.871 0.722 0.529",
   burlywood1 => "1.000 0.827 0.608",
   burlywood2 => "0.933 0.773 0.569",
   burlywood3 => "0.804 0.667 0.490",
   burlywood4 => "0.545 0.451 0.333",
   chocolate => "0.824 0.412 0.118",
   chocolate1 => "1.000 0.498 0.141",
   chocolate2 => "0.933 0.463 0.129",
   chocolate3 => "0.804 0.400 0.114",
   chocolate4 => "0.545 0.271 0.075",
   peru => "0.804 0.522 0.247",
   tan => "0.824 0.706 0.549",
   tan1 => "1.000 0.647 0.310",
   tan2 => "0.933 0.604 0.286",
   tan3 => "0.804 0.522 0.247",
   tan4 => "0.545 0.353 0.169",
   DarkGreen => "0.000 0.392 0.000",
   DarkKhaki => "0.741 0.718 0.420",
   DarkOliveGreen => "0.333 0.420 0.184",
   DarkOliveGreen1 => "0.792 1.000 0.439",
   DarkOliveGreen2 => "0.737 0.933 0.408",
   DarkOliveGreen3 => "0.635 0.804 0.353",
   DarkOliveGreen4 => "0.431 0.545 0.239",
   DarkSeaGreen => "0.561 0.737 0.561",
   DarkSeaGreen1 => "0.757 1.000 0.757",
   DarkSeaGreen2 => "0.706 0.933 0.706",
   DarkSeaGreen3 => "0.608 0.804 0.608",
   DarkSeaGreen4 => "0.412 0.545 0.412",
   ForestGreen => "0.133 0.545 0.133",
   GreenYellow => "0.678 1.000 0.184",
   LawnGreen => "0.486 0.988 0.000",
   LightSeaGreen => "0.125 0.698 0.667",
   LimeGreen => "0.196 0.804 0.196",
   MediumSeaGreen => "0.235 0.702 0.443",
   MediumSpringGreen => "0.000 0.980 0.604",
   MintCream => "0.961 1.000 0.980",
   OliveDrab => "0.420 0.557 0.137",
   OliveDrab1 => "0.753 1.000 0.243",
   OliveDrab2 => "0.702 0.933 0.227",
   OliveDrab3 => "0.604 0.804 0.196",
   OliveDrab4 => "0.412 0.545 0.133",
   PaleGreen => "0.596 0.984 0.596",
   PaleGreen1 => "0.604 1.000 0.604",
   PaleGreen2 => "0.565 0.933 0.565",
   PaleGreen3 => "0.486 0.804 0.486",
   PaleGreen4 => "0.329 0.545 0.329",
   SeaGreen => "0.180 0.545 0.341",
   SeaGreen1 => "0.329 1.000 0.624",
   SeaGreen2 => "0.306 0.933 0.580",
   SeaGreen3 => "0.263 0.804 0.502",
   SeaGreen4 => "0.180 0.545 0.341",
   SpringGreen => "0.000 1.000 0.498",
   SpringGreen1 => "0.000 1.000 0.498",
   SpringGreen2 => "0.000 0.933 0.463",
   SpringGreen3 => "0.000 0.804 0.400",
   SpringGreen4 => "0.000 0.545 0.271",
   YellowGreen => "0.604 0.804 0.196",
   chartreuse => "0.498 1.000 0.000",
   chartreuse1 => "0.498 1.000 0.000",
   chartreuse2 => "0.463 0.933 0.000",
   chartreuse3 => "0.400 0.804 0.000",
   chartreuse4 => "0.271 0.545 0.000",
   green => "0.000 1.000 0.000",
   green1 => "0.000 1.000 0.000",
   green2 => "0.000 0.933 0.000",
   green3 => "0.000 0.804 0.000",
   green4 => "0.000 0.545 0.000",
   khaki => "0.941 0.902 0.549",
   khaki1 => "1.000 0.965 0.561",
   khaki2 => "0.933 0.902 0.522",
   khaki3 => "0.804 0.776 0.451",
   khaki4 => "0.545 0.525 0.306",
   DarkOrange => "1.000 0.549 0.000",
   DarkOrange1 => "1.000 0.498 0.000",
   DarkOrange2 => "0.933 0.463 0.000",
   DarkOrange3 => "0.804 0.400 0.000",
   DarkOrange4 => "0.545 0.271 0.000",
   DarkSalmon => "0.914 0.588 0.478",
   LightCoral => "0.941 0.502 0.502",
   LightSalmon => "1.000 0.627 0.478",
   LightSalmon1 => "1.000 0.627 0.478",
   LightSalmon2 => "0.933 0.584 0.447",
   LightSalmon3 => "0.804 0.506 0.384",
   LightSalmon4 => "0.545 0.341 0.259",
   PeachPuff => "1.000 0.855 0.725",
   PeachPuff1 => "1.000 0.855 0.725",
   PeachPuff2 => "0.933 0.796 0.678",
   PeachPuff3 => "0.804 0.686 0.584",
   PeachPuff4 => "0.545 0.467 0.396",
   bisque => "1.000 0.894 0.769",
   bisque1 => "1.000 0.894 0.769",
   bisque2 => "0.933 0.835 0.718",
   bisque3 => "0.804 0.718 0.620",
   bisque4 => "0.545 0.490 0.420",
   coral => "1.000 0.498 0.314",
   coral1 => "1.000 0.447 0.337",
   coral2 => "0.933 0.416 0.314",
   coral3 => "0.804 0.357 0.271",
   coral4 => "0.545 0.243 0.184",
   honeydew => "0.941 1.000 0.941",
   honeydew1 => "0.941 1.000 0.941",
   honeydew2 => "0.878 0.933 0.878",
   honeydew3 => "0.757 0.804 0.757",
   honeydew4 => "0.514 0.545 0.514",
   orange => "1.000 0.647 0.000",
   orange1 => "1.000 0.647 0.000",
   orange2 => "0.933 0.604 0.000",
   orange3 => "0.804 0.522 0.000",
   orange4 => "0.545 0.353 0.000",
   salmon => "0.980 0.502 0.447",
   salmon1 => "1.000 0.549 0.412",
   salmon2 => "0.933 0.510 0.384",
   salmon3 => "0.804 0.439 0.329",
   salmon4 => "0.545 0.298 0.224",
   sienna => "0.627 0.322 0.176",
   sienna1 => "1.000 0.510 0.278",
   sienna2 => "0.933 0.475 0.259",
   sienna3 => "0.804 0.408 0.224",
   sienna4 => "0.545 0.278 0.149",
   Color => "0.000 0.000 0.000",
   DeepPink => "1.000 0.078 0.576",
   DeepPink1 => "1.000 0.078 0.576",
   DeepPink2 => "0.933 0.071 0.537",
   DeepPink3 => "0.804 0.063 0.463",
   DeepPink4 => "0.545 0.039 0.314",
   HotPink => "1.000 0.412 0.706",
   HotPink1 => "1.000 0.431 0.706",
   HotPink2 => "0.933 0.416 0.655",
   HotPink3 => "0.804 0.376 0.565",
   HotPink4 => "0.545 0.227 0.384",
   IndianRed => "0.804 0.361 0.361",
   IndianRed1 => "1.000 0.416 0.416",
   IndianRed2 => "0.933 0.388 0.388",
   IndianRed3 => "0.804 0.333 0.333",
   IndianRed4 => "0.545 0.227 0.227",
   LightPink => "1.000 0.714 0.757",
   LightPink1 => "1.000 0.682 0.725",
   LightPink2 => "0.933 0.635 0.678",
   LightPink3 => "0.804 0.549 0.584",
   LightPink4 => "0.545 0.373 0.396",
   MediumVioletRed => "0.780 0.082 0.522",
   MistyRose => "1.000 0.894 0.882",
   MistyRose1 => "1.000 0.894 0.882",
   MistyRose2 => "0.933 0.835 0.824",
   MistyRose3 => "0.804 0.718 0.710",
   MistyRose4 => "0.545 0.490 0.482",
   OrangeRed => "1.000 0.271 0.000",
   OrangeRed1 => "1.000 0.271 0.000",
   OrangeRed2 => "0.933 0.251 0.000",
   OrangeRed3 => "0.804 0.216 0.000",
   OrangeRed4 => "0.545 0.145 0.000",
   PaleVioletRed => "0.859 0.439 0.576",
   PaleVioletRed1 => "1.000 0.510 0.671",
   PaleVioletRed2 => "0.933 0.475 0.624",
   PaleVioletRed3 => "0.804 0.408 0.537",
   PaleVioletRed4 => "0.545 0.278 0.365",
   VioletRed => "0.816 0.125 0.565",
   VioletRed1 => "1.000 0.243 0.588",
   VioletRed2 => "0.933 0.227 0.549",
   VioletRed3 => "0.804 0.196 0.471",
   VioletRed4 => "0.545 0.133 0.322",
   firebrick => "0.698 0.133 0.133",
   firebrick1 => "1.000 0.188 0.188",
   firebrick2 => "0.933 0.173 0.173",
   firebrick3 => "0.804 0.149 0.149",
   firebrick4 => "0.545 0.102 0.102",
   pink => "1.000 0.753 0.796",
   pink1 => "1.000 0.710 0.773",
   pink2 => "0.933 0.663 0.722",
   pink3 => "0.804 0.569 0.620",
   pink4 => "0.545 0.388 0.424",
   red => "1.000 0.000 0.000",
   red1 => "1.000 0.000 0.000",
   red2 => "0.933 0.000 0.000",
   red3 => "0.804 0.000 0.000",
   red4 => "0.545 0.000 0.000",
   tomato => "1.000 0.388 0.278",
   tomato1 => "1.000 0.388 0.278",
   tomato2 => "0.933 0.361 0.259",
   tomato3 => "0.804 0.310 0.224",
   tomato4 => "0.545 0.212 0.149",
   DarkOrchid => "0.600 0.196 0.800",
   DarkOrchid1 => "0.749 0.243 1.000",
   DarkOrchid2 => "0.698 0.227 0.933",
   DarkOrchid3 => "0.604 0.196 0.804",
   DarkOrchid4 => "0.408 0.133 0.545",
   DarkViolet => "0.580 0.000 0.827",
   LavenderBlush => "1.000 0.941 0.961",
   LavenderBlush1 => "1.000 0.941 0.961",
   LavenderBlush2 => "0.933 0.878 0.898",
   LavenderBlush3 => "0.804 0.757 0.773",
   LavenderBlush4 => "0.545 0.514 0.525",
   MediumOrchid => "0.729 0.333 0.827",
   MediumOrchid1 => "0.878 0.400 1.000",
   MediumOrchid2 => "0.820 0.373 0.933",
   MediumOrchid3 => "0.706 0.322 0.804",
   MediumOrchid4 => "0.478 0.216 0.545",
   MediumPurple => "0.576 0.439 0.859",
   MediumPurple1 => "0.671 0.510 1.000",
   MediumPurple2 => "0.624 0.475 0.933",
   MediumPurple3 => "0.537 0.408 0.804",
   MediumPurple4 => "0.365 0.278 0.545",
   lavender => "0.902 0.902 0.980",
   magenta => "1.000 0.000 1.000",
   magenta1 => "1.000 0.000 1.000",
   magenta2 => "0.933 0.000 0.933",
   magenta3 => "0.804 0.000 0.804",
   magenta4 => "0.545 0.000 0.545",
   maroon => "0.690 0.188 0.376",
   maroon1 => "1.000 0.204 0.702",
   maroon2 => "0.933 0.188 0.655",
   maroon3 => "0.804 0.161 0.565",
   maroon4 => "0.545 0.110 0.384",
   orchid => "0.855 0.439 0.839",
   orchid1 => "1.000 0.514 0.980",
   orchid2 => "0.933 0.478 0.914",
   orchid3 => "0.804 0.412 0.788",
   orchid4 => "0.545 0.278 0.537",
   plum => "0.867 0.627 0.867",
   plum1 => "1.000 0.733 1.000",
   plum2 => "0.933 0.682 0.933",
   plum3 => "0.804 0.588 0.804",
   plum4 => "0.545 0.400 0.545",
   purple => "0.627 0.125 0.941",
   purple1 => "0.608 0.188 1.000",
   purple2 => "0.569 0.173 0.933",
   purple3 => "0.490 0.149 0.804",
   purple4 => "0.333 0.102 0.545",
   thistle => "0.847 0.749 0.847",
   thistle1 => "1.000 0.882 1.000",
   thistle2 => "0.933 0.824 0.933",
   thistle3 => "0.804 0.710 0.804",
   thistle4 => "0.545 0.482 0.545",
   violet => "0.933 0.510 0.933",
   AntiqueWhite => "0.980 0.922 0.843",
   AntiqueWhite1 => "1.000 0.937 0.859",
   AntiqueWhite2 => "0.933 0.875 0.800",
   AntiqueWhite3 => "0.804 0.753 0.690",
   AntiqueWhite4 => "0.545 0.514 0.471",
   FloralWhite => "1.000 0.980 0.941",
   GhostWhite => "0.973 0.973 1.000",
   NavajoWhite => "1.000 0.871 0.678",
   NavajoWhite1 => "1.000 0.871 0.678",
   NavajoWhite2 => "0.933 0.812 0.631",
   NavajoWhite3 => "0.804 0.702 0.545",
   NavajoWhite4 => "0.545 0.475 0.369",
   OldLace => "0.992 0.961 0.902",
   WhiteSmoke => "0.961 0.961 0.961",
   gainsboro => "0.863 0.863 0.863",
   ivory => "1.000 1.000 0.941",
   ivory1 => "1.000 1.000 0.941",
   ivory2 => "0.933 0.933 0.878",
   ivory3 => "0.804 0.804 0.757",
   ivory4 => "0.545 0.545 0.514",
   linen => "0.980 0.941 0.902",
   seashell => "1.000 0.961 0.933",
   seashell1 => "1.000 0.961 0.933",
   seashell2 => "0.933 0.898 0.871",
   seashell3 => "0.804 0.773 0.749",
   seashell4 => "0.545 0.525 0.510",
   snow => "1.000 0.980 0.980",
   snow1 => "1.000 0.980 0.980",
   snow2 => "0.933 0.914 0.914",
   snow3 => "0.804 0.788 0.788",
   snow4 => "0.545 0.537 0.537",
   wheat => "0.961 0.871 0.702",
   wheat1 => "1.000 0.906 0.729",
   wheat2 => "0.933 0.847 0.682",
   wheat3 => "0.804 0.729 0.588",
   wheat4 => "0.545 0.494 0.400",
   white => "1.000 1.000 1.000",
   BlanchedAlmond => "1.000 0.922 0.804",
   DarkGoldenrod => "0.722 0.525 0.043",
   DarkGoldenrod1 => "1.000 0.725 0.059",
   DarkGoldenrod2 => "0.933 0.678 0.055",
   DarkGoldenrod3 => "0.804 0.584 0.047",
   DarkGoldenrod4 => "0.545 0.396 0.031",
   LemonChiffon => "1.000 0.980 0.804",
   LemonChiffon1 => "1.000 0.980 0.804",
   LemonChiffon2 => "0.933 0.914 0.749",
   LemonChiffon3 => "0.804 0.788 0.647",
   LemonChiffon4 => "0.545 0.537 0.439",
   LightGoldenrod => "0.933 0.867 0.510",
   LightGoldenrod1 => "1.000 0.925 0.545",
   LightGoldenrod2 => "0.933 0.863 0.510",
   LightGoldenrod3 => "0.804 0.745 0.439",
   LightGoldenrod4 => "0.545 0.506 0.298",
   LightGoldenrodYellow => "0.980 0.980 0.824",
   LightYellow => "1.000 1.000 0.878",
   LightYellow1 => "1.000 1.000 0.878",
   LightYellow2 => "0.933 0.933 0.820",
   LightYellow3 => "0.804 0.804 0.706",
   LightYellow4 => "0.545 0.545 0.478",
   PaleGoldenrod => "0.933 0.910 0.667",
   PapayaWhip => "1.000 0.937 0.835",
   cornsilk => "1.000 0.973 0.863",
   cornsilk1 => "1.000 0.973 0.863",
   cornsilk2 => "0.933 0.910 0.804",
   cornsilk3 => "0.804 0.784 0.694",
   cornsilk4 => "0.545 0.533 0.471",
   gold => "1.000 0.843 0.000",
   gold1 => "1.000 0.843 0.000",
   gold2 => "0.933 0.788 0.000",
   gold3 => "0.804 0.678 0.000",
   gold4 => "0.545 0.459 0.000",
   goldenrod => "0.855 0.647 0.125",
   goldenrod1 => "1.000 0.757 0.145",
   goldenrod2 => "0.933 0.706 0.133",
   goldenrod3 => "0.804 0.608 0.114",
   goldenrod4 => "0.545 0.412 0.078",
   moccasin => "1.000 0.894 0.710",
   yellow => "1.000 1.000 0.000",
   yellow1 => "1.000 1.000 0.000",
   yellow2 => "0.933 0.933 0.000",
   yellow3 => "0.804 0.804 0.000",
   yellow4 => "0.545 0.545 0.000"
} ;

   if (!exists $color2rgb->{$color_name}) {
      return {
         rgb => "0.000 0.000 0.000",
         hexcode => "000000"
      } ;
   }

   my @rgb = split(' ', $color2rgb->{$color_name}) ;

   my $rgb255 = join(",",POSIX::floor($rgb[0] * 255),
                         POSIX::floor($rgb[1] * 255),
                         POSIX::floor($rgb[2] * 255));

   my $hexcode  ;
   foreach my $val (@rgb) {
      my $t_hex = sprintf("%x",POSIX::floor($val * 255)) ;
      if (length($t_hex) == 1) {$t_hex = '0'.$t_hex;}
      $hexcode .= $t_hex ;
   }


   return {
      rgb => $color2rgb->{$color_name},
      rgb255 => $rgb255,
      hexcode => $hexcode
   } ;

}

=head2 pngconvert_bdp_interaction_topology_graph()

   Title:       pngconvert_bdp_interaction_topology_graph()
   Function:    Iterates over all eps format topology graphs and converts to
                 PNG using ImageMagick
   Args:        $_->{pibase_specs} = pibase_specs [optional - if not get_specs()]
   Returns:     Nothing

=cut

sub pngconvert_bdp_interaction_topology_graph {

   my $in = shift ;
   my $pibase_specs ;
   if (!exists $in->{pibase_specs}) {
      $pibase_specs = pibase::get_specs() ;
   } else {
      $pibase_specs = $in->{pibase_specs};
   }
   my $binaries = pibase::locate_binaries() ;

   my ($dbh) = pibase::connect_pibase() ;
   my ($bdp_id, $subset_source_id, $eps_file) = pibase::mysql_fetchcols($dbh,
      "SELECT bdp_id, subset_source_id, file_path ".
      "FROM bdp_interaction_topology_graph") ;
   print STDERR "File to process: ".($#{$bdp_id} + 1)."\n" ;
   print STDERR "Now on file: 0";
   foreach my $j ( 0 .. $#{$bdp_id}) {
      my $counter = $j + 1 ;
      print STDERR "\b"x(length($j)).($j + 1) ;
      my $old_fn = $eps_file->[$j] ;
      my $new_fn = $eps_file->[$j] ; $new_fn =~ s/eps$/png/ ;
      system($binaries->{imagemagick_convert}." ".$old_fn." ".$new_fn) ;
   }
   print STDERR "\n" ;

}

=head2 view_object()

   Title:       view_object()
   Function:    Retrieves and displays to STDOUT PIBASE structure file for
                 domain, interface, or complex
   Args:        $_->{subset_id} = subset_id of domain
                $_->{subset_id_1} = subset_id of domain 1 of an interface
                $_->{subset_id_2} = subset_id of domain 2 of an interface
                $_->{bdp_id} = bdp_id of a complex
                $_->{subset_source_id} = domain classification system of a complex
   Returns:     STDOUT html page of search results

=cut


sub view_object { #routine to retrieve structure files for
#  (1) domains, (2) interactions, (3) complexes w domain defs

   my $in = shift ;
   my $q = new CGI ;

   my $pb_specs = pibase::get_specs() ;

   # get all options into a hash
   my $input = $q->Vars ;
   my $object_type ;
   if (exists $input->{bdp_id} && exists $input->{subset_source_id}) {
      $object_type = 'complex' ;
   } elsif (exists $input->{subset_id_2}) {
      $object_type = 'interface' ;
   } elsif (exists $input->{subset_id}) {
      $object_type = 'domain' ;
   }

   my ($dbh) = pibase::connect_pibase ;

   if ($object_type eq 'complex') {

      print $q->header('chemical/x-ras') ;
      my $bdp_id = $input->{bdp_id} ;
      my $subset_source_id = $input->{subset_source_id} ;

      print "load inline\n" ;
      print "wireframe off\n" ;
      print "set background lightgrey\n" ;
      print "set specular on\n" ;
      my ($subset_id, $class, $subset_fn) = pibase::mysql_fetchcols($dbh,
         "SELECT a.subset_id, a.class, b.file_path FROM subsets as a, ".
         "subsets_files as b WHERE a.bdp_id = \"$bdp_id\" AND ".
         "a.subset_source_id = \"$subset_source_id\" AND ".
         "a.subset_id = b.subset_id") ;

      my $class2color = {};
      my $class_counter = 0;
      foreach my $t_class (sort @{$class}) {
         if (exists $class2color->{$t_class}) {next;}
         $class2color->{$t_class} =
            $pb_specs->{web}->{domain_colors}->[$class_counter] ;
         print "# color for class $t_class is ".$class2color->{$t_class}."\n";
         $class_counter++ ;
      }

      my $cur_colors ; my $sid2domnum ;
      foreach my $j ( 0 .. $#{$subset_id}) {
         $cur_colors->[$j] = $class2color->{$class->[$j]} ;
         $sid2domnum->{$subset_id->[$j]} = $j ;
         print "# color for d$j is ".$cur_colors->[$j]."\n ";
      }

      _view_interface_subset2rasmol({
         dbh => $dbh,
         subsets => $subset_id,
         colors => $cur_colors
      }) ;

      my ($subset_id_1, $subset_id_2) = pibase::mysql_fetchcols($dbh,
         "SELECT a.subset_id_1, a.subset_id_2 FROM intersubset_contacts as a, ".
         "subsets as b WHERE a.bdp_id = \"$bdp_id\" AND ".
         "a.subset_id_1 = b.subset_id AND ".
         "b.subset_source_id = \"$subset_source_id\"") ;
      print "select all\n" ;
      print "ribbons 200\n" ;
      my @interfaces ;
      foreach my $j (0 .. $#{$subset_id_1}) {
         my $s12sig = "i$j" ; push @interfaces, $s12sig ;
         my $s1sig = "d".$sid2domnum->{$subset_id_1->[$j]} ;
         my $s2sig = "d".$sid2domnum->{$subset_id_2->[$j]} ;
         print "define $s12sig ((within(".$pb_specs->{interface_detect_calc}->{dist_cutoff}.", $s1sig) and $s2sig) or (within(".$pb_specs->{interface_detect_calc}->{dist_cutoff}.", $s2sig) and $s1sig))\n" ;
         print "select $s12sig\n" ;
         print "backbone off\n" ;
         print "ribbons\n" ;
      }
      print "define allints (".join(',', @interfaces).")\n" ;
      print "select allints\ndots 200\n" ;
      print "exit\n" ;
      foreach my $j ( 0 .. $#{$subset_fn}) {
         my $t_1 = $pb_specs->{subsets_dir} ;
         my $t_2 = $pb_specs->{web}->{subsets_files_basedir} ;
         $subset_fn->[$j] =~ s/$t_1/$t_2/ ;
         if ($subset_fn->[$j] =~ /\.gz$/) {
            system($pb_specs->{binaries}->{zcat}." ".$subset_fn->[$j]) ;
         } else {
            system("cat ".$subset_fn->[$j]) ;
         }
      }

   } elsif ($object_type eq 'interface') {

      print $q->header('chemical/x-ras') ;

#      print $q->header('text/html') ;
#      print $q->start_html(-title => 'PIBASE view results',
#             -style=>{'src'=>'http://salilab.org/~fred/style/styles.css'}) ;
#      my $greeting = pibase::web::greeting() ;
#      print $greeting ;

      my $subset_id_1 = $input->{subset_id_1} ;
      my $subset_id_2 = $input->{subset_id_2} ;

      my ($subset1_fn) = pibase::mysql_singleval($dbh,
      "SELECT file_path FROM subsets_files WHERE subset_id = \"$subset_id_1\"");
      my ($subset2_fn) = pibase::mysql_singleval($dbh,
      "SELECT file_path FROM subsets_files WHERE subset_id = \"$subset_id_2\"");
      my $t_1 = $pb_specs->{subsets_dir} ;
      my $t_2 = $pb_specs->{web}->{subsets_files_basedir} ;
      $subset1_fn =~ s/$t_1/$t_2/ ;
      $subset2_fn =~ s/$t_1/$t_2/ ;

      my @cur_colors = @{$pb_specs->{web}->{domain_colors}}[0..1] ;

      print "load inline\n" ;
      print "wireframe off\n" ;
      print "set background lightgrey\n" ;
      print "set specular on\n" ;
      _view_interface_subset2rasmol({
         dbh => $dbh,
         subsets => [$subset_id_1, $subset_id_2],
         colors => \@cur_colors
      }) ;
      print "define interface ((within(".$pb_specs->{interface_detect_calc}->{dist_cutoff}.", d1) and not d1) or (within(".$pb_specs->{interface_detect_calc}->{dist_cutoff}.", d0) and not d0))\n" ;
      print "select all\n" ;
      print "ribbons 200\n" ;
      print "select interface\n" ;
      print "backbone off\n" ;
      print "ribbons\n" ;
      print "dots 200\n" ;
#      print "# am now going to print $subset1_fn\n" ;
#      print "# am now going to print $subset2_fn\n" ;
      print "exit\n" ;
      foreach my $t_subset_fn ($subset1_fn, $subset2_fn) {
         if ($t_subset_fn =~ /\.gz$/) {
            system($pb_specs->{binaries}->{zcat}." ".$t_subset_fn) ;
         } else {
            system("cat ".$t_subset_fn) ;
         }
      } 

   } elsif ($object_type eq 'domain') {

      print $q->header('chemical/x-pdb') ;
      my $subset_id = $input->{subset_id} ;
      my ($subset_fn) = pibase::mysql_singleval($dbh,
       "SELECT file_path FROM subsets_files WHERE subset_id = \"$subset_id\"") ;
      my $t_1 = $pb_specs->{subsets_dir} ;
      my $t_2 = $pb_specs->{web}->{subsets_files_basedir} ;
      $subset_fn =~ s/$t_1/$t_2/ ;
      if ($subset_fn =~ /\.gz$/) {
         system($pb_specs->{binaries}->{zcat}." ".$subset_fn) ;
      } else {
         system("cat ".$subset_fn) ;
      }

# make an array of color values based on the classificaitons ala topo graph
   }

}


=head2 view_interface_subset2rasmol()

   Title:       view_object()
   Function:    Returns rasmol script commands to define a series of domains
   Args:        $_->{dbh} = PIBASE database handle
                $_->{subsets} = [subset_id_1, 2, ...]
                $_->{colors} = [name of color to use for domain 1, 2, ...]
   Returns:     STDOUT prints rasmol script commands

=cut

sub _view_interface_subset2rasmol {

   my $in = shift ;
   my $dbh = $in->{dbh} ;
   my $subsets = $in->{subsets} ;
   my $colors = $in->{colors} ;

   foreach my $j ( 0 .. $#{$subsets}) {
      my ($chain, $start, $end) = pibase::mysql_fetchcols($dbh,
         "SELECT chain_id, start_resno, end_resno FROM subsets_details ".
         "WHERE subset_id = \"$subsets->[$j]\"") ;

      my $define_comm ;
      $define_comm = "define d$j ( " ;
      my $lastnum = $#{$chain} ;
      foreach my $k ( 0 .. ($#{$chain} - 1)) {
         $define_comm .= '(' ;
         if (defined $start->[$k]) {
            $define_comm .= " resno >= $start->[$k] " ; }

         if (defined $end->[$j]) {
            $define_comm .= " and resno >= $end->[$k] " ; }

         if (defined $chain->[$j]) {
            $define_comm .= " and *$chain->[$k]" ; }

         $define_comm .= ') OR' ;
      }


      $define_comm .= '(' ;
      if (defined $start->[$lastnum]) {
         $define_comm .= " resno >= $start->[$lastnum] " ; }

      if (defined $end->[$lastnum]) {
         $define_comm .= " and resno <= $end->[$lastnum] " ; }

      if (defined $chain->[$lastnum]) {
         $define_comm .= " and *$chain->[$lastnum]" ; }

      $define_comm .= ')' ;
      $define_comm .= ')' ;
      print $define_comm."\n" ;
   }

   foreach my $j ( 0 .. $#{$subsets}) {
      print  "select d$j\n" ;
      my $t_color = pibase::web::get_color_codes({name => $colors->[$j]}) ;
      print "color [".$t_color->{rgb255}."]\n" ; }

}

1 ;
