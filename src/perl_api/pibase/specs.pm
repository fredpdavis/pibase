=head1 NAME

pibase::specs - Perl package containing definition of pibase table structures

=head1 DESCRIPTION

Contains pibase table structure definitions

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

package pibase::specs ;
use strict;
use warnings;
use Exporter;
my @ISA = qw/Exporter/ ;
my @EXPORT_OK = qw/table_spec sql_table_spec/ ;

use pibase::raw_table_specs qw/full_table_specs/ ;

=head2 table_spec(@tablelist)

   Title:       table_spec()
   Function:    Return mysql DDL format table specs.
   Returns:     $_ hashref pointing from tablename to specs
                  $_->{i} = specs for ith table

=cut

sub table_spec {

   my @tables = @_ ;
   my $all_tables = pibase::raw_table_specs::full_table_specs() ;

   my $ans ;

   foreach my $j ( 0 ..$#tables) {
      if (exists $all_tables->{$tables[$j]}) {
         push @{$ans}, {name => $tables[$j], specs => $all_tables->{$tables[$j]}} ;
      }
   }

   return $ans ;

}



=head2 SUB sql_table_spec(@tablelist)

   Title:       sql_table_spec()
   Function:    Return mysql DDL format table specs.
   Args:        $_ hashref pointing from tablename to specs
                $_->{i} = specs for ith table

=cut

sub sql_table_spec {

   my @tables = @_ ;
   my $all_tables = pibase::raw_table_specs::full_table_specs() ;

   my $ans ;

   foreach my $j ( 0 ..$#tables) {
      if (exists $all_tables->{$tables[$j]}) {
         my $cursql ;
         $cursql = '('."\n" ;

	 foreach my $k ( 0 .. $#{$all_tables->{$tables[$j]}->{field_name}}) {
	    $cursql .= "   ".$all_tables->{$tables[$j]}->{field_name}->[$k]." ".$all_tables->{$tables[$j]}->{field_spec}->[$k] ;
	    if ( $k < $#{$all_tables->{$tables[$j]}->{field_name}} ||
	         exists $all_tables->{$tables[$j]}->{prikey} ) {
	       $cursql .= ",\n" ;
	    } else {
	       $cursql .= "\n" ;
	    }
	 }

	 if (exists $all_tables->{$tables[$j]}->{prikey}) {
	    $cursql .= '   PRIMARY KEY '.$all_tables->{$tables[$j]}->{prikey}."\n" ; }

	 $cursql .= ")" ;
         $ans->{$tables[$j]} = $cursql ;
      }
   }

   return $ans ;

}

1;
