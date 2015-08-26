=head1 NAME

pibase::tables_on_disk - Perl module that provides an SQL like query interface to tables stored on disk.

=head1 DESCRIPTION

Perl module to interface with tables stored on disk.

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

package pibase::tables_on_disk ;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/contacts_select/ ;

use File::Temp qw/tempfile/ ;
use Sys::Hostname ;
use pibase::specs ;


=head2 select_tod()

   Title:       select_tod()
   Function:    provides a pseudo-mysql select like interface to [gzipped]
                contacts file in pibase.interatomic_contacts_prototype format

   Args:        $_[0]- source file
                $_[1]- arrayref ['bdp_id', 'resno_1', 'resno_2', 'distance']
                $_[2]- table name
                $_[3]- where clause

   Returns:     @_ - array of arrays
                $_[i][j] = jth field of ith result

=cut

sub select_tod {

   my ($source_file, $fields, $table_name, $whereclause) = @_ ;

# Set uncompression binary

   my $hostname = hostname() ;
   my $zcat_bin = 'gzcat' ;
   if ($hostname =~ /^node/) {
      $zcat_bin = 'zcat' ;
   } elsif ( ($hostname =~ /alto/) || ($hostname =~ /diva/)) {
      $zcat_bin = 'gzcat' ;
   }

# Determine table specs.
# Buid a hash pointing from field name to field number

   my $table_spec = table_spec($table_name) ;
   my @format ;
   push @format, @{$table_spec->[0]->{specs}->{field_name}} ;

   my $format ;
   foreach my $j (0 .. $#format) {
      $format->{$format[$j]} = $j + 1 ; }

   my $gz_fl = 0 ;
   if ($source_file =~ /gz$/) {
      $gz_fl = 1; }

   my @extract_fields ;
   my $count_fl = 0 ;
   foreach my $j ( 0 .. $#{$fields}) {
      if (exists $format->{$fields->[$j]}) {
         push @extract_fields, '$'.$format->{$fields->[$j]} ;
      } elsif ($fields->[$j] eq 'count') {
         $count_fl = 1 ; }
   }

   my $awk_where = '';
   if (defined $whereclause) {
      my @raw_wheres = split(/\b[Aa][Nn][Dd]\b/, $whereclause) ;
      my @awk_where ;
      foreach my $j ( 0 .. $#raw_wheres) {
         $raw_wheres[$j] =~ s/^ *//g ;
         $raw_wheres[$j] =~ s/ *$//g ;
         my ($tfield, $trelop, $tval) = split(/\s+/, $raw_wheres[$j]) ;
         $tfield =~ s/ //g ;
         if ($trelop eq '=') {
            $trelop = '==' ; }
         push @awk_where, '( $'.$format->{$tfield}.$trelop.$tval.' )' ;
      }
      $awk_where = 'if ('.join('&&', @awk_where).')' ;
   }

   my $tcom ;
   if ($gz_fl) {
      $tcom = "$zcat_bin " ; }
   else {
      $tcom = "cat " ; }
   $tcom .= $source_file.' | ' ;

   $tcom .= "awk -F\'	\' \'{" ;
   $tcom .= " $awk_where " ;
   $tcom .= "print ".join("\"	\"", @extract_fields)."}\' " ;

   $tcom .= ' | sort | uniq ' ;
   if ($count_fl) {
      $tcom .= "-c " }

   my @results ;
   open(RESULTS, "$tcom | ") ;
   while (my $line = <RESULTS>) {
      if ($line !~ /^#/) {

         chomp $line ;
         my $curcount ;
         if ($count_fl) {
	    ($curcount) = ($line =~ /^ *([0-9]+)/) ;
	    $line =~ s/ *[0-9]* *// ;
	    $line =~ s/^\s+// ;
	 }

         my @t = split(/\t/, $line) ;
         foreach my $j ( 0 .. $#t) {
            push @{$results[$j]}, $t[$j] ; }

         if ($count_fl) {
            push @{$results[($#t + 1)]}, $curcount; }

      }
   }
   close(RESULTS) ;

   return @results ;

}

1;
