=head1 NAME

pibase::interatomic_contacts - perl module that deals with interatomic contacts

=head1 DESCRIPTION

The interatomic_contacts.pm module deals with queries involving interatomic contacts.

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

package pibase::interatomic_contacts ;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/contacts_select_inter contacts_select special_contact raw_contacts_select/ ;

use File::Temp qw/tempfile/ ;
use Sys::Hostname ;
use pibase qw/complete_pibase_specs/;
use pibase::specs qw/table_spec/ ;


=head2 contacts_select()

   Function:    provides a pseudo-mysql select like interface to [gzipped]
                  contacts file in pibase.interatomic_contacts_prototype format
   Args:        $fields - arrayref: ['bdp_id', 'resno_1', 'resno_2', 'distance']
   Return:      results

=cut

sub contacts_select {

   my ($fields, $source_file, $whereclause) = @_ ;

#Set uncompression binary

   my $hostname = hostname() ;
   my $zcat_bin = 'gzcat' ;
   if ($hostname =~ /^node/) {
      $zcat_bin = 'zcat' ;
   } elsif ($hostname =~ /^tombak/) {
      $zcat_bin = 'zcat' ;
   } elsif ($hostname =~ /^daf/) {
      $zcat_bin = 'zcat' ;
   } elsif ( ($hostname =~ /alto/) || ($hostname =~ /diva/)) {
      $zcat_bin = 'gzcat' ;
   }

#Set the interatomic contacts source file format.

   my $specs = pibase::complete_pibase_specs() ;
   my $tablespecs = table_spec('interatomic_contacts_prototype') ;
   my @format = @{$tablespecs->[0]->{specs}->{field_name}} ;


#   my @format = qw/bdp_id/ ;
#   push @format, qw/chain_id_1 resno_1 resno_1_int resna_1 atomna_1/ ;
#   push @format, qw/chain_id_2 resno_2 resno_2_int resna_2 atomna_2/ ;
#   push @format, qw/distance/ ;

#Buid a hash pointing from field name to field number

   my $format ;
   foreach my $j (0 .. $#format) {fjj
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


=head2 raw_contacts_select()

   Title:       raw_contacts_select()
   Function:    Allows SQL-like SELECT from the raw interatomic contacts
                 tables stored on disk
   Args:        $_[0] = source_file
                $_[1] = SQL-like SELECT query
                $_[2] = params
                  $_[2]->{maxdist} - distance cutoff
   Returns:     $_ - filehandle to of results file

=cut

sub raw_contacts_select {

   my ($source_file, $select_sql, $params) = @_ ;

   $select_sql =~ s/^SELECT // ;
   $select_sql =~ s/FROM.*$// ;
   $select_sql =~ s/ //g ;
   my @fields = split(/\,/, $select_sql) ;

#Set uncompression binary

   my ($t0, $tsec) ;

   my $hostname = hostname() ;
   my $zcat_bin = 'gzcat' ;
   if ($hostname =~ /^node/) {
      $zcat_bin = 'zcat' ;
   } elsif ($hostname =~ /^tombak/) {
      $zcat_bin = 'zcat' ;
   } elsif ($hostname =~ /^daf/) {
      $zcat_bin = 'zcat' ;
   } elsif ( ($hostname =~ /alto/) || ($hostname =~ /diva/)) {
      $zcat_bin = 'gzcat' ;
   }

#Set the interatomic contacts source file format.

   my $specs = pibase::complete_pibase_specs() ;
   my $tablespecs = table_spec('interatomic_contacts_prototype') ;
   my @format = @{$tablespecs->[0]->{specs}->{field_name}} ;

#Buid a hash pointing from field name to field number

   my $format ;
   foreach my $j (0 .. $#format) {
      $format->{$format[$j]} = $j + 1 ; }

   my $gz_fl = 0 ;
   if ($source_file =~ /gz$/) {
      $gz_fl = 1; }

   my @extract_fields ;
   my $count_fl = 0 ;
   foreach my $j ( 0 .. $#fields) {
      if (exists $format->{$fields[$j]}) {
         push @extract_fields, '$'.$format->{$fields[$j]} ; } }

   push @extract_fields, '$'.$format->{chain_id_1} ;
   push @extract_fields, '$'.$format->{resno_1} ;
   push @extract_fields, '$'.$format->{chain_id_2} ;
   push @extract_fields, '$'.$format->{resno_2} ;

   my $awk_where = '';
   if (defined $params->{maxdist}) {
      $awk_where = ' if ($'.$format->{distance}." <= $params->{maxdist} ) " ; }

   my $tcom ;
   if ($gz_fl) {
      $tcom = "$zcat_bin " ;
   } else {
      $tcom = "cat " ; }
   $tcom .= $source_file.' | ' ;

   $tcom .= "awk -F\'	\' \'{" ;
   $tcom .= $awk_where ;
   $tcom .= "print ".join("\"	\"", @extract_fields)."}\' " ;

   my @results ;
   my $contacts_fh ;
   open ($contacts_fh, "$tcom | ") ;

   return $contacts_fh ;
}


=head2 contacts_select_inter()

   Title:       contacts_select_inter()
   Function:    Returns inter-domain interatomic-contacts as specified by an
                  SQL-like SELECT query from the raw interatomic contacts tables
                  stored on disk
   Args:        $_[0] = source_file
                $_[1] = SQL-like SELECT query
                $_[2] = resno_2_subset - hash from residue number to domain
                  identifier
                $_[3] = params
                  $_[3]->{maxdist} - distance cutoff
   Returns:     @_ = array of results 

=cut

sub contacts_select_inter {

   my ($source_file, $select_sql, $resno_2_subset, $params) = @_ ;

   $select_sql =~ s/^SELECT // ;
   $select_sql =~ s/FROM.*$// ;
   $select_sql =~ s/ //g ;
   my @fields = split(/\,/, $select_sql) ;

#Set uncompression binary

   my ($t0, $tsec) ;
#   print STDERR "   misc prepartory shit\t" ;

   my $hostname = hostname() ;
   my $zcat_bin = 'gzcat' ;
   if ($hostname =~ /^node/) {
      $zcat_bin = 'zcat' ;
   } elsif ($hostname =~ /^tombak/) {
      $zcat_bin = 'zcat' ;
   } elsif ($hostname =~ /^daf/) {
      $zcat_bin = 'zcat' ;
   } elsif ( ($hostname =~ /alto/) || ($hostname =~ /diva/)) {
      $zcat_bin = 'gzcat' ;
   }

#Set the interatomic contacts source file format.

   my $specs = pibase::complete_pibase_specs() ;
   my $tablespecs = table_spec('interatomic_contacts_prototype') ;
   my @format = @{$tablespecs->[0]->{specs}->{field_name}} ;

#Buid a hash pointing from field name to field number

   my $format ;
   foreach my $j (0 .. $#format) {
      $format->{$format[$j]} = $j + 1 ; }

   my $gz_fl = 0 ;
   if ($source_file =~ /gz$/) {
      $gz_fl = 1; }

   my @extract_fields ;
   my $count_fl = 0 ;
   foreach my $j ( 0 .. $#fields) {
      if (exists $format->{$fields[$j]}) {
         push @extract_fields, '$'.$format->{$fields[$j]} ; } }

   push @extract_fields, '$'.$format->{chain_id_1} ;
   push @extract_fields, '$'.$format->{resno_1} ;
   push @extract_fields, '$'.$format->{chain_id_2} ;
   push @extract_fields, '$'.$format->{resno_2} ;

   my $awk_where = '';
   if (defined $params->{maxdist}) {
      $awk_where = ' if ($'.$format->{distance}." <= $params->{maxdist} ) " ; }

   my $tcom ;
   if ($gz_fl) {
      $tcom = "$zcat_bin " ;
   } else {
      $tcom = "cat " ; }
   $tcom .= $source_file.' | ' ;

   $tcom .= "awk -F\'	\' \'{" ;
   $tcom .= $awk_where ;
   $tcom .= "print ".join("\"	\"", @extract_fields)."}\' " ;

   my @results ;
#   print STDERR "   open interatomic files:" ;
   open(RESULTS, "$tcom | ") ;
#   print STDERR "   parse in interdomain:" ;
   while (my $line = <RESULTS>) {
      if ($line !~ /^#/) {

         chomp $line ;
         my $curcount ;

         my @t = split(/\t/, $line) ;

	 my $sig1 = $t[-3]."\n".$t[-4] ;
	 my $sig2 = $t[-1]."\n".$t[-2] ;

	 if (!(exists $resno_2_subset->{$sig1} &&
	       exists $resno_2_subset->{$sig2}) ||
	     ($resno_2_subset->{$sig1} ne $resno_2_subset->{$sig2})) {
            foreach my $j ( 0 .. ($#t - 4)) {
               push @{$results[$j]}, $t[$j] ; }
         }
      }
   }
   close(RESULTS) ;

   return @results ;
}


=head2 special_params()

   Title:       special_params()
   Function:    Specifies the parameters (like distance thresholds) for the
                  ``special'' contacts (salt bridges, hydrogen bonds, strong
                  hydrogen bonds, disulfide bonds)
   Args:        $_[0] = source_file
                $_[1] = SQL-like SELECT query
                $_[2] = resno_2_subset - hash from residue number to domain
                  identifier
                $_[3] = params
                  $_[3]->{maxdist} - distance cutoff
   Returns:     @_ = array of results 

=cut

sub special_params {

   my $res ;
   my ($thresh, $hbres, $sbres) ;
   $thresh->{salt} = 4 ;
   $thresh->{hbond} = 3.5 ;
   $thresh->{hbond_s} = 4.0 ;
   $thresh->{ssbond} = 3.0 ;

   $res->{thresh} = $thresh ;


# 0 = HB donora hnd acceptor
# 1 = HB donor
# -1 = HB acceptor

   $hbres->{' N  '} = 1 ;
   $hbres->{'HIS'}->{' ND1'} = 1 ;
   $hbres->{'HIS'}->{' NE2'} = 1 ;
   $hbres->{'ASN'}->{' ND2'} = 1 ; #SC flip degeneracy
   $hbres->{'GLN'}->{' NE2'} = 1 ; #SC flip degeneracy
   $hbres->{'ARG'}->{' NE '} = 1 ;
   $hbres->{'ARG'}->{' NH1'} = 1 ;
   $hbres->{'ARG'}->{' NH2'} = 1 ;
   $hbres->{'LYS'}->{' NZ '} = 1 ;
   $hbres->{'TRP'}->{' NE1'} = 1 ;

   $hbres->{'SER'}->{' OG '} = 0 ; # both acceptor and donor
   $hbres->{'THR'}->{' OG1'} = 0 ;
   $hbres->{'TYR'}->{' OH '} = 0 ;

   $hbres->{' O  '} = -1 ;
   $hbres->{' OXT'} = -1 ;
   $hbres->{'ASP'}->{' OD1'} = -1 ;
   $hbres->{'ASP'}->{' OD2'} = -1 ;
   $hbres->{'GLU'}->{' OE1'} = -1 ;
   $hbres->{'GLU'}->{' OE2'} = -1 ;
   $hbres->{'ASN'}->{' OD1'} = 0 ; #SC flip degeneracy
   $hbres->{'GLN'}->{' OE1'} = 0 ; #SC flip degeneracy
   $hbres->{'CYS'}->{' SG '} = 0 ;
   $hbres->{'MET'}->{' SD '} = -1 ;

   $sbres->{' OXT'} = -1 ;
   $sbres->{'ASP'}->{' OD1'} = -1 ;
   $sbres->{'ASP'}->{' OD2'} = -1 ;
   $sbres->{'GLU'}->{' OE1'} = -1 ;
   $sbres->{'GLU'}->{' OE2'} = -1 ;
#   $sbres->{' N  '} = 1 ;
   $sbres->{'HIS'}->{' ND1'} = 1 ;
   $sbres->{'HIS'}->{' NE2'} = 1 ;
   $sbres->{'LYS'}->{' NZ '} = 1 ;
   $sbres->{'ARG'}->{' NH1'} = 1 ;

   $res->{hbond} = $hbres ;
   $res->{salt} = $sbres ;

   return $res ;

}


=head2 special_contact()

   Title:       special_contact()
   Function:    given the distance between a pair of atom/residues, decide
                  whether it meets dist requirements for a hbond, ssbond,
                  or saltbridge
   Args:        $_[0] = contacts information
                 ->{resna1} = name of residue 1
                 ->{atomna1} = name of atom 1
                 ->{resna2} = name of residue 2
                 ->{atomna2} = name of atom 2
                 ->{dist} = distance
                $_[1] = $t - Time::Benchmark timer handle
   Returns:     $_ = contact category (none, salt, ssbond, hbond)

=cut

sub special_contact {

   my $val = shift ;
#   my $t = shift ;

#   $t->start('special_contact preset') ;
   my $thresh ;
   $thresh->{salt} = 4 ;
   $thresh->{hbond} = 3.5 ;
   $thresh->{ssbond} = 3.0 ;

   if ( ($val->{resna1} eq 'CYS') &&
        ($val->{resna2} eq 'CYS') ) {
      if ( ($val->{atomna1} =~ /SG/) &&
           ($val->{atomna2} =~ /SG/) &&
	   ($val->{dist} <= $thresh->{ssbond} ) ) {
         return 'ssbond' ;
      } else {
         return 'none' ;
      }
   }

   if ($val->{atomna1} eq ' SG ' ||
       $val->{atomna2} eq ' SG ' ||
       $val->{atomna1} eq ' SD ' ||
       $val->{atomna2} eq ' SD ' ) {
      $thresh->{hbond} = 4 ;
   }

   my ($resspec, $resind) ;
   $resind->{' N  '} = 'HBD' ;
   $resspec->{'HIS'}->{' ND1'} = 'HBD' ;
   $resspec->{'HIS'}->{' NE2'} = 'HBD' ;
   $resspec->{'ASN'}->{' ND2'} = 'HBb' ; #SC flip degeneracy
   $resspec->{'GLN'}->{' NE2'} = 'HBb' ; #SC flip degeneracy
   $resspec->{'ARG'}->{' NE '} = 'HBD' ;
   $resspec->{'ARG'}->{' NH1'} = 'HBD' ;
   $resspec->{'ARG'}->{' NH2'} = 'HBD' ;
   $resspec->{'LYS'}->{' NZ '} = 'HBD' ;
   $resspec->{'TRP'}->{' NE1'} = 'HBD' ;

   $resspec->{'SER'}->{' OG '} = 'HBb' ; # both acceptor and donor
   $resspec->{'THR'}->{' OG1'} = 'HBb' ;
   $resspec->{'TYR'}->{' OH '} = 'HBb' ;

   $resind->{' O  '} = 'HBA' ;
   $resind->{' OXT'} = 'HBA' ;
   $resspec->{'ASP'}->{' OD1'} = 'HBA' ;
   $resspec->{'ASP'}->{' OD2'} = 'HBA' ;
   $resspec->{'GLU'}->{' OE1'} = 'HBA' ;
   $resspec->{'GLU'}->{' OE2'} = 'HBA' ;
   $resspec->{'ASN'}->{' OD1'} = 'HBb' ; #SC flip degeneracy
   $resspec->{'GLN'}->{' OE1'} = 'HBb' ; #SC flip degeneracy

   $resspec->{'CYS'}->{' SG '} = 'HBb' ;
   $resspec->{'MET'}->{' SD '} = 'HBA' ;
#   $t->stop('special_contact preset') ;


   my $flag ;

#   $t->start('special_contact assign') ;
   my @atomna = ($val->{atomna1}, $val->{atomna2}) ;
   my @resna = ($val->{resna1}, $val->{resna2}) ;
   my $sbposs ;
   foreach my $j ( 0 .. $#atomna) {
      if (exists $resind->{$atomna[$j]}) {
         $flag->[$j]->{$resind->{$atomna[$j]}}++ ; }

      if (exists $resspec->{$resna[$j]}->{$atomna[$j]}) {
         $flag->[$j]->{$resspec->{$resna[$j]}->{$atomna[$j]}}++ ; }

      if ($resna[$j] eq 'ASP' ||
           $resna[$j] eq 'GLU') {
         $sbposs->[$j] = 1 ;
      } elsif ( $resna[$j] eq 'ARG' ||
                $resna[$j] eq 'LYS' ||
                $resna[$j] eq 'HIS' ) {
         $sbposs->[$j] = -1 ;
      } else {
         $sbposs->[$j] = 0 ;
      }
   }
#   $t->stop('special_contact assign') ;


#   $t->start('special_contact assign pair') ;
   if (((exists $flag->[0]->{HBA} && $sbposs->[0] == 1 &&
         exists $flag->[1]->{HBD} && $sbposs->[1] == -1 ) ||
        (exists $flag->[1]->{HBA} && $sbposs->[1] == 1 &&
         exists $flag->[0]->{HBD} && $sbposs->[0] == -1 )) &&
       ($val->{dist} <= $thresh->{salt}) ) {
      return 'salt' ;
   }

   if ( ((($flag->[0]->{HBA} || $flag->[0]->{HBb}) &&
          ($flag->[1]->{HBD} || $flag->[1]->{HBb}) ) ||
         (($flag->[0]->{HBD} || $flag->[0]->{HBb}) &&
          ($flag->[1]->{HBA} || $flag->[1]->{HBb}) ) ) &&
        ($val->{dist} <= $thresh->{hbond})) {
      return 'hbond' ;
   }
#   $t->stop('special_contact assign pair') ;

   return 'none' ;
}

1;
