=head1 NAME

pibase::benchmark - Perl package for benchmarking pibase

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

package pibase::benchmark;
require Exporter;
@ISA = qw/Exporter/ ;
@EXPORT = qw/memusage/ ;

use strict;
use warnings;

=head2 memusage()

   Title:       memusage()
   Function:    Returns the VmSize of the process
      works on linux /proc/$$/status

=cut

sub memusage {
   my $memusage = `grep ^VmSize /proc/$$/status` ;
   print STDERR "memsusage is $memusage\n" ;
   chop $memusage ; $memusage =~ s/^VmSize:\s+//g ;
   return $memusage ;
}

1;
