#!/usr/local/bin/perl
=head1 NAME

setuprun.pl - script to setup an SGE cluster run

=head1 DESCRIPTION

Splits an input file into a specified number of tasks
and displays a tasklist for inclusion in an SGE job script

=head1 SYNOPSIS

setuprun.pl inputfile output_tasklist

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

use strict;
use warnings;
use pibase::auxil ;

main() ;

sub main {

   pibase::auxil::setup_clusterrun() ;

}
