#!/usr/local/bin/perl
=head1 NAME

sqlcreate.parser.pl

=head1 DESCRIPTION

Parses (My)SQL CREATE TABLE statements and displays perl code that
defines table structure

=head1 SYNOPSIS

sqlcreate.parser.pl < createtables.SQL

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
use pibase::create_raw_table_specs qw/create_raw_table_specs/ ;
main() ;

sub main {

   pibase::create_raw_table_specs::create_raw_table_specs() ;

}
