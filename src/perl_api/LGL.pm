=head1 NAME

LGL.pm - perl interface to Alex Adai's Large Graph Layout (LGL) program

=head1 DESCRIPTION

The LGL.pm package provides a perl interface to LGL so that edge/node lists can
be converted into postscript or png images with user-configurable edge and node
colors, shapes, and sizes.

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


=head1 SYNOPSIS

=head2 Example 1

   my $edges = {
      a => { 'b' => 1} ,
      a => { 'd' => 1} ,
      a => { 'e' => 1} ,
      b => { 'c' => 1} ,
      c => { 'e' => 1} ,
   } ;
   my $coords = LGL::edges2coords({edges => $edges}) ;
   coords2ps({
      edges => $edges,
      coords => $coords,
   }) ;


=head2 Example 2

   my $edges ;
   $edges->{a}->{b} = 1 ;
   $edges->{a}->{c} = 1 ;
   $edges->{b}->{d} = 1 ;

   my $ncols;
   $ncols->{"a"} = "black" ;
   $ncols->{"b"} = "yellow";
   $ncols->{"c"} = "cyan";
   $ncols->{"d"} = "red";

   my $ecols ;
   $ecols->{a}->{b} = "black" ;
   $ecols->{a}->{c} = "purple" ;
   $ecols->{b}->{d} = "brown" ;

   my $coords = LGL::edges2coords({edges => $edges}) ;
   LGL::coords2ps({coords => $coords,
                   edges => $edges,
                   mag => 100,
                   ethick => 1,
                   nrad => 5,
                   nshape => 'fbox',
                   ncol => "0.5 0.3 0.2",
                   ncols => $ncols,
                   ecols => $ecols,
                   ecol => "0.2 0.7 0.4"}) ;


=head1 SUBROUTINES

=cut

package LGL ;
use Exporter ;
use strict;
use warnings;

use Sys::Hostname ;
use File::Copy ;
use File::Temp qw/tempfile tempdir/;
use File::Path qw/rmtree/ ;

BEGIN{
   our $VERSION = "1.200712" ;
}

our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/edges2coords coords2ps/ ;

my $lgl_specs = { 
   loc => "/groups/eddy/home/davisf/work/pibase/pibase200708/auxil/lgl",
   progs => ["lglayout2D", "LGLFormatHandler.pm", "lgl.pl",
             "ParseConfigFile.pm", 'lglbreakup', 'lglrebuild']
} ;

# don't remember where these came from - think it is from X.Org
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


=head2 edges2coords()

   Title:       edges2coords()
   Function:    Performs graph layout using LGL, returning coordinate info
   Args:        $_->{edges}->{node1}->{node2} ;
   Returns:     $_->{node} = [$x, $y] ;

=cut

sub edges2coords {

   my $params = shift ;
   my $edges = $params->{edges} ;

   my $specs = $lgl_specs ;

   my $edgeseen ;
   my $node2ind ;my @nodes; my $ncount = 0 ;

   my $orig_dir = `pwd` ; chomp $orig_dir ;
   my $run_dir = tempdir("lglrunXXXXX") ;
#   my $run_dir = "lglrunXXXXX$$" ; mkdir $run_dir ;
   chdir $run_dir ;

   my ($ncol_fh, $ncol_fn) = tempfile("edgesXXXXX", SUFFIX => ".ncol") ;

   if (!exists $params->{type} || $params->{type} eq 'hash') {
      foreach my $n1 (keys %{$edges}) {
         if (!exists $node2ind->{$n1}) {
            push @nodes, $n1 ; $node2ind->{$n1} = $#nodes ; }

         foreach my $n2 (keys %{$edges->{$n1}}) {
            if (!exists $node2ind->{$n2}) {
               push @nodes, $n2 ; $node2ind->{$n2} = $#nodes ; }

            if ($n1 eq $n2) {next; }
            if (!exists $edgeseen->{$n1}->{$n2}) {
               print $ncol_fh "$node2ind->{$n1} $node2ind->{$n2}\n" ;
               $edgeseen->{$n1}->{$n2}++ ;
               $edgeseen->{$n2}->{$n2}++ ;
            }
         }
      }
   }
   close($ncol_fh) ;


   my ($cnfg_fh, $cnfg_fn) =
      tempfile("lglconfigXXXXX",SUFFIX => ".config") ;
   my ($coords_fh, $coords_fn) =
      tempfile("finalcoordsXXXXX",SUFFIX => ".coords") ;
   close ($coords_fh) ; unlink $coords_fn ;
   my $tmp_dir = tempdir("lgltempXXXXX") ;
#   my $tmp_dir = "lgltempXXXXX$$" ; mkdir $tmp_dir ;

   my $lglparams = {
      ncol_fn => $ncol_fn,
      cnfg_fh => $cnfg_fh,
      tmp_dir => $tmp_dir,
      coords_fn => $coords_fn
   } ;

   _make_lglconfig($lglparams) ;
   close($cnfg_fh) ;

   foreach my $prog (@{$specs->{progs}}) {
      system("cp $specs->{loc}/$prog ./") ; }

   system("perl lgl.pl -c $cnfg_fn 2>/dev/null >/dev/null") ;

   my $coords ;
   if (-s $cnfg_fn) {
      open (COORDS, "$tmp_dir/$coords_fn") ;
      while (my $line = <COORDS>) {
         chomp $line;
         my ($nind, $x, $y) = split(/ /, $line) ;
         $coords->{$nodes[$nind]} = [$x, $y] ;
      }
      close(COORDS) ;
   } else {
      print STDERR "BUG: no final coords output\n" ;
   }

   chdir $orig_dir ;
   rmtree($run_dir) ;
   return $coords ;
}


=head2 _make_lglconfig()

   Title:       _make_lglconfig
   Function:    Makes a config file for an LGL run.
   Returns:     nothing
   Args:        $_->{cnfg_fh} (file handle to display the configuration file to)

=cut

sub _make_lglconfig {

   my $params = shift ;

   print {$params->{cnfg_fh}} "tmpdir = '$params->{tmp_dir}'\n" ;
   print {$params->{cnfg_fh}} "inputfile = '$params->{ncol_fn}'\n" ;
   print {$params->{cnfg_fh}} "finaloutcoords = '$params->{coords_fn}'\n" ;
   print {$params->{cnfg_fh}} "treelayout = '0'\n" ;
   print {$params->{cnfg_fh}} "useoriginalweights = '0'\n" ;
   print {$params->{cnfg_fh}} "edgelevelmap = '1'\n" ;
   print {$params->{cnfg_fh}} "outputmst = '0'\n" ;
   print {$params->{cnfg_fh}} "threadcount = '1'\n" ;
   print {$params->{cnfg_fh}} "dimension = '2'\n" ;
   print {$params->{cnfg_fh}} "cutoff = ''\n" ;
   print {$params->{cnfg_fh}} "usemst = '0'\n" ;
   print {$params->{cnfg_fh}} "issilent = '0'\n" ;
   print {$params->{cnfg_fh}} "pickupdir = ''\n" ;
   print {$params->{cnfg_fh}} "integratetype = ''\n" ;
   print {$params->{cnfg_fh}} "placmentdistance = ''\n" ;
   print {$params->{cnfg_fh}} "placementradius = ''\n" ;
   print {$params->{cnfg_fh}} "placeleafsclose = '0'\n" ;

   return ;

}


=head2 coords2pngmap()

   Title:	coords2pngmap
   Function:	Creates a client-side png image map given coordinate and
                  edge information. 
   Returns:	nothing
   Args:
      - $_->{ncol}	node color
      - $_->{ecol}	edge color
      - $_->{estyle}	edge style (solid, dashed)
      - $_->{nshape}	node shape (triangle, square, circle)
      - $_->{nrad}	node radius
      - $_->{ethick}	edge thickness
      - $_->{mag}	magnification parameter
      - $_->{coords}	node coordinate information
      - $_->{edges}->{n1}->{n2}	edge list (nhash)
      - $_->{map_fh}	output imagemap file name
      - $_->{nodenames}	node names
      - $_->{nodeurls}	node URLs

=cut

sub coords2pngmap {

   my $input = shift ;

   my $params = {
      ncol => "0 0 0",
      ecol => "0 0 0",
      estyle => "solid",
      nshape => "fcircle",
      nrad => .65,
      ethick => 0.3,
      mag => 7
   } ;

   foreach my $prop (keys %{$input}) {
      $params->{$prop} = $input->{$prop} ; }

   my $coords = $input->{coords};
   my $edges = $input->{edges};

   my $bbox ;
   my $new_coords ;

   foreach my $node (keys %{$coords}) {
      $new_coords->{$node} = [ $coords->{$node}->[0] * $params->{mag},
                              $coords->{$node}->[1] * $params->{mag} ] ;
      if (!exists $bbox->{xmax}) {
         $bbox->{xmin} = $new_coords->{$node}->[0] ;
         $bbox->{xmax} = $new_coords->{$node}->[0] ;
         $bbox->{ymin} = $new_coords->{$node}->[1] ;
         $bbox->{ymax} = $new_coords->{$node}->[1] ;
      } else {
         if ($new_coords->{$node}->[0] < $bbox->{xmin}) {
            $bbox->{xmin} = $new_coords->{$node}->[0] ; }
         if ($new_coords->{$node}->[0] > $bbox->{xmax}) {
            $bbox->{xmax} = $new_coords->{$node}->[0] ; }

         if ($new_coords->{$node}->[1] < $bbox->{ymin}) {
            $bbox->{ymin} = $new_coords->{$node}->[1] ; }

         if ($new_coords->{$node}->[1] > $bbox->{ymax}) {
            $bbox->{ymax} = $new_coords->{$node}->[1] ; }
      }
   }
   $coords = $new_coords ;
   $bbox->{xmin} -= $params->{nrad} * 2;
   $bbox->{ymin} -= $params->{nrad} * 2;
   $bbox->{ymax} += $params->{nrad} * 2;
   $bbox->{xmax} += $params->{nrad} * 2;

   $bbox->{ymax} -= $bbox->{ymin} ;
   $bbox->{xmax} -= $bbox->{xmin} ;

   my $mapname = $input->{map_name} ;
   if (!defined $mapname) {$mapname = 'lglpmmap';}
   print {$input->{map_fh}} "<html><body bgcolor='black'>\n" ;
   print {$input->{map_fh}} "<MAP NAME=$mapname>\n";
   foreach my $point (keys %{$coords}) {
      $coords->{$point}->[0] -= $bbox->{xmin} ;
      $coords->{$point}->[1] -= $bbox->{ymin} ;

      my $nodename = $point;
      if (exists $input->{nodenames}->{$point}) {
         $nodename = $input->{nodenames}->{$point} ; }

      if (exists $input->{nodeurls}->{$point}) {
         print {$input->{map_fh}} '<AREA HREF="'.$input->{nodeurls}->{$point}.
               '" shape="circle" title="'.$nodename.'" COORDS="'.
               $coords->{$point}->[0].','.$coords->{$point}->[1].
               ','.$params->{nrad}.'">'."\n";
      } else {
         print {$input->{map_fh}} '<AREA NOHREF="" shape="circle" title="'.$nodename.'" COORDS="'.$coords->{$point}->[0].','.$coords->{$point}->[1].','.$params->{nrad}.'">'."\n";
      }
   }
   print {$input->{map_fh}} "</MAP>\n";
   print {$input->{map_fh}} '<IMG SRC="'.$input->{img_fn}.'" USEMAP="#'.$mapname.'"'."\n" ;
   print {$input->{map_fh}} "</html></body>\n" ;

}


=head2 coords2png()

   Title:	coords2png
   Function:	Creates a png image of a graph, given layout info
   Returns:	GD::image object
   Args:
      $_->{coords}	node coordinate information
      $_->{edges}->{n1}->{n2}	edge list (nhash)
      $_->{ncol}	default node color
      $_->{ncols}->{n}	node-specific color
      $_->{ecol}	default edge color
      $_->{ecols}->{n1}->{n2}	edge-specific color
      $_->{estyle}	default edge style - solid or dashed
      $_->{estyles}->{n1}->{n2}	edge-specific style
      $_->{nshape}	node shape ([f]{triangle|box|circle})
      $_->{nshapes}->{n}	node-specific shape
      $_->{nrad}	node radius
      $_->{ethick}	edge thickness
      $_->{mag}	magnification parameter

=cut

sub coords2png {

   my $input = shift ;

   my $params = {
      ncol => "0 0 0",
      ecol => "0 0 0",
      estyle => "solid",
      nshape => "fcircle",
      nrad => .65,
      ethick => 0.3,
      mag => 7
   } ;

   foreach my $prop (keys %{$input}) {
      $params->{$prop} = $input->{$prop} ; }

   my $coords = $input->{coords};
   my $edges = $input->{edges};

   my $bbox ;
   my $new_coords ;
   foreach my $node (keys %{$coords}) {
      $new_coords->{$node} = [ $coords->{$node}->[0] * $params->{mag},
                              $coords->{$node}->[1] * $params->{mag} ] ;
      if (!exists $bbox->{xmax}) {
         $bbox->{xmin} = $new_coords->{$node}->[0] ;
         $bbox->{xmax} = $new_coords->{$node}->[0] ;
         $bbox->{ymin} = $new_coords->{$node}->[1] ;
         $bbox->{ymax} = $new_coords->{$node}->[1] ;
      } else {
         if ($new_coords->{$node}->[0] < $bbox->{xmin}) {
            $bbox->{xmin} = $new_coords->{$node}->[0] ; }
         if ($new_coords->{$node}->[0] > $bbox->{xmax}) {
            $bbox->{xmax} = $new_coords->{$node}->[0] ; }

         if ($new_coords->{$node}->[1] < $bbox->{ymin}) {
            $bbox->{ymin} = $new_coords->{$node}->[1] ; }

         if ($new_coords->{$node}->[1] > $bbox->{ymax}) {
            $bbox->{ymax} = $new_coords->{$node}->[1] ; }
      }
   }
   my $old_coords = $coords;
   $coords = $new_coords ;
   $bbox->{xmin} -= $params->{nrad} * 2;
   $bbox->{ymin} -= $params->{nrad} * 2;
   $bbox->{ymax} += $params->{nrad} * 2;
   $bbox->{xmax} += $params->{nrad} * 2;

   $bbox->{ymax} -= $bbox->{ymin} ;
   $bbox->{xmax} -= $bbox->{xmin} ;

   foreach my $point (keys %{$coords}) {
      $coords->{$point}->[0] -= $bbox->{xmin} ;
      $coords->{$point}->[1] -= $bbox->{ymin} ;
   }


   my $im = new GD::Image($bbox->{xmax}, $bbox->{ymax}) ;
   my $white = $im->colorAllocate(255,255,255) ;
   $im->transparent($white) ;

   my $gdcolors = {} ;
   foreach my $n1 (keys %{$params->{edges}}) {
      foreach my $n2 (keys %{$params->{edges}->{$n1}}) {
         if ((exists $params->{ecols}->{$n1}->{$n2}) &&
              $params->{ecols}->{$n1}->{$n2} !~ / /) {
            my $col = $params->{ecols}->{$n1}->{$n2} ;
            if ($col eq 'invisible') {
               $params->{ecols}->{$n1}->{$n2} = $col ;
            } else {
               if (exists $color2rgb->{$col}) {
                  $params->{ecols}->{$n1}->{$n2} = $color2rgb->{$col} ;
               } else {
                  $params->{ecols}->{$n1}->{$n2} = $params->{ecol} ;
               }
            }
         } else {
            $params->{ecols}->{$n1}->{$n2} = $params->{ecol};
         }

         if ($params->{ecols}->{$n1}->{$n2} ne 'invisible' &&
             !exists $gdcolors->{$params->{ecols}->{$n1}->{$n2}}) {
                  my @colrgb = split(/ /, $params->{ecols}->{$n1}->{$n2}) ;
                  $gdcolors->{$params->{ecols}->{$n1}->{$n2}} =
                     $im->colorAllocate((255 * $colrgb[0]),
                                          (255 * $colrgb[1]),
                                          (255 * $colrgb[2])) ;
         }

       }
    }

    foreach my $n (keys %{$coords}) {
         if (exists $params->{ncols}->{$n} &&
             $params->{ncols}->{$n} !~ / /) {

            my $col = $params->{ncols}->{$n} ;
            if ($col eq 'invisible') {
               $params->{ncols}->{$n} = 'invisible' ;
            } else {
               if (exists $color2rgb->{$col}) {
                  $params->{ncols}->{$n} = $color2rgb->{$col} ;
               } else {
                  $params->{ncols}->{$n} = $params->{ncol} ; }

               if (!exists $gdcolors->{$params->{ncols}->{$n}}) {
                  my @colrgb = split(/ /, $params->{ncols}->{$n}) ;
                  $gdcolors->{$params->{ncols}->{$n}} = $im->colorAllocate(
                                   (255 * $colrgb[0]),
                                   (255 * $colrgb[1]),
                                   (255 * $colrgb[2])) ;
               }
            }
         } else {
            $params->{ncols}->{$n} = $params->{ncol} ;
         }

         if ($params->{ncols}->{$n} ne 'invisible' &&
             !exists $gdcolors->{$params->{ncols}->{$n}}) {
                  my @colrgb = split(/ /, $params->{ncols}->{$n}) ;
                  $gdcolors->{$params->{ncols}->{$n}} =
                     $im->colorAllocate((255 * $colrgb[0]),
                                          (255 * $colrgb[1]),
                                          (255 * $colrgb[2])) ;
         }
   }


   my $date = `date` ; chomp $date ;

   my $props = 0 ;
   $im->setThickness($params->{ethick}) ;
   foreach my $node1 (keys %{$edges}) {
      foreach my $node2 (keys %{$edges->{$node1}}) {
         my @colrgb ;
         if ($params->{ecols}->{$node1}->{$node2} eq 'invisible') {next;}
         my $curestyle = $params->{estyle} ;
         if (exists $params->{estyles}->{$node1}->{$node2}) {
         $curestyle = $params->{estyles}->{$node1}->{$node2} ; }

         if ($curestyle eq 'solid') {
            $im->line( @{$coords->{$node1}}, @{$coords->{$node2}},
                       $gdcolors->{$params->{ecols}->{$node1}->{$node2}} ) ;
         } elsif ($curestyle eq 'dashed') {
            $im->dashedLine(@{$coords->{$node1}},@{$coords->{$node2}},
                            $gdcolors->{$params->{ecols}->{$node1}->{$node2}}) ;
         }
      }
   }

   foreach my $node (keys %{$coords}) {
      if ($params->{ncols}->{$node} eq 'invisible') { next; }

      my $tshape ;
      if (exists $params->{nshapes}->{$node}) {
         $tshape = $params->{nshapes}->{$node} ;
      } else {
         $tshape = $params->{nshape} ;
      }

      my $fill = 0 ;
      if ($tshape =~ /^f/) {$fill = 1; $tshape =~ s/^f// ;}

      if ($tshape eq 'circle') {
         if (! $fill) {
            $im->arc($coords->{$node}->[0], $coords->{$node}->[1],
                     ($params->{nrad} * 2), ($params->{nrad} * 2),
                     0, 360, $gdcolors->{$params->{ncols}->{$node}}) ;
         } else {
            $im->filledArc($coords->{$node}->[0], $coords->{$node}->[1],
                     ($params->{nrad} * 2), ($params->{nrad} * 2),
                     0, 360, $gdcolors->{$params->{ncols}->{$node}}) ;
         }

      } elsif ($tshape eq 'triangle') {

         my $p1 = [$coords->{$node}->[0],
                   $coords->{$node}->[1] + $params->{nrad}] ;
         my $p2 = [$coords->{$node}->[0] - $params->{nrad},
                   $coords->{$node}->[1] - $params->{nrad}] ;
         my $p3 = [$coords->{$node}->[0] + $params->{nrad},
                   $coords->{$node}->[1] - $params->{nrad}] ;
         my $poly = new GD::Polygon;
         $poly->addPt($p1->[0], $p1->[1]) ;
         $poly->addPt($p2->[0], $p2->[1]) ;
         $poly->addPt($p3->[0], $p3->[1]) ;
         if ($fill)  {
            $im->filledPolygon($poly,$gdcolors->{$params->{ncols}->{$node}}) ;
         } else {
            $im->polygon($poly,$gdcolors->{$params->{ncols}->{$node}}) ;
         }

      } elsif ($tshape eq 'box') {

         my $p1 = [$coords->{$node}->[0] + $params->{nrad},
                   $coords->{$node}->[1] + $params->{nrad}] ;
         my $p2 = [$coords->{$node}->[0] + $params->{nrad},
                   $coords->{$node}->[1] - $params->{nrad}] ;
         my $p3 = [$coords->{$node}->[0] - $params->{nrad},
                   $coords->{$node}->[1] - $params->{nrad}] ;
         my $p4 = [$coords->{$node}->[0] - $params->{nrad},
                   $coords->{$node}->[1] + $params->{nrad}] ;
         my $poly = new GD::Polygon;
         $poly->addPt($p1->[0], $p1->[1]) ;
         $poly->addPt($p2->[0], $p2->[1]) ;
         $poly->addPt($p3->[0], $p3->[1]) ;
         $poly->addPt($p4->[0], $p4->[1]) ;
         if ($fill) {
            $im->filledPolygon($poly,$gdcolors->{$params->{ncols}->{$node}}) ;
         } else {
            $im->polygon($poly,$gdcolors->{$params->{ncols}->{$node}}) ;
         }
      }
   }

   if (exists $input->{out_fn}) {
      open(OUTPNG, ">".$input->{out_fn}) ;
      binmode OUTPNG ;
      print OUTPNG $im->png() ;
      close(OUTPNG) ;
   } else {
      return $im ;
   }

}


=head2 coords2ps()

   Title:	coords2ps
   Function:	Creates a postscript image of a graph, given layout info
   Returns:	nothing
   Args:
      $_->{coords}	node coordinate information
      $_->{edges}->{n1}->{n2}	edge list (nhash)
      $_->{ncol}	default node color
      $_->{ncols}->{n}	node-specific color
      $_->{ecol}	default edge color
      $_->{ecols}->{n1}->{n2}	edge-specific color
      $_->{estyle}	default edge style - solid or dashed
      $_->{estyles}->{n1}->{n2}	edge-specific style
      $_->{nshape}	node shape ([f]{triangle|box|circle})
      $_->{nshapes}->{n}	node-specific shape
      $_->{nrad}	node radius
      $_->{ethick}	edge thickness
      $_->{mag}	magnification parameter

=cut

sub coords2ps {

   my $input = shift ;

   my $params = {
      ncol => "0 0 0",
      ecol => "0 0 0",
      estyle => " ",
      nshape => "fcircle",
      nrad => '.65',
      ethick => '0.3',
      mag => 7
   } ;

   foreach my $prop (keys %{$input}) {
      $params->{$prop} = $input->{$prop} ; }

   if ($params->{estyle} eq 'solid') { $params->{estyle} = " ";}
   if ($params->{estyle} eq 'dashed') { $params->{estyle} = "3 3";}

   if (exists $params->{out_fn}) {
      open($params->{fh}, ">".$params->{out_fn}) ; }
   if (!exists $params->{fh}) {
      open($params->{fh}, ">-") ; }

   my $coords = $input->{coords};
   my $edges = $input->{edges};

   if (exists $params->{ecols}) {
      foreach my $n1 (keys %{$params->{ecols}}) {
       foreach my $n2 (keys %{$params->{ecols}->{$n1}}) {
         my $col = $params->{ecols}->{$n1}->{$n2} ;
         if ($col !~ / /) {
            if (exists $color2rgb->{$col}) {
               $params->{ecols}->{$n1}->{$n2} = $color2rgb->{$col} ;
            } else {
               $params->{ecols}->{$n1}->{$n2} = $params->{ecol} ;
            }
         }
       }
      }
   }

   if (exists $params->{ncols}) {
      foreach my $n (keys %{$params->{ncols}}) {
         my $col = $params->{ncols}->{$n} ;
         if ($col !~ / /) {
            if (exists $color2rgb->{$col}) {
               $params->{ncols}->{$n} = $color2rgb->{$col} ;
            } else {
               $params->{ncols}->{$n} = $params->{ncol} ;
            }
         }
      }
   }

   my $bbox ;
   my $new_coords ;
   foreach my $node (keys %{$coords}) {
      $new_coords->{$node} = [ $coords->{$node}->[0] * $params->{mag},
                              $coords->{$node}->[1] * $params->{mag} ] ;
      if (!exists $bbox->{xmax}) {
         $bbox->{xmin} = $new_coords->{$node}->[0] ;
         $bbox->{xmax} = $new_coords->{$node}->[0] ;
         $bbox->{ymin} = $new_coords->{$node}->[1] ;
         $bbox->{ymax} = $new_coords->{$node}->[1] ;
      } else {
         if ($new_coords->{$node}->[0] < $bbox->{xmin}) {
            $bbox->{xmin} = $new_coords->{$node}->[0] ; }
         if ($new_coords->{$node}->[0] > $bbox->{xmax}) {
            $bbox->{xmax} = $new_coords->{$node}->[0] ; }

         if ($new_coords->{$node}->[1] < $bbox->{ymin}) {
            $bbox->{ymin} = $new_coords->{$node}->[1] ; }

         if ($new_coords->{$node}->[1] > $bbox->{ymax}) {
            $bbox->{ymax} = $new_coords->{$node}->[1] ; }
      }
   }
   my $old_coords = $coords;
   $coords = $new_coords ;
   $bbox->{xmin} -= $params->{nrad} * 2;
   $bbox->{ymin} -= $params->{nrad} * 2;
   $bbox->{ymax} += $params->{nrad} * 2;
   $bbox->{xmax} += $params->{nrad} * 2;

   my $date = `date` ; chomp $date ;
   print {$params->{fh}} '%!PS-Adobe-3.0 EPSF-3.0'."\n";
   print {$params->{fh}} '%%Creator: graph2ps.pl'."\n";
   print {$params->{fh}} '%%Title: graph layout'."\n";
   print {$params->{fh}} '%%CreationDate: '.$date."\n";
   print {$params->{fh}} '%%DocumentData: Clean7Bit'."\n";
   print {$params->{fh}} '%%Origin: [0 0]'."\n";
   print {$params->{fh}} '%%BoundingBox: '.join(' ',($bbox->{xmin},
         $bbox->{ymin}, $bbox->{xmax}, $bbox->{ymax}))."\n";
   print {$params->{fh}} '%%LanguageLevel: 2'."\n";
   print {$params->{fh}} '%%Pages: 1'."\n";
   print {$params->{fh}} '%%Page: 1 1'."\n";
   print {$params->{fh}} '%%EOF'."\n";

   if (exists $params->{label_fl}) {
      print {$params->{fh}} "/Times-Roman findfont\n" ;
      print {$params->{fh}} "10 scalefont\n" ;
      print {$params->{fh}} "setfont\n" ;
   }

   if (exists $params->{estyles}) {
      foreach my $n1 (keys %{$params->{estyles}}) {
         foreach my $n2 (keys %{$params->{estyles}->{$n1}}) {
            if ($params->{estyles}->{$n1}->{$n2} eq 'solid') {
               $params->{estyles}->{$n1}->{$n2} = " ";
            } elsif ($params->{estyles}->{$n1}->{$n2} eq 'dashed') {
               $params->{estyles}->{$n1}->{$n2} = "3 3";
            }
         }
      }
   }

   my $props = 0 ;
   foreach my $node1 (keys %{$edges}) {
      foreach my $node2 (keys %{$edges->{$node1}}) {
         if (exists $params->{estyles}->{$node1}->{$node2}) {
            print {$params->{fh}} "[$params->{estyles}->{$node1}->{$node2}]".
                  " 0 setdash\n" ;
         } else {
            print {$params->{fh}} "[$params->{estyle}] 0 setdash\n" ;
         }
         if (exists $params->{ecols}->{$node1}->{$node2}) {
            print {$params->{fh}} $params->{ecols}->{$node1}->{$node2}.
                  " setrgbcolor\n" ;
         } else {
            print {$params->{fh}} "$params->{ecol} setrgbcolor\n" ;
         }
         print {$params->{fh}} "$params->{ethick} setlinewidth\n" ;
         print {$params->{fh}} join(" ", @{$coords->{$node1}}, "moveto")."\n" ;
         print {$params->{fh}} join(" ", @{$coords->{$node2}}, "lineto")."\n" ;
         print {$params->{fh}} "stroke\n" ;
      }
   }
   print {$params->{fh}} "[] 0 setdash\n" ;

   foreach my $node (keys %{$coords}) {
      if (exists $params->{ncols}->{$node}) {
         print {$params->{fh}} $params->{ncols}->{$node}." setrgbcolor\n" ;
      } else {
         print {$params->{fh}} "$params->{ncol} setrgbcolor\n" ;
      }

      my $tshape ;
      if (exists $params->{nshapes}->{$node}) {
         $tshape = $params->{nshapes}->{$node} ;
      } else {
         $tshape = $params->{nshape} ;
      }

      my $fill = 0 ;
      if ($tshape =~ /^f/) {$fill = 1; $tshape =~ s/^f// ;}

      if (exists $params->{label_fl}) {
         print {$params->{fh}} "newpath\n" ;
         print {$params->{fh}} "$coords->{$node}->[0] ".
                                 "$coords->{$node}->[1] moveto\n" ;
         print {$params->{fh}} "($node) show\n" ;
      }

      if ($tshape eq 'circle') {
         print {$params->{fh}} join(" ", ($coords->{$node}->[0],
                                          $coords->{$node}->[1],
                                          $params->{nrad}, 0, 360,
                                          "arc"))."\n" ;
      } elsif ($tshape eq 'triangle') {
         my $p1 = [$coords->{$node}->[0],
                   $coords->{$node}->[1] + $params->{nrad}] ;
         my $p2 = [$coords->{$node}->[0] - $params->{nrad},
                   $coords->{$node}->[1] - $params->{nrad}] ;
         my $p3 = [$coords->{$node}->[0] + $params->{nrad},
                   $coords->{$node}->[1] - $params->{nrad}] ;
         print {$params->{fh}} "newpath\n" ;
         print {$params->{fh}} join(" ", ($p1->[0], $p1->[1], "moveto"))."\n" ;
         print {$params->{fh}} join(" ", ($p2->[0], $p2->[1], "lineto"))."\n" ;
         print {$params->{fh}} join(" ", ($p3->[0], $p3->[1], "lineto"))."\n" ;
      } elsif ($tshape eq 'box') {
         my $p1 = [$coords->{$node}->[0] + $params->{nrad},
                   $coords->{$node}->[1] + $params->{nrad}] ;
         my $p2 = [$coords->{$node}->[0] + $params->{nrad},
                   $coords->{$node}->[1] - $params->{nrad}] ;
         my $p3 = [$coords->{$node}->[0] - $params->{nrad},
                   $coords->{$node}->[1] - $params->{nrad}] ;
         my $p4 = [$coords->{$node}->[0] - $params->{nrad},
                   $coords->{$node}->[1] + $params->{nrad}] ;
         print {$params->{fh}} "newpath\n" ;
         print {$params->{fh}} join(" ", ($p1->[0], $p1->[1], "moveto"))."\n" ;
         print {$params->{fh}} join(" ", ($p2->[0], $p2->[1], "lineto"))."\n" ;
         print {$params->{fh}} join(" ", ($p3->[0], $p3->[1], "lineto"))."\n" ;
         print {$params->{fh}} join(" ", ($p4->[0], $p4->[1], "lineto"))."\n" ;
      }
      if ($fill) { print {$params->{fh}} "fill\n" ;}
      print {$params->{fh}} "closepath\nstroke\n" ;
   }

}


=head2 get_color2rgb()

   Title:       get_color2rgb()
   Function:    Returns a hash pointing from color name to rgb (0-1) values
   Args:        Nothing
   Returns:     $_->{color} = "$r $g $b" ;

=cut

sub get_color2rgb {

   my %a = %{$color2rgb} ;
   return \%a ;

}


1 ;
