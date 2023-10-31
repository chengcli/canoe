## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2020 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## Adapted to Athena++/Atmosphere by Cheng Li
## ---------------------------------------------------------------------

if ($#ARGV != 2) {
  print "\nUsage: make_step.pl step root_dir id\n";
  exit;
}

$step=$ARGV[0];
$step_underscore=$step;
$step_underscore=~ s/-/_/;

$root_dir=$ARGV[1];
$count=$ARGV[2];

print
"/**
\@page $step_underscore Example #$count : The $step Problem
";

#open BF, "$root_dir/examples/$step/builds-on"
#    or die "Can't open builds-on file $root_dir/examples/$step/builds-on";
#my $buildson = <BF>;
#close BF;
#chop $buildson;
$buildson = "";

# At the very top, print which other programs this one builds on. The
# filter script will replace occurrences of step-XX by the appropriate
# links.
if ($buildson ne "")
{
    $buildson =~ s/ /, /g;
    print "This tutorial depends on $buildson.\n\n";
}

# then show the table of contents
print
"\@htmlonly
<div class=\"contents\">
<div class=\"toc\">
<p><b>Table of contents</b></p>
<ul>
  <li><a href=\"#Intro\" class=bold>Introduction</a></li>
";

system $^X, "$root_dir/doc/scripts/intro2toc", "$root_dir/examples/$step/intro.dox";

print "  <li><a href=\"#CommProg\" class=bold>Commented program</a></li>\n";

my $file_extension;

my $pgen = substr $step, index($step, ".")+1;

if (-f "$root_dir/drum/pgen/$pgen.cpp")
{
  $file_extension = cpp;
}

if (-f "$root_dir/drum/pgen/$pgen.cu")
{
  $file_extension = cu;
}

system $^X, "$root_dir/doc/scripts/program2toc", "$root_dir/drum/pgen/$pgen.$file_extension";

print
"  <li><a href=\"#Results\" class=bold>Results</a></li>
";

system $^X, "$root_dir/doc/scripts/intro2toc", "$root_dir/examples/$step/results.dox";

print
"  <li><a href=\"#PlainProg\" class=bold>Plain program</a></li>
</ul></div></div>
\@endhtmlonly
";

system $^X, "$root_dir/doc/scripts/create_anchors", "$root_dir/examples/$step/intro.dox";


# Start the commented program by writing two empty lines. We have had
# cases where the end of the intro.dox was missing a newline, and in
# that case doxygen might get confused about what is being added here
# to the end of an existing line. So add a newline.
#
# But then we also had a situation where doxygen was confused about a
# line starting with an anchor (see #9357). It's not clear what the
# cause is, but making sure that there is an empty line in between
# solves the problem -- so a second newline character.
print " *\n";
print " *\n";
print " * <a name=\"CommProg\"></a>\n";
print " * <h1> The commented program</h1>\n";

system $^X, "$root_dir/doc/scripts/program2doxygen", "$root_dir/drum/pgen/$pgen.$file_extension";

system $^X, "$root_dir/doc/scripts/create_anchors", "$root_dir/examples/$step/results.dox";


# Move to the stripped, plain program. The same principle as above
# applies for newlines.
print " *\n";
print " *\n";
print
"<a name=\"PlainProg\"></a>
<h1> The plain program</h1>
\@include \"${pgen}_plain.txt\"
*/
";
