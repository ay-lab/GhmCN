#!/usr/bin/env perl
use strict;
use Getopt::Long;

my ($reference,$query,%Reference,$out);

# >>> GetOpts
GetOptions("r=s"=>\$reference,
"q=s"=>\$query,
"o=s"=>\$out) or die("Error in command line arguments");

# >>> Build Reference Hash
open(FILE,"$reference");
 while(<FILE>){
  if(/(chr.*)\t(\d*)\t(\d*)\t(\d*)/){ #$2=C; $3=CT
   $Reference{$4}="$1\t$2";
  }
 }
close(FILE);

# >>> Process query with Hashed Reference:
open(FILE,"$query");
 open(OUT,"> $out");
 while(<FILE>){
  if(/(\d*)\t(\d*)\t(.*)/){ #$2=C; $3=CT
   print OUT "$Reference{$1}\t$Reference{$2}\t$3\n";
  }
 }
 close(OUT);
close(FILE);

