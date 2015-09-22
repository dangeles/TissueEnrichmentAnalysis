#!/usr/bin/perl

use strict;
open(my $fh, '>', 'output.txt');
my (@files) = <../input/genesets_golden/*>;
foreach my $file (@files) {
  print $fh qq($file\n);
  my $pyreturn = `python hypergeometricTests.py $file 2>&1`;
  print $fh qq($pyreturn\n);
  print $fh qq(#####\n);
  print qq($file\n);
}

close $fh;
