#!/usr/bin/perl

# given input TSV file of gene + anatomy term + something, group by something and anatomy term all genes with both in common and output into a file for each pair.  output to genesets/<file>  2015 09 10

use strict;

my $outdir = '/home/raymond/local/src/git/tissue_enrichment_tool_hypergeometric_test/input/genesets/';
my $infile = '/home/raymond/local/src/git/tissue_enrichment_tool_hypergeometric_test/input/WS250_ExpressionClusterWithAnatomy.TSV';

my %hash;
open (IN, "<$infile") or die "Cannot open $infile : $!";
while (my $line = <IN>) {
  chomp $line;
  $line =~ s/"//g;
  my ($gene, $wbbt, $third) = split/\t/, $line;
  $third =~ s/:/_/g;
  $hash{$third}{$wbbt}{$gene}++;
} # while (my $line = <IN>)
close (IN) or die "Cannot close $infile : $!";

foreach my $third (sort keys %hash) {
  foreach my $wbbt (sort keys %{ $hash{$third} }) {
    my (@genes) = sort keys %{ $hash{$third}{$wbbt} };
    my $count   = scalar @genes;
#     my $outfile = $outdir . $wbbt . '_' . $count . '_' . $third;
    my $outfile = $outdir . $third . '_' . $wbbt . '_' . $count;
    open (OUT, ">$outfile") or die "Cannot create $outfile : $!";
    print OUT "gene,reads\n";
    foreach my $gene (@genes) { print OUT "$gene\n"; }
    close (OUT) or die "Cannot close $outfile : $!";
  } # foreach my $wbbt (sort keys %{ $hash{$third} })
} # foreach my $third (sort keys %hash)

__END__

"WBGene00006528"	"WBbt:0003666"	"WBPaper00045974:NSM_enriched_totalRNA_RNAseq"
"WBGene00001648"	"WBbt:0003666"	"WBPaper00045974:NSM_enriched_totalRNA_RNAseq"
"WBGene00001130"	"WBbt:0003666"	"WBPaper00045974:NSM_enriched_totalRNA_RNAseq"
