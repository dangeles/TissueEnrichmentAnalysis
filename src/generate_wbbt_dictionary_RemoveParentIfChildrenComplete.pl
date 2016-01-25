#!/usr/bin/perl


# generate dictionary of genes to anatomy terms based on
# http://wormbase.caltech.edu:8080/wormbase/manual/geneenrichment/davidahypergeometrictest
# eliminates anatomy terms that have all children (with some annotation) already in the list to remove some redundancy




# we're generating KEEP and using that.  We have GOOD and we're never
# adding nor removing from it.  DISCARD we're keeping just in case we
# want to look at it later.
# Terms only go in KEEP if it either has a transitive term not in GOOD,
# or no transitive terms at all.
# Reading the regulates_closure counts into a %regulatesCount hash.
# Then when looping through GOOD, if an edge has a sub with a count of
# zero, it gets ignored (the edge gets ignored as if it was never
# there).

use CGI;
use strict;
use LWP::Simple;
use JSON;

my $solr_url = 'http://wobr.caltech.edu:8082/solr/anatomy/';

my $json = JSON->new->allow_nonref;

my $annot_cutoff = 100;
if ($ARGV[0]) {$annot_cutoff = $ARGV[0]}
my $url = $solr_url . "select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=0&q=document_category:bioentity&facet=true&facet.field=regulates_closure&facet.limit=-1&facet.mincount=$annot_cutoff&facet.sort=count&fq=source:%22WB%22&fq=-qualifier:%22not%22";




my $page_data = get $url;

my $perl_scalar = $json->decode( $page_data );
my %jsonHash = %$perl_scalar;

## Fetch anatomy terms that have sufficient number (specified in solr query above, $url) of annotated genes
my $arrRef = $jsonHash{"facet_counts"}{"facet_fields"}{"regulates_closure"} ;
my @array = @$arrRef;
my %wbbt;
my %good;
my %keep;
my %discard;
while (@array) {
  my $wbbt = shift @array;
  my $count = shift @array;
  $wbbt{$wbbt} = $count;
  $good{$wbbt} = $count;
}

## Fetch 
my %regulatesCount;
my $urlRegulates = $solr_url . "select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=0&q=document_category:bioentity&facet=true&facet.field=regulates_closure&facet.limit=-1&facet.sort=count&fq=source:%22WB%22&fq=-qualifier:%22not%22";
my $pageRegulates_data = get $urlRegulates;
my $perlRegulates_scalar = $json->decode( $pageRegulates_data );
my %jsonHashRegulates = %$perlRegulates_scalar;
my $arrRefRegulates = $jsonHashRegulates{'facet_counts'}{'facet_fields'}{'regulates_closure'};
my @arrayRegulates = @$arrRefRegulates;


while (scalar @arrayRegulates > 0) {
   my $id    = shift @arrayRegulates;
   my $count = shift @arrayRegulates;
   $regulatesCount{$id} = $count;
}

foreach my $wbbt (sort { $wbbt{$a} <=> $wbbt{$b} } keys %wbbt) { 
  my $url = $solr_url . "select?qt=standard&fl=topology_graph_json&version=2.2&wt=json&indent=on&rows=1&q=id:%22" . $wbbt . "%22&fq=document_category:%22ontology_class%22";
  my $page_data = get $url;
  my $perl_scalar = $json->decode( $page_data );
  my %jsonHash = %$perl_scalar;
  my $topo_data = $json->decode( $jsonHash{"response"}{"docs"}[0]{"topology_graph_json"} );
  my %topo = %$topo_data;
  my (@edges) = @{ $topo{"edges"} };
  my $hasTransChild = 0;
  my $hasTransChildNotGood = 0;
  for my $index (0 .. @edges) {
    my ($sub, $obj, $pred) = ('', '', '');
    my $isTrans = 0;
    if ($edges[$index]{'sub'}) {  $sub  = $edges[$index]{'sub'};  }
    if ($edges[$index]{'obj'}) {  $obj  = $edges[$index]{'obj'};  }
    if ($edges[$index]{'pred'}) { $pred = $edges[$index]{'pred'}; }
    next unless ($regulatesCount{$sub});
    if ( ($pred eq 'is_a') || ($pred eq 'part_of') || ($pred eq 'occurs_in') || ($pred =~ m/.*regulates/) ) { $isTrans++; }
    if ( ($isTrans) && ($obj eq $wbbt) ) { 
      unless ($good{$sub}) { $hasTransChildNotGood++; }
      $hasTransChild++; }
  }
  if ($hasTransChildNotGood) { $keep{$wbbt}++;    }
    elsif ($hasTransChild) {   $discard{$wbbt}++; }
    else {                     $keep{$wbbt}++;    }
}

# =for testing
# foreach my $wbbt (sort keys %discard) { print qq(DISCARD $wbbt DISC\n); }
# foreach my $wbbt (sort keys %keep) { print qq(KEEP $wbbt KEEP\n); }
# foreach my $wbbt (sort keys %wbbt) { print qq(GOOD $wbbt GOOD\n); }
# =cut


## Fetch genes annotated to each KEEP anatomy term
my %genes;
foreach my $wbbt (sort keys %keep) {
  my $url = $solr_url . "select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=10000&q=document_category:bioentity&fq=source:%22WB%22&fq=-qualifier:%22not%22&fq=regulates_closure:%22" . $wbbt . "%22";
  my $page_data = get $url;
  my $perl_scalar = $json->decode( $page_data );
  my %jsonHash = %$perl_scalar;
  my $arrRef = $jsonHash{'response'}{'docs'};
  my @array = @$arrRef;
  foreach my $hashRef (@array) {
    my %hash = %$hashRef;
    my $gene = $hash{'id'};
    $gene =~ s/^WB://;
#     print "$wbbt OTHER $hash{'id'}\n";
    $genes{$gene}{$wbbt}++; } }

## Fetch human-readable names for anatomy terms
my %wbbtName;
my $url = $solr_url . "select?qt=standard&fl=id,annotation_class_label&version=2.2&wt=json&indent=on&rows=100000&q=id:*&fq=document_category:ontology_class&fq=-is_obsolete:true";
my $page_data = get $url;
my $perl_scalar = $json->decode( $page_data );
my %jsonHash = %$perl_scalar;
my $arrRef = $jsonHash{'response'}{'docs'};
my @array = @$arrRef;
foreach my $hashRef (@array) {
  my %hash = %$hashRef;
  my $id = $hash{'id'};
  my $name = $hash{'annotation_class_label'};
  $name =~ tr/,/comma/;
  $wbbtName{$id} = $name; }


## Print dictionary in csv
my @wbbts_header;
foreach my $wbbt (sort keys %keep) {
  push @wbbts_header, qq($wbbtName{$wbbt}($wbbt));
}
my $wbbts_header = join",", @wbbts_header;
# my $wbbts_header = join",", sort keys %wbbt;
print qq(wbid,$wbbts_header\n);
foreach my $gene (sort keys %genes) {
  my @out;
  foreach my $wbbt (sort keys %keep) {
    if ($genes{$gene}{$wbbt}) { push @out, '1'; } else { push @out, '0'; } }
  my $out = join",", @out;
  print qq($gene,$out\n);
} # foreach my $gene (sort keys %genes)


# print qq($reg_clos\n);


__END__

http://131.215.12.204:8080/solr/anatomy/select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=0&q=document_category:bioentity&facet=true&facet.field=regulates_closure&facet.limit=-1&facet.mincount=100&facet.sort=count&fq=source:%22WB%22&fq=-qualifier:%22not%22

http://131.215.12.204:8080/solr/anatomy/select?qt=standard&fl=regulates_closure&version=2.2&wt=json&indent=on&rows=1&q=id:%22WBbt:0005373%22&fq=document_category:%22ontology_class%22

http://131.215.12.204:8080/solr/anatomy/select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=10000&q=document_category:bioentity&fq=source:%22WB%22&fq=-qualifier:%22not%22&fq=regulates_closure:%22WBbt:0005373%22
