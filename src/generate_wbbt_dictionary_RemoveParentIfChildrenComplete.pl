#!/usr/bin/perl


# generate dictionary of genes to anatomy terms based on
# http://wormbase.caltech.edu:8080/wormbase/manual/geneenrichment/davidahypergeometrictest
## variant that eliminates anatomy terms that have all children already in the list to remove some redundancy

use CGI;
use strict;
use LWP::Simple;
use JSON;

my $solr_url = 'http://wobr.caltech.edu:8082/solr/anatomy/';

my $json = JSON->new->allow_nonref;

my $url = $solr_url . "select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=0&q=document_category:bioentity&facet=true&facet.field=regulates_closure&facet.limit=-1&facet.mincount=100&facet.sort=count&fq=source:%22WB%22&fq=-qualifier:%22not%22";

my $page_data = get $url;

my $perl_scalar = $json->decode( $page_data );
my %jsonHash = %$perl_scalar;

## Fetch anatomy terms that have sufficient number (specified in solr query, $url) of annotated genes
my $arrRef = $jsonHash{"facet_counts"}{"facet_fields"}{"regulates_closure"} ;
my @array = @$arrRef;
my %wbbt;
while (@array) {
  my $wbbt = shift @array;
  my $count = shift @array;
  $wbbt{$wbbt} = $count;
}

## Remove terms from list if any of its descendants is in the list
=comment out
foreach my $wbbt (sort { $wbbt{$a} <=> $wbbt{$b} } keys %wbbt) { 
  next unless ($wbbt{$wbbt});				# terms removed while looping, so skip those already removed
#  print qq($wbbt\t$wbbt{$wbbt}\n);
  my $url = $solr_url . "select?qt=standard&fl=regulates_closure&version=2.2&wt=json&indent=on&rows=1&q=id:%22" . $wbbt . "%22&fq=document_category:%22ontology_class%22";
  my $page_data = get $url;
  my $perl_scalar = $json->decode( $page_data );
  my %jsonHash = %$perl_scalar;
  my $arrRef = $jsonHash{'response'}{'docs'}[0]{'regulates_closure'};
  my @array = @$arrRef;
  foreach my $other (@array) {
    if ( ($other ne $wbbt) && ($wbbt{$other}) ) { delete $wbbt{$other}; } }	# if other term is in list, remove from list
} 
=cut

## Fetch genes annotated to each anatomy term
my %genes;
foreach my $wbbt (sort keys %wbbt) {
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
foreach my $wbbt (sort keys %wbbt) {
  push @wbbts_header, qq($wbbtName{$wbbt}($wbbt));
}
my $wbbts_header = join",", @wbbts_header;
# my $wbbts_header = join",", sort keys %wbbt;
print qq(wbid,$wbbts_header\n);
foreach my $gene (sort keys %genes) {
  my @out;
  foreach my $wbbt (sort keys %wbbt) {
    if ($genes{$gene}{$wbbt}) { push @out, '1'; } else { push @out, '0'; } }
  my $out = join",", @out;
  print qq($gene,$out\n);
} # foreach my $gene (sort keys %genes)


# print qq($reg_clos\n);


__END__

http://131.215.12.204:8080/solr/anatomy/select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=0&q=document_category:bioentity&facet=true&facet.field=regulates_closure&facet.limit=-1&facet.mincount=100&facet.sort=count&fq=source:%22WB%22&fq=-qualifier:%22not%22

http://131.215.12.204:8080/solr/anatomy/select?qt=standard&fl=regulates_closure&version=2.2&wt=json&indent=on&rows=1&q=id:%22WBbt:0005373%22&fq=document_category:%22ontology_class%22

http://131.215.12.204:8080/solr/anatomy/select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=10000&q=document_category:bioentity&fq=source:%22WB%22&fq=-qualifier:%22not%22&fq=regulates_closure:%22WBbt:0005373%22
