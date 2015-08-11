#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Text::CSV;
use Bio::LITE::Taxonomy::NCBI;

my $options = {};
GetOptions($options, 
    'blast_tsv|b=s',
    'taxonomy_csv|t=s',
    );

$options->{blast_tsv} and $options->{taxonomy_csv} or die "Usage: $0 -b <blast tsv> -t <taxonomy csv>\n";

my $csv = Text::CSV->new ({ binary => 1, eol => "\n" });

# blast tab output headers (-outfmt 6)
my @blast_headers = qw(query_id subject_id pct_identity alignment_length mismatches 
    gap_opens q_start q_end s_start s_end evalue bit_score);
# taxonomic ranks to report
my @ranks = qw /superkingdom kingdom phylum class order/;
# headers to include in output taxonomy file.
my @taxa_headers = (@ranks, "ID", "Rank", "Scientific_Name");
my @output_headers = (@taxa_headers, @blast_headers);

my $taxDB = Bio::LITE::Taxonomy::NCBI->new (
    db => "NCBI",
    names => "/isilon/biodiversity/reference/ncbi_taxonomy/names.dmp",
    nodes => "/isilon/biodiversity/reference/ncbi_taxonomy/nodes.dmp",
    dict => "gi_taxid_nucl.bin",
    );

#my @tax = $taxDB->get_taxonomy_from_gi(690029403);
#print join("\t", @tax) . "\n";

# Open the blast output file. Assume -outfmt 6 was used.
open (my $fb, '<:encoding(utf8)', $options->{blast_tsv}) or die "Error: couldn't open file " . $options->{blast_tsv} . "\n";
open (my $ftax, '>:encoding(utf8)', $options->{taxonomy_csv}) or die "Error: couldn't open file" . $options->{taxonomy_csv} . "\n";
$csv->print($ftax, \@output_headers);

while (my $line = <$fb>) {
    my @fields = split (/\t/, $line);
    my $ranks_found = ['','','','','','','',''];
    my $subject_id = $fields[1];
    my $gi;
    if ($subject_id =~ /gi\|([0-9]+)\|/) {
	$gi = $1;
	my @tax = $taxDB->get_taxonomy_from_gi($gi);
	my $combined_fields = [@tax, @fields];
	$csv->print($ftax, $combined_fields);
    } else {
	print "Warning: could not parse gi at line $.\n";
    }
}

close ($fb);
close ($ftax);
