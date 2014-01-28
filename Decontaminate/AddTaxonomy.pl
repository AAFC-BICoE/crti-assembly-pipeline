#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Text::CSV;

# Take parsed blast output file, line by line
# Get the species name
# get taxonomy output
# Add this to the front of the blast line
# Add the full header
# Write CSV and/or tab-delimited output

my $parsed_blast = '';
my $taxa_tsv = '';
my $taxa_csv = '';

my @parsed_blast_headers = qw(Query_Name Query_Description Hit_Name Hit_Description 
    HSP_Num Evalue Score Bits Pct_Identical Pct_Conserved Q_Start Q_End Q_Match 
    Q_Length Q_Strand S_Start S_End S_Match S_Length S_Strand);
my @taxa_headers = qw(Superkingdom Kingdom Phylum);
my @output_headers = ();

my $csv = Text::CSV->new ({ binary => 1, eol => "\n" });
my $tsv = Text::CSV->new ({ binary => 1, sep_char => "\t", eol => "\n" });
$tsv->column_names(@parsed_blast_headers);

GetOptions ('parsed_blast|p=s' => \$parsed_blast,
            'taxa_csv|c=s' => \$taxa_csv,);
$parsed_blast and $taxa_csv or die "Usage: $0 -p <parsed blast file> -c <taxa csv file>\n";

open (my $fpb, '<:encoding(utf8)', $parsed_blast) or die "Error: couldn't open file $parsed_blast\n";
open (my $ftax, '>:encoding(utf8)', $taxa_csv) or die "Error: couldn't open file $taxa_csv\n";

sub create_taxonomy_object
{
    my $taxa_dir = "/isilon/biodiversity/reference/ncbi_taxonomy";
    my $nodes = $taxa_dir . "/nodes.dmp";
    my $names = $taxa_dir . "/names.dmp";
    my $dbh = new Bio::DB::Taxonomy (-source => "flatfile",
				 -directory => $taxa_dir,
				 -nodesfile => $nodes,
				 -namesfile => $names,
				 # -force => $options->{'reindex'},
				);
	return $dbh;
}

sub get_taxon_rank
{
    my $dbh = shift;
    my $species = shift;
    my @ranks = qw /superkingdom kingdom phylum class order/;
    my $node = $dbh->get_taxon (-name => $species);
    my $ranks_found = [];
    if (defined $node) {
        if (scalar @ranks > 0) {
            my $rank_header = join "\t", @ranks;
            #print STDERR $rank_header."\n";
            my $ref = $node;
            my $ranks = {};
            while (defined $ref) {
                $ranks->{$ref->rank()} = $ref->scientific_name();
                $ref = $ref->ancestor();
            }
            foreach my $rank (@ranks) {
                if (defined $ranks->{$rank}) {
                    push (@$ranks_found, $ranks->{$rank});
                } else {
                    push (@$ranks_found, "-");
                }
            }
        }
        push (@$ranks_found, $node->id());
        push (@$ranks_found, $node->rank());
        push (@$ranks_found, $node->scientific_name());
    } else {
        push (@$ranks_found, "Name \"$species\" not found\n");
    }
    return $ranks_found;
}

my $dbh = create_taxonomy_object;

while (my $row = $tsv->getline_hr($fpb)) {
    my $hit = $row->{"Hit_Description"};
    $hit =~ s/^PREDICTED:\s+//;
    my $species_name = join(' ', (split (/\s+/, $hit))[0..1]); # Get first two words as species name
    my $ranks_found = get_taxon_rank ($dbh, $species_name);
    print join ("\t", @$ranks_found) . "\n";
}
    
    