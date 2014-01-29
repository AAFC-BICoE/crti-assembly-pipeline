#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Text::CSV;
use Bio::DB::Taxonomy;

# Take parsed blast output file, line by line
# Get the species name
# get taxonomy output
# Add this to the front of the blast line
# Add the full header
# Write CSV and/or tab-delimited output

my @parsed_blast_headers = qw(Query_Name Query_Description Hit_Name Hit_Description 
    HSP_Num Evalue Score Bits Pct_Identical Pct_Conserved Q_Start Q_End Q_Match 
    Q_Length Q_Strand S_Start S_End S_Match S_Length S_Strand);
my @ranks = qw /superkingdom kingdom phylum class order/;
my @taxa_headers = (@ranks, "ID", "Rank", "Scientific_Name");
my @output_headers = (@taxa_headers, "Target_Rank_Distance", @parsed_blast_headers);

my $csv = Text::CSV->new ({ binary => 1, eol => "\n" });
my $tsv = Text::CSV->new ({ binary => 1, sep_char => "\t", eol => "\n" });
$tsv->column_names(@parsed_blast_headers);

my $parsed_blast = '';
my $taxa_tsv = '';
my $taxa_csv = '';
my $target_taxon = '';
GetOptions ('parsed_blast|p=s' => \$parsed_blast,
            'taxa_csv|c=s' => \$taxa_csv,
            'target_taxon|t=s' => \$target_taxon,);
$parsed_blast and $taxa_csv or die "Usage: $0 -p <parsed blast file> -c <taxa csv file>\n";

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
        push (@$ranks_found, $node->rank()); #Seems this is always species...
        push (@$ranks_found, $node->scientific_name());
    } else {
        push (@$ranks_found, "Name \"$species\" not found\n");
    }
    while (scalar @$ranks_found < scalar @taxa_headers) {
        push (@$ranks_found, '');
    }
    return $ranks_found;
}

# Calculate the rank distance between this blast hit and the target.
# This function is likely buggy
sub calc_target_distance
{
    my $target_ranks = shift;
    my $hit_ranks = shift;
    my $h = scalar @$hit_ranks;
    my $t = scalar @$target_ranks;
    my $dist = 0;
    my $divergent = 0;
    if ($h == $t) {
        my $i = 0;
        while ($i < $h) {
            my $hr = $hit_ranks->[$i];
            my $tr = $target_ranks->[$i];
            if ($hr ne $tr) {
                $divergent = 1;
                $dist += 2;
            } elsif ($divergent) {
                # assert hr == tr == '' or '-'
                # assuming rank names etc are always unique across different lineages
                $dist += 2;
            }
            $i++;
        }
    }
    return $dist;
}


open (my $fpb, '<:encoding(utf8)', $parsed_blast) or die "Error: couldn't open file $parsed_blast\n";
open (my $ftax, '>:encoding(utf8)', $taxa_csv) or die "Error: couldn't open file $taxa_csv\n";
my $dbh = create_taxonomy_object;

my $target_ranks = [];
if ($target_taxon) {
    $target_ranks = get_taxon_rank ($dbh, $target_taxon);
}

$csv->print($ftax, \@output_headers);
#while (my $row = $tsv->getline_hr($fpb)) {
#    my $hit = $row->{"Hit_Description"};
while (my $line = <$fpb>) {
    my @fields = split (/\t/, $line);
    my $ranks_found = ['','','','','','','',''];
    my $hit = $fields[3];
    my $target_distance = "N/A";
    if ($hit) {
        $hit =~ s/^PREDICTED:\s+//;
        my @words = split (/\s+/, $hit);
        if (scalar @words >= 2) {
            my $species_name = join(' ', @words[0..1]); # Get first two words as species name
            $ranks_found = get_taxon_rank ($dbh, $species_name);
            
            if ($target_taxon) {
                $target_distance = calc_target_distance ($target_ranks, $ranks_found);
            }
        }
    }
    my $combined_fields = [@$ranks_found, $target_distance, @fields];
    $csv->print($ftax, $combined_fields);
}

close ($fpb);
close ($ftax);
    