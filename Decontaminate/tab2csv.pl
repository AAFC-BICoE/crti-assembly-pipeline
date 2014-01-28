use strict;
use warnings;
use autodie;
use Text::CSV;
use Getopt::Long;

my ($tabfile, $csvfile) = ('','');
GetOptions('tabfile|t=s' => \$tabfile,
        'csvfile|c=s' => \$csvfile,);
unless ($tabfile and $csvfile) {
    die "Usage: $0 -t <tab file> -c <csvfile>\n";
}

my $csv = Text::CSV->new ({ binary => 1 });
my $tsv = Text::CSV->new ({ binary => 1, sep_char => "\t", eol => "\n" });

open (my $infh,  "<:encoding(utf8)", $csvfile) or die "Error: couldn't open CSV file $csvfile\n";
open (my $outfh, ">:encoding(utf8)", $tabfile) or die "Error: coulnd't open tab file $tabfile\n";

while (my $row = $tsv->getline ($infh)) {
    $csv->print ($outfh, $row);
}
