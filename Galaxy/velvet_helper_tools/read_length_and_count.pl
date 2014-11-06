#!/usr/bin/env perl
use strict;
use warnings;

my $file = $ARGV[0];
my $read = $ARGV[1];

print `head -n 2 $file | tail -n 1 | wc -m `-1;
print "\n";
print `grep -c ^$read $file`;
