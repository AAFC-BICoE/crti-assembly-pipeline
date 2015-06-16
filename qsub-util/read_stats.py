#!/usr/bin/env python
import itertools, sys
import HTSeq
import argparse

key_order = ['filename', 'numreads', 'avgreadlen', 'minlen', 'maxlen', 'ns', 'numbases']

parser = argparse.ArgumentParser()
parser.add_argument('--reads', dest='reads', nargs='+',
            help='whitespace-separated list of reads files.')
parser.add_argument('--stats_file', dest='stats_file',
            help='output stats file name.')
parser.add_argument('--stats_header', action='store_true',
            help='write header to output stats file')

args = parser.parse_args()

# function to get the numreads, avg readlen, min, max, number of N's, total bp in a reads file.
def read_stats(filename):
    rfq = HTSeq.FastqReader( filename )
    stats = { 'filename': filename, 'numreads': 0, 'numbases': 0, 'ns': 0, 'minlen': 1e9, 'maxlen': 0 }
    for read in rfq:
        rlen = len(read.seq)
        stats['numreads'] = stats['numreads'] + 1
        stats['numbases'] = stats['numbases'] + rlen
        if rlen < stats['minlen']:
            stats['minlen'] = rlen
        if rlen > stats['maxlen']:
            stats['maxlen'] = rlen
        stats['ns'] = stats['ns'] + read.seq.count('N')

    stats['avgreadlen'] = float (stats['numbases'] / stats['numreads'])
    y = '\t'.join(map(str, [stats[x] for x in key_order])) + '\n'
    return y

sf = open(args.stats_file, 'w')
if (args.stats_header):
    sf.write('#' + '\t'.join(key_order) + '\n')

for rf in args.reads:
    stats_str = read_stats(rf)
    sf.write(stats_str)

sf.close()
