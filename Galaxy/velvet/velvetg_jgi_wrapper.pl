#!/usr/bin/env perl

# Conventience wrapper for velvetg; copies outfiles to galaxy-specified destinations.
# Please email bugs/feature requests to Edward Kirton (ESKirton@LBL.gov)
#
# History:
# - 2010/03/04 : file created
# - 2001/02/05 : added new options, outfiles; renamed to velvetg_jgi to avoid collision with the other velvetg tool

use strict;
use warnings;
use File::Copy;

# shift wrapper args
my $velveth_path=shift @ARGV or die;
my $velvetg_path=shift @ARGV or die;
my $velvetg_outfile=shift @ARGV or die;
my $contigs_outfile=shift @ARGV or die;
my $stats_outfile=shift @ARGV or die;
my $lastgraph_outfile=shift @ARGV or die;
my $unused_reads_outfile=shift @ARGV or die;
my $amos_outfile=shift @ARGV or die;

# setup velvetg folder

my $exp_cov;
my $size= @ARGV;

for (my $i = 0; $i < $size ; $i ++)
{

	if ($ARGV[$i] eq "-exp_cov") {

		$exp_cov = $ARGV[$i+1];

	}

}

if ($exp_cov eq "fun_time")
{

    my $exp_cov_file=shift@ARGV or die;
    open FH, $exp_cov_file;
    my @array = ();
    while (<FH>)
    {

        push @array , $_ ;

    }

    print @array;
    print @array[0];
    my $size_arr = @array;
    my $r =0;
    my @new_A = ();
    for ($r = 0; $r < $size_arr ; $r ++)
    {



       @{$array[$r]} = split "\t",  $array[$r];
    }
    print "a @array\n";
    for (my $u = 0; $u < $size_arr ; $u ++)
    {
        print $array[$u];
    }


    close FH;


}else {



    if(!(-d $velveth_path))
    {

       $velveth_path = $velveth_path . "*";

       $velveth_path = ` ls -d $velveth_path | head -n 1`;

    }
    chomp $velveth_path;
    #print "My velvetg_path is $velvetg_path\n";
    #print "My velveth path is $velveth_path\n";

    -d $velvetg_path or mkdir($velvetg_path) or die("Unable to create output folder, $velvetg_path: $!\n");

    # run command (remaining args, starting with exe path)
    print "The args @ARGV \n";
    print `cp -r $velveth_path/* $velvetg_path`;

    open (VELVETG, "@ARGV|") or die("Unable to run velvetg\n");
    open (OUT, ">$velvetg_outfile") or die("Unable to open outfile, $velvetg_outfile: $!\n");
    while (<VELVETG>) {
      print OUT $_;
     print if /^Final graph/;
    }
    close VELVETG;
    close OUT;

    # process output
    unlink($contigs_outfile);
    move("$velvetg_path/contigs.fa", $contigs_outfile);
    unlink($stats_outfile);
    move("$velvetg_path/stats.txt", $stats_outfile);

    unlink($lastgraph_outfile);
    if ( -f "$velvetg_path/LastGraph") {
     move("$velvetg_path/LastGraph", $lastgraph_outfile);
    } elsif ( -f "$velvetg_path/Graph2") {
     move("$velvetg_path/Graph2", $lastgraph_outfile);
    } else {
     open(OUT, ">$lastgraph_outfile") or die($!);
     print OUT "ERROR: $velvetg_path/LastGraph not found!\n";
     close OUT;
   }
  unlink($unused_reads_outfile);
  move("$velvetg_path/UnusedReads.fa", $unused_reads_outfile);
    if ( $amos_outfile ne 'None' ) {
      unlink($amos_outfile);
      move("$velvetg_path/velvet_asm.afg", $amos_outfile);
    }
  exit;



}

