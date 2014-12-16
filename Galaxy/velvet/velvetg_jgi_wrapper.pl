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

sub getParam # (my $lookFor , my @arr) 
{

	my $lookFor = shift @_;
	my @arr = @_;
	my $ret = -1;
	
	for(my $i = 0; $i < @arr ; $i ++)
	{
		
		if ($arr[$i] eq $lookFor)
		{
			$ret = $i + 1;
		}
		
	}
	
	return $ret;
	
}






# shift wrapper args
my $velveth_path=shift @ARGV or die;
my $velvetg_path=shift @ARGV or die;
my $velvetg_outfile=shift @ARGV or die;
#my $contigs_outfile=shift @ARGV or die;
#my $stats_outfile=shift @ARGV or die;
#my $lastgraph_outfile=shift @ARGV or die;
#my $velvetg_log_file= shift @ARGV or die;

my $kmer = shift @ARGV or die;

my $kmerBool = shift @ARGV or die;

my $outputID = shift @ARGV or die;

my $outpath = shift @ARGV or die;

my $extension = shift @ARGV or die;




my $exp_cov_file ;


if (($kmerBool  eq "false"))
{
	$exp_cov_file = shift @ARGV or die;
}


#my $unused_reads_outfile=shift @ARGV or die;
#my $amos_outfile=shift @ARGV or die;

# setup velvetg folder


chomp $velveth_path;



open (OUT, ">$velvetg_outfile") or die("Unable to open outfile, $velvetg_outfile: $!\n");


if ($kmerBool eq "true" )
{
	print "Such TRUE\n";
	$velveth_path = $velveth_path . "/kmer";
	if (!( -d $velveth_path))
	{
	
	    $velveth_path = $velveth_path . "_" . $kmer ;
	    print "AJLKDFS";
	}

	-d $velvetg_path or mkdir($velvetg_path) or die("Unable to create output folder, $velvetg_path: $!\n");

 	# run command (remaining args, starting with exe path)
#	$ARGV[1] = "$velvetg_path" . 
#	print  "Velveth path is $velveth_path \n";
#	print  "The args @ARGV \n";
	print  `cp -r $velveth_path/* $velvetg_path`;
	if ($kmer <= 31)
	{
		$ARGV[0] = "velvetg31";
	}
	elsif ($kmer <= 63  )
	{
		$ARGV[0] = "velvetg63";
	}
	elsif ($kmer <= 127  )
	{
		$ARGV[0] = "velvetg127";
	}
	elsif ($kmer <= 145 )
	{
		$ARGV[0] = "velvetg145";
	}
	elsif ($kmer <= 195  )
	{
		$ARGV[0] = "velvetg195";
	}
	else
	{
		$ARGV[0] = "velvetg245";
	}
	print "THE COMMAND @ARGV";
	open (VELVETG, "@ARGV|") or die("Unable to run velvetg\n");
	#open (OUT, ">$velvetg_outfile") or die("Unable to open outfile, $velvetg_outfile: $!\n");

	




	while (<VELVETG>) {
		print OUT $_;
		print if /^Final graph/;
	}
	close VELVETG;
	
	chomp $outpath;
	chomp $outputID;
	chomp $kmer;
	chomp $extension;
	my $loc = $outpath . "/primary_" . "$outputID" . "_contigsKmer$kmer" . "_visible_" . "$extension"      ;
		
	
	chomp $velvetg_path;

	


#	print OUT "Loc $velvetg_path/contigs.fa";
	open FH, "< $velvetg_path/contigs.fa" or die "Cannot open file $velvetg_path/contigs.fa : $!";
	chomp $loc;
	print $loc;
	open (O , ">$loc");
		
	while (my $line = <FH> )
	{
		chomp $line;
		print O $line, "\n";
	}

	close O;
	
	
	


	$loc = $outpath . "/primary_" . "$outputID" . "_statsKmer$kmer" . "_visible_txt"      ;
	chomp $velvetg_path;
#	print OUT "Loc $velvetg_path/stats.txt.fa";
	open FH, "< $velvetg_path/stats.txt" or die "Cannot open file: $!";
	chomp $loc;
	print $loc;
	open (O , ">$loc");
		
	while (my $line = <FH> )
	{
		chomp $line;
		print O $line, "\n";
	}

	close O;


	
	$loc = $outpath . "/primary_" . "$outputID" . "_PreGraphKmer$kmer" . "_visible_txt";
	chomp $velvetg_path;
	print OUT "Loc $velvetg_path/PreGraph";
	open FH, "< $velvetg_path/PreGraph" or die "Cannot open file: $!";
	chomp $loc;
	print $loc;
	open (O , ">$loc");
		
	while (my $line = <FH> )
	{
		chomp $line;
		print O $line, "\n";
	}

	close O;	
	
	print OUT "getParam: ",  getParam("-amos_file", @ARGV); 
	if (getParam ("-amos_file", @ARGV) != -1 && $ARGV[ getParam ("-amos_file", @ARGV) ] eq "yes" )
		{
			$loc = $outpath . "/primary_" . "$outputID" . "_VelvetAsmKmer$kmer" . "_visible_txt";
			chomp $velvetg_path;
			print OUT "Loc $velvetg_path/velvet_asm.afg";
			open FH, "< $velvetg_path/velvet_asm.afg" or die "Cannot open file: $!";
			chomp $loc;
			print $loc;
			open (O , ">$loc");
		
			while (my $line = <FH> )
			{
				chomp $line;
				print O $line, "\n";
			}
	
			close O;	
		}
		
	print OUT "getParam: ",  getParam("-unused_reads", @ARGV); 
	if (getParam ("-unused_reads", @ARGV) != -1 && $ARGV[ getParam ("-unused_reads", @ARGV) ] eq "yes" )
	{
		$loc = $outpath . "/primary_" . "$outputID" . "_UnusedReadsKmer$kmer" . "_visible_fa";
		chomp $velvetg_path;
		print OUT "Loc $velvetg_path/UnusedReads.fa";
		open FH, "< $velvetg_path/UnusedReads.fa" or die "Cannot open file: $!";
		chomp $loc;
		print $loc;
		open (O , ">$loc");
	
		while (my $line = <FH> )
		{
			chomp $line;
			print O $line, "\n";
		}
		close O;	
	}
		
	
	

	
	print `rm $velvetg_path/*`;

}else
{
	
	my @kmers = ` cut -f1 $exp_cov_file `;
	my @expected = ` cut -f2 $exp_cov_file`;
#	print "expected coverage file: $exp_cov_file \n";
	print OUT @kmers;
	print OUT @expected;
	my $expectedCov = getParam ("-exp_cov" , @ARGV) ;
#	for (my $a = 0; $a < @ARGV ; $a ++)
#	{
#		if ($ARGV[$a] eq "-exp_cov")
#		{
#			$expectedCov = $a + 1;
#		}
#	}
	-d $velvetg_path or mkdir($velvetg_path) or die("Unable to create output folder, $velvetg_path: $!\n");
#	print `mkdir $velvetg_path/temp`;
	
#	print `cp -r $velveth_path/* $velvetg_path/temp`;

#	print $velveth_path;
	print `mkdir $velvetg_path/temp`;
	
	for (my $i = 0; $i < @kmers -1; $i ++)
	{
		chomp $velveth_path;
		chomp $velvetg_path;
		chomp $kmers[$i];
		my $rty = $velveth_path . "/kmer_$kmers[$i]" ;
		
		print OUT $rty;		
		
		
		if ($kmers[$i] <= 31)
		{
			$ARGV[0] = "velvetg31";
		}
		elsif ($kmers[$i] <= 63  )
		{
			$ARGV[0] = "velvetg63";
		}
		elsif ($kmers[$i] <= 127  )
		{
			$ARGV[0] = "velvetg127";
		}
		elsif ($kmers[$i] <= 145 )
		{
			$ARGV[0] = "velvetg145";
		}
		elsif ($kmers[$i] <= 195  )
		{
			$ARGV[0] = "velvetg195";
		}
		else
		{
			$ARGV[0] = "velvetg245";
		}

		chomp $rty;
		print OUT `cp $rty/*  $velvetg_path/temp `;
		chomp $expected[$i];
		$ARGV[1] = "$velvetg_path" . "/temp" ;
		$ARGV[$expectedCov] = $expected[$i];
		print OUT `ls $ARGV[1]`;
		print OUT "The A: @ARGV \n";
		print OUT "AAA\n";

	#	print "echo $expected[$i] > $velvetg_path/temp/contigs.fa";
	#	print `echo $expected[$i] > $velvetg_path/temp/contigs.fa`;
		
		open (VELVETG, "@ARGV|") or die("Unable to run velvetg\n");

		while (<VELVETG>) {
			print OUT $_;
			print if /^Final graph/;
		}
		close VELVETG;
		print "ls $ARGV[1]";
		print OUT`ls $ARGV[1]`;
		chomp $outpath;
		chomp $outputID;
		chomp $kmers[$i];
		chomp $extension;
		my $loc = $outpath . "/primary_" . "$outputID" . "_contigs Kmer $kmers[$i]" . "_visible_" . "$extension"      ;
		
		# $outpath ."/primary_".$outid."_$data"."_visible_$intype";
#		open ( OUT , ">$loc");
#		print OUT $expected[$i], "\n";

		chomp $velvetg_path;
		print OUT "Loc $velvetg_path/temp/contigs.fa";
		open FH, "< $velvetg_path/temp/contigs.fa" or die "Cannot open file: $!";
		chomp $loc;
		print $loc;
		open (O , ">$loc");
		
		while (my $line = <FH> )
		{
			chomp $line;
			print O $line, "\n";
		}

		close O;
		close FH;


		$loc = $outpath . "/primary_" . "$outputID" . "_statsKmer$kmers[$i]" . "_visible_txt"      ;

		chomp $velvetg_path;
		print OUT "Loc $velvetg_path/temp/stats.txt";
		open FH, "< $velvetg_path/temp/stats.txt" or die "Cannot open file: $!";
		chomp $loc;
		print $loc;
		open (O , ">$loc");
		
		while (my $line = <FH> )
		{
			chomp $line;
			print O $line, "\n";
		}

		close O;
		close FH;

		$loc = $outpath . "/primary_" . "$outputID" . "_PreGraphKmer$kmers[$i]" . "_visible_txt";
		chomp $velvetg_path;
		print OUT "Loc $velvetg_path/temp/PreGraph";
		open FH, "< $velvetg_path/temp/PreGraph" or die "Cannot open file: $!";
		chomp $loc;
		print $loc;
		open (O , ">$loc");
		
		while (my $line = <FH> )
		{
			chomp $line;
			print O $line, "\n";
		}

		close O;
		close FH;	
	
		print OUT "\ngetParam: ",  getParam("-amos_file", @ARGV); 
		
		if (getParam ("-amos_file", @ARGV) != -1 && $ARGV[ getParam ("-amos_file", @ARGV) ] eq "yes" )
		{
			$loc = $outpath . "/primary_" . "$outputID" . "_VelvetAsmKmer$kmers[$i]" . "_visible_txt";
			chomp $velvetg_path;
			print OUT "Loc $velvetg_path/temp/velvet_asm.afg";
			open FH, "< $velvetg_path/temp/velvet_asm.afg" or die "Cannot open file: $!";
			chomp $loc;
			print $loc;
			open (O , ">$loc");
		
			while (my $line = <FH> )
			{
				chomp $line;
				print O $line, "\n";
			}
	
			close O;
			close FH;	
		}
		
		print OUT "getParam: ",  getParam("-unused_reads", @ARGV); 
		
		if (getParam ("-unused_reads", @ARGV) != -1 && $ARGV[ getParam ("-unused_reads", @ARGV) ] eq "yes" )
		{
			$loc = $outpath . "/primary_" . "$outputID" . "_UnusedReadsKmer$kmers[$i]" . "_visible_fa";
			chomp $velvetg_path;
			print OUT "Loc $velvetg_path/temp/UnusedReads.fa";
			open FH, "< $velvetg_path/temp/UnusedReads.fa" or die "Cannot open file: $!";
			chomp $loc;
			print $loc;
			open (O , ">$loc");
		
			while (my $line = <FH> )
			{
				chomp $line;
				print O $line, "\n";
			}
	
			close O;	
			close FH;
		}
		
		
		
		
		
		
		
	
		print `rm $velvetg_path/temp/*`;
		
	#	print "mv $velvetg_path/temp/contigs.fa $loc";
	#	print `cp $velvetg_path/temp/contigs.fa $loc`;
	#	print "\n LS \n";
	#	print `ls $outpath`;
	#	print "\n END LS \n";
	}
	
		
}

close OUT;



#print ` mkdir $velvetg_path/apps`;


#`cp $velvetg_path/contigs.fa $velvetg_path/apps`;
#`cp $velvetg_path/stats.txt $velvetg_path/apps`;





    # process output
  #unlink($velvetg_log_file);
#  move("$velvetg_path/Log", $velvetg_log_file);
#  unlink($contigs_outfile);
#  move("$velvetg_path/contigs.fa", $contigs_outfile);
#  unlink($stats_outfile);
#  move("$velvetg_path/stats.txt", $stats_outfile);
#
#  unlink($lastgraph_outfile);
#  if ( -f "$velvetg_path/PreGraph") {
 #   move("$velvetg_path/PreGraph", $lastgraph_outfile);
#  } elsif ( -f "$velvetg_path/Graph2") {
#    move("$velvetg_path/Graph2", $lastgraph_outfile);
#  } else {
 #   open(OUT, ">$lastgraph_outfile") or die($!);
#    print OUT "ERROR: $velvetg_path/LastGraph not found!\n";
#    close OUT;
#  }
#  unlink($unused_reads_outfile);
#  move("$velvetg_path/UnusedReads.fa", $unused_reads_outfile);
#    if ( $amos_outfile ne 'None' ) {
#      unlink($amos_outfile);
#      move("$velvetg_path/velvet_asm.afg", $amos_outfile);
#   }




exit;
