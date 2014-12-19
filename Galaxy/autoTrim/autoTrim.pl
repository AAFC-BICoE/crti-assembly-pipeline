#!/usr/bin/env perl
use strict;
use warnings;

my $file = shift @ARGV;
my $output = shift @ARGV;
my $input = shift @ARGV;
my $quality = shift @ARGV;
my $cat = shift @ARGV;
print "$file/fastqc_data.txt";




no warnings 'numeric';


#my $first = `grep -n ">>Per sequence quality scores" $file/fastqc_data.txt | cut  -d':' -f 1`;

#print  `tail -n $first $file/fastqc_data.txt | grep -n ">>" | cut -d':' -f 1 | head -n 1`;

#chomp $first;

#chomp $last;

#print "\n$first"; #,$last";

open R , ">$cat";

open FH, "$file/fastqc_data.txt";

my $score = "False";

open OUT1, ">pbsq.txt";

#while (my $line = <FH>)
#{
#	chomp $line;
#	print "L$line\n";
#	>>Per base sequence quality	pass
#	if ( index($line,">>Per base sequence quality") != -1 && $line ne ">>Per base sequence quality	pass")
#	{
#		$score = "True";
#		next;	
#	}	
#	if ($score eq "True")
#	{
#		if ($line eq ">>END_MODULE")
#		{
#			last;	
#		}	
#		print OUT1 "$line\n";
#		print "$line\n";
#	}
#}

while (my $line = <FH>)
{
	chomp $line;
#	>>Per base sequence quality	pass
	if ( index($line,">>Per base sequence quality") != -1 && $line ne ">>Per base sequence content	pass")
	{
		$score = "True";
		next;	
	}

	if ($score eq "True")
	{
		if ($line eq ">>END_MODULE")
		{
			last;	
		}
		print OUT1 "$line\n";
	#	print "$line\n";
	}
}




close FH;

close OUT1;
#print "DONE\n";

open T, "pbsq.txt";


my $first = 1;
my $last = `tail -n 1 pbsq.txt | cut -d'\t' -f 1 | cut -d'-' -f 2`;


#my $first2;
#my $last2 = `tail -n 1 pbsq.txt | cut -d'\t' -f 1 | cut -d'-' -f 2`;

my $bool = 0;
my @fail;
my @list;
my $i = 0;
while (my $l = <T>)
{
	chomp $l;
	if ($l eq "#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile")
	{
		next;
	}
	my @A = split (/\t/, $l);
#	print "A@A\n";
	if ($A[2] < 25 || $A[3] < 10 )
	{
#		if(index($A[0],"-") != -1)
#		{
#			chomp $A[0];
#			push @fail, `echo $A[0] |cut -d'-' -f 2`;
#		}
#		else
#		{
#			chomp $A[0];
#			push @fail, $A[0];		
#		}
#		print "f $A[0]\n";

		push @fail, $i;

	}


	


	if(index($A[0],"-") != -1)
	{
		chomp $A[0];
		my $b = `echo $A[0] |cut -d'-' -f 2`;
		chomp $b;
		push @list, $b ;
	}
	else
	{
		chomp $A[0];
		push @list, $A[0];		
	}
	$i ++;
}
close T;

##print "LL @list LL\n";
#print `wc pbsq.txt`;
#print `cat pbsq.txt `;
#print
#print `tail -n 1 pbsq.txt `;
#print "\n";
#print `tail -n 1 pbsq.txt | cut -d'\t' -f 1`;
#print "\n";
#print "LASSSSSSSSSSSSSSSSSSSSSSSST $last \n";
open FH, "$file/fastqc_data.txt";

$score = "False";

open OUT1, ">args.txt";

while (my $line = <FH>)
{
	chomp $line;
#	>>Per base sequence quality	pass
	if ( index($line,">>Per base sequence content") != -1 && $line ne ">>Per base sequence content	pass")
	{
		$score = "True";
		next;	
	}
	
	if ($score eq "True")
	{
		if ($line eq ">>END_MODULE")
		{
			last;	
		}
		print OUT1 "$line\n";
	#	print "$line\n";
	}
}

close FH;

close OUT1;
#print `wc args.txt`;
open T, "args.txt";



my @B;

while(my $l = <T>)
{
	chomp $l;
	if ($l eq "#Base	G	A	T	C")
	{
		next;
	}
	push @B, $l;
}
close T;
my $size = @B;
my @fail2;

for (my $i = 0; $i < $size ; $i = $i + 1)
{
	
	#
	# Check that we are not in a point that causes a warn or a fail
	#
	
	my @C = split /\t/, $B[$i];
	print R "$C[0] ";
	
#	if ( abs($C[1] - $C[4] )/ (($C[1] + $C[4])/2) > 0.1)
#	{
	print R "GC ",abs($C[1] - $C[4] ) / (($C[1] + $C[4])/2);	
#	}

#	if (abs($C[2] - $C[3] )/ (($C[2] + $C[3])/2) > 0.1 )
#	{
	print R " AT ", abs($C[2] - $C[3] )/ (($C[2] + $C[3])/2) , "\n";
#	}
	
	if (abs($C[1] - $C[4] )/ (($C[1] + $C[4])/2) > 0.1 || abs($C[2] - $C[3] )/ (($C[2] + $C[3])/2) > 0.1 )
	{
#		print "bb $B[$i] bb \n";

#		if ( abs($C[1] - $C[4] )/ (($C[1] + $C[4])/2) > 0.1)
#		{
#			print "GC ",abs($C[1] - $C[4] / (($C[1] + $C[4])/2)), "\n";	
#		}

#		if (abs($C[2] - $C[3]) / (($C[2] + $C[3])/2) > 0.1 )
#		{
#			print "AT ", abs($C[2] - $C[3] / (($C[2] + $C[3])/2)) , "\n";
#		}

		push @fail2 , $i;

#		if(index($C[0],"-") != -1)
#		{
#			push @fail2, `echo $C[0] |cut -d'-' -f 2`;	
#		}
#		else
#		{
#			push @fail2, $C[0];	
#		}
	}
}

my $j = 0;
my $break = 0;

#print "F1 @fail \n";
#print "F2 @fail2 \n";

# combining the two arrays

for (my $k = 0; $k < @fail ; $k ++)
{
	if ($fail[$k] <= $fail2[$j])
	{
		splice @fail2 , $j , 0 , $fail[$k] ;
		splice @fail , $k;
	}else
	{
		while ($fail[$k] > $fail2[$j] )
		{
			$j++;
			## Break out of for loop here
			if ($j >= @fail2)
			{
				$break = 1;
				last;
			}
		}
		
		if ($break == 1)
		{
			last;
		}
		
		splice @fail2 , $j , 0 , $fail[$k] ;
		splice @fail , $k;
		
	}
}

my $firstFail = 0;
my $lastFail = 0;

print "The list values: $list[$fail2[0]] $list[$fail2[-1]] \n";

if ( @fail2 && $first != $list[$fail2[0]])
{
	$firstFail = 1;
}


if ( @fail2 && $last != $list[$fail2[-1]])
{
	$lastFail = 1;
}




###########################################################
#							  # 
#	Add first and last to the end and start of array  #
#	Also remove duplicates				  #
#	HERE						  #
#							  #
###########################################################


my $prev = -1;
#unshift (@fail2 , "1");
splice @fail2, 0, 0, 0;
my $t = @list - 1 ;
push (@fail2, $t);
#my $s = 0;#@fail2;
my $r = 0;
#$s = @fail2;
#chomp $s;
print "FYUESES @fail2 \n";
#print @list;
#print "S$s\n";

$prev = 0;
for (my $t = 1; $t < @fail2; $t ++)
{
	if ( $prev == $fail2[$t])
	{
#		print $fail2[$t];
		splice @fail2 , $t , 1;
		$t = $t - 1;
	}else
	{
		$prev = $fail2[$t];
	}
}






print "FAIL @fail2";

my $leng = -1;
my $start = -1;
my $end = -1;


for (my $j = 0; $j < @fail2 - 1 ; $j ++)
{	
#	print "\nstart $start , end $end\n";
	if ($list[$fail2[$j+1]] - $list[$fail2[$j]] >=  $leng)
	{
		$leng = $fail2[$j+1] - $fail2[$j];
		$start = $fail2[$j];
		$end = $fail2[$j+1];
	}
	
}

#print "\nstart $start , end $end\n";

chomp $start;
chomp $end;

## moving them up and down by 1

print "start: $list[$start] end: $list[$end] \n";


print "FF $firstFail";

if ((@list && $list[$start] != 1 )|| $firstFail == 0)
{
	$start = $start + 1;
}

if (@list &&  $list[$end] != 1 && $lastFail == 0)
{
	$end = $end - 1;
}

#33 and 64
print "fastx_trimmer -Q $quality -f $list[$start] -l $list[$end] -i $input -o $output";
print `fastx_trimmer -Q $quality -f $list[$start] -l $list[$end] -i $input -o '$output'`;

close R;
print "\n $start , $end\n";

#die;
