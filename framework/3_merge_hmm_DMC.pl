#!/usr/bin/perl
#该程序用于合并HMM状态和DMC为一个矩阵
#输入1：hmm/mutstates.txt， hmm/wtstates.txt
#输入2：refumr_CG_DMC_state_stat.txt
#输出：refumr_CG_DMC_stat_HMM.txt
use warnings;
use strict;

open HMM1,"<","hmm/mutstates.txt" or die $!;
my %mut_state_hs = ();
while(<HMM1>){
	s/[\r\n]//g;
	my @arr = split /\t/;
	$mut_state_hs{"$arr[0],$arr[1]"} = $arr[9]; #chr CG -> value (new)
}
close HMM1;

open HMM2,"<","hmm/wtstates.txt" or die $!;
my %wt_state_hs = ();
while(<HMM2>){
	s/[\r\n]//g;
	my @arr = split /\t/;
	$wt_state_hs{"$arr[0],$arr[1]"} = $arr[9];
}
close HMM2;

open DMC,"<","refumr_CG_DMC_state_stat.txt" or die $!;
open OUT,">","refumr_CG_DMC_stat_HMM.txt" or die $!;

my $header = <DMC>;
$header =~ s/[\r\n]//g;
print OUT $header."\tMutState\tWtState\n";
while(<DMC>){
	s/[\r\n]//g;
	my @arr = split /\t/;
	
	if(exists $mut_state_hs{"$arr[0],$arr[1]"} && exists $wt_state_hs{"$arr[0],$arr[1]"}){
		print OUT "$_\t".$mut_state_hs{"$arr[0],$arr[1]"}."\t".$wt_state_hs{"$arr[0],$arr[1]"}."\n";
	}
	else{
		print OUT "$_\t2\t2\n"; #CG not in HMM, No diff.
	}
}
close DMC;
close OUT;

