#!/usr/bin/perl
#获取CRs的TF interaction network
use warnings;
use strict;

open CR,"<","chromatin_regulators.txt" or die $!;
my %CR_hs = ();
while(<CR>){
	s/[\r\n]//g;
	$CR_hs{$_} = 1;
}
close CR;

my %TF_hs = ();
open TF,"<","TFs.txt" or die $!;
while(<TF>){
	s/[\r\n]//g;
	my @arr = split /\t/;
	$TF_hs{uc($arr[0])} = 1;
}
close TF;

open BIOGRID,"<","sort_uniq_BIOGRID_intersect.txt" or die $!;
while(<BIOGRID>){
	s/[\r\n]//g;
	my @arr = split /\t/;
	if( (exists $CR_hs{$arr[0]} && exists $TF_hs{$arr[1]}) ||
		(exists $CR_hs{$arr[0]} && exists $TF_hs{$arr[1]}) ){
		print $_."\n";
	}
}
close BIOGRID;

