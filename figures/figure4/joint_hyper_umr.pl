#!/usr/bin/perl
use warnings;
use strict;

die "
perl $0 H3K4me3_hyper_matrix.gz H3K4me3_umr_matrix.gz H3K4me3_joint_matrix.gz
" unless @ARGV;

my ($hyper_file, $umr_file, $joint_file) = @ARGV;
open IN1,"zcat $hyper_file |" or die $!;
open IN2,"zcat $umr_file |" or die $!;
open OUT,"| gzip >$joint_file" or die $!;

#0,1000,2000,3000
my $header = <IN1>;
$header =~ s/2000/4000/g; #body 2000 to 4000
$header =~ s/0,800,1600/0,1000,2000/g;

print OUT $header;
<IN2>;
while(my $line1 = <IN1>){
	my $line2 = <IN2>;
	my @arr1 = split /\t/, $line1; #hyper
	my @arr2 = split /\t/, $line2; #umr
	
	my @s1_umr_add_arr = @arr2[306..805];
	my @s2_umr_add_arr = @arr2[1106..1605];
	
	my @new_arr = (@arr1[0..505], @s1_umr_add_arr, @arr1[806..1305], @s2_umr_add_arr);
	
	print OUT join "\t", @new_arr;
}

close IN1;
close IN2;
close OUT;
