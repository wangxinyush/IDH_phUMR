#!/usr/bin/perl
#这个程序用来获取所有在refumrs里的CG位点甲基化水平，样本为Brain Normal,IDH Glioma, WT Glioma
use warnings;
use strict;
use File::Basename;

warn call_time()."Start\n";
#1.get all CGs in refumrs
my $reference = "/share/pub/wangxy/IDH_Glioma/1_Identificaiton_IDH_Hyper/202111_new_refumrs/refumr_CG_methy/mr200_refumr_0.6_IGV.bed"; #21716
#my $wig_list = "test_wig_list.txt";
my $wig_list = "Brain_GBM_wig_list.txt";

my $ref_head = 0;				#reference has a head?, default: 0
my $chr_ref_hs_ref = bed_chr_split( $reference , $ref_head );  # chr => [start , end]

my $chr_site_hs_ref = CG_chr_split("/share/pub/wangxy/software/genome/ucsc/hg19/CG/hg19_CpG.txt");

#my $region_CG_ref = ""; #"chr\tcg" => ""
my %CG_region_hs = (); #"chr\tcg" => "chr\tstart\tend" (region, eg. refumr)
foreach my $chr(sort keys %{ $chr_ref_hs_ref } ){ #get CG in refumrs for each chromosome.
	my ($chr_region_CG_ref, $chr_CG_region_ref) = get_region_cg( $chr_ref_hs_ref -> {$chr} , $chr_site_hs_ref -> {$chr}, $chr);
	foreach my $cg(keys %{ $chr_CG_region_ref }){
		$CG_region_hs{$cg} = $chr_CG_region_ref -> {$cg};
	}
}

#methy matrix (merge)
my @sample_wig_files = get_sample_wigs($wig_list);
my %methy_matrix_hs = ();
my $header = "chr\tCG\tchr\tstart\tend\t";
foreach my $sample_i(0..$#sample_wig_files){  #for each sample. 
	open IN,"<",$sample_wig_files[$sample_i] or die $!;
	my $chr = "0";
	
	my %cg_in_refumr_hs = (); #"chr\tcg" => methy_level
	while(<IN>){
		next if(/track type/);
		s/[\r\n]$//g;
		if(/chrom/){ 	#declaration line
			my @declaration_arrs = split /\s+/;
			$chr = substr( $declaration_arrs[1], 6 ); #get chr from declaration line, eg. variableStep chrom=chr1
		}
		else{
			my @arr = split ;
			if(exists $CG_region_hs{"$chr\t$arr[0]"}){ #TCGA, Feinberg
				$cg_in_refumr_hs{"$chr\t$arr[0]"} = $arr[1];
			}
			elsif( exists $CG_region_hs{ "$chr\t".($arr[0]+1) }){ #DKFZ, Roadmap, need to add 1
				$cg_in_refumr_hs{ "$chr\t".($arr[0]+1) } = $arr[1];
			}
		}
	}
	#merge into matrix
	foreach my $cg(keys %CG_region_hs){
		if(exists $cg_in_refumr_hs{$cg}){
			$CG_region_hs{$cg} .= "\t".$cg_in_refumr_hs{$cg};
		}
		else{
			$CG_region_hs{$cg} .= "\tNA";
		}
	}
	
	close IN;
}

#print into matrix
#open OUT,">","test4_refumr_CG_methy_matrix.txt" or die $!;
open OUT,">","refumr_CG_methy_matrix.txt" or die $!;
print OUT "chr\tCG\tchr\tstart\tend";
foreach my $sample_path(@sample_wig_files){
	my $sample_name = fileparse($sample_path,".wig");
	print OUT "\t".$sample_name;
}
print OUT "\n";


foreach my $cg(sort keys %CG_region_hs){
	print OUT "$cg\t".$CG_region_hs{$cg}."\n";
}
close OUT;

warn call_time()."End\n";

=DESCRIPTION CG_chr_split  	[ INDEPENDENT ]
    * Name: CG_chr_split
    * Function: split file into chr arrs.
    * Params: $infile
    * Return: \@chr_methy_arrs , \@chr_site_arrs
    * Independence: [ INDEPENDENT ]
=cut

sub CG_chr_split{
	my ( $infile ) = @_;

	open IN,"<",$infile;
	
	my %chr_site_hs  = ();

	my $chr = "0";
	while( <IN> ){
		s/[\r\n]$//g;
		my ($chr, $site) = split /\t/;
		if(exists $chr_site_hs{$chr}){
			push @{ $chr_site_hs{$chr} } , $site;
		}
		else{
			$chr_site_hs{$chr}  = [];
		}
	}
	close IN;
	return ( \%chr_site_hs );
}

=Data Description
@hg19 CG (sorted)
chr1    10469
chr1    10471
chr1    10484
@refumr (sorted)
chr1    713543  714842
chr1    762105  762804
chr1    804991  805540
=cut
=Usage
my $ref = get_region_cg(\@region,\@site);
foreach my $i(@$ref){
	print $i."\n";
}
#date: 10/10/2021
=cut
#get_region_cg(\@region,\@site)
sub get_region_cg{
	my ( $region_ref, $site_ref, $cur_chr ) = @_;
	
	my $site_index = 0;
	my $first_right_site_index =0;
	
	my %region_CG_hs = ();
	my %CG_region_hs = ();
	
	foreach my $i(0..$#$region_ref){
		my ( $start , $end ) = ( $$region_ref[$i][0] , $$region_ref[$i][1] );
		
		my $num = 0;

		while( $$site_ref[$site_index] < $start ){
			last if ( $site_index >= $#$site_ref );
			$site_index ++;
		}#end: maybe the first site in region. (>= start)
		$first_right_site_index = $site_index ; #mark first mapping site

		while( $$site_ref[$site_index] <= $end ){	#all sites in region. ( >= start and <= end )
			last if ( $site_index >= $#$site_ref );
			
			$region_CG_hs{"$cur_chr\t".$$site_ref[$site_index]} = "";
			$CG_region_hs{"$cur_chr\t".$$site_ref[$site_index]} = "$cur_chr\t$start\t$end";
			
			$num ++;
			$site_index ++;
		}#end: out region
		
		$site_index = $first_right_site_index ; #to return first mapping site
	}
	
	return (\%region_CG_hs, \%CG_region_hs);
}


sub get_sample_wigs{
	my ( $wig_list ) = @_ ;

	my @sample_wig_files = ();
	if( $wig_list =~ /,/ ){ #mutiple wig file
		@sample_wig_files = split /,/,$wig_list;
	}
	elsif( $wig_list =~ /wig$/ or $wig_list =~ /wig\.gz$/){ #single wig file
		push @sample_wig_files, $wig_list;
	}
	elsif( -e $wig_list ){ #file exists, and not a wig file.
		open FL,"<",$wig_list or die "$wig_list \n wig file open error, Please check up this path!";
		@sample_wig_files = <FL>;

		map s/[\r\n]$//g,@sample_wig_files;

		close FL;
	}
	else{
		die "wig lists must be a file or a string seperated by comma.\n";
	}
	return @sample_wig_files;
}

=DESCRIPTION bed_chr_split  	[ INDEPENDENT ]
    * Name: bed_chr_split
    * Function: split bed file into chr arrs.
    * Params: $infile , $has_head
    * Return: \@chr_bed_hs
    * Independence: [ INDEPENDENT ]
=cut

sub bed_chr_split{
	my ( $infile , $has_head ) = @_;
	open my $in_fh,"<",$infile or die $infile."\nBed region open error,Please check up this path!\n";
	
	#delete the head 
	if($has_head){
		<$in_fh>;
	}
	
	my %chr_bed_hs  = ();

	my $chr = "";
	while( <$in_fh> ){
		s/[\r\n]$//g;
		my @arr = split;
		my ( $curChr, $start, $end ) = @arr[0, 1, 2];

		if($curChr ne $chr){  	# new chromosome
			$chr = $curChr;
			$chr_bed_hs{$chr}  = [];
			push @{ $chr_bed_hs{$chr} } , [ $arr[1], $arr[2] ];
		}
		else{
			push @{ $chr_bed_hs{$chr} } , [ $arr[1], $arr[2] ];
		}
	}
	close $in_fh;
	return \%chr_bed_hs ;
}

=DESCRIPTION wig_chr_split  	[ INDEPENDENT ]
    * Name: wig_chr_split
    * Function: split file into chr arrs.
    * Params: $infile
    * Return: \@chr_methy_arrs , \@chr_site_arrs
    * Independence: [ INDEPENDENT ]
=cut

sub wig_chr_split{
	my ( $infile ) = @_;

	my $in_fh;
	if( $infile =~ /gz$/ ){
		open $in_fh,"zcat $infile |" or die $infile."\nwig file open error, Please check up this path!\n";
	}
	elsif( $infile =~ /wig$/ ){
		open $in_fh,"<",$infile or die $infile."\nwig file open error, Please check up this path!\n";
	}
	else{
		die "ERROR: $infile is not .wig or .gz file.\n";
	}
	
	my %chr_site_hs  = ();
	my %chr_methy_hs = ();

	my $chr = "0";
	while( <$in_fh> ){
		next if(/track type/);
		s/[\r\n]$//g;
		if(/chrom/){ 	#declaration line
			my @declaration_arrs = split /\s+/;
			$chr = substr( $declaration_arrs[1], 6 ); #get chr from declaration line, eg. variableStep chrom=chr1
			$chr_site_hs{$chr}  = [];
			$chr_methy_hs{$chr} = [];
		}
		else{
			my @arr = split ;
			push @{ $chr_site_hs{$chr}  } , $arr[0];
			push @{ $chr_methy_hs{$chr} } , $arr[1];
		}
	}
	close $in_fh;
	return ( \%chr_site_hs , \%chr_methy_hs );
}

#get current time [hour:min:sec month,mday]
sub call_time{
	my $tab_nums = 0; #DEFAULT: has no tab 
	$tab_nums = $_[0] if ($_[0]);

	my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime;
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	#my @weeks = qw(Sun Mon Tue Wed Thu Fri Sat);
	#$year += 1900;
	
	return "[$hour:$min:$sec $months[$mon],$mday] ".("   " x $tab_nums);
}

