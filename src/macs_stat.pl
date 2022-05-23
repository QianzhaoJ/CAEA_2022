#!/usr/bin/perl
use warnings;
use strict;
die "perl $0 macs_peaks.bed" unless @ARGV==1;
open IN, "<", "$ARGV[0]" || die "$!: $ARGV[0]\n";  ## input is name_peak.bed
my ($read_line, @line_arr, $length, $sum_length, $sum_peaks, $average_length,
    $coverage, @length_arr, $median, $median_length, $row, $temp, $gene_size);
$gene_size = 3090000000;
my ($min_length, $max_length) = (500,0);
my $sample_name = (split(/_/,$ARGV[0]))[0];

print STDERR "---Program\t$0\tstarts --> ".localtime()."\n";
while(<IN>){
	chomp;
	$read_line = $_;
	@line_arr = split(/\s+/,$read_line);
	$length = $line_arr[2] - $line_arr[1];
	push( @length_arr, $length);
	++$sum_peaks;
	$sum_length += $length;
	$min_length = $length if($min_length>$length);
	$max_length = $length if($max_length<$length);
}
$average_length = $sum_length/$sum_peaks;
$coverage = sprintf("%.2f", ($sum_length/$gene_size)*100);
@length_arr = sort {$a <=> $b} @length_arr;
$temp=int($sum_peaks/2);

print join "\t",'Sample','total peaks','total length','mean length','median length','min length','max length','Coverage';
print "\n";
print "$sample_name\t";
print $sum_peaks."\t";
print $sum_length."\t";
print $average_length."\t";
if($sum_peaks/2>$temp){
    $median=($sum_peaks+1)/2;
    $median_length=$length_arr[$median];
    print $median_length."\t";
}
else{
    $median=$sum_peaks/2;
    $median_length=($length_arr[$median]+$length_arr[$median+1])/2;
    print $median_length."\t";
}
print $min_length."\t";
print $max_length."\t";
print $coverage."\%\n";

print STDERR "---Program\t$0\tends  --> ".localtime()."\n";
