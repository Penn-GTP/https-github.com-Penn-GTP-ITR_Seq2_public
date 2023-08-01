#!/bin/env perl
# This script is used to relabel paired-end FASTQ files to add the UMI information
# UMIs are the 3' end bases of the I2 sequences
use strict;
use warnings;

my $usage = "Usage: $0 FWD-TRIMFILE REV-TRIMFILE OUTFILE";
my $in1 = shift or die $usage;
my $in2 = shift or die $usage;
my $out = shift or die $usage;

# check options
unless(defined $in1 && defined $in2) {
	print STDERR "$usage\n";
	exit;
}

# open inputs
if($in1 =~ /\.gz$/) {
	open(IN1, "zcat $in1 |") || die "Unable to open $in1: $!";
}
elsif($in1 =~ /\.bz2$/) {
	open(IN1, "bzcat $in1 |") || die "Unable to open $in1: $!";
}
else {
	open(IN1, "<$in1") || die "Unable to open $in1: $!";
}

if($in2 =~ /\.gz$/) {
	open(IN2, "zcat $in2 |") || die "Unable to open $in2 $!";
}
elsif($in2 =~ /\.bz2$/) {
	open(IN2, "bzcat $in2 |") || die "Unable to open $in2: $!";
}
else {
	open(IN2, "<$in2") || die "Unable to open $in2: $!";
}

# open outputs
open(OUT, ">$out") || die "Unable to write to $out: $!";

# Scan R1/R2 masked lowercase bases
while(!eof(IN1) && !eof(IN2)) {
	my $def1 = <IN1>; chomp $def1;
	my $def2 = <IN2>; chomp $def2;
	my $seq1 = <IN1>; chomp $seq1;
	my $seq2 = <IN2>; chomp $seq2;
	my $sep1 = <IN1>;
	my $sep2 = <IN2>;
	my $qual1 = <IN1>;
	my $qual2 = <IN2>;

	my ($name1, $desc1) = $def1 =~ /^@(\S+)(.*)/;
	my ($name2, $desc2) = $def2 =~ /^@(\S+)(.*)/;
	$name1 .= '/1';
	$name2 .= '/2';

# search lowercase base range
  my $qname = 'primer';
	my $strand1 = '-';
	my $strand2 = '+';
	my ($start1, $end1) = (-1, -1);
	my ($start2, $end2) = (-1, -1);
	my $score = 0;
	if($seq1 =~ /[a-z]+/g) {
		$end1 = pos($seq1);
		$start1 = $end1 - length($&);
		$score += $end1 - $start1;
	}
	if($seq2 =~ /[a-z]+/g) {
		$end2 = pos($seq2);
		$start2 = $end2 - length($&);
		$score += $end2 - $start2;
	}
	print OUT "$name1\t$start1\t$end1\t$name2\t$start2\t$end2\t$qname\t$score\t$strand1\t$strand2\n";
}


close(IN1);
close(IN2);
close(OUT);
