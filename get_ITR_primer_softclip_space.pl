#!/usr/bin/env perl
# get ITR read space between primer and softclip regions, or negative values if they overlap
use strict;
use warnings;

my $usage = "Usage: $0 TRIM-INFILE SOFTCLIP-INFILE OUTFILE";

my $trim_in = shift or die $usage;
my $clip_in = shift or die $usage;
my $out = shift or die $usage;

open(TRIM_IN, "<$trim_in") || die "Unable to open $trim_in: $!";
open(CLIP_IN, "<$clip_in") || die "Unable to open $clip_in: $!";
open(OUT, ">$out") || die "Unable to write to $out: $!";

# read in primer pos
my %qname2primer_pos;
while(my $line = <TRIM_IN>) {
	chomp $line;
	my ($qname1, $from1, $to1, $qname2, $from2, $to2) = split(/\t/, $line);
	if($from1 != -1 && $to1 != -1) {
		$qname2primer_pos{$qname1} = [$from1, $to1];
	}
	if($from2 != -1 && $to2 != -1) {
		$qname2primer_pos{$qname2} = [$from2, $to2];
	}
}

# read softclip pos and output
print OUT "qname\tprimer_from\tprimer_to\tsoftclip_from\tsoftclip_to\tref_strand\tspace\n";
while(my $line = <CLIP_IN>) {
	chomp $line;
	my ($qname, $clip_from, $clip_to, $score, $ref_strand) = split(/\t/, $line);
	if(exists $qname2primer_pos{$qname}) {
		my ($primer_from, $primer_to) = @{$qname2primer_pos{$qname}};
		my ($mate) = $qname =~ /\/(\d)$/;
		my $space = 'NA';
		if($primer_from < $clip_to && $primer_to > $clip_from) { # primer and clip regions overlap
			$space = $mate == 1 ? $primer_from - $clip_from : $clip_to - $primer_to;
		}
		print OUT "$qname\t$primer_from\t$primer_to\t$clip_from\t$clip_to\t$ref_strand\t$space\n";
	}
}

close(TRIM_IN);
close(CLIP_IN);
close(OUT);
