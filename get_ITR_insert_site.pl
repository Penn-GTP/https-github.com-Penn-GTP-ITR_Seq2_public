#!/usr/bin/env perl
# This script is ued to get ITR insertion site (of given flanking length) in BED format from dedup-novec alignments
use strict;
use warnings;

my $insert_size = 2;
my $min_softclip = 29;
my $keep_dup = 0;

my $usage = "Usage: $0 INFILE LOCFILE POSFILE [--insert-size $insert_size] [--min-softclip $min_softclip] [--keep-dup]";
my $infile = shift or die $usage;
my $loc_outfile = shift or die $usage;
my $pos_outfile = shift or die $usage;

# parse options
for(my $i = 0; $i < @ARGV; $i++) {
	if($ARGV[$i] eq '--insert-size') {
		$insert_size = $ARGV[++$i];
	}
	elsif($ARGV[$i] eq '--min-softclip') {
		$min_softclip = $ARGV[++$i];
	}
	elsif($ARGV[$i] eq '--keep-dup') {
		$keep_dup = $ARGV[++$i];
	}
	else {
		print STDERR "Error: unknown option $ARGV[$i]\n";
		exit;
	}
}

if(!($insert_size > 0)) {
	print STDERR "--insert-size must be positive\n";
	exit;
}

if(!($min_softclip > 0)) {
	print STDERR "--min-softclip must be positive\n";
	exit;
}


# open BAM input
my $flags = $keep_dup ? "" : "-F 0x400";
open(IN, "samtools view -h $flags $infile | ") || die "Unable to open $infile: $!";
open(LOC, ">$loc_outfile") || die "Unable to write to $loc_outfile: $!";
open(POS, ">$pos_outfile") || die "Unable to write to $pos_outfile: $!";

# read and output
my %chr2len;
while(my $line = <IN>) {
	chomp $line;
	if($line =~ /^@/) {
		if($line =~ /^\@SQ\s+SN:(\S+)\s+LN:(\d+)/) {
			$chr2len{$1} = $2;
		}
	}
	else {
		my ($qname, $flag, $chr, $start, $mapQ, $cigar, $rnext, $pnext, $tlen, $seq, $qual) = split(/\t/, $line);
		my $qlen = length($seq);
		my $strand = ($flag & 0x10) ? '-' : '+';
		my $mate = ($flag & 0x40) ? 1 : 2;
		$start--; # use 0-based start
			my $clip_end;
		if($mate == 1) {
			$clip_end = $strand eq '+' ? 3 : 5;
		}
		else {
			$clip_end = $strand eq '-' ? 3 : 5;
		}

		my $align_len = get_align_len_from_cigar($cigar);
		my $clip_len = get_clip_len_from_cigar($cigar, $clip_end);
		my $end = $start + $align_len;

		my $clip_from = -1;
		my $clip_to = -1;
		if($clip_len >= $min_softclip) {
			if($clip_end == 5 && $strand eq '+' || $clip_end == 3 && $strand eq '-') {
				$clip_from = 0;
				$clip_to = $clip_from + $clip_len;
			}
			else {
				$clip_to = $qlen;
				$clip_from = $clip_to - $clip_len;
			}

			my $insert_pos;
			if($mate == 1) {
				$insert_pos = $strand eq '+' ? $end : $start;
			}
			else {
				$insert_pos = $strand eq '-' ? $end : $start;
			}
			my $insert_start = $insert_pos - $insert_size / 2;
			my $insert_end = $insert_pos + $insert_size / 2;
			if(0 <= $insert_start && $insert_end <= $chr2len{$chr}) {
				print LOC "$chr\t$insert_start\t$insert_end\t$qname/$mate\t$mapQ\t$strand\n";
				print POS "$qname/$mate\t$clip_from\t$clip_to\t$mapQ\t$strand\n";
			}
		}
	}
}

close(IN);
close(LOC);
close(POS);

sub get_align_len_from_cigar {
	my $cigar = shift;
	my $align_len = 0;
	while($cigar =~ /(\d+)([MIDNSHPX=])/g) {
    if($2 eq 'M' || $2 eq 'D' || $2 eq 'N' || $2 eq '=' || $2 eq 'X') {
			$align_len += $1;
		}
	}
	return $align_len;
}

sub get_clip_len_from_cigar {
	my ($cigar, $clip_end) = @_;
	my $clip_regex = $clip_end == 3 ? qr/(\d+)S$/ : qr/^(\d+)S/;
	my $clip_len = 0;
	if($cigar =~ /$clip_regex/) {
		$clip_len = $1;
	}
	return $clip_len;
}
