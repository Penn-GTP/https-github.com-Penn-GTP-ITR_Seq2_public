#!/bin/env perl
# Prepare sh script for filtering reference mapping files
our $VERSION = 'v2.2.4';

use strict;
use warnings;
use File::Basename;
use lib dirname (__FILE__);
use ITRSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
#my $sh_path = '/bin/bash';
my $samtools = 'samtools';
my $bedtools = 'bedtools';
my $cmd = "$0 " . join(" ", @ARGV);
my @comments = qq(Sample name\tTotal reads\tITR-containing reads\tHost-mapped reads\tHost-mapped deduplexed reads\tVector-mapped reads\tHost-mapped deduplexed non-vector reads\tInsert sites\tInsert sites filtered\tUnique insert sites\tUnique insert sites filtered\tMerged insert peaks\tRead abundance of insert peaks\tOn-target insert peaks\tRead abundance of on-target insert peaks\tOff-target insert peaks\tRead abundance of off-target insert peaks\tClonal insert sites\tUMI-locus abundance of clonal insert sites\tFrequency of UMI locus abundance of clonal insert sites);
my @headers = qw(sample_name total_read trimmed_read ref_mapped ref_mapped_dedup vec_mapped ref_mapped_dedup_novec
insert_site insert_site_filtered insert_site_uniq insert_site_uniq_filtered
peak_count peak_clone ontarget_peak_count ontarget_peak_clone offtarget_peak_count offtarget_peak_clone
clonal_count clonal_loc_count clonal_loc_freq);

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new ITRSeqExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $BASE_DIR = $design->get_global_opt('BASE_DIR');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
my $DEMUX_DIR = $design->get_global_opt('DEMUX_DIR');
my $WORK_DIR = $design->get_global_opt('WORK_DIR');
#my $UMI_LEN = $design->get_global_opt('UMI_LEN');
my $DEFAULT_MIN_CLONE_EXP = 2;

# check required directories
if(!(-e $BASE_DIR && -d $BASE_DIR)) {
	print STDERR "Error: BASE_DIR $BASE_DIR not exists\n";
	exit;
}

if(!(-e $SCRIPT_DIR && -d $SCRIPT_DIR)) {
	print STDERR "Error: SCRIPT_DIR $SCRIPT_DIR not exists\n";
	exit;
}

if(!(-e $DEMUX_DIR)) {
	  print STDERR "Error: DEMUX $DEMUX_DIR not exists\n";
		  exit;
}

if(!(-e $WORK_DIR)) {
	print STDERR "Error: WORK_DIR $WORK_DIR not exists\n";
	exit;
}

open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
# write header
print OUT qq(# CMD:"$cmd"\n# VER:$VERSION\n);
print OUT "# ", join("\t", @comments), "\n";
print OUT join("\t", @headers), "\n";

foreach my $sample ($design->get_sample_names()) {
	print STDERR "getting stats for $sample\n";
# get total read
  my $total_read;
	{
		my $in = $design->sample_opt($sample, 'fastq_R1');
		$total_read = $in =~ /\.gz$/ ? `zcat $DEMUX_DIR/$in | wc -l` : `cat $DEMUX_DIR/$in | wc -l`;
		chomp $total_read;
		$total_read /= 4;
	}

# get trimmed read
  my $trimmed_read;
  {
		my $in = $design->get_sample_fwd_ITRtrim_file($sample);
		$trimmed_read = $in =~ /\.gz$/ ? `zcat $WORK_DIR/$in | wc -l` : `cat $DEMUX_DIR/$in | wc -l`;
		chomp $trimmed_read;
		$trimmed_read /= 4;
	}

# get ref mapped
  my $ref_mapped;
	{
		my $in = $design->get_sample_ref_filtered_sorted_file($sample);
		$ref_mapped = `samtools view $WORK_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $ref_mapped;
	}

# get ref dedup
  my $ref_dedup;
	{
		my $in = $design->get_sample_ref_dedup_file($sample);
		$ref_dedup = `samtools view -F 0x400 $WORK_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $ref_dedup;
	}

# get vec mapped
	my $vec_mapped;
	{
		my $in = $design->get_sample_vec_sorted_file($sample);
		$vec_mapped = `samtools view $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $vec_mapped;
	}

# get ref_novec
  my $ref_novec;
	{
		my $in = $design->get_sample_ref_novec_file($sample);
		$ref_novec = `samtools view -F 0x400 $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $ref_novec;
	}

# get insert site info
	my ($site, $site_filtered) = (0, 0);
	{
		my $in = $design->get_sample_ref_insert_site($sample);
		$site = `cat $WORK_DIR/$in | wc -l`; chomp $site;

		$in = $design->get_sample_ref_insert_site_filtered($sample);
		$site_filtered = `cat $WORK_DIR/$in | wc -l`; chomp $site_filtered;
	}

# get insert site uniq info
	my ($site_uniq, $site_uniq_filtered) = (0, 0);
	{
		my $in = $design->get_sample_ref_insert_site_uniq($sample);
		$site_uniq = `cat $WORK_DIR/$in | wc -l`; chomp $site_uniq;

		$in = $design->get_sample_ref_insert_site_filtered_uniq($sample);
		$site_uniq_filtered = `cat $WORK_DIR/$in | wc -l`; chomp $site_uniq_filtered;
	}

# get peak info
  my ($peak_count, $peak_clone) = (0, 0);
	{
		my $in = $design->get_sample_ref_peak_track($sample);
		open(BED, "<$BASE_DIR/$in") || die "Unable to open $in: $!";
		while(my $line = <BED>) {
			next if($line =~ /^(?:#|track)/);
			chomp $line;
			$peak_count++;
			my $peak_name = (split(/\t/, $line))[3];
			my ($dedup_count) = $peak_name =~ /ReadCount=(\d+)/;
			$peak_clone += $dedup_count;
		}
		close(BED);
	}

# get on/off-target peak info
	my ($ontarget_count, $ontarget_clone, $offtarget_count, $offtarget_clone) = (0, 0, 0, 0);
	{
		my $on_in = $design->get_sample_ref_peak_track_ontarget($sample);
		my $off_in = $design->get_sample_ref_peak_track_offtarget($sample);

		open(ON, "<$BASE_DIR/$on_in") || die "Unable to open $on_in: $!";
		while(my $line = <ON>) {
			next if($line =~ /^(?:#|track)/);
			chomp $line;
			$ontarget_count++;
			my $peak_name = (split(/\t/, $line))[3];
			my ($dedup_count) = $peak_name =~ /ReadCount=(\d+)/;
			$ontarget_clone += $dedup_count;
		}
		close(ON);
		open(OFF, "<$BASE_DIR/$off_in") || die "Unable to open $off_in: $!";
		while(my $line = <OFF>) {
			next if($line =~ /^(?:#|track)/);
			chomp $line;
			$offtarget_count++;
			my $peak_name = (split(/\t/, $line))[3];
			my ($dedup_count) = $peak_name =~ /ReadCount=(\d+)/;
			$offtarget_clone += $dedup_count;
		}
		close(OFF);
	}

# get clone info
	my ($clone_count, $clone_loc_count) = (0, 0, 0);
	my %clone_loc_freq;
	{
		my $in = $design->get_sample_ref_clone_track($sample);
		open(BED, "<$BASE_DIR/$in") || die "Unable to open $in: $!";
		while(my $line = <BED>) {
			next if($line =~ /^(?:#|track)/);
			chomp $line;
			my ($clone_name) = (split(/\t/, $line))[3];
			my ($loc_count) = $clone_name =~ /LocCount=(\d+)/;
			$clone_count++;
			$clone_loc_count += $loc_count;
			$clone_loc_freq{$loc_count}++;
		}
		close(BED);
	}

# output
  print OUT "$sample\t$total_read\t$trimmed_read\t$ref_mapped\t$ref_dedup\t$vec_mapped\t$ref_novec\t",
	"$site\t$site_filtered\t$site_uniq\t$site_uniq_filtered\t",
	"$peak_count\t$peak_clone\t$ontarget_count\t$ontarget_clone\t$offtarget_count\t$offtarget_clone\t",
  "$clone_count\t$clone_loc_count\t", get_freq_str(%clone_loc_freq), "\n";
}

close(OUT);

# subroutine definitions
sub get_freq_str {
  return 'NA' if(!@_);
  my %freq = @_;
  return join(',', map { "$_:$freq{$_}" } sort {$a <=> $b} keys %freq);
}
