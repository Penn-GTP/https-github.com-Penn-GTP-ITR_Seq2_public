#!/bin/env perl
# Prepare bash script for annotating ITR-peaks
our $VERSION = 'v2.2.1';
our $ENV_FILE = 'set_annotate_env.sh';

use strict;
use warnings;
use File::Basename;
use lib dirname (__FILE__);
use ITRSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $track_script = 'get_peak_track.pl';
my $bedtools = 'bedtools';
my $anno_script = 'get_peak_anno.pl';
my $sample_stats_script = 'get_sample_stats.pl';
#my $clone_loc_script = 'show_clone_loc_distrib.R';
my $cmd = "$0 " . join(" ", @ARGV);

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new ITRSeqExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $BASE_DIR = $design->get_global_opt('BASE_DIR');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
#my $DEMUX_DIR = $design->get_global_opt('DEMUX_DIR');
my $WORK_DIR = $design->get_global_opt('WORK_DIR');
#my $UMI_LEN = $design->get_global_opt('UMI_LEN');

# check required directories
if(!(-e $BASE_DIR && -d $BASE_DIR)) {
	print STDERR "Error: BASE_DIR $BASE_DIR not exists\n";
	exit;
}

if(!(-e $SCRIPT_DIR && -d $SCRIPT_DIR)) {
	print STDERR "Error: SCRIPT_DIR $SCRIPT_DIR not exists\n";
	exit;
}

if(!(-e $WORK_DIR)) {
	print STDERR "Error: WORK_DIR $WORK_DIR not exists, creating now\n";
	exit;
}

open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
# write header
print OUT "#!$sh_path\n";
print OUT qq(# CMD:"$cmd"\n# VER:$VERSION\n);
# set env
print OUT "source $SCRIPT_DIR/$ENV_FILE\n\n";

foreach my $sample ($design->get_sample_names()) {
# prepare peak track cmd
  {
		my $in = $design->get_sample_ref_peak($sample);
		my $out = $design->get_sample_ref_peak_track($sample);

		my $cmd = "$SCRIPT_DIR/$track_script $BASE_DIR/$in $BASE_DIR/$out --name $sample";

		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: annotated ref map file already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare annotate peak cmd
  {
		my $in = $design->get_sample_ref_peak_track($sample);
		my $gff = $design->sample_opt($sample, 'ref_gff');
		my $out = $design->get_sample_ref_peak_anno($sample);
		my $opts = $design->sample_opt($sample, 'anno_opts');

		my $cmd = "$SCRIPT_DIR/$anno_script $gff $BASE_DIR/$in $BASE_DIR/$out $opts";

		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: peak annotation already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare clone track cmd
	{
		my $in = $design->get_sample_ref_clone($sample);
		my $out = $design->get_sample_ref_clone_track($sample);

		my $cmd = "$SCRIPT_DIR/$track_script $BASE_DIR/$in $BASE_DIR/$out --name $sample";

		if(!-e "$BASE_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare annotate clone cmd
	{
		my $in = $design->get_sample_ref_clone_track($sample);
		my $gff = $design->sample_opt($sample, 'ref_gff');
		my $out = $design->get_sample_ref_clone_anno($sample);
		my $opts = $design->sample_opt($sample, 'anno_opts');

		my $cmd = "$SCRIPT_DIR/$anno_script $gff $BASE_DIR/$in $BASE_DIR/$out $opts";

		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: peak annotation already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

	print OUT "\n";
}

# prepare sample stat cmd
{
	my $out = $design->get_exp_stats_file($infile);

	my $cmd = "$SCRIPT_DIR/$sample_stats_script $infile $BASE_DIR/$out";
	if(!(-e "$BASE_DIR/$out")) {
		print OUT "$cmd\n";
	}
	else {
		print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
		print OUT "# $cmd\n";
	}
}

close(OUT);
# change to exacutable
chmod 0750, $outfile;
