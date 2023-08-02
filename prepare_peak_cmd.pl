#!/bin/env perl
# Prepare sh script for getting ITR peaks from filtered alignments
our $VERSION = 'v2.2.2';
our $ENV_FILE = 'set_peak_env.sh';

use strict;
use warnings;
use Bio::SeqIO;
use File::Basename;
use lib dirname (__FILE__);
use ITRSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $samtools = 'samtools';
my $bedtools = 'bedtools';
#my $picard = 'picard.jar';
my $insert_script = 'get_ITR_insert_site.pl';
my $space_script = 'get_ITR_primer_softclip_space.pl';
my $peak_script = 'get_ITR_peak.pl';
my $clone_script = 'get_ITR_clone.pl';
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
my $INSERT_SIZE = $design->get_global_opt('INSERT_SIZE');
my $KEEP_UNPAIR = $design->get_global_opt('KEEP_UNPAIR');
my $KEEP_STRAND = $design->get_global_opt('KEEP_STRAND');
my $MAX_PEAK_DIST = $design->get_global_opt('MAX_PEAK_DIST');
my $MIN_SPACE = $design->get_global_opt('MIN_SPACE');
#my $ONTARGET_FLANK = $design->get_global_opt('ONTARGET_FLANK');

my $DEFAULT_MIN_CLONE_LOC = 2;

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
	print STDERR "Error: WORK_DIR $WORK_DIR not exists\n";
	exit;
}

open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
# write header
print OUT "#!$sh_path\n";
print OUT qq(# CMD:"$cmd"\n# VER:$VERSION\n);
# set env
print OUT "source $SCRIPT_DIR/$ENV_FILE\n\n";

foreach my $sample ($design->get_sample_names()) {
# prepare get ITR insert site cmd
	{
		my $in = $design->get_sample_ref_novec_file($sample);
		my $loc_out = $design->get_sample_ref_insert_site($sample); 
		my $pos_out = $design->get_sample_ref_softclip_pos($sample); 
		my $primer_file = $design->get_global_primer_fwd();
		my $clip_len = 0;
		my $seq_in = new Bio::SeqIO(-file => "<$BASE_DIR/$primer_file", -format => 'fasta');
		while(my $seq_obj = $seq_in->next_seq()) {
			if($clip_len == 0 || $seq_obj->length() < $clip_len) {
				$clip_len = $seq_obj->length();
			}
		}

		my $cmd = "$SCRIPT_DIR/$insert_script $BASE_DIR/$in $WORK_DIR/$loc_out $WORK_DIR/$pos_out --insert-size $INSERT_SIZE --min-softclip $clip_len";

		if(!(-e "$WORK_DIR/$loc_out" && -e "$WORK_DIR/$pos_out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$loc_out, $WORK_DIR/$pos_out exist, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare primer softclip space cmd
  {
		my $trim_in = $design->get_sample_ITRtrim_pos($sample);
		my $clip_in = $design->get_sample_ref_softclip_pos($sample);

		my $out = $design->get_sample_ref_primer_softclip_space($sample);
		my $cmd = "$SCRIPT_DIR/$space_script $WORK_DIR/$trim_in $WORK_DIR/$clip_in $WORK_DIR/$out";
		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare filter insert site cmd
  {
		my $in = $design->get_sample_ref_insert_site($sample);
		my $space = $design->get_sample_ref_primer_softclip_space($sample);
		my $out = $design->get_sample_ref_insert_site_filtered($sample);

		my $cmd = "awk 'BEGIN{FS=\"\\t\"} (NR == FNR && FNR > 1 && \$NF != \"NA\" && \$NF >= $MIN_SPACE) {ID[\$1]; next} (\$4 in ID)' $WORK_DIR/$space $WORK_DIR/$in > $WORK_DIR/$out";
		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare insert site uniq cmd
  {
		my $in = $design->get_sample_ref_insert_site($sample);
		my $out = $design->get_sample_ref_insert_site_uniq($sample);

		my $cmd = "if [ -s $WORK_DIR/$in ]; then $bedtools sort -i $WORK_DIR/$in | $bedtools merge -d -$INSERT_SIZE -c 4,5,6 -o collapse,sum,collapse -i - > $WORK_DIR/$out ; else cp $WORK_DIR/$in $WORK_DIR/$out ; fi;";

		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare insert site filtered uniq cmd
  {
		my $in = $design->get_sample_ref_insert_site_filtered($sample);
		my $out = $design->get_sample_ref_insert_site_filtered_uniq($sample);

		my $cmd = "if [ -s $WORK_DIR/$in ]; then $bedtools sort -i $WORK_DIR/$in | $bedtools merge -d -$INSERT_SIZE -c 4,5,6 -o collapse,sum,collapse -i - > $WORK_DIR/$out ; else cp $WORK_DIR/$in $WORK_DIR/$out ; fi;";

		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare insert site merged cmd
  {
		my $in = $design->get_sample_ref_insert_site_filtered_uniq($sample);
		my $out = $design->get_sample_ref_insert_site_merged($sample);

		my $cmd = "if [ -s $WORK_DIR/$in ]; then $bedtools merge -d $MAX_PEAK_DIST -c 4,5,6 -o collapse,sum,collapse -i $WORK_DIR/$in > $WORK_DIR/$out ; else cp $WORK_DIR/$in $WORK_DIR/$out ; fi;";

		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare call peak cmd
	{
		my $in = $design->get_sample_ref_insert_site_merged($sample);
		my $out = $design->get_sample_ref_peak($sample);
		my $cmd = "$SCRIPT_DIR/$peak_script $WORK_DIR/$in $WORK_DIR/$out --keep-strand $KEEP_STRAND";

		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare call clone cmd
  {
		my $site_in = $design->get_sample_ref_insert_site_filtered_uniq($sample);
		my $aln_in = $design->get_sample_ref_novec_file($sample);
		my $out = $design->get_sample_ref_clone($sample);
		my $min_clone_loc = $design->sample_opt($sample, 'min_clone_loc') ? $design->sample_opt($sample, 'min_clone_loc') : $DEFAULT_MIN_CLONE_LOC;
		my $cmd = "$SCRIPT_DIR/$clone_script $WORK_DIR/$site_in $BASE_DIR/$aln_in $WORK_DIR/$out --min-loc $min_clone_loc";

		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

	print OUT "\n";
}

close(OUT);
# change to exacutable
chmod 0750, $outfile;
