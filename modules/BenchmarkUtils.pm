#!/usr/bin/env perl
package BenchmarkUtils;

use strict;
use warnings;
use diagnostics;
use File::Basename;

##########################
sub normalize_file {
    # Input is a caller result file, output is a uniform file
	my $raw_file = shift;
	my $outDir   = shift;
	my $bed      = shift;
    my $threshold = shift; # For cnvkit


	my $normalized_file = basename($raw_file);

	$normalized_file =~s/(.bed|.txt)/.normalized.txt/;
	$normalized_file = "$outDir/$normalized_file";


	#chr16	2097454	2098078	CA17296.rg.calls	-1.15206

	my $cmd = "cat $raw_file | sort -V | uniq | cut -f 1,2,3 | grep -v 'START' | ";
	$cmd.= "awk '{ if (\$2 > \$3){ print \$1\"\t\"\$3\"\t\"\$2 } else {print \$0} }'| ";
	
	# Desplegar els exons que solapen la CNV
	$cmd .= "bedtools intersect -a $bed -b stdin -wa | sort -V | uniq > $normalized_file";
	system $cmd;

	return $normalized_file;
}

##########################
sub write_bed {

	my $chr    = shift;
	my $start  = shift;
	my $end    = shift;
	my $outDir = shift;

	my $outBed = "$outDir/$chr.$start.$end.bed";

	open (BED, ">",$outBed);
	print BED "$chr\t$start\t$end";
	close BED;

	return $outBed;
}

##########################
sub write_bed2 {

	my $href   = shift;
	my $sample = shift;
	my $outDir = shift;
	my %SampleCnv = %$href;

	my $outBed = "$outDir/$sample.bed";
	open (BED, ">",$outBed);

	foreach my $call (@{ $SampleCnv{$sample} }) {
		my ($chr, $start, $end, $info) = split(/\t/, $call);
		print BED "$chr\t$start\t$end\n";
	}
	close BED;

	return $outBed;
}


 ##########################
sub true_positives {	
	
	# Get numer of True positives rois that span the CNV call
	my $caller_file = shift;
	my $variant_bed = shift;
	my $roi_bed     = shift;
	my $do_roi      = shift;
	my $all_vars    = shift;
	my $outDir      = shift;

	my $tp = 0;
	my $str;

	my @TruePositives = ();
	if ($do_roi) {

		my $cmd = "bedtools intersect -a $caller_file -b $variant_bed -wa | sort -V | uniq > $outDir/tp.bed";
		system $cmd;

		open (TMP, "<", "$outDir/tp.bed");
		while (my $line=<TMP>) {
			chomp $line;
			next if $line eq "";
			push @TruePositives, $line;
			$tp++;
		}
		close TMP;
	}
	else {
		$tp = `bedtools intersect -a $roi_bed -b $caller_file | sort -V | uniq | wc -l`;
		chomp $tp;
	}
	my $tp_str = join("\n", @TruePositives);
	unlink "$outDir/tp.bed";

	return $tp, $tp_str;
}

##########################
sub false_positives {
	
	# Get numer of True positives rois that span the CNV call
	my $caller_file = shift;
	my $variant_bed = shift;
	my $roi_bed     = shift;
	my $do_roi      = shift;
	my $all_vars    = shift;
	my $outDir      = shift;

	my $str;
	my $fp = 0;

	my @FalsePositives = ();
	if ($do_roi) {

		#print "$all_vars\n";
		my $cmd = "bedtools intersect -a $caller_file -b $variant_bed -v | ";
		$cmd.= " bedtools intersect -a stdin -b $all_vars -v " if -s $all_vars;
		$cmd.= " | sort -V | uniq > $outDir/fp.bed";
		system $cmd;

		open (TMP, "<", "$outDir/fp.bed");
		while (my $line=<TMP>) {
			chomp $line;
			next if $line eq "";
			push @FalsePositives, $line;
			$fp++;
		}
		close TMP;
	}
	else {
		$fp = `bedtools intersect -a $roi_bed -b $caller_file -v | sort -V | uniq | wc -l`;
		chomp $fp;
	}
	my $fp_str = join("\n", @FalsePositives);
	unlink "$outDir/fp.bed";
	return $fp, $fp_str;

}




return 1;