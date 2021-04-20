#!/usr/bin/env perl
package GRAPES;

use strict;
use Getopt::Long;
use YAML::Tiny;
use File::Basename;

##############
sub callCNV {

    my $inDir  = shift;
    my $bed    = shift;
    my $genome = shift;
    my $outDir = shift;

    # Open the config
    #my $yaml = YAML::Tiny->read( "$::dirname/modules/grapes.optimal.yaml" );
    my $yaml = YAML::Tiny->read( "$::dirname/modules/grapes.yaml" );

    # Get a reference to the first document
    my $grapesConfig = $yaml->[0];

    # Get a hash of params
    my %params = %{ $grapesConfig };

    $outDir = "$outDir/GRAPES";

    GRAPES::cleanData("$outDir/GRAPES", $bed);

    my $cmd = "$::Callers{GRAPES} wes -all --noofftarget --lowerdupcutoff $params{lowerdupcutoff}"
    . " --lowerdelcutoff $params{lowerdelcutoff} --upperdelcutoff $params{upperdelcutoff}"
    . " --mincorr $params{mincorr} --minsvsize 15 --minzscore $params{minzscore} "
    . " --pooled $inDir --bed $bed --genome $genome --outdir $outDir -t 4";
    print "$cmd\n";
    system($cmd);

}
################
sub benchmark {
    my $analysisDir = shift;
    my $bed         = shift;
    my $known_vars  = shift;
    my $outDir      = shift;

    my $TP = 0;
    my $FP = 0;

	# Step2 Benchmark itself
	open IN, "<", $known_vars || die " ERROR: Unable to open $known_vars\n";

	while (my $line=<IN>) {
		chomp $line;
        #print "$line\n";

		my ($sample, $chr, $start, $end, $CNV) = split(/\t/, $line);

        # Get sample results
		my ($grapesResult) = grep ($_=~/.CNV.bed/, glob("$analysisDir/GRAPES/$sample*"));

		# Write coordinates to a tmp bed file
		my $variant_bed = BenchmarkUtils::write_bed($chr, $start, $end, $outDir);

		# All variants per sample 
		my $all_variants_bed = BenchmarkUtils::write_bed2(\%::SampleCnv, $sample, $outDir);

		# Reformat variant files for intersecting purposes
		my $normalized_file = BenchmarkUtils::normalize_file($grapesResult, $outDir, $bed);

		if ($chr ne 'None') {

			# Get true positives, return number of rois if roi analysis 
			my ( $true_positives, $calls_tp ) = BenchmarkUtils::true_positives($normalized_file, $variant_bed, $bed, 1, $all_variants_bed, $outDir);
			my @tpCalls = split (/\n/, $calls_tp);

			# Adding true-positives if found
			$TP += $true_positives;

			# Get false positives, return number of rois
			my ( $false_positives, $calls_fp ) = BenchmarkUtils::false_positives($normalized_file, $variant_bed, $bed, 1, $all_variants_bed, $outDir);
			my @fpCalls = split (/\n/, $calls_fp);

            $FP += $false_positives;
        }
        unlink $variant_bed;
        unlink $all_variants_bed;
        unlink $normalized_file;
    }

    close IN;

    my $recall = sprintf "%.3f", $TP/$::total_cnv_rois;
    my $precision = sprintf "%.3f", $TP/($TP+$FP);

    return $recall, $precision;
}

###################
sub cleanData {
    # Since GRAPES concatenates some results when running multiple times o
    # over the same directory, we will remove them
    my $outDir = shift;
    my $bed    = shift;

    my $name = basename($bed);
    $name =~s/.bed//;

    my @bed = glob ("$outDir/*.CNV.bed");
    unlink @bed;

    @bed = glob ("$outDir/ON_TARGET/*.CNV.bed");
    unlink @bed;

    my @segmented = glob ("$outDir/ON_TARGET/SEGMENT_DATA/*");
    unlink @segmented;

    @segmented = glob ("$outDir/SEGMENT_DATA/*");
    unlink @segmented;

    my @ratios = glob ("$outDir/ON_TARGET/*ratios.txt.gz");
    unlink @ratios;
    
    my ($Ratios) = glob ("$outDir/GRAPES.Ratios.bed.gz");
    unlink ($Ratios);

    my $db = "$::Callers{GRAPES}/db/sqlite/$name.db";
    unlink $db;

}
###################
sub eraseResults  {

}

return 1;