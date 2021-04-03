#!/usr/bin/env perl
package CONTRA;

use strict;
use Getopt::Long;
use File::Basename;
use YAML::Tiny;

sub callCNV {
    
 my $inDir   = shift;
 my $control = shift;
 my $bed     = shift;
 my $outDir  = shift;

 my @bams = glob ("$inDir/*.bam");

 # Create CONTRA output directory
 mkdir "$outDir/CONTRA";

# Open the config
my $yaml = YAML::Tiny->read( "$::dirname/modules/contra.optimal.yaml" );
 
# Get a reference to the first document
my $contraConfig = $yaml->[0];

# Get a hash of params
my %params = %{ $contraConfig };

my $pm = Parallel::ForkManager->new($::threads);

 foreach my $bam(@bams) {
    my $pid = $pm -> start() and next; 

	my $name = basename($bam);
	$name =~s/.bam//;

	my $cmd = "python $::Callers{CONTRA} -t $bed -s $bam -c $control -o $outDir/CONTRA/$name --pval $params{pval} -l --numBin 500 --minExon 10 $::devNull";
	system($cmd);
	#print "$cmd\n";

    my ($allCNV)   = glob ("$outDir/CONTRA/$name/table/*FILTERED.txt");
    my ($largeCNV) = glob ("$outDir/CONTRA/$name/table/*LargeDeletion.txt");

    #print "$allCNV\t$largeCNV\n"; 

	if ($allCNV) {

  	    $cmd = "cat $allCNV | awk '{ print \$4\"\t\"\$5\"\t\"\$6\"\t\"\$3\";\"\$14\";\"\$8}'| tail -n +2 -f -> $outDir/CONTRA/$name.all_cnvs.bed";
		system $cmd;

		$cmd = "cat $largeCNV | awk '{ print \$1\"\t\"\$5\"\t\"\$6\"\t\"\$1\":\"\$5\"-\"\$6\";\"\$7\";\"\$8}'| tail -n +2 -f -> $outDir/CONTRA/$name.large_cnvs.bed";
		system $cmd;
				
		my $str = `cat $outDir/CONTRA/$name.large_cnvs.bed | grep 'NA'`;
		chomp $str;

		if (!$str) {
			$cmd = "bedtools intersect -a $outDir/CONTRA/$name.all_cnvs.bed -b $outDir/CONTRA/$name.large_cnvs.bed -v > $outDir/CONTRA/$name.single_cnvs.bed";
            system $cmd;

            # Merging small and large CNVs that overlap
			$cmd = "cat $outDir/CONTRA/$name.large_cnvs.bed $outDir/CONTRA/$name.single_cnvs.bed | sort -V > $outDir/CONTRA/$name.CONTRA.bed";
            system $cmd;
		}
		else {
		    $cmd = "cat $outDir/CONTRA/$name.all_cnvs.bed | sort -V > $outDir/CONTRA/$name.CONTRA.bed";
            system $cmd;
		}
		unlink ("$outDir/CONTRA/$name.single_cnvs.bed");
		unlink ("$outDir/CONTRA/$name.all_cnvs.bed");
		unlink ("$outDir/CONTRA/$name.large_cnvs.bed");
	}
    $pm->finish;
  }
  $pm->wait_all_children; 
}
################
sub cleanData {
	my $inDir = shift;

	my @files = glob ("$inDir/*");
	foreach my $file (@files) {
		print "$file\n";
		my @thrash = glob("$file/table/*");
		unlink @thrash;
	}

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

		my ($sample, $chr, $start, $end, $CNV) = split(/\t/, $line);

        # Get sample results
		my ($contraResult)= grep ($_=~/CONTRA.bed/, glob("$outDir/CONTRA/$sample*"));

		# Write coordinates to a tmp bed file
		my $variant_bed = BenchmarkUtils::write_bed($chr, $start, $end, $outDir);

		# All variants per sample 
		my $all_variants_bed = BenchmarkUtils::write_bed2(\%::SampleCnv, $sample, $outDir);

		# Reformat variant files for intersecting purposes
		my $normalized_file = BenchmarkUtils::normalize_file($contraResult, $outDir, $bed);

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


return 1;