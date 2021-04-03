#!/usr/bin/env perl
package CNVKIT;

use strict;
use Getopt::Long;
use File::Basename;
use YAML::Tiny;
use Cwd;

#############
sub cleanData {
    my $inDir = shift;

    my @segmentedCalls = glob ("$inDir/*calls.cns");
    unlink @segmentedCalls;

    my @calls = glob ("$inDir/*CNVKIT.bed");
    unlink @calls;

    my @zscores = glob ("$inDir/*zscore.txt");
    unlink @zscores;
}

#############
sub callCNV {

    # Reference will be created dynamically
    my $inDir = shift;
    my $bed   = shift;
    my $genome= shift;
    my $outDir= shift;

    # Open the config
    #my $yaml = YAML::Tiny->read( "$::dirname/modules/cnvkit.optimal.yaml" );
    my $yaml = YAML::Tiny->read( "$::dirname/modules/cnvkit.yaml" );

    # Get a reference to the first document
    my $cnvkitConfig = $yaml->[0];

    # Get a hash of params
    my %params = %{ $cnvkitConfig };
    my $currentDir = cwd();

    $outDir = "$outDir/CNVKIT";
    mkdir $outDir;
    my $target_bed = basename($bed);
    $target_bed=~s/.bed/.target.bed/;    
    my $nameTargetBed = $currentDir . "/" . $target_bed;
    $target_bed = "$outDir/$target_bed";

    my $antitarget_bed = basename($bed);
    $antitarget_bed=~s/.bed/.antitarget.bed/;
    my $nameAntiTargetBed = $currentDir . "/" . $antitarget_bed;
    $antitarget_bed = "$outDir/$antitarget_bed";

    # Creating accessible genomic regions
    if (!-e  "$outDir/access.hg19.bed") {
        my $cmd = "python3 $::Callers{CNVKIT} access $genome -o $outDir/access.hg19.bed > /dev/null 2>&1";
        system $cmd;
    }

    # Creating bins
    my $cmd = "python3 $::Callers{CNVKIT} autobin $inDir/*.bam -t $bed -g $outDir/access.hg19.bed";
    system $cmd;

    $cmd = "mv $nameTargetBed $target_bed";
    system $cmd;

    $cmd = "mv $nameAntiTargetBed $antitarget_bed";
    system $cmd;

    my @bams = glob "$inDir/*.bam";
    my $pm = Parallel::ForkManager->new($::threads);

    # Exctracting coverage
    foreach my $bam (@bams) {
        #my $pid = $pm -> start() and next; 
        my $name = basename($bam);
        $name =~s/.bam//;

        $cmd = "python3 $::Callers{CNVKIT} coverage $bam $target_bed -o $outDir/$name.targetcoverage.cnn";
        system $cmd if !-e "$outDir/$name.targetcoverage.cnn";

        #$cmd = "python3 $::Callers{CNVKIT} coverage $bam $antitarget_bed -o $outDir/$name.antitargetcoverage.cnn";
        #system $cmd if !-e "$outDir/$name.antitargetcoverage.cnn";
        #$pm->finish;
    }
    #$pm->wait_all_children; 

    # Create a pooled reference for all samples
    foreach my $bam (@bams) {

        my $name = basename($bam);
        $name =~s/.bam//;
        
        # Restrict analysis to target regions only
        my @coverages = glob ("$outDir/*.targetcoverage.cnn");
        @coverages = grep ($_ !~/$name/, @coverages);

        my $samples_2ref = join (" ", @coverages);

        $cmd = "python3 $::Callers{CNVKIT} reference $samples_2ref --fasta $genome -o $outDir/$name.reference.cnn > /dev/null 2>&1";
        system $cmd if !-e "$outDir/$name.reference.cnn";

        $cmd = "python3 $::Callers{CNVKIT} fix $outDir/$name.targetcoverage.cnn $outDir/$name.targetcoverage.cnn $outDir/$name.reference.cnn -o $outDir/$name.cnr > /dev/null 2>&1";
        system $cmd if !-e  "$outDir/$name.cnr";

        $cmd = "python3 $::Callers{CNVKIT} segment $outDir/$name.cnr -o $outDir/$name.cns > /dev/null 2>&1";
        system $cmd if !-e "$outDir/$name.cns";

        $cmd = "python3 $::Callers{CNVKIT} segmetrics -s $outDir/$name.cns --sem $outDir/$name.cnr -o $outDir/$name.metrics.cns";
        system $cmd if !-e "$outDir/$name.metrics.cns";

        $cmd = "python3 $::Callers{CNVKIT} call $outDir/$name.metrics.cns --filter sem -o $outDir/$name.calls.cns "
         . " -m threshold -t=$params{lowerDelCutoff},$params{upperDelCutoff},$params{lowerDupCutoff},$params{upperDupCutoff}";
         print"$cmd\n";
        system $cmd;

        $cmd = "python3 $::Callers{CNVKIT} export bed  $outDir/$name.calls.cns -o $outDir/$name.CNVKIT.bed > /dev/null 2>&1";
        system $cmd;

    } 
    foreach my $bam (@bams) {
        my $pid = $pm -> start() and next; 

        my $name = basename($bam);
        $name =~s/.bam//;
        $cmd = "python3 $::Utils{CNVKIT_ZTEST} -a $params{alpha} $outDir/$name.cnr -t | grep -v 'ztest' > $outDir/$name.zscore.txt";
        system $cmd;

        if (-e "$outDir/$name.zscore.txt"){ 
            
            my $cmd = "cat $outDir/$name.zscore.txt $outDir/$name.CNVKIT.bed | sort -V | uniq  > $outDir/$name.tmp.CNVKIT.bed";
            system $cmd;

            unlink "$outDir/$name.CNVKIT.bed";
            rename "$outDir/$name.tmp.CNVKIT.bed", "$outDir/$name.CNVKIT.bed";

            addRatioInfo("$outDir/$name.CNVKIT.bed", "$outDir/$name.calls.cns");
        }
        $pm->finish;
    }
    $pm->wait_all_children; 
}


################
sub benchmark {
    my $analysisDir = shift;
    my $bed         = shift;
    my $known_vars  = shift;
    my $outDir      = shift;
    my $threshold   = shift;

    my $TP = 0;
    my $FP = 0;

	# Step2 Benchmark itself
	open IN, "<", $known_vars || die " ERROR: Unable to open $known_vars\n";

	while (my $line=<IN>) {
		chomp $line;

		my ($sample, $chr, $start, $end, $CNV) = split(/\t/, $line);

        # Get sample results
		my ($cnvkitResult) = grep ($_=~/.CNVKIT.bed/, glob("$analysisDir/CNVKIT/$sample*"));

		# Write coordinates to a tmp bed file
		my $variant_bed = BenchmarkUtils::write_bed($chr, $start, $end, $outDir);

		# All variants per sample 
		my $all_variants_bed = BenchmarkUtils::write_bed2(\%::SampleCnv, $sample, $outDir);

		# Reformat variant files for intersecting purposes
		my $normalized_file = BenchmarkUtils::normalize_file($cnvkitResult, $outDir, $bed, $threshold);

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
sub addRatioInfo {

    my $cnvkitBed = shift;
    my $cnvkitSeg = shift;

    my $newCnvkitBed = $cnvkitBed;
    $newCnvkitBed=~s/.bed/.tmp.bed/;
    open (OUT, ">", $newCnvkitBed) || die " ERROR: Unable to open $newCnvkitBed\n";
    open (IN, "<", $cnvkitBed) || die " ERROR: Unable to open $cnvkitBed\n";
    while (my $line=<IN>) {
        chomp $line;
        my @tmp = split("\t", $line);
        my $value;
        #chr19	1206912	51364624	CA17296.rg.calls	1
        #chr19	33792551	33792705	chr19_33792243_33793321;CEBPA	-5.48163	0.37013	0.28419	4.16671e-07
        if (scalar @tmp < 6) {
            #chr16	2097454	2098078	chr16_2097454_2098078;TSC2	-1.15206	0	106.694	4	2.68125
            my $result = `grep -P '$tmp[0]\t$tmp[1]\t$tmp[2]' $cnvkitSeg`;
            my @tmpResult = split("\t", $result);
            $value = $tmpResult[4];
        }
        else {
            $value = $tmp[4];
        }
        print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$value\n";
        #print OUT "$line\t$value\n";
    }
    close IN;
    close OUT;
    unlink $cnvkitBed;
    rename $newCnvkitBed, $cnvkitBed;
}
###################

return 1;