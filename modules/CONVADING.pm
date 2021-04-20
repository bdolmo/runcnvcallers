#!/usr/bin/env perl
package CONVADING;

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use YAML::Tiny;


sub callCNV {
 # Subroutine that does all calling commands

 my $inDir  = shift;
 my $bed    = shift;
 my $genome = shift;
 my $outDir = shift;
 my $controlDir = shift;

# Open the config
my $yaml = YAML::Tiny->read( "$::dirname/modules/convading.yaml" );
 
# optimized parameters
#my $yaml = YAML::Tiny->read("$::dirname/modules/convading.optimal.yaml");

# Get a reference to the first document
my $convadingConfig = $yaml->[0];

# Get a hash of params
my %params = %{ $convadingConfig };

 my $bed_name = basename($bed);
 my $convading_name = "$outDir/$bed_name";
 $convading_name =~s/.bed/.convading.bed/;
 
 # We need to modify initial bed file to fulfill convading requirements
 my $cmd = "cat $bed | awk '{ split(\$4, a, \";\"); gene=a[2]; print \$1\"\t\"\$2\"\t\"\$3\"\t\"gene }' > $convading_name";
 system $cmd;
  
 my $convadingOut = "$outDir/CONVADING";
 mkdir $convadingOut;

 # Extract counts and normalize
 my ($normalized) = glob ("$convadingOut/StartWithBam/*normalized.coverage.txt");
 if (!-e $normalized) {
	$cmd = "perl $::Callers{CONVADING} -mode StartWithBam -inputDir $inDir -controlsDir $convadingOut/controls" 
	." -useSampleAsControl -outputDir $convadingOut/StartWithBam -bed $convading_name";
	print "$cmd\n";
	system($cmd);
 }

 # Select most informative samples
 $cmd = "perl $::Callers{CONVADING} -mode StartWithMatchScore -inputDir $convadingOut/StartWithBam -controlsDir "
 . " $convadingOut/controls -outputDir $convadingOut/StartWithMatchScore";
 print "$cmd\n";
 system($cmd);
 
 # CNV detection
 $cmd = "perl $::Callers{CONVADING} -mode StartWithBestScore -regionThreshold $params{regionThreshold} -ratioCutOffLow $params{ratioCutOffLow} -ratioCutOffHigh $params{ratioCutOffHigh}" 
 . " -zScoreCutOffHigh $params{zScoreCutOffHigh} -zScoreCutOffLow $params{zScoreCutOffLow} -inputDir $convadingOut/StartWithMatchScore -outputDir $convadingOut/StartWithBestScore $::devNull";
 system($cmd); 
 print "$cmd\n";

 # GenerateTargetQcList
 $cmd = "perl $::Callers{CONVADING} -mode GenerateTargetQcList -inputDir $convadingOut/StartWithBam "
 . " -controlsDir $convadingOut/controls -outputDir $convadingOut $::devNull";
 print "$cmd\n";
 system($cmd);

 # CreateFinalList
 $cmd = "perl $::Callers{CONVADING} -mode CreateFinalList -percentageLessReliableTargets 20 -inputDir $convadingOut/StartWithBestScore "
 . " -targetQcList $convadingOut/targetQcList.txt -outputDir $convadingOut";
 print "$cmd\n";
 system($cmd);

 my @totalListFiles = glob("$convadingOut/StartWithBestScore/*totallist.txt");

 foreach my $tFile (@totalListFiles) {
	 my $convadingBed = $convadingOut . "/" . basename($tFile);
	 $convadingBed =~s/.best.score.totallist.txt/.CONVADING.bed/;
	 open (OUTPUT, ">", $convadingBed);
	 open (INPUT, "<", $tFile) || die " ERROR: Unable to open $tFile\n";
	 while (my $line=<INPUT>) {
		chomp $line;
		my @tmp = split (/\t/, $line);

		# If no DUP or DEL found in line, variant is filtered		
		if ($line !~/DUP/ && $line !~/DEL/) {
			next;
		}
		print OUTPUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t2\t2\t$tmp[10]\t$tmp[4]\n";
		#CHR	START	STOP	GENE	AUTO_RATIO	AUTO_ZSCORE	AUTO_VC	GENE_RATIO	GENE_ZSCORE	GENE_VC	ABBERATION	QUALITY	SHAPIRO-WILK
		#chr8	145742433	145742575	RECQL4	1.0298847443299	0.346739533956389	0.0861878770756533	0.998158829493921	-0.0318553529542735	0.0577978372652784	DUP	.	0.92626
		#chr8	145742433	145743169	RECQL4	4	2	DUP
	 }
	 close INPUT;
	 close OUTPUT;
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
        print "$line\n";

		my ($sample, $chr, $start, $end, $CNV) = split(/\t/, $line);

        # Get sample results
		my ($convadingResult)= grep ($_=~/best.score.totallist.txt/, glob("$outDir/CONVADING/StartWithBestScore/$sample*"));

		# Remove LOW QUAL cnv calls
		$convadingResult = selectValidHitsConvading($convadingResult);

		# Reformat variant files for intersecting purposes
		my $normalized_file = BenchmarkUtils::normalize_file($convadingResult, $outDir, $bed);

		# Write coordinates of the known CNV to a tmp bed file
		my $variant_bed = BenchmarkUtils::write_bed($chr, $start, $end, $outDir);

		# Write all validated CNV-exons for the sample being tested
		my $all_variants_bed = BenchmarkUtils::write_bed2(\%::SampleCnv, $sample, $outDir);

		if ($chr ne 'None') {

			# Get true positives, return number of rois if roi analysis 
			my ( $true_positives, $calls_tp ) = BenchmarkUtils::true_positives($normalized_file, $variant_bed, $bed, 1, $all_variants_bed, $outDir);
			my @tpCalls = split (/\n/, $calls_tp);

			foreach my $tpcall (@tpCalls) {
				print "$tpcall\n";
			}

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
	print "$::total_cnv_rois\t$TP\n";

    my $recall = sprintf "%.3f", $TP/$::total_cnv_rois;
    my $precision = sprintf "%.3f", $TP/($TP+$FP);

    return $recall, $precision;
}

##########################
sub selectValidHitsConvading {
	my $input = shift;
	my $output = $input;

	$output =~s/.totallist.txt/.filtered.totallist.txt/;
	open (OUTPUT, ">", $output);
	open (INPUT, "<", $input) || die " ERROR: Unable to open $input\n";
	while (my $line=<INPUT>) {
		chomp $line;
		my @tmp = split (/\t/, $line);

		# If no DUP or DEL found in line, variant is filtered		
		if ($line !~/DUP/ && $line !~/DEL/) {
			next;
		}
		print OUTPUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t2\t2\t$tmp[10]\t$tmp[4]\n";
		#CHR	START	STOP	GENE	AUTO_RATIO	AUTO_ZSCORE	AUTO_VC	GENE_RATIO	GENE_ZSCORE	GENE_VC	ABBERATION	QUALITY	SHAPIRO-WILK
		#chr8	145742433	145742575	RECQL4	1.0298847443299	0.346739533956389	0.0861878770756533	0.998158829493921	-0.0318553529542735	0.0577978372652784	DUP	.	0.92626
		#chr8	145742433	145743169	RECQL4	4	2	DUP
	}
	close INPUT;
	close OUTPUT;
	return $output;
}



return 1;