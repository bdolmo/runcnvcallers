#!/usr/bin/env perl
package VISCAP;

use strict;
use Getopt::Long;
use File::Basename;
use Sort::Key::Natural qw(natsort);
use YAML::Tiny;


#####################################
sub callCNV {

 my $inDir = shift;
 my $bed   = shift;
 my $genome= shift;
 my $outDir= shift;

 $outDir = "$outDir/VISCAP";
 mkdir $outDir;
 mkdir "$outDir/DoC";

 # Open the config
 my $yaml = YAML::Tiny->read( "$::dirname/modules/viscap.optimal.yaml" );
 
 # Get a reference to the first document
 my $viscapConfig = $yaml->[0];

 # Get a hash of params
 my %params = %{ $viscapConfig };

 createConfig(\%params);

 VISCAP::cleanData("$outDir/VISCAP", $bed);

 # Edit bed file to be compatible with viscap
 my $bedVisCap = createTargetRegionsVisCap( $bed );

 my $samplesFile;
 my @bams = glob "$inDir/*.bam";
 my $pm = Parallel::ForkManager->new($::threads);

 # Exctracting coverage with GATK
 foreach my $bam (@bams) {
	my $name = basename($bam);
	$name =~s/.bam//;

    my $RG = `samtools view $bam -H | grep 'RG'`;
    chomp $RG;
	my $addedRG = $bam;
	$addedRG =~s/.bam/.rg.bam/ if $bam !~/rg/;

    if (!$RG && !-e $addedRG) {

		print " INFO: Adding read-group into $bam\n";
        my $cmd = "$::samtools addreplacerg -r 'ID:$name-1' -r 'SM:$name' -r 'LB:$bed' -r 'PL:ILLUMINA' -o $addedRG -O BAM $bam";
		system $cmd;

		$cmd = "$::samtools index $addedRG";
		system $cmd;

		#rename bam

		unlink($bam);
		unlink("$bam.bai");
		rename $addedRG, $bam;

		$cmd = "$::samtools index $bam";
		system $cmd;

    }
	if (!-e "$outDir/DoC/$name.DoC") {
		print " INFO: $bam has read group\n";
		my $cmd = "$::samtools index $bam";
		system $cmd;

		$cmd = "java -jar $::Utils{GATK} -T DepthOfCoverage -R $genome -o $outDir/DoC/$name.DoC -I $bam -L $bed\n";
		system $cmd;
	}
	$samplesFile = "$outDir/DoC/$name.DoC.DATA.sample_interval_summary";

	if (!-e $samplesFile) {
		$samplesFile =~s/.DATA//;
	}

	# Remove non-consitent lines between sample_interval_summary and Viscap BED ROI
	# This happens because GATK merges overlapping coordinates
 	makeConsistentFilesVisCap( $samplesFile, $bedVisCap ); 
 }

 # Run VisCap
 my $cmd = "Rscript $::Callers{VISCAP} $outDir/DoC $outDir $::Utils{VISCAP_CONFIG}";
 #system $cmd;

 # Transform XLS results to BED format
 VisCap2Bed($outDir);
}
##############
sub cleanData {
	# Remove all resulting files
	my $outDir = shift;

	my @resultsBed = glob ("$outDir/*VisCap.bed");
	unlink @resultsBed;

	my @runs = glob ("$outDir/VISCAP_*");
	foreach my $run (@runs) {
		print "$run\n";
		my @data = glob ("$run/*");
		unlink @data;
	}
}

###################
sub createConfig {

	my $href = shift;
	my %params = %$href;

    my $baseBed   = dirname($::bed);
    my $viscapBed = basename($::bed);
    $viscapBed =~s/.bed/.viscap.bed/;

    open (CONFIG, ">", $::Utils{VISCAP_CONFIG}) || die " ERROR: Unable to open $::Utils{VISCAP_CONFIG}\n";
    print CONFIG "interval_list_dir           <- \"$baseBed\"\n";
    #explorer_file               <- "C:\\Windows\\explorer.exe" 
    # this is the extension of the coverage file from DepthOfCoverage
    #print CONFIG "cov_file_pattern            <- \".DoC.DATA.sample_interval_summary\$\"\n";
    print CONFIG "cov_file_pattern            <- \".DoC.sample_interval_summary\$\"\n";

    print CONFIG "cov_field                   <- \"_total_cvg\"\n";
    print CONFIG "interval_file_pattern       <- \"$viscapBed\$\"\n"; # interval filename"
    print CONFIG "ylimits                     <- c(-2, 2)\n"; # interval filename"
    print CONFIG "iqr_multiplier              <- 3\n";
    print CONFIG "threshold.min_exons         <- 1\n";
    print CONFIG "threshold.cnv_log2_cutoffs  <- c($params{lowerCutoff}, $params{upperCutoff})\n";
    print CONFIG "iterative.calling.limit     <- 0\n";
    print CONFIG "infer.batch.for.sub.out_dir <- TRUE\n";
    print CONFIG "clobber.output.directory    <- FALSE\n";
    print CONFIG "dev_dir     <- FALSE\n";
    close CONFIG;
}

#####################################
sub createTargetRegionsVisCap {

	my $bed = shift;
	my $bedVisCap = $bed;
	$bedVisCap =~s/.bed/.viscap.bed/;
	open (BED, "<", $bed ) || die " ERROR: Unable to open $bed\n";
	open (VISCAP, ">", $bedVisCap ) || die " ERROR: Unable to open $bedVisCap\n";

	while (my $line=<BED>) {
		chomp $line;
		my @tmp = split (/\t/, $line);
		my @info = split (/;/, $tmp[3]);
		$tmp[1] = $tmp[1]+1;
		print VISCAP "$tmp[0]\t$tmp[1]\t$tmp[2]\t$info[1]\n";
	}
	close BED;
	close VISCAP;

	return $bedVisCap;
}

##############################################
sub VisCap2Bed {
    
	my $outDir = shift;

	my @runDirs = glob ("$outDir/VISCAP_run*");
	@runDirs = grep ( -d $_, @runDirs );
	
	# Hash of samples
	my %HoS = ();

	my %Samples = ();

	foreach my $runDir ( @runDirs ) {

		my @xls = glob ("$runDir/*cnvs.xls");
		foreach my $file ( @xls ) {

			my $sample = basename ($file);
			$sample=~s/.cnvs.xls//;
			
			$Samples{$sample}++;

			my $csvfile = $file;
			$csvfile =~s/.xls/.csv/;

			if (!-e $csvfile) {
				# Converting a CSV file to TXT file
				my $cmd = "ssconvert $file $csvfile $::devNull";
				system $cmd;
			}

			open (IN, "<", $csvfile) || die " ERROR: Unable to open $csvfile\n";
			while (my $line=<IN>) {
				chomp $line;
				next if $line =~/CNV_id/;
				next if $line !~/Gain|Loss/;
				my @tmp = split (",", $line);

				# 5' coordinates
				my ($chr1, $start1, $end1) = split (/[:-]/, $tmp[3]);

				# 3' coordinates
				my ($chr2, $start2, $end2) = split (/[:-]/, $tmp[4]);
				my $svtype = lossGain2DelDup($tmp[2]);

				if ($start1 > $end2) {
					push @{ $HoS{$sample} }, "$chr1\t$end2\t$start1\t$tmp[5]\t$svtype\t$tmp[8]";
				}
				else {
					push @{ $HoS{$sample} }, "$chr1\t$start1\t$end2\t$tmp[5]\t$svtype\t$tmp[8]";
				}
				#CA17296.rg	5	0	chr1:45794978-45795110	chr1:45805891-45805927	MUTYH	MUTYH	16	-342	39	151	-0.75	0.4	96
			}
			close IN;
		}
	}
	foreach my $sample ( natsort keys %Samples) {
		open (OUT, ">", "$outDir/$sample.VisCap.bed");
		foreach my $call ( natsort @{$HoS{$sample}} ) {
    		print OUT "$call\n";
		}  
		close OUT;
	}
}

###############################
sub makeConsistentFilesVisCap { 

	my $file1  = shift;
	my $file2  = shift;

	my %hash1;
	my %hash2;

	open (IN, "<", $file1);
	while (my $line=<IN>) {
		chomp $line;
		my @tmp = split (/[\t:-]/, $line);
		next if @tmp < 3;

		my $coord = "$tmp[0]\t$tmp[1]\t$tmp[2]";
		$hash1{$coord}++;
	}
	close IN;

	open (IN, "<", $file2);
	while (my $line=<IN>) {
		chomp $line;
		my @tmp = split (/\t/, $line);
		my $coord = "$tmp[0]\t$tmp[1]\t$tmp[2]";
		$hash2{$coord}++;
	}
	close IN;

	my %merged;
	foreach my $line (keys %hash2) {
		if (!exists $hash1{$line} ) {
			$merged{$line}++;
		}
	}

	foreach my $line (keys %hash1) {
		if (!exists $hash2{$line} ) {
			$merged{$line}++;
		}
	}

	# Rewrite files
	my $tmpFile = $file1;
	$tmpFile =~s/.sample_interval_summary/.tmp.sample_interval_summary/;
	open (OUT, ">", $tmpFile);
	open (IN, "<", $file1);
	while ( my $line=<IN>) {
		chomp $line;
		if ($line =~/^Target/) {
			print OUT "$line\n";
			next;
		}
		my @tmp = split (/[\t:-]/, $line);
		next if @tmp < 3;
		my $coord = "$tmp[0]\t$tmp[1]\t$tmp[2]";
		next if exists $merged{$coord};
		print OUT "$line\n";
	}
	close IN;
	close OUT;
	unlink($file1);
	rename $tmpFile, $file1;

	$tmpFile = $file2;	
	$tmpFile =~s/.bed/.tmp.bed/;
	open (OUT, ">", $tmpFile);
	open (IN, "<", $file2);
	while ( my $line=<IN>) {
		chomp $line;
		my @tmp = split (/\t/, $line);
		my $coord = "$tmp[0]\t$tmp[1]\t$tmp[2]";
		next if exists $merged{$coord};
		print OUT "$line\n";
	}
	close IN;
	close OUT;
	unlink($file2);
	rename $tmpFile, $file2;
}

###################
sub lossGain2DelDup {
	my $lg = shift;
	if ($lg eq 'Loss') {
		return 'DEL';
	}
	if ($lg eq 'Gain') {
		return 'DUP';
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
        #print "$line\n";

		my ($sample, $chr, $start, $end, $CNV) = split(/\t/, $line);

        # Get sample results
		my ($viscapResult) = grep ($_=~/.VisCap.bed/, glob("$analysisDir/VISCAP/$sample*"));

        #print "$analysisDir\t$sample\t$cnvkitResult\n"; exit;

		# Write coordinates to a tmp bed file
		my $variant_bed = BenchmarkUtils::write_bed($chr, $start, $end, $outDir);

		# All variants per sample 
		my $all_variants_bed = BenchmarkUtils::write_bed2(\%::SampleCnv, $sample, $outDir);

		# Reformat variant files for intersecting purposes
		my $normalized_file = BenchmarkUtils::normalize_file($viscapResult, $outDir, $bed);

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