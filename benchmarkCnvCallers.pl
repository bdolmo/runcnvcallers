#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use List::Util qw/shuffle/;
use List::Util qw( min max );
use Sort::Key::Natural qw(natsort);
use Excel::Writer::XLSX;

my $dirname = dirname (__FILE__);

my $roiLevel  = "";
my $geneLevel = "";
my $known_variants;
my $breakmer_output;
my $cnvkit_output;
my $grapes_output;
my $convading_output;
my $contra_output;
my $viscap_output;
my $useMlpa;
my $bed;
my $mlpa_baseline;
my $inDir;
my $outDir;
my $bySvType;
my $programs;

 my $bedtools = `which bedtools`;
 chomp $bedtools;
 if (!$bedtools) {
	 print " ERROR: BEDtools was not found on PATH\n";
 } 


Help () if (@ARGV < 8 or !GetOptions(
	'indir|i=s' =>\$inDir, 
	'known=s'   =>\$known_variants,
    'programs=s'=>\$programs,
	'bed=s'     =>\$bed,
	#'useMlpa'=>\$useMlpa,
	'mlpa=s'    =>\$mlpa_baseline,
	'by_svtype' =>\$bySvType,
	'outdir|o=s'=>\$outDir,
	)
);	

 # Validated genes by MLPA
 my $MLPAresults = "$dirname/ICR96_DATA/CANCER_ICR96_MLPA_BASELINE.txt";
 if (!-e $MLPAresults) {
	print " ERROR: $MLPAresults does not exist\n";
	exit;
 }

 my @arrProg = ( "breakmer", "contra", "grapes", "convading", "cnvkit", "viscap", "decon", "panelcnmops" );
 my %ProgHash = map { $_ => 1 } @arrProg;

 # Hash of callers 
 my %HoC = ();


 # Check program inputs  
 if ($programs) {
     # Parse programs
     my @parsed = split(",", lc $programs);
     foreach my $input(@parsed) {
         if (!exists $ProgHash{$input}) {
             print " ERROR: program $input is not accepted\n"; 
             Help();
         }
		 if ($input eq 'breakmer') {
			$HoC{BREAKMER}{INDIR} = "$inDir/BREAKMER";
			if (!-e $HoC{BREAKMER}{INDIR}){
				print " ERROR: $HoC{BREAKMER}{INDIR}\n"; 
			}
		 }
		 if ($input eq 'contra') {
			$HoC{CONTRA}{INDIR} = "$inDir/CONTRA";
			if (!-e $HoC{CONTRA}{INDIR}){
				print " ERROR: $HoC{CONTRA}{INDIR}\n"; 
			}
		 }
		 if ($input eq 'grapes'){
			 $HoC{GRAPES}{INDIR} = "$inDir/GRAPES/ON_TARGET";
			if (!-e $HoC{GRAPES}{INDIR}){
				print " ERROR: $HoC{GRAPES}{INDIR}\n"; 
			}
		 }
		 if ($input eq 'convading'){
			$HoC{CONVADING}{INDIR} = "$inDir/CONVADING";
			if (!-e $HoC{CONVADING}{INDIR}){
				print " ERROR: $HoC{CONVADING}{INDIR}\n"; 
			}
		 } 
		 if ($input eq 'cnvkit'){
			$HoC{CNVKIT}{INDIR} = "$inDir/CNVKIT";
			if (!-e $HoC{CNVKIT}{INDIR}){
				print " ERROR: $HoC{CNVKIT}{INDIR}\n"; 
			}
		 } 
		 if ($input eq 'viscap') {
			$HoC{VISCAP}{INDIR} = "$inDir/VISCAP";
			if (!-e $HoC{VISCAP}{INDIR}){
				print " ERROR: $HoC{VISCAP}{INDIR}\n"; 
			}
		 }
		 if ($input eq 'decon'){
			$HoC{DECON}{INDIR} = "$inDir/DECON";
			if (!-e $HoC{DECON}{INDIR}){
				print " ERROR: $HoC{DECON}{INDIR}\n"; 
			}
		 }
		 if ($input eq 'panelcnmops'){
			$HoC{PANELCNMOPS}{INDIR} = "$inDir/PANELCNMOPS";
			if (!-e $HoC{PANELCNMOPS}{INDIR}){
				print " ERROR: $HoC{PANELCNMOPS}{INDIR}\n"; 
			}
		 } 
	 }
 }

# Hash that will store all CNVs per sample
my %SampleCnv = ();

mkdir $outDir;

my %TP = ();
my %FP = ();

main();

###########################
sub main {

	# Initialize caller TP , FP
	foreach my $caller (sort keys %HoC){
		$HoC{$caller}{FP} = 0;
		$HoC{$caller}{TP} = 0;
	}  

	# Find how many exons have been validated by MLPA across all the dataset
	my ($total_validated_genes, $total_validated_rois) 
		= getMlpaInfo($MLPAresults, $bed);

	# rewrite bed as a sub-bed file with only ROIs validated by MLPA
	$bed = getEffectiveRois($MLPAresults, $bed);

	# Find how many rois are affected by CNVs
	my ($total_cnv_rois, $total_del_rois, $total_dup_rois, $refSampleCNV) 
		= getSampleKnownCnvs($known_variants, $bed);

	# De-referencing array 
	my %sampleKnownCnv = %$refSampleCNV;

	# THis hash will be used to create an output XLSX 
	my %tpResults = ();
	my %fpResults = ();
	foreach my $sample ( natsort keys %sampleKnownCnv ) {
		foreach my $call (@{ $sampleKnownCnv{$sample} }) {
			next if $call eq 'None';
			foreach my $caller ( natsort keys %HoC ){
				$tpResults{$sample}{$call}{$caller}  = "NO";  
			}	
		} 
	} 

	# Print summary validated information 
	print " INFO: MLPA_genes:$total_validated_genes\t"
		. "MLPA_ROIs:$total_validated_rois\tAffected_ROIs:$total_cnv_rois "
		. "(Deletions=$total_del_rois, Duplications=$total_dup_rois)\n";

	foreach my $caller ( natsort keys %HoC ){
		my @cnvFiles = getCnvsFromCaller($caller, $HoC{$caller}{INDIR});
		foreach my $cnvf (natsort @cnvFiles) {
			my @tmpFile = split(/\./, basename($cnvf));
			my $sample  = $tmpFile[0]; 
			$HoC{$caller}{$sample} = $cnvf; 
		}
	}

	foreach my $sample ( natsort keys %sampleKnownCnv ) {
		foreach my $caller ( natsort keys %HoC ){
			if (!exists $HoC{$caller}{$sample}){
				print "$sample\n";
				next;
			}
			if (-e $HoC{$caller}{$sample}){
				my ($tp, $fp) = doBenchmarkt($sample, $caller, $HoC{$caller}{$sample}, $outDir, \%tpResults, \%fpResults );
				$HoC{$caller}{TP}+= $tp;
				$HoC{$caller}{FP}+= $fp;
			} 
		} 
	} 
	my %HoW = ();
    my $workbook  = Excel::Writer::XLSX->new( $outDir . "/" . "Benchmarkt.xlsx" );
	foreach my $caller ( natsort keys %HoC ){
		$HoW{$caller} = $workbook->add_worksheet($caller);
        my $row = 0;
        $HoW{$caller}->write( $row, 0, 'Chromosome');
        $HoW{$caller}->write( $row, 1, 'Start'); 
        $HoW{$caller}->write( $row, 2, 'End'); 
        $HoW{$caller}->write( $row, 3, 'Exon'); 
        $HoW{$caller}->write( $row, 4, 'Sample'); 
        $HoW{$caller}->write( $row, 5, 'Detected'); 
        $HoW{$caller}->write( $row, 7, 'Chromosome');
        $HoW{$caller}->write( $row, 8, 'Start'); 
        $HoW{$caller}->write( $row, 9, 'End'); 
        $HoW{$caller}->write( $row, 10, 'Exon'); 
        $HoW{$caller}->write( $row, 11, 'Sample'); 
        $HoW{$caller}->write( $row, 12, 'Class'); 

		$row++;

   		foreach my $sample ( natsort keys %sampleKnownCnv ) {

			foreach my $call ( natsort keys %{$tpResults{$sample}} ){ 

				my ($chr, $start, $end, $exon) = split("\t", $call);
				my $class = $tpResults{$sample}{$call}{$caller};
				next if !$class;
				$HoW{$caller}->write( $row, 0, $chr);
				$HoW{$caller}->write( $row, 1, $start); 
				$HoW{$caller}->write( $row, 2, $end); 
				$HoW{$caller}->write( $row, 3, $exon); 
				$HoW{$caller}->write( $row, 4, $sample); 
				$HoW{$caller}->write( $row, 5, $class); 
            	$row++;
        	}
		}
		$row = 1;
   		foreach my $sample ( natsort keys %sampleKnownCnv ) {

			foreach my $call ( natsort keys %{$fpResults{$sample}{$caller} } ){ 

				my ($chr, $start, $end, $exon) = split("\t", $call);
				my $class = "FP";

				$HoW{$caller}->write( $row, 7, $chr);
				$HoW{$caller}->write( $row, 8, $start); 
				$HoW{$caller}->write( $row, 9, $end); 
				$HoW{$caller}->write( $row, 10, $exon); 
				$HoW{$caller}->write( $row, 11, $sample); 
				$HoW{$caller}->write( $row, 12, $class); 
            	$row++;
        	}
		}		

    }
    $workbook->close();

	open DATA, ">", "$outDir/data.txt";
	print DATA "Caller\tMetric\tValue\n";

	foreach my $caller ( natsort keys %HoC ) {
		my $precision = sprintf "%.2f", 100*($HoC{$caller}{TP}/($HoC{$caller}{TP}+$HoC{$caller}{FP}));
		$HoC{$caller}{PRECISION}   = $precision;
		$HoC{$caller}{SENSITIVITY} = sprintf "%.2f", 100*($HoC{$caller}{TP}/$total_cnv_rois);
		print " INFO: $caller\tTP:$HoC{$caller}{TP}\tFP:$HoC{$caller}{FP}\tPrecision:$precision\n"; 
		print DATA "$caller\tRecall\t$HoC{$caller}{SENSITIVITY}\n";
		print DATA "$caller\tPrecision\t$HoC{$caller}{PRECISION}\n";
	}
	close DATA;

	plotBench("$outDir/data.txt");
}

##############################
sub plotBench {
	my $data = shift;
	open (R, ">", "$outDir/benchplot.R");
	print R "library(ggplot2)\n";
	print R "library(ggpattern)\n";
	print R "library(ggsci)\n";
	print R "mydata<-read.table(file=\"$data\", sep=\"\t\", header=TRUE)\n";
	print R "attach(mydata)\n";
	print R "png(\"$outDir/out.png\", width=1400, height=800, res=200)\n";
	print R "yvalues<-seq(from=0, to=100, by=20)\n"; 
	print R "myplot<-ggplot(mydata, aes(x=Caller,y=Value,fill=Caller, alpha=Metric, pattern=Metric))+geom_bar_pattern(color = \"black\", pattern_fill = \"black\",pattern_angle = 45,pattern_density = 0.05,pattern_spacing = 0.025,pattern_key_scale_factor = 0.6
, stat=\"identity\", position=position_dodge()) + scale_pattern_manual(values = c(Recall = \"none\", Precision = \"stripe\")) + scale_fill_npg()\n";
	print R "myplot+theme_bw()+scale_alpha_manual(values = c(Recall = 1, Precision = 0.4))\n";
	print R "dev.off()\n"; 
	close R;
	`Rscript $outDir/benchplot.R`;
} 


##############################
sub doBenchmarkt {

	my $sample          = shift;
	my $caller          = shift;
	my $caller_results  = shift;
	my $outDir          = shift;
	my $tpResults       = shift;
	my $fpResults       = shift;

	my $tp = 0;
	my $fp = 0;

	# Sample validated ROIs (Normal + CNV) from MLPA
	my $validated_rois = "$outDir/$sample.all_rois.bed";
	my $cmd = "grep '$sample' $MLPAresults | cut -f 1,2,3 | $bedtools intersect -a $bed -b stdin > $validated_rois";
	system $cmd;

	# Normal ROIs from MPLA
	my $normal_rois = "$outDir/$sample.normal_rois.bed";
	$cmd = "grep '$sample' $MLPAresults | grep 'Normal' | cut -f 1,2,3 | $bedtools intersect -a $bed -b stdin > $normal_rois";
	system $cmd;

    # CNV ROIs from MPLA
	my $cnv_rois    = "$outDir/$sample.cnv_rois.bed";
	$cmd = "grep '$sample' $MLPAresults | grep 'Exon' | cut -f 1,2,3 | $bedtools intersect -a $bed -b stdin > $cnv_rois";
	system $cmd;

	# Total number of affected ROIs from MLPA
	my $n_cnv_rois = `cat $cnv_rois | wc -l`;
	chomp $n_cnv_rois;

	my $keepSegmentCalls = "";
	if ($caller eq 'CNVKIT') {
		#$keepSegmentCalls = " | awk '{ if (\$0 ~/calls/){ print \$0} }' ";
		$keepSegmentCalls = "";
	} 

	if ($n_cnv_rois > 0) {

		my $normalized_calls = "$outDir/$sample.normalized.calls.bed";
		my $cmd = "cat $caller_results $keepSegmentCalls | $bedtools intersect -a $validated_rois -b stdin -wa | sort -V | uniq  > $normalized_calls";
		system $cmd;

		my $tpStr = `$bedtools intersect -a $normalized_calls -b $cnv_rois -wa | sort -V | uniq | wc -l`;
		chomp $tpStr;
		$tp+=$tpStr;

		open (TP, "$bedtools intersect -a $normalized_calls -b $cnv_rois -wa | sort -V | uniq |");
		while (my $call=<TP>) {
			chomp $call;
			$tpResults->{$sample}->{$call}->{$caller} = "YES";
		} 
		close TP;

		my $fpStr = `$bedtools intersect -a $normalized_calls -b $normal_rois -wa | sort -V | uniq | wc -l`;
		chomp $fpStr;
		$fp+=$fpStr;
		open (FP, "$bedtools intersect -a $normalized_calls -b $normal_rois -wa | sort -V | uniq |" );
		while (my $call=<FP>) {
			chomp $call;
			$fpResults->{$sample}->{$caller}->{$call}++;
		} 
		close FP;		
	}
	else {

		my $normalized_calls = "$outDir/$sample.normalized.calls.bed";
		my $cmd = "cat $caller_results $keepSegmentCalls | $bedtools intersect -a $validated_rois -b stdin -wa | sort -V | uniq > $normalized_calls";
		system $cmd;

		my $fpStr = `$bedtools intersect -a $normalized_calls -b $normal_rois  -wa | sort -V | uniq | wc -l`;
		chomp $fpStr;
		$fp+=$fpStr;	
		open (FP, "$bedtools intersect -a $normalized_calls -b $normal_rois -wa | sort -V | uniq |" );
		while (my $call=<FP>) {
			chomp $call;
			$fpResults->{$sample}->{$caller}->{$call}++;
		} 
		close FP;
	} 

	return ($tp, $fp);
} 

##############################
sub getCnvsFromCaller  {

	my $caller = shift;
	my $inDir  = shift;

	my @cnvResults = ();
	if ($caller eq "GRAPES"){
		@cnvResults = glob("$inDir/*.CNV.bed");  
	}
	if ($caller eq "DECON"){
		@cnvResults = glob("$inDir/*.DECON.bed");  
	} 
	if ($caller eq "PANELCNMOPS"){
		@cnvResults = glob("$inDir/*.PANELCNMOPS.bed");  
	}
	if ($caller eq "CNVKIT"){
		@cnvResults = glob("$inDir/*.CNVKIT.bed");  
	}
	if ($caller eq "CONVADING"){
		@cnvResults = glob("$inDir/*.CONVADING.bed");  
	}
	if ($caller eq "VISCAP"){
		@cnvResults = glob("$inDir/*.VisCap.bed");  
	}
	return @cnvResults; 
} 

##############################
sub selectValidHitsConvading {

	my $input = shift;
	my $output = $input;

	$output =~s/.totallist.txt/.filtered.totallist.txt/;

	open (OUTFILE, ">", $output)  || die " ERROR: Unable to open $output\n";
	open (INFILE, "<", $input) || die " ERROR: Unable to open $input\n";
	while (my $line=<INFILE>) {
		chomp $line;
		my @tmp = split (/\t/, $line);
		next if $line =~/LOW_QUALITY/;
		if ($line !~/DUP/ && $line !~/DEL/) {
			next;
		}
		print OUTFILE "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t2\t2\t$tmp[10]\t$tmp[4]\n";
		#CHR	START	STOP	GENE	AUTO_RATIO	AUTO_ZSCORE	AUTO_VC	GENE_RATIO	GENE_ZSCORE	GENE_VC	ABBERATION	QUALITY	SHAPIRO-WILK
		#chr8	145742433	145742575	RECQL4	1.0298847443299	0.346739533956389	0.0861878770756533	0.998158829493921	-0.0318553529542735	0.0577978372652784	DUP	.	0.92626
	}
	close INFILE;
	close OUTFILE;
	return $output;
}
##########################
sub updateBenchmark {

	# This module actually does the benchmarkt
	my $chr         = shift;
	my $start       = shift;
	my $end         = shift;
	my $svtype      = shift;
	my $href        = shift;
	my $sample      = shift;
	my $bed         = shift;
	my $MLPAresults = shift;
	my $bySvType    = shift;
	my $del_threshold = shift;

	# Skipping sample without CNVs
	if ($chr eq 'None') {
		return;
	}

	# De-reference hash of files
	my %HoFiles = %$href;

	# Write coordinates to a tmp bed file
	my $variant_bed = write_bed($chr, $start, $end);

	foreach my $caller ( sort keys %HoFiles) {

		# All variants per sample 
		my $all_variants_bed = write_bed2(\%SampleCnv, $sample);

		# Reformat variant files for intersecting purposes
		my $normalized_file = normalize_file($HoFiles{$caller}, $svtype, $caller, $del_threshold, $del_threshold);

		# Get true positives
		my ( $true_positives, $calls_tp ) = true_positives($normalized_file, $variant_bed, $bed, 1, $all_variants_bed);
		
		my @tpCalls = split (/\n/, $calls_tp);
		foreach my $tpcall(@tpCalls) {
			my($chrom, $st, $ed, $ratioValue)=split(/\t/, $tpcall);
			#print ALL_CALLS "$sample\t$chrom\t$st\t$ed\t$svtype\t$caller\t$ratioValue\tTP\n";
		}		

		# Summing true-positives
		$TP{$caller} += $true_positives;

		# Get false positives, return number of rois
		my ( $false_positives, $calls_fp ) = false_positives($normalized_file, $variant_bed, $bed, 1, $all_variants_bed);

		my @fpCalls = split (/\n/, $calls_fp);
		foreach my $fpcall(@fpCalls) {			
			my($chrom, $st, $ed, $ratioValue) = split(/\t/, $fpcall);
			#print ALL_CALLS "$sample\t$chrom\t$st\t$ed\t$svtype\t$caller\t$ratioValue\tFP\n";
		}

		# Adding false-positives if found
		$FP{$caller} += $false_positives;
		if ( $true_positives > 0 ) {
			print GLOBAL_FN "\tYES";
		}
		else {
			print GLOBAL_FN "\tNO";
		}
		
		unlink($all_variants_bed);
		unlink($normalized_file);
	}		
	unlink($variant_bed);
}
 ##########################
sub true_positives {	
	
	# Get numer of True positives rois that span the CNV call
	my $caller_file = shift;
	my $variant_bed = shift;
	my $roi_bed     = shift;
	my $do_roi      = shift;
	my $all_vars    = shift;

	my $tp = 0;
	my $str;

	my @TruePositives = ();

	if ($do_roi) {

		my $cmd = "bedtools intersect -a $caller_file -b $variant_bed -wa | sort -V | uniq > $outDir/tp.bed";
		system $cmd;

		open (TMP, "<", "$outDir/tp.bed");
		my %seen = ();
		while (my $line=<TMP>) {
			chomp $line;
			next if $line eq "";
			my @tmp = split (/\t/, $line);
			$seen{"$tmp[0]\t$tmp[1]\t$tmp[2]"}++;
			if ($seen{"$tmp[0]\t$tmp[1]\t$tmp[2]"} > 1) {
				next;
			}

			push @TruePositives, $line;
			$tp++;
		}
		close TMP;
		unlink("$outDir/tp.bed");
	}
	else {
		$tp = `bedtools intersect -a $bed -b $caller_file | sort -V | uniq | wc -l`;
		chomp $tp;
	}
	my $tp_str = join("\n", @TruePositives);
	unlink("$outDir/tp.bed");

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

	my $str;
	my $fp = 0;

	my @FalsePositives = ();
	if ($do_roi) {

		my $cmd = "bedtools intersect -a $caller_file -b $variant_bed -v | ";
		$cmd.= " bedtools intersect -a stdin -b $all_vars -v " if -s $all_vars;
		$cmd.= " | sort -V | uniq > $outDir/fp.bed";
		system $cmd;
		my %seen = ();

		open (TMP, "<", "$outDir/fp.bed");
		while (my $line=<TMP>) {
			chomp $line;
			next if $line eq "";
			my @tmp = split (/\t/, $line);

			$seen{"$tmp[0]\t$tmp[1]\t$tmp[2]"}++;
			if ($seen{"$tmp[0]\t$tmp[1]\t$tmp[2]"} > 1) {
				next;
			}

			push @FalsePositives, $line;
			$fp++;
		}
		close TMP;
	}
	else {
		$fp = `bedtools intersect -a $bed -b $caller_file -v | sort -V | uniq | wc -l`;
		chomp $fp;
	}
	my $fp_str = join("\n", @FalsePositives);
	unlink("$outDir/fp.bed");

	return $fp, $fp_str;
}
##########################
sub normalize_file {

	my $raw_file = shift;
	my $svtype   = shift,
	my $algorithm = shift;
	my $del_threshold = shift;
	my $dup_threshold = shift;

	my $addedRatioFile = $raw_file;
	$addedRatioFile =~s/.bed/.added.ratio.bed/ if $raw_file =~/.bed/;
	$addedRatioFile =~s/.txt/.added.ratio.txt/ if $raw_file =~/.txt/;

	my $ratioValue;
	open (ADDED, ">", $addedRatioFile);
	open (INPUT, "<", $raw_file);

	while (my $line=<INPUT>) {
		chomp $line;

		next if $line =~/START/;
		next if $svtype eq 'DUP';
		my @tmp = split(/\t/, $line);
		my @info = split(";", $tmp[3]);

		if ($algorithm eq 'CONTRA') {
			if ($svtype eq 'DEL') {
				next if $line !~/loss/;
			}
			if ($svtype eq 'DUP') {
				next if $line !~/gain/;
			}
			$ratioValue = $info[-1];			
		}
		if ($algorithm eq 'CNVKIT') {
			$ratioValue = $tmp[4];
			my $cnvkit_svtype;
			if ($ratioValue > 0) {
				$cnvkit_svtype = "DUP";
			}
			else {
				$cnvkit_svtype = "DEL";
			}
			next if $cnvkit_svtype ne $svtype;
		}
		if ($algorithm eq 'GRAPES') {
			next if $line !~/$svtype/;

			($ratioValue) = grep ($_=~/RRD=/, @info);
			$ratioValue =~s/RRD=//;
			$ratioValue = 0.54 if $ratioValue eq '.';
		}
		if ($algorithm eq 'CONVADING') {
			next if $line !~/$svtype/;
			$ratioValue = $tmp[-1];
		}
		if ($algorithm eq 'VISCAP') {
			next if $line !~/$svtype/;
			$ratioValue = $tmp[5];
		}
		if ($svtype eq 'DEL') {
			if ($algorithm eq 'CONTRA' or $algorithm eq 'CNVKIT' or $algorithm eq 'VISCAP')  {
				if ( abs($ratioValue) > 100 && abs($ratioValue) < 1000 ) {
					$ratioValue = $ratioValue/1000;
				}
				elsif ( abs($ratioValue) >= 1000) {
					$ratioValue = $ratioValue/1000;
				}
				my $threshold = log2($del_threshold);
				#print"$algorithm\t$ratioValue\t$threshold\n";
				next if $ratioValue > $threshold;
			}
			else {
				next if $ratioValue > $del_threshold;
			}
		}
		if ($svtype eq 'DUP') {
			next if $ratioValue > $dup_threshold;
		}

		print ADDED "$tmp[0]\t$tmp[1]\t$tmp[2]\t$ratioValue\n";
	}
	close INPUT;
	close ADDED;
	my $normalized_file = basename($raw_file);

	$normalized_file =~s/.bed/.normalized.txt/ if $normalized_file =~/.bed/;
	$normalized_file =~s/.txt/.normalized.txt/ if $normalized_file =~/.txt/;

	$normalized_file = "$outDir/$normalized_file";

	my $cmd = "cat $addedRatioFile | sort -V | uniq | cut -f 1,2,3,4 | grep -v 'START' | ";
	$cmd.= "awk '{ if (\$2 > \$3){ print \$1\"\t\"\$3\"\t\"\$2\"\t\"\$4 } else {print \$0} }'| ";
	
	# Desplegar els exons que solapen la CNV
	$cmd .= "bedtools intersect -a $bed -b stdin -wa -wb | sort -V | uniq | cut -f 1,2,3,8 > $normalized_file";
	system $cmd;

	unlink($addedRatioFile);
	return $normalized_file;
}

##########################
sub write_bed {

	my $chr = shift;
	my $start = shift;
	my $end = shift;

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
	my %SampleCnv = %$href;

	my $outBed = "$outDir/$sample.bed";
	open (BED, ">",$outBed);

	foreach my $call (@{ $SampleCnv{$sample} }) {
		my ($chr, $start, $end, $info) = split(/\t/, $call);
		print BED "$chr\t$start\t$end\t$info\n";
	}
	close BED;
	return $outBed;
}

##########################
sub getEffectiveRois {

	my $MLPAresults = shift;
	my $rawBed      = shift;

	my $effective_bed = basename($rawBed);
	$effective_bed =~s/.bed/.validated.bed/;
	$effective_bed = "$outDir/$effective_bed";

	my $cmd = "bedtools intersect -a $rawBed -b $MLPAresults -wa";
	$cmd.= " | sort -V | uniq  | grep -v 'rs'> $effective_bed";
	system $cmd;
	
	return $effective_bed;
}
###########################
sub getMlpaInfo {

	my $MLPAresults = shift;
	my $ROIbed      = shift;

	my $validated_rois = 0;
	my $num_genes = `cat $MLPAresults | wc -l`;
	chomp $num_genes;

	open IN, "$bedtools intersect -a $ROIbed -b $MLPAresults -wa -wb |"
		|| die " ERROR: Cannot intersect\n";
	while (my $line =<IN>) {
		chomp $line;
		$validated_rois++;
	}
	close IN;
	return $num_genes, $validated_rois;
}
############################
sub getSampleKnownCnvs {

	my $known_variants = shift;
	my $ROIbed      = shift;

	my $num_cnv_rois = 0;
	my $num_del_rois = 0;
	my $num_dup_rois = 0;

	my %sampleCNV = ();

	open ( IN, "<", $known_variants) || die " ERROR: Unable to open $known_variants\n";
	while (my $line=<IN>) {

		chomp $line;
		my ($sample, $chr, $start, $end, $CNV) = split(/\t/, $line);
		my $call = $line =~/None/ ? "None" : "$chr\t$start\t$end"; 
		#print "$line\n";
		#push @{ $sampleCNV{$sample} }, $call;
		if ($line =~/None/) {
			push @{ $sampleCNV{$sample} }, $call;
			next;
		} 

		my $tmpBed = write_bed($chr, $start, $end);
		my $x = `$bedtools intersect -a $ROIbed -b $tmpBed | wc -l`;
		chomp $x;

		my $exons = `$bedtools intersect -a $ROIbed -b $tmpBed`;
		chomp $exons;

		my @tmpExons = split(/\n/, $exons);
		foreach my $exon(@tmpExons) {
			push @{ $sampleCNV{$sample} }, $exon;
		}

		if ($line =~/deletion/i) {
			$num_del_rois+=$x;
		}
		if ($line =~/duplication/i) {
			$num_dup_rois+=$x;
		}

		$num_cnv_rois+=$x;
		unlink($tmpBed);
	}
	close IN;
	
	return $num_cnv_rois, $num_del_rois, $num_dup_rois, \%sampleCNV;

}
############################
sub shuffle_array {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}
############################
sub log2 {
    my $n = shift;
    return log($n)/log(2);
}

###########################
sub Help {

  print "\n Usage: $0 [PARAMS]

 Params:  
  -i,--indir  STRING Input directory with caller results
  -o,--outdir      STRING 	 Output directory
  --programs    STRING   Programs separated by comma: breakmer,contra,convading,cnvkit,decon,grapes,viscap
  --bed		STRING	 BED regions file
  --mlpa	STRING	 mlpa validation file
  --known	STRING	 ICR96 validated CNVs
  --useMlpa     Use validated exons only\n\n";
	exit;
}