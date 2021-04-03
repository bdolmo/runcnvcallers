#!/usr/bin/env perl
package PANELCNMOPS;

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

 $outDir = "$outDir/PANELCNMOPS";
 mkdir $outDir;

# Open the config
 my $yaml = YAML::Tiny->read( "$::dirname/modules/panelcnmops.yaml" );
 
 # Get a reference to the first document
 my $panelcnmopsConfig = $yaml->[0];

 # Get a hash of params
 my %params = %{ $panelcnmopsConfig };

 my $panelcnR = "$outDir/panelcnmops.R";
 my $countWindowsRdata = "$outDir/countWindows.rds";
 my $countRdata = "$outDir/counts.rds";

 open (R, ">", $panelcnR) || die " ERROR: Unable to open $panelcnR\n";
 print R "library(panelcn.mops)\n";
 print R "library(plyr)\n";
 print R "readLength<-100\n";

 print R "auxCNname <- function(x) {
  if (x \%in\% c(\"CN0\", \"CN1\")) return(\"deletion\") 
  else if (x \%in\% c(\"CN3\", \"CN4\")) return(\"duplication\")
}\n";

 print R "allbams <- list.files(path=\"$inDir\", pattern=\"*.bam\$\", full.names = TRUE)\n";

 if (-e $countWindowsRdata) {
    print R "countWindows <- readRDS(\"$countWindowsRdata\")\n"; 
 }
 else {
    print R "countWindows <- getWindows(\"$bed\")\n";
    print R "saveRDS(countWindows, file.path(\"$outDir\", \"countWindows.rds\"))\n"; 
 } 
 if (-e $countRdata){
    print R "counts <- readRDS(\"$countRdata\")\n";
 }
 else{
    print R "counts <- countBamListInGRanges(countWindows = countWindows, bam.files = allbams, read.width = readLength)\n";
    print R "saveRDS(counts, file.path(\"$outDir\", \"counts.rds\"))\n"; 
 } 

 print R "XandCB <- counts\n";
 print R "elementMetadata(XandCB) <- cbind(elementMetadata(XandCB), elementMetadata(XandCB))\n";
 print R "classes <- c($params{CN0}, $params{CN1}, 1, $params{CN3}, $params{CN4})\n";
 print R "resultList <- runPanelcnMops(XandCB, 1:ncol(elementMetadata(counts)),countWindows = countWindows, I = classes, sizeFactor = \"$params{sizeFactor}\" , norm = $params{norm}, normType = \"$params{normType}\", qu = $params{qu}, quSizeFactor = $params{quSizeFactor}, priorImpact = $params{priorImpact}, minMedianRC = $params{minMedianRC}, maxControls = $params{maxControls})\n";

 print R "sampleNames <- colnames(elementMetadata(counts))\n";
 print R "finalResultsTable <- createResultTable(resultlist = resultList, XandCB = XandCB, countWindows = countWindows,sampleNames = sampleNames)\n"; 
 print R "allResults <- ldply(finalResultsTable, data.frame)\n";
 print R "outputFile <- file.path(\"$outDir\", \"panelcnmops.results.txt\")\n";

 # Build output file
 print R "colNames <- c(\"Sample\", \"Gene\", \"Chr\", \"Start\", \"End\", \"lowQual\", \"CN\")\n"; 
 print R "filteredResults <- allResults[(allResults\$CN != \"CN2\") & (allResults\$lowQual != \"lowQual\"),colNames]\n"; 
 print R "filteredResults\$CNV.type <- lapply(filteredResults\$CN, function(x) sapply(x, auxCNname))\n";
 print R "filteredResults\$CNV.type <- as.factor(unlist(filteredResults\$CNV.type))\n";
 print R "write.table(filteredResults, outputFile, sep=\'\t\', row.names=FALSE, quote = FALSE)\n"; 
 close R;
 if (!-e "$outDir/panelcnmops.results.txt"){
     `Rscript $panelcnR`;
 } 

  outputCNV($outDir);

}

##########################
sub outputCNV {
    my $inDir = shift;

    my $panelCNresults = "$inDir/panelcnmops.results.txt";
    open (IN, "<", $panelCNresults) || die " ERROR: Unable to open $panelCNresults\n";

    my %sampleCalls = ();
    while (my $line=<IN>){

        chomp $line;
        next if $line =~/^Sample/;

        my @tmp = split("\t", $line);

        my $sample = $tmp[0];
        $sample =~s/.bam//;

        my $chr    = $tmp[2];
        $chr = "chr" . $chr if $chr !~/^chr/;
        my $pos    = $tmp[3];
        my $end    = $tmp[4];
        my $svtype = $tmp[7] eq "deletion" ? "DEL" : "DUP";

        my $call = "$chr\t$pos\t$end\t$svtype";
        push @{ $sampleCalls{$sample} }, $call; 
    } 
    close IN;

    foreach my $sample ( natsort keys %sampleCalls ) {
        open (OUT, ">", "$inDir/$sample.PANELCNMOPS.bed") || die " ERROR: Unable to open $inDir/$sample.PANELCNMOPS.bed\n";
        foreach my $call (@{$sampleCalls{$sample}}){
            print OUT "$call\n";
        } 
        close OUT;  
    } 
}




return 1;