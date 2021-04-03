#!/usr/bin/env perl
package DECON;

use strict;
use Getopt::Long;
use File::Basename;
use Sort::Key::Natural qw(natsort);
use YAML::Tiny;
use List::MoreUtils qw(uniq);

#####################################
sub callCNV {

 my $inDir = shift;
 my $bed   = shift;
 my $genome= shift;
 my $outDir= shift;

 $outDir = "$outDir/DECON";
 mkdir $outDir;

 # Open the config
 my $yaml = YAML::Tiny->read( "$::dirname/modules/decon.yaml" );
 
 # Get a reference to the first document
 my $deconConfig = $yaml->[0];

 # Get a hash of params
 my %params = %{ $deconConfig };

 chdir $::Utils{DECON_DIR}; 


 my $outputRdata = "$outDir/OUTPUT.RData";
 my $failedRois  = "$outDir/FAILED";
 my $failedRoisWhole = $failedRois . "_Failures.txt";
 my $cnvCalls    = "$outDir/CNV";
 my $cnvResults  = "$outDir/CNV_all.txt";

 # Get count data
 my $cmd = "Rscript $::Callers{DECON_READBAMS} --bams $inDir --bed $bed --fasta $genome --out $outDir/OUTPUT";
 print "$cmd\n";
 system $cmd if !-e $outputRdata;

 # Call part 2
 $cmd = "Rscript $::Callers{DECON_FAILEDROIS} --Rdata $outputRdata --mincorr $params{mincorr} --mincov $params{mincov} --out $failedRois";
 print "$cmd\n";
 system $cmd if !-e $failedRoisWhole;

 # Call part 3
 $cmd = "Rscript $::Callers{DECON_CALLCNV} --Rdata $outputRdata --transProb $params{transProb} --out $cnvCalls";
 print "$cmd\n";
 system $cmd if !-e $cnvResults;

 # Write sample level calls to BED
 outputCNV($outDir);

}
##########################
sub outputCNV {
    my $inDir = shift;

    my $deconCNV = "$inDir/CNV_all.txt";
    open (IN, "<", $deconCNV) || die " ERROR: Unable to open $deconCNV\n";

    my %sampleCalls = ();
    while (my $line=<IN>){
        chomp $line;
        my @tmp = split("\t", $line);
        next if $line =~/^CNV/;

        my $sample = $tmp[1]; 
        my $chr    = $tmp[10];
        $chr = "chr" . $chr if $chr !~/^chr/;
        my $pos    = $tmp[8];
        my $end    = $tmp[9];
        my $svtype = $tmp[6] eq "deletion" ? "DEL" : "DUP";

        my $call = "$chr\t$pos\t$end\t$svtype";
        push @{ $sampleCalls{$sample} }, $call; 
    } 
    close IN;

    foreach my $sample ( natsort keys %sampleCalls ) {
        open (OUT, ">", "$inDir/$sample.DECON.bed") || die " ERROR: Unable to open $inDir/$sample.DECON.bed\n";
        foreach my $call (uniq@{$sampleCalls{$sample}}){
            print OUT "$call\n";
        } 
        close OUT;  
    } 
}

return 1; 