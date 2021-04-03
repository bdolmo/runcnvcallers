#!/usr/bin/env perl
package BREAKMER;

use strict;
use Getopt::Long;
use YAML::Tiny;
use File::Basename;

##############
sub callSV {

    my $inDir  = shift;
    my $bed    = shift;
    my $genome = shift;
    my $outDir = shift;
    
    my @bams = glob ("$inDir/*.bam");

    foreach my $bam (@bams) {
        print "$bam\n";
        my $name = basename($bam);
        $name =~s/.bam//;
        my $bam_basename = basename($bam);

        mkdir  "$outDir/BREAKMER";
        my $sampleDir = "$outDir/BREAKMER/$name";
        mkdir  "$outDir/BREAKMER/$name";

        # Heuristics to analyze only the affected exon, Breakmer is slow
        my $known_partial =  "/home/bernat/Escritorio/SIMULATIONS_TARGET/config_partial_deletions.txt";
        if (-e $known_partial) {
            my $cmd = "grep -P \'$name\t' $known_partial | cut -f 2,3,4  | $::bedtools intersect -a $bed -b stdin -wa |  sed \'s/;/_/\'> $sampleDir/$name.bed";
            system $cmd;
        }
        exit;
        # Create BREAKMER output directory
        my $tmp_config = $::Utils{BREAKMER_CONFIG};
        $tmp_config =~s/.cfg/.tmp.cfg/;

        # Rewrite config file
        open (CFG_IN, "<", $::Utils{BREAKMER_CONFIG}) or die " ERROR: Unable to open $::Utils{BREAKMER_CONFIG}\n";
        open (CFG_OUT, ">", $tmp_config) or die " ERROR: Unable to open $tmp_config\n";
        while (my $line =<CFG_IN>) {
            chomp $line;
            my @tmp = split("=", $line);
            if ($line =~/targets_bed_file/) {
                print CFG_OUT "$tmp[0]" . "=$sampleDir/$name.bed\n";
            }
            elsif ($line =~/sample_bam_file/) {
                print CFG_OUT "$tmp[0]" . "=$inDir/$bam_basename\n";
            }
            elsif ($line =~/analysis_dir/) {
                print CFG_OUT "$tmp[0]" . "=$sampleDir\n";
            }
            elsif ($line =~/reference_fasta/) {
                print CFG_OUT "$tmp[0]" . "=$genome\n";
            } 
            else {
                print CFG_OUT "$line\n";
            }
        }
        close CFG_IN;
        close CFG_OUT;

        unlink ($::Utils{BREAKMER_CONFIG});
        rename $tmp_config, $::Utils{BREAKMER_CONFIG};

        my $cmd = "python $::Callers{BREAKMER} run -c $::Utils{BREAKMER_CONFIG}";
        system $cmd;

        # Convert BreaKmer output to Bed format
        my ($breakmerResults) = glob ("$sampleDir/output/*.out");

        my $breakmerBed = "$outDir/BREAKMER/$name.BreaKmer.bed";
        open (OUT, ">", $breakmerBed)  || die " ERROR: Unable to open $breakmerBed\n";

        if ($breakmerResults) {
            open (IN, "<", $breakmerResults) || die " ERROR: Unable to open $breakmerResults\n";
            while (my $line=<IN>) {
                chomp $line;
                next if $line =~/^genes/;
                my @tmp = split ("\t", $line);

                my ($chr, $start, $end) = split(/[:\-\s]/, $tmp[1]);
                $chr = "chr$chr" if $chr !~/chr/;
                print OUT "$chr\t$start\t$end\t$tmp[0]\n";
                #genes	target_breakpoints	mismatches	strands	total_matching	sv_type	sv_subtype	split_read_count	disc_read_count	breakpoint_coverages	contig_id	contig_seq
                #NM_006939_10_11_SOS2	chr14:50625292-50625425 (D133)	0	+	138	indel	None	56	0	0,0	NM_006939_10_11_SOS2-contig1	AAATACAGTTAATTTTTCACAAAATGTCAATATAACTGAAATACATATTTTTCTTACCGTTCAATCAGTACCATATTTGGAGGCCTTCCTCCCTTTTTCTTTTTACTTTGGCTCCTCAGATTGTGGCTTTTCCTAACA
            }
            close IN;
        }
        close OUT;
    }
}

return 1;