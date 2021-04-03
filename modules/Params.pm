#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

package Params;

sub loadAndCheck {

    my $checkCnvkit = `which cnvkit`;
    chomp $checkCnvkit;
    if (!$checkCnvkit) {
        $checkCnvkit = "$::dirname/algorithms/cnvkit/cnvkit.py";
    }

    my $cnvkit_ztest = `which cnv_ztest.py`;
    chomp $cnvkit_ztest;
    if (!$cnvkit_ztest) {
        $cnvkit_ztest = "$::dirname/algorithms/cnvkit/scripts/cnv_ztest.py";
    }

    # my $docker = `which docker`;
    # chomp $docker;

    # if (!$docker) {
    #     print " ERROR: docker was not found on PATH\n";
    #     exit;
    # } 

    # my %dockerEnv = ();

    # # Check that DECON is available
    # my $deconImage = `$docker image ls nderoo324/decon| grep -v 'REPOSITORY' `;
    # chomp $deconImage;
 
    # if (!$deconImage){
    #     print " ERROR: decon image was not found\n";
    #     exit;
    # }
    # else {
    #     $deconImage = "nderoo324/decon:1.0.2";
    # }  

    # All executables under the same directory
    %::Callers = (
        ATLASCNV  => "$::dirname/algorithms/Atlas-CNV/atlas_cnv.pl",
        CONVADING => "$::dirname/algorithms/CoNVaDING/CoNVaDING.pl",
        CONTRA    => "$::dirname/algorithms/CONTRA.v2.0.8/contra.py",
        CNVKIT    => $checkCnvkit,
        GRAPES    => "$::dirname/algorithms/GRAPES/GRAPES",
        VISCAP    => "$::dirname/algorithms/VisCap/VisCap.R",
        BREAKMER  => "$::dirname/algorithms/BreaKmer/breakmer.py",
        DECON_READBAMS   => "ReadInBams.R",
        DECON_CALLCNV    => "makeCNVcalls.R",
        DECON_FAILEDROIS => "IdentifyFailures.R"
    );

    %::Utils = (
        VISCAP_CONFIG => "$::dirname/algorithms/VisCap/VisCap.cfg",
        ATLAS_CONFIG  => "$::dirname/algorithms/Atlas-CNV/config",
        GATK          => "$::dirname/Utils/GenomeAnalysisTK.jar",
        CNVKIT_ACCESSIBLE_REGIONS => "$::dirname/algorithms/cnvkit/data/access-5k-mappable.hg19.bed",
        CNVKIT_ZTEST  => $cnvkit_ztest,
        BREAKMER_CONFIG => "$::dirname/algorithms/BreaKmer/breakmer.cfg",
        #DECON_CONFIG => "$::dirname/algorithms/DECoN-1.0.2/decon.cfg",
        DECON_DIR => "$::dirname/algorithms/DECoN-1.0.2/Linux",
    );

    # ICR96 specific data
    %::ICR_DATA = (
        # Known CNVs
        KNOWN_VARS=> "$::dirname/ICR96_DATA/CANCER_ICR96_known_variants.txt",

        # Sample distribution in two lanes
        LANE_INFO => "$::dirname/ICR96_DATA/CANCER_ICR96_LANES.txt",

        # MLPA validation data
        MLPA      => "$::dirname/ICR96_DATA/CANCER_ICR96_MLPA_BASELINE.txt"
    );

    # Validate all executables and data is available
    foreach my $caller ( sort keys %::Callers) {
        next if $caller =~/^DECON/;  
        if (!-e $::Callers{$caller}) {
            print " ERROR: $caller is not present\n";
            exit;
        }
    }
    foreach my $util ( sort keys %::Utils) {
        if (!-e $::Utils{$util}) {
            print " ERROR: $util is not present\n";
            exit;
        }
    }
    foreach my $data ( sort keys %::ICR_DATA) {
        if (!-e $::ICR_DATA{$data}) {
            print " ERROR: $data is not present\n";
            exit;
        }
    }

    print " INFO: All executables are available\n";
}


return 1;