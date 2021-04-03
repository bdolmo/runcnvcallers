#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use File::Basename;
use Parallel::ForkManager;
use Sort::Key::Natural qw(natsort);
use Spreadsheet::ParseExcel;
use Getopt::Long;

our $dirname = dirname (__FILE__);
use lib (dirname (__FILE__));
use modules::Params;
use modules::SysUtils;
use modules::BREAKMER;
use modules::CONTRA;
use modules::CONVADING;
use modules::CNVKIT;
use modules::DECON;
use modules::GRAPES;
use modules::VISCAP;
use modules::PANELCNMOPS;

 our $threads = 3;
 our %Callers;
 our %Utils;
 our %ICR_DATA;
 our $programs;

 Params::loadAndCheck();

 my $inDir;
 my $outDir;
 our $bed;
 my $genome;
 my $controlDir;
 my $controlSampleCONTRA;
 our $devNull  = ">/dev/null 2>&1";

our $samtools = `which samtools`;
chomp $samtools;

if (!$samtools) {
    print " ERROR: samtools was not found on path\n";
    exit;
}

our $bedtools = `which bedtools`;
chomp $bedtools;

if (!$bedtools) {
    print " ERROR: bedtools was not found on path\n";
    exit;
}

Help () if (@ARGV < 1 or !GetOptions(
	'input|i=s'=>\$inDir,
	'output|=s'=>\$outDir,
	'genome|g=s'=>\$genome,
	'bed|b=s'=>\$bed,
	'c|controldir =s'=>\$controlDir,
    'x=s'=>\$programs,
	'contra=s'=>\$controlSampleCONTRA,
));

 my $doCNVKIT   = "";
 my $doGRAPES   = "";
 my $doCONVADING= "";
 my $doVISCAP   = "";
 my $doCONTRA   = "";
 my $doBREAKMER = "";
 my $doDECON    = "";
 my $doPANELCNMOPS = "";

 my @arrProg = ( "breakmer", "contra", "grapes", "convading", "cnvkit", "viscap", "decon", "panelcnmops" );
 my %ProgHash = map { $_ => 1 } @arrProg;

 if ($programs) {
     # Parse programs
     my @parsed = split(",", lc $programs);
     foreach my $input(@parsed) {
         if (!exists $ProgHash{$input}) {
             print " ERROR: program $input is not accepted\n"; 
             Help();
         }
         if ($input eq 'contra') {
             $doCONTRA = 1;
         }
         if ($input eq 'grapes') {
             $doGRAPES = 1;
         }
         if ($input eq 'convading') {
             $doCONVADING = 1;
         }
         if ($input eq 'viscap') {
             $doVISCAP = 1;
         }
         if ($input eq 'cnvkit') {
             $doCNVKIT = 1;
         }
         if ($input eq 'breakmer') {
             $doBREAKMER = 1;
         }
         if ($input eq 'decon') {
             $doDECON = 1;
         }
         if ($input eq 'panelcnmops'){
             $doPANELCNMOPS = 1;
         } 
     }
 }

 if (!$outDir)  {
    print " ERROR: output directory (-o) was not specified\n";
    Help();
    exit;
 }
 mkdir $outDir;
 if ($doDECON) {
    DECON::callCNV($inDir, $bed, $genome, $outDir);
 }
 if ($doCNVKIT) {
    CNVKIT::callCNV($inDir, $bed, $genome, $outDir);
 }
 if ($doGRAPES) {
    GRAPES::callCNV($inDir, $bed, $genome, $outDir);
 }
 if ($doCONVADING) {
    CONVADING::callCNV($inDir, $bed, $genome, $outDir, $controlDir);
 }
 if ($doCONTRA) {
    CONTRA::callCNV($inDir, $controlSampleCONTRA, $bed, $outDir);
 }
 if ($doVISCAP) {
    VISCAP::callCNV($inDir, $bed, $genome, $outDir);
 }
 if ($doBREAKMER) {
    BREAKMER::callSV($inDir, $bed, $genome, $outDir);
 }
 if ($doPANELCNMOPS) {
    PANELCNMOPS::callCNV($inDir, $bed, $genome, $outDir);
 } 
  #doAtlasCNV($inDir, $bed, $genome, $outDir);


##############################################
 sub Help {
     
  print "\n Usage: $0 <PARAMS>

  -i,--input        STRING     Input BAM directory
  -o,--output       STRING     Output directory
  -g,--genome       STRING     Genome in FASTA format
  -b,--bed          STRING     BED roi file
  -c,--controldir   STRING     Directory with controls
  -x                STRING     Programs separated by comma: breakmer,contra,convading,cnvkit,decon,grapes,viscap
  -contra           STRING     Control sample for CONTRA\n\n";
	exit;
 }