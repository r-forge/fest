#!/usr/bin/perl
use strict;
my $usage = "Usage: thinning.pl \<ped-file-name\> \<dat-file-name\>  \<map-file-name\> <freq-file-name\> \<nNotMarker\> \<limitCentiMorgan\> \<freqThreshold\> \<nmarker\> \<suffix\>\n";
if (scalar @ARGV < 9) {
    print $usage;
    print "Thin input files to Merlin\n";
    die;
}
##my $narg = scalar @ARGV;
##print min(1, 2, 3), "\n";

my $inpedfile = $ARGV[0];
my $indatfile = $ARGV[1];
my $inmapfile = $ARGV[2];
my $infreqfile = $ARGV[3];
my $nNotMarker = $ARGV[4]; # number of columns in pedfile before 
my $limitCentiMorgan = $ARGV[5];
my $freqThreshold = $ARGV[6];
my $nmarker = $ARGV[7];
my $suffix = $ARGV[8];

##print $limitCentiMorgan, "\n";
##print $suffix, "\n";

my $outpedfile = $inpedfile . $suffix;
my $outdatfile = $indatfile . $suffix;
my $outmapfile = $inmapfile . $suffix;
my $outfreqfile = $infreqfile . $suffix;

##print $outfreqfile, "\n";

open INMAPFILE, "$inmapfile";
open INDATFILE, "$indatfile";
open INFREQFILE, "$infreqfile";

open OUTMAPFILE, ">$outmapfile";
open OUTDATFILE, ">$outdatfile";
open OUTFREQFILE, ">$outfreqfile";

my @indThinned;
my @indThinnedPed;
my $index = 0;

my $firstLine = 1;
my $firstChr = 1;
my $prevChr = -1;
my @Fprev;
while (<INMAPFILE>) {

    chomp;
    my @F = split;


    if ($firstLine) { # print header
	print OUTMAPFILE join(" ", @F), "\n";
	$firstLine = 0;
    }
    else {
##	my $tmp1 = <INFREQFILE>;
##	my @freq1 = split(' ', $tmp1);
	my @freq1 = split(' ', <INFREQFILE>);
	my $tmp2 = <INFREQFILE>;
	my @freq2 = split(' ', $tmp2);
	my $m = scalar @freq2 - 1;
	my @freqSort = sort @freq2[1..$m];
	my $maf = $freqSort[0];
	if ($F[0] != $prevChr) {
	    $firstChr = 1;
	    $prevChr = $F[0];
	}
	
	if (($firstChr || ($F[2] >= $Fprev[2] + $limitCentiMorgan)) && ($maf > $freqThreshold)) {
	    push @indThinned, $index;
	    push @indThinnedPed, (2*$index, 2*$index + 1);
	    print OUTMAPFILE join(" ", @F), "\n";
	    @Fprev = @F[0..2];
	    $firstChr = 0;
	}
	
	$index = $index + 1;
    }
}

close INFREQFILE;
open INFREQFILE, "$infreqfile";

##print "DEBUG", join(" ", @indThinned), "\n";
##print "DEBUG", join(" ", @indThinnedPed), "\n";

$index = 0;
my $indT = 0;
my $indF = 0; # 0 if M, 1 if F
while (<INFREQFILE>) {
    chomp;
    my @F1 = split;
    
    if ($index == $indThinned[$indT]) {
	print OUTFREQFILE join(" ", @F1), "\n";
	if ($indF == 1) {
	    $indT = $indT + 1;
	}
    }
    if ($indF == 1) {
	$index = $index + 1;
	$indF = 0;
    }
    else {
	$indF = 1;
    }
}
close INFREQFILE;
close OUTFREQFILE;


$index = 0;
$indT = 0;
while (<INDATFILE>) {
    chomp;
    my @F = split;

    if ($index == $indThinned[$indT]) {
	print OUTDATFILE join(" ", @F), "\n";
	$indT = $indT + 1;
    }
    $index = $index + 1;
}
close INDATFILE;
close OUTDATFILE;


open INPEDFILE, "$inpedfile";
open OUTPEDFILE, ">$outpedfile";
while (<INPEDFILE>) {

    chomp;
    my @F = split;
    my $n = scalar @F;
    my @markers = @F[$nNotMarker..$n];

##    print join(" ", @markers[@indThinnedPed]), "\n";

    print OUTPEDFILE join(" ", @F[0..4]), "  ", join(" ", @markers[@indThinnedPed]), "\n";
}

close INPEDFILE;
close OUTPEDFILE;
