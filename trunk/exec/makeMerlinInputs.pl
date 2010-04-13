#!/usr/bin/perl
use strict;
my $usage = "Usage: makeMerlinInputs.pl \<prefix-file-name\> \<chr-length\> \<n.each\> \<n.chr\> \<n.alleles\> [\<frequencies\>]\n";

if (scalar @ARGV < 5) {
    print $usage;
    print "Make input files to Merlin\n";
    die;
}

my $narg = scalar @ARGV;

my $prefix = $ARGV[0];
my $chrLength = $ARGV[1];
my $neach = $ARGV[2];
my $nchr = $ARGV[3];
my $nall = $ARGV[4];

my $i;
my @freq;
if ($narg > 5) {
    for $i (1..$nall) {
	push @freq, ($ARGV[$i+4]);
    }
}
else
{
    for $i (1..$nall) {
	push @freq, (1.0/$nall);
    }
}
  


my $freqfile = $prefix . ".freq";
my $datfile = $prefix . ".dat";
my $mapfile = $prefix . ".map";

#print "neach = $neach, nchr = $nchr, nalleles = $nall\n";

open FREQ, ">$freqfile";
open MAP, ">$mapfile";
open DAT, ">$datfile";

## Write ped file
my $nMarkers = $neach * $nchr;
my $nMarkers2 = 2* $nMarkers;
my @markers0 = (0);
my @markers1 = (1);
for $i (2..$nMarkers2) {
    push @markers0, (0);
    push @markers1, (1);
}

my $j;
#my $freq = 1.0/$nall;
for $i (1..$nMarkers) {
    print FREQ "M ", "M".$i, "\n";
    for $j (0..$nall-1) {
	print FREQ "F ", $freq[$j], "\n";
    }
}

close FREQ;

#my @lim = (0,190);
my @lim = (0,$chrLength);
my $deltapos = 0;
if ($neach > 1) {
    $deltapos = ($lim[1] - $lim[0])/($neach-1);
}
my $k = 1;
print MAP "CHR MARKER POS\n";
for $i (1..$nchr) {
    my $pos = 0.0;
    for $j (1..$neach) {
	print MAP $i, " ", "M".$k, " ", $pos, "\n";
	$k = $k + 1;
	$pos = $pos + $deltapos;
    }
}

close MAP;

for $i (1..$nMarkers) {
    print DAT "M ", "M".$i, "\n";
}
print DAT "A locus1\n";

close DAT;
