#!/usr/bin/perl
use strict;
my $usage = "Usage: makepedigree.pl \<inped-file-name\> \<prefix-file-name\> \<n.marker\> \<typed\>\n";

if (scalar @ARGV < 4) {
    print $usage;
    print "Make pedigree input file to Merlin\n";
    die;
}

my $narg = scalar @ARGV;

my $inped = $ARGV[0];
my $prefix = $ARGV[1];
my $nMarkers = $ARGV[2];
my %typed = ();
#{@ARGV[3..$narg]};
foreach (@ARGV[3..$narg]) {
  $typed{ $_ } = 1;
}


my $pedfile = $prefix . ".ped";

open(INPED, "$inped") || die "Could not open $inped";
open PED, ">$pedfile";

## Write ped file
my $nMarkers2 = 2* $nMarkers;
my @markers0 = (0);
my @markers1 = (1);
my $i;
for $i (2..$nMarkers2) {
    push @markers0, (0);
    push @markers1, (1);
}

while (<INPED>) {
    chomp;
    my @F = split;
    my @F2;
    if (exists  $typed{$F[1]}) {
	@F2 = (@F[0..4], @markers1, $F[5]);
    }
    else {
	@F2 = (@F[0..4], @markers0, $F[5]);
    }

    print PED join(" ", @F2), "\n";
}
close PED;
