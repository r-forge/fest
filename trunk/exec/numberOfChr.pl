#!/usr/bin/perl
use strict;
my $usage = "Usage: numberOfChr.pl \<map-file-name\>\n";
if (scalar @ARGV != 1) {
    print $usage;
    print "Thin input files to Merlin\n";
    die;
}
##my $narg = scalar @ARGV;
##print min(1, 2, 3), "\n";

my $inmapfile = $ARGV[0];

open INMAPFILE, "$inmapfile";

my %chr;
my $firstLine = 1;
while (<INMAPFILE>) {

    chomp;
    my @F = split;

    if ($firstLine) { # header
	$firstLine = 0;
    }
    else {
	$chr{$F[0]} = 1;
    }
}

close INMAPFILE;

my $key;
my $value;;
my $nchr = 0;
while (($key, $value) = each(%chr)){
    $nchr = $nchr + 1
}

print $nchr, "\n";
