#!/usr/bin/perl
use strict;
my $usage = "Usage: hapmapLD_divide.pl \<ld-file-name\> \<chunk-size\> \n";
if (scalar @ARGV < 2) {
    print $usage;
    print "Divide hapmap linkage disequilibrium input file in to independent chunks\n";
    die;
}

my $ldfile = $ARGV[0];
my $chunkSize = $ARGV[1];

my $nmarker = 1;
my %marker;

open LDFILE, "$ldfile";
while (<LDFILE>) {

    chomp;
    my @F = split;

    my $snp = $F[3];
    my $snp2 = $F[4];

    if (exists($marker{$snp})) {
	##$marker{$snp} = $marker{$snp} + 1;
    }
    else {
	$marker{$snp} = $nmarker;
	$nmarker = $nmarker + 1;
    }

##    $marker{$snp} = 1;
}
##print "nmarker", $nmarker, "\n";
close LDFILE;

my $nchunk = int($nmarker / $chunkSize) + 1;
my @ldfiles_out;
for my $i (0..($nchunk-1)) {
    my $outfile = $ldfile . ".chunk" . ($i+1);
    open $ldfiles_out[$i], ">$outfile";
    my $out = $ldfiles_out[$i];
    print $out "snp1 snp2 R2\n";
}

open LDFILE, "$ldfile";
while (<LDFILE>) {

    chomp;
    my @F = split;

    my $snp = $F[3];
    my $snp2 = $F[4];

    my $ichunk1 = int($marker{$snp}/$chunkSize);
    my $ichunk2 = int($marker{$snp2}/$chunkSize);

    if ($ichunk1 >= $nchunk) {
	print "$ichunk1 $nchunk\n";
	die;
    }
    if ($ichunk1 == $ichunk2) {
	my $out = $ldfiles_out[$ichunk1];
	print $out join(" ", @F[3,4,6]), "\n";
    }
##    $marker{$snp} = 1;
}
close LDFILE;
for my $i (0..($nchunk-1)) {
    close $ldfiles_out[$i];
}
print "$nchunk $nmarker\n";
