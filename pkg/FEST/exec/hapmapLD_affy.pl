#!/usr/bin/perl
use strict;
my $usage = "Usage: hapmapLD_affy.pl \<ld-file-name\> \<freq-file-name\> \n";
if (scalar @ARGV < 2) {
    print $usage;
    print "Summary of linkage disequilibrium input file\n";
    print "Write Hapmap LD for Affymetrix subset to file\n";
    die;
}

my $ldfile = $ARGV[0];
my $freqfile = $ARGV[1];

my $logfile = "hapmapLD_affy.log";
open LOGFILE, ">$logfile";

##my $freqfile = "../Data/Affy500K_Allele_Frequency_Files/chr22.freq.txt";
open FREQFILE, "$freqfile";
my %freqmarker;
my $nfreqmarker = 1;
print LOGFILE "Reads frequency file for SNP information\n";
while (<FREQFILE>) {

    chomp;
    my @F = split;

    my $snp = $F[0];

    if ($snp ne "---") {
	$freqmarker{$snp} = 1;
	$nfreqmarker = $nfreqmarker + 1
    }
}
close FREQFILE;
print LOGFILE "Number of markers in Affymetrix file: $nfreqmarker\n";

##my $narg = scalar @ARGV;
##print min(1, 2, 3), "\n";
open LDFILE, "$ldfile";
my $ldfile_out = $ARGV[0] . ".affy";
open LDFILEOUT, ">$ldfile_out";

my $nmarker = 1;
my %marker;

print LOGFILE "Hapmap LD file with frequency file SNPs written to file $ldfile_out\n";
while (<LDFILE>) {

    chomp;
    my @F = split;

    my $snp = $F[3];
    my $snp2 = $F[4];

    if (exists($marker{$snp})) {
	$marker{$snp} = $marker{$snp} + 1;
    }
    else {
	$marker{$snp} = 1;
    }
    if (exists($freqmarker{$snp}) && exists($freqmarker{$snp2})) {
	print LDFILEOUT join(" ", @F), "\n";
    }
##    $marker{$snp} = 1;
}
close LDFILE;
close LDFILEOUT;

my $nmarker=1;
my $key;
my $value;
$ldfile_out = $ARGV[0] . ".snplist";
print LOGFILE "Writing list of hapmap SNPs to file $ldfile_out\n";
open LDFILEOUT, ">$ldfile_out";
while (($key, $value) = each(%marker)){
    print LDFILEOUT "$key $value\n";
    $nmarker = $nmarker + 1
}
close LDFILEOUT;
print LOGFILE "Number of hapmap markers: $nmarker\n";

$nmarker=1;
my $ldfile_out2 = $ARGV[0] . ".snplist.affy";
print LOGFILE "Writing list of hapmap SNPs that are in frequency file to $ldfile_out2\n";
my %marker_sub;
open LDFILEOUT2, ">$ldfile_out2";
while ($key = each(%freqmarker)){
    if (exists($marker{$key})) {
	print LDFILEOUT2 "$key\n";
	$nmarker = $nmarker + 1
    }
}
close LDFILEOUT2;
print LOGFILE "Number of hapmap markers in Affymetrix freqfile: $nmarker\n";
close LOGFILE;
