#!/usr/bin/perl
use strict;
my $usage = "Usage: merge.pl \<prefix-file-name\> \<nrows\>  \<nmarkers\> \<nsim-each-pedigree\> \[alternative-pedfiles-in\] \[merged-alternative-pedfiles-out\]\n";

if (scalar @ARGV < 6) { # should at least have one altenative pedfile + merged alternativ
    print $usage;
    print "For each hypothesized pedigree relation: Merge input pedigree files, put typed individuals at\n";
    print "correct place in pedigree, and change family id such that each simulated pedigree file\n";
    print "has its own unique family id.\n";
    die;
}

# Arguments
my $prefix = $ARGV[0]; # prefix for input pedigree files
my $nrows = $ARGV[1];
my $nMarkers = $ARGV[2];
my $nSimInEachPedigree = $ARGV[3];

##my $outfilePrefix1 = $ARGV[4]; # prefix for merged true pedigree file
my $nAlt = ((scalar @ARGV) - 4)/2; # number of alternative pedigrees
my @inPrefixAlt = @ARGV[4..(4+$nAlt-1)]; # prefix for input alternative pedigree files
my @outPrefixAlt = @ARGV[(4+$nAlt)..$#ARGV]; # prefix for merged alternative pedigree files

##my $outfileTrue;
my @infileAlt;
my @outfileAlt;

##$outfileTrue = $outfilePrefix1 . ".ped";
for my $i (0..($nAlt-1)) {
    $infileAlt[$i] = $inPrefixAlt[$i] . ".ped";
    $outfileAlt[$i] = $outPrefixAlt[$i] . ".ped";
}

#print "$prefix $nrows $nSimInEachPedigree $nAlt $outfilePrefix1\n@infileAlt\n@outfileAlt\n";

## Merge all simulated pedigree files in to a single temporary file
my $tmpfile = "merged.ped.tmp";
##system("cat $prefix*.ped > $tmpfile");
my @pedfiles = <$prefix*.ped>;
open MERGEDPEDFILE, ">$tmpfile";
foreach my $pedfile (@pedfiles) {
    open INPEDFILE, "$pedfile";
    while (<INPEDFILE>) {
	print MERGEDPEDFILE $_;
    }
    close INPEDFILE;
}
close MERGEDPEDFILE;

open(INFILE, "$tmpfile") || die "Could not open $tmpfile";
##open OUTFILETRUE, ">$outfileTrue";

## Open the alternative pedigree output files
my @outAlt;
for my $i (0..($nAlt-1)) {
    open $outAlt[$i], ">$outfileAlt[$i]";
}

## Null markers
my $nMarkers2 = 2*$nMarkers;
my @markers0 = (0);
for my $i (0..($nMarkers2-1)) {
    push @markers0, (0);
}

## Read alternative input pedigrees
## 6 first columns stored in @Falt. Number of subjects stored in @nAlt.
## Store information about which subjects are typed (@typedRow)
my @Falt;
my @nAltNotTyped;
my @nAlt;
my @typedRow;
for my $i (0..($nAlt-1)) {
    open INFILEALT, "$infileAlt[$i]";
    my @Fmatrix;
    my @typed;
    my $indtyped = 0;
    my $n=0;
    while (<INFILEALT>) {
	chomp;
	my @F = split;
	
	push @Fmatrix, [ @F[0..5] ];
	if ($F[5]) {
	    $typed[$indtyped] = $n;
	    $indtyped = $indtyped + 1;
	}
	$n = $n+1;
    }
    push @Falt, [@Fmatrix];
    push @typedRow, [@typed];
    $nAltNotTyped[$i] = $n-2;
    $nAlt[$i] = $n;
    close INFILEALT;
}


## Read merged input pedigree file, insert markers of typed individuals
## into the alternative pedigrees and write this to
## the alternative output pedigrees
my @indtyped;
my $famid = 1;
my $i = 1;
my $isplit = 1;
my $newsplit = 0;
my $first;
while (<INFILE>) {

    chomp;
    my @F = split;

    if ($F[0] ne "end") {
	if ($newsplit) {
	    $newsplit = 0;
##	    close OUTFILETRUE;
	    for my $j (0..($nAlt-1)) {
		close $outAlt[$i];
	    }
	    $isplit = $isplit + 1;
##	    $outfileTrue = $outfilePrefix1 . $isplit . ".ped";
##	    open OUTFILETRUE, ">$outfileTrue";
	    for my $j (0..($nAlt-1)) {
		$outfileAlt[$j] = $outPrefixAlt[$j] . $isplit . ".ped";
		open $outAlt[$j], ">$outfileAlt[$j]";
	    }
	}

	$F[0] = $famid;
	if ($i == 1) {
	    $first = 1;
	}
	else {
	    $first = 0;
	}
##	print OUTFILETRUE join(" ", @F), "\n";
##	@F[2..3] = (0,0);

	## First time: Write all rows that are not typed to the output pedigrees
	if ($first) {
	    for my $j (0..($nAlt-1)) {
		my $out = $outAlt[$j];
		for my $k (0..($nAlt[$j]-1)) {
		    if (!$Falt[$j]->[$k]->[5]) { # !typed
			print $out "$famid ";
			for my $l (1..4) {
			    print $out "$Falt[$j]->[$k]->[$l] ";
			}
			print $out join(" ", @markers0), "\n";
		    }
		}
		$indtyped[$j] = 0;
#		$typedRow[$j] = $nAltNotTyped[$j];
#		$typedRow[$j] = $typedRow[$j] + 1;
	    }
	}

	## If individual is typed: write data to the output pedigrees
##	my $typed = ($F[5] != "0/");
	my $notTyped = ($F[5] == "0/") || ($F[5] == "0"); ## different formats
	my $typed = !$notTyped;
	if ($typed) {
	    for my $j (0..($nAlt-1)) {
		my $out = $outAlt[$j];
		print $out "$famid ";
		for my $l (1..4) {
		    my $row = $typedRow[$j][$indtyped[$j]];
		    print $out "$Falt[$j]->[$row]->[$l] ";
		}
		print $out join(" ", @F[5..$#F]), "\n";
#		print join(" ", @F[5..$#F]), "\n";
#		$typedRow[$j] = $typedRow[$j] + 1;
		$indtyped[$j] = $indtyped[$j] + 1;
	    }
	}

	if ($i < $nrows) {
	    $i = $i + 1;
	}
	else {
	    $i = 1;
	    $famid = $famid + 1;
	    if ($nSimInEachPedigree > 0) {
		if ($famid % $nSimInEachPedigree == 1) {
		    $newsplit = 1;
		}
	    }
	}
    }
}
close INFILE;
for my $i (0..($nAlt-1)) {
    close $outAlt[$i];
}

#system("rm $tmpfile")
unlink "$tmpfile";
