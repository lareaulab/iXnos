#!/usr/bin/perl
use strict;

# May 2017: fixed quality string bug

# Trim off end sequence matching the barcode and linker, AND random nt
# (5 at end - this experiment used the RT primer with no random nt)

# this version does not allow de-duplication

# full linker sequence:
# NI-814 GATCA 5’-/5Phos/NNNNNGATCAAGATCGGAAGAGCACACGTCTGAA/3ddC/
# NI-815 GATCA 5’-/5Phos/NNNNNGCATAAGATCGGAAGAGCACACGTCTGAA/3ddC/
# NNNNNiiiiiAGATCGGAAGAGCACACGTCTGAAC

# Untreated: NI-814 5' GATCA 3'
# Cycloheximide treated: NI-815 5' GCATA 3'
# Library reverse indexing primer: NI-799 
# RT primer: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGnnnnGTGACTGGAGTTCAGACGTGTGCTC

my $barcode = shift;

my $printstring = "";

while(<>)
{
    my $header;
    my $seq;
    my $qheader;
    my $qual;
    
    if (/^@/)
    { 
	print $printstring;
	$printstring = "";

	$header = $_; chomp $header;
	$seq = <>; chomp $seq;
	$qheader = <>; chomp $qheader;
	$qual = <>; chomp $qual;
    }

    my $trimmed;

    if ( $seq =~ /^[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTCTGAAC/ )
    {
	next;
    }

    # no random barcodes at beginning
    elsif ( ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTCTGAAC/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTCTGAAC$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTCTGAA$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTCTGA$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTCTG$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTCT$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTC$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGT$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACG$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACAC$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACA$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCAC$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCA$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGC$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAG$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGA$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAG$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAA$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGA$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGG$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCG$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATC$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGAT$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGA$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AG$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}A$/ ) or
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}$/ ))
    {
	
	my $length = length($trimmed);

	if ($length > 10)
	{
	    $printstring = $header . "\n" . $trimmed . "\n" . $qheader . "\n";
	    # no random barcodes at beginning
	    $printstring .= substr( $qual, 0, $length );
	    $printstring .= "\n";
	}
    }
    print $printstring;
    $printstring = "";
}


