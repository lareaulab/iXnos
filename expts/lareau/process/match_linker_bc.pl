#!/usr/bin/perl
use strict;

# Match sequences ending in 5 random nt + 5-nt barcode + linker (any prefix of linker)
# (this experiment used the RT primer with no random nt)

# this version does not allow de-duplication

# full linker sequence:
# NNNNNiiiiiAGATCGGAAGAGCACACGTCTGAAC

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

    # no random barcodes at beginning
    if ( ( $seq =~ /^([ACGTN]*)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTCTGAAC/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTCTGAAC$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTCTGAA$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTCTGA$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTCTG$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTCT$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGTC$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACGT$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACACG$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACAC$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCACA$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCAC$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGCA$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAGC$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGAG$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAGA$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAAG$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGAA$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGGA$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCGG$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATCG$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGATC$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGAT$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AGA$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}AG$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}A$/ ) or
	 ( $seq =~ /^([ACGTN]+)[ACGTN]{5}${barcode}$/ ))
    {
	$printstring = $header . "\n" . $seq . "\n" . $qheader . "\n";
	# no random barcodes at beginning
	$printstring .= $qual;
	$printstring .= "\n";
    }
    print $printstring;
    $printstring = "";
}


