#!/usr/bin/perl
use strict;

# Trim off end sequence matching the linker from Rachel Green's experiments:
# CTGTAGGCACCATCAAT... (it can be longer than this)

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

    # no random barcodes at beginning
    if ( ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGTAGGCACCATCAAT/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGTAGGCACCATCAA$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGTAGGCACCATCA$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGTAGGCACCATC$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGTAGGCACCAT$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGTAGGCACCA$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGTAGGCACC$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGTAGGCAC$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGTAGGCA$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGTAGGC$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGTAGG$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGTAG$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGTA$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGT$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CTG$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)CT$/ ) or
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)C$/ ) or ## ...eats up anything ending in a C... but that's ok
	 ( ($trimmed) = $seq =~ /^([ACGTN]+)$/ )  ## ...no evidence of linker..
	)
    {
    
	my $length = length($trimmed);

	if ($length > 10)
	{
	    $printstring = $header . "\n" . $trimmed . "\n" . $qheader . "\n";
	    $printstring .= substr( $qual, 0, $length ); ## fix this bug in other versions of the script!!
	    $printstring .= "\n";
	}
    }
    
    print $printstring;
    $printstring = "";
}
