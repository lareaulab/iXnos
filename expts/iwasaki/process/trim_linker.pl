#!/usr/bin/perl
use strict;

#take fastq files 
#trim off any end sequence matching the linker: CTGTAGGCACCATCAAT
#(so, any terminal C, CT, CTG, CTGT, ... , up to CTGTAGGCACCATCAAT*)
#print out trimmed fastq


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
    if ( $seq =~ /^CTGTAGGCACCATCAAT/ )
    {
	next;
    }
    elsif ( ( ($trimmed) = $seq =~ /^([ACGTN]+)CTGTAGGCACCATCAAT/ ) or
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
	    ( ($trimmed) = $seq =~ /^([ACGTN]+)C$/ ) )
    {

	my $length = length($trimmed);
	if ($length > 10)
	{
	    $printstring = $header . "\n" . $trimmed . "\n" . $qheader . "\n";
	    $printstring .= substr( $qual, 0, $length );
	    $printstring .= "\n";
	}
    }
}

print $printstring;
