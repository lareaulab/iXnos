#!/usr/bin/perl

# /mnt/lareaulab/rtunney/Bias/RawData/SRR1049521.fastq
    
# 3' adapter sequence:
# 5' AppTCGTATGCCGTCTTCTGCTTGidT 3'
# 8 random nt at the beginning of the read

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
    if ( $seq =~ /^[ACGTN]{8}TCGTATGCCGTCTTCTGCTTG/ )
    {
	next;
    }
    elsif ( ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGCCGTCTTCTGCTTG/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGCCGTCTTCTGCTTG$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGCCGTCTTCTGCTT$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGCCGTCTTCTGCT$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGCCGTCTTCTGC$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGCCGTCTTCTG$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGCCGTCTTCT$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGCCGTCTTC$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGCCGTCTT$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGCCGTCT$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGCCGTC$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGCCGT$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGCCG$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGCC$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATGC$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTATG$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTAT$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGTA$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCGT$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TCG$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)TC$/ ) or
	    ( ($trimmed) = $seq =~ /^[ACGTN]{8}([ACGTN]+)T$/ ) )
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

