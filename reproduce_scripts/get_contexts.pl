#!/usr/bin/perl

use strict;
use Bio::Seq;
use Bio::SeqIO;

my $codonlist = shift; # high_asite.csv
my $seqfile = shift; # iXnos/genome_data/scer.transcripts.13cds10.fa

my %codons;
my %testseqs;


open(CODONS, $codonlist) or die;
my $header = <CODONS>;
while(<CODONS>) {

    chomp;
    (my $gene, my $idx, my $codon) = split(',');

    $codons{$gene}{$idx} = $codon;
}

my $seq_io = Bio::SeqIO->new(
	-file   => $seqfile,
	-format => 'fasta',
    );

while (my $seq_obj = $seq_io->next_seq() ) {

    my $id = $seq_obj->display_id();
    my $seq = $seq_obj->seq();
    
    if ( $codons{$id} ) {
	$testseqs{$id} = $seq;
    }
}
    
for my $gene ( keys %codons ) {

    my $seq = $testseqs{$gene};
    
    for my $idx ( sort {$a <=> $b} keys %{$codons{$gene}} ) {

	my $start = ($idx * 3) + 13;

	# -3, -2, -1, A, 1, 2
	# $start - 9

	my @neighborhood;
	for my $i ( -3 .. 2 ) {
	    push( @neighborhood, substr( $seq, $start - $i * 3, 3 ) );
	}

	$neighborhood[3] = lc $neighborhood[3];
	
	print join " ", @neighborhood;
	print "\n";
    }
}
