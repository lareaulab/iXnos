GENOME_DIR=$1
WEINBERG_PROC_DIR=$2

if [ ! -f $GENOME_DIR/scer.transcripts.13cds10.1.ebwt ]; then
bowtie-build $GENOME_DIR/scer.transcripts.13cds10.fa $GENOME_DIR/scer.transcripts.13cds10
fi
if [ ! -f $GENOME_DIR/scer.transcripts.13cds10.idx.fa ]; then
rsem-prepare-reference $GENOME_DIR/scer.transcripts.13cds10.fa $GENOME_DIR/scer.transcripts.13cds10
fi

cat $WEINBERG_PROC_DIR/SRR1049521.fastq | $WEINBERG_PROC_DIR/trim_linker_weinberg.pl > $WEINBERG_PROC_DIR/SRR1049521.trimmed.fastq

bowtie -v 2 -p 36 -S --un $WEINBERG_PROC_DIR/weinberg.not_rrna.fastq \
	$GENOME_DIR/ScerRRNA \
	$WEINBERG_PROC_DIR/SRR1049521.trimmed.fastq > $WEINBERG_PROC_DIR/weinberg.rrna.sam 2> $WEINBERG_PROC_DIR/weinberg.rrna.bowtiestats

bowtie -v 2 -p 36 -S --un $WEINBERG_PROC_DIR/weinberg.not_ncrna.fastq \
       $GENOME_DIR/rna_coding \
       $WEINBERG_PROC_DIR/weinberg.not_rrna.fastq > $WEINBERG_PROC_DIR/weinberg.ncrna.sam 2> $WEINBERG_PROC_DIR/weinberg.ncrna.bowtiestats

bowtie -a --norc -v 2 -p 36 -S --un $WEINBERG_PROC_DIR/weinberg.unmapped.fastq \
       $GENOME_DIR/scer.transcripts.13cds10 \
       $WEINBERG_PROC_DIR/weinberg.not_ncrna.fastq > $WEINBERG_PROC_DIR/weinberg.footprints.sam 2> $WEINBERG_PROC_DIR/weinberg.footprints.bowtiestats

rsem-calculate-expression --sam $WEINBERG_PROC_DIR/weinberg.footprints.sam $GENOME_DIR/scer.transcripts.13cds10 $WEINBERG_PROC_DIR/weinberg 2> $WEINBERG_PROC_DIR/weinberg.rsem.stderr

samtools view -h $WEINBERG_PROC_DIR/weinberg.transcript.bam > $WEINBERG_PROC_DIR/weinberg.transcript.sam
