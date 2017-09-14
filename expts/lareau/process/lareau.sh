# bowtie-build
# rsem-prepare-reference

GENOME_DIR=$1
LAREAU_PROC_DIR=$2

if [ ! -f $GENOME_DIR/scer.transcripts.13cds10.1.ebwt ]; then
	bowtie-build $GENOME_DIR/scer.transcripts.13cds10.fa $GENOME_DIR/scer.transcripts.13cds10
fi
if [ ! -f $GENOME_DIR/scer.transcripts.13cds10.idx.fa ]; then
	rsem-prepare-reference $GENOME_DIR/scer.transcripts.13cds10.fa $GENOME_DIR/scer.transcripts.13cds10
fi
		
cat $LAREAU_PROC_DIR/LLMG004_S31_L007_R1_001.fastq $LAREAU_PROC_DIR/LLMG005_S1_L001_R1_001.fastq | $LAREAU_PROC_DIR/trim_linker_lareau.pl > $LAREAU_PROC_DIR/LLMG004_005.trimmed.fastq

~/bowtie-1.2.1.1/bowtie -v 2 -p 36 -S --un $LAREAU_PROC_DIR/lareau.not_rrna.fastq \
	$GENOME_DIR/ScerRRNA \
	$LAREAU_PROC_DIR/LLMG004_005.trimmed.fastq > $LAREAU_PROC_DIR/lareau.rrna.sam 2> $LAREAU_PROC_DIR/lareau.rrna.bowtiestats

~/bowtie-1.2.1.1/bowtie -v 2 -p 36 -S --un $LAREAU_PROC_DIR/lareau.not_ncrna.fastq \
       $GENOME_DIR/rna_coding \
       $LAREAU_PROC_DIR/lareau.not_rrna.fastq > $LAREAU_PROC_DIR/lareau.ncrna.sam 2> $LAREAU_PROC_DIR/lareau.ncrna.bowtiestats

~/bowtie-1.2.1.1/bowtie -a --norc -v 2 -p 36 -S --un $LAREAU_PROC_DIR/lareau.unmapped.fastq \
       $GENOME_DIR/scer.transcripts.13cds10 \
       $LAREAU_PROC_DIR/lareau.not_ncrna.fastq > $LAREAU_PROC_DIR/lareau.footprints.sam 2> $LAREAU_PROC_DIR/lareau.footprints.bowtiestats

rsem-calculate-expression --sam $LAREAU_PROC_DIR/lareau.footprints.sam $GENOME_DIR/scer.transcripts.13cds10 $LAREAU_PROC_DIR/lareau 2> $LAREAU_PROC_DIR/lareau.rsem.stderr

samtools view -h $LAREAU_PROC_DIR/lareau.transcript.bam > $LAREAU_PROC_DIR/lareau.transcript.sam
