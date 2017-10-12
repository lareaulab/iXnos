# bowtie-build
# rsem-prepare-reference

GENOME_DIR=$1
GREEN_PROC_DIR=$2

if [ ! -f $GENOME_DIR/scer.transcripts.13cds10.1.ebwt ]; then
	bowtie-build $GENOME_DIR/scer.transcripts.13cds10.fa $GENOME_DIR/scer.transcripts.13cds10
fi
if [ ! -f $GENOME_DIR/scer.transcripts.13cds10.idx.fa ]; then
	rsem-prepare-reference $GENOME_DIR/scer.transcripts.13cds10.fa $GENOME_DIR/scer.transcripts.13cds10
fi
		
cat $GREEN_PROC_DIR/SRR5008134.fastq $GREEN_PROC_DIR/SRR5008135.fastq | $GREEN_PROC_DIR/trim_linker_green.pl > $GREEN_PROC_DIR/SRR5008134_SRR5008135.trimmed.fastq

bowtie -v 2 -p 36 -S --un $GREEN_PROC_DIR/green.not_rrna.fastq \
	$GENOME_DIR/ScerRRNA \
	$GREEN_PROC_DIR/SRR5008134_SRR5008135.trimmed.fastq > $GREEN_PROC_DIR/green.rrna.sam 2> $GREEN_PROC_DIR/green.rrna.bowtiestats

bowtie -v 2 -p 36 -S --un $GREEN_PROC_DIR/green.not_ncrna.fastq \
       $GENOME_DIR/rna_coding \
       $GREEN_PROC_DIR/green.not_rrna.fastq > $GREEN_PROC_DIR/green.ncrna.sam 2> $GREEN_PROC_DIR/green.ncrna.bowtiestats

bowtie -a --norc -v 2 -p 36 -S --un $GREEN_PROC_DIR/green.unmapped.fastq \
       $GENOME_DIR/scer.transcripts.13cds10 \
       $GREEN_PROC_DIR/green.not_ncrna.fastq > $GREEN_PROC_DIR/green.footprints.sam 2> $GREEN_PROC_DIR/green.footprints.bowtiestats

rsem-calculate-expression --sam $GREEN_PROC_DIR/green.footprints.sam $GENOME_DIR/scer.transcripts.13cds10 $GREEN_PROC_DIR/green 2> $GREEN_PROC_DIR/green.rsem.stderr

samtools view -h $GREEN_PROC_DIR/green.transcript.bam > $GREEN_PROC_DIR/green.transcript.sam
