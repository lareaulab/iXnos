#bowtie-build $GENOME_DIR/gencode.v22.transcript.13cds10.fa $GENOME_DIR/human.transcripts.13cds10
#rsem-prepare-reference $GENOME_DIR/gencode.v22.transcript.13cds10.fa $GENOME_DIR/human.transcripts.13cds10

GENOME_DIR=$1
IWASAKI_PROC_DIR=$2

echo $GENOME_DIR

if [ ! -f $GENOME_DIR/human.transcripts.13cds10.1.ebwt ]; then
	bowtie-build $GENOME_DIR/gencode.v22.transcript.13cds10.fa $GENOME_DIR/human.transcripts.13cds10
fi
if [ ! -f $GENOME_DIR/human.transcripts.13cds10.idx.fa ]; then
	echo "I'm making a human index!!!"
	rsem-prepare-reference $GENOME_DIR/gencode.v22.transcript.13cds10.fa $GENOME_DIR/human.transcripts.13cds10
fi

cat $IWASAKI_PROC_DIR/SRR2075925.fastq $IWASAKI_PROC_DIR/SRR2075925.fastq | $IWASAKI_PROC_DIR/trim_linker.pl > $IWASAKI_PROC_DIR/SRR2075925_SRR2075926.trimmed.fastq

~/bowtie-1.2.1.1/bowtie -v 2 -p 36 -S --un $IWASAKI_PROC_DIR/iwasaki.not_rrna_trna.fastq \
	$GENOME_DIR/human_rrna_trna \
	$IWASAKI_PROC_DIR/SRR2075925_SRR2075926.trimmed.fastq > $IWASAKI_PROC_DIR/iwasaki.rrna_trna.sam 2> $IWASAKI_PROC_DIR/iwasaki.rrna_trna.bowtiestats

~/bowtie-1.2.1.1/bowtie -a --norc -v 2 -p 36 -S --un $IWASAKI_PROC_DIR/iwasaki.unmapped.fastq \
	$GENOME_DIR/human.transcripts.13cds10 \
	$IWASAKI_PROC_DIR/iwasaki.not_rrna_trna.fastq > $IWASAKI_PROC_DIR/iwasaki.footprints.sam 2> $IWASAKI_PROC_DIR/iwasaki.footprints.bowtiestats

rsem-calculate-expression --sam $IWASAKI_PROC_DIR/iwasaki.footprints.sam $GENOME_DIR/human.transcripts.13cds10 $IWASAKI_PROC_DIR/iwasaki 2> $IWASAKI_PROC_DIR/iwasaki.rsem.stderr

samtools view -h $IWASAKI_PROC_DIR/iwasaki.transcript.bam > $IWASAKI_PROC_DIR/iwasaki.transcript.sam

