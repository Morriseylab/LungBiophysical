


# trimmomatic
trimmomatic SE -threads 3 -trimlog AT1_1_Input.txt AT1_1_Input.fq AT1_1_Input.trimmom.fastq ILLUMINACLIP:illumina_adapters.fa:2:30:10

# mapping using mm10
bwa aln -q 5 -l 32 -k 2 /Users/*/Desktop/mm10/Sequence/BWAIndex/genome.fa AT1_1_Input.trimmom.fastq > bwa/AT1Input1.sai
bwa samse /Users/*/Desktop/mm10/Sequence/BWAIndex/genome.fa bwa/AT1Input1.sai AT1_1_Input.trimmom.fastq > bwa/AT1Input1.sam

# Make BAM files
samtools view -bhS -F 1804 -q 30 bwa/AT1LB1-1.sam | samtools sort -T bwa/AT1LB1-1.sam - > bwa/AT1-LB1-1.bam
samtools index bwa/AT1-Input1.bam

# Remove Duplicates
picard MarkDuplicates I=AT1-LB1-1.bam O=AT1-LB1-1.dupl.bam M=AT1-LB1-1.bam.log.txt REMOVE_DUPLICATES=true
samtools index AT1-LB1-1.dupl.bam

# make bigwig files from BAM
bamCoverage -b AT2-input.bam -o AT2-input.bw -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 -e 200 -bl mm10-blacklist.v2.bed

