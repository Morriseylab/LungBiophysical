~/bin/STAR --soloType CB_UMI_Simple \
--soloUMIlen 12 \
--soloFeatures Gene Velocyto \
--soloCellFilter EmptyDrops_CR \
--soloUMIfiltering MultiGeneUMI \
--soloCBmatchWLtype 1MM_multi_pseudocounts \
--soloCBwhitelist ~/NGSshare/3M-february-2018.txt \
--readFilesIn ~/NGS/30-624603005/00_fastq/EEM-scRNA-113_S1_L001_R2_001.fastq.gz ~/NGS/30-624603005/00_fastq/EEM-scRNA-113_S1_L001_R1_001.fastq.gz \
--genomeDir ~/NGSshare/mm39_tdTomato_STAR \
--runThreadN 64 \
--readFilesCommand zcat \
--soloBarcodeReadLength 0 \
--outSAMtype BAM Unsorted \
--outFileNamePrefix STARsolo/

gzip STARsolo/Solo.out/Gene/filtered/*

