~/bin/STAR --soloType CB_UMI_Simple \
--soloUMIlen 12 \
--soloFeatures Gene Velocyto \
--soloCellFilter EmptyDrops_CR \
--soloUMIfiltering MultiGeneUMI \
--soloCBmatchWLtype 1MM_multi_pseudocounts \
--soloCBwhitelist ~/NGSshare/3M-february-2018.txt \
--readFilesIn R2.fastq.gz R1.fastq.gz \
--genomeDir ~/NGSshare/mm39_STAR \
--runThreadN 64 \
--readFilesCommand zcat \
--soloBarcodeReadLength 0 \
--outSAMtype BAM Unsorted \
--outFileNamePrefix STARsolo/

gzip STARsolo/Solo.out/Gene/filtered/*

