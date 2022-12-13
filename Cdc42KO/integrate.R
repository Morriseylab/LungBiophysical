library(Seurat)
library(tidyverse)
library(scExtras)
library(clustree)
source('~/dsdata/lungmap/bin/PaperPlotFunctions.R')

outdir <- 'Integrate'
dir.create(outdir)
plotdir <- paste0(outdir,'/Plots/')
dir.create(plotdir)

m2h <- readr::read_csv('~/dsdata/NGSshare/homologs/mouse_human.csv')
cc.genes$s.genes <- m2h %>% filter(human_name %in% cc.genes$s.genes) %>% pull(mouse_name)
cc.genes$g2m.genes <- m2h %>% filter(human_name %in% cc.genes$g2m.genes) %>% pull(mouse_name)

scrna.list <- list()
scrna.list$lig = RunQC(dir=outdir,org='mouse',name='ligation',files='Mouse_Kazu_ligation_cd45neg//STARsolo/Solo.out/Gene/filtered/',filter=T, doubletdetection = T,UpperMitoCutoff=5)
scrna.list$wt <- RunQC(dir=outdir,org='mouse',name='wt',files='../Mouse_Kazu_Ctrl_cd45neg/STARsolo/Solo.out/Gene/filtered/',filter=T, doubletdetection = T,UpperMitoCutoff=5)


scrna.list <- scrna.list  %>% map(~CellCycleScoring(.x, s.features  = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes) %>%
                                    SCTransform(vst.flavor = "v2", verbose = FALSE) %>%
                                    RunPCA(npcs = 30, verbose = FALSE)
)


features <- SelectIntegrationFeatures(object.list = scrna.list,nfeatures = 3000)
scrna.list <- PrepSCTIntegration(object.list = scrna.list , anchor.features = features)


anchors <- FindIntegrationAnchors(object.list = scrna.list , normalization.method = "SCT",
                                  anchor.features = features)
scrna <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

scrna$sample <- scrna$orig.ident



scrna <- RunPCA(scrna, verbose = FALSE)
ElbowPlot(scrna,ndims = 50)
dims <- 1:30

DefaultAssay(scrna) <- 'integrated'
scrna <- RunUMAP(scrna, reduction = "pca", dims = dims, n.neighbors = 15,verbose = FALSE)
scrna <- FindNeighbors(scrna, reduction = "pca", dims = dims)
scrna <- FindClusters(scrna, resolution = 0.3)


scrna <- PrepSCTFindMarkers(scrna)


genes <- c('Hopx','Lamp3','Lyz1','Cldn4','Scgb1a1','Krt8','Foxj1','Malat1','hybrid_score')

DefaultAssay(scrna) <- 'SCT'

f1 <- FeaturePlot(scrna,genes, order=T) * coord_equal() * theme_void() * NoLegend()
d1 <- DimPlot(scrna,label = T,cols=cpallette,label.size = 8,repel = T) + theme_void() + coord_equal()
b <- sampleBarGraph(scrna,group.by = 'seurat_clusters',col = c('red','blue')) + ylab('% Cells')
#loadded from lungmap Paperfunctions code
s1 <- samplePlot(scrna,colorpal =c('red','blue'),nrow = 2)
(d1/s1)|(f1/b)

ggsave(paste0(plotdir,'Overplot1.png'),width = 12,height = 12)

genes2 <- c('Hopx','Foxj1','Sftpc','Sftpb','Lamp3','Lyz1')

v1 <- VlnPlot(scrna, genes2,split.by = 'sample',cols =  c('red','blue'))

(d1|s1)/v1 
ggsave(paste0(plotdir,'Overplot2.png'),width = 12,height = 12)

scrna
scrna@misc$findallmarkers <- FindAllMarkers(scrna,assay = 'SCT',only.pos = T)

d2 <- DimPlot(scrna, group.by = 'Phase') + theme_void() + coord_equal()
v2 <- VlnPlot(scrna, c('Pcna','Mki67','Top2a','Mcm6'),split.by = 'sample',cols =  c('red','blue'))
(d1|d2)/v2
ggsave(paste0('~/dsdata/projects/Morrisey/Kazu/ligation/Integrate/Overviewplot_CCgenes.png'),width = 12,height = 12)

saveRDS(scrna,'~/dsdata/projects/Morrisey/Kazu/ligation/Integrate/Seurat.RDS')




