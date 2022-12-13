library(Seurat)
library(tidyverse)
library(scExtras)
library(clustree)
source('~/dsdata/lungmap/bin/PaperPlotFunctions.R')

outdir <- 'Integrate'
dir.create(outdir)
plotdir <- paste0(outdir,'/Plots/')
dir.create(plotdir)


all <- readRDS(paste0(outdir,'/Seurat.RDS'))
d1 <- DimPlot(all,label = T,cols=cpallette,label.size = 8,repel = T) + theme_void() + coord_equal()
d1

#Subset
DefaultAssay(all) <- 'RNA'


#scrna <- subset(all, subset= (Lamp3 > 2 | Hopx>2) & Foxj1 < 1 & Mgp < 1)
scrna <- subset(all, idents = c(0,2,3,4))%>% subset(subset = nCount_RNA>1000)


#Re UMAP and cluster

DefaultAssay(scrna) <- 'integrated'
scrna <- RunPCA(scrna)
ElbowPlot(scrna,ndims = 50)
dims <- 1:20
scrna <- RunUMAP(scrna, reduction = "pca", dims = dims, n.neighbors = 10,min.dist = .4,verbose = FALSE)
scrna <- FindNeighbors(scrna, reduction = "pca", dims = dims)
scrna <- FindClusters(scrna, algorithm = 2,resolution = 0.3)


genes <- c('Hopx','Lamp3','Lyz1','Cldn4','Scgb1a1','Krt8','Foxj1','Malat1','hybrid_score')

DefaultAssay(scrna) <- 'SCT'


f1 <- FeaturePlot(scrna,genes, order=T) * coord_equal() * theme_void() * NoLegend()
d1 <- DimPlot(scrna,label = T,cols=cpallette,label.size = 8,repel = T) + theme_void() + coord_equal()
b <- sampleBarGraph(scrna,group.by = 'seurat_clusters',col = c('red','blue')) + ylab('% Cells')
#loadded from lungmap Paperfunctions code
s1 <- samplePlot(scrna,colorpal =c('red','blue'),nrow = 2)
(d1/s1)|(f1/b)
ggsave(paste0(plotdir,'Subset_overviewplot2.png'),width = 12,height = 12)



scrna$celltype_l2 <- plyr::mapvalues(scrna$seurat_clusters,seq(0,4),c('AT1a','AT2a','AT1b','AT2b', 'AT2c'))
scrna$celltype_l1 <- plyr::mapvalues(scrna$seurat_clusters,seq(0,4),c('AT1','AT2','AT1','AT2','AT2'))
scrna$celltypel1_sample <- paste0(scrna$celltype_l1, '_',scrna$sample)
scrna$celltypel2_sample <- paste0(scrna$celltype_l2, '_',scrna$sample)
  
Idents(scrna) <- 'celltypel1_sample'


saveRDS(scrna,paste0(outdir,'/Seurat_subset.RDS'))

scrna$celltypel1_sample

# Lvl 1 test --------------------------------------------------------------


tests <- list(
  'AT2ko_WT' = list(A = 'AT2_Cdc42KO',B='AT2_wt'),
  'AT1ko_WT' = list(A = 'AT1_Cdc42KO',B='AT1_wt')
)


Idents(scrna) <- 'celltypel1_sample'
dge <- map(tests,function(c){
  FindMarkers(scrna, ident.1 = c$A, ident.2 = c$B, verbose = FALSE,recorrect_umi = FALSE) %>% rownames_to_column('gene') %>% 
    mutate(diffpct=pct.1-pct.2)
}
)

writexl::write_xlsx(dge,paste0(outdir,'/DiffGeneExpress_subset.xlsx'))





#Look at the diff between the Lyz1+ Lyz1 cells. 
Idents(scrna) <- 'celltype'

at2lyz1 <- FindMarkers(scrna, ident.1 = 'AT2b', ident.2 = 'AT2a', verbose = FALSE,recorrect_umi = FALSE) %>% rownames_to_column('gene') %>% 
  mutate(diffpct=pct.1-pct.2)
writexl::write_xlsx(at2lyz1,paste0(outdir,'/AT2Lyz1_diffexp.xlsx'))

