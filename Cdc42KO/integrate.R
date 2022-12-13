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

saveRDS(scrna,'Seurat.RDS')




