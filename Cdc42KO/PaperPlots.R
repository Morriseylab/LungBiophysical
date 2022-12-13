library(Seurat)
library(tidyverse)
library(scExtras)
library(corrplot)


source('~/dsdata/lungmap/bin/PaperPlotFunctions.R')

outdir <- 'Integrate'
dir.create(outdir)
plotdir <- paste0(outdir,'/PaperPlots/')
dir.create(plotdir)


scrna <- readRDS("Seurat_subset.RDS")

scrna$cellcolorsl1 <-  plyr::mapvalues(scrna$celltypel1_sample,
                                       c('AT1_Cdc42KO','AT1_wt','AT2_Cdc42KO','AT2_wt'), 
                                       c('#FB9E9E','#730707','#76B4F8','#355271')
                                       )
scrna$celltypel1_sample <- factor(scrna$celltypel1_sample, c("AT2_wt","AT2_Cdc42KO","AT1_wt","AT1_Cdc42KO"))
scrna$sample <- factor(scrna$sample, c("wt","Cdc42KO"))

cellcolors <-  scrna@meta.data %>% select(celltypel1_sample,cellcolorsl1) %>% 
  distinct() %>% 
  arrange(celltypel1_sample) %>% 
  mutate(cellcolor=as.character(cellcolorsl1)) %>% pull(cellcolor)



d1 <- DimPlot(scrna, group.by = 'celltypel1_sample',cols=cellcolors,order = T) + coord_equal() + theme_void() + ggtitle('')
d1
ggsave(paste0(plotdir,'UMAP_celltype.png'),d1,width = 6,height = 6)


scrna$sample <- factor(scrna$sample, levels=c('wt','Cdc42KO'))



s1 <- samplePlot(scrna,colorpal =c('red','blue'),nrow = 2)
s1
ggsave(paste0(plotdir,'UMAP_Sample.png'),s1,width = 6,height = 6)



# FeaturePlot -------------------------------------------------------------

DefaultAssay(scrna) <- "SCT"
genes <- c('Sftpc','Sftpb','Ager','Hopx')

f1 <- FeaturePlot(scrna,genes,order = T,min.cutoff = 'q01') * coord_equal() *  theme_void() * theme(legend.position = 'none')
f1
ggsave(paste0(plotdir,'MarkerGene_featurePlots.png'),f1,width = 6,height = 6)


# Dotplot of Markergenes --------------------------------------------------


at1genes <-  c('Ager', 'Hopx', 'Rtkn2', 'Aqp5', 'Cav1', 'Cav2', 'Wasl', 'Igfbp2', 'Itgb6', 'Sema3a', 'Sema3e', 'Col4a3')
at2genes <- c('Sftpc', 'Sftpb', 'Sftpd', 'Sftpa1', 'Lamp3', 'Abca3', 'Slc34a2', 'Lpcat1', 'Etv5', 'Il33', 'Cpm')


d1 <- DotPlot(scrna,assay='SCT',features = c(rev(at1genes),rev(at2genes)),group.by = "celltypel1_sample",cols='OrRd') + coord_flip() + xlab('') + ylab('') +
  scale_color_distiller(palette = 'OrRd',direction = 1)
d1
ggsave(paste0(plotdir,'MarkerGene_DotPlots.png'),d1,width = 8,height = 6,bg = 'white')


# Heatmap -----------------------------------------------------------------

#Just get At2


Idents(scrna) <- 'celltype_l1'
DefaultAssay(scrna) <- "RNA"
DiffExpress <- FindConservedMarkers(scrna,ident.1 = 'AT1',ident.2 = 'AT2',grouping.var = 'sample') 

HMgenes <- DiffExpress %>% rownames_to_column('gene') %>% 
  filter(wt_avg_log2FC < 0) %>% 
  #top_n(100,abs(wt_avg_log2FC))  #%>%
  head(100) #%>%
#mutate(cell=if_else(wt_avg_log2FC>0,'AT1','AT2'))
#pull(gene) 



data <- FetchData(object=scrna, vars=unique(HMgenes$gene),slot='data') %>% scale() %>% t()
annD <- scrna@meta.data %>% rownames_to_column('cellid')  %>% 
  select(cellid,sample,celltypel1_sample)


#randomize the samples 
library(ComplexHeatmap)
ha = HeatmapAnnotation(sample=annD$celltypel1_sample,
                       col=list(
                         sample=setNames(c('#FB9E9E','#730707','#76B4F8','#355271'),
                                         c('AT1_Cdc42KO','AT1_wt','AT2_Cdc42KO','AT2_wt'))
                       )
)

f1=circlize::colorRamp2(c(-2,0,2), c('skyblue1', "grey10","yellow"))

ht <-Heatmap(data,
             col=f1,
             show_row_dend = F,
             row_names_side='left',
             show_column_names = F,
             show_row_names=T,
             # right_annotation = ra,
             #added back in row names and makde the font smaller
             cluster_columns = F,
             row_names_gp = grid::gpar(fontsize = 8),
             show_column_dend = F,
             cluster_rows = F,
             column_split=annD$celltypel1_sample,
             #             use_raster=FALSE,
             #             row_split = HMgenes$cell,
             top_annotation = ha
)

ht


png(filename = paste0(plotdir,'/AT2_Heatmap.png'),width = 12,height=12,units = 'in',res=300)
draw(ht)
dev.off()







# CorrPlot ----------------------------------------------------------------

scrna <- FindVariableFeatures(scrna, assay = 'RNA')

varfeatures <- scrna@assays$integrated@var.features
mat <- AverageExpression(scrna, group.by='celltypel1_sample',features = varfeatures, assays = 'SCT', slot = 'data')

sp <- cor(mat$SCT,method='spearman')
png(filename=paste0(plotdir,'Corrplot_upper_circle.png'),width = 10,height = 8,units = 'in',res=300)
corrplot(sp, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,col.lim = c(0,1),is.corr=F,col=COL2('RdBu', 200) )
dev.off()

png(filename=paste0(plotdir,'Corrplot_upper_circle_label.png'),width = 10,height = 8,units = 'in',res=300)
corrplot(sp, method='circle',type = "upper", order = "hclust", tl.col = "black", tl.srt = 45,addCoef.col = 'white')
dev.off()




