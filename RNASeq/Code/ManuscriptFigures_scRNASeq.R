###########################################
# Project: bulkRNA-seq validation in scRNA-Seq
# Author: Cheyret Wood
# #########################################

## linux code I want to keep handy 
#PATH="$HOME/miniconda/bin:$PATH"
# screen -S mlmseurat

rm(list = ls())
# Load Libraries
# Be sure to install all of these
library(Seurat) #Seurat_4.3.0.1
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)

# project directory 
dir <- "/shared/moorelab/maier/proj/maierLiaoMould_scRNASeq/"
out <- "/shared/moorelab/maier/proj/maier_bulkRNASeq/maier_bulkRNASeq_paperAnalysis/Dissemination/"

# Read in data
nonObject <- readRDS(paste0(dir, 'DataProcessed/MaierLiaoMould_scRNASeq_non_macrophage_recluster_hybsc_05.RDS'))
load(paste0(dir, "DataProcessed/sample_phenotypes.Rdata")) # all_pheno

all_pheno <- all_pheno %>% mutate(phenotype = case_when(phenotype == "R turned P" ~ "P", TRUE ~ phenotype))

# labels
phCols <- c("steelblue1", "lightpink1", "red") # control, NP, P
dxCols <- c("steelblue1", "darkorange3") # control, Sarc

#### assigning out the cell types and ordering phenotype ####
# rename the re-clusters as per Ivana and Shu-Yi
nonObject@meta.data$celltype <- case_when(nonObject@meta.data$seurat_clusters %in% c(1,5,7,11) ~ "CD4+ T-cells",
                                          nonObject@meta.data$seurat_clusters %in% c(4,8,13) ~ "Dendritic Cells",
                                          nonObject@meta.data$seurat_clusters %in% c(0,2)  ~ "CD8+ T-cells",
                                          nonObject@meta.data$seurat_clusters == 3  ~ "T-cells", 
                                          nonObject@meta.data$seurat_clusters == 6  ~ "B Cells", 
                                          nonObject@meta.data$seurat_clusters == 9  ~ "Natural Killer Cells",
                                          nonObject@meta.data$seurat_clusters == 12 ~ "Epithelial Cells", 
                                          nonObject@meta.data$seurat_clusters == 15 ~ "Plasmacytoid Dendritic Cells", 
                                          TRUE ~ "Exclude")

nonObject@meta.data$phenotype <- case_when(nonObject@meta.data$orig.ident %in% all_pheno$id[which(all_pheno$phenotype == "C")] ~ "Control", 
                                           nonObject@meta.data$orig.ident %in% all_pheno$id[which(all_pheno$phenotype == "P")] ~ "SarcP", 
                                           nonObject@meta.data$orig.ident %in% all_pheno$id[which(all_pheno$phenotype == "R")] ~ "SarcNP", 
                                           TRUE ~ "exclude")

nonObject@meta.data$phenotype <- factor(nonObject@meta.data$phenotype, levels = c("Control", "SarcNP", "SarcP", "exclude"), ordered = T)

# order the subjects
nonObject@meta.data$sampleID <- factor(nonObject@meta.data$sampleID, levels = c("mould_1", "mould_2", "mould_3", "mould_4", "mould_5", "mould_6", "mould_7", "mould_8",
                                                                                "maier_7", "maier_8", "maier_10", "maier_14", "mould_18", "mould_19", 
                                                                                "maier_2", 'maier_5', "maier_12", 'maier_11',  "maier_15", "maier_17", "liao_29715", "liao_42787", 
                                                                                "maier_1", 'maier_3', "maier_4", "maier_9", "maier_13","maier_18", 'liao_4276', "liao_24703"), ordered = T)
# subset just to what we are interested in
cd4 <- subset(nonObject, subset = celltype %in% c("CD4+ T-cells"))
# make sure we use the assay we want
DefaultAssay(cd4) <- "SCT"

#CHN1, IL18RAP, NCR3
genes <- c("CHN1", "IL18RAP", "NCR3")

VlnPlot(cd4, features = genes, group.by = "sampleID", pt.size = 0, log = T)
ggsave(filename = paste0(out, "selectGenes_violin_CD4Tcells.png"), width = 11, height = 3.5)
ggsave(filename = paste0(out, "selectGenes_violin_CD4Tcells.pdf"), width = 11, height = 3.5)

# calculate the pseudobulk data
sample_list <- unique(cd4@meta.data$orig.ident)
avg_gene <- matrix(rep(NA, length(sample_list)*length(genes)), nrow = length(genes))
colnames(avg_gene) <- sample_list
rownames(avg_gene) <- genes

for (sample in sample_list){# within each sample 
  print(sample)
  
  if((cd4@meta.data %>% filter(orig.ident == sample) %>% nrow())> 0){
    tmp <- subset(x = cd4, orig.ident == sample)
    avg_gene[genes, sample] <- rowMeans(tmp@assays$SCT@counts[genes,])
  }else{ # if that sample doesn't have any cells in a cluster, filled with zeros 
    x <- rep(0, length(genes))
    names(x) <- genes 
    avg_gene[genes, sample] <- x
  }
}

plotDat <- merge(all_pheno, t(avg_gene), by.x = "id", by.y = 0)
plotDat <- plotDat %>% pivot_longer(CHN1:NCR3)
plotDat$id_order <- factor(plotDat$id, ordered = T, 
                            levels = c("mould_1", "mould_2", "mould_3", "mould_4", "mould_5", "mould_6", "mould_7", "mould_8",
                                       "maier_7", "maier_8", "maier_10", "maier_14", "mould_18", "mould_19", 
                                       "maier_2", 'maier_5', "maier_12", 'maier_11', "maier_15", "maier_17", "liao_29715", "liao_42787", 
                                       "maier_1", 'maier_3', "maier_4", "maier_9", "maier_13", "maier_18",'liao_4276', "liao_24703"))
plotDat$phenotype <- case_when(plotDat$phenotype == "C" ~ "Control", 
                               plotDat$phenotype == "P" ~ "SarcP", 
                               plotDat$phenotype == "R" ~ "SarcR")
plotDat$phenotype <- factor(plotDat$phenotype, levels = c("Control", "SarcR", "SarcP"), ordered = T)


plotDat %>% ggplot(aes(x = phenotype, y = value, fill = phenotype)) + 
  geom_violin() + geom_beeswarm() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = phCols) + 
  labs(title = "", y = "Avgerage Gene Expression Across Cells", x = "", fill = "Phenotype") + 
  facet_wrap(~name, scales = "free_y")
ggsave(filename = paste0(out, "violin_CD4Tcells_phenotype.png"), width = 10, height = 4)
ggsave(filename = paste0(out, "violin_CD4Tcells_phenotype.pdf"), width = 10, height = 4)



