###########################################
# Project: bulkRNA-seq sensitivity 
# Author: Cheyret Wood
# #########################################

rm(list = ls())

library(tidyverse)
library(openxlsx)
library(edgeR)
library(DESeq2)

`%notin%` <- Negate(`%in%`)

options(future.globals.maxSize = 100*1000 * 1024^2)

source("processingData.R")
source("XX_formattingFunctions.R")
source("XX_enrichmentFunctions.R")

dataDir <- "P:/BRANCHES/NJHMoore_CBD_Sarcoidosis/maier_bulkRNASeq_paperAnalysis/DataRaw/"
outDir  <- "P:/BRANCHES/NJHMoore_CBD_Sarcoidosis/maier_bulkRNASeq_paperAnalysis/"

# get the prepped data
outDat <- getDataStarted(dataDir, 6, outliers = T)
counts.filtered <- outDat[[1]]
all_md <- outDat[[2]]
fdata <- outDat[[3]]
genes <- outDat[[4]]

# reorder the phenotype so we can get the comparsions we want easily 
all_md$Phenotype <- case_when(all_md$Phenotype == "control" ~ "cntl",
                              all_md$Phenotype == "NP" ~ "aNP",
                              all_md$Phenotype == "P" ~ "zP")

# functions for running models
runDDS_DX <- function(samples, sex = T, age = F){
  dsgn <- as.formula(ifelse(sex, "~ Birth.sex + Dx", "~ Dx"))
  dsgn <- as.formula(ifelse(age, "~ Age.at.sampling + Dx", "~ Dx"))
  
  # Create DESeq2 Dataset and Model
  dds <- DESeqDataSetFromMatrix(countData = counts.filtered[,samples], colData = all_md[samples,],
                                design = dsgn, rowData = fdata)
  dds <- DESeq(dds)
  ddr <- results(dds, contrast = c("Dx", "Sarc", "Control"))
  # Order results by p-value
  ddr <- ddr[order(ddr$padj), ]
  # Add in the gene symbols
  ddr$gene <- genes[rownames(ddr)]
  
  return(as.data.frame(ddr))
}

runDDS_pNp <- function(samples, sex = T){
  dsgn <- as.formula(ifelse(sex, "~ Birth.sex + Phenotype", "~ Phenotype"))
  
  # Create DESeq2 Dataset and Model
  dds <- DESeqDataSetFromMatrix(countData = counts.filtered[,samples], colData = all_md[samples,],
                                design = dsgn, rowData = fdata)
  dds <- DESeq(dds)
  # Code contrasts to compare groups
  ddr <- results(dds)
  # Order results by p-value
  ddr <- ddr[order(ddr$padj), ]
  # Add in the gene symbols
  ddr$gene <- genes[rownames(ddr)]
  
  return(as.data.frame(ddr))
}

plotLog <- function(dat, y = "_males"){
  dat <- dat %>% mutate(threshold = case_when(padj_all < 0.1  & get(paste0("padj", y)) < 0.1  ~ "FDR < 0.1 Both",
                                              padj_all < 0.1  & get(paste0("padj", y)) >= 0.1 ~ "FDR < 0.1 Full Group Only",
                                              TRUE ~ "Out")) %>% 
    filter(threshold != "Out")
  
  
  return(dat %>%  
           ggplot(aes(x = log2FoldChange_all, y = get(paste0("log2FoldChange", y)), color = threshold, shape = threshold)) + 
           geom_point(alpha = 0.6, size = 2) + geom_abline(slope=1, intercept = 0) +
           xlab("Log2 Fold-Change in Full Group") + theme_bw() + labs(color = "", shape = "") +
           scale_color_manual(values = c("hotpink3", "lightskyblue3")) +
           scale_shape_manual(values = c(15,17)) +
           theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
                 axis.title = element_text(size = rel(1.25))))
}


# list of samples for the sensitvity
nosmoke <- all_md$Sample.ID[which(all_md$Smoking.Status == "Never")]
males   <- all_md$Sample.ID[which(all_md$Birth.sex == "Male")]
females <- all_md$Sample.ID[which(all_md$Birth.sex != "Male")]
oldies  <- all_md$Sample.ID[which(all_md$Age.at.sampling >= 38)]


# dx analyses 
all_ddr <- runDDS_DX(all_md$Sample.ID)
# sensitivity for nonsmokers
nos_ddr <- runDDS_DX(nosmoke)
# sensitivity for male only
mal_ddr <- runDDS_DX(males, sex = F) 
# sensitivity for female only
fem_ddr <- runDDS_DX(females, sex = F) 
# sensitivity for oldies 
old_ddr <- runDDS_DX(oldies, age = T) 
# run with age in the model
age_ddr <- runDDS_DX(all_md$Sample.ID, age = T)


# sensitivity for comparing male sarc NP to male sarc P
all_p_ddr <- runDDS_pNp(all_md$Sample.ID)
mal_p_ddr <- runDDS_pNp(males, sex = F) 

# combine together the sensitivities for plotting
dx_sen <- merge(nos_ddr, mal_ddr, by = "gene", suffixes = c("_noSmoke", "_males"))
dx_sen <- merge(dx_sen,  fem_ddr, by = "gene")
dx_sen <- merge(dx_sen,  all_ddr, by = "gene", suffixes = c("_females", "_all"))
dx_sen <- merge(dx_sen,  old_ddr, by = "gene")
dx_sen <- merge(dx_sen,  age_ddr, by = "gene", suffixes = c("_old", "_age"))
#colnames(dx_sen)[26:31] <- paste0(colnames(dx_sen)[26:31], "_old")


pNP_sen <- merge(all_p_ddr, mal_p_ddr, by = "gene", suffixes = c("_all", "_males"))


plotLog(dx_sen, "_noSmoke") + ylab("Log2 Fold-Change in NonSmokers") 
ggsave(paste0(outDir, "Dissemination/sens_noSmokers_dx.png"), width = 6, height = 4)
ggsave(paste0(outDir, "Dissemination/sens_noSmokers_dx.pdf"), width = 6, height = 4)

plotLog(dx_sen) + ylab("Log2 Fold-Change in Males") 
ggsave(paste0(outDir, "Dissemination/sens_malesOnly_dx.png"), width = 6, height = 4)
ggsave(paste0(outDir, "Dissemination/sens_malesOnly_dx.pdf"), width = 6, height = 4)

plotLog(dx_sen, "_females") + ylab("Log2 Fold-Change in Females") 
ggsave(paste0(outDir, "Dissemination/sens_femalesOnly_dx.png"), width = 6, height = 4)
ggsave(paste0(outDir, "Dissemination/sens_femalesOnly_dx.pdf"), width = 6, height = 4)

plotLog(dx_sen, "_old") + ylab("Log2 Fold-Change for 38+ years") 
ggsave(paste0(outDir, "Dissemination/sens_38plus_dx.png"), width = 6, height = 4)
ggsave(paste0(outDir, "Dissemination/sens_38plus_dx.pdf"), width = 6, height = 4)

dx_sen %>% dplyr::select("gene", ends_with("_all"), ends_with("_old")) %>% filter(padj_all < 0.1 & padj_old > 0.1) %>%
  filter(log2FoldChange_all < -2) %>% arrange(pvalue_all) %>% 
  dplyr::select(gene, log2FoldChange_all, pvalue_all, padj_all, log2FoldChange_old, pvalue_old, padj_old)

plotLog(dx_sen, "_age") + ylab("Log2 Fold-Change with Age in Model") 


plotLog(pNP_sen) + ylab("Log2 Fold-Change in Males") 
ggsave(paste0(outDir, "Dissemination/sens_malesOnly_pNP.png"), width = 6, height = 4)
ggsave(paste0(outDir, "Dissemination/sens_malesOnly_pNP.pdf"), width = 6, height = 4)

