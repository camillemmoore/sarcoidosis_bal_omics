###########################################
# Project: bulkRNA-seq validation in scRNA-Seq
# Author: Cheyret Wood
# #########################################

rm(list = ls())

library(tidyverse)
library(edgeR)
library(DESeq2)
library(openxlsx)
library(readxl)
library(ggrepel)
library(arsenal)

`%notin%` <- Negate(`%in%`)

options(future.globals.maxSize = 100*1000 * 1024^2)

setwd("P:/BRANCHES/NJHMoore_CBD_Sarcoidosis/maier_bulkRNASeq_paperAnalysis/Code")

source("processingData.R")
source("XX_formattingFunctions.R")

nSmallGroup <- 13
phCols <- c("steelblue1", "lightpink1", "red") # control, NP, P
dxCols <- c("steelblue1", "darkorange3") # control, Sarc

dataDir <- "P:/BRANCHES/NJHMoore_CBD_Sarcoidosis/maier_bulkRNASeq_paperAnalysis/DataRaw/"
outDir  <- "P:/BRANCHES/NJHMoore_CBD_Sarcoidosis/maier_bulkRNASeq_paperAnalysis/"

# get the rnaSeq data ready 
outDat <- getDataStarted(dataDir, nSmallGroup, outliers = T)
counts.filtered <- outDat[[1]]
all_md <- outDat[[2]]
fdata <- outDat[[3]]
genes <- outDat[[4]]

# get the methylation ids 
sMeth0  <- read_excel(paste0(dataDir, "sarc meth data with phenotype.xlsx"))
# get age at sampling
sMeth1 <- sMeth0 %>% mutate(Age.at.sample = difftime(`date of sample`, `Date of Birth`, unit = "weeks")/52.25) 
sMeth2 <- sMeth1 %>% filter(Phenotype %in% c("control", "NP", "P"))

sMeth <- sMeth2

tab1 <- tableby(Phenotype ~ Age.at.sample + `Birth sex` + Race + Hispanic + `Smoking Status` + `FVC pre % pred` + `FEV1 pre % pred` + `DLCO uncorr % pred`, 
                data = sMeth, test = F, total = F, numeric.stats = c("medianq1q3", "meansd"), cat.stats = c("countpct"))

# used this to get the Tab1_control_NP_P.csv
summary(tab1, digits = 0)

liliesNP <- c("29715", "42787", "34350", "35612", "34703", "31498", "20711", "21553", "30347", "13727", "14691", "42503", "39810", "39583", "33042")

sMeth[,1:3] %>% filter(`Sample ID` %in% liliesNP)
# 39583 is labeled as a P not NP
liliesNP[which(liliesNP %notin% sMeth$`Sample ID`)]
# 14691 is missing 

sMeth0[,1:3] %>% filter(`Sample ID` %in% liliesNP) %>% distinct(`Sample ID`, .keep_all = T)

cardwellIDs <- read_excel("C:/Users/woodche/Dropbox/Sarcoidosis&CBD Analysis Group/Sarc BAL CD4+ T cell manuscript/Cardwell methylation/44 sample analysis/samples_used_in_analysis.xlsx")
dim(sMeth %>% filter(`Sample ID` %in% cardwellIDs$`Sample ID`) %>% distinct(`Sample ID`)) 
dim(cardwellIDs)
test <- sMeth %>% filter(`Sample ID` %in% cardwellIDs$`Sample ID`) %>% distinct(`Sample ID`, .keep_all = T) 
table(test$Phenotype)

# get together a table one ####
sMeth3 <- sMeth2 %>% mutate(Race = factor(Race, levels = c("White", "Black"), ordered = T),
                            `Birth sex` = factor(`Birth sex`, levels = c("Male", "Female"), ordered = T),
                            Hispanic = factor(Hispanic, levels = c("Non-Hispanic", "Hispanic"), ordered = T),
                            `Smoking Status` = factor(`Smoking Status`, levels = c("Previous", "Never"), ordered = T),
                            flowData = ifelse(`Sample ID` %in% liliesNP, "included", "not included"), 
                            flowData = ifelse(Phenotype %in% c("P", "control") | flowData == "included", "included", "not included"), # all of lilies P and Control are accounted for from what we understand 
                            methylation = ifelse(`Sample ID` %in% cardwellIDs$`Sample ID`, "included", "not included"), 
                            rnaSeq   = ifelse(`Sample ID` %in% all_md$Sample.ID, "included", "not included"),
                            rna_and_meth = ifelse(rnaSeq == "included" & methylation == "included", "included", "not included")) %>%
  filter(Phenotype %in% c("P", "NP", "control"))


tab2 <- tableby(Phenotype ~ Age.at.sample + `Birth sex` + Race + Hispanic + `Smoking Status` + `FVC pre % pred` + `FEV1 pre % pred` + `DLCO uncorr % pred` + flowData + methylation, 
                data = sMeth3, test = F, total = F, numeric.stats = c("medianq1q3", "meansd"), cat.stats = c("countpct"))

# used this to get the Tab1_control_NP_P.csv
summary(tab2, digits = 0, text = T)


sMeth3 %>% select('Sample ID', Age.at.sample, Dx, Phenotype, flowData:rna_and_meth) %>%
  pivot_longer(cols = flowData:rna_and_meth) %>% 
  filter(value == "included") %>%
  mutate(name = case_when(name == "flowData" ~ "Flow Data", 
                          name == "methylation" ~ "Methylation", 
                          name == "rnaSeq" ~ "Bulk RNA-Seq", 
                          name == "rna_and_meth" ~ "In RNA-Seq and Methylation"), 
         name = factor(name, levels = c("Flow Data",  "Methylation", "Bulk RNA-Seq", "In RNA-Seq and Methylation"), ordered = T)) %>%
  ggplot(aes(x = Age.at.sample, fill = Phenotype)) + 
  geom_histogram(alpha = 0.5, position = "identity", bins = 15) + theme_bw() + 
  # use some slightly different colors
  scale_fill_manual(values = c(phCols[1], "lightsalmon3", "lightpink3")) +
  facet_wrap(~name) + labs(x = "Age at Sample (yrs)")
ggsave("../Dissemination/hist_age_pheno.png")
ggsave("../Dissemination/hist_age_pheno.pdf")

sMeth3 %>% select('Sample ID', Age.at.sample, Dx, Phenotype, flowData:rna_and_meth) %>%
  pivot_longer(cols = flowData:rna_and_meth) %>% 
  filter(value == "included") %>%
  mutate(name = case_when(name == "flowData" ~ "Flow Data", 
                          name == "methylation" ~ "Methylation", 
                          name == "rnaSeq" ~ "Bulk RNA-Seq", 
                          name == "rna_and_meth" ~ "In RNA-Seq and Methylation"), 
         name = factor(name, levels = c("Flow Data",  "Methylation", "Bulk RNA-Seq", "In RNA-Seq and Methylation"), ordered = T)) %>%
  ggplot(aes(x = Age.at.sample, fill = Dx)) + 
  geom_histogram(alpha = 0.5, position = "identity", bins = 15) + theme_bw() + 
  scale_fill_manual(values = dxCols) +
  facet_wrap(~name) + labs(x = "Age at Sample (yrs)")
ggsave("../Dissemination/hist_age_dx.png")
ggsave("../Dissemination/hist_age_dx.pdf")

sMeth3 %>% dplyr::select('Sample ID', Age.at.sample, Dx, Phenotype, rnaSeq) %>%
  filter(rnaSeq == "included") %>%
  ggplot(aes(x = Age.at.sample, fill = Dx)) + 
  geom_histogram(alpha = 0.5, position = "identity", bins = 15) + theme_bw() + 
  scale_fill_manual(values = dxCols) + labs(x = "Age at Sample (yrs)")
ggsave("../Dissemination/hist_rna_age_dx.png")
ggsave("../Dissemination/hist_rna_age_dx.pdf")


# Create DESeq2 Dataset and Model
dds <- DESeqDataSetFromMatrix(countData = counts.filtered, colData = all_md,
                              design = ~ Birth.sex + Dx, rowData = fdata)
dds <- DESeq(dds)

# Simple DESeq2 Differential Expression Analysis #####

# Get results for the comparison of the groups
Sarc_v_Control <- results(dds)
Sarc_v_Control <- Sarc_v_Control[order(Sarc_v_Control$padj),]
# Add in the gene symbols
Sarc_v_Control$gene <- genes[rownames(Sarc_v_Control)]
Sarc_v_Control_df <- formatDat(Sarc_v_Control)

sarcVCtrlPlot <- Sarc_v_Control_df %>% 
  mutate(threshold = ifelse(padj < 0.05, "FDR < 0.05", "FDR > 0.05")) %>% 
  filter(is.na(threshold) == F) 

probPoints <- sarcVCtrlPlot %>% filter(padj < 0.05 & log2FoldChange > 0 & log2FoldChange < 1 )

sigColor <- "black"

sarcVCtrlPlot %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue), color = threshold)) + 
  geom_point() + theme_bw() + xlim(min(sarcVCtrlPlot$log2FoldChange), min(sarcVCtrlPlot$log2FoldChange)*-1) +
  geom_text_repel(data = sarcVCtrlPlot %>% filter(padj < 0.05 & gene %notin% probPoints$gene & log2FoldChange >= 1),
                  aes(label = gene), xlim = c(2, 3.7), color = sigColor,
                  max.overlaps = 20, seed = 400, min.segment.length = 0.1, size = 3.3) +
  geom_text_repel(data = sarcVCtrlPlot %>% filter(padj < 0.05 & gene %notin% probPoints$gene & log2FoldChange <= 0),
                  aes(label = gene), xlim = c(-3.7, -1), color = sigColor,
                  max.overlaps = 20, seed = 400, min.segment.length = 0.1, size = 3.3) +
  geom_text_repel(data = sarcVCtrlPlot %>% filter(padj < 0.05 & gene %in% probPoints$gene & -log10(pvalue) < 4.5),
                  aes(label = gene), xlim = c(-1, 2), ylim = c(0, 4), color = sigColor, force = 5, 
                  max.overlaps = 20, seed = 400, min.segment.length = 0.1, size = 3.3) +
  geom_text_repel(data = sarcVCtrlPlot %>% filter(padj < 0.05 & gene %in% probPoints$gene & -log10(pvalue) >= 4.5),
                  aes(label = gene), xlim = c(-2, 2), ylim = c(5, 8), color = sigColor, force = 4, 
                  max.overlaps = 20, seed = 400, min.segment.length = 0.1, size = 3.3) +
  labs(x = "log2(Fold Change)", y = "-log10(p-value)") +
  scale_color_manual(values = c(alpha(sigColor, 0.5), alpha("grey", 0.5))) +
  theme(legend.position = "none") 

ggsave(filename = "../Dissemination/SarcVControl_volcanoPlot_bw.png", width = 8, height = 5, units = "in",  dpi = "retina")
ggsave(filename = "../Dissemination/SarcVControl_volcanoPlot_bw.pdf", width = 8, height = 5, units = "in",  dpi = "retina")
