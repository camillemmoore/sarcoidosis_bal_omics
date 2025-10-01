###########################################
# Project: bulkRNA-seq validation in scRNA-Seq
# Author: Cheyret Wood
# #########################################

rm(list = ls())

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(arsenal))

`%notin%` <- Negate(`%in%`)

options(future.globals.maxSize = 100*1000 * 1024^2)

setwd("P:/BRANCHES/NJHMoore_CBD_Sarcoidosis/maier_bulkRNASeq_paperAnalysis/Code")

source("processingData.R")
source("XX_formattingFunctions.R")

nSmallGroup <- 13

dataDir <- "P:/BRANCHES/NJHMoore_CBD_Sarcoidosis/maier_bulkRNASeq_paperAnalysis/DataRaw/"
outDir  <- "P:/BRANCHES/NJHMoore_CBD_Sarcoidosis/maier_bulkRNASeq_paperAnalysis/"
cardwellIDs <- read_excel("C:/Users/woodche/Dropbox/Sarcoidosis&CBD Analysis Group/Sarc BAL CD4+ T cell manuscript/Cardwell methylation/44 sample analysis/samples_used_in_analysis.xlsx")
sMeth0  <- read_excel(paste0(dataDir, "sarc meth data with phenotype.xlsx")) 
sarc0   <- read_excel(paste0(dataDir, "sarc data with PFTs corrected.xlsx")) 
rnaSeq  <- read.xlsx(paste0(dataDir, "bal_cd4+_rnaseqc_rnaseq_reads.xlsx"))
# all of lilies controls and P are accounted for from what we understand 
liliesNP <- c("29715", "42787", "34350", "35612", "34703", "31498", "20711", "21553", "30347", "13727", "14691", "42503", "39810", "39583", "33042")

# add the additional methylation sample
cardwellIDs <- rbind(cardwellIDs, c("49618", NA))

# get the RNA-Seq clinical data 
source("processingData.R")
dataDir <- "P:/BRANCHES/NJHMoore_CBD_Sarcoidosis/maier_bulkRNASeq_paperAnalysis/DataRaw/"
nSmallGroup <- 13
outDat <- getDataStarted(dataDir, nSmallGroup, outliers = T)
all_md <- outDat[[2]]

# reduce to columns of interest 
sMeth1 <- sMeth0 %>% select(`Sample ID`, Dx, Phenotype, `date of sample`, `Date of Birth`, `Birth sex`, Race, Hispanic, `Smoking Status`, ends_with("pred"),  `study #`,  `study # 2`, `Sarc Proj`)
sarc1  <- sarc0  %>% select(`Sample ID`, `Dx of data for Ivana`, `date of sample`, `Date of Birth`, `Birth sex`, Race, Hispanic, `Smoking Status`, ends_with("pred"))

# add in age at sample
sMeth2 <- sMeth1 %>% mutate(Age.at.sample = difftime(`date of sample`, `Date of Birth`, unit = "weeks")/52.25) 
sarc2  <- sarc1  %>% mutate(Age.at.sample = difftime(`date of sample`, `Date of Birth`, unit = "weeks")/52.25) 

dat <- merge(sMeth2, sarc2, by = colnames(sMeth2)[which(colnames(sMeth2) %in% colnames(sarc2))], all = T)
# reduce to unique rows
dat <- dat %>% distinct(.keep_all = T)

# check the overlap with lily ####
#dat[,c("Sample ID", "Dx", "Phenotype")] %>% filter(`Sample ID` %in% liliesNP)
# 39583 is labeled as a P not NP
#liliesNP[which(liliesNP %notin% dat$`Sample ID`)]
# 14691 is missing 
# add in fix for lilies ID
dat %>% filter(`Sample ID` %in% "16462")

#sMeth0[,1:3] %>% filter(`Sample ID` %in% liliesNP) %>% distinct(`Sample ID`, .keep_all = T)

# Methylation Breakdown ####
#dim(dat %>% filter(`Sample ID` %in% cardwellIDs$`Sample ID`) %>% distinct(`Sample ID`)) 
#dim(cardwellIDs)
#test <- sMeth %>% filter(`Sample ID` %in% cardwellIDs$`Sample ID`) %>% distinct(`Sample ID`, .keep_all = T) 
#table(test$Phenotype)

# get the rnaseq ids ####
rnaSubId <- all_md$Sample.ID

# get together a table one ####
sMeth3 <- dat %>% mutate(Race = factor(Race, levels = c("White", "Black"), ordered = T),
                         `Birth sex` = factor(`Birth sex`, levels = c("Male", "Female"), ordered = T),
                         Hispanic = factor(Hispanic, levels = c("Non-Hispanic", "Hispanic"), ordered = T),
                         `Smoking Status` = factor(`Smoking Status`, levels = c("Previous", "Never"), ordered = T),
                         flowData = ifelse(`Sample ID` %in% liliesNP, "included", "not included"), 
                         flowData = ifelse(Phenotype %in% c("P", "control") | flowData == "included", "included", "not included"), # all of lilies P and Control are accounted for from what we understand 
                         methylation = ifelse(`Sample ID` %in% cardwellIDs$`Sample ID`, "included", "not included"), 
                         rnaSeq   = ifelse(`Sample ID` %in% rnaSubId, "included", "not included"),
                         rna_and_meth = ifelse(rnaSeq == "included" & methylation == "included", "included", "not included")) %>%
  filter(Phenotype %in% c("P", "NP", "control"))

# save out the sample IDs
write.csv(sMeth3[,c("Sample ID", "methylation", "rnaSeq")], row.names = F, file = paste0(outDir, "DataProcessed/Sarc_BAL_CD4+Tcell_sampleID_by_datatype.csv"))


#temp <- tableby(Phenotype ~ Age.at.sample + `Birth sex` + Race + Hispanic + `Smoking Status` + `FVC pre % pred` + `FEV1 pre % pred` + `DLCO uncorr % pred` + methylation + rnaSeq + rna_and_meth, 
#                data = sMeth3, test = T, total = F)
#summary(temp, digits = 1, text = T, pfootnote = T)
 
# the overall table 1 
tab2 <- tableby(Phenotype ~ anova(Age.at.sample, "medianq1q3") + fe(`Birth sex`) + fe(Race) + fe(Hispanic) + fe(`Smoking Status`) + `FVC pre % pred` + `FEV1 pre % pred` + `DLCO uncorr % pred` + 
                  methylation + rnaSeq + rna_and_meth, 
                data = sMeth3, test = T, total = F, numeric.stats = c("meansd"), cat.stats = c("countpct"))
summary(tab2, digits = 1, pfootnote = T)

# get the p-values for the lung functions 
t.test(`FVC pre % pred` ~ Phenotype, data = sMeth3[which(sMeth3$Phenotype != "control"),])
t.test(`FEV1 pre % pred` ~ Phenotype, data = sMeth3[which(sMeth3$Phenotype != "control"),])
t.test(`DLCO uncorr % pred` ~ Phenotype, data = sMeth3[which(sMeth3$Phenotype != "control"),])

# really getting into the sample numbers 
#tab3 <- tableby(Phenotype ~ methylation + rnaSeq + rna_and_meth, #Age.at.sample + `Birth sex` + Race + Hispanic + `Smoking Status` + `FVC pre % pred` + `FEV1 pre % pred` + `DLCO uncorr % pred` 
#                data = sMeth3 %>% filter(methylation == "included"), test = F, total = T, numeric.stats = c("medianq1q3", "meansd"), cat.stats = c("countpct"))
#summary(tab3, digits = 0, text = T)

#tab4 <- tableby(Phenotype ~ Age.at.sample + `Birth sex` + methylation + rnaSeq + rna_and_meth, # + Race + Hispanic + `Smoking Status` + `FVC pre % pred` + `FEV1 pre % pred` + `DLCO uncorr % pred` 
#                data = sMeth3, test = F, total = T, numeric.stats = c("medianq1q3"), cat.stats = c("countpct"))
#summary(tab4, digits = 0, text = T)

# code for original table 1 ####

#sMeth2 <- sMeth2 %>% filter(Phenotype %in% c("control", "NP", "P"))
#sMeth <- sMeth2

#tab1 <- tableby(Phenotype ~ Age.at.sample + `Birth sex` + Race + Hispanic + `Smoking Status` + `FVC pre % pred` + `FEV1 pre % pred` + `DLCO uncorr % pred`, 
#                data = sMeth, test = F, total = F, numeric.stats = c("medianq1q3", "meansd"), cat.stats = c("countpct"))

# used this to get the Tab1_control_NP_P.csv
#summary(tab1, digits = 0)


