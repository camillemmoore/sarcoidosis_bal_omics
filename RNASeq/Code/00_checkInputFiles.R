###########################################
# Project: maier bulk rnaseq redo
# Author: Cheyret Wood
# Date: 11/30/2023
# #########################################

## linux code I want to keep handy 
#PATH="$HOME/miniconda/bin:$PATH"
#srun --pty /bin/bash

rm(list = ls())

# Load Libraries
library(tidyverse)
library(skimr)

`%notin%` <- Negate(`%in%`)

options(future.globals.maxSize = 100*1000 * 1024^2)

dataDir <- "/Moore/proj/maier_bulkRNASeq/Maier_Repro_Bulk_Results/Data/"

rnaSeq <- read.xlsx(paste0(dataDir, "bal_cd4+_rnaseqc_rnaseq_reads.xlsx"))
sData  <- read.xlsx(paste0(dataDir, "sarc data with PFTs corrected.xlsx"))
sMeth  <- read.xlsx(paste0(dataDir, "sarc meth data with phenotype.xlsx"))

cat("\n dim(rnaSeq) \n")
dim(rnaSeq) # 59,841 X 37
cat("\n")

cat("\n head(rnaSeq[,1:15]) \n ")
head(rnaSeq[,1:15])
cat("\n")

cat("\n dim(sData) \n")
dim(sData) # 52 X 34 
cat("\n")

cat("\n head(sData) \n")
head(sData)
cat("\n")

cat("\n dim(sMeth) \n")
dim(sMeth) # 52 X 36
cat("\n")

cat("\n head(sMeth) \n")
head(sMeth)
cat("\n")

# check what's unique to sMeth ####

# same samples 
identical(sData$Sample.ID[order(sData$Sample.ID)], sMeth$Sample.ID[order(sMeth$Sample.ID)])

# unique column names ("Dx", "Phenotype", "study.#.2")
colnames(sMeth)[which(colnames(sMeth) %notin% colnames(sData))]
# unique column names (Dx.of.data.for.Ivana)
colnames(sData)[which(colnames(sData) %notin% colnames(sMeth))]

cat("\n summary(as.factor(sMeth$Dx)) \n")
summary(as.factor(sMeth$Dx))
cat("\n")
# ? sarc  Control     Sarc Sarc/CBD
#      1       20       30        1

cat("\n summary(as.factor(sMeth$Phenotype)) \n")
summary(as.factor(sMeth$Phenotype))
cat("\n")
# CBD     control   non bx proven sarc sarc   NP    P
#   1          20                         1   14   16

cat("\n summary(as.factor(sMeth$`study.#.2`)) \n")
summary(as.factor(sMeth$`study.#.2`))
cat("\n")
#  P   NA's
#  3     49

cat("\n summary(as.factor(sData$Dx.of.data.for.Ivana)) \n")
summary(as.factor(sData$Dx.of.data.for.Ivana))
cat("\n")
# Control   Sarc  Sarc/CBD
#      20     31          1

# see if I can merge them together 
sMethDat <- merge(sMeth, sData, by = colnames(sMeth)[which(colnames(sMeth) %in% colnames(sData))])

cat("\n dim(sMethDat) \n")
dim(sMethDat) # 52 X 37
cat("\n")

# now it's all together

# look at the skim results
cat("\n skim(sMethDat) \n")
skim(sMethDat)
cat("\n")


# RNA-Seq ####

# get the list of colnames
rnaSamp <- colnames(rnaSeq[,3:ncol(rnaSeq)])
# remove the extra labeling
rnaSamp1 <- str_split(rnaSamp, pattern = "_", simplify = T)[,2]
# check they are unique ids
identical(length(rnaSamp1), length(unique(rnaSamp1)))

#fix colnames 
colnames(rnaSeq)[3:ncol(rnaSeq)] <- rnaSamp1

# now pull colnames to just the subject id
rnaSubId <- gsub("B", "", rnaSamp1)

# we have the information for all 37 samples
cat("\n dim(sMethDat[which(sMethDat$Sample.ID %in% rnaSubId),]) \n")
dim(sMethDat[which(sMethDat$Sample.ID %in% rnaSubId),])
cat("\n")

# number of sarc subject with RNA-Seq information
summary(as.factor(sMethDat[which(sMethDat$Sample.ID %in% rnaSubId),"Dx"]))
#? sarc  Control   Sarc  Sarc/CBD
#     1       15     20         1
