library(sesame) #DNAm preprocessing
library(dplyr) #tidying
library(matrixStats) #matrix manipulation
library(stringr) #string manipulation
library(limma) #linear regression
library(ggplot2) #visualization
library(ggrepel)
library(EnhancedVolcano) #visualization
library(bacon) #p-value adjustment for bias & inflation
library(ENmix) #principal component regression
library(qqman) #qq and manhattan plotting
library(gprofiler2) #enrichment analysis
library(randomForest) #prediction
library(data.table)
library(ExperimentHub)
set.seed(11) #set seed for reproducibility

setwd(".")
#cache relevant annotation data from ExperimentHub
sesameDataCache() #this is only necessary once per sesame installation

#read in phenotype data
pd <- read.csv("../meth_phenotype.csv")
par(mfrow=c(1,1))
boxplot(age~DX,data=pd, main="Distribution of Age by DX")

#prepare array annotation data 
ann = as.data.frame(fread("EPIC.hg38.manifest.tsv.gz"))

#detect idat files in your working directory
files <- searchIDATprefixes(getwd())

#run QC on idats
qc_pre <- openSesame(files, prep="", func=sesameQC_calcStats, BPPARAM = BiocParallel::SnowParam(2))
qc_pre[1]
saveRDS(qc_pre, "qc_pre.rds")

#preprocess idats
sdf <- openSesame(files, func = NULL, BPPARAM = BiocParallel::SnowParam(2))
saveRDS(sdf, "sdf.rds")

#run QC on preprocessed data
qc_post <- openSesame(sdf, func = sesameQC_calcStats, BPPARAM = BiocParallel::SnowParam(2))
saveRDS(qc_post, "qc_post.rds")

#add additional QC metrics: bisulfite conversion, predicted sex, karyotype, & ethnicity
qc_pre_df <- do.call(rbind, lapply(qc_pre, as.data.frame))
qc_post_df <- do.call(rbind, lapply(qc_post, as.data.frame))

#Calculate bisulfite conversion score
qc_pre_df$gct <- sapply(sdf, bisConversionControl)
qc_post_df$gct <- qc_pre_df$gct

#infer sex
qc_pre_df$sex_predicted <- sapply(sdf, inferSex)
qc_post_df$sex_predicted <- qc_pre_df$sex_predicted

#infer karyotype
qc_pre_df$karyotype_predicted <- sapply(sdf, inferSexKaryotypes)
qc_post_df$karyotype_predicted <- qc_pre_df$karyotype_predicted

#infer ethnicity
qc_pre_df$ethnicity_predicted <- sapply(sdf, inferEthnicity)
qc_post_df$ethnicity_predicted <- qc_pre_df$ethnicity_predicted

write.csv(qc_pre_df, "qc_pre.csv")
write.csv(qc_post_df, "qc_post.csv")

betas <- do.call(cbind, lapply(sdf, getBetas))
saveRDS(betas, "betas_raw.rds")

#remove non-cpg probes
betas <- betas[-grep("ch", rownames(betas)),]
betas <- betas[-grep("rs", rownames(betas)),]
betas <- betas[-grep("ctl", rownames(betas)),]

#remove probes flagged in all samples
missing <- apply(betas, 1, function(x){sum(is.na(x))})
betas <- betas[missing == 0,] 
nrow(betas)

saveRDS(betas, "betas_cleaned.rds")

#convert to M values
Ms <- BetaValueToMValue(betas)
saveRDS(Ms, "Ms.rds")
Ms = readRDS("Ms.rds")

#remove sex chromosome probes. Otherwise PCA will strongly seperate by sex. 
Ms_sub <- Ms[!rownames(Ms) %in% ann[ann$CpG_chrm == "chrX",]$Probe_ID & !rownames(Ms) %in% ann[ann$CpG_chrm == "chrX",]$Probe_ID, ]
pca <- prcomp(t(Ms_sub))
summary(pca)

#add PCs to phenotype data
pd <- merge(pd, pca$x, by.x = "sentrix", by.y = "row.names")
rownames(pd) <- pd$sentrix

#plot top PCs
library(MetBrewer)
pal=met.brewer(name="Austria",n=7,type="discrete")
ggplot(pd, aes(PC1, PC2)) + geom_point() + theme_bw()
labplot = ggplot(pd, aes(PC1, PC2),label=sentrix) + geom_point(color = dplyr::case_when(pd$DX == "Control" ~ 'blue', pd$DX == "Sarc" ~ pal[1]), size = 3, alpha = 0.8) + theme_bw()
### geom_label_repel
labplot + 
  geom_label_repel(aes(label = sentrix),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   ) +
  theme_classic()

#save current phenotype data
write.csv(pd, "WorkingPD.csv")

#build a dataframe of relevant covariates. We have to convert everything to numeric variables. 
covs = pd %>% dplyr::select(DX, Progressor, Age.at.diagnosis, Corrected_sex, Race, Smoking.Status, Sentrix.Barcode,Sentrix.Position)
covs$sarc <- 0
covs[covs$DX != "Control", ]$sarc <- 1

covs$prog <- 2
covs[covs$Progressor == "C", ]$prog <- 0
covs[covs$Progressor == "NP", ]$prog <- 1

covs[covs$Corrected_sex == "Male", ]$Corrected_sex <- 1
covs[covs$Corrected_sex == "Female", ]$Corrected_sex <- 0
covs$Sentrix.Barcode = as.factor(covs$Sentrix.Barcode)
covs$Sentrix.Position = as.factor(covs$Sentrix.Position)

covs[covs$Smoking.Status == "Previous", ]$Smoking.Status <- 1
covs[covs$Smoking.Status == "Never", ]$Smoking.Status <- 0

covs[covs$Race == "Black", ]$Race <- 1
covs[covs$Race == "White", ]$Race <- 0

#covs <- covs %>% dplyr::select(sarc, age, sex, smoke_current, smoke_former)

# There are more arrays than I have clinical data for
Ms_sub_filtered = Ms_sub[ , colnames(Ms_sub) %in% rownames(covs) ]
Ms_filtered = Ms[ , colnames(Ms) %in% rownames(covs) ]

# Write out M values for PEER
write.table(Ms_filtered, "Final_Ms_sarc.csv", sep = ",", row.names = T,quote=F)
write.table(covs, "meth_cov.csv", sep = ",", row.names = T,quote=F)

pcrplot(Ms_sub_filtered, covs, npc = 10, subset = F) #plots to pcr_diag.jpg. Note we are using data without sex chromosome probes for this qc step; sex must always be adjusted for if sex chromosomes are included in analysis. 

#specify model
#design <- model.matrix( ~ DX + age + Birth.sex + Smoking.Status + Race, pd)
design <- model.matrix( ~ DX + Sentrix.Barcode + Sentrix.Position + Corrected_sex + Race, covs)
#design <- model.matrix( ~ DX + Sentrix.Barcode + Sentrix.Position + Birth.sex + Race + age, covs)

#fit model to data
fit <- lmFit(Ms_filtered, design)
colnames(design)
contrasts <- makeContrasts(diagnosis = "DXSarc", levels = design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)
summary(decideTests(fit))

res <- limma::topTable(fit, num = Inf, coef = "diagnosis")
bake <- bacon(res$t)
#plot bacon results. Histogram is of test statistics. Black line represent normal distribution. Red line represents bacon-derived null distribution. Confirm that the red line better fits the test statistic distribution. 
plot(bake)
traces(bake) #parameter estimates should stabilize over iterations; high stochasticity can indicate bacon not performing well. 

#update results with bacon-adjusted test statistics
res$bacon_p <- bacon::pval(bake)
res$bacon_q <- p.adjust(res$bacon_p, method = "BH")
nrow(res[res$bacon_q <= .05, ])

medianbeta_sarc <- rowMedians(betas[,colnames(betas) %in% rownames(pd[pd$DX == "Sarc", ])])
medianbeta_control <- rowMedians(betas[,colnames(betas) %in% rownames(pd[pd$DX == "Control", ])])
medianbeta_dif <- medianbeta_sarc - medianbeta_control
beta_df <- cbind(medianbeta_sarc, medianbeta_control, medianbeta_dif)
rownames(beta_df) <- rownames(betas)
beta_df <- beta_df[rownames(res),]
res <- cbind(res, beta_df)

res <- merge(res,ann, by = "row.names")
rownames(res) <- res$Row.names
res$Row.names <- NULL

res <- res[order(res$bacon_p),]
write.csv(res, "DifferentialMethylation_FullResults.csv")
write.csv(res[res$bacon_p <= .05,], "DifferentialMethylation_SignificantResults.csv")

## Visualization
EnhancedVolcano(res, lab = res$gene, x = "logFC", y = "bacon_q", pCutoff = 0.05) 

# ======================================== #
# Get PEER factors
# ======================================== #
library(peer)
library(qtl)

print("starting")
#read in sample info (samples in rows)
cov=read.csv("meth_cov.csv",stringsAsFactors=F)
cov=cov[order(rownames(cov)),]
print(head(cov))

#read in methylation data (samples in columns)
dat=read.csv("Final_Ms_sarc.csv",row.names=1)
colnames(dat)=gsub("X","",colnames(dat))
dat=dat[,order(colnames(dat))]
print(dim(dat))
print(dat[c(1:10),c(1:10)])

#transpose data
dat=as.matrix(dat)
dat=t(dat)

#set up model for PEER with covariates of interest. 
cov_mod=as.matrix(model.matrix(~Age.at.diagnosis+Corrected_sex+Race,cov))
print(dim(dat))
print(head(cov_mod))

#IMPORTANT NOTE: PEER outputs lack column and row names. For the residuals file, these will correspond with names of the final methylation data matrix. 
#I recommend saving the colnames and rownames of dat here so you can easily incorporate with the residuals file. 
write.csv(rownames(dat), "peer_20k_rownames.csv",row.names=TRUE,quote=FALSE)
write.csv(dat, "peer_20k_dat.csv",row.names=TRUE,quote=FALSE)
write.csv(colnames(dat), "peer_20k_colnames.csv",row.names=TRUE,quote=FALSE)
write.csv(cov_mod, "peer_20k_cov_mod.csv",row.names=TRUE,quote=FALSE)

#set up PEER model
model = PEER()
Nmax_iterations = 750 #100 iteration convergence cutoff
PEER_setPhenoMean(model,as.matrix((dat)))
dim(PEER_getPhenoMean(model))
PEER_setPriorAlpha(model,0.001,0.1)
PEER_setPriorEps(model,0.1,10.)
PEER_setCovariates(model, cov_mod)
PEER_setNmax_iterations(model,Nmax_iterations)
PEER_setNk(model,20) #Calculate 20 PEER factors
PEER_update(model)
factors = PEER_getX(model)
Yc = PEER_getResiduals(model)
Alpha = PEER_getAlpha(model)
weights=PEER_getW(model)
write.csv(factors, "peer_20k_factors.csv",row.names=TRUE,quote=FALSE) #PEER factors
Yc=t(Yc) #will read and write much faster if loci in rows instead of columns. 
write.csv(Yc, "peer_20k_residuals.csv",row.names=TRUE,quote=FALSE) #data residuals after PEER and covariate regression
write.csv(Alpha, "peer_20k_alpha.csv",row.names=TRUE,quote=FALSE)
write.csv(weights, "peer_20k_weights.csv",row.names=TRUE,quote=FALSE)


peer_factors = read.csv(file="redo_peer_20k_factors.csv",head=TRUE,sep=",",check.names=FALSE,stringsAsFactors = FALSE)
actual_pf = peer_factors[,c(6:25)]

# Make a scree plot
pc.cr <- princomp(actual_pf, cor = TRUE)
screeplot(pc.cr,type="line",main="Scree Plot of 20 PEER Factors",npcs=20)

# PCRA
covs=read.csv("meth_cov.csv",stringsAsFactors=F)
test_pheno = covs[,c(1:8)]
test_pheno$Sentrix.Barcode = as.factor(test_pheno$Sentrix.Barcode)

rownames(actual_pf) = rownames(test_pheno)
npc = 100
p <- lmmatrix(as.matrix(actual_pf),test_pheno)
plotp("pcr_diag_peer_factors.jpg", 
      p, colnames(p), length(colnames(actual_pf)), title="Principal Component Regression Analysis on 20 PEER factors") 

# ======================================== #
# Back to Methylation workflow
# ======================================== #
dat = Ms_filtered
# read file with covariates combined with PEER factors
meth_cov_peer <- read.csv("redo_meth_cov_peer.csv",stringsAsFactors = TRUE,row.names=1)
ann <- as.data.frame(sesameAnno_get("Anno/EPIC/EPIC.hg38.manifest.rds"))
ann = as.data.frame(fread("EPIC.hg38.manifest.tsv.gz"))
betas = readRDS("betas_cleaned.rds")

pd <- read.csv("../meth_phenotype.csv")
Ms = readRDS("Ms.rds")
Ms_sub <- Ms[!rownames(Ms) %in% ann[ann$CpG_chrm == "chrX",]$Probe_ID & !rownames(Ms) %in% ann[ann$CpG_chrm == "chrX",]$Probe_ID, ]
pca <- prcomp(t(Ms_sub))
summary(pca)
pd <- merge(pd, pca$x, by.x = "sentrix", by.y = "row.names")
rownames(pd) <- pd$sentrix

#### specify model - Sarc vs. Control ####
design <- model.matrix( ~ DX + Corrected_sex + Race + P2, meth_cov_peer)

#fit model to data
fit <- lmFit(dat, design)
colnames(design)
contrasts <- makeContrasts(diagnosis = "DXSarc", levels = design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)
summary(decideTests(fit))

res <- limma::topTable(fit, num = Inf, coef = "diagnosis")
bake <- bacon(res$t)
#plot bacon results. Histogram is of test statistics. Black line represent normal distribution. Red line represents bacon-derived null distribution. Confirm that the red line better fits the test statistic distribution. 
plot(bake)
traces(bake) #parameter estimates should stabilize over iterations; high stochasticity can indicate bacon not performing well. 

#update results with bacon-adjusted test statistics
res$bacon_p <- bacon::pval(bake)
res$bacon_q <- p.adjust(res$bacon_p, method = "BH")
nrow(res[res$bacon_q <= .05, ])

# QQ plots before and after bacon
qq(res$P.Value)
qq(res$bacon_p)

qq(res$adj.P.Val)
qq(res$bacon_q)

medianbeta_sarc <- rowMedians(betas[,colnames(betas) %in% rownames(pd[pd$DX == "Sarc", ])])
medianbeta_control <- rowMedians(betas[,colnames(betas) %in% rownames(pd[pd$DX == "Control", ])])
medianbeta_dif <- medianbeta_sarc - medianbeta_control
beta_df <- cbind(medianbeta_sarc, medianbeta_control, medianbeta_dif)
rownames(beta_df) <- rownames(betas)
beta_df <- beta_df[rownames(res),]
res <- cbind(res, beta_df)

res <- merge(res,ann, by.x = "row.names",by.y = "Probe_ID")
rownames(res) <- res$Row.names
res$Row.names <- NULL

res <- res[order(res$bacon_p),]
head(res,7)
write.csv(res, "redo_SvC_DifferentialMethylation_FullResults_peer_correctedsex.csv")
write.csv(res[res$bacon_p <= .05,], "redo_SvC_DifferentialMethylation_SignificantResults_peer_correctedsex.csv")

## Visualization
par(mfrow=c(1,1))
EnhancedVolcano(res, lab = res$gene, x = "logFC", y = "bacon_q", pCutoff = 0.05) 

res$chrom <- gsub("chr", "", res$seqnames)
res[res$chrom == "X", ]$chrom <- "23"
res[res$chrom == "Y", ]$chrom <- "24"
res$chrom <- as.numeric(res$chrom)

manhattan(res, chr = "chrom", bp = "start", p = "bacon_q", snp = "gene", suggestiveline = -log10(.05), genomewideline = 0, annotatePval = 0.05)

# Pathway Analysis
#subset to significant associations
sig <- res[res$bacon_q <= .05, ]

#perform enrichment on hypermethylated and hypomethylated probes using sesame
enrich_hyper <- testEnrichment(rownames(sig[sig$logFC > 0, ]))
enrich_hypo <- testEnrichment(rownames(sig[sig$logFC < 0, ]))

#plot significant associations for hypermethylated hits
KYCG_plotEnrichAll(enrich_hyper)

KYCG_plotDot(enrich_hyper, n_max=20)

KYCG_plotDot(enrich_hyper[enrich_hyper$group == "KYCG.EPIC.chromHMM.20211020",], n_max = 20) 
#gprofiler2 enrichment
genes_hyper <- sesameData_getGenesByProbes(rownames(sig[sig$logFC > 0, ]))
## Platform set to: EPIC
genes_hypo <- sesameData_getGenesByProbes(rownames(sig[sig$logFC < 0, ]))
## Platform set to: EPIC
enrich_dbs <- gost(query = list("hypo" = genes_hypo$gene_name, "hyper" = genes_hyper$gene_name))

#plot associations
gostplot(enrich_dbs)

#### specify model - Progressive vs. Remitting (P vs. NP) #### 
P_vs_NP = subset(meth_cov_peer, Progressor %in% c("NP","P"))
P_vs_NP$Progressor = factor(P_vs_NP$Progressor,levels=c("NP","P"))
test_dat = dat[ , colnames(dat) %in% rownames(P_vs_NP) ]
design <- model.matrix( ~ Progressor + Corrected_sex + Race + P2, P_vs_NP)

write.csv(P_vs_NP, "PvNP_metadata.csv")
#only 1 female progressor)

#fit model to data
fit <- lmFit(test_dat, design)
colnames(design)
contrasts <- makeContrasts(diagnosis = "ProgressorP", levels = design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)
summary(decideTests(fit))

res <- limma::topTable(fit, num = Inf, coef = "diagnosis")
bake <- bacon(res$t)
#plot bacon results. Histogram is of test statistics. Black line represent normal distribution. Red line represents bacon-derived null distribution. Confirm that the red line better fits the test statistic distribution. 
plot(bake)
traces(bake) #parameter estimates should stabilize over iterations; high stochasticity can indicate bacon not performing well. 

#update results with bacon-adjusted test statistics
res$bacon_p <- bacon::pval(bake)
res$bacon_q <- p.adjust(res$bacon_p, method = "BH")
nrow(res[res$bacon_q <= .05, ])

medianbeta_sarc <- rowMedians(betas[,colnames(betas) %in% rownames(P_vs_NP[P_vs_NP$Progressor == "P", ])])
medianbeta_control <- rowMedians(betas[,colnames(betas) %in% rownames(P_vs_NP[P_vs_NP$Progressor == "NP", ])])
medianbeta_dif <- medianbeta_sarc - medianbeta_control
beta_df <- cbind(medianbeta_sarc, medianbeta_control, medianbeta_dif)
rownames(beta_df) <- rownames(betas)
beta_df <- beta_df[rownames(res),]
res <- cbind(res, beta_df)

res <- merge(res,ann, by = "row.names")
rownames(res) <- res$Row.names
res$Row.names <- NULL

res <- res[order(res$bacon_p),]
head(res,7)

par(mfrow=c(1,1))
EnhancedVolcano(res, lab = res$gene, x = "logFC", y = "bacon_q", pCutoff = 0.05) 
res$chrom <- gsub("chr", "", res$seqnames)
res[res$chrom == "X", ]$chrom <- "23"
res[res$chrom == "Y", ]$chrom <- "24"
res$chrom <- as.numeric(res$chrom)
manhattan(res, chr = "chrom", bp = "start", p = "bacon_q", snp = "gene", suggestiveline = -log10(.05), genomewideline = 0, annotatePval = 0.05)

write.csv(res, "PvNP_DifferentialMethylation_FullResults_peer_correctedsex.csv")
write.csv(res[res$bacon_p <= .05,], "PvNP_DifferentialMethylation_SignificantResults_peer_correctedsex.csv")

# Boxplots of random top P vs. NP cpgs
random_top_cpgs = c("cg23406740","cg15509024","cg16218368","cg04935278","cg21016177","cg16970087")
test = subset(dat,rownames(dat) %in% random_top_cpgs)
cpgs = t(test)
boxplot_dat = merge(pd,cpgs,how='left',by.x=0,by.y=0,all.x=T)
par(mfrow=c(3,2))
boxplot(cg23406740~Progressor,data=boxplot_dat)
boxplot(cg15509024~Progressor,data=boxplot_dat)
boxplot(cg16218368~Progressor,data=boxplot_dat)
boxplot(cg04935278~Progressor,data=boxplot_dat)
boxplot(cg21016177~Progressor,data=boxplot_dat)
boxplot(cg16970087~Progressor,data=boxplot_dat)

#### specify model - Progressive vs. Remitting (NP vs. C) #### 
NP_vs_C = subset(meth_cov_peer, Progressor %in% c("C","NP"))
NP_vs_C$Progressor = factor(NP_vs_C$Progressor,levels=c("C","NP"))
test_dat = dat[ , colnames(dat) %in% rownames(NP_vs_C) ]
design <- model.matrix( ~ Progressor + Corrected_sex + Race + P2, NP_vs_C)

#fit model to data
fit <- lmFit(test_dat, design)
colnames(design)
contrasts <- makeContrasts(diagnosis = "ProgressorNP", levels = design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)
summary(decideTests(fit))

res <- limma::topTable(fit, num = Inf, coef = "diagnosis")
bake <- bacon(res$t)
#plot bacon results. Histogram is of test statistics. Black line represent normal distribution. Red line represents bacon-derived null distribution. Confirm that the red line better fits the test statistic distribution. 
plot(bake)
traces(bake) #parameter estimates should stabilize over iterations; high stochasticity can indicate bacon not performing well. 

#update results with bacon-adjusted test statistics
res$bacon_p <- bacon::pval(bake)
res$bacon_q <- p.adjust(res$bacon_p, method = "BH")
nrow(res[res$bacon_q <= .05, ])

medianbeta_sarc <- rowMedians(betas[,colnames(betas) %in% rownames(NP_vs_C[NP_vs_C$Progressor == "NP", ])])
medianbeta_control <- rowMedians(betas[,colnames(betas) %in% rownames(NP_vs_C[NP_vs_C$Progressor == "C", ])])
medianbeta_dif <- medianbeta_sarc - medianbeta_control
beta_df <- cbind(medianbeta_sarc, medianbeta_control, medianbeta_dif)
rownames(beta_df) <- rownames(betas)
beta_df <- beta_df[rownames(res),]
res <- cbind(res, beta_df)

res <- merge(res,ann, by = "row.names")
rownames(res) <- res$Row.names
res$Row.names <- NULL

res <- res[order(res$bacon_p),]
head(res,7)
write.csv(res, "NPvC_DifferentialMethylation_FullResults_peer_correctedsex.csv")
write.csv(res[res$bacon_p <= .05,], "NPvC_DifferentialMethylation_SignificantResults_peer_correctedsex.csv")

#### specify model - Progressive vs. Remitting (P vs. C) #### 
P_vs_C = subset(meth_cov_peer, Progressor %in% c("C","P"))
P_vs_C$Progressor = factor(P_vs_C$Progressor,levels=c("C","P"))
test_dat = dat[ , colnames(dat) %in% rownames(P_vs_C) ]
design <- model.matrix( ~ Progressor + Corrected_sex + Race + P2, P_vs_C)

#fit model to data
fit <- lmFit(test_dat, design)
colnames(design)
contrasts <- makeContrasts(diagnosis = "ProgressorP", levels = design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)
summary(decideTests(fit))

res <- limma::topTable(fit, num = Inf, coef = "diagnosis")
bake <- bacon(res$t)
#plot bacon results. Histogram is of test statistics. Black line represent normal distribution. Red line represents bacon-derived null distribution. Confirm that the red line better fits the test statistic distribution. 
plot(bake)
traces(bake) #parameter estimates should stabilize over iterations; high stochasticity can indicate bacon not performing well. 

#update results with bacon-adjusted test statistics
res$bacon_p <- bacon::pval(bake)
res$bacon_q <- p.adjust(res$bacon_p, method = "BH")
nrow(res[res$bacon_q <= .05, ])

medianbeta_sarc <- rowMedians(betas[,colnames(betas) %in% rownames(P_vs_C[P_vs_C$Progressor == "P", ])])
medianbeta_control <- rowMedians(betas[,colnames(betas) %in% rownames(P_vs_C[P_vs_C$Progressor == "C", ])])
medianbeta_dif <- medianbeta_sarc - medianbeta_control
beta_df <- cbind(medianbeta_sarc, medianbeta_control, medianbeta_dif)
rownames(beta_df) <- rownames(betas)
beta_df <- beta_df[rownames(res),]
res <- cbind(res, beta_df)

res <- merge(res,ann, by = "row.names")
rownames(res) <- res$Row.names
res$Row.names <- NULL

res <- res[order(res$bacon_p),]
head(res,7)
write.csv(res, "PvC_DifferentialMethylation_FullResults_peer_correctedsex.csv")
write.csv(res[res$bacon_p <= .05,], "PvC_DifferentialMethylation_SignificantResults_peer_correctedsex.csv")

#### Working with only P vs. NP in PEER - 05/22/23 ####
dat=read.csv("Final_Ms_sarc.csv",row.names=1,check.names=FALSE)
meth_cov_peer <- read.csv("meth_cov_peer.csv",stringsAsFactors = TRUE,row.names=1)
ann <- as.data.frame(sesameAnno_get("Anno/EPIC/EPIC.hg38.manifest.rds"))
betas = readRDS("betas_cleaned.rds")

pd <- read.csv("../meth_phenotype.csv")
Ms = readRDS("Ms.rds")
Ms_sub <- Ms[!rownames(Ms) %in% ann[ann$CpG_chrm == "chrX",]$Probe_ID & !rownames(Ms) %in% ann[ann$CpG_chrm == "chrX",]$Probe_ID, ]
pca <- prcomp(t(Ms_sub))
summary(pca)
pd <- merge(pd, pca$x, by.x = "sentrix", by.y = "row.names")
rownames(pd) <- pd$sentrix

P_vs_NP = subset(meth_cov_peer, Progressor %in% c("NP","P"))
P_vs_NP$Progressor = factor(P_vs_NP$Progressor,levels=c("NP","P"))
test_dat = dat[ , colnames(dat) %in% rownames(P_vs_NP) ]
design <- model.matrix( ~ Progressor + Corrected_sex + Race + P1, P_vs_NP)

write.csv(P_vs_NP, "PvNP_metadata.csv")

#fit model to data
fit <- lmFit(test_dat, design)
colnames(design)
contrasts <- makeContrasts(diagnosis = "ProgressorP", levels = design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)
summary(decideTests(fit))

res <- limma::topTable(fit, num = Inf, coef = "diagnosis")
bake <- bacon(res$t)
#plot bacon results. Histogram is of test statistics. Black line represent normal distribution. Red line represents bacon-derived null distribution. Confirm that the red line better fits the test statistic distribution. 
plot(bake)
traces(bake) #parameter estimates should stabilize over iterations; high stochasticity can indicate bacon not performing well. 

#update results with bacon-adjusted test statistics
res$bacon_p <- bacon::pval(bake)
res$bacon_q <- p.adjust(res$bacon_p, method = "BH")
nrow(res[res$bacon_q <= .05, ])

medianbeta_sarc <- rowMedians(betas[,colnames(betas) %in% rownames(P_vs_NP[P_vs_NP$Progressor == "P", ])])
medianbeta_control <- rowMedians(betas[,colnames(betas) %in% rownames(P_vs_NP[P_vs_NP$Progressor == "NP", ])])
medianbeta_dif <- medianbeta_sarc - medianbeta_control
beta_df <- cbind(medianbeta_sarc, medianbeta_control, medianbeta_dif)
rownames(beta_df) <- rownames(betas)
beta_df <- beta_df[rownames(res),]
res <- cbind(res, beta_df)

res <- merge(res,ann, by = "row.names")
rownames(res) <- res$Row.names
res$Row.names <- NULL

res <- res[order(res$bacon_p),]
head(res,7)

par(mfrow=c(1,1))
EnhancedVolcano(res, lab = res$gene, x = "logFC", y = "bacon_q", pCutoff = 0.05) 
res$chrom <- gsub("chr", "", res$seqnames)
res[res$chrom == "X", ]$chrom <- "23"
res[res$chrom == "Y", ]$chrom <- "24"
res$chrom <- as.numeric(res$chrom)
manhattan(res, chr = "chrom", bp = "start", p = "bacon_q", snp = "gene", suggestiveline = -log10(.05), genomewideline = 0, annotatePval = 0.05)

write.csv(res, "PvNP_DifferentialMethylation_FullResults_peer_correctedsex.csv")
write.csv(res[res$bacon_p <= .05,], "PvNP_DifferentialMethylation_SignificantResults_peer_correctedsex.csv")

#### Male vs. Female ####
peer_factors = read.csv(file="peer_20k_factors.csv",head=TRUE,sep=",",check.names=FALSE,stringsAsFactors = FALSE)
actual_pf = peer_factors[,c(6:25)]

# PCRA
#peer_factors_mod = peer_factors[,c(2:ncol(peer_factors))]
covs=read.csv("meth_cov.csv",stringsAsFactors=F)
test_pheno = covs[,c(1:8)]
test_pheno$Sentrix.Barcode = as.factor(test_pheno$Sentrix.Barcode)

rownames(actual_pf) = rownames(test_pheno)

#### trying with peer factors
dat=read.csv("Final_Ms_sarc.csv",row.names=1,check.names=FALSE)
#dat = Ms_filtered
meth_cov_peer <- read.csv("meth_cov_peer.csv",stringsAsFactors = TRUE,row.names=1)
ann <- as.data.frame(sesameAnno_get("Anno/EPIC/EPIC.hg38.manifest.rds"))
betas = readRDS("betas_cleaned.rds")

pd <- read.csv("../meth_phenotype.csv")
Ms = readRDS("Ms.rds")
Ms_sub <- Ms[!rownames(Ms) %in% ann[ann$CpG_chrm == "chrX",]$Probe_ID & !rownames(Ms) %in% ann[ann$CpG_chrm == "chrX",]$Probe_ID, ]
pca <- prcomp(t(Ms_sub))
summary(pca)
pd <- merge(pd, pca$x, by.x = "sentrix", by.y = "row.names")
rownames(pd) <- pd$sentrix

##specify model - Male  vs. Female
design <- model.matrix( ~ Corrected_sex + P2, meth_cov_peer)

#fit model to data
fit <- lmFit(dat, design)
colnames(design)
contrasts <- makeContrasts(Corrected_sex = 1, levels = design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)
summary(decideTests(fit))

res <- limma::topTable(fit, num = Inf, coef = "Corrected_sex")
bake <- bacon(res$t)
#plot bacon results. Histogram is of test statistics. Black line represent normal distribution. Red line represents bacon-derived null distribution. Confirm that the red line better fits the test statistic distribution. 
plot(bake)
traces(bake) #parameter estimates should stabilize over iterations; high stochasticity can indicate bacon not performing well. 

#update results with bacon-adjusted test statistics
res$bacon_p <- bacon::pval(bake)
res$bacon_q <- p.adjust(res$bacon_p, method = "BH")
nrow(res[res$bacon_q <= .05, ])

medianbeta_sarc <- rowMedians(betas[,colnames(betas) %in% rownames(pd[pd$DX == "Sarc", ])])
medianbeta_control <- rowMedians(betas[,colnames(betas) %in% rownames(pd[pd$DX == "Control", ])])
medianbeta_dif <- medianbeta_sarc - medianbeta_control
beta_df <- cbind(medianbeta_sarc, medianbeta_control, medianbeta_dif)
rownames(beta_df) <- rownames(betas)
beta_df <- beta_df[rownames(res),]
res <- cbind(res, beta_df)

res <- merge(res,ann, by = "row.names")
rownames(res) <- res$Row.names
res$Row.names <- NULL

res <- res[order(res$bacon_p),]
head(res,7)
write.csv(res[res$bacon_p <= .05,], "MaleVsFemale_DifferentialMethylation_SignificantResults_peer_correctedsex.csv")



#### functions ####
# Try making a PCRA using non-DESeq2 approach
lmmatrix <- function(y,cov)
{
  p <- matrix(NA,ncol=ncol(cov),nrow=ncol(y))
  for(j in 1:ncol(cov))
  {
    x <- cov[,j]
    for(i in 1:ncol(y))
    {
      fit <- summary(lm(y[,i]~x,na.action=na.omit))
      f <- fit$fstatistic
      p[i,j] <- pf(f["value"],f["numdf"],f["dendf"],lower.tail=FALSE)
    }
  }
  colnames(p) <- names(cov)
  p
}

plotp<-function(jpgfile,p,yaxis,xmax,title)
{
  jpeg(filename=jpgfile,width=1000,height=700,quality = 100)
  par(mar=c(5, 4.2, 4, 2) + 0.1)
  plot(1,xlim=c(0,xmax),ylim=c(0,length(yaxis)+1),type="n",bty="n",
       axes=FALSE,xlab="Principal Component",ylab="",main=title)
  axis(1,at=c(1:xmax),pos=0.5,las=1,lwd=3)
  for(i in 1:length(yaxis)){text(0.3,i, yaxis[i],xpd=TRUE,adj=1)}
  for(i in 1:ncol(p)){
    for(j in 1:nrow(p)){
      pp <- p[j,i]
      colcode <- "white"
      if(pp <= 10e-10){colcode="darkred"}
      else if (pp <= 10e-5){colcode="red"}
      else if (pp <= 0.01){colcode="orange"}
      else if(pp <= 0.05){colcode="pink"}
      polygon(c(j-0.5,j-0.5,j+0.5,j+0.5), c(i-0.5,i+0.5,i+0.5,i-0.5),col=colcode,
              border=NA)
    }}
  legend("topright",c("<0.05","<0.01","<10E-5","<10E-10"),
         col=c("pink","orange","red","darkred"),
         pch=15,pt.cex=2,bty="o",horiz=TRUE,xpd=TRUE)
  dev.off()
}

 
