# MULTIOMICS ANALYSIS WITH DIABLO


library(openxlsx)
library(tidyverse)
library(DESeq2)
library(edgeR)
library(mixOmics)

dataDir <- './BAL_TCELL_RNASEQ/'
setwd(dataDir)

##############################################################
# RNASEQ DATA MANAGEMENT
##############################################################
# rnaseq residuals (sex regressed out)
result.counts <- read.xlsx(paste0(dataDir, "Sex_regressed_residual_counts.xlsx"), rowNames = T)

# rnaseq metadata
all_md <- read.xlsx("rnaseq_md.xlsx")
rownames(all_md) <- all_md$Sample.ID

# gene names
genes_tab <- read.xlsx("gene_names.xlsx", rowNames = T)
genes <- genes_tab$genes
names(genes) <- rownames(genes_tab)
table(rownames(all_md)==rownames(result.counts))

#### Select top n variable genes to include in DIABLO
n <- 15000
bulk.res.sort.var <- colVars(as.matrix(result.counts))
names(bulk.res.sort.var) <- colnames(result.counts)
bulk.res.sort.var <- bulk.res.sort.var[order(bulk.res.sort.var, decreasing = TRUE)]
top.n.genes <- names(bulk.res.sort.var)[1:n]


##############################################################
# METHYLATION DATA
##############################################################

# methylation residuals
meth.counts <- read.csv("methylation_residuals_063023.csv", row.names = 1)

# methylation metadata
meth.md <- read.csv("meth_phenotype.csv")
colnames(meth.md)[1] <- "Sample.ID"
meth.md$sentrix <- paste0("X", meth.md$sentrix)

# only keeping subjects with methylation data
meth.md <- meth.md[meth.md$sentrix %in% colnames(meth.counts), ]

# matching order of columns in counts to that with md to get the sample id with same order
meth.md <- meth.md[match(colnames(meth.counts), meth.md$sentrix), ]

# replacing the sentrix names with sample id in the count data
colnames(meth.counts) <- meth.md$Sample.ID

# replacing the rownames of md with sample id
rownames(meth.md) <- meth.md$Sample.ID

# Top n variable cpgs for DIABLO
n <- 15000
meth.res.sort.var <- meth.counts
meth.res.sort.var$vars <- rowVars(as.matrix(meth.res.sort.var))
meth.res.sort.var <- meth.res.sort.var[order(meth.res.sort.var$vars, decreasing = TRUE), ]
top.n.cpg <- rownames(meth.res.sort.var)[1:n]


##############################################################
# Get data for multi-omics networks
##############################################################
both.md <- meth.md[rownames(meth.md) %in% rownames(result.counts), ]
meth.res <- meth.counts[, rownames(both.md)]
bulk.res <- result.counts[rownames(both.md),]
bulk.res <- bulk.res[,colnames(bulk.res) %in% top.n.genes ]
meth.res <- meth.res[rownames(meth.res) %in% top.n.cpg,]
table(rownames(both.md)==rownames(bulk.res))
table(rownames(both.md)==colnames(meth.res))

write.xlsx(list(md = both.md, rnaseq = bulk.res, meth =meth.res),
           'data_for_multiomics.xlsx', rowNames=T)

##############################################################
# DIABLO
##############################################################
both.md <- read.xlsx('data_for_multiomics.xlsx', sheet = 'md', rowNames = T)
meth.res <- read.xlsx('data_for_multiomics.xlsx', sheet = 'meth', rowNames = T)
bulk.res <- read.xlsx('data_for_multiomics.xlsx', sheet = 'rnaseq', rowNames = T)

bulk.top.genes.res <- t(bulk.res[,colnames(bulk.res) %in% bulk.top.genes ])
meth.top.cpgs.res <- meth.res[rownames(meth.res) %in% meth.top.cpgs,]
all.dx <- factor(both.md$DX, levels = c("Control", "Sarc"))
names(all.dx) <- rownames(both.md)

# RUN DIABLO MODELS FOR 3 DIFFERENT DESIGN VALUES
for(i in c(0.25, 0.5, 0.75)){

design_value <- i

## wrangling the data in the required format
X <- list(bulk = t(bulk.top.genes.res),
          meth = t(meth.top.cpgs.res))

Y <- all.dx

## Design matrix
design <- matrix(design_value, ncol = length(X), nrow = length(X), 
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

## fitting the PLSDA model
diablo.mod <- block.plsda(X, Y, ncomp = 5, design = design)

## obtain the performance of the model with respect to the different prediction distances
set.seed(1) # For reproducibility
perf.diablo.mod <- perf(diablo.mod, validation = 'Mfold', folds = 5, nrepeat = 10)

## choosing number of components
ncomp <- perf.diablo.mod$choice.ncomp$WeightedVote[2,3]

## choosing number of variables
# First round of tuning: grid search from 5 to 55 by 5
test.keepX <- list(bulk = seq(5, 55, 5),
                   meth = seq(5, 55, 5))
tune.diablo.mod <- tune.block.splsda(X, Y, ncomp = ncomp, 
                                     test.keepX = test.keepX, design = design,
                                     validation = 'Mfold', folds = 5, nrepeat = 10, 
                                     BPPARAM = BiocParallel::SnowParam(),
                                     dist = "mahalanobis.dist")

list.keepX <- tune.diablo.mod$choice.keepX

# Second round of tuning using a more fine grid
temp <- unique(list.keepX$bulk)
temp2 <- unique(list.keepX$meth)
bulk.temp <- meth.temp <- NULL
for(j in 1:length(temp)){bulk.temp <- c(bulk.temp, max(temp[j]-5, 5):(temp[j]+5))
meth.temp <- c(meth.temp, max(temp2[j]-5, 5):(temp2[j]+5))}
bulk.temp <- unique(bulk.temp)
meth.temp <- unique(meth.temp)
  
test.keepX <- list(bulk = c(bulk.temp),
                     meth = c(meth.temp))

                   
                   
tune.diablo.mod <- tune.block.splsda(X, Y, ncomp = ncomp, 
                                     test.keepX = test.keepX, design = design,
                                     validation = 'Mfold', folds = 5, nrepeat = 10, 
                                     BPPARAM = BiocParallel::SnowParam(),
                                     dist = "mahalanobis.dist")

list.keepX <- tune.diablo.mod$choice.keepX


## final model
diablo.mod.final <- block.splsda(X, Y, ncomp = ncomp, 
                                 keepX = list.keepX, design = design)


Y_df <- as.data.frame(diablo.mod.final$Y)
Y_merged <- merge(Y_df, both.md, by = "row.names", sort = FALSE)[, c("Sample.ID", "Progressor")]

# Create visualizations
pdf(paste0(dataDir, 'diablo_plots_set_', "15kvariable_design_",design_value, '.pdf'), width=11, height=8.5)

for(comp in 1:ncomp){
plotDiablo(diablo.mod.final, ncomp = comp)
  
df_1 <- data.frame(x = diablo.mod.final$variates$meth[,comp], 
                   y = diablo.mod.final$variates$bulk[,comp],
                   group = diablo.mod.final$Y, 
                   id = names(diablo.mod.final$variates$meth[,comp]))

df_1 <- merge(df_1, Y_merged, by.x = "id", by.y = "Sample.ID", sort = FALSE)

p_1 <- ggplot(df_1, aes(x = x, y = y, color = Progressor, label = id)) +
  geom_point() + ggtitle(paste0('Component ', comp))
  theme_classic()
print(p_1)

cimDiablo(diablo.mod.final, color.blocks = c('darkorchid', 'lightgreen'),
          comp = comp, margin=c(20,20), legend.position = "right")

}
dev.off()

saveRDS(diablo.mod.final, paste0(dataDir, 'diablo_model_set_15k_variable_design', design_value, '.rds'))


# Save components with annotation
epic <- readRDS('./EPIC.hg38.manifest.rds')
epic <- as.data.frame(epic)

component_elements <- list()
nms <- paste0(c(rep("bulk", ncomp), rep('meth', ncomp)), rep(1:ncomp, 2))

for(comp in 1:ncomp){
  temp_nms <- nms[grep(comp, nms)]
component_elements[[temp_nms[1]]] <- genes[rownames(selectVar(diablo.mod.final, block = 'bulk', comp = comp)$bulk$value)]
component_elements[[temp_nms[2]]] <- as.data.frame(epic[rownames(selectVar(diablo.mod.final, block = 'meth', comp = comp)$meth$value),])
                           
}

write.xlsx(component_elements, paste0(dataDir, 'diablo_components_set_', '15kvariable_design_', design_value,'.xlsx'), rowNames=T)

}



############### ADDITIONAL VISUALIZATIONS
library(pheatmap)
library(biomaRt)

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# HEATMAP WITH ANNOTATIONS
for(design_value in c(0.25, 0.5, 0.75)){
  sns <- getSheetNames(paste0(dataDir, 'diablo_components_set_', '15kvariable_design_', design_value,'.xlsx'))
  ncomp <- length(sns)/2
  pdf(paste0(dataDir, 'diablo_heatmaps_', '15kvariable_design_', design_value,'.pdf'), height=11, width=8.5)
  
  for(j in 1:ncomp){
    rna_features <- read.xlsx(paste0(dataDir, 'diablo_components_set_', '15kvariable_design_', design_value,'.xlsx'), sheet=paste0('bulk',j), colNames = F)
    ensemblids <- names(genes[genes %in% rna_features$X1])
    ensemblids <- unlist(lapply(strsplit(ensemblids,split='[.]'), function(x) x[1]))
    rna_annot <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                                      'start_position', 'end_position'),
                       filters = 'ensembl_gene_id', 
                       values = ensemblids, 
                       mart = ensembl)
    
    rna_annot <- merge(rna_annot, data.frame(ensembl_gene_id = ensemblids, 
                                             orig_id = names(genes[genes %in% rna_features$X1])), all.y=T)
    
    rownames(rna_annot) <- rna_annot$orig_id
    
    rna_annot$chromosome_name <- as.character(rna_annot$chromosome_name )
    
    meth_features <- read.xlsx(paste0(dataDir, 'diablo_components_set_', '15kvariable_design_', design_value,'.xlsx'),sheet=paste0('meth',j))
    colnames(meth_features)[1] <- 'id'
    meth_annot <- meth_features[,c("id", "seqnames", "start", "end", 'gene_HGNC')]
    
    meth_annot$hgnc_symbol <- meth_annot$gene_HGNC
    meth_annot$start_position <- meth_annot$start
    meth_annot$end_position <- meth_annot$end
    meth_annot$chromosome_name <- (gsub('chr', '',meth_annot$seqnames))
    meth_annot$orig_id <- meth_annot$id
    meth_annot$Type <- 'Methylation'
    
    rna_annot$id <- rna_annot$ensembl_gene_id
    rna_annot$Type <- "Expression"
    
    feature_annot <- rbind(meth_annot[,c('hgnc_symbol', 'chromosome_name', 'orig_id', 'id', 'Type')], 
                           rna_annot[,c('hgnc_symbol', 'chromosome_name', 'orig_id', 'id', 'Type')])
    
    rownames(feature_annot) <- feature_annot$orig_id
    feature_annot$display <- paste0(feature_annot$id,ifelse(is.na(feature_annot$hgnc_symbol), '',
                                                            ' - ') , 
                                    ifelse(is.na(feature_annot$hgnc_symbol), '',
                                           feature_annot$hgnc_symbol))
    
    feature_annot$Chr <- factor(feature_annot$chromosome_name)
    
    
    feature_mat <- rbind(bulk.top.genes.res[rownames(feature_annot[feature_annot$Type=='Expression',]),], 
                         meth.top.cpgs.res[rownames(feature_annot[feature_annot$Type=='Methylation',]),])
    
    annot_col <- as.data.frame(both.md[,c('Progressor')])
    rownames(annot_col) <- rownames(both.md)
    colnames(annot_col) <- "Diagnosis"
    
    
    annot_colors <- list(DX=c('Control' = '#a6cee3', 'Sarc'='#fb6a4a'),
                         Diagnosis = c('C' = '#a6cee3', 
                                       'NP' = '#fcbba1', 
                                       'P'='#e31a1c'), 
                         Type = c('Expression' = '#6a3d9a', 'Methylation' = '#cab2d6'))
    
    
    
    # Nice color scale from RColorBrewer
    mycols = colorRampPalette(c("blue","white","red"))(1000)
    breakscale <- c(-4,seq(-3,3, length.out=length(mycols)-1),5)
    
    display_nms <- feature_annot$display
    names(display_nms) <- rownames(feature_annot)
    display_nms <- display_nms[rownames(feature_mat)]
    
    pheatmap(feature_mat, scale='row',
             annotation_col = annot_col, 
             annotation_colors = annot_colors, 
             annotation_row = feature_annot[,c('Type', 'Chr')],
             color=mycols, fontsize_row = 6,
             breaks=breakscale,
             main = paste0('Component ', j), 
             labels_row = display_nms, 
             clustering_method = "ward.D2")
    
    
  }
  
  dev.off()
}

# SCATTER PLOTS OF EXPRESSION AND METHYLATION COMPONENTS
for(design_value in c(0.25, 0.5, 0.75)){
diablo.mod.final <- readRDS(paste0(dataDir, 'diablo_model_set_15k_variable_design', design_value, '.rds'))

Y_df <- as.data.frame(diablo.mod.final$Y)
Y_merged <- merge(Y_df, both.md, by = "row.names", sort = FALSE)[, c("Sample.ID", "Progressor")]

colnames(Y_merged) <- c('Sample.ID', 'Diagnosis')

pdf(paste0(dataDir, 'diablo_scatter_', "15kvariable_design_",design_value, '.pdf'), width=5, height=4)

for(comp in 1:diablo.mod.final$ncomp[1]){
  #plotDiablo(diablo.mod.final, ncomp = comp)
  
  df_1 <- data.frame(x = diablo.mod.final$variates$meth[,comp], 
                     y = diablo.mod.final$variates$bulk[,comp],
                     group = diablo.mod.final$Y, 
                     id = names(diablo.mod.final$variates$meth[,comp]))
  
  df_1 <- merge(df_1, Y_merged, by.x = "id", by.y = "Sample.ID", sort = FALSE)
  
  p_1 <- ggplot(df_1, aes(x = x, y = y, color = Diagnosis, label = id)) +
    geom_point() + ggtitle(paste0('Component ', comp))+annotate("text",x=max(df_1$x),y=min(df_1$y),hjust=1,label=paste0('r = ', round(cor(df_1$x, df_1$y),2)))+
    scale_color_manual(values = c('C' = '#a6cee3', 
                                  'NP' = '#fcbba1', 
                                  'P'='#e31a1c'))+ylab(paste0('Expression Component ', comp))+xlab(paste0('Methylation Component ', comp))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5))
  print(p_1)
}
dev.off()
}
