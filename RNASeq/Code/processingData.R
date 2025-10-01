
getDataStarted <- function(dataDir, nSmallGroup, outliers = F){
  pheno.data <- read.xlsx(paste0(dataDir, "sarc meth data with phenotype.xlsx"), detectDates = T)
  pheno.data.fltrd <- pheno.data[pheno.data$Phenotype %in% c("P", "NP", "control"), ]
  pheno.data.fltrd$Date.of.Birth <- as.Date(pheno.data.fltrd$Date.of.Birth)
  pheno.data.fltrd$date.of.sample <- as.Date(pheno.data.fltrd$date.of.sample)
  pheno.data.fltrd$Age.at.sampling <- trunc((pheno.data.fltrd$Date.of.Birth %--%
                                               pheno.data.fltrd$date.of.sample) / years(1))
  raw.counts <- read.xlsx(paste0(dataDir, "bal_cd4+_rnaseqc_rnaseq_reads.xlsx"))
  all_md <- pheno.data.fltrd
  
  if(outliers == T){# Removing outliers
  outliers <- c("39810", "33042", "16092")
  all_md <- all_md[!c(all_md$Sample.ID %in% outliers), ]
  }
  
  counts <- raw.counts
  sample.id.raw <- colnames(counts)[-c(1:2)]
  sample.id <- gsub(".{1}$", "", (sapply(strsplit(sample.id.raw, "_"), "[[", 2)))
  all_md <- all_md[all_md$Sample.ID %in% sample.id, ]
  all_md <- all_md %>% distinct(Sample.ID, .keep_all = TRUE)
  colnames(counts)[-c(1:2)] <- sample.id
  rownames(all_md) <- all_md$Sample.ID
  all_md$Sample.ID <- as.character(all_md$Sample.ID)
  
  # set the rownmaes of the count matrix to the ensembl gene ID
  rownames(counts) <- counts$Name
  
  # save the gene names in a vector to help label results later
  genes <- counts$Description
  names(genes) <- counts$Name
  
  # Subset the count matrix to columns for samples in the meta data file
  counts <- counts[,colnames(counts) %in% rownames(all_md)]
  
  # Similarly, subset MD to only those samples in the RNAseq data
  all_md <- all_md[colnames(counts),]
  
  # Filter out lowly expressed genes####
  dgList <- DGEList(counts = counts)
  
  countsPerMillion <- cpm(dgList)
  countCheck <- countsPerMillion > 1
  keep <- which(rowSums(countCheck) >= nSmallGroup) # I usually include genes with at least 1 CPM in N samples where N is the sample size of the smallest group
  dgList <- dgList[keep,]
  
  counts.filtered <- counts[names(keep),]
  
  # gene or feature data
  fdata <- data.frame(geneid = names(genes) , gene_name = genes)
  rownames(fdata) <- names(genes)
  fdata <- fdata[rownames(counts.filtered),]
  
  
  out <- list(counts.filtered, all_md, fdata, genes)
  return(out)
}
