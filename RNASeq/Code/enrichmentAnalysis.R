##########################################
# Gene set enrichment analysis
##########################################

doEnrichment <- function(dat){
  setEnrichrSite("Enrichr") # Human genes
  
  dbs <- c("BioPlanet_2019", "Reactome_2016", "GO_Biological_Process_2021", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018", 
           "ARCHS4_Tissues", "KEGG_2019_Human", "MSigDB_Hallmark_2020", "WikiPathways_2019_Human")
  
  enrichement_list_up <- list()
  enrichement_list_down <- list()
  
  #nms <- dat
  degs <- dat
  
 # for (i in 1:length(nms)){
    #degs <- read.xlsx("H:/Output files/Bulk RNASeq outputs/P_v_NP_v_C_outliers_removed_Maier_BlkRNASeq.xlsx",
    #                  sheet = i)
    geneset_up <- as.character(degs[is.na(degs$padj)==F & degs$log2FoldChange>0, 'gene'])[1:100]
    geneset_down <- as.character(degs[is.na(degs$padj)==F & degs$log2FoldChange<0, 'gene'])[1:100]
        
    #enriched <- enrichr(geneset, dbs)
    enriched_up <- enrichr(geneset_up, dbs)
    enriched_down <- enrichr(geneset_down, dbs)
        
    enriched_up <- enriched_up[lapply(enriched_up, function(x) nrow(x)) > 0]
    enriched_down <- enriched_down[lapply(enriched_down, function(x) nrow(x)) > 0]
    
    for (j in 1:length(enriched_up)){
      #enriched_up[[j]]$comparison <- nms[i]
      enriched_up[[j]]$database <- dbs[j]
      enriched_up[[j]] <- enriched_up[[j]][enriched_up[[j]]$Adjusted.P.value < 0.05,]
    }
    
    for (j in 1:length(enriched_down)){
      #enriched_down[[j]]$comparison <- nms[i]
      enriched_down[[j]]$database <- dbs[j]
      enriched_down[[j]] <- enriched_down[[j]][enriched_down[[j]]$Adjusted.P.value < 0.05,]
    }
    
    enrichement_list_up <- do.call(rbind, enriched_up)
    enrichement_list_down <- do.call(rbind, enriched_down)
            
 # }
  
#  names(enrichement_list_up) <- paste0(nms, '_up')
#  names(enrichement_list_down) <-paste0(nms, '_down')
  
 out <- list(geneset_down, geneset_up, enrichement_list_down, enrichement_list_up) 
 return(out)
}
