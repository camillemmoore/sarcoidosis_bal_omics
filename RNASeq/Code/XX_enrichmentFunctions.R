
# function to run analysis
runEnrichment <- function(dat){
  geneset_up   <- as.character(dat[is.na(dat$padj)==F & dat$log2FoldChange>0, 'gene'])[1:100]
  geneset_down <- as.character(dat[is.na(dat$padj)==F & dat$log2FoldChange<0, 'gene'])[1:100]
  
  enriched_up   <- enrichr(geneset_up, dbs)
  enriched_down <- enrichr(geneset_down, dbs)
  
  enriched_up   <- enriched_up[lapply(enriched_up, function(x) nrow(x)) > 0] 
  enriched_down <- enriched_down[lapply(enriched_down, function(x) nrow(x)) > 0] 
  
  for(j in 1:length(enriched_up)){
    enriched_up[[j]]$database <- dbs[j]
    enriched_up[[j]] <- enriched_up[[j]][enriched_up[[j]]$Adjusted.P.value < 0.05,]
    enriched_up[[j]] <- enriched_up[[j]] %>% arrange(P.value)
  }
  
  for(j in 1:length(enriched_down)){
    enriched_down[[j]]$database <- dbs[j]
    enriched_down[[j]] <- enriched_down[[j]][enriched_down[[j]]$Adjusted.P.value < 0.05,]
    enriched_down[[j]] <- enriched_down[[j]] %>% arrange(P.value)
  }
  
  enrichement_list_up <- do.call(rbind, enriched_up)
  enrichement_list_down <- do.call(rbind, enriched_down)
  
  out <- list(enrichement_list_up,  enrichement_list_down)
  
  return(out)
}








