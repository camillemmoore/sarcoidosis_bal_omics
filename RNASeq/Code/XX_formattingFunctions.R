
# get the numbers ready for the table
formatP <- function(x){return(ifelse(x < 0.0001, formatC(x, format = "E", digits = 1), round(x, 4)))}

# get the data formatted for output
formatDat <- function(dat){
  dat <- data.frame(dat)
  dat <- dat %>% relocate(gene)
  dat$foldChange <- 2^(dat$log2FoldChange)
  non_pval_col <- colnames(dat)[!colnames(dat) %in% c("gene", "pvalue", "padj")]
  dat[, non_pval_col] <- apply(dat[, non_pval_col], 2, function(x) round(x, 4))
  pval_col <- c("pvalue", "padj")
  dat[, pval_col] <- apply(dat[, pval_col], 2, function(x) round(x, 8))
  return(dat)
}

# put together the results table
getDEResultsTab <- function(dat){
  dat <- dat %>% select(gene, baseMean, log2FoldChange, stat, pvalue, padj) %>% mutate_if(is.numeric, formatP)
  return(datatable(dat, extensions = c('Buttons', 'KeyTable', "FixedHeader"),
                   options = list(dom = 'lBftrip',
                                  buttons = c('copy', 'csv', 'excel', 'print'),
                                  keys = TRUE,
                                  searchHighlight = TRUE,
                                  pageLength = 10,
                                  lengthMenu = c("10", "25", "50", "100"),
                                  scrollX  = TRUE), rownames = FALSE))
}

# plot volcano 
plotVolcano <- function(dat, cutoff = 0.05){
  return(dat %>% mutate(threshold = ifelse(padj < cutoff, "FDR < 0.05", "FDR > 0.05")) %>% 
           filter(is.na(threshold) == F) %>%
           ggplot(aes(x = log2FoldChange, y = -log10(pvalue), color = threshold, label = gene)) + 
           geom_point() +
           labs(x = "log2 fold change", y = "-log10 p-value", color = "") + 
           theme(legend.position = "none", legend.title = element_blank(), 
                 plot.title = element_text(size = rel(1.5), hjust = 0.5), 
                 axis.title = element_text(size = rel(1.25))) +
           theme_bw()
  )
}

# get together enrichment tab 
getEnResutltsTab <- function(dat){
  return(datatable(dat %>% #arrange(database, P.value) %>% 
              select(database, Term, Overlap, P.value, Adjusted.P.value, Odds.Ratio, Genes) %>% 
              mutate_if(is.numeric, formatP),
            extensions = c('Buttons', 'KeyTable', "FixedHeader"),
            options = list(dom = 'lBftrip',
                           buttons = c('copy', 'csv', 'excel', 'print'),
                           keys = TRUE,
                           searchHighlight = TRUE,
                           pageLength = 10,
                           lengthMenu = c("10", "25", "50", "100"),
                           scrollX = TRUE),
            rownames = FALSE)
  )
}