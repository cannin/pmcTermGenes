library(org.Hs.eg.db)
library(magrittr)

#' Retrieve Pubtator Data
#' 
#' @param term the search term; there needs to be a corresponding TERM_genes.rds created by getPubtator() (e.g., sarcoma)
#' @param minGenes keep result entries with this minimum count (Default: 0)
#' @param pmcidsToKeep vector of PMC IDs to keep and discard others; ignored if NULL (Default: NULL)
#' @param verbose show debugging output (Default: TRUE)
#' 
#' @return a data.frame with gene and pmcid_cnt columns
#' 
#' @note Current code only matches human genes. Publications with genes exclusively from other organisms will return errors
#' 
#' @examples 
#' #getPubtator(term="sarcoma", pmcids, sleepTime=0, pmcidsFilter="intersect", forceFlag=FALSE, randomizePmcids=FALSE)
#'
#' @export
getPmcGeneCounts <- function(term, minGenes=0, pmcidsToKeep=NULL, verbose=FALSE) {
  pmcGenes <- readRDS(paste0(term, "_genes.rds"))
  
  allGenes <- list()
  
  cat("TOTAL: ", length(pmcGenes), "\n")    
  
  for(i in 1:length(pmcGenes)) {
    # i <- 1
    if(verbose) {
      cat("I: ", i, "\n")      
    }

    tryCatch({
      pmcid <- names(pmcGenes)[i]
      
      if(!is.null(pmcidsToKeep) && !(pmcid %in% pmcidsToKeep)) {
        next
      }
      
      entrezids <- pmcGenes[[pmcid]]
      symbols <- suppressMessages(mapIds(org.Hs.eg.db, keys=entrezids, column="SYMBOL", keytype="ENTREZID", multiVals="first"))
      allGenes[[pmcid]] <- symbols %>% unname %>% sort
      
      if(verbose) {
        cat("I: ", i, " SUCCESS: PMCID: ", pmcid, "\n")
      }
    }, error=function(e) {
      cat("ERROR: PMCID: ", pmcid, "\n")
    })
  }
  
  tmp <- allGenes %>% unlist %>% table %>% sort(., decreasing = TRUE)
  termGenes <- tmp[which(tmp >= minGenes)]

  termGenesDf <- as.data.frame(termGenes)
  colnames(termGenesDf) <- c("gene", "pmcid_cnt")
  write.table(termGenesDf, paste0(term, "_genes_symbols_cnt.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE)

  #allGenesInverted <- inverseListOfVectors(allGenes)
  #saveRDS(allGenesInverted, paste0(term, "_pmcids_by_gene.rds"))
  
  return(termGenesDf)
}


