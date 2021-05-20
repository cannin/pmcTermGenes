library(xml2)
library(httr)
library(magrittr)
library(jsonlite)
library(simpleRCache)
library(jqr)

#' Retrieve Pubtator Data
#' 
#' @param term the search term (e.g., sarcoma)
#' @param cacheDir path to a cache directory to save intermediate files so they do not need to be re-downloaded
#' @param pmcids a vector of pmcids
#' @param sleepTime seconds to sleep between PMC ID processing; this to deal 
#'   with rate limits by PubTator (if they exist) (Default: 5) 
#' @param pmcidsFilter two options: "setdiff" download/process articles missing 
#'   from cacheDir OR "intersect" to only process articles already downloaded
#' @param forceFlag creates an infinite loop to force code to ignore all errors (Default: FALSE)
#' @param randomizePmcids randomize the list of PMC IDs so if there connection 
#'   issues, new articles will be run each time (Default: TRUE) 
#' @param mainSectionGenesOnly Extract genes from anywhere in Pubtator results (FALSE) OR 
#'   exclude genes extracted from References, Supplemental, Methods (TRUE) (Default: TRUE) 
#' @param abstractOnly Extract genes only from abstract; requires mainSectionGenesOnly=TRUE (Default: FALSE)
#' @param verbose show debugging output (Default: TRUE)
#' 
#' @note Not all PMC ID files may appear in cacheDir; this may be due to API error or empty file from API
#' 
#' @return a list with PMC IDs as names and genes as the entries 
#' 
#' @examples 
#' getPubtator(term="sarcoma", pmcids, sleepTime=0, pmcidsFilter="intersect", forceFlag=FALSE, randomizePmcids=FALSE)
#'
#' @export
getPubtator <- function(term, 
                        cacheDir, 
                        pmcids,
                        sleepTime=5, 
                        pmcidsFilter=c("setdiff", "intersect"), 
                        forceFlag=FALSE, 
                        randomizePmcids=FALSE, 
                        mainSectionGenesOnly=TRUE, 
                        abstractOnly=FALSE, 
                        verbose=TRUE) {
  # SETUP ----
  # From: https://github.com/jeroen/curl/issues/156 for 'Error in the HTTP2 framing layer' (or install curl with homebrew and setup access in .Renviron)
  httr::set_config(httr::config(http_version=0))
  
  # SETUP CACHING ----
  Sys.setenv("PREFIX_SIMPLERCACHE"="pubtator")
  Sys.setenv("DEBUG_SIMPLERCACHE"="TRUE")
  
  setCacheRootPath(cacheDir)
  GETCached <- addMemoization(GET)
  
  cacheFiles <- dir(cacheDir, "PMC" )
  cachePmcids <- gsub("_.*", "", cacheFiles)
  
  # PUBTATOR ----
  #Sys.setenv("DEBUG_COUNT"="TRUE")
  # while(Sys.getenv("DEBUG_COUNT") == "TRUE") {
  #   cacheFiles <- dir(cacheDir, "PMC" )
  #   cachePmcids <- gsub("_.*", "", cacheFiles)
  #   
  #   tmp <- pmcids %in% cachePmcids %>% which %>% length
  #   cat("T: ", tmp, "\n")
  #   
  #   Sys.sleep(5)
  # }
  
  # pmcidsFilter intersect to process what is available, setdiff to download what is missing
  if(!is.null(pmcidsFilter) && pmcidsFilter == "intersect") {
    pmcids <- intersect(pmcids, cachePmcids)
  }
  
  # This is changed by forceFlag after 1 run
  forceState <- TRUE 
  
  while(forceState) { tryCatch({
    if(!is.null(pmcidsFilter) && pmcidsFilter == "setdiff") {
      pmcids <- setdiff(pmcids, cachePmcids)
    }
    
    maxPmcids <- length(pmcids)
    allGenes <- list()
    idxRange <- 1:maxPmcids
    
    if(randomizePmcids) {
      idxRange <- sample(idxRange)
    }
    
    for(i in idxRange) {
      # i <- 1
      cat("I: ", i, "\n")
      query <- pmcids[i]
      
      Sys.setenv("PREFIX_SIMPLERCACHE"=query)
      
      q2 <- url_escape(query)
      url <- paste0('https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson?pmcids=',q2)
      url
      cat("URL: ", url, "\n")
      
      ignoreQueriesFilename <- "empty_pmcids.txt"
      
      if(file.exists(ignoreQueriesFilename)) {
        ignoreQueries <- readLines(ignoreQueriesFilename)
      } else {
        ignoreQueries <- NULL
      }
      
      if(query %in% ignoreQueries) {
        next
      }
      
      # Get hash
      resultsHash <- capture.output({
        req <- url %>% GETCached() 
        status <- status_code(req)
        tmpPubtatorDat <- httr::content(req, as="text", encoding="UTF-8")
      })
      cat("DEBUG: ", paste(resultsHash, collapse="\n"), "\n")
      
      resultsHash <- resultsHash[1] %>% trimws %>% sub("keyHash:  ", "", .)
      cat("HASH: ", resultsHash, "\n")
      
      writeLines(tmpPubtatorDat, "pubtator.txt")
      tmpPubtatorDat <- readLines("pubtator.txt")
      #tmpPubtatorDat
      
      tryCatch({
        if(tmpPubtatorDat == "" && status == 200) {
          ignoreQueries <- c(ignoreQueries, query) %>% unique 
          writeLines(ignoreQueries, ignoreQueriesFilename)
          
          badCacheFile <- dir(cacheDir, pattern=resultsHash)
          file.remove(file.path(cacheDir, badCacheFile))
        }
        
        if(length(tmpPubtatorDat) > 0 && tmpPubtatorDat != "") {
          
          if(mainSectionGenesOnly) {
            
            if(abstractOnly) {
              jq_mainsection_str <- '.passages|map(select((.infons."section_type"=="ABSTRACT")))|map(.annotations[])|map(select(.infons.type=="Gene"))|map(.infons.identifier)|unique'
            } else {
              jq_mainsection_str <- '.passages|map(select((.infons."section_type"!="REF") and (.infons."section_type"!="SUPPL") and (.infons."section_type"!="METHODS")))|map(.annotations[])|map(select(.infons.type=="Gene"))|map(.infons.identifier)|unique'
            }
            
            jq_not_mainsection_str <- '.passages|map(select((.infons."section_type"=="REF") or (.infons."section_type"=="SUPPL") or (.infons."section_type"=="METHODS")))|map(.annotations[])|map(select(.infons.type=="Gene"))|map(.infons.identifier)|unique' 

            genes <- jq(tmpPubtatorDat[1], jq_not_mainsection_str) %>% fromJSON
            cat("DEBUG: GENE COUNT OUTSIDE MAIN SECTIONS: ", length(genes), "\n")   
            
            genes <- jq(tmpPubtatorDat[1], jq_mainsection_str) %>% fromJSON
            cat("DEBUG: GENE COUNT INSIDE MAIN SECTIONS ", length(genes), "\n")
          } else {
            pubtator <- fromJSON(tmpPubtatorDat[1])
            pubtator$accessions
            idx <- grepl("gene", pubtator$accessions)
            
            genes <- pubtator$accessions[idx] %>%
              gsub("gene@", "", .) %>%
              strsplit(., ";" ) %>%
              unlist %>%
              unique            
          }
        } else {
          genes <- NULL
          badCacheFile <- dir(cacheDir, pattern=resultsHash)
          file.remove(file.path(cacheDir, badCacheFile))
        }
        
        cat("DEBUG: FINAL GENE COUNT: ", length(genes), "\n")
        genes
      }, error = function(e) {
        cat("ERROR: Retrieval issue:", e$message, "\n")
      })
      
      if(length(genes) > 0) {
        allGenes[[query]] <- genes        
      }
      
      Sys.sleep(sleepTime)
    }
    
    if(!forceFlag) {
      forceState <- FALSE
    }
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
  }) }
  
  saveRDS(allGenes, paste0(term, "_genes.rds"))
  #allGenes <- readRDS(paste0(term, "_genes.rds"))  
  
  return(allGenes)
}
