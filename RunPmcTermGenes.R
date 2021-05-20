library(magrittr)
library(httr)
library(jsonlite)
library(xml2)
library(simpleRCache)
library(jqr)
library(org.Hs.eg.db)

source("getPmcIds.R")
source("getPubtator.R")
source("getPmcGeneCounts.R")

term <- "osteosarcoma"
cacheDir <- "cache_osteosarcoma"

# GET PMCIDS ----
pmcids <- getPmcIds(term)
tmp_pmcids <- pmcids[1:100] # Subset for example 

# GET PUBTATOR ----
pmcid_gene_lst <- getPubtator(term=term, 
                              cacheDir=cacheDir, 
                              pmcids=tmp_pmcids, 
                              sleepTime=3, 
                              pmcidsFilter="intersect", # change from "setdiff" to "intersect" if files downloaded previously
                              forceFlag=FALSE, 
                              randomizePmcids=TRUE, 
                              mainSectionGenesOnly=TRUE,
                              abstractOnly=FALSE)

# GET GENE COUNTS ----
gene_counts <- getPmcGeneCounts(term=term, pmcidsToKeep=NULL, verbose=FALSE)

