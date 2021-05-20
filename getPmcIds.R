library(xml2)
library(httr)
library(magrittr)
library(jsonlite)

#' Retrieve PMC IDs for a search term
#' 
#' @param term the search term (e.g., sarcoma)
#' @param maxIds max results to return
#' 
#' @return a vector of PMC IDs
#' 
#' @examples 
#' #getPubtator(term="sarcoma", pmcids, sleepTime=0, pmcidsFilter="intersect", forceFlag=FALSE, randomizePmcids=FALSE)
#'
#' @export
getPmcIds <- function(term, maxIds=10000) {
  # TODO maxIds not working
  
  #query <- paste0('("', term, '") AND (pubmed pmc [sb]) AND hasabstract')
  query <- paste0('(', term, ') AND (pubmed pmc [sb]) AND hasabstract')
  q2 <- url_escape(query)
  url <- paste0('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&usehistory=y&datetype=crdt&retmax=', maxIds, '&rettype=XML&term=',q2)
  url
  cat("URL: ", url, "\n")
  
  tmpEsearchDat <- url %>% GET() %>% httr::content("text")
  tmpEsearchDat
  
  esearchDat <- read_xml(tmpEsearchDat, verbose=TRUE)
  write_xml(esearchDat, paste0(term, "_esearch.xml"))
  
  ## Use PubMed WebEnv to retrieve results
  t1 <- xml_find_all(esearchDat, "//WebEnv")
  webenv <- xml_text(t1)
  webenv
  
  # Fetch PubMed content
  url <- paste0('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&query_key=1&rettype=xml&WebEnv=',url_escape(webenv))
  url
  cat("URL: ", url, "\n")
  
  tmpEfetchDat <- url %>% GET() %>% httr::content("text")
  #tmpEfetchDat
  #efetchDat <- read_xml("efetch.xml")
  efetchDat <- read_xml(tmpEfetchDat, verbose=TRUE)
  write_xml(efetchDat, paste0(term, "_efetch.xml"))
  
  pmcids <- xml_text(xml_find_all(efetchDat, ".//ArticleId[@IdType='pmc']"))
  pmcids
  writeLines(pmcids, paste0(term, "_pmcids.txt"))
  
  return(pmcids)
}
