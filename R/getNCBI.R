##' NCBI Database API - Get NCBI taxonomy information from given NCBI taxonomy IDs
##'
##' Get NCBI taxonomy information. 
##' @title Get NCBI taxonomy information
##' @param NCBITaxoIDs A vector of NCBI taxonomy IDs.
##' @return A list containing taxonomy information for each ID.
##' @examples
##' threeTax <- getNCBITaxo(c('9606', '511145', '797302'))
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl postForm
##' @importFrom xml2 read_xml xml_children xml_contents
##' @references Entrez Programming Utilities Help \url{http://www.ncbi.nlm.nih.gov/books/NBK25499/}
##' @export
##'
getNCBITaxo <- function(NCBITaxoIDs) {

  ## compress taxonomy IDs
  taxoIDs <- paste(NCBITaxoIDs, collapse = ',')

  ## read in xml content with HTTP POST method
  urlBase <- 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
  xmlStr <- postForm(urlBase,
                     db = 'taxonomy',
                     id = taxoIDs,
                     retmode = 'xml')
  taxXML <- read_xml(xmlStr)

  ## extract child node of each tax
  taxChildren <- xml_children(taxXML)
  taxList <- lapply(taxChildren, function(x){

    ## all TaxId, all ScientificName, and all Rank
    xPathVec <- c('.//TaxId', './/ScientificName', './/Rank')
    allTax <- sapply(xPathVec, function(eachXPath) {
      eachNodeSet <- xml_find_all(x, eachXPath)
      eachContent <- xml_contents(eachNodeSet)
      eachContent <- as.character(eachContent)
    })

    ## move input tax to the last
    allTax <- rbind(allTax[-1, ], allTax[1, ])
    colnames(allTax) <- c('TaxId', 'ScientificName', 'Rank')

    return(allTax)
  })

  names(taxList) <- NCBITaxoIDs

  return(taxList)
}


##' NCBI Database API - Get NCBI gene information from given NCBI gene IDs
##'
##' Get NCBI gene information, including gene name, description, genetic source, aliases, gene location.
##' To retrieve thousands of proteins, use EPost to post record into the web server and then retrieve data using ESummary.
##' @title Get NCBI gene information
##' @param NCBIGeneIDs A vector of NCBI gene IDs.
##' @inheritParams getKEGGGeneSeq
##' @return A list containing gene information for each ID. "NA" will be returned for the items if the contents are not found.
##' @examples
##' twoGenes <- getNCBIGeneInfo(c('948242', '15486644'))
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl postForm
##' @importFrom xml2 read_xml xml_children xml_contents
##' @references Entrez Programming Utilities Help \url{http://www.ncbi.nlm.nih.gov/books/NBK25499/}
##' @export
##'
##' 
getNCBIGeneInfo <- function(NCBIGeneIDs, n = 4) {

  ##~~~~~~~~~~~~~~~~~~~~~~~~~EPost~~~~~~~~~~~~~~~~~~~~~~~
  ## compress gene IDs
  geneIDs <- paste(NCBIGeneIDs, collapse = ',')
  infoPostPara <- list(db = 'gene', id = geneIDs)
  infoPost <- EPostNCBI(infoPostPara)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ##~~~~~~~~~~~~~~~~~~~~~~ESummary~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  return(infoPost)

}


##' NCBI Database API - Directly use NCBI EPost API
##'
##' NCBI EPost provide thousands of queries in one HTTP POST.
##' @title NCBI EPost API
##' @param postPara A named list of HTTP POST terms. For example, "db" is E-utility Database Name, see \url{http://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly}.
##' @param ... Parameters inherited from the "postForm()" function in the "RCurl" package.
##' @return A names list containing "QueryKey", "Count", "WebEnv".
##' @examples
##' genePara <- list(db = "gene", id="948242,15486644")
##' genePost <- EPostNCBI(genePara)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl postForm
##' @importFrom xml2 read_xml xml_children xml_contents
##' @references Entrez Programming Utilities Help \url{http://www.ncbi.nlm.nih.gov/books/NBK25499/}
##' @export
##'
##' 
EPostNCBI <- function(postPara, ...) {
  
  ## NCBI EPost url
  urlBase <- 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi'
  
  ## HTTP POST with EPost
  xmlPostStr <- postForm(uri = urlBase, .params = postPara)
  xmlPost <- read_xml(xmlPostStr)

  ## retrieve key and webenv
  xPathPost <- c('.//QueryKey', './/WebEnv')
  infoPost <- lapply(xPathPost, function(eachXPath) {
    eachNodeSet <- xml_find_all(xmlPost, eachXPath)
    eachContent <- xml_contents(eachNodeSet)
    eachContent <- as.character(eachContent)
  })
  names(infoPost) <- c('QueryKey', 'WebEnv')

  return(infoPost)
}
