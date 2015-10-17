##' NCBI Database Additional API - Get NCBI taxonomy information from a given NCBI taxonomy IDs
##'
##' Get NCBI taxonomy information. 
##' @title Get NCBI taxonomy information
##' @param NCBITaxoIDs A vector of NCBI taxonomy IDs.
##' @return A list containing taxonomy information for each ID.
##' @examples
##' threeTax <- getNCBITaxo(c("9606", "511145", "797302"))
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl postForm
##' @importFrom xml2 read_xml xml_children xml_contents
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
