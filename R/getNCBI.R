##' NCBI Database Additional API - Get NCBI taxonomy information from a given NCBI taxonomy IDs
##'
##' Get NCBI taxonomy information. 
##' @title Get NCBI taxonomy information
##' @param NCBITaxoIDs A vector of NCBI taxonomy IDs.
##' @return A list containing taxonomy information for each ID.
##' @examples
##' threeTax <- getNCBITaxo(c("9606", "511145", "797302"))
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom xml2 read_xml xml_children xml_contents
##' @export
##'
getNCBITaxo <- function(NCBITaxoIDs) {

  ## compress taxonomy IDs
  taxoIDs <- paste(NCBITaxoIDs, collapse = ',')

  ## read in xml content
  url <- paste0('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=', taxoIDs, '&retmode=xml')
  taxXML <- read_xml(url)

  taxChildren <- xml_children(taxXML)

  taxList <- lapply(taxChildren, function(x){
    ## input TaxId is the 1st child, ScientificName is 2nd, and Rank is the 5th
    inputTax <- xml_contents(xml_children(x)[c(1, 2, 5)])
    inputTax <- as.character(inputTax)

    ## higher tax, LineageEx is 10th child
    higherEle <- xml_children(xml_children(x)[10])
    higherTax <- sapply(higherEle, function(y){
      higherEach <- xml_contents(xml_children(y))
      higherEach <- as.character(higherEach)

      return(higherEach)
    })
    higherTax <- t(higherTax)

    ## combine higher and input tax
    allTax <- rbind(higherTax, inputTax, deparse.level = 0)

    return(allTax)
  })

  return(taxList)
}
