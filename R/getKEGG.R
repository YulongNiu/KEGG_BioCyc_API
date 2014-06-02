##' Get species KEGG/NCBI ID.
##'
##' Get the phylogenetic information of given species.
##' It supports both batch input and regular expression search.
##'
##' @title Get species list from KEGG.
##' @name getSpePhylo
##' @param speList The species list that is a vector like 'c("hsa", "eco")'. The input 'speList' should be consistent with the parameter 'speType'.
##' @param speType It supports five types: 'KEGG', 'Tnum', 'regexpr', and 'phylo'.
##' KEGG type is a three or four letters, for exmaple 'hsa' is the KEGG ID for Homo sapiens,
##' while the corresponding T number is 'T01001'.
##' The 'regexpr' is used for regulare expression search with the Latin name ('Escherichia coli'), sub-species name ('K-12 MG1655'), and common name ('human').
##' The 'phylo' uses phylogentic orders for search, and it supports 'Domain' (either 'Eukaryotes' or 'Prokaryotes'), 'Kingdom' ('Animals'), 'phylum' ('Vertebrates'), and 'class' ('Mammals').
##' But it does not support mixed 'KEGG'.
##' Attention: Mutiple KEGG species ID may correspond to one taxonomy ID, for exmaple 'lph' and 'lpo' to '91891'
##' @param whole Whether or not get the whole KEGG species list,
##' and the default value is FALSE.
##' @return Matrix of species information.
##' @examples
##' # search species list from KEGG ID
##' getSpePhylo(c('hsa', 'eco'))
##' # search species whose names include 'Escherichia coli'
##' getSpePhylo('Escherichia coli', speType = 'regexpr')
##' # search species whose class is 'Mammals'
##' getSpePhylo('Mammals', speType = 'phylo')
##' # get whole KEGG species information table
##' getSpePhylo(whole = TRUE)
##'
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
##'
getSpePhylo <- function(speList, speType = 'KEGG', whole = FALSE){

  if (!(speType %in% c('KEGG', 'regexpr', 'phylo', 'Tnum'))) {
    stop('"speType" now only supports "NCBI", "KEGG", "Tnum", "regexpr", and "phylo".')
  } else {}

  # get whole KEGG species information
  speAnno <- webTable('http://rest.kegg.jp/list/organism', ncol = 4)

  colnames(speAnno) <- c('TID', 'KEGGID', 'LatinName', 'Phylo')

  if (whole) {
    return(speAnno)
  } else {
    if (speType == 'Tnum') {
      selSpeAnno <- speAnno[speAnno[, 1] %in% speList, , drop = FALSE]
    }
    else if (speType == 'KEGG') {
      selSpeAnno <- speAnno[speAnno[, 2] %in% speList, , drop = FALSE]
    }
    else if (speType == 'regexpr') {
      selSpeAnno <- speAnno[grep(speList, speAnno[, 3]), , drop = FALSE]
    }
    else if (speType == 'phylo') {
      selSpeAnno <- speAnno[grep(speList, speAnno[, 4]), , drop = FALSE]
    }
  }

  return(selSpeAnno)

}



##' Get the NCBI taxonomy ID from a given KEGG ID
##'
##' NCBI taxonomy ID is used as unique ID accoss KEGG and BioCyc databases. This functions is used to get the corresponding NCBI Taxonomy ID from KEGG.
##' @title Get NCBI Taxonomy ID From KEGG ID
##' @param KEGGID The KEGG support multiple species ID, for example c('hsa', 'eco').
##' @param n The number of CPUs or processors, and the default value is 4.
##' @return The corresponding NCBI Taxonomy ID in character vector.
##' @examples
##' # get human and Ecoli NCBI taxonomy ID
##' KEGG2Tax(c('hsa', 'eco'))
##' # transfer all KEGG species ID to NCBI taxonomy ID
##' \dontrun{
##' wKEGGSpe <- getSpePhylo(whole = TRUE)
##' wNCBISpe <- KEGG2Tax(wKEGGSpe[, 2])
##' }
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl getURL
##' @importFrom doMC registerDoMC
##' @importFrom foreach foreach %dopar%
##' @export
##'
KEGG2Tax <- function(KEGGID, n = 4){

  require(RCurl)
  require(foreach)
  require(doMC)
  registerDoMC(n)

  getcontent <- function(s,g) {
    substring(s,g,g+attr(g,'match.length')-1)
  }

  getSingleTax <- function (KEGGspeID) {
    # USE: get KEGGSpeID webpage
    # INPUT: 'KEGGID' is the KEGG species ID.
    # OUTPUT: The NCBI taxonomy ID.
    KEGGLink <- paste('http://www.genome.jp/kegg-bin/show_organism?org=', KEGGspeID, sep = '')
    KEGGWeb <- getURL(KEGGLink)

    # get Taxonomy ID. The taxonomy ID is in the web-link like 'http://www.ncbi.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=593907'
    taxIDLink <- gregexpr('wwwtax\\.cgi\\?id=\\d+', KEGGWeb)
    taxIDLink <- getcontent(KEGGWeb, taxIDLink[[1]])
    taxID <- gregexpr('\\d+', taxIDLink)
    taxID <- getcontent(taxIDLink, taxID[[1]])

    return(taxID)
  }

  NCBITax <- foreach(i = 1:length(KEGGID), .combine = c) %dopar% {
    print(paste('It is running ', i, '.', sep = ''))
    taxID <- getSingleTax(KEGGID[i])
    names(taxID) <- KEGGID[i]
    return(taxID)
  }

  return(NCBITax)

}

##' Get the KEGG orthology list.
##'
##' Get the KEGG orthology list by a given KEGG KO ID.
##' @title Get KEEE orthology.
##' @param KOID The KEGG orthology ID.
##' @return A matrix of genes under the given orthology ID.
##' @examples getKEGGKO('K02110')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
getKEGGKO <- function(KOID){

  # get KO webpage
  url <- paste('http://rest.kegg.jp/link/genes/', KOID, sep = '')

  # get KO list
  KOWebMat <- webTable(url, ncol = 2)

  # seperate the species and genes
  KOMat <- unlist(strsplit(KOWebMat[, 2], split = ':', fixed = TRUE))
  KOMat <- matrix(KOMat, ncol = 2, byrow = TRUE)
  colnames(KOMat) <- c('speID', 'geneID')

  return(KOMat)

}


##' Get the whole pathway ID from KEGG database.
##'
##' Get the pathway ID and annoation of a given KEGG species ID.
##' @title List pathway of a given species ID
##' @param KEGGspec The KEGG species ID. Only one species ID once input.
##' @return A matrix of pathway ID and annotation.
##' @examples getKEGGPathAnno('hsa')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
getKEGGPathAnno <- function(KEGGspec){

  # get KEGG pathway annotation list
  url <- paste('http://rest.kegg.jp/list/pathway/', KEGGspec, sep = '')
  pathAnno <- webTable(url, ncol = 2)

  # the transfer the pathname of 'path:hsa00010' to 'hsa00010'.
  pathID <- pathAnno[, 1]
  pathID <- sapply(strsplit(pathID, split = ':', fixed = TRUE), '[', 2)
  pathAnno[, 1] <- pathID

  colnames(pathAnno) <- c('pathID', 'Annotation')

  return(pathAnno)
}


##' Get the pathway and genes.
##'
##' Get the pathway and genes according to KEGG species ID.
##' @title List pathways and genes of a given KEGG species ID
##' @param KEGGspec The KEGG species ID. Only one species ID once input.
##' @return A List named with KEGG pathway IDs, and each element of the list contains the KEGG gene IDs.
##' @examples getKEGGPathGenes('hsa')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
getKEGGPathGenes <- function(KEGGspec){

  # get KEGG pathway list
  url <- paste('http://rest.kegg.jp/link/', KEGGspec, '/pathway',  sep = '')

  pathAnno <- webTable(url, ncol = 2)

  pathInfo <- aggregate(pathAnno[, 2], list(pathAnno[, 1]), '[')
  pathList <- list()
  # transfer into list
  for (i in 1:nrow(pathInfo)){
    pathList[[i]] <- as.character(pathInfo[i, 2][[1]])
  }
  names(pathList) <- pathInfo[, 1]

  return(pathList)
}


##' Get a R matrix object if the weblink returned as a matrix.
##'
##' If the web return a matrix, use this function to extract it as a R matrix object.
##' @title Get R matrix from weblink
##' @param url The weblink.
##' @param ncol The column number of the matrix.
##' @return A R matrix
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl getURL
webTable <- function(url, ncol) {

  require(RCurl)

  webPage <-getURL(url)

  # transfer webpage into a matrix
  webMat <- unlist(strsplit(webPage, split = '\n', fixed = TRUE))
  webMat <- sapply(webMat, strsplit, split = '\t', fixed = TRUE)
  webMat <- matrix(unlist(webMat), ncol = ncol, byrow = TRUE)

  return(webMat)
}
