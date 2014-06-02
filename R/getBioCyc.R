##' Get species information for BioCyc.
##'
##' Get the BioCyc species information including the BioCyc species ID and the Latin name.
##' @title Get species from BioCyc
##' @param speList The species list that is a vector like 'c("HUMAN", "ECOLI", "ZMOB579138")'. The input speList should be consistent with the parameter 'speType'.
##' @param speType It supports two types: 'BioCyc' and 'regexpr'.
##' BioCyc type is the BioCyc species ID, for exmaple 'HUMAN' is the BioCyc ID for the Homo sapiens.The 'regexpr' is used for regulare expression search with the Latin name for example 'Escherichia coli'.
##' @param whole Whether or not get the whole BioCyc species list,
##' and the default value is FALSE.
##' @return Matrix of species information.
##' @examples
##' # search species list from BioCyc ID
##' getPhyloCyc(c("HUMAN", "ECOLI", "ZMOB579138"), speType = 'BioCyc')
##' # search species whose names include 'Escherichia coli'
##' getPhyloCyc('Escherichia coli', speType = 'regexpr')
##' # get whole BioCyc species information table
##' getPhyloCyc(whole = TRUE)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML xmlRoot xmlTreeParse getNodeSet
##' @export
##'
getPhyloCyc <- function(speList, speType = 'BioCyc', whole = FALSE) {

  require(XML)

  if (!(speType %in% c('BioCyc', 'regexpr'))) {
    stop('"speType" now only supports "BioCyc" and "regexpr".')
  } else {}

  # read in the whole biocyc XML file
  cycSpeXML <- xmlRoot(xmlTreeParse('http://biocyc.org/xmlquery?dbs'))

  # get each species
  cycSpe <- getNodeSet(cycSpeXML, '//PGDB')

  # get biocyc ID for each species
  cycSpeID <- t(sapply(cycSpe, xmlAttrs))

  # get latin name of each species
  cycLatin <- sapply(cycSpe, function(x) {
    sapply(getNodeSet(x, '//PGDB/*'), xmlValue)
  })
  cycLatin <- sapply(cycLatin, paste, collapse = ' ')
  cycSpeMat <- cbind(cycSpeID[, 1], cycLatin, cycSpeID[, 2])
  colnames(cycSpeMat) <- c('BioCycID', 'LatinName', 'Version')

  if (whole) {
    return(cycSpeMat)
  } else {
    if (speType == 'BioCyc') {
      cycSpeMat <- cycSpeMat[cycSpeMat[, 1] %in% speList, , drop = FALSE]
    }
    else if (speType == 'regexpr') {
      cycSpeMat <- cycSpeMat[grep(speList, cycSpeMat[, 2]), , drop = FALSE]
    }
  }

  return(cycSpeMat)

}



##' Get the NCBI taxonomy ID from a given BioCyc ID
##'
##' NCBI taxonomy ID is used as unique ID accoss cyc and BioCyc databases. This functions is used to get the corresponding NCBI Taxonomy ID from BioCyc. Only one 'cycID' should be input once. It is easy to batch input 'cycID' by using the function sapply().
##' @title Get NCBI Taxonomy ID From BioCyc ID
##' @param cycID The cyc species ID. The KEGG support multiple species ID, for example c('HUMAN', 'MOUSE', 'ECOLI').
##' @param n The number of CPUs or processors, and the default value is 4.
##' @return The corresponding NCBI Taxonomy ID in character vector.
##' @examples
##' # get human, mouse, and Ecoli NCBI taxonomy ID
##' cyc2Tax(c('HUMAN', 'MOUSE', 'ECOLI'))
##' # transfer all BioCyc species ID to NCBI taxonomy ID
##' \dontrun{
##' wBioCyc <- getPhyloCyc(whole = TRUE)
##' wNCBISpe <- cyc2Tax(wBioCyc[, 1])
##' }
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl getURL
##' @importFrom doMC registerDoMC
##' @importFrom foreach foreach
##' @export
##'
cyc2Tax <- function(cycID, n = 4){

  require(RCurl)
  require(foreach)
  require(doMC)
  registerDoMC(n)

  getcontent <- function(s,g) {
    substring(s,g,g+attr(g,'match.length')-1)
  }

  getSingleTax <- function(cycSpeID) {
    # USE: get the NCBI taxonomy ID from one BioCyc species ID
    # INPUT: 'cycSpeID' is one BioCyc species ID
    # OUTPUT: the NCBI taxnomy ID

    # get cycID webpage
    cycLink <- paste('http://biocyc.org/', cycSpeID, '/organism-summary?object=', cycSpeID, sep = '')
    cycWeb <- getURL(cycLink)

    # get Taxonomy ID. The taxonomy ID is in the web-link like 'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=593907'
    taxIDLink <- gregexpr('wwwtax\\.cgi\\?mode=Info\\&id=\\d+', cycWeb)
    taxIDLink <- getcontent(cycWeb, taxIDLink[[1]])
    taxID <- gregexpr('\\d+', taxIDLink)
    taxID <- getcontent(taxIDLink, taxID[[1]])

    return(taxID)
  }

  cycTax <- foreach(i = 1:length(cycID), .combine = c) %dopar% {
    print(paste('It is running ', i, '.', sep = ''))
    singleTaxID <- getSingleTax(cycID[i])
    names(singleTaxID) <- cycID[i]
    return(singleTaxID)
  }

  return(cycTax)

}


##' Get the whole genes/proteins list from a given species
##'
##' Get the BioCyc gene ID or protein ID list of a given species.
##' @title Get whole genes/proteins list
##' @param speID The BioCyc species ID, for example 'ECOLI' is for 'Escherichia coli K-12 substr. MG1655'.
##' @param type Get the 'genes' or 'proteins', and the default value is 'genes'.
##' @return A vector of genes of proteins with BioCyc ID.
##' @examples getCycGenesList('ECOLI')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML xmlRoot xmlTreeParse getNodeSet xmlGetAttr
##' @export
##'
getCycGenesList <- function(speID, type = 'genes'){

  require(XML)

  # read in the whole biocyc XML file
  url <- paste('http://biocyc.org/xmlquery?query=[x:x<-', speID, '^^', type, ']&detail=none', sep = '')
  cycListXML <- xmlRoot(xmlTreeParse(url))

  # get gene/protein information
  if (type == 'genes') {
    cycList <- getNodeSet(cycListXML, '//Gene')
  } else if (type == 'proteins') {
    cycList <- getNodeSet(cycListXML, '//Protein')
  }
  # delect the ones with 'class = "true"', for example "<Gene resource="getxml?HUMAN:BC-1.7.29" orgid="HUMAN" frameid="BC-1.7.29" class="true"/>".
  cycList <- sapply(cycList, function(x){
    if (is.null(xmlGetAttr(x, 'class'))) {
      cycID <- xmlGetAttr(x, 'frameid')
    } else {
      cycID <- NULL
    }
    return(cycID)
  })

  return(unlist(cycList))

}

##' Get gene information from BioCyc database.
##'
##' The gene information from BioCyc including genome location and gene name information. Some genes in BioCyc may do not have common names or accession names In this circumstance, 'NULL' will be return
##' @title Get one gene information from BioCyc.
##' @param geneID A BioCyc gene ID with the length of 1.
##' @param speID The BioCyc species ID, for example 'ECOLI' is for 'Escherichia coli K-12 substr. MG1655'.
##' @return A list.
##' @examples
##' # get 'atpE' gene information from Ecoli K-12 MG1655 strain.
##' getCycGeneInfo('EG10102', 'ECOLI')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML xmlRoot xmlTreeParse getNodeSet
##' @export
##'
getCycGeneInfo <- function(geneID, speID){

  require(XML)

  # read in gene information XML
  url <- paste('http://biocyc.org/getxml?', speID, ':',geneID, '&detail=full', sep = '')
  geneInfoXML <- xmlRoot(xmlTreeParse(url))

  # location in genome
  # direction
  dirTrans <- xmlNodeVal(geneInfoXML, '//transcription-direction')
  # right end
  rightPos <- xmlNodeVal(geneInfoXML, '//right-end-position')
  # left end
  leftPos <- xmlNodeVal(geneInfoXML, '//left-end-position')
  locTrans <- list(dirTrans = dirTrans,
                   rightPos = rightPos,
                   leftPos = leftPos)

  # gene names, some genes may not have names
  # common name
  comName <- xmlNodeVal(geneInfoXML, '//common-name')
  if (length(comName) == 0) {
    comName = NULL
  }
  # accession, the accession name of node is like 'accession-1', 'accession-2'
  allNodeName <- sapply(getNodeSet(geneInfoXML, '//*'), xmlName)
  accNodeName <- allNodeName[grepl('^accession-\\d+', allNodeName)]
  # some genes may have no accession name
  if (length(accNodeName) == 0) {
    accName = NULL
  } else {
    accNodePath <- paste('//', accNodeName, sep = '')
    accName <- unname(sapply(accNodePath, xmlNodeVal, xmlFile = geneInfoXML))
  }
  geneName <- list(comName = comName,
                   accName = accName)

  # merge
  cycGene <- list(locTrans = locTrans,
                  geneName = geneName)

  return(cycGene)

}



##' Get transcription unit (TU) from gene
##'
##' Get TU from a given BioCyc gene ID.
##' @title Get TU from ID. If the given gene has no TU, 'NULL' will be returned. If the 'evidence' is set to TRUE, a list will return.
##' @param geneID A BioCyc gene.
##' @param speID The BioCyc species ID, for example 'ECOLI' is for 'Escherichia coli K-12 substr. MG1655'.
##' @param evidence  Logical value indicates whether to return the evidence value.
##' @return A vector contains TUs or NULL
##' @examples getCycTUfGene('EG10102', 'ECOLI')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML
##' @export
##'
getCycTUfGene <- function(geneID, speID, evidence = FALSE) {

  require(XML)

  # read in TU information XML
  url <- paste('http://biocyc.org/apixml?fn=transcription-units-of-gene&id=', speID, ':', geneID, '&detail=low', sep = '')
  TUInfoXML <- xmlRoot(xmlTreeParse(url))

  # get TU
  TUID <- xmlNodeAttr(TUInfoXML, '/ptools-xml/Transcription-Unit', 'frameid')
  # may have no TU
  if (length(TUID) == 0) {
    TUID = NULL
  } else {}

  return(TUID)
}



##' Get TU information from BioCyc database.
##'
##' Get TU information including genes in TU and evidence.
##' There is another way to get the genes 'http://biocyc.org/apixml?fn=transcription-unit-genes&id=ECOLI:TU0-42328&detail=full'.
##' @title Get one TU information
##' @param TUID A BioCyc TU ID with the length of 1.
##' @param speID The BioCyc species ID, for example 'ECOLI' is for 'Escherichia coli K-12 substr. MG1655'.
##' @return A list of TU information.
##' @examples getCycTUInfo('TU0-6636', 'ECOLI')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML xmlRoot xmlTreeParse
##' @export
##'
getCycTUInfo <- function(TUID, speID) {

  require(XML)

  # read in gene information XML
  url <- paste('http://biocyc.org/apixml?fn=transcription-units-of-gene&id=', speID, ':', TUID, '&detail=full', sep = '')
  TUInfoXML <- xmlRoot(xmlTreeParse(url))

  # genes in the TU
  geneIDs <- xmlNodeAttr(TUInfoXML, '//component/Gene', 'frameid')

  # evidence
  TUEv <- xmlNodeAttr(TUInfoXML, '//evidence/Evidence-Code', 'frameid')

  # merge
  cycTU <- list(geneIDs = geneIDs,
                TUEv = TUEv)

  return(cycTU)

}



##' Translate KEGG ID to BioCyc ID.
##'
##' Translate the KEGG gene ID to BioCyc gene ID is tricky. In BioCyc, the gene names is not in a uniform; some of them use symbol like 'dnaK', but some use the KEGGID. At first, transfer the KEGG ID to symbol. Secondly, If symbol is used, we extract the information from a website like 'http://websvc.biocyc.org/ECOLI/foreignid?ids=b3734'. Otherwise, we use 'http://biocyc.org/xmlquery?query=[x:x<-AACT754507^^genes,x^name="ANH9381_0646"]&detail=full'. There are two circumstances that will return "0": one is that  BioCyc database may marker some genes as "Pseudo-Genes", and the other is different gene symbols in KEGG and BioCyc.
##' @title Transfer KEGG ID to BioCyc ID.
##' @param KEGGID Only one KEGG ID
##' @param speKEGGID Species BioCyc ID.
##' @param speCycID Species KEGG ID
##' @param type 'gene' or 'protein'
##' @return The BioCyc gene ID or "0", if gene is not found.
##' @examples
##' KEGGID2CycID('b3732', 'eco', 'ECOLI')
##' KEGGID2CycID('ANH9381_0646', 'aao', 'AACT754507', type = 'protein')
##' # return '0' because of different gene symbol ('atpE/H' in BioCyc, but 'atpE-H'in KEGG), and the actural gene ID is 'GSGK-10' in BioCyc.
##' KEGGID2CycID('Bd0010', 'bba', 'BBAC264462-WGS')
##' # return '0' because of 'Pseudo-Genes', and the actural gene ID is 'GJTC-3643' in BioCyc
##' KEGGID2CycID('LIV_2438', 'liv', 'LIVA881621')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML xmlRoot xmlTreeParse
##' @export
##'
KEGGID2CycID <- function(KEGGID, speKEGGID, speCycID, type = 'gene') {

  require(XML)

  # transfer KEGG ID to symbol, if it has one; otherwise, we just use the KEGGID
  KEGGsymTable <- webTable(paste('http://rest.kegg.jp/list/', speKEGGID, ':', KEGGID, sep = ''), n = 2)
  KEGGsym <- KEGGsymTable[1, 2]
  KEGGsym <- unlist(strsplit(KEGGsym, split = ';', fixed = TRUE))
  if (!grepl(' ', KEGGsym)[1]) {
    # The first element is gene symbol, because the gene symbol has not text space
    KEGGID <- KEGGsym[1]
  } else {}

  # try symbol
  url <- paste('http://websvc.biocyc.org/', speCycID, '/foreignid?ids=', KEGGID, sep = '')
  genePage <- webTable(url, n = 1)

  if (genePage[2, 1] == '1') {
    if (type == 'gene') {
      cycID <- genePage[3, 1]
    }
    else if (type == 'protein') {
      cycID <- genePage[6, 1]
    }
  } else if (genePage[2, 1] == '0'){
    url <- paste('http://biocyc.org/xmlquery?query=[x:x<-', speCycID, '^^genes,x^name=', '"', KEGGID, '"', ']&detail=full', sep = '')
    geneXML <- xmlRoot(xmlTreeParse(url))
    if (type == 'gene') {
      cycID <- xmlNodeAttr(geneXML, '/ptools-xml/Gene', 'frameid')
    }
    else if (type == 'protein') {
      cycID <- xmlNodeAttr(geneXML, '//Protein', 'frameid')
    }
  }

  if (length(cycID) == 0) {
    cycID = '0'
  } else {}

  return(cycID)

}

##' Get node values.
##'
##' Get node values from the XML file.
##' @title Get all nodes value directly from the XML file.
##' @param xmlFile XML file.
##' @param nodePath The XPath of nodeset (one or mutiple nodes).
##' @return Nodeset value
##' @examples
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML getNodeSet xmlValue
xmlNodeVal <- function(xmlFile, nodePath){

  require(XML)
  nodeSet <- getNodeSet(xmlFile, nodePath)
  nodeValue <- sapply(nodeSet, xmlValue)

  return(nodeValue)

}

##' Get attribute values.
##'
##' Get attribute values from the XML file.
##' @title Get all nodes attributes directly from the XML file.
##' @param xmlFile XML file.
##' @param nodePath The XPath of nodeset (one or mutiple nodes).
##' @param attrName Attributes name
##' @return Nodeset attributes
##' @examples
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML getNodeSet xmlValue
xmlNodeAttr <- function(xmlFile, nodePath, attrName){

  require(XML)
  nodeSet <- getNodeSet(xmlFile, nodePath)
  nodeAttr <- sapply(nodeSet, xmlGetAttr, name = attrName)

  return(nodeAttr)

}




## ##' Transfer KEGG ID to BioCyc ID.
## ##'
## ##' Tranfser KEGG ID to BioCyc ID by the route, whole species gene list --> unit gene information --> match KEGG ID
## ##' @title Transfer genes ID from KEGG to BioCyc
## ##' @param KEGGID A vector of KEGG gene ID.
## ##' @param speID Species BioCyc ID.
## ##' @param n The number of CPUs or processors, and the default value is 4.
## ##' @param type Only support genes
## ##' @return A vector of BioCycID
## ##' # EG10098 EG10101
## ##' # "b3734" "b3732"
## ##' @examples KEGGID2CycID(c('b3734', 'b3732'), 'ECOLI')
## ##' @author Yulong Niu \email{niuylscu@@gmail.com}
## ##' @importFrom doMC registerDoMC
## ##' @importFrom foreach foreach
## ##' @export
## ##'
## KEGGID2CycID <- function(KEGGID, speID, n = 4, type = 'genes'){

##   require(foreach)
##   require(doMC)

##   registerDoMC(n)

##   # whole list of given species
##   geneList <- getCycGenesList(speID)

##   # whole accession number
##   accList <- foreach(i = 1:length(geneList)) %dopar% {
##     print(paste('It is running ', i, '.', sep = ''))
##     geneInfo <- getCycGeneInfo(geneList[i], speID)
##     accInfo <- geneInfo$geneName$accName
##     # some genes have no accession names
##     if (!is.null(accInfo)) {
##       names(accInfo) <- rep(geneList[i], length(accInfo))
##     } else {}
##     return(accInfo)
##   }

##   accList <- unlist(accList)

##   KEGGList <- accList[accList %in% KEGGID]

##   return(KEGGList)


## }
