##' BioCyc Database API - Get species information from BioCyc.
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
##' getCycPhylo(c('HUMAN', 'ECOLI', 'ZMOB579138'), speType = 'BioCyc')
##' 
##' # search species whose names include 'Escherichia coli'
##' getCycPhylo('Escherichia coli', speType = 'regexpr')
##' 
##' \dontrun{
##' # get whole BioCyc species information table
##' getCycPhylo(whole = TRUE)}
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom xml2 read_xml xml_find_all xml_attrs xml_text xml_children
##' @export
##'
getCycPhylo <- function(speList, speType = 'BioCyc', whole = FALSE) {

  if (!(speType %in% c('BioCyc', 'regexpr'))) {
    stop('"speType" now only supports "BioCyc" and "regexpr".')
  } else {}

  ## read in the whole biocyc XML file
  cycSpeXML <- read_xml('http://biocyc.org/xmlquery?dbs')

  ## get each species
  cycSpe <- xml_find_all(cycSpeXML, '//PGDB')

  ## get biocyc ID for each species
  cycSpeID <- t(sapply(cycSpe, xml_attrs))

  ## get latin name of each species
  cycLatin <- lapply(cycSpe, function(x) {
    eachLatin <- xml_text(xml_children(x))
    return(eachLatin)
  })
  
  cycLatin <- sapply(cycLatin, paste, collapse = ' ')
  cycSpeMat <- cbind(cycSpeID[, 1], cycLatin, cycSpeID[, 2])
  colnames(cycSpeMat) <- c('BioCycID', 'LatinName', 'Version')

  if(!whole) {
    if (speType == 'BioCyc') {
      cycSpeMat <- cycSpeMat[cycSpeMat[, 1] %in% speList, , drop = FALSE]
    }
    else if (speType == 'regexpr') {
      cycSpeMat <- cycSpeMat[grep(speList, cycSpeMat[, 2]), , drop = FALSE]
    }
  } else {}

  return(cycSpeMat)
}


##' BioCyc Database API - Get the whole genes/proteins list from a given species
##'
##' Get the BioCyc gene ID or protein ID list of a given species.
##' @title Get whole genes/proteins list
##' @param speID The BioCyc species ID, for example "ECOLI" is for "Escherichia coli K-12 substr. MG1655".
##' @param type Get the "genes" or "proteins", and the default value is "genes".
##' @return A vector of genes of proteins with BioCyc ID.
##' @examples getCycGenesList('ECOLI')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML xmlRoot xmlTreeParse getNodeSet xmlGetAttr
##' @export
##'
getCycGenesList <- function(speID, type = 'genes'){

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

##' BioCyc Database API - Get gene information from BioCyc database.
##'
##' The gene information from BioCyc including genome location and gene name information. Some genes in BioCyc may do not have common names or accession names In this circumstance, "NULL" will be return
##' @title Get one gene information from BioCyc.
##' @param geneID A BioCyc gene ID with the length of 1.
##' @param speID The BioCyc species ID, for example "ECOLI" is for "Escherichia coli K-12 substr. MG1655".
##' @return A list.
##' @examples
##' # get 'atpE' gene information from Ecoli K-12 MG1655 strain.
##' getCycGeneInfo('EG10102', 'ECOLI')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML xmlRoot xmlTreeParse getNodeSet xmlName
##' @export
##'
getCycGeneInfo <- function(geneID, speID){

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
  comName <- testLen(comName, NULL, comName)

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
                  geneName = geneName,
                  url = url)

  return(cycGene)

}



##' BioCyc Database API - Get transcription unit (TU) from gene
##'
##' Get TU from a given BioCyc gene ID.
##' If the given gene has no TU, "NULL" will be returned. If the "evidence" is set to TRUE, a list will return.
##' @title Get TU from ID.
##' @param geneID A BioCyc gene.
##' @param speID The BioCyc species ID, for example "ECOLI" is for "Escherichia coli K-12 substr. MG1655".
##' @param evidence  Logical value indicates whether to return the evidence value.
##' @return A vector contains TUs or NULL
##' @examples getCycTUfGene('EG10102', 'ECOLI')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML xmlRoot xmlTreeParse
##' @export
##'
getCycTUfGene <- function(geneID, speID, evidence = FALSE) {

  # read in TU information XML
  url <- paste0('http://biocyc.org/apixml?fn=transcription-units-of-gene&id=', speID, ':', geneID, '&detail=low')
  TUInfoXML <- xmlRoot(xmlTreeParse(url))

  # get TU
  TUID <- xmlNodeAttr(TUInfoXML, '/ptools-xml/Transcription-Unit', 'frameid')
  TUID <- testLen(TUID, NULL, TUID)

  return(TUID)
}



##' BioCyc Database API - Get TU information from BioCyc database.
##'
##' Get TU information including genes in TU and evidence.
##' There is another way to get the genes 'http://biocyc.org/apixml?fn=transcription-unit-genes&id=ECOLI:TU0-42328&detail=full'.
##' @title Get one TU information
##' @param TUID A BioCyc TU ID with the length of 1.
##' @param speID The BioCyc species ID, for example 'ECOLI' is for 'Escherichia coli K-12 substr. MG1655'.
##' @return A list of TU information.
##' @examples
##' getCycTUInfo('TU0-6636', 'ECOLI')
##' getCycTUInfo('TU00260', 'ECOLI')
##' getCycTUInfo('TUC7Z-43', 'SMUT210007')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML xmlRoot xmlTreeParse
##' @export
##'
getCycTUInfo <- function(TUID, speID) {

  # read in gene information XML
  url <- paste('http://biocyc.org/apixml?fn=transcription-units-of-gene&id=', speID, ':', TUID, '&detail=full', sep = '')
  TUInfoXML <- xmlRoot(xmlTreeParse(url))

  # genes in the TU
  geneIDs <- xmlNodeAttr(TUInfoXML, '//component/Gene', 'frameid')

  # terminator
  terminatorName <- xmlNodeAttr(TUInfoXML, '//component/Terminator', 'frameid')
  terminatorName <- testLen(terminatorName, NULL, terminatorName)
  # right end
  rightPos <- xmlNodeVal(TUInfoXML, '//component/Terminator/right-end-position')
  rightPos <- testLen(rightPos, NULL, rightPos)
  # left end
  leftPos <- xmlNodeVal(TUInfoXML, '//component/Terminator/left-end-position')
  leftPos <- testLen(leftPos, NULL, leftPos)
  # teminator evidence
  TMEv <- xmlNodeAttr(TUInfoXML, '//component/Terminator/evidence/Evidence-Code', 'frameid')
  TMEv <- testLen(TMEv, NULL, TMEv)
  terminator <- list(name = terminatorName,
                     rightPos = rightPos,
                     leftPos = leftPos,
                     Ev = TMEv)

  # promoter
  promoterName <- xmlNodeAttr(TUInfoXML, '//component/Promoter', 'frameid')
  promoterName <- testLen(promoterName, NULL, promoterName)
  # promoter common name
  promoterComName <- xmlNodeVal(TUInfoXML, '//component/Promoter/common-name')
  promoterComName <- testLen(promoterComName, NULL, promoterComName)
  # right end
  rightPos <- xmlNodeVal(TUInfoXML, '//component/Promoter/right-end-position')
  rightPos <- testLen(rightPos, NULL, rightPos)
  # left end
  leftPos <- xmlNodeVal(TUInfoXML, '//component/Promoter/left-end-position')
  leftPos <- testLen(leftPos, NULL, leftPos)
  # promoter evidence
  PMEv <- xmlNodeAttr(TUInfoXML, '//component/Promoter/evidence/Evidence-Code', 'frameid')
  PMEv <- testLen(PMEv, NULL, PMEv)
  promoter <- list(name = promoterName,
                   comName = promoterComName,
                   rightPos = rightPos,
                   leftPos = leftPos,
                   Ev = PMEv)
  
  # TU evidence
  TUEv <- xmlNodeAttr(TUInfoXML, '//Transcription-Unit/evidence/Evidence-Code', 'frameid')
  TUEv <- testLen(TUEv, NULL, TUEv)

  # merge
  cycTU <- list(geneIDs = geneIDs,
                terminator = terminator,
                promoter = promoter,
                TUEv = TUEv,
                url = url)

  return(cycTU)

}


##' BioCyc Database API - Get whole TU list of a given species from BioCyc database.
##' Get the whole transcription units from a given species. It may take more than 10 minutes to retrieve the xml file.
##' 
##' @title Get whole TU information from a given species
##' @param speID The BioCyc species ID, for example "ECOLI" is for "Escherichia coli K-12 substr. MG1655".
##' @return The TU id vector
##' @examples
##' \dontrun{
##' # get Streptococcus mutans UA159 TU
##' SMTU <- getCycSpeTU('SMUT210007')
##' }
##' @author Yulong Niu \email{niuylscu@@gmail.com}
getCycSpeTU <- function(speID){
  
  url <- paste('http://websvc.biocyc.org/xmlquery?[x:x<-', speID, '^^Transcription-Units]', sep = '')

  TUxml <- xmlRoot(xmlTreeParse(url))
  TUvec <- xmlNodeAttr(TUxml, '//Transcription-Unit', 'frameid')
}




##' Translate KEGG ID to BioCyc ID.
##'
##' Translate the KEGG gene ID to BioCyc gene ID is tricky. In BioCyc, the gene names is not in a uniform; some of them use symbol like 'dnaK', but some use the KEGGID. For genes, if symbol is given, we use 'http://biocyc.org/xmlquery?query=[x:x<-ECOLI^^genes,x^name="atpA"]&detail=full'. There are two circumstances that will return "0": one is that  BioCyc database may marker some genes as "Pseudo-Genes", and the other is different gene symbols in KEGG and BioCyc. For proteins, we at first transfer KEGG gene IDs to UniProt IDs, and then to BioCyc gene IDs.
##' 
##' @title Transfer KEGG ID to BioCyc ID.
##' @param KEGGID Only one KEGG ID
##' @param speKEGGID Species BioCyc ID.
##' @param speCycID Species KEGG ID
##' @param type 'gene' or 'protein'
##' @return The BioCyc gene ID or "0", if gene is not found.
##' @examples
##' # symbol is "atpD"
##' transGeneIDKEGG2Cyc('b3732', 'eco', 'ECOLI')
##' 
##' # symbol is "SMU_408" but the first annotation word is "permease"
##' transGeneIDKEGG2Cyc('SMU_408', 'smu', 'SMUT210007')
##' 
##' # It will return "0" because of the symbol 'atpE_H' from KEGG.
##' # The symbol in BioCyc is 'atpE/H'.
##' transGeneIDKEGG2Cyc('Bd0010', 'bba', 'BBAC264462')
##' 
##' # retrieve protein
##' transGeneIDKEGG2Cyc('b0001', 'eco', 'ECOLI', type = 'protein')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom xml2 read_xml xml_children xml_attr
##' @export
##'
transGeneIDKEGG2Cyc <- function(KEGGID, speKEGGID, speCycID, type = 'gene') {


  KEGGID2symbol <- function(KEGGIDtry, speKEGGIDtry) {
    ## transfer KEGG ID to symbol, if it has one; otherwise, we just use the KEGGID
    ## To save the gene symbol at most, the return gene information is split by ";" and ",". All elements, expect one contains space, are returned.
    ## USE: try to convert KEGG gene ID to symbole from KEGG gene annoation. If the first 'proper word' is all in lower case, return it plus KEGGID. If no 'proper word', KEGGID returns.
    ## INPUT: 'KEGGIDtry' is the KEGG gene ID. 'speKEGGIDtry' is the KEGG species ID.
    ## OUTPUT: A vector (length may be bigger than 1)
    ## EXAMPLE: KEGGID2symbol('SMU_23', 'smu')
    ## EXAMPLE: KEGGID2symbol('SMU_t02', 'smu')
    ## EXAMPLE: KEGGID2symbol('SMU_24', 'smu')
    ## EXAMPLE: KEGGID2symbol('SMU_20', 'smu')
    ## EXAMPLE: KEGGID2symbol('SMU_408', 'smu')
    KEGGsymTable <- webTable(paste('http://rest.kegg.jp/list/', speKEGGIDtry, ':', KEGGIDtry, sep = ''), ncol = 2)
    KEGGsym <- KEGGsymTable[1, 2]

    ## split by '; ' and ', '
    ## if KEGGsym is '', then NULL will return. c('test', NULL) is equal to 'test'
    KEGGEle <- unlist(strsplit(KEGGsym, split = '; ', fixed = TRUE))
    KEGGEle <- unlist(strsplit(KEGGEle, split = ', ', fixed = TRUE))

    ## remove elements with space
    KEGGEle <- KEGGEle[!grepl(' ', KEGGEle)]

    ## combine KEGGIDtry
    KEGGEle <- c(KEGGIDtry, KEGGEle)

    return(KEGGEle)
  }

  TrySymConv <- function(symboltry, speCycIDtry) {
    ## USE: try to convert KEGG symbol to Biocyc gene ID
    ## INPUT: 'symboltry' is the gene symbol extract from KEGG. 'speCycIDtry' is Biocyc species ID.
    ## OUTPU: converted BioCyc gene ID. '"0"' will return, if not found.
    ## EXAMPLE: TrySymConv('dnaA', 'SMUT210007')
    ## EXAMPLE: TrySymConv('SMU_06', 'SMUT210007')

    ## ## old code for xml, but will not work because "%" --> "%25"
    ## url <- paste0('http://biocyc.org/xmlquery?query=[x:x<-', speCycIDtry, '^^genes,x^name','="', symboltry, '"]&detail=full')
    ## geneXML <- xmlRoot(xmlTreeParse(url))
    ## cycID <- xmlNodeAttr(geneXML, '/ptools-xml/Gene', 'frameid')
    ## cycID <- testLen(cycID, '0', cycID)
    
    ## try symbol
    url <- paste0('http://websvc.biocyc.org/xmlquery?query=[x:x%3C-',speCycIDtry, '^^genes,x^name%3D%22', symboltry, '%22]&detail=full')
    geneXML <- read_xml(url)
    geneXMLChild <- xml_children(geneXML)
    
    ## check if xml returns geneID
    if (length(geneXMLChild) < 2) {
      cycID = '0'
    } else {
      geneXMLChildCont <- geneXMLChild[[2]]
      cycID <- xml_attr(geneXMLChildCont, 'frameid', default = '0')
    }
    
    cycIDList <- list(cycID = cycID, url = url)
    ## # also get protein
    ## cycID <- xmlNodeAttr(geneXML, '//Protein', 'frameid')
    
    return(cycIDList)
  }

  if (type == 'gene') {
    ## convert to symbol
    symbol <- KEGGID2symbol(KEGGID, speKEGGID)
    for (i in 1:length(symbol)) {
      cycIDList <- TrySymConv(symbol[i], speCycID)
      if (cycIDList$cycID != '0') {
        break
      } else {}
    }
  }
  else if (type == 'protein') {
    # transfer KEGG ID to unipro ID
    standKEGGID <- paste(speKEGGID, KEGGID, sep = ':')
    uniproID <- convKEGG('uniprot', standKEGGID, convertType = 'identity', n = 1)
    uniproID <- uniproID[1, 2]
    
    # 'up:Q8DWN9' --> 'Q8DWN9'
    uniproID <- sapply(strsplit(uniproID, split = ':', fixed = TRUE), '[[', 2)
    url <- paste0('http://websvc.biocyc.org/', speCycID, '/foreignid?ids=Uniprot:', uniproID)
    cycIDMat <- webTable(url, ncol = 1)
    if (cycIDMat[2, 1] == 1) {
      cycID <- cycIDMat[3, 1]
    }
    else if (cycIDMat[2, 1] == 0) {
      cycID <- NULL
    }

    cycID <- testLen(cycID, '0', cycID)
    cycIDList <- list(cycID = cycID, url = url)
  }
  
  return(cycIDList)
}

##' Get node values.
##'
##' Get node values from the XML file.
##' @title Get all nodes value directly from the XML file.
##' @param xmlFile XML file.
##' @param nodePath The XPath of nodeset (one or mutiple nodes).
##' @return Nodeset value
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML getNodeSet xmlValue
##' @keywords internal
##' 
xmlNodeVal <- function(xmlFile, nodePath){

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
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML getNodeSet xmlValue
##' @keywords internal
##' 
xmlNodeAttr <- function(xmlFile, nodePath, attrName){

  nodeSet <- getNodeSet(xmlFile, nodePath)
  nodeAttr <- sapply(nodeSet, xmlGetAttr, name = attrName)

  return(nodeAttr)

}


##' Test if the input vector's length is 0
##'
##' To test the length of input vector, if the length is 0, return 'trueVal', else return 'falseval'
##' @title Test length is 0 or not 
##' @param inputVal vector
##' @param trueVal return this value, if the length of 'inputVal' is 0.
##' @param falseVal return this value, if the length of 'falseVal' is not 0.
##' @return 'trueVal' or 'falseVal'
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @keywords internal
##' 
testLen <- function(inputVal, trueVal, falseVal) {

  if (length(inputVal) == 0) {
    return(trueVal)
  } else {
    return(falseVal)
  }
}


## Transfer KEGG ID to BioCyc ID.
##
## Tranfser KEGG ID to BioCyc ID by the route, whole species gene list --> unit gene information --> match KEGG ID
## @title Transfer genes ID from KEGG to BioCyc
## @param KEGGID A vector of KEGG gene ID.
## @param speID Species BioCyc ID.
## @param n The number of CPUs or processors, and the default value is 4.
## @param type Only support genes
##  @return A vector of BioCycID
## # EG10098 EG10101
## # "b3734" "b3732"
## # @examples transGeneIDKEGG2Cyc(c('b3734', 'b3732'), 'ECOLI')
## # @author Yulong Niu \email{niuylscu@@gmail.com}
##  @importFrom doMC registerDoMC
##  @importFrom foreach foreach
##  @export
##
## transGeneIDKEGG2Cyc <- function(KEGGID, speID, n = 4, type = 'genes'){

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

