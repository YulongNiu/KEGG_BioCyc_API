##' KEGG Database API - Get species KEGG/NCBI ID.
##'
##' Get the phylogenetic information of given species.
##' It supports both batch input and regular expression search.
##'
##' @title Get species list from KEGG.
##' @name getKEGGPhylo
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
##' getKEGGPhylo(c('hsa', 'eco'))
##' # search species whose names include 'Escherichia coli'
##' getKEGGPhylo('Escherichia coli', speType = 'regexpr')
##' # search species whose class is 'Mammals'
##' getKEGGPhylo('Mammals', speType = 'phylo')
##' # get whole KEGG species information table
##' getKEGGPhylo(whole = TRUE)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
getKEGGPhylo <- function(speList, speType = 'KEGG', whole = FALSE){

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


##' KEGG Database API - Get the KEGG orthology list.
##'
##' Get the KEGG orthology list by a given KEGG KO ID.
##' @title Get KEGG orthology.
##' @param KOID The KEGG orthology ID.
##' @return A character vector of KEGG gene IDs
##' @examples getKEGGKO('K02110')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
getKEGGKO <- function(KOID){

  # get KO webpage
  url <- paste('http://rest.kegg.jp/link/genes/', KOID, sep = '')

  # get KO list
  KOWebMat <- webTable(url, ncol = 2)

  geneIDs <- KOWebMat[, 2]
  ## # seperate the species and genes
  ## KOMat <- unlist(strsplit(KOWebMat[, 2], split = ':', fixed = TRUE))
  ## KOMat <- matrix(KOMat, ncol = 2, byrow = TRUE)
  ## colnames(KOMat) <- c('speID', 'geneID')

  return(geneIDs)

}


##' KEGG Database API - Get the whole pathway ID from KEGG database.
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

  ## # the transfer the pathname of 'path:hsa00010' to 'hsa00010'.
  ## pathID <- pathAnno[, 1]
  ## pathID <- sapply(strsplit(pathID, split = ':', fixed = TRUE), '[', 2)
  ## pathAnno[, 1] <- pathID

  colnames(pathAnno) <- c('pathID', 'Annotation')

  return(pathAnno)
}


##' KEGG Database API - Get the pathway and genes.
##'
##' Get the pathway and genes according to KEGG species ID.
##' @title List pathways and genes of a given KEGG species ID
##' @param KEGGspec The KEGG species ID. Only one species ID once input.
##' @return A List named with KEGG pathway IDs, and each element of the list contains the KEGG gene IDs.
##' @examples
##' \dontrun{
##' getKEGGPathGenes('hsa')}
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

##' KEGG Database API - Get the whole KEGG IDs from one species.
##'
##' Get the KEGG protein ID list and annotation.
##' @title Get whole KEGG IDs and annotation
##' @param KEGGspec KEGSS species org code or T number , for example 'hsa' or 'T01001'.
##' @return A matrix of KEGG IDs and annotation
##' @examples
##' # KEGG org cord
##' getProID('eco')
##'
##' # KEGG T number
##' getProID('T00007')
##'
##' # KEGG T number with empty elements
##' getProID('T10004')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
##' 
getProID <- function(KEGGspec){
  
  # get KEGG ID annotation list
  url <- paste('http://rest.kegg.jp/list/', KEGGspec, sep = '')

  # transfer webpage into a matrix
  speIDAnno <- webTable(url, ncol = 2)

  return(speIDAnno)
}



##' KEGG Database API - Get the nucleotide acid and amino acid sequences 
##'
##' Get the protein and gene sequences in fasta format. This function support mutiple querys.
##' @title Get protein and gene sequences
##' @param KEGGID A vector of KEGG IDs. Seqences from different species could be combined together.
##' @param seqType Choose nucleotide acid (ntseq) or amino acid (aaseq) seqences, and the default is amino acid sequences.
##' @param n The number of CPUs or processors, and the default value is 4.
##' @return A BStringSet 
##' @examples
##' # two amino acid seqences from different sepecies with 2 threads.
##' twoAASeqs <- getKEGGGeneSeq(c('mja:MJ_0011', 'hsa:10458'), n = 2)
##' \dontrun{
##' # export fasta format files
##' require('Biostrings')
##' writeXStringSet(twoAASeqs, 'twoAASeqs.fasta')}
##'
##' \dontrun{
##' # more examples
##' twoNTSeqs <- getKEGGGeneSeq(c('shy:SHJG_7159', 'shy:SHJG_7160'), 'ntseq')
##' mutilAASeqs <- getKEGGGeneSeq(c('eco:b0202', 'eco:b0203', 'eco:b0204',
##' 'eco:b0205', 'eco:b0206', 'eco:b0216', 'eco:b0244',
##' 'eco:b4626', 'eco:b3796', 'eco:b3797', 'eco:b3296',
##' 'eco:b3297'))}
##' 
##' \dontrun{
##' # get the whole E.coli genome protein seqences
##' ecoProIDs <- getProID('eco')
##' ecoGenomePro <- getKEGGGeneSeq(ecoProIDs[, 1])}
##' @importFrom RCurl getURL
##' @importFrom doMC registerDoMC
##' @importFrom foreach foreach %dopar%
##' @importFrom Biostrings BStringSet
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @references \url{http://www.kegg.jp/kegg/rest/keggapi.html}
##' @seealso getKEGGTIDGeneSeq
##' @export
##'
##' 
getKEGGGeneSeq <- function(KEGGID, seqType = 'aaseq', n = 4){

  # register mutiple cores
  registerDoMC(n)

  doTen <- function(tenWebSeq, seqTypeBio = seqType){
    # USE: a temprary function to deal with the return web "10 seqence", and it also for less ten or without coding sequence. The basic idea to split sequences is by the marker of '(A)' for amino acid sequences and '(N)' for nucleotide sequences.
    # INPUT: return results from function getURL()
    # OUTPUT: BStingSet or NA (for the case one gene has no cording protein sequence)
    if (nchar(tenWebSeq) != 0){
      splitTen <- unlist(strsplit(tenWebSeq, split = '\n', fixed = TRUE))
      if (seqTypeBio == 'aaseq') {
        namePoint <- which(grepl(' \\(A\\)', splitTen))
      }
      else if (seqTypeBio == 'ntseq') {
        namePoint <- which(grepl(' \\(N\\)', splitTen))
      }
      # sequence start point
      startPoint <- namePoint + 1
      endPoint <- c(namePoint[-1] - 1, length(splitTen))
    } else {
      # gene without coding sequence
      return(NULL)
    }

    tenSeqBS <- foreach(i = 1:length(namePoint), .combine = append) %dopar% {
      proSeqName <- splitTen[namePoint[i]]
      proSeqName <- substring(proSeqName, 2)
      proSeq <- paste(splitTen[startPoint[i] : endPoint[i]], collapse = '')
      # to BStringSet
      proSeqBS <- BStringSet(proSeq)
      names(proSeqBS) <- proSeqName
      return(proSeqBS)
    }

    return(tenSeqBS)
  }

  # cut the input 'KEGGID' into 10
  cutMat <- CutSeqEqu(length(KEGGID), 10)

  # deal with ten sequences each time
  print(paste('The input ID length is ', length(KEGGID), '.', sep = ''))
  proSeq <- foreach(i = cutMat[1, ], j = cutMat[2, ], .combine = append) %dopar% {
    print(paste('Get ', j , ' proteins.'))
    if (i == j){
      # only one input KEGG ID
      linkKEGGPro <- paste('http://rest.kegg.jp/get/', KEGGID[i], '/', seqType, sep = '')
    } else {
      mergeID <- paste(KEGGID[i:j], collapse = "+")
      linkKEGGPro <- paste('http://rest.kegg.jp/get/', mergeID, '/', seqType, sep = '')
    }
    webProSeq <- getURL(linkKEGGPro)

    doTen(webProSeq)
  }

  return(proSeq)

}



##' KEGG Database API - Convert IDs between KEGG databases and outside databases
##'
##' Convert gene identifiers or chemical substance identifiers between KEGG databases and oursite outside databases. For gene identifiers, this API provides functions to convert IDs/databases between KEGG databases and ncbi-gi/ncbi-geneid/uniprot. For chemical substance identifiers, it converts IDs/databases between drug/compound/glycan and pubchem/chebi. This API doesn't convert IDs between outside database. For example, the convert between pubchem and chebi IDs is not allowed. Try to use mutiple CPUs when conver large number of IDs. 
##'
##' The IDs and database convert is controlled by the argument "convertType".
##' @title KEGG convert function
##' @param targetDB "targetDB" and "convertType" are set correspondingly.
##' "convertType" --> "database"
##' For gene convert from KEGG to outside databases.
##' "sourceEntry" --> KEGG organism code, or T number.
##' "targetDB" --> "ncbi-gi", "ncbi-geneid", or "uniprot".
##' For gene convert from outside databases to KEGG.
##' "sourceEntry" --> "ncbi-gi", "ncbi-geneid", or "uniprot".
##' "targetDB" --> KEGG organism code, or T number.
##' For chemical substance convert from KEGG to outside databases.
##' "sourceEntry" --> "drug", "compound", or "glycan".
##' "targetDB" --> "pubchem", or "chebi".
##' For chemical substance convert from outside databases to KEGG.
##' "sourceEntry" --> "pubchem", or "chebi".
##' "targetDB" --> "drug", "compound", or "glycan".
##' 
##' 
##' "convertType" --> "identity"
##' For gene convert
##' "sourceEntry" --> KEGG organism code, T number, "genes", "ncbi-gi", "ncbi-geneid", or "uniprot". "genes" is set to convert outside identities to KEGG when the organism code is not known.
##' "targetDB" --> A vector (length can be bigger than 1) of identities. 
##' For chemical substance convert
##' "sourceEntry" --> "drug", "compound", "glycan", "pubchem", or "chebi".
##' "targetDB" --> A vector (length can be bigger than 1) of identities.
##' 
##' @param sourceEntry see "targetDB"
##' @param convertType set to be "database" or "identity".
##' @inheritParams getKEGGGeneSeq
##' @return A matrix that the first column is "targetDB"
##' @examples
##' # convert database from KEGG to outside databases.
##' convKEGG('ncbi-geneid', 'eco')
##' convKEGG('pubchem', 'drug')
##' 
##' \dontrun{
##' # convert database from outside databases to KEGG.
##' convKEGG('smu', 'uniprot')
##' convKEGG('glycan', 'chebi')}
##'
##' # convert identities from KEGG to outside database.
##' # mutiple organism convert.
##' convKEGG('ncbi-gi', c('hsa:10458', 'ece:Z5100'), convertType = 'identity', n = 2)
##' convKEGG('pubchem', 'cpd:C00004', convertType = 'identity', n = 2)
##'
##' # convert identities from outside databases to KEGG.
##' # the organism code is unknown.
##' convKEGG('genes', 'ncbi-geneid:3113320', convertType = 'identity', n = 2)
##' convKEGG('genes', 'ncbi-gi:54293358', convertType = 'identity', n = 2)
##' @importFrom foreach foreach %dopar%
##' @importFrom doMC registerDoMC
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @references \url{http://www.kegg.jp/kegg/rest/keggapi.html}
##' @export
##' 
convKEGG <- function(targetDB, sourceEntry, convertType = 'database', n = 4) {

  if (convertType == 'identity') {
    
    registerDoMC(n)
    
    ## cut the every 10 sourceEntry
    cutMat <- CutSeqEqu(length(sourceEntry), 10)
    sourceSeq <- apply(cutMat, 2, function(x) {
      eachSeq <- paste(sourceEntry[x[1] : x[2]], collapse = '+')
      return(eachSeq)
    })

    urlSeq <- paste0('http://rest.kegg.jp/conv/', targetDB, '/', sourceSeq)
    
    convRes <- foreach(i = 1:length(urlSeq), .combine = rbind) %dopar% {
      convRes <- webTable(urlSeq[i], ncol = 2)
      ## deal with NA
      ## rbind(NULL, NULL) is NULL
      if (is.na(convRes[1, 1])) {
        convRes <- NULL
      } else {}
      
      return(convRes)
    }

    ## remove no result
    if (is.null(convRes)) {
      convRes <- matrix(rep(NA, 2), nrow = 1)
    } else {}
    
  } else {
    url <- paste0('http://rest.kegg.jp/conv/', targetDB, '/', sourceEntry)
    convRes <- webTable(url, ncol = 2)
  }

  return(convRes)
}


##' Get a R matrix object if the weblink returned as a matrix.
##'
##' If the web return a matrix, use this function to extract it as a R matrix object. An empty element is used to represent a "NA" in KEGG. So the empty web element is set to "NA" according the "ncol" parameter.
##' @title Get R matrix from weblink
##' @param url The weblink.
##' @param ncol The column number of the matrix.
##' @return A R matrix
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl getURL
##' @keywords internal
##' 
webTable <- function(url, ncol) {

  webPage <-getURL(url)

  ## transfer webpage into a matrix
  webMat <- unlist(strsplit(webPage, split = '\n', fixed = TRUE))
  webMat <- sapply(webMat, strsplit, split = '\t', fixed = TRUE)

  ## transfer empty web elements to NAs
  webMat <- sapply(webMat, function(x) {
    lenSub <- ncol - length(x)

    if (lenSub > 0) {
      x <- c(x, rep(NA, lenSub))
    } else {}

    return(x)
  })

  webMat <- matrix(unlist(webMat), ncol = ncol, byrow = TRUE)

  return(webMat)
}


##' Cut vectors with invervals
##'
##' CutSeq() is used to cut a vector with different invervals. CutSeqEqu() is used to cut a vector with same invervals.
##' @title Cut vectors
##' @param cutSeq The inverals vector. The length of cutSeq could be more than 1, and 0 will be automatically excluded. 
##' @return A cut matrix in which the first row is the start point and second row is the end point.
##' @examples
##' # with one interval
##' CutSeq(10)
##' # with multiple interval
##' CutSeq(c(2, 3, 5))
##' # exclude 0
##' @rdname CutSeqInterval
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
##' 
CutSeq <- function(cutSeq){

  # remove 0, because we cannot cut a sequence by the internal of 0.
  cutSeq <- cutSeq[cutSeq != 0]
  vecCutseq <- length(cutSeq)

  if (vecCutseq == 1) {
    headCut <- 1
    endCut <- cutSeq
  } else {
    # loopCutSeq is the circle of vecCutseq
    loopCutSeq <- list()
    for(i in 1:vecCutseq) {
      loopCutSeq[[i]] <- cutSeq[1:i]
    }

    loopSumCutSeq <-  sapply(loopCutSeq, sum)

    # the head and tail sequence
    headCut <- c(1,loopSumCutSeq[1:(vecCutseq-1)]+1)
    endCut <- loopSumCutSeq
  }

  cutMat <- matrix(c(headCut, endCut), 2, byrow=TRUE)

  return(cutMat)

}


##' @param vecLen The length of vector used to cut.
##' @param equNum The equal internal.
##' @examples
##' # equal interval is the same as the length of vector
##' CutSeqEqu(10, equNum = 10)
##' CutSeqEqu(21, equNum = 10)
##' # euqal interval is larger than the length of vector
##' CutSeqEqu(10, equNum = 20)
##' @rdname CutSeqInterval
##' @export
##'
##' 
CutSeqEqu <- function(vecLen, equNum){
  
  if (equNum > vecLen){
    # the internal is bigger than the length of vecLen. So we use the full vecLen.
    cutMat <- matrix(c(1, vecLen))
  } else {
    timeNum <- vecLen %/% equNum
    remainer <- vecLen %% equNum
    cutSeq <- c(rep(equNum, timeNum), remainer)
    cutMat <- CutSeq(cutSeq)
  }

  return(cutMat)
}

