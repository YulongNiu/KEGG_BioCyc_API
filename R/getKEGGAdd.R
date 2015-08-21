##' Get the gene motif table
##'
##' Get the gene motif tables and additional information.
##' 
##' @title Get motif information from KEGG
##' @param geneID A single KEGG gene ID
##' @param hasAddInfo A logic element whether to show the additional information. The default value is "FALSE".
##' @return A list if "hasAddInfo" is TRUE, else a matrix. A special case retuns "NULL" instead of a matrix (other information is not effected), if there is no motif information.
##' @examples
##' getKEGGGeneMotif('brp:103873230')
##' 
##' # no motif information
##' getKEGGGeneMotif('hsa:4558', hasAddInfo = TRUE)
##' @importFrom XML readHTMLTable
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
getKEGGGeneMotif <- function(geneID, hasAddInfo = FALSE) {
  ## motif webpage URL
  url <- paste0('http://www.kegg.jp/ssdb-bin/ssdb_motif?kid=', geneID)

  ## read merged TableList
  readMotifList <- readHTMLTable(url)

  ## gene information
  Organism <- colnames(readMotifList[[1]])[2]

  ## GeneName
  GeneName <- geneID

  ## Definition
  Definition <- as.character(readMotifList[[1]][2, 2])

  ## motif Table
  motifTable <- readMotifList[[2]]

  if (!is.null(motifTable)) {
    motifTable <- apply(motifTable, 1:2, as.character)
  } else {}

  if (hasAddInfo) {
    motifList <- list(Organism = Organism,
                      GeneName = GeneName,
                      Definition = Definition,
                      motifTable = motifTable)
  } else {
    motifList <- motifTable
  }

  return(motifList)
}

