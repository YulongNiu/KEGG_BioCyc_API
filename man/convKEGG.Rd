% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/getKEGG.R
\name{convKEGG}
\alias{convKEGG}
\title{KEGG convert function}
\usage{
convKEGG(targetDB, sourceEntry, convertType = "database", n = 4)
}
\arguments{
\item{targetDB}{"targetDB" and "convertType" are set correspondingly.
"convertType" --> "database"
For gene convert from KEGG to outside databases.
"sourceEntry" --> KEGG organism code, or T number.
"targetDB" --> "ncbi-gi", "ncbi-geneid", or "uniprot".
For gene convert from outside databases to KEGG.
"sourceEntry" --> "ncbi-gi", "ncbi-geneid", or "uniprot".
"targetDB" --> KEGG organism code, or T number.
For chemical substance convert from KEGG to outside databases.
"sourceEntry" --> "drug", "compound", or "glycan".
"targetDB" --> "pubchem", or "chebi".
For chemical substance convert from outside databases to KEGG.
"sourceEntry" --> "pubchem", or "chebi".
"targetDB" --> "drug", "compound", or "glycan".


"convertType" --> "identity"
For gene convert
"sourceEntry" --> KEGG organism code, T number, "genes", "ncbi-gi", "ncbi-geneid", or "uniprot". "genes" is set to convert outside identities to KEGG when the organism code is not known.
"targetDB" --> A vector (length can be bigger than 1) of identities.
For chemical substance convert
"sourceEntry" --> "drug", "compound", "glycan", "pubchem", or "chebi".
"targetDB" --> A vector (length can be bigger than 1) of identities.}

\item{sourceEntry}{see "targetDB"}

\item{convertType}{set to be "database" or "identity".}

\item{n}{The number of CPUs or processors, and the default value is 4.}
}
\value{
A matrix that the first column is "targetDB"
}
\description{
KEGG Database API - Convert IDs between KEGG databases and outside databases
}
\details{
Convert gene identifiers or chemical substance identifiers between KEGG databases and oursite outside databases. For gene identifiers, this API provides functions to convert IDs/databases between KEGG databases and ncbi-gi/ncbi-geneid/uniprot. For chemical substance identifiers, it converts IDs/databases between drug/compound/glycan and pubchem/chebi. This API doesn't convert IDs between outside database. For example, the convert between pubchem and chebi IDs is not allowed. Try to use mutiple CPUs when conver large number of IDs.

The IDs and database convert is controlled by the argument "convertType".
}
\examples{
# convert database from KEGG to outside databases.
convKEGG('ncbi-geneid', 'eco')
convKEGG('pubchem', 'drug')

\dontrun{
# convert database from outside databases to KEGG.
convKEGG('smu', 'uniprot')
convKEGG('glycan', 'chebi')}

# convert identities from KEGG to outside database.
# mutiple organism convert.
convKEGG('ncbi-gi', c('hsa:10458', 'ece:Z5100'), convertType = 'identity', n = 2)
convKEGG('pubchem', 'cpd:C00004', convertType = 'identity', n = 2)

# convert identities from outside databases to KEGG.
# the organism code is unknown.
convKEGG('genes', 'ncbi-geneid:3113320', convertType = 'identity', n = 2)
convKEGG('genes', 'ncbi-gi:54293358', convertType = 'identity', n = 2)
}
\author{
Yulong Niu \email{niuylscu@gmail.com}
}
\references{
\url{http://www.kegg.jp/kegg/rest/keggapi.html}
}

