% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/getNCBI.R
\name{getNCBITaxo}
\alias{getNCBITaxo}
\title{Get NCBI taxonomy information}
\usage{
getNCBITaxo(NCBITaxoIDs)
}
\arguments{
\item{NCBITaxoIDs}{A vector of NCBI taxonomy IDs.}
}
\value{
A list containing taxonomy information for each ID.
}
\description{
NCBI Database Additional API - Get NCBI taxonomy information from a given NCBI taxonomy IDs
}
\details{
Get NCBI taxonomy information.
}
\examples{
threeTax <- getNCBITaxo(c("9606", "511145", "797302"))
}
\author{
Yulong Niu \email{niuylscu@gmail.com}
}

