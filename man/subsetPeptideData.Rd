\name{subsetPeptideData}
\alias{subsetPeptideData}
\title{subsetPeptideData}
\description{Given a matrix of peptide data, omit columns with excess missing data, specified by NAs.}
\usage{subsetPeptideData(pepdat, numNAsAllowed = NULL, percentageNAsAllowed = 0.05)}
\arguments{
  \item{pepdat}{The peptide matrix, with peptides in columns and samples in rows.}
  \item{numNAsAllowed}{The maximum count of missing values for each peptide (counts NAs).}
  \item{percentageNAsAllowed}{The percentage of missing data allowed for each peptide over samples.}
}

\value{Returns a matrix.}
\examples{
data(ProCoNA_Data)
subsetPeptideData(peptideMatrix, percentageNAsAllowed=0.2)
}
\author{David L Gibbs}





