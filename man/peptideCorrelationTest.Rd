\name{peptideCorrelationTest}
\alias{peptideCorrelationTest}
\title{peptideCorrelationTest}
\description{Take the data, and a mapping of peptides to
proteins, and compute the correlation between
peptides linked to a given protein. Compare
to random.}
\usage{peptideCorrelationTest(dat, pepinfo, pepCol, protCol)}
\arguments{
  \item{dat}{The data with samples as rows and peptides as columns}
  \item{pepinfo}{The mapping of peptides to proteins as a data frame}
  \item{pepCol}{The column name of peptide info table containing peptide IDs}
  \item{protCol}{The column name of pepinfo info table containing protein IDs}
}

\value{return a list of protein correlations and random correlations}

\author{David L Gibbs}





