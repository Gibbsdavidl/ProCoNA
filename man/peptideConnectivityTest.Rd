\name{peptideConnectivityTest}
\alias{peptideConnectivityTest}
\title{peptideConnectivityTest}
\description{This function will compare the connectivity between
peptides linked to a given protein, against a randomly
drawn, similarly sized, selection of peptides.
The hypothesis is that peptides from a given protein
should be more connected than random.}
\usage{
peptideConnectivityTest(pnet,
pepInfo,
pepCol,
protCol)
}
\arguments{
  \item{pnet}{The peptide net object}
  \item{pepInfo}{The peptide information table, mapping peptides to proteins}
  \item{pepCol}{The string identifying the column in the pepInfo table with peptide ID}
  \item{protCol}{String identifying column in pepInfo with Protein ID.}
}

\value{Returns a list of the connected peptides and the random samples.}

\author{David L Gibbs}





