\name{peptideConnectivityTest}
\alias{peptideConnectivityTest}
\title{peptideConnectivityTest}

\description{This function will compare the connectivity between
peptides mapped to a given protein, against a randomly
drawn, similarly sized, selection of peptides.
The hypothesis is that peptides from a given protein
should be more connected than random.}

\usage{peptideConnectivityTest(pnet, pepInfo, pepCol, protCol, repsPerProt)}

\arguments{
  \item{pnet}{The peptide net object}
  \item{pepInfo}{The peptide information table, mapping peptides to proteins}
  \item{pepCol}{The string identifying the column in the pepInfo table with peptide ID}
  \item{protCol}{String identifying column in pepInfo with Protein ID.}
  \item{repsPerProt}{number of repetitions for the null}
}

\value{Returns a list of the connected peptides and the random samples.}

\examples{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData, pow=12)
p <- peptideConnectivityTest(net1, masstagdb, "Mass_Tag_ID", "Reference", 200)
}

\author{David L Gibbs}





