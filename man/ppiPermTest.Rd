\name{ppiPermTest}
\alias{ppiPermTest}
\title{ppiPermTest}

\description{Performs a permutation test for enrichment of PPI edges given a database. Peptides are selected from each module and mapped to potential protein parents in the mass tag database. We check if these proteins are found in the PPI network, and record any edges between them. This is compared to edges found using randomly selected proteins (taken from the mass tag database). A p-value is computed as the number of times the randomly sampled proteins incurred more edges than the observed proteins, divided by the number of iterations.}

\usage{ppiPermTest(pnet, pepdat, pepinfo, pepColName, pi_colName, pi_edges,  
   			   threshold, iterations)}

\arguments{
  \item{pnet}{procona network object}
  \item{pepdat}{the data matrix with peptides as columns.}
  \item{pepinfo}{Maps peptides to proteins ... same format as in ppiTable}
  \item{pepColName}{The column in pepinfo with peptide IDs... as in pepdat (the peptide data matrix)}
  \item{pi_colName}{The column in pepinfo that maps peptides to unit found in pi_edges}
  \item{pi_edges}{Must be two columns A-B ... sort out evidence levels (in vivo or in vitro) in advance}
  \item{threshold}{Minimum peptide correlation with module eigenvector.}
  \item{iterations}{Number of repititions}
}

\value{returns list of test results.}

\examples{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData, pow=12)
ppis <- data.frame(A=sample(masstagdb$Reference, 50), B=sample(masstagdb$Reference, 50))
ppiPermTest(net1, peptideData, masstagdb, "Mass_Tag_ID", "Reference", ppis, 0.33, 100)
}

\author{David L Gibbs}





