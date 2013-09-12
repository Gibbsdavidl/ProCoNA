\name{peptideData}
\alias{peptideData}
\title{Simulated data matrix with peptides in columns and samples in rows.}
\description{
  This simulated peptide dataset was generated using OpenMS's MSSimulator. A set
  of proteins was randomly sampled, and used to generate a likely set of observed
  peptides. Then data for a co-expression network was simulated with WGCNA's 
  simulation function, and columns were named with simulated peptides.

}
\usage{data(ProCoNA_Data)}
\format{A matrix with peptide columns and samples in rows.}
\keyword{datasets}
