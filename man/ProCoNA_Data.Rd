\name{ProCoNA-Data}
\alias{masstagdb}
\alias{peptideData}
\alias{universe}
\alias{phenotypes}
\alias{net1}
\alias{net2}


\title{A simulated mass tag data base}
\description{
  The mass tag database, which would be used to identify peptides, simply
  maps peptide IDs to peptide sequences and protein matches.

  
  This simulated peptide dataset was generated using OpenMS's MSSimulator. A set
  of proteins was randomly sampled, and used to generate a likely set of observed
  peptides. Then data for a co-expression network was simulated with WGCNA's 
  simulation function, and columns were named with simulated peptides.

  
  The mass tag database, which would be used to identify peptides, simply
  maps peptide IDs to peptide sequences and protein matches. This represents 
  a mapping to Entrez IDs.

  
    The matrix annotates the biological samples according to
  ... phenotypic observations!


  The two network objects are included to avoid rebuilding them in the
  other man page examples.
}

\keyword{datasets}



