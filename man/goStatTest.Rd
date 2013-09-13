\name{goStatTest}
\alias{goStatTest}
\alias{goStatTest,proconaNet,numeric,data.frame,character,character,character,character,character,numeric,logical-method}
\title{goStatTest}

\description{Wrapper function to run the hyperGTest from package GOstats, after mapping each peptide to an entrez ID.}

\usage{goStatTest(pnet, module, pepinfo, pepColName, protColName, universe, 
    onto, annot, pvalue, cond)}

\arguments{
  \item{pnet}{Procona network object.}
  \item{module}{Module of interest (numeric)}
  \item{pepinfo}{The mass tag info, mapping peptides to proteins.}
  \item{pepColName}{Column name in mass tag info for peptides}
  \item{protColName}{Column name in mass tag info for proteins}
  \item{universe}{Table mapping protein IDs to entrez IDs}
  \item{onto}{The ontology catagory (bp etc)..}
  \item{annot}{The annotation database to use}
  \item{pvalue}{pvalue cutoff}
  \item{cond}{conditional parameter, see GOstats.}
}

\value{Returns the results of the hyper geometric test.}

\examples{
\dontrun{
data(ProCoNA_Data)
#net1 <- buildProconaNetwork("peptide network", peptideData, pow=12)
goStatTest(net1, 1, masstagdb, "Mass_Tag_ID", "Reference", universe, "BP", "org.Mm.eg.db", 0.005, FALSE)
}
}
\author{David L Gibbs}





