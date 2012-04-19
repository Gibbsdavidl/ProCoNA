\name{goStatTest}
\alias{goStatTest}
\title{goStatTest}
\description{Wrapper function to run the hyperGTest from package GOstats}
\usage{
goStatTest(pnet,
module,
pepinfo,
pepColName,
protColName,
universe,
onto,
annot,
pvalue,
cond)
}
\arguments{
  \item{pnet}{the procona network}
  \item{module}{module of interest}
  \item{pepinfo}{the mass tag info}
  \item{pepColName}{column in mass tag info for peptides}
  \item{protColName}{column in mass tag info for proteins}
  \item{universe}{the universe to consider, list of proteins}
  \item{onto}{the ontology catagory (bp etc)..}
  \item{annot}{the annotation database to use}
  \item{pvalue}{pvalue cutoff}
  \item{cond}{conditional parameter, see GOstats.}
}

\value{returns the results of the hyper geometric test.}

\author{David L Gibbs}





