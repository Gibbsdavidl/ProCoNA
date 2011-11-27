\name{tukeyAdjacency}
\alias{tukeyAdjacency}
\title{tukeyAdjacency}
\description{Adjacency matrix computed with Tukey's bi-weight correlation.}
\usage{
tukeyAdjacency(pepdat,
type = "unsigned",
power = 1,
cores = 1,
median = T,
full.init = T)
}
\arguments{
  \item{pepdat}{peptide data, rows = samples, cols = peptides}
  \item{type}{type of adjacency matrix, unsigned, signed, positiveOnly}
  \item{power}{scaling power (beta)}
  \item{cores}{number of cores to use}
  \item{median}{median or mcv init?}
  \item{full.init}{init from full data? or pairwise?}
}

\value{returns a scaled adjacency matrix}

\author{David L Gibbs}





