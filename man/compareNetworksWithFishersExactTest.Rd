\name{compareNetworksWithFishersExactTest}
\alias{compareNetworksWithFishersExactTest}
\title{compareNetworksWithFishersExactTest}
\description{This makes the heatmap of agreement betewen networks, module-wise,
by comparing each using Fisher's exact test where...
n == number of entities in the network
m == number of entities in intersection of two modules
d1 == number of entities in module A but not in module B
d2 == number of entities in module B but not in module A
2x2 matrix for the test is then:
m  d1
d2 n-d1-d2-m}
\usage{
compareNetworksWithFishersExactTest(net1,
net2,
colors1,
colors2,
title = "",
net1label = "",
net2label = "")
}
\arguments{
  \item{net1}{Net 1 entities, character vector}
  \item{net2}{Net 2 entiites, character vector}
  \item{colors1}{modules for net 1}
  \item{colors2}{modules for net 2}
  \item{title}{Plot title}
  \item{net1label}{xlabel}
  \item{net2label}{ylabel}
}

\value{returns fishers exact test -log pvalues and overlap matrix}

\author{David L Gibbs}





