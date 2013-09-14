\name{Accessors}
\Rdversion{1.1}
\docType{methods}
\alias{TOM}
\alias{TOM,proconaNet-method}
\alias{TOM<-}
\alias{TOM<-,proconaNet-method}
\alias{proconaVersion}
\alias{proconaVersion,proconaNet-method}
\alias{proconaVersion<-}
\alias{proconaVersion<-,proconaNet-method}
\alias{networkName}
\alias{networkName,proconaNet-method}
\alias{networkName<-}
\alias{networkName<-,proconaNet-method}
\alias{samples}
\alias{samples,proconaNet-methods}
\alias{samples<-}
\alias{samples<-,proconaNet-methods}
\alias{adj}
\alias{adj,proconaNet-method}
\alias{adj<-}
\alias{adj,proconaNet-method}
\alias{peptides}
\alias{peptides,proconaNet-method}
\alias{peptides<-}
\alias{peptides<-,proconaNet-method}
\alias{pepTree}
\alias{pepTree,proconaNet-method}
\alias{pepTree<-}
\alias{pepTree<-,proconaNet-method}
\alias{dynamicColors}
\alias{dynamicColors,proconaNet-method}
\alias{dynamicColors<-}
\alias{dynamicColors<-,proconaNet-method}
\alias{mergedColors}
\alias{mergedColors,proconaNet-method}
\alias{mergedColors<-}
\alias{mergedColors<-,proconaNet-method}
\alias{MEs}
\alias{MEs,proconaNet-method}
\alias{MEs<-}
\alias{MEs<-,proconaNet-method}
\alias{colorOrder}
\alias{colorOrder,proconaNet-method}
\alias{colorOrder<-}
\alias{colorOrder<-,proconaNet-method}
\alias{samples}
\alias{samples,proconaNet-method}
\alias{samples<-}
\alias{samples<-,proconaNet-method}
\alias{mergedMEs}
\alias{mergedMEs,proconaNet-method}
\alias{mergedMEs<-}
\alias{mergedMEs<-,proconaNet-method}
\alias{colorOrder}
\alias{power}
\alias{networkType}
\alias{permtest}
\alias{pepnet-class}
\alias{pepnet}
\alias{networkPower}
\alias{colorOrder,proconaNet-method}
\alias{power,proconaNet-method}
\alias{networkType,proconaNet-method}
\alias{permtest,proconaNet-method}
\alias{networkPower,proconaNet-method}
\alias{colorOrder<-}
\alias{power<-}
\alias{networkType<-}
\alias{permtest<-}
\alias{networkPower<-}
\alias{colorOrder<-,proconaNet-method}
\alias{power<-,proconaNet-method}
\alias{networkType<-,proconaNet-method}
\alias{permtest<-,proconaNet-method}
\alias{networkPower<-,proconaNet-method}


\title{Accessors for the proconaNet S4 class}
\description{Accessor functions allow access to the object data.}
\section{Methods}{
  \describe{
    \item{\code{TOM}:}{The topological overlap matrix or TOM. \code{"matrix"}}
    \item{\code{TOM<-}:}{The topological overlap matrix or TOM. \code{"matrix"}}

    \item{\code{networkName}:}{
      A name describing the data or experiment used to build the
      network. \code{"character"}
    }
    \item{\code{samples}:}{The names of samples used in building the network. \code{"character"}}
    \item{\code{adj}:}{The adjacency matrix. \code{"matrix"}}
    \item{\code{peptides}:}{The names of peptides used in the network,
      also the node names. \code{"character"}}
    \item{\code{pepTree}:}{The network dendrogram. \code{"hclust"}}
    \item{\code{dynamicColors}:}{The module labels on each node (or peptide). \code{"numeric"}}
    \item{\code{MEs}:}{The module eigenvectors (or eigen-peptides). \code{"data.frame"}}
    \item{\code{mergedMEs}:}{The module eigenvectors after merging
      similar modules. \code{"data.frame"}}
    \item{\code{mergedColors}:}{The module labels after merging similar modules. \code{"numeric"}}
    \item{\code{colorOrder}:}{Modules are ordered by size, these labels
      correspond to that order. \code{"character"}}
    \item{\code{power}:}{The soft thresholding power used in scaling the
      adjacency matrix. \code{"numeric"}}
    \item{\code{networkType}:}{Either a signed or unsigned network
      regarding the method used in computing the initial correlations
      between nodes. \code{"character"}}
    \item{\code{permtest}:}{The results of the permutation test on
      signficance of topological overlap within modules. \code{"matrix"}}
    \item{\code{proconaVersion}:}{
      Returns the version number of the software that built the object. \code{"character"}
    }
  }
}
\author{David L Gibbs}
\examples{
data(ProCoNA_Data)
tomMatrix <- TOM(net1)
}


