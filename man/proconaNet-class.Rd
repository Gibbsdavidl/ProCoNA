\name{proconaNet-class}

\Rdversion{1.1}
\docType{class}
\alias{proconaNet-class}

\alias{proconaNet}
\alias{proconaVersion}
\alias{networkName}
\alias{samples}
\alias{adj}
\alias{TOM}
\alias{peptides}
\alias{pepTree}
\alias{dynamicColors}
\alias{mergedColors}
\alias{MEs}
\alias{colorOrder}
\alias{mergedMEs}
\alias{colorOrder}
\alias{power}
\alias{networkType}
\alias{permtest}
\alias{pepnet-class}
\alias{pepnet}
\alias{networkPower}


\title{proconaNet S4 class}
\description{NA}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new(proconaNet ...)}}
\section{Slots}{
  \describe{
    \item{\code{networkName}:}{
      A name describing the data or experiment used to build the
      network. \code{"character"}
    }
    \item{\code{samples}:}{The names of samples used in building the network. \code{"character"}}
    \item{\code{adj}:}{The adjacency matrix. \code{"matrix"}}
    \item{\code{TOM}:}{The topological overlap matrix or TOM. \code{"matrix"}}
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
\section{Methods}{
  \describe{
    \item{\code{modulePhenotypeCorrelations}:}{Correlating the module
      summaries to phenotypes.}
    \item{\code{moduleMemberCorrelations}:}{Correlations between modules.}
    \item{\code{compareNetworksWithFishersExactTestProcona}:}{Comparing
      two networks by module contents, pairwise between modules. }
    \item{\code{correlationWithPhenotypesHeatMap}:}{Heatmap of
      associations between modules and phenotypes.}
    \item{\code{MMvsPS}:}{Module members vs (a given) Phenotypes }
    \item{\code{MMvsPSallModules}:}{Module members vs all phenotypes.}
    \item{\code{toPermTest}:}{The permutation test for topological connectivity.}
    \item{\code{peptideConnectivityTest}:}{Test for connectivity among
      peptides of a similar protein.}
    \item{\code{peptideCorrelationTest}:}{Are peptides of a given protein
      more correlated than peptides drawn at random?}
    \item{\code{goStatTest}:}{GO term enrichment on modules.}
    \item{\code{keggStatTest}:}{Kegg pathway enrichment on modules.}
    \item{\code{ppiPermTest}:}{Protein-protein interaction enrichment on modules.}
  }
}

\author{David L Gibbs}





