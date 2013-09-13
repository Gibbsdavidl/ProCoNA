\name{corBootstrap}
\alias{corBootstrap}
\title{corBootstrap}

\description{Boostraps a correlation matrix. In order to bootstrap a large correlation matrix,
several thousand samplings may be necessary. To avoid storing thousands of matrices, a running 
mean is kept for each pairwise correlation. In addition, a running standard deviation is computed
so that for each pairwise correlation, we can estimate the distribution of values across resamplings. 
After each resampling, a new correlation matrix is computed. A difference is taken between this new
matrix and the running mean. If all differences are less than the specified threshold, then the
bootstrapped matrix has converged to a final state. 
}

\usage{corBootstrap(dataMatrix, networkType = "signed", threshold = 1e-04, 
    tmpSaveFile = TRUE)}

\arguments{
  \item{dataMatrix}{Matrix with samples in rows and peptides (or other data type) in columns.}
  \item{networkType}{Whether the sign is considered in constructing adjacency and TOM}
  \item{threshold}{Maximum difference allowed between running mean boostrap correlation matrix, and new resampled cor matrix.  Defines how soon we consider the bootstrap to have converged. }
  \item{tmpSaveFile}{Should temporary saves be done?}
}

\value{Returns a list of the bootstrapped matrix, standard deviation matrix, and the number of resamplings done.}

\examples{
data(ProCoNA_Data)
x <- peptideData[,1:10]
y <- corBootstrap(dataMatrix=x, networkType="unsigned", threshold=0.1, tmpSaveFile=FALSE)
}
\author{David L Gibbs}





