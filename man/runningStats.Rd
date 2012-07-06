\name{runningStats}
\alias{runningStats}
\title{runningStats}
\description{Computing the running mean and variance }
\usage{runningStats(newMat, runningMean, Mk1, Sk1, k)}
\arguments{
  \item{newMat}{The matrix from resampled data}
  \item{runningMean}{The running mean matrix}
  \item{Mk1}{Matrix used in calculation of mean}
  \item{Sk1}{Matrix used in calculation of sd}
  \item{k}{Current resampling iteration}
}

\value{returns the list of runningMean, runningSD, Mk, Sk}

\author{David L Gibbs}





