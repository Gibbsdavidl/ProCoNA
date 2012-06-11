\name{toPermTest}
\alias{toPermTest}
\title{toPermTest}
\description{Uses the procona network object,
the data with peptides as columns, samples in rows.
And the power that the net was built at
the number of permutations to do...
Modules are permuted and mean topological overlap
is recorded, constructing the null.  The number
of random permutations with mean TO greater than
observed provides the p-value.}
\usage{toPermTest(pnet, dat, numPermutes = 100)}
\arguments{
  \item{pnet}{to get the colors}
  \item{dat}{the data used to build the model}
  \item{numPermutes}{The number of permutations to perform}
}

\value{returns the network obj with the perm test}

\author{David L Gibbs}





