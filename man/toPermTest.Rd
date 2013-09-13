\name{toPermTest}
\alias{toPermTest}
\alias{toPermTest,proconaNet,numeric-method}
\title{toPermTest}
\description{Uses the procona network object,
the data with peptides as columns, samples in rows.
And the power that the net was built at
the number of permutations to do...
Modules are permuted and mean topological overlap
is recorded, constructing the null.  The number
of random permutations with mean TO greater than
observed provides the p-value.}
\usage{toPermTest(pnet, numPermutes)}
\arguments{
  \item{pnet}{ProCoNA network object.}
  \item{numPermutes}{The number of permutations to perform}
}

\value{returns the network obj with the perm test}
\examples{
data(ProCoNA_Data)
#net1 <- buildProconaNetwork("peptide network", peptideData, pow=12)
toPermTest(net1, 100)
}

\author{David L Gibbs}





