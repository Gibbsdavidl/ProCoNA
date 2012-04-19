
##|| The pepcor class will hold all the results from peptide network calculations
##|| the name of the network
##|| the samples used to make the network
##|| the adjacency matrix
##|| the topological overlap matrix
##|| the names of the peptides, usually AMT IDs
##|| the original dendrogram
##|| The cuts from dynamic tree cut
##|| the module eigenvectors from dynamic cut
##|| module eigen-peptides
##|| the merged colors
##|| the color order
##|| the power used in the model
##|| the results of the permutation test

setClass("pepnet",
         representation(networkName="character",
                        samples="character",
                        adjacency = "matrix",
                        TOM = "matrix",
                        peptides = "character", 
                        pepTree = "hclust", 
                        dynamicColors = "numeric", 
                        MEs = "data.frame", 
                        mergedMEs = "data.frame", 
                        mergedColors = "numeric", 
                        colorOrder = "character", 
                        power= "numeric", 
                        permtest="matrix" 
                        )
        
         )


proconaPeptides <- function
### this function returns the names of the peptides in the pepnet object
(p  ##<< the procona obj
 ) {
  return(p@peptides)
### returns the vector of peptide IDs
}

#setMethod("peptides", signature(pnet="pepnet"),
#          function (pnet) peptides(pnet)
#          )

proconaName <- function
### this function returns the names of the peptides in the pepnet object
(p  ##<< the procona obj
 ) {
  return(p@proconaName)
### returns the vector of peptide IDs
}


neworkSamples <- function
### this function returns the names of the samples
(p  ##<< the procona obj
 ) {
  return(p@samples)
### returns the vector of sample IDs
}

proconaAdj <- function
### this function returns the adjacency matrix
(p  ##<< the procona obj
 ) {
  return(p@adjacency)
### returns the adjacency matrix
}


proconaTOM <- function
### this function returns the topological overlap matrix
(p  ##<< the procona obj
 ) {
  return(p@TOM)
### returns the TOM
}


proconaPepTree <- function
### this function returns the peptide dendrogram
(p  ##<< the procona obj
 ) {
  return(p@peptree)
### returns the peptide dendrogram
}


proconaMergedMEs <- function
### this function returns the merged module eigenvectors
(p  ##<< the procona obj
 ) {
  return(p@mergedMEs)
### returns the matrix of module eigenvectors in columns
}


proconaMergedColors <- function
### this function returns the vector of module identifiers
(p  ##<< the procona obj
 ) {
  return(p@mergedColors)
### returns the vector module names for each peptide
}



                                        #setMethod("mergedColors", signature(pnet="pepnet"),
#          function(pnet) mergedColors(pnet) 
#          )


proconaPower <- function
### this function returns the names of the peptides in the pepnet object
(p  ##<< the procona obj
 ) {
  return(p@power)
### returns the vector of peptide IDs
}


proconaPermtest <- function
### this function returns the names of the peptides in the pepnet object
(p  ##<< the procona obj
 ) {
  return(p@permtest)
### returns the vector of peptide IDs
}

