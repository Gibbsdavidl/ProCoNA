
##|| The pepcor class will hold all the results from peptide network calculations
##|| the name of the network
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
                        adjacency = "matrix",
                        TOM = "matrix",
                        peptides = "character", 
                        geneTree = "hclust", 
                        dynamicColors = "numeric", 
                        MEs = "data.frame", 
                        mergedMEs = "data.frame", 
                        mergedColors = "numeric", 
                        colorOrder = "character", 
                        power= "numeric", 
                        permtest="matrix" 
                        )
        
         )


peptides <- function
### this function returns the names of the peptides in the pepnet object
(p  ##<< the procona obj
 ) {
  return(p@peptides)
### returns the vector of peptide IDs
}

#setMethod("peptides", signature(pnet="pepnet"),
#          function (pnet) peptides(pnet)
#          )

mergedColors <- function
### this function returns the names of the peptides in the pepnet object
(p  ##<< the procona obj
 ) {
  return(p@mergedColors)
### returns the vector of peptide IDs
}

#setMethod("mergedColors", signature(pnet="pepnet"),
#          function(pnet) mergedColors(pnet) 
#          )
