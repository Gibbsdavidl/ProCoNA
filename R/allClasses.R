

# Some code from the OneHandClapping package
# to enable the use of hclust in a S4 class.
# dummy class because 'hclust' is S3 class
setClass("hclust",contains="list")


##|| The pepcor class will hold all the results from peptide network calculations
##|| Version of procona that built the network
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
##|| the adjacency matrix is signed or unsigned
##|| the results of the permutation test

setClass("proconaNet",
         representation(
           proconaVersion = "character",
           networkName="character",
           samples="character",
           adj = "matrix",
           TOM = "matrix",
           peptides = "character", 
           pepTree = "hclust", 
           dynamicColors = "numeric", 
           MEs = "data.frame", 
           mergedMEs = "data.frame", 
           mergedColors = "numeric", 
           colorOrder = "character", 
           power= "numeric",
           networkType= "character",
           permtest="matrix" 
           ),
         prototype(
           proconaVersion="0",
           networkName="ProCoNA Network",
           samples=c("a","b","c"),
           adj=matrix(data=0, ncol=3, nrow=3),
           TOM=matrix(data=0, ncol=3, nrow=3),
           peptides=c("d","e","f"),
           pepTree= hclust(d=dist(matrix(c(1,1,1,1)))),
           dynamicColors=c(-1,-1,-1),
           MEs=data.frame(A=c(0,0,0), B=c(0,0,0), C=c(0,0,0)),
           mergedMEs=data.frame(A=c(0,0,0), B=c(0,0,0), C=c(0,0,0)),
           mergedColors=c(-1,-1,-1),
           colorOrder=c(""),
           power=1,
           networkType="signed",
           permtest=matrix()
           )
         #validity = check_proconaNet
         )
