

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
           samples="",
           adj=matrix(),
           TOM=matrix(),
           peptides="",
           pepTree= hclust(d=dist(matrix(c(1,1,1,1)))),
           dynamicColors=c(0),
           MEs=data.frame(),
           mergedMEs=data.frame(),
           mergedColors=c(0),
           colorOrder="",
           power=1,
           networkType="",
           permtest=matrix()
           )
         #validity = check_proconaNet
         )
