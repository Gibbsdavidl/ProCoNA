# Building the network by pieces...


newProconaObj <- function
### This function returns a peptide co-expression network object.
(networkName,             ##<< Name of this network
 pepdat,                  ##<< This variable is the data set with rows as samples and cols as peptides
 pow=NULL,                ##<< The scaling power, NULL if unknown
 signed=F,                ##<< Whether the sign is considered in constructing adjacency and TOM
 scaleFreeThreshold=0.8,  ##<< The threshold for fitting to scale-free topology.. will use closest power.
 deepSplit = 2,           ##<< Course grain control of module size
 minModuleSize=30,        ##<< The minimum module size allowed
 mergeThreshold=0.1,      ##<< Below this threshold, modules are merged.
 clusterType="average",   ##<< Clustering option
 performTOPermtest=TRUE,  ##<< Performs permutation testing on modules
 toPermTestPermutes=10000 ##<< Number of permutations to do.
 ){

  require(WGCNA)
  pnet = new("pepnet")
  pnet@networkName=networkName
  if (is.null(pow)) {
    print("Computing soft threshold power")
    softThreshold <- pickSoftThreshold(pepdat, powerVector=1:20, RsquaredCut=scaleFreeThreshold)
    pnet@power=softThreshold$powerEstimate
  } else {
    pnet@power=pow
  }
  cat("Using power: ", pnet@power, "\n")
  pnet@peptides=colnames(pepdat)
  print("Computing adjacency")
  pnet@adjacency <- adjacency(datExpr=pepdat, power=pnet@power, type=signed)
  print("Computing TOM")
  pnet@TOM = TOMsimilarity(pnet@adjacency, TOMType=signed);
  dissTOM = 1-pnet@TOM
  print("Clustering")
  pnet@geneTree = flashClust(as.dist(dissTOM), method = clusterType);
  pnet@dynamicColors = cutreeDynamic(dendro = pnet@geneTree, distM = dissTOM,
    deepSplit = 2, pamRespectsDendro = FALSE,
    minClusterSize = minModuleSize);

  print(table(pnet@dynamicColors))
                                        #merging modules
  print("Merging modules")
  MEDissThres = mergeThreshold
  MEList = moduleEigengenes(pepdat, colors = pnet@dynamicColors)
  pnet@MEs = MEList$eigengenes
                                        # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(pnet@MEs);
                                        # Cluster module eigengenes
  METree = flashClust(as.dist(MEDiss), method = clusterType );
                                        # Call an automatic merging function
  merge = mergeCloseModules(pepdat, pnet@dynamicColors,
    cutHeight = MEDissThres, verbose = 4, relabel=T)
                                        # The merged module colors
  pnet@mergedColors = merge$colors;
                                        # Eigengenes of the new merged modules:
  pnet@mergedMEs = merge$newMEs;  
                                        # Construct numerical labels corresponding to the colors
  pnet@colorOrder = c("grey", standardColors(50));
                                        # Permutation test for module topological overlap
  pnet <- toPermTest(pnet, pepdat, toPermTestPermutes)
  
  print("DONE!")
  return(pnet)
  ### returns the procona network object
}
  


newRobustProconaObj <- function
### This function returns a peptide co-expression network object.
(networkName,             ##<< Name of this network
 pepdat,                  ##<< This variable is the data set with rows as samples and cols as peptides
 pow=NULL,                ##<< The scaling power, NULL if unknown
 signed="unsigned",       ##<< Whether the sign is considered in constructing adjacency and TOM
 scaleFreeThreshold=0.8,  ##<< The threshold for fitting to scale-free topology.. will use closest power.
 deepSplit = 2,           ##<< Course grain control of module size
 minModuleSize=30,        ##<< The minimum module size allowed
 mergeThreshold=0.2,      ##<< Below this threshold, modules are merged.
 clusterType="average",   ##<< Clustering option
 performTOPermtest=TRUE,  ##<< Performs permutation testing on modules
 toPermTestPermutes=10000, ##<< Number of permutations to do.
 cores=1                   ##<< Number of cores to use in tukey computation
 ){

  require(WGCNA)
  require(fastukeycor)
  pnet = new("pepnet")
  pnet@networkName=networkName
  if (is.null(pow)) {
    print("Computing soft threshold power")
    softThreshold <- pickSoftThreshold(pepdat, powerVector=1:20, RsquaredCut=scaleFreeThreshold)
    pnet@power=softThreshold$powerEstimate
  } else {
    pnet@power=pow
  }
  cat("Using power: ", pnet@power, "\n")
  pnet@peptides=colnames(pepdat)
  print("Computing adjacency")
  pnet@adjacency <- tukeyAdjacency(pepdat=pepdat, power=pnet@power, type=signed, cores=cores)
  print("Computing TOM")
  pnet@TOM = TOMsimilarity(pnet@adjacency, TOMType=signed);
  dissTOM = 1-pnet@TOM
  print("Clustering")
  pnet@geneTree = flashClust(as.dist(dissTOM), method = clusterType);
  pnet@dynamicColors = cutreeDynamic(dendro = pnet@geneTree, distM = dissTOM,
    deepSplit = 2, pamRespectsDendro = FALSE,
    minClusterSize = minModuleSize);

  print(table(pnet@dynamicColors))
                                        #merging modules
  print("Merging modules")
  MEDissThres = mergeThreshold
  MEList = moduleEigengenes(pepdat, colors = pnet@dynamicColors)
  pnet@MEs = MEList$eigengenes
                                        # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(pnet@MEs);
                                        # Cluster module eigengenes
  METree = flashClust(as.dist(MEDiss), method = clusterType );
                                        # Call an automatic merging function
  merge = mergeCloseModules(pepdat, pnet@dynamicColors,
    cutHeight = MEDissThres, verbose = 4, relabel=T)
                                        # The merged module colors
  pnet@mergedColors = merge$colors;
                                        # Eigengenes of the new merged modules:
  pnet@mergedMEs = merge$newMEs;  
                                        # Construct numerical labels corresponding to the colors
  pnet@colorOrder = c("grey", standardColors(50));
                                        # Permutation test for module topological overlap
  pnet <- toPermTest(pnet, pepdat, toPermTestPermutes)
  
  print("DONE!")
  return(pnet)
  ### returns the procona network object
}
  

