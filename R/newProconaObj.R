

newProconaObj <- function
### This function returns a peptide co-expression network object.
(networkName="procona",   ##<< Name of this network
 pepdat=NULL,             ##<< This variable is the data set with rows as samples and cols as peptides
 pow=NULL,                ##<< The scaling power, NULL if unknown
 powMax=20,               ##<< The maximum power to be searched.
 signed="unsigned",       ##<< Whether the sign is considered in constructing adjacency and TOM
 pearson=T,               ##<< use Pearson's cor or the robust bi-weight correlation
 scaleFreeThreshold=0.8,  ##<< The threshold for fitting to scale-free topology.. will use closest power.
 deepSplit = 2,           ##<< Course grain control of module size
 minModuleSize=30,        ##<< The minimum module size allowed
 mergeThreshold=0.1,      ##<< Below this threshold, modules are merged.
 clusterType="average",   ##<< Clustering option
 performTOPermtest=TRUE,  ##<< Performs permutation testing on modules
 toPermTestPermutes=100   ##<< Number of permutations to do.
 ){

  print("Constructing New ProCoNA Object")
  pnet = new("pepnet")
  pnet@networkName=networkName
  pnet@samples=rownames(pepdat)
  if (is.null(pow)) {
    print("Computing soft threshold power")
    softThreshold <- pickSoftThreshold(pepdat, powerVector=1:powMax,
                                       RsquaredCut=scaleFreeThreshold,
                                       networkType=signed)
    pnet@power=softThreshold$powerEstimate
  } else {
    pnet@power=pow
  }
  cat("Using power: ", pnet@power, "\n")
  pnet@peptides=colnames(pepdat)

  print("Computing adjacency")
  if (pearson) {
    pnet@adjacency <- adjacency(datExpr=pepdat,
                                power=pnet@power,
                                type=signed,
                                corOptions="use='pairwise.complete.obs'")    
  } else {
    pnet@adjacency <- adjacency(datExpr=pepdat,
                                power=pnet@power,
                                type=signed,
                                corFnc="bicor",
                                corOptions="use='pairwise.complete.obs'")
  }

  print("Computing TOM")
  pnet@TOM = TOMsimilarity(pnet@adjacency, TOMType=signed);
  rownames(pnet@TOM) <- pnet@peptides
  colnames(pnet@TOM) <- pnet@peptides
  rownames(pnet@adjacency) <- pnet@peptides
  colnames(pnet@adjacency) <- pnet@peptides
  dissTOM = 1-pnet@TOM

  print("Clustering")
  pnet@pepTree = flashClust(as.dist(dissTOM), method = clusterType);
  pnet@dynamicColors = cutreeDynamic(dendro = pnet@pepTree,
    distM = dissTOM,
    deepSplit = deepSplit, pamRespectsDendro = TRUE,
    minClusterSize = minModuleSize);
  print(table(pnet@dynamicColors))
                                        #merging modules
  print("Merging modules")
  MEDissThres = mergeThreshold
  MEList = moduleEigengenes(pepdat, colors = pnet@dynamicColors)
  pnet@MEs = MEList$eigengenes   # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(pnet@MEs)       # Cluster module eigengenes
  METree = flashClust(as.dist(MEDiss),
    method = clusterType );      # Call an automatic merging function
  merge = mergeCloseModules(pepdat, pnet@dynamicColors,
    cutHeight = MEDissThres, verbose = 4, relabel=T)    # The merged module colors
  pnet@mergedColors = merge$colors;     # Eigengenes of the new merged modules:
  pnet@mergedMEs = merge$newMEs;        # Construct numerical labels corresponding to the colors
  pnet@colorOrder = c("grey", standardColors(50));

  print("Topological Overlap Permutation Test On Modules")
  if(performTOPermtest) {
    pnet <- toPermTest(pnet, pepdat, toPermTestPermutes)
  }

  print("DONE!")
  return(pnet)
  ### returns the procona network object
}
