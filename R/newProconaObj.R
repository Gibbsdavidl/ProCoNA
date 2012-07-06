
proconaVersion <- function
### Keep track of the version of procona used to build networks.
() {
  return("0.13")
}

newProconaObj <- function
### This function returns a peptide co-expression network object.
(networkName="procona",   ##<< Name of this network
 pepdat=NULL,             ##<< This variable is the data set with rows as samples and cols as peptides
 pow=NULL,                ##<< The scaling power, NULL if unknown
 powMax=20,               ##<< The maximum power to be searched.
 networkType="signed",    ##<< Whether the sign is considered in constructing adjacency and TOM
 pearson=F,               ##<< use Pearson's cor or the robust bi-weight correlation
 scaleFreeThreshold=0.8,  ##<< The threshold for fitting to scale-free topology.. will use closest power.
 deepSplit = 2,           ##<< Course grain control of module size
 minModuleSize=30,        ##<< The minimum module size allowed
 mergeThreshold=0.1,      ##<< Below this threshold, modules are merged.
 clusterType="average",   ##<< Clustering option
 pamRespectsDendro=T,     ##<< When cutting the dendrogram, pay attention to branch membership.
 performTOPermtest=TRUE,  ##<< Performs permutation testing on modules
 toPermTestPermutes=100   ##<< Number of permutations to do.
 ){

  print("Constructing New ProCoNA Object")
  pnet = new("pepnet")
  pnet@proconaVersion = proconaVersion()
  pnet@networkName=networkName
  pnet@networkType=networkType;
  pnet@samples=rownames(pepdat)
  if (is.null(pow)) {
    print("Computing soft threshold power")
    if(pearson) {
      softThreshold <- pickSoftThreshold(pepdat, powerVector=1:powMax,
                                         RsquaredCut=scaleFreeThreshold,
                                         networkType=networkType)
    } else {
      softThreshold <- pickSoftThreshold(pepdat, powerVector=1:powMax,
                                         RsquaredCut=scaleFreeThreshold,
                                         networkType=networkType, corFnc="bicor")
    }
    pnet@power=softThreshold$powerEstimate
    if(is.na(pnet@power)) {
      print("Power Not Found!")
      return(NULL)
    }
  } else {
    pnet@power=pow
  }
  cat("Using power: ", pnet@power, "\n")
  pnet@peptides=colnames(pepdat)

  print("Computing adjacency")
  if (pearson) {
    pnet@adjacency <- adjacency(datExpr=pepdat,
                                power=pnet@power,
                                type=networkType,
                                corOptions="use='pairwise.complete.obs'")    
  } else {
    pnet@adjacency <- adjacency(datExpr=pepdat,
                                power=pnet@power,
                                type=networkType,
                                corFnc="bicor",
                                corOptions="use='pairwise.complete.obs'")
  }

  print("Computing TOM")
  pnet@TOM = TOMsimilarity(pnet@adjacency, TOMType=networkType);
  rownames(pnet@TOM) <- pnet@peptides
  colnames(pnet@TOM) <- pnet@peptides
  rownames(pnet@adjacency) <- pnet@peptides
  colnames(pnet@adjacency) <- pnet@peptides
  dissTOM = 1-pnet@TOM

  print("Clustering")
  pnet@pepTree = flashClust(as.dist(dissTOM), method = clusterType);
  pnet@dynamicColors = cutreeDynamic(dendro = pnet@pepTree,
    distM = dissTOM,
    deepSplit = deepSplit,
    pamRespectsDendro = pamRespectsDendro,
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
