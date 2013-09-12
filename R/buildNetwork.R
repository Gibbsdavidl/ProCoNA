
proconaVersionFun <- function
### Keep track of the version of procona used to build networks.
() {
  return("0.99.1")
}

### This function returns a peptide co-expression network object.

setGeneric("buildProconaNetwork",
           valueClass = "proconaNet",
           function(networkName, pepdat, pow, powMax, networkType, pearson, scaleFreeThreshold, deepSplit, minModuleSize,
                   mergeThreshold, clusterType, pamRespectsDendro, performTOPermtest, toPermTestPermutes) {
               standardGeneric("buildProconaNetwork")
           })

setMethod("buildProconaNetwork",
          signature(networkName="character",  ##<< Name of this network
                    pepdat="matrix",           ##<< This variable is the data set with rows as samples and cols as peptides
                    pow="numeric",             ##<< The scaling power, NULL if unknown
                    powMax="numeric",          ##<< The maximum power to be searched.
                    networkType="character",   ##<< Whether the sign is considered in constructing adjacency and TOM
                    pearson="logical",         ##<< use Pearson's cor or the robust bi-weight correlation
                    scaleFreeThreshold="numeric",  ##<< The threshold for fitting to scale-free topology.. will use closest power.
                    deepSplit="numeric",       ##<< Course grain control of module size
                    minModuleSize="numeric",   ##<< The minimum module size allowed
                    mergeThreshold="numeric",  ##<< Below this threshold, modules are merged.
                    clusterType="character",   ##<< Clustering option
                    pamRespectsDendro="logical", ##<< When cutting the dendrogram, pay attention to branch membership.
                    performTOPermtest="logical", ##<< Performs permutation testing on modules
                    toPermTestPermutes="numeric"), ##<< Number of permutations to do.
          function(networkName, pepdat, pow, powMax, networkType, pearson, scaleFreeThreshold, deepSplit, minModuleSize,
                   mergeThreshold, clusterType, pamRespectsDendro, performTOPermtest, toPermTestPermutes) {

              #error checking
              #args <- as.list(match.call(expand.dots = TRUE)[-1])
              #prebuild_check(args,pepdat)
              
              print("Constructing New ProCoNA Network Object")
              pnet = new("proconaNet")
              proconaVersion(pnet) = proconaVersionFun()
              networkName(pnet)=networkName
              networkType(pnet)=networkType
              samples(pnet)=rownames(pepdat)
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
                  networkPower(pnet)=softThreshold$powerEstimate
                  if(is.na(networkPower(pnet))) {
                      print("Power Not Found!")
                      return(NULL)
                  }
              } else {
                  networkPower(pnet)=pow
              }
              cat("Using power: ", networkPower(pnet), "\n")
              peptides(pnet)=colnames(pepdat)
              
              print("Computing adjacency")
              if (pearson) {
                  adj(pnet) <- adjacency(datExpr=pepdat,
                                         power=networkPower(pnet),
                                         type=networkType,
                                         corOptions="use='pairwise.complete.obs'")    
              } else {
                  adj(pnet) <- adjacency(datExpr=pepdat,
                                         power=networkPower(pnet),
                                         type=networkType,
                                         corFnc="bicor",
                                         corOptions="use='pairwise.complete.obs'")
              }
              
              print("Computing TOM")
              TOM(pnet) = TOMsimilarity(adj(pnet), TOMType=networkType);
              rownames(TOM(pnet)) <- peptides(pnet)
              colnames(TOM(pnet)) <- peptides(pnet)
              rownames(adj(pnet)) <- peptides(pnet)
              colnames(adj(pnet)) <- peptides(pnet)
              dissTOM = 1-TOM(pnet)
              
              print("Clustering")
              pepTree(pnet) = flashClust(as.dist(dissTOM), method = clusterType);
              dynamicColors(pnet) = cutreeDynamic(dendro = pepTree(pnet),
                               distM = dissTOM,
                               deepSplit = deepSplit,
                               pamRespectsDendro = pamRespectsDendro,
                               minClusterSize = minModuleSize);
              print(table(dynamicColors(pnet)))
                                        #merging modules
              print("Merging modules")
              MEDissThres = mergeThreshold
              MEList = moduleEigengenes(pepdat, colors = dynamicColors(pnet))
              MEs(pnet) = MEList$eigengenes   # Calculate dissimilarity of module eigengenes
              MEDiss = 1-cor(MEs(pnet))       # Cluster module eigengenes
              METree = flashClust(as.dist(MEDiss),
                  method = clusterType );      # Call an automatic merging function
              merge = mergeCloseModules(pepdat, dynamicColors(pnet),
                  cutHeight = MEDissThres, verbose = 4, relabel=T)    # The merged module colors
              mergedColors(pnet) = merge$colors;     # Eigengenes of the new merged modules:
              mergedMEs(pnet) = merge$newMEs;        # Construct numerical labels corresponding to the colors
              colorOrder(pnet) = c("grey", standardColors(50));
              
              print("Topological Overlap Permutation Test On Modules")
              if(performTOPermtest) {
                  pnet <- toPermTest(pnet, toPermTestPermutes)
              }
              
              print("DONE!")
              validObject(pnet)
              return(pnet)
              ### returns the procona network object
          }
          )
