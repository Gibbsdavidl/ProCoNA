
proconaVersionFun <- function
### Keep track of the version of procona used to build networks.
() {
  return("0.99.1")
}

### This function returns a peptide co-expression network object.


### This function returns a peptide co-expression network object.
buildProconaNetwork <- function
(networkName="ProCoNA",   ##<< Name of this network
 pepdat,                  ##<< This variable is the data set with rows as samples and cols as peptides
 pow=1,                   ##<< The scaling power, NULL if unknown
 powMax=20,               ##<< The maximum power to be searched.
 networkType="signed",    ##<< Whether the sign is considered in constructing adjacency and TOM
 pearson=FALSE,           ##<< use Pearson's cor or the robust bi-weight correlation
 scaleFreeThreshold=0.8,  ##<< The threshold for fitting to scale-free topology.. will use closest power.
 deepSplit = 2,           ##<< Course grain control of module size
 minModuleSize=30,        ##<< The minimum module size allowed
 mergeThreshold=0.1,      ##<< Below this threshold, modules are merged.
 clusterType="average",   ##<< Clustering option
 pamRespectsDendro=TRUE,     ##<< When cutting the dendrogram, pay attention to branch membership.
 performTOPermtest=TRUE,  ##<< Performs permutation testing on modules
 toPermTestPermutes=100   ##<< Number of permutations to do.
 ) {
                                        #error checking
                                        #args <- as.list(match.call(expand.dots = TRUE)[-1])
                                        #prebuild_check(args,pepdat)
    
    message("Constructing New ProCoNA Network Object")

    if (inherits(pepdat, "MSnSet")) {
        pepdat <- t(exprs(pepdat))
    }

    pnet = new("proconaNet")
    proconaVersion(pnet) = proconaVersionFun()
    networkName(pnet)=networkName
    networkType(pnet)=networkType
    samples(pnet)=rownames(pepdat)
    if (is.null(pow)) {
        message("Computing soft threshold power")
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
            error("Power Not Found!")
        }
    } else {
        networkPower(pnet)=pow
    }
    cat("Using power: ", networkPower(pnet), "\n")
    peptides(pnet)=colnames(pepdat)
    
    message("Computing adjacency")
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
    
    message("Computing TOM")
    TOM(pnet) = TOMsimilarity(adj(pnet), TOMType=networkType);
    rownames(TOM(pnet)) <- peptides(pnet)
    colnames(TOM(pnet)) <- peptides(pnet)
    rownames(adj(pnet)) <- peptides(pnet)
    colnames(adj(pnet)) <- peptides(pnet)
    dissTOM = 1-TOM(pnet)
    
    message("Clustering")
    pepTree(pnet) = flashClust(as.dist(dissTOM), method = clusterType);
    dynamicColors(pnet) = cutreeDynamic(dendro = pepTree(pnet),
                     distM = dissTOM,
                     deepSplit = deepSplit,
                     pamRespectsDendro = pamRespectsDendro,
                     minClusterSize = minModuleSize);
    print(table(dynamicColors(pnet)))
                                        #merging modules
    message("Merging modules")
    MEDissThres = mergeThreshold
    MEList = moduleEigengenes(pepdat, colors = dynamicColors(pnet),
        softPower=networkPower(pnet))
    MEs(pnet) = MEList$eigengenes   # Calculate dissimilarity of module eigengenes
    MEDiss = 1-cor(MEs(pnet))       # Cluster module eigengenes
    METree = flashClust(as.dist(MEDiss),
        method = clusterType );      # Call an automatic merging function
    merge = mergeCloseModules(pepdat, dynamicColors(pnet),
        cutHeight = MEDissThres, verbose = 4, relabel=TRUE)    # The merged module colors
    mergedColors(pnet) = merge$colors;     # Eigengenes of the new merged modules:
    mergedMEs(pnet) = merge$newMEs;        # Construct numerical labels corresponding to the colors
    colorOrder(pnet) = c("grey", standardColors(50));
    
    message("Topological Overlap Permutation Test On Modules")
    if(performTOPermtest) {
        pnet <- toPermTest(pnet, toPermTestPermutes)
    }
    
    message("DONE!")
    validObject(pnet)
    return(pnet)
### returns the procona network object
}


