
#
# These functions should be considered experimental!! #


runningStats <- function
### Computing the running mean and variance 
(newMat,        ##<< The matrix from resampled data
 runningMean,   ##<< The running mean matrix
 Mk1,           ##<< Matrix used in calculation of mean
 Sk1,           ##<< Matrix used in calculation of sd
 k              ##<< Current resampling iteration
 ) {
  Mk <- Mk1 + (newMat-Mk1)/k
  Sk <- Sk1 + (newMat-Mk1)*(newMat-Mk)
  runningSD <- Sk/(k-1)
  runningMean <- (newMat + k*runningMean)/(k+1)
  gc()
  list(runningMean,runningSD,Mk,Sk)
  ### returns the list of runningMean, runningSD, Mk, Sk
}


corBootstrap <- function
### This function returns a bootstrapped correlation matrix, std dev of correlations matrix, and number
### of samplings.
(dataMatrix,              ##<< This variable is the data set with rows as samples and cols as peptides
 networkType="signed",    ##<< Whether the sign is considered in constructing adjacency and TOM
 threshold = 0.0001,      ##<< When to stop running
 tmpSaveFile = T          ##<< Should temporary saves be done?
 ){
  
  #check_corBootstrap(tmpSaveFile,networkType,threshold,dataMatrix)
  
  n           <- ncol(dataMatrix)
  Mk          <- mat.or.vec(n,n)
  Sk          <- mat.or.vec(n,n)
  runningSD   <- mat.or.vec(n,n)
  runningMean <- mat.or.vec(n,n)
  oldSD       <- mat.or.vec(n,n)
  oldMean     <- mat.or.vec(n,n)
  oldMean[1,1] <- 1
  k           <- 1

  print("Resampling is ... GO!")  

  while(any(abs(runningMean - oldMean) > threshold)) {

    # save the results from the last run #
    oldMean <- runningMean
    oldSD <- runningSD

    # RESAMPLING #
    resamp <- sample(1:nrow(dataMatrix), size=nrow(dataMatrix), replace=T)
    dat <- dataMatrix[resamp,]

    # This Adjacency #    
    adjMat <- adjacency(datExpr=dat, power=1, type=networkType,
                        corFnc="bicor", corOptions="use='pairwise.complete.obs'")
    if (k==1) { # first run...
      Mk <- adjMat  #Sk is zero here
      runningMean <- adjMat  #runningSD still 0
    } else { # k > 1
      x <- runningStats(adjMat, runningMean, Mk, Sk, k)
      runningMean <- x[[1]]
      runningSD <- x[[2]]
      Mk <- x[[3]]
      Sk <- x[[4]]
    }

    if(tmpSaveFile == T && k%%500 == 0) {
      cat("On resampling: ", k, "\n")
      save(list(runningMean, runningSD, (k-1)),
           file="Resampling_Temp_Save_File.rda")
    }
    
    k <- k+1
    gc();
  }
  cat("All Done!", "Performed: ", k-1, " resamplings\n")
  list(runningMean, runningSD, k-1)
  ### Returns a list of the mean matrix, sd matrix, and the number of resamplings done.
}




bootstrapProconaNetwork <- function
### This function returns a peptide co-expression network object.
(networkName="bootstrap procona",   ##<< Name of this network
 pepdat=NULL,             ##<< This variable is the data set with rows as samples and cols as peptides
 pow=NULL,                ##<< The scaling power, NULL if unknown
 powMax=20,               ##<< The maximum power to be searched.
 networkType="signed",    ##<< Whether the sign is considered in constructing adjacency and TOM
 scaleFreeThreshold=0.8,  ##<< The threshold for fitting to scale-free topology.. will use closest power.
 deepSplit = 2,           ##<< Course grain control of module size
 minModuleSize=30,        ##<< The minimum module size allowed
 mergeThreshold=0.1,      ##<< Below this threshold, modules are merged.
 clusterType="average",   ##<< Clustering option
 pamRespectsDendro=T,     ##<< When cutting the dendrogram, pay attention to branch membership.
 performTOPermtest=TRUE,  ##<< Performs permutation testing on modules
 toPermTestPermutes=100,  ##<< Number of permutations to do.
 bootstrapThreshold=0.0001  ##<< When to stop resampling...
 ){

                                        #args <- as.list(match.call(expand.dots = TRUE)[-1])
                                        #prebuild_check(args,pepdat)
    
    print("Constructing New ProCoNA Object")
    pnet <- new("proconaNet")
    proconaVersion(pnet) <- proconaVersionFun()
    networkName(pnet) <- networkName
    networkType(pnet) <- networkType;
    samples(pnet) <- rownames(pepdat)
    
    print("Computing adjacency")
    bootstrapCor <- corBootstrap(pepdat, networkType, bootstrapThreshold)
    adj(pnet) <- bootstrapCor[[1]]
    
    if (is.null(pow)) {
        print("Computing soft threshold power")
        softThreshold <- pickSoftThreshold.fromSimilarity(adj(pnet),
                                                          powerVector=1:powMax,
                                                          RsquaredCut=scaleFreeThreshold,
                                                          networkType=networkType)
        networkPower(pnet)=softThreshold$powerEstimate
    } else {
        networkPower(pnet)=pow
    }
    
    cat("Using power: ", networkPower(pnet), "\n")
    peptides(pnet)=colnames(pepdat)
    adj(pnet) <- adj(pnet)^networkPower(pnet)
    
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
    
    print("Validating Object...")
                                        #check_proconaNet(pnet)
    
    print("DONE!")
    return(pnet)
### returns the procona network object
}
