


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
### This function returns a peptide co-expression network object.
(pepdat,                  ##<< This variable is the data set with rows as samples and cols as peptides
 networkType="signed",    ##<< Whether the sign is considered in constructing adjacency and TOM
 threshold = 0.0001,      ##<< When to stop running
 tmpSaveFile = T          ##<< Should temporary saves be done?
 ){
  
  n           <- ncol(pepdat)
  Mk          <- mat.or.vec(n,n)
  Sk          <- mat.or.vec(n,n)
  runningSD   <- mat.or.vec(n,n)
  runningMean <- mat.or.vec(n,n)
  oldSD       <- mat.or.vec(n,n)
  oldMean     <- mat.or.vec(n,n)
  oldMean[1,1] <- 1
  k           <- 1
  pnet = new("pepnet")  
  print("Resampling is ... GO!")  

  while(any(abs(runningMean - oldMean) > threshold)) {

    # save the results from the last run #
    oldMean <- runningMean
    oldSD <- runningSD

    # RESAMPLING #
    resamp <- sample(1:nrow(pepdat), size=nrow(pepdat), replace=T)
    dat <- pepdat[resamp,]

    # This Adjacency #    
    pnet@adjacency <- adjacency(datExpr=dat, power=pnet@power, type=networkType,
                                corFnc="bicor", corOptions="use='pairwise.complete.obs'")
    if (k==1) { # first run...
      Mk <- pnet@adjacency  #Sk is zero here
      runningMean <- pnet@adjacency  #runningSD still 0
    } else { # k > 1
      x <- runningStats(pnet@adjacency, runningMean, Mk, Sk, k)
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


parallelCorBoostrap <- function
### This function returns a peptide co-expression network object.
(pepdat,                  ##<< This variable is the data set with rows as samples and cols as peptides
 networkType="signed",    ##<< Whether the sign is considered in constructing adjacency and TOM
 threads=1,               ##<< NUMBER OF Threads
 threshold=0.0001,        ##<< when to stop each thread.
 tmpSaveFiles=T
 ){
  stopifnot(library(multicore, logical.return=T))
    
  ff <- function(i) {
    i <- i+1
    corBootstrap(
      pepdat,   
      networkType,
      threshold)
  }
  
  x <- mclapply(X=1:threads, FUN=ff)
  x
### returns the bootstrapped adjacency matrix and sd matrix
}



newBootstrapProconaObj <- function
### This function returns a peptide co-expression network object.
(networkName="procona",   ##<< Name of this network
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
 bootstrapThreshold=0.0001,  ##<< When to stop resampling...
 parallel=F,              ##<< Should we compute the bootstrap in parallel?
 parallelThreads=1        ##<< How many threads?
 ){

  print("Constructing New ProCoNA Object")
  pnet = new("pepnet")
  pnet@proconaVersion = proconaVersion()
  pnet@networkName=networkName
  pnet@networkType=networkType;
  pnet@samples=rownames(pepdat)

  print("Computing adjacency")
  bootstrapCor <- corBootstrap(pepdat, networkType, bootstrapThreshold)
  pnet@adjacency <- bootstrapCor[[1]]
  
  if (is.null(pow)) {
    print("Computing soft threshold power")
    softThreshold <- pickSoftThreshold.fromSimilarity(pnet@adjacency,
                                                      powerVector=1:powMax,
                                                      RsquaredCut=scaleFreeThreshold,
                                                      networkType=networkType)
    pnet@power=softThreshold$powerEstimate
  } else {
    pnet@power=pow
  }
  
  cat("Using power: ", pnet@power, "\n")
  pnet@peptides=colnames(pepdat)
  pnet@adjacency <- pnet@adjacency^pnet@power

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
