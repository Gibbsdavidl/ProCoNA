
toPermTest <- function
### Uses the procona network object,
### the data with peptides as columns, samples in rows.
### And the power that the net was built at
### the number of permutations to do...
### Modules are permuted and mean topological overlap
### is recorded, constructing the null.  The number
### of random permutations with mean TO greater than
### observed provides the p-value.
(
 pnet, ##<< to get the colors
 dat,  ##<< the data used to build the model
 numPermutes=100  ##<< The number of permutations to perform
 ) {
  colors = pnet@mergedColors
  tom = pnet@TOM
  u_colors = unique(colors)
  u_colors = u_colors[u_colors != "NA" & u_colors != "0"]
  sizes   = vector("numeric", length(u_colors))
  pvalues = vector("numeric", length(u_colors))
  meantom = vector("numeric", length(u_colors))
  meanran = vector("numeric", length(u_colors))
  
  for(i in 1:length(u_colors)) {  # for each unique color
    color = u_colors[i]           #   it's this color

    idx = colors == color 
    tomSubset = tom[idx,idx]
    meanTO = mean(tomSubset)
    size = sum(idx)           # idx is boolean with length == # peptides
    x = 1:length(idx)         # indices to sample from

    cat("Permuting module: ", color, "\n")
    cat("   dim of tom subset: ", dim(tomSubset), "\n")

                                        # create null distribution
    dist = vector("numeric", numPermutes+1)
    dist[numPermutes+1] = meanTO
    for(j in 1:numPermutes) {
      pdx = sample(x, size, replace=F)
      ptom = tom[pdx,pdx]
      pmean = mean(ptom)
      #print(paste(i, ": ", round(meanTO,3), ", ", j, ": ", round(pmean,3)))
      dist[j] = pmean
    }

    # how many times did the permuted mean TO measure
    # greater than the real module mean TO?
    pvalues[i] = sum(meanTO <= dist)/numPermutes
    sizes[i]   = size
    meantom[i] = meanTO
    meanran[i] = mean(dist)
  }

  to_test = cbind(u_colors, sizes, meantom, meanran, pvalues)
  colnames(to_test) = c("module", "module size", "TOMmean", "PermMean", "TO p-value")
  pnet@permtest = to_test
  return(pnet)
  ### returns the network obj with the perm test
}


peptideConnectivityTest <- function
### This function will compare the connectivity between
### peptides linked to a given protein, against a randomly
### drawn, similarly sized, selection of peptides.
### The hypothesis is that peptides from a given protein
### should be more connected than random.
(pnet,    ##<< The peptide net object
 pepInfo, ##<< The peptide information table, mapping peptides to proteins
 pepCol,   ##<< The string identifying the column in the pepInfo table with peptide ID
 protCol   ##<< String identifying column in pepInfo with Protein ID.
 ) {

  pepInfo = pepInfo[which(pepInfo[,pepCol] %in% pnet@peptides),] 
  uProts  = unique(pepInfo[,protCol])  # The unique proteins in the network
  cat(length(uProts), " number of proteins being investigated.\n")
  connPeps = vector("numeric", length(uProts))  # peptides associated with uProt
  randPeps = vector("numeric", length(uProts))  # random peptide connections
  i = 1
  numPeps = 0
  for (p in uProts) {
    pepIdx = which(pepInfo[,protCol] == p)
    ids    = as.character(pepInfo[pepIdx,pepCol])
    x      = which(pnet@peptides %in% ids) 
    if (length(x) > 1) {
      numPeps = numPeps + length(x)
      m = pnet@TOM[x,x]
      connPeps[i] = mean(m[upper.tri(m)])
      l = sample(1:nrow(pnet@TOM), size=length(x), replace=F)
      randM = pnet@TOM[l,l]
      randPeps[i] = mean(randM[upper.tri(randM)])
    } else {
      connPeps[i] = NA
      randPeps[i] = NA
    }
    i = i+1
  }
  cat(numPeps, " number of peptides associated with these proteins.\n")
  print(t.test(connPeps, randPeps))
  return(list(connPeps, randPeps))
  ### Returns a list of the connected peptides and the random samples.
}
    

peptideCorrelationTest <- function
### Take the data, and a mapping of peptides to
### proteins, and compute the correlation between
### peptides linked to a given protein. Compare
### to random.
(dat,     ##<< The data with samples as rows and peptides as columns
 pepinfo, ##<< The mapping of peptides to proteins as a data frame
 pepCol,  ##<< The column name of peptide info table containing peptide IDs
 protCol  ##<< The column name of pepinfo info table containing protein IDs
 ) {

  uprots <- unique(pepinfo[,protCol]) 
  protcors <- c()
  randcors <- c()
  
  for(p in uprots) {
    # The peps of interest for protein p
    peps     <- pepinfo[which(pepinfo[,protCol] == p), pepCol]
    pepidx   <- which(colnames(dat) %in% peps)
    cat("protein: ", p, " has ", length(peps), " peptides... ", length(pepidx), "are in the data.\n")

    if (length(pepidx) > 1) {
      pcors    <- cor(dat[, pepidx], use="pairwise.complete.obs")
      protcors <- c(protcors, pcors[upper.tri(pcors)])
      rcors    <- cor(dat[,sample(1:ncol(dat), size=length(pepidx))], use="pairwise.complete.obs")
      randcors <- c(randcors, rcors[upper.tri(rcors)])
    }
  }

  print(t.test(protcors, randcors, na.rm=T))
  
  return(list(protcors, randcors))
  ### return a list of protein correlations and random correlations
}


