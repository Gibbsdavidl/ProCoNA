
plotNet <- function
(object
 ) {
  plotDendroAndColors(dendro=pepTree(object), colors=mergedColors(object), dendroLabels=F)
}


printNet <- function(object) {
  cat("Network Name: ", networkName(object), "\n")
  cat("  Number of samples : ", length(samples(object)), "\n")
  cat("  Number of peptides: ", length(peptides(object)), "\n")
  cat("  Number of modules : ", length(unique(mergedColors(object))), "\n")
  cat("  Power used        : ", networkPower(object), "\n")
  cat("  Network Type      : ", networkType(object), "\n")
  cat("  ProCoNA version   : ", proconaVersion(object), "\n")
}
 
