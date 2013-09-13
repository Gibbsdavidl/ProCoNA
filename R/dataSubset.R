

getPeptideNAs <- function
### This function returns the number of NAs for each peptide.
(pepdat ##<< the peptide data.
 ) {
  apply(pepdat, MARGIN=2, FUN= function(x) {length(which(is.na(x)))})
### returns a list of counts of NAs for each peptide.
}


subsetPeptideData <- function
### Return the smaller peptide table, omitting NAs
(pepdat, ##<< the peptide information
 numNAsAllowed=NULL, ##<< the number of NAs allowed for each peptide
 percentageNAsAllowed = 0.05 ##<< the percentage of NAs allowed for each peptide.
 ) {
  naCount = getPeptideNAs(pepdat)
  if (is.null(numNAsAllowed)) {
    nasAllowed <- as.integer(percentageNAsAllowed * nrow(pepdat))
    littlePepdat <- pepdat[,which(naCount <= nasAllowed)]
    message("Removing peptides with more than ", nasAllowed, " missing data points\n")
  } else {
    littlePepdat <- pepdat[,which(naCount <= numNAsAllowed)]
    message("Removing peptides with more than ", numNAsAllowed, " missing data points\n")
  }
  return(littlePepdat)
### Returns the subset table.
}
