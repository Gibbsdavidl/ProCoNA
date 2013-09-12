

# unit test for the topological overlap permutation test

test_peptideCorrelationTest <- function() {
  data(ProCoNA_Data)
  
  #random test#
  randomdat <- matrix(runif(length(peptideData)), ncol=ncol(peptideData))
  colnames(randomdat) <- colnames(peptideData)
  res0 <- peptideCorrelationTest(randomdat, masstagdb, "Mass_Tag_ID", "Reference")

  #perfect correlation test#
  prots <- unique(masstagdb$Reference[masstagdb$Mass_Tag_ID %in% colnames(peptideData)])
  dat <- peptideData
  for (p in prots) {
    idx <- which(colnames(peptideData) %in% masstagdb$Mass_Tag_ID[masstagdb$Reference == p])
    central <- rnorm(nrow(dat))
    for (i in idx) {
      dat[,i] <- central + rnorm(mean=0.01, sd=0.05, n=nrow(dat))
    }
  }
  res1 <- peptideCorrelationTest(dat, masstagdb, "Mass_Tag_ID", "Reference")

  checkTrue(res0$p.value > 0.05)
  checkTrue(res1$p.value < 0.05)
}
