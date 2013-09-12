

# unit test for the topological overlap permutation test

test_ppiPermTest <- function() {
  options(stringsAsFactors=F)
  #sink(file="/dev/null",type="output")
  # random ppi test
  data(ProCoNA_Data)
  net1 <- buildProconaNetwork("peptide network", peptideData)
  ppis <- data.frame(A=sample(masstagdb$Reference, 50), B=sample(masstagdb$Reference, 50))
  x <- ppiPermTest(net1, peptideData, masstagdb, "Mass_Tag_ID", "Reference", ppis, 0.33, 1000)

  # with mulitple testing, p-values should be above 0.006
  # for non-significance...
  checkTrue(all(x$pval > 0.006))


  # perfect ppi enrichment for a single module
  mod <- which(mergedColors(net1) == 2)
  modprots <- as.character(unique(masstagdb$Reference[masstagdb$Mass_Tag_ID %in% peptides(net1)[mod]]))
  a1 <- modprots[1:10]
  b1 <- modprots[3:12]
  background1 <- paste(sample(masstagdb$Reference,200), sample(masstagdb$Reference,200), sep="")
  background2 <- paste(sample(masstagdb$Reference,200), sample(masstagdb$Reference,200), sep="")
  ppi2 <- data.frame(A=c(a1,background1),B=c(b1,background2))
  x <- ppiPermTest(net1, peptideData, masstagdb, "Mass_Tag_ID", "Reference", ppi2, 0.3, 1000)
  #sink(file=NULL)
  checkTrue(x$pval[names(x$pval) == "2"] < 0.05)
}
