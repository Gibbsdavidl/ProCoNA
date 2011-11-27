
# this is true of the WGCNA function adjacency.
# z <- matrix(rnorm(100), ncol=10)
# all(adjacency(datExpr=z, power=1, type="unsigned") == abs(cor(z)))

tukeyAdjacency <- function
### Adjacency matrix computed with Tukey's bi-weight correlation.
(pepdat,                 ##<< peptide data, rows = samples, cols = peptides
 type = "unsigned",      ##<< type of adjacency matrix, unsigned, signed, positiveOnly
 power = 1,              ##<< scaling power (beta)
 cores=1,                ##<< number of cores to use
 median=T,               ##<< median or mcv init?
 full.init=T             ##<< init from full data? or pairwise?           
 ) {
  require(biwt)
  require(fastukeycor)
  cormat <- partukeycor(x=t(pepdat),
                        median=median,
                        full.init=full.init,
                        output="matrix",
                        cores=cores)
  
  if (type == "unsigned") {
    cormat = abs(cormat)
  }
  else if (type == "signed") {
    cormat = (1 + cormat)/2
  }
  cormat^power
  ### returns a scaled adjacency matrix
}

tukeyDataQuality <- function 
### Check data quality by comparing Tukey and Pearson.
(pepdat ##<< peptide data
 ) {
  require(biwt)
  tukeyCorMat <- biwt.cor(x=pepdat)
  pearsonCorMat <- cor(pepdat)

  tukey <- tukeyCorMat[upper.tri(tukeyCorMat)]
  pearson <- pearsonCorMat[upper.tri(pearsonCorMat)]

  flaggedHigh <- which(abs(tukey) > 0.85 | abs(pearson) > 0.85)
  pearsonHigh <- pearson[flaggedHigh]
  tukeyHigh <- tukey[flaggedHigh]
  
  flagged <- which(abs(pearsonHigh - tukeyHigh) > 1)
  cat("Number of peptides with questionable data: ", length(flagged), "\n")
  unique(colnames(pepdat)[i2c(nrow(pepdat),flagged)])
  ### returns names of peptides that are questionable.
}
