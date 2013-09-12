

# unit test for the topological overlap permutation test

test_toPermTest <- function() {

  tom <- matrix(runif(500*500), ncol=500)
  i <- 1
  for (j in seq(100, 500, 100)) {
    for (k in i:j) {
      tom[i:j, k] <- 0.9
    }
    i <- j+1
  }
  ms <- c(rep(1,100), rep(2,100), rep(3,100), rep(4,100), rep(5,100))
  net <- new("proconaNet")
  net@mergedColors <- ms
  net@TOM <- tom
  res0 <- toPermTest(net,1000)
  res0 <- res0@permtest

  checkEquals(length(res0[,1]), 5)
  checkEqualsNumeric(mean(res0[,3]), 0.9, tolerance=1.0e-4)
  checkTrue(all(res0[,5] <= 1/1000))
}
