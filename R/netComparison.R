

orderMatrixIndex <- function
### Order the the matrix by upper diag in a greedy fashion
(
 mat ##<< A matrix
 ) {
  d1 <- dim(mat)[1] #rows
  d2 <- dim(mat)[2] #cols
  thisrow <- d1 # start at last row
  thiscol <- 1  # and first col
  neworder <- 1:d2

  # If there are more rows than columns (and vice versa)
  # order the diagonal, then tack on the left overs.
  if (d2 > d1) {
    for (i in 1:d1) {
      colmax <- max(mat[thisrow,])  # max in bottom row
      imax <- which(mat[thisrow,] == colmax)
                                        # swap the columns in new order
      x <- which(neworder == imax[1])
      y <- neworder[thiscol]
      neworder[thiscol] <- imax[1]
      neworder[x] <- y
      mat[,imax[1]] <- -1
      thisrow <- thisrow-1
      thiscol <- thiscol+1
    }
  } else {
    for (i in 1:d2) {
      print(neworder)
      colmax <- max(mat[thisrow,])  # max in bottom row
      imax <- which(mat[thisrow,] == colmax)
      x <- which(neworder == imax[1])
      y <- neworder[thiscol]
      neworder[thiscol] <- imax[1]
      neworder[x] <- y
      mat[,imax[1]] <- -1
      thisrow <- thisrow-1
      thiscol <- thiscol+1
    }
  }
  return(neworder)
  ### returns a matrix in order of greatest in upper diagonal direction.
}


setGeneric("getFisherMatrix",
           function(peps1,peps2,colors1,colors2) {
               standardGeneric("getFisherMatrix")
           })

setMethod("getFisherMatrix",
          signature(peps1="character",
                    peps2="character",
                    colors1="numeric",
                    colors2="numeric"),
          
          function ### Returns the matrix of fisher pvalues
          (peps1,    ##<<  Peps1 list of entities in the pepswork (nodes of pepswork 1)
           peps2,    ##<<  Peps2 list of entities in the pepswork (nodes of pepswork 2)
           colors1, ##<<  the module assignments for pepswork 1
           colors2  ##<<  the module assignments for pepswork 2
           ) {

              modules1 <- unique(colors1)
              modules2 <- unique(colors2)
              fishermatrix  <- mat.or.vec(length(modules1), length(modules2))
              overlapmatrix <- mat.or.vec(length(modules1), length(modules2))
              n <- length(intersect(peps1, peps2))
              
              for (i in 1:length(modules1)) {
                  for (j in 1:length(modules2)) {
                      module1 <- peps1[which(colors1 == modules1[i])]
                      module2 <- peps2[which(colors2 == modules2[j])]
                      m <- length(intersect(module1, module2))
                      n1 <- length(setdiff(module1, module2))
                      n2 <- length(setdiff(module2, module1))
                      u <- n-n1-n2-m
                      if (m == 0) {
                          fishermatrix[i,j] <- 1 # tank the result if no overlap
                      } else {
                          fisherMat <- matrix(c(m, n1, n2, u), byrow=T, nrow=2)
                          fout <- fisher.test(fisherMat, conf.int=F)
                          fishermatrix[i,j] <- fout$p.value
                      }
                      overlapmatrix[i,j] <- m
                  }
              }
              rownames(fishermatrix) <- modules1
              colnames(fishermatrix) <- modules2
              rownames(overlapmatrix) <- modules1
              colnames(overlapmatrix) <- modules2
              
              return(list(fishermatrix, overlapmatrix))
### returns the fisher test pvalues and overlaps 
})


#           signature=c("character","character","numeric","numeric",
#               "character","character","character"),
         


setGeneric("compareNetworksWithFishersExactTest",
           function
### This makes the heatmap of agreement betewen networks, module-wise,
### by comparing each using Fisher's exact test where...
### n == number of entities in the network
### m == number of entities in intersection of two modules
### d1 == number of entities in module A but not in module B
### d2 == number of entities in module B but not in module A
### 2x2 matrix for the test is then:
### m  d1
### d2 n-d1-d2-m
          (peps1, ##<< Net 1 entities, character vector
           peps2, ##<< Net 2 entiites, character vector
           colors1, ##<< modules for net 1
           colors2, ##<< modules for net 2
           title="", ##<< Plot title
           net1label="", ##<< xlabel
           net2label=""  ##<< ylabel
           ) {
  
                                        # first get the matrix of fisher test pvalues and overlaps
              m <- getFisherMatrix(peps1, peps2, colors1, colors2)
              fishers <- m[[1]]
              overlap <- m[[2]]
              
                                        # order overlaps in a greedy way, try to make nice diagonal
              idx <- orderMatrixIndex(overlap)
              fishers <- fishers[,idx]
              overlap <- overlap[,idx]
              
                                        # will use the -log of hochberg adjusted p-values
              fishLog <- matrix(p.adjust(fishers, "hochberg"), nrow=nrow(m[[1]]), byrow=F)
              fishLog <- -log10(fishLog)
              colnames(fishLog) <- colnames(fishers); rownames(fishLog) <- rownames(fishers)
              logLevels <- sort(cut(fishLog, 64)) # factors 1 - 64 for legend
              fishLogMatrix <- fishLog
              
                                        # normalize for nice colors
              fishLog[which(is.infinite(fishLog))] <- NA
              fishLog[which(is.na(fishLog))] <- max(fishLog, na.rm=T)
              fishLog <- 1-(fishLog/max(fishLog))
              d <- dim(fishLog)
              colnames(fishLog) <- colnames(fishers)
              rownames(fishLog) <- rownames(fishers)
              
                                        # the main image
              b <- seq(from=0, to=1, length=65)
              par( mar=c(5,5,5,7), xpd=T)
              image(x=1:d[1], y=1:d[2], z=fishLog, col=heat.colors(64), breaks=b, main=title,
                    xlab = net1label, ylab=net2label, axes=F)
              
                                        # add a legend with 5 levels ... legendLable == "legend label"
              legendLable <- names(table(logLevels))
              legendLable <- legendLable[quantile(1:length(legendLable))]
              legend(x=d[1]+1, y=10,legend=legendLable, cex=0.6, title="-log adj p-values",
                     fill=c(heat.colors(64)[64],
                         heat.colors(64)[48], heat.colors(64)[32],
                         heat.colors(64)[16], heat.colors(64)[1]))
              axis(1, at=1:d[1], (rownames(fishers)), cex.axis=0.65)
              axis(2, at=1:d[2], (colnames(fishers)), cex.axis=0.65)
              
                                        # add the overlaps at each intersection of modules
              for (i in 1:d[1]) {
                  for (j in 1:d[2]) {
                      text(x=i, y=j, labels = overlap[i,j], cex = 0.65)
                  }
              }
              return(list(fishLogMatrix, overlap))
### returns fishers exact test -log pvalues and overlap matrix
          })
          


setGeneric("compareNetworksWithFishersExactTestProcona",
 #          signature=c("proconaNet",
 #              "proconaNet",
 #              "character"),
           function(net1,net2,title) {
           #    standardGeneric("compareNetworksWithFishersExactTestProcona")
               compareNetworksWithFishersExactTest(peptides(net1), peptides(net2),
                                                   mergedColors(net1), mergedColors(net2),
                                                   title, networkName(net1), networkName(net2))

           })

