i2c <- function
### Index to coordinates
(nrows, ##<< number of rows in the matrix
 i      ##<< the index into the matrix
 ) {
  y <- ceiling(i/nrows)
  x <- i%%nrows
  if (x == 0) {x <- nrows}
  c(x,y)
### the row col coordinates into the matrix
}

c2i <- function
### coordinates to index
(nrows, ##<< number of rows in the matrix
 x,     ##<< the row coordinate
 y      ##<< the col coordinate
 ) {
  (y-1) * nrows + x
### the index into the matrix
}
  
i2col <- function
### index to column
(nrows, ##<< number of rows in matrix
 i      ##<< the index
 ) {
  y <- ceiling(i/nrows)
### returns the column of the matrix
}

utri <- function
### The upper triangle of a matrix
(mat  ##<< A matrix
 ){
  mat[upper.tri(mat)]
  ### Returns a vector
}
  
