
check_proconaNet <- function
### Check the validity of the procona object.
(object  ##<< The procona network object.
 ){

  errors <- character()

  sample_len <- length(samples(object))
  if(sample_len < 2) {
    msg <- paste("Number of samples is: ", sample_len, "\n", sep="")
    errors <- c(errors, msg)
  }

  pep_len <- length(peptides(object))
  if(pep_len < 2) {
    msg <- paste("Number of peptides is: ", pep_len, "\n", sep="")
    errors <- c(errors, msg)
  }
  
  if(networkPower(object) < 1) {
    msg <- paste("The network power is: ", networkPower(object), ".\n", sep="")
    errors <- c(errors, msg)
  }

  if (! (networkType(object) == "signed" || networkType(object) == "unsigned")) {
      msg <-paste(networkType(object), ":  networkType must be either: \"signed\" or \"unsigned\"")
      errors <- c(errors,msg)
  }

  if(ncol(TOM(object)) < 2) {
      msg <- paste("TOM size is ", ncol(TOM(object)), "while number of peptides is ", pep_len, "\n", sep="")
      errors <- c(errors, msg)
  }

  if(ncol(adj(object)) < 2) {
      msg <- paste("Adjacency matrix size is ", ncol(adj(object)), "\n", sep="")
      errors <- c(errors, msg)
  }

  if(length(mergedColors(object)) < 2) {
    msg <- paste("The number of mergedColors is less than 2.\n", sep="")
    errors <- c(errors, msg)
  }

  if(length(dynamicColors(object)) < 2) {
    msg <- paste("The number of dynamicColors is less than 2.\n", sep="")
    errors <- c(errors, msg)
  }

  if(ncol(mergedMEs(object)) < 2) {
    msg <- "The network module eigenvectors are too few.\n"
    errors <- c(errors, msg)
}
    
  if(length(errors) == 0)
    {TRUE}
  else
    {errors}
  ### Returns TRUE if no errors detected, and a character vector of errors otherwise.
}


setValidity("proconaNet", check_proconaNet)


