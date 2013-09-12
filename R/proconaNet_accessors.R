
# Accessor functions for the proconaNet object. #

# the peptides
setGeneric("peptides", function(x) standardGeneric("peptides"))
setMethod("peptides", signature(x="proconaNet"), function(x) slot(x, "peptides"))

setGeneric("peptides<-", function(x, value) standardGeneric("peptides<-"))
setReplaceMethod("peptides", "proconaNet",
                 function(x, value) {x@peptides <- value; validObject(x); x})


### The names of the peptides used to build the peptideNet object.
setGeneric("proconaVersion", function(x) standardGeneric("proconaVersion"))
setMethod("proconaVersion", "proconaNet", function(x) x@proconaVersion)

setGeneric("proconaVersion<-", function(x, value) standardGeneric("proconaVersion<-"))
setReplaceMethod("proconaVersion", "proconaNet",
                 function(x, value) {x@proconaVersion <- value; validObject(x); x})


### The networkName of the proconaNet.
setGeneric("networkName", function(x) standardGeneric("networkName"))
setMethod("networkName", signature(x="proconaNet"), function(x) slot(x, "networkName"))

setGeneric("networkName<-", function(x, value) standardGeneric("networkName<-"))
setReplaceMethod("networkName", "proconaNet",
                 function(x, value) {x@networkName <- value; validObject(x); x})


### The type of correlations used in proconaNet.
setGeneric("networkType", function(x) standardGeneric("networkType"))
setMethod("networkType", signature(x="proconaNet"), function(x) slot(x, "networkType"))

setGeneric("networkType<-", function(x, value) standardGeneric("networkType<-"))
setReplaceMethod("networkType", "proconaNet",
                 function(x, value) {x@networkType <- value; validObject(x); x})


### The sample names used to build the network. Useful for checking phenotype ordering.
setGeneric("samples", function(x) standardGeneric("samples"))
setMethod("samples", signature(x="proconaNet"), function(x) slot(x, "samples"))

setGeneric("samples<-", function(x, value) standardGeneric("samples<-"))
setReplaceMethod("samples", "proconaNet",
                 function(x, value) {x@samples <- value; validObject(x); x})


### this function returns the adjacency matrix
setGeneric("adj", function(x) standardGeneric("adj"))
setMethod("adj", signature(x="proconaNet"), function(x) slot(x, "adj"))

setGeneric("adj<-", function(x, value) standardGeneric("adj<-"))
setReplaceMethod("adj", "proconaNet",
                 function(x, value) {x@adj <- value; validObject(x); x})


### The topological overlap matrix.
setGeneric("TOM", function(x) standardGeneric("TOM"))
setMethod("TOM", signature(x="proconaNet"), function(x) slot(x, "TOM"))

setGeneric("TOM<-", function(x, value) standardGeneric("TOM<-"))
setReplaceMethod("TOM", "proconaNet",
                 function(x, value) {x@TOM <- value; validObject(x); x})


### The peptide dendrogram.
setGeneric("pepTree", function(x) standardGeneric("pepTree"))
setMethod("pepTree", signature(x="proconaNet"), function(x) slot(x, "pepTree"))

setGeneric("pepTree<-", function(x, value) standardGeneric("pepTree<-"))
setReplaceMethod("pepTree", "proconaNet",
                 function(x, value) {x@pepTree <- value; validObject(x); x})


### The module eigenvectors after the module merging step.
setGeneric("mergedMEs", function(x) standardGeneric("mergedMEs"))
setMethod("mergedMEs", signature(x="proconaNet"), function(x) slot(x, "mergedMEs"))

setGeneric("mergedMEs<-", function(x, value) standardGeneric("mergedMEs<-"))
setReplaceMethod("mergedMEs", "proconaNet",
                 function(x, value) {x@mergedMEs <- value; validObject(x); x})


### The vector of module assignments after the merging step.
setGeneric("mergedColors", function(x) standardGeneric("mergedColors"))
setMethod("mergedColors", signature(x="proconaNet"), function(x) slot(x, "mergedColors"))

setGeneric("mergedColors<-", function(x, value) standardGeneric("mergedColors<-"))
setReplaceMethod("mergedColors", "proconaNet",
                 function(x, value) {x@mergedColors <- value; validObject(x); x})


### The module eigenvectors after the module merging step.
setGeneric("MEs", function(x) standardGeneric("MEs"))
setMethod("MEs", signature(x="proconaNet"), function(x) slot(x, "MEs"))

setGeneric("MEs<-", function(x, value) standardGeneric("MEs<-"))
setReplaceMethod("MEs", "proconaNet",
                 function(x, value) {x@MEs <- value; validObject(x); x})


### The vector of module assignments after the merging step.
setGeneric("dynamicColors", function(x) standardGeneric("dynamicColors"))
setMethod("dynamicColors", signature(x="proconaNet"), function(x) slot(x, "dynamicColors"))

setGeneric("dynamicColors<-", function(x, value) standardGeneric("dynamicColors<-"))
setReplaceMethod("dynamicColors", "proconaNet",
                 function(x, value) {x@dynamicColors <- value; validObject(x); x})


### this function returns the names of the peptides in the pepnet object
setGeneric("networkPower", function(x) standardGeneric("networkPower"))
setMethod("networkPower", signature(x="proconaNet"), function(x) slot(x, "power"))

setGeneric("networkPower<-", function(x, value) standardGeneric("networkPower<-"))
setReplaceMethod("networkPower", "proconaNet",
                 function(x, value) {x@power <- value; validObject(x); x})


### The results of the module-wise topological overlap permutation test.
setGeneric("permtest", function(x) standardGeneric("permtest"))
setMethod("permtest", signature(x="proconaNet"), function(x) slot(x, "permtest"))

setGeneric("permtest<-", function(x, value) standardGeneric("permtest<-"))
setReplaceMethod("permtest", "proconaNet",
                 function(x, value) {x@permtest <- value; validObject(x); x})


### The color order..
setGeneric("colorOrder", function(x) standardGeneric("colorOrder"))
setMethod("colorOrder", signature(x="proconaNet"), function(x) slot(x, "colorOrder"))

setGeneric("colorOrder<-", function(x, value) standardGeneric("colorOrder<-"))
setReplaceMethod("colorOrder", "proconaNet",
                 function(x, value) {x@colorOrder <- value; validObject(x); x})
