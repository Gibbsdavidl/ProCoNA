
subsetModCors <- function
### subsets the module-phenotype correlation matrix which has funny rownames
(modCors,  ##<< The matrix of module-phenotype correlations
 modules   ##<< Which modules are desired.
 ) {
  rows <- rownames(modCors)  # like ME3
  ms <- sapply(modules, function(x) paste("ME", x, sep=""))
  modCors[which(rows %in% ms),]
}



setGeneric("correlationWithPhenotypesHeatMap",
           function(net,phenotypes,modules,plotName,title,textSize) {
               standardGeneric("correlationWithPhenotypesHeatMap")
           })

setMethod("correlationWithPhenotypesHeatMap",
          signature(net="proconaNet",
                    phenotype="matrix",
                    modules="numeric",
                    plotName="character",
                    title="character",
                    textSize="numeric"),
         
          function
### The heatmap showing the relation of the
### modules with the phenotypes computed by
### using the function relateModules(...)

          (net,               ##<< the peptide network object
           phenotypes,       ##<< matrix of phenotypic traits
           modules = c(),   ##<< the vector of modules to plot
           plotName="",    ##<< the name of the saved plot, NULL to show on screen
           title="Module-trait relationships", ##<< plot main title
           textSize=0.4
           ) {
              
              if (nrow(phenotypes) != nrow(mergedMEs(net))) {
                  stop("Net must have merged eigenvectors with dimensions matching the phenotype data! See mergedMEs(net).")
              }
              
              if (is.null(modules)) {modules <- unique(mergedColors(net))}
              modCors <- modulePhenotypeCorrelations(net, phenotypes)
              modCors <- subsetModCors(modCors, modules)
              mms <- 1:(ncol(modCors)/2)
              mmps <- ((ncol(modCors)/2)+1):ncol(modCors)
              modCorMM <- as.matrix(modCors[,mms])
              modCorPS <- as.matrix(modCors[,mmps])
              
                                        # pastes into a vector down columns
              textMatrix <- paste(signif(modCorMM, 2), "\n(", sep="")
              textMatrix <- paste(textMatrix, paste(signif(modCorPS, 1), ")", sep = ""), sep="")
              textMatrix <- matrix(textMatrix, byrow=F, nrow=nrow(modCorMM))  
              
              par(mar = c(6, 8.5, 3, 3));
                                        # Display the correlation values within a heatmap plot
              if(!is.null(plotName)) {pdf(plotName)};
              
              labeledHeatmap(Matrix = modCors[,mms],
                             xLabels = names(phenotypes),
                             yLabels = rownames(modCors),
                             ySymbols = rownames(modCors),
                             colorLabels = FALSE,
                             colors = blueWhiteRed(50),
                             textMatrix = textMatrix,
                             setStdMargins = FALSE,
                             cex.text = textSize,
                             zlim = c(-1,1),
                             main = title)
              
              if(!is.null(plotName)) {dev.off()}
              
              return(modCors)
### the module eigenvector correlations
          })


setGeneric("MMvsPS",
           function(pnet,pepdat, phenoVec, mod) {
               standardGeneric("MMvsPS")
           })

setMethod("MMvsPS",
          signature(pnet="proconaNet",
                    pepdat="matrix",
                    phenoVec="numeric",
                    mod="numeric"),
          function
### Plots the module membership against the peptide significance
### for a given trait and module
          (pnet,      ##<< The procona network
           pepdat,   ##<< the peptide data, with rows as samples and columns as peptides
           phenoVec, ##<< the phenotypic trait, vector
           mod       ##<< the module of interest
           ) {
              
              modme <- paste("ME",mod,sep="")
              
              thesePeps <- which(mergedColors(pnet) == mod)
              pepnum <- length(thesePeps)
              message(pepnum, " peptides in module ", mod, "\n")
              
              mmCors <- vector("numeric", pepnum)
              psCors <- vector("numeric", pepnum)
              
              dat <- pepdat[,thesePeps]
              
              for (i in 1:pepnum) {
                  
                  mmCors[i] <- abs(cor(dat[,i], mergedMEs(pnet)[,modme], use="pairwise.complete.obs"))
                  psCors[i] <- abs(cor(dat[,i], phenoVec, use="pairwise.complete.obs"))
                  
              }
              
              plot(x=mmCors, y=psCors, main="", xlab="Module Membership", ylab="Peptide Significance")
              
              return(list(mmCors,psCors))
### returns a list of module memberships and peptide significances.
          })



setGeneric("MMvsPSallModules",
           function(net,peptable,phenoVec,prefixName) {
               standardGeneric("MMvsPSallModules")
           })

setMethod("MMvsPSallModules",
          signature(net="proconaNet",
                    peptable="matrix",
                    phenoVec="numeric",
                    prefixName="character"),
          
          function
### Produce pdfs for all modules.
          (net,       ##<< the procona network object
           peptable,  ##<< the peptide data
           phenoVec,  ##<< the phenotypic trait, as a numeric vector
           prefixName="mm_vs_ps_" ##<< the  
           ){
              
              mods <- unique(mergedColors(net))
              
              for (m in mods) {
                  message("Printing module: ", m, "\n")
                  pdfName <- paste(prefixName, m, ".pdf", sep="")
                  pdf(pdfName)
                  MMvsPS(net, peptable, phenoVec, m)
                  dev.off()
              }
### nothing returned
          })
