
subsetModCors <- function
### subsets the module-phenotype correlation matrix which has funny rownames
(modCors,  ##<< The matrix of module-phenotype correlations
 modules   ##<< Which modules are desired.
 ) {
  rows <- rownames(modCors)  # like ME3
  ms <- sapply(modules, function(x) paste("ME", x, sep=""))
  modCors[which(rows %in% ms),]
}



correlationWithPhenotypesHeatMap <- function
### The heatmap showing the relation of the
### modules with the phenotypes computed by
### using the function relateModules(...)

(pnet,             ##<< the peptide network object
 pdat,             ##<< the peptide data
 phenotypes,       ##<< matrix of phenotypic traits
 modules,          ##<< the vector of modules to plot
 plotName=NULL,    ##<< the name of the saved plot, NULL to show on screen
 title="Module-trait relationships", ##<< plot main title
 textSize=0.4
 ) {

  modCors <- modulePhenotypeCorrelations(pnet,pdat,phenotypes)
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
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = textSize,
                 zlim = c(-1,1),
                 main = title)
  if(!is.null(plotName)) {dev.off()}

  return(modCors)
  ### the module eigenvector correlations
}



MMvsPS <- function
### Plots the module membership against the peptide significance
### for a given trait and module
(net,      ##<< The procona network
 pepdat,   ##<< the peptide data, with rows as samples and columns as peptides
 phenoVec, ##<< the phenotypic trait, vector
 mod       ##<< the module of interest
 ) {

  modme <- paste("ME",mod,sep="")
  
  thesePeps <- net@peptides[which(net@mergedColors == mod)]
  pepnum <- length(thesePeps)
  cat(pepnum, " number of peptides in module ", mod, "\n")

  mmCors <- vector("numeric", pepnum)
  psCors <- vector("numeric", pepnum)
  
  dat <- peptable[,thesePeps]
  
  for (i in 1:pepnum) {

    mmCors[i] <- abs(cor(dat[,i], net@mergedMEs[,modme], use="pairwise.complete.obs"))
    psCors[i] <- abs(cor(dat[,i], phenoVec, use="pairwise.complete.obs"))
    
  }

  plot(x=mmCors, y=psCors, main="", xlab="Module Membership", ylab="Peptide Significance")
  
  return(list(mmCors,psCors))
  ### returns a list of module memberships and peptide significances.
}



MMvsPSallModules <- function
### Produce pdfs for all modules.
(net,       ##<< the procona network object
 peptable,  ##<< the peptide data
 phenoVec,  ##<< the phenotypic trait, as a numeric vector
 prefixName="mm_vs_ps_" ##<< the  
 ){

  mods <- unique(net@mergedColors)

  for (m in mods) {
    cat("Printing module: ", m, "\n")
    pdfName <- paste(prefixName, m, ".pdf", sep="")
    pdf(pdfName)
    MMvsPS(net, peptable, phenoVec, m)
    dev.off()
  }
  ### nothing returned
}
