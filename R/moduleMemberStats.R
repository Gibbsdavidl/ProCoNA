

setGeneric("modulePhenotypeCorrelations",
           function(pnet, phenotypes) {
               standardGeneric("modulePhenotypeCorrelations")
           })

setMethod("modulePhenotypeCorrelations",
          signature(pnet="proconaNet",
                    phenotypes="matrix"),
          
          function(pnet, phenotypes) {
                                       ### Computes the relation between the modules and
                                       ### the phenotypes.              
                                        # Code borrowed/repurposed from WGCNA tutorials
                                        # The following setting is important, do not omit.
              options(stringsAsFactors = FALSE);
              for (i in 1:ncol(phenotypes)) {
                  x <- phenotypes[,i]
                  phenotypes[,i] <- if(typeof(x) != "double") {as.numeric(factor(x))} else x
              }
              
              nSamples = length(samples(pnet));
                                        # Recalculate MEs with color labels
              MEs0 = mergedMEs(pnet)
              MEs = orderMEs(MEs0)
              moduleTraitCor = cor(MEs, phenotypes, use = "pairwise.complete.obs");
              moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
              pvalueNames <- sapply(colnames(moduleTraitCor), function(x){paste(x,".pvalue",sep="")}) 
              colnames(moduleTraitPvalue) <- pvalueNames
              
              moduleCors <- data.frame(moduleTraitCor,
                                       moduleTraitPvalue)
              
              return(moduleCors)
### returns a list of module member stats.
          })



setGeneric("moduleMemberCorrelations",
           function(pnet, pepdat, phenotypes) {
               standardGeneric("moduleMemberCorrelations")
           })

setMethod("moduleMemberCorrelations",
          signature(pnet="proconaNet",
                    pepdat="matrix",
                    phenotypes="matrix"),
          
          function(pnet, pepdat, phenotypes) {
                                        # Code borrowed/repurposed from WGCNA tutorials
                                        # The following setting is important, do not omit.
              options(stringsAsFactors = FALSE);
              for (i in 1:ncol(phenotypes)) {
                  x <- phenotypes[,i]
                  phenotypes[,i] <- if(typeof(x) != "double") {as.numeric(factor(x))} else x
              }
                                        # Define numbers of genes and samples
              nPeptides = ncol(pepdat);
              nSamples = nrow(pepdat);
                                        # Recalculate MEs with color labels
              MEs0 = mergedMEs(pnet)
              MEs = orderMEs(MEs0)
              
              modNames = substring(names(MEs), 3)
              
              peptideModuleMembership = as.data.frame(cor(pepdat, MEs, use = "pairwise.complete.obs"));
              
              MMPvalue = as.data.frame(corPvalueStudent(
                  as.matrix(peptideModuleMembership), nSamples));
              
              names(peptideModuleMembership) = paste("MM", modNames, sep="");
              names(MMPvalue) = paste("p.MM", modNames, sep="");
              
              peptideTraitSignificance = as.data.frame(cor(pepdat, phenotypes, use = "pairwise.complete.obs"));
              
              PSPvalue = as.data.frame(corPvalueStudent(
                  as.matrix(peptideTraitSignificance), nSamples));
              
              names(peptideTraitSignificance) = paste("PS.", names(phenotypes), sep="");
              
              names(PSPvalue) = paste("p.PS.", names(phenotypes), sep="");
              
              memberCors <- data.frame(Peptide=peptides(pnet),
                                       Module=mergedColors(pnet),
                                       peptideModuleMembership,
                                       MMPvalue,
                                       peptideTraitSignificance,
                                       PSPvalue)
              return(memberCors)
### returns a data frame of peptide correlation stats.
          })
