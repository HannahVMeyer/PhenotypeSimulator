## ---- echo = FALSE, message=FALSE----------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library("PhenotypeSimulator")
library("ggplot2")

## ---- eval=FALSE---------------------------------------------------------
#  # load kinship data, add small value to diagonal for numerical stability and
#  # do an eigendecomposition
#  indir <- "~/data/hmeyer/Supplementary/CEU.0908.impute.files"
#  
#  kinship <- read.table(paste(indir,
#                              "/genotypes_genome_hapgen.controls.grm.rel",sep=""),
#                        sep="\t", header=FALSE)
#  
#  kinship <- as.matrix(kinship)
#  diag(kinship) <- diag(kinship) + 1e-4
#  
#  kinship_decomposed <- eigen(kinship)
#  write.table(kinship,
#  			paste(indir, "/genotypes_genome_hapgen.controls.grm.rel",
#  				  sep=""),
#  			sep="\t", col.names=FALSE, row.names=FALSE)
#  
#  write.table(kinship_decomposed$values,
#  			paste(indir, "/genotypes_eigenvalues",
#  				  sep=""),
#  			sep="\t", col.names=FALSE, row.names=FALSE)
#  
#  write.table(kinship_decomposed$vectors,
#  			paste(indir, "/genotypes_eigenvectors",
#  				  sep=""),
#  			sep="\t", col.names=FALSE, row.names=FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  indir <- "/homes/hannah/data/hmeyer/Supplementary/CEU.0908.impute.files"
#  
#  # specify directory to save data; if it doesn't exist yet, create i.e. needs
#  # writing permissions
#  datadir <- '/homes/hannah/tmp/phenotypeSimulator'
#  if (!dir.exists(datadir)) dir.create(datadir, recursive=TRUE)
#  
#  # specify filenames and parameters
#  totalGeneticVar <- 0.4
#  totalSNPeffect <- 0.1
#  h2s <- totalSNPeffect/totalGeneticVar
#  kinshipfile <- paste(indir, "/genotypes_genome_hapgen.controls.grm.rel",
#                                  sep="")
#  
#  genoFilePrefix <- paste(indir, "/genotypes_", sep="")
#  genoFileSuffix <- "_hapgen.controls.gen"
#  
#  # simulate phenotype with three phenotype components
#  simulation <- runSimulation(N = 1000, P = 3, cNrSNP=10, seed=43,
#                             kinshipfile = kinshipfile,
#                             oxgen = TRUE,
#                             genoFilePrefix = genoFilePrefix,
#                             genoFileSuffix = genoFileSuffix,
#                             chr = c(5,7,11),
#                             mBetaGenetic = 0, sdBetaGenetic = 0.2,
#                             theta=1,
#                             NrSNPsOnChromosome=c(533966, 503129, 428863),
#                             genVar = totalGeneticVar, h2s = h2s,
#                             phi = 0.6, delta = 0.2, rho=0.2,
#                             NrFixedEffects = 2, NrConfounders = c(2, 2),
#                             distConfounders = c("bin", "norm"),
#                             probConfounders = 0.2,
#                             genoDelimiter=" ",
#                             kinshipDelimiter="\t",
#                             kinshipHeader=FALSE,
#                             verbose = TRUE )

## ---- eval=FALSE, echo=FALSE---------------------------------------------
#  savedir <- paste(system.file("extdata", package="PhenotypeSimulator"),
#                              "/resultsSimulationAndLinearModel", sep="")
#  saveRDS(simulation$phenoComponentsFinal,
#          paste(savedir, "/simulation_phenoComponentsFinal.rds", sep=""))

## ---- echo=FALSE---------------------------------------------------------
savedir <- paste(system.file("extdata", package="PhenotypeSimulator"), 
                            "/resultsSimulationAndLinearModel", sep="")
simulation <- list(
    phenoComponentsFinal=readRDS(paste(savedir, 
                                       "/simulation_phenoComponentsFinal.rds", 
                            sep="")))

## ----heatmaps, fig.keep='all', fig.height=3, fig.width=7, eval=TRUE, fig.cap="\\label{fig:heatmaps}Heatmaps of the trait-by-trait correlation (Pearson correlation) of the simulated phenotype and its five phenotype components."----
cor_phenotype <- reshape2::melt(cor(simulation$phenoComponentsFinal$Y))
cor_phenotype$type <- "Y" 

cor_genBg <- reshape2::melt(cor(simulation$phenoComponentsFinal$Y_genBg))
cor_genBg$type <- "U" 

cor_genFixed <- reshape2::melt(cor(simulation$phenoComponentsFinal$Y_genFixed))          
cor_genFixed$type <- "XB" 

cor_noiseFixed <- reshape2::melt(cor(simulation$phenoComponentsFinal$Y_noiseFixed))          
cor_noiseFixed$type <- "WA" 

cor_noiseBg <- reshape2::melt(cor(simulation$phenoComponentsFinal$Y_noiseBg))          
cor_noiseBg$type <- "Psi" 

cor_correlatedBg <- reshape2::melt(cor(simulation$phenoComponentsFinal$Y_correlatedBg))          
cor_correlatedBg$type <- "T" 

cor_components <- rbind(cor_phenotype, cor_genBg, cor_genFixed, cor_noiseBg, 
                        cor_noiseFixed, cor_correlatedBg)
cor_components$type <- factor(cor_components$type, levels=c("Y", 
                                                            "WA",
                                                            "XB", 
                                                            "U", 
                                                            "T",
															"Psi"),
													labels=c("bold(Y)", 
                                                            "bold(WA)",
                                                            "bold(XB)", 
                                                            "bold(U)", 
                                                            "bold(T)",
															"bold(Psi)"))
colnames(cor_components) <- c("TraitA", "TraitB", "Correlation", "type")

# For prettier plotting
cor_components$TraitA <- gsub("_", " ", cor_components$TraitA)
cor_components$TraitB <- gsub("_", " ", cor_components$TraitB)

p_corr <- ggplot(data=cor_components, aes(x=TraitA, y=TraitB, fill=Correlation))  
p_corr <- p_corr + geom_tile() +
    facet_grid(~ type, labeller = label_parsed) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradient(limits = c(-1,1), 
                        guide=guide_colorbar(direction = "horizontal",
                                             title.position="top",
                                             title.hjust =0.5)) +
    xlab(" ") +
    ylab(" ") +
    theme_bw() +
    coord_fixed() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5),
			legend.position = "bottom",
			strip.background = element_rect(fill="white", color="white"),
            axis.text  =element_text(size=8),
            legend.title  =element_text(size=8),
            legend.text  =element_text(size=8),
            strip.text  =element_text(size=8))
print(p_corr)

## ---- echo=FALSE, eval=FALSE---------------------------------------------
#  ggsave(plot=p_corr,
#          filename=paste(savedir, "/correlation_phenotype.pdf", sep=""),
#          units="mm", width=86, height=60)

## ---- eval=FALSE---------------------------------------------------------
#  # Save phenotypes, kinship and non-genetic covariates in csv and
#  # GEMMA-specific format
#  outdirectory <- savePheno(simulation, directory = datadir,
#                              intercept_gemma=TRUE, format=c("csv", "gemma"),
#                              saveIntermediate=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  datadir="/homes/hannah/tmp/phenotypeSimulator"
#  
#  causalSNPs <- fread(paste(datadir, "/SNP_NrSNP10.csv", sep=""),
#  					header=TRUE, sep=",", stringsAsFactors=FALSE,
#  					data.table=FALSE)
#  causalSNPnames <- colnames(causalSNPs)
#  
#  LMM_mt <- fread(paste(datadir, "/GWAS_gemma_LMM_mt_pvalues.txt", sep=""),
#  				data.table=FALSE, header=FALSE, stringsAsFactors=FALSE)
#  LMM_st_trait1 <- fread(paste(datadir, "/GWAS_gemma_LMM_st_trait1_pvalues.txt",
#                              sep=""),
#  				data.table=FALSE, header=FALSE, stringsAsFactors=FALSE)
#  LMM_st_trait2 <- fread(paste(datadir, "/GWAS_gemma_LMM_st_trait2_pvalues.txt",
#                              sep=""),
#  				data.table=FALSE, header=FALSE, stringsAsFactors=FALSE)
#  LMM_st_trait3 <- fread(paste(datadir, "/GWAS_gemma_LMM_st_trait3_pvalues.txt",
#                              sep=""),
#  				data.table=FALSE, header=FALSE, stringsAsFactors=FALSE)
#  
#  LMM_mt <- LMM_mt[order(LMM_mt[,2]),]
#  LMM_mt$expected <- stats::ppoints(nrow(LMM_mt))
#  LMM_mt$model <- "mvLMM"
#  LMM_mt$trait <- "all"
#  colnames(LMM_mt)[1:2] <- c("rsID", "observed")
#  
#  LMM_st_trait1 <- LMM_st_trait1[order(LMM_st_trait1[,2]),]
#  LMM_st_trait1$expected <- LMM_mt$expected
#  LMM_st_trait1$model <- "uvLMM"
#  LMM_st_trait1$trait <- "trait1"
#  colnames(LMM_st_trait1)[1:2] <- c("rsID", "observed")
#  
#  LMM_st_trait2 <- LMM_st_trait2[order(LMM_st_trait2[,2]),]
#  LMM_st_trait2$expected <- LMM_mt$expected
#  LMM_st_trait2$model <- "uvLMM"
#  LMM_st_trait2$trait <- "trait2"
#  colnames(LMM_st_trait2)[1:2] <- c("rsID", "observed")
#  
#  LMM_st_trait3 <- LMM_st_trait3[order(LMM_st_trait3[,2]),]
#  LMM_st_trait3$expected <- LMM_mt$expected
#  LMM_st_trait3$model <- "uvLMM"
#  LMM_st_trait3$trait <- "trait3"
#  colnames(LMM_st_trait3)[1:2] <- c("rsID", "observed")
#  
#  pvalues <- rbind(LMM_mt, LMM_st_trait1, LMM_st_trait2, LMM_st_trait3)
#  pvalues$causal <- "all"
#  pvalues$causal[which(pvalues$rsID %in% causalSNPnames)] <- "simulated causal"
#  pvalues$causal <- factor(pvalues$causal, levels=c("all", "simulated causal"))
#  
#  p_qq <- ggplot(data=pvalues,
#  				aes(x=-log10(expected), y=-log10(observed)))
#  p_qq <- p_qq + geom_point(aes(color=causal, shape=trait)) +
#      scale_color_manual(values=c('darkgrey', '#1b9e77'),
#                         name="SNPs",
#  					guide=guide_legend(nrow = 2, title.position="top",
#                                         title.hjust =0.5,
#                                         override.aes = list(shape=21))) +
#      scale_shape_manual(values=c(18, 15:17),
#                         name="Traits",
#  					guide=guide_legend(nrow = 2,
#                                         override.aes = list(colour='darkgrey'),
#                                         title.position="top",
#                                         title.hjust =0.5)) +
#      geom_abline(intercept=0, slope=1, col="black") +
#  	geom_point(data=dplyr::filter(pvalues, causal == "simulated causal"),
#  				color='#1b9e77', aes(shape=trait)) +
#      xlab(expression(Expected~~-log[10](italic(p))) ) +
#      ylab(expression(Observed~~-log[10](italic(p)))) +
#      facet_grid(~ model) +
#      theme_bw() +
#  	theme(strip.background = element_rect(fill="white"),
#  			axis.title = element_text(size=12),
#  			axis.text  =element_text(size=10),
#  			legend.title  =element_text(size=12),
#  			legend.text  =element_text(size=10),
#  			strip.text  =element_text(size=10),
#  			legend.position = "bottom",
#              legend.direction = "horizontal")

## ----echo=FALSE, eval=FALSE, out.width='100%'----------------------------
#  ggsave(plot=p_qq,
#  		filename=paste(savedir, "/gwas_results_gemma_qqplot.png", sep=""),
#  		units="mm", width=86, height=100)

## ---- echo=FALSE, fig.cap="\\label{fig:qqplot}Quantile-quantile plots of p-values observed from the multivariate linear mixed model (mvLMM, traits:all) and the univariate linear mixed models (uvLMM, traits: trait1/trait2/trait3) fitted to each of about eight million genome-wide SNPs (grey), including the ten SNPs for which a phenotype effect was modelled (green)."----
knitr::include_graphics(paste(savedir, "/gwas_results_gemma_qqplot.png", sep=""))

## ---- eval=FALSE---------------------------------------------------------
#  sigSNPs <- dplyr::filter(pvalues, causal == "simulated causal",
#                           observed < 5*10^-8)
#  sigSNPs$observed_adjusted <- sigSNPs$observed
#  sigSNPs$observed_adjusted[which(sigSNPs$trait != "all")] <-
#      3*sigSNPs$observed[which(sigSNPs$trait != "all")]
#  sigSNPs <- dplyr::filter(sigSNPs, observed_adjusted < 5*10^-8)

## ---- echo=FALSE, eval=FALSE---------------------------------------------
#  write.table(sigSNPs,
#              paste(savedir, "/sigSNPs.csv", sep=""),
#              col.names=TRUE, row.names=FALSE, quote=FALSE, sep=",")

## ----sigSNPs, echo=FALSE-------------------------------------------------
                        
sigSNPs <- read.table(paste(savedir, "/sigSNPs.csv", sep=""), header=TRUE, 
                      sep=",")
knitr::kable(sigSNPs[, c(1,4,5,7)], 
             caption="\\label{tab:sigSNPs}P-values from multi-trait GWAS (model == mvLMM) and 
                        single-trait GWAS (model == uvLMM)",
             table.attr = "id=\"sigSNPs\"",
             digits=40)

## ---- eval=FALSE---------------------------------------------------------
#  # Read the SNPs simulated to have an effect on the phenotype and their SNP
#  # effects
#  SNPs <- read.csv(paste(datadir,"/SNP_NrSNP10.csv", sep=""),
#                   row.names=1)
#  SNPeffect <- read.csv(paste(datadir,"/SNP_effects_NrSNP10.csv", sep=""),
#                        row.names=1)
#  
#  
#  # HapGen simulated SNPs might be encoded as 7-25481146 , which R converts to
#  # X7.25481146 via make.names in read.table; substitute X and . to get original
#  # names
#  SNPnames <- gsub("\\.", "-", gsub("X", "", colnames(SNPs)))
#  SNPnames <- SNPnames[order(SNPnames)]
#  simulated_withEffect <- dplyr::filter(pvalues, rsID %in% SNPnames)

## ---- echo=FALSE, eval=FALSE---------------------------------------------
#  write.table(simulated_withEffect, paste(savedir,"/simulated_withEffect.csv", sep=""), col.names=TRUE, quote=FALSE,row.names=FALSE)

## ---- echo=FALSE, eval=TRUE----------------------------------------------
SNPeffect <- read.csv(paste(savedir,"/SNP_effects_NrSNP10.csv", sep=""),
                      row.names=1)
SNPs <- read.csv(paste(savedir,"/SNP_NrSNP10.csv", sep=""),
                 row.names=1)
SNPnames <- gsub("\\.", "-", gsub("X", "", colnames(SNPs)))
SNPnames <- SNPnames[order(SNPnames)]
simulated_withEffect <- read.csv(paste(savedir,"/simulated_withEffect.csv", sep=""))

## ----freq-beta, fig.cap="\\label{fig:freq-beta}P-values, allele frequencies and simulated effect sizes of the ten SNPs simulated to have an effect on phenotype."----
# get the allele frequencies of the SNPs simulated to have an effect on the 
# phenotype
freq <- apply(SNPs, 2, getAlleleFrequencies)
minor <- apply(freq, 2, min)
minor <- minor[order(gsub("X", "", names(minor)))]

# format the effect size table and combine with p-value and frequency
# information
effects <- reshape2::melt(SNPeffect, value.name="beta", variable.name="name")
effects$type <- gsub('\\..*', '', effects$name)
effects$rsID <- gsub("\\.", "-", gsub('.*_', '', effects$name))
effects$trait <- rep(paste("trait", 1:3, sep=""), 10)
effects <- effects[order(effects$rsID),]
effects <- dplyr::filter(effects, rsID %in% simulated_withEffect$rsID)

simulated_withEffect <- simulated_withEffect[order(simulated_withEffect[, 1]),]
simulated_withEffect$beta <- NA
simulated_withEffect$beta[simulated_withEffect$trait != "all"] <- effects$beta
simulated_withEffect$freq <- rep(minor[which(SNPnames %in% effects$rsID)], each=4)

# plot association p-values in relation to the absolute value of the simulated 
# effect sizes and the SNPs allele frequencies
p <- ggplot(dplyr::filter(simulated_withEffect, trait != "all"))
p + geom_point(aes(x=freq, y=-log10(observed), color=abs(beta), shape=trait)) +
    geom_hline(yintercept=-log10(5*10^(-8)), color="darkgrey") +
    scale_color_gradientn(colours=c('#f0f9e8','#ccebc5','#a8ddb5','#7bccc4',
                                    '#4eb3d3','#2b8cbe','#08589e'),
                          name=expression("| simulated" ~beta~ "|")) +
    scale_shape_discrete(name="Phenotype in GWAS") +
    xlab("Allele frequencies") +
    ylab(expression(-log[10](p-value))) +
    theme_bw()

