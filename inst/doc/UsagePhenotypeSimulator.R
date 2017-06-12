## ---- echo = FALSE, message=FALSE----------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library("PhenotypeSimulator")

## ------------------------------------------------------------------------
### step-by-step

# set genetic and noise models
modelGenetic <- "geneticBgOnly"
modelNoise <- "noiseBgOnly"

# simulate genotypes and estimate kinship
genotypes <- simulateGenotypes(N = 100, NrSNP = 10000, 
                               frequencies = c(0.05, 0.1, 0.3, 0.4), 
                               verbose = FALSE)
kinship <- getKinship(genotypes$X_sd, norm=TRUE, verbose = FALSE)

# simulate phenotype components
genBg <- geneticBgEffects(kinship = kinship, P = 15)
noiseBg <- noiseBgEffects(N = 100, P = 15)

# combine components into final phenotype with genetic variance component 
# explaining 40% of total variance
phenotype <- createPheno(N = 100, P = 15, noiseBg = noiseBg, genBg = genBg, 
                         modelNoise = modelNoise, modelGenetic = modelGenetic, 
                         genVar = 0.4, verbose = FALSE)

## ------------------------------------------------------------------------
### via `runSimulation`

# simulate phenotype with genetic and noise random effects only
# genetic variance
genVar <- 0.4
# random genetic variance: h2b 
phenotype <- runSimulation(N = 100, P = 15,  tNrSNP = 10000, 
                           SNPfrequencies = c(0.05, 0.1,0.3,0.4), 
                           normalise = TRUE, genVar = 0.4, h2bg = 1, phi = 1, 
                           verbose = FALSE)

## ------------------------------------------------------------------------
### step-by- step

# set genetic and noise models
modelGenetic <- "geneticFixedAndBg"
modelNoise <- "noiseFixedAndBgAndCorrelated"

# simulate genotypes and estimate kinship
genotypes <- simulateGenotypes(N = 100, NrSNP = 10000, 
                               frequencies = c(0.05, 0.1,0.3,0.4), 
                               verbose = FALSE)
# kinship estimate based on standardised SNPs (as described in )
kinship <- getKinship(X=genotypes$X_sd, norm=TRUE, verbose = FALSE)

# simulate 30 fixed genetic effects (from non-standardised SNP genotypes)
causalSNPs <- getCausalSNPs(genotypes = genotypes, NrCausalSNPs = 30, 
                            standardise = FALSE, verbose = FALSE)
genFixed <- geneticFixedEffects(N = 100, P = 15, X_causal = causalSNPs)  

# simulate random genetic effects
genBg <- geneticBgEffects(kinship = kinship, P = 15)

# simulate 4 fixed noise effects:
# * 1 binomial fixed noise effect shared across all traits
# * 2 categorical (3 categories) independent fixed noise traits
# * 1 categorical (4 categories) independent fixed noise traits
# * 2 normally distributed independent and shared fixed noise traits
noiseFixed <- noiseFixedEffects(N = 100, P = 15, NrFixedEffects = 4, 
                                NrConfounders = c(1, 2, 1, 2),
                                pIndependentConfounders = c(0, 1, 1, 0.5),  
                                distConfounders = c("bin", "cat_norm", 
                                                    "cat_unif", "norm"),
                                probConfounders = 0.2, 
                                catConfounders = c(0, 3, 4, 0))

# simulate correlated noise effects with max correlation of 0.8
correlatedBg <- correlatedBgEffects(N = 100, P = 15, pcorr = 0.8)

# simulate random noise effects
noiseBg <- noiseBgEffects(N = 100, P = 15)

# total SNP effect on phenotype: 0.01
totalGeneticVar <- 0.4
totalSNPeffect <- 0.01
h2s <- totalSNPeffect/totalGeneticVar

# combine components into final phenotype with genetic variance component 
# explaining 40% of total variance
phenotype <- createPheno(N = 100, P = 15, noiseBg = noiseBg, 
                         noiseFixed = noiseFixed, correlatedBg = correlatedBg, 
                         genFixed = genFixed, genBg = genBg, 
                         modelNoise = modelNoise, modelGenetic = modelGenetic, 
                         genVar = totalGeneticVar, h2s = h2s, phi = 0.6, 
                         rho = 0.1, delta = 0.3, gamma = 1,  verbose = FALSE)

## ---- tidy=TRUE, tidy.opts = list(width.cutoff = 60)---------------------
### via `runSimulation`

# simulate phenotype with the same five phenotype components and settings as 
# above; display progress via verbose=TRUE
phenotype <- runSimulation(N = 100, P = 15,  tNrSNP = 10000,  cNrSNP=30, 
                           SNPfrequencies = c(0.05, 0.1,0.3,0.4), 
                           normalise = TRUE, genVar = totalGeneticVar, 
                           h2s = h2s, phi = 0.6, delta = 0.3, gamma = 1,
                           NrFixedEffects = 4, NrConfounders = c(1, 2, 1, 2),
                           pIndependentConfounders = c(0, 1, 1, 0.5),  
                           distConfounders = c("bin", "cat_norm", 
                                               "cat_unif", "norm"), 
                           probConfounders = 0.2, 
                           catConfounders = c(0, 3, 4, 0),
                           pcorr = 0.8,
                           verbose = TRUE )


## ---- echo=FALSE---------------------------------------------------------
# show proportion of variance of the different phenotype components in the 
# final phenotype
varGenFixed <- t(phenotype$varComponents
                 [grepl("var_genFix", names(phenotype$varComponents))])
varGenBg <- t(phenotype$varComponents
              [grepl("var_genBg", names(phenotype$varComponents))])

varNoiseFixed <- t(phenotype$varComponents
                   [grepl("var_noiseFixed", names(phenotype$varComponents))])
varNoiseBg <- t(phenotype$varComponents
                [grepl("var_noiseBg", names(phenotype$varComponents))])
varNoiseCorr <- t(phenotype$varComponents
                  [grepl("var_noiseCor", names(phenotype$varComponents))])

totalPropVariance <- as.matrix(t(data.frame(varGenFixed, 
                                            varGenBg, 
                                            varNoiseFixed,
                                            varNoiseBg, 
                                            varNoiseCorr=c(varNoiseCorr, 0))))
totalPropVariance <- cbind(totalPropVariance, rowSums(totalPropVariance))
totalPropVariance <- rbind(totalPropVariance, 
                           sumProportion=colSums(totalPropVariance))

colnames(totalPropVariance) <- c("shared effect", "independent effect", 
                                 "total effect")

knitr::kable(totalPropVariance, caption="Proportion of variance explained
             by the different phenotype components")

## ---- fig.show='hold', fig.height=3.4, fig.width=3.4---------------------
### show 'image' of phenotype and correlation between phenotypic traits
image(t(phenotype$phenoComponents$Y), main="Phenotype: [samples x traits]", 
      axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)
image(cor(phenotype$phenoComponents$Y), 
      main="Correlation of traits [traits x traits]", axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Traits", line = 1)
mtext(side = 2, text = "Traits", line = 1)

## ---- eval = FALSE-------------------------------------------------------
#  out <- savePheno(phenotype, directoryGeno="/tmp/genotypes",
#            directoryPheno="/tmp/phenotypes", outstring="test_simulation",
#            sample_subset_vec = c(50, 70), pheno_subset_vec = c(5, 10),
#            saveAsTable=TRUE, saveAsPlink=TRUE, verbose=FALSE)

## ------------------------------------------------------------------------
## a) Draw cuasal SNPs from a simulated genotype matrix
# simulate 10,000 bi-allelic SNP genotypes for 100 samples with randomly drawn 
# allele frequencies of 0.05, 0.1, 0.3 and 0.4. 
genotypes <- simulateGenotypes(N = 100, NrSNP = 10000, 
                               frequencies = c(0.05, 0.1, 0.3,0.4), 
                               verbose = FALSE)

# draw 10 causal SNPs from the genotype matrix (use non-standardised allele 
# codes i.e. (0,1,2))
causalSNPs <- getCausalSNPs(NrCausalSNPs = 10, genotypes = genotypes, 
                            standardise = FALSE)

## ------------------------------------------------------------------------
## b) draw 10 causal SNPs from external genotype files: sample 10 SNPs from 
## chromsome 22
# use sample genotype file provided in the extdata/genotypes folder
genotypeFile <- system.file("extdata/genotypes/", "genotypes_chr22.csv", 
                            package = "PhenotypeSimulator")
genoFilePrefix <- gsub("chr.*", "", genotypeFile) 
genoFileSuffix <- ".csv" 

causalSNPsFromFile <- getCausalSNPs(NrCausalSNPs = 10, chr = 22, 
                                    genoFilePrefix = genoFilePrefix, 
                                    genoFileSuffix = genoFileSuffix,  
                                    standardise = FALSE, 
                                    genoFileDelimiter = ",", verbose=FALSE)

## ---- fig.show='hold'----------------------------------------------------
# create genetic fixed effects with 20% of SNPs having a specific effect, 
# affecting 40% of all simulated traits
fixedGenetic <- geneticFixedEffects(X_causal = causalSNPs, N = 100, P = 10, 
                                    pIndependentGenetic = 0.2, 
                                    pTraitIndependentGenetic = 0.4)

## ---- fig.show='hold', echo=FALSE, fig.height=3.4, fig.width=3.4---------
image(fixedGenetic$shared, main="Shared fixed genetic effects", axes=FALSE, 
      cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)
image(fixedGenetic$independent, main="Independent fixed genetic effects", 
      axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)


## ---- fig.show='hold'----------------------------------------------------
## a) Estimate kinship from simulated genotypes
kinship <- getKinship(genotypes$X_sd, norm=TRUE, verbose = FALSE)

## b) Read kinship from external kinship file
kinshipFile <- system.file("extdata/kinship/", "kinship.csv", 
                           package = "PhenotypeSimulator")
kinshipFromFile <- getKinship(kinshipfile = kinshipFile, norm=FALSE, 
                              verbose = FALSE)

genBg <- geneticBgEffects(kinship = kinship, P = 15)

## ---- fig.show='hold', echo=FALSE, fig.height=3.4, fig.width=3.4---------
image(genBg$shared, main="Shared random genetic effects", axes=FALSE, 
      cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)
image(genBg$independent, main="Independent random genetic effects",  
      axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)

## ---- fig.show='hold'----------------------------------------------------
# create 1 noise fixed effects affecting 30% of all simulated traits. The effect 
# follows a uniform distribution between 30 and 40  (resembling for instance age 
# in a study cohort).
fixedNoiseUniform <- noiseFixedEffects(N = 100, P = 10, NrConfounders = 1, 
                                       pIndependentConfounders = 1, 
                                       pTraitIndependentConfounders = 0.3, 
                                       distConfounders = "unif", 
                                       mConfounders = 35, sdConfounders = 5)

# create 2 noise fixed effects with 1 specific confounder affecting 20% of all 
# simulated traits. The effects follow a normal distribution
fixedNoiseNormal <- noiseFixedEffects(N = 100, P = 10, NrConfounders = 2, 
                                      pIndependentConfounders = 0.5, 
                                      pTraitIndependentConfounders = 0.2, 
                                      distConfounders = "norm", 
                                      mConfounders = 0, sdConfounders = 1)

# create 1 noise fixed effects affecting  all simulated traits. The effect 
# follows a binomial distribution with probability 0.5 (resembling for instance 
# sex in a study cohort).
fixedNoiseBinomial <- noiseFixedEffects(N = 100, P = 10, NrConfounders = 1, 
                                        pIndependentConfounders = 0, 
                                        distConfounders = "bin", 
                                        probConfounders = 0.5)

## ---- fig.show='hold', echo=FALSE, fig.height=3.5, fig.width=4-----------
image(fixedNoiseUniform$independent, 
      main="Independent fixed noise effects\n(uniform confounder dist)", 
      axes=FALSE,  cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)

image(fixedNoiseNormal$shared, 
      main="Shared fixed noise effects\n(normal confounder dist)", 
      axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)

image(fixedNoiseNormal$independent, 
      main="Independent fixed noise effects\n(normal confounder dist)", 
      axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)

image(fixedNoiseBinomial$shared, 
     main="Shared fixed noise effects\n(binomial confounder dist)",  
     axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)

## ---- fig.show='hold', fig.height=3.4, fig.width=3.4---------------------
# simulate correlated noise effect for 10 traits with top-level 
# correlation of 0.8
correlatedNoise <- correlatedBgEffects(N = 100, P = 10, pcorr = 0.8 )

# correlation structure of the traits: strong the closer to the diagonal, 
# little correlation at the furthest distance to the diagonal 
furthestDistCorr <- 0.4^(10-1)
pairs(correlatedNoise, pch = ".", 
      main=paste("Correlation at furthest distance to diagonal:\n",
                 furthestDistCorr), cex.main=0.8)
image(correlatedNoise, main="Correlated noise effects",  axes=FALSE, 
      cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)

## ----fig.show='hold'-----------------------------------------------------
# simulate a noise random effect for 10 traits
noiseBg <- noiseBgEffects(N = 100, P = 10, mean = 0, sd = 1)

## ---- fig.show='hold', echo=FALSE, fig.height=3.4, fig.width=3.4---------
image(noiseBg$shared, main="Shared random noise effects", axes=FALSE, 
      cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)
image(noiseBg$independent, main="Independent random noise effects", 
      axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)

