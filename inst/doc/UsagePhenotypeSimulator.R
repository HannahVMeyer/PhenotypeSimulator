## ---- echo = FALSE, message=FALSE----------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
##devtools::load_all()
library("PhenotypeSimulator")

## ------------------------------------------------------------------------
# Set parameters
genVar <- 0.4
noiseVar <- 1- genVar
shared <- 0.6
independent <- 1 - shared

# simulate simple bi-allelic genotypes and estimate kinship
genotypes <- simulateGenotypes(N = 100, NrSNP = 10000, 
                               frequencies = c(0.05, 0.1, 0.3, 0.4), 
                               verbose = FALSE)
genotypes_sd <- standardiseGenotypes(genotypes$genotypes)
kinship <- getKinship(N=N, X=genotypes_sd, verbose = FALSE)

# simulate phenotype components
genBg <- geneticBgEffects(N = 100, kinship = kinship, P = 15)
noiseBg <- noiseBgEffects(N = 100, P = 15)

# rescale phenotype components
genBg_shared_scaled <- rescaleVariance(genBg$shared, shared * genVar)
genBg_independent_scaled <- rescaleVariance(genBg$independent, 
                                            independent * genVar)
noiseBg_shared_scaled <- rescaleVariance(noiseBg$shared, shared * noiseVar)
noiseBg_independent_scaled <- rescaleVariance(noiseBg$independent, 
                                              independent *noiseVar)

# Total variance proportion shave to add up yo 1
total <- independent * noiseVar + independent * genVar + 
    shared * noiseVar + shared * genVar

total == 1

# combine components into final phenotype
Y <- scale(genBg_shared_scaled$component + noiseBg_shared_scaled$component +
    genBg_independent_scaled$component + noiseBg_independent_scaled$component)

## ------------------------------------------------------------------------
# simulate phenotype with population structure and observational noise effects 
# only

# genetic variance
genVar <- 0.4

# random genetic variance: h2b 
phenotype <- runSimulation(N = 100, P = 15,  tNrSNP = 10000, 
                           SNPfrequencies = c(0.05, 0.1,0.3,0.4), 
                           genVar = 0.4, h2bg = 1, phi = 1, 
                           verbose = TRUE)

## ------------------------------------------------------------------------
# read genotypes from external file
# use one of the sample genotype file provided in the 
# extdata/genotypes/subfolders (e.g.extdata/genotypes/hapgen )
genotypefile <- system.file("extdata/genotypes/hapgen", 
                            "genotypes_hapgen.controls.gen", 
                            package = "PhenotypeSimulator")
# remove the .gen ending (oxgen specific endings .gen and .sample are added 
# automatically )
genotypefile <- gsub("\\.gen","", genotypefile)

genotypes <- readStandardGenotypes(N=100, filename = genotypefile,
                                    format="oxgen",
                                    delimiter = ",", verbose=TRUE)

genotypes_sd <- standardiseGenotypes(genotypes$genotypes)

# kinship estimate based on standardised SNPs 
kinship <- getKinship(N=100, X=genotypes_sd, verbose = FALSE)

# simulate 30 genetic variant effects (from non-standardised SNP genotypes)
causalSNPs <- getCausalSNPs(N=100, genotypes = genotypes$genotypes, 
                            NrCausalSNPs = 30, verbose = FALSE)
genFixed <- geneticFixedEffects(N = 100, P = 15, X_causal = causalSNPs)  

# simulate infinitesimal genetic effects
genBg <- geneticBgEffects(N=100, kinship = kinship, P = 15)

# simulate 4 different confounder effects:
# * 1 binomial covariate effect shared across all traits
# * 2 categorical (3 categories) independent covariate effects
# * 1 categorical (4 categories) independent  covariate effect
# * 2 normally distributed independent and shared covariate effects
noiseFixed <- noiseFixedEffects(N = 100, P = 15, NrFixedEffects = 4, 
                                NrConfounders = c(1, 2, 1, 2),
                                pIndependentConfounders = c(0, 1, 1, 0.5),  
                                distConfounders = c("bin", "cat_norm", 
                                                    "cat_unif", "norm"),
                                probConfounders = 0.2, 
                                catConfounders = c(0, 3, 4, 0))

# simulate correlated effects with max correlation of 0.8
correlatedBg <- correlatedBgEffects(N = 100, P = 15, pcorr = 0.8)

# simulate observational noise effects
noiseBg <- noiseBgEffects(N = 100, P = 15)

# total SNP effect on phenotype: 0.01
genVar <- 0.6
noiseVar <- 1 - genVar
totalSNPeffect <- 0.01
h2s <- totalSNPeffect/genVar
phi <- 0.6 
rho <- 0.1
delta <- 0.3
shared <- 0.8
independent <- 1 - shared

# rescale phenotype components
genFixed_shared_scaled <- rescaleVariance(genFixed$shared, shared * h2s *genVar)
genFixed_independent_scaled <- rescaleVariance(genFixed$independent, 
                                            independent * h2s *genVar)
genBg_shared_scaled <- rescaleVariance(genBg$shared, shared * (1-h2s) *genVar)
genBg_independent_scaled <- rescaleVariance(genBg$independent, 
                                            independent * (1-h2s) * genVar)

noiseBg_shared_scaled <- rescaleVariance(noiseBg$shared, shared * phi* noiseVar)
noiseBg_independent_scaled <- rescaleVariance(noiseBg$independent, 
                                              independent * phi* noiseVar)
correlatedBg_scaled <- rescaleVariance(correlatedBg$correlatedBg, 
                                       shared * rho * noiseVar)

noiseFixed_shared_scaled <- rescaleVariance(noiseFixed$shared, shared * delta * 
                                                noiseVar)
noiseFixed_independent_scaled <- rescaleVariance(noiseFixed$independent, 
                                              independent * delta * noiseVar)

# Total variance proportions have to add up yo 1
total <- shared * h2s *genVar +  independent * h2s * genVar +
    shared * (1-h2s) * genVar +   independent * (1-h2s) * genVar +
    shared * phi* noiseVar +  independent * phi* noiseVar +
    rho * noiseVar +
    shared * delta * noiseVar +  independent * delta * noiseVar

total == 1

# combine components into final phenotype
Y <- scale(genBg_shared_scaled$component + noiseBg_shared_scaled$component +
        genBg_independent_scaled$component + noiseBg_independent_scaled$component +
        genFixed_shared_scaled$component + noiseFixed_shared_scaled$component +
        genFixed_independent_scaled$component + noiseFixed_independent_scaled$component +
        correlatedBg_scaled$component)


## ---- tidy=TRUE, tidy.opts = list(width.cutoff = 60)---------------------
# simulate phenotype with the same five phenotype components and settings as 
# above; display progress via verbose=TRUE
phenotype <- runSimulation(N = 100, P = 15, genotypefile=genotypefile, 
                           format ="oxgen",
                           cNrSNP=30, 
                           genVar = genVar, 
                           h2s = h2s, phi = 0.6, delta = 0.3,
                           distBetaGenetic = "unif", mBetaGenetic = 0.5, 
                           sdBetaGenetic = 1,
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
image(t(phenotype$phenoComponentsFinal$Y), main="Phenotype: [samples x traits]", 
      axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)
image(cor(phenotype$phenoComponentsFinal$Y), 
      main="Correlation of traits [traits x traits]", axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Traits", line = 1)
mtext(side = 2, text = "Traits", line = 1)

## ---- eval = FALSE-------------------------------------------------------
#  out <- savePheno(phenotype, directory="/tmp",
#                   outstring="test_simulation",
#                   format=c("csv", "plink"), verbose=FALSE)

## ------------------------------------------------------------------------
## a) Draw cuasal SNPs from a simulated genotype matrix
# simulate 10,000 bi-allelic SNP genotypes for 100 samples with randomly drawn 
# allele frequencies of 0.05, 0.1, 0.3 and 0.4. 
genotypes <- simulateGenotypes(N = 100, NrSNP = 10000, 
                               frequencies = c(0.05, 0.1, 0.3,0.4), 
                               verbose = FALSE)

# draw 10 causal SNPs from the genotype matrix (use non-standardised allele 
# codes i.e. (0,1,2))
causalSNPs <- getCausalSNPs(N=100, NrCausalSNPs = 10, 
                            genotypes = genotypes$genotypes)

## ------------------------------------------------------------------------
## b) Draw SNPs from external genotype files: 
# read genotypes from external file
# use one of the sample genotype file provided in the 
# extdata/genotypes/subfolders (e.g.extdata/genotypes/hapgen )
genotypefile <- system.file("extdata/genotypes/hapgen", 
                            "genotypes_hapgen.controls.gen", 
                            package = "PhenotypeSimulator")
# remove the .gen ending (oxbgen specific endings .gen and .sample are added 
# automatically )
genotypefile <- gsub("\\.gen","", genotypefile)

genotypes <- readStandardGenotypes(N = 100, filename = genotypefile,
                                    format="oxgen",
                                    delimiter = ",", verbose=TRUE)

causalSNPsFromFile <- getCausalSNPs(N = 100, NrCausalSNPs = 10, 
                                    genotypes=genotypes$genotypes)
                                    

## ------------------------------------------------------------------------
## c) draw 10 causal SNPs from external genotype files: sample 10 SNPs from 
## chromosome 22
# use sample genotype file provided in the extdata/genotypes folder
genotypeFile <- system.file("extdata/genotypes/", "genotypes_chr22.csv", 
                            package = "PhenotypeSimulator")
genoFilePrefix <- gsub("chr.*", "", genotypeFile) 
genoFileSuffix <- ".csv" 

causalSNPsSampledFromFile <- getCausalSNPs(N = 10, NrCausalSNPs = 10, chr = 22, 
                                    genoFilePrefix = genoFilePrefix, 
                                    genoFileSuffix = genoFileSuffix,  
                                    delimiter = ",", verbose=TRUE)

## ---- fig.show='hold'----------------------------------------------------
# create genetic variant effects with 20% of SNPs having a specific effect, 
# affecting 40% of all simulated traits
fixedGenetic <- geneticFixedEffects(X_causal = causalSNPsFromFile, 
                                    N = 100, P = 10, 
                                    pIndependentGenetic = 0.2, 
                                    pTraitIndependentGenetic = 0.4)

## ---- fig.show='hold', echo=FALSE, fig.height=3.4, fig.width=3.4---------
image(fixedGenetic$shared, main="Shared genetic variant effects", axes=FALSE, 
      cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)
image(fixedGenetic$independent, main="Independent genetic variant effects", 
      axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)


## ---- fig.show='hold'----------------------------------------------------
## a) Estimate kinship from simulated genotypes
genotypes <- simulateGenotypes(N = 100, NrSNP = 10000, 
                               frequencies = c(0.05, 0.1, 0.3,0.4), 
                               verbose = FALSE)

genotypes_sd <- standardiseGenotypes(genotypes$genotypes)
kinship <- getKinship(N=100, X=genotypes_sd, verbose = FALSE)

## b) Read kinship from external kinship file
kinshipFile <- system.file("extdata/kinship/", "kinship.csv", 
                           package = "PhenotypeSimulator")
kinshipFromFile <- getKinship(N=50, kinshipfile = kinshipFile, 
                              verbose = FALSE)

genBg <- geneticBgEffects(N=100, kinship = kinship, P = 15)

## ---- fig.show='hold', echo=FALSE, fig.height=3.4, fig.width=3.4---------
image(genBg$shared, main="Shared infinitesimal genetic effects", axes=FALSE, 
      cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)
image(genBg$independent, main="Independent infinitesimal genetic effects",  
      axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)

## ---- fig.show='hold'----------------------------------------------------
# create 1 non-genetic covariate effect affecting 30% of all simulated traits. The effect 
# follows a uniform distribution between 30 and 40  (resembling for instance age 
# in a study cohort).
fixedNoiseUniform <- noiseFixedEffects(N = 100, P = 10, NrConfounders = 1, 
                                       pIndependentConfounders = 1, 
                                       pTraitIndependentConfounders = 0.3, 
                                       distConfounders = "unif", 
                                       mConfounders = 35, sdConfounders = 5)

# create 2 non-genetic covariate effects with 1 specific confounder affecting 
# 20% of all simulated traits. The effects follow a normal distribution
fixedNoiseNormal <- noiseFixedEffects(N = 100, P = 10, NrConfounders = 2, 
                                      pIndependentConfounders = 0.5, 
                                      pTraitIndependentConfounders = 0.2, 
                                      distConfounders = "norm", 
                                      mConfounders = 0, sdConfounders = 1)

# create 1 non-genetic covariate effects affecting  all simulated traits. The 
# effect follows a binomial distribution with probability 0.5 (resembling for
# instance sex in a study cohort).
fixedNoiseBinomial <- noiseFixedEffects(N = 100, P = 10, NrConfounders = 1, 
                                        pIndependentConfounders = 0, 
                                        distConfounders = "bin", 
                                        probConfounders = 0.5)

## ---- fig.show='hold', echo=FALSE, fig.height=3.5, fig.width=4-----------
image(fixedNoiseUniform$independent, 
      main="Independent non-genetic covariate effects\n(uniform confounder dist)", 
      axes=FALSE,  cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)

image(fixedNoiseNormal$shared, 
      main="Shared non-genetic covariate effects\n(normal confounder dist)", 
      axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)

image(fixedNoiseNormal$independent, 
      main="Independent non-genetic covariate effects\n(normal confounder dist)", 
      axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)

image(fixedNoiseBinomial$shared, 
     main="Shared non-genetic covariate effects\n(binomial confounder dist)",  
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
pairs(correlatedNoise$correlatedBg, pch = ".", 
      main=paste("Correlation at furthest distance to diagonal:\n",
                 furthestDistCorr), cex.main=0.8)
image(correlatedNoise$correlatedBg, main="Correlated effects",  axes=FALSE, 
      cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)

## ----fig.show='hold'-----------------------------------------------------
# simulate a noise random effect for 10 traits
noiseBg <- noiseBgEffects(N = 100, P = 10, mean = 0, sd = 1)

## ---- fig.show='hold', echo=FALSE, fig.height=3.4, fig.width=3.4---------
image(noiseBg$shared, main="Shared observational noise effects", axes=FALSE, 
      cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)
image(noiseBg$independent, main="Independent observational noise effects", 
      axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)

