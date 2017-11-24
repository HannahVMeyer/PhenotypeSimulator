#' Scale phenotype component.
#'
#' The function scales the specified phenotype component such that the average 
#' column variance is equal
#' to the user-specified proportion of variance. 
#'
#' @param component numeric [N x P] phenotype matrix where N are the number of 
#' observations and P numer of phenotypes
#' @param propvar numeric specifying the proportion of variance that should be 
#' explained by this phenotype component
#' @return component_scaled numeric [N x P] phenotype matrix if propvar != 0 or 
#' NULL
#' @export
#' @examples
#' x <- matrix(rnorm(100), nc=10)
#' x_scaled <- rescaleVariance(x, propvar=0.4)
rescaleVariance <- function(component, propvar) {
    if (propvar != 0) {
        var_component <- var(component)
        mean_var <- mean(diag(var_component))
        scale_factor <- mean_var/propvar
        component_scaled <- component/sqrt(scale_factor)
        return(component_scaled)
    } else {
        return(NULL)
    }
}


#' Set simulation model.
#'
#' Based on parameters provided, this function sets the name for the phenotype 
#' simulation. The model name is needed downstream for phenotype component 
#' simulations, but can also be set manually.
#'
#' @param genVar Total genetic variance [double] 
#' @param h2s Proportion [double] of variance of fixed genetic effects
#' @param h2bg Proportion [double] of variance of random genetic effects
#' @param theta Proportion [double] of variance of shared fixed genetic effects
#' @param eta Proportion [double] of variance of shared bg genetic effects
#' @param noiseVar Total genetic variance [double] 
#' @param rho Proportion [double] of variance of correlated noise effects
#' @param pcorr Correlation [double] between phenotypes
#' @param delta Proportion [double] of fixed noise variance
#' @param gamma Proportion [double] of variance of shared fixed noise effects
#' @param phi Proportion [double] of variance of background noise effects
#' @param alpha Proportion [double] of Variance of shared bg noise effect
#' @param pIndependentConfounders Proportion [double] of noise effects 
#' (confounders) to have a trait-independent effect
#' @param pTraitIndependentConfounders Proportion [double] of traits influenced  
#' by independent fixed noise effects
#' @param pIndependentGenetic Proportion [double] of genetic effects (SNPs) to 
#' have a trait-independent fixed effect
#' @param pTraitIndependentGenetic Proportion [double] of traits influenced by 
#' independent fixed genetic effects
#' @param v [boolean]; if TRUE, progress info is printed to standard out
#' @return named list containing the genetic model (modelGenetic) and the noise 
#' model (modelNoise). Options are:
#' modelNoise: "nNoise", "noiseFixedOnly", "noiseBgOnly", "noiseCorrelatedOnly",
#'  "noiseFixedAndBg","noiseCorrelatedAndBg", "noiseFixedAndCorrelated",
#'  "noiseFixedAndBgAndCorrelated"
#' modelGenetic: "noGenetic","geneticBgOnly", "geneticFixedOnly",
#' "geneticFixedAndBg"
#' @export
#' @examples
#' #genetic fixed effects only
#' model <- setModel(genVar=1, h2s=1)
#' 
#' #genetic fixed and bg effects
#' model <- setModel(genVar=1, h2s=0.01)
#' 
#' #genetic and noise fixed effects only
#' model <- setModel(genVar=0.4, h2s=1, delta=1)
setModel <- function(genVar=NULL, h2s=NULL, theta=0.8, h2bg=NULL, eta=0.8, 
                     noiseVar=NULL, delta=NULL, gamma=0.8, rho=NULL, phi=NULL, 
                     alpha=0.8, pcorr=0.6, pIndependentConfounders=0.4,  
                     pTraitIndependentConfounders=0.2,  pIndependentGenetic=0.4, 
                     pTraitIndependentGenetic=0.2, v=TRUE)  {
    if (is.null(c(genVar, noiseVar, h2bg, h2s, delta, rho, phi))) {
        stop("No variance components specified")
    }
    if (is.null(genVar)) {
        if (is.null(noiseVar)) {
            stop(paste("Neither genVar nor noiseVar are provided, thus",
                      "proportion of variance from genetics cannot be deduced")
                 )
        } else {
            genVar <- 1 - noiseVar
        } 
    }
    if (is.null(noiseVar)) {
        noiseVar <- 1 - genVar
    }
    if (noiseVar + genVar > 1) {
        stop("Sum of genetic variance and noise variance is greater than 1")
    }
    if (noiseVar == 0 && any(!is.null(phi), !is.null(rho), !is.null(delta))) {
        stop(paste("The noise variance is set to 0 (or genetic variance set to",
                    " 1) but noise variance components are supplied")
        )
    }
    if (genVar == 0 && any(!is.null(h2s), !is.null(h2bg))) {
        stop(paste("The genetic variance is set to 0 (or noise variance set to",
                   " 1) but genetic variance components are supplied")
        )
    }
    vmessage(c("The total noise variance (noiseVar) is:", noiseVar), verbose=v)
    if ((noiseVar) == 0 ) {
        modelNoise="noNoise"
        vmessage(c("The noise model is:", modelNoise), verbose=v)
    } else {
        if (all(c(is.null(delta), is.null(rho), is.null(phi)))) {
            stop(paste("Neither delta (fixed noise effect variance) nor rho",
                       "(correlated noise effect variance) or phi (random",
                       "noise effect variance) are provided, at least", 
                       "one is required")
            )
        } 
        if (length(c(delta, rho, phi)) >=2) {
            if (sum(c(delta, rho, phi)) > 1 ) {
                stop(paste("Sum of the proportion of the variance of noise",
                           "effects is greater than 1; change noiseVar, delta",
                           "(fixed noise effect variance), rho (correlated",
                           " noise effect variance) or phi (random noise effect"
                           , "variance) such that delta + rho + phi = 1")
                )
            }
            if (length(c(delta, rho, phi)) == 3 && 
                sum(c(delta, rho, phi)) < 1) {
                stop(paste("Sum of the proportion of the variance of noise",
                           "effects is less than 1; change noiseVar, delta",
                           "(fixed noise effect variance), rho (correlated",
                           " noise effect variance) or phi (random noise effect"
                           , "variance) such that delta + rho + phi = 1")
                )
            }
        } 
        if (is.null(phi) && all(c(!is.null(delta), !is.null(rho)))) { 
            phi <- 1 - delta - rho
        } 
        if (is.null(rho) && all(c(!is.null(delta), !is.null(phi)))) { 
            rho <- 1 - delta - phi
        } 
        if (is.null(delta) && all(c(!is.null(rho), !is.null(phi)))) { 
            delta <- 1 - rho - phi
        } 
        if (is.null(delta)) { 
            delta <- 0
        } 
        if (is.null(phi)) { 
            phi <- 0
        } 
        if (is.null(rho)) { 
            rho  <- 0
        }
        if(delta + rho + phi != 1) {
            stop(paste("Not enough components provided to set proportions of",
                       "noise variance correctly; if noise variance is only", 
                       "explained by one component, its proportion of variance",
                       "needs to be set to 1; otherwise, the proportion of",
                       "variance of at least 2 components need to be specified")
            )
        }
        if (phi == 1) {
            modelNoise="noiseBgOnly"
            vmessage(c("The noise model is:", modelNoise), verbose=v)
            vmessage(c("Proportion of random noise variance (phi):", phi), 
                     verbose=v)
            vmessage(c("Variance of shared random noise effect (alpha):", 
                       alpha), verbose=v)
        } else if (delta == 1) {
            modelNoise="noiseFixedOnly"
            vmessage(c("The noise model is:", modelNoise), verbose=v)
            vmessage(c("Proportion of fixed noise variance (delta):", 
                       delta), verbose=v)
            vmessage(c("Proportion of variance of shared fixed noise",
                           "effects (gamma):", gamma), verbose=v)
            vmessage(c("Proportion of noise effects (confounders) to have",
                           "a trait-independent effect (pIndependentConfounders"
                            , "):", pIndependentConfounders), 
                     verbose=v)
            vmessage(c("Proportion of traits influenced by independent",
                           "fixed noise effects (pTraitIndependentConfounders):"
                           , pTraitIndependentConfounders), 
                     verbose=v)
        } else if (rho == 1) {
            modelNoise="noiseCorrelatedOnly"
            vmessage(c("The noise model is:", modelNoise), verbose=v)
            vmessage(c("Proportion of variance of correlated noise effects",
                           "(rho):", rho), verbose=v)
            vmessage(c("Correlation between phenotypes (pcorr):", pcorr), 
                     verbose=v)
        } else if (1 - rho == 1) {
            modelNoise="noiseFixedAndBg"
            vmessage(c("The noise model is:", modelNoise), verbose=v)
            vmessage(c("Proportion of fixed noise variance (delta):", 
                       delta), verbose=v)
            vmessage(c("Proportion of variance of shared fixed noise",
                           "effects (gamma):", gamma), verbose=v)
            vmessage(c("Proportion of fixed noise effects to have",
                           "a trait-independent effect (pIndependentConfounders"
                           , "):", pIndependentConfounders), 
                     verbose=v)
            vmessage(c("Proportion of traits influenced by independent",
                           "fixed noise effects (pTraitIndependentConfounders):"
                           , pTraitIndependentConfounders), 
                     verbose=v)            
            vmessage(c("Proportion of random noise variance (phi):", phi), 
                     verbose=v)
            vmessage(c("Variance of shared random noise effect (alpha):", 
                       alpha), verbose=v)
        } else if (1 - phi == 1) {
            modelNoise="noiseFixedAndCorrelated"
            vmessage(c("The noise model is:", modelNoise), verbose=v)
            vmessage(c("Proportion of fixed noise variance (delta):", 
                       delta), verbose=v)
            vmessage(c("Proportion of variance of shared fixed noise",
                           "effects (gamma):", gamma), verbose=v)
            vmessage(c("Proportion of fixed  noise effects to have",
                           "a trait-independent effect (pIndependentConfounders"
                           , "):", pIndependentConfounders), 
                     verbose=v)
            vmessage(c("Proportion of traits influenced by independent",
                           "fixed noise effects (pTraitIndependentConfounders):"
                           , pTraitIndependentConfounders), 
                     verbose=v)
            vmessage(c("Proportion of variance of correlated noise effects",
                           "(rho):", rho), verbose=v)
            vmessage(c("Correlation between phenotypes (pcorr):", pcorr), 
                     verbose=v)
        } else if (1 - delta == 1 ) {
            modelNoise="noiseCorrelatedAndBg"
            vmessage(c("The noise model is:", modelNoise), verbose=v)
            vmessage(paste("Proportion of variance of correlated noise effects",
                           "(rho):", rho), verbose=v)
            vmessage(c("Correlation between phenotypes (pcorr):", pcorr), 
                     verbose=v)
            vmessage(c("Proportion of random noise variance (phi):", phi), 
                     verbose=v)
            vmessage(c("Variance of shared random noise effect (alpha):", 
                       alpha), verbose=v)
        } else {
            modelNoise="noiseFixedAndBgAndCorrelated"
            vmessage(c("The noise model is:", modelNoise), verbose=v)
            vmessage(c("Proportion of fixed noise variance (delta):", 
                       delta), verbose=v)
            vmessage(c("Proportion of variance of shared fixed noise",
                           "effects (gamma):", gamma), verbose=v)
            vmessage(c("Proportion of fixed noise effects to have",
                           "a trait-independent effect (pIndependentConfounders"
                           , "):", pIndependentConfounders), 
                     verbose=v)
            vmessage(c("Proportion of traits influenced by independent",
                           "fixed noise effects (pTraitIndependentConfounders):"
                           , pTraitIndependentConfounders), 
                     verbose=v)
            vmessage(c("Proportion of variance of correlated noise effects",
                           "(rho):", rho), verbose=v)
            vmessage(c("Correlation between phenotypes (pcorr):", pcorr), 
                     verbose=v)
            vmessage(c("Proportion of random noise variance (phi):", phi), 
                     verbose=v)
            vmessage(c("Variance of shared random noise effect (alpha):", 
                       alpha), verbose=v)
        }
        vmessage("\n", verbose=v)
    }

    vmessage(c("The total genetic variance (genVar) is:", genVar), verbose=v)
    if ( genVar == 0 ) {
        modelGenetic="noGenetic"
        vmessage(c("The genetic model is:", modelGenetic), verbose=v)
    } else {   
        if (all(c(is.null(h2bg), is.null(h2s)))) {
            stop(paste("Neither fixed genetic effect variance (h2s) nor random"
                    , "genetic effect variance (h2bg) provided, at least one is"
                    , "required")
            )
        }        
        if (length(c(h2bg, h2s)) == 2 && sum(c(h2bg, h2s)) != 1) {
            stop(paste("Sum of the proportion of the variance of genetic",
                       "effects is not equal to 1; change h2s (fixed effect",
                       "variance) or h2bg (random effect variance) such that",
                       "h2s + h2bg = 1")
            )
        } 
        if (is.null(h2s)) {
            h2s <- 1 - h2bg
        } else {
            h2bg <- 1 - h2s
        }

        if ( h2s == 0) {
            modelGenetic="geneticBgOnly"
            vmessage(c("The genetic model is:", modelGenetic), verbose=v)
            vmessage(c("Proportion of variance of random genetic effects",
                           "(h2bg):", h2bg), verbose=v)
            vmessage(c("Proportion of variance of shared random genetic",
                           "effects (eta):", eta), verbose=v)
        } else if ( h2s == 1) {
            modelGenetic="geneticFixedOnly"
            vmessage(c("The genetic model is:", modelGenetic), verbose=v)
            vmessage(c("Proportion of variance of fixed genetic effects",
                           "(h2s):", h2s), verbose=v)
            vmessage(c("Proportion of variance of shared fixed genetic",
                           "effects (theta):", theta), verbose=v)
            vmessage(c("Proportion of fixed genetic effects to have a", 
                           "trait-independent fixed effect", 
                           "(pIndependentGenetic):", pIndependentGenetic), 
                     verbose=v)
            vmessage(c("Proportion of traits influenced by independent",
                           "fixed genetic effects (pTraitIndependentGenetic):",
                           pTraitIndependentGenetic), 
                     verbose=v)
        } else {
            modelGenetic="geneticFixedAndBg"
            vmessage(c("The genetic model is:", modelGenetic), verbose=v)
            vmessage(c("Proportion of variance of fixed genetic effects",
                           "(h2s):", h2s), verbose=v)
            vmessage(c("Proportion of variance of shared fixed genetic",
                           "effects (theta):", theta), verbose=v)
            vmessage(c("Proportion of fixed genetic effects to have a", 
                           "trait-independent fixed effect", 
                           "(pIndependentGenetic):", pIndependentGenetic), 
                     verbose=v)
            vmessage(c("Proportion of traits influenced by independent",
                           "fixed genetic effects (pTraitIndependentGenetic):",
                           pTraitIndependentGenetic), 
                     verbose=v)
            vmessage(c("Proportion of variance of random genetic effects",
                           "(h2bg):", h2bg), verbose=v)
            vmessage(c("Proportion of variance of shared random genetic",
                           "effects (eta):", eta), verbose=v)
        }
        vmessage("\n", verbose=v)
    }
    return(list(modelGenetic=modelGenetic, modelNoise=modelNoise,
                genVar=genVar, h2s=h2s, h2bg=h2bg, noiseVar=noiseVar, rho=rho, 
                delta=delta, phi=phi, gamma=gamma, theta=theta, eta=eta, 
                alpha=alpha, pcorr=pcorr))
}


#' Run phenotype simulation.
#'
#' runSimulation wraps around the phenotype component functions (genFixedEffects
#' , genBgEffects, noiseBgEffects, noiseFixedEffects and correlatedBgEffects),
#' rescales each component and combines them into the final phenotype.
#'
#' @param P number [integer] of phenotypes to simulate 
#' @param N number [integer] of samples to simulate
#' @param tNrSNP total number [integer] of SNPs to simulate; these SNPs are used
#'  for kinship estimation
#' @param cNrSNP number [integer] of causal SNPs; used as genetic fixed effects
#' @param NrConfounders number [integer] of confounders; used as noise fixed 
#' effects
#' @param NrFixedEffects number [integer] of different fixed effects to simulate
#' ; allows to simulate fixed effects from different distributions or with 
#'  different parameters
#' @param chr numeric vector of chromosomes to chose NrCausalSNPs from; only 
#' used when external genotype data is provided i.e. is.null(genoFilePrefix) == 
#' FALSE
#' @param NrChrCausal Number [integer] of causal chromosomes to  chose 
#' NrCausalSNPs from (as opposed to the actual chromosomes to chose from via chr
#' );  only used when external genotype data is provided i.e. 
#' is.null(genoFilePrefix) == FALSE. 
#' @param SNPfrequencies vector of allele frequencies [double] from which to 
#' sample
#' @param genotypefile needed when loading entire genotype set into memory, 
#' path/to/genotype file [string] in format specified by \link{format}
#' @param format needed when loading entire genotype set into memory, specifies 
#' the format of the genotype data; has to be one of  plink, oxgen, genome, 
#' bimbam and delim; for details see Details in \link{readStandardGenotypes}
#' @param genoFilePrefix full path/to/chromosome-wise-genotype-file-ending-
#' before-"chrChromosomeNumber" (no '~' expansion!) [string]
#' @param genoFileSuffix [string] following chromosome number including 
#' .fileformat (e.g. ".csv"); has to be a text format i.e. comma/tab/space
#' separated
#' @param genoDelimiter field separator [string] of genotypefile or genoFile 
#' @param kinshipfile path/to/kinshipfile [string] to be read; either X or 
#' kinshipfile must be provided
#' @param kinshipDelimiter field separator [string] of kinship file 
#' @param kinshipHeader [boolean], if TRUE kinship file has header information 
#' @param standardise [boolean], if TRUE standardised genotypes will be returned
#' @param distConfounders name [string] of distribution to use to simulate 
#' confounders; one of "unif", "norm", "bin", "cat_norm", "cat_unif"
#' @param mConfounders mean/midpoint [double] of normal/uniform distribution for
#'  confounders
#' @param sdConfounders standard deviation/extension from midpoint [double] of 
#' normal/uniform distribution for confounders
#' @param catConfounders confounder categories [factor]; required if 
#' distConfounders "cat_norm" or "cat_unif" 
#' @param probConfounders probability [double] of binomial confounders (0/1); 
#' required if distConfounders "bin" 
#' @param distBetaConfounders name [string] of distribution to use to simulate 
#' effect sizes of confounders; one of "unif" or "norm"
#' @param mBetaConfounders mean/midpoint [double] of normal/uniform distribution 
#' for effect sizes of confounders
#' @param sdBetaConfounders standard deviation/extension from midpoint [double] 
#' of normal/uniform distribution for effect sizes of confounders
#' @param pIndependentConfounders Proportion [double] of noise effects 
#' (confounders) to have a trait-independent effect
#' @param pTraitIndependentConfounders Proportion [double] of traits influenced 
#' by independent fixed noise effects
#' @param distBetaGenetic name [string] of distribution to use to simulate 
#' effect sizes of SNPs; one of "unif" or "norm"
#' @param mBetaGenetic mean/midpoint [double] of normal/uniform distribution 
#' for effect sizes of SNPs
#' @param sdBetaGenetic standard deviation/extension from midpoint [double] 
#' of normal/uniform distribution for effect sizes of SNPs
#' @param pIndependentGenetic Proportion [double] of genetic effects (SNPs) to 
#' have a trait-independent fixed effect
#' @param pTraitIndependentGenetic Proportion [double] of traits influenced by 
#' independent fixed genetic effects
#' @param meanNoiseBg mean [double] of the normal distribution for noise bg 
#' effects
#' @param sdNoiseBg standard deviation [double] of the normal distribution for 
#' noise bg effects
#' @param sampleID prefix [string] for naming samples (followed by sample number
#'  from 1 to N)
#' @param phenoID prefix [string] for naming traits (followed by trait number 
#' from 1 to P)
#' @param snpID prefix [string] for naming SNPs (followed by SNP number from 1 
#' to NrSNP)
#' @param genVar Proportion [double] of total genetic variance
#' @param h2s Proportion [double] of gentic variance of fixed effects 
#' @param h2bg Proportion [double] of genetic variance of background effects; 
#' either h2s or h2bg have to be specified and h2s + h2bg = 1
#' @param theta Proportion [double] of variance of shared fixed genetic effects
#' @param eta Proportion [double] of variance of shared bg genetic effects
#' @param noiseVar Proportion [double] of total noise variance
#' @param rho Proportion [double] of noise variance of correlated effects; sum 
#' of rho, delta and phi has to be equal 1
#' @param pcorr Correlation [double] between phenotypes
#' @param delta Proportion [double] of noise variance of fixed effects; sum of 
#' rho, delta and phi  has to be equal 1
#' @param gamma Proportion [double] of variance of shared fixed noise effects
#' @param phi Proportion [double] of noise variance of background effects; sum 
#' of rho, delta and phi  has to be equal 1
#' @param alpha Variance [double] of shared bg noise effect
#' @param seed seed [integer] to initiate random number generation
#' @param verbose [boolean]; if TRUE, progress info is printed to standard out
#' @return named list of absolute levels of variance explained by each phenotype 
#' component (varComponents), 
#' the scaled phenotype components (phenoComponents) and and parameters used for 
#' phenotype set-up and labeling (setup)
#' @seealso \code{\link{createPheno}}, which this function wraps
#' @export
#' @examples
#' # simulate phenotype of 100 samples, 10 traits from genetic and noise 
#' # background effects, with variance explained of 0.2 and 0.8 respectively
#' genVar = 0.2
#' simulatedPhenotype <- runSimulation(N=100, P=5, cNrSNP=10,
#' genVar=genVar, h2s=1, phi=1)
runSimulation <- function(N=1000, P=10, tNrSNP=5000, cNrSNP=20, 
                          NrConfounders=10, seed=219453, 
                          chr=NULL, NrChrCausal=NULL,
                          genVar=NULL, h2s=NULL, theta=0.8, h2bg=NULL, eta=0.8, 
                          noiseVar=NULL, rho=NULL, delta=NULL, gamma=0.8, 
                          phi=NULL, alpha=0.8, 
                          sampleID="ID_", phenoID="Trait_", snpID="SNP_",
                          genotypefile=NULL, format=NULL,
                          genoFilePrefix=NULL, genoFileSuffix=NULL, 
                          genoDelimiter=",", kinshipfile=NULL, 
                          kinshipHeader=TRUE, kinshipDelimiter=",", 
                          standardise=TRUE,
                          NrFixedEffects=1, distConfounders="norm",
                          mConfounders=0, sdConfounders=1,
                          catConfounders=NULL, probConfounders=NULL,
                          distBetaConfounders="norm", mBetaConfounders=0, 
                          sdBetaConfounders=1,
                          distBetaGenetic="norm", mBetaGenetic=0, 
                          sdBetaGenetic=1,
                          pIndependentConfounders=0.4, 
                          pTraitIndependentConfounders=0.2, 
                          pcorr=0.8, meanNoiseBg=0, sdNoiseBg=1, 
                          SNPfrequencies=c(0.1, 0.2, 0.4), 
                          pIndependentGenetic=0.4, pTraitIndependentGenetic=0.2,
                          verbose=TRUE) {

    vmessage(c("Set seed:", seed), verbose=verbose)
    set.seed(seed)
       # find model
    model <- setModel(genVar=genVar, h2s=h2s, h2bg=h2bg, theta=theta, eta=eta, 
                      noiseVar=noiseVar, rho=rho, delta=delta, gamma=gamma, 
                      phi=phi, alpha=alpha, pcorr=pcorr, 
                      pIndependentConfounders=pIndependentConfounders,  
                      pTraitIndependentConfounders=pTraitIndependentConfounders, 
                      pIndependentGenetic=pIndependentGenetic, 
                      pTraitIndependentGenetic=pTraitIndependentGenetic, 
                      v=verbose)
    
    id_traits <- paste(phenoID, seq(1, P, 1), sep="")
    id_samples <- paste(sampleID, seq(1, N, 1), sep="")

    ### create simulated phenotypes
    # 1. Simulate genetic terms
    vmessage(c("Simulate genetic effects (genetic model:", 
               model$modelGenetic,")"), verbose=verbose)
    if (grepl('Fixed', model$modelGenetic)) {
        if (is.null(genoFilePrefix) && is.null(genotypefile)) {
            if (!grepl('Bg', model$modelGenetic) && tNrSNP != cNrSNP) {
                warning(paste("The genetic model does not contain random",
                            "effects but the total number of SNPs to",
                            "simulate (tNrSNP:",
                            tNrSNP, ") is larger than the causal number of",
                            "SNPs (cNrSNP:", cNrSNP, "). If genotypes are",
                            "not needed, consider setting tNrSNP=cNrSNP",
                            "to speed up computation"))
            }
            genotypes <- simulateGenotypes(N=N, NrSNP=tNrSNP, 
                                           frequencies=SNPfrequencies, 
                                           sampleID=sampleID, 
                                           snpID=snpID, 
                                           verbose=verbose)
        } else if (! is.null(genotypefile)) {
            genotypes <- readStandardGenotypes(genotypefile, format=format,
                                               verbose=verbose, 
                                               sampleID = sampleID, 
                                               snpID = snpID, 
                                               delimiter = genoDelimiter)
            id_samples <- genotypes$id_samples
        } else {
            genotypes <- NULL
        }
        causalSNPs <- getCausalSNPs(NrCausalSNPs=cNrSNP, chr=chr, 
                                    NrChrCausal=NrChrCausal,
                                    genotypes=genotypes$genotypes,
                                    genoFilePrefix=genoFilePrefix, 
                                    genoFileSuffix=genoFileSuffix, 
                                    delimiter=genoDelimiter, 
                                    sampleID=sampleID, 
                                    verbose=verbose)
        id_snps <- genotypes$id_snps
        
        vmessage("Simulate genetic fixed effects", verbose=verbose)
        genFixed <- geneticFixedEffects(X_causal=causalSNPs, N=N, P=P, 
                                        pIndependentGenetic=
                                            pIndependentGenetic, 
                                        pTraitIndependentGenetic=
                                            pTraitIndependentGenetic,
                                        distBeta=distBetaGenetic, 
                                        mBeta=mBetaGenetic, 
                                        sdBeta=sdBetaGenetic)
        
        var_genFixed_shared <- model$theta * model$h2s * model$genVar
        var_genFixed_independent <- (1 - model$theta) * model$h2s * model$genVar
        
        genFixed_shared_rescaled <- rescaleVariance(genFixed$shared, 
                                                    var_genFixed_shared)
        genFixed_independent_rescaled <- rescaleVariance(genFixed$independent, 
                                                var_genFixed_independent)
        
        Y_genFixed <- addNonNulls(list(genFixed_shared_rescaled, 
                                       genFixed_independent_rescaled))
        colnames(Y_genFixed) <- id_traits
        rownames(Y_genFixed) <- id_samples
        
    } else {
        genFixed <- NULL
        genotypes <- NULL
        var_genFixed_shared <- 0
        var_genFixed_independent <- 0
        Y_genFixed <- NULL
        cNrSNP <- 0
    }
    if (grepl('Bg', model$modelGenetic)) {
        if (is.null(kinshipfile)) {
            if (is.null(genotypes) && !is.null(genotypefile)){
                genotypes <- readStandardGenotypes(genotypefile, format=format,
                                                   verbose=verbose, 
                                                   sampleID = sampleID, 
                                                   snpID = snpID, 
                                                   delimiter = delimiter)
                id_samples <- genotypes$id_samples
                id_snps <- genotypes$id_snps
            }
            if (is.null(genotypes)) {
                genotypes <- simulateGenotypes(N=N, NrSNP=tNrSNP, 
                                               frequencies=SNPfrequencies, 
                                               sampleID=sampleID, 
                                               verbose=verbose)
            }
            kinship <- getKinship(X=genotypes$genotypes, standardise=standardise, 
                                  sampleID=sampleID, 
                                  verbose=verbose)
        } else {
            vmessage("Read kinship from file", verbose=verbose)
            kinship <- getKinship(kinshipfile=kinshipfile, sampleID=sampleID, 
                                  sep=kinshipDelimiter, 
                                  header=kinshipHeader, verbose=verbose)
        }
        
        vmessage("Simulate genetic background effects", verbose=verbose)
        genBg <- geneticBgEffects(P=P, kinship=kinship)
        
        var_genBg_shared <- model$eta * model$h2bg * model$genVar
        var_genBg_independent <- (1 - model$eta) * model$h2bg * model$genVar
        
        genBg_shared_rescaled <- rescaleVariance(genBg$shared, var_genBg_shared)
        genBg_independent_rescaled <- rescaleVariance(genBg$independent,
                                                      var_genBg_independent)
        
        Y_genBg <- addNonNulls(list(genBg_shared_rescaled, 
                                    genBg_independent_rescaled))
        colnames(Y_genBg) <- id_traits
        rownames(Y_genBg) <- id_samples
        
        cov_Y_genBg <- t(Y_genBg) %*% Y_genBg
        cov_Y_genBg <-  cov_Y_genBg/mean(diag(cov_Y_genBg))
        diag(cov_Y_genBg) <- diag(cov_Y_genBg) + 1e-4
    } else {
        genBg <- NULL
        kinship <- NULL
        var_genBg_shared <- 0
        var_genBg_independent <- 0
        Y_genBg <- NULL
        cov_Y_genBg <- NULL
    }
    # 1. Simulate noise terms
    vmessage(c("Simulate noise terms (noise model:", model$modelNoise, ")"),
             verbose=verbose)
    if (grepl('Correlated', model$modelNoise))  {
        vmessage("Simulate correlated background effects", verbose=verbose)
        correlatedBg <- correlatedBgEffects(N=N, P=P, pcorr=model$pcorr)
        var_noiseCorrelated <- model$rho *  model$noiseVar
        correlatedBg_rescaled <- rescaleVariance(correlatedBg, 
                                                 var_noiseCorrelated)
        
        Y_correlatedBg <- correlatedBg_rescaled
        colnames(Y_correlatedBg) <- id_traits
        rownames(Y_correlatedBg) <- id_samples
    } else {
        correlatedBg <- NULL
        var_noiseCorrelated <- 0
        Y_correlatedBg <- NULL 
    }
    if (grepl('Bg', model$modelNoise)) {
        vmessage("Simulate noise background effects", verbose=verbose)
        noiseBg <- noiseBgEffects(N=N, P=P, mean=meanNoiseBg, sd=sdNoiseBg)
        
        var_noiseBg_shared <- model$alpha * model$phi * model$noiseVar
        var_noiseBg_independent <- (1 - model$alpha) * model$phi * model$noiseVar
        
        noiseBg_shared_rescaled <- rescaleVariance(noiseBg$shared, 
                                                   var_noiseBg_shared)
        noiseBg_independent_rescaled <- rescaleVariance(noiseBg$independent, 
                                                        var_noiseBg_independent)
        
        Y_noiseBg <- addNonNulls(list(noiseBg_shared_rescaled, 
                                      noiseBg_independent_rescaled))
        colnames(Y_noiseBg) <- id_traits
        rownames(Y_noiseBg) <- id_samples
        cov_Y_noiseBg <- t(Y_noiseBg) %*% Y_noiseBg
        cov_Y_noiseBg <-  cov_Y_noiseBg/mean(diag(cov_Y_noiseBg))
        diag(cov_Y_noiseBg) <- diag(cov_Y_noiseBg) + 1e-4
    } else {
        noiseBg <- NULL
        var_noiseBg_shared <- 0
        var_noiseBg_independent <- 0
        Y_noiseBg <- NULL
        cov_Y_noiseBg <- NULL
    }
    if (grepl('Fixed', model$modelNoise)) {
        vmessage("Simulate noise fixed effects", verbose=verbose)
        noiseFixed <- noiseFixedEffects(P=P, N=N, 
                                        NrFixedEffects = NrFixedEffects,
                                        NrConfounders=NrConfounders,
                                        pIndependentConfounders=
                                            pIndependentConfounders,
                                        pTraitIndependentConfounders=
                                            pTraitIndependentConfounders,
                                        distConfounders=distConfounders, 
                                        mConfounders=mConfounders, 
                                        sdConfounders=sdConfounders, 
                                        catConfounders=catConfounders, 
                                        probConfounders = probConfounders,
                                        distBeta=distBetaConfounders, 
                                        mBeta=mBetaConfounders, 
                                        sdBeta=sdBetaConfounders)
        
        
        var_noiseFixed_shared <- model$gamma * model$delta * model$noiseVar
        var_noiseFixed_independent <- (1 - model$gamma) * model$delta * 
            model$noiseVar
        
        noiseFixed_shared_rescaled <- rescaleVariance(noiseFixed$shared, 
                                                      var_noiseFixed_shared)
        noiseFixed_independent_rescaled <- rescaleVariance(
            noiseFixed$independent, 
            var_noiseFixed_independent)
        
        Y_noiseFixed <- addNonNulls(list(noiseFixed_shared_rescaled, 
                                         noiseFixed_independent_rescaled))
        colnames(Y_noiseFixed) <- id_traits
        rownames(Y_noiseFixed) <- id_samples
    } else {
        noiseFixed <- NULL
        var_noiseFixed_shared <- 0
        var_noiseFixed_independent <- 0
        Y_noiseFixed <- NULL
    }

    
    # 3. Construct final simulated phenotype 
    vmessage("Construct final simulated phenotype"
             , verbose=verbose)
    
    components <- list(Y_genFixed, Y_genBg, Y_noiseFixed, Y_noiseBg, 
                       Y_correlatedBg)
    Y <-  addNonNulls(components)
    Y <- scale(Y)
    colnames(Y) <- id_traits
    rownames(Y) <- id_samples
    
    varComponents <- data.frame(genVar=model$genVar, h2s=model$h2s, 
                                h2bg=model$h2bg, 
                                var_genFixed_shared=var_genFixed_shared, 
                                var_genFixed_independent=
                                    var_genFixed_independent, 
                                var_genBg_shared=var_genBg_shared, 
                                var_genBg_independent=var_genBg_independent, 
                                noiseVar=model$noiseVar, 
                                var_noiseFixed_shared=var_noiseFixed_shared, 
                                var_noiseFixed_independent=
                                    var_noiseFixed_independent, 
                                var_noiseBg_shared=var_noiseBg_shared, 
                                var_noiseBg_independent=var_noiseBg_independent, 
                                var_noiseCorrelated=var_noiseCorrelated)
    phenoComponents <- list(Y=Y, Y_genFixed=Y_genFixed, Y_genBg=Y_genBg, 
                            Y_noiseFixed=Y_noiseFixed, Y_noiseBg=Y_noiseBg, 
                            Y_correlatedBG=Y_correlatedBg, 
                            cov_Y_genBg=cov_Y_genBg, 
                            cov_Y_noiseBg=cov_Y_noiseBg, genFixed=genFixed, 
                            noiseFixed=noiseFixed)
    setup <- list(P=P, N=N, NrCausalSNPs=cNrSNP, 
                  modelGenetic=model$modelGenetic, 
                  modelNoise=model$modelNoise, 
                  id_samples=id_samples, id_traits=id_traits, id_snps=id_snps)
    rawComponents <- list(kinship=kinship, genotypes=genotypes)

    return(list(varComponents=varComponents, phenoComponents=phenoComponents, 
                setup=setup, rawComponents=rawComponents))
}



