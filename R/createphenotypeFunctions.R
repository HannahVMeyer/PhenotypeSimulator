#' Set simulation model.
#'
#' Based on parameters provided, this function sets the name for the phenotype 
#' simulation. The model name is needed downstream for phenotype component 
#' simulations, but can also be set manually.
#'
#' @param genVar Total genetic variance [double] 
#' @param h2s Proportion [double] of variance of fixed genetic effects
#' @param h2bg Proportion [double] of variance of random genetic effects
#' @param theta Proportion [double] of variance of common fixed genetic effects
#' @param eta Proportion [double] of variance of common bg genetic effects
#' @param noiseVar Total genetic variance [double] 
#' @param rho Proportion [double] of variance of correlated noise effects
#' @param pcorr Correlation [double] between phenotypes
#' @param delta Proportion [double] of fixed noise variance
#' @param gamma Proportion [double] of variance of common fixed noise effects
#' @param phi Proportion [double] of variance of background noise effects
#' @param alpha Proportion [double] of Variance of common bg noise effect
#' @param pSpecificConfounders Proportion [double] of noise effects 
#' (confounders) to have a trait-specific effect
#' @param pTraitSpecificConfounders Proportion [double] of traits influenced by 
#' specific fixed noise effects
#' @param pSpecificGenetic Proportion [double] of genetic effects (SNPs) to have
#'  a trait-specific fixed effect
#' @param pTraitSpecificGenetic Proportion [double] of traits influenced by 
#' specific fixed genetic effects
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
                     alpha=0.8, pcorr=0.6, pSpecificConfounders=0.4,  
                     pTraitSpecificConfounders=0.2,  pSpecificGenetic=0.4, 
                     pTraitSpecificGenetic=0.2, v=TRUE)  {
    if (is.null(c(genVar, noiseVar, h2bg, h2s, delta, rho, phi))) {
        stop("No variance components specified")
    }
    if (is.null(genVar)) {
        if (is.null(noiseVar)) {
            stop("Neither genVar nor noiseVar are provided, thus proportion of 
                 variance from genetics cannot be deduced")
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
    if ((noiseVar) == 0 ) {
        modelNoise="noNoise"
        vmessage(c("The noise model is:", modelNoise), verbose=v)
    } else {
        if (all(c(is.null(delta), is.null(rho), is.null(phi)))) {
            stop("Neither delta nor rho or phi are provided, at least one is 
                 required")
        } else if (all(c(!is.null(delta), !is.null(rho), !is.null(phi)))) {
            if(delta + rho + phi != 1) {
                stop("Sum of the proportion of the variance of noise effects is 
                     greater than 1; change noiseVar, delta (fixed effect 
                     variance), rho (correlated background variance) or phi 
                     (random effect variance)\n delta + rho + phi = 1")
            }
            } else if (all(c(!is.null(delta), !is.null(rho)))) { 
                phi <- 1 - delta - rho
            } else if (all(c(!is.null(delta), !is.null(phi)))) { 
                rho <- 1 - delta - phi
            } else if (all(c(!is.null(rho), !is.null(phi)))) { 
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
        
        if (phi == 1) {
            modelNoise="noiseBgOnly"
            vmessage(c("The noise model is:", modelNoise), verbose=v)
            vmessage(c("Variance of common bg noise effect:", alpha), verbose=v)
        } else if (delta == 1) {
            modelNoise="noiseFixedOnly"
            vmessage(c("The noise model is:", modelNoise), verbose=v)
            vmessage(c("Proportion of fixed noise variance:", delta), verbose=v)
            vmessage(c("Proportion of variance of common fixed noise effects:",
                       gamma), verbose=v)
            vmessage(c("Proportion of noise effects (confounders) to have a 
                     trait-specific effect :", pSpecificConfounders), verbose=v)
            vmessage(c("Proportion of traits influenced by specific fixed noise 
                       effects:", pTraitSpecificConfounders), verbose=v)
        } else if (rho == 1) {
            modelNoise="noiseCorrelatedOnly"
            vmessage(c("The noise model is:", modelNoise), verbose=v)
            vmessage(c("Correlation between phenotypes:", pcorr), verbose=v)
            vmessage(c("Proportion of variance of correlated bg noise effects:", 
                       rho), verbose=v)
        } else if (1 - rho == 1) {
            modelNoise="noiseFixedAndBg"
            vmessage(c("The noise model is:", modelNoise), verbose=v)
            vmessage(c("Proportion of fixed noise variance:", delta), verbose=v)
            vmessage(c("Proportion of variance of common fixed noise effects:", 
                       gamma), verbose=v)
            vmessage(c("Proportion of noise effects (confounders) to have a 
                       trait-specific effect :", pSpecificConfounders), 
                     verbose=v)
            vmessage(c("Proportion of traits influenced by specific fixed noise 
                       effects:", pTraitSpecificConfounders), verbose=v)
            vmessage(c("Variance of common bg noise effect:", alpha), verbose=v)
        } else if (1 - phi == 1) {
            modelNoise="noiseFixedAndCorrelated"
            vmessage(c("The noise model is:", modelNoise), verbose=v)
            vmessage(c("Proportion of fixed noise variance:", delta), verbose=v)
            vmessage(c("Proportion of variance of common fixed noise effects:", 
                       gamma), verbose=v)
            vmessage(c("Proportion of noise effects (confounders) to have a 
                       trait-specific effect :", pSpecificConfounders), 
                     verbose=v)
            vmessage(c("Proportion of traits influenced by specific fixed noise 
                       effects:", pTraitSpecificConfounders), verbose=v)
            vmessage(c("Correlation between phenotypes:", pcorr), verbose=v)
            vmessage(c("Proportion of variance of correlated bg noise effects:", 
                       rho), verbose=v)
        } else if (1 - delta == 1 ) {
            modelNoise="noiseCorrelatedAndBg"
            vmessage(c("The noise model is:", modelNoise), verbose=v)
            vmessage(c("Correlation between phenotypes:", pcorr), verbose=v)
            vmessage(c("Proportion of variance of correlated bg noise effects:",
                       rho), verbose=v)
            vmessage(c("Variance of common bg noise effect:", alpha), verbose=v)
        } else {
            modelNoise="noiseFixedAndBgAndCorrelated"
            vmessage(c("The noise model is:", modelNoise), verbose=v)
            vmessage(c("Proportion of fixed noise variance:", delta), verbose=v)
            vmessage(c("Proportion of variance of common fixed noise effects:",
                       gamma), verbose=v)
            vmessage(c("Proportion of noise effects (confounders) to have a 
                       trait-specific effect :", pSpecificConfounders), 
                     verbose=v)
            vmessage(c("Proportion of traits influenced by specific fixed noise 
                       effects:", pTraitSpecificConfounders), verbose=v)
            vmessage(c("Correlation between phenotypes:", pcorr), verbose=v)
            vmessage(c("Proportion of variance of correlated bg noise effects:", 
                       rho), verbose=v)
            vmessage(c("Variance of common bg noise effect:", alpha), verbose=v)
        }
    }

    if ( genVar == 0 ) {
        modelGenetic="noGenetic"
        vmessage(c("The genetic model is:", modelGenetic), verbose=v)
    } else {   
        if (all(c(is.null(h2bg), is.null(h2s)))) {
            stop("Neither h2bg nor h2s provided, at least one is required")
        }
        if (is.null(h2s)) {
            h2s <- 1 - h2bg
        } else {
            h2bg <- 1 - h2s
        }
        if (h2bg + h2s != 1) {
            stop("Sum of the proportion of the variance of genetic effects is 
                 greater than 1; change h2s (fixed effect variance) or h2bg 
                 (random effect variance)\n h2s + h2bg = 1")
        } 
        if ( h2s == 0) {
            modelGenetic="geneticBgOnly"
            vmessage(c("The genetic model is:", modelGenetic), verbose=v)
            vmessage(c("Proportion of variance of common bg genetic effects:", 
                       eta), verbose=v)
        } else if ( h2s == 1) {
            modelGenetic="geneticFixedOnly"
            vmessage(c("The genetic model is:", modelGenetic), verbose=v)
            vmessage(c("Proportion of variance of fixed genetic effects:", h2s), 
                     verbose=v)
            vmessage(c("Proportion of variance of common fixed genetic effects:"
                       , theta), verbose=v)
            vmessage(c("Proportion of genetic effects (SNPs) to have a 
                       trait-specific fixed effect:", pSpecificGenetic), 
                     verbose=v)
            vmessage(c("Proportion of traits influenced by specific fixed 
                        genetic effects:", pTraitSpecificGenetic), verbose=v)
        } else {
            modelGenetic="geneticFixedAndBg"
            vmessage(c("The genetic model is:", modelGenetic), verbose=v)
            vmessage(c("Proportion of variance of fixed genetic effects:", h2s), 
                     verbose=v)
            vmessage(c("Proportion of variance of common fixed genetic effects:"
                       , theta), verbose=v)
            vmessage(c("Proportion of genetic effects (SNPs) to have a 
                       trait-specific fixed effect:", pSpecificGenetic), 
                     verbose=v)
            vmessage(c("Proportion of traits influenced by specific fixed 
                        genetic effects:", pTraitSpecificGenetic), verbose=v)
        }
    }
    return(list(modelGenetic=modelGenetic, modelNoise=modelNoise))
}

### phenotype function

#' Combine different phenotype components.
#'
#' createPheno takes precomputed phenotype components and rescales them 
#' according to the specified proportion of variance they should explain in the 
#' final phenotype.
#'
#' @param P number [integer] of phenotypes to simulate 
#' @param N number [integer] of samples to simulate
#' @param sampleID prefix [string] for naming samples (followed by sample number
#'  from 1 to N)
#' @param phenoID prefix [string] for naming traits (followed by trait number 
#' from 1 to P)
#' @param genBg list of specific and common genetic background effects as 
#' obtained by \code{\link{geneticBgEffects}}
#' @param genFixed list of specific and common genetic fixed effects as 
#' obtained by \code{\link{geneticFixedEffects}}
#' @param noiseBg list of specific and common noise background effects as 
#' obtained by \code{\link{noiseBgEffects}}
#' @param noiseFixed list of specific and common noise fixed effects as 
#' obtained by \code{\link{noiseFixedEffects}}
#' @param correlatedBg list of correlated background effects as obtained 
#' by \code{\link{correlatedBgEffects}}
#' @param genVar Proportion [double] of total genetic variance
#' @param h2s Proportion [double] of variance of fixed genetic effects
#' @param theta Proportion [double] of variance of common fixed genetic effects
#' @param h2bg Proportion [double] of variance of background genetic effects; 
#' either h2s or h2b have to be specified and 
#' h2s + h2b = genVar
#' @param eta Proportion [double] of variance of common bg genetic effects
#' @param noiseVar Proportion [double] of total noise variance
#' @param rho Proportion [double] of variance of correlated noise effects
#' @param pcorr Correlation [double] between phenotypes
#' @param delta Proportion [double] of fixed noise variance
#' @param gamma Proportion [double] of variance of common fixed noise effects
#' @param phi Proportion [double] of variance of background noise effects
#' @param alpha Proportion [double] of variance of common bg noise effect
#' @param modelNoise name [string] of noise model for the phenotype simulation; 
#' based on model specification, phenotype components will be added to the final 
#' phenotype. Possible models: "noNoise", "noiseFixedOnly", "noiseBgOnly", 
#' "noiseCorrelatedOnly", "noiseFixedAndBg","noiseCorrelatedAndBg", 
#' "noiseFixedAndCorrelated", "noiseFixedAndBgAndCorrelated"
#' @param modelGenetic name [string] of genetic model for the phenotype 
#' simulation; based on model specification, phenotype components will be added 
#' to the final phenotype. Possible models: "noGenetic","geneticBgOnly", 
#' "geneticFixedOnly","geneticFixedAndBg"
#' @param verbose [boolean]; if TRUE, progress info is printed to standard out
#' @return named list of absolute levels of variance explained by each phenotype 
#' component (varComponents), 
#' the scaled phenotype components (phenoComponents) and and parameters used for 
#' phenotype set-up and labeling (setup)
#' @seealso \code{\link{rescaleVariance}} for the function used to rescale 
#' phenotype components
#' @export
#' @examples
#' noiseBg <- noiseBgEffects(N=100, P=50)
#' phenotypeNoiseBgOnly <- createPheno(P=50, N=100, noiseBg=noiseBg, 
#' modelNoise="noiseBgOnly", modelGenetic="noGenetic")
#'
#' genotypes <- simulateGenotypes(N=100, NrSNP=50)
#' causalSNPs <- getCausalSNPs(genotypes=genotypes, standardise=FALSE)
#' genFixed <- geneticFixedEffects(N=100, P=50, X_causal=causalSNPs)
#' phenotypeNoiseBgOnlyGenFixed <- createPheno(P=50, N=100, noiseBg=noiseBg,
#' genFixed=genFixed, modelNoise="noiseBgOnly", 
#' modelGenetic="geneticFixedOnly", genVar=0.05)
createPheno <- function(P, N, sampleID="ID_", phenoID="trait_", 
                        correlatedBg=NULL, genFixed=NULL, genBg=NULL, 
                        noiseFixed=NULL, noiseBg=NULL, genVar=NULL, h2s=NULL, 
                        h2bg=NULL, noiseVar=NULL, rho=NULL, delta=NULL, phi=NULL
                        , gamma=0.8, theta=0.8, eta=0.8, alpha=0.8, pcorr=0.6, 
                        modelNoise="noNoise",  modelGenetic="noGenetic", 
                        verbose=TRUE) {
    if (is.null(c(correlatedBg, genFixed, genBg, noiseFixed, noiseBg))) {
        stop("No phenotype components provided, at least one is required")
    }
    if (!is.null(genFixed) && !grepl("Fixed", modelGenetic)) {
        stop(paste("Genetic fixed effect is provided, but the genetic model is", 
        modelGenetic))
    }
    if (!is.null(genBg) && !grepl("Bg", modelGenetic)) {
        stop(paste("Genetic bg effect is provided, but the genetic model is", 
                   modelGenetic))
    }
    if (!is.null(noiseFixed) && !grepl("Fixed", modelNoise)) {
        stop(paste("Noise fixed effect is provided, but the noise model is", 
                   modelNoise))
    }
    if (!is.null(correlatedBg) && !grepl("Correlated", modelNoise)) {
        stop(paste("Correlated noise effect is provided, but the noise model is"
                   , modelNoise))
    }
    if (!is.null(noiseBg) && !grepl("Bg", modelNoise)) {
        stop(paste("Noise background effect is provided, but the noise model is"
                   , modelNoise))
    }
    ### i) rescale different components based on their proportional contribution
    # 2 main variance components; genVar + noiseVar =1
    if (modelGenetic == "noGenetic") {
        if (!is.null(genVar) && genVar != 0) {
            stop("geneticModel is noGenetic but genVar is unequal to zero")
        } else {
            genVar <- 0
        }
    } else {
        if (is.null(genVar)) {
            if (is.null(noiseVar)) {
                stop(paste("geneticModel is ", modelGenetic, "but neither genVar
                nor noiseVar are provided, thus proportion of variance from
                genetics cannot be deduced"))
            } else {
                genVar <- 1 - noiseVar
            } 
        } 
        if (is.null(c(h2s, h2bg))){
            if(!is.null(genFixed)) {
                h2s <- 1
                h2bg <- 0
            } else if (!is.null(genBg)) {
                h2bg <- 1
                h2s <- 0
            } else {
                stop("Genetic variance is unequal to zero, but no genetic 
                phenotype components are provided")
            }
        } else if (is.null(h2s)) {
            h2s <- 1 - h2bg
        } else {
            h2bg <- 1 - h2s
        }
        if (h2bg + h2s != 1) {
            stop("Sum of the proportion of the variance of genetic effects is 
                 greater than 1; change h2s (fixed effect variance) or h2bg 
                 (random effect variance)\n h2s + h2bg = 1")
        }
    }
    if (modelNoise == "noNoise") {
        if (!is.null(noiseVar) && modelNoise == "noNoise") {
            stop("modelNoise is noNoise but noiseVar is unequal to zero")
        } else {
            noiseVar <- 0
        }
    } else {
        if (is.null(noiseVar)) {
            noiseVar <- 1 - genVar
        }
        if (all(c(is.null(delta), is.null(rho), is.null(phi)))) {
            if (modelNoise == "noiseFixedOnly") {
                delta <- 1
                phi <- 0
                rho <- 0
            } else if (modelNoise == "noiseBgOnly") {
                phi <- 1
                delta <- 0
                rho <- 0
            } else if (modelNoise == "noiseCorrelatedOnly") {
                rho <- 1
                phi <- 0
                delta <- 0
            } else {
                stop(paste("modelNoise is ", modelNoise, "but neither delta nor 
                           rho or phi are provided, at least one is required"))
            }
        } else if (all(c(!is.null(delta), !is.null(rho), !is.null(phi)))) {
            if(delta + rho + phi != 1) {
                stop("Sum of the proportion of the variance of noise effects is 
                     greater than 1; change delta (fixed effect variance), rho 
                     (correlated background variance) or phi (random effect 
                     variance)\n delta + rho + phi = 1")
            }
        } else if (all(c(!is.null(delta), !is.null(rho)))) { 
            phi <- 1 - delta - rho
        } else if (all(c(!is.null(delta), !is.null(phi)))) { 
            rho <- 1 - delta - phi
        } else if (all(c(!is.null(rho), !is.null(phi)))) { 
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
    }
    # proportions of genetic var: snp variance h2s and background variance h2bg
    # h2s + h2b = 1
    # genetic fixed
    if (grepl("Fixed", modelGenetic)) {
        if (is.null(genFixed)) stop("Genetic model includes fixed effects, but
                                    genFixed was not provided")
        if (! all(dim(genFixed$common) == c(N, P))) {
            dim_genFixed <- paste(dim(genFixed$common), collapse=",")
            stop("Dimensions of the genetic fixed effect (", dim_genFixed,") are
                different from the specified dimensions: number of columns P: ",
                 P, ", number of rows N: ", N) 
        }
        if (modelGenetic == "geneticFixedOnly") {
            if (h2s != 1) stop("Genetic model 'geneticFixedOnly' implies all 
                               genetic variance is explained by the fixed effect
                               , however h2s != 1")
        }
        if (is.null(genFixed$common) && theta != 0) {
            warning("No common fixed specific effect chosen, setting theta to 
                    0")
            theta <- 0
        }
        if (is.null(genFixed$specific && theta != 1)) {
            warning("No specific fixed genetic effect chosen, setting theta to 
                    1")
            theta <- 1
        }
        var_genFixed_comm <- theta * h2s * genVar
        var_genFixed_spec <- (1 - theta) * h2s * genVar

        genFixed_comm_rescaled <- rescaleVariance(genFixed$common, 
                                                  var_genFixed_comm)
        genFixed_spec_rescaled <- rescaleVariance(genFixed$specific, 
                                                  var_genFixed_spec)
    
        Y_genFixed <- addNonNulls(list(genFixed_comm_rescaled, 
                                       genFixed_spec_rescaled))
        colnames(Y_genFixed) <- paste(phenoID, seq(1, P, 1), sep="")
        rownames(Y_genFixed) <- paste(sampleID, seq(1, N, 1), sep="")
    } else {
        h2s <- 0
        var_genFixed_comm <- 0
        var_genFixed_spec <- 0
        Y_genFixed <- NULL
    }
    # genetic bg
    if (grepl("Bg", modelGenetic)) {
        if (is.null(genBg)) stop("Genetic model includes bg effects, but genBg 
                                 was not provided")
        if (! all(dim(genBg$common) == c(N, P))) {
            dim_genBg <- paste(dim(genBg$common), collapse=",")
            stop("Dimensions of the genetic bg effect (", dim_genBg,") are 
                 different from the specified dimensions: number of columns P: "
                 , P, ", number of rows N: ", N) 
        }
        if (modelGenetic == "geneticBgOnly" && h2bg != 1) {
            stop("Genetic model 'geneticBgOnly' implies all genetic variance is 
                 explained by the bg effect, however genVar != h2bg")
        }
        if (is.null(genBg$common) && eta != 0) {
            warning("No common background genetic effect chosen, setting eta to 
                    0")
            eta <- 0
        }
        if (is.null(genBg$specific) && eta != 1) {
            warning("No specific background genetic effect chosen, setting eta 
                    to 1")
            eta <- 1
        }
        var_genBg_comm <- eta * h2bg * genVar
        var_genBg_spec <- (1 - eta) * h2bg * genVar
        
        genBg_comm_rescaled <- rescaleVariance(genBg$common, var_genBg_comm)
        genBg_spec_rescaled <- rescaleVariance(genBg$specific, var_genBg_spec)
        
        Y_genBg <- addNonNulls(list(genBg_comm_rescaled, genBg_spec_rescaled))
        colnames(Y_genBg) <- paste(phenoID, seq(1, P, 1), sep="")
        rownames(Y_genBg) <- paste(sampleID, seq(1, N, 1), sep="")

        cov_Y_genBg <- t(Y_genBg) %*% Y_genBg
        cov_Y_genBg <-  cov_Y_genBg/mean(diag(cov_Y_genBg))
        diag(cov_Y_genBg) <- diag(cov_Y_genBg) + 1e-4
    } else {
        h2bg <- 0
        var_genBg_comm <- 0
        var_genBg_spec <- 0
        Y_genBg <- NULL
        cov_Y_genBg <- NULL
    }
    # proportions of noise variance: fixed delta, correlated rho and noise bg
    # noiseVar = rho + delta + rest
    # noise fixed 
    if (grepl("Fixed", modelNoise)) {
        if (is.null(noiseFixed)) stop("Noise model includes fixed effects, 
                                      but noiseFixed was not provided")
        if (! all(dim(noiseFixed$common) == c(N, P))) {
            dim_noiseFixed <- paste(dim(noiseFixed$common), collapse=",")
            stop("Dimensions of the noise fixed effect (", dim_noiseFixed,") are
                 different from the specified dimensions: number of columns P: "
                 , P, ", number of rows N: ", N) 
        }
        if (modelNoise == "noiseFixedOnly" && delta != 1) {
             stop("Noise model 'noiseFixedOnly' implies all noise variance is 
                  explained by the fixed effect, however noiseVar != delta")
        }
        if ((modelNoise == "noiseFixedAndBg" && delta + phi != 1) || 
            (modelNoise == "noiseFixedAndCorrelated" && delta + rho != 1)) {
            stop(paste("Noise model is", modelNoise, "and sum of noise variance 
                       effects is greater > 0"))
        }
        if (is.null(noiseFixed$common) && gamma != 0) {
            warning("No common fixed noise effect chosen, setting gamma to 0")
            gamma <- 0
        }
        if (is.null(noiseFixed$specific) && gamma != 1) {
            warning("No specific fixed noise effect chosen, setting gamma to 1")
            gamma <- 1
        }
        var_noiseFixed_comm <- gamma * delta * noiseVar
        var_noiseFixed_spec <- (1 - gamma) * delta * noiseVar
        
        noiseFixed_comm_rescaled <- rescaleVariance(noiseFixed$common, 
                                                    var_noiseFixed_comm)
        noiseFixed_spec_rescaled <- rescaleVariance(noiseFixed$specific, 
                                                    var_noiseFixed_spec)
        
        Y_noiseFixed <- addNonNulls(list(noiseFixed_comm_rescaled, 
                                         noiseFixed_spec_rescaled))
        colnames(Y_noiseFixed) <- paste(phenoID, seq(1, P, 1), sep="")
        rownames(Y_noiseFixed) <- paste(sampleID, seq(1, N, 1), sep="")
    } else {
        delta <- 0
        var_noiseFixed_comm <- 0
        var_noiseFixed_spec <- 0
        Y_noiseFixed <- NULL
    }

    # correlated bg
    if (grepl("Correlated", modelNoise)) {
        if (is.null(correlatedBg)) {
            stop("Noise model includes correlated bg effects, but correlatedBg 
                 was not provided")
        }
        if (! all(dim(correlatedBg) == c(N, P))) {
            dim_correlatedBg <- paste(dim(correlatedBg$common), collapse=",")
            stop("Dimensions of the correlated bg effect (", dim_correlatedBg,")
                 are different from the specified dimensions: number of columns 
                 P: ", P, ", number of rows N: ", N) 
        }
        if (modelNoise == "noiseCorrelatedOnly" && rho != 1) {
            stop("Noise model 'noiseCorrelatedOnly' implies all noise variance 
                 is explained by correlated bg effects, but noiseVar != rho")
        }
        if (modelNoise == "noiseCorrelatedAndBg" && rho + phi != 1) {
            stop(paste("Noise model is", modelNoise, "and sum of noise variance 
                       effects is greater > 0"))
        }
        var_noiseCorrelated <- rho *  noiseVar
        correlatedBg_rescaled <- rescaleVariance(correlatedBg, 
                                                 var_noiseCorrelated)
        
        Y_correlatedBg <- correlatedBg_rescaled
        colnames(Y_correlatedBg) <- paste(phenoID, seq(1, P, 1), sep="")
        rownames(Y_correlatedBg) <- paste(sampleID, seq(1, N, 1), sep="")
    } else {
        rho <- 0
        var_noiseCorrelated <- 0
        Y_correlatedBg <- NULL 
    }

    # noise bg
    if (grepl("Bg", modelNoise)) {
        if (is.null(noiseBg)) {
            stop("Noise model includes bg effects, but noiseBg not provided")
        }
        if (! all(dim(noiseBg$common) == c(N, P))) {
            dim_noiseBg <- paste(dim(noiseBg$common), collapse=",")
            stop("Dimensions of the noise bg effect (", dim_noiseBg,") are 
                 different from the specified dimensions: number of columns P: "
                 , P, ", number of rows N: ", N) 
        }
        if (modelNoise == "noiseBgOnly" && phi != 1) {
            stop("Noise model 'noiseBgOnly' implies all noise variance is 
                 explained by the bg effect, however noiseVar != phi")
        }
        if (is.null(noiseBg$common) && alpha != 0) {
            warning("No common background noise effect chosen, setting alpha to 
                    0")
            alpha <- 0
        }
        if (is.null(noiseBg$specific) && alpha != 1) {
            warning("No specific background noise effect chosen, setting alpha 
                    to 1")
            alpha <- 1
        }
        var_noiseBg_comm <- alpha * phi * noiseVar
        var_noiseBg_spec <- (1 - alpha) * phi * noiseVar
        
        noiseBg_comm_rescaled <- rescaleVariance(noiseBg$common, 
                                                 var_noiseBg_comm)
        noiseBg_spec_rescaled <- rescaleVariance(noiseBg$specific, 
                                                 var_noiseBg_spec)
        
        Y_noiseBg <- addNonNulls(list(noiseBg_comm_rescaled, 
                                      noiseBg_spec_rescaled))
        colnames(Y_noiseBg) <- paste(phenoID, seq(1, P, 1), sep="")
        rownames(Y_noiseBg) <- paste(sampleID, seq(1, N, 1), sep="")
        cov_Y_noiseBg <- t(Y_noiseBg) %*% Y_noiseBg
        cov_Y_noiseBg <-  cov_Y_noiseBg/mean(diag(cov_Y_noiseBg))
        diag(cov_Y_noiseBg) <- diag(cov_Y_noiseBg) + 1e-4
    } else {
        var_noiseBg_comm <- 0
        var_noiseBg_spec <- 0
        Y_noiseBg <- NULL
        cov_Y_noiseBg <- NULL
    }
    
    vmessage(c("Put all phenotype components together..."), verbose=verbose)
    ### ii) put all components together and scale to mean=0, sd=1
    components <- list(Y_genFixed, Y_genBg, Y_noiseFixed, Y_noiseBg, 
                       Y_correlatedBg)
    Y <-  addNonNulls(components)
    Y <- scale(Y)
    colnames(Y) <- paste(phenoID, seq(1, P, 1), sep="")
    rownames(Y) <- paste(sampleID, seq(1, N, 1), sep="")

    varComponents <- data.frame(genVar=genVar, h2s=h2s, h2bg=h2bg, 
                                var_genFixed_comm=var_genFixed_comm, 
                                var_genFixed_spec=var_genFixed_spec, 
                                var_genBg_comm=var_genBg_comm, 
                                var_genBg_spec=var_genBg_spec, 
                                noiseVar=noiseVar, 
                                var_noiseFixed_comm=var_noiseFixed_comm, 
                                var_noiseFixed_spec=var_noiseFixed_spec, 
                                var_noiseBg_comm=var_noiseBg_comm, 
                                var_noiseBg_spec=var_noiseBg_spec, 
                                var_noiseCorrelated=var_noiseCorrelated)
    phenoComponents <- list(Y=Y, Y_genFixed=Y_genFixed, Y_genBg=Y_genBg, 
                            Y_noiseFixed=Y_noiseFixed, Y_noiseBg=Y_noiseBg, 
                            Y_correlatedBG=Y_correlatedBg, 
                            cov_Y_genBg=cov_Y_genBg, 
                            cov_Y_noiseBg=cov_Y_noiseBg, genFixed=genFixed, 
                            noiseFixed=noiseFixed)
    setup <- list(P=P, N=N, modelGenetic=modelGenetic, modelNoise=modelNoise, 
                  sampleID=sampleID, phenoID=phenoID)
    
    return(list(varComponents=varComponents, phenoComponents=phenoComponents, 
                setup=setup))
}

#' Run phenotype simulation.
#'
#' runSimulation wraps around the phenotype component functions (genFixedEffects
#' , genBgEffects, noiseBgEffects, noiseFixedEffects and correlatedBgEffects)
#' and combines the simulated phenotype components via create phenotype.
#'
#' @param P number [integer] of phenotypes to simulate 
#' @param N number [integer] of samples to simulate
#' @param tNrSNP total number [integer] of SNPs to simulate; these SNPs are used
#'  for kinship estimation
#' @param cNrSNP number [integer] of causal SNPs; used as genetic fixed effects
#' @param NrConfounders number [integer] of confounders; used as noise fixed 
#' effects
#' @param NrFixedEffects number [integer] of different fixed effects to simulate;
#' allows to simulate fixed effects from different distributions or with 
#'  different parameters
#' @param chr numeric vector of chromosomes to chose NrCausalSNPs from; only 
#' used when external genotype data is provided i.e. is.null(genoFilePrefix) == 
#' FALSE
#' @param chr_string alternative to chr, a [string] with chromosomes separated 
#' by comma; most often used when run as a command line application.
#' @param NrChrCausal Number [integer] of causal chromosomes to  chose 
#' NrCausalSNPs from (as opposed to the actual chromosomes to chose from via chr
#' 'chr_string);  only used when external genotype data is provided i.e. 
#' is.null(genoFilePrefix) == FALSE. 
#' @param SNPfrequencies vector of allele frequencies [double] from which to 
#' sample
#' @param SNPfrequencyString alternative to a frequencies vector, a [string] 
#' with frequencies separated by comma
#' @param genoFilePrefix path/to/chromosome-wise-genotype-file-ending-before-
#' "chrChromosomeNumber" [string]
#' @param genoFileSuffix [string] following chromosome number including 
#' .fileformat (e.g. ".csv")
#' @param genoFileDelimiter field separator [string] of genotype file
#' @param kinshipfile path/to/kinshipfile [string] to be read; either X or 
#' kinshipfile must be provided
#' @param kinshipDelimiter field separator [string] of kinship file 
#' @param kinshipHeader [boolean], if TRUE kinship file has header information 
#' @param normalise [boolean], if TRUE kinship matrix will be normalised by the 
#' mean of its diagonal elements
#' and 1e-4 added to the diagonal for numerical stability
#' @param standardise [boolean], if TRUE standardised genotypes will be returned
#' @param distConfounders name [string] of distribution to use to simulate 
#' confounders; one of "unif", "norm", "bin", "cat_norm", "cat_unif"
#' @param distBeta name [string] of distribution to use to simulate effect sizes
#'  of confounders; one of "unif" or "norm"
#' @param mConfounders mean/midpoint [double] of normal/uniform distribution for
#'  confounders
#' @param sdConfounders standard deviation/extension from midpoint [double] of 
#' normal/uniform distribution for confounders
#' @param catConfounders confounder categories [factor]; required if 
#' distConfounders "cat_norm" or "cat_unif" 
#' @param probConfounders probability [double] of binomial confounders (0/1); 
#' required if distConfounders "bin" 
#' @param mBeta mean/midpoint [double] of normal/uniform distribution for effect
#'  sizes of confounders
#' @param sdBeta standard deviation/extension from midpoint [double] of normal/
#' uniform distribution for effect sizes of confounders
#' @param pSpecificConfounders Proportion [double] of noise effects 
#' (confounders) to have a trait-specific effect
#' @param pTraitSpecificConfounders Proportion [double] of traits influenced by 
#' specific fixed noise effects
#' @param pSpecificGenetic Proportion [double] of genetic effects (SNPs) to have
#'  a trait-specific fixed effect
#' @param pTraitSpecificGenetic Proportion [double] of traits influenced by 
#' specific fixed genetic effects
#' @param NrConfoundersString alternative to NrConfounder, a comma-separated
#' [string] with number(s) [integer] of confounders to simulate; typically used 
#' when run as command line application
#' @param pSpecificConfoundersString alternative to pSpecificConfounders, a 
#' comma-separated [string] with proportion(s) [double] of noise effects 
#' (confounders) to have a trait-specific effect; typically used when run as 
#' command line application
#' @param pTraitSpecificConfoundersString alternative to 
#' pTraitSpecificConfounders,  a comma-separated [string] with proportion(s) 
#' [double] of traits influenced by specific fixed noise effects; typically used
#'  when run as command line application
#' @param distConfoundersString alternative to distConfounders, a comma-
#' separated [string] with name(s) [string] of distributions to use to 
#' simulate confounders; one of "unif", "norm", "bin", "cat_norm", "cat_unif";
#' typically used when run as command line application 
#' @param distBetaString alternative to distBeta, a comma-separated [string] 
#' with name(s) [string] of distribution to use to simulate effect sizes of 
#' confounders; one of "unif" or "norm"; typically used when run as command line 
#' application
#' @param mConfoundersString alternative to mConfounders, a comma-
#' separated [string] with of mean/midpoint(s) [double] of normal/uniform 
#' distribution for confounders; typically used when run as command line 
#' application
#' @param sdConfoundersString alternative to sdConfounders, a comma-
#' separated [string] with standard deviation(s)/distance from 
#' midpoint(s) [double] of normal/uniform distribution for confounders; 
#' typically used when run as command line application
#' @param catConfoundersString alternative to catConfounders, a comma-
#' separated [string] with the number of confounder categories [integer]; 
#' required if distConfounders "cat_norm" or "cat_unif"; 
#' typically used when run as command line application
#' @param probConfoundersString alternative to probConfounders, a comma-
#' separated [string] with probability(s) [double] of binomial 
#' confounders (0/1); required if distConfounders "bin"; typically used
#'  when run as command line application
#' @param mBetaString  alternative to mBeta, a comma- separated [string] with 
#' means/midpoints [double] of normal/uniform distribution for effect sizes of 
#' confounders; typically used when run as command line application
#' @param sdBetaString alternative to sdBeta, a comma- separated [string] with 
#' standard deviation/distance from midpoint [double] of normal/uniform 
#' distribution for effect sizes of confounders; typically used when run as 
#' command line application

#' @param meanNoiseBg mean [double] of the normal distribution for noise bg 
#' effects
#' @param sdNoiseBg standard deviation [double] of the normal distribution for 
#' noise bg effects
#' @param sampleID prefix [string] for naming samples (followed by sample number
#'  from 1 to N)
#' @param phenoID prefix [string] for naming traits (followed by trait number 
#' from 1 to P)
#' @param genVar Proportion [double] of total genetic variance
#' @param h2s Proportion [double] of gentic variance of fixed effects 
#' @param h2bg Proportion [double] of genetic variance of background effects; 
#' either h2s or h2bg have to be specified and 
#' h2s + h2bg = 1
#' @param theta Proportion [double] of variance of common fixed genetic effects
#' @param eta Proportion [doubl:w
#' e] of variance of common bg genetic effects
#' @param noiseVar Proportion [double] of total noise variance
#' @param rho Proportion [double] of noise variance of correlated effects; sum 
#' of rho, delta and phi cannot be greater than 1
#' @param pcorr Correlation [double] between phenotypes
#' @param delta Proportion [double] of noise variance of fixed effects; sum of 
#' rho, delta and phi cannot be greater than 1
#' @param gamma Proportion [double] of variance of common fixed noise effects
#' @param phi Proportion [double] of noise variance of background effects; sum 
#' of rho, delta and phi cannot be greater than 1
#' @param alpha Variance [double] of common bg noise effect
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
#' simulatedPhenotype <- runSimulation(N=100, P=10, genVar=genVar, h2bg=genVar, 
#' phi=1)
runSimulation <- function(N=1000, P=10, tNrSNP=5000, cNrSNP=20, 
                          NrConfounders=10, seed=219453, chr_string=NULL, 
                          chr=20, NrChrCausal=NULL,
                          genVar=NULL, h2s=NULL, theta=0.8, h2bg=NULL, eta=0.8, 
                          noiseVar=NULL, rho=NULL, delta=NULL, gamma=0.8, 
                          phi=NULL, alpha=0.8, sampleID="ID_", phenoID="Trait_", 
                          genoFilePrefix=NULL, genoFileSuffix=NULL, 
                          genoFileDelimiter=",", kinshipfile=NULL, 
                          kinshipHeader=TRUE, kinshipDelimiter=",", 
                          normalise=TRUE, standardise=TRUE,
                          NrFixedEffects=1, distConfounders="norm",
                          mConfounders=0, sdConfounders=1,
                          catConfounders=NULL, probConfounders=NULL,
                          distBeta="norm", mBeta=0, sdBeta=1,
                          pSpecificConfounders=0.4, 
                          pTraitSpecificConfounders=0.2, 
                          pcorr=0.8, meanNoiseBg=0, sdNoiseBg=1, 
                          SNPfrequencies=c(0.1, 0.2, 0.4), 
                          SNPfrequencyString=NULL,
                          pSpecificGenetic=0.4, pTraitSpecificGenetic=0.2,
                          NrConfoundersString=NULL, 
                          pSpecificConfoundersString=NULL, 
                          pTraitSpecificConfoundersString=NULL, 
                          distConfoundersString=NULL, 
                          mConfoundersString=NULL, 
                          sdConfoundersString=NULL, 
                          catConfoundersString=NULL, 
                          probConfoundersString=NULL, 
                          distBetaString=NULL, 
                          mBetaString=NULL, 
                          sdBetaString=NULL,
                          verbose=TRUE) {

    vmessage(c("Set seed:", seed))
    set.seed(seed)
       # find model
    model <- setModel(genVar=genVar, h2s=h2s, h2bg=h2bg, theta=theta, eta=eta, 
                      noiseVar=noiseVar, rho=rho, delta=delta, gamma=gamma, 
                      phi=phi, alpha=alpha, pcorr=pcorr, 
                      pSpecificConfounders=pSpecificConfounders,  
                      pTraitSpecificConfounders=pTraitSpecificConfounders, 
                      pSpecificGenetic=pSpecificGenetic, 
                      pTraitSpecificGenetic=pTraitSpecificGenetic, v=verbose)
    modelNoise <- model$modelNoise
    modelGenetic <- model$modelGenetic

    ### create simulated phenotypes
    # 1. Simulate noise terms
    if (modelNoise != 'noNoise') {
        vmessage(c("Simulate noise terms (noise model:", modelNoise, ")"),
                   verbose=verbose)
        if (grepl('Correlated', modelNoise))  {
            vmessage("Simulate correlated background effects", verbose=verbose)
            correlatedBg <- correlatedBgEffects(N=N, P=P, pcorr=pcorr)
        } else {
            correlatedBg <- NULL
        }
        if (grepl('Bg', modelNoise)) {
            vmessage("Simulate noise background effects", verbose=verbose)
            noiseBg <- noiseBgEffects(N=N, P=P, mean=meanNoiseBg, sd=sdNoiseBg)
        } else {
            noiseBg <- NULL
        }

        if (grepl('Fixed', modelNoise)) {
            vmessage("Simulate noise fixed effects", verbose=verbose)
            noiseFixed <- noiseFixedEffects(P=P, N=N, 
                                            NrFixedEffects = NrFixedEffects,
                                            NrConfounders=NrConfounders,
                                            pSpecificConfounders=
                                                pSpecificConfounders,
                                            pTraitSpecificConfounders=
                                                pTraitSpecificConfounders,
                                            distConfounders=distConfounders, 
                                            mConfounders=mConfounders, 
                                            sdConfounders=sdConfounders, 
                                            catConfounders=catConfounders, 
                                            probConfounders = probConfounders,
                                            distBeta=distBeta, mBeta=mBeta, 
                                            sdBeta=sdBeta,
                                            NrConfoundersString=
                                                NrConfoundersString,
                                            pSpecificConfoundersString=
                                                pSpecificConfoundersString,
                                            pTraitSpecificConfoundersString=
                                                pTraitSpecificConfoundersString,
                                            distConfoundersString=
                                                distConfoundersString, 
                                            mConfoundersString=
                                                mConfoundersString, 
                                            sdConfoundersString=
                                                sdConfoundersString, 
                                            catConfoundersString=
                                                catConfoundersString, 
                                            probConfoundersString=
                                                probConfoundersString,
                                            distBetaString=distBetaString, 
                                            mBetaString=mBetaString, 
                                            sdBetaString=sdBetaString)
        } else {
            noiseFixed <- NULL
        }
    }
    
    # 2. Simulate genetic terms
    if ( modelGenetic != 'noGenetic') {
        vmessage(c("Simulate genetic effects (genetic model:", modelGenetic,")")
                 , verbose=verbose)
        if (is.null(genoFilePrefix)) {
            genotypes <- simulateGenotypes(N=N, NrSNP=tNrSNP, 
                                           frequencies=SNPfrequencies, 
                                           sampleID=sampleID, 
                                           frequencyString=SNPfrequencyString, 
                                           verbose=verbose)
            kinship <- getKinship(X=genotypes$X_sd, sampleID=sampleID, 
                                  norm=normalise, verbose=verbose)
            causalSNPs <- getCausalSNPs(NrCausalSNPs=cNrSNP, genotypes=genotypes,
                                        sampleID=sampleID, 
                                        standardise=standardise, 
                                        verbose=verbose)
        } else {
            kinship <- getKinship(kinshipfile=kinshipfile, sampleID=sampleID, 
                                  norm=normalise, sep=kinshipDelimiter, 
                                  header=kinshipHeader, verbose=verbose)
            causalSNPs <- getCausalSNPs(NrCausalSNPs=cNrSNP, chr=chr, 
                                        chr_string=chr_string, 
                                        genoFilePrefix=genoFilePrefix, 
                                        genoFileSuffix=genoFileSuffix, 
                                        genoFileDelimiter=genoFileDelimiter, 
                                        sampleID=sampleID, 
                                        standardise=standardise, 
                                        verbose=verbose)
        }   

        if (grepl('Fixed', modelGenetic)) {
            vmessage("Simulate genetic fixed effects", verbose=verbose)
            genFixed <- geneticFixedEffects(X_causal=causalSNPs, N=N, P=P, 
                                            pSpecificGenetic=pSpecificGenetic, 
                                            pTraitSpecificGenetic=
                                                pTraitSpecificGenetic)
        } else {
            genFixed <- NULL
        }

        if (grepl('Bg', modelGenetic)) {
            vmessage("Simulate genetic background effects (kinship-based)"
                     , verbose=verbose)
            genBg <- geneticBgEffects(N=N, P=P, kinship=kinship)
        } else {
            genBg <- NULL
        }
    } else {
        genotypes <- NULL
        kinship <- NULL
        causalSNPs <- NULL
    }
    
    # 3. Construct final simulated phenotype 
    vmessage("Construct final simulated phenotype"
             , verbose=verbose)
    finalPheno <- createPheno(N=N, P=P, 
                            correlatedBg=correlatedBg, genFixed=genFixed, 
                            genBg=genBg, noiseFixed=noiseFixed, noiseBg=noiseBg, 
                            genVar=genVar, h2s=h2s, h2bg=h2bg, theta=theta, 
                            eta=eta, noiseVar=noiseVar, delta=delta, 
                            gamma=gamma, phi=phi, alpha=alpha, rho=rho,  
                            modelNoise=modelNoise, modelGenetic=modelGenetic, 
                            verbose=verbose ) 
}

#' Save final phenotype and phenotype components.
#'
#' savePheno saves model setup parameters and simulated genotypes to the 
#' specified directories. Requires a simulatedData list which is the output of 
#' either \link{runSimulation} or \link{createPheno}
#'
#' @param simulatedData named list of i) dataframe of proportion of variance 
#' explained for each component (varComponents), 
#' ii) a named list with the simulated phenotype components (phenoComponents) 
#' and iii) a named list of parameters describing the model setup (setup);
#' obtained from either \link{runSimulation} or \link{createPheno} 
#' @param directoryGeno name of parent directory [string] where genotypes from 
#' simulations should be saved [needs user writing permission]
#' @param directoryPheno name of parent directory [string] where final phenotype 
#' and phenotype components from simulations should be saved [needs user writing 
#' permission]
#' @param outstring optional name [string] of subdirectory (in relation to 
#' directoryPheno/directoryGeno) to save set-up
#' specific simulation results
#' @param kinship [N x N] matrix of kinship estimate; optional: if provided, 
#' kinship estimate will be saved to file
#' @param genotypes [N x totalNrSNPs] matrix of all simulated SNPs; optional: if
#'  provided, simulated SNPs will be saved to file 
#' @param causalSNPs [N x NrCausalSNPs] matrix of causal SNPs; optional: if 
#' provided, causalSNPs will be saved to file
#' @param sample_subset_vec optional vector of sample subset sizes [integer] e.g.
#' if provided, draws subsets of samples out of the total simulated dataset and 
#' saves them separately 
#' @param pheno_subset_vec optional vector of phenotype subset sizes [integer] 
#' e.g.if provided, draws subsets of traits out of the total simulated dataset 
#' and saves them separately 
#' @param sample_subset_string optional [string] of comma-separated numbers e.g.
#'  "50,100,500"; alternative to \code{sample_subset_vec}, typically used when
#'  run as command line application
#' @param pheno_subset_string optional [string] of comma-separated numbers e.g. 
#' "10,50,80"; alternative to \code{pheno_subset_vec}, typically used when
#'  run as command line application
#' @param saveAsTable [boolean] if TRUE, data will be saved as .csv 
#' @param saveAsRDS [boolean] if TRUE, data will be saved as .RDS; at least one 
#' of 'saveAsTable' or 'saveAsRDS' has to be TRUE, both can be TRUE
#' @param verbose [boolean]; if TRUE, progress info is printed to standard out
#' @return no return value, sole purpose of function is writing output files
#' @export
#' @examples
#' simulatedPhenotype <- runSimulation(N=100, P=10, genVar=0.4, h2bg=1, 
#' phi=1, verbose=FALSE)
#' #not run
#' #savePheno(simulatedPhenotype, directoryGeno="/path/to/dir/",  
#' #directoryPheno="/path/to/dir/", outstring="Date_simulation")
savePheno <- function(simulatedData, directoryGeno, directoryPheno, 
                      kinship=NULL, genotypes=NULL, causalSNPs=NULL,
                      sample_subset_vec=NULL, pheno_subset_vec=NULL, 
                      sample_subset_string=NULL, pheno_subset_string=NULL, 
                      outstring=NULL, saveAsTable=TRUE, saveAsRDS=FALSE, 
                      verbose=TRUE) {
    if (!saveAsTable && !saveAsRDS) {
        stop("Either one of saveAsTable or saveAsRDS must be true in order to 
             save output")
    }
    modelGenetic <- simulatedData$setup$modelGenetic
    modelNoise <- simulatedData$setup$modelNoise
    N <- simulatedData$setup$N
    P <- simulatedData$setup$P
    sampleID <-  simulatedData$setup$sampleID
    phenoID <-  simulatedData$setup$phenoID
    NrSNP <- ncol(simulatedData$phenoComponents$causalSNPs)
    genVar <- simulatedData$varComponents$genVar

    if (! is.null(sample_subset_string) || ! is.null(sample_subset_vec)) {
        vmessage(c("Create sample subsets:", sample_subset_string), 
                 verbose=verbose)
        if (! is.null(sample_subset_string)) {
            ssample <- commaList2vector(sample_subset_string)
        } else {
            ssample <- sample_subset_vec
        }
        sample_subset <- sapply(ssample, function(s) {
            tmp <- sample(1:N, s, replace=FALSE)
            names(tmp) <- paste(sampleID, tmp, sep="")
            return(tmp)
        })
    } else {
        sample_subset <- list(set=seq(1, N, 1))
    }

    if (! is.null(pheno_subset_string) || ! is.null(pheno_subset_vec)) {
        vmessage(c("Create pheno subsets:", pheno_subset_string), 
                 verbose=verbose)
        if (! is.null(pheno_subset_string)) {
            psample <- commaList2vector(pheno_subset_string)
        } else {
            psample <- pheno_subset_vec
        }
        pheno_subset <- sapply(psample, function(s) {
            tmp <- sample(1:P, s, replace=FALSE)
            names(tmp) <- paste(phenoID, tmp, sep="")
            return(tmp)
        })
    } else {
        pheno_subset <- list(set=seq(1, P, 1))
    }
    
    ### set-up directories
    if (is.null(outstring)) {
        outstring=paste("samples", N, "_NrSNP", NrSNP, "_Cg", genVar, "_model", 
                        modelNoise, modelGenetic, sep="")
    }

    directoryGeno <- file.path(directoryGeno, outstring)
    ifelse(!dir.exists(directoryGeno), 
           dir.create(directoryGeno, recursive=TRUE), FALSE)

    directoryPheno <- file.path(directoryPheno, outstring)
    ifelse(!dir.exists(directoryPheno), 
           dir.create(directoryPheno, recursive=TRUE), FALSE)
    
    vmessage(c("Save simulation results"), verbose=verbose)
    out <- l_ply(sample_subset, function(ss) {
            l_ply(pheno_subset, function(sp, ss) {
    
            nrsamples <- length(ss)
            nrpheno <- length(sp)

            outstring=paste("samples", nrsamples, "_traits", nrpheno, "_NrSNP",
                            NrSNP, "_Cg", genVar, "_model", modelNoise, 
                            modelGenetic, sep="")
            directoryPheno = paste(directoryPheno,"/", "samples", nrsamples, 
                                   "_NrSNP",NrSNP, "_Cg", genVar, "_model", 
                                   modelNoise, modelGenetic, sep="")
            
            ifelse(!dir.exists(directoryPheno), 
                   dir.create(directoryPheno, recursive=TRUE), FALSE)

            vmessage(c("Save phenotype to ", directoryPheno, "/Y..."), 
                     verbose=verbose, sep="")
            subset_Y <- simulatedData$phenoComponents$Y[ss,sp]
            if (saveAsRDS) {
                saveRDS(subset_Y, 
                        paste(directoryPheno, "/Ysim_", outstring ,".rds", 
                              sep=""))
            }
            if (saveAsTable) {
                write.table(subset_Y, paste(directoryPheno, "/Ysim_",outstring,
                            ".csv", sep=""), sep=",", quote=FALSE,
                            col.names=NA, row.names=TRUE)
            }
            if (grepl("Bg", modelGenetic)) {
                subset_genBg <- simulatedData$phenoComponents$Y_genBg[ss,sp]
                subset_cov_Y_genBg <- 
                    simulatedData$phenoComponents$cov_Y_genBg[sp,sp]
                vmessage(c("Save genetic background to ", directoryPheno, 
                           "/Y_genBg..."), verbose=verbose, sep="")
                if (saveAsRDS) saveRDS(subset_genBg, paste(directoryPheno, 
                                       "/Y_genBg_", outstring,".rds", sep=""))
                if (saveAsRDS) saveRDS(subset_cov_Y_genBg, paste(directoryPheno, 
                                       "/cov_Y_genBg_", outstring,".rds", 
                                       sep=""))
                if (saveAsTable) write.table(subset_genBg, paste(directoryPheno, 
                                             "/Y_genBg_", outstring,".csv", 
                                             sep=""),
                                             sep=",",quote=FALSE, col.names=NA, 
                                             row.names=TRUE)
                if (saveAsTable) write.table(subset_cov_Y_genBg, 
                                 paste(directoryPheno, "/cov_Y_genBg_", 
                                       outstring, ".csv", sep=""), sep=",",
                                 quote=FALSE, col.names=FALSE, row.names=FALSE)

                if(!is.null(kinship)) {
                    vmessage(c("Save kinship to", directoryGeno), 
                             verbose=verbose)
                    write.table(kinship, paste(directoryGeno, "/samples", N, 
                              "_NrSNP", NrSNP, "_kinship.csv", sep=""), sep=",",
                              col.names=TRUE, row.names=FALSE)
                }
            }
            if (grepl("Fixed", modelGenetic)) {
                subset_genFixed <- 
                    simulatedData$phenoComponents$Y_genFixed[ss,sp]
                vmessage(c("Save genetic fixed effects to ", directoryPheno, 
                           "/Y_genFixed..."), verbose=verbose, sep="")
                if (saveAsRDS) {
                    saveRDS(subset_genFixed, 
                            paste(directoryPheno, "/Y_genFixed_", outstring,
                                  ".rds", sep=""))
                }
                if (saveAsTable) {
                    write.table(subset_genFixed, 
                                paste(directoryPheno, "/Y_genFixed_", outstring,
                                      ".csv", sep=""), sep=",", quote=FALSE, 
                                col.names=NA, row.names=TRUE)
                }

                vmessage(c("Save SNPs and effect sizes to ", directoryGeno, 
                           "/SNP..."), verbose=verbose, sep="")
                SNP <- simulatedData$phenoComponents$genFixed$cov[,ss]
                colnames(SNP) <- rownames(subset_Y)
                SNP_effect <- 
                    simulatedData$phenoComponents$genFixed$cov_effect[sp,]
                colnames(SNP_effect) <- colnames(subset_Y)
                if (saveAsRDS) {
                    saveRDS(SNP,  
                            paste(directoryGeno,"/SNP_NrSNP", NrSNP, "_",  
                                outstring,".rds",sep=""))
                }
                if (saveAsRDS) {
                    saveRDS(SNP_effect,  
                            paste(directoryGeno, "/SNP_effects_NrSNP", NrSNP, 
                                  "_", outstring, ".rds",sep=""))
                }
                if (saveAsTable) {
                    write.table(SNP,
                                paste(directoryGeno, "/SNP_NrSNP", NrSNP, "_",  
                                        outstring, ".csv",sep=""), sep=",", 
                                        col.names=TRUE, row.names=FALSE, 
                                        quote=FALSE)
                }
                if (saveAsTable) {
                    write.table(SNP_effect,  
                                paste(directoryGeno, "/SNP_effects_NrSNP", 
                                      NrSNP, "_", outstring, ".csv",sep=""), 
                                sep=",", col.names=TRUE, row.names=FALSE, 
                                quote=FALSE)
                }

                if (!is.null(genotypes)) {
                    vmessage(c("Save genotypes to", directoryGeno), 
                             verbose=verbose)
                    samples <-paste(sampleID, seq(1,N,1), sep="")
                    X_id <- data.frame(FID=samples, IID=samples, PAT=rep(0,N), 
                                       MAT=rep(0,N), SEX=rep(0,N), 
                                       PHENOTYPE=rep(-9, N))
                    write.table(cbind(X_id, X_id), paste(directoryGeno, 
                                "/samples", N, "_NrSNP", NrSNP, 
                                "_genotypes_ID.txt", sep=""), sep="\t", 
                                col.names=TRUE, row.names=FALSE)
                    write.table(genotypes, paste(directoryGeno, "/samples", N, 
                              "_NrSNP", NrSNP, "_genotypes.txt", sep=""), 
                              sep="\t", col.names=TRUE, row.names=FALSE)
                }
                if (!is.null(causalSNPs)) {
                    write.table(causalSNPs, paste(directoryGeno, "/NrTraits", P,
                                "_causalSNPs.csv", sep=""), sep=",", 
                                col.names=NA, row.names=TRUE, quote=FALSE)
                }
            }

            if (grepl("Correlatd", modelNoise)) {
                subset_correlatedBg <- 
                    simulatedData$phenoComponents$Y_correlatedBg[ss,sp]
                vmessage(c("Save correlated background to ", directoryPheno, 
                           "/Y_correlatedBg..."), verbose=verbose, sep="")
                if (saveAsRDS) saveRDS(subset_correlatedBg, paste(directoryPheno
                                       , "/Y_correlatedBg_", outstring,".rds", 
                                       sep=""))
                if (saveAsTable) write.table(subset_correlatedBg,
                                             paste(directoryPheno, 
                                                   "/Y_correlatedBg_",
                                                   outstring, ".csv", sep=""), 
                                             sep=",", quote=FALSE, col.names=NA, 
                                             row.names=TRUE)
            }
            if (grepl("Bg", modelNoise)) {
                subset_noiseBg <- simulatedData$phenoComponents$Y_noiseBg[ss,sp]
                subset_cov_Y_noiseBg <- 
                    simulatedData$phenoComponents$cov_Y_noiseBg[sp,sp]
                vmessage(c("Save noise background to ", directoryPheno, 
                           "/Y_noiseBg..."), verbose=verbose, sep="")
                if (saveAsRDS) saveRDS(subset_noiseBg, paste(directoryPheno, 
                                       "/Y_noiseBg_", outstring,".rds", sep=""))
                if (saveAsRDS) saveRDS(subset_cov_Y_noiseBg, 
                                       paste(directoryPheno, "/cov_Y_noiseBg_", 
                                             outstring,".rds", sep=""))
                if (saveAsTable) write.table(subset_noiseBg, 
                                       paste(directoryPheno, "/Y_noiseBg_", 
                                             outstring,".csv", sep=""), sep=",",
                                       quote=FALSE, col.names=NA, 
                                       row.names=TRUE)
                if (saveAsTable) write.table(subset_cov_Y_noiseBg, 
                                    paste(directoryPheno, "/cov_Y_noiseBg_",
                                           outstring,".csv", sep=""), sep=",",
                                    quote=FALSE, col.names=FALSE, 
                                    row.names=FALSE)
            }

            if (grepl("Fixed", modelNoise)) {
                subset_noiseFixed <- 
                    simulatedData$phenoComponents$Y_noiseFixed[ss,sp]
                vmessage(c("Save noise fixed effects to ", directoryPheno, 
                           "/Y_noiseFixed..."), verbose=verbose, sep="")
                if (saveAsRDS) {
                    saveRDS(subset_noiseFixed, 
                            paste(directoryPheno, "/Y_noiseFixed_", outstring,
                                  ".rds", sep=""))
                }
                if (saveAsTable) {
                    write.table(subset_noiseFixed, 
                                paste(directoryPheno, "/Y_noiseFixed_", 
                                      outstring,".csv", sep=""), sep=",", 
                                quote=FALSE, col.names=NA, row.names=TRUE)
                }
                vmessage(c("Save covariates and effect sizes to ", 
                           directoryPheno, "/Covs..."), verbose=verbose, sep="")
                cov <- simulatedData$phenoComponents$noiseFixed$cov[,ss]
                colnames(cov) <- rownames(subset_Y)
                cov_effect <- 
                    simulatedData$phenoComponents$noiseFixed$cov_effect[sp,]
                if (saveAsRDS) {
                    saveRDS(cov, paste(directoryPheno, "/Covs_", 
                                        outstring, ".rds", sep=""))
                }
                if (saveAsRDS) {
                    saveRDS(simulatedData$phenoComponents$noiseFixed$cov_effects,
                            paste(directoryPheno, "/Covs_effect_", outstring,
                                  ".rds", sep=""))
                }
                if (saveAsTable) {
                    write.table(cov, paste(directoryPheno, "/Covs_",
                                            outstring, ".csv", sep=""), sep=",",
                                            quote=FALSE, col.names=NA, 
                                            row.names=TRUE)
                }
                if (saveAsTable) {
                    write.table(
                        simulatedData$phenoComponents$noiseFixed$cov_effects, 
                        paste(directoryPheno, "/Covs_effect_", outstring,".csv",
                              sep=""), sep=",", quote=FALSE, col.names=TRUE, 
                        row.names=FALSE)
                }
            }

            if (saveAsRDS) {
                saveRDS(simulatedData$varComponents, 
                        paste(directoryPheno, "/varComponents_",  outstring,
                              ".rds", sep=""))
            }
            if (saveAsTable) {
                write.table(simulatedData$varComponents, 
                         paste(directoryPheno, "/varComponents_", outstring,
                               ".csv", sep=""), sep=",", quote=FALSE,
                         col.names=TRUE, row.names=FALSE)
            }
        }, ss =ss)
    })
}


