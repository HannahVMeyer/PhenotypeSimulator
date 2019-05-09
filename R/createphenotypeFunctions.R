#' Scale phenotype component.
#'
#' The function scales the specified component such that the average column 
#' variance is equal to the user-specified proportion of variance. 
#'
#' @param component [N x P] Phenotype matrix [double] where [N] are the number 
#' of samples and P the number of phenotypes 
#' @param propvar Number [double] specifying the proportion of variance that 
#' should be explained by this phenotype component
#' @return If propvar != 0, a named list with the [N x P] matrix of the scaled 
#' component (component) and its scale factor [double] (scale_factor) else 
#' returns NULL
#' @export
#' @examples
#' x <- matrix(rnorm(100), nc=10)
#' x_scaled <- rescaleVariance(x, propvar=0.4)
rescaleVariance <- function(component, propvar) {
    if (!is.numeric(propvar)) {
        stop("propvar needs to be numeric")
    }
    if (propvar < 0 || propvar > 1) {
        stop("propvar cannot be less than 0 and or greater than 1")
    }
    if (!is.null(component) && !is.matrix(component)){
        stop("component needs to be a matrix")
    }
    if (!is.null(component) && propvar != 0) {
        var_component <- var(component)
        mean_var <- mean(diag(var_component))
        scale_factor <- sqrt(propvar/mean_var)
        component_scaled <- component * scale_factor
        colnames(component_scaled) <- colnames(component)
        rownames(component_scaled) <- rownames(component)
        return(list(component=component_scaled,
                    scale_factor=scale_factor))
    } else {
        return(NULL)
    }
}


#' Phenotype transformation.
#' 
#' Transformation of phenotype component by applying a user-specified 
#' non-linear transformation to the phenotype component.
#' 
#' @param component [N x P] Phenotype matrix [double] where [N] are the number 
#' of samples and P the number of phenotypes 
#' @param method [string] one of exp (exponential), log (logarithm), poly
#' (polynomial), sqrt (squareroot) or custom (user-supplied function)
#' @param transformNeg [string] one of abs (absolute value) or set0 (set all 
#' negative values to zero). If method==log and transformNeg==set0, negative
#' values set to 1e-5
#' @param alpha [double] weighting scalar for non-linearity: alpha==0 fully
#' linear phenotype, alpha==1 fully non-linear phenotype. See @details.
#' @param logbase [int] when method==log, sets the log base for transformation
#' @param expbase [double] when method==exp, sets the exp base for
#' transformation.
#' @param power [double] when method==poly, sets the power to raise to.
#' @param f [function] function accepting component as a single argument.
#' @param verbose [boolean]; If TRUE, progress info is printed to standard out.
#' @details transformNonlinear takes a phenotype component as input and 
#' transforms it according to the specified transformation method. The user can 
#' choose how strongly non-linear the resulting phenotype component should be, 
#' by specifying the weighting parameter alpha:
#' component_transformed = (1 - alpha) \* component + 
#' alpha \* transformfunction(component)
#' @return [N x P] transformed phenotype matrix [double]
#' @export
#' @examples 
#' # Simulate non-genetic covariate effects 
#' cov_effects <- noiseFixedEffects(N=100, P=5)
#' # Transform logarithmically
#' covs_log <- transformNonlinear(cov_effects$shared, alpha=0.5, method="log",
#' transformNeg="abs")
#' # Transform custom
#' f_custom <- function(x) {x^2 + 3*x}
#' covs_custom <- transformNonlinear(cov_effects$shared, alpha=0.5, 
#' method="custom", f=f_custom)

transformNonlinear <- function(component, alpha, method, logbase=10, power=2, 
                               expbase=NULL, transformNeg="abs", f=NULL,
                               verbose=TRUE) {
    testNumerics(numbers=c(alpha, expbase, logbase), proportions=alpha, 
                 positives=logbase)
    nonlinear <- function(x, method, expbase, logbase, power, f) {
        if (!is.null(transformNeg)) {
            if (transformNeg == "abs") {
                x <- abs(x)
            } else if (transformNeg == "set0") {
                if (method == "log") {
                    x[x < 0] <- 1e-5
                } else {
                    x[x < 0] <- 0
                }
            } else {
                stop("Negative transformation method not known")
            }
        }
        vmessage(c("Use", method, "as transformation method"), verbose=verbose)
        if (method == "exp") {
            if (is.null(expbase)) y <- exp(x)
            if (!is.null(expbase)) y <- expbase^x
        } else if (method == "log") {
            if (any(x <= 0)) {
                stop(paste("For transformation log, all values have to", 
                     " be greater than zero"))
            }
            y <- log(x, base=logbase)
        } else if (method == "poly") {
            y <- x^power
        } else if (method == "sqrt") {
            if (any(x < 0)) {
                stop(paste("For transformation square root, all values have to", 
                     "be greater or equal to zero"))
            }
            y <- sqrt(x)
        } else if (method == "custom") {
            y <- f(x)
        } else {
            stop("Transformation method not known")
        }
        return(y)
    }
    component_trans <- (1 - alpha) * component + 
        alpha * nonlinear(component, method=method, logbase=logbase, 
                          expbase=expbase, power=power, f=f)
    return(component_trans)
}

#' Set simulation model.
#'
#' Based on parameters provided, this function sets the name for the phenotype 
#' simulation. It carries out compatibiltiy checks of the specifie parameters 
#' and checks for any missing information. 
#' 
#' @param genVar Total genetic variance [double].
#' @param h2s Proportion [double] of variance of genetic variant effects.
#' @param h2bg Proportion [double] of variance of infinitesimal genetic effects
#' i.e. correlation introduced by sample kinship).
#' @param theta Proportion [double] of variance of shared genetic variant 
#' effects.
#' @param eta Proportion [double] of variance of shared infinitesimal genetic 
#' effects.
#' @param noiseVar Total noise variance [double]. 
#' @param rho Proportion [double] of variance of correlated noise effects.
#' @param pcorr Correlation [double] between phenotypes.
#' @param delta Proportion [double] of variance of non-genetic covariate effect.
#' @param gamma Proportion [double] of variance of shared non-genetic covariate 
#' effects.
#' @param phi Proportion [double] of variance of observational noise effects.
#' @param alpha Proportion [double] of variance of shared observational noise 
#' effect.
#' @param pIndependentConfounders Proportion [double] of non-genetic covariate 
#' to have a trait-independent effect.
#' @param pTraitIndependentConfounders Proportion [double] of traits influenced  
#' by independent non-genetic covariate effects.
#' @param pIndependentGenetic Proportion [double] of genetic variant effects to 
#' have a trait-independent fixed effect.
#' @param pTraitIndependentGenetic Proportion [double] of traits influenced by 
#' independent genetic variant effects.
#' @param proportionNonlinear [double] proportion of the phenotype to be non-
#' linear
#' @param cNrSNP Number [integer] of causal SNPs; used as genetic variant 
#' effects.
#' @param NrConfounders Number [integer] of non-genetic covariates; used as 
#' non-genetic covariate effects.
#' @param verbose [boolean]; If TRUE, progress info is printed to standard out.
#' @return Named list containing the genetic model (modelGenetic), the noise 
#' model (modelNoise) and the input parameters (h2s, h2bg, noiseVar, rho, delta, 
#' phi, gamma, theta, eta, alpha, pcorr, proportionNonlinear). Model options 
#' are: modelNoise: "noNoise", "noiseFixedOnly", "noiseBgOnly", 
#' "noiseCorrelatedOnly", "noiseFixedAndBg","noiseCorrelatedAndBg", 
#' "noiseFixedAndCorrelated", "noiseFixedAndBgAndCorrelated"
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
                     pTraitIndependentGenetic=0.2, proportionNonlinear=0,
                     cNrSNP=NULL, NrConfounders=10,
                     verbose=TRUE)  {
    if (is.null(c(genVar, noiseVar, h2bg, h2s, delta, rho, phi))) {
        stop("No variance components specified")
    }
    if (!is.null(NrConfounders) && all(NrConfounders != 0)) {
        # for testing purposes set to anything other integer than 0
        NrConfounders=10
    }
    numbers <- list(genVar=genVar, h2s=h2s, h2bg=h2bg, theta=theta,
                    eta=eta, noiseVar=noiseVar, delta=delta, gamma=gamma, 
                    rho=rho, phi=phi, alpha=alpha, pcorr=pcorr, 
                    pIndependentConfounders=pIndependentConfounders, 
                    pTraitIndependentConfounders=
                        pTraitIndependentConfounders, 
                    pIndependentGenetic=pIndependentGenetic, 
                    pTraitIndependentGenetic=pTraitIndependentGenetic,
                    proportionNonlinear=proportionNonlinear,
                    cNrSNP=cNrSNP, NrConfounders=NrConfounders)
    proportions <- list(genVar=genVar, h2s=h2s, h2bg=h2bg, theta=theta,
                        eta=eta, noiseVar=noiseVar, delta=delta, gamma=gamma, 
                        rho=rho, phi=phi, alpha=alpha, pcorr=pcorr, 
                        pIndependentConfounders=pIndependentConfounders, 
                        pTraitIndependentConfounders=
                            pTraitIndependentConfounders, 
                        pIndependentGenetic=pIndependentGenetic, 
                        pTraitIndependentGenetic=pTraitIndependentGenetic,
                        proportionNonlinear=proportionNonlinear)
    testNumerics(numbers=numbers, proportions=proportions)
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
    vmessage(c("The total noise variance (noiseVar) is:", noiseVar), 
             verbose=verbose)
    if ((noiseVar) == 0 ) {
        modelNoise="noNoise"
        delta <- 0 
        rho <- 0
        phi <- 0
        vmessage(c("The noise model is:", modelNoise), verbose=verbose)
    } else {
        if (all(c(is.null(delta), is.null(rho), is.null(phi)))) {
            stop(paste("Neither delta (non-genetic covariate effect variance)", 
                       "nor rho (correlated noise effect variance) or phi", 
                       "(observational noise effect variance) are provided, at", 
                       "least one is required")
            )
        } 
        if (length(c(delta, rho, phi)) >=2) {
            if (sum(c(delta, rho, phi)) > 1 ) {
                stop(paste("Sum of the proportion of the variance of noise",
                           "effects is greater than 1; change noiseVar, delta",
                           "(non-genetic covariate effect variance),",  
                           "rho (correlated noise effect variance) or phi", 
                           "(observational noise effect variance) such that", 
                           "delta + rho + phi = 1")
                )
            }
            if (length(c(delta, rho, phi)) == 3 && 
                sum(c(delta, rho, phi)) < 1) {
                stop(paste("Sum of the proportion of the variance of noise",
                           "effects is less than 1; change noiseVar, delta",
                           "(non-genetic covariate effect variance),",  
                           "rho (correlated noise effect variance) or phi", 
                           "(observational noise effect variance) such that", 
                           "delta + rho + phi = 1")
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
                       "needs to be set to 1; otherwise, the proportions of",
                       "variance of at least 2 components need to be specified")
            )
        }
        if (delta != 0 && NrConfounders == 0) {
            stop(paste("Proportion of of non-genetic covariate variance ",
                       "(delta) is", delta, "but number of",
                       "NrConfounder is set to zero"))
        }
        if (gamma == 1) {
            pIndependentConfounders=0
            pTraitIndependentConfounders=0
        }
        if (gamma == 0) {
            pIndependentConfounders=1
            pTraitIndependentConfounders=1
        }
        if (phi == 1) {
            modelNoise="noiseBgOnly"
            vmessage(c("The noise model is:", modelNoise), verbose=verbose)
            vmessage(c("Proportion of observational noise variance (phi):", 
                       phi), verbose=verbose)
            vmessage(c("Variance of shared observational noise effect (alpha):", 
                       alpha), verbose=verbose)
        } else if (delta == 1) {
            modelNoise="noiseFixedOnly"
            vmessage(c("The noise model is:", modelNoise), verbose=verbose)
            vmessage(c("Proportion of non-genetic covariate variance (delta):", 
                       delta), verbose=verbose)
            vmessage(c("Proportion of variance of shared non-genetic covariate",
                       "effects (gamma):", gamma), verbose=verbose)
            vmessage(c("Proportion of non-genetic covariates to have",
                       "a trait-independent effect (pIndependentConfounders"
                       , "):", pIndependentConfounders), 
                     verbose=verbose)
            vmessage(c("Proportion of traits influenced by independent",
                       "non-genetic covariate effects", 
                       "(pTraitIndependentConfounders):", 
                       pTraitIndependentConfounders), verbose=verbose)
        } else if (rho == 1) {
            modelNoise="noiseCorrelatedOnly"
            vmessage(c("The noise model is:", modelNoise), verbose=verbose)
            vmessage(c("Proportion of variance of correlated noise effects",
                       "(rho):", rho), verbose=verbose)
        } else if (rho == 0) {
            modelNoise="noiseFixedAndBg"
            vmessage(c("The noise model is:", modelNoise), verbose=verbose)
            vmessage(c("Proportion of non-genetic covariate variance (delta):", 
                       delta), verbose=verbose)
            vmessage(c("Proportion of variance of shared non-genetic covariate",
                       "effects (gamma):", gamma), verbose=verbose)
            vmessage(c("Proportion of non-genetic covariates to have",
                       "a trait-independent effect (pIndependentConfounders"
                       , "):", pIndependentConfounders), 
                     verbose=verbose)
            vmessage(c("Proportion of traits influenced by independent",
                       "non-genetic covariate effects", 
                       "(pTraitIndependentConfounders):", 
                       pTraitIndependentConfounders), verbose=verbose)
            vmessage(c("Proportion of observational noise variance (phi):", 
                       phi), verbose=verbose)
            vmessage(c("Variance of shared observational noise effect (alpha):", 
                       alpha), verbose=verbose)
        } else if (phi == 0) {
            modelNoise="noiseFixedAndCorrelated"
            vmessage(c("The noise model is:", modelNoise), verbose=verbose)
            vmessage(c("Proportion of non-genetic covariate variance (delta):", 
                       delta), verbose=verbose)
            vmessage(c("Proportion of variance of shared non-genetic covariate",
                       "effects (gamma):", gamma), verbose=verbose)
            vmessage(c("Proportion of non-genetic covariates to have",
                       "a trait-independent effect (pIndependentConfounders"
                       , "):", pIndependentConfounders), 
                     verbose=verbose)
            vmessage(c("Proportion of traits influenced by independent",
                       "non-genetic covariate effects", 
                       "(pTraitIndependentConfounders):", 
                       pTraitIndependentConfounders), verbose=verbose)
            vmessage(c("Proportion of variance of correlated noise effects",
                       "(rho):", rho), verbose=verbose)
            vmessage(c("Correlation between phenotypes (pcorr):", pcorr), 
                     verbose=verbose)
        } else if (delta == 0 ) {
            modelNoise="noiseCorrelatedAndBg"
            vmessage(c("The noise model is:", modelNoise), verbose=verbose)
            vmessage(paste("Proportion of variance of correlated noise effects",
                           "(rho):", rho), verbose=verbose)
            vmessage(c("Proportion of observational noise variance (phi):", 
                       phi), verbose=verbose)
            vmessage(c("Variance of shared observational noise effect (alpha):", 
                       alpha), verbose=verbose)
        } else {
            modelNoise="noiseFixedAndBgAndCorrelated"
            vmessage(c("The noise model is:", modelNoise), verbose=verbose)
            vmessage(c("Proportion of non-genetic covariate variance (delta):", 
                       delta), verbose=verbose)
            vmessage(c("Proportion of variance of shared non-genetic covariate",
                       "effects (gamma):", gamma), verbose=verbose)
            vmessage(c("Proportion of non-genetic covariates to have",
                       "a trait-independent effect (pIndependentConfounders"
                       , "):", pIndependentConfounders), 
                     verbose=verbose)
            vmessage(c("Proportion of traits influenced by independent",
                       "non-genetic covariate effects", 
                       "(pTraitIndependentConfounders):", 
                       pTraitIndependentConfounders), verbose=verbose)
            vmessage(c("Proportion of variance of correlated noise effects",
                       "(rho):", rho), verbose=verbose)
            vmessage(c("Proportion of observational noise variance (phi):", 
                       phi), verbose=verbose)
            vmessage(c("Variance of shared observational noise effect (alpha):", 
                       alpha), verbose=verbose)
        }
        vmessage("\n", verbose=verbose)
    }
    
    vmessage(c("The total genetic variance (genVar) is:", genVar), 
             verbose=verbose)
    if ( genVar == 0 ) {
        modelGenetic="noGenetic"
        h2s <- 0
        h2bg <- 0
        vmessage(c("The genetic model is:", modelGenetic), verbose=verbose)
    } else {   
        if (all(c(is.null(h2bg), is.null(h2s)))) {
            stop(paste("Neither genetic variant effect variance (h2s) nor",
                       "infinitesimal genetic effect variance (h2bg) provided",
                       ", at least one is required")
            )
        }
        if (length(c(h2bg, h2s)) == 2 && sum(c(h2bg, h2s)) != 1) {
            stop(paste("Sum of the proportion of the variance of genetic",
                       "effects is not equal to 1; change h2s (genetic variant", 
                       "effect variance) or h2bg (infinitesimal genetic effect", 
                       "variance) such that h2s + h2bg = 1")
            )
        } 
        if (is.null(h2s)) {
            h2s <- 1 - h2bg
        } else {
            h2bg <- 1 - h2s
        }
        if (!is.null(cNrSNP)) {
            if (h2s != 0 && cNrSNP == 0) {
                stop(paste("Proportion of variants of genetic variant effects ",
                           "(h2s) is", h2s, "but number of",
                           "cNrSNP is set to zero"))
            }
        }
        if (theta == 1) {
            pIndependentGenetic=0
            pTraitIndependentGenetic=0
        }
        if (theta == 0) {
            pIndependentGenetic=1
            pTraitIndependentGenetic=1
        }
        if ( h2s == 0) {
            modelGenetic="geneticBgOnly"
            vmessage(c("The genetic model is:", modelGenetic), verbose=verbose)
            vmessage(c("Proportion of variance of infinitesimal genetic", 
                       "effects (h2bg):", h2bg), verbose=verbose)
            vmessage(c("Proportion of variance of shared infinitesimal genetic",
                       "effects (eta):", eta), verbose=verbose)
        } else if ( h2s == 1) {
            modelGenetic="geneticFixedOnly"
            vmessage(c("The genetic model is:", modelGenetic), verbose=verbose)
            vmessage(c("Proportion of variance of genetic variant effects",
                       "(h2s):", h2s), verbose=verbose)
            vmessage(c("Proportion of variance of shared genetic variant",
                       "effects (theta):", theta), verbose=verbose)
            vmessage(c("Proportion of genetic variant effects to have a", 
                       "trait-independent fixed effect", 
                       "(pIndependentGenetic):", pIndependentGenetic), 
                     verbose=verbose)
            vmessage(c("Proportion of traits influenced by independent",
                       "genetic variant effects (pTraitIndependentGenetic):",
                       pTraitIndependentGenetic), 
                     verbose=verbose)
        } else {
            modelGenetic="geneticFixedAndBg"
            vmessage(c("The genetic model is:", modelGenetic), verbose=verbose)
            vmessage(c("Proportion of variance of genetic variant effects",
                       "(h2s):", h2s), verbose=verbose)
            vmessage(c("Proportion of variance of shared genetic variant",
                       "effects (theta):", theta), verbose=verbose)
            vmessage(c("Proportion of genetic variant effects to have a", 
                       "trait-independent fixed effect", 
                       "(pIndependentGenetic):", pIndependentGenetic), 
                     verbose=verbose)
            vmessage(c("Proportion of traits influenced by independent",
                       "genetic variant effects (pTraitIndependentGenetic):",
                       pTraitIndependentGenetic), 
                     verbose=verbose)
            vmessage(c("Proportion of variance of infinitesimal genetic", 
                       "effects (h2bg):", h2bg), verbose=verbose)
            vmessage(c("Proportion of variance of shared infinitesimal genetic",
                       "effects (eta):", eta), verbose=verbose)
        }
        vmessage(c("Proportion of non-linear phenotype transformation",
                   "(proportionNonlinear):", proportionNonlinear), 
                 verbose=verbose)
        vmessage("\n", verbose=verbose)
    }
    return(list(modelGenetic=modelGenetic, modelNoise=modelNoise,
                genVar=genVar, h2s=h2s, h2bg=h2bg, noiseVar=noiseVar, rho=rho, 
                delta=delta, phi=phi, gamma=gamma, theta=theta, eta=eta, 
                alpha=alpha, pcorr=pcorr,
                pTraitIndependentGenetic=pTraitIndependentGenetic,
                pIndependentGenetic=pIndependentGenetic,
                pTraitIndependentConfounders=pTraitIndependentConfounders,
                pIndependentConfounders=pIndependentConfounders,
                proportionNonlinear=proportionNonlinear
    ))
}


#' Run phenotype simulation.
#'
#' runSimulation wraps around setModel, the phenotype component functions 
#' (genFixedEffects, genBgEffects, noiseBgEffects, noiseFixedEffects and 
#' correlatedBgEffects), rescales each component and combines them into the 
#' final phenotype. For details to all parameters, see the respective functions.
#'
#' @param P Number [integer] of phenotypes to simulate. 
#' @param N Number [integer] of samples to simulate.
#' @param genVar Proportion [double] of total genetic variance.
#' @param h2s Proportion [double] of genetic variance of genetic variant effects. 
#' @param h2bg Proportion [double] of genetic variance of infinitesimal genetic 
#' effects; either h2s or h2bg have to be specified and h2s + h2bg = 1.
#' @param theta Proportion [double] of variance of shared genetic variant 
#' effects.
#' @param eta Proportion [double] of variance of shared infinitesimal genetic 
#' effects.
#' @param noiseVar Proportion [double] of total noise variance.
#' @param rho Proportion [double] of noise variance of correlated effects; sum 
#' of rho, delta and phi has to be equal 1.
#' @param delta Proportion [double] of noise variance of non-genetic covariate 
#' effects; sum of rho, delta and phi  has to be equal 1.
#' @param gamma Proportion [double] of variance of shared non-genetic covariate 
#' effects.
#' @param phi Proportion [double] of noise variance of observational noise 
#' effects; sum of rho, delta and phi has to be equal 1.
#' @param alpha Variance [double] of shared observational noise effect.
#' @param tNrSNP Total number [integer] of SNPs to simulate; these SNPs are used
#' for kinship estimation.
#' @param cNrSNP Number [integer] of causal SNPs; used as genetic variant 
#' effects.
#' @param SNPfrequencies Vector of allele frequencies [double] from which to 
#' sample.
#' @param genotypefile Needed when reading external genotypes (into memory), 
#' path/to/genotype file [string] in format specified by \link{format}.
#' @param format Needed when reading external genotypes, specifies 
#' the format of the genotype data; has to be one of plink, oxgen, genome, 
#' bimbam and delim when reading files into memory, or one of oxgen, bimbam or
#' delim if sampling genetic variants from file; for details see
#' \link{readStandardGenotypes} and \link{getCausalSNPs}.
#' @param genoFilePrefix Needed when sampling cuasal SNPs from file, full 
#' path/to/chromosome-wise-genotype-file-ending-before-"chrChromosomeNumber" 
#' (no '~' expansion!) [string]
#' @param genoFileSuffix Needed when sampling causal SNPs from file, 
#' following chromosome number including fileformat (e.g. ".csv") [string]
#' @param genoDelimiter Field separator [string] of genotypefile or genoFile if
#' format == delim. 
#' @param skipFields Number [integer] of fields (columns) in to skip in 
#' genoFilePrefix-genoFileSuffix-file. See details in \link{getCausalSNPs} if
#' format == delim. 
#' @param header [logical] Can be set to indicate if
#' genoFilePrefix-genoFileSuffix file has a header for format == 'delim'. 
#' See details in \link{getCausalSNPs}.
#' @param probabilities [bool]. If set to TRUE, the genotypes in the files 
#' described by genoFilePrefix and genoFileSuffix are provided as triplets of 
#' probablities (p(AA), p(Aa), p(aa)) and are converted into their expected 
#' genotype frequencies by 0*p(AA) + p(Aa) + 2p(aa) via \link{probGen2expGen}.
#' @param chr Numeric vector of chromosomes [integer] to chose NrCausalSNPs 
#' from; only used when external genotype data is sampled i.e. 
#' !is.null(genoFilePrefix) 
#' @param NrSNPsOnChromosome Specifies the number of SNPs [integer] per entry in 
#' chr (see above); has to be the same length as chr. If not provided, lines in 
#' genoFilePrefix-genoFileSuffix file will be counted (which can be slow for 
#' large files).
#' @param NrChrCausal Number [integer] of causal chromosomes to chose 
#' NrCausalSNPs from (as opposed to the actual chromosomes to chose from via chr
#' );  only used when external genotype data is sampled i.e. 
#' !is.null(genoFilePrefix).
#' @param kinshipfile path/to/kinshipfile [string]; if provided, 
#' kinship for simulation of genetic backgound effect will be read from file.
#' @param kinshipHeader [boolean] If TRUE kinship file has header information. 
#' @param kinshipDelimiter Field separator [string] of kinship file. 
#' @param standardise [boolean] If TRUE genotypes will be standardised for 
#' kinship estimation (recommended).
#' @param distBetaGenetic Name [string] of distribution to use to simulate 
#' effect sizes of genetic variants; one of "unif" or "norm".
#' @param mBetaGenetic Mean/midpoint [double] of normal/uniform distribution 
#' for effect sizes of genetic variants.
#' @param sdBetaGenetic Standard deviation/extension from midpoint [double] 
#' of normal/uniform distribution for effect sizes of genetic variants.
#' @param pIndependentGenetic Proportion [double] of genetic variant effects to 
#' have a trait-independent fixed effect.
#' @param pTraitIndependentGenetic Proportion [double] of traits influenced by 
#' independent genetic variant effects.
#' @param keepSameIndependentSNPs [boolean] If set to TRUE, the 
#' independent SNPs effects always influence the same subset of traits.
#' @param pTraitsAffectedGenetics Proportion [double] of traits affected by the 
#' genetic variant effect. For non-integer results of pTraitsAffected*P, the 
#' ceiling of the result is used. Allows to simulate for instance different 
#' levels of pleiotropy.
#' @param NrFixedEffects Number [integer] of different non-genetic covariate 
#' effects to simulate; allows to simulate non-genetic covariate effects from 
#' different distributions or with different parameters.
#' @param NrConfounders Number [integer] of non-genetic covariates; used as 
#' non-genetic covariate effects.
#' @param distConfounders Vector of name(s) [string] of distributions to use to 
#' simulate confounders; one of "unif", "norm", "bin", "cat_norm", "cat_unif".
#' @param mConfounders Vector of mean(s)/midpoint(s) [double] of 
#' normal/uniform distribution for confounders.
#' @param sdConfounders Vector of standard deviation(s)/extension from 
#' midpoint(s) [double] of normal/uniform distribution for confounders.
#' @param catConfounders Vector of confounder categories [factor]; required if 
#' distConfounders "cat_norm" or "cat_unif".
#' @param probConfounders Vector of probability(ies) [double] of binomial 
#' confounders (0/1); required if distConfounders "bin". 
#' @param distBetaConfounders Vector of name(s) [string] of distribution to use 
#' to simulate effect sizes of confounders; one of "unif" or "norm".
#' @param mBetaConfounders Vector of mean(s)/midpoint(s) [double] of 
#' normal/uniform distribution for effect sizes of confounders.
#' @param sdBetaConfounders Vector of standard deviation(s)/extension from 
#' midpoint(s) [double] of normal/uniform distribution for effect sizes of 
#' confounders.
#' @param pIndependentConfounders Vector of proportion(s) [double] of 
#' non-genetic covariate effects to have a trait-independent effect.
#' @param pTraitIndependentConfounders Vector of proportion(s) [double] of 
#' traits influenced by independent non-genetic covariate effects.
#' @param keepSameIndependentConfounders [boolean] If set to TRUE, the 
#' independent confounder effects always influence the same subset of traits.
#' @param pTraitsAffectedConfounders Proportion(s) [double] of traits 
#' affected by the non-genetic covariates. For non-integer results of 
#' pTraitsAffected*P, the ceiling of the result is used.
#' @param meanNoiseBg Mean [double] of the normal distributions for the 
#' simulation observational noise effects.
#' @param sdNoiseBg Standard deviation [double] of the normal distributions for 
#' the simulations of the observational noise effects.
#' @param pcorr Correlation [double] between phenotypes.
#' @param corrmatfile path/to/corrmatfile.csv [string] with comma-separated 
#' [P x P] numeric [double] correlation matrix; if provided,  correlation matrix 
#' for simulation of correlated backgound effect will be read from file; 
#' file should NOT contain an index or header column.
#' @param nonlinear nonlinear transformation method [string]; one exp 
#' (exponential), log (logarithm), poly (polynomial), sqrt (squareroot) or 
#' custom (user-supplied function); if log or exp, base can be specified; if 
#' poly, power can be specified; if custom, a custom function (see for details). 
#' Non-linear transformation is optional, default is NULL ie no transformation
#' (see details).
#' @param logbase [int] base of logarithm for non-linear phenotype 
#' transformation (see details).
#' @param expbase [int] base of exponential function for non-linear phenotype 
#' transformation (see details).
#' @param power [double] power of polynomial function for non-linear phenotype 
#' transformation.
#' @param transformNeg [string] transformation method for negative values in non
#' linear phenotype transformation. One of abs (absolute value) or set0 (set all 
#' negative values to zero). If nonlinear==log and transformNeg==set0, negative
#' values set to 1e-5
#' @param customTransform [function] custom transformation function accepting 
#' a single argument.
#' @param proportionNonlinear [double] proportion of the phenotype to be non-
#' linear (see details)
#' @param sampleID Prefix [string] for naming samples (will be followed by 
#' sample number from 1 to N when constructing sample IDs); only used if 
#' genotypes/kinship are simulated/do not have sample IDs.
#' @param phenoID Prefix [string] for naming traits (will be followed by 
#' phenotypes number from 1 to P when constructing phenotype IDs).
#' @param snpID Prefix [string] for naming SNPs (will be followed by 
#' SNP number from 1 to NrSNP when constructing SNP IDs).
#' @param seed Seed [integer] to initiate random number generation.
#' @param verbose [boolean]; If TRUE, progress info is printed to standard out
#' @return Named list of i) dataframe of proportion of variance 
#' explained for each component (varComponents), 
#' ii) a named list with the final simulated phenotype components 
#' (phenoComponentsFinal), iii) a named list with the intermediate simulated 
#' phenotype components (phenoComponentsIntermediate), iv) a named list of 
#' parameters describing the model setup (setup) and v) a named list of raw 
#' components (rawComponents) used for genetic effect simulation (genotypes 
#' and/or kinship, eigenvalues and eigenvectors of kinship)
#' @details Phenotypes are modeled under a linear additive model where
#' Y = WA + BX + G + C + Phi, with WA the non-genetic covariates, BX the genetic
#' variant effects, G the infinitesimal genetic effects, C the correlated 
#' background effects and the Phi the observational noise. For more information
#' on these components look at the respective function descriptions (see also)
#' Optionally the phenotypes can be non-linearly transformed via:
#' Y_trans = (1-alpha) x Y + alpha x f(Y). Alpha is the proportion of non-
#' linearity of the phenotype and f is a non-linear transformation, and one of
#' exp, log or sqrt. 
#' @export
#' @seealso \link{setModel}, \link{geneticFixedEffects},
#'  \link{geneticBgEffects}, \link{noiseBgEffects}, \link{noiseFixedEffects},
#' \link{correlatedBgEffects} and \link{rescaleVariance}.
#' @examples
#' # simulate phenotype of 100 samples, 10 traits from genetic and noise 
#' # background effects, with variance explained of 0.2 and 0.8 respectively
#' genVar = 0.2
#' simulatedPhenotype <- runSimulation(N=100, P=5, cNrSNP=10,
#' genVar=genVar, h2s=1, phi=1)
runSimulation <- function(N, P,
                          genVar=NULL, h2s=NULL, theta=0.8, h2bg=NULL, eta=0.8, 
                          noiseVar=NULL, rho=NULL, delta=NULL, gamma=0.8, 
                          phi=NULL, alpha=0.8, 
                          tNrSNP=5000, cNrSNP=20, 
                          SNPfrequencies=c(0.1, 0.2, 0.4), 
                          genotypefile=NULL, format='delim',
                          genoFilePrefix=NULL, genoFileSuffix=NULL, 
                          genoDelimiter=",", skipFields=NULL, 
                          header=FALSE,
                          probabilities=FALSE,
                          chr=NULL, NrSNPsOnChromosome=NULL, 
                          NrChrCausal=NULL,
                          kinshipfile=NULL, 
                          kinshipHeader=FALSE, kinshipDelimiter=",", 
                          standardise=TRUE,
                          distBetaGenetic="norm", mBetaGenetic=0, 
                          sdBetaGenetic=1,pTraitsAffectedGenetics=1,
                          pIndependentGenetic=0.4, pTraitIndependentGenetic=0.2,
                          keepSameIndependentSNPs=FALSE,
                          NrFixedEffects=1, NrConfounders=10, 
                          distConfounders="norm", mConfounders=0, 
                          sdConfounders=1,catConfounders=NULL, 
                          probConfounders=NULL, distBetaConfounders="norm", 
                          mBetaConfounders=0, sdBetaConfounders=1,
                          pTraitsAffectedConfounders=1,
                          pIndependentConfounders=0.4, 
                          pTraitIndependentConfounders=0.2,
                          keepSameIndependentConfounders=FALSE,
                          pcorr=0.8, corrmatfile=NULL,
                          meanNoiseBg=0, sdNoiseBg=1, 
                          nonlinear=NULL, logbase=10, expbase=NULL, power=NULL,
                          customTransform=NULL, transformNeg="abs",
                          proportionNonlinear=0,
                          sampleID="ID_", phenoID="Trait_", snpID="SNP_",
                          seed=219453, verbose=FALSE) {

    vmessage(c("Set seed:", seed), verbose=verbose)
    set.seed(seed)
    model <- setModel(genVar=genVar, h2s=h2s, h2bg=h2bg, theta=theta, eta=eta, 
                      noiseVar=noiseVar, rho=rho, delta=delta, gamma=gamma, 
                      phi=phi, alpha=alpha, pcorr=pcorr, 
                      pIndependentConfounders=pIndependentConfounders,  
                      pTraitIndependentConfounders=pTraitIndependentConfounders, 
                      pIndependentGenetic=pIndependentGenetic, 
                      pTraitIndependentGenetic=pTraitIndependentGenetic, 
                      proportionNonlinear=proportionNonlinear,
                      cNrSNP=cNrSNP, NrConfounders=NrConfounders,
                      verbose=verbose)
    id_snps <- NULL
    id_samples <- NULL
    id_phenos <- paste(phenoID, 1:P, sep="")

    ### create simulated phenotypes
    # 1. Simulate genetic terms
    vmessage(c("Simulate genetic effects (genetic model: ", 
               model$modelGenetic,")"), sep="", verbose=verbose)
    if (grepl('Fixed', model$modelGenetic)) {
        if (is.null(genoFilePrefix) && is.null(genotypefile)) {
            if (!grepl('Bg', model$modelGenetic) && tNrSNP != cNrSNP) {
                warning(paste("The genetic model does not contain infinitesimal",
                              "genetic effects but the total number of SNPs to",
                              "simulate (tNrSNP:",
                              tNrSNP, ") is larger than the number of genetic 
                            variant effects SNPs (cNrSNP:", cNrSNP, 
                              "). If genotypes are not needed, consider setting 
                            tNrSNPs=cNrSNPs to speed up computation"))
            }
            genotypes <- simulateGenotypes(N=N, NrSNP=tNrSNP, 
                                           frequencies=SNPfrequencies, 
                                           sampleID=sampleID, 
                                           snpID=snpID, 
                                           verbose=verbose)
            id_samples <- genotypes$id_samples
            id_snps <- genotypes$id_snps
        } else if (! is.null(genotypefile)) {
            genotypes <- readStandardGenotypes(N=N, filename=genotypefile, 
                                               format=format,
                                               verbose=verbose, 
                                               sampleID=sampleID, 
                                               snpID=snpID, 
                                               delimiter=genoDelimiter)
            id_samples <- genotypes$id_samples
            id_snps <- genotypes$id_snps
        } else {
            genotypes <- NULL
        }
        causalSNPs <- getCausalSNPs(N=N, NrCausalSNPs=cNrSNP, chr=chr, 
                                    NrChrCausal=NrChrCausal,
                                    NrSNPsOnChromosome=NrSNPsOnChromosome,
                                    genotypes=genotypes$genotypes,
                                    genoFilePrefix=genoFilePrefix, 
                                    genoFileSuffix=genoFileSuffix, 
                                    format=format,
                                    probabilities=probabilities,
                                    skipFields=skipFields,
                                    header=header,
                                    delimiter=genoDelimiter, 
                                    sampleID=sampleID, 
                                    verbose=verbose)
        if (is.null(id_snps)) {
            id_snps <- colnames(causalSNPs)
        }
        if (is.null(id_samples)) {
            id_samples <- rownames(causalSNPs)
        }
        vmessage("Simulate genetic variant effects", verbose=verbose)
        genFixed <- geneticFixedEffects(X_causal=causalSNPs, N=N, P=P, 
                                        pTraitsAffected=
                                            pTraitsAffectedGenetics,
                                        pIndependentGenetic=
                                            model$pIndependentGenetic, 
                                        pTraitIndependentGenetic=
                                            model$pTraitIndependentGenetic,
                                        keepSameIndependent=
                                            keepSameIndependentSNPs,
                                        distBeta=distBetaGenetic, 
                                        mBeta=mBetaGenetic, 
                                        sdBeta=sdBetaGenetic,
                                        id_phenos=id_phenos,
                                        id_samples=id_samples,
                                        phenoID=phenoID,
                                        verbose=verbose)

        var_genFixed_shared <- model$theta * model$h2s * model$genVar
        var_genFixed_independent <- (1 - model$theta) * model$h2s * model$genVar

        genFixed_shared_rescaled <- rescaleVariance(genFixed$shared, 
                                                    var_genFixed_shared)
        genFixed_independent_rescaled <- 
            rescaleVariance(genFixed$independent, var_genFixed_independent)

        Y_genFixed <- addNonNulls(list(genFixed_shared_rescaled$component, 
                                       genFixed_independent_rescaled$component))
    } else {
        genFixed <- NULL
        genFixed_shared_rescaled <- NULL
        genFixed_independent_rescaled <- NULL
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
                                                   delimiter = genoDelimiter)
            }
            if (is.null(genotypes)) {
                genotypes <- simulateGenotypes(N=N, NrSNP=tNrSNP, 
                                               frequencies=SNPfrequencies, 
                                               sampleID=sampleID, 
                                               verbose=verbose)
            }
            id_samples <- genotypes$id_samples
            id_snps <- genotypes$id_snps
            
            kinship <- getKinship(X=genotypes$genotypes, 
                                  N=N, standardise=standardise, 
                                  id_samples=id_samples,
                                  sampleID=sampleID,
                                  verbose=verbose)
        } else {
            vmessage("Read kinship from file", verbose=verbose)
            kinship <- getKinship(kinshipfile=kinshipfile, 
                                  sep=kinshipDelimiter, 
                                  N=N, header=kinshipHeader, verbose=verbose,
                                  sampleID=sampleID, id_samples=id_samples)
            if(!is.null(genotypes)) {
                if(colnames(kinship) != genotypes$id_samples) {
                    stop(paste("Sample names in kinship file do not match", 
                               "sample names of genotypes", sep=""))
                }
            } else {
                id_samples <- colnames(kinship)
            }
        }

        if (eta == 1) {
            genBgShared <- TRUE
            genBgIndependent <- FALSE
        }
        else if (eta == 0) {
            genBgShared <- FALSE
            genBgIndependent <- TRUE
        } else {
            genBgShared <- TRUE
            genBgIndependent <- TRUE
        }

        vmessage("Simulate infinitesimal genetic effects", verbose=verbose)
        genBg <- geneticBgEffects(N=N, P=P, kinship=kinship, 
                                  shared=genBgShared, 
                                  independent=genBgIndependent,
                                  id_phenos=id_phenos)
        eval_kinship <- genBg$eval_kinship
        evec_kinship <- genBg$evec_kinship
        
        var_genBg_shared <- model$eta * model$h2bg * model$genVar
        var_genBg_independent <- (1 - model$eta) * model$h2bg * model$genVar

        genBg_shared_rescaled <- rescaleVariance(genBg$shared, var_genBg_shared)
        genBg_independent_rescaled <- rescaleVariance(genBg$independent,
                                                      var_genBg_independent)

        Y_genBg <- addNonNulls(list(genBg_shared_rescaled$component, 
                                    genBg_independent_rescaled$component))

        cov_genBg_shared <- genBg$cov_shared
        cov_genBg_shared_rescaled <- cov_genBg_shared * 
            genBg_shared_rescaled$scale_factor^2

        cov_genBg_independent <- genBg$cov_independent
        cov_genBg_independent_rescaled <- cov_genBg_independent * 
            genBg_independent_rescaled$scale_factor^2

        cov_genBg <- cov_genBg_shared_rescaled + cov_genBg_independent_rescaled
    } else {
        genBg <- NULL
        genBg_shared_rescaled <- NULL
        genBg_independent_rescaled <- NULL
        kinship <- NULL
        var_genBg_shared <- 0
        var_genBg_independent <- 0
        Y_genBg <- NULL
        cov_genBg <- NULL
        cov_genBg_shared <- NULL
        cov_genBg_independent <- NULL
        eval_kinship <- NULL
        evec_kinship <- NULL
    }
    # 1. Simulate noise terms
    vmessage(c("Simulate noise terms (noise model: ", model$modelNoise, ")"),
             sep="", verbose=verbose)
    if (grepl('Correlated', model$modelNoise))  {
        corr_mat <- NULL
        vmessage("Simulate correlated background effects", verbose=verbose)
        if (!is.null(corrmatfile)){
            vmessage("Read file with correlation matrix for correlated",  
                     "background effect", verbose=verbose)
            corr_mat <- as.matrix(data.table::fread(corrmatfile, sep=",",
                                                    data.table=FALSE,
                                                    stringsAsFactors=FALSE))
        }
        correlatedBg <- correlatedBgEffects(N=N, P=P, corr_mat=corr_mat, 
                                            pcorr=model$pcorr, 
                                            id_phenos=id_phenos,
                                            id_samples=id_samples,
                                            sampleID=sampleID,
                                            phenoID=phenoID)
        var_noiseCorrelated <- model$rho *  model$noiseVar
        correlatedBg_rescaled <- rescaleVariance(correlatedBg$correlatedBg, 
                                                 var_noiseCorrelated)

        Y_correlatedBg <- correlatedBg_rescaled$component

        cov_correlatedBg <- correlatedBg$cov_correlated * 
            correlatedBg_rescaled$scale_factor^2
    } else {
        correlatedBg <- NULL
        var_noiseCorrelated <- 0
        Y_correlatedBg <- NULL 
        cov_correlatedBg <- NULL
    }
    if (grepl('Bg', model$modelNoise)) {
        if (alpha == 1) {
            noiseBgShared <- TRUE
            noiseBgIndependent <- FALSE
        }
        else if (alpha == 0) {
            noiseBgShared <- FALSE
            noiseBgIndependent <- TRUE
        } else {
            noiseBgShared <- TRUE
            noiseBgIndependent <- TRUE
        }
        vmessage("Simulate observational noise effects", verbose=verbose)
        noiseBg <- noiseBgEffects(N=N, P=P, mean=meanNoiseBg, sd=sdNoiseBg,
                                  shared=noiseBgShared, 
                                  independent=noiseBgIndependent,
                                  id_phenos=id_phenos, id_samples=id_samples,
                                  sampleID=sampleID, phenoID=phenoID)

        var_noiseBg_shared <- model$alpha * model$phi * model$noiseVar
        var_noiseBg_independent <- (1 - model$alpha) * model$phi * model$noiseVar

        noiseBg_shared_rescaled <- rescaleVariance(noiseBg$shared, 
                                                   var_noiseBg_shared)
        noiseBg_independent_rescaled <- rescaleVariance(noiseBg$independent, 
                                                        var_noiseBg_independent)

        Y_noiseBg <- addNonNulls(list(noiseBg_shared_rescaled$component, 
                                      noiseBg_independent_rescaled$component))

        cov_noiseBg_shared <- noiseBg$cov_shared
        cov_noiseBg_shared_rescaled <- cov_noiseBg_shared * 
            noiseBg_shared_rescaled$scale_factor^2

        cov_noiseBg_independent <- noiseBg$cov_independent
        cov_noiseBg_independent_rescaled <- cov_noiseBg_independent * 
            noiseBg_independent_rescaled$scale_factor^2

        cov_noiseBg <- cov_noiseBg_shared_rescaled + 
            cov_noiseBg_independent_rescaled

    } else {
        noiseBg <- NULL
        noiseBg_shared_rescaled <- NULL
        noiseBg_independent_rescaled <- NULL
        var_noiseBg_shared <- 0
        var_noiseBg_independent <- 0
        Y_noiseBg <- NULL
        cov_noiseBg <- NULL
        cov_noiseBg_shared <- NULL
        cov_noiseBg_independent <- NULL
    }
    if (grepl('Fixed', model$modelNoise)) {
        vmessage("Simulate confounder effects", verbose=verbose)
        noiseFixed <- noiseFixedEffects(P=P, N=N, 
                                        NrFixedEffects = NrFixedEffects,
                                        NrConfounders=NrConfounders,
                                        pTraitsAffected=
                                            pTraitsAffectedConfounders,
                                        pIndependentConfounders=
                                            model$pIndependentConfounders,
                                        pTraitIndependentConfounders=
                                            model$pTraitIndependentConfounders,
                                        keepSameIndependent=
                                            keepSameIndependentConfounders,
                                        distConfounders=distConfounders, 
                                        mConfounders=mConfounders, 
                                        sdConfounders=sdConfounders, 
                                        catConfounders=catConfounders, 
                                        probConfounders = probConfounders,
                                        distBeta=distBetaConfounders, 
                                        mBeta=mBetaConfounders, 
                                        sdBeta=sdBetaConfounders,
                                        id_phenos=id_phenos,
                                        id_samples=id_samples,
                                        sampleID=sampleID,
                                        phenoID=phenoID,
                                        verbose=verbose)

        var_noiseFixed_shared <- model$gamma * model$delta * model$noiseVar
        var_noiseFixed_independent <- (1 - model$gamma) * model$delta * 
            model$noiseVar

        noiseFixed_shared_rescaled <- rescaleVariance(noiseFixed$shared, 
                                                      var_noiseFixed_shared)
        noiseFixed_independent_rescaled <- rescaleVariance(
            noiseFixed$independent, 
            var_noiseFixed_independent)

        Y_noiseFixed <- 
            addNonNulls(list(noiseFixed_shared_rescaled$component, 
                            noiseFixed_independent_rescaled$component))
    } else {
        noiseFixed <- NULL
        noiseFixed_shared_rescaled <- NULL
        noiseFixed_independent_rescaled <- NULL
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

    # 4. Transformation
    if (!is.null(nonlinear)) {
        vmessage("Transform phenotypes", verbose=verbose)
        
        Y_transformed <- transformNonlinear(Y, method=nonlinear, 
                                          alpha=proportionNonlinear, 
                                          logbase=logbase, expbase=expbase,
                                          power=power, f=customTransform,
                                          transformNeg=transformNeg)
        Y_transformed <- scale(Y_transformed)
    } else {
        Y_transformed <- NULL
    }
    Y <- scale(Y)

    varComponents <- data.frame(genVar=model$genVar, h2s=model$h2s, 
                                h2bg=model$h2bg, 
                                proportionNonlinear=model$proportionNonlinear,
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
    phenoComponentsFinal <- list(Y=Y, Y_genFixed=Y_genFixed, Y_genBg=Y_genBg, 
                                 Y_noiseFixed=Y_noiseFixed, Y_noiseBg=Y_noiseBg, 
                                 Y_correlatedBg=Y_correlatedBg, 
                                 Y_transformed=Y_transformed,
                                 cov_genBg=cov_genBg, 
                                 cov_noiseBg=cov_noiseBg,
                                 cov_correlatedBg = cov_correlatedBg)
    phenoComponentsIntermediate <- list(
        Y_genFixed_shared=
            genFixed_shared_rescaled$component,
        Y_genFixed_independent=
            genFixed_independent_rescaled$component,
        Y_noiseFixed_shared=
            noiseFixed_shared_rescaled$component,
        Y_noiseFixed_independent=
            noiseFixed_independent_rescaled$component,
        Y_genBg_shared=
            genBg_shared_rescaled$component,
        Y_genBg_independent=
            genBg_independent_rescaled$component,
        Y_noiseBg_shared=
            noiseBg_shared_rescaled$component,
        Y_noiseBg_independent=
            noiseBg_independent_rescaled$component,
        cov_genBg_shared=cov_genBg_shared, 
        cov_genBg_independent=cov_genBg_independent, 
        cov_noiseBg_shared=cov_noiseBg_shared,
        cov_noiseBg_independent=cov_noiseBg_independent, 
        genFixed=genFixed, 
        noiseFixed=noiseFixed)

    setup <- list(P=P, N=N, NrCausalSNPs=cNrSNP, 
                  modelGenetic=model$modelGenetic, 
                  modelNoise=model$modelNoise, 
                  id_samples=id_samples, id_phenos=id_phenos, id_snps=id_snps)
    rawComponents <- list(kinship=kinship, genotypes=genotypes,
                          eval_kinship=eval_kinship, evec_kinship=evec_kinship)

    return(list(varComponents=varComponents, 
                phenoComponentsFinal=phenoComponentsFinal, 
                phenoComponentsIntermediate=phenoComponentsIntermediate, 
                setup=setup, rawComponents=rawComponents))
}



