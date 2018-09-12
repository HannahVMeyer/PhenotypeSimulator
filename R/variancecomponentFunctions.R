#' Simulate genetic variant effects.
#'
#' geneticFixedEffects takes genetic variants which should be added as genetic 
#' variant effects to the phenotype. These variants can have the same effects 
#' across all traits (shared) or can be independent across traits (independent); 
#' in addition, only a certain proportion of traits can be affected by the 
#' genetic variants.
#' 
#' @param X_causal [N x NrCausalSNPs] Matrix of [NrCausalSNPs] SNPs from [N] 
#' samples.
#' @param N Number [integer] of samples to simulate; has to be provided as a 
#' dimnesionality check for X_causal and downstream analyses; nrow(X_causal) has 
#' to be equal to N.
#' @param P Number [integer] of phenotypes to simulate.
#' @param phenoID Prefix [string] for naming traits.
#' @param id_samples Vector of [NrSamples] sample IDs [string]; if not provided
#' colnames(X_causal) used.
#' @param id_phenos Vector of [NrTraits] phenotype IDs [string]; if not provided
#' constructed by paste(phenoID, 1:P, sep="").
#' @param pTraitsAffected Proportion [double] of traits affected by the genetic
#' effect. For non-integer results of pTraitsAffected*P, the ceiling of the 
#' result is used. Allows to simulate for instance different levels of 
#' pleiotropy.
#' @param pIndependentGenetic Proportion [double] of genetic effects (SNPs) to
#' have an independent fixed effect.
#' @param pTraitIndependentGenetic Proportion [double] of traits influenced by 
#' independent fixed genetic effects.
#' @param keepSameIndependent [boolean] If set to TRUE, the independent genetic 
#' effects always influence the same subset of traits.
#' @param distBeta Vector of name(s) [string] of distribution to use to simulate 
#' effect sizes of SNPs; one of "unif" or "norm".
#' @param mBeta Vector of mean/midpoint(s) [double] of normal/uniform 
#' distribution for effect sizes of SNPs.
#' @param sdBeta Vector of standard deviation/distance from midpoint [double] 
#' of normal/uniform distribution for effect sizes of SNPs.
#' @param verbose [boolean] If TRUE, progress info is printed to standard out
#' @return Named list of shared fixed genetic effects (shared: [N x P] matrix), 
#' independent fixed genetic effects (independent: [N x P] matrix), 
#' the causal SNPs labeled as shared or independent effect 
#' (cov: [NrCausalSNPs x N] matrix) and the simulated effect sizes of the causal 
#' SNPs (cov_effect: [P x NrCausalSNPs] dataframe).
#' @export
#' @examples
#' genotypes <- simulateGenotypes(N=100, NrSNP=20, verbose=FALSE)
#' causalSNPs <- getCausalSNPs(N=100, genotypes=genotypes$genotypes)
#' geneticFixed <- geneticFixedEffects(N=100, X_causal=causalSNPs, 
#' P=10)
geneticFixedEffects <- function(X_causal, P, N, phenoID="Trait_",
                                id_samples = rownames(X_causal),
                                id_phenos = NULL,
                                pTraitsAffected=1,
                                pIndependentGenetic=0.4, 
                                pTraitIndependentGenetic=0.2, 
                                keepSameIndependent=FALSE,
                                distBeta="norm", mBeta=0, sdBeta=1, 
                                verbose=FALSE) {
    numbers <- list(P=P, N=N, mBeta=mBeta, sdBeta=sdBeta,
                    pIndependentGenetic=pIndependentGenetic, 
                    pTraitIndependentGenetic=pTraitIndependentGenetic, 
                    pTraitsAffected=pTraitsAffected)
    positives <- list(P=P, N=N, sdBeta=sdBeta)
    proportions <- list(pIndependentGenetic=pIndependentGenetic, 
                        pTraitIndependentGenetic=pTraitIndependentGenetic, 
                        pTraitsAffected=pTraitsAffected)
    testNumerics(numbers=numbers, positives=positives, proportions=proportions)
    if (!is.numeric(X_causal)) {
        stop("Genetic variant matrix to simulate genetic variant effects from is
              not numeric. Check your genotype simulation or the parameters
             specified for reading the genotypes from file. Did you provide the
             correct format information?")
    }
    if (nrow(X_causal) != N){
        stop("Number of samples in SNP matrix (", nrow(X_causal), ") is 
             different from number of samples to be simulated")
    }
    if (length(id_samples) !=  nrow(X_causal)) {
        stop("Length of id_samples (", length(id_samples), ") is different ",
             "from number of samples in X_causal (", nrow(X_causal), "). Does ",
             "your X_causal have rownames (default to retrieve id_samples if ",
             "not provided)?")
    }
    if (!(is.character(phenoID) && length(phenoID) == 1)) {
        stop("phenoID has to be of length 1 and of type character")
    }
    if (!is.null(id_phenos) && length(id_phenos) !=  P) {
        stop("Length of id_phenos (", length(id_phenos), ") is different ",
             "from P (", P, ")")
    }
    if (is.null(id_phenos)) id_phenos <- paste(phenoID, 1:P, sep="")
    NrCausalSNPs <- ncol(X_causal)
    traitsAffected <- ceiling(P*pTraitsAffected)
    
    NrIndependentSNPs <- round(pIndependentGenetic * NrCausalSNPs)
    NrSharedSNPs <- NrCausalSNPs - NrIndependentSNPs
    
    vmessage(c("Out of", P, "total phenotypes,", traitsAffected, "traits will",
               "be affected by genetic variant effects"), verbose=verbose)
    
    if (traitsAffected == 1 && NrIndependentSNPs != 0) {
        vmessage(c("The total number of traits affected by genetic variant",
                 "effects is 1, so pTraitIndependentGenetic",
                 " will automatically be set to 1."), verbose=verbose)
        pTraitIndependentGenetic <- 1
    }
    
    if (NrIndependentSNPs != 0) {
        vmessage(c("Out of these affected traits (", traitsAffected, "), ", 
                   ceiling(pTraitIndependentGenetic * traitsAffected), 
                   " trait(s)",
                   " will have independent genetic variant effects"), sep="",
                 verbose=verbose)
    }
    
    Gshared <- NULL
    Gindependent <- NULL
    
    if (NrSharedSNPs != 0) {
        if (NrIndependentSNPs != 0) {
            shared <- sample(c(rep(TRUE, NrSharedSNPs), 
                               rep(FALSE, NrIndependentSNPs)), 
                             replace=FALSE)
            X_shared <-  X_causal[,shared]
            snpIDshared <- colnames(X_causal)[shared]
        } else {
            X_shared <- X_causal
            snpIDshared <- colnames(X_causal)
        }
        
        if (distBeta == "unif") {
            betaX_exp <- rexp(NrSharedSNPs) %*% t(rexp(traitsAffected))
            ## transfrom to uniform with mean=m, with m*x=0.5
            betaX_unif <- exp(-betaX_exp/mean(betaX_exp))
            multiplicative = 2*sdBeta/(max(betaX_unif)-min(betaX_unif))
            additive = mBeta + sdBeta - multiplicative * max(betaX_unif)
            betaX_shared  = multiplicative * betaX_unif + additive
        }
        if (distBeta == "norm") {
            betaX_shared <- simulateDist(NrSharedSNPs,  dist=distBeta, 
                                         m=mBeta, std=sdBeta) %*% 
                t(abs(simulateDist(traitsAffected, dist=distBeta, m=0, 
                                   std=1)))
        }
        
        if (P != traitsAffected) {
            betaX_shared <- cbind(betaX_shared, 
                                  matrix(0, ncol=P-traitsAffected, 
                                         nrow=NrSharedSNPs))
        }
        rownames(betaX_shared) <- paste("sharedEffect", 
                                        1:nrow(betaX_shared), sep="")
        
        cov <- data.frame(X_shared)
        colnames(cov) <- snpIDshared
        rownames(cov) <- id_samples
        
        cov_effect <- data.frame(t(betaX_shared))
        colnames(cov_effect) <- paste(rownames(betaX_shared), "_", 
                                      colnames(cov), sep="")
        rownames(cov_effect) <- id_phenos
        
        Gshared = X_shared %*% betaX_shared
        colnames(Gshared) <- id_phenos
        rownames(Gshared) <- id_samples
    }
    
    if (NrIndependentSNPs != 0) {
        if (NrSharedSNPs != 0) {
            independent <- !shared
            X_independent <- X_causal[,independent]
            snpIDindependent <- colnames(X_causal)[independent]
        } else {
            X_independent <- X_causal
            snpIDindependent <- colnames(X_causal)
        }
        
        betaX_independent <- matrix(simulateDist(traitsAffected * 
                                                     NrIndependentSNPs, 
                                                 dist=distBeta,
                                                 m=mBeta, std=sdBeta), 
                                    ncol=traitsAffected)
        
        TraitIndependentGenetic <- ceiling(pTraitIndependentGenetic * 
                                               traitsAffected)
        if (keepSameIndependent) {
            p_nongenetic <- sample(
                c(rep(FALSE, TraitIndependentGenetic), 
                  rep(TRUE, (traitsAffected - TraitIndependentGenetic))), 
                replace=FALSE)
            p_nongenetic <- matrix(rep(p_nongenetic, NrIndependentSNPs), 
                                   NrIndependentSNPs, byrow = TRUE)
        } else {
            p_nongenetic <- t(sapply(1:ncol(X_independent), function(x) {
                sample(c(rep(FALSE, TraitIndependentGenetic), 
                         rep(TRUE, (traitsAffected - TraitIndependentGenetic))), 
                       replace=FALSE)
            }))
        }
        
        betaX_independent[p_nongenetic] <- 0
        
        if (P != traitsAffected) {
            betaX_independent <- cbind(betaX_independent, 
                                       matrix(0, ncol=P-traitsAffected, 
                                              nrow=NrIndependentSNPs))
        }
        rownames(betaX_independent) <- paste("independentEffect", 
                                             1:nrow(betaX_independent), sep="") 
        cov <- data.frame(X_independent)
        colnames(cov) <- snpIDindependent
        rownames(cov) <- id_samples
        
        cov_effect <- data.frame(t(betaX_independent))
        colnames(cov_effect) <- paste(rownames(betaX_independent), "_", 
                                      colnames(cov), sep="")
        rownames(cov_effect) <- id_phenos
        
        Gindependent = X_independent %*% betaX_independent
        colnames(Gindependent) <- id_phenos
        rownames(Gindependent) <- id_samples
    }
    if (NrSharedSNPs != 0 && NrIndependentSNPs != 0) {
        cov = cbind(X_shared, X_independent)
        colnames(cov) <- c( snpIDshared, snpIDindependent)
        cov_effect = data.frame(betaX_shared=t(betaX_shared), 
                                betaX_independent=t(betaX_independent))
        colnames(cov_effect) <- paste(colnames(cov_effect), "_", 
                                      colnames(cov), sep="")
        rownames(cov_effect) <- id_phenos
    }
    
    return(list(shared=Gshared, 
                independent=Gindependent, 
                cov=as.matrix(cov),
                cov_effect=as.matrix(cov_effect)))
}

#' Simulate noise fixed effects.
#'
#' noiseFixedEffects simulates a number of non-genetic covariate effects 
#' (confounders). Confounders can have effects across all traits (shared) 
#' or to a number of traits only (independent); in addition, only a certain 
#' proportion of traits can be affected by the confounders.
#' Confounders can be simulated as categorical variables or following a binomial 
#' , uniform or normal distribution. Effect sizes for the noise effects can be 
#' simulated from a uniform or normal distribution. Multiple confounder sets 
#' drawn from different distributions/different parameters of the same 
#' distribution can be simulated by specifying NrFixedEffects and supplying the 
#' respective distribution parameters. 
#'
#' @param N Number [integer] of samples to simulate.
#' @param P Number [integer] of phenotypes to simulate.
#' @param sampleID Prefix [string] for naming samples.
#' @param phenoID Prefix [string] for naming traits.
#' @param id_samples Vector of [NrSamples] sample IDs [string]; if not provided
#' constructed by paste(sampleID, 1:N, sep="").
#' @param id_phenos Vector of [NrTraits] phenotype IDs [string]; if not provided
#' constructed by paste(phenoID, 1:P, sep="").
#' @param pTraitsAffected Vector of proportion(s) [double] of traits affected by 
#' the confounders. For non-integer results of pTraitsAffected*P, the ceiling of 
#' the result is used.
#' @param NrFixedEffects Number [integer] of different confounder effects to 
#' simulate; allows to simulate fixed effects from different distributions or 
#' with different parameters; if only one type of confounder distribution is 
#' wanted, set NrFixedEffects=1 and choose the number of confounders with eg 
#' NrConfounders=10.
#' @param NrConfounders Vector of number(s) [integer] of confounders from a 
#' specified distribution to simulate.
#' @param pIndependentConfounders Vector of proportion(s) [double] of 
#' confounders to have a trait-independent effect.
#' @param pTraitIndependentConfounders Vector of proportion(s) [double] of 
#' traits influenced by independent confounder effects.
#' @param keepSameIndependent [boolean] If set to TRUE, the independent genetic 
#' effects always influence the same subset of traits.
#' @param distConfounders Vector of name(s) [string] of distribution to use to 
#' simulate confounders; one of "unif", "norm", "bin", "cat_norm", "cat_unif".
#' @param distBeta Vector of name(s) [string] of distribution to use to simulate 
#' effect sizes of confounders; one of "unif" or "norm".
#' @param mConfounders Vector of mean/midpoint(s) [double] of normal/uniform 
#' distribution for confounders.
#' @param sdConfounders Vector of standard deviation(s)/distance from 
#' midpoint(s) [double] of normal/uniform distribution for confounders.
#' @param catConfounders Vector of number(s) of confounder categories [integer]; 
#' required if distConfounders "cat_norm" or "cat_unif".
#' @param probConfounders Vector of probability(s) [double] of binomial 
#' confounders (0/1); required if distConfounders "bin".
#' @param mBeta Vector of mean/midpoint [double] of normal/uniform distribution 
#' for effect sizes of confounders.
#' @param sdBeta Vector of standard deviation/distance from midpoint [double] 
#' of normal/uniform distribution for effect sizes of confounders.
#' @param verbose [boolean] If TRUE, progress info is printed to standard out
#' @return Named list of shared confounder effects (shared: [N x P] matrix), 
#' independent confoudner effects (independent: [N x P] matrix), 
#' the confounders labeled as shared or independent effect 
#' (cov: [NrConfounders x N] matrix) and the simulated effect sizes of the 
#' confounders (cov_effect: [P x NrConfounders] dataframe).
#' @seealso \code{\link{simulateDist}} 
#' @export
#' @examples
#' # fixed noise effect with default setting
#' noiseFE <- noiseFixedEffects(P=5, N=20)
#' 
#' # 1 categorical fixed noise effect with uniform distribution of the 
#' # categories
#' noiseFE_catUnif <- noiseFixedEffects(P=10, N=20, NrConfounders=1, 
#' distConfounders="cat_unif", catConfounders=3)
#' 
#' # 10 fixed noise effect with uniform distribution between 1 and 5 (3 +/- 2) 
#' # categories
#' noiseFE_uniformConfounders_normBetas <- noiseFixedEffects(P=10, N=20, 
#' NrConfounders=10, distConfounders="unif", mConfounders=3, sdConfounders=2, 
#' distBeta="norm",  sdBeta=2)
#' 
#'  # 4 fixed noise effect with binomial distribution with p=0.2 
#' noiseFE_binomialConfounders_uniformBetas <- noiseFixedEffects(P=10, N=20, 
#' NrConfounders=4, distConfounders="bin", probConfounders=0.2, distBeta="norm", 
#' sdBeta=2)
#' 
#'  # 2 fixed noise effect with 1 binomial confounders and 1 normally 
#'  # distributed confounder; the latter only affects 2 traits 
#'  noiseFE_binomialandNormalConfounders <- noiseFixedEffects(P=10, N=20, 
#'  NrFixedEffects=2, pTraitsAffected =c (1,0.2), NrConfounders=c(2,2), 
#'  distConfounders=c("bin", "norm"),  probConfounders=0.2)
noiseFixedEffects <- function(N, P, NrConfounders=10, sampleID="ID_",
                              phenoID="Trait_", 
                              id_samples=NULL,
                              id_phenos=NULL,
                              pTraitsAffected=1,
                              NrFixedEffects=1, 
                              pIndependentConfounders=0.4, 
                              pTraitIndependentConfounders=0.2, 
                              keepSameIndependent=FALSE,
                              distConfounders="norm", mConfounders=0, 
                              sdConfounders=1, catConfounders=NULL, 
                              probConfounders=NULL, distBeta="norm", mBeta=0, 
                              sdBeta=1,
                              verbose=FALSE) {
    numbers <- list(N=N, P=P, NrFixedEffects=NrFixedEffects, 
                    NrConfounders=NrConfounders)
    positives <- list(P=P, N=N, NrFixedEffects=NrFixedEffects, 
                      NrConfounders=NrConfounders)
    integers <- list(P=P, N=N, NrFixedEffects=NrFixedEffects, 
                     NrConfounders=NrConfounders)
    testNumerics(numbers=numbers, positives=positives, integers=integers)
    if (!is.null(id_samples) && length(id_samples) !=  N) {
        stop("Length of id_samples (", length(id_samples), ") is different ",
             "from N (", N, ")")
    }
    if (!is.null(id_phenos) && length(id_phenos) !=  P) {
        stop("Length of id_phenos (", length(id_phenos), ") is different ",
             "from P (", P, ")")
    }
    if (!(is.character(phenoID) && length(phenoID) == 1)) {
        stop("phenoID has to be of length 1 and of type character")
    }
    if (!(is.character(sampleID) && length(sampleID) == 1)) {
        stop("sampleID has to be of length 1 and of type character")
    }
    if (is.null(id_samples)) id_samples <- paste(sampleID, 1:N, sep="")
    if (is.null(id_phenos)) id_phenos <- paste(phenoID, 1:P, sep="")
    
    
    oneFixedEffectComponent <- function(N, P, NrOfEffect, NrConfounders, 
                                        pTraitsAffected,
                                        pIndependentConfounders, 
                                        pTraitIndependentConfounders,
                                        keepSameIndependent,
                                        distConfounders, mConfounders, 
                                        sdConfounders, catConfounders, 
                                        probConfounders, distBeta, mBeta, 
                                        sdBeta) {
        numbers <- list(mConfounders=mConfounders,
                        sdConfounders=sdConfounders, mBeta=mBeta, 
                        sdBeta=sdBeta,
                        pIndependentConfounders=pIndependentConfounders, 
                        pTraitIndependentConfounders=
                            pTraitIndependentConfounders, 
                        probConfounders=probConfounders,
                        pTraitsAffected=pTraitsAffected)
        positives <- list(sdBeta=sdBeta, sdConfounders=sdConfounders)
        proportions <- list(pIndependentConfounders=pIndependentConfounders, 
                            pTraitIndependentConfounders=
                                pTraitIndependentConfounders, 
                            probConfounders=probConfounders,
                            pTraitsAffected=pTraitsAffected)
        testNumerics(numbers=numbers, positives=positives, 
                     proportions=proportions)
        
        if (is.null(catConfounders) && grepl("cat", distConfounders)) {
            stop("Confounder distribution set to ", distConfounders, " but ",
                 "no categories provided")
        }
        if (is.null(probConfounders) && grepl("bin", distConfounders)) {
            stop("Confounder distribution set to ", distConfounders, " but ",
                 "no probabilities provided")
        }
        
        traitsAffected <- ceiling(P*pTraitsAffected)
        
        if(is.null(NrOfEffect)) {
            vmessage(c("Out of", P, "total phenotypes,", traitsAffected, 
                       "trait(s) will",
                       "be affected by the covariate effect"), 
                     verbose=verbose)
        } else {
            vmessage(c("Out of", P, "total phenotypes,", traitsAffected, 
                       "trait(s) will",
                       "be affected by the ", NrOfEffect ," covariate effect"), 
                     verbose=verbose)
        }
        
        NrIndependentConfounders <- round(pIndependentConfounders * 
                                              NrConfounders)
        NrSharedConfounders <- NrConfounders - NrIndependentConfounders
        
        if (traitsAffected == 1 && NrIndependentConfounders != 0) {
            vmessage(c("The total number of traits affected by covariate ",
                     "effects is 1, so pTraitIndependentConfounders",
                     " will automatically be set to 1."), verbose=verbose)
            pTraitIndependentConfounders <- 1
        }
        
        if (NrIndependentConfounders != 0) {
            vmessage(c("Out of these affected traits (", traitsAffected, "), ", 
                       ceiling(pTraitIndependentConfounders * traitsAffected), 
                       " trait(s)",
                       " will have independent covariate effects"), sep="",
                     verbose=verbose)
        }
        
        Cshared <- NULL
        Cindependent <- NULL
        
        if (NrSharedConfounders != 0) {
            shared <- matrix(simulateDist(N * NrSharedConfounders, 
                                          dist=distConfounders, m=mConfounders, 
                                          std=sdConfounders, 
                                          categories=catConfounders, 
                                          prob=probConfounders), 
                             ncol=NrSharedConfounders)
            colnames(shared) <- paste("sharedConfounder", NrOfEffect, "_", 
                                      distConfounders, 
                                      seq(1, NrSharedConfounders, 1), sep="")
            
            if (distBeta == "unif") {
                beta_exp <- rexp(NrSharedConfounders) %*% t(rexp(traitsAffected))
                ## transfrom to uniform with mean=m, with m*x=0.5
                beta_unif <- exp(-beta_exp/mean(beta_exp))
                multiplicative = 2*sdBeta/(max(beta_unif)-min(beta_unif))
                additive = mBeta + sdBeta - multiplicative * max(beta_unif)
                beta_shared  = multiplicative * beta_unif + additive
            }
            if (distBeta == "norm") {
                beta_shared <- simulateDist(NrSharedConfounders,  dist=distBeta, 
                                            m=mBeta, std=sdBeta) %*% 
                    t(abs(simulateDist(traitsAffected, dist=distBeta, m=0, 
                                       std=1)))
            }
            rownames(beta_shared) <- paste("sharedConfounder", NrOfEffect, "_", 
                                           distConfounders,
                                           "_Beta_", distBeta, 
                                           seq(1, NrSharedConfounders, 1), sep="")
            
            if (P != traitsAffected) {
                beta_shared <- cbind(beta_shared, 
                                     matrix(0, ncol=P-traitsAffected, 
                                            nrow=NrSharedConfounders))
            }
            
            Cshared <- shared %*% beta_shared
            colnames(Cshared) <- id_phenos
            rownames(Cshared) <- id_samples
            
            cov <- data.frame(shared)
            rownames(cov) <- id_samples
            cov_effect <- data.frame(t(beta_shared))
            rownames(cov_effect) <- id_phenos
        } 
        if (NrIndependentConfounders != 0) {
            independent <- matrix(simulateDist(N * NrIndependentConfounders, 
                                               dist=distConfounders, 
                                               m=mConfounders, 
                                               std=sdConfounders, 
                                               categories=catConfounders, 
                                               prob=probConfounders), 
                                  ncol=NrIndependentConfounders)
            colnames(independent) <- paste("independentConfounder", NrOfEffect, 
                                           "_", distConfounders, 
                                           seq(1, NrIndependentConfounders, 1), 
                                           sep="")
            beta_independent <- matrix(simulateDist(traitsAffected * 
                                                        NrIndependentConfounders, 
                                                    dist=distBeta, m=mBeta, 
                                                    std=sdBeta), 
                                       ncol=traitsAffected)
            rownames(beta_independent) <- paste("independentConfounder", 
                                                NrOfEffect, "_", 
                                                distConfounders,
                                                "_Beta_", distBeta, 
                                                seq(1, NrIndependentConfounders, 
                                                    1), 
                                                sep="")
            
            TraitIndependentConfounders <- ceiling(pTraitIndependentConfounders* 
                                                       traitsAffected)
            if (keepSameIndependent) {
                p_nonconfounders <- sample(
                    c(rep(FALSE, TraitIndependentConfounders), 
                      rep(TRUE, (traitsAffected - TraitIndependentConfounders))), 
                    replace=FALSE)
                p_nonconfounders <- matrix(rep(p_nonconfounders, 
                                               NrIndependentConfounders), 
                                           NrIndependentConfounders, 
                                           byrow = TRUE)
            } else {
                p_nonconfounders <- t(sapply(1:ncol(independent), function(x) {
                    sample(c(rep(FALSE, TraitIndependentConfounders), 
                             rep(TRUE, 
                                 (traitsAffected - TraitIndependentConfounders))
                    ), 
                    replace=FALSE)
                }))
            }
            
            beta_independent[p_nonconfounders] <- 0
            
            if (P != traitsAffected) {
                beta_independent <- cbind(beta_independent, 
                                          matrix(0, ncol=P-traitsAffected, 
                                                 nrow=NrIndependentConfounders))
            }
            
            
            Cindependent <- independent %*% beta_independent
            colnames(Cindependent) <- id_phenos
            rownames(Cindependent) <- id_samples
            
            cov <- data.frame(independent)
            rownames(cov) <- id_samples
            cov_effect <- data.frame(t(beta_independent))
            rownames(cov_effect) <- id_phenos
        }
        if (NrSharedConfounders != 0 && NrIndependentConfounders != 0) {
            cov <- data.frame(shared, independent)
            rownames(cov) <- id_samples
            cov_effect <- data.frame(t(beta_shared), t(beta_independent))
            rownames(cov_effect) <- id_phenos
        }
        return(list(shared=Cshared, independent=Cindependent, cov=cov, 
                    cov_effect=cov_effect))
    }
    
    param_length <- c(pTraitsAffected=length(pTraitsAffected), 
                      pIndependentConfounders=
                          length(pIndependentConfounders),
                      pTraitIndependentConfounders=
                          length(pTraitIndependentConfounders),
                      keepSameIndependent=length(keepSameIndependent),
                      distConfounders=length(distConfounders),
                      mConfounders=length(mConfounders),
                      sdConfounders=length(sdConfounders),
                      probConfounders=length(probConfounders),
                      catConfounders=length(catConfounders),
                      distBeta=length(distBeta),
                      mBeta=length(mBeta),
                      sdBeta=length(sdBeta))
    
    common <- c(pTraitsAffected=length(pTraitsAffected), 
                pIndependentConfounders=
                    length(pIndependentConfounders),
                pTraitIndependentConfounders=
                    length(pTraitIndependentConfounders),
                keepSameIndependent=length(keepSameIndependent),
                distConfounders=length(distConfounders),
                distBeta=length(distBeta),
                mBeta=length(mBeta),
                sdBeta=length(sdBeta))
    
    norm_unif <- c(mConfounders=length(mConfounders),
                   sdConfounders=length(sdConfounders))
    norm_unif_dist <- length(which(grepl("^unif", distConfounders))) + 
        length(which(grepl("^norm",distConfounders)))
    
    bin <- c(probConfounders=length(probConfounders))
    bin_dist <- length(which(grepl("bin", distConfounders)))
    
    categ <- c(catConfounders=length(catConfounders))
    categ_dist <- length(which(grepl("cat", distConfounders)))
    
    if (NrFixedEffects != 1 && all(param_length <= 1)) {
        stop("NrFixedEffects specified to greater than 1 (i.e. more than one ",
             "type of fixed effects), but not enough parameters specified for ",
             "multiple fixed effects")
    }
    check_common <- sapply(seq_along(common), function(x){
        if (common[x]  > 1 && common[x] != NrFixedEffects) {
            stop("Length of ", names(common)[x], " (", common[x], 
                 ") doesn't match NrFixedEffects (", NrFixedEffects, ")")
        }
    })
    check_nu <- sapply(seq_along(norm_unif), function(x){
        if (norm_unif[x]  > 1 && norm_unif[x] != norm_unif_dist) {
            stop("Length of ", names(norm_unif)[x], " (", norm_unif[x], 
                 ") doesn't match unif/norm distConfounder length (", 
                 norm_unif_dist, ")")
        }
    })
    if (bin > 1 && bin != bin_dist) {
        stop("Length of ", names(bin), " (", bin, 
             ") doesn't match bin distConfounder length (", 
             bin_dist, ")")
    }   
    if (categ > 1 && categ != categ_dist) {
        stop("Length of ", names(categ), " (", categ, 
             ") doesn't match cat_unif/cat_norm distConfounder length (", 
             categ_dist, ")")
    } 
    
    if (NrFixedEffects != 1)  { 
        un <- intersect(which(!grepl("cat", distConfounders)),
                        which(!grepl("bin", distConfounders)))
        if (length(un) != 0) {
            tmp_m_un <- rep(9, length(distConfounders))
            tmp_sd_un <- rep(9, length(distConfounders))
            tmp_m_un[un] <- mConfounders
            tmp_sd_un[un] <- sdConfounders
            mConfounders <- tmp_m_un
            sdConfounders <- tmp_sd_un
        }
        if (!is.null(catConfounders)) {
            tmp_cat <- rep(9, length(distConfounders))
            tmp_cat[which(grepl("cat", distConfounders))] <- catConfounders
            catConfounders <- tmp_cat
        }
        if (!is.null(probConfounders)) {
            tmp_bin <- rep(0, length(distConfounders))
            tmp_bin[which(grepl("bin", distConfounders))] <- probConfounders
            probConfounders <- tmp_bin
        }
        
        if (length(pTraitsAffected) == 1) {
            pTraitsAffected <- rep(pTraitsAffected, 
                                   NrFixedEffects)
        }
        if (length(pIndependentConfounders) == 1) {
            pIndependentConfounders <- rep(pIndependentConfounders, 
                                           NrFixedEffects)
        }
        if (length(pTraitIndependentConfounders) == 1) {
            pTraitIndependentConfounders <- rep(pTraitIndependentConfounders, 
                                                NrFixedEffects)
        }
        if (length(keepSameIndependent) == 1) {
            keepSameIndependent <- rep(keepSameIndependent, NrFixedEffects)
        }
        if (length(NrConfounders) == 1) {
            NrConfounders <- rep(NrConfounders, NrFixedEffects)
        }
        if (length(distConfounders) == 1) {
            distConfounders <- rep(distConfounders, NrFixedEffects)
        }
        if (length(mConfounders) == 1) {
            mConfounders <- rep(mConfounders, NrFixedEffects)
        }
        if (length(sdConfounders) == 1) {
            sdConfounders <- rep(sdConfounders, NrFixedEffects)
        }
        if (length(catConfounders) <= 1) {
            catConfounders <- rep(catConfounders, NrFixedEffects)
        }
        if (length(probConfounders) <= 1) {
            probConfounders <- rep(probConfounders, NrFixedEffects)
        }
        if (length(distBeta) == 1) {
            distBeta <- rep(distBeta, NrFixedEffects)
        }
        if (length(mBeta) == 1) {
            mBeta <- rep(mBeta, NrFixedEffects)
        }
        if (length(sdBeta) == 1) {
            sdBeta <- rep(sdBeta, NrFixedEffects)
        }
        
        tmp <- lapply(1:NrFixedEffects, function(x) {
            oneFixedEffectComponent(N=N, P=P, NrOfEffect=x, 
                                    NrConfounders=NrConfounders[x], 
                                    pTraitsAffected=pTraitsAffected[x],
                                    pIndependentConfounders=
                                        pIndependentConfounders[x], 
                                    pTraitIndependentConfounders=
                                        pTraitIndependentConfounders[x], 
                                    keepSameIndependent=
                                        keepSameIndependent[x],
                                    distConfounders=distConfounders[x], 
                                    mConfounders=mConfounders[x], 
                                    sdConfounders=sdConfounders[x], 
                                    catConfounders=catConfounders[x], 
                                    probConfounders=probConfounders[x], 
                                    distBeta=distBeta[x], mBeta=mBeta[x], 
                                    sdBeta=sdBeta[x])
        })
        shared <- addNonNulls(lapply(tmp, function(x) x$shared))
        independent <- addNonNulls(lapply(tmp, function(x) x$independent))
        cov <- do.call(cbind, lapply(tmp, function(x) x$cov))
        cov_effect <- do.call(cbind, lapply(tmp, function(x) x$cov_effect))
        return(list(shared=shared, independent=independent, cov=cov, 
                    cov_effect=cov_effect))
    } else {
        oneFixedEffectComponent(N=N, P=P, NrOfEffect=NULL, 
                                NrConfounders=NrConfounders, 
                                pTraitsAffected=pTraitsAffected,
                                pIndependentConfounders=
                                    pIndependentConfounders, 
                                pTraitIndependentConfounders=
                                    pTraitIndependentConfounders,
                                keepSameIndependent=keepSameIndependent,
                                distConfounders=distConfounders, 
                                mConfounders=mConfounders, 
                                sdConfounders=sdConfounders, 
                                catConfounders=catConfounders, 
                                probConfounders=probConfounders, 
                                distBeta=distBeta, mBeta=mBeta, sdBeta=sdBeta) 
    } 
}


#' Simulate infinitesimal genetic effects (reflecting sample kinship).
#'
#' geneticBgEffects simulates an infinitesimal genetic effects with a proportion 
#' of the effect shared across samples and a proportion independent across 
#' samples; they are based on the kinship estimates of the (simulated) samples.
#'
#' @param P Number [integer] of phenotypes to simulate .
#' @param N Number [integer] of samples to simulate; has to be provided as a 
#' dimnesionality check for kinship and downstream analyses; nrow(kinship) has 
#' to be equal to N.
#' @param kinship [N x N] Matrix of kinship estimates [double].
#' @param independent [bool] independent effect simulated if set to TRUE.
#' @param shared [bool] shared effect simulated if set to TRUE; at least one of 
#' shared or independent has to be set to TRUE.
#' @param phenoID Prefix [string] for naming traits.
#' @param id_samples Vector of [NrSamples] sample IDs [string]; if not provided
#' colnames(kinship) are used.
#' @param id_phenos Vector of [NrTraits] phenotype IDs [string]; if not provided
#' constructed by paste(phenoID, 1:P, sep="").
#' @return Named list of shared infinitesimal genetic effects (shared: [N x P] 
#' matrix) and independent infinitesimal genetic effects (independent: [N x P] 
#' matrix), the covariance term of the shared effect (cov_shared: [P x P] 
#' matrix), the covariance term of the independent effect (cov_independent: 
#' [P x P] matrix), the eigenvectors (eigenvec_kinship: [N x N]) and eigenvalues
#'  (eigenval_kinship: [N]) of the kinship matrix.
#' @details For the simulation of the infinitesimal genetic effects, three 
#' matrix components are used: i) the kinship matrix K [N x N] which is treated 
#' as the sample design matrix, ii) matrix B [N x P] with vec(B) drawn from a 
#' normal distribution and iii) the trait design matrix A [P x P]. For the
#' independent effect, A is a diagonal matrix with normally distributed values.  
#' A for the shared effect is a matrix of rowrank one, with normally distributed 
#' entries in row 1 and zeros elsewhere. To construct the final effects, the 
#' three matrices are multiplied as: E = cholesky(K)BA^T. 
#' @export
#' @examples
#' genotypes <- simulateGenotypes(N=100, NrSNP=400, verbose=FALSE)
#' kinship <- getKinship(N=100, X=genotypes$genotypes, standardise=TRUE, 
#' verbose=FALSE)
#' geneticBg <- geneticBgEffects(N=100, P=10, kinship=kinship)
geneticBgEffects <- function(P, N, kinship, phenoID="Trait_",
                             id_samples = colnames(kinship), 
                             shared=TRUE, independent=TRUE,
                             id_phenos = NULL) {
    numbers <- list(N=N, P=P)
    positives <- list(P=P, N=N)
    integers <- list(P=P, N=N)
    testNumerics(numbers=numbers, positives=positives, integers=integers)
    if (diff(dim(kinship)) !=0) {
        stop ("Kinship matrix needs to be a square matrix ,",
              "however it has ", nrow(kinship), " rows and ", ncol(kinship), 
              " columns")
    }
    if (N != ncol(kinship)) {
        stop ("The number of samples specified (", N, ") and the number of ",
              "samples in the kinship matrix (", ncol(kinship), 
              ") are different.")
    }
    if (length(id_samples) !=  N) {
        stop("Length of id_samples (", length(id_samples), ") is different ",
             "from N (", N, "). Does you kinship matrix contain column names ",
             "(default for retrieving id_samples)?")
    }
    if (!is.null(id_phenos) && length(id_phenos) !=  P) {
        stop("Length of id_phenos (", length(id_phenos), ") is different ",
             "from P (", P, ")")
    }
    if (!(is.character(phenoID) && length(phenoID) == 1)) {
        stop("phenoID has to be of length 1 and of type character")
    }
    if (is.null(id_phenos)) id_phenos <- paste(phenoID, 1:P, sep="")
    
    kin_eigen <- eigen(kinship, symmetric=TRUE)
    if (any(Re(as.complex(kin_eigen$values)) <= 0)) {
        stop("Kinship matrix is not positive semi-definite")
    }
    kinship_chol <- t(chol(kinship))
    N <- ncol(kinship)
    
    if (!shared && !independent) {
        stop("At least one geneticBgEffect has to be specified; set either ", 
             "shared or independent or both to TRUE")
    }
    # shared effect
    if (shared) {
        B <- matrix(rnorm(N * P), ncol=P)
        A <- matrix(rep(0, P * P), ncol=P)
        A[,1] <- rnorm(P)
        genBgShared <- kinship_chol %*% (B %*% t(A))
        colnames(genBgShared) <- id_phenos
        rownames(genBgShared) <- id_samples
        
        cov_shared <- tcrossprod(A)
        diag(cov_shared) <- diag(cov_shared) + 1e-4
    } else {
        genBgShared <- NULL
        cov_shared <- NULL
    }
    
    # independent effect
    if (independent) {
        D <- matrix(rnorm(N * P), ncol=P)
        C <- matrix(rep(0, P * P), ncol=P)
        diag(C) <- rnorm(P)
        genBgIndependent <- kinship_chol %*% (D %*% C)
        colnames(genBgIndependent) <- id_phenos
        rownames(genBgIndependent) <- id_samples
        
        cov_independent <- tcrossprod(C)
        diag(cov_independent) <- diag(cov_independent) + 1e-4
    } else {  
        genBgIndependent <- NULL
        cov_independent <- NULL
    }
    return(list(shared=genBgShared, independent=genBgIndependent, 
                cov_shared=cov_shared, cov_independent=cov_independent,
                evec_kinship=kin_eigen$vectors,
                eval_kinship=kin_eigen$values))
}

#' Simulate observational noise effects.
#'
#' noiseBgEffects simulates observational noise with a proportion of the effect 
#' shared across samples and a proportion independent across samples.
#'
#' @param P Number [integer] of phenotypes to simulate. 
#' @param N Number [integer] of samples to simulate.
#' @param mean Mean [double] of the normal distribution.
#' @param sd Standard deviation [double] of the normal distribution.
#' @param independent [bool] independent effect simulated if set to TRUE.
#' @param shared [bool] shared effect simulated if set to TRUE; at least one of 
#' shared or independent has to be set to TRUE.
#' @param sampleID Prefix [string] for naming samples.
#' @param phenoID Prefix [string] for naming traits.
#' @param id_samples Vector of [NrSamples] sample IDs [string]; if not provided
#' constructed by paste(sampleID, 1:N, sep="").
#' @param id_phenos Vector of [NrTraits] phenotype IDs [string]; if not provided
#' constructed by paste(phenoID, 1:P, sep="").
#' @return Named list of shared noise effects (shared: [N x P] matrix) and 
#' independent noise effects (independent: [N x P] matrix), the covariance term 
#' of the shared effect (cov_shared: [P x P] matrix) and the covariance term of 
#' the independent effect (cov_independent: [P x P] matrix).
#' @details For the simulation of the observational noise effects, two 
#' components are used: i) matrix B [N x P] with vec(B) drawn from a normal 
#' distribution with mean=mean and sd=sd and ii) the trait design matrix A 
#' [P x P]. For the independent effect, A is a diagonal matrix with normally 
#' distributed values.  A for the shared effect is a matrix of rowrank one, with 
#' normally distributed entries in row 1 and zeros elsewhere.  To construct the 
#' final effects, the two matrices are multiplied as: E = BA^T. 
#' @export
#' @examples
#' noiseBG <- noiseBgEffects(N=100, P=20, mean=2)
noiseBgEffects <- function(N, P, mean=0, sd=1, sampleID="ID_", phenoID="Trait_",
                           shared=TRUE, independent=TRUE,
                           id_samples = NULL,
                           id_phenos = NULL) {
    numbers <- list(N=N, P=P, sd)
    positives <- list(P=P, N=N, sd=sd)
    integers <- list(P=P, N=N)
    testNumerics(numbers=numbers, positives=positives, integers=integers)
    
    if (!is.null(id_samples) && length(id_samples) !=  N) {
        stop("Length of id_samples (", length(id_samples), ") is different ",
             "from N (", N, ")")
    }
    if (!is.null(id_phenos) && length(id_phenos) !=  P) {
        stop("Length of id_phenos (", length(id_phenos), ") is different ",
             "from P (", P, ")")
    }
    if (!(is.character(phenoID) && length(phenoID) == 1)) {
        stop("phenoID has to be of length 1 and of type character")
    }
    if (!(is.character(sampleID) && length(sampleID) == 1)) {
        stop("sampleID has to be of length 1 and of type character")
    }
    if (is.null(id_samples)) id_samples <- paste(sampleID, 1:N, sep="")
    if (is.null(id_phenos)) id_phenos <- paste(phenoID, 1:P, sep="")
    
    if (!shared && !independent) {
        stop("At least one noiseBgEffect has to be specified; set either ", 
             "shared or independent or both to TRUE")
    }
    # shared effect
    if (shared) {
        B <- matrix(rnorm(N * P, mean=mean, sd=sd), ncol=P)
        A <- matrix(rep(0, P * P), ncol=P)
        A[,1] <- rnorm(P, mean=mean, sd=sd)
        noiseBgShared <-  B %*% t(A)
        colnames(noiseBgShared) <- id_phenos
        rownames(noiseBgShared) <- id_samples
        
        cov_shared <- tcrossprod(A)
        diag(cov_shared) <- diag(cov_shared) + 1e-4
    } else {
        noiseBgShared <- NULL
        cov_shared <- NULL
    }
    
    # independent effect
    if (independent) {
        D <- matrix(rnorm(N * P, mean=mean, sd=sd), ncol=P)
        C <- matrix(rep(0, P * P), ncol=P)
        diag(C) <- rnorm(P, mean=mean, sd=sd)
        noiseBgIndependent <- D %*% C
        colnames(noiseBgIndependent) <- id_phenos
        rownames(noiseBgIndependent) <- id_samples
        
        cov_independent <- tcrossprod(C)
        diag(cov_independent) <- diag(cov_independent) + 1e-4
    } else {
        noiseBgIndependent <- NULL
        cov_independent <- NULL
    }
    return(list(shared=noiseBgShared, independent=noiseBgIndependent,
                cov_shared=cov_shared, cov_independent=cov_independent))
}

#' Simulate correlated background effects.
#'
#' correlatedBgEffects computes a background effect that simulates structured 
#' correlation between the phenotypes.
#' 
#' @param N Number [integer] of samples to simulate.
#' @param P Number [integer] of phenotypes to simulate. 
#' @param corr_mat [P x P] correlation matrix [double] as covariance 
#' component for the multivariate normal distribution. If not provided, pcorr is
#' used to construct the correlation matrix.
#' @param pcorr Initial strength of correlation [double] between neighbouring 
#' traits. Decreases by pcorr^(distance); distance from 0 to P-1. See details.
#' @param sampleID Prefix [string] for naming samples.
#' @param phenoID Prefix [string] for naming traits.
#' @param id_samples Vector of [NrSamples] sample IDs [string]; if not provided
#' constructed by paste(sampleID, 1:N, sep="").
#' @param id_phenos Vector of [NrTraits] phenotype IDs [string]; if not provided
#' constructed by paste(phenoID, 1:P, sep="").
#' @param verbose [boolean] If TRUE, progress info is printed to standard out.
#' @return Named list with [N x P] matrix of correlated background effects (
#' correlatedBg) and the correlation matrix (cov_correlated). If corr_mat 
#' provided corr_mat == cov_correlated.
#' @seealso \code{\link[mvtnorm]{rmvnorm}} which is used to simulate the 
#' multivariate normal distribution
#' @details correlatedBgEffects can be used to simulate phenotypes with a 
#' defined level of correlation between traits. If the corr_mat is not provided,
#' a simple correlation structure based on the distance of the traits will be 
#' constructed. Traits of distance d=1 (adjacent columns) will have correlation 
#' cor=\eqn{pcorr^1}{pcorr^1}, traits with d=2 have cor=\eqn{pcorr^2}{pcorr^2} 
#' up to traits with d=(P-1) cor=\eqn{pcorr^{(P-1)}}{pcorr^{(P-1)}} and 
#' 0 < pcorr < 1. The correlated background effect correlated is simulated based 
#' on this correlation structure C: 
#' \eqn{correlated ~ N_{NP}(0,C)}{correlated ~ N_{NP}(0,C)}.  
#' @export
#' @examples
#' correlatedBg <- correlatedBgEffects(N=100, P=20, pcorr=0.4)
correlatedBgEffects <- function(N, P, pcorr=NULL, corr_mat=NULL,
                                sampleID="ID_", phenoID="Trait_",
                                id_samples = NULL,
                                id_phenos = NULL,
                                verbose=FALSE) {
    numbers <- list(N=N, P=P)
    positives <- list(P=P, N=N)
    testNumerics(numbers=numbers, positives=positives)
    if (!is.null(id_samples) && length(id_samples) !=  N) {
        stop("Length of id_samples (", length(id_samples), ") is different ",
             "from N (", N, ")")
    }
    if (!is.null(id_phenos) && length(id_phenos) !=  P) {
        stop("Length of id_phenos (", length(id_phenos), ") is different ",
             "from P (", P, ")")
    }
    if (!(is.character(phenoID) && length(phenoID) == 1)) {
        stop("phenoID has to be of length 1 and of type character")
    }
    if (!(is.character(sampleID) && length(sampleID) == 1)) {
        stop("sampleID has to be of length 1 and of type character")
    }
    if (is.null(id_samples)) id_samples <- paste(sampleID, 1:N, sep="")
    if (is.null(id_phenos)) id_phenos <- paste(phenoID, 1:P, sep="")
    
    if(!is.null(corr_mat) && !is.null(pcorr)) {
        vmessage(c("Both pcorr (", pcorr, ") and corr_mat provided; corr_mat",
                   "will be used to simulate the correlatedBgEffects"), 
                 verbose=verbose)
    }
    if(is.null(corr_mat)) {
        if (is.null(pcorr)) {
            stop("At least one of pcorr or corr_mat have to be provided to", 
                 " simulated the correlatedBgEffects")
        }
        if (!is.numeric(pcorr)) {
            stop("pcorr has to be of type numeric")
        }
        if (length(pcorr) > 1) {
            stop("only a single pcorr value (double) can be specified")
        }
        if (pcorr <= 0 || pcorr >= 1) {
            stop("pcorr has to be greater than zero and less than 1")
        }
        corr_vec <- cumprod(rep(pcorr, P - 1))
        tri_corr_vec <- unlist(sapply(0:(length(corr_vec) -1), function(pos) {
            return(corr_vec[1:(length(corr_vec) - pos)])
        }))

        corr_mat <- diag(P)
        corr_mat[lower.tri(corr_mat, diag=FALSE)] <- tri_corr_vec
        corr_mat <- corr_mat + t(corr_mat) - diag(P)
    } else {
        corr_mat <- as.matrix(corr_mat)
        if (diff(dim(corr_mat)) != 0) {
            stop(paste("Correlation matrix needs to be a square matrix,",
                       "however it has", nrow(corr_mat), "rows and", 
                       ncol(corr_mat), "columns"))
        }
        if (ncol(corr_mat) != P) {
            stop("Dimensions of correlation matrix (", 
                 paste(dim(corr_mat), collapse=","), ") do ", 
                 "not match the number phenotypes (", P, ") to be simulated")
        }
    }
    pheno <- mvtnorm::rmvnorm(N, rep(0,P), corr_mat)
    colnames(pheno) <- id_phenos
    rownames(pheno) <- id_samples
    
    return(list(correlatedBg=pheno, cov_correlated=corr_mat))
}


