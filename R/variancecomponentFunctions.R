### 1. fixed effects

#' Simulate genetic fixed effects.
#'
#' geneticFixedEffects takes genetic variants which should be added as fixed 
#' effect to the phenotype. These variants can have the same effects across all 
#' traits (shared) or can be independent across traits (independent); in 
#' addition, only a certain proportion of traits can be affected by the 
#' genetic variants.
#' 
#' @param X_causal [N x NrCausalSNPs] matrix of [NrCausalSNPs] SNPs from [N] 
#' samples.
#' @param N number [integer] of samples to simulate; has to be provided as a 
#' dimnesionality check for X_causal and downstream analyses; nrow(X_causal) has 
#' to be equal to N.
#' @param P number [integer] of phenotypes to simulate.
#' @param phenoID prefix [string] for naming traits.
#' @param id_samples vector of [NrSamples] sample IDs [string]; if not provided
#' colnames(X_causal) used.
#' @param id_phenos vector of [NrTraits] phenotype IDs [string]; if not provided
#' constructed by paste(phenoID, 1:P, sep="").
#' @param pTraitsAffected proportion [double] of traits affected by the genetic
#' effect. For non-integer results of pTraitsAffected*P, the ceiling of the 
#' result is used. Allows to simulate for instance different levels of 
#' pleiotropy.
#' @param pIndependentGenetic Proportion [double] of genetic effects (SNPs) to
#' have an independent fixed effect.
#' @param pTraitIndependentGenetic Proportion [double] of traits influenced by 
#' independent fixed genetic effects
#' @param distBeta vector of name(s) [string] of distribution to use to simulate 
#' effect sizes of SNPs; one of "unif" or "norm".
#' @param mBeta vector of mean/midpoint [double] of normal/uniform distribution 
#' for effect sizes of SNPs.
#' @param sdBeta vector of standard deviation/distance from midpoint [double] 
#' of normal/uniform distribution for effect sizes of SNPs.
#' @param verbose [boolean]; if TRUE, progress info is printed to standard out
#' @return named list of shared fixed genetic effects (shared: [N x P] matrix), 
#' independent fixed genetic effects (independent: [N x P] matrix), 
#' the causal SNPs labeled as shared or independent effect 
#' (cov: [NrCausalSNPs x N] matrix) and the simulated effect sizes of the causal 
#' SNPs (cov_effect: [P x NrCausalSNPs] dataframe).
#' @export
#' @examples
#' genotypes <- simulateGenotypes(N=100, NrSNP=20, verbose=FALSE)
#' causalSNPs <- getCausalSNPs(N=100, genotypes=genotypes$genotypes)
#' geneticFixed <- geneticFixedEffects(N=100, X_causal=causalSNPs$causalSNPs, 
#' P=10)
geneticFixedEffects <- function(X_causal, P, N, phenoID="Trait_",
                                id_samples = rownames(X_causal),
                                id_phenos = paste(phenoID, 1:P, sep=""),
                                pTraitsAffected=1,
                                pIndependentGenetic=0.4, 
                                pTraitIndependentGenetic=0.2, 
                                distBeta="norm", mBeta=0, sdBeta=1, 
                                verbose=FALSE) {

    if(nrow(X_causal) != N){
        stop("Number of samples in SNP matrix (", nrow(X_causal), ") is 
             different from number of samples to be simulated")
    }
    
    NrCausalSNPs <- ncol(X_causal)
    traitsAffected <- ceiling(P*pTraitsAffected)
    
    if (traitsAffected == 1) {
        NrIndependentSNPs <- NrCausalSNPs
    } else {
        NrIndependentSNPs <- round(pIndependentGenetic * NrCausalSNPs)
    }
    NrSharedSNPs <- NrCausalSNPs - NrIndependentSNPs
    
    vmessage(c("Out of", P, "total phenotypes, ", traitsAffected, "traits 
                will be affected by the fixed genetic effects. "))
    
    if (NrIndependentSNPs != 0) {
        vmessage(c("Out of the these affected traits (", traitsAffected, ")", 
                ceiling(pTraitIndependentGenetic * traitsAffected), "trait(s) 
                will have independent genetic effects"))
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

        betaX_shared <- simulateDist(NrSharedSNPs, dist=distBeta, m=mBeta, 
                                     std=sdBeta) %*% 
                        t(simulateDist(traitsAffected, dist=distBeta, m=mBeta, 
                                       std=sdBeta))
        if (P != traitsAffected) {
            betaX_shared <- cbind(betaX_shared, 
                                  matrix(0, ncol=P-traitsAffected, 
                                         nrow=NrSharedSNPs))
        }
        cov <- data.frame(X_shared)
        colnames(cov) <- snpIDshared
        rownames(cov) <- id_samples
        
        cov_effect <- data.frame(t(betaX_shared))
        colnames(cov_effect) <- paste(colnames(cov_effect), "_", 
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
        p_nongenetic <- sample(c(rep(FALSE, TraitIndependentGenetic), 
                                 rep(TRUE, 
                                     (traitsAffected - TraitIndependentGenetic))), 
                               replace=FALSE)
        betaX_independent[,p_nongenetic] <- 
            matrix(rep(0,  length(which(p_nongenetic)) * NrIndependentSNPs), 
                                                nrow=NrIndependentSNPs)
        
        if (P != traitsAffected) {
            betaX_independent <- cbind(betaX_independent, 
                                       matrix(0, ncol=P-traitsAffected, 
                                              nrow=NrIndependentSNPs))
        }
        
        cov <- data.frame(X_independent)
        colnames(cov) <- snpIDindependent
        rownames(cov) <- id_samples
        
        cov_effect <- data.frame(t(betaX_independent))
        colnames(cov_effect) <- paste(colnames(cov_effect), "_", 
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
                cov=cov, 
                cov_effect=cov_effect))
}

#' Simulate noise fixed effects.
#'
#' noiseFixedEffects simulates a number of fixed noise effects 
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
#' @param N number [integer] of samples to simulate.
#' @param P number [integer] of phenotypes to simulate.
#' @param sampleID prefix [string] for naming samples.
#' @param phenoID prefix [string] for naming traits.
#' @param id_samples vector of [NrSamples] sample IDs [string]; if not provided
#' constructed by paste(sampleID, 1:N, sep="").
#' @param id_phenos vector of [NrTraits] phenotype IDs [string]; if not provided
#' constructed by paste(phenoID, 1:P, sep="").
#' @param pTraitsAffected vector of proportion(s) [double] of traits affected by 
#' the confounders. For non-integer results of pTraitsAffected*P, the ceiling of 
#' the result is used.
#' @param NrFixedEffects number [integer] of different fixed effects to simulate
#' ; allows to simulate fixed effects from different distributions or with 
#'  different parameters; if only one type of confounder distribution is wanted,
#'  set NrFixedEffects=1 and choose the number of confounders with eg 
#'  NrConfounders=10.
#' @param NrConfounders vector of number(s) [integer] of confounders from a 
#' specified distribution to simulate.
#' @param pIndependentConfounders vector of proportion(s) [double] of noise 
#' effects (confounders) to have a trait-independent effect.
#' @param pTraitIndependentConfounders vector of proportion(s) [double] of 
#' traits influenced by independent fixed noise effects
#' @param distConfounders vector of name(s) [string] of distribution to use to 
#' simulate confounders; one of "unif", "norm", "bin", "cat_norm", "cat_unif".
#' @param distBeta vector of name(s) [string] of distribution to use to simulate 
#' effect sizes of confounders; one of "unif" or "norm".
#' @param mConfounders vector of mean/midpoint(s) [double] of normal/uniform 
#' distribution for confounders.
#' @param sdConfounders vector of standard deviation(s)/distance from 
#' midpoint(s) [double] of normal/uniform distribution for confounders.
#' @param catConfounders vector of number(s) of confounder categories [integer]; 
#' required if distConfounders "cat_norm" or "cat_unif".
#' @param probConfounders vector of probability(s) [double] of binomial 
#' confounders (0/1); required if distConfounders "bin".
#' @param mBeta vector of mean/midpoint [double] of normal/uniform distribution 
#' for effect sizes of confounders.
#' @param sdBeta vector of standard deviation/distance from midpoint [double] 
#' of normal/uniform distribution for effect sizes of confounders.
#' @return named list of shared fixed noise effects (shared: [N x P] matrix), 
#' independent fixed noise effects (independent: [N x P] matrix), 
#' the causal SNPs labeled shared or independent effect 
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
noiseFixedEffects <- function(N, P, sampleID="ID_",phenoID="Trait_", 
                              id_samples = paste(sampleID, 1:N, sep=""),
                              id_phenos = paste(phenoID, 1:P, sep=""),
                              pTraitsAffected=1,
                              NrFixedEffects=1, NrConfounders=10, 
                              pIndependentConfounders=0.4, 
                              pTraitIndependentConfounders=0.2, 
                              distConfounders="norm", mConfounders=0, 
                              sdConfounders=1, catConfounders=NULL, 
                              probConfounders=NULL, distBeta="norm", mBeta=0, 
                              sdBeta=1) {
    
    oneFixedEffectComponent <- function(N, P, NrConfounders, 
                                        pTraitsAffected,
                                        pIndependentConfounders, 
                                        pTraitIndependentConfounders, 
                                        distConfounders, mConfounders, 
                                        sdConfounders, catConfounders, 
                                        probConfounders, distBeta, mBeta, 
                                        sdBeta) {
        if (is.null(catConfounders) && grepl("cat", distConfounders)) {
            stop(paste("Confounder distribution set to", distConfounders, "but",
                       "no categories provided"))
        }
        if (is.null(probConfounders) && grepl("bin", distConfounders)) {
            stop(paste("Confounder distribution set to", distConfounders, "but",
                       "no probabilities provided"))
        }
        
        traitsAffected <- ceiling(P*pTraitsAffected)
        
        Cshared <- NULL
        Cindependent <- NULL
        if (P == 1) {
            NrIndependentConfounders <- 0
        } else {
            NrIndependentConfounders <- round(pIndependentConfounders * 
                                                  NrConfounders)
        }
        NrSharedConfounders <- NrConfounders - NrIndependentConfounders
        
        if (NrSharedConfounders != 0) {
            shared <- matrix(simulateDist(N * NrSharedConfounders, 
                                        dist=distConfounders, m=mConfounders, 
                                        std=sdConfounders, 
                                        categories=catConfounders, 
                                        prob=probConfounders), 
                           ncol=NrSharedConfounders)
            colnames(shared) <- paste("sharedConfounder_", distConfounders, 
                                    seq(1, NrSharedConfounders, 1), sep="")
        
            beta_shared <- simulateDist(NrSharedConfounders,  dist=distBeta, 
                                      m=mBeta, std=sdBeta) %*% 
                        t(simulateDist(traitsAffected, dist=distBeta, m=mBeta, 
                                       std=sdBeta))
            
            rownames(beta_shared) <- paste("sharedConfounder_", distConfounders,
                                        "_Beta_", distBeta, 
                                         seq(1, NrSharedConfounders, 1), sep="")
            if (P != traitsAffected) {
                beta_shared <- cbind(beta_shared, 
                                      matrix(0, ncol=P-traitsAffected, 
                                             nrow=NrSharedConfounders))
            }
            
            Cshared <- shared %*% beta_shared
            cov <- data.frame(shared)
            rownames(cov) <- id_samples
            cov_effect <- data.frame(t(beta_shared))
            rownames(cov_effect) <- id_phenos
        } 
        if (NrIndependentConfounders != 0) {
            independent <- matrix(simulateDist(N * NrIndependentConfounders, 
                                        dist=distConfounders, m=mConfounders, 
                                        std=sdConfounders, 
                                        categories=catConfounders, 
                                        prob=probConfounders), 
                           ncol=NrIndependentConfounders)
            colnames(independent) <- paste("independentConfounder_", 
                                           distConfounders, 
                                           seq(1, NrIndependentConfounders, 1), 
                                           sep="")
            beta_independent <- matrix(simulateDist(traitsAffected * 
                                                        NrIndependentConfounders, 
                                             dist=distBeta, m=mBeta, 
                                             std=sdBeta), ncol=traitsAffected)
            rownames(beta_independent) <- paste("independentConfounder_", 
                                                distConfounders,
                                                "_Beta_", distBeta, 
                                                seq(1, NrIndependentConfounders, 
                                                    1), 
                                                sep="")
            
            TraitIndependentConfounders <- ceiling(pTraitIndependentConfounders* 
                                                   traitsAffected)
            p_nonconfounders <- sample(
                c(rep(FALSE, TraitIndependentConfounders), 
                rep(TRUE, (traitsAffected - TraitIndependentConfounders))), 
                replace=FALSE)
            beta_independent[,p_nonconfounders] <- matrix(
                rep(0, 
                    length(which(p_nonconfounders)) * NrIndependentConfounders),
                nrow=NrIndependentConfounders)
            
            if (P != traitsAffected) {
                beta_independent <- cbind(beta_independent, 
                                          matrix(0, ncol=P-traitsAffected, 
                                                 nrow=NrIndependentConfounders))
            }
            
            
            Cindependent <- independent %*% beta_independent
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
    
    if (length(pTraitsAffected) != 1 && 
        length(pTraitsAffected) != NrFixedEffects ) {
        stop(paste("Length of pTraitsAffected (", 
                   length(pTraitsAffected), ")
                   doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (length(pIndependentConfounders) != 1 && 
        length(pIndependentConfounders) != NrFixedEffects ) {
        stop(paste("Length of pIndependentConfounders (", 
                    length(pIndependentConfounders), ")
                    doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (length(pTraitIndependentConfounders) != 1 && 
        length(pTraitIndependentConfounders) != NrFixedEffects ) {
            stop(paste("Length of pTraitIndependentConfounders (", 
                       length(pTraitIndependentConfounders), ") doesn't match 
                       NrFixedEffects (", NrFixedEffects, ")"))
    }        
    if (length(NrConfounders) != 1 && 
        length(NrConfounders) != NrFixedEffects ) {
        stop(paste("Length of NrConfounders (", length(NrConfounders), ")
                   doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (length(distConfounders) != 1 && 
        length(distConfounders) != NrFixedEffects ) {
        stop(paste("Length of distConfounders (", length(distConfounders)
            , ") doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (length(mConfounders) != 1 && 
        length(mConfounders) != NrFixedEffects ) {
        stop(paste("Length of mConfounders (", length(mConfounders), ")
                   doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (length(sdConfounders) != 1 && 
        length(sdConfounders) != NrFixedEffects ) {
        stop(paste("Length of sdConfounders (", length(sdConfounders), ")
                   doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (length(probConfounders) > 1 && 
        length(probConfounders) != NrFixedEffects ) {
        stop(paste("Length of probConfounders (", length(probConfounders),")
                   doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (length(catConfounders) > 1 && 
        length(catConfounders) != NrFixedEffects ) {
        stop(paste("Length of catConfounders (", length(catConfounders), ")
                   doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (length(distBeta) != 1 && 
        length(distBeta) != NrFixedEffects ) {
        stop(paste("Length of distBeta (", length(distBeta), ")
                   doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (length(mBeta) != 1 && 
        length(mBeta) != NrFixedEffects ) {
        stop(paste("Length of mBeta (", length(mBeta), ")
                   doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (length(sdBeta) != 1 && 
        length(sdBeta) != NrFixedEffects ) {
        stop(paste("Length of sdBeta (", length(sdBeta), ")
                   doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (NrFixedEffects != 1)  { 
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
            oneFixedEffectComponent(N, P, NrConfounders[x], 
                                    pTraitsAffected[x],
                                    pIndependentConfounders[x], 
                                    pTraitIndependentConfounders[x], 
                                    distConfounders[x], 
                                    mConfounders[x], sdConfounders[x], 
                                    catConfounders[x], probConfounders[x], 
                                    distBeta[x], mBeta[x], sdBeta[x])
        })
        shared <- addNonNulls(lapply(tmp, function(x) x$shared))
        independent <- addNonNulls(lapply(tmp, function(x) x$independent))
        cov <- do.call(cbind, lapply(tmp, function(x) x$cov))
        cov_effect <- do.call(cbind, lapply(tmp, function(x) x$cov_effect))
        return(list(shared=shared, independent=independent, cov=cov, 
                    cov_effect=cov_effect))
    } else {
        oneFixedEffectComponent(N, P, NrConfounders, pTraitsAffected,
                                pIndependentConfounders, 
                                pTraitIndependentConfounders, distConfounders, 
                                mConfounders, sdConfounders, catConfounders, 
                                probConfounders, distBeta, mBeta, sdBeta) 
    } 
}


### 2. Background effects

#' Simulate infinitesimal genetic effects (reflecting sample kinship).
#'
#' geneticBgEffects simulates an infinitesimal genetic effects with a proportion 
#' of the effect shared across samples and a proportion independent across 
#' samples; they are based on the kinship estimates of the (simulated) samples.
#'
#' @param P number [integer] of phenotypes to simulate .
#' @param N number [integer] of samples to simulate; has to be provided as a 
#' dimnesionality check for kinship and downstream analyses; nrow(kinship) has 
#' to be equal to N.
#' @param kinship [N x N] matrix of kinship estimates [double]
#' @param phenoID prefix [string] for naming traits.
#' @param id_samples vector of [NrSamples] sample IDs [string]; if not provided
#' colnames(kinship) are used.
#' @param id_phenos vector of [NrTraits] phenotype IDs [string]; if not provided
#' constructed by paste(phenoID, 1:P, sep="").
#' @return named list of shared infinitesimal genetic effects (shared: [N x P] 
#' matrix) and independent infinitesimal genetic effects (independent: [N x P] 
#' matrix), the covariance term of the shared effect (cov_shared: [P x P] 
#' matrix) and the covariance term of the independent effect (cov_independent: 
#' [P x P] matrix).
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
                             id_phenos = paste(phenoID, 1:P, sep="")) {
    kinship_chol <- t(chol(kinship))
    N <- ncol(kinship)
   
    # shared effect
    B <- matrix(rnorm(N * P), ncol=P)
    A <- matrix(rep(0, P * P), ncol=P)
    A[,1] <- rnorm(P)
    genBgShared <- kinship_chol %*% (B %*% t(A))
    colnames(genBgShared) <- id_phenos
    rownames(genBgShared) <- id_samples
    
    # independent effect
    D <- matrix(rnorm(N * P), ncol=P)
    C <- matrix(rep(0, P * P), ncol=P)
    diag(C) <- rnorm(P)
    genBgIndependent <- kinship_chol %*% (D %*% C)
    colnames(genBgIndependent) <- id_phenos
    rownames(genBgIndependent) <- id_samples
    
    cov_shared <- A%*%t(A)
    cov_independent <- C%*%t(C)

    return(list(shared=genBgShared, independent=genBgIndependent, 
                cov_shared=cov_shared, cov_independent=cov_independent))
}

#' Simulate observational noise effects.
#'
#' noiseBgEffects simulates observational noise with a proportion of the effect 
#' shared across samples and a proportion independent across samples.
#'
#' @param P number [integer] of phenotypes to simulate. 
#' @param N number [integer] of samples to simulate.
#' @param mean mean [double] of the normal distribution.
#' @param sd standard deviation [double] of the normal distribution.
#' #' @param sampleID prefix [string] for naming samples.
#' @param phenoID prefix [string] for naming traits.
#' @param id_samples vector of [NrSamples] sample IDs [string]; if not provided
#' constructed by paste(sampleID, 1:N, sep="").
#' @param id_phenos vector of [NrTraits] phenotype IDs [string]; if not provided
#' constructed by paste(phenoID, 1:P, sep="").
#' @return named list of shared noise effects (shared: [N x P] matrix) and 
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
                           id_samples = paste(sampleID, 1:N, sep=""),
                           id_phenos = paste(phenoID, 1:P, sep="")) {
    # shared effect
    B <- matrix(rnorm(N * P, mean=mean, sd=sd), ncol=P)
    A <- matrix(rep(0, P * P), ncol=P)
    A[,1] <- rnorm(P, mean=mean, sd=sd)
    noiseBgShared <-  B %*% t(A)
    colnames(noiseBgShared) <- id_phenos
    rownames(noiseBgShared) <- id_samples
    
    # independent effect
    D <- matrix(rnorm(N * P, mean=mean, sd=sd), ncol=P)
    C <- matrix(rep(0, P * P), ncol=P)
    diag(C) <- rnorm(P, mean=mean, sd=sd)
    noiseBgIndependent <- D %*% C
    colnames(noiseBgIndependent) <- id_phenos
    rownames(noiseBgIndependent) <- id_samples
    
    cov_shared <- A %*% t(A)
    cov_independent <- C %*% t(C)
    
    return(list(shared=noiseBgShared, independent=noiseBgIndependent,
                cov_shared=cov_shared, cov_independent=cov_independent))
}

#' Simulate correlated background effects.
#'
#' correlatedBgEffects computes a background effect that simulates structured 
#' correlation between the phenotypes.
#' 
#' @param N number [integer] of samples to simulate
#' @param P number [integer] of phenotypes to simulate 
#' @param pcorr initial strength of correlation [double] between neighbouring 
#' traits. Decreases by pcorr^(distance); distance from 0 to P-1. See details.
#' @param corr_mat optional [P x P] correlation matrix [double] as covariance 
#' component for the multivariate normal distribution. If not provided, pcorr is
#' used to construct the correlation matrix.
#' @param sampleID prefix [string] for naming samples.
#' @param phenoID prefix [string] for naming traits.
#' @param id_samples vector of [NrSamples] sample IDs [string]; if not provided
#' constructed by paste(sampleID, 1:N, sep="").
#' @param id_phenos vector of [NrTraits] phenotype IDs [string]; if not provided
#' constructed by paste(phenoID, 1:P, sep="").
#' @return named list with [N x P] matrix of correlated background effects (
#' correlatedBg) and the correlation matrix (cov_correlated). If corr_mat 
#' provided corr_mat == cov_correlated.
#' @seealso \code{\link[mvtnorm]{rmvnorm}} which is used to simulate the 
#' multivariate normal distribution
#' @details correlatedBgEffects can be used to simulate phenotypes with a 
#' defined level of correlation between traits. If the corr_mat is not provided,
#' a simple correlation structure based on the distance of the traits will be 
#' constructed. Traits of distance d=1 (adjacent columns) will have correlation 
#' cor=\eqn{pcorr^1}{pcorr^1}, traits with d=2 have cor=\eqn{pcorr^2}{pcorr^2} 
#' up to traits with d=(P-1) cor=\eqn{pcorr^{(P-1)}}{pcorr^{(P-1)}}. The 
#' correlated background effect correlated is simulated based on this 
#' correlation structure C: 
#' \eqn{correlated ~ N_{NP}(0,C)}{correlated ~ N_{NP}(0,C)}.  
#' @export
#' @examples
#' correlatedBg <- correlatedBgEffects(N=100, P=20, pcorr=0.4)
correlatedBgEffects <- function(N, P, pcorr, corr_mat=NULL,
                                sampleID="ID_", phenoID="Trait_",
                                id_samples = paste(sampleID, 1:N, sep=""),
                                id_phenos = paste(phenoID, 1:P, sep="")) {
        
    if(is.null(corr_mat)) {
        corr_vec <- cumprod(rep(pcorr, P - 1))
        tri_corr_vec <- unlist(sapply(0:(length(corr_vec) -1), function(pos) {
            return(corr_vec[1:(length(corr_vec) - pos)])
        }))

        corr_mat <- diag(P)
        corr_mat[lower.tri(corr_mat, diag=FALSE)] <- tri_corr_vec
        corr_mat <- corr_mat + t(corr_mat) - diag(P)
    }
    pheno <- mvtnorm::rmvnorm(N, rep(0,P), corr_mat)
    colnames(pheno) <- id_phenos
    rownames(pheno) <- id_samples
    
    return(list(correlatedBg=pheno, cov_correlated=corr_mat))
}


