### 1. fixed effects

#' Simulate genetic fixed effects.
#'
#' geneticFixedEffects takes a matrix of SNPs which should be added as fixed 
#' effect to the phenotype. SNPs can have effects across all traits (shared) 
#' or to a number of traits only (independent); the proportion of independent 
#' SNPs from the total input SNPs can be chosen via pIndependentGenetic. The 
#' number of traits that are associated with independent genetic effects can be 
#' chosen via pTraitIndependentGenetic. 
#'
#' @param X_causal [N x NrCausalSNPs] matrix of standardised (depending on 
#' standardise option) SNPs 
#' @param N number [integer] of samples to simulate; has to be less than or 
#' equal to the number of samples in X_causal (rows); if less than  number of 
#' samples in X_causal, N random rows of X_causal will be drawn; if not 
#' specified, N assumed to be equal to samples in X_causal
#' @param P number [integer] of phenotypes to simulate 
#' @param pIndependentGenetic Proportion [double] of genetic effects (SNPs) to
#'  have an independent fixed effect
#' @param pTraitIndependentGenetic Proportion [double] of traits influenced by 
#' independent fixed genetic effects
#' @param verbose [boolean]; if TRUE, progress info is printed to standard out
#' @return named list of shared fixed genetic effects (shared: [N x P] matrix), 
#' independent fixed genetic effects (independent: [N x P] matrix), 
#' the causal SNPs named by having a shared or independent effect 
#' (cov: [NrCausalSNPs x N] matrix) and the simulated effect sizes of the causal 
#' SNPs (cov_effect: [P x NrCausalSNPs] dataframe)
#' @export
#' @examples
#' genotypes <- simulateGenotypes(N=100, NrSNP=20, verbose=FALSE)
#' causalSNPs <- getCausalSNPs(genotypes=genotypes)
#' geneticFixed <- geneticFixedEffects(X_causal=causalSNPs, P=10, N=100)
geneticFixedEffects <- function(X_causal, P, N=NULL, pIndependentGenetic=0.4, 
                                pTraitIndependentGenetic=0.2, verbose=TRUE) {
    NrGenotypeSamples <- nrow(X_causal) 
    if (!is.null(N)) {
        if (N > NrGenotypeSamples) {
            stop("Sample number specified exceeds number of genotypes provided")
        }
        if (N < NrGenotypeSamples) {
            vmessage(c("Sampling", N, "samples from", NrGenotypeSamples,
                       "genotypes provided"))
            X_causal <- X_causal[sample(NrGenotypeSamples, N),]
        }
    }
    NrCausalSNPs <- ncol(X_causal)
    if (P == 1) {
        NrIndependentSNPs <- 0
    } else {
        NrIndependentSNPs <- round(pIndependentGenetic * NrCausalSNPs)
    }
    NrSharedSNPs <- NrCausalSNPs - NrIndependentSNPs
    
    Gshared <- NULL
    Gindependent <- NULL
   
    if (NrSharedSNPs != 0) {
        # shared
        shared <- sample(c(rep(TRUE, NrSharedSNPs), 
                           rep(FALSE, NrIndependentSNPs)), 
                   replace=FALSE)
        X_shared <-  X_causal[,shared]
        betaX_shared <- rnorm(NrSharedSNPs) %*% t(rnorm(P))
        cov <- t(data.frame(X_shared))
        cov_effect <- data.frame(t(betaX_shared))
        colnames(cov_effect) <- paste(colnames(cov_effect), "_", 
                                      rownames(cov), sep="")
        Gshared = X_shared %*% betaX_shared
    }
   
    if (NrIndependentSNPs != 0) {
        # independent
        independent <- !shared
        X_independent <- X_causal[,independent]
        betaX_independent <- matrix(rnorm(P * NrIndependentSNPs), ncol=P)
        TraitIndependentGenetic <- ceiling(pTraitIndependentGenetic * P)
        p_nongenetic <- sample(c(rep(TRUE, TraitIndependentGenetic), 
                                 rep(FALSE, 
                                     (P - TraitIndependentGenetic))), 
                               replace=FALSE)
        betaX_independent[,p_nongenetic] <- 
            matrix(rep(0,  TraitIndependentGenetic * NrIndependentSNPs), 
                                                nrow=NrIndependentSNPs)
        cov <- t(data.frame(independent))
        cov_effect <- data.frame(t(betaX_independent))
        colnames(cov_effect) <- paste(colnames(cov_effect), "_", 
                                      rownames(cov), sep="")
    }
    if (NrSharedSNPs != 0 && NrIndependentSNPs != 0) {
        cov = rbind(t(X_shared), t(X_independent))
        cov_effect = data.frame(betaX_shared=t(betaX_shared), 
                                betaX_independent=t(betaX_independent))
        colnames(cov_effect) <- paste(colnames(cov_effect), "_", 
                                      rownames(cov), sep="")
        
        Gindependent = X_independent %*% betaX_independent
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
#' or to a number of traits only (independent); the proportion of independent 
#' confounders from the total of simulated confounders can be chosen via 
#' pIndependentConfounders. The number of traits that are associated with 
#' independent noise effects can be chosen via  pTraitIndependentConfounders. 
#' Confounders can be simulated as categorical variables or following a binomial 
#' , uniform or normal distribution. Effect sizes for the noise effects can be 
#' simulated from a uniform or normal distribution. Multiple confounder sets 
#' drawn from different distributions/different parameters of the same 
#' distribution can be simulated by specifying NrFixedEffects and supplying the 
#' respective distribution parameters (*Confounders and *Beta) explained below. 
#' If a model parameter and its "String' version are provided, i.e. 
#' NrConfounders and NrConfoundersStrings, the "String' specification will be 
#' used!
#'
#' @param N number [integer] of samples to simulate
#' @param P number [integer] of phenotypes to simulate
#' @param NrFixedEffects number [integer] of different fixed effects to simulate
#' ; allows to simulate fixed effects from different distributions or with 
#'  differen parameters
#' @param NrConfounders vector of number(s) [integer] of confounders to simulate
#' @param pIndependentConfounders vector of proportion(s) [double] of noise 
#' effects (confounders) to have a trait-independent effect
#' @param pTraitIndependentConfounders vector of proportion(s) [double] of 
#' traits influenced by independent fixed noise effects
#' @param distConfounders vector of name(s) [string] of distribution to use to 
#' simulate confounders; one of "unif", "norm", "bin", "cat_norm", "cat_unif"
#' @param distBeta vector of name(s) [string] of distribution to use to simulate 
#' effect sizes of confounders; one of "unif" or "norm"
#' @param mConfounders vector of mean/midpoint(s) [double] of normal/uniform 
#' distribution for confounders
#' @param sdConfounders vector of standard deviation(s)/distance from 
#' midpoint(s) [double] of normal/uniform distribution for confounders
#' @param catConfounders vector of number(s) of confounder categories [integer]; 
#' required if distConfounders "cat_norm" or "cat_unif" 
#' @param probConfounders vector of probability(s) [double] of binomial 
#' confounders (0/1); required if distConfounders "bin" 
#' @param mBeta vector of mean/midpoint [double] of normal/uniform distribution 
#' for effect sizes of confounders
#' @param sdBeta vector of standard deviation/distance from midpoint [double] 
#' of normal/uniform distribution for effect sizes of confounders
#' @param NrConfoundersString alternative to NrConfounder, a comma-separated
#' [string] with number(s) [integer] of confounders to simulate; typically used 
#' when run as command line application
#' @param pIndependentConfoundersString alternative to pIndependentConfounders,  
#' a comma-separated [string] with proportion(s) [double] of noise effects 
#' (confounders) to have a trait-independent effect; typically used when run as 
#' command line application
#' @param pTraitIndependentConfoundersString alternative to 
#' pTraitIndependentConfounders,  a comma-separated [string] with proportion(s) 
#' [double] of traits influenced by independent fixed noise effects; typically 
#' used when run as command line application
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
#' @return named list of shared fixed noise effects (shared: [N x P] matrix), 
#' independent fixed noise effects (independent: [N x P] matrix), 
#' the causal SNPs named by having a shared or independent effect 
#' (cov: [NrConfounders x N] matrix) and the simulated effect sizes of the 
#' confounders
#' (cov_effect: [P x NrConfounders] dataframe)
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
#'  # 4 fixed noise effect with 2 binomial confounders and 2 normally 
#'  # distributed confounders 
#'  noiseFE_binomialandNormalConfounders <- noiseFixedEffects(P=10, N=20, 
#'  NrFixedEffects=2, NrConfounders=c(2,2), distConfounders=c("bin", "norm"), 
#'  probConfounders=0.2)
noiseFixedEffects <- function(N, P, NrFixedEffects=1, NrConfounders=10, 
                              pIndependentConfounders=0.4, 
                              pTraitIndependentConfounders=0.2, 
                              distConfounders="norm", mConfounders=0, 
                              sdConfounders=1, catConfounders=NULL, 
                              probConfounders=NULL, distBeta="norm", mBeta=0, 
                              sdBeta=1, 
                              NrConfoundersString=NULL, 
                              pIndependentConfoundersString=NULL, 
                              pTraitIndependentConfoundersString=NULL, 
                              distConfoundersString=NULL, 
                              mConfoundersString=NULL, 
                              sdConfoundersString=NULL, 
                              catConfoundersString=NULL, 
                              probConfoundersString=NULL, 
                              distBetaString=NULL, 
                              mBetaString=NULL, 
                              sdBetaString=NULL) {
    
    oneFixedEffectComponent <- function(N, P, NrConfounders, 
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
                        t(simulateDist(P, dist=distBeta, m=mBeta, std=sdBeta))
            rownames(beta_shared) <- paste("sharedConfounder_", distConfounders,
                                        "_Beta_", distBeta, 
                                         seq(1, NrSharedConfounders, 1), sep="")
            
            Cshared <- shared %*% beta_shared
            cov <- t(data.frame(shared))
            cov_effect <- data.frame(t(beta_shared))
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
            beta_independent <- matrix(simulateDist(P * NrIndependentConfounders, 
                                             dist=distBeta, m=mBeta, 
                                             std=sdBeta), ncol=P)
            rownames(beta_independent) <- paste("independentConfounder_", 
                                                distConfounders,
                                                "_Beta_", distBeta, 
                                                seq(1, NrIndependentConfounders, 
                                                    1), 
                                                sep="")
            p_nonconfounders <- sample(
                c(rep(TRUE, pTraitIndependentConfounders * P), 
                rep(FALSE, (1 - pTraitIndependentConfounders) * P)), 
                replace=FALSE)
            beta_independent[,p_nonconfounders] <- matrix(
                rep(0, 
                    length(which(p_nonconfounders)) * NrIndependentConfounders),
                nrow=NrIndependentConfounders)
        
            Cindependent <- independent %*% beta_independent
            cov <- t(data.frame(independent))
            cov_effect <- data.frame(t(beta_independent))
        }
        if (NrSharedConfounders != 0 && NrIndependentConfounders != 0) {
            cov <- t(data.frame(shared, independent))
            cov_effect <- data.frame(t(beta_shared), t(beta_independent))
        }
        return(list(shared=Cshared, independent=Cindependent, cov=cov, 
                    cov_effect=cov_effect))
    }
    
    if (!is.null(pIndependentConfoundersString)) {
        pIndependentConfounders <- commaList2vector(
            pIndependentConfoundersString)
    }
    if (length(pIndependentConfounders) != 1 && 
        length(pIndependentConfounders) != NrFixedEffects ) {
        stop(paste("Length of pIndependentConfounders (", 
                    length(pIndependentConfounders), ")
                    doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (!is.null(pTraitIndependentConfoundersString)) {
        pTraitIndependentConfounders <- 
            commaList2vector(pTraitIndependentConfoundersString)
    }
    if (length(pTraitIndependentConfounders) != 1 && 
        length(pTraitIndependentConfounders) != NrFixedEffects ) {
            stop(paste("Length of pTraitIndependentConfounders (", 
                       length(pTraitIndependentConfounders), ") doesn't match 
                       NrFixedEffects (", NrFixedEffects, ")"))
    }        
    if (!is.null(NrConfoundersString)) {
        NrConfounders <- commaList2vector(NrConfoundersString)
    }
    if (length(NrConfounders) != 1 && 
        length(NrConfounders) != NrFixedEffects ) {
        stop(paste("Length of NrConfounders (", length(NrConfounders), ")
                   doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (!is.null(distConfoundersString)) {
        distConfounders <- 
            unlist(strsplit(distConfoundersString, split=","))
    }
    if (length(distConfounders) != 1 && 
        length(distConfounders) != NrFixedEffects ) {
        stop(paste("Length of distConfounders (", length(distConfounders)
            , ") doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (!is.null(mConfoundersString)) {
        mConfounders <- commaList2vector(mConfoundersString)
    }
    if (length(mConfounders) != 1 && 
        length(mConfounders) != NrFixedEffects ) {
        stop(paste("Length of mConfounders (", length(mConfounders), ")
                   doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (!is.null(sdConfoundersString)) {
        sdConfounders <- commaList2vector(sdConfoundersString)
    }
    if (length(sdConfounders) != 1 && 
        length(sdConfounders) != NrFixedEffects ) {
        stop(paste("Length of sdConfounders (", length(sdConfounders), ")
                   doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (!is.null(probConfoundersString)) {
        probConfounders <- commaList2vector(probConfoundersString)
    }
    if (length(probConfounders) > 1 && 
        length(probConfounders) != NrFixedEffects ) {
        stop(paste("Length of probConfounders (", length(probConfounders),")
                   doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (!is.null(catConfoundersString)) {
        catConfounders <- commaList2vector(catConfoundersString)
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
                                    pIndependentConfounders[x], 
                                    pTraitIndependentConfounders[x], 
                                    distConfounders[x], 
                                    mConfounders[x], sdConfounders[x], 
                                    catConfounders[x], probConfounders[x], 
                                    distBeta[x], mBeta[x], sdBeta[x])
        })
        shared <- addNonNulls(lapply(tmp, function(x) x$shared))
        independent <- addNonNulls(lapply(tmp, function(x) x$independent))
        cov <- do.call(rbind, lapply(tmp, function(x) x$cov))
        cov_effect <- do.call(cbind, lapply(tmp, function(x) x$cov_effect))
        return(list(shared=shared, independent=independent, cov=cov, 
                    cov_effect=cov_effect))
    } else {
        oneFixedEffectComponent(N, P, NrConfounders, pIndependentConfounders, 
                            pTraitIndependentConfounders, distConfounders, 
                            mConfounders, sdConfounders, catConfounders, 
                            probConfounders, distBeta, mBeta, sdBeta) 
    } 
}


### 2. Background effects

#' Simulate genetic background effects.
#'
#' geneticBgEffects simulates two random genetic effects (shared and 
#' independent) based on the kinship estimates of the (simulated) samples.
#'
#' @param P number [integer] of phenotypes to simulate 
#' @param kinship [N x N] matrix of kinship estimates [double]
#' @return named list of shared background genetic effects (shared: [N x P] 
#' matrix) and independent background genetic effects (independent: [N x P] 
#' matrix) 
#' @details For the simulation of the genetic background effects, three matrix 
#' components are used: i) the kinship matrix K [N x N] which is treated as the 
#' sample-design matrix
#' (the genetic profile of the samples), ii) an effect-size matrix B [N x P] 
#' with vec(B) drawn from a normal distribution and iii) the trait design matrix
#'  A [P x P]. For the
#' independent effect, A is a diagonal matrix with normally distributed values.  
#' A for the shared effect is a matrix of rowrank one, with normally distributed 
#' entries in row 1 and
#' zeros elsewhere. The final effect E, the three matrices are multiplied as: 
#' E = KBA
#' @export
#' @examples
#' genotypes <- simulateGenotypes(N=100, NrSNP=400, verbose=FALSE)
#' kinship <- getKinship(genotypes$X_sd, norm=TRUE, verbose=FALSE)
#' geneticBg <- geneticBgEffects(P=10, kinship=kinship)
geneticBgEffects <- function(P, kinship) {
    N <- ncol(kinship)
    kinship_chol <- t(chol(kinship))
   
    # shared effect
    B <- matrix(rnorm(N * P), ncol=P)
    A <- matrix(rep(0, P * P), ncol=P)
    A[,1] <- rnorm(P)
    genBgShared <- kinship_chol %*% (B %*% t(A))
    
    # independent effect
    B <- matrix(rnorm(N * P), ncol=P)
    A <- matrix(rep(0, P * P), ncol=P)
    diag(A) <- rnorm(P)
    genBgIndependent <- kinship_chol %*% (B %*% A)

    return(list(shared=genBgShared, independent=genBgIndependent))
}

#' Simulate noise background effects.
#'
#' noiseBgEffects simulates two random noise effects (shared and independent).
#'
#' @param P number [integer] of phenotypes to simulate 
#' @param N number [integer] of samples to simulate
#' @param mean mean [double] of the normal distribution
#' @param sd standard deviation [double] of the normal distribution
#' @return named list of shared background noise effects (shared: [N x P] 
#' matrix) and independent background noise effects (independent: [N x P] 
#' matrix) 
#' @details The independent background effect is simulated as vec(shared) ~ 
#' N(mean,sd). The shared background effect is simulated as the 
#' matrix product of two normal
#' distributions A [N x 1] and B [P x 1]: AB^t
#' @export
#' @examples
#' noiseBG <- noiseBgEffects(N=100, P=20, mean=2)
noiseBgEffects <- function(N, P, mean=0, sd=1) {
    # shared effect
    noiseBgShared <- rnorm(n=N, mean=mean, sd=sd) %*% t(rnorm(n=P, 
                                                                   mean=mean, 
                                                                   sd=sd))
    
    # independent effect
    noiseBgIndependent <- matrix(rnorm(n= N * P, mean=mean, sd=sd), ncol=P) 
   
    return(list(shared=noiseBgShared, independent=noiseBgIndependent))
}

#' Simulate correlated background effects.
#'
#' correlatedBgEffects computes a background effect that simulates structured 
#' correlation between the phenotypes.
#' 
#' @param N number [integer] of samples to simulate
#' @param P number [integer] of phenotypes to simulate 
#' @param pcorr initial strength of correlation [double] between neighbouring 
#' traits
#' @return [N x P] matrix of correlated background effects
#' @seealso \code{\link[mvtnorm]{rmvnorm}} which is used to simulate the 
#' multivariate normal distribution
#' @details correlatedBgEffects can be used to simulate phenotypes with a 
#' defined level of correlation between traits. The level of correlation depends
#'  on the distance
#' of the traits. Traits of distance d=1 (adjacent columns) will have 
#' correlation cor=\eqn{pcorr^1}{pcorr^1}, traits with d=2 have 
#' cor=\eqn{pcorr^2}{pcorr^2} up to traits with d=(P-1) 
#' cor=\eqn{pcorr^{(P-1)}}{pcorr^{(P-1)}}. The correlated background effect 
#' correlated is simulated based on this correlation structure C: 
#' \eqn{correlated ~ N_{NP}(0,C)}{correlated ~ N_{NP}(0,C)}.  
#' @export
#' @examples
#' correlatedBg <- correlatedBgEffects(N=100, P=20, pcorr=0.4)
correlatedBgEffects <- function(N, P, pcorr) {
    corr_vec <- cumprod(rep(pcorr, P - 1))
    tri_corr_vec <- unlist(sapply(0:(length(corr_vec) -1), function(pos) {
        return(corr_vec[1:(length(corr_vec) - pos)])
    }))

    Sigma <- diag(P)
    Sigma[lower.tri(Sigma, diag=FALSE)] <- tri_corr_vec
    Sigma <- Sigma + t(Sigma) - diag(P)

    pheno <- mvtnorm::rmvnorm(N, rep(0,P), Sigma)

    return(pheno)
}


