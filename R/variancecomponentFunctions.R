### 1. fixed effects

#' Simulate genetic fixed effects.
#'
#' geneticFixedEffects takes a matrix of SNPs which should be added as fixed 
#' effect to the phenotype. SNPs can have effects across all traits (common) 
#' or to a number of traits only (specific); the proportion of specific SNPs 
#' from the total input SNPs can be chosen via pSpecificGenetic. The number of 
#' traits 
#' that are associated with specific genetic effects can be chosen via 
#' pTraitSpecificGenetic. 
#'
#' @param X_causal [N x NrCausalSNPs] matrix of standardised (depending on 
#' standardise option) SNPs 
#' @param N number [integer] of samples to simulate
#' @param P number [integer] of phenotypes to simulate 
#' @param pSpecificGenetic Proportion [double] of genetic effects (SNPs) to have
#'  a trait-specific fixed effect
#' @param pTraitSpecificGenetic Proportion [double] of traits influenced by 
#' specific fixed genetic effects
#' @return named list of common fixed genetic effects (common: [N x P] matrix), 
#' specific fixed genetic effects (specific: [N x P] matrix), 
#' the causal SNPs named by having a common or specific effect 
#' (cov: [NrCausalSNPs x N] matrix) and the simulated effect sizes of the causal 
#' SNPs (cov_effect: [P x NrCausalSNPs] dataframe)
#' @export
#' @examples
#' genotypes <- simulateGenotypes(N=100, NrSNP=20, verbose=FALSE)
#' causalSNPs <- getCausalSNPs(genotypes=genotypes)
#' geneticFixed <- geneticFixedEffects(X_causal=causalSNPs, P=10, N=100)
geneticFixedEffects <- function(X_causal, N, P, pSpecificGenetic=0.4, 
                                pTraitSpecificGenetic=0.2) {
    NrCausalSNPs <- ncol(X_causal)
    NrSpecificSNPs <- round(pSpecificGenetic * NrCausalSNPs)
    NrCommonSNPs <- NrCausalSNPs - NrSpecificSNPs
    
    ### i) common
    comm <- sample(c(rep(TRUE, NrCommonSNPs), rep(FALSE, NrSpecificSNPs)), 
                   replace=FALSE)
    X_comm <-  X_causal[,comm]
    betaX_comm <- simulateDist(NrCommonSNPs, dist="norm") %*% 
        t(simulateDist(P, dist="norm"))
   
    ### ii) specific
    spec <- !comm
    X_specific <- X_causal[,spec]
    betaX_specific <- matrix(simulateDist(P * NrSpecificSNPs, dist="norm"), 
                             ncol=P)
    p_nongenetic <- sample(c(rep(TRUE, pTraitSpecificGenetic * P), 
                             rep(FALSE, (1 - pTraitSpecificGenetic) * P)), 
                           replace=FALSE)
    betaX_specific[,p_nongenetic] <- matrix(rep(0, length(which(p_nongenetic)) 
                                                * NrSpecificSNPs), 
                                            nrow=NrSpecificSNPs)

    cov = rbind(t(X_comm), t(X_specific))
    cov_effect = data.frame(betaX_comm=t(betaX_comm), 
                            betaX_spec=t(betaX_specific))
    colnames(cov_effect) <- paste(colnames(cov_effect),"_", 
                                  rownames(cov), sep="")

    Gcomm = X_comm %*% betaX_comm
    Gspec = X_specific %*% betaX_specific
    
    return(list(common=Gcomm, specific=Gspec, cov=cov, cov_effect=cov_effect))
}

#' Simulate noise fixed effects.
#'
#' noiseFixedEffects simulates a specific number of fixed noise effects 
#' (confounders). Confounders can have effects across all traits (common) 
#' or to a number of traits only (specific); the proportion of specific 
#' confounders from the total of simulated confounders can be chosen via 
#' pSpecificConfounders. The number of traits that are associated with specific 
#' noise effects can be chosen via  pTraitSpecificConfounders. Confounders can 
#' be simulated as categorical variables or following a binomial, uniform or 
#' normal distribution. Effect sizes for the noise effects can be simulated from
#' a uniform or normal distribution. Multiple confounder sets drawn from 
#' different distributions/different parameters of the same distribution can be 
#' simulated by specifying NrFixedEffects and supplying the respective
#' distribution parameters (*Confounders and *Beta) explained below. If a model
#' parameter and its "String' version are provided, i.e. NrConfounders and 
#' NrConfoundersStrings, the "String' specification will be used!
#'
#' @param N number [integer] of samples to simulate
#' @param P number [integer] of phenotypes to simulate
#' @param NrFixedEffects number [integer] of different fixed effects to simulate;
#' allows to simulate fixed effects from different distributions or with 
#'  differen parameters
#' @param NrConfounders vector of number(s) [integer] of confounders to simulate
#' @param pSpecificConfounders vector of proportion(s) [double] of noise effects 
#' (confounders) to have a trait-specific effect
#' @param pTraitSpecificConfounders vector of proportion(s) [double] of traits 
#' influenced by specific fixed noise effects
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
#' @return named list of common fixed noise effects (common: [N x P] matrix), 
#' specific fixed noise effects (specific: [N x P] matrix), 
#' the causal SNPs named by having a common or specific effect 
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
#'  distributed confounders 
#'  noiseFE_binomialandNormalConfounders <- noiseFixedEffects(P=10, N=20, 
#'  NrFixedEffects=2, NrConfounders=c(2,2), distConfounders=c("bin", "norm"), 
#'  probConfounders=0.2)
noiseFixedEffects <- function(N, P, NrFixedEffects=1, NrConfounders=10, 
                              pSpecificConfounders=0.4, 
                              pTraitSpecificConfounders=0.2, 
                              distConfounders="norm", mConfounders=0, 
                              sdConfounders=1, catConfounders=NULL, 
                              probConfounders=NULL, distBeta="norm", mBeta=0, 
                              sdBeta=1, 
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
                              sdBetaString=NULL) {
    
    oneFixedEffectComponent <- function(N, P, NrConfounders, 
                                        pSpecificConfounders, 
                                        pTraitSpecificConfounders, 
                                        distConfounders, mConfounders, 
                                        sdConfounders, catConfounders, 
                                        probConfounders, distBeta, mBeta, 
                                        sdBeta) {
        Ccomm <- NULL
        Cspec <- NULL
        
        NrSpecificConfounders <- round(pSpecificConfounders * NrConfounders)
        NrCommonConfounders <- NrConfounders - NrSpecificConfounders
    
        if (NrCommonConfounders != 0) {
            comm <- matrix(simulateDist(N * NrCommonConfounders, 
                                        dist=distConfounders, m=mConfounders, 
                                        std=sdConfounders, 
                                        categories=catConfounders, 
                                        prob=probConfounders), 
                           ncol=NrCommonConfounders)
            colnames(comm) <- paste("commonConfounder_", distConfounders, 
                                    seq(1, NrCommonConfounders, 1), sep="")
        
            beta_comm <- simulateDist(NrCommonConfounders,  dist=distBeta, 
                                      m=mBeta, std=sdBeta) %*% 
                        t(simulateDist(P, dist=distBeta, m=mBeta, std=sdBeta))
            rownames(beta_comm) <- paste("commonConfounder_", distConfounders,
                                        "_Beta_", distBeta, 
                                         seq(1, NrCommonConfounders, 1), sep="")
            
            Ccomm <- comm %*% beta_comm
            cov <- t(data.frame(comm))
            cov_effect <- data.frame(t(beta_comm))
        } 
        if (NrSpecificConfounders != 0) {
            spec <- matrix(simulateDist(N * NrSpecificConfounders, 
                                        dist=distConfounders, m=mConfounders, 
                                        std=sdConfounders, 
                                        categories=catConfounders, 
                                        prob=probConfounders), 
                           ncol=NrSpecificConfounders)
            colnames(spec) <- paste("specificConfounder_", distConfounders, 
                                    seq(1, NrSpecificConfounders, 1), sep="")
            beta_spec <- matrix(simulateDist(P * NrSpecificConfounders, 
                                             dist=distBeta, m=mBeta, 
                                             std=sdBeta), ncol=P)
            rownames(beta_spec) <- paste("specificConfounder_", distConfounders,
                                        "_Beta_", distBeta, 
                                         seq(1, NrSpecificConfounders, 1), 
                                         sep="")
            p_nonconfounders <- sample(
                c(rep(TRUE, pTraitSpecificConfounders * P), 
                rep(FALSE, (1 - pTraitSpecificConfounders) * P)), replace=FALSE)
            beta_spec[,p_nonconfounders] <- matrix(
                rep(0, length(which(p_nonconfounders)) * NrSpecificConfounders),
                nrow=NrSpecificConfounders)
        
            Cspec <- spec %*% beta_spec
            cov <- t(data.frame(spec))
            cov_effect <- data.frame(t(beta_spec))
        }
        if (NrCommonConfounders != 0 && NrSpecificConfounders != 0) {
            cov <- t(data.frame(comm, spec))
            cov_effect <- data.frame(t(beta_comm), t(beta_spec))
        }
        return(list(common=Ccomm, specific=Cspec, cov=cov, 
                    cov_effect=cov_effect))
    }
    if (!is.null(pSpecificConfoundersString)) {
        pSpecificConfounders <- commaList2vector(pSpecificConfoundersString)
    }
    if (length(pSpecificConfounders) != 1 && 
        length(pSpecificConfounders) != NrFixedEffects ) {
        stop(paste("Length of pSpecificConfounders (", 
                    length(pSpecificConfounders), ")
                    doesn't match NrFixedEffects (", NrFixedEffects, ")"))
    }
    if (!is.null(pTraitSpecificConfoundersString)) {
        pTraitSpecificConfounders <- 
            commaList2vector(pTraitSpecificConfoundersString)
    }
    if (length(pTraitSpecificConfounders) != 1 && 
        length(pTraitSpecificConfounders) != NrFixedEffects ) {
            stop(paste("Length of pTraitSpecificConfounders (", 
                       length(pTraitSpecificConfounders), ") doesn't match 
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
        if (length(pSpecificConfounders) == 1) {
            pSpecificConfounders <- rep(pSpecificConfounders, NrFixedEffects)
        }
        if (length(pTraitSpecificConfounders) == 1) {
            pTraitSpecificConfounders <- rep(pTraitSpecificConfounders, 
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
                                    pSpecificConfounders[x], 
                                    pTraitSpecificConfounders[x], 
                                    distConfounders[x], 
                                    mConfounders[x], sdConfounders[x], 
                                    catConfounders[x], probConfounders[x], 
                                    distBeta[x], mBeta[x], sdBeta[x])
        })
        common <- addNonNulls(lapply(tmp, function(x) x$common))
        specific <- addNonNulls(lapply(tmp, function(x) x$specific))
        cov <- do.call(rbind, lapply(tmp, function(x) x$cov))
        cov_effect <- do.call(cbind, lapply(tmp, function(x) x$cov_effect))
        return(list(common=common, specific=specific, cov=cov, 
                    cov_effect=cov_effect))
    } else {
        oneFixedEffectComponent(N, P, NrConfounders, pSpecificConfounders, 
                            pTraitSpecificConfounders, distConfounders, 
                            mConfounders, sdConfounders, catConfounders, 
                            probConfounders, distBeta, mBeta, sdBeta) 
    } 
}


### 2. Background effects

#' Simulate genetic background effects.
#'
#' geneticBgEffects simulates two random genetic effects (common and specific) 
#' based on the kinship estimates of the (simulated) samples.
#'
#' @param N number [integer] of samples to simulate
#' @param P number [integer] of phenotypes to simulate 
#' @param kinship [N x N] matrix of kinship estimates [double]
#' @return named list of common background genetic effects (common: [N x P] 
#' matrix) and specific background genetic effects (specific: [N x P] matrix) 
#' @details For the simulation of the genetic background effects, three matrix 
#' components are used: i) the kinship matrix K [N x N] which is treated as the 
#' sample-design matrix
#' (the genetic profile of the samples), ii) an effect-size matrix B [N x P] 
#' with vec(B) drawn from a normal distribution and iii) the trait design matrix
#'  A [P x P]. For the
#' specific effect, A is a diagonal matrix with normally distributed values. A 
#' for the common effect is a matrix of rowrank one, with normally distributed 
#' entries in row 1 and
#' zeros elsewhere. The final effect E, the three matrices are multiplied as: 
#' E = KBA
#' @export
#' @examples
#' genotypes <- simulateGenotypes(N=100, NrSNP=400, verbose=FALSE)
#' kinship <- getKinship(genotypes$X_sd, norm=TRUE, verbose=FALSE)
#' geneticBg <- geneticBgEffects(N=100, P=10, kinship=kinship)
geneticBgEffects <- function(N, P, kinship) {
    kinship_chol <- t(chol(kinship))
    
    # specific effect
    B <- matrix(simulateDist(x=N * P, dist="norm"), ncol=P)
    A <- matrix(rep(0, P * P), ncol=P)
    diag(A) <- simulateDist(x=P, dist="norm")
    genBgSpecific <- kinship_chol %*% (B %*% A)
    
    # common effect
    B <- matrix(simulateDist(x=N * P, dist="norm"), ncol=P)
    A <- matrix(rep(0, P * P), ncol=P)
    A[,1] <- simulateDist(x=P, dist="norm")
    genBgCommon <- kinship_chol %*% (B %*% t(A))
    
    return(list(common=genBgCommon, specific=genBgSpecific))
}

#' Simulate noise background effects.
#'
#' noiseBgEffects simulates two random noise effects (common and specific).
#'
#' @param P number [integer] of phenotypes to simulate 
#' @param N number [integer] of samples to simulate
#' @param mean mean [double] of the normal distribution
#' @param sd standard deviation [double] of the normal distribution
#' @return named list of common background noise effects (common: [N x P] 
#' matrix) and specific background noise effects (specific: [N x P] matrix) 
#' @details The common background effect common as simulated as vec(common) ~ 
#' N(mean,sd). The specific background effect sepcific is simulated as the 
#' matrix product of two normal
#' distributions A [N x 1] and B [P x 1]: specific At(B)
#' @export
#' @examples
#' noiseBG <- noiseBgEffects(N=100, P=20, mean=2)
noiseBgEffects <- function(N, P, mean=0, sd=1) {
    # common effect
    noiseBgCommon <- matrix(simulateDist(x= N * P, dist="norm", m=mean, std=sd), 
                            ncol=P) 
    
    # specific effect
    noiseBgSpecific <- simulateDist(x=N, dist="norm", m=mean, std=sd) %*%  
        t(simulateDist(x=P, dist="norm", m=mean, std=sd))
    return(list(common=noiseBgCommon, specific=noiseBgSpecific))
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


