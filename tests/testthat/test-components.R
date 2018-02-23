context('Test variance component functions')
geno <- simulateGenotypes(N=10, NrSNP=100, verbose=FALSE)

context('Test geneticFixedEffects')
test_that('geneticFixedEffects throws sample error',{
    expect_error(geneticFixedEffects(X_causal = geno$genotypes,
                                     P=10, N=20, verbose=FALSE),
                 "Number of samples in SNP matrix")
})

test_that('geneticFixedEffects throws id_samples error when no rownames',{
    rownames(geno$genotypes) <- NULL
    expect_error(geneticFixedEffects(X_causal = geno$genotypes,
                                     P=10, N=10,
                                     verbose=FALSE),
                 "Length of id_samples")
})

test_that('geneticFixedEffects throws id_samples error',{
    expect_error(geneticFixedEffects(X_causal = geno$genotypes,
                                     P=10, N=10, id_samples="sample1",
                                     verbose=FALSE),
                 "Length of id_samples")
})

test_that('geneticFixedEffects throws id_phenos error',{
    expect_error(geneticFixedEffects(X_causal = geno$genotypes,
                                     P=10, N=10, id_phenos="pheno1",
                                     verbose=FALSE),
                 "Length of id_phenos")
})

test_that('geneticFixedEffects throws sample error',{
    expect_error(geneticFixedEffects(X_causal = geno$genotypes,
                                     P=10, N=20, verbose=FALSE),
                 "Number of samples in SNP matrix")
})

test_that('geneticFixedEffects throws length error',{
    expect_error(geneticFixedEffects(X_causal=geno$genotypes, N=10, P=2, 
                                     phenoID=c("a", "b")),
                 "phenoID has to be of length 1 and of type character")
})

test_that('geneticFixedEffects throws type error',{
    expect_error(geneticFixedEffects(X_causal=geno$genotypes, N=10, P=2, 
                                     phenoID=2),
                "phenoID has to be of length 1 and of type character")
})

test_that('geneticFixedEffects throws pTraitAffected proportion error',{
    expect_error(geneticFixedEffects(X_causal=geno$genotypes, N=10, P=2, 
                                     pTraitsAffected = 4),
                 "Proportions have to be specified between 0 and 1:")
})

test_that('geneticFixedEffects throws pTraitAffected type error',{
    expect_error(geneticFixedEffects(X_causal=geno$genotypes, N=10, P=2, 
                                     pTraitsAffected = "all"),
                "pTraitsAffected is/are not numeric")
})

test_that('geneticFixedEffects throws pIndependentGenetic type error',{
    expect_error(geneticFixedEffects(X_causal=geno$genotypes, N=10, P=2, 
                                     pIndependentGenetic = "all"),
                "pIndependentGenetic is/are not numeric")
})

test_that('geneticFixedEffects throws pIndependentGenetic proportion error',{
    expect_error(geneticFixedEffects(X_causal=geno$genotypes, N=10, P=2, 
                                     pIndependentGenetic = -0.2),
                 "Proportions have to be specified between 0 and 1")
})

test_that('geneticFixedEffects throws pTraitIndependentGenetic proportion error',{
    expect_error(geneticFixedEffects(X_causal=geno$genotypes, N=10, P=2, 
                                     pTraitIndependentGenetic = 1.1), 
                "Proportions have to be specified between 0 and 1")
})

test_that('geneticFixedEffects throws pTraitIndependentGenetic type error',{
    expect_error(geneticFixedEffects(X_causal=geno$genotypes, N=10, P=2, 
                                     pTraitIndependentGenetic = "none"),
                 "pTraitIndependentGenetic is/are not numeric")
})

test_that('geneticFixedEffects throws type error',{
    expect_error(geneticFixedEffects(X_causal=geno$genotypes, N='a', P='b'),
                 "P,N is/are not numeric")
})

test_that('geneticFixedEffects throws sdBeta out of range error',{
    expect_error(geneticFixedEffects(X_causal=geno$genotypes, N=10, P=2, 
                                     sdBeta=-.05),
                "sdBeta has/have to be greater than zero")
})

test_that('geneticFixedEffects returns correct effect dimensions',{
    geno <- simulateGenotypes(N=100, NrSNP=50, verbose=FALSE)
    f <- geneticFixedEffects(X_causal = geno$genotypes,
                                     P=10, N=100, verbose=FALSE)
    expect_equal(dim(f$shared), c(100,10))
    expect_equal(dim(f$independent), c(100,10))
})

test_that('geneticFixedEffects returns correct cov and coveffectdimensions',{
    geno <- simulateGenotypes(N=100, NrSNP=50, verbose=FALSE)
    f <- geneticFixedEffects(X_causal = geno$genotypes,
                             P=10, N=100, verbose=FALSE)
    expect_equal(dim(f$cov), c(100,50))
    expect_equal(dim(f$cov_effect), c(10,50))
})

test_that('geneticFixedEffects returns NULL for non-specified components',{
    geno <- simulateGenotypes(N=100, NrSNP=50, verbose=FALSE)
    f <- geneticFixedEffects(X_causal = geno$genotypes,
                             P=10, N=100, verbose=FALSE,
                             pIndependentGenetic = 0)
    expect_equal(dim(f$independent), NULL)
    f <- geneticFixedEffects(X_causal = geno$genotypes,
                             P=10, N=100, verbose=FALSE,
                             pIndependentGenetic = 1)
    expect_equal(dim(f$shared), NULL)
})

test_that('geneticFixedEffects shared effect is correlated',{
    geno <- simulateGenotypes(N=100, NrSNP=50, verbose=FALSE)
    f <- geneticFixedEffects(X_causal = geno$genotypes,
                             P=10, N=100, verbose=FALSE)
    expect_true(all(abs(cor(f$shared)) > 0.9))
})

test_that('geneticFixedEffects independent effect is not perfectly correlated',{
    geno <- simulateGenotypes(N=100, NrSNP=50, verbose=FALSE)
    f <- geneticFixedEffects(X_causal = geno$genotypes,
                             P=10, N=100, verbose=FALSE,
                             pIndependentGenetic = 0.9)
    expect_false(all(cor(f$independent) == 1))
})

test_that("geneticFixedEffects returns correct number of independent effects",{
    pIndependentGenetic <- 0.5
    NrSNP <- 10
    geno <- simulateGenotypes(N=100, NrSNP=NrSNP, verbose=FALSE)
    f <- geneticFixedEffects(X_causal = geno$genotypes,
                             P=10, N=100, 
                             pIndependentGenetic = pIndependentGenetic, 
                             verbose=FALSE)
    expect_equal(length(which(grepl("independent", colnames(f$cov_effect)))),
                 pIndependentGenetic * NrSNP)
})

test_that("geneticFixedEffects fails with proportions error",{
    pIndependentGenetic <- -0.5
    NrSNP <- 10
    geno <- simulateGenotypes(N=100, NrSNP=NrSNP, verbose=FALSE)
    expect_error(geneticFixedEffects(X_causal = geno$genotypes,
                                     P=10, N=100, 
                                     pIndependentGenetic = pIndependentGenetic, 
                                     verbose=FALSE),
                 "Proportions have to be specified between 0 and 1")
})

test_that("geneticFixedEffects fails with phenoID error", {
    NrSNP <- 10
    geno <- simulateGenotypes(N=100, NrSNP=NrSNP, verbose=FALSE)
    expect_error(geneticFixedEffects(X_causal = geno$genotypes,
                                     P=10, N=100, phenoID=c("ID", 2),
                                     verbose=FALSE),
                 "phenoID has to be of length 1 and of type character")
})
test_that("noiseFixedEffects fails with out of range error",{
    pIndependentGenetic <- -0.5
    NrSNP <- 10
    geno <- simulateGenotypes(N=100, NrSNP=NrSNP, verbose=FALSE)
    expect_error(geneticFixedEffects(X_causal = geno$genotypes,
                             P=10, N=100, 
                             pIndependentGenetic = pIndependentGenetic, 
                             verbose=FALSE),
                 "Proportions have to be specified between")
})

context('Test noiseFixedEffects')
test_that("noiseFixedEffects fails with sample error", {
    expect_error(noiseFixedEffects(N=-10, P=10),
                 "N has/have to be greater than zero")
})

test_that("noiseFixedEffects fails with pheno error", {
    expect_error(noiseFixedEffects(N=10, P=-2),
                 "P has/have to be greater than zero")
})

test_that("noiseFixedEffects fails with integer error", {
    expect_error(noiseFixedEffects(N=10, P=2, NrConfounders=4.5),
                 "NrConfounders has/have to be integers, given")
})


test_that("noiseFixedEffects fails with missing category information",{
    expect_error(noiseFixedEffects(N=100, P=10, distConfounders ="cat_unif", 
                                   catConfounders = NULL),
                 "Confounder distribution set to ")
})

test_that("noiseFixedEffects fails with sampleID type error", {
    expect_error(noiseFixedEffects(N=10, P=2, NrConfounders=4, sampleID=2),
                 "sampleID has to be of length 1 and of type character")
})
test_that("noiseFixedEffects fails with sampleID length error", {
    expect_error(noiseFixedEffects(N=10, P=2, NrConfounders=4, 
                                   sampleID=c("p1", "p2")),
                "sampleID has to be of length 1 and of type character")
})

test_that("noiseFixedEffects fails with id_phenos length error", {
    expect_error(noiseFixedEffects(N=10, P=2, NrConfounders=4, id_phenos="p1"), 
                 "Length of id_phenos")
})

test_that("noiseFixedEffects fails with id_samples length error", {
    expect_error(noiseFixedEffects(N=10, P=2, NrConfounders=4, 
                                   id_samples=paste("ID", 1:15)), 
                 "Length of id_samples")
})

test_that("noiseFixedEffects fails with pTraitsAffects proportion error", {
    expect_error(noiseFixedEffects(N=10, P=2, NrConfounders=4, 
                                   pTraitsAffected=1.2), 
                 "Proportions have to be specified between 0 and 1")
})

test_that("noiseFixedEffects fails with NrFixedEffect number error", {
    expect_error(noiseFixedEffects(N=10, P=2, NrConfounders=4, 
                                   pTraitsAffected=0.8, NrFixedEffects=2), 
                 "NrFixedEffects specified to greater than 1")
})

test_that("noiseFixedEffects fails with parameter length mismatch error", {
    expect_error(noiseFixedEffects(N=10, P=2, NrConfounders=4, 
                                   pTraitsAffected=0.8, NrFixedEffects=3, 
                                   sdBeta=c(1,2)),
                 "Length of sdBeta")
})

tmp <- noiseFixedEffects(N=10000, P=2, NrConfounders=1, pTraitsAffected=0.8, 
                         NrFixedEffects=4, mBeta=1, sdBeta=c(2,1,1,3), 
                         catConfounders=3,  
                         distConfounders = c("cat_unif", "norm", "unif", "bin"), 
                         mConfounders=c(2,3), probConfounders=0.5)

test_that("noiseFixedEffects returns expected output for normal confounders", {
    expect_equal(mean(tmp$cov$sharedConfounder1_cat_unif1), 2, tolerance=0.05)
})
test_that("noiseFixedEffects returns expected output for normal confounders", {
    expect_equal(mean(tmp$cov$sharedConfounder2_norm1), 2, tolerance=0.05)
})
test_that("noiseFixedEffects returns expected output for normal confounders", {
    expect_equal(mean(tmp$cov$sharedConfounder3_unif1), 3, tolerance=0.05)
})
test_that("noiseFixedEffects returns expected output for binomial confounders", 
          {
    expect_equal(mean(tmp$cov$sharedConfounder4_bin1), 0.5, tolerance=0.05)
})


test_that("noiseFixedEffects fails with missing probabbility information",{
    expect_error(noiseFixedEffects(N=100, P=10, distConfounders ="bin", 
                                   probConfounders = NULL),
                 "Confounder distribution set to ")
})

test_that("noiseFixedEffects fails with length mismatch of confounder 
          information",{
    expect_error(noiseFixedEffects(N=100, P=10, 
                                   distConfounders =c("bin", "norm"),
                                   NrFixedEffects = 3),
                 "Length of distConfounders")
})

test_that("noiseFixedEffects fails with out of range error",{
    pIndependentConfounders <- 6
    expect_error(noiseFixedEffects(N=100, P=10, NrConfounders = 10,
                                pIndependentConfounders = 
                                    pIndependentConfounders),
                 "Proportions have to be specified between")
})

test_that("noiseFixedEffects returns correct number of independent effects",{
    pIndependentConfounders <- 0.6
    NrConfounders <- 10
    f <- noiseFixedEffects(N=100, P=10, NrConfounders = NrConfounders,
                           pIndependentConfounders = pIndependentConfounders)
    expect_equal(length(which(grepl("independent", colnames(f$cov_effect)))),
                 pIndependentConfounders * NrConfounders)
})

test_that("noiseFixedEffects fails with out of range error for probConfounders",
          {
    expect_error(noiseFixedEffects(N=100, P=10, NrConfounders = 10,
                                   distConfounders="bin", probConfounders=-0.5),
                 "Proportions have to be specified between 0 and 1")
})

test_that("noiseFixedEffects fails with categorical distribution error", {
    expect_error(noiseFixedEffects(N=10, P=2, NrConfounders=1, 
                                   catConfounders="small", 
                                   distConfounders = "cat_unif"),
    "Simulating categorical distribution: categories has to be numeric")
})

context('Test geneticBgEffects')
test_that('geneticBgEffects fails with sample input error', {
    X <- rnorm(10)
    kinship <- tcrossprod(X)
    colnames(kinship) <- paste("ID_", 1:10, sep="")
    expect_error(geneticBgEffects(P=10, N="all", kinship),
                 "N is/are not numeric")
})
test_that('geneticBgEffects fails with pheno input error', {
    X <- rnorm(10)
    kinship <- tcrossprod(X)
    diag(kinship) <- diag(kinship) + 1e-4
    colnames(kinship) <- paste("ID_", 1:10, sep="")
    expect_error(geneticBgEffects(P=2.5, N=10, kinship),
                 "P has/have to be integers")
})
test_that('geneticBgEffects fails with dimension error', {
    X <- rnorm(10)
    kinship <- tcrossprod(X)
    expect_error(geneticBgEffects(P=10, N=10, kinship[,-1]),
                 "Kinship matrix needs to be a square matrix")
})

test_that('geneticBgEffects fails with sample dimension error', {
    X <- rnorm(10)
    kinship <- tcrossprod(X)
    diag(kinship) <- diag(kinship) + 1e-4
    colnames(kinship) <- paste("ID_", 1:10, sep="")
    expect_error(geneticBgEffects(P=5, N=50, kinship, shared = TRUE,
                                  independent = FALSE),
                 "The number of samples specified")
})

test_that('geneticBgEffects fails with missing colnames', {
    X <- rnorm(10)
    kinship <- tcrossprod(X)
    diag(kinship) <- diag(kinship) + 1e-4
    expect_error(geneticBgEffects(P=5, N=10, kinship, shared = TRUE,
                                  independent = FALSE),
                 "Length of id_samples")
})

test_that('geneticBgEffects fails with id_samples dimension', {
    X <- rnorm(10)
    kinship <- tcrossprod(X)
    diag(kinship) <- diag(kinship) + 1e-4
    expect_error(geneticBgEffects(P=5, N=10, kinship, shared = TRUE,
                                  independent=FALSE,
                                  id_samples=paste("ID_", 1:20, sep="")),
                 "Length of id_samples")
})


test_that('geneticBgEffects fails with id_phenos dimension', {
    X <- rnorm(10)
    kinship <- tcrossprod(X)
    diag(kinship) <- diag(kinship) + 1e-4
    colnames(kinship) <- paste("ID_", 1:10, sep="")
    expect_error(geneticBgEffects(P=5, N=10, kinship, shared = TRUE,
                                  independent = FALSE,
                                  id_phenos=paste("P_", 1:2, sep="")),
                 "Length of id_phenos")
})

test_that('geneticBgEffects fails with phenoID type', {
    X <- rnorm(10)
    kinship <- tcrossprod(X)
    diag(kinship) <- diag(kinship) + 1e-4
    colnames(kinship) <- paste("ID_", 1:10, sep="")
    expect_error(geneticBgEffects(P=5, N=10, kinship, shared = TRUE,
                                  independent = FALSE,
                                  phenoID=5),
                 "phenoID has to be of length")
})
test_that('geneticBgEffects fails with positive semi-definite error', {
    X <- rnorm(10)
    kinship <- tcrossprod(X)
    diag(kinship) <- diag(kinship) + 1e-4
    kinship[1,] <- 0
    colnames(kinship) <- paste("ID_", 1:10, sep="")
    expect_error(geneticBgEffects(P=10, N=10, kinship),
                 "Kinship matrix is not positive semi-definite")
    
})

test_that('geneticBgEffects fails with input error', {
    X <- rnorm(10)
    kinship <- tcrossprod(X)
    diag(kinship) <- diag(kinship) + 1e-4
    colnames(kinship) <- paste("ID_", 1:10, sep="")
    expect_error(geneticBgEffects(P=10, N=10, kinship, shared = FALSE,
                                  independent = FALSE),
                 "At least one geneticBgEffect has to be specified")
})




test_that('geneticBgEffects returns correct dimensions', {
    NrSamples <- 10
    NrPheno <- 5
    X <- rnorm(NrSamples)
    kinship <- tcrossprod(X)
    diag(kinship) <- diag(kinship) + 1e-4
    colnames(kinship) <- paste("ID_", 1:10, sep="")
    b <- geneticBgEffects(P=NrPheno, N=NrSamples, kinship)
    expect_equal(dim(b$cov_shared), c(NrPheno, NrPheno))
    expect_equal(dim(b$cov_independent), c(NrPheno, NrPheno))
    expect_equal(dim(b$shared), c(NrSamples, NrPheno))
    expect_equal(dim(b$independent), c(NrSamples, NrPheno))
})

context('Test noiseBgEffects')
test_that('noiseBgEffects fails with sample input error', {
    expect_error(noiseBgEffects(P=10, N="cohort"),
                 "N is/are not numeric")
})
test_that('noiseBgEffects fails with pheno input error', {
    expect_error(geneticBgEffects(P=-3, N=10),
                 "P has/have to be greater")
})
test_that('noiseBgEffects fails with id_samples dimension', {
    expect_error(noiseBgEffects(P=5, N=10, shared = TRUE,
                                  independent=FALSE,
                                  id_samples=paste("ID_", 1:20, sep="")),
                 "Length of id_samples")
})


test_that('noiseBgEffects fails with id_phenos dimension', {
    expect_error(noiseBgEffects(P=5, N=10, shared = TRUE,
                                  independent = FALSE,
                                  id_phenos=paste("P_", 1:7, sep="")),
                 "Length of id_phenos")
})

test_that('noiseBgEffects fails with phenoID type', {
    expect_error(noiseBgEffects(P=5, N=10, shared = TRUE,
                                  independent = FALSE,
                                  phenoID=5),
                 "phenoID has to be of length")
})

test_that('noiseBgEffects fails with sampleID type', {
    expect_error(noiseBgEffects(P=5, N=10, shared = TRUE,
                                independent = FALSE,
                                sampleID=paste("ID", 1:10, sep="")),
                 "sampleID has to be of length")
})
test_that('noiseBgEffects fails with input error', {
    expect_error(noiseBgEffects(P=10, N=10, shared = FALSE,
                                  independent = FALSE),
                 "At least one noiseBgEffect has to be specified")
    
})

test_that('noiseBgEffects returns correct dimensions', {
    NrSamples <- 20
    NrPheno <- 7
    b <- noiseBgEffects(P=NrPheno, N=NrSamples)
    expect_equal(dim(b$cov_shared), c(NrPheno, NrPheno))
    expect_equal(dim(b$cov_independent), c(NrPheno, NrPheno))
    expect_equal(dim(b$shared), c(NrSamples, NrPheno))
    expect_equal(dim(b$independent), c(NrSamples, NrPheno))
})

context('Test correlatedBgEffects')
test_that('correlatedBgEffects fails with sample input error', {
    expect_error(correlatedBgEffects(P=10, N="all", pcorr=0.9),
                 "N is/are not numeric")
})
test_that('correlatedBgEffects fails with pheno input error', {
    expect_error(correlatedBgEffects(P=0, N=10, pcorr=0.9),
                 "P has/have to be greater")
})
test_that('correlatedBgEffects fails with id_samples dimension', {
    expect_error(correlatedBgEffects(P=5, N=10, pcorr=0.8,
                                id_samples=paste("ID_", 1:20, sep="")),
                 "Length of id_samples")
})


test_that('correlatedBgEffects fails with id_phenos dimension', {
    expect_error(correlatedBgEffects(P=5, N=10, pcorr=0.7,
                                id_phenos=paste("P_", 1:7, sep="")),
                 "Length of id_phenos")
})

test_that('correlatedBgEffects fails with phenoID type', {
    expect_error(correlatedBgEffects(P=5, N=10, pcorr=-0.8,
                                phenoID=5),
                 "phenoID has to be of length")
})

test_that('correlatedBgEffects fails with sampleID type', {
    expect_error(correlatedBgEffects(P=5, N=10, pcorr=0.1,
                                sampleID=paste("ID", 1:10, sep="")),
                 "sampleID has to be of length")
})

test_that('correlatedBgEffects fails with pcorr range error', {
    expect_error(correlatedBgEffects(P=5, N=10, pcorr=1),
                 "pcorr has to be greater than zero and less than 1")
    expect_error(correlatedBgEffects(P=5, N=10, pcorr=0),
                 "pcorr has to be greater than zero and less than 1")
})

test_that("correlatedBgEffects fails with not square matrix error",{
    X <- rnorm(20)
    corr_mat <- tcrossprod(X)
    expect_error(correlatedBgEffects(N=100, P=20, corr_mat = corr_mat[,-1]),
                 "Correlation matrix needs to be a square matrix")
})

test_that("correlatedBgEffects fails with dimension mismatch error",{
    X <- rnorm(20)
    corr_mat <- tcrossprod(X)
    expect_error(correlatedBgEffects(N=100, P=10, corr_mat = corr_mat),
                 "Dimensions of correlation matrix")
})
      
test_that("correlatedBgEffects returns correct dimensions",{
    b <- correlatedBgEffects(N=100, P=10, pcorr=0.6)
    expect_equal(dim(b$correlatedBg), c(100,10))
    expect_equal(dim(b$cov_correlated), c(10,10))
})           