context('Test variance component functions')

context('Test geneticFixedEffects')
test_that('geneticFixedEffects throws sample error',{
    geno <- simulateGenotypes(N=10, NrSNP=100, verbose=FALSE)
    expect_error(geneticFixedEffects(X_causal = geno$genotypes,
                                     P=10, N=20, verbose=FALSE),
                 "Number of samples in SNP matrix")
})

test_that('geneticFixedEffects throws sample error',{
    geno <- simulateGenotypes(N=10, NrSNP=100, verbose=FALSE)
    expect_error(geneticFixedEffects(X_causal = geno$genotypes,
                                     P=10, N=20, verbose=FALSE),
                 "Number of samples in SNP matrix")
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
test_that("noiseFixedEffects fails with missing category information",{
    expect_error(noiseFixedEffects(N=100, P=10, distConfounders ="cat_unif", 
                                   catConfounders = NULL),
                 "Confounder distribution set to ")
    
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

context('Test geneticBgEffects')
test_that('geneticBgEffects fails with dimension error', {
    X <- rnorm(10)
    kinship <- tcrossprod(X)
    expect_error(geneticBgEffects(P=10, N=10, kinship[,-1]),
                 "Kinship matrix needs to be a square matrix")
                 
})
test_that('geneticBgEffects fails with positive semi-definite error', {
    X <- rnorm(10)
    kinship <- tcrossprod(X)
    diag(kinship) <- diag(kinship) + 1e-4
    kinship[1,] <- 0
    expect_error(geneticBgEffects(P=10, N=10, kinship),
                 "Kinship matrix is not positive semi-definite")
    
})

test_that('geneticBgEffects fails with input error', {
    X <- rnorm(10)
    kinship <- tcrossprod(X)
    diag(kinship) <- diag(kinship) + 1e-4
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
    b <- geneticBgEffects(P=NrPheno, N=NrSamples, kinship)
    expect_equal(dim(b$cov_shared), c(NrPheno, NrPheno))
    expect_equal(dim(b$cov_independent), c(NrPheno, NrPheno))
    expect_equal(dim(b$shared), c(NrSamples, NrPheno))
    expect_equal(dim(b$independent), c(NrSamples, NrPheno))
})

context('Test noiseBgEffects')
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