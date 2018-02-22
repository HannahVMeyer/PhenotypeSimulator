context('Test createphenotype functions')

noiseBg <- noiseBgEffects(N=100, P=10)
noiseFixed <- noiseFixedEffects(N=100, P=10)
noiseCorrelated <- correlatedBgEffects(N=100, P=10, pcorr=0.6)
genotypes <- simulateGenotypes(N=100, NrSNP=20, verbose=FALSE)
kinship <- getKinship(N=100, X=genotypes$genotypes, standardise=TRUE, 
                      verbose=FALSE)
causalSNPs <- getCausalSNPs(N=100, genotypes=genotypes$genotypes, verbose=FALSE)
geneticFixed <- geneticFixedEffects(X_causal=causalSNPs, P=10, N=100, 
                                    verbose=FALSE)
geneticBg <- geneticBgEffects(P=10, N=100, kinship=kinship)

context("Tests rescaleVariance")
test_that("rescaleVariance fails when propvar is not numeric", {
    expect_error(rescaleVariance(noiseBg$shared, propvar = 'a'),
                 "propvar needs to be numeric")
})
test_that("rescaleVariance fails when propvar is out of range", {
    expect_error(rescaleVariance(noiseBg$shared, propvar = 1.1),
                 "propvar cannot be less than 0 and or greater than 1")
})
test_that("rescaleVariance fails when component is not a matrix", {
    expect_error(rescaleVariance(data.frame(noiseBg$shared), propvar = 0.1),
                 "component needs to be a matrix")
})
test_that("rescaleVariance returns NULL when no propvar is specified", {
    expect_equal(rescaleVariance(noiseBg$shared, propvar = 0), NULL)
})

test_that("rescaleVariance returns correclty scaled average column variance", {
    tmp <- rescaleVariance(noiseBg$shared, propvar = 0.1)
    expect_equal(mean(diag(var(tmp$component))), 0.1, tolerance=5e-5)
})


context("Tests setModel")
test_that("setModel fails with input error", {
    expect_error(setModel(), "No variance components specified")
})

test_that("setModel fails with out of range error", {
    expect_error(setModel(genVar=0.4, h2s=1.2, theta=0.8,  eta=0.8, 
                          delta=0.2, gamma=0.8, rho=0.2, phi=0.6, 
                          alpha=0.8, pcorr=0.6, pIndependentConfounders=0.4,  
                          pTraitIndependentConfounders=0.2,  
                          pIndependentGenetic=0.4, 
                          pTraitIndependentGenetic=0.2, verbose=FALSE),
                 "Proportions have to be specified between 0 and 1:")
})
test_that("setModel fails with missing variance error", {
    expect_error(setModel(h2s=1, theta=0.8,  eta=0.8, 
                          delta=0.2, gamma=0.8, rho=0.2, phi=0.6, 
                          alpha=0.8, pcorr=0.6, pIndependentConfounders=0.4,  
                          pTraitIndependentConfounders=0.2,  
                          pIndependentGenetic=0.4, 
                          pTraitIndependentGenetic=0.2, verbose=FALSE),
                 "Neither genVar nor noiseVar are provided, thus")
})

test_that("setModel fails with sum of variance error", {
    expect_error(setModel(genVar=0.2, noiseVar=0.9, h2s=1, theta=0.8,  eta=0.8, 
                          delta=0.2, gamma=0.8, rho=0.2, phi=0.6, 
                          alpha=0.8, pcorr=0.6, pIndependentConfounders=0.4,  
                          pTraitIndependentConfounders=0.2,  
                          pIndependentGenetic=0.4, 
                          pTraitIndependentGenetic=0.2, verbose=FALSE),
                 "Sum of genetic variance and noise variance is greater than 1")
})

test_that("setModel fails with extra noise component error", {
    expect_error(setModel(noiseVar=0, h2s=1, theta=0.8,  eta=0.8, 
                          delta=0.2, gamma=0.8, rho=0.2, phi=0.6, 
                          alpha=0.8, pcorr=0.6, pIndependentConfounders=0.4,  
                          pTraitIndependentConfounders=0.2,  
                          pIndependentGenetic=0.4, 
                          pTraitIndependentGenetic=0.2, verbose=FALSE),
                 "The noise variance is set to 0 ")
})

test_that("setModel fails with extra genetic component error", {
    expect_error(setModel(genVar=0, h2s=1, theta=0.8,  eta=0.8, 
                          delta=0.2, gamma=0.8, rho=0.2, phi=0.6, 
                          alpha=0.8, pcorr=0.6, pIndependentConfounders=0.4,  
                          pTraitIndependentConfounders=0.2,  
                          pIndependentGenetic=0.4, 
                          pTraitIndependentGenetic=0.2, verbose=FALSE),
                 "The genetic variance is set to 0 ")
})

test_that("setModel fails with missing noise component error", {
    expect_error(setModel(noiseVar=0.8, h2s=1, theta=0.8,  eta=0.8, 
                          pcorr=0.6, pIndependentConfounders=0.4,  
                          pTraitIndependentConfounders=0.2,  
                          pIndependentGenetic=0.4, 
                          pTraitIndependentGenetic=0.2,
                          verbose=FALSE),
                 "Neither delta ")
})

test_that("setModel fails with proportions larger than 1 error", {
    expect_error(setModel(noiseVar=0.8, h2s=1, theta=0.8, delta=0.8, phi=0.2, 
                          rho=0.2, eta=0.8, 
                          pcorr=0.6, pIndependentConfounders=0.4,  
                          pTraitIndependentConfounders=0.2,  
                          pIndependentGenetic=0.4, 
                          pTraitIndependentGenetic=0.2,
                          verbose=FALSE),
                 "Sum of the proportion of the variance of noise")
})

test_that("setModel fails with proportions less than 1 error", {
    expect_error(setModel(noiseVar=0.8, h2s=1, theta=0.8, delta=0.2, phi=0.2, 
                          rho=0.2, eta=0.8, 
                          pcorr=0.6, pIndependentConfounders=0.4,  
                          pTraitIndependentConfounders=0.2,  
                          pIndependentGenetic=0.4, 
                          pTraitIndependentGenetic=0.2,
                          verbose=FALSE),
                 "Sum of the proportion of the variance of noise")
})

test_that("setModel fails with insufficient noise component error", {
    expect_error(setModel(noiseVar=0.8, h2s=1, theta=0.8,  eta=0.8, 
                          delta=0.8, pcorr=0.6, pIndependentConfounders=0.4,  
                          pTraitIndependentConfounders=0.2,  
                          pIndependentGenetic=0.4, 
                          pTraitIndependentGenetic=0.2,
                          verbose=FALSE),
                 "Not enough components provided to set proportions of")
})

test_that("setModel fails with insufficient genetic component error", {
    expect_error(setModel(genVar=0.8, theta=0.8,  eta=0.8, 
                          delta=0.2, gamma=0.8, rho=0.2, phi=0.6, 
                          alpha=0.8, pcorr=0.6, pIndependentConfounders=0.4,  
                          pTraitIndependentConfounders=0.2,  
                          pIndependentGenetic=0.4, 
                          pTraitIndependentGenetic=0.2, verbose=FALSE),
                 "Neither genetic variant effect variance")
})

test_that("setModel fails with sum of genetic component error", {
    expect_error(setModel(genVar=0.8, h2s=0.4, h2bg=0.3, theta=0.8,  eta=0.8, 
                          delta=0.2, gamma=0.8, rho=0.2, phi=0.6, 
                          alpha=0.8, pcorr=0.6, pIndependentConfounders=0.4,  
                          pTraitIndependentConfounders=0.2,  
                          pIndependentGenetic=0.4, 
                          pTraitIndependentGenetic=0.2, verbose=FALSE),
                 "Sum of the proportion of the variance of genetic")
})

context("Tests runSimulation")
test_that("runSimulation returns correct output", {
    expect_equal(length(runSimulation(N=10, P=5, genVar=1, h2s=0.2)), 5)
})