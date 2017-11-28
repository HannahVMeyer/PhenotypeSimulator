context('Test phenotype functions')

noiseBg <- noiseBgEffects(N=100, P=10)
noiseFixed <- noiseFixedEffects(N=100, P=10)
noiseCorrelated <- correlatedBgEffects(N=100, P=10, pcorr=0.6)
genotypes <- simulateGenotypes(N=100, NrSNP=20, verbose=FALSE)
kinship <- getKinship(genotypes$genotypes, standardise=TRUE, verbose=FALSE)
causalSNPs <- getCausalSNPs(genotypes=genotypes$genotypes)
geneticFixed <- geneticFixedEffects(X_causal=causalSNPs, P=10, N=100)
geneticBg <- geneticBgEffects(P=10, kinship=kinship)

test_that('getAlleleFrequencies returns a vector of allele frequencies equal to 
          the simulated frequencies', {
    snp <- rbinom(1000, 2, 0.3)
    allelefreq <- getAlleleFrequencies(snp)
    expect_equal(allelefreq[2], 0.3, tolerance = 0.1, scale = 0.3) 
})




