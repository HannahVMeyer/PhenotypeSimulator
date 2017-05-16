context('Test phenotype functions')

noiseBg <- noiseBgEffects(N=100, P=10)
noiseFixed <- noiseFixedEffects(N=100, P=10)
noiseCorrelated <- correlatedBgEffects(N=100, P=10, pcorr=0.6)
genotypes <- simulateGenotypes(N=100, NrSNP=20, verbose=FALSE)
kinship <- getKinship(genotypes$X_sd, norm=TRUE, verbose=FALSE)
causalSNPs <- getCausalSNPs(genotypes=genotypes)
geneticFixed <- geneticFixedEffects(X_causal=causalSNPs, P=10, N=100)
geneticBg <- geneticBgEffects(P=10, kinship=kinship)

test_that('getAlleleFrequencies returns a vector of allele frequencies equal to 
          the simulated frequencies', {
    snp <- rbinom(1000, 2, 0.3)
    allelefreq <- getAlleleFrequencies(snp)
    expect_equal(allelefreq[1], 0.3, tolerance = 0.1, scale = 0.3) 
})

test_that('createPheno fails if no phenotype components are provided',{
    expect_error(createPheno(P=10, N=100), 
                 "No phenotype components provided, at least one is required")
})

test_that('createPheno computes variance components such that the sum of the 
          components adds to 1', {
    phenotype <- createPheno(P=10, N=100, noiseBg=noiseBg, 
                             genFixed=geneticFixed, genBg=geneticBg, 
                             noiseFixed=noiseFixed, 
                             correlatedBg=noiseCorrelated, 
                             modelNoise="noiseFixedAndBgAndCorrelated", 
                             modelGenetic="geneticFixedAndBg", 
                             genVar=0.4, h2s=0.025, delta=0.4, rho=0.1)
    expect_equal( sum(phenotype$varComponents
                      [grepl("Var", names(phenotype$varComponents))]), 1)
    expect_equal( sum(phenotype$varComponents
                      [grepl("h2", names(phenotype$varComponents))]), 1)
    expect_equal( sum(phenotype$varComponents
                      [grepl("var_gen", names(phenotype$varComponents))]),  
                  as.numeric(phenotype$varComponents
                             [grepl("genVar", names(phenotype$varComponents))]))
    expect_equal( sum(phenotype$varComponents
                      [grepl("var_noise", names(phenotype$varComponents))]),  
                  as.numeric(phenotype$varComponents
                      [grepl("noiseVar", names(phenotype$varComponents))]))
})

