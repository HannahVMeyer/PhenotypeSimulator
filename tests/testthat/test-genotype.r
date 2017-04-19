context('Test genotype functions')
test_that('standardiseGenotypes returns output matrix/dataframe of same 
          dimension as input matrix/dataframe', {
    geno <- cbind(rbinom(100, 2, 0.3), rbinom(100, 2, 0.4))
    geno_sd <- standardiseGenotypes(geno)
    expect_equal(dim(geno), dim(geno_sd))
})


test_that('getAlleleFrequencies returns a vector of allele frequencies equal to 
          the simulated frequencies', {
    snp <- rbinom(1000, 2, 0.3)
    allelefreq <- getAlleleFrequencies(snp)
    expect_equal(allelefreq[1], 0.3, tolerance = 0.1, scale = 0.3) 
})

test_that('getKinshio fails if neither kinshipfile nor genotype matrix are 
          provided',{
    expect_error(getKinship(), "Either X or kinshipfile must be provided")
})

