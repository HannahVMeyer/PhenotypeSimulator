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
    expect_equal(allelefreq[2], 0.3, tolerance = 0.1, scale = 0.3) 
})

test_that('getAlleleFrequencies fails with encoding error', {
              snp <- rep(c(1,2,3), 4)
              expect_error(getAlleleFrequencies(snp), 
                          "SNP vector contains alleles not encoded as 0, 1 or 2"
              )
})

test_that('getAlleleFrequencies fails with dimension error', {
    snp <- matrix(rep(c(1,2,3), 4), nrow=2)
    expect_error(getAlleleFrequencies(snp), 
                 paste("getAllelelFrequencies .*", dim(snp)[1], " x ",
                       dim(snp)[2], "]", sep=""))
})

test_that('simulateGenotypes fails with range error', {
    expect_error(simulateGenotypes(N=10, NrSNP=10, freq=c(-0.1, 0.4)), 
                 paste("Allele frequencies must be between 0 and 1"))
})

test_that('simulateGenotypes simulates correct frequencies', {
    x <- simulateGenotypes(N=1000, NrSNP=1, freq=0.4)
    f <- min(getAlleleFrequencies(x$genotypes))
    expect_equal(f, 0.4, tolerance = 0.1, scale = 0.4) 
})

test_that('readStandardGenotypes fails with format error', {
    expect_error(readStandardGenotypes(N=100, filename="test", format="GCTA"),
                 ".*is not a supported genotype format.*")
})

test_that('readStandardGenotypes fails with sample error', {
    filename_plink  <- system.file("extdata/genotypes/plink/",
                                   "genotypes_plink.bed",
                                   package = "PhenotypeSimulator") 
    filename_plink <- gsub("\\.bed", "", filename_plink)
    expect_error(readStandardGenotypes(N=10000, filename=filename_plink, 
                                       format="plink"),
                 "Sample number specified.*")
})

test_that('readStandardGenotypes subsampling works', {
    filename_plink  <- system.file("extdata/genotypes/plink/",
                                   "genotypes_plink.bed",
                                   package = "PhenotypeSimulator") 
    filename_plink <- gsub("\\.bed", "", filename_plink)
    geno <- readStandardGenotypes(N=50, filename=filename_plink, 
                                       format="plink")
    expect_equal(nrow(geno$genotypes), 50)
    expect_equal(length(geno$id_samples), 50)
})


test_that('getCausalSNPs fails with missing input error', {
    expect_error(getCausalSNPs(N=100, NrCausalSNPs = 20, verbose=FALSE),
                 "No genotype information provided,.*")
})

test_that('getCausalSNPs fails with incorrect input error', {
    expect_error(getCausalSNPs(N=100, NrCausalSNPs = 20,
                               genoFilePrefix = "~/test", verbose=FALSE),
                 "genoFilePrefix contains ~.*")
})

test_that('getCausalSNPs fails with missing chromosome error', {
    expect_error(getCausalSNPs(N=100, NrCausalSNPs = 20,
                               genoFilePrefix = "/test", verbose=FALSE),
                 "No information about chromosomes to.*")
})

test_that('getCausalSNPs fails with too much chromosome information error', {
    expect_error(getCausalSNPs(N=100, NrCausalSNPs = 20,
                               genoFilePrefix = "/test", NrChrCausal = 20,
                               chr=10, verbose=FALSE),
                 "Too much information for sampling .*")
})
test_that('getCausalSNPs fails with wrong chromosome information error', {
    expect_error(getCausalSNPs(N=100, NrCausalSNPs = 20,
                               genoFilePrefix = "/test", NrChrCausal = 3,
                               NrSNPsOnChromosome = c(4,5), verbose=FALSE),
                 "Not enough information about numbers .*")
})

test_that('getCausalSNPs fails with wrong file error', {
    expect_error(getCausalSNPs(N=100, NrCausalSNPs = 20,
                               genoFilePrefix = "/test", 
                               genoFileSuffix = ".csv", 
                               NrChrCausal = 3,
                               NrSNPsOnChromosome = c(4,5,6), verbose=FALSE),
                 "does not exist.*")
})

test_that('getCausalSNPs fails with SNP error', {
    geno <- simulateGenotypes(N=100, NrSNP=10, verbose=FALSE)
    expect_error(getCausalSNPs(N=100, NrCausalSNPs = 20, 
                               genotypes = geno$genotypes,
                               verbose=FALSE),
                 "Number of genotypes is less than number of causal SNPs")
})

test_that('getCausalSNPs fails with sample error', {
    geno <- simulateGenotypes(N=100, NrSNP=10, verbose=FALSE)
    expect_error(getCausalSNPs(N=1000, NrCausalSNPs = 10, 
                               genotypes = geno$genotypes,
                               verbose=FALSE),
                 "Sample number specified exceeds number.*")
})

test_that('getCausalSNPs returns correct output dimensions', {
    geno <- simulateGenotypes(N=100, NrSNP=10, verbose=FALSE)
    expect_equal(dim(getCausalSNPs(N=100, NrCausalSNPs = 10, 
                               genotypes = geno$genotypes,
                               verbose=FALSE)),
                 c(100,10))
})

test_that('getKinship fails if neither kinshipfile nor genotype matrix are 
          provided',{
    expect_error(getKinship(N=100), "Either X or kinshipfile must be provided")
})

test_that('getKinship returns positive semi-definite matrix',{
    geno <- simulateGenotypes(N=100, NrSNP=10, verbose=FALSE)
    kin <- getKinship(X=geno$genotypes, verbose=FALSE)
    kin_eigen <- as.complex(eigen(kin, only.values = TRUE)$values)
    expect_true(all(Re(kin_eigen) > 0))
})

test_that('getKinship fails with sample error', {
    kinshipfile <- system.file("extdata/kinship", "kinship.csv",
                               package = "PhenotypeSimulator")
    expect_error(getKinship(N=1000, kinshipfile=kinshipfile, verbose=FALSE),
                 "Number of samples specifid is greater than number")
})

test_that('getKinship returns kinship of correct dimensions', {
    kinshipfile <- system.file("extdata/kinship", "kinship.csv",
                               package = "PhenotypeSimulator")
    K <- getKinship(N=50, kinshipfile=kinshipfile, verbose=FALSE)
    expect_equal(dim(K), c(50,50))
})


test_that('getKinship sampling by sample ID works', {
    kinshipfile <- system.file("extdata/kinship", "kinship.csv",
                               package = "PhenotypeSimulator")
    K <- getKinship(N=10, kinshipfile=kinshipfile, 
                    id_samples = paste("ID_", 1:10, sep=""),
                    verbose=FALSE)
    expect_equal(dim(K), c(10,10))
})

test_that('getKinship fails with sample ID error', {
    kinshipfile <- system.file("extdata/kinship", "kinship.csv",
                               package = "PhenotypeSimulator")
    expect_error(getKinship(N=10, kinshipfile=kinshipfile, 
                    id_samples = paste("Sample_", 1:10, sep=""),
                    verbose=FALSE),
                 "Not all id_samples are present")
})