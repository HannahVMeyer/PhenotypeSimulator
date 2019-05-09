context('Test genotype functions')

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
                 paste("getAllelelFrequencies .*", sep=""))
})

test_that('standardiseGenotypes returns output matrix/dataframe of same 
          dimension as input matrix/dataframe', {
    geno <- cbind(rbinom(100, 2, 0.3), rbinom(100, 2, 0.4))
    geno_sd <- standardiseGenotypes(geno)
    expect_equal(dim(geno), dim(geno_sd))
})

test_that('simulateGenotypes fails with range error', {
    expect_error(simulateGenotypes(N=10, NrSNP=10, freq=c(-0.1, 0.4)), 
                 paste("Allele frequencies must be between 0 and 1"))
})

test_that('simulateGenotypes fails with snpID error', {
    expect_error(simulateGenotypes(N=10, NrSNP=10, snpID=2, freq=c(0.1, 0.4)), 
                 "snpID has to be of length 1 and of type character")
})

test_that('simulateGenotypes fails with sampleID error', {
    expect_error(simulateGenotypes(N=10, NrSNP=10, sampleID=c("ID1", "ID2"), 
                                   freq=c(0.1, 0.4)), 
                 "sampleID has to be of length 1 and of type character")
})

test_that('simulateGenotypes simulates correct frequencies', {
    x <- simulateGenotypes(N=1000, NrSNP=1, freq=0.4)
    f <- min(getAlleleFrequencies(x$genotypes))
    expect_equal(f, 0.4, tolerance = 0.1, scale = 0.4) 
})

test_that('readStandardGenotypes fails with missing format error', {
    expect_error(readStandardGenotypes(N=100, filename="test"),
                 "Genotypefile format has to be specified, supported")
})

test_that('readStandardGenotypes fails with format error', {
    expect_error(readStandardGenotypes(N=100, filename="test", format="GCTA"),
                 ".*is not a supported genotype format.*")
})

test_that('readStandardGenotypes returns correctly formated genotypes GENOME', {
    filename_genome  <- system.file("extdata/genotypes/genome/",
                                    "genotypes_genome.txt", 
                                    package = "PhenotypeSimulator") 
    data_genome <- readStandardGenotypes(N=100, filename_genome, 
                                         format ="genome")
    expect_equal(dim(data_genome$genotypes)[1], length(data_genome$id_samples))
    expect_equal(dim(data_genome$genotypes)[2], length(data_genome$id_snps))
    expect_true(is.null(data_genome$format_files))
    expect_true(all(data_genome$genotypes %in% c(0,1,2)))
})


test_that('readStandardGenotypes returns correctly formated genotypes plink', {
    filename_plink  <- system.file("extdata/genotypes/plink/",
                                   "genotypes_plink.bed", 
                                   package = "PhenotypeSimulator") 
    filename_plink <- gsub("\\.bed", "", filename_plink)
    data_plink <- readStandardGenotypes(N=100, filename=filename_plink, 
                                        format="plink")
    expect_equal(dim(data_plink$genotypes)[1], length(data_plink$id_samples))
    expect_equal(dim(data_plink$genotypes)[2], length(data_plink$id_snps))
    expect_identical(data_plink$id_snps[1:3], c("SNP_0", "SNP_1","SNP_2"))
    expect_identical(data_plink$id_samples[1:3], c("per0", "per1","per2"))
    expect_false(is.null(data_plink$format_files))
    expect_true(all(data_plink$genotypes %in% c(0,1,2)))
})


test_that('readStandardGenotypes returns correctly formated genotypes HAPGEN2', {
    filename_hapgen  <- system.file("extdata/genotypes/hapgen/",
                                    "genotypes_hapgen.controls.gen",
                                    package = "PhenotypeSimulator") 
    filename_hapgen <- gsub("\\.gen", "", filename_hapgen)
    data_hapgen <- readStandardGenotypes(N=100, filename_hapgen, format='oxgen')
    expect_equal(dim(data_hapgen$genotypes)[1], length(data_hapgen$id_samples))
    expect_equal(dim(data_hapgen$genotypes)[2], length(data_hapgen$id_snps))
    expect_identical(data_hapgen$id_snps[1:2], c("rs10399749", "rs4030303"))
    expect_identical(data_hapgen$id_samples[1:3], c("id1_0", "id1_1","id1_2"))
    expect_false(is.null(data_hapgen$format_files))
    expect_true(all(data_hapgen$genotypes %in% c(0,1,2)))
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

test_that('getCausalSNPs fails with genotype input error', {
    geno <- simulateGenotypes(N=10, NrSNP=10)
    expect_error(getCausalSNPs(N=10, NrCausalSNPs=10, genotypes=geno),
                 "Genotypes have to be provided as")
})

test_that('getCausalSNPs fails with genotype input error', {
    geno <- simulateGenotypes(N=10, NrSNP=10)
    expect_error(getCausalSNPs(N=10, NrCausalSNPs=50, genotypes=geno$genotypes),
                 "Number of genotypes is less than number of causal SNPs")
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

test_that('getCausalSNPs fails with SNP number error', {
    genotypeFile <- system.file("extdata/genotypes/",
                                "genotypes_chr22.csv",
                                package = "PhenotypeSimulator")
    genoFilePrefix <- gsub("chr.*", "", genotypeFile)
    genoFileSuffix <- ".csv"
    expect_error(getCausalSNPs(N=50, NrCausalSNPs=2000, chr=22, 
                               genoFilePrefix=genoFilePrefix, 
                               genoFileSuffix=genoFileSuffix),
                 "Number of causal SNPs to be chosen from chromosome")
})

test_that('getCausalSNPs fails with delimiter error', {
    genotypeFile <- system.file("extdata/genotypes/",
                                "genotypes_chr22.csv",
                                package = "PhenotypeSimulator")
    genoFilePrefix <- gsub("chr.*", "", genotypeFile)
    genoFileSuffix <- ".csv"
    expect_error(getCausalSNPs(N=50, NrCausalSNPs=10, chr=22, 
                               genoFilePrefix=genoFilePrefix, 
                               genoFileSuffix=genoFileSuffix,
                               delim="\t", verbose=FALSE),
                 "Delimiter specified for genoFilePrefix-genoFileSuffix")
})

test_that('getCausalSNPs samples from oxgen formated genotypes', {
    genotypeFile  <- system.file("extdata/genotypes/hapgen/",
                                    "genotypes_hapgen_chr22.gen",
                                    package = "PhenotypeSimulator") 
    genoFilePrefix <- gsub("chr.*", "", genotypeFile)
    genoFileSuffix <- ".gen"
    tmp <- getCausalSNPs(N=50, NrCausalSNPs=10, chr=22, 
                               genoFilePrefix=genoFilePrefix, 
                               genoFileSuffix=genoFileSuffix,
                               format='oxgen', verbose=FALSE)
    expect_equal(dim(tmp)[1], 50)
    expect_equal(dim(tmp)[2], 10)

})

test_that('getCausalSNPs samples from delim formated genotypes', {
    genotypeFile <- system.file("extdata/genotypes/",
                                "genotypes_chr22.csv",
                                package = "PhenotypeSimulator")
    genoFilePrefix <- gsub("chr.*", "", genotypeFile)
    genoFileSuffix <- ".csv"
    tmp <- getCausalSNPs(N=50, NrCausalSNPs=10, chr=22, 
                         genoFilePrefix=genoFilePrefix, 
                         genoFileSuffix=genoFileSuffix,
                         format='delim', delim=",", verbose=FALSE)
    expect_equal(dim(tmp)[1], 50)
    expect_equal(dim(tmp)[2], 10)
})

test_that('getCausalSNPs samples from bimbam formated genotypes', {
    genotypeFile <- system.file("extdata/genotypes/",
                                "genotypes_chr22.bimbam",
                                package = "PhenotypeSimulator")
    genoFilePrefix <- gsub("chr.*", "", genotypeFile)
    genoFileSuffix <- ".bimbam"
    tmp <- getCausalSNPs(N=50, NrCausalSNPs=10, chr=22, 
                         genoFilePrefix=genoFilePrefix, 
                         genoFileSuffix=genoFileSuffix,
                         format='bimbam', delim=",", verbose=FALSE)
    expect_equal(dim(tmp)[1], 50)
    expect_equal(dim(tmp)[2], 10)
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
              expect_error(getKinship(N=100),
                           "Either X or kinshipfile must be provided")
          })

test_that('getKinship fails with non-existing kinship file',{
    expect_error(getKinship(N=100, kinshipfile="~/test/kinship"), 
                 "does not exist")
})


test_that('getKinship returns positive semi-definite matrix',{
    geno <- simulateGenotypes(N=100, NrSNP=10, verbose=FALSE)
    kin <- getKinship(N=100, X=geno$genotypes, verbose=FALSE)
    kin_eigen <- as.complex(eigen(kin, only.values = TRUE)$values)
    expect_true(all(Re(kin_eigen) > 0))
})

test_that('getKinship fails with sample error', {
    kinshipfile <- system.file("extdata/kinship", "kinship.csv",
                               package = "PhenotypeSimulator")
    expect_error(getKinship(N=1000, kinshipfile=kinshipfile, verbose=FALSE),
                 "Number of samples specifid is greater than number")
})

test_that('getKinship fails with sample ID error', {
    kinshipfile <- system.file("extdata/kinship", "kinship.csv",
                               package = "PhenotypeSimulator")
    expect_error(getKinship(N=10, kinshipfile=kinshipfile, 
                            id_samples = paste("Sample_", 1:10, sep=""),
                            verbose=FALSE),
                 "Not all id_samples are present")
})

test_that('getKinship fails with sample ID to sample mismatch error', {
    kinshipfile <- system.file("extdata/kinship", "kinship.csv",
                               package = "PhenotypeSimulator")
    expect_error(getKinship(N=10, kinshipfile=kinshipfile, 
                            id_samples = paste("Sample_", 1:12, sep=""),
                            verbose=FALSE),
                 "Length of id_samples")
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

