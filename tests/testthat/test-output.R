context('Test output functions')

tmp <- './tmpdata'
simulation <- runSimulation(N=10, P=2, genVar=0.4, h2s=0.2, phi=0.8, delta=0.2)
genotypes <- simulation$rawComponents$genotypes
kinship <-  simulation$rawComponents$kinship
phenotypes <- simulation$phenoComponentsFinal$Y
covariates <- simulation$phenoComponentsIntermediate$noiseFixed$cov

context('Test writeStandardOutput')
test_that('writeStandardOutput fails with missing output format', {
    expect_error(writeStandardOutput(directory=tempdir(), 
                                     genotypes=genotypes$genotypes, 
                                     phenotypes=phenotypes, 
                                     id_samples = genotypes$id_samples, 
                                     id_snps = genotypes$id_snps, 
                                     id_phenos = colnames(phenotypes), 
                                     kinship=kinship),
                 "Output format has to be specified") 
          })

test_that('writeStandardOutput fails with path error', {
    expect_error(writeStandardOutput(directory="~/test", 
                                     genotypes=genotypes$genotypes, 
                                     phenotypes=phenotypes, 
                                     id_samples = genotypes$id_samples, 
                                     id_snps = genotypes$id_snps, 
                                     id_phenos = colnames(phenotypes), 
                                     kinship=kinship, 
                                     format="plink"),
                 "directory contains ~:") 
})

test_that('writeStandardOutput writes relevant plink files', {
    write_out <- writeStandardOutput(directory=tmp, 
                                     genotypes=genotypes$genotypes, 
                                     phenotypes=phenotypes, 
                                     id_samples = genotypes$id_samples, 
                                     id_snps = genotypes$id_snps, 
                                     id_phenos = colnames(phenotypes), 
                                     kinship=kinship, 
                                     covariates=covariates,
                                     format="plink", verbose=FALSE)
    expect_true(file.exists(paste(tmp, "/Ysim_plink.txt", sep="")))
    expect_true(file.exists(paste(tmp, "/Covs_plink.txt", sep="")))
    expect_false(file.exists(paste(tmp, "/Kinship_plink.txt", sep="")))
    expect_true(file.exists(paste(tmp, "/Genotypes.fam", sep="")))
    expect_true(file.exists(paste(tmp, "/Genotypes.bim", sep="")))
    expect_true(file.exists(paste(tmp, "/Genotypes.bed", sep="")))
})

test_that('writeStandardOutput writes relevant oxgen files', {
    write_out <- writeStandardOutput(directory=tmp, 
                                     genotypes=genotypes$genotypes, 
                                     phenotypes=phenotypes, 
                                     id_samples = genotypes$id_samples, 
                                     id_snps = genotypes$id_snps, 
                                     id_phenos = colnames(phenotypes), 
                                     kinship=kinship, 
                                     covariates=covariates,
                                     format="snptest", verbose=FALSE)
    expect_true(file.exists(paste(tmp, "/Ysim_snptest.sample", sep="")))
    expect_true(file.exists(paste(tmp, "/Genotypes_snptest.gen", sep="")))
    expect_false(file.exists(paste(tmp, "/Kinship_snptest.txt", sep="")))
})


test_that('writeStandardOutput writes relevant gemma files', {
    write_out <- writeStandardOutput(directory=tmp, 
                                     genotypes=genotypes$genotypes, 
                                     phenotypes=phenotypes, 
                                     id_samples = genotypes$id_samples, 
                                     id_snps = genotypes$id_snps, 
                                     id_phenos = colnames(phenotypes), 
                                     kinship=kinship, 
                                     covariates=covariates,
                                     format="gemma", verbose=FALSE)
    expect_true(file.exists(paste(tmp, "/Ysim_gemma.txt", sep="")))
    expect_true(file.exists(paste(tmp, "/Kinship_gemma.txt", sep="")))
    expect_true(file.exists(paste(tmp, "/Covs_gemma.txt", sep="")))
    expect_true(file.exists(paste(tmp, "/Genotypes_gemma.txt", sep="")))
})


test_that('writeStandardOutput writes relevant bimbam files', {
    write_out <- writeStandardOutput(directory=tmp, 
                                     genotypes=genotypes$genotypes, 
                                     phenotypes=phenotypes, 
                                     id_samples = genotypes$id_samples, 
                                     id_snps = genotypes$id_snps, 
                                     id_phenos = colnames(phenotypes), 
                                     kinship=kinship, 
                                     covariates=covariates,
                                     format="bimbam", verbose=FALSE)
    expect_true(file.exists(paste(tmp, "/Ysim_bimbam.txt", sep="")))
    expect_false(file.exists(paste(tmp, "/Covs_bimbam.txt", sep="")))
    expect_false(file.exists(paste(tmp, "/Kinship_bimbam.txt", sep="")))
    expect_true(file.exists(paste(tmp, "/Genotypes_bimbam.csv", sep="")))
})

context('Test savePheno')

test_that('writeStandardOutput fails with path error', {
    expect_error(savePheno(simulation, directory="~/fail",  
                       outstring="Data_simulation", format=c("csv", "plink")),
                 "directory contains ~")
})

test_that('writeStandardOutput writes relevant main csv output', {
    testdir <- savePheno(simulation, directory=tmp,
                         outstring="Data_simulation", format="csv",
                         saveIntermediate=FALSE,
                         verbose=FALSE)
    expect_true(file.exists(paste(testdir, "/Ysim_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Y_genBg_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Y_genFixed_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Y_noiseBg_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Y_noiseFixed_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Genotypes_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/SNP_NrSNP20_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/SNP_effects_NrSNP20_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Covs_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Covs_effect_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Kinship_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/varComponents_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/cov_noiseBg_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/cov_genBg_Data_simulation.csv", 
                                  sep="")))
    unlink(paste(testdir, "/*", sep=""))
})

test_that('writeStandardOutput writes additional csv output', {
    testdir <- savePheno(simulation, directory=tmp,
                         outstring="Data_simulation", format="csv",
                         saveIntermediate=TRUE,
                         verbose=FALSE)
    expect_true(file.exists(paste(testdir, 
                                  "/Y_genBg_independent_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/Y_genBg_shared_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/Y_genFixed_independent_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/Y_genFixed_shared_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/Y_noiseBg_independent_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/Y_noiseBg_shared_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/Y_noiseFixed_independent_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/Y_noiseFixed_shared_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/cov_genBg_independent_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/cov_genBg_shared_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/cov_noiseBg_independent_Data_simulation.csv", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/cov_noiseBg_shared_Data_simulation.csv", 
                                  sep="")))
    unlink(paste(testdir, "/*", sep=""))
})

test_that('writeStandardOutput writes relevant main rds output', {
    testdir <- savePheno(simulation, directory=tmp,
                         outstring="Data_simulation", format="rds",
                         saveIntermediate=FALSE,
                         verbose=FALSE)
    expect_true(file.exists(paste(testdir, "/Ysim_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Y_genBg_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Y_genFixed_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Y_noiseBg_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Y_noiseFixed_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/SNP_NrSNP20_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/SNP_effects_NrSNP20_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Covs_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Covs_effect_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Kinship_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/varComponents_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/cov_noiseBg_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/cov_genBg_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, "/Genotypes_Data_simulation.rds", 
                                  sep="")))
    unlink(paste(testdir, "/*", sep=""))
})

test_that('writeStandardOutput writes additional rds output', {
    testdir <- savePheno(simulation, directory=tmp,
                         outstring="Data_simulation", format="rds",
                         saveIntermediate=TRUE,
                         verbose=FALSE)
    expect_true(file.exists(paste(testdir, 
                                  "/Y_genBg_independent_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/Y_genBg_shared_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/Y_genFixed_independent_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/Y_genFixed_shared_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/Y_noiseBg_independent_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/Y_noiseBg_shared_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/Y_noiseFixed_independent_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/Y_noiseFixed_shared_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/cov_genBg_independent_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/cov_genBg_shared_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/cov_noiseBg_independent_Data_simulation.rds", 
                                  sep="")))
    expect_true(file.exists(paste(testdir, 
                                  "/cov_noiseBg_shared_Data_simulation.rds", 
                                  sep="")))
    unlink(paste(testdir, "/*", sep=""))
})

unlink(paste(tmp, "/*", sep=""))
