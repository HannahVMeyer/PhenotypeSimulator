#' Write simulated data into formats used by standard GWAS software  
#'
#' writeStandardOutput can write genotypes and phenotypes as well as possible
#' covariates and kinship matrices into a number of formats for standard GWAS: 
#' plink, snptest, bimbam, gemma. Alternatively, simple text files (
#' with specified delimiter) can be written. For more information on the 
#' different file formats see \emph{External formats}.
#'
#' @param genotypes [NrSample x NrSNP] data.frame/matrix of genotypes [integers]
#' /[doubles]
#' @param phenotypes [NrSample x NrTrait] data.frame/matrix of phenotypes 
#' [doubles]
#' @param covariates [NrSample x NrCovariates] data.frame/matrix of covariates
#' [integers]/[doubles]
#' @param kinship [NrSample x NrSamples] data.frame/matrix of kinship estimates
#' [doubles]
#' @param format name [string] of genotype file format, options are: "plink", 
#' "snptest", "gemma", "bimbam", "delim". For details on the file formats see 
#' \emph{External formats}.
#' @param id_samples vector of [NrSamples] sample IDs [string] of simulated 
#' phenotypes, genotypes and covariates
#' @param id_snps vector of [NrSNPs] snp IDs [string] of (simulated) 
#' genotypes
#' @param id_traits vector of [NrTraits] phenotype IDs [string] of 
#' simulated phenotypes
#' @param standardInput_samples data.frame of sample information obtained when 
#' reading genotypes from plink, oxgen or genome file
#' @param standardInput_genotypes data.frame of genotypes obtained when 
#' reading genotypes from plink, oxgen, or genome file
#' @param directory absolute path (no tilde expansion) to parent directory 
#' [string] where the data should be saved [needs user writing permission]
#' @param outstring optional name [string] of subdirectory (in relation to 
#' directoryPheno/directoryGeno) to save set-up
#' independent simulation results
#' @param intercept_gemma [boolean] When modeling an intercept term in gemma, a 
#' column of 1's have to be appended to the covariate files. Set intercept_gemma
#' to TRUE to include a column of 1's in the output.
#' @param delimiter field separator [string] of genotype file
#' @param verbose [boolean]; if TRUE, progress info is printed to standard out
#' @return List of [NrSamples X NrSNPs] genotypes, their [NrSNPs] ID (id_snps), 
#' their [NrSamples] IDs (id_samples) and format specific additional files (
#' might be used for output writing).
#' @section External formats:
#' \itemize{
#' \item plink format: consists of three files, .bed, .bim and .fam. 
#' From \url{https://www.cog-genomics.org/plink/1.9/formats}:  The .bed 
#' files contain the  primary representation of genotype calls at biallelic 
#' variants in a binary format. The .bim is a text file with no header 
#' line, and one line  per variant with the following six fields: i) Chromosome 
#' code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or 
#' name, ii) Variant identifier, iii) Position in morgans or centimorgans (safe 
#' to use dummy value of '0'), iv) Base-pair coordinate (normally 1-based, but 
#' 0 ok; limited to 231-2), v) Allele 1 (corresponding to clear bits in .bed; 
#' usually minor), vi) Allele 2 (corresponding to set bits in .bed; usually 
#' major). The .fam file is a text file with no header line, and one line per 
#' sample with the following six fields: i) Family ID ('FID'), ii), Within-
#' family ID ('IID'; cannot be '0'), iii) Within-family ID of father ('0' if 
#' father isn't in dataset, iv) within-family ID of mother ('0' if mother isn't 
#' in dataset), v) sex code ('1' = male, '2' = female, '0' = unknown), vi) 
#' Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing 
#' data if case/control)
#' \item snptest format: consists of two files, the genotype file ending in .gen 
#' and the sample file ending in .sample. From
#'  \url{http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html}:
#' The genotype file stores data on a one-line-per-SNP format. The first 5 
#' entries of each line should be the SNP ID, RS ID of the SNP, base-pair 
#' position of the SNP, the allele coded A and the allele coded B. The SNP ID 
#' can be used to denote the chromosome number of each SNP. The next three 
#' numbers on the line should be the probabilities of the three genotypes AA, 
#' AB and BB at the SNP for the first individual in the cohort. The next three 
#' numbers should be the genotype probabilities for the second individual in the 
#' cohort. The next three numbers are for the third individual and so on. The 
#' order of individuals in the genotype file should match the order of the 
#' individuals in the sample file. The sample file has three parts (a) a header 
#' line detailing the names of the columns in the file, (b) a line detailing the 
#' types of variables stored in each column, and (c) a line for each individual 
#' detailing the information for that individual. 
#' a) The header line needs a minimum of three entries. The first three entries 
#' should always be ID_1, ID_2 and missing. They denote that the first three 
#' columns contain the first ID, second ID and missing data proportion of each 
#' individual. Additional entries on this line should be the names of covariates 
#' or phenotypes that are included in the file. In the above example, there are 
#' 4 covariates named cov_1, cov_2, cov_3, cov_4, a continuous phenotype named 
#' pheno1 and a binary phenotype named bin1. All phenotypes should appear after 
#' the covariates in this file. b) The second line (the variable type line) 
#' details the type of variables included in each column. The first three 
#' entries of this line should be set to 0. Subsequent entries in this line for 
#' covariates and phenotypes should be specified by the following rules: D for
#' Discrete covariates (coded using positive integers), C for Continuous 
#' covariates, P for Continuous Phenotype, B for Binary Phenotype (0 = Controls, 
#' 1 = Cases). c) Individual information: one line for each individual 
#' containing the information specified by the entries of the header line.
#' Entries of the sample file are separated by spaces.
#' \item bimbam format: consists of a) a simple, tab-separated phenotype file 
#' without sample or phenotype header/index and b) the mean genotype file format 
#' which is a single file, without information on individuals. From 
#' \url{http://www.haplotype.org/bimbam.html}: The first column of the mean 
#' genotype files is the SNP ID, the second and third columns are allele types 
#' with minor allele first. The rest columns are the mean genotypes of different 
#' individuals – numbers between 0 and 2 that represents the (posterior) mean 
#' genotype, or dosage of the minor allele.
#' \item gemma format: consists of a) a simple, tab-separated phenotype file 
#' without sample or phenotype header/index and b) the mean genotype file format 
#' which is a single file, without information on individuals (both the same as 
#' above for bimbam format). In addition and if applicable, c) a kinship file 
#' and d) covariate file. From 
#' \url{http://www.xzlab.org/software/GEMMAmanual.pdf}: The kinship file 
#' contains a NrSample × NrSample matrix, where each row and each column 
#' corresponds to individuals in the same order as in the mean genotype file, 
#' and ith row and jth column is a number indicating the relatedness value 
#' between ith and jth individuals. The covariates file has the same format as 
#' the phenotype file dsecribed above and must contain a column of 1’s if one 
#' wants to include an intercept term (set parameter intercept_gemma=TRUE).
#' }
#' @export
#' @seealso \link{readStandardGenotypes}
#' @examples 
#' simulation <- runSimulation(N=10, P=2, genVar=0.4, h2s=0.2, phi=1)
#' genotypes <- simulation$rawComponents$genotypes
#' kinship <-  simulation$rawComponents$kinship
#' phenotypes <- simulation$phenoComponents$Y
#' 
#' \dontrun{
#' # Save in plink format (.bed, .bim, .fam, Y_sim_plink.txt)
#' writeStandardOutput(directory="/path/to/output", genotypes=genotypes$genotypes, 
#' phenotypes=phenotypes, id_samples = genotypes$id_samples, 
#' id_snps = genotypes$id_snps, id_traits = colnames(phenotypes), 
#' format="plink")
#' 
#' # Save in gemma and snptest format
#' writeStandardOutput(directory="/path/to/output", genotypes=genotypes$genotypes, 
#' phenotypes=phenotypes, id_samples = genotypes$id_samples, 
#' id_snps = genotypes$id_snps, id_traits = colnames(phenotypes),
#'  kinship=kinship, format=c("snptest", "gemma"))
#' }
writeStandardOutput <- function(directory, 
                                genotypes=NULL, phenotypes=NULL, 
                                covariates=NULL, kinship=NULL,
                                id_samples, id_snps, id_traits, 
                                outstring=NULL, 
                                standardInput_samples=NULL,
                                standardInput_genotypes=NULL,
                                format = NULL,
                                intercept_gemma=FALSE, 
                                verbose=TRUE, delimiter = ",") {
    if (is.null(format)) {
        stop("Output format has to be specified, supported formats are plink", 
             "snptest, gemma, bimbam, delim (where the delimiter is specified", "
             via 'delimiter=')")
    }
    if ("plink" %in% format) {
        if (!is.null(phenotypes)) {
            if (is.null(standardInput_samples)) {
                pheno_format <- cbind(id_samples, id_samples, phenotypes)
            } else {
                pheno_format <- cbind(standardInput_samples[,1:2], phenotypes)
            }
            write.table(pheno_format, paste(directory, "/Ysim", outstring,
                                     "_plink.txt", sep=""), 
                        sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
        }
        if (!is.null(genotypes)) {
            if (is.null(standardInput_genotypes)) {
                N <- nrow(genotypes)
                X_id <- data.frame(FID=id_samples, IID=id_samples, 
                                   PAT=rep(0, N), MAT=rep(0, N), SEX=rep(0, N), 
                                   PHENOTYPE=rep(-9, N))
                rownames(X_id) <- id_samples
            } else {
                X_id <- standardInput_genotypes
                colnames(X_id) <- c("FID", "IID", "PAT", "MAT", "SEX", 
                                    "PHENOTYPE")
            }
            plink.out <- snpStats::write.plink(file.base=paste(directory,
                                                               "/genotypes", 
                                                               outstring, 
                                                               sep=""), 
                                         snps=as(genotypes, "SnpMatrix"), 
                                         sex=X_id$SEX, 
                                         father=X_id$PAT, 
                                         mother=X_id$MAT, 
                                         pedigree=X_id$FID, 
                                         id=X_id$IID, 
                                         phenotype=X_id$PHENOTYPE)
        }
        if (!is.null(covariates)) {
            if (is.null(standardInput_samples)) {
                covs_format <- cbind(id_samples, id_samples, covariates)
            } else {
                covs_format <- cbind(standardInput_samples[,1:2], covariates)
            }
            colnames(covs_format) <- c("FID", "IID", colnames(covariates))
            write.table(covs_format, paste(directory, "/Covs", outstring,
                                           "_plink.txt", sep=""), 
                        sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
        }
    }
    if ("bimbam" %in% format) {
        if (!is.null(phenotypes)) {
            write.table(phenotypes, paste(directory, "/Ysim", outstring, 
                                "_bimbam.txt", sep=""), sep="\t", 
                        quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
        if (!is.null(genotypes)) {
            if (is.null(standardInput_genotypes)) {
                geno_format <- cbind(id_snps, rep("A", length(id_snps)),
                                     rep("B", length(id_snps)), t(genotypes))
            } else {
                geno_format <- standardInput_genotypes
            }
            write.table(geno_format, paste(directory, "/genotypes", outstring, 
                                           ".bimbam", sep=""),
                        sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)
        }
        
    }
    if ("gemma" %in% format) {
        if (!is.null(phenotypes)) {
            write.table(phenotypes, paste(directory, "/Ysim", outstring, 
                                    "_gemma.txt", sep=""), sep="\t", 
                        quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
        if (!is.null(genotypes)) {
            if (is.null(standardInput_genotypes)) {
                geno_format <- cbind(id_snps, rep("A", length(id_snps)),
                                     rep("B", length(id_snps)), t(genotypes))
            } else {
                geno_format <- standardInput_genotypes
            }
            write.table(geno_format, paste(directory, "/genotypes", outstring, 
                                           ".gemma", sep=""),
                        sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)
            
        }
        if (!is.null(kinship)) {
            write.table(kinship, paste(directory, "/Kinship", outstring,
                                    "_gemma.txt", sep=""), sep="\t", 
                        quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
        if (!is.null(covariates)) {
            if (intercept_gemma) {
                covariates <- cbind(rep(1, length(id_samples)), covariates)
            }
            write.table(covariates, paste(directory, "/Covs", outstring,
                                    "_gemma.txt", sep=""), sep="\t", 
                        quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
    }
    if ("snptest" %in% format) {
        if (!is.null(phenotypes)) {
            line2 <- rep("P", length(id_traits))
            pheno_tmp <- rbind(line2, phenotypes)
            
            if (is.null(standardInput_samples)) {
                line2 <- c(0, 0, 0)
                samples_tmp <- cbind(ID_1=id_samples, 
                             ID_2=id_samples, 
                             missing=rep(0, length(id_samples)))
                samples_tmp <- rbind(line2 , samples_tmp)
                pheno_format <- cbind(samples_tmp, pheno_tmp)
            } else {
                pheno_format <- cbind(standardInput_samples[,1:3], pheno_tmp)
            }
            rownames(pheno_format) <- 1:(length(id_samples) + 1)
            if (is.null(covariates)) {
                write.table(pheno_format, paste(directory, "/Ysim", outstring, 
                                            "_snptest.samples", 
                                            sep=""), 
                        sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)
            }
        }
        if (!is.null(genotypes)) {
            if (is.null(standardInput_genotypes)) {
                probGen <- t(apply(genotypes, 2, expGen2probGen))
                geno_format <- cbind(id_snps, id_snps, 1:length(id_snps),
                                     rep("A", length(id_snps)),
                                     rep("B", length(id_snps)), probGen)
            } else {
                geno_format <- standardInput_genotypes
            }
            write.table(geno_format, paste(directory, "/genotypes", outstring, 
                                           ".snptest", sep=""),
                        col.names=FALSE, row.names=FALSE, quote=FALSE)
        }
        if (!is.null(covariates)) {
            line2 <- rep("C", length(colnames(covariates)))
            line2[grepl("cat",colnames(covariates))] <- "D"                
            line2[grepl("bin",colnames(covariates))] <- "D" 
            covs <- rbind(line2, covariates, make.row.names=FALSE)
            ids <- pheno_format[, 1:3]
            pheno <- pheno_format[, -c(1:3)]
            covs_format <- cbind(ids, covs, pheno)
            write.table(covs_format, paste(directory, "/Ysim", outstring, 
                                           "_snptest.samples", sep=""),
                        sep=" ", col.names=TRUE, row.names=FALSE, quote=FALSE)
        }
    }
}
        
#' Save final phenotype and phenotype components.
#'
#' savePheno saves model setup parameters and simulated genotypes to the 
#' specified directories. Requires a simulatedData list which is the output of 
#' \link{runSimulation} .
#'
#' @param simulatedData named list of i) dataframe of proportion of variance 
#' explained for each component (varComponents), 
#' ii) a named list with the simulated phenotype components (phenoComponents), 
#' iii) a named list of parameters describing the model setup (setup) and iv) 
#' a named list of raw components (rawComponents) used for genetic effect 
#' simulation (genotypes and/or kinship); obtained from \link{runSimulation} 
#' @param directory absolute path (no tilde expansion) to parent directory 
#' [string] where simulated data should be saved [needs user writing 
#' permission]
#' @param outstring optional name [string] of subdirectory (in relation to 
#' directoryPheno/directoryGeno) to save set-up dependent simulation results; if
#' set to NULL, subdirectory named by NrSamples, NrSNPs, genetic Model and 
#' noise Model and genVar is created.
#' @param intercept_gemma [boolean] When modeling an intercept term in gemma, a 
#' column of 1's have to be appended to the covariate files. Set intercept_gemma
#' to TRUE to include a column of 1's in the output.
#' @param format vector of format names [string] specifying the output format;
#' multiple output formats can be requested. Options are: plink, bimbam, 
#' snptest, gemma, csv or rds. For information on format see details. In order
#' to save intermediate phenotype components, at least one of csv or rds need to 
#' be specified. plink/bimbam/snptest will only save final phenotype/genotype, 
#' kinship and covariate data.
#' @param verbose [boolean]; if TRUE, progress info is printed to standard out
#' @return list of paths [strings] to final output phenotype (directoryPheno) 
#' and genotype (directoryGeno) directories. If outstring is NULL, these 
#' directories will be subdirectories of the input phenotype and genotype directories.
#' @export
#' @examples
#' simulatedPhenotype <- runSimulation(N=100, P=5, cNrSNP=10,
#' genVar=0.2, h2s=0.2, phi=1)
#' \dontrun{
#' outputdir <- savePheno(simulatedPhenotype, directoryGeno="/path/to/dir/",  
#' directoryPheno="/path/to/dir/", outstring="Data_simulation", 
#' format=c("csv", "plink"))}
savePheno <- function(simulatedData, directoryGeno, directoryPheno, 
                      format=".csv",
                      outstring="",  intercept_gemma = TRUE, verbose=TRUE) {
    if (grepl("~", directoryGeno)) {
        stop("directoryGeno contains ~: path expansion not guaranteed on 
             every platform (see path.expand{base}), please provide full file
             path to the genotype directory")
    }
    if (grepl("~", directoryPheno)) {
        stop("directoryPheno contains ~: path expansion not guaranteed on 
             every platform (see path.expand{base}), please provide full file
             path to the phenotype directory")
    }
    modelGenetic <- simulatedData$setup$modelGenetic
    modelNoise <- simulatedData$setup$modelNoise
    nrsamples <- simulatedData$setup$N
    nrpheno <- simulatedData$setup$P
    id_samples <-  simulatedData$setup$id_samples
    id_traits <-  simulatedData$setup$id_traits
    id_snps <-  simulatedData$setup$id_snps
    NrSNP <-simulatedData$setup$NrCausalSNPs
    genVar <- simulatedData$varComponents$genVar
    rawComponents <-  simulatedData$rawComponents
    phenoComponents <- simulatedData$phenoComponents
    
    
    ### set-up directories
    if (is.null(outstring)) {
        outstring=paste("samples", nrsamples, "_NrSNP", NrSNP, "_Cg", genVar, 
                        "_model", modelNoise, modelGenetic, sep="")
    }
    
    directoryGeno <- file.path(directoryGeno, outstring)
    if (!dir.exists(directoryGeno)) dir.create(directoryGeno, recursive=TRUE)
    
    directoryPheno <- file.path(directoryPheno, outstring)
    if (!dir.exists(directoryPheno)) dir.create(directoryPheno, recursive=TRUE)
    
    vmessage(c("Save simulation results"), verbose=verbose)
           
    vmessage(c("Save phenotype to ", directoryPheno, "/Y..."), verbose=verbose, 
             sep="")

    if (outstring != "") {
        outstring <- paste("_", outstring, sep="")
    }
    if ("rds" %in% format) {
        saveRDS(phenoComponents$Y, 
                paste(directoryPheno, "/Ysim", outstring ,".rds", sep=""))
    }
    if ("csv" %in% format) {
        write.table(phenoComponents$Y, 
                    paste(directoryPheno, "/Ysim", outstring, ".csv", sep=""), 
                    sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
    }
    if ("plink" %in% format) {
        plink <- writeStandardOutput(phenotypes=phenoComponents$Y, 
                                     genotypes=rawComponents$genotypes$genotypes,
                                     covariates=phenoComponents$noiseFixed$cov,
                                     kinship=rawComponents$kinship,
                                     format="plink", 
                                     standardInput_samples = 
                                         rawComponents$genotypes$format_files$plink_fam,
                                     directory=directoryPheno, 
                                     id_samples=id_samples,
                                     id_snps=id_snps,
                                     id_traits=id_traits)
    }
    if ("snptest" %in% format) {
        snptest <- writeStandardOutput(phenotypes=phenoComponents$Y, 
                                       genotypes=rawComponents$genotypes$genotypes,
                                       covariates=phenoComponents$noiseFixed$cov,
                                       kinship=rawComponents$kinship,
                                       format="snptest", 
                                       standardInput_genotypes = 
                                           rawComponents$genotypes$format_files$oxgen_genotypes,
                                       standardInput_samples = 
                                           rawComponents$genotypes$format_files$snptest_samples,
                                       directory=directoryPheno, 
                                       id_samples=id_samples,
                                       id_snps=id_snps,
                                       id_traits=id_traits)
    }
    if ("bimbam" %in% format) {
        bimbam <- writeStandardOutput(phenotypes=phenoComponents$Y, 
                                      genotypes=rawComponents$genotypes$genotypes,
                                      covariates=phenoComponents$noiseFixed$cov,
                                      kinship=rawComponents$kinship,
                                      format="bimbam", 
                                      standardInput_genotypes = 
                                          rawComponents$genotypes$format_files$bimbam_snp_info,
                                      directory=directoryPheno, 
                                      id_samples=id_samples,
                                      id_snps=id_snps,
                                      id_traits=id_traits)
    }
    if ("gemma" %in% format) {
        gemma <- writeStandardOutput(phenotypes=phenoComponents$Y, 
                                     genotypes=rawComponents$genotypes$genotypes,
                                     covariates=phenoComponents$noiseFixed$cov,
                                     kinship=rawComponents$kinship,
                                     format="gemma", 
                                     standardInput_genotypes = 
                                         rawComponents$genotypes$format_files$bimbam_snp_info,
                                     directory=directoryPheno, 
                                     id_samples=id_samples,
                                     id_snps=id_snps,
                                     id_traits=id_traits)
    }
    
    if (grepl("Bg", modelGenetic)) {
        vmessage(c("Save genetic background component to ", directoryPheno,
                   "/Y_genBg..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_genBg, 
                    paste(directoryPheno, "/Y_genBg", outstring,".rds",sep=""))
            saveRDS(phenoComponents$cov_Y_genBg, 
                    paste(directoryPheno, "/cov_Y_genBg", outstring,".rds", 
                          sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_genBg, 
                        paste(directoryPheno, "/Y_genBg", outstring,".csv", 
                              sep=""),
                        sep=",",quote=FALSE, col.names=NA, row.names=TRUE)
            write.table(phenoComponents$cov_Y_genBg, 
                        paste(directoryPheno, "/cov_Y_genBg", outstring,".csv", 
                              sep=""), 
                        sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
        vmessage(c("Save kinship to", directoryGeno), verbose=verbose)
        if ("rds" %in% format) {
            saveRDS(rawComponents$kinship, 
                    paste(directoryPheno, "/kinship",outstring,".csv", sep="")) 
        }
        if ("csv" %in% format) {
            write.table(rawComponents$kinship, 
                        paste(directoryPheno, "/kinship",outstring,".csv", 
                              sep=""), 
                        sep=",", col.names=TRUE, row.names=FALSE)
        }
    }
        
    if (grepl("Fixed", modelGenetic)) {
        vmessage(c("Save genetic fixed effects to ", directoryPheno, 
                   "/Y_genFixed..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_genFixed, 
                    paste(directoryPheno, "/Y_genFixed", outstring, ".rds", 
                          sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_genFixed, 
                        paste(directoryPheno, "/Y_genFixed", outstring, ".csv", 
                              sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
        }
        
        vmessage(c("Save causal SNPs and their effect sizes to ", directoryGeno, 
                   "/SNP..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$genFixed$cov,  
                    paste(directoryGeno,"/SNP_NrSNP", NrSNP, outstring,
                          ".rds",sep=""))
            saveRDS(phenoComponents$genFixed$cov_effect,  
                    paste(directoryGeno, "/SNP_effects_NrSNP", NrSNP, 
                          outstring, ".rds",sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$genFixed$cov,
                        paste(directoryGeno, "/SNP_NrSNP", NrSNP,  
                              outstring, ".csv",sep=""), 
                        sep=",", col.names=NA, row.names=TRUE, quote=FALSE)
            write.table(phenoComponents$genFixed$cov_effect,  
                        paste(directoryGeno, "/SNP_effects_NrSNP", NrSNP, 
                              outstring, ".csv",sep=""), 
                        sep=",", col.names=NA, row.names=TRUE, quote=FALSE)
        }
    }
    if (!is.null(rawComponents$genotypes)) {
        vmessage(c("Save genotypes to", directoryGeno), verbose=verbose)
        if ("csv" %in% format) {
            write.table(t(rawComponents$genotypes$genotypes), 
                paste(directoryGeno, "/genotypes", outstring, ".csv", sep=""), 
                sep=",", col.names=NA, row.names=TRUE, quote=FALSE)
        }
        if ("rds" %in% format) {
            saveRDS(t(rawComponents$genotypes$genotypes), 
                    paste(directoryGeno, "/genotypes", outstring, ".rds", 
                          sep=""))
        }
    }
    
    if (grepl("Correlated", modelNoise)) {
        vmessage(c("Save correlated background to ", directoryPheno, 
                   "/Y_correlatedBg..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_correlatedBg, 
                    paste(directoryPheno, "/Y_correlatedBg", outstring,".rds", 
                          sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_correlatedBg,
                        paste(directoryPheno, "/Y_correlatedBg", outstring, 
                              ".csv", sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
        }
    }
    if (grepl("Bg", modelNoise)) {
        vmessage(c("Save noise background to ", directoryPheno, "/Y_noiseBg..."
                   ), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_noiseBg, 
                    paste(directoryPheno, "/Y_noiseBg", outstring,".rds", 
                          sep=""))
            saveRDS(phenoComponents$cov_Y_noiseBg, 
                    paste(directoryPheno, "/cov_Y_noiseBg", outstring,".rds", 
                          sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_noiseBg, 
                        paste(directoryPheno, "/Y_noiseBg", outstring,".csv", 
                              sep=""), sep=",",
                        quote=FALSE, col.names=NA, row.names=TRUE)
            write.table(phenoComponents$cov_Y_noiseBg, 
                        paste(directoryPheno, "/cov_Y_noiseBg",
                              outstring,".csv", sep=""), sep=",",
                        quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
    }
    
    if (grepl("Fixed", modelNoise)) {
        vmessage(c("Save noise fixed effects to ", directoryPheno, 
                   "/Y_noiseFixed..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_noiseFixed, 
                    paste(directoryPheno, "/Y_noiseFixed", outstring, ".rds", 
                          sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_noiseFixed, 
                        paste(directoryPheno, "/Y_noiseFixed", outstring,
                              ".csv", sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
        }
        
        vmessage(c("Save covariates and their effect sizes to ", 
                   directoryPheno, "/Covs..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$noiseFixed$cov,
                    paste(directoryPheno, "/Covs", outstring, ".rds", sep=""))
            saveRDS(phenoComponents$noiseFixed$cov_effects,
                    paste(directoryPheno, "/Covs_effect", outstring, ".rds", 
                          sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$noiseFixed$cov, 
                        paste(directoryPheno, "/Covs", outstring, ".csv", 
                              sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
            write.table(phenoComponents$noiseFixed$cov_effect, 
                        paste(directoryPheno, "/Covs_effect", outstring,".csv",
                              sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
        }
    }
    if ("rds" %in% format) {
        saveRDS(simulatedData$varComponents, 
                paste(directoryPheno, "/varComponents",  outstring, ".rds", 
                      sep=""))
    }
    if ("csv" %in% format) {
        write.table(simulatedData$varComponents, 
                    paste(directoryPheno, "/varComponents", outstring, ".csv",
                          sep=""), 
                    sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)
    }
    return(list(directoryPheno=directoryPheno, directoryGeno=directoryGeno))
}
