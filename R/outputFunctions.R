#' Write simulated data into formats used by standard GWAS software. 
#'
#' writeStandardOutput can write genotypes and phenotypes as well as possible
#' covariates and kinship matrices into a number of formats for standard GWAS 
#' software: plink, snptest, bimbam, gemma, limmbo. For more information on the 
#' different file formats see \emph{External formats}.
#'
#' @param genotypes [NrSamples x NrSNP] Data.frame/matrix of genotypes 
#' [integers]/[doubles].
#' @param phenotypes [NrSamples x NrTrait] Data.frame/matrix of phenotypes 
#' [doubles].
#' @param additionalPhenotypes [NrSamples x NrTrait] Data.frame/matrix of 
#' additional phenotypes (for instance non-linearly tranformed orginal 
#  phenotypes via runSimulation) [doubles].
#' @param covariates [NrSamples x NrCovariates] Data.frame/matrix of covariates
#' [integers]/[doubles].
#' @param kinship [NrSamples x NrSamples] Data.frame/matrix of kinship estimates
#' [doubles].
#' @param eval_kinship [NrSamples] vector with eigenvalues of kinship matrix
#' [doubles].
#' @param evec_kinship [NrSamples x NrSamples] Data.frame/matrix with
#' eigenvectors of kinship matrix [doubles].
#' @param format Vector of name(s) [string] of file formats, options are: 
#' "plink", "snptest", "gemma", "bimbam", "delim". For details on the file 
#' formats see \emph{External formats}.
#' @param id_samples Vector of [NrSamples] sample IDs [string] of simulated 
#' phenotypes, genotypes and covariates.
#' @param id_snps Vector of [NrSNPs] SNP IDs [string] of (simulated) 
#' genotypes.
#' @param id_phenos Vector of [NrTraits] phenotype IDs [string] of 
#' simulated phenotypes.
#' @param standardInput_samples (optional) Data.frame of sample information 
#' obtained when genotypes were read from plink, oxgen or genome file.
#' @param standardInput_genotypes (optional) Data.frame of genotypes obtained 
#' when reading genotypes from plink, oxgen, or genome file.
#' @param directory Absolute path (no tilde expansion) to parent directory 
#' [string] where the data should be saved [needs user writing permission]
#' @param outstring (optional) Name [string] of subdirectory (in relation to 
#' directory) to save set-up independent simulation results.
#' @param intercept_gemma [boolean] When modeling an intercept term in gemma, a 
#' column of 1's have to be appended to the covariate files. Set intercept_gemma
#' to TRUE to include a column of 1's in the output.
#' @param nameAdditional name [string] of additonal phenotypes to be appended to 
#' filename. 
#' @param verbose [boolean]; If TRUE, progress info is printed to standard out
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
#' (genotypes_snptest.gen) and the sample file ending in .sample 
#' (Ysim_snptest.sample). From
#'  \url{https://www.well.ox.ac.uk/~gav/snptest/#input_file_formats}:
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
#' without sample or phenotype header/index (Ysim_bimbam.txt) and b) the mean 
#' genotype file format which is a single file, without information on individuals: 
#' (genotypes.bimbam). From \url{http://www.haplotype.org/bimbam.html}: 
#' The first column of the mean genotype files is the SNP ID, the second and 
#' third columns are allele types with minor allele first. The rest columns are 
#' the mean genotypes of different individuals – numbers between 0 and 2 that 
#' represents the (posterior) mean genotype, or dosage of the minor allele.
#' \item gemma format: consists of a) a simple, tab-separated phenotype file 
#' without sample or phenotype header/index (Ysim_gemma.txt) and b) the mean 
#' genotype file format which is a single file, without information on 
#' individuals(genotypes.gemma); a) and b) both the same as above for bimbam 
#' format). In addition and if applicable, c) a kinship file (kinship_gemma.txt)
#' and d) covariate file (Covs_gemma.txt). From 
#' \url{http://www.xzlab.org/software/GEMMAmanual.pdf}: The kinship file 
#' contains a NrSample × NrSample matrix, where each row and each column 
#' corresponds to individuals in the same order as in the mean genotype file, 
#' and ith row and jth column is a number indicating the relatedness value 
#' between ith and jth individuals. The covariates file has the same format as 
#' the phenotype file dsecribed above and must contain a column of 1’s if one 
#' wants to include an intercept term (set parameter intercept_gemma=TRUE).
#' \item limmbo format: consists of a) a comma-separated phenotype file 
#' without sample IDs as index and phenotype IDs as header (Ysim_limmbo.csv),
#' b) the mean genotype file format with one genetic variant per line. The
#' first column contains the variant ID, column 2-N+1 contain the genotype code
#' (numbers between 0 and 2 that represent the (posterior) mean genotype/dosage
#' of the minor allele) for N samples, c)  a kinship file (kinship_limmbo.csv)
#' and d) covariate file (covs_limmbo.csv). From
#' 
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
#' writeStandardOutput(directory=tempdir(), 
#' genotypes=genotypes$genotypes, phenotypes=phenotypes, 
#' id_samples = genotypes$id_samples, id_snps = genotypes$id_snps, 
#' id_phenos = colnames(phenotypes), format="plink")
#' 
#' # Save in gemma and snptest format
#' writeStandardOutput(directory=tempdir(), 
#' genotypes=genotypes$genotypes, phenotypes=phenotypes, 
#' id_samples = genotypes$id_samples, id_snps = genotypes$id_snps, 
#' id_phenos = colnames(phenotypes), kinship=kinship, 
#' format=c("snptest", "gemma"))
#' }
writeStandardOutput <- function(directory, phenotypes=NULL,
                                genotypes=NULL, 
                                additionalPhenotypes=NULL,
                                covariates=NULL, kinship=NULL,
                                eval_kinship=NULL, evec_kinship=NULL,
                                id_samples, id_snps, id_phenos, 
                                outstring=NULL, 
                                standardInput_samples=NULL,
                                standardInput_genotypes=NULL,
                                format = NULL,
                                intercept_gemma=FALSE, 
                                nameAdditional="_nonLinear",
                                verbose=TRUE) {
    if (is.null(format)) {
        stop("Output format has to be specified, supported formats are plink ",
             "snptest, gemma, limmbo and bimbam")
    }
    if (is.null(genotypes)) {
        vmessage(c("No genotypes provided. Remaining data still saved in ",
                   paste(format, collapse=","), " format."), verbose=verbose)
    }
    if (grepl("~", directory)) {
        stop("directory contains ~: path expansion not guaranteed on 
             every platform (see path.expand{base}), please provide full file
             path to the directory")
    }
    if (!dir.exists(directory)) dir.create(directory, recursive=TRUE)

    if ("plink" %in% format) {
        if (!is.null(phenotypes)) {
            vmessage(c("Save phenotypes in PLINK format"), verbose=verbose)
            if (is.null(standardInput_samples)) {
                pheno_format <- cbind(id_samples, id_samples, phenotypes)
            } else {
                pheno_format <- cbind(standardInput_samples[,1:2], phenotypes)
            }
            write.table(pheno_format, paste(directory, "/Ysim", outstring,
                                     "_plink.txt", sep=""), 
                        sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
            if (!is.null(additionalPhenotypes)) {
                vmessage(c("Save additional phenotypes in PLINK format"), 
                         verbose=verbose)
                if (is.null(standardInput_samples)) {
                    pheno_format <- cbind(id_samples, id_samples, 
                                          additionalPhenotypes)
                } else {
                    pheno_format <- cbind(standardInput_samples[,1:2], 
                                          additionalPhenotypes)
                }
                write.table(pheno_format, paste(directory, "/Ysim",
                                                nameAdditional, outstring,
                                                "_plink.txt", sep=""), 
                            sep="\t", quote=FALSE, col.names=FALSE,
                            row.names=FALSE)
            }
        }
        if (!is.null(genotypes)) {
            vmessage(c("Save genotypes in PLINK format"), verbose=verbose)
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
                                                               "/Genotypes", 
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
            vmessage(c("Save covariates in PLINK format"), verbose=verbose)
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
        if (!is.null(kinship)) {
            vmessage(c("PLINK format does not support a kinship, supplied ",
                       "kinship ignored"), verbose=verbose)
        }
    }
    if ("bimbam" %in% format) {
        if (!is.null(phenotypes)) {
        vmessage(c("Save phenotypes in BIMBAM format"), verbose=verbose)
            write.table(phenotypes, 
                        paste(directory, "/Ysim", outstring, "_bimbam.txt",
                              sep=""),
                        sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
            if (!is.null(additionalPhenotypes)) {
                vmessage(c("Save additional phenotypes in BIMBAM format"), 
                         verbose=verbose)
                write.table(additionalPhenotypes, 
                            paste(directory, "/Ysim", nameAdditional, outstring, 
                                  "_bimbam.txt",
                                  sep=""),
                            sep="\t", quote=FALSE, col.names=FALSE,
                            row.names=FALSE)
            }
        }
        if (!is.null(genotypes)) {
            vmessage(c("Save genotypes in BIMBAM format"), verbose=verbose)
            if (is.null(standardInput_genotypes)) {
                geno_format <- cbind(id_snps, rep("A", length(id_snps)),
                                     rep("B", length(id_snps)), t(genotypes))
            } else {
                geno_format <- standardInput_genotypes
            }
            write.table(geno_format, 
                        paste(directory, "/Genotypes", outstring, "_bimbam.csv", 
                              sep=""),
                        sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)
        }
        if (!is.null(covariates)) {
            vmessage(c("BIMBAM format does not support covariates, supplied ",
                       "covariates ignored"), verbose=verbose)
        }
        if (!is.null(kinship)) {
            vmessage(c("BIMBAM format does not support a kinship, supplied ",
                       "kinship ignored"), verbose=verbose)
        }
    }
    if ("gemma" %in% format) {
        if (!is.null(phenotypes)) {
            vmessage(c("Save phenotypes in GEMMA format"), verbose=verbose)
            write.table(phenotypes, 
                        paste(directory, "/Ysim", outstring, "_gemma.txt", 
                              sep=""), 
                        sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
            if (!is.null(additionalPhenotypes)) {
                vmessage(c("Save additional phenotypes in GEMMA format"), 
                         verbose=verbose)
                write.table(additionalPhenotypes, 
                            paste(directory, "/Ysim", nameAdditional, outstring, 
                                  "_gemma.txt", 
                                  sep=""), 
                            sep="\t", quote=FALSE, col.names=FALSE,
                            row.names=FALSE)
            }
        }
        if (!is.null(genotypes)) {
            vmessage(c("Save genotypes in GEMMA format"), verbose=verbose)
            if (is.null(standardInput_genotypes)) {
                geno_format <- cbind(id_snps, rep("A", length(id_snps)),
                                     rep("B", length(id_snps)), t(genotypes))
            } else {
                geno_format <- standardInput_genotypes
            }
            write.table(geno_format,
                        paste(directory, "/Genotypes", outstring, "_gemma.txt", 
                              sep=""),
                        sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)

        }
        if (!is.null(kinship)) {
            vmessage(c("Save kinship in GEMMA format"), verbose=verbose)
            write.table(kinship, 
                        paste(directory, "/Kinship", outstring, "_gemma.txt", 
                              sep=""), 
                        sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
            write.table(evec_kinship, paste(directory,
                                            "/Kinship_eigenvec_gemma.txt",
                                            sep=""),
                        sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE)
            write.table(eval_kinship, paste(directory,
                                            "/Kinship_eigenval_gemma.txt",
                                            sep=""),
                        sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
        if (!is.null(covariates)) {
            vmessage(c("Save covariates in GEMMA format"), verbose=verbose)
            if (intercept_gemma) {
                covariates <- cbind(rep(1, length(id_samples)), covariates)
            }
            write.table(covariates, 
                        paste(directory, "/Covs", outstring, "_gemma.txt", 
                              sep=""), 
                        sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
    }
    if ("limmbo" %in% format) {
        if (!is.null(phenotypes)) {
            vmessage(c("Save phenotypes in LiMMBo format"), verbose=verbose)
            write.table(phenotypes, 
                        paste(directory, "/Ysim", outstring, "_limmbo.csv", 
                              sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
            if (!is.null(additionalPhenotypes)) {
                vmessage(c("Save additional phenotypes in LiMMBo format"), 
                         verbose=verbose)
                write.table(additionalPhenotypes, 
                            paste(directory, "/Ysim", nameAdditional, outstring, 
                                  "_limmbo.csv", 
                                  sep=""), 
                            sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
            }
        }
        if (!is.null(genotypes)) {
            vmessage(c("Save genotypes in LiMMBo format"), verbose=verbose)
            geno_format <- cbind(id_snps, t(genotypes))
            write.table(geno_format,
                        paste(directory, "/Genotypes", outstring, "_limmbo.csv", 
                              sep=""),
                        sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

        }
        if (!is.null(kinship)) {
            vmessage(c("Save kinship in LiMMBo format"), verbose=verbose)
            write.table(kinship, 
                        paste(directory, "/Kinship", outstring, "_limmbo.csv", 
                              sep=""), 
                        sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)
            write.table(evec_kinship, paste(directory, 
                                            "/Kinship_eigenvec_limmbo.csv",
                                            sep=""),
                        sep=",", quote=FALSE, col.names=colnames(kinship), 
                        row.names=FALSE)
            write.table(eval_kinship, paste(directory, 
                                            "/Kinship_eigenval_limmbo.csv",
                                            sep=""),
                        sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
        if (!is.null(covariates)) {
            vmessage(c("Save covariates in LiMMBo format"), verbose=verbose)
            write.table(covariates, 
                        paste(directory, "/Covs", outstring, "_limmbo.csv", 
                              sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
        }
    }
    if ("snptest" %in% format) {
        if (!is.null(phenotypes)) {
            line2 <- rep("P", length(id_phenos))
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
                vmessage(c("Save phenotypes in SNPTEST format"), 
                         verbose=verbose)
                write.table(pheno_format, 
                            paste(directory, "/Ysim", outstring, 
                                  "_snptest.sample", sep=""), 
                            sep=" ", quote=FALSE, col.names=TRUE, 
                            row.names=FALSE)
            }
            if (!is.null(additionalPhenotypes)) {
                line2 <- rep("P", length(id_phenos))
                additional_pheno_tmp <- rbind(line2, additionalPhenotypes)
    
                if (is.null(standardInput_samples)) {
                    line2 <- c(0, 0, 0)
                    samples_tmp <- cbind(ID_1=id_samples, 
                                         ID_2=id_samples, 
                                         missing=rep(0, length(id_samples)))
                    samples_tmp <- rbind(line2 , samples_tmp)
                    additional_pheno_format <- cbind(samples_tmp, 
                                                     additional_pheno_tmp)
                } else {
                    additional_pheno_format <-
                        cbind(standardInput_samples[,1:3], additional_pheno_tmp)
                }
                rownames(pheno_format) <- 1:(length(id_samples) + 1)
                if (is.null(covariates)) {
                    vmessage(c("Save additional phenotypes in SNPTEST format"), 
                             verbose=verbose)
                    write.table(additional_pheno_format, 
                                paste(directory, "/Ysim", nameAdditional,
                                      outstring, "_snptest.sample", sep=""), 
                                sep=" ", quote=FALSE, col.names=TRUE, 
                                row.names=FALSE)
                }
            }
        }
        if (!is.null(genotypes)) {
            vmessage(c("Save genotypes in SNPTEST format"), verbose=verbose)
            if (is.null(standardInput_genotypes)) {
                probGen <- t(apply(genotypes, 2, expGen2probGen))
                geno_format <- cbind(id_snps, id_snps, 1:length(id_snps),
                                     rep("A", length(id_snps)),
                                     rep("B", length(id_snps)), probGen)
            } else {
                geno_format <- standardInput_genotypes
            }
            write.table(geno_format, 
                        paste(directory, "/Genotypes", outstring,
                              "_snptest.gen", sep=""),
                        col.names=FALSE, row.names=FALSE, quote=FALSE)
        }
        if (!is.null(covariates) && !is.null(phenotypes)) {
            vmessage(c("Save phenotypes and covariates in SNPTEST format"), 
                         verbose=verbose)
            line2 <- rep("C", length(colnames(covariates)))
            line2[grepl("cat",colnames(covariates))] <- "D"                
            line2[grepl("bin",colnames(covariates))] <- "D" 
            covs <- rbind(line2, covariates, make.row.names=FALSE)
            ids <- pheno_format[, 1:3]
            pheno <- pheno_format[, -c(1:3)]
            covs_format <- cbind(ids, covs, pheno)
            write.table(covs_format, 
                        paste(directory, "/Ysim", outstring, "_snptest.sample", 
                              sep=""),
                        sep=" ", col.names=TRUE, row.names=FALSE, quote=FALSE)
            if (!is.null(additionalPhenotypes)) {
                additional_pheno <- additional_pheno_format[, -c(1:3)]
                covs_format <- cbind(ids, covs, additional_pheno)
                write.table(covs_format, 
                            paste(directory, "/Ysim", nameAdditional, 
                                  outstring, "_snptest.sample", 
                                  sep=""),
                            sep=" ", col.names=TRUE, row.names=FALSE, 
                            quote=FALSE)
            }
        }
        if (!is.null(kinship)) {
            vmessage(c("SNPTEST format does not support a kinship, supplied ",
                       "kinship ignored"), verbose=verbose)
        }
    }
}

#' Save final phenotype and phenotype components.
#'
#' savePheno saves simulated phenotypes and their components, model setup 
#' parameters and variance components to the specified directories. Requires a 
#' simulatedData list which is the output of \link{runSimulation}.
#'
#' @param simulatedData Named list of i) dataframe of proportion of variance 
#' explained for each component (varComponents), 
#' ii) a named list with the final simulated phenotype components 
#' (phenoComponentsFinal), iii) a named list with the intermediate simulated 
#' phenotype components (phenoComponentsIntermediate), iv) a named list of 
#' parameters describing the model setup (setup) and v) a named list of raw 
#' components (rawComponents) used for genetic effect simulation (genotypes 
#' and/or kinship); obtained from \link{runSimulation} 
#' @param directory Absolute path (no tilde expansion) to parent directory 
#' [string] where simulated data should be saved [needs user writing 
#' permission]
#' @param outstring Optional name [string] of subdirectory (in relation to 
#' directory) to save set-up dependent simulation results; if
#' set to NULL, subdirectory named by NrSamples, NrSNPs, genetic Model and 
#' noise Model and genVar is created.
#' @param intercept_gemma [boolean] When modeling an intercept term in gemma, a 
#' column of 1's have to be appended to the covariate files. Set intercept_gemma
#' to TRUE to include a column of 1's in the output; only used when 
#' "gemma" \%in\% format
#' @param format Vector of format name(s) [string] specifying the output format;
#' multiple output formats can be requested. Options are: plink, bimbam, 
#' snptest, gemma, limmbo, csv or rds. For information on format see details. In
#' orde to save intermediate phenotype components, at least one of csv or rds
#' need to be specified. plink/bimbam/snptest will only save final
#' phenotype/genotype, kinship and covariate data.
#' @param saveIntermediate [bool] If TRUE, intermediate phenotype components
#' such as shared and independent effect components are saved.
#' @param verbose [boolean]; If TRUE, progress info is printed to standard out
#' @return Path [string] to final output directory. If outstring is 
#' NULL, this directory will be a subdirectory of the input directory.
#' @export
#' @examples
#' simulatedPhenotype <- runSimulation(N=100, P=5, cNrSNP=10,
#' genVar=0.2, h2s=0.2, phi=1)
#' \dontrun{
#' outputdir <- savePheno(simulatedPhenotype, directory=tempdir(),  
#' outstring="Data_simulation", format=c("csv", "plink"))}
savePheno <- function(simulatedData, directory, format=".csv",
                      outstring="", saveIntermediate=TRUE, 
                      intercept_gemma = TRUE, verbose=TRUE) {
    if (grepl("~", directory)) {
        stop("directory contains ~: path expansion not guaranteed on 
             every platform (see path.expand{base}), please provide full file
             path to the directory")
    }
    modelGenetic <- simulatedData$setup$modelGenetic
    modelNoise <- simulatedData$setup$modelNoise
    nrsamples <- simulatedData$setup$N
    nrpheno <- simulatedData$setup$P
    id_samples <-  simulatedData$setup$id_samples
    id_phenos <-  simulatedData$setup$id_phenos
    id_snps <-  simulatedData$setup$id_snps
    NrSNP <- simulatedData$setup$NrCausalSNPs
    genVar <- simulatedData$varComponents$genVar
    rawComponents <-  simulatedData$rawComponents
    phenoComponents <- simulatedData$phenoComponentsFinal
    phenoIntermediate <- simulatedData$phenoComponentsIntermediate
    
    
    ### set-up directories
    if (is.null(outstring)) {
        outstring=paste("samples", nrsamples, "_NrSNP", NrSNP, "_Cg", genVar, 
                        "_model", modelNoise, modelGenetic, sep="")
    }
    
    directory <- file.path(directory, outstring)
    if (!dir.exists(directory)) dir.create(directory, recursive=TRUE)
    
    vmessage(c("Save simulation results"), verbose=verbose)
           
    vmessage(c("Save phenotype to ", directory, "Y..."), verbose=verbose, 
             sep="")

    if (outstring != "") {
        outstring <- paste("_", outstring, sep="")
    }
    if ("rds" %in% format) {
        saveRDS(phenoComponents$Y, 
                paste(directory, "/Ysim", outstring ,".rds", sep=""))
        if (!is.null(phenoComponents$Y_nonLinear)) {
            saveRDS(phenoComponents$Y_nonLinear, 
                    paste(directory, "/Ysim_nonLinear", outstring ,".rds", 
                          sep=""))
        }
    }
    if ("csv" %in% format) {
        write.table(phenoComponents$Y, 
                    paste(directory, "/Ysim", outstring, ".csv", sep=""), 
                    sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
        if (!is.null(phenoComponents$Y_nonLinear)) {
            write.table(phenoComponents$Y_nonLinear, 
                        paste(directory, "/Ysim_nonLinear", outstring, ".csv", 
                              sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
        }
    }
    if (grepl("Bg", modelGenetic)) {
        vmessage(c("Save infinitesimal genetic component to ", directory,
                   "Y_genBg..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_genBg, 
                    paste(directory, "/Y_genBg", outstring,".rds",sep=""))
            saveRDS(phenoComponents$cov_genBg, 
                    paste(directory, "/cov_genBg", outstring,".rds", 
                          sep=""))
            if (saveIntermediate){
                saveRDS(phenoIntermediate$Y_genBg_shared, 
                        paste(directory, "/Y_genBg_shared", outstring, ".rds",
                              sep=""))
                saveRDS(phenoIntermediate$cov_genBg_shared, 
                        paste(directory, "/cov_genBg_shared", outstring, ".rds", 
                              sep=""))
                saveRDS(phenoIntermediate$Y_genBg_independent, 
                        paste(directory, "/Y_genBg_independent", outstring,
                              ".rds", sep=""))
                saveRDS(phenoIntermediate$cov_genBg_independent, 
                        paste(directory, "/cov_genBg_independent", outstring,
                              ".rds", sep=""))
            }
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_genBg, 
                        paste(directory, "/Y_genBg", outstring,".csv", 
                              sep=""),
                        sep=",",quote=FALSE, col.names=NA, row.names=TRUE)
            write.table(phenoComponents$cov_Y_genBg, 
                        paste(directory, "/cov_genBg", outstring,".csv", 
                              sep=""), 
                        sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)
            if (saveIntermediate){
                write.table(phenoIntermediate$Y_genBg_shared, 
                            paste(directory, "/Y_genBg_shared", outstring,
                                  ".csv", sep=""),
                            sep=",",quote=FALSE, col.names=NA, row.names=TRUE)
                write.table(phenoIntermediate$cov_genBg_shared, 
                            paste(directory, "/cov_genBg_shared", outstring,
                                  ".csv", sep=""), 
                            sep=",", quote=FALSE, col.names=FALSE,
                            row.names=FALSE)
                write.table(phenoIntermediate$Y_genBg_independent, 
                            paste(directory, "/Y_genBg_independent", outstring,
                                  ".csv", sep=""),
                            sep=",",quote=FALSE, col.names=NA, row.names=TRUE)
                write.table(phenoIntermediate$cov_genBg_independent, 
                            paste(directory, "/cov_genBg_independent", 
                                  outstring,".csv", sep=""), 
                            sep=",", quote=FALSE, col.names=FALSE, 
                            row.names=FALSE)
            }
        }
        vmessage(c("Save kinship to", directory), verbose=verbose)
        if ("rds" %in% format) {
            saveRDS(rawComponents$kinship, 
                    paste(directory, "/Kinship", outstring,".rds", sep=""))
            saveRDS(rawComponents$evec_kinship, paste(directory, 
                                            "/Kinship_eigenvec", outstring,
                                            ".rds", sep=""))
            saveRDS(rawComponents$eval_kinship, paste(directory, 
                                            "/Kinship_eigenval", outstring,
                                            ".rds", sep=""))
        }
        if ("csv" %in% format) {
            write.table(rawComponents$kinship, 
                        paste(directory, "/Kinship",outstring,".csv", 
                              sep=""), 
                        sep=",", col.names=TRUE, row.names=FALSE)
            colnames(rawComponents$evec_kinship) <- 
                colnames(rawComponents$kinship)
            write.table(rawComponents$evec_kinship, paste(directory, 
                                            "/Kinship_eigenvec", outstring,
                                            ".csv", sep=""),
                        sep=",", quote=FALSE,
                        col.names=TRUE, 
                        row.names=FALSE)
            write.table(rawComponents$eval_kinship, paste(directory, 
                                            "/Kinship_eigenval", outstring,
                                            ".csv", sep=""),
                        sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
    }
        
    if (grepl("Fixed", modelGenetic)) {
        vmessage(c("Save genetic variant effects to ", directory, 
                   "Y_genFixed..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_genFixed, 
                    paste(directory, "/Y_genFixed", outstring, ".rds", 
                          sep=""))
            if (saveIntermediate){
                saveRDS(phenoIntermediate$Y_genFixed_shared, 
                        paste(directory, "/Y_genFixed_shared", outstring,
                              ".rds", sep=""))
                saveRDS(phenoIntermediate$Y_genFixed_independent, 
                        paste(directory, "/Y_genFixed_independent", outstring,
                              ".rds", sep=""))
            }
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_genFixed, 
                        paste(directory, "/Y_genFixed", outstring, ".csv", 
                              sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
            if (saveIntermediate){
                write.table(phenoIntermediate$Y_genFixed_shared, 
                            paste(directory, "/Y_genFixed_shared", outstring,
                                  ".csv", sep=""),
                            sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
                write.table(phenoIntermediate$Y_genFixed_independent, 
                            paste(directory, "/Y_genFixed_independent", 
                                  outstring, ".csv", sep=""),
                            sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
            }
        }
        
        vmessage(c("Save genetic variants and their effect sizes to ", 
                   directory, "SNP..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoIntermediate$genFixed$cov,  
                    paste(directory,"/SNP_NrSNP", NrSNP, outstring,
                          ".rds",sep=""))
            saveRDS(phenoIntermediate$genFixed$cov_effect,  
                    paste(directory, "/SNP_effects_NrSNP", NrSNP, 
                          outstring, ".rds",sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoIntermediate$genFixed$cov,
                        paste(directory, "/SNP_NrSNP", NrSNP,  
                              outstring, ".csv",sep=""), 
                        sep=",", col.names=NA, row.names=TRUE, quote=FALSE)
            write.table(phenoIntermediate$genFixed$cov_effect,  
                        paste(directory, "/SNP_effects_NrSNP", NrSNP, 
                              outstring, ".csv",sep=""), 
                        sep=",", col.names=NA, row.names=TRUE, quote=FALSE)
        }
        if (outstring == "") {
            outstringGenFixed <- "_genFixed"
        } else {
            outstringGenFixed <- paste("_genFixed", outstring, sep="")
        }
        if ("plink" %in% format) {
            writeStandardOutput(directory=directory,
                                genotypes=phenoIntermediate$genFixed$cov,
                                outstring=outstringGenFixed,
                                format="plink",
                                standardInput_samples = 
                                    rawComponents$genotypes$format_files$plink_fam,
                                id_samples=id_samples,
                                id_snps=
                                    colnames(phenoIntermediate$genFixed$cov),
                                verbose=FALSE)
        }
        if ("gemma" %in% format) {
            writeStandardOutput(directory=directory,
                                genotypes=phenoIntermediate$genFixed$cov,
                                outstring=outstringGenFixed,
                                format="gemma",                                
                                id_samples=id_samples,
                                id_snps=
                                    colnames(phenoIntermediate$genFixed$cov),
                                verbose=FALSE)
        }
        if ("limmbo" %in% format) {
            writeStandardOutput(directory=directory,
                                genotypes=phenoIntermediate$genFixed$cov,
                                outstring=outstringGenFixed,
                                format="limmbo",
                                id_samples=id_samples,
                                id_snps=
                                    colnames(phenoIntermediate$genFixed$cov),
                                verbose=FALSE)
        }
        if ("snptest" %in% format) {
            writeStandardOutput(directory=directory,
                                genotypes=phenoIntermediate$genFixed$cov,
                                outstring=outstringGenFixed,
                                format="snptest",
                                standardInput_samples = 
                                    rawComponents$genotypes$format_files$snptest_samples,
                                id_samples=id_samples,
                                id_snps=
                                    colnames(phenoIntermediate$genFixed$cov),
                                verbose=FALSE)
        }
        if ("bimbam" %in% format) {
            writeStandardOutput(directory=directory,
                                genotypes=phenoIntermediate$genFixed$cov,
                                outstring=outstringGenFixed,
                                format="bimbam",                                
                                id_samples=id_samples,
                                id_snps=
                                    colnames(phenoIntermediate$genFixed$cov),
                                verbose=FALSE)
        }
    }
    if (!is.null(rawComponents$genotypes)) {
        vmessage(c("Save genotypes to", directory), verbose=verbose)
        if ("csv" %in% format) {
            write.table(t(rawComponents$genotypes$genotypes), 
                paste(directory, "/Genotypes", outstring, ".csv", sep=""), 
                sep=",", col.names=NA, row.names=TRUE, quote=FALSE)
        }
        if ("rds" %in% format) {
            saveRDS(t(rawComponents$genotypes$genotypes), 
                    paste(directory, "/Genotypes", outstring, ".rds", 
                          sep=""))
        }
    }
    if (is.null(rawComponents$genotypes) && grepl("Fixed", modelGenetic)) {
        vmessage(c("Only causal SNPs sampled from genoFilePrefix-genoFileSuffix",
                   "file were read into memory, so no new genotypes files are",
                   "written"), verbose=verbose)
    }
    if (grepl("Correlated", modelNoise)) {
        vmessage(c("Save correlated background to ", directory, 
                   "Y_correlatedBg..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_correlatedBg, 
                    paste(directory, "/Y_correlatedBg", outstring,".rds", 
                          sep=""))
            saveRDS(phenoComponents$cov_correlatedBg, 
                    paste(directory, "/cov_correlatedBg", outstring,".rds", 
                          sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_correlatedBg,
                        paste(directory, "/Y_correlatedBg", outstring, 
                              ".csv", sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
            write.table(phenoComponents$cov_correlatedBg,
                        paste(directory, "/cov_correlatedBg", outstring, 
                              ".csv", sep=""), 
                        sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
    }
    if (grepl("Bg", modelNoise)) {
        vmessage(c("Save observational noise to ", directory, "Y_noiseBg..."
                   ), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_noiseBg, 
                    paste(directory, "/Y_noiseBg", outstring,".rds", 
                          sep=""))
            saveRDS(phenoComponents$cov_noiseBg, 
                    paste(directory, "/cov_noiseBg", outstring,".rds", 
                          sep=""))
            if (saveIntermediate){
                saveRDS(phenoIntermediate$Y_noiseBg_shared, 
                        paste(directory, "/Y_noiseBg_shared", outstring, ".rds",
                              sep=""))
                saveRDS(phenoIntermediate$cov_noiseBg_shared, 
                        paste(directory, "/cov_noiseBg_shared", outstring, 
                              ".rds", sep=""))
                saveRDS(phenoIntermediate$Y_noiseBg_independent, 
                        paste(directory, "/Y_noiseBg_independent", outstring,
                              ".rds", sep=""))
                saveRDS(phenoIntermediate$cov_noiseBg_independent, 
                        paste(directory, "/cov_noiseBg_independent", outstring,
                              ".rds", sep=""))
            }
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_noiseBg, 
                        paste(directory, "/Y_noiseBg", outstring,".csv", 
                              sep=""), sep=",",
                        quote=FALSE, col.names=NA, row.names=TRUE)
            write.table(phenoComponents$cov_Y_noiseBg, 
                        paste(directory, "/cov_noiseBg",
                              outstring,".csv", sep=""), sep=",",
                        quote=FALSE, col.names=FALSE, row.names=FALSE)
            if (saveIntermediate){
                write.table(phenoIntermediate$Y_noiseBg_shared, 
                            paste(directory, "/Y_noiseBg_shared", outstring,
                                  ".csv", sep=""),
                            sep=",",quote=FALSE, col.names=NA, row.names=TRUE)
                write.table(phenoIntermediate$cov_noiseBg_shared, 
                            paste(directory, "/cov_noiseBg_shared", outstring,
                                  ".csv", sep=""), 
                            sep=",", quote=FALSE, col.names=FALSE,
                            row.names=FALSE)
                write.table(phenoIntermediate$Y_noiseBg_independent, 
                            paste(directory, "/Y_noiseBg_independent", 
                                  outstring, ".csv", sep=""),
                            sep=",",quote=FALSE, col.names=NA, row.names=TRUE)
                write.table(phenoIntermediate$cov_noiseBg_independent, 
                            paste(directory, "/cov_noiseBg_independent", 
                                  outstring,".csv", sep=""), 
                            sep=",", quote=FALSE, col.names=FALSE, 
                            row.names=FALSE)
            }
        }
    }
    
    if (grepl("Fixed", modelNoise)) {
        vmessage(c("Save confounder effects to ", directory, 
                   "Y_noiseFixed..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_noiseFixed, 
                    paste(directory, "/Y_noiseFixed", outstring, ".rds", 
                          sep=""))
            if (saveIntermediate){
                saveRDS(phenoIntermediate$Y_noiseFixed_shared, 
                        paste(directory, "/Y_noiseFixed_shared", outstring,
                              ".rds", sep=""))
                saveRDS(phenoIntermediate$Y_genFixed_independent, 
                        paste(directory, "/Y_noiseFixed_independent", outstring,
                              ".rds", sep=""))
            }
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_noiseFixed, 
                        paste(directory, "/Y_noiseFixed", outstring,
                              ".csv", sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
            if (saveIntermediate){
                write.table(phenoIntermediate$Y_noiseFixed_shared, 
                            paste(directory, "/Y_noiseFixed_shared", outstring,
                                  ".csv", sep=""),
                            sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
                write.table(phenoIntermediate$Y_noiseFixed_independent, 
                            paste(directory, "/Y_noiseFixed_independent", 
                                  outstring, ".csv", sep=""),
                            sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
            }
        }
        
        vmessage(c("Save confounders and their effect sizes to ", 
                   directory, "Covs..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoIntermediate$noiseFixed$cov,
                    paste(directory, "/Covs", outstring, ".rds", sep=""))
            saveRDS(phenoIntermediate$noiseFixed$cov_effect,
                    paste(directory, "/Covs_effect", outstring, ".rds", 
                          sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoIntermediate$noiseFixed$cov, 
                        paste(directory, "/Covs", outstring, ".csv", 
                              sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
            write.table(phenoIntermediate$noiseFixed$cov_effect, 
                        paste(directory, "/Covs_effect", outstring,".csv",
                              sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
        }
    }
    if ("rds" %in% format) {
        saveRDS(simulatedData$varComponents, 
                paste(directory, "/varComponents",  outstring, ".rds", 
                      sep=""))
    }
    if ("csv" %in% format) {
        write.table(simulatedData$varComponents, 
                    paste(directory, "/varComponents", outstring, ".csv",
                          sep=""), 
                    sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)
    }
    if ("plink" %in% format) {
        plink <- writeStandardOutput(phenotypes=phenoComponents$Y, 
                                     additionalPhenotypes=
                                         phenoComponents$Y_nonLinear,
                                     genotypes=rawComponents$genotypes$genotypes,
                                     covariates=phenoIntermediate$noiseFixed$cov,
                                     kinship=rawComponents$kinship,
                                     format="plink", 
                                     standardInput_samples = 
                                         rawComponents$genotypes$format_files$plink_fam,
                                     directory=directory, 
                                     id_samples=id_samples,
                                     id_snps=id_snps,
                                     id_phenos=id_phenos)
    }
    if ("snptest" %in% format) {
        snptest <- writeStandardOutput(phenotypes=phenoComponents$Y, 
                                       additionalPhenotypes=
                                         phenoComponents$Y_nonLinear,
                                       genotypes=rawComponents$genotypes$genotypes,
                                       covariates=phenoIntermediate$noiseFixed$cov,
                                       format="snptest", 
                                       standardInput_genotypes = 
                                           rawComponents$genotypes$format_files$oxgen_genotypes,
                                       standardInput_samples = 
                                           rawComponents$genotypes$format_files$snptest_samples,
                                       directory=directory, 
                                       id_samples=id_samples,
                                       id_snps=id_snps,
                                       id_phenos=id_phenos)
    }
    if ("bimbam" %in% format) {
        bimbam <- writeStandardOutput(phenotypes=phenoComponents$Y, 
                                      additionalPhenotypes=
                                         phenoComponents$Y_nonLinear,
                                      genotypes=rawComponents$genotypes$genotypes,
                                      covariates=phenoIntermediate$noiseFixed$cov,
                                      format="bimbam", 
                                      standardInput_genotypes = 
                                          rawComponents$genotypes$format_files$bimbam_snp_info,
                                      directory=directory, 
                                      id_samples=id_samples,
                                      id_snps=id_snps,
                                      id_phenos=id_phenos)
    }
    if ("gemma" %in% format) {
        gemma <- writeStandardOutput(phenotypes=phenoComponents$Y, 
                                     additionalPhenotypes=
                                         phenoComponents$Y_nonLinear,
                                     genotypes=rawComponents$genotypes$genotypes,
                                     covariates=phenoIntermediate$noiseFixed$cov,
                                     kinship=rawComponents$kinship,
                                     eval_kinship=rawComponents$eval_kinship,
                                     evec_kinship=rawComponents$evec_kinship,
                                     format="gemma", 
                                     intercept_gemma=intercept_gemma, 
                                     standardInput_genotypes = 
                                         rawComponents$genotypes$format_files$bimbam_snp_info,
                                     directory=directory, 
                                     id_samples=id_samples,
                                     id_snps=id_snps,
                                     id_phenos=id_phenos)
    }
    if ("limmbo" %in% format) {
        limmbo <- writeStandardOutput(phenotypes=phenoComponents$Y, 
                                      additionalPhenotypes=
                                          phenoComponents$Y_nonLinear,
                                      genotypes=rawComponents$genotypes$genotypes,
                                      covariates=phenoIntermediate$noiseFixed$cov,
                                      kinship=rawComponents$kinship,
                                      eval_kinship=rawComponents$eval_kinship,
                                      evec_kinship=rawComponents$evec_kinship,
                                      format="limmbo", 
                                      directory=directory, 
                                      id_samples=id_samples,
                                      id_snps=id_snps,
                                      id_phenos=id_phenos)    
    }
    return(directory=directory)
}
