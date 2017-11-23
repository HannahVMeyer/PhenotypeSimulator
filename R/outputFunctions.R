#' Write simulated data into formats used by standard GWAS software  
#'
#' writeStandardOutput can write genotypes and phenotypes as well as possible
#' covariates and kinship matrices into a number of formats for standard GWAS: 
#' plink, snptest, bimbam, gemma. Alternatively, simple text files (
#' with specified delimiter) can be read. For more information on the different 
#' file formats see \emph{Externa formats}.
#'
#' @param filename path/to/genotypefile [string] in plink, oxgen 
#' (impute2/snptest/hapgen2), genome, bimbam or [delimiter]-delimited format (
#' for format information see \emph{External genotype software and formats}).
#' @param format name [string] of genotype file format
#' @param id_samples  [string] 
#' @param id_snps [string] 
#' @param intercept_gemma [boolean] When modeling an intercept term in gemma, a 
#' column of 1's have to be appended to the covariate files. Set intercept_gemma
#' to TRUE to include a column of 1's in the output.
#' @param pheno_snptest The snptest .samples file contains sample IDs, 
#' covariates and phenotypes. Thus, when writing the covariates file, the 
#' formated sample ID and phenotype data have to be provided 
#' (output from type = "pheno", format ="snptest")
#' @param delimiter field separator [string] of genotype file
#' @param verbose [boolean]; if TRUE, progress info is printed to standard out
#' @return List of [NrSamples X NrSNPs] genotypes, their [NrSNPs] ID (id_snps), 
#' their [NrSamples] IDs (id_samples) and format specific additional files (
#' might be used for output writing).
#' @section External formats:
#' \itemize{
#' \item PLINK: consists of three files, .bed, .bim and .fam. 
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
#' \item snptest: consists of two files, the genotype file ending in .gen 
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
#' detailing the information for that individual. For more information on the 
#' sampple file visit the above url
#' \item bimbam: Mean genotype file format of bimbam which is a single file, 
#' without information on individuals. From 
#' \url{http://www.haplotype.org/bimbam.html}: The first column of the mean 
#' genotype files is the SNP ID, the second and third columns are allele types 
#' with minor allele first. The rest columns are the mean genotypes of different 
#' individuals – numbers between 0 and 2 that represents the (posterior) mean 
#' genotype, or dosage of the minor allele.
#' }
#' @export
#' @examples 
writeStandardOutput <- function(data, type, directory, 
                                id_samples, id_snps, id_traits, 
                                outstring=NULL, 
                                standardInput=NULL,
                                format = c("plink", "snptest", "gemma", 
                                           "bimbam", "delim"),
                                intercept_gemma=FALSE, 
                                pheno_snptest=NULL,
                                verbose=TRUE, delimiter = ",", ...) {
    if (is.null(format)) {
        stop("Output format has to be specified, supported formats are plink", 
             "snptest, gemma, bimbam, delim (where the delimiter is specified", "
             via 'delimiter=')")
    }
    if (format == "plink") {
        if (type == "pheno") {
            if (is.null(standardInput)) {
                data_format <- cbind(id_samples, id_samples, data)
            } else {
                data_format <- cbind(standardInput[,1:2], data)
            }
            write.table(data_format, paste(directory, "/Ysim_", outstring,
                                     "_plink.txt", sep=""), 
                        sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
        }
        if (type == "geno") {
            if (is.null(standardInput)) {
                X_id <- data.frame(FID=id_samples, IID=id_samples, 
                                   PAT=rep(0, N), MAT=rep(0, N), SEX=rep(0, N), 
                                   PHENOTYPE=rep(-9, N))
                rownames(X_id) <- id_samples
            } else {
                X_id <- standardInput
            }
            plink.out <- write.plink(file.base=paste(directory, "/genotypes_", 
                                                     outstring, sep=""), 
                                         snps=as(geno, "SnpMatrix"), 
                                         sex=X_id$SEX, 
                                         father=X_id$PAT, 
                                         mother=X_id$MAT, 
                                         pedigree=X_id$FID, 
                                         id=X_id$IID, 
                                         phenotype=X_id$PHENOTYPE)
        }
        if (type == "covs") {
            if (is.null(standardInput)) {
                data_format <- cbind(id_samples, id_samples, data)
            } else {
                data_format <- cbind(standardInput[,1:2], data)
            }
            colnames(data_format) <- c("FID", "IID", colnames(data))
            write.table(data_format, paste(directoryPheno, "/Covs_", outstring,
                                           "_plink.txt", sep=""), 
                        sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
        }
    }
    if (format == "bimbam") {
        if (type == "pheno") {
            write.table(data, paste(directory, "/Ysim_", outstring, 
                                "_bimbam.txt", sep=""), sep="\t", 
                        quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
        if (type == "geno") {
            if (is.null(standardInput)) {
                data_format <- cbind(id_snps, rep("A", length(id_snps)),
                                     rep("B", length(id_snps)), data)
            } else {
                data_format <- standardInput
            }
            write.table(data_format, paste(directory, "/genotypes_", outstring, 
                                           ".bimbam", sep=""),
                        col.names=FALSE, row.names=FALSE, quote=FALSE)
        }
        
    }
    if (format == "gemma") {
        if (type == "pheno") {
            write.table(data, paste(directory, "/Ysim_", outstring, 
                                    "_gemma.txt", sep=""), sep="\t", 
                        quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
        if (type == "geno") {
            if (is.null(standardInput)) {
                data_format <- cbind(id_snps, rep("A", length(id_snps)),
                                     rep("B", length(id_snps)), data)
            } else {
                data_format <- standardInput
            }
            write.table(data_format, paste(directory, "/genotypes_", outstring, 
                                           ".gemma", sep=""),
                        col.names=FALSE, row.names=FALSE, quote=FALSE)
            
        }
        if (type == "kinship") {
            write.table(data, paste(directory, "/Kinship_", outstring,
                                    "_gemma.txt", sep=""), sep="\t", 
                        quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
        if (type == "covs") {
            if (intercept) {
                data_format <- cbind(rep(1, length(id_samples)), data)
            }
            write.table(data_format, paste(directory, "/Covs_", outstring,
                                    "_gemma.txt", sep=""), sep="\t", 
                        quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
    }
    if (format == "snptest") {
        if (type == "pheno") {
            line2pheno <- rep("P", length(id_traits))
            pheno_tmp <- rbind(line2pheno, data)
            
            if (is.null(standardInput)) {
                line2samples <- c(0, 0, 0)
                samples_tmp <- cbind(ID_1=id_samples, 
                             ID_2=id_samples, 
                             missing=rep(0, length(id_samples)))
                samples_tmp <- rbind(line2samples , samples_tmp)
                data_format <- cbind(samples_tmp, pheno_tmp)
            } else {
                data_format <- cbind(standardInput[,1:3], pheno_tmp)
            }
            write.table(data_format, paste(directory, "/Ysim_", outstring, 
                                            "_pheno_only_snptest.samples", 
                                            sep=""), 
                        sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)
        }
        if (type == "geno") {
            if (is.null(standardInput)) {
                probGen <- apply(data, 1, expGen2probGen)
                data_format <- cbind(id_snps, id_snps, 1:length(id_snps),
                                     rep("A", length(id_snps)),
                                     rep("B", length(id_snps)), probGen)
            } else {
                data_format <- standardInput
            }
            write.table(data_format, paste(directory, "/genotypes_", outstring, 
                                           ".snptest", sep=""),
                        col.names=FALSE, row.names=FALSE, quote=FALSE)
        }
        if (type == "covs") {
            line2covs <- rep("C", length(colnames(data)))
            line2covs[grepl("cat",colnames(data))] <- "D"                
            line2covs[grepl("bin",colnames(data))] <- "D" 
            covs <- rbind(line2covs, data)
            ids <- pheno_snptest[, 1:3]
            pheno <- pheno_snptest[, -c(1:3)]
            data_format <- cbind(ids, covs, pheno)
            write.table(data_format, paste(directory, "/Ysim_", outstring, 
                                           "_snptest.samples", sep=""),
                        col.names=TRUE, row.names=FALSE, quote=FALSE)
        }
    }
    return(data_format)
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
#' @param directoryGeno absolute path (no tilde expansion) to parent directory 
#' [string] where genotypes from simulations should be saved [needs user writing 
#' permission]
#' @param directoryPheno absolute path (no tilde expansion) to parent directory 
#' [string] where final phenotype and phenotype components from simulations 
#' should be saved [needs user writing permission]
#' @param outstring optional name [string] of subdirectory (in relation to 
#' directoryPheno/directoryGeno) to save set-up
#' independent simulation results
#' @param sample_subset_vec optional vector of sample subset sizes [integer];
#' if provided, draws subsets of samples out of the total simulated dataset and 
#' saves them separately 
#' @param pheno_subset_vec optional vector of phenotype subset sizes [integer] 
#' if provided, draws subsets of traits out of the total simulated dataset 
#' and saves them separately 
#' @param format vector of format names [string] specifying the output format;
#' multiple output formats can be requested. Options are: plink, bimbam, 
#' snptest, csv or rds. For information on format see details. In order to 
#' save intermediate phenotype components, at least one of csv or rds need to 
#' be specified. plink/bimbam/snptest will only save final phenotype/genotype 
#' and covariate data.
#' @param verbose [boolean]; if TRUE, progress info is printed to standard out
#' @return list of paths [strings] to final output phenotype (directoryPheno) 
#' and genotype (directoryGeno) directories. If outstring or subset settings not
#' NULL, these directories will be subdirectories of the input phenotype and 
#' genotype directories.
#' @export
#' @examples
#' simulatedPhenotype <- runSimulation(N=100, P=5, cNrSNP=10,
#' genVar=0.2, h2s=1, phi=1)
#' #not run
#' #outputdir <- savePheno(simulatedPhenotype, directoryGeno="/path/to/dir/",  
#' #directoryPheno="/path/to/dir/", outstring="Date_simulation", 
#' #saveAsPlink=TRUE)
savePheno <- function(simulatedData, directoryGeno, directoryPheno, 
                      format=".csv",
                      sample_subset_vec=NULL, pheno_subset_vec=NULL, 
                      outstring=NULL,  verbose=TRUE) {
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
    sampleID <-  simulatedData$setup$sampleID
    phenoID <-  simulatedData$setup$phenoID
    NrSNP <-simulatedData$setup$NrCausalSNPs
    genVar <- simulatedData$varComponents$genVar
    rawComponents <-  simulatedData$rawComponents
    phenoComponents <- simulatedData$phenoComponents
    
    
    ### set-up directories
    if (is.null(outstring)) {
        outstring=paste("samples", N, "_NrSNP", NrSNP, "_Cg", genVar, "_model", 
                        modelNoise, modelGenetic, sep="")
    }
    
    directoryGeno <- file.path(directoryGeno, outstring)
    ifelse(!dir.exists(directoryGeno), 
           dir.create(directoryGeno, recursive=TRUE), FALSE)
    
    directoryPheno <- file.path(directoryPheno, outstring)
    ifelse(!dir.exists(directoryPheno), 
           dir.create(directoryPheno, recursive=TRUE), FALSE)
    
    vmessage(c("Save simulation results"), verbose=verbose)
           
    vmessage(c("Save phenotype to ", directoryPheno, "/Y..."), verbose=verbose, 
             sep="")

    if ("rds" %in% format) {
        saveRDS(phenoComponents$Y, 
                paste(directoryPheno, "/Ysim_", outstring ,".rds", sep=""))
    if ("csv" %in% format) {
        write.table(phenoComponents$Y, 
                    paste(directoryPheno, "/Ysim_", outstring, ".csv", sep=""), 
                    sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
    }
    if ("plink" %in% format) {
        plink_pheno <- writeStandardOutput(phenoComponents$Y, 
                                          type="pheno", format="plink", 
                                          directory=directoryGeno, 
                                          id_samples=sampleID,
                                          id_snps=snpID,
                                          id_traits=phenoID)
    }
    if ("snptest" %in% format) {
        snptest_pheno <- writeStandardOutput(phenoComponents$Y, 
                                            type="pheno", format="snptest", 
                                            directory=directoryGeno, 
                                            id_samples=sampleID,
                                            id_snps=snpID,
                                            id_traits=phenoID)
    }
    if ("bimbam" %in% format) {
        bimbam_pheno <- writeStandardOutput(phenoComponents$Y, 
                                           type="pheno", format="bimbam", 
                                           directory=directoryGeno, 
                                           id_samples=sampleID,
                                           id_snps=snpID,
                                           id_traits=phenoID)
    }
    if ("gemma" %in% format) {
        gemma_pheno <- writeStandardOutput(phenoComponents$Y, 
                                          type="pheno", format="gemma", 
                                          directory=directoryGeno, 
                                          id_samples=sampleID,
                                          id_snps=snpID,
                                          id_traits=phenoID)
    }
    
    if (grepl("Bg", modelGenetic)) {
        vmessage(c("Save genetic background component to ", directoryPheno,
                   "/Y_genBg..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_genBg, 
                    paste(directoryPheno, "/Y_genBg_", outstring,".rds",sep=""))
            saveRDS(phenoComponents$cov_Y_genBg, 
                    paste(directoryPheno, "/cov_Y_genBg_", outstring,".rds", 
                          sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_genBg, 
                        paste(directoryPheno, "/Y_genBg_", outstring,".csv", 
                              sep=""),
                        sep=",",quote=FALSE, col.names=NA, row.names=TRUE)
            write.table(phenoComponents$cov_Y_genBg, 
                        paste(directoryPheno, "/cov_Y_genBg_", outstring,".csv", 
                              sep=""), 
                        sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
        vmessage(c("Save kinship to", directoryGeno), verbose=verbose)
        if ("rds" %in% format) {
            saveRDS(rawComponents$kinship, 
                    paste(directoryPheno, "/kinship_",outstring,".csv", sep=""), 
        }
        if ("csv" %in% format) {
            write.table(rawComponents$kinship, 
                        paste(directoryPheno, "/kinship_",outstring,".csv", 
                              sep=""), 
                        sep=",", col.names=TRUE, row.names=FALSE)
        }
        if ("gemma" %in% format) {
            gemma_geno <- writeStandardOutput(rawComponents$kinship, 
                                              type="kinship", format="gemma", 
                                              directory=directoryGeno, 
                                              id_samples=sampleID,
                                              id_snps=snpID,
                                              id_traits=phenoID)
        }
    }
        
    if (grepl("Fixed", modelGenetic)) {
        vmessage(c("Save genetic fixed effects to ", directoryPheno, 
                   "/Y_genFixed..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_genFixed, 
                    paste(directoryPheno, "/Y_genFixed_", outstring, ".rds", 
                          sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_genFixed, 
                        paste(directoryPheno, "/Y_genFixed_", outstring, ".csv", 
                              sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
        }
        
        vmessage(c("Save causal SNPs and their effect sizes to ", directoryGeno, 
                   "/SNP..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$genFixed$cov,  
                    paste(directoryGeno,"/SNP_NrSNP", NrSNP, "_", outstring,
                          ".rds",sep=""))
            saveRDS(phenoComponents$genFixed$cov_effect,  
                    paste(directoryGeno, "/SNP_effects_NrSNP", NrSNP, "_", 
                          outstring, ".rds",sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$genFixed$cov,
                        paste(directoryGeno, "/SNP_NrSNP", NrSNP, "_",  
                              outstring, ".csv",sep=""), 
                        sep=",", col.names=NA, row.names=TRUE, quote=FALSE)
            write.table(phenoComponents$genFixed$cov_effect,  
                        paste(directoryGeno, "/SNP_effects_NrSNP", NrSNP, "_", 
                              outstring, ".csv",sep=""), 
                        sep=",", col.names=NA, row.names=TRUE, quote=FALSE)
        }
    }
    if (!is.null(rawComponents$genotypes)) {
        vmessage(c("Save genotypes to", directoryGeno), verbose=verbose)
        if ("csv" %in% format) {
            write.table(t(rawComponents$genotypes$genotypes), 
                paste(directoryGeno, "/genotypes_", outstring, ".csv", sep=""), 
                sep=",", col.names=NA, row.names=TRUE, quote=FALSE)
        }
        if ("rds" %in% format) {
            saveRDS(t(rawComponents$genotypes$genotypes), 
                    paste(directoryGeno, "/genotypes_", outstring, ".rds", 
                          sep=""))
        }
        if ("plink" %in% format) {
            plink_geno <- writeStandardOutput(rawComponents$genotypes$genotypes, 
                                              type="geno", format="plink", 
                                              standardInput = 
                            rawComponents$genotypes$format_files$plink_fam,
                                              directory=directoryGeno, 
                                              id_samples=sampleID,
                                              id_snps=snpID,
                                              id_traits=phenoID)
        }
        if ("snptest" %in% format) {
            snptest_geno <- writeStandardOutput(rawComponents$genotypes$genotypes, 
                                              type="geno", format="snptest", 
                                              directory=directoryGeno, 
                                              standardInput = 
                             rawComponents$genotypes$format_files$oxgen_samples,
                                              id_samples=sampleID,
                                              id_snps=snpID,
                                              id_traits=phenoID)
        }
        if ("bimbam" %in% format) {
            bimbam_geno <- writeStandardOutput(rawComponents$genotypes$genotypes, 
                                                type="geno", format="bimbam", 
                                                directory=directoryGeno, 
                                                standardInput = 
                           rawComponents$genotypes$format_files$bimbam_snp_info,
                                                id_samples=sampleID,
                                                id_snps=snpID,
                                                id_traits=phenoID)
        }
        if ("gemma" %in% format) {
            gemma_geno <- writeStandardOutput(rawComponents$genotypes$genotypes, 
                                                type="geno", format="gemma", 
                                                directory=directoryGeno, 
                                                standardInput = 
                           rawComponents$genotypes$format_files$bimbam_snp_info,
                                                id_samples=sampleID,
                                                id_snps=snpID,
                                                id_traits=phenoID)
        }
    }
    
    if (grepl("Correlated", modelNoise)) {
        vmessage(c("Save correlated background to ", directoryPheno, 
                   "/Y_correlatedBg..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_correlatedBg, 
                    paste(directoryPheno, "/Y_correlatedBg_", outstring,".rds", 
                          sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_correlatedBg,
                        paste(directoryPheno, "/Y_correlatedBg_", outstring, 
                              ".csv", sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
        }
    }
    if (grepl("Bg", modelNoise)) {
        vmessage(c("Save noise background to ", directoryPheno, "/Y_noiseBg..."
                   ), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_noiseBg, 
                    paste(directoryPheno, "/Y_noiseBg_", outstring,".rds", 
                          sep=""))
            saveRDS(phenoComponents$cov_Y_noiseBg, 
                    paste(directoryPheno, "/cov_Y_noiseBg_", outstring,".rds", 
                          sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_noiseBg, 
                        paste(directoryPheno, "/Y_noiseBg_", outstring,".csv", 
                              sep=""), sep=",",
                        quote=FALSE, col.names=NA, row.names=TRUE)
            write.table(phenoComponents$cov_Y_noiseBg, 
                        paste(directoryPheno, "/cov_Y_noiseBg_",
                              outstring,".csv", sep=""), sep=",",
                        quote=FALSE, col.names=FALSE, row.names=FALSE)
        }
    }
    
    if (grepl("Fixed", modelNoise)) {
        vmessage(c("Save noise fixed effects to ", directoryPheno, 
                   "/Y_noiseFixed..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$Y_noiseFixed, 
                    paste(directoryPheno, "/Y_noiseFixed_", outstring, ".rds", 
                          sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$Y_noiseFixed, 
                        paste(directoryPheno, "/Y_noiseFixed_", outstring,
                              ".csv", sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
        }
        
        vmessage(c("Save covariates and their effect sizes to ", 
                   directoryPheno, "/Covs..."), verbose=verbose, sep="")
        if ("rds" %in% format) {
            saveRDS(phenoComponents$noiseFixed$cov,
                    paste(directoryPheno, "/Covs_", outstring, ".rds", sep=""))
            saveRDS(phenoComponents$noiseFixed$cov_effects,
                    paste(directoryPheno, "/Covs_effect_", outstring, ".rds", 
                          sep=""))
        }
        if ("csv" %in% format) {
            write.table(phenoComponents$noiseFixed$cov, 
                        paste(directoryPheno, "/Covs_", outstring, ".csv", 
                              sep=""), 
                        sep=",", quote=FALSE, col.names=NA, row.names=TRUE)
            write.table(phenoComponents$noiseFixed$cov_effects, 
                        paste(directoryPheno, "/Covs_effect_", outstring,".csv",
                              sep=""), 
                        sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)
        }
        if ("plink" %in% format) {
            plink_cov <- writeStandardOutput(phenoComponents$noiseFixed$cov, 
                                              type="covs", format="plink", 
                                              directory=directoryGeno, 
                                              standardInput = 
                            rawComponents$genotypes$format_files$plink_fam,
                                              id_samples=sampleID,
                                              id_snps=snpID,
                                              id_traits=phenoID)
        }
        if ("snptest" %in% format) {
            snptest_cov <- writeStandardOutput(phenoComponents$noiseFixed$cov, 
                                                type="covs", format="snptest", 
                                                directory=directoryGeno, 
                                                standardInput = 
                                                    rawComponents$genotypes$format_files$snptest_samples,
                                                pheno_snptest=snptest_pheno,
                                                id_samples=sampleID,
                                                id_snps=snpID,
                                                id_traits=phenoID)
        }
        if ("gemma" %in% format) {
            gemma_cov <- writeStandardOutput(phenoComponents$noiseFixed$cov, 
                                              type="covs", format="gemma",
                                              intercept_gemma=intercept_gemma, 
                                              directory=directoryGeno, 
                                              standardInput = 
                                                  rawComponents$genotypes$format_files$bimbam_snp_info,
                                              id_samples=sampleID,
                                              id_snps=snpID,
                                              id_traits=phenoID)
        }
        
    }
    
    if ("rds" %in% format) {
        saveRDS(simulatedData$varComponents, 
                paste(directoryPheno, "/varComponents_",  outstring, ".rds", 
                      sep=""))
    }
    if ("csv" %in% format) {
        write.table(simulatedData$varComponents, 
                    paste(directoryPheno, "/varComponents_", outstring, ".csv",
                          sep=""), 
                    sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)
    }
    return(list(directoryPheno=directoryPheno, directoryGeno=directoryGeno))
    }

        
        
        
        