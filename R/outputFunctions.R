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
writeStandardOutput <- function(data, type, directory, outstring,
                                id_samples, id_snps, id_traits, 
                                standardInput=NULL,
                                format = c("plink", "snptest", "gemma", 
                                           "bimbam", "delim"),
                                intercept_gemma=FALSE, covs_snptest=FALSE,
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
            write.table(data_format, paste(directoryPheno, "/Ysim_", outstring,
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
            write.table(data_format, paste(directoryPheno, "/Ysim_", outstring, 
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
        
        
        
        
        