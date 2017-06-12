#' Compute allele frequencies from genotype data.
#'
#' @param snp vector of length N samples with genotypes encoded as 0,1 and 2
#' @return vector with minor and major allele frequency
#' @export
#' @examples
#' # create snp vector with minor allele frequency 0.3
#' snp <- rbinom(200, 2, 0.3)
#' allelefreq <- getAlleleFrequencies(snp)
getAlleleFrequencies <- function(snp) {
    counts <- as.data.frame(table(snp))[,2]
    frequencies <- counts/length(snp)
    major_a <- sqrt(max(frequencies))
    minor_a <- 1 - major_a
    return(c(minor_a, major_a))
}


#' Simulate bi-allelic genotypes.
#'
#' @param N number of samples for which to simulate bi-allelelic genotypes
#' @param NrSNP number of SNPs to simulate
#' @param frequencies vector of allele frequencies [double] from which to sample
#' @param frequencyString alternative to a frequencies vector, a [string] with 
#' frequencies separated by comma
#' can be supplied; most often used when run as a command line application
#' @param sampleID prefix for naming samples (followed by sample number from 1 
#' to N)
#' @param verbose boolean; if TRUE, progress info is printed to standard out
#' @return list of a [N x NrSNP] matrix of simulated genotypes, [N x NrSNP] 
#' matrix of standardised simulated genotypes and vector of sample IDs
#' @seealso \code{\link{standardiseGenotypes}}
#' @export
#' @examples
#' N10NrSNP10 <- simulateGenotypes(N=10, NrSNP=10)
#' N10NrSNP10 <- simulateGenotypes(N=10, NrSNP=10,
#' frequencyString="0.2,0.3,0.4")
simulateGenotypes <- function(N, NrSNP=5000, frequencies=c(0.1, 0.2, 0.4), 
                              sampleID="ID_", verbose=TRUE, 
                              frequencyString=NULL) {
	if (! is.null(frequencyString)) {
		frequencies=commaList2vector(frequencyString)
	}
    if (any(frequencies < 0) || any(frequencies > 1)) {
        stop ("Allele frequencies must be between 0 and 1")
    }
    samples <-paste(sampleID, seq(1, N, 1), sep="")
	vmessage(c("Simulate", NrSNP, "SNPs..."), verbose=verbose)
    X <- sapply(1:NrSNP, function(x) rbinom(N, 2, sample(frequencies, 1)))
    colnames(X) <- paste(rep(1, ncol(X)), "-", 1:ncol(X), "-SNPID", 1:ncol(X), 
                         sep="")
    rownames(X) <- samples
    X_sd <- standardiseGenotypes(X)
	return(list(X=X, X_sd=X_sd, samples=samples))
}

#' Standardise genotypes 
#'
#' Genotypes are standardised as described in Yang et al:
#' snp_standardised = (snp - 2 * minor_allele_freq)/
#' sqrt(2 * minor_allele_freq * major_allele_freq)
#'
#' @param geno [N x NrSNP] matrix/dataframe of genotypes
#' @return [N x NrSNP] matrix of standardised genotypes
#' @seealso \code{\link{getAlleleFrequencies}}
#' @export
#' @references Yang, J., Lee, S.H., Goddard, M.E., Visscher, P.M. (2011) GCTA: 
#' a tool for genome-wide complex trait analysis, AJHG: 88
#' @examples
#' geno <- cbind(rbinom(100, 2, 0.3), rbinom(100, 2, 0.4))
#' geno_sd <- standardiseGenotypes(geno)
standardiseGenotypes <- function(geno) {
    apply(geno, 2, function(snp) {
        allele_freq <- getAlleleFrequencies(snp)
        var_snp <- sqrt(2*allele_freq[1]*allele_freq[2])
        if (var_snp == 0) {
            snp_sd <- snp
        } else {
            snp_sd <- sapply(snp, function(s) (s - 2*allele_freq[1])/var_snp)
        }
        return (snp_sd)
    })
}

#' Draw random SNPs from genotypes.  
#'
#' Draw random SNPs from either simulated genotypes or external genotype files.  
#'
#' @param NrCausalSNPs number [integer] of SNPs to chose at random
#' @param genotypes [NrSamples x totalNrSNPs] matrix of genotypes
#' @param chr numeric vector of chromosomes to chose NrCausalSNPs from; only 
#' used when external genotype data is provided i.e. is.null(genoFilePrefix) 
#' == FALSE
#' @param chr_string [string] alternative to chr, a string with chromosomes 
#' separated by comma; most often used when run as a command line application
#' @param NrChrCausal Number [integer] of causal chromosomes to  chose 
#' NrCausalSNPs from (as opposed to the actual chromosomes to chose from via chr
#' 'chr_string);  only used when external genotype data is provided i.e. 
#' is.null(genoFilePrefix) == FALSE. 
#' @param genoFilePrefix full path/to/chromosome-wise-genotype-file-ending-
#' before-"chrChromosomeNumber" (no '~' expansion!) [string]
#' @param genoFileSuffix [string] following chromosome number including 
#' .fileformat (e.g. ".csv"); has to be a text format i.e. comma/tab/space
#' separated
#' @param genoFileDelimiter field separator [string] of genotype file
#' @param sampleID prefix [string] for naming samples (followed by sample number
#'  from 1 to NrSamples)
#' @param standardise [boolean]; if TRUE standardised genotypes will be returned
#' @return [N x NrCausalSNPs] matrix of randomly drawn, standardised (depending 
#' on standardise option) SNPs 
#' @param verbose [boolean]; if TRUE, progress info is printed to standard out
#' @details In order to chose SNPs from external genotype files without reading 
#' them into memory, genotypes for each chromosome need to be accesible as 
#' [SNPs x samples] in a separate file, containing "chrChromosomenumber" (e.g 
#' chr22) in the file name (e.g. /path/to/dir/related_nopopstructure_chr22.csv).  
#' All genotype files need to be saved in the same directory. genoFilePrefix 
#' (/path/to/dir/related_nopopstructure_) and genoFileSuffix (.csv) specify the 
#' strings leading and following the "chrChromosomenumber". The first column in 
#' each file needs to be the SNP_ID and files cannot contain a header. 
#' getCausalSNPs generates a vector of chromosomes from which to sample the SNPs
#' . For each of the chromosomes, it counts the number of SNPs in the chromosome
#'  file and creates vectors of random numbers ranging from 1:NrSNPSinFile. Only 
#'  the lines corresponding to these numbers are then read into R. The example 
#'  data provided for chromosome 22 contains genotypes (50 samples) of the first 
#'  500 SNPs on chromosome 22 with a minor allele frequency of greater than 2% 
#'  from the European populations of the the 1000 Genomes project.
#'  
#' @seealso \code{\link{standardiseGenotypes}} 
#' @export
#' @examples
#' geno <- simulateGenotypes(N=10, NrSNP=10)
#' causalSNPsFromSimulatedGenoStandardised <- getCausalSNPs(NrCausalSNPs=10,
#' genotypes=geno, 
#' standardise=TRUE)
#'
#' genotypeFile <- system.file("extdata/genotypes/",
#' "genotypes_chr22.csv",
#' package = "PhenotypeSimulator")
#' genoFilePrefix <- gsub("chr.*", "", genotypeFile) 
#' genoFileSuffix <- ".csv" 
#' causalSNPsFromFile <- getCausalSNPs(NrCausalSNPs=10, chr=22, 
#' genoFilePrefix=genoFilePrefix, 
#' genoFileSuffix=genoFileSuffix)
getCausalSNPs <- function(NrCausalSNPs=20,  genotypes=NULL, chr=NULL, 
                          chr_string=NULL, NrChrCausal=NULL,
                          genoFilePrefix=NULL, genoFileSuffix=NULL, 
                          genoFileDelimiter=",", 
                          sampleID="ID_", standardise=FALSE, verbose=TRUE) {
	if (! is.null(genotypes)) {
        if (standardise) genotypes <- genotypes$X_sd
        if (! standardise) genotypes <- genotypes$X
        N <- ncol(genotypes)
        if ( N < NrCausalSNPs) {
            stop(paste("Number of genotypes is less than number of causal SNPs." 
                 , "Increase number of simulated genotypes in simulateGenotypes"
                 , "or decrease number of causal SNPs"))
        }
		causalSNPs <- genotypes[, sample(1:ncol(genotypes), NrCausalSNPs)]
	} else {
	    if (grepl("~", genoFilePrefix)) {
	        stop(paste("genoFilePrefix contains ~: path expansion not", 
	                   "guaranteed on every platform (see path.expand{base}),",
	                   "please provide full file path to genotype files"))
	    }
	    if (all(c( is.null(chr_string), is.null(chr), is.null(NrChrCausal)))) {
	        stop(paste("No information about chromosomes to sample from", 
	                   "provided; please specify either chr_string, chr or", 
                       "NrChrCausal"))
	    }
	    if (all(c( !is.null(chr_string), !is.null(chr))) ||
	        all(c( !is.null(chr_string), !is.null(NrChrCausal))) ||
	        all(c( !is.null(NrChrCausal), !is.null(chr)))) {
	        stop(paste("Too much information for sampling chromosomes provided,"
	                   , "please specifiy only either chr_string, chr or",
                        "NrChrCausal"))
	    }
		if (! is.null(chr_string)) {
			ChrCausal <- commaList2vector(chr_string)
		} else if (! is.null(chr)) {
			ChrCausal <- chr
		} else {
			ChrCausal <- sample(1:22, NrChrCausal)
		}
		NrChrCausal <- length(ChrCausal)
		vmessage(c("Draw", NrCausalSNPs, "causal SNPs from", NrChrCausal, 
		           "chromosomes..."), verbose=verbose)
		NrCausalSNPsChr <- rep(floor(NrCausalSNPs/NrChrCausal), NrChrCausal)
		if ( NrCausalSNPs %% NrChrCausal != 0) {
			addSNP <- sample(NrChrCausal, NrCausalSNPs %% NrChrCausal)
			NrCausalSNPsChr[addSNP] <- NrCausalSNPsChr[addSNP] + 1
		}
		vmessage(c("Causal chromosomes:", ChrCausal), verbose=verbose)
		vmessage(c("Get causal SNPs from chromsome-wide SNP files (", 
		           genoFilePrefix, "...)", sep=""), verbose=verbose)
		
		causalSNPs <- lapply(seq_along(ChrCausal), function(chrom) {
			chromosomefile <- paste(genoFilePrefix, "chr", ChrCausal[chrom], 
			                        genoFileSuffix, sep="")
			SNPsOnChromosome <- R.utils::countLines(chromosomefile) - 1
			if (SNPsOnChromosome <  NrCausalSNPsChr[chrom]) {
			    stop(paste("Number of causal SNPs to be chosen from chromosome", 
			               chr, "is larger than actual number of SNPs provided",
			               "in chromosome file"))
			}
			randomSNPindex <- sample(1:SNPsOnChromosome, NrCausalSNPsChr[chrom])
            randomSNPindex <- randomSNPindex[order(randomSNPindex, 
                                                   decreasing=FALSE)]
 			causalSNPsChr <- read.table(
 			    text=read_lines(chromosomefile, randomSNPindex, sep="\n"), 
 			    sep=genoFileDelimiter, row.names=1)
		})
		
		causalSNPs <- t(do.call(rbind, causalSNPs))
        N <- nrow(causalSNPs)
        rownames(causalSNPs) <- paste(sampleID, seq(1, N, 1), sep="")
        
        if (standardise) {
            vmessage("Standardise SNPs...", verbose=verbose)
            causalSNPs <- standardiseGenotypes(causalSNPs)
        }
	}
	return(causalSNPs)	
}


#' Get genetic kinship.
#'
#' Estimate kinship from standardised genotypes or read pre-computed kinship 
#' file. Standardised genotypes can be obtained via 
#' \code{\link{standardiseGenotypes}} or directly from the 
#' \code{\link{simulateGenotypes}} return object. 
#' Kinship matrices can optionally be normalised.
#'
#' @param X [NrSamples x totalNrSNPs] matrix of standardised genotypes
#' @param sampleID prefix [string] for naming samples (followed by sample number 
#' from 1 to NrSamples)
#' @param norm [boolean], if TRUE kinship matrix will be normalised by the mean 
#' of its diagonal elements and 1e-4 added to the diagonal for numerical 
#' stability
#' @param standardise [boolean], if TRUE genotypes will be standardised before
#' kinship estimation 
#' @param kinshipfile path/to/kinshipfile [string] to be read; either X or 
#' kinshipfile must be provided
#' @param sep field separator [string] of kinship file 
#' @param header [boolean], if TRUE kinship file has header information 
#' @param verbose [boolean]; if TRUE, progress info is printed to standard out
#' @return [NrSamples x NrSamples] matrix of kinship estimate 
#' @details The kinship is estimated as \eqn{K = XX_T}, with X the standardised
#' genotypes of the samples. When estimating the kinship from the provided 
#' genotypes, the kinship should be normalised by the mean of its diagonal 
#' elements and 1e-4 added to the diagonal for numerical stability via norm=TRUE.
#' If a kinship file is provided, normalising can optionally be chosen. For the 
#' provided kinship file, normalisation has already been done a priori and norm 
#' should be set to FALSE. The provided kinship contains estimates for 50 
#' samples across the entire genome. 
#' @export
#' @examples
#' geno <- simulateGenotypes(N=10, NrSNP=50)
#' K_fromGenotypesNormalised <- getKinship(geno$X_sd, norm=TRUE)
#'
#' kinshipfile <- system.file("extdata/kinship", 
#' "kinship.csv",
#' package = "PhenotypeSimulator")
#' K_fromFile <- getKinship(kinshipfile=kinshipfile, norm=FALSE)
getKinship <- function(X=NULL, kinshipfile=NULL, sampleID="ID_", norm=TRUE, 
                       standardise=TRUE, sep=",", header=TRUE, verbose=TRUE) {
    if (!is.null(X)) {
        if (abs(mean(X)) > 0.2 && (sd(X) > 1.2 || sd(X) < 0.8 ) 
            && ! standardise) {
            warning(paste("It seems like genotypes are not standardised, set", 
                    "standardise=TRUE to estimate kinship from standardised",
                    "genotypes (recommended)"))
        }
        if (standardise) {
            X <- standardiseGenotypes(X)
        }
        N <- nrow(X)
        NrSNP <- ncol(X)
        vmessage(c("Estimating kinship from", NrSNP, "SNPs provided"), 
                 verbose=verbose)
        kinship <- X %*% t(X)
        colnames(kinship) <- paste(sampleID, seq(1, N, 1), sep="")
    } else if (!is.null(kinshipfile)) {
        vmessage(c("Reading kinship file from", kinshipfile), verbose=verbose)
        kinship <- as.matrix(read.table(kinshipfile, sep=sep, header=header, 
                                        stringsAsFactors=FALSE))
        if (diff(dim(kinship)) !=0) {
            if (abs(diff(dim(kinship))) == 1) {
                stop (paste("Kinship matrix needs to be a square matrix,",
                "however it has", nrow(kinship), "rows and", ncol(kinship), 
                "columns, did you specify the kinship header information",
                "correctly?"))
            } else {
                stop (paste("Kinship matrix needs to be a square matrix,",
                "however it has", nrow(kinship), "rows and", ncol(kinship), 
                "columns"))
            }
        }
    } else {
        stop ("Either X or kinshipfile must be provided")
    }
    if (norm) {
        vmessage("Normalising kinship", verbose=verbose)
        kinship <- kinship/mean(diag(kinship))
        diag(kinship) <- diag(kinship) + 1e-4
    }
    return(kinship)
}
