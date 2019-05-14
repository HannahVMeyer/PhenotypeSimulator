#' Compute allele frequencies from genotype data.
#'
#' @param snp [N x 1] Vector of length [N] samples with genotypes of a single
#' bi-allelic genetic variant/SNP encoded as 0,1 and 2.
#' @return Vector with ref (0-encoded) and alt (1-encoded) allele frequencies.
#' @export
#' @examples
#' # create snp vector with minor allele frequency 0.3
#' snp <- rbinom(200, 2, 0.3)
#' allelefreq <- getAlleleFrequencies(snp)
getAlleleFrequencies <- function(snp) {
    if (dim(data.frame(snp))[2] != 1) {
        stop ("getAllelelFrequencies takes a [N x 1] vector of genotypes, but ",
              "provided genotype dimensions are: [", dim(snp)[1], " x ",
              dim(snp)[2], "]")
    }
    pp <- length(which(snp == 0))
    pq2 <- length(which(snp == 1))
    qq <- length(which(snp == 2))

    if ((pp + pq2 + qq != length(snp))) {
        stop ("SNP vector contains alleles not encoded as 0, 1 or 2")
    }
    p <- (2*pp +  pq2)/(2*length(snp))
    return(c(p, 1-p))
}

#' Standardise genotypes.
#'
#' Genotypes are standardised as described in Yang et al:
#' snp_standardised = (snp - 2 * ref_allele_freq)/
#' sqrt(2 * ref_allele_freq * alt_allele_freq).
#'
#' @param geno [N x NrSNP] Matrix/dataframe of genotypes [integer]/[double].
#' @return [N x NrSNP] Matrix of standardised genotypes [double].
#' @seealso \code{\link{getAlleleFrequencies}}
#' @export
#' @references Yang, J., Lee, S.H., Goddard, M.E., Visscher, P.M. (2011) GCTA:
#' a tool for genome-wide complex trait analysis, AJHG: 88
#' @examples
#' geno <- cbind(rbinom(2000, 2, 0.3), rbinom(2000, 2, 0.4),rbinom(2000, 2, 0.5))
#' geno_sd <- standardiseGenotypes(geno)
standardiseGenotypes <- function(geno) {
    allele_freq <-  sapply(data.frame(geno),  getAlleleFrequencies)
    var_geno <- sqrt(2*allele_freq[1,]*allele_freq[2,])
    var_geno[var_geno == 0] <- 1
    geno_mean <- sweep(geno, 2, 2*allele_freq[1,], "-")
    geno_sd <- sweep(geno_mean, 2, var_geno, "/")
    return (geno_sd)
}

#' Simulate bi-allelic genotypes.
#'
#' @param N Number of samples for which to simulate bi-allelic genotypes.
#' @param NrSNP Number of SNPs to simulate.
#' @param frequencies Vector of allele frequencies [double] from which to
#' sample.
#' @param sampleID Prefix [string] for naming samples (will be followed by
#' sample number from 1 to N when constructing id_samples).
#' @param snpID Prefix [string] for naming SNPs (will be followed by SNP number
#' from 1 to NrSNP when constructing id_snps).
#' @param verbose [boolean] If TRUE, progress info is printed to standard out.
#' @return Named list with [N x NrSNP] matrix of simulated genotypes
#' (genotypes), their SNP frequencies (freq), a vector of sample IDs
#' (id_samples) and a vector of SNP IDs (id_snps).
#' @seealso \code{\link{standardiseGenotypes}}
#' @export
#' @examples
#' N10NrSNP10 <- simulateGenotypes(N=10, NrSNP=10)
#' N10NrSNP10 <- simulateGenotypes(N=10, NrSNP=10,
#' frequencies=c(0.2,0.3,0.4))
simulateGenotypes <- function(N, NrSNP=5000, frequencies=c(0.1, 0.2, 0.4),
                              sampleID="ID_", snpID="SNP_", verbose=TRUE) {
    numbers <- list(N=N, NrSNP=NrSNP, frequencies=frequencies)
    integers <- list(N=N, NrSNP=NrSNP)
    testNumerics(numbers=numbers, integers=integers)
    if (any(frequencies < 0) || any(frequencies > 1)) {
        stop ("Allele frequencies must be between 0 and 1")
    }
    if (!(is.character(snpID) && length(snpID) == 1)) {
        stop("snpID has to be of length 1 and of type character")
    }
    if (!(is.character(sampleID) && length(sampleID) == 1)) {
        stop("sampleID has to be of length 1 and of type character")
    }
    samples <-paste(sampleID, 1:N, sep="")
    snps <- paste(snpID, 1:NrSNP, sep="")
    vmessage(c("Simulate", NrSNP, "SNPs..."), verbose=verbose)
    freq <- sample(frequencies, NrSNP, replace=TRUE)
    X <- sapply(1:NrSNP, function(x) rbinom(N, 2, freq[x]))
    colnames(X) <- snps
    rownames(X) <- samples
    return(list(genotypes=X, freq = freq, id_snps=snps, id_samples=samples))
}


#' Read genotypes from file.
#'
#' readStandardGenotypes can read genotypes from a number of input
#' formats for standard GWAS (binary plink, snptest, bimbam, gemma) or
#' simulation software (binary plink, hapgen2, genome). Alternatively, simple
#' text files (with specified delimiter) can be read. For more information on
#' the different file formats see \emph{External genotype software and formats}.
#'
#' @param filename path/to/genotypefile [string] in plink, oxgen
#' (impute2/snptest/hapgen2), genome, bimbam or [delimiter]-delimited format (
#' for format information see \emph{External genotype software and formats}).
#' @param format Name [string] of genotype file format. Options are: "plink",
#' "oxgen", "genome", "bimbam" or "delim".
#' @param N Number [integer] of samples to simulate.
#' @param sampleID Prefix [string] for naming samples (will be followed by
#' sample number from 1 to N when constructing id_samples).
#' @param snpID Prefix [string] for naming SNPs (will be followed by SNP number
#' from 1 to NrSNP when constructing id_snps).
#' @param delimiter Field separator [string] of genotype file when
#' format == "delim".
#' @param verbose [boolean] If TRUE, progress info is printed to standard out.
#' @return Named list of [NrSamples X NrSNPs] genotypes (genotypes), their
#' [NrSNPs] SNP IDs (id_snps), their [NrSamples] sample IDs (id_samples) and
#' format-specific additional files (such as format-specific genotypes encoding
#' or sample information; might be used for output writing).
#' @details The file formats and software formates supported are described
#' below. For large external genotypes, consider the option to only read
#' randomly selected, causal SNPs into memory via \link{getCausalSNPs}.
#' @section External genotype software and formats:
#' \subsection{Formats:}{
#' \itemize{
#' \item PLINK: consists of three files: .bed, .bim and .fam. When specifying
#' the filepath, only the core of the name without the ending should be
#' specified (i.e. for geno.bed, geno.bim and geno.fam, geno should be
#' specified). When reading data from plink files, the absolute path to the
#' plink-format file has to be provided, tilde expansion not provided.
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
#' \item oxgen: consists of two files: the space-separated genotype file ending
#' in .gen and the space-separated sample file ending in .sample. When
#' specifying the filepath, only the core of the name without the ending should
#' be specified (i.e. for geno.gen and geno.sample, geno should be specified).
#' From
#'  \url{https://www.well.ox.ac.uk/~gav/snptest/#input_file_formats}:
#' The genotype file stores data on a one-line-per-SNP format. The first five
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
#' sample file visit the above url or see \link{writeStandardOutput}.
#' \item genome: The entire output of genome can be saved via `genome -options >
#' outputfile`. The /path/to/outputfile should be provided and this function
#' extracts the relevant genotype information from this output file.
#' \url{http://csg.sph.umich.edu/liang/genome}
#' \item bimbam: Mean genotype file format of bimbam which is a single, comma-
#' separated file, without information on individuals. From
#' \url{http://www.haplotype.org/bimbam.html}: the first column of the mean
#' genotype files is the SNP ID, the second and third columns are allele types
#' with minor allele first. The remaining columns are the mean genotypes of
#' different individuals – numbers between 0 and 2 that represents the
#' (posterior) mean genotype, or dosage of the minor allele.
#' \item delim: a [delimter]-delimited file of [NrSNPs x NrSamples] genotypes
#' with the snpIDs in the first column and the sampleIDs in the first row and
#' genotypes encoded as numbers between 0 and 2 representing the (posterior)
#' mean genotype, or dosage of the minor allele. Can be user-genotypes or
#' genotypes simulated with foward-time algorithms such as simupop
#' (\url{http://simupop.sourceforge.net/Main/HomePage}) or MetaSim
#' (\url{project.org/web/packages/rmetasim/vignettes/CreatingLandscapes.html}),
#' that allow for user-specified output formats.
#' }}
#'
#' \subsection{Genotype simulation characteristics:}{
#' \itemize{
#' \item PLINK: simple, bi-allelelic genotypes without LD structure, details can
#' be found at \url{https://www.cog-genomics.org/plink/1.9/input#simulate}.
#' \item Hapgen2: resampling-based genotype simulation, details can be found at
#' \url{http://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html}.
#' \item Genome: coalescent-based genotype simulation, details can be found at
#' \url{http://csg.sph.umich.edu/liang/genome/GENOME-manual.pdf}.
#' }
#' Examples on how to call these genotype simulation tools can be found in the
#' vignette \emph{sample-scripts-external-genotype-simulation}.}
#' @export
#' @examples
#' # Genome format
#' filename_genome  <- system.file("extdata/genotypes/genome/",
#' "genotypes_genome.txt",
#' package = "PhenotypeSimulator")
#' data_genome <- readStandardGenotypes(N=100, filename_genome, format ="genome")
#'
#' filename_hapgen  <- system.file("extdata/genotypes/hapgen/",
#' "genotypes_hapgen.controls.gen",
#' package = "PhenotypeSimulator")
#' filename_hapgen <- gsub("\\.gen", "", filename_hapgen)
#' data_hapgen <- readStandardGenotypes(N=100, filename_hapgen, format='oxgen')
#'
#' filename_plink  <- system.file("extdata/genotypes/plink/",
#' "genotypes_plink.bed",
#' package = "PhenotypeSimulator")
#' filename_plink <- gsub("\\.bed", "", filename_plink)
#' data_plink <- readStandardGenotypes(N=100, filename=filename_plink,
#' format="plink")
readStandardGenotypes <- function(N, filename, format = NULL,
                                  verbose=TRUE, sampleID = "ID_",
                                  snpID = "SNP_", delimiter = ",") {
    if (is.null(format)) {
        stop("Genotypefile format has to be specified, supported",
             "formats are plink, oxgen, genome, bimbam and delim (where the
             delimiter is specified via 'delimiter=')  ")
    }
    if (format == "plink") {
        data <- snpStats::read.plink(bed=filename)
        genotypes <- as(data$genotypes, 'numeric')
        if (N > nrow(genotypes)) {
            stop("Sample number specified (", N, ") exceeds number of genotypes
                 provided (",nrow(genotypes), ")")
        }
        if (N < (nrow(genotypes))) {
            Nsamples <- sort(sample(nrow(genotypes), N))
            genotypes <- genotypes[Nsamples,]
            data$fam <-  data$fam[Nsamples,]
        }
        id_samples <- rownames(genotypes)
        id_snps <- colnames(genotypes)
        format_files <- list(plink_fam = data$fam, plink_map = data$map)
    } else if (format == "bimbam") {
        genotypes_raw <- data.table::fread(paste(filename, ".txt", sep=""),
                                           data.table=FALSE, sep=delimiter)
        genotypes <- t(genotypes_raw[, -c(1:3)])
        if (N > nrow(genotypes)) {
            stop("Sample number specified (", N, ") exceeds number of genotypes
                 provided (",nrow(genotypes), ")")
        }
        if (N < (nrow(genotypes))) {
            Nsamples <- sort(sample(nrow(genotypes), N))
            genotypes <- genotypes[Nsamples,]
        }
        bimbam_snp_info <- genotypes_raw[, c(1:3)]
        bimbam_snp_postion <- data.table::fread(paste(filename, ".position.txt",
                                                      sep=""),
                                                data.table=FALSE, sep=delimiter)
        id_samples <- paste(sampleID, nrow(genotypes), sep="")
        id_snps <- genotypes_raw[,1]
        format_files <- list(bimbam_snp_postion=bimbam_snp_postion,
                             bimbam_snp_info=bimbam_snp_info)
    } else if (format == "oxgen") {
        genotypes_raw <- data.table::fread(paste(filename, ".gen", sep=""),
                                           data.table=FALSE)
        samples<- data.table::fread(paste(filename, ".sample", sep=""),
                                    data.table=FALSE)
        genotypes <- apply(genotypes_raw[, -c(1:5)], 1, probGen2expGen)
        if (N > nrow(genotypes)) {
            stop("Sample number specified (", N, ") exceeds number of genotypes
                 provided (",nrow(genotypes), ")")
        }
        if (N < (nrow(genotypes))) {
            allsamples <- nrow(genotypes)
            Nsamples <- sort(sample(allsamples, N))
            genotypes <- genotypes[Nsamples,]
            samples <- samples[c(1,(Nsamples+1)),]

            Nsamples_oxgen <- 5 + sort(c(Nsamples*3, (Nsamples*3)-1,
                                         (Nsamples*3)-2))
            genotypes_raw <- genotypes_raw[ ,c(1:5, Nsamples_oxgen)]
        }
        id_samples <- samples$ID_1[-1]
        id_snps <- genotypes_raw$V2
        format_files <- list(oxgen_samples = samples,
                             oxgen_genotypes=genotypes_raw)
    } else if (format == "genome") {
        data <- data.table::fread( filename, skip="POP1:", sep=" ",
                                   colClasses="character", header=FALSE)
        genotypes <- matrix(as.numeric(unlist(strsplit(data$V2, split=""))),
                            nrow= nrow(data), byrow=TRUE)
        if (N > nrow(genotypes)) {
            stop("Sample number specified (", N, ") exceeds number of genotypes
                 provided (",nrow(genotypes), ")")
        }
        if (N < nrow(genotypes)) {
            Nsamples <- sort(sample(nrow(genotypes), N))
            genotypes <- genotypes[Nsamples,]
        }
        id_samples <- paste(sampleID, 1:N, "_", gsub(":", "", data$V1), sep="")
        id_snps <- paste(snpID, 0:(ncol(genotypes) -1), sep="")
        format_files = NULL
    } else if (format == "delim") {
        data <- data.table::fread(filename, data.table=FALSE,
                                  sep=delimiter)
        id_snps <- data$V1
        genotypes <- t(data[,-1])
        colnames(genotypes) <- id_snps
        if (N > nrow(genotypes)) {
            stop("Sample number specified (", N, ") exceeds number of genotypes
                 provided (",nrow(genotypes), ")")
        }
        if (N < nrow(genotypes)) {
            Nsamples <- sort(sample(nrow(genotypes), N))
            genotypes <- genotypes[Nsamples,]
        }
        id_samples <- rownames(genotypes)
        format_files <- NULL
    } else {
        stop(format, " is not a supported genotype format. Supported",
             "formats are plink, oxgen, genome, bimbam and delim (where the
                 delimiter is specified via 'delimiter=') ")
    }
    colnames(genotypes) <- id_snps
    rownames(genotypes) <- id_samples
    return(list(genotypes=genotypes, id_snps=id_snps, id_samples=id_samples,
                format_files = format_files))
}

#' Draw random SNPs from genotypes.
#'
#' Draw random SNPs from genotypes provided or external genotype files.
#' When drawing from external genotype files, only lines of randomly
#' chosen SNPs are read, which is recommended for large genotype files. See
#' details for more information. The latter option currently supports file in
#' simple delim-formats (with specified delimiter and optional number of fields
#' to skip) and the bimbam and the oxgen format.
#'
#' @param N Number [integer] of samples to simulate.
#' @param NrCausalSNPs Number [integer] of SNPs to chose at random.
#' @param genotypes [NrSamples x totalNrSNPs] Matrix of genotypes [integer]/
#' [double].
#' @param chr Vector of chromosome(s) [integer] to chose NrCausalSNPs from; only
#' used when external genotype data is provided i.e. !is.null(genoFilePrefix).
#' @param NrSNPsOnChromosome Vector of number(s) of SNPs [integer] per entry in
#' chr (see above); has to be the same length as chr. If not provided, number of
#' SNPS in file will be determined from line count (which can be slow for large
#' files); (optional) header lines will be ignored, so accurate number of SNPs
#' not lines in file should be specified.
#' @param NrChrCausal Number [integer] of causal chromosomes to sample
#' NrCausalSNPs from (as opposed to the actual chromosomes to chose from via chr
#' ); only used when external genotype data is provided i.e.
#' !is.null(genoFilePrefix).
#' @param genoFilePrefix full path/to/chromosome-wise-genotype-file-ending-
#' before-"chrChromosomeNumber" (no '~' expansion!) [string].
#' @param genoFileSuffix [string] Following chromosome number including
#' .fileformat (e.g. ".csv"); File described by genoFilePrefix-genoFileSuffix
#' has to be a text format i.e. comma/tab/space separated.
#' @param format Name [string] of genotype file format. Options are:
#' "oxgen", "bimbam" or "delim". See \link{readStandardGenotypes} for details.
#' @param delimiter Field separator [string] of genotypefile or
#' genoFilePrefix-genoFileSuffix file if format == 'delim'.
#' @param skipFields Number [integer] of fields (columns) to skip in
#' genoFilePrefix-genoFileSuffix file if format == 'delim'. See details.
#' @param header [logical] Can be set to indicate if
#' genoFilePrefix-genoFileSuffix file has a header for format == 'delim'. 
#' See details.
#' @param probabilities [boolean]. If set to TRUE, the genotypes in the files
#' described by genoFilePrefix-genoFileSuffix are provided as triplets of
#' probabilities (p(AA), p(Aa), p(aa)) and are converted into their expected
#' genotype frequencies by 0*p(AA) + p(Aa) + 2p(aa) via \link{probGen2expGen}.
#' @param sampleID Prefix [string] for naming samples (will be followed by
#' sample number from 1 to N when constructing id_samples)
#' @param verbose [boolean] If TRUE, progress info is printed to standard out
#' @return [N x NrCausalSNPs] Matrix of randomly drawn genotypes [integer]/
#' [double]
#' @details In order to chose SNPs from external genotype files without reading
#' them into memory, genotypes for each chromosome need to be accesible as
#' [SNPs x samples] in a separate file, containing "chrChromosomenumber" (e.g
#' chr22) in the file name (e.g. /path/to/dir/related_nopopstructure_chr22.csv).
#' All genotype files need to be saved in the same directory. genoFilePrefix
#' (/path/to/dir/related_nopopstructure_) and genoFileSuffix (.csv) specify the
#' strings leading and following the "chrChromosomenumber". If format== delim,
#' the first column in each file needs to be the SNP_ID, the first row can
#' either contain sample IDs or the first row of genotypes (specified with
#' header). Subsequent columns containing additional SNP information can be
#' skipped by setting skipFields. If format==oxgen or bimbam, files need to be
#' in the oxgen or bimbam format (see \link{readStandardGenotypes} for details)
#' and no additional information about delim, header or skipFields will be
#' considered.
#' getCausalSNPs generates a vector of chromosomes from which to sample the
#' SNPs. For each of the chromosomes, it counts the number of SNPs in the
#' chromosome file and creates vectors of random numbers ranging
#' from 1:NrSNPSinFile. Only the lines corresponding to these numbers are then
#' read into R. The example data provided for chromosome 22 contains genotypes
#' (50 samples) of the first 500 SNPs on chromosome 22 with a minor allele
#' frequency of greater than 2% from the European populations of the the 1000
#' Genomes project.
#'
#' @seealso \code{\link{standardiseGenotypes}}
#' @export
#' @examples
#' # get causal SNPs from genotypes simulated within PhenotypeSimulator
#' geno <- simulateGenotypes(N=10, NrSNP=10)
#' causalSNPsFromSimulatedGenoStandardised <- getCausalSNPs(N=10,
#' NrCausalSNPs=10, genotypes=geno$genotypes)
#'
#'# Get causal SNPs by sampling lines from large SNP files
#' genotypeFile <- system.file("extdata/genotypes/",
#' "genotypes_chr22.csv",
#' package = "PhenotypeSimulator")
#' genoFilePrefix <- gsub("chr.*", "", genotypeFile)
#' genoFileSuffix <- ".csv"
#' causalSNPsFromLines <- getCausalSNPs(N=50, NrCausalSNPs=10, chr=22,
#' genoFilePrefix=genoFilePrefix,
#' genoFileSuffix=genoFileSuffix)
getCausalSNPs <- function(N, NrCausalSNPs=20,  genotypes=NULL,
                          chr=NULL, NrSNPsOnChromosome=NULL,
                          NrChrCausal=NULL, genoFilePrefix=NULL,
                          genoFileSuffix=NULL, format='delim',
                          delimiter=",", header=FALSE,
                          skipFields=NULL, probabilities=FALSE,
                          sampleID="ID_", verbose=TRUE) {
    if (! is.null(genotypes)) {
        if (!(is.matrix(genotypes) || is.data.frame(genotypes))) {
            stop("Genotypes have to be provided as a [NrSamples x NrSNP] matrix",
                 " or data.frame but are provided as ", typeof(genotypes))
        }
        if ( ncol(genotypes) < NrCausalSNPs) {
            stop(paste("Number of genotypes is less than number of causal SNPs."
                       , "Increase number of simulated genotypes in",
                       "simulateGenotypes or decrease number of causal SNPs"))
        }
        drawSNPs <- sort(sample(ncol(genotypes), NrCausalSNPs))
        causalSNPs <- as.matrix(genotypes[,drawSNPs])
        colnames(causalSNPs) <-colnames(genotypes)[drawSNPs]
        if (!is.numeric(causalSNPs)) {
            stop(paste("Provided genotypes are not numeric."))
        }
    } else if (!is.null(genoFilePrefix)) {
        if (grepl("~", genoFilePrefix)) {
            stop(paste("genoFilePrefix contains ~: path expansion not",
                       "guaranteed on every platform (see path.expand{base}),",
                       "please provide full file path to genotype files"))
        }
        if (all(c(is.null(chr), is.null(NrChrCausal)))) {
            stop(paste("No information about chromosomes to sample from",
                       "provided; please specify either chr or",
                       "NrChrCausal"))
        }
        if (all(c( !is.null(NrChrCausal), !is.null(chr)))) {
            stop(paste("Too much information for sampling chromosomes provided,"
                       , "please specifiy only either chr or",
                       "NrChrCausal"))
        }
        if (! is.null(chr)) {
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
        if (!is.null(NrSNPsOnChromosome ) &&
            length(ChrCausal) != length(NrSNPsOnChromosome)) {
            stop("Not enough information about numbers of SNPs per chromosome",
                 "provided. Number of causal chromosomes: ", length(ChrCausal),
                 ", Information about SNPs on these chromosomes given for ",
                 length(NrSNPsOnChromosome), ".")
        }
        skipRows <- 0
        if (format == "oxgen") {
            delimiter <- " "
            skipFields <- 5
            probabilities <- TRUE
        } else if (format == "bimbam") {
            delimiter <- ","
            skipFields <- 3
            probabilities <- FALSE
        } else if (format == "delim") {
            if (is.null(skipFields)) skipFields <- 1
            if (header) skipRows <- 1
        } else {
            stop("Format: ", format, "not supported for sampling genotypes
                     from file.")
        }
        
        vmessage(c("Get causal SNPs from chromsome-wide SNP files (",
                   genoFilePrefix, "...)", sep=""), verbose=verbose)

        causalSNPs <- lapply(seq_along(ChrCausal), function(chrom) {
            vmessage(c("Get", NrCausalSNPsChr[chrom],
                       "causal SNPs from chromsome", ChrCausal[chrom], "..."),
                     verbose=verbose)
            chromosomefile <- paste(genoFilePrefix, "chr", ChrCausal[chrom],
                                    genoFileSuffix, sep="")
            if (!file.exists(chromosomefile)) {
                stop(chromosomefile , " does not exist, did you specify the ",
                    "correct genoFilePrefix and genoFileSuffix?") 
            }
            if(is.null(NrSNPsOnChromosome)) {
                vmessage(c("Count number of SNPs on chromosome",
                           ChrCausal[chrom], "...", sep=""), verbose=verbose)
                SNPsOnChromosome <- R.utils::countLines(chromosomefile) -
                    skipRows
            } else {
                SNPsOnChromosome <- NrSNPsOnChromosome[chrom]
            }
            if (SNPsOnChromosome <  NrCausalSNPsChr[chrom]) {
                stop(paste("Number of causal SNPs to be chosen from chromosome", 
                           ChrCausal[chrom], "is larger than actual number of", 
                           "SNPs provided In chromosome file"))
            }
            vmessage(c("Sample SNPs on chromosome", ChrCausal[chrom], "..."),
                     verbose=verbose)
            randomSNPindex <- sample((skipRows + 1):SNPsOnChromosome,
                                     NrCausalSNPsChr[chrom])
            randomSNPindex <- randomSNPindex[order(randomSNPindex,
                                                   decreasing=FALSE)]
            text <- read_lines(chromosomefile, randomSNPindex, sep="\n")

            if (! grepl(delimiter, text[1])) {
                stop("Delimiter specified for genoFilePrefix-genoFileSuffix (",
                     delimiter, ") cannot be found in file. Did you specify",
                    " the correct delimiter and/or file format?")
            }
            causalSNPsChr <- read.table(text=text, sep=delimiter,
                                        stringsAsFactors=FALSE)
            
            if (format == "oxgen") {
                causalSNPsChr[,2] <- gsub("-", ":", causalSNPsChr[,2])
                causalSNPsChr[,3] <- gsub("", "", causalSNPsChr[,3])
                rownames(causalSNPsChr) <- apply(causalSNPsChr[,1:5], 1, paste,
                                                 collapse= "-")

                samplefile <- paste(genoFilePrefix, "chr", ChrCausal[chrom], 
                                    gsub("gen", "sample", genoFileSuffix),
                                    sep="")
                if (!file.exists(samplefile)) {
                    stop(samplefile , " which is required in oxgen format, ",
                         "does not exist. Did you specify the",
                         "correct genoFilePrefix and genoFileSuffix?") 
                }
                samples <- read.table(samplefile, sep=delimiter,
                                            stringsAsFactors=FALSE, skip=2)[,1]
            } else if (format == "bimbam") {
                rownames(causalSNPsChr) <- apply(causalSNPsChr[,1:3], 1, paste,
                                                 collapse= "-")
                samples <- paste(sampleID, 1:nrow(causalSNPsChr), sep="")
            } else {
                rownames(causalSNPsChr) <- causalSNPsChr[,1]
                if (header) {
                    samples <- read.table(text=read_lines(chromosomefile, 1,
                                                          sep="\n"),
                                          sep=delimiter,
                                          stringsAsFactors=FALSE)[,-1]
                } else {
                    samples <- paste(sampleID, 1:nrow(causalSNPsChr), sep="")
                }

            } 
            if (!is.null(skipFields)) {
                causalSNPsChr <- causalSNPsChr[,-c(1:skipFields)]
            }
            if (probabilities) {
                causalSNPsChr <- t(apply(causalSNPsChr, 1, probGen2expGen))
            }
            colnames(causalSNPsChr) <- samples
            return(causalSNPsChr)
        })
        causalSNPs <- t(do.call(rbind, causalSNPs))
        if (!is.numeric(causalSNPs)) {
            stop(paste("Sampled genotypes are not numeric. Did you specify the",
                       "correct format and/or number of skipFields?"))
        }
    } else {
        stop(paste("No genotype information provided, please specify either",
                   "genotypefile to read genotypes from file, or",
                   "genoFilePrefix and genoFileSuffix to sample lines from",
                   "file or directly provide [N x S] matrix of [S] genotypes",
                   "from [N] samples"))
    }
    NrGenotypeSamples <- nrow(causalSNPs)
    if (N > NrGenotypeSamples) {
        stop("Sample number specified exceeds number of genotypes provided")
    }
    if (N < NrGenotypeSamples) {
        vmessage(c("Sampling", N, "samples from", NrGenotypeSamples,
                   "genotypes provided"), verbose = verbose)
        Nsamples <- sort(sample(NrGenotypeSamples, N))
        causalSNPs <- causalSNPs[Nsamples,]
    }
    return(causalSNPs)
}


#' Get genetic kinship.
#'
#' Estimate kinship from standardised genotypes or read pre-computed kinship
#' file. Standardised genotypes can be obtained via
#' \code{\link{standardiseGenotypes}}.
#'
#' @param N Number [integer] of samples to simulate.
#' @param X [NrSamples x totalNrSNPs] Matrix of (standardised) genotypes.
#' @param sampleID Prefix [string] for naming samples (will be followed by
#' sample number from 1 to N when constructing id_samples).
#' @param id_samples Vector of [NrSamples] sample IDs [string]; if not provided
#' constructed by paste(sampleID, 1:N, sep="").
#' @param standardise [boolean] If TRUE genotypes will be standardised before
#' kinship estimation.
#' @param kinshipfile path/to/kinshipfile [string] to be read; either X or
#' kinshipfile must be provided.
#' @param sep Field separator [string] of kinship file.
#' @param header [boolean], If TRUE kinship file has header information.
#' @param verbose [boolean]; If TRUE, progress info is printed to standard out
#' @return [NrSamples x NrSamples] Matrix of kinship estimate.
#' @details The kinship is estimated as \eqn{K = XX_T}, with X the standardised
#' genotypes of the samples. When estimating the kinship from the provided
#' genotypes, the kinship is normalised by the mean of its diagonal
#' elements and 1e-4 added to the diagonal for numerical stability.
#' @export
#' @examples
#' geno <- simulateGenotypes(N=10, NrSNP=50)
#' K_fromGenotypesNormalised <- getKinship(N=10, X=geno$genotypes,
#' standardise=TRUE)
#'
#' kinshipfile <- system.file("extdata/kinship",
#' "kinship.csv",
#' package = "PhenotypeSimulator")
#' K_fromFile <- getKinship(N=50, kinshipfile=kinshipfile)
getKinship <- function(N, sampleID="ID_", X=NULL, kinshipfile=NULL,
                       id_samples=NULL,
                       standardise=FALSE, sep=",", header=TRUE, verbose=TRUE) {
    numbers <- list(N=N)
    positives <- list(N=N)
    integers <- list(N=N)
    testNumerics(numbers=numbers, positives=positives, integers=integers)

    if(is.null(id_samples)) {
        if (!(is.character(sampleID) && length(sampleID) == 1)) {
            stop("sampleID has to be of length 1 and of type character")
        }
        id_samples <- paste(sampleID, 1:N, sep="")
    }
    if (!is.null(X)) {
        N <- nrow(X)
        NrSNP <- ncol(X)
        if (standardise) {
            vmessage(c("Standardising the", NrSNP, "SNPs provided"),
                     verbose=verbose)
            X <- standardiseGenotypes(X)
        }
        vmessage(c("Estimating kinship from", NrSNP, "SNPs provided"),
                 verbose=verbose)
        kinship <- tcrossprod(X)
        vmessage("Normalising kinship", verbose=verbose)
        kinship <- kinship/mean(diag(kinship))
        # Make numerically stable
        diag(kinship) <- diag(kinship) + 1e-4
        colnames(kinship) <- id_samples
    } else if (!is.null(kinshipfile)) {
        if (!file.exists(kinshipfile)) {
            stop(kinshipfile , " does not exist.")
        }
        vmessage(c("Reading kinship file from", kinshipfile), verbose=verbose)
        kinship <- as.matrix(data.table::fread(kinshipfile, sep=sep,
                                               header=header,
                                               data.table=FALSE,
                                               stringsAsFactors=FALSE))
        if (diff(dim(kinship)) !=0) {
            if (abs(diff(dim(kinship))) == 1) {
                stop ("Kinship matrix needs to be a square matrix,",
                      "however it has ", nrow(kinship), " rows and ",
                      ncol(kinship), " columns, did you specify the kinship ",
                      "header information correctly?")
            } else {
                stop ("Kinship matrix needs to be a square matrix ,",
                      "however it has ", nrow(kinship), " rows and ",
                      ncol(kinship), " columns")
            }
        }

        if (!header) {
            vmessage(c("Kinship does not have column names, sample names",
                       "created by paste(sampleID, 1:N, sep='')"), verbose=
                         verbose)
            colnames(kinship) <- paste(sampleID, 1:ncol(kinship), sep='')
        }

        NrKinshipSamples <- ncol(kinship)
        if (N > NrKinshipSamples) {
            stop("Number of samples specifid is greater than number of ",
                 "samples in kinship matrix")
        }
        if (N < NrKinshipSamples) {
            vmessage(c("Number of samples specifid is smaller than number of",
                       "samples in kinship matrix"), verbose=verbose)
            if (length(id_samples) == N) {
                if (all(id_samples %in% colnames(kinship))) {
                    vmessage(c("Extract kinship samples based on provided",
                               "id_samples"), verbose=verbose)
                    kinship <- kinship[which(colnames(kinship) %in% id_samples),
                                       which(colnames(kinship) %in% id_samples)]
                } else {
                    stop("Not all id_samples are present in sample names of ",
                         "kinship (colnames), subsampling not possible. Please",
                         " check dimensions and/or names of kinship")
                }
            }
            else if (length(id_samples) == NrKinshipSamples) {
                vmessage(c("Sampling", N, "samples from", NrKinshipSamples,
                           "samples provided in kinship"), verbose=verbose)
                Nsamples <- sort(sample(NrKinshipSamples, N))
                kinship <- kinship[Nsamples,Nsamples]
                id_samples <- colnames(kinship)
            } else {
                stop("Length of id_samples (", length(id_samples), ") matches
                    neither the number of samples in the kinship (",
                     NrKinshipSamples , ") nor the specified sample number (",
                     N, "). Please check dimensions and/or names of kinship")
            }
        }
    } else {
        stop ("Either X or kinshipfile must be provided")
    }
    return(kinship)
}
