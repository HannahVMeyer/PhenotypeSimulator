#' Command line execution for PhenotypeSimulator
#' 
#' simulatePhenotypes runs without arguments. Upon call, it reads command-line
#' parameters and supplies these to \code{\link{runSimulation}} and 
#' \code{\link{savePheno}}. For details on input to \code{\link{runSimulation}} 
#' and \code{\link{savePheno}}, please refer to their help pages. For help on 
#' the command line arguments that can be passed, see examples below. From the 
#' command line, the help function can be called via 
#' 'Rscript -e "PhenotypeSimulator::simulatePhenotypes()" --args --help 
#'  
#' @export
#' 
#' @examples
#' # (not run)
#' # Simulate simple phenotype of genetic and noise background effects only:
#' # (not run)
#' # Rscript -e "PhenotypeSimulator::simulatePhenotypes()" \
#' #--args \ 
#' #--NrSamples=100 --NrPhenotypes=15 \
#' #--tNrSNP=10000 --cNrSNP=30 \
#' #--SNPfrequencies=0.05,0.1,0.3,0.4 \
#' #--genVar=0.4 --h2s=0.025 --phi=0.6 --delta=0.3 --gamma=1 \
#' #--pcorr=0.8 \
#' #--NrFixedEffects=4 --NrConfounders=1,2,1,2 \
#' #--pIndependentConfounders=0,1,1,0.5 \
#' #--distConfounders=bin,cat_norm,cat_unif,norm \
#' #--probConfounders=0.2 \
#' #--catConfounders=0,3,4,0 \
#' #--directoryGeno=~/tmp/genotypes \
#' #--directoryPheno=~/tmp/phenotypes \
#' #--sampleSubset=50,70 \
#' #--phenoSubset=5,10 \
#' #--normalise \
#' #--showProgress \
#' #--saveTable \
#' #--savePlink


simulatePhenotypes <- function() {
    option_list <- list(
        make_option(c("-P", "--NrPhenotypes"), action="store", dest="P", 
                    default=10, type="integer", help="Number of phenotypes to 
                    simulate [default: %default]"),
        make_option(c("-N", "--NrSamples"), action="store", dest="N", 
                    default=1000, type="integer", help="Number of samples to 
                    simulate [default: %default]"),
        make_option(c("-C", "--NrConfounders"), action="store",
                    dest="NrConfoundersString", default=10, type="character",
                    help="Number of confounders (fixed noise effects) per set 
                    (default 1 set of fixedEffects) to simulate 
                    [default: %default]"),
        make_option(c("-F", "--NrFixedEffects"), action="store", 
                    dest="NrFixedEffects", default=1, type="integer", 
                    help="Number of different fixed noise effects to simulate, 
                    allows to simulate fixed effects from different 
                    distributions or with different parameters 
                    [default: %default]"),
        make_option(c("-tS", "--tNrSNPs"), action="store", dest="tNrSNPs", 
                    default=5000, type="integer", help="Number of total SNPs to 
                    simulate [default: %default]"),
        make_option(c("-cS", "--cNrSNPs"), action="store", dest="cNrSNPs", 
                    default=20, type="integer", help="Number of causal SNPs to 
                    draw from total SNPs [default: %default]"),
        
        make_option(c("-SNPfrequencies", "--SNPfrequencies"), action="store", 
                    dest="SNPfrequencyString", default="0.1,0.2,0.4", 
                    type="character", help="Comma-separated list of allele 
                    frequencies from which to sample to simulate genotypes 
                    [default: %default]"),
        make_option(c("-psg", "--pIndependentGenetic"), action="store", 
                    dest="pIndependentGenetic", default=0.4, type="double", 
                    help="Proportion of variance of fixed genetic effects to 
                    have a trait-independent effect [default: %default]"),
        make_option(c("-ptsg", "--pTraitIndependentGenetic"), action="store", 
                    dest="pTraitIndependentGenetic", default=0.2, type="double", 
                    help="Proportion of traits influenced by independent fixed 
                    genetic effects [default: %default]"),
        make_option(c("-psn", "--pIndependentConfounders"), action="store", 
                    dest="pIndependentConfounders", default=0.4, 
                    type="character", 
                    help="Proportion of variance of fixed noise effects to have 
                    a trait-independent effect [default: %default]"),
        make_option(c("-ptsn", "--pTraitIndependentConfounders"),
                    action="store", 
                    dest="pTraitIndependentConfounders", default=0.2, 
                    type="character", help="Proportion of traits influenced by 
                    independent fixed noise effects [default: %default]"),
        
        make_option(c("-genVar", "--genVar"), action="store", dest="genVar", 
                    default=NULL, type="double", help="Total genetic variance 
                    [default: %default]"),
        make_option(c("-noiseVar", "--noiseVar"), action="store", 
                    dest="noiseVar", default=NULL, type="double", 
                    help="Total noise variance [default: %default]"),
        make_option(c("-h2s", "--h2s"), action="store", dest="h2s", default=NULL
                    , type="double", help="Proportion of genetic variance of 
                    fixed genetic effects [default: %default]"),
        make_option(c("-h2bg", "--h2bg"), action="store", dest="h2bg", 
                    default=NULL, type="double", help="Proportion of genetic 
                    variance of background genetic effects [default: %default]"),
        make_option(c("-theta", "--theta"), action="store", dest="theta", 
                    default=0.8, type="double", help="Proportion of genetic 
                    variance of shared fixed genetic effects [default: %default]"
        ),
        make_option(c("-eta", "--eta"), action="store", dest="eta", default=0.8, 
                    type="double", help="Proportion of genetic variance of 
                    shared genetic background effects  [default: %default]"),
        
        make_option(c("-delta", "--delta"), action="store", dest="delta", 
                    default=NULL, type="double", help="Proportion of variance of 
                    fixed noise effects [default: %default]"),
        make_option(c("-gamma", "--gamma"), action="store", dest="gamma", 
                    default=0.8, type="double", help="Proportion of noise 
                    variance of shared fixed noise effects [default: %default]"),
        make_option(c("-alpha", "--alpha"), action="store", dest="alpha", 
                    default=0.8, type="double", help="Proportion of noise 
                    variance of shared background noise effect 
                    [default: %default]"),
        make_option(c("-rho", "--rho"), action="store", dest="rho", default=NULL
                    , type="double", help="Proportion of noise variance of 
                    correlated noise effects [default: %default]"),
        make_option(c("-phi", "--phi"), action="store", dest="phi", default=NULL
                    , type="double", help="Proportion of noise variance of 
                    background noise effects [default: %default]"),
        make_option(c("-pcorr", "--pcorr"), action="store", dest="pcorr", 
                    default=0.6, type="double", help="Correlation strength of
                    correlated noise effects [default: %default]"),
        make_option(c("--distConfounders"), action="store", 
                    dest="distConfoundersString", default="norm", 
                    type="character", help="[distribution to simulate 
                    confounders; one of 'unif', 'norm', 'bin', 'cat_norm', 
                    'cat_unif', default: %default]"),
        make_option(c("--mConfounders"), action="store", 
                    dest="mConfounderString", default="0", type="character", 
                    help="Mean/midpoint(s) of normal/uniform distribution for 
                    confounders [default: %default]"),
        make_option(c("--sdConfounders"), action="store", 
                    dest="sdConfounderString", default="1", type="character", 
                    help="standard deviation(s)/distance from midpoint(s) of 
                    normal/uniform distribution for confounders 
                    [default: %default]"),
        make_option(c("--catConfounders"), action="store", 
                    dest="catConfoundersString", default=NULL, type="character", 
                    help="number(s) of confounder categories; required if 
                    distConfounders 'cat_norm' or 'cat_unif' 
                    [default: %default]"),
        make_option(c("--probConfounders"), action="store", 
                    dest="probConfoundersString", default=0, type="character", 
                    help="probability(s) of binomial confounders (0/1); required 
                    if distConfounders 'bin' [default: %default]"),
        make_option(c("--distBetaConfounders"), action="store", 
                    dest="distBetaConfoundersString", default="norm", 
                    type="character", help="Name(s) of distribution to use to 
                    simulate effect sizes of confounders; one of 'unif' or 
                    'norm' [default: %default]"),
        make_option(c("--mBetaConfounders"), action="store", 
                    dest="mBetaConfoundersString", 
                    default=0, type="character", help="Mean/midpoint of normal
                    /uniform distribution for effect sizes of confounders 
                    [default: %default]"),
        make_option(c("--sdBetaConfounders"), action="store", 
                    dest="sdBetaConfoundersString", 
                    default=1, type="character", help="Standard deviation/
                    distance from midpoint of normal/uniform distribution for 
                    effect sizes of confounders [default: %default]"),
        make_option(c("--distBetaGenetic"), action="store", 
                    dest="distBetaGeneticString", default="norm", 
                    type="character", help="Name(s) of distribution to use to 
                    simulate effect sizes of SNPs; one of 'unif' or 'norm' 
                    [default: %default]"),
        make_option(c("--mBetaGenetic"), action="store", 
                    dest="mBetaGeneticString", default=0, type="character", 
                    help="Mean/midpoint of normal/uniform distribution for 
                    effect sizes of SNPs [default: %default]"),
        make_option(c("--sdBetaGenetic"), action="store", 
                    dest="sdBetaGeneticString", default=1, type="character", 
                    help="Standard deviation/distance from midpoint of 
                    normal/uniform distribution for effect sizes of SNPs 
                    [default: %default]"),
        
        make_option(c("--meanNoiseBg"), action="store", dest="meanNoiseBg", 
                    default=0, type="double", help="Mean of the normal 
                    distribution [default: %default]"),
        make_option(c("--sdNoiseBg"), action="store", dest="sdNoiseBg", 
                    default=1, type="double", help="Standard deviation of the 
                    normal distribution [default: %default]"),
        
        make_option(c("--seed"), action="store", dest="seed", 
                    default=219453, type="integer", help="Seed to initialise 
                    random number generator [default: %default]"),
        make_option(c("--showProgress"), action="store_true", dest="verbose", 
                    default=FALSE, type="logical", help=" [default: %default]"),
        make_option(c("-norm", "--normalise"), action="store_true", 
                    default=FALSE,
                    dest="normalise", type="logical", help="Should user-supplied 
                    kinship be normalised [default: %default]"),
        make_option(c("-stand", "--standardise"), action="store_true", 
                    dest="standardise", default=FALSE, type="logical", 
                    help="Should genotypes be standardised for simulation of 
                    fixed effects [default: %default]"),
        
        make_option(c("-chrom", "--chromosomes"), action="store", 
                    dest="chr_string", default=NULL, type="character", 
                    help="Comma-separated list of chromosomes to draw causal 
                    SNPs from [default: %default]"),
        make_option(c("--NrCausalChrom"), action="store", 
                    dest="NrChrCausal", default=1, type="integer", 
                    help="Number of chromosomes to draw causal 
                    SNPs from (as opposed to a independent list of chromosomes 
                    to draw from via --chrom) [default: %default]"),
        
        make_option(c("-dg", "--directoryGeno"), action="store", 
                    dest="directoryGeno", default=NULL, type="character", help=
                    "Absolute path (no tilde expansion) to parent directory
                    for genotypes; [needs user writing permission,
                    [default: %default]"),
        make_option(c("-dp", "--directoryPheno"), action="store", 
                    dest="directoryPheno", default=NULL, type="character", 
                    help="Absolute path (no tilde expansion) to parent directory
                    for phenotypes; needs user writing permission,
                    [default: %default]"),
        make_option(c("-ds", "--subdirectory"), action="store", dest="outstring"
                    , default=NULL, type="character", help="Name of subdirectory
                    to be created within dg/dp [default: %default]"),
        make_option(c("-kf", "--kinshipfile"), action="store", 
                    dest="kinshipfile", default=NULL, type="character", 
                    help="Path to pre-computed, comma-separated kinshipfile 
                    (header=sample IDs) [default: %default]"),
        make_option(c("--kinshipfileHasHeader"), action="store_true", 
                    dest="kinshipHeader", default=FALSE, type="logical", 
                    help="Does kinship have a header line (e.g. header=sample 
                    IDs)? [default: %default]"),
        make_option(c("--kinshipfileDelimiter"), action="store", 
                    dest="kinshipDelimiter", default=",", type="character", 
                    help="Field separator of kinship file (e.g. `,`) 
                    [default: %default]"),
        make_option(c("--genoFilePrefix"), action="store", 
                    dest="genoFilePrefix", default=NULL, type="character", 
                    help="Path to and prefix of per-chromosome comma-separated 
                    genotypes file [default: %default]"),
        make_option(c("--genoFileSuffix"), action="store", 
                    dest="genoFileSuffix", default="", type="character", 
                    help="Optional string of genotype file file ending including 
                    format indication (e.g. '.csv') [default: %default]"),
        make_option(c("--genoFileDelimiter"), action="store", 
                    dest="genoFileDelimiter", default=",", type="character", 
                    help="Field separator of genotype file [default: %default]"),
        make_option(c("--sampleID"), action="store", dest="sampleID", 
                    default="ID_", type="character", help="Prefix for naming 
                    simulated samples [default: %default]"),
        make_option(c("--phenoID"), action="store", dest="phenoID", 
                    default="Trait_", type="character", help="Prefix for naming 
                    simulated phenotypes [default: %default]"),
        
        make_option(c("-sSubset", "--sampleSubset"), action="store", 
                    dest="sample_subset_string", default=NULL, type="character", 
                    help="Comma-separated list of samples sizes to be drawn and 
                    saved from simulation [default: %default]"),
        make_option(c("-pSubset", "--phenoSubset"), action="store", 
                    dest="pheno_subset_string", default=NULL, type="character", 
                    help="Comma-separated list of phenotype sizes to be drawn 
                    and saved from simulation [default: %default]"),
        
        make_option(c("-saveTable", "--saveTable"), action="store_true", 
                    dest="saveAsTable", default=FALSE, type="logical", 
                    help="Output format of results: when flag set, output saved 
                    as .csv; at least one of -saveTable or -saveRDS needs to be 
                    set [default: %default]"),
        make_option(c("-saveRDS", "--saveRDS"), action="store_true", 
                    dest="saveAsRDS", default=FALSE, type="logical", 
                    help="Output format of results: when flag set, output saved 
                    as .rds;  at least one of -saveTable or -saveRDS needs to be 
                    set  [default: %default]"),
        make_option(c("-savePlink", "--savePlink"), action="store_true", 
                    dest="saveAsPlink", default=FALSE, type="logical", 
                    help="When flag is set, simulated genotypes can additionally
                    be saved in the binary plink format, i.e. .bed, .bim and 
                    .fam files. [default: %default]")
        )
    args <- parse_args(OptionParser(option_list=option_list))
    if(!args$saveAsRDS && !args$saveAsTable) {
        stop("At least one of --saveRDS or --saveTable need to be set")
    }
    
    if (args$verbose) {
        message("Phenotype directory: ", args$directoryPheno)
        message("Genotype directory: ", args$directoryGeno)
        message("Subdirectory: ", args$outstring)
        if (!is.null(args$kinshipfile)) {
            message("Kinship file: ", args$kinshipfile)
        }
        if (!is.null(args$genoFilePrefix)) {
            message("SNP file prefix: ", args$genoFilePrefix)
            message("SNP file suffix: ", args$genoFileSuffix)
        } else {
            message(paste("Number of SNPs to simulate (for kinship estimation",
                          "and drawing of causal SNPs): "), args$tNrSNP)
        }
        message("Number of phenotypes: ", args$P)
        message("Number of samples: ", args$N)
        message("Number of causal SNPs: ", args$cNrSNP)
        message("Number of confounders: ", args$NrConfoundersString)
        
        message("Total genetic variance: ", args$genVar)
        message("Total noise variance: ", 1 - args$genVar)
    }
    
    simulatedPheno <- runSimulation( N=args$N, P=args$P, seed=args$seed,
                                     tNrSNP=args$tNrSNP, cNrSNP=args$cNrSNP,
                                     NrConfoundersString=
                                         args$NrConfoundersStrings, 
                                     NrFixedEffects=args$NrFixedEffects,
                                     chr_string=args$chr_string, 
                                     NrChrCausal=args$NrChrCausal,
                                     sampleID=args$sampleID, 
                                     phenoID=args$phenoID,
                                     genoFilePrefix=args$genoFilePrefix, 
                                     genoFileSuffix=args$genoFileSuffix, 
                                     SNPfrequencyString=args$SNPfrequencyString,
                                     genoFileDelimiter=args$genoFileDelimiter,
                                     kinshipfile=args$kinshipfile,
                                     kinshipHeader=args$kinshipHeader,
                                     kinshipDelimiter=args$kinshipDelimiter,
                                     normalise=args$normalise, 
                                     standardise=args$standardise, 
                                     genVar=args$genVar, 
                                     h2s=args$h2s, h2bg=args$h2bg,
                                     theta=args$theta, eta=args$eta,
                                     noiseVar=args$noiseVar, 
                                     delta=args$delta, rho=args$rho, 
                                     phi=args$phi,
                                     alpha=args$alpha,
                                     gamma=args$gamma, 
                                     pcorr=args$pcorr,
                                     pIndependentConfoundersString=
                                         args$pIndependentConfoundersString, 
                                     pTraitIndependentConfoundersString=
                                        args$pTraitIndependentConfoundersString, 
                                     pIndependentGenetic=
                                         args$pIndependentGenetic, 
                                     pTraitIndependentGenetic=
                                         args$pTraitIndependentGenetic, 
                                     distBetaGeneticString=
                                         args$distBetaGeneticString,
                                     mBetaGeneticString=args$mBetaGeneticString, 
                                     sdBetaGeneticString=
                                         args$sdBetaGeneticString,
                                     distConfoundersString=
                                         args$distConfoundersString, 
                                     mConfoundersString=args$mConfoundersString, 
                                     sdConfoundersString=
                                         args$sdConfoundersString,
                                     catConfoundersString=
                                         args$catConfoundersString,
                                     probConfoundersString=
                                         args$probConfoundersString,
                                     distBetaConfoundersString=
                                         args$distBetaConfoundersString,
                                     mBetaConfoundersString=
                                         args$mBetaConfoundersString, 
                                     sdBetaConfoundersString=
                                         args$sdBetaConfoundersString,
                                     meanNoiseBg=args$meanNoiseBg, 
                                     sdNoiseBg=args$sdNoiseBg,
                                     verbose=args$verbose)
    
    outdir <- savePheno(simulatedPheno, 
                            sample_subset_string=args$sample_subset_string, 
                            pheno_subset_string=args$pheno_subset_string,  
                            outstring=args$outstring, 
                            directoryGeno=args$directoryGeno, 
                            directoryPheno=args$directoryPheno,
                            verbose=args$verbose,
                            saveAsTable=args$saveAsTable,
                            saveAsRDS=args$saveAsRDS,
                            saveAsPlink=args$saveAsPlink)
}
