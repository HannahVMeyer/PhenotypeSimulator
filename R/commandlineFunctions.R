#' Command line execution for PhenotypeSimulator.
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
#' #--tNrSNPs=10000 --cNrSNPs=30 \
#' #--SNPfrequencies=0.05,0.1,0.3,0.4 \
#' #--genVar=0.4 --h2s=0.025 --phi=0.6 --delta=0.3 --gamma=1 \
#' #--pcorr=0.8 \
#' #--NrFixedEffects=4 --NrConfounders=1,2,1,2 \
#' #--pIndependentConfounders=0,1,1,0.5 \
#' #--distConfounders=bin,cat_norm,cat_unif,norm \
#' #--probConfounders=0.2 \
#' #--catConfounders=0,3,4,0 \
#' #--directory=/tmp \
#' #--showProgress \


simulatePhenotypes <- function() {
    option_list <- list(
        make_option(c("--NrSamples"), action="store", dest="NrSamples", 
                    type="integer", help="Number of samples to 
                    simulate [default: %default]."),
        make_option(c("--NrPhenotypes"), action="store", dest="NrPhenotypes", 
                    type="integer", help="Number of phenotypes to 
                    simulate [default: %default]."),

        make_option(c("--genVar"), action="store", dest="genVar", 
                    default=NULL, type="double", help="Total genetic variance 
                    [default: %default]."),
        make_option(c("--noiseVar"), action="store", 
                    dest="noiseVar", default=NULL, type="double", 
                    help="Total noise variance [default: %default]."),
        make_option(c("--h2s"), action="store", dest="h2s", default=NULL
                    , type="double", help="Proportion of genetic variant effects 
                    from total genetic variance [default: %default]."),
        make_option(c("--h2bg"), action="store", dest="h2bg", 
                    default=NULL, type="double", help="Proportion of 
                    infinitesimal genetic effects from total genetic variance 
                    [default: %default]."),
        make_option(c("--theta"), action="store", dest="theta", 
                    default=0.8, type="double", help="Proportion of shared 
                    genetic variant effects from total genetic variant variance
                    [default: %default]."),
        make_option(c("--eta"), action="store", dest="eta", default=0.8, 
                    type="double", help="Proportion of shared infinitesimal 
                    genetic effects from total infinitesimal genetic variance  
                    [default: %default]."),
        make_option(c("--delta"), action="store", dest="delta", 
                    default=NULL, type="double", help="Proportion of variance of 
                    confounder effects [default: %default]."),
        make_option(c("--gamma"), action="store", dest="gamma", 
                    default=0.8, type="double", help="Proportion of shared 
                    confounder variance from total confounder variance 
                    [default: %default]."),
        make_option(c("--alpha"), action="store", dest="alpha", 
                    default=0.8, type="double", help="Proportion of shared
                    obervational noise variance from total observational noise 
                    variance [default: %default]."),
        make_option(c("--rho"), action="store", dest="rho", default=NULL
                    , type="double", help="Proportion of correlated noise 
                    effects from total noise variance [default: %default]."),
        make_option(c("--phi"), action="store", dest="phi", default=NULL
                    , type="double", help="Proportion of observational noise 
                    variance from total noise variance [default: %default]."),

        make_option(c("--tNrSNP"), action="store", dest="tNrSNP", 
                    default=5000, type="integer", help="Total number of SNPs to 
                    simulate [default: %default]."),
        make_option(c("--cNrSNP"), action="store", dest="cNrSNP", 
                    default=20, type="integer", help="Number of causal SNPs to 
                    draw from total number of SNPs [default: %default]."),
        make_option(c("--SNPfrequencies"), action="store", 
                    dest="SNPfrequencyString", default="0.1,0.2,0.4", 
                    type="character", help="Comma-separated list of allele 
                    frequencies from which to sample when simulating genotypes 
                    [default: %default]."),

        make_option(c("--genotypefile"), action="store", 
                    dest="genotypefile", default=NULL, type="character", 
                    help="Path to external genotype file (to be fully read into 
                    memory) in format specified by --format 
                    [default: %default]."),
        make_option(c("--format"), action="store", 
                    dest="format", default='delim', type="character", 
                    help="Needed when --genotypefile or
                    --genoFilePrefix/--genoFileSuffix are specified; specifies  
                    format of the genotype data; if --genotypefile: has to be 
                    one of plink, oxgen, genome, bimbam or delim, if 
                    --genoFilePrefix/--genoFileSuffix has to be 
                    one of oxgen, bimbam or delim [default: %default]"),

        make_option(c("--genoFilePrefix"), action="store", 
                    dest="genoFilePrefix", default=NULL, type="character", 
                    help="Full path to file (no tilde-expansion) and prefix of 
                    per-chromosome comma-separated genotypes file e.g. 
                    '/tmp/genotypes_' [default: %default]."),
        make_option(c("--genoFileSuffix"), action="store", 
                    dest="genoFileSuffix", default=NULL, type="character", 
                    help="Optional string of genotype file ending including 
                    format indication (e.g. '.csv') [default: %default]."),
        make_option(c("--genoFileDelimiter"), action="store", 
                    dest="genoFileDelimiter", default=",", type="character", 
                    help="Field separator of --genotypefile or
                    --genoFilePrefix/--genoFileSuffix; if file is 
                    tab-separated, please specify 'tab' [default: %default]."),
        make_option(c("--probabilities"), action="store_true", 
                    dest="probabilities", default=FALSE, type="logical", 
                    help="Set this flag if the genotypes in genoFilePrefix-
                    genoFileSuffix are provided as triplets of probablities 
                    (p(AA), p(Aa), p(aa)); they will be converted into their 
                    expected genotype frequencies by 0*p(AA) + p(Aa) + 2*p(aa) 
                    [default: %default]."),
        make_option(c("--skipFields"), action="store", 
                    dest="skipFields", default=NULL, type="integer", 
                    help="Number of fields (columns) to skip in genoFilePrefix-
                    genoFileSuffix file [default: %default]."),
        make_option(c("--genoFileHeader"), action="store_true", 
                    dest="header", default=FALSE,
                    help="Use flag to indicate that genoFilePrefix-
                    genoFileSuffix file has a header for format == 'delim'. 
                    [default: %default]."),
        make_option(c("--chr"), action="store", 
                    dest="chr_string", default=NULL, type="character", 
                    help="Comma-separated list of chromosomes to draw causal 
                    SNPs from [default: %default]."),
        make_option(c("--NrSNPsOnChromosome"), action="store", 
                    dest="NrSNPsOnChromosomeString", default=NULL, 
                    type="character", help="Comma-separated list of the number 
                    of SNPs per chromosome specified in --chr (see above); has 
                    to be the same length as --chr. If not provided, lines in 
                    file will be counted (which can be slow for large files) 
                    [default: %default]."),
        make_option(c("--NrChrCausal"), action="store", 
                    dest="NrChrCausal", default=NULL, type="integer", 
                    help="Number of chromosomes to randomly draw causal SNPs 
                    from (as opposed to a specific list of chromosomes 
                    to draw from via --chr) [default: %default]."),

        make_option(c("--kinshipFile"), action="store", 
                    dest="kinshipfile", default=NULL, type="character", 
                    help="Path to pre-computed, comma-separated kinshipfile 
                    [default: %default]."),
        make_option(c("--kinshipFileHasNoHeader"), action="store_false", 
                    dest="kinshipHeader", default=TRUE, type="logical", 
                    help="Set to specify that kinship file does not have a 
                    header line [default: header=sample IDs]."),
        make_option(c("--kinshipFileDelimiter"), action="store", 
                    dest="kinshipFileDelimiter", default=",", type="character", 
                    help="Field separator of kinship file (e.g. `,`); if 
                    kinship file is tab-separated, please specify 'tab'
                    [default: %default]."),
        make_option(c("--noStandardise"), action="store_false", 
                    dest="standardise", default=TRUE, type="logical", 
                    help="If set, genotypes will not be standardised for 
                    kinship estimation (standardissation is recommended). 
                    [default: %default]."),
        make_option(c("--pIndependentGenetic"), 
                    action="store", 
                    dest="pIndependentGenetic", default=0.4, type="double", 
                    help="Proportion of variance of genetic variant effects to 
                    have a trait-independent effect [default: %default]."),
        make_option(c("--pTraitIndependentGenetic"), action="store", 
                    dest="pTraitIndependentGenetic", default=0.2, type="double", 
                    help="Proportion of traits influenced by independent 
                    genetic variant effects [default: %default]."),
        make_option(c("--keepSameIndependentSNPs"), 
                    action="store_true", dest="keepSameIndependentSNPs", 
                    default=FALSE, type="logical", help="If this flag is set, 
                    the independent genetic variant effects always influence the 
                    same subset of traits [default: %default]."),
        make_option(c("--pTraitsAffectedGenetics"), action="store", 
                    dest="pTraitsAffectedGeneticsString", default=1, 
                    type="double", help="Proportion of traits affected by 
                    the genetic variant effects [default: %default]."),
        make_option(c("--distBetaGenetic"), action="store", 
                    dest="distBetaGenetic", default="norm", 
                    type="character", help="Distribution to use for 
                    simulating the effect sizes of the genetic variant effects; 
                    one of 'unif' or 'norm' [default: %default]."),
        make_option(c("--mBetaGenetic"), action="store", 
                    dest="mBetaGenetic", default=0, type="double", 
                    help="Mean/midpoint of normal/uniform distribution for 
                    effect sizes of genetic variant effects [default: 
                    %default]."),
        make_option(c("--sdBetaGenetic"), action="store", 
                    dest="sdBetaGenetic", default=1, type="double", 
                    help="Standard deviation/distance from midpoint of 
                    normal/uniform distribution for effect sizes of genetic 
                    variant effects [default: %default]."),

        make_option(c("--NrFixedEffects"), action="store", 
                    dest="NrFixedEffects", default=1, type="integer", 
                    help="Number of different confounder effects to simulate: 
                    allows to simulate confounder effects from different 
                    distributions or with different parameters 
                    [default: %default]."),        
        make_option(c("--NrConfounders"), action="store",
                    dest="NrConfoundersString", default=10, type="character",
                    help="Number of confounders to simulate per set 
                    (default 1 set of confounder effects --NrFixedEffects=1) 
                    [default: %default]."),

        make_option(c("--pTraitsAffectedConfounders"), action="store", 
                    dest="pTraitsAffectedConfoundersString", default="1", 
                    type="character", help="Proportion(s) of traits affected by 
                    the confounder(s); for more than one type of confounders (
                    i.e. --NrFixedEffects > 1), provide values separated by 
                    commas, e.g. '0.2,0.5,1', [default: %default]."),
        make_option(c("--pIndependentConfounders"),
                    action="store", dest="pIndependentConfoundersString", 
                    default=0.4, type="character", 
                    help="Proportion(s) of variance of confounder effects to 
                    have a trait-independent effect; for more than one type of 
                    confounders, provide values separated by commas, e.g. 
                    '0.3,0.4,0.8' [default: %default]."),
        make_option(c("--pTraitIndependentConfounders"), action="store", 
                    dest="pTraitIndependentConfoundersString", default=0.2, 
                    type="character", help="Proportion(s) of traits influenced 
                    by independent confounder effects; for more than one type of 
                    confounders, provide values separated by commas, e.g. 
                    '0.8,0.5,0.9' [default: %default]."),
        make_option(c("--keepSameIndependentConfounders"), 
                    action="store", dest="keepSameIndependentConfoundersString", 
                    default="FALSE", type="character", 
                    help="Specifies if the independent confounder effects should 
                    always influence; for more than one set of confounders, 
                    provide a comma-separated list if TRUE and FALSE, e.g.
                    'TRUE,TRUE,FALSE' [default: %default]."),
        make_option(c("--distConfounders"), action="store", 
                    dest="distConfoundersString", default="norm", 
                    type="character", help="Distribution(s) for simulating 
                    confounders; one of 'unif', 'norm', 'bin', 'cat_norm', 
                    'cat_unif'; for more than one type of confounders, provide
                    values separated by commas, e.g. 'norm,cat_unif', 
                    [default: %default]."),
        make_option(c("--mConfounders"), action="store", 
                    dest="mConfoundersString", default="0", type="character", 
                    help="Mean/midpoint(s) of normal/uniform distribution for 
                    confounders; for more than one type of confounders, provide
                    values separated by commas, e.g. '0,2,1', 
                    [default: %default]."),
        make_option(c("--sdConfounders"), action="store", 
                    dest="sdConfoundersString", default="1", type="character", 
                    help="Standard deviation(s)/distance from midpoint(s) of 
                    normal/uniform distribution for confounders; for more than
                    one type of confounders, provide values separated by commas, 
                    e.g. '1,2,1', [default: %default]."),
        make_option(c("--catConfounders"), action="store", 
                    dest="catConfoundersString", default=NULL, type="character", 
                    help="Number(s) of confounder categories; required if 
                    distConfounders 'cat_norm' or 'cat_unif'; for more than one 
                    type of confounders, provide values separated by commas, 
                    e.g. '2,0,5' (second confounder not categorical),
                    [default: %default]."),
        make_option(c("--probConfounders"), action="store", 
                    dest="probConfoundersString", default=NULL, type="character", 
                    help="Probability(s) of binomial confounders (0/1); required 
                    if distConfounders 'bin', for more than one 
                    type of confounders, provide values separated by commas, 
                    e.g. '0.2,0.6,0' (third confounder not binomial), 
                    [default: %default]."),
        make_option(c("--distBetaConfounders"), action="store", 
                    dest="distBetaConfoundersString", default="norm", 
                    type="character", help="Name(s) of distribution to use to 
                    simulate effect sizes of confounders; one of 'unif' or 
                    'norm'; for more than one type of confounders, provide 
                    values separated by commas, e.g. 'norm,unif', 
                    [default: %default]."),
        make_option(c("--mBetaConfounders"), action="store", 
                    dest="mBetaConfoundersString", 
                    default=0, type="character", help="Mean/midpoint of normal
                    /uniform distribution for effect sizes of confounders;  for 
                    more than one type of confounders, provide values separated 
                    by commas, e.g. '0,0.5', [default: %default]."),
        make_option(c("--sdBetaConfounders"), action="store", 
                    dest="sdBetaConfoundersString", 
                    default=1, type="character", help="Standard deviation/
                    distance from midpoint of normal/uniform distribution for 
                    effect sizes of confounders;  for more than one type of 
                    confounders, provide values separated by commas, e.g. 
                    '1,2', [default: %default]"),

        make_option(c("--pcorr"), action="store", dest="pcorr", 
                    default=0.8, type="double", help="Correlation strength of
                    correlated noise effects [default: %default]."),
        make_option(c("--corrmatFile"), action="store", 
                    dest="corrmatfile", default=NULL, type="character", 
                    help="path/to/corrmatfile.csv [string] with comma-separated
                    [P x P] numeric [double] correlation matrix; if provided,  
                    correlation matrix for simulation of correlated backgound 
                    effect will be read from file; file should NOT contain an 
                    index or header column [default: %default]."),

        make_option(c("--meanNoiseBg"), action="store", dest="meanNoiseBg", 
                    default=0, type="double", help="Mean of the normal 
                    distribution for simulating observational noise effects 
                    [default: %default]."),
        make_option(c("--sdNoiseBg"), action="store", dest="sdNoiseBg", 
                    default=1, type="double", help="Standard deviation of the 
                    normal distribution for simulating observational noise
                    effects [default: %default]."),

        make_option(c("--sampleID"), action="store", dest="sampleID", 
                    default="ID_", type="character", help="Prefix for naming 
                    simulated samples; will be followed by sample number from 
                    1 to --NrSamples when constructing sample IDs; only used if
                    genotypes/kinship are simulated or provided data does not 
                    contain sample IDs [default: %default]."),
        make_option(c("--phenoID"), action="store", dest="phenoID", 
                    default="Trait_", type="character", help="Prefix for naming 
                    simulated phenotypes; will be followed by trait number from 
                    1 to --NrPhenotypes when constructing trait IDs
                    [default: %default]."),
        make_option(c("--snpID"), action="store", dest="snpID", 
                    default="SNP_", type="character", help="Prefix for naming 
                    simulated snps; will be followed by SNP number from 1 to 
                    --tNrSNPs when constructing SNP IDs [default: %default]."),

        make_option(c("--seed"), action="store", dest="seed", 
                    default=219453, type="integer", help="Seed to initialise 
                    random number generator [default: %default]."),
        make_option(c("--showProgress"), action="store_true", dest="verbose", 
                    default=FALSE, type="logical", help="If set, progress 
                    messages about simulation steps are printed to standard out
                    [default: %default]."),

        make_option(c("--nonlinear"), action="store", dest="nonlinear", 
                    default=NULL, type="character", 
                    help="Nonlinear transformation method; one of exp, log, or 
                    sqrt; if log, base can be specified; non-linear 
                    transformation is optional, default is NULL ie no 
                    transformation"),
        make_option(c("--logbase"), action="store", dest="logbase", default=10,
                    type="integer", help="Base of logarithm for non-linear phenotype
                    transformation"),
        make_option(c("--expbase"), action="store", dest="expbase", 
                    default=NULL,
                    type="double", help="Base of exponential function for
                    non-linear phenotype transformation; if non given, Euler's
                    number is used (exp)"),
        make_option(c("--power"), action="store", dest="power", default=2,
                    type="double", help="Power of polynomial non-linear phenotype
                    transformation"),
        make_option(c("--proportionNonlinear"), action="store", default=0,
                    dest="proportionNonlinear", type="double", help="proportion
                    of the phenotype to be non-linear"),
        make_option(c("--transfromNegNonlinear"), action="store", default='abs',
                    dest="transformNeg", type="character", help="transformation
                    method for negative values in non linear phenotype
                    transformation. One of abs (absolute value) or set0 (set all
                    negative values to zero). If nonlinear==log and 
                    transformNegNonlinear==set0, negative values set to 1e-5."),

        make_option(c("--directory"), action="store", 
                    dest="directory", default=NULL, type="character", help=
                    "Absolute path (no tilde expansion) to parent directory
                    where simulation results should be saved; needs user 
                    writing permission [default: %default]."),
        make_option(c("--subdirectory"), action="store", dest="outstring"
                    , default=NULL, type="character", help="Optional name of 
                    subdirectory (in relation to --directory) to save set-up 
                    dependent simulation results; if not specified, subdirectory 
                    named by NrSamples, NrSNPs, genetic Model, noise Model and 
                    genVar is created; if no subdirectory is required specify
                    '' [default: %default]."),

        make_option(c("--saveTable"), action="store_true", 
                    dest="saveAsTable", default=FALSE, type="logical", 
                    help="Output format of results: when flag set, output saved 
                    as .csv; at least one of --saveTable or --saveRDS needs to 
                    be set [default: %default]."),
        make_option(c("--saveRDS"), action="store_true", 
                    dest="saveAsRDS", default=FALSE, type="logical", 
                    help="Output format of results: when flag set, output saved 
                    as .rds;  at least one of -saveTable or -saveRDS needs to be 
                    set [default: %default]."),
        make_option(c("--saveLIMMBO"), action="store_true", 
                    dest="saveAsLIMMBO", default=FALSE, type="logical",
                    help="When flag is set, simulated data is saved in
                    LiMMBo format, i.e. Ysim_limmbo.csv (phenotype file),
                    Covs_limmbo.csv (covariates file), Kinship_limmbo.csv
                    (kinship file) and genotypes_limmbo.csv (genotype file)
                    [default: %default]."),
        make_option(c("--savePLINK"), action="store_true", 
                    dest="saveAsPLINK", default=FALSE, type="logical", 
                    help="When flag is set, simulated data is saved in  
                    binary plink format, i.e. .bed, .bim and .fam (genotypes) 
                    files, Ysim_plink.txt (phenotypes) and Covs_plink.txt (
                    covariates) [default: %default]."),
        make_option(c("--saveGEMMA"), action="store_true", 
                    dest="saveAsGEMMA", default=FALSE, type="logical", 
                    help="When flag is set, simulated data is saved in
                    gemma format, i.e. Ysim_gemma.txt (phenotype file), 
                    Covs_gemma.txt (covariates file), Kinship_gemma.txt 
                    (kinship file) and genotypes.gemma (genotype file)
                    [default: %default]."),
        make_option(c("--noGemmaIntercept"),
                    action="store_false", dest="intercept_gemma", default=TRUE, 
                    type="logical", help ="if --saveGEMMA: when modeling an 
                    intercept term in gemma, a column of 1's has to be appended 
                    to the covariate files. Set noGemmaIntercept when appending 
                    a column of 1's is not desired [default: %default"),
        make_option(c("--saveBIMBAM"), action="store_true", 
                    dest="saveAsBIMBAM", default=FALSE, type="logical", 
                    help="When flag is set, simulated genotypes are saved in the 
                    BIMBAM format, i.e. Ysim_bimbam.txt (phenotype file) and 
                    genotypes.bimbam (genotype file) [default: %default]."),
        make_option(c("--saveSNPTEST"), action="store_true", 
                    dest="saveAsSNPTEST", default=FALSE, type="logical", 
                    help="When flag is set, simulated genotypes are saved in  
                    snptest format, i.e.  Ysim_snptest.sample (phenotype and 
                    covariate file) and genotypes_snptest.gen (genotype file)
                    [default: %default]."),
        make_option(c("--dontSaveIntermediate"), 
                    action="store_false", dest="saveIntermediate", 
                    default=TRUE, 
                    type="logical", help="Set to specify that intermediate 
                    phenotype components such as shared and independent effects 
                    components should not be saved; [default: %default]")
        )
    args <- parse_args(OptionParser(option_list=option_list))
    if(!args$NrSamples || !args$NrPhenotypes) {
        stop("Number of samples and number of phenotypes have to be specified 
             via --NrSamples and --NrPhenotypes; 
             if you just want usage information please use the --help flag")
    }
    if(!args$saveAsRDS && !args$saveAsTable) {
        stop("At least one of --saveRDS or --saveTable need to be set")
    }

    if (grepl("~", args$directory)) {
        stop("directory contains ~: path expansion not guaranteed on 
             every platform (see path.expand{base}), please provide full file
             path to the output directory")
    }
    NrConfounders <- commaList2vector(args$NrConfoundersString)
    SNPfrequencies <- commaList2vector(args$SNPfrequencyString)
    pIndependentConfounders <- 
        commaList2vector(args$pIndependentConfoundersString)
    pTraitIndependentConfounders <-
        commaList2vector(args$pTraitIndependentConfoundersString)
    pTraitsAffectedConfounders <-
        commaList2vector(args$pTraitsAffectedConfoundersString)
    keepSameIndependentConfounders <-
        commaList2vector(args$keepSameIndependentConfoundersString,
                                            type="logical")
    distConfounders <- commaList2vector(args$distConfoundersString,
                                        type="character")
    mConfounders <- commaList2vector(args$mConfoundersString)
    sdConfounders <- commaList2vector(args$sdConfoundersString)
    catConfounders <- commaList2vector(args$catConfoundersString)
    probConfounders <- commaList2vector(args$probConfoundersString)
    distBetaConfounders <- commaList2vector(args$distBetaConfoundersString,
                                            type="character")
    mBetaConfounders <- commaList2vector(args$mBetaConfoundersString)
    sdBetaConfounders <- commaList2vector(args$sdBetaConfoundersString)

    chr <- commaList2vector(args$chr_string)
    NrSNPsOnChromosome <- commaList2vector(args$NrSNPsOnChromosomeString)

    if (args$verbose) {
        message("Output directory: ", args$directory)
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
        message("Number of phenotypes: ", args$NrPhenotypes)
        message("Number of samples: ", args$NrSamples)
        message("Number of causal SNPs: ", args$cNrSNP)
        message("Number of confounders: ", sum(NrConfounders))
    }

    format <- NULL
    if (args$saveAsRDS) format <- c(format, "rds")
    if (args$saveAsTable) format <- c(format, "csv")
    if (args$saveAsPLINK) format <- c(format, "plink")
    if (args$saveAsGEMMA) format <- c(format, "gemma")
    if (args$saveAsBIMBAM) format <- c(format, "bimbam")
    if (args$saveAsSNPTEST) format <- c(format, "snptest")
    if (args$saveAsLIMMBO) format <- c(format, "limmbo")

    if(tolower(args$kinshipFileDelimiter) == "tab") {
        args$kinshipFileDelimiter="\t"
    }
    if(tolower(args$genoFileDelimiter) == "tab") {
        args$genoFileDelimiter="\t"
    }
    simulatedPheno <- runSimulation( N=args$NrSamples, P=args$NrPhenotypes,
                                     seed=args$seed,
                                     tNrSNP=args$tNrSNP, cNrSNP=args$cNrSNP,
                                     pTraitsAffectedGenetics=
                                         args$pTraitsAffectedGenetics,
                                     pTraitsAffectedConfounders=
                                         pTraitsAffectedConfounders,
                                     NrConfounders=NrConfounders,
                                     NrFixedEffects=args$NrFixedEffects,
                                     chr=chr,
                                     NrChrCausal=args$NrChrCausal,
                                     NrSNPsOnChromosome=NrSNPsOnChromosome,
                                     probabilities = args$probabilities,
                                     skipFields=args$skipFields,
                                     header=args$header,
                                     sampleID=args$sampleID,
                                     phenoID=args$phenoID,
                                     genoFilePrefix=args$genoFilePrefix,
                                     genoFileSuffix=args$genoFileSuffix,
                                     SNPfrequencies=SNPfrequencies,
                                     genoDelimiter=args$genoFileDelimiter,
                                     kinshipfile=args$kinshipfile,
                                     kinshipHeader=args$kinshipHeader,
                                     kinshipDelimiter=args$kinshipFileDelimiter,
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
                                     corrmatfile=args$corrmatfile,
                                     pIndependentConfounders=
                                         pIndependentConfounders,
                                     pTraitIndependentConfounders=
                                        pTraitIndependentConfounders,
                                     keepSameIndependentConfounders=
                                         keepSameIndependentConfounders,
                                     pIndependentGenetic=
                                         args$pIndependentGenetic,
                                     pTraitIndependentGenetic=
                                         args$pTraitIndependentGenetic,
                                     keepSameIndependentSNPs=
                                         args$keepSameIndependentSNPs,
                                     distBetaGenetic=args$distBetaGenetic,
                                     mBetaGenetic=args$mBetaGenetic,
                                     sdBetaGenetic=
                                         args$sdBetaGenetic,
                                     distConfounders=distConfounders,
                                     mConfounders=mConfounders,
                                     sdConfounders=sdConfounders,
                                     catConfounders=catConfounders,
                                     probConfounders=probConfounders,
                                     distBetaConfounders=distBetaConfounders,
                                     mBetaConfounders=mBetaConfounders,
                                     sdBetaConfounders=sdBetaConfounders,
                                     meanNoiseBg=args$meanNoiseBg,
                                     sdNoiseBg=args$sdNoiseBg,
                                     nonlinear=args$nonlinear,
                                     logbase=args$logbase,
                                     expbase=args$expbase,
                                     power=args$power,
                                     transformNeg=
                                         args$transformNegNonlinear,
                                     proportionNonlinear=
                                         args$proportionNonlinear,
                                     verbose=args$verbose)

    outdir <- savePheno(simulatedPheno,
                        format=format,
                        saveIntermediate=args$saveIntermediate,
                        intercept_gemma=args$intercept_gemma,
                        outstring=args$outstring,
                        directory=args$directory,
                        verbose=args$verbose)
}
