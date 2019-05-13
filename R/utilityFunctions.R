#' Print userinfo.
#'
#' Wrapper function around \code{\link{message}} that allows to turn the 
#' printing of messages to standard out.
#' on or off
#'
#' @param userinfo Vector of [string] element(s) and variables
#' @param verbose [boolean] If TRUE message is displayed on standard out, if 
#' FALSE, message is suppressed.
#' @param sep Delimiter [string] to separate message elements when userinfo 
#' given as vector.
#' @seealso \code{\link{message}} which this function wraps
vmessage <- function(userinfo, verbose=TRUE, sep=" ") {
    if (!is.character(sep)) {
        stop("Separator passed to vmessage must be of type character")
    }
    if (verbose) {
        message(paste(userinfo, collapse=sep))
    }
}

#' Comma-separated string to numeric vector.
#'
#' Split input of comma-separated string into vector of numeric, logical or 
#' character values.
#'
#' @param commastring Input [character] vector containing numbers separated by 
#' commas.
#' @param type Name [string] of the type of input variables; one of numeric, 
#' logical or character
#' @return Numeric vector of values extracted from commastring.
commaList2vector <- function(commastring=NULL, type="numeric") {
    if (is.null(commastring)) {
        return(NULL)
    }
    commastring <- gsub(" ", "", commastring)
    if (type == "numeric") {
        tmp <- as.numeric(unlist(strsplit(commastring, ",")))
    } else if (type == "logical") {
        tmp <- as.logical(unlist(strsplit(commastring, ",")))
    } else if (type == "character") {
        tmp <- unlist(strsplit(commastring, ","))
    } else {
        stop("Unknown type of comma-separated list elements")
    }
    return(tmp)
}

#' Add all non-NULL elements of list.
#'
#' @param compList List of numeric matrices or data.frames of the equal 
#' dimensions.
#' @return Matrix or data.frame containing sum of all list elements where 
#' \code{is.null} is FALSE.
addNonNulls <- function(compList) {
            if (!is.list(compList)) {
                stop("addNonNulls expects input of type list")
            }
            nonNulls <- compList[!sapply(compList, is.null)]
            if (length(nonNulls) == 0) return(NULL)
            if (length(unique(sapply(nonNulls, ncol))) != 1) {
                stop("Column dimensions of list elements are different")
            }
            if (length(unique(sapply(nonNulls, nrow))) != 1) {
                stop("Row dimensions of list elements are different")
            }
            Reduce("+", nonNulls)
}

#' Data simulation for different distributions.
#'
#' Wrapper function to simulate data from different distribution with different 
#' parameter settings.
#'
#' @param x The number [integer] of observations to simulate.
#' @param dist Name of distribution [string] from which the observations are 
#' drawn. 'norm' is the normal distribution, 'unif' the uniform distribution 
#' 'bin' the binomial distribution, "cat_norm" samples categorical variables 
#' according to a normal distribution and "cat_unif" 
#' according to a uniform distribution. For "cat_norm", length(category)/2 is 
#' used mean for the normal distribution unless
#' specified otherwise.
#' @param m Mean of the normal distribution [double]/the mean between min 
#' and max for the uniform distribution [double]/ 
#' the rank of the category to be used as mean for "cat_norm" [integer].
#' @param std Standard deviation of the normal distribution or the distance 
#' of min/max from the mean for the uniform distribution [double].
#' @param categories Number of categories [integer] for simulating categorical 
#' variables (for distr="cat_norm" or "cat_unif").
#' @param prob Probability [double] of success for each trial 
#' (for distr="bin").
#' @return Numeric vector of length [x] with the sampled values
#' @seealso \code{\link{runif}}, \code{\link{rnorm}}, \code{\link{rbinom}} for 
#' documentation of the underlying distributions.
#' @export
#' @examples
#' normal <- simulateDist(x=10, dist="norm", m=2, std=4)
#' cat_normal <- simulateDist(x=10, dist="cat_norm", categories=5)
#' cat_uniform <- simulateDist(x=10, dist="cat_unif", categories=5)
#' uniform <- simulateDist(x=10, dist="unif", m=4, std=1)
#' binomial <- simulateDist(x=10, dist="bin", prob=0.4)
simulateDist <- function(x, 
                         dist=c("unif", "norm", "bin", "cat_norm", "cat_unif"), 
                         m=NULL, std=1, categories=NULL, prob=NULL) {
    if (x < 0) {
        stop(paste("The number of observations to simulate (", x, ") has ",
                   "to be greater than 0"))
    }
    if (length(dist) > 1) {
        stop("Please specify exactly one distribution to sample from,",
             "currently ", length(dist), " provided.")
    }
    if (dist == "unif") {
        if (!is.null(m) && !is.numeric(m)) {
            stop("mean between min and max for the uniform distribution is ", 
                 "not numeric")
        }
        if (!is.numeric(m)) {
            stop("distance of min/max from the mean for the uniform ", 
                 "distribution is not numeric")
        }
        if (is.null(m)) m <- 0
        d <- runif(n=x, min=m - std, max=m + std)
    }
    else if (dist == "norm") {
        if (!is.null(m) && !is.numeric(m)) {
            stop("Mean of normal distribution is not numeric")
        }
        if (!is.numeric(std)) {
            stop("Standard deviation of normal distribution is not numeric")
        }
        if (is.null(m)) m <- 0
        if (std < 0) {
            stop(paste("Simulating normal distribution: standard deviation ",
                       "to be greater than 0"))
        }
        d <- rnorm(n=x, mean=m, sd=std)
    }
    else if (dist == "bin") {
        if (is.null(prob)) {
            stop(paste("Simulating binomial distribution: Probability has to ",
                       "be specified"))
        }
        if (!is.numeric(prob)) {
            stop(paste("Simulating binomial distribution: Probability has to",
                       "be numeric"))
        }
        if (prob < 0) {
            stop(paste("Simulating binomial distribution: Probability has to",
                "be between 0 and 1"))
        }
        d <- rbinom(n=x, size=1, prob=prob)
    }
    else if (grepl("cat", dist)) {
        if (is.null(categories)) {
            stop(paste("Simulating categorical distribution: number of ",
                       "categories has to be specified"))
        }
        if (!is.numeric(categories)) {
            stop(paste("Simulating categorical distribution: categories has to",
                       "be numeric"))
        }
        if (categories <= 0) {
            stop(paste("Simulating categorical distribution: number of ",
                       "categories has to be greater than 0"))
        }
        if (dist == "cat_norm") {
            # generate probabilities for categories
            if (!is.null(m) && !is.numeric(m)) {
                stop("Median of cat_norm distribution is not numeric")
            }
            if (is.null(m)) m <- median(1:categories)
            ptmp=sapply(1:categories, pnorm, mean=m, sd=std)
            prob = c(ptmp[1],diff(ptmp))/sum(c(ptmp[1],diff(ptmp)))
        }
        else if (dist == "cat_unif") {
            prob=rep(1/categories, categories)
        } else {
            stop("Unknown distribution")
        }
        d = sample(1:categories, x, replace=TRUE, prob=prob)
        d <- as.numeric(d)
    } else {
        stop("Unknown distribution")
    }
    return(d)
}

#' Compute expected genotypes from genotype probabilities.
#'
#' Convert genotypes encoded as triplets of probablities (p(AA), p(Aa), p(aa))
#' into their expected genotype frequencies by 0*p(AA) + p(Aa) + 2p(aa).
#'
#' @param probGeno Vector [numeric] with genotype probabilites; has to be a
#' multiple of 3.
#' @return Numeric vector of length [length(probGeno)/3] with the expected 
#' genotype value per individual.
#' @export
#' @examples
#' nrSamples <- 10
#' # Construct genotype probability vector (usually from external input)
#' # First, assign zero probabilty of AA, Aa and aa for all samples
#' genotype_prob <- rep(0, 3*nrSamples)
#' # Second, for each sample draw one of 0,1,2 (corresponding to AA, Aa and aa)
#' genotype_prob[seq(1, nrSamples*3, 3) + sample(0:2, 10, replace=TRUE)] <- 1
#' genotype_exp <- probGen2expGen(genotype_prob)
probGen2expGen <- function(probGeno) {
    if (is.list(probGeno) || !is.vector(probGeno, mode="numeric")) {
        stop("probGen2expGen takes a vector of genotype probabilities as input,
             but ",  typeof(probGeno), " provided.")
    }
    if( length(probGeno) %% 3 != 0) {
        stop("Length of genotype probabilty vector (",  length(probGeno), 
             ") is not a multiple of three")
    }
    if( any(is.na(probGeno))) {
        stop("Samples need to be fully genotyped, but missing genotype data
             detected. Consider imputing missing data, or excluding missing 
             genotypes")
    }
    multiples <- seq(1, length(probGeno), 3)
    testProb <- zoo::rollapply(unlist(probGeno), 3, by = 3, sum, 
                              partial = FALSE)
    if (any(testProb != 1)) {
        stop("Genotype probabilities do not sum to 1")
    }
    probGeno[multiples] <- 0
    probGeno[(multiples + 2)] <- 2 *  probGeno[(multiples + 2)]
    expGeno <- zoo::rollapply(unlist(probGeno), 3, by = 3, sum, 
                              partial = FALSE)
    return(expGeno )
}

#' Rewrite expected genotypes into genotype probabilities.
#'
#' Convert genotype frequencies to genotypes encoded as triplets of probablities 
#' (p(AA), p(Aa), p(aa)).
#'
#' @param geno Vector [numeric] with genotypes 
#' @return Numeric vector of length [length(geno)*3] with the genotype encoded 
#' as probabbilities (p(AA), p(Aa), p(aa)).
#' @export
#' @examples
#' nrSamples <- 10
#' # Simulate binomial SNP with 0.2 allele frequency
#' geno <- rbinom(nrSamples, 2, p=0.2)
#' geno_prob<- expGen2probGen(geno)
expGen2probGen <- function(geno) {
    if (is.list(geno)) {
        stop("expGen2probGen takes a vector of genotypes as input,
             but list provided.")
    } 
    if (length(geno) <= 1) {
        if (!is.na(geno) && !is.numeric(geno) ) {
            stop("expGen2probGen takes a vector of genotypes as input,
             but ",  typeof(geno), " provided.")
        }
    } else if (!is.vector(geno, mode="numeric")) {
        stop("expGen2probGen takes a vector of genotypes as input,
             but ",  typeof(geno), " provided.")
    }
    unlist(lapply(geno, function(g) {
        if (is.na(g)) return( c(NA,NA,NA))
        else if (g == 0) return( c(1,0,0))
        else if (g == 1) return( c(0,1,0))
        else if (g == 2) return( c(0,0,1))
        else {
            stop("Genotypes can only be encoded as 0,1,2 or denoted as missing", 
                "via NA")
        }
    }))
}
#' Test lists for different properties of numerics.
#' 
#' Test all elements of a list if they are numeric, positive numbers, integers
#' or proportions (range 0-1)
#' 
#' @param numbers List whose elements are tested for being numeric
#' @param positives List whose elements are tested for being positive numbers
#' @param integers List whose elements are tested for being integers
#' @param proportions List whose elements are tested for being proportions 
#' between 0 and 1
#' @return NULL
testNumerics <- function(numbers, positives=NULL, integers=NULL, 
                         proportions=NULL) {
    notNullnotNumeric <- function(x) !is.null(x) && !is.numeric(x)
    notNullnotPositive <- function(x) !is.null(x)  && any(x <= 0)
    notNullnotInteger <- function(x) !is.null(x)  && any(x%%1!=0) 
    notNullNotInRange <- function(x) !is.null(x)  && (any(x < 0) || any(x > 1)) 
    if (any(sapply(numbers, notNullnotNumeric))) {
        notNumbers <- which(sapply(numbers, notNullnotNumeric))
        stop(paste(names(numbers)[notNumbers], collapse=","), 
             " is/are not numeric (",  
             paste(numbers[notNumbers], collapse=","), ")", sep="")
    }
    if (!is.null(positives)) {
        if (any(sapply(positives, notNullnotPositive))) {
            notPositives <- which(sapply(positives, notNullnotPositive))
            stop(paste(names(positives)[notPositives], collapse=","), 
                 " has/have to be greater than zero (",  
                 paste(positives[notPositives], collapse=","), ")", 
                 sep="")
        }  
    }
    if (!is.null(integers)) {
        if (any(sapply(integers, notNullnotInteger))) {
            notIntegers <- which(sapply(integers, notNullnotInteger))
            stop(paste(names(integers)[notIntegers], collapse=","), 
                 " has/have to be integers, given ",  
                 paste(integers[notIntegers], collapse=","), 
                 sep="")
        }  
    }
    if (!is.null(proportions)) {
        if (any(sapply(proportions, notNullNotInRange))) {
            outOfRange <- which(sapply(proportions, notNullNotInRange))
            stop("Proportions have to be specified between 0 and 1: ", 
                 paste(names(proportions)[outOfRange], collapse=","), 
                 " are outside of this range (",  
                 paste(proportions[outOfRange], collapse=","), ")", sep="")
        }  
    }
}
