#' Print userinfo.
#'
#' Wrapper function around \code{\link{message}} that allows to turn the 
#' printing of messages to standard out.
#' on or off
#'
#' @param userinfo string or vector of string elements and variables
#' @param verbose [boolean], if TRUE message is displayed on standard out, if 
#' FALSE, message is suppressed
#' @param sep delimiter [string] to separate message elements when userinfo 
#' given as vector
#' @seealso \code{\link{message}} which this function wraps
vmessage <- function(userinfo, verbose=TRUE, sep=" ") {
    if (verbose) {
        message(paste(userinfo, collapse=sep))
    }
}

#' Comma-separated string to numeric vector.
#'
#' Split input of comma-separated string into vector of numeric values.
#'
#' @param commastring input character vector containing numbers separated by 
#' commas
#' @return numeric vector of values extracted from commastring
commaList2vector <- function(commastring) {
    as.numeric(unlist(strsplit(commastring, ",")))
}

#' Add all non-NULL elements of list.
#'
#' @param compList list of numeric matrices or data.frames of the equal 
#' dimensions
#' @return matrix or data.frame containing sum of all list elements where 
#' \code{is.null} is FALSE
addNonNulls <- function(compList) {
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
#' @param x the number [integer] of observations to simulate
#' @param dist name of distribution [string] from which the observations are 
#' drawn. 'norm' is the normal distribution, 'unif' the uniform distribution 
#' 'bin' the binomial distribution, "cat_norm" samples categorical variables 
#' according to a normal distribution and "cat_unif" 
#' according to a uniform distribution. For "cat_norm", length(category)/2 is 
#' used mean for the normal distribution unless
#' specified otherwise.
#' @param m the mean of the normal distribution [double]/the mean between min 
#' and max for the uniform distribution [double]/ 
#' the rank of the category to be used as mean for "cat_norm" [integer]
#' @param std the standard deviation of the normal distribution or the distance 
#' of min/max from the mean for the uniform distribution [double]
#' @param categories number of categories [integer] for simulating categorical 
#' variables (for distr="cat_norm" or "cat_unif")
#' @param prob the probability [double] of success for each trial 
#' (for distr="bin")
#' @return numeric vector of length [x] with the sampled values
#' @seealso \code{\link{runif}}, \code{\link{rnorm}}, \code{\link{rbinom}} for 
#' documentation of the underlying distributions
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
    if (dist == "unif") {
        if (is.null(m)) m <- 0
        d <- runif(n=x, min=m - std, max=m + std)
    }
    else if (dist == "norm") {
        if (is.null(m)) m <- 0
        d <- rnorm(n=x, mean=m, sd=std)
    }
    else if (dist == "bin") {
        if (prob < 0) {
            stop(paste("Simulating binomial distribution: Probability has to",
                "be between 0 and 1"))
        }
        d <- rbinom(n=x, size=1, prob=prob)
    }
    else if (grepl("cat", dist)) {
        if (dist == "cat_norm") {
            # generate probabilities for categories
            if (is.null(m)) m <- median(1:categories)
            ptmp=sapply(1:categories, pnorm, mean=m, sd=std)
            prob = c(ptmp[1],diff(ptmp))/sum(c(ptmp[1],diff(ptmp)))
        }
        else if (dist == "cat_unif") {
            prob=rep(1/categories, categories)
        }
        d = sample(1:categories, x, replace=TRUE, prob=prob)
        d <- as.numeric(d)
    } else {
        stop("Unknown distribution")
    }
    return(d)
}
