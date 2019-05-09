#' PhenotypeSimulator: A package for simulating phenotypes from different
#' genetic and noise components
#'
#' @docType package
#' @name PhenotypeSimulator
#' @import snpStats
#' @import ggplot2 
#' @importFrom cowplot plot_grid
#' @importFrom methods getPackageName as
#' @importFrom stats median pnorm rbinom rnorm runif sd var na.omit rexp
#' @importFrom utils read.table write.table 
#' @importFrom optparse make_option parse_args OptionParser
#' @importFrom dplyr filter
#' @importFrom reshape2 melt
#' @importFrom snpStats write.plink
#' @importClassesFrom snpStats SnpMatrix
#' @importFrom R.utils countLines
#' @keywords internal
NULL
