#' @useDynLib PhenotypeSimulator
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
    library.dynam.unload(getPackageName(), libpath)
}
