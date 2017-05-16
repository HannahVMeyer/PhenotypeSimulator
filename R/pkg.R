#' @useDynLib PhenotypeSimulator, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
    library.dynam.unload(getPackageName(), libpath)
}
