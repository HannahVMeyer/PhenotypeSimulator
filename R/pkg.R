#' @useDynLib PhenotypeSimulator
NULL

.onUnload <- function (libpath) {
    library.dynam.unload(getPackageName(), libpath)
}
