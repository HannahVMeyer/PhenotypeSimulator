% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genotypeFunctions.R
\name{getKinship}
\alias{getKinship}
\title{Get genetic kinship.}
\usage{
getKinship(
  N,
  sampleID = "ID_",
  X = NULL,
  kinshipfile = NULL,
  id_samples = NULL,
  standardise = FALSE,
  sep = ",",
  header = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{N}{Number [integer] of samples to simulate.}

\item{sampleID}{Prefix [string] for naming samples (will be followed by
sample number from 1 to N when constructing id_samples).}

\item{X}{[NrSamples x totalNrSNPs] Matrix of (standardised) genotypes.}

\item{kinshipfile}{path/to/kinshipfile [string] to be read; either X or
kinshipfile must be provided.}

\item{id_samples}{Vector of [NrSamples] sample IDs [string]; if not provided
constructed by paste(sampleID, 1:N, sep="").}

\item{standardise}{[boolean] If TRUE genotypes will be standardised before
kinship estimation.}

\item{sep}{Field separator [string] of kinship file.}

\item{header}{[boolean], If TRUE kinship file has header information.}

\item{verbose}{[boolean]; If TRUE, progress info is printed to standard out}
}
\value{
[NrSamples x NrSamples] Matrix of kinship estimate.
}
\description{
Estimate kinship from standardised genotypes or read pre-computed kinship
file. Standardised genotypes can be obtained via
\code{\link{standardiseGenotypes}}.
}
\details{
The kinship is estimated as \eqn{K = XX_T}, with X the standardised
genotypes of the samples. When estimating the kinship from the provided
genotypes, the kinship is normalised by the mean of its diagonal
elements and 1e-4 added to the diagonal for numerical stability.
}
\examples{
geno <- simulateGenotypes(N=10, NrSNP=50)
K_fromGenotypesNormalised <- getKinship(N=10, X=geno$genotypes,
standardise=TRUE)

kinshipfile <- system.file("extdata/kinship",
"kinship.csv",
package = "PhenotypeSimulator")
K_fromFile <- getKinship(N=50, kinshipfile=kinshipfile)
}
