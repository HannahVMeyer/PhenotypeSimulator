% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genotypeFunctions.R
\name{readStandardGenotypes}
\alias{readStandardGenotypes}
\title{Read genotypes from file.}
\usage{
readStandardGenotypes(
  N,
  filename,
  format = NULL,
  verbose = TRUE,
  sampleID = "ID_",
  snpID = "SNP_",
  delimiter = ","
)
}
\arguments{
\item{N}{Number [integer] of samples to simulate.}

\item{filename}{path/to/genotypefile [string] in plink, oxgen
(impute2/snptest/hapgen2), genome, bimbam or [delimiter]-delimited format (
for format information see \emph{External genotype software and formats}).}

\item{format}{Name [string] of genotype file format. Options are: "plink",
"oxgen", "genome", "bimbam" or "delim".}

\item{verbose}{[boolean] If TRUE, progress info is printed to standard out.}

\item{sampleID}{Prefix [string] for naming samples (will be followed by
sample number from 1 to N when constructing id_samples).}

\item{snpID}{Prefix [string] for naming SNPs (will be followed by SNP number
from 1 to NrSNP when constructing id_snps).}

\item{delimiter}{Field separator [string] of genotype file when
format == "delim".}
}
\value{
Named list of [NrSamples X NrSNPs] genotypes (genotypes), their
[NrSNPs] SNP IDs (id_snps), their [NrSamples] sample IDs (id_samples) and
format-specific additional files (such as format-specific genotypes encoding
or sample information; might be used for output writing).
}
\description{
readStandardGenotypes can read genotypes from a number of input
formats for standard GWAS (binary plink, snptest, bimbam, gemma) or
simulation software (binary plink, hapgen2, genome). Alternatively, simple
text files (with specified delimiter) can be read. For more information on
the different file formats see \emph{External genotype software and formats}.
}
\details{
The file formats and software formates supported are described
below. For large external genotypes, consider the option to only read
randomly selected, causal SNPs into memory via \link{getCausalSNPs}.
}
\section{External genotype software and formats}{

\subsection{Formats:}{
\itemize{
\item PLINK: consists of three files: .bed, .bim and .fam. When specifying
the filepath, only the core of the name without the ending should be
specified (i.e. for geno.bed, geno.bim and geno.fam, geno should be
specified). When reading data from plink files, the absolute path to the
plink-format file has to be provided, tilde expansion not provided.
From \url{https://www.cog-genomics.org/plink/1.9/formats}:  The .bed
files contain the  primary representation of genotype calls at biallelic
variants in a binary format. The .bim is a text file with no header
line, and one line  per variant with the following six fields: i) Chromosome
code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or
name, ii) Variant identifier, iii) Position in morgans or centimorgans (safe
to use dummy value of '0'), iv) Base-pair coordinate (normally 1-based, but
0 ok; limited to 231-2), v) Allele 1 (corresponding to clear bits in .bed;
usually minor), vi) Allele 2 (corresponding to set bits in .bed; usually
major). The .fam file is a text file with no header line, and one line per
sample with the following six fields: i) Family ID ('FID'), ii), Within-
family ID ('IID'; cannot be '0'), iii) Within-family ID of father ('0' if
father isn't in dataset, iv) within-family ID of mother ('0' if mother isn't
in dataset), v) sex code ('1' = male, '2' = female, '0' = unknown), vi)
Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing
data if case/control)
\item oxgen: consists of two files: the space-separated genotype file ending
in .gen and the space-separated sample file ending in .sample. When
specifying the filepath, only the core of the name without the ending should
be specified (i.e. for geno.gen and geno.sample, geno should be specified).
From
 \url{https://www.well.ox.ac.uk/~gav/snptest/#input_file_formats}:
The genotype file stores data on a one-line-per-SNP format. The first five
entries of each line should be the SNP ID, RS ID of the SNP, base-pair
position of the SNP, the allele coded A and the allele coded B. The SNP ID
can be used to denote the chromosome number of each SNP. The next three
numbers on the line should be the probabilities of the three genotypes AA,
AB and BB at the SNP for the first individual in the cohort. The next three
numbers should be the genotype probabilities for the second individual in the
cohort. The next three numbers are for the third individual and so on. The
order of individuals in the genotype file should match the order of the
individuals in the sample file. The sample file has three parts (a) a header
line detailing the names of the columns in the file, (b) a line detailing the
types of variables stored in each column, and (c) a line for each individual
detailing the information for that individual. For more information on the
sample file visit the above url or see \link{writeStandardOutput}.
\item genome: The entire output of genome can be saved via `genome -options >
outputfile`. The /path/to/outputfile should be provided and this function
extracts the relevant genotype information from this output file.
\url{http://csg.sph.umich.edu/liang/genome/}
\item bimbam: Mean genotype file format of bimbam which is a single, comma-
separated file, without information on individuals. From the documentation
for bimbam at \url{http://www.haplotype.org/software.html}: the first column
of the mean genotype files is the SNP ID, the second and third columns are
allele types with minor allele first. The remaining columns are the mean
genotypes of different individuals – numbers between 0 and 2 that represents
the (posterior) mean genotype, or dosage of the minor allele. \item delim: a
[delimter]-delimited file of [(NrSNPs+1) x (NrSamples+1)] genotypes with the
snpIDs in the first column and the sampleIDs in the first row and genotypes
encoded as numbers between 0 and 2 representing the (posterior) mean
genotype, or dosage of the minor allele. Can be user-genotypes or genotypes
simulated with foward-time algorithms such as simupop
(\url{http://simupop.sourceforge.net/Main/HomePage})
that allow for user-specified output formats.
}}

\subsection{Genotype simulation characteristics:}{
\itemize{
\item PLINK: simple, bi-allelelic genotypes without LD structure, details can
be found at \url{https://www.cog-genomics.org/plink/1.9/input#simulate}.
\item Hapgen2: resampling-based genotype simulation, details can be found at
\url{http://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html}.
\item Genome: coalescent-based genotype simulation, details can be found at
\url{http://csg.sph.umich.edu/liang/genome/GENOME-manual.pdf}.
}
Examples on how to call these genotype simulation tools can be found in the
vignette \emph{sample-scripts-external-genotype-simulation}.}
}

\examples{
# Genome format
filename_genome  <- system.file("extdata/genotypes/genome/",
"genotypes_genome.txt",
package = "PhenotypeSimulator")
data_genome <- readStandardGenotypes(N=100, filename_genome, format ="genome")

filename_hapgen  <- system.file("extdata/genotypes/hapgen/",
"genotypes_hapgen.controls.gen",
package = "PhenotypeSimulator")
filename_hapgen <- gsub("\\\\.gen", "", filename_hapgen)
data_hapgen <- readStandardGenotypes(N=100, filename_hapgen, format='oxgen')

filename_plink  <- system.file("extdata/genotypes/plink/",
"genotypes_plink.bed",
package = "PhenotypeSimulator")
filename_plink <- gsub("\\\\.bed", "", filename_plink)
data_plink <- readStandardGenotypes(N=100, filename=filename_plink,
format="plink")

filename_delim  <- system.file("extdata/genotypes/",
"genotypes_chr22.csv",
package = "PhenotypeSimulator")
data_delim <- readStandardGenotypes(N=50, filename=filename_delim,
format="delim")
}
