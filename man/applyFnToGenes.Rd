% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mega2riterate.R
\name{applyFnToGenes}
\alias{applyFnToGenes}
\title{apply a function to the genotypes (markers) in each gene transcript and/or base pair range}
\usage{
applyFnToGenes(op           = function (markers, range, envir) {},
               genes_arg    = NULL,
               ranges_arg   = matrix(ncol = 3, nrow = 0),
               chrs_arg     = vector("integer", 0),
               markers_arg  = vector("character", 0),
               type_arg     = "TX",
               fuzz_arg     = 0,
               envir        = ENV)
}
\arguments{
\item{op}{Is a function of three arguments.  It will be called repeatedly by
\code{applyFnToGenes} in a try/catch context.  The arguments are:
\describe{
\item{markers}{Marker data for each marker selected.  A marker is a data frame with the following 5 observations:
\describe{
\item{locus_link}{is the ordinal ranking of this marker among all loci}
\item{locus_link_fill}{is the position of corresponding marker genotype data in the
\emph{unified_genotype_table}}
\item{MarkerName}{is the text name of the marker}
\item{chromosome}{is the integer chromosome number}
\item{position}{is the integer base pair position of marker}
 }
}
\item{range}{An indicator of which range argument these markers correspond to.}
\item{envir}{An 'environment' holding Mega2R data frames and state data.}
}}

\item{genes_arg}{a character vector of gene names.
All the transcripts identified with the specified gene in BioConductor Annotation,\cr
\bold{TxDb.Hsapiens.UCSC.hg19.knownGene}, are selected.  This produces multiple "range"
elements containing chromosome, start base pair, end base pair.  (If the gene name is "*",
all the transcript will be selected.) Note: BioCoductor Annotation
\bold{org.Hs.eg.db} is used to convert from gene name to ENTREZ gene id.}

\item{ranges_arg}{an integer matrix of three columns.  The columns define a range:
a chromosome number, a start base pair value, and an end base pair value.}

\item{chrs_arg}{an integer vector of chromosome numbers.  All of the base pairs on each
chromosomes will be selected as a single range.}

\item{markers_arg}{a data frame with the following 5 observations:
\describe{
\item{locus_link}{is the ordinal ranking of this marker among all loci}
\item{locus_link_fill}{is the position of corresponding marker genotype data in the\cr
\emph{unified_genotype_table}}
\item{MarkerName}{is the text name of the marker}
\item{chromosome}{is the integer chromosome number}
\item{position}{is the integer base pair position of marker}
}}

\item{type_arg}{a character vector of length 1 that contains \bold{"TX"} or does not.  If it is
\bold{"TX"}, which is the default, the \bold{TX} fields of BioConductor Annotation,\cr
\bold{TxDb.Hsapiens.UCSC.hg19.knownGene} are used to define the base pair ranges and chromosome.
Otherwise, the \bold{CDS} fields are used.}

\item{fuzz_arg}{is an integer vector of length one or two.  The first argument is used to reduce
the start base pair selected from each transcript and the second to increase the end base pair
position.  (If only one value is present, it is used for both adjustments.)  Note: The values
can be positive or negative.}

\item{envir}{an 'environment' that contains all the data frames created from the SQLite database.}
}
\value{
None
}
\description{
This function generates base pair ranges from its input arguments.
 Each range specifies a chromosome, a start
 base pair and end base pair.  Typically, a range could be a gene transcript, though
 it could be a whole chromosome, or a run of base pairs on a chromosome.  Once the
 ranges are generated, \code{applyFnToRanges} is called to find all the
 rows (i.e. markers) from the \emph{markers} data frame that fall in each range.  For these
 markers, a matrix of the genotypes is generated.  Finally, the \code{op} function is called for
 each range with the arguments: markers, range, and 'environment'.
}
\note{
If you want subsequent calls to \code{op} to share information, data can be placed in
 a data frame that is added to the 'environment'.
}
\examples{
  db = system.file("exdata", "seqsimm.db", package="Mega2R")
  ENV = read.Mega2DB(db)

  show = function(m, r, e) {
      print(r)
      print(m)
      print(head(getgenotypes(m, envir = e)))
  }

   # apply function "show" to all transcripts on genes ELL2 and CARD15
\donttest{
    # donttestcheck: time
    applyFnToGenes(show, genes_arg = c("CEP104"))
}

   # apply function "show" to all genotypes on chromosomes 11 for two base
   # pair ranges
   applyFnToGenes(show, ranges_arg = matrix(c(1, 5000000, 10000000,
                  1, 10000000, 15000000), ncol = 3, nrow = 2, byrow = TRUE))

   # apply function "show" to all genotypes for first marker in each chromosome
   applyFnToGenes(show, markers_arg = ENV$markers[! duplicated(ENV$markers$chromosome), 3])

   # apply function "show" to all genotypes on chromosomes 24 and 26
   applyFnToGenes(show, chrs_arg=c(24, 26))


}
