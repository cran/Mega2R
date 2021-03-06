% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mega2rcreate.R
\name{getgenotypesraw}
\alias{getgenotypesraw}
\title{fetch genotype integer matrix for specified markers}
\usage{
getgenotypesraw(markers_arg, envir = ENV)
}
\arguments{
\item{markers_arg}{a data.frame with the following 5 observations:
\describe{
\item{locus_link}{is the ordinal ranking of this marker among all loci}
\item{locus_link_fill}{is the position of corresponding genotype data in the
\emph{unified_genotype_table}}
\item{MarkerName}{is the text name of the marker}
\item{chromosome}{is the integer chromosome number}
\item{position}{is the integer base pair position of marker}
 }}

\item{envir}{an environment that contains all the data frames created from the SQLite database.}
}
\value{
a matrix of genotypes represented as integers.  Each 32 bit integer represents contains
 two allele values: the high 16 bits contains the index of allele1 and the low 16 bits contains
 the index of allele2.  In the matrix, there is one column for each
 marker in the \emph{markers_arg} argument.  There is one row for each person in the family
 (\emph{fam}) table.
}
\description{
This function calls a C++ function that does all the heavy lifting.  It passes the arguments
 necessary for the C++ function: some from the caller's arguments and some from data frames
 that are in the "global" environment, \bold{envir}.  From its markers_arg argument, it gets
 the locus_index and the index in the \emph{unified_genotype_table}.
 From the "global" environment, \bold{envir}, it gets a bit vector of compressed genotype information,
 and some bookkeeping related data.
 Note: This function also contains a dispatch/switch on the type of compression in the genotype
 vector.  A different C++ function is called when there is compression versus when there is no
 compression.
}
\details{
The \emph{unified_genotype_table} contains one raw vector for each person.  In the vector,
 there are two bits for each genotype.  This function creates an output matrix by fixing
 the marker and collecting genotype information for each person and then repeating for
 all the needed markers.
}
\examples{
db = system.file("exdata", "seqsimm.db", package="Mega2R")
ENV = read.Mega2DB(db)

# two ints in upper/lower half integer representing allele # for all persons in chromosome 1
getgenotypesraw(ENV$markers[ENV$markers$chromosome == 1,])

}
