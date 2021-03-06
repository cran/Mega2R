% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mega2riterate.R
\name{applyFnToMarkers}
\alias{applyFnToMarkers}
\title{apply a function to the genotypes from a set of markers}
\usage{
applyFnToMarkers(op      = function (markers, range, envir) {},
                markers_arg,
                envir = ENV)
}
\arguments{
\item{op}{Is a function of three arguments.  It will be called once by
\code{applyFnToMarkers} in a try/catch context.  The arguments are:
\describe{
\item{markers}{Marker data for each marker in \strong{geno}.  A marker is a data frame with the following 5 observations:
\describe{
\item{locus_link}{is the ordinal ranking of this marker among all loci}
\item{locus_link_fill}{is the position of corresponding marker genotype data in the
\emph{unified_genotype_table}}
\item{MarkerName}{is the text name of the marker}
\item{chromosome}{is the integer chromosome number}
\item{position}{is the integer base pair position of marker}
 }
}
\item{range}{NULL: to indicate no explicit range was specified.}
\item{envir}{An 'environment' holding Mega2R data frames and state data.}
}}

\item{markers_arg}{a data frame with the following 5 observations:
\describe{
\item{locus_link}{is the ordinal ranking of this marker among all loci}
\item{locus_link_fill}{is the position of corresponding marker genotype data in the
\emph{unified_genotype_table}}
\item{MarkerName}{is the text name of the marker}
\item{chromosome}{is the integer chromosome number}
\item{position}{is the integer base pair position of marker}
 }}

\item{envir}{an 'environment' that contains all the data frames created from the SQLite database.}
}
\value{
None
}
\description{
A matrix of the genotypes for all the specified markers is generated.  Then, the call back function, \code{op},
  is called with the markers, NULL (for the range), and the 'environment'.
}
\examples{
  db = system.file("exdata", "seqsimm.db", package="Mega2R")
  ENV = read.Mega2DB(db)
  show = function(m, r, e) {
      print(r)
      print(m)
      print(head(getgenotypes(m, envir = e)))
  }

   # apply function "show" to all genotypes > 5,000,000 bp
   applyFnToMarkers(show, ENV$markers[ENV$markers$position > 5000000,])


}
