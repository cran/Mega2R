% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mega2gdsfmt.R
\name{mkAlleles}
\alias{mkAlleles}
\title{generate allele pairs in with MAJ(or) allele first}
\usage{
mkAlleles(markers = NULL, separator = "/", envir = ENV)
}
\arguments{
\item{markers}{data frame of markers to be processed}

\item{envir}{'environment' containing SQLite database and other globals}
}
\value{
None
}
\description{
The genotypes are ordered so the allele with the greater frequency appears first.
}
\examples{
\dontrun{
mkAlleles(envir)
}
}
\keyword{internal}
