% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mega2vcf.R
\name{mkVCFpen}
\alias{mkVCFpen}
\title{generate required Mega2 penetrance (.pen) file}
\usage{
mkVCFpen(prefix, envir)
}
\arguments{
\item{prefix}{prefix for .pen file (penetrance)}

\item{envir}{"environment" containing SQLite database and other globals}
}
\value{
None
}
\description{
Generate the initial boiler plate VCF, then generate ##INFO entries for each entry tag.
 Finally, generate the ##contig entries for each chromosome.
}
\examples{
\dontrun{
mkVCFpen(prefix, envir)
}
}
\keyword{internal}
