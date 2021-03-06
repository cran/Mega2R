% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mega2rtutorial.R
\name{dump_mega2rtutorial_data}
\alias{dump_mega2rtutorial_data}
\title{dump tutorial data}
\usage{
dump_mega2rtutorial_data(dir = file.path(tempdir(), "Mega2Rtutorial"))
}
\arguments{
\item{dir}{The directory to store the tutorial data to.  By default, this is
tempdir()/Mega2Rtutorial}
}
\value{
None
}
\description{
This function retrieves data stored in the Mega2rtutorial (inst/exdata).  It
dumps them in the specified directory.
}
\examples{
dump_mega2rtutorial_data()

}
