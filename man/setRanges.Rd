% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mega2riterate.R
\name{setRanges}
\alias{setRanges}
\title{set default range data: chromosome and start/end base pair}
\usage{
setRanges(ranges, indices, envir = ENV)
}
\arguments{
\item{ranges}{a data frame that contains at least 4 observations: a name, a chromosome, a start
base pair position and an end base pair position.}

\item{indices}{a vector of 3 or 4 integers that specify the chromosome column, start base pair,
column and end base pair column of range data frame and lastly the name column.  If the vector
only contains 3 integers, a name will be generated from the three range elements and it will be
appended to the ranges and the last range column will be added to the indices.}

\item{envir}{an 'environment' that contains all the data frames created from the SQLite database.}
}
\value{
None
}
\description{
This function sets the default list of ranges used by \code{applyFnToRanges}.  \code{applyFnToRanges}
 examines each range and the set of markers that fall within the range will be
 processed.
}
\examples{
db = system.file("exdata", "seqsimm.db", package="Mega2R")
ENV = read.Mega2DB(db)

ranges = matrix(c(1, 2240000, 2245000,
                  1, 2245000, 2250000,
                  1, 3760000, 3761000,
                  1, 3761000, 3762000,
                  1, 3762000, 3763000,
                  1, 3763000, 3764000,
                  1, 3764000, 3765000,
                  1, 3765000, 3763760,
                  1, 3763760, 3767000,
                  1, 3767000, 3768000,
                  1, 3768000, 3769000,
                  1, 3769000, 3770000),
                 ncol = 3, nrow = 12, byrow = TRUE)

setRanges(ranges, 1:3)

ENV$refRanges

ranges = matrix(c(1, 2240000, 2245000,
                  1, 2245000, 2250000,
                  1, 3760000, 3761000,
                  1, 3761000, 3762000,
                  1, 3762000, 3763000,
                  1, 3763000, 3764000,
                  1, 3764000, 3765000,
                  1, 3765000, 3763760,
                  1, 3763760, 3767000,
                  1, 3767000, 3768000,
                  1, 3768000, 3769000,
                  1, 3769000, 3770000),
                 ncol = 3, nrow = 12, byrow = TRUE)
ranges = data.frame(ranges)
ranges$name = LETTERS[1:12]
names(ranges) = c("chr", "start", "end", "name")

setRanges(ranges, 1:4)

ENV$refRanges
}
