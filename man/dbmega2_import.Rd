% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mega2rcreate.R
\name{dbmega2_import}
\alias{dbmega2_import}
\title{read Mega2 SQLite database into R}
\usage{
dbmega2_import(dbname,
               bpPosMap = NULL,
               verbose = FALSE)
}
\arguments{
\item{dbname}{file path to SQLite database.}

\item{bpPosMap}{index that specifies which map in the map_table should be used for
marker chromosome/position.  If it is NULL, the internal variable
\emph{base_pair_position_index} is used instead.
\code{showMapNames()} shows the association between map name and map number.}

\item{verbose}{print out statistics on the name/size of each table read and show column headers.
Also, save the verbose value for use by other Mega2R functions.}
}
\value{
envir an environment that contains all the data frames made from the SQLite database.
}
\description{
Read the fields of SQLite data base tables that are required for Mega2R into
 data frames.
 These data frames are stored in an 'environment' which is returned.
 This function also adds some state data, extra data frames, and computed data frames
 to the 'environment'.
}
\examples{
db = system.file("exdata", "seqsimm.db", package="Mega2R")
ENV = dbmega2_import(db, verbose = TRUE)

ENV = dbmega2_import(db)
}
