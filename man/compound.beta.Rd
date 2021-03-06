% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.functions.R
\name{compound.beta}
\alias{compound.beta}
\title{Prepares beta time-course parameters for mb.write()}
\usage{
compound.beta(beta.1)
}
\arguments{
\item{beta.1}{A two-element list whose elements have the following names:
\itemize{
\item \code{pool} Can take either \code{"rel"}, \code{"arm"} or \code{"const"}
\item \code{method} Can take either \code{"common"}, \code{"random"}, or be assigned a numeric value
}}
}
\description{
Checks that beta time-course parameters have been specified correctly and converts them
to the correct format for \code{mb.write()} and other subsequent functions.
}
