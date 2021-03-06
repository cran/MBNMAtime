% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rank.functions.R
\name{rankauc}
\alias{rankauc}
\title{Calculates ranking probabilities for AUC from a time-course MBNMA}
\usage{
rankauc(mbnma, decreasing = FALSE, treats = NULL, int.range = NULL,
  n.iter = mbnma$BUGSoutput$n.sims, ...)
}
\arguments{
\item{mbnma}{An S3 object of class \code{"mbnma"} generated by running
a time-course MBNMA model}

\item{decreasing}{A boolean object to indicate whether higher values are better (\code{decreasing=TRUE}) or
worse (\code{decreasing=FALSE}).}

\item{treats}{A character vector of treatment/class names (depending on the value of \code{level}). If left \code{NULL`` then rankings will be calculated for all treatments/classes. Note that unlike }rank.mbnma()` this argument
cannot take a numeric vector.}

\item{int.range}{A numeric vector with two elements that indicates the range
over which to calculate AUC. Takes the form c(lower bound, upper bound). If left
as \code{NULL} (the default) then the range will be between zero and the maximum follow-up
time in the data for the treatments specified in \code{treats}.}

\item{n.iter}{The number of iterations for which to calculate AUC (if \code{"auc"} is included in \code{params}).
Must be a positive integer. Default is the value used in \code{mbnma}.}

\item{...}{Arguments to be sent to \code{integrate()}}
}
\value{
A named list whose elements correspond to parameters given in
\code{params}. Each element contains:
\itemize{
\item \code{summary.rank} A data frame containing
mean, sd, and quantiles for the ranks of each treatment given in \code{treats}
\item \code{prob.matrix} A matrix of the proportions of MCMC results for which each
treatment in \code{treats} ranked in which position for the given parameter
\item \code{rank.matrix} A matrix of the ranks of MCMC results for each treatment in
\code{treats} for the given parameter.
}
}
\description{
Calculates ranking probabilities for AUC from a time-course MBNMA
}
\details{
\code{"auc"} can be included in \code{params} to rank treatments based on
Area Under the Curve (AUC). This accounts for the effect of multiple
time-course parameters simultaneously on the treatment response, but will
be impacted by the range of time over which AUC is calculated (\code{int.range}).
Currently \code{"auc"} cannot be ranked for class effect models.

As with other post-estimation functions, \code{rank()} should only be performed on
models which have successfully converged. Note that rankings can be very sensitive to
even small changes in treatment effects and therefore failure to converge in only
one parameter may have substantial impact on rankings.
}
