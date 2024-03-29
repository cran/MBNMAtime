% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mbnma-class.R
\name{rank.mbnma}
\alias{rank.mbnma}
\title{Rank parameters from a time-course MBNMA}
\usage{
\method{rank}{mbnma}(
  x,
  param = "auc",
  lower_better = FALSE,
  treats = NULL,
  int.range = NULL,
  n.iter = x$BUGSoutput$n.sims,
  ...
)
}
\arguments{
\item{x}{An S3 object of \code{class("mbnma")} generated by running
a time-course MBNMA model}

\item{param}{A character object containing any model parameter monitored
in \code{mbnma} for which ranking is desired (e.g. \code{"beta.1"}, \code{"emax"}).
Parameters must vary by treatment for ranking to be possible. Can also be specified as
\code{"auc"} (see details).}

\item{lower_better}{Indicates whether negative responses are better (\code{lower_better=TRUE}) or
positive responses are better (\code{lower_better=FALSE})}

\item{treats}{A character vector of treatment/class names (depending on the parameter to be ranked) or
a numeric vector of treatment/class codes (as coded in \code{mbnma})
that indicate which treatments/classes to calculate rankings for. If left `NULL``
then rankings will be calculated for all treatments/classes.}

\item{int.range}{A numeric vector with two elements that indicates the range
over which to calculate AUC. Takes the form c(lower bound, upper bound). If left
as \code{NULL} (the default) then the range will be between zero and the maximum follow-up
time in the dataset.}

\item{n.iter}{The number of iterations for which to calculate AUC (if \code{"auc"} is included in \code{params}).
Must be a positive integer. Default is the value used in \code{mbnma}.}

\item{...}{Arguments to be sent to \code{integrate()}}
}
\value{
A named list whose elements include:
\itemize{
\item \code{summary.rank} A data frame containing
mean, sd, and quantiles for the ranks of each treatment given in \code{treats}
\item \code{prob.matrix} A matrix of the proportions of MCMC results for which each
treatment/class in \code{treats} ranked in which position for the given parameter
\item \code{rank.matrix} A matrix of the ranks of MCMC results for each treatment/class in
\code{treats} for the given parameter.
}
}
\description{
Ranks desired parameters saved from a time-course MBNMA model from "best" to
"worst".
}
\details{
\code{"auc"} can be specified in \code{param} to rank treatments based on
Area Under the Curve (AUC). This accounts for the effect of multiple
time-course parameters simultaneously on the treatment response, but will
be impacted by the range of time over which AUC is calculated (\code{int.range}).
This requires integration over \code{int.range} and can take some time to run (particularly)
for spline functions as this uses the trapezoid method rather than adaptive quadrature).
Note that \code{"auc"} can only be calculated at the treatment-level in class effect models.

As with other post-estimation functions, \code{rank()} should only be performed on
models which have successfully converged. Note that rankings can be very sensitive to
even small changes in treatment effects and therefore failure to converge in only
one parameter may have substantial impact on rankings.
}
\examples{
\donttest{
# Create an mb.network object from a dataset
network <- mb.network(alog_pcfb)

# Run an MBNMA model with an Emax time-course
emax <- mb.run(network,
  fun=temax(pool.emax="rel", method.emax="common",
            pool.et50="rel", method.et50="random"),
  intercept=FALSE)

# Rank treatments by time-course parameter from the model with lower scores being better
rank(emax, param=c("emax"), lower_better=TRUE)

# Rank treatments 1-3 by AUC
rank(emax, param="auc", treats=c(1:3), lower_better=TRUE,
  int.range=c(0,20))
}

}
