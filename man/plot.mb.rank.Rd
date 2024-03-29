% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mb.rank-class.R
\name{plot.mb.rank}
\alias{plot.mb.rank}
\title{Plot histograms of rankings from MBNMA models}
\usage{
\method{plot}{mb.rank}(x, treat.labs = NULL, ...)
}
\arguments{
\item{x}{An object of class \code{"mb.rank"} generated by \code{rank.mbnma()}}

\item{treat.labs}{A vector of treatment labels in the same order as treatment codes.
Easiest to use treatment labels stored by \code{mb.network()}}

\item{...}{Arguments to be sent to \code{ggplot2::ggplot()}}
}
\value{
A histogram that shows rankings for each treatment/agent/prediction.
The object returned is an object of class \code{c("gg", "ggplot")}.
}
\description{
Plot histograms of rankings from MBNMA models
}
\examples{
\donttest{
# Create an mb.network object from a dataset
painnet <- mb.network(osteopain)

# Run an MBNMA model with an Emax time-course
emax <- mb.run(painnet,
  fun=temax(pool.emax="rel", method.emax="common",
    pool.et50="abs", method.et50="random"),
  positive.scale=TRUE)

# Calculate treatment rankings for AUC and emax
ranks <- rank(emax,
  param=c("auc"),
  int.range=c(0,15), n.iter=500)

# Plot histograms for ranking by AUC
plot(ranks)
}

}
