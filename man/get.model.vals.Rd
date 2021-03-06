% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.functions.R
\name{get.model.vals}
\alias{get.model.vals}
\title{Get MBNMA model values}
\usage{
get.model.vals(mbnma, timecourse, beta.incl, E0 = 0)
}
\arguments{
\item{mbnma}{An S3 object of class \code{"mbnma"} generated by running
a time-course MBNMA model}

\item{timecourse}{A character object that specifies the time-course used in
\code{mbnma} (in terms of alpha, beta, and time), as generated by
\code{init.predict()}}

\item{beta.incl}{A numeric vector that indicates the time-course parameters
that were included in \code{mbnma}, as generated by \code{init.predict()}}

\item{E0}{An object to indicate the value(s) to use for the response at time = 0
in the prediction. This can take a number of different formats depending
on how it will be used/calculated. The default is \code{0} but this may lead
to non-sensical predictions.
\itemize{
\item \code{numeric()} A single numeric value representing the deterministic response at time = 0,
given.
\item \code{character()} A single string representing a stochastic distribution for the response
at time = 0. This is specified as a random number generator
(RNG) given as a string, and can take any RNG distribution for which a function exists
in R. For example: \code{"rnorm(n, 7, 0.5)"}.
}}
}
\value{
A list containing named elements that correspond to different
time-course parameters in \code{mbnma}. These elements contain MCMC results
either taken directly from \code{mbnma} or (in the case of random time-course
parameters specified as \code{method="random"}) randomly
generated using parameter values estimated in \code{mbnma}.

Additional elements contain the following values:
\itemize{
\item \code{timecourse} A character object that specifies the time-course used in \code{mbnma} in terms of
alpha, beta, mu, d and time. Consistency relative time-course parameters
are specified in terms of mu and d.
\item \code{mu.prior} A character vector that
indicates for which time-course parameters a network reference treatment
effect will be required for prediction.
\item \code{time.params} A character vector
that indicates the different time-course parameters that are required for
the prediction
}

@noRd
}
\description{
Extracts specific information required for prediction from a time-course
MBNMA model
}
