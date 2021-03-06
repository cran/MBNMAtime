% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write.functions.R
\name{write.cor}
\alias{write.cor}
\title{Adds correlation between time-course relative effects}
\usage{
write.cor(model, var.scale = NULL, class.effect = list())
}
\arguments{
\item{model}{A character object of JAGS MBNMA model code}

\item{var.scale}{A numeric vector indicating the relative scale of variances between
correlated time-course parameters when relative effects are modelled on more than
one time-course parameter(see Details LINK). Each element of
the vector refers to the relative scale of each of the time-course parameters that is
modelled using relative effects.}

\item{class.effect}{A list of named strings that determines which time-course
parameters to model with a class effect and what that effect should be
(\code{"common"} or \code{"random"}). For example: \code{list("beta.2"="common", "beta.3"="random")}.}
}
\description{
This uses a Wishart prior as default for modelling the correlation
}
