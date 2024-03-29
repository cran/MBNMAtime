% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mb.predict-class.R
\name{plot.mb.predict}
\alias{plot.mb.predict}
\title{Plots predicted responses from a time-course MBNMA model}
\usage{
\method{plot}{mb.predict}(
  x,
  disp.obs = FALSE,
  overlay.ref = TRUE,
  overlay.nma = NULL,
  method = "random",
  col = "blue",
  max.col.scale = NULL,
  treat.labs = NULL,
  plot.bins = TRUE,
  ...
)
}
\arguments{
\item{x}{An object of class \code{"mb.predict"} generated by
\code{predict("mbnma")}}

\item{disp.obs}{A boolean object to indicate whether to show shaded sections
of the plot for where there is observed data (\code{TRUE}) or not (\code{FALSE})}

\item{overlay.ref}{A boolean object indicating whether to overlay a line
showing the median network reference treatment response over time on the
plot (\code{TRUE}) or not (\code{FALSE}). The network reference treatment (treatment
\enumerate{
\item must be included in \code{predict} for this to display the network reference
treatment properly.
}}

\item{overlay.nma}{Numeric vector used to overlay the results from a standard NMA model that
"lumps" time-points together within the time bin ranges specified in \code{overlay.nma}.
The numbers in \code{overlay.nma} define the boundaries of the time bins within which to perform
a standard NMA. Length must be >=2, or can be left as \code{NULL} (the default) to indicate that no NMA
should be perfomed. \code{overlay.nma} can only be specified if \code{overlay.ref==TRUE}.
See Details for further information.}

\item{method}{Can take \code{"common"} or \code{"random"} to indicate the type of NMA model used to synthesise data points
given in \code{overlay.nma}. The default is \code{"random"} since this assumes different
time-points in \code{overlay.nma} have been lumped together to estimate the NMA.}

\item{col}{A character indicating the colour to use for shading if \code{disp.obs}
is set to \code{TRUE}. Can be either \code{"blue"}, \code{"green"}, or \code{"red"}}

\item{max.col.scale}{Rarely requires adjustment. The maximum count of
observations (therefore the darkest shaded color) only used if \code{disp.obs} is
used. This allows consistency of shading between multiple plotted graphs.
It should always be at least as high as the maximum count of observations
plotted}

\item{treat.labs}{A vector of treatment labels in the same order as treatment codes.
Easiest to use treatment labels stored by \code{mb.network()}}

\item{plot.bins}{Plot time bin boundaries as vertical dashed lines. Setting \code{plot.bins=TRUE} if \code{overlay.nma}
is specified also sets x-axis ticks to time bin boundaries automatically.}

\item{...}{Arguments for \code{ggplot()} or \code{R2jags()}}
}
\description{
Plots predicted responses from a time-course MBNMA model
}
\details{
For the S3 method \code{plot()}, if \code{disp.obs} is set to \code{TRUE} it is
advisable to ensure predictions in \code{predict} are estimated using an even
sequence of time points to avoid misrepresentation of shaded densities.
Shaded counts of observations will be relative to the treatment plotted in
each panel rather than to the network reference treatment if \code{disp.obs} is
set to \code{TRUE}.
}
\section{Overlaying NMA results}{


\code{overlay.nma} indicates regions of the data (defined as "time bins") over which it may be reasonable to "lump" different
follow-up times from different studies together and assume a standard NMA model. For example:
\itemize{
\item \code{overlay.nma=c(5,10)} indicates a single NMA of studies with follow-up times \verb{>5} and \verb{<=10}
\item \code{overlay.nma=c(5,10,15)} indicates two NMAs should be performed of studies with follow-up times \verb{>5} and \verb{<=10}
of studies with follow-up times \verb{>10} and \verb{<=15}
}

When used with MBNMA (via \code{predict.mbnma()}) this allows comparison to MBNMA results over a specific range of time within each time bin.
It can be useful to assess which time-course function might be suitable when using \code{binplot()}, or to
to assess if the MBNMA predictions are in agreement with predictions from an NMA model when using \code{plot.mb.predict()}
for a specific range of time-points. This can be a general indicator of the fit of the time-course model.

However, it is important to note that the wider the range specified in \code{overlay.nma}, the more likely it is that different time-points
are included, and therefore that there is greater heterogeneity/inconsistency in the NMA model. If \code{overlay.nma} includes
several follow-up times for any study then only a single time-point will be taken (the one closest to \code{mean(overlay.nma)}).
The NMA predictions are plotted over the range specified in \code{overlay.nma} as a horizontal line, with the 95\%CrI shown by a grey
rectangle. The NMA predictions represent those for \emph{any time-points within this range} since they lump together data at
all these time-points. Predictions for treatments that are disconnected from
the network reference treatment at data points specified within \code{overlay.nma} cannot be estimated so are not included.

It is important to note that the NMA model is not necessarily the "correct" model, since it "lumps" different time-points
together and ignores potential differences in treatment effects that may arise from this. The wider the range specified in
\code{overlay.nma}, the greater the effect of "lumping" and the stronger the assumption of similarity between studies.

For an NMA model to be estimated and a corresponding prediction to be made from it, \strong{each} time bin
must include the network reference treatment (treatment=1) evaluated in at least 1 connected study in the time bin.
If a given time bin does not meet this criteria then an NMA will not be calculated for it.
}

\examples{
\donttest{
# Create an mb.network object from a dataset
copdnet <- mb.network(copd)

# Run an MBNMA model with a log-linear time-course
loglin <- mb.run(copdnet,
  fun=tloglin(pool.rate="rel", method.rate="common"),
  rho="dunif(0,1)", covar="varadj")

# Predict responses using the original dataset to estimate the network reference
#treatment response
df.ref <- copd[copd$treatment=="Placebo",]
predict <- predict(loglin, times=c(0:20), E0=0, ref.resp=df.ref)

# Plot the predicted responses with observations displayed on plot as green shading
plot(predict, disp.obs=TRUE, overlay.ref=FALSE, col="green")

# Plot the predicted responses with the median network reference treatment response overlayed
#on the plot
plot(predict, disp.obs=FALSE, overlay.ref=TRUE)

# Plot predictions from NMAs calculated between different time-points
plot(predict, overlay.nma=c(5,10), n.iter=20000)
plot(predict, overlay.nma=c(5,10,15,20), n.iter=20000)
# Time-course fit may be less well at 15-20 weeks follow-up
}

}
