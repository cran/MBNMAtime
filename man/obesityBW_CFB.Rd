% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{obesityBW_CFB}
\alias{obesityBW_CFB}
\title{Studies of treatments for reducing body weight in patients with obesity}
\format{
A data frame with 710 rows and 7 variables:
\itemize{
\item \code{studyID} Study identifiers
\item \code{time} Numeric data indicating follow-up times
\item \code{y} Numeric data indicating the mean response for a given observation
\item \code{se} Numeric data indicating the standard error for a given observation
\item \code{n} Numeric data indicating the number of participants used to calculate means for each observation
\item \code{treatment} Treatment identifiers as factors. Labels are shortened treatment names.
\item \code{treatname} Character data giving the full names of each treatment in the format agent_dose
\item \code{agent} Agent (drug) names stored as characters
\item \code{class} The drug class of the agent (a broader category than \code{agent}) stored as characters
}
}
\source{
Pfizer Ltd.
}
\usage{
obesityBW_CFB
}
\description{
A dataset from a systematic review of pharmacological treatments for reducing body weight in patients with
obesity. The outcome is continuous, and aggregate data responses are given as mean change from baseline in
body weight (KG). Overall there are 35 RCTs investigating
26 treatments of 16 agents (/combinations of agents) in the network. Standard
deviations have been imputed for 421 observations.
}
\details{
\code{obesityBW_CFB} is a data frame in long format (one row per observation, arm and study),
with the variables \code{studyID}, \code{time}, \code{y}, \code{se}, \code{n}, \code{treatment}, \code{treatname}, \code{agent} and \code{class}.
}
\keyword{datasets}
