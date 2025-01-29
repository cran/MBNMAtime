## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  eval=rmarkdown::pandoc_available("1.12.3")
)

library(MBNMAtime)
library(rmarkdown)
library(knitr)
library(dplyr)
library(scales)
library(RColorBrewer)
library(ggplot2)
library(ggdist)
#load(system.file("extdata", "vignettedata.rda", package="MBNMAtime", mustWork = TRUE))

## ----warning=FALSE------------------------------------------------------------
# Loops of evidence within the alogliptin dataset
network.alog <- mb.network(alog_pcfb)
splits.alog <- mb.nodesplit.comparisons(network.alog)
print(splits.alog)

## ----eval=FALSE, results="hide"-----------------------------------------------
#  # Fit a B-spline MBNMA with a knot at 2.5 weeks and
#  #common relative effects on slope.1 and slope.2
#  mbnma <- mb.run(network.pain,
#                  fun=tspline(type="bs", knots=2.5,
#                              pool.1 = "rel", method.1="common",
#                              pool.2 = "rel", method.2="common"
#                              ))
#  
#  # Fit a UME model on both spline coefficients simultaneously
#  ume <- mb.run(network.pain,
#                  fun=tspline(type="bs", knots=2.5,
#                              pool.1 = "rel", method.1="common",
#                              pool.2 = "rel", method.2="common"
#                              ),
#                UME=TRUE)
#  
#  # Fit a UME model on the 1nd coefficient only
#  ume.slope.1 <- mb.run(network.pain,
#                  fun=tspline(type="bs", knots=2.5,
#                              pool.1 = "rel", method.1="common",
#                              pool.2 = "rel", method.2="common"
#                              ),
#                UME="beta.1")
#  
#  # Fit a UME model on the 2nd coefficient only
#  ume.slope.2 <- mb.run(network.pain,
#                  fun=tspline(type="bs", knots=2.5,
#                              pool.1 = "rel", method.1="common",
#                              pool.2 = "rel", method.2="common"
#                              ),
#                UME="beta.2")

## ----echo=FALSE---------------------------------------------------------------
print("Deviance for mbnma: 397.7")
print("Deviance for ume on beta.1 and beta.2: 386.0")
print("Deviance for ume on beta.1: 385.2")
print("Deviance for ume on beta.2: 390.1")

## ----include=FALSE------------------------------------------------------------
load(system.file("extdata", "nodesplit.rda", package="MBNMAtime", mustWork = TRUE))
load(system.file("extdata", "ns.itp.rda", package="MBNMAtime", mustWork = TRUE))

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  # Nodesplit using an Emax MBNMA
#  nodesplit <- mb.nodesplit(network.pain,
#                            fun=temax(pool.emax="rel", method.emax = "random",
#                                      pool.et50="abs", method.et50 = "common"),
#                            nodesplit.parameters="all"
#                            )

## ----echo=FALSE, eval=FALSE, include=FALSE------------------------------------
#  save(nodesplit, file="inst/extdata/nodesplit.rda")

## -----------------------------------------------------------------------------
print(nodesplit)

## ----fig.height=2.5, fig.show="hold"------------------------------------------
# Plot forest plots of direct and indirect results for each node-split comparison
plot(nodesplit, plot.type="forest")

# Plot posterior densities of direct and indirect results for each node-split comparisons
plot(nodesplit, plot.type="density")

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  # Nodesplit on emax of 1-parameter ITP MBNMA
#  ns.itp <- mb.nodesplit(network.pain,
#                         fun=titp(pool.emax = "rel", method.emax="common"),
#                         nodesplit.parameters="all")

## ----echo=FALSE, eval=FALSE, include=FALSE------------------------------------
#  save(ns.itp, file="inst/extdata/ns.itp.rda")

## -----------------------------------------------------------------------------
print(ns.itp)

plot(ns.itp, plot.type="forest")

