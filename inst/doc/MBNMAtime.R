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
#load(system.file("extdata", "vignettedata.rda", package="MBNMAtime", mustWork = TRUE))

## ---- echo=FALSE--------------------------------------------------------------
kable(head(osteopain), digits=2) 

## ---- echo=FALSE--------------------------------------------------------------
kable(head(alog_pcfb), digits=2) 

## ---- echo=FALSE--------------------------------------------------------------
kable(head(copd), digits=2) 

## ---- echo=FALSE--------------------------------------------------------------
kable(head(obesityBW_CFB), digits=2) 

## ---- echo=FALSE--------------------------------------------------------------
kable(head(goutSUA_CFB), digits=2) 

## ----network.pain-------------------------------------------------------------
# Using the pain dataset
network.pain <- mb.network(osteopain, reference = "Pl_0")
print(network.pain)

## -----------------------------------------------------------------------------
# Prepare data using the alogliptin dataset
network.alog <- mb.network(alog_pcfb, reference = "placebo")

# Plot network
plot(network.alog)

## ---- eval=FALSE--------------------------------------------------------------
#  plot(network.gout, level = "class", remove.loops = TRUE, label.distance = 5)

## ---- echo=FALSE--------------------------------------------------------------
suppressWarnings(plot(network.gout, level = "class", remove.loops = TRUE, label.distance = 5))

## ---- eval=FALSE--------------------------------------------------------------
#  plot(network.gout, level = "treatment", v.color = "class", label.distance = 5)

## ---- echo=FALSE--------------------------------------------------------------
suppressWarnings(plot(network.gout, level = "treatment", v.color = "class", label.distance = 5))

## ----pain.time----------------------------------------------------------------
# Prepare data using the pain dataset
network.pain <- mb.network(osteopain, reference="Pl_0")

# Draw plot of raw study responses over time
timeplot(network.pain)

## ----obese.time, message=FALSE------------------------------------------------
# Draw plot of within-study relative effects over time grouped by class
network.gout <- mb.network(goutSUA_CFBcomb)
timeplot(network.gout, level="class", plotby="rel")

## ---- results="hide"----------------------------------------------------------
# Run a linear time-course MBNMA
mbnma <- mb.run(network.alog, fun=tpoly(degree=1, pool.1="rel", method.1="common"))

## -----------------------------------------------------------------------------
summary(mbnma)

## ---- results="hide"----------------------------------------------------------
# Run an Emax time-course MBNMA with two parameters
mbnma <- mb.run(network.alog, fun=temax(
  pool.emax = "rel", method.emax="common",
  pool.et50 = "abs", method.et50="common"
))

## -----------------------------------------------------------------------------
summary(mbnma)

## ---- eval=FALSE--------------------------------------------------------------
#  # Using the COPD dataset
#  network.copd <- mb.network(copd)
#  
#  # Run an log-linear time-course MBNMA
#  # that accounts for correlation between time points using variance adjustment
#  mbnma <- mb.run(network.copd,
#                  fun=tloglin(pool.rate="rel", method.rate="random"),
#                  rho="dunif(0,1)", covar="varadj")

## ---- results="hide", message=FALSE, warning=FALSE----------------------------
# Run a B-spline time-course MBNMA with a knot at 0.2 times the max follow-up
# Common class effect on beta.2, the 2nd spline coefficient
mbnma <- mb.run(network.gout, 
                fun=tspline(type="bs", knots=c(0.2),
                            pool.1 = "rel", method.1="common",
                            pool.2="rel", method.2="random"),
                class.effect = list(beta.2="common"))


## -----------------------------------------------------------------------------
summary(mbnma)

## ---- eval=FALSE--------------------------------------------------------------
#  mbnma <- mb.run(network.copd,
#                  fun=tloglin(pool.rate="rel", method.rate="random"),
#                  priors=list(rate="dnorm(0,2) T(0,)"))

## ---- eval=FALSE--------------------------------------------------------------
#  mbnma <- mb.run(network.pain,
#                  fun=tspline(type="ls", knots=1,
#                              pool.1="rel", method.1="random",
#                              pool.2="rel", method.2="common"),
#                  omega=matrix(c(10,3,3,1), nrow=2))

## -----------------------------------------------------------------------------
print(mbnma)

## ---- eval=FALSE--------------------------------------------------------------
#  # Traceplots
#  mcmcplots::traplot(mbnma, "sd.beta.1")
#  
#  # Running mean plots
#  mcmcplots::rmeanplot(mbnma, "sd.beta.1")
#  
#  # Posterior densities
#  mcmcplots::denplot(mbnma, "sd.beta.1")
#  
#  # Autocorrelation plots
#  coda::autocorr.plot(mbnma)

## ---- results="hide", fig.show="hold", eval=FALSE-----------------------------
#  # Run a first-order fractional polynomial time-course MBNMA
#  mbnma <- mb.run(network.pain,
#                  fun=tfpoly(degree=1,
#                            pool.1="rel", method.1="random",
#                            method.power1="common"))
#  
#  # Plot a box-plot of deviance contributions (the default)
#  devplot(mbnma, n.iter=1000)

## ---- echo=FALSE, results="hide", fig.show="hold"-----------------------------
# Run a first-order fractional polynomial time-course MBNMA
mbnma <- mb.run(network.pain, 
                fun=tfpoly(degree=1,
                          pool.1="rel", method.1="random",
                          method.power1="common"), n.iter=5000)

# Plot a box-plot of deviance contributions (the default)
devplot(mbnma, n.iter=500)

## ---- eval=FALSE--------------------------------------------------------------
#  # Plot fitted and observed values with treatment labels
#  fitplot(mbnma, n.iter=1000)

## ---- results="hide"----------------------------------------------------------
# Run a quadratic time-course MBNMA using the alogliptin dataset
mbnma <- mb.run(network.alog, 
                fun=tpoly(degree=2,
                          pool.1="rel", method.1="random",
                          pool.2="rel", method.2="common"
                          )
)

plot(mbnma)

## -----------------------------------------------------------------------------
allres <- get.relative(mbnma, time=20,
                       treats = c("alog_100", "alog_50", "placebo"))
print(allres)

## ---- include=FALSE, eval=rmarkdown::pandoc_available("1.12.3")---------------
load(system.file("extdata", "ranks.rda", package="MBNMAtime", mustWork = TRUE))

## ---- results="hide", eval=rmarkdown::pandoc_available("1.12.3")--------------
# Identify quantile for knot at 1 week
timequant <- 1/max(network.pain$data.ab$time)

# Run a piecewise linear time-course MBNMA with a knot at 1 week
mbnma <- mb.run(network.pain,
                fun=tspline(type="ls", knots = timequant,
                            pool.1 = "rel", method.1="common",
                            pool.2 = "rel", method.2="common"))


# Rank results based on AUC (calculated 0-10 weeks), more negative slopes considered to be "better"
ranks <- rank(mbnma, params=c("auc", "d.2"), 
                    int.range=c(0,10),  lower_better = TRUE, n.iter=1000)

## ---- echo=FALSE, eval=FALSE, include=FALSE-----------------------------------
#  save(ranks, file="inst/extdata/ranks.rda")

## ---- eval=rmarkdown::pandoc_available("1.12.3")------------------------------
print(ranks)

## ---- eval=rmarkdown::pandoc_available("1.12.3")------------------------------
# Ranking histograms for AUC
plot(ranks, params = "auc")

## ---- eval=rmarkdown::pandoc_available("1.12.3")------------------------------
# Cumulative ranking for all ranked parameters
cumrank(ranks)

## ---- results="hide", message=FALSE, eval=FALSE-------------------------------
#  # Run an Emax time-course MBNMA using the osteoarthritis dataset
#  mbnma <- mb.run(network.pain,
#                  fun=temax(pool.emax="rel", method.emax="common",
#                            pool.et50="abs", method.et50="common"),
#                  rho="dunif(0,1)", covar="varadj")

## ---- results="hide", message=FALSE, echo=FALSE-------------------------------
# Run an Emax time-course MBNMA using the osteoarthritis dataset
mbnma <- mb.run(network.pain,
                fun=temax(pool.emax="rel", method.emax="common",
                          pool.et50="abs", method.et50="common"),
                rho="dunif(0,1)", covar="varadj", n.iter=3000)

## ---- results="hide", message=FALSE, eval=rmarkdown::pandoc_available("1.12.3")----
# Specify placebo time-course parameters
ref.params <- list(emax=-2)

# Predict responses for a selection of treatments using a stochastic E0 and
# placebo parameters defined in ref.params to estimate the network reference treatment effect
pred <- predict(mbnma, treats=c("Pl_0", "Ce_200", "Du_90", "Et_60", 
                                        "Lu_400", "Na_1000", "Ox_44", "Ro_25",
                                        "Tr_300", "Va_20"),
                        E0=~rnorm(n, 8, 0.5), ref.resp=ref.params)

print(pred)

## ---- results="hide", message=FALSE, eval=rmarkdown::pandoc_available("1.12.3")----
# Generate a dataset of network reference treatment responses over time
placebo.df <- network.pain$data.ab[network.pain$data.ab$treatment==1,]

# Predict responses for a selection of treatments using a deterministic E0 and 
#placebo.df to model the network reference treatment effect
pred <- predict(mbnma, treats=c("Pl_0", "Ce_200", "Du_90", "Et_60", 
                                        "Lu_400", "Na_1000", "Ox_44", "Ro_25",
                                        "Tr_300", "Va_20"),
                        E0=10, ref.resp=placebo.df)

print(pred)

## ---- message=FALSE, eval=rmarkdown::pandoc_available("1.12.3")---------------
plot(pred, overlay.ref=TRUE, disp.obs=TRUE)

## ---- fig.height=3, results="hide", eval=FALSE--------------------------------
#  # Fit a quadratic time-course MBNMA to the Obesity dataset
#  network.obese <- mb.network(obesityBW_CFB, reference = "plac")
#  
#  mbnma <- mb.run(network.obese,
#                  fun=tpoly(degree=2,
#                            pool.1 = "rel", method.1="common",
#                            pool.2="rel", method.2="common"))
#  
#  # Define stochastic values centred at zero for network reference treatment
#  ref.params <- list(beta.1=~rnorm(n, 0, 0.05), beta.2=~rnorm(n, 0, 0.0001))
#  
#  # Predict responses over the
#  pred.obese <- predict(mbnma, times=c(0:50), E0=100, treats = c(1,4,15),
#                          ref.resp=ref.params)
#  
#  # Plot predictions
#  plot(pred.obese, disp.obs = TRUE)

## ---- fig.height=3, results="hide", echo=FALSE, message=FALSE-----------------
# Fit a quadratic time-course MBNMA to the Obesity dataset
network.obese <- mb.network(obesityBW_CFB, reference = "plac")

mbnma <- mb.run(network.obese,
                fun=tpoly(degree=2,
                          pool.1 = "rel", method.1="common",
                          pool.2="rel", method.2="common"), n.iter=3000)

# Define stochastic values centred at zero for network reference treatment
ref.params <- list(beta.1=~rnorm(n, 0, 0.05), beta.2=~rnorm(n, 0, 0.0001))

# Predict responses over the
pred.obese <- predict(mbnma, times=c(0:50), E0=100, treats = c(1,4,15),
                        ref.resp=ref.params)

# Plot predictions
plot(pred.obese, disp.obs = TRUE)

## ---- results="hide", warning=FALSE-------------------------------------------
# Overlay predictions from lumped NMA between 8-10 weeks follow-up
plot(pred, overlay.nma=c(8,10), n.iter=20000)

## ---- warning=FALSE-----------------------------------------------------------
# Loops of evidence within the alogliptin dataset
splits.alog <- mb.nodesplit.comparisons(network.alog)
print(splits.alog)

## ---- eval=FALSE, results="hide"----------------------------------------------
#  # Identify quantile for knot at 0.5 weeks
#  timequant <- 0.5/max(network.pain$data.ab$time)
#  
#  # Fit a B-spline MBNMA with common relative effects on slope.1 and slope.2
#  mbnma <- mb.run(network.pain,
#                  fun=tspline(type="bs", knots=timequant,
#                              pool.1 = "rel", method.1="common",
#                              pool.2 = "rel", method.2="common"
#                              ))
#  
#  # Fit a UME model on both spline coefficients simultaneously
#  ume <- mb.run(network.pain,
#                  fun=tspline(type="bs", knots=timequant,
#                              pool.1 = "rel", method.1="common",
#                              pool.2 = "rel", method.2="common"
#                              ),
#                UME=TRUE)
#  
#  # Fit a UME model on the 1nd coefficient only
#  ume.slope.1 <- mb.run(network.pain,
#                  fun=tspline(type="bs", knots=timequant,
#                              pool.1 = "rel", method.1="common",
#                              pool.2 = "rel", method.2="common"
#                              ),
#                UME="beta.1")
#  
#  # Fit a UME model on the 2nd coefficient only
#  ume.slope.2 <- mb.run(network.pain,
#                  fun=tspline(type="bs", knots=timequant,
#                              pool.1 = "rel", method.1="common",
#                              pool.2 = "rel", method.2="common"
#                              ),
#                UME="beta.2")

## ---- echo=FALSE--------------------------------------------------------------
print("Deviance for mbnma: -110.54")
print("Deviance for ume on beta.1 and beta.2: -118.16")
print("Deviance for ume on beta.1: -117.51")
print("Deviance for uyme on beta.2: -118.04")

## ---- include=FALSE-----------------------------------------------------------
load(system.file("extdata", "nodesplit.rda", package="MBNMAtime", mustWork = TRUE))
load(system.file("extdata", "ns.exp.rda", package="MBNMAtime", mustWork = TRUE))

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  # Nodesplit using an Emax MBNMA
#  nodesplit <- mb.nodesplit(network.pain,
#                            fun=temax(pool.emax="rel", method.emax = "random",
#                                      pool.et50="abs", method.et50 = "common"),
#                            nodesplit.parameters="all"
#                            )

## ---- echo=FALSE, eval=FALSE, include=FALSE-----------------------------------
#  save(nodesplit, file="inst/extdata/nodesplit.rda")

## -----------------------------------------------------------------------------
print(nodesplit)

## ---- fig.height=2.5, fig.show="hold"-----------------------------------------
# Plot forest plots of direct and indirect results for each node-split comparison
plot(nodesplit, plot.type="forest")

# Plot posterior densities of direct and indirect results for each node-split comparisons
plot(nodesplit, plot.type="density")

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  # Nodesplit on emax of 1-parameter exponential MBNMA
#  ns.exp <- mb.nodesplit(network.pain,
#                         fun=texp(pool.emax = "rel", method.emax="common"),
#                         nodesplit.parameters="all")

## ---- echo=FALSE, eval=FALSE, include=FALSE-----------------------------------
#  save(ns.exp, file="inst/extdata/ns.exp.rda")

## ---- fig.height=2.5----------------------------------------------------------
print(ns.exp)

plot(ns.exp, plot.type="forest")

