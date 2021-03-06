# Functions for manipulating/preparing MBNMA datasets
# Author: Hugo Pedder
# Date created: 2018-09-10


#' Create an mb.network object
#'
#' Creates an object of class `mb.network`. Various MBNMA functions can subsequently be applied
#' to this object.
#'
#' @param data.ab A data frame of arm-level data in "long" format containing the columns:
#' * `studyID` Study identifiers
#' * `time` Numeric data indicating follow-up times
#' * `y` Numeric data indicating the aggregate response for a given observation (e.g. mean)
#' * `se` Numeric data indicating the standard error for a given observation
#' * `treatment` Treatment identifiers (can be numeric, factor or character)
#' * `class` An optional column indicating a particular class identifier. Observations with the same treatment
#' identifier must also have the same class identifier.
#' * `N` An optional column indicating the number of participants used to calculate the
#' response at a given observation
#' @param reference A number or character (depending on the format of `treatment` within `data.ab`)
#' indicating the reference treatment in the network (i.e. those for which estimated relative treatment
#' effects estimated by the model will be compared to).
#' @param description Optional. Short description of the network.
#'
#' @details Missing values (`NA`) cannot be included in the dataset. Studies must have a baseline
#' measurement and more than a single follow-up time (unless change from baseline data are being used).
#' Data must be present for all arms within a study at each follow-up time.
#'
#' @return An object of class `mb.network` which is a list containing:
#' * `description` A short description of the network
#' * `data.ab` A data frame containing the arm-level network data (treatment identifiers will have
#' been recoded to a sequential numeric code)
#' * `treatments` A character vector indicating the treatment identifiers that correspond to the
#' new treatment codes.
#' * `classes` A character vector indicating the class identifiers (if included in the original data)
#' that correspond to the new class codes.
#'
#' @examples
#' # Using the osteoarthritis dataset
#' print(osteopain)
#'
#' # Define network
#' network <- mb.network(osteopain, description="Osteoarthritis Dataset")
#'
#' # Define network with different network reference treatment
#' network <- mb.network(osteopain, reference="Ce_200")
#'
#'
#' # Using the alogliptin dataset
#' network <- mb.network(alog_pcfb, description="Alogliptin Dataset")
#'
#' # Examine networks
#' print(network)
#' plot(network)
#'
#' @export
mb.network <- function(data.ab, reference=1, description="Network") {

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertCharacter(description, len=1, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  mb.validate.data(data.ab)

  index.data <- add_index(data.ab=data.ab, reference=reference)

  network <- index.data
  network <- c(list("description" = description), network)

  #mb.validate.network(network) # need to write a script to validate MBNMA network...should include sorting

  class(network) <- "mb.network"
  return(network)

}




#' Add follow-up time and arm indices to a dataset
#'
#' Adds follow-up time (`fups`, `fupcount`) and arm (`arms`, `narms`) indices to a dataset.
#'
#' @param data.ab A data frame of arm-level data in "long" format containing the columns:
#' * `studyID` Study identifiers
#' * `time` Numeric data indicating follow-up times
#' * `y` Numeric data indicating the mean response for a given observation
#' * `se` Numeric data indicating the standard error for a given observation
#' * `treatment` Treatment identifiers (can be numeric, factor or character)
#' * `class` An optional column indicating a particular class code. Treatments with the same identifier
#' must also have the same class code.
#'
#' @inheritParams mb.network
#'
#' @return A data frame similar to `data.ab` but with additional columns:
#' * `arm` Arm identifiers coded for each study
#' * `fupcount` Follow-up identifiers coded for each study
#' * `fups` The total number of follow-up measurements in each study
#' * `narm` The total number of arms in each study
#'
#' If `treatment` or `class` are non-numeric or non-sequential (i.e. with missing numeric codes),
#' treatments/classes in the returned data frame will be numbered and recoded to enforce sequential
#' numbering (a warning will be shown stating this).
#'
#' @examples
#' # Add indices to osteoarthritis pain dataset
#' data.ab <- add_index(osteopain)
#'
#' # Add indices to dataset using different network reference treatment
#' data.ab <- add_index(osteopain, reference=3)
#' @export
add_index <- function(data.ab, reference=1) {

  # Run Checks
  checkmate::assertDataFrame(data.ab)

  if (length(reference)>1) {
    stop("reference must be an object of length 1")
  }

  if (is.factor(data.ab$treatment)) {
    if (is.null(reference)) {
      reference <- levels(data.ab$treatment)[1]
      message(paste0("Reference treatment has automatically been set to `", reference, "`"))
    } else if (is.numeric(reference)) {
      reference <- as.character(data.ab$treatment[as.numeric(data.ab$treatment)==reference])[1]
      message(paste0("Reference treatment is `", reference, "`"))
    } else if (is.character(reference)) {
      if (is.na(match(reference, as.character(data.ab$treatment)))) {
        stop("Reference treatment specified is not a treatment given in the data")
      }
    }
  }

  if (is.numeric(data.ab$treatment)) {
    if (is.null(reference)) {
      reference <- sort(data.ab$treatment)[1]
      message(paste0("Reference treatment has automatically been set to `", reference, "`"))
    } else if (is.character(reference)) {
      stop("Reference treatment must correspond to format of treatments provided: a number corresponding to a treatment code within the data")
    } else if (is.numeric(reference)) {
      if (is.na(match(reference, data.ab$treatment))) {
        stop("Reference treatment specified is not a treatment given in the data")
      }
    }

    if (max(data.ab$treatment) != length(unique(data.ab$treatment)) |
        !all.equal(data.ab$treatment, as.integer(data.ab$treatment))
        ) {
      message("Treatments are being recoded to enforce sequential numbering")
    }
  }

  if (is.character(data.ab$treatment)) {
    if (is.null(reference)) {
      stop("Reference treatment must be specified if treatments are given as characters")
    } else if (is.numeric(reference)) {
      stop("Reference treatment must correspond to format of treatments provided: a character corresponding to a named treatment within the data")
    } else if (is.character(reference)) {
      if (is.na(match(reference, data.ab$treatment))) {
        stop("Reference treatment specified is not a treatment given in the data")
      }
    }
  }


  # Numeric data must be checked that sequence is consistent for sequential numbering
  # Factor data must be allocated codes based on factor levels
  # Character data must be allocated codes automatically (reference is only one that matters)


  #message("Treatments are being recoded to enforce sequential numbering")
  if (is.numeric(data.ab$treatment)) {
    activelist <- unique(data.ab$treatment)[-match(reference, unique(data.ab$treatment))]

  } else if (is.factor(data.ab$treatment)) {
    labs <- levels(data.ab$treatment)[levels(data.ab$treatment) %in% unique(data.ab$treatment)]
    data.ab$treatment <- factor(data.ab$treatment, labels=labs)

    activelist <- levels(data.ab$treatment)[levels(data.ab$treatment) %in% unique(data.ab$treatment)]
    activelist <- activelist[levels(data.ab$treatment)!=reference]

  } else if (is.character(data.ab$treatment)) {
    activelist <- sort(unique(data.ab$treatment))
    activelist <- activelist[activelist!=reference]
  }

  orderlist <- c(reference, as.character(activelist))
  #data.ab$treatname <- data.ab$treatment


  # Must be numeric for mb.run
  data.ab$treatment <- as.numeric(factor(data.ab$treatment,
                                    levels=orderlist)) # provide factor for sorting so that reference is #1

  #data.ab$treatment <- as.numeric(factor(data.ab$treatment,
  #                                       levels=orderlist)) # provide factor for sorting so that reference is #1


  #### Check for multiple arms of same treatment in same study and throw an error if found
  check <- data.ab %>%
    dplyr::group_by(studyID, time, treatment) %>%
    dplyr::mutate(duplicate=dplyr::n())

  duplicate.data.ab <- data.frame(check$studyID, check$duplicate)

  if (length(check$duplicate) != sum(check$duplicate)) {
    duplicateID <- unique(as.character(check$studyID[check$duplicate>1]))
    duplicateID <- paste(duplicateID, collapse="\n")
    msg <- paste0("Studies have multiple arms of the same treatment - multiple arms of the same treatment must be pooled into a single arm for studyID:\n",
                  duplicateID)
    stop(msg)
  }


  #### Add indices

  data.ab <- dplyr::arrange(data.ab, studyID, time, treatment)
  #data.ab <- arrange(data.ab, studyID, time, desc(treatment))

  #data.ab <- transform(data.ab,studyID=as.numeric(factor(studyID)))
  #data.ab <- arrange(data.ab, studyID, time, treatment)   NOT SURE THIS LINE IS NECESSARY

  # Do not run this function with pylr loaded!!
  data.ab <- data.ab %>%
    dplyr::group_by(studyID, time) %>%
    dplyr::mutate(arm = sequence(dplyr::n()))

  data.ab <- data.ab %>%
    dplyr::group_by(studyID, treatment) %>%
    dplyr::mutate(fupcount = sequence(dplyr::n()))

  data.ab <- data.ab %>%
    dplyr::group_by(studyID, treatment) %>%
    dplyr::mutate(fups=dplyr::n())

  data.ab <- data.ab %>%
    dplyr::group_by(studyID, time) %>%
    dplyr::mutate(narm=dplyr::n())


  # Store class labels and recode (if they exist in data.ab)
  if ("class" %in% names(data.ab)) {
    # Create class labels
    classdata <- data.ab[, names(data.ab) %in% c("treatment", "class")]
    classdata <- dplyr::arrange(classdata, treatment)
    classes <- as.character(unique(classdata$class))

    # Recode classes
    data.ab$class <- as.numeric(factor(data.ab$class, levels=classes))

    # Generate class key
    classkey <- unique(classdata)
    classkey$treatment <- factor(classkey$treatment, labels=orderlist)
    classkey$class <- factor(classkey$class, labels=classes)

    return(list("data.ab"=data.ab, "treatments"=orderlist, "classes"=classes, "classkey"=classkey))

  } else {
    return(list("data.ab"=data.ab, "treatments"=orderlist))
  }
}






#' Prepares data for JAGS
#'
#' Converts MBNMA data frame to a list for use in JAGS model
#'
#' @inheritParams mb.run
#' @inheritParams mb.network
#' @param class A boolean object indicating whether or not `data.ab` contains
#'   information on different classes of treatments
#' @param covstruct A character to indicate the covariance structure required for modelling correlation between
#' time points (if any), since
#' this determines some of the data. Can be either `"CS"` (compound symmetry) or `"AR1"` (autoregressive AR1).
#'
#' @return A named list of numbers, vector, matrices and arrays to be sent to
#'   JAGS. List elements are:
#'   * `y` An array of mean responses for each observation in each arm within each study
#'   * `se` An array of standard errors for each observation in each arm within each study
#'   * `time` A matrix of follow-up times within each study
#'   * `fups` A numeric vector with the number of follow-up measurements per study
#'   * `narm` A numeric vector with the number of arms per study
#'   * `NS` The total number of studies in the dataset
#'   * `NT` The total number of treatments in the dataset
#'   * `treat` A matrix of treatment codes within each study
#'   * `Nclass` Optional. The total number of classes in the dataset
#'   * `class` Optional. A matrix of class codes within each study
#'   * `classkey` Optional. A vector of class codes that correspond to treatment codes.
#'   Same length as the number of treatment codes.
#'   * `mat.triangle` Optional. A matrix with number indicating how to fill covariance
#'   matrices within the JAGS code.
#'   * `mat.order` Optional. A matrix with number indicating what order to fill
#'   covariance matrices within the JAGS code.
#'   * `timedif.0` Optional. A vector of the difference in times between the first and second
#'   follow-up time in each study.
#'
#' @examples
#' # Using the alogliptin dataset
#' network <- mb.network(alog_pcfb)
#' jagsdat <- getjagsdata(network$data.ab)
#'
#'
#' # Get JAGS data with class
#' netclass <- mb.network(goutSUA_CFBcomb)
#' jagsdat <- getjagsdata(netclass$data.ab, class=TRUE)
#'
#'
#' # Get JAGS data that allows for modelling correlation between time points
#' painnet <- mb.network(osteopain)
#' jagsdat <- getjagsdata(painnet$data.ab, rho="estimate", covstruct="AR1")
#'
#' @export
getjagsdata <- function(data.ab, class=FALSE, rho=NULL, covstruct="CS") {

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertLogical(class, len=1, null.ok=FALSE, add=argcheck)
  checkmate::assertChoice(covstruct, choices=c("CS", "AR1"), null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  df <- data.ab

  # dataset = factor indicating which data is required for the analysis (treatment univariate, treatment multivariate, class univariate, class multivariate)
  if (class==FALSE & is.null(rho)) {
    dataset <- 1
  } else if (class==TRUE & is.null(rho)) {
    dataset <- 2
  } else if (class==FALSE & !is.null(rho)) {
    dataset <- 3
  } else if (class==TRUE & !is.null(rho)) {
    dataset <- 4
  }

  if (dataset==1 | dataset==3) {
    varnames <- c("studyID", "y", "se", "treatment", "time", "arm", "fupcount")
  } else if (dataset==2 | dataset==4) {
    varnames <- c("studyID", "y", "se", "treatment", "time", "arm", "fupcount", "class")
  }
  nvar <- length(varnames)

  # Check correct variables are present
  if (!identical((varnames %in% names(df)), rep(TRUE, nvar))) {
    msg <- paste0("Variables are missing from dataset:\n",
                 paste(varnames[!(varnames %in% names(df))], collapse="\n"))
    stop(msg)
  }

  #df <- add_index(df, reference=reference)[["data"]]

  df <- dplyr::arrange(df, dplyr::desc(narm), dplyr::desc(fups), studyID, arm, time)

  df <- transform(df,studyID=as.numeric(factor(studyID, levels=as.character(unique(df$studyID)))))

  y <- array(rep(NA, max(as.numeric(df$studyID))*max(df$arm)*max(df$fupcount)),
             dim=c(max(as.numeric(df$studyID)),
                   max(df$arm),
                   max(df$fupcount)))

  se <- array(rep(NA, max(as.numeric(df$studyID))*max(df$arm)*max(df$fupcount)),
              dim=c(max(as.numeric(df$studyID)),
                    max(df$arm),
                    max(df$fupcount)))

  time <- array(rep(NA, max(as.numeric(df$studyID))*max(df$arm)*max(df$fupcount)),
                dim=c(max(as.numeric(df$studyID)),
                      max(df$arm),
                      max(df$fupcount)))

  narm <- vector()
  fups <- vector()
  NS <- max(as.numeric(df$studyID))
  NT <- max(df$treatment)

  treat <- matrix(rep(NA, max(as.numeric(df$studyID))*max(df$arm)),
                  nrow = max(as.numeric(df$studyID)), ncol = max(df$arm)
  )

  if (dataset==2 | dataset==4) {
    Nclass <- max(df$class)
    # class <- matrix(rep(NA, max(as.numeric(df$studyID))*max(df$arm)),
    #                 nrow = max(as.numeric(df$studyID)), ncol = max(df$arm)
    # )

    codes <- data.frame(df$treatment, df$class)
    codes <- dplyr::arrange(codes, df$treatment)
    classcode <- unique(codes)$df.class
  }

  for (i in 1:max(as.numeric(df$studyID))) {
    for (k in 1:max(df$arm[df$studyID==i])) {
      for (m in 1:max(df$fupcount[df$studyID==i])) {

        y[i,k,m] <- df$y[as.numeric(df$studyID)==i &
                           df$arm==k &
                           df$fupcount==m]

        se[i,k,m] <- df$se[as.numeric(df$studyID)==i &
                             df$arm==k &
                             df$fupcount==m]

        time[i,k,m] <- df$time[as.numeric(df$studyID)==i &
                                 df$arm==k &
                                 df$fupcount==m]
      }

      treat[i,k] <- max(df$treatment[as.numeric(df$studyID)==i &
                                       df$arm==k])

      #if (dataset==2 | dataset==4) {
      #  class[i,k] <- max(df$class[as.numeric(df$studyID)==i &
      #                               df$arm==k])
      #}

    }
    narm <- append(narm, max(df$arm[as.numeric(df$studyID)==i]))
    fups <- append(fups, max(df$fupcount[as.numeric(df$studyID)==i]))
  }

  time <- time[1:max(as.numeric(df$studyID)),1,]

  mat.triangle <- vector()
  p <- 1
  for (i in 2:max(fups)) {
    mat.triangle[fups==i] <- p
    p <- p+i
  }


  mat.order <- array(dim=c(max(fups), max(fups), NS))
  for (i in 1:NS) {
    p <- 1
    for (c in 1:(fups[i]-1)) {
      for (r in (c+1):fups[i]) {
        mat.order[r,c,i] <- p
        mat.order[c,r,i] <- p
        p <- p+1
      }
    }
  }


  # Return data array
  datalist <- list(y=y, se=se, time=time, fups=fups, narm=narm, NS=NS)

  # Add treatment details only if multiple treatments available
  #(allows for jags model of reference treatment)
  if (NT>1) {
    datalist[["NT"]] <- NT
    datalist[["treat"]] <- treat
  }


  if (dataset==2) {
    datalist[["Nclass"]] <- Nclass
    datalist[["class"]] <- classcode
    #datalist[["classkey"]] <- classcode
  } else if (dataset==3) {
    #datalist[["mat.triangle"]] <- mat.triangle
    #datalist[["mat.order"]] <- mat.order
  } else if (dataset==4) {
    #datalist[["mat.triangle"]] <- mat.triangle
    #datalist[["mat.order"]] <- mat.order
    datalist[["Nclass"]] <- Nclass
    datalist[["class"]] <- classcode
    #datalist[["classkey"]] <- classcode
  }

  if (!is.null(rho) & !is.null(covstruct)) {
    if (covstruct=="AR1") {
      #datalist[["mat.triangle"]] <- mat.triangle
      #datalist[["mat.order"]] <- mat.order
      datalist[["timedif.0"]] <- time[,2]-time[,1]
    }
  }

  return(datalist)

}


#' Create a dataset with the latest time point only
#'
#' Takes the latest time point from each arm in each study within an
#' `mb.network` object. Useful for network plots.
#'
#' @inheritParams mb.run
#'
#' @return A data frame in long format of responses at the latest time point in
#'   each arm of each study.
#'
#' @examples
#' # Using the alogliptin dataset
#' network <- mb.network(alog_pcfb)
#'
#' # Generate a data frame with only the latest time point included in each study
#' get.latest.time(network)
#'
#' @export
get.latest.time <- function(network) {
  # Create dataset with latest time point in studies only (useful for network plots, checking indirect evidence, etc.)
  # network is class "mb.network"

  checkmate::assertClass(network, "mb.network", null.ok=FALSE)

  data <- do.call("rbind",
                  by(network[["data.ab"]], INDICES=list(network[["data.ab"]]$studyID, network[["data.ab"]]$arm),
                     FUN=function(DF) DF[which.max(DF$time), ]))

  data <- dplyr::arrange(data, studyID, arm, time)

  return(data)
}






#' Create a dataset with the earliest time point only
#'
#' Takes the earliest time point from each arm in each study within an
#' `mb.network` object. Useful for network plots.
#'
#' @inheritParams mb.run
#'
#' @return A data frame in long format of responses at the earliest time point in
#'   each arm of each study.
get.earliest.time <- function(network) {

  checkmate::assertClass(network, "mb.network", null.ok=FALSE)

  data <- do.call("rbind",
                  by(network[["data.ab"]], INDICES=list(network[["data.ab"]]$studyID, network[["data.ab"]]$arm),
                     FUN=function(DF) DF[which.min(DF$time), ]))

  data <- dplyr::arrange(data, studyID, arm, time)

  return(data)
}









# FUNCTION IS DEPRACATED!!!
makecontrast <- function(id, treatdose, y=NULL, sd=NULL, n, r=NULL) {
  ######################################################################
  #### Make contrast converts arm-level data to contrast-level data ####
  ######################################################################

  # Variables must be sorted in terms of study ID and reference treatment

  #### Parameters:
  # id = study ID
  # treatdose = treatment
  # y = arm mean
  # sd = arm SD
  # n = arm N
  # r = arm response

  #### Returns:
  # treat1 = treatment in arm1 of contrast
  # treat2 = treatment in arm2 of contrast
  # TE = relative treatment effect (on log scale for binary data)
  # seTE = SE of relative treatment effect (on log scale for binary data)
  # n1 = N in arm1
  # n2 = N in arm2
  # studLab = study ID


  TE <- vector()
  seTE <- vector()
  treat1 <- vector()
  treat2 <- vector()
  n1 <- vector()
  n2 <- vector()
  studLab <- vector()

  for (i in seq_along(id)) {

    #print(i)

    k <- i+1

    while (k<=length(id) & id[k] == id[i] & !is.null(id[k])) {

      treat1 <- append(treat1, treatdose[i])
      treat2 <- append(treat2, treatdose[k])


      if (!is.null(y) & is.null(r)) {
        TE <- append(TE, (y[k] - y[i]))

        sigma <- ((n[i]-1)*sd[i] + (n[k]-1)*sd[k]) / (n[i] + n[k] - 2)
        se <- sigma / ((n[i] + n[k])^0.5)
        seTE <- append(seTE, se)

      } else if (is.null(y) & !is.null(r)) {
        TE <- append(TE, log(
          (r[k]/(n[k]-r[k])) / (r[i]/(n[i]-r[i]))
        ))
        variance <- (1/r[k] + 1/(n[k]-r[k]) + 1/r[i] + 1/(n[i]-r[i]))
        seTE <- append(seTE, variance^0.5)

      } else {stop("Response information is not given for either continuous or binary data")}


      n1 <- append(n1, n[i])
      n2 <- append(n2, n[k])

      studLab <- append(studLab, id[i])

      k = k+1



    }

  }

  result <- data.frame(treat1, treat2, TE, seTE, n1, n2, studLab)

  return(result)

}







#' Convert arm-based MBNMA data to contrast data
#'
#' Converts an object of class `mb.network` from arm-based long MBNMA data to a data frame with
#' contrast data (a separate contrast for each treatment comparison at each time point within each
#' study). Data can be either long or wide.
#'
#' @param network An object of class `mb.network`
#' @param datatype A string indicating the data type. Can be `binomial` or `normal`
#' @param format A string indicating the data format. Can be `wide` (two additional columns for each
#' variable - contrast arms) or `long`.
#'
#' @return A data frame with the following columns. In `wide` format, some columns are given the indices
#' 1 and 2 to indicate each arm in a given treatment comparison.:
#' * `t` The treatment in each arm
#' * `TE` The treatment effect (mean difference, log-odds) for the treatment in arm 1 versus the treatment
#' in arm 2
#' * `seTE` The standard error for the treatment effect (mean difference, log-odds) for the treatment in
#' arm 1 versus the treatment in arm 2
#' * `y` The mean response in each arm
#' * `se` The standard error of the mean in each arm
#' * `r` The number of responders in each arm
#' * `n` The total number of participants in each arm
#' * `fupcount` Follow-up identifier
#' * `time` The time the data are reported
#' * `studyID` Study identifier
#'
#' @examples
#' # Create mb.network
#' network <- mb.network(osteopain)
#'
#' # Convert to wide contrast data
#' mb.make.contrast(network, format="wide")
#'
#' # Convert to long contrast data
#' mb.make.contrast(network, format="long")
#' @export
mb.make.contrast <- function(network, datatype=NULL, format="wide") {
  # network is an object of class mb.network
  # datatype can be "binomial" or "normal"
  # format can be "wide" (two additional columns for each variable - contrast arms) or "long"
  #Wide format still has to be implemented!!

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mb.network", null.ok=FALSE, add=argcheck)
  checkmate::assertChoice(datatype, choices=c("normal", "binomial"), null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(format, choices=c("wide", "long"), null.ok=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Possibly need to add something to ensure data is sorted..?

  data <- network[["data.ab"]]

  # Checks
  if (is.null(datatype) & ("n" %in% names(data))) {
    datatype <- "binomial"
    warning("Data type not specified and data frame contains n - data are assumed to be binary")
  } else if (is.null(datatype) & ("y" %in% names(data))) {
    datatype <- "normal"
    warning("Data type not specified and data frame contains y - data are assumed to be normal")
  } else if (is.null(datatype) & ("y" %in% names(data)) & "n" %in% names(data)) {
    stop("Data type not specified and data frame contains n and y - unclear whether to use binomial or normal")
  } else if (is.null(datatype) & !("y" %in% names(data)) & !("n" %in% names(data))) {
    stop("Data type not specified and data frame does not contain n or y - unclear whether to use binomial or normal")
  }

  TE <- vector()
  seTE <- vector()
  t1 <- vector()
  t2 <- vector()
  class1 <- vector()
  class2 <- vector()
  y1 <- vector()
  y2 <- vector()
  se1 <- vector()
  se2 <- vector()
  n1 <- vector()
  n2 <- vector()
  N <- vector()
  studyID <- vector()
  fupcount <- vector()
  time <- vector()

  if (format=="long" & datatype=="normal") {
    warning("studyID has been changed to allow separate ID for each contrast rather than study")
    #longdata <- data[0,]
    longdata <- data.frame("studyID"=NA, "treatment"=NA, "time"=NA, "y"=NA, "se"=NA, "arm"=NA, "fupcount"=NA)

    if ("class" %in% names(data)) {
      longdata$class <- NA
    }

    #studyID.long <- 1
    row <- 1
  } else if (format=="long" & datatype=="binomial") {
    stop("Long format not currently supported with binomial data")
  }

  for (i in seq_along(data[["studyID"]])) {

    #print(i)

    k <- i+1

    while (k<=length(data[["studyID"]]) &
           data[["studyID"]][k] == data[["studyID"]][i] &
           data[["fupcount"]][k] == data[["fupcount"]][i] &
           !is.null(data[["studyID"]][k])) {

      # if fupcount changes then don't increase studyID.long and use studyID.long from previously
      # if fupcount remains the same then increase studyID.long

      t1 <- append(t1, data[["treatment"]][i])
      t2 <- append(t2, data[["treatment"]][k])

      if ("class" %in% names(data)) {
        class1 <- append(class1, data[["class"]][i])
        class2 <- append(class2, data[["class"]][k])
      }

      TE <- append(TE, (data[["y"]][k] - data[["y"]][i]))

      if (datatype=="binomial") {

        sigma <- (((data[["n"]][i]-1)*data[["se"]][i]*(data[["n"]][i]^0.5)) +
                    ((data[["n"]][k]-1)*data[["se"]][k]*(data[["n"]][k]^0.5))) /
          (data[["n"]][i] + data[["n"]][k] - 2)
        se <- sigma / ((data[["n"]][i] + data[["n"]][k])^0.5)
        n1 <- append(n1, data[["n"]][i])
        n2 <- append(n2, data[["n"]][k])
        N <- append(N, n1+n2)
      } else if (datatype=="normal") {
        se <- (data[["se"]][i] + data[["se"]][k])/2

        y1 <- append(y1, data[["y"]][i])
        y2 <- append(y2, data[["y"]][k])
        se1 <- append(se1, data[["se"]][i])
        se2 <- append(se2, data[["se"]][k])
      }

      seTE <- append(seTE, se)

      studyID <- append(studyID, data[["studyID"]][i])
      fupcount <- append(fupcount, data[["fupcount"]][i])
      time <- append(time, data[["time"]][i])

      if (datatype=="normal" & format=="long") {

        longdata[row,] <- rep(NA, ncol(longdata))
        longdata$y[row] <- data[["y"]][i]
        longdata$se[row] <- data[["se"]][i]
        longdata$arm[row] <- 1
        longdata$treatment[row] <- data[["treatment"]][i]
        longdata$fupcount[row] <- data[["fupcount"]][i]
        longdata$time[row] <- data[["time"]][i]
        longdata$studyID[row] <- data[["studyID"]][i]
        #longdata$test[row] <- studyID.long

        if ("class" %in% names(data)) {
          longdata$class[row] <- data[["class"]][i]
        }
        row <- row + 1

        longdata[row,] <- rep(NA, ncol(longdata))
        longdata$y[row] <- data[["y"]][k]
        longdata$se[row] <- data[["se"]][k]
        longdata$arm[row] <- 2
        longdata$treatment[row] <- data[["treatment"]][k]
        longdata$fupcount[row] <- data[["fupcount"]][k]
        longdata$time[row] <- data[["time"]][k]
        longdata$studyID[row] <- data[["studyID"]][k]
        #longdata$test[row] <- studyID.long

        if ("class" %in% names(data)) {
          longdata$class[row] <- data[["class"]][k]
        }
        row <- row + 1

        #studyID.long <- studyID.long + 1

      }

      k = k+1

    }


    #k = k+1





  }
  #print(studyID)
  #print(TE)

  if (format=="wide") {
    if (datatype=="binomial") {
      result <- data.frame("t1"=t1, "t2"=t2, "TE"=TE, "seTE"=seTE, "n1"=n1, "n2"=n2, "fupcount"=fupcount, "time"=time, "studyID"=studyID)
    } else if (datatype=="normal") {
      result <- data.frame("t1"=t1, "t2"=t2, "TE"=TE, "seTE"=seTE, "y1"=y1, "y2"=y2, "se1"=se1, "se2"=se2, "fupcount"=fupcount, "time"=time, "studyID"=studyID)
    }

    if ("class" %in% names(data)) {
      result$class1 <- class1
      result$class2 <- class2
    }
  } else if (format=="long" & datatype=="normal") {

    longdata <- longdata %>%
      dplyr::group_by(studyID, fupcount) %>%
      dplyr::mutate(contrastID = c(0,0) + rep(seq(1:(dplyr::n()/2)), each=2))

    longdata$studyID <- as.numeric(factor(paste(longdata$studyID, longdata$contrastID, sep="_"),
                                          levels=unique(paste(longdata$studyID, longdata$contrastID, sep="_"))
    ))

    result <- longdata[,!(names(longdata) %in% c("contrastID"))]
  }

  return(result)

}







#' Validates that a dataset fulfills requirements for MBNMA
#'
#' @inheritParams mb.network
#' @param CFB A boolean object to indicate if the dataset is composed of studies measuring change from
#' baseline (`TRUE`) or not (`FALSE`). It is not essential to specify this correctly but failing to do so
#' may lead to warnings.
#' @param single.arm A boolean object to indicate whether or not function should allow singe arm studies to
#' be allowed in the network without returning an error. Default is not to allow their inclusion (`single.arm=FALSE`)
#'
#'
#' @details Checks done within the validation:
#' * Checks data.ab has required column names
#' * Checks there are no NAs
#' * Checks that all SEs are positive
#' * Checks that studies have baseline measurement (unless change from baseline data is being used)
#' * Checks that arms are balanced at each time point
#' * Checks that class codes are consistent within each treatment
#' * Checks that treatment codes are consistent across different time points within a study
#' * Checks that studies have at least two arms (if `single.arm = FALSE`)
#'
#' @return An error or warnings if checks are not passed. Runs silently if checks are passed
#'
mb.validate.data <- function(data.ab, single.arm=FALSE, CFB=TRUE) {
  # data.ab must have columns c("studyID", "time", "y", "se", "treatment")
  # optional column of class

  # Checks data.ab has required column names
  # Checks there are no NAs
  # Checks that all SEs are positive
  # Checks that studies have baseline measurement (if non-CFB data is being used)
  # Checks that arms are balanced at each time point
  # Checks that studies have more than a single time point - NO LONGER INCLUDED
  # Checks that class codes are consistent within each treatment
  # Checks that studies have more than one arm if single.arm==FALSE

  varnames <- c("studyID", "time", "y", "se", "treatment")
  data.ab <- dplyr::arrange(data.ab, studyID, time, treatment)

  # Check data.ab has required column names
  if (all(varnames %in% names(data.ab))==FALSE) {
    stop("Required variable names are: 'studyID', 'time', 'y', 'se', 'treatment'")
  }

  # Check there are no NAs
  na.vars <- vector()
  for (i in seq_along(varnames)) {
    if (anyNA(data.ab[[varnames[i]]])) {
      na.vars <- append(na.vars, varnames[i])
    }
  }
  if (length(na.vars)>0) {
    stop(paste0("NA values in:\n", paste(na.vars, collapse="\n")))
  }


  # Check that all SEs are positive
  if (!all(data.ab[["se"]]>0)) {
    stop("All SEs must be >0")
  }

  # Generate narms index for checking if studies are only single-arm
  if (single.arm==FALSE) {
    data.ab <- data.ab %>%
      dplyr::group_by(studyID, time) %>%
      dplyr::mutate(narms = dplyr::n())
  }

  singlearm.studyID <- vector()
  nozero.studyID <- vector()
  onetime.studyID <- vector()
  unbalance.studyID <- vector()
  for (i in seq_along(unique(data.ab$studyID))) {
    subset <- data.ab[data.ab$studyID==unique(data.ab$studyID)[i],]

    if (CFB==FALSE) {
      # Check that studies have baseline measurement (if non-CFB data is being used)
      if (subset$time[1] != 0) {
        nozero.studyID <- append(nozero.studyID, as.character(subset$studyID[1]))
      }

      # Check that studies have more than a single time point
      if (length(unique(subset$time))<2) {  # This value can be changed (or section removed) if using CFB data
        onetime.studyID <- append(onetime.studyID, as.character(subset$studyID[1]))
      }
    }

    # Check that no studies are single arm
    if (all(subset$narms<2)) {
      singlearm.studyID <- append(singlearm.studyID, as.character(subset$studyID[1]))
    }

    # Check that arms are balanced at each time point
    studyIDtreat <- subset$treatment[subset$time == subset$time[1]]
    for (m in seq_along(unique(subset$time))) {
      if (identical(studyIDtreat, subset$treatment[subset$time == unique(subset$time)[m]]) == FALSE) {
        unbalance.studyID <- append(unbalance.studyID, as.character(subset$studyID[1]))
      }
    }



  }

  # Return printed errors
  if (CFB==FALSE) {
    if (length(nozero.studyID) >0) {
      warning(paste0("The following studies are missing a baseline measurement and\ndata have not been specified as change from baseline (`CFB=TRUE`):\n",
                  paste(nozero.studyID, collapse="\n")))
    }

    if (length(onetime.studyID)>0) {
      warning(paste0("The following studies only have data at a single time point and\ndata have not been specified as change from baseline (`CFB=TRUE`):\n",
                  paste(onetime.studyID, collapse="\n")))
    }
  }

  if (length(singlearm.studyID) >0) {
    stop(paste0("The following studies do not contain more than a single study arm:\n",
                paste(unique(singlearm.studyID), collapse="\n")))
  }
  if (length(unique(unbalance.studyID))>0) {
    stop(paste0("The following studies either do not report data for each treatment at each timepoint (data are unbalanced),\nor they do not have consistent treatment codes at different time points within the study:\n",
               paste(unique(unbalance.studyID), collapse="\n")))
  }



  # Check that class codes are consistent within each treatment
  if ("class" %in% names(data.ab)) {
    class.mismatch <- vector()
    for (i in seq_along(unique(data.ab$treatment))) {
      match <- data.ab$class[data.ab$treatment==unique(data.ab$treatment)[i]]
      if (length(unique(match)) > 1) {
        class.mismatch <- append(class.mismatch, as.character(unique(data.ab$treatment)[i]))
      }
    }
    if (length(unique(class.mismatch))>0) {
      stop(paste0("Class codes are different within the same treatment for the following treatments:\n",
                 paste(unique(class.mismatch), collapse="\n")))
    }
  }

}




#' #' Converts JAGS data held in `class("mbnma")` object to a data.frame
#' #'
#' #' @return A data frame similar in format to the data frame in an `"mb.network"` object.
#' #' @noRd
#' jagstonetwork <- function(mbnma) {
#'   checkmate::assertClass(mbnma, "mbnma")
#'
#'   jagsdat <- mbnma$model$data()
#'
#'   df <- data.frame(matrix(nrow=0, ncol=9))
#'
#'   class <- vector()
#'   for (i in 1:jagsdat$NS) {
#'     for (k in 1:jagsdat$narm[i]) {
#'       for (m in 1:jagsdat$fups[i]) {
#'         row <- c(i, jagsdat$treat[i,k], jagsdat$time[i,m], jagsdat$y[i,k,m], jagsdat$se[i,k,m],
#'                 k, m, jagsdat$fups[i], jagsdat$narm[i])
#'         df <- rbind(df, row)
#'       }
#'     }
#'   }
#'
#'   names(df) <- c("studyID", "treatment", "time", "y", "se", "arm", "fupcount", "fups", "narm")
#'
#'   if ("class" %in% names(jagsdat)) {
#'     df$class <- jagsdat$class[df$treatment]
#'   }
#'
#'   dplyr::arrange(df, fupcount, arm)
#'
#'   return(df)
#' }
