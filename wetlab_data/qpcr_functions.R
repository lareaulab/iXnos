#' Extracts melt curve data from QPCR results files
#'
#' Reads in and reorganises the *Melt Curve Derivative Results.csv" produced by
#' the QPCR machine. If a sample info file is described, the function will
#' attempt to merge this with the resulting dataframe.
#' @param sampleInfo A dataframe describing the samples in the analysis.
#' Must include a column called "Well".
#' @return A dataframe describing for each sample/well the second derivative of
#' fluorescence with respect to temperature at each temperature in the melt
#' analysis.
#' @usage mdat <- getMeltData(sampleInfo = sample.info)
#' @export
getMeltData <- function(sampleInfo){
  suppressPackageStartupMessages(require(tidyverse))
  # This throws an error, but it's not a big deal
  md <- read_csv(
    file = list.files(
      pattern = "*Melt Curve Derivative Results*"
    ),
    col_names = TRUE,
    col_types = cols()
  )
  mdl <- gather(
    md,
    well,
    Dev2,
    A1:H12
  )
  mdl <- mdl[,-1]
  if (missing(sampleInfo)) {
    print(
      "You didn't define a sampleInfo object, so I'm just returning all the melt curves."
    )
    return(mdl)
  } else {
    return(
      merge(
        mdl,
        sample.info,
        by.x = "well",
        by.y = "Well"
        )
    )
  }
  }

#' Convenience function for low-level processing of QPCR data
#'
#' Reads in the "*Quantification Amplification Results_SYBR.csv" file output
#' by the QPCR machine and reformats it. After excluding any requested wells,
#' it then uses functions from the libraries "qpcR" & dpcR to
#' automatically find the best-fitting model from its selection, and use this
#' to derive a per-reaction/well/sample measure of the take-off value
#' (Cy0 method) and the PCR efficiency.
#'
#'
#' @param sampleInfo A dataframe describing the samples in the analysis.
#' Must include a column called "Well".
#' @param bad.ones A character vector of wells known to be empty or bad.
#' @param verbose Logical, default is TRUE. Should the function print data on
#' the per-sample curve fitting to the terminal?
#' @return 1. A plot of the qpcr model fits for each sample.
#' 2. A list containing two dataframes: "fluorescence_data", that
#' describes the raw fluorescence at each cycle for each well/sample;
#' "Cy0_and_efficiency", that describes the take-off cycle (according to the
#' Cy0 method) and efficiency for each well/sample.
#' @usage mdat <- getQPCRdata(sampleInfo = sample.info, bad.ones,
#' verbose = TRUE)
#' @export
getQPCRdata <- function(
  sampleInfo,
  bad.ones,
  verbose = TRUE
) {
  # sampleInfo object must have a col called "Well"
  suppressPackageStartupMessages(require("tidyverse"))
  results.list <- list()
  # Raw fluorescence
  qdat <- read_csv(
    file = list.files(
      pattern = "*Quantification Amplification Results_SYBR.csv"
    ),
    col_names = TRUE,
    col_types = cols()
  )[,-1]
  qdatl <- gather(
    qdat,
    unique.id,
    intensity,
    A1:H12
  )
  if(missing(sampleInfo)) {
    results.list[[1]] <- qdatl
  } else {
    results.list[[1]] <- merge(
      sampleInfo,
      qdatl,
      by.x = "Well",
      by.y = "unique.id"
      )
  }
  names(results.list)[1] <- "fluorescence_data"
  # QPCR model fitting
  require("qpcR")
  require("dpcR")
  #
  if (missing(bad.ones)) {
    l <- modlist(
      qdat,
      cyc = 1,
      fluo = NULL,
      # baseline = "median",
      remove = "none",
      opt = TRUE,
      verbose = verbose
    )
  } else {
    l <- modlist(
      qdat %>% dplyr::select(
        match(
          setdiff(
            names(qdat),
            bad.ones
          ),
          colnames(qdat)
        )
      ),
      cyc = 1,
      fluo = NULL,
      # baseline = "median",
      remove = "none",
      opt = TRUE,
      verbose = verbose
    )
  }
  #
  # png(
  #   filename = "pcr-fits.png",
  #   width = 7,
  #   height = 7,
  #   units = "in",
  #   res = 400)
  # plot(
  #   l,
  #   which = "single",
  #   confband = "confidence"
  # )
  # dev.off()
  #
  dq <- data.frame(
    well = colnames(
      qdat %>% dplyr::select(
        match(
          setdiff(
            names(qdat),
            bad.ones
          ),
          colnames(qdat)
        )
      )
    )[-1],
    qpcr_analyser(
      l,
      takeoff = FALSE
    )
  )
  if(missing(sampleInfo)) {
    results.list[[2]] <- dq
  } else {
    results.list[[2]] <- merge(
      sampleInfo,
      dq,
      by.x = "Well",
      by.y = "well"
      )
  }
  names(results.list)[2] <- "Cy0_and_efficiency"
  return(results.list)
}
