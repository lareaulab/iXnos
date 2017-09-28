#' Gate flow cytometry data
#'
#' Picks out the data from the most dense region of fsc-ssc space.
#' Requires a csv in the working directory called sample-info.csv that has the
#' field "uid"; and that all the fcs files to be analysed are in the working
#' directory.
#' @param var.names A character vector of the variables collected by the flow
#' cytometer; defaults to: "fsc", "ssc", "green", "red" & "time".
#' @param flow.plots Logical - plot fsc-ssc or not
#' @return A dataframe with all the flow data for each event
#' @usage fdat <- gateFlowData()
#' @export
gateFlowData <- function(
  var.names = c("fsc", "ssc", "green", "red", "time"),
  flow.plots = TRUE
){
  # Load pertinent libraries.
  require(flowCore)
  require(flowStats)
  require(flowViz)
  # Read sample info file
  # should be a csv with a field called "uid".
  sample.info <- read.csv(
    file = "sample-info.csv",
    header = TRUE,
    stringsAsFactors = FALSE
  )
  # Read in flow data.
  flowData <- read.flowSet(
    files = list.files(pattern = "*.fcs"),
    transformation = FALSE
  )
  # The number and order of the columns is set when you export the fcs files
  # if they're different to below, this won't work.
  if (length(var.names) != length(flowCore::colnames(flowData))) {
    print(
      "The fcs files contain more variables than you've described, please supply a different vector of names to var.names."
    )
  } else {
    flowCore::colnames(flowData) <- var.names
    sampleNames(flowData) <- sample.info$uid
    print("I've loaded the data!")
    # Gating
    # For the red and green channels, if an event has a value < 0, set it to 1
    truncTrans <- truncateTransform(
      transformationId = "Truncate-transformation",
      a = 1
    )
    myTrans <- transformList(
      c(
        var.names[3],
        var.names[4]
      ),
      truncTrans
    )
    flowData.tt <- flowCore::transform(
      flowData,
      myTrans
    )
    print("I transformed the data!")
    # curv2filter gate based on fsc & ssc.
    # bwfacValue is a tuning factor
    bwfacValue <- 2.2
    c2f <- curv2Filter(
      "fsc",
      "ssc",
      bwFac = bwfacValue
    )
    c2f.results <- flowCore::filter(
      flowData.tt,
      c2f
    )
    c2f.split <- flowCore::split(
      flowData.tt,
      c2f.results
    )
    # Make some diagnostic fsc v ssc plots with the gate drawn
    if (flow.plots == TRUE){
      system("mkdir fsc-ssc-plots")
      FscSscGraphs <- function(x) {
        png(
          filename = paste(
            "fsc-ssc-plots/",
            paste(
              description(x)$GUID,
              "-fsc-ssc-curv2Filter.png",
              sep = ""
            ),
            sep = ""),
          width = 7,
          height = 7,
          units = "in",
          res = 400
        )
        print(
          xyplot(
            `ssc` ~ `fsc`,
            data = x,
            smooth = FALSE,
            filter = c2f
          )
        )
        dev.off()
      }
      fsApply(
        flowData.tt,
        FscSscGraphs
      )
      # If there are < 20 samples, make a facetted plot for ease of comparison
      if (length(list.files(pattern = "*.fcs")) < 20) {
        png(
          filename = "fsc-ssc-c2f.png",
          width = 7,
          height = 7,
          units = "in",
          res = 400
        )
        print(
          xyplot(
            `ssc` ~ `fsc`,
            data = flowData.tt,
            smooth = FALSE,
            filter = c2f
          )
        )
        dev.off()
      }
      print("I gated the data and make some graphs!")
    }
    # Take the curv2filter object and for each sample extract the most populous
    # area
    # Load some more libraries
    require(parallel)
    require(plyr)
    require(tidyverse)
    foo <- base::Filter(
      function(x) !(nrow(x) < 2),
      mclapply(
        mclapply(
          mclapply(
            c2f.split[2:length(c2f.split)],
            fsApply,
            exprs,
            simplify = FALSE,
            use.exprs = FALSE
          ),
          lapply,
          data.frame
        ),
        ldply
      )
    )
    bar <- ldply(
      Map(
        function(x, y) dplyr::mutate(x, area = y),
        foo,
        names(foo)
      )
    )
    bar$.id <- as.factor(bar$.id)
    bar$area <- as.factor(bar$area)
    write_csv(
      x = bar %>% dplyr::count(.id, area),
      path = "events-by-sample-by-area.csv"
    )
    monkey <- ldply(
      Map(
        function(x, y) dplyr::filter(x, area == y),
        dlply(
          bar,
          .(.id)
        ),
        (bar %>%
           dplyr::count(.id, area) %>%
           dplyr::filter(n==max(n)))$area
      )
    )
    write_csv(
      x = monkey %>% dplyr::count(.id, area),
      path = "post-filter-summary.csv"
    )
    # Add the sample.info
    print("I've finished filtering the data")
    fDat <- merge(
      monkey,
      sample.info,
      by.x = ".id",
      by.y = "uid"
    )
    write_csv(
      x = fDat,
      path = "gated-facs-data.csv"
    )
    library("devtools")
    capture.output(
      session_info(),
      file = "session-info-gateFlowData.txt"
    )
    print("Finished!")
    return(fDat)
  }
}
