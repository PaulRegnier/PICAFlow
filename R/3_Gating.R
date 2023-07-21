#' Gate on cells of interest
#'
#' This function allows to gate on cells of interest.
#'
#' @param flowset A string defining the `flowSet` to use. Defaults to `NULL`.
#'
#' @param sampleToPlot A numeric vector defining which samples should be plotted for gating visualization. Defaults to `NULL`.
#'
#' @param xParameter A string defining which parameter to use as x. Defaults to `NULL`.
#'
#' @param yParameter A string defining which parameter to use as y. Defaults to `NULL`.
#'
#' @param xlim A numeric vector of length 2 defining the minimum and maximum display values for the x axis. Defaults to `NULL`.
#'
#' @param ylim A numeric vector of length 2 defining the minimum and maximum display values for the y axis. Defaults to `NULL`.
#'
#' @param xGateLim A numeric vector of length 2 defining the minimum and maximum limits for the gate x axis. `-Inf` and `+Inf` values can be used. Defaults to `NULL`.
#'
#' @param yGateLim A numeric vector of length 2 defining the minimum and maximum limits for the gate y axis. `-Inf` and `+Inf` values can be used. Defaults to `NULL`.
#'
#' @param subset A boolean defining if gated cells should be actually extracted from the `flowSet`. Defaults to `FALSE`.
#'
#' @param gateName A string defining a user-friendly name for the current gate. Defaults to `NULL`.
#'
#' @param exportAllPlots A boolean defining if generated plots should be exported as PDF. If set to `TRUE`, all samples will be plotted according to the `samplesPerPage` argument. Defaults to `FALSE`.
#'
#' @param samplesPerPage An integer defining the number of samples to be plotted per PDF page. This argument is only used when `exportAllPlots = TRUE`. Defaults to `6`.
#'
#' @param recursivity A boolean defining if generated plots should be simply printed (if set to `FALSE`) or returned by the function (if set to `TRUE`). This argument is automatically set to `TRUE` when `exportAllPlots` is also set to `TRUE`. Please note that this argument should not be changed manually. Defaults to `FALSE`.
#'
#' @param inverseGating A boolean defining if the gate should rather consider cells that are outside the gate (if set to `TRUE`). Defaults to `FALSE`.
#'
#' @return If `exportAllPlots = TRUE`, the function will output PDF files containing the gating plots for each sample of the dataset. If `subset = TRUE`, the function will return a list of 2 elements named `flowset` (containing the actual gated `flowSet`) and `summary` (containing basic statistics about the gating such as the number of cells before gating, the number of cells gated and the proportion of cells gated).
#'
#' @importFrom foreach %do%
#'
#' @export

gateData = function(flowset = NULL, sampleToPlot = NULL, xParameter = NULL, yParameter = NULL, xlim = NULL, ylim = NULL, xGateLim = NULL, yGateLim = NULL, subset = FALSE, gateName = NULL, exportAllPlots = FALSE, samplesPerPage = 6, recursivity = FALSE, inverseGating = FALSE)
{
  a = NULL
  p = NULL

  xParameterOriginal = xParameter
  yParameterOriginal = yParameter

  if (length(sampleToPlot) > 1)
  {
    xParameter = flowCore::getChannelMarker(flowset[[sampleToPlot[1]]], xParameter)$desc
    yParameter = flowCore::getChannelMarker(flowset[[sampleToPlot[1]]], yParameter)$desc
  } else
  {
    xParameter = flowCore::getChannelMarker(flowset[[sampleToPlot]], xParameter)$desc
    yParameter = flowCore::getChannelMarker(flowset[[sampleToPlot]], yParameter)$desc
  }

  plot = ggcyto::ggcyto(flowset[sampleToPlot], ggplot2::aes(x = !!rlang::sym(xParameter), y = !!rlang::sym(yParameter)))
  plot = plot + ggplot2::geom_bin_2d(bins = 256)  # geom_hex() is bugging right now (28/11/2022)

  if (is.null(xlim) == FALSE & is.null(ylim) == FALSE)
  {
    myPars = ggcyto::ggcyto_par_set(limits = list(x = xlim, y = ylim))
    plot = plot + myPars
  }

  if (length(as.numeric(which(flowCore::markernames(flowset) == xParameter))) > 0)
  {
    matchingID = as.numeric(which(flowCore::markernames(flowset) == xParameter))
    xParameter = names(flowCore::markernames(flowset))[matchingID]
  }

  if (length(as.numeric(which(flowCore::markernames(flowset) == yParameter))) > 0)
  {
    matchingID = as.numeric(which(flowCore::markernames(flowset) == yParameter))
    yParameter = names(flowCore::markernames(flowset))[matchingID]
  }

  if (is.null(xGateLim) == FALSE & is.null(yGateLim) == FALSE)
  {
    namedParametersList = list(xParameter = xGateLim, yParameter = yGateLim)
    names(namedParametersList) = c(xParameter, yParameter)

    rect.g = flowCore::rectangleGate(namedParametersList, filterId = "Gate")

    if (inverseGating == TRUE)
    {
      plot = plot + ggcyto::geom_gate(rect.g) + ggcyto::geom_stats(size = 3, color = "red")

      rect.g = !rect.g

      plot = plot + ggplot2::annotate("text", x = xlim[1], y = ylim[1], hjust = 0, vjust = "bottom", label = "Complementary gate!", color = "red", size = 3)

    } else
    {
      plot = plot + ggcyto::geom_gate(rect.g) + ggcyto::geom_stats(size = 3)
    }

    if (subset == TRUE)
    {
	    # For a very curious and not understandable reason, the "filter" function from "flowCore" does not output the results in the same format when called directly (within a script) or within a function (like here in the "gateData" function). So, the "summary" function cannot work anymore in our case. I had to directly specify which slot to use ("subSet") from the filter output ("resultFilter") then manually apply "summary" on it to make it work. This should be further investigated by the flowCore package maintainer.

	    flowset = methods::as(flowset, "flowSet")

      resultFilter = flowCore::filter(flowset, rect.g)

	    resultFilter = as.list(resultFilter)

      resultFilterStats = data.frame(matrix(0, nrow = length(resultFilter), ncol = 3), stringsAsFactors = FALSE)

      foreach::foreach(a = 1:length(resultFilter)) %do%
      {
        currentData = resultFilter[[a]]
        resultFilterStats[a, ] = c(length(currentData@subSet), as.numeric(summary(currentData@subSet)["TRUE"]), as.numeric(summary(currentData@subSet)["TRUE"])/length(currentData@subSet))
      }

      colnames(resultFilterStats) = c("n", "true", "p")
      rownames(resultFilterStats) = flowWorkspace::sampleNames(flowset)

      result = flowCore::Subset(flowset, resultFilter)

      return(list(flowset = result, summary = resultFilterStats))
    }
  }

  plot = plot + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 6, color = "black"), strip.text.y = ggplot2::element_text(size = 6, color = "black"), axis.text.x = ggplot2::element_text(size = 6), axis.text.y = ggplot2::element_text(size = 6))

  if (recursivity == FALSE)
  {
    print(plot)
  } else
  {
    return(plot)
  }

  if (exportAllPlots == TRUE & samplesPerPage > 0)
  {
    pagesToExport = ceiling(length(flowset)/samplesPerPage)

    pagesSamplesList = list()
    foreach::foreach(p = 1:pagesToExport) %do%
    {
      if (p == pagesToExport)
      {
        pagesSamplesList[[p]] = c(((p - 1) * samplesPerPage) + 1, length(flowset))
      } else
      {
        pagesSamplesList[[p]] = c(((p - 1) * samplesPerPage) + 1, p * samplesPerPage)
      }
    }

    foreach::foreach(q = 1:length(pagesSamplesList)) %do%
    {
      gateData(flowset = flowset, sampleToPlot = seq(pagesSamplesList[[q]][1], pagesSamplesList[[q]][2]), xParameter = xParameterOriginal, yParameter = yParameterOriginal, xlim = xlim, ylim = ylim, xGateLim = xGateLim, yGateLim = yGateLim, exportAllPlots = FALSE, recursivity = TRUE, inverseGating = inverseGating)

      if (dir.exists(file.path("output", "3_Gating", paste("x=", xParameterOriginal, "_y=", yParameterOriginal, sep = ""))))
      {
        unlink(file.path("output", "3_Gating", paste("x=", xParameterOriginal, "_y=", yParameterOriginal, "*.*", sep = "")))
      } else
      {
        dir.create(file.path("output", "3_Gating", paste("x=", xParameterOriginal, "_y=", yParameterOriginal, sep = "")))
      }

      ggplot2::ggsave(filename = paste("x=", xParameterOriginal, "_y=", yParameterOriginal, "_samples-", pagesSamplesList[[q]][1], "-to-", pagesSamplesList[[q]][2], ".pdf", sep = ""), device = "pdf", path = file.path("output", "3_Gating", paste("x=", xParameterOriginal, "_y=", yParameterOriginal, sep = "")))
    }
  }
}

#' Export gating parameters
#'
#' This function allows to export a text file summarizing all the parameters used for the gating process.
#'
#' @param gatingParameters A list defining the parameters to export. It should be a list generated using `gateData()` function. Defaults to `NULL`.
#'
#' @return Generated text file is saved to `output > 3_Gating` directory.
#'
#' @importFrom foreach %do%
#'
#' @export

exportGatingParameters = function(gatingParameters = NULL)
{
  a = NULL

  finalGatingParameters = data.frame(matrix("", nrow = length(gatingParameters), ncol = 11))
  colnames(finalGatingParameters) = c("GateNumber", "GateName", "xParameter", "yParameter", "xlim", "ylim", "xGateLim", "yGateLim", "samplesToUse", "samplesPerPage", "inverseGating")

  foreach::foreach(a = 1:length(gatingParameters)) %do%
  {
    currentListElement = gatingParameters[[a]]

    finalGatingParameters[a, ] = c(as.character(a), names(gatingParameters)[a], currentListElement$xParameter_value, currentListElement$yParameter_value, paste(currentListElement$xlim_value, collapse = ", "), paste(currentListElement$ylim_value, collapse = ", "), paste(currentListElement$xGateLim_value, collapse = ", "), paste(currentListElement$yGateLim_value, collapse = ", "), paste(currentListElement$samplesToUse_value, collapse = ", "), as.character(currentListElement$samplesPerPage_value), as.character(currentListElement$inverseGating_value))
  }

  utils::write.table(finalGatingParameters, file.path("output", "3_Gating", "gatingParameters.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}

#' Export gating statistics
#'
#' This function allows to export all the gating statistics.
#'
#' @param totalStats A list of gating statistics generated with the `gateData()` method. Defaults to `NULL`.
#'
#' @return Generated text file is saved to `output > 3_Gating` directory.
#'
#' @export

exportGatingStatistics = function(totalStats = NULL)
{
  totalStats = as.data.frame(totalStats)
  totalStats = rbind(colnames(totalStats), totalStats)
  rownames(totalStats)[1] = "Samples"
  utils::write.table(totalStats, file.path("output", "3_Gating", "gatingStatistics.txt"), quote = FALSE, col.names = FALSE, row.names = TRUE, sep = "\t")
}
