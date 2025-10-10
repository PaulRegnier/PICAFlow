#' Plot densities for each parameter
#'
#' This function plots, for each parameter, the stacked densities of fluorescence of several samples on a single page.
#'
#' @param parametersToPlot A vector of all parameters to plot. Defaults to `NULL`.
#'
#' @param maxSamplesNbPerPage An integer defining the maximum number of samples to plot on each page. Defaults to `10`.
#'
#' @param folder A string defining the folder to create for the current graphs generation. Defaults to `NULL`.
#'
#' @param suffix A string defining the pattern of which `rds` files to use. Defaults to `_raw`.
#'
#' @param downsample_n An integer defining the maximum number of cells to use for each plotted sample. Defaults to `50000`.
#'
#' @param nCoresToExploit An integer defining the number of cores to use to export the graphs. It corresponds to the number of parameters to be treated concomitantly. Defaults to `NULL`, which translates to the actual number of cores of the processor minus 1.
#'
#' @return Generated plots are saved to `output > 2_Normalization > folder` directory.
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#'
#' @export

plotFacets = function(parametersToPlot = NULL, maxSamplesNbPerPage = 10, folder = NULL, suffix = "_raw", downsample_n = 50000, nCoresToExploit = NULL)
{
  p = NULL
  n = NULL
  h = NULL
  name = NULL
  i = NULL
  .data = NULL

  if (dir.exists(file.path("output", "2_Normalization", folder)))
  {
    unlink(file.path("output", "2_Normalization", folder), recursive = TRUE)

  }

  dir.create(file.path("output", "2_Normalization", folder))

  filesToOpen = dir(file.path("rds"), pattern = suffix, full.names = TRUE)

  if (is.null(nCoresToExploit) == FALSE)
  {
    coresNumber = nCoresToExploit

  } else
  {
    coresNumber = parallel::detectCores() - 1
  }

  cl = parallel::makeCluster(coresNumber, type = "PSOCK")
  doSNOW::registerDoSNOW(cl)

  pb = utils::txtProgressBar(min = 0, max = length(filesToOpen), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)

  foreach::foreach(p = filesToOpen, .packages = c("foreach", "flowCore", "ggplot2"), .options.snow = opts) %dopar%
  {
    currentData = readRDS(p)

    currentParameter = colnames(currentData[[1]]@exprs)

    currentParameter_numberID = which(parametersToPlot %in% currentParameter)

    peaksToPlot = 0

    if(file.exists(file.path("output", "2_Normalization", "peaks", paste("peaksAnalysis_Parameter-", currentParameter_numberID, "_", currentParameter, ".txt", sep = ""))) == TRUE)
    {

     currentParameter_peaksValues = utils::read.table(file = file.path("output", "2_Normalization", "peaks", paste("peaksAnalysis_Parameter-", currentParameter_numberID, "_", currentParameter, ".txt", sep = "")), sep = "\t", header = TRUE)

     rownames(currentParameter_peaksValues) = currentParameter_peaksValues$Samples
     currentParameter_peaksValues$Samples = NULL

     peaksToPlot = ncol(currentParameter_peaksValues)
    }

    nbOfPages = ceiling(length(currentData)/maxSamplesNbPerPage)

    pagesSamplesList = list()

    foreach::foreach(n = 1:nbOfPages) %do%
    {
      if (n == nbOfPages)
      {
        pagesSamplesList[[n]] = c(((n - 1) * maxSamplesNbPerPage) + 1, length(currentData))
      } else
      {
        pagesSamplesList[[n]] = c(((n - 1) * maxSamplesNbPerPage) + 1, n * maxSamplesNbPerPage)
      }
    }

    currentParameterTotalData = list()
    foreach::foreach(h = 1:length(currentData)) %do%
    {
      currentSampleData = as.data.frame(currentData[[h]]@exprs[, currentParameter])
      currentSampleName = flowWorkspace::pData(Biobase::phenoData(currentData))[h, "name"]
      currentSampleData$name = currentSampleName
      colnames(currentSampleData)[1] = currentParameter

      if (downsample_n > 0)
      {
        if (downsample_n > nrow(currentSampleData))
        {
          downsample_n = nrow(currentSampleData)
        }

        linesToKeep = sample(1:nrow(currentSampleData), downsample_n)
        currentSampleData = currentSampleData[linesToKeep, ]
      }

      currentParameterTotalData[[h]] = currentSampleData
    }

    currentParameterTotalData = do.call(rbind.data.frame, currentParameterTotalData)

    currentParameterQuantiles = stats::quantile(currentParameterTotalData[, 1], probs = seq(0, 1, 0.05))

    currentPlot_minScale = floor(as.numeric(currentParameterQuantiles[2])) - 1

    currentPlot_maxScale = ceiling(as.numeric(currentParameterQuantiles[20])) + 1

    foreach::foreach(q = 1:length(pagesSamplesList)) %do%
    {
      currentPage_lowLimit = pagesSamplesList[[q]][1]
      currentPage_highLimit = pagesSamplesList[[q]][2]
      currentPage_sampleIndices = c(currentPage_lowLimit:currentPage_highLimit)

      currentPage_sampleName = flowWorkspace::pData(Biobase::phenoData(currentData))[, "name"][currentPage_sampleIndices]



        currentParameterTotalDataTemp = currentParameterTotalData[currentParameterTotalData$name %in% currentPage_sampleName, ]



    if(peaksToPlot > 0)
    {

      currentParameterTotalDataTemp[, colnames(currentParameter_peaksValues)] = currentParameter_peaksValues[currentParameterTotalDataTemp$name, ]
    }

      grDevices::pdf(file.path("output", "2_Normalization", folder, paste("Parameter-", currentParameter_numberID, "_", currentParameter, "_samples-", currentPage_lowLimit, "-to-", currentPage_highLimit, ".pdf", sep = "")))


      plot1 = ggplot2::ggplot(currentParameterTotalDataTemp, ggplot2::aes(currentParameterTotalDataTemp[, 1], fill = name)) + ggplot2::geom_density(alpha = 1, linewidth = 0.5, na.rm = TRUE) + ggplot2::scale_x_continuous(limits = c(currentPlot_minScale, currentPlot_maxScale)) + ggplot2::facet_grid(name ~ .) + ggplot2::theme(legend.text=ggplot2::element_text(size = 6))

      if(peaksToPlot > 0)
      {
        foreach::foreach(i = 1:peaksToPlot) %do%
          {
            plot1 = plot1 + ggplot2::geom_vline(ggplot2::aes(xintercept = .data[[paste("Peak", i, sep = "")]]), color = "black", linewidth = 0.5)

          }
      }

      print(plot1)
      grDevices::dev.off()
    }
  }
}

#' Analyze the peaks for each sample and parameter
#'
#' This function generates one text file per parameter which contains all the determined peaks for every sample.
#'
#' @param parametersToAnalyze A vector of all parameters to analyze. Defaults to `NULL`.
#'
#' @param max.lms.sequence A named numeric list containing the maximum expected/desired number of peaks for each parameter. Defaults to `NULL`.
#'
#' @param suffix A string defining the pattern of which `rds` files to use. Defaults to `_raw`.
#'
#' @param samplesToDelete A character vector defining the samples to remove from the peaks analysis. Defaults to `NULL`.
#'
#' @param nCoresToExploit An integer defining the number of cores to use to analyze the peaks. It corresponds to the number of parameters to be treated concomitantly. Defaults to `NULL`, which translates to the actual number of cores of the processor minus 1.
#'
#' @param minpeakdistance An integer defining the minimum possible distance (in indices) between two peaks. The user should be careful if changing this value (see `findpeaks()` method from `pracma` package for more information). Defaults to `150`.
#'
#' @param minpeakheight An integer defining the minimum possible absolute height for a peak. The user should be careful when changing this value (see `findpeaks()` method from `pracma` package for more information). Defaults to `0.1`.
#'
#' @return Generated plots are saved to `output > 2_Normalization > peaks` directory. Generated text file is saved to `output > 2_Normalization` directory.
#'
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#'
#' @export

analyzePeaks = function(parametersToAnalyze = NULL, max.lms.sequence = NULL, suffix = NULL, samplesToDelete = NULL, nCoresToExploit = NULL, minpeakdistance = 150, minpeakheight = 0.1)
{
  a = NULL
  d = NULL
  n = NULL

  filesToOpen = dir(file.path("rds"), pattern = suffix, full.names = TRUE)

  if (is.null(nCoresToExploit) == FALSE)
  {
    coresNumber = nCoresToExploit
  } else
  {
    coresNumber = parallel::detectCores() - 1
  }

  cl = parallel::makeCluster(coresNumber, type = "PSOCK")
  doSNOW::registerDoSNOW(cl)

  pb = utils::txtProgressBar(min = 0, max = length(filesToOpen), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)

  totalPeaksList = foreach::foreach(a = filesToOpen, .packages = c("foreach", "flowCore", "pracma"), .combine = "c", .options.snow = opts) %dopar%
  {
    currentPeaksList = list()

    # Debug only: a = rds.files[1:coresNumber] a = filesToOpen[1]

    currentData = readRDS(a)

    if (length(samplesToDelete) > 0)
    {
      samplesToDeleteID = which(as.vector(flowWorkspace::pData(currentData)$name) %in% samplesToDelete)
      currentData = currentData[-samplesToDeleteID]
    }

    currentParameter = colnames(currentData[[1]]@exprs)

    currentParameterRefs = sort(currentData@phenoData[[1]])

    currentPeaksList[[currentParameter]] = data.frame(matrix(0, ncol = as.numeric(max.lms.sequence[currentParameter]), nrow = length(currentParameterRefs)), stringsAsFactors = FALSE)

    colnames(currentPeaksList[[currentParameter]]) = paste("Peak", seq(1:ncol(currentPeaksList[[currentParameter]])), sep = "")

    rownames(currentPeaksList[[currentParameter]]) = currentParameterRefs

    foreach::foreach(d = 1:length(currentParameterRefs)) %do%
    {
      currentSampleName = currentParameterRefs[d]

      currentSampleData = flowCore::exprs(currentData[[currentSampleName]])

      currentParameterDensity = stats::density(currentSampleData, n = 2^12)$y

      currentParameterResults = pracma::findpeaks(currentParameterDensity, minpeakdistance = minpeakdistance, minpeakheight = minpeakheight, sortstr = TRUE, npeaks = as.numeric(max.lms.sequence[currentParameter]))

      currentParameterPeaksIndexes = currentParameterResults[, 2]
      currentParameterPeaks = sort(stats::density(currentSampleData, n = 2^12)$x[currentParameterPeaksIndexes])

      if (length(currentParameterPeaks) != as.numeric(max.lms.sequence[currentParameter]))
      {
        currentPeaksList[[currentParameter]][currentSampleName, ] = rep(NA, as.numeric(max.lms.sequence[currentParameter]))
        currentPeaksList[[currentParameter]][currentSampleName, 1:(length(currentParameterPeaks))] = as.numeric(currentParameterPeaks)
      } else
      {
        currentPeaksList[[currentParameter]][currentSampleName, ] = as.numeric(currentParameterPeaks)
      }
    }

    return(currentPeaksList)

    gc()
  }

  parallel::stopCluster(cl)

  if (dir.exists(file.path("output", "2_Normalization", "peaks")))
  {
    unlink(file.path("output", "2_Normalization", "peaks"), recursive = TRUE)
  }

  dir.create(file.path("output", "2_Normalization", "peaks"))

  foreach::foreach(n = 1:length(totalPeaksList)) %do%
  {
    currentPeaksListOut = totalPeaksList[[n]]

    currentPeaksListOut = cbind(rownames(currentPeaksListOut), currentPeaksListOut)
    colnames(currentPeaksListOut)[1] = "Samples"
    currentPeaksName = names(totalPeaksList)[n]

    currentParameter_numberID = which(parametersToAnalyze %in% currentPeaksName)

    utils::write.table(currentPeaksListOut, file.path("output", "2_Normalization", "peaks", paste("peaksAnalysis_Parameter-", currentParameter_numberID, "_", currentPeaksName, ".txt", sep = "")), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  }
}

#' Generate new peak files using previously generated peak files
#'
#' This function generates new peak files from previously generated peak files but keeps the data already present in the old peak files. Typically, this function is used when new samples are added to the dataset but peak information from the already present samples should be kept. To use this function, please create a `refs` directory within the `output > 2_Normalization > peaks` directory, then copy the old peak files within this `refs` directory. Afterwards, generate the peak files for the current dataset using `analyzePeaks()` function, then run the `keepPreviousPeaksData()` function. Newly generated peak files will overwrite the ones present in the `output > 2_Normalization > peaks` directory. Basically, if a sample is found in both new and old peak files, corresponding data from old peak files will overwrite the corresponding data from new peak files, thus keeping intact the already determined information in the old dataset. If needed, a backup (before edition) of the peak files located in the `output > 2_Normalization > peaks` can be found in the `output > 2_Normalization > peaks > backup` directory.
#'
#' @importFrom foreach %do%
#'
#' @export

keepPreviousPeaksData = function()
{
  filesToOpen = dir(file.path("output", "2_Normalization", "peaks"), pattern = "peaksAnalysis", full.names = TRUE)

  filesRefsToOpen = dir(file.path("output", "2_Normalization", "peaks", "refs"), pattern = "peaksAnalysis", full.names = TRUE)

  if (dir.exists(file.path("output", "2_Normalization", "peaks", "backup")) == TRUE)
  {
    unlink(file.path("output", "2_Normalization", "peaks", "backup"), recursive = TRUE)
  }

  dir.create(file.path("output", "2_Normalization", "peaks", "backup"))

  totalPeaksList = list()

  foreach::foreach(q = 1:length(filesToOpen)) %do%
  {
    currentFileToOpen = filesToOpen[q]
    currentFileRefToOpen = filesRefsToOpen[q]

    currentParameterName = gsub("(.+)_Parameter-([0-9]+)_(.+)\\.txt", "\\3", currentFileToOpen)
    currentParameterNumber = gsub("(.+)_Parameter-([0-9]+)_(.+)\\.txt", "\\2", currentFileToOpen)

    currentPeakData = utils::read.table(currentFileToOpen, header = TRUE, sep = "\t")
    rownames(currentPeakData) = currentPeakData$Samples

    currentPeakDataRef = utils::read.table(currentFileRefToOpen, header = TRUE)

    newRowsID = which((currentPeakData$Samples %in% currentPeakDataRef$Samples) == FALSE)
    newRowsSamples = currentPeakData[newRowsID, "Samples"]
    refRowsID = which((currentPeakData$Samples %in% currentPeakDataRef$Samples) == TRUE)
    refRowsSamples = currentPeakData[refRowsID, "Samples"]

    mixedData_ref = currentPeakDataRef[currentPeakDataRef$Samples %in% refRowsSamples, ]
    mixedData_new = currentPeakData[currentPeakData$Samples %in% newRowsSamples, ]

    newMixedData = rbind(mixedData_ref, mixedData_new)

    currentFileToOpen_backup = file.path("output", "2_Normalization", "peaks", "backup", paste("peaksAnalysis_Parameter-", currentParameterNumber, "_", currentParameterName, ".txt", sep = ""))

    file.copy(currentFileToOpen, currentFileToOpen_backup)

    utils::write.table(newMixedData, currentFileToOpen, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  }
}

#' Synthesize peaks information into a single object
#'
#' This function takes all peak files from the `output > 2_Normalization > peaks` directory and returns a single object which contains all the useful elements for the next steps.
#'
#' @return A named list of size 3 (`info`, `best` and `raw`) containing all the summarized elements of the peak files. The `info` element contains general information about the peaks, organized as a table with the following columns: `Parameter` (parameter name), `PeaksNb` (number of peaks specified), `PeakType` (which precise peak is described) and `Mean` (the mean of all the peaks for every sample). The `best` element is a named list containing the mean peaks found for each parameter. The `raw` element stores the total peaks list for every parameter and sample.
#'
#' @importFrom foreach %do%
#'
#' @export

synthesizePeaks = function()
{
  o = NULL
  p = NULL

  filesToOpen = dir(file.path("output", "2_Normalization", "peaks"), pattern = "peaksAnalysis", full.names = TRUE)

  totalPeaksList = list()

  foreach::foreach(q = 1:length(filesToOpen)) %do%
  {
    currentFileToOpen = filesToOpen[q]
    currentParameterName = gsub("(.+)_Parameter-([0-9]+)_(.+)\\.txt", "\\3", currentFileToOpen)

    currentPeakData = utils::read.table(currentFileToOpen, header = TRUE, sep = "\t")
    rownames(currentPeakData) = currentPeakData$Samples

    currentPeakData = currentPeakData[sort(rownames(currentPeakData)), ]

    currentPeakData[, 1] = NULL

    totalPeaksList[[currentParameterName]] = currentPeakData
  }

  totalPeaksList = totalPeaksList[sort(names(totalPeaksList))]

  base.lms.list = list()

  totalPeaksInfo = data.frame(matrix(0, ncol = 4), stringsAsFactors = FALSE)

  foreach::foreach(o = 1:length(totalPeaksList)) %do%
  {
    currentParameterData = totalPeaksList[[o]]
    currentDataRownames = rownames(currentParameterData)
    currentParameterData = apply(currentParameterData, 2, as.numeric)
    rownames(currentParameterData) = currentDataRownames

    currentParameterName = names(totalPeaksList)[o]

    currentParameterPeaksInfo = data.frame(matrix(0, ncol = 4), stringsAsFactors = FALSE)

    currentPeaksList = NULL

    currentParameterMeans = colMeans(currentParameterData, na.rm = TRUE)

    foreach::foreach(p = 1:ncol(currentParameterData)) %do%
    {
      peakText = paste("Peak", p, sep = "")

      currentPeakData = currentParameterData[, p]
      currentPeakData_metric = mean(currentPeakData, na.rm = TRUE)

      currentPeaksList = c(currentPeaksList, currentPeakData_metric)

      currentParameterPeaksInfo = rbind(currentParameterPeaksInfo, c(currentParameterName, ncol(currentParameterData), peakText, currentPeakData_metric))

      rowsNA_ID = which(is.na(totalPeaksList[[o]][, p]))
      if (length(rowsNA_ID) > 0)
      {
        totalPeaksList[[o]][rowsNA_ID, p] = currentParameterMeans[p]
      }
    }

    currentParameterPeaksInfo = currentParameterPeaksInfo[-1, ]

    totalPeaksInfo = rbind(totalPeaksInfo, currentParameterPeaksInfo)

    base.lms.list[[currentParameterName]] = currentPeaksList
  }

  totalPeaksInfo = totalPeaksInfo[-1, ]

  colnames(totalPeaksInfo) = c("Parameter", "PeaksNb", "PeakType", "Mean")

  peaksAnalysis = list(info = totalPeaksInfo, best = base.lms.list, raw = totalPeaksList)

  return(peaksAnalysis)
}

# This function is slightly modified from `gaussNorm()` method (from `flowStats` package). It now supports the specification of a custom set of peaks for each parameter instead of being directly generated from this function without further control. Please note that this function should not be used as is, as it is directly called in the `normalizeData()` function.

gaussNorm2 = function(flowset, channel.names, max.lms = 2, base.lms = NULL, peak.density.thr = 0.05, peak.distance.thr = 0.05, debug = FALSE, fname = "", custom.lms = NULL)
{
  iter = NULL
  iter2 = NULL

  remb.flowset = remBoundary(flowset, channel.names)
  expr.list = c(1:length(flowset))

  if (length(max.lms) == 1)
  {
    max.lms = rep(max.lms, times = length(channel.names))
  }
  if (length(max.lms) != length(channel.names))
  {
    cat("Error: length of max.lms and channel.names doesn't match\n")
    return(NULL)
  }

  names(max.lms) = channel.names

  lms = extract.landmarks2(remb.flowset$res, expr.list, channel.names, max.lms, peak.density.thr, peak.distance.thr, custom.lms = custom.lms)
  if (is.null(base.lms))
    base.lms = extract.base.landmarks(lms$filter, channel.names, max.lms)

  ## finds the best matching between landmarks and base landmarks.  the score of a matching is defined as a function of the distance between landmark and the base landmark and the score of the landmark itself.
  if (is.null(custom.lms) == TRUE)
  {
    matched.lms = match.all.lms(lms, base.lms, channel.names, max.lms)
  } else
  {
    matched.lms = lms$filter

    foreach::foreach(iter = 1:length(matched.lms)) %do%
    {
      currentParameterCustomLms = custom.lms[[channel.names]]

      current.matched.lms = NULL

      foreach::foreach(iter2 = 1:ncol(currentParameterCustomLms)) %do%
      {
        tempDataFramePeak = c(currentParameterCustomLms[iter, iter2], base.lms[[channel.names]][iter2])

        current.matched.lms = c(current.matched.lms, tempDataFramePeak)
      }

      matched.lms[iter][[1]] = current.matched.lms
    }
  }

  confidence = compute.confidence(matched.lms, base.lms)
  cat("\nAdjusting the distance between landmarks\n")
  newRange = matrix(ncol = length(channel.names), nrow = 2)
  colnames(newRange) = channel.names
  newRange[, channel.names] = c(Inf, -Inf)

  for (i in expr.list)
  {
    cat(".")
    if (fname != "")
      file.name = paste(fname, as.character(i), sep = ".") else file.name = ""

      ## normalizing sample i.
      flowCore::exprs(remb.flowset$res[[i]]) = normalize.one.expr(flowCore::exprs(remb.flowset$res[[i]]), base.lms, lms$filter[, i], lms$original[, i], matched.lms[, i], channel.names, max.lms, file.name, debug)
      for (p in channel.names)
      {
        newRange[1, p] = min(newRange[1, p], min(flowCore::exprs(remb.flowset$res[[i]])[, p], na.rm = TRUE))
        newRange[2, p] = max(newRange[2, p], max(flowCore::exprs(remb.flowset$res[[i]])[, p], na.rm = TRUE))
      }
  }

  cat("\n")
  restoreBoundary(flowset, remb.flowset, channel.names)

  # updating the new ranges.
  for (i in expr.list)
  {
    for (p in channel.names)
    {
      ip = match(p, flowWorkspace::pData(flowCore::parameters(remb.flowset$res[[i]]))$name)
      tmp = flowCore::parameters(remb.flowset$res[[i]])
      oldRanges = unlist(range(flowset[[i]])[, p])
      flowWorkspace::pData(tmp)[ip, c("minRange", "maxRange")] = c(min(oldRanges[1], newRange[1, p]), max(oldRanges[2], newRange[2, p]))
      remb.flowset$res[[i]]@parameters = tmp

    }
  }

  list(flowset = remb.flowset$res, confidence = confidence)
}

# This function is slightly modified from `extract.landmarksNC()` method (from `flowStats` package). It now supports the specification of a custom set of peaks for each parameter instead of being directly generated from this function without further control. Please note that this function should not be used as is, as it is directly called in the `normalizeData()` function.

extract.landmarksNC2 = function(ncflowset, indices, expr.list, channel.names, max.lms, peak.density.thr, peak.distance.thr, custom.lms)
{
  ## defining output variables.
  lms.list = list()
  lms.list$original = matrix(vector("list"), length(channel.names), length(expr.list))
  lms.list$score = matrix(vector("list"), length(channel.names), length(expr.list))
  lms.list$filter = matrix(vector("list"), length(channel.names), length(expr.list))
  row.names(lms.list$original) = channel.names
  row.names(lms.list$score) = channel.names
  row.names(lms.list$filter) = channel.names
  max.lms.num = 0
  c = 1

  ## iterating over channels.
  for (p in channel.names)
  {
    if (max.lms[c] != 0)
    {
      j = 1
      for (i in expr.list)
      {
        data = flowCore::exprs(ncflowset[[i]])

        ## finding the landmarks of channel p of sample i.
        data[indices$index[[i]], ] = NA

        if (is.null(custom.lms) == TRUE)
        {
          lms = landmarker(data, p, max.lms[c])
        } else
        {
          lms = as.numeric(custom.lms[[p]][i, ])
        }

        lms.list$original[p, j] = list(lms)

        ## returns the max.lms[c] top score landmarks.
        filtered = filter.lms2(lms, data[, p], max.lms[c], peak.density.thr, peak.distance.thr, custom.lms)
        lms.list$filter[p, j] = list(filtered$lms)
        lms.list$score[p, j] = list(filtered$score)
        j = j + 1
      }
    }

    c = c + 1
  }

  return(lms.list)
}

# This function is slightly modified from `extract.landmarks()` method (from `flowStats` package). It now supports the specification of a custom set of peaks for each parameter instead of being directly generated from this function without further control. Please note that this function should not be used as is, as it is directly called in the `normalizeData()` function.

extract.landmarks2 = function(flowset, expr.list, channel.names, max.lms, peak.density.thr, peak.distance.thr, custom.lms)
{
  ## defining output variables.
  lms.list = list()
  lms.list$original = matrix(vector("list"), length(channel.names), length(expr.list))
  lms.list$score = matrix(vector("list"), length(channel.names), length(expr.list))
  lms.list$filter = matrix(vector("list"), length(channel.names), length(expr.list))
  row.names(lms.list$original) = channel.names
  row.names(lms.list$score) = channel.names
  row.names(lms.list$filter) = channel.names
  max.lms.num = 0
  c = 1

  ## iterating over channels.
  for (p in channel.names)
  {
    if (max.lms[c] != 0)
    {
      j = 1

      for (i in expr.list)
      {
        data = flowCore::exprs(flowset[[i]])

        ## finding the landmarks of channel p of sample i.
        if (is.null(custom.lms) == TRUE)
        {
          lms = landmarker(data, p, max.lms[c])
        } else
        {
          lms = as.numeric(custom.lms[[p]][i, ])
        }

        lms.list$original[p, j] = list(lms)

        ## returns the max.lms[c] top score landmarks.
        filtered = filter.lms2(lms, data[, p], max.lms[c], peak.density.thr, peak.distance.thr, custom.lms)
        lms.list$filter[p, j] = list(filtered$lms)
        lms.list$score[p, j] = list(filtered$score)
        j = j + 1
      }
    }

    c = c + 1
  }

  return(lms.list)
}


# This function is slightly modified from `filter.lms()` method (from `flowStats` package). It now supports the specification of a custom set of peaks for each parameter instead of being directly generated from this function without further control. Please note that this function should not be used as is, as it is directly called in the `normalizeData()` function.

filter.lms2 = function(lms, data, max.lms, peak.density.thr, peak.distance.thr, custom.lms)
{
  filtered = list()
  if (length(lms) == 0)
  {
    filtered$lms = vector()
    filtered$score = vector()
    return(filtered)
  }

  filtered$score = score.lms(lms, data, max.lms, peak.density.thr, peak.distance.thr)
  lms.score = data.frame(score = filtered$score, ind = c(1:length(lms)))
  lms.score = lms.score[do.call(order, c(lms.score["score"], decreasing = T)), ]

  if (is.null(custom.lms) == TRUE)
  {
    ind = which(lms.score$score > 0)
  } else
  {
    ind = which(lms.score$score >= 0)
  }

  if (length(ind) == 0)
  {
    filtered$lms = vector()
    filtered$score = vector()
    return(filtered)
  }

  lms.score.ind = lms.score$ind[ind]
  if (length(lms.score.ind) < max.lms)
    filtered$lms = sort(lms[lms.score.ind], decreasing = F) else filtered$lms = sort(lms[lms.score.ind[c(1:max.lms)]], decreasing = F)
  return(filtered)
}

#' Normalize peaks for each sample and parameter
#'
#' This function normalizes peaks for each sample and parameter, using a slightly modified version of the `gaussNorm()` method from the `flowStats` package or the unmodified `warpSet()` method from `flowStats` package. This version allows to take the peaks individually determined for each sample into account, in order to make them reach the target values for each peak of each parameter.
#'
#' @param try A boolean defining if the normalization should be performed for real or as a test. Defaults to `TRUE`.
#'
#' @param max.lms.sequence A named numeric list containing the maximum expected/desired number of peaks for each parameter. Defaults to `NULL`.
#'
#' @param suffix A string defining the `rds` files to use. Defaults to `_raw`.
#'
#' @param base.lms.list A list containing the target peaks values to reach. Should be the `best` list element from the `synthesizePeaks()` method output. Please note that this argument is only used when the normalization method is set to `gaussNorm`.
#'
#' @param warpSet A boolean defining if the normalization should use the `warpSet()` method from `flowStats` package. Defaults to `FALSE`.
#'
#' @param gaussNorm A boolean defining if the normalization should use the `gaussNorm2()` method adapted from `flowStats` package. Defaults to `TRUE`.
#'
#' @param samplesToDelete A character vector defining the samples to remove before proceeding to the normalization. Defaults to `NULL`.
#'
#' @param nCoresToExploit An integer defining the number of cores to use to export the `rds` files. Defaults to `NULL` which translates to the actual number of cores of the processor minus 1.
#'
#' @param custom.lms.list A list containing the actual peaks values for each sample and parameter. Should be the `raw` list element from the `synthesizePeaks()` method output. Please note that this argument is only used when the normalization method is set to `gaussNorm`. Defaults to `NULL`.
#'
#' @return If `try = TRUE`, the function will only return the output messages from the normalization attempt. If `try = FALSE`, generated `rds` files are saved to `rds` directory.
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#'
#' @export

normalizeData = function(try = TRUE, max.lms.sequence = NULL, suffix = NULL, base.lms.list = NULL, warpSet = FALSE, gaussNorm = TRUE, samplesToDelete = NULL, nCoresToExploit = NULL, custom.lms.list = NULL)
{
  a = NULL

  filesToOpen = dir(file.path("rds"), pattern = suffix, full.names = TRUE)

  if (is.null(nCoresToExploit) == FALSE)
  {
    coresNumber = nCoresToExploit
  } else
  {
    coresNumber = parallel::detectCores() - 1
  }

  cl = parallel::makeCluster(coresNumber, type = "PSOCK")
  doSNOW::registerDoSNOW(cl)

  pb = utils::txtProgressBar(min = 0, max = length(filesToOpen), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)

  normalizationResult = foreach::foreach(a = filesToOpen, .export = c("gaussNorm2", "extract.landmarksNC2", "extract.landmarks2", "filter.lms2", "compute.confidence", "extract.base.landmarks", "landmarker", "match.all.lms", "normalize.one.expr", "remBoundary", "restoreBoundary", "score.lms", "best.match", "register.channel", "returny", "choose.nk", "combinations.itr", "register.function", "match.score", "gau"), .packages = c("foreach", "flowCore", "attempt", "flowStats"), .combine = "c", .options.snow = opts) %dopar%
  {
    # Debug only: a = filesToOpen[1]

    currentData = readRDS(a)
    currentParameter = colnames(currentData[[1]]@exprs)

    if (length(samplesToDelete) > 0)
    {
      samplesToDeleteID = which(as.vector(flowWorkspace::pData(currentData)$name) %in% samplesToDelete)
      currentData = currentData[-samplesToDeleteID]
    }

    if (gaussNorm == TRUE)
    {
      if (try == TRUE)
      {
        print("## Testing gaussNorm normalization ##")
      } else
      {
        print("## Applying gaussNorm normalization ##")
      }

      totalErrors = 0
      totalErrorsParameters = NULL

      base.lms.used = NULL

      if (is.null(base.lms.list) == FALSE)
      {
        base.lms.used = base.lms.list[currentParameter]
      }

      custom.lms.used = NULL

      if (is.null(custom.lms.list) == FALSE)
      {
        custom.lms.used = custom.lms.list[currentParameter]

        if (length(samplesToDelete) > 0)
        {
          samplesToDeleteID = which(as.vector(rownames(custom.lms.used[[currentParameter]])) %in% samplesToDelete)
          custom.lms.used[[currentParameter]] = custom.lms.used[[currentParameter]][-samplesToDeleteID, ]
        }
      }

      if (try == TRUE)
      {
        sink("temp")

        res = attempt::silent_attempt(gaussNorm2(currentData, currentParameter, max.lms = as.numeric(max.lms.sequence[currentParameter]), base.lms = base.lms.used, custom.lms = custom.lms.used)$flowset, msg = "")
        sink()

        if (attempt::is_try_error(res) == TRUE)
        {
          textRes = paste("Parameter: ", currentParameter, " => Error during normalization! Please edit the max.lms value accordingly.", sep = "")
        } else
        {
          textRes = paste("Parameter: ", currentParameter, " => No error during normalization.", sep = "")
        }

        return(textRes)

      } else
      {
        sink("temp")
        currentData = gaussNorm2(currentData, currentParameter, max.lms = as.numeric(max.lms.sequence[currentParameter]), base.lms = base.lms.used, custom.lms = custom.lms.used)$flowset
        sink()
        saveRDS(currentData, file.path("rds", paste("step2_", currentParameter, "_normalized.rds", sep = "")))
      }

      gc()
    } else if (warpSet == TRUE)
    {
      print("## Applying warpSet normalization ##")

      sink("temp")
      currentData = flowStats::warpSet(currentData, currentParameter, peakNr = as.numeric(max.lms.sequence[currentParameter]), nbreaks = 20)
      sink()
      saveRDS(currentData, file.path("rds", paste("step2_", currentParameter, "_normalized.rds", sep = "")))

      gc()
    } else
    {
    }

    gc()
  }

  parallel::stopCluster(cl)

  if (try == TRUE & gaussNorm == TRUE)
  {
    return(normalizationResult)
  }
}

#' Merge all parameters into a single file
#'
#' This function merges all the parameter-restricted `rds` files into a single `flowSet` contained in a `rds` file. Please note that once the merging process finished, input `rds` files will be deleted.
#'
#' @param suffix A string defining the `rds` files to use. Defaults to `_normalized`.
#'
#' @return Generated `rds` file is saved to `rds` directory.
#'
#' @importFrom foreach %do%
#'
#' @export

mergeParameters = function(suffix = "_normalized")
{
  a = NULL
  b = NULL

  filesToOpen = dir(file.path("rds"), pattern = suffix, full.names = TRUE)

  print("Merging data...")

  totalParameters = NULL
  pooledData = list()
  foreach::foreach(a = 1:length(filesToOpen)) %do%
  {
    currentData = readRDS(filesToOpen[a])

    currentData = flowCore::flowSet_to_list(currentData)

    currentParameter = colnames(currentData[[1]]@exprs)

    if (a == 1)
    {
      pooledData = currentData
    } else
    {
      foreach::foreach(b = 1:length(pooledData)) %do%
      {
        pooledData[[b]] = flowCore::fr_append_cols(pooledData[[b]], currentData[[b]]@exprs)
      }
    }

    print(paste("File ", a, "/", length(filesToOpen), sep = ""))

    totalParameters = c(totalParameters, currentData[[1]]@parameters@data$desc)
  }

  foreach::foreach(b = 1:length(pooledData)) %do%
  {
    pooledData[[b]]@parameters@data$desc = totalParameters
  }

  pooledData = methods::as(pooledData, "flowSet")
  saveRDS(pooledData, file.path("rds", paste("2_Normalized.rds", sep = "")))

  unlink(dir(file.path("rds"), pattern = "step1_", full.names = TRUE))
  unlink(dir(file.path("rds"), pattern = "step2_", full.names = TRUE))

  gc()
}

# The following functions are unmodified from flowStats package (gaussNorm.R more precisely), but they need to be integrated to the package because import of unexported functions from another package is deprecated. I actually did not find any other way to include these functions without triggering any note.

## computes the confidence value based on the matched landmarks.
compute.confidence = function(matched.lms, base.lms)
{
  confidence = rep(1, times = length(matched.lms[1, ]))
  # for each channel c.
  for (c in 1:length(base.lms))
  {
    for (l in 1:length(base.lms[[c]]))
    {
      lms.list = 0
      for (i in 1:length(matched.lms[1, ]))
      {
        ind = which(matched.lms[c, i][[1]] == base.lms[[c]][l])
        if (length(ind) == 0)
        {
          lms.list[i] = NA
        } else
        {
          if (ind[1]%%2 != 0)
            ind[1] = ind[1] + 1
          lms.list[i] = matched.lms[c, i][[1]][ind[1] - 1]
        }
      }
      d = stats::density(stats::na.omit(lms.list))
      den.lms = 0
      for (i in 1:length(lms.list))
      {
        if (is.na(lms.list[i]))
          den.lms[i] = NA else den.lms[i] = d$y[which(abs(d$x - lms.list[i]) == min(abs(d$x - lms.list[i])))]
      }
      m.den.lms = mean(den.lms, na.rm = T)
      ind1 = which(den.lms > m.den.lms/4 & den.lms < m.den.lms/3)
      ind2 = which(den.lms < m.den.lms/4)
      if (length(ind1) != 0)
        confidence[ind1] = confidence[ind1] * 0.8
      if (length(ind2) != 0)
        confidence[ind2] = confidence[ind2] * 0.7
      if (length(which(is.na(lms.list))) != 0)
        confidence[which(is.na(lms.list))] = confidence[which(is.na(lms.list))] * 0.6

      # lms.frame=data.frame(lms=lms.list, ind=c(1:length(lms.list))) lms.frame=lms.frame[do.call(order, c(lms.frame['lms'], decreasing=F)), ] diff=unlist(lms.frame['lms'])[2:length(lms.list)]-unlist(lms.frame['lms'])[1:(length(lms.list)-1)] bw=2*mean(na.omit(diff)) den=0 for(i in 1:length(lms.list)){ if(is.na(lms.list[i])) den[i]=NA else den[i]=length(which(lms.list>lms.list[i]-bw & lms.list<lms.list[i]+bw))/length(lms.list) }
    }
  }
  confidence
}

## computes the base landmarks.  output: a vector of base landmarks for each channel p.
extract.base.landmarks = function(lms, channel.names, max.lms)
{
  lms.list = list()
  max.lms.num = 0
  for (p in channel.names)
  {
    if (max.lms[p] != 0)
    {
      ## if the total number of landmarks for channel p is less than max.lms[p] return them as base landmarks.
      if (length(unlist(lms[p, ])) <= max.lms[p])
      {
        lms.list[[p]] = sort(unlist(lms[p, ]))
      } else
      {

        if (max.lms[p] == 1)
        {
          lms.list[[p]] = stats::median(unlist(lms[p, ]))
        } else
        {
          ## first identify samples that have exactly max.lms[p] landmarks on their channel p. These landmarks are labeled from 1 to max.lms the landmarks samples with less than max.lms[p] landmarks are stored in short.lms vector.
          short.lms = vector()
          lms.class = list()
          for (jj in 1:max.lms[p]) lms.class[[jj]] = vector()
          for (ii in 1:length(lms[p, ]))
          {
            lms.p.ii = lms[p, ii][[1]]
            if (length(lms.p.ii) == max.lms[p])
            {
              for (jj in 1:max.lms[p]) lms.class[[jj]] = append(lms.class[[jj]], lms.p.ii[jj])
            } else
            {
              short.lms = append(short.lms, lms.p.ii)
            }
          }

          if (length(lms.class[[1]]) == 0)
          {
            cat("No curve with ", max.lms[p], " landmarks was found for channel ", p, ". Decrease max.lms[", p, "] and try again.\n")
            stop()
          }
          mean.lms.class = 0
          for (jj in 1:max.lms[p])
          {
            mean.lms.class[jj] = mean(lms.class[[jj]])
          }

          ## assigning short.lms landmarks to the class of labeled landmarks wich are closer to.
          if (length(short.lms) != 0)
          {
            for (jj in 1:length(short.lms))
            {
              kk = which(abs(mean.lms.class - short.lms[jj]) == min(abs(mean.lms.class - short.lms[jj])))
              kk = kk[1]
              lms.class[[kk]] = append(lms.class[[kk]], short.lms[jj])
            }
          }
          lms.class.len = 0
          lms.class.med = 0
          for (jj in 1:max.lms[p])
          {
            lms.class.len[jj] = length(lms.class[[jj]])
            lms.class.med[jj] = stats::median(lms.class[[jj]], na.rm = T)
          }
          s = stats::sd(lms.class.len, na.rm = T)
          m = mean(lms.class.len, na.rm = T)
          ll = which(m - lms.class.len > 2 * s)
          ## if a class of landmarks has too few landmarks in it just ignore this class.
          if (length(ll) != 0)
          {
            cat("warning: fewer landmark classes found in channel ", p, "\n")
            tmp.lms = list()
            kk = 1
            for (jj in (1:max.lms[p])[-ll])
            {
              tmp.lms[[kk]] = lms.class[[jj]]
              kk = kk + 1
            }
            lms.class = tmp.lms

          }

          ## returning the median of each class as the base landmark.
          for (jj in 1:length(lms.class))
          {
            lms.class.med[jj] = stats::median(lms.class[[jj]], na.rm = T)
          }

          lms.list[[p]] = lms.class.med

        }
      }
    }
  }
  return(lms.list)
}

## returns the peaks (local maxima's) in the kernel density estimate of data.
landmarker = function(data, channel.name, max.lms, span = 3)
{
  d = data[, channel.name]
  A = stats::density(stats::na.omit(d))
  d = A$y
  pks = c()
  ## sliding a window of size span and returns locations where the middle point in the window is maximum.
  for (i in 1:(length(d) - span))
  {
    if (!is.na(d[i + span%/%2]) & (d[i + span%/%2] == max(d[c(i:(i + span))], na.rm = T)))
      if (!is.na(d[i]) & !is.na(d[i + span]) & d[i + span%/%2] != d[i] & d[i + span%/%2] != d[i + span])
        pks = append(pks, i + span%/%2)

  }
  return(A$x[pks])
}

## finds the best matching between landmarks (lms) and the base landmarks (base.lms).  the score of a matching is defined as a function of the distance between the base landmark and the landmark and the score of the landmark itself.  returns: matched.lms-for each channel p of a sample i matched.lms[p, i] is a list of pairs of the form (lms, base.lms)
match.all.lms = function(lms, base.lms, channel.names, max.lms)
{
  n = length(lms$filter[channel.names[1], ])  ##number of samples.
  matched.lms = matrix(vector("list"), length(channel.names), n)
  row.names(matched.lms) = channel.names
  lms.class = list()
  lms.class.median = list()
  for (p in channel.names)
  {
    lms.class[[p]] = list()
    lms.class.median[[p]] = list()
  }
  for (p in channel.names)
  {
    if (max.lms[[p]] == 0)
      next
    for (i in 1:n)
    {
      matched.lms[p, i][[1]] = best.match(lms$original[p, i][[1]], lms$score[p, i][[1]], base.lms[[p]], max.lms[[p]])
    }
  }
  return(matched.lms)
}

## shifts the data in such a way that the peak at matched.lms[i] is moved to matched.lms[i+1] for each i.  base.lms, lms.lis and lms.original are passed for the debug purposes only.
normalize.one.expr = function(data, base.lms, lms.list, lms.original, matched.lms, channel.names, max.lms, fname = "", debug = FALSE)
{
  norm.data = data
  if (debug == TRUE)
  {
    if (fname == "")
      grDevices::dev.new() else grDevices::png(filename = fname, bg = "white", width = 1800, height = 1000)
    graphics::par(mfrow = c(1, length(which(max.lms != 0))))
  }
  cn = 1
  i = 0
  ## normalizing each channel.
  for (p in channel.names)
  {
    i = i + 1
    if (max.lms[i] == 0)
    {
      norm.data[, p] = data[, p]
    } else
    {
      ## registering the data of channel p given the matching landmarks.
      norm.data[, p] = register.channel(data[, p], matched.lms[[p]])
      ## drawing the pre- and post- normalization density curves.
      if (debug == TRUE)
      {
        if (length(lms.list[[p]]) == 0)
        {
          A1 = stats::density(stats::na.omit(data[, p]))
          xlim = c(min(A1$x, na.rm = T), max(A1$x, na.rm = T))
          ylim = c(min(A1$y, na.rm = T), max(A1$y, na.rm = T))
          graphics::par(mfg = c(1, cn))
          plot(A1, type = "l", col = "black", xlim = xlim, ylim = ylim, xlab = p, main = "No lms is found")
        } else
        {
          A1 = stats::density(stats::na.omit(data[, p]))
          A2 = stats::density(stats::na.omit(norm.data[, p]))
          xlim = c(min(min(A1$x, na.rm = T), min(A2$x, na.rm = T)), max(max(A1$x, na.rm = T), max(A2$x, na.rm = T)))
          ylim = c(min(min(A1$y, na.rm = T), min(A2$y, na.rm = T)), max(max(A1$y, na.rm = T), max(A2$y, na.rm = T)))
          graphics::par(mfg = c(1, cn))
          plot(A1, type = "l", col = "blue", xlim = xlim, ylim = ylim, xlab = p, main = "")
          Lms = matched.lms[[p]][seq(1, length(matched.lms[[p]]), by = 2)]
          M.Lms = matched.lms[[p]][seq(2, length(matched.lms[[p]]), by = 2)]
          M.Lms.index = 1
          for (k in 1:(length(M.Lms))) M.Lms.index[k] = which(base.lms[[p]] == M.Lms[k])
          graphics::points(Lms, returny(A1, Lms), pch = 19, col = "red")
          graphics::text(Lms, returny(A1, Lms), as.character(M.Lms.index), pos = 3, col = "red")
          graphics::points(lms.original[[p]], returny(A1, lms.original[[p]]), pch = 21, col = "black")  ##all the peaks
          graphics::par(mfg = c(1, cn))
          plot(A2, type = "l", col = "red", xlim = xlim, ylim = ylim, xlab = p)
          graphics::points(M.Lms, returny(A2, M.Lms), pch = 15, col = "blue")
          graphics::text(M.Lms, returny(A2, M.Lms), as.character(M.Lms.index), pos = 3, col = "blue")
        }
        cn = cn + 1
      }
    }
  }
  if (debug)
  {
    if (fname == "")
      readline()
    grDevices::dev.off()
  }
  return(norm.data)

}

# setting the boundary values to NA. Also returns the indices of the NA values for recovery.
remBoundary = function(flowset, channel.names)
{
  ranges = flowCore::fsApply(flowset, range)
  index = list()
  res = list()
  for (i in 1:length(flowset))
  {
    d = Biobase::exprs(flowset[[i]])
    index[[i]] = vector()
    for (p in channel.names)
    {
      ind = which(d[, p] == ranges[[i]][1, p])
      index[[i]] = append(index[[i]], ind)
      d[ind, p] = NA
      ind = which(d[, p] == ranges[[i]][2, p])
      index[[i]] = append(index[[i]], ind)
      d[ind, p] = NA
    }
    res[[i]] = methods::new("flowFrame", exprs = d)
    flowCore::parameters(res[[i]]) = flowCore::parameters(flowset[[i]])
  }
  names(res) = flowWorkspace::sampleNames(flowset)
  res = methods::as(res, "flowSet")
  Biobase::phenoData(res) = Biobase::phenoData(flowset)
  return(list(index = index, res = res))
}

# Restore the boundary values that were set to NA.
restoreBoundary = function(org.flowset, remb.flowset, channel.names)
{
  for (i in 1:length(org.flowset))
  {
    for (p in channel.names)
    {
      Biobase::exprs(remb.flowset$res[[i]])[remb.flowset$index[[i]], p] = Biobase::exprs(org.flowset[[i]])[remb.flowset$index[[i]], p]
    }
  }
}

## assigns a score to each landmark. the score of a landmarks is a funciton of its sharpness and density value.  the peaks with density value less than peak.density.thr*maximum.peak.density are discarded.  of the peaks with distance less than peak.distance.thr*range.data only one is considered.
score.lms = function(lms, data, max.lms, peak.density.thr, peak.distance.thr)
{
  bw = 64
  score = vector()
  height.cutoff = peak.density.thr
  if (length(lms) == 0)
    return(score)
  A = stats::density(stats::na.omit(data))
  bw = min(64, length(A$x)/10)
  lms.max.height = max(returny(A, lms), na.rm = T)
  MIN.LMS.DIST = (max(A$x, na.rm = T) - min(A$x, na.rm = T)) * peak.distance.thr
  last.lms = -1
  last.lms.i = -1
  last.lms.score = 0
  for (i in 1:length(lms))
  {
    lms.ind = which(stats::na.omit(abs(A$x - lms[i])) == min(stats::na.omit(abs(A$x - lms[i]))))
    ind = (max(lms.ind - bw%/%2, 1)):(min(lms.ind + bw%/%2, length(A$x)))
    if (length(ind) == 0)
      ind = 1

    if (A$y[lms.ind] < height.cutoff * lms.max.height)
    {
      w = 0
    } else
    {
      ## computing the sharpness
      w = A$y[lms.ind] - A$y[ind]
      w[which(w < 0)] = 3 * w[which(w < 0)]
    }
    ## computing final score
    score[i] = sum(w, na.rm = T) * A$y[lms.ind]
    if (score[i] < 0)
      score[i] = 0
    if (last.lms < 0)
    {
      last.lms = lms[i]
      last.lms.i = i

    } else
    {
      ## If two lms's are very close only choose one with the better score.
      if (lms[i] - last.lms < MIN.LMS.DIST)
      {
        if (score[i] > score[last.lms.i])
        {
          last.lms = lms[i]
          score[last.lms.i] = 0
          last.lms.i = i
        } else
        {
          score[i] = 0
        }
      } else
      {
        last.lms = lms[i]
        last.lms.i = i
      }
    }
  }
  return(score)
}

best.match = function(lms, lms.score, base.lms, max.lms)
{
  comb = choose.nk(max(length(lms), length(base.lms)), min(length(lms), length(base.lms)))
  sc = list()
  k = 1
  max.s = -1
  if (length(lms) < length(base.lms))
  {
    for (i in 1:(dim(comb)[1]))
    {
      if (length(which(comb[i, ] == 0)) != 0)
        c = comb[i, ][-which(comb[i, ] == 0)] else c = comb[i, ]
        d = combinations.itr(length(comb[i, ]), length(c))
        for (j in 1:(dim(d)[1]))
        {
          s = match.score(lms[d[j, ]], base.lms[c], lms.score[d[j, ]])
          if (max.s < s)
          {
            lms.index = d[j, ]
            base.lms.index = c
            max.s = s
          }
          sc[[k]] = list(score = s, lms.index = d[j, ], base.lms.index = c)
          k = k + 1
        }
    }

  } else
  {
    for (i in 1:(dim(comb)[1]))
    {
      if (length(which(comb[i, ] == 0)) != 0)
        c = comb[i, ][-which(comb[i, ] == 0)] else c = comb[i, ]
        d = combinations.itr(length(comb[i, ]), length(c))
        for (j in 1:(dim(d)[1]))
        {
          s = match.score(lms[c], base.lms[d[j, ]], lms.score[c])
          if (max.s < s)
          {
            lms.index = c
            base.lms.index = d[j, ]
            max.s = s
          }

          k = k + 1
        }
    }
  }
  res = 1:(2 * length(lms.index))
  res[seq(1, 2 * length(lms.index), by = 2)] = lms[lms.index]
  res[seq(2, 2 * length(lms.index), by = 2)] = base.lms[base.lms.index]
  return(res)
}

## manipulates the data in such a way that the landmark at matched.lms[i] is moved to matched.lms[i+1] for each i.
register.channel = function(data, matched.lms)
{
  if (length(matched.lms) == 0)
  {
    cat("*")
    return(data)
  }
  s = m = shift = vector()
  lms = vector()
  for (i in seq(1, length(matched.lms), by = 2))
  {
    shift = append(shift, matched.lms[i + 1] - matched.lms[i])
    lms = append(lms, matched.lms[i])
    s = append(s, stats::sd(stats::na.omit(data)))
  }
  r.data = register.function(data, s, lms, shift)
  return(r.data)
}

returny = function(A, X)
{
  y = vector()
  i = 1
  for (x in X)
  {
    y[i] = A$y[which(abs(A$x - x) == min(abs(A$x - x)))[1]]
    i = i + 1
  }
  return(y)
}

choose.nk = function(n, k)
{
  v = matrix(0, ncol = k, nrow = 0)
  for (j in 1:k)
  {
    c = combinations.itr(n, j)
    if (j < k)
      c = cbind(c, matrix(0, nrow = dim(c)[1], ncol = k - j))
    v = rbind(v, c)
  }
  return(v)
}

combinations.itr = function(n, k)
{
  v = 1
  res = NULL
  while (length(v) != 0)
  {
    if (v[length(v)] > n)
    {
      v = v[-length(v)]
      v[length(v)] = v[length(v)] + 1
    } else if (length(v) == k)
    {

      res = rbind(res, v)
      v[length(v)] = v[length(v)] + 1
    } else
    {
      v = append(v, v[length(v)] + 1)
    }
  }
  return(res)
}

register.function = function(data, s, m, shift)
{
  sum = 0
  if (length(m) == 1)
  {
    return(data + shift)
  }
  if (length(m) == 2)
  {
    sh = (shift[1] - shift[2])
    data = data + gau(data, s[1], m[1]) * (sh/2)
    data = data - gau(data, s[2], m[2]) * (sh/2)
    return(data + shift[1] - sh/2)
  }
  max.shift = which(abs(shift) == max(abs(shift)))[1]
  if (shift[max.shift] > 0)
  {
    sh = (shift[max.shift] - (shift[max.shift] - min(shift[-max.shift]))/2)
  } else
  {
    sh = (shift[max.shift] - (shift[max.shift] - max(shift[-max.shift]))/2)
  }
  data = data + sh
  shift = shift - sh

  for (i in 1:length(m)) data = data + gau(data, s[i], m[i]) * shift[i]
  return(data)
}

## defines the score of a matching between landmarks (lms) and base landmarks (base.lms)
match.score = function(lms, base.lms, lms.score)
{
  d = abs(lms - base.lms)
  if (length(which(d == 0)) != 0)
    d[which(d == 0)] = 0.01
  return(sum(1/d * lms.score * lms.score))
}

## gaussian function used in shifting the data.
gau = function(d, s, m)
{
  return(2.7182^(-((d - m)^2)/(2 * s * s)))
}

#' Manually fine tune peaks
#'
#' This function allows to manually fine tune the pre-determined peaks for each sample and parameter using a dedicated Shiny application.
#'
#' @importFrom foreach %do%
#' @import flowCore
#'
#' @export

fineTunePeaks = function()
{
  peaks = NULL
  p = NULL
  x = NULL
  y = NULL

  # RDS file to store the whole peaks information

  RDS_filePath = file.path("output", "2_Normalization", "ShinyApp_TunePeaks.rds")

  # Helper function for null-coalescing

  `%||%` = function(a, b) if (is.null(a) || is.na(a)) b else a

  # Creating the persistent data structure to store peak information

  firstParameterData = readRDS(file.path("rds", list.files("rds", pattern = "*.rds")[1]))
  totalSampleNames = firstParameterData@phenoData@data$name
  parametersList = list.files(file.path("output", "2_Normalization", "peaks"), pattern = "*.txt")

  totalPeaksValues = vector(mode = "list", length = length(totalSampleNames))
  names(totalPeaksValues) = totalSampleNames

  foreach::foreach(sample = 1:length(totalSampleNames)) %do%
    {
      currentSample = totalSampleNames[sample]

      currentSample_parametersList = vector(mode = "list", length = length(parametersList))
      names(currentSample_parametersList) = parametersList

      foreach::foreach(p = 1:length(parametersList)) %do%
        {
          currentParameterID = grep(paste("Parameter-", p, "_", sep = ""), parametersList)
          currentParameterName = parametersList[currentParameterID]

          peaksData = utils::read.csv(file = file.path("output", "2_Normalization", "peaks", parametersList[currentParameterID]), sep = "\t")
          associatedPeaksInfos = peaksData[peaksData$Samples == currentSample, ]
          associatedPeaksInfos = as.numeric(associatedPeaksInfos[1, -1])

          peaksList = vector(mode = "list", length = 6)
          foreach::foreach(peaks = 1:3) %do%
            {
              if(peaks <= length(associatedPeaksInfos))
              {
                peaksList[[peaks]] = associatedPeaksInfos[peaks]
                peaksList[[peaks + 3]] = FALSE
              } else
              {
                peaksList[[peaks]] = 0
                peaksList[[peaks + 3]] = TRUE
              }

              names(peaksList)[peaks] = paste("x", peaks, sep = "")
              names(peaksList)[peaks + 3] = paste("c", peaks, sep = "")
            }

          peaksList = append(list(num_peaks = length(associatedPeaksInfos)), peaksList)

          currentSample_parametersList[[currentParameterName]] = peaksList
        }
      totalPeaksValues[[currentSample]] = currentSample_parametersList
    }


  namedListParameters = as.list(names(totalPeaksValues[[1]]))
  names(namedListParameters) = gsub("^peaksAnalysis_Parameter-([0-9]+)_(.+)\\.txt$", "Parameter-\\1 (\\2)", names(totalPeaksValues[[1]]))

  associatedIDs = as.numeric(gsub("^peaksAnalysis_Parameter-([0-9]+)_(.+)\\.txt$", "\\1", names(totalPeaksValues[[1]])))

  namedListParameters = namedListParameters[order(associatedIDs)]

  # Save data

  save_positions_data = function(data, file_path) {

    saveRDS(data, file_path)

  }

  # Shiny app: UI section

  ui = shiny::fluidPage(

    shiny::titlePanel("Peaks fine editor"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(

        shiny::h4(shiny::strong("Samples")),

        shiny::fluidRow(
          shiny::column(6, shiny::actionButton("prev_sample", "< Previous", width = '100%')),
          shiny::column(6, shiny::actionButton("next_sample", "Next >", width = '100%'))
        ),
        shiny::br(),

        shiny::selectInput("sample_choice",
                           label = NULL,
                           choices = totalSampleNames,
                           selected = totalSampleNames[1]),

        shiny::hr(),

        shiny::h4(shiny::strong("Parameters")),

        shiny::fluidRow(
          shiny::column(6, shiny::actionButton("prev_parameter", "< Previous", width = '100%')),
          shiny::column(6, shiny::actionButton("next_parameter", "Next >", width = '100%'))
        ),

        shiny::br(),

        shiny::selectInput("parameter_choice",
                           label = NULL,
                           choices = namedListParameters,
                           selected = namedListParameters[1]),

        shiny::hr(),

        shiny::h4(shiny::strong("Number of peaks")),
        shiny::verbatimTextOutput("num_peaks_display"),

        shiny::hr(),

        shiny::h4(shiny::strong("Peaks tuning")),

        shiny::uiOutput("peak_controls"),

        shiny::h4(shiny::strong("Current position of peaks")),
        shiny::verbatimTextOutput("positions_x"),
      ),

      shiny::mainPanel(
        shiny::plotOutput("density_plot", width = "100%")
      )
    )
  )


  # Shiny app: server section

  server = function(input, output, session) {

    # -----------------------------------------------------------
    # 1. INITIALIZING PERSISTENT DATA
    # -----------------------------------------------------------

    # We create a reactive app_data object based on totalPeaksValues

    app_data = shiny::reactiveValues()
    for (name in totalSampleNames) {
      app_data[[name]] = totalPeaksValues[[name]]
    }

    # Initially displayed sample and parameter are the first available ones

    initial_sample_name = totalSampleNames[1]
    initial_parameter_name = names(totalPeaksValues[[1]])[1]

    # We grab the associated peaks and checkboxes values from the totalPeaksValues

    initial_x1 = totalPeaksValues[[initial_sample_name]][[initial_parameter_name]]$x1
    initial_x2 = totalPeaksValues[[initial_sample_name]][[initial_parameter_name]]$x2
    initial_x3 = totalPeaksValues[[initial_sample_name]][[initial_parameter_name]]$x3

    initial_c1 = totalPeaksValues[[initial_sample_name]][[initial_parameter_name]]$c1
    initial_c2 = totalPeaksValues[[initial_sample_name]][[initial_parameter_name]]$c2
    initial_c3 = totalPeaksValues[[initial_sample_name]][[initial_parameter_name]]$c3

    # QC of initial peaks values (set them to their actual values if they are ok, otherwise set them to 0)

    init_x1_val = if(length(initial_x1) == 1 && !is.na(initial_x1)) initial_x1 else 0
    init_x2_val = if(length(initial_x2) == 1 && !is.na(initial_x2)) initial_x2 else 0
    init_x3_val = if(length(initial_x3) == 1 && !is.na(initial_x3)) initial_x3 else 0

    # Updating peaks values displaying if the checkboxes are checked or not

    current_x = shiny::reactiveValues(
      x1 = if(initial_c1) NA_real_ else init_x1_val,
      x2 = if(initial_c2) NA_real_ else init_x2_val,
      x3 = if(initial_c3) NA_real_ else init_x3_val
    )

    # QC of initial checkboxes values (set them to their actual values if they are ok, otherwise set them to FALSE)

    current_c = shiny::reactiveValues(
      c1 = if(length(initial_c1) == 1 && !is.na(initial_c1)) initial_c1 else FALSE,
      c2 = if(length(initial_c2) == 1 && !is.na(initial_c2)) initial_c2 else FALSE,
      c3 = if(length(initial_c3) == 1 && !is.na(initial_c3)) initial_c3 else FALSE
    )

    # Store the last value from slider to put it back afterwards

    last_valid_x = shiny::reactiveValues(
      x1 = init_x1_val,
      x2 = init_x2_val,
      x3 = init_x3_val
    )

    # Reactive variables which store the inital number of peaks

    initial_num_peaks = totalPeaksValues[[initial_sample_name]][[initial_parameter_name]]$num_peaks
    initial_num_peaks_val = max(1, min(3, as.integer(initial_num_peaks %||% 1)))

    num_peaks_r = shiny::reactiveVal(initial_num_peaks_val)
    prev_num_peaks = shiny::reactiveVal(initial_num_peaks_val)

    # Reactive variables to store the min/max limits of the sliders
    slider_min = shiny::reactiveVal(0)
    slider_max = shiny::reactiveVal(1)

    # Colors and styles of displayed peak lines

    peak_colors = c("red", "blue", "darkgreen")
    peak_lineTypes = c("dashed", "dashed", "dashed")

    # -----------------------------------------------------------
    # 2. AUTOMATIC SAVE LOGIC
    # -----------------------------------------------------------

    shiny::observe({
      dummy_trigger = shiny::reactiveValuesToList(app_data)
      shiny::isolate({
        save_positions_data(dummy_trigger, RDS_filePath)
      })
    })

    # -----------------------------------------------------------
    # 3. LOGIC TO MANAGE THE NUMBER OF PEAKS
    # -----------------------------------------------------------

    # Observe any modification of num_peaks_r variable

    shiny::observeEvent(num_peaks_r(), {
      shiny::req(input$sample_choice, input$parameter_choice)
      new_n = num_peaks_r()
      old_n = prev_num_peaks()

      sample = input$sample_choice
      parameter = input$parameter_choice

      # Save the new number of peaks in app_data

      app_data[[sample]][[parameter]]$num_peaks = new_n

      # Deletion logic

      if (new_n < old_n) {
        for (i in (new_n + 1):old_n) {
          pos_key = paste0("x", i)
          check_key = paste0("c", i)

          app_data[[sample]][[parameter]][[pos_key]] = NULL
          app_data[[sample]][[parameter]][[check_key]] = NULL

          current_x[[pos_key]] = 0
          current_c[[check_key]] = FALSE
          last_valid_x[[pos_key]] = 0
        }
      }

      # Update the previous value for the next call

      prev_num_peaks(new_n)
    }, ignoreInit = TRUE)


    # -----------------------------------------------------------
    # 4. LOADING LOGIC (CHANGING SAMPLE OR PARAMETER)
    # -----------------------------------------------------------

    # Getting the desired peak value

    get_pos = function(key, sample, parameter) {
      val = app_data[[sample]][[parameter]][[key]]
      if (is.null(val)) {
        val = totalPeaksValues[[sample]][[parameter]][[key]]
      }
      if (length(val) != 1 || !is.numeric(val)) {
        return(0)
      }
      return(val)
    }

    # Getting the desired checkbox value

    get_check = function(key, sample, parameter) {
      val = app_data[[sample]][[parameter]][[key]]
      if (is.null(val)) {
        val = totalPeaksValues[[sample]][[parameter]][[key]]
      }
      if (length(val) != 1 || !is.logical(val)) {
        return(FALSE)
      }
      return(val)
    }

    # Observer which triggers when the sample or the peak changes

    shiny::observeEvent({
      input$sample_choice
      input$parameter_choice
    }, {
      shiny::req(input$sample_choice, input$parameter_choice)

      sample = input$sample_choice
      parameter = input$parameter_choice

      # Computing density to get the limits

      currentDensity = current_density_function()
      new_min_x = min(currentDensity$x)
      new_max_x = max(currentDensity$x)

      if (new_min_x == new_max_x) {
        new_min_x = new_min_x - 1
        new_max_x = new_max_x + 1
      }

      new_min_x = floor(new_min_x)
      new_max_x = ceiling(new_max_x)

      # Update the reactive variables for limits

      slider_min(new_min_x)
      slider_max(new_max_x)

      # Get the number of peaks

      current_n = app_data[[sample]][[parameter]]$num_peaks %||% 1
      current_n = max(1, min(3, as.integer(current_n)))

      # Update the reactive variable

      num_peaks_r(current_n)

      # Update prev_num_peaks

      shiny::isolate({
        prev_num_peaks(current_n)
      })

      # Get peaks values and checkboxes states (+ some QC)

      current_c$c1 = get_check("c1", sample, parameter)
      current_c$c2 = get_check("c2", sample, parameter)
      current_c$c3 = get_check("c3", sample, parameter)

      x1_val = get_pos("x1", sample, parameter)
      x2_val = get_pos("x2", sample, parameter)
      x3_val = get_pos("x3", sample, parameter)

      last_valid_x$x1 = if(is.na(x1_val)) last_valid_x$x1 else x1_val
      last_valid_x$x2 = if(is.na(x2_val)) last_valid_x$x2 else x2_val
      last_valid_x$x3 = if(is.na(x3_val)) last_valid_x$x3 else x3_val

      current_x$x1 = if(current_c$c1) NA_real_ else x1_val
      current_x$x2 = if(current_c$c2) NA_real_ else x2_val
      current_x$x3 = if(current_c$c3) NA_real_ else x3_val

      # Update limits and input controls

      n = current_n

      shiny::isolate({
        for (i in 1:3) {
          pos_key = paste0("x", i)
          check_key = paste0("c", i)

          slider_val = last_valid_x[[pos_key]]

          if (slider_val < new_min_x || slider_val > new_max_x) {
            slider_val = (new_min_x + new_max_x) / 2
          }

          shiny::updateSliderInput(session,
                                   paste0(pos_key, "_slider"),
                                   value = slider_val,
                                   min = new_min_x,
                                   max = new_max_x)

          shiny::updateCheckboxInput(session,
                                     paste0(check_key, "_checkbox"),
                                     value = current_c[[check_key]])
        }
      })
    }, ignoreInit = FALSE)

    # Observers for previous/next buttons

    parameter_choices_vec = names(totalPeaksValues[[1]])

    shiny::observeEvent(input$next_parameter, {
      current_index = match(input$parameter_choice, parameter_choices_vec)
      new_index = if (current_index == length(parameter_choices_vec)) 1 else current_index + 1
      shiny::updateSelectInput(session, "parameter_choice", selected = parameter_choices_vec[new_index])
    })

    shiny::observeEvent(input$prev_parameter, {
      current_index = match(input$parameter_choice, parameter_choices_vec)
      new_index = if (current_index == 1) length(parameter_choices_vec) else current_index - 1
      shiny::updateSelectInput(session, "parameter_choice", selected = parameter_choices_vec[new_index])
    })

    sample_choices_vec = totalSampleNames

    shiny::observeEvent(input$next_sample, {
      current_index = match(input$sample_choice, sample_choices_vec)
      new_index = if (current_index == length(sample_choices_vec)) 1 else current_index + 1
      shiny::updateSelectInput(session, "sample_choice", selected = sample_choices_vec[new_index])
    })

    shiny::observeEvent(input$prev_sample, {
      current_index = match(input$sample_choice, sample_choices_vec)
      new_index = if (current_index == 1) length(sample_choices_vec) else current_index - 1
      shiny::updateSelectInput(session, "sample_choice", selected = sample_choices_vec[new_index])
    })

    # -----------------------------------------------------------
    # 5. DYNAMIC RENDERING OF UI CONTROLS
    # -----------------------------------------------------------

    # Rendering the non-editable number of peaks value

    output$num_peaks_display = shiny::renderText({
      return(num_peaks_r())
    })

    # Rendering the dynamic UI

    output$peak_controls = shiny::renderUI({

      shiny::req(num_peaks_r())
      n = num_peaks_r()

      min_x = slider_min()
      max_x = slider_max()

      control_list = lapply(1:n, function(i) {
        pos_key = paste0("x", i)
        check_key = paste0("c", i)

        initial_check_value = FALSE
        initial_slider_value = (min_x + max_x) / 2

        shiny::tagList(
          shiny::h4(shiny::span(paste("Peak ", i), style = paste0("color: ", peak_colors[i], ";"))),

          shiny::checkboxInput(paste0(check_key, "_checkbox"),
                               label = "Unset (set to NA)",
                               value = initial_check_value),

          shiny::sliderInput(paste0(pos_key, "_slider"),
                             paste0("Peak x position:"),
                             min = min_x, max = max_x,
                             value = initial_slider_value, step = 0.01),
          shiny::hr()
        )
      })

      do.call(shiny::tagList, control_list)
    })

    # Observer to deactivate the slider if necessary

    shiny::observe({
      shiny::req(num_peaks_r())
      for (i in 1:num_peaks_r()) {
        check_key = paste0("c", i)
        pos_key = paste0("x", i, "_slider")

        is_disabled = shiny::isolate(current_c[[check_key]])

        session$sendInputMessage(pos_key, list(disabled = is_disabled))
      }
    })

    # -----------------------------------------------------------
    # 6. SYNCHRONISATION : WRITING IN THE PERSISTENT DATA STRUCTURE
    # -----------------------------------------------------------

    update_peak_data = function(index) {

      shiny::req(num_peaks_r() >= index)
      pos_key = paste0("x", index)
      check_key = paste0("c", index)
      slider_key = paste0(pos_key, "_slider")
      check_box_key = paste0(check_key, "_checkbox")

      sample = input$sample_choice
      parameter = input$parameter_choice

      is_checked = input[[check_box_key]]

      # Update of the checkbox state

      app_data[[sample]][[parameter]][[check_key]] = is_checked
      current_c[[check_key]] = is_checked

      shiny::isolate({
        # Update the peak x value

        if (is_checked) {
          app_data[[sample]][[parameter]][[pos_key]] = NA_real_
          current_x[[pos_key]] = NA_real_
        } else {
          val = shiny::isolate(input[[slider_key]])

          app_data[[sample]][[parameter]][[pos_key]] = val
          current_x[[pos_key]] = val
          last_valid_x[[pos_key]] = val
        }
      })
    }

    # Observer for each slider and checkbox

    shiny::observeEvent(input$x1_slider, update_peak_data(1), ignoreInit = TRUE)
    shiny::observeEvent(input$c1_checkbox, update_peak_data(1), ignoreInit = TRUE)
    shiny::observeEvent(input$x2_slider, update_peak_data(2), ignoreInit = TRUE)
    shiny::observeEvent(input$c2_checkbox, update_peak_data(2), ignoreInit = TRUE)
    shiny::observeEvent(input$x3_slider, update_peak_data(3), ignoreInit = TRUE)
    shiny::observeEvent(input$c3_checkbox, update_peak_data(3), ignoreInit = TRUE)

    # -----------------------------------------------------------
    # 7. DENSITY FUNCITONS AND RENDERING
    # -----------------------------------------------------------

    current_density_function = shiny::reactive({

      sample_name = input$sample_choice
      parameter_name = input$parameter_choice

      shiny::req(sample_name, parameter_name)

      parameterName = gsub("^(.+)_(.+)_(.+)$", "\\3", parameter_name)
      parameterName = gsub("\\.txt", "", parameterName)

      RDSFilesList = list.files("rds", pattern = "*.rds")
      RDSFileToOpenID = grep(paste("(.+)_", parameterName, "_(.+)$", sep = ""), RDSFilesList)

      if(length(RDSFileToOpenID) == 0) {
        return(list(x = c(-1, 1), y = c(0, 0)))
      }

      RDSFile = readRDS(file.path("rds", RDSFilesList[RDSFileToOpenID]))

      # Use of isolated local values

      patientData = RDSFile[RDSFile@phenoData@data$name == sample_name, ]

      if(length(patientData) == 0 || nrow(patientData[[1]]@exprs) == 0) {
        return(list(x = c(-1, 1), y = c(0, 0)))
      }

      data_for_density = as.numeric(patientData[[1]]@exprs[, 1])

      if(length(data_for_density) < 2) {
        return(list(x = c(min(data_for_density)-1, max(data_for_density)+1), y = c(0, 0)))
      }

      density = density(data_for_density)

      outputValues = list(x = density$x, y = density$y)
      return(outputValues)
    })

    # Displaying peaks x positions

    output$positions_x = shiny::renderText({
      shiny::req(num_peaks_r())
      n = num_peaks_r()

      text_output = ""
      for (i in 1:n) {
        pos_key = paste0("x", i)

        val = current_x[[pos_key]]
        display_val = if (is.na(val)) "NA (unset)" else format(val, digits = 4)

        text_output = paste0(text_output,
                             paste0("Peak ", i, " (", peak_colors[i], "): x = ", display_val, "\n"))
      }
      return(text_output)
    })

    # Rendering density graph

    output$density_plot = shiny::renderPlot({

      shiny::req(num_peaks_r(), input$sample_choice, input$parameter_choice)

      currentDensity = current_density_function()

      data_density = data.frame(
        x = currentDensity$x,
        y = currentDensity$y
      )

      plot_title = NULL

      x_min = floor(min(data_density$x))
      x_max = ceiling(max(data_density$x))
      y_max = ceiling(max(data_density$y) * 1.05)

      p = ggplot2::ggplot(data_density, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line(color = "blue", linewidth = 1) +
        ggplot2::geom_area(fill = "lightblue", alpha = 0.5) +
        ggplot2::labs(title = plot_title,
             x = "Fluorescence intensity (x)", y = "Density (y)") +
        ggplot2::theme_minimal() +
        ggplot2::coord_cartesian(xlim = c(x_min, x_max),
                        ylim = c(0, max(y_max, 0.1)))

      n = num_peaks_r()

      for (i in 1:n) {
        pos_key = paste0("x", i)

        if (!is.na(current_x[[pos_key]])) {
          p = p + ggplot2::geom_vline(
            xintercept = current_x[[pos_key]],
            color = peak_colors[i],
            linetype = peak_lineTypes[i],
            size = 1.2
          )
        }
      }

      print(p)


    },

    # Defining the height of the graph as its current width

    height = function() {
      session$clientData$output_density_plot_width
    })
  }

  # Launch the Shiny app

  shiny::shinyApp(ui = ui, server = server)

}

#' Reexport peak tabular text files from RDS peaks data structure
#'
#' This function allows to reexport all the tabular text files containing all the corrected peaks for every parameter and sample after peaks have been corrected through the dedicatedShiny application.
#'
#' @importFrom foreach %do%
#'
#' @export

exportRDSPeaksDataToTabularTextFiles = function()
{
  sample = NULL
  parameter = NULL
  peak = NULL

  RDSDataToExport = readRDS(file = file.path("output", "2_Normalization", "ShinyApp_TunePeaks.rds"))

  formattedData = data.frame()

  foreach::foreach(sample = 1:length(RDSDataToExport)) %do%
    {
      currentSampleData = RDSDataToExport[[sample]]
      currentSampleName = names(RDSDataToExport)[sample]

      foreach::foreach(parameter = 1:length(currentSampleData)) %do%
        {
          currentParameterData = currentSampleData[[parameter]]
          currentParameterName = names(currentSampleData)[parameter]

          nbOfPeaks = as.numeric(currentParameterData$num_peaks)

          formattedData = rbind(formattedData, c(currentSampleName, currentParameterName, nbOfPeaks, currentParameterData[[2]], currentParameterData[[3]], currentParameterData[[4]]))

        }

    }

  colnames(formattedData) = c("Samples", "Parameter", "NbOfPeaks", "Peak1", "Peak2", "Peak3")

  uniqueParameters = unique(formattedData$Parameter)
  uniqueSamples = unique(formattedData$Sample)

  foreach::foreach(parameter = 1:length(uniqueParameters)) %do%
    {
      currentParameter = uniqueParameters[parameter]
      currentParameterData = data.frame()

      foreach::foreach(sample = 1:length(uniqueSamples)) %do%
        {
          currentSample = uniqueSamples[sample]

          currentData = formattedData[formattedData$Sample == currentSample & formattedData$Parameter == currentParameter, ]

          currentParameterData = rbind(currentParameterData, currentData)

        }

      currentParameter_peaksNb = as.numeric(unique(currentParameterData$NbOfPeaks))

      currentParameter_peaksColumnsToKeep = NULL

      foreach::foreach(peak = 1:currentParameter_peaksNb) %do%
        {
          currentParameter_peaksColumnsToKeep = c(currentParameter_peaksColumnsToKeep, paste("Peak", peak, sep = ""))

        }

      dataFrameToModify = as.data.frame(currentParameterData[, currentParameter_peaksColumnsToKeep])
      colnames(dataFrameToModify) = currentParameter_peaksColumnsToKeep

      currentParameterData[, currentParameter_peaksColumnsToKeep] = apply(dataFrameToModify, 2, as.numeric)

      currentParameterData = currentParameterData[, colnames(currentParameterData) %in% c("Samples", currentParameter_peaksColumnsToKeep)]

      utils::write.table(currentParameterData, file = file.path("output", "2_Normalization", "peaks", currentParameter), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    }
}
