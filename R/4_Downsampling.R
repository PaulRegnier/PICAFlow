#' Reduce and downsample dataset
#'
#' This function allows to reduce the dataset by removing useless paramaters as well as tag cells for the construction of a downsampled dataset in which all groups contribute equally to the dataset and in which each sample contributes equally within a given group. Please note that from this step, `rds` files will contain table-like data and not `flowFrame` nor `flowSet` objects anymore.
#'
#' @param flowSet A string defining the `flowSet` to use. Defaults to `NULL`.
#'
#' @param groupVector A character vector defining the groups for each sample. Defaults to `NULL`.
#'
#' @param parametersToKeep A character vector defining the parameters to keep in the subsequent dataset. Typically, one wants to remove parameters that were already used to gate cells. Defaults to `NULL`.
#'
#' @param downsampleMinEvents An integer defining the minimum number of cells per sample to be included in the downsampled dataset. Defaults to `10000`.
#'
#' @param rescale A boolean defining if the data should be rescaled parameter-by-parameter. Defaults to `TRUE`.
#'
#' @param rescale_min An integer defining the minimum value for the rescaled dataset. Used only if `rescale` is set to `TRUE`. Defaults to `1`.
#'
#' @param rescale_max An integer defining the maximum value for the rescaled dataset. Used only if `rescale` is set to `TRUE`. Defaults to `10`.
#'
#' @param maxCellsNb An integer defining the maximum number of cells for the downsampled dataset. If the actual number of cells is greater than `maxCellsNb`, a global downsampling factor will be applied to keep the number of cells of the downsampled dataset below `maxCellsNb`. Defaults to `500000`.
#'
#' @param estimateThreshold A boolean defining if a given `maxCellsNb` threshold should rather be replaced by the test of all thresholds possible according to the number of cells in each sample (if set to `TRUE`). Be careful, this function is expensive in terms of computational resources, so do not use it if your dataset is rather large. Defaults to `FALSE`.
#'
#' @param coresNumber An integer defining the number of cores to use to analyze the peaks. Be careful, this function is expensive in terms of computational resources, so do not increase this value too much (typically 1 to 4 maximum). Defaults to `1`.
#'
#' @return If `estimateThreshold = FALSE`, the function will return a list of 2 elements named `data` (containing the actual downsampled data) and `message` (containing the log messages generated during the downsampling process). If `estimateThreshold = TRUE`, the function will rather return a plot in PDF format showing the relation between a given threshold and the final number of cells after downsampling.
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#'
#' @export

poolData = function(flowSet = NULL, groupVector = NULL, parametersToKeep = NULL, downsampleMinEvents = 10000, rescale = TRUE, rescale_min = 1, rescale_max = 10, maxCellsNb = 5e+05, estimateThreshold = FALSE, coresNumber = 1)
{
  g = NULL
  a = NULL
  z = NULL
  w = NULL
  u = NULL
  v = NULL

  outputMessage = NULL
  samplesToRemoveNames = NULL

  uniqueGroups = unique(groupVector)
  totalDataList = rep(list(NULL), length(uniqueGroups))
  discardedSamplesWhenSampling = NULL
  

  totalSamplesSizes = NULL
  totalSamplesNames = NULL
  totalSamplesGroups = NULL

  pb = tcltk::tkProgressBar("Busy...", paste("Group 0/", length(uniqueGroups), sep = ""), 0, length(uniqueGroups), 200)

  foreach::foreach(g = 1:length(uniqueGroups)) %do%
  {
    currentGroupID = uniqueGroups[g]

    currentGroupSamplesID = which(groupVector == currentGroupID)

    currentGroupDataList = rep(list(NULL), length(currentGroupSamplesID))

    foreach::foreach(a = 1:length(currentGroupSamplesID)) %do%
    {
      currentSampleID = currentGroupSamplesID[a]
      currentFileData = flowSet[[currentSampleID]]
      currentSampleName = gsub("(.+)_(.+)_(.+)", "\\3", flowSet@phenoData@data$name[currentSampleID])

      res = data.frame(currentFileData@exprs, stringsAsFactors = FALSE)

		if(nrow(res) > 0)
		{
		  res$group = currentGroupID
      res$name = currentSampleName

      currentGroupDataList[[a]] = res
		} else
		{
		currentGroupDataList[[a]] = NULL

		}


      totalSamplesSizes = c(totalSamplesSizes, nrow(res))

      totalSamplesNames = c(totalSamplesNames, currentSampleName)

      totalSamplesGroups = c(totalSamplesGroups, currentGroupID)
     
    }
	
	samplesToRemoveID = which(totalSamplesSizes == 0)
	samplesToRemoveNames = c(samplesToRemoveNames, totalSamplesNames[samplesToRemoveID])
	
	if(length(samplesToRemoveID) > 0)
	{
	
	
	
	currentGroupDataList = currentGroupDataList[-c(samplesToRemoveID)]
	totalSamplesSizes = totalSamplesSizes[-c(samplesToRemoveID)]

      totalSamplesNames = totalSamplesNames[-c(samplesToRemoveID)]

      totalSamplesGroups = totalSamplesGroups[-c(samplesToRemoveID)]
	}

    gc()

    currentGroupData = do.call(rbind.data.frame, currentGroupDataList)

    totalDataList[[g]] = currentGroupData
    names(totalDataList)[g] = currentGroupID
    gc()

    tcltk::setTkProgressBar(pb, g, label = paste("Group ", g, "/", length(uniqueGroups), sep = ""))
  }

  close(pb)
  
  
  if(length(samplesToRemoveNames) == 1)
	{
		outputMessage = c(outputMessage, paste("## Prior to any downsampling, the following sample was discarded because there were no cells remaining in the associated flowSet: ", samplesToRemoveNames, " ##", sep = ""))

	}
	
	if(length(samplesToRemoveNames) > 1)
	{
		outputMessage = c(outputMessage, paste("## Prior to any downsampling, the following samples were discarded because there were no cells remaining in the associated flowSets: ", paste(samplesToRemoveNames, collapse = ", "), " ##", sep = ""))

	}
  

  if (estimateThreshold == TRUE)
  {
    thresholdsToTest = sort(totalSamplesSizes)
  } else
  {
    thresholdsToTest = downsampleMinEvents
  }

  cl = parallel::makeCluster(coresNumber, type = "PSOCK")
  doSNOW::registerDoSNOW(cl)

  pb = utils::txtProgressBar(min = 0, max = length(thresholdsToTest), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)

  estimateThresholdTable = foreach::foreach(z = thresholdsToTest, .packages = c("foreach", "tcltk"), .combine = "c", .options.snow = opts) %dopar%
  {
    editedTotalDataList = totalDataList

    totalDeletedDataList = list()

    downsampleMinEvents = z

    samplesToDiscardID = which(totalSamplesSizes < downsampleMinEvents)
    samplesToDiscardNames = totalSamplesNames[samplesToDiscardID]
    samplesToDiscardGroups = totalSamplesGroups[samplesToDiscardID]

    if (length(samplesToDiscardNames) > 0)
    {
      foreach::foreach(w = 1:length(samplesToDiscardID)) %do%
      {
        currentSampleName = samplesToDiscardNames[w]
        currentSampleGroup = samplesToDiscardGroups[w]

        currentSampleNameForList = paste(currentSampleGroup, "_", currentSampleName, sep = "")

        currentSampleGroupData = editedTotalDataList[[currentSampleGroup]]
        currentSampleGroupDataDeleted = currentSampleGroupData[currentSampleGroupData$name == currentSampleName, ]
        currentSampleGroupData = currentSampleGroupData[currentSampleGroupData$name != currentSampleName, ]

        editedTotalDataList[[currentSampleGroup]] = currentSampleGroupData
        totalDeletedDataList[[currentSampleNameForList]] = currentSampleGroupDataDeleted
      }
    }

    totalSamplesMinSizes = vector(mode = "list", length = length(editedTotalDataList))
    totalGroupsMinSizes = vector(mode = "list", length = length(editedTotalDataList))
    foreach::foreach(x = 1:length(editedTotalDataList)) %do%
    {
      currentGroupData = editedTotalDataList[[x]]

      currentGroupSamplesUnique = unique(currentGroupData$name)

      currentGroupLengths = vector(mode = "list", length = length(currentGroupSamplesUnique))
      foreach::foreach(y = 1:length(currentGroupSamplesUnique)) %do%
      {
        currentSampleName = currentGroupSamplesUnique[y]

        currentSampleAssociatedLength = nrow(currentGroupData[currentGroupData$name == currentSampleName, ])
        currentGroupLengths[[y]] = currentSampleAssociatedLength
      }

      totalSamplesMinSizes[[x]] = min(unlist(currentGroupLengths))
      totalGroupsMinSizes[[x]] = min(unlist(currentGroupLengths)) * length(currentGroupSamplesUnique)
    }

    totalSamplesMinSizes = as.vector(unlist(totalSamplesMinSizes))
    totalGroupsMinSizes = as.vector(unlist(totalGroupsMinSizes))

    if (all(totalGroupsMinSizes > 0) == TRUE)
    {
      if (estimateThreshold == FALSE)
      {
        outputMessage = c(outputMessage, paste("## Downsampling activated (minimum events number = ", downsampleMinEvents, ") ##", sep = ""))

        if (length(samplesToDiscardNames) == 0)
        {
          outputMessage = c(outputMessage, paste("## No samples were discarded ##", sep = ""))
        } else if (length(samplesToDiscardNames) == 1)
        {
          outputMessage = c(outputMessage, paste("## ", length(samplesToDiscardNames), " sample was discarded: ", samplesToDiscardNames, " ##", sep = ""))
        } else
        {
          outputMessage = c(outputMessage, paste("## ", length(samplesToDiscardNames), " samples were discarded: ", paste(samplesToDiscardNames, collapse = ", "), " ##", sep = ""))
        }
      }

      downsamplingFactors = 1/(totalGroupsMinSizes/min(totalGroupsMinSizes))

      finalGroupsRealNumbers = 0
      downsampledDataset = list()

      foreach::foreach(u = 1:length(editedTotalDataList)) %do%
      {
        currentGroupData = editedTotalDataList[[u]]
        currentGroupDownsamplingFactor = downsamplingFactors[[u]]
        currentGroupMinSampleSize = totalSamplesMinSizes[[u]]
        currentGroupName = unique(currentGroupData$group)

        currentGroupLength = length(unique(currentGroupData$name))

        currentGroupFinalCellsNumber = round(currentGroupMinSampleSize * currentGroupDownsamplingFactor, 0) * currentGroupLength

        if ((currentGroupFinalCellsNumber * length(editedTotalDataList)) > maxCellsNb)
        {
          if (estimateThreshold == FALSE)
          {
            outputMessage = c(outputMessage, paste("## Group ", u, " (", currentGroupName, "): Warning! the total number of cells exceeds the defined maximum threshold of ", as.numeric(maxCellsNb), " cells. The associated downsampling factor of ", currentGroupDownsamplingFactor, " was corrected to ", (currentGroupDownsamplingFactor * (maxCellsNb/length(editedTotalDataList)))/currentGroupFinalCellsNumber, " ##", sep = ""))

          }

          currentGroupDownsamplingFactor = (currentGroupDownsamplingFactor * (maxCellsNb/length(editedTotalDataList)))/currentGroupFinalCellsNumber
        }

        foreach::foreach(v = 1:currentGroupLength) %do%
        {
          currentSampleName = unique(currentGroupData$name)[[v]]
          currentSampleData = currentGroupData[currentGroupData$name == currentSampleName, ]

          currentSampleNameForList = paste(currentGroupName, "_", currentSampleName, sep = "")

          samplePool = 1:nrow(currentSampleData)
          linesToKeep = sample(samplePool, round(currentGroupMinSampleSize * currentGroupDownsamplingFactor, 0))

          linesToDeleteID = which(samplePool %in% linesToKeep == FALSE)
          linesToDelete = samplePool[linesToDeleteID]

          downsampledDataset[[currentSampleNameForList]] = currentSampleData[linesToKeep, ]

          totalDeletedDataList[[currentSampleNameForList]] = currentSampleData[linesToDelete, ]
        }

        if (estimateThreshold == FALSE)
        {
          outputMessage = c(outputMessage, paste("## Group ", u, " (", currentGroupName, "): ", round(currentGroupMinSampleSize * currentGroupDownsamplingFactor, 0), " cells were taken in each of the ", currentGroupLength, " samples, totalling ", currentGroupFinalCellsNumber, " cells for this group (downsampling factor = ", currentGroupDownsamplingFactor, ") ##", sep = ""))
        }

        finalGroupsRealNumbers = finalGroupsRealNumbers + currentGroupFinalCellsNumber
      }

      editedTotalDataList = downsampledDataset
      deletedDataNotSampledList = totalDeletedDataList

      currentThresholdResult = list(c(downsampleMinEvents, finalGroupsRealNumbers), editedTotalDataList, deletedDataNotSampledList)
    } else
    {
      currentThresholdResult = list(c(downsampleMinEvents, NULL), NULL, NULL)
    }

    return(list(data = currentThresholdResult, message = outputMessage))
  }

  close(pb)

  parallel::stopCluster(cl)

  outputMessage = estimateThresholdTable$message
  estimateThresholdTable = estimateThresholdTable$data

  editedTotalDataList_ID = seq(1, length(estimateThresholdTable), 3) + 1
  editedTotalDataList = estimateThresholdTable[[editedTotalDataList_ID]]

  deletedDataNotSampledList_ID = seq(1, length(estimateThresholdTable), 3) + 2
  deletedDataNotSampledList = estimateThresholdTable[[deletedDataNotSampledList_ID]]

  estimateThresholdTable_ID = seq(1, length(estimateThresholdTable), 3)
  estimateThresholdTable = do.call(rbind.data.frame, estimateThresholdTable[estimateThresholdTable_ID])

  colnames(estimateThresholdTable) = c("Threshold", "TotalCellsNumbersRetained")

  if (estimateThreshold == FALSE)
  {
    totalDataSampled = do.call(rbind.data.frame, editedTotalDataList)
    colnames(totalDataSampled) = as.vector(colnames(editedTotalDataList[[1]]))
    rownames(totalDataSampled) = NULL
    totalDataSampled$state = "sampled"

    totalDataNotSampled = do.call(rbind.data.frame, deletedDataNotSampledList)
    colnames(totalDataNotSampled) = as.vector(colnames(deletedDataNotSampledList[[1]]))
    rownames(totalDataNotSampled) = NULL

    totalDataNotSampled$state = "notSampled"

    totalData = rbind(totalDataSampled, totalDataNotSampled)

    totalData = data.frame(totalData, stringsAsFactors = FALSE)

    groups = totalData$group
    names = totalData$name
    states = totalData$state

    colnames(totalData) = c(flowSet[[1]]@parameters@data$desc, "group", "name", "state")

    totalData = totalData[, -c((ncol(totalData) - 2):ncol(totalData))]
    totalData = totalData[, order(colnames(totalData))]

    totalData$group = groups
    totalData$name = names
    totalData$state = states

    totalData = totalData[sample(nrow(totalData)), ]

    rownames(totalData) = NULL

    foreach::foreach(c = 1:(ncol(totalData) - 3)) %do%
    {
      if (rescale == TRUE)
      {
        currentData = totalData[, c]

        max_old = max(currentData)
        min_old = min(currentData)

        currentData = (((rescale_max - rescale_min)/(max_old - min_old)) * (currentData - max_old)) + rescale_max

        totalData[, c] = currentData
      }
    }

    if (is.null(parametersToKeep) == FALSE)
    {
      totalData = totalData[, c(parametersToKeep, "group", "name", "state")]
    }
  }

  gc()

  if (estimateThreshold == TRUE)
  {
    y = estimateThresholdTable$TotalCellsNumbersRetained
    x = estimateThresholdTable$Threshold
    mod = stats::loess(y ~ x, span = 0.25)
    yfit = stats::predict(mod, newdata = x)

    grDevices::pdf(file.path("output", "3_Gating", "CutThreshold-vs-FinalDatasetCellNumber.pdf"))

    plot(x, y, xlab = "Cut threshold", ylab = "Final dataset cell number")
    graphics::grid(NULL, NULL)
    graphics::lines(x, yfit, col = "blue", lty = 3, lwd = 2)

    grDevices::dev.off()

  } else
  {
    return(list(data = totalData, message = outputMessage))
  }
}

#' Export all outputs from downsampling
#'
#' This function allows to export data and log from the downsampling step.
#'
#' @param poolData A list of length 2 produced by `poolData()` method. The `data` element will be exported as a new `rds` file, whereas the `message` element will be exported as a text file. Defaults to `NULL`.
#'
#' @return Generated `rds` file is saved to `rds` directory. Generated text file is saved to `output > 4_Downsampling` directory.
#'
#' @export

exportDownsamplingOutput = function(poolData = NULL)
{
  saveRDS(poolData$data, file.path("rds", "4_Downsampled_full.rds"))

  downsamplingLog = poolData$message
  outputLog = file(file.path("output", "4_Downsampling", "downsamplingLog.txt"))
  writeLines(downsamplingLog, outputLog)
  close(outputLog)
}


#' Generate two subdatasets from downsampled dataset
#'
#' This function allows to export two different subdatasets from the original downsampled dataset. This is typically used when one wants to generate training and validation datasets. Please note that this function do not directly split the samples between the two subdatasets, but rather the cells themselves. It means that at the end, each subdataset will have the same number of samples as compared to the original dataset, but will not contain the same cells: some will be attributed to the training subdataset, whereas the others will be attributed to the validation subdataset.
#'
#' @param trainingDatasetProportion An decimal numeric value (between 0 and 1) defining the frequency of cells that will go to the training subdataset. Defaults to `0.75`.
#'
#' @return Generated `rds` files are saved to `rds` directory.
#'
#' @export

splitDataset = function(trainingDatasetProportion = 0.75)
{
  data = readRDS(file.path("rds", "4_Downsampled_full.rds"))

  validationDatasetProportion = 1 - trainingDatasetProportion

  sampledData = data[data$state == "sampled", ]

  notSampledData = data[data$state == "notSampled", ]

  sampledData_trainingRowsToKeepID = sample(1:nrow(sampledData), round(nrow(sampledData) * trainingDatasetProportion))
  notSampledData_trainingRowsToKeepID = sample(1:nrow(notSampledData), round(nrow(notSampledData) * trainingDatasetProportion))

  sampledData_validationRowsToKeepID = which((1:nrow(sampledData)) %in% sampledData_trainingRowsToKeepID == FALSE)
  notSampledData_validationRowsToKeepID = which((1:nrow(notSampledData)) %in% notSampledData_trainingRowsToKeepID == FALSE)

  trainingDataset_sampled = sampledData[sampledData_trainingRowsToKeepID, ]
  validationDataset_sampled = sampledData[sampledData_validationRowsToKeepID, ]

  trainingDataset_notSampled = notSampledData[notSampledData_trainingRowsToKeepID, ]
  validationDataset_notSampled = notSampledData[notSampledData_validationRowsToKeepID, ]

  trainingDataset = rbind(trainingDataset_sampled, trainingDataset_notSampled)
  validationDataset = rbind(validationDataset_sampled, validationDataset_notSampled)

  saveRDS(trainingDataset, file.path("rds", "4_Downsampled_training.rds"))
  saveRDS(validationDataset, file.path("rds", "4_Downsampled_validation.rds"))

  rm(trainingDataset_sampled)
  rm(validationDataset_sampled)
  rm(trainingDataset_notSampled)
  rm(validationDataset_notSampled)
  rm(trainingDataset)
  rm(validationDataset)
  gc()
}
