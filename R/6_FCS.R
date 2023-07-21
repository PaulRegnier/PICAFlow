#' Export FCS files from flowSet
#'
#' This function allows to export new FCS files from a data frame.
#'
#' @param data A data frame from which to generate the FCS files. It must originate from `openDownsampledData()` method call and must contain the UMAP parameters for each cell. Defaults to `NULL`.
#'
#' @param datasetFolder A string defining in which folder to save the results. It should match the `datasetToUse` value. Can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. The value must match the origin of the `data` argument. Defaults to `full`.
#'
#' @return Generated FCS and text files are saved to `output > 6_FCS` directory.
#'
#' @importFrom foreach %do%
#'
#' @export

exportFCS = function(data = NULL, datasetFolder = "full")
{
  a = NULL
  d = NULL
  workingDirectory = getwd()

  if (dir.exists(file.path(workingDirectory, "output", "6_FCS", "samples")) == FALSE)
  {
    dir.create(file.path(workingDirectory, "output", "6_FCS", "samples"))
  } else
  {
    unlink(file.path(workingDirectory, "output", "6_FCS", "samples", "*.*"))
  }

  if (dir.exists(file.path(workingDirectory, "output", "6_FCS", "groups")) == FALSE)
  {
    dir.create(file.path(workingDirectory, "output", "6_FCS", "groups"))
  } else
  {
    unlink(file.path(workingDirectory, "output", "6_FCS", "groups", "*.*"))
  }

  if (dir.exists(file.path(workingDirectory, "output", "6_FCS", "samples", datasetFolder)) == FALSE)
  {
    dir.create(file.path(workingDirectory, "output", "6_FCS", "samples", datasetFolder))
  } else
  {
    unlink(file.path(workingDirectory, "output", "6_FCS", "samples", datasetFolder, "*.*"))
  }

  if (dir.exists(file.path(workingDirectory, "output", "6_FCS", "groups", datasetFolder)) == FALSE)
  {
    dir.create(file.path(workingDirectory, "output", "6_FCS", "groups", datasetFolder))
  } else
  {
    unlink(file.path(workingDirectory, "output", "6_FCS", "groups", datasetFolder, "*.*"))
  }

  groupsToPlot = length(unique(data$group))

  totalCorrespondenceSamplesID = data.frame(matrix(NA, nrow = 0, ncol = 5), stringsAsFactors = FALSE)

  pb = tcltk::tkProgressBar("Exporting FCS files...", paste("Group 0/", groupsToPlot, sep = ""), 0, groupsToPlot, 0, 300)

  foreach::foreach(a = 1:groupsToPlot) %do%
  {
    currentGroup = unique(data$group)[a]

    subsetData = data[data$group == currentGroup, ]

    rowsToReplace = which(subsetData$group == currentGroup)
    subsetData[rowsToReplace, "group"] = a

    subsetDataFCS = subsetData

    finalCorrespondenceSamplesID = vector(mode = "list", 0)

    currentGroupSamplesList = unique(subsetDataFCS$name)

    foreach::foreach(q = 1:length(currentGroupSamplesList)) %do%
    {
      currentGroupCurrentSample = currentGroupSamplesList[q]

      rowsToEdit = which(subsetDataFCS$name == currentGroupCurrentSample)

      if (a == 1)
      {
        startPoint = 0
      } else
      {
        startPoint = max(as.numeric(totalCorrespondenceSamplesID[, 5]))
      }

      subsetDataFCS[rowsToEdit, "name"] = startPoint + q
      finalCorrespondenceSamplesID[[q]] = c(a, currentGroup, q, currentGroupCurrentSample, (startPoint + q))
    }

    finalCorrespondenceSamplesID = t(as.data.frame(finalCorrespondenceSamplesID))

    totalCorrespondenceSamplesID = rbind(totalCorrespondenceSamplesID, finalCorrespondenceSamplesID)

    colsToDeleteID = which(colnames(subsetDataFCS) %in% c("group"))
    subsetDataFCS = subsetDataFCS[, -colsToDeleteID]
    colnameToReplaceID = which(colnames(subsetDataFCS) == "name")
    colnames(subsetDataFCS)[colnameToReplaceID] = "ID"

    colnameToReplaceID = which(colnames(subsetDataFCS) == "state")
    colnames(subsetDataFCS)[colnameToReplaceID] = "UMAP_Train"

    subsetDataFCS$UMAP_Train[subsetDataFCS$UMAP_Train != "sampled"] = 0

    subsetDataFCS$UMAP_Train[subsetDataFCS$UMAP_Train == "sampled"] = 1

    if (length(unique(subsetDataFCS$clusteringGroup)) > 0)
    {
      subsetDataFCS$clusteringGroup[subsetDataFCS$clusteringGroup != "training"] = 1

      subsetDataFCS$clusteringGroup[subsetDataFCS$clusteringGroup == "training"] = 0
    }

    subsetDataFCS = as.matrix(apply(subsetDataFCS, 2, as.numeric))

    if (a == 1)
    {
      fullTotalData = subsetDataFCS
    } else
    {
      fullTotalData = rbind(fullTotalData, subsetDataFCS)
    }

    subsetDataFCS = flowCore::flowFrame(subsetDataFCS)

    flowCore::write.FCS(subsetDataFCS, filename = file.path(workingDirectory, "output", "6_FCS", "groups", datasetFolder, paste("Group-", currentGroup, "_AllSamples_", datasetFolder, ".fcs", sep = "")), delimiter = "|")

    uniqueSamples = unique(subsetData$name)

    foreach::foreach(d = 1:length(uniqueSamples)) %do%
    {
      currentSample = uniqueSamples[d]

      subsetDataSample = subsetData

      subsetDataSample = subsetDataSample[subsetDataSample$name == currentSample, ]

      subsetDataSampleFCS = subsetDataSample

      colnameToReplaceID = which(colnames(subsetDataSampleFCS) == "state")
      colnames(subsetDataSampleFCS)[colnameToReplaceID] = "UMAP_Train"

      subsetDataSampleFCS$UMAP_Train[subsetDataSampleFCS$UMAP_Train != "sampled"] = 0

      subsetDataSampleFCS$UMAP_Train[subsetDataSampleFCS$UMAP_Train == "sampled"] = 1

      if (length(unique(subsetDataSampleFCS$clusteringGroup)) > 0)
      {
        subsetDataSampleFCS$clusteringGroup[subsetDataSampleFCS$clusteringGroup != "training"] = 1

        subsetDataSampleFCS$clusteringGroup[subsetDataSampleFCS$clusteringGroup == "training"] = 0
      }

      subsetDataSampleFCS$name = as.numeric(totalCorrespondenceSamplesID[totalCorrespondenceSamplesID[, 1] == as.numeric(unique(subsetDataSample$group)) & totalCorrespondenceSamplesID[, 4] == currentSample, 5])
      colnameToReplaceID = which(colnames(subsetDataFCS) == "name")
      colnames(subsetDataSampleFCS)[colnameToReplaceID] = "ID"

	  if(nrow(subsetDataSampleFCS) > 1)
	  {
	    subsetDataSampleFCS = as.matrix(apply(subsetDataSampleFCS, 2, as.numeric))
	  } else
	  {
	  	subsetDataSampleFCS = t(as.matrix(apply(subsetDataSampleFCS, 2, as.numeric)))
	  }

      subsetDataSampleFCS = flowCore::flowFrame(subsetDataSampleFCS)

      flowCore::write.FCS(subsetDataSampleFCS, filename = file.path(workingDirectory, "output", "6_FCS", "samples", datasetFolder, paste("Group-", currentGroup, "_Sample-", currentSample, "_", datasetFolder, ".fcs", sep = "")), delimiter = "|")
    }

    tcltk::setTkProgressBar(pb, a, label = paste("Group ", a, "/", groupsToPlot, sep = ""))

    gc()
  }

  close(pb)

  fullTotalData = flowCore::flowFrame(fullTotalData)

  flowCore::write.FCS(fullTotalData, filename = file.path(workingDirectory, "output", "6_FCS", paste("FullData_", datasetFolder, ".fcs", sep = "")), delimiter = "|")

  colnames(totalCorrespondenceSamplesID) = c("GroupNumber", "GroupName", "SampleNumber", "SampleName", "SampleID")
  rownames(totalCorrespondenceSamplesID) = NULL

  utils::write.table(totalCorrespondenceSamplesID, file.path("output", "6_FCS", paste("FCS_Groups&SamplesIdentificationIDs_", datasetFolder, ".txt", sep = "")), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}
