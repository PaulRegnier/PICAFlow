#' Open clusters abundance files
#'
#' This function allows to open one or several clusters abundance files. As this function is called by `mergeData()` function, it should not be called as is.
#'
#' @param pattern A string defining the pattern used to determine which files should be opened. Defaults to `NULL`.
#'
#' @param sampleNamesColumn A string defining the column name the be used as common comparator for each file to merge. This column should be present in each file to be opened and should present identical values (for instance, individualized sample names). Defaults to `NULL`.
#'
#' @return Return value is directly passed to `mergeData()` function.

openFiles = function(pattern = NULL, sampleNamesColumn = NULL)
{
  a = NULL

  columnsToRemove = c("Sample")

  filesToMergeNames = list.files(file.path("output", "8_Analysis"), pattern = pattern)
  totalRownames = NULL
  totalColnames = NULL
  totalData = list()

  foreach::foreach(a = 1:length(filesToMergeNames)) %do%
  {
    currentFileToOpen = filesToMergeNames[a]

    currentFileData = as.data.frame(utils::read.csv(file.path("output", "8_Analysis", currentFileToOpen), header = TRUE, sep = "\t"))

    currentRownames = currentFileData[, sampleNamesColumn]
    currentColnames = colnames(currentFileData)

    currentColnamesID = which(currentColnames %in% columnsToRemove == FALSE)

    if (length(currentColnamesID) > 0)
    {
      currentFileData = currentFileData[, currentColnamesID]
    }

    totalRownames = c(totalRownames, currentFileData[, sampleNamesColumn])

    currentColnamesToAddID = which(colnames(currentFileData) != sampleNamesColumn)
    currentColnamesToAdd = colnames(currentFileData)[currentColnamesToAddID]

    totalColnames = c(totalColnames, currentColnamesToAdd)

    rownames(currentFileData) = currentFileData[, sampleNamesColumn]
    currentFileData = currentFileData[, which(colnames(currentFileData) != sampleNamesColumn)]
    totalData[[currentFileToOpen]] = currentFileData
  }

  totalColnames = unique(totalColnames)
  totalColnames = totalColnames[order(totalColnames)]

  totalRownames = unique(totalRownames)
  totalRownames = totalRownames[order(totalRownames)]

  return(list(totalData = totalData, totalColnames = totalColnames, totalRownames = totalRownames))
}

#' Open and merge clusters abundance files
#'
#' This function allows to open and merge one or several clusters abundance files. Please note that the input text file(s) should be generated with the `exportClustersStatsAndPlots()` function and placed in the `output > 8_Analysis` directory.
#'
#' @param pattern A string defining the pattern used to determine which files should be opened. Defaults to `NULL`.
#'
#' @param sampleNamesColumn A string defining the column name the be used as common comparator for each file to merge. This column should be present in each file to be opened and should present identical values (for instance, individualized sample names). Defaults to `NULL`.
#'
#' @return A data frame containing the merged data.
#'
#' @importFrom foreach %do%
#'
#' @export

mergeData = function(pattern = NULL, sampleNamesColumn = NULL)
{
  b = NULL
  d = NULL
  c = NULL

  dataToMerge = openFiles(pattern = pattern, sampleNamesColumn = sampleNamesColumn)

  totalColnames = dataToMerge$totalColnames
  totalRownames = dataToMerge$totalRownames
  totalData = dataToMerge$totalData

  mergedData = data.frame(matrix(NA, ncol = length(totalColnames), nrow = length(totalRownames)), stringsAsFactors = FALSE)
  rownames(mergedData) = totalRownames
  colnames(mergedData) = totalColnames

  foreach::foreach(b = 1:length(totalRownames)) %do%
  {
    currentPatient = totalRownames[b]
    currentPatientData = NULL

    foreach::foreach(c = 1:length(totalData)) %do%
    {
      currentListData = totalData[[c]]
      currentListPatientData = currentListData[rownames(currentListData) == currentPatient, ]

      foreach::foreach(d = 1:ncol(currentListPatientData)) %do%
      {
        currentColumnName = colnames(currentListPatientData)[d]
        currentColumnPatientData = as.numeric(currentListPatientData[currentColumnName])

        mergedData[currentPatient, currentColumnName] = currentColumnPatientData
      }
    }
  }

  return(mergedData)
}

#' Open and merge metadata to dataset
#'
#' This function allows to open and merge metadata (such as clinical data) to the dataset. Please note that the metadata file should be an Excel spreadsheet named `metadata.xslx` and placed in the `output > 8_Analysis` directory.
#'
#' @param data A data frame containing the original dataset to which metadata are to be applied. This argument should be the output of the `mergeData()` function. Defaults to `NULL`.
#'
#' @param mergeColumn A string defining the column name the be used as comparator for the metadata. This column should be present in the metadata file and contain values that match the rownames of `data` (for instance, individualized sample names). Defaults to `NULL`.
#'
#' @param replaceColumn A string defining the column name the be used as the final name for samples. This column should be present in the metadata file. Defaults to `NULL`.
#'
#' @return A data frame containing the merged data.
#'
#' @importFrom foreach %do%
#'
#' @export

mergeMetadata = function(data = NULL, mergeColumn = NULL, replaceColumn = NULL)
{
  a = NULL

  clinical = as.data.frame(readxl::read_excel(file.path("output", "8_Analysis", "metadata.xlsx")))

  uniquePatientsCodes = unique(rownames(data))

  totalAssociatedClinicalData = data.frame(matrix(NA, ncol = ncol(clinical), nrow = length(uniquePatientsCodes)), stringsAsFactors = FALSE)
  colnames(totalAssociatedClinicalData) = colnames(clinical)

  patientsNotFound = NULL

  foreach::foreach(a = 1:length(uniquePatientsCodes)) %do%
  {
    currentPatient = uniquePatientsCodes[a]

    currentPatient_matchingID = which(clinical[, mergeColumn] == currentPatient)

    if (length(currentPatient_matchingID) > 0)
    {
      associatedClinicalData = clinical[currentPatient_matchingID, ]
    } else
    {
      associatedClinicalData = rep(NA, ncol(clinical))
      patientsNotFound = c(patientsNotFound, currentPatient)
    }

    totalAssociatedClinicalData[a, ] = associatedClinicalData
  }

  if (length(patientsNotFound) > 0)
  {
    out = print(paste(length(patientsNotFound), " samples (", paste(patientsNotFound, collapse = ", "), ") were not found in the metadata file, please check!", sep = ""))
  } else
  {
    out = cbind(data, totalAssociatedClinicalData)
  }

  if(is.null(replaceColumn) == FALSE)
  {
    rownames(out) = out[, replaceColumn]
  }

  return(out)
}

#' Construct rearranged datasets ready for further analyses
#'
#' This function allows to subset and rearrange a dataset file for subsequent analyses. If needed, desired column can also be extracted. Please note that no data are removed from the dataset, but are rather rearranged. Columns of interest will be put in the `dataSubset` element of the final output list, whereas other columns will be put in the `dataRemoved` element of the final output list.
#'
#' @param data A data frame containing the original dataset to rearrange. This argument should be the output of the `mergeMetadata()` function. Defaults to `NULL`.
#'
#' @param columnsToKeepData A character vector defining the columns of the dataset to keep. This argument was rather designed for actual clusters data, but can be used for other kind of columns. Defaults to `NULL`. Please note that this argument can either be considered as a string (if `isColumnsToKeepDataRegex = FALSE`) or as a regex (if `isColumnsToKeepDataRegex = TRUE`).
#'
#' @param columnsToKeepMetadata A character vector defining other columns of the dataset to keep. This argument was rather designed for metadata, but can be used for other kind of columns. Defaults to `NULL`. Please note that this argument can only be considered as a string but never as a regex.
#'
#' @param isColumnsToKeepDataRegex A boolean defining if the columnsToKeepData argument should be considered as a string (if `isColumnsToKeepDataRegex = FALSE`) or as a regex (if `isColumnsToKeepDataRegex = TRUE`). Defaults to `TRUE`.
#'
#' @return A list of length 2: the `dataSubset` element which contains a data frame of the desired columns of interest and the `dataRemoved` element which contains the remaining columns of the dataset.
#'
#' @importFrom foreach %do%
#'
#' @export

subsetDataUMAP = function(data = NULL, columnsToKeepData = NULL, columnsToKeepMetadata = NULL, isColumnsToKeepDataRegex = TRUE)
{
  a = NULL

  if (isColumnsToKeepDataRegex == TRUE)
  {
    dataColsToKeepID = grep(columnsToKeepData, colnames(data))
  } else
  {
    dataColsToKeepID = which(colnames(data) %in% columnsToKeepData)
  }

  if (length(dataColsToKeepID) == 0)
  {
    dataColsToKeepID = NULL
  }

  clinicalColsToKeepID = which(colnames(data) %in% columnsToKeepMetadata)

  if (length(clinicalColsToKeepID) == 0)
  {
    clinicalColsToKeepID = NULL
  }

  totalColsToKeepID = c(dataColsToKeepID, clinicalColsToKeepID)
  totalColsToRemoveID = which((1:ncol(data) %in% totalColsToKeepID) == FALSE)

  if (length(totalColsToRemoveID) > 0)
  {
    dataRemoved = data[, totalColsToRemoveID]
  }

  if (length(totalColsToKeepID) > 0)
  {
    data = data[, totalColsToKeepID]
  }

  rownamesData = rownames(data)
  data = apply(data, 2, as.numeric)
  data = as.data.frame(data)
  rownames(data) = rownamesData

  linesToRemove = NULL

  foreach::foreach(a = 1:nrow(data)) %do%
  {
    currentLineData = data[a, ]
    currentLineNA = length(which(is.na(currentLineData)))

    if (currentLineNA > 0)
    {
      linesToRemove = c(linesToRemove, a)
    }
  }

  if (length(linesToRemove) > 0)
  {
    data = data[-linesToRemove, ]
    dataRemoved = dataRemoved[-linesToRemove, ]
  }

  dataOut = list(dataSubset = data, dataRemoved = dataRemoved)

  return(dataOut)
}

#' Perform UMAP analysis of columns of interest from specified dataset
#'
#' This function allows to perform UMAP analysis of the columns of interest from specified dataset.
#'
#' @param data A list containing the rearranged dataset. This argument should be the output of the `subsetDataUMAP()` function. Defaults to `NULL`.
#'
#' @param n_neighbors_UMAP An integer defining the number of neighbors to use for UMAP analysis. Defaults to `NULL`.
#'
#' @param min_dist_UMAP An integer defining the minimum distance to use for UMAP analysis. Defaults to `NULL`.
#'
#' @param feature A string defining the feature used for overlay and ellipses construction to the UMAP plot. This feature should be a column name of the `dataRemoved` element of `data`. Defaults to `NULL`.
#'
#' @param suffix A string defining a custom suffix to add to the plot name. Defaults to `NULL`.
#'
#' @param returnUMAPData A boolean defining if UMAP parameters (`UMAP_1` and `UMAP_2`) should be added to the `dataRemoved` element of `data`. Defaults to `FALSE`.
#'
#' @param computeUMAP A boolean defining if actual UMAP computation should be performed or not. Defaults to `TRUE`.
#'
#' @param seedValue An integer defining the seed value to use for UMAP computation. Setting a seed helps to obtain reproducible UMAP plots. Defaults to `42`.
#'
#' @param plotHighlightedFeatureItems A boolean defining if UMAP plots displaying highlighted items of the selected `feature` should also be generated or not. Defaults to `TRUE`.
#'
#' @param drawEllipses A boolean defining if ellipses on UMAP plots should be drawn or not. Defaults to `TRUE`.
#'
#' @return A list of length 2: the `dataSubset` element containing a data frame of the desired columns of interest and the `dataRemoved` element containing the remaining columns of the dataset. Generated PDF files are saved to `output > 8_Analysis` directory.
#'
#' @export

UMAP_clusters = function(data = NULL, n_neighbors_UMAP = NULL, min_dist_UMAP = NULL, feature = NULL, suffix = NULL, returnUMAPData = FALSE, computeUMAP = TRUE, seedValue = 42, plotHighlightedFeatureItems = TRUE, drawEllipses = TRUE)
{
  set.seed(seedValue)

  UMAP_1 = NULL
  UMAP_2 = NULL
  a = NULL

  dataUMAP = data$dataSubset

  if (computeUMAP == TRUE)
  {
    UMAP_out = uwot::umap(data$dataSubset, n_neighbors = n_neighbors_UMAP, min_dist = min_dist_UMAP, metric = "euclidean", init = "random", scale = FALSE, verbose = TRUE)

    dataUMAP$UMAP_1 = UMAP_out[, 1]
    dataUMAP$UMAP_2 = UMAP_out[, 2]
  } else
  {
    dataUMAP$UMAP_1 = data$dataRemoved[, "UMAP_1"]
    dataUMAP$UMAP_2 = data$dataRemoved[, "UMAP_2"]
  }

  dataUMAP$sample = gsub("Group-(.+)_Sample-(.+)", "\\2", rownames(data$dataSubset))

  dataUMAP[, feature] = data$dataRemoved[, feature]


	if(drawEllipses == TRUE)
	{

	ggplot2::ggplot(dataUMAP, ggplot2::aes(x = UMAP_1, y = UMAP_2, linewidth = 8, color = get(feature))) +
    ggplot2::geom_point() +
    ggplot2::geom_text(ggplot2::aes(label = sample), hjust = -0.4, vjust = 0.4, size = 2.5, color = "black") +
    ggplot2::stat_ellipse(ggplot2::aes(x = UMAP_1, y = UMAP_2, color = get(feature)), type = "norm", level = 0.5, linewidth = 1.5)

	} else
	{
	ggplot2::ggplot(dataUMAP, ggplot2::aes(x = UMAP_1, y = UMAP_2, linewidth = 8, color = get(feature))) +
    ggplot2::geom_point() +
    ggplot2::geom_text(ggplot2::aes(label = sample), hjust = -0.4, vjust = 0.4, size = 2.5, color = "black")
	}

  currentPlotName = paste("UMAP_Overlay-", feature, suffix, ".pdf", sep = "")
  totalPlotFileName = file.path("output", "8_Analysis", currentPlotName)
  ggplot2::ggsave(totalPlotFileName, width = 10, height = 10)

  if (plotHighlightedFeatureItems == TRUE)
  {
	  featureItems = unique(dataUMAP[, feature])

  	foreach::foreach(a = 1:length(featureItems)) %do%
  	{
  		currentFeatureItem = featureItems[a]
  		dataUMAP_temp = dataUMAP

  		dataUMAP_temp[dataUMAP_temp[, feature] != currentFeatureItem, feature] = "Other"


if(drawEllipses == TRUE)
	{
	ggplot2::ggplot(dataUMAP_temp, ggplot2::aes(x = UMAP_1, y = UMAP_2, linewidth = 8, color = get(feature))) +
  		ggplot2::geom_point() +
  		ggplot2::geom_text(ggplot2::aes(label = sample), hjust = -0.4, vjust = 0.4, size = 2.5, color = "black") +
  	  ggplot2::stat_ellipse(ggplot2::aes(x = UMAP_1, y = UMAP_2, color = get(feature)), type = "norm", level = 0.5, linewidth = 1.5)


	} else
	{
	ggplot2::ggplot(dataUMAP_temp, ggplot2::aes(x = UMAP_1, y = UMAP_2, linewidth = 8, color = get(feature))) +
  		ggplot2::geom_point() +
  		ggplot2::geom_text(ggplot2::aes(label = sample), hjust = -0.4, vjust = 0.4, size = 2.5, color = "black")
	}

  	  currentPlotName_outlined = paste("UMAP_Overlay-", feature, suffix, "_", currentFeatureItem, "_highlighted.pdf", sep = "")
  	  totalPlotFileName_outlined = file.path("output", "8_Analysis", currentPlotName_outlined)
  	  ggplot2::ggsave(totalPlotFileName_outlined, width = 10, height = 10)
  	}
  }

  if (returnUMAPData == TRUE)
  {
    data$dataRemoved[, "UMAP_1"] = UMAP_out[, 1]
    data$dataRemoved[, "UMAP_2"] = UMAP_out[, 2]

    return(data)
  }
}

#' Remove outliers from dataset
#'
#' This function allows to remove from the dataset samples that are considered as outliers. Please note that these samples will be permanently removed and cannot be recovered afterwards.
#'
#' @param data A list containing the merged dataset. This argument should be the output of the `mergeMetadata()` function. Defaults to `NULL`.
#'
#' @param outliers An characted vector defining the outliers to remove from the dataset. Defaults to `NULL`.
#'
#' @return A data frame containing `data` minus the outliers.
#'
#' @export

removeOutliers = function(data = NULL, outliers = NULL)
{
  outliersToRemoveID = which(rownames(data) %in% outliers)

  if (length(outliersToRemoveID) > 0)
  {
    data = data[-outliersToRemoveID, ]
  }

  return(data)
}

#' Perform hierarchical clustering of UMAP parameters
#'
#' This function allows to perform hierarchical clustering of UMAP parameters. Please note that the determined clusters for each cell will be saved as a new column named `cluster` in the `dataRemoved` element of `data`.
#'
#' @param data A list containing the rearranged dataset. This argument should be the output of the `UMAP_clusters()` function with `addUMAPOutputToData = TRUE`. Defaults to `NULL`.
#'
#' @param clustersNb An integer defining the number expected/observable/supposed clusters from the UMAP plot. Defaults to `NULL`.
#'
#' @return A list of length 2: the `dataSubset` element containing a data frame of the desired columns of interest and the `dataRemoved` element containing the remaining columns of the dataset.
#'
#' @export

hierarchicalClusteringData = function(data = NULL, clustersNb = NULL)
{
  clusteringData = data$dataRemoved[, c("UMAP_1", "UMAP_2")]

  d = stats::dist(clusteringData, method = "euclidean")
  fit = stats::hclust(d, method = "complete")

  groups = stats::cutree(fit, k = clustersNb)

  grDevices::pdf(file = file.path("output", "8_Analysis", "clusteringDendrogram.pdf"), paper = "a4r", width = 25, height = 15)
  plot(fit, cex = 0.25)

  graphics::par(mar = c(0, 0, 0, 0))
  stats::rect.hclust(fit, k = clustersNb, border = "red")
  grDevices::dev.off()
  data$dataRemoved[, "cluster"] = as.character(groups)

  return(data)
}

#' Bind data to a single data frame
#'
#' This function allows transform the rearranged dataset to simple data frame, to make subsequent data exportation and analysis easier.
#'
#' @param data A list containing the rearranged dataset. This argument should be the output of the `hierarchicalClusteringData()`. Defaults to `NULL`.
#'
#' @return A data frame containing the whole columns from `data`.
#'
#' @export

bindData = function(data = NULL)
{
  data = cbind(data$dataSubset, data$dataRemoved)

  return(data)
}


#' Construct abundance plots for desired features
#'
#' This function allows to generate and export plots showing the cluster abundance for one or several desired feature(s) of interest, both using boxplots and UMAP overlays. A feature can be for instance a group or a metadata element.
#'
#' @param data A data frame containing the dataset to use. This argument should be the output of the `bindData()` function. Defaults to `NULL`.
#'
#' @param columnsToPlot A character vector defining the columns to be used as plotted elements (typically cell clusters). Defaults to `NULL`.
#'
#' @param features A character vector defining the columns to be used as features (typically group, cluster or metadata). Defaults to `NULL`.
#'
#' @param plotUMAPOverlays A boolean defining if the UMAP overlays should be plotted (`TRUE`) or not (`FALSE`). This is typically used if the user does not intend to use the metadata addition and analysis, but only wants the boxplots to be produced. Defaults to `TRUE`.
#'
#' @return Generated PDF files are saved to `output > 8_Analysis > columnsToPlot` subdirectories.
#'
#' @importFrom foreach %do%
#'
#' @export

constructPlots = function(data = NULL, columnsToPlot = NULL, features = NULL, plotUMAPOverlays = TRUE)
{
  UMAP_1 = NULL
  UMAP_2 = NULL
  mean_se = NULL
  b = NULL

  foreach::foreach(b = 1:length(columnsToPlot)) %do%
  {
    currentColumn = columnsToPlot[b]

    dir.create(file.path("output", "8_Analysis", currentColumn))

    if(plotUMAPOverlays == TRUE)
    {
      currentMidpointValue = (min(data[, currentColumn]) + max(data[, currentColumn]))/2
      ggplot2::ggplot(data, ggplot2::aes(x = UMAP_1, y = UMAP_2, size = 8, color = get(currentColumn))) +
        ggplot2::scale_colour_gradient2(low = "blue", high = "red", mid = "yellow", midpoint = currentMidpointValue) +
        ggplot2::geom_point()

      currentFileName_plot1 = file.path("output", "8_Analysis", currentColumn, paste("UMAP_Overlay-", currentColumn, ".pdf", sep = ""))
      ggplot2::ggsave(currentFileName_plot1, width = 10, height = 10)


    }


    foreach::foreach(c = 1:length(features)) %do%
    {
      currentFeature = features[c]

      currentFeatureDataToRemoveID = which(data[, currentFeature] == "Excluded")

      currentFeatureData = data

      if (length(currentFeatureDataToRemoveID) > 0)
      {
        currentFeatureData = data[-currentFeatureDataToRemoveID, ]
      }

      ggplot2::ggplot(currentFeatureData, ggplot2::aes(x = get(currentFeature), y = get(currentColumn), group = get(currentFeature))) +
        ggplot2::geom_boxplot(outlier.shape = NA, col = "darkblue") +
        ggplot2::stat_boxplot(geom = "errorbar", width = 0.25, col = "darkblue") +
        ggplot2::stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.5, linewidth = 1.5, col = "red") +
        ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 0.25, linewidth = 0.75, col = "red") +
        ggplot2::geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, binwidth = diff(range(currentFeatureData[, currentColumn]))/80, method = "histodot")

      currentFileName_plot2 = file.path("output", "8_Analysis", currentColumn, paste("Boxplots_", currentFeature, "_", currentColumn, ".pdf", sep = ""))
      ggplot2::ggsave(currentFileName_plot2, width = 10, height = 10)
    }
  }
}

#' Export data bind
#'
#' This function allows to export the final data bind to a tabular text file.
#'
#' @param data A data frame containing the dataset to export. This argument should be the output of the `bindData()` function. Defaults to `NULL`.
#'
#' @return Generated text file is saved to `output > 8_Analysis` directory.
#'
#' @export

exportDataBind = function(data = NULL)
{
  dataUMAP_export = cbind(rownames(data), data)
  colnames(dataUMAP_export)[1] = "Samples"

  utils::write.table(dataUMAP_export, file = file.path("output", "8_Analysis", paste("FullData.txt", sep = "")), sep = "\t", col.names = TRUE, row.names = FALSE)
}

#' Export global tables for defined features
#'
#' This function allows to export text files which contain several statistics about columns of interest using selected features of interest as data delimiters.
#'
#' @param data A data frame containing the dataset to export. This argument should be the output of the `bindData()` function. Defaults to `NULL`.
#'
#' @param columnsToUse A character vector defining the columns considered as the variable of interest (typically cell clusters). Defaults to `NULL`.
#'
#' @param featuresToUse A character vector defining the columns to be used as features (typically group, cluster or metadata). Defaults to `NULL`.
#'
#' @param metricToUse A string defining the metric to use for statistics computation. It can be either `mean` or `median`. Defaults to `median`.
#'
#' @return Generated text files are saved to `output > 8_Analysis` directory.
#'
#' @importFrom foreach %do%
#'
#' @export

exportFeaturesTables = function(data = NULL, columnsToUse = NULL, featuresToUse = NULL, metricToUse = "median")
{
  f = NULL

  foreach::foreach(f = 1:length(featuresToUse)) %do%
  {
    currentFeature = featuresToUse[f]

    currentFeatureAnalysis = analyzeFeature(data = data, columnsToUse = columnsToUse, feature = currentFeature, metricToUse = metricToUse)

    currentFeatureAnalysis_export = cbind(rownames(currentFeatureAnalysis), currentFeatureAnalysis)
    currentFeatureAnalysis_export = rbind(colnames(currentFeatureAnalysis_export), currentFeatureAnalysis_export)
    currentFeatureAnalysis_export[1, 1] = paste("Feature_", currentFeature, sep = "")

    utils::write.table(currentFeatureAnalysis_export, file = file.path("output", "8_Analysis", paste("Feature_", currentFeature, "_analysisOut.txt", sep = "")), sep = "\t", col.names = FALSE, row.names = FALSE)
  }
}

#' Analyze a given feature
#'
#' This function allows to analyze a given feature. As this function is called by `exportFeaturesTables()` function, it should not be called as is.
#'
#' @param data A data frame containing the dataset to export. This argument should be the output of the `bindData()` function. Defaults to `NULL`.
#'
#' @param columnsToUse A character vector defining the columns considered as the variable of interest (typically cell clusters). Defaults to `NULL`.
#'
#' @param feature A character vector defining the column to be used as a feature (typically group, cluster or metadata). Defaults to `NULL`.
#'
#' @param metricToUse A string defining the metric to use for statistics computation. It can be either `mean` or `median`. Defaults to `median`.
#'
#' @return Return value is directly passed to `exportFeaturesTables()` function.

analyzeFeature = function(data = NULL, columnsToUse = NULL, feature = NULL, metricToUse = "median")
{
  a = NULL
  b = NULL
  d = NULL
  e = NULL
  f = NULL

  rowsToRemoveID = which(data[, feature] == "Excluded")

  if (length(rowsToRemoveID) > 0)
  {
    data = data[-rowsToRemoveID, ]
  }

  uniqueFeatures = unique(data[, feature])

  totalPairsList = list()

  foreach::foreach(a = 1:length(uniqueFeatures)) %do%
  {
    currentGroupName_1 = uniqueFeatures[a]

    foreach::foreach(b = 1:length(uniqueFeatures)) %do%
    {
      currentGroupName_2 = uniqueFeatures[b]
      totalPairsList[[paste(currentGroupName_1, currentGroupName_2, sep = "_")]] = sort(c(currentGroupName_1, currentGroupName_2))
    }
  }

  elementsToRemoveID = which(duplicated(totalPairsList))

  foreach::foreach(c = 1:length(totalPairsList)) %do%
  {
    currentElementToTest = totalPairsList[[c]]

    if (length(unique(currentElementToTest)) == 1)
    {
      elementsToRemoveID = c(elementsToRemoveID, c)
    }
  }

  totalPairsList = totalPairsList[-elementsToRemoveID]

  totalColumnsResults_pVal = data.frame(matrix(NA, ncol = length(totalPairsList), nrow = length(columnsToUse)), stringsAsFactors = FALSE)
  colnames(totalColumnsResults_pVal) = names(totalPairsList)
  rownames(totalColumnsResults_pVal) = columnsToUse

  totalColumnsResults_diffMetric = totalColumnsResults_pVal

  foreach::foreach(d = 1:length(columnsToUse)) %do%
  {
    currentColumnToTest = columnsToUse[d]

    currentColumnsTestResults_pVal = totalPairsList
    currentColumnsTestResults_diffMetric = totalPairsList

    foreach::foreach(e = 1:length(totalPairsList)) %do%
    {
      currentPairToTest = totalPairsList[[e]]

      currentPairToTest_1_ID = which(data[, feature] == currentPairToTest[1])
      currentPairToTest_2_ID = which(data[, feature] == currentPairToTest[2])

      dataset_x = data[currentPairToTest_1_ID, currentColumnToTest]
      dataset_y = data[currentPairToTest_2_ID, currentColumnToTest]

      if (length(dataset_x) >= 3 & length(dataset_y) >= 3)
      {
        currentTestResult = stats::t.test(x = dataset_x, y = dataset_y)
        currentColumnsTestResults_pVal[[e]] = currentTestResult$p.value
      } else
      {
        currentColumnsTestResults_pVal[[e]] = NA
      }

      if (metricToUse == "mean")
      {
        currentColumnsTestResults_diffMetric[[e]] = mean(dataset_x) - mean(dataset_y)
      } else if (metricToUse == "median")
      {
        currentColumnsTestResults_diffMetric[[e]] = stats::median(dataset_x) - stats::median(dataset_y)
      } else
      {
        currentColumnsTestResults_diffMetric[[e]] = mean(dataset_x) - mean(dataset_y)
      }
    }

    totalColumnsResults_pVal[currentColumnToTest, ] = currentColumnsTestResults_pVal
    totalColumnsResults_diffMetric[currentColumnToTest, ] = currentColumnsTestResults_diffMetric
  }

  uncorrectedPValues = as.vector(as.matrix(totalColumnsResults_pVal))

  correctedPValues = stats::p.adjust(uncorrectedPValues, method = "BH")

  totalColumnsResults_pValuesCorrected = data.frame(matrix(correctedPValues, ncol = length(names(totalPairsList))), stringsAsFactors = FALSE)
  rownames(totalColumnsResults_pValuesCorrected) = rownames(totalColumnsResults_pVal)
  colnames(totalColumnsResults_pValuesCorrected) = paste("adjPVal_Comp_", feature, "_", colnames(totalColumnsResults_pVal), sep = "")

  colnames(totalColumnsResults_diffMetric) = paste("DiffMetric_Comp_", feature, "_", colnames(totalColumnsResults_pVal), sep = "")

  totalColumnsResults_score = totalColumnsResults_diffMetric
  colnames(totalColumnsResults_score) = gsub("adjPVal", "Score", colnames(totalColumnsResults_pValuesCorrected))
  totalColumnsResults_score$sumScore = NA

  foreach::foreach(f = 1:nrow(totalColumnsResults_pValuesCorrected)) %do%
  {
    currentRowData_pVal = totalColumnsResults_pValuesCorrected[f, ]

    currentRowData_diffMetric = abs(totalColumnsResults_diffMetric[f, ])

    currentRowData_score = currentRowData_diffMetric * (1 - currentRowData_pVal)

    totalColumnsResults_score[f, ] = cbind(currentRowData_score, sum(currentRowData_score, na.rm = TRUE))
  }

  currentFeatureRawMetrics = computeMetrics(data = data, columnsToUse = columnsToUse, featureToUse = feature, metricToUse = metricToUse)

  totalFeatureResults = cbind(currentFeatureRawMetrics, totalColumnsResults_pValuesCorrected, totalColumnsResults_diffMetric, totalColumnsResults_score)
  totalFeatureResults = totalFeatureResults[order(totalFeatureResults$sumScore, decreasing = TRUE), ]

  return(totalFeatureResults)
}

#' Compute metrics for the given feature
#'
#' This function allows to compute metrics for the given feature. As this function is called by `exportFeaturesTables()` function, it should not be called as is.
#'
#' @param data A data frame containing the dataset to export. This argument should be the output of the `bindData()` function. Defaults to `NULL`.
#'
#' @param columnsToUse A character vector defining the columns considered as the variable of interest (typically cell clusters). Defaults to `NULL`.
#'
#' @param featureToUse A character vector defining the column to be used as feature (typically group, cluster or metadata). Defaults to `NULL`.
#'
#' @param metricToUse A string defining the metric to use for statistics computation. It can be either `mean` or `median`. Defaults to `median`.
#'
#' @return Return value is directly passed to `exportFeaturesTables()` function.

computeMetrics = function(data = NULL, columnsToUse = NULL, featureToUse = NULL, metricToUse = "median")
{
  a = NULL
  b = NULL

  uniqueFeatures = unique(data[, featureToUse])

  dataMetrics = data.frame(matrix(NA, ncol = length(uniqueFeatures), nrow = length(columnsToUse)), stringsAsFactors = FALSE)
  colnames(dataMetrics) = uniqueFeatures
  rownames(dataMetrics) = columnsToUse

  foreach::foreach(a = 1:length(columnsToUse)) %do%
  {
    currentColumn = columnsToUse[a]
    currentAssociatedData = data[, c(currentColumn, featureToUse)]

    currentColumnMetricsPerFeature = NULL

    foreach::foreach(b = 1:length(uniqueFeatures)) %do%
    {
      currentFeature = uniqueFeatures[b]

      currentFeatureAssociatedDataID = which(currentAssociatedData[, featureToUse] == currentFeature)

      if (metricToUse == "mean")
      {
        currentColumnMetricsPerFeature = c(currentColumnMetricsPerFeature, mean(currentAssociatedData[currentFeatureAssociatedDataID, currentColumn]))
      } else if (metricToUse == "median")
      {
        currentColumnMetricsPerFeature = c(currentColumnMetricsPerFeature, stats::median(currentAssociatedData[currentFeatureAssociatedDataID, currentColumn]))
      } else
      {
        currentColumnMetricsPerFeature = c(currentColumnMetricsPerFeature, mean(currentAssociatedData[currentFeatureAssociatedDataID, currentColumn]))
      }
    }

    dataMetrics[currentColumn, ] = currentColumnMetricsPerFeature
  }

  colnames(dataMetrics) = paste("metric_", featureToUse, "_", colnames(dataMetrics), sep = "")

  return(dataMetrics)
}

#' Export clusters abundance heatmap using a given feature of interest
#'
#' This function allows to export a heatmap showing clusters abundance using a given feature of interest as delimiter.
#'
#' @param feature A string defining the column to be used a feature (typically group, cluster or metadata). Defaults to `NULL`.
#'
#' @param clustersToKeepRegex A string defining the regex to use for the identification of columns associated to cell clusters. Defaults to `NULL`.
#'
#' @param metricToUse A string defining the metric to use for statistics computation. It can be either `mean` or `median`. Defaults to `median`.
#'
#' @return Generated PDF file is saved to `output > 8_Analysis` directory.
#'
#' @export

heatmapAbundancesGroups = function(feature = NULL, clustersToKeepRegex = NULL, metricToUse = "median")
{
  data_abundance = utils::read.table(file.path("output", "8_Analysis", paste("Feature_", feature, "_analysisOut.txt", sep = "")), header = TRUE, sep = "\t")

  rownames(data_abundance) = data_abundance[, paste("Feature_", feature, sep = "")]
  data_abundance[, paste("Feature_", feature, sep = "")] = NULL

  clustersToKeepID = grep(clustersToKeepRegex, rownames(data_abundance))
  data_abundance = data_abundance[clustersToKeepID, ]

  dataColsID = grep(("metric_"), colnames(data_abundance))
  data_abundance = as.matrix(t(data_abundance[, dataColsID]))

  data_abundance = data_abundance[, order(colnames(data_abundance))]

  columnsToRemove = as.numeric(which(is.na(colSums(data_abundance))))

  if(length(columnsToRemove) > 0)
  {
    samplesToRemove = colnames(data_abundance)[columnsToRemove]
	  data_abundance = data_abundance[, -c(columnsToRemove)]

	  print(paste("Warning! The following cluster(s) was/were removed from the heatmap representation, mainly because all the metrics inside were identical for all the groups:", paste(samplesToRemove, collapse = ", ", sep = "")))
  }

  columnsToRemove = as.numeric(which(colSums(data_abundance) == 0))

  if(length(columnsToRemove) > 0)
  {
    samplesToRemove = colnames(data_abundance)[columnsToRemove]
	  data_abundance = data_abundance[, -c(columnsToRemove)]

	  print(paste("Warning! The following cluster(s) was/were removed from the heatmap representation, mainly because all the metrics inside were equal to 0 for all the groups:", paste(samplesToRemove, collapse = ", ", sep = "")))
  }

  data_abundance = scale(data_abundance)

  min = min(data_abundance, na.rm = TRUE)
  max = max(data_abundance, na.rm = TRUE)

  if (metricToUse == "mean")
  {
    middle = mean(data_abundance, na.rm = TRUE)
  } else if (metricToUse == "median")
  {
    middle = stats::median(data_abundance, na.rm = TRUE)
  } else
  {
    middle = mean(data_abundance, na.rm = TRUE)
  }

	separator = min(abs(min-middle), abs(max-middle))/2

  colors = c(seq(min, (middle - separator) - 0.001, length = 25), seq(middle - separator, middle + separator, length = 50), seq((middle + separator) + 0.001, max, length = 25))
  colors = sort(colors)
  customPalette = (grDevices::colorRampPalette(c("lightblue", "blue", "black", "yellow", "orange")))(n = 99)

  grDevices::pdf(file = file.path("output", "8_Analysis", paste("heatmap_", feature, "_", clustersToKeepRegex, "_unscaled.pdf", sep = "")), bg = "transparent", width = 20, height = 10, paper = "a4r")

  gplots::heatmap.2(data_abundance, trace = "none", scale = "none", key = TRUE, keysize = 1, density.info = "none", col = customPalette, breaks = colors, cexRow = 0.5, cexCol = 0.8, margins = c(10, 10), dendrogram = "both")

  grDevices::dev.off()
}

#' Perform ROC analysis
#'
#' This function allows to perform ROC analysis on the full data exported using the `exportDataBind()` function. To use this function, please create a Excel copy (`xslx` format) of the `fullData.txt` file. One can use this Excel file to easily create new features, notably by combining several other ones (for instance by summing clusters).
#'
#' @param dataFile A string defining the Excel file to use. Defaults to `NULL`.
#'
#' @param predictors A list defining the predictors to use. Each element of the list is considered as a predictor. If an element of the `predictors` list is itself a list, the final associated predictor will be the sum of the elements inside this list. Defaults to `NULL`.
#'
#' @param pairsToAnalyze A list defining the pairs to use for the ROC analysis. Each element of the `pairsToAnalyze` list should also be a list of length 2 where both sub-elements are a vector representing the paired groups. Each of these sub-elements can be a vector of the desired length. In the case the vector length of a sub-element is greater than 1, the final group will consist of the pool of every sub-element put inside. Defaults to `NULL`.
#'
#' @param pairs_columnToCheck A string defining the column to use for pairs identification. Defaults to `NULL`.
#'
#' @return Generated text and PDF files are saved to `output > 8_Analysis > ROC` directory.
#'
#' @importFrom foreach %do%
#'
#' @export

ROCanalysis = function(dataFile = NULL, predictors = NULL, pairsToAnalyze = NULL, pairs_columnToCheck = NULL)
{
  a = NULL
  i = NULL
  b = NULL

  if (dir.exists(file.path("output", "8_Analysis", "ROC")) == FALSE)
  {
    dir.create(file.path("output", "8_Analysis", "ROC"))
  } else
  {
    unlink(file.path("output", "8_Analysis", "ROC", "*.*"))

  }

  data = data.frame(readxl::read_excel(file.path("output", "8_Analysis", dataFile)))

  finalData = data.frame(matrix(ncol = 6))

  foreach::foreach(i = 1:length(predictors)) %do%
  {
    currentPredictor = predictors[[i]]

    grDevices::pdf(file.path("output", "8_Analysis", "ROC", paste("ROC_", paste(currentPredictor, collapse = "+"), ".pdf", sep = "")))

    graphics::par(mfrow = c(2, 2))

    foreach::foreach(a = 1:length(pairsToAnalyze)) %do%
    {
      currentPair = pairsToAnalyze[[a]]

      if (is.list(currentPair) == FALSE)
      {
        currentData = data[which((data[, pairs_columnToCheck] %in% currentPair) == TRUE), ]
        tagsToUse = currentPair
      } else
      {
        currentData = NULL
        totalDataTags = NULL

        foreach::foreach(b = 1:length(currentPair)) %do%
        {
          currentPair_elements = currentPair[[b]]
          currentPair_tag = paste(currentPair[[b]], collapse = "+")

          foreach::foreach(c = 1:length(currentPair_elements)) %do%
          {
            currentPair_subelement = currentPair_elements[c]
            currentData_subelement = data[which((data[, pairs_columnToCheck] %in% currentPair_subelement) == TRUE), ]

            currentData = rbind(currentData, currentData_subelement)

            totalDataTags = c(totalDataTags, rep(currentPair_tag, nrow(currentData_subelement)))
          }
        }

        currentData[, pairs_columnToCheck] = totalDataTags
        tagsToUse = unique(totalDataTags)
      }

      if (length(currentPredictor) > 1)
      {
        scoreUsed = rowSums(apply(currentData[, currentPredictor], 2, as.numeric), na.rm = TRUE)
      } else
      {
        scoreUsed = currentData[, currentPredictor]
      }

      roc_empirical = ROCit::rocit(score = as.numeric(scoreUsed), class = currentData[, pairs_columnToCheck], negref = tagsToUse[1], method = "empirical")

      plot = graphics::plot(roc_empirical, values = TRUE)
      graphics::text(1, 0.3, paste("Pair = ", paste(tagsToUse, collapse = " VS "), sep = ""), col = "black", cex = 0.35, adj = 1)
      graphics::text(1, 0.25, paste("Predictors = ", paste(currentPredictor, collapse = "+"), sep = ""), col = "black", cex = 0.35, adj = 1)

      associatedData = c(paste(currentPredictor, collapse = "+"), paste(tagsToUse, collapse = "-VS-"), as.numeric(plot$AUC), as.numeric(plot$"optimal Youden Index point"["cutoff"]), as.numeric(plot$"optimal Youden Index point"["FPR"]), as.numeric(plot$"optimal Youden Index point"["TPR"]))

      finalData = rbind(finalData, associatedData)
    }

    grDevices::dev.off()
  }

  finalData = finalData[-1, ]
  colnames(finalData) = c("Predictor", "Pair", "AUC", "OptimalYoudenIndexPoint_cutoff", "OptimalYoudenIndexPoint_FPR", "OptimalYoudenIndexPoint_TPR")

  finalData = rbind(colnames(finalData), finalData)

  utils::write.table(finalData, file = file.path("output", "8_Analysis", "ROC", "ROC_values.txt"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}


#' Export all parameters used
#'
#' This function allows to export all the parameters that were used throughout the entire analysis. Some parameters will not be exported here because they were already exported in a special format when the associated functions were used.
#'
#' @param parametersToExport A list containing all the parameters to export. Defaults to `NULL`.
#'
#' @return Generated text file is saved to `output` directory.
#'
#' @export

exportParametersUsed = function(parametersToExport = NULL)
{
  options(scipen = 999)

  sink(file.path("output", paste("analysisParametersUsed_", parametersToExport$datasetToUse_value, ".txt", sep = "")))

  cat("# Global analysis seed value")
  cat("\n")
  cat(paste("seed_value = ", parametersToExport$seed_value, sep = ""))
  cat("\n\n")

  cat("# Initial parameters to keep in the dataset")
  cat("\n")
  cat(paste("parametersToKeep = c(\"", gsub(", ", "\", \"", paste(parametersToExport$parametersToKeep, collapse = ", ")), "\")", sep = ""))
  cat("\n")
  cat(paste("customNames = c(\"", gsub(", ", "\", \"", paste(parametersToExport$customNames, collapse = ", ")), "\")", sep = ""))
  cat("\n\n")

  cat("# 2_Normalization")
  cat("\n")
  cat(paste("warpSet_value = ", parametersToExport$warpSet_value, sep = ""))
  cat("\n")
  cat(paste("gaussNorm_value = ", parametersToExport$gaussNorm_value, sep = ""))
  cat("\n")
  cat(paste("max.lms.sequence = list(", paste(parametersToExport$max.lms.sequence, collapse = ", "), ")", sep = ""))
  cat("\n")

  if (is.null(parametersToExport$samplesToDelete) == TRUE)
  {
    cat(paste("samplesToDelete = NULL", sep = ""))
  } else
  {
    cat(paste("samplesToDelete = c(", paste(parametersToExport$samplesToDelete, collapse = ", "), ")", sep = ""))
  }

  cat("\n\n")

  cat("# 4_Downsampling")
  cat("\n")
  cat(paste("parametersToKeepFinal = c(\"", gsub(", ", "\", \"", paste(parametersToExport$parametersToKeepFinal, collapse = ", ")), "\")", sep = ""))
  cat("\n")
  cat(paste("downsampleMinEvents_value = ", parametersToExport$downsampleMinEvents_value, sep = ""))
  cat("\n")
  cat(paste("maxCellsNb = ", parametersToExport$maxCellsNb, sep = ""))
  cat("\n")
  cat(paste("estimateThreshold_value = ", parametersToExport$estimateThreshold_value, sep = ""))
  cat("\n")

  if (is.null(parametersToExport$trainingDatasetProportion_value) == TRUE)
  {
    cat(paste("trainingDatasetProportion_value = NULL", sep = ""))
  } else
  {
    cat(paste("trainingDatasetProportion_value = ", parametersToExport$trainingDatasetProportion_value, sep = ""))
  }

  cat("\n")
  cat(paste("datasetToUse_value = \"", parametersToExport$datasetToUse_value, "\"", sep = ""))
  cat("\n\n")

  cat("# 5_UMAP")
  cat("\n")
  cat(paste("min_dist_value = ", parametersToExport$min_dist_value, sep = ""))
  cat("\n")
  cat(paste("n_neighbors_value = ", parametersToExport$n_neighbors_value, sep = ""))
  cat("\n\n")

  cat("# 7_Clustering")
  cat("\n")
  cat(paste("subsetDownsampled_value = ", parametersToExport$subsetDownsampled_value, sep = ""))
  cat("\n")
  cat(paste("clusterMinPercentage_value = ", parametersToExport$clusterMinPercentage_value, sep = ""))
  cat("\n")
  cat(paste("metricToUse_value = \"", parametersToExport$metricToUse_value, "\"", sep = ""))
  cat("\n")
  cat(paste("cutoff_value = ", parametersToExport$cutoff_value, sep = ""))
  cat("\n")
  cat(paste("thresholds = c(", paste(parametersToExport$parametersToUseThresholds$threshold, collapse = ", "), ")", sep = ""))
  cat("\n")
  cat(paste("folder_value = \"", parametersToExport$folder_value, "\"", sep = ""))
  cat("\n")
  cat(paste("prefix_value = \"", parametersToExport$prefix_value, "\"", sep = ""))
  cat("\n")
  cat(paste("maxCellsPerPlot_value = ", parametersToExport$maxCellsPerPlot_value, sep = ""))
  cat("\n\n")

  cat("# 8_Analysis")
  cat("\n")
  cat(paste("sampleNamesColumn = \"", parametersToExport$sampleNamesColumn_value, "\"", sep = ""))
  cat("\n")
  cat(paste("mergeColumn = \"", parametersToExport$mergeColumn_value, "\"", sep = ""))
  cat("\n")
  cat(paste("columnsToKeepData = \"", parametersToExport$columnsToKeepData_value, "\"", sep = ""))
  cat("\n")

  if (is.null(parametersToExport$columnsToKeepClinical_value) == TRUE)
  {
    cat(paste("columnsToKeepMetadata = NULL", sep = ""))
  } else
  {
    cat(paste("columnsToKeepMetadata = c(", gsub(", ", "\", \"", paste(parametersToExport$columnsToKeepClinical_value, collapse = ", ")), ")", sep = ""))
  }

  cat("\n")
  cat(paste("n_neighbors_value_UMAP = ", parametersToExport$n_neighbors_value_UMAP, sep = ""))
  cat("\n")
  cat(paste("min_dist_value_UMAP = ", parametersToExport$min_dist_value_UMAP, sep = ""))
  cat("\n")
  cat(paste("outliersToRemove = c(\"", gsub(", ", "\", \"", paste(parametersToExport$outliersToRemove_value, collapse = ", ")), "\")", sep = ""))
  cat("\n")
  cat(paste("clustersNb = ", parametersToExport$clustersNb_value, sep = ""))
  cat("\n")
  cat(paste("columnsToPlot = c(\"", gsub(", ", "\", \"", paste(parametersToExport$columnsToPlot_value, collapse = ", ")), "\")", sep = ""))
  cat("\n")
  cat(paste("featuresToUse = c(\"", gsub(", ", "\", \"", paste(parametersToExport$featuresToUse_value, collapse = ", ")), "\")", sep = ""))

  sink()
}
