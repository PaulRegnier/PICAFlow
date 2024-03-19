#' Test normality of dataset parameters
#'
#' This function allows to test the normality of each parameter of the dataset using Shapiro-Wilk test. This is useful to determine whether mean or median metric should be used in later analyses. Typically, a skewed distribution will show poorly informative mean values, but much more informative median values.
#'
#' @param data A string defining the dataset to use. Defaults to `NULL`.
#'
#' @param parametersToAnalyze A character vector defining the parameters to analyze in the subsequent dataset. Defaults to `NULL`.
#'
#' @param datasetFolder A string defining which dataset to open. It can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to `full`.
#'
#' @return Generated text file is saved to `output > 7_Clustering > datasetFolder` directory.
#'
#' @importFrom foreach %do%
#'
#' @export

testDatasetNormality = function(data = NULL, parametersToAnalyze = NULL, datasetFolder = "full")
{
  a = NULL

  if (dir.exists(file.path("output", "7_Clustering", datasetFolder)))
  {
    unlink(file.path("output", "7_Clustering", datasetFolder), recursive = TRUE)
  }

  dir.create(file.path("output", "7_Clustering", datasetFolder))

  totalParametersPvaluesShapiro = data.frame(matrix(0, ncol = 3, nrow = length(parametersToAnalyze)))

  foreach::foreach(a = 1:length(parametersToAnalyze)) %do%
  {
    currentParameter = parametersToAnalyze[a]
    currentShapiroTest = stats::shapiro.test(sample(data[, currentParameter], 5000))
    totalParametersPvaluesShapiro[a, 1] = currentShapiroTest$statistic
    totalParametersPvaluesShapiro[a, 2] = currentShapiroTest$p.value

    if (currentShapiroTest$p.value < 0.05)
    {
      totalParametersPvaluesShapiro[a, 3] = "No"
    } else
    {
      totalParametersPvaluesShapiro[a, 3] = "Yes"
    }
  }

  rownames(totalParametersPvaluesShapiro) = parametersToAnalyze
  colnames(totalParametersPvaluesShapiro) = c("ShapiroWilk_statistic", "ShapiroWilk_pValue", "isDistributionNormal")

  print(totalParametersPvaluesShapiro)

  totalParametersPvaluesShapiro = cbind(rownames(totalParametersPvaluesShapiro), totalParametersPvaluesShapiro)
  colnames(totalParametersPvaluesShapiro)[1] = "parameter"

  utils::write.table(totalParametersPvaluesShapiro, file = file.path("output", "7_Clustering", datasetFolder, paste("ShapiroWilkNormalityTest_Summary.txt", sep = "")), quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)
}

#' Perform initial cell clustering training
#'
#' This function allows perform the initial step of the cell clustering training.
#'
#' @param data A data frame containing the dataset to use. Defaults to `NULL`.
#'
#' @param parametersToAnalyze A character vector defining the parameters to analyze in the subsequent dataset. Defaults to `NULL`.
#'
#' @param subsetDownsampled An integer defining the maximum number of cells to use for the analysis. Defaults to `25000`.
#'
#' @param clusterMinPercentage An integer defining the minimum frequency of a cluster to be kept as is. Values represent percentages and range from `0` to `100`. Defaults to `0.5`.
#'
#' @return A list of 4 elements: `dataTraining` which contains the actual data used for training, `clustering` which contains the clusters associated to each cell, `treeCuts` which contains the number of clusters for different hierarchical tree cutting and `hSequence` which contains the different thresholds used to cut the hierarchical tree.
#'
#' @importFrom foreach %do%
#'
#' @export

initialClusteringTraining = function(data = NULL, parametersToAnalyze = NULL, subsetDownsampled = 25000, clusterMinPercentage = 0.5)
{
  a = NULL

  dataTraining = data[data$state == "sampled", ]

  if (subsetDownsampled > nrow(dataTraining))
  {
    subsetDownsampled = nrow(dataTraining)
  }

  rowsToKeep_dataTraining = sample(as.numeric(rownames(dataTraining)), subsetDownsampled)

  dataTraining = dataTraining[as.numeric(rownames(dataTraining)) %in% rowsToKeep_dataTraining, ]

  dataTraining_restricted = dataTraining[, parametersToAnalyze]

  dataTraining_distanceMatrix = stats::dist(dataTraining_restricted, method = "euclidean")
  dataTraining_clustering = stats::hclust(dataTraining_distanceMatrix, method = "complete")

  heightsSequence = seq(1, floor(max(dataTraining_clustering$height)), 0.05)
  treeCuttingResults = data.frame(matrix(NA, ncol = 3, nrow = length(heightsSequence)))
  colnames(treeCuttingResults) = c("h", "clustersNb", "clustersNb.SupMin")

  foreach::foreach(a = 1:length(heightsSequence)) %do%
  {
    currentHeight = heightsSequence[a]
    currentHeight_clustering = stats::cutree(dataTraining_clustering, h = currentHeight)
    currentHeight_clustersNb = length(table(currentHeight_clustering))

    clustersSupMin_ID = which((table(currentHeight_clustering)/subsetDownsampled) * 100 >= clusterMinPercentage)
    currentHeight_clustersNbSupMin = length(clustersSupMin_ID)

    treeCuttingResults[a, ] = c(currentHeight, currentHeight_clustersNb, currentHeight_clustersNbSupMin)
  }

  treeCuttingResults$ID = rownames(treeCuttingResults)

  treeCuttingResults$clustersPercentage.SupMin = treeCuttingResults$clustersNb.SupMin/treeCuttingResults$clustersNb

  return(list(dataTraining = dataTraining, clustering = dataTraining_clustering, treeCuts = treeCuttingResults, hSequence = heightsSequence))

  gc()
}

#' Display an interactive plot for choosing an optimal cutoff value
#'
#' This function allows to display an interactive `ggplotly` decision plot to help the user to choose the optimal cutoff for the hierarchical tree cutting.
#'
#' @param initialClustering An object generated using the `initialClusteringTraining()` function. Defaults to `NULL`.
#'
#' @return A `ggplotly` plot. Please note that after the cutoff selection is performed, the graph should be close with `dev.off()` call.
#'
#' @export

plotInitialClusteringTraining = function(initialClustering = NULL)
{
  clustersNb.SupMin = NULL
  clustersPercentage.SupMin = NULL
  clustersNb = NULL
  ID = NULL
  h = NULL

  treeCuts = initialClustering$treeCuts

  p = ggplot2::ggplot(treeCuts, ggplot2::aes(clustersNb.SupMin, clustersPercentage.SupMin, text = paste("clustersNb: ", clustersNb, "\nh: ", h, "\nID: ", ID, sep = ""))) +
    ggplot2::geom_point() +
    ggplot2::theme_bw()

  plotly::ggplotly(p)
}

#' Perform final cell clustering training
#'
#' This function allows to apply the selected cutoff (for hierarchical clustering tree cutting) to the training dataset.
#'
#' @param initialClusteringData An object generated using the `initialClusteringTraining()` function. Defaults to `NULL`.
#'
#' @param clusterMinPercentage An integer defining the minimum frequency of a cluster to be kept as is. Values represent percentages and range from `0` to `100`. Defaults to `0.5`.
#'
#' @param cutoff An integer defining the actual cutoff to use. Defaults to `NULL`.
#'
#' @return A list of 2 elements: `dataTraining` which contains the actual data used for training and `cutoff` which contains the cutoff value used for the tree cutting.
#'
#' @export

finalClusteringTraining = function(initialClusteringData = NULL, clusterMinPercentage = NULL, cutoff = NULL)
{
  hSequence = initialClusteringData$hSequence
  clustering = initialClusteringData$clustering
  dataTraining = initialClusteringData$dataTraining

  associatedHeight = hSequence[cutoff]

  dataTraining_clustering_cut = stats::cutree(clustering, h = associatedHeight)

  clustersToDelete = which((table(dataTraining_clustering_cut)/nrow(dataTraining)) * 100 < clusterMinPercentage)

  dataTraining$cluster = dataTraining_clustering_cut

  if (length(clustersToDelete) > 0)
  {
    cellsToDeclusterize = which(dataTraining$cluster %in% clustersToDelete)
    dataTraining = dataTraining[-cellsToDeclusterize, ]
  }

  return(list(dataTraining = dataTraining, cutoff = cutoff))
}

#' Export plot showing the cutoff on the previous decision plot
#'
#' This function allows to export the previous decision plot on which is displayed the chosen cutoff for hierarchical clustering tree cutting.
#'
#' @param initialClustering An object generated using the `initialClusteringTraining()` function. Defaults to `NULL`.
#'
#' @param finalClustering An object generated using the `finalClusteringTraining()` function. Defaults to `NULL`.
#'
#' @param datasetFolder A string defining in which folder to save the results. It can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to `full`.
#'
#' @return A `ggplotly` plot. Please note that after the cutoff selection is performed, the graph should be close with `dev.off()` call.
#'
#' @export

plotFinalClusteringTraining = function(initialClustering = NULL, finalClustering = NULL, datasetFolder = "full")
{
  clustersNb.SupMin = NULL
  clustersPercentage.SupMin = NULL
  clustersNb = NULL
  ID = NULL
  h = NULL

  treeCuts = initialClustering$treeCuts

  treeCuts_cutoff = treeCuts[treeCuts$ID == finalClustering$cutoff, ]

  p = ggplot2::ggplot(treeCuts, ggplot2::aes(clustersNb.SupMin, clustersPercentage.SupMin, text = paste("clustersNb: ", clustersNb, "\nh: ", h, "\nID: ", ID, sep = ""))) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::geom_point(data = treeCuts_cutoff, colour = "red") +
    ggplot2::geom_text(data = treeCuts_cutoff, label = paste("ID = ", treeCuts_cutoff$ID, "\nclustersNb = ", treeCuts_cutoff$clustersNb, "\nclustersNb.SupMin = ", treeCuts_cutoff$clustersNb.SupMin, "\nclustersPercentages.SupMin = ", treeCuts_cutoff$clustersPercentage.SupMin, sep = ""), hjust = 1, vjust = 0.5, color = "red", size = 2.5)

  grDevices::pdf(file.path("output", "7_Clustering", datasetFolder, paste("clusterNb.SupMin_VS_clustersPercentage.Supmin_training.pdf", sep = "")))
  print(p)
  grDevices::dev.off()
}

#' Export UMAP plot with clusters determined in training dataset
#'
#' This function allows to export a PDF file showing the UMAP parameters of the involved training dataset and their colored determined clusters.
#'
#' @param dataTrainingClustered A data frame containing the clustered training dataset. Defaults to `NULL`.
#'
#' @param datasetFolder A string defining in which folder to save the results. It can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to `full`.
#'
#' @return Generated PDF plot is saved in the `output > 7_Clustering` directory.
#'
#' @export

plotUMAP_projectionTraining = function(dataTrainingClustered = NULL, datasetFolder = "full")
{
  UMAP_1 = NULL
  UMAP_2 = NULL
  cluster = NULL

  p = ggplot2::ggplot(dataTrainingClustered, ggplot2::aes(x = UMAP_2, y = UMAP_1, color = as.character(cluster))) +
    ggplot2::geom_point(size = 0.75, alpha = 0.5, show.legend = FALSE)

  grDevices::pdf(file.path("output", "7_Clustering", datasetFolder, paste("clustersProjectionUMAP_training.pdf", sep = "")))
  print(p)
  grDevices::dev.off()
}

#' Apply the training clustering on the whole dataset
#'
#' This function allows to apply the previously generated clustering on the remaining dataset not used in the training (also called validation dataset).
#'
#' @param dataTraining A data frame containing the training dataset. Defaults to `NULL`.
#'
#' @param dataValidation A data frame containing the validation dataset. Defaults to `NULL`.
#'
#' @param parametersToUse A vector defining the parameters to use for the model application. Defaults to `NULL`.
#'
#' @param coresNumber An integer defining the number of cores to use to apply the clustering model. Defaults to `1`.
#'
#' @param chunksMaxSize An integer defining the maximum cell number of a chunk to be processed. Defaults to `100000`.
#'
#' @return A data frame containing the clustered cells of the `dataValidation` dataset.
#'
#' @importFrom foreach %do%
#'
#' @export

applyClusterModel = function(dataTraining = dataTraining, dataValidation = dataValidation, parametersToUse = NULL, coresNumber = 1, chunksMaxSize = 1e+05)
{
  a = NULL
  b = NULL

  chunksNb = ceiling(nrow(dataValidation)/chunksMaxSize)

  chunksLimits = list()
  options(scipen = 999)

  foreach::foreach(a = 1:chunksNb) %do%
  {
    currentChunk_lowerLimit = ((a - 1) * chunksMaxSize) + 1

    if (a == chunksNb)
    {
      currentchunk_upperLimit = nrow(dataValidation)
    } else
    {
      currentchunk_upperLimit = a * chunksMaxSize
    }

    chunksLimits[[a]] = c(currentChunk_lowerLimit, currentchunk_upperLimit)
  }

  cl = parallel::makeCluster(coresNumber, type = "PSOCK")
  doSNOW::registerDoSNOW(cl)

  pb = utils::txtProgressBar(min = 0, max = length(chunksLimits), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)

  totalValidationClusters = foreach::foreach(b = 1:length(chunksLimits), .packages = c("foreach", "tcltk", "class"), .options.snow = opts) %dopar%
  {
    lowLimitChunk = chunksLimits[[b]][1]
    highLimitChunk = chunksLimits[[b]][2]

    knnClustTemp = class::knn(train = dataTraining[, parametersToUse], test = dataValidation[c(lowLimitChunk:highLimitChunk), parametersToUse], k = 10, cl = factor(dataTraining$cluster), prob = FALSE)

    currentChunkOutData = data.frame(cluster = knnClustTemp, rowID = as.numeric(rownames(dataValidation[c(lowLimitChunk:highLimitChunk), ])))
    return(currentChunkOutData)
  }

  close(pb)

  parallel::stopCluster(cl)

  dataValidation_clusters = do.call(rbind, totalValidationClusters)

  sort(unique(dataTraining$cluster)) == sort(unique(dataValidation_clusters$cluster))

  dataValidation$cluster = dataValidation_clusters$cluster
  dataValidation$clusteringGroup = "validation"

  dataTraining$clusteringGroup = "training"

  totalData = rbind(dataTraining, dataValidation)

  totalData$cluster = as.numeric(totalData$cluster)

  return(totalData)
}


#' Open the desired UMAP dataset
#'
#' This function allows to open the desired UMAP dataset for parameter threshold generation.
#'
#' @param datasetToUse A string defining which dataset to open. It can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to `full`.
#'
#' @return A data frame containing the selected dataset to use.
#'
#' @export

openUMAPData = function(datasetToUse)
{
  if (datasetToUse == "training")
  {
    data = readRDS(file.path("rds", "5_UMAP_training.rds"))
  } else if (datasetToUse == "validation")
  {
    data = readRDS(file.path("rds", "5_UMAP_validation.rds"))
  } else
  {
    datasetToUse = "full"
    data = readRDS(file.path("rds", "5_UMAP_full.rds"))
  }

  return(data)
}

#' Plot density of a given parameter on an interactive plot
#'
#' This function allows to display the density of a desired parameter on an interactive `ggplotly` plot to choose the best threshold to use for negative and positive cells separation.
#'
#' @param data A string defining which dataset to use. It should be the output of the `openUMAPData()` function. Defaults to `NULL`.
#'
#' @param parameter A string defining which parameter to plot Defaults to `NULL`.
#'
#' @param displayedCells A string defining the maximum number of cells to be plotted. Defaults to `100000`.
#'
#' @return A `ggplotly` plot. Please note that after the threshold selection is performed, the graph should be close with `dev.off()` call.
#'
#' @export

determineParameterThreshold = function(data = NULL, parameter = NULL, displayedCells = 1e+05)
{
  sampledDataThresholds = data[sample(1:nrow(data), displayedCells), ]
  p = ggplot2::ggplot(sampledDataThresholds, ggplot2::aes(x = get(parameter))) +
    ggplot2::geom_density(fill = "#69b3a2", color = "#e9ecef", alpha = 0.8)

  plotly::ggplotly(p)
}

#' Collapse the phenotypically close cell clusters
#'
#' This function allows to collapse the phenotypically close cell clusters in order to reduce their number and ease data interpretation. Results are saved in a new `rds` file.
#'
#' @param data A data frame containing the dataset to use. Defaults to `NULL`.
#'
#' @param parametersToUse A character vector defining the parameters to analyze in the subsequent dataset. Defaults to `NULL`.
#'
#' @param parametersToUseThresholds A data frame defining the thresholds for each parameter Defaults to `NULL`.
#'
#' @param metricUsed A string defining the metric to be used. It can be either `mean` or `median`. Defaults to `median`.
#'
#' @param datasetFolder A string defining in which folder to save the results. It can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to `full`.
#'
#' @param customClusteringCut A numeric value defining which `h` value to use for the dendrogram tree cutting. This determines how strong the cluster collapsing will be (higher value means stronger collapsing, thus resulting in fewer clusters at the end). Defaults to `0`, which corresponds to no cut.
#'
#' @return Generated rds file is saved to `rds` directory.
#'
#' @importFrom foreach %do%
#'
#' @export

collapseCloseClusters = function(data = NULL, parametersToUse = NULL, parametersToUseThresholds = NULL, metricUsed = "median", datasetFolder = "full", customClusteringCut = 0)
{
  a = NULL
  b = NULL

  clustersToBrowse = unique(data$cluster)

  clustersBinary = data.frame(matrix(NA, nrow = length(clustersToBrowse), ncol = length(parametersToUse)))

  foreach::foreach(a = 1:length(clustersToBrowse)) %do%
  {
    currentClusterID = clustersToBrowse[a]

    currentClusterAssociatedData = data[data$cluster == currentClusterID, ]
    currentClusterAssociatedData = currentClusterAssociatedData[, parametersToUse]

    if (metricUsed == "mean")
    {
      currentCluster_colMetric = colMeans(currentClusterAssociatedData, na.rm = TRUE)
    } else if (metricUsed == "median")
    {
      currentCluster_colMetric = array(matrixStats::colMedians(as.matrix(currentClusterAssociatedData), na.rm = TRUE), dimnames = list(colnames(currentClusterAssociatedData)))
    } else
    {
      currentCluster_colMetric = colMeans(currentClusterAssociatedData, na.rm = TRUE)
    }

    foreach::foreach(b = 1:length(parametersToUse)) %do%
    {
      currentParameter = parametersToUse[b]

      currentParameterThreshold = as.numeric(parametersToUseThresholds[parametersToUseThresholds$parameter == currentParameter, "threshold"])

      if (currentCluster_colMetric[currentParameter] < currentParameterThreshold)
      {
        currentCluster_colMetric[currentParameter] = 0
      } else
      {
        currentCluster_colMetric[currentParameter] = 1
      }
    }

    clustersBinary[a, ] = currentCluster_colMetric
    rownames(clustersBinary)[a] = paste("C", currentClusterID, sep = "")
  }

  colnames(clustersBinary) = parametersToUse

  distance = stats::dist(clustersBinary, method = "euclidean")
  clustering = stats::hclust(distance, method = "complete")

  clustering_cut = stats::cutree(clustering, h = customClusteringCut)

  clustersToPoolUnique = as.numeric(which(table(clustering_cut) > 1))

  maxOriginalClustersNb = max(unique(data$cluster))

  foreach::foreach(c = 1:length(clustersToPoolUnique)) %do%
  {
    currentClusterID = clustersToPoolUnique[c]
    clusterToPoolDetails = as.numeric(gsub("C", "", names(which(clustering_cut == currentClusterID))))
    currentClusterNewName = maxOriginalClustersNb + c
    data[data$cluster %in% clusterToPoolDetails, "cluster"] = currentClusterNewName
  }

  grDevices::pdf(file.path("output", "7_Clustering", datasetFolder, paste("clustersCollapsing_dendrogramForTreeCutting.pdf", sep = "")), paper = "a4r", width = 20)

  plot(clustering, cex = 0.5)
  abline(h = customClusteringCut, lwd = 2, col = "red")

  grDevices::dev.off()

  print(paste("Before collapsing: ", length(clustersToBrowse), " clusters", sep = ""))
  print(paste("After collapsing: ", length(unique(data$cluster)), " clusters", sep = ""))

  saveRDS(data, file.path("rds", paste("7_Clustered_", datasetFolder, ".rds", sep = "")))

  return(data)
}

#' Open the desired clustered dataset
#'
#' This function allows to open the desired clustered dataset for subsequent analysis.
#'
#' @param datasetToUse A string defining which dataset to open. It can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to `full`.
#'
#' @return A data frame containing the selected dataset to use.
#'
#' @export

openClusteredFullDataCollapsed = function(datasetToUse)
{
  if (datasetToUse == "training")
  {
    data = readRDS(file.path("rds", "7_Clustered_training.rds"))
  } else if (datasetToUse == "validation")
  {
    data = readRDS(file.path("rds", "7_Clustered_validation.rds"))
  } else if (datasetToUse == "full")
  {
    data = readRDS(file.path("rds", "7_Clustered_full.rds"))
  } else
  {
    data = readRDS(file.path("rds", "7_Clustered_full.rds"))
  }

  return(data)
}

#' Export final clusters statistics and plots
#'
#' This function allows to export statistics and plots for the final determined clusters.
#'
#' @param data A data frame containing the clustered dataset. It should be the output of `openClusteredFullDataCollapsed()` function. Defaults to `NULL`.
#'
#' @param folder A string defining the directory in which the cluster plots will be saved. Defaults to `NULL`.
#'
#' @param parametersToUse A character vector defining the parameters to use in the exported statistics. Defaults to `NULL`.
#'
#' @param coresNumber An integer defining the number of processor cores to be used for data exportation. Defaults to `2`.
#'
#' @param prefix A string defining the prefix to adjunct to each cluster name in each plot exported. Defaults to `NULL`.
#'
#' @param maxCellsPerPlot An integer defining the maximum cell number displayed per plot. Defaults to `1000`.
#'
#' @param metricUsed A string defining the metric to use for statistics computation. It can be either `mean` or `median`. Defaults to `median`.
#'
#' @param datasetFolder A string defining in which folder to save the results. It can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to `full`.
#'
#' @return Generated PDF plots are saved to `output > 7_Clustering > datasetFolder > folder` directory. Generated text files are saved to `output > 7_Clustering > datasetFolder` directory.
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#'
#' @export

exportClustersStatsAndPlots = function(data = NULL, folder = NULL, parametersToUse = NULL, coresNumber = 2, prefix = NULL, maxCellsPerPlot = 1000, metricUsed = "median", datasetFolder = "full")
{
  a = NULL
  UMAP_1 = NULL
  UMAP_2 = NULL
  b = NULL
  e = NULL
  z = NULL

  uniqueClusters = unique(data$cluster)

  minX = min(data$UMAP_2)
  maxX = max(data$UMAP_2)

  minY = min(data$UMAP_1)
  maxY = max(data$UMAP_1)

  if (dir.exists(file.path("output", "7_Clustering", datasetFolder, folder)))
  {
    unlink(file.path("output", "7_Clustering", datasetFolder, folder), recursive = TRUE)
  }

  dir.create(file.path("output", "7_Clustering", datasetFolder, folder))

  cl = parallel::makeCluster(coresNumber, type = "PSOCK")
  doSNOW::registerDoSNOW(cl)

  pb = utils::txtProgressBar(min = 0, max = length(uniqueClusters), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)

  uniqueGroups = unique(data$group)

  out = foreach::foreach(a = 1:length(uniqueClusters), .packages = c("foreach", "ggplot2", "tcltk", "matrixStats"), .options.snow = opts) %dopar%
  {
    currentCluster = uniqueClusters[a]

    dir.create(file.path("output", "7_Clustering", datasetFolder, folder, paste(prefix, "_C", currentCluster, sep = "")))

    currentCluster_data = data[data$cluster == currentCluster, ]

    uniquePatients = unique(paste("Group-", data$group, "_Sample-", data$name, sep = ""))

    if (nrow(currentCluster_data) > maxCellsPerPlot)
    {
      currentCluster_data_plotOverall = currentCluster_data[sample(rownames(currentCluster_data), maxCellsPerPlot), ]
    } else
    {
      currentCluster_data_plotOverall = currentCluster_data
    }

    grDevices::pdf(file.path("output", "7_Clustering", datasetFolder, folder, paste(prefix, "_C", currentCluster, sep = ""), paste("Cluster-", paste(prefix, "_C", currentCluster, sep = ""), "_Overall.pdf", sep = "")))

    plot = ggplot2::ggplot(currentCluster_data_plotOverall, ggplot2::aes(x = UMAP_2, y = UMAP_1)) +
      ggplot2::geom_point(alpha = 0.25) +
      ggplot2::xlim(minX, maxX) +
      ggplot2::ylim(minY, maxY)
    print(plot)

    grDevices::dev.off()

    if(nrow(currentCluster_data[currentCluster_data$clusteringGroup == "training", ]) > 0)
    {

      currentCluster_data_plot_training = currentCluster_data[currentCluster_data$clusteringGroup == "training", ]

      if (nrow(currentCluster_data_plot_training) > maxCellsPerPlot)
      {
        currentCluster_data_plot_training = currentCluster_data_plot_training[sample(rownames(currentCluster_data_plot_training), maxCellsPerPlot), ]
      }

      grDevices::pdf(file.path("output", "7_Clustering", datasetFolder, folder, paste(prefix, "_C", currentCluster, sep = ""), paste("Cluster-", paste(prefix, "_C", currentCluster, sep = ""), "_Training.pdf", sep = "")))

      plot = ggplot2::ggplot(currentCluster_data_plot_training, ggplot2::aes(x = UMAP_2, y = UMAP_1)) + ggplot2::geom_point(alpha = 0.25) + ggplot2::xlim(minX, maxX) + ggplot2::ylim(minY, maxY)
      print(plot)

      grDevices::dev.off()

    }


    if(nrow(currentCluster_data[currentCluster_data$clusteringGroup == "validation", ]) > 0)
    {
      currentCluster_data_plot_validation = currentCluster_data[currentCluster_data$clusteringGroup == "validation", ]

      if (nrow(currentCluster_data_plot_validation) > maxCellsPerPlot)
      {
        currentCluster_data_plot_validation = currentCluster_data_plot_validation[sample(rownames(currentCluster_data_plot_validation), maxCellsPerPlot), ]
      }



      grDevices::pdf(file.path("output", "7_Clustering", datasetFolder, folder, paste(prefix, "_C", currentCluster, sep = ""), paste("Cluster-", paste(prefix, "_C", currentCluster, sep = ""), "_Validation.pdf", sep = "")))

      plot = ggplot2::ggplot(currentCluster_data_plot_validation, ggplot2::aes(x = UMAP_2, y = UMAP_1)) + ggplot2::geom_point(alpha = 0.25) + ggplot2::xlim(minX, maxX) + ggplot2::ylim(minY, maxY)
      print(plot)

      grDevices::dev.off()

    }




    currentCluster_percentages = NULL
    currentCluster_samples = NULL

    foreach::foreach(b = 1:length(uniqueGroups)) %do%
    {
      currentGroup = uniqueGroups[b]

      currentCluster_data_currentGroup = currentCluster_data[currentCluster_data$group == currentGroup, ]

      currentGroup_data = data[data$group == currentGroup, ]

      uniqueGroupPatientsID = grep(paste("Group-", currentGroup, "_", sep = ""), uniquePatients)
      currentCluster_currentGroup_data_patients = gsub(paste("Group-", currentGroup, "_Sample-", sep = ""), "", uniquePatients[uniqueGroupPatientsID])

      currentCluster_data_currentGroup_plot = currentCluster_data_currentGroup

      if (nrow(currentCluster_data_currentGroup_plot) > maxCellsPerPlot)
      {
        currentCluster_data_currentGroup_plot = currentCluster_data_currentGroup_plot[sample(rownames(currentCluster_data_currentGroup_plot), maxCellsPerPlot), ]
      }

      foreach::foreach(c = 1:length(currentCluster_currentGroup_data_patients)) %do%
      {
        currentPatient = currentCluster_currentGroup_data_patients[c]

        currentCluster_data_currentGroup_currentPatient = currentCluster_data_currentGroup[currentCluster_data_currentGroup$name == currentPatient, ]
        currentPatient_data = currentGroup_data[currentGroup_data$name == currentPatient, ]

        currentCluster_patientPercentage = (nrow(currentCluster_data_currentGroup_currentPatient)/nrow(currentPatient_data)) * 100

        currentCluster_percentages = c(currentCluster_percentages, currentCluster_patientPercentage)
        currentCluster_samples = c(currentCluster_samples, paste("Group-", currentGroup, "_Sample-", currentPatient, sep = ""))

        currentCluster_data_currentGroup_currentPatient_plot = currentCluster_data_currentGroup_currentPatient

        if (nrow(currentCluster_data_currentGroup_currentPatient_plot) > maxCellsPerPlot)
        {
          currentCluster_data_currentGroup_currentPatient_plot = currentCluster_data_currentGroup_currentPatient_plot[sample(rownames(currentCluster_data_currentGroup_currentPatient_plot), maxCellsPerPlot), ]
        }
      }
    }

    currentCluster_samplesPercentages = data.frame(sample = currentCluster_samples, percentage = currentCluster_percentages)

    colnames(currentCluster_samplesPercentages)[2] = paste(prefix, "_C", currentCluster, sep = "")

    currentCluster_samplesPercentages = currentCluster_samplesPercentages[order(currentCluster_samplesPercentages$sample), ]

    rownames(currentCluster_samplesPercentages) = currentCluster_samplesPercentages$sample

    currentCluster_samplesPercentages$sample = NULL

    currentCluster_data_parametersOnly = currentCluster_data[, parametersToUse]

    if (metricUsed == "mean")
    {
      currentCluster_data_colMetric = data.frame(colMeans(currentCluster_data_parametersOnly))
      colnames(currentCluster_data_colMetric) = paste(prefix, "_C", currentCluster, sep = "")
    } else if (metricUsed == "median")
    {
      currentCluster_data_colMetric = data.frame(matrixStats::colMedians(as.matrix(currentCluster_data_parametersOnly)))
      rownames(currentCluster_data_colMetric) = colnames(currentCluster_data_parametersOnly)
      colnames(currentCluster_data_colMetric) = paste(prefix, "_C", currentCluster, sep = "")
    } else
    {
      currentCluster_data_colMetric = data.frame(colMeans(currentCluster_data_parametersOnly))
      colnames(currentCluster_data_colMetric) = paste(prefix, "_C", currentCluster, sep = "")
    }

    return(list(clustersPercentages = currentCluster_samplesPercentages, clustersPhenotypes = currentCluster_data_colMetric))
  }

  parallel::stopCluster(cl)

  foreach::foreach(e = 1:length(out)) %do%
  {
    currentData = out[[e]]

    currentData_clustersPhenotypes = currentData$clustersPhenotypes
    currentData_clustersPercentages = currentData$clustersPercentages

    if (e == 1)
    {
      totalClustersPhenotypes = currentData_clustersPhenotypes
      totalClustersPercentages = currentData_clustersPercentages

    } else
    {
      totalClustersPhenotypes = cbind(totalClustersPhenotypes, currentData_clustersPhenotypes)
      totalClustersPercentages = cbind(totalClustersPercentages, currentData_clustersPercentages)
    }
  }

  totalClustersPhenotypes_export = cbind(rownames(totalClustersPhenotypes), totalClustersPhenotypes)
  colnames(totalClustersPhenotypes_export)[1] = "Clusters"

  utils::write.table(totalClustersPhenotypes_export, file.path("output", "7_Clustering", datasetFolder, paste(prefix, "_clustersPhenotypes.txt", sep = "")), sep = "\t", row.names = FALSE, col.names = TRUE)

  totalClustersPercentages = cbind(rownames(totalClustersPercentages), totalClustersPercentages)
  colnames(totalClustersPercentages)[1] = "Sample"

  totalClustersPercentages = cbind(rownames(totalClustersPercentages), totalClustersPercentages)
  colnames(totalClustersPercentages)[1] = "SampleCorrected"

  utils::write.table(totalClustersPercentages, file.path("output", "7_Clustering", datasetFolder, paste(prefix, "_clustersPercentages.txt", sep = "")), sep = "\t", row.names = FALSE, col.names = TRUE)

  uniqueGroups = unique(gsub("Group-(.+)_Sample-(.+)", "\\1", totalClustersPercentages$Sample))

  newGroupedTotalData = data.frame(matrix(0, nrow = length(uniqueGroups), ncol = ncol(totalClustersPercentages) - 1), stringsAsFactors = FALSE)

  foreach::foreach(z = 1:length(uniqueGroups)) %do%
  {
    currentGroup = uniqueGroups[z]

    associatedSamples = grep(paste("Group-", currentGroup, "_Sample-(.+)", sep = ""), totalClustersPercentages$Sample)
    currentData = totalClustersPercentages[associatedSamples, ]

    if (metricUsed == "mean")
    {
      currentDataGrouped = colMeans(currentData[, -c(1:2)])

    } else if (metricUsed == "median")
    {
      currentDataGrouped = matrixStats::colMedians(as.matrix(currentData[, -c(1:2)]))

    } else
    {
      currentDataGrouped = colMeans(currentData[, -c(1:2)])
    }

    newGroupedTotalData[z, ] = c(currentGroup, currentDataGrouped)
  }

  colnames(newGroupedTotalData) = colnames(totalClustersPercentages)[-1]
  colnames(newGroupedTotalData)[1] = "Group"

  utils::write.table(newGroupedTotalData, file.path("output", "7_Clustering", datasetFolder, paste(prefix, "_clustersPercentages_perGroup.txt", sep = "")), sep = "\t", row.names = FALSE, col.names = TRUE)
}

#' Export UMAP plot with final clusters
#'
#' This function allows to export a PDF file showing the UMAP parameters of the involved dataset and their colored determined clusters.
#'
#' @param data A data frame containing the clustered dataset. It should be the output of `openClusteredFullDataCollapsed()` function. Defaults to `NULL`.
#'
#' @param displayedCells An integer defining the maximum number of cells to display on the plot. Defaults to `100000`.

#' @param datasetFolder A string defining in which folder to save the results. It can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to `full`.
#'
#' @return Generated PDF plot is saved to `output > 7_Clustering` directory.
#'
#' @export

plotUMAP_projectionFinalClusters = function(data = NULL, displayedCells = 1e+05, datasetFolder = "full")
{
  UMAP_1 = NULL
  UMAP_2 = NULL
  cluster = NULL

  if (dir.exists(file.path("output", "7_Clustering", datasetFolder)) == FALSE)
  {
    dir.create(file.path("output", "7_Clustering", datasetFolder))
  }


  clusterizedFullData_collapsed_plot = data[sample(1:nrow(data), displayedCells), ]
  p = ggplot2::ggplot(clusterizedFullData_collapsed_plot, ggplot2::aes(x = UMAP_2, y = UMAP_1, color = as.character(cluster))) +
    ggplot2::geom_point(size = 0.75, alpha = 0.5, show.legend = FALSE)

  grDevices::pdf(file.path("output", "7_Clustering", datasetFolder, paste("clustersProjectionUMAP_allData_final.pdf", sep = "")))
  print(p)
  grDevices::dev.off()
}

#' Perform FlowSOM clustering
#'
#' This function allows to perform FlowSOM clustering algorithm on the dataset and parameters of interest.
#'
#' @param data A data frame containing the dataset to use. Defaults to `NULL`.
#'
#' @param parametersToUse A vector defining the parameters to use for the clustering. Defaults to `NULL`.
#'
#' @param seed An integer defining the seed to use to obtain reproducible results. Defaults to `NULL`.
#'
#' @param maxMeta An integer defining the maximum number of clusters to generate. Please note that FlowSOM cannot algorithm cannot exceed 90 clusters maximum. However, FlowSOM will surely output way less than this number, depending on the dataset of interest. Defaults to `90`.
#'
#' @param datasetFolder A string defining in which folder to save the results. It can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to `full`.
#'
#' @return Generated rds file is saved to `rds` directory.
#'
#' @export

FlowSOM_clustering = function(data = NULL, parametersToUse = NULL, seed = NULL, maxMeta = 90, datasetFolder = "full")
{

  FlowSOM_output = FlowSOM::FlowSOM(input = as.matrix(data[, parametersToUse]), compensate = FALSE, transform = FALSE, scale = FALSE, silent = FALSE, seed = seed, maxMeta = maxMeta)


  FlowSOM_clustering = as.numeric(FlowSOM::GetMetaclusters(FlowSOM_output))

  print(paste("FlowSOM identified ", length(unique(FlowSOM_clustering)), " clusters within the dataset.", sep = ""))

  data$cluster = FlowSOM_clustering
  data$clusteringGroup = "validation"

  saveRDS(data, file.path("rds", paste("7_Clustered_", datasetFolder, ".rds", sep = "")))

  return(data)
}

#' Perform PhenoGraph clustering
#'
#' This function allows to perform PhenoGraph clustering algorithm on the dataset and parameters of interest.
#'
#' @param data A data frame containing the dataset to use. Defaults to `NULL`.
#'
#' @param parametersToUse A vector defining the parameters to use for the clustering. Defaults to `NULL`.
#'
#' @param k An integer defining the maximum number of clusters to generate. PhenoGraph will surely output way less than this number, depending on the dataset of interest. Defaults to `100`.
#'
#' @param coresNumber An integer defining the number of cores to use to apply the clustering model. Defaults to `1`.
#'
#' @param datasetFolder A string defining in which folder to save the results. It can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to `full`.
#'
#' @return Generated rds file is saved to `rds` directory.
#'
#' @export

FastPhenoGraph_clustering = function(data = NULL, parametersToUse = NULL, k = 100, coresNumber = 1, datasetFolder = "full")
{


  FastPhenoGraph_output = FastPG::fastCluster(as.matrix(data[, parametersToUse]), k, num_threads = coresNumber, verbose = TRUE)

  FastPhenoGraph_clustering = as.numeric(FastPhenoGraph_output$communities)

  data$cluster = FastPhenoGraph_clustering
  data$clusteringGroup = "validation"


  saveRDS(data, file.path("rds", paste("7_Clustered_", datasetFolder, ".rds", sep = "")))

  return(data)

}

#' Export information about clusters phenotype
#'
#' This function allows to export heatmaps and text files about clusters phenotype. It also summarizes the phenotypes by using binary mapping of the data using the pre-determined thresholds for each parameter.
#'
#' @param prefix A string defining the prefix of the clusters phenotypes file to use. Defaults to `NULL`.
#'
#' @param thresholds (Optional) A data frame defining the thresholds for each parameter. If set to `NULL`, then binary heatmaps will not be generated. Thresholds values for each parameter are typically obtained when the user chooses the hierarchical clustering + k-nearest neighbors approach for cell clustering. however, the user can still provide these thresholds by following the subsection called "Determine binary thresholds" detailed in the tutorial. Defaults to `NULL`.
#'
#' @param metricUsed A string defining the metric to use for statistics computation. It can be either `mean` or `median`. Defaults to `median`.
#'
#' @param datasetFolder A string defining in which folder to save the results. It can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to `full`.
#'
#' @return Generated PDF files are saved to `output > 7_Clustering > datasetFolder` directory. Generated text files are saved to `output > 7_Clustering > datasetFolder` directory.
#'
#' @export

clustersPhenotypesHeatmap = function(prefix = NULL, thresholds = NULL, metricUsed = "median", datasetFolder = "full")
{
  a = NULL

  totalClustersPhenotypes = utils::read.csv(file = file.path("output", "7_Clustering", datasetFolder, paste(prefix, "_clustersPhenotypes.txt", sep = "")), sep = "\t")

  rownames(totalClustersPhenotypes) = totalClustersPhenotypes$Clusters
  totalClustersPhenotypes$Clusters = NULL

  markersNames = colnames(totalClustersPhenotypes)
  clustersNames = rownames(totalClustersPhenotypes)

  totalClustersPhenotypes = apply(totalClustersPhenotypes, 2, as.numeric)
  colnames(totalClustersPhenotypes) = markersNames
  rownames(totalClustersPhenotypes) = clustersNames

  clustersToKeep = clustersNames

  totalClustersPhenotypes = totalClustersPhenotypes[clustersToKeep, ]

  featuresToDelete = NULL

  featuresToDeleteID = which(colnames(totalClustersPhenotypes) %in% featuresToDelete)

  if (length(featuresToDeleteID) > 0)
  {
    totalClustersPhenotypes = totalClustersPhenotypes[, -featuresToDeleteID]
  }

  totalClustersPhenotypes_unscaled = totalClustersPhenotypes

  totalClustersPhenotypes = scale(totalClustersPhenotypes)

  min = min(totalClustersPhenotypes, na.rm = TRUE)
  max = max(totalClustersPhenotypes, na.rm = TRUE)

  if (metricUsed == "mean")
  {
    middle = mean(totalClustersPhenotypes, na.rm = TRUE)
  } else if (metricUsed == "median")
  {
    middle = stats::median(totalClustersPhenotypes, na.rm = TRUE)
  } else
  {
    middle = mean(totalClustersPhenotypes, na.rm = TRUE)
  }


	separator = min(abs(min-middle), abs(max-middle))/2

  colors = sort(c(seq(min, (middle - separator) - 0.001, length = 25), seq(middle - separator, middle + separator, length = 50), seq((middle + separator) + 0.001, max, length = 25)))

  customPalette = (grDevices::colorRampPalette(c("lightblue", "blue", "black", "yellow", "orange")))(n = 99)

  grDevices::pdf(file = file.path("output", "7_Clustering", datasetFolder, paste("heatmap_clustersPhenotypes_rowScaled.pdf", sep = "")), bg = "transparent", width = 20, height = 10, paper = "a4r")

  gplots::heatmap.2(totalClustersPhenotypes, trace = "none", scale = "row", key = TRUE, keysize = 1, density.info = "none", col = customPalette, cexRow = 0.35, cexCol = 0.35, margins = c(6, 8), dendrogram = "both")

  grDevices::dev.off()

  grDevices::pdf(file = file.path("output", "7_Clustering", datasetFolder, paste("heatmap_clustersPhenotypes_rowUnscaled.pdf", sep = "")), bg = "transparent", width = 20, height = 10, paper = "a4r")

  gplots::heatmap.2(totalClustersPhenotypes, trace = "none", scale = "none", key = TRUE, keysize = 1, density.info = "none", col = customPalette, breaks = colors, cexRow = 0.35, cexCol = 0.35, margins = c(6, 8), dendrogram = "both")

  grDevices::dev.off()

  totalClustersPhenotypes_exported = totalClustersPhenotypes

  totalClustersPhenotypes_exported = cbind(rownames(totalClustersPhenotypes_exported), totalClustersPhenotypes_exported)
  colnames(totalClustersPhenotypes_exported)[1] = "Markers"

  utils::write.table(totalClustersPhenotypes_exported, file = file.path("output", "7_Clustering", datasetFolder, paste("heatmap_clustersPhenotypes_data.txt", sep = "")), quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)

  if(is.null(thresholds) == FALSE)
  {

    totalClustersPhenotypes_binary = totalClustersPhenotypes_unscaled

    foreach::foreach(a = 1:nrow(totalClustersPhenotypes_binary)) %do%
      {
        currentRowData = totalClustersPhenotypes_binary[a, ]
        currentThreshold = thresholds[a]

        currentRowData[currentRowData < currentThreshold] = 0
        currentRowData[currentRowData >= currentThreshold] = 1

        totalClustersPhenotypes_binary[a, ] = currentRowData
      }

    min = min(totalClustersPhenotypes_binary, na.rm = TRUE)
    max = max(totalClustersPhenotypes_binary, na.rm = TRUE)

    if (metricUsed == "mean")
    {
      middle = mean(totalClustersPhenotypes_binary, na.rm = TRUE)
    } else if (metricUsed == "median")
    {
      middle = stats::median(totalClustersPhenotypes_binary, na.rm = TRUE)
    } else
    {
      middle = mean(totalClustersPhenotypes_binary, na.rm = TRUE)
    }

    grDevices::pdf(file = file.path("output", "7_Clustering", datasetFolder, paste("heatmap_clustersPhenotypes_binary.pdf", sep = "")), bg = "transparent", width = 20, height = 10, paper = "a4r")

    gplots::heatmap.2(totalClustersPhenotypes_binary, trace = "none", scale = "none", key = TRUE, keysize = 1, density.info = "none", col = c("blue", "yellow"), cexRow = 0.35, cexCol = 0.35, margins = c(6, 8), dendrogram = "both")

    grDevices::dev.off()

    totalClustersPhenotypes_binary_exported = totalClustersPhenotypes_binary

    totalClustersPhenotypes_binary_exported = cbind(rownames(totalClustersPhenotypes_binary_exported), totalClustersPhenotypes_binary_exported)
    colnames(totalClustersPhenotypes_binary_exported)[1] = "Markers"

    utils::write.table(totalClustersPhenotypes_binary_exported, file = file.path("output", "7_Clustering", datasetFolder, paste("heatmap_clustersPhenotypes_binary_data.txt", sep = "")), quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)
  }


}

#' Export information about clusters abundance (group- or sample-wise)
#'
#' This function allows to export heatmaps about clusters abundance per group or per sample.
#'
#' @param prefix A string defining the prefix of the clusters abundance file to use. Defaults to `NULL`.
#'
#' @param metricUsed A string defining the metric to use for statistics computation. It can be either `mean` or `median`. Defaults to `median`.
#'
#' @param datasetFolder A string defining in which folder to save the results. It can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to `full`.
#'
#' @param mode A string defining the mode to use for heatmap generation. It can be either `group` or `sample`. Defaults to `group`.
#'
#' @return Generated PDF files are saved to `output > 7_Clustering > datasetFolder`.
#'
#' @export

clustersPercentagesHeatmap = function(prefix = NULL, metricUsed = "median", datasetFolder = "full", mode = "group")
{
  d = NULL

  if(mode == "group")
  {
    totalClustersPercentages = utils::read.csv(file = file.path("output", "7_Clustering", datasetFolder, paste(prefix, "_clustersPercentages_perGroup.txt", sep = "")), sep = "\t")


    rownames(totalClustersPercentages) = totalClustersPercentages$Group
    totalClustersPercentages$Group = NULL

    clustersNames = colnames(totalClustersPercentages)
    groupsNames = rownames(totalClustersPercentages)

    totalClustersPercentages = apply(totalClustersPercentages, 2, as.numeric)
    colnames(totalClustersPercentages) = clustersNames
    rownames(totalClustersPercentages) = groupsNames
  } else if(mode == "sample")
  {
    totalClustersPercentages = utils::read.csv(file = file.path("output", "7_Clustering", datasetFolder, paste(prefix, "_clustersPercentages.txt", sep = "")), sep = "\t")

    rownames(totalClustersPercentages) = totalClustersPercentages$Sample
    totalClustersPercentages$Sample = NULL
    totalClustersPercentages$SampleCorrected = NULL

    clustersNames = colnames(totalClustersPercentages)
    sampleNames = rownames(totalClustersPercentages)

    totalClustersPercentages = apply(totalClustersPercentages, 2, as.numeric)
    colnames(totalClustersPercentages) = clustersNames
    rownames(totalClustersPercentages) = sampleNames

  } else
  {
    totalClustersPercentages = utils::read.csv(file = file.path("output", "7_Clustering", datasetFolder, paste(prefix, "_clustersPercentages_perGroup.txt", sep = "")), sep = "\t")


    rownames(totalClustersPercentages) = totalClustersPercentages$Group
    totalClustersPercentages$Group = NULL

    clustersNames = colnames(totalClustersPercentages)
    groupsNames = rownames(totalClustersPercentages)

    totalClustersPercentages = apply(totalClustersPercentages, 2, as.numeric)
    colnames(totalClustersPercentages) = clustersNames
    rownames(totalClustersPercentages) = groupsNames
  }



  totalClustersPercentages = scale(totalClustersPercentages)

  if(length(as.numeric(which(colSums(apply(totalClustersPercentages, 2, is.na)) > 0))) > 0)
  {
    colsToReplace = as.numeric(which(colSums(apply(totalClustersPercentages, 2, is.na)) > 0))
    totalClustersPercentages[, colsToReplace] = 0
  }

  min = min(totalClustersPercentages, na.rm = TRUE)
  max = max(totalClustersPercentages, na.rm = TRUE)

  if (metricUsed == "mean")
  {
    middle = mean(totalClustersPercentages, na.rm = TRUE)
  } else if (metricUsed == "median")
  {
    middle = stats::median(totalClustersPercentages, na.rm = TRUE)
  } else
  {
    middle = mean(totalClustersPercentages, na.rm = TRUE)
  }

    	separator = min(abs(min-middle), abs(max-middle))/2


  colors = sort(c(seq(min, (middle - separator) - 0.001, length = 25), seq(middle - separator, middle + separator, length = 50), seq((middle + separator) + 0.001, max, length = 25)))



  customPalette = (grDevices::colorRampPalette(c("lightblue", "blue", "black", "yellow", "orange")))(n = 99)

  if(mode == "group")
  {

    grDevices::pdf(file = file.path("output", "7_Clustering", datasetFolder, paste("heatmap_clustersPercentages_perGroup_columnScaled.pdf", sep = "")), bg = "transparent", width = 20, height = 10, paper = "a4r")

    gplots::heatmap.2(totalClustersPercentages, trace = "none", scale = "none", key = TRUE, keysize = 1, density.info = "none", col = customPalette, breaks = colors, cexRow = 0.35, cexCol = 0.35, margins = c(20, 8), dendrogram = "both")


  } else if(mode == "sample")
  {

    groupsPerSample = gsub("Group-(.+)_Sample-(.+)", "\\1", rownames(totalClustersPercentages))

    uniqueGroups = unique(groupsPerSample)

    groupsPalette = grDevices::rainbow(length(uniqueGroups))

    foreach::foreach(d = 1:length(uniqueGroups)) %do%
      {

        currentGroup = uniqueGroups[d]

        groupsPerSample[groupsPerSample == currentGroup] = groupsPalette[d]
      }

    grDevices::pdf(file = file.path("output", "7_Clustering", datasetFolder, paste("heatmap_clustersPercentages_perSample_columnScaled.pdf", sep = "")), bg = "transparent", width = 20, height = 10, paper = "a4r")

    gplots::heatmap.2(totalClustersPercentages, trace = "none", scale = "none", key = TRUE, keysize = 1, density.info = "none", col = customPalette, breaks = colors, cexRow = 0.35, cexCol = 0.35, margins = c(18, 14), dendrogram = "both", RowSideColors = groupsPerSample)

  } else
  {
    grDevices::pdf(file = file.path("output", "7_Clustering", datasetFolder, paste("heatmap_clustersPercentages_perGroup_columnScaled.pdf", sep = "")), bg = "transparent", width = 20, height = 10, paper = "a4r")

    gplots::heatmap.2(totalClustersPercentages, trace = "none", scale = "none", key = TRUE, keysize = 1, density.info = "none", col = customPalette, breaks = colors, cexRow = 0.35, cexCol = 0.35, margins = c(20, 8), dendrogram = "both")

  }



  grDevices::dev.off()

  # grDevices::pdf(file = file.path("output", "7_Clustering", datasetFolder, paste("heatmap_clustersPercentages_columnUnscaled.pdf", sep = "")), bg = "transparent", width = 20, height = 10, paper = "a4r")
  #
  # gplots::heatmap.2(totalClustersPercentages, trace = "none", scale = "none", key = FALSE, keysize = 0.5, col = customPalette, cexRow = 0.8, cexCol = 0.8, margins = c(20, 8), dendrogram = "both")
  #
  # grDevices::dev.off()
}
