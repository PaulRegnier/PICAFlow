#' Open the full, training or validation downsampled dataset
#'
#' This function allows to open one of the three possible downsampled datasets: the full dataset, the training subdataset or the validation subdataset.
#'
#' @param datasetToUse A string defining which dataset to open. Can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. Defaults to `full`.
#'
#' @return Returns the associated `rds` data.
#'
#' @export

openDownsampledData = function(datasetToUse = "full")
{
  if (datasetToUse == "training")
  {
    data = readRDS(file.path("rds", "4_Downsampled_training.rds"))
  } else if (datasetToUse == "validation")
  {
    data = readRDS(file.path("rds", "4_Downsampled_validation.rds"))
  } else
  {
    datasetToUse = "full"
    data = readRDS(file.path("rds", "4_Downsampled_full.rds"))
  }

  return(data)
}

#' Determine optimal UMAP parameters
#'
#' This function computes UMAP dimensionality reduction using different combinations of `n_neighbors` and `min_dist` parameters to allow the user to choose the combination that best suits the current dataset. Usually, this function is called on a small subdataset, typically the cells that are tagged as downsampled.
#'
#' @param data A data frame on which to perform the analysis. It must originate from `openDownsampledData()` method call. Defaults to `NULL`.
#'
#' @param nNeighborsToTest A numeric vector defining the different `n_neighbors` to test. Defaults to `NULL`.
#'
#' @param minDistToTest A numeric vector defining the different `min_dist` to test. Defaults to `NULL`.
#'
#' @param downsample An integer defining the maximum number of cells to use for the analysis. Defaults to `50000`.
#'
#' @param datasetFolder A string defining which dataset to open. Can be either `full` for the full downsampled dataset, `training` for the training downsampled subdataset or `validation` for the validation downsampled subdataset. The value must match the origin of the `data` argument. Defaults to `full`.
#'
#' @return Generated PDF file is saved to `output > 5_UMAP` directory.
#'
#' @importFrom foreach %do%
#'
#' @export

determineOptimalUMAPParameters = function(data = NULL, nNeighborsToTest = NULL, minDistToTest = NULL, downsample = 50000, datasetFolder = "full")
{
  a = NULL
  b = NULL
  UMAP1 = NULL
  UMAP2 = NULL

  if (downsample > 0)
  {
    if (downsample > nrow(data))
    {
      downsample = nrow(data)
    }

    data = data[sample(nrow(data), downsample), ]
  }

  coresNumber = parallel::detectCores() - 1
  set.seed(42)

  plotsList = list()

  print("Applying UMAP algorithm with various parameters...")

  foreach::foreach(a = 1:length(nNeighborsToTest)) %do%
  {
    currentNNeighborsToTest = nNeighborsToTest[a]

    foreach::foreach(b = 1:length(minDistToTest)) %do%
    {
      currentMinDistToTest = minDistToTest[b]

      currentParameters_umap = uwot::umap(data[, -c((ncol(data) - 2):ncol(data))], init = "random", scale = FALSE, verbose = TRUE, n_threads = coresNumber, n_components = 2, min_dist = currentMinDistToTest, n_neighbors = currentNNeighborsToTest, fast_sgd = TRUE)

      currentParameters_umapEmbeddings = as.data.frame(currentParameters_umap, stringsAsFactors = FALSE)
      colnames(currentParameters_umapEmbeddings) = c("UMAP1", "UMAP2")

      currentPlot = ggplot2::ggplot(currentParameters_umapEmbeddings, ggplot2::aes(x = UMAP1, y = UMAP2)) +
        ggplot2::geom_bin2d(bins = 250) +
        ggplot2::scale_fill_continuous(type = "viridis") +
        ggplot2::ggtitle(paste("n_neighbors = ", currentNNeighborsToTest, "\nmin_dist = ", currentMinDistToTest, sep = "")) +
        ggplot2::theme_bw() +
        ggplot2::theme(text = ggplot2::element_text(size = 6), legend.position = "none", plot.title = ggplot2::element_text(size = 6, hjust = 0.5), axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 6), axis.text.y = ggplot2::element_text(size = 6), panel.grid.major = ggplot2::element_line(linewidth = 0.5), panel.grid.minor = ggplot2::element_line(linewidth = 0.5))

      plotsList[[length(plotsList) + 1]] = cowplot::as_grob(currentPlot)
    }

    print(paste("Done: ", round(a/length(nNeighborsToTest) * 100, 0), "%", sep = ""))
  }

  outputPlots = cowplot::plot_grid(plotlist = plotsList, ncol = length(minDistToTest), nrow = length(nNeighborsToTest))

  rm(currentParameters_umap)
  rm(plotsList)
  gc()

  if (dir.exists(file.path("output", "5_UMAP", datasetFolder)) == TRUE)
  {
    unlink(file.path("output", "5_UMAP", datasetFolder), recursive = TRUE)
  }

  dir.create(file.path("output", "5_UMAP", datasetFolder))

  grDevices::pdf(file = file.path("output", "5_UMAP", datasetFolder, "facets_nNeighbors-vs-minDist.pdf"), width = length(minDistToTest) * 2, height = length(nNeighborsToTest) * 2)
  print(outputPlots)
  grDevices::dev.off()
}

#' Compute UMAP on the whole downsampled dataset
#'
#' This function computes UMAP dimensionality reduction on the whole downsampled dataset and returns the associated UMAP embeddings as well as the model used.
#'
#' @param data A data frame on which to perform the analysis. It must originate from `openDownsampledData()` method call. Defaults to `NULL`.
#'
#' @param n_threads An integer defining the number of cores to use for UMAP computation. Be careful, this function is expensive in terms of computational resources, so do not increase this value too much. Defaults to `1`.
#'
#' @param min_dist A numeric vector defining the different `min_dist` to test. Defaults to `NULL`.
#'
#' @param n_neighbors A numeric vector defining the different `n_neighbors` to test. Defaults to `NULL`.
#'
#' @return A full list of UMAP information computed for the provided data.
#'
#' @export

UMAP_downsampledDataset = function(data = NULL, n_threads = 1, min_dist = NULL, n_neighbors = NULL)
{
  data_umap_out = uwot::umap(data[, -c((ncol(data)-2):ncol(data))], init = "random", scale = FALSE, verbose = TRUE, n_threads = n_threads, ret_model = TRUE, n_components = 2, min_dist = min_dist, n_neighbors = n_neighbors)

  return(data_umap_out)
}

#' Apply UMAP model on the remaining cells
#'
#' This function allows to apply the UMAP model generated on downsampled cells to all the other cells. Usually, this function is called on a large subdataset, typically the cells that were not tagged as downsampled.
#'
#' @param data A data frame on which to perform the analysis. It must originate from `openDownsampledData()` method call. Defaults to `NULL`.
#'
#' @param model An object obtained after a call to umap() method from uwot package with ret_model argument set to TRUE. Defaults to `NULL`.
#'
#' @param chunksMaxSize An integer defining the maximum cell number of each chunk for UMAP computation. Defaults to `1000000`.
#'
#' @param n_threads An integer defining the number of cores to use for UMAP computation. Be careful, this function is expensive in terms of computational resources, so do not increase this value too much. Defaults to `1`.
#'
#' @return Returns the data and associated UMAP parameters.
#'
#' @export

UMAPFlowset = function(data = NULL, model = NULL, chunksMaxSize = 1e+06, n_threads = 1)
{
  a = NULL
  b = NULL

  chunksNb = ceiling(nrow(data)/chunksMaxSize)
  chunksLimits = list()
  options(scipen = 999)

  foreach::foreach(a = 1:chunksNb) %do%
  {
    currentChunk_lowerLimit = ((a - 1) * chunksMaxSize) + 1

    if (a == chunksNb)
    {
      currentchunk_upperLimit = nrow(data)
    } else
    {
      currentchunk_upperLimit = a * chunksMaxSize
    }

    chunksLimits[[a]] = c(currentChunk_lowerLimit, currentchunk_upperLimit)
  }

  dataNotSampled_UMAPIndexes_list = list()

  print("Applying previous UMAP model to not sampled data...")

  foreach::foreach(b = 1:length(chunksLimits)) %do%
  {
    currentChunk_lowerLimit = chunksLimits[[b]][1]
    currentchunk_upperLimit = chunksLimits[[b]][2]

    currentChunkData = data[currentChunk_lowerLimit:currentchunk_upperLimit, ]

    currentChunk_dataNotSampled_umap_out = uwot::umap_transform(currentChunkData[, -c((ncol(currentChunkData) - 2):ncol(currentChunkData))], model, verbose = TRUE, n_threads = n_threads)

    dataNotSampled_UMAPIndexes_list[[b]] = data.frame(UMAP_1 = currentChunk_dataNotSampled_umap_out[, 1], UMAP_2 = currentChunk_dataNotSampled_umap_out[, 2])
    rm(currentChunkData)
    rm(currentChunk_dataNotSampled_umap_out)
    gc()

    paste("Done: ", round(b/length(chunksLimits) * 100, 0), "%", sep = "")
  }

  dataNotSampled_UMAPIndexes_list = do.call(rbind.data.frame, dataNotSampled_UMAPIndexes_list)
  data$UMAP_1 = dataNotSampled_UMAPIndexes_list[, 1]
  data$UMAP_2 = dataNotSampled_UMAPIndexes_list[, 2]

  rm(model)
  rm(dataNotSampled_UMAPIndexes_list)
  gc()

  return(data)
}
