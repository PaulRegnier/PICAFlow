#' Check new version
#'
#' This function checks whether there is a new available version for `PICAFlow` in the GitHub repository. If this is the case, users are notified.
#' @param ... Useful to transfer other arguments to the function.
#'
#' @export

checkNewVersion = function(...)
{
  latestRelease = gh::gh("GET /repos/{owner}/{repo}/releases", owner = "PaulRegnier", repo = "PICAFlow")
  latestRelease = as.character(latestRelease[[1]][["tag_name"]])
  installedRelease = as.character(utils::packageVersion("PICAFlow"))

  if(installedRelease != latestRelease)
  {
    print(paste("[WARNING] The version of PICAFlow which is currently installed on this computer ('", installedRelease, "') does not match with the latest version available on GitHub ('", latestRelease, "'). We kindly invite you to update to the latest version to avoid any troubles related to bugs and/or missing features. You can launch the update with the following command: 'devtools::install_github('PaulRegnier/PICAFlow', force = TRUE)'. Please see the tutorial for more information. Thank you for your understanding.", sep = ""))
  } else
  {
    print(paste("[NOTE] The version of PICAFlow which is currently installed on this computer (", installedRelease, ") is up-to-date.", sep = ""))
  }
}

#' Setup the working directory
#'
#' This function setups the working directory. It does not take any input argument but assumes that the `workingDirectory` variable is set to the path of interest using the base R `setwd()` function. From this path, everything will be wiped and `input`, `output` and `rds` directories will be created. The `output`` directory will also contain several other specific directories that will be used in the subsequent analyses.
#'
#' @export

setupWorkingDirectory = function()
{
  workingDirectory = getwd()

  if(dir.exists(file.path(workingDirectory)) == TRUE && length(grep("(input)|(output)|(rds)", list.dirs(file.path(workingDirectory)))) > 0)
  {
    print(paste("The '", file.path(workingDirectory), "' working directory already exists and seems to contain folders ususally created by PICAFlow. To avoid any unwanted deletion of data, please check first that you already backed up its content to another place then delete the content of the '", file.path(workingDirectory), "' directory before running this function again.", sep = ""))

  } else
  {
    unlink(file.path(workingDirectory, "input"), recursive = TRUE)
    dir.create(file.path(workingDirectory, "input"))

    unlink(file.path(workingDirectory, "output"), recursive = TRUE)
    dir.create(file.path(workingDirectory, "output"))
    dir.create(file.path(workingDirectory, "output", "1_Transformation"))
    dir.create(file.path(workingDirectory, "output", "2_Normalization"))
    dir.create(file.path(workingDirectory, "output", "3_Gating"))
    dir.create(file.path(workingDirectory, "output", "4_Downsampling"))
    dir.create(file.path(workingDirectory, "output", "5_UMAP"))
    dir.create(file.path(workingDirectory, "output", "6_FCS"))
    dir.create(file.path(workingDirectory, "output", "7_Clustering"))
    dir.create(file.path(workingDirectory, "output", "8_Analysis"))

    unlink(file.path(workingDirectory, "rds"), recursive = TRUE)
    dir.create(file.path(workingDirectory, "rds"))

  }


}

#' Convert `fcs` files to `rds` files
#'
#' This function converts every `fcs` file in the `input` directory to a `flowFrame` object (from `flowCore` package) encapsulated into a `rds` file. It eventually renames the parameters used in the dataset if desired. This step helps to decrease the overall computing time and complexity of several next steps.
#'
#' @param conversionTable A tabular-delimited text file containing the appropriate information to convert one or several channels to other ones. This is used when the overall staining mix presents the same specificities but slightly different fluorophores depending on the samples. This typically occurs when the cytometer configuration is modified or when the cytometer is different from one batch to another. Defaults to `NULL`. We typically advise users to launch this function the first time with this argument set to `NULL`, in order to clearly see where the differences for names/descriptions are located. Then, users can construct and specify a conversion table with the `conversionTable` parameter which allows to resolve the differences observed in the dataset.
#'
#' The `conversionTable` table follows a pre-defined format: 4 columns in any order (`from_desc`, `to_desc`, `from_name` and `to_name`) and any given number of line, each line referring to a specific matching to be treated. For instance, if a line has the values `from_desc = CXCR5 B610-ECD-A`, `to_desc = CXCR5 ECD-A`, `from_name = FL2-A` and `to_name = FL11-A`, it means that any occurrence of a parameter (in any `rds` file) named `FL2-A` which also matches the description `CXCR5 B610-ECD-A` will see its values respectively replaced with `FL11-A` and `CXCR5 ECD-A`.
#'
#' Please note that the renaming both affects the parameters of each file AND the parameters in each self-contained compensation matrix.
#'
#' If needed, it is also possible to ask for the deletion of one or several channels in the dataset by specifying the word `DELETE` in the `to_desc` and `to_name` columns of the `conversionTable` table. In this case, it (or they) will simply be deleted from the samples which match with the associated pair of `from_desc` and `from_name` parameters used in the table.
#'
#' @return Generated `rds` files are saved to `rds` directory. The function also returns a matrix of all unique parameters (once correctly renamed) used in the dataset.
#'
#' @importFrom foreach %dopar%
#'
#' @export

convertToRDS = function(conversionTable = NULL)
{
  a = NULL
  b = NULL

  workingDirectory = getwd()

  fcs.dir = file.path(workingDirectory, "input")
  filesToOpen = dir(fcs.dir, full.names = TRUE)

  if (parallel::detectCores() == 1)
  {
    coresNumber = 1
  } else
  {
    coresNumber = parallel::detectCores() - 1
  }

  cl = parallel::makeCluster(coresNumber, type = "PSOCK")
  doSNOW::registerDoSNOW(cl)

  pb = utils::txtProgressBar(min = 0, max = length(filesToOpen), style = 3)

  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)

  totalParametersInfos = foreach::foreach(a = filesToOpen, .packages = c("foreach", "flowCore", "tcltk"), .combine = "c", .options.snow = opts) %dopar%
    {
      currentData = flowCore::read.FCS(a, transformation = FALSE, truncate_max_range = FALSE, ignore.text.offset = TRUE)
      currentFilename = gsub("(.+)/(.+)\\.fcs", "\\2", a)

      if (is.null(conversionTable) == FALSE)
      {

        foreach::foreach(t = 1:ncol(conversionTable)) %do%
          {
            currentColumnData = conversionTable[, t]

            currentColumnData_isNA = which(is.na(currentColumnData))

            if(length(currentColumnData_isNA) > 0)
            {
              conversionTable[currentColumnData_isNA, t] = "NA"

            }

          }


        currentFileParameterDescriptions = as.character(as.vector(currentData@parameters@data[, "desc"]))
        currentFileParameterDescriptions[is.na(currentFileParameterDescriptions)] = "NA"

        currentFileParameterNames = as.character(as.vector(currentData@parameters@data[, "name"]))
        currentFileParameterNames[is.na(currentFileParameterNames)] = "NA"

        isCompensationMatricesPresent = try(flowStats::spillover(currentData), silent = TRUE)

        if(class(isCompensationMatricesPresent) != "try-error")
        {
          compensationMatricesSlot = as.numeric(which(lengths(flowStats::spillover(currentData)) > 0))
          compensationMatricesSlotName = names(flowStats::spillover(currentData))[compensationMatricesSlot]

        } else
        {
          compensationMatricesSlotName = NULL

        }

        parametersToRemoveID = NULL

        foreach::foreach(b = 1:length(currentFileParameterDescriptions)) %do% {
          currentDescriptionFrom = currentFileParameterDescriptions[b]
          currentDescriptionName = currentFileParameterNames[b]

          currentNameFromMatching = which((conversionTable$from_desc == currentDescriptionFrom) & conversionTable$from_name == currentDescriptionName)

          if (length(currentNameFromMatching) > 0)
          {
            currentDescriptionTo = conversionTable[currentNameFromMatching, "to_desc"]
            currentNameTo = conversionTable[currentNameFromMatching, "to_name"]
            currentNameFrom = conversionTable[currentNameFromMatching, "from_name"]

            rowToReplaceID = as.numeric(which(currentFileParameterDescriptions == currentDescriptionFrom & currentFileParameterNames == currentDescriptionName))

            if(currentDescriptionTo == "DELETE" & currentNameTo == "DELETE")
            {
              parametersToRemoveID = c(parametersToRemoveID, rowToReplaceID)
              # currentData = currentData[, -rowToReplaceID]
            } else
            {
              currentData@parameters@data[rowToReplaceID, "desc"] = paste(currentDescriptionTo, "_replaced", sep = "")
              currentData@parameters@data[rowToReplaceID, "name"] = paste(currentNameTo, "_replaced", sep = "")

              colToReplaceID = as.numeric(which(colnames(currentData@exprs) == currentNameFrom))
              colnames(currentData@exprs)[colToReplaceID] = paste(currentNameTo, "_replaced", sep = "")


              if (length(compensationMatricesSlotName) > 0)
              {
                matchingCompensationNameID = which(colnames(currentData@description[compensationMatricesSlotName][[1]]) == currentNameFrom)

                colnames(currentData@description[compensationMatricesSlotName][[1]])[matchingCompensationNameID] = paste(currentNameTo, "_replaced", sep = "")
              }

            }


          }
        }

        if(is.null(parametersToRemoveID) == FALSE & length(parametersToRemoveID) > 0)
        {
          currentData = currentData[, -parametersToRemoveID]
        }

        currentData@parameters@data[, "desc"] = gsub("_replaced", "", currentData@parameters@data[, "desc"])
        currentData@parameters@data[, "name"] = gsub("_replaced", "", currentData@parameters@data[, "name"])
        colnames(currentData@exprs) = gsub("_replaced", "", colnames(currentData@exprs))


        if (length(compensationMatricesSlotName) > 0)
        {
          colnames(currentData@description[compensationMatricesSlotName][[1]]) = gsub("_replaced", "", colnames(currentData@description[compensationMatricesSlotName][[1]]))

        }

      }

      currentFileDescriptions = as.character(as.vector(currentData@parameters@data[, "desc"]))

      currentFileDescriptions[is.na(currentFileDescriptions)] = "NA"

      currentSampleOrder = order(as.vector(currentData@parameters@data[, "name"]))

      currentFileNames = as.vector(currentData@parameters@data[, "name"])[currentSampleOrder]
      currentFileDescriptions = currentFileDescriptions[currentSampleOrder]

      currentParametersInfos = list(currentFileNames, currentFileDescriptions)

      names(currentParametersInfos) = c(paste("names_", currentFilename, sep = ""), paste("descriptions_", currentFilename, sep = ""))

      saveRDS(currentData, file.path("rds", paste(currentFilename, ".rds", sep = "")))

      gc()

      return(currentParametersInfos)
    }

  close(pb)

  parallel::stopCluster(cl)

  totalParametersNamesID = grep("names", names(totalParametersInfos))
  totalParametersNames = totalParametersInfos[totalParametersNamesID]


  # Generate a list of unique sets of parameters names

  n = NULL

  totalParametersNames_table = t(data.frame(totalParametersNames))

  totalParametersNames_tableDuplicatedIDs = as.numeric(which(duplicated(totalParametersNames_table)))

  if(length(totalParametersNames_tableDuplicatedIDs) > 0)
  {
    totalParametersNames_table = data.frame(totalParametersNames_table[-totalParametersNames_tableDuplicatedIDs, ])
  }

  totalParameterNames_synthetic = NULL

  if(ncol(totalParametersNames_table) == 1)
  {

    totalParameterNames_synthetic = totalParametersNames_table[, 1]
  } else
  {
    foreach::foreach(n = 1:ncol(totalParametersNames_table)) %do%
      {
        currentParameterNames = as.character(totalParametersNames_table[, n])

        if(length(currentParameterNames) > 1)
        {
          currentParameterNames_synthetic = paste("'", paste(unique(currentParameterNames), collapse = "' or '"), "'", sep = "")

        } else
        {
          currentParameterNames_synthetic = currentParameterNames

        }

        totalParameterNames_synthetic = c(totalParameterNames_synthetic, currentParameterNames_synthetic)

      }

  }




  totalParametersDescriptionsID = grep("descriptions", names(totalParametersInfos))
  totalParametersDescriptions = totalParametersInfos[totalParametersDescriptionsID]

  totalParametersDescriptions_table = t(data.frame(totalParametersDescriptions))

  totalParametersDescriptions_tableDuplicatedIDs = as.numeric(which(duplicated(totalParametersDescriptions_table)))

  if(length(totalParametersDescriptions_tableDuplicatedIDs) > 0)
  {
    totalParametersDescriptions_table = data.frame(totalParametersDescriptions_table[-totalParametersDescriptions_tableDuplicatedIDs, ])
  }


  totalParameterDescriptions_synthetic = NULL

  if(ncol(totalParametersDescriptions_table) == 1)
  {

    totalParameterDescriptions_synthetic = totalParametersDescriptions_table[, 1]
  } else
  {
    foreach::foreach(n = 1:ncol(totalParametersDescriptions_table)) %do%
      {
        currentParameterDescriptions = as.character(totalParametersDescriptions_table[, n])

        if(length(currentParameterDescriptions) > 1)
        {
          currentParameterDescriptions_synthetic = paste("'", paste(unique(currentParameterDescriptions), collapse = "' or '"), "'", sep = "")

        } else
        {
          currentParameterDescriptions_synthetic = currentParameterDescriptions

        }

        totalParameterDescriptions_synthetic = c(totalParameterDescriptions_synthetic, currentParameterDescriptions_synthetic)

      }

  }




  totalParametersData = data.frame(totalParameterNames_synthetic, totalParameterDescriptions_synthetic, stringsAsFactors = FALSE)
  totalParametersData = cbind(1:nrow(totalParametersData), totalParametersData)
  colnames(totalParametersData) = c("Parameter_ID", "Parameter_Name", "Parameter_Description")

  return(totalParametersData)
}

#' Export new `rds` files from a previously generated pooled version
#'
#' This function exports new `rds` files (one per sample) from a previously generated pooled version
#'
#' @param RDSFileToUse A vector of all parameters to plot. Defaults to `NULL`.
#'
#' @param coresNumber An integer defining the number of cores to use to export the files. It corresponds to the number of samples to be treated concomitantly. Defaults to `2`.
#'
#' @return Generated `rds` files are saved to `rds` directory. Please note that they will overwrite any preexisting `rds` file with the same name.
#'
#' @importFrom foreach %dopar%
#'
#' @export

exportRDSFilesFromPool = function(RDSFileToUse = NULL, coresNumber = 2)
{
  a = NULL

  flowSetToSplit = readRDS(file.path("rds", paste(RDSFileToUse, ".rds", sep = "")))

  samplesToExport = flowCore::sampleNames(flowSetToSplit)


  cl = parallel::makeCluster(coresNumber, type = "PSOCK")
  doSNOW::registerDoSNOW(cl)

  pb = utils::txtProgressBar(min = 0, max = length(samplesToExport), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)


  foreach::foreach(a = samplesToExport, .packages = c("foreach", "flowCore", "tcltk"), .options.snow = opts) %dopar%
    {

      flowSetToSplit = readRDS(file.path("rds", paste(RDSFileToUse, ".rds", sep = "")))
      currentFlowFrame = flowSetToSplit[[a]]
      currentFlowFrameName = a

      saveRDS(currentFlowFrame, file.path("rds", paste(currentFlowFrameName, ".rds", sep = "")))

      gc()

    }

  parallel::stopCluster(cl)

}


#' Subset `rds` files
#'
#' This function subsets every `rds` file with the specified parameters of interest. It can also rename each parameter in a more user-friendly way.
#'
#' @param parametersToKeep A vector of all parameters to keep in each `rds` file. Defaults to `NULL`.
#'
#' @param customNames A vector of the same length as `parametersToKeep` containing the user-friendly parameter names. Defaults to `NULL`.
#'
#' @return Generated `rds` files are saved to `rds` directory, overwriting the previous ones.
#'
#' @importFrom foreach %dopar%
#'
#' @export

subsetData = function(parametersToKeep = NULL, customNames = NULL)
{
  a = NULL

  filesToOpen = dir(file.path("rds"), full.names = TRUE)

  coresNumber = parallel::detectCores() - 1

  cl = parallel::makeCluster(coresNumber, type = "PSOCK")
  doSNOW::registerDoSNOW(cl)

  pb = utils::txtProgressBar(min = 0, max = length(filesToOpen), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)

  foreach::foreach(a = filesToOpen, .packages = c("foreach", "flowCore", "tcltk"), .options.snow = opts) %dopar%
  {
    currentData = readRDS(a)
    currentFilename = gsub("(.+)/(.+)\\.rds", "\\2", a)

    currentData = currentData[, parametersToKeep]

    currentData@parameters@data[, "desc"] = customNames


  # If the sample does not have any spillover matrices embedded, flowStats::spillover() will throw an error and crash

  currentDataCompensationExists = try(flowStats::spillover(currentData), silent = TRUE)

  if((class(currentDataCompensationExists) == "try-error") == FALSE)
  {

    compensationMatricesSlot = as.numeric(which(lengths(flowStats::spillover(currentData)) > 0))
    compensationMatricesSlotName = names(flowStats::spillover(currentData))[compensationMatricesSlot]

    # Removing useless parameters from the embedded compensation matrix

    matchingCompensationNameID = which(colnames(currentData@description[compensationMatricesSlotName][[1]]) %in% parametersToKeep)

    if (length(matchingCompensationNameID) > 0)
    {

      currentData@description[compensationMatricesSlotName][[1]] =   currentData@description[compensationMatricesSlotName][[1]][matchingCompensationNameID, matchingCompensationNameID]

      }
  }

    saveRDS(currentData, file.path("rds", paste(currentFilename, ".rds", sep = "")))

    rm(currentData)
    rm(currentFilename)
    gc()
  }

  parallel::stopCluster(cl)
}


#' Launch the R Shiny interactive application for edition of data compensations
#'
#' This function launches a R Shiny interactive application which helps to visually adjust data compensations.
#'
#' @param fs_shiny A flowSet generated using `mergeData()` function. Defaults to `NULL`.
#'
#' @param maxEventsNumber The maximum number of events to use for data displaying. Defaults to `10000`.
#'
#' @param options A list containing any desired parameter to pass to the Shiny application
#'
#' @importFrom foreach %do%
#'
#' @export

launchCompensationTuningShinyApp = function(fs_shiny = NULL, maxEventsNumber = 10000, options = list())
{
  o = NULL
  a = NULL
  b = NULL
  m = NULL

  samplesValuesList = as.list(seq(1, length(fs_shiny), 1))
  names(samplesValuesList) = rownames(Biobase::phenoData(fs_shiny))

  totalSampleLengths = NULL
  foreach::foreach(o = 1:length(samplesValuesList)) %do%
    {
      totalSampleLengths = c(totalSampleLengths, nrow(fs_shiny[[o]]@exprs))
    }

  pooledValuesShiny = list()

  foreach::foreach(o = 1:length(samplesValuesList)) %do%
    {
      names(samplesValuesList)[o] = paste(names(samplesValuesList)[o], " (#", nrow(fs_shiny[[o]]@exprs), ")", sep = "")

      pooledValuesShiny[[o]] = flowCore::exprs(fs_shiny[[o]])[sample(nrow(fs_shiny[[o]]), round(mean(totalSampleLengths)/length(totalSampleLengths))), ]
    }

  pooledValuesShiny = do.call(rbind.data.frame, pooledValuesShiny)
  pooledValuesShiny = as.matrix(pooledValuesShiny)

  fs_shiny = methods::rbind2(fs_shiny, fs_shiny[[1]])

  flowCore::exprs(fs_shiny[[(length(samplesValuesList) + 1)]]) = pooledValuesShiny

  flowWorkspace::sampleNames(fs_shiny)[(length(samplesValuesList) + 1)] = "All_Samples"

  samplesValuesList[paste("All_Samples (#", nrow(fs_shiny[[(length(samplesValuesList) + 1)]]@exprs), ")", sep = "")] = length(samplesValuesList) + 1

  parametersNamesList = as.list(seq(1, length(as.vector(colnames(fs_shiny[[1]]@exprs))), 1))
  names(parametersNamesList) = as.vector(colnames(fs_shiny[[1]]@exprs))

  foreach::foreach(o = 1:length(parametersNamesList)) %do%
    {
      parameterName = names(parametersNamesList)[o]

      if (parameterName %in% names(flowCore::markernames(fs_shiny)) == FALSE)
      {
        names(parametersNamesList)[o] = paste(parameterName, sep = "")
        parametersNamesList[o] = parameterName
      } else
      {
        associatedParameterNameId = which(names(flowCore::markernames(fs_shiny)) %in% parameterName)
        associatedParameterName = as.character(flowCore::markernames(fs_shiny)[associatedParameterNameId])
        names(parametersNamesList)[o] = paste(parameterName, " (", associatedParameterName, ")", sep = "")
        parametersNamesList[o] = parameterName
      }
    }





  # Define UI
  ui = shiny::pageWithSidebar(

    shiny::titlePanel(shiny::h1(shiny::strong("Tuning of compensation"), align = "center")),

    shiny::sidebarPanel(

      shiny::selectInput("sample", "Sample", samplesValuesList),
      shiny::selectInput("parameterX", "Parameter x", parametersNamesList),
      shiny::selectInput("parameterY", "Parameter y", parametersNamesList),
      shiny::sliderInput("parameterX_compensation", shiny::strong("Parameter x compensation"), min = -2, max = 2, value = 0, step = 0.01),
      shiny::sliderInput("parameterY_compensation", shiny::strong("Parameter y compensation"), min = -2, max = 2, value = 0, step = 0.01),


      shiny::actionButton('export_rds', 'Export rds'),
      shiny::actionButton('clear_all', 'Clear exports'),
      width = 4),

    shiny::mainPanel(
      shiny::plotOutput(outputId = "scatterPlot", height = "900px")
    )
  )

  # Define server logic
  server = function(input, output, session)
  {


    shiny::observe({

      if(file.exists(file.path("rds", "compensationMatrix.rds")) == TRUE)
      {

        compensationMatrix = readRDS(file.path("rds", "compensationMatrix.rds"))
      } else
      {
        # Instead of taking the compensation matrix from the first sample of the dataset, it is better to extract the first available one which is not empty
        totalCompensationMatrices = list()

        foreach::foreach(m = 1:length(fs_shiny)) %do%
          {
            currentFlowFrame = fs_shiny[[m]]

            compensationMatricesSlot = as.numeric(which(lengths(flowStats::spillover(currentFlowFrame)) > 0))
            compensationMatrix = flowStats::spillover(currentFlowFrame)[[compensationMatricesSlot]]
            rownames(compensationMatrix) = colnames(compensationMatrix)

            if(all(colSums(compensationMatrix) == rowSums(compensationMatrix)) == TRUE)
            {

              totalCompensationMatrices[[m]] = "empty"
            } else
            {
              totalCompensationMatrices[[m]] = "notEmpty"


            }
          }

        notEmptyCompensationMatricesIDs = which(totalCompensationMatrices == "notEmpty")
        referenceSampleIDForCompensationMatrix = notEmptyCompensationMatricesIDs[1]


        compensationMatricesSlot = as.numeric(which(lengths(flowStats::spillover(fs_shiny[[referenceSampleIDForCompensationMatrix]])) > 0))
        compensationMatrix = flowStats::spillover(fs_shiny[[referenceSampleIDForCompensationMatrix]])[[compensationMatricesSlot]]
        rownames(compensationMatrix) = colnames(compensationMatrix)

      }


      currentPairCompensations_1 = compensationMatrix[input$parameterX, input$parameterY]
      currentPairCompensations_2 = compensationMatrix[input$parameterY, input$parameterX]

      shiny::updateSliderInput(session, "parameterX_compensation", value = currentPairCompensations_1)

      shiny::updateSliderInput(session, "parameterY_compensation", value = currentPairCompensations_2)





    })

    output$scatterPlot = shiny::renderPlot({

      if(file.exists(file.path("rds", "compensationMatrix.rds")) == TRUE)
      {

        compensationMatrix = readRDS(file.path("rds", "compensationMatrix.rds"))
      } else
      {


        # Instead of taking the compensation matrix from the first sample of the dataset, it is better to extract the first available one which is not empty
        totalCompensationMatrices = list()

        foreach::foreach(m = 1:length(fs_shiny)) %do%
          {
            currentFlowFrame = fs_shiny[[m]]

            compensationMatricesSlot = as.numeric(which(lengths(flowStats::spillover(currentFlowFrame)) > 0))
            compensationMatrix = flowStats::spillover(currentFlowFrame)[[compensationMatricesSlot]]
            rownames(compensationMatrix) = colnames(compensationMatrix)

            if(all(colSums(compensationMatrix) == rowSums(compensationMatrix)) == TRUE)
            {

              totalCompensationMatrices[[m]] = "empty"
            } else
            {
              totalCompensationMatrices[[m]] = "notEmpty"


            }
          }

        notEmptyCompensationMatricesIDs = which(totalCompensationMatrices == "notEmpty")
        referenceSampleIDForCompensationMatrix = notEmptyCompensationMatricesIDs[1]


        compensationMatricesSlot = as.numeric(which(lengths(flowStats::spillover(fs_shiny[[referenceSampleIDForCompensationMatrix]])) > 0))
        compensationMatrix = flowStats::spillover(fs_shiny[[referenceSampleIDForCompensationMatrix]])[[compensationMatricesSlot]]
        rownames(compensationMatrix) = colnames(compensationMatrix)

      }

      compensationMatrix[input$parameterX, input$parameterY] = input$parameterX_compensation
      compensationMatrix[input$parameterY, input$parameterX] = input$parameterY_compensation


      data = fs_shiny[[as.numeric(input$sample)]]

      parametersListName = flowCore::markernames(fs_shiny)

      parameterXNameId = which(input$parameterX == names(flowCore::markernames(fs_shiny)))
      parameterYNameId = which(input$parameterY == names(flowCore::markernames(fs_shiny)))

      parameterXName = as.character(flowCore::markernames(fs_shiny)[parameterXNameId])
      parameterYName = as.character(flowCore::markernames(fs_shiny)[parameterYNameId])

      data = flowCore::compensate(data, compensationMatrix)

      saveRDS(compensationMatrix, file.path("rds", "compensationMatrix.rds"))


      resPlot = as.data.frame(data@exprs, stringsAsFactors = FALSE)[1:maxEventsNumber,]


      plotMainTitle = paste("Dot plot for parameter x = ", input$parameterX, " and parameter y = ", input$parameterY, sep = "")


      plot = ggplot2::ggplot(resPlot, ggplot2::aes(x = resPlot[, input$parameterX], y = resPlot[, input$parameterY])) +

        ggplot2::geom_bin_2d(bins = 256) +
        ggplot2::scale_fill_continuous(type = "viridis") +

        ggplot2::xlim(floor(min(resPlot[, input$parameterX])) , ceiling(max(resPlot[, input$parameterX]))) +
        ggplot2::ylim(floor(min(resPlot[, input$parameterY])) , ceiling(max(resPlot[, input$parameterY]))) +

        ggplot2::labs(title = plotMainTitle, x = paste("Intensity of fluorescence (", parameterXName, ")", sep = ""), y = paste("Intensity of fluorescence (", parameterYName, ")", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 14), axis.title = ggplot2::element_text(size = 16, face = "bold"), plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5))
      plot
    })




    shiny::observeEvent(input$export_rds, {
      exportedValues = shiny::reactiveValuesToList(input)




      if(file.exists(file.path("rds", "compensationMatrix.rds")) == TRUE)
      {

        compensationMatrix = readRDS(file.path("rds", "compensationMatrix.rds"))
      } else
      {


        # Instead of taking the compensation matrix from the first sample of the dataset, it is better to extract the first available one which is not empty
        totalCompensationMatrices = list()

        foreach::foreach(m = 1:length(fs_shiny)) %do%
          {
            currentFlowFrame = fs_shiny[[m]]

            compensationMatricesSlot = as.numeric(which(lengths(flowStats::spillover(currentFlowFrame)) > 0))
            compensationMatrix = flowStats::spillover(currentFlowFrame)[[compensationMatricesSlot]]
            rownames(compensationMatrix) = colnames(compensationMatrix)

            if(all(colSums(compensationMatrix) == rowSums(compensationMatrix)) == TRUE)
            {

              totalCompensationMatrices[[m]] = "empty"
            } else
            {
              totalCompensationMatrices[[m]] = "notEmpty"


            }
          }

        notEmptyCompensationMatricesIDs = which(totalCompensationMatrices == "notEmpty")
        referenceSampleIDForCompensationMatrix = notEmptyCompensationMatricesIDs[1]


        compensationMatricesSlot = as.numeric(which(lengths(flowStats::spillover(fs_shiny[[referenceSampleIDForCompensationMatrix]])) > 0))
        compensationMatrix = flowStats::spillover(fs_shiny[[referenceSampleIDForCompensationMatrix]])[[compensationMatricesSlot]]
        rownames(compensationMatrix) = colnames(compensationMatrix)


      }



      saveRDS(compensationMatrix, file.path("output", "1_Transformation", "compensationMatrix.rds"))



    })

    shiny::observeEvent(input$clear_all, {
      exportedValues = shiny::reactiveValuesToList(input)


      if (file.exists(file.path("output", "1_Transformation", "compensationMatrix.rds")))
      {
        unlink(file.path("output", "1_Transformation", "compensationMatrix.rds"))
      }

      if (file.exists(file.path("rds", "compensationMatrix.rds")))
      {
        unlink(file.path("rds", "compensationMatrix.rds"))
      }


    })
  }

  shiny::shinyApp(ui = ui, server = server, options = options)
}



#' Compensate `rds` files
#'
#' This function applies the self-contained compensation matrix of each `rds` file to every desired parameter.
#'
#' @param parametersToCompensate A vector of all parameters to compensate in each `rds` file. Defaults to `NULL`.
#'
#' @param useCustomCompensationMatrix A boolean defining if untouched embedded (if set to `FALSE`) or tuned (if set to `TRUE`) compensation matrix should be used for data compensation. Defaults to `FALSE`.
#'
#' @return Generated `rds` files are saved to `rds` directory, overwriting the previous ones.
#'
#' @importFrom foreach %dopar%
#'
#' @export

compensateData = function(parametersToCompensate = NULL, useCustomCompensationMatrix = FALSE)
{
  a = NULL

  filesToOpen = dir(file.path("rds"), full.names = TRUE)

  coresNumber = parallel::detectCores() - 1

  cl = parallel::makeCluster(coresNumber, type = "PSOCK")
  doSNOW::registerDoSNOW(cl)

  pb = utils::txtProgressBar(min = 0, max = length(filesToOpen), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)

  foreach::foreach(a = filesToOpen, .packages = c("foreach", "flowCore", "tcltk"), .options.snow = opts) %dopar%
  {
    # Debug only : a = rds.files[1:coresNumber]

    currentData = readRDS(a)
    currentFilename = gsub("(.+)/(.+)\\.rds", "\\2", a)

    if(useCustomCompensationMatrix == TRUE)
    {
    currentSampleComp = readRDS(file.path("output", "1_Transformation", "compensationMatrix.rds"))

    } else
    {
      compensationMatricesSlot = as.numeric(which(lengths(flowStats::spillover(currentData)) > 0))
      compensationMatrix = flowStats::spillover(currentData)[[compensationMatricesSlot]]
      parametersToKeepID = which(colnames(compensationMatrix) %in% parametersToCompensate)
      currentSampleComp = compensationMatrix[parametersToKeepID, parametersToKeepID]


    }


    currentData = flowCore::compensate(currentData, currentSampleComp)

    saveRDS(currentData, file.path("rds", paste(currentFilename, ".rds", sep = "")))

    rm(currentData)
    rm(currentSampleComp)
    rm(parametersToKeepID)
    rm(compensationMatrix)
    rm(compensationMatricesSlot)
    gc()
  }

  parallel::stopCluster(cl)
}

#' Get all channels information
#'
#' This function gets the names and descriptions of all the parameters contained in the first available `rds` file in the `rds` directory.
#'
#' @return A dataframe containing the ID, name and description of all the found parameters.
#'
#' @export

getAllChannelsInformation = function()
{
  listRDSfiles = list.files(file.path("rds"), pattern = ".rds")
  fileToOpen = listRDSfiles[1]

  fileContent = readRDS(file.path("rds", fileToOpen))

  if(length(fileContent) > 1)
  {
    fileContent = fileContent[[1]]

  }

  allParametersDescriptions = as.vector(flowCore::parameters(fileContent)$desc)
  allParametersNames = as.vector(flowCore::parameters(fileContent)$name)

  allParametersInformation = data.frame(matrix(0, nrow = length(allParametersDescriptions), ncol = 3))
  colnames(allParametersInformation) = c("Parameter_ID", "Parameter_Name", "Parameter_Description")

  allParametersInformation$Parameter_ID = 1:length(allParametersDescriptions)
  allParametersInformation$Parameter_Name = allParametersNames
  allParametersInformation$Parameter_Description = allParametersDescriptions

  return(allParametersInformation)
}

