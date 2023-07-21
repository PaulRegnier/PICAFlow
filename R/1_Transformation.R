#' Merge `rds` files into a single `rds` file
#'
#' This function merges all the `flowFrame` objects from each `rds` file into a single `flowSet` object encapsulated in a single R-specific `rds` file.
#'
#' @param suffix A character string eventually used to select some `rds` files to merge. Defaults to `NULL`.
#'
#' @return Generated `rds` file is saved to `rds` directory.
#'
#' @importFrom foreach %do%
#'
#' @export

mergeSamples = function(suffix = NULL)
{
  a = NULL

  filesToOpen = dir(file.path("rds"), pattern = suffix, full.names = TRUE)

  pb = tcltk::tkProgressBar("Merging data...", paste("File 0/", length(filesToOpen), sep = ""), 0, length(filesToOpen), 200)

  totalParameters = NULL
  pooledData = list()
  foreach::foreach(a = 1:length(filesToOpen)) %do%
  {
    currentData = readRDS(filesToOpen[a])

    currentSample = gsub(paste("(.+)/(.+)", suffix, "\\.rds", sep = ""), "\\2", filesToOpen[a])

    currentData@description$GUID = currentSample

    pooledData[[currentSample]] = currentData

    tcltk::setTkProgressBar(pb, a, label = paste("File ", a, "/", length(filesToOpen), sep = ""))
  }

  close(pb)

  pooledData = methods::as(pooledData, "flowSet")
  saveRDS(pooledData, file.path("rds", paste("pooledSamples.rds", sep = "")))

  gc()

  return(pooledData)
}

#' Launch the R Shiny interactive application for visualization of data logicle transformation
#'
#' This function launches a R Shiny interactive application which helps to visually determine and adjust the parameters for logicle transformation of data.
#'
#' @param fs_shiny A flowSet generated using `mergeData()` function. Defaults to `NULL`.
#'
#' @importFrom foreach %do%
#'
#' @export

launchManualLogicleShinyApp = function(fs_shiny = NULL)
{
  o = NULL
  a = NULL
  b = NULL

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

    shiny::titlePanel(shiny::h1(shiny::strong("Tuning of logicle parameters"), align = "center")),

    shiny::sidebarPanel(

      shiny::selectInput("sample", "Sample", samplesValuesList),
      shiny::selectInput("parameter", "Parameter", parametersNamesList),
      shiny::verbatimTextOutput("auto_parameters"),
      shiny::textOutput("status_parameter"),
      shiny::sliderInput("logicle_t", shiny::strong("t = highest value of the dataset"), min = 1, max = 20, value = 1, step = 0.1),
      shiny::textOutput("logicle_t_lin"),
      shiny::sliderInput("logicle_w", shiny::strong("w = linearization width (slope at 0)"), min = 0, max = 6, value = 0.2, step = 0.1),
      shiny::sliderInput("logicle_m", shiny::strong("m = number of decades of transformed data"), min = 0, max = 10, value = 1, step = 0.1),
      shiny::sliderInput("logicle_a", shiny::strong("a = constant to add to data"), min = 0, max = 2, value = 0, step = 0.1),

      shiny::actionButton('save_inputs', 'Save current parameter'),
      shiny::actionButton('export_txt', 'Export txt'),
      shiny::actionButton('clear_all', 'Clear exports'),
      shiny::actionButton('apply_autologicle', 'Apply autologicle'),
      width = 4),

    shiny::mainPanel(
      shiny::plotOutput(outputId = "histoPlot", height = "900px"),
    )
  )

  # Define server logic
  server = function(input, output, session)
  {
    output$logicle_t_lin = shiny::renderText({
      paste("log10(t) = ", format(round(10^input$logicle_t, 0), big.mark = " ", trim = TRUE), sep = "")
    })

    output$auto_parameters = shiny::renderText({

      data = fs_shiny[[as.numeric(input$sample)]]
      dataToUseForTransformationEstimation = data@exprs[, input$parameter]

      currentParameterRange = range(dataToUseForTransformationEstimation)
      currentParameter_t = currentParameterRange[2]
      currentParameter_m = log10(currentParameter_t)
      currentParameter_r = currentParameterRange[1]
      currentParameter_w = (currentParameter_m - log10(currentParameter_t/abs(currentParameter_r)))/2

      currentParameter_max = max(dataToUseForTransformationEstimation)
      currentParameter_min = min(dataToUseForTransformationEstimation)

      currentParameter_mean = mean(dataToUseForTransformationEstimation)

      currentParameter_a = 0

      paste("Data info:\nmin = ", format(round(currentParameter_min, 0), big.mark = " ", trim = TRUE), "\nmax = ", format(round(currentParameter_max, 0), big.mark = " ", trim = TRUE), "\nmean = ", format(round(currentParameter_mean, 0), big.mark = " ", trim = TRUE), "\n\nAutologicle:\nt = ", format(round(currentParameter_t, 0), big.mark = " ", trim = TRUE), " or log10(t) = ", round(log10(currentParameter_t), 1), "\nw = ", round(currentParameter_w, 1), "\nm = ", round(currentParameter_m, 1), "\na = ",
            round(currentParameter_a, 1), sep = "")
    })

    shiny::observe({

      if (file.exists(file.path("output", "1_Transformation", "parametersLogicle.txt")))
      {
        importedTxt = utils::read.table(file.path("output", "1_Transformation", "parametersLogicle.txt"), header = FALSE)
        colnames(importedTxt) = importedTxt[1, ]
        importedTxt = importedTxt[-1, ]

        if (input$parameter %in% colnames(importedTxt))
        {
          parameterDataID = which(colnames(importedTxt) == input$parameter)
          associatedParameterData = importedTxt[, c(1, parameterDataID)]

          shiny::updateSliderInput(session, "logicle_w", value = as.numeric(associatedParameterData[associatedParameterData$parameter == "w", 2]))
          shiny::updateSliderInput(session, "logicle_t", value = log10(as.numeric(associatedParameterData[associatedParameterData$parameter == "t", 2])))
          shiny::updateSliderInput(session, "logicle_a", value = as.numeric(associatedParameterData[associatedParameterData$parameter == "a", 2]))
          shiny::updateSliderInput(session, "logicle_m", value = as.numeric(associatedParameterData[associatedParameterData$parameter == "m", 2]))
        } else
        {
          data = fs_shiny[[as.numeric(input$sample)]]

          dataToUseForTransformationEstimation = data@exprs[, input$parameter]

          currentParameterRange = range(dataToUseForTransformationEstimation)
          currentParameter_t = currentParameterRange[2]
          currentParameter_m = log10(currentParameter_t)
          currentParameter_r = currentParameterRange[1]
          currentParameter_w = (currentParameter_m - log10(currentParameter_t/abs(currentParameter_r)))/2
          currentParameter_a = 0

          shiny::updateSliderInput(session, "logicle_w", value = currentParameter_w)
          shiny::updateSliderInput(session, "logicle_t", value = log10(currentParameter_t))
          shiny::updateSliderInput(session, "logicle_a", value = currentParameter_a)
          shiny::updateSliderInput(session, "logicle_m", value = currentParameter_m)

          output$status_parameter = shiny::renderText("Status: unsaved")
        }
      } else if (file.exists(file.path("rds", "parameters.rds")))
      {
        importedRDS = readRDS(file.path("rds", "parameters.rds"))

        if (input$parameter %in% names(importedRDS))
        {
          parameterDataID = which(names(importedRDS) == input$parameter)
          associatedParameterData = importedRDS[[parameterDataID]]

          shiny::updateSliderInput(session, "logicle_w", value = associatedParameterData$w)
          shiny::updateSliderInput(session, "logicle_t", value = log10(associatedParameterData$t))
          shiny::updateSliderInput(session, "logicle_a", value = associatedParameterData$a)
          shiny::updateSliderInput(session, "logicle_m", value = associatedParameterData$m)
        } else
        {
          data = fs_shiny[[as.numeric(input$sample)]]

          dataToUseForTransformationEstimation = data@exprs[, input$parameter]

          currentParameterRange = range(dataToUseForTransformationEstimation)
          currentParameter_t = currentParameterRange[2]
          currentParameter_m = log10(currentParameter_t)
          currentParameter_r = currentParameterRange[1]
          currentParameter_w = (currentParameter_m - log10(currentParameter_t/abs(currentParameter_r)))/2
          currentParameter_a = 0

          shiny::updateSliderInput(session, "logicle_w", value = currentParameter_w)
          shiny::updateSliderInput(session, "logicle_t", value = log10(currentParameter_t))
          shiny::updateSliderInput(session, "logicle_a", value = currentParameter_a)
          shiny::updateSliderInput(session, "logicle_m", value = currentParameter_m)

          output$status_parameter = shiny::renderText("Status: unsaved")
        }
      } else
      {
        data = fs_shiny[[as.numeric(input$sample)]]

        dataToUseForTransformationEstimation = data@exprs[, input$parameter]

        currentParameterRange = range(dataToUseForTransformationEstimation)
        currentParameter_t = currentParameterRange[2]
        currentParameter_m = log10(currentParameter_t)
        currentParameter_r = currentParameterRange[1]
        currentParameter_w = (currentParameter_m - log10(currentParameter_t/abs(currentParameter_r)))/2
        currentParameter_a = 0

        shiny::updateSliderInput(session, "logicle_w", value = currentParameter_w)
        shiny::updateSliderInput(session, "logicle_t", value = log10(currentParameter_t))
        shiny::updateSliderInput(session, "logicle_a", value = currentParameter_a)
        shiny::updateSliderInput(session, "logicle_m", value = currentParameter_m)

        output$status_parameter = shiny::renderText("Status: unsaved")
      }
    })

    output$histoPlot = shiny::renderPlot({

      data = fs_shiny[[as.numeric(input$sample)]]

      parametersListName = flowCore::markernames(fs_shiny)

      parameterNameId = which(input$parameter == names(flowCore::markernames(fs_shiny)))
      parameterName = as.character(flowCore::markernames(fs_shiny)[parameterNameId])

      decidedTransform = flowCore::logicleTransform(w = input$logicle_w, t = 10^input$logicle_t, m = input$logicle_m, a = input$logicle_a)
      decidedTransformList = flowCore::transformList(input$parameter, decidedTransform)

      res = flowCore::transform(data, decidedTransformList)

      resPlot = as.data.frame(res@exprs, stringsAsFactors = FALSE)

      if (length(parameterName) == 0)
      {
        plotMainTitle = paste("Density curve for parameter ", input$parameter, sep = "")
      } else
      {
        plotMainTitle = paste("Density curve for parameter ", input$parameter, " (", parameterName, ")", sep = "")
      }

      plot = ggplot2::ggplot(resPlot, ggplot2::aes(x = resPlot[, input$parameter])) +
        ggplot2::geom_density(fill = "grey", linewidth = 1) +
        ggplot2::xlim(-1, 8) +
        ggplot2::labs(title = plotMainTitle, x = "Intensity of fluorescence (logicle-transformed)", y = "Density") +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 14), axis.title = ggplot2::element_text(size = 16, face = "bold"), plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5))
      plot
    })

    totalParametersToExport = list()

    shiny::observeEvent(input$save_inputs, {
      exportedValues = shiny::reactiveValuesToList(input)

      if (file.exists(file.path("rds", "parameters.rds")))
      {
        totalParametersToExport = readRDS(file.path("rds", "parameters.rds"))
      }

      totalParametersToExport[[exportedValues$parameter]] = list(w = input$logicle_w, t = 10^input$logicle_t, m = input$logicle_m, a = input$logicle_a)

      saveRDS(totalParametersToExport, file.path("rds", "parameters.rds"))

      output$status_parameter = shiny::renderText("Status: saved")

      print(length(totalParametersToExport))
      print(totalParametersToExport)
    })

    shiny::observeEvent(input$apply_autologicle, {
      data = fs_shiny[[as.numeric(input$sample)]]

      dataToUseForTransformationEstimation = data@exprs[, input$parameter]

      currentParameterRange = range(dataToUseForTransformationEstimation)
      currentParameter_t = currentParameterRange[2]
      currentParameter_m = log10(currentParameter_t)
      currentParameter_r = currentParameterRange[1]
      currentParameter_w = (currentParameter_m - log10(currentParameter_t/abs(currentParameter_r)))/2
      currentParameter_a = 0

      shiny::updateSliderInput(session, "logicle_w", value = currentParameter_w)
      shiny::updateSliderInput(session, "logicle_t", value = log10(currentParameter_t))
      shiny::updateSliderInput(session, "logicle_a", value = currentParameter_a)
      shiny::updateSliderInput(session, "logicle_m", value = currentParameter_m)

      output$status_parameter = shiny::renderText("Status: unsaved")
    })

    shiny::observeEvent(input$export_txt, {
      exportedValues = shiny::reactiveValuesToList(input)

      if (file.exists(file.path("rds", "parameters.rds")))
      {
        totalParametersToExport = readRDS(file.path("rds", "parameters.rds"))

        res = data.frame(matrix(0, ncol = length(totalParametersToExport), nrow = 4), stringsAsFactors = FALSE)
        rownames(res) = names(totalParametersToExport[[1]])
        colnames(res) = names(totalParametersToExport)

        foreach::foreach(a = 1:length(totalParametersToExport)) %do%
        {
          currentParameterDataName = names(totalParametersToExport)[a]

          foreach::foreach(b = 1:length(totalParametersToExport[[a]])) %do%
          {
            currentParameterSubdataValue = totalParametersToExport[[a]][b]
            currentParameterSubdataName = names(totalParametersToExport[[a]])[b]
            res[currentParameterSubdataName, currentParameterDataName] = as.numeric(currentParameterSubdataValue)
          }
        }

        print(res)

        res = cbind(rownames(res), res)
        colnames(res)[1] = "parameter"
        utils::write.table(res, file.path("output", "1_Transformation", "parametersLogicle.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
      }
    })

    shiny::observeEvent(input$clear_all, {
      exportedValues = shiny::reactiveValuesToList(input)

      if (file.exists(file.path("rds", "parameters.rds")))
      {
        unlink(file.path("rds", "parameters.rds"))
      }

      if (file.exists(file.path("output", "1_Transformation", "parametersLogicle.txt")))
      {
        unlink(file.path("output", "1_Transformation", "parametersLogicle.txt"))
      }

      output$status_parameter = shiny::renderText("Status: unsaved")
    })
  }

  shiny::shinyApp(ui = ui, server = server)
}

#' Transform `rds` files
#'
#' This function applies a mathematical transformation to data from each `rds` file.
#'
#' @param parametersToTransform A vector of all parameters to transform in each `rds` file. Defaults to `NULL`.
#'
#' @param transformationMethod A string defining the transformation method to apply. It can be either `logicle` or `biexp`. Defaults to `logicle`.
#'
#' Please note that the `logicle` transformation needs a file generated using the `manualLogicle` R Shiny application which contains the parameters of the `logicle` transformation to be used for each parameter. The `biexp` transformation does not need this file and can be directly applied.
#'
#' @return Generated `rds` files are saved to `rds` directory, overwriting the previous ones.
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#'
#' @export

transformData = function(parametersToTransform = NULL, transformationMethod = "logicle")
{
  a = NULL
  b = NULL
  customWidthBasis = NULL

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

    totalErrors = 0
    totalErrorsParameters = NULL

    foreach::foreach(b = 1:length(parametersToTransform)) %do%
    {
      currentParameterToTransform = parametersToTransform[b]

      if (transformationMethod == "logicle")
      {
        parametersLogicle = utils::read.table(file.path("output", "1_Transformation", "parametersLogicle.txt"), header = FALSE, as.is = TRUE, sep = "\t")
        rownames(parametersLogicle) = parametersLogicle[, 1]
        colnames(parametersLogicle) = parametersLogicle[1, ]
        parametersLogicle = parametersLogicle[-1, -1]
        parametersLogicle = parametersLogicle[, order(colnames(parametersLogicle))]

        currentParameter_w = as.numeric(parametersLogicle["w", currentParameterToTransform])
        currentParameter_m = as.numeric(parametersLogicle["m", currentParameterToTransform])
        currentParameter_t = as.numeric(parametersLogicle["t", currentParameterToTransform])
        currentParameter_a = as.numeric(parametersLogicle["a", currentParameterToTransform])

        decidedTransform = flowCore::logicleTransform(w = currentParameter_w, t = currentParameter_t, m = currentParameter_m, a = currentParameter_a)
        decidedTransformList = flowCore::transformList(currentParameterToTransform, decidedTransform)
      } else if (transformationMethod == "biexp")
      {
        decidedTransform = flowWorkspace::flowjo_biexp(channelRange = 1e+06, maxValue = 262144, pos = 6, neg = 0.25, widthBasis = customWidthBasis)
        decidedTransformList = flowCore::transformList(currentParameterToTransform, decidedTransform)
      }

      currentData = flowCore::transform(currentData, translist = decidedTransformList)
    }

    saveRDS(currentData, file.path("rds", paste(currentFilename, ".rds", sep = "")))

    rm(currentData)
    rm(decidedTransform)
    rm(decidedTransformList)
    gc()
  }

  parallel::stopCluster(cl)
}

#' Export a single `rds` file per parameter
#'
#' This function reconstructs the dataset by exporting a single `rds` file per parameter instead of one `rds` file per sample. Please note that after the completion of the exportation, the sample-specific `rds` file will be deleted.
#'
#' @param parametersToExport A vector of all parameters to export. Defaults to `NULL`.
#'
#' @param nCoresToExploit An integer defining the number of cores to use to export the `rds` files. Defaults to `NULL` which translates to the actual number of cores of the processor minus 1.
#'
#' @return Generated `rds` files are saved to `rds` directory.
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#'
#' @export

exportPerParameter = function(parametersToExport = NULL, nCoresToExploit = NULL)
{
  a = NULL
  b = NULL

  if (is.null(nCoresToExploit) == FALSE)
  {
    coresNumber = nCoresToExploit
  } else
  {
    coresNumber = parallel::detectCores() - 1
  }

  cl = parallel::makeCluster(coresNumber, type = "PSOCK")
  doSNOW::registerDoSNOW(cl)

  pb = utils::txtProgressBar(min = 0, max = length(parametersToExport), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)

  filesToOpen = dir(file.path("rds"), full.names = TRUE)

  foreach::foreach(a = parametersToExport, .packages = c("foreach", "flowCore", "tcltk"), .options.snow = opts) %dopar%
  {
    # Debug only : a = parametersToKeep[1]

    currentParameter = a

    currentParameterTotalData = list()
    currentParameterTotalNames = NULL

    foreach::foreach(b = 1:length(filesToOpen)) %do%
    {
      currentData = readRDS(filesToOpen[b])
      currentFilename = gsub("(.+)/(.+)\\.rds", "\\2", filesToOpen[b])

      currentData = currentData[, currentParameter]
      currentParameterTotalData[[b]] = currentData
      currentParameterTotalNames = c(currentParameterTotalNames, currentFilename)
    }

    names(currentParameterTotalData) = currentParameterTotalNames

    currentParameterFlowset = methods::as(currentParameterTotalData, "flowSet")

    saveRDS(currentParameterFlowset, file.path("rds", paste("step1_", currentParameter, "_raw.rds", sep = "")))

    rm(currentData)
    rm(currentParameterTotalData)
    rm(currentFilename)
    gc()
  }

  unlink(filesToOpen)

  parallel::stopCluster(cl)
}
