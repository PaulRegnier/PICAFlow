#' Merge `rds` files into a single `rds` file
#'
#' This function merges all the `flowFrame` objects from each `rds` file into a single `flowSet` object encapsulated in a single R-specific `rds` file.
#'
#' @param suffix A character string eventually used to select some `rds` files to merge. Defaults to `NULL`.
#'
#' @param useStructureFromReferenceSample A numeric specifying the sample number to use as reference for channel order and description. Setting this to a value > 0 will reorder each channel and overwrite the descriptions for each flowFrame according to the reference sample. Defaults to `0`, which prevents any reordering and overwriting.
#'
#' @return Generated `rds` file is saved to `rds` directory.
#'
#' @importFrom foreach %do%
#'
#' @export

mergeSamples = function(suffix = NULL, useStructureFromReferenceSample = 0)
{
  a = NULL

  if(file.exists(file.path("rds", "pooledSamples.rds")) == TRUE)
  {

    unlink(file.path("rds", "pooledSamples.rds"))
  }

  filesToOpen = dir(file.path("rds"), pattern = suffix, full.names = TRUE)

  pb = tcltk::tkProgressBar("Merging data...", paste("File 0/", length(filesToOpen), sep = ""), 0, length(filesToOpen), 200)

  pooledData = list()
  foreach::foreach(a = 1:length(filesToOpen)) %do%
  {
    currentData = readRDS(filesToOpen[a])

    currentSample = gsub(paste("(.+)/(.+)", suffix, "\\.rds", sep = ""), "\\2", filesToOpen[a])

    currentData@description$GUID = currentSample

    pooledData[[currentSample]] = currentData


    descriptionsToReplaceID = as.vector(which(is.na(flowCore::parameters(pooledData[[currentSample]])$desc)))

    if(length(descriptionsToReplaceID) > 0)
    {

      flowCore::parameters(pooledData[[currentSample]])$desc[descriptionsToReplaceID] = as.vector(flowCore::parameters(pooledData[[currentSample]])$name[descriptionsToReplaceID])
    }

    if(useStructureFromReferenceSample > 0)
    {


      flowCore::exprs(pooledData[[currentSample]]) = flowCore::exprs(pooledData[[currentSample]])[, colnames(flowCore::exprs(pooledData[[useStructureFromReferenceSample]]))]

      flowCore::parameters(pooledData[[currentSample]])$desc = flowCore::parameters(pooledData[[useStructureFromReferenceSample]])$desc



    }

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

launchTransformationTuningShinyApp = function(fs_shiny = NULL)
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

    pooledValuesShiny[[o]] = flowCore::exprs(fs_shiny[[o]])[sample(nrow(fs_shiny[[o]]), round(min(totalSampleLengths))), ]
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

  transformationsList = as.vector(c("Logicle", "Biexponential", "Arcsinh"))

  # Define UI
  ui = shiny::pageWithSidebar(

    shiny::titlePanel(shiny::h1(shiny::strong("Tuning of transformation parameters"), align = "center")),

    shiny::sidebarPanel(

      shiny::selectInput("sample", "Sample", samplesValuesList),
      shiny::selectInput("parameter", "Parameter", parametersNamesList),
      shiny::selectInput("transformation", "Transformation", transformationsList),
      shiny::verbatimTextOutput("auto_parameters"),
      shiny::textOutput("status_parameter"),


      shiny::tabsetPanel(
        id = "params",
        type = "hidden",
        shiny::tabPanel("Logicle",
                 shiny::actionButton('apply_autologicle', 'Apply autologicle'),

                 shiny::sliderInput("logicle_t", shiny::strong("t = Highest value of the dataset"), min = 1, max = 20, value = 1, step = 0.1),
                 shiny::textOutput("logicle_t_lin"),
                 shiny::sliderInput("logicle_w", shiny::strong("w = Linearization width (slope at 0)"), min = 0, max = 6, value = 0.2, step = 0.1),
                 shiny::sliderInput("logicle_m", shiny::strong("m = Number of decades of transformed data"), min = 0, max = 10, value = 1, step = 0.1),
                 shiny::sliderInput("logicle_a", shiny::strong("a = Constant to add to data"), min = 0, max = 2, value = 0, step = 0.1),





        ),
        shiny::tabPanel("Biexponential",
                 shiny::sliderInput("biexp_a", shiny::strong("a"), min = 0, max = 10, value = 0.5, step = 0.1),
                 shiny::sliderInput("biexp_b", shiny::strong("b"), min = 0, max = 10, value = 1, step = 0.1),
                 shiny::sliderInput("biexp_c", shiny::strong("c"), min = 0, max = 10, value = 0.5, step = 0.1),
                 shiny::sliderInput("biexp_d", shiny::strong("d"), min = 0, max = 10, value = 1, step = 0.1),
                 shiny::sliderInput("biexp_f", shiny::strong("f = Constant bias for the intercept"), min = 0, max = 50, value = 0, step = 1),
                 shiny::sliderInput("biexp_w", shiny::strong("w = Constant bias for the 0 point of the data"), min = 0, max = 50, value = 0, step = 1),




        ),
        shiny::tabPanel("Arcsinh",
                 shiny::sliderInput("arcsinh_a", shiny::strong("a = Shift about 0"), min = 0, max = 1, value = 1, step = 0.01),
                 shiny::sliderInput("arcsinh_b", shiny::strong("b = Scale factor"), min = 0, max = 1, value = 1, step = 0.01),
                 shiny::sliderInput("arcsinh_c", shiny::strong("c = Constant to add to data"), min = 0, max = 50, value = 0, step = 1),



        )
      ),

      shiny::actionButton('save_inputs', 'Save current parameter'),
      shiny::actionButton('export_rds', 'Export rds'),
      shiny::actionButton('clear_all', 'Clear exports'),
      width = 4),

    shiny::mainPanel(
      shiny::plotOutput(outputId = "histoPlot", height = "900px")
    )
  )

  # Define server logic
  server = function(input, output, session)
  {



      output$logicle_t_lin = shiny::renderText({

        if(input$transformation == "Logicle")
        {

          paste("log10(t) = ", format(round(10^input$logicle_t, 0), big.mark = " ", trim = TRUE), sep = "")

        }
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

      if(input$transformation == "Logicle")
      {

        nextText = paste("\n\nAutologicle:\nt = ", format(round(currentParameter_t, 0), big.mark = " ", trim = TRUE), " or log10(t) = ", round(log10(currentParameter_t), 1), "\nw = ", round(currentParameter_w, 1), "\nm = ", round(currentParameter_m, 1), "\na = ", round(currentParameter_a, 1), sep = "")
      } else if(input$transformation == "Biexponential")
      {

        nextText = "\n\nFormula applied: biexp(x) = a*exp(b*(x-w))-c*exp(-d*(x-w))+f"
      } else if(input$transformation == "Arcsinh")
      {
        nextText = "\n\nFormula applied: arcsinh(x) = asinh(a+b*x)+c"

      } else
      {

        nextText = paste("\n\nAutologicle:\nt = ", format(round(currentParameter_t, 0), big.mark = " ", trim = TRUE), " or log10(t) = ", round(log10(currentParameter_t), 1), "\nw = ", round(currentParameter_w, 1), "\nm = ", round(currentParameter_m, 1), "\na = ", round(currentParameter_a, 1), sep = "")
      }

      paste("Data info:\nmin = ", format(round(currentParameter_min, 0), big.mark = " ", trim = TRUE), "\nmax = ", format(round(currentParameter_max, 0), big.mark = " ", trim = TRUE), "\nmean = ", format(round(currentParameter_mean, 0), big.mark = " ", trim = TRUE), nextText, sep = "")
    })

    shiny::observeEvent(input$transformation, {
      shiny::updateTabsetPanel(inputId = "params", selected = input$transformation)
    })

    shiny::observe({


      fileExistsToCheckRds = file.path("rds", "parametersTransformations.rds")


     if (file.exists(fileExistsToCheckRds)) # If a rds file only has already been exported
      {
        importedRDS = readRDS(fileExistsToCheckRds)

        if (input$parameter %in% names(importedRDS)) # If the selected parameter is present in the file, we use this value
        {
          parameterDataID = which(names(importedRDS) == input$parameter)
          associatedParameterData = importedRDS[[parameterDataID]]



          if(associatedParameterData$transformation == "Logicle")
          {
            shiny::updateSelectInput(session, "transformation", selected = associatedParameterData$transformation)

            shiny::updateSliderInput(session, "logicle_w", value = associatedParameterData$w)
            shiny::updateSliderInput(session, "logicle_t", value = log10(associatedParameterData$t))
            shiny::updateSliderInput(session, "logicle_a", value = associatedParameterData$a)
            shiny::updateSliderInput(session, "logicle_m", value = associatedParameterData$m)

          } else if(associatedParameterData$transformation == "Biexponential")
          {
            shiny::updateSelectInput(session, "transformation", selected = associatedParameterData$transformation)

            shiny::updateSliderInput(session, "biexp_a", value = associatedParameterData$a)
            shiny::updateSliderInput(session, "biexp_b", value = associatedParameterData$b)
            shiny::updateSliderInput(session, "biexp_c", value = associatedParameterData$c)
            shiny::updateSliderInput(session, "biexp_d", value = associatedParameterData$d)
            shiny::updateSliderInput(session, "biexp_f", value = associatedParameterData$f)
            shiny::updateSliderInput(session, "biexp_w", value = associatedParameterData$w)

          } else if(associatedParameterData$transformation == "Arcsinh")
          {
            shiny::updateSelectInput(session, "transformation", selected = associatedParameterData$transformation)

            shiny::updateSliderInput(session, "arcsinh_a", value = associatedParameterData$a)
            shiny::updateSliderInput(session, "arcsinh_b", value = associatedParameterData$b)
            shiny::updateSliderInput(session, "arcsinh_c", value = associatedParameterData$c)

          } else
          {

            shiny::updateSelectInput(session, "transformation", selected = "Logicle")

            shiny::updateSliderInput(session, "logicle_w", value = associatedParameterData$w)
            shiny::updateSliderInput(session, "logicle_t", value = log10(associatedParameterData$t))
            shiny::updateSliderInput(session, "logicle_a", value = associatedParameterData$a)
            shiny::updateSliderInput(session, "logicle_m", value = associatedParameterData$m)
          }





        } else # If the selected parameter is not present in the file, we give it default values
        {
          data = fs_shiny[[as.numeric(input$sample)]]

          dataToUseForTransformationEstimation = data@exprs[, input$parameter]



          if(input$transformation == "Logicle")
          {

            currentParameterRange = range(dataToUseForTransformationEstimation)
            currentParameter_t = currentParameterRange[2]
            currentParameter_m = log10(currentParameter_t)
            currentParameter_r = currentParameterRange[1]
            currentParameter_w = (currentParameter_m - log10(currentParameter_t/abs(currentParameter_r)))/2
            currentParameter_a = 0

            shiny::updateSelectInput(session, "transformation", selected = input$transformation)


            shiny::updateSliderInput(session, "logicle_w", value = currentParameter_w)
            shiny::updateSliderInput(session, "logicle_t", value = log10(currentParameter_t))
            shiny::updateSliderInput(session, "logicle_a", value = currentParameter_a)
            shiny::updateSliderInput(session, "logicle_m", value = currentParameter_m)

          } else if(input$transformation == "Biexponential")
          {


            currentParameter_a = 0.5
            currentParameter_b = 1
            currentParameter_c = currentParameter_a
            currentParameter_d = currentParameter_b
            currentParameter_f = 0
            currentParameter_w = 0

            shiny::updateSelectInput(session, "transformation", selected = input$transformation)


            shiny::updateSliderInput(session, "biexp_a", value = currentParameter_a)
            shiny::updateSliderInput(session, "biexp_b", value = currentParameter_b)
            shiny::updateSliderInput(session, "biexp_c", value = currentParameter_c)
            shiny::updateSliderInput(session, "biexp_d", value = currentParameter_d)
            shiny::updateSliderInput(session, "biexp_f", value = currentParameter_f)
            shiny::updateSliderInput(session, "biexp_w", value = currentParameter_w)

          } else if(input$transformation == "Arcsinh")
          {


            currentParameter_a = 1
            currentParameter_b = 1
            currentParameter_c = 0

            shiny::updateSelectInput(session, "transformation", selected = input$transformation)


            shiny::updateSliderInput(session, "arcsinh_a", value = currentParameter_a)
            shiny::updateSliderInput(session, "arcsinh_b", value = currentParameter_b)
            shiny::updateSliderInput(session, "arcsinh_c", value = currentParameter_c)
          } else
          {


            currentParameterRange = range(dataToUseForTransformationEstimation)
            currentParameter_t = currentParameterRange[2]
            currentParameter_m = log10(currentParameter_t)
            currentParameter_r = currentParameterRange[1]
            currentParameter_w = (currentParameter_m - log10(currentParameter_t/abs(currentParameter_r)))/2
            currentParameter_a = 0

            shiny::updateSelectInput(session, "transformation", selected = "Logicle")

            shiny::updateSliderInput(session, "logicle_w", value = currentParameter_w)
            shiny::updateSliderInput(session, "logicle_t", value = log10(currentParameter_t))
            shiny::updateSliderInput(session, "logicle_a", value = currentParameter_a)
            shiny::updateSliderInput(session, "logicle_m", value = currentParameter_m)
          }






          output$status_parameter = shiny::renderText("Status: unsaved")
        }
      } else # If no text or rds files have been exported, we only can give default values
      {
        data = fs_shiny[[as.numeric(input$sample)]]

        dataToUseForTransformationEstimation = data@exprs[, input$parameter]



        if(input$transformation == "Logicle")
        {

          currentParameterRange = range(dataToUseForTransformationEstimation)
          currentParameter_t = currentParameterRange[2]
          currentParameter_m = log10(currentParameter_t)
          currentParameter_r = currentParameterRange[1]
          currentParameter_w = (currentParameter_m - log10(currentParameter_t/abs(currentParameter_r)))/2
          currentParameter_a = 0

          shiny::updateSelectInput(session, "transformation", selected = input$transformation)


          shiny::updateSliderInput(session, "logicle_w", value = currentParameter_w)
          shiny::updateSliderInput(session, "logicle_t", value = log10(currentParameter_t))
          shiny::updateSliderInput(session, "logicle_a", value = currentParameter_a)
          shiny::updateSliderInput(session, "logicle_m", value = currentParameter_m)

        } else if(input$transformation == "Biexponential")
        {


          currentParameter_a = 0.5
          currentParameter_b = 1
          currentParameter_c = currentParameter_a
          currentParameter_d = currentParameter_b
          currentParameter_f = 0
          currentParameter_w = 0

          shiny::updateSelectInput(session, "transformation", selected = input$transformation)


          shiny::updateSliderInput(session, "biexp_a", value = currentParameter_a)
          shiny::updateSliderInput(session, "biexp_b", value = currentParameter_b)
          shiny::updateSliderInput(session, "biexp_c", value = currentParameter_c)
          shiny::updateSliderInput(session, "biexp_d", value = currentParameter_d)
          shiny::updateSliderInput(session, "biexp_f", value = currentParameter_f)
          shiny::updateSliderInput(session, "biexp_w", value = currentParameter_w)

        } else if(input$transformation == "Arcsinh")
        {

          currentParameter_a = 1
          currentParameter_b = 1
          currentParameter_c = 0

          shiny::updateSelectInput(session, "transformation", selected = input$transformation)


          shiny::updateSliderInput(session, "arcsinh_a", value = currentParameter_a)
          shiny::updateSliderInput(session, "arcsinh_b", value = currentParameter_b)
          shiny::updateSliderInput(session, "arcsinh_c", value = currentParameter_c)

        } else
        {

          currentParameterRange = range(dataToUseForTransformationEstimation)
          currentParameter_t = currentParameterRange[2]
          currentParameter_m = log10(currentParameter_t)
          currentParameter_r = currentParameterRange[1]
          currentParameter_w = (currentParameter_m - log10(currentParameter_t/abs(currentParameter_r)))/2
          currentParameter_a = 0

          shiny::updateSelectInput(session, "transformation", selected = "Logicle")


          shiny::updateSliderInput(session, "logicle_w", value = currentParameter_w)
          shiny::updateSliderInput(session, "logicle_t", value = log10(currentParameter_t))
          shiny::updateSliderInput(session, "logicle_a", value = currentParameter_a)
          shiny::updateSliderInput(session, "logicle_m", value = currentParameter_m)

        }



        output$status_parameter = shiny::renderText("Status: unsaved")
      }
    })

    output$histoPlot = shiny::renderPlot({





      data = fs_shiny[[as.numeric(input$sample)]]

      parametersListName = flowCore::markernames(fs_shiny)

      parameterNameId = which(input$parameter == names(flowCore::markernames(fs_shiny)))
      parameterName = as.character(flowCore::markernames(fs_shiny)[parameterNameId])

      if(input$transformation == "Logicle")
      {

        decidedTransform = flowCore::logicleTransform(w = input$logicle_w, t = 10^input$logicle_t, m = input$logicle_m, a = input$logicle_a)


      } else if(input$transformation == "Biexponential")
      {
        decidedTransform = flowCore::biexponentialTransform(a = input$biexp_a, b = input$biexp_b, c = input$biexp_c, d = input$biexp_d, f = input$biexp_f, w = input$biexp_w)


      } else if(input$transformation == "Arcsinh")
      {
        decidedTransform = flowCore::arcsinhTransform(a = input$arcsinh_a, b = input$arcsinh_b, c = input$arcsinh_c)

      } else
      {
        decidedTransform = flowCore::logicleTransform(w = input$logicle_w, t = 10^input$logicle_t, m = input$logicle_m, a = input$logicle_a)



      }





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
        ggplot2::xlim(floor(min(resPlot[, input$parameter])) , ceiling(max(resPlot[, input$parameter]))) +
        ggplot2::labs(title = plotMainTitle, x = paste("Intensity of fluorescence (using ", input$transformation, " transformation)", sep = ""), y = "Density") +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 14), axis.title = ggplot2::element_text(size = 16, face = "bold"), plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5))
      plot
    })




    shiny::observeEvent(input$save_inputs, {
      exportedValues = shiny::reactiveValuesToList(input)


      totalParametersToExport = list()


      fileExistsToCheckRds = file.path("rds", "parametersTransformations.rds")


      if (file.exists(fileExistsToCheckRds))
      {
        totalParametersToExport = readRDS(fileExistsToCheckRds)
      }

      if(input$transformation == "Logicle")
      {

        totalParametersToExport[[exportedValues$parameter]] = list(transformation = input$transformation, w = input$logicle_w, t = 10^input$logicle_t, m = input$logicle_m, a = input$logicle_a)

      } else if(input$transformation == "Biexponential")
      {

        totalParametersToExport[[exportedValues$parameter]] = list(transformation = input$transformation, a = input$biexp_a, b = input$biexp_b, c = input$biexp_c, d = input$biexp_d, f = input$biexp_f, w = input$biexp_w)

      } else if(input$transformation == "Arcsinh")
      {

        totalParametersToExport[[exportedValues$parameter]] = list(transformation = input$transformation, a = input$arcsinh_a, b = input$arcsinh_b, c = input$arcsinh_c)

      } else
      {

        totalParametersToExport[[exportedValues$parameter]] = list(transformation = input$transformation, w = input$logicle_w, t = 10^input$logicle_t, m = input$logicle_m, a = input$logicle_a)

      }


      saveRDS(totalParametersToExport, fileExistsToCheckRds)

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

    shiny::observeEvent(input$export_rds, {
      exportedValues = shiny::reactiveValuesToList(input)




        fileExistsToCheckRds = file.path("rds", "parametersTransformations.rds")



      if (file.exists(fileExistsToCheckRds))
      {
        totalParametersToExport = readRDS(fileExistsToCheckRds)

        saveRDS(totalParametersToExport, file.path("output", "1_Transformation", "parametersTransformations.rds"))


        # res = data.frame(matrix(0, ncol = length(totalParametersToExport), nrow = length(totalParametersToExport[[1]])), stringsAsFactors = FALSE)
        # rownames(res) = names(totalParametersToExport[[1]])
        # colnames(res) = names(totalParametersToExport)
        #
        # foreach::foreach(a = 1:length(totalParametersToExport)) %do%
        # {
        #   currentParameterDataName = names(totalParametersToExport)[a]
        #
        #   foreach::foreach(b = 1:length(totalParametersToExport[[a]])) %do%
        #   {
        #     currentParameterSubdataValue = totalParametersToExport[[a]][b]
        #     currentParameterSubdataName = names(totalParametersToExport[[a]])[b]
        #     res[currentParameterSubdataName, currentParameterDataName] = as.numeric(currentParameterSubdataValue)
        #   }
        # }
        #
        # print(res)

        # res = cbind(rownames(res), res)
        # colnames(res)[1] = "parameter"
        # utils::write.table(res, fileExistsToCheckTxt, quote = TRUE, col.names = TRUE, row.names = FALSE, sep = "\t")
      }
    })

    shiny::observeEvent(input$clear_all, {
      exportedValues = shiny::reactiveValuesToList(input)


      if (file.exists(file.path("output", "1_Transformation", "parametersTransformations.rds")))
      {
        unlink(file.path("output", "1_Transformation", "parametersTransformations.rds"))
      }

      if (file.exists(file.path("rds", "parametersTransformations.rds")))
      {
        unlink(file.path("rds", "parametersTransformations.rds"))
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
#' Please note that this function needs the `parametersTransformations.rds` file located in the `output > 1_Transformation` directory which is generated using the `transformationTuning` R Shiny application. This file contains the type as well as the parameters of the transformation of interest to be used for each flow cytometry parameter in the dataset.
#'
#' @return Generated `rds` files are saved to `rds` directory, overwriting the previous ones.
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#'
#' @export

transformData = function(parametersToTransform = NULL)
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

    transformationsParameters = readRDS(file.path("output", "1_Transformation", "parametersTransformations.rds"))

    foreach::foreach(b = 1:length(parametersToTransform)) %do%
    {
      currentParameterToTransform = parametersToTransform[b]

      currentParametersTransformationParameters = transformationsParameters[[currentParameterToTransform]]
      currentParameterTransformationMethod = currentParametersTransformationParameters$transformation

      if (currentParameterTransformationMethod == "Logicle")
      {

        currentParameter_w = as.numeric(currentParametersTransformationParameters$"w")
        currentParameter_m = as.numeric(currentParametersTransformationParameters$"m")
        currentParameter_t = as.numeric(currentParametersTransformationParameters$"t")
        currentParameter_a = as.numeric(currentParametersTransformationParameters$"a")

        decidedTransform = flowCore::logicleTransform(w = currentParameter_w, t = currentParameter_t, m = currentParameter_m, a = currentParameter_a)
        decidedTransformList = flowCore::transformList(currentParameterToTransform, decidedTransform)
      } else if (currentParameterTransformationMethod == "Biexponential")
      {
        currentParameter_a = as.numeric(currentParametersTransformationParameters$"a")
        currentParameter_b = as.numeric(currentParametersTransformationParameters$"b")
        currentParameter_c = as.numeric(currentParametersTransformationParameters$"c")
        currentParameter_d = as.numeric(currentParametersTransformationParameters$"d")
        currentParameter_f = as.numeric(currentParametersTransformationParameters$"f")
        currentParameter_w = as.numeric(currentParametersTransformationParameters$"w")

        decidedTransform = flowCore::biexponentialTransform(a = currentParameter_a, b = currentParameter_b, c = currentParameter_c, d = currentParameter_d, f = currentParameter_f, w = currentParameter_w)
        decidedTransformList = flowCore::transformList(currentParameterToTransform, decidedTransform)

      } else if (currentParameterTransformationMethod == "Arcsinh")
      {

        currentParameter_a = as.numeric(currentParametersTransformationParameters$"a")
        currentParameter_b = as.numeric(currentParametersTransformationParameters$"b")
        currentParameter_c = as.numeric(currentParametersTransformationParameters$"c")

        decidedTransform = flowCore::arcsinhTransform(a = currentParameter_a, b = currentParameter_b, c = currentParameter_c)
        decidedTransformList = flowCore::transformList(currentParameterToTransform, decidedTransform)

      } else
      {

        currentParameter_w = as.numeric(currentParametersTransformationParameters$"w")
        currentParameter_m = as.numeric(currentParametersTransformationParameters$"m")
        currentParameter_t = as.numeric(currentParametersTransformationParameters$"t")
        currentParameter_a = as.numeric(currentParametersTransformationParameters$"a")

        decidedTransform = flowCore::logicleTransform(w = currentParameter_w, t = currentParameter_t, m = currentParameter_m, a = currentParameter_a)
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
