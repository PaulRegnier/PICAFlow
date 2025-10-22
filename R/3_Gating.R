#' Get 2nd and 99th percentiles from parameter distribution
#'
#' This function allows to extract the 2nd and 99th percentiles for any given parameter distribution contained in the dataset.
#'
#' @param flowset A string defining the `flowSet` to use. Defaults to `NULL`.
#'
#' @param sample A numeric vector defining which sample should be used to get the limits. Defaults to `1`.
#'
#' @param parameter A string defining which parameter to extract limits from. Defaults to `NULL`.
#'
#' @return A vector of length 2 containing the 2nd and 99th percentiles, in the first and second position, respectively.
#'
#' @export

getParameterLimits = function(flowset = NULL, sample = 1, parameter = NULL)
{
  limits = c(
    round(as.numeric(quantile(flowset[[sample]]@exprs[, parameter], probs = seq(0, 1, 0.01))[2])),
    round(as.numeric(quantile(flowset[[sample]]@exprs[, parameter], probs = seq(0, 1, 0.01))[100]))
  )

  return(limits)
}

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
#' @param specificGatesSampleIDs A numeric vector defining the sample IDs to use for the generation of sample-wise adapted gates. Defaults to `NULL`
#'
#' @param redrawGate A boolean defining if the gates should be manualy redrawn or not. Defaults to `TRUE`
#'
#' @param gatingset A `gatingset` object used internally to store data in a compatible way with `flowGate` package. This value should not be manually provided apart from what is described in the tutorial. Defaults to `NULL`.
#'
#' @param generatedGates A named list containing the gates that were generated during the gating process. This value should not be manually provided apart from what is described in the tutorial. Defaults to `NULL`.
#'
#' @param customBinWidth A numeric vector defining the width of the bins to use for data plotting after interactive gating. Defaults to `5000`. High values (over `100` or `1000`) are better suited for linear parameters with very broad range (such FSC, SSC, etc.) whereas low values (under `10` or `1`) are better suited for log-transformed parameters with a small range. Defaults to `5000`.
#'
#' @return If `exportAllPlots = TRUE`, the function will output PDF files containing the gating plots for each sample of the dataset. If `subset = TRUE`, the function will return a list of 2 elements named `flowset` (containing the actual gated `flowSet`) and `summary` (containing basic statistics about the gating such as the number of cells before gating, the number of cells gated and the proportion of cells gated).
#'
#' @importFrom foreach %do%
#'
#' @export

gateData = function(flowset = NULL, sampleToPlot = NULL, xParameter = NULL, yParameter = NULL, xlim = NULL, ylim = NULL, subset = FALSE, gateName = NULL, exportAllPlots = FALSE, samplesPerPage = 6, recursivity = FALSE, inverseGating = FALSE, specificGatesSampleIDs = NULL, redrawGate = TRUE, gatingset = NULL, generatedGates = NULL, customBinWidth = 5000)
{
  a = NULL
  p = NULL
  g = NULL
  h = NULL
  i = NULL

  xParameterOriginal = xParameter
  yParameterOriginal = yParameter
  gateNameOriginal = gateName

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

  plot = plot + ggplot2::geom_bin_2d(binwidth = customBinWidth)  # geom_hex() is bugging right now (28/11/2022)

  if (length(sampleToPlot) > 1)
  {
    plot = plot + ggplot2::ggtitle(paste("Displaying samples number: ", paste(sampleToPlot, collapse = ", "), sep = ""))

  } else
  {

    plot = plot + ggplot2::ggtitle(paste("Displaying sample number: ", paste(sampleToPlot, collapse = ", "), sep = ""))

  }


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


  if(is.null(gatingset) == TRUE)
  {

    gatingset = flowWorkspace::GatingSet(flowset)

    shiny::runApp(shinyAppGating(flowset = flowset, sample = 1, param_x = xParameter, param_y = yParameter, gateKind = "globalGate", param_x_minLim = xlim[1], param_x_maxLim = xlim[2], param_y_minLim = ylim[1], param_y_maxLim = ylim[2], bins = customBinWidth))

    globalGate_coordinates = readRDS(file.path("rds", "globalGate_coordinates.rds"))

    globalGate = flowCore::polygonGate(.gate = globalGate_coordinates, filterId = "defaultPolygonGate")


    generatedGates = vector(mode = "list", length = length(flowset))

    foreach::foreach(g = 1:length(flowset)) %do%
      {
        generatedGates[[g]] = globalGate
      }

    names(generatedGates) = Biobase::phenoData(flowset)$name



    flowWorkspace::gs_pop_add(gatingset, generatedGates, parent = "root")


  } else
  {

    if(is.null(generatedGates) == TRUE)
    {

      nodes = flowWorkspace::gs_get_pop_paths(gatingset, path = "auto")


      generatedGates = flowWorkspace::gs_pop_get_gate(gatingset, nodes[2])
    }


    if(redrawGate == TRUE & is.null(specificGatesSampleIDs) == FALSE)
    {




      foreach::foreach(h = 1:length(specificGatesSampleIDs)) %do%
        {
          currentSampleSpecificGate = specificGatesSampleIDs[h]

          currentSampleName = rownames(Biobase::phenoData(flowset))[currentSampleSpecificGate]

          nodes_temp = flowWorkspace::gs_get_pop_paths(gatingset, path = "auto")

          print(paste("Opening sample ", h, "/", length(specificGatesSampleIDs), " for special gating: ", "ID #", currentSampleSpecificGate, " (", currentSampleName, ")", sep = ""))

          shiny::runApp(shinyAppGating(flowset = flowset, sample = currentSampleSpecificGate, param_x = xParameter, param_y = yParameter, gateKind = "specialGate", param_x_minLim = xlim[1], param_x_maxLim = xlim[2], param_y_minLim = ylim[1], param_y_maxLim = ylim[2], bins = customBinWidth))

          currentSpecificGate_coordinates = readRDS(file.path("rds", "specialGate_coordinates.rds"))

          currentSpecificGate = flowCore::polygonGate(.gate = currentSpecificGate_coordinates, filterId = "specialGate")

          generatedGates[[currentSampleName]] = currentSpecificGate

        }


    }


  }


  if (inverseGating == TRUE)
  {

    plot = plot + ggcyto::geom_gate(generatedGates) + ggcyto::geom_stats(size = 3, color = "red")

    plot = plot + ggplot2::annotate("text", x = xlim[1], y = ylim[1], hjust = 0, vjust = "bottom", label = "Complementary gate!", color = "red", size = 3)


  } else
  {
    plot = plot + ggcyto::geom_gate(generatedGates) + ggcyto::geom_stats(size = 3)
  }

  if (subset == TRUE)
  {
    # For a very curious and not understandable reason, the "filter" function from "flowCore" does not output the results in the same format when called directly (within a script) or within a function (like here in the "gateData" function). So, the "summary" function cannot work anymore in our case. I had to directly specify which slot to use ("subSet") from the filter output ("resultFilter") then manually apply "summary" on it to make it work. This should be further investigated by the flowCore package maintainer.

    flowset = methods::as(flowset, "flowSet")


    if (inverseGating == TRUE)
    {
      foreach::foreach(i = 1:length(generatedGates)) %do%
        {
          generatedGates[[i]] = !generatedGates[[i]]

        }

    }



    resultFilter = flowCore::filter(flowset, generatedGates)

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

    return(list(flowset = result, summary = resultFilterStats, generatedGates = generatedGates))

    if (inverseGating == TRUE)
    {
      foreach::foreach(i = 1:length(generatedGates)) %do%
        {
          generatedGates[[i]] = !generatedGates[[i]]

        }

    }

  }


  plot = plot + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 6, color = "black"), strip.text.y = ggplot2::element_text(size = 6, color = "black"), axis.text.x = ggplot2::element_text(size = 6), axis.text.y = ggplot2::element_text(size = 6))

  if (recursivity == FALSE)
  {

    ggplot2::ggsave(filename = paste("name=", gateNameOriginal, "_x=", xParameterOriginal, "_y=", yParameterOriginal, "_testGating.pdf", sep = ""), device = "pdf", path = file.path("output", "3_Gating", paste("name=", gateNameOriginal, "_x=", xParameterOriginal, "_y=", yParameterOriginal, sep = "")), plot = plot, width = 29.7, height = 21, units = "cm")

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
        plot = gateData(flowset = flowset, sampleToPlot = seq(pagesSamplesList[[q]][1], pagesSamplesList[[q]][2]), xParameter = xParameterOriginal, yParameter = yParameterOriginal, xlim = xlim, ylim = ylim, exportAllPlots = FALSE, recursivity = TRUE, inverseGating = inverseGating, specificGatesSampleIDs = specificGatesSampleIDs, redrawGate = FALSE, gatingset = gatingset, generatedGates = generatedGates, customBinWidth = customBinWidth, gateName = gateNameOriginal)

        if (dir.exists(file.path("output", "3_Gating", paste("name=", gateNameOriginal, "_x=", xParameterOriginal, "_y=", yParameterOriginal, sep = ""))))
        {
          unlink(file.path("output", "3_Gating", paste("name=", gateNameOriginal, "_x=", xParameterOriginal, "_y=", yParameterOriginal, "*.*", sep = "")))
        } else
        {
          dir.create(file.path("output", "3_Gating", paste("name=", gateNameOriginal, "_x=", xParameterOriginal, "_y=", yParameterOriginal, sep = "")))
        }

        ggplot2::ggsave(filename = paste("name=", gateNameOriginal, "_x=", xParameterOriginal, "_y=", yParameterOriginal, "_samples-", pagesSamplesList[[q]][1], "-to-", pagesSamplesList[[q]][2], ".pdf", sep = ""), device = "pdf", path = file.path("output", "3_Gating", paste("name=", gateNameOriginal, "_x=", xParameterOriginal, "_y=", yParameterOriginal, sep = "")), plot = plot, width = 29.7, height = 21, units = "cm")
      }
  }


  return(list(gatingset = gatingset, generatedGates = generatedGates))




}

#' Shiny application to gate on cells of interest
#'
#' This Shiny application allows to interactively gate on cells of interest. Warning: this function is called when using the `gateData()` mother function, so it should not be called directly.
#'
#' @param flowset A string defining the `flowSet` to use. Defaults to `NULL`.
#'
#' @param sample A numeric vector defining which sample should be plotted for gating visualization. Defaults to `1`.
#'
#' @param param_x A string defining which parameter to use as x. Defaults to `NULL`.
#'
#' @param param_y A string defining which parameter to use as y. Defaults to `NULL`.
#'
#' @param gateKind A string defining the kind of gate to define (either `globalGate` or `specialGate`). Defaults to `globalGate`.
#'
#' @param param_x_minLim A numeric defining the minimum graph limit to use for the parameter `param_x` on the x-axis. Defaults to `NULL`.
#'
#' @param param_x_maxLim A numeric defining the maximum graph limit to use for the parameter `param_x` on the x-axis. Defaults to `NULL`.
#'
#' @param param_y_minLim A numeric defining the maximum graph limit to use for the parameter `param_y` on the y-axis. Defaults to `NULL`.
#'
#' @param param_y_maxLim A numeric defining the maximum graph limit to use for the parameter `param_y` on the y-axis. Defaults to `NULL`.
#'
#' @param bins A numeric vector defining the width of the bins to use for data plotting after interactive gating. Defaults to `5000` (value passed by the mother `gateData()` function). High values (over `100` or `1000`) are better suited for linear parameters with very broad range (such FSC, SSC, etc.) whereas low values (under `10` or `1`) are better suited for log-transformed parameters with a small range.
#'
#' @return The function will not directly return any value. Instead, it will write a RDS file containing the necessary gate information which can be loaded and used afterwards.
#'
#' @importFrom foreach %do%
#' @importFrom dplyr %>%
#'
#' @export

shinyAppGating = function(flowset = NULL, sample = 1, param_x = NULL, param_y = NULL, gateKind = "globalGate", param_x_minLim = NULL, param_x_maxLim = NULL, param_y_minLim = NULL, param_y_maxLim = NULL, bins = 5000)
{
  # Variables declaration

  x = NULL
  y = NULL
  X = NULL
  Y = NULL
  xend = NULL
  yend = NULL

  # Definition of UI

  ui = shiny::fluidPage(
    shiny::titlePanel(
      paste("Gating on parameters X = ", param_x, " and Y = ", param_y, sep = "")),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::h4("Graph controls"),

        # X and Y parameters limits

        shiny::numericInput("x_min", paste0(param_x, " (min):"), value = param_x_minLim),
        shiny::numericInput("x_max", paste0(param_x, " (max):"), value = param_x_maxLim),
        shiny::numericInput("y_min", paste0(param_y, " (min):"), value = param_y_minLim),
        shiny::numericInput("y_max", paste0(param_y, " (max):"), value = param_y_maxLim),

        shiny::hr(),

        # Slider to tune bin width

        shiny::sliderInput("bin_width", "Bin size",
                    min = 1, max = 10000, value = bins, step = 1),

        shiny::hr(),

        shiny::h4("Step 1: simplified calibration"),
        shiny::p("Calibration must be done on the graph corners that are currently displayed in the box below **like this**."),

        shiny::actionButton("start_cal", "Start calibration"),
        shiny::actionButton("reset_cal", "Reinitialize calibration", class = "btn-danger"),

        shiny::verbatimTextOutput("calibration_status"),

        shiny::hr(),
        shiny::h4("Step 2: gate construction and export"),
        shiny::p("Click on the graph to add vertices to the gate polygon."),

        # Button to reinitialize the polygon only

        shiny::actionButton("reset_polygon", "Reinitialize polygon only", class = "btn-warning"),

        shiny::actionButton("close_polygon", "Close polygon", class = "btn-success"),

        # Export button

        shiny::actionButton("download_coords", "Export gate"),

        shiny::hr(),
        shiny::h5("Polygon status:"),
        shiny::verbatimTextOutput("polygon_status_info"),

        shiny::hr(),
        shiny::h5("Last clicked coordinates:"),
        shiny::verbatimTextOutput("click_info")
      ),

      shiny::mainPanel(
        shiny::plotOutput("dot_plot",
                   click = "plot_click",
                   width = "800px",
                   height = "800px"),
        shiny::hr(),
        shiny::h4("Polygon vertices:"),
        shiny::tableOutput("selected_data")
      )
    )
  )

  # Server definition

  server = function(input, output, session) {

    poly_coords = shiny::reactiveValues(x = NULL, y = NULL, closed = FALSE)

    # Calibration points of the graph
    # BL = Bottom-left
    # BR = Bottom-right
    # TL = Top-left

    calibration_points = shiny::reactiveValues(
      calibration_BL = NULL,
      calibration_BR = NULL,
      calibration_TL = NULL,
      step = 0
    )

    # Reactive element for graph limits set by user

    current_limits = shiny::reactive({
      list(
        x_min = input$x_min,
        x_max = input$x_max,
        y_min = input$y_min,
        y_max = input$y_max,
        x_range = input$x_max - input$x_min,
        y_range = input$y_max - input$y_min
      )
    })

    # ---------------------------------------
    # Calibration and polygon control logics
    # ---------------------------------------

    shiny::observeEvent(input$start_cal, {
      if (calibration_points$step < 4) {
        calibration_points$step = 1
        poly_coords$x = poly_coords$y = NULL
        poly_coords$closed = FALSE
        calibration_points$calibration_BL = calibration_points$calibration_BR = calibration_points$calibration_TL = NULL
      }
    })

    shiny::observeEvent(input$reset_cal, {
      calibration_points$calibration_BL = calibration_points$calibration_BR = calibration_points$calibration_TL = NULL
      calibration_points$step = 0
      poly_coords$x = poly_coords$y = NULL
      poly_coords$closed = FALSE
      shiny::showNotification("Calibration and polygon reinitialized.", type = "warning", duration = 2)
    })

    shiny::observeEvent(input$reset_polygon, {
      poly_coords$x = poly_coords$y = NULL
      poly_coords$closed = FALSE
      shiny::showNotification("Polygon erased, but calibration preserved.", type = "warning", duration = 2)
    })

    shiny::observeEvent(input$close_polygon, {
      shiny::req(length(poly_coords$x) >= 3)
      poly_coords$closed = TRUE
      shiny::showNotification("Polygon fnished and closed.", type = "message", duration = 3)
    })

    # ---------------
    # Plot click observer
    # ---------------

    shiny::observeEvent(input$plot_click, {

      shiny::req(input$plot_click, input$plot_click$coords_css)
      x_pixel = input$plot_click$coords_css$x
      y_pixel = input$plot_click$coords_css$y

      limits = current_limits()

      # --- A. calibration phase (step 1, 2 or 3) ---

      if (calibration_points$step > 0 && calibration_points$step < 4) {

        current_step = calibration_points$step

        if (current_step == 1) {
          calibration_points$calibration_BL = c(x_pixel, y_pixel)
          message = paste0("Bottom-left calibration point (", x_pixel, ", ", y_pixel, ") saved. Next calibration click expected: (", limits$x_max, ", ", limits$y_min, ")")
        } else if (current_step == 2) {
          calibration_points$calibration_BR = c(x_pixel, y_pixel)
          message = paste0("Bottom-right calibration point (", x_pixel, ", ", y_pixel, ") saved. Next calibration click expected: (", limits$x_min, ", ", limits$y_max, ")")
        } else if (current_step == 3) {
          calibration_points$calibration_TL = c(x_pixel, y_pixel)
          message = paste0("Top-left calibration point (", x_pixel, ", ", y_pixel, ") saved. Calibration done.")
        }

        calibration_points$step = current_step + 1
        shiny::showNotification(message, type = "message", duration = 2)
        return()
      }

      # --- B. Polygon construction phase (step 4) ---

      if (calibration_points$step == 4 && !poly_coords$closed) {

        factors = calibration_factors()

        # Compute the distance in pixels relative to the bottom-left calibration point (calibration_BL), for the X axis

        pixel_dist_x = x_pixel - factors$calibration_BL_pixel[1]

        # Compute the distance in pixels relative to the bottom-left calibration point (calibration_BL), for the Y axis

        pixel_dist_y = factors$calibration_BL_pixel[2] - y_pixel

        # Conversion in absolute coordinates (relative to limits$x_min and limits$y_min)

        final_x = limits$x_min + (pixel_dist_x * factors$x_scale_factor)
        final_y = limits$y_min + (pixel_dist_y * factors$y_scale_factor)

        # Make sure that the point stays within the user limits (zoom)

        final_x = max(limits$x_min, min(limits$x_max, final_x))
        final_y = max(limits$y_min, min(limits$y_max, final_y))

        poly_coords$x = c(poly_coords$x, final_x)
        poly_coords$y = c(poly_coords$y, final_y)
        shiny::showNotification(paste("Vertice added: (", round(final_x, 0), ", ", round(final_y, 0), ")"), type = "default", duration = 1)
      }
    })

    # ------------------------------
    # Compute the rescaling factors
    # ------------------------------

    calibration_factors = shiny::reactive({
      shiny::req(calibration_points$step == 4)
      limits = current_limits()

      # Coordinates in pixels of calibration points

      calibration_BL_x = calibration_points$calibration_BL[1]
      calibration_BL_y = calibration_points$calibration_BL[2]
      calibration_BR_x = calibration_points$calibration_BR[1]
      calibration_TL_y = calibration_points$calibration_TL[2]

      # Distances in pixels for the display span of visible data

      x_pixel_span = calibration_BR_x - calibration_BL_x
      y_pixel_span = calibration_BL_y - calibration_TL_y

      # Scaling factor: (data span) / (pixels span)

      x_scale_factor = limits$x_range / x_pixel_span
      y_scale_factor = limits$y_range / y_pixel_span

      return(list(
        x_scale_factor = x_scale_factor, # Data Units per Pixel
        y_scale_factor = y_scale_factor, # Data Units per Pixel
        calibration_BL_pixel = calibration_points$calibration_BL # Bottom-left corner in pixels
      ))
    })

    # Creation of the reactive data frame for display and export

    polygon_data_table = shiny::reactive({
      if (length(poly_coords$x) > 0) {
        data.frame(
          ID = 1:length(poly_coords$x),
          X = poly_coords$x,
          Y = poly_coords$y
        )
      } else {
        NULL
      }
    })

    # -------------------
    # Plot display logic
    # -------------------

    output$dot_plot = shiny::renderPlot({

      # 1. Determine graph limits from user inputs

      current_x_min = input$x_min
      current_x_max = input$x_max
      current_y_min = input$y_min
      current_y_max = input$y_max

      # 2. Create the graph with ggcyto

      p = ggcyto::ggcyto(flowset[sample], ggplot2::aes(x = !!rlang::sym(param_x), y = !!rlang::sym(param_y))) +
        ggplot2::geom_bin_2d(binwidth = input$bin_width)

      # 3. Apply the dynamic limits with ggcyto_par_set

      myPars = ggcyto::ggcyto_par_set(limits = list(x = c(current_x_min, current_x_max), y = c(current_y_min, current_y_max)))
      p = p + myPars

      # 4. Force the aspect ratio to 1:1 via ggplot2 theme function (important for calibration)

      p = p + ggplot2::theme(aspect.ratio = 1)

      # 5. Showing the points and polygone

      if (length(poly_coords$x) > 0) {
        poly_data = data.frame(x = as.numeric(poly_coords$x), y = as.numeric(poly_coords$y), stringsAsFactors = FALSE)
        colnames(poly_data) = c(param_x, param_y)

        if (length(poly_coords$x) > 1) {

          p = p + ggplot2::geom_path(data = poly_data, color = "blue", linewidth = 1, group = 1)

          if (poly_coords$closed && length(poly_coords$x) > 2) {
            closing_segment = data.frame(
              x = utils::tail(poly_data[, param_x], 1),
              xend = poly_data[1, param_x],
              y = utils::tail(poly_data[, param_y], 1),
              yend = poly_data[1, param_y]
            )

            p = p + ggplot2::geom_segment(data = closing_segment, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, group = 1), color = "blue", linetype = "solid", linewidth = 1)

          } else {
            NULL
          }
        } else {
          NULL
        }

        p = p + ggplot2::geom_point(data = poly_data, color = "red", size = 5, shape = 19, group = 1) +
          ggplot2::geom_text(data = poly_data, ggplot2::aes(label = seq_along(poly_coords$x)), vjust = -1, color = "red")
      }

      p
    })

    # -----------------------
    # Statuses display logic
    # -----------------------

    output$calibration_status = shiny::renderText({
      limits = current_limits()
      if (calibration_points$step == 0) return("Status: click on \"Start calibration\".")
      if (calibration_points$step == 1) return(paste0("Status: Click expected on the **bottom-left** calibration point with coordinates (", limits$x_min, ", ", limits$y_min, ")."))
      if (calibration_points$step == 2) return(paste0("Status: Click expected on the **bottom-right** calibration point with coordinates (", limits$x_max, ", ", limits$y_min, ")."))
      if (calibration_points$step == 3) return(paste0("Status: Click expected on the **top-left** calibration point with coordinates (", limits$x_min, ", ", limits$y_max, ")."))
      if (calibration_points$step == 4) return("Status: calibration done. You can now proceed to add vertices to the gate.")
      return("")
    })

    output$polygon_status_info = shiny::renderText({
      if (poly_coords$closed) {
        "Status: closed (polygone complete)"
      } else if (length(poly_coords$x) >= 1) {
        paste0("Status : open (", length(poly_coords$x), " vertices)")
      } else {
        "Statut: inactive"
      }
    })

    output$click_info = shiny::renderText({
      if (length(poly_coords$x) == 0) return("No saved points.")

      paste("Last clicked coordinates: (", round(utils::tail(poly_coords$x, 1), 0), ", ", round(utils::tail(poly_coords$y, 1), 0), ")", sep = "")
    })

    output$selected_data = shiny::renderTable({
      if (length(poly_coords$x) > 0) {
        polygon_data_table() %>%
          dplyr::mutate(X = round(X, 0), Y = round(Y, 0)) # Round the coordinates for displaying
      } else {
        data.frame(Message = "Click on the graph after calibration.")
      }
    }, rownames = FALSE)

    # ----------------------------------------------------------------------------------
    # Export + end of application logic
    # ----------------------------------------------------------------------------------

    shiny::observeEvent(input$download_coords, {

      exportTable = polygon_data_table()

      if (is.null(exportTable)) {
        warning("No polygon point saved. The application is stopping without coordinates saved.")
      } else {
        # Make sure the column names are correct for RDS export

        colnames(exportTable) = c("ID", param_x, param_y)

        saveRDS(exportTable[, -1], file.path("rds", paste(gateKind, "_coordinates.rds", sep = "")))
        message(paste0("Polygon coordinates saved in: rds > ", gateKind, "_coordinates.rds"))
      }

      shiny::stopApp(output)
    })
  }

  # --- Launch the Shiny application ---

  app = shiny::shinyApp(ui = ui, server = server)

  return(app)
}

#' Export gating statistics
#'
#' This function allows to export all the gating statistics.
#'
#' @param totalStats A list of gating statistics generated with the `gateData()` method. Defaults to `NULL`.
#'
#' @param filename A character vector defining the custom file name to use for the gating statistics to export. Defaults to `gatingStatistics`.
#'
#' @return Generated text file is saved to `output > 3_Gating` directory.
#'
#' @export

exportGatingStatistics = function(totalStats = NULL, filename = "gatingStatistics")
{
  totalStats = as.data.frame(totalStats)
  totalStats = rbind(colnames(totalStats), totalStats)
  rownames(totalStats)[1] = "Samples"
  utils::write.table(totalStats, file.path("output", "3_Gating", paste(filename, ".txt", sep = "")), quote = FALSE, col.names = FALSE, row.names = TRUE, sep = "\t")
}
