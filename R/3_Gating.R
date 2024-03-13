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
#' @param customBins A numeric vector defining the number of bins to use for data plotting after interactive gating. Defaults to `256`.
#'
#' @return If `exportAllPlots = TRUE`, the function will output PDF files containing the gating plots for each sample of the dataset. If `subset = TRUE`, the function will return a list of 2 elements named `flowset` (containing the actual gated `flowSet`) and `summary` (containing basic statistics about the gating such as the number of cells before gating, the number of cells gated and the proportion of cells gated).
#'
#' @importFrom foreach %do%
#'
#' @export

gateData = function(flowset = NULL, sampleToPlot = NULL, xParameter = NULL, yParameter = NULL, xlim = NULL, ylim = NULL, subset = FALSE, gateName = NULL, exportAllPlots = FALSE, samplesPerPage = 6, recursivity = FALSE, inverseGating = FALSE, specificGatesSampleIDs = NULL, redrawGate = TRUE, gatingset = NULL, generatedGates = NULL)
{
  a = NULL
  p = NULL
  g = NULL
  h = NULL
  i = NULL

  xParameterOriginal = xParameter
  yParameterOriginal = yParameter

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

  plot = plot + ggplot2::geom_bin_2d(bins = customBins)  # geom_hex() is bugging right now (28/11/2022)

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



  ### TEMPORARY FIX FOR FLOWGATE ###

  library("ggplot2")

  #' Apply gate from gs_gate_interactive
  #'
  #' @param coords The coordinates of the interactively drawn gate.
  #' @param gateType The selected type of gate (from UI).
  #' @param filterId The gate name specified by the user.
  #' @param gg The plot object from vars$plot.
  #'
  #' @return The original GatingSet with the newly drawn gate applied.
  #'
  #' @noRd
  #'
  applyGateClose <- function(
    gs, subset, coords, gateType, filterId, gg, useBiex, bins, xMax, xWidth,
    xPos, xNeg, yMax, yWidth, yPos, yNeg){
    if(gateType == "polygonGate"){
      names(coords) <- c(names(gg[[1]])[[3]], names(gg[[1]])[[4]])
      coords <- as.matrix(coords)
      gate <- flowCore::polygonGate(coords, filterId = filterId)
    } else if(gateType == "spanGate"){
      names(coords) <- c(names(gg[[1]])[[3]])
      gate <- flowCore::rectangleGate(coords, filterId = filterId)
    } else if(gateType == "quadGate"){
      names(coords) <- c(names(gg[[1]])[[3]], names(gg[[1]])[[4]])
      gate <- flowCore::quadGate(coords, filterId = filterId)
    } else if(gateType == "rectangleGate"){
      names(coords) <- c(names(gg[[1]])[[3]], names(gg[[1]])[[4]])
      gate <- flowCore::rectangleGate(coords, filterId = filterId)
    }
    gs_pop_add(gs, gate, parent = subset)
    recompute(gs)

    if(useBiex){
      varsBiex <- list(X = list(
        maxValue = xMax, widthBasis = xWidth, pos = xPos, neg = xNeg),
        Y = list(
          maxValue = yMax, widthBasis = yWidth, pos = yPos, neg = yNeg))
    } else{
      varsBiex <- "unused"
    }
    output <- list("Gate" = gate, "Bins" = bins, "Scaling" = varsBiex)
    return(output)
  }


  #' Gate Handlers for gs_gate_interactive
  #'
  #' These two functions convert the brushes and clicks on the plot generated in
  #' gs_gate_interactive into usable dimensions for flowCore's filter definition
  #' functions
  #'
  #' @param brush a reactive from the shiny app referring to brushes on the plot.
  #' @param click a reactive from the shiny app referring to clicks on the plot.
  #' @param gateType an input from the shiny app showing the type of gate the user
  #'   wishes to draw.
  #'
  #' @return A set of gate coordinates formatted according to the gate type chosen
  #'   by the user.
  #'
  #' @noRd
  #'
  coordBrush <- function(brush, gateType, useBiex, transX, transY){
    if(useBiex){
      brush$xmin <- transX(brush$xmin)
      brush$xmax <- transX(brush$xmax)
      brush$ymin <- transY(brush$ymin)
      brush$ymax <- transY(brush$ymax)
    }

    if(gateType == "rectangleGate"){
      res <- list(
        "X" = c(brush$xmin, brush$xmax), "Y" = c(brush$ymin, brush$ymax))
    } else if(gateType == "spanGate"){
      res <- list("X" = c(brush$xmin, brush$xmax))
    }
    return(res)
  }

  coordClick <- function(click, gateType, useBiex, transX, transY){
    if(useBiex){
      click$x <- transX(click$x)
      click$y <- transY(click$y)
    }
    if(gateType == "polygonGate"){
      res <- data.frame("x" = click$x, "y" = click$y)
    } else if(gateType == "quadGate"){
      res <- list("X" = click$x, "Y" = click$y)
    }
    return(res)
  }

  gateHandler <- function(gateType, brush, click, biex, tX, tY, append){
    switch(gateType,
           rectangleGate = ,
           spanGate = coordBrush(brush, gateType, biex, tX, tY),
           quadGate = coordClick(click, gateType, biex, tX, tY),
           polygonGate = dplyr::bind_rows(
             append, coordClick(click, gateType, biex, tX, tY)),
           stop("Invalid gate type"))
  }


  #' Sequentially apply a manual gating strategy to a GatingSet or list
  #'
  #' This function allows for a "semi-automatic" approach to interactive manual
  #' gating of flow cytometry data. It leverages the purrr package to let you
  #' easily define a gating strategy and then apply it sequentially to a
  #' GatingSet. This will call gs_gate_interactive() once for each line in your
  #' gating strategy, applying it to your GatingSet as soon as you draw each
  #' prompted gate.
  #'
  #' The gating strategy should be a tibble, with each column name corresponding
  #' to one parameter from gs_gate_interactive. Any parameters not specified in
  #' this tibble will either use their defaults from gs_gate_interactive or can be
  #' specified directly in the function call to gs_apply_gating_strategy.
  #' Typically, this gating strategy will have a column for 'filterId', 'dims',
  #' 'subset', and 'coords', but techinicaly only filterId is required. See
  #' examples below for an easy way to construct this strategy using tribble().
  #'
  #'
  #' @param gs A GatingSet or list of GatingSets.
  #' @param gating_strategy A tibble-formatted gating strategy (see examples below)
  #' @param ... Other parameters to pass to gs_gate_interactive(). Note that only
  #'   constant parameters should be supplied here---anything that varies should
  #'   be included in the gating_strategy tibble.
  #'
  #' @return the GatingSet or list of GatingSets with the gates in gating_strategy
  #'   applied as specified.
  #'
  #' @importFrom tibble tribble
  #' @importFrom methods is
  #'
  #' @examples
  #'
  #' fs <- flowCore::read.flowSet(
  #'   path = system.file("extdata", package = "flowGate"), pattern = ".FCS$")
  #'
  #' gs <- flowWorkspace::GatingSet(fs)
  #'
  #' # Note - this is a very rudamentary GatingSet for example purposes only.
  #' # Please see the vignette accompanying this package or the flowWorkspace
  #' # documentation # for a complete look at creating a GatingSet.
  #'
  #' gating_strategy <- tibble::tribble(
  #' ~filterId,      ~dims,                       ~subset,        ~coords,
  #' "Lymphocytes", list("FSC-H", "SSC-H"),       "root",         list(c(0, 3e5), c(0, 3e5)),
  #' "CD45 CD15",   list("CD45 PE", "CD15 FITC"), "Lymphocytes",  list(c(0, 3e5), c(0, 2e5)),
  #' )
  #'
  #'
  #' if(interactive()){
  #' gs_apply_gating_strategy(gs,
  #' gating_strategy = gating_strategy,
  #' bins = 512) # note that extra args for gs_gate_interactive can be supplied.
  #' }
  #' @export
  gs_apply_gating_strategy <- function(gs, gating_strategy, ...){
    if(methods::is(gs, "GatingSet")){
      purrr::pmap(gating_strategy,
                  gs_gate_interactive, gs = gs, ...)
    } else {
      stop("'gs' must be a GatingSet")
    }
  }


  #' Interactive Manual Gating
  #'
  #' \code{gs_gate_interactive} opens a new graphical window where you can draw
  #' rectangle, polygon, 1-D span, or 2-D quadrant gates that will be applied to
  #' an entire GateSet (see the flowWorkspace package for complete information
  #' about GateSets).
  #'
  #' @param gs The GateSet that will be gated on.
  #' @param filterId String that gives the name of the new gate. Must be unique
  #'   (can specify parent gates to aid in this).
  #' @param sample Numeric specifying which of the GatingHierarchy objects (i.e.
  #'   which FCS file/flow sample) that make up the GateSet do you want to use to
  #'   draw the gate? Note that the gate you draw will be applied to all
  #'   GatingHierarchy objects in the GateSet. Defaults to the first
  #'   GatingHierarchy object in the GateSet.
  #' @param dims A list of strings, length-1 or length-2, that specifies the x-
  #'   and y- parameters that you will be gating on. Giving a length-1 list will
  #'   result in a histogram, while a length-2 list will result in a dot-plot.
  #'   Giving a length-3 or longer list will result in only the first two
  #'   dimensions being used, and will generate a warning to say as much. Defaults
  #'   to forward scatter ("FSC-A") and side scatter ("SSC-A").
  #' @param subset String that gives the name of the parent gate that you want to
  #'   sample from. For example, if you wanted to gate all live cells out of a
  #'   previously drawn "lymphocytes" gate, you would specify "lymphocytes" here.
  #'   Defaults to "root" (ungated).
  #' @param regate A boolean specifying whether all gates with a name matching
  #'   \code{filterId} should first be deleted before being re-drawn. Attempting
  #'   to draw a gate with a non-unique \code{filterId} without specifying
  #'   \code{regate = TRUE} will result in an error. Defaults to \code{FALSE}
  #' @param overlayGates List of strings giving the \code{filterId}s of other
  #'   gates to draw on the example plot when gating. Useful for drawing multiple
  #'   gates on the same population (for example, after specifying a marker-low
  #'   population, you can overlay the marker-low gate to aid in drawing a
  #'   marker-high gate). Defaults to \code{NULL} (no overlaid gates).
  #'
  #' @examples
  #'
  #' path_to_fcs <- system.file("extdata", package = "flowGate")
  #' fs <- read.flowSet(path = path_to_fcs,
  #'                    pattern = ".FCS$",
  #'                    full.names = TRUE)
  #' gs <- GatingSet(fs)
  #'
  #' if(interactive()) { # only run in interactive sessions
  #' gs_gate_interactive(gs,
  #'                     filterId = "Lymphocytes",
  #'                     dims = list("FSC-H", "SSC-H"))
  #' }
  #'
  #' # returns gs with the same "Lymphocytes" gate on FSC-H and SSC-H applied to
  #' # the root node (all events) of each sample in the GateSet.
  #'
  #' if(interactive()) {
  #' gs_gate_interactive(gs,
  #'                     filterId = "Live cells",
  #'                     dims = "Viability",
  #'                     subset = "Lymphocytes")
  #' }
  #'
  #' # returns gs with a "Live cells" gate drawn on all cells included in the
  #' # parent "Lymphocytes" gate. This gate would be based on a histogram of a
  #' # marker called Viability, using the first GatingHierarchy sample as an
  #' # example.
  #'
  #' if(interactive()){
  #' gs_gate_interactive(gs,
  #'                     filterId = "Live cells",
  #'                     dims = list("Viability", "SSC-A"),
  #'                     subset = "Lymphocytes",
  #'                     regate = TRUE)
  #' }
  #'
  #' # first deletes the "Live cells" gate drawn above, then adds a new "Live
  #' # cells" gate to the set, this time based on a dot plot of Viability by
  #' # side-scatter.
  #'
  #' if(interactive()){
  #' gs_gate_interactive(gs,
  #'                     filterId = "Dead cells",
  #'                     dims = list("Viability", "SSC-A"),
  #'                     subset = "Lymphocytes",
  #'                     overlayGates = "Live cells")
  #' }
  #'
  #' # returns gs with a "Dead cells" gate drawn on the same example graph that
  #' # was used to draw the "Live cells" gate above. Overlays the "Live cells"
  #' # gate on top of this graph to aid in drawing the "Dead cells" gate.
  #'
  #' @return A list of the interactively-specified parameters, including the drawn
  #'   gate's coordinates, plot bins, and any flowjo biex coefs used to calculate
  #'   those transforms.
  #'
  #' @import flowWorkspace
  #' @import ggcyto
  #' @import BiocManager
  #' @importFrom ggplot2 aes_ aes geom_density scale_x_continuous
  #' @importFrom ggplot2 scale_y_continuous geom_path geom_hex
  #' @importFrom ggplot2 theme element_blank coord_cartesian
  #' @importFrom rlang .data
  #' @importFrom shiny updateTabsetPanel reactive
  #'
  #' @export
  gs_gate_interactive <- function(
    gs, filterId, sample = 1, dims = list("FSC-A", "SSC-A"), subset = "root",
    regate = FALSE, overlayGates = NULL){
    # Delete gate if regating ==================================================
    if(regate == TRUE){gs_pop_remove(gs, filterId)}
    # Server Function ==========================================================
    server <- function(input, output, session) {
      vals <- shiny::reactiveValues(gateCoords = data.frame(
        "x" = numeric(), "y" = numeric()))
      # Biex Handling --------------------------------------------------------
      shiny::observeEvent(input$useBiex, {
        if(input$useBiex){
          updateTabsetPanel(inputId = "biexTab", selected = "biexPanel")
        }else{
          updateTabsetPanel(inputId = "biexTab", selected = "blankPanel")
        }})
      transX <- reactive(flowjo_biexp(
        maxValue = input$xMaxVal, pos = input$xPos, neg = input$xNeg,
        widthBasis = input$xWidth, inverse = TRUE))
      transY <- reactive(flowjo_biexp(
        maxValue = input$yMaxVal, pos = input$yPos, neg = input$yNeg,
        widthBasis = input$yWidth, inverse = TRUE))
      # Prepare main panel plot ----------------------------------------------
      FPlot <- reactive(preparePlot(
        gs, sample, dims, subset, input$bins, input$useCoords,
        c(input$XMin, input$XMax, input$YMin, input$YMax), overlayGates,
        input$gateType, vals$gateCoords, input$useBiex, input$xMaxVal,
        input$xWidth, input$xPos, input$xNeg, input$yMaxVal, input$yWidth,
        input$yPos, input$yNeg))
      output$plot1 <- shiny::renderPlot(FPlot(), height = function() {
        session$clientData$output_plot1_width})
      output$filterId <- shiny::renderText({paste0("Gate Name: ", filterId)})
      output$subset <- shiny::renderText({paste0("subset of: ", subset)})
      # Gate Handling --------------------------------------------------------
      gateH <- reactive(gateHandler(
        input$gateType, input$plot1_brush, input$plot1_click, input$useBiex,
        transX(), transY(), vals$gateCoords))
      shiny::observeEvent(input$plot1_brush, {vals$gateCoords <- gateH()})
      shiny::observeEvent(input$plot1_click, {vals$gateCoords <- gateH()})
      shiny::observeEvent(input$reset, {
        vals$gateCoords <- data.frame("x" = numeric(), "y" = numeric())})
      # Apply gate and close -------------------------------------------------
      shiny::observeEvent(input$done, {
        output <- applyGateClose(
          gs, subset, vals$gateCoords, input$gateType, filterId, FPlot(),
          input$useBiex, input$bins, input$xMaxVal, input$xWidth,
          input$xPos, input$xNeg, input$yMaxVal, input$yWidth, input$yPos,
          input$yNeg)
        shiny::stopApp(output)})}
    shiny::runApp(shiny::shinyApp(ui, server))
  }


  #' Interactively adjust a gate from a GatingSet
  #'
  #' CAUTION: Experimental Function. Still probably has bugs.
  #'
  #' Call gs_gate_transform_interactive to open a small Shiny app to allow for
  #' manual, interactive adjustments to gates. Currently only supports
  #' rectangleGates and polygonGates.
  #'
  #' @param gs The GatingSet containing the gate you want to adjust
  #' @param node String specifying the (unambiguous) name of the node to adjust
  #' @param sample Numeric specifying which sample in the GatingSet to use for
  #'   example purposes. Note that the adjusted gate will be applied to ALL
  #'   samples, not just this one.
  #' @param dims List of characters specifying channel names or marker names to
  #'   plot on x and y axis. Defaults to list("FSC-A", "SSC-A") mostly just to
  #'   make the format clear.
  #' @param overlayGates (optional) string or character vector specifying names of
  #'   gates to draw on the plot but NOT adjust, for ease of adjusting a gate in
  #'   the vicinity of other gates. Leave NULL to not overlay any gates.
  #'
  #' @return NULL, but silently deletes the old gate, adds the new one, and
  #'   recomputes the GatingSet.
  #'
  #' @examples
  #' path_to_fcs <- system.file("extdata", package = "flowGate")
  #' fs <- read.flowSet(path = path_to_fcs,
  #'                    pattern = ".FCS$",
  #'                    full.names = TRUE)
  #' gs <- GatingSet(fs)
  #'
  #' if(interactive()) { # only run in interactive sessions
  #' gs_gate_interactive(gs,
  #'                     filterId = "Lymphocytes",
  #'                     dims = list("FSC-H", "SSC-H"))
  #'
  #' # Adds a lymphocytes gate to the GatingSet (exactly as in gs_gate_interactive)
  #'
  #' gs_gate_transform_interactive(gs,
  #'                               filterId = "Lymphocytes",
  #'                               dims = list("FSC-H", "SSC-H"))
  #' }
  #'
  #' # Opens a window to adjust the gate manually
  #'
  #' @import flowWorkspace
  #' @import ggcyto
  #' @import BiocManager
  #' @importFrom ggplot2 aes_ aes geom_density scale_x_continuous
  #' @importFrom ggplot2 scale_y_continuous geom_path geom_hex
  #' @importFrom ggplot2 theme element_blank coord_cartesian geom_vline
  #' @importFrom rlang .data
  #' @importFrom shiny reactive
  #'
  #' @export
  gs_gate_transform_interactive <- function(
    gs, node, sample = 1, dims = list("FSC-A", "SSC-A"), overlayGates = NULL){
    # Server Function ==========================================================
    transServer <- function(input, output, session) {
      # Prepare main panel plot ----------------------------------------------
      TPlot <- shiny::reactive(prepareTransPlot(
        gs, sample, dims, node, input$transBins, input$transUseCoords,
        c(input$transXMin, input$transXMax, input$transYMin,
          input$transYMax), input$transUseBiex, overlayGates,
        input$transScaleToggle, c(input$transScaleX, input$transScaleY),
        input$transRotate, input$transShiftX, input$transShiftY))

      output$transPlot <- shiny::renderPlot(TPlot(), height = function() {
        session$clientData$output_transPlot_width})
      # Apply gate and close -------------------------------------------------
      shiny::observeEvent(input$transDone, {
        if(input$transScaleToggle == "uniform"){
          scaleDims <- 1
        }else if(input$transScaleToggle == "independent"){
          scaleDims <- 2
        }
        updateGate(
          gs, node, scaleDims,
          scale = c(input$transScaleX, input$transScaleY),
          deg = input$transRotate, dx = input$transShiftX,
          dy = input$transShiftY
        )
        shiny::stopApp()})}
    shiny::runApp(shiny::shinyApp(uiTransform, transServer))
  }

  # Helpers ----------------------------------------------------------------------

  updateGate <- function(gs, node, scaleDims, scale, deg, dx, dy){
    gate <- flowWorkspace::gh_pop_get_gate(gs[[1]], node)
    if(is(gate, "rectangleGate")){
      deg <- NULL
      if(length(purrr::pluck(gate, "parameters")) == 1){
        dy <- NULL
      }
    }
    if(scaleDims == 1){
      scale <- scale[[1]]
    }

    flowCore::transform_gate(gs, y = node, scale = scale, deg = deg, dx = dx,
                             dy = dy)

    flowWorkspace::recompute(gs)
  }


  prepareTransPlot <- function(gs, sample, dims, node, bins, useCoords, coords,
                               useBiex, overlayGates, scaleMode, scale,
                               deg, dx, dy){

    sample.gs <- gs[[sample]]

    gg <- prepTransPlot(sample.gs, dims, node, bins, useCoords, coords, useBiex)

    if(!is.null(overlayGates)){gg <- gg + geom_gate(overlayGates)}

    gate <- flowWorkspace::gh_pop_get_gate(sample.gs, node)

    if(is(gate, "rectangleGate")){
      deg <- NULL
      if(length(purrr::pluck(gate, "parameters")) == 1){
        dy <- NULL
      }
    }
    if(scaleMode == "uniform"){
      scale <- scale[[1]]
    }

    newGate <- flowCore::transform_gate(gate, scale = scale, deg = deg, dx = dx,
                                        dy = dy)

    if(is(newGate, "rectangleGate")){
      newGate <- ggcyto:::fortify.rectangleGate(newGate)
    }else if(is(newGate, "polygonGate")){
      newGate <- ggcyto:::fortify.polygonGate(newGate)
    }

    if(length(purrr::pluck(gate, "parameters")) == 1 & is(gate, "rectangleGate")){
      gg <- gg +
        geom_vline(xintercept = min(newGate[[1]]), colour = "firebrick") +
        geom_vline(xintercept = max(newGate[[1]]), colour = "firebrick")
    }else{
      gg <- gg + geom_path(data = newGate, colour = "firebrick")
    }

    gg <- ggcyto::as.ggplot(gg)

    return(gg)
  }

  prepTransPlot <- function(sample.gs, dims, node, bins, useCoords, coords, useBiex){
    if(length(dims) > 2){
      warning("The first two dims will be used, the others discarded.")
    }

    if(length(dims) == 1){
      gg <- ggcyto::ggcyto(sample.gs, aes(!!dims[[1]])) +
        geom_density() +
        geom_gate(node, colour = "grey50") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme_flowGate
      if(useCoords){
        gg <- gg + coord_cartesian(xlim = c(coords[[1]], coords[[2]]))
      }
    } else {
      gg <- ggcyto::ggcyto(sample.gs, aes(!!dims[[1]], !!dims[[2]])) +
        geom_hex(bins = bins) +
        geom_gate(node, colour = "grey50") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme_flowGate
      if(useCoords){
        gg <- gg + coord_cartesian(xlim = c(coords[[1]], coords[[2]]),
                                   ylim = c(coords[[3]], coords[[4]]))
      }
    }

    if(useBiex){
      suppressMessages(if(length(dims) == 1){
        gg <- gg + ggcyto::scale_x_flowjo_biexp()
      }else{
        gg <- gg + ggcyto::scale_x_flowjo_biexp() +
          ggcyto::scale_y_flowjo_biexp()
      })
    }
    return(gg)
  }


  #' Prepare Plot for Interactive Gating
  #'
  #' A helper function for gs_gate_interactive to prepare the click-able flow
  #' plots from ggcyto. All params are pulled from the parent call to
  #' gs_gate_interactive.
  #'
  #' @param gs The GatingSet to use.
  #' @param filterId Name of the gate to be drawn.
  #' @param sample Which sample within the GatingSet to use for the plot (gates
  #'   will be applied to all samples).
  #' @param dims List of X and (optionally) Y parameters to plot.
  #' @param subset The parent gate from which to draw the current plot.
  #' @param bins Passed to geom_hex; how fine resolution the plot should be.
  #' @param coords Passed to coord_cartesian; list of length-2 specifying the plot
  #'   dimensions.
  #' @param regate Boolean; should the gs first be stripped of gates matching the
  #'   filterId?
  #' @param overlayGates List of filterIds to plot on top of the current plot.
  #'
  #' @return A ggplot object ready to pass into the shiny app.
  #' @noRd
  #'
  preparePlot <- function(gs, sample, dims, subset, bins, useCoords, coords, overlayGates,
                          addGateType, addCoords, useBiex, x_max, x_wide, x_pos,
                          x_neg, y_max, y_wide, y_pos, y_neg){

    sample.gs <- gs[[sample]]

    gg <- prepMainPlot(sample.gs, dims, subset, bins, useCoords, coords)

    if(!is.null(overlayGates)){gg <- gg + geom_gate(overlayGates)}

    gg <- drawGates(gg, addGateType, addCoords)

    gg <- biexAdjust(gg, useBiex, dims, x_max, x_wide, x_pos, x_neg, y_max,
                     y_wide, y_pos, y_neg)

    gg <- ggcyto::as.ggplot(gg)

    return(gg)
  }

  # Helper Functions -------------------------------------------------------------

  prepMainPlot <- function(sample.gs, dims, subset, bins, useCoords, coords){
    if(length(dims) > 2){
      warning("gs_gate_interactive can only handle one or two dims.
                The first two dims will be used, the others discarded.")
    }

    if(length(dims) == 1){
      gg <- ggcyto::ggcyto(sample.gs, aes(!!dims[[1]]), subset = subset) +
        geom_density() + scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) + theme_flowGate
      if(useCoords){
        gg <- gg + coord_cartesian(xlim = c(coords[[1]], coords[[2]]))
      }
    } else {
      gg <- ggcyto::ggcyto(
        sample.gs, aes(!!dims[[1]], !!dims[[2]]), subset = subset) +
        geom_hex(bins = bins) + scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) + theme_flowGate
      if(useCoords){
        gg <- gg + coord_cartesian(xlim = c(coords[[1]], coords[[2]]),
                                   ylim = c(coords[[3]], coords[[4]]))
      }
    }

    return(gg)
  }

  theme_flowGate <- ggplot2::theme_gray() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(hjust = 0.5))

  drawGates <- function(gg, addGateType, addCoords){
    if(addGateType == "polygonGate"){
      if(!is.null(addCoords) & nrow(addCoords) > 1){
        gg <- gg + geom_path(
          data = addCoords, aes(.data$x, .data$y), inherit.aes = FALSE)
      }
    }else if(addGateType == "quadGate"){
      gg <- gg + ggplot2::geom_vline(xintercept = addCoords$X) +
        ggplot2::geom_hline(yintercept = addCoords$Y)
    }

    return(gg)
  }

  biexAdjust <- function(gg, useBiex, dims, x_max, x_wide, x_pos, x_neg, y_max,
                         y_wide, y_pos, y_neg){
    if(useBiex){
      suppressMessages(if(length(dims)==1){
        gg <- gg + ggcyto::scale_x_flowjo_biexp(
          maxValue = x_max, widthBasis = x_wide, pos = x_pos, neg = x_neg)
      } else{
        gg <- gg + ggcyto::scale_x_flowjo_biexp(
          maxValue = x_max, widthBasis = x_wide, pos=x_pos, neg=x_neg) +
          ggcyto::scale_y_flowjo_biexp(
            maxValue = y_max, widthBasis = y_wide, pos = y_pos, neg = y_neg)
      })
    }

    return(gg)
  }


  sliderInput_MaxVal <- function(id, tag){
    shiny::sliderInput(id, tag, min = -1000, max = 300000, value = 250000)
  }

  sliderInput_Width <- function(id, tag){
    shiny::sliderInput(id, tag, min = -1000, max = -1, value = -10)
  }

  sliderInput_Neg <- function(id, tag){
    shiny::sliderInput(id, tag, min = 0, max = 1, value = 0)
  }

  sliderInput_Pos <- function(id, tag){
    shiny::sliderInput(id, tag, min = 2, max = 7, value = 4, step = 0.1)
  }

  ui <- shiny::fluidPage(
    shiny::titlePanel("Draw your gate"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::actionButton("reset", "Reset"),
        shiny::actionButton("done", "Done"),
        shiny::sliderInput("bins", "Bins", min = 2, max = 2048, value = 256),
        shiny::radioButtons("gateType", "Gate Type:", c(
          "Rectangle" = "rectangleGate", "Polygon" = "polygonGate",
          "Span" = "spanGate", "Quadrant" = "quadGate"),
          selected = "rectangleGate"),
        shiny::checkboxInput("useCoords", "Enable Manual Coords?"),
        shiny::numericInput("XMin", "X Minimum", -1000),
        shiny::numericInput("XMax", "X Maximum", 50000),
        shiny::numericInput("YMin", "Y Minimum", -1000),
        shiny::numericInput("YMax", "Y Maximum", 50000),
        shiny::checkboxInput("useBiex", "Use FlowJo Biex?"),
        shiny::tabsetPanel(id = "biexTab", type = "hidden",
                           shiny::tabPanel("blankPanel", " "),
                           shiny::tabPanel(
                             "biexPanel",
                             sliderInput_MaxVal("xMaxVal", "X Max Value"),
                             sliderInput_Width("xWidth", "X Width Basis"),
                             sliderInput_Neg("xNeg", "X Extra Negative Decades"),
                             sliderInput_Pos("xPos", "X Positive Decades"),
                             sliderInput_MaxVal("yMaxVal", "Y Max Value"),
                             sliderInput_Width("yWidth", "Y Width Basis"),
                             sliderInput_Neg("yNeg", "Y Extra Negative Decades"),
                             sliderInput_Pos("yPos", "Y Positive Decades")))),
      # Main panel for displaying outputs ----------------------------------
      shiny::mainPanel(
        shiny::textOutput("filterId"),
        shiny::textOutput("subset"),
        shiny::plotOutput("plot1",
                          height="auto",
                          click = "plot1_click",
                          brush = shiny::brushOpts(id = "plot1_brush")))))


  uiTransform <- shiny::fluidPage(
    shiny::titlePanel("Update a Gate"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::actionButton("transReset", "Reset"),
        shiny::actionButton("transDone", "Done"),
        shiny::radioButtons("transScaleToggle", "Scale Type:",
                            c("Uniform" = "uniform",
                              "Independent" = "independent")),
        shiny::numericInput("transScaleX", "Scale Factor (Uniform/X)", 1, step = 0.05),
        shiny::numericInput("transScaleY", "Scale Factor (Y)", 1, step = 0.05),
        shiny::sliderInput("transRotate", "Rotation (deg)", value = 0,
                           min = -360, max = 360),
        shiny::numericInput("transShiftX", "Shift Gate (X)", 0, step = 10),
        shiny::numericInput("transShiftY", "Shift Gate (Y)", 0, step = 10),
        shiny::sliderInput("transBins", "Bins", min = 2, max = 2048, value = 256),
        shiny::checkboxInput("transUseBiex", "Use FlowJo Biex?"),
        shiny::checkboxInput("transUseCoords", "Enable Manual Coords?"),
        shiny::numericInput("transXMin", "X Minimum", -1000),
        shiny::numericInput("transXMax", "X Maximum", 20000),
        shiny::numericInput("transYMin", "Y Minimum", -1000),
        shiny::numericInput("transYMax", "Y Maximum", 20000)),
      # Main panel for displaying outputs ----------------------------------
      shiny::mainPanel(
        shiny::plotOutput("transPlot",
                          height="auto"))))


  ###


  ## Put again flowGate:: namespace when the fix is incorporated into flowGate officially

  globalGate = gs_gate_interactive(gatingset, filterId = "globalGate", dims = list(xParameter, yParameter))



  globalGate = globalGate$Gate



  generatedGates = vector(mode = "list", length = length(flowset))

  foreach::foreach(g = 1:length(flowset)) %do%
    {
      generatedGates[[g]] = globalGate
    }

  names(generatedGates) = Biobase::phenoData(flowset)$name








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

if(length(grep("specialGate", nodes_temp)) == 0)
{

  currentSpecificGate = gs_gate_interactive(gatingset, sample = currentSampleSpecificGate, filterId = "specialGate", dims = list(xParameter, yParameter), overlayGates = "globalGate")
} else
{

  currentSpecificGate = gs_gate_interactive(gatingset, sample = currentSampleSpecificGate, filterId = "specialGate", dims = list(xParameter, yParameter), overlayGates = "globalGate", regate = TRUE)

}







          # currentSpecificGate = gs_gate_interactive(gatingset[[currentSampleSpecificGate]], filterId = rownames(phenoData(flowset))[currentSampleSpecificGate], dims = list(xParameter, yParameter), overlayGates = rownames(phenoData(flowset))[currentSampleSpecificGate], regate = TRUE)


          generatedGates[[currentSampleName]] = currentSpecificGate$Gate

        }




  }





#
# nodes = gs_get_pop_paths(gatingset, path = "auto")
#
#
# generatedGates = gs_pop_get_gate(gatingset, nodes[2])


}





    if (inverseGating == TRUE)
    {
      plot = plot + ggcyto::geom_gate(generatedGates) + ggcyto::geom_stats(size = 3, color = "red")



      plot = plot + ggplot2::annotate("text", x = xlim[1], y = ylim[1], hjust = 0, vjust = "bottom", label = "Complementary gate!", color = "red", size = 3)

foreach::foreach(i = 1:length(generatedGates)) %do%
  {
  generatedGates[[i]] = !generatedGates[[i]]

}


    } else
    {
      plot = plot + ggcyto::geom_gate(generatedGates) + ggcyto::geom_stats(size = 3)
    }

    if (subset == TRUE)
    {
	    # For a very curious and not understandable reason, the "filter" function from "flowCore" does not output the results in the same format when called directly (within a script) or within a function (like here in the "gateData" function). So, the "summary" function cannot work anymore in our case. I had to directly specify which slot to use ("subSet") from the filter output ("resultFilter") then manually apply "summary" on it to make it work. This should be further investigated by the flowCore package maintainer.

	    flowset = methods::as(flowset, "flowSet")

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
    }


  plot = plot + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 6, color = "black"), strip.text.y = ggplot2::element_text(size = 6, color = "black"), axis.text.x = ggplot2::element_text(size = 6), axis.text.y = ggplot2::element_text(size = 6))

  if (recursivity == FALSE)
  {
    print(plot)
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
      gateData(flowset = flowset, sampleToPlot = seq(pagesSamplesList[[q]][1], pagesSamplesList[[q]][2]), xParameter = xParameterOriginal, yParameter = yParameterOriginal, xlim = xlim, ylim = ylim, exportAllPlots = FALSE, recursivity = TRUE, inverseGating = inverseGating, specificGatesSampleIDs = specificGatesSampleIDs, redrawGate = FALSE, gatingset = gatingset, generatedGates = generatedGates)

      if (dir.exists(file.path("output", "3_Gating", paste("x=", xParameterOriginal, "_y=", yParameterOriginal, sep = ""))))
      {
        unlink(file.path("output", "3_Gating", paste("x=", xParameterOriginal, "_y=", yParameterOriginal, "*.*", sep = "")))
      } else
      {
        dir.create(file.path("output", "3_Gating", paste("x=", xParameterOriginal, "_y=", yParameterOriginal, sep = "")))
      }

      ggplot2::ggsave(filename = paste("x=", xParameterOriginal, "_y=", yParameterOriginal, "_samples-", pagesSamplesList[[q]][1], "-to-", pagesSamplesList[[q]][2], ".pdf", sep = ""), device = "pdf", path = file.path("output", "3_Gating", paste("x=", xParameterOriginal, "_y=", yParameterOriginal, sep = "")))
    }
  }


    return(list(gatingset = gatingset, generatedGates = generatedGates))




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
