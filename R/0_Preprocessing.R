#' Setup the working directory
#'
#' This function setups the working directory. It does not take any input argument but assumes that the `workingDirectory` variable is set to the path of interest using the base R `setwd()` function. From this path, everything will be wiped and `input`, `output` and `rds` directories will be created. The `output`` directory will also contain several other specific directories that will be used in the subsequent analyses.
#'
#' @export

setupWorkingDirectory = function()
{
  workingDirectory = getwd()

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

#' Convert `fcs` files to `rds` files
#'
#' This function converts every `fcs` file in the `input` directory to a `flowFrame` object (from `flowCore` package) encapsulated into a `rds` file. It eventually renames the parameters used in the dataset if desired. This step helps to decrease the overall computing time and complexity of several next steps.
#'
#' @param conversionTable A tabular-delimited text file containing the appropriate information to convert one or several channels to other ones. This is used when the overall staining mix presents the same specificities but slightly different fluorophores depending on the samples. This typically occurs when the cytometer configuration is modified or when the cytometer is different from one batch to another. Defaults to `NULL`.
#'
#' The `conversionTable` table follows a pre-defined format: 4 columns in any order (`from_desc`, `to_desc`, `from_name` and `to_name`) and any given number of line, each line referring to a specific matching to be treated. For instance, if a line has the values `from_desc = CXCR5 B610-ECD-A`, `to_desc = CXCR5 ECD-A`, `from_name = FL2-A` and `to_name = FL11-A`, it means that any occurrence of a parameter (in any `rds` file) named `FL2-A` which also matches the description `CXCR5 B610-ECD-A` will see its values respectively replaced with `FL11-A` and `CXCR5 ECD-A`.
#'
#' Please note that the renaming both affects the parameters of each file AND the parameters in each self-contained compensation matrix.
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
      currentFileParameterDescriptions = as.vector(currentData@parameters@data[, "desc"])
      currentFileParameterNames = as.vector(currentData@parameters@data[, "name"])

      compensationMatricesSlot = as.numeric(which(lengths(flowStats::spillover(currentData)) > 0))
      compensationMatricesSlotName = names(flowStats::spillover(currentData))[compensationMatricesSlot]

      foreach::foreach(b = 1:length(currentFileParameterDescriptions)) %do% {
        currentDescriptionFrom = currentFileParameterDescriptions[b]
        currentDescriptionName = currentFileParameterNames[b]

        currentNameFromMatching = which((conversionTable$from_desc == currentDescriptionFrom) & conversionTable$from_name == currentDescriptionName)

        if (length(currentNameFromMatching) > 0)
        {
          currentDescriptionTo = conversionTable[conversionTable$from_desc == currentDescriptionFrom, "to_desc"]
          currentNameTo = conversionTable[conversionTable$from_desc == currentDescriptionFrom, "to_name"]
          currentNameFrom = conversionTable[conversionTable$from_desc == currentDescriptionFrom, "from_name"]

          rowToReplaceID = which(currentData@parameters@data$desc == currentDescriptionFrom)

          currentData@parameters@data[rowToReplaceID, "desc"] = paste(currentDescriptionTo, "_replaced", sep = "")
          currentData@parameters@data[rowToReplaceID, "name"] = paste(currentNameTo, "_replaced", sep = "")

          colToReplaceID = as.numeric(which(colnames(currentData@exprs) == currentNameFrom))
          colnames(currentData@exprs)[colToReplaceID] = paste(currentNameTo, "_replaced", sep = "")

          matchingCompensationNameID = which(colnames(currentData@description[compensationMatricesSlotName][[1]]) == currentNameFrom)

          if (length(matchingCompensationNameID) > 0)
          {
            colnames(currentData@description[compensationMatricesSlotName][[1]])[matchingCompensationNameID] = paste(currentNameTo, "_replaced", sep = "")
          }
        }
      }

      currentData@parameters@data[, "desc"] = gsub("_replaced", "", currentData@parameters@data[, "desc"])
      currentData@parameters@data[, "name"] = gsub("_replaced", "", currentData@parameters@data[, "name"])
      colnames(currentData@exprs) = gsub("_replaced", "", colnames(currentData@exprs))

      colnames(currentData@description[compensationMatricesSlotName][[1]]) = gsub("_replaced", "", colnames(currentData@description[compensationMatricesSlotName][[1]]))
    }

    currentFileDescriptions = as.vector(currentData@parameters@data[, "desc"])
    currentFileDescriptions[is.na(currentFileDescriptions)] = paste("Empty Description #", c(1:length(which(is.na(currentFileDescriptions)))), sep = "")

    currentParametersInfos = list(as.vector(currentData@parameters@data[, "name"]), currentFileDescriptions)

    names(currentParametersInfos) = c(paste("names_", currentFilename, sep = ""), paste("descriptions_", currentFilename, sep = ""))

    saveRDS(currentData, file.path("rds", paste(currentFilename, ".rds", sep = "")))

    gc()

    return(currentParametersInfos)
  }

  close(pb)

  parallel::stopCluster(cl)

  totalParametersNamesID = grep("names", names(totalParametersInfos))
  totalParametersNames = totalParametersInfos[totalParametersNamesID]
  totalParametersNames = totalParametersNames[[1]]

  totalParametersDescriptionsID = grep("descriptions", names(totalParametersInfos))
  totalParametersDescriptions = totalParametersInfos[totalParametersDescriptionsID]
  totalParametersDescriptions = totalParametersDescriptions[[1]]

  totalParametersData = data.frame(totalParametersNames, totalParametersDescriptions, stringsAsFactors = FALSE)
  totalParametersData = cbind(1:nrow(totalParametersData), totalParametersData)
  colnames(totalParametersData) = c("Parameter_ID", "Parameter_Name", "Parameter_Description")

  return(totalParametersData)
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

    saveRDS(currentData, file.path("rds", paste(currentFilename, ".rds", sep = "")))

    rm(currentData)
    rm(currentFilename)
    gc()
  }

  parallel::stopCluster(cl)
}

#' Compensate `rds` files
#'
#' This function applies the self-contained compensation matrix of each `rds` file to every desired parameter.
#'
#' @param parametersToCompensate A vector of all parameters to compensate in each `rds` file. Defaults to `NULL`.
#'
#' @return Generated `rds` files are saved to `rds` directory, overwriting the previous ones.
#'
#' @importFrom foreach %dopar%
#'
#' @export

compensateData = function(parametersToCompensate = NULL)
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

    compensationMatricesSlot = as.numeric(which(lengths(flowStats::spillover(currentData)) > 0))
    compensationMatrix = flowStats::spillover(currentData)[[compensationMatricesSlot]]
    parametersToKeepID = which(colnames(compensationMatrix) %in% parametersToCompensate)
    currentSampleComp = compensationMatrix[parametersToKeepID, parametersToKeepID]

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
