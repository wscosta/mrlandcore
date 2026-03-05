#' @title calcLUH3
#' @description Prepares the LUH3 historic landuse-dataset for usage in MAgPIE.
#'
#' @param landuseTypes magpie: magpie landuse classes,
#'                     LUH3: original landuse classes
#' @param irrigation    if true: areas are returned separated by irrigated and rainfed,
#'                      if false: total areas (irrigated + rainfed)
#' @param cellular      if true: dataset is returned on 0.5 degree resolution,
#'                      if false: return country-level data
#' @param yrs           years to be returned
#'
#' @return magpie object with land data in Mha
#'
#' @author Wanderson Costa, Alexandre Koberle, Pascal Sauer, Miodrag Stevanovic, 
#' @seealso
#' [calcLanduseInitialisation()]
#' @examples
#' \dontrun{
#' calcOutput("LUH3")
#' }
calcLUH3 <- function(landuseTypes = "magpie", irrigation = FALSE,
                     cellular = FALSE, yrs = seq(1965, 2020, 5)) {
  yrs <- as.integer(gsub("y", "", yrs))

  .aggregateWithMapping <- function(x) {
    mapping <- calcOutput("ResolutionMapping", input = "magpie", target = "luh3", aggregate = FALSE)
    mapping$x.y.iso <- paste0(mapping$cellOriginal, ".", mapping$country)
    mapping <- mapping[, c("cell", "x.y.iso")]

    x <- toolAggregate(x, mapping)
    names(dimnames(x))[1] <- "x.y.iso"

    return(x)
  }

  .ensureAllCells <- function(x, clustermap) {
    missingCells <- setdiff(clustermap$cell, getItems(x, 1))
    x <- magclass::add_columns(x, missingCells, dim = 1, fill = 0)
    stopifnot(setequal(clustermap$cell, getItems(x, 1)))
    x <- x[clustermap$cell, , ]
    return(x)
  }

  clustermap <- readSource("MagpieFulldataGdx", subtype = "clustermap")

  if (!landuseTypes %in% c("magpie", "LUH3")) {
    stop("Unknown landuseTypes = \"", landuseTypes, "\", allowed are: magpie, LUH3")
  }

  states <- readSource("LUH3", "states", yrs)
  stopifnot(terra::units(states) == "Mha")
  # skip plantations (forestry) as it's all zeros in LUH3 at the moment
  states <- states[[grep("pltns", names(states), invert = TRUE)]]
  states <- as.magpie(states)

  x <- .aggregateWithMapping(states)
  getItems(x, 2) <- sub("-.+", "", getItems(x, 2))
  x <- .ensureAllCells(x, clustermap)

  names(dimnames(x)) <- c("x.y.iso", "t", "landuse")
  getSets(x, fulldim = FALSE)[3] <- "landuse"

  # adjusts a single cell so that the sum of its classes does not exceed the ceiling (max cell area in the Equator).
  # the largest class is reduced first, and reductions continue until the total is <= ceil.
  .adjustCell <- function(cellValues, ceil = 0.30914) {
    total <- sum(cellValues)

    while (total > ceil) {
      # find the index of the class with maximum value
      maxClassIdx <- which.max(cellValues)

      # calculate excess
      excess <- total - ceil

      # reduce the maximum class by the excess area, never going below zero
      reduction <- min(excess, cellValues[maxClassIdx])
      cellValues[maxClassIdx] <- cellValues[maxClassIdx] - reduction

      # recalculate the total after reduction
      total <- sum(cellValues)
    }
    return(cellValues)
  }

  # applies .adjustCell to an entire grid (cells x years x land use classes)
  # loops over all cells and years, ensuring that no cell exceeds the ceiling
  .adjustGrid <- function(gridArray, ceil = 0.30914) {
    dims <- dim(gridArray)
    # loop over cells and years
    for (i in 1:dims[1]) {
      for (j in 1:dims[2]) {
        gridArray[i, j, ] <- .adjustCell(gridArray[i, j, ], ceil)
      }
    }
    return(gridArray)
  }

  # select Brazil cells
  brazilCells <- getCells(x)[grepl("\\.BRA$", getCells(x))]
  xBra <- x[brazilCells, , ]

  # apply adjustment to not exceed max area
  dataAdjusted <- .adjustGrid(as.array(xBra))

  # put adjusted data back into magpie object
  xBra[, ] <- dataAdjusted

  # check sums per cell
  cellSums <- dimSums(xBra, dim = 3)
  range(as.vector(cellSums))

  cropIBGE <- readSource("IBGE", "Cropland")

  for (yearInt in yrs) {
    year <- paste0("y", yearInt)

    cropPastOrig <- xBra[, year, c("c3ann", "c4ann", "c3per", "c4per", "c3nfx", "pastr", "range")]
    # Crop = soma das culturas anuais e perenes
    cropOrig <- magclass::dimSums(
      xBra[, year, c("c3ann", "c4ann", "c3per", "c4per", "c3nfx")],
      dim = 3
    )
    dimnames(cropOrig)[[3]] <- "crop.past"
    names(dimnames(cropOrig))[3] <- "landuse"

    # Past = soma de pastagem manejada + natural
    pastOrig <- magclass::dimSums(
      xBra[, year, c("pastr", "range")],
      dim = 3
    )
    dimnames(pastOrig)[[3]] <- "crop.past"
    names(dimnames(pastOrig))[3] <- "landuse"

    # compute farm = crop + past
    farm <- cropOrig + pastOrig
    dimnames(farm)[[3]] <- "crop.past"
    names(dimnames(farm))[3] <- "landuse"

    # replace crop with IBGE
    cropNew <- cropIBGE[, year, c("c3ann", "c3per", "c4ann", "c4per", "c3nfx")]
    names(dimnames(cropNew))[3] <- "landuse"
    idx <- match(dimnames(farm)[[1]], dimnames(cropNew)[[1]])
    cropNew <- cropNew[idx, , drop = FALSE]
    idx <- match(dimnames(farm)[[1]], dimnames(cropNew)[[1]])

    # Crop Total = soma das culturas anuais e perenes
    cropNewTotal <- magclass::dimSums(
      cropNew[, year, c("c3ann", "c4ann", "c3per", "c4per", "c3nfx")],
      dim = 3
    )
    dimnames(cropNewTotal)[[3]] <- "crop.past"
    names(dimnames(cropNewTotal))[3] <- "landuse"

    # realign cropNew with farm
    idx <- match(dimnames(farm)[[1]], dimnames(cropNewTotal)[[1]])
    cropNewTotal <- cropNewTotal[idx, , drop = FALSE]
    idx <- match(dimnames(farm)[[1]], dimnames(cropNewTotal)[[1]])

    # compute new past = farming LUH3 - crop IBGE
    pastNewTotal <- farm - cropNewTotal
    dimnames(pastNewTotal)[[3]] <- "past"
    names(dimnames(pastNewTotal))[3] <- "landuse"

    # Seleciona pastr e range para o ano
    pastNewSplit <- xBra[, year, c("pastr", "range")]

    # Redistribui a diferença usando subsetting direto
    pastNewSplit[, , "pastr"] <- pastNewTotal[, , 1] * 0.35
    pastNewSplit[, , "range"] <- pastNewTotal[, , 1] * 0.65

    sum(pastNewSplit, na.rm = TRUE)

    .adjustPastCellsLUH3 <- function(pastrValues, rangeValues,
                                     primnValues, secdnValues,
                                     primfValues, secdfValues,
                                     urbanValues) {
      # soma das pastagens
      pastTotal <- pastrValues + rangeValues

      # deficit = valores negativos da soma de pastr + range
      deficit <- -pmin(pastTotal, 0)

      # realocação para primn
      adjustPrimn <- pmin(deficit, primnValues)
      primnValues <- primnValues - adjustPrimn
      deficit <- deficit - adjustPrimn

      # realocação para secdn
      idx <- deficit > 0
      if (any(idx)) {
        adjustSecdn <- pmin(deficit[idx], secdnValues[idx])
        secdnValues[idx] <- secdnValues[idx] - adjustSecdn
        deficit[idx] <- deficit[idx] - adjustSecdn
      }

      # realocação para primf
      idx <- deficit > 0
      if (any(idx)) {
        adjustPrimf <- pmin(deficit[idx], primfValues[idx])
        primfValues[idx] <- primfValues[idx] - adjustPrimf
        deficit[idx] <- deficit[idx] - adjustPrimf
      }

      # realocação para secdf
      idx <- deficit > 0
      if (any(idx)) {
        adjustSecdf <- pmin(deficit[idx], secdfValues[idx])
        secdfValues[idx] <- secdfValues[idx] - adjustSecdf
        deficit[idx] <- deficit[idx] - adjustSecdf
      }

      # realocação para urban
      idx <- deficit > 0
      if (any(idx)) {
        adjustUrban <- pmin(deficit[idx], urbanValues[idx])
        urbanValues[idx] <- urbanValues[idx] - adjustUrban
        deficit[idx] <- deficit[idx] - adjustUrban
      }

      pastrValues[deficit == 0] <- 0
      rangeValues[deficit == 0] <- 0

      list(
        pastr = pastrValues,
        range = rangeValues,
        primn = primnValues,
        secdn = secdnValues,
        primf = primfValues,
        secdf = secdfValues,
        urban = urbanValues,
        residual = deficit
      )
    }


    negCells <- (pastNewSplit[, year, "pastr"] + pastNewSplit[, year, "range"]) < 0

    # number of negative cells
    numNegCells <- sum(as.vector(negCells))
    numNegCells

    # sum of negative values (pastNewSplit < 0)
    totalNeg <- sum(as.vector(pastNewSplit[negCells]))
    totalNeg

    if (any(negCells)) {
      primnOrig <- xBra[, year, "primn"]
      secdnOrig <- xBra[, year, "secdn"]
      primfOrig <- xBra[, year, "primf"]
      secdfOrig <- xBra[, year, "secdf"]
      urbanOrig <- xBra[, year, "urban"]

      adjusted <- .adjustPastCellsLUH3(
        pastrValues = pastNewSplit[negCells, year, "pastr"],
        rangeValues = pastNewSplit[negCells, year, "range"],
        primnValues = primnOrig[negCells],
        secdnValues = secdnOrig[negCells],
        primfValues = primfOrig[negCells],
        secdfValues = secdfOrig[negCells],
        urbanValues = urbanOrig[negCells]
      )

      # atualiza o xBra
      pastNewSplit[negCells, year, "pastr"] <- adjusted$pastr
      pastNewSplit[negCells, year, "range"] <- adjusted$range
      xBra[negCells, year, "primn"] <- adjusted$primn
      xBra[negCells, year, "secdn"] <- adjusted$secdn
      xBra[negCells, year, "primf"] <- adjusted$primf
      xBra[negCells, year, "secdf"] <- adjusted$secdf
      xBra[negCells, year, "urban"] <- adjusted$urban

      negCellsAfter <- (pastNewSplit[, year, "pastr"] + pastNewSplit[, year, "range"]) < 0
      sum(as.vector(pastNewSplit[, year, "pastr"] < 0 & pastNewSplit[, year, "range"] < 0))
      sum(negCellsAfter)

      residual <- adjusted$residual
      # If local crop exceeds total available cell area, crop is locally reduced
      # to preserve physical feasibility and total area conservation.
      if (any(residual > 0)) {
        # Calculate total available area in each cell
        cropPastTotal <- dimSums(cropPastOrig, dim = 3)
        totalAvailable <- as.vector(
          cropPastTotal[negCells] +
            primnOrig[negCells] +
            secdnOrig[negCells] +
            primfOrig[negCells] +
            secdfOrig[negCells] +
            urbanOrig[negCells]
        )


        .redistributeCrop <- function(cropVec, deficit) {
          if (deficit <= 0 || sum(cropVec) <= 0) {
            return(cropVec)
          }

          remaining <- cropVec
          rest <- deficit

          while (rest > 1e-12 && any(remaining > 0)) {
            w <- remaining / sum(remaining)
            delta <- w * rest
            delta <- pmin(delta, remaining)

            remaining <- remaining - delta
            rest <- rest - sum(delta)
          }

          return(remaining)
        }


        # 1️st Case: total > 0 and residual positive → reduce crop by residual
        idxPositive <- which(residual > 0 & totalAvailable > 0)
        if (length(idxPositive) > 0) {
          cells <- which(negCells)

          for (k in idxPositive) {
            cell <- cells[k]

            cropVec <- as.numeric(cropNew[cell, , drop = TRUE])
            newCrop <- .redistributeCrop(cropVec, residual[k])

            cropNew[cell, , ] <- newCrop
            pastNewSplit[cell, , ] <- 0
          }
        }

        # 2️nd Case: total == 0 → set crop and past to zero
        idxZero <- which(totalAvailable == 0)
        if (length(idxZero) > 0) {
          cells <- which(negCells)

          for (k in idxZero) {
            cell <- cells[k]

            cropNew[cell, , ]      <- 0
            pastNewSplit[cell, , ] <- 0
          }
        }
      }

      xBra[, year, getItems(cropNew, 3)]      <- cropNew
      xBra[, year, getItems(pastNewSplit, 3)] <- pastNewSplit


      # identify cells with negative past
      negxBraCells <- xBra < 0

      # number of negative cells
      numNegXCells <- sum(as.vector(negxBraCells))
      numNegXCells
    }
  }

  x[brazilCells, , ] <- xBra

  if (isTRUE(irrigation)) {
    crops <- c("c3ann", "c3per", "c4ann", "c4per", "c3nfx")

    # identify cells who are from Brazil
    brazilCells <- getCells(x)[grepl("\\.BRA$", getCells(x))]
    nonBrazilCells <- getCells(x)[!grepl("\\.BRA$", getCells(x))]

    irrigLUH <- readSource("LUH3", "management", yrs, convert = FALSE)
    irrigLUH <- as.magpie(irrigLUH)
    irrigLUH <- irrigLUH[, , paste0("irrig_", crops)]
    stopifnot(0 <= irrigLUH, irrigLUH <= 1)

    getNames(irrigLUH) <- crops
    # convert to Mha by multiplying with cropland in Mha
    irrigLUH <- irrigLUH * states[, , crops]

    irrigLUH <- .aggregateWithMapping(irrigLUH)
    irrigLUH <- .ensureAllCells(irrigLUH, clustermap)

    names(dimnames(irrigLUH)) <- c("x.y.iso", "t", "data")

    x <- add_dimension(x, dim = 3.2, add = "irrigation", nm = "total")
    getItems(x, 3.2, full = TRUE)[getItems(x, 3.1, full = TRUE) %in% crops] <- "rainfed"
    x <- add_columns(x, dim = 3, addnm = paste0(crops, ".irrigated"))

    irrigLUH <- add_dimension(irrigLUH, dim = 3.2, add = "irrigation", nm = "irrigated")
    x[nonBrazilCells, , getItems(irrigLUH, 3)] <- irrigLUH[nonBrazilCells, , ]
    x[nonBrazilCells, , "rainfed"] <- x[nonBrazilCells, , "rainfed"] - collapseNames(x[nonBrazilCells, , "irrigated"])

    yearLabels <- paste0("y", yrs)

    # set irrigated = 0 for Brazil
    x[brazilCells, yearLabels, paste0(crops, ".irrigated")] <- 0

    stopifnot(min(x[, , "rainfed"]) >= 0)

    #code logic after having the irrigation map for Brazil (irrigatedBRA)
    # x[brazilCells, , "*.irrigated"] <- irrigatedBRA
    # x[brazilCells, , "*.rainfed"]  <- x[brazilCells, , "*.rainfed"] - irrigatedBRA
  }

  if (landuseTypes == "magpie") {
    mapping <- toolGetMapping("LUH3.csv", where = "mrlandcore")
    if (isTRUE(irrigation)) {
      mapping <- data.frame(luh3 = getItems(x, 3),
                            land = paste0(mapping[match(getItems(x, 3.1, full = TRUE), mapping$luh3), ]$land,
                                          ".", getItems(x, 3.2, full = TRUE)))
    }
    stopifnot(setequal(getItems(x, 3), mapping$luh3))
    x <- toolAggregate(x, mapping, dim = 3, from = "luh3", to = "land")
  }

  if (!cellular) {
    x <- mstools::toolConv2CountryByCelltype(x, cells = "lpjcell")
  }

  return(list(x            = x,
              weight       = NULL,
              unit         = "Mha",
              description  = "land area for different land use types",
              isocountries = !cellular))
}
