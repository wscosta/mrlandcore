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
#' @author Wanderson Costa, Pascal Sauer, Miodrag Stevanovic, Alexandre Koberle
#' @seealso
#' [calcLanduseInitialisation()]
#' @examples
#' \dontrun{
#' calcOutput("LUH3")
#' }
calcLUH3 <- function(landuseTypes = "magpie", irrigation = FALSE,
                     cellular = FALSE, yrs = seq(1965, 2020, 5)) {
  # landuseTypes <- "magpie"
  # irrigation <- FALSE
  # cellular <- FALSE
  # yrs <- seq(1965, 2020, 5)
  # yrs <- as.integer(gsub("y", "", yrs))

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

  cropIBGE <- readSource("IBGE", "Cropland")

  # # total per class (sum over all Brazil cells)
  # totalPerClassBefore <- as.data.frame(dimSums(xBra, dim = 1))
  # totalPerClassBefore

  # totalCountry <- dimSums(xBra, dim = c(1, 3))
  # totalCountry


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

  year <- "y1995"

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
      #idx <- which(residual > 0)
      # cellsResidual <- getCells(xBra)[negCells][idx]

      # debugTable <- data.frame(
      #   cell     = cellsResidual,
      #   c3ann    = as.vector(cropPastOrig[negCells, year, "c3ann"][idx]),
      #   c4ann    = as.vector(cropPastOrig[negCells, year, "c4ann"][idx]),
      #   c3per    = as.vector(cropPastOrig[negCells, year, "c3per"][idx]),
      #   c4per    = as.vector(cropPastOrig[negCells, year, "c4per"][idx]),
      #   c3nfx    = as.vector(cropPastOrig[negCells, year, "c3nfx"][idx]),
      #   pastr    = as.vector(cropPastOrig[negCells, year, "pastr"][idx]),
      #   range    = as.vector(cropPastOrig[negCells, year, "range"][idx]),
      #   primn    = as.vector(primnOrig[negCells][idx]),
      #   secdn    = as.vector(secdnOrig[negCells][idx]),
      #   primf    = as.vector(primfOrig[negCells][idx]),
      #   secdf    = as.vector(secdfOrig[negCells][idx]),
      #   urban    = as.vector(urbanOrig[negCells][idx]),
      #   farm     = as.vector(farm[negCells][idx]),
      #   c3annIBGE = as.vector(cropNew[negCells, year, "c3ann"][idx]),
      #   c4annIBGE = as.vector(cropNew[negCells, year, "c4ann"][idx]),
      #   c3perIBGE = as.vector(cropNew[negCells, year, "c3per"][idx]),
      #   c4perIBGE = as.vector(cropNew[negCells, year, "c4per"][idx]),
      #   c3nfxIBGE = as.vector(cropNew[negCells, year, "c3nfx"][idx]),
      #   pastrNew  = as.vector(pastNewSplit[negCells, year, "pastr"][idx]),
      #   rangeNew  = as.vector(pastNewSplit[negCells, year, "range"][idx]),
      #   primnNew  = as.vector(xBra[negCells, year, "primn"][idx]),
      #   secdnNew  = as.vector(xBra[negCells, year, "secdn"][idx]),
      #   primfNew  = as.vector(xBra[negCells, year, "primf"][idx]),
      #   secdfNew  = as.vector(xBra[negCells, year, "secdf"][idx]),
      #   urbanNew  = as.vector(xBra[negCells, year, "urban"][idx]),
      #   residual  = as.vector(residual[idx])
      # )
      # debugTable$cropTotal <- rowSums(debugTable[, c("c3ann", "c3per", "c4ann", "c4per", "c3nfx")])
      # debugTable$pastTotal <- rowSums(debugTable[, c("pastr", "range")])
      # debugTable$cropIBGE <- rowSums(debugTable[, c("c3annIBGE", "c3perIBGE", "c4annIBGE", "c4perIBGE", "c3nfxIBGE")])
      # debugTable$pastNew <- rowSums(debugTable[, c("pastrNew", "rangeNew")])
      # debugTable$total <- rowSums(debugTable[, c("c3ann", "c3per", "c4ann", "c4per", "c3nfx", "pastr", "range",
      #                                            "primn", "secdn", "primf", "secdf", "urban")])
      # debugTable$totalNew <- rowSums(debugTable[, c("cropIBGE", "pastNew", "primnNew", "secdnNew",
      #                                               "primfNew", "secdfNew", "urbanNew")])
      # debugTable$total <- formatC(debugTable$total, format = "f", digits = 8)
      # debugTable$totalNew <- formatC(debugTable$totalNew, format = "f", digits = 8)


      # print(
      #   cbind(
      #     debugTable["cell"],
      #     lapply(debugTable[ , -1], formatC, format = "f", digits = 8)
      #   )
      # )

      # write.table(
      #   debugTable,
      #   file = "debugTable.csv",   # nome do arquivo
      #   sep = ";",                 # separador ponto e vírgula
      #   dec = ".",                 # decimal com ponto
      #   row.names = FALSE,         # não salva índice das linhas
      #   quote = FALSE              # não coloca aspas em strings
      # )

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

  # # total per class (sum over all Brazil cells)
  # totalPerClass <- as.data.frame(dimSums(xBra, dim = 1))
  # totalPerClass

  # totalCountry <- dimSums(xBra, dim = c(1, 3))
  # totalCountry

  x[brazilCells, , ] <- xBra

  if (isTRUE(irrigation)) {
    crops <- c("c3ann", "c3per", "c4ann", "c4per", "c3nfx")
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
    x[, , getItems(irrigLUH, 3)] <- irrigLUH

    x[, , "rainfed"] <- x[, , "rainfed"] - collapseNames(x[, , "irrigated"])
    stopifnot(min(x[, , "rainfed"]) >= 0)
  }

  # # adjust past to avoid negative values
  # .adjustCellsMagpie <- function(pastValues, otherValues, forestValues, urbanValues) {
  #   deficit <- -pastValues

  #   adjustOther <- pmin(deficit, otherValues)
  #   otherValues <- otherValues - adjustOther
  #   deficit <- deficit - adjustOther

  #   idx <- deficit > 0
  #   if (any(idx)) {
  #     adjustForest <- pmin(deficit[idx], forestValues[idx])
  #     forestValues[idx] <- forestValues[idx] - adjustForest
  #     deficit[idx] <- deficit[idx] - adjustForest
  #   }

  #   idx <- deficit > 0
  #   if (any(idx)) {
  #     adjustUrban <- pmin(deficit[idx], urbanValues[idx])
  #     urbanValues[idx] <- urbanValues[idx] - adjustUrban
  #     deficit[idx] <- deficit[idx] - adjustUrban
  #   }

  #   pastValues[deficit == 0] <- 0

  #   list(past = pastValues,
  #        other = otherValues,
  #        forest = forestValues,
  #        urban = urbanValues,
  #        residual = deficit)
  # }

  if (landuseTypes == "magpie") {
    mapping <- toolGetMapping("LUH3.csv", where = "mrlandcore")
    if (isTRUE(irrigation)) {
      mapping <- data.frame(luh3 = getItems(x, 3),
                            land = paste0(mapping[match(getItems(x, 3.1, full = TRUE), mapping$luh3), ]$land,
                                          ".", getItems(x, 3.2, full = TRUE)))
    }
    stopifnot(setequal(getItems(x, 3), mapping$luh3))
    x <- toolAggregate(x, mapping, dim = 3, from = "luh3", to = "land")

    # cropIBGE <- readSource("IBGE", "Cropland")

    # # select Brazil cells
    # brazilCells <- getCells(x)[grepl("\\.BRA$", getCells(x))]
    # xBra <- x[brazilCells, , ]

    # # total per class (sum over all Brazil cells)
    # totalPerClassBefore <- as.data.frame(dimSums(xBra, dim = 1))
    # totalPerClassBefore

    # totalCountry <- dimSums(xBra, dim = c(1, 3))
    # totalCountry

    # # check sums per cell
    # cellSums <- dimSums(xBra, dim = 3)
    # range(as.vector(cellSums))

    # yearsToAdjust <- c("y1995")
    # year <- "y1995"

    # for (year in yearsToAdjust) {
    #   cropOrig <- xBra[, year, "crop"]
    #   pastOrig <- xBra[, year, "past"]

    #   # compute farm = crop + past
    #   farm <- cropOrig + pastOrig
    #   dimnames(farm)[[3]] <- "crop.past"
    #   names(dimnames(farm))[3] <- "landuse"


    #   # replace crop with IBGE
    #   cropNew <- cropIBGE[, year, "cropland"]
    #   dimnames(cropNew)[[3]] <- "crop.past"
    #   names(dimnames(cropNew))[3] <- "landuse"

    #   # realign cropNew with farm
    #   idx <- match(dimnames(farm)[[1]], dimnames(cropNew)[[1]])
    #   cropNew <- cropNew[idx, , drop = FALSE]
    #   idx <- match(dimnames(farm)[[1]], dimnames(cropNew)[[1]])

    #   # compute new past = farming LUH3 - crop IBGE
    #   pastNew <- farm - cropNew
    #   dimnames(pastNew)[[3]] <- "past"
    #   names(dimnames(pastNew))[3] <- "landuse"

    #   # identify cells with negative past
    #   negCells <- pastNew < 0

    #   # number of negative cells
    #   numNegCells <- sum(as.vector(negCells))
    #   numNegCells

    #   # sum of negative values (pastNew < 0)
    #   totalNeg <- sum(as.vector(pastNew[negCells]))
    #   totalNeg

    #   if (any(negCells)) {
    #     otherOrig <- xBra[, year, "other"]
    #     forestOrig <- xBra[, year, "forest"]
    #     urbanOrig <- xBra[, year, "urban"]

    #     adjustedCells <- .adjustCellsMagpie(
    #       pastNew[negCells],
    #       otherOrig[negCells],
    #       forestOrig[negCells],
    #       urbanOrig[negCells]
    #     )

    #     pastNew[negCells] <- adjustedCells$past
    #     xBra[negCells, year, "other"] <- adjustedCells$other
    #     xBra[negCells, year, "forest"] <- adjustedCells$forest
    #     xBra[negCells, year, "urban"] <- adjustedCells$urban

    #     sum(as.vector(pastNew < 0))
    #     residual <- adjustedCells$residual
    #     # If local crop exceeds total available cell area, crop is locally reduced
    #     # to preserve physical feasibility and total area conservation.
    #     if (any(residual > 0)) {
    #       idx <- which(residual > 0)
    #       cellsResidual <- getCells(xBra)[negCells][idx]

    #       debugTable <- data.frame(
    #         cell   = cellsResidual,
    #         crop   = as.vector(cropOrig[negCells][idx]),
    #         past   = as.vector(pastOrig[negCells][idx]),
    #         forest = as.vector(forestOrig[negCells][idx]),
    #         other  = as.vector(otherOrig[negCells][idx]),
    #         urban  = as.vector(urbanOrig[negCells][idx]),
    #         farm   = as.vector(farm[negCells][idx]),
    #         cropIBGE = as.vector(cropNew[negCells][idx]),
    #         pastNew  = as.vector(pastNew[negCells][idx]),
    #         residual = as.vector(residual[idx])
    #       )
    #       debugTable$total <- rowSums(debugTable[, c("crop", "past", "forest", "other", "urban")])
    #       debugTable$total <- formatC(debugTable$total, format = "f", digits = 8)

    #       print(
    #         cbind(
    #           debugTable["cell"],
    #           lapply(debugTable[ , -1], formatC, format = "f", digits = 8)
    #         )
    #       )

    #       # Calculate total available area in each cell
    #       totalAvailable <- as.vector(
    #         cropOrig[negCells] + pastOrig[negCells] + forestOrig[negCells] + otherOrig[negCells] + urbanOrig[negCells]
    #       )

    #       # 1️st Case: total > 0 and residual positive → reduce crop by residual
    #       idxPositive <- which(residual > 0 & totalAvailable > 0)
    #       if (length(idxPositive) > 0) {
    #         cropNew[negCells][idxPositive] <- cropNew[negCells][idxPositive] - residual[idxPositive]
    #         pastNew[negCells][idxPositive] <- 0
    #       }

    #       # 2️nd Case: total == 0 → set crop and past to zero
    #       idxZero <- which(totalAvailable == 0)
    #       if (length(idxZero) > 0) {
    #         cropNew[negCells][idxZero] <- 0
    #         pastNew[negCells][idxZero] <- 0
    #       }

    #     }

    #     # update xBra with new crop and past
    #     xBra[, year, "crop"] <- cropNew
    #     xBra[, year, "past"] <- pastNew

    #     # identify cells with negative past
    #     negxBraCells <- xBra < 0

    #     # number of negative cells
    #     numNegXCells <- sum(as.vector(negxBraCells))
    #     numNegXCells
    #   }
    # }

    # # apply adjustment to not exceed max area
    # dataAdjusted <- .adjustGrid(as.array(xBra))

    # # put adjusted data back into magpie object
    # xBra[, ] <- dataAdjusted
    # # check sums per cell
    # cellSums <- dimSums(xBra, dim = 3)
    # range(as.vector(cellSums))

    # # total per class (sum over all Brazil cells)
    # totalPerClass <- as.data.frame(dimSums(xBra, dim = 1))
    # totalPerClass

    # totalCountry <- dimSums(xBra, dim = c(1, 3))
    # totalCountry

    # x[brazilCells, , ] <- xBra
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
