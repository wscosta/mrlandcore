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

  # substitute Brazilian cells with MapBiomas landcover (years before 1990 held constant at 1990)
  # time_interpolate restructures dim-3 incorrectly for objects read via dimnames bypass;
  # avoid [<-] operator for the same reason (dim mismatch). Use @.Data directly.
  lcBRA <- readSource("MapBiomas", "LandCover")
  getSets(lcBRA, fulldim = FALSE)[3] <- "landuse"
  brazilCells <- getCells(x)[grepl("\\.BRA$", getCells(x))]
  commonCells <- intersect(brazilCells, getCells(lcBRA))
  if (length(commonCells) > 0) {
    cellX      <- match(commonCells, getCells(x))
    cellLc     <- match(commonCells, getCells(lcBRA))
    yrLc       <- getYears(lcBRA)
    yrIdxLc    <- match(ifelse(getYears(x) %in% yrLc, getYears(x), yrLc[1]), yrLc)
    classIdxLc <- match(getItems(x, 3), getItems(lcBRA, 3))
    stopifnot(!any(is.na(classIdxLc)))
    x@.Data[cellX, , ] <- lcBRA@.Data[cellLc, yrIdxLc, classIdxLc]
  }

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

    # Brazil: rainfed = MapBiomas total - MapBiomas irrigated (explicit, no LUH3 management)
    # irrigated is capped at class total to handle edge cases in stage-4 weight distribution
    irrBRA <- readSource("MapBiomas", "Irrigation")
    irrBRA <- irrBRA[, , paste0(crops, ".irrigated")]
    brazilIrrCells <- intersect(brazilCells, getCells(irrBRA))
    stopifnot(length(brazilIrrCells) == length(brazilCells))

    # time_interpolate restructures irrBRA's dim-3 into a wrong shape (290 x 1 x 125)
    # when irrBRA has dot-names from the readMapBiomas dimnames bypass. Avoid it:
    # map x years to irrBRA years manually (constant extrapolation for pre-1990 years).
    cellX     <- match(brazilIrrCells, getCells(x))
    cellIrr   <- match(brazilIrrCells, getCells(irrBRA))
    yrIrr     <- getYears(irrBRA)
    yrIdxIrr <- match(ifelse(getYears(x) %in% yrIrr, getYears(x), yrIrr[1]), yrIrr)
    for (crop in crops) {
      rainfedClass <- paste0(crop, ".rainfed")
      irrigClass   <- paste0(crop, ".irrigated")
      d3Irr  <- which(getItems(x,      3) == irrigClass)
      d3Rain <- which(getItems(x,      3) == rainfedClass)
      d3Src  <- which(getItems(irrBRA, 3) == irrigClass)
      newIrrig  <- irrBRA@.Data[cellIrr, yrIdxIrr, d3Src]
      totalCrop <- x@.Data[cellX, , d3Rain]
      x@.Data[cellX, , d3Irr]  <- newIrrig
      x@.Data[cellX, , d3Rain] <- totalCrop - newIrrig
    }

    # rest of world: standard LUH3 management approach
    nonBrazilCells <- getCells(x)[!grepl("\\.BRA$", getCells(x))]
    x[nonBrazilCells, , "rainfed"] <- x[nonBrazilCells, , "rainfed"] -
      collapseNames(x[nonBrazilCells, , "irrigated"])
    stopifnot(min(x[nonBrazilCells, , "rainfed"]) >= 0)
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
