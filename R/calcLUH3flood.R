#' @title calcLUH3flood
#' @description Prepares the LUH3 historic flood data for usage in MAgPIE,
#' by calculating the share of c3ann with c3ann area in Mha to obtain flood area in Mha.
#'
#' @param cellular      if true: dataset is returned on 0.5 degree resolution,
#'                      if false: return country-level data
#' @param yrs           years to be returned
#'
#' @return magpie object with flood data in Mha
#'
#' @author Wanderson Costa, Pascal Sauer, Miodrag Stevanovic, Alexandre Koberle
#' @seealso
#' [calcLanduseInitialisation()]
#' @examples
#' \dontrun{
#' calcOutput("LUH3flood")
#' }
calcLUH3flood <- function(cellular = FALSE, yrs = seq(1965, 2020, 5)) {
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

  management <- readSource("LUH3", "management", yrs, convert = FALSE)
  states <- readSource("LUH3", "states", yrs, convert = TRUE) # convert to Mha
  # convert from shares to Mha, by multiplying flood share with c3ann in Mha
  x <- as.magpie(management[[paste0("y", yrs, "..flood")]] * states[[paste0("y", yrs, "..c3ann")]])

  x <- .aggregateWithMapping(x)
  x <- .ensureAllCells(x, clustermap)

  names(dimnames(x)) <- c("x.y.iso", "t", "data")

  # Brazil: substitute with MapBiomas flood data (Inundacao, capped at c3ann in script 16)
  # Same pattern as calcLUH3.R irrigation block:
  # - avoid time_interpolate (restructures dim-3 of non-standard magpie from dimnames bypass)
  # - use manual year mapping for constant extrapolation of pre-1990 years
  # - use @.Data direct assignment to bypass magclass [<- dim-3 mismatch
  brazilCells    <- getCells(x)[grepl("\\.BRA$", getCells(x))]
  irrBRA         <- readSource("MapBiomas", "Irrigation")
  irrBRA         <- irrBRA[, , "flood"]
  brazilIrrCells <- intersect(brazilCells, getCells(irrBRA))
  stopifnot(length(brazilIrrCells) == length(brazilCells))

  cellX     <- match(brazilIrrCells, getCells(x))
  cellIrr   <- match(brazilIrrCells, getCells(irrBRA))
  yrIrr     <- getYears(irrBRA)
  yrIdxIrr <- match(ifelse(getYears(x) %in% yrIrr, getYears(x), yrIrr[1]), yrIrr)
  x@.Data[cellX, , 1] <- irrBRA@.Data[cellIrr, yrIdxIrr, 1]

  if (!cellular) {
    x <- mstools::toolConv2CountryByCelltype(x, cells = "lpjcell")
  }

  return(list(x            = x,
              weight       = NULL,
              unit         = "Mha",
              description  = "flood area",
              isocountries = !cellular))
}
