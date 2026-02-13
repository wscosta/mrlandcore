#' @title readIBGE
#' @description
#' Reads IBGE cropland (planted area) time series data for Brazil and converts
#' it into a magpie object.
#' Data was collected from the SIDRA IBGE Platform
#' (https://sidra.ibge.gov.br//).
#'
#' @param subtype Data subtype (currently only "Cropland")
#'
#' @return A magpie object with cropland (planted area) data for Brazil (Mha)
#' @author Wanderson Costa, Alex Koberle, Miodrag Stevanovic
#' @examples
#' \dontrun{
#' readSource("IBGE", subtype="Cropland")
#' }
#'
#' @importFrom utils read.csv

readIBGE <- function(subtype = "Cropland") {

  files <- c(
    #Cropland = "/p/projects/rd3mod/inputdata/sources/IBGE/cropland_1995_all.csv"
    #Cropland = "cropland_1995_all.csv"
    #Cropland = "/p/projects/rd3mod/inputdata/sources/IBGE/crop_planted_area_1995_luh3_all.csv"
    Cropland = "crop_planted_area_1995_luh3_all.csv"

  )

  if (!subtype %in% names(files)) {
    stop("Unknown subtype: ", subtype)
  }

  .capValues <- function(data, valueCol, ceil = 0.30914) {

    cellsModified <- sum(data[[valueCol]] > ceil, na.rm = TRUE)

    idx <- data[[valueCol]] > ceil
    data[[valueCol]][idx] <- ceil

    message(cellsModified, " cells were capped at ", ceil)

    return(data)
  }

  dat <- read.csv(
    file = files[subtype],
    stringsAsFactors = FALSE
  )

  # Separate column "x.y.iso.t.kcr.value"
  parts <- strsplit(dat[["x.y.iso.t.kcr.value"]], ";", fixed = TRUE)
  parts <- do.call(rbind, parts)

  colnames(parts) <- c("x.y.iso", "t", "kcr", "value")

  dat <- cbind(
    dat[, setdiff(names(dat), "x.y.iso.t.kcr.value"), drop = FALSE],
    as.data.frame(parts, stringsAsFactors = FALSE)
  )

  dat[["t"]]     <- as.integer(dat[["t"]])
  dat[["value"]] <- as.numeric(dat[["value"]])

  # Conversion (ha to Mha)
  dat[["value"]] <- dat[["value"]] / 1e6

  # Reverse rows to standardize the x.y.iso order with LUH3 data
  #dat <- dat[rev(seq_len(nrow(dat))), ]

  # Cap values
  dat <- .capValues(dat, "value")

  mag <- magclass::as.magpie(
    dat,
    spatial = "x.y.iso",
    tidy = TRUE
  )

  dimnames(mag)[[1]] <- unique(dat[["x.y.iso"]])
  #magclass::getItems(mag, 1, raw = TRUE) <- dat[["x.y.iso"]]

  sum(mag, na.rm = TRUE)

  return(mag)
}
