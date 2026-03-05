#' @title readIBGE
#' @description
#' Reads IBGE cropland (planted area) time series data for Brazil and converts
#' it into a magpie object.
#' Data was collected from the SIDRA IBGE Platform
#' (https://sidra.ibge.gov.br//).
#'
#' @param subtype Data subtype (currently only "Cropland")
#' @param yrs     years to be returned
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

  yrs <- seq(1965, 2020, 5)

  files <- c(
    Cropland = "crop_planted_area_1995_to_2024_luh3_all.csv"
  )

  if (!subtype %in% names(files)) {
    stop("Unknown subtype: ", subtype)
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

  # If years before 1995 are requested, replicate 1995 values
  earlyYears <- yrs[yrs < 1995]

  if (length(earlyYears) > 0) {
    base1995 <- dat[dat$t == 1995, ]

    extra <- do.call(rbind, lapply(earlyYears, function(y) {
      tmp <- base1995
      tmp$t <- y
      tmp
    }))

    dat <- rbind(dat, extra)
  }

  # Filter selected years
  dat <- dat[dat$t %in% yrs, ]

  # Conversion (ha to Mha)
  dat[["value"]] <- dat[["value"]] / 1e6

  dat <- dat[order(dat$x.y.iso, dat$t, dat$kcr), ]

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