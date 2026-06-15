#' @title readMapBiomas
#' @description
#' Reads MapBiomas land cover and irrigation data for Brazil and converts
#' into magpie objects. Data was produced from MapBiomas Collection 10
#' combined with IBGE planted area statistics, covering 12 LUH3 land use
#' classes (land cover) and 6 irrigation classes.
#'
#' @param subtype Data subtype: "LandCover" or "Irrigation"
#'
#' @return A magpie object with land use data for Brazil (Mha)
#' @author Wanderson Costa, Alexandre Koberle, Miodrag Stevanovic
#' @examples
#' \dontrun{
#' readSource("MapBiomas", subtype = "LandCover")
#' readSource("MapBiomas", subtype = "Irrigation")
#' }
#'
#' @importFrom utils read.csv

readMapBiomas <- function(subtype = "LandCover") {

  files <- c(
    LandCover  = "mapbiomas_luh3_landcover.csv",
    Irrigation = "mapbiomas_irrigation.csv"
  )

  if (!subtype %in% names(files)) {
    stop("Unknown subtype: ", subtype, ". Available: ", paste(names(files), collapse = ", "))
  }

  dat <- read.csv(
    file             = files[[subtype]],
    sep              = ";",
    stringsAsFactors = FALSE
  )

  names(dat)[names(dat) == "year"] <- "t"

  # ha to Mha
  dat[["value"]] <- dat[["value"]] / 1e6

  # explicit column order required by as.magpie tidy = TRUE
  dat <- dat[, c("x.y.iso", "t", "data", "value")]

  # as.magpie converts "." to "p" in class names (e.g. "c3ann.irrigated" → "c3annpirrigated")
  # save original names and expected transformed names before conversion
  originalClasses    <- unique(dat[["data"]])
  transformedClasses <- gsub(".", "p", originalClasses, fixed = TRUE)

  mag <- magclass::as.magpie(dat, spatial = "x.y.iso", tidy = TRUE)
  dimnames(mag)[[1]] <- unique(dat[["x.y.iso"]])

  # verify order and transformation match before restoring — catches any as.magpie behaviour change
  stopifnot(identical(getItems(mag, 3), transformedClasses))
  getItems(mag, 3) <- originalClasses

  return(mag)
}
