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

  # as.magpie interprets dots as sub-dimension separators (e.g. "c3ann.irrigated" splits into
  # dim 3.1="c3ann" / dim 3.2="irrigated" and collapses back as "c3ann_irrigated").
  # Use a placeholder to preserve literal dots in class names.
  dat[["data"]] <- gsub(".", "DOTPLACEHOLDER", dat[["data"]], fixed = TRUE)
  mag <- magclass::as.magpie(dat, spatial = "x.y.iso", tidy = TRUE)
  dimnames(mag)[[1]] <- unique(dat[["x.y.iso"]])
  getItems(mag, 3) <- gsub("DOTPLACEHOLDER", ".", getItems(mag, 3), fixed = TRUE)

  return(mag)
}
