Rcpp::loadModule("mod", TRUE)

#' Create a BoostLines Object
#'
#' Creates a BoostLines object.
#'
#' BoostLine objects may not contain multilines.
#'
#' @param x A SpatialLines object or a Rcpp_line_collection object.
#' @param proj4string A proj4string, ignored if \code{x} is SpatialLines.
#'
#' @return A BoostLine S3 object.
#'
#' @import Rcpp
#'
#' @export
BoostLines <- function(x, proj4string='') {

  if (inherits(x, 'SpatialLines')) {
    sl <- x

    ## Make lc object
    coords_ls <- lapply(sl@lines, {
      function(x) {
        if (length(x@Lines) > 1) {
          stop('No multilines!')
        }
        x@Lines[[1]]@coords
      }
    })
    lc <- make_lc(coords_ls)

    ## Make and return class
    me <- list(
      lc = lc,
      proj4string = sl@proj4string
    )
  } else if (inherits(x, 'Rcpp_line_collection')) {

    me <- list(
      lc = x,
      proj4string = proj4string
    )
  } else {
    stop('x of invalid class.')
  }

  class(me) <- append(class(me), "BoostLines")

  return(me)
}


#' Cast BoostLines as SpatialLines
#'
#' Creates a SpatialLines object.
#'
#' @param bl A BoostLines object.
#'
#' @return A SpatialLines object.
#'
#' @export
unboost <- function(bl) {
  UseMethod("unboost", bl)
}

unboost.BoostLines <- function(bl) {

  coords_ls <- unpack_lc(bl$lc)
  cntr <- 1
  lns_ls <- lapply(coords_ls, {
    function(x) {
      ln <- Lines(slinelist = list(Line(coords = x)), ID = as.character(cntr))
      cntr <<- cntr + 1
      return(ln)
    }
  })
  sl <- SpatialLines(lns_ls, proj4string = bl$proj4string)

  return(sl)
}

#' Node BoostLines object
#'
#' Breaks lines at all proper intersections.
#'
#' @param bl A BoostLines object.
#'
#' @return A list of (1) the noded BoostLines object and (2) an integer vector of line IDs.
#'
#' @export
bNode <- function(bl) {
  UseMethod("bNode", bl)
}

bNode.BoostLines <- function(bl) {

  broken_ls <- node_lc(bl$lc)
  broken_lc <- broken_ls[[1]]
  broken_ids <- broken_ls[[2]]
  broken_bl <- BoostLines(broken_lc, bl$proj4string)
  out_ls <- list(broken_bl, broken_ids)

  return(out_ls)
}

#' Test BoostLines Intersections
#'
#' Checks whether BoostLines intersect.
#'
#' @param bl A BoostLines object.
#'
#' @return A logical matrix.
#'
#' @export
bIntersects <- function(bl) {
  UseMethod("bIntersects", bl)
}

bIntersects.BoostLines <- function(bl) {

  lm <- intersects_lc(bl$lc)

  return(lm)
}

#' Calculate BoostLines Minimal Distances
#'
#' Calculates minimal distances between BoostLines.
#'
#' @param bl A BoostLines object.
#'
#' @return A numeric matrix.
#'
#' @export
bDistance <- function(bl) {
  UseMethod("bDistance", bl)
}

bDistance.BoostLines <- function(bl) {

  dm <- distance_lc(bl$lc)

  return(dm)
}
