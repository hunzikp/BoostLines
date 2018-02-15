Rcpp::loadModule("mod", TRUE)

#' Create a BoostLines Object
#'
#' Creates a BoostLines object.
#'
#' BoostLine objects may not contain multilines.
#'
#' @param x A SpatialLines object or a Rcpp_boost_line_collection object.
#' @param proj4string A proj4string, ignored if \code{x} is SpatialLines.
#'
#' @return A BoostLine S3 object.
#'
#' @import Rcpp
#' @import sp
#'
#' @export
BoostLines <- function(x, proj4string=CRS('')) {

  if (inherits(x, 'SpatialLines')) {
    sl <- x

    ## Make blc object
    coords_ls <- lapply(sl@lines, {
      function(x) {
        if (length(x@Lines) > 1) {
          stop('No multilines!')
        }
        x@Lines[[1]]@coords
      }
    })
    blc <- make_blc(coords_ls)

    ## Make and return class
    me <- list(
      blc = blc,
      proj4string = sl@proj4string
    )
  } else if (inherits(x, 'Rcpp_boost_line_collection')) {

    me <- list(
      blc = x,
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
#' @import sp
#'
#' @export
unboost <- function(bl) {
  UseMethod("unboost", bl)
}

unboost.BoostLines <- function(bl) {

  coords_ls <- unpack_blc(bl$blc)
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
#' @return A list of (1) a new, noded BoostLines object and (2) an integer vector of line IDs.
#'
#' @export
bNode <- function(bl) {
  UseMethod("bNode", bl)
}

bNode.BoostLines <- function(bl) {

  broken_ls <- node_blc(bl$blc)
  broken_blc <- broken_ls[[1]]
  broken_ids <- broken_ls[[2]]
  broken_bl <- BoostLines(broken_blc, bl$proj4string)
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

  lm <- intersects_blc(bl$blc)

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

  dm <- distance_blc(bl$blc)

  return(dm)
}

#' Calculate BoostLines Lengths
#'
#' Calculates line lengths.
#'
#' @param bl A BoostLines object.
#' @param lonlat A boolean; are line coordinates in lon/lat? If TRUE, returns lenghts in meters.
#'
#' @return A numeric vector.
#'
#' @export
bLength <- function(bl, lonlat) {
  UseMethod("bLength", bl)
}

bLength.BoostLines <- function(bl, lonlat) {

  lengths <- get_lengths(bl$blc, lonlat)

  return(lengths)
}

#' Create a Gridded Line Collection
#'
#' Creates a BoostLines object with gridded lines.
#'
#'
#' @param x Upper left coordinate pair OR a RasterLayer object OR a VeloxRaster object.
#' @param nrow Number of rows. Required if x is a coordinate pair.
#' @param ncol Number of columns. Required if x is a coordinate pair.
#' @param res Resolution (x and y). Required if x is a coordinate pair.
#' @param n8 Boolean. Whether 8 neighborhood should be connected. If FALSE, connects 4 neighborhood.
#' @param proj4string A proj4string. Ignored if x is a RasterLayer OR a VeloxRaster object.
#'
#' @return A BoostLine S3 object.
#'
#' @import Rcpp
#' @import sp
#'
#' @export
bGrid <- function(x, nrow=NULL, ncol=NULL, res=NULL, n8 = TRUE, proj4string = CRS('')) {

  if (inherits(x, what = 'numeric')) {
    if (length(x) == 2) {

      if (!is.null(nrow) & !is.null(ncol) & !is.null(res)) {
        x_upper_left <- x[1]
        y_upper_left <- x[2]
      } else {
        stop('Arguments nrow, ncol, and res are required if x is a coordinate pair.')
      }

    } else {
      stop('x is numeric but not a coordinate pair.')
    }
  } else if (inherits(x, what = 'RasterLayer')) {

    if (res(x)[1] != res(x)[2]) {
      warning('x and y resolution differ. Using x resolution')
    }
    res <- res(x)[1]
    nrow <- nrow(x)
    ncol <- ncol(x)
    x_upper_left <- xmin(x) + res/2
    y_upper_left <- ymax(x) - res/2
    proj4string <- CRS(proj4string(x))

  } else if (inherits(x, what = 'VeloxRaster')) {

    vxres <- x$res
    if (vxres(x)[1] != vxres(x)[2]) {
      warning('x and y resolution differ. Using x resolution')
    }

    res <- vxres[1]
    nrow <- x$dim[1]
    ncol <- x$dim[2]
    x_upper_left <- x$extent[1] + res/2
    y_upper_left <- x$extent[4] - res/2
    proj4string <- CRS(x$crs)

  } else {
    stop('x is an unsupported object.')
  }

  ## Make gridded boost_line_collection
  blc <- make_gridded_blc(x = x_upper_left, y = y_upper_left, nrow = nrow, ncol = ncol, res = res, n8 = n8)

  ## Make BoostLines object
  bl <- BoostLines(blc, proj4string)

  return(bl)
}


#' Concatenate a List of BoostLines Objects
#'
#' Creates a single BoostLines object from multiple BoostLines objects.
#'
#'
#' @param x A List of BoostLines objects.
#'
#' @return A BoostLine S3 object.
#'
#' @import Rcpp
#'
#' @export
bAppend <- function(x) {

  ## Sanity checks
  if (!inherits(x, 'list')) {
    stop('x must be a list of BoostLines objects.')
  }
  if (length(x) < 2) {
    stop('x must be of length > 1.')
  }
  bl.class <- unlist(lapply(x, function(x) inherits(x, 'BoostLines')))
  if (!all(bl.class)) {
    stop('x must be a list of BoostLines objects.')
  }
  proj.ls <- lapply(x, function(x) x$proj4string)
  if (length(unique(proj.ls)) > 1) {
    stop('BoostLine objects have different proj4strings.')
  }

  ## Append
  new_blc = x[[1]]$blc
  for (i in 2:length(x)) {
    new_blc = append_blc(new_blc, x[[i]]$blc)
  }
  new_bl <- BoostLines(new_blc, proj.ls[[1]])

  return(new_bl)
}


#' Eliminate Multilines
#'
#' Removes all multilines from a SpatialLinesDataFrame by replacing them with simple lines.
#' Attribute data for multilines is duplicated accordingly.
#'
#' @param s A SpatialLinesDataFrame or SpatialLines object.
#'
#' @return A SpatialLinesDataFrame or SpatialLines object.
#'
#' @import sp
#'
#' @export
remove_multilines <- function(s) {

  if (inherits(s, 'SpatialLines') | inherits(s, 'SpatialLinesDataFrame')) {
    lines.ls <- s@lines
    multiline.count <- unlist(lapply(lines.ls, function(x) length(x@Lines)))
    exploded_lines.id <- rep(1:length(lines.ls), times=multiline.count)
    exploded_lines.ls <- lapply(lines.ls, function(x) {
      this.Lines.ls <- x@Lines
      this.id <- x@ID
      cnt <- 0
      this.lines.ls <- lapply(this.Lines.ls, function(x) {
        cnt <<- cnt + 1;
        Lines(list(x), paste0(this.id, '.', cnt));
      })
      return(this.lines.ls);
    })
    exploded_lines.ls <- unlist(exploded_lines.ls)
    exploded.sl <- SpatialLines(exploded_lines.ls, proj4string = s@proj4string)

    if (inherits(s, 'SpatialLinesDataFrame')) {
      exploded.sldf <- SpatialLinesDataFrame(exploded.sl, s@data[exploded_lines.id,], FALSE)
      out <- exploded.sldf
    } else {
      out <- exploded.sl
    }

  } else {
    stop('s is not a supported object.')
  }

  return(out)
}

length.BoostLines <- function(bl) {

  l <- get_line_count(bl$blc)

  return(l)
}

#' @title
#' Fast boost2graph
#'
#' @description
#' Generates an igraph graph from a BoostLine object and a DataFrame.
#'
#' @details
#' Generates an igraph graph from a BoostLine object and a DataFrame.
#' Line attribute data is transformed into edge attribute data, and an extra 'length' edge attribute with
#' Euclidean or geodesic (see attribute \code{lonlat}) line lengths is added.
#'
#'
#' @param x A BoostLine object.
#' @param df A DataFrame with as many rows as there are lines in \code{x}.
#' @param lonlat A boolean indiciating whether BoostLines are in lat/lon format.
#' If TRUE, the length edge attribute added to the igraph object is measured in meters.
#' @param plot.result Should resulting graph be plotted?
#'
#' @return A \code{igraph::graph} object.
#'
#' @import igraph
#'
boost2graph <- function(x, df, lonlat=FALSE, plot.result=FALSE) {

  if (!inherits(x, 'BoostLines')) {
    stop('x is not a BoostLines object.')
  }
  if (nrow(df) != length(x)) {
    stop('df has an incorrect number of rows.')
  }

  ## Make vertex DF
  coord.ls <- get_line_coords(x$blc)
  endcoords.ls <- lapply(coord.ls, function(x) x[c(1,nrow(x)),])
  endcoord.mat <- unique(do.call("rbind", endcoords.ls))
  vertex.df <- data.frame(vid=1:nrow(endcoord.mat), x=endcoord.mat[,1], y=endcoord.mat[,2])

  ## Make edgelist DF
  endcoords.el.ls <- lapply(endcoords.ls, function(x) matrix(c(x[1,], x[2,]), ncol=4))
  el.mat <- do.call("rbind", endcoords.el.ls)
  el.df <- as.data.frame(el.mat)
  names(el.df) <- c('from.x', 'from.y', 'to.x', 'to.y')
  el.df$order <- 1:nrow(el.df)
  fromvertex.df <- vertex.df
  names(fromvertex.df) <- c('from.vid', 'from.x', 'from.y')
  el.df <- merge(el.df, fromvertex.df, all.x=TRUE, all.y=FALSE)
  tovertex <- vertex.df
  names(tovertex) <- c('to.vid', 'to.x', 'to.y')
  el.df <- merge(el.df, tovertex, all.x=TRUE, all.y=FALSE)
  el.df <- el.df[order(el.df$order),]
  el.df <- el.df[,c('from.vid', 'to.vid', 'order')]

  ## Add attribute data, euclidean length
  el.df <- cbind(el.df, df)
  if (!('length' %in% names(el.df))){
    el.df$length <- bLength(x, lonlat = lonlat)
  }

  ## Make graph
  graph <- graph_from_data_frame(el.df, FALSE, vertex.df)

  ## Optional plotting
  if (plot.result) {
    asp <- (max(endcoord.mat[,2])-min(endcoord.mat[,2]))/(max(endcoord.mat[,1])-min(endcoord.mat[,1]))
    plot(graph, vertex.size=0.25, vertex.label=NA, layout=as.matrix(vertex.df[,2:3]), asp=asp)
  }

  return(graph)
}


#' Split Lines
#'
#' Splits lines into shorter lines of length maxlen.
#'
#' Input lines are split into shorter lines of length maxlen, plus a remainder of length smaller than maxlength.
#' Input lines shorter than maxlen are left untouched.
#'
#'
#' @param bl A BoostLines object.
#' @param maxlen The maximum length of the split lines.
#'
#' @return A list of (1) a new, split BoostLines object and (2) an integer vector of line IDs.
#'
#' @export
bSplit <- function(bl, maxlen) {
  UseMethod("bSplit", bl)
}

bSplit.BoostLines <- function(bl, maxlen) {

  splitted_ls <- split_blc(bl$blc, maxlen)
  splitted_blc <- splitted_ls[[1]]
  splitted_ids <- splitted_ls[[2]]
  splitted_bl <- BoostLines(splitted_blc, bl$proj4string)
  out_ls <- list(splitted_bl, splitted_ids)

  return(out_ls)
}



