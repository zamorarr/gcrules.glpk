#' Minimum Set Cover for a Binary Matrix
#'
#' @param x binary matrix of features
#' @param y binary vector of predictor
#' @param xe binary vector of example observation
#' @param ye binary value of example predictor
#' @param complement whether to add the complement of x as features
#' @export
#' @importClassesFrom Matrix TsparseMatrix
min_set_cover <- function(x, y, xe, ye, complement = TRUE) {
  # add complement variables
  if (complement) {
    x <- cbind(x, 1L - x)
    xe <- c(xe, 1L - xe)
  }

  # check dimensions
  n <- nrow(x)
  d <- ncol(x)
  stopifnot(length(y) == n)
  stopifnot(length(xe) == d)
  stopifnot(length(ye) == 1)

  # setup logfile
  #logfile <- tempfile(pattern = "gcr-", fileext = ".log")
  #logfile <- normalizePath(logfile, mustWork = FALSE)
  #cat(sprintf("writing solver log to %s\n", logfile))

  # build constraint matrix
  idx <- which(y != ye)
  s <- length(idx)

  # initialize empty contraint matrix
  m <- matrix(0L, nrow = 2 + s, ncol = 1 + d)

  # add "consistent" constraints
  m[1:s,1:d] <- 1L - x[idx,]

  # add "relevance" constraint
  m[1 + s, 1:d] <- xe
  m[1 + s, 1 + d] <- -1L

  # add auxillary variable constraint
  m[2 + s, 1:d] <- 1L
  m[2 + s, 1 + d] <- -1L


  # convert contraint matrix to sparse i,j,v triplet
  m <- methods::as(m, "TsparseMatrix")
  mi <- m@i + 1L
  mj <- m@j + 1L
  mv <- as.integer(m@x)

  # run min_set_cover algorithm
  r <- min_set_cover_cpp(mi, mj, mv, n, d, s)

  # add some elements to return variable
  alpha_idx <- which(r$alpha > 0)
  if (complement) {
    alpha_idx <- ifelse(alpha_idx > d/2, alpha_idx - d/2, alpha_idx)
  }
  r$alpha_idx <- alpha_idx
  r
}
