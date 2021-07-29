#' Find Max Support for a Binary Matrix
#'
#' @param x binary matrix of features
#' @param y binary vector of predictor
#' @param xe binary vector of example observation
#' @param ye binary value of example predictor
#' @param complement whether to add the complement of x as features
#' @export
#' @importClassesFrom Matrix TsparseMatrix
max_support <- function(x, y, xe, ye, R_max = ncol(x), ws = 1, complement = TRUE) {
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

  # build constraint matrix
  idx <- which(y != ye)
  n_neg <- length(idx)
  n_pos <- length(which(y == ye))

  # initialize empty contraint matrix
  m <- matrix(0L, nrow = 4 + n_neg + n, ncol = 2 + d + n)

  # add "consistent" constraints
  m[1:n_neg,1:d] <- 1L - x[idx,]

  # add "rule" constraints
  m[n_neg + (1:n), 1:d] <- 1 - x[1:n,1:d]
  for (i in 1:n) {
    m[n_neg + i, d + i] <- 1L + d
  }
  #m[s + (1:n), d + (1:n)] <- 1L + d

  # add "relevance" constraint
  m[1 + n_neg + n, 1:d] <- xe
  m[1 + n_neg + n, 1 + d + n] <- -1L

  # add auxillary (sum(alpha) == R)
  m[2 + n_neg + n, 1:d] <- 1L
  m[2 + n_neg + n, 1 + d + n] <- -1L

  # max sparsity
  m[3 + n_neg + n, 1 + d + n] <- 1L

  # add auxillary (sum(r) == S)
  m[4 + n_neg + n, d + (1:n)] <- 1L
  m[4 + n_neg + n, 2 + d + n] <- -1L

  #return(m)

  # convert contraint matrix to sparse i,j,v triplet
  m <- methods::as(m, "TsparseMatrix")
  mi <- m@i + 1L
  mj <- m@j + 1L
  mv <- as.integer(m@x)

  # run min_set_cover algorithm
  r <- max_support_cpp(mi, mj, mv, n, d, n_neg, R_max, ws)

  # remove complement variables
  if (complement) {
    alpha <- r$alpha
    r$alpha <- ifelse(alpha > d/2, alpha - d/2, alpha)
  }

  r
}
