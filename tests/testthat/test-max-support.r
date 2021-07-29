test_that("max_support works", {
  x <- matrix(c(1,1,1,1,1,1,1,0,0,1,0,0,1,1,0,1), nrow = 4)
  y <- c(1,1,1,0)

  e <- 1
  xe <- x[e,]
  ye <- y[e]

  r <- max_support(x, y, xe, ye, complement = FALSE, ws = 1)

  expect_equal(r$sparsity, 1)
  expect_equal(r$alpha, 2)
  expect_equal(r$support, 3)
  expect_equal(r$r, c(1,2,3))

  x2 <- rbind(x, c(0,1,1,1))
  y2 <- c(y,0)
  r2 <- max_support(x2, y2, xe, ye, complement = FALSE, ws = 1)

  expect_equal(r2$sparsity, 2)
  expect_equal(r2$alpha, c(1,2))
  expect_equal(r$support, 3)
  expect_equal(r$r, c(1,2,3))
})
