test_that("min_set_cover works", {
  x <- matrix(c(1,1,1,1,1,1,1,0,0,1,0,0,1,1,0,1), nrow = 4)
  y <- c(1,1,1,0)

  e <- 1
  xe <- x[e,]
  ye <- y[e]

  r <- min_set_cover(x, y, xe, ye, complement = FALSE)

  expect_equal(r$R, 1)
  expect_equal(r$alpha, c(0,1,0,0))

  x2 <- rbind(x, c(0,1,1,1))
  y2 <- c(y,0)
  r2 <- min_set_cover(x2, y2, xe, ye, complement = FALSE)

  expect_equal(r2$R, 2)
  expect_equal(r2$alpha, c(1,1,0,0))
})
