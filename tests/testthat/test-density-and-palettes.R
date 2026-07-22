test_that("named and custom palettes resolve correctly", {
  expect_identical(resolve_palette("refined"), refined_density_palette())
  expect_identical(resolve_palette("flowjo"), flowjo_palette())
  expect_identical(resolve_palette(c("#000000", "#FFFFFF")), c("#000000", "#FFFFFF"))
  expect_error(resolve_palette("unknown"), "Unknown palette")
})

test_that("point density preserves missing positions", {
  set.seed(1)
  x <- c(stats::rnorm(20), NA_real_)
  y <- c(stats::rnorm(20), 1)
  density <- compute_point_density(x, y, n = 40)
  expect_length(density, length(x))
  expect_true(all(is.finite(density[seq_len(20)])))
  expect_true(is.na(density[[21]]))
  expect_error(compute_point_density(1:9, 1:9), "Too few valid events")
})

test_that("density colors are clipped and gamma-scaled to zero through one", {
  color <- prepare_density_color(1:100, 0.1, 0.9, gamma = 1)
  expect_equal(min(color), 0)
  expect_equal(max(color), 1)
  expect_true(all(color >= 0 & color <= 1))
  expect_error(prepare_density_color(1:10, 0.9, 0.1), "must satisfy")
  expect_error(prepare_density_color(1:10, gamma = 0), "positive finite")
})
