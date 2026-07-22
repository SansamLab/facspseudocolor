test_that("current EdU CSV workflow preserves fitted models and normalization", {
  cfg <- read_example_config("config_edu.yml")
  cfg$data_dir <- example_data_dir("example")
  settings <- current_settings_from_config(cfg)
  manifest <- make_sample_manifest(replicates = cfg$replicates)

  models <- fit_replicate_reference_models(
    manifest, cfg$data_dir, cfg$suffixes, settings
  )
  results <- lapply(seq_len(nrow(manifest)), function(i) {
    read_and_normalize_sample(
      prefix = manifest$prefix[[i]],
      condition_label = manifest$condition[[i]],
      data_dir = cfg$data_dir,
      file_suffixes = cfg$suffixes,
      settings = settings,
      baseline_slope = models[[as.character(manifest$replicate_index[[i]])]]$slope
    )
  })

  expect_equal(
    vapply(models, function(x) x$slope, numeric(1)),
    c(`1` = 6.63379258818644, `2` = 4.66104681005732),
    tolerance = 1e-8
  )
  expect_equal(
    vapply(models, function(x) x$negative_intercept, numeric(1)),
    c(`1` = 2499.02152598892, `2` = 3475.09405621298),
    tolerance = 1e-8
  )
  expect_equal(
    vapply(models, function(x) x$negative_event_n, numeric(1)),
    c(`1` = 5564, `2` = 5967)
  )
  expect_equal(
    vapply(results, function(x) stats::median(x$data$target_norm, na.rm = TRUE),
           numeric(1)),
    c(1752.82011390766, 1481.56260982336, 1353.90517673248,
      2349.03387094802, 1664.84907317370, 1587.49882678020),
    tolerance = 1e-8
  )
  expect_identical(
    vapply(results, function(x) nrow(x$data), integer(1)),
    c(18724L, 18412L, 18546L, 18223L, 18569L, 17902L)
  )
})

test_that("current POI CSV workflow preserves fitted models and normalization", {
  cfg <- read_example_config("config_poi.yml")
  cfg$data_dir <- example_data_dir("example_poi")
  settings <- current_settings_from_config(cfg)
  manifest <- make_sample_manifest(replicates = cfg$replicates)

  models <- fit_replicate_background_models(
    manifest, cfg$data_dir, cfg$suffixes, settings
  )
  results <- lapply(seq_len(nrow(manifest)), function(i) {
    read_and_normalize_sample(
      prefix = manifest$prefix[[i]],
      condition_label = manifest$condition[[i]],
      data_dir = cfg$data_dir,
      file_suffixes = cfg$suffixes,
      settings = settings,
      background_model = models[[as.character(manifest$replicate_index[[i]])]]
    )
  })

  expect_equal(
    vapply(models, function(x) x$intercept, numeric(1)),
    c(`1` = 1047.32666698722, `2` = 774.471642274585),
    tolerance = 1e-8
  )
  expect_equal(
    vapply(models, function(x) x$slope, numeric(1)),
    c(`1` = 7.89440695503372, `2` = 7.98133225253185),
    tolerance = 1e-8
  )
  expect_equal(
    vapply(models, function(x) x$cutoff, numeric(1)),
    c(`1` = 1292.26160094065, `2` = 1286.36360141968),
    tolerance = 1e-8
  )
  expect_equal(
    vapply(results, function(x) stats::median(x$data$target_norm, na.rm = TRUE),
           numeric(1)),
    c(956.742596862670, 4991.65220845589, 1638.55373071234,
      1826.48029208043, 968.992295650424, 4098.33453848040,
      1648.79195919279, 1689.28981136966),
    tolerance = 1e-8
  )
  expect_identical(
    vapply(results, function(x) nrow(x$data), integer(1)),
    c(9837L, 9818L, 9725L, 14503L, 9461L, 8792L, 9607L, 9619L)
  )
})
