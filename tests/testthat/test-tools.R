test_that("build_adjacency_matrix", {
  # test asserts
  ## test data
  msg  <-  "data must be a dataframe that contains at least two columns"
  expect_error(build_adjacency_matrix(),
               regexp = msg)
  expect_error(build_adjacency_matrix(NULL),
               regexp = msg)
  expect_error(build_adjacency_matrix(NA),
               regexp = msg)
  expect_error(build_adjacency_matrix("data"),
               regexp = msg)
  expect_error(build_adjacency_matrix(12),
               regexp = msg)
  expect_error(build_adjacency_matrix(drugs_targets_matrix[, 1]),
               regexp = msg)

  ## test source
  msg <- "source must be a character value reprsents a valid data column name"
  expect_error(build_adjacency_matrix(drugs_targets_matrix),
               regexp = msg)
  expect_error(build_adjacency_matrix(drugs_targets_matrix, NULL),
               regexp = msg)
  expect_error(build_adjacency_matrix(drugs_targets_matrix, NA),
               regexp = msg)
  expect_error(build_adjacency_matrix(drugs_targets_matrix, "Not_column"),
               regexp = msg)
  expect_error(build_adjacency_matrix(drugs_targets_matrix, 123),
               regexp = msg)
  expect_error(build_adjacency_matrix(drugs_targets_matrix, c("parent_key",
                                                              "ID")),
               regexp = msg)

  ## test target
  msg <- "target must be a character value reprsents a valid data column name"
  expect_error(build_adjacency_matrix(drugs_targets_matrix,
                                      source = "parent_key"),
               regexp = msg)
  expect_error(build_adjacency_matrix(drugs_targets_matrix,
                                      source = "parent_key",
                                      NULL),
               regexp = msg)
  expect_error(build_adjacency_matrix(drugs_targets_matrix,
                                      source = "parent_key",
                                      NA),
               regexp = msg)
  expect_error(build_adjacency_matrix(drugs_targets_matrix,
                                      source = "parent_key",
                                      "Not_column"),
               regexp = msg)
  expect_error(build_adjacency_matrix(drugs_targets_matrix,
                                      source = "parent_key",
                                      123),
               regexp = msg)
  expect_error(build_adjacency_matrix(drugs_targets_matrix,
                                      source = "parent_key",
                                      c("parent_key", "ID")),
               regexp = msg)
})
