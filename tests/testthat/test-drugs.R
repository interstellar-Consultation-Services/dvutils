test_that("drug_descriptors_categories", {
  expect_equal(drug_descriptors_categories(),
               c("hybrid",
                 "constitutional",
                 "topological",
                 "electronic",
                 "geometrical"))
})

test_that("drug_descriptors_names", {
  expect_equal(drug_descriptors_names(drug_descriptors_categories()[1]),
               c("org.openscience.cdk.qsar.descriptors.molecular.WHIMDescriptor",
                 "org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor"))
  #Test asserts
  msg <- paste("category must be 'all' or value(s) from",
               "drug_descriptors_categories() returned values.",
               "Setting default value 'all'")
  expect_warning(names_len <- length(drug_descriptors_names(NULL)),
                 regexp = msg,
                 fixed = TRUE)
  expect_equal(names_len, 51)

  expect_warning(names_len <- length(drug_descriptors_names(NA)),
                 regexp = msg,
                 fixed = TRUE)
  expect_equal(names_len, 51)

  expect_warning(names_len <- length(drug_descriptors_names(123)),
                 regexp = msg,
                 fixed = TRUE)
  expect_equal(names_len, 51)

  expect_warning(names_len <- length(drug_descriptors_names("abc")),
                 regexp = msg,
                 fixed = TRUE)
  expect_equal(names_len, 51)

  expect_warning(names_len <- length(drug_descriptors_names(c("abc", "electronic"))),
                 regexp = msg,
                 fixed = TRUE)
  expect_equal(names_len, 51)

})

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
