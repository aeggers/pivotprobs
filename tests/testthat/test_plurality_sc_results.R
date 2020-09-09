context("Plurality results via SC")

test_that("Pivot probs in neutral case are equal across candidates", {

  plurality_sc <- plurality_event_list(k = 3) %>%
    event_probabilities_from_event_list(method = "sc", alpha = c(5,5,5), tol = .1)

  expect_equal(plurality_sc[["a_"]]$integral, 1/3, tolerance = .01)
  expect_equal(plurality_sc[["b_"]]$integral, 1/3, tolerance = .01)
  expect_equal(plurality_sc[["c_"]]$integral, 1/3, tolerance = .01)

  expect_equal(plurality_sc[["a_b"]]$integral, plurality_sc[["b_a"]]$integral, tolerance = .00001) # will actually be exactly the same, because same calculation

  expect_equal(plurality_sc[["a_bc"]]$integral, plurality_sc[["b_ac"]]$integral, tolerance = .00001) # will actually be exactly the same, because same calculation

})

test_that("Pivot probs in non-neutral case have reasonable ordering", {

  plurality_sc <- plurality_event_list(k = 3) %>%
    event_probabilities_from_event_list(method = "sc", alpha = c(10,7,4), tol = .1)

  expect_true(plurality_sc[["a_"]]$integral > plurality_sc[["b_"]]$integral)
  expect_true(plurality_sc[["b_"]]$integral > plurality_sc[["c_"]]$integral)
  expect_true(plurality_sc[["a_b"]]$integral > plurality_sc[["a_c"]]$integral)
  expect_true(plurality_sc[["b_a"]]$integral > plurality_sc[["b_c"]]$integral)

})

test_that("Pivot probs in non-neutral case have reasonable ordering when drop_dimension", {

  plurality_sc <- plurality_event_list(k = 3) %>%
    event_probabilities_from_event_list(method = "sc", drop_dimension = T, alpha = c(10,7,4), tol = .1)

  expect_true(plurality_sc[["a_"]]$integral > plurality_sc[["b_"]]$integral)
  expect_true(plurality_sc[["b_"]]$integral > plurality_sc[["c_"]]$integral)
  expect_true(plurality_sc[["a_b"]]$integral > plurality_sc[["a_c"]]$integral)
  expect_true(plurality_sc[["b_a"]]$integral > plurality_sc[["b_c"]]$integral)

})


test_that("Pivot probs in non-neutral case have reasonable ordering when drop_dimension and skip_non_pivot_events", {

  plurality_sc <- plurality_event_list(k = 3) %>%
    event_probabilities_from_event_list(method = "sc", drop_dimension = T, skip_non_pivot_events = T, alpha = c(10,7,4), tol = .01)

  expect_true(is.null(plurality_sc[["a_"]]))
  expect_true(plurality_sc[["a_b"]]$integral > plurality_sc[["a_c"]]$integral)
  expect_true(plurality_sc[["b_a"]]$integral > plurality_sc[["b_c"]]$integral)
  # not sure if we have enough precision for this
  expect_true(plurality_sc[["a_b"]]$integral > plurality_sc[["b_a"]]$integral)

})


test_that("Adjacent pivot probs in non-neutral case are equal when merge_adjacent_pivot_events", {

  plurality_sc <- plurality_event_list(k = 3) %>%
    event_probabilities_from_event_list(method = "sc", drop_dimension = T, skip_non_pivot_events = T, merge_adjacent_pivot_events = T, alpha = c(10,7,4), tol = .01)

  expect_true(is.null(plurality_sc[["a_"]]))
  expect_equal(plurality_sc[["a_b"]]$integral, plurality_sc[["a_b"]]$integral)
  expect_equal(plurality_sc[["a_c"]]$integral, plurality_sc[["c_a"]]$integral)
})
