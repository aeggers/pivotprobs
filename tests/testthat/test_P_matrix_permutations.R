context("Permuting the P matrix")

# I only do this for MC, EV and EN here because it's a common structure for others and they take a long time.

test_that("main function yields correct P matrices for plurality via EV", {

  plurality_result <- plurality_election(k = 4) %>%
    election_event_probs(method = "ev", alpha = c(12, 8, 5, 2.5))

  expect_true(is.null(plurality_result[["a_"]]))
  expect_equal(plurality_result[["a_b"]]$P, rbind(c(1,0,1,1), c(0,1,0,0), 0, 0))
  expect_equal(plurality_result[["b_a"]]$P, rbind(c(1,0,0,0), c(0,1,1,1), 0, 0))
  expect_equal(plurality_result[["c_a"]]$P, rbind(c(1,0,0,0), 0, c(0,1,1,1), 0))

})

test_that("main function yields correct P matrices for plurality via MC", {

  plurality_result <- plurality_election(k = 4) %>%
    election_event_probs(method = "mc", alpha = c(12, 8, 5, 2.5), num_sims = 10001)

  expect_equal(plurality_result[["a_"]]$P, rbind(rep(1, 4), 0, 0, 0))
  expect_equal(plurality_result[["a_b"]]$P, rbind(c(1,0,1,1), c(0,1,0,0), 0, 0))
  expect_equal(plurality_result[["b_a"]]$P, rbind(c(1,0,0,0), c(0,1,1,1), 0, 0))
  expect_equal(plurality_result[["c_a"]]$P, rbind(c(1,0,0,0), 0, c(0,1,1,1), 0))

})


test_that("main function yields correct P matrices for irv second round pivot events via EN", {

  irv_result <- irv_election() %>%
    election_event_probs(method = "en", alpha = c(12, 8, 5, 2.5, 7, 8))

  expect_equal(irv_result[["a_b"]]$P, rbind(c(1,1,0,0,1,0), c(0,0,1,1,0,1), 0))
  expect_equal(irv_result[["b_a"]]$P, rbind(c(1,1,0,0,1,0), c(0,0,1,1,0,1), 0))
  expect_equal(irv_result[["c_a"]]$P, rbind(c(1,1,1,0,0,0), 0, c(0,0,0,1,1,1)))

})

test_that("main function yields correct P matrices for irv first round pivot events via EN", {

  irv_result <- irv_election() %>%
    election_event_probs(method = "en", alpha = c(12, 8, 5, 2.5, 7, 8))

  expect_equal(irv_result[["a_b|ab"]]$P, rbind(c(1,1,0,0,1,1), c(0,0,1,1,0,0), 0))
  expect_equal(irv_result[["b_a|ba"]]$P, rbind(c(1,1,0,0,0,0), c(0,0,1,1,1,1), 0))
  expect_equal(irv_result[["c_a|ba"]]$P, rbind(c(1,1,0,0,0,0), c(0,0,1,1,1,1), 0))
  expect_equal(irv_result[["c_b|ab"]]$P, rbind(c(1,1,0,0,1,1), c(0,0,1,1,0,0), 0))
  expect_equal(irv_result[["c_a|ca"]]$P, rbind(c(1,1,0,0,0,0), 0, c(0,0,1,1,1,1)))
  expect_equal(irv_result[["c_a|cb"]]$P, rbind(0, c(1,1,0,0,0,0), c(0,0,1,1,1,1)))

})

test_that("main function yields correct P matrices for all pivot events via MC", {

  irv_result <- irv_election() %>%
    election_event_probs(method = "mc", alpha = c(12, 8, 5, 2.5, 7, 8), num_sims = 10001)

  expect_equal(irv_result[["a__b"]]$P, rbind(c(1,1,1,1,1,1), 0, 0))
  expect_equal(irv_result[["b__a"]]$P, rbind(0, c(1,1,1,1,1,1), 0))
  expect_equal(irv_result[["a_b"]]$P, rbind(c(1,1,0,0,1,0), c(0,0,1,1,0,1), 0))
  expect_equal(irv_result[["b_a"]]$P, rbind(c(1,1,0,0,1,0), c(0,0,1,1,0,1), 0))
  expect_equal(irv_result[["c_a"]]$P, rbind(c(1,1,1,0,0,0), 0, c(0,0,0,1,1,1)))
  expect_equal(irv_result[["a_b|ab"]]$P, rbind(c(1,1,0,0,1,1), c(0,0,1,1,0,0), 0))
  expect_equal(irv_result[["b_a|ba"]]$P, rbind(c(1,1,0,0,0,0), c(0,0,1,1,1,1), 0))
  expect_equal(irv_result[["c_a|ba"]]$P, rbind(c(1,1,0,0,0,0), c(0,0,1,1,1,1), 0))
  expect_equal(irv_result[["c_b|ab"]]$P, rbind(c(1,1,0,0,1,1), c(0,0,1,1,0,0), 0))
  expect_equal(irv_result[["c_a|ca"]]$P, rbind(c(1,1,0,0,0,0), 0, c(0,0,1,1,1,1)))
  expect_equal(irv_result[["c_a|cb"]]$P, rbind(0, c(1,1,0,0,0,0), c(0,0,1,1,1,1)))

})



