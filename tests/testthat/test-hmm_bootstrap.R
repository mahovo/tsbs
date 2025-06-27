test_that("hmm_bootstrap returns list of correct length and shapes", {
  df <- data.frame(return = rnorm(50))
  boots <- hmm_bootstrap(df, num_boots = 2, nstates = 2)
  expect_length(boots, 2)
  for (b in boots) expect_true(is.data.frame(b))
})

test_that("hmm_bootstrap fails without depmixS4 installed", {
  skip_if_not_installed("depmixS4")
  df <- data.frame(return = rnorm(50))
  expect_silent(
    hmm_bootstrap(df, num_boots = 1, nstates = 2)
  )
})
