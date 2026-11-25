# Author: Davide Delfino
# Input: none.
# Output: replication configuration, shared across models.

replication_config <- list(
  countries = c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT"),
  aggregate = "EA",
  start_est = as.Date("2000-04-01"),
  start_eval = as.Date("2017-01-01"),
  end_eval = as.Date("2026-02-01"),
  covid_start = as.Date("2020-03-01"),
  covid_end = as.Date("2021-07-01"),
  information_set_size = 50L,
  initial_monthly = 25L,
  initial_quarterly = 5L,
  selection_end = as.Date("2026-02-01"),
  dfm_idiosyncratic = "q0",
  monte_carlo_draws = 500L,
  an_draws = 2000L,
  seed = 1712L
)
