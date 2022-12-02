library(POMADE)

mdes_MADE(
  J = seq(10,100,5),
  tau = c(0.05,0.10,0.20),
  omega = 0,
  rho = c(0.5,0.7,0.9),
  target_power = .9,
  alpha = 0.05,
  sigma2_dist = \(x) rgamma(x, shape = 5, rate = 10),
  n_ES_dist = \(x) 1 + rpois(x, 5.5 - 1)
) |>
  plot_MADE(expected_studies = c(30, 50), numbers = FALSE)

ggsave("header.svg")
