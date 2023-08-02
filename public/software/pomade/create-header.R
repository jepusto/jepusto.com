library(POMADE)
library(ggplot2)
library(MetBrewer)
library(future)
library(furrr)

plan(multisession)

mdes_res <- 
  mdes_MADE(
    J = seq(10,100,5),
    tau = c(0.1,0.2,0.3),
    omega = 0.05,
    rho = c(0.5,0.7,0.9),
    target_power = .9,
    alpha = 0.05,
    sigma2_dist = \(x) rgamma(x, shape = 5, rate = 10),
    n_ES_dist = \(x) 1 + rpois(x, 5.5 - 1)
  )

plot_MADE(
  mdes_res, expected_studies = c(30, 50), numbers = FALSE,
  y_limits = c(0,0.90), y_breaks = seq(0,0.9,0.1),
  x_breaks = seq(10,100,10)
) + 
  scale_color_manual(values = met.brewer("Demuth", n = 3))

ggsave("static/img/headers/MDES-MADE.png", width = 10, height = 5)

mdes_res |>
  subset(tau == 0.2) |>
plot_MADE(
  expected_studies = c(30, 50), numbers = FALSE,
  y_limits = c(0,0.90), y_breaks = seq(0,0.9,0.1),
  x_breaks = seq(10,100,10)
) + 
  scale_color_manual(values = met.brewer("Demuth", n = 3))

ggsave("content/software/POMADE/featured.png", width = 5, height = 5)
