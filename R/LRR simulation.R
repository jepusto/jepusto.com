library(tidyverse)

# simulate AR(1) poisson series

r_pois_AR1 <- function(mu_vec, ar_vec, 
                       len = length(mu_vec), 
                       lambda = pmax(0, mu_vec[-1] - ar_vec * mu_vec[-len])) {
  
  # binomial thinning
  y <- vector(mode = "numeric", length = len)
  y[1] <- rpois(1, lambda = mu_vec[1])
  for (i in 1:(len-1)) 
    y[i+1] <- rbinom(1, size = y[i], prob = ar_vec[i]) + rpois(1, lambda = lambda[i])
  
  y
}

# estimate LRR three ways

calc_LRRs <- function(x, y) {

  n <- table(x)
  M <- pmax(tapply(y, x, mean), 1 / (2 * n))
  V <- pmax(tapply(y, x, var), 1 / n^3)
  
  LRR1 <- as.numeric(diff(log(M)))
  BC2 <- log(M) + V / (2 * n * M^2)
  LRR2 <- BC2[[2]] - BC2[[1]]
  BC3 <- log(M) + log(n + (n + 1) * sqrt(1 + ((2 * n + 1) * (n - 1) * V) / (n * (n + 1)^2 * M^2))) - log(2 * n + 1)
  LRR3 <- BC3[[2]] - BC3[[1]]
  V_LRR <- sum(V / (n * M^2))
  
  data.frame(LRR1 = LRR1, LRR2 = LRR2, LRR3 = LRR3, V_LRR = V_LRR)
}

sim_LRRs <- function(theta, m, n, mu_A, phi, iterations, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  x <- c(rep(0, m), rep(1, n))
  mu_vec <- mu_A * exp(x * theta)
  len <- m + n
  ar_vec <- ifelse(mu_vec[-1] - phi * mu_vec[-len] >= 0, phi, mu_vec[-1] / mu_vec[-len])
  lambda <- pmax(0, mu_vec[-1] - ar_vec * mu_vec[-len])
  
  LRRs_df <- 
    rerun(iterations, {
        y <- r_pois_AR1(mu_vec, ar_vec, len = len, lambda = lambda)
        calc_LRRs(x = x, y = y)
      }) %>%
    bind_rows()
  
  LRRs_df %>%
    summarize(
      across(c(LRR1, LRR2, LRR3), list(E = mean, V = var, cor = ~ cor(., V_LRR))),
      E_V = mean(V_LRR),
    ) %>%
    pivot_longer(-E_V, names_to = c("estimator", ".value"), names_pattern = "(LRR.)_(.+)") %>%
    mutate(
      bias = E - theta,
      V_RB = E_V / V,
    ) %>%
    select(estimator, bias, var = V, cor, V_RB)
}


# design factors
design_factors <- list(
  theta = seq(-2.5, 2.5, 0.25),
  m = c(5, 10, 15),
  n = c(5, 10, 15),
  mu_A = 10,
  phi = c(0, 0.2, 0.4)
)


params <- 
  cross_df(design_factors) %>%
  mutate(
    iterations = 10000,
    seed = round(runif(1) * 2^30) + 1:n()
  )

library(tictoc)
library(future)
library(furrr)
plan(multisession)

tic()

res <-
  params %>%
  mutate(res = future_pmap(., .f = sim_LRRs)) %>%
  unnest(cols = res)

toc()

saveRDS(res, "R/LRR-simulation-results.rds")