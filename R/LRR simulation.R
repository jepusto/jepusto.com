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
  BC3 <- log(M) + log(n + (n + 1) * sqrt(1 + ((2 * n + 1) * V) / (n * (n + 1) * M^2))) - log(2 * n + 1)
  LRR3 <- BC3[[2]] - BC3[[1]]
  V_LRR <- sum(V / (n * M^2))
  
  data.frame(estimator = paste0("LRR", 1:3), LRR = c(LRR1, LRR2, LRR3), V_LRR = V_LRR)
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
    group_by(estimator) %>%
    summarize(
      bias = mean(LRR) - theta,
      var = var(LRR),
      rmse = mean((LRR - theta)^2),
      cov = cov(LRR, V_LRR),
      V_RB = mean(V_LRR) / var(LRR)
    )
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
    iterations = 5000,
    seed = round(runif(1) * 2^30) + 1:n()
  )

nrow(params)

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

res <- readRDS("R/LRR-simulation-results.rds")

# Bias

res %>%
  filter(phi == 0) %>%
  ggplot(aes(theta, bias, color = estimator, shape = estimator)) + 
  facet_grid(m ~ n, labeller = "label_both") + 
  geom_point() + 
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) + 
  theme_minimal() +
  labs(x = "True LRR", y = "Bias")
  
res %>%
  filter(phi == 0.2) %>%
  ggplot(aes(theta, bias, color = estimator, shape = estimator)) + 
  facet_grid(m ~ n, labeller = "label_both") + 
  geom_point() + 
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) + 
  theme_minimal() +
  labs(x = "True LRR", y = "Bias")

res %>%
  filter(phi == 0.4) %>%
  ggplot(aes(theta, bias, color = estimator, shape = estimator)) + 
  facet_grid(m ~ n, labeller = "label_both") + 
  geom_point() + 
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) + 
  theme_minimal() +
  labs(x = "True LRR", y = "Bias")


# RMSE

res %>%
  filter(phi == 0) %>%
  ggplot(aes(theta, rmse, color = estimator, shape = estimator)) + 
  facet_grid(m ~ n, labeller = "label_both") + 
  geom_point() + 
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) + 
  theme_minimal() +
  labs(x = "True LRR", y = "RMSE")

res %>%
  filter(phi == 0.2) %>%
  ggplot(aes(theta, rmse, color = estimator, shape = estimator)) + 
  facet_grid(m ~ n, labeller = "label_both") + 
  geom_point() + 
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) + 
  theme_minimal() +
  labs(x = "True LRR", y = "RMSE")

res %>%
  filter(phi == 0.4) %>%
  ggplot(aes(theta, rmse, color = estimator, shape = estimator)) + 
  facet_grid(m ~ n, labeller = "label_both") + 
  geom_point() + 
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) + 
  theme_minimal() +
  labs(x = "True LRR", y = "RMSE")


# Variance

res %>%
  filter(phi == 0) %>%
  ggplot(aes(theta, var, color = estimator, shape = estimator)) + 
  facet_grid(m ~ n, labeller = "label_both") + 
  geom_point() + 
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) + 
  theme_minimal() +
  labs(x = "True LRR", y = "Variance")

res %>%
  filter(phi == 0.2) %>%
  ggplot(aes(theta, var, color = estimator, shape = estimator)) + 
  facet_grid(m ~ n, labeller = "label_both") + 
  geom_point() + 
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) + 
  theme_minimal() +
  labs(x = "True LRR", y = "Variance")

res %>%
  filter(phi == 0.4) %>%
  ggplot(aes(theta, var, color = estimator, shape = estimator)) + 
  facet_grid(m ~ n, labeller = "label_both") + 
  geom_point() + 
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) + 
  theme_minimal() +
  labs(x = "True LRR", y = "Variance")


# Covariance with V_LRR

res %>%
  filter(phi == 0) %>%
  ggplot(aes(theta, cov, color = estimator, shape = estimator)) + 
  facet_grid(m ~ n, labeller = "label_both") + 
  geom_point() + 
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) + 
  theme_minimal() +
  labs(x = "True LRR", y = "Covariance with V")

res %>%
  filter(phi == 0.2) %>%
  ggplot(aes(theta, cov, color = estimator, shape = estimator)) + 
  facet_grid(m ~ n, labeller = "label_both") + 
  geom_point() + 
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) + 
  theme_minimal() +
  labs(x = "True LRR", y = "Covariance with V")

res %>%
  filter(phi == 0.4) %>%
  ggplot(aes(theta, cov, color = estimator, shape = estimator)) + 
  facet_grid(m ~ n, labeller = "label_both") + 
  geom_point() + 
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) + 
  theme_minimal() +
  labs(x = "True LRR", y = "Covariance with V")

