library(tidyverse)

f <- function(phi, mu, tau, omega, sigma, rho, k, s, alpha = .025) {
  crit <- qnorm(1 - alpha)
  Z <- (phi - crit * sigma) / sqrt(omega^2 + (1 - rho) * sigma^2)
  prob <- pnorm(Z)
  phi_sd <- sqrt(tau^2 + rho * sigma^2)
  if (phi_sd > 0) {
    dbinom(s, size = k, prob = prob) * dnorm(phi, mean = mu, sd = phi_sd)
  } else {
    dbinom(s, size = k, prob = prob)
  }
}

binom_norm_density <- function(s, k, mu, tau, omega, sigma, rho, alpha = .025) {
  phi_sd <- sqrt(tau^2 + rho * sigma^2)
  if (phi_sd > 0) {
    integrate(
      f, 
      lower = mu - 6 * phi_sd, 
      upper = mu + 6 * phi_sd,
      mu = mu, tau = tau, omega = omega, sigma = sigma, 
      rho = rho, k = k, s = s, alpha = alpha
    )$value
  } else {
    f(phi = mu, mu = mu, tau = tau, omega = omega,
      sigma = sigma, rho = rho, k = k, s = s, alpha = alpha)
  }
}

find_dist <- function(k, ESS, mu, tau, omega, rho, alpha = .025) {
  s <- 0:k
  sigma <- 2 / sqrt(ESS)
  probs <- sapply(
    s, 
    binom_norm_density, 
    k = k, mu = mu,
    tau = tau, omega = omega,
    sigma = sigma, rho = rho,
    alpha = alpha
  )
  data.frame(
    s = s, 
    p = probs
  )
}

find_dist_GH <- function(k, ESS, mu, tau, omega, rho, qp = 5, alpha = .025) {
  quad_points <- rmutil::gauss.hermite(qp)
  sigma <- 2 / sqrt(ESS)
  zeta_sd <- sqrt(tau^2 + rho * sigma^2)
  zeta <- zeta_sd * (quad_points[,"Points"]) + mu
  crit <- qnorm(1 - alpha)
  Z <- (zeta - crit * sigma) / sqrt(omega^2 + (1 - rho) * sigma^2)
  probs <- pnorm(Z)
  s <- 0:k
  h_qp <- sapply(probs, \(x) dbinom(s, size = k, prob = x))
  p <- as.vector(h_qp %*% (quad_points[,"Weights"]))
  data.frame(s = s, p = p)
}

k <- 10
ESS <- 20
mu <- 0.5
tau <- 0.1
omega <- 0.05
rho <- 0.5

res <- 
  inner_join(
    find_dist(k, ESS, mu, tau, omega, rho),
    find_dist_GH(k, ESS, mu, tau, omega, rho, qp = 21),
    by = "s"
  )
res
ggplot(res) + 
  geom_point(aes(x = s, y = p.x), color = "blue") + 
  geom_point(aes(x = s, y = p.y), color = "green")

