library(brms)
library(gamlss.dist)

# Beta-binomial model

data("cbpp", package = "lme4")
head(cbpp)


beta_binomial2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "vint1[n]"
)

stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
"

stanvars <- stanvar(scode = stan_funs, block = "functions")

bf(
  incidence | vint(size) ~ period + (1|herd), 
  data = cbpp, 
  family = beta_binomial2, 
  stanvars = stanvars
)

bbin_fit <- brm(
  incidence | vint(size) ~ period + (1|herd), data = cbpp, 
  family = beta_binomial2, stanvars = stanvars
)

summary(bbin_fit)


# Double-Poisson model

set.seed(20230913)
N <- 500
X <- rnorm(N)
mu <- exp(2 + 0.3 * X)
theta_inv <- 0.7
Y <- rDPO(N, mu = mu, sigma = theta_inv)
dat <- data.frame(X = X, Y = Y, maxval = 1000)


stancode_lpmf <- "
  real double_Poisson_lpmf(int X, real mu, real phi) {
    real ans;
    real A = inv(2) * log(phi) - phi * mu;
    if (X == 0)
      ans = A;
    else
      ans = A + X * (phi * (1 + log(mu)) - 1) - lgamma(X + 1) + (1 - phi) * X * log(X);
    return ans;
  }
"

double_Poisson <- custom_family(
  "double_Poisson", dpars = c("mu","phi"),
  links = c("log","log"),
  # lb = c(0, 0), ub = c(NA, NA),
  type = "int"
)


phi_prior <- prior(normal(0,10), class = "phi")

double_Poisson_stanvars <- stanvar(scode = stancode_lpmf, block = "functions")

get_prior(
  Y ~ X,
  family = double_Poisson,
  stanvars = double_Poisson_stanvars,
  prior = phi_prior,
  data = dat
)

dpo_spec <- brmsformula(
  Y ~ X,
  family = double_Poisson,
  stanvars = double_Poisson_stanvars,
  prior = phi_prior,
  data = dat
)

dpo_fit <- brm(
  Y ~ X,
  data = dat, 
  family = double_Poisson,
  prior = phi_prior,
  stanvars = double_Poisson_stanvars
)
