# Numerical validation

library(tidyverse)
library(metadat)
library(metafor)
library(clubSandwich)

dat.tannersmith2016$ESid <- 1:nrow(dat.tannersmith2016)
rho <- 0.7
Vmat <- with(dat.tannersmith2016, impute_covariance_matrix(vi = vi, cluster = studyid, r = rho, smooth_vi = TRUE))

TS_fit <- rma.mv(yi ~ 0 + sexmix, V = Vmat, 
                 random = ~ 1 | studyid / ESid,
                 data = dat.tannersmith2016,
                 method = "REML", sparse = TRUE)
TS_cat <- conf_int(TS_fit, vcov = "CR2")

Cmat <- constrain_equal("^sexmix", reg_ex = TRUE, 
                        coefs = coef(TS_fit))
Cmat_pivot <- vcov(TS_fit)[-1,-1] %*% Cmat

Wald_cat <- Wald_test(TS_fit, constraints = Cmat, vcov = "CR2")
Wald_test(TS_fit, constraints = Cmat_pivot, vcov = "CR2")

# calculate category-specific degrees of freedom
tau_sq <- TS_fit$sigma2[1]
omega_sq <- TS_fit$sigma2[2]

cat_calcs <- 
  dat.tannersmith2016 %>%
  group_by(sexmix, studyid) %>%
  summarise(
    ybar = mean(yi),
    kj = n(),
    sigma2j = mean(vi),
    .groups = "drop_last"
  ) %>%
  mutate(
    wj = kj/ (kj * tau_sq + (kj - 1) * rho * sigma2j + omega_sq + sigma2j)
  ) %>%
  summarise(
    ybar_c = weighted.mean(ybar, wj),
    W_c = sum(wj),
    E_VR = 1 / W_c,
    nu_VR = 1 / (sum(wj^2 / (W_c - wj)^2) - (2 / W_c) * sum(wj^3 / (W_c - wj)^2) + (1 / W_c^2) * sum(wj^2 / (W_c - wj))^2)
  )

df_Z <- 
  cat_calcs %>%
  summarize(
    q = n() - 1,
    num = q * (q + 1),
    den = 2 * sum((1 - 1 / (E_VR * sum(W_c)))^2 / nu_VR),
    df_Z = num / den
  ) %>%
  pull(df_Z)

mu_hat <- coef(TS_fit)
VR <- vcovCR(TS_fit, type = "CR2") %>% as.matrix()
q <- nrow(Cmat)
Q <- as.numeric(t(Cmat %*% mu_hat) %*% solve(Cmat %*% VR %*% t(Cmat)) %*% (Cmat %*% mu_hat))
delta <- (df_Z - q + 1) / (df_Z)
Fstat <- delta * Q / q
df_num <- q
df_den <- df_Z - q + 1
all.equal(df_num, Wald_cat$df_num)
all.equal(df_den, Wald_cat$df_denom)
all.equal(delta, Wald_cat$delta)
all.equal(Fstat, Wald_cat$Fstat)

