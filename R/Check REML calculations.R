library(metafor)
library(clubSandwich)
tausq <- 0.1^2
omegasq <- 0.05^2
phi <- 0.4
J <- 6
k_j <- 1 + rpois(J, 3)
sigmasq_j <- 4 / pmax(rgamma(J, 2, scale = 30), 20)
sID <- rep(LETTERS[1:J], k_j)
sigmasq_ij <- rep(sigmasq_j, k_j)

make_omega_mat <- function(sID, sigmasq_ij, tausq, omegasq, rho = 0) {
  tausq_mat <- tapply(sigmasq_ij, sID, 
                      function(x) matrix(tausq, nrow = length(x), ncol = length(x)), 
                      simplify = FALSE)
  omegasq_mat <- tapply(sigmasq_ij, sID, 
                        function(x) diag(omegasq, nrow = length(x)), 
                        simplify = FALSE)
  Omega_rho <- impute_covariance_matrix(vi = sigmasq_ij, cluster = sID, r = rho)
  mapply(function(a,b,c) a + b + c, 
         a = tausq_mat, b = omegasq_mat, c = Omega_rho)
  
}

Phi_mat <- make_omega_mat(sID = sID, sigmasq_ij = sigmasq_ij, 
                          tausq = tausq, omegasq = omegasq, rho = phi)

tausq_t <- 0.2^2
omegasq_t <- 0.02^2
rho <- 0.8
wj_tilde <- k_j / (k_j * tausq_t + k_j * rho * sigmasq_j + omegasq_t + (1 - rho) * sigmasq_j)
W_tilde <- sum(wj_tilde)

Omega_mat <- make_omega_mat(sID = sID, sigmasq_ij = sigmasq_ij, 
                            tausq = tausq_t, omegasq = omegasq_t, rho = rho)

Omega_inv <- lapply(Omega_mat, function(x) chol2inv(chol(x)))
sapply(Omega_inv, sum)
all.equal(sapply(Omega_inv, sum), wj_tilde, check.attributes = FALSE)

tr2 <- wj_tilde * (1 + (k_j - 1) * (tausq_t + rho * sigmasq_j) / (omegasq_t + (1 - rho) * sigmasq_j))
tr2
sapply(Omega_inv, function(x) sum(diag(x)))
all.equal(sapply(Omega_inv, function(x) sum(diag(x))), tr2, check.attributes = FALSE)

Jmat <- tapply(sigmasq_ij, sID, 
               function(x) matrix(1L, nrow = length(x), ncol = length(x)), 
               simplify = FALSE)

mapply(function(a,b) sum(a %*% b %*% a), a = Omega_inv, b = Jmat)
wj_tilde^2
all.equal(mapply(function(a,b) sum(a %*% b %*% a), a = Omega_inv, b = Jmat), wj_tilde^2, check.attributes = FALSE)

sapply(Omega_inv, function(a) sum(a %*% a))
wj_tilde^2 / k_j
all.equal(sapply(Omega_inv, function(a) sum(a %*% a)), wj_tilde^2 / k_j, check.attributes = FALSE)

wj_star <- k_j / (k_j * tausq + k_j * phi * sigmasq_j + omegasq + (1 - phi) * sigmasq_j)
Omega_inv_big <- bldiag(Omega_inv)
Omega_inv_1 <- rowSums(Omega_inv_big)
Qmat <- Omega_inv_big - tcrossprod(Omega_inv_1) / W_tilde
Phi_big <- bldiag(Phi_mat)
QJQPhi <- Qmat %*% bldiag(Jmat) %*% Qmat %*% Phi_big
sum(diag(QJQPhi))

Qform1 <- sum(wj_tilde^2 / wj_star) - 
  2 * sum(wj_tilde^3 / wj_star) / W_tilde + 
  sum(wj_tilde^2) * sum(wj_tilde^2 / wj_star) / W_tilde^2
Qform1
all.equal(sum(diag(QJQPhi)), Qform1)

Oi_Phi_j <- mapply(function(a,b) a %*% b, a = Omega_inv, b = Phi_mat)
Oi_Phi_j_ <- 
  mapply(function(a, b, c, d) a / b - c * d, 
         a = Phi_mat, 
         b = omegasq_t + (1 - rho) * sigmasq_j,
         c = Jmat, 
         d = wj_tilde * (tausq_t + rho * sigmasq_j) / wj_star / (omegasq_t + (1 - rho) * sigmasq_j))
all.equal(Oi_Phi_j, Oi_Phi_j_)

Oi_Phi_Oi_j <- mapply(function(a,b) a %*% b %*% a, a = Omega_inv, b = Phi_mat)
Oi_Phi_Oi_j_ <- 
  mapply(function(a, b, c, d, e) a / b^2 - c * d / b^2 - c * e / b, 
         a = Phi_mat, 
         b = omegasq_t + (1 - rho) * sigmasq_j,
         c = Jmat, 
         d = wj_tilde * (tausq_t + rho * sigmasq_j) / wj_star,
         e = wj_tilde^2 * (tausq_t + rho * sigmasq_j) / (wj_star * k_j))
all.equal(Oi_Phi_Oi_j, Oi_Phi_Oi_j_)

tr_Phi_j <- sapply(Phi_mat, function(a) sum(diag(a)))
tr_Phi_j_ <- k_j / wj_star + (k_j - 1) * (omegasq + (1 - phi) * sigmasq_j)
all.equal(tr_Phi_j, tr_Phi_j_, check.attributes = FALSE)

tr_Oi_Phi_Oi <- mapply(function(a,b) sum(diag(a %*% b %*% a)),
                       a = Omega_inv, b = Phi_mat)

Qform2a <- (k_j / wj_star + (k_j - 1) * (omegasq + (1 - phi) * sigmasq_j) 
            - (tausq_t + rho * sigmasq_j) * k_j * wj_tilde / wj_star
            - (tausq_t + rho * sigmasq_j) * (omegasq_t + (1 - rho) * sigmasq_j) * wj_tilde^2 / wj_star ) / 
                 (omegasq_t + (1 - rho) * sigmasq_j)^2
all.equal(tr_Oi_Phi_Oi, Qform2a, check.attributes = FALSE)

Oi_Phi_Oi_Oi <- mapply(function(a,b) 2 * sum(a %*% b %*% a %*% a) / W_tilde,
                       a = Omega_inv, b = Phi_mat)
Qform2b <- 2 * wj_tilde^3 / (k_j * wj_star) / W_tilde
all.equal(Oi_Phi_Oi_Oi, Qform2b, check.attributes = FALSE)

sum_Oi_Phi_Oi <- mapply(function(a,b) sum(a %*% b %*% a),
                            a = Omega_inv, b = Phi_mat)
sum_Oi_Oi <- sapply(Omega_inv, function(a) sum(a %*% a))

Qform2c <- sum(wj_tilde^2 / k_j) * sum(wj_tilde^2 / wj_star) / W_tilde^2
all.equal(sum(sum_Oi_Oi) * sum(sum_Oi_Phi_Oi) / W_tilde^2, Qform2c)

PhiQ <- Qmat %*% Phi_big %*% Qmat
sum(diag(QPhiQ))
