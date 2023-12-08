functions {

int gpo_quantile(real p, real mu, real phi) {
  int q = 0;
  real phi_sqrt = sqrt(phi);
  real mu_phi_sqrt = mu * phi_sqrt;
  real m;
  if (phi > 1) 
    m = mu / (1 - inv(phi_sqrt));
  else 
    m = positive_infinity();
  real lpmf = - mu_phi_sqrt;
  real cdf = exp(lpmf);
  real ln_inc;
  while (cdf < p && q < m) {
    q += 1;
    ln_inc = (q - 1) * log(mu_phi_sqrt + q * (1 - phi_sqrt)) - (q - 2) * log(mu_phi_sqrt + (q - 1) * (1 - phi_sqrt)) + (phi_sqrt - 1) - log(q);
    lpmf += ln_inc;    
    cdf += exp(lpmf);
  }
  return q;
}

int gpo_rng(real mu, real phi) {
  real p = uniform_rng(0,1);
  int x = gpo_quantile(p, mu, phi);
  return x;
}

}
