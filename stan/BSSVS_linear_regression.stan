functions{
  real coefficient_lpdf(vector b, real[] p, real[] t, real[] c){
    int K = size(p);
    real ldens[K];
    for (i in 1:K){
      ldens[i] = log_mix(p[i], normal_lpdf(b[i] | 0, t[i]), normal_lpdf(b[i] | 0, c[i]*t[i]));
    }
    return(log_sum_exp(ldens));
  }
}
data{
  int<lower=0> N;
  int<lower=0> P;
  matrix[N, P] X;
  vector[N] Y;
  real<lower=0> tau[P];
  real<lower=0> c[P];
  real<lower=0> a_p;
  real<lower=0> b_p;
} 
parameters{
  real<lower=0, upper=1> p[P];
  real alpha;
  vector[P] beta;
  real<lower=0> sigma_Y;
}
transformed parameters{
}
model{
  alpha ~ normal(0, 10);
  sigma_Y ~ cauchy(0, 2.5);
  p ~ beta(a_p, b_p);
  beta ~ coefficient(p, tau, c);
  Y ~ normal(X * beta + alpha, sigma_Y);
}
generated quantities{
}
