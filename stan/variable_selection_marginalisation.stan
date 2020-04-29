data{
  int<lower=0> N;
  int<lower=0> P;
  matrix[N, P] X;
  vector[N] Y;
  real<lower=0> tauB0;
  real<lower=0> tau;
  int<lower=0> K;
  int indicators[K, P];
  real<lower=0> w;
  real<lower=0> a_g;
  real<lower=0> b_g;
} 
parameters{
  real<lower=0, upper=1> gamma[P];
  real alpha;
  vector[P] beta;
  real<lower=0> sigma_Y;
}
transformed parameters{
  vector[K] lp;
  for(k in 1:K){
      lp[k] = normal_lpdf(Y |X * (beta .* to_vector(indicators[k, ])) + alpha, sigma_Y);
      lp[k] += bernoulli_lpmf(indicators[k, ] | gamma);
  }
}
model{
  target += log_sum_exp(lp);
  alpha ~ normal(0, 1/sqrt(tauB0));
  beta ~ normal(0, 1/sqrt(tau));
  sigma_Y ~ cauchy(0, 2.5);
  gamma ~ beta(a_g, b_g);
}
generated quantities{
  int model_index;
  model_index =  categorical_logit_rng(lp);
}
