data{
  int<lower=0> N;
  int<lower=0> P;
  matrix[N, P] X;
  vector[N] Y;
  int<lower=0> K;
  int indicators[K, P];
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
  vector[K] lp;
  for(k in 1:K){
      lp[k] = normal_lpdf(Y |X * (beta .* to_vector(indicators[k, ])) + alpha, sigma_Y);
      lp[k] += bernoulli_lpmf(indicators[k, ] | p);
  }
}
model{
  target += log_sum_exp(lp);
  alpha ~ normal(0, 10);
  beta ~ normal(0, 5);
  sigma_Y ~ cauchy(0, 2.5);
  p ~ beta(a_p, b_p);
}
generated quantities{
  int model_index;
  model_index =  categorical_logit_rng(lp);
}