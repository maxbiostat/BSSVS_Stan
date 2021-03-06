source("generate_BSSVS_data.r")
true.p <- c(0, 1, 0, 1, 1, 1)
which(true.p == 1)
K <- length(true.p)
simu_data <- generate_BSSVS_data(n = 500, K = K, block_size = 2, rho = 0.0,
                                 p_inc = true.p, seed = 666)
pairs(simu_data$X)
ols <- lm(simu_data$Y ~ simu_data$X)
summary(ols)
c(1.5, simu_data$true.betas * true.p)

######################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

probZero <- 1/2
Q <- 1-((probZero)^(1/simu_data$P))
devtools::source_url("https://raw.githubusercontent.com/maxbiostat/opinion_pooling/master/code/beta_elicitator.r")
probPars <- elicit_beta(m0 = Q, c = 1) ## elicit a prior with mean Q and coefficient of variation 1.
curve(dbeta(x, probPars$a, probPars$b), ylab = "density", xlab = expression(p[i]), main = "Probability of inclusion")
BSSVS.data <- list(N = length(simu_data$Y),
                   Y = simu_data$Y,
                   X = simu_data$X,
                   P = simu_data$P,
                   tau = rep(1, simu_data$P),
                   c = rep(5, simu_data$P),
                   a_p = probPars$a,
                   b_p = probPars$b
)

marg_BSSVS.model <- stan_model(file = "stan/BSSVS_linear_regression.stan",
                       save_dso = TRUE)
system.time(
  MAP <- optimizing(marg_BSSVS.model, data = BSSVS.data,
                    hessian = TRUE)
)
get_confidences <- function(map_est){
  ## For MAP estimates with a Hessian matrix, return confidence intervals
  ## assuming estimates are asymptotically normal
  p <- ncol(map_est$hessian)
  mus <- map_est$par[1:p]
  sds <- sqrt(diag(solve(-map_est$hessian)))
  return(data.frame(lwr = qnorm(c(.025), mean = mus, sd = sds),
                    mean = mus,
                    upr = qnorm(c(.975), mean = mus, sd = sds),
                    sd = sds))
}
get_confidences(MAP)

marg_BSSVS.posterior <- sampling(marg_BSSVS.model, 
                                 data = BSSVS.data, 
                                 control = list(adapt_delta = 0.95, max_treedepth = 10))
print(marg_BSSVS.posterior, pars = c("p", "beta", "alpha", "sigma_Y", "model_index"))
check_hmc_diagnostics(marg_BSSVS.posterior)

beta.hat <- colMeans(extract(marg_BSSVS.posterior, 'beta')$beta)

p.hat <- colMeans(extract(marg_BSSVS.posterior, 'p')$p)

getBF <- function(p, q = .5){
  {p*(1-q)}/{q*(1-p)}
}
sapply(p.hat, getBF, q = Q) ## weird Bayes factors

models <- extract(marg_BSSVS.posterior, 'model_index')$model_index
which(BSSVS.data$indicators[median(models), ] == 1)
simu_data$included
simu_data$true.betas[simu_data$included]
beta.hat[simu_data$included]

barplot(table(models))