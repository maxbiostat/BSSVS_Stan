source("generate_BSSVS_data.r")
simu_data <- generate_BSSVS_data(n = 500, K = 5, block_size = 3, rho = 0.6,
                                 p_inc = c(0, 1, 0, 1, 1), seed = 666)
lm(simu_data$Y ~ simu_data$X)
c(1.5, simu_data$true.betas)

######################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

probZero <- 1/2
Q <- 1-((probZero)^(1/simu_data$P))
devtools::source_url("https://raw.githubusercontent.com/maxbiostat/opinion_pooling/master/code/beta_elicitator.r")
# gammaPars <- elicit_beta(m0 = Q, c = 1) ## elicit a prior with mean Q and coefficient of variation 1.
gammaPars <- list(a = 1, b = 1)
curve(dbeta(x, gammaPars$a, gammaPars$b), ylab = "density", xlab = expression(gamma[i]), main = "Probability of inclusion")
Ind <- expand.grid(replicate(simu_data$P, 0:1, simplify = FALSE))
BSSVS.data <- list(N = length(simu_data$Y),
                   P = simu_data$P,
                   tau = 4,
                   tauB0 = 10,
                   Y = simu_data$Y,
                   X = simu_data$X,
                   indicators = Ind,
                   K = nrow(Ind),
                   w = Q,
                   a_g = gammaPars$a,
                   b_g = gammaPars$b
)

marg_BSSVS.model <- stan_model(file = "stan/variable_selection_marginalisation.stan",
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
print(marg_BSSVS.posterior, pars = c("gamma", "beta", "alpha", "sigma_Y", "model_index"))
check_hmc_diagnostics(marg_BSSVS.posterior)

beta.hat <- colMeans(extract(marg_BSSVS.posterior, 'beta')$beta)

gammas.hat <- colMeans(extract(marg_BSSVS.posterior, 'gamma')$gamma)

getBF <- function(p, q = .5){
  {p*(1-q)}/{q*(1-p)}
}
sapply(gammas.hat, getBF) ## weird Bayes factors

models <- extract(marg_BSSVS.posterior, 'model_index')$model_index
which(BSSVS.data$indicators[median(models), ] == 1)
simu_data$included
simu_data$true.betas[simu_data$included]
beta.hat[simu_data$included]