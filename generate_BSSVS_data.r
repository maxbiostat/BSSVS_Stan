block_allocation <- function(K, block_size){
  ## works out how many blocks of what size are needed
  nblocks <- floor(K/block_size)
  remain <- K - (nblocks * block_size)
  if(remain == 0){
    sizes <- rep(block_size, nblocks)
  }else{
    sizes <- c(rep(block_size, nblocks), remain)
  }
  if(sum(sizes) != K) stop("Something is wrong, block sizes don't sum up to K")
  ans <- list(
    nblocks = length(sizes),
    sizes = sizes
  )
  return(ans)
}
##
make_corr_matrix <- function(rho, dim){
  ans <- matrix(rho, nrow = dim, ncol = dim)
  diag(ans) <- 1
  return(ans)
}
##
generate_predictors <- 
##
generate_BSSVS_data <- function(n, K, block_size = NULL,
                                rho = 0, intercept = 1.5,  p_inc, sd = 2, seed = NULL){
  ## Generates a design matrix with blocks of correlated covariates, each block of size 'block_size'
  ###  within each block the correlation is 'rho'.
  ## For an 'independent' design matrix, make 'block_size = 1'.
  ## 'p_inc' is the probability of inclusion of each predictor. For a deterministic inclusion (easier to debug) set them to 0 or 1.
  if(length(p_inc) != K) stop("p_inc should be of dimension K")
  if(!is.null(seed)) set.seed(seed)
  ## Assumes all blocks of covariates have the same internal correlation, easily relaxed
  if(is.null(block_size)) block_size <- 1
  block.specs <- block_allocation(K = K, block_size = block_size)
  ## generate design matrix
  J <- block.specs$nblocks
  Xs <- vector(J, mode = "list")
  corrMats <- vector(J, mode = "list")
  for(j in 1:J){
    corrMats[[j]] <- make_corr_matrix(rho = rho, dim = block.specs$sizes[j])
    Xs[[j]] <- MASS::mvrnorm(n, mu = rep(0, block.specs$sizes[j]), Sigma = corrMats[[j]])
  }
  ## notice the code below only works because I chose predictors to be "standardised" (sd = 1)
  full_X <- do.call(cbind, Xs)
  # if(K != block_size) full_X <- cbind(full_X, rnorm(n))
  colnames(full_X) <- paste("X", 1:K, sep = "")
  betas <- rnorm(K, mean = 0, sd = 2)
  inds <- rbinom(K, size = 1, prob = p_inc)
  Mu <- full_X%*% (betas * inds)
  Y <- rnorm(n, mean = Mu + intercept, sd = sd)
  return(
    list(
      Y = Y,
      X = full_X,
      P = K,
      included = which(as.logical(inds)),
      true.betas = betas
    )
  )
}