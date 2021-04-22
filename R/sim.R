#' Simulate (uncorrelated) multivariate normal data of size N x p
#'
#' @param N number of observations (will adjust if N not a multiple of clusters)
#' @param clusters number of clusters
#' @param p size of N x p matrix corresponding to number of predictors
#' @param ranef.sd standard deviation of simulated cluster random intercepts, assuming normal distribution
#' @param eff.sd standard deviation of distribution of effect sizes (difference between treat = 1 and treat = 0), assuming normal distribution with mean = 0
#' @param resid.sd optional standard deviation of residuals (inferred from icc if not supplied)
#' @param icc optional intra-class correlation (inferred if resid.sd supplied)
#' @param z current z of R bootstrap replicate, provided by boot_glmnet()
#' @param sparsity 1-proportion of non-zero coefficients (equivalently proportion of coefficients that are exactly zero)
#' @param treat vector of treatment group coefficients (defaults to sample(c(0,1), clusters, replace = TRUE))
#' @return list with $data data frame of simulated data, $eff.size data frame of effect sizes and standardized effects sizes for p variables, $icc, $resid.sd

sim <- function(N = 100,
                clusters = 20,
                p = 10,
                ranef.sd = .1,
                eff.sd = .1,
                resid.sd = NULL,
                icc = .5,
                sparsity = .5,
                treat = NULL){
  if(is.null(treat)){
    treat <- sample(c(0,1), clusters, replace = TRUE)
  }
  if(N %% clusters){
    N <- floor(N/clusters)*clusters
    warning('\nnumber of clusters not a multiple of N, using N =', N, '\n')

  }
  treat <- rep(treat, each = N/clusters)
  ranef <- rnorm(clusters, 0, ranef.sd)
  ranef.x <- rep(ranef, each = N/clusters)
  id <- rep(1:clusters, each = N/clusters)
  eff <- rnorm(p, 0, eff.sd)
  eff[sample(1:p, floor(p*sparsity))] <- 0
  if(!is.null(resid.sd)){
    icc <- ranef.sd^2/(ranef.sd^2+resid.sd^2)
    cat('\nICC =', icc)
  }
  if(is.null(resid.sd)){
    resid.sd <- sqrt((ranef.sd^2*(icc-1))/(-1*icc))
  }

  X <- treat%*%t(eff) + ranef.x + matrix(rnorm(N*p, 0, resid.sd), ncol = p)
  data = data.frame(id = id,
                    ranef.x = ranef.x,
                    treat = treat,
                    X)
  eff.size = data.frame(variable = paste0('X', 1:length(eff)),
                        eff.size = eff,
                        z = eff/resid.sd)

  return(list(data = data,
              eff.size = eff.size,
              icc = icc,
              resid.sd = resid.sd))

}
