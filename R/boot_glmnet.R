usethis::use_pipe()

#' Fit a bootstrap replicate
#'
#' @param x Model matrix or sparse model matrix
#' @param y Dependent or outcome variable
#' @param cluster optional vector of length(y) labels denoting clusters
#' @param alpha mixing parameter ranging from 0 to 1, with 0 being ridge regression, 1 being the LASSO and values in-between being elastic net regression
#' @param lambda penalty term
#' @param t0 starting time, typically provided by boot_glmnet()
#' @param R number of bootstrap replicates
#' @param z current z of R bootstrap replicate, provided by boot_glmnet()
#' @param verbose print replicate number, elapsed time, and expected time to complete if TRUE
#' @param file if a character string is supplied, mc_/boot_glmnet will attempt to find a file with that name
#' @return RMSE, r-squared, and a matrix of penalized coefficients
bootfit <- function(x, y, cluster, alpha, lambda, t0, R, z, verbose, ...) {
  if (verbose) {
    cat('\rFitting bootstrap replicate',
        z,
        'of',
        R)
  }

  if(is.null(cluster)){
    k <- sample(1:length(y), length(y), replace = TRUE)
  } else {
    # very basic check of supplied cluster vector
    if(length(cluster) != length(y)){
      stop('\nLength of cluster vector != length of y vector\n')
    }
    # sample from clusters
    j <- sample(unique(cluster), length(unique(cluster)), replace = TRUE)
    # index rows of data then subset only those rows from sampled clusters
    deck <- 1:length(y)
    deck <- deck[cluster %in% j]
    # sample n = length(y) rows from clusters with replacement
    k <- sample(deck, length(y), replace = TRUE)
  }
  y <- y[k]
  x <- x[k, ]
  m <-
    glmnet::glmnet(
      x = x,
      y = y,
      alpha = alpha,
      lambda = lambda,
      family = "gaussian"
    )
  pred <- predict(m, x, lambda = lambda, type = "response")
  t1 <- Sys.time()
  elapsed <- as.numeric(t1) - as.numeric(t0)
  tms <- ifelse(elapsed / z * R < 120, 1, 60)
  tms.label <- ifelse(tms == 1, 'sec', 'min')
  if (verbose) {
    cat(
      '\b . Time elapsed',
      round(elapsed / tms, 1),
      tms.label,
      '\b. Expected time to complete',
      round((elapsed / z) / tms * R, 1),
      tms.label,
      '\b.             '
    )
  }
  return(list(
    RMSE = caret::RMSE(pred, y),
    R2 = caret::R2(pred, y)[1, 1],
    coefs = coef(m, lambda = lambda)
  ))
}

# this is the superceded bootfit function which does not allow user to
# account for any clustering in the data
bootfit0 <- function(x, y, alpha, lambda, t0, R, z, verbose, ...) {
  if (verbose) {
    cat('\rFitting bootstrap replicate',
        z,
        'of',
        R)
  }

  k <- sample(1:length(y), length(y), replace = TRUE)
  y <- y[k]
  x <- x[k, ]
  m <-
    glmnet::glmnet(
      x = x,
      y = y,
      alpha = alpha,
      lambda = lambda,
      family = "gaussian"
    )
  pred <- predict(m, x, lambda = lambda, type = "response")
  t1 <- Sys.time()
  elapsed <- as.numeric(t1 - t0)
  tms <- ifelse(elapsed / z * R < 120, 1, 60)
  tms.label <- ifelse(elapsed / z * R < 120, 'sec', 'min')
  if (verbose) {
    cat(
      '\b . Time elapsed',
      round(elapsed / tms, 1),
      tms.label,
      '\b. Expected time to complete',
      round((elapsed / z) / tms * R, 1),
      tms.label,
      '\b.             '
    )
  }
  return(list(
    RMSE = caret::RMSE(pred, y),
    R2 = caret::R2(pred, y)[1, 1],
    coefs = coef(m, lambda = lambda)
  ))
}


#' Fit bootstrap replicates using glmnet::glmnet
#'
#' @param x a model matrix or sparse representation of the model matrix
#' @param y a vector of the dependent or outcome observations
#' @param alpha mixing parameter which ranges from 0 to 1, where 0 is ridge regression and 1 is LASSO
#' @param lambda penalty term
#' @param verbose prints current progress and expected time to completion if TRUE
#' @param R number of bootstrap replicates or bootstrap replicates per block if computing in parallel
#' @return a boot_glmnet object
#' @export
boot_glmnet <-
  function(x,
           y,
           alpha,
           lambda,
           cluster = NULL,
           verbose = TRUE,
           R = 100,
           file = NULL,
           ...) {
    t0 <- Sys.time()
    if(is.character(file)){
      if(file.exists(file)){
        obj <- readRDS(file = file)
          if(class(obj) %in% c('boot_glmnet', 'boot_glmet_q')){
            return(obj)
          }
        warning('\nfile supplied was not boot_glmnet object\n')
      }
    }
    if (verbose) {
      cat('\n')
    }
    out <- lapply(1:R, function(z)
      bootfit(
        y = y,
        x = x,
        cluster = cluster,
        alpha = alpha,
        lambda = lambda,
        t0 = t0,
        R = R,
        z = z,
        verbose = verbose,
        ...
      ))
    out <- simplify2array(out)

    r2vec <- unlist(out['R2', ])
    rmsevec <- unlist(out['RMSE', ])
    coefs <- out['coefs', ]
    coefmat <- do.call(cbind, coefs)
    obj <- list(
      alpha = alpha,
      lambda = lambda,
      R = R,
      r2vec = r2vec,
      rmsevec = rmsevec,
      coefmat = coefmat
    )
    class(obj) <- 'boot_glmnet'
    if(is.character(file)){
      saveRDS(obj, file = file)
    }
    return(obj)
  }

#' Combine two boot_glmnet objects
#'
#' @param fit1 a boot_glmnet object
#' @param fit2 a boot_glmnet object, which must have the same alpha, lambda, and number of rows as fit1, but may have any number of columns > 1
#' @return a boot_glmnet object
#' @export
boot_combine <- function(fit1, fit2) {
  if (fit1$alpha != fit2$alpha) {
    warning('\nModel alphas not equal, cannot combine')
    return(NULL)
  }
  if (fit1$lambda != fit2$lambda) {
    warning('\nModel lambdas not equal, cannot combine')
    return(NULL)
  }
  if (nrow(fit1$coefmat) != nrow(fit2$coefmat)) {
    warning('\nModel matrixes not equal, cannot combine')
    return(NULL)
  }
  r2vec <- c(fit1$r2vec, fit2$r2vec)
  rmsevec <- c(fit1$rmsevec, fit2$rmsevec)
  coefmat <- cbind(fit1$coefmat, fit2$coefmat)
  obj <- list(
    alpha = fit1$alpha,
    lambda = fit1$lambda,
    R = fit1$R + fit2$R,
    r2vec = r2vec,
    rmsevec = rmsevec,
    coefmat = coefmat
  )
  class(obj) <- 'boot_glmnet'
  return(obj)

}

#' Compute bootstrap statistics
#'
#' @param obj a boot_glmnet object
#' @param ci specified confidence interval between 0 and 1
#' @return a boot_glmnet_q object
#' @export
quant.boot_glmnet <- function(obj, ci = .95) {
  if (ci > 1) {
    ci <- ci / 100
    warning('\nci > 1, assumed to be percentage (divided by 100)\n')
  }
  ci.lwr <- (1 - abs(ci)) / 2
  ci.upr <- 1 - ((1 - abs(ci)) / 2)
  R <- obj$R
  coefmat <- obj$coefmat
  r2vec <- obj$r2vec
  rmsevec <- obj$rmsevec

  coefs.lwr <-
    apply(coefmat, 1, function(z)
      quantile(z, ci.lwr, na.rm = TRUE))
  coefs.upr <-
    apply(coefmat, 1, function(z)
      quantile(z, ci.upr, na.rm = TRUE))
  coefs.m <- apply(coefmat, 1, function(z)
    mean(z, na.rm = TRUE))
  # p(selection|alpha, lambda)
  p.lasso <- apply(coefmat, 1, function(z)
    sum(abs(z) > 0) / R)
  # Approx CIs for p.lasso given R replicates
  p.lasso.lwr <- p.lasso - qnorm(ci.upr) * sqrt(p.lasso * (1 - p.lasso) /
                                                  R)
  p.lasso.upr <- p.lasso + qnorm(ci.upr) * sqrt(p.lasso * (1 - p.lasso) /
                                                  R)
  p.lasso.upr[p.lasso.upr > 1] <- 1
  p.lasso.lwr[p.lasso.lwr < 0] <- 0
  r2 <-
    c(mean(r2vec, na.rm = TRUE), quantile(r2vec, c(ci.lwr, ci.upr), na.rm = TRUE))
  rmse <-
    c(mean(rmsevec, na.rm = TRUE), quantile(rmsevec, c(ci.lwr, ci.upr), na.rm = TRUE))
  coefs <- data.frame(coefs.m,
                      coefs.lwr,
                      coefs.upr,
                      p.lasso,
                      p.lasso.lwr,
                      p.lasso.upr)
  obj <- list(
    r2 = r2,
    rmse = rmse,
    coefs = coefs,
    alpha = obj$alpha,
    lambda = obj$lambda,
    R = R,
    ci = ci
  )
  class(obj) <- 'boot_glmnet_q'
  return(obj)
}

#' Calculate proportion of bootstrap replicates where a predictor has been selected
#'
#' @param obj a boot_glmnet or boot_glmnet_q object
#' @return a data.frame with proportion of replicates where a predictor has been selected and an approximate confidence interval
#' @export
p.lasso <- function(obj, ...) {
  if (!class(obj) %in% c('boot_glmnet', 'boot_glmnet_q')) {
    warning('\nNot a boot_glmnet or boot_glmnet_q object')
  }
  x <- obj
  if (class(obj) == 'boot_glmnet') {
    x <- quant.boot_glmnet(obj, ...)
  }

  y <- x$coefs[, c(4:6)]
  colnames(y) <-
    c('p.lasso',
      paste0('lwr.', x$ci * 100),
      paste0('upr.', x$ci * 100))
  attr(y, 'ci') <- x$ci
  attr(y, "R") <- x$R
  attr(y, 'alpha') <- x$alpha
  attr(y, 'lambda') <- x$lambda
  return(y)
}

#' Return bootstrap coefficients and confidence intervals
#'
#' @param obj a boot_glmnet object
#' @return a data.frame of bootstrap coefficients
#' @export
coef.boot_glmnet <- function(obj, ...) {
  x <- quant.boot_glmnet(obj, ...)
  y <- x$coefs[, c(1:3)]
  colnames(y) <-
    c('Estimate',
      paste0('lwr.', x$ci * 100),
      paste0('upr.', x$ci * 100))
  attr(y, 'ci') <- x$ci
  attr(y, "R") <- x$R
  attr(y, 'alpha') <- x$alpha
  attr(y, 'lambda') <- x$lambda
  return(y)
}

#' Return bootstrap coefficients and confidence intervals
#'
#' @param obj a boot_glmnet_q object
#' @return a data.frame of bootstrap coefficients
#' @export
coef.boot_glmnet_q <- function(obj, ...) {
  x <- obj
  y <- x$coefs[, c(1:3)]
  colnames(y) <-
    c('Estimate',
      paste0('lwr.', x$ci * 100),
      paste0('upr.', x$ci * 100))
  attr(y, 'ci') <- x$ci
  attr(y, "R") <- x$R
  attr(y, 'alpha') <- x$alpha
  attr(y, 'lambda') <- x$lambda
  return(y)
}


sig_coefs <- function(obj,
                      p.lasso.thresh = 1,
                      ...) {
  if (class(obj) == 'boot_glmnet') {
    obj <- quant.boot_glmnet(obj, ...)
  }
  x <- obj$coefs
  flag <- sign(x[, 2]) == sign(x[, 3])
  flag[x[, 2] == 0] <- FALSE
  plt <- x$p.lasso > p.lasso.thresh

  names(x) <- c(
    'Estimate',
    paste0('lwr'),
    paste0('upr'),
    'p.lasso',
    paste0('p.lasso.lwr', obj$ci * 100),
    paste0('p.lasso.upr', obj$ci * 100)
  )
  return(data.frame(x, flag = flag, p.lasso.thresh = plt))
}

#' @export
summary.boot_glmnet <-
  function(obj,
           digits = 3,
           p.lasso.thresh = 1,
           silent = FALSE,
           ...) {
    obj <- quant.boot_glmnet(obj, ...)
    if(!silent){
      cat('\nalpha: ')
      print(obj$alpha)
      cat('\nlambda: ')
      print(obj$lambda)


      cat('\nRMSE: \n')
      print(obj$rmse, digits = digits)
      cat('\nR-squared:  \n')
      print(obj$r2, digits = digits)
      ci <- sig_coefs(obj,
                      p.lasso.thresh = p.lasso.thresh,
                      ...)
      cat(
        '\nCoefficients with significant coefficients at p <',
        1 - obj$ci,
        'or P(selection with LASSO) >',
        p.lasso.thresh,
        ':  \n'
      )
      ci$sig <- ifelse(ci$flag, '*', ' ')
      print(ci[ci$flag | ci$p.lasso.thresh, -c(7, 8)], digits = digits)
      cat('\n\n')
    }
    if(silent){
      ci <- sig_coefs(obj,
                      p.lasso.thresh = p.lasso.thresh,
                      ...)
      ci$sig <- ifelse(ci$flag, '*', ' ')
      ci[,-c('flag')]
    }
  }


#' @export
summary.boot_glmnet_q <-
  function(obj,
           p.lasso.thresh = 1,
           digits = 3,
           ...) {
    cat('\nalpha: ')
    print(obj$alpha)
    cat('\nlambda: ')
    print(obj$lambda)


    cat('\nRMSE: \n')
    print(obj$rmse, digits = digits)
    cat('\nR-squared:  \n')
    print(obj$r2, digits = digits)
    ci <- sig_coefs(obj,
                    p.lasso.thresh = p.lasso.thresh,
                    ...)
    cat(
      '\nCoefficients with significant coefficients at p <',
      1 - obj$ci,
      'or P(selection with LASSO) >',
      p.lasso.thresh,
      ':  \n'
    )
    ci$sig <- ifelse(ci$flag, '*', ' ')
    print(ci[ci$flag | ci$p.lasso.thresh, -c(7, 8)], digits = digits)
    cat('\n\n')
  }

#' @export
plot.boot_glmnet <- function(obj,
                             keep = NULL,
                             ...) {
  s <- sig_coefs(obj, ...) %>% tibble::rownames_to_column() %>%
    dplyr::rename(coef = rowname) %>%
    dplyr::filter(flag |
                    coef %in% keep | p.lasso.thresh, coef != '(Intercept)') %>%
    dplyr::mutate(coef = forcats::fct_reorder(coef, Estimate))

  g <-
    ggplot2::ggplot(s, aes(x = coef, y = Estimate, colour = flag)) +
    geom_errorbar(width = .2, aes(ymin = lwr, ymax = upr)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_colour_manual(values = c("black", "red")) +
    xlab("") +
    ylab("[Penalized] coefficient (in standardized units, SD)") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none")
  g
}

#' @export
plot.boot_glmnet_q <- function(obj,
                               keep = NULL,
                               ...) {
  s <- sig_coefs(obj, ...) %>% tibble::rownames_to_column() %>%
    dplyr::rename(coef = rowname) %>%
    dplyr::filter(flag |
                    coef %in% keep | p.lasso.thresh, coef != '(Intercept)') %>%
    dplyr::mutate(coef = forcats::fct_reorder(coef, Estimate))

  g <- ggplot(s, aes(x = coef, y = Estimate, colour = flag)) +
    geom_errorbar(width = .2, aes(ymin = lwr, ymax = upr)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_colour_manual(values = c("black", "red")) +
    xlab("") +
    ylab("[Penalized] coefficient (in standardized units, SD)") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none")
  g
}

#' Fit n x R blocks of bootstrap replicates in parallel using future.apply
#'
#' @param n number of blocks to compute in parallel
#' @param x Model matrix or sparse model matrix
#' @param y Dependent or outcome variable
#' @param alpha mixing parameter ranging from 0 to 1, with 0 being ridge regression, 1 being the LASSO and values in-between being elastic net regression
#' @param lambda penalty term
#' @param t0 starting time
#' @param R number of bootstrap replicates
#' @param z current z or R bootstrap replicate
#' @param verbose print replicate number, elapsed time, and expected time to complete if TRUE
#' @param file if a character string is supplied, mc_boot_glmnet will attempt to find a file with that name
#' @return RMSE, r-squared, and a matrix of penalized coefficients
#' @export
mc_boot_glmnet <- function(n = 4L,
                           x,
                           y,
                           alpha = 1,
                           lambda,
                           R = 100,
                           verbose = TRUE,
                           file = NULL,
                           ...) {
  if (n < 2) {
    warning('\n number of blocks n must be >=2, increasing n to 2\n')
    n = 2L
  }
  if(is.character(file)){
    if(file.exists(file)){
      obj <- readRDS(file = file)
      if(class(obj) %in% c('boot_glmnet', 'boot_glmet_q')){
        return(obj)
      }
      warning('\nfile supplied was not boot_glmnet object\n')
    }
  }
  if(!class(x) %in% 'dgCMatrix'){
    x <- Matrix::Matrix(x, sparse = TRUE)
  }
  future::plan("future::multisession")
  k <- future.apply::future_replicate(
    n,
    boot_glmnet(
      x = x,
      y = y,
      alpha = alpha,
      lambda = lambda,
      R = R,
      verbose = verbose,
      ...
    ),
    simplify = "list",
    future.seed = TRUE
  )
  future::plan("future::sequential")
  j <- k[, 1]
  for (i in 2:n) {
    j <- boot_combine(j, k[, i])
  }
  if(is.character(file)){
    saveRDS(j, file = file)
  }
  return(j)
}
