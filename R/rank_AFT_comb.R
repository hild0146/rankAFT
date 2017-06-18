#' Parameter Estimation for Rank Based Accelerated Failure Time Model
#'
#' This function can be used to estimate the parameters associated with the rank accelerated failure time model as described by Heller (2007)
#'
#' @param x.mat design matrix containing parameters of interest
#' @param surv.time survival time
#' @param surv.cens survival censoring variable where 1 represents a failure, and 0 right censoring
#'
#' @return a list containing
#'  \item{est}{a vector of parameter estimates}
#'  \item{var}{the individual variance estimate or covariance matrix}
#'  \item{opt.val}{the initial convergence tolerance (should be near zero)}
#'
#' @references Heller G. (2007) Smoothed Rank Regression With Censored Data
#'             Journal of the American Statistical Association, 102 (478), 552 - 559.

#' @examples
#'
#' n.pts <- 100
#'
#' betas <- c(-0.2, 0.5, 1)
#'
#' x.mat <- cbind(rnorm(n.pts, mean = -2, sd = 1),
#'                rnorm(n.pts, mean = 1, sd = 3),
#'                rbinom(n.pts, 1, 0.5))
#'
#' cens.time <- log(runif(n = n.pts, min = 0, max = 50))
#'
#' eps.i <- rnorm(n = n.pts, mean = 0, sd = 1)
#'
#' fail.time <- x.mat %*% betas + eps.i
#'
#' surv.time <- exp(pmin(cens.time, fail.time))
#'
#' surv.cens <- 1 * (fail.time < cens.time)
#'
#'rankAFT(x.mat, surv.time, surv.cens)
#'
#' @export rankAFT

rankAFT <- function(x.mat, surv.time, surv.cens, method, maxit = 2000) {

  # ----- initial parameter estimates and bandwidth -----

  bandwidth.optim <- function(beta.est, x.mat, surv.time, surv.cens){

    beta.res <- log(surv.time) - x.mat %*% beta.est

    comb <- fyg_rit_orig(surv_cens = surv.cens,
                         beta_res = beta.res,
                         x_mat = x.mat)

    comb <- comb * nrow(x.mat) ^ (-2)

    sum(abs(comb))
  }

  # use either nelder-mead or brent optimization depending on the number of
  # parameters

  if(ncol(x.mat) > 1) {

    opt.prelim <- optim(par = rep(0, ncol(x.mat)),
                        method = 'Nelder-Mead',
                        fn = bandwidth.optim,
                        x.mat = x.mat,
                        surv.time = surv.time,
                        surv.cens = surv.cens)

  } else {

    opt.prelim <- optim(par = 0,
                     lower = - 10,
                     upper = 10,
                     method = 'Brent',
                     fn = bandwidth.optim,
                     x.mat = x.mat,
                     surv.time = surv.time,
                     surv.cens = surv.cens)
  }

  beta.prelim <- opt.prelim$par

  res.prelim <- log(surv.time) - x.mat %*% beta.prelim

  h.prelim <- sd(res.prelim[surv.cens == 1]) * nrow(x.mat) ^ (-0.26)

  # --- final parameter estimates -----

  aft.optim <- function(beta.est, x.mat, surv.time, surv.cens, h){

    beta.res <- log(surv.time) - x.mat %*% beta.est

    comb <- fyg_rit_smooth(surv_cens = surv.cens,
                           beta_res = beta.res,
                           x_mat = x.mat,
                           h = h)

    comb <- comb * nrow(x.mat) ^ (-2)

    sum(abs(comb))
  }

  opt.final <- optim(par = beta.prelim,
                     fn = aft.optim,
                     method = method,
                     x.mat = x.mat,
                     surv.time = surv.time,
                     surv.cens = surv.cens,
                     h = h.prelim,
                     control = list(maxit = maxit))

  beta.final <- opt.final$par

  res.final <- log(surv.time) - x.mat %*% beta.final

  h.final <- sd(res.final[surv.cens == 1]) * nrow(x.mat) ^ (-0.26)

  # --- variance calculation

  A.n <- A_n_calc(surv_cens = surv.cens,
                  beta_res = res.final,
                  x_mat = x.mat,
                  h = h.final)

  A.n[upper.tri(A.n)] <- A.n[lower.tri(A.n)]

  A.n <- A.n * h.final ^ (-1) * nrow(x.mat) ^ (-3/2)


  V.n <- V_n_calc(surv_cens = surv.cens,
                  beta_res = res.final,
                  x_mat = x.mat,
                  h = h.final)

  V.n[upper.tri(V.n)] <- V.n[lower.tri(V.n)]

  V.n <- V.n * nrow(x.mat) ^ (-3)

  var.beta <- solve(A.n) %*% V.n %*% solve(A.n)

  out <- list(est = beta.final,
              var = var.beta,
              opt.prelim = opt.prelim,
              h.prelim = h.prelim,
              opt.final = opt.final,
              h.final = h.final)

  return(out)

}




