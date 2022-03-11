#' Estimating marginalized zero-inflated Poisson model
#'
#' Function to estimate a marginalized zero-inflated Poisson model
#'
#' @references Long, D. L., Preisser, J. S., Herring, A. H., & Golin, C. E. (2014). A marginalized zero‐inflated Poisson regression model with overall exposure effects. Statistics in Medicine, 33(29), 5151-5165.
#' @param formula an object of class "\code{formula}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. A typical formula has the form \code{response ~ terms} where \code{response} is the count response vector and \code{terms} is a series of terms that predict response. For example, \code{formula = y ~ x1 + x2 + x3}. Do not write intercept in the \code{formula}; intercept will be automatically added in model fitting.
#' @param data a data frame containing variables in the model.
#' @importFrom bbmle parnames mle2
#' @details Function returns an object of class "\code{mle2}" from \pkg{bbmle} R package. Apply \code{summary} function to the resulting object from the function to obtain more estimation information.
#' @return
#' Suffix \code{_zero} corresponds to the parameters associated with the structrual zero rate part of a model. \cr \cr
#' Suffix \code{_mean} corresponds to the parameters associated with the overall mean, which evaluate the effects of covariates on the overall mean. \cr
#' @examples
#' head(dat.pfi)
#'
#' #Fit a marginalized zero-inflated Poisson model
#' res = mzip(formula = y ~ m0 + int_PF + year_new + race_new, data = dat.pfi)
#'
#' #Obtain estimation results
#' bbmle::summary(res)
#' @export

mzip <- function(formula, data)
{

  covs = data[, labels(stats::terms(formula))]
  y = data[, all.vars(formula)[1]]
  intercept = 1

  dat = data.frame(intercept, covs, y)

  names = c("intercept", labels(stats::terms(formula)))


  numPar <- ncol(dat)-1
  covnames <- names #colnames(dat)[1:numPar]

  minus_ll_function = function(par)
  {
    dat2 <- dat
    logit_psi <- rep(0,nrow(dat2))
    log_mui <- rep(0,nrow(dat2))
    for(k in 1:numPar)
    {
      logit_psi = logit_psi + par[k]*dat2[,k]
      log_mui = log_mui + par[k+numPar]*dat2[,k]
    }
    psi = exp(logit_psi)/(1+exp(logit_psi))
    mui = exp(log_mui)
    lambdai = mui/(1-psi)
    dat2$psi <- psi
    dat2$mui <- mui
    dat2$lambdai <- lambdai
    dat2$iszero <- 0
    dat2$iszero[which(dat2$y == 0)] <- 1
    dat2$ll1 <- log(dat2$psi+exp(-dat2$lambdai)*(1-dat2$psi))*dat2$iszero
    dat2$ll2 <- log((1-dat2$psi)*exp(-dat2$lambdai)*dat2$lambdai^(dat2$y)/factorial(dat2$y))*(1-dat2$iszero)
    dat2$ll <- dat2$ll1+dat2$ll2
    output = -sum(dat2$ll)
    return(output)
  }

  bbmle::parnames(minus_ll_function) = c(paste0(covnames, "_zero"), paste0(covnames, "_mean"))


  start_list = rep(0, 2*numPar)

  names(start_list) = c(paste0(covnames, "_zero"), paste0(covnames, "_mean"))


  est_pars <- bbmle::mle2(minuslogl = minus_ll_function, start  = start_list, method = "BFGS")

  return(est_pars)
}
