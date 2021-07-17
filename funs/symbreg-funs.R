#' functions for symbolic regression
#' author: mwelz
#' last changed: Jul 17, 2021


# void function for input checks
input.check.X <- function(data.lower, data.upper){

  if(any(dim(data.lower) != dim(data.upper))) stop("data need to be of same dimension")
  if(any(data.lower > data.upper)) stop("lower bound cannot exceed upper bound")
  
} # FUN


# void function for input checks
input.check.Y <- function(data.lower, data.upper){
  
  if(!is.vector(data.lower) | !is.vector(data.upper)) stop("data need to be a vector")
  if(any(data.lower > data.upper)) stop("lower bound cannot exceed upper bound")
  
} # FUN


# helper function to get X'X with symbolic data
get.XX <- function(X.lower, X.upper){
  
  # sum() isn't vectorized, so we need to vectorize its mapping function, see here: https://stackoverflow.com/questions/18110397/why-doesnt-outer-work-the-way-i-think-it-should-in-r
  outer(1:ncol(X.lower), 1:ncol(X.upper),
        FUN = Vectorize(function(i,j){
          sum((X.upper[,i] + X.lower[,i]) * (X.upper[,j] + X.lower[,j]))
        })) * 0.25 / nrow(X.lower)
} # FUN


# helper function to get X'Y with symbolic data
get.XY <- function(X.lower, X.upper, Y.lower, Y.upper){
  
  Ysum <- Y.upper + Y.lower
  
  sapply(1:ncol(X.lower), FUN = function(i){
    sum((X.upper[,i] + X.lower[,i]) * Ysum)
  }) * 0.25 / nrow(X.lower)
 
} # FUN


#' returns the symbolic OLS estimate (with an intercept)
#' 
#' @param X.lower design matrix of lower bounds for the explanatory variables
#' @param X.upper design matrix of upper bounds for the explanatory variables
#' @param Y.lower vector of lower bounds for the dependent variable
#' @param Y.upper vector of upper bounds for the dependent variable
#' 
#' @note Unlike in Billard and Diday (2000), X'X and X'Y are corrected by dividing by the number of observations in order to meet the symbolic cross-product expressions derived there. This obviously doesn't change the symbolic OLS estimate.
#' 
#' @export
symbolic.regression <- function(X.lower, X.upper, Y.lower, Y.upper){
  
  # input checks
  input.check.X(data.lower = X.lower, data.upper = X.upper)
  input.check.Y(data.lower = Y.lower, data.upper = Y.upper)
  
  # get X'X matrix (with intercepts)
  XX <- get.XX(X.lower = cbind(1, X.lower), 
               X.upper = cbind(1, X.upper))
  
  # get X'Y matrix
  XY <- get.XY(X.lower = cbind(1, X.lower), 
               X.upper = cbind(1, X.upper),
               Y.lower = Y.lower,
               Y.upper = Y.upper)
  
  # return the symbolic OLS estimate 
  return(as.numeric(solve(XX) %*% XY))
  
} # FUN



#' returns the symbolic OLS estimate (with an intercept) and its standard errors by means of bootstrapping
#' 
#' @param X.lower design matrix of lower bounds for the explanatory variables
#' @param X.upper design matrix of upper bounds for the explanatory variables
#' @param Y.lower vector of lower bounds for the dependent variable
#' @param Y.upper vector of upper bounds for the dependent variable
#' @param num.bootstrap.samples number of bootstrap samples. Default is 2000.
#' 
#' @note Unlike in Billard and Diday (2000), X'X and X'Y are corrected by dividing by the number of observations in order to meet the symbolic cross-product expressions derived there. This obviously doesn't change the symbolic OLS estimate.
#' 
#' @export
symbolic.regression_bootstrap <- function(X.lower, X.upper, Y.lower, Y.upper,
                                          num.bootstrap.samples = 2000){
  
  # initialize
  n            <- length(Y.lower)
  estimates_bs <- matrix(NA_real_, num.bootstrap.samples, ncol(X.lower) + 1)
  
  for(b in 1:num.bootstrap.samples){
    
    # reshuffle the indices to generate bootstrap draws
    idx.reshuffled   <- sample(1:n, size = n, replace = TRUE)
    
    # obtain bootstrap coefficient estimates
    estimates_bs[b,] <- symbolic.regression(X.lower = X.lower[idx.reshuffled,],
                                            X.upper = X.upper[idx.reshuffled,],
                                            Y.lower = Y.lower[idx.reshuffled],
                                            Y.upper = Y.upper[idx.reshuffled])
    
  } # FOR
  
  # return bootstrap estimates of location and scale of the coefficients
  return(list(bootstrap.coefficients = apply(estimates_bs, MARGIN = 2, FUN = mean),
              boostrap.standarderrors = sqrt(apply(estimates_bs, MARGIN = 2, FUN = var))))
  
} # FUN
