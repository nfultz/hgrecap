#' Bayesian Changepoint for Binomial Data
#'
#' This package provides a binomial variant of bayesian changepoint analysis.
#'
#' @author
#' Neal Fultz \email{njf@@zestfinance.com},
#'
#' @name hgrecap-package
#' @docType package
#' @useDynLib hgrecap
#' @import Rcpp
NULL

#' Hypergeometric Recapture Sampler
#'
#' Uses MCMC to sample the size.
#'
#' @param s a list of observation ids
#' @param prior a list specifying the lambda
#' @param mcmc.control see \code{\link{hgrecap.control}}
#'
#' @return a list containing samples and transition probabilities
#'
#' @export
hgrecap <- function(s, prior=list(lambda=0), mcmc.control=hgrecap.control()) {

  r <- length(s) + 1
  marked <- Reduce(union, s, c(), accumulate = TRUE)

  N_min <- length(marked[[r]])

  n_drawn   <-  sapply(s, length)
  K_marked   <-  sapply(marked[-r], length)
  d_unmarked <-  sapply(Map(setdiff, s, marked[-r]), length)

  rcppmm <- rcpp_hgrecap_metropolis(n_drawn, K_marked, d_unmarked, N_min, prior, mcmc.control)
}

#' MCMC control parameters
#'
#' @param mcmc.iterations Number of iterations to run
#' @param mcmc.thin an integer > 0, Number of iterations to thin between samples
#' @param mcmd.burnin an integer >= 0
#' @param verbose - Not used right now.
#'
#'
#' @export
hgrecap.control <- function(mcmc.iterations=500, mcmc.thin=5, mcmc.burnin=50, verbose=0){
  if(verbose > 0) message("Verbose > 0 not implemented yet")
  as.list(environment())
}


