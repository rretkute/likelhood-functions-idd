#' Negative Binomial likelihood
#'
#' @param sim 
#' @param obs 
#' @param k 
#'
#' @return log-likelihood value
#' @export
#'
#' 
lll_obs_NB <- function(sim, obs, k){
  L <- sum(dnbinom(obs, mu=sim, size=k,  log = TRUE))
  return(L)
}
