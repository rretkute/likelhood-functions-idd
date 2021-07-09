#' Poisson likelihood
#'
#' @param sim 
#' @param obs 
#'
#' @return log-likelihood value
#' @export
#'
#' 
lll_obs_pois <- function(sim, obs){
  L <- sum(dpois(obs, sim, log = TRUE))
  return(L)
}