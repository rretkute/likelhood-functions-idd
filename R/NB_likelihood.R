#' Negative Binomial likelihood
#'
#' @param sim 
#' @param obs 
#'
#' @return
#' @export
#'
#' @examples
lll_obs_NB<-function(sim, obs, k){
#  L<-sum(dpois(obs, sim, log = TRUE))
  L<-sum(dnbinom(obs, mu=sim, size=k,  log = TRUE))
  return(L)
}
