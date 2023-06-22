#' Calculate likelihood for Poisson distributed observations
#'
#' @param sim 
#' @param obs 
#'
#' @return
#' @export
#'
#' @examples
lll_obs_pois<-function(sim, obs){
  L<-sum(dpois(obs[which(sim>=0)], sim[which(sim>=0)], log = TRUE))
  return(L)
}
