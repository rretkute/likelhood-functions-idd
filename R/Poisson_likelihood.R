#' Poisson likelihood
#'
#' @param sim 
#' @param obs 
#'
#' @return
#' @export
#'
#' @examples
lll_obs_pois<-function(sim, obs){
  L<-sum(dpois(obs, sim, log = TRUE))
  return(L)
}