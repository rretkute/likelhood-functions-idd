#' find_NB_likelihoods
#'
#' @param obs 
#' @param model 
#' @param beta.cand 
#' @param gamma.cand 
#' @param alpha.cand 
#'
#' @return the log-likelihood values, range and maximum
#' @export
#'
#' @examples
find_NB_likelihoods <- function(obs, k, model, initial_state, tmax, dt,
                                     beta.cand, gamma.cand = NULL, alpha.cand = NULL,
                                     N=NULL, other_pars = NULL){
  
  if(model=="flu"){
    ans<-data.frame(beta=c(), gamma=c(), ll=c())
    for(j in 1:length(beta.cand)){
      for(k in 1:length(gamma.cand)){
        parameters <- c(beta =  beta.cand[j], gamma = gamma.cand[k], N=N)
        sim <- run_ODE_model(parameters, initial_state, tmax, dt,
                             model = Flu_ODE)
        sims<-pmax(0, diff(sim$R)) 
        ll<- lll_obs_NB(sims, diff(obs$R), k)
        ans<-rbind(ans, data.frame(beta=beta.cand[j], gamma=gamma.cand[k], ll=ll))
      }
    }
  }
  if(model=="COVID19"){
    
    ans<-data.frame(beta=c(), gamma=c(), ll=c())
    for(j in 1:length(beta.cand)){
      for(k in 1:length(alpha.cand)){
        parameters <- c(beta =  beta.cand[j], p = other_pars$p, gamma_1 = other_pars$gamma_1,
                        gamma_2 = other_pars$gamma_2, alpha = alpha.cand[k])
        sim <- run_ODE_model(parameters, initial_state, tmax, dt,
                             model = COVID_ODE)
        sims<-pmax(0,diff(sim$D)); 
        ll<- lll_obs_NB(sims, diff(obs$D), k)
        ans<-rbind(ans, data.frame(beta=beta.cand[j], alpha=alpha.cand[k], ll=ll))
      }
    }
    
  }
  # create indicator of finite ll
  ind <- which(ans$ll>-Inf)
  range_ll <- c(min(ans$ll[ind]), max(ans$ll[ind])) # range
  max_ll <- ans[ans$ll==max(ans$ll),] # maximum
  
  return(list(ans = ans, range = range_ll, max_ll = max_ll))
}
