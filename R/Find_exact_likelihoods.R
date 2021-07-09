#' find exact likelihoods
#'
#' @param events 
#' @param model 
#' @param beta.cand 
#' @param gamma.cand 
#' @param alpha.cand 
#'
#' @return
#' @export
#'
#' @examples
find_exact_likelihoods <- function(events, model, 
                                   beta.cand, gamma.cand = NULL, alpha.cand = NULL, N=NULL,
                                   other_pars = NULL){
  if(model=="flu"){
    ans<-data.frame(beta=c(), gamma=c(), ll=c())
    for(j in 1:length(beta.cand)){
      for(k in 1:length(gamma.cand)){
        ll=lll_exact_flu(events, data.frame(beta=beta.cand[j], gamma=gamma.cand[k], N=N))
        ans<-rbind(ans, data.frame(beta=beta.cand[j], gamma=gamma.cand[k], ll=ll))
      }
    }
  }
  if(model=="COVID19"){
  
  ans<-data.frame(beta=c(), alpha=c(), ll=c())
  for(j in 1:length(beta.cand)){
    for(k in 1:length(alpha.cand)){
      parameters <- data.frame(beta =  beta.cand[j], p = other_pars$p, gamma_1 = other_pars$gamma_1,
                      gamma_2 = other_pars$gamma_2, alpha = alpha.cand[k], N=N)
      ll<- lll_exact_COVID19(events, parameters)
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
