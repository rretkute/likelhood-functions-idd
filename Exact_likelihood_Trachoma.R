#' Process based likelihood for trachoma
#'
#' @param events observed infection and removal times, indicators, number of susceptible and infectious at time of obervation 
#' @param parameters
#'
#' @return
#' @export
#'
#' @import statip
#' @examples
lll_exact_Trachoma<-function(events, parameters){
  events<-events[order(events$time),]
  L<-0
  k<-nrow(events) #
  for(i in 2:k){
    SS<-events$S[i-1]
    II<-events$I[i-1]
    ll<-log( (parameters$beta*SS*II/parameters$N)^events$eta1[i] * 
        (parameters$gamma*II)^(1-events$eta1[i]))-                
       ((parameters$beta)*SS*II/parameters$N + parameters$gamma*II)*
      (events$time[i]-events$time[i-1])
    L<-L+ll
  }
  return(L)
}
