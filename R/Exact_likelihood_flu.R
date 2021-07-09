#' Exact likelihood for flu
#'
#' @param events observed infection and removal times, indicators, number of susceptible and infectious at time of observation 
#' @param parameters
#'
#' @return log-likelihood value
#' @export
#'
#'
lll_exact_flu<-function(events, parameters){
  events<-events[order(events$time),]
  L<-0
  k<-nrow(events)
  for(i in 2:k){
    I<-events$I[i-1]
    S<-events$S[i-1]
    L<-L+log(((parameters$beta)*S*I/parameters$N)^events$eta[i] * (parameters$gamma*I)^(1-events$eta[i]))-
      ((parameters$beta)*S*I/parameters$N + parameters$gamma*I)*(events$time[i]-events$time[i-1])
  }
  return(L)
}