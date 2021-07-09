#' Exact likelihood for flu
#'
#' @param events observed infection and removal times, indicators, number of susceptible and infectious at time of obervation 
#' @param parameters
#'
#' @return
#' @export
#'
#' @examples
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