#' Exact likelihood for flu
#'
#' @param events observed infection and removal times, indicators, number of susceptible and infectious at time of observation 
#' @param parameters
#'
#' @return log-likelihood value
#' @export
#'
#' @import statip
#' 
lll_exact_COVID19<-function(events, parameters){
  events<-events[order(events$time),]
  L<-0
  k<-nrow(events) # 
  for(i in 2:k){
    S<-events$S[i-1]
    I1<-events$I1[i-1]
    I2<-events$I2[i-1]
    L<-L+log( ((parameters$beta)*S*(I1+I2)/parameters$N)^events$eta1[i] * 
                (dbern(events$eta2[i], parameters$p, log=FALSE))^events$eta1[i]*
                (parameters$gamma_1*I1)^((1-events$eta2[i])*(1-events$eta1[i]))*
                (parameters$gamma_2*I2)^((1-events$eta1[i])*events$eta2[i]*(1-events$eta3[i]))*
                (parameters$alpha*I2)^((1-events$eta1[i])*events$eta2[i]*events$eta3[i]))-
      ((parameters$beta)*S*(I1+I2)/parameters$N+  parameters$gamma_1*I1 + 
         parameters$gamma_2*I2 +parameters$alpha*I2)*(events$time[i]-events$time[i-1])
  }
  return(L)
}