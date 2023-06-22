#' Exact likelihood for leprosy
#'
#' @param events observed infection and removal times, indicators, number of susceptible and infectious at time of obervation 
#' @param parameters
#'
#' @return
#' @export
#'
#' @import statip
#' @examples
lll_exact_Leprosy<-function(events, parameters){
  events<-events[order(events$time),]
  L<-0
  k<-nrow(events) #
  for(i in 2:k){
    SS<-events$S[i-1]
    EE<-events$E[i-1]
    II<-events$I[i-1]
    DD<-events$D[i-1]
    RR<-events$R[i-1]
    ll<-log( (parameters$beta*SS*(II+DD)/parameters$N)^(events$eta1[i]*(1-events$eta4[i])) * 
                (dbern(events$eta2[i], parameters$pE, log=FALSE))^((1-events$eta1[i])*(1-events$eta3[i])*(1-events$eta4[i]))*
                (parameters$sigma*EE)^((1-events$eta1[i])*(1-events$eta3[i])*(1-events$eta4[i]))*
                ((parameters$rho+parameters$gamma)*DD)^((1-events$eta1[i])*events$eta2[i]*events$eta3[i]*(1-events$eta4[i]))*
                (parameters$gamma*II)^((1-events$eta1[i])*(1-events$eta2[i])*events$eta3[i]*(1-events$eta4[i]))*
                (parameters$mu*EE)^(events$eta1[i]*(1-events$eta2[i])*(1-events$eta3[i])*events$eta4[i])*
                (parameters$mu*II)^((1-events$eta1[i])*(1-events$eta2[i])*(1-events$eta3[i])*events$eta4[i])*
                (parameters$mu*DD)^((1-events$eta1[i])*events$eta2[i]*(1-events$eta3[i])*events$eta4[i])*
                (parameters$mu*RR)^((1-events$eta1[i])*(1-events$eta2[i])*events$eta3[i]*events$eta4[i])
              )-                
       ((parameters$beta)*SS*(II+DD)/parameters$N+  
          parameters$sigma*EE + 
          (parameters$rho+parameters$gamma)*DD+
          parameters$gamma*II+
          parameters$mu*(EE+II+DD+RR))*
      (events$time[i]-events$time[i-1])
    L<-L+ll
  }
  return(L)
}


