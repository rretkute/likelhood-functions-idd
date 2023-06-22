#' Trachoma Stochastic model
#'
#' @param parameters 
#' @param initial_state 
#' @param tmax 
#'
#' @return dataframe with event timing, indicators and sizes of compartments
#' @export
#'
#' @examples
Trachoma_stochastic_model<-function(parameters, initial_state, t0, tmax){
  
  ## All events are coded using 3 indexes:
  # eta1: 1 if S->E, 0 otherwise

  # assign numbers in each state for first time step
  nS <- initial_state$S
  nI <- initial_state$I

  
  N<-nS+nI
 
  #initialise time
  time <- t0
  
  if(nI>0) {
    events<-data.frame(time=t0, eta1=1, S=nS, I=nI, nI=1)
  }
  
  while((time <= tmax & (nI>0))){ 
    r <- c(parameters$beta * nS * nI /N, # S-> I
           parameters$gamma * nI # I-> R 
    )
    rtotal <- sum(r)
    r<-cumsum(r)
    delta <- (-1/rtotal)*log(runif(1)) 
    time <- time + delta # Time to the next event
    rnd<-runif(1)*rtotal
    if (rnd <= r[1]) { 
      # S->E
      nS <- max(0, nS - 1) 
      nI<-nI+1
      events<-rbind(events, data.frame(time=time, eta1=1, S=nS, I=nI, nI=1))
    }   else {
      nS <- nS + 1
      nI<-max(0, nI-1)
      events<-rbind(events, data.frame(time=time, eta1=0, S=nS, I=nI, nI=0))
    }
  }  ##  END of while
  
  return(events)
}
