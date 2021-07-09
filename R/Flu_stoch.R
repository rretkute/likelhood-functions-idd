#' Flu Stochastic model
#'
#' @param parameters 
#' @param initial_state 
#' @param tmax 
#'
#' @return
#' @export
#'
#' @examples
Flu_stochastic_model<-function(parameters, initial_state, tmax){
  # eta: 1 if S->I, 0 if I->R
  
  # assign numbers in each state for first time step
  nS <- initial_state$S
  nI <- initial_state$I
  nR<- initial_state$R

  
  #initialise time
  time <- 0
  
  events<-data.frame(time=0, eta=1, S=nS, I=nI)
  
  
  while(time <= tmax & nI>0){ 
    r <- c(parameters$beta * nS * nI /parameters$N, parameters$gamma * nI)
    rtotal <- sum(r)
    delta <- (-1/rtotal)*log(runif(1)) 
    time <- time + delta # Time to the next event
    if (runif(1)*rtotal <= r[1]) { # Transmission event
      nS <- max(0, nS - 1)
      nI<-nI+1
      events<-rbind(events, data.frame(time=time, eta=1, S=nS, I=nI))
    } else {
      nI <- max(0, nI - 1)  # Recovery event
      nR <- nR + 1
      events<-rbind(events, data.frame(time=time, eta=0, S=nS, I=nI))
    }
  }  ##  END of while
  
  return(events=events)
}