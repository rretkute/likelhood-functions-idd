#' COVID-19 Stochastic model
#'
#' @param parameters 
#' @param initial_state 
#' @param tmax 
#'
#' @return dataframe with event timing, indicators and sizes of compartments
#' @export
#'
#'
COVID19_stochastic_model<-function(parameters, initial_state, tmax){
  
  ## All events are coded using three indexes:
  # eta1: 1 if S->I, 0 if I->R/D
  # eta2: 1 if infection is severe, 0 is mild
  # eta3: 1 if individual died, 0 if removed
  
  # assign numbers in each state for first time step
  nS <- initial_state$S
  nI1 <- initial_state$I1
  nI2 <- initial_state$I2
  nR<- initial_state$R
  nD<- initial_state$D
  #####  RR: otherwise makes these as data.frame 
  
  
  #initialise time
  time <- 0
  
  if(nI1>0) {
    events<-data.frame(time=0, eta1=1, eta2=0, eta3=0,  S=nS, I1=nI1, I2=nI2, R=nR, D=nD)
  }
  if(nI2>0) {
     events<-data.frame(time=0, eta1=1, eta2=1, eta3=0, eta4=0, S=nS, I1=nI1, I2=nI2, R=nR, D=nD)
  }
  
  while((time <= tmax & (nI1+nI2>0))){ 
    r <- c(parameters$beta * nS * (nI1 + nI2)/parameters$N, # S-> I1 or I2
           parameters$gamma_1 * nI1,  # I1-> R 
           parameters$gamma_2 * nI2, # I2-> R 
           parameters$alpha * nI2) # I2->D
    rtotal <- sum(r)
    r<-cumsum(r)
    delta <- (-1/rtotal)*log(runif(1)) 
    time <- time + delta # Time to the next event
    if (runif(1)*rtotal <= r[1]) { # Transmission event
      nS <- max(0, nS - 1)  # Transmission event.
      if(runif(1)>=parameters$p){ # Mild infection
        nI1<- nI1 + 1
        events<-rbind(events, data.frame(time=time, eta1=1, eta2=0, eta3=0, S=nS, I1=nI1, I2=nI2, R=nR, D=nD))
      }  else {
        nI2<-nI2+1
        events<-rbind(events, data.frame(time=time, eta1=1, eta2=1, eta3=0, S=nS, I1=nI1, I2=nI2, R=nR, D=nD))
      }
    }  
    if ((runif(1)*rtotal > r[1]) & (runif(1)*rtotal <= r[2]) & nI1>0) { # I1-> R
      nI1 <- max(0, nI1 - 1)  # Recovery event
      nR <- nR + 1
      events<-rbind(events, data.frame(time=time, eta1=0, eta2=0, eta3=0, S=nS, I1=nI1, I2=nI2, R=nR, D=nD))
    }
    if ((runif(1)*rtotal > r[2]) & (runif(1)*rtotal <= r[3]) & nI2>0) { # I2-> R
      nI2 <- max(0, nI2 - 1)  # Recovery event
      nR <- nR + 1
      events<-rbind(events, data.frame(time=time, eta1=0, eta2=1, eta3=0, S=nS, I1=nI1, I2=nI2, R=nR, D=nD))
    }
    if ((runif(1)*rtotal > r[3]) &  nI2>0) { # I2-> D
      nI2 <- max(0, nI2 - 1)  # Recovery event
      nD <- nD + 1
      events<-rbind(events, data.frame(time=time, eta1=0, eta2=1, eta3=1, S=nS, I1=nI1, I2=nI2, R=nR, D=nD))
    } 
  }  ##  END of while
  
  return(events)
}
