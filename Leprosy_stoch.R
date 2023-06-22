#' Leprocy Stochastic model
#'
#' @param parameters 
#' @param initial_state 
#' @param tmax 
#'
#' @return dataframe with event timing, indicators and sizes of compartments
#' @export
#'
#' @examples
Leprocy_stochastic_model<-function(parameters, initial_state, tmax){
  
  ## All events are coded using 3 indexes:
  # eta1: 1 if S->E, 0 otherwise
  # eta2: 1 if becomes infectious and detected, 0 if becomes infectious
  # eta3: 1 if recovers, 0 otherwise
  # eta4: 1 if death, 0 otherwise
  
  # assign numbers in each state for first time step
  nS <- initial_state$S
  nE <- initial_state$E
  nI <- initial_state$I
  nD<- initial_state$D
  nR<- initial_state$R
  N<-nS+nE+nI+nD+nR
 
  #initialise time
  time <- 0
  
  if(nI>0) {
    events<-data.frame(time=0, eta1=1, eta2=0, eta3=0, eta4=0, 
                       S=nS, E=nE, I=nI, R=nR, D=nD, nD=0)
  }
  
  while((time <= tmax & (nI+nD+nE>0))){ 
    r <- c(parameters$beta * nS * (nI + nD)/N, # S-> E
           parameters$sigma * nE,  # E->D/I 
           (parameters$rho+parameters$gamma) * nD, # D->R
           parameters$gamma * nI, # I-> R 
           parameters$mu * nE, # Death
           parameters$mu * nI, # Death
           parameters$mu * nD, # Death
           parameters$mu * nR # Death
    )
    rtotal <- sum(r)
    r<-cumsum(r)
    delta <- (-1/rtotal)*log(runif(1)) 
    time <- time + delta # Time to the next event
    rnd<-runif(1)*rtotal
    if (rnd <= r[1]) { 
      # S->E
      nS <- max(0, nS - 1) 
      nE<-nE+1
      events<-rbind(events, data.frame(time=time, eta1=1, eta2=0, eta3=0, eta4=0,
                                       S=nS, E=nE, I=nI, R=nR, D=nD, nD=0))
    }  
    if ((rnd > r[1]) & (rnd <= r[2]) & nE>0) { 
      # E-> D/I
      nE <- max(0, nE - 1)  
      if(runif(1)<=parameters$pE){ 
        # E->D
        nD <- nD + 1
        events<-rbind(events, data.frame(time=time, eta1=0, eta2=1, eta3=0, eta4=0,
                                         S=nS, E=nE, I=nI, R=nR, D=nD, nD=1))
      } else { 
       # E-> I
        nI <- nI + 1
        events<-rbind(events, data.frame(time=time, eta1=0, eta2=0, eta3=0, eta4=0,
                                         S=nS, E=nE, I=nI, R=nR, D=nD, nD=0))
      }
    }
    if ((rnd > r[2]) & (rnd <= r[3]) &  nD>0) { 
      #  D -> R
      nD <- max(0, nD - 1)  
      nR <- nR + 1
      events<-rbind(events, data.frame(time=time, eta1=0, eta2=1, eta3=1, eta4=0,
                                       S=nS, E=nE, I=nI, R=nR, D=nD, nD=0))
    }
    if ((rnd > r[3]) & (rnd <= r[4]) &  nI>0) { 
      #  I -> R
      nI <- max(0, nI - 1)  
      nR <- nR + 1
      events<-rbind(events, data.frame(time=time, eta1=0, eta2=0, eta3=1,eta4=0,
                                       S=nS, E=nE, I=nI, R=nR, D=nD, nD=0))
    }
    if ((rnd > r[4]) & (rnd <= r[5]) &  nE>0) { 
      # Death E -> Birth S
      nE <- max(0, nE - 1)  
      nS <- nS + 1
      events<-rbind(events, data.frame(time=time, eta1=1, eta2=0, eta3=0, eta4=1, 
                                       S=nS, E=nE, I=nI, R=nR, D=nD, nD=0))
    } 
    if ((rnd > r[5]) & (rnd <= r[6]) &  nD>0) { 
      # Death D -> Birth S
      nD <- max(0, nD - 1)  
      nS <- nS + 1
      events<-rbind(events, data.frame(time=time, eta1=0, eta2=1, eta3=0, eta4=1,
                                       S=nS, E=nE, I=nI, R=nR, D=nD, nD=0))
    } 
    if ((rnd > r[6]) & (rnd <= r[7]) &  nI>0) { 
      # Death I -> Birth S
      nI <- max(0, nI - 1)  
      nS <- nS + 1
      events<-rbind(events, data.frame(time=time, eta1=0, eta2=0, eta3=0,
                  eta4=1, S=nS, E=nE, I=nI, R=nR, D=nD, nD=0))
    } 
    if ((rnd > r[7]) &  nR>0) { 
      # Death R -> Birth S
      nR <- max(0, nR - 1)  
      nS <- nS + 1
      events<-rbind(events, data.frame(time=time, eta1=0, eta2=0, eta3=1, eta4=1,
                                       S=nS, E=nE, I=nI, R=nR, D=nD, nD=0))
    } 
  }  ##  END of while
  
  return(events)
}


