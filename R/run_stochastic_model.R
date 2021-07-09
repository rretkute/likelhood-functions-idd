#' Run stochastic model
#'
#' @param parameters 
#' @param initial_state 
#' @param tmax Outbreak observation time (days)
#' @param dt Time aggregation interval (days)
#' @param model 
#'
#' @return Infectious events and data aggregated at a population level
#' @export
#'
#' @examples
run_stochastic_model <- function(parameters, initial_state, tmax, dt, model){
  N<-parameters$N
  times<-seq(0,tmax, by=dt)
  
  if(model=="flu"){
    # eta: 1 if S->I, 0 if I->R
    out <- Flu_stochastic_model(parameters, initial_state,  tmax)
    # Number of S and I at dt day interval
    ans.S<-sapply(1:length(times), function(a) N-sum(out$eta[which(out$time<=times[a])]))
    ans.I<-sapply(1:length(times), function(a) sum(out$eta[which(out$time<=times[a])])-
                    sum(1-out$eta[which(out$time<=times[a])]))
    ans.R<-sapply(1:length(times), function(a) sum(1-out$eta[which(out$time<=times[a])]))
    
    return(list(events=out, pop=data.frame(time=times, S=ans.S, I=ans.I, R=ans.R)))
  }
  
  if(model=="COVID19"){
    # eta1: 1 if S->I, 0 if I->R/D
    # eta2: 1 if infection is severe, 0 is mild
    # eta3: 1 if individual died, 0 if removed
    out <- COVID19_stochastic_model(parameters, initial_state,  tmax)
    # Number of S, I1, I2, R and D at dt day interval
    ans.S<-sapply(1:length(times), function(a) N-sum(out$eta1[which(out$time<=times[a])]))
    ans.I1<-sapply(1:length(times), function(a) sum(out$eta1[which(out$time<=times[a])]*
                                              (1-out$eta2[which(out$time<=times[a])]))-
    sum((1-out$eta1[which(out$time<=times[a])])* (1-out$eta2[which(out$time<=times[a])])))
    ans.I2<-sapply(1:length(times), function(a) sum(out$eta1[which(out$time<=times[a])]*
                                                    out$eta2[which(out$time<=times[a])])-
                   sum((1-out$eta1[which(out$time<=times[a])])* out$eta2[which(out$time<=times[a])]))
    ans.R<-sapply(1:length(times), function(a) sum((1-out$eta1[which(out$time<=times[a])])*
                                                   (1-out$eta3[which(out$time<=times[a])])))
    ans.D<-sapply(1:length(times), function(a) sum((1-out$eta1[which(out$time<=times[a])])*
                                                   out$eta3[which(out$time<=times[a])]))
    return(list(events=out, pop=data.frame(time=times, S=ans.S, I1=ans.I1, I2=ans.I2, 
                             R=ans.R, D= ans.D)))
  }
   if(!(model %in% c("COVID19", "flu"))) print("Please specify correct model")
}
