library(deSolve)
library(ggplot2)
library(gridExtra)
library(viridis)
library(matrixStats)
library(statip)

source("R/run_stochastic_model.R")
source("R/Flu_stoch.R")
source("R/Find_exact_likelihoods.R")
source("R/Exact_likelihood_flu.R")

source('R/find_likelihoods.R')
source('R/Poisson_likelihood.R')
source('R/Flu_ODE.R')
source('R/run_ODE_model.R')


pal <- rev(viridis_pal(option="C")(100))

dt<-1 # Aggregate data to a daily cases


# 1. Make synthetic data
I0<-1
N<-763
beta<-1.6
gamma<-1/2

initial_state <- data.frame(S = N - I0, I = I0)
parameters <- data.frame(beta=beta, gamma=gamma, N=N)

tmax <- 25
set.seed(100)
sim<-run_stochastic_model (parameters, initial_state, tmax, dt, "flu")

events<-sim$events
output<-sim$pop


# 2. Plot synthetic data 

id<-1
tmp<-data.frame(id=c(), x=c(), y=c(), val=c())
for(i in 1:nrow(events)){
  if(events$eta[i]==1){
    tmp<-rbind(tmp, data.frame(id=id, x=events$time[i], y=0, val=1))
    tmp<-rbind(tmp, data.frame(id=id, x=events$time[i],
                               y=0.5+sum(events$eta[events$time<=events$time[i]])/(sum(events$eta)), val=1))
    id<-id+1
  } else {
    tmp<-rbind(tmp, data.frame(id=id, x=events$time[i], y=0, val=0))
    tmp<-rbind(tmp, data.frame(id=id, x=events$time[i], 
                               y=-0.5-sum(1-events$eta[events$time<=events$time[i]])/(sum(1-events$eta)), val=0))
    id<-id+1
  }
}
tmp$val<-as.factor(tmp$val)

f1<-ggplot(tmp) + geom_line(aes(x=x, y=y, group=id, col=val),show.legend = FALSE) +
  scale_color_manual(values=c("red",  "blue"))+  theme_bw() +
  xlab("time") + ylab("Event") + ggtitle("(a)") +
  scale_y_continuous(breaks = c(-1.5, 1.5), labels=c("Recovery", "Infection")) 

ans<-rbind(data.frame(time=output$time, Number=output$S, Compartment="S"),
           data.frame(time=output$time, Number=output$I, Compartment="I"),
           data.frame(time=output$time, Number=output$R, Compartment="R"))

f2<-ggplot(ans) + geom_line( aes(x=time, y=Number, col=Compartment), size=1) +
  theme_bw() + ggtitle("(b)") +  
  scale_color_manual(values=c("green", "red", "blue"))


# 3. Exact likelihood
true.param <- data.frame(beta = beta, gamma = gamma)

beta.cand<-seq(0.1, 10, by=0.1)
inf.per.cand<-seq(0.1, 10, by=0.1)

# Total number of grid points
gamma.cand<-1/inf.per.cand

start_time = Sys.time()
ll_exact<-find_exact_likelihoods(events, 'flu', beta.cand=beta.cand, gamma.cand = gamma.cand, N=N)
end_time = Sys.time()
print(end_time - start_time)

ll_exact$range

#MLE
ll_exact$ans$beta[which(ll_exact$ans$ll==max(ll_exact$ans$ll))]*N

as.numeric(quantile(ll_exact$ans$ll, probs=seq(0,1,0.1)))

f3<-ggplot(data = ll_exact$ans) +
  geom_tile(aes(x =beta, y = 1/gamma, fill = ll)) +
  scale_fill_gradientn(colors = pal, limits=c(4464.371, ll_exact$range[2]), 
                       na.value = "gray") +
  geom_point(data = true.param, 
             mapping = aes(x = beta, y = 1/gamma), shape = 8, col='white', size=4) +
  xlab(expression(beta)) + ylab("Recovery time (days)") +
  labs(fill='log likelihood') + ggtitle("(c)")

# 4. Poisson likelihood
initial_state <- c(S = N - I0, I = I0, R = 0)

start_time = Sys.time()
ll_pois <- find_poisson_likelihoods(output, 'flu', initial_state, tmax, dt,
                                    beta.cand=beta.cand, gamma.cand = gamma.cand, N=N)
end_time = Sys.time()
print(end_time - start_time)

ll_pois$range

#MLE
ll_pois$ans$beta[which(-ll_pois$ans$ll==sort(-ll_pois$ans$ll)[1])]*N
1/(ll_pois$ans$gamma[which(-ll_pois$ans$ll==sort(-ll_pois$ans$ll)[1])])
as.numeric(quantile(ll_pois$ans$ll[ll_pois$ans$ll>-Inf], probs=seq(0,1,0.1), na.rm=TRUE))


f4<-ggplot(data = ll_pois$ans) +
  geom_tile(aes(x =beta, y = 1/gamma, fill = ll)) +
  scale_fill_gradientn(colors = pal, limits=c(-177.10694, ll_pois$range[2]), 
                       na.value = "gray") +
  geom_point(data = true.param, 
             mapping = aes(x = beta, y = 1/gamma), shape = 8, col='white', size=4) +
  xlab(expression(beta)) + ylab("Recovery time (days)") +
  labs(fill='log likelihood') + ggtitle("(d)")

f <- grid.arrange(f1, f2, f3, f4, ncol = 2)

ggsave("flu.png",plot = f, 
       width = 12, height = 8)