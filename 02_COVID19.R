library(deSolve)
library(ggplot2)
library(gridExtra)
library(viridis)
library(matrixStats)
library(statip)

source("R/run_stochastic_model.R")
source("R/COVID19_stoch.R")
source("R/Find_exact_likelihoods.R")
source("R/Exact_likelihood_COVID19.R")

source('R/find_likelihoods.R')
source('R/Poisson_likelihood.R')
source('R/COVID_ODE.R')
source('R/run_ODE_model.R')


pal <- rev(viridis_pal(option="C")(100))

dt<-1 # Aggregate data to a daily cases


# 1. Make synthetic data
I0 <- 1
N <- 3700
initial_state <- data.frame(S = N - I0, I1 = I0, I2 = 0, R = 0, D = 0)
parameters <- data.frame(beta = 0.21, p = 20/712, gamma_1 = 1/10.63,
                         gamma_2 = 1/18.70, alpha = 0.007923033,
                         N = N)
tmax <- 180
set.seed(100)
sim<-run_stochastic_model (parameters, initial_state, tmax, dt, "COVID19")

events<-sim$events
output<-sim$pop


# 2. Plot synthetic data
# eta1: 1 if S->I, 0 if I->R/D
# eta2: 1 if infection is severe, 0 is mild
# eta3: 1 if individual died, 0 if removed
id<-1
tmp<-data.frame(id=c(), x=c(), y=c(), val=c())
for(i in 1:nrow(events)){
  if(events$eta1[i]==1 & events$eta2[i]==0){ # mild infection
    h<-0.5+sum(events$eta1[events$time<=events$time[i]])/(sum(events$eta1))
    tmp<-rbind(tmp, data.frame(id=id, x=events$time[i], y=0, val=1))
    tmp<-rbind(tmp, data.frame(id=id, x=events$time[i], y=h, val=1))
    id<-id+1
  }
  if(events$eta1[i]==1 & events$eta2[i]==1){ # severe infection
    h<-0.5+sum(events$eta1[events$time<=events$time[i]])/(sum(events$eta1))
    tmp<-rbind(tmp, data.frame(id=id, x=events$time[i], y=0, val=2))
    tmp<-rbind(tmp, data.frame(id=id, x=events$time[i], y=h, val=2))
    id<-id+1
  }
  if(events$eta1[i]==0 & events$eta3[i]==0){ # recovery
    h<-0.5+sum(1-events$eta1[events$time<=events$time[i]])/(sum(1-events$eta1))
    tmp<-rbind(tmp, data.frame(id=id, x=events$time[i], y=0, val=3))
    tmp<-rbind(tmp, data.frame(id=id, x=events$time[i], y=-h, val=3))
    id<-id+1
  }
  if(events$eta1[i]==0 & events$eta3[i]==1){ # death
    h<-0.5+sum(1-events$eta1[events$time<=events$time[i]])/(sum(1-events$eta1))
    tmp<-rbind(tmp, data.frame(id=id, x=events$time[i], y=0, val=4))
    tmp<-rbind(tmp, data.frame(id=id, x=events$time[i], y=-h, val=4))
    id<-id+1
  }
}

tmp$val<-as.factor(tmp$val)

f1<-ggplot(tmp) + geom_line(aes(x=x, y=y, group=id, col=val),show.legend = FALSE) +
  scale_color_manual(values=c("red", "yellow", "blue", "black"))+
  geom_line(data=tmp[which(tmp$val==4),],aes(x=x, y=y, group=id), col='black')+
  theme(legend.position = "none") + theme_bw() +
  xlab("time") + ylab("Event") + ggtitle("(a)") +
  scale_y_continuous(breaks = c(-1.5, 1.5), labels=c("Recovery/Death", "Infection")) 


ans<-rbind(data.frame(time=output$time, Number=output$S, Compartment="S"),
           data.frame(time=output$time, Number=output$I1, Compartment="I1"),
           data.frame(time=output$time, Number=output$I2, Compartment="I2"),
           data.frame(time=output$time, Number=output$R, Compartment="R"),
           data.frame(time=output$time, Number=output$D, Compartment="D"))

f2<-ggplot(ans) + geom_line( aes(x=time, y=Number, col=Compartment), size=1) +
  theme_bw() + ggtitle("(b)") +
  scale_color_manual(values=c("green", "red", "yellow","blue", "black"))


# 3. Exact likelihood
true.param <- data.frame(beta = 0.21, alpha = 0.007923033)

beta.cand<-seq(0.01, 1, l = 100)
duration.cand<-seq(1, 200, l = 100)
alpha.cand<-1/duration.cand

start_time = Sys.time()
ll_exact<-find_exact_likelihoods(events, 'COVID19', beta.cand, 
                                 alpha.cand = alpha.cand, N=N,
                                 other_pars = data.frame(p = 20/712, gamma_1 = 1/10.63, gamma_2 = 1/18.70))

end_time = Sys.time()
print(end_time - start_time)

ll_exact$range

#MLE
ll_exact$ans$beta[which(ll_exact$ans$ll==max(ll_exact$ans$ll))]*N
1/ll_exact$ans$alpha[which(ll_exact$ans$ll==max(ll_exact$ans$ll))]

as.numeric(quantile(ll_exact$ans$ll, probs=seq(0,1,0.1)))

f3<-ggplot(data = ll_exact$ans) +
  geom_tile(aes(x = beta, y = 1/alpha, fill = ll)) +
  scale_fill_gradientn(colors = pal, limits=c(15658.361, ll_exact$range[2]), 
                       na.value = "gray") +
  geom_point(data = true.param, 
             mapping = aes(x = beta, y = 1/alpha), shape = 8, col='white', size=4) +
  xlab(expression(beta)) + ylab("Average time to death (days)") +
  labs(fill='log likelihood') + ggtitle("(c)")

# 4. Poisson likelihood
initial_state <- c(S = N - I0, I1 = I0, I2 = 0, R = 0, D = 0)

start_time = Sys.time()
ll_pois <- find_poisson_likelihoods(output, 'COVID19', initial_state, tmax, dt,
                                    beta.cand, alpha.cand = alpha.cand, N=N,
                                    other_pars = data.frame(p = 20/712, gamma_1 = 1/10.63, gamma_2 = 1/18.70))
end_time = Sys.time()
print(end_time - start_time)

ll_pois$range

#MLE
ll_pois$ans$beta[which(-ll_pois$ans$ll==sort(-ll_pois$ans$ll)[1])]
1/(ll_pois$ans$alpha[which(-ll_pois$ans$ll==sort(-ll_pois$ans$ll)[1])])
as.numeric(quantile(ll_pois$ans$ll[ll_pois$ans$ll>-Inf], probs=seq(0,1,0.1), na.rm=TRUE))


f4<-ggplot(data = ll_pois$ans) +
  geom_tile(aes(x = beta, y = 1/alpha, fill = ll)) +
  scale_fill_gradientn(colors = pal, limits=c(-46.12182, ll_pois$range[2]), 
                       na.value = "gray") +
  geom_point(data = true.param, 
             mapping = aes(x = beta, y = 1/alpha), shape = 8, col='white', size=4) +
  xlab(expression(beta)) + ylab("Average time to death (days)") +
  labs(fill='log likelihood') + ggtitle("(d)")



f <- grid.arrange(f1, f2, f3, f4,  ncol = 2)

ggsave("covid_1.png",plot = f, 
       width = 12, height = 8)