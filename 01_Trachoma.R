library(deSolve)
library(ggplot2)
library(gridExtra)
library(viridis)
library(matrixStats)
library(statip)

source("Trachoma_ODE.R")
source("Trachoma_stoch.R")
source("Exact_likelihood_Trachoma.R")
source('Poisson_likelihood.R')

# Set up
pal <- rev(viridis_pal(option="C")(100))

N0 <- 100
I0 <- 1
parameters <- data.frame(beta = 0.044, gamma = 0.017, N=N0)
tmax <- 52*6 # 6 years

beta.cand<-seq(0.01, 0.11, l = 100)
gamma.cand<-seq(0.005, 0.08, l = 100)


###################################
#   Figure 1
###################################
start_time = Sys.time()

# Run stochastic model (to produce observation data)
initial_state <- data.frame(S=N0-I0, I = I0)
set.seed(50556)
out_stoch<-Trachoma_stochastic_model(parameters, initial_state, 0, tmax)

# Epidemic curves
ans<-rbind(data.frame(time=out_stoch$time, Number=out_stoch$S, Compartment="S"),
           data.frame(time=out_stoch$time, Number=out_stoch$I, Compartment="I"))

# Individual events
id<-1
tmp<-data.frame(id=c(), x=c(), y=c(), val=c())
for(i in 1:nrow(out_stoch)){
  if(out_stoch$eta[i]==1){
    tmp<-rbind(tmp, data.frame(id=id, x=out_stoch$time[i], y=0, val=1))
    tmp<-rbind(tmp, data.frame(id=id, x=out_stoch$time[i],
                y=0.5+sum(out_stoch$eta[out_stoch$time<=out_stoch$time[i]])/(sum(out_stoch$eta)), val=1))
    id<-id+1
  } else {
    tmp<-rbind(tmp, data.frame(id=id, x=out_stoch$time[i], y=0, val=0))
    tmp<-rbind(tmp, data.frame(id=id, x=out_stoch$time[i], 
                               y=-0.5-sum(1-out_stoch$eta[out_stoch$time<=out_stoch$time[i]])/(sum(1-out_stoch$eta)), val=0))
    id<-id+1
  }
}
tmp$val<-as.factor(tmp$val)

fig1<-ggplot(tmp) + geom_line(aes(x=x/52, y=y, group=id, col=val),show.legend = FALSE) +
  scale_color_manual(values=c("red",  "blue"))+  theme_bw() +
  xlab("Time (years)") + ylab("Event") + ggtitle("(a)") +
  scale_y_continuous(breaks = c(-1.5, 1.5), labels=c("Recovery", "Infection")) 

# Run deterministic model (to compare output)
xstart <- c(I = I0)
times <- seq(from = 0,to = tmax, by = 1) 
out_det<- as.data.frame(ode(xstart, times, Trachoma_ODE_model, 
                            parameters, method=rk4))

fig2<-ggplot(ans) + 
  geom_line( aes(x=time/52, y=Number, col=Compartment), size=1) +
  geom_line(data=out_det, aes(x=time/52, y=I), col="red", size=1)+
  geom_line(data=out_det, aes(x=time/52, y=N0-I), col="green", size=1)+
  theme_bw() + ggtitle("(b)") +  
  xlab("Time (years)")+
  scale_color_manual(values=c("darkred", "darkgreen"), name="")

# Process based likelihood
ans1<-data.frame(beta=c(), beta=c(), ll=c())
for(j in 1:length(beta.cand)){
  for(k in 1:length(gamma.cand)){
    prm <- data.frame(beta = beta.cand[j], gamma = gamma.cand[k], N=N0)
    ll<-lll_exact_Trachoma(out_stoch, prm)
    ans1<-rbind(ans1, data.frame(beta=beta.cand[j],
                                 gamma=gamma.cand[k], ll=ll))
  }
}

fig3<-ggplot(ans1) +
  geom_tile(aes(x = beta, y = gamma, fill = ll)) +
  scale_fill_gradientn(colors = pal, 
                       limits=c(as.numeric(quantile(ans1$ll, 0.95)) , max(ans1$ll)), 
                       na.value = "gray") +
  geom_point(data = parameters, 
             mapping = aes(x = beta, y = gamma),
             shape = 8, col='white', size=2) +
  xlab(expression(beta)) + ylab(expression(gamma)) +
  labs(fill='log likelihood') + ggtitle("(c)") +
  theme_bw()

#  Observation based likelihood
obs.dat<-out_stoch$I  
ans2<-data.frame(beta=c(), pE=c(), ll=c())
for(j in 1:length(beta.cand)){
  for(k in 1:length(gamma.cand)){
    prm <- data.frame(beta = beta.cand[j], gamma = gamma.cand[k], N=N0)
    out_det<- as.data.frame(ode(xstart, out_stoch$time, Trachoma_ODE_model, 
                                prm, method=rk4))
    sims<-out_det$I
    ll<- lll_obs_pois(sims, obs.dat)
    ans2<-rbind(ans2, data.frame(beta=beta.cand[j],
                                 gamma=gamma.cand[k], ll=ll))
  }
}

fig4<-ggplot(ans2) +
  geom_tile(aes(x = beta, y =gamma, fill = ll)) +
  scale_fill_gradientn(colors = pal, 
                       limits=c(as.numeric(quantile(ans2$ll, 0.95)) , max(ans2$ll)), 
                       na.value = "gray") +
  geom_point(data = parameters, 
             mapping = aes(x = beta, y = gamma),
             shape = 8, col='white', size=2) +
  xlab(expression(beta)) + ylab(expression(gamma)) +
  labs(fill='log likelihood') + ggtitle("(d)") +
  theme_bw()

f<-grid.arrange(fig1, fig2, fig3, fig4, ncol = 2)

ggsave("Figure_1.tiff",plot = f, 
       width = 8, height = 6)

end_time = Sys.time()
print(end_time - start_time)
# Time difference of 4.472006 mins

###################################
#   Figure 2
###################################
start_time = Sys.time()

# Effectiveness of MDA
eff<-0.8

# Process based likelihood estimates
elm.prm1<-data.frame(id=c(), prev=c(), elim=c())
wh<-which(ans1$ll>=as.numeric(quantile(ans1$ll, 0.99)))
for(ii in wh){
  prm <- data.frame(beta = ans1$beta[ii], gamma = ans1$gamma[ii], N=N0)
  Prev<-c()
  Elim.time<-c()
  for(jj in 1:100){
    tmax <- 52*6
    initial_state <- data.frame(S=N0-I0, I = I0)
    out_stoch <- Trachoma_stochastic_model(prm, initial_state,0,tmax)
    prev<-tail(out_stoch$I,1)
    if(prev>20){
      tk<- seq(max(out_stoch$time), 52*100, 26)
      i<-2
      elim<-0
      elim.time<- 100
      while(elim==0){
        new.inf<- length(which(runif(tail(out_stoch$I,1))<=(1-eff)))
        initial_state <- data.frame(S=N0-new.inf, I = new.inf)
        out_stoch<- rbind(out_stoch,
                          as.data.frame(Trachoma_stochastic_model(prm, initial_state, tk[i-1], tk[i])))
        if(tail(out_stoch$I,1)<1) {
          elim<-1
          elim.time<- (i-1)*0.5
        }
        if(i>=100){
          elim<-1
        }
        i<-i+1
      }
      Prev<-c(Prev, prev)
      Elim.time<-c(Elim.time, elim.time)
    }
  }
  elm.prm1<-rbind(elm.prm1, data.frame(id=ii, prev=mean(Prev), 
                                       elim=mean(Elim.time)))
}

# Observation based likelihood estimates
elm.prm2<-data.frame(id=c(), prev=c(), elim=c())
wh<-which(ans2$ll>=as.numeric(quantile(ans2$ll, 0.99)))

for(ii in wh){
  tmax <- 52*6
  prm <- data.frame(beta = ans2$beta[ii], gamma = ans2$gamma[ii], N=N0)
  times <- seq(from = 0,to = tmax, by = 1) 
  xstart <- c(I = I0)
  out_det<- as.data.frame(ode(xstart, times, Trachoma_ODE_model, 
                              prm, method = "rk4"))
  prev<-tail(out_det$I,1)
  tk<- seq(52*6, 52*100, 26)
  i<-2
  elim<-0
  while(elim==0){
    times <- seq(from = tk[i-1],to = tk[i], by = 1) 
    xstart <- c(I = (1-eff)*tail(out_det$I,1))
    out_det<- rbind(out_det,
                    as.data.frame(ode(xstart, times, Trachoma_ODE_model, 
                                      prm, method = "rk4")))
    if(tail(out_det$I,1)<0.001) {
      elim<-1
      elim.time<- i*0.5
    }
    i<-i+1
  }
  elm.prm2<-rbind(elm.prm2, data.frame(id=ii, prev=prev, elim=elim.time))
}

elm.prm1$Model<-"Stochastic"
elm.prm2$Model<-"Deterministic"

elm.prm<-rbind(elm.prm1, elm.prm2)

f<-ggplot(elm.prm, aes(x=prev, y=elim))+
  geom_point(pch=21, fill="#43a2ca", size=2) +
  facet_grid(.~Model) +
  xlab("Pre-treatment prevalence (%)") + ylab("Years until elimination")+
  theme_bw()

ggsave("Figure_2.tiff",plot = f, 
       width = 8, height = 3)

end_time = Sys.time()
print(end_time - start_time)
# Time difference of 23.39428 mins
