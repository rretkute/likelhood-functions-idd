library(deSolve)
library(ggplot2)
library(gridExtra)
library(viridis)
library(matrixStats)
library(statip)

source("Leprosy_ODE.R")
source("Leprosy_stoch.R")
source("Exact_likelihood_Leprosy.R")
source('Poisson_likelihood.R')

# Set up
pal <- rev(viridis_pal(option="C")(100))

N0 <- 1000
I0 <- 1
parameters <- data.frame(beta = 0.1, mu = 1/(60*52),
                         gamma = 1/20, pE = 0.01, sigma = 0.6,
                         rho = 1/12, theta = 0, N=N0)
true.param <- data.frame(beta = 0.1, pE = -2)
tmax <- 52*8 # 8 years

beta.cand<-seq(0.01, 0.25, l = 100)
pE.cand<-seq(-5, 0, l = 100)  # log10 scale

###################################
#   Figure 3
###################################
start_time = Sys.time()

# Run deterministic model
xstart <- c(S = N0-I0, E = 0, I = I0, R = 0, D = 0)
times <- seq(from = 0,to = tmax, by = 1) 
out_det<- as.data.frame(ode(xstart, times, Leprocy_ODE_model, 
                            parameters, method=rk4))

ans<-rbind(data.frame(time=out_det$time, Number=out_det$S, Compartment="S"),
           data.frame(time=out_det$time, Number=out_det$E, Compartment="E"),
           data.frame(time=out_det$time, Number=out_det$I, Compartment="I"),
           data.frame(time=out_det$time, Number=out_det$D, Compartment="D"),
           data.frame(time=out_det$time, Number=out_det$R, Compartment="R")
)
ans$Compartment<-factor(ans$Compartment, levels = c("S", "E", "I", "D", "R"))
fig1<-ggplot(ans) + 
  geom_line( aes(x=time/52, y=Number, col=Compartment), size=1) +
  theme_bw() + ggtitle("(a)") +  
  xlab("Time (years)")+
  scale_color_manual(values=c("darkgreen", "gray", "darkred", "darkblue", "black"), name="")

# Run stochastic model (to produce observation data)
initial_state <- data.frame(S = N0 - I0, I = I0, E = 0, D = 0, R = 0)
set.seed(1)
out_stoch<-Leprocy_stochastic_model(parameters, initial_state,  tmax)

ans<-rbind(data.frame(time=out_stoch$time, Number=out_stoch$S, Compartment="S"),
           data.frame(time=out_stoch$time, Number=out_stoch$E, Compartment="E"),
           data.frame(time=out_stoch$time, Number=out_stoch$I, Compartment="I"),
           data.frame(time=out_stoch$time, Number=out_stoch$D, Compartment="D"),
           data.frame(time=out_stoch$time, Number=out_stoch$R, Compartment="R")
)
ans$Compartment<-factor(ans$Compartment, levels = c("S", "E", "I", "D", "R"))
fig2<-ggplot(ans) + 
  geom_line( aes(x=time/52, y=Number, col=Compartment), size=1) +
  theme_bw() + ggtitle("(b)") +  
  xlab("Time (years)")+
  scale_color_manual(values=c("darkgreen", "gray", "darkred", "darkblue", "black"), name="")
fig2A<-ggplot(ans[ans$Compartment=="D",]) + 
  geom_line( aes(x=time/52, y=Number, col=Compartment), size=1, col="darkblue") +
  theme_bw() + xlab("") +ylab("")
fig2<-fig2 + 
  annotation_custom(
    ggplotGrob(fig2A), 
    xmin = 3, xmax = 8, ymin = 650, ymax = 1060
  )
fig2<-fig2+
  geom_line(data=ans[ans$Compartment=="R",], aes(x=time/52, y=Number),
            size=1, col="black") +
  theme_bw() + xlab("") +ylab("")

# Process based likelihood
ans1<-data.frame(beta=c(), pE=c(), ll=c())
for(j in 1:length(beta.cand)){
  for(k in 1:length(pE.cand)){
    prm <- data.frame(beta = beta.cand[j], mu = 1/(60*52),
                      gamma = 1/20, pE = 10^pE.cand[k], sigma = 0.6,
                      rho = 1/12, theta = 0, N=N0)
    ll<-lll_exact_Leprosy(out_stoch, prm)
    ans1<-rbind(ans1, data.frame(beta=beta.cand[j],
                                 pE=pE.cand[k], ll=ll))
  }
}

fig3<-ggplot(ans1) +
  geom_tile(aes(x = beta, y = pE, fill = ll)) +
  scale_fill_gradientn(colors = pal, 
                       limits=c(as.numeric(quantile(ans1$ll, 0.95)) , max(ans1$ll)), 
                       na.value = "gray") +
  geom_point(data = true.param, 
             mapping = aes(x = beta, y = pE),
             shape = 8, col='white', size=2) +
  xlab(expression(beta)) + ylab("log10(pE)") +
  labs(fill='log likelihood') + ggtitle("(c)") +
  theme_bw()

#  Observation based likelihood
obs.dat<-out_stoch$D  
ans2<-data.frame(beta=c(), pE=c(), ll=c())
for(j in 1:length(beta.cand)){
  for(k in 1:length(pE.cand)){
    prm <- data.frame(beta = beta.cand[j], mu = 1/(60*52),
                      gamma = 1/20, pE = 10^pE.cand[k], sigma = 0.6,
                      rho = 1/12, theta = 0, N=N0)
    out_det<- as.data.frame(ode(xstart, out_stoch$time, Leprocy_ODE_model, 
                                prm, method=rk4))
    sims<-out_det$D
    ll<- lll_obs_pois(sims, obs.dat)
    ans2<-rbind(ans2, data.frame(beta=beta.cand[j],
                                 pE=pE.cand[k], ll=ll))
  }
}

fig4<-ggplot(ans2) +
  geom_tile(aes(x = beta, y =pE, fill = ll)) +
  scale_fill_gradientn(colors = pal, 
                       limits=c(as.numeric(quantile(ans2$ll, 0.95)) , max(ans2$ll)), 
                       na.value = "gray") +
  geom_point(data = true.param, 
             mapping = aes(x = beta, y = pE),
             shape = 8, col='white', size=2) +
  xlab(expression(beta)) + ylab("log10(pE)") +
  labs(fill='log likelihood') + ggtitle("(d)") +
  theme_bw()

f<-grid.arrange(fig1, fig2, fig3, fig4, ncol = 2)

ggsave("Figure_3.tiff",plot = f, 
       width = 8, height = 6)

end_time = Sys.time()
print(end_time - start_time)
# Time difference of 1.046644 hours
