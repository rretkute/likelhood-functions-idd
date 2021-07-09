library(deSolve)
library(ggplot2)
library(gridExtra)
library(viridis)
library(matrixStats)
library(statip)

source("R/run_stochastic_model.R")
source("R/Flu_stoch.R")
source("R/COVID19_stoch.R")

source('R/NegBinom_likelihood.R')
source('R/NB_likelihood.R')
source('R/Flu_ODE.R')
source('R/COVID_ODE.R')
source('R/run_ODE_model.R')


pal <- rev(viridis_pal(option="C")(100))

dt<-1 # Aggregate data to a daily cases


###########################
###        Flu   ######
##########################

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

output<-sim$pop

true.param <- data.frame(beta = beta, gamma = gamma)

beta.cand<-seq(0.1, 10, by=0.1)
inf.per.cand<-seq(0.1, 10, by=0.1)


gamma.cand<-1/inf.per.cand

initial_state <- c(S = N - I0, I = I0, R = 0)



start_time = Sys.time()
ll_NB <- find_NB_likelihoods(output, k=0.57, 'flu', initial_state, tmax, dt,
                                    beta.cand=beta.cand, gamma.cand = gamma.cand, N=N)
end_time = Sys.time()
print(end_time - start_time)

ll_NB$range

#MLE
ll_NB$ans$beta[which(-ll_NB$ans$ll==sort(-ll_NB$ans$ll)[1])]
1/(ll_NB$ans$gamma[which(-ll_NB$ans$ll==sort(-ll_NB$ans$ll)[1])])
as.numeric(quantile(ll_NB$ans$ll[ll_NB$ans$ll>-Inf], probs=seq(0,1,0.1), na.rm=TRUE))

fS1<-ggplot(data = ll_NB$ans) +
  geom_tile(aes(x =beta, y = 1/gamma, fill = ll)) +
  scale_fill_gradientn(colors = pal, limits=c(-128.22653, ll_NB$range[2]), 
                       na.value = "gray") +
  geom_point(data = true.param, 
             mapping = aes(x = beta, y = 1/gamma), shape = 8, col='white', size=4) +
  xlab(expression(beta)) + ylab("Recovery time (days)") +
  labs(fill='log likelihood') + ggtitle("(a)")


###########################
###        COVID19  
##########################

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

output<-sim$pop

true.param <- data.frame(beta = 0.21, alpha = 0.007923033)

beta.cand<-seq(0.01, 1, l = 100) 
duration.cand<-seq(1, 200, l = 100)
alpha.cand<-1/duration.cand


initial_state <- c(S = N - I0, I1 = I0, I2 = 0, R = 0, D = 0)
start_time = Sys.time()
ll_NB <- find_NB_likelihoods(output, k=0.57, 'COVID19', initial_state, tmax, dt,
                                    beta.cand, alpha.cand = alpha.cand,
                                    other_pars = data.frame(p = 0.01, gamma_1 = 1/12.56,
                                                            gamma_2 = 1/21.2))
end_time = Sys.time()
print(end_time - start_time)

#MLE
ll_NB$ans$beta[which(-ll_NB$ans$ll==sort(-ll_NB$ans$ll)[1])]
1/(ll_NB$ans$alpha[which(-ll_NB$ans$ll==sort(-ll_NB$ans$ll)[1])])

as.numeric(quantile(ll_NB$ans$ll[ll_NB$ans$ll>-Inf], probs=seq(0,1,0.1), na.rm=TRUE))


fS2<-ggplot(data = ll_NB$ans) +
  geom_tile(aes(x = beta, y = 1/alpha, fill = ll)) +
  scale_fill_gradientn(colors = pal, limits=c( -46.66632, ll_NB$range[2]), 
                       na.value = "gray") +
  geom_point(data = true.param, 
             mapping = aes(x = beta, y = 1/alpha), shape = 8, col='white', size=4) +
  xlab(expression(beta)) + ylab("Average time to death (days)") +
  labs(fill='log likelihood') + ggtitle("(b)")

f <- grid.arrange(fS1, fS2, ncol = 2)

