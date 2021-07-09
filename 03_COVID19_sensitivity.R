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

###########################
###        COVID19  - low values
##########################

# 1. Make synthetic data
I0 <- 1
N <- 3700
initial_state <- data.frame(S = N - I0, I1 = I0, I2 = 0, R = 0, D = 0)
parameters <- data.frame(beta = 0.21, p = 0.01, gamma_1 = 1/12.56,
                         gamma_2 = 1/21.2, alpha = 0.007923033,
                         N = N)
tmax <- 180
set.seed(100)
sim<-run_stochastic_model (parameters, initial_state, tmax, dt, "COVID19")

events<-sim$events
output<-sim$pop


# 3. Exact likelihood
true.param <- data.frame(beta = 0.21, alpha = 0.007923033)

beta.cand<-seq(0.01, 1, l = 100) 
duration.cand<-seq(1, 200, l = 100)
alpha.cand<-1/duration.cand

ll_exact<-find_exact_likelihoods(events, 'COVID19', beta.cand, 
                                 alpha.cand = alpha.cand, N=N,
                                 other_pars = data.frame(p = 0.01, gamma_1 = 1/12.56,
                                                         gamma_2 = 1/21.2))

ll_exact$range

as.numeric(quantile(ll_exact$ans$ll, probs=seq(0,1,0.1)))

f3A<-ggplot(data = ll_exact$ans) +
  geom_tile(aes(x = beta, y = 1/alpha, fill = ll)) +
  scale_fill_gradientn(colors = pal, limits=c(18399.39, ll_exact$range[2]), 
                       na.value = "gray") +
  geom_point(data = true.param, 
             mapping = aes(x = beta, y = 1/alpha), shape = 8, col='white', size=4) +
  xlab(expression(beta)) + ylab("Average time to death (days)") +
  labs(fill='log likelihood') + ggtitle("(a)")

# 4. Poisson likelihood
initial_state <- c(S = N - I0, I1 = I0, I2 = 0, R = 0, D = 0)
ll_pois <- find_poisson_likelihoods(output, 'COVID19', initial_state, tmax, dt,
                                    beta.cand, alpha.cand = alpha.cand, N=N,
                                    other_pars = data.frame(p = 0.01, gamma_1 = 1/12.56,
                                                            gamma_2 = 1/21.2))


ll_pois$range
as.numeric(quantile(ll_pois$ans$ll[ll_pois$ans$ll>-Inf], probs=seq(0,1,0.1), na.rm=TRUE))


f4A<-ggplot(data = ll_pois$ans) +
  geom_tile(aes(x = beta, y = 1/alpha, fill = ll)) +
  scale_fill_gradientn(colors = pal, limits=c(-28.31010, ll_pois$range[2]), 
                       na.value = "gray") +
  geom_point(data = true.param, 
             mapping = aes(x = beta, y = 1/alpha), shape = 8, col='white', size=4) +
  xlab(expression(beta)) + ylab("Average time to death (days)") +
  labs(fill='log likelihood') + ggtitle("(b)")


###########################
###        COVID19  - high values
##########################

# 1. Make synthetic data
I0 <- 1
N <- 3700
initial_state <- data.frame(S = N - I0, I1 = I0, I2 = 0, R = 0, D = 0)
parameters <- data.frame(beta = 0.21, p = 0.184, gamma_1 = 1/8.7,
                         gamma_2 = 1/16.2, alpha = 0.007923033,
                         N = N)
tmax <- 180
set.seed(100)
sim<-run_stochastic_model (parameters, initial_state, tmax, dt, "COVID19")

events<-sim$events
output<-sim$pop


# 3. Exact likelihood
true.param <- data.frame(beta = 0.21, alpha = 0.007923033)

beta.cand<-seq(0.01, 1, l = 100)
duration.cand<-seq(1, 200, l = 100)
alpha.cand<-1/duration.cand

ll_exact<-find_exact_likelihoods(events, 'COVID19', beta.cand, 
                                 alpha.cand = alpha.cand, N=N,
                                 other_pars = data.frame(p = 0.184, gamma_1 = 1/8.7,
                                                         gamma_2 = 1/16.2))

ll_exact$range
as.numeric(quantile(ll_exact$ans$ll, probs=seq(0,1,0.1)))

f3B<-ggplot(data = ll_exact$ans) +
  geom_tile(aes(x = beta, y = 1/alpha, fill = ll)) +
  scale_fill_gradientn(colors = pal, limits=c(10064.004, ll_exact$range[2]), 
                       na.value = "gray") +
  geom_point(data = true.param, 
             mapping = aes(x = beta, y = 1/alpha), shape = 8, col='white', size=4) +
  xlab(expression(beta)) + ylab("Average time to death (days)") +
  labs(fill='log likelihood') + ggtitle("(c)")

# 4. Poisson likelihood
initial_state <- c(S = N - I0, I1 = I0, I2 = 0, R = 0, D = 0)
ll_pois <- find_poisson_likelihoods(output, 'COVID19', initial_state, tmax, dt,
                                    beta.cand, alpha.cand = alpha.cand, N=N,
                                    other_pars = data.frame(p = 0.184, gamma_1 = 1/8.7,
                                                            gamma_2 = 1/16.2))


ll_pois$range
as.numeric(quantile(ll_pois$ans$ll[ll_pois$ans$ll>-Inf], probs=seq(0,1,0.1), na.rm=TRUE))


f4B<-ggplot(data = ll_pois$ans) +
  geom_tile(aes(x = beta, y = 1/alpha, fill = ll)) +
  scale_fill_gradientn(colors = pal, limits=c( -166.0913, ll_pois$range[2]), 
                       na.value = "gray") +
  geom_point(data = true.param, 
             mapping = aes(x = beta, y = 1/alpha), shape = 8, col='white', size=4) +
  xlab(expression(beta)) + ylab("Average time to death (days)") +
  labs(fill='log likelihood') + ggtitle("(d)")


####  

f <- grid.arrange(f3A, f4A, f3B, f4B,  ncol = 2)

ggsave("covid_2.png",plot = f, 
       width = 12, height = 8)