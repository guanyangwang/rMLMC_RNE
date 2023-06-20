##  reproduce the experiments and plots in Appendix F
rm(list = ls())
source('vanilla_NMC.R')
source('mlmc.R')
source('sigmoid_gaussian_setup.R')
library(hrbrthemes)
library(viridis)
# some auxiliary functions
cummean <- function(x){cumsum(x)/1:length(x)}
cumsd <- function(x){
  cmean_x <- cummean(x)
  cmean_xx <- cummean(x*x)
  cvar_x <- cummean(x*x) - cummean(x)^2
  csd_x <-sqrt(cvar_x)
  return(csd_x)
}

##################
### READ estimator
###################
# READ estimator, generates 10^6 estimators
# parameter r_1 = 0.74, r_2 = 0.6

set.seed(1)

# points we choose to draw plots
n_plt <- seq(1e3,2e6,2e4)

res_READ <- rep(0,length(n_plt))
for(i in 1:length(n_plt)){
  print(i)
  res_READ[i] <-  mean(unlist(replicate(n_plt[i], recursive_MLMC(simulator = simulator_gaussian_martingale, trajectory = c(), 
                                                          params = c(0.74,0.6), funs = funs, time_horizon = 3), simplify=FALSE)))
}

# avg sample size for calling recursive_MLMC once
avg_sample <- 0.74/(2 * 0.74 - 1) * 0.6/(2* 0.6 - 1)




#########################
### NMC1 in Rainforth et al. (2018)
#########################
set.seed(1)
m_max = 400
m_rep = 1
#NMC_1 uses N0=N1=N2 = m (our parameter here)
# total sample cost is m^3
m <- seq(10,m_max,3)

res_NMC1 <- rep(0, length(m))

# for each m[j] we run the NMC
for(j in 1:length(m)){
    print(j)
    res_NMC1[j] <- mean(NMC(simulator = simulator_gaussian_martingale,
                              trajectory = c(), samples_vec = rep(m[j],3),funs, time_horizon = 3))}


#########################
### NMC2 in Rainforth et al. (2018)
#########################
set.seed(1)
l_max = 100
l_rep = 1
#NMC2 uses N0 = N1^2 = N2^2 = l
# total sample cost is N0N1N2 = l^4

l <- seq(4,l_max,1)
res_NMC2 <- rep(0, length(l))
for(j in 1:length(l)){
    print(j)
    res_NMC2[j] <- mean(NMC(simulator = simulator_gaussian_martingale,
                              trajectory = c(), samples_vec = c(l[j]^2,l[j],l[j] ),funs, time_horizon = 3)) 
}


###############################
# plots here
###############################


##############
# Estimation plot (Figure 3)
##############

# the following reproduces Figure 3 in  Appendix F
library(ggplot2)
method = c(rep('READ', length(n_plt)), rep('NMC1', length(m)), rep('NMC2', length(l)))

df <- data.frame(Log_Samples = c(log10(n_plt* avg_sample), log10(m^3), log10(l^4)), 
                 Estimation = c(res_READ, res_NMC1, res_NMC2), method = method)

ground <- mean(res_READ)
p <- ggplot(data = df, aes(x = Log_Samples, y = Estimation, color = method)) + geom_point(aes(shape = method),size = 2, alpha = 0.8)+ xlab('log(Total Sample Cost)') + ylab('Estimation') + xlim(4, 7.5) + ylim(0.58,0.65)
p <- p + geom_hline(yintercept = 0.612, color = 'black', linetype = 'dashed', lwd = 0.3)
p



##############
# Heatmap (Figure 4-5)
##############


## calculate the standard deviation and normalized sd
## when r_0 and r_1 (p_1 and p_2 here) are within the range
## of theorem 2.2 
set.seed(1)
p_1 <- seq(0.6,0.74,0.014)
p_2 <- seq(0.55,0.60,0.005)
res <- matrix(0, length(p_1),length(p_2))
for(i in 1:length(p_1)){
  for(j in 1:length(p_2)){
    print(p_1[i])
    print(p_2[j])
    res[i,j] <- sd(unlist(replicate(1e6, recursive_MLMC(simulator = simulator_gaussian_martingale, trajectory = c(), 
                                                          params = c(p_1[i],p_2[j]), funs = funs, time_horizon = 3), simplify=FALSE)))
  }
}
res_normalize <- matrix(0,length(p_1),length(p_2))

for(i in 1:length(p_1)){
  for(j in 1:length(p_2)){
    res_normalize[i,j] <- res[i,j] * sqrt(p_1[i] * p_2[j]/((2 * p_1[i] - 1) * (2 * p_2[j] - 1)))
}}


## calculate the standard deviation and normalized sd
## when r_0 and r_1 (p_1_new and p_2_new here) are out of the range
## of theorem 2.2 

set.seed(1)
p_1_new <- seq(0.8,0.9,0.01)
p_2_new <- seq(0.7,0.8,0.01)
res_new <- matrix(0, length(p_1_new),length(p_2_new))
res_new_normalize <- matrix(0, length(p_1_new),length(p_2_new))
for(i in 1:length(p_1)){
  for(j in 1:length(p_2)){
    print(p_1_new[i])
    print(p_2_new[j])
    res_new[i,j] <- sd(unlist(replicate(1e6, recursive_MLMC(simulator = simulator_gaussian_martingale, trajectory = c(), 
                                                        params = c(p_1_new[i],p_2_new[j]), funs = funs, time_horizon = 3), simplify=FALSE)))
  }
}
for(i in 1:length(p_1)){
  for(j in 1:length(p_2)){
    res_new_normalize[i,j] <- res_new[i,j] * sqrt(p_1_new[i] * p_2_new[j]/((2 * p_1_new[i] - 1) * (2 * p_2_new[j] - 1)))
  }}


# the following reproduces Figure 4-5 in  Appendix F


data <- expand.grid(r0=p_1, r1=p_2)
data$sd <- matrix(res,ncol = 1, byrow = TRUE)
data$normalize_sd <- matrix(res_normalize,ncol = 1, byrow = TRUE)
ggplot(data, aes(r0, r1, fill= sd)) + 
  geom_tile() +
  scale_fill_gradient(low="azure2", high="steelblue4") +
  theme_ipsum()


ggplot(data, aes(r0, r1, fill= normalize_sd)) + 
  geom_tile() +
  scale_fill_gradient(low="azure2", high="steelblue4") +
  theme_ipsum()

data_2<- expand.grid(r0=p_1_new, r1=p_2_new)
data_2$sd <- matrix(res_new,ncol = 1, byrow = TRUE)
data_2$normalize_sd <- matrix(res_new_normalize, ncol = 1, byrow = TRUE)
ggplot(data_2, aes(r0, r1, fill= sd)) + 
  geom_tile() +
  scale_fill_gradient(low="azure2", high="steelblue4") +
  theme_ipsum()
ggplot(data_2, aes(r0, r1, fill= normalize_sd)) + 
  geom_tile() +
  scale_fill_gradient(low="azure2", high="steelblue4") +
  theme_ipsum()


