rm(list = ls())
source('vanilla_MC2.R')
source('mlmc.R')
source('sin_t_setup.R')
# some auxiliary functions
cummean <- function(x){cumsum(x)/1:length(x)}
cumsd <- function(x){
  cmean_x <- cummean(x)
  cmean_xx <- cummean(x*x)
  cvar_x <- cummean(x*x) - cummean(x)^2
  csd_x <-sqrt(cvar_x)
  return(csd_x)
}
set.seed(1)

n_plt <- seq(1e3,2e6,4e4)

res_READ <- rep(0,length(n_plt))
for(i in 1:length(n_plt)){
  print(i)
  res_READ[i] <-  mean(unlist(replicate(n_plt[i], recursive_MLMC(simulator = simulator_t_martingale, trajectory = c(), 
                                                                 params = c(0.74,0.6), funs = funs, time_horizon = 3), simplify=FALSE)))
  print(res_READ[i])
}

# avg sample size for calling recursive_MLMC once
avg_sample <- 0.74/(2 * 0.74 - 1) * 0.6/(2* 0.6 - 1)

#READ_estimate <- cummean(res_READ[1,])[n_plt]



set.seed(1)
m_max = 400
m_rep = 1
m <- seq(10,m_max,3)
res_NMC1 <- rep(0, length(m))
for(j in 1:length(m)){
  print(j)
  res_NMC1[j] <- mean(NMC(simulator = simulator_t_martingale,
                          trajectory = c(), samples_vec = rep(m[j],3),funs, time_horizon = 3))
  print(res_NMC1[j])}



set.seed(1)
l_max = 100
l_rep = 1
l <- seq(4,l_max,1)
res_NMC2 <- rep(0, length(l))
for(j in 1:length(l)){
  print(j)
  res_NMC2[j] <- mean(NMC(simulator = simulator_t_martingale,
                          trajectory = c(), samples_vec = c(l[j]^2,l[j],l[j] ),funs, time_horizon = 3)) 
  print(res_NMC2[j])
}


library(ggplot2)
method = c(rep('READ', length(n_plt)), rep('NMC1', length(m)), rep('NMC2', length(l)))

df <- data.frame(Log_Samples = c(log10(n_plt* avg_sample), log10(m^3), log10(l^4)), 
                 Estimation = c(res_READ, res_NMC1, res_NMC2), method = method)

ground <- mean(res_READ)
p <- ggplot(data = df, aes(x = Log_Samples, y = Estimation, color = method)) + geom_point(aes(shape = method),size = 2, alpha = 0.8)+ xlab('log(Total Sample Cost)') + ylab('Estimation') + xlim(4, 7.5) + ylim(0.2,0.35)
p <- p + geom_hline(yintercept = 0.261, color = 'black', linetype = 'dashed', lwd = 0.3)
p

set.seed(1)
READ_estimate <- unlist(replicate(1e6, recursive_MLMC(simulator = simulator_t_martingale, trajectory = c(), 
                                                      params = c(0.74,0.6), funs = funs, time_horizon = 3), simplify=FALSE))

y = READ_estimate
estimate_y <- cummean(y)
sd_y <-cumsd(y)
points <- seq(4e4,1e6,1e3)
l = length(points)

df_y <- data.frame(estimate = estimate_y[points], sd = sd_y[points], no_estimators = points)
df_y$upper <- df_y$estimate + 1.96 * df_y$sd/sqrt(df_y$no_estimators)
df_y$lower <- df_y$estimate - 1.96 * df_y$sd/sqrt(df_y$no_estimators)

p2 <- ggplot(data = df_y, aes(x = no_estimators)) + geom_line(aes(y = estimate), color = 'steelblue4') 

p2 <- p2 + geom_line(aes(y = upper),linetype = "dotted", color = 'steelblue4') + 
  geom_line(aes(y = lower),linetype = "dotted", color = 'steelblue4') 
p2 <- p2 + xlab('Number of estimators') + ylab('Running Average') + ylim(0.2,0.3) 
p2

