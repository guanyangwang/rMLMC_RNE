##  reproduce the experiments and plots in Section 3
rm(list = ls())
source('vanilla_NMC.R')
source('mlmc.R')
source('sin_gaussian_setup.R')

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
# READ estimator, each repetition generates 10^6 estimators, 20 repetitions
# parameter r_1 = 0.74, r_2 = 0.6
set.seed(1)
n_max = 1e6
n_rep = 20
res_READ <- matrix(0,n_rep,n_max)
for(i in 1:n_rep){
  print(i)
  res_READ[i,] <-  unlist(replicate(n_max, recursive_MLMC(simulator = simulator_gaussian_martingale, trajectory = c(), 
                                                         params = c(0.74,0.6), funs = funs, time_horizon = 3), simplify=FALSE))
  }

# average sample size for each call of recursive_MLMC
avg_sample <- 0.74/(2 * 0.74 - 1) * 0.6/(2* 0.6 - 1)


temp <- t(apply(res_READ, 1, cummean))
MSE_READ <- (colMeans(temp) - exp(-0.5))^2

# points we choose (to make plot)
n_plt <- floor(10^seq(0.5,6,0.2)) 
#plot(log10(n_plt* avg_sample), log10(MSE_READ[n_plt]))

# the corresponding log(MSE)
Log_READ = log10(MSE_READ[n_plt])

# fit a regression between log(MSE) and log(total sample cost)
# slope corresponds to (empirical) convergence rate
lm(log10(MSE_READ[n_plt])~(log10(n_plt* avg_sample)))


#########################
### NMC1 in Rainforth et al. (2018)
#########################

set.seed(1)

#NMC_1 uses N0=N1=N2 = m (our parameter here)
m_max = 400
m_rep = 20
m <- seq(10,m_max,15)
# total sample cost is m^3
res_NMC1 <- matrix(0, m_rep, length(m))

# for each m, repeat 20 times and calculate the MSE
for(i in 1:m_rep){
  print(i)
  for(j in 1:length(m)){
  print(j)
  res_NMC1[i,j] <- mean(NMC(simulator = simulator_gaussian_martingale,
                            trajectory = c(), samples_vec = rep(m[j],3),funs, time_horizon = 3)) 
  res_NMC1[i,j] <- (res_NMC1[i,j] - exp(-0.5))^2}
}
#plot(log10(m^3), log10(colMeans(res_NMC1)))

# fit linear regression to get the convergence rate
lm(log10(colMeans(res_NMC1))~log10(m^3))

#log(MSE) of NMC1
Log_NMC1 <- log10(colMeans(res_NMC1))



#########################
### NMC2 in Rainforth et al. (2018)
#########################
set.seed(1)

#NMC2 uses N0 = N1^2 = N2^2 = l

# total sample cost is N0N1N2 = l^4
l_max = 100
l_rep = 20
l <- seq(4,l_max,4)
res_NMC2 <- matrix(0, m_rep, length(l))
for(i in 1:l_rep){
  print(i)
  for(j in 1:length(l)){
    print(j)
    res_NMC2[i,j] <- mean(NMC(simulator = simulator_gaussian_martingale,
                              trajectory = c(), samples_vec = c(l[j]^2,l[j],l[j] ),funs, time_horizon = 3)) 
    res_NMC2[i,j] <- (res_NMC2[i,j] - exp(-0.5))^2}
}

#plot(log10(l^4), log10(colMeans(res_NMC2)))

Log_NMC2 <- log10(colMeans(res_NMC2))
lm(log10(colMeans(res_NMC2))~log10(l^4))


#####################
### count (clock) time
#####################
set.seed(1)

# READ, 10^5 repetitions

t_start_READ <- Sys.time()

r_READ<- mean(unlist(replicate(1e5, recursive_MLMC(simulator = simulator_gaussian_martingale, trajectory = c(), 
                                                      params = c(0.74,0.6), funs = funs, time_horizon = 3), simplify=FALSE)))
# 11.86271 secs
time_READ <- Sys.time() - t_start_READ

# err = 3.186* 10^{-6}
sq_error_READ <- (r_READ - exp(-0.5))^2


# NMC1, N_0 = 400, total sample cost = N_0^3 = 6.4 * 10^7
set.seed(1)
t_start_NMC1 <- Sys.time()
r_NMC1 <- mean(NMC(simulator = simulator_gaussian_martingale,
                   trajectory = c(), samples_vec = c(400,400,400),funs, time_horizon = 3))
# time = 43.792s
time_NMC1 <- Sys.time() - t_start_NMC1

# error = 3.51 * 10^{-4}
sq_error_NMC1 <- (r_NMC1 - exp(-0.5))^2

time_NMC1 * sq_error_NMC1
set.seed(1)

t_start_NMC2 <- Sys.time()

r_NMC2 <- mean(NMC(simulator = simulator_gaussian_martingale,
                        trajectory = c(), samples_vec = c(1e4,100,100),funs, time_horizon = 3))
# time = 1.21 mins = 72.6 seconds
time_NMC2 <- Sys.time() - t_start_NMC2

# error = 6.80 * 10^{-5}

sq_error_NMC2 <- (r_NMC2 - exp(-0.5))^2
###############################
# plots here
###############################


###############################
# The following reproduces 
# Figure 1 in Section 3
###############################
library(ggplot2)
# Figure 1: plot the log(MSE) vs log(Total Sample Size) for all three methods
method = c(rep('READ', length(n_plt)), rep('NMC1', length(m)), rep('NMC2', length(l)))

df <- data.frame(Log_Samples = c(log10(n_plt* avg_sample), log10(m^3), log10(l^4)), 
                 Log_MSE = c(Log_READ, Log_NMC1, Log_NMC2), method = method)
p <- ggplot(data = df, aes(x = Log_Samples, y = Log_MSE, color = method)) + geom_point(aes(shape = method),size = 0.8) + 
      stat_smooth(method = lm, se = FALSE, size = 0.3, fullrange = TRUE) + xlab('log(Total Sample Cost)') + ylab('log(MSE)')
p

###############################
# The following reproduces 
# Figure 2 in Section 3
###############################
# Figure 2: plot the running everage with confidence interval
y = res_READ[1,]
estimate_y <- cummean(y)
sd_y <-cumsd(y)
points <- seq(4e4,1e6,1e3)
l = length(points)

df_y <- data.frame(estimate = estimate_y[points], sd = sd_y[points], no_estimators = points)
df_y$upper <- df_y$estimate + 1.96 * df_y$sd/sqrt(df_y$no_estimators)
df_y$lower <- df_y$estimate - 1.96 * df_y$sd/sqrt(df_y$no_estimators)

p2 <- ggplot(data = df_y, aes(x = no_estimators)) + geom_line(aes(y = estimate), color = 'steelblue4') + ylim(0.57,0.65)

p2 <- p2 + geom_line(aes(y = upper),linetype = "dotted", color = 'steelblue4') + 
    geom_line(aes(y = lower),linetype = "dotted", color = 'steelblue4') 

p2 <- p2 + geom_hline(yintercept = exp(-1/2), color = 'red', linetype = 'dashed')
p2 <- p2 + xlab('Number of estimators') + ylab('Running Average')
p2


