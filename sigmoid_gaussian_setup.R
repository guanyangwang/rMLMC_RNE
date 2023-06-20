# The function generates random samples from a normal distribution
# It takes in two arguments: trajectory and num_samples
# If the length of trajectory is equal to zero, it generates num_samples of random samples from a normal distribution with mean 0 and sd 1
# If trajectory is not empty, the samples are generated with mean equal to the last element of the trajectory vector and sd of 1
# The generated samples are then returned by the function
simulator_gaussian_martingale <- function(trajectory, num_samples){
  if(length(trajectory) == 0) samples <- rnorm(num_samples)
  else  samples <- rnorm(num_samples, mean = trajectory[length(trajectory)], sd = 1)
  return(samples)
} 


#### the following three functions are g_0, g_1, g_2 in Appendix F ###
#### g_0 (y_0, z) = sigmoid(y_0 + z) #################################
#### g_1 (y_0, y_1, z) = sigmoid(y_1 + z) ############################
#### g_2 (y_0, y_1, y_2) = sigmoid(y_2)   ############################



# The function calculates the sigmoid function of a given input x
# The function takes in a single argument x and returns the sigmoid of x, which is calculated as 1/(1 + e^-x)
sigmoid <- function(x){
  return(1/(1 + exp(-x)))
}

# The function takes in two arguments: traj and sample
# traj is a vector of previous values and sample is a scalar value
# The function returns the sigmoid of the last element of the traj vector plus the sample argument
func1 <- function(traj, sample){
  return(sigmoid(traj[length(traj)] + sample))
}

# The function takes in two arguments: traj and sample
# traj is a vector of previous values and sample is a scalar value
# The function returns the sigmoid of the last element of the traj vector plus the sample argument
func2 <- function(traj, sample){
  return(sigmoid(traj[length(traj)] + sample))
}

# The function takes in one argument: traj
# traj is a vector of previous values
# The function returns the sigmoid of the last element of the traj vector
func3 <- function(traj){
  return(sigmoid(traj[length(traj)]))
}

# The last line of code creates a list named 'funs' containing the three functions defined above
# This list can be used to store and call the functions
funs <- list(func1, func2, func3)
