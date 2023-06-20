
# Function to simulate a Gaussian martingale
simulator_gaussian_martingale <- function(trajectory, num_samples, mean = pi/2){
  # Check if trajectory is empty
  if(length(trajectory) == 0) {
    # Generate samples from normal distribution with mean = mean and sd = 1
    samples <- rnorm(num_samples, mean = mean, sd = 1)
  } else {
    # Generate samples from normal distribution with mean = last element of trajectory and sd = 1
    samples <- rnorm(num_samples, mean = trajectory[length(trajectory)], sd = 1)
  }
  # Return the generated samples
  return(samples)
} 


#### the following three functions are g_0, g_1, g_2 in Section 3 ####
#### g_0 (y_0, z) = sin(y_0 + z) #####################################
#### g_1 (y_0, y_1, z) = sin(y_1 - z) ################################
#### g_2 (y_0, y_1, y_2) = y_2 #######################################

# Function 1: takes in traj and sample, returns sin(traj[last element] + sample)
func1 <- function(traj, sample){
  return(sin(traj[length(traj)] + sample))
}

# Function 2: takes in traj and sample, returns sin(traj[last element] - sample)
func2 <- function(traj, sample){
  return(sin(traj[length(traj)] - sample))
}

# Function 3: takes in traj, returns the last element
func3 <- function(traj){
  return(traj[length(traj)])
}

# A list containing all three functions
funs <- list(func1, func2, func3)
