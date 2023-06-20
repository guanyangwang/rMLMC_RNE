
simulator_t_martingale <- function(trajectory, num_samples, location = pi/2){
  if(length(trajectory) == 0) samples <- rt(num_samples, df = 10, ncp = 0.5)
  else  samples <- rt(num_samples, df= 10, ncp = 0.5)
  return(samples)
} 



func1 <- function(traj, sample){
  return(sin(traj[length(traj)] + sample))
}

func2 <- function(traj, sample){
  return(sin(traj[length(traj)] - sample))
}

func3 <- function(traj){
  return(traj[length(traj)])
}

funs <- list(func1, func2, func3)