vanilla_MC2 <- function(simulator, trajectory = c(), samples_vec,funs, time_horizon){
  #' vanilla Monte Carlo estimator by path simulation
  #' @param simulator: function, simulatior to generate new samples, 
  #' when trajectory = empty, simulator generates y^{(0)}, otherwise,
  #' it generates y^{(d)} given the previous trajectory (y^{(0)}, ... , y^{(d-1)})
  #' @param trajectory: vector, the trajectory of process history, default = empty
  #' @param num_paths: integer, the number of the paths sampled
  #' @return scalar,which is a biased estimator of the RNE
  #' @noRd
  # @export
  
  func = funs[[1]]
  samples <- samples_vec[1]
  x <- simulator(trajectory, num_samples = 1)
  trajectory = c(trajectory, x)
  if(isTRUE(time_horizon == 1)) return(func(trajectory))
  samples_vec <- samples_vec[-1]
  funs <- funs[-1]
  samples <- mean(unlist(replicate(samples, vanilla_MC2(simulator, trajectory, samples_vec, funs,time_horizon - 1), 
                              simplify=FALSE)))
  return(func(trajectory, samples))
  }


## the following version is a faster NMC estimator, but only works when 
## 1. every function depends on only the last element of trajectory, 2. D = 2
## this estimator is used in our paper to compare with our READ estimator
NMC <- function(simulator, trajectory = c(), samples_vec,funs, time_horizon = 3){
  #' vanilla Monte Carlo estimator by path simulation
  #' @param function, simulatior to generate new samples, 
  #' when trajectory = empty, simulator generates y^{(0)}, otherwise,
  #' it generates y^{(d)} given the previous trajectory (y^{(0)}, ... , y^{(d-1)})
  #' @param trajectory: vector, the trajectory of process history, default = empty
  #' @param num_paths: integer, the number of the paths sampled
  #' @return scalar,which is a biased estimator of the RNE
  #' @noRd
  # @export
  
  x_0 <- simulator(trajectory, samples_vec[1])
  x_1 <- lapply(x_0, function(x) simulator(trajectory = c(x), samples_vec[2]))
  x_1 <- matrix(unlist(x_1), ncol = samples_vec[2], byrow = TRUE)
  func <- funs[[3]]
  y <- function(x) mean(vapply(x, func, numeric(1)))
  res_1 <- structure(vapply(x_1, function(x) y(simulator(x, samples_vec[3])), numeric(1)), dim = dim(x_1))
  res_0 <- matrix(0, dim(res_1)[1], dim(res_1)[2])
  for(i in 1: dim(x_1)[1]){
    for(j in 1: dim(x_1)[2]){
      res_0[i,j] = funs[[2]](x_1[i,j], res_1[i,j])
    }
  }
  gamma_1 = rowMeans(res_0)
  res = rep(0,length(x_0))
  for(i in 1:length(x_0)){
    res[i] <- funs[[1]](x_0[i], gamma_1[i])
  }
  return(mean(res))
}

