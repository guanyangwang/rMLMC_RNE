
recursive_MLMC <- function(simulator, trajectory = c(), params,funs, time_horizon){
  #' unbiased estimator of estimating the repeatedly nested expectations
  #' @param simulator: function, simulatior to generate new samples, 
  #' when trajectory = empty, simulator generates y^{(0)}, otherwise,
  #' it generates y^{(d)} given the previous trajectory (y^{(0)}, ... , y^{(d-1)})
  #' @param trajectory: vector, the trajectory of process history, default = empty
  #' @param params: vector, the parameters of the geometric distribution at each level
  #' @param funs: list, the functions of interest at each level
  #' @param time_horizon: integer, equals D + 1 in our READ paper
  #' @return scalar, an unbiased estimator of gamma_{D - time_horizon}
  #' @noRd
  #' Target: We want to get ONE unbiased estimator of gamma_0 
  #' How to get the estimator: 
  #'   When time_horizon = 1, we generate one sample from the simulator conditioning on the trajectory.
  #'   When time_horizon > 1, we call the algorithm itself 2^N times with one less time_horizon,
  #'   and apply the rMLMC method
  # @export
  
 # get the first function 
  func = funs[[1]]
 # append x to the trajectory
  x <- simulator(trajectory, num_samples = 1)
  curr = x
  trajectory = c(trajectory, x)
  
  # if already reach the last depth (i.e., horizon = 1)
  if(isTRUE(time_horizon == 1)) {
    res <- func(trajectory)
    return(res)
  }
  N <- rgeom(n = 1, prob = params[1]) # number of samples used -- warning: in R, rgeom starts with 0
  num_samples <- 2^N
  pN <- dgeom(N, prob = params[1]) # the pmf of N
  params = params[-1] # remove the first parameter from the list
  funs = funs[-1] # remvoe the first function from the list
  
  # the recursive with new list of functions, new parameters, and one less time horizon
  samples <- unlist(replicate(num_samples, recursive_MLMC(simulator, trajectory, params, funs, time_horizon - 1), 
                              simplify=FALSE))
  if (isTRUE(N == 0)) return(func(trajectory, samples)/pN) # edge case where N = 0
  else{
  # split samples into even and odd terms
    samples_odd <- samples[seq(1,num_samples,2)]
    samples_even <- samples[-seq(1,num_samples,2)]
    mean_odd <- mean(samples_odd)
    mean_even <- mean(samples_even)
    mean_all <- mean_odd * 0.5 + mean_even * 0.5
    estimator_all <- func(trajectory, mean_all)
    estimator_odd <- func(trajectory, mean_odd)
    estimator_even <- func(trajectory, mean_even)
  # return the ratio estimator
    Delta = estimator_all - (estimator_odd + estimator_even)/2
    return(Delta/pN)
  }
}





