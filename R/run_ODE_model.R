#' Run ODE model
#'
#' @param parameters 
#' @param initial_state 
#' @param tmax 
#' @param dt Time aggregation interval (days)
#' @param model 
#'
#' @return solution to ODE model
#' @export
#'
#' @examples
run_ODE_model <- function(parameters, initial_state, tmax, dt, model){
  times <- seq(from = 0,to = tmax, by = dt) 
  output_det<- as.data.frame(ode(initial_state, times, model, parameters, method=rk4))
  return(output_det)
}
