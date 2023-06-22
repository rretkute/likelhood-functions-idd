#' Leprocy ODE model
#'
#' @param t 
#' @param x 
#' @param parameters 
#'
#' @return
#' @export
#'
#' @examples

Leprocy_ODE_model <- function(t, x, parameters){ 
  S <- x[1]
  E <- x[2]
  I <- x[3]
  R <- x[4]
  D <- x[5]
  
  N <- S + E + I + R + D
  
  with( as.list(parameters), { 
    dS <- mu * N - beta * S * (I + D)/N - mu * S
    dE <- beta * S * (I + D)/N - sigma * E - mu * E
    dI <- (1 - pE) * sigma * E - gamma * I - theta * I - mu * I
    dR <- gamma * (I + D) + rho * D - mu * R
    dD <-  pE * sigma * E + theta * I - (rho + gamma) * D - mu * D
    
    res <- c(dS,dE, dI, dR, dD)
    list(res)
  }
  )
}
