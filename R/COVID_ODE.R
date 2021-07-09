#' COVID-19 ODE model
#'
#' @param t 
#' @param x 
#' @param parameters 
#'
#' @return list of derivatives
#' @export
#'
#
COVID_ODE <- function(t, x, parameters){ 
  
  S <- x[1] # Susceptible
  I_1 <- x[2] # Mild infection
  I_2 <- x[3] # Severe infection
  R <- x[4] # Recovered
  D <- x[5] # Deaths
  
  with( as.list(parameters), { 
    dS <- - (beta * S * (I_1 + I_2)) / N
    dI_1 <- (1 - p) * (beta * S * (I_1 + I_2)) / N - gamma_1 * I_1
    dI_2 <- p* (beta * S * (I_1 + I_2)) / N - (gamma_2 + alpha) * I_2
    dR <- gamma_1 * I_1 + gamma_2 * I_2
    dD <-  alpha * I_2
    
    res <- c(dS, dI_1, dI_2, dR, dD)
    list(res)
  }
  )
}
