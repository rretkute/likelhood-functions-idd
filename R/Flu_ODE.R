#' Flu ODE 
#'
#' @param t 
#' @param x 
#' @param parameters 
#'
#' @return
#' @export
#'
#' @examples
Flu_ODE <- function(t, x, parameters){
  S <- x[1]
  I <- x[2]
  R <- x[3]
  N <- S + I + R

  with( as.list(parameters), {
    dS <- -beta * (S * I) / N
    dI <- beta * (S * I) / N  - gamma * I
    dR <- gamma * I
    res <- c(dS, dI, dR)
    list(res)
  }
  )
}