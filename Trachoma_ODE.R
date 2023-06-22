#' Deterministic Trachoma model
#'
#' @param t 
#' @param x 
#' @param parameters 
#'
#' @return
#' @export
#'
#' @examples

Trachoma_ODE_model <- function(t, x, parameters){ 
  I <- x[1]
  with( as.list(parameters), { 
    dI <- beta * (N-I)*I/N - gamma*I
    res <- c(dI)
    list(res)
  }
  )
}
