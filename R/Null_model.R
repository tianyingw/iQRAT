#' Fit the null model
#' This function provides null model fit. The default option is to use quantreg since it is faster. When
#' dataset is too large, quantreg will have memory limitation issue, and we will use the self-implemented
#' algorithm.
#' @param Y A vector of response, size \eqn{n}.
#' @param C A vector or matrix of covariates.
#' @import quantreg Matrix MASS
#' @export
#'
Null_model <- function(Y,C){
  

   null.fit =rq(Y~C, tau = -1)
   
  
  
  return(null.fit)
}
