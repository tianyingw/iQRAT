#' Function for meta-analysis using Fisher's method
#'
#' This funciton uses Fisher's method to combine a list of independent p values and return the p value for mets-analysis.
#'
#' @param pval.list A list of independent p values. Each element is a vector with same length.
#'
#' @return The p value for meta-analysis.
#' @export
#'
meta_fisher <- function(pval.list){
  nstudy = length(pval.list)
  if(is.vector((pval.list[[1]]))){
    mattemp = rep(0, length(pval.list[[1]]))
    matind = 1
  }else{
    ngene = nrow(pval.list[[1]])
    nmethod  = ncol(pval.list[[1]])
    mattemp = matrix(0, nrow = ngene, ncol = nmethod)
    matind = 2
  }


  for(i in 1:nstudy){

    A = as.numeric(pval.list[[i]])
    A = (A == 0)*1e-15 + (A>0)*A

    if(matind == 2){
      A = matrix(A, ncol = nmethod, nrow = ngene)
    }

    mattemp = mattemp - 2*log(A)
  }
  pval_meta = pchisq( mattemp, df = 2*nstudy, lower.tail = FALSE)


  return(pval_meta)
}
  return(pval_meta)
}


