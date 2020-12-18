
QS_Single <- function(null_Sn, Xstar, method.p = "davies", acc =1e-9){
  #### Common set
  Q.test = sum((null_Sn)^2)
  CovS = crossprod(Xstar)
  if(method.p == "liu"){
    pval_S_c = Liu.pval(CovS, Q.test)
  }
  if(method.p == "davies"){
    pval_S_c = davies.pval(CovS, Q.test, acc =acc)
  }


  return( pval_S_c)

}



QB_Single <- function(null_Sn,  Xstar, ind1, method.p = "davies", acc =1e-9){
  ### Commmon
  Q.test = (sum(null_Sn ))^2

  if(length(ind1) ==1){
    Xstar = matrix(Xstar, nrow=length(Xstar), ncol = 1)
  }
  CovS = crossprod(Xstar%*%rep(1, length(ind1)))

  if(method.p == "liu"){
    pval_B = Liu.pval(CovS, Q.test)
  }
  if(method.p == "davies"){
    pval_B = davies.pval(CovS, Q.test)
  }

  return(pval_B)
}



