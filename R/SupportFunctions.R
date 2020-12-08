
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




#' Generate a simulation sample
#' @export

GenerateSample <- function(hp, temp_list, temp_list2, pnum, local, error_id, Standardize = TRUE, local.tau = NULL){


    n = length(temp_list$x_index)/2
    x1 = hp[temp_list$x_index[1:n], temp_list$indexSelect:(temp_list$indexSelect + pnum -1)]
    x2 = hp[temp_list$x_index[(n+1):(2*n)], temp_list$indexSelect:(temp_list$indexSelect + pnum -1)]
    x = x1 + x2

    if(local == TRUE){
        u=seq(0,1,length=50)
        if(local.tau < 0.5){
          betax=function(u, u0) pmin( 0, 5*(u-u0)/(u0))
        }else{
          betax=function(u, u0) pmax( 0, 5*(u-u0)/(1-u0))
        }

        uu = runif(n)
        slopex = betax(uu,u0=local.tau)


        if(error_id == 1){
            error = qnorm(uu)
        }
        if(error_id == 2){
            error = qchisq(uu, df = 2)
        }
        if(error_id == 3){
            error = qcauchy(uu)
        }
        if(error_id == 4){
            error = qt(uu, 2)
        }
        if(error_id == 5){
            error = qt(uu, 3)
        }
        y_a=1+1.2*temp_list2$C+x%*%temp_list$beta*slopex+error
        y_0=1+1.2*temp_list2$C+error


    }else{

        y_a=1+1.2*temp_list2$C+x%*%temp_list$beta+(1+x%*%temp_list$gamma)*temp_list2$error
        y_0=1+1.2*temp_list2$C+temp_list2$error}

    if(Standardize == TRUE){
        percentile<-ecdf(y_0)
        y_0<-qnorm(percentile(y_0)-min(percentile(y_0))/2)

        percentile<-ecdf(y_a)
        y_a<-qnorm(percentile(y_a)-min(percentile(y_a))/2)}

    return(Data = list(y_0 = y_0, y_a = as.vector(y_a), x = x, c = temp_list2$C, error = temp_list2$error))

}




