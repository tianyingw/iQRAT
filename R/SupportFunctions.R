

QS_Common <- function(null_Sn_0_c, null_Sn_1_c, null_Sn_2_c, Xstar_c){
  #### Common set
  QSKAT_W = sum((null_Sn_0_c)^2)
  QSKAT_L = sum((null_Sn_1_c)^2)
  QSKAT_N = sum((null_Sn_2_c)^2)


  Q.test = c(QSKAT_W, QSKAT_L, QSKAT_N)
  CovS = crossprod(Xstar_c)
  pval_S_c = Liu.pval(CovS, Q.test)

  #### Cauchy test
  CT_S_c = mean(tan((0.5-pval_S_c[1:2])*pi))
  pval_S_CT2_c = pcauchy(CT_S_c, lower.tail = FALSE)
  CT_S_c = mean(tan((0.5-pval_S_c)*pi))
  pval_S_CT3_c = pcauchy(CT_S_c, lower.tail = FALSE)

  return(list(pval_S_CT3_c = pval_S_CT3_c, pval_S_CT2_c = pval_S_CT2_c))

}



QS_Rare <- function( null_Sn_0_r, null_Sn_2_r,  Xstar_r){
  #### Rare Set
  QSKAT_W = sum((null_Sn_0_r)^2)
  QSKAT_N = sum((null_Sn_2_r)^2)


  Q.test = c(QSKAT_W, QSKAT_N)
  CovS = crossprod(Xstar_r)
  pval_S_r = Liu.pval(CovS, Q.test)
  pval_S_CT2_r = pval_S_r[1]
  #### Cauchy test
  CT_S_r = mean(tan((0.5-pval_S_r)*pi))
  pval_S_CT3_r = pcauchy(CT_S_r, lower.tail = FALSE)



  return(list(pval_S_CT3_r = pval_S_CT3_r, pval_S_CT2_r =pval_S_CT2_r ))
}




QB_Common <- function(null_Sn_0_c, null_Sn_1_c, null_Sn_2_c,  Xstar_c, ind_common){
  ### Commmon
  QB_W = (sum(null_Sn_0_c))^2
  QB_L = (sum(null_Sn_1_c))^2
  QB_N = (sum(null_Sn_2_c))^2

  Q.test = c(QB_W, QB_L, QB_N)
  CovS = crossprod(Xstar_c%*%rep(1, length(ind_common)))
  pval_B_c = Liu.pval(CovS, Q.test)

  CT_B_c = mean(tan((0.5-pval_B_c[1:2])*pi))
  pval_B_CT2_c = pcauchy(CT_B_c, lower.tail = FALSE)
  CT_B_c = mean(tan((0.5-pval_B_c)*pi))
  pval_B_CT3_c = pcauchy(CT_B_c, lower.tail = FALSE)

  return(list(pval_B_CT3_c = pval_B_CT3_c, pval_B_CT2_c = pval_B_CT2_c))
}



QB_Rare <- function(null_Sn_0_r, null_Sn_2_r,  Xstar_r, ind_rare){

  QB_W = (sum(null_Sn_0_r))^2
  QB_N = (sum(null_Sn_2_r))^2


  Q.test = c(QB_W,QB_N)
  CovS = crossprod(Xstar_r%*%rep(1, length(ind_rare)))
  pval_B_r = Liu.pval(CovS, Q.test)
  pval_B_CT2_r = pval_B_r[1]
  ### Cauchy test
  CT_B_r = mean(tan((0.5-pval_B_r)*pi))
  pval_B_CT3_r = pcauchy(CT_B_r, lower.tail = FALSE)

  return(list(pval_B_CT3_r = pval_B_CT3_r, pval_B_CT2_r = pval_B_CT2_r))
}



Q_Combine <- function(out_c, out_r){
 out_c = as.numeric(out_c)
 out_r = as.numeric(out_r)

  ###### Combine p values


  Fisher_S = -2*log(out_c[2]) - 2*log(out_r[2])
  pval_S2 = pchisq(Fisher_S, df = 4, lower.tail = FALSE)

  Fisher_S = -2*log(out_c[1]) - 2*log(out_r[1])
  pval_S3 = pchisq(Fisher_S, df = 4, lower.tail = FALSE)

  return(list(pval3 = pval_S3, pval2 = pval_S2))
}

Liu.pval <- function(CovS, Q.all){
    A1<-CovS
    A2<-A1 %*% A1

    c1<-rep(0,4)
    c1[1]<-sum(diag(A1))
    c1[2]<-sum(diag(A2))
    c1[3]<-sum(A1*t(A2))
    c1[4]<-sum(A2*t(A2))
    param<-Get_Liu_Params_Mod(c1)

    Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
    Q.Norm1<-Q.Norm * param$sigmaX + param$muX
    p.value<- pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE)

    return(p.value)
}

Get_Liu_Params_Mod<-function(c1){
    ## Helper function for getting the parameters for the null approximation

    muQ<-c1[1]
    sigmaQ<-sqrt(2 *c1[2])
    s1 = c1[3] / c1[2]^(3/2)
    s2 = c1[4] / c1[2]^2

    beta1<-sqrt(8)*s1
    beta2<-12*s2
    type1<-0

    #print(c(s1^2,s2))
    if(s1^2 > s2){
        a = 1/(s1 - sqrt(s1^2 - s2))
        d = s1 *a^3 - a^2
        l = a^2 - 2*d
    } else {
        type1<-1
        l = 1/s2
        a = sqrt(l)
        d = 0
    }
    muX <-l+d
    sigmaX<-sqrt(2) *a

    re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
    return(re)
}



GenerateSample <- function(hp, temp_list, temp_list2, pnum, local, error_id, Standardize = TRUE){


    n = length(temp_list$x_index)/2
    x1 = hp[temp_list$x_index[1:n], temp_list$indexSelect:(temp_list$indexSelect + pnum -1)]
    x2 = hp[temp_list$x_index[(n+1):(2*n)], temp_list$indexSelect:(temp_list$indexSelect + pnum -1)]
    x = x1 + x2

    if(local == TRUE){
        u=seq(0,1,length=50)
        betax=function(u, u0) pmax( 0, 5*(u-u0)/(1-u0))
        uu = runif(n); slopex = betax(uu,u0=0.7)


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
        y_a=1+1.2*temp_list2$C+(x%*%temp_list$beta)*slopex+error
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


