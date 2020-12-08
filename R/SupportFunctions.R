
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



##################################################
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



###############################
# The following two functions related to davies method are from SKAT package
# Reference: @Manual{SKATr,
# title = {SKAT: SNP-Set (Sequence) Kernel Association Test},
# author = {Seunggeun Lee and with contributions from Larisa Miropolsky and Michael Wu},
# year = {2017},
# note = {R package version 1.3.2.1},
# url = {https://CRAN.R-project.org/package=SKAT},
# }
###############################
davies.pval <- function(CovS, QSKAT, acc =1e-9){
  eig1 = eigen(CovS)$values
  acc1 = acc
  temp_davies = Get_PValue.Lambda(eig1, QSKAT)
  return(pval_davies = temp_davies$p.value)
}
davies.pval2 <- function(CovS, QSKAT, acc =1e-9){
  eig1 = eigen(CovS)$values
  eta = eig1[which(eig1 > 1e-05)]
  eta1 = eta[which(eta > mean(eta)/1e5)]
  acc1 = acc
  temp_davies = SKAT_davies(QSKAT, eta1, acc = acc1)
  while(temp_davies$ifault == 1){
    acc1 = acc1*10
    warning(paste0("use accuracy ", acc1))
    temp_davies = SKAT_davies(QSKAT, eta1, acc = acc1)
  }
  
  pval_davies = max(temp_davies$Qq,0)
  return(pval_davies)
}


SKAT_davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=1e6,acc=1e-9) {
  
  r <- length(lambda)
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")
  
  out <- .C("qfc",lambdas=as.double(lambda),noncentral=as.double(delta),df=as.integer(h),r=as.integer(r),sigma=as.double(sigma),q=as.double(q),lim=as.integer(lim),acc=as.double(acc),trace=as.double(rep(0,7)),ifault=as.integer(0),res=as.double(0),PACKAGE="SKAT")
  
  out$res <- 1 - out$res
  
  return(list(trace=out$trace,ifault=out$ifault,Qq=out$res))
  
}

