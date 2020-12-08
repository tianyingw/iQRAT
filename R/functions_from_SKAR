###############################
# The following functions are from SKAT package
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



SKAT_davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=1e6,acc=1e-9) {

  r <- length(lambda)
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")

  out <- .C("qfc",lambdas=as.double(lambda),noncentral=as.double(delta),df=as.integer(h),r=as.integer(r),sigma=as.double(sigma),q=as.double(q),lim=as.integer(lim),acc=as.double(acc),trace=as.double(rep(0,7)),ifault=as.integer(0),res=as.double(0),PACKAGE="SKAT")

  out$res <- 1 - out$res

  return(list(trace=out$trace,ifault=out$ifault,Qq=out$res))

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
Get_PValue.Lambda<-function(lambda,Q){

  #print(lambda)
  n1<-length(Q)

  p.val<-rep(0,n1)
  p.val.liu<-rep(0,n1)
  is_converge<-rep(0,n1)
  p.val.liu<-Get_Liu_PVal.MOD.Lambda(Q, lambda)

  for(i in 1:n1){
    out<-SKAT_davies(Q[i],lambda,acc=10^(-6))

    p.val[i]<-out$Qq
    #p.val.liu[i]<-SKAT_liu(Q[i],lambda)

    is_converge[i]<-1

    # check convergence
    if(length(lambda) == 1){
      p.val[i]<-p.val.liu[i]
    } else if(out$ifault != 0){
      is_converge[i]<-0
    }

    # check p-value
    if(p.val[i] > 1 || p.val[i] <= 0 ){
      is_converge[i]<-0
      p.val[i]<-p.val.liu[i]
    }
  }

  p.val.msg = NULL
  p.val.log=NULL
  #cat(p.val[1])
  if(p.val[1] == 0){

    param<-Get_Liu_Params_Mod_Lambda(lambda)
    p.val.msg<-Get_Liu_PVal.MOD.Lambda.Zero(Q[1], param$muQ, param$muX, param$sigmaQ, param$sigmaX, param$l, param$d)
    p.val.log<-Get_Liu_PVal.MOD.Lambda(Q[1], lambda, log.p=TRUE)[1]

  }

  return(list(p.value=p.val, p.val.liu=p.val.liu, is_converge=is_converge, p.val.log=p.val.log, pval.zero.msg=p.val.msg))

}

Get_Liu_PVal.MOD.Lambda<-function(Q.all, lambda, log.p=FALSE){

  param<-Get_Liu_Params_Mod_Lambda(lambda)

  Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
  Q.Norm1<-Q.Norm * param$sigmaX + param$muX
  p.value<- pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE, log.p=log.p)

  return(p.value)

}


Get_Liu_PVal.MOD.Lambda.Zero<-function(Q, muQ, muX, sigmaQ, sigmaX, l, d){


  Q.Norm<-(Q - muQ)/sigmaQ
  Q.Norm1<-Q.Norm * sigmaX + muX

  temp<-c(0.05,10^-10, 10^-20,10^-30,10^-40,10^-50, 10^-60, 10^-70, 10^-80, 10^-90, 10^-100)
  #qchisq(temp, df=1000000000,lower.tail=FALSE)
  out<-qchisq(temp,df = l,ncp=d, lower.tail=FALSE)
  #cat(c(Q.Norm1,l,d, out))
  #cat("\n")
  IDX<-max(which(out < Q.Norm1))

  pval.msg<-sprintf("Pvalue < %e", temp[IDX])
  return(pval.msg)

}
Get_Liu_Params_Mod_Lambda<-function(lambda){
  ## Helper function for getting the parameters for the null approximation

  c1<-rep(0,4)
  for(i in 1:4){
    c1[i]<-sum(lambda^i)
  }

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


