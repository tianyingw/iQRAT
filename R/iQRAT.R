#' Integrated Quantile Rank Test for group-wise joint effect of rare and common variants in sequencing study
#'
#' This function implements the efficient quantile rank group-wise test, especially designed for testing the joint effect of rare
#' and common variants in sequencing study. Comparing to existing methods, it tests integrated associations across the entire response distribution, rather than
#' only focusing on mean or variability. This method does not have assumptions for error distributions. Thus, no quantile and rank normalization
#' is required for data pre-processing. Note that when applying to univariate, no weights will be assigned.
#'
#' @param Y A vector of response, size \eqn{n}.
#' @param X A n by p matrix of genotype.
#' @param C A vector or matrix of covariates.
#' @param cutoff A number for separating common and rare variants. If not assigned, the default cutoff is \eqn{1/\sqrt(2*n)}.
#' @param weight A logicalvariable indicating if different weights will be assigned to common and rare variants. The default is TRUE.
#' @param method.type A character, either "S" for \eqn{Q_S} or "B" for \eqn{Q_B}.
#'
#' @return A list contains p value, case and method type. Case indicates whether \code{X} is common, rare or mix.
#'
#' @import quantreg stats
#'
#' @section References:
#' Wang, T., Ionita-Laza, I. and Wei, Y. Efficient Quantile Rank Test for heterogeneous joint effect of rare and common variants in sequencing studies.
#'
#' @example man/examples/iQRAT_eg.R
#' @export
iQRAT <- function(Y, X, C, cutoff = NULL, weight = TRUE, method.type = "S"){

  if(is.matrix(X)){
    out = iQRAT_mat(Y, X, C, cutoff = cutoff, weight, method.type = method.type)
  }else{
    out = iQRAT_vec(Y, X, C, cutoff = cutoff, method.type = method.type)
  }

  return(out)
}



iQRAT_vec<- function(Y, X, C, cutoff = NULL, method.type = "S"){


  if(sum(X) == 0){
    stop("No variations!")
  }else{
    x_new = X
  }


  n = length(Y)
  if(is.null(cutoff)){
    cutoff = 1/sqrt(2*n)
  }
  maf = mean(x_new)/2
  ind.err = which(maf>0.5)
  if(length(ind.err) > 0){
    msg<-sprintf("Genotypes of some variants are not the number of minor alleles! These genotypes are flipped!")
    warning(msg,call.=FALSE)

    x_new[,ind.err] = 2 - x_new[, ind.err]

    maf[ind.err] = 1-maf[ind.err]
  }

  v = rq(Y~C,tau = -1)


  Xstar = lm(x_new~C)$residuals




  rr = ranks( v, score = "wilcoxon")
  null_phi_0 = rr$ranks

  sd_phi0 = sqrt(rr$A2)


  null_Sn_0= (null_phi_0/sd_phi0)%*%Xstar


  rr = ranks( v, score = "lehmann")


  null_phi_1 = rr$ranks

  sd_phi1 = sqrt(rr$A2)


  null_Sn_1 = (null_phi_1/sd_phi1)%*%Xstar


  rr = ranks( v, score = "normal")


  null_phi_2 = rr$ranks

  sd_phi2 = sqrt(rr$A2)

  null_Sn_2 = (null_phi_2/sd_phi2)%*%Xstar

  if(maf>cutoff){
    out_QS_c = QS_Common(null_Sn_0, null_Sn_1, null_Sn_2, Xstar)
    out_QS_c = as.numeric(out_QS_c)


    out = list(pval = out_QS_c[1], case = "Common1", method.type = method.type )


  }else{

    out_QS_r = QS_Rare(null_Sn_0, null_Sn_2, Xstar)
    out_QS_r = as.numeric(out_QS_r)

    out = list( pval = out_QS_r[1], case = "Rare1", method.type = method.type )



  }


  return(out)
}



iQRAT_mat <- function(Y, X, C, cutoff = NULL, weight = TRUE, method.type = "S"){


  x_colsum = colSums(X)
  trim_index = which(x_colsum == 0)
  if(length(trim_index) > 0){x_new = X[, -trim_index]
  }else {x_new = X}


  n = length(Y)

  if(is.null(cutoff)){
    cutoff = 1/sqrt(2*n)
  }
  maf = colMeans(x_new)/2
  ind.err = which(maf>0.5)
  if(length(ind.err) > 0){
    msg<-sprintf("Genotypes of some variants are not the number of minor alleles! These genotypes are flipped!")
    warning(msg,call.=FALSE)

    x_new[,ind.err] = 2 - x_new[, ind.err]

    maf[ind.err] = 1-maf[ind.err]
  }

  v = rq(Y~C,tau = -1)

  ind_common = which(maf>cutoff)
  ind_rare = which(maf <= cutoff)

  Xstar1 = lm(x_new~C)$residuals

  if(length(ind_common) > 0 & length(ind_rare) > 0){



    Xstar_c = Xstar1[, ind_common]
    Xstar_r = Xstar1[, ind_rare]
    Xstar_r = lm(Xstar_r~Xstar_c)$residuals


    if(weight == TRUE){
      W_c = dbeta(maf[ind_common], 0.5,0.5)
      W_r = dbeta(maf[ind_rare], 1,25)

      if(length(ind_common) == 1){
        Xstar_c = Xstar_c*W_c
      }else{
        Xstar_c = t(t(Xstar_c)*W_c)
      }
      if(length(ind_rare) == 1){
        Xstar_r = Xstar_r*W_r
      }else{
        Xstar_r = t(t(Xstar_r)*W_r)
      }

    }

    rr = ranks( v, score = "wilcoxon")
    null_phi_0 = rr$ranks

    sd_phi0 = sqrt(rr$A2)


    null_Sn_0_c = (null_phi_0/sd_phi0)%*%Xstar_c
    null_Sn_0_r = (null_phi_0/sd_phi0)%*%Xstar_r


    rr = ranks( v, score = "lehmann")


    null_phi_1 = rr$ranks

    sd_phi1 = sqrt(rr$A2)


    null_Sn_1_c = (null_phi_1/sd_phi1)%*%Xstar_c


    rr = ranks( v, score = "normal")


    null_phi_2 = rr$ranks

    sd_phi2 = sqrt(rr$A2)

    null_Sn_2_c = (null_phi_2/sd_phi2)%*%Xstar_c
    null_Sn_2_r = (null_phi_2/sd_phi2)%*%Xstar_r



    ### Multiple common and multiple rare
    if(method.type == "S"){
      out_QS_c = QS_Common(null_Sn_0_c, null_Sn_1_c, null_Sn_2_c, Xstar_c)
      out_QS_r = QS_Rare(null_Sn_0_r, null_Sn_2_r, Xstar_r)

      out_QS_Combo = Q_Combine(out_QS_c, out_QS_r)
      out = list(pval = out_QS_Combo$pval3,  case = "Mix", method.type = method.type)
    }

    if(method.type == "B"){

      out_QB_c = QB_Common(null_Sn_0_c, null_Sn_1_c, null_Sn_2_c, Xstar_c, ind_common)
      out_QB_r = QB_Rare(null_Sn_0_r, null_Sn_2_r, Xstar_r, ind_rare)

      out_QB_Combo = Q_Combine(out_QB_c, out_QB_r)
      out = list(pval = out_QB_Combo$pval3,  case = "Mix", method.type = method.type)
    }



  }else if(length(ind_common) == 0 & length(ind_rare) > 0){

    ##### Rare only case

    Xstar = Xstar1
    if(weight == TRUE){

      W = dbeta(maf, 1,25)
      Xstar = t(t(Xstar)*W)
    }

    rr = ranks( v, score = "wilcoxon")
    null_phi_0 = rr$ranks

    sd_phi0 = sqrt(rr$A2)
    null_Sn_0 = (null_phi_0/sd_phi0)%*%Xstar


    rr = ranks( v, score = "normal")
    null_phi_2 = rr$ranks

    sd_phi2 = sqrt(rr$A2)
    null_Sn_2 = (null_phi_2/sd_phi2)%*%Xstar

    if(method.type == "S"){
      out_QS_r = as.numeric(QS_Rare(null_Sn_0, null_Sn_2, Xstar))
      out = list(pval = out_QS_r[1], case = "RareOnly", method.type = method.type )

    }
    if(method.type == "B"){
      out_QB_r = as.numeric(QB_Rare(null_Sn_0, null_Sn_2, Xstar, ind_rare))
      out = list(pval = out_QB_r[1],case = "RareOnly", method.type = method.type)

    }





  }else{
    ### Common only case

    Xstar_c = Xstar1



    if(weight == TRUE){

      W_c = dbeta(maf, 0.5,0.5)
      Xstar_c = t(t(Xstar_c)*W_c)
    }

    rr = ranks( v, score = "wilcoxon")
    null_phi_0 = rr$ranks

    sd_phi0 = sqrt(rr$A2)



    null_Sn_0_c = (null_phi_0/sd_phi0)%*%Xstar_c

    rr = ranks( v, score = "lehmann")
    null_phi_1 = rr$ranks
    sd_phi1 = sqrt(rr$A2)
    null_Sn_1_c = (null_phi_1/sd_phi1)%*%Xstar_c


    rr = ranks( v, score = "normal")


    null_phi_2 = rr$ranks

    sd_phi2 = sqrt(rr$A2)

    null_Sn_2_c = (null_phi_2/sd_phi2)%*%Xstar_c

    if(method.type == "S"){
      out_QS_c = as.numeric(QS_Common(null_Sn_0_c, null_Sn_1_c, null_Sn_2_c, Xstar_c))
      out = list(pval = out_QS_c[1],case = "CommonOnly", method.type = method.type )

    }

    if(method.type ==  "B"){
      out_QB_c = as.numeric(QB_Common(null_Sn_0_c, null_Sn_1_c, null_Sn_2_c, Xstar_c, ind_common))
      out = list(pval = out_QB_c[1],case = "CommonOnly", method.type = method.type )

    }

  }




  return(out)
}
