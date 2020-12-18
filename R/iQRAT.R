#' Integrated Quantile Rank Test for group-wise joint effect of rare and common variants in sequencing study
#'
#' This function implements the efficient quantile rank group-wise test, especially designed for testing the joint effect of rare
#' and common variants in sequencing study. Comparing to existing methods, it tests integrated associations across the entire response distribution, rather than
#' only focusing on mean or variability. This method does not have assumptions for error distributions. Thus, no quantile and rank normalization
#' is required for data pre-processing. Note that when applying to univariate, no weights will be assigned.
#'
#' @param X A n by p matrix of genotype.
#' @param C A vector or matrix of covariates.
#' @param v Null model fitted based on quantile process.
#' @param cutoff A number for separating common and rare variants. If not assigned, the default cutoff is \eqn{1/\sqrt(2*n)}.
#' @param weight A logical variable indicating if different weights will be assigned to common and rare variants. The default is TRUE.
#' @param method.type A character, either "S" for \eqn{Q_S} or "B" for \eqn{Q_B}.
#' @param w A vector with length p. It could be functional annotation or other user-specified weights for variants. If w is NULL, the default beta density weight will apply.
#' @param SingleScore Either "Wilcoxon", "Normal", "Lehmann" or "InverseLehmann" if you want to use only one weight function; otherwise, "NULL" represents a combination of four weight functions will be used.
#' @param RareOnly A logical variable indicating if only consider rare variants effect. The default is FALSE. Once RareOnly is set to be TRUE, the cutoff value will be ignored.
#' @param w.dbeta Two parameters for the beta density. The default value is NULL.
#'
#' @return A list contains p value, case and method type. Case indicates whether \code{X} is common, rare or mix.
#'
#' @import quantreg stats
#'
#' @section References:
#' Wang, T., Ionita-Laza, I. and Wei, Y. Integrated Quantile RAnk Test (iQRAT) for gene-level associations.
#' Seunggeun Lee and with contributions from Larisa Miropolsky and Michael Wu.(2017) SKAT: SNP-Set (Sequence) Kernel Association Test.
#' @example man/examples/iQRAT_eg.R
#' @export
iQRAT <- function(X, C, v, cutoff = NULL, weight = TRUE, method.type = "S", w = NULL, SingleScore = NULL, RareOnly = FALSE, w.dbeta = NULL, acc = 1e-9){

  if(is.matrix(X)){
    x_colsum = colSums(X)
    trim_index = which(x_colsum == 0)
    if(length(trim_index) > 0){
      X = as.matrix(X[, -trim_index])
    if(!is.null(w)){
      w = as.vector(w)
      w = w[-trim_index]
    }}

    if(ncol(X) > 1){
      out = iQRAT_mat4(X, C, v, cutoff = cutoff, weight, method.type = method.type, w = w, SingleScore = SingleScore, RareOnly = RareOnly, w.dbeta =w.dbeta, acc =acc)

    }else{
      out = iQRAT_vec4(X, C, v, cutoff = cutoff, method.type = method.type, SingleScore = SingleScore, acc =acc)

    }
  }else{
    out = iQRAT_vec4(X, C, v, cutoff = cutoff, method.type = method.type, SingleScore = SingleScore, acc =acc)
  }

  return(out)
}



iQRAT_vec4<- function(X, C, v, cutoff = NULL, method.type = "S", SingleScore = "NULL", RareOnly = FALSE, acc =1e-9){

  # Test if v comes from quantreg or qr_tau_diff
  if(length(v$diff_dsol) > 0){
    diff_dsol = TRUE
  }else{
    diff_dsol = FALSE
  }
  v = transfer_v(v, diff_dsol = diff_dsol)
  if(sum(X) == 0){
    stop("No variations!")
  }else{
    x_new = X
  }


  n = length(X)
  if(is.null(cutoff)){
    cutoff = 1/sqrt(2*n)
  }

  if(RareOnly){
    cutoff = 999 # we consider all variants as rare, hence only the weights for rare variants will be applied.
  }

  maf = mean(x_new)/2
  ind.err = which(maf>0.5)
  if(length(ind.err) > 0){
    msg<-sprintf("Genotypes of some variants are not the number of minor alleles! These genotypes are flipped!")
    warning(msg,call.=FALSE)

    x_new[,ind.err] = 2 - x_new[, ind.err]

    maf[ind.err] = 1-maf[ind.err]
  }




  Xstar = lm(x_new~C)$residuals




  if(!is.null(SingleScore)){

    if(SingleScore == "Wilcoxon"){
      rr0 = NewRanks( v, score = "wilcoxon")
      null_phi_0 = rr0$ranks
      sd_phi0 = sqrt(rr0$A2)
      null_Sn = (null_phi_0/sd_phi0)%*%Xstar

    }

    if(SingleScore == "Lehmann"){
      if(maf > cutoff){
        rr1 = NewRanks( v, score = "lehmann")
      }else{
        stop("Cannot use Lehmann only for one rare variant!")

      }
      null_phi_1 = rr1$ranks
      sd_phi1 = sqrt(rr1$A2)
      null_Sn  = (null_phi_1/sd_phi1)%*%Xstar

    }
    if(SingleScore == "InverseLehmann"){
      if(maf > cutoff){
        rr1 = NewRanks( v, score = "inverselehmann")
      }else{
        stop("Cannot use InverseLehmann only for one rare variant!")
      }

      null_phi_1 = rr1$ranks
      sd_phi1 = sqrt(rr1$A2)
      null_Sn  = (null_phi_1/sd_phi1)%*%Xstar
    }

    if(SingleScore == "Normal"){
      rr2 = NewRanks( v, score = "normal")
      null_phi_2 = rr2$ranks
      sd_phi2 = sqrt(rr2$A2)
      null_Sn  = (null_phi_2/sd_phi2)%*%Xstar
    }



    if(maf>cutoff){
      out_QS = QS_Single(null_Sn, Xstar, method.p = "davies", acc = acc)
      out = list(pval = out_QS[1], case = "Common1", method.type = method.type )


    }else{

      out_QS = QS_Single(null_Sn, Xstar, method.p = "davies", acc = acc)
      if(SingleScore != "Lehmann" & SingleScore != "InverseLehmann"){
        out = list(pval = out_QS[1], case = "Rare1", method.type = method.type )
      }else{

          out = list(pval = NA, case = "Rare1", method.type = method.type )


      }




    }

  }else{

    # wilcoxon
    rr0 = NewRanks( v, score = "wilcoxon")
    null_phi_0 = rr0$ranks
    sd_phi0 = sqrt(rr0$A2)
    null_Sn_0= (null_phi_0/sd_phi0)%*%Xstar

    # lehmann

    if(maf > cutoff){
      # lehmann
      rr1 = NewRanks( v, score = "lehmann")
      null_phi_1 = rr1$ranks
      sd_phi1 = sqrt(rr1$A2)
      null_Sn_1 = (null_phi_1/sd_phi1)%*%Xstar
    }



    rr2 = NewRanks( v, score = "normal")
    null_phi_2 = rr2$ranks
    sd_phi2 = sqrt(rr2$A2)
    null_Sn_2 = (null_phi_2/sd_phi2)%*%Xstar

    # inverse lehmann
    if(maf > cutoff){
      # lehmann
      rr1 = NewRanks( v, score = "inverselehmann")
      null_phi_1 = rr1$ranks
      sd_phi1 = sqrt(rr1$A2)
      null_Sn_3 = (null_phi_1/sd_phi1)%*%Xstar
    }


    QSKAT_W = sum((null_Sn_0)^2)
    QSKAT_N = sum((null_Sn_2)^2)

    if(maf > cutoff){
      QSKAT_L = sum((null_Sn_1)^2)
      QSKAT_IL = sum((null_Sn_3)^2)
    }


    CovS = crossprod(Xstar)

    if(maf>cutoff){
      pval_S_c = c(davies.pval(CovS, QSKAT_W, acc =acc), davies.pval(CovS, QSKAT_L, acc =acc), davies.pval(CovS, QSKAT_N, acc =acc), davies.pval(CovS, QSKAT_IL, acc =acc))#davies.pval(CovS, c(QSKAT_W, QSKAT_L, QSKAT_N, QSKAT_IL), acc =acc)
      pval_S_CT3_c = pcauchy(mean(tan((0.5-pval_S_c)*pi)), lower.tail = FALSE)

      out = list(pval = pval_S_CT3_c, p_Wilcoxon = pval_S_c[1], p_Normal = pval_S_c[3], p_Lehmann = pval_S_c[2], p_InverseLehmann = pval_S_c[4], case = "Common1", method.type = method.type )
    }else{

      pval_S_c = c(davies.pval(CovS, QSKAT_W, acc =acc),  davies.pval(CovS, QSKAT_N, acc =acc))#davies.pval(CovS, c(QSKAT_W, QSKAT_L, QSKAT_N, QSKAT_IL), acc =acc)
      pval_S_CT3_c = pcauchy(mean(tan((0.5-pval_S_c)*pi)), lower.tail = FALSE)

      out = list(pval = pval_S_CT3_c, p_Wilcoxon = pval_S_c[1], p_Normal = pval_S_c[2], case = "Rare1", method.type = method.type )

    }
  }





  return(out)
}



iQRAT_mat4 <- function(X, C, v, cutoff = NULL, weight = TRUE, method.type = "S", w = NULL, SingleScore = "NULL", RareOnly = FALSE, w.dbeta = NULL, acc =1e-9){


  # Test if v comes from quantreg or qr_tau_diff
  if(length(v$diff_dsol) > 0){
    diff_dsol = TRUE
  }else{
    diff_dsol = FALSE
  }

  v = transfer_v(v, diff_dsol = diff_dsol)

  x_colsum = colSums(X)
  trim_index = which(x_colsum == 0)
  if(length(trim_index) > 0){x_new = as.matrix(X[, -trim_index])
  if(!is.null(w)){
    w = as.vector(w)
    w = w[-trim_index]
  }
  }else {x_new = X}


  n = nrow(X)

  if(is.null(cutoff)){
    cutoff = 1/sqrt(2*n)
  }

  if(RareOnly){
    cutoff = 999 # consider all variants as rare variants
  }

  maf = colMeans(x_new)/2
  ind.err = which(maf>0.5)
  if(length(ind.err) > 0){
    msg<-sprintf("Genotypes of some variants are not the number of minor alleles! These genotypes are flipped!")
    warning(msg,call.=FALSE)

    x_new[,ind.err] = 2 - x_new[, ind.err]

    maf[ind.err] = 1-maf[ind.err]
  }



  ind_common = which(maf>cutoff)
  ind_rare = which(maf <= cutoff)

  Xstar1 = lm(x_new~C)$residuals

  ### common and rare combined case
  if(length(ind_common) > 0 & length(ind_rare) > 0){

    Xstar_c = Xstar1[, ind_common]
    Xstar_r = Xstar1[, ind_rare]


    if(weight == TRUE){

      if(is.null(w)){
        if(!is.null(w.dbeta)){
          W_c = dbeta(maf[ind_common], w.dbeta[1],w.dbeta[2])
          W_r = dbeta(maf[ind_rare], w.dbeta[1],w.dbeta[2])
        }else{
          W_c = dbeta(maf[ind_common], 0.5,0.5)
          W_r = dbeta(maf[ind_rare], 1,25)
        }

      }else{
        W_c = as.vector(w[ind_common])
        W_r = as.vector(w[ind_rare])
      }


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

    if(!is.null(SingleScore)){
      if(SingleScore == "Wilcoxon"){
        # wilcoxon
        rr0 = NewRanks( v, score = "wilcoxon")
        null_phi_0 = rr0$ranks

        sd_phi0 = sqrt(rr0$A2)
        null_Sn_c = (null_phi_0/sd_phi0)%*%Xstar_c
        null_Sn_r = (null_phi_0/sd_phi0)%*%Xstar_r
      }

      if(SingleScore == "Lehmann"){
        # lehmann
        rr1 = NewRanks( v, score = "lehmann")
        null_phi_1 = rr1$ranks
        sd_phi1 = sqrt(rr1$A2)
        null_Sn_c = (null_phi_1/sd_phi1)%*%Xstar_c


      }
      if(SingleScore == "InverseLehmann"){
        rr1 = NewRanks( v, score = "inverselehmann")
        null_phi_1 = rr1$ranks
        sd_phi1 = sqrt(rr1$A2)
        null_Sn_c = (null_phi_1/sd_phi1)%*%Xstar_c



      }

      if(SingleScore == "Normal"){
        # normal
        rr2 = NewRanks( v, score = "normal")
        null_phi_2 = rr2$ranks
        sd_phi2 = sqrt(rr2$A2)
        null_Sn_c = (null_phi_2/sd_phi2)%*%Xstar_c
        null_Sn_r = (null_phi_2/sd_phi2)%*%Xstar_r
      }

      # calculate common and rare p values
      if(method.type == "S"){
        out_c = QS_Single(null_Sn = null_Sn_c, Xstar = Xstar_c, method.p = "davies", acc=acc)
        if(SingleScore == "Normal" | SingleScore == "Wilcoxon"){
          out_r = QS_Single(null_Sn = null_Sn_r, Xstar = Xstar_r, method.p = "davies", acc=acc)
        }else{
          out_r = NULL
        }

      }

      if(method.type == "B"){
        out_c = QB_Single(null_Sn_c,  Xstar_c, ind_common, method.p = "davies", acc=acc)
        if(SingleScore == "Normal" | SingleScore == "Wilcoxon"){
          out_r = QB_Single(null_Sn_r,  Xstar_r, ind_rare, method.p = "davies", acc=acc)
        }else{
          out_r = NULL
        }

      }


      # use cauchy to combine common and rare p values
      if(!is.null(out_r)){
        p_vec = c(out_c, out_r)
      }else{
        p_vec = out_c
      }

      # use Cauchy to combine p value
      pval_all = pcauchy(mean(tan((0.5-p_vec)*pi)), lower.tail = FALSE)

      out = list(pval = pval_all,  case = "Mix", method.type = method.type)


    }else{

      # wilcoxon
      rr0 = NewRanks( v, score = "wilcoxon")
      null_phi_0 = rr0$ranks
      sd_phi0 = sqrt(rr0$A2)

      null_Sn_0_c = (null_phi_0/sd_phi0)%*%Xstar_c
      null_Sn_0_r = (null_phi_0/sd_phi0)%*%Xstar_r

      # lehmann
      rr1 = NewRanks( v, score = "lehmann")
      null_phi_1 = rr1$ranks
      sd_phi1 = sqrt(rr1$A2)

      null_Sn_1_c = (null_phi_1/sd_phi1)%*%Xstar_c



      # normal
      rr2 = NewRanks( v, score = "normal")
      null_phi_2 = rr2$ranks
      sd_phi2 = sqrt(rr2$A2)

      null_Sn_2_c = (null_phi_2/sd_phi2)%*%Xstar_c
      null_Sn_2_r = (null_phi_2/sd_phi2)%*%Xstar_r


      # inverse lehmann
      rr1 = NewRanks( v, score = "inverselehmann")
      null_phi_1 = rr1$ranks
      sd_phi1 = sqrt(rr1$A2)


      null_Sn_3_c = (null_phi_1/sd_phi1)%*%Xstar_c



      ### Multiple common and multiple rare
      if(method.type == "S"){

        QSKAT_W = sum((null_Sn_0_c)^2)
        QSKAT_L = sum((null_Sn_1_c)^2)
        QSKAT_N = sum((null_Sn_2_c)^2)
        QSKAT_IL = sum((null_Sn_3_c)^2)

        CovS = crossprod(Xstar_c)
        pval_S_c = c(davies.pval(CovS, QSKAT_W, acc =acc), davies.pval(CovS, QSKAT_L, acc =acc), davies.pval(CovS, QSKAT_N, acc =acc), davies.pval(CovS, QSKAT_IL, acc =acc))#davies.pval(CovS, c(QSKAT_W, QSKAT_L, QSKAT_N, QSKAT_IL), acc =acc)
        pval_S_CT3_c = pcauchy(mean(tan((0.5-pval_S_c)*pi)), lower.tail = FALSE)


        QSKAT_W = sum((null_Sn_0_r)^2)

        QSKAT_N = sum((null_Sn_2_r)^2)


        CovS = crossprod(Xstar_r)
        pval_S_r = c(davies.pval(CovS, QSKAT_W, acc =acc), davies.pval(CovS, QSKAT_N, acc =acc))
        pval_S_CT3_r = pcauchy(mean(tan((0.5-pval_S_r)*pi)), lower.tail = FALSE)




      }

      if(method.type == "B"){

        ### Commmon
        QSKAT_W = (sum(null_Sn_0_c))^2
        QSKAT_L = (sum(null_Sn_1_c))^2
        QSKAT_N = (sum(null_Sn_2_c))^2
        QSKAT_IL = (sum(null_Sn_3_c))^2


        if(length(ind_common) ==1){
          Xstar_c = matrix(Xstar_c, nrow=length(Xstar_c), ncol = 1)
        }
        CovS = crossprod(Xstar_c%*%rep(1, length(ind_common)))
        pval_S_c = c(davies.pval(CovS, QSKAT_W, acc =acc), davies.pval(CovS, QSKAT_L, acc =acc), davies.pval(CovS, QSKAT_N, acc =acc), davies.pval(CovS, QSKAT_IL, acc =acc))#davies.pval(CovS, c(QSKAT_W, QSKAT_L, QSKAT_N, QSKAT_IL), acc =acc)
        pval_S_CT3_c = pcauchy(mean(tan((0.5-pval_S_c)*pi)), lower.tail = FALSE)

        ### Rare
        QSKAT_W = (sum(null_Sn_0_r))^2

        QSKAT_N = (sum(null_Sn_2_r))^2



        if(length(ind_rare) ==1){
          Xstar_r = matrix(Xstar_r, nrow=length(Xstar_r), ncol = 1)
        }
        CovS = crossprod(Xstar_r%*%rep(1, length(ind_rare)))
        pval_S_r = c(davies.pval(CovS, QSKAT_W, acc =acc), davies.pval(CovS, QSKAT_N, acc =acc))#davies.pval(CovS, c(QSKAT_W, QSKAT_L, QSKAT_N, QSKAT_IL), acc =acc)
        pval_S_CT3_r = pcauchy(mean(tan((0.5-pval_S_r)*pi)), lower.tail = FALSE)
      }
      #

      p_vec = as.numeric(c(pval_S_CT3_r, pval_S_CT3_c))
      # use Cauchy to combine p value
      pval_all = pcauchy(mean(tan((0.5-p_vec)*pi)), lower.tail = FALSE)
      out = list(pval = pval_all,  p_Wilcoxon =  pcauchy(mean(tan((0.5-c(pval_S_c[1],pval_S_r[1]))*pi)), lower.tail = FALSE),
                 p_Normal =  pcauchy(mean(tan((0.5-c(pval_S_c[3],pval_S_r[2]))*pi)),lower.tail = FALSE),
                 p_Lehmann =  pval_S_c[2],
                 p_InverseLehamnn =  pval_S_c[4],case = "Mix", method.type = method.type)

    }

  }else if(length(ind_common) == 0 & length(ind_rare) > 0){

    ##### Rare only case

    Xstar = Xstar1
    if(weight == TRUE){
      if(is.null(w)){
        W = dbeta(maf, 1,25)
      }else{
        W = as.vector(w)
      }


      Xstar = t(t(Xstar)*W)
    }



    if(!is.null(SingleScore)){

      if(SingleScore == "Wilcoxon"){
        rr0 = NewRanks( v, score = "wilcoxon")
        null_phi_0 = rr0$ranks
        sd_phi0 = sqrt(rr0$A2)
        null_Sn = (null_phi_0/sd_phi0)%*%Xstar

        if(method.type == "S"){

          out_r = QS_Single(null_Sn, Xstar, method.p = "davies", acc =acc)
        }

        if(method.type == "B"){

          out_r = QB_Single(null_Sn,  Xstar, ind_rare, method.p = "davies", acc =acc)
        }
      }




      if(SingleScore == "Normal"){
        rr2 = NewRanks( v, score = "normal")
        null_phi_2 = rr2$ranks
        sd_phi2 = sqrt(rr2$A2)
        null_Sn = (null_phi_2/sd_phi2)%*%Xstar

        if(method.type == "S"){

          out_r = QS_Single(null_Sn, Xstar, method.p = "davies", acc =acc)
        }

        if(method.type == "B"){

          out_r = QB_Single(null_Sn,  Xstar, ind_rare, method.p = "davies", acc =acc)
        }
      }

      if(SingleScore == "Wilcoxon" | SingleScore == "Normal"){
        out = list(pval = out_r[1], case = "RareOnly", method.type = method.type )
      }else{
        out = list(pval = NA, case = "RareOnly", method.type = method.type )
      }



    }else{

      # wilcoxon
      rr0 = NewRanks( v, score = "wilcoxon")
      null_phi_0 = rr0$ranks
      sd_phi0 = sqrt(rr0$A2)
      null_Sn_0 = (null_phi_0/sd_phi0)%*%Xstar


      # normal
      rr2 = NewRanks( v, score = "normal")
      null_phi_2 = rr2$ranks
      sd_phi2 = sqrt(rr2$A2)
      null_Sn_2 = (null_phi_2/sd_phi2)%*%Xstar


      if(method.type == "S"){
        QSKAT_W = sum((null_Sn_0)^2)
        QSKAT_N = sum((null_Sn_2)^2)

        CovS = crossprod(Xstar)
        pval_S_r = c(davies.pval(CovS, QSKAT_W, acc =acc), davies.pval(CovS, QSKAT_N, acc =acc))#davies.pval(CovS, c(QSKAT_W, QSKAT_L, QSKAT_N, QSKAT_IL), acc =acc)
        pval_S_CT3_r = pcauchy(mean(tan((0.5-pval_S_r)*pi)), lower.tail = FALSE)
      }
      if(method.type == "B"){
        QSKAT_W = (sum(null_Sn_0))^2
        QSKAT_N = (sum(null_Sn_2))^2


        if(length(ind_rare) ==1){
          Xstar_r = matrix(Xstar, nrow=length(Xstar), ncol = 1)
        }
        CovS = crossprod(Xstar%*%rep(1, length(ind_rare)))
        pval_S_r = c(davies.pval(CovS, QSKAT_W, acc =acc),  davies.pval(CovS, QSKAT_N, acc =acc))#davies.pval(CovS, c(QSKAT_W, QSKAT_L, QSKAT_N, QSKAT_IL), acc =acc)
        pval_S_CT3_r = pcauchy(mean(tan((0.5-pval_S_r)*pi)), lower.tail = FALSE)


      }
      out = list(pval = pval_S_CT3_r, p_Wilcoxon = pval_S_r[1], p_Normal = pval_S_r[2], case = "RareOnly", method.type = method.type )

    }



  }else{
    ### Common only case

    Xstar_c = Xstar1



    if(weight == TRUE){
      if(is.null(w)){
        W_c = dbeta(maf, 0.5,0.5)
      }else{
        W_c = as.vector(w)
      }

      Xstar_c = t(t(Xstar_c)*W_c)
    }


    if(!is.null(SingleScore)){

      if(SingleScore == "Wilcoxon"){
        rr0 = NewRanks( v, score = "wilcoxon")
        null_phi_0 = rr0$ranks
        sd_phi0 = sqrt(rr0$A2)
        null_Sn = (null_phi_0/sd_phi0)%*%Xstar_c

      }

      if(SingleScore == "Lehmann"){
        rr1 = NewRanks( v, score = "lehmann")
        null_phi_1 = rr1$ranks
        sd_phi1 = sqrt(rr1$A2)
        null_Sn = (null_phi_1/sd_phi1)%*%Xstar_c

      }

      if(SingleScore == "InverseLehmann"){
        rr1 = NewRanks( v, score = "inverselehmann")
        null_phi_1 = rr1$ranks
        sd_phi1 = sqrt(rr1$A2)
        null_Sn = (null_phi_1/sd_phi1)%*%Xstar_c

      }

      if(SingleScore == "Normal"){
        rr2 = NewRanks( v, score = "normal")
        null_phi_2 = rr2$ranks
        sd_phi2 = sqrt(rr2$A2)
        null_Sn = (null_phi_2/sd_phi2)%*%Xstar_c

      }
      if(method.type == "S"){
        out_c = QS_Single(null_Sn, Xstar_c, method.p = "davies", acc =acc)

      }

      if(method.type == "B"){
        out_c = QB_Single(null_Sn,  Xstar_c, ind_common, method.p = "davies", acc =acc)

      }

      out = list(pval = out_c[1],case = "CommonOnly", method.type = method.type )

    }else{
      rr0 = NewRanks( v, score = "wilcoxon")
      null_phi_0 = rr0$ranks
      sd_phi0 = sqrt(rr0$A2)
      null_Sn_0 = (null_phi_0/sd_phi0)%*%Xstar_c

      rr1 = NewRanks( v, score = "lehmann")
      null_phi_1 = rr1$ranks
      sd_phi1 = sqrt(rr1$A2)
      null_Sn_1 = (null_phi_1/sd_phi1)%*%Xstar_c

      rr2 = NewRanks( v, score = "normal")
      null_phi_2 = rr2$ranks
      sd_phi2 = sqrt(rr2$A2)
      null_Sn_2 = (null_phi_2/sd_phi2)%*%Xstar_c

      rr3 = NewRanks( v, score = "inverselehmann")
      null_phi_3 = rr3$ranks
      sd_phi3 = sqrt(rr3$A2)
      null_Sn_3 = (null_phi_3/sd_phi3)%*%Xstar_c

      if(method.type == "S"){
        QSKAT_W = sum((null_Sn_0)^2)
        QSKAT_L = sum((null_Sn_1)^2)
        QSKAT_N = sum((null_Sn_2)^2)
        QSKAT_IL = sum((null_Sn_3)^2)

        CovS = crossprod(Xstar_c)
        pval_S_c = c(davies.pval(CovS, QSKAT_W, acc =acc), davies.pval(CovS, QSKAT_L, acc =acc), davies.pval(CovS, QSKAT_N, acc =acc), davies.pval(CovS, QSKAT_IL, acc =acc))#davies.pval(CovS, c(QSKAT_W, QSKAT_L, QSKAT_N, QSKAT_IL), acc =acc)
        pval_S_CT3_c = pcauchy(mean(tan((0.5-pval_S_c)*pi)), lower.tail = FALSE)

      }

      if(method.type ==  "B"){
        QSKAT_W = (sum(null_Sn_0))^2
        QSKAT_L = (sum(null_Sn_1))^2
        QSKAT_N = (sum(null_Sn_2))^2
        QSKAT_IL = (sum(null_Sn_3))^2


        if(length(ind_common) ==1){
          Xstar_c = matrix(Xstar_c, nrow=length(Xstar_c), ncol = 1)
        }
        CovS = crossprod(Xstar_c%*%rep(1, length(ind_common)))
        pval_S_c = c(davies.pval(CovS, QSKAT_W, acc =acc), davies.pval(CovS, QSKAT_L, acc =acc), davies.pval(CovS, QSKAT_N, acc =acc), davies.pval(CovS, QSKAT_IL, acc =acc))#davies.pval(CovS, c(QSKAT_W, QSKAT_L, QSKAT_N, QSKAT_IL), acc =acc)
        pval_S_CT3_c = pcauchy(mean(tan((0.5-pval_S_c)*pi)), lower.tail = FALSE)

      }


      out = list(pval = pval_S_CT3_c, p_Wilcoxon = pval_S_c[1], p_Normal = pval_S_c[3], p_Lehmann = pval_S_c[2], p_InverseLehmann = pval_S_c[4], case = "CommonOnly", method.type = method.type )


    }

  }


  return(out)
}
