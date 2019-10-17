#' Visualize quantile effect for the entire process
#'
#' This function implements the visualization tool as an extension of the main approach.
#'
#' @param Y A vector of response, size \eqn{n}.
#' @param X A n by p matrix of genotype.
#' @param C A vector or matrix of covariates.
#' @return A list contains quantile levels and corresponding p values.
#'
#' @example man/examples/plotiQRAT_eg.R
#' @export

plotiQRAT<- function(Y, X, C){
  if(is.matrix(X)){
    x_colsum = colSums(X)
    trim_index = which(x_colsum == 0)
    if(length(trim_index) > 0){x_new = X[, -trim_index]
    }else {x_new = X}
    #maf = colMeans(x_new)/2
    #tr1 = which(maf>0.05)
    # x_new = x_new[,tr1]
    maf = colMeans(x_new)/2
    ind.err = which(maf>0.5)
  }else{
    x_new = X
    maf = mean(x_new)/2
    ind.err = which(maf>0.5)
  }


  n = length(Y)
cutoff = 1/sqrt(2*n)
  if(length(ind.err) > 0){
    msg<-sprintf("Genotypes of some variants are not the number of minor alleles! These genotypes are flipped!")
    warning(msg,call.=FALSE)

    x_new[,ind.err] = 2 - x_new[, ind.err]

    maf[ind.err] = 1-maf[ind.err]
  }

  v = rq(Y~C,tau = -1)



  Xstar1 = lm(x_new~C)$residuals


  if(is.matrix(x_new)){


    ind_common = which(maf>cutoff)
    ind_rare = which(maf <= cutoff)
    Xstar_c = Xstar1[, ind_common]
    Xstar_r = Xstar1[, ind_rare]
    Xstar_r = lm(Xstar_r~Xstar_c)$residuals


    W_c = dbeta(maf[ind_common], 0.5,0.5)
    W_r = dbeta(maf[ind_rare], 1,25)
    Xstar_c = t(t(Xstar_c)*W_c)
    Xstar_r = t(t(Xstar_r)*W_r)

    taus <- v$sol[1, ]
    J = length(taus)
    taumat = matrix(rep(taus, n), nrow = n, ncol = length(taus), byrow = TRUE)
    ranksmat = taumat - (1 - v$dsol)
    ranks_1 = apply(ranksmat, 2, sd)
    #ranks_1 = sqrt(taus*(1-taus))

    ranks_sd = ranks_1[2:(J-1)]
    ranks_int = ranksmat[,2:(J-1)]
    null_phi_0 = (ranks_int/ranks_sd)
    #null_phi_0 = cbind(ranksmat[,1], null_phi_0, ranksmat[,J])
    null_Sn_0_c =t(Xstar_c)%*% null_phi_0
    null_Sn_0_r =t(Xstar_r)%*% null_phi_0




    V0 = colSums((null_Sn_0_c)^2)
    CovS = crossprod(Xstar_c)
    pval_c = Liu.pval(CovS, V0)

    V0 = colSums((null_Sn_0_r)^2)
    CovS = crossprod(Xstar_r)
    pval_r = Liu.pval(CovS, V0)


    TFisher = -2*log(pval_c)-2*log(pval_r)
    pval_Fisher = pchisq(TFisher, df = 4, lower.tail = FALSE)
  }else{


    taus <- v$sol[1, ]
    J = length(taus)
    taumat = matrix(rep(taus, n), nrow = n, ncol = length(taus), byrow = TRUE)
    ranksmat = taumat - (1 - v$dsol)
    ranks_1 = apply(ranksmat, 2, sd)
    #ranks_1 = sqrt(taus*(1-taus))

    ranks_sd = ranks_1[2:(J-1)]
    ranks_int = ranksmat[,2:(J-1)]
    null_phi_0 = (ranks_int/ranks_sd)
    null_Sn_0 =t(Xstar1)%*% null_phi_0




    V0 = colSums((null_Sn_0)^2)
    CovS = crossprod(Xstar1)
    pval_Fisher = Liu.pval(CovS, V0)
  }


  return(list(pval = pval_Fisher, taus = taus[2:(J-1)]))
}
