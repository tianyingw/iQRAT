#' Function to generate a test data
#'
#' @param error_id A number indicates error distribution: 1 for N(0,1); 2 for Chi square with df=2; 3 for Cauchy (0,1); 4 for t(df = 2); 5 for t(df = 3).
#' @param n Sample size.
#' @param pnum Number of variants.
#' @param beta_star coefficients for \code{X}.
#' @param gamma_star Control heteroscedastic case: \code{error = (1+gamma_star*X)*e}.
#' @param causalrare Causal ratio for rare variants.
#' @param causalcommon Causal ratio for common variants.
#' @param SKAT.haplotypes Data from SKAT pacakge.
#'
#' @return A list contains response \code{y} under null distribution, \code{y} under alternative distrubtion, \code{x} genotype, \code{c} covariates and \code{error}.
#'
#' @import SKAT
#' @export
#'
GenerateTestData <- function(error_id, n, pnum, beta_star, gamma_star, causalrare, causalcommon, SKAT.haplotypes){


  hp = SKAT.haplotypes$Haplotype
  nhp = nrow(hp)
  maf = colSums(hp)/nrow(hp)
  CHROM_POS = SKAT.haplotypes$SNPInfo$CHROM_POS
  rm(SKAT.haplotypes)


  set.seed(12)
  index50 = sample(ncol(hp)-pnum, 1)
  CHROM_select = CHROM_POS[index50:(index50 + pnum -1)]
  selected.region = range(CHROM_select)
  region.size = selected.region[2] - selected.region[1]

  x_index = sample(nhp, 2*n)


  hp50 = hp[x_index, index50:(index50 + pnum -1)]
  maf50 = colMeans(hp50)
  mafrare = which(maf50 <= 0.05 & maf50 > 0)
  mafcom = which(maf50 > 0.05)
  pnum1 = length(mafrare)
  pnum2 = length(mafcom)

  while(pnum2 < 10){
    index50 = sample(ncol(hp)-pnum, 1)
    CHROM_select = CHROM_POS[index50:(index50 + pnum -1)]
    selected.region = range(CHROM_select)
    region.size = selected.region[2] - selected.region[1]

    if(n>5000){
      x_index = sample(nhp, 2*n, replace = TRUE)

    }else{
      x_index = sample(nhp, 2*n)
    }
    hp50 = hp[x_index, index50:(index50 + pnum -1)]
    maf50 = colMeans(hp50)
    mafrare = which(maf50 <= 0.05 & maf50>0)
    mafcom = which(maf50 > 0.05)
    pnum1 = length(mafrare)
    pnum2 = length(mafcom)
  }


  causal_index1 = sample(mafcom, max(pnum2*causalcommon, 1))
  causal_index2 = sample(mafrare, max(pnum1*causalrare, 1))
  causal_index = c(causal_index1, causal_index2)

  beta = gamma = rep(0, pnum)
  beta[causal_index]= beta_star * abs(log10(maf50[causal_index]))
  gamma[causal_index]= gamma_star * abs(log10(maf50[causal_index]))



  output_all = list()
  output_all$region.size = region.size
  output_all$indexSelect = index50
  output_all$beta = beta
  output_all$gamma = gamma
  output_all$x_index = x_index
  output_all$causal_index = causal_index
  output_all$pnum_rare = pnum1
  output_all$pnum_common = pnum2

  output_all$causal_index1 = causal_index1
  output_all$causal_index2 = causal_index2

  temp_list = output_all

  C=rnorm(n,mean=4,sd=1)
  if (error_id == 1){
    error = rnorm(n)
  }
  if (error_id == 2){
    error = rchisq(n, 2)
  }
  if (error_id == 3){
    error = rcauchy(n)
  }
  if (error_id == 4){
    error = rt(n, 2)
  }
  if (error_id ==5){
    error = rt(n, 3)
  }

  temp_list2 = list()
  temp_list2$C = C
  temp_list2$error = error

  Data = GenerateSample(hp, temp_list, temp_list2, pnum, local = FALSE, error_id)

  return(Data)
}
