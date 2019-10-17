library(quantreg)
library(SKAT)

error_id = 1

n= 500
pnum = 70
beta_star = 0.3
gamma_star = 0.1

causalrare = 0.3
causalcommon = 0.2

# select region


data("SKAT.haplotypes")
Data = GenerateTestData(error_id, n, pnum, beta_star, gamma_star,
                        causalrare, causalcommon, SKAT.haplotypes)
out = plotiQRAT(Data$y_a, Data$x, Data$c)
plot(out$taus, -log10(out$pval), type = "l", xlab = "quantiles", ylab = "-log10(p)")
