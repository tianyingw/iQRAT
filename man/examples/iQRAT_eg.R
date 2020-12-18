
rm(list=ls())
library(quantreg)
library(SKAT)
data("SampleData")

# Step 1: fit null model
null.fit = Null_model(Y = SampleData$y, C = SampleData$c)
# Step 2: run the test, p value will return

# SKAT version iQRAT
test.iQRAT1 = iQRAT(X = SampleData$x, C = SampleData$c, v =  null.fit, method.type = "S")

# If you want to specify weights
w = dbeta(colMeans(SampleData$x)/2,0.5,0.5) # we use beta density as an example
test.iQRAT2 = iQRAT(X = SampleData$x, C = SampleData$c, v =  null.fit, method.type = "S", w = w)

# Burden version iQRAT
test.iQRAT3 = iQRAT(X = SampleData$x, C = SampleData$c, v =  null.fit, method.type = "B")


