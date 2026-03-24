library(lattice)
library(MASS)
library(nnet)
library(mice)

Part <- All_dataframe[, c(1:31)]

data_mice <- mice(Part, m = 5, maxit = 50, meth = 'pmm', seed = 500)
summary(data_mice)
result <- complete(data_mice, action = 3)
sum(is.na(result))

All_dataframe2 <- left_join(result, All_dataframe[, c(1, 32:52)])