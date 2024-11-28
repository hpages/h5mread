### This is the code that was used to produce the test.h5 file located
### in the extdata/ folder of the h5mread package.

library(HDF5Array)

set.seed(123)

m1 <- matrix(101:160, ncol=5)

m2 <- matrix((runif(360000) - 0.5) * 10, ncol=90)
m2[1:2, 2:3] <- c(NA, Inf, NaN, -Inf)
rownames(m2) <- sprintf("row%04d", 1:4000)
colnames(m2) <- sprintf("col%02d", 1:90)

a3 <- poissonSparseArray(c(180, 75, 4), density=0.15)

m4 <- matrix(as.raw(sample(0:255, 112000, replace=TRUE)), nrow=28)

## 30k random words:
rwords <- sapply(1:3e4, function(i) paste(sample(letters, sample(10, 1)), collapse=""))

writeHDF5Array(m1, "test.h5", "m1")
writeHDF5Array(m2, "test.h5", "m2", chunkdim=c(50, 40))
writeHDF5Array(a3, "test.h5", "a3", chunkdim=c(20, 20, 4))
writeHDF5Array(m4, "test.h5", "m4", chunkdim=c(28, 250))
h5write(rwords, "test.h5", "rwords")

h5ls("test.h5")

