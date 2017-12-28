context('Test utility functions')

test_that('commaList2vector returns vector from comma-separated inputlist', {
    expect_equal(commaList2vector("1,2,3"), c(1,2,3))
})

test_that('addNonNulls fails when dimensions of list elements are not the same', 
          {
    expect_error(addNonNulls(list(matrix(1:10, ncol=2), matrix(11:20, ncol=2), 
                                  matrix(21:30, ncol=5), NULL, NULL)),
                 "Column dimensions of list elements are different")
    expect_error(addNonNulls(list(matrix(1:10, ncol=2), matrix(11:20, ncol=2), 
                                  matrix(21:40, ncol=2), NULL, NULL)),
                 "Row dimensions of list elements are different")
})

test_that('vmessage returns NULL when verbose set to FALSE', {
    expect_equal(vmessage("Hello world", verbose=FALSE), NULL)
})

test_that('simulateDist fails when distr provided is not one of "unif", "norm", 
          "bin", "cat_norm" or "cat_unif" ', {
    expect_error(simulateDist(x=10, dist="gamma"),"Unknown distribution")
})

test_that('commaList2vector fails with unknown type', {
    expect_error(commaList2vector(type='matrix', commastring="1,2,2,3,4"), 
                 "Unknown type of comma-separated list elements")
})

test_that('probGen2expGen and expGen2probGen conversions work',{
    nrSamples <- 10
    genotype_prob <- rep(0, 3*nrSamples)
    genotype_prob[seq(1, nrSamples*3, 3) + sample(0:2, 10, replace=TRUE)] <- 1
    expect_equal(genotype_prob, expGen2probGen(probGen2expGen(genotype_prob)))
})

test_that('simulateDist fails with unknown distribution', {
    expect_error(simulateDist(dist="beta"), "Unknown distribution")
})

test_that('simulateDist fails with length of distribution argument', {
    dist=c("unif", "norm", "bin", "cat_norm", "cat_unif")
    expect_error(simulateDist(), 
                 paste("Please specify exactly one distribution to sample from,",
                 "currently ", length(dist), " provided.", sep=""))
})