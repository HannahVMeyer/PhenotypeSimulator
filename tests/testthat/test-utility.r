context('Test utility functions')

test_that("vmessage fails when non-character specified as separator", {
    expect_error(vmessage(c("hello", "world"), sep=9))
})

test_that("vmessage returns correct messsage", {
    expect_that(vmessage(c("Hello", "World"), sep="Little"), 
                shows_message("HelloLittleWorld"))
})

test_that('vmessage returns NULL when verbose set to FALSE', {
    expect_equal(vmessage("Hello world", verbose=FALSE), NULL)
})

test_that('commaList2vector returns vector from comma-separated inputlist', {
    expect_equal(commaList2vector("1,2,3"), c(1,2,3))
})

test_that('commaList2vector returns NULL', {
    expect_equal(commaList2vector(), NULL)
})

test_that('commaList2vector fails with unknown type', {
    expect_error(commaList2vector(type='matrix', commastring="1,2,2,3,4"), 
                 "Unknown type of comma-separated list elements")
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

test_that('addNonNulls fails when non-list supplied as input', {
              expect_error(addNonNulls(matrix(1:10, ncol=2)),
                           "addNonNulls expects input of type list")
})

test_that('simulateDist fails when distr provided is not one of "unif", "norm", 
          "bin", "cat_norm" or "cat_unif" ', {
    expect_error(simulateDist(x=10, dist="gamma"),"Unknown distribution")
})

test_that('simulateDist fails when number of observations is less than 0', {
    expect_error(simulateDist(x=-10, dist="bin"),
                 "The number of observations to simulate")
          })

test_that('simulateDist fails when standard deviation for normal distribution
          less than 0', {
    expect_error(simulateDist(x=10, dist="norm", std=-0.3),
                 "Simulating normal distribution")
})

test_that('simulateDist fails when probability for binomial dist is not
          provided', {
              expect_error(simulateDist(x=10, dist="bin"),
                           "Simulating binomial distribution")
          })

test_that('simulateDist fails when probability for binomial dist is less
          than zero', {
              expect_error(simulateDist(x=10, dist="bin", prob=-0.02),
                           "Simulating binomial distribution")
          })

test_that('simulateDist fails with unknown distribution', {
    expect_error(simulateDist(x=10, dist="beta"), "Unknown distribution")
})

test_that('simulateDist fails with length of distribution argument', {
    dist=c("unif", "norm", "bin", "cat_norm", "cat_unif")
    expect_error(simulateDist(x=10), 
                 paste("Please specify exactly one distribution to sample from,",
                       "currently ", length(dist), " provided.", sep=""))
})

test_that('simulateDist fails when number of categories for categorical dist
          is not provided', {
              expect_error(simulateDist(x=10, dist="cat_unif"),
                           "Simulating categorical distribution")
          })

test_that('simulateDist fails when  number of categories for categorical dist
          is  less than zero', {
              expect_error(simulateDist(x=10, dist="cat_unif", categories=-2),
                           "Simulating categorical distribution")
          })

test_that("probGen2expGen fails when input is not a numeric vector", {
    expect_error( probGen2expGen(list(c(0.1, 0.9, 0))), 
                  "probGen2expGen takes")
})

test_that("probGen2expGen fails when input vector is not a multiple of three", {
    expect_error( probGen2expGen(c(0.1, 0.9)), 
                  "Length of genotype probabilty vector")
})

test_that("probGen2expGen fails with NA", {
    expect_error( probGen2expGen(c(0.1, 0.8, NA)), 
                  "Samples need to be fully genotyped")
})

test_that("probGen2expGen fails with incorrect probabilities", {
    expect_error( probGen2expGen(c(0.1, 0.8, 0)), 
                  "Genotype probabilities do not sum to 1")
})

test_that('expGen2probGen fails when wrong encoding is given', {
    expect_error(expGen2probGen(c(0,2,3)),
                 "Genotypes can only be encoded as 0,1,2")
})

test_that('expGen2probGen fails when wrong input format is given', {
    expect_error(expGen2probGen(list(c(0,2,1))),
                 "expGen2probGen takes a vector of ")
})

test_that('expGen2probGen returns triple NA', {
    expect_equal(expGen2probGen(c(NA)), c(NA,NA,NA))
})

test_that('expGen2probGen returns expected output', {
    expect_equal(expGen2probGen(c(2)), c(0,0,1))
}) 
test_that('probGen2expGen and expGen2probGen conversions work',{
    nrSamples <- 10
    genotype_prob <- rep(0, 3*nrSamples)
    genotype_prob[seq(1, nrSamples*3, 3) + sample(0:2, 10, replace=TRUE)] <- 1
    expect_equal(genotype_prob, expGen2probGen(probGen2expGen(genotype_prob)))
})

