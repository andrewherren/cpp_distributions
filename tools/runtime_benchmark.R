library(microbenchmark)

n <- 100000

# Case 1
lambda <- 1.2
a <- 1.7
b <- 1.5
microbenchmark(x <- sample_gig(n, lambda, a, b), 
               y <- rgig(n, lambda, b, a))

# Case 2
lambda <- 0.5
a <- 0.5
b <- 0.75
microbenchmark(x <- sample_gig(n, lambda, a, b), 
               y <- rgig(n, lambda, b, a))

# Case 3
lambda <- 0.5
a <- 0.5
b <- 0.25
microbenchmark(x <- sample_gig(n, lambda, a, b), 
               y <- rgig(n, lambda, b, a))
