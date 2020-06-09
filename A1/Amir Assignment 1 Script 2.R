##### ---------- Assignment 1 ---------- #####
rm(list=ls()) #remove objects
library(spuRs)

### ---------- Task 1: Compute Volume ---------- ###

# Defining volume function that uses Simpson's method

simpson <- function(ftn, a, b, tol = 1e-8, verbose = FALSE) {
  # numerical integral of ftn from a to b
  # using Simpson's rule with tolerance tol
  #
  # ftn is a function of a single variable and a < b
  # if verbose is TRUE then n is printed to the screen
  
  # initialise
  n <- 4
  h <- (b - a)/4
  fx <- sapply(seq(a, b, by = h), ftn)
  S <- sum(fx*c(1, 4, 2, 4, 1))*h/3
  S.diff <- tol + 1  # ensures we loop at least once
  
  # increase n until S changes by less than tol
  while (S.diff > tol) {
    # cat('n =', n, 'S =', S, '\n')  # diagnostic
    S.old <- S
    n <- 2*n
    h <- h/2
    fx[seq(1, n+1, by = 2)] <- fx  # reuse old ftn values
    fx[seq(2, n, by = 2)] <- sapply(seq(a+h, b-h, by = 2*h), ftn)
    S <- h/3*(fx[1] + fx[n+1] + 4*sum(fx[seq(2, n, by = 2)]) +
                2*sum(fx[seq(3, n-1, by = 2)]))
    S.diff <- abs(S - S.old)
  }
  if (verbose) cat('partition size', n, '\n')
  return(S)
}


volume <- function(h, hmax, ftn) {
  # h is assumed to be the water level of the dam
  # hmax is the height of the dam
  # ftn is the area, a single-variate function of h; A by default for later use
  # ep will be our default tolerance, so that we don't have to use
  # the default tolerance of 1e-8 for simpson.r
  
  # if h is out of bounds, the function should return 'NA'
  if (h < 0) { 
    return(NA)
  } else if(h > hmax) {
    return(NA)
  } else {
    # force simpson to use a tolerance of 1e-15
    # this will allow the volume function to use its own tolerance when nested inside another function
    # with a possibly different tolerance
    ep <- 1e-9
    V <- simpson(ftn, 0, h, tol = ep)
    return(V)
  }
}


### ---------- Task 2: Compute Height from Volume ---------- ###

# Defining height function using the bisection method as the root-finding algorithm

# if using bisection method:
#source('spuRs/resources/scripts/bisection.r')
# if using Newton-Raphson method:
source('spuRs/resources/scripts/newtonraphson.r')

height <- function(h, hmax, v, ftn) {
  # defining vmax and w in terms of our inputs: h, hmax, and v
  vmax <- volume(hmax, hmax, ftn)
  if (h > hmax) {
    h <- hmax
  } else if (h < 0) {
    h <- 0
  }
  w <- as.numeric(volume(h, hmax, ftn) + v)
  if (w > vmax){
    u <- hmax
  } else if (w < 0){
    u <- 0
  } else {
    # defining a single-variate function for Newton-Raphson, that returns the vector (f(u), f'(u))
    # note that this function calls on the volume function,
    # and depends on the 'hmax' and 'ftn' arguments of height(...)
    eta <- function(x){
      m <- hmax
      y <- ftn
      f <- volume(x, m, y) - w
      d <- ftn(x)
      return(c(f,d))
    }
    u <- newtonraphson(eta, hmax, tol = 1e-9, max.iter = 100)
  }
  return(u)
}

### ---------- Task 3: Testing Height function ---------- ###
for (i in 1:examplem) {
  utest[i] <- suppressWarnings(height(exampleheight[i], hmax = 4, examplevols[i], A))
  show(utest)
}
# Given equations for area A(h) = h,
# and H(h,v) = sqrt(h^2 + 2v)
# where hmax = 4

hmax <- 4
A <- function(x) { return(x)}

# should simulate the height function, but using a precise calculation
# in place of a root-finding algorithm
exampleformula <- function(h, v) {
  # First we condition that our value for h makes sense,
  # since the water level will always be in [0, hmax],
  # even if the volume function isn't defined outside
  # this interval
  
  # We therefore 'make sense' of the input, h
  # by re-interpreting it
  if (h < 0) {
    h <- 0
  } else if (h > hmax){
    h <- hmax
  }
  # now that we've made sure our h is sensible
  # we wish to calculate the new height, u
  # using the precise formula given
  u2 <- h^2 + 2*v
  
  # since we want to return a real and practical
  # value for u,
  # we check that u^2 is inside our interval
  if (u2 < 0){
    u2 <- 0
  } else if (u2 > hmax^2){
    u2 <- hmax^2
  }
  # thus we have ensured that our interpretation
  # of the new height, u, makes sense
  # i.e., is inside our interval
  u <- sqrt(u2)
  return( u )
  }

# Initializing a vector for the correct test-heights

# Taking test data for test-height and test-volume from the assignment
exampleheight <- c(-1, .5, .5, 3.5, 3.5, 3.5, 5)
examplevols <- c(1, 1, -1, 1, 2, -1, 1)

# Finding the length of our test vectors
examplem <- length(exampleheight)

# Initializing vectors to calculate the test heights
# And test our function
ucheck <- c(rep(0, examplem))
utest <- c(rep(0, examplem))

# Calculating test heights via formula given
for (i in 1:examplem ) {
  ucheck[i] <- suppressWarnings(exampleformula(exampleheight[i], examplevols[i]))
  show(ucheck)
  }

# Testing our function with the same data
for (i in 1:examplem) {
  utest[i] <- suppressWarnings(height(exampleheight[i], hmax = 4, examplevols[i], A))
  show(utest)
  }

# A successful test suggests that our height function is ready
# to be used on our data

### ---------- Task 4: Tracking the height of the dam ---------- ###
# Reading the data
catch_path = 'C:/Users/amalt/Documents/UniMelb Semester 1 (2020)/4 - Systems Modelling and Simulation/Assignment/catchment_b-1.txt'

catch_data <- scan(file = catch_path)

# defining the length of our data vector
n <- length(catch_data)

# setting our parameters

a <- 2      # alpha is our usage volume
b <- 0.1    # beta is our evaporation coefficient
m <- 3      # max height

# initializing our height vector
lvl <- c(1, rep(0,n))

# we're now ready to calculate h[t+1]
# for t = 1, ..., n
# we will do this by analogy from the test above

# We have two possible area functions
A <- function(x) {
  if (x <= 2) { area = 100*(x^2)
  } else {
    area = 400*(x-1)
  }
  return(area)
}

for (t in 1:n){
  # first we want to make sure to interpret the level correctly
  # by ensuring that it is in our interval [0, hmax]
  # we don't actually need these conditionals since
  # this is ensured from within the height function
  if (lvl[t] < 0) {
    lvl[t] <- 0
  } else if (lvl[t] > m){
    lvl[t] <- m
  }
  # now we invoke our glorious height function
  # which calls on our wonderful volume function
  # and uses a root-finding algorithm
  v <- (catch_data[t] - a - b*A(lvl[t]))
  u <- height(lvl[t], m, v, A)
  
  # lastly, we again want to interpret
  # the new level correctly, meaning we want
  # to ensure its value is within our interval
  # this conditional should also be unnecessary
  # if the height function works properly
  if (u < 0){
    u <- 0
  } else if (u > m){
    u <- m
  }
  lvl[t+1] <- u
}







# Testing why NR gets hung:

volume(lvl[40],3,A.test)
w = catch_data[40] - 2 - 0.1*A.test(lvl[40])
eta.test <- function(x){
  Z = volume(x, 3, A) - w
  return( Z )
}
zeta.test <- function(x){
  L <- c(eta.test(x), A(x))
  return( L )
}
newtonraphson(zeta.test, 3)
