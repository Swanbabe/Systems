##### ---------- Assignment 1 ---------- #####
rm(list=ls()) #remove objects
library(spuRs)

### ---------- Task 1: Compute Volume ---------- ###

# Defining volume function that uses Simpson's method

#source('spuRs/resources/scripts/simpson.r')

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
    ep <- 1e-15
    V <- simpson(ftn, 0, h, tol = ep)
    return(V)
  }
}


### ---------- Task 2: Compute Height from Volume ---------- ###

# Defining height function using the bisection method as the root-finding algorithm

# if using bisection method:
#source('spuRs/resources/scripts/bisection.r')
# if using Newton-Raphson method:
# program spuRs/resources/scripts/newtonraphson.r
# loadable spuRs function

newtonraphson <- function(ftn, x0, tol = 1e-9, max.iter = 100) {
  # Newton_Raphson algorithm for solving ftn(x)[1] == 0
  # we assume that ftn is a function of a single variable that returns
  # the function value and the first derivative as a vector of length 2
  #
  # x0 is the initial guess at the root
  # the algorithm terminates when the function value is within distance
  # tol of 0, or the number of iterations exceeds max.iter
  
  # initialise
  x <- x0
  fx <- ftn(x)
  iter <-  0
  
  # continue iterating until stopping conditions are met
  while ((abs(fx[1]) > tol) && (iter < max.iter)) {
    x <- x - fx[1]/fx[2]
    fx <- ftn(x)
    iter <- iter + 1
    cat("At iteration", iter, "value of x is:", x, "\n")
  }
  
  # output depends on success of algorithm
  if (abs(fx[1]) > tol) {
    cat("Algorithm failed to converge\n")
    return(NULL)
  } else {
    cat("Algorithm converged\n")
    return(x)
  }
}


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
      return(volume(x, m, y) - w)
    }
    zeta <- function(x){
      z <- c(eta(x), ftn(x))
      return(z)
    }
    u <- newtonraphson(zeta, hmax, tol = 1e-9, max.iter = 100)
  }
  return(u)
}

### ---------- Task 3: Testing Height function ---------- ###

# Given equations for area A(h) = h,
# and H(h,v) = sqrt(h^2 + 2v)

A <- function(x) { return(x)}

uformula <- function(h, v) { return(sqrt(h^2 + 2*v))}

# Initializing a vector for the correct test-heights

# Taking test data for test-height and test-volume from the assignment
htest <- c(-1, .5, .5, 3.5, 3.5, 3.5, 5)
vtest <- c(1, 1, -1, 1, 2, -1, 1)

# Finding the length of our test vectors
m <- length(htest)

# Initializing vectors to calculate the test heights
# And test our function
ucheck <- c(rep(0,m))
utest <- c(rep(0, m))

# Calculating test heights via formula given
for (i in 1:m ) {
  ucheck[i] <- suppressWarnings(uformula(htest[i], vtest[i]))
}

# Testing our function with the same data
for (i in 1:m) {
  utest[i] <- suppressWarnings(height(htest[i], hmax = 4, vtest[i], A))
}
show(utest)

### ---------- Task 4: Tracking the height of the dam ---------- ###
# Reading the data
catch_path = 'C:/Users/amalt/Documents/UniMelb Semester 1 (2020)/4 - Systems Modelling and Simulation/Assignment/catchment_b-1.txt'

catch_data <- scan(file = catch_path)

n <- length(catch_data)



# Calculating h(t+1)

lvl <- c(1, rep(0,n))
a <- 2    # alpha is our usage volume
b <- 0.1    # beta is our evaporation coefficient
m <- 3    # max height

for (i in 1:n){
  if (lvl[i] <= 2){
    A <- function(x) {return(100*(x^2))}
  } else {
    A <- function(x) {return(400*(x-1))}
  }
  v <- catch_data[i] - a - b*A(lvl[i])
  lvl[i+1] <- height(lvl[i], m, v, A)
}
