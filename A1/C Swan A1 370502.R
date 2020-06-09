rm(list=ls())

#spuRs alogorithms

# Simpson's Rule, taken from spuRs:
# program spuRs/resources/scripts/simpson.r

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


######################################################################################################

#volume
#we have function 'ftn' as the cross-sectional area of at height h 

volume = function(h,hmax,ftn) {
  if (h < 0) {
    return(NA)
  } else if (h > hmax) {
    return(NA)
  } else {
    V = simpson(ftn,0,h, tol = 1e-15)
    return(V)
  }
}

#height


height = function(h,hmax,v,ftn) {
  if (h < 0) {
    return(NA)
  } else if (h > hmax) {
    return(NA)
  } else {
    new_volume = as.numeric(volume(h, hmax, ftn) + v)
    max_volume = volume(hmax,hmax,ftn)
    if (new_volume > max_volume) {
      dam_level = hmax
      print(cat('The dam has reached the maximum level of', dam_level, '\n'))
    } else if (new_volume < 0) {
      dam_level = 0
      print(cat('The dam has reached the minimum level of', dam_level, '\n'))
    } else { 
      newt_input = function(x) {
        N = volume(x,hmax,ftn) - new_volume
        return(c(N,ftn(x)))
      }
    dam_level = newtonraphson(newt_input,hmax, tol = 1e-15, max.iter = 1000)
    }
  return(dam_level)
  }
}


#translating our variables to the assignment:

# dam_level = u = H(h,v)
# volume (function) = V(h)
# new_volume = V(u)
# ftn (function) = A(h)

######################################################################################################

#test case




heights = c(-1,0.5,0.5,3.5,3.5,3.5,5)
volumes = c(1,1,-1,1,2,-1,1)

# Let's evaluate our test criteria:

A = function(h) return(h)
V_h = function(k) return(k^2/2)
H_h_v = function(i,j) return(sqrt(i^2 + 2*j))

for (i in seq(1,7)) {
      print(cat("The Area of A", i," is ", A(heights[i]), "\n", 
                'The Volume of V', i," is ", V_h(heights[i]), "\n",
                'The Dam Level of H',i, " is ", H_h_v(heights[i],volumes[i]), "\n",
                'The New Height of h', i, " is ", height(heights[i], hmax = 4, volumes[i], ftn = A), "\n"))
}

######################################################################################################


#tracking height over time

catchment = scan(file = "/School/Mast of Sci/MAST90045/A1/catchment_b-1.txt")

dam_calc = function(volume_data, h_vector = c(1), hmax = 3, alpha = 2, beta = 0.1) {
  len = length(catchment)
  Area = function(h) {
    if (0 <= h & h <= 2) {
      return(100*h^2)
    } else if (2 <= h & h <= 3){
      return(400*(h - 1))
    } else (
      return(NA)
    )
  }
  for (i in seq(1,len)) {
    if (h_vector[i] < 0) {
      h_vector[i] = 0
    } else if (h_vector[i] > hmax) {
      h_vector[i] = hmax
    } 
    vol = volume_data[i] - alpha - beta*Area(h_vector[i])  
    next_h = height(h_vector[i], hmax, vol, Area)
       
    h_vector = append(h_vector, next_h)
    
  } 
  return(h_vector)
}

rain_volumes = append(catchment, NA)
dam_heights = dam_calc(catchment)
Results = cbind(rain_volumes,dam_heights)
write.csv(Results, file = "/School/Mast of Sci/MAST90045/A1/results.csv")



