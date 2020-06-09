### Assignment 3 Code ###

library(spuRs)
rm( list = ls( ) ) # this removes everything in the R workspace
# so that you can start with a clean slate

###########################################

days <- c(31,28,31,30,31,30,31,31,30,31,30,31)

rainfall <- read.csv(file = "D:/School/Mast of Sci/MAST90045/A3/Melbourne_rainfall_mm-1.csv")

seasame <- 1 # seed for pseudorandom number generator

set.seed(seasame)
###########################################
##### Task 1: Y_sim and F_Y Functions #####

Ysim <- function(p, m, la, n) {
  A <- rbinom(n, size = 1, prob = p)
  B <- rgamma(n, shape = m, rate = la)
  return(sum(A*B))
}

F_Y <- function(x, p, m, la, n) {
  df.vec <- c(rep(NA, n)) # initializing a vector of probabilities to be added
  pN0 <- dbinom(0, n, p)
  for (k in 1:n) { # summing Pr(C_k <= x)*Pr(N = k), each term will be in df.vec
    kay <- k # to pass on k as a parameter in the distribution functions
    df.vec[k] <- dbinom(kay, size = n, prob = p)*pgamma(x, shape = kay*m, rate = la)
  }
  return(pN0 + sum(df.vec))
}


ptest <- 0.27
mtest <- 0.24
latest <- 0.043
xtest <- 46
ntest <- 31
cdftest <- c(rep(NA, ntest))


##### Plotting F_Y Function #####

xvars <- c(seq(0, 400, 1))
dfY <- c(rep(NA, length(xvars)))

for (i in 1:length(xvars)) {
  dfY[i] <- F_Y(xvars[i], ptest, mtest, latest, ntest)
}


plot(xvars, dfY)



##### Using Ysim to check F_Y #####

# Simulating Ysim
yCheck <- c(rep(NA, 10^6))
for (i in 1:length(yCheck)) {
  yCheck[i] <- Ysim(ptest, mtest, latest, ntest)
}

# Defining the bin sizes
bins = seq(0,((max(yCheck)%/%20)*20+20),20)

# Plotting a histogram of Ysim
hist(yCheck, breaks = bins)

# Finding mean and SD for each bin

# To calculate expected value of each bin more easily
meanf_Y <- function(x) {
  fYa <- F_Y(x, ptest, mtest, latest, ntest)
  fYb <- F_Y(x + 20, ptest, mtest, latest, ntest)
  return(fYb - fYa)
}

# To calculate SD of each bin more easily
sdf_Y <- function(x) {
  fYa <- F_Y(x, ptest, mtest, latest, ntest)
  fYb <- F_Y(x + 20, ptest, mtest, latest, ntest)
  return((fYb - fYa)*(1 - fYb + fYa))
}

# Initializing vectors to store the mean and SD for each bin
yBins_mean <- c(rep(NA, length(bins)))
yBins_sd <- c(rep(NA, length(bins)))

# Calculating each bin's mean and SD with the functions defined above
for (i in 1:length(bins)) {
  yBins_mean[i] <- length(yCheck)*meanf_Y(bins[i])
  yBins_sd[i] <- sqrt(length(yCheck)*sdf_Y(bins[i]))
}


#############################################
##### Task 2: Estimating the Parameters #####

# Some quick estimates
y.hat <- rainfall$mean_rainfall
p.hat <- rainfall$mean_days_rain/days

# Defining the deciles (for legibility in defining the loss function)
d1 <- rainfall$decile_1_rainfall
d5 <- rainfall$median_rainfall
d9 <- rainfall$decile_9_rainfall

# Creating the loss function
L = function(x) {
  w1 <- 0.09
  w2 <- 0.25
  loss1 <- ((F_Y(d1[i], p.hat[i], (x*y.hat[i])/(days[i]*p.hat[i]), x, days[i]) - 0.1)^2)/w1
  loss5 <- ((F_Y(d5[i], p.hat[i], (x*y.hat[i])/(days[i]*p.hat[i]), x, days[i]) - 0.5)^2)/w2
  loss9 <- ((F_Y(d9[i], p.hat[i], (x*y.hat[i])/(days[i]*p.hat[i]), x, days[i]) - 0.9)^2)/w1
  return(loss1 + loss5 + loss9)
}

# Initializing vectors for lambda and m estimates
w <- length(rainfall$mean_days_rain)
la.hat <- c(rep(NA, w))
la.hatsd <- c(rep(NA, w))


# Minimizing the loss function to estimate lambda
for (i in 1:w) {
  la.hat[i] <- optim(latest, L)$par
  la.hatsd[i] <- optim(latest, L)$value
}

# Estimating m
m.hat <- (la.hat*y.hat)/(p.hat*days)


##########################################
##### Task 3: Simulating Water Usage #####

# Let's define a function to return the month depending on the day

month <- function(day) {
  if (day <= 31) {
    m <- 1
  } else if (day <= 59) {
    m <- 2
  } else if (day <= 90) {
    m <- 3
  } else if (day <= 120) {
    m <- 4
  } else if (day <= 151) {
    m <- 5
  } else if (day <= 181) {
    m <- 6
  } else if (day <= 212) {
    m <- 7
  } else if (day <= 243) {
    m <- 8
  } else if (day <= 273) {
    m <- 9
  } else if (day <= 304) {
    m <- 10
  } else if (day <= 334) {
    m <- 11
  } else m <- 12
  return(m)
  }





tanksim = function(num_years, maxraintank, maxgreytank, phat, mhat, lahat, plotflag = F) {
  # Constants
  days <- 1:365
  days_in_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  mean_max_temp <- c(25.9, 25.8, 23.9, 20.3, 16.7, 14.1, 13.5, 15.0, 17.2, 19.7, 22.0, 24.2)
  roofarea <- 100 # m^2
  gardenarea <- 200 # m^2
  flushsize <- 5 # litres per flush
  numflush <- function() rbinom(1, 15, .8) # flushes per day for four people, at home half the day
  showersize <- 35 # in litres
  numshower <- 4 # showers per day for four people
  washsize <- 35 # litres per load
  numwash <- function() rbinom(1, 8, .125) # washes per day for four people
  
  Xsim <- function(month) {
    # simulate rainfall for a day in given month
    rbinom(1, 1, phat[month])*rgamma(1, mhat[month], lahat[month])
  }
  
  # initializing total data
  yearly_saved <- c(rep(NA, length(num_years)))
  avg_raintank <- c(rep(0, 365))
  avg_greytank <- c(rep(0, 365))
  yearly_raintank <- matrix(NA, nrow = 365, ncol = num_years)
  yearly_greytank <- matrix(NA, nrow = 365, ncol = num_years)

  
  for (year in 1:num_years) {
    # initializing yearly vectors; each year these vectors reset
    raintank <- c(rep(0, 365))
    greytank <- c(rep(0, 365))
    daily_saved <- c(rep(0, 365))
    gardenwater <- c(rep(0, 365))
    for (i in 1:365) {
      M <- month(i) # determining the month
      rain <- Xsim(M) # simulating daily rain based on month
      flushwater <- numflush() * flushsize
      
      #updating tank sizes and adjusting for the maximum sizes
      raintank[i] <-  raintank[i] + (roofarea * rain) # updating the rainwater tank based on the rain
      greytank[i] <- greytank[i] + (numshower * showersize) + (numwash() * washsize)  # updating the greywater tank
      if (raintank[i] > maxraintank) {
        raintank[i] <- maxraintank
      }
      if (greytank[i] > maxgreytank) {
        greytank[i] <- maxgreytank
      }
      
      # defining garden variables
      mmT <- mean_max_temp[M]/15 # finding the required garden water
      gardenwater[i] <- rain*gardenarea # how much the garden was watered by the rain
      if (i < 3) {
        avggardenwater <- sum(gardenwater[1:i])/i
      } else {
        avggardenwater <- mean(gardenwater[(i-2):i])
      }
      
      if (flushwater <= greytank[i]) {
        greytank[i] <- greytank[i] - flushwater
      } else if (flushwater <= (raintank[i] + greytank[i])) {
        raintank[i] <- raintank[i] + greytank[i] - flushwater
        greytank[i] <- 0
      } else {
        greytank[i] <- 0
        raintank[i] <- 0
      }
      
      # finding the required amount of garden water
      if (avggardenwater < mmT) {
        gardenreq <- mmT - avggardenwater
      } else {
        gardenreq <- 0
      }
      
      if (gardenreq*gardenarea <= raintank[i]) {
        raintank[i] <- raintank[i] - gardenarea*gardenreq
      } else if (gardenreq*gardenarea <= (raintank[i] + greytank[i])) {
        greytank[i] <- greytank[i] + raintank[i] - gardenreq*gardenarea
        raintank[i] <- 0
      } else {
        greytank[i] <- 0
        raintank[i] <- 0
      }
      daily_saved[i] <- raintank[i] + greytank[i]
      # updating next day's tank levels
      if (i < 365){
        raintank[i+1] <- raintank[i]
        greytank[i+1] <- greytank[i]
      }
    }
    # Summing up the water saved this year
    yearly_saved[year] <- sum(daily_saved)
    yearly_greytank[,year] <- greytank
    yearly_raintank[,year] <- raintank
  }
  # recording average grey/rain levels
  for (n in days) {
    avg_greytank[n] <- mean(yearly_greytank[n,])
    avg_raintank[n] <- mean(yearly_raintank[n,])
  }
  
  # plot flag
  if (plotflag == T) {
    plot(x = days, y = avg_raintank,
         type = "l", col = "blue", xlab = "Days", ylab = "Litres",
         lwd = 2)
    lines(x = days, y = avg_greytank,
          type = "l", col = "red", lwd = 2 )
    if (num_years > 1) {
      for (y in 1:num_years){
        lines(x = days, y = yearly_raintank[,y],
              type = "l", lty = 3, lwd = 0.75,
              col = "lightblue")
        lines(x = days, y = yearly_greytank[,y],
              type = "l", lty = 3, lwd = 0.75,
              col = "lightpink")
      }
    }
  }
  return(yearly_saved)
}


tanksim(5, 3000, 1000, p.hat, m.hat, la.hat, T)