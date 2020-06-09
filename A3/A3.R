
rm(list=ls())

# Part 1


dm = c(31,28,31,30,31,30,31,31,30,31,30,31) # days in a month

# testing distributions in R

set.seed(13)

rbinom(n = 31, size = 1, prob = 0.27)*rgamma(31,0.24,0.043)


# Ysim function

Ysim = function(p, m, la, n) {
  
  A = rbinom(n = n, size = 1, prob = p)
  
  B = rgamma(n = n, shape = m, rate = la)
  
  return(sum(A*B))

}



# F_Y function - we need to use dbinom, NOT pbinom!

F_Y = function(x, p, m, la, n) {
  my_list = c()
  for(k in 1:n) {my_list = 
    append(my_list, 
           pgamma(x, shape = k*m, rate = la)*dbinom(k, size = n, p))}
  
  cd_Y = dbinom(0, n, p) + sum(my_list)
  return(cd_Y)
}


# CDF of F_Y

xs = 1:500

fYplt = c()
set.seed(13)
for(a in 1:length(xs)) fYplt = append(fYplt, F_Y(a, 0.27, 0.24, 0.043, 31))

plot(xs, fYplt)


# Import rain data

data = read.csv(file = "D:/School/Mast of Sci/MAST90045/A3/Melbourne_rainfall_mm-1.csv")


# Using Ysim to check F_Y is correct
set.seed(13)
Ysim(0.27,0.24,0.043,31)

# checking this will work...
lh = c(rep(NA,10^4))
for(k in 1:10^4) lh[k] = Ysim(0.27,0.24,0.043,31)

max(lh)


# making our bins
binsize = c(seq(0,360,20))

Bins = findInterval(lh, binsize)
Bins2 = as.data.frame(table(Bins))

Bins2


mod_fY = function(x) return(F_Y(x, 0.27,0.24,0.043,31))

set.seed(13)
(mod_fY(20) - mod_fY(0))*10^4

#Standard deviation
sqrt((mod_fY(20) - mod_fY(0))*10^4*(1- mod_fY(20) + mod_fY(0)))

# Generate the frequencies
Fmean = c()
Fsd = c()
for(i in c(seq(0,340,20))) {
  Fmean = append(Fmean, (mod_fY(i+20) - mod_fY(i))*10^4)
  Fsd = append(Fsd, sqrt((mod_fY(i + 20) - 
                            mod_fY(i))*10^4*(1- mod_fY(i + 20) + mod_fY(i))))
}

freq = Bins2$Freq
FRE = c(freq[1:14],0,0,freq[15:16])
FRE

#Getting everything into a nice dataframe for plotting
BINS = data.frame('Bins' = seq(20,360,20), 'Freq' = FRE, 'Mean' = Fmean, 'SD' = Fsd )
BINS




library(ggplot2)

H = ggplot(BINS, aes(factor(Bins), Freq, group=Bins)) + geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymin = Freq-SD, ymax = Freq+SD), width =.3, colour = 'red') +
  geom_point(aes(y=Mean), colour = 'blue')
  

print(H)


# Part 2

#defining parameters
p_est = c()
x_est = c()
dec_1 = c()
dec_5 = c()
dec_9 = c()
y_est = c()
for(i in 1:12) {
  p_est = append(p_est, data$mean_days_rain[i]/dm[i])
  x_est = append(x_est, data$mean_rainfall[i]/dm[i])
  dec_1 = append(dec_1, data$decile_1_rainfall[i])
  dec_5 = append(dec_5, data$median_rainfall[i])
  dec_9 = append(dec_9, data$decile_9_rainfall[i])
  y_est = append(y_est, data$mean_rainfall[i])
}

# lambda function test?

L = function(lambda) {
  jingo = (100/9)*(F_Y(x = dec_1[1], 
                       p = p_est[1], 
                       m =( y_est[1]*lambda)/(dm[1]*p_est[1]), 
                       la = lambda, n = dm[1] )-0.1)^2 + 
    4*(F_Y(x = dec_5[1], 
           p = p_est[1], 
           m =( y_est[1]*lambda)/(dm[1]*p_est[1]), 
           la = lambda, n = dm[1] )- 0.5)^2 + 
    (100/9)*(F_Y(x = dec_9[1], 
                 p = p_est[1], 
                 m =( y_est[1]*lambda)/(dm[1]*p_est[1]), 
                 la = lambda, 
                 n = dm[1] )-0.9)^2
  return(jingo)
}

#testing optimsation
optim(par = 0.5, fn = L)$par

# actual lambda function

L = function(lambda) {
  jingo = (100/9)*(F_Y(x = dec_1[i], 
                       p = p_est[i], 
                       m =( y_est[i]*lambda)/(dm[i]*p_est[i]), 
                       la = lambda, 
                       n = dm[i] )-0.1)^2 + 
    4*(F_Y(x = dec_5[i], 
           p = p_est[i], 
           m =( y_est[i]*lambda)/(dm[i]*p_est[i]), 
           la = lambda, n = dm[i] )- 0.5)^2 + 
    (100/9)*(F_Y(x = dec_9[i], 
                 p = p_est[i], 
                 m =( y_est[i]*lambda)/(dm[i]*p_est[i]), 
                 la = lambda, 
                 n = dm[i] )-0.9)^2
  return(jingo)
}

lam_est = c()
lam_sd = c(rep(NA,12))
for(i in 1:12) {
  lam_est = append(lam_est, optim(par = 0.5, fn = L)$par)
  lam_sd[i] = optim(par = 0.5, fn = L)$value
}

m_est = c(rep(NA,12))
for(i in 1:12) m_est[i] = (y_est[i]*lam_est[i])/(dm[i]*p_est[i])

print(lam_est)
print(m_est)
print(p_est)

library(stats)
sd(lam_est)
sd(p_est)
sd(m_est)

# Part 3 


month = function(day){
  dys = c(31,59,90,120,151,181,212,243,273,304,334,365)
  m = c(1:12)
  for (i in 1:length(dys)) {
    if (day <= dys[i]) {
      return(m[i])
    }
  }
}

tanksim = function(num_years, maxraintank, maxgreytank, 
                   phat = p_est, mhat = m_est, lahat = lam_est, plotflag = F) {
  # Constants
  days = 1:365
  days_in_month = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  mean_max_temp = c(25.9, 25.8, 23.9, 20.3, 16.7, 14.1, 13.5, 15.0, 
                    17.2, 19.7, 22.0, 24.2)
  roofarea = 100 # m^2
  gardenarea = 200 # m^2
  flushsize = 5 # litres per flush
  numflush = function() rbinom(1, 15, .8) # flushes per day 
                                          #for four people, at home half the day
  showersize = 35 # in litres
  numshower = 4 # showers per day for four people
  washsize = 35 # litres per load
  numwash = function() rbinom(1, 8, .125) # washes per day for four people
  
  Xsim = function(month) {
    # simulate rainfall for a day in given month
    rbinom(1, 1, phat[month])*rgamma(1, mhat[month], lahat[month])
  }
  
  save_year = c(rep(NA, length(num_years)))
  tank_rain_av = c(rep(0, 365))
  grey_rain_av = c(rep(0, 365))
  yearly_raintank = matrix(NA, nrow = 365, ncol = num_years)
  yearly_greytank = matrix(NA, nrow = 365, ncol = num_years)
  
  
  for (year in 1:num_years) {
    raintank = c(rep(0, 365))
    greytank = c(rep(0, 365))
    daily_saved = c(rep(0, 365))
    gardenwater = c(rep(0, 365))
    for (i in 1:365) {
      M = month(i) 
      rain = Xsim(M) 
      flushwater = numflush() * flushsize
      
      raintank[i] =  raintank[i] + (roofarea * rain)
      greytank[i] = greytank[i] + (numshower * showersize) +
        (numwash() * washsize)
      if (raintank[i] > maxraintank) raintank[i] = maxraintank
      if (greytank[i] > maxgreytank) greytank[i] = maxgreytank
      
      avMAXtemp = mean_max_temp[M]/15
      gardenwater[i] = rain*gardenarea
      if (i < 3) avggardenwater = sum(gardenwater[1:i])/i
      else avggardenwater = mean(gardenwater[(i-2):i])
      
      if (flushwater <= greytank[i]) greytank[i] = greytank[i] - flushwater
      else if (flushwater <= (raintank[i] + greytank[i])) {
        raintank[i] = raintank[i] + greytank[i] - flushwater
        greytank[i] = 0
      } 
      else {
        greytank[i] = 0
        raintank[i] = 0
      }
      
      if (avggardenwater < avMAXtemp) gardenreq = avMAXtemp - avggardenwater 
      else gardenreq = 0
      
      if (gardenreq*gardenarea <= raintank[i]) {raintank[i] = raintank[i] - 
        gardenarea*gardenreq}
      else if (gardenreq*gardenarea <= (raintank[i] + greytank[i])) {
        greytank[i] = greytank[i] + raintank[i] - gardenreq*gardenarea
        raintank[i] = 0
      } 
      else {
        greytank[i] = 0
        raintank[i] = 0
      }
      daily_saved[i] = raintank[i] + greytank[i]
      if (i < 365){
        raintank[i+1] = raintank[i]
        greytank[i+1] = greytank[i]
      }
    }
    save_year[year] = sum(daily_saved)
    yearly_greytank[,year] = greytank
    yearly_raintank[,year] = raintank
  }
  total_grey = c()
  total_rain = c()
  for (n in 1:num_years) {
    
    total_grey = c(total_grey,yearly_greytank[,n])
    total_rain = c(total_rain, yearly_raintank[,n])
    
  }
  
  if (plotflag == T) {

    plot(x = 1:length(total_rain), y = total_rain,
         type = "l", col = "blue", xlab = "Days", ylab = "Litres",
         lwd = 0.5)
    lines(x = 1:length(total_grey), y = total_grey,
          type = "l", col = "red", lwd = 0.5 )
    }
  total_save = sum(save_year)
  return(total_save)
}


tanksim(5, 3000, 1000, p_est, m_est, lam_est, T)


#prelim exploration of most water saved for different tank sizes for a 5 year period

I = c(seq(2000,10000,500))
J = c(seq(500,8500,500))
B = 0
for (i in I) {
  for (j in J) {
    A = tanksim(5,i,j,p_est, m_est,lam_est,F)
    if (A > B) {
      B = A
      best_i = i
      best_j = j
    }
  }
}
cat('Best Rainwater Tank Size: ', best_i, '\n Best Greywater Tank Size: ', best_j,
    '\n Water Saved: ', B)
tanksim(5,10000,8000,plotflag = T)

tanksim(10,10000,8000,plotflag = T)


I = c(seq(2000,50000,500))
J = c(seq(500,48500,500))
B = 0
for (i in I) {
  for (j in J) {
    A = tanksim(10,i,j)
    if (A > B) {
      B = A
      best_i = i
      best_j = j
    }
  }
}
cat('Best Rainwater Tank Size: ', best_i, '\n Best Greywater Tank Size: ', best_j,
    '\n Water Saved: ', B)

tanksim(10,best_i, best_j, plotflag = T)
