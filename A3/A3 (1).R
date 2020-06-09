
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
  for(k in 1:n) {my_list = append(my_list, pgamma(x, shape = k*m, rate = la)*dbinom(k, size = n, p))}
  
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
  Fsd = append(Fsd, sqrt((mod_fY(i + 20) - mod_fY(i))*10^4*(1- mod_fY(i + 20) + mod_fY(i))))
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
  jingo = (100/9)*(F_Y(x = dec_1[1], p = p_est[1], m =( y_est[1]*lambda)/(dm[1]*p_est[1]), la = lambda, n = dm[1] )-0.1)^2 + 
    4*(F_Y(x = dec_5[1], p = p_est[1], m =( y_est[1]*lambda)/(dm[1]*p_est[1]), la = lambda, n = dm[1] )- 0.5)^2 + 
    (100/9)*(F_Y(x = dec_9[1], p = p_est[1], m =( y_est[1]*lambda)/(dm[1]*p_est[1]), la = lambda, n = dm[1] )-0.9)^2
  return(jingo)
}

#testing optimsation
optim(par = 0.5, fn = L)$par

# actual lambda function

L = function(lambda) {
  jingo = (100/9)*(F_Y(x = dec_1[i], p = p_est[i], m =( y_est[i]*lambda)/(dm[i]*p_est[i]), la = lambda, n = dm[i] )-0.1)^2 + 
    4*(F_Y(x = dec_5[i], p = p_est[i], m =( y_est[i]*lambda)/(dm[i]*p_est[i]), la = lambda, n = dm[i] )- 0.5)^2 + 
    (100/9)*(F_Y(x = dec_9[i], p = p_est[i], m =( y_est[i]*lambda)/(dm[i]*p_est[i]), la = lambda, n = dm[i] )-0.9)^2
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
print(x)


# Part 3 

tanksim = function(num_years, maxraintank, maxgreytank, phat, mhat, lahat, plotflag = F) {

  # Constants
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

  rainwatertank = 0
  greywatertank = 0
  
  Xsim <- function(month) {
    # simulate rainfall for a day in given month
    rbinom(1, 1, phat[month])*rgamma(1, mhat[month], lahat[month])
  }
  total_raintank = c()
  total_grey = c()
  total_saved = c()
  total_gardenwater = c()
  
  for(i in 1:num_years) {
    year_raintank = c()
    year_grey = c()
    year_saved = c()
    year_gardenwater = c()
      
    for(j in 1:12){
      month_raintank = c()
      month_grey = c()
      month_saved = c()
      month_gardenwater = c()
      
      for(k in 1:days_in_month[j]){
        water_saved = 0
        rain = Xsim(j)
        rainwatertank = rainwatertank + rain*roofarea/(10^4)
        greywatertank = greywatertank + showersize*numshower + washsize*numwash()
        toil_flush = numflush()*flushsize
      
        #flushing the toilet
        
        if(greywatertank > 0 & greywatertank >= toil_flush) {
          greywatertank = greywatertank - toil_flush
          water_saved = water_saved + toil_flush
        }
        else if (greywatertank > 0 & toil_flush > greywatertank){
            left_over_flush = toil_flush - greywatertank
            water_saved = water_saved + greywatertank
            greywatertank = 0
            if(rainwatertank > 0 & rainwatertank >= left_over_flush) {
              rainwatertank = rainwatertank - left_over_flush
              water_saved = water_saved + left_over_flush}
            else if (rainwatertank > 0 & left_over_flush >= rainwatertank) {
              #mains are used
              water_saved = water_saved + rainwatertank
              rainwatertank = 0
            }
        }
        
        else if (rainwatertank > 0 & rainwatertank >= toil_flush) {
          rainwatertank = rainwatertank - toil_flush
          water_saved = water_saved + toil_flush
        }
        else if (rainwatertank > 0 & toil_flush >= rainwatertank) {
          water_saved = water_saved + rainwatertank
          rainwatertank = 0
          # the rest of the flush comes from the mains
          
        }
        #else no water is saved
        
        
        #garden water
        
        water_use = function() {

          if(water_need > ((rain + TWOdaysRain)/3)) {
            Need = gardenarea*(water_need - (rain + TWOdaysRain))
            
            if(rainwatertank >= Need) {
              rainwatertank = rainwatertank - Need
              water_saved = water_saved + Need
            }
            if(rainwatertank > 0 & Need > rainwatertank) {
              left_over_garden = Need - rainwatertank
              water_saved = water_saved + rainwatertank
              rainwatertank = 0
              
              if(greywatertank > 0 & greywatertank >= left_over_garden) {
                greywatertank = greywatertank - left_over_garden
                water_saved = water_saved + left_over_garden}
              
              else if (greywatertank > 0 & left_over_garden >= greywatertank) {
                #mains are used
                water_saved = water_saved + greywatertank
                greywatertank = 0
              }
            }
          }
        else{
          Need = 0
        }
        gardenwater = rain + Need
        out = c(rainwatertank, greywatertank, water_saved, gardenwater)
        return(out)
        }
        
        water_need = mean_max_temp[j]/15
        
        #starting out
        if( i > 1) {
          if( j == 1) {
            if(k == 1) {
              TWOdaysRain = total_gardenwater[length(total_gardenwater)-2] + 
                total_gardenwater[length(total_gardenwater) - 1]
              
              rainwatertank = water_use()[1]
              greywatertank = water_use()[2]
              water_saved = water_use()[3]
              gardenwater = water_use()[4]
              }
            if(k==2) {
              TWOdaysRain = total_gardenwater[length(total_gardenwater)-1] + 
                month_gardenwater[length(month_gardenwater) - 1]
              
              rainwatertank = water_use()[1]
              greywatertank = water_use()[2]
              water_saved = water_use()[3]
              gardenwater = water_use()[4]
              }
            }
          }
        if(i ==1) {
          if(j ==1) {
            if(k==1) {
              TWOdaysRain = 0
              rainwatertank = water_use()[1]
              greywatertank = water_use()[2]
              water_saved = water_use()[3]
              gardenwater = water_use()[4]
            }
            if(k==2) {
              TWOdaysRain = month_gardenwater[length(month_gardenwater)-1]
              rainwatertank = water_use()[1]
              greywatertank = water_use()[2]
              water_saved = water_use()[3]
              gardenwater = water_use()[4]
            }
          }
        }
        if(k == 1) {
          TWOdaysRain = year_gardenwater[length(year_gardenwater) -2] +
            year_gardenwater[length(year_gardenwater)-1]
          rainwatertank = water_use()[1]
          greywatertank = water_use()[2]
          water_saved = water_use()[3]
          gardenwater = water_use()[4]
        }
        if(k==2) {
          TWOdaysRain = year_gardenwater[length(year_gardenwater)-1] + 
            month_gardenwater[length(month_gardenwater)-1]
          rainwatertank = water_use()[1]
          greywatertank = water_use()[2]
          water_saved = water_use()[3]
          gardenwater = water_use()[4]
        }
        else {
          TWOdaysRain = month_gardenwater[length(month_gardenwater)-2] +
            month_gardenwater[length(month_gardenwater)-1]
          rainwatertank = water_use()[1]
          greywatertank = water_use()[2]
          water_saved = water_use()[3]
          gardenwater = water_use()[4]
        }
        
        month_saved = append(month_saved, water_saved)
        month_gardenwater = append(month_gardenwater, gardenwater)
        month_grey = append(month_grey, greywatertank)
        month_raintank = append(month_raintank, rainwatertank)
        
      
        
        }#end of days loop 
    
      year_saved = append(year_saved, month_saved)
      year_gardenwater = append(year_gardenwater, month_gardenwater)
      year_grey = append(year_grey, month_grey)
      year_raintank = append(year_raintank, month_raintank)
      
      } #end of months loop
    
    total_saved = append(total_saved, year_saved)
    total_gardenwater = append(total_gardenwater, year_gardenwater)
    total_grey = append(total_grey, year_grey)
    total_raintank = append(total_raintank, year_raintank)
      
    
      
    } #end of years loop
    
  if(plotflag == T) {
    plot(x = seq(1,num_years*365,1), y = total_grey,
         type = l, col = "red", ylab = "Litres",
         xlab = "Days")
    plot(x = seq(1,num_years*365,1), y = total_raintank,
         type = l, col = "blue" )
  }
  
    #insert plot here
  return(sum(total_saved))
    
}

tanksim(1,3000, 1000, p_est, m_est, lam_est, T)

