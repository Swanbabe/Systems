## Assignment 2: A Game of Squash


rm(list=ls())

##########################################################################################################

## ___ Status of the game ___


status = function(x,y) {
  #' for integer inputs x and y, returns:
  #' 
  #' 'unfinished' if the game has not yet finished; 
  #'  
  #' 'player 1 win' if player 1 has won the game; 
  #' 
  #' 'player 2 win' if player 2 has won the game; 
  #' 
  #' 'impossible' if x and y are impossible scores. 
  
  im = 'impossible'
  p1 = 'player 1 win'
  p2 = 'player 2 win'
  un = 'unfinished'
  
  if (x < 0 | y < 0) return(im)
  
  if (x == 9 & x - y > 1) return(p1)
  else if (x > 9 & x - y == 2) return (p1)
  
  if (y == 9 & y - x > 1) return(p2)
  else if (y > 9 & y - x == 2) return(p2)
  
  if (x > 9 & x - y > 2) return(im)
  if (y > 9 & y - x > 2) return(im)
  
  else return(un)
}


# From the spuRs library:

# Program spuRs/resources/scripts/status.test.r

status.test <- function(s.ftn) {
  
  x.vec <- (-1):11
  y.vec <- (-1):11
  plot(x.vec, y.vec, type = "n", xlab = "player 1", ylab = "player 2")
  
  for (x in x.vec) {
    for (y in y.vec) {
      s <- s.ftn(x, y)
      if (s == "impossible") text(x, y, "X", col = "red")
      else if (s == "unfinished") text(x, y, "?", col = "blue")
      else if (s == "player 1 win") text(x, y, "1", col = "green")
      else if (s == "player 2 win") text(x, y, "2", col = "green")
    }
  }
  return(invisible(NULL))
}



status.test(status)

##########################################################################################################

## ___ Simulating a game ___


play_point = function(state, a, b) {
  #' We simulate a point played based on a binomial distrubition with parameters of 1 observation
  #' and 1 trial, with probabilities
  #'  'a': the probability that player 1 wins a point if player 1 serves,  or
  #'  'b': the probability that player 1 wins a point if player 2 serves
  #' 
  #' state will be the vector (x,y,z) where x: number of points won by player 1, y: number of points won by player 2
  #' z = 1 if player is serving and z=2 if player 2 is serving
  
  if (state[3] == 1) {
    bin1 = rbinom(1,1,a)
    new_state = c( bin1, 0, 1 - bin1 )
  }
  if (state[3] == 2) {
    bin2 = rbinom(1,1,b)
    new_state = c( 0, 1 - bin2, - bin2)
  }
  state = state + new_state
}


# Program spuRs/resources/scripts/play_game.r

play_game <- function(a, b) {
  state <- c(0, 0, 1)
  while (status(state[1], state[2]) == "unfinished") {
    # show(state)
    state <- play_point(state, a, b)
  }
  if (status(state[1], state[2]) == "player 1 win") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


MySeed = 3133

squash_prob = function(a,b,n,seed_gen){
  #' For producing the average probability of player 1 wining a game based on probabilities 
  #' 'a' and 'b' as defined above, and on 'n games played
  
  game = c()
  set.seed(seed_gen)
  
  for (k in 1:n) { game[k] = play_game(a,b) }
  
  return(mean(game))
}


p_hat = c()

for (i in 1:12) { p_hat = append(p_hat, squash_prob(0.55,0.45,2^i,MySeed)) }

p_hat

plot(1:12, p_hat, type = 'p', main = "Estimated probability for sample size 'n'", xlab = "log(n)/log(2)", ylab = "p_hat")


##########################################################################################################

## ___ Probability of Winning ___


biggest_n = 10^4

A = seq(0.1,0.9,0.1)
B = seq(0.1,0.9,0.1)

est_p_ab = matrix( nrow = 9, ncol = 9 )

for (x in 1:9) { for (y in 1:9) { est_p_ab[x,y] = squash_prob(A[x], B[y], biggest_n, MySeed) }}

Est = round(est_p_ab,2)

View(Est)

print(xtable(Est))

##########################################################################################################

## ___ Length of a Game ___


game_length <- function(a, b) {
  #'the length of a game played (until state not 'unfinished'), with 'a','b' as defined above
  #'
  state <- c(0, 0, 1)
  points_played = 0
  
  while (status(state[1], state[2]) == "unfinished") {
    # show(state)
    
    state <- play_point(state, a, b)
    
    points_played = points_played + 1
    
  }
  return(points_played)
  
}

squash_length_prob = function(a,b,n,seed_gen){
  #' the estimated average game length for 'n' games, with 'a','b' as defined above
  #' 
  game = c()
  set.seed(seed_gen)
  
  for (k in 1:n) { game[k] = game_length(a,b) }
  
  return(mean(game))
}


length_p_ab = matrix( nrow = 9, ncol = 9)

for (x in 1:9) { for (y in 1:9) { length_p_ab[x,y] = squash_length_prob(A[x], B[y], biggest_n, MySeed)}}

Len = round(length_p_ab,2)

View(Len)

print(xtable(Len))
