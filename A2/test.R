
play_point = function(state, a, b) {
  #' We simulate a point played based on a binomial distrubition with parameters of 1 observation
  #' and 1 trial, with probability 'a' or 'b'
  #' 
  #' state will be the vector (x,y,z) where x: number of points won by player 1, y: number of points won by player 2
  #' z = 1 if player is serving and z=2 if player 2 is serving
  
  if (state[3] == 1) {
    bin1 = rbinom(1,1,a)
    new_state = c( bin1, 0, 1 - bin1 )
  }
  if (state[3] == 2) {
    bin2 = rbinom(1,1,b)
    new_state = c( 0, bin2, bin2 - 1)
  }
  state = state + new_state
}
