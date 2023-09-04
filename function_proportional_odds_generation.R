get_prob = function(theta, beta){
  A = length(theta) + 1
  
  pi_1 = NULL
  pi_0 = NULL
  
  for(i in 1:A){
    if(i==1){
      inc_1 = exp(theta[i] + beta) / (1 + exp(theta[i] + beta))
      inc_0 = exp(theta[i]) / (1 + exp(theta[i]))
      pi_1 = append(pi_1, inc_1)
      pi_0 = append(pi_0, inc_0)
    }else if(i!=A){
      inc_1 = (exp(theta[i] + beta) - exp(theta[i-1] + beta)) / ((1 + exp(theta[i] + beta)) * (1 + exp(theta[i-1] + beta)))
      inc_0 = (exp(theta[i]) - exp(theta[i-1])) / ((1 + exp(theta[i])) * (1 + exp(theta[i-1])))
      pi_1 = append(pi_1, inc_1)
      pi_0 = append(pi_0, inc_0)
    }else{
      inc_1 = 1 - exp(theta[i-1] + beta) / (1 + exp(theta[i-1] + beta))
      inc_0 = 1 - exp(theta[i-1]) / (1 + exp(theta[i-1]))
      pi_1 = append(pi_1, inc_1)
      pi_0 = append(pi_0, inc_0)
    }
  }
  pi = data.frame(pi_1 = pi_1, pi_0 = pi_0)
  return(pi)
}