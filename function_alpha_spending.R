pocock = function(alpha, t){
  f = NULL
  for(i in 1:length(t)){
    f = append(f, alpha * min(log(1 + (exp(1) - 1) * t[i]), 1))
  }
  return(f)
}

obf = function(alpha, t){
  f = NULL
  for(i in 1:length(t)){
    if(t[i] != 0){
      f = append(f, min(2 * (1 - pnorm(qnorm(1 - alpha / 2) / sqrt(t[i]))), alpha))
    }
    else{
      f = append(f, 0)
    }
  }
  return(f)
}

gamma_family = function(alpha, gamma, t){
  f = NULL
  if(gamma == 0){
    for(i in 1:length(t)){
      f = append(f, alpha * min(t[i], 1))
    }
  }
  else{
    for(i in 1:length(t)){
      f = append(f, alpha * min(1, (1 - exp(-gamma * t[i])) / (1 - exp(-gamma)), 1))
    }
  }
  return(f)
}
