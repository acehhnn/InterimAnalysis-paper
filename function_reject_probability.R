library(MASS)
source('function_alpha_spending.R')

# critical values under normality assumption
get_boundary = function(K, alpha = 0.05, t = NULL, func = 'obf', gamma = 0){
  if(is.null(t)){
    t = seq(1/K, 1, 1/K)
  }
  if(length(t) != K){
    stop('Length of t must be equal to K.')
  }
  
  c = NULL
  S = matrix(rep(0,K*K), nrow = K) # Sigma
  for(i in 1:K){
    for(j in 1:K){
      S[i,j] = min(t[i],t[j]) / max(t[i],t[j])
    }
  }
  
  mean = rep(0,K)
  sigma = S
  n = 1e+4
  set.seed(20230810)
  simu = mvrnorm(n, mean, sigma) # generate random variables
  
  if(func == 'obf'){
    a = obf(alpha, t)
  }else if(func == 'pocock'){
    a = pocock(alpha, t)
  }else{
    a = gamma_family(alpha, gamma, t)
  }
  
  for(k in 1:length(t)){ # get d
    if(k == 1){
      c0 = qnorm(1 - a[k])
      c = append(c, c0)
    }
    else{
      delta = a[k] - a[k - 1]
      itr = seq(c[k - 1], 0, -0.01) # set of critical values
      for(d in itr){
        s = 0 # number of OK points
        for(t in 1:n){
          flag = 1
          for(j in 1:(k - 1)){
            if(simu[t,j] >= c[j]){
              flag = 0
              break
            }
          }
          if(flag == 1){
            if(simu[t,k] <= d){
              flag = 0
            }else{
              s = s + 1
            }
          }
        }
        prob = s / n
        if(prob >= delta){
          c = append(c,d)
          break
        }
      }
    }
  }
  return(c)
}

# proportional odds model
test_prop_p = function(K, y1, y0, c, t = NULL){
  if(is.null(t)){
    t = seq(1/K, 1, 1/K)
  }
  if(length(t) != K){
    stop('Length of t must be equal to K.')
  }
  
  flag = rep(0, K)
  
  n1 = length(y1)
  n0 = length(y0)
  
  for(k in 1:K){
    n1_k = min(ceiling(n1 * t[k]), n1)
    n0_k = min(ceiling(n0 * t[k]), n0)
    
    y1_k = y1[1:n1_k]
    y0_k = y0[1:n0_k]
    
    df_k = data.frame(
      y_k = c(y1_k, y0_k),
      x_k = c(rep(1, n1_k), rep(0, n0_k))
    )
    df_k$y_k = as.factor(df_k$y_k)
    
    fit_k = polr(y_k ~ x_k, df_k, method = "logistic", Hess = TRUE)
    beta = fit_k$coefficients
    sigma = (solve(fit_k$Hessian))[1,1]
    z = beta / sqrt(sigma) 
    if(z >= c[k]){
      flag[k] = 1
      break
    }
  }
  return(flag)
}

# Using rank sum statistics
test_rank_p = function(K, y1, y0, c, t = NULL){
  if(is.null(t)){
    t = seq(1/K, 1, 1/K)
  }
  if(length(t) != K){
    stop('Length of t must be equal to K.')
  }
  
  flag = rep(0, K)
  
  n1 = length(y1)
  n0 = length(y0) 
  
  for(k in 1:K){
    n1_k = min(ceiling(n1 * t[k]), n1)
    n0_k = min(ceiling(n0 * t[k]), n0)
    
    y_k = c(y1[1:n1_k], y0[1:n0_k])
    d_k = NULL
    for(a in 1:A){
      d_k = append(d_k, sum(y_k == a))
    }
    r_k = rank(y_k)
    w_k = sum(r_k[1:n1_k])
    mu_k = n1_k * (n1_k + n0_k +1) / 2
    var_k = n1_k * n0_k * (n1_k + n0_k + 1) / 12 - n1_k * n0_k * sum(d_k^3 - d_k) / (12 * (n1_k + n0_k) * (n1_k + n0_k - 1))
    z = (w_k - mu_k) / sqrt(var_k)
  
    if(z >= c[k]){
      flag[k] = 1
      break
    }
  }
  return(flag)
}