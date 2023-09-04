library(MASS)

# proportional odds model
stand_prop = function(y1, y0, K, t = NULL){
  if(is.null(t)){
    t = seq(1/K, 1, 1/K)
  }
  if(length(t) != K){
    stop('Length of t must be equal to K.')
  }
  
  n1 = length(y1) # number of patients in the treatment group
  n0 = length(y0) # number of patients in the control group
  
  z_seq = NULL # a sequence of standardized statistics
  
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
    
    fit_k = polr(y_k ~ x_k, df_k, method = "logistic", Hess = TRUE) # fit the proportional odds model
    beta = fit_k$coefficients # estimate beta
    sigma = (solve(fit_k$Hessian))[1,1] # estimate variance of beta
    z = beta / sqrt(sigma) # standardized statistics
    
    z_seq = append(z_seq, z)
  }
  return(z_seq)
}

# Wilcoxon rank sum test
stand_rank = function(y1, y0, K, A, t = NULL){
  if(is.null(t)){
    t = seq(1/K, 1, 1/K)
  }
  if(length(t) != K){
    stop('Length of t must be equal to K.')
  }
  
  n1 = length(y1) # number of patients in the treatment group
  n0 = length(y0) # number of patients in the control group
  
  z_seq = NULL # a sequence of standardized statistics
  
  for(k in 1:K){
    n1_k = min(ceiling(n1 * t[k]), n1)
    n0_k = min(ceiling(n0 * t[k]), n0)
    
    y_k = c(y1[1:n1_k], y0[1:n0_k])
    
    d_k = NULL # number of patients in each level
    for(a in 1:A){
      d_k = append(d_k, sum(y_k == a))
    }
    
    r_k = rank(y_k)
    w_k = sum(r_k[1:n1_k]) # rank sum
    mu_k = n1_k * (n1_k + n0_k +1) / 2 # mean
    var_k = n1_k * n0_k * (n1_k + n0_k + 1) / 12 - n1_k * n0_k * sum(d_k^3 - d_k) / (12 * (n1_k + n0_k) * (n1_k + n0_k - 1)) # adjusted variance
    z = (w_k - mu_k) / sqrt(var_k) # standardized statistics
    z_seq = append(z_seq, z)
  }
  return(z_seq)
}