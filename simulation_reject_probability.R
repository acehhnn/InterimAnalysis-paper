source('function_reject_probability.R')
source('function_proportional_odds_generation.R')

##### Simulation function #####
p_simu = function(K, A, c, n_seq, pi_1, pi_0){
  p_prop_seq = NULL
  p_rank_seq = NULL
  
  for(n in n_seq){
    a = sequence(A)
    
    # simulation start!
    itr = 1e+4
    flag_prop = NULL
    flag_rank = NULL
    
    for(seed in 1:itr){
      set.seed(seed)
      y1 = sample(a, n, pi_1, replace = T)
      y0 = sample(a, n, pi_0, replace = T)
      
      flag_prop = rbind(flag_prop, test_prop_p(K, y1, y0, c))
      flag_rank = rbind(flag_rank, test_rank_p(K, y1, y0, c))
    }
    
    n_prop = colSums(flag_prop)
    n_rank = colSums(flag_rank)
    
    p_prop = n_prop / itr
    p_rank = n_rank / itr
    
    p_prop_seq = c(p_prop_seq, p_prop, sum(p_prop))
    p_rank_seq = c(p_rank_seq, p_rank, sum(p_rank))
  }
  df = data.frame(
    prop = p_prop_seq,
    rank = p_rank_seq,
    unit = rep(n_seq, each = (K + 1))
  )
  return(df)
}


##### Simulation start #####
####### Simulation I #######

K = 3
A = 6
c = get_boundary(K = K)
n_seq = c(100,1000)

## Setting-2-1
theta = c(.1,.2,.4,.8,1.6)
beta = .1

pi = get_prob(theta, beta)
df = p_simu(K, A, c, n_seq, pi$pi_1, pi$pi_0)


## Setting-2-2
theta = c(.1,.2,.4,.8,1.6)
beta = .5

pi = get_prob(theta, beta)
df = p_simu(K, A, c, n_seq, pi$pi_1, pi$pi_0)

## Setting-2-3
theta = c(.01,.04,.16,.64,1.28)
beta = .01

pi = get_prob(theta, beta)
df = p_simu(K, A, c, n_seq, pi$pi_1, pi$pi_0)

## Setting-2-4
theta = c(.01,.04,.16,.64,1.28)
beta = .05

pi = get_prob(theta, beta)
df = p_simu(K, A, c, n_seq, pi$pi_1, pi$pi_0)

## Setting-1-1
theta = c(.1,.2,.4,.8,1.6)
beta = .1

pi_1 = c(.05,.05,.05,.05,.05,.75)
pi_0 = c(.75,.05,.05,.05,.05,.05)

df = p_simu(K, A, c, n_seq, pi_1, pi_0)

## Setting-1-2
theta = c(.1,.2,.4,.8,1.6)
beta = .5

pi_1 = c(.05,.05,.05,.05,.4,.4)
pi_0 = c(.4,.4,.05,.05,.05,.05)

df = p_simu(K, A, c, n_seq, pi_1, pi_0)

## Setting-1-3
theta = c(.01,.04,.16,.64,1.28)
beta = .01

pi_1 = c(.1,.1,.2,.2,.2,.2)
pi_0 = c(.2,.2,.2,.2,.1,.1)

df = p_simu(K, A, c, n_seq, pi_1, pi_0)

## Setting-1-4
theta = c(.01,.04,.16,.64,1.28)
beta = .05

pi_1 = c(.16,.16,.16,.2,.16,.16)
pi_0 = c(.2,.16,.16,.16,.16,.16)

df = p_simu(K, A, c, n_seq, pi_1, pi_0)


##### Simulation start #####
###### Simulation II #######

K = 1
A = 6
c = get_boundary(K = K)
n_seq = c(100,1000)

## Setting-2-1
theta = c(.1,.2,.4,.8,1.6)
beta = .1

pi = get_prob(theta, beta)
df = p_simu(K, A, c, n_seq, pi$pi_1, pi$pi_0)

## Setting-2-2
theta = c(.1,.2,.4,.8,1.6)
beta = .5

pi = get_prob(theta, beta)
df = p_simu(K, A, c, n_seq, pi$pi_1, pi$pi_0)

## Setting-2-3
theta = c(.01,.04,.16,.64,1.28)
beta = .01

pi = get_prob(theta, beta)
df = p_simu(K, A, c, n_seq, pi$pi_1, pi$pi_0)

## Setting-2-4
theta = c(.01,.04,.16,.64,1.28)
beta = .05

pi = get_prob(theta, beta)
df = p_simu(K, A, c, n_seq, pi$pi_1, pi$pi_0)

## Setting-1-1
theta = c(.1,.2,.4,.8,1.6)
beta = .1

pi_1 = c(.05,.05,.05,.05,.05,.75)
pi_0 = c(.75,.05,.05,.05,.05,.05)

df = p_simu(K, A, c, n_seq, pi_1, pi_0)


## Setting-1-2
theta = c(.1,.2,.4,.8,1.6)
beta = .5

pi_1 = c(.05,.05,.05,.05,.4,.4)
pi_0 = c(.4,.4,.05,.05,.05,.05)

df = p_simu(K, A, c, n_seq, pi_1, pi_0)

## Setting-1-3
theta = c(.01,.04,.16,.64,1.28)
beta = .01

pi_1 = c(.1,.1,.2,.2,.2,.2)
pi_0 = c(.2,.2,.2,.2,.1,.1)

df = p_simu(K, A, c, n_seq, pi_1, pi_0)

## Setting-1-4
theta = c(.01,.04,.16,.64,1.28)
beta = .05

pi_1 = c(.16,.16,.16,.2,.16,.16)
pi_0 = c(.2,.16,.16,.16,.16,.16)

df = p_simu(K, A, c, n_seq, pi_1, pi_0)
