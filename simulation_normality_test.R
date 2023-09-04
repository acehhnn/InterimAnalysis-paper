library(ggplot2)
library(cowplot)
source('function_proportional_odds_generation.R')
source('function_standardized_statistics_generation.R')

##### Scenario-1 #####
# Proportional odds model holds

theta = c(.01,.2,.4,.8,1.6)
beta = 0
A = 6

pi = get_prob(theta, beta)
pi_1 = pi$pi_1
pi_0 = pi$pi_0

K = 3

seed_base = 20230822


##### start simulation #####

k = 2
a = 4

seed_seq = seq(1, 1000, 1) + seed_base
expo_seq = seq(2,4,.1) # a sequence of sample sizes
ks_seq_prop = NULL # a sequence of Kolmogorov-Smirnov statistics of proportional odds model approach
ks_seq_rank = NULL # a sequence of Kolmogorov-Smirnov statistics of Wilcoxon rank sum approach

for(ex in expo_seq){
  n = ceiling(10^ex)
  z_seq_prop = NULL # a sequence of standardized statistics of proportional odds model approach
  z_seq_rank = NULL # a sequence of standardized statistics of Wilcoxon rank sum approach
  
  for(seed in seed_seq){
    set.seed(seed)
    y1 = sample(seq(1,A,1), size = n, replace = TRUE, prob = pi_1)
    y0 = sample(seq(1,A,1), size = n, replace = TRUE, prob = pi_0)
    
    z_prop = stand_prop(y1, y0, K)
    z_seq_prop = append(z_seq_prop, z_prop[k])
    
    z_rank = stand_rank(y1, y0, K, A)
    z_seq_rank = append(z_seq_rank, z_rank[k])
  }
  
  ks_prop = ks.test(z_seq_prop, 'pnorm')
  ks_seq_prop = append(ks_seq_prop, ks_prop$statistic)
  
  ks_rank = ks.test(z_seq_rank, 'pnorm')
  ks_seq_rank = append(ks_seq_rank, ks_rank$statistic)
}

plot_df_1_normal = data.frame(
  log_units = rep(expo_seq, 2),
  test_statistics = c(ks_seq_prop, ks_seq_rank),
  methods = rep(c('Proportional', 'Wilcoxon'), each = length(expo_seq))
)


##### Scenario-2 #####
# pi1 = pi0 = (.1,.1,.2,.2,.2,.2)

pi_1 = c(.1,.1,.2,.2,.2,.2)
pi_0 = c(.1,.1,.2,.2,.2,.2)
A = 6

K = 3

seed_base = 20230822


##### start simulation #####

k = 2
a = 4

seed_seq = seq(1, 1000, 1) + seed_base
expo_seq = seq(2,4,.1)
ks_seq_prop = NULL
ks_seq_rank = NULL

for(ex in expo_seq){
  n = ceiling(10^ex)
  z_seq_prop = NULL
  z_seq_rank = NULL

  for(seed in seed_seq){
    set.seed(seed)
    y1 = sample(seq(1,A,1), size = n, replace = TRUE, prob = pi_1)
    y0 = sample(seq(1,A,1), size = n, replace = TRUE, prob = pi_0)

    z_prop = stand_prop(y1, y0, K)
    z_seq_prop = append(z_seq_prop, z_prop[k])

    z_rank = stand_rank(y1, y0, K, A)
    z_seq_rank = append(z_seq_rank, z_rank[k])
  }
  ks_prop = ks.test(z_seq_prop, 'pnorm')
  ks_seq_prop = append(ks_seq_prop, ks_prop$statistic)

  ks_rank = ks.test(z_seq_rank, 'pnorm')
  ks_seq_rank = append(ks_seq_rank, ks_rank$statistic)
}

plot_df_2_normal = data.frame(
  log_units = rep(expo_seq, 2),
  test_statistics = c(ks_seq_prop, ks_seq_rank),
  methods = rep(c('Proportional', 'Wilcoxon'), each = length(expo_seq))
)


p1 = ggplot(plot_df_1_normal, aes(x = log_units, y = test_statistics, color = methods)) +
  geom_line() + theme_minimal() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Log_number', y = 'Kolmogorov-Smirnov statistics') + coord_cartesian(ylim = c(0.015, 0.05))

p2 = ggplot(plot_df_2_normal, aes(x = log_units, y = test_statistics, color = methods)) + 
  geom_line() + theme_minimal() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Log_number', y = 'Kolmogorov-Smirnov statistics') + coord_cartesian(ylim = c(0.015, 0.05))

plot_grid(p1, p2)
