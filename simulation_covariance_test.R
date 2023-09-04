library(ggplot2)
library(cowplot)
source('function_proportional_odds_generation.R')
source('function_standardized_statistics_generation.R')

##### Scenario-1 #####
# proportional odds model holds

theta = c(.01,.2,.4,.8,1.6)
beta = 0
A = 6

pi = get_prob(theta, beta)
pi_1 = pi$pi_1
pi_0 = pi$pi_0

K = 3
mat_index = c(1,2,3,5,6,9) # indices of elements need in the covariance matrix of z

seed_base = 20230822


##### start simulation #####

a = 4

seed_seq = seq(1, 1000, 1) + seed_base
expo_seq = seq(2,4,.1)

var_mat_prop = NULL
var_mat_rank = NULL

for(ex in expo_seq){
  n = ceiling(10^ex)
  z_mat_prop = NULL
  z_mat_rank = NULL
  
  for(seed in seed_seq){
    set.seed(seed)
    y1 = sample(seq(1,A,1), size = n, replace = TRUE, prob = pi_1)
    y0 = sample(seq(1,A,1), size = n, replace = TRUE, prob = pi_0)

    z_prop = stand_prop(y1, y0, K)
    z_mat_prop = rbind(z_mat_prop, z_prop)

    z_rank = stand_rank(y1, y0, K, A)
    z_mat_rank = rbind(z_mat_rank, z_rank)
  }
  var_prop = as.vector(var(z_mat_prop))[mat_index]
  var_mat_prop = rbind(var_mat_prop, var_prop)
  
  var_rank = as.vector(var(z_mat_rank))[mat_index]
  var_mat_rank = rbind(var_mat_rank, var_rank)
}

plot_df_1_var = data.frame(
  log_units = rep(expo_seq, 2),
  var_1 = c(var_mat_prop[,1], var_mat_rank[,1]),
  cov_1_2 = c(var_mat_prop[,2], var_mat_rank[,2]),
  cov_1_3 = c(var_mat_prop[,3], var_mat_rank[,3]),
  var_2 = c(var_mat_prop[,4], var_mat_rank[,4]),
  cov_2_3 = c(var_mat_prop[,5], var_mat_rank[,5]), 
  var_3 = c(var_mat_prop[,6], var_mat_rank[,6]),
  methods = rep(c('Proportional', 'Wilcoxon'), each = length(expo_seq))
)


p_1_1 = ggplot(plot_df_1_var, aes(x = log_units, y = var_1, color = methods)) +
  geom_line(linewidth = .8) + theme_minimal() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Log_number', y = 'Var (Z1)') + coord_cartesian(ylim = c(0.8, 1.2))
p_1_2 = ggplot(plot_df_1_var, aes(x = log_units, y = var_2, color = methods)) +
  geom_line(linewidth = .8) + theme_minimal() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Log_number', y = 'Var (Z2)') + coord_cartesian(ylim = c(0.8, 1.2))
p_1_3 = ggplot(plot_df_1_var, aes(x = log_units, y = var_3, color = methods)) +
  geom_line(linewidth = .8) + theme_minimal() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Log_number', y = 'Var (Z3)') + coord_cartesian(ylim = c(0.8, 1.2))
p_1_4 = ggplot(plot_df_1_var, aes(x = log_units, y = cov_1_2, color = methods)) +
  geom_line(linewidth = .8) + theme_minimal() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Log_number', y = 'Cov (Z1, Z2)') + coord_cartesian(ylim = c(0.4, 0.9))
p_1_5 = ggplot(plot_df_1_var, aes(x = log_units, y = cov_1_3, color = methods)) +
  geom_line(linewidth = .8) + theme_minimal() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Log_number', y = 'Cov (Z1, Z3)') + coord_cartesian(ylim = c(0.4, 0.9))
p_1_6 = ggplot(plot_df_1_var, aes(x = log_units, y = cov_2_3, color = methods)) +
  geom_line(linewidth = .8) + theme_minimal() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Log_number', y = 'Cov (Z2, Z3)') + coord_cartesian(ylim = c(0.4, 0.9))

plot_grid(p_1_1, p_1_2, p_1_3, p_1_4, p_1_5, p_1_6, nrow = 2, byrow = TRUE)



##### Scenario-2 #####
# pi1 = pi0 = (.1,.1,.2,.2,.2,.2)

pi_1 = c(.1,.1,.2,.2,.2,.2)
pi_0 = c(.1,.1,.2,.2,.2,.2)
A = 6

K = 3
mat_index = c(1,2,3,5,6,9)

seed_base = 20230822


##### start simulation #####

a = 4

seed_seq = seq(1, 1000, 1) + seed_base
expo_seq = seq(2,4,.1)
var_mat_prop = NULL
var_mat_rank = NULL

for(ex in expo_seq){
  n = ceiling(10^ex)
  z_mat_prop = NULL
  z_mat_rank = NULL

  for(seed in seed_seq){
    set.seed(seed)
    y1 = sample(seq(1,A,1), size = n, replace = TRUE, prob = pi_1)
    y0 = sample(seq(1,A,1), size = n, replace = TRUE, prob = pi_0)

    z_prop = stand_prop(y1, y0, K)
    z_mat_prop = rbind(z_mat_prop, z_prop)

    z_rank = stand_rank(y1, y0, K, A)
    z_mat_rank = rbind(z_mat_rank, z_rank[k])
  }
  var_prop = as.vector(var(z_mat_prop))[mat_index]
  var_mat_prop = rbind(var_mat_prop, var_prop)

  var_rank = as.vector(var(z_mat_rank))[mat_index]
  var_mat_rank = rbind(var_mat_rank, var_rank)
}

plot_df_2_var = data.frame(
  log_units = rep(expo_seq, 2),
  var_1 = c(var_mat_prop[,1], var_mat_rank[,1]),
  cov_1_2 = c(var_mat_prop[,2], var_mat_rank[,2]),
  cov_1_3 = c(var_mat_prop[,3], var_mat_rank[,3]),
  var_2 = c(var_mat_prop[,4], var_mat_rank[,4]),
  cov_2_3 = c(var_mat_prop[,5], var_mat_rank[,5]),
  var_3 = c(var_mat_prop[,6], var_mat_rank[,6]),
  methods = rep(c('Proportional', 'Wilcoxon'), each = length(expo_seq))
)


p_2_1 = ggplot(plot_df_2_var, aes(x = log_units, y = var_1, color = methods)) +
  geom_line(linewidth = .8) + theme_minimal() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Log_number', y = 'Var (Z1)') + coord_cartesian(ylim = c(0.8, 1.2))
p_2_2 = ggplot(plot_df_2_var, aes(x = log_units, y = var_2, color = methods)) +
  geom_line(linewidth = .8) + theme_minimal() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Log_number', y = 'Var (Z2)') + coord_cartesian(ylim = c(0.8, 1.2))
p_2_3 = ggplot(plot_df_2_var, aes(x = log_units, y = var_3, color = methods)) +
  geom_line(linewidth = .8) + theme_minimal() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Log_number', y = 'Var (Z3)') + coord_cartesian(ylim = c(0.8, 1.2))
p_2_4 = ggplot(plot_df_2_var, aes(x = log_units, y = cov_1_2, color = methods)) +
  geom_line(linewidth = .8) + theme_minimal() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Log_number', y = 'Cov (Z1, Z2)') + coord_cartesian(ylim = c(0.4, 0.9))
p_2_5 = ggplot(plot_df_2_var, aes(x = log_units, y = cov_1_3, color = methods)) +
  geom_line(linewidth = .8) + theme_minimal() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Log_number', y = 'Cov (Z1, Z3)') + coord_cartesian(ylim = c(0.4, 0.9))
p_2_6 = ggplot(plot_df_2_var, aes(x = log_units, y = cov_2_3, color = methods)) +
  geom_line(linewidth = .8) + theme_minimal() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Log_number', y = 'Cov (Z2, Z3)') + coord_cartesian(ylim = c(0.4, 0.9))

plot_grid(p_2_1, p_2_2, p_2_3, p_2_4, p_2_5, p_2_6, nrow = 2, byrow = TRUE)
