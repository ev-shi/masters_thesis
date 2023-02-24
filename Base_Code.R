## Generate test negative data
## Allow for different levels of severity

library(survival)
library(tidyverse)
library(gam)
library(Metrics)

# Number of simulation runs
num_runs <- 10

# population size
N <- 10000

# study end time
tau <- 1000

################################################################################

# vaccination coverage
pct_vax <- 0.5

# vaccine efficacy against infection (1-hazard ratio)
VE_s <- 0.6
# vaccine efficacy against progression to symptoms given infection
VE_p <- 0.4
# vaccine efficacy against progression to severe disease given symptoms
VE_h <- 0.2

# true value of VE_sp
VE_sp <- 0.76
# true value of VE_sph
VE_sph <- 0.808

# percent of population that is older
pct_older <- 0.5
# hazard ratio relating older age and testing positive
hr_older_tp <- 1.0
# hazard ratio relating older age and testing negative
hr_older_tn <- 2.0

# can further add another layer that allows different vaccination probability

################################################################################

# define hazard function for infection with target pathogen in reference group
# start out with a weibull
shape_tp <- 1.2
rate_tp0 <- 0.01

# baseline probability that a younger unvax person progresses to symptoms given infection
prob_dx_tp_young <- 0.8
# baseline probability that a younger unvax person progresses to severe disease given symptoms
prob_hos_tp_young <- 0.45

# baseline probability that an older unvax person progresses to symptoms given infection
prob_dx_tp_old <- 0.8
# baseline probability that an older unvax person progresses to severe disease given symptoms
prob_hos_tp_old <- 0.45

################################################################################

# define hazard function for testing negative in reference group
# start out with a weibull
shape_tn <- 0.8
rate_tn0 <- 0.01

# baseline probability that a younger unvax person progresses to symptoms given infection
prob_dx_tn_young <- 0.5
# baseline probability that a younger unvax person progresses to severe disease given symptoms
prob_hos_tn_young <- 0.05

# baseline probability that an older unvax person progresses to symptoms given infection
prob_dx_tn_old <- 0.5
# baseline probability that an older unvax person progresses to severe disease given symptoms
prob_hos_tn_old <- 0.05

################################################################################

# probability that a vaccinated person with no symptoms seeks testing
prob_test_1_vax <- 0.2
# probability that a vaccinated person with mild symptoms seeks testing
prob_test_2_vax <- 0.7
# probability that a vaccinated person with severe symptoms seeks testing
prob_test_3_vax <- 0.9

# probability that an unvaccinated person with no symptoms seeks testing
prob_test_1_unvax <- 0.2
# probability that an unvaccinated person with mild symptoms seeks testing
prob_test_2_unvax <- 0.7
# probability that an unvaccinated person with severe symptoms seeks testing
prob_test_3_unvax <- 0.9

## can further modify this to depend on age

################################################################################

# run the simulation
set.seed(12345)

# create dataframe for bias output
cols <- c("est_VE_s", "est_VE_sp", "est_VE_sph")
VE_df <- data.frame(matrix(ncol=length(cols),nrow=0))
colnames(VE_df) <- cols

for(i in 1:num_runs){
  source("Functions.R")
  
  df_keep <- df[df$keep==1,]
  df_s <- df_keep
  df_sp <- df_keep[df_keep$severity==2 | df_keep$severity==3,]
  df_sph <- df_keep[df_keep$severity==3,]
  
  # fit logistic regression model adjusting for age and calendar time
  # finding VE_s
  logfit_s <- glm(result ~ vax + age + s(time,3), data=df_s, family=binomial)
  est_VE_s <- 1-exp(logfit_s$coefficients[2])
  
  # finding VE_sp
  logfit_sp <- glm(result ~ vax + age + s(time,3), data=df_sp, family=binomial)
  est_VE_sp <- 1-exp(logfit_sp$coefficients[2])
  
  # finding VE_sph
  logfit_sph <- glm(result ~ vax + age + s(time,3), data=df_sph, family=binomial)
  est_VE_sph <- 1-exp(logfit_sph$coefficients[2])
  
  VE_df <- VE_df %>% add_row(est_VE_s = est_VE_s, est_VE_sp = est_VE_sp, 
                             est_VE_sph = est_VE_sph)
}

view(VE_df)

### Results
cols <- c("avg_est_VE_s", "avg_est_VE_sp", "avg_est_VE_sph",
          "mean_bias_s", "mean_bias_sp", "mean_bias_sph", 
          "mcse_bias_s", "mcse_bias_sp", "mcse_bias_sph",
          "root_mse_s", "root_mse_sp", "root_mse_sph")
results_df <- data.frame(matrix(ncol=length(cols),nrow=0))
colnames(results_df) <- cols

# Averages
avg_est_VE_s <- mean(est_VE_s)
avg_est_VE_sp <- mean(est_VE_sp)
avg_est_VE_sph <- mean(est_VE_sph)

# Calculating mean biases
mean_bias_s <- mean(VE_df$est_VE_s - VE_s)
mean_bias_sp <- mean(VE_df$est_VE_sp - VE_sp)
mean_bias_sph <- mean(VE_df$est_VE_sph - VE_sph)

# Monte Carlo Standard Error of Bias
mcse_bias_s <- sqrt((1/(num_runs * (num_runs - 1))) * sum((VE_df$est_VE_s - VE_s)^2))
mcse_bias_sp <- sqrt((1/(num_runs * (num_runs - 1))) * sum((VE_df$est_VE_sp - VE_sp)^2))
mcse_bias_sph <- sqrt((1/(num_runs * (num_runs - 1))) * sum((VE_df$est_VE_sph - VE_sph)^2))

# Calculating Square Root MSE
root_MSE_s <- sqrt(mean((VE_s - VE_df$est_VE_s)^2))
root_MSE_sp <- sqrt(mean((VE_sp - VE_df$est_VE_sp)^2))
root_MSE_sph <- sqrt(mean((VE_sph - VE_df$est_VE_sph)^2))


# Mean Absolute Error
#MAE_s <- mae(VE_df$est_VE_s, 0.85)
#MAE_sp <- mae(VE_df$est_VE_sp, 0.85)
#MAE_sph <- mae(0.85, VE_df$est_VE_sph)

results_df <- results_df %>% add_row(avg_est_VE_s = avg_est_VE_s, 
                                     avg_est_VE_sp = avg_est_VE_sp, 
                                     avg_est_VE_sph = avg_est_VE_sph,
                                     mean_bias_s = mean_bias_s, 
                                     mean_bias_sp = mean_bias_sp, 
                                     mean_bias_sph = mean_bias_sph, 
                                     mcse_bias_s = mcse_bias_s, 
                                     mcse_bias_sp = mcse_bias_sp, 
                                     mcse_bias_sph = mcse_bias_sph,
                                     root_mse_s = root_MSE_s, 
                                     root_mse_sp = root_MSE_sp, 
                                     root_mse_sph = root_MSE_sph)
view(results_df)

# Monte Carlo Standard Error
# average standard error 

