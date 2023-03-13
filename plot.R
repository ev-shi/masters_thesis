# plotting

library(survival)
library(tidyverse)
library(gam)
library(Metrics)

# population size
N <- 10000

# study end time
tau <- 1000

################################################################################

# vaccination coverage
pct_vax <- 0.5

# vaccine efficacy against infection (1-hazard ratio)
VE_s <- 0.7
# vaccine efficacy against progression to symptoms given infection
VE_p_vec <- seq(0, 0.9, 0.1)
# vaccine efficacy against progression to severe disease given symptoms
VE_h_vec <- seq(0, 0.9, 0.1)

################################################################################

# percent of population that is older
pct_older <- 0.5
# hazard ratio relating older age and testing positive
hr_older_tp <- 1.2
# hazard ratio relating older age and testing negative
hr_older_tn <- 2.0

# define hazard function for infection with target pathogen in reference group
# start out with a weibull
shape_tp <- 1.2
rate_tp0 <- 0.0015

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
rate_tn0 <- 0.015

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
prob_test_1_vax <- 0.7
# probability that a vaccinated person with mild symptoms seeks testing
prob_test_2_vax <- 0.7
# probability that a vaccinated person with severe symptoms seeks testing
prob_test_3_vax <- 0.7

# probability that an unvaccinated person with no symptoms seeks testing
prob_test_1_unvax <- 0.7
# probability that an unvaccinated person with mild symptoms seeks testing
prob_test_2_unvax <- 0.7
# probability that an unvaccinated person with severe symptoms seeks testing
prob_test_3_unvax <- 0.7

################################################################################

# run the simulation
set.seed(12345)

# create dataframe for bias output
cols <- c("est_VE_s", "est_VE_sp", "est_VE_sph")
VE_p_df <- data.frame(matrix(ncol=length(cols),nrow=0))
colnames(VE_p_df) <- cols

# changing VE_p
for(i in 1:length(VE_p_vec)){
  VE_p <- VE_p_vec[i]
  VE_h <- 0.2
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
  
  VE_p_df <- VE_p_df %>% add_row(est_VE_s = est_VE_s, est_VE_sp = est_VE_sp, 
                             est_VE_sph = est_VE_sph)
}

#saveRDS(VE_p_df, file = "estimates_from_changing_VE_p")

#jpeg(file = "VE_s vs VE_p.jpeg")
plot(VE_p_vec, VE_p_df$est_VE_s, main = "VE_s estimate vs. VE_p",
     xlab = "VE_p", ylab = "VE_s estimate")
#dev.off()

#jpeg(file = "VE_sp vs VE_p.jpeg")
plot(VE_p_vec, VE_p_df$est_VE_sp, main = "VE_sp estimate vs. VE_p",
     xlab = "VE_p", ylab = "VE_sp estimate")
#dev.off()

# jpeg(file = "VE_sph vs VE_p.jpeg")
plot(VE_p_vec, VE_p_df$est_VE_sph, main = "VE_sph estimate vs. VE_p",
     xlab = "VE_p", ylab = "VE_sph estimate")
# dev.off()

# changing VE_h
cols <- c("est_VE_s", "est_VE_sp", "est_VE_sph")
VE_h_df <- data.frame(matrix(ncol=length(cols),nrow=0))
colnames(VE_h_df) <- cols

for(i in 1:length(VE_h_vec)){
  VE_h <- VE_h_vec[i]
  VE_p <- 0.2
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
  
  VE_h_df <- VE_h_df %>% add_row(est_VE_s = est_VE_s, est_VE_sp = est_VE_sp, 
                             est_VE_sph = est_VE_sph)
}

# saveRDS(VE_h_df, file = "estimates_from_changing_VE_h")
# 
# jpeg(file = "VE_s vs VE_h.jpeg")
plot(VE_p_vec, VE_h_df$est_VE_s, main = "VE_s estimate vs. VE_h",
     xlab = "VE_h", ylab = "VE_s estimate")
# dev.off()
# 
# jpeg(file = "VE_sp vs VE_h.jpeg")
plot(VE_p_vec, VE_h_df$est_VE_sp, main = "VE_sp estimate vs. VE_h",
     xlab = "VE_h", ylab = "VE_sp estimate")
# dev.off()
# 
# jpeg(file = "VE_sph vs VE_h.jpeg")
plot(VE_p_vec, VE_h_df$est_VE_sph, main = "VE_sph estimate vs. VE_h",
     xlab = "VE_h", ylab = "VE_sph estimate")
# dev.off()