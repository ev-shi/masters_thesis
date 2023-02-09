## Generate test negative data
## Allow for different levels of severity

library(survival)
library(tidyverse)
library(gam)

# population size
N <- 20000

# study end time
tau <- 1000

################################################################################

# vaccination coverage
pct_vax <- 0.50

# vaccine efficacy against infection (1-hazard ratio)
VE_sph <- 0.5
# vaccine efficacy against progression to symptoms given infection
VE_ph <- 0.67
# vaccine efficacy against progression to severe disease given symptoms
VE_h <- 0.775

################################################################################

# define hazard function for infection with target pathogen in reference group
# start out with a weibull
shape_tp <- 1.2
rate_tp0 <- 0.01

# baseline probability that an unvaccinated person progresses to symptoms given infection
prob_dx_tp <- 0.7
# baseline probability that an unvaccinated person progresses to severe disease given symptoms
prob_hos_tp <- 0.05

################################################################################

# define hazard function for testing negative in reference group
# start out with a weibull
shape_tn <- 0.8
rate_tn0 <- 0.01

# baseline probability that an unvaccinated person progresses to symptoms given infection
prob_dx_tp <- 0.5
# baseline probability that an unvaccinated person progresses to severe disease given symptoms
prob_hos_tp <- 0.01

################################################################################

# probability that a vaccinated person with no symptoms seeks testing
prob_test_s_vax <- 0.10
# probability that a vaccinated person with mild symptoms seeks testing
prob_test_p_vax <- 0.50
# probability that a vaccinated person with severe symptoms seeks testing
prob_test_h_vax <- 1.00

# probability that an unvaccinated person with no symptoms seeks testing
prob_test_s_unvax <- 0.05
# probability that an unvaccinated person with mild symptoms seeks testing
prob_test_p_unvax <- 0.20
# probability that an unvaccinated person with severe symptoms seeks testing
prob_test_h_unvax <- 0.90

################################################################################

# run the simulation
source("TNDfunctions.R")
View(df)

split_df <- split(df, df$severity)

# fit logistic regression model adjusting for age and calendar time
logfit <- glm(result ~ vax + age + s(time,3), data=split_df$`1`, family=binomial)
summary(logfit)
1-exp(logfit$coefficients[2])

logfit <- glm(result ~ vax + age + s(time,3), data=split_df$`2`, family=binomial)
summary(logfit)
1-exp(logfit$coefficients[2])

logfit <- glm(result ~ vax + age + s(time,3), data=split_df$`3`, family=binomial)
summary(logfit)
1-exp(logfit$coefficients[2])
