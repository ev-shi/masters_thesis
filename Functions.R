# supporting functions

################################################################################

# create data frame
cols <- c("id", "vax", "age", "time", "severity", "result", "keep")
df <- data.frame(matrix(ncol=length(cols),nrow=0))
colnames(df) <- cols

# Generate testing data for each person
for(i in 1:N){
  
  # Generate vaccination status
  vax <- rbinom(1,1,pct_vax)
  
  # Generate age
  age <- rbinom(1,1,pct_older)
  
  ## Generate test positive from a Weibull proportional hazards model
  
  # Calculate person's TP hazard ratio as a function of their covariates
  hr_tp <- ((1-VE_s)^vax)*((hr_older_tp)^age)
  
  # Updated rate parameter is a function of hazard ratio, baseline rate, shape
  rate_tp <- (hr_tp)^(1/shape_tp)*rate_tp0
  
  # Generate a Weibull by generating a uniform random variable
  Ti <- (-log(1-runif(1)))^(1/shape_tp)/rate_tp
  # this is equivalent to qweibull()
  # qweibull(runif(1),shape=shape_tp,scale=rate_tp^(-1))
  
  # If time is before the study end time, determine severity
  # Store the test result
  if(Ti <= tau){
    # set probability of progression based on age
    if(age==1){prob_dx_tp <- prob_dx_tp_old; prob_hos_tp <- prob_hos_tp_old}
    if(age==0){prob_dx_tp <- prob_dx_tp_young; prob_hos_tp <- prob_hos_tp_young}
    
    # this determines whether the person develops symptoms
    if(rbinom(1,1,prob_dx_tp*(1-VE_p)^vax)==1){
      # this determines whether the person develops severe disease given symptoms
      if(rbinom(1,1,prob_hos_tp*(1-VE_h)^vax)==1){
        severity <- 3
      }else{
        severity <- 2
      }
    }else{
      severity <- 1
    }
    
    # adjust the probability that someone seeks testing
    if(vax==1){
      if(severity==1){prob_test <- prob_test_1_vax}
      if(severity==2){prob_test <- prob_test_2_vax}
      if(severity==3){prob_test <- prob_test_3_vax}
    }else{
      if(severity==1){prob_test <- prob_test_1_unvax}
      if(severity==2){prob_test <- prob_test_2_unvax}
      if(severity==3){prob_test <- prob_test_3_unvax}        
    }
    # store the result only if they seek testing
    if(rbinom(1,1,prob_test)){
      df <- df %>% add_row(id=i, vax=vax, age=age, time=round(Ti), 
                            severity=severity, result = 1, keep = 1)
    } else {
      df <- df %>% add_row(id=i, vax=vax, age=age, time=round(Ti), 
                           severity=severity, result = 1, keep = 0)  
      }
  }
  
  ## Generate recurrent test negatives from a Weibull proportional hazards model
  
  # Make a test to control the while loop
  z=TRUE
  
  # Set the last test negative time as 0 to start
  last_tn <- 0
  
  while(z){
    
    # Calculate person's TN hazard ratio as a function of their covariates
    hr_tn <- ((hr_older_tn)^age)
    
    # Updated rate parameter is a function of hazard ratio, baseline rate, shape
    rate_tn <- (hr_tn)^(1/shape_tn)*rate_tn0
    
    # Generate a Weibull that conditions on the previous event time
    Ui <- last_tn + (-log(1-runif(1)) + (rate_tn*last_tn)^shape_tn )^(1/shape_tn)/rate_tn
    
    # Store the test if it is less than the test positive time Ti and the study end time
    if(Ui <= Ti & Ui <= tau){
      
      # set probability of progression based on age
      if(age==1){prob_dx_tn <- prob_dx_tn_old; prob_hos_tn <- prob_hos_tn_old}
      if(age==0){prob_dx_tn <- prob_dx_tn_young; prob_hos_tn <- prob_hos_tn_young}
      
      # this determines whether the person develops symptoms
      if(rbinom(1,1,prob_dx_tn)==1){
        # this determines whether the person develops severe disease given symptoms
        if(rbinom(1,1,prob_hos_tp)==1){
          severity <- 3
        }else{
          severity <- 2
        }
      }else{
        severity <- 1
      }
      
      # adjust the probability that someone seeks testing
      if(vax==1){
        if(severity==1){prob_test <- prob_test_1_vax}
        if(severity==2){prob_test <- prob_test_2_vax}
        if(severity==3){prob_test <- prob_test_3_vax}
      }else{
        if(severity==1){prob_test <- prob_test_1_unvax}
        if(severity==2){prob_test <- prob_test_2_unvax}
        if(severity==3){prob_test <- prob_test_3_unvax}        
      }
      # store the result only if they seek testing
      if(rbinom(1,1,prob_test)){
        df <- df %>% add_row(id=i, vax=vax, age=age, time=round(Ui), 
                             severity=severity,result = 0, keep = 1)
      } else {
        df <- df %>% add_row(id=i, vax=vax, age=age, time=round(Ui), 
                             severity=severity,result = 0, keep = 0)
      }      
      
      last_tn <- round(Ui)
      
    }else{
      z=FALSE # otherwise, stop the loop
    }
  }
}