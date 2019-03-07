#!/usr/bin/env Rscript
# Author: Eva Linehan
# Date: February 2019
# Desc: Fitting models for Miniproject

#clear environments
rm(list=ls())
graphics.off()
library(minpack.lm)
library(R.utils)
library(ggplot2)
library(plyr)


initialtime = proc.time()[3] # Start timer to time script execution
DF <- read.csv("../Data/Ready_to_fit.csv", header = TRUE) # Read in and observe data
#values <- subset(DF, DF$FinalID == "MTD3375" , select = c(OriginalTraitValue, ConTemp, Kelvin))



############################################## Model Functions ######################################################


############################################## Briere's Model #######################################################



briere<- function(Temp, T0, Tm, c){
  # T = temperature values in TPC
  # T0 = lower thermal limit
  # Tl = upper thermal limit
  # c = constant
  return(c*Temp*(Temp-T0)*(abs(Tm-Temp)^(1/2))*as.numeric(Temp<Tm)*as.numeric(Temp>T0))
  # Function forced to be 0 above and below max and min temp
}



############################################# Schoofield Models #####################################################


# Global Boltzmann constant 
assign("k", 8.617 * 10^-5, envir = .GlobalEnv) 


# k    : Boltzmann constant
# B0   : Normalization constant (trait value), initialized as the trait value that corresponds to the closest recorded temperature at 283.15 K.
# E    : Activation energy (slope right hand side of curve)
# E_h  : High temp deactivation energy, i.e peak trait value achieved (slope left hand side of curve)
# T_h  : Temperture at which 50% high temp deactivation occurs
# E_l  : Low temp deactivation energy (slope right hand side of the curve)
# T_l  : Temperature estimated at which the trait value was halved from the peak (50% low temp deactivation occurs) (linear regression of right hand side of peak)
# Temp : Vector of temperatures in K


# 278.15K (5 degrees C).taken as reference temperature



# Equation 7 in Schoolfield et al.,1991 (no low-temp inactivation)


Linear_Schoolnolow <- function(B0, E, E_h, T_h, Temp){
  
  return(B0 * exp((-E/k)*((1/Temp) - (1/278.15)))/(1  + exp(E_h/k * ((1/T_h) - (1 / Temp)))))
  
}

# Equation 6 in Schoolfield et al.,1991 (no high-temp interaction)
Linear_Schoolnohigh <- function(B0, E, E_l, T_l, Temp){
  
  return(B0 * exp((-E/k)*((1/Temp) - (1/278.15)))/(1 + exp(E_l/k * ((1/T_l) - (1 / Temp)))))
  
}


# Equation 4 in Schoolfield et al.,1991 (full model)

Linear_Schoolfull<- function(B0, E_h, E_l, E, T_l, T_h, Temp){
  
  
  return((B0 * exp((-E/k)*((1/Temp) - (1/278.15))))/(1 + exp((E_l/k) * ((1/T_l) - (1 / Temp))) + exp((E_h/k) * ((1/T_h) - (1 / Temp)))))
}




#################### Estimating Schoolfield starting parameters #######################


Startvals <- function(values){
  # Calculate starting values from log-linear tranformation of below and above peak values (including temperature and trait values)
  Peak <- max(values[which(values[, "OriginalTraitValue"] == max(values[, "OriginalTraitValue"])), "Kelvin"]) # Temperature of highest trait value
  
  # Identify left and right had side of the curve
  
  # Left handside
  tryCatch({Below_peak <<- subset(values, values[,"Kelvin"] <= Peak)
  },
  error = function(e){
    # If there are no values below the peak, initialize the following
    E <<- -0.65
    T_l <<- min(values$Kelvin)
  })
  # Right hand side 
  tryCatch({Above_peak <<- subset(values, values[,"Kelvin"] >= Peak)
  },
  error = function(e){
    # If there are no values above the peak, initialize the following
    T_h <<- max(values$Kelvin)
    E_h <<- abs(E) *3
  })
  
 # For below peak subsets containing more than 1 value, perform log-linear tranformation
  if(length(Below_peak$OriginalTraitValue)>1){
    x <<- 1/(k * Below_peak$Kelvin)
    y <<- log(Below_peak$OriginalTraitValue)
    model <<- lm(y~x)
    fit <<- summary(model)
    E <<- fit$coefficients[2] # Initialize E as slope
    Half_val <<- log(mean(Below_peak$OriginalTraitValue)) # Find mean log Trait value on y axis
    One_over_kT_l <<- ((Half_val - coef(fit)[1])/coef(fit)[2]) # Calculate the 1/kt value (x) at mean log Trait value (y) according to linear model
    T_l <<- 1/(One_over_kT_l*k) # Convert 1/kt value to Kelvin and initialize as temp at 50% enzyme low deactivation
  } else {E <- -0.65; T_l <- min(Below_peak$Kelvin)}
  
  # For above peak subsets containing more than 1 value, perform log-linear tranformation
  if(length(Above_peak$OriginalTraitValue)>1){
    x <- 1/(k * Above_peak$Kelvin)
    y <- log(Above_peak$OriginalTraitValue)
    model2<-lm(y~x)
    fit <- summary(model2)
    E_h <- fit$coefficients[2]  # Initialize E_h as slope
    Half_val <- log(mean(Above_peak$OriginalTraitValue))
    One_over_kT_h <- ((Half_val - coef(fit)[1])/coef(fit)[2]) # Calculate the 1/kt value (x) at mean log Trait value (y) according to linear model
    T_h <- 1/(One_over_kT_h*k) # Convert 1/kt value to Kelvin and initialize as temp at 50% enzyme high deactivation
  } else {T_h<- max(Above_peak$Kelvin); E_h <- abs(E) *3}
  
  
  B0 <- with(values, if(which(abs(Kelvin-278.15) == min(abs(Kelvin-278.15)))) return(OriginalTraitValue[which.min(abs(Kelvin-278.15))]))
  # Trait value closest to 283.15 K
  
  E_l <- E * 0.5 # Low temp deactivation energy
  
  return(list(B0 = B0, E_h = E_h, E_l = abs(E_l), E = abs(E), T_l = T_l, T_h = T_h, Peak = Peak))
}
  
#Startvals(values)


#################################################################### Fitting ########################################################################


Gaussian_fluctuation<- function(param_calc){
  # Randomly sample new parameter estimate around a normal distribution
  mu = param_calc # Mean set as calculated parameter estimate
  new_param <- rnorm(1, mean = mu, sd = 0.5)
  return(new_param)
}

AICc <- function(aic, K, n){
  # Calculate AICc for fitted model
  val <- pmax(K,n) # Select the highest value among k and n
  result <- aic + ((2*K*(K+1)*(K+2))/max(val+3) - K - 2)
  # n is sample size
  # K is number of model parameters
  return(result)
}



cubic <- function(values){
  # Fit cubic model
  K <- 4 # Number of parameters in model
  n <- length(values$OriginalTraitValue) # Sample size
  fit <- lm(values$OriginalTraitValue ~ poly(values$ConTemp, 3, raw = TRUE))
  aicc <- AICc(AIC(fit), K, n) # Calculate AICcc
  return(cAICc = aicc)
}
#cubic(values)


try_briere <- function(values){
  # Fit Briere model
  K <- 3 # Number of parameters for AICC calculation
  n <- length(values$OriginalTraitValue) # Sample size
  tryCatch({fit <<- nlsLM(OriginalTraitValue ~ briere(ConTemp, T0, Tm, c), start = list(T0=min(values$ConTemp), Tm=max(values$ConTemp), c=1), 
                          data = values, control = nls.control(maxiter = 1000), 
                          lower = c(0,0,0), upper = c(20, 40, Inf))
  aicc <<- AICc(AIC(fit), K, n) # Calculate AICc of fitted model
  return(list(T0 = coef(fit)["T0"], Tm = coef(fit)["Tm"], c = coef(fit)["c"], bAICc = aicc)) # Return best model starting parameters and AICc
  rm(fit, aicc) # Prevent variables declared in trycatch to be left in global environment
  },
  error = function(e){ 
    # Return NAs if model does not converge
    return(list(T0 = NA, Tm = NA, c = NA, bAICc = NA)) 
  })
}

#try_briere(values)


try_lin_schoolnl <- function(values){ 
  # Fit simplified Schoolfied model with high temp deactivation, S1
  param_calc <- Startvals(values) # Calculate starting parameters
  params_required <- c(param_calc[1], param_calc[4], param_calc[2], param_calc[6]) # Append required starting estimates to vector
  Peak <- as.numeric(param_calc[7]) # Peak to bound T_h
  counter = 0
  error_counter = 0
  lowest_aicc <- 1000
  K <- 4 # Number of parameters for aicc function
  n <- length(values$OriginalTraitValue) # sample size, n, for aicc function
  while (counter < 20) {
    counter <- counter + 1
    if(counter == 1){
      tryCatch({mod <<- nlsLM(OriginalTraitValue ~ Linear_Schoolnolow(B0, E, E_h, T_h, values$Kelvin), 
                              start = list(B0 = params_required[1], E = params_required[2], E_h = params_required[3], 
                                           T_h = params_required[4]), data = values, control = nls.control(maxiter = 1000),
                              lower = c(0,0,0,Peak),  upper = c(Inf, 3, Inf, 273.15 + 45))
      new_aicc <<- AICc(AIC(mod), K, n )}, # Calculate AICc of model
      error = function(e){
        error_counter <<- error_counter + 1
        new_aicc <<- 1005 # If error set an AICc value higher than lowest AICc
      })}
    else{
      parameters <- lapply(params_required, Gaussian_fluctuation) # Randomize parameters
      start_val <- c(parameters[1], parameters[2], parameters[3], parameters[4]) # Append new starting estimates to vector
      tryCatch({mod <<- nlsLM(OriginalTraitValue ~ Linear_Schoolnolow(B0, E, E_h, T_h, values$Kelvin), 
                              start = list(B0 = start_val[1], E = start_val[2], E_h = start_val[3], 
                                           T_h = start_val[4]), data = values, control = nls.control(maxiter = 1000),
                              lower = c(0,0,0,Peak),  upper = c(Inf, 3, Inf, 273.15 + 45))
      new_aicc <<- AICc(AIC(mod), K, n)},  # Calculate AICc of model
      error = function(e){
        error_counter <<- error_counter + 1
        new_aicc <<- 1005 # If error set an AICc value higher than lowest AICc
      }
      )}
    if (new_aicc < lowest_aicc){
      # Compare each new model AICc to the lowest AICc
      lowest_aicc <- new_aicc # Set as new lowest AICc
      best_mod1 <- mod # Include model 
    }
  }
  Convergence_prop <- 1 - (sum(error_counter)/20) # Proportion of times models converged 
  while(counter == 20){
    if(lowest_aicc == 1000){
      # If no models converged return NA values
      return(list(B0 = NA, E = NA, E_h = NA, T_h = NA, Lin_Snl_AICc = NA, Lin_schoolnl_convergence = Convergence_prop))
      rm(mod, new_aicc, error_counter) # Prevent variables declared in TryCatch hanging in global environment
    } 
    else{
      # Return best model starting parameters, AICc and convergence proportion
      return(list(B0 = summary(best_mod1)$coef[1,1], E = summary(best_mod1)$coef[2,1], E_h = summary(best_mod1)$coef[3,1], T_h = summary(best_mod1)$coef[4,1], Lin_Snl_AICc = lowest_aicc, Lin_schoolnl_convergence = Convergence_prop))
      rm(mod, new_aicc, error_counter) # Prevent variables declared in TryCatch hanging in global environment
    }
  } 
}
#try_lin_schoolnl(values)


try_lin_schoolnh <- function(values){
  # Fit simplified Schoolfied model with low temp deactivation, S2
  param_calc <- Startvals(values)  # Calculate starting parameters
  params_required <- c(param_calc[1], param_calc[4], param_calc[3], param_calc[5]) # Append required starting estimates to vector
  Peak <- as.numeric(param_calc[7]) # Peak to bound T_l
  counter = 0
  error_counter = 0
  K <- 4 # Number of parameters for aicc function
  n <- length(values$OriginalTraitValue) # sample size, n, for aicc function
  Lowest_aicc <- 1000
  while (counter < 20) {
    counter <- counter + 1
    if(counter == 1){
      tryCatch({mod <<- nlsLM(OriginalTraitValue ~ Linear_Schoolnohigh(B0, E, E_l, T_l, values$Kelvin), 
                              start = list(B0 = params_required[1], E = params_required[2], E_l = params_required[3], 
                                           T_l = params_required[4]), data = values, control = nls.control(maxiter = 1000),
                              lower = c(0,0,0,0), upper = c(Inf, 3, E, Peak))
      new_aicc <<- AICc(AIC(mod), K, n)}, # Calculate AICc of model
      error = function(e){
        error_counter <<- error_counter + 1
        new_aicc <<- 1005 # If error set an AICc value higher than lowest AICc
      })}
    else{
      parameters <- lapply(params_required, Gaussian_fluctuation) # Randomize parameter
      start_val <- c(parameters[1], parameters[2], parameters[3], parameters[4]) # Append new starting estimates to vector
      tryCatch({mod <<- nlsLM(OriginalTraitValue ~ Linear_Schoolnohigh(B0, E, E_l, T_l, values$Kelvin), 
                              start = list(B0 = start_val[1], E = start_val[2], E_l = start_val[3], 
                                           T_l = start_val[4]), data = values, control = nls.control(maxiter = 1000),
                              lower = c(0,0,0,0), upper = c(Inf, 3, E, Peak))
      new_aicc <<- AICc(AIC(mod), K, n)}, # Calculate AICc of fitted model
      error = function(e){
        error_counter <<- error_counter + 1
        new_aicc <<- 1005} # If error set an AICc value higher than lowest AICc
      )}
    if (new_aicc < Lowest_aicc){
      # Compare each new model AICc to the lowest AICc
      Lowest_aicc <- new_aicc # Set as new lowest AICc
      best_mod2 <- mod # Include model 
    } 
  }
  Convergence_prop <- 1 - (sum(error_counter)/20) # Proportion of times models converged 
  while(counter == 20){
    if(Lowest_aicc == 1000){
      # If no models converged return NA values
      return(list(B0 = NA, E = NA, E_l = NA, T_l = NA, Lin_Snh_AICc = NA, Lin_schoolnh_convergence = Convergence_prop))
      rm(mod, new_aicc, error_counter) # Prevent variables declared in TryCatch hanging in global environment
    } 
    else{
      # Return best model starting parameters, AICc and convergence proportion
      return(list(B0 = summary(best_mod2)$coef[1,1], E = summary(best_mod2)$coef[2,1], E_l = summary(best_mod2)$coef[3,1], T_l = summary(best_mod2)$coef[4,1], Lin_Snh_AICc = Lowest_aicc, Lin_schoolnh_convergence = Convergence_prop))
      rm(mod, new_aicc, error_counter) # Prevent variables declared in TryCatch hanging in global environment
    }
  } 
}
#try_lin_schoolnh(values)




try_lin_schoolfull <- function(values){
  # Fit full Schoolfied model with low temp  and high temp deactivation, S3
  param_calc <- Startvals(values) # Calculate starting parameters
  params_required <- c(param_calc[1], param_calc[2], param_calc[3], param_calc[4], param_calc[5], param_calc[6]) # Append required starting estimates to vector
  counter = 0
  Peak = as.numeric(param_calc[7]) # Peak to bound T_l and T_h
  error_counter = 0
  K <- 6 # Number of parameters for aicc function
  n <- length(values$OriginalTraitValue)  # sample size, n, for aicc function
  Lowest_aicc <- 1000
  while (counter < 20) {
    counter <- counter + 1
    if(counter == 1){
      tryCatch({mod1 <<- nlsLM(OriginalTraitValue ~ Linear_Schoolfull(B0, E_h, E_l, E, T_l, T_h, values$Kelvin),
                               start = list(B0 = params_required[1], E_h = params_required[2], E_l = params_required[3], 
                                            E = params_required[4], T_l = params_required[5], T_h = params_required[6]), 
                               data = values, control = nls.control(maxiter = 1000),
                               lower = c(0,0,0,0,273.15,Peak), upper = c(Inf, Inf, E , 3, Peak, 273.15 + 45))
      new_aicc <<- AICc(AIC(mod1), K, n)}, # Calculate AICc of fitted model
      error = function(e){
        error_counter <<- error_counter + 1
        new_aicc <<- 1005 # If error set an AICc value higher than lowest AICc
      },
      warning = function(w){
        #print(paste("NaNs produced with these parameters in Final ID loop", counter))
      })}
    else{
      parameters <- lapply(params_required, Gaussian_fluctuation)
      start_val <- c(parameters[1], parameters[2], parameters[3], parameters[4], parameters[5], param_calc[6])
      tryCatch({mod1 <<- nlsLM(OriginalTraitValue ~ Linear_Schoolfull(B0, E_h, E_l, E, T_l, T_h, values$Kelvin),
                               start = list(B0 = start_val[1], E_h = start_val[2], E_l = start_val[3], 
                                            E = start_val[4], T_l = start_val[5], T_h = start_val[6]), 
                               data = values, control = nls.control(maxiter = 1000),
                               lower = c(0,0,0,0,273.15,Peak), upper = c(Inf, Inf, E , 3, Peak, 273.15 + 45))
      new_aicc <<- AICc(AIC(mod1), K, n)}, # Calculate AICc of fitted model
      error = function(e){
        error_counter <<- error_counter + 1
        new_aicc <<- 1005}, # If error set an AICc value higher than lowest AICc
      warning = function(w){
        #print(paste("NaNs produced with these parameters in Final ID loop", counter))
      }
      )}
    if (new_aicc < Lowest_aicc){
      # Compare each new model AICc to the lowest AICc
      Lowest_aicc <- new_aicc # Set as new lowest AICc
      best_mod3 <- mod1 # Include model 
    } 
  }
  Convergence_prop <- 1 - (sum(error_counter)/20) # Proportion of times models converged 
  while(counter == 20){
    if(Lowest_aicc == 1000){
      # If no models converged return NA values
      return(list(B0 = NA, E_h = NA, E_l = NA, E = NA, T_l = NA, T_h = NA, Lin_Sfull_AICc = NA, Lin_Sfull_convergence = Convergence_prop))
      rm(mod1, new_aicc, error_counter)  # Prevent variables declared in TryCatch hanging in global environment
    } 
    else{
      # Return best model starting parameters, AICc and convergence proportion
      return(list(B0 = summary(best_mod3)$coef[1,1], E_h = summary(best_mod3)$coef[2,1], E_l = summary(best_mod3)$coef[3,1], E = summary(best_mod3)$coef[4,1], T_l = summary(best_mod3)$coef[5,1], T_h = summary(best_mod3)$coef[6, 1], Lin_Sfull_AICc = Lowest_aicc, Lin_Sfull_convergence = Convergence_prop))
      rm(mod1, new_aicc, error_counter)  # Prevent variables declared in TryCatch hanging in global environment
    }
  } 
}

#try_lin_schoolfull(values)


###################################### Fit Models and Output csv ###########################################


IDlist = unique(DF$FinalID) # List of individual IDs
IDs <- data.frame(NULL) # Pre-allocate a list to append each individual ID to add to final table
Mod_fit_table <- matrix(ncol = 26, nrow = length(IDlist)) # Pre-allocate output matrix


#initialtime = proc.time()[3]
for (i in (1:length(IDlist))){
  subset_id<- IDlist[i] # Take out the ID for each iteration
  IDs[i,1] <- subset_id # Append ID to list
  values = subset(DF, DF$FinalID == subset_id, select = c("OriginalTraitValue", "Kelvin", "ConTemp")) # Subset out_table by ID in unique list
  # Run cubic
  Mod_fit_table[i,1] <- cubic(values)
  # Run briere and append to matrix
  model2<- try_briere(values)
  Mod_fit_table[i,2] <- model2$T0
  Mod_fit_table[i,3] <- model2$Tm
  Mod_fit_table[i,4] <- model2$c
  Mod_fit_table[i,5] <- model2$bAICc
  # Run Schoolfield models and append to matrix
  # Schoolnolow (S1)
  model3 <- try_lin_schoolnl(values)
  Mod_fit_table[i,6] <- model3$B0
  Mod_fit_table[i,7] <- model3$E
  Mod_fit_table[i,8] <- model3$E_h
  Mod_fit_table[i,9] <- model3$T_h
  Mod_fit_table[i,10] <- model3$Lin_Snl_AICc
  Mod_fit_table[i,11] <- model3$Lin_schoolnl_convergence
  #Schoolnohigh (S2)
  model4 <- try_lin_schoolnh(values)
  Mod_fit_table[i,12] <- model4$B0
  Mod_fit_table[i,13] <- model4$E
  Mod_fit_table[i,14] <- model4$E_l
  Mod_fit_table[i,15] <- model4$T_l
  Mod_fit_table[i,16] <- model4$Lin_Snh_AICc
  Mod_fit_table[i,17] <- model4$Lin_schoolnh_convergence
  #Schoolfull (S3)
  model5 <- try_lin_schoolfull(values)
  Mod_fit_table[i,18] <- model5$B0
  Mod_fit_table[i,19] <- model5$E_h
  Mod_fit_table[i,20] <- model5$E_l
  Mod_fit_table[i,21] <- model5$E
  Mod_fit_table[i,22] <- model5$T_l
  Mod_fit_table[i,23] <- model5$T_h
  Mod_fit_table[i,24] <- model5$Lin_Sfull_AICc
  Mod_fit_table[i,25] <- model5$Lin_Sfull_convergence
  # Number of datapoints
  Mod_fit_table[i,26] <- length(values$OriginalTraitValue)
}

#finish_time = proc.time()[3] - initialtime  

#Mod_fit <- list()

#call_models <- function (values){
#  output = list()
#  output[1] <- list(cubic(values))
#  output[2]<- list(try_briere(values))
#  output[3] <- list(try_lin_schoolnl(values))
#  output[4] <- list(try_lin_schoolnh(values))
#  output[5] <- list(try_lin_schoolfull(values))
  
#  return(output)
 
#}
#call_models(values)

#finallist <- list()

#initialtime = proc.time()[3]
#sapply(IDlist, function(id) {
#  values <- subset(DF, DF$FinalID == id, select = c("OriginalTraitValue", "Kelvin", "ConTemp")) # Subset out_table by ID in unique list
#  row_out <- call_models(values)
#  out <- list.append(finallist, row_out)
#  do.call(rbind, out)
#  }
#)


#finish_time = proc.time()[3] - initialtime  
  
  
  
#output_table <- as.data.frame(t(as.data.frame(finallist))) 


Mod_fit_table <- as.data.frame(Mod_fit_table)
colnames(Mod_fit_table) <- c("C_AICC", "T0", "Tm", "c", "bAICC",
                         "Snl_B0", "Snl_E", "Snl_E_h", "Snl_T_h", "Snl_AICC", "Snl_Converge",
                         "Snh_B0", "Snh_E", "Snh_E_l", "Snh_T_l", "Snh_AICC", "Snh_Converge",
                         "Full_B0", "Full_E_h", "Full_E_l", "Full_E", "Full_T_l", "Full_T_h", "Full_AICC", "Full_Converge", "Sample_size")
Out_table<- cbind(IDs, Mod_fit_table) # Append infividual IDs to final dataframe
write.csv(Out_table, file = "../Data/Model_Results.csv")
finish_time = (proc.time()[3] - initialtime)/60  
print(paste("Script took", finish_time, "minutes to run"))
