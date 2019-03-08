#!/usr/bin/env Rscript
# Author: Eva Linehan
# Date: February 2019
# Desc: PLotting TPC curves

#clear environments
rm(list=ls())
graphics.off()
library(ggplot2)
library(minpack.lm)
library(plyr)
suppressMessages(library(dplyr))
suppressMessages(library(gridExtra))
library(xtable)
library(RColorBrewer)
library(stringr)
suppressMessages(library(Hmisc))
library(reshape2)
library(compiler)

enableJIT(3)


# Temporarily turn off warnings
options(warn = -1)

initialtime = proc.time()[3]

Data <- read.csv("../Data/Ready_to_fit.csv", header = TRUE) # Read in wrangled Biotraits data
Models <- read.csv("../Data/Model_Results.csv", header = TRUE) # Read in output from model fitting


########################################### Model Funtions #####################################################


# Global Boltzmann constant 
assign("k", 8.617 * 10^-5, envir = .GlobalEnv) 


briere<- function(Temp, T0, Tm, c){
  return(c*Temp*(Temp-T0)*(abs(Tm-Temp)^(1/2))*as.numeric(Temp<Tm)*as.numeric(Temp>T0))
  # Function forced to be 0 above and below max and min temp
}

briere <- cmpfun(briere)

# Equation 7 in Schoolfield et al.,1991 (no low-temp inactivation)


Linear_Schoolnolow <- function(B0, E, E_h, T_h, Temp){
  
  return(B0 * exp((-E/k)*((1/Temp) - (1/278.15)))/(1  + exp(E_h/k * ((1/T_h) - (1 / Temp)))))
  
}

Linear_Schoolnolow <- cmpfun(Linear_Schoolnolow)


# Equation 6 in Schoolfield et al.,1991 (no high-temp interaction)

Linear_Schoolnohigh <- function(B0, E, E_l, T_l, Temp){
  
  return(B0 * exp((-E/k)*((1/Temp) - (1/278.15)))/(1 + exp(E_l/k * ((1/T_l) - (1 / Temp)))))
  
}

Linear_Schoolnohigh <- cmpfun(Linear_Schoolnohigh)


# Equation 4 in Schoolfield et al.,1991 (full model)

Linear_Schoolfull<- function(B0, E_h, E_l, E, T_l, T_h, Temp){
  
  
  return((B0 * exp((-E/k)*((1/Temp) - (1/278.15))))/(1 + exp((E_l/k) * ((1/T_l) - (1 / Temp))) + exp((E_h/k) * ((1/T_h) - (1 / Temp)))))
}


Linear_Schoolfull <- cmpfun(Linear_Schoolfull)


################################################## Plot fits #######################################################



plot_cubic <- function(traits_data){
  # Predict trait values according to the cubic model from sequence temperature data
  new_data_c <- seq(ceiling(min(traits_data$ConTemp)), ceiling(max(traits_data$ConTemp)), length.out = length(traits_data$OriginalTraitValue))
  cube_fit <- lm(traits_data$OriginalTraitValue ~ poly(traits_data$ConTemp, 3, raw = TRUE)) # Fit sequence data to model
  Model_data_c <- data.frame(TraitValue = predict(cube_fit), Temperature = new_data_c) # Output dataframe of temperatures and predicted trait values
  return(Model_data_c)
}

plot_cubic <- cmpfun(plot_cubic)


plot_briere <- function(best_params, new_data){
  # Predict trait values according to the briere model from sequence temperature data
  pred.b <- briere(new_data, best_params$T0, best_params$Tm, best_params$c) # Fit sequence data to model
  Model_data_b <- data.frame(TraitValue = pred.b, Temperature = new_data) # Output dataframe of temperatures and predicted trait values
  return(Model_data_b)
}

plot_briere <- cmpfun(plot_briere)


plot_school_nl <- function(best_params, new_data, new_data_school){
  # Predict trait values according to the schoolfield high temp deactivation model (S1) from sequence temperature data
  snl.fit <- Linear_Schoolnolow(best_params$Snl_B0, best_params$Snl_E, best_params$Snl_E_h, best_params$Snl_T_h, new_data_school) # Fit sequence data to model
  Model_data_snl <- data.frame(TraitValue = snl.fit, Temperature = new_data) # Output dataframe of temperatures and predicted trait values
  return(Model_data_snl)
}

plot_school_nl <- cmpfun(plot_school_nl)


plot_school_nh <- function(best_params, new_data, new_data_school){
  # Predict trait values according to the schoolfield low temp deactivation model (S2) from sequence temperature data
  pred.snh <- Linear_Schoolnohigh(best_params$Snh_B0, best_params$Snh_E, best_params$Snh_E_l, best_params$Snh_T_l, new_data_school) # Fit sequence data to model
  Model_data_snh <- data.frame(TraitValue = pred.snh, Temperature = new_data) # Output dataframe of temperatures and predicted trait values
  return(Model_data_snh)
}

plot_school_nh <- cmpfun(plot_school_nh)


plot_school_full<- function(best_params, new_data, new_data_school){
  # Predict trait values according to the full schoolfield model (S3) from sequence temperature data
  pred.full <- Linear_Schoolfull(best_params$Full_B0, best_params$Full_E_h, best_params$Full_E_l, best_params$Full_E, best_params$Full_T_l, best_params$Full_T_h, new_data_school) # Fit sequence data to model
  Model_data_full <- data.frame(TraitValue = pred.full, Temperature = new_data) # Output dataframe of temperatures and predicted trait values
  return(Model_data_full)
}

plot_school_full <- cmpfun(plot_school_full)


Build_plot <- function(id){
  # Function to plot all models on original dataset for input ID 
  traits_data <- subset(Data, Data$FinalID == id, select = c("ConTemp", "OriginalTraitValue", "Kelvin")) # Original trait data for ID
  best_params <- Models [i,] # Best parameter data for ID
  new_data<- seq(ceiling(min(traits_data$ConTemp)), ceiling(max(traits_data$ConTemp)), length.out = 100) # Temperature data to be predicted in Celsius
  new_data_school<- seq(ceiling(min(traits_data$Kelvin)), ceiling(max(traits_data$Kelvin)), length.out = 100) # Temperature data to be predicted in Kelvin
  # Plot
  fit_plot <- ggplot(traits_data, aes(x = ConTemp, y = OriginalTraitValue)) + 
          xlab(expression(paste("Temperature (", degree, C, ")"))) + 
          ylab("Trait Value") + xlim (min(traits_data$ConTemp), max(traits_data$ConTemp)) +
          geom_line(data = plot_cubic(traits_data), aes(x = Temperature, y = TraitValue, colour = "deepskyblue2"), lwd = 1) +
          geom_line(data = plot_briere(best_params, new_data), aes(x = Temperature, y = TraitValue, colour = "chartreuse2"), lwd = 1) +
          geom_line(data = plot_school_nl(best_params, new_data, new_data_school), aes(x = Temperature, y = TraitValue, colour = "darkorchid2"), lwd = 1) +
          geom_line(data = plot_school_nh(best_params, new_data, new_data_school), aes(x = Temperature, y = TraitValue, colour = "gold1"), lwd = 1) +
          geom_line(data = plot_school_full(best_params, new_data, new_data_school), aes(x = Temperature, y = TraitValue, colour = "brown2"), lwd = 1) +
          geom_point(size = 3, col = "black", alpha = 0.9, pch = 20) +
          scale_color_identity(name = "Model fit",
                               breaks = c("deepskyblue2", "chartreuse2", "darkorchid2", "gold1", "brown2" ),
                               labels = c("Cubic", "Briere", "S1", "S2", "S3" ),guide = "legend") +
          theme_light(base_size = 16)
  return(fit_plot) # Output plot created
}

Build_plot <- cmpfun(Build_plot)


IDlist = unique(Models$V1) # Unique list of IDs
for (i in (1:length(IDlist))){
  #Plot created using Build_plot fuction and saved in results folder for each ID
  ID <- IDlist[i]
  pdf(file = paste("../Results/Models_fitted/", ID,".pdf", sep = ""))
  print(Build_plot(ID))
  dev.off()
}
  
  

##################################################### Analysis #############################################################


######### Range of sample sizes per ID in dataset ##########

Observations <- as.data.frame(count(Data, FinalID)) # Number of observations per ID
#quantile(Observations[,2], c(.25, .50, .75 ))
#min(Observations[,2])
#max(Observations[,2])
#median(Observations[,2])

df <- Observations %>% group_by(n) %>% count() # Tally the number of observations for each ID


# Create and save a boxplot to illustrate range of sample sizes for each TPC
pdf(file = paste("../Results/Samples_boxplot.pdf"))
print(ggplot(df, aes(x = n, y = nn, group = 1)) + xlab("Frequency") + ylab("Sample Size") + 
  geom_boxplot()+ 
  geom_text(data = subset(df, nn > 131), 
            mapping = aes(label = nn),
            nudge_x = 240, size = 5, check_overlap = T) +
    theme_minimal(base_size = 15) +
  theme_light())
invisible(dev.off())


######### Starting params and AICc for Fig. 1 ################
# Figure 1 generated in LaTeX

best_paramsA <- subset(Models, Models$V1 == "MTD2100") # Model starting values for ID in subfigure 1
A <- best_paramsA[,grepl("AICC", names(best_paramsA))] # Find all columns with AICC values
best_paramsB <- subset(Models, Models$V1 == "MTD2170") # Model starting values for ID in subfigure 2
B <- best_paramsB[,grepl("AICC", names(best_paramsB))] # Find all columns with AICC values
fit_table2 <- data.frame(rbind(as.numeric(A), as.numeric(B))) # Generate table of AICC values for each model
colnames(fit_table2)<- c("Cubic", "Briere", "S1", "S2", "S3")
TPC <- c("a) Green Seaweed", "b) Proteobacteria")
fit_table2 <- cbind(TPC, fit_table2)

# Output created table as .tex file
table2 <- xtable(fit_table2)
align(table2) <- "cc|ccccc"
print.xtable(table2, type = "latex", floating = F, caption.placement = "top", file = "../Results/Table2.tex",
             include.rownames = F)



################## Model Performance #######################

############## Summary of fits  ################

aicc <- as.matrix(Models[,grepl("AICC", names(Models))]) # Create a matrix of all aicc values
converge <- as.matrix(Models[,grepl("Converge", names(Models))]) # Seperate convergence proportions

Number <- apply(aicc, 2, function(x) 1577 - sum(is.na(x))) # Calculate number of model fits for the entire data set
Proportion <- apply(aicc, 2, function(x) paste(format(round(((1577 - sum(is.na(x)))/1577)*100), digits = 2), "%", sep = "")) # Calculate percentage proportion
Conv <- apply(converge, 2, function(x) paste(format(round(mean(x)*100), digits = 2), "%", sep = "")) # Calculate average convergence rate for each schoolfield model
Convergence <- c("-", "-", Conv[[1]], Conv[[2]], Conv[[3]]) # Convergence row, dashes indicating both cubic and briere which were fit once
fit_table <- as.matrix(rbind(Number, Proportion, Convergence))
colnames(fit_table)<- c("Cubic", "Briere", "S1", "S2", "S3")


# Output created table as .tex file
table1 <- xtable(fit_table, align = "l|ccccc")
print.xtable(table1, type = "latex", floating = F, caption.placement = "top", 
             file = "../Results/Table1.tex", scalebox = '1.5')


################ Delta AICC ####################

min_val <- apply(aicc, 1, function(x) min(x, na.rm = T)) # For each TPC, select model with lowest AICC
delta_mat <- matrix(ncol = 5, nrow = length(min_val)) # Matrix to fill with delta AICC values
for(i in 1:length(min_val)){
  # Create matric converting AICCs into delta AICCs by comparing distance with the best model
 delta_mat[i,1] <- abs(min_val[i] - aicc[i,1])
 delta_mat[i,2] <- abs(min_val[i] - aicc[i,2])
 delta_mat[i,3] <- abs(min_val[i] - aicc[i,3])
 delta_mat[i,4] <- abs(min_val[i] - aicc[i,4])
 delta_mat[i,5] <- abs(min_val[i] - aicc[i,5])
}

# Categorize Delta AICCs based on guidelines
less_2 <- apply(delta_mat, 2, function(x) sum(x <=2, na.rm = T)) # Accumulate those less than equal to 2
two_four <- apply(delta_mat, 2, function(x) sum(x > 2 & x <= 4, na.rm = T)) # Accumulate those greater than 2, less than equal to 4
four_seven <- apply(delta_mat, 2, function(x) sum(x > 4 & x <= 7, na.rm = T)) # Accumulate those greater than 4, less than equal to 7
seven_ten <- apply(delta_mat, 2, function(x) sum(x > 7 & x <= 10, na.rm = T)) # Accumulate those greater than 7, less than equal to 10
ten <- apply(delta_mat, 2, function(x) sum(x > 10, na.rm = T)) # Accumulate those greater than 2
Model <- c("Cubic", "Briere", "S1", "S2", "S3") 
delta_table <- data.frame(cbind(Model,less_2, two_four, four_seven, seven_ten, ten)) # Generate table
colnames(delta_table) <- c("Model", "$\\Delta < 2$", "$2 < \\Delta \\leq 4$", "$4 < \\Delta \\leq 7$", "$ 7 < \\Delta \\leq 10$", "$\\Delta > 10$") # Column names in LaTeX syntax

# Output created table as .tex file
table3 <- xtable(delta_table, type ="latex")
align(table3) <- "lc|ccccc"
print(table3, floating = F, caption.placement = "top", file = "../Results/Table3.tex",
             include.rownames = F, scalebox = '0.6', sanitize.text.function = function(x){x}, width = 6)

################ AICC Weightage #####################

delta_weight <- na.omit(delta_mat) # Omit NAs from weightage analysis
  
Weight1<- function(id, col){
# Calculate weight for each model in row i
  return(exp(-(delta_weight[id,col])/2))
}

Weight2 <- function(weight, sum_weight){
# Calculate relative weight among models
  return(weight/sum_weight)
}
weightages <- matrix(ncol = 5, nrow = length(delta_weight[,1])) # Weightage matrix
for(i in 1:length(delta_weight[,1])){
  # Calculate weightage of each AICC
  # Calculate independent weight for each model
  c <- Weight1(i, 1)
  b <- Weight1(i, 2)
  snl <- Weight1(i, 3)
  snh <- Weight1(i, 4)
  sf <- Weight1(i, 5)
  sum_weight <- sum(c+b+snl+snh+sf) # Sum weightages
  # Calculate proportional/ actual weights
  weightages[i,1] <- Weight2(c, sum_weight) 
  weightages[i,2] <- Weight2(b, sum_weight) 
  weightages[i,3] <- Weight2(snl, sum_weight) 
  weightages[i,4] <- Weight2(snh, sum_weight) 
  weightages[i,5] <- Weight2(sf, sum_weight) 
}

weightages <- as.data.frame(weightages)
colnames(weightages) <- c("Cubic", "Briere", "S1", "S2", "S3")
melt_weightages <- melt(weightages, measure.vars = c("Cubic", "Briere", "S1", "S2", "S3")) # Melt dataframe for figure


# Boxplot of AICC distributions
pdf(file = paste("../Results/Fig_5.pdf"))
print(ggplot(melt_weightages, aes(factor(variable), value, color = variable)) + 
  geom_boxplot() +
  xlab("Models") + ylab ("Akaike Weight") +
  scale_color_brewer(palette="Dark2") +
  theme_minimal(base_size = 15))
invisible(dev.off())


########## Model performance bar plots ##############

best_mod <- apply(aicc, 1, function(x) colnames(aicc)[which.min(x)]) # Select minimum aic from each row of matrix
aicc_table <- cbind(as.character(Models$V1), best_mod) # Create a table of each ID and best model
# Count the number of times each model was selected as the best model
one <- sum(aicc_table[,2] == "C_AICC")
two <- sum(aicc_table[,2] == "bAICC")
three <- sum(aicc_table[,2] == "Snl_AICC")
four <- sum(aicc_table[,2] == "Snh_AICC")
five <- sum(aicc_table[,2] == "Full_AICC")
#sum(c_best + b_best + snl_best + snh_best + full_best)
fig2 <- data.frame(AICC = c("Cubic", "Briere", "S1", "S2", "S3"), Count = c(one, two, three, four, five))

# Barplot of best models 
pdf(file = paste("../Results/Fig_2.pdf"))
print(ggplot(fig2, aes(x = AICC, y = Count, fill = AICC)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Count), vjust=1.6, color="white", size=3.5) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal(base_size = 16))
invisible(dev.off())

### Compare traits

aicc_df <- data.frame(aicc_table) # Best model per ID
traits <- unique(Data %>% group_by(FinalID) %>% select(StandardisedTraitName))  # For each ID select it's Trait name
trait_type <- cbind(aicc_df, traits[,2]) # Trait type per best selected model
# Rename traits as either respiration, photosynthesis, growth
#trait_type$StandardisedTraitName <- gsub("^.?photosynth", "P", trait_type$StandardisedTraitName)
#trait_type$StandardisedTraitName <- gsub("^.?Growth", "G", trait_type$StandardisedTraitName)
#trait_type$StandardisedTraitName <- str_replace_all(trait_type, 'Respiration', 'R')
# Create vectors of related trait names 
R <- c("Mass-Specific Respiration Rate", "respiration rate", "Surface Area-Specific Dark Respiration Rate",
       "Surface Area-Specific Mitochondrial Respiration Rate")
P <- c("Surface Area-Specific Maximum Photosynthesis Rate", "net photosynthesis", "gross photosynthesis",
       "Mass-Specific Photosynthetic Oxygen Production Rate", "Surface Area-Specific Photosynthetic Oxygen Production Rate",
       "POC Photosynthetic Oxygen Production Rate", "Petic 14-NaHCO3 Incorperation", "Chlorophyll-a-Specific Carbon Production Rate", "photosynthetic 14-NaHCO3 Incorperation", "Rate of photosynthesis")
G <- c("Population Growth Rate", "Individual Mass Growth Rate", "Radial Growth Rate", "Individual Length Growth Rate", "Specific Growth Rate")

# Rename trait levels into categories
levels(trait_type$StandardisedTraitName)[levels(trait_type$StandardisedTraitName) %in% R] <-"Respiration"
levels(trait_type$StandardisedTraitName)[levels(trait_type$StandardisedTraitName) %in% P] <-"Photosynthesis"
levels(trait_type$StandardisedTraitName)[levels(trait_type$StandardisedTraitName) %in% G] <-"Growth"

# Generate proportion plot
pdf(file = paste("../Results/Fig_3.pdf"))
print(ggplot(trait_type, aes(x = StandardisedTraitName, fill = best_mod)) +
  geom_bar(position = "fill", width = 4) +
  ylab("Proportion") +
  xlab("Trait Type") +
  facet_grid(.~StandardisedTraitName) +
  scale_fill_brewer("Best Model", labels = c("Briere", "Cubic", "S3", "S2", "S1"), palette = "Dark2") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())) 
invisible(dev.off())

### Compare Kingdoms
King <- unique(Data %>% group_by(FinalID)  %>% select(ConKingdom)) # For each ID select consumer Kingdom
Kingdom <- cbind(aicc_df, King[,2]) # Kingdom per best selected model
Kingdom <- Kingdom[!(Kingdom$ConKingdom == ""),] # ID for which Kingdom data is not present


# Generate proportion plot
#pdf(file = paste("../Results/Fig_4.pdf"))
ggplot(Kingdom, aes(x = ConKingdom, fill = best_mod)) +
  geom_bar(position = "fill", width = 12) +
  ylab("Proportion") +
  xlab("Kingdom") +
  facet_grid(.~ConKingdom) +
  scale_fill_brewer("Best Model", labels = c("Briere", "Cubic", "S3", "S2", "S1"), palette = "Dark2") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())
ggsave(filename = paste("../Results/Fig_4.pdf"), width = 300, height = 150, units = "mm")
invisible(dev.off())

# Calculate mean activation energies
#E_vals <- as.data.frame(Models[,grepl("E", names(Models))]) # Create a matrix of all aicc values
#Average <- apply(na.omit(E_vals), 2, function(x) mean(x))

# Appendix 2
aicc_tally <- cbind(Observations, aicc) # Number of observations per ID and selected score
aicc_tally <- na.omit(aicc_tally) # Omit NA's
aicc_tally <- aicc_tally[(aicc_tally$n < 13),] # Remove values greater than 13

# Print plot
pdf(file = paste("../Results/Appendix2.pdf"))
print(ggplot(data = aicc_tally) + 
  ylab("AICC score") +
  xlab("Sample Size") +
  geom_point(aes(x = n, y = C_AICC, colour = "Cubic")) +
  geom_point(aes(x = n, y = bAICC, color = "Briere")) +
  geom_point(aes(x = n, y = Snl_AICC, color = "S1")) +
  geom_point(aes(x = n, y = Snh_AICC, color = "S2")) +
  geom_point(aes(x = n, y = Full_AICC, color = "S3")) +
  labs(color = "Models") +
  theme_minimal(base_size = 16))
invisible(dev.off())

finish_time2 = (proc.time()[3] - initialtime)/60  
print(paste("Script took", finish_time2, "minutes to run"))