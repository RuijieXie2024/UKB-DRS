# ===============================================================
# Title:     Novel type 2 diabetes prediction score based on traditional risk factors and circulating metabolites: Model derivation and validation in two large cohort studies
# Author:    [Ruijie Xie]
# Date:      [30/09/2024]
# ===============================================================

# ---------------------------
# Load Required Packages
# ---------------------------

library(survival)
library(glmnet)
library(survIDINRI)
library(nricens)
library(pec)
library(dplyr)
library(tidyr)
library(doParallel)
library(missForest)
library(ggplot2)
library(runway)

# ---------------------------
# Data Preparation
# ---------------------------

# Note: Replace 'UKB_data.csv' and 'ESTHER_data.csv'
# with the actual file names when using real data.

# Load derivation cohort
UKB_data <- read.csv('UKB_data.csv')

# Load validation cohort
ESTHER_data <- read.csv('ESTHER_data.csv')


# ---------------------------
# Variable Imputation
# ---------------------------

# Impute missing values in metabolites using missForest
missForest <- function(xmis, maxiter = 10, ntree = 100, variablewise = FALSE,
                       decreasing = FALSE, verbose = FALSE,
                       mtry = floor(sqrt(ncol(xmis))), replace = TRUE,
                       classwt = NULL, cutoff = NULL, strata = NULL,
                       sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                       xtrue = NA, parallelize = c('no', 'variables', 'forests'))
  

# ---------------------------
# Variable Preparation
# ---------------------------

# List of variables with missing values (e.g., metabolite variables)
metabolite_vars <- c('Metabolite1', 'Metabolite2', 'Metabolite3', 'Metabolite4',
                     'Metabolite5', 'Metabolite6', 'Metabolite7', 'Metabolite8',
                     'Metabolite9', 'Metabolite10', '.........' , 'Metabolite249')


# Log-transform and standardize metabolite variables
for (var in metabolite_vars) {
  training_data_imputed[[var]] <- scale(log(training_data_imputed[[var]] + 1))
  validation_data_imputed[[var]] <- scale(log(validation_data_imputed[[var]] + 1))
  external_data_imputed[[var]] <- scale(log(external_data_imputed[[var]] + 1))
}




# ---------------------------
# co-correlation matrix of selected metabolites
# ---------------------------
# Define and Order Metabolite List with Standardized Names
# Create a named vector where the keys are original names and the values are standardized names
metabolite_mapping <- c(
  "bOHbutyrate" = "3-Hydroxybutyrate",
  "Acetate" = "Acetate",
  "Citrate" = "Citrate",
  "Gln" = "Glutamine",
  "Glucose" = "Glucose",
  "IDL_CE_pct" = "IDL-CE-pct",
  "LA_pct" = "LA-pct",
  "Lactate" = "Lactate",
  "M_LDL_TG_pct" = "M-LDL-TG-pct",
  "Pyruvate" = "Pyruvate",
  "Tyr" = "Tyrosine"
)

# Specify the order for sorting
metabolites_ordered <- c(
  "3-Hydroxybutyrate",
  "Acetate",
  "Citrate",
  "Glutamine",
  "Glucose",
  "IDL-CE-pct",
  "LA-pct",
  "Lactate",
  "M-LDL-TG-pct",
  "Pyruvate",
  "Tyrosine"
)

# Rename Metabolite Columns in the Datasets
# Assume that the datasets 'training_data' and 'validation_data' have been loaded into the R environment

# Define a function to rename columns in a dataframe based on a mapping
rename_metabolites <- function(df, mapping) {
  # Check if all original names exist in the dataframe
  if(!all(names(mapping) %in% colnames(df))){
    missing <- names(mapping)[!(names(mapping) %in% colnames(df))]
    stop(paste("The following metabolites are missing in the dataset:", paste(missing, collapse = ", ")))
  }
  
  # Rename the columns
  df_renamed <- df
  colnames(df_renamed)[match(names(mapping), colnames(df_renamed))] <- mapping
  return(df_renamed)
}

# Apply the renaming function to both datasets
training_data_renamed <- rename_metabolites(training_data, metabolite_mapping)
validation_data_renamed <- rename_metabolites(validation_data, metabolite_mapping)

# Extract Relevant Data
# Select metabolites in the specified order
training_metabolites <- training_data_renamed[, metabolites_ordered]
validation_metabolites <- validation_data_renamed[, metabolites_ordered]

# Calculate Correlation Matrices
# Use Pearson correlation coefficient
cor_training <- cor(training_metabolites, use = "pairwise.complete.obs", method = "pearson")
cor_validation <- cor(validation_metabolites, use = "pairwise.complete.obs", method = "pearson")

# Generate Aesthetic Correlation Matrix Plots
# Define a unified theme
theme_set(theme_minimal())

# Create correlation matrix plot for the training dataset
p1 <- ggcorrplot(cor_training, 
                 hc.order = FALSE,  # Do not reorder based on hierarchical clustering
                 type = "upper",
                 lab = TRUE,
                 lab_size = 3,
                 method = "square",
                 colors = c("blue", "white", "red"),
                 title = "Training Dataset",
                 ggtheme = theme_minimal()) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Create correlation matrix plot for the validation dataset
p2 <- ggcorrplot(cor_validation, 
                 hc.order = FALSE,  # Do not reorder based on hierarchical clustering
                 type = "upper",
                 lab = TRUE,
                 lab_size = 3,
                 method = "square",
                 colors = c("blue", "white", "red"),
                 title = "Validation Dataset",
                 ggtheme = theme_minimal()) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Use gridExtra to display both plots side by side
combined_plot <- grid.arrange(p1, p2, ncol = 2)




# ---------------------------
# Variable Selection with Cox-LASSO
# ---------------------------

# Combine clinical variables and metabolite variables
predictor_vars <- c('Age', 'Sex', 'BMI', 'FamilyHistory', 'SmokingStatus', 'BloodPressure', 'SteroidUse', 'HbA1c', metabolite_vars)

# Prepare data for variable selection
X <- as.matrix(training_data_imputed[, predictor_vars])
y <- Surv(training_data_imputed$FollowUpTime, training_data_imputed$Event == 1)

# Set up parallel processing
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Number of bootstrap samples
n_bootstrap <- 1000

# Store selected variables from each bootstrap sample
selected_variables_list <- foreach(i = 1:n_bootstrap, .packages = c('glmnet', 'survival')) %dopar% {
  # Bootstrap sample
  bootstrap_indices <- sample(1:nrow(training_data_imputed), replace = TRUE)
  bootstrap_sample <- training_data_imputed[bootstrap_indices, ]
  
  # Prepare data
  X_boot <- as.matrix(bootstrap_sample[, predictor_vars])
  y_boot <- Surv(bootstrap_sample$FollowUpTime, bootstrap_sample$Event == 1)
  
  # Perform LASSO with cross-validation
  cv.lasso_boot <- cv.glmnet(X_boot, y_boot, family = 'cox', alpha = 1, nfolds = 10)
  best_lambda_boot <- cv.lasso_boot$lambda.min
  lasso_model_boot <- glmnet(X_boot, y_boot, family = 'cox', alpha = 1, lambda = best_lambda_boot)
  selected_vars_boot <- rownames(coef(lasso_model_boot))[which(coef(lasso_model_boot) != 0)]
  selected_vars_boot[selected_vars_boot != '(Intercept)']
}

# Stop parallel processing
stopCluster(cl)

# Aggregate selection frequency
variable_selection_counts <- table(unlist(selected_variables_list))
variable_selection_percentages <- (variable_selection_counts / n_bootstrap) * 100

# Select variables appearing in at least 95% of bootstrap samples
stable_variables <- names(variable_selection_percentages[variable_selection_percentages >= 95])
print('Variables selected in at least 95% of bootstrap samples:')
print(stable_variables)


# ---------------------------
# Model Development
# ---------------------------

# Fit the old model (clinical variables only)
old_model_formula <- as.formula('Surv(FollowUpTime, Event == 1) ~ Age + Sex + BMI + FamilyHistory + SmokingStatus + BloodPressure + SteroidUse + HbA1c')
old_model <- coxph(old_model_formula, data = training_data_imputed)

# Fit the new model (clinical variables + selected metabolites)
new_model_formula <- as.formula(paste('Surv(FollowUpTime, Event == 1) ~', paste(stable_variables, collapse = ' + ')))
new_model <- coxph(new_model_formula, data = training_data_imputed)


# ---------------------------
# Model Validation
# ---------------------------

# Function to calculate Harrell's C-index
calculate_cindex <- function(model, data) {
  surv_obj <- Surv(data$FollowUpTime, data$Event == 1)
  cindex <- concordance(model, newdata = data)$concordance
  return(cindex)
}

# Calculate C-index for validation data
cindex_old_validation <- calculate_cindex(old_model, validation_data_imputed)
cindex_new_validation <- calculate_cindex(new_model, validation_data_imputed)

# Calculate C-index for external data
cindex_old_external <- calculate_cindex(old_model, external_data_imputed)
cindex_new_external <- calculate_cindex(new_model, external_data_imputed)

# Print C-index results
cat('C-index for Old Model (Validation Data):', round(cindex_old_validation, 3), '\n')
cat('C-index for New Model (Validation Data):', round(cindex_new_validation, 3), '\n')
cat('C-index for Old Model (External Data):', round(cindex_old_external, 3), '\n')
cat('C-index for New Model (External Data):', round(cindex_new_external, 3), '\n')


# ---------------------------
# Compare C-index in diffrent models
# ---------------------------
####################CompareC##########################
result <- compareC(timeX = your_data$dur_diab, statusX = your_data$i_diab, scoreY = your_data$risk_old, scoreZ = your_data$risk_new)

# print result
print(result)

est_diff_c <- 0.0126978
est_vardiff_c <- 2.541948e-06




# ---------------------------
# NRI and IDI Calculations
# ---------------------------

#--- sample data ---
D=subset(pbc, select=c('Age', 'Sex', 'BMI', 'FamilyHistory', 'SmokingStatus', 'BloodPressure', 'SteroidUse', 'HbA1c', metabolite_vars))
D$status=as.numeric(D$status==2)
D=D[!is.na(apply(D,1,mean)),] ; dim(D)
mydata=D[1:100,]

t0=365.25*10
indata1=mydata;
indata0=mydata[,-7] ; n=nrow(D) ;
covs1<-as.matrix(indata1[,c(-1,-2)])
covs0<-as.matrix(indata0[,c(-1,-2)])

#--- inference ---
x<-IDI.INF(mydata[,1:2], covs0, covs1, t0, npert=500) ;

#--- results ---
IDI.INF.OUT(x) ;

#--- Graphical presentaion of the estimates ---
# IDI.INF.GRAPH(x) ;


indat <- your_data[,c("dur_diab","i_diab")]

# Risk of model variable prediction for old and new models
x_old <- as.matrix(your_data[, "risk_old"])
x_new <- as.matrix(your_data[, "risk_new"])


# Perform 500 bootstrap calculations of IDI and NRI in 3652.5 days
result <- IDI.INF(indat, covs0 = x_old, covs1 = x_new, t0 = 3652.5, npert = 500)

# output result
print(IDI.INF.OUT(result))
IDI.INF.GRAPH(result)



# ---------------------------
# Forest plot
# ---------------------------

# Extract the biomarker names
selected_bmrs <- df %>% pull(name)

# Join the association data frame with group data
df <-
  df_logodds_associations %>%
  # Set the study variable to a factor to preserve order of appearance
  mutate(
    study = factor(
      study,
      levels = c("training set", "test set", "external validation set")
    )
  ) %>%
  # use right_join, with df_grouping on the right, to preserve the order of
  # biomarkers it specifies.
  dplyr::right_join(., df_grouping, by = "name") %>%
  tidyr::drop_na(.data$beta)

# Draw a forestplot of odds ratios
ggforestplot::forestplot(
  df = df,
  name = name,
  estimate = beta,
  se = se,
  pvalue = pvalue,
  psignif = 0.002,
  colour = study,
  xlab = "Odds ratio for incident type 2 diabetes (95% CI)\nper 1−SD increment in metabolite concentration",
  title = "Metabolites",
  logodds = TRUE
)



# ---------------------------
# ROC curves
# ---------------------------

# Temporary Dataset 1 - clinical CDRS
temp_dataset1 <- UKB_data[, c("i_diab", "risk1")]
temp_dataset1$outcomes <- temp_dataset1$i_diab
temp_dataset1$predictions <- temp_dataset1$risk1
temp_dataset1$model_name <- 'Base model'
temp_dataset1 <- temp_dataset1[, c("outcomes", "predictions", "model_name")]

# Temporary Dataset 2 - Concise UKB-DRS
temp_dataset2 <- UKB_data[, c("i_diab", "risk2")]
temp_dataset2$outcomes <- temp_dataset2$i_diab
temp_dataset2$predictions <- temp_dataset2$risk2
temp_dataset2$model_name <- 'Base model + Metabolites'
temp_dataset2 <- temp_dataset2[, c("outcomes", "predictions", "model_name")]

# Temporary Dataset 3 - Concise UKB-DRS
temp_dataset3 <- UKB_data[, c("i_diab", "risk3")]
temp_dataset3$outcomes <- temp_dataset3$i_diab
temp_dataset3$predictions <- temp_dataset3$risk3
temp_dataset3$model_name <- 'Base model + Metabolites'
temp_dataset3 <- temp_dataset3[, c("outcomes", "predictions", "model_name")]

# Combine the three datasets
runway_UKB <- rbind(temp_dataset1, temp_dataset2, temp_dataset3)


roc_plot_multi(runway_UKB, 
               outcome = 'outcomes', 
               positive = '1',
               prediction = 'predictions', 
               model = 'model_name',
               ci = TRUE,
               plot_title = 'ROC curves')



# ---------------------------
# Calibration curves
# ---------------------------

# Temporary Dataset 1 - clinical CDRS
temp_dataset1 <- UKB_data[, c("i_diab", "risk1")]
temp_dataset1$outcomes <- temp_dataset1$i_diab
temp_dataset1$predictions <- temp_dataset1$risk1
temp_dataset1$model_name <- 'Base model'
temp_dataset1 <- temp_dataset1[, c("outcomes", "predictions", "model_name")]

# Temporary Dataset 2 - Concise UKB-DRS
temp_dataset2 <- UKB_data[, c("i_diab", "risk2")]
temp_dataset2$outcomes <- temp_dataset2$i_diab
temp_dataset2$predictions <- temp_dataset2$risk2
temp_dataset2$model_name <- 'Base model + Metabolites'
temp_dataset2 <- temp_dataset2[, c("outcomes", "predictions", "model_name")]

# Temporary Dataset 3 - Concise UKB-DRS
temp_dataset3 <- UKB_data[, c("i_diab", "risk3")]
temp_dataset3$outcomes <- temp_dataset3$i_diab
temp_dataset3$predictions <- temp_dataset3$risk3
temp_dataset3$model_name <- 'Base model + Metabolites'
temp_dataset3 <- temp_dataset3[, c("outcomes", "predictions", "model_name")]

# Combine the three datasets
runway_UKB <- rbind(temp_dataset1, temp_dataset2, temp_dataset3)


cal_plot_multi(runway_UKB,
               outcome = 'outcomes',
               positive = '1',
               prediction = 'predictions',
               model = 'model_name',
               n_bins = 10,
               plot_title = 'Calibration plot')
