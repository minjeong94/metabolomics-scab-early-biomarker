library(glmnet)
library(tidyverse)
library(e1071)
# Set working directory
setwd("C:/Users/Min Jeong Kang/Desktop")
# Read the data
data <- read_csv("scab_data_final.csv")

# Add prefix "meta_" to metabolite columns
new_names <- paste0("meta_", names(data[, c(6:ncol(data))]))
colnames(data)[6:ncol(data)] <- new_names
glimpse(data)
# Normalize metabolite levels
data <- data %>%
  mutate_at(vars(matches("meta_")), ~scale(.x)) %>%
  mutate_at(vars(matches("meta_")), ~as.numeric(.x))

# Convert hpi to a factor, For 3-4DPI and 5-7DPI, change the filter()
data_1 <- data %>%
  mutate(hpi = factor(hpi)) %>% filter(hpi %in% c("1day","2day")) %>% 
  mutate(hpi = factor(as.character(hpi)))

#-----------LASSO REGRESSION + LOGISTIC CONTROL x INFECTED--------------------------

#--------Control VS Susceptible-------------

# Divide the data into resistant and susceptible)
data_con_sus <- data_1 %>% filter(treatment %in% c("control", "susceptible")) 

data_con_sus <- data_con_sus %>%
  mutate(infected = ifelse(treatment == "control", 0, 1))

glimpse(data_con_sus)
#---------------------------------------

lasso_model_infected <- function(data) {  
  sample_rows <- sample(1:nrow(data), size = 8, replace = FALSE)
  test_data <- data[sample_rows, ]
  train_data <- data[-sample_rows, ]
  
  # Subset the outcome column using the column name
  train_Y <- as.numeric(train_data$infected) 
  test_Y <- as.numeric(test_data$infected) 
  
  # Select predictor columns for training and testing data (C30)/ For HILIC, use c(77:140)
  train_X <- train_data[, c(6:76)]
  test_X <- test_data[, c(6:76)]   
  
  # Convert to model.matrix
  train_X <- model.matrix(~ . - 1, data = train_X)
  test_X <- model.matrix(~ . - 1, data = test_X)
  
  # Perform cross-validation to find the optimal lambda
  cvfit <- cv.glmnet(
    train_X,
    train_Y,
    family = "binomial",         # Specify logistic regression
    type.measure = "class",      # Use "class" for classification error rate
    nfolds = 10,                 # Number of cross-validation folds
    alpha = 0.999                # Alpha close to 1 for Lasso
  )  
  
  lambda_min <- cvfit$lambda.min
  
  # Fit the final model using the optimal lambda
  lasso_reg <- glmnet(
    train_X,
    train_Y,
    lambda = lambda_min,
    family = "binomial"
  )
  prediction <- predict(lasso_reg, newx = test_X)
  # Convert predictions to binary by rounding to nearest integer (you may adjust threshold if needed)
  predicted_class <- ifelse(prediction >= 0.5, 1, 0)
  
  # Calculate true positives, false positives, false negatives, true negatives
  tp <- sum(predicted_class == 1 & test_Y == 1)
  fp <- sum(predicted_class == 1 & test_Y == 0)
  fn <- sum(predicted_class == 0 & test_Y == 1)
  tn <- sum(predicted_class == 0 & test_Y == 0)
  
  # Calculate accuracy, sensitivity, specificity, and precision
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  precision <- tp / (tp + fp)
  
  f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
  
  # Calculate Matthews correlation coefficient (MCC)
  mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  
  # Extract coefficients for the optimal lambda
  coefficients <- coef(lasso_reg)
  
  coefficients_df <-tibble(
    metabolite = rownames(coefficients),
    coefficient = coefficients[,1],
  )
  
  return(list(
    tibble(
      accuracy = accuracy,
      sensitivity = sensitivity,
      specificity = specificity,
      precision = precision,
      tp = tp, fp = fp, fn = fn, tn = tn,
      f1_score = f1_score,
      mcc = mcc
    ), coefficients = coefficients_df))
}

# Replicate 20 times
lasso_con_sus <- map(1:20, ~lasso_model_infected(data = data_con_sus))

lasso_cont_sus <- bind_rows(lapply(lasso_con_sus, function(sublist) sublist[[1]]))
print(lasso_cont_sus)

lasso_sus_av <- lasso_cont_sus %>% summarise_all(mean)*100 
lasso_sus_sd <- lasso_cont_sus %>%summarise_all(sd)*100
lasso_sus_results <- bind_rows(lasso_sus_av,lasso_sus_sd) %>%
  mutate(
    tp = tp/100,
    fp = fp/100,
    fn = fn/100,
    tn = tn/100
  ) %>%
  mutate_all(~ round(., 2))
lasso_sus_results

# List of feature importance scores for each cross-validation run
lass_sus_coeff <- lapply(1:20, function(sublist) lasso_con_sus[[sublist]]$coefficients)
lasso_sus_coeffc<- bind_rows(lass_sus_coeff) %>% 
  group_by(metabolite) %>% 
  summarise_all(mean)%>%
  filter(coefficient != 0) %>%
  filter(metabolite != "(Intercept)")

print(lasso_sus_coeffc, n=31)

# Export the data into xlsx, Save the file accordingly
library(writexl)
library(openxlsx)

write_xlsx(lasso_sus_coeffc, path = ".xlsx")
write.xlsx(lasso_sus_coeffc, file = ".xlsx")


#--------Control VS Resistant-------------

# Divide the data into resistant and susceptible)
data_con_res <- data_1 %>% filter(treatment %in% c("control", "resistant")) 

data_con_res <- data_con_res %>%
  mutate(infected = ifelse(treatment == "control", 0, 1))

glimpse(data_con_res)
#---------------------------------------

lasso_model_infected <- function(data) {  
  sample_rows <- sample(1:nrow(data), size = 8, replace = FALSE)
  test_data <- data[sample_rows, ]
  train_data <- data[-sample_rows, ]
  
  # Subset the outcome column using the column name
  train_Y <- as.numeric(train_data$infected) 
  test_Y <- as.numeric(test_data$infected) 
  
  # Select predictor columns for training and testing data(C30)/ For HILIC, use c(77:140)
  train_X <- train_data[, c(6:76)]
  test_X <- test_data[, c(6:76)]
  
  # Convert to model.matrix
  train_X <- model.matrix(~ . - 1, data = train_X)
  test_X <- model.matrix(~ . - 1, data = test_X)
  
  # Perform cross-validation to find the optimal lambda
  cvfit <- cv.glmnet(
    train_X,
    train_Y,
    family = "binomial",         # Specify logistic regression
    type.measure = "class",      # Use "class" for classification error rate
    nfolds = 10,                 # Number of cross-validation folds
    alpha = 0.999                # Alpha close to 1 for Lasso
  )  
  
  lambda_min <- cvfit$lambda.min
  
  # Fit the final model using the optimal lambda
  lasso_reg <- glmnet(
    train_X,
    train_Y,
    lambda = lambda_min,
    family = "binomial"
  )
  prediction <- predict(lasso_reg, newx = test_X)
  # Convert predictions to binary by rounding to nearest integer (you may adjust threshold if needed)
  predicted_class <- ifelse(prediction >= 0.5, 1, 0)
  
  # Calculate true positives, false positives, false negatives, true negatives
  tp <- sum(predicted_class == 1 & test_Y == 1)
  fp <- sum(predicted_class == 1 & test_Y == 0)
  fn <- sum(predicted_class == 0 & test_Y == 1)
  tn <- sum(predicted_class == 0 & test_Y == 0)
  
  # Calculate accuracy, sensitivity, specificity, and precision
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  precision <- tp / (tp + fp)
  
  f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
  
  # Calculate Matthews correlation coefficient (MCC)
  mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  
  # Extract coefficients for the optimal lambda
  coefficients <- coef(lasso_reg)
  
  coefficients_df <-tibble(
    metabolite = rownames(coefficients),
    coefficient = coefficients[,1],
  )
  
  return(list(
    tibble(
      accuracy = accuracy,
      sensitivity = sensitivity,
      specificity = specificity,
      precision = precision,
      tp = tp, fp = fp, fn = fn, tn = tn,
      f1_score = f1_score,
      mcc = mcc
    ), coefficients = coefficients_df))
}

# Replicate 20 times
lasso_con_res <- map(1:20, ~lasso_model_infected(data = data_con_res))

lasso_cont_res <- bind_rows(lapply(lasso_con_res, function(sublist) sublist[[1]]))
print(lasso_cont_res)

lasso_res_av <- lasso_cont_res %>% summarise_all(mean)*100 
lasso_res_sd <- lasso_cont_res %>%summarise_all(sd)*100
lasso_res_results <- bind_rows(lasso_res_av,lasso_res_sd) %>%
  mutate(
    tp = tp/100,
    fp = fp/100,
    fn = fn/100,
    tn = tn/100
  ) %>%
  mutate_all(~ round(., 2))
lasso_res_results

# List of feature importance scores for each cross-validation run
lass_res_coeff <- lapply(1:20, function(sublist) lasso_con_res[[sublist]]$coefficients)
lasso_res_coeffc<- bind_rows(lass_res_coeff) %>% 
  group_by(metabolite) %>% 
  summarise_all(mean)%>%
  filter(coefficient != 0) %>%
  filter(metabolite != "(Intercept)")

print(lasso_res_coeffc, n=31)

LR_performance_1 <- bind_rows(
  lasso_sus_results %>% mutate(source = "cont_vs_sus"),
  lasso_res_results %>% mutate(source = "cont_vs_res")
)


LR_coeff_res_sus_withControl <- bind_rows(
  lasso_sus_coeffc %>% mutate(source = "cont_vs_sus"),
  lasso_res_coeffc %>% mutate(source = "cont_vs_res")
)

# Export the data into xlsx, Save the file accordingly
library(writexl)
library(openxlsx)

write_xlsx(LR_coeff_res_sus_withControl, path = ".xlsx")
write.xlsx(LR_coeff_res_sus_withControl, file = ".xlsx")

write_xlsx(LR_performance_1, path = ".xlsx")
write.xlsx(LR_performance_1, file = ".xlsx")


