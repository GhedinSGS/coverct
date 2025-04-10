library(pacman)
p_load(tidyverse, duckdb, dbplyr, ranger)

cs = function(fraction_below) {
  return(log2(fraction_below + 1e-4))
}

coverage_featurize = function(cov_ct_nods) {
  cov_ct_nods %>%
    group_by(PatientID, replicate, CT, downsample) %>%
    summarize(
      mean_dp = mean(COV_DP),
      sd_dp = sd(COV_DP),#,
      # Fraction below values
    fraction1x    = mean(if_else(COV_DP < 1, 1L, 0L), na.rm = T),
    fraction10x   = mean(if_else(COV_DP < 10, 1L, 0L), na.rm = T),
    fraction20x   = mean(if_else(COV_DP < 20, 1L, 0L), na.rm = T),
    fraction50x   = mean(if_else(COV_DP < 50, 1L, 0L), na.rm = T),
    fraction100x  = mean(if_else(COV_DP < 100, 1L, 0L), na.rm = T),
    fraction200x  = mean(if_else(COV_DP < 200, 1L, 0L), na.rm = T),
    fraction500x  = mean(if_else(COV_DP < 500, 1L, 0L), na.rm = T),
    fraction1000x = mean(if_else(COV_DP < 1000, 1L, 0L), na.rm = T),
    fraction1500x = mean(if_else(COV_DP < 1500, 1L, 0L), na.rm = T),
    fraction2000x = mean(if_else(COV_DP < 2000, 1L, 0L), na.rm = T),
    fraction3000x = mean(if_else(COV_DP < 3000, 1L, 0L), na.rm = T),
    fraction4000x = mean(if_else(COV_DP < 4000, 1L, 0L), na.rm = T),
    fraction5000x = mean(if_else(COV_DP < 5000, 1L, 0L), na.rm = T),
    fraction10000x = mean(if_else(COV_DP < 10000, 1L, 0L), na.rm = T)
    ) %>%
    ungroup %>% distinct() %>%
    collect() %>%
    mutate(
    CHAOS = -mean_dp/sd_dp,
    #   # Calculate cs values using the cs function
      cs1    = cs(fraction1x),
      cs10   = cs(fraction10x),
      cs20   = cs(fraction20x),
      cs50   = cs(fraction50x),
      cs100  = cs(fraction100x),
      cs200  = cs(fraction200x),
      cs500  = cs(fraction500x),
      cs1000 = cs(fraction1000x),
      cs1500 = cs(fraction1500x),
      cs2000 = cs(fraction2000x),
      cs3000 = cs(fraction3000x),
      cs4000 = cs(fraction4000x),
      cs5000 = cs(fraction5000x),
      cs10000 = cs(fraction10000x)
    )
}

wrap_rf = function(hmh_train, features, size = 0.8, 
                   test_data = NULL, num_trees = 500,
                   f1_threshold = 22, model_type = "rf") {
  data_split <- split_data(hmh_train, size = size)
  rf_model <- fit_model(data_split$train_data, model_type, features = features, num_trees = num_trees)
  # test data
  if (is.null(test_data)) test_data = data_split$test_data
  list_out = evaluate_ranger(rf_model, test_data, f1_threshold = f1_threshold, model_type)
  return(list_out)
}

split_data <- function(data, run_col = "run", ct_col = "CT", size = 0.8, 
                       seed = 42, n_bins = 5) {
  set.seed(seed)
  
  # Check if the necessary columns exist
  if (run_col %in% names(data) && ct_col %in% names(data)) {
    # Add CT bins within each run
    data <- data %>%
      group_by(!!sym(run_col)) %>%
      dplyr::mutate(CT_bin = cut(!!sym(ct_col), breaks = n_bins, 
                                 include.lowest = TRUE, labels = FALSE)) %>%
      ungroup()
    
    # Perform stratified sampling within each run and CT bin
    train_data <- data %>%
      group_by(!!sym(run_col), CT_bin) %>%
      sample_frac(size, replace = FALSE) %>%
      ungroup()  # Always ungroup after operations to avoid unexpected behavior
    
    # Determine the unique identifiers of train data assuming 'mcov_id' is the unique ID
    train_ids <- train_data %>% pull(mcov_id)
    
    # Create test data by excluding train_ids
    test_data <- data %>% 
      filter(!mcov_id %in% train_ids) %>%
      ungroup()  # Ensure data is ungrouped for further operations
  } else {
    # If one of the columns doesn't exist, revert to random sampling
    train_data <- data %>% sample_frac(size, replace = FALSE)
    test_data <- data %>% anti_join(train_data, by = "mcov_id")
  }
  
  return(list(train_data = train_data, test_data = test_data))
}


fit_model <- function(train_data, model_type = "rf", num_trees = 500,
                      features = "DADI") {
  set.seed(42)
  train_data <- train_data %>% select(one_of(c("CT", features)))
  if (model_type == "lm") {
    return(lm(CT ~ ., 
              data = train_data))
  } else {
    return(ranger(CT ~ ., data = train_data, 
                  num.trees = num_trees, importance = 'impurity'))
  }
}

evaluate_ranger = function(rf_model, test_data, f1_threshold = 22, model_type = "rf") {
  rf_predictions <- make_predictions(rf_model, test_data, model_type)
  rf_results <- calculate_r_squared(test_data$CT, rf_predictions, 
                                    length(rf_model$variable.importance))
  rf_mape = calculate_mae_mape(test_data$CT, rf_predictions)
  rf_f1 = calculate_classification_metrics(test_data$CT, 
                                           rf_predictions, threshold = f1_threshold)$f1_score
  
  # Calculate additional performance metrics
  pearson_corr <- cor(test_data$CT, rf_predictions, method = "pearson")
  spearman_corr <- cor(test_data$CT, rf_predictions, method = "spearman")
  
  if (model_type == "rf") {
    feature_importances <- data.frame(Variable = names(rf_model$variable.importance), 
                                      Importance = rf_model$variable.importance) %>% 
      arrange(-Importance) %>% mutate(Variable = factor(Variable, levels = Variable))
    out_list = list(model = rf_model, pred = rf_predictions, actual = test_data, 
                    performance = list(f1 = rf_f1, r2 = rf_results, mape = rf_mape,
                                       pearson = pearson_corr, spearman = spearman_corr), 
                    importance = feature_importances, f1_threshold = f1_threshold)
  } else {
    out_list = list(model = rf_model, pred = rf_predictions, actual = test_data, 
                    performance = list(f1 = rf_f1, r2 = rf_results, mape = rf_mape,
                                       pearson = pearson_corr, spearman = spearman_corr),
                    f1_threshold = f1_threshold)
  }
  
  return(out_list)
}


# Function to calculate R-squared and Adjusted R-squared
calculate_r_squared <- function(actuals, predictions, num_predictors) {
  sst <- sum((actuals - mean(actuals))^2)
  ssr <- sum((actuals - predictions)^2)
  r_squared <- 1 - ssr / sst
  n <- length(actuals)
  adjusted_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - num_predictors - 1))
  return(list(r_squared = r_squared, adjusted_r_squared = adjusted_r_squared))
}

calculate_mae_mape <- function(actuals, predictions) {
  # Calculate Mean Absolute Error (MAE)
  mae <- mean(abs(actuals - predictions))
  
  # Calculate Mean Absolute Percentage Error (MAPE)
  mape <- mean(abs((actuals - predictions) / actuals))
  
  # Return both MAE and MAPE
  return(list(MAE = mae, MAPE = mape))
}

calculate_classification_metrics <- function(actuals, predictions, threshold) {
  # Convert continuous predictions into binary class based on the threshold
  predicted_classes <- ifelse(predictions < threshold, 1L, 0L)
  actual_classes <- ifelse(actuals < threshold, 1L, 0L)
  
  # Calculate True Positives (TP), False Positives (FP), True Negatives (TN), and False Negatives (FN)
  TP <- sum(predicted_classes == 1 & actual_classes == 1)
  FP <- sum(predicted_classes == 1 & actual_classes == 0)
  TN <- sum(predicted_classes == 0 & actual_classes == 0)
  FN <- sum(predicted_classes == 0 & actual_classes == 1)
  
  # Calculate Precision and Recall
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  # Avoid division by zero
  if (is.nan(precision) | is.nan(recall)) {
    precision <- 0
    recall <- 0
  }
  # Calculate F1 Score
  f1_score <- 
    ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / 
             (precision + recall))
  # Return a list with all metrics
  return(list(precision = precision, recall = recall, f1_score = f1_score))
}


make_predictions <- function(model, test_data, model_type = "lm") {
  if (model_type == "lm") {
    return(predict(model, newdata = test_data))
  } else {
    return(predict(model, data = test_data)$predictions)
  }
}

summarize_replicate_samples = function(summarized_patient_true_replicate_raw) {
  tally_variable_maf_counts_1 <- summarized_patient_true_replicate_raw %>%
    group_by(PatientID) %>%
    summarize(
      count_above_0.001 = sum(if_else(ALT_FREQ_1 > 0.001, 1L, 0L), na.rm = TRUE),
      count_above_0.01  = sum(if_else(ALT_FREQ_1 > 0.01,  1L, 0L), na.rm = TRUE),
      count_above_0.03  = sum(if_else(ALT_FREQ_1 > 0.03,  1L, 0L), na.rm = TRUE),
      ALT_FREQ = mean(ALT_FREQ_1)
    ) %>% mutate(run = 1) %>% as.data.frame()
  
  tally_variable_maf_counts_2 = summarized_patient_true_replicate_raw %>%
    group_by(PatientID) %>%
    summarize(
      count_above_0.001 = sum(if_else(ALT_FREQ_2 > 0.001, 1L, 0L), na.rm = TRUE),
      count_above_0.01  = sum(if_else(ALT_FREQ_2 > 0.01,  1L, 0L), na.rm = TRUE),
      count_above_0.03  = sum(if_else(ALT_FREQ_2 > 0.03,  1L, 0L), na.rm = TRUE),
      ALT_FREQ = mean(ALT_FREQ_2)
    ) %>% mutate(run = 2) %>% as.data.frame()
  
  return(list(rep1 = tally_variable_maf_counts_1, rep2 = tally_variable_maf_counts_2))
}

