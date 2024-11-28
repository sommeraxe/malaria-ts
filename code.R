# Load necessary libraries
library(stats)
library(forecast)
library(Metrics)
library(tseries)
library(imputeTS)
library(tidyr)
library(dplyr)
library(MASS)  # For robust regression
options(max.print = 2000)
options(warn = -1)
# Read the data
data <- read.csv("ncsu_colab/statewise_monthly.csv")
names = names(data)
preprocess_data <- function(data, state) {
  # Convert Date column to Date type and sort data by date
  data$Date <- as.Date(data$Date, format="%Y-%m-%d")
  data <- data[order(data$Date),]
  
  # Extract the month as a number and then convert it to a month name
  data$month_num <- format(data$Date, "%m")
  data$month_name <- months(data$Date)
  
  # Create dummy variables for each month
  month_dummies <- model.matrix(~factor(data$month_name) - 1, data)
  colnames(month_dummies) <- levels(factor(data$month_name))
  data <- cbind(data, month_dummies)
  data <- data[, !(colnames(data) %in% "January")]
  
  # Add a time column as a numeric representation of the Date
  data$time <- as.numeric(data$Date - min(data$Date))
  data$cases <- data[[state]]
  data <- data %>% drop_na(cases)
  data$cases <- 2 * sqrt(data$cases + 3/8)
  data$cases_lag12 <- dplyr::lag(data$cases, 12)
  data$diff <- data$cases - data$cases_lag12
  return(data)
}

rolling_forecast <- function(data, method, predict_ahead = 3) {
  results <- data.frame(MAPE1=numeric(), MAPE2=numeric(), MAPE3=numeric(), MAPE4=numeric(), MAPE5=numeric())
  residuals <- list()
  predictions <- list()
  end <- length(data$cases)
  limit <- end-predict_ahead-35
  v <- 1:limit
  for (i in 1:limit) {
    train <- data[1:(i+35), ]
    test <- data[(i+36):(i+35+ predict_ahead), ]
    
    if (method == "lm") {
      model1 <- lm(cases ~ poly(time, 2) + factor(month_name) - 1, data = train)
      model2 <- lm(log(cases) ~ poly(time, 2) + factor(month_name) - 1, data = train)
    } else if (method == "rlm") {
      model1 <- rlm(cases ~ poly(time, 2) + factor(month_name) - 1, data = train)
      model2 <- rlm(log(cases) ~ poly(time, 2) + factor(month_name) - 1, data = train)
    }
    
    predictions1 <- predict(model1, newdata = test)
    log_predictions <- predict(model2, newdata = test)
    v_squared <- summary(model2)$sigma^2
    v[i] <- v_squared
    predictions2 <- exp(log_predictions + (1/2 * v_squared))
    
    train <- na.omit(train)
    if (method == "lm") {
      model3 <- lm(diff ~ poly(time, 2), data = train)
    } else if (method == "rlm") {
      model3 <- rlm(diff ~ poly(time, 2), data = train)
    }
    predictions3 <- predict(model3, newdata = test)
    predictions3 <- predictions3 + test$cases_lag12
    
    if (method == "lm") {
      model4 <- lm(cases ~ cases_lag12 + poly(time, 2), data = train)
    } else if (method == "rlm") {
      model4 <- rlm(cases ~ cases_lag12 + poly(time, 2), data = train)
    }
    predictions4 <- predict(model4, newdata = test)
    
    
    models <- list(model1, model2, model3, model4)
    residuals[[i]] <- list(resid(model1), resid(model2), resid(model3), resid(model4))
    predictions[[i]] <- list(predictions1, log_predictions, predictions3, predictions4)
    mape1 <- mean(abs((test$cases - predictions1) / test$cases))
    mape2 <- mean(abs((test$cases - predictions2) / test$cases))
    mape3 <- mean(abs((test$cases - predictions3) / test$cases))
    mape4 <- mean(abs((test$cases - predictions4) / test$cases))
    
    results[i, 1:4] <- c(mape1, mape2, mape3, mape4)
  }
  # Calculate the mean of the first four columns
  column_means <- sapply(results[, 1:4], mean)
 
  # Find the index of the column with the lowest average value
  min_index <- which.min(column_means)
  
  for (i in 1:limit) {
    test <- data[(i+36):(i+35+ predict_ahead), ]
    model5 <- auto.arima(residuals[[i]][[min_index]])
    predictions5 <- forecast(model5, predict_ahead)$mean+predictions[[i]][[min_index]]
    if (min_index == 2){
      predictions5 <- exp(predictions5 + (1/2 * v[i]))
    }
    mape5 <- mean(abs((test$cases - predictions5) / test$cases))
    results[i, 5] <- mape5
  }
  
  return (list(min_index, results))
}

insample_fit <- function(data, min_index, method, state){
  end <- length(data$cases)
  # Fit models
  if (method == "lm") {
    model_with_dummies1 <- lm(cases ~ poly(time, 2) + factor(month_name) - 1, data = data)
    model_with_dummies_log <- lm(log(cases) ~ poly(time, 2) + factor(month_name) - 1, data = data)
  } else if (method == "rlm") {
    model_with_dummies1 <- rlm(cases ~ poly(time, 2) + factor(month_name) - 1, data = data)
    model_with_dummies_log <- rlm(log(cases) ~ poly(time, 2) + factor(month_name) - 1, data = data)
  }
  
  summary(model_with_dummies1)
  summary(model_with_dummies_log)
  predict_ployseas <- predict(model_with_dummies1)
  
  # Residuals for Model 1
  residuals1 <- resid(model_with_dummies1)
  plot(residuals1, main=paste("Residuals with Seasonality (Month) + Quadratic Trend for", state), ylab="Residuals", lag.max=35)
  acf(residuals1, main=paste("ACF of Residuals for", state, "- Quadratic Trend"))
  abline(h=0)
  
  # Log-transformed predictions
  v_squared <- summary(model_with_dummies_log)$sigma^2
  predicted_log_cases <- predict(model_with_dummies_log)
  predicted_log_cases <- exp(predicted_log_cases + (1/2 * v_squared))
  residuals_log <- resid(model_with_dummies_log)
  plot(residuals_log, main=paste("Residuals from Log Model with Seasonality (Month) + Quadratic for", state), ylab="Residuals")
  acf(residuals_log, main=paste("ACF of Residuals for", state, "- Log Model with Seasonality (Month) + Quadratic"), lag.max = 35)
  abline(h = 0)
  
  # Model 3: Lag with trend
  
  if (method == "lm") {
    model_with_lag_and_trend <- lm(data$diff ~ poly(time, 2), data = data)
  } else if (method == "rlm") {
    model_with_lag_and_trend <- rlm(data$diff ~ poly(time, 2), data = data)
  }
  
  summary(model_with_lag_and_trend)
  predict_lag_trend <- predict(model_with_lag_and_trend)
  predict_lag_trend <- predict_lag_trend + data$cases_lag12[13:end]
  residuals_lag <- data$cases[13:end] - predict_lag_trend
  #residuals_lag <- na.omit(residuals_lag)
  plot(residuals_lag, main=paste("Residuals from Model with Lag 12 for", state), ylab="Residuals")
  acf(residuals_lag, main=paste("ACF of Residuals for", state, "- Model with Lag 12"), lag.max = 35)
  abline(h = 0)
  
  # Model 4: Lag with polynomial trend
  if (method == "lm") {
    model_with_lag_and_polytrend <- lm(cases ~ cases_lag12 + poly(time, 2), data = data)
  } else if (method == "rlm") {
    model_with_lag_and_polytrend <- rlm(cases ~ cases_lag12 + poly(time, 2), data = data)
  }
  
  summary(model_with_lag_and_polytrend)
  predict_lag_polytrend <- predict(model_with_lag_and_polytrend)
  rmse(data$cases[13:end], predict_lag_polytrend)
  plot(predict_lag_polytrend, main=paste("Residuals from Model with Lag 12 for", state), ylab="Residuals")
  
  # Handle missing values before ACF and PACF
  residuals_lag_polytrend <- data$cases[13:end] - predict_lag_polytrend
  
  acf(residuals_lag_polytrend, main=paste("ACF of Residuals for", state, "- Model with Lag 12"), lag.max = 35)
  abline(h = 0)
  
  allpredictions = list(predict_ployseas, predicted_log_cases,predict_lag_trend, predict_lag_polytrend)
  allresiduals = list(residuals1, residuals_log, residuals_lag, residuals_lag_polytrend)
  # ARIMA Model
  
  arima_model <- auto.arima((allresiduals[[min_index]]))
  if (min_index == 1) {
    arima_fitted <- fitted(arima_model)+allpredictions[[min_index]]
    arima_residuals <- arima_fitted-data$cases
  } else if (min_index == 2) {
    arima_fitted <- exp(predict(model_with_dummies_log)+fitted(arima_model) + (1/2 * v_squared))
    arima_residuals <- arima_fitted-data$cases
  } else {
    arima_fitted <- fitted(arima_model)+allpredictions[[min_index]]
    arima_residuals <- arima_fitted-data$cases[13:end]
    arima_fitted <- c(rep(NA, 12), arima_fitted)
  }
  
  # Function to calculate metrics
  calculate_metrics <- function(actual, predicted, model = NULL) {
    mse_val <- mse(actual, predicted)
    rmse_val <- rmse(actual, predicted)
    mape_val <- mape(actual, predicted)
    rse_val <- if (!is.null(model)) rmse_val * sqrt((length(actual)) / summary(model)$df[2]) else NA
    
    return(list(MSE = mse_val, RMSE = rmse_val, RSE = rse_val, MAPE = mape_val))
  }
  
  results1 <- calculate_metrics(data$cases, predict_ployseas, model_with_dummies1)
  results2 <- calculate_metrics(data$cases, predicted_log_cases, model_with_dummies_log)
  results3 <- calculate_metrics(data$cases[13:end], predict_lag_trend, model_with_lag_and_trend)
  results4 <- calculate_metrics(data$cases[13:end], predict_lag_polytrend, model_with_lag_and_polytrend)
 
  if (min_index == 3 | min_index == 4) {
    results5 <- calculate_metrics(data$cases[13:end], arima_fitted[13:end], arima_model)
  } else {
    results5 <- calculate_metrics(data$cases, arima_fitted, arima_model)
  }
  # Ensure all data frames have the same structure before combining them
  df1 <- as.data.frame(t(unlist(results1)), stringsAsFactors = FALSE)
  df2 <- as.data.frame(t(unlist(results2)), stringsAsFactors = FALSE)
  df3 <- as.data.frame(t(unlist(results3)), stringsAsFactors = FALSE)
  df4 <- as.data.frame(t(unlist(results4)), stringsAsFactors = FALSE)
  df5 <- as.data.frame(t(unlist(results5)), stringsAsFactors = FALSE)
  
  # Ensure all data frames have the same columns
  columns <- union(union(union(colnames(df1), colnames(df2)), colnames(df3)), union(colnames(df4), colnames(df5)))
  
  # Function to add missing columns
  add_missing_columns <- function(df, columns) {
    for (col in columns) {
      if (!col %in% colnames(df)) {
        df[[col]] <- NA
      }
    }
    return(df)
  }
  
  df1 <- add_missing_columns(df1, columns)
  df2 <- add_missing_columns(df2, columns)
  df3 <- add_missing_columns(df3, columns)
  df4 <- add_missing_columns(df4, columns)
  df5 <- add_missing_columns(df5, columns)
  
  # Combine all data frames into one
  results_df <- rbind(df1, df2, df3, df4, df5)
  rownames(results_df) <- c("Model 1: Poly + Season", "Model 2: Log Transformed", "Model 3: Lag + Trend", "Model 4: Lag + Poly Trend", "Model 5: ARIMA")
  #print(results_df)
  
  png(paste0("ncsu_colab/plots/AvsP.", state, ".png"), width = 960, height = 540)
  # Plotting actual vs predicted for all models including ARIMA
  predict_lag_polytrend_final <- c(rep(NA, 12), predict_lag_polytrend)
  predict_lag_trend_final <- c(rep(NA, 12), predict_lag_trend)
  ylim_range <- range(c(data$cases, predict_ployseas, predicted_log_cases, predict_lag_trend, predict_lag_polytrend, arima_fitted), na.rm = TRUE)
  plot(data$cases, type = 'l', col = 'black', ylim = ylim_range, ylab=paste( "Cases"), xlab='Time', main=paste('Original vs. Predicted', 'Cases'))
  lines(predict_ployseas, col = 'blue', lty = 2)
  lines(predicted_log_cases, col = 'red', lty = 3)
  lines(predict_lag_trend_final, col = 'green', lty = 4)
  lines(predict_lag_polytrend_final, col = 'orange', lty = 5)
  lines(arima_fitted, col = 'purple', lty = 6)
  
  legend("topleft", legend = c(paste(state, "Actual"), "Predicted (Poly+Season)", "Predicted (Log+Poly+Season)", "Predicted (Lag 12)", "Predicted (Lag + Poly Trend)", "ARIMA Fitted"), col = c("black", "blue", "red", "green", "orange", "purple"), lty = 1:6, cex = 0.75, bty = "n")  
  dev.off()
  # Display p, d, q values
  arima_order <- arima_model$arma
  cat("ARIMA Model Order: (p, d, q) = (", arima_order[1], ", ", arima_order[6], ", ", arima_order[2], ")\n")
  
  return(list(residuals1, residuals_log, residuals_lag, residuals_lag_polytrend, arima_residuals))
}

plot_all <- function(data, mapes_all, residuals, state){
  mapes_all$MAPE1 <- as.numeric(as.character(mapes_all$MAPE1))
  mapes_all$MAPE2 <- as.numeric(as.character(mapes_all$MAPE2))
  mapes_all$MAPE3 <- as.numeric(as.character(mapes_all$MAPE3))
  mapes_all$MAPE4 <- as.numeric(as.character(mapes_all$MAPE4))
  mapes_all$MAPE5 <- as.numeric(as.character(mapes_all$MAPE5))
  
  png(paste0("ncsu_colab/plots/MAPE.", state, ".png"), width = 960, height = 540)  
  plot(mapes_all$MAPE1, type = 'l', col = 'blue', ylim = range(c(mapes_all$MAPE1, mapes_all$MAPE2, mapes_all$MAPE3, mapes_all$MAPE4, mapes_all$MAPE5), na.rm = TRUE), ylab = "MAPE", xlab = "Iteration", main = paste("MAPE for Rolling Forecasts for", state))
  lines(mapes_all$MAPE2, col = 'red')
  lines(mapes_all$MAPE3, col = 'green')
  lines(mapes_all$MAPE4, col = 'orange')
  lines(mapes_all$MAPE5, col = 'purple')
  legend("topright", legend = c("Model 1: Poly+Season", "Model 2: Log+Poly+Season", "Model 3: Lag 12 Diff", "Model 4: Lag 12_diff_poly", "Model 5: ARIMA"), col = c("blue", "red", "green", "orange", "purple"), lty = 1)
  dev.off()
  
  # Residuals, ACF, and PACF plots for each model
  par(mfrow = c(3, 1))
  
  # Model 1: Poly + Season
  plot(residuals[[1]], main = paste("Residuals from Model 1: Poly + Season for", state), ylab = "Residuals")
  acf(residuals[[1]], main = paste("ACF of Residuals from Model 1 for", state), lag.max = 35)
  pacf(residuals[[1]], main = paste("PACF of Residuals from Model 1 for", state), lag.max = 35)
  
  # Model 2: Log + Poly + Season
  plot(residuals[[2]], main = paste("Residuals from Model 2: Log + Poly + Season for", state), ylab = "Residuals")
  acf(residuals[[2]], main = paste("ACF of Residuals from Model 2 for", state), lag.max = 35)
  pacf(residuals[[2]], main = paste("PACF of Residuals from Model 2 for", state), lag.max = 35)
  
  # Model 3: Lag + Trend
  plot(residuals[[3]], main = paste("Residuals from Model 3: Lag + Trend for", state), ylab = "Residuals")
  acf(residuals[[3]], main = paste("ACF of Residuals from Model 3 for", state), lag.max = 35)
  pacf(residuals[[3]], main = paste("PACF of Residuals from Model 3 for", state), lag.max = 35)
  
  # Model 4: Lag + Poly Trend
  plot(residuals[[4]], main = paste("Residuals from Model 4: Lag + Poly Trend for", state), ylab = "Residuals")
  acf(residuals[[4]], main = paste("ACF of Residuals from Model 4 for", state), lag.max = 35)
  pacf(residuals[[4]], main = paste("PACF of Residuals from Model 4 for", state), lag.max = 35)
  
  # Model 5: ARIMA
  plot(residuals[[5]], main = paste("Residuals from Model 5: ARIMA for", state), ylab = "Residuals")
  acf(residuals[[5]], main = paste("ACF of Residuals from Model 5 for", state), lag.max = 35)
  pacf(residuals[[5]], main = paste("PACF of Residuals from Model 5 for", state), lag.max = 35)
  
  par(mfrow = c(1, 1))
  
}

run <- function(state, method = c("lm", "rlm")){
  method <- match.arg(method)
  
  data = preprocess_data(data, state)
  
  print(length(data$cases - 3 - 35))
  
  x <- rolling_forecast(data, method)
  
  print(x[[1]])
  
  y <- insample_fit(data, x[[1]], method, state)
  
  plot_all(data, x[[2]], y, state)
  
  print(x[[2]])
}

for (state in names[2:29]){
  print(state)
  run(state, "lm")
}
