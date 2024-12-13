# Load necessary libraries
packages <- c("data.table",
              "ggplot2",
              "reshape2",
              "pracma",   # For matrix operations
              "vars",
              "MTS",  
              "stats"    
)

# Install missing packages
packages_to_install <- packages[!packages %in% installed.packages()[,"Package"]]
if(length(packages_to_install)) install.packages(packages_to_install, repos = "https://cloud.r-project.org/")

# Load libraries
lapply(packages, require, character.only = TRUE)

# Set working directory
setwd()

###--------------------------------------------------------------------------###
### Data Preparation
df <- fread("owid-covid-data.csv")
sum(is.na(df))

# Data selection and imputation
# Remove columns with > 50% NA values in data
total_rows <- nrow(df)
columns_to_keep <- sapply(df, function(col) sum(is.na(col)) <= (0.5 * total_rows))
df_filtered <- df[, names(df)[columns_to_keep], with = FALSE]

# Filter for Vietnam
df_vietnam <- df_filtered[location == "Vietnam"]
sum(is.na(df_vietnam))

# Select numeric columns and the 'date' column
covid_df <- df_vietnam[, c(names(df_vietnam)[sapply(df_vietnam, is.numeric)], "date"), with = FALSE]

# Apply nafill to numeric columns only (Forward and Backward Filling)
numeric_cols <- names(covid_df)[sapply(covid_df, is.numeric)]
covid_imputed_fw <- covid_df[, (numeric_cols) := lapply(.SD, function(col) nafill(col, type = "locf")), .SDcols = numeric_cols]
covid_imputed <- covid_imputed_fw[, (numeric_cols) := lapply(.SD, function(col) nafill(col, type = "nocb")), .SDcols = numeric_cols]

# Check if any missing values remain, should be 0
cat("Remaining missing values:", sum(is.na(covid_imputed)), "\n")

# Convert 'date' column to Date format
covid_imputed$date <- as.Date(covid_imputed$date, format = "%m/%d/%Y")

# Remove columns that do not change over time (zero variance columns)
no_variance_cols <- sapply(covid_imputed[, -which(names(covid_imputed) == "date"), with = FALSE], function(x) var(x, na.rm = TRUE) == 0)
covid_imputed <- covid_imputed[, !names(covid_imputed) %in% names(no_variance_cols[no_variance_cols]), with = FALSE]

### Add Total Cases for China and the US to the Vietnam dataset

# Filter data for China and the US
df_china <- df_filtered[location == "China"]
df_us <- df_filtered[location == "United States"]

# Extract total cases data for China and the US, ensure 'date' is Date type
df_china_cases <- df_china[, .(date, total_cases_china = total_cases)]
df_us_cases <- df_us[, .(date, total_cases_us = total_cases)]

# Ensure the 'date' column is of Date type in all datasets
covid_imputed$date <- as.Date(covid_imputed$date, format = "%m/%d/%Y")
df_china_cases$date <- as.Date(df_china_cases$date, format = "%m/%d/%Y")
df_us_cases$date <- as.Date(df_us_cases$date, format = "%m/%d/%Y")

# Merge the data for China and the US into the Vietnam dataset by 'date'
covid_imputed <- merge(covid_imputed, df_china_cases, by = "date", all.x = TRUE)
covid_imputed <- merge(covid_imputed, df_us_cases, by = "date", all.x = TRUE)

# Check the final dataset with the added columns
head(covid_imputed)

# Transform total cases in Vietnam to log_scale
# Add a log-transformed column for Vietnam's total cases
covid_imputed[, total_cases := log(total_cases)]

###----------------------------------------------------------------------------###

### Correlation matrix

# Remove the 'date' column and select numeric columns only
numeric_data <- covid_imputed[, lapply(.SD, as.numeric), .SDcols = setdiff(names(covid_imputed), c("date", "total_cases"))]

# Compute the correlation matrix
variable_mapping <- data.table(Variable = names(numeric_data), Index = as.character(seq_along(numeric_data)))
numeric_data <- setnames(numeric_data, old = names(numeric_data), new = variable_mapping$Index)
cor_matrix <- cor(numeric_data, use = "complete.obs")
cor_melt <- as.data.table(as.table(cor_matrix))
setnames(cor_melt, old = c("V1", "V2", "N"), new = c("Var1", "Var2", "value"))

# Plot the correlation heatmap with values displayed
ggplot(cor_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10)) +
  labs(title = "Correlation Heatmap", x = "Variable Index", y = "Variable Index")

# Print the mapping of indexes to variable names
cat("Variable Index Mapping:\n")
print(variable_mapping)

# Identify highly correlated variables
high_corr <- cor_melt[abs(value) > 0.9 & Var1 != Var2]
print(high_corr)

# Keep only one variable from each pair
to_remove <- unique(high_corr[, ifelse(Var1 < Var2, Var2, Var1)])  # Remove the "second" variable alphabetically

# Remove the highly correlated variables
numeric_data_clean <- numeric_data[, !to_remove, with = FALSE]

# Verify remaining variables
cat("Remaining variables after removing high correlations:", names(numeric_data_clean), "\n")

cor_matrix_clean <- cor(numeric_data_clean, use = "complete.obs")
numeric_data_clean$total_cases <- covid_imputed$total_cases
numeric_data_clean$date <- covid_imputed$date
covid_imputed <- numeric_data_clean

# Plot the correlation heatmap

# Convert the correlation matrix to long format
cor_matrix_melt <- melt(cor_matrix_clean)
colnames(cor_matrix_melt) <- c("Var1", "Var2", "value")

###---------------------------------------------------------------------------###
### Data Exploration
# Plot time series for total cases in Vietnam
ggplot(data = covid_imputed, aes(x = date, y = total_cases)) +
  geom_line(color = "blue") + # Correct: Line graph
  labs(
    title = "Time Series of Total COVID-19 Cases in Vietnam (log scale)",
    x = "Date",
    y = "Total Cases"
  )

###---------------------------------------------------------------------------###
# Train and Test Data
# Define the cutoff date
cutoff_date <- as.Date("2022-03-31")

# Split the data into training and testing sets
train_data <- covid_imputed[date <= cutoff_date]
test_data <- covid_imputed[date > cutoff_date]

# Check the sizes of the training and testing sets
cat("Training data rows:", nrow(train_data), "\n")
cat("Testing data rows:", nrow(test_data), "\n")

# View a summary of the split datasets
summary(train_data)
summary(test_data)

###---------------------------------------------------------------------------###
# Choose K using eigenvalue ratio .. (2013)
# Exclude both 'total_cases' (dependent variable) and 'date' column
numeric_cols <- names(train_data)[sapply(train_data, is.numeric)]
independent_vars <- setdiff(numeric_cols, c("total_cases", "date"))  # Exclude 'total_cases' and 'date'
numeric_train_data <- train_data[, ..independent_vars]

# Perform PCA on the independent variables
pca_result <- prcomp(numeric_train_data, center = TRUE, scale. = TRUE)
eigenvalues <- pca_result$sdev^2  # Extract eigenvalues

# Calculate eigenvalue ratios
eigen_ratios <- eigenvalues[-length(eigenvalues)] / eigenvalues[-1]

# Find the optimal number of principal components (K) using the maximum eigenvalue ratio
optimal_K_eigen <- which.max(eigen_ratios)
cat("Optimal number of principal components (K) using Eigenvalue Ratios:", optimal_K_eigen, "\n")

# Plot the eigenvalues to visualize their contribution
eigen_df <- data.table(
  Component = 1:length(eigenvalues),
  Eigenvalue = eigenvalues
)

ggplot(eigen_df, aes(x = Component, y = Eigenvalue)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(
    title = "Scree Plot of Eigenvalues",
    x = "Principal Component",
    y = "Eigenvalue"
  ) +
  theme_minimal()

# STEP 1: Compute K first principal components
# Perform PCA
pca_result <- prcomp(numeric_train_data, center = TRUE, scale. = TRUE)

# Extract the first K principal components
K <- 4 # Replace this with the optimal K obtained earlier
principal_components <- as.data.table(pca_result$x[, 1:K])

# STEP 2: Estimate the coefficients lambda and VAR coefficient matrices
# Add the dependent variable 'total_cases' to the principal components data
dependent_variable <- train_data$total_cases
regression_data <- cbind(principal_components, total_cases = dependent_variable)

# Perform linear regression to estimate lambda_y
formula <- as.formula(paste("total_cases ~", paste(colnames(principal_components), collapse = " + ")))
linear_model <- lm(formula, data = regression_data)

# Output the estimated coefficients (lambda_y)
lambda_y <- coef(linear_model)
cat("Estimated coefficients (lambda_y):\n")
print(lambda_y)

# Use VARorder() to determine the optimal lag
principal_components_ts <- ts(principal_components)
var_order <- VARorder(principal_components_ts, maxp = 10)
print(var_order)

# Fit a VAR model on the principal components
principal_components_ts <- ts(principal_components, start = 1, frequency = 1)
optimal_lag <- 10  # Select the optimal lag (e.g., based on AIC)
var_model <- MTS::VAR(principal_components_ts, p = optimal_lag)

# Print the summary of the VAR model
summary(var_model)
###---------------------------------------------------------------------------###
## Forecast 
# Set the forecast horizon
h <- 14

# Forecast the principal components using the VAR model
forecasted_components <- VARpred(var_model, h )

forecasted_principal_components <- forecasted_components$pred

# Extract the intercept and slopes from lambda_y
intercept <- lambda_y[1]
slopes <- lambda_y[-1]

# Compute the forecast for total_cases
forecasted_total_cases <- intercept + as.matrix(forecasted_principal_components) %*% slopes
cat("Forecasted total_cases for horizon h:\n")
print(forecasted_total_cases)

# Combine historical and forecasted values
historical_values <- train_data$total_cases
forecasted_values <- c(historical_values, forecasted_total_cases)

# Create a time-series object
forecasted_ts <- ts(forecasted_values, start = 1)

plot(forecasted_ts, type = "o", col = "blue", main = "Forecasted Total Cases",
     xlab = "Time", ylab = "Total Cases")
abline(v = length(historical_values), col = "red", lty = 2)  # Mark forecast start

###---------------------------------------------------------------------------###
# ERROR METRICS
# Plot the historical and forecasted values

# Extract the actual values for the same forecast horizon
actual_total_cases <- test_data$total_cases[1:h]  # First `h` values from the test set

# Compute the residuals (differences between actual and forecasted values)
residuals <- actual_total_cases - forecasted_total_cases

# Compute RMSE
rmse <- sqrt(mean(residuals^2))
cat("RMSE of the forecast:", rmse, "\n")

# View comparison of forecasted vs actual values
comparison <- data.frame(
  Date = test_data$date[1:h],
  Actual = actual_total_cases,
  Forecasted = forecasted_total_cases
)
print(comparison)
