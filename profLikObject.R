## 10/01/2025
## some experimental code to handle profile likelihoods

#' Create a Log-Likelihood Function with Interpolation and Extrapolation
#'
#' This function takes a data frame representing a profile log-likelihood and
#' returns a function that can calculate the log-likelihood and its gradient
#' at any given parameter value. It uses cubic spline interpolation for points
#' within the data range and quadratic extrapolation for points outside the range.
#'
#' @param lik A data frame with two columns:
#'            - `point`: The parameter value (beta).
#'            - `value`: The corresponding log-likelihood value.
#'
#' @return A function `loglikelihood(beta)` that takes a numeric value or vector `beta`
#'         and returns a list containing two vectors of the same length as `beta`:
#'         - `value`: The calculated log-likelihood values.
#'         - `gradient`: The calculated gradient (first derivative) values.
#'
#' @examples
#' # 1. Generate sample profile log-likelihood data
#' set.seed(42)
#' beta_points <- seq(0, 10, by = 1)
#' # A simple quadratic function centered at beta = 5, with some noise
#' loglik_values <- -((beta_points - 5)^2) + rnorm(length(beta_points), 0, 0.1)
#' lik_data <- data.frame(point = beta_points, value = loglik_values)
#'
#' # 2. Create the loglikelihood function from the data
#' loglikelihood_func <- create_loglikelihood_function(lik_data)
#'
#' # 3. Test the function with a single value and a vector
#' cat("Test with a single value (beta = 4.5):\n")
#' print(loglikelihood_func(4.5))
#'
#' cat("\nTest with a vector of values:\n")
#' print(loglikelihood_func(c(-2, 4.5, 12)))
#'
create_loglikelihood_function <- function(lik) {
  
  # --- 1. Input Validation ---
  if (!is.data.frame(lik) || !all(c("point", "value") %in% names(lik))) {
    stop("Input 'lik' must be a data frame with 'point' and 'value' columns.")
  }
  if (nrow(lik) < 4) {
    stop("Input 'lik' must have at least 4 rows for stable interpolation and extrapolation.")
  }
  
  # --- 2. Prepare Data ---
  # Sort data by the parameter point to ensure correct interpolation
  lik_sorted <- lik[order(lik$point), ]
  
  # Define the boundaries for interpolation vs. extrapolation
  min_beta <- min(lik_sorted$point)
  max_beta <- max(lik_sorted$point)
  
  # --- 3. Set up Interpolation ---
  # Create a cubic spline interpolation function. This object can also compute derivatives.
  spline_interpolator <- splinefun(lik_sorted$point, lik_sorted$value)
  
  # --- 4. Set up Extrapolation ---
  # Fit a quadratic model (y = c0 + c1*x + c2*x^2) to the first 3 points for left extrapolation
  extrap_data_left <- head(lik_sorted, 3)
  model_left <- lm(value ~ point + I(point^2), data = extrap_data_left)
  coefs_left <- coef(model_left)
  
  # Fit a quadratic model to the last 3 points for right extrapolation
  extrap_data_right <- tail(lik_sorted, 3)
  model_right <- lm(value ~ point + I(point^2), data = extrap_data_right)
  coefs_right <- coef(model_right)
  
  
  # --- 5. Define and Return the Main Function ---
  # This is the function that will be returned to the user. It "closes over"
  # the prepared data and models from the parent environment.
  loglikelihood <- function(beta) {
    
    # Ensure beta is a numeric vector
    if (!is.numeric(beta)) {
      stop("Input 'beta' must be a numeric value or vector.")
    }
    
    # Initialize result vectors
    n <- length(beta)
    ll_value <- numeric(n)
    ll_gradient <- numeric(n)
    
    # Create logical indices for each case
    idx_interp <- beta >= min_beta & beta <= max_beta
    idx_extrap_left <- beta < min_beta
    idx_extrap_right <- beta > max_beta
    
    # --- Interpolation ---
    # Process values within the original data range
    if (any(idx_interp)) {
      beta_interp <- beta[idx_interp]
      ll_value[idx_interp] <- spline_interpolator(beta_interp)
      ll_gradient[idx_interp] <- spline_interpolator(beta_interp, deriv = 1)
    }
    
    # --- Left Extrapolation ---
    # Process values below the original data range
    if (any(idx_extrap_left)) {
      beta_extrap_left <- beta[idx_extrap_left]
      new_data_left <- data.frame(point = beta_extrap_left)
      ll_value[idx_extrap_left] <- predict(model_left, newdata = new_data_left)
      ll_gradient[idx_extrap_left] <- coefs_left[2] + 2 * coefs_left[3] * beta_extrap_left
    }
    
    # --- Right Extrapolation ---
    # Process values above the original data range
    if (any(idx_extrap_right)) {
      beta_extrap_right <- beta[idx_extrap_right]
      new_data_right <- data.frame(point = beta_extrap_right)
      ll_value[idx_extrap_right] <- predict(model_right, newdata = new_data_right)
      ll_gradient[idx_extrap_right] <- coefs_right[2] + 2 * coefs_right[3] * beta_extrap_right
    }
    
    # Return a named list with the result vectors
    return(list(value = ll_value, gradient = ll_gradient))
  }
  
  return(loglikelihood)
}


#### --- Example Usage --- ####

# 1. Generate sample profile log-likelihood data
set.seed(42)
beta_points <- seq(0, 10, by = 1)
loglik_values <- -((beta_points - 5)^2) + rnorm(length(beta_points), 0, 0.1)
lik_data <- data.frame(point = beta_points, value = loglik_values)

# 2. Create the loglikelihood function from the data
loglikelihood_func <- create_loglikelihood_function(lik_data)

# 3. Test the function
cat("--- Testing the Generated Function ---\n")
cat("Test with a single value (beta = 4.5):\n")
print(loglikelihood_func(4.5))

cat("\nTest with a vector of values:\n")
print(loglikelihood_func(c(-2, 4.5, 12)))


# 4. Visualize the results
# Generate a sequence of beta values for a smooth plot
plot_betas <- seq(-5, 15, by = 0.1)

# Calculate log-likelihood and gradient for all betas in one vectorized call
plot_results <- loglikelihood_func(plot_betas)

# Combine results into a data frame for plotting
plot_df <- data.frame(
  beta = plot_betas,
  value = plot_results$value,
  gradient = plot_results$gradient
)

# Set up plot layout
par(mfrow = c(2, 1), mar = c(4.1, 4.1, 3.1, 1.1))

# Plot the log-likelihood function
plot(lik_data$point, lik_data$value,
     pch = 19, col = "black", cex = 1.2,
     xlab = "Parameter (beta)", ylab = "Log-Likelihood",
     main = "Interpolated and Extrapolated Log-Likelihood",
     xlim = range(plot_df$beta),
     ylim = range(plot_df$value, lik_data$value, na.rm = TRUE)
)
lines(plot_df$beta, plot_df$value, col = "dodgerblue", lwd = 2)
abline(v = c(min(lik_data$point), max(lik_data$point)), lty = 2, col = "firebrick")
legend("bottomright", bty = "n",
       legend = c("Original Data Points", "Fitted Function", "Extrapolation Boundary"),
       col = c("black", "dodgerblue", "firebrick"),
       lty = c(NA, 1, 2), pch = c(19, NA, NA))

# Plot the gradient
plot(plot_df$beta, plot_df$gradient,
     type = "l", col = "darkgreen", lwd = 2,
     xlab = "Parameter (beta)", ylab = "Gradient",
     main = "Gradient of the Log-Likelihood Function")
abline(v = c(min(lik_data$point), max(lik_data$point)), lty = 2, col = "firebrick")
abline(h = 0, lty = 1, col = "grey50")
legend("topright", bty = "n",
       legend = c("Gradient", "Extrapolation Boundary"),
       col = c("darkgreen", "firebrick"),
       lty = c(1, 2))

