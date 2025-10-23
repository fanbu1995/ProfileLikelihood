## 10/01/2025
## try to handle batch profile likelihoods, as downloaded from OHDSI server

#' Create a Batch Log-Likelihood Function from a data frame of profiles
#'
#' This function processes a data frame where each row contains a complete
#' profile log-likelihood curve, identified by a combination of IDs. The point and
#' value data are stored as semicolon-separated strings. It returns a single
#' function that can evaluate all profiles for a given beta value or vector.
#'
#' @param lik_batch_df A data frame with the following columns:
#'        - At least one of `exposure_id`, `outcome_id`, `period_id` to identify profiles.
#'        - `point`: A string of semicolon-separated numeric parameter values.
#'        - `value`: A string of semicolon-separated numeric log-likelihood values.
#'
#' @return A function `loglikelihoodBatch(beta)` that takes a numeric value or vector `beta`
#'         and returns a named list. Each element of the list corresponds to one
#'         profile and contains the log-likelihood `value` and `gradient`.
#'
#' @examples
#' # 1. Create sample batch data (note 'exposure_id' is missing)
#' batch_data <- data.frame(
#'   outcome_id = c("outA", "outB", "outC"),
#'   period_id = c(2020, 2021, 2020),
#'   point = c("0;1;2;3;4", "5;6;7;8;9", "2;3;4;5;6"),
#'   value = c("-8;-3;-1;-3;-8", "-9;-4;-2;-4;-9", "-4;-1;0;-1;-4")
#' )
#'
#' # 2. Create the batch processing function
#' loglik_batch_func <- create_loglikelihood_batch_function(batch_data)
#'
#' # 3. Test with a single beta value
#' cat("\n--- Testing Batch Function (Single Beta) ---\n")
#' print(loglik_batch_func(2.5))
#'
#' # 4. Test with a vector of beta values
#' cat("\n--- Testing Batch Function (Vector of Betas) ---\n")
#' print(loglik_batch_func(c(0.5, 7.5)))
#'
create_loglikelihood_batch_function <- function(lik_batch_df) {
  
  # --- 1. Input Validation ---
  potential_id_cols <- c("exposure_id", "outcome_id", "period_id")
  id_cols_present <- intersect(potential_id_cols, names(lik_batch_df))
  
  if (length(id_cols_present) == 0) {
    stop(paste("Input data frame must contain at least one of:", paste(potential_id_cols, collapse = ", ")))
  }
  
  required_data_cols <- c("point", "value")
  if (!all(required_data_cols %in% names(lik_batch_df))) {
    stop("Input data frame is missing 'point' and/or 'value' columns.")
  }
  
  # --- 2. Pre-process and Create All Interpolator Functions (Efficiently) ---
  
  # Vectorized creation of unique keys for each profile
  keys <- do.call(paste, c(lik_batch_df[id_cols_present], sep = "_"))
  
  # Vectorized parsing of string data into lists of numeric vectors
  points_list <- lapply(strsplit(lik_batch_df$point, ";"), as.numeric)
  values_list <- lapply(strsplit(lik_batch_df$value, ";"), as.numeric)
  
  # Use mapply to iterate over the parsed lists and create an interpolator for each profile
  interpolator_list <- mapply(function(p, v, key) {
    lik_single <- data.frame(point = p, value = v)
    
    # Check for parsing errors or insufficient data for this specific profile
    if(any(is.na(p)) || any(is.na(v)) || nrow(lik_single) < 4) {
      warning(paste("Skipping profile", key, "due to invalid data or insufficient points."))
      return(NULL) # Return NULL for invalid entries that will be filtered out
    }
    
    # If data is valid, create the log-likelihood function
    return(create_loglikelihood_function(lik_single))
    
  }, p = points_list, v = values_list, key = keys, SIMPLIFY = FALSE) # SIMPLIFY=FALSE ensures a list is returned
  
  # Assign the generated keys as names for the list of functions
  names(interpolator_list) <- keys
  
  # Remove any NULL entries that resulted from skipped profiles
  interpolator_list <- interpolator_list[!sapply(interpolator_list, is.null)]
  
  # --- 3. Define and Return the Main Batch Function ---
  loglikelihoodBatch <- function(beta) {
    
    # Apply each pre-computed function to the input beta
    results_list <- lapply(interpolator_list, function(func) {
      func(beta)
    })
    
    return(results_list)
  }
  
  return(loglikelihoodBatch)
}

# --- Example Usage for the Batch Function ---

# 1. Create sample batch data (note 'exposure_id' is missing)
batch_data <- data.frame(
  outcome_id = c("outA", "outB", "outC"),
  period_id = c(2020, 2021, 2020),
  point = c("0;1;2;3;4", "5;6;7;8;9", "2;3;4;5;6"),
  value = c("-8;-3;-1;-3;-8", "-9;-4;-2;-4;-9", "-4;-1;0;-1;-4")
)

# 1b Load example LPs from EUMAEUS study
load("data/CCAE_HC_profiles.RData")
batch_data <- exampleLPs[1:5, c("exposure_id", "outcome_id", "period_id", "point", "value")]

# 2. Create the batch processing function
loglik_batch_func <- create_loglikelihood_batch_function(batch_data)

# 3. Test with a single beta value
cat("\n\n--- Testing Batch Function (Single Beta) ---\n")
print(loglik_batch_func(2.5))

# 4. Test with a vector of beta values
cat("\n--- Testing Batch Function (Vector of Betas) ---\n")
print(loglik_batch_func(c(0.5, 7.5)))

