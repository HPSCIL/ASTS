#######################################################################
##'# R version 4.3.1 (2023-06-16 ucrt) -- "Beagle Scouts"
# Platform: x86_64-w64-mingw32/x64 (64-bit)
##'
##'This code aims to provide the source code for a spatiotemporal adaptive sampling algorithm model.
##'The model is applicable to research aimed at minimizing  cumulative prediction errors.
##'The minimum requirements for model construction are: projection coordinates utmX, utmY (used for building the inla-spde mesh) and the target variable.
##'If there are covariates involved in the model building, they can also be applied (just modify the expression as needed).
##'
##'This software package is compiled and executed in a Windows environment. 
##'The Mac environment requires independent configuration of packages such as INLA 
##'The INLA package and related papers referenced for use are as follows.'https://www.r-inla.org/
##' @author JunfengGu  v1.0

#######################################################################
#############  cube1 load library                                  
#######################################################################
library(ggplot2)
library(tidyr)
library(dplyr)
library(INLA)
library(geosample)
library(magrittr)
library(purrr)
library("GA")
#######################################################################
#############  cube2 load data                                   
#######################################################################
load("origin.data.rdata")

head(data.border)

# UTM_X    UTM_Y
# 1 839.406 3434.387
# 2 841.059 3320.620
# 3 837.062 3445.284
# 4 838.012 3339.756
# 5 881.109 3430.434
# 6 802.687 3435.190

head(data.coordinates.utm)

# Station_ID   UTM_X    UTM_Y
# 1          1 808.460 3469.927
# 2          2 812.737 3470.054
# 3          3 808.585 3465.670
# 4          4 812.863 3465.797
# 5          5 817.142 3465.925
# 6          6 821.420 3466.056

head(data.coordinates.sf)

# Station_ID                 geometry
# 1          1  POINT (808.46 3469.927)
# 2          2 POINT (812.737 3470.054)
# 3          3  POINT (808.585 3465.67)
# 4          4 POINT (812.863 3465.797)
# 5          5 POINT (817.142 3465.925)
# 6          6  POINT (821.42 3466.056)


head(data.pm.all)
# Station_ID      Lon     Lat  UTM_X    UTM_Y month PM2017 PM2018 PM2019 PM2020 PM2021
# 1          1 114.2413 31.3229 808.46 3469.927     1   78.4   75.6   86.2   58.4   60.3
# 2          1 114.2413 31.3229 808.46 3469.927     2   62.9   63.2   67.3   39.2   44.6
# 3          1 114.2413 31.3229 808.46 3469.927     3   56.2   45.6   43.3   35.6   36.3
# 4          1 114.2413 31.3229 808.46 3469.927     4   39.9   39.8   34.7   32.8   30.4
# 5          1 114.2413 31.3229 808.46 3469.927     5   40.1   33.1   30.9   28.4   26.8


#######################################################################
#############  cube3 parameter combination                                  
#######################################################################
#'The settings of the parameters influence the computational efficiency of the model and the final sampling results. 
#'The representations of each parameter are as follows:
#' ini_sample_num:
#' The initial sampling quantity serves as prior information for initializing the overall sampling scheme. 
#' The initial sampling quantity influences the final sampling results and should be adjusted based on the target number of samples. 
#' For example, if the ultimate goal is to have 10 sampling points, setting the initial points to 5 or 6 is more appropriate.
#' While 3 or 4 sampling points can still allow the model to function properly, the lack of sufficient prior information may lead to biased results. 
#' However, this is not a significant issue, as the model's use of multiple parameters will eventually yield optimal results, which can help address some of these concerns from certain perspectives.
#'ini_sample_border_dis: Initialize the distance between sampling points and boundaries to avoid the influence of INLA-SPDE boundary effects on prior information.
#'ini_sample_dis:ini_sample_border_dis: Initializes the minimum distance between sampling points. 
#'               This parameter allows the initial model to incorporate better prior information, preventing the sampling set from being overly dense.
#'
#'NOTE:
#'(ini_sample_num,ini_sample_border_dis,ini_sample_dis)
#'These three parameters are necessary for the initial sampling process and mainly influence the prior information.
#' Better prior information will inevitably yield a superior initial plan. 
#' The greater the difference between the initial sampling points and the target sampling points, the less impact the initial sampling points will have on the final sampling scheme. 
#' Additionally, adjustments can be made based on the research objectives.
#' 
#' 
#' The prior.range_max and alpha parameters are used to construct the spde component in INLA-SPDE.
#' INLA-SPDE create an inla.spde2 model object for a Matern model, using a PC prior for the parameters.
#' prior.range_max : the max of spatial range of the random field
#' alpha:Fractional operator order, 0<alpha≤2 supported, for v=alpha-d/2>0
#' 
#' sampleDis:Minimum distance between two sampling points during the ASTS process.
#' borderDis:Minimum distance from sampling points to boundaries in ASTS
#' 
#' 
#' You may not know exactly how to set the parameters to achieve optimal results from the model. 
#' Below is a method for establishing a set of parameters that can yield the best sampling scheme and optimal parameters. 
#' However, it is important to note that this is just an example; specific boundary values need to be established based on the scale of the research.


ini_sample_num <- 5
target_samplesize <- 10

data.pm <- data.pm.all[, c("Station_ID", 
                           "Lon", 
                           "Lat", 
                           "UTM_X", 
                           "UTM_Y", 
                           "month", 
                           "PM2017")]


alpha_value <- c("1","3/2","2")
prior.range_max_value <- seq(20, 30, by = 5)
sampleDis_value <- seq(10, 20, by = 2)
borderDis_value <- seq(10,15,by=1)


#NOte!!! :Too many combinations mean it will be more time-consuming.
para.combinations <- expand.grid(alpha = alpha_value,
                            prior.range_max = prior.range_max_value,
                            sampleDis = sampleDis_value,
                            borderDis=borderDis_value)


#######################################################################
#############  cube4 load ASTS                                  
#######################################################################
n_stations <- length(data.coordinates.utm$Station_ID)
n_data <- length(data.pm.2017$Station_ID) 
n_time <- as.integer(n_data/n_stations)
data.stationID <- unique(data.coordinates.utm$Station_ID)

#Each parameter's sampling scheme under different iteration requirements is directly stored in a TXT file. 
# For example, if the target sample size is 10 and there are 5 parameter sets, with each set having 5 iterations, then the file name could be sample10_3.5.
# Here, 10_3 indicates the third parameter space, while 5 denotes the fifth iterations within that space.
filename <-"sample_10_" 

# The variable 'outpath' represents the folder where all the txt files from this experiment are stored.
outpath <- "txtOut/sample10"


# Start the timer
start_time <- Sys.time()
# Execute the main function with specified parameters
asts_start(para.combinations, iterations = 5)
# Capture the end time right after the function execution
end_time <- Sys.time()
# Calculate the elapsed execution time
execution_time <- end_time - start_time
# Print the execution time in a readable format
cat("Execution Time:", execution_time, "\n")


#######################################################################
#############  cube6 Extract the optimal solution                           
#######################################################################

# Extracting the optimal results requires only the input of the results folder, regardless of how many txt files are within it; 
# the system will automatically find the best solution.

extract_result("txtOut/sample10")
# The function returns content as follows：
# File with minimum RMSE: txtOut/sample10/sample10_1.2.txt 
# RMSE = 0.2212 
# R² = 0.7302 
# MAE = 0.149 
# MSE = 0.0489 
# SampleSize = 10 
# SampleIDs = SampleID1=31, SampleID2=136, SampleID3=182, SampleID4=192, SampleID5=299, SampleID6=440, SampleID7=437, SampleID8=398, SampleID9=2, SampleID10=386 

# In addition to finding the optimal sampling scheme, 
# sample10_1.2 represents the parameter set where the best results are located, which is the first group.


################################################################################
#############              function collection    ##############################                         
################################################################################ 
#' @author JunFengGu
#' @description asts function collection  
#' @v1.0

# Function: mse
# Description: Calculates the Mean Squared Error (MSE) between actual and predicted numeric vectors.
#
# Parameters:
# - actual: Numeric vector. The actual observed values.
# - predicted: Numeric vector. The predicted values corresponding to the actual values.
#
# Returns: Numeric. The MSE value calculated as the mean of squared differences between actual and predicted values.

mse <- function(actual, predicted) {
  tryCatch({
    if (!is.numeric(actual) || !is.numeric(predicted)) {
      stop("Error: Both 'actual' and 'predicted' must be numeric vectors.")
    }
    if (length(actual) == 0 || length(predicted) == 0) {
      stop("Error: Both 'actual' and 'predicted' cannot be empty.")
    }
    if (length(actual) != length(predicted)) {
      stop("Error: Vectors 'actual' and 'predicted' must be of the same length.")
    }
    
    cat("Calculating MSE for provided actual and predicted vectors.\n")
    mse_value <- mean((actual - predicted)^2)
    cat("MSE calculated successfully.\n")
    return(mse_value)
  }, error = function(e) {
    cat("An error occurred in MSE calculation: ", e$message, "\n")
    return(NA)
  })
}

# Function: mae
# Description: Calculates the Mean Absolute Error (MAE) between actual and predicted numeric vectors.
#
# Parameters:
# - actual: Numeric vector. The actual observed values.
# - predicted: Numeric vector. The predicted values corresponding to the actual values.
#
# Returns: Numeric. The MAE value calculated as the mean of absolute differences between actual and predicted values.

mae <- function(actual, predicted) {
  tryCatch({
    if (!is.numeric(actual) || !is.numeric(predicted)) {
      stop("Error: Both 'actual' and 'predicted' must be numeric vectors.")
    }
    if (length(actual) == 0 || length(predicted) == 0) {
      stop("Error: Both 'actual' and 'predicted' cannot be empty.")
    }
    if (length(actual) != length(predicted)) {
      stop("Error: Vectors 'actual' and 'predicted' must be of the same length.")
    }
    
    cat("Calculating MAE for provided actual and predicted vectors.\n")
    mae_value <- mean(abs(actual - predicted))
    cat("MAE calculated successfully.\n")
    return(mae_value)
  }, error = function(e) {
    cat("An error occurred in MAE calculation: ", e$message, "\n")
    return(NA)
  })
}

# Function: calculate_r2_score
# Description: Computes the R² (coefficient of determination) score for actual and predicted values.
#
# Parameters:
# - y_actual: Numeric vector. The actual observed values.
# - y_predicted: Numeric vector. The predicted values corresponding to the actual values.
#
# Returns: Numeric. The R² score indicating the proportion of variance in the dependent variable explained by the independent variable(s).

calculate_r2_score <- function(y_actual, y_predicted) {
  tryCatch({
    if (!is.numeric(y_actual) || !is.numeric(y_predicted)) {
      stop("Error: Both 'y_actual' and 'y_predicted' must be numeric vectors.")
    }
    if (length(y_actual) == 0 || length(y_predicted) == 0) {
      stop("Error: Both 'y_actual' and 'y_predicted' cannot be empty.")
    }
    if (length(y_actual) != length(y_predicted)) {
      stop("Error: Vectors 'y_actual' and 'y_predicted' must be of the same length.")
    }
    
    avr_y_actual <- mean(y_actual)
    ss_total <- sum((y_actual - avr_y_actual)^2)
    ss_residuals <- sum((y_actual - y_predicted)^2)
    
    if (ss_total == 0) {
      stop("Total sum of squares is zero; cannot calculate R².")
    }
    
    cat("Calculating R² score for provided actual and predicted values.\n")
    r2 <- 1 - ss_residuals / ss_total
    cat("R² calculated successfully.\n")
    return(r2)
  }, error = function(e) {
    cat("An error occurred in R² calculation: ", e$message, "\n")
    return(NA)
  })
}

# Function: cal_dis
# Description: Calculates the Euclidean distance between two points in 2D space.
#
# Parameters:
# - x1, y1: Numeric. The coordinates of the first point.
# - x2, y2: Numeric. The coordinates of the second point.
#
# Returns: Numeric. The Euclidean distance between the two points.

cal_dis <- function(x1, y1, x2, y2) {
  tryCatch({
    if (!is.numeric(x1) || !is.numeric(y1) || !is.numeric(x2) || !is.numeric(y2)) {
      stop("Error: All inputs to the distance function must be numeric.")
    }
    if (length(x1) != 1 || length(y1) != 1 || length(x2) != 1 || length(y2) != 1) {
      stop("Error: All inputs must be scalars (single numeric values).")
    }
    
    cat("Calculating Euclidean distance between points.\n")
    distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
    cat("Distance calculated successfully.\n")
    return(distance)
  }, error = function(e) {
    cat("An error occurred in distance calculation: ", e$message, "\n")
    return(NA)
  })
}

# Function: calculate_validation_metrics
# Description: Computes various validation metrics for model performance based on predictions from an INLA model.
#
# Parameters:
# - stack: Object. The INLA stack object containing prediction data.
# - result.1: Object. The result object from an INLA model fit, containing predictive results.
# - val.data: Data frame. Contains validation data, including actual observations to compare against.
#
# Returns: List. A named list of validation metrics including DIC, RMSE, MAE, MSE, R², coverage probability, and residuals.

calculate_validation_metrics <- function(stack, result.1, val.data) {
  tryCatch({
    validation.res <- list()
    
    # Extract indices and predictor summaries  
    index_val <- inla.stack.index(stack, "val")$data
    
    predictor_summary <- result.1$summary.linear.predictor  
    hyper_summary <- result.1$summary.hyperpar 
    
    tmp_val.mean <- result.1$summary.linear.predictor[index_val, "mean"]
    tmp_val.sd <- result.1$summary.linear.predictor[index_val, "sd"]
    
    if (length(index_val) != nrow(val.data)) {  
      stop("Mismatch between validation index and provided validation data.")  
    }  
    
    pred_mean <- predictor_summary[index_val, "mean"]  
    pred_sd <- predictor_summary[index_val, "sd"]  
    
    
    # Assign predictions and residuals  
    val.data$pm_val <- pred_mean  
    val.data$res <- abs(val.data$logPM - pred_mean)  
    
    
    # Compute standardised residuals  
    validation.res$res.std <- (val.data$logPM - tmp_val.mean) / sqrt(tmp_val.sd^2 + 1/result.1$summary.hyperpar[1, "mean"])
    validation.res$p <- pnorm(validation.res$res.std)
    validation.res$cover <- mean((validation.res$p > 0.025) & (validation.res$p < 0.975), na.rm = TRUE)
    
    
    
    # Classic metrics  
    validation.res$dic <- result.1[["dic"]][["dic"]]
    validation.res$rmse <- sqrt(mean(val.data$res^2, na.rm = TRUE))
    validation.res$mae <- mae(val.data$logPM, val.data$pm_val)
    validation.res$mse <- mse(val.data$logPM, val.data$pm_val)
    validation.res$r2 <- calculate_r2_score(val.data$logPM, val.data$pm_val)
    validation.res$GroupRho <- result.1$summary.hyperpar["GroupRho for field", ][[1]]
    validation.res$cor <- cor(val.data$logPM, val.data$pm_val, use = "pairwise.complete.obs", method = "pearson")
    
    
    validation.res$val.data <- val.data
    message("Validation metrics calculated successfully.")  

    return(validation.res)
    
  }, error = function(e) {
    cat("An error occurred in calculate_validation_metrics: ", e$message, "\n")
    return(NULL)
  })
}



# Function: analyze_inhibit_sample
# Description: This function attempts to generate an inhibitory sample of spatial points 
# within a given border, ensuring that the generated points are at least a specified 
# minimum distance from the border. It uses a discrete inhibition sampling method.
#
# Parameters:
# - data_border: Data frame. Contains the spatial coordinates (UTM_X, UTM_Y) defining the border.
# - data_utm: Data frame. Contains the spatial coordinates (UTM_X, UTM_Y) of the points to sample from.
# - ini_sampleNUM: Numeric. The initial number of samples to attempt to generate.
# - del: Numeric. The inhibition distance parameter used in the sampling process.
# - dis: Numeric. The minimum distance required between generated samples and the border.
#
# Returns: 
# - A vector of Station_IDs representing the successful inhibitory sample if generated, 
#   or NULL if the process fails after a specified number of attempts.
#
# Notes:
# - The function will attempt to generate a valid sample up to 'maxTryTimes' times. 
#   If unsuccessful, it suggests adjusting the 'dis' parameter.
# - The function includes input validation and error handling to ensure robustness.
analyze_inhibit_sample <- function(ini_sampleNUM, del, dis) {
      tryTimes <- 0
      maxTryTimes <- 1000
      
      tryCatch({
        # Validate inputs
        if (!is.numeric(ini_sampleNUM) || ini_sampleNUM <= 0) {
          stop("Error: 'ini_sampleNUM' must be a positive number.")
        }
        if (!is.numeric(del) || del <= 0) {
          stop("Error: 'del' must be a positive number.")
        }
        if (!is.numeric(dis) || dis <= 0) {
          stop("Error: 'dis' must be a positive number.")
        }
        while (tryTimes < maxTryTimes) {
          # Increment tryTimes
          tryTimes <- tryTimes + 1
          cat("Attempting to generate inhibit sample, try #", tryTimes, "\n")
          # Attempt to generate inhibit sample
          
          inhibit_sample <- discrete.inhibit.sample(obj = data.coordinates.sf, size = ini_sampleNUM, delta = del, plotit = FALSE)
          cat("================inhibit_sample success================ \n")
          
          # Check if inhibit_sample is valid
          if (is.null(inhibit_sample) || length(inhibit_sample[[4]]) == 0) {
            cat("Warning: Inhibit sample generation failed, retrying...\n")
            next
          }
          inhibit_sample_selected_data <- data.coordinates.utm[data.coordinates.utm$Station_ID %in% inhibit_sample[[4]][[1]], ]
          cat("================inhibit_sample_selected_data success================ \n")
          
          # Calculate minimum distances
          min_distances <- sapply(1:nrow(inhibit_sample_selected_data), function(i) {
            point <- inhibit_sample_selected_data[i, ]
            distances <- sqrt((data.border$UTM_X - point$UTM_X)^2 + (data.border$UTM_Y - point$UTM_Y)^2)
            min(distances)
          })
          
          # Check distances
          check_distances <- function(min_distances) {
            if (any(min_distances < dis)) {
              return(FALSE)
            } else {
              return(TRUE)
            }
          }
          result <- check_distances(min_distances)
          if (result) {
            cat("Inhibit sample generation successful.\n")
            plot(data.border)
            points(st_coordinates(inhibit_sample[[4]])[,1], st_coordinates(inhibit_sample[[4]])[,2], col = "blue", pch = 19)
            return(inhibit_sample[[4]][[1]])
          } else {
            cat("Minimum distance requirement not met, retrying...\n")
          }
        }
        
        stop("\nReached maximum tryTimes. Please change 'dis'.")
      }, error = function(e) {
        cat("An error occurred in analyze_inhibit_sample: ", e$message, "\n")
        return(NULL)
      })
    }


# Function: dis_Jug  
# Description:  
#   This function filters out (excludes) validation stations from a dataset if they are within  
#   a specified minimum distance from any estimation station. It is commonly used in spatial   
#   validation processes to ensure that the validation set is spatially independent from the   
#   estimation set by removing validation points that are too close to estimation points.  
#  
# Parameters:  
# - data_val: Data frame. Contains all validation station information.  
#             Must include at least the columns: Station_ID, UTM_X, UTM_Y.  
# - distance: Numeric. The minimum allowable distance (in the same unit as UTM_X/UTM_Y, e.g., meters).  
#             Validation stations within this distance from any estimation station will be excluded.  
# - est_station: Data frame. Contains the spatial coordinates for the estimation stations.  
#                Must include Station_ID, UTM_X, UTM_Y.  
# - val_station: Data frame. Contains the spatial coordinates for the validation stations.  
#                Must include Station_ID, UTM_X, UTM_Y.  
#  
# Returns:  
# - A data frame with the remaining validation stations that are at least 'distance' away  
#   from all estimation stations. The filtered data has the same structure as 'data_val'.  
#  
# Notes:  
# - The function iterates through each estimation station and compares the distance to every  
#   validation station. If the distance is less than the threshold, the validation station is excluded.  
# - The function prints out the number of remaining stations and total exclusions for monitoring.  
# - If, after exclusion, the number of remaining validation stations is negative, a warning is triggered.  
# - It is assumed that a helper function 'cal_dis(x1, y1, x2, y2)' exists for Euclidean distance calculation.  
# - This function has no side effects and does not modify the original data frames.
dis_Jug <- function(data_val, distance, est_station, val_station) {  
  # Generate all pairs of estimation and validation stations  
  idx_est <- rep(1:nrow(est_station), each = nrow(val_station))  
  idx_val <- rep(1:nrow(val_station), times = nrow(est_station))  
  
  # Calculate distances for all pairs  
  dists <- mapply(  
    function(i, j) cal_dis(est_station$UTM_X[i], est_station$UTM_Y[i],  
                           val_station$UTM_X[j], val_station$UTM_Y[j]),  
    idx_est, idx_val  
  )  
  
  # Identify validation stations to exclude  
  to_exclude <- unique(val_station$Station_ID[idx_val[dists < distance]])  
  
  residual <- nrow(val_station) - length(to_exclude)  
  if (residual < 0) {  
    warning("Negative number of residual stations after exclusion.")  
  }  
  
  message("**----residual number ", residual, " Station ----**")  
  message("**----total exclude ", length(to_exclude), " Station ----**")  
  
  # Filter out excluded validation stations  
  new_data <- data_val[!data_val$Station_ID %in% to_exclude, ]  
  return(new_data)  
}  
    


# Function: asts_function
# Description: The function processes geographic and environmental data, supporting iterative sample selection based on an adaptive spatiotemporal sampling strategy, 
# where the integrated nested Laplace approximation (INLA) method is used to fit spatiotemporal models.
# Parameters:
# - txt_name: Character. The name of the text file where results will be written.
# - data.coordinates: Data frame. Contains the spatial coordinates (UTM_X, UTM_Y) of the stations.
# - data.pm: Data frame. Contains particulate matter measurements (PM) and other relevant data.
# - target_samplesize: Numeric. The target sample size for model estimation.
# - inhibit_samplesize: Numeric. The size of the inhibitory sample to be used initially.
# - alpha: Numeric. The range parameter for the SPDE model, controlling spatial smoothness.
# - prior.range_max: Numeric. The maximum prior range value for the SPDE model.
#
# Returns: NULL. The function primarily writes data to a specified text file.
# Each loop iteration prints messages to indicate progress and status of the process.


asts_function <- function(file_name,
                          alpha_value,
                          prior.range_max_value,
                          inhibit_sample) {
  tryCatch({
    # Validate inputs
    stopifnot(  
      is.character(file_name) && nchar(file_name) > 0,  
      is.data.frame(data.coordinates.utm),  
      is.data.frame(data.pm),  
      all(c('UTM_X', 'UTM_Y', 'Station_ID') %in% names(data.coordinates.utm)),  
      is.numeric(target_samplesize) && target_samplesize > 0,  
      is.numeric(ini_sample_num) && ini_sample_num > 0,  
      is.numeric(alpha_value) && alpha_value > 0,  
      is.numeric(prior.range_max_value) && prior.range_max_value > 0  
    ) 
    
    iteration  <- 0
    max_iteration  <- target_samplesize - ini_sample_num + 1
    est_ID <- c()
    
    while (iteration < max_iteration) {
      message("**---- Iteration I =", I, "----**\n")
      if (!is.null(inhibit_sample) && length(inhibit_sample) == ini_sample_num) {
        est_ID <- inhibit_sample
        val_ID <- data.stationID[!data.stationID %in% est_ID]
        inhibit_sample <- NULL
        cat("**---- Inhibit sample success ---**\n")
      } else {
        est_ID <- append(est_ID, result_max2$StationID)
        cat("**---- Append success ---**\n")
      }

      # Working with estimation and validation sets
      est_station <- filter(data.coordinates.utm, Station_ID %in% est_ID)  
      est_data    <- filter(data.pm, Station_ID %in% est_ID)  
      val_station <- filter(data.coordinates.utm, Station_ID %in% val_ID)  
      val_data    <- filter(data.pm, Station_ID %in% val_ID)

      # Normalize covariates
      mean_covariates <- colMeans(data.pm[,4:5], na.rm = TRUE)
      sd_covariates <- apply(data.pm[,4:5], 2, sd, na.rm = TRUE)
      est_data[,4:5] <- scale(est_data[,4:5], mean_covariates, sd_covariates)
      val_data[,4:5] <- scale(val_data[,4:5], mean_covariates, sd_covariates)
      
      # Log-transform PM values
      est_data$logPM <- log(est_data$PM)
      val_data$logPM <- log(val_data$PM)
      
      # Create mesh
      mesh <- inla.mesh.2d(loc = cbind(est_station$UTM_X, est_station$UTM_Y),
                           loc.domain = data.border,
                           max.edge = c(15, 100),
                           min.angle = c(26, 21),
                           cutoff = 5, 
                           plot.delay = NULL)
      cat("Created mesh with", mesh$n, "vertices\n")
      plot(mesh)
      points(x = data.border$UTM_X, y = data.border$UTM_Y, cex = 0.1, col = 'red')
      points(x = est_station$UTM_X, y = est_station$UTM_Y, pch = 17, cex = 1, col = "blue")
      
      # SPDE model
      spde <- inla.spde2.pcmatern(mesh = mesh, 
                                  alpha = alpha_value, 
                                  constr = TRUE,
                                  prior.range = c(prior.range_max_value, 0.01),
                                  prior.sigma = c(3, 0.01))
      
      # Field indices
      field.indices <- inla.spde.make.index("field",
                                            n.spde = spde$n.spde,
                                            n.group = n_time)
      cat("**---- SPDE model and field indices setup success! ---**\n")
      
      # Projection matrices
      A.est <- inla.spde.make.A(mesh, 
                                loc =as.matrix(data.coordinates.utm[est_data$Station_ID, c("UTM_X", "UTM_Y")]),
                                group = est_data$month, n.group = n_time)
      A.val <- inla.spde.make.A(mesh, 
                                loc =as.matrix(data.coordinates.utm[val_data$Station_ID, c("UTM_X", "UTM_Y")]),
                                group = val_data$month, n.group = n_time)
      
      # Stacks
      stack.est <- inla.stack(data = list(logPM = est_data$logPM), 
                              A = list(A.est, 1),
                              effects = list(c(field.indices, list(Intercept = 1)),list(est_data[,4:5])),
                              tag = "est")
      
      stack.val <- inla.stack(data = list(logPM = NA), 
                              A = list(A.val, 1),
                              effects = list(c(field.indices, list(Intercept = 1)),list(val_data[,4:5])),
                              tag = "val")
      
      stack <- inla.stack(stack.est, stack.val)
 
     
      rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.9)))
      
      # Model formula
      formula <- (logPM ~ -1 + Intercept + UTM_X + UTM_Y +
                    f(field, model=spde, group=field.group, 
                      control.group=list(model="ar1", hyper = rprior)))
      
      cat("**---- formula setup success! ---**\n")
      # INLA model fit
      result.1 <- inla(formula, data = inla.stack.data(stack, spde = spde),
                       family = "gaussian",
                       control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
                       control.compute = list(cpo = FALSE, dic = TRUE, config = TRUE, return.marginals.predictor = TRUE),
                       control.inla = list(reordering = "metis", strategy = 'laplace'),
                       keep = FALSE, verbose = TRUE)
      
      cat("INLA model fit summary:\n")
      print(summary(result.1))
      cat("**---- Model fitting success! ---**\n")
      
      # Calculate validation metrics
      validation.res <- calculate_validation_metrics(stack, result.1, val_data)
      validation.res$sampleSize <- length(est_ID)
      validation.res$SampleID <- est_ID 
      file_path <- file.path(outpath, file_name)
      
      output_dir <- dirname(file_path)  
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)  
      if (!file.exists(file_path)) file.create(file_path)

      selected_items <- round(unlist(validation.res[c("sampleSize", "cover", "dic", "rmse", "mae", "mse", "r2", "GroupRho", "cor", "SampleID")]), 4)
      write.table(selected_items, file_path, append = TRUE, col.names = FALSE)
      cat("\n", file = file_path, append = TRUE)
      cat("**---- Results written to file ----**\n")
      # Update inhibition sample
      new_data <- dis_Jug(validation.res$val.data, sampleDis_value, est_station, val_station)
      
      result_max2 <- new_data %>%
        select(Station_ID, month, res) %>%
        group_by(Station_ID) %>%
        summarise(total_res = sum(res)) %>%
        arrange(desc(total_res)) %>%
        slice(1) %>%
        select(StationID = Station_ID)
      I <- I + 1
      cat("**---- Iteration update success! ---**\n")
    }
    cat("**---- Calibration finished successfully ---**\n")
  }, error = function(e) {
    cat("An error occurred in asts_function: ", e$message, "\n")
    return(NULL)
  })
}


# 
# Function: load_asts
# Description
# The load_asts function is designed to perform a series of operations across multiple iterations, generating filenames and processing data samples. 
# It includes input validation to ensure that the parameters are within expected ranges and types. 
# During each iteration, it calls external functions to analyze samples and perform additional computations.
# 
# Parameters
# alpha: Numeric. A weight or proportion used for certain calculations, expected to be between 0 and 1.
# prior_range_max: Numeric. The maximum value for a prior range, used for statistical or probabilistic analysis. Must be positive.
# sampleDis: Numeric. Represents the distance or distribution of the sample, affecting data sampling methods. Must be non-negative.
# borderDis: Numeric. The boundary distance, used to define the limits or constraints of sampling. Must be non-negative.# iterations: Numeric. 
#            The number of attempts or iterations the loop should execute. Must be a positive integer.
# ini_sample_num: Numeric. The initial number of samples, used to determine the starting quantity of data. Must be positive.
# target_samplesize: Numeric. The target sample size, indicating the desired number of samples in the dataset. Must be positive.

load_asts <- function(ini_sample_num,alpha, prior_range_max, sampleDis, borderDis, iterations) {
      # Basic input validation
      if (!is.numeric(prior_range_max) || prior_range_max <= 0) {
        stop("prior_range_max should be a positive numeric value.")
      }
      
      if (!is.numeric(sampleDis) || sampleDis < 0) {
        stop("sampleDis should be a non-negative numeric value.")
      }
      
      if (!is.numeric(borderDis) || borderDis < 0) {
        stop("borderDis should be a non-negative numeric value.")
      }
      
      if (!is.numeric(iterations) || iterations <= 0) {
        stop("iterations should be a positive integer.")
      }
      
      if (!is.numeric(ini_sample_num) || ini_sample_num <= 0) {
        stop("ini_sample_num should be a positive numeric value.")
      }
      
      if (!is.numeric(target_samplesize) || target_samplesize <= 0) {
        stop("target_samplesize should be a positive numeric value.")
      }
      
      # Start the main loop
      for (i in 1:as.integer(iterations)) {
        
        # Generate file name
        file_name <- paste0(filename, i, ".txt")
        
        # Assume analyze_inhibit_sample is a predefined function
        inhibit_sample <- analyze_inhibit_sample(ini_sample_num, sampleDis, as.integer(borderDis))
        
        # Condition to check, e.g. whether inhibit_sample meets criteria
        if (is.null(inhibit_sample)) {
          message("Inhibit sample returned NULL for iteration ", i, ". Skipping.")
          next
        }
        
        asts_function(file_name,
                      as.numeric(alpha),
                      as.numeric(prior_range_max),
                      inhibit_sample)

        # Print the processed file name for confirmation
        print(paste("Processed:", file_name))
      }
  }

# Function: asts_start
# Description: Iterates through combinations of parameters to execute the `load_asts` function multiple times.
#
# Parameters:
# - para.combinations: A data frame. Each row contains a set of parameters required for `load_asts`.
#                      Expected columns include "alpha", "prior.range_max", "sampleDis", and "borderDis".
# - iterations: Numeric. The number of times `load_asts` should be executed for each parameter combination.
#
# Dependencies: Requires the `load_asts` function to be defined elsewhere.
#
# Returns: None. This function performs its operations via side effects, calling `load_asts`.

asts_start <- function(para.combinations, iterations) {
  # Check that iterations is a numeric value
  if (!is.numeric(iterations) || iterations <= 0) {
    stop("Iterations must be a positive numeric value.")
  }

  # Iterate through each parameter combination
  for (i in 1:nrow(para.combinations)) {
    # Extract each set of parameters, ensuring the expected column names are present
    if (!all(c("alpha", "prior.range_max", "sampleDis", "borderDis") %in% names(para.combinations))) {
      stop("para.combinations must contain columns: alpha, prior.range_max, sampleDis, and borderDis.")
    }
    
    alpha <- as.numeric(para.combinations[i, "alpha"])
    prior_range_max <- as.numeric(para.combinations[i, "prior.range_max"])
    sampleDis <- as.numeric(para.combinations[i, "sampleDis"])
    borderDis <- as.numeric(para.combinations[i, "borderDis"])
    
    # Call load_asts function with extracted parameters
    load_asts(ini_sample_num,alpha, prior_range_max , sampleDis, borderDis, iterations)
  }
}


# Function to process text files in a directory and find the file with the smallest RMSE
# Arguments:
#   directory_path: A string indicating the path to the directory containing text files
# Returns:
#   Outputs details of the file with the minimum RMSE, including RMSE, R², MAE, MSE, and sample information
extract_result <- function(directory_path) {
  # Validate that the directory exists
  if (!dir.exists(directory_path)) {
    stop("The specified directory does not exist.")
  }
  
  # List all txt files in the directory
  file_list <- list.files(path = directory_path, pattern = "\\.txt$", full.names = TRUE)
  
  # Check if there are any txt files in the directory
  if (length(file_list) == 0) {
    stop("No .txt files found in the specified directory.")
  }
  
  # Function to read and process each individual file
  process_file <- function(file_path) {
    # Attempt to read the file, handle any reading errors
    file_lines <- tryCatch(readLines(file_path), error = function(e) {
      warning(paste("Error reading file:", file_path, "Skipped."))
      return(NULL)
    })
    
    # Return NULL if file reading fails
    if (is.null(file_lines)) return(NULL)
    
    # Initialize variables to store file data
    data_list <- list()          # Holds all parsed entries
    current_entry <- list()      # Holds current parsing entry
    sample_size_current <- NA    # Tracks current sample size
    
    # Process each line in the file
    for (line in file_lines) {
      line <- trimws(line)  # Remove any leading/trailing whitespace
      
      # Only process non-empty lines
      if (nchar(line) > 0) {
        key_value <- unlist(strsplit(line, " "))  # Split the line into key-value pairs
        key <- key_value[1]                        # Extract the key
        value <- as.numeric(key_value[2])          # Convert value to numeric
        
        # Check for sampleSize key indicating a new entry
        if (key == "\"sampleSize\"") {
          # Save current entry if it exists
          if (!is.na(sample_size_current)) {
            data_list[[as.character(sample_size_current)]] <- current_entry
          }
          # Update to new sample size and reset the current entry
          sample_size_current <- value
          current_entry <- list()
        }
        
        # Add the key-value pair to the current entry
        current_entry[[key]] <- value
      }
    }
    
    # Ensure the last entry is added to the data list
    if (!is.na(sample_size_current)) {
      data_list[[as.character(sample_size_current)]] <- current_entry
    }
    
    # Determine the last entry for the file processed
    last_sample_size <- as.character(max(as.numeric(names(data_list))))
    last_entry <- data_list[[last_sample_size]]
    
    # Return the collected data
    list(
      rmse = last_entry[["\"rmse\""]],
      entry = last_entry,
      sample_size = last_sample_size,
      file_path = file_path
    )
  }
  
  # Process all files in the directory
  results <- map(file_list, process_file)
  
  # Filter out any NULL results due to failed file processing
  results <- results[!sapply(results, is.null)]
  
  # Check if there are valid results
  if (length(results) == 0) {
    stop("No valid data processed from files.")
  }
  
  # Extract RMSE values from the results and identify the file with the minimum RMSE
  rmse_values <- map_dbl(results, "rmse")
  min_rmse_index <- which.min(rmse_values)
  min_rmse_scenario <- results[[min_rmse_index]]
  
  # Output details of the file with the minimum RMSE0
  cat(
    "File with minimum RMSE:", min_rmse_scenario$file_path, "\n",
    "RMSE =", min_rmse_scenario$rmse, "\n",
    "R² =", min_rmse_scenario$entry[["\"r2\""]], "\n",
    "MAE =", min_rmse_scenario$entry[["\"mae\""]], "\n",
    "MSE =", min_rmse_scenario$entry[["\"mse\""]], "\n",
    "SampleSize =", min_rmse_scenario$entry[["\"sampleSize\""]], "\n",
    "SampleIDs =", paste0("SampleID", 1:min_rmse_scenario$entry[["\"sampleSize\""]],
                          "=", unlist(min_rmse_scenario$entry[
                            paste0("\"SampleID", 1:min_rmse_scenario$entry[["\"sampleSize\""]], "\"")]),
                          collapse=", "), "\n"
  )
}







#==============================================================================#
#========================= GA Combined with INLA===============================#
#===== The following is the full code for GA Combined with INLA.===============#
#==============================================================================#


################################################################################
###########       Part 1.    Predefinition and Invocation           ############
################################################################################

n_data <- length(data.pm.2017$Station_ID) #5244
n_time <- as.integer(n_data/n_stations)#12month
all_ids <- data.coordinates.utm$Station_ID  


# Run the GA to select the optimal combination of 10 station IDs.  
# - n_select = 10: Number of stations to select  
# - popSize = 5: Number of individuals per population in GA  
# - maxiter = 10: Maximum iterations for the GA  
# - run = 5: Early stopping if no improvement for 5 generations  
# The function uses global 'all_ids' and 'fitness_function'.  
# The result is a list with:  
#   best_id: The optimal station ID combination  
#   lowest_rmse: The corresponding lowest RMSE  
#   GA_res: Full GA result object  
result <- run_GA_station_selection(n_select = 10, 
                                   popSize = 5, 
                                   maxiter = 10,
                                   run = 5)
# NOTICE!!
# The final result is the best ID and RMSE. 
# If you need other metrics, input the ID sequence into the GA_INLA function. 
# All result metrics, such as MAE, COR, and R2, are available in validation_val.
# Extract and return them as needed.
# The log file, documenting the optimization process of each generation in the 
# genetic algorithm, is stored in GA_ID_RMSE_log.txt.

################################################################################
###########       Part 2.    GA-INLA Function collection            ############
################################################################################

# Function: run_GA_station_selection  
# Description: Runs the GA optimization and prints/returns best station IDs and lowest RMSE.  
# Uses external global variables: all_ids, fitness_function.  
#  
# Parameters:  
# - n_select: Integer. Number of stations to select.  
# - popSize: Integer. Population size for GA.  
# - maxiter: Integer. Maximum number of iterations for GA.  
# - run: Integer. Early stopping rounds for GA.  
#  
# Returns: A list with best_id, lowest_rmse, and GA_res (GA object).  

run_GA_station_selection <- function(  
    n_select,  
    popSize = 20,  
    maxiter = 30,  
    run = 10  
) {  
  n_stations <- length(all_ids)  
  lower <- rep(1, n_select)  
  upper <- rep(n_stations, n_select)  
  GA_res <- ga(  
    type = "permutation",  
    fitness = fitness_function,  
    min = lower,  
    max = upper,  
    popSize = popSize,  
    maxiter = maxiter,  
    run = run  
  )  
  best_x <- as.integer(round(GA_res@solution[1, 1:n_select]))  
  best_id <- all_ids[best_x]  
  best_rmse <- -GA_res@fitnessValue  
  cat("Best ID combination:", best_id, "\n")  
  cat("Lowest RMSE:", best_rmse, "\n")  
  return(list(best_id = best_id, lowest_rmse = best_rmse, GA_res = GA_res))  
}

# Function: fitness_function  
# Description: Fitness function for the genetic algorithm (GA). Given a candidate solution (vector of indices),  
#              selects corresponding station IDs, evaluates them using GA_Inla (INLA spatial validation),  
#              logs the selected IDs and RMSE to a file, and returns the negative RMSE as the fitness value.  
#  
# Parameters:  
# - x: Vector. Indices representing a candidate solution in GA; used to select n_select station IDs from all_ids.  
#  
# Returns: Numeric. The negative RMSE (so that GA maximizes fitness == minimizes RMSE).  
#  
# Side effects:  
# - Appends the selected IDs and corresponding RMSE to the file "GA_ID_RMSE_log2.txt".  
#  
# Dependencies:  
# - Variables: all_ids, n_select  
# - Functions: GA_Inla 
fitness_function <- function(x) {  
  sel_id <- all_ids[x[1:n_select]]  
  rmse <- tryCatch({  
    GA_Inla(sel_id)
  }, error=function(e) 1e6)
  cat("select ID =", paste(sel_id, collapse = ","), "| RMSE =", rmse, "\n",   
      file = "GA_ID_RMSE_log.txt", append = TRUE)
  -rmse  
}

# Function: GA_Inla  
# Description: Given a vector of station IDs for model estimation, this function partitions data,  
#              performs spatial triangulation, constructs the SPDE model, fits an INLA spatial model,  
#              and computes validation metrics on the hold-out stations. Returns the RMSE on validation data.  
#  
# Parameters:  
# - x: Vector. A vector of station IDs to be used for model estimation (training set).  
#  
# Returns: Numeric. The RMSE (Root Mean Squared Error) of the model predictions on the validation stations.  
#  
# Dependencies:  
# - Uses global variables/data: data.stationID, data.coordinates.utm, data.pm.2017, data.border, n_time  
# - Requires: dplyr, INLA, and custom function calculate_validation_metrics.
GA_Inla <- function(x){
  
  ##PART1 DATA Partition and Standardize
  est_ID <- x 
  val_ID <- setdiff(data.stationID, est_ID)  
  
  est_station <- data.coordinates.utm %>% filter(Station_ID %in% est_ID)  
  val_station <- data.coordinates.utm %>% filter(Station_ID %in% val_ID)  
  est_data <- data.pm.2017 %>% filter(Station_ID %in% est_ID)  
  val_data <- data.pm.2017 %>% filter(Station_ID %in% val_ID)  
  cat("**----GA_Inla DATA Partition Success!---**\n")  
  
  mean_cov <- colMeans(data.pm.2017[,2:3], na.rm=TRUE)  
  sd_cov <- apply(data.pm.2017[,2:3], 2, sd, na.rm=TRUE)  
  names(mean_cov) <- names(sd_cov) <- names(data.pm.2017)[2:3]  
  
  est_data <- est_data %>%  
    mutate(across(all_of(names(mean_cov)),  
                  ~(. - mean_cov[cur_column()]) / sd_cov[cur_column()])) %>%  
    mutate(logPM = log(PM))  
  
  val_data <- val_data %>%  
    mutate(across(all_of(names(mean_cov)),  
                  ~(. - mean_cov[cur_column()]) / sd_cov[cur_column()])) %>%  
    mutate(logPM = log(PM))  
  
  
  ###########PART2 Triangulation using borders
  mesh =
    inla.mesh.2d(loc=cbind(est_station$UTM_X,
                           est_station$UTM_Y),
                 loc.domain=data.border,
                 max.edge=c(15, 100),
                 min.angle=c(26, 21),
                 cutoff=5, 
                 plot.delay=NULL
    )
  # plot(mesh)
  # points(x=data_WH_border$UTM_X,y=data_WH_border$UTM_Y, cex=0.1,col='red')
  # points(x=est_station$UTM_X, y=est_station$UTM_Y, pch =17,cex=1, col="blue")
  # 
  spde <- inla.spde2.pcmatern(
    mesh = mesh, alpha = 3/2, constr = TRUE,
    prior.range = c(25, 0.01), # P(range < 10000) = 0.01
    prior.sigma = c(3, 0.01) # P(sigma > 3) = 0.01
  )
  
  field.indices =
    inla.spde.make.index("field",
                         n.spde=spde$n.spde,
                         n.group=n_time)
  lengths(field.indices)
  print("**----GA_INLA  Triangulation Success!---**")
  ###########      PART3 Make the SPDE object and the formula    #############
  A.est =
    inla.spde.make.A(mesh,
                     loc=
                       as.matrix(data.coordinates.utm[est_data$Station_ID,
                                                      c("UTM_X","UTM_Y")]),
                     group=est_data$month,
                     n.group=n_time
    )
  A.val =
    inla.spde.make.A(mesh,
                     loc=
                       as.matrix(data.coordinates.utm[val_data$Station_ID,
                                                      c("UTM_X","UTM_Y")]),
                     group=val_data$month,
                     n.group=n_time)
  stack.est =
    inla.stack(data=list(logPM=est_data$logPM),
               A=list(A.est, 1),
               effects=
                 list(c(field.indices,
                        list(Intercept=1)),
                      list(est_data[,2:3])),
               tag="est")
  stack.val =
    inla.stack(data=list(logPM=NA),
               A=list(A.val, 1),
               effects=
                 list(c(field.indices,
                        list(Intercept=1)),
                      list(val_data[,2:3])),
               tag="val")
  stack = inla.stack(stack.est, stack.val)
  print("**----PART4 Make the SPDE object and the formula  Success!---**")
  
  ###########      PART4 Call INLA and get results               #############
  rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.9)))
  
  formula <- (logPM ~ -1 + Intercept + UTM_X + UTM_Y +
                f(field, model=spde, group=field.group, 
                  control.group=list(model="ar1", hyper = rprior)))
  result.1 =
    inla(formula,
         data=inla.stack.data(stack, spde=spde),
         family="gaussian",
         control.predictor=list(A=inla.stack.A(stack), compute=TRUE),
         control.compute=list(cpo=FALSE,dic = TRUE,config = TRUE,
                              return.marginals.predictor=TRUE),
         control.inla = list(reordering = "metis",strategy='laplace'),
         keep=FALSE, verbose=TRUE)
  
  print(summary(result.1))
  print("**----cube 6 Success!---**")
  
  #############  cube 7 Extract results                             
  
  validation_val <- calculate_validation_metrics(stack,result.1,val_data)
  return(validation_val$rmse) 
}









