#==============================================================================#
#========================= GA Combined with INLA===============================#
#===== The following is the full code for GA Combined with INLA.===============#
#==============================================================================#

################################################################################
###########      PART 1. Load data and library            ######################
################################################################################
install.packages("GA")
library("GA")
library(INLA)
library(dplyr)  
library(tidyr)  
load("0RData/origin.data.rdata")

################################################################################
###########      PART 2. Predefinition and Invocation     ######################
################################################################################
n_data <- length(data.pm.2017$Station_ID) #5244
all_ids <- data.coordinates.utm$Station_ID  
n_time <- as.integer(n_data/length(all_ids))#12month



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

################################################################################
###########   Part 3.  GA-INLA Function collection           ###################
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