#' predict_out_of_sample
#' 
#' @details 
#' given an independent data set predict the model at these locations and time points. Can think of this is as the predict.glm function with newdata specified
#' 
#' A note on factors when predicting. if out of sample data set is doesn't have all the levels of a factor that are estimated
#' in the observed_df use the levels() function to set the levels of a factor so there are the same levels identified as in the observed_df. This is because
#' we construct a model matrix with coeffecients based on levels for a factor. So to correctly predict with the estimated coeffecients we would need the same structured model matricies
#' 
#' @param conf_obj CPUEspatial object, created from configure_obj and the obj component has been optimised
#' @param out_of_sample_df SpatialPointsDataFrame data set which should be consistent with observed_df used in the configure_obj() function in terms of variable names and classes
#' @param observed_df SpatialPointsDataFrame, which contains response variable and covariates for glmm analysis mut contain column 'area'
#' @param simulate whether to simulate parameter uncertainty. Wont represent prediction intervals as we have not incorporated variation in the response variable
#' @param mesh an inla.mesh object that has been created before this function is applied
#' @param n_sims number of simulated parameters for predictions
#' @return list model predictions (fitted) Squared prediction errro (SPE)
#' @importFrom INLA inla.spde.make.A
#' @importFrom stats formula
#' @importFrom TMB sdreport
#' @importFrom MASS mvrnorm
#' @export
predict_out_of_sample = function(conf_obj, out_of_sample_df, observed_df, mesh, simulate = FALSE, n_sims = 100) {
  if(class(out_of_sample_df) != "SpatialPointsDataFrame")
    stop(paste0("out_of_sample_df needs to be of class SpatialPointsDataFrame, see the example for more information, can be converted by using coordinates(out_of_sample_df) <- ~ x + y"))
  response_variable_label = eval(conf_obj$Call$func_call$response_variable_label)
  time_variable_label = eval(conf_obj$Call$func_call$time_variable_label)
  catchability_covariates = eval(conf_obj$Call$func_call$catchability_covariates)
  spatial_covariates = eval( conf_obj$Call$func_call$spatial_covariates)
  spline_spatial_covariates = eval( conf_obj$Call$func_call$spline_spatial_covariates)
  spline_catchability_covariates = eval( conf_obj$Call$func_call$spline_catchability_covariates)
  linear_basis = eval( conf_obj$Call$func_call$linear_basis)
  linear_basis = ifelse(is.null(linear_basis), 0, linear_basis)
  n_pred = nrow(out_of_sample_df@data)
  vars_needed = c(response_variable_label, time_variable_label, catchability_covariates, spatial_covariates, spline_spatial_covariates, spline_catchability_covariates, "area")
  if(!all(vars_needed %in% colnames(out_of_sample_df@data)))
    stop(paste0("out_of_sample_df: needs colnames ", paste(vars_needed, collapse = ", ")))
  ## original time-steps
  original_time_steps = unique(observed_df@data[, time_variable_label])
  n_t = conf_obj$tmb_data$n_t
  ## build model matricies
  year_factor = factor(out_of_sample_df@data[, time_variable_label], levels = levels(factor(observed_df@data[,time_variable_label])))
  y_pred_raw = out_of_sample_df@data[, response_variable_label]
  time_model_matrix <- model.matrix(formula("y_pred_raw ~ 0 + year_factor"), drop.unused.levels = F)
  check_col = check_dims(time_model_matrix, conf_obj$tmb_data$time_model_matrix, row = FALSE)
  if(!check_col$passed_check)
    stop(paste0("When constructing time_model_matrix, ", check_col$msg))
  
  model_matrix = matrix(1, nrow = n_pred, ncol = 1);
  if(length(catchability_covariates) != 0) {
    ff = paste0(response_variable_label," ~ ",paste(catchability_covariates, collapse = " + "))
    m <- model.frame(formula(conf_obj$Call$catchability), out_of_sample_df@data, drop.unused.levels = F)
    model_matrix <- model.matrix(formula(ff), m)
  }
  check_col = check_dims(model_matrix, conf_obj$tmb_data$model_matrix, row = FALSE)
  if(!check_col$passed_check)
    stop(paste0("When constructing model_matrix, ", check_col$msg))
  
  data_spatial_model_matrix = matrix(1, ncol = 2, nrow = n_pred) 
  if(length(spatial_covariates) != 0) {
    m <- model.frame(formula(conf_obj$Call$spatial), out_of_sample_df@data, drop.unused.levels = F)
    data_spatial_model_matrix <- model.matrix(formula(conf_obj$Call$spatial), m)
  }
  
  check_col = check_dims(data_spatial_model_matrix, conf_obj$tmb_data$X_spatial_ipt[,,1], row = FALSE)
  if(!check_col$passed_check)
    stop(paste0("When constructing spatial_model_matrix, ", check_col$msg))
  
  ## model matrix for spline covariates
  ## Spline based stuff
  spline_ <- NULL
  spline_spatial_ <- NULL
  for_plotting <- NULL
  S_catchability_list <- list()
  S_spatial_list <- list()
  if(length(spline_catchability_covariates) > 0) {
    ff = formula(paste0(response_variable_label," ~ ", paste("s(", spline_catchability_covariates,", bs = 'cs')",collapse = " + ")))
    spline_ <- mgcv::gam(ff, data = out_of_sample_df@data, fit = F)
    
    for(i in 1:length(spline_$smooth)) {
      S_null <- spline_$smooth[[i]]$S[[1]]
      S_catchability_list[[i]] <- S_null
    }
  } else {
    ## create a dummy variable
    ff = formula(paste0(response_variable_label," ~ s(area, bs = 'cs')"))
    spline_ <- mgcv::gam(ff, data = out_of_sample_df@data, fit = F)
    S_null <- spline_$smooth[[1]]$S[[1]]
    S_catchability_list[[1]] <- S_null
  }

  if(length(spline_spatial_covariates) > 0) {
    ff = formula(paste0(response_variable_label," ~ ", paste("s(", spline_spatial_covariates,", bs = 'cs')",collapse = " + ")))
    spline_spatial_ <- mgcv::gam(ff, data = out_of_sample_df@data, fit = F)
    for(i in 1:length(spline_spatial_$smooth)) {
      S_null <- spline_spatial_$smooth[[i]]$S[[1]]
      S_spatial_list[[i]] <- S_null
    }
  } else {
    ## create a dummy variable
    ff = formula(paste0(response_variable_label," ~ s(area, bs = 'cs')"))
    spline_spatial_ <- mgcv::gam(ff, data = observed_df@data, fit = F)
    S_null <- spline_spatial_$smooth[[1]]$S[[1]]
    S_spatial_list[[1]] <- S_null
  }

  S_catchability_combined <- .bdiag(S_catchability_list)         # join S's in sparse matrix
  S_catchability_dims <- unlist(lapply(S_catchability_list, nrow)) # Find dimension of each S
  S_spatial_combined <- .bdiag(S_spatial_list)         # join S's in sparse matrix
  S_spatial_dims <- unlist(lapply(S_spatial_list, nrow)) # Find dimension of each S
  
  ## build mesh interpolation
  A = inla.spde.make.A(mesh, loc = cbind(coordinates(out_of_sample_df[, 1]), coordinates(out_of_sample_df[, 2])))
  data_vertex_ndx = nn2(mesh$loc[,c(1,2)], coordinates(out_of_sample_df), k = 1)$nn.idx
  ## deal with epsilon 
  time_variable = out_of_sample_df@data[, time_variable_label]
  time_variable = as.integer(as.character(time_variable))
  obs_t = table(match(x = time_variable, table = original_time_steps))
  t_i = time_variable - min(original_time_steps)
  obs_lab = as.integer(names(obs_t))
  year_ndx_for_each_obs = matrix(-99, nrow = n_t, ncol = max(obs_t))
  obs_ndx = 1
  for(t in 1:n_t) {
    if(t %in% obs_lab) {
      year_ndx_for_each_obs[t, 1:obs_t[obs_ndx]] = which(t_i == (t - 1))
      obs_ndx = obs_ndx + 1
    }
  }
  
  ## Calculate fitted value
  ran_effects = FALSE
  if(length(conf_obj$obj$env$random) > 0)
    ran_effects = TRUE
  
  MLE_pars = conf_obj$obj$env$last.par.best
  fitted = NULL
  obs = out_of_sample_df@data[, response_variable_label]
  SPE = NULL
  if(simulate) {
    sd_rep = sdreport(conf_obj$obj, getJointPrecision = T)
    if(ran_effects) {
      sim_pars = t(rmvnorm_prec(MLE_pars, sd_rep$jointPrecision, n_sims)) ## redraw parameters based on hessian
    } else {
      sim_pars = mvrnorm(n_sims, MLE_pars, sd_rep$cov.fixed) ## redraw parameters based on hessian
    }
    for(j in 1:n_sims) {
      this_rep = conf_obj$obj$report(sim_pars[j, ])
      this_eta = time_model_matrix %*%  this_rep$time_betas +  data_spatial_model_matrix %*% this_rep$spatial_betas + model_matrix %*% this_rep$betas
      if(!is.null(spline_spatial_covariates)) 
        this_eta = this_eta + spline_spatial_$X[,-1] %*% this_rep$gammas_spatial 
      if(!is.null(spline_catchability_covariates))
        this_eta = this_eta + spline_$X[,-1] %*% this_rep$gammas 
      if(conf_obj$tmb_data$omega_indicator) {
        if(linear_basis == 0) {
          this_eta = this_eta + A %*% this_rep$omega_input / this_rep$tau_omega
        } else {
          this_eta = this_eta + this_rep$omega_input[data_vertex_ndx]
        }
      }
      if(conf_obj$tmb_data$epsilon_indicator) {
        for(t in 1:ncol(this_rep$epsilon_input)) {
          this_eps = NULL
          if(linear_basis == 0) {
            this_eps = A %*% this_rep$epsilon_input[,t] / this_rep$tau_epsilon
          } else {
            this_eps = this_rep$epsilon_input[data_vertex_ndx, t]
          }
          obs_this_year = sum(year_ndx_for_each_obs[t,] != -99)
          if(obs_this_year != 0) {
            for(i in 1:obs_this_year)
              this_eta[year_ndx_for_each_obs[t,i]] = this_eta[year_ndx_for_each_obs[t,i]] + this_eps[year_ndx_for_each_obs[t,i]]
          }
        }
      }
      fitted = cbind(fitted,  inverse_link(this_eta, conf_obj$Call$link))
    }
    
    SPE =  apply(fitted, MARGIN = 2, FUN = function(fit, obs) {(obs - fit)^2}, obs = obs)
    
  } else {
    this_rep = conf_obj$obj$report(MLE_pars)
    this_eta = time_model_matrix %*%  this_rep$time_betas +  data_spatial_model_matrix %*% this_rep$spatial_betas + model_matrix %*% this_rep$betas
    if(!is.null(spline_spatial_covariates)) 
      this_eta = this_eta + spline_spatial_$X[,-1] %*% this_rep$gammas_spatial 
    if(!is.null(spline_catchability_covariates))
      this_eta = this_eta + spline_$X[,-1] %*% this_rep$gammas 
    if(conf_obj$tmb_data$omega_indicator) {
      if(linear_basis == 0) {
        this_eta = this_eta + A %*% this_rep$omega_input / this_rep$tau_omega
      } else {
        this_eta = this_eta + this_rep$omega_input[data_vertex_ndx]
      }
    }
    if(conf_obj$tmb_data$epsilon_indicator) {
      for(t in 1:ncol(this_rep$epsilon_input)) {
        this_eps = NULL
        if(linear_basis == 0) {
          this_eps = A %*% this_rep$epsilon_input[,t] / this_rep$tau_epsilon
        } else {
          this_eps = this_rep$epsilon_input[data_vertex_ndx, t]
        }
        obs_this_year = sum(year_ndx_for_each_obs[t,] != -99)
        if(obs_this_year != 0) {
          for(i in 1:obs_this_year)
            this_eta[year_ndx_for_each_obs[t,i]] = this_eta[year_ndx_for_each_obs[t,i]] + this_eps[year_ndx_for_each_obs[t,i]]
        }
      }
    }
    ## invert to go to Y space
    fitted = inverse_link(this_eta, conf_obj$Call$link)
    if(conf_obj$Call$family != 1)
      fitted = out_of_sample_df$area * inverse_link(this_eta, conf_obj$Call$link)
    ## do some metrics Squared prediction error
    SPE = (obs - fitted)^2
  }
  return(list(fitted = fitted, obs = obs, SPE = SPE))
}

