#' configure_obj
#' 
#' @details 
#' given a data frame and some user defined settings will return an TMB object that represents either a GLM, or geo-statistical model.
#' is implied in this model, otherwise you could just use the standard GLM approach
#' @param data data.frame, assumes it has an column name 'area'
#' @param include_epsilon time-invariant spatial GF
#' @param include_omega time-varying spatial GF
#' @param response_variable_label character 
#' @param time_variable_label character assumed to be an integer variable 
#' @param catchability_covariates string vector 
#' @param catchability_covariate_type string vector indicating if the factor, numeric
#' @param spatial_covariates vector of strings
#' @param spatial_covariate_type string vector indicating if the factor, numeric
#' @param spline_catchability_covariates string vector 
#' @param spline_spatial_covariates vector of strings
#' @param projection_df a data.frame needs to have the same variable names (colnames) as data. Should supply variable values for all projection cells over all time steps
#' @param mesh an inla.mesh object that has been created before this function is applied
#' @param family 0 = Poisson, 1 = Negative Binomial, 2 = Gaussian, 3 = Gamma 
#' @param link link function 0 = log, 1 = logit, 2 = probit, 3 = inverse, 4 = identity
#' @param trace_level 'none' don't print any information, 'low' print steps in the function 'medium' print gradients of TMB optimisation, 'high' print parameter candidates as well as gradients during oprimisation. 
#' @export
#' @importFrom sp coordinates
#' @importFrom INLA inla.mesh.projector inla.spde2.matern inla.spde.make.A
#' @importFrom RANN nn2
#' @importFrom TMB MakeADFun
#' @importFrom mgcv gam PredictMat
#' @importFrom Matrix .bdiag
#' @importFrom stats model.matrix rnorm sd terms terms.formula
#' @return: list of estimated objects and data objects
configure_obj = function(data, projection_df, mesh, family, link, include_omega, include_epsilon, response_variable_label, time_variable_label, catchability_covariates = NULL, catchability_covariate_type = NULL, spatial_covariates = NULL, spatial_covariate_type = NULL, spline_catchability_covariates = NULL,
                         spline_spatial_covariates = NULL, trace_level = "none") {
  if(!trace_level %in% c("none", "low", "medium","high"))
    stop(paste0("trace_level needs to be 'none', 'low', 'medium','high'"))
  if(class(data) != "SpatialPointsDataFrame")
    stop(paste0("data needs to be of class SpatialPointsDataFrame, see the example for more information, can be converted by using coordinates(data) <- ~ x + y"))
  if(class(projection_df) != "SpatialPointsDataFrame")
    stop(paste0("projection_df needs to be of class SpatialPointsDataFrame, see the example for more information, can be converted by using coordinates(projection_df) <- ~ x + y"))
  if(class(mesh) != "inla.mesh")
    stop(paste0("mesh needs to be of class inla.mesh, see the example for more information"))
  if(!family %in% c(0:3))
    stop(paste0("family needs to be a value from 0 to 3, for a valid distribution"))
  if(!link %in% c(0:4))
    stop(paste0("link needs to be a value from 0 to 4, for a valid link function"))
  if(length(catchability_covariates) != length(catchability_covariate_type))
    stop(paste0("catchability_covariates needs to be length catchability_covariate_type"))
  if(length(spatial_covariates) != length(spatial_covariate_type))
    stop(paste0("spatial_covariates needs to be length spatial_covariate_type"))
  vars_needed = c(response_variable_label, time_variable_label, catchability_covariates, spatial_covariates, "area")
  proj_vars_needed = c(response_variable_label, time_variable_label, spatial_covariates, "area")
  # get response variable
  if(!all(proj_vars_needed %in% colnames(projection_df@data))) 
    stop(paste0("projection_df: needs colnames ", paste(proj_vars_needed, collapse = ", ")))
  if(!all(vars_needed %in% colnames(data@data)))
    stop(paste0("data: needs colnames ", paste(vars_needed, collapse = ", ")))
  
  if(length(catchability_covariates) != length(catchability_covariate_type)) 
    stop(paste0("catchability_covariates needs to be the same length as catchability_covariate_type"))
  if(length(spatial_covariates) != length(spatial_covariate_type)) 
    stop(paste0("spatial_covariates needs to be the same length as spatial_covariate_type"))
  
  # time step 
  time_variable = data@data[,time_variable_label]
  if(class(time_variable) != "integer") 
    stop(paste0("time_variable_name needs to be an integer this can be achieved by as.integer(data@data$", time_variable_label,")"))
  
  time_levels = sort(unique(time_variable))
  
  if(length(time_levels) < 2) 
    stop(paste0("Need at least two time-steps to run this type of model, found ", length(time_levels)))
    
  if(trace_level != "none")
    print(paste0("Passed initial input checks"))
  
  ## map mesh to observations
  A = inla.spde.make.A(mesh, loc = cbind(coordinates(data[, 1]), coordinates(data[, 2])))
  ## Create Sparse Matern objects to pass to TMB
  spde = inla.spde2.matern(mesh, alpha = 2)
  ## get some index
  n_t = length(unique(time_variable))
  # number of spatial cells in projection matrix
  n_z = length(unique(coordinates(projection_df)[,1])) * length(unique(coordinates(projection_df)[,2]))
  # number of vertices ie. random effects for each GF
  n_v = mesh$n
  n = nrow(data)
  
  ## model matrix for linear effects
  time_model_matrix = evalit(paste0("model.matrix(",response_variable_label," ~ 0 + factor(",time_variable_label,"),  data = data@data)"))
  model_matrix = matrix(1, nrow = n, ncol = 1);
  if(length(catchability_covariates) != 0)
    model_matrix = evalit(paste0("model.matrix(",response_variable_label," ~ ",paste(ifelse(catchability_covariate_type == "factor", paste0("factor(",catchability_covariates,")"), catchability_covariates), collapse = " + "),",  data = data@data)"))
  
  if(trace_level != "none")
    print(paste0("built model_matrix"))
  
  data_spatial_model_matrix = matrix(1, ncol = 2, nrow = n) 
  ## needs to be 2 columns as we estimate n-1 coeffectients so we need at least 2 cols so we have 1 transformed parameter
  ## will be ignored as, the transformed parameters will be fixed at 0
  if(length(spatial_covariates) != 0) 
    data_spatial_model_matrix = evalit(paste0("model.matrix(",response_variable_label," ~ 0 + ",paste(ifelse(spatial_covariate_type == "factor", paste0("factor(",spatial_covariates,")"), spatial_covariates), collapse = " + "),",  data = data@data)"))
  
  ## model matrix for spline covariates
  ## Spline based stuff
  spline_ = spline_spatial_ = NULL
  S_catchability_list = list()
  S_catchability_reporting_list = list()
  S_spatial_list = list()
  S_spatial_reporting_list = list()
  if(length(spline_catchability_covariates) > 0) {
    spline_ = evalit(paste0("mgcv::gam(",response_variable_label," ~ ", paste("s(", spline_catchability_covariates,", bs = 'cs')",collapse = " + "),", data = data@data, fit = F)"))
    if(trace_level == "high")
      print(paste0("length spline = ", length( spline_$smooth)))
    
    for(i in 1:length(spline_$smooth)) {
      S_null = spline_$smooth[[i]]$S[[1]]
      for_plotting = seq(min(data@data[,spline_catchability_covariates[i]]),max(data@data[,spline_catchability_covariates[i]]),by = diff(range(data@data[,spline_catchability_covariates[i]])) / 50)
      forReport = evalit(paste0("mgcv::PredictMat(spline_$smooth[[i]], data = data.frame(",spline_catchability_covariates[i]," = for_plotting))"))
      S_catchability_list[[i]] = S_null
      S_catchability_reporting_list[[i]] = forReport
    }
  } else {
    spline_ = evalit(paste0("mgcv::gam(",response_variable_label," ~ s(area, bs = 'cs'), data = data@data, fit = F)"))
    if(trace_level == "high")
      print(paste0("length spline = ", length( spline_$smooth)))
    S_null = spline_$smooth[[1]]$S[[1]]
    for_plotting = seq(min(data@data[,"area"]),max(data@data[,"area"]),by = diff(range(data@data[,"area"])) / 50)
    forReport = mgcv::PredictMat(spline_$smooth[[1]], data = data.frame(area = for_plotting))
    S_catchability_list[[1]] = S_null
    S_catchability_reporting_list[[1]] = forReport
  }

  if(trace_level != "none")
    print(paste0("Passed catchability spline section"))
  
  if(length(spline_spatial_covariates) > 0) {
    spline_spatial_ = evalit(paste0("mgcv::gam(",response_variable_label," ~ ", paste("s(", spline_spatial_covariates,", bs = 'cs')",collapse = " + "),", data = data@data, fit = F)"))
    if(trace_level == "high")
      print(paste0("length spline_spatial_ = ", length( spline_spatial_$smooth)))
    for(i in 1:length(spline_spatial_$smooth)) {
      S_null = spline_spatial_$smooth[[i]]$S[[1]]
      for_plotting = seq(min(data@data[,spline_spatial_covariates[i]]), max(data@data[,spline_spatial_covariates[i]]), by = diff(range(data@data[,spline_spatial_covariates[i]])) / 50)
      forReport = evalit(paste0("mgcv::PredictMat(spline_spatial_$smooth[[i]], data = data.frame(",spline_spatial_covariates[i]," = for_plotting))"))
      S_spatial_list[[i]] = S_null
      S_spatial_reporting_list[[i]] = forReport
    }
  } else {
    spline_spatial_ = evalit(paste0("mgcv::gam(",response_variable_label," ~ s(area, bs = 'cs'), data = data@data, fit = F)"))
    if(trace_level == "high")
      print(paste0("length spline_spatial_ = ", length( spline_spatial_$smooth)))
    
    S_null = spline_spatial_$smooth[[1]]$S[[1]]
    for_plotting = seq(min(data@data[,"area"]),max(data@data[,"area"]),by = diff(range(data@data[,"area"])) / 50)
    forReport = mgcv::PredictMat(spline_spatial_$smooth[[1]], data = data.frame(area = for_plotting))
    S_spatial_list[[1]] = S_null
    S_spatial_reporting_list[[1]] = forReport
  }
  if(trace_level != "none")
    print(paste0("Passed spatial spline section"))
  
  S_catchability_combined = .bdiag(S_catchability_list)         # join S's in sparse matrix
  S_catchability_dims = unlist(lapply(S_catchability_list, nrow)) # Find dimension of each S
  S_catchability_design_matrix = .bdiag(S_catchability_reporting_list)
  S_spatial_combined = .bdiag(S_spatial_list)         # join S's in sparse matrix
  S_spatial_dims = unlist(lapply(S_spatial_list, nrow)) # Find dimension of each S
  S_spatial_design_matrix = .bdiag(S_spatial_reporting_list)
  
  
  if(trace_level != "none")
    print(paste0("Passed: model matrix configurations and spline configurations"))
  
  ## Projection model matrix stuff
  X_spatial_proj_zpt = array(1, dim = c(n_z, ncol(data_spatial_model_matrix), n_t))
  X_spatial_ipt = array(1, dim = c(n, ncol(data_spatial_model_matrix), n_t))
  spline_spatial_model_matrix_proj_zpt = array(0, dim = c(n_z, ncol(spline_spatial_$X[,-1]), n_t))
  spline_spatial_model_matrix_ipt = array(spline_spatial_$X[,-1], dim = c(n, ncol(spline_spatial_$X[,-1]), n_t))
  
  
  ## map mesh to extrapolation grid. Assumes projection_df has the same spatial cell for all time-cells
  ## so need to define the projector just for a single time-step
  proj_df_t = subset(projection_df, subset = projection_df@data[,time_variable_label] == time_levels[1])
  Proj <- inla.mesh.projector(mesh, loc = coordinates(proj_df_t))
  P = Proj$proj$A
  Proj_area = proj_df_t@data$area
  
  if(trace_level != "none")
    print(paste0("Passed: mapping projection grid to mesh"))
  
  for(t in 1:n_t) {
    proj_df_subset = subset(projection_df, subset = projection_df@data[,time_variable_label] == time_levels[t])

    if(length(spatial_covariates) > 0) {
      proj_spatial_model_matrix = evalit(paste0("model.matrix(formula = ",response_variable_label," ~ 0 + ", paste(ifelse(spatial_covariate_type == "factor", paste0("factor(",spatial_covariates,")"), spatial_covariates), collapse = " + "),",  data = proj_df_subset@data)"))
      data_spatial_model_matrix = evalit(paste0("model.matrix(formula = ",response_variable_label," ~ 0 + ",paste(ifelse(spatial_covariate_type == "factor", paste0("factor(",spatial_covariates,")"), spatial_covariates), collapse = " + "),",  data = data@data)"))
      X_spatial_proj_zpt[,,t] = proj_spatial_model_matrix
      X_spatial_ipt[,,t] = data_spatial_model_matrix
    }
    if(length(spline_spatial_covariates) > 0) {
      ## for the spatial spline
      S_spatial_proj_model_mat = NULL
      for(i in 1:length(spline_spatial_covariates)) {
        spatial_spline_proj = PredictMat(spline_spatial_$smooth[[i]], data = proj_df_subset@data)
        S_spatial_proj_model_mat = cbind(S_spatial_proj_model_mat, spatial_spline_proj)
      }
      spline_spatial_model_matrix_proj_zpt[,,t] = S_spatial_proj_model_mat
    }
  }
  if(trace_level != "none")
    print(paste0("Passed: projection model matrix construction"))
  
  ## Set up TMB object.
  tmb_data <- list(
               model = "SpatialTemporalCPUE",
               n_i = n,
               n_t = n_t,
               y_i = data@data[,response_variable_label],
               t_i = time_variable - min(time_variable), # index for C++ language could be 1990 1991 etc,
               area_i = data@data$area,
               obs_t = as.numeric(table(time_variable)),
               A = A,
               spde = spde$param.inla[c("M0","M1","M2")],
               Proj = P,
               Proj_Area = Proj_area,
               family = family,
               link = link,
               pref_coef_bounds = c(-0.001, 5),
               model_matrix = model_matrix,
               time_model_matrix = time_model_matrix,
               X_spatial_ipt = X_spatial_ipt,
               X_spatial_proj_zpt = X_spatial_proj_zpt,
               apply_pref = 0,
               omega_indicator = ifelse(include_omega, 1, 0),
               epsilon_indicator = ifelse(include_epsilon, 1, 0),
               spline_flag = c(ifelse(length(spline_catchability_covariates) > 0, 1, 0), ifelse(length(spline_spatial_covariates) > 0, 1, 0)),
               spline_model_matrix = spline_$X[,-1],
               spline_spatial_model_matrix_ipt = spline_spatial_model_matrix_ipt,
               spline_spatial_model_matrix_proj_zpt = spline_spatial_model_matrix_proj_zpt,
               S = S_catchability_combined,
               Sdims = S_catchability_dims,
               designMatrixForReport = S_catchability_design_matrix,
               S_spatial = S_spatial_combined,
               Sdims_spatial = S_spatial_dims,
               designMatrixForReport_spatial = S_spatial_design_matrix
  )
  year_ndx_for_each_obs = matrix(-99, nrow = tmb_data$n_t, ncol = max(tmb_data$obs_t))
  for(t in 1:tmb_data$n_t) {
    year_ndx_for_each_obs[t, 1:tmb_data$obs_t[t]] = which(tmb_data$t_i == (t - 1)) - 1
  }
  tmb_data$year_ndx_for_each_obs = year_ndx_for_each_obs
  
  if(trace_level != "none")
    print(paste0("Passed: TMB data list construction"))
  
  # parameters for TMB
  dis_cor_0.1 = mean(0.5 * c(diff(range(mesh$loc[,1])), diff(range(mesh$loc[,2])))) # halfway 
  sigma_guess = 1 ## more intuative to think of marginal variance when setting starting values for params

  params = list(
    betas = rep(0, ncol(tmb_data$model_matrix)),
    constrained_spatial_betas = rep(0, max(dim(tmb_data$X_spatial_ipt)[2] - 1,1)),
    constrained_time_betas = rep(0, max(dim(tmb_data$time_model_matrix)[2] - 1,1)),
    ln_kappa_epsilon = log(sqrt(8) / dis_cor_0.1),#sqrt(8) / dis_cor_0.1,
    ln_kappa_omega = log(sqrt(8) / dis_cor_0.1),#sqrt(8) / dis_cor_0.1,
    ln_tau_epsilon =  log(1/(sigma_guess * sqrt(4*pi * sqrt(8) / dis_cor_0.1^2))),
    ln_tau_omega =  log(1/(sigma_guess * sqrt(4*pi * sqrt(8) / dis_cor_0.1^2))),
    logit_pref_coef = logit_general(0, tmb_data$pref_coef_bounds[1], tmb_data$pref_coef_bounds[2]),
    logit_eps_rho = logit_general(0, -0.99, 0.99),
    ln_phi = 0,
    omega_input = rep(0,spde$n.spde),
    epsilon_input = array(0,dim = c(spde$n.spde, tmb_data$n_t)),
    gammas = rep(0,sum(tmb_data$Sdims)),  # Spline coefficients
    ln_lambda = rep(0,length(tmb_data$Sdims)), #Log spline penalization coefficients
    gammas_spatial = rep(0,sum(tmb_data$Sdims_spatial)),  # Spline coefficients
    ln_lambda_spatial = rep(0,length(tmb_data$Sdims_spatial)) #Log spline penalization coefficients
  )
  
  if(trace_level != "none")
    print(paste0("Passed: TMB parameter list construction"))
  
  any_errors = check_inputs(tmb_data = tmb_data, tmb_params = params)
  if(!any_errors$result) {
    print("Error found incompatible data and parameter lists. I am returning the errors to help debug this function.")
    return(any_errors$errors)
  }
  ## set which parameters are being estimated and which are held constant.drr
  pars_to_fix = c()
  if(!include_epsilon)
    pars_to_fix = c(pars_to_fix, "epsilon_input")
  if(!include_omega)
    pars_to_fix = c(pars_to_fix, "omega_input")
  if(length(spatial_covariates) == 0)
    pars_to_fix = c(pars_to_fix, "constrained_spatial_betas")
  if(length(spline_catchability_covariates) == 0)
    pars_to_fix = c(pars_to_fix, "ln_lambda", "gammas")
  if(length(spline_spatial_covariates) == 0)
    pars_to_fix = c(pars_to_fix, "ln_lambda_spatial", "gammas_spatial")
  ## do we need to estimate a dispersion parameter
  if(family %in% c(0))
    pars_to_fix = c(pars_to_fix, "ln_phi")

  fixed_pars = list()
  if(length(pars_to_fix) > 0)
    fixed_pars = fix_pars(par_list = params, pars_to_exclude = pars_to_fix)

  ## compile model
  #compile(file.path("src","debug_standalone_version.cpp"))
  #file.exists(file.path("src","debug_standalone_version.dll"))
  #dyn.load(dynlib(file.path("src","debug_standalone_version")))
  obj = NULL
  ## create obj
  if(include_epsilon & include_omega) {
    obj <- MakeADFun(tmb_data, params, random = c("epsilon_input","omega_input"), map = fixed_pars, DLL = "CPUEspatial_TMBExports", method = "nlminb", hessian = T, silent=T)
  } else if(include_epsilon & !include_omega) {
    obj <- MakeADFun(tmb_data, params, random = c("epsilon_input"), map = fixed_pars, DLL = "CPUEspatial_TMBExports", method = "nlminb", hessian = T, silent=T)
  } else if(!include_epsilon & include_omega) {
    obj <- MakeADFun(tmb_data, params, random = c("omega_input"), map = fixed_pars, DLL = "CPUEspatial_TMBExports", method = "nlminb", hessian = T, silent=T)
  }
  
  if(trace_level != "none")
    print(paste0("Passed: successfully built obj"))
  
  return(list(obj = obj, tmb_pars = params, tmb_data = tmb_data))
}


