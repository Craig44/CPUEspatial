#' configure_obj
#' 
#' @details 
#' given a data frame and some user defined settings will return an TMB object that represents either a GLM, or geo-statistical model.
#' is implied in this model, otherwise you could just use the standard GLM approach. 
#' pref_hyper_distribution integer specifies whether preference coeffecient (the logit transformed parameter) is time-varying (and thus treated as random effect) if time-varying pref, logit_pref ~ N(mu_pref, cv_pref), 
#'  == 0: Not time-varying (estimated)
#'  == < 0: Not time-varying (Fixed). Used for profiling, diagnosing identfiability
#'  == 1: mu_pref (est) & sd_pref (est), 
#'  == 2: mu_pref (fixed)  & sd_pref(fixed)
#'  == 3: mu_pref (est)  & sd_pref(fixed)
#'  == 4: mu_pref (fixed)  & sd_pref(est)

#' @param observed_df SpatialPointsDataFrame, which contains response variable and covariates for glmm analysis mut contain column 'area'
#' @param projection_df SpatialPointsDataFrameneeds to have the same variable names (colnames) as observed_df. Should supply variable values for all projection cells over all time steps
#' @param include_epsilon boolean time-varying spatial GF
#' @param include_omega boolean time-invariant spatial GF 
#' @param epsilon_structure character either "iid" or "ar1"
#' @param response_variable_label character relating to a column name in observed_df
#' @param time_variable_label character relating to a column names needs to be an integer variable 
#' @param catchability_covariates string vector of column labels corresponding for catchability realted covariates, if type not numeric it is treated as categorical
#' @param spatial_covariates string vector of column labels corresponding for spatial habitat realted covariates, if type not numeric it is treated as categorical
#' @param spline_catchability_covariates string vector 
#' @param spline_spatial_covariates vector of strings indicating column name
#' @param mesh an inla.mesh object that has been created before this function is applied
#' @param family 0 = Gaussian, 1 = Binomial, 2 = Gamma, 3 = Poisson, 4 = Negative Binomial
#' @param link link function 0 = log, 1 = logit, 2 = probit, 3 = inverse, 4 = identity, 5 = inverse squared (1/x^2)
#' @param linear_basis 0 = apply triangulation sparse matrix approach, 1 = Nearest Neighbour
#' @param apply_preferential_sampling whether to jointly model observation location 
#' @param preference_model_type integer 0 = Dinsdale approach, 1 = LGCP lattice approach
#' @param pref_hyper_distribution integer #specifies whether preference coeffecient (the logit transformed parameter) is time-varying (and thus treated as random effect) if time-varying pref, logit_pref = N(mu_pref, sd_pref). See details for more information
#' @param logit_pref_hyper_prior_vals vector<double> specifing the mean and sd (note sd is estimated interanlly as ln_sd to constrain sd > 0, so this value is logged in the model). if estimated as specified by pref_hyper_distribution, these are the starting values, otherwise they are the values fixed during estimation
#' @param trace_level 'none' don't print any information, 'low' print steps in the function 'medium' print gradients of TMB optimisation, 'high' print parameter candidates as well as gradients during oprimisation. 
#' @param projection_raster_layer a RasterLayer object only required if apply_preferential_sampling = TRUE, and preference_model_type == 1. Should be the same resolution as projection_df. Used to collate sample locations. The observed_df slot should have vaules 0 = cell not in projection grid or 1 = active projection cell
#' @param pref_bounds lower and upper bounds
#' @param kappa_bounds lower and upper bounds for kappa parameter
#' @param init_vals list of of parameter values passed to tmb::MakeADFun, useful for investigating different starting values.
#' TODO add whether we want to use Nearest Neighbour approach NN.
#' @export
#' @importFrom sp coordinates
#' @importFrom INLA inla.mesh.projector inla.spde2.matern inla.spde.make.A
#' @importFrom RANN nn2
#' @importFrom TMB MakeADFun
#' @importFrom mgcv gam PredictMat
#' @importFrom Matrix .bdiag
#' @importFrom raster raster rasterize crs
#' @importFrom stats model.matrix rnorm sd terms terms.formula dist
#' @return: list of estimated objects and data objects
configure_obj = function(observed_df, projection_df, mesh, family, link, include_omega, include_epsilon, epsilon_structure = "iid", response_variable_label, time_variable_label, catchability_covariates = NULL, spatial_covariates = NULL, spline_catchability_covariates = NULL,
                         spline_spatial_covariates = NULL, linear_basis = 0, apply_preferential_sampling = FALSE, preference_model_type = 1, pref_hyper_distribution = 0, logit_pref_hyper_prior_vals = c(0,1), projection_raster_layer = NULL, trace_level = "none", pref_bounds= c(-10, 10), kappa_bounds = NULL, init_vals = NULL) {
  Call = list()
  Call$func_call <- match.call()
  if(!trace_level %in% c("none", "low", "medium","high"))
    stop(paste0("trace_level needs to be 'none', 'low', 'medium','high'"))
  if(class(observed_df) != "SpatialPointsDataFrame")
    stop(paste0("observed_df needs to be of class SpatialPointsDataFrame, see the example for more information, can be converted by using coordinates(observed_df) <- ~ x + y"))
  if(class(projection_df) != "SpatialPointsDataFrame")
    stop(paste0("projection_df needs to be of class SpatialPointsDataFrame, see the example for more information, can be converted by using coordinates(projection_df) <- ~ x + y"))
  if(class(mesh) != "inla.mesh")
    stop(paste0("mesh needs to be of class inla.mesh, see the example for more information"))
  if(!family %in% c(0:4))
    stop(paste0("family needs to be a value from 0 to 3, for a valid distribution"))
  if(!link %in% c(0:4))
    stop(paste0("link needs to be a value from 0 to 4, for a valid link function"))
  if(!epsilon_structure %in% c("iid", "ar1"))
    stop("epsilon_structure, needs to either be 'iid' or 'ar1'")
  set_up_dummy_proj = TRUE
  if(apply_preferential_sampling & preference_model_type == 1) {
    if(is.null(projection_raster_layer))
      stop(paste0("If applying preferential sampling with LGCP approach you need to supply projection_raster_layer"))
    set_up_dummy_proj = FALSE
  }
  
  if(family == 1 & link != 1)
    stop("Currently binomial family is only parameterised for link = logit.")
  # check projections are the same.
  #if(is.na(attr(crs(projection_df), "projargs") != attr(crs(observed_df), "projargs")) attr(crs(projection_df), "projargs") != attr(crs(observed_df), "projargs"))
  #  stop(paste0("projection_df and observed_df need to have the same crs projection system, please check raster::crs() for both these objects"))
  
  catchability_covariate_type = NULL
  if(length(catchability_covariates) > 0) {
    for(i in 1:length(catchability_covariates)) {
      this_type = class(observed_df@data[, catchability_covariates[i]])
      catchability_covariate_type = c(catchability_covariate_type, this_type)
    }
  }
  
  spatial_covariate_type = NULL
  if(length(spatial_covariates) > 0) {
    for(i in 1:length(spatial_covariates)) {
      this_type = class(observed_df@data[, spatial_covariates[i]])
      if(this_type == "integer")
        stop(paste0("currently the variable ", spatial_covariates[i], " is of type integer. Change this to either a factor or numeric so this function can correctly configure the model"))
      spatial_covariate_type = c(spatial_covariate_type, this_type)
    }
  }
  vars_needed = c(response_variable_label, time_variable_label, catchability_covariates, spatial_covariates, spline_spatial_covariates, "area")
  proj_vars_needed = c(response_variable_label, time_variable_label, spatial_covariates, "area")
  # get response variable
  if(!all(proj_vars_needed %in% colnames(projection_df@data))) 
    stop(paste0("projection_df: needs colnames ", paste(proj_vars_needed, collapse = ", ")))
  if(!all(vars_needed %in% colnames(observed_df@data)))
    stop(paste0("observed_df: needs colnames ", paste(vars_needed, collapse = ", ")))
  
  if(length(catchability_covariates) != length(catchability_covariate_type)) 
    stop(paste0("catchability_covariates needs to be the same length as catchability_covariate_type"))
  if(length(spatial_covariates) != length(spatial_covariate_type)) 
    stop(paste0("spatial_covariates needs to be the same length as spatial_covariate_type"))
  
  # time step 
  time_variable = observed_df@data[,time_variable_label]
  if(class(time_variable) != "integer") 
    stop(paste0("time_variable_name needs to be an integer this can be achieved by as.integer(observed_df@data$", time_variable_label,")"))
    ## get some index
  time_levels = min(time_variable):max(time_variable) 
  n_t = length(time_levels)
  
  if(length(time_levels) < 2) 
    stop(paste0("Need at least two time-steps to run this type of model, found ", length(time_levels)))
  
  if(trace_level != "none")
    print(paste0("Passed initial input checks"))
  
  ## map mesh to observations
  A = inla.spde.make.A(mesh, loc = cbind(coordinates(observed_df[, 1]), coordinates(observed_df[, 2])))
  ## Create Sparse Matern objects to pass to TMB
  spde = inla.spde2.matern(mesh, alpha = 2)

  
  # number of vertices ie. random effects for each GF
  n_v = mesh$n
  n = nrow(observed_df)
  
  ## model matrix for linear effects
  ff = paste0(response_variable_label," ~ 0 + factor(",time_variable_label,")")
  m <- model.frame(formula(ff), observed_df@data, drop.unused.levels = T)
  time_model_matrix <- model.matrix(formula(ff), data = m)
  Call$time = ff
  if(trace_level != "none")
    print(paste0("built time_model_matrix"))
  
  model_matrix = matrix(1, nrow = n, ncol = 1);
  if(length(catchability_covariates) != 0) {
    ff = paste0(response_variable_label," ~ ",paste(catchability_covariates, collapse = " + "))
    Call$catchability = ff
    m <- model.frame(formula(ff), observed_df@data, drop.unused.levels = T)
    model_matrix <- model.matrix(formula(ff), m)
  }
  
  if(trace_level != "none")
    print(paste0("built model_matrix"))
  data_spatial_model_matrix = matrix(1, ncol = 2, nrow = n) 
  spatial_constrained_coeff_ndx = matrix(0, nrow = 1, ncol = 1);
  spatial_covar_type = c(1) # defualt to numeric
  ## needs to be 2 columns as we estimate n-1 coeffectients so we need at least 2 cols so we have 1 transformed parameter
  ## will be ignored as, the transformed parameters will be fixed at 0
  if(length(spatial_covariates) != 0) {
    ff = paste0(response_variable_label," ~ 0 +",paste(spatial_covariates, collapse = " + "))
    Call$spatial = ff
    m <- model.frame(formula(ff), observed_df@data, drop.unused.levels = T)
    data_spatial_model_matrix <- model.matrix(formula(ff), m)
    p_s = max(attributes(data_spatial_model_matrix)$assign) # should be the same length as spatial covariate
    if(attributes(data_spatial_model_matrix)$dim[2] == 1) {
      spatial_constrained_coeff_ndx = matrix(-99, nrow = p_s, ncol = 1)
    } else {
      spatial_constrained_coeff_ndx = matrix(-99, nrow = p_s, ncol = max(table(attributes(data_spatial_model_matrix)$assign)) - 1)
    }
    counter = 0;
    spatial_covar_type = vector();
    for(i in 1:p_s) {
      if(spatial_covariate_type[i] == "numeric") {
        spatial_covar_type[i] = 1
        spatial_constrained_coeff_ndx[i, 1] = counter
        counter = counter + 1
      } else {
        spatial_covar_type[i] = 0
        for(j in 1:(sum(attributes(data_spatial_model_matrix)$assign == i) - 1)) {
          spatial_constrained_coeff_ndx[i, j] = counter
          counter = counter + 1
        }
      }
    }
  }
  if(trace_level != "none")
    print(paste0("Configured spatial coeffecients"))
  ## model matrix for spline covariates
  ## Spline based stuff
  spline_ <- NULL
  spline_spatial_ <- NULL
  for_plotting <- NULL
  S_catchability_list <- list()
  S_catchability_reporting_list <- list()
  S_spatial_list <- list()
  S_spatial_reporting_list <- list()
  if(length(spline_catchability_covariates) > 0) {
    ff = formula(paste0(response_variable_label," ~ ", paste("s(", spline_catchability_covariates,", bs = 'ts')",collapse = " + ")))
    spline_ <- mgcv::gam(ff, data = observed_df@data, fit = F)
    #attach(spline_)
    if(trace_level == "high") {
      print(paste0("length spline = ", length( spline_$smooth)))
      print(paste0("formula = mgcv::gam(",response_variable_label," ~ ", paste("s(", spline_catchability_covariates,", bs = 'ts')",collapse = " + ")))
    }
    
    for(i in 1:length(spline_$smooth)) {
      S_null <- spline_$smooth[[i]]$S[[1]]
      for_plotting <- data.frame("temp_name" = seq(min(observed_df@data[,spline_catchability_covariates[i]]),max(observed_df@data[,spline_catchability_covariates[i]]),by = diff(range(observed_df@data[,spline_catchability_covariates[i]])) / 50));
      colnames(for_plotting) = spline_catchability_covariates[i]
      forReport <- mgcv::PredictMat(spline_$smooth[[i]], data = for_plotting)
      S_catchability_list[[i]] <- S_null
      S_catchability_reporting_list[[i]] <- forReport
    }
  } else {
    ## create a dummy variable
    ff = formula(paste0(response_variable_label," ~ s(",response_variable_label,", bs = 'ts')"))
    spline_ <-  mgcv::gam(ff, data = observed_df@data, fit = F)
    if(trace_level == "high")
      print(paste0("length spline = ", length( spline_$smooth)))
    S_null <- spline_$smooth[[1]]$S[[1]]
    for_plotting <- seq(min(observed_df@data[,response_variable_label]),max(observed_df@data[,response_variable_label]),by = diff(range(observed_df@data[,response_variable_label])) / 50)
    plot_data = data.frame(x = for_plotting)
    colnames(plot_data) = response_variable_label
    forReport <- mgcv::PredictMat(spline_$smooth[[1]], data = plot_data)
    S_catchability_list[[1]] <- S_null
    S_catchability_reporting_list[[1]] <- forReport
  }
  
  if(trace_level != "none")
    print(paste0("Passed catchability spline section"))
  
  if(length(spline_spatial_covariates) > 0) {
    ff = formula(paste0(response_variable_label," ~ ", paste("s(", spline_spatial_covariates,", bs = 'ts')",collapse = " + ")))
    spline_spatial_ <- mgcv::gam(ff, data = observed_df@data, fit = F)
    if(trace_level == "high") {
      print(paste0("formula = mgcv::gam(",response_variable_label," ~ ", paste("s(", spline_spatial_covariates,", bs = 'ts')",collapse = " + ")))
      print(paste0("length spline_spatial_ = ", length( spline_spatial_$smooth)))
    }
    for(i in 1:length(spline_spatial_$smooth)) {
      S_null <- spline_spatial_$smooth[[i]]$S[[1]]
      for_plotting <- data.frame("temp_name" = seq(min(observed_df@data[,spline_spatial_covariates[i]]),max(observed_df@data[,spline_spatial_covariates[i]]),by = diff(range(observed_df@data[,spline_spatial_covariates[i]])) / 50));
      colnames(for_plotting) = spline_spatial_covariates[i]
      forReport <- mgcv::PredictMat(spline_spatial_$smooth[[i]], data = for_plotting)
      S_spatial_list[[i]] <- S_null
      S_spatial_reporting_list[[i]] <- forReport
    }
  } else {
    ## create a dummy variable
    ff = formula(paste0(response_variable_label," ~ s(",response_variable_label,", bs = 'ts')"))
    spline_spatial_ <- mgcv::gam(ff, data = observed_df@data, fit = F)
    if(trace_level == "high")
      print(paste0("length spline_spatial_ = ", length( spline_spatial_$smooth)))
    
    S_null <- spline_spatial_$smooth[[1]]$S[[1]]
    for_plotting <- seq(min(observed_df@data[,response_variable_label]),max(observed_df@data[,response_variable_label]),by = diff(range(observed_df@data[,response_variable_label])) / 50)
    plot_data = data.frame(x = for_plotting)
    colnames(plot_data) = response_variable_label
    
    forReport <- mgcv::PredictMat(spline_spatial_$smooth[[1]], data = plot_data)
    S_spatial_list[[1]] <- S_null
    S_spatial_reporting_list[[1]] <- forReport
  }
  if(trace_level != "none")
    print(paste0("Passed spatial spline section"))
  
  S_catchability_combined <- .bdiag(S_catchability_list)         # join S's in sparse matrix
  S_catchability_dims <- unlist(lapply(S_catchability_list, nrow)) # Find dimension of each S
  S_catchability_design_matrix <- .bdiag(S_catchability_reporting_list)
  S_spatial_combined <- .bdiag(S_spatial_list)         # join S's in sparse matrix
  S_spatial_dims <- unlist(lapply(S_spatial_list, nrow)) # Find dimension of each S
  S_spatial_design_matrix <- .bdiag(S_spatial_reporting_list)
  
  auxillary_objects = list(splinereportinglist = S_catchability_reporting_list, splinespatialreportinglist = S_spatial_reporting_list)
  
  if(trace_level != "none")
    print(paste0("Passed: model matrix configurations and spline configurations"))
  
  ## Projection model matrix stuff
  proj_df_subset = subset(projection_df, subset = projection_df@data[,time_variable_label] == time_levels[1])
  # number of spatial cells in projection matrix
  n_z = nrow(proj_df_subset)
  if(!set_up_dummy_proj) {
    if(apply_preferential_sampling & preference_model_type == 1) {
      if(sum(projection_raster_layer@data@values == 1) != n_z)
        stop(paste0("the object projection_raster_layer, needs to have the data values == 1, for each projection cell. Found sum(projection_raster_layer@data@values == 1) = ", sum(projection_raster_layer@data@values == 1) , " and n_z projection cells = ", n_z))
    }
  }
  
  
  X_spatial_proj_zpt = array(1, dim = c(n_z, ncol(data_spatial_model_matrix), n_t), dimnames = list(NULL, colnames(data_spatial_model_matrix), NULL))
  X_spatial_ipt = array(1, dim = c(n, ncol(data_spatial_model_matrix), n_t), dimnames = list(NULL, colnames(data_spatial_model_matrix), NULL))
  spline_spatial_model_matrix_proj_zpt = array(0, dim = c(n_z, ncol(spline_spatial_$X[,-1]), n_t))
  spline_spatial_model_matrix_ipt = array(spline_spatial_$X[,-1], dim = c(n, ncol(spline_spatial_$X[,-1]), n_t))
  Nij = array(0, dim = c(n_z, n_t))
  
  ## map mesh to extrapolation grid. Assumes projection_df has the same spatial cell for all time-cells
  Proj <- inla.mesh.projector(mesh, loc = coordinates(proj_df_subset))
  P = Proj$proj$A
  Proj_area = proj_df_subset@data$area
  
  proj_vertex_ndx = nn2(mesh$loc[,c(1,2)], coordinates(proj_df_subset), k = 1)$nn.idx
  data_vertex_ndx = nn2(mesh$loc[,c(1,2)], coordinates(observed_df), k = 1)$nn.idx
  
  
  if(trace_level != "none")
    print(paste0("Passed: mapping projection grid to mesh"))
  
  for(t in 1:n_t) {
    proj_df_subset <- subset(projection_df@data, subset = projection_df@data[,time_variable_label] == time_levels[t])
    data_df_subset <- subset(observed_df, subset = observed_df@data[,time_variable_label] == time_levels[t])
    data_df_subset@data$indicator = 1
    # validate binning
    #plot(proj_count)
    #points(coordinates(data_df_subset), pch = 16, col = adjustcolor(col = "red", alpha.f = 0.2))
    if(!set_up_dummy_proj) {
      proj_count <- rasterize(data_df_subset, projection_raster_layer, field = data_df_subset$indicator, sum, na.rm = T)
      Nij[,t] = proj_count$layer@data@values[projection_raster_layer@data@values == 1]
      if(sum(proj_count$layer@data@values[projection_raster_layer@data@values == 1], na.rm = T) != nrow(data_df_subset))
        stop(paste0("in time-step ", time_levels[t], " there were ",  nrow(data_df_subset), " observatons, but when we rasterised this data only calculated ", sum(proj_count$layer@data@values[projection_raster_layer@data@values == 1], na.rm = T) , " check the raster layer is correct dimensions"))
      
      #plot(proj_count)
      #points(data_df_subset, pch = 16, cex = 0.3)
      
      Nij[,t][is.na(Nij[,t])] = 0
    }
    
    
    if(length(spatial_covariates) > 0) {
      ##  build data spatial model matrix
      #data_spatial_model_matrix <- evalit(paste0("model.matrix(",response_variable_label," ~ 0 + ",paste(spatial_covariates, collapse = " + "),",  data = data@data)"), e1)
      ff = paste0(response_variable_label," ~ 0 + ",paste(spatial_covariates, collapse = " + "))
      m <- model.frame(formula(ff), observed_df@data, drop.unused.levels = F)
      data_spatial_model_matrix <- model.matrix(formula(ff), m)
      
      ##  build projection spatial model matrix
      ff <- paste0(response_variable_label," ~ 0 + ", paste(spatial_covariates, collapse = " + "))
      m <- model.frame(formula(ff), proj_df_subset, drop.unused.levels = F)
      proj_spatial_model_matrix <- model.matrix(formula(ff), m)
      if(ncol(proj_spatial_model_matrix) != ncol(data_spatial_model_matrix))
        stop("when building model matrix for spatial_covariates. The columns of projection model matrix differed from the columns of observed_df model matrix. This can happen when there are levels in one data set that are not in the other. you might want to investigate these covariates.")
      X_spatial_proj_zpt[,,t] <- proj_spatial_model_matrix
      X_spatial_ipt[,,t] <- data_spatial_model_matrix
    }
    if(length(spline_spatial_covariates) > 0) {
      ## for the spatial spline
      S_spatial_proj_model_mat = NULL
      for(i in 1:length(spline_spatial_covariates)) {
        spatial_spline_proj <- PredictMat(spline_spatial_$smooth[[i]], data = proj_df_subset)
        S_spatial_proj_model_mat <- cbind(S_spatial_proj_model_mat, spatial_spline_proj)
      }
      spline_spatial_model_matrix_proj_zpt[,,t] <- S_spatial_proj_model_mat
    }
  }
  if(trace_level != "none")
    print(paste0("Passed: projection model matrix construction"))
  
  if(is.null(kappa_bounds)) {
    Dist = stats::dist(mesh$loc)
    kappa_bounds = c(sqrt(8)/max(Dist),  sqrt(8)/min(Dist)) # Range = nu*sqrt(8)/kappa
  }
  
  ## Set up TMB object.
  tmb_data <- list(
    model = "SpatialTemporalCPUE",
    n_i = n,
    n_t = n_t,
    y_i = observed_df@data[,response_variable_label],
    t_i = time_variable - min(time_variable), # index for C++ language could be 1990 1991 etc,
    area_i = observed_df@data$area,
    obs_t = as.numeric(table(time_variable)),
    A = A,
    spde = spde$param.inla[c("M0","M1","M2")],
    Proj = P,
    Proj_Area = Proj_area,
    family = family,
    link = link,
    pref_coef_bounds = pref_bounds, ## perhaps make a user input, can be a difficult parameter to estimate
    kappa_bounds = kappa_bounds,
    Nij = Nij,
    apply_pref = ifelse(apply_preferential_sampling, 1, 0),
    LCGP_approach = preference_model_type, ## 0 = dinsdale, 1 = LGCP Lattice
    model_matrix = model_matrix,
    time_model_matrix = time_model_matrix,
    X_spatial_ipt = X_spatial_ipt,
    X_spatial_proj_zpt = X_spatial_proj_zpt,
    omega_indicator = ifelse(include_omega, 1, 0),
    epsilon_indicator = ifelse(include_epsilon, 1, 0),
    epsilon_ar1 = ifelse(epsilon_structure == "iid", 0, 1),
    spline_flag = c(ifelse(length(spline_catchability_covariates) > 0, 1, 0), ifelse(length(spline_spatial_covariates) > 0, 1, 0)),
    spline_model_matrix = spline_$X[,-1],
    spline_spatial_model_matrix_ipt = spline_spatial_model_matrix_ipt,
    spline_spatial_model_matrix_proj_zpt = spline_spatial_model_matrix_proj_zpt,
    S = S_catchability_combined,
    Sdims = S_catchability_dims,
    designMatrixForReport = S_catchability_design_matrix,
    S_spatial = S_spatial_combined,
    Sdims_spatial = S_spatial_dims,
    designMatrixForReport_spatial = S_spatial_design_matrix,
    spatial_constrained_coeff_ndx = spatial_constrained_coeff_ndx,
    spatial_covar_type = spatial_covar_type,
    simulate_GF = 0
    
  )
  year_ndx_for_each_obs = matrix(-99, nrow = tmb_data$n_t, ncol = max(tmb_data$obs_t))
  for(t in 1:tmb_data$n_t) {
    year_ndx_for_each_obs[t, 1:tmb_data$obs_t[t]] = which(tmb_data$t_i == (t - 1)) - 1
  }
  tmb_data$year_ndx_for_each_obs = year_ndx_for_each_obs
  
  if(linear_basis == 0) {
    tmb_data$A = A
    tmb_data$Proj = P
    tmb_data$model = "SpatialTemporalCPUE"
  } else if(linear_basis == 1) {
    tmb_data$index_proj_vertex = proj_vertex_ndx - 1 # R index -> C++ index
    tmb_data$index_data_vertex = data_vertex_ndx - 1 # R index -> C++ index
    tmb_data$model = "SpatialTemporalCPUENN"
    
  } else if(linear_basis == 2) {
    tmb_data$index_proj_vertex = proj_vertex_ndx - 1 # R index -> C++ index
    tmb_data$index_data_vertex = data_vertex_ndx - 1 # R index -> C++ index
    tmb_data$model = "SpatialTemporalCPUEVAST"
  }
  
  if(trace_level != "none")
    print(paste0("Passed: TMB data list construction"))
  
  # parameters for TMB
  dis_cor_0.1 = mean(0.5 * c(diff(range(mesh$loc[,1])), diff(range(mesh$loc[,2])))) # halfway 
  sigma_guess = 1 ## more intuative to think of marginal variance when setting starting values for params
  
  params = list(
    betas = rep(0, ncol(tmb_data$model_matrix)),
    constrained_spatial_betas = rep(0, max(tmb_data$spatial_constrained_coeff_ndx) + 1),
    constrained_time_betas = rep(0, max(dim(tmb_data$time_model_matrix)[2] - 1,1)),
    ln_phi = 0,
    logit_pref_coef = logit_general(0, tmb_data$pref_coef_bounds[1], tmb_data$pref_coef_bounds[2]),
    logit_pref_hyper_params = c(logit_pref_hyper_prior_vals[1], log(logit_pref_hyper_prior_vals[2])),
    lgcp_intercept = 0,
    ln_kappa_omega = log(sqrt(8) / dis_cor_0.1),#sqrt(8) / dis_cor_0.1,
    ln_tau_omega =  log(1/(sigma_guess * sqrt(4*pi * sqrt(8) / dis_cor_0.1^2))),
    ln_kappa_epsilon = log(sqrt(8) / dis_cor_0.1),#sqrt(8) / dis_cor_0.1,
    ln_tau_epsilon =  log(1/(sigma_guess * sqrt(4*pi * sqrt(8) / dis_cor_0.1^2))),
    trans_eps_rho = 0,
    omega_input = rep(0,spde$n.spde),
    epsilon_input = array(0,dim = c(spde$n.spde, tmb_data$n_t)),
    gammas = rep(0,sum(tmb_data$Sdims)),  # Spline coefficients
    ln_lambda = rep(0,length(tmb_data$Sdims)), #Log spline penalization coefficients
    gammas_spatial = rep(0,sum(tmb_data$Sdims_spatial)),  # Spline coefficients
    ln_lambda_spatial = rep(0,length(tmb_data$Sdims_spatial)) #Log spline penalization coefficients
  )
  params_vast = list(
    betas = rep(0, ncol(tmb_data$model_matrix)),
    constrained_spatial_betas = rep(0, max(tmb_data$spatial_constrained_coeff_ndx) + 1),
    constrained_time_betas = rep(0, max(dim(tmb_data$time_model_matrix)[2] - 1,1)),
    ln_phi = 0,
    logit_pref_coef = logit_general(0, tmb_data$pref_coef_bounds[1], tmb_data$pref_coef_bounds[2]),
    logit_pref_hyper_params = c(logit_pref_hyper_prior_vals[1], log(logit_pref_hyper_prior_vals[2])),
    lgcp_intercept = 0,
    logit_kappa = 0,#sqrt(8) / dis_cor_0.1,
    omega_eta = 1,
    epsilon_eta = 1,
    trans_eps_rho = 0,
    omega_input = rep(0,spde$n.spde),
    epsilon_input = array(0,dim = c(spde$n.spde, tmb_data$n_t)),
    gammas = rep(0,sum(tmb_data$Sdims)),  # Spline coefficients
    ln_lambda = rep(0,length(tmb_data$Sdims)), #Log spline penalization coefficients
    gammas_spatial = rep(0,sum(tmb_data$Sdims_spatial)),  # Spline coefficients
    ln_lambda_spatial = rep(0,length(tmb_data$Sdims_spatial)) #Log spline penalization coefficients
  )
  if(pref_hyper_distribution > 0) {
    params$logit_pref_coef =  rep(logit_general(0, tmb_data$pref_coef_bounds[1], tmb_data$pref_coef_bounds[2]), n_t)
	params_vast$logit_pref_coef =  rep(logit_general(0, tmb_data$pref_coef_bounds[1], tmb_data$pref_coef_bounds[2]), n_t)
  }
  
  if(trace_level != "none")
    print(paste0("Passed: TMB parameter list construction"))
  
  any_errors = check_inputs(tmb_data = tmb_data, tmb_params = params)
  if(!any_errors$result) {
    print("Error found incompatible data and parameter lists. I am returning the errors to help debug this function.")
    return(list(errors = any_errors$errors, tmb_pars = params, tmb_data = tmb_data))
  }
  ## set which parameters are being estimated and which are held constant.drr
  pars_to_fix = NULL;
  if(linear_basis == 2) {
    if(epsilon_structure == "iid")
      pars_to_fix = c(pars_to_fix, "trans_eps_rho")
    if(!include_epsilon & !include_omega)
      pars_to_fix = c(pars_to_fix, "epsilon_input", "logit_kappa", "omega_input", "omega_eta", "epsilon_eta")
    if(!include_omega)
      pars_to_fix = c(pars_to_fix, "omega_input", "omega_eta")
    if(!include_epsilon)
      pars_to_fix = c(pars_to_fix, "epsilon_input", "epsilon_eta")
    if(length(spatial_covariates) == 0)
      pars_to_fix = c(pars_to_fix, "constrained_spatial_betas")
    if(length(spline_catchability_covariates) == 0)
      pars_to_fix = c(pars_to_fix, "ln_lambda", "gammas")
    if(length(spline_spatial_covariates) == 0)
      pars_to_fix = c(pars_to_fix, "ln_lambda_spatial", "gammas_spatial")
    ## do we need to estimate a dispersion parameter
    if(family %in% c(1,3)) ## binomial and Poisson
      pars_to_fix = c(pars_to_fix, "ln_phi")
  } else {
    if(epsilon_structure == "iid")
      pars_to_fix = c(pars_to_fix, "trans_eps_rho")
    if(!include_epsilon)
      pars_to_fix = c(pars_to_fix, "epsilon_input", "ln_kappa_epsilon", "ln_tau_epsilon")
    if(!include_omega)
      pars_to_fix = c(pars_to_fix, "omega_input", "ln_kappa_omega", "ln_tau_omega")
    if(length(spatial_covariates) == 0)
      pars_to_fix = c(pars_to_fix, "constrained_spatial_betas")
    if(length(spline_catchability_covariates) == 0)
      pars_to_fix = c(pars_to_fix, "ln_lambda", "gammas")
    if(length(spline_spatial_covariates) == 0)
      pars_to_fix = c(pars_to_fix, "ln_lambda_spatial", "gammas_spatial")
    ## do we need to estimate a dispersion parameter
    if(family %in% c(1,3)) ## binomial and Poisson
      pars_to_fix = c(pars_to_fix, "ln_phi")
  }

  vec_pars_to_adjust = NULL
  vec_elements_to_exclude = NULL
  if (apply_preferential_sampling) {
    if(preference_model_type == 0)
      pars_to_fix = c(pars_to_fix, "lgcp_intercept")
    if(pref_hyper_distribution < 0) {
      # not time-varying
      pars_to_fix = c(pars_to_fix, "logit_pref_coef", "logit_pref_hyper_params")
    } else if(pref_hyper_distribution == 0) {
      # not time-varying
      pars_to_fix = c(pars_to_fix, "logit_pref_hyper_params")
    } else if(pref_hyper_distribution == 2) {
      ## both fixed
      pars_to_fix = c(pars_to_fix, "logit_pref_hyper_params")
    } else if (pref_hyper_distribution == 3) {
      ## estimate hyper mean
      pars_to_fix = c(pars_to_fix, "logit_pref_hyper_params")
      vec_pars_to_adjust = c(vec_pars_to_adjust,"logit_pref_hyper_params")
      vec_elements_to_exclude = list(logit_pref_hyper_params = c(1))
    } else if(pref_hyper_distribution == 4) {
      ## estimate hyper sd
      pars_to_fix = c(pars_to_fix, "logit_pref_hyper_params")
      vec_pars_to_adjust = c(vec_pars_to_adjust,"logit_pref_hyper_params")
      vec_elements_to_exclude = list(logit_pref_hyper_params = c(2))
    }
  } else {
    pars_to_fix = c(pars_to_fix, "logit_pref_coef","lgcp_intercept", "logit_pref_hyper_params")
  }
  
  
  fixed_pars = list()
  if(linear_basis == 2) {
    if(length(pars_to_fix) > 0)
      fixed_pars = fix_pars(par_list = params_vast, pars_to_exclude = pars_to_fix, vec_pars_to_adjust = vec_pars_to_adjust, vec_elements_to_exclude = vec_elements_to_exclude)
  } else {
    if(length(pars_to_fix) > 0)
      fixed_pars = fix_pars(par_list = params, pars_to_exclude = pars_to_fix, vec_pars_to_adjust = vec_pars_to_adjust, vec_elements_to_exclude = vec_elements_to_exclude)
  }
  
  ## compile model
  #compile(file.path("src","debug_standalone_version.cpp"))
  #file.exists(file.path("src","debug_standalone_version.dll"))
  #dyn.load(dynlib(file.path("src","debug_standalone_version")))
  ## create ob
  to_be_silent = ifelse(trace_level == "none", TRUE, FALSE)
  random_pars = NULL
  if(include_epsilon)
    random_pars = c(random_pars, "epsilon_input")
  if(include_omega)
    random_pars = c(random_pars, "omega_input")
  if(apply_preferential_sampling & pref_hyper_distribution != 0)
    random_pars = c(random_pars, "logit_pref_coef")
  
  if(length(spline_catchability_covariates) != 0)
    random_pars = c(random_pars, "gammas")
  if(length(spline_spatial_covariates) == 0)
    random_pars = c(random_pars, "gammas_spatial")
  
  obj = NULL;
  tmb_pars = NULL
  if(linear_basis == 2) {
	tmb_pars = params_vast
    if(!is.null(init_vals)) {
      if(!all(unique(names(init_vals)) %in% unique(names(params_vast)))) {
        stop(paste0("init_vals: doesn't have all the parameter labels expected. expected ", paste(unique(names(params_vast)), collapse = "\n")))
      } else {
        params_vast = init_vals
      }
    }
    obj <- tryCatch(expr = MakeADFun(tmb_data, params_vast, random = random_pars, map = fixed_pars, DLL = "CPUEspatial_TMBExports", method = "nlminb", hessian = T, silent = to_be_silent), error = function(e){e})
    if(inherits(obj, "error")) {
      stop("An error occured in MakeADFun, probably and incompatiability between data and parameters")
    }
	
  } else {
  	tmb_pars = params
	if(!is.null(init_vals)) {
	  if(!all(unique(names(init_vals)) %in% unique(names(params)))) {
		stop(paste0("init_vals: doesn't have all the parameter labels expected. expected ", paste(unique(names(params)), collapse = "\n")))
	  } else {
		params = init_vals
	  }
	}
    obj <- tryCatch(expr = MakeADFun(tmb_data, params, random = random_pars, map = fixed_pars, DLL = "CPUEspatial_TMBExports", method = "nlminb", hessian = T, silent = to_be_silent), error = function(e){e})
    if(inherits(obj, "error")) { 
      stop("An error occured in MakeADFun, probably and incompatiability between data and parameters")
    }
  }
  
  if(trace_level != "none")
    print(paste0("Passed: successfully built obj"))
  Call$link = link
  Call$family = family
  return(list(obj = obj, tmb_pars = tmb_pars, tmb_data = tmb_data, Call = Call, auxillary_objects = auxillary_objects))
}


