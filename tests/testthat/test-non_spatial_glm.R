#'
#' @details The purpose of these test are to test likelihood consistency and also estimation is consistent with other packages i.e glm
#'



#' A non spatial GLM model with  Gamma log link response variable
test_that("gamma_glm", {
  load(system.file("testdata", "non_spatial_glm.RData",package="CPUEspatial"))
  shape = 50 ### for gamma response variable
  
  ## simulate data 
  data@data$y_i = rgamma(n = nrow(data), shape = shape, scale = (data$area * exp(data$eta)) / shape)
  
  ## build model
  simple_model = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                               spatial_covariates = NULL, spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none", linear_basis = 0)
  simple_model_NN = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                               spatial_covariates = NULL, spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none", linear_basis = 1)
  
  ## estimate model
  opt = nlminb(simple_model$obj$par, simple_model$obj$fn, simple_model$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  opt_NN = nlminb(simple_model_NN$obj$par, simple_model_NN$obj$fn, simple_model_NN$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  #opt$convergence
  ## get derived quantities
  rep_CPUE_spatial_NN = simple_model_NN$obj$report(simple_model_NN$obj$env$last.par.best)
  rep_CPUE_spatial = simple_model$obj$report(simple_model$obj$env$last.par.best)
  simple_glm = glm(y_i ~  factor(year) + factor(fleet_ndx), offset = log(area), data = data, family = Gamma(link = log))
  expect_equal(object = deviance_calc(y = data$y_i, mu = rep_CPUE_spatial$mu, dist = "gamma"),
               expected = summary(simple_glm)$deviance, 
               tolerance = 0.001)
  expect_equal(object = opt$objective, opt_NN$objective,  tolerance = 0.001)
  expect_equal(object = -1 * opt$objective, as.numeric(logLik(simple_glm)),  tolerance = 0.001)
  
})


#' A non spatial GLM model with  Poisson log link response variable
test_that("poisson_glm", {
  load(system.file("testdata", "non_spatial_glm.RData",package="CPUEspatial"))
  ## simulate data 
  data$area = rlnorm(nrow(data), 0, 0.1)
  data@data$y_i = rpois(n = nrow(data), lambda = (data$area *3 * exp(data$eta)) )
  
  ## build model
  simple_model = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx",
                               spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none")
  simple_model_NN = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                  response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                                  spatial_covariates = NULL, spline_catchability_covariates = NULL,
                                  spline_spatial_covariates = NULL, trace_level = "none", linear_basis = 1)
  
  ## estimate model
  opt = nlminb(simple_model$obj$par, simple_model$obj$fn, simple_model$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  opt_NN = nlminb(simple_model_NN$obj$par, simple_model_NN$obj$fn, simple_model_NN$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  #opt$convergence
  ## get derived quantities
  rep_CPUE_spatial = simple_model_NN$obj$report(simple_model_NN$obj$env$last.par.best)
  rep_CPUE_spatial_NN = simple_model$obj$report(simple_model$obj$env$last.par.best)
  simple_glm = glm(y_i ~  factor(year) + factor(fleet_ndx), offset = log(area), data = data, family = poisson(link = log))
  
  expect_equal(object = opt$objective, opt_NN$objective,  tolerance = 0.001)
  expect_equal(object = -1 * opt$objective, as.numeric(logLik(simple_glm)),  tolerance = 0.001)
  
})


#' A non spatial GLM model with  Normal response variable
test_that("normal_glm", {
  load(system.file("testdata", "non_spatial_glm.RData",package="CPUEspatial"))
  std_dev = 0.0001 ### for gamma response variable
  ## simulate data 
  data@data$y_i = rnorm(n = nrow(data), mean = data$area * data$eta, sd = std_dev)
  data@data$y_per_area =  data@data$y_i /  data@data$area
  full_proj_df@data$y_per_area = 1
  glm_version_data = data
  glm_version_data@data$area = rnorm(n = nrow(data), 1, 0.01)
  ## build model
  simple_model = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 0, link = 4, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                               spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none")
  simple_model_NN = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 0, link = 4, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx",
                               spatial_covariates = NULL, spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none", linear_basis = 1)
  
  simple_glm_model = configure_obj(observed_df = glm_version_data, projection_df = full_proj_df, mesh = mesh, family = 0, link = 4, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_per_area", time_variable_label = "year", catchability_covariates = "fleet_ndx",
                               spatial_covariates = NULL, spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none")
  ## estimate model
  opt = nlminb(simple_model$obj$par, simple_model$obj$fn, simple_model$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  opt_NN = nlminb(simple_model_NN$obj$par, simple_model_NN$obj$fn, simple_model_NN$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  opt_glm = nlminb(simple_glm_model$obj$par, simple_glm_model$obj$fn, simple_glm_model$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  #opt$convergence
  ## get derived quantities
  rep_CPUE_spatial = simple_model$obj$report(simple_model$obj$env$last.par.best)
  rep_CPUE_spatial_glm = simple_glm_model$obj$report(simple_glm_model$obj$env$last.par.best)
  
  simple_glm = glm(y_per_area ~  factor(year) + factor(fleet_ndx), data = data, family = gaussian(link = "identity"))
  #expect_equal(object = -1 * opt_glm$objective, expected = as.numeric(logLik(simple_glm)),  tolerance = 0.001)
  expect_equal(object = opt$objective, expected = opt_NN$objective,  tolerance = 0.001)
  
  ## we dont' expect the log likeilihood, deviance and estimated standard devaition.
  ## when area is treated as a multiplier on the response variable in CPUEspatial compared with the catch/area in the GLM
  ## look at coeffecients particularly conanical index
  res <- as.data.frame(coefficients(simple_glm))
  index <- regexpr("year", row.names(res)) > 0
  X <- res[index, 1]
  V <- summary(simple_glm)$cov.unscaled[index, , drop = FALSE][, index]
  n <- length(X) + 1
  A <- matrix(-1/n, n, n - 1)
  A[-1, ] <- A[-1, ] + diag(rep(1, n - 1))
  CPUE <- A %*% X
  COV <- sqrt(diag((A %*% V) %*% t(A)))
  index <- exp(CPUE)
  geometric_mean_CPUE_spatial = rep_CPUE_spatial$relative_index / exp(sum(log(rep_CPUE_spatial$relative_index[rep_CPUE_spatial$relative_index > 0]), na.rm=T) / length(rep_CPUE_spatial$relative_index)) 
  for(i in 1:length(geometric_mean_CPUE_spatial))
    expect_equal(object = geometric_mean_CPUE_spatial[i], expected = index[i],  tolerance = 0.1)
  
})

#' A non spatial GLM model with  normal - log response variable
test_that("lognormal_glm", {
  load(system.file("testdata", "non_spatial_glm.RData",package="CPUEspatial"))
  std_dev = 0.15
  data@data$y_i = rlnorm(n = nrow(data), meanlog = log(data$area * exp(data$eta)), std_dev)
  simple_model = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 0, link = 0, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                               spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none")
  simple_model_NN = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 0, link = 0, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx",  
                               spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none", linear_basis = 1)
  ## estimate model
  opt = nlminb(simple_model$obj$par, simple_model$obj$fn, simple_model$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  opt_NN = nlminb(simple_model_NN$obj$par, simple_model_NN$obj$fn, simple_model_NN$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  #opt$convergence
  ## get derived quantities
  rep_CPUE_spatial = simple_model$obj$report(simple_model$obj$env$last.par.best)
  simple_glm = glm(y_i ~  factor(year) + factor(fleet_ndx), data = data, offset = log(area), family = gaussian(link = "log"))
  ## compare log likelihoods from MLEs
  expect_equal(object = -1 * opt$objective, expected = as.numeric(logLik(simple_glm)),  tolerance = 0.001)
  expect_equal(object = opt$objective, expected = opt_NN$objective,  tolerance = 0.001)
  
})


#' A non spatial GLM model with  binomial response variable
test_that("binomial_glm", {
  load(system.file("testdata", "non_spatial_glm.RData",package="CPUEspatial"))
  data@data$area = rnorm(nrow(data))
  data@data$y_i = rbinom(n = nrow(data), size = 1, prob = plogis(data$eta))
  #table(data@data$y_i) 
  simple_model = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 1, link = 1, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx",
                               spatial_covariates = NULL, spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none")
  simple_model_NN = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 1, link = 1, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                               spatial_covariates = NULL, spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none", linear_basis = 1)
  ## estimate model
  opt = nlminb(simple_model$obj$par, simple_model$obj$fn, simple_model$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  opt_NN = nlminb(simple_model_NN$obj$par, simple_model_NN$obj$fn, simple_model_NN$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  #opt$convergence
  ## get derived quantities
  rep_CPUE_spatial = simple_model$obj$report(simple_model$obj$env$last.par.best)
  rep_CPUE_spatial_NN = simple_model_NN$obj$report(simple_model_NN$obj$env$last.par.best)
  
  simple_glm = glm(y_i ~  factor(year) + factor(fleet_ndx), data = data, family = binomial(link = "logit"))
  ## compare log likelihoods from MLEs
  expect_equal(object = -1 * opt$objective, expected = as.numeric(logLik(simple_glm)),  tolerance = 0.001)
  
  #sum(dbinom(x = data@data$y_i, size = 1, prob = rep_CPUE_spatial$mu, log =T))
  #sum(dbinom(x = data@data$y_i, size = 1, prob = fitted(simple_glm), log =T))
  expect_equal(object = opt$objective, expected = opt_NN$objective,  tolerance = 0.001)
})
  
  
#' A non spatial GLM model with  binomial response variable
test_that("negative_binomial_glm", {
  library(MASS)
  load(system.file("testdata", "non_spatial_glm.RData",package="CPUEspatial"))
  #data@data$area = rnorm(nrow(data))
  phi = 1 # 0 implies poisson no overdispersion
  mu = data$area * exp(data$eta)
  var = mu + (mu^2/phi)
  data@data$y_i = rnbinom(nrow(data), size = phi, mu = mu)
  
  #table(data@data$y_i) 
  simple_model = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 4, link = 0, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx",
                               spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none")
  simple_model_NN = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 4, link = 0, include_omega = F, include_epsilon = F, 
                                  response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                                  spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                                  spline_spatial_covariates = NULL, trace_level = "none", linear_basis = 1)
  ## estimate model
  opt = suppressWarnings(nlminb(simple_model$obj$par, simple_model$obj$fn, simple_model$obj$gr, control = list(eval.max = 10000, iter.max = 10000)))
  opt_NN = suppressWarnings( nlminb(simple_model_NN$obj$par, simple_model_NN$obj$fn, simple_model_NN$obj$gr, control = list(eval.max = 10000, iter.max = 10000)))
  #opt$convergence
  ## get derived quantities
  rep_CPUE_spatial = simple_model$obj$report(simple_model$obj$env$last.par.best)
  rep_CPUE_spatial_NN = simple_model_NN$obj$report(simple_model_NN$obj$env$last.par.best)
  
  simple_glm = suppressWarnings(glm.nb(y_i ~  factor(year) + factor(fleet_ndx) + offset(log(area)), data = data, link = "log"))
  ## compare log likelihoods from MLEs
  expect_equal(object = -1 * opt$objective, expected = as.numeric(logLik(simple_glm)),  tolerance = 0.001)
  
  #sum(dbinom(x = data@data$y_i, size = 1, prob = rep_CPUE_spatial$mu, log =T))
  #sum(dbinom(x = data@data$y_i, size = 1, prob = fitted(simple_glm), log =T))
  expect_equal(object = opt$objective, expected = opt_NN$objective,  tolerance = 0.001)
})
