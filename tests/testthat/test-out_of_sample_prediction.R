#'
#' @details Test out of sample prediction functionality
#'

#' A non spatial GLM model with  Gamma log link response variable
test_that("gamma_glm", {
  load(system.file("testdata", "non_spatial_glm.RData",package="CPUEspatial"))
  shape = 50 ### for gamma response variable
  
  ## simulate data 
  data@data$y_i = rgamma(n = nrow(data), shape = shape, scale = (data$area * exp(data$eta)) / shape)
  n = nrow(data@data)
  train_ndx = sample(1:n, size = round(n * 0.8), replace = F)
  train_df = subset(data, subset = 1:n %in% train_ndx)
  test_df = subset(data, subset =!1:n %in% train_ndx)
  
  ## build model
  simple_model = configure_obj(observed_df = train_df, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                               spatial_covariates = NULL, spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none")

  ## estimate model
  opt = nlminb(simple_model$obj$par, simple_model$obj$fn, simple_model$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  MLE_ = simple_model$obj$report(simple_model$obj$env$last.par.best)
  ## 
  self_test = predict_out_of_sample(conf_obj = simple_model, out_of_sample_df = train_df, observed_df = train_df, mesh = mesh, simulate = F)
  expect_equal(object =sum(MLE_$mu), sum(self_test$fitted),  tolerance = 0.001)
  test_pred = predict_out_of_sample(conf_obj = simple_model, out_of_sample_df = test_df, observed_df = train_df, mesh = mesh, simulate = F)
  expect_true(length(test_pred) == 3)
})


#' A non spatial GLM model with  Poisson log link response variable
test_that("poisson_glm", {
  load(system.file("testdata", "non_spatial_glm.RData",package="CPUEspatial"))
  ## simulate data 
  data$area = rlnorm(nrow(data), 0, 0.1)
  data@data$y_i = rpois(n = nrow(data), lambda = (data$area *3 * exp(data$eta)) )
  n = nrow(data@data)
  train_ndx = sample(1:n, size = round(n * 0.8), replace = F)
  train_df = subset(data, subset = 1:n %in% train_ndx)
  test_df = subset(data, subset =!1:n %in% train_ndx)
  
  ## build model
  simple_model = configure_obj(observed_df = train_df, projection_df = full_proj_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx",
                               spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none")
  ## estimate model
  opt = nlminb(simple_model$obj$par, simple_model$obj$fn, simple_model$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  ## get derived quantities
  rep_CPUE_spatial = simple_model$obj$report(simple_model$obj$env$last.par.best)
  ## 
  self_test = predict_out_of_sample(conf_obj = simple_model, out_of_sample_df = train_df, observed_df = train_df, mesh = mesh, simulate = F)
  expect_equal(object =sum(rep_CPUE_spatial$mu), sum(self_test$fitted),  tolerance = 0.001)
  test_pred = predict_out_of_sample(conf_obj = simple_model, out_of_sample_df = test_df, observed_df = train_df, mesh = mesh, simulate = F)
  expect_true(length(test_pred) == 3)
  
})


#' A non spatial GLM model with  Normal response variable
test_that("normal_glm", {
  load(system.file("testdata", "non_spatial_glm.RData",package="CPUEspatial"))
  std_dev = 0.0001 ### for gamma response variable
  ## simulate data 
  data@data$y_i = rnorm(n = nrow(data), mean = data$area * data$eta, sd = std_dev)
  n = nrow(data@data)
  train_ndx = sample(1:n, size = round(n * 0.8), replace = F)
  train_df = subset(data, subset = 1:n %in% train_ndx)
  test_df = subset(data, subset =!1:n %in% train_ndx)
  
  ## build model
  simple_model = configure_obj(observed_df = train_df, projection_df = full_proj_df, mesh = mesh, family = 0, link = 4, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                               spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none")
## estimate model
  opt = nlminb(simple_model$obj$par, simple_model$obj$fn, simple_model$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  #opt$convergence
  ## get derived quantities
  rep_CPUE_spatial = simple_model$obj$report(simple_model$obj$env$last.par.best)
  ## 
  self_test = predict_out_of_sample(conf_obj = simple_model, out_of_sample_df = train_df, observed_df = train_df, mesh = mesh, simulate = F)
  expect_equal(object =sum(rep_CPUE_spatial$mu), sum(self_test$fitted),  tolerance = 0.001)
  test_pred = predict_out_of_sample(conf_obj = simple_model, out_of_sample_df = test_df, observed_df = train_df, mesh = mesh, simulate = F)
  expect_true(length(test_pred) == 3)
})
