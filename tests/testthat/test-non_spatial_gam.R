#'
#' @details The purpose of these test are to test likelihood consistency and also estimation is consistent with other packages i.e glm
#'



#' A non spatial GAM model with  Gamma log link response variable
test_that("normal_gam", {
  load(system.file("testdata", "non_spatial_glm.RData",package="CPUEspatial"))
  library(mgcv)
  f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^3 * (1 - x)^10
  set.seed(123)
  ysim <- function(n = 500, scale = 2) {
    x <- runif(n)
    t = sample(x = c(1,2), replace = T, size = n, prob = c(0.3,0.7))
    t_coef = c(0.8,0.2)
    e <- rnorm(n, 0, scale)
    f <- f2(x)
    y <- t_coef[t] + f + e
    data.frame(y_i = y, x2 = x, f2 = f, t = t)
  }
  
  df <- ysim()
  df$t_factor = factor(df$t)
  df$t = as.integer(df$t)
  m <- gam(y_i ~ t_factor + s(x2, bs = "ts"), data = df, method="REML")
  df$long = rnorm(nrow(df), 180, 3)
  df$lat = rnorm(nrow(df), -44, 3)
  df$area = 1
  coordinates(df) <- ~ long + lat
  ## build model
  head(full_proj_df@data)
  full_proj_df = subset(full_proj_df, subset = full_proj_df$year %in% c(1,2))
  full_proj_df$t = full_proj_df$year
  full_proj_df$x2 = 1
  simple_model = configure_obj(observed_df = df, projection_df = full_proj_df, mesh = mesh, family = 0, link = 4, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "t", catchability_covariates = NULL, 
                               spatial_covariates = NULL, spline_catchability_covariates = c("x2"),
                               spline_spatial_covariates = NULL, trace_level = "none", linear_basis = 0)
  
  ## estimate model
  opt = nlminb(simple_model$obj$par, simple_model$obj$fn, simple_model$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  ## 
  rep_CPUE_spatial = simple_model$obj$report(simple_model$obj$env$last.par.best)
  ## compare fitted values
  expect_equal(rep_CPUE_spatial$mu, fitted(m), tolerance = 0.1)
  
  ## estimate phi a little off but close
  ## summary(m)$dispersion rep_CPUE_spatial$phi
  
  ## test predict out of sample function
  test_predict = predict_out_of_sample(conf_obj = simple_model, out_of_sample_df = df, observed_df = df, mesh = mesh)
  expect_equal(rep_CPUE_spatial$mu, as.numeric(test_predict$fitted), tolerance = 0.001)
  
  ## predict on new dataset
  test_df <- ysim(n = 10)
  test_df$t_factor = factor(test_df$t)
  test_df$t = as.integer(test_df$t)
  test_df$long = rnorm(nrow(test_df), 180, 3)
  test_df$lat = rnorm(nrow(test_df), -44, 3)
  test_df$area = 1
  coordinates(test_df) <- ~ long + lat
  test_predict = predict_out_of_sample(conf_obj = simple_model, out_of_sample_df = test_df, observed_df = df, mesh = mesh)
  gam_mgcv_predict = predict(m, type = "response", newdata = test_df)
  expect_equal(as.numeric(test_predict$fitted), as.numeric(gam_mgcv_predict), tolerance = 0.001)
  
})


#' A non spatial GAM model with normal error multiple splines for catchability
test_that("normal_gam_multiple_smoothers", {
  load(system.file("testdata", "non_spatial_glm.RData",package="CPUEspatial"))
  library(mgcv)
  f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^3 * (1 - x)^10
  f3 = function(x) {sin(x)}
  set.seed(123)
  ysim <- function(n = 500, scale = 2) {
    x <- runif(n)
    x3 <- runif(n)
    t = sample(x = c(1,2), replace = T, size = n, prob = c(0.3,0.7))
    t_coef = c(0.8,0.2)
    e <- rnorm(n, 0, scale)
    f <- f2(x)
    f_3 = f3(x3)
    y <- t_coef[t] + f + f_3 + e
    data.frame(y_i = y, x2 = x, f2 = f, t = t, x3)
  }
  
  df <- ysim()
  df$t_factor = factor(df$t)
  df$t = as.integer(df$t)
  m <- gam(y_i ~ t_factor + s(x2, bs = "ts")+ s(x3, bs = "ts"), data = df, method="REML")
  df$long = rnorm(nrow(df), 180, 3)
  df$lat = rnorm(nrow(df), -44, 3)
  df$area = 1
  coordinates(df) <- ~ long + lat
  ## build model
  head(full_proj_df@data)
  full_proj_df = subset(full_proj_df, subset = full_proj_df$year %in% c(1,2))
  full_proj_df$t = full_proj_df$year
  full_proj_df$x2 = 1
  simple_model = configure_obj(observed_df = df, projection_df = full_proj_df, mesh = mesh, family = 0, link = 4, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "t", catchability_covariates = NULL, 
                               spatial_covariates = NULL, spline_catchability_covariates = c("x2", "x3"),
                               spline_spatial_covariates = NULL, trace_level = "none", linear_basis = 0)
  
  ## estimate model
  opt = nlminb(simple_model$obj$par, simple_model$obj$fn, simple_model$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  ## 
  rep_CPUE_spatial = simple_model$obj$report(simple_model$obj$env$last.par.best)
  ## compare fitted values
  expect_equal(rep_CPUE_spatial$mu, fitted(m), tolerance = 0.1)
  
  ## estimate phi a little off but close
  ## summary(m)$dispersion rep_CPUE_spatial$phi
  
  ## test predict out of sample function
  test_predict = predict_out_of_sample(conf_obj = simple_model, out_of_sample_df = df, observed_df = df, mesh = mesh)
  expect_equal(rep_CPUE_spatial$mu, as.numeric(test_predict$fitted), tolerance = 0.001)
  
  
  ## predict on new dataset
  test_df <- ysim(n = 10)
  test_df$t_factor = factor(test_df$t)
  test_df$t = as.integer(test_df$t)
  test_df$long = rnorm(nrow(test_df), 180, 3)
  test_df$lat = rnorm(nrow(test_df), -44, 3)
  test_df$area = 1
  coordinates(test_df) <- ~ long + lat
  test_predict = predict_out_of_sample(conf_obj = simple_model, out_of_sample_df = test_df, observed_df = df, mesh = mesh)
  gam_mgcv_predict = predict(m, type = "response", newdata = test_df)
  expect_equal(as.numeric(test_predict$fitted), as.numeric(gam_mgcv_predict), tolerance = 0.001)
  
})
#' A non spatial GLM model with  Poisson log link response variable
test_that("gamma_gam", {
  load(system.file("testdata", "non_spatial_glm.RData",package="CPUEspatial"))
  shape = 50 ### for gamma response variable
  library(mgcv)
  f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^2 * (1 - x)^10
  set.seed(123)
  ysim <- function(n = 500, scale = 2) {
    x <- runif(n)
    t = sample(x = c(1,2), replace = T, size = n, prob = c(0.3,0.7))
    t_coef = c(0.8,0.2)
    area = rlnorm(n, log(1), 0.5)
    f <- f2(x)
    eta = t_coef[t] + f 
    y_i = rgamma(n = n, shape = shape, scale = (area * exp(eta)) / shape)
    data.frame(y_i = y_i, x2 = x, f2 = f, t = t, eta = eta, area)
  }
  
  df <- ysim()
  df$t_factor = factor(df$t)
  df$t = as.integer(df$t)
  gam_mod <- gam(y_i ~ offset(log(area)) + t_factor + s(x2, bs = "ts"), data = df, family = Gamma(link = log), method="REML", control  = gam.control(maxit = 10000))
  
  df$long = rnorm(nrow(df), 180, 3)
  df$lat = rnorm(nrow(df), -44, 3)
  coordinates(df) <- ~ long + lat
  ## build model
  full_proj_df = subset(full_proj_df, subset = full_proj_df$year %in% c(1,2))
  full_proj_df$t = full_proj_df$year
  full_proj_df$x2 = 1
  ## build model
  simple_model = configure_obj(observed_df = df, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = F, 
                               response_variable_label = "y_i", time_variable_label = "t", catchability_covariates = NULL,
                               spatial_covariates = NULL,  spline_catchability_covariates = c("x2"),
                               spline_spatial_covariates = NULL, trace_level = "none")
  
  
  ## estimate model
  opt = nlminb(simple_model$obj$par, simple_model$obj$fn, simple_model$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  #opt$convergence
  ## get derived quantities
  rep_CPUE_spatial = simple_model$obj$report(simple_model$obj$env$last.par.best)
  expect_equal(rep_CPUE_spatial$mu, fitted(gam_mod), tolerance = 0.1)
  #expect_equal(object = -1 * opt$objective, as.numeric(logLik(gam_mod)),  tolerance = 0.001)
  ## 1/summary(gam_mod)$dispersion rep_CPUE_spatial$phi
  
  ## test predict function
  test_predict = predict_out_of_sample(conf_obj = simple_model, out_of_sample_df = df, observed_df = df, mesh = mesh)
  expect_equal(rep_CPUE_spatial$mu, as.numeric(test_predict$fitted), tolerance = 0.001)
  
})
