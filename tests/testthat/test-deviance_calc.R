#'
#' @details The purpose are test configuration with a range of
#' covariates configurations.
#'

#' A non spatial GLM model with  Gamma log link response variable
test_that("deviance", {
  ## gamma
  load(system.file("testdata", "non_spatial_glm.RData",package="CPUEspatial"))
  shape = 50 ### for gamma response variable
  data@data$y_i = rgamma(n = nrow(data), shape = shape, scale = (data$area * exp(data$eta)) / shape)
  simple_glm = glm(y_i ~  factor(year) + factor(fleet_ndx), offset = log(area), data = data, family = Gamma(link = log))
  pkg_dev = deviance_calc(y =  data@data$y_i, mu = fitted(simple_glm), dist = "gamma")
  glm_dev = deviance(simple_glm) 
  expect_equal(object = pkg_dev, expected = glm_dev,  tolerance = 0.001)
  ## Gaussian
  std_dev = 0.0001 ### for gamma response variable
  data@data$y_i = rnorm(n = nrow(data), mean = data$area * data$eta, sd = std_dev)
  data@data$y_per_area =  data@data$y_i /  data@data$area
  simple_glm = glm(y_per_area ~  factor(year) + factor(fleet_ndx), data = data, family = gaussian(link = "identity"))
  pkg_dev = deviance_calc(y =  data@data$y_per_area, mu = fitted(simple_glm), dist = "gaussian")
  glm_dev = deviance(simple_glm) 
  expect_equal(object = pkg_dev, expected = glm_dev,  tolerance = 0.001)
  ## negative binomial
  library(MASS)
  phi = 1
  mu = data$area * exp(data$eta)
  var = mu + (mu^2/phi)
  data@data$y_i = rnbinom(nrow(data), size = phi, mu = mu)
  simple_glm = suppressWarnings(glm.nb(y_i ~  factor(year) + factor(fleet_ndx) + offset(log(area)), data = data, link = "log"))

  pkg_dev = deviance_calc(y =  data@data$y_i, mu = fitted(simple_glm), phi = simple_glm$theta, dist = "neg_binomial")
  glm_dev = deviance(simple_glm) 
  expect_equal(object = pkg_dev, expected = glm_dev,  tolerance = 0.001)
  ## poisson
  data$area = rlnorm(nrow(data), 0, 0.1)
  data@data$y_i = rpois(n = nrow(data), lambda = (data$area *3 * exp(data$eta)) )
  simple_glm = glm(y_i ~  factor(year) + factor(fleet_ndx), offset = log(area), data = data, family = poisson(link = log))
  pkg_dev = deviance_calc(y =  data@data$y_i, mu = fitted(simple_glm), dist = "poisson")
  glm_dev = deviance(simple_glm) 
  expect_equal(object = pkg_dev, expected = glm_dev,  tolerance = 0.001)
  ## Binomial
  data@data$area = rnorm(nrow(data))
  data@data$y_i = rbinom(n = nrow(data), size = 1, prob = plogis(data$eta))
  simple_glm = glm(y_i ~  factor(year) + factor(fleet_ndx), data = data, family = binomial(link = "logit"))
  pkg_dev = deviance_calc(y =  data@data$y_i, mu = fitted(simple_glm), dist = "binomial")
  glm_dev = deviance(simple_glm) 
  expect_equal(object = pkg_dev, expected = glm_dev,  tolerance = 0.001)
  })
