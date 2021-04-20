#'
#' @details The purpose are test configuration with a range of
#' covariates configurations.
#'

#' A non spatial GLM model with  Gamma log link response variable
test_that("config", {
  load(system.file("testdata", "config.RData",package="CPUEspatial"))
  simple_model = configure_obj(observed_df = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = T, include_epsilon = T, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, 
                               spatial_covariates = NULL, spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none")
  
  single_catch_model = configure_obj(observed_df = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                     response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx",
                                     spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                                     spline_spatial_covariates = NULL, trace_level = "none")
  
  single_catch_sptial_model = configure_obj(observed_df = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                            response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                                            spatial_covariates = "habitat",  spline_catchability_covariates = NULL,
                                            spline_spatial_covariates = NULL, trace_level = "none")
  
  single_catch_sptial_model_num = configure_obj(observed_df = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                                response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                                                spatial_covariates = "omega", spline_catchability_covariates = NULL,
                                                spline_spatial_covariates = NULL, trace_level = "none")
  
  double_catch_model = configure_obj(observed_df = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                     response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = c("fleet_ndx","spatial_component"),
                                     spatial_covariates = NULL, spline_catchability_covariates = NULL,
                                     spline_spatial_covariates = NULL, trace_level = "none")

  double_catch_sptial_model = configure_obj(observed_df = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                            response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = c("fleet_ndx","spatial_component"), 
                                            spatial_covariates = c("habitat", "spatial_component"), spline_catchability_covariates = NULL,
                                            spline_spatial_covariates = NULL, trace_level = "none")
  ## check spline congigurations
  spline_catch_model = configure_obj(observed_df = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                     response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, 
                                     spatial_covariates = NULL,  spline_catchability_covariates = "omega",
                                     spline_spatial_covariates = NULL, trace_level = "none")
  
  spline_spatial_model = configure_obj(observed_df = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                       response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, 
                                       spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                                       spline_spatial_covariates = "omega", trace_level = "none")
  
  spline_both_model = configure_obj(observed_df = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                    response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, 
                                    spatial_covariates = NULL, spline_catchability_covariates = "omega",
                                    spline_spatial_covariates = "epsilon", trace_level = "none")
  
  multi_spline_both_model = configure_obj(observed_df = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                          response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, 
                                          spatial_covariates = NULL, spline_catchability_covariates = c("omega","epsilon"),
                                          spline_spatial_covariates = c("omega","epsilon"), trace_level = "none")
  

  ## validate by looking at log likelihood is value
  expect_true(length(simple_model) == 4)
  
  })


#' A spatial GLM with preferential sampling, checks the configure_obj() function works for expected inputs.
test_that("config_with_pref", {
  load(system.file("testdata", "config.RData",package="CPUEspatial"))

  eps_model = configure_obj(observed_df = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = T, 
                                       response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, 
                                       spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                                       spline_spatial_covariates = NULL, apply_preferential_sampling = TRUE, preference_model_type = 0, pref_hyper_distribution = 0,
                                       trace_level = "none")
  
  eps_model_time_vary = configure_obj(observed_df = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = T, 
                                      response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, 
                                      spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                                      spline_spatial_covariates = NULL, apply_preferential_sampling = TRUE, preference_model_type = 0, pref_hyper_distribution = 1,
                                      logit_pref_hyper_prior_vals = c(0,0.5), trace_level = "none")  
  
  eps_model_time_vary_2 = configure_obj(observed_df = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = T, 
                                      response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, 
                                      spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                                      spline_spatial_covariates = NULL, apply_preferential_sampling = TRUE, preference_model_type = 0, pref_hyper_distribution = 2,
                                      logit_pref_hyper_prior_vals = c(0,0.5), trace_level = "high")  
  eps_model_time_vary_3 = configure_obj(observed_df = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = T, 
                                      response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, 
                                      spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                                      spline_spatial_covariates = NULL, apply_preferential_sampling = TRUE, preference_model_type = 0, pref_hyper_distribution = 3,
                                      logit_pref_hyper_prior_vals = c(0,0.5), trace_level = "none")  
  
  
  
  ## validate by looking at log likelihood is value
  expect_true(length(eps_model_time_vary_3) == 4)
  
})