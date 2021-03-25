#'
#' @details The purpose are test configuration with a range of
#' covariates configurations.
#'

#' A non spatial GLM model with  Gamma log link response variable
test_that("config", {
  load(system.file("testdata", "config.RData",package="CPUEspatial"))
  simple_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = T, include_epsilon = T, 
                               response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, 
                               spatial_covariates = NULL, spline_catchability_covariates = NULL,
                               spline_spatial_covariates = NULL, trace_level = "none")
  
  single_catch_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                     response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx",
                                     spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                                     spline_spatial_covariates = NULL, trace_level = "none")
  
  single_catch_sptial_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                            response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                                            spatial_covariates = "habitat",  spline_catchability_covariates = NULL,
                                            spline_spatial_covariates = NULL, trace_level = "none")
  
  single_catch_sptial_model_num = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                                response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                                                spatial_covariates = "omega", spline_catchability_covariates = NULL,
                                                spline_spatial_covariates = NULL, trace_level = "none")
  
  double_catch_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                     response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = c("fleet_ndx","spatial_component"),
                                     spatial_covariates = NULL, spline_catchability_covariates = NULL,
                                     spline_spatial_covariates = NULL, trace_level = "none")

  double_catch_sptial_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                            response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = c("fleet_ndx","spatial_component"), 
                                            spatial_covariates = c("habitat", "spatial_component"), spline_catchability_covariates = NULL,
                                            spline_spatial_covariates = NULL, trace_level = "none")
  ## check spline congigurations
  spline_catch_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                     response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, 
                                     spatial_covariates = NULL,  spline_catchability_covariates = "omega",
                                     spline_spatial_covariates = NULL, trace_level = "none")
  
  spline_spatial_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                       response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, 
                                       spatial_covariates = NULL,  spline_catchability_covariates = NULL,
                                       spline_spatial_covariates = "omega", trace_level = "none")
  
  spline_both_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                    response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, 
                                    spatial_covariates = NULL, spline_catchability_covariates = "omega",
                                    spline_spatial_covariates = "epsilon", trace_level = "none")
  
  multi_spline_both_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = F, include_epsilon = F, 
                                          response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, 
                                          spatial_covariates = NULL, spline_catchability_covariates = c("omega","epsilon"),
                                          spline_spatial_covariates = c("omega","epsilon"), trace_level = "none")
  

  ## validate by looking at log likelihood is value
  expect_true(length(simple_model) == 4)
  
  })
