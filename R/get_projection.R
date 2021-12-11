#' get_projection
#' 
#' @details used to get the projected fitted values from a given model.
#' 
#' @param obj TMB object which has been optimised and checked for convergence
#' @param data the tmb_data element that is returned from configure_obj() function
#' @param projection_df spatial data frame, should be the same one given to configure_obj()
#' @param time_variable_label character
#' @export
#' @return: projection_df with a column called predicted_y which has all the spatial components plus the time_coeffecient
get_projection <- function(obj, data, projection_df, time_variable_label) {
  model_type = obj$env$data$model
  time_variable = projection_df@data[,time_variable_label]
  time_levels = sort(unique(time_variable))
  MLE_report = obj$report(obj$env$last.par.best)
  projection_df@data$predicted_y = NA
  for(t in 1:data$n_t) {
    time_ndx = projection_df@data[,time_variable_label] == time_levels[t]
    sub_proj = subset(projection_df, subset = time_ndx)
    spatial_prediction = NULL
    spatial_betas = MLE_report$spatial_betas;
    if(length(spatial_betas) == 1)
      spatial_betas = matrix(spatial_betas, nrow = 1)
    if(model_type == "SpatialTemporalCPUE") {
      spatial_prediction =  MLE_report$time_betas[t] + MLE_report$omega_proj + data$X_spatial_proj_zpt[,,t] %*% spatial_betas + (data$Proj %*% MLE_report$epsilon_input[,t])  / MLE_report$tau_epsilon
    } else if (model_type == "SpatialTemporalCPUENN") {
      spatial_prediction =  MLE_report$time_betas[t] + MLE_report$omega_proj + data$X_spatial_proj_zpt[,,t] %*% spatial_betas + MLE_report$epsilon_input[data$index_proj_vertex + 1, t]
    } else if (model_type == "SpatialTemporalCPUEVAST") {
      spatial_prediction =  MLE_report$time_betas[t] + MLE_report$omega_proj + data$X_spatial_proj_zpt[,,t] %*% spatial_betas + MLE_report$epsilon_v[data$index_proj_vertex + 1, t]
    }
    projection_df@data$predicted_y[time_ndx] = inverse_link(spatial_prediction, data$link)
  }
  return(projection_df)
}
