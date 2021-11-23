#' get_correlation
#' 
#' @details Calculate the correlation between number of samples in a projected cell with predicted spatial abundance from a model
#' with no Preferential sampling. Used as an exploratory tool to see if Preferential sampling could be a useful model addition
#' 
#' @param observed_df SpatialPointsDataFrame, which contains response variable and covariates for glmm analysis mut contain column 'area'
#' @param projection_df SpatialPointsDataFrameneeds to have the same variable names (colnames) as observed_df. Should supply variable values for all projection cells over all time steps
#' @param proj_variable_label character relating to a column name for the projected variable, usually derived from get_projection()
#' @param time_variable_label character relating to a column names needs to be an integer variable 
#' @param projection_raster_layer a RasterLayer object only required if apply_preferential_sampling = TRUE, and preference_model_type == 1. Should be the same resolution as projection_df. Used to collate sample locations. The observed_df slot should have vaules 0 = cell not in projection grid or 1 = active projection cell
#' @export
#' @importFrom raster rasterize
#' @return: projection_df with a column called predicted_y which has all the spatial components plus the time_coeffecient
get_correlation <- function(observed_df, projection_df, time_variable_label, proj_variable_label, projection_raster_layer) {
  time_variable = observed_df@data[,time_variable_label]
  time_levels = sort(unique(time_variable))
  n_t = length(unique(time_variable))
  n_z = sum(projection_df@data[,time_variable_label] == time_levels[1])
  Nij = array(NA, dim = c(n_z, n_t))
  y_hat = array(NA, dim = c(n_z, n_t))
  
  for(t in 1:n_t) {
    proj_df_subset <- subset(projection_df, subset = projection_df@data[,time_variable_label] == time_levels[t])
    data_df_subset <- subset(observed_df, subset = observed_df@data[,time_variable_label] == time_levels[t])
    data_df_subset@data$indicator = 1
    samples_count <- rasterize(data_df_subset, projection_raster_layer, field = data_df_subset$indicator, sum, na.rm = T)
    proj_count <- rasterize(proj_df_subset, projection_raster_layer, field = proj_df_subset@data[,proj_variable_label], sum, na.rm = T)
    
    Nij[,t] = samples_count$layer@data@values[projection_raster_layer@data@values == 1]
    y_hat[which(!is.na(Nij[,t])),t] = proj_count$layer@data@values[which(!is.na(Nij[,t]))]
  }
  result = list(samples_per_cell = Nij, projected_value_per_cell = y_hat, pairwise_correlation = cor(as.vector(Nij), as.vector(y_hat), use = "pairwise.complete.obs"))
  return(result)
}

