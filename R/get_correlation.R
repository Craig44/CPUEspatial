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
#' @return: a list with the following elements
#' \itemize{
#'   \item samples_rasters list of raster layers for each time-step with the accumulated observations in each projection grid
#'   \item proj_rasters list of raster layers for each time-step with the projected value
#'   \item correlation_by_time_step Pearsons correlation between complete obs of the two raster layers by time-step. non-sampled areas not included
#'   \item correlation_by_time_step_alt Pearsons correlation where non-sampled cells are set = 0 and used in correlation.
#'   \item overall_correlation Pearsons correlation for complete obs over all time steps
#' }
#' @examples
#'\dontrun{
#' you can plot the layers side by side using the following R-code
#' corr_obj = get_correlation(...)
#' plot for the first time-step
#' par(mfrow = c(2,1))
#' plot(corr_obj$samples_rasters[[1]])
#' plot(proj_rasters$samples_rasters[[1]])
#'}
get_correlation <- function(observed_df, projection_df, time_variable_label, proj_variable_label, projection_raster_layer) {
  time_variable = observed_df@data[,time_variable_label]
  time_levels = sort(unique(time_variable))
  time_levels = min(time_levels):max(time_levels)
  n_t = length(unique(time_variable))
  n_z = length(projection_raster_layer$layer@data@values)
  sample_raster = proj_raster = list()
  correlation_by_time_step = correlation_by_time_step_alt = Nij = y_hat = Nij_alt = y_hat_alt = vector();
  ## for each year
  for(t in 1:n_t) {
    proj_df_subset <- subset(projection_df, subset = projection_df@data[,time_variable_label] == time_levels[t])
    data_df_subset <- subset(observed_df, subset = observed_df@data[,time_variable_label] == time_levels[t])
    data_df_subset@data$indicator = 1
    samples_count <- rasterize(data_df_subset, projection_raster_layer, field = data_df_subset$indicator, sum, na.rm = T)
    proj_count <- rasterize(proj_df_subset, projection_raster_layer, field = proj_df_subset@data[,proj_variable_label], sum, na.rm = T)
    sample_raster[[t]] = samples_count
    proj_raster[[t]] = proj_count
    Nij = c(Nij, samples_count$layer@data@values / mean(samples_count$layer@data@values, na.rm = T))
    y_hat = c(y_hat, proj_count$layer@data@values / mean(proj_count$layer@data@values, na.rm = T))
    correlation_by_time_step[t] = cor(samples_count$layer@data@values, proj_count$layer@data@values, use = "pairwise.complete.obs")
    ## include correlation with non-sampled cells
    obs_with_zero = samples_count$layer@data@values
    obs_with_zero[is.na(obs_with_zero)] = 0
    correlation_by_time_step_alt[t] = cor(obs_with_zero, proj_count$layer@data@values, use = "pairwise.complete.obs")
  }
  result = list(samples_rasters = sample_raster, proj_rasters = proj_raster, correlation_by_time_step_alt = correlation_by_time_step_alt, correlation_by_time_step = correlation_by_time_step, overall_correlation = cor(Nij, y_hat, use = "pairwise.complete.obs"))
  return(result)
}

