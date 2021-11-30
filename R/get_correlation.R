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
#'   \item correlation_by_time_step Pearsons correlation between complete obs of the two raster layers by time-step
#'   \item overall_correlation Pearsons correlation for complete obs over all time steps
#' }
#' @details the overall correlation uses standardised values for a given year. if n = obs in a year and p = projected values in a year then this returns cor(n / sum(n), p / sum(p))
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
  n_t = length(unique(time_variable))
  n_z = length(projection_raster_layer$layer@data@values)
  sample_raster = proj_raster = list()
  correlation_by_time_step =  nij_stand = y_hat_stand = vector();
  ## for each year
  for(t in 1:n_t) {
    proj_df_subset <- subset(projection_df, subset = projection_df@data[,time_variable_label] == time_levels[t])
    data_df_subset <- subset(observed_df, subset = observed_df@data[,time_variable_label] == time_levels[t])
    data_df_subset@data$indicator = 1
    samples_count <- rasterize(data_df_subset, projection_raster_layer, field = data_df_subset$indicator, sum, na.rm = T)
    proj_count <- rasterize(proj_df_subset, projection_raster_layer, field = proj_df_subset@data[,proj_variable_label], sum, na.rm = T)
    sample_raster[[t]] = samples_count
    proj_raster[[t]] = proj_count
    nij_stand = c(nij_stand, samples_count$layer@data@values / sum(samples_count$layer@data@values, na.rm  =T))
    ## standardised in each year
    ## otherwise get correlations that
    y_hat_stand = c(y_hat_stand, proj_count$layer@data@values/ sum(proj_count$layer@data@values, na.rm = T))
    correlation_by_time_step[t] = cor(samples_count$layer@data@values, proj_count$layer@data@values, use = "pairwise.complete.obs")
  }
  result = list(samples_rasters = sample_raster, proj_rasters = proj_raster, correlation_by_time_step = correlation_by_time_step, overall_correlation = cor(nij_stand, y_hat_stand, use = "pairwise.complete.obs"))
  return(result)
}

