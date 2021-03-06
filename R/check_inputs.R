#' check_inputs
#' @details a wrapper function to check data is compatible with parameters, so we don't get a segment fault from the TMB MakeAdFun
#' @param tmb_data a list of data
#' @param tmb_params a list of parameters
#' @keywords internal
#' @return list with two elements a boolean called result (True if everything is okay and False if flagged warnigns) and string vector of errors if result == FALSE
check_inputs = function (tmb_data, tmb_params)  {
  error_msgs = c();
  ## to fill this out with more sanity tests
  ## linear covariates
  if(sum(tmb_data$year_ndx_for_each_obs >= 0) != tmb_data$n_i)
    error_msgs = c(error_msgs, "year_ndx_for_each_obs has to have the same number of indices as there are data values, this will cause a crash.")

  test = tryCatch(expr = tmb_data$spline_spatial_model_matrix_ipt[,,1] %*% tmb_params$gammas_spatial, error = function(e) {e})
  if(inherits(test, "error")) 
    error_msgs = c(error_msgs, "spline_spatial_model_matrix_ipt not compatible with gammas_spatial")

  
  ## spline dimensions
  test = tryCatch(expr = tmb_data$spline_spatial_model_matrix_ipt[,,1] %*% tmb_params$gammas_spatial, error = function(e) {e})
  if(inherits(test, "error")) 
    error_msgs = c(error_msgs, "spline_spatial_model_matrix_ipt not compatible with gammas_spatial")
  
  test = tryCatch(expr = tmb_data$spline_spatial_model_matrix_proj_zpt[,,1] %*% tmb_params$gammas_spatial, error = function(e) {e})
  if(inherits(test, "error")) 
    error_msgs = c(error_msgs, "spline_habitat_model_matrix_proj_zpt not compatible with gammas_spatial")
  
  test = tryCatch(expr = tmb_data$designMatrixForReport_spatial %*% tmb_params$gammas_spatial, error = function(e) {e})
  if(inherits(test, "error")) 
    error_msgs = c(error_msgs, "designMatrixForReport_spatial not compatible with gammas_spatial")
  
  test = tryCatch(expr = tmb_data$spline_model_matrix %*% tmb_params$gammas, error = function(e) {e})
  if(inherits(test, "error")) 
    error_msgs = c(error_msgs, "spline_model_matrix not compatible with gammas")
  
  test = tryCatch(expr = tmb_data$designMatrixForReport %*% tmb_params$gammas, error = function(e) {e})
  if(inherits(test, "error")) 
    error_msgs = c(error_msgs, "designMatrixForReport not compatible with gammas")
  
  if(any(tmb_data$t_i > tmb_data$n_t)) 
    error_msgs = c(error_msgs, "t_i has a value equal to or larger than n_t. This will cause an out of memory access fault as C++ starts indexing from 0!!!")
  
  
  if(length(error_msgs) != 0) 
    return(list(result = FALSE, errors = error_msgs))
  
  return(list(result = TRUE, errors = error_msgs))
}
