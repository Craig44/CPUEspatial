#' calculate_influence
#' 
#' @details 
#' given an optimsied model that has been output from configure_obj() function, we will generate a bunch of influence statistics following the work from Nokome's influence plot metrics.
#' Bentley et.al. 2011 Influence plots and metrics: tools for better understanding fisheries catch-per-unit-effort standardizations
#' @param conf_obj object for which the obj element has been optimsied
#' @param data_df data.frame of original data set
#' @return list objects that summarise of influence statistics that can then be passed to a plotting function
#' @importFrom stats formula get_all_vars
#' @importFrom TMB sdreport
#' @importFrom dplyr group_by mutate summarise
#' @export
calculate_influence = function(conf_obj, data_df) {
  ## do some checks and balances
  if(!conf_obj$tmb_data$link %in% c(0))
    warning("This is functionality is made for multiplicative models. Usually associated with log-link models, could be the identity link if you have logged the response variable previously. See the original Bentley 2011 paper for more information")
  
  ## get standard errors
  sd_rep = sdreport(conf_obj$obj)
  sd_rep_labs = names(sd_rep$value)
  ## get the terms
  y_lab = eval(conf_obj$Call$func_call$response_variable_label)
  time_lab = eval(conf_obj$Call$func_call$time_variable_label)
  catch_terms = eval(conf_obj$Call$func_call$catchability_covariates)
  spatial_terms = eval( conf_obj$Call$func_call$spatial_covariates)
  all_terms = c(catch_terms, spatial_terms)
  ## get term type
  catchability_type = spatial_type = vector()
  for(i in 1:length(spatial_terms))
    spatial_type[i] = class(data_df@data[, spatial_terms[i]])
  for(i in 1:length(catch_terms))
    catchability_type[i] = class(data_df@data[, catch_terms[i]])
  ##  Create model matricies for the two components
  catchability_mod_mat = conf_obj$tmb_data$model_matrix
  spatial_mod_mat = conf_obj$tmb_data$X_spatial_ipt[,,1]
  catchability_coef_labs = colnames(catchability_mod_mat)
  spatial_coef_labs = colnames(spatial_mod_mat)

  ## create a data frame that has each variable for each each term
  catchability_lab_df = get_all_vars(formula(conf_obj$Call$catchability), data = data_df@data)
  catchability_lab_df$time = data_df@data[, time_lab]
  spatial_lab_df = get_all_vars(formula(conf_obj$Call$spatial) , data = data_df@data)
  spatial_lab_df$time = data_df@data[, time_lab]
  
  ## get coeffecients for each term and estimated standand error
  terms_ls = list()
  for(i in 1:length(catch_terms)) {
    coeff_ndx = grepl(catchability_coef_labs, pattern = catch_terms[i])
    labs = substring(catchability_coef_labs[coeff_ndx], first = nchar(catch_terms[i]) + 1)
    if(catchability_type[i] == "factor") { #include intercept
      coeff_ndx[1] = TRUE
      MLE = sd_rep$value[sd_rep_labs %in% "betas_w_intercept"][coeff_ndx]
      se = sd_rep$sd[sd_rep_labs %in% "betas_w_intercept"][coeff_ndx]
      ## find levels
      data_labs = unique(as.character(catchability_lab_df[,catch_terms[i]]))
      intercept = data_labs[!data_labs %in% labs]
      labs = c(intercept, labs)
      labs = factor(labs, levels = levels(data_df@data[, catch_terms[i]]))
      terms_ls[[catch_terms[i]]] = data.frame(labs = labs, MLE = MLE, SE = se, lower = MLE - 2*se, upper = MLE + 2 * se)
      ## attach coeffecients catchability_lab_df
      orig_colnames = colnames(catchability_lab_df)
      data_terms = get(catch_terms[i], catchability_lab_df)
      ndx = match(data_terms, catch_terms_ls[[catch_terms[i]]]$labs)
      coeff_terms = catch_terms_ls[[catch_terms[i]]]$MLE[ndx]
      catchability_lab_df = cbind(catchability_lab_df, coeff_terms)
      colnames(catchability_lab_df) = c(orig_colnames, paste0(catch_terms[i], ".fit"))
    } else {
      ## just slope parameter
      MLE = sd_rep$value[sd_rep_labs %in% "betas"][coeff_ndx]
      se = sd_rep$sd[sd_rep_labs %in% "betas"][coeff_ndx]
      terms_ls[[catch_terms[i]]] = data.frame(labs = labs, MLE = MLE, SE = se, lower = MLE - 2*se, upper = MLE + 2 * se)
      orig_colnames = colnames(catchability_lab_df)
      data_terms = get(catch_terms[i], catchability_lab_df)
      catchability_lab_df = cbind(catchability_lab_df, MLE * data_terms)
      colnames(catchability_lab_df) = c(orig_colnames, paste0(catch_terms[i], ".fit"))
    }
  }
  ## spatial factors are treated differently
  for(i in 1:length(spatial_terms)) {
    coeff_ndx = grepl(spatial_coef_labs, pattern = spatial_terms[i])
    labs = substring(spatial_coef_labs[coeff_ndx], first = nchar(spatial_terms[i]) + 1)
    if(spatial_type[i] == "factor") { #include intercept
      MLE = sd_rep$value[sd_rep_labs %in% "spatial_betas"][coeff_ndx]
      se = sd_rep$sd[sd_rep_labs %in% "spatial_betas"][coeff_ndx]
      labs = factor(labs)
      ## find levels
      terms_ls[[spatial_terms[i]]] = data.frame(labs = labs, MLE = MLE, SE = se, lower = MLE - 2*se, upper = MLE + 2 * se)
      ## attach coeffecients catchability_lab_df
      orig_colnames = colnames(spatial_lab_df)
      data_terms = get(spatial_terms[i], spatial_lab_df)
      ndx = match(data_terms, spatial_terms_ls[[spatial_terms[i]]]$labs)
      coeff_terms = spatial_terms_ls[[spatial_terms[i]]]$MLE[ndx]
      spatial_lab_df = cbind(spatial_lab_df, coeff_terms)
      colnames(spatial_lab_df) = c(orig_colnames, paste0(spatial_terms[i], ".fit"))
    } else {
      ## just slope parameter
      MLE = sd_rep$value[sd_rep_labs %in% "spatial_betas"][coeff_ndx]
      se = sd_rep$sd[sd_rep_labs %in% "spatial_betas"][coeff_ndx]
      terms_ls[[spatial_terms[i]]] = data.frame(labs = labs, MLE = MLE, SE = se, lower = MLE - 2*se, upper = MLE + 2 * se)
      orig_colnames = colnames(spatial_lab_df)
      data_terms = get(spatial_terms[i], spatial_lab_df)
      spatial_lab_df = cbind(spatial_lab_df, MLE * data_terms)
      colnames(spatial_lab_df) = c(orig_colnames, paste0(spatial_terms[i], ".fit"))
    }
  }
  ## get distribution of observations
  ## number of observatons over time among the levels in a variable
  all_distr = list()
  for(i in 1:length(catch_terms)) {
    distr = suppressMessages(catchability_lab_df %>% 
      group_by(get(catch_terms[i]), time) %>%
      summarise(n = length(time)))
    distr = distr %>% group_by(time) %>% mutate(n_total = sum(n))
    distr$n_prop = distr$n / distr$n_total
    distr$time = factor(distr$time, levels = min(distr$time):max(distr$time))
    colnames(distr) = c("term", conf_obj$Call$func_call$time_variable_label, "n","n_total","n_prop")
    all_distr[[catch_terms[i]]] = distr
  }
  for(i in 1:length(spatial_terms)) {
    distr = suppressMessages(spatial_lab_df %>% 
      group_by(get(spatial_terms[i]), time) %>%
      summarise(n = length(time)))
    distr = distr %>% group_by(time) %>% mutate(n_total = sum(n))
    distr$n_prop = distr$n / distr$n_total
    distr$time = factor(distr$time, levels = min(distr$time):max(distr$time))
    colnames(distr) = c("term", conf_obj$Call$func_call$time_variable_label, "n","n_total","n_prop")
    all_distr[[spatial_terms[i]]] = distr
  }
  ## now influence metrics
  overall_influence = vector()
  influence_df = NULL;
  counter = 1;
  for(i in 1:length(catch_terms)) {
    mean_term = mean(get(paste0(catch_terms[i], ".fit"), catchability_lab_df), i = 1)
    influence_ = suppressMessages(catchability_lab_df %>%
      group_by(time) %>%
      summarise(delta = mean(get(paste0(catch_terms[i], ".fit")) - mean_term)))
    if(i == 1) {
      influence_df = influence_
    } else {
      influence_df = cbind(influence_df, influence_$delta)
    }
    overall_influence[counter] = exp(mean(abs(influence_$delta))) - 1
    counter = counter + 1
  }
  for(i in 1:length(spatial_terms)) {
    mean_term = mean(get(paste0(spatial_terms[i], ".fit"), spatial_lab_df), i = 1)
    influence_ = suppressMessages(spatial_lab_df %>%
      group_by(time) %>%
      summarise(delta = mean(get(paste0(spatial_terms[i], ".fit")) - mean_term)))
    if(is.null(influence_df)) {
      influence_df = influence_
    } else {
      influence_df = cbind(influence_df, influence_$delta)
    }
    overall_influence[counter] = exp(mean(abs(influence_$delta))) - 1
    counter = counter + 1
  }
  names(overall_influence) = c(catch_terms, spatial_terms)
  colnames(influence_df) = c("time", catch_terms, spatial_terms)
  ## overall influence
  
  return(list(influence_df = influence_df, all_distr = all_distr, terms_ls = terms_ls, catchability_lab_df = catchability_lab_df, spatial_lab_df = spatial_lab_df, overall_influence = overall_influence))
}


