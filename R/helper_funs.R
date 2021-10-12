#' evalit
#' @details evaluates a character string as if you typed it into the console
#' @param x a string to 
#' @keywords internal
#' @examples 
#' evalit("print('hello')")
#' df = data.frame(x = c(1,32,4,3), y = c(23,32,1,1))
#' var_label = "x"
#' evalit(paste0("df$",var_label)) # should return df$x
#' @return what ever the evaluation string defines
evalit = function (x, ...)  {
  eval(parse(text = x), ...)
}

#' check_dims
#' @details check a dimension is consistent between two matricies
#' @param m1 a matrix
#' @param m2 a matrix
#' @param row boolean check row if false check columns
#' @keywords internal
#' 
#' @return list with bool and message
check_dims = function (m1, m2, row)  {
  passed_check = TRUE
  msg = "success"
  if(row) {
    if(nrow(m1) != nrow(m2)) {
      passed_check = FALSE
      msg = paste0("the first matrix had '", nrow(m1), "' rows and the second had, '", nrow(m2), "'.")
    } 
  } else {
    if(ncol(m1) != ncol(m2)) {
      passed_check = FALSE
      msg = paste0("the first matrix had '", ncol(m1), "' columns and the second had, '", ncol(m2), "'.")
    } 
  }
  return(list(passed_check = passed_check, msg = msg))
}

#' check_gradients
#' @details checks a TMB object for fixed effect parameters that have 0 gradients, suggesting they don't contribute the log likelihood. And you should 
#' look into these parameters.
#' @param obj A TMB list that has been built by MakeAdFun
#' @export
#' @return character of good news or labels of problem parameters
#' 
check_gradients = function(obj) {
  if(sum(obj$gr() == 0)) {
    return("no Non-zero gradients, give estimating a go")
  } else {
    return(names(obj$par)[which(obj$gr() == 0)])
  }
  return(NULL)
}



#' inverse_link
#' @details applies the inverse of the link function
#' log_link                 = 0
#' logit_link               = 1
#' probit_link              = 2
#' inverse_link             = 3
#' identity_link            = 4
#' @param x a vector of values to transform
#' @param link_fun integer based on the above
#' @keywords internal
#' @return 
inverse_link = function (x, link_fun) {
  return_val = x; ## default identity link
  if(link_fun == 0) {
    return_val = exp(x)
  } else if(link_fun == 1) {
    return_val = plogis(x)
  } else if(link_fun == 2) {
    return_val = pnorm(x)
  } else if(link_fun == 3) {
    return_val = pnorm(x)
  } else if(link_fun == 4) {
    # do nothing identiy link
  }
  return(return_val);
}

#' logit_general 
#' @details bounds X which is between [lb,ub] to -inf -> inf based on the logit transformation
#' @param X scalar range [lb,ub]
#' @param ub upper bound for X
#' @param lb lower bound for X
#' @return Y to be between [-inf, inf]
logit_general = function(X, lb, ub) {
  X1 = (X - lb) / (ub - lb)
  log(X1/(1 - X1))
}

#' rmvnorm_prec 
#' @details
#' simualte parameters from the joint precision matrix derived from a TMB objects
#' @param mu vector of MLE both fixed and random effect parameters
#' @param prec precision matrix, derived from sdreport(obj, getJointPrecision = T)
#' @param n.sims integer number of simulations
#' @param random_seed integer seed
#' @importFrom stats rnorm
#' @importFrom Matrix solve Cholesky
#' @export
#' @return matrix of simulated parameter values with dimensions [n.sims x length(mu)]
rmvnorm_prec <- function(mu, prec, n.sims, random_seed = trunc(runif(1, 1, 1e5))) {
  set.seed( random_seed )
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L = Matrix::Cholesky(prec, super=TRUE)
  z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z = as.matrix(z)
  return(mu + z)
}

#' get_tmb_fixed_effects
#' @details returns MLE estimates for fixed effects
#' @param obj An optimised list that has been build by MakeAdFun
#' 
get_tmb_fixed_effects = function(obj) {
  if (length(obj$env$random) == 0) {
    return(obj$env$last.par.best)
  }
  return(obj$env$last.par.best[-obj$env$random])
}


#' check_tmb_convergence
#' @author C.Marsh
#' @details use TMB object to check gradients and if the hessian is definite positive
#' @param obj An optimised list that has been build by MakeAdFun
#' @param delta Gradient threshold for defining a converged model
#' @importFrom TMB sdreport 
#' @importFrom stats optimHess 
#' 
check_tmb_convergence = function(obj, delta = 0.001) {
  best_fixed_eff_pars = get_tmb_fixed_effects(obj)
  grads = tryCatch( expr = obj$gr(best_fixed_eff_pars));
  if (inherits(grads, "error"))
    return("Could not get gradients with an error. generally happens from unconverged models please check")
  if (any(is.na(grads)))
    return("Found gradients with NaN, conclusion = unconverged model")
  labs = names(obj$par)
  if(max(abs(grads)) > delta)
    return(paste0("param with label = ", labs[which.max(abs(grads) > delta)], " has gradient > ", delta, " |grad| = ", round(max(abs(grads)), 5)))
  hess = optimHess(best_fixed_eff_pars, fn = obj$fn, gr = obj$gr)
  sdr = sdreport(obj, getJointPrecision = TRUE)
  if (is.character(try(chol(sdr$cov.fixed), silent = TRUE)))
    return("Covariance of fixed effect no positive Definitive")
  if ("jointPrecision" %in% names(sdr)) {
    if (is.character(try(chol(sdr$jointPrecision), silent = TRUE)))
      return("Joint Precision not positive Definitive")
  }
  return("No evidence of non convergence")
}
  
#' TMB helper function
#' @author C.Marsh
#' @description this function returns a list of factors used in the map argument of the MakeADFun function
#' @param par_list a named list that you give to the par argument in the MakeADFun
#' @param pars_to_exclude a vector of strings with names of parmeters you want to FIX in the objective object.
#' @param vec_pars_to_adjust a vector string of parameter labels that we want to exclude certain elements.
#' @param vec_elements_to_exclude a named list (names same as par_list) with number of elements = length(vec_pars_to_adjust). each list element 
#' contains a vector of elements that we want to exclude from estimation.
#' @export
#' @return a list of factors used in the MakeADFun function
fix_pars = function(par_list, pars_to_exclude, vec_pars_to_adjust = NULL, vec_elements_to_exclude = NULL) {
  if (!any(pars_to_exclude %in% names(par_list))) {
    stop(paste0("The parameters ", paste(pars_to_exclude[!pars_to_exclude %in% names(par_list)],collapse = " ")," in exclusion parameters could not be found in the 'par_list', please sort this out"))
  }
  pars = names(par_list)
  mapped_pars = list();
  tailor_vectors = FALSE
  if (!is.null(vec_pars_to_adjust)) {
    tailor_vectors = TRUE
    if (!all(vec_pars_to_adjust %in% pars_to_exclude))
      stop("parmaeters noted in vec_pars_to_adjust, need to also be in pars_to_exclude")
  }
  param_factor = 1;
  for(i in 1:length(pars)) {
    if (pars[i] %in% pars_to_exclude) {
      params_in_this_par = par_list[[pars[i]]];
      if (tailor_vectors & (pars[i] %in% vec_pars_to_adjust)) {
        include_element_index = c(1:length(params_in_this_par))[-vec_elements_to_exclude[[which(pars[i] %in% names(vec_elements_to_exclude))]]]
        params_vals = factor(rep(NA, length(params_in_this_par)), levels = factor(param_factor:(param_factor + length(include_element_index) - 1)))
        params_vals[include_element_index] = factor(param_factor:(param_factor + length(include_element_index) - 1))#, levels = factor(include_element_index))
        param_factor = param_factor + length(include_element_index)
        mapped_pars[[pars[i]]] = params_vals;
      } else {
        mapped_pars[[pars[i]]] = rep(factor(NA),length(params_in_this_par));
      }
    } else {
      params_in_this_par = par_list[[pars[i]]];
      params_vals = factor(param_factor:(param_factor + length(params_in_this_par) - 1))
      param_factor = param_factor + length(params_in_this_par)
      mapped_pars[[pars[i]]] = params_vals
    }
  }
  return(mapped_pars);
}
#' eigen_decomp_covariance Do an eigen decomposition to look at poorly estimated parameters from MLE fit
#' @param covariance_matrix symetric covariance matrix
#' @param param_labels vector of param labels (optional)
#' @param delta a cut off value for 'poorly' defined parameters.
#' @export
#' @return: data frame of eiegen values for the matrix and index of good and bad parameters based on delta
#'
eigen_decomp_covariance = function(covariance_matrix, param_labels = NULL, delta = .Machine$double.eps) {
  ## check covariance is invertable
  if (!isSymmetric(covariance_matrix))
    stop("covariance matrix is not symetric something is wrong here.")
  ## check positive semi defintie matrix
  if(class(try(solve(covariance_matrix),silent=T)) != "matrix")
    stop("covariance not invertible")
  ## calculate hessian
  hess = solve(covariance_matrix)
  ## eigen decomposition
  Eig = eigen(hess)
  WhichBad = which(Eig$values < sqrt(delta))
  df = NULL;
  if (is.null(param_labels)) 
    param_labels = as.character(1:ncol(covariance_matrix))
  
   
  if (length(WhichBad) == 0) {    
    message( "All parameters are identifiable" )
  } else {
    # Check for parameters
    RowMax = apply( Eig$vectors[,WhichBad,drop=FALSE], MARGIN=1, FUN=function(vec){max(abs(vec))} )
    df = data.frame("Param"=param_labels, "eigenvalues", Eig$values ,"Param_check"=ifelse(RowMax>0.1, "Bad","OK"))
  }
  return(df)
}

#' cumulant_fun
#' @details This is the cumulative function for the exponential families table 5.1 of GLM with examples in R DUnn 2018
#' @param dist character values options availble {gaussian, gamma,binomial, neg_binomial, inverse_gaussian, poisson}
#' @param theta canonical form of the distirbution
#' @keywords internal
#' @return returns the evaluation of the cumulative function given a theta
cumulant_fun = function(theta, dist) {
  result = NULL
  if(dist == "gaussian") {
    result = theta^2 /2
  } else if (dist == "gamma") {
    result = -log(-theta)
  } else if (dist == "binomial") {
    result = exp(theta) / (1 + exp(theta))
  } else if (dist == "neg_binomial") {
    result = - log(1 - exp(theta))
  } else if (dist == "poisson") {
    result = exp(theta)
  } else if (dist == "inverse_gaussian") {
    result = -sqrt(-2*theta)
  } else {
    stop("unknown distribution")
  }
  return(result)
}
#' canonical_fun
#' @details returns the canonical form (often denoted as theta in the literature) given mu
#' @param mu fitted values from the model
#' @param dist character values options availble {gaussian, gamma,binomial, neg_binomial, inverse_gaussian, poisson}
#' @param k the k parameter needed for the negative binomial evaluation
#' @keywords internal
#' @return canonical given fitted value
canonical_fun = function(mu, dist, k = NULL) {
  theta = NULL
  if(dist == "gaussian") {
    theta = mu
  } else if (dist == "gamma") {
    theta = - 1 / mu
  } else if (dist == "binomial") {
    theta = log(mu / (1 - mu))
  } else if (dist == "neg_binomial") {
    theta = log(mu / (mu + k))
  } else if (dist == "poisson") {
    theta = log(mu)
  } else if (dist == "inverse_gaussian") {
    theta = - 1 / (2*mu^2)
  } else {
    stop("unknown distribution")
  }
  return(theta)
}

#' deviance_calc
#' @details returns the residual deviance based on the exponential family. That is 2 time fitted compared with the saturated model
#' @param y response variable values
#' @param mu fitted values from the model
#' @param dist character values options availble {gaussian, gamma,binomial, neg_binomial, inverse_gaussian, poisson}
#' @param k the k parameter needed for the negative binomial evaluation
#' @importFrom MASS neg.bin

#' @export
#' @return residual deviance
deviance_calc = function(y, mu, dist, phi = NULL) {
  result = NULL
  if(dist == "poisson") {
    result = sum(poisson()$dev.resids(y, mu, wt = rep(1, length(y))))
  } else if(dist == "binomial") {
    result = sum(binomial()$dev.resid(y, mu, wt = rep(1, length(y))))
  } else if(dist == "neg_binomial") {
    result = sum(neg.bin(theta = phi)$dev.resids(y, mu, wt = rep(1, length(y))))
  } else {
    ## saturated t(y,y)
    theta_s = canonical_fun(y, dist = dist, k = k)
    ## fitted t(y,mu)
    theta_fit = canonical_fun(mu, dist = dist, k = k)
    result = sum(2 * ((y * theta_s - cumulant_fun(theta_s, dist))
                      - (y * theta_fit - cumulant_fun(theta_fit, dist))))
  }

  

  return(result)
}
