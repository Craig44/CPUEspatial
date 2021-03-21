#'
#' A script to validate functionality with GLM and GAM
#' gamma response variable with discrete spatail area plus numeric linear spatial variable fleet and year effects
#' with time-invariante GF where observation locations are tied to high abundance
#'
library(TMB)
library(CPUEspatial)
library(INLA)
library(RandomFields)
library(RANN)
library(ggplot2)
library(gridExtra)
library(mgcv) 
library(spatstat)
library(DHARMa)
library(dplyr)

source(file.path("inst","examples","canonical.index.R"))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
set.seed(123)
# spatial domain even though we are not doing a spatial model
# this model requires it.
grid_dim = c("x"=100, "y"=100)
max.edge = c(4, 15)
cutoff = c(6)
# marginal variance of the spatial-temporal field 
epsilon_sigma2 <- 1  # marginal variance of the eta field
epsilon_range <- 30     # range-decorrelation parameter
nu = 1             # smoothness parameter  
#omega_model = RMwhittle(nu = nu, var = omega_sigma2, scale = omega_range )
epsilon_model = RMmatern(nu = nu, var = epsilon_sigma2, scale = epsilon_range / sqrt(8))

# sample size
#n = 1000 # per year n is generated from a LGCP so is a random variable
n_years = 10
year_coef = c(rnorm(n_years / 2, 0.5, 0.6),rnorm(n_years / 2, 0, 0.6))
plot(1:n_years, year_coef, type = "o")
normalised_year_coef = year_coef - mean(year_coef)
shape = 50 ### for gamma response variable

region_coef = rnorm(4,0, 0.6) # region coeffecients
region_lab = c("A","B","C","D")
normalised_region_coef = region_coef - mean(region_coef)
#
fleet_coef = c(-0.7, 0.31, -0.2) # simple vessel coeffecients
normalised_fleet_coef = fleet_coef - mean(fleet_coef)
## change fleet over time
prob_fleet1_per_year = seq(from = 0.6, to = 0.2, length.out = n_years)
prob_fleet2_per_year = c(rep(0, n_years - (n_years - 6)), seq(from = 0.2, to = 0.4, length.out = n_years - 6))
prob_fleet3_per_year = 1 - (prob_fleet1_per_year + prob_fleet2_per_year)
intercept = mean(year_coef) + fleet_coef[1] + mean(region_coef)

## projection space
predxseq = seq(0, grid_dim["x"],length.out = 50)
unique(round(diff(predxseq),4))
predyseq = seq(0, grid_dim["y"],length.out = 50)
unique(round(diff(predyseq),4))
proj_df = expand.grid(predxseq, predyseq)
colnames(proj_df) = c("x","y")
proj_df$region = NA
proj_df$region[proj_df$x < 50 & proj_df$y > 50] = 1
proj_df$region[proj_df$x > 50 & proj_df$y > 50] = 2
proj_df$region[proj_df$x < 50 & proj_df$y < 50] = 3
proj_df$region[proj_df$x > 50 & proj_df$y < 50] = 4


## area cells
plot(1, 1, xlim = c(0, grid_dim["x"]), ylim = c(0, grid_dim["y"]), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i")
abline(v = c(seq(0,grid_dim["x"], by = 50)), lwd = 3)
abline(h = c(seq(0,grid_dim["y"], by = 50)), lwd = 3)
text(x = 25, y  = 75, labels = "A", cex = 3)
text(x = 75, y  = 75, labels = "B", cex = 3)
text(x = 25, y  = 25, labels = "C", cex = 3)
text(x = 75, y  = 25, labels = "D", cex = 3)
#year_samples = sample(1:n_years, size = n * n_years, replace = T,  prob = rep(1,n_years) / n_years)


sampData = fleet_ndx = NULL;
full_proj_df = NULL
for(i in 1:n_years) {
  ## Simulate sample locations based on population spatial distribution
  proj_df$epsilon <- RFsimulate(epsilon_model, x=as.matrix(proj_df[,c("x","y")]), exactness=TRUE)$variable1
  ggplot(proj_df, aes(x = x, y = y, fill = epsilon)) +
    geom_tile()
  
  ## Point process
  intercept_pp = -2.3 # why an intercept model doesn't work, number of observations
  # would be tied to abundance.
  proj_df$log_lambda = intercept_pp + proj_df$epsilon + normalised_region_coef[proj_df$region]
  lam_Im = im(t(matrix(exp(proj_df$log_lambda), nrow=length(predxseq), ncol=length(predyseq))), xcol = predxseq, yrow = predyseq)
  integral(lam_Im)
  
  #ggplot(proj_df, aes(x = x, y = y, fill = omega)) + 
  #  geom_tile()
  # samples_this_year = sum(year_samples == i)
  ## simualte from a point process
  ppp <- rpoispp(lam_Im, win = win)
  
  fleet_ndx = c(fleet_ndx, sample(1:length(fleet_coef), size = ppp$n, replace = T,  prob = c(prob_fleet1_per_year[i], prob_fleet2_per_year[i],prob_fleet3_per_year[i])))
  ##
  ## match the closest depth variable
  ndx = nn2(proj_df[ ,c("x","y")], cbind(ppp$x, ppp$y), k = 1)$nn.idx
  sampData = rbind(sampData, data.frame(year = i, x = ppp$x, y = ppp$y, epsilon = proj_df$epsilon[ndx]))
  
  ##
  proj_df$year = i
  full_proj_df = rbind(full_proj_df, proj_df)
}
n = nrow(sampData)
sampData$region = NA
sampData$region[sampData$x < 50 & sampData$y > 50] = 1
sampData$region[sampData$x > 50 & sampData$y > 50] = 2
sampData$region[sampData$x < 50 & sampData$y < 50] = 3
sampData$region[sampData$x > 50 & sampData$y < 50] = 4
table(sampData$region)
table(sampData$year)

log_intensity = ggplot(proj_df, aes(x = x, y = y, fill = log_lambda)) + 
  geom_tile() + 
  labs(fill = expression(log(lambda)))

log_intensity_w_p  = log_intensity +
  geom_point(data = sampData, aes(x = x, y = y), size = 0.2, alpha = 0.3, inherit.aes = F)
grid.arrange(log_intensity, log_intensity_w_p, ncol = 2)



sampData$fleet_ndx = fleet_ndx 
head(sampData)

mesh = inla.mesh.2d(max.edge = max.edge, n =10, cutoff = cutoff, loc.domain = SpatialPoints( data.frame(x = c(0,0,grid_dim["x"],grid_dim["x"]), y= c(0,grid_dim["y"],grid_dim["y"],0))))
plot(mesh)
points(sampData$x, sampData$y, pch = 16, col = adjustcolor(col = "red", alpha.f = 0.3))

A = inla.spde.make.A(mesh, loc = cbind(sampData$x, sampData$y))
Proj <- inla.mesh.projector(mesh, loc = cbind(proj_df$x, proj_df$y))
P = Proj$proj$A
spde = inla.spde2.matern(mesh, alpha = 2)

sampData$area = rlnorm(n, log(0.01 * 0.02), 0.4)
sampData$eta =  intercept + normalised_year_coef[sampData$year] + normalised_region_coef[sampData$region] + ifelse(sampData$fleet_ndx == 1, 0, fleet_coef[2] - fleet_coef[1]) + sampData$epsilon
sampData$y_i = rgamma(n = nrow(sampData), shape = shape, scale = (sampData$area * exp(sampData$eta)) / shape)
sampData$catch__per_km_2 = sampData$y_i / sampData$area

hist(sampData$y_i)
summary(sampData$y_i)

full_proj_df$y_i = 1
unique(diff(proj_df$x))
unique(diff(proj_df$y))
full_proj_df$area = 2.04 * 2.04 ## equal area

data = sampData
coordinates(data) <- ~ x + y
coordinates(full_proj_df) <- ~ x + y

## check they all configure correclty
spatial_model_w_epsilon = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = T, 
                                      response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                      spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                      spline_spatial_covariates = NULL, trace_level = "high")

spatial_model_w_epsilon_NN = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = T, 
                                           response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                         spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                         spline_spatial_covariates = NULL, linear_basis = 1, trace_level = "none")

spatial_model_w_epsilon_w_pref_din = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = T, 
                                                 response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                                 spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                                 spline_spatial_covariates = NULL, apply_preferential_sampling = T, preference_model_type = 0,trace_level = "none")

spatial_model_w_epsilon_NN_w_pref_din = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = T, 
                                                    response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                                    spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                                    spline_spatial_covariates = NULL, linear_basis = 1, apply_preferential_sampling = T, preference_model_type = 0, trace_level = "none")

spatial_model_w_epsilon_w_pref_lgcp = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = T, 
                                                  response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                                  spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                                  spline_spatial_covariates = NULL, apply_preferential_sampling = T, preference_model_type = 1,trace_level = "none")

spatial_model_w_epsilon_NN_w_pref_lgcp = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = T, 
                                                     response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                                     spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                                     spline_spatial_covariates = NULL, linear_basis = 1, apply_preferential_sampling = T, preference_model_type = 1, trace_level = "none")

## trace fixed effect 
spatial_model_w_epsilon$obj$env$tracepar = T
spatial_model_w_epsilon_NN$obj$env$tracepar = T
spatial_model_w_epsilon_w_pref_din$obj$env$tracepar = T
spatial_model_w_epsilon_NN_w_pref_din$obj$env$tracepar = T
spatial_model_w_epsilon_w_pref_lgcp$obj$env$tracepar = T
spatial_model_w_epsilon_NN_w_pref_lgcp$obj$env$tracepar = T

## estimate models no preference
start_time = Sys.time()
opt_epsilon = nlminb(spatial_model_w_epsilon$obj$par, spatial_model_w_epsilon$obj$fn, spatial_model_w_epsilon$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
opt_epsilon$opt_time = Sys.time() - start_time

start_time = Sys.time()
opt_epsilon_NN = nlminb(spatial_model_w_epsilon_NN$obj$par, spatial_model_w_epsilon_NN$obj$fn, spatial_model_w_epsilon_NN$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
opt_epsilon_NN$opt_time = Sys.time() - start_time

start_time = Sys.time()
opt_epsilon_NN_din = nlminb(spatial_model_w_epsilon_NN_w_pref_din$obj$par, spatial_model_w_epsilon_NN_w_pref_din$obj$fn, spatial_model_w_epsilon_NN_w_pref_din$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
opt_epsilon_NN_din$opt_time = Sys.time() - start_time

start_time = Sys.time()
opt_epsilon_NN_lgcp = nlminb(spatial_model_w_epsilon_NN_w_pref_lgcp$obj$par, spatial_model_w_epsilon_NN_w_pref_lgcp$obj$fn, spatial_model_w_epsilon_NN_w_pref_lgcp$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
opt_epsilon_NN_lgcp$opt_time = Sys.time() - start_time


rep_epsilon = spatial_model_w_epsilon$obj$report(spatial_model_w_epsilon$obj$env$last.par.best)
rep_epsilon_NN = spatial_model_w_epsilon_NN$obj$report(spatial_model_w_epsilon_NN$obj$env$last.par.best)
rep_epsilon_NN_pref_lgcp = spatial_model_w_epsilon_NN_w_pref_lgcp$obj$report(spatial_model_w_epsilon_NN_w_pref_lgcp$obj$env$last.par.best)
rep_epsilon_NN_pref_din = spatial_model_w_epsilon_NN_w_pref_din$obj$report(spatial_model_w_epsilon_NN_w_pref_din$obj$env$last.par.best)
## estimated marginal standard ev
rep_epsilon$MargSD_epsilon
rep_epsilon_NN$MargSD_epsilon
rep_epsilon_NN_pref_lgcp$MargSD_epsilon
rep_epsilon_NN_pref_din$MargSD_epsilon
## estimated range
rep_epsilon$Range_epsilon
rep_epsilon_NN$Range_epsilon
rep_epsilon_NN_pref_lgcp$Range_epsilon
rep_epsilon_NN_pref_din$Range_epsilon

## spatial plots
proj_epsilon = get_projection(obj = spatial_model_w_epsilon$obj, data = spatial_model_w_epsilon$tmb_data, projection_df = full_proj_df, time_variable_label = "year")
proj_epsilon_NN = get_projection(obj = spatial_model_w_epsilon_NN$obj, data = spatial_model_w_epsilon_NN$tmb_data, projection_df = full_proj_df, time_variable_label = "year")
proj_epsilon_NN_pref_lgcp = get_projection(obj = spatial_model_w_epsilon_NN$obj, data = spatial_model_w_epsilon_NN$tmb_data, projection_df = full_proj_df, time_variable_label = "year")
proj_epsilon_NN_pref_din = get_projection(obj = spatial_model_w_epsilon_NN$obj, data = spatial_model_w_epsilon_NN$tmb_data, projection_df = full_proj_df, time_variable_label = "year")

z_range = range(c(proj_epsilon$predicted_y, proj_epsilon_NN$predicted_y))

## 
years = unique(data@data$year)
n_t = length(years)

for(t in 1:n_t) {
  print(years[t])
  df = data.frame(subset(proj_epsilon, proj_epsilon$year == years[t]))
  plt = ggplot(df, aes(x = x, y = y, fill = predicted_y)) +
    geom_tile() +
    ggtitle("Triangulation") +
    scale_fill_gradientn(colours=c("white","orange","red","dark red"), limits = z_range)
  df = data.frame(subset(proj_epsilon_NN, proj_epsilon_NN$year == years[t]))
  
  plt1 = ggplot(df, aes(x = x, y = y, fill = predicted_y)) +
    geom_tile() +
    ggtitle("nearest Neighbour")       +
    scale_fill_gradientn(colours=c("white","orange","red","dark red"), limits = z_range)
  
  joint_plt = grid.arrange(plt, plt1, ncol = 2)
  ## if you want to save it
  ## ggsave(filename = , plot = joint_plt...)
  Sys.sleep(2)
}


## goodness of fits
rep_opt_epsilon_pref_din = spatial_model_w_epsilon_w_pref_din$obj$report()
rep_opt_epsilon_NN_pref_din = spatial_model_w_epsilon_NN_w_pref_din$obj$report()
rep_opt_epsilon_pref_lgcp = spatial_model_w_epsilon_w_pref_lgcp$obj$report()
rep_opt_epsilon_NN_pref_lgcp = spatial_model_w_epsilon_NN_w_pref_lgcp$obj$report()

sim_data_epsilon_w_pref_din = simulate_data(spatial_model_w_epsilon_w_pref_din$obj, n_sims = 100, fixed_effect = 1, random_effect = 0)
sim_data_epsilon_NN_w_pref_din = simulate_data(spatial_model_w_epsilon_NN_w_pref_din$obj, n_sims = 100, fixed_effect = 1, random_effect = 0)
sim_data_epsilon_w_pref_lgcp = simulate_data(spatial_model_w_epsilon_w_pref_lgcp$obj, n_sims = 100, fixed_effect = 1, random_effect = 0)
sim_data_epsilon_NN_w_pref_lgcp = simulate_data(spatial_model_w_epsilon_NN_w_pref_lgcp$obj, n_sims = 100, fixed_effect = 1, random_effect = 0)
## use DHARMa package for test-statistics

## simulation
N_sims = 100
set.seed(123)
data_ls = no_pref_ls = no_pref_NN_ls =  pref_din = pref_lgcp = list()
for(sim in 1:N_sims) {
  if(sim  %% 5 == 0)
    cat("sim = ", sim, " ")
  sampData = fleet_ndx = NULL;
  full_proj_df = NULL
  for(i in 1:n_years) {
    ## Simulate sample locations based on population spatial distribution
    proj_df$epsilon <- RFsimulate(epsilon_model, x=as.matrix(proj_df[,c("x","y")]), exactness=TRUE)$variable1
    ggplot(proj_df, aes(x = x, y = y, fill = epsilon)) +
      geom_tile()
    
    ## Point process
    intercept_pp = -2.3 # why an intercept model doesn't work, number of observations
    # would be tied to abundance.
    proj_df$log_lambda = intercept_pp + proj_df$epsilon + normalised_region_coef[proj_df$region]
    lam_Im = im(t(matrix(exp(proj_df$log_lambda), nrow=length(predxseq), ncol=length(predyseq))), xcol = predxseq, yrow = predyseq)
    integral(lam_Im)
    
    #ggplot(proj_df, aes(x = x, y = y, fill = omega)) + 
    #  geom_tile()
    # samples_this_year = sum(year_samples == i)
    ## simualte from a point process
    ppp <- rpoispp(lam_Im, win = win)
    
    fleet_ndx = c(fleet_ndx, sample(1:length(fleet_coef), size = ppp$n, replace = T,  prob = c(prob_fleet1_per_year[i], prob_fleet2_per_year[i],prob_fleet3_per_year[i])))
    ##
    ## match the closest depth variable
    ndx = nn2(proj_df[ ,c("x","y")], cbind(ppp$x, ppp$y), k = 1)$nn.idx
    sampData = rbind(sampData, data.frame(year = i, x = ppp$x, y = ppp$y, epsilon = proj_df$epsilon[ndx]))
    
    ##
    proj_df$year = i
    full_proj_df = rbind(full_proj_df, proj_df)
  }
  n = nrow(sampData)
  sampData$region = NA
  sampData$region[sampData$x < 50 & sampData$y > 50] = 1
  sampData$region[sampData$x > 50 & sampData$y > 50] = 2
  sampData$region[sampData$x < 50 & sampData$y < 50] = 3
  sampData$region[sampData$x > 50 & sampData$y < 50] = 4
  table(sampData$region)
  table(sampData$year)
  
  log_intensity = ggplot(proj_df, aes(x = x, y = y, fill = log_lambda)) + 
    geom_tile() + 
    labs(fill = expression(log(lambda)))
  
  log_intensity_w_p  = log_intensity +
    geom_point(data = sampData, aes(x = x, y = y), size = 0.2, alpha = 0.3, inherit.aes = F)
  grid.arrange(log_intensity, log_intensity_w_p, ncol = 2)
  
  
  
  sampData$fleet_ndx = fleet_ndx 
  head(sampData)
  
  mesh = inla.mesh.2d(max.edge = max.edge, n =10, cutoff = cutoff, loc.domain = SpatialPoints( data.frame(x = c(0,0,grid_dim["x"],grid_dim["x"]), y= c(0,grid_dim["y"],grid_dim["y"],0))))
  plot(mesh)
  points(sampData$x, sampData$y, pch = 16, col = adjustcolor(col = "red", alpha.f = 0.3))
  
  A = inla.spde.make.A(mesh, loc = cbind(sampData$x, sampData$y))
  Proj <- inla.mesh.projector(mesh, loc = cbind(proj_df$x, proj_df$y))
  P = Proj$proj$A
  spde = inla.spde2.matern(mesh, alpha = 2)
  
  sampData$area = rlnorm(n, log(0.01 * 0.02), 0.4)
  sampData$eta =  intercept + normalised_year_coef[sampData$year] + normalised_region_coef[sampData$region] + ifelse(sampData$fleet_ndx == 1, 0, fleet_coef[2] - fleet_coef[1]) + sampData$epsilon
  sampData$y_i = rgamma(n = nrow(sampData), shape = shape, scale = (sampData$area * exp(sampData$eta)) / shape)
  sampData$catch__per_km_2 = sampData$y_i / sampData$area
  
  hist(sampData$y_i)
  summary(sampData$y_i)
  
  full_proj_df$y_i = 1
  unique(diff(proj_df$x))
  unique(diff(proj_df$y))
  full_proj_df$area = 2.04 * 2.04 ## equal area
  
  data = sampData
  coordinates(data) <- ~ x + y
  coordinates(full_proj_df) <- ~ x + y
  
  ## check they all configure correclty
  spatial_model_w_epsilon = configure_obj(data = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = T, 
                                          response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                          spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                          spline_spatial_covariates = NULL, trace_level = "none")
  
  spatial_model_w_epsilon_NN = configure_obj(data = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = T, 
                                             response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                             spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                             spline_spatial_covariates = NULL, linear_basis = 1, trace_level = "none")
  
  spatial_model_w_epsilon_w_pref_din = configure_obj(data = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = T, 
                                                     response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                                     spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                                     spline_spatial_covariates = NULL, apply_preferential_sampling = T, preference_model_type = 0,trace_level = "none")
  
  spatial_model_w_epsilon_NN_w_pref_din = configure_obj(data = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = T, 
                                                        response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                                        spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                                        spline_spatial_covariates = NULL, linear_basis = 1, apply_preferential_sampling = T, preference_model_type = 0, trace_level = "none")
  
  spatial_model_w_epsilon_w_pref_lgcp = configure_obj(data = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = T, 
                                                      response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                                      spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                                      spline_spatial_covariates = NULL, apply_preferential_sampling = T, preference_model_type = 1,trace_level = "none")
  
  spatial_model_w_epsilon_NN_w_pref_lgcp = configure_obj(data = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = F, include_epsilon = T, 
                                                         response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                                         spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                                         spline_spatial_covariates = NULL, linear_basis = 1, apply_preferential_sampling = T, preference_model_type = 1, trace_level = "none")
  
  ## estimate models no preference
  start_time = Sys.time()
  opt_epsilon = nlminb(spatial_model_w_epsilon$obj$par, spatial_model_w_epsilon$obj$fn, spatial_model_w_epsilon$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  opt_epsilon$opt_time = Sys.time() - start_time
  
  start_time = Sys.time()
  opt_epsilon_NN = nlminb(spatial_model_w_epsilon_NN$obj$par, spatial_model_w_epsilon_NN$obj$fn, spatial_model_w_epsilon_NN$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  opt_epsilon_NN$opt_time = Sys.time() - start_time
  
  start_time = Sys.time()
  opt_epsilon_NN_din = nlminb(spatial_model_w_epsilon_NN_w_pref_din$obj$par, spatial_model_w_epsilon_NN_w_pref_din$obj$fn, spatial_model_w_epsilon_NN_w_pref_din$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  opt_epsilon_NN_din$opt_time = Sys.time() - start_time
  
  start_time = Sys.time()
  opt_epsilon_NN_lgcp = nlminb(spatial_model_w_epsilon_NN_w_pref_lgcp$obj$par, spatial_model_w_epsilon_NN_w_pref_lgcp$obj$fn, spatial_model_w_epsilon_NN_w_pref_lgcp$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  opt_epsilon_NN_lgcp$opt_time = Sys.time() - start_time
  
  
  rep_epsilon = spatial_model_w_epsilon$obj$report(spatial_model_w_epsilon$obj$env$last.par.best)
  rep_epsilon_NN = spatial_model_w_epsilon_NN$obj$report(spatial_model_w_epsilon_NN$obj$env$last.par.best)
  rep_epsilon_NN_pref_lgcp = spatial_model_w_epsilon_NN_w_pref_lgcp$obj$report(spatial_model_w_epsilon_NN_w_pref_lgcp$obj$env$last.par.best)
  rep_epsilon_NN_pref_din = spatial_model_w_epsilon_NN_w_pref_din$obj$report(spatial_model_w_epsilon_NN_w_pref_din$obj$env$last.par.best)
  
  data_ls[[sim]] = sampData
  no_pref_ls[[sim]] = rep_epsilon
  no_pref_NN_ls[[sim]] = rep_epsilon_NN
  pref_din[[sim]] = rep_epsilon_NN_pref_din
  pref_lgcp[[sim]] = rep_epsilon_NN_pref_lgcp
  
  
  
}

##
# data_ls = no_pref_ls = no_pref_NN_ls =  pref_din = pref_lgcp = list()
## look at parameters
## marginal standard deviation
sd_df = extract_and_merge(par_label = "MargSD_epsilon", reports = list(no_pref_ls, no_pref_NN_ls, pref_din, pref_lgcp),
                          report_labels = c("no_pref","no_pref_NN","pref_din","pref_lgcp"))

ggplot(sd_df, aes(y = value, x = lab, color = lab)) +
  geom_boxplot() +
  geom_abline(intercept = sqrt(epsilon_sigma2), slope = 0) +
  ylim(0,2)

## range
range_df = extract_and_merge(par_label = "Range_epsilon", reports = list(no_pref_ls, no_pref_NN_ls, pref_din, pref_lgcp),
                          report_labels = c("no_pref","no_pref_NN","pref_din","pref_lgcp"))

ggplot(range_df, aes(y = value, x = lab, color = lab)) +
  geom_boxplot() +
  geom_abline(intercept = (epsilon_range), slope = 0) +
  ylim(0,50)

## time-coeffecients
time_df = extract_and_merge(par_label = "time_betas", reports = list(no_pref_ls, no_pref_NN_ls, pref_din, pref_lgcp),
                             report_labels = c("no_pref","no_pref_NN","pref_din","pref_lgcp"))
## get quantiles
sum_time_df = time_df  %>% 
  group_by(lab, element) %>%
  summarize(lower = quantile(value, probs = c(0.025)), mid = quantile(value, probs = c(0.5)), upper = quantile(value, probs = c(0.975)))

true_vals = data.frame(x = 1:n_years, y = normalised_year_coef)

ggplot(sum_time_df, aes(y = mid, x = element, color = lab, fill = lab,  ymin = lower, ymax = upper, alpha = 0.1)) +
  geom_ribbon() +
  ylim(-5,5) +
  geom_line(data = true_vals, inherit.aes = F, aes(x = x, y = y), size = 2) +
  guides(alpha = FALSE) +
  labs(y = "Relative year index") +
  scale_x_continuous(labels = as.character(1:10), breaks = 1:10, name = "Year")

