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
library(mgcv) 
library(spatstat)
library(DHARMa)
source(file.path("inst","examples","canonical.index.R"))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
set.seed(123)
# spatial domain even though we are not doing a spatial model
# this model requires it.
grid_dim = c("x"=100, "y"=100)
max.edge = c(4, 10)
cutoff = c(6)
# marginal variance of the spatial-temporal field 
omega_sigma2 <- 1  # marginal variance of the eta field
omega_range <- 30     # range-decorrelation parameter
nu = 1             # smoothness parameter  
#omega_model = RMwhittle(nu = nu, var = omega_sigma2, scale = omega_range )
omega_model = RMmatern(nu = nu, var = omega_sigma2, scale = omega_range / sqrt(8))

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

rep_ls = opt_ls = rep_NN_ls = opt_NN_ls = list()
## Simulate sample locations based on population spatial distribution
proj_df$omega <- RFsimulate(omega_model, x=as.matrix(proj_df[,c("x","y")]), exactness=TRUE)$variable1
## Point process
intercept_pp = -2.3 # why an intercept model doesn't work, number of observations
# would be tied to abundance.
proj_df$log_lambda = intercept_pp + proj_df$omega + normalised_region_coef[proj_df$region]
lam_Im = im(t(matrix(exp(proj_df$log_lambda), nrow=length(predxseq), ncol=length(predyseq))), xcol = predxseq, yrow = predyseq)
integral(lam_Im)

sampData = fleet_ndx = NULL;
full_proj_df = NULL
for(i in 1:n_years) {
  #ggplot(proj_df, aes(x = x, y = y, fill = omega)) + 
  #  geom_tile()
  # samples_this_year = sum(year_samples == i)
  ## simualte from a point process
  ppp <- rpoispp(lam_Im, win = win)

  fleet_ndx = c(fleet_ndx, sample(1:length(fleet_coef), size = ppp$n, replace = T,  prob = c(prob_fleet1_per_year[i], prob_fleet2_per_year[i],prob_fleet3_per_year[i])))
  ##
  sampData = rbind(sampData, data.frame(year = i, x = ppp$x, y = ppp$y))
  
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

## match the closest depth variable
ndx = nn2(proj_df[ ,c("x","y")], sampData[,c("x","y")], k = 1)$nn.idx
sampData$omega = proj_df$omega[ndx]

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
sampData$eta =  intercept + normalised_year_coef[sampData$year] + normalised_region_coef[sampData$region] + ifelse(sampData$fleet_ndx == 1, 0, fleet_coef[2] - fleet_coef[1]) + sampData$omega
sampData$y_i = rgamma(n = nrow(sampData), shape = shape, scale = (sampData$area * exp(sampData$eta)) / shape)
sampData$catch__per_km_2 = sampData$y_i / sampData$area

hist(sampData$y_i)
summary(sampData$y_i)

full_proj_df$y_i = 1
full_proj_df$area = 1 ## equal area

data = sampData
coordinates(data) <- ~ x + y
coordinates(full_proj_df) <- ~ x + y

## check they all configure correclty
spatial_model_w_omega = configure_obj(data = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = T, include_epsilon = F, 
                                     response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                     spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                     spline_spatial_covariates = NULL, trace_level = "none")

spatial_model_w_omega_NN = configure_obj(data = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = T, include_epsilon = F, 
                                      response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                      spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                      spline_spatial_covariates = NULL, linear_basis = 1, trace_level = "none")

spatial_model_w_omega_w_pref_din = configure_obj(data = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = T, include_epsilon = F, 
                                      response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                      spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                      spline_spatial_covariates = NULL, apply_preferential_sampling = T, preference_model_type = 0,trace_level = "none")

spatial_model_w_omega_NN_w_pref_din = configure_obj(data = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = T, include_epsilon = F, 
                                         response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                         spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                         spline_spatial_covariates = NULL, linear_basis = 1, apply_preferential_sampling = T, preference_model_type = 0, trace_level = "none")

spatial_model_w_omega_w_pref_lgcp = configure_obj(data = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = T, include_epsilon = F, 
                                                 response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                                 spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                                 spline_spatial_covariates = NULL, apply_preferential_sampling = T, preference_model_type = 1,trace_level = "none")

spatial_model_w_omega_NN_w_pref_lgcp = configure_obj(data = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = T, include_epsilon = F, 
                                                    response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                                    spatial_covariates = c("region"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                                                    spline_spatial_covariates = NULL, linear_basis = 1, apply_preferential_sampling = T, preference_model_type = 1, trace_level = "none")

spatial_model_w_omega_w_pref_din$obj$gr()
spatial_model_w_omega_NN_w_pref_din$obj$gr()
spatial_model_w_omega_w_pref_lgcp$obj$gr()
spatial_model_w_omega_NN_w_pref_lgcp$obj$gr()

## trace fixed effect 
spatial_model_w_omega_w_pref_din$obj$env$tracepar = T
spatial_model_w_omega_NN_w_pref_din$obj$env$tracepar = T
spatial_model_w_omega_w_pref_lgcp$obj$env$tracepar = T
spatial_model_w_omega_NN_w_pref_lgcp$obj$env$tracepar = T

## estimate models
start_time = Sys.time()
opt_omega_pref_din = nlminb(spatial_model_w_omega_w_pref_din$obj$par, spatial_model_w_omega_w_pref_din$obj$fn, spatial_model_w_omega_w_pref_din$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
opt_omega_pref_din$opt_time = Sys.time() - start_time
start_time = Sys.time()
opt_omega_NN_pref_din = nlminb(spatial_model_w_omega_NN_w_pref_din$obj$par, spatial_model_w_omega_w_pref_din$obj$fn, spatial_model_w_omega_NN_w_pref_din$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
opt_omega_NN_pref_din$opt_time = Sys.time() - start_time
start_time = Sys.time()
opt_omega_pref_lgcp = nlminb(spatial_model_w_omega_w_pref_lgcp$obj$par, spatial_model_w_omega_w_pref_lgcp$obj$fn, spatial_model_w_omega_w_pref_lgcp$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
opt_omega_pref_lgcp$opt_time = Sys.time() - start_time
start_time = Sys.time()
opt_omega_NN_pref_lgcp = nlminb(spatial_model_w_omega_NN_w_pref_lgcp$obj$par, spatial_model_w_omega_NN_w_pref_lgcp$obj$fn, spatial_model_w_omega_NN_w_pref_lgcp$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
opt_omega_NN_pref_lgcp$opt_time = Sys.time() - start_time

opt_omega_pref_din$opt_time
opt_omega_NN_pref_din$opt_time
opt_omega_pref_lgcp$opt_time 
opt_omega_NN_pref_lgcp$opt_time
## estimate preference parameter
rep_opt_omega_pref_din = spatial_model_w_omega_w_pref_din$obj$report()
rep_opt_omega_NN_pref_din = spatial_model_w_omega_NN_w_pref_din$obj$report()
rep_opt_omega_pref_lgcp = spatial_model_w_omega_w_pref_lgcp$obj$report()
rep_opt_omega_NN_pref_lgcp = spatial_model_w_omega_NN_w_pref_lgcp$obj$report()

rep_opt_omega_pref_din$pref_coef
rep_opt_omega_NN_pref_din$pref_coef
rep_opt_omega_pref_lgcp$pref_coef
rep_opt_omega_NN_pref_lgcp$pref_coef
opt_omega_pref_lgcp$par["lgcp_intercept"]
opt_omega_NN_pref_lgcp$par["lgcp_intercept"]

#opt_spa_NN = nlminb(spatial_model_w_omega_NN$obj$par, spatial_model_w_omega_NN$obj$fn, spatial_model_w_omega_NN$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
#rep_omega = spatial_model_w_omega$obj$report()
#rep_omega_NN = spatial_model_w_omega_NN$obj$report()

#sd_rep = sdreport(simple_spatial_model$obj)
#sd_rep_NN = sdreport(spatial_model_w_omega_NN$obj)
rep_opt_omega_pref_din$MargSD_omega
rep_opt_omega_NN_pref_din$MargSD_omega
rep_opt_omega_pref_lgcp$MargSD_omega
rep_opt_omega_NN_pref_lgcp$MargSD_omega

rep_opt_omega_pref_din$Range_omega
rep_opt_omega_NN_pref_din$Range_omega
rep_opt_omega_pref_lgcp$Range_omega
rep_opt_omega_NN_pref_lgcp$Range_omega

## spatial plots
proj_omega_din = get_projection(obj = spatial_model_w_omega_w_pref_din$obj, data = spatial_model_w_omega_w_pref_din$tmb_data, projection_df = full_proj_df, time_variable_label = "year")
proj_omega_NN_din = get_projection(obj = spatial_model_w_omega_NN_w_pref_din$obj, data = spatial_model_w_omega_NN_w_pref_din$tmb_data, projection_df = full_proj_df, time_variable_label = "year")
z_range = range(proj_omega_din$predicted_y)
z_range_NN = range(proj_omega_NN_din$predicted_y)

## goodness of fits
rep_opt_omega_pref_din = spatial_model_w_omega_w_pref_din$obj$report()
rep_opt_omega_NN_pref_din = spatial_model_w_omega_NN_w_pref_din$obj$report()
rep_opt_omega_pref_lgcp = spatial_model_w_omega_w_pref_lgcp$obj$report()
rep_opt_omega_NN_pref_lgcp = spatial_model_w_omega_NN_w_pref_lgcp$obj$report()

sim_data_omega_w_pref_din = simulate_data(spatial_model_w_omega_w_pref_din$obj, n_sims = 100, fixed_effect = 1, random_effect = 0)
sim_data_omega_NN_w_pref_din = simulate_data(spatial_model_w_omega_NN_w_pref_din$obj, n_sims = 100, fixed_effect = 1, random_effect = 0)
sim_data_omega_w_pref_lgcp = simulate_data(spatial_model_w_omega_w_pref_lgcp$obj, n_sims = 100, fixed_effect = 1, random_effect = 0)
sim_data_omega_NN_w_pref_lgcp = simulate_data(spatial_model_w_omega_NN_w_pref_lgcp$obj, n_sims = 100, fixed_effect = 1, random_effect = 0)
## use DHARMa package for test-statistics


