#'
#' A simulation script to look at preference vs non-preference sampling
#' with just time-varying GF.
#' gamma response variable with discrete spatail area plus numeric linear spatial variable fleet and year effects
#'
library(TMB)
library(CPUEspatial)
library(INLA)
library(RandomFields)
library(RANN)
library(ggplot2)
library(mgcv) 
library(reshape2)
library(dplyr)
library(gratia) # useful for plotting splines/smoother functions https://github.com/gavinsimpson/gratia
source(file.path("inst","examples","canonical.index.R"))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
set.seed(123)
# spatial domain even though we are not doing a spatial model
# this model requires it.
grid_dim = c("x"=100, "y"=100)
max.edge = c(10)
cutoff = c(6)
# marginal variance of the spatial-temporal field 
epsilon_sigma2 <- 1  # marginal variance of the eta field
epsilon_range <- 30     # range-decorrelation parameter
nu = 1             # smoothness parameter  
#epsilon_model = RMwhittle(nu = nu, var = epsilon_sigma2, scale = epsilon_range )
epsilon_model = RMmatern(nu = nu, var = epsilon_sigma2, scale = epsilon_range / sqrt(8))

# sample size
n = 1000 # per year
n_years = 10
#pref_beta = c(rep(3, 5), rep(1.5, 5)) ## constant preference coeffecient
pref_beta = rep(1.5, n_years) ## constant preference coeffecient

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

## look at the coefs and expectation
sys_range = range(rowSums(expand.grid(year_coef,region_coef, fleet_coef)))
area_range = range(rlnorm(10000, log(0.01 * 0.02), 0.4))
## expected
shape * area_range * exp(sys_range)
scale_range =  1/shape * area_range * exp(sys_range)
var_range = shape * scale_range^2
cv = 1/sqrt(shape)
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
proj_df$area = unique(round(diff(predxseq),4)) * unique(round(diff(predyseq),4))

## area cells
plot(1, 1, xlim = c(0, grid_dim["x"]), ylim = c(0, grid_dim["y"]), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i")
abline(v = c(seq(0,grid_dim["x"], by = 50)), lwd = 3)
abline(h = c(seq(0,grid_dim["y"], by = 50)), lwd = 3)
text(x = 25, y  = 75, labels = "A", cex = 3)
text(x = 75, y  = 75, labels = "B", cex = 3)
text(x = 25, y  = 25, labels = "C", cex = 3)
text(x = 75, y  = 25, labels = "D", cex = 3)
year_samples = sample(1:n_years, size = n * n_years, replace = T,  prob = rep(1,n_years) / n_years)

rep_ls = opt_ls = sd_ls = data_ls = list()
true_index_ls = list();
rep_pref_ls = opt_pref_ls = sd_pref_ls  = list()
N_sims = 100
################################################
# This simulation can take some time FYI
# load(file = "GF_omega_pref.RData")
#
################################################
for(sim in 1:N_sims) {
  if(sim %% 5 == 0)
    cat("sim = ", sim, "\n")
  ## Simulate sample locations based on population spatial distribution
  
  sampData = fleet_ndx = NULL;
  full_proj_df = NULL
  for(i in 1:n_years) {
    proj_df$epsilon <- RFsimulate(epsilon_model, x=as.matrix(proj_df[,c("x","y")]), exactness=TRUE)$variable1
    ## spatial expected catch rate
    proj_df$eta =  intercept + normalised_year_coef[i] + normalised_region_coef[proj_df$region] +  proj_df$epsilon
    
    #ggplot(proj_df, aes(x = x, y = y, fill = epsilon)) + 
    #  geom_tile()
    samples_this_year = sum(year_samples == i)
    ind <- sample(1:nrow(proj_df), size = samples_this_year, replace = T, prob = exp(pref_beta[i] * proj_df$eta))
    this_samp =proj_df[ind, ]
    fleet_ndx = c(fleet_ndx, sample(1:length(fleet_coef), size = samples_this_year, replace = T,  prob = c(prob_fleet1_per_year[i], prob_fleet2_per_year[i],prob_fleet3_per_year[i])))
    ##
    sampData = rbind(sampData, data.frame(year = i, x = this_samp$x, y = this_samp$y, region = this_samp$region, epsilon = this_samp$epsilon))
    
    ##
    proj_df$year = i
    full_proj_df = rbind(full_proj_df, proj_df)
  }
  
  table(sampData$region)
  
  sampData$fleet_ndx = fleet_ndx 
  head(sampData)
  
  mesh = inla.mesh.2d(loc = cbind(sampData$x, sampData$y), max.edge = max.edge, n =10, cutoff = cutoff, loc.domain = SpatialPoints( data.frame(x = c(0,0,grid_dim["x"],grid_dim["x"]), y= c(0,grid_dim["y"],grid_dim["y"],0))))
  plot(mesh)
  points(sampData$x, sampData$y, pch = 16, col = adjustcolor(col = "red", alpha.f = 0.1))
  
  A = inla.spde.make.A(mesh, loc = cbind(sampData$x, sampData$y))
  Proj <- inla.mesh.projector(mesh, loc = cbind(proj_df$x, proj_df$y))
  P = Proj$proj$A
  spde = inla.spde2.matern(mesh, alpha = 2)
  
  sampData$area = rlnorm(nrow(sampData), log(0.01 * 0.02), 0.4)
  sampData$eta =  intercept + normalised_year_coef[sampData$year] + normalised_region_coef[sampData$region] + ifelse(sampData$fleet_ndx == 1, 0, fleet_coef[2] - fleet_coef[1]) + sampData$epsilon
  sampData$y_i = rgamma(n = nrow(sampData), shape = shape, scale = (sampData$area * exp(sampData$eta)) / shape)
  sampData$catch__per_km_2 = sampData$y_i / sampData$area
  
  hist(sampData$y_i)
  summary(sampData$y_i)
  
  full_proj_df$y_i = 1
  #full_proj_df$area = 1 ## equal area
  
  sampData$region = factor(sampData$region)
  sampData$fleet_ndx = factor(sampData$fleet_ndx)
  full_proj_df$region = factor(full_proj_df$region)
  
  data = sampData
  coordinates(data) <- ~ x + y
  coordinates(full_proj_df) <- ~ x + y
  data_ls[[sim]] = sampData
  ## check they all configure correclty
  spatial_model_w_omega = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = T, include_epsilon = F, 
                                        response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                                        spatial_covariates = c("region"),  spline_catchability_covariates = NULL,
                                        spline_spatial_covariates = NULL, trace_level = "none")
  # tips on debugging
  opt_spa = nlminb(spatial_model_w_omega$obj$par, spatial_model_w_omega$obj$fn, spatial_model_w_omega$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  rep_epsilon = spatial_model_w_omega$obj$report(spatial_model_w_omega$obj$env$last.par.best)
  sd_rep = sdreport(spatial_model_w_omega$obj)
  
  opt_ls[[sim]] = opt_spa
  rep_ls[[sim]] = rep_epsilon
  sd_ls[[sim]] = sd_rep
  
  spatial_model_w_omega_pref = configure_obj(observed_df = data, projection_df = full_proj_df, mesh = mesh, family = 2, link = 0, include_omega = T, include_epsilon = F, 
                                             response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", 
                                             spatial_covariates = c("region"),  spline_catchability_covariates = NULL,
                                             spline_spatial_covariates = NULL, apply_preferential_sampling = T, preference_model_type = 0, trace_level = "none")
  
  opt_pref_spa = nlminb(spatial_model_w_omega_pref$obj$par, spatial_model_w_omega_pref$obj$fn, spatial_model_w_omega_pref$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  rep_pref_epsilon = spatial_model_w_omega_pref$obj$report(spatial_model_w_omega_pref$obj$env$last.par.best)
  sd_pref_rep = sdreport(spatial_model_w_omega_pref$obj)
  
  rep_pref_ls[[sim]] = rep_pref_epsilon
  opt_pref_ls[[sim]] = opt_pref_spa
  sd_pref_ls[[sim]] = sd_pref_rep
  
  true_ndx = tapply(full_proj_df$area * exp(year_coef[full_proj_df$year] + region_coef[full_proj_df$region] + full_proj_df$epsilon), INDEX = full_proj_df$year, FUN = sum)
  true_index_ls[[sim]] = true_ndx
  ## Compare indices
  par(mfrow = c(1,2))
  plot(1:n_years, normalised_year_coef, lwd= 3, lty = 1, col = "black", type  ="l", main = "Year coeffecients", ylab = "", xlab = "time")
  lines(1:n_years, rep_pref_epsilon$time_betas, lwd = 3, lty = 2, col = "blue")
  lines(1:n_years, rep_epsilon$time_betas, lwd = 3, lty = 3, col = "red")
  
  plot(1:n_years, true_ndx / gm_mean(true_ndx), lwd = 3, lty = 1, type = "l", main = "Relative index", ylab = "", xlab = "time")
  lines(1:n_years, rep_pref_epsilon$relative_index / gm_mean(rep_pref_epsilon$relative_index), lty = 2, col = "blue", lwd =3)
  lines(1:n_years, rep_epsilon$relative_index / gm_mean(rep_epsilon$relative_index), lty = 3, col = "red", lwd =3)
  legend("topright", col = c("blue","red","black"), legend = c("Preference","No Preference", "OM"), lwd =3, cex = 0.8)
  ##
}
#
save(rep_ls, opt_ls, sd_ls, data_ls, rep_pref_ls, opt_pref_ls, sd_pref_ls, true_index_ls, file = "GF_TV_pref_1.5.RData")

est_sd_omega = Reduce(c, lapply(rep_ls, FUN = function(x) {x$MargSD_omega}))
est_range_omega = Reduce(c, lapply(rep_ls, FUN = function(x) {x$Range_omega}))

est_sd_omega_pref = Reduce(c, lapply(rep_pref_ls, FUN = function(x) {x$MargSD_omega}))
est_range_omega_pref = Reduce(c, lapply(rep_pref_ls, FUN = function(x) {x$Range_omega}))
est_omega_pref_pref = Reduce(c, lapply(rep_pref_ls, FUN = function(x) {x$pref_coef}))

boxplot(est_omega_pref_pref, ylim = c(0,2))
abline(h = unique(pref_beta), lwd = 3, lty = 2, col = "red")

boxplot(cbind(est_sd_omega, est_sd_omega_pref))
abline(h = sqrt(epsilon_sigma2), lwd = 3, lty = 2, col = "red")

boxplot(cbind(est_range_omega, est_range_omega_pref))
abline(h = (epsilon_range), lwd = 3, lty = 2, col = "red")

## get coeffecients and relative index
relative_index = Reduce(rbind, lapply(rep_ls, FUN = function(x) {x$relative_index}))
relative_index_pref = Reduce(rbind, lapply(rep_pref_ls, FUN = function(x) {x$relative_index}))
standardised_index = Reduce(rbind, lapply(rep_ls, FUN = function(x) {x$relative_index / gm_mean(x$relative_index)}))
standardised_index_pref = Reduce(rbind, lapply(rep_pref_ls, FUN = function(x) {x$relative_index / gm_mean(x$relative_index)}))
year_index = Reduce(rbind, lapply(rep_ls, FUN = function(x) {x$time_betas}))
year_index_pref = Reduce(rbind, lapply(rep_pref_ls, FUN = function(x) {x$time_betas}))
## plot relative index, geometric index
rownames(year_index) = rownames(year_index_pref) = rownames(standardised_index) = rownames(relative_index)  = rownames(relative_index_pref) = rownames(standardised_index_pref) = paste0("sim_", 1:nrow(relative_index))

melt_year = melt(year_index)
melt_year_pref = melt(year_index_pref)
melt_year$model = "No Preference"
melt_year_pref$model = "Preference"
year_df = rbind(melt_year, melt_year_pref)
colnames(year_df) = c("sim", "time", "value", "model")


melt_stand = melt(standardised_index)
melt_stand_pref = melt(standardised_index_pref)
melt_stand$model = "No Preference"
melt_stand_pref$model = "Preference"
stand_df = rbind(melt_stand, melt_stand_pref)
colnames(stand_df) = c("sim", "time", "value", "model")

melt_relative = melt(relative_index)
melt_relative_pref = melt(relative_index_pref)
melt_relative$model = "No Preference"
melt_relative_pref$model = "Preference"
relative_df = rbind(melt_relative, melt_relative_pref)
colnames(relative_df) = c("sim", "time", "value", "model")

## get quanitles
relative_summary = relative_df %>% 
  group_by(time, model) %>% 
  summarize(lower = quantile(value, 0.025), upper = quantile(value, 0.975), mid = quantile(value, 0.5))

standardise_summary = stand_df %>% 
  group_by(time, model) %>% 
  summarize(lower = quantile(value, 0.025), upper = quantile(value, 0.975), mid = quantile(value, 0.5))

year_summary = year_df %>% 
  group_by(time, model) %>% 
  summarize(lower = quantile(value, 0.025), upper = quantile(value, 0.975), mid = quantile(value, 0.5))

## visualise
ggplot(relative_summary, aes(x = time, y = mid, col = model, alpha = 0.1, fill = model)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper)) +
  geom_line(size = 2) +
  scale_alpha(guide = "none") +
  ylab("Relative index")

ggplot(standardise_summary, aes(x = time, y = mid, col = model, alpha = 0.1, fill = model)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper)) +
  geom_line(size = 2) + 
  scale_alpha(guide = "none") +
  ylab("Geometric index")

true_coef_df = data.frame(time = 1:n_years, true_year_coef = normalised_year_coef)
ggplot(year_summary, aes(x = time, y = mid, col = model, alpha = 0.1, fill = model)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper)) +
  geom_line(size = 2) + 
  scale_alpha(guide = "none") + 
  geom_line(data = true_coef_df, aes(x = time, y = true_year_coef), inherit.aes = F, size = 1.2, linetype = "dashed") + 
  ylab("Normalised year coeffecients")

