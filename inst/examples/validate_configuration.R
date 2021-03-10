library(RandomFields)
library(geoR)
library(fields)
library(TMB)
# load INLA for mesh + sparse matrices
library(INLA)
## What about when there is a spatially varying covariate.
library(DSpat)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
#library(mgcv)
source(file.path("R","helper_funs.R"))
source(file.path("R","check_inputs.R"))

library(CPUEspatial)

set.seed(123)
# sample size
n = 1000 # per year
###########################################
# Set field parameters ####################
###########################################

# spatial domain
grid_dim = c("x"=100, "y"=100)
# mesh configuration values
# these are relative to the domain units.
max.edge = c(4, 10)
cutoff = c(6)

# marginal variance of the spatial field
# marginal variance of the spatial-temporal field 
epsilon_sigma2 <- 1   # marginal variance of the eta field
omega_sigma2 <- 1   # marginal variance of the eta field
epsilon_range <- 30     # range-decorrelation parameter
omega_range <- 80   # range-decorrelation parameter
nu = 1             # smoothness parameter  
omega_model = RMmatern(nu = nu, var = omega_sigma2, scale = omega_range / sqrt(8))
epsilon = RMmatern(nu = nu, var = epsilon_sigma2, scale = epsilon_range / sqrt(8))

## linear component
# choose preferential parameter (beta=0 is uniform random)
beta = 1.5
# covariates
hab_coef = c(1, -1.02, -0.5,  1.2) ## contrasts
normalise_hab_coef = hab_coef - mean(hab_coef)

habitat_probs = c(1/3, 2/3, 7/8)
n_years = 10
year_coef = c(rnorm(n_years / 2, 0.5, 0.6),rnorm(n_years / 2, 0, 0.6))
normalised_year_coef = year_coef - mean(year_coef)

#
fleet_coef = c(-0.7, 0.31, -0.2) # simple vessel coeffecients
normalised_fleet_coef = fleet_coef - mean(fleet_coef)
## change fleet over time
prob_fleet1_per_year = seq(from = 0.6, to = 0.2, length.out = n_years)
prob_fleet2_per_year = c(rep(0, n_years - (n_years - 6)), seq(from = 0.2, to = 0.4, length.out = n_years - 6))
prob_fleet3_per_year = 1 - (prob_fleet1_per_year + prob_fleet2_per_year)

intercept = mean(hab_coef) + mean(year_coef) + fleet_coef[1]
## fleet dynamics change over time, one fleet for another
## for habitiat layer
hab.range = 30 ## low number = patch, large = diffuse
# define covariance model for the field GF
shape = 50 ### for gamma response variable
n_h = length(habitat_probs) + 1
model_hab = RMgauss(var = 1, scale = hab.range/sqrt(-log(0.05)))

x_seq = seq(0, grid_dim["x"], 1)
y_seq = seq(0, grid_dim["y"], 1)
coords = cbind(x = rep(x_seq, each = length(y_seq)), y = rep(y_seq, length(x_seq)))

## projection space
predxseq = seq(0, grid_dim["x"],length.out = 50)
unique(round(diff(predxseq),4))
predyseq = seq(0, grid_dim["y"],length.out = 50)
unique(round(diff(predyseq),4))
predgridFull = expand.grid(predxseq, predyseq)

## 
# simulate habitat
habRaw = RFsimulate(model = model_hab, x = coords[, "x"], y = coords[, "y"])@data[, 1]
quan = c(-Inf, quantile(habRaw, habitat_probs), Inf)
habitat = as.numeric(cut(habRaw, quan))
covariates = data.frame(x = coords[,"x"], y = coords[,"y"], habitat = habitat)
covariates$x = covariates$x
covariates$y = covariates$y
covariates$omega = RFsimulate(omega_model, x = as.matrix(covariates[,c("x","y")]))$variable1
popdata = covariates
popdata$hab_factor = hab_coef[covariates$habitat]
popdata$norm_hab_factor = normalise_hab_coef[covariates$habitat]
## set up pop Dataframe and simulate epsilon
time_popdata = do.call("rbind", replicate(n_years, popdata, simplify = FALSE))
time_popdata$year = sort(rep(1:n_years, nrow(popdata)))
time_popdata$space_time = NA
dim(time_popdata)
# Simulate Population Epsilon assumed iid
for(i in 1:n_years) {
  this_year_GF <- RFsimulate(epsilon, x = as.matrix(covariates[,c("x","y")]))$variable1
  time_popdata$epsilon[time_popdata$year == i] = this_year_GF
}
time_popdata$spatial_component = time_popdata$epsilon + time_popdata$norm_hab_factor + time_popdata$omega
# Randomly select sample size per year roughly equal sampling in each year
year_samples = sample(1:n_years, size = n, replace = T,  prob = rep(1,n_years) / n_years)
sampData = NULL;
fleet_ndx = NULL;
## Simulate sample locations based on population spatial distribution
projection_df = NULL;
for(i in 1:n_years) {
  samples_this_year = sum(year_samples == i)
  pop_ndx = time_popdata$year == i
  ind <- sample(which(pop_ndx), size = samples_this_year, replace = T, prob = exp(beta * time_popdata$spatial_component[pop_ndx]))
  sampData = rbind(sampData, time_popdata[ind, ])
  fleet_ndx = c(fleet_ndx, sample(1:length(fleet_coef), size = samples_this_year, replace = T,  prob = c(prob_fleet1_per_year[i], prob_fleet2_per_year[i],prob_fleet3_per_year[i])))
  
  ##
  this_df = time_popdata[pop_ndx, ]
  pop_red_ndx = RANN::nn2(this_df[,c(1,2)], predgridFull[ ,c(1,2)], k = 1)$nn.idx
  temp_df = this_df[pop_red_ndx,]
  temp_df$year = i
  projection_df = rbind(projection_df, this_df[pop_red_ndx,])

}
projection_df$area = 2.0408 * 2.0408
projection_df$y_i = 1
sampData$fleet_ndx = fleet_ndx 
time_popdata$relative_spatial = exp(normalised_year_coef[time_popdata$year] + time_popdata$spatial_component)

pop_over_space = ggplot(data = time_popdata, aes(x = x, y = y, fill = relative_spatial)) +
  geom_tile() + 
  labs(fill = "Relative\nCatch rate") +
  facet_wrap( ~ year, nrow = 2, ncol = 5 )
#pop_over_space
#ggsave(plot = pop_over_space, filename = file.path(FIG_DIR, "population", paste0("Relative_spatial_population_",sim,".png")), width = 12, height = 9)

## attach pop locations with Prediction locations for covariates
pop_red_ndx = as.data.frame(RANN::nn2(popdata[,c(1,2)], predgridFull[ ,c(1,2)], k = 1))
## Create GMRF objects
mesh = inla.mesh.2d(max.edge = max.edge, n =10, cutoff = cutoff, loc.domain = SpatialPoints( data.frame(x = c(0,0,grid_dim["x"],grid_dim["x"]), y= c(0,grid_dim["y"],grid_dim["y"],0))))
A = inla.spde.make.A(mesh, loc = cbind(sampData$x, sampData$y))
Proj <- inla.mesh.projector(mesh, loc = cbind(predgridFull$Var1, predgridFull$Var2))
P = Proj$proj$A
spde = inla.spde2.matern(mesh, alpha = 2)
## visualise the mesh wtih respect to sample data
plot(mesh)
points(sampData$x, sampData$y, col = adjustcolor(col = "black", alpha = 0.2), pch = 16)

## create projection model matrices for spatial covariates
temp_pred_df = predgridFull
X_spatial_proj_zpt = array(0, dim = c(nrow(predgridFull), length(hab_coef), n_years))
X_spatial_ipt = array(0, dim = c(nrow(sampData), length(hab_coef), n_years))
for(t in 1:n_years) {
  pop_ndx = time_popdata$year == t
  temp_pred_df$hab = time_popdata[pop_ndx, ][pop_red_ndx[,1], "habitat"]
  X_spatial_proj_zpt[,,t] = model.matrix(Var1 ~ 0 + factor(hab), data = temp_pred_df)
  ## sample one
  sub_obs = sampData[sampData$year == t,]
  mod_mat =  model.matrix(y ~ 0 + factor(habitat), data = sampData)
  X_spatial_ipt[,,t] = mod_mat
  ## check correct projectiosn
  #hab_plot = ggplot(data = temp_pred_df, aes(x = Var1, y = Var2, fill = factor(hab))) +
  #  geom_tile() + 
  #  labs(fill = "Habitat Type")
  #print(hab_plot)
  #Sys.sleep(3)
}

## model matricies
time_model_matrix = model.matrix(x ~ 0 + factor(year),  data = sampData)
model_matrix = model.matrix(x ~ factor(fleet_ndx),  data = sampData)
# randomly assing different area swept
sampData$area = rlnorm(n, log(0.01 * 0.02), 0.4)

## add observation error to Y's
## simualte observations based on a Gamma distribution
sampData$eta =  intercept + sampData$spatial_component + normalised_year_coef[sampData$year] + ifelse(sampData$fleet_ndx == 1, 0, fleet_coef[2] - fleet_coef[1])
sampData$y_i = rgamma(n = nrow(sampData), shape = shape, scale = (sampData$area * exp(sampData$eta)) / shape)

data = sampData
coordinates(data) <- ~ x + y
coordinates(projection_df) <- ~ x + y

## debug configure_obj.R
data = data
include_epsilon = T
include_omega = T
response_variable_label = "y_i" 
time_variable_label = "year"
catchability_covariates = NULL
catchability_covariate_type =NULL
spatial_covariates = NULL
spatial_covariate_type = NULL
spline_catchability_covariates = NULL
spline_spatial_covariates = NULL
projection_df = projection_df
mesh = mesh
family = 3
link = 0
trace_level = "none"
linear_basis = 0
warnings()

detach("package:CPUEspatial", unload=TRUE)
library(CPUEspatial)


## check they all configure correclty
simple_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = T, include_epsilon = T, 
                      response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, catchability_covariate_type = NULL, 
                      spatial_covariates = NULL, spatial_covariate_type = NULL, spline_catchability_covariates = NULL,
                      spline_spatial_covariates = NULL, trace_level = "high")


#obj_simple <- MakeADFun(simple_model$tmb_data, simple_model$tmb_pars, random = c("epsilon_input","omega_input"), DLL = "CPUEspatial_TMBExports", method = "nlminb", hessian = T, silent=T)

single_catch_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = T, include_epsilon = T, 
                             response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                             spatial_covariates = NULL, spatial_covariate_type = NULL, spline_catchability_covariates = NULL,
                             spline_spatial_covariates = NULL, trace_level = "high")
#obj_single_catch <- MakeADFun(single_catch_model$tmb_data, single_catch_model$tmb_pars, random = c("epsilon_input","omega_input"), DLL = "CPUEspatial_TMBExports", method = "nlminb", hessian = T, silent=T)

single_catch_sptial_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = T, include_epsilon = T, 
                                   response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                   spatial_covariates = "habitat", spatial_covariate_type = "factor", spline_catchability_covariates = NULL,
                                   spline_spatial_covariates = NULL, trace_level = "high")
#obj_single_catch_sptial <- MakeADFun(single_catch_sptial_model$tmb_data, single_catch_sptial_model$tmb_pars, random = c("epsilon_input","omega_input"), DLL = "CPUEspatial_TMBExports", method = "nlminb", hessian = T, silent=T)

single_catch_sptial_model_num = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = T, include_epsilon = T, 
                                          response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = "fleet_ndx", catchability_covariate_type = "factor", 
                                          spatial_covariates = "omega", spatial_covariate_type = "numeric", spline_catchability_covariates = NULL,
                                          spline_spatial_covariates = NULL, trace_level = "high")
#obj_single_catch_sptial <- MakeADFun(single_catch_sptial_model$tmb_data, single_catch_sptial_model$tmb_pars, random = c("epsilon_input","omega_input"), DLL = "CPUEspatial_TMBExports", method = "nlminb", hessian = T, silent=T)

double_catch_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = T, include_epsilon = T, 
                                   response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = c("fleet_ndx","spatial_component"), catchability_covariate_type = c("factor","numeric"), 
                                   spatial_covariates = NULL, spatial_covariate_type = NULL, spline_catchability_covariates = NULL,
                                   spline_spatial_covariates = NULL, trace_level = "high")
#obj_double_catch <- MakeADFun(double_catch_model$tmb_data, double_catch_model$tmb_pars, random = c("epsilon_input","omega_input"), DLL = "CPUEspatial_TMBExports", method = "nlminb", hessian = T, silent=T)

double_catch_sptial_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = T, include_epsilon = T, 
                                          response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = c("fleet_ndx","spatial_component"), catchability_covariate_type = c("factor","numeric"), 
                                          spatial_covariates = c("habitat", "spatial_component"), spatial_covariate_type = c("factor","numeric"), spline_catchability_covariates = NULL,
                                          spline_spatial_covariates = NULL, trace_level = "high")
#obj_double_catch_sptial <- MakeADFun(double_catch_sptial_model$tmb_data, double_catch_sptial_model$tmb_pars, random = c("epsilon_input","omega_input"), DLL = "CPUEspatial_TMBExports", method = "nlminb", hessian = T, silent=T)

## marginal log-likelihood
simple_model$obj$fn()
single_catch_model$obj$fn()
single_catch_sptial_model$obj$fn()
single_catch_sptial_model_num$obj$fn()
double_catch_model$obj$fn()
double_catch_sptial_model$obj$fn()
## gradient
simple_model$obj$gr()
single_catch_model$obj$gr()
single_catch_sptial_model$obj$gr()
single_catch_sptial_model_num$obj$gr()
double_catch_model$obj$gr()
double_catch_sptial_model$obj$fn()
## look at the pars
simple_model$obj$par
single_catch_model$obj$par
single_catch_sptial_model$obj$par
single_catch_sptial_model_num$obj$par
double_catch_model$obj$par
double_catch_sptial_model$obj$par

## check spline congigurations
spline_catch_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = T, include_epsilon = T, 
                             response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, catchability_covariate_type = NULL, 
                             spatial_covariates = NULL, spatial_covariate_type = NULL, spline_catchability_covariates = "omega",
                             spline_spatial_covariates = NULL, trace_level = "high")

spline_spatial_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = T, include_epsilon = T, 
                                   response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, catchability_covariate_type = NULL, 
                                   spatial_covariates = NULL, spatial_covariate_type = NULL, spline_catchability_covariates = NULL,
                                   spline_spatial_covariates = "omega", trace_level = "high")

spline_both_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = T, include_epsilon = T, 
                                     response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, catchability_covariate_type = NULL, 
                                     spatial_covariates = NULL, spatial_covariate_type = NULL, spline_catchability_covariates = "omega",
                                     spline_spatial_covariates = "epsilon", trace_level = "high")

multi_spline_both_model = configure_obj(data = data, projection_df = projection_df, mesh = mesh, family = 3, link = 0, include_omega = T, include_epsilon = T, 
                                  response_variable_label = "y_i", time_variable_label = "year", catchability_covariates = NULL, catchability_covariate_type = NULL, 
                                  spatial_covariates = NULL, spatial_covariate_type = NULL, spline_catchability_covariates = c("omega","epsilon"),
                                  spline_spatial_covariates = c("omega","epsilon"), trace_level = "high")



simple_model$obj$fn()
single_catch_model$obj$fn()
single_catch_sptial_model$obj$fn()
double_catch_model$obj$fn()
double_catch_sptial_model$obj$fn()



compile(file.path("src","debug_standalone_version.cpp"))
file.exists(file.path("src","debug_standalone_version.dll"))
dyn.load(dynlib(file.path("src","debug_standalone_version")))


