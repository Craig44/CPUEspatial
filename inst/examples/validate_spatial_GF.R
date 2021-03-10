#'
#' A script to validate functionality with GLM and GAM
#' gamma response variable with discrete spatail area plus numeric linear spatial variable fleet and year effects
#' with time-invariante GF 
#'
library(TMB)
library(CPUEspatial)
library(INLA)
library(RandomFields)
library(RANN)
library(ggplot2)
library(mgcv) 
library(gratia) # useful for plotting splines/smoother functions https://github.com/gavinsimpson/gratia
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
n = 1000 # per year
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
year_samples = sample(1:n_years, size = n * n_years, replace = T,  prob = rep(1,n_years) / n_years)

rep_ls = opt_ls = rep_NN_ls = opt_NN_ls = list()
#N_sims = 5
#for(sim in 1:N_sims) {
sim = 1
  if(sim %% 5 == 0)
    cat("sim = ", sim, "\n")
  ## Simulate sample locations based on population spatial distribution
  proj_df$omega <- RFsimulate(omega_model, x=as.matrix(proj_df[,c("x","y")]), exactness=TRUE)$variable1
  
  sampData = fleet_ndx = NULL;
  full_proj_df = NULL
  for(i in 1:n_years) {
    #ggplot(proj_df, aes(x = x, y = y, fill = omega)) + 
    #  geom_tile()
    samples_this_year = sum(year_samples == i)
    x_loc = runif(samples_this_year, 0, grid_dim["x"])
    y_loc = runif(samples_this_year, 0, grid_dim["y"])
    fleet_ndx = c(fleet_ndx, sample(1:length(fleet_coef), size = samples_this_year, replace = T,  prob = c(prob_fleet1_per_year[i], prob_fleet2_per_year[i],prob_fleet3_per_year[i])))
    ##
    sampData = rbind(sampData, data.frame(year = i, x = x_loc, y = y_loc))
    
    ##
    proj_df$year = i
    full_proj_df = rbind(full_proj_df, proj_df)
  }
  sampData$region = NA
  sampData$region[sampData$x < 50 & sampData$y > 50] = 1
  sampData$region[sampData$x > 50 & sampData$y > 50] = 2
  sampData$region[sampData$x < 50 & sampData$y < 50] = 3
  sampData$region[sampData$x > 50 & sampData$y < 50] = 4
  table(sampData$region)
  
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
  
  # tips on debugging
  # track fixed effect parameters
  spatial_model_w_omega$obj$env$tracepar = T
  spatial_model_w_omega_NN$obj$env$tracepar = T
  ## check no zero gradients or NA evaluations
  spatial_model_w_omega$obj$fn()
  spatial_model_w_omega_NN$obj$fn()
  spatial_model_w_omega$obj$gr()
  spatial_model_w_omega_NN$obj$gr()
  
  names(spatial_model_w_omega$obj$par)[spatial_model_w_omega$obj$gr() == 0]
  which(spatial_model_w_omega$obj$gr() == 0)
  spatial_model_w_omega$obj$par
  opt_spa = nlminb(spatial_model_w_omega$obj$par, spatial_model_w_omega$obj$fn, spatial_model_w_omega$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  opt_spa_NN = nlminb(spatial_model_w_omega_NN$obj$par, spatial_model_w_omega_NN$obj$fn, spatial_model_w_omega_NN$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
  rep_omega = spatial_model_w_omega$obj$report()
  rep_omega_NN = spatial_model_w_omega_NN$obj$report()
  
  #sd_rep = sdreport(simple_spatial_model$obj)
  #sd_rep_NN = sdreport(spatial_model_w_omega_NN$obj)
  
  opt_ls[[sim]] = opt_spa
  rep_ls[[sim]] = rep_omega
  
  opt_NN_ls[[sim]] = opt_spa_NN
  rep_NN_ls[[sim]] = rep_omega_NN
#}

est_sd_omega = Reduce(c, lapply(rep_ls, FUN = function(x) {x$MargSD_omega}))
est_sd_omega_NN = Reduce(c, lapply(rep_NN_ls, FUN = function(x) {x$MargSD_omega}))
est_range_omega = Reduce(c, lapply(rep_ls, FUN = function(x) {x$Range_omega}))
est_range_omega_NN = Reduce(c, lapply(rep_NN_ls, FUN = function(x) {x$Range_omega}))

boxplot(cbind(est_sd_omega, est_sd_omega_NN), names = c("Triangulation", "Nearest Neighbour"))
abline(h = sqrt(omega_sigma2), lwd = 3, lty = 2, col = "red")
boxplot(cbind(est_range_omega, est_range_omega_NN), names = c("Triangulation", "Nearest Neighbour"))
abline(h = (omega_range), lwd = 3, lty = 2, col = "red")

glm_fit_spa = glm(y_i ~ factor(year) + factor(fleet_ndx) + factor(region), offset = log(area), data = sampData, family = Gamma(link = log))

## Compare the fits
opt_spa$objective 
logLik(glm_fit_spa)
## get year indes
con_ndx_spa_glm = canonical.index(GLM = glm_fit_spa, year = 1990:1999, base  = 1, year.name = "year")
geo_index = con_ndx_spa_glm$index
plot(geo_index$year, geo_index$index, type = "l", lwd = 3, xlab = "Year", ylab = "Relative Index", ylim = c(0, 4))
lines(geo_index$year, rep_omega$relative_index / gm_mean(rep_omega$relative_index), lwd = 3, lty = 2, col = "red")
lines(geo_index$year, rep_omega_NN$relative_index / gm_mean(rep_omega_NN$relative_index), lwd = 3, lty = 2, col = "blue")
arrows(x0 = geo_index$year, x1 = geo_index$year, y0 = geo_index$lower.CI, 
       y1 = geo_index$upper.CI, lwd = 3, angle = 90, length = 0.05, code = 3)
legend('topright', legend = c("CPUEspatial (NN)", "CPUEspatial (tri)","GLM"), col = c("blue","red","black"), lwd = 3)
## manually calculate the dispersion of Gamma log GLM
1 / summary(glm_fit_spa)$dispersion
rep$phi
AIC = 2*opt_spa$objective + 2*length(simple_spatial_model$obj$par)

## influence plots
#myInfl = Influence$new(glm_fit_spa)
#myInfl$calc()
attr(glm_fit_spa, "family")
attributes(glm_fit_spa)


# Diagnostics
obs_fit_df = data.frame(obs = sampData$y_i, glm = fitted.values(glm_fit_spa), CPUEspatial = rep$mu)
round(head(obs_fit_df),3)

plot(obs_fit_df$obs, obs_fit_df$glm , pch = 16, xlab = "Observed", ylab = "Fitted")
points(obs_fit_df$obs, obs_fit_df$CPUEspatial, pch = 16, col = "red", cex = 0.6)
legend('topright', legend = c("CPUEspatial","GLM"), col = c("red","black"), pch = 16)

## compare estiamted coeffecients
glm_coefs = coefficients(glm_fit_spa)
contrast_region = glm_coefs[grepl(names(glm_coefs), pattern = "region")]
reg_coefs = c(glm_coefs[grepl(names(glm_coefs), pattern = "Intercept")], glm_coefs[grepl(names(glm_coefs), pattern = "Intercept")] + contrast_region)
# transform so comparible to zero sum coeffecients
comparible_coefs = rbind(reg_coefs - mean(reg_coefs)
                         ,rep$spatial_betas[1:(length(rep$spatial_betas) )])

dimnames(comparible_coefs) = list(c("GLM","CPUEspatial"), paste0("region ", 1:4))
comparible_coefs
# depth 
rep$spatial_betas[length(rep$spatial_betas)]
glm_coefs[grepl(names(glm_coefs), pattern = "depth")]

# Year coeffecients
contrast_year = glm_coefs[grepl(names(glm_coefs), pattern = "year")]
year_coefs = c(glm_coefs[grepl(names(glm_coefs), pattern = "Intercept")], glm_coefs[grepl(names(glm_coefs), pattern = "Intercept")] + contrast_year)
comparible_year_coefs = rbind(year_coefs - mean(year_coefs)
                              ,rep$time_betas)
dimnames(comparible_year_coefs) = list(c("GLM","CPUEspatial"), paste0("year ", 1:10))
comparible_year_coefs

# fleet This has the intercept term. GLM has intercept based on reference levels
# where as we have it the mean of factors.
contrast_catch = glm_coefs[grepl(names(glm_coefs), pattern = "fleet_ndx")]
catch_coefs = c(glm_coefs[grepl(names(glm_coefs), pattern = "Intercept")], glm_coefs[grepl(names(glm_coefs), pattern = "Intercept")] + contrast_catch)
## 
cpue_int = catch_coefs[1] + (mean(year_coefs) - year_coefs[1]) + (mean(reg_coefs) - reg_coefs[1])
## all the other catchabilit coeffecients should be the same
comparible_fleet_coefs = rbind(c(cpue_int, contrast_catch)
                               ,rep$betas)

dimnames(comparible_fleet_coefs) = list(c("GLM","CPUEspatial"), paste0("Fleet ", 1:3))
comparible_fleet_coefs

## plot spline effects
sm <- evaluate_smooth(glm_fit_spa, "s(spline_var)")
sm[["upper"]] <- sm[["est"]] + (2 * sm[["se"]])
sm[["lower"]] <- sm[["est"]] - (2 * sm[["se"]])
sm$true = sin(sm$spline_var)
plt <- ggplot(sm, aes_string(x = "spline_var", y = "est")) +
  geom_line(size = 2) +
  geom_ribbon(mapping = aes_string(ymin = "lower", ymax = "upper"), alpha = 0.2) + 
  geom_line(aes(y = true), color="steelblue", size = 2, linetype="twodash") 
plt
DF = data.frame(spline_var = seq(min(sampData[,"spline_var"]), max(sampData[,"spline_var"]), by = diff(range(sampData[,"spline_var"])) / 50));
DF$CPUEspatial = rep$splineForReport_spatial
plt_alt = plt +
  geom_line(data = DF, aes(x = spline_var, y = CPUEspatial), color="red", linetype="twodash", size = 2, inherit.aes = F) 
plt_alt
## TODO add a legend


## look at th GF parameters
rep_omega$MargSD_omega
rep_omega$Range_omega
