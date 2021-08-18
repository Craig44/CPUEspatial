#'
#' Real example loosely based on BAR 5 & 6
 # 603 = 602
## 604 = 602
## 303 = 504

nz.polygon(proj_poly, border = "blue")
## remove observations outside of this polugon
loc_ndx = point.in.polygon(point.x = cpue_df$start_longitude, point.y = cpue_df$start_latitude, pol.x = proj_poly$x, pol.y = proj_poly$y)
data2use2 = subset(x = cpue_df, subset = loc_ndx == 1)
table(loc_ndx)

# convert factors to characters
data2use2$fish_year = as.integer(as.character(data2use2$fish_year))
data2use2$fish_month = (as.character(data2use2$fish_month))
data2use2$primary_method = (as.character(data2use2$primary_method))
data2use2$vessel_key = (as.character(data2use2$vessel_key))
data2use2$target_species = (as.character(data2use2$target_species))
data2use2$start_stats_area_code = (as.character(data2use2$start_stats_area_code))

# convert to Spatial point DF
coordinates(data2use2) <- ~ start_longitude + start_latitude
proj4string(data2use2) <- unproj
# convert to UTM, better for distance based analysis than Degrees
data2use2_utm <- spTransform(data2use2, proj)
## land polygons nz.coast_sp
## Build Mesh with islands in it
plot(data2use2_utm$start_longitude, y=data2use2_utm$start_latitude, pch=20, asp=1)
plot(nz.coast_utm_sp, add=T, border="black", col="grey")

bndary_utm = inla.nonconvex.hull(cbind(data2use2_utm$start_longitude, data2use2_utm$start_latitude), resolution = 50, convex = -0.10)
max.edge = 35 # km
mesh_utm = inla.mesh.2d(loc=cbind(data2use2_utm$start_longitude, data2use2_utm$start_latitude),
                        max.edge = c(1) * max.edge,
                        boundary = bndary_utm,
                        # - use 5 times max.edge in the outer extension/offset/boundary
                        cutoff = max.edge/3)
## visualise for UTM
par(mar = c(1,1,1,1))
plot(mesh_utm, main= "With Boundary"); 
points(data2use2_utm$start_longitude, data2use2_utm$start_latitude, col="red", cex = 0.4)
plot(mesh_utm, add = T); 
dim(mesh_utm$loc)

### barier mesh
## two sources
# https://github.com/skaug/tmb-case-studies/blob/master/spdeBarrier/spdeBarrier.cpp
# https://haakonbakkagit.github.io/btopic107.html

tl = length(mesh$graph$tv[,1])
# - the number of triangles in the mesh
posTri = matrix(0, tl, 2)
for (t in 1:tl){
  temp = mesh$loc[mesh$graph$tv[t, ], ]
  posTri[t,] = colMeans(temp)[c(1,2)] 
}
posTri = SpatialPoints(posTri)

proj4string(posTri) <- proj

# - compute the triangle positions
normal = over(nz.coast_utm_sp, posTri, returnList=T)
# - checking which mesh triangles are inside the normal area
normal = unlist(normal)
barrier.triangles = setdiff(1:tl, normal)


fem = INLA:::inla.barrier.fem(mesh = mesh, barrier.triangles = barrier.triangles)
spdeMatricesBarrier = list(C0 = fem$C[[1]],C1 =fem$C[[2]] ,D0 = fem$D[[1]],D1 = fem$D[[2]],I = fem$I )
barrier = 1
c = c(1,0.2)
#poly.barrier = inla.barrier.polygon(mesh, barrier.triangles)
#barrier.model = inla.barrier.pcmatern(mesh, barrier.triangles = barrier.triangles, prior.range = c(6, .5), prior.sigma = c(3, 0.01))


########################
## Exploratory analysis
########################
## rasterise our variabels to this resolttion
ylim_utm = range(data2use2_utm$start_latitude)
xlim_utm = range(data2use2_utm$start_longitude)
ylim_dec = range(data2use2$start_latitude)
xlim_dec = range(data2use2$start_longitude)

#train_df_utm$catch_rate = train_df_utm$SCI / train_df_utm$area_swept
dec_raster <- raster(resolution = 5, xmn = xlim_utm[1], xmx = xlim_utm[2], ymn = ylim_utm[1], ymx = ylim_utm[2], crs = proj)
dec_raster_dec <- raster(resolution = 0.2, xmn = xlim_dec[1], xmx = xlim_dec[2], ymn = ylim_dec[1], ymx = ylim_dec[2], crs = unproj)
## aggregrate over all years
mean_density <- rasterize(data2use2, dec_raster_dec, field = data2use2$BAR_catch, median)
total <- rasterize(data2use2, dec_raster_dec, field = data2use2$BAR_catch, sum)
n_tows <- rasterize(data2use2, dec_raster_dec, field = data2use2$BAR_catch, fun = function(x, na.rm = T){length(as.vector(x))})
year_plot = data.frame(x = coordinates(mean_density)[,1], y = coordinates(mean_density)[,2], mean_density = mean_density@data@values)
year_plot$tows = n_tows@data@values
year_plot$Sum = total@data@values

###################
## Create Projection data frame
###################
# define lattice discretisation for TMB
# these are mid points
Y_proj <- seq(min(bndary_utm$loc[, 2]),max(bndary_utm$loc[, 2]),length.out=100)
X_proj <- seq(min(bndary_utm$loc[, 1]),max(bndary_utm$loc[, 1]),length.out=100)
x_proj_res = unique(round(diff(X_proj),2))
y_proj_res = unique(round(diff(Y_proj),2))
## create a raster layer. more for jointly modelling sample locations, will be ignored for most models
proj_grid <- expand.grid(X_proj, Y_proj)
x_temp = matrix(proj_grid$Var1, nrow = length(X_proj), ncol = length(Y_proj), byrow = T)
x_temp[1:5,1:5]
y_temp = matrix(proj_grid$Var2, nrow = length(X_proj), ncol = length(Y_proj), byrow = T)
y_temp[1:5,1:5]

loc_ndx = point.in.polygon(point.x = proj_grid[,1], point.y = proj_grid[,2], 
                           pol.x = bndary_utm$loc[,1], pol.y = bndary_utm$loc[, 2])
## if an irregular projection grid then need to supply a raster layer
active_proj_cell = matrix(loc_ndx, nrow = length(X_proj), ncol = length(Y_proj), byrow = T)
active_proj_cell = active_proj_cell[nrow(active_proj_cell):1, ]

proj_raster_with_active_cells = raster(x = (active_proj_cell), xmn = min(X_proj) - x_proj_res/2, xmx = max(X_proj) + x_proj_res/2, ymn = min(Y_proj) - y_proj_res/2, ymx = max(Y_proj) + y_proj_res/2, crs = proj)
#proj_raster_with_active_cells = raster(x = active_proj_cell, xmn = min(X_proj), xmx = max(X_proj), ymn = min(Y_proj), ymx = max(Y_proj), crs = proj)
## check 
par(mfrow = c(1,1))
plot(proj_raster_with_active_cells)
plot(mesh_utm, main= "With Boundary", add = T); 
points(data2use2_utm, pch = 16, col = "red")


data2use2_utm$indicator = 1
proj_count <- rasterize(data2use2_utm, proj_raster_with_active_cells, field = data2use2_utm$indicator, sum, na.rm = T)
nrow(data2use2_utm)
sum(proj_count$layer@data@values[proj_raster_with_active_cells@data@values == 1], na.rm = T)


## not compatible
plot(proj_raster_with_active_cells)
points(data2use2_utm, pch = 16, col = "red")

proj_grid = proj_grid[loc_ndx == 1, ]


head(proj_grid)

nrow(proj_grid) # number of projection cells
## get the size of a single cell for area projections
approX_proj_grid = cbind(x = c(X_proj[1], X_proj[1], X_proj[2], X_proj[2]),
                         y = c(Y_proj[1], Y_proj[2], Y_proj[2], Y_proj[1]))
sp_poly_grid = SpatialPolygons(list(Polygons(list(Polygon(coords = coordinates(approX_proj_grid))) ,ID = "grid")), proj4string = proj)
## convert to meters
raster::crs(sp_poly_grid)
sp_poly_grid = spTransform(sp_poly_grid, CRSobj = proj)
raster::crs(sp_poly_grid)
grid_area = rgeos::gArea(sp_poly_grid) #* 1e-6 # convert from m^2 * 1e-6 = km^2, changed to 9 as it is a constant
grid_area 

## convert utm - > dec so we can attach bathymetry
boundary_dec = SpatialPolygons(list(Polygons(list(Polygon(coords = bndary_utm$loc)) ,ID = "bndry")), proj4string = proj)
boundary_dec = spTransform(boundary_dec, CRSobj = unproj)

coordinates(proj_grid) <- ~ Var1 + Var2
proj4string(proj_grid) <- proj

## assign start_stats_area_code
proj_utm_spdf = proj_grid

ndx = over(proj_utm_spdf, nz.stat_area_utm)
proj_utm_spdf$stat_area = names(nz.stat_area_utm)[ndx]
plot(proj_utm_spdf)
points(proj_utm_spdf[is.na(ndx),] , col = "red")
## land!!! associate closest stat area

## Set up containers for results
pts = proj_utm_spdf[is.na(ndx),]
n <- length(pts)
nearest_stat_area <- character(n)
distToNearestCanton <- numeric(n)
## For each point, find name of nearest polygon (in this case, Belgian cantons)
for (i in seq_along(nearest_stat_area)) {
  gDists <- gDistance(pts[i,], nz.stat_area_utm, byid=TRUE)
  nearest_stat_area[i] <- names(nz.stat_area_utm)[which.min(gDists)]
  distToNearestCanton[i] <- min(gDists)
}
## visually inspect this
plot(nz.stat_area_utm)
table(is.na(proj_utm_spdf$stat_area))
proj_utm_spdf$stat_area[is.na(ndx)] = nearest_stat_area

points(proj_utm_spdf, col = proj_utm_spdf$stat_area, pch = 16, cex = 0.3)
plot(nz.stat_area_utm, add = T)
points(data2use2_utm, pch = 16, cex = 0.1)
## no spatial temporal  varying covariates so, we
## can just append with a year variable to get full proj_df
years = unique(data2use2$fish_year)
full_proj_df = NULL
for(i in 1:length(years)) {
  df = data.frame(x = coordinates(proj_utm_spdf)[,1], y = coordinates(proj_utm_spdf)[,2], start_stats_area_code = proj_utm_spdf$stat_area, area =  rgeos::gArea(sp_poly_grid), fish_year = as.character(years[i]))
  full_proj_df = rbind(full_proj_df, df)
}
dim(full_proj_df)
full_proj_df$log_catch = 1
full_proj_df$BAR_catch = 1

coordinates(full_proj_df) <- ~ x + y
proj4string(full_proj_df) <- proj

## now we can fit the model I think
data = data2use2_utm
projection_df = full_proj_df

## check there are not levels of factors in projection grid that
## are not in observed data.
table(full_proj_df$start_stats_area_code)
table(data2use2_utm$start_stats_area_code)
table(is.na(data2use2_utm$start_stats_area_code))
areas_not_in_data = !names(table(full_proj_df$start_stats_area_code)) %in% names(table(data2use2_utm$start_stats_area_code))
names(table(full_proj_df$start_stats_area_code))[areas_not_in_data]
## 503 = 29
## 603 = 602
## 604 = 602
## 303 = 504
## 26 = 35
neighbour_match = data.frame(area = (c("503", "603", "604", "303", "026")), neighbour_area = (c("029", "602", "602", "504", "030")))
ndx = match(full_proj_df$start_stats_area_code, table = neighbour_match$area)
full_proj_df$start_stats_area_code[which(!is.na(ndx))] = neighbour_match$neighbour_area[ndx[!is.na(ndx)]]
areas_not_in_data = !names(table(full_proj_df$start_stats_area_code)) %in% names(table(data2use2_utm$start_stats_area_code))
names(table(full_proj_df$start_stats_area_code))[areas_not_in_data]
table(is.na(full_proj_df$start_stats_area_code))
na_ndx = which(is.na(full_proj_df$start_stats_area_code))
plot(nz.stat_area_utm)
points(coordinates(full_proj_df[na_ndx,]), col = "red", cex =2)


## deal with covaraiate types
data2use2_utm$fish_year = as.integer(data2use2_utm$fish_year)
data2use2_utm$start_stats_area_code = as.factor(data2use2_utm$start_stats_area_code)
data2use2_utm$fish_month = factor(data2use2_utm$fish_month, levels = month.abb)
data2use2_utm$vessel_key = as.factor(data2use2_utm$vessel_key)
data2use2_utm$target_species = as.factor(data2use2_utm$target_species)

non_spatial = configure_obj(data = data2use2_utm, projection_df = full_proj_df, mesh = mesh_utm, family = 0, link = 0, include_omega = F, include_epsilon = F, 
                           response_variable_label = "BAR_catch", time_variable_label = "fish_year", catchability_covariates = c("vessel_key","target_species","fish_month"), 
                           spatial_covariates = c("start_stats_area_code"), spline_catchability_covariates = NULL,
                           spline_spatial_covariates = NULL, trace_level = "high")

check_gradients(non_spatial$obj)


opt = nlminb(non_spatial$obj$par, non_spatial$obj$fn, non_spatial$obj$gr, control = list(eval.max = 10000, iter.max = 10000))
opt$convergence

sd_rep = sdreport(non_spatial$obj)

non_spatial_rep = non_spatial$obj$report(non_spatial$obj$env$last.par.best)
## glm approach
glm_data = data2use2_utm@data
glm_data$fish_year = as.factor(glm_data$fish_year)
#glm_data$vessel_key
log_positive = glm((BAR_catch) ~ fish_year + vessel_key + target_species + start_stats_area_code +fish_month, offset = log(area), data = glm_data, family = gaussian(link = "log"))
summary(log_positive)$dispersion
non_spatial_rep$phi
## check they give the same estimates
logLik(log_positive)
opt$objective

par_labs = names(sd_rep$value)

sd_rep$sd[par_labs %in% "betas"]
sd_rep$sd[par_labs %in% "betas_w_intercept"]
sd_rep$value[par_labs %in% "betas"]
sd_rep$value[par_labs %in% "betas_w_intercept"]


## compare estiamted coeffecients
term = "fish_month"
glm_coefs = coefficients(log_positive)
contrast_month = glm_coefs[grepl(names(glm_coefs), pattern = term)]

catchability_mod_mat = non_spatial$tmb_data$model_matrix
catchability_coef_labs = colnames(catchability_mod_mat)
catchability_intercept = non_spatial_rep$betas[1]
mod_ndx = grepl(catchability_coef_labs, pattern = term)
catchability_coeffs = non_spatial_rep$betas[mod_ndx]

labs = substring(catchability_coef_labs[mod_ndx], first = nchar(term) + 1)


ff = formula(non_spatial$Call$catchability)
factor_lab = get_all_vars(ff, data = data2use2_utm@data, drop.unused.levels = T)
factor_lab$time = data2use2_utm@data[, "fish_year"]
data_labs = unique(as.character(factor_lab[,term]))
intercept = data_labs[!data_labs %in% labs]
labs = c(intercept, labs)
mod_ndx[1] = T
coeffs_ = data.frame(lab = labs, MLE = sd_rep$value[par_labs %in% "betas_w_intercept"][mod_ndx],
                     SE = sd_rep$sd[par_labs %in% "betas_w_intercept"][mod_ndx])
coeffs_$upper = coeffs_$MLE + 2* coeffs_$SE
coeffs_$lower = coeffs_$MLE - 2* coeffs_$SE
coeffs_$lab = factor(coeffs_$lab,levels=month.abb,ordered=T)

contrast_month
catchability_coeffs
# looks good
##################################################
## Influence Plots calculations
##################################################
myInfl = Influence$new(model = log_positive)
myInfl$calc()
myInfl$summary

myInfl_cpue = calculate_influence(non_spatial, data2use2_utm)
myInfl$summary$overall
myInfl_cpue$overall_influence

plot_influence(myInfl_cpue, "fish_month")
myInfl$cdiPlot('fish_month')

plot_influence(myInfl_cpue, "target_species")
myInfl$cdiPlot('target_species')

plot_influence(myInfl_cpue, "vessel_key")
myInfl$cdiPlot('vessel_key')

plot_influence(myInfl_cpue, "start_stats_area_code")
myInfl$cdiPlot('start_stats_area_code')



tab = kable(myInfl$summary,
            row.names = FALSE,
            caption = "Influence metrics.",
            format = "latex", escape=F, digits = 2, booktabs = T)
## if you want the float [H] so table stays in a latex document use latex_options = "HOLD_position"
tab = tab %>%
  kable_styling(position = "center", latex_options = c("HOLD_position", "scale_down"))

print(tab) # copy this into the latex document
#save to file called "car_table.tex"
#cat(tab, file = file.path("..","writing", "influence.tex"))

#png(filename = file.path(DIR$fig, "stan_vs_unstan.png"), res = 150, height = 7, width = 10, units = "in")
myInfl$stanPlot()
#dev.off()



labs = names(sd_rep$value)

## get levels
## if numeric cut() into break factors
## Reorder levels according to coefficients if necessary
## distrs = aggregate(.$preds[,1],list(levels,.$preds[,.$focus]),length)
## 
## 
. = myInfl
term = "fish_month"
variable=NULL


#Define levels of term on which coefficient and distribution plots will be based
#This is done by search for each column name in the data as a whole word in the
#each term. This allows for matching of terms like 'poly(log(var),3)' with 'var'
if(is.null(variable)){
  for(name in names(.$preds)){
    match = grep(paste('([^[:alnum:]_])+',name,'([^[:alnum:]_])+',sep=''),paste('(',term,')')) # ([^\\w])+ Means ignore any 'word' character (alphanumerics plus underscore)
    if(length(match)>0){
      variable = name
      break
    }
  }
}
if(is.null(variable)) stop('Unable to find a matching variable for term "',term,'". Please specify in the argument "variable".')
levels = .$preds[,variable]

#Numeric terms are cut into factors
if(is.numeric(levels)){
  breaks = pretty(levels,30)
  step = breaks[2]-breaks[1]
  labels = breaks+step/2
  breaks = c(breaks,breaks[length(breaks)]+step)
  levels = cut(levels,breaks,labels=labels,include.lowest=T)
}
#Reorder levels according to coefficients if necessary
if(.$orders[[term]]=='coef'){
  coeffs = aggregate(.$preds[,paste('fit',term,sep='.')],list(levels),mean)
  names(coeffs) = c('term','coeff')
  coeffs = coeffs[order(coeffs$coeff),]
  levels = factor(levels,levels=coeffs$term,ordered=T)
}


## get the terms equivalant from CPUEspatial
catchability_mod_mat = non_spatial$tmb_data$model_matrix
catchability_coef_labs = colnames(catchability_mod_mat)
mod_ndx = grepl(catchability_coef_labs, pattern = term)
mod_ndx[1] = TRUE

catchability_intercept = non_spatial_rep$betas[1]
catchability_coeffs = non_spatial_rep$betas[mod_ndx]
spatial_mod_mat = non_spatial$tmb_data$X_spatial_ipt[,,1]
catch_terms = sweep(catchability_mod_mat ,MARGIN=2,non_spatial_rep$betas,`*`)  

spatial_terms = sweep(spatial_mod_mat ,MARGIN=2,non_spatial_rep$spatial_betas,`*`)  

this_term_cpue = c(catchability_intercept, catchability_intercept + catchability_coeffs)
names(this_term_cpue) = c(catchability_coef_labs[1], catchability_coef_labs[mod_ndx])

## get the GLM equivalent
raw_coeffs = coef(log_positive)
glm_raw_coeffs = grepl(names(raw_coeffs), pattern = term)
raw_coeffs[glm_raw_coeffs]

this_term_glm = c(raw_coeffs[1], raw_coeffs[1] + raw_coeffs[glm_raw_coeffs])

this_term_glm - mean(this_term_glm)
this_term_cpue - mean(this_term_cpue)

plot(this_term_glm)
points(this_term_cpue, col = "red")

# Estiamted Coefficients
coeffs = aggregate(.$preds[,paste(c('fit','se.fit'),term,sep='.')],list(levels),mean)
this_term_glm
## the same to a constant
coeffs$fit.fish_month - min(coeffs$fit.fish_month)
this_term_glm - min(this_term_glm)





names(coeffs) = c('term','coeff','se')
coeffs = within(coeffs,{
  lower = coeff-se
  upper = coeff+se
})

par(mar=c(0,5,3,0),las=1)
with(coeffs,{
  xs = 1:max(as.integer(term))
  ylim = c(min(exp(lower)),max(exp(upper)))
  if(ylim[1]<0.5*min(exp(coeff))) ylim[1] = 0.5*min(exp(coeff))
  if(ylim[2]>2*max(exp(coeff))) ylim[2] = 2*max(exp(coeff))
  plot(as.integer(term),exp(coeff),xlim=range(xs),ylim=ylim,pch=2,cex=1.5,xaxt='n',ylab='',log='y')
  mtext('Coefficient',side=2,line=4,las=0,cex=0.8)
  abline(h=1,lty=2)
  abline(v=xs,lty=1,col='grey')
  segments(as.integer(term),exp(lower),as.integer(term),exp(upper))
  arrows(as.integer(term),exp(lower),as.integer(term),exp(upper),angle=90,code=3,length=0.05)
  axis(3,at=xs,labels=levels(term)[xs])
})


par(mar=c(0,5,3,0),las=1)
with(coeffs_,{
  xs = 1:max(as.integer(lab))
  ylim = c(min((lower)),max((upper)))
  if(ylim[1]<0.5*min(( MLE))) ylim[1] = 0.5*min(( MLE))
  if(ylim[2]>2*max(( MLE))) ylim[2] = 2*max(( MLE))
  plot(as.integer(lab), (MLE), xlim = range(xs),ylim=ylim,pch=2,cex=1.5,xaxt='n',ylab='',log='y')
  mtext('Coefficients',side=2,line=4,las=0,cex=0.8)
  abline(h=1,lty=2)
  abline(v=xs,lty=1,col='grey')
  segments(as.integer(lab),(lower),as.integer(lab),(upper))
  arrows(as.integer(lab),(lower),as.integer(lab),(upper),angle=90,code=3,length=0.05)
  axis(3,at=xs,labels=levels(lab))
})

## number of observatons over time among the levels in a variable
distr = factor_lab %>% 
  group_by(get(term), time) %>%
  summarise(n = length(time))
distr = distr %>% group_by(time) %>% mutate(n_total = sum(n))
distr$n_prop = distr$n / distr$n_total
colnames(distr) = c("term", non_spatial$Call$func_call$time_variable_label, "n","n_total","n_prop")
distr$fish_year = factor(distr$fish_year, levels = min(distr$fish_year):max(distr$fish_year))
#Distribution
distrs = aggregate(.$preds[,1],list(levels,.$preds[,.$focus]),length)
names(distrs) = c('term','focus','count')
distrs = merge(distrs,aggregate(list(total=distrs$count),list(focus=distrs$focus),sum))
distrs$prop = with(distrs,count/total)
par(mar=c(5,5,0,0),las=1)
xlab = .$labels[[variable]]

table(distrs$total)
table(distr$n_total)


par(mar=c(5,5,0,0),las=1)
xlab = term
if(is.null(xlab)) xlab = variable
ylab= .$labels[[.$focus]]
if(is.null(ylab)) ylab = .$focus
with(distrs,{
  xs = 1:max(as.integer(term))
  ys = 1:max(as.integer(focus))
  plot(NA,xlim=range(xs),ylim=range(ys),xaxt='n',yaxt='n',xlab=xlab,ylab='')
  abline(v=xs,lty=1,col='grey')
  axis(1,at=xs,labels=levels(term)[xs])
  abline(h=ys,lty=1,col='grey')
  axis(2,at=ys,labels=levels(focus)[ys])
  mtext(ylab,side=2,line=4,las=0,cex=0.8)
  points(as.integer(term),as.integer(focus),cex=sqrt(prop)*12)
})

## distribution of events (observatons) not catch
par(mar=c(5,5,0,0),las=1)
xlab = term
if(is.null(xlab)) xlab = variable
ylab= non_spatial$Call$func_call$time_variable_label
if(is.null(ylab)) ylab = .$focus
with(distr,{
  xs = 1:max(as.integer(term))
  ys = 1:max(as.integer(get(ylab)))
  plot(NA,xlim=range(xs),ylim=range(ys),xaxt='n',yaxt='n',xlab=xlab,ylab='')
  abline(v=xs,lty=1,col='grey')
  axis(1,at=xs,labels=levels(term)[xs])
  abline(h=ys,lty=1,col='grey')
  axis(2,at=ys,labels=levels(get(ylab))[ys])
  mtext(ylab,side=2,line=4,las=0,cex=0.8)
  points(as.integer(term),as.integer(get(ylab)),cex=sqrt(n_prop)*12)
})



#Calculate influences and statisitcs
.$influences = data.frame(level=levels(.$model$model[,.$focus]))
overall = c(NA,NA) # NAs for null model and for focus term
trend = c(NA,NA)
for(term in .$terms){
  if(term != .$focus){
    infl = aggregate(
      list(value = .$preds[,paste('fit',term,sep='.')]),
      list(level = .$preds[,.$focus]),
      mean
    )
    overall = c(overall,with(infl,exp(mean(abs(value)))-1))
    trend = c(trend,with(infl,exp(cov(1:length(value),value)/var(1:length(value)))-1))
    names(infl) = c('level',term)
    .$influences = merge(.$influences,infl,all.x=T,by='level')
  }
}



