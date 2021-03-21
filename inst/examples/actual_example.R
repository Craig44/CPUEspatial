#'
#' Real example
#'

#source("influ.R")
library(CPUEspatial)
library(kableExtra)
library(nzPlot)
library(geoR)
library(rgdal)
library(rgeos)
library(fields)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(raster)
library(INLA)
source(file.path("inst", "examples","read_NZ_polys.R"))
head(cpue_df)
cpue_df$area = cpue_df$fishing_distance * cpue_df$effort_width/1000 # km^2
cpue_df$log_catch = log(cpue_df$BAR_catch)
bad_ndx = names(table(cpue_df$vessel_key))[table(cpue_df$vessel_key) < 5]
dim(cpue_df)
cpue_df = subset(cpue_df, subset = !cpue_df$vessel_key %in% bad_ndx )
dim(cpue_df)
## fit glm model
log_positive = glm(log_catch ~ fish_year + vessel_key + target_species + start_stats_area_code + fish_month, offset = area, data = cpue_df)
log_positive_alt = glm(BAR_catch ~ fish_year + vessel_key + target_species + start_stats_area_code +fish_month, offset = area, data = cpue_df, family = gaussian(link = "log"))
gamma_glm = glm(BAR_catch ~ fish_year + vessel_key + target_species + start_stats_area_code + fish_month, offset = area, data = cpue_df, family = Gamma(link = "log"))
summary(log_positive)
summary(gamma_glm)
summary(log_positive_alt)

##Geostatistical model set up
nz(ylim = c(-45, -53))
nz.polygon(nz.statarea(labels = T))
nz.points(cpue_df$start_longitude, cpue_df$start_latitude, pch = 16)
#proj_poly = nz.locator()
proj_poly =  list(x = c(165.4688, 165.6495, 165.6947, 167.9993, 167.9993, 168.9935, 170.1684, 170.1684, 168.2253),
                  y = c(-46.32088, -50.02528, -52.04248, -51.37043, -48.99842, -49.02806, -48.52180, -46.66310, -46.03929))
## 503 = 29
## 603 = 602
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
  df = data.frame(x = coordinates(proj_utm_spdf)[,1], y = coordinates(proj_utm_spdf)[,1], start_stats_area_code = proj_utm_spdf$stat_area, area =  rgeos::gArea(sp_poly_grid), fish_year = as.character(years[i]))
  full_proj_df = rbind(full_proj_df, df)
}
dim(full_proj_df)
full_proj_df$log_catch = 1
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


data = data2use2_utm
projection_df = full_proj_df
mesh = mesh_utm
family = 0
link = 4
include_omega = F
include_epsilon = F 
response_variable_label = "log_catch"
time_variable_label = "fish_year"
catchability_covariates = c("vessel_key","target_species","fish_month")
catchability_covariate_type = c("factor", "factor","factor")
spatial_covariates = c("start_stats_area_code")
spatial_covariate_type = c("factor")
spline_catchability_covariates = NULL
spline_spatial_covariates = NULL
trace_level = "high"
apply_preferential_sampling = T
preference_model_type = 1
projection_raster_layer = proj_raster_with_active_cells

simple_glm = configure_obj(data = data2use2_utm, projection_df = full_proj_df, mesh = mesh_utm, family = 0, link = 4, include_omega = F, include_epsilon = F, 
                           response_variable_label = "log_catch", time_variable_label = "fish_year", catchability_covariates = c("vessel_key","target_species","fish_month"), catchability_covariate_type = c("factor", "factor","factor"), 
                           spatial_covariates = c("start_stats_area_code"), spatial_covariate_type = c("factor"), spline_catchability_covariates = NULL,
                           spline_spatial_covariates = NULL, trace_level = "high")

simple_glm$obj$fn()

