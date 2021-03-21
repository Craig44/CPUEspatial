## Convert nzPlot objects to spatialDF's so we can use multiple projections 
## for coast polygons, depth lines
## nz objects are either vector or data.frames that create polygons when there is a row of NA's
library(nzPlot)
library(sp)
library(ggplot2)
## convert lat long to UTM
unproj_args = " +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 +lon_wrap=180"
unproj <- CRS(unproj_args) #default WGS84 projection
zone = 60 ## change this based on data, for NZ good enough
proj <- CRS(paste0("+proj=utm +zone=", zone, "+datum=WGS84 +units=km"))  
## create a polygon class for stat areas
nz_stat_area_polys = list()
areas = unique(nz.statarea.lines$area)
for(i in 1:length(areas)) {
  ndx = nz.statarea.lines$area == areas[i]
  nz.statarea.lines$area[ndx]
  coords = cbind(nz.statarea.lines$long[ndx], nz.statarea.lines$lat[ndx])
  second_ndx = complete.cases(coords)
  nz_stat_area_polys[[i]] = Polygons(list(Polygon(coords[second_ndx,])), ID = as.character(sprintf("%03d", areas[i])))
}
nz.stat_area = SpatialPolygons(Srl = nz_stat_area_polys, proj4string = unproj)
nz.stat_area_utm = spTransform(nz.stat_area, proj)

#nz.stat_area_df = SpatialPolygonsDataFrame(nz.stat_area)

poly_ndx = c(0, which(is.na(nz.coast$long)), nrow(nz.coast) + 1) # the first and last don't have NA's so just add them
nz_coast_polys = list()
for(i in 2:length(poly_ndx)) {
  nz_coast_polys[[i - 1]] = Polygon(as.matrix(nz.coast[(poly_ndx[i - 1] + 1):(poly_ndx[i] - 1),]))
}
nz.coast_sp = SpatialPolygons(Srl = list(Polygons(nz_coast_polys, "coast")), proj4string = unproj)
nz.coast_utm_sp = spTransform(nz.coast_sp, proj)

# in order to plot polygons, first fortify the data
# create a data.frame from our spatial object
nz.coast_sp_df <- fortify(nz.coast_sp, region = "coast")
nz.coast_utm_sp = spTransform(nz.coast_sp, proj)
nz.coast_sp_utm_df <- fortify(nz.coast_utm_sp, region = "coast")
## do the same but with nz.islands
head(which(is.na(nz.islands[["long"]])))
head(which(is.na(nz.islands[["lat"]])))
# A nuiance below, the first entry is NA unlike the above case see: head(nz.islands[["long"]])
poly_ndx = c(which(is.na(nz.islands[["long"]])), nrow(nz.islands) + 1) # the first and last don't have NA's so just add them
nz_island_polys = list()
for(i in 2:length(poly_ndx)) {
  nz_island_polys[[i - 1]] = Polygon(cbind(nz.islands[["long"]][(poly_ndx[i - 1] + 1):(poly_ndx[i] - 1)], nz.islands[["lat"]][(poly_ndx[i - 1] + 1):(poly_ndx[i] - 1)]))
}
nz.island_sp = SpatialPolygons(Srl = list(Polygons(nz_island_polys, "coast")), proj4string = unproj)
nz.island_sp_df <- fortify(nz.island_sp, region = "coast")
nz.island_utm_sp = spTransform(nz.island_sp, proj)
nz.island_sp_utm_df <- fortify(nz.island_utm_sp, region = "coast")

# in order to plot polygons, first fortify the data
# create a data.frame from our spatial object

## Do it for lines
## Do it for lines
line_ndx = c(which(is.na(nzPlot::nz.200m[,"long"])), nrow(nzPlot::nz.200m) + 1) # the first and last don't have NA's so just add them
nz_200 = list()
for(i in 2:length(line_ndx)) {
  nz_200[[i - 1]] = Line(nzPlot::nz.200m[(line_ndx[i - 1] + 1):(line_ndx[i] - 1),])
}
unproj = CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 +lon_wrap=180")
nz_200_sp = SpatialLines(LinesList = list(Lines(nz_200, "200m")), proj4string = unproj)
nz_200_sp_df = SpatialLinesDataFrame(nz_200_sp, data = data.frame(z = 1), match.ID = F)
nz_200_utm_sp = spTransform(nz_200_sp_df, proj)
nz_200_sp_df = fortify(nz_200_sp_df, region = "200m")
nz_200_sp_utm_df <- fortify(nz_200_utm_sp, region = "200m")


line_ndx = c(which(is.na(nzPlot::nz.150m[,"long"])), nrow(nzPlot::nz.150m) + 1) # the first and last don't have NA's so just add them
nz_150 = list()
for(i in 2:length(line_ndx)) {
  nz_150[[i - 1]] = Line(nzPlot::nz.150m[(line_ndx[i - 1] + 1):(line_ndx[i] - 1),])
}
nz_150_sp = SpatialLines(LinesList = list(Lines(nz_150, "150m")), proj4string = unproj)
nz_150_sp_df = SpatialLinesDataFrame(nz_150_sp, data = data.frame(z = 1), match.ID = F)
nz_150_utm_sp = spTransform(nz_150_sp_df, proj)
nz_150_sp_df = fortify(nz_150_sp_df, region = "150m")
nz_150_sp_utm_df <- fortify(nz_150_utm_sp, region = "150m")

line_ndx = c(which(is.na(nzPlot::nz.250m[,"long"])), nrow(nzPlot::nz.250m) + 1) # the first and last don't have NA's so just add them
nz_250 = list()
for(i in 2:length(line_ndx)) {
  nz_250[[i - 1]] = Line(nzPlot::nz.250m[(line_ndx[i - 1] + 1):(line_ndx[i] - 1),])
}
nz_250_sp = SpatialLines(LinesList = list(Lines(nz_250, "250m")), proj4string = unproj)
nz_250_sp_df = SpatialLinesDataFrame(nz_250_sp, data = data.frame(z = 1), match.ID = F)
nz_250_utm_sp = spTransform(nz_250_sp_df, proj)
nz_250_sp_df = fortify(nz_250_sp_df, region = "250m")
nz_250_sp_utm_df <- fortify(nz_250_utm_sp, region = "250m")


line_ndx = c(which(is.na(nzPlot::nz.100m[,"long"])), nrow(nzPlot::nz.100m) + 1) # the first and last don't have NA's so just add them
nz_100 = list()
for(i in 2:length(line_ndx)) {
  nz_100[[i - 1]] = Line(nzPlot::nz.100m[(line_ndx[i - 1] + 1):(line_ndx[i] - 1),])
}
nz_100_sp = SpatialLines(LinesList = list(Lines(nz_100, "100m")), proj4string = unproj)
nz_100_sp_df = SpatialLinesDataFrame(nz_100_sp, data = data.frame(z = 1), match.ID = F)
nz_100_utm_sp = spTransform(nz_100_sp_df, proj)
nz_100_sp_df = fortify(nz_100_sp_df, region = "100m")
nz_100_sp_utm_df <- fortify(nz_100_utm_sp, region = "100m")


line_ndx = c(which(is.na(nzPlot::nz.500m[,"long"])), nrow(nzPlot::nz.500m) + 1) # the first and last don't have NA's so just add them
nz_500 = list()
for(i in 2:length(line_ndx)) {
  nz_500[[i - 1]] = Line(nzPlot::nz.500m[(line_ndx[i - 1] + 1):(line_ndx[i] - 1),])
}
nz_500_sp = SpatialLines(LinesList = list(Lines(nz_500, "500m")), proj4string = unproj)
nz_500_sp_df = SpatialLinesDataFrame(nz_500_sp, data = data.frame(z = 1), match.ID = F)
nz_500_utm_sp = spTransform(nz_500_sp_df, proj)
nz_500_sp_df = fortify(nz_500_sp_df, region = "500m")
nz_500_sp_utm_df <- fortify(nz_500_utm_sp, region = "500m")
