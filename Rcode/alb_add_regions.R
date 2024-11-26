# Add regional component
# Modify a_el in Extrapolation_list with regional assignments

add.region <- function(data, Extrapolation_List, strata.limits, grid.df){

library(sp)

load("E:/ALB_CPUE/2024/Background_Data/alb_regions_poly_2024_shp.Rdata")

r1.coords = regions.shp@polygons[[1]]@Polygons[[1]]@coords
r2.coords = regions.shp@polygons[[2]]@Polygons[[1]]@coords

extrap.cells = as.data.frame(Extrapolation_List$Data_Extrap)
coordinates(extrap.cells) = ~Lon + Lat
proj4string(extrap.cells) = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

extrap.cells$r1 <- point.in.polygon(extrap.cells$Lon, extrap.cells$Lat, r1.coords[,1], r1.coords[,2])
extrap.cells$r2 <- point.in.polygon(extrap.cells$Lon, extrap.cells$Lat, r2.coords[,1], r2.coords[,2])

# Create a_el for Extrapolation_List 
# Is cell is in region - assign 1, otherwise 0
a_el = as.data.frame(extrap.cells[,c('r1','r2')])[,1:2]
a_el %<>% mutate(r1 = ifelse(r1==0,0,1), r2 = ifelse(r2==0,0,1))

# Area of each cell in extrapolation grid, by region
a_el = a_el*grid.df$Area_km2
names(a_el) = strata.limits[1:2,1]

# Add Region to data set
data$Region = NA
for(i in 1:2){
  coords = regions.shp@polygons[[i]]@Polygons[[1]]@coords
  data$Region = ifelse(point.in.polygon(data$lond, data$latd, coords[,1], coords[,2]) %in% c(1,2),i,data$Region)
}

# Adjust the area for the regional indices
Extrapolation_List$a_el = a_el

return(list(data, Extrapolation_List))
}
