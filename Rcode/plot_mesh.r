
## Plot spatial mesh for albacore to match the 2021 knot locations

##_____________________________________________________________________________________________________________________________
## load packages
library(data.table)
library(dplyr)
library(magrittr)
library(sdmTMB)
library(ggplot2)
library(sdmTMBextra)

##_____________________________________________________________________________________________________________________________
## set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list=ls())

##_____________________________________________________________________________________________________________________________
## bring in region and coast shapes
load(file="../Background_Data/alb_regions_poly_2024_shp.RData")
load(file = "../Background_Data/SP.ALB.coast.shp.RData")
##_____________________________________________________________________________________________________________________________

## load the global mesh
load(file = "../Data/global.mesh.boundary.Rdata")
load(file = "../Data/global.mesh.Rdata")

## define equidistant coordinates
crs_eqd = sp::CRS("+proj=tpeqd +lat_1=-25 +lon_1=170 +lat_2=-25 +lon_2=260 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
coast_eqd = sp::spTransform(SP.ALB.coast.shp, crs_eqd)
regions_eqd = sp::spTransform(regions.shp, crs_eqd)

# plotting
# add plot to show the mesh structure
# pull-out needed quantities for plotting
t.sub=1:nrow(mesh$mesh$graph$tv)
idx = cbind(mesh$mesh$graph$tv[t.sub,c(1:3,1), drop=FALSE], NA)
x = mesh$mesh$loc[t(idx), 1]
y = mesh$mesh$loc[t(idx), 2]
idx_normal = idx[-mesh$barrier_triangles,]
x_normal = mesh$mesh$loc[t(idx_normal), 1]
y_normal = mesh$mesh$loc[t(idx_normal), 2]
idx_boundary = mesh_boundary$idx
# add NA at each non-consecutive index
non_consec_idx = c(1,which(diff(idx_boundary[,2])<0)+1)
new_idx_boundary = matrix(c(NA,1),nrow=1,ncol=2)
for(i in 1:(length(non_consec_idx)-1)){
  new_idx_boundary = rbind(new_idx_boundary,
                           idx_boundary[non_consec_idx[i]:non_consec_idx[i+1],],
                           c(idx_boundary[non_consec_idx[i+1],1],NA))
}
# remove duplicate indices for plotting
new_idx_boundary = unique(new_idx_boundary)
x_boundary = mesh_boundary$loc[t(new_idx_boundary), 1]
y_boundary = mesh_boundary$loc[t(new_idx_boundary), 2] 

normal_col = "black"
  barrier_col = "gray70"
    boundary_col = "hotpink"
      land_col = "gray90"
        region_col = "blue"
          
        png(filename = paste0("../Figures/global/mesh_eqd.png"),width = 12,
            height = 8, units = "in", res=300)
        par(mar=c(5,5,1,1))
        plot(x,y,type="n",axes=TRUE,xlab="Eastings (km)",ylab="Northings (km)",
             cex=1.5,cex.axis=1.5,cex.lab=1.5,las=1)
        sp::plot(coast_eqd,add=TRUE,col=land_col)
        lines(x,y,col=barrier_col)
        points(x,y,pch=16,cex=0.5,col=barrier_col)
        lines(x_normal,y_normal,col=normal_col)
        lines(x_boundary,y_boundary,col=boundary_col,lwd=3.5)
        points(x_normal,y_normal,pch=16,cex=0.5,col=normal_col)
        sp::plot(regions_eqd,border=region_col,lwd=2,add=TRUE)
        legend("bottomleft",legend=c("normal knot","barrier knot",
                                     "normal edate","barrier edate",
                                     "knot boundary","region boundary","land"),
               lwd=c(NA,NA,1,1,3,3,1),pch=c(16,16,NA,NA,NA,NA,NA),
               col=c(normal_col,barrier_col,normal_col,barrier_col,
                     boundary_col,region_col,NA),
               fill=c(NA,NA,NA,NA,NA,NA,land_col),
               border=c(NA,NA,NA,NA,NA,NA,"black"),bty="n",cex=0.75)
        dev.off()

        