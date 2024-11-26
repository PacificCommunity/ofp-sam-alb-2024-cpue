
rm(list=ls())

library(dplyr)
library(magrittr)
library(ggplot2)
library(sdmTMB)
library(sdmTMBextra)
library(INLA)
library(sp)
library(data.table)
library(scales)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## plotting format
###################################################################
gg.theme <- theme_bw() +
  theme(axis.line = element_line(color="black"),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=9),
        text = element_text(size=12),
        axis.title = element_text(size=20),
        strip.background =element_rect(fill='white'),
        strip.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"))  

############################################################################

ltb.cpue.col = c("royalblue3","deepskyblue1","gold","orange1","indianred1",
                             "firebrick2","#AC2020")

load("../Background_Data/SP.ALB.coast.shp.RData")
load("../Background_Data/alb_regions_poly_2024_shp.RData")
load("../Background_Data/coast.RData")

model.name <- paste0("global")

load(file = "../Data/subsample_1954.Rdata")
DG <- sub; rm(sub)
DG %<>% mutate(Response_variable = alb_cpue, month = as.numeric(month))

## convert to equal distant projection
crs_eqd = sp::CRS("+proj=tpeqd +lat_1=-25 +lon_1=170 +lat_2=-25 +lon_2=260 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
crs_ll = sp::CRS(sp::proj4string(SP.ALB.coast.shp))
cpue_coords = DG[,c("lond","latd")]
cpue_sp = sp::SpatialPoints(cpue_coords)
sp::proj4string(cpue_sp) = crs_ll
cpue_eqd = sp::spTransform(cpue_sp, crs_eqd)@coords
colnames(cpue_eqd) = c("lon_eqd","lat_eqd")
DG = cbind(DG,cpue_eqd)
coast_eqd = sp::spTransform(SP.ALB.coast.shp, crs_eqd)
coast_ll = sp::spTransform(SP.ALB.coast.shp, crs_ll)
regions_eqd = sp::spTransform(regions.shp, crs_eqd)
coast_eqd_sf = sf::st_as_sf(coast_eqd) 
coast_sf = sf::st_as_sf(coast_ll) 

## create unique points of sample data for plotting
unique.points <- DG %>% select(c(lond, latd)) %>% distinct()

## remove samples that are on land
pts <- sf::st_as_sf(unique.points, coords=1:2, crs=crs_ll)
## Find which points fall over land
on.land <- !is.na(as.numeric(sf::st_intersects(pts, coast_sf)))

## see what is on land
plot(sf::st_geometry(coast_sf))
plot(pts, col=1+on.land, pch=16, add=TRUE)

## remove land from sampling locations
unique.points <- unique.points[!on.land,]
DG %<>% left_join(unique.points %>% mutate(keep = 1)) %>%
  filter(keep == 1) %>% select(-keep)

## define the mesh
grid_coords = expand.grid(lon=(130.5:287.5),lat=(0.5:-49.5))
grid_sp = sp::SpatialPoints(grid_coords)
sp::proj4string(grid_sp) = crs_ll
## over lap with land
over_index = which(is.na(sp::over(grid_sp,SP.ALB.coast.shp)))
grid_sp = grid_sp[over_index]
plot(grid_sp)

## convert to equi-distant coordinates
grid_sp_eqd = sp::spTransform(grid_sp, crs_eqd)
grid_table = as.data.table(cbind(grid_sp@coords,grid_sp_eqd@coords))
colnames(grid_table) = c("lon","lat","lon_eqd","lat_eqd")

## remove cells that have mean sst < 16C
load(file = "../Background_Data/enviro.RData")
mean.sst <- enviro.dat %>% mutate(lat = lat - 0.5, lon = lon - 0.5) %>%
  filter(lat <= 0) %>% select(lat, lon, sst) %>%
  group_by(lat, lon) %>% summarise(mean.sst = mean(sst, na.rm = T)) %>%
  ungroup() %>% group_by(lat) %>%
  mutate(mean.sst = ifelse(is.na(mean.sst), mean(mean.sst, na.rm = T),
                           mean.sst))
grid_table %<>% left_join(mean.sst) %>% filter(mean.sst >= 16) %>%
  select(-mean.sst)

plot(grid_table$lon, grid_table$lat)

## define mesh
mesh_boundary = INLA::inla.nonconvex.hull(cbind(grid_table$lon_eqd,
                                                grid_table$lat_eqd),
                                          convex = -0.005,resolution=205)
#save(mesh_boundary, file = "../Data/global.mesh.boundary.Rdata")
mesh_inla = INLA::inla.mesh.2d(loc = grid_sp_eqd,boundary = mesh_boundary,
                               max.n.strict=c(500,16),cutoff=770)
mesh_sdmTMB = make_mesh(DG, c("lon_eqd", "lat_eqd"), mesh = mesh_inla)
mesh = add_barrier_mesh(mesh_sdmTMB,coast_eqd_sf,
                        range_fraction = 0.01,plot=T)
#save(mesh, file = "../Data/global.mesh.Rdata")

rm(mesh_sdmTMB)

plot(mesh)
mesh$mesh$n

fmla <- as.formula("Response_variable ~ as.factor(year) + as.factor(flag_id) + s(hbf, bs='bs', m=c(3,2)) + s(month, bs = 'cc', k = 12)")

fmla

graphics.off()

A <- Sys.time()
fit<- sdmTMB(
  fmla,
  time = "year",
  data = DG,
  mesh = mesh,
  spatial = "on",
  spatiotemporal = "iid",
  family = delta_gamma(link1 = "logit", link2 = "log"),
  control = sdmTMBcontrol(newton_loops = 1, parallel = 13L))

B <- Sys.time()
fit.time <- B - A
fit <- c(fit, fit.time = list(fit.time))
class(fit) <- "sdmTMB"

save(fit, file = paste0("../Results/fit.", model.name, ".RData"))
load(file = paste0("../Results/fit.", model.name, ".RData"))

max(fit$gradients)

fit$fit.time
fit$formula

AIC(fit)
logLik(fit)
## extract parameters as a df
fit$formula
print(tidy(fit, conf.int = TRUE))
tidy(fit, effects = "ran_pars", conf.int = TRUE)
## check if realistic results
sanity(fit)

## make predictions
replicate_df <- function(dat, time_name, time_values) {
  nd <- do.call("rbind", replicate(length(time_values), dat, simplify = FALSE))
  nd[[time_name]] <- rep(time_values, each = nrow(dat))
  nd
}

## get the extrapolation grid from the 2021 VAST model
load("../Background_Data/extrapolation.grid.Rdata")

## create new data to predict on
new.data <- extrap.grid

## create unique points of sample data for plotting
unique.points <- new.data %>% select(c(lond, latd)) %>% distinct()

## remove samples that are on land
pts <- sf::st_as_sf(unique.points, coords=1:2, crs=crs_ll)
## Find which points fall over land
on.land <- !is.na(as.numeric(sf::st_intersects(pts, coast_sf)))

## see what is on land
plot(sf::st_geometry(coast_sf))
plot(pts, col=1+on.land, pch=16, add=TRUE)

## remove land from sampling locations
unique.points <- unique.points[!on.land,]
new.data %<>% left_join(unique.points %>% mutate(keep = 1)) %>%
  filter(keep == 1) %>% select(-keep)

## identify regions of interest
num.regions <- 2
new.data$Region = NA
for(i in 1:num.regions){
  coords = regions.shp@polygons[[i]]@Polygons[[1]]@coords
  new.data$Region = ifelse(sp::point.in.polygon(new.data$latd,
                                                new.data$lond, coords[,2],
                                                coords[,1]) %in%
                             c(1,2),i,new.data$Region)
}

## separate sub-regions in WCPO
new.data %<>% mutate(Region = ifelse(Region == 1 & latd <= -25, "1-EF",
                                     ifelse(Region == 1 & latd > -25,
                                            "1-ABCD", Region)))
new.data %>% filter(!is.na(Region)) %>% ggplot() +
  geom_point(aes(x = lond, y = latd, color = as.factor(Region))) +
  geom_sf(data=coast_sf,alpha=0.5)
graphics.off()

fit$formula
## most common flag
new.data[,"flag_id"] <- names(sort(table(DG$flag_id),decreasing=TRUE)[1])
new.data[,"hbf"] <- median(DG$hbf, na.rm = T)
new.data <- data.frame(new.data)

## replicate every location for every time step
new.data %<>% filter(!is.na(Region))
new.data <- replicate_df(new.data, "year", unique(DG$year))
new.data %<>% mutate(month = 7) # treat month as catchability covariate

## predict by region-area
num.regions <- 3
regs <- unique(new.data$Region)
A <- Sys.time()
for(i in 1:num.regions){
  new.df <- new.data %>% filter(Region == regs[i])
  p.tmp <- predict(fit, newdata = new.df, return_tmb_object = TRUE)
  if(!exists("pred.df")){
    pred.df <- p.tmp$data
  } else {
    pred.df %<>% rbind(p.tmp$data)
  }
  save(pred.df, file = paste0("../Results/pred.", model.name, ".RData"))
  rm(p.tmp, new.df)
}
B <- Sys.time()
pred.time <- B - A
pred.time

load(file = paste0("../Results/pred.", model.name, ".RData"))

## make sure the link functions are appropriate to your model
index.final <- pred.df %>%
  mutate(prob1 = boot::inv.logit(est1), prob2 = exp(est2)) %>%
  mutate(rel_abund = prob1 * prob2 * Area_km2) %>%
  group_by(year, Region) %>% summarise(rel_abund = sum(rel_abund)) %>%
  ungroup()

save(index.final, file = paste0("../Results/index.", model.name, ".RData"))
load(file = paste0("../Results/index.", model.name, ".RData"))

index.final %<>% mutate(Model = "global") %>%
  select(year, Index = rel_abund, Region, Model) %>%
  mutate(YrQtr = year, Region = paste("Region", Region))

# Make plots
dir.create(paste0("../Figures/", model.name),
           recursive = T)
setwd(paste0("../Figures/", model.name))

## mean-centered
index.final %>%
  group_by(Region, Model) %>%
  mutate(std = Index / mean(Index)) %>%
ggplot(aes(YrQtr, std, group = as.factor(Model))) +
  geom_line(linewidth = 0.75) +
  geom_point(size=1) +
  labs(y = 'Mean-Centered Standardized Indices', x ='Year') +
  theme_bw() + gg.theme + theme(legend.position = "none") + 
  facet_grid(rows = vars(Region), scales='free')

ggsave(paste0('MC.Std.Indices.', model.name, '.png'), width=10, height=7,
       units='in', dpi=300)

## get index in 5x5 cells for Tom Peatman
ji <- function(xy, origin=c(0,0), cellsize=c(5,5)) {
  t(apply(xy, 1,
          function(z) cellsize/2+origin+cellsize*(floor((z - origin)/cellsize))))
}
index.1x1 <- pred.df %>%
  mutate(prob1 = boot::inv.logit(est1), prob2 = exp(est2)) %>%
  mutate(rel_abund = prob1 * prob2 * Area_km2)
JI <- ji(cbind(index.1x1$lond, index.1x1$latd), cellsize = c(5, 5))
index.1x1$lon_5 <- JI[, 1]
index.1x1$lat_5 <- JI[, 2]

index.5x5 <- index.1x1 %>%
  group_by(lat_5, lon_5, year) %>%
  summarise(rel_abund = sum(rel_abund)) %>% ungroup()
rm(JI)

## check against the previous estimates by plotting
## assign regions
index.5x5$Region = NA
for(i in 1:2){
  coords = regions.shp@polygons[[i]]@Polygons[[1]]@coords
  index.5x5$Region = ifelse(sp::point.in.polygon(index.5x5$lon_5,
                                                  index.5x5$lat_5, coords[,1],
                                                  coords[,2]) %in%
                               c(1,2),i,index.5x5$Region)
}

## separate sub-regions in WCPO
index.5x5%<>% mutate(Region = ifelse(Region == 1 & lat_5 <= -25, "1-EF",
                                     ifelse(Region == 1 & lat_5 > -25,
                                            "1-ABCD", Region)))
index.final %<>%
  rbind(index.5x5 %>% filter(!is.na(Region)) %>%
          mutate(Region = paste("Region", Region)) %>%
          group_by(year, Region) %>%
          summarise(Index = sum(rel_abund)) %>% ungroup() %>%
          mutate(Model = "5x5", YrQtr = year) %>%
          select(year, YrQtr, Region, Index, Model))

index.final %>% filter(Model %in% c("5x5", "global")) %>%
  group_by(Region, Model) %>%
  mutate(std = Index / mean(Index)) %>%
  ggplot(aes(YrQtr, std, group = as.factor(Model))) +
  geom_line(aes(col=as.factor(Model)), linewidth = 0.75) +
  geom_point(size=1, aes(col=as.factor(Model))) +
  labs(y = 'Mean-Centered Std. Indices', x ='YrQtr') +
  theme_bw() + gg.theme + labs(col = "Model") + 
  facet_grid(rows = vars(Region), scales='free')

index.5x5 %<>% mutate(Qtr = 3)

yrs <- sort(unique(index.5x5$year))
index.5x5 %<>% left_join(data.frame(year = yrs, ts = seq(1, length(yrs), 1)))

write.csv(index.5x5, file = paste0("../../Results/", model.name,
                                   ".index.5x5.csv"), row.names = F)
index.5x5 <- read.csv(file = paste0("../../Results/", model.name,
                                   ".index.5x5.csv"))

## compare with nominal
## add area to data
library(raster)
r <- raster(ymn = -51, ymx = 1, res = 1)
a <- area(r)
lat <- yFromRow(r, 1:nrow(r))
area <- a[,1]
plot(lat, area, type = "l"); points(lat, area, pch = 16)
area.df <- data.frame(latd = lat,
                      Area_km2 = area)

DG %<>% left_join(area.df)
detach('package:raster')

nom.idx <- DG %>% filter(!is.na(Region)) %>%
  mutate(Region = ifelse(Region == 1 & latd <= -25, "1-EF",
                         ifelse(Region == 1 & latd > -25,
                                "1-ABCD", Region))) %>%
  group_by(year, Region, latd, lond) %>%
  summarise(CPUE = mean(Response_variable),
            area.weight = sum(Area_km2)) %>% ungroup() %>%
  mutate(CPUE.weighted = CPUE * area.weight) %>% group_by(year, Region) %>%
  summarise(Index = sum(CPUE.weighted)) %>% ungroup()

index.comb <- index.final %>% filter(Model == "global") %>%
  mutate(Model = "Standardized") %>%
  select(Region, Index, year, Model) %>%
  rbind(nom.idx %>% mutate(Model = "Nominal",
                           Region = paste("Region", Region)))

## nominal
index.comb %>%
  group_by(Region, Model) %>%
  mutate(std = Index / mean(Index)) %>%
  ggplot(aes(year, std, group = as.factor(Model))) +
  geom_line(aes(col=as.factor(Model)), linewidth = 0.5) +
  scale_x_continuous(breaks = seq(1950, 2030, 5)) +
  scale_color_manual(values = c("green", "black", "grey")) +
  labs(y = 'Mean-Centered Indices', x ='Year') +
  theme_bw() + gg.theme + labs(col = "Model") + 
  facet_grid(rows = vars(Region), scales='free')

ggsave(paste0('MC.Std.Nominal.Indices.', model.name, '.png'), width=10,
       height=7, units='in', dpi=300)

## plot sub-sampled data spatially by decade
DG %>% mutate(decade = floor(YrQtr / 10) * 10) %>%
  group_by(latd, lond, decade) %>%
  summarise(samp.size = n()) %>% ungroup() %>%
  ggplot() +
  facet_wrap(~decade) +
  geom_tile(aes(x = lond, y = latd, fill = samp.size)) +
  geom_polygon(data = SP.ALB.coast.shp,
               aes(x = long, y = lat, group = group),
               colour = "black", fill = "bisque3") +
  geom_polygon(data = regions.shp,
               aes(x = long, y = lat, group = group),
               colour = "firebrick4", fill = NA) +
  scale_fill_gradientn("# Observations", colours = ltb.cpue.col) +
  theme_void() + gg.theme +
  labs(x = "Longitude", y = "Latitude")

ggsave(paste0('subsample.distribution.', model.name, '.png'), width=10,
       height=3, units='in', dpi=300)

## plot CPUE spatially by 5-year intervals
spat.cpue.df <- pred.df %>%
  mutate(prob1 = boot::inv.logit(est1), prob2 = exp(est2)) %>%
  mutate(rel_abund = prob1 * prob2 * Area_km2) %>%
  mutate(year = 5 * floor(year / 5)) %>% group_by(year, latd, lond) %>%
  summarise(rel_abund = sum(rel_abund)) %>% ungroup()

spat.cpue.df %>%
  ggplot() +
  facet_wrap(~year) +
  geom_tile(aes(x = lond, y = latd, fill = rel_abund)) +
  geom_polygon(data = SP.ALB.coast.shp,
               aes(x = long, y = lat, group = group),
               colour = "black", fill = "bisque3") +
  geom_polygon(data = regions.shp,
               aes(x = long, y = lat, group = group),
               colour = "firebrick4", fill = NA) +
  scale_fill_gradient2("Relative Abundance", high = muted("blue"),
                       low = muted("white")) +
  theme_void() + gg.theme +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_blank()) +
  labs(x = "Longitude", y = "Latitude")

ggsave(paste0('prop.cpue.by.year.', model.name, '.png'), width=20,
       height=9, units='in', dpi=300)

## Model diagnostics & plot output
## need to simulate responses to calculate DHARMa style residuals
## can use built in functions from sdmTMB package but DHARMa package has
## some built in tests/plots
set.seed(1134)
load(file = paste0("../../Results/fit.", model.name, ".RData"))

sim_fit = simulate(fit, nsim = 500)
pred.DG <- predict(fit)

pred.DG %<>% left_join(area.df) %>%
  mutate(prob1 = boot::inv.logit(est1), prob2 = exp(est2)) %>%
  mutate(rel_abund = prob1 * prob2 * Area_km2)

save(pred.DG, file = paste0("../../Results/pred.DG.", model.name,
                            ".RData"))
load(file = paste0("../../Results/pred.DG.", model.name,
                   ".RData"))

residuals_dharma = DHARMa::createDHARMa(simulatedResponse = sim_fit,
                                        observedResponse = pred.DG$Response_variable)	

save(residuals_dharma, file = paste0("../../Results/resid.dharma.",
                                     model.name, ".RData"))
load(file = paste0("../../Results/resid.dharma.",
                   model.name, ".RData"))

## basic QQ & residual v. predicted plot
plot(residuals_dharma)
# residual v. predicted by model covariate
DHARMa::plotResiduals(residuals_dharma, form = as.factor(pred.DG$year))
DHARMa::plotResiduals(residuals_dharma, form = pred.DG$lond)
DHARMa::plotResiduals(residuals_dharma, form = pred.DG$latd)
DHARMa::plotResiduals(residuals_dharma, form = pred.DG$hbf)
DHARMa::plotResiduals(residuals_dharma, form = pred.DG$month)
DHARMa::plotResiduals(residuals_dharma, form = pred.DG$flag_id)
# test for uniformity in residuals, overdispersion, outliers
DHARMa::testResiduals(residuals_dharma)
# test for zero inflation
DHARMa::testZeroInflation(residuals_dharma)
# test for spatial autocorrelation
# DHARMa::testSpatialAutocorrelation needs unique locations
pred.DG$spatial_group = as.numeric(as.factor(paste0(pred.DG$lond,"_",
                                                    pred.DG$latd)))
spatial_group_dt = data.table(spatial_group=pred.DG$spatial_group,
                              x=pred.DG$lond,y=pred.DG$latd) %>%
  .[,.(x=mean(x), y=mean(y)), by=spatial_group]
residuals_spatial_group = DHARMa::recalculateResiduals(residuals_dharma,
                                                       group = pred.DG$spatial_group)	
DHARMa::testSpatialAutocorrelation(residuals_spatial_group,
                                   x = spatial_group_dt$x,
                                   y = spatial_group_dt$y)
# test for temporal autocorrelation
pred.DG$temporal_group = factor(pred.DG$year,levels=as.character(sort(unique(pred.DG$year))))
#pred.DG$temporal_group = factor(pred.DG$YrQtr,levels=as.character(sort(unique(pred.DG$YrQtr))))
residuals_temporal_group = DHARMa::recalculateResiduals(residuals_dharma, group = pred.DG$temporal_group)	
DHARMa::testTemporalAutocorrelation(residuals_temporal_group, time=levels(pred.DG$temporal_group))

pred.DG$temporal_group =factor(pred.DG$year,levels=as.character(sort(unique(pred.DG$year))))
residuals_temporal_group = DHARMa::recalculateResiduals(residuals_dharma, group = pred.DG$temporal_group)	
DHARMa::testTemporalAutocorrelation(residuals_temporal_group, time=levels(pred.DG$temporal_group))

## spatial plots of residuals
## assign data points to 5x5 cells
JI <- ji(cbind(pred.DG$lond, pred.DG$latd))
pred.DG$lon_5 <- JI[, 1]
pred.DG$lat_5 <- JI[, 2]
pred.DG %<>% mutate(cell.5x5 = paste0(lat_5, "-", lon_5))
rm(JI)

# Histogram of residuals aggregated to annual time scale
pred.DG %<>% mutate(scaled_res = residuals_dharma$scaledResiduals)
res.agg = pred.DG %>%
  group_by(year, cell.5x5) %>%
  summarise(res = mean(scaled_res, na.rm=T))
res.agg %<>% left_join(pred.DG %>% select(cell.5x5, lat_5, lon_5) %>%
                         distinct())

newmap = fortify(maps::map(wrap=c(0,360), plot=FALSE, fill=TRUE))

ggplot() +   
  geom_tile(data=res.agg, 
            aes(lon_5, lat_5, fill=res)) +
  xlab("Longitude") + ylab("Latitude") + 
  scale_x_continuous(limits=c(min(res.agg$lon_5), max(res.agg$lon_5)),
                     breaks = seq(min(res.agg$lon_5), max(res.agg$lon_5), 10)) +
  coord_equal()+
  scale_y_continuous(limits=c(-50, 5)) +
  labs(fill='Residual')+ 
  scale_fill_gradient2(low='firebrick', high='navy', midpoint=0.5)   +
  geom_map(data=newmap, map=newmap,
           aes(map_id=region), fill='gray90', col='gray60') +  gg.theme

ggsave(paste0('Dharma_res_', model.name, '.png'), width=12, height=8,
       units='in', dpi=300)

ggplot(res.agg, aes(res)) + geom_histogram(col='white') +
  gg.theme + labs(y = 'Frequency', x = 'PIT residual')

ggsave(paste0('Dharma_res_hist_', model.name, '.png'), width=12, height=8,
       units='in', dpi=300)

##_____________________________________________________________________________________________________________________________
land_col = "gray90"
## estimate uncertainty
## conduct boot-strapping from joint precision matrix
predict_sim = predict(fit, newdata = new.data, nsim=500)

save(predict_sim, file = paste0("../../Results/predict_sim.",
                                     model.name, ".RData"))
load(file = paste0("../../Results/predict_sim.",
                                model.name, ".RData"))

predict_nd_sim = cbind(new.data,predict_sim) %>% as.data.table() %>%
  melt(.,id.vars=colnames(new.data))

save(predict_nd_sim, file = paste0("../../Results/predict_nd_sim.",
                                model.name, ".RData"))
load(file = paste0("../../Results/predict_nd_sim.",
                   model.name, ".RData"))

index.cv.5x5 <- predict_nd_sim %>% .[,.(mean=mean(exp(value)), sd=sd(exp(value))),
                                 by=.(year,lon_eqd,lat_eqd, latd, lond)] %>%
  mutate(cv = sd / mean) %>%
  select(year, latd, lond, cv) %>%
  left_join(index.5x5 %>%
              select(year, latd = lat_5, lond = lon_5, Index = rel_abund))

## assign regions
index.cv.5x5$Region = NA
for(i in 1:2){
  coords = regions.shp@polygons[[i]]@Polygons[[1]]@coords
  index.cv.5x5$Region = ifelse(sp::point.in.polygon(index.cv.5x5$lond,
                                                 index.cv.5x5$latd, coords[,1],
                                                 coords[,2]) %in%
                              c(1,2),i,index.cv.5x5$Region)
}
## separate sub-regions in WCPO
index.cv.5x5 %<>% mutate(Region = ifelse(Region == 1 & latd <= -25,
                                                    "1-EF",
                                     ifelse(Region == 1 & latd > -25,
                                            "1-ABCD", Region)))

index.cv.5x5 %>% ggplot() +
  geom_point(aes(x = lond, y = latd, color = as.factor(Region))) +
  geom_sf(data=coast_sf,fill=land_col,alpha=0.5)
  
graphics.off()

index.cv <- index.cv.5x5 %>% filter(!is.na(Region)) %>%
  group_by(Region, year) %>%
  summarise(Index = sum(Index), cv = mean(cv)) %>%
  ungroup()

write.csv(index.cv, file = paste0("../../Results/index.cv.",
                                  model.name, ".csv"), row.names = F)
index.cv <- read.csv(file = paste0("../../Results/index.cv.",
                                   model.name, ".csv"))

## plot estimated response uncertainty
predict_nd_sim %>% .[,.(mean=mean(exp(value)),sd=sd(exp(value))),
                     by=.(year,lon_eqd,lat_eqd, latd, lond)] %>%
  .[,cv:=sd/mean] %>%
  data.frame() %>% mutate(decade = floor(year / 10) * 10) %>%
  group_by(lat_eqd, lon_eqd, decade, latd, lond) %>%
  summarise(cv = mean(cv)) %>%
  ggplot() + 
  facet_wrap(~decade) +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_point(aes(x = lon_eqd, y = lat_eqd, color = cv),size=3) +
  geom_sf(data=coast_eqd_sf,fill=land_col,alpha=0.5) +
  gg.theme + theme(axis.text.x=element_blank(),
                   axis.text.y=element_blank()) +
  viridis::scale_color_viridis("CV (response)",begin = 0.1,end = 0.8,
                               direction = 1,option = "H",trans="sqrt",
                               limits = c(0, 8))

ggsave(paste0('Uncertainty.spatial.', model.name, '.png'), width=12,
       height=8, units='in', dpi=300)

## get region specific indices with uncertainty
A <- Sys.time()
for(i in 1:num.regions){
  new.df <- new.data %>% filter(Region == regs[i])
  p.tmp <- predict(fit, newdata = new.df, return_tmb_object = TRUE)
  area.vec <- new.df$Area_km2 %>% unlist()
  index.tmp <- get_index(p.tmp, area = area.vec,
                         bias_correct = T)
  index.tmp <- index.tmp %>% mutate(Region = regs[i])
  
  if(!exists("index.df")){
    index.df<- index.tmp
  } else {
    index.df <- rbind(index.df, index.tmp)
  }
  save(index.df, file = paste0("../../Results/index.conf.intervs.",
                               model.name, ".RData"))
  rm(index.tmp, p.tmp, area.vec, new.df)
}
B <- Sys.time()
idx.time <- B - A
idx.time

## combine the various model results
load(file = paste0("../../Results/index.conf.intervs.",
                             model.name, ".RData"))
NZ.cpue <- read.csv(file = "../../Results/alb1_daily_full_indices.csv")
NZ.cpue %<>% mutate(year = fyear, est = index_Positive, Region = "1-EF") %>%
  mutate(lwr = est - 2 * sd, upr = est + 2 * sd, log_est = log(est)) %>%
  select(year, est, lwr, upr, log_est, se = cv, Region) %>%
  mutate(model = "NZ Troll")
idx <- index.df %>% mutate(model = "Global") %>% rbind(NZ.cpue)
load(file = paste0("../../Results/index.conf.intervs.spawn.RData"))

idx %<>% rbind(index.df %>% mutate(model = "Spawn", Region = "1-ABCD"))

## mean-centered with uncertainty
idx %>% filter(model == "Global" | model == "Spawn") %>%
  group_by(model, Region) %>%
  mutate(est = est / mean(est)) %>% ungroup() %>%
  mutate(sd = se * est) %>%
  mutate(lwr = est - 2 * sd, upr = est + 2 * sd) %>%
  mutate(lwr = ifelse(lwr < 0, 0, lwr)) %>%
  ggplot() +
  geom_line(aes(year, est, col = model), linewidth = 0.85) +
  geom_ribbon(aes(x = year, ymin=lwr,ymax=upr, fill = model), alpha=0.3) +
  scale_color_manual(values = ltb.cpue.col[c(1,9)]) +
  scale_fill_manual(values = ltb.cpue.col[c(1,9)]) +
  scale_x_continuous(breaks = seq(1950, 2030, 10)) +
  labs(y = 'Mean-centered Std. Indices (Confedence Interval)', x ='Year',
       col = "Index Est.", fill = "Index CI") +
  theme_bw() + gg.theme +
  facet_grid(rows = vars(Region), scales='free')

ggsave(paste0('Std.Indices.with.uncertainty.', model.name, '.png'), width=10,
       height=7, units='in', dpi=300)

## Estimated fixed effects
fit$formula

pred.DG %>% 
  ggplot() +
  geom_boxplot(aes(x = as.factor(hbf), y = rel_abund)) +
  gg.theme + labs(x = "Hooks between Floats",
                  y = "Predicted CPUE")

ggsave(paste0('hbf.', model.name, '.png'), width=10,
       height=7, units='in', dpi=300)

pred.DG %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(flag_id), y = rel_abund)) +
  gg.theme + labs(x = "FLAG", y = "Predicted CPUE")

ggsave(paste0('FLAG.', model.name, '.png'), width=10, height=7,
       units='in', dpi=300)

pred.DG %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(month), y = rel_abund)) +
  gg.theme + labs(x = "Month", y = "Predicted CPUE")

ggsave(paste0('Month.', model.name, '.png'), width=10, height=7,
       units='in', dpi=300)

## encounter prob
pred.DG %>% filter(!is.na(Region)) %>%
  mutate(Region = ifelse(Region == 1, "WCPO",
                         ifelse(Region == 2, "EPO", NA))) %>%
  mutate(Region = factor(Region, levels = c("WCPO", "EPO"))) %>%
  group_by(year, Region) %>%
  summarise(prob1 = mean(prob1)) %>% filter(!is.na(Region)) %>%
  ggplot() +
  facet_grid(rows = vars(Region)) +
  geom_point(aes(x = year, y = prob1)) +
  geom_smooth(aes(x = year, y = prob1), size = 0.2) +
  scale_x_continuous(breaks = seq(1950, 2030, 10)) +
  gg.theme + labs(x = "Year", y = "Encounter Probability")

ggsave(paste0('encounter.prob.', model.name, '.png'), width=8, height=8,
       units='in', dpi=300)

# positive component
pred.DG %>% filter(!is.na(Region)) %>%
  mutate(Region = ifelse(Region == 1, "WCPO",
                         ifelse(Region == 2, "EPO", NA))) %>%
  mutate(Region = factor(Region, levels = c("WCPO", "EPO"))) %>%
  group_by(year, Region) %>%
  summarise(prob2 = mean(prob2)) %>% filter(!is.na(Region)) %>%
  ggplot() +
  facet_grid(rows = vars(Region)) +
  geom_point(aes(x = year, y = prob2)) +
  geom_smooth(aes(x = year, y = prob2), size = 0.2) +
  scale_x_continuous(breaks = seq(1950, 2030, 10)) +
  gg.theme + labs(x = "Year", y = "Positive Component")

ggsave(paste0('positive.comp.', model.name, '.png'), width=8, height=8,
       units='in', dpi=300)
