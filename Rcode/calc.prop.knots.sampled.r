rm(list=ls())

library(dplyr)
library(data.table)
library(magrittr)
library(ggplot2)
library(ggthemes)
library(sp)
library(sf)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Formatting
###################################################################
newmap = fortify(maps::map(wrap=c(0,360), plot=FALSE, fill=TRUE))


colpal <- function(n=12) {
  fb35 <- "#AC2020" 
    colorRampPalette(c("royalblue3","deepskyblue1","gold", "orange1",
                                   "indianred1","firebrick2",fb35))(n)
}

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                        "#0072B2", "#D55E00", "#CC79A7")
                        
gg.theme <- theme_bw() + theme(axis.line = element_line(color="black"),
                               #panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               text = element_text(size=14),
                               axis.title = element_text(size=16),
                               strip.background =element_rect(fill="white"),
                               strip.text = element_text(size=14),
                               legend.box.background = element_blank(),
                               legend.background = element_blank(),
                               panel.background = element_rect(fill = NA,
                                                               color = "black"))

axis_labels <- rbind(data.frame(long = rep(120,6), lat = seq(-50,0,10),
                                labels = seq(-50,0,10)),
                     data.frame(long = seq(140,280,20), lat = rep(-55,8),
                                labels = seq(140,280,20)),
                     data.frame(long = c(110, 200), lat = c(-25, -60),
                                labels = c("Lat", "Long")))

############################################################################

## Define regions
load("../Background_Data/SP.ALB.coast.shp.RData")
load("../Background_Data/alb_regions_poly_2024_shp.RData")
coast.sf <- st_as_sf(SP.ALB.coast.shp)
regions.sf <- st_as_sf(regions.shp)

## load dataset and sdmTMB global run data
load(file = "../Data/subsample_1954.Rdata")
DG <- sub; rm(sub)
load(file = "../Results/fit.global.RData")

## knot locations
knot.locs.eqd <- data.frame(fit$spde$mesh$loc) %>% select(-X3)
colnames(knot.locs.eqd) <- c("eqd_lon", "eqd_lat")

## change to lat lons
crs_ll = sp::CRS(sp::proj4string(SP.ALB.coast.shp))
crs_eqd = sp::CRS("+proj=tpeqd +lat_1=-25 +lon_1=170 +lat_2=-25 +lon_2=260 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
knot.locs.eqd.sp = sp::SpatialPoints(knot.locs.eqd)
sp::proj4string(knot.locs.eqd.sp) = crs_eqd
knot.locs.sp = sp::spTransform(knot.locs.eqd.sp, crs_ll)@coords
colnames(knot.locs.sp) = c("Lon","Lat")
knot.locs <-  cbind(knot.locs.eqd, knot.locs.sp)
knot.locs %<>% mutate(Lon =  ifelse(Lon < 0, 360 + Lon, Lon))

## convert original data points to sf
orig.dat.sf <- DG %>% select(latd, lond)
coordinates(orig.dat.sf) <- ~lond+latd
orig.dat.sf@proj4string <- crs_ll
orig.dat.sf %<>% st_as_sf(orig.dat.sf)

## identify regions of interest
num.regions <- 2
knot.locs$Region = NA
for(i in 1:num.regions){
  coords = regions.shp@polygons[[i]]@Polygons[[1]]@coords
  knot.locs$Region = ifelse(sp::point.in.polygon(knot.locs$Lat,
                                                 knot.locs$Lon, coords[,2],
                                                 coords[,1]) %in%
                              c(1,2),i,knot.locs$Region)
}

knot.locs.sf <- knot.locs %>% filter(!is.na(Region))
coordinates(knot.locs.sf) <- ~Lon+Lat
knot.locs.sf@proj4string <- crs_ll
knot.locs.sf %<>% st_as_sf(knot.locs.sf)

## check plot
ggplot() +
  geom_sf(data = coast.sf) +
  geom_text(data = axis_labels, aes(x = long, y = lat, label = labels),
            size = 5) +
  geom_segment(aes(x = 130 , y = seq(-50, 0, 10), xend = 287.5,
                   yend = seq(-50, 0, 10)), col = "lightgrey") +
  geom_segment(aes(x = seq(140, 280, 20), y = -50, xend = seq(140, 280, 20),
                   yend = 0), col = "lightgrey") +
  geom_sf(data = regions.sf, fill = NA, col = "blue") +
  geom_sf(data = orig.dat.sf, col = "cyan") +
  geom_sf(data = knot.locs.sf %>% filter(!is.na(Region)), size = 2, col = "darkred") +
  theme(axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = NA, color = "white"))

knot.dt <- as.data.frame(summary(fit$spde$A_st)) %>%
  group_by(i) %>% summarise(x = max(x)) %>% ungroup() %>%
  left_join(as.data.frame(summary(fit$spde$A_st))) %>% select(-x) %>%
  left_join(data.frame(fit$data) %>% select(yyqq = YrQtr) %>%
              mutate(i = seq(1, nrow(fit$data)))) %>% select(-i) %>%
  left_join(knot.locs %>% mutate(j = seq(1, nrow(knot.locs), 1))) %>%
  rename(knot = j) %>% filter(!is.na(Region))

knots_sampled = knot.dt %>% group_by(yyqq, Region) %>%
  summarise(n_sampled = n_distinct(knot), N = n()) %>% ungroup()

knot.locs %<>% mutate(knot = seq(1, nrow(fit$spde$mesh$loc), 1))

knots_per_region <- knot.locs %>% filter(!is.na(Region)) %>%
  group_by(Region) %>%
  summarise(total_knots = n_distinct(knot))

knots_sampled %<>% left_join(knots_per_region) %>%
  mutate(prop_sampled = n_sampled / total_knots) %>%
  mutate(Region = as.factor(Region))

## fill in missing rows
knots_sampled %<>%
  right_join(expand.grid(yyqq = seq(min(knots_sampled$yyqq),
                                    max(knots_sampled$yyqq), 0.25),
                         Region = as.factor(seq(1, num.regions, 1))) %>%
               left_join(knots_sampled %>% select(Region, total_knots) %>%
                           distinct()))
knots_sampled[is.na(knots_sampled)] <- 0
knots_sampled %<>% mutate(Qtr = as.factor((yyqq - floor(yyqq) + 0.25) * 4))

# plot results
p = knots_sampled %>%
  ggplot() +
  xlab("Year") + ylab("Proportion of Spatial Knots Sampled") +
  scale_x_continuous(breaks=c(1960,1990,2020)) +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  facet_wrap(~Region) +
  geom_hline(yintercept=c(0),size=1) +
  geom_hline(yintercept=c(0,0.2,0.4),col="gray70") +
  geom_point(aes(x=yyqq,y=prop_sampled,fill=Region,size=N),
             shape=21,color="black",alpha=0.4) +
  geom_smooth(aes(x=yyqq,y=prop_sampled),linewidth=2,se=FALSE,color="black") +
  geom_smooth(aes(x=yyqq,y=prop_sampled,color=Region),linewidth=1.25) +
  theme_few(base_size=20)
p

ggsave(filename=paste0("../Figures/Catch_effort_summaries/prop_knots_sampled.png"), plot = p, device = "png",
       width = 16, height = 9, units = c("in"),
       dpi = 300, limitsize = TRUE)

## check spatiotemporal coverage
yrqtr.cells <- DG %>% mutate(yrqtr.cell.5x5 = paste0(YrQtr, "-", cell.5x5)) %>%
  select(yrqtr.cell.5x5) %>% mutate(freq = 1) %>% group_by(yrqtr.cell.5x5) %>%
  summarise(freq = sum(freq))
summary(yrqtr.cells$freq)
table(yrqtr.cells$freq)
