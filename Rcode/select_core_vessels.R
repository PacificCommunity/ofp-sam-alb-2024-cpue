rm(list=ls())

library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(sf)
library(sp)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read in filtered data
load('../Data/ops.model.Rdata')
load(file="../Background_Data/SP.ALB.coast.shp.RData")
load(file="../Background_Data/alb_regions_poly_2021_shp.RData")

min.threshold <- 3 ## minimum number of samples per year-quarter per cell.5x5

crs_ll <- SP.ALB.coast.shp@proj4string
coast.sf <- st_as_sf(SP.ALB.coast.shp)
regions.sf <- st_as_sf(regions.shp)

axis_labels <- rbind(data.frame(long = rep(120,6), lat = seq(-50,0,10),
                                labels = seq(-50,0,10)),
                     data.frame(long = seq(140,280,20), lat = rep(-55,8),
                                labels = seq(140,280,20)),
                     data.frame(long = c(110, 200), lat = c(-25, -60),
                                labels = c("Lat", "Long")))


df1 <- ops.filt %>% filter(vessel_id == "5982") %>%
  select(latd, lond) %>% data.frame()

centroid.df1 <- df1 %>% summarise(latd = mean(latd), lond = mean(lond))

coordinates(df1) <- ~lond+latd
df1@proj4string <- crs_ll
df1 %<>% st_as_sf(df1)

coordinates(centroid.df1) <- ~lond+latd
centroid.df1@proj4string <- crs_ll
centroid.df1 %<>% st_as_sf(centroid.df1)

#unique(ops.filt$vessel_id)
df2 <- ops.filt %>% filter(vessel_id == "2000") %>%
  select(latd, lond) %>% data.frame()

centroid.df2 <- df2 %>% summarise(latd = mean(latd), lond = mean(lond))

coordinates(df2) <- ~lond+latd
df2@proj4string <- crs_ll
df2 %<>% st_as_sf(df2)

coordinates(centroid.df2) <- ~lond+latd
centroid.df2@proj4string <- crs_ll
centroid.df2 %<>% st_as_sf(centroid.df2)

dir.create("../Figures/vessel/")

ggplot() +
  geom_sf(data = coast.sf) +
  geom_text(data = axis_labels, aes(x = long, y = lat, label = labels),
            size = 5) +
  geom_segment(aes(x = 130 , y = seq(-50, 0, 10), xend = 287.5,
                   yend = seq(-50, 0, 10)), col = "lightgrey") +
  geom_segment(aes(x = seq(140, 280, 20), y = -50, xend = seq(140, 280, 20),
                   yend = 0), col = "lightgrey") +
  geom_sf(data = regions.sf, fill = NA, col = "blue") +
  geom_sf(data = df1, size = 1, col = "deeppink") +
  geom_sf(data = centroid.df1, size = 4, col = "darkred") +
  geom_sf(data = df2, size = 1, col = "cyan") +
  geom_sf(data = centroid.df2, size = 4, col = "darkblue") +
  theme(axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = NA, color = "white"))

ggsave('../Figures/vessel/example.of.spatial.variability.png', width=10,
       height=5, units='in', dpi=300)

rm(df1, df2, centroid.df1, centroid.df2)

## evaluate longevity of vessels
temporal.df <- ops.filt %>% filter(!is.na(vessel_id)) %>%
  select(vessel_id, YrQtr) %>%
  group_by(vessel_id) %>%
  mutate(yrs.fished = max(YrQtr) - min(YrQtr)) %>%
  summarise(mean.yrqtr = median(YrQtr),
            yrs.fished = mean(yrs.fished)) %>% ungroup()

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

temporal.df %>% filter(yrs.fished > 7) %>%
  ggplot() +
  geom_point(aes(x = mean.yrqtr, y = yrs.fished)) +
  labs(y = "Years Fished", x = "Midpoint of Career") +
  gg.theme

## select vessels with > 7 years of fishing to get more obs per vessel
hist(temporal.df$yrs.fished)
core.dat <- ops.filt %>% left_join(temporal.df)
## classify flags into DW, PICT, or AU/NZ
DW.flags <- c("BZ", "ES", "SU", "TW", "CN", "KR", "JP", "US", "VN", "VU")
PICT.flags <- c("AU", "NZ", "NC", "FJ", "PF", "SB", "TO", "PG", "CK", "WS",
                "NU", "KI", "TV", "ID", "FM")
core.dat %<>% mutate(Fleet.type = ifelse(flag_id %in% DW.flags, "DW",
                                         ifelse(flag_id %in% PICT.flags,
                                                "PICT", NA)))
table(core.dat$Fleet.type)

## create yrqtr.cell.5x5
core.dat %<>%
  mutate(yrqtr.cell.5x5 = paste0(YrQtr, "-", cell.5x5)) 

## see just how poor the sampling coverage is
yrqtr.cells <- core.dat %>% filter(!is.na(vessel_id)) %>%
  select(yrqtr.cell.5x5) %>%
  mutate(freq = 1) %>%
  group_by(yrqtr.cell.5x5) %>% summarise(freq = sum(freq))

ts <- length(unique(core.dat$YrQtr)) ## number of time steps

## plot all the cells in the dataset
unique.pts <- core.dat %>% select(cell.5x5) %>% distinct() %>%
  tidyr::separate(cell.5x5, into =  c("lon", "lat"), sep = "x", remove = F) %>%
  mutate(lat = as.numeric(lat) + 2.5, lon = as.numeric(lon) + 2.5)
coordinates(unique.pts) <- ~lon+lat
unique.pts@proj4string <- crs_ll
unique.pts %<>% st_as_sf(unique.pts)

## remove cells that get less than an average of min threshold per time step
cells <- core.dat %>% filter(!is.na(vessel_id)) %>%
  select(cell.5x5) %>% mutate(freq = 1) %>%
  group_by(cell.5x5) %>% summarise(freq = sum(freq)) %>%
  ungroup() %>% filter(freq >= min.threshold * ts) %>% select(-freq) %>%
  tidyr::separate(cell.5x5, into =  c("lon", "lat"), sep = "x", remove = F) %>%
  mutate(lat = as.numeric(lat) + 2.5, lon = as.numeric(lon) + 2.5)

new.cells <- cells %>% select(lat, lon) %>% distinct()
coordinates(new.cells) <- ~lon+lat
new.cells@proj4string <- crs_ll
new.cells %<>% st_as_sf(new.cells)

ggplot() +
  geom_sf(data = coast.sf) +
  geom_text(data = axis_labels, aes(x = long, y = lat, label = labels),
            size = 5) +
  geom_segment(aes(x = 130 , y = seq(-50, 0, 10), xend = 287.5,
                   yend = seq(-50, 0, 10)), col = "lightgrey") +
  geom_segment(aes(x = seq(140, 280, 20), y = -50, xend = seq(140, 280, 20),
                   yend = 0), col = "lightgrey") +
  geom_sf(data = regions.sf, fill = NA, col = "blue") +
  geom_sf(data = unique.pts, shape = 15, size = 7.5, col = "cyan") +
  geom_sf(data = new.cells, shape = 15, size = 7, col = "darkred") +
  theme(axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = NA, color = "white"))

#ggsave('../Figures/all.cells.in.data.png', width=10, height=5,
#       units='in', dpi=300)
graphics.off()

## remove cells that are not sampled adequately
## over time (on average, min.threshold for each time step)
yrqtr.cells <- core.dat %>% select(yrqtr.cell.5x5, cell.5x5) %>%
  left_join(cells %>% select(cell.5x5) %>% mutate(freq = 1)) %>%
  filter(freq == 1) %>% distinct() %>%
  select(yrqtr.cell.5x5)
## ensure the yrqtr.cells have reasonable number of obs to select from in data
yrqtr.cells <- core.dat %>% select(yrqtr.cell.5x5, cell.5x5) %>%
  left_join(yrqtr.cells %>% mutate(freq = 1)) %>% filter(freq == 1) %>%
  group_by(yrqtr.cell.5x5, cell.5x5) %>% summarise(freq = sum(freq)) %>%
  filter(freq >= 10) %>% ungroup() %>% select(yrqtr.cell.5x5)  %>% distinct()
  
num.cells <- nrow(yrqtr.cells)
vessels <- core.dat %>% filter(yrs.fished >= min.yrs.fished) %>%
  select(vessel_id) %>% distinct() %>% unlist()
vessels <- vessels[!is.na(vessels)]
sample.coverage.final <- 0
vessels.final <- NA

max.vessels <- as.integer(0.1 * length(vessels)) ## max number of vessels in core fleet

seeds <- c(44505843, 88892444, 88147626, 96279512, 14096433, 15551159, 81062538,
           12216060, 32127134, 61515742)

dir.create("../Results/vessel/")

for(j in 1:length(seeds)){
  vessels.final <- NA
  sample.coverage.final <- 0
  for(i in 1:length(vessels)){
    vess <- sample(vessels, size = 1)
    vessels <- vessels[vessels != vess]
    if(i == 1){
      vessels.final <- vess
    } else {
      vessels.final <- append(vessels.final, vess)
    }
    
    df <- core.dat %>% filter(vessel_id %in% vessels.final) %>%
      mutate(freq = 1) %>%
      group_by(yrqtr.cell.5x5) %>% summarise(freq = sum(freq))
    df %<>% right_join(yrqtr.cells) %>%
      mutate(freq = ifelse(is.na(freq), 0, freq)) %>%
      mutate(sample.num = ifelse(freq >= (min.threshold),
                                 (min.threshold), freq))
    
    sample.coverage <- (sum(df$sample.num) / ((min.threshold) * num.cells))
    if(sample.coverage > sample.coverage.final){
      sample.coverage.final <- sample.coverage
    } else {
      vessels.final <- setdiff(vessels.final, vess)
    }
    
    num.vessels <- length(vessels.final)
    
    if(sample.coverage.final == 1 | num.vessels == max.vessels){break}
  }

  vessels.final %<>% data.frame()
  colnames(vessels.final) <- "vessel_id"
  
  save(vessels.final, file = paste0("../Results/vessel/core_vessels", j,
                                    ".RData"))
}
