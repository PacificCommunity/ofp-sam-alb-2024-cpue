library(VAST)
library(tidyverse)
library(magrittr)

# Remove all previous objects
rm(list=ls())

gg.theme <- theme_bw() +
  theme(axis.line = element_line(color="black"),
        panel.grid.minor = element_blank(),
        text = element_text(size=18),
        axis.title = element_text(size=20),
        strip.background =element_rect(fill="slategray"),
        strip.text = element_text(size=12, colour='white'),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        panel.background = element_rect(fill = NA, color = "black")) 

setwd('E:/ALB_CPUE/2024')

## load 2021 regional structure
load(file = "Background_Data/alb_regions_poly_2024_shp.RData")
regions.2024 <- regions.shp; rm(regions.shp)

# Read in filtered data - from filter_data.r
load('Data/ops.model.Rdata')

# Correct logdate month, quarter - as format changed for original file
ops.filt %<>% mutate(month = as.numeric(substr(logdate,4,5)),
                     quarter=(ceiling(month/3)/4)-0.25, YrQtr = year+quarter )
ops.filt %<>% filter(year<=2022)

## summarise catch by region (2021 4 region structure)
catch.summary <- ops.filt %>% select(latd, lond, alb_n, YrQtr)

## assign regions based on 2021 regional structure
num.regions <- 2 ## number of regions in 2024
catch.summary$Region = NA
for(i in 1:num.regions){
  coords = regions.2024@polygons[[i]]@Polygons[[1]]@coords
  catch.summary$Region = ifelse(sp::point.in.polygon(catch.summary$latd,
                                                      catch.summary$lond, coords[,2],
                                                      coords[,1]) %in%
                                   c(1,2),i,catch.summary$Region)
}
catch.summary %<>% mutate(Region = ifelse(Region == 1 & latd <= -25, "1-CD",
                                          ifelse(Region == 1 & latd > -25,
                                                 "1-AB", Region))) %>%
  mutate(Region = paste("Region", Region)) %>% mutate(year = floor(YrQtr)) %>%
  mutate(qtr =  4 * (YrQtr - year) + 1) %>% group_by(Region, year, qtr) %>%
  summarise(alb_n = sum(alb_n)) %>% ungroup() %>% filter(Region != "Region NA")

## plot catch time series by region
catch.summary %>% mutate(alb_n = alb_n / 10000) %>% ggplot() +
  geom_line(aes(year, alb_n, col = as.factor(qtr)), size = 1) +
  labs(y = 'Catch (x 10,000)', x = 'Year', col = "Quarter")+
  theme_bw() + gg.theme +
  facet_grid(Region~.)

ggsave('Figures/catch.time.series.by.quarter.png',
       width=10, height=7, units='in', dpi=300)

## plot effort by quarter
## summarise effort by region (2021 4 region structure)
effort.summary <- ops.filt %>% select(latd, lond, hhook, YrQtr)

## assign regions based on 2021 regional structure
num.regions <- 2 ## number of regions in 2024
effort.summary$Region = NA
for(i in 1:num.regions){
  coords = regions.2024@polygons[[i]]@Polygons[[1]]@coords
  effort.summary$Region = ifelse(sp::point.in.polygon(effort.summary$latd,
                                                     effort.summary$lond, coords[,2],
                                                     coords[,1]) %in%
                                  c(1,2),i,effort.summary$Region)
}
effort.summary %<>% mutate(Region = ifelse(Region == 1 & latd <= -25, "1-CD",
                                          ifelse(Region == 1 & latd > -25,
                                                 "1-AB", Region))) %>%
  mutate(Region = paste("Region", Region)) %>% mutate(year = floor(YrQtr)) %>%
  mutate(qtr =  4 * (YrQtr - year) + 1) %>% group_by(Region, year, qtr) %>%
  summarise(hhook = sum(hhook)) %>% ungroup() %>% filter(Region != "Region NA")

## plot effort time series by region
effort.summary %>% mutate(hhook = hhook / 10000) %>%
  ggplot() +
  geom_line(aes(year, hhook, col = as.factor(qtr)), size = 1) +
  labs(y = 'Hooks (x 1,000,000)', x = 'Year', col = "Quarter")+
  theme_bw() + gg.theme +
  facet_grid(Region~.)

ggsave('Figures/effort.time.series.by.quarter.png',
       width=10, height=7, units='in', dpi=300)
