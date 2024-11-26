## compare indices from separate model runs

library(tidyverse)
library(magrittr)

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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

## need to set year-quarter time step
yr.start <- 1960 ## start year-quarter in model
yr.end <- 2022.75 ## last year-quarter in model
time.step <- 0.25 ## quarterly time step
ts.df <- data.frame(YrQtr = seq(yr.start, yr.end, time.step))
ts.df %<>% mutate(ts =  seq(1, nrow(ts.df), 1))

## load 2021 regional structure
load(file = "../Background_Data/alb_regions_poly_2021_shp.RData")
regions.2021 <- regions.shp
## load 2024 regional structure
load(file = "../Background_Data/alb_regions_poly_2024_shp.RData")
regions.2024 <- regions.shp
rm(regions.shp)
############################
## Load results
############################
setwd("E:/2021_SC/CPUE/ALB/Results/")

## Load 2021 results
load("seed123_234760_knots_150_1960_2019_targ3H__flag_17082021.Rdata")
orig.2021 = out[[1]]
orig.2021.dat = out[[2]]

## Get 2021 index
Index.2021 =summary(orig.2021$parameter_estimates$SD,
                    select='report')[grep('^Index_ctl$',
                                          rownames(summary(orig.2021$parameter_estimates$SD,
                                                           select='report')) ), ] 
Index.2021 %<>% data.frame() %>%
  mutate(Region = paste0('Region ',rep(c(1:4, "All"), each = nrow(Index.2021)/5)),
         YrQtr = rep(sort(unique(orig.2021.dat$YrQtr)),times=5)) %>%
  rename(SE=Std..Error)

## Index 2021
ggplot(Index.2021  %>%  mutate(lwr = Estimate-SE, upr = Estimate + SE),
       aes(YrQtr, Estimate)) + geom_line(col='lightgray') +
  geom_ribbon(aes(ymin=lwr,ymax=upr),col="white",alpha=0.3) +
  geom_point(size=1) +
  ylab('Relative abundance (1000 mt)') + xlab('Year-quarter')+
  theme_bw() + gg.theme + theme(legend.position='none') + 
  facet_grid(Region~., scales='free')

## load VAST results for updated data
setwd("E:/ALB_CPUE/2024/")
load('Results/2024_alb_fit_1954_annual.Rdata')
new.2024 = out[[1]]
new.2024.dat = out[[2]]

## Get index
Index.2024 =summary(new.2024$parameter_estimates$SD,
                    select='report')[grep('^Index_ctl$',
                                          rownames(summary(new.2024$parameter_estimates$SD,
                                                           select='report')) ), ] 
Index.2024 %<>% data.frame() %>%
  mutate(Region = paste0('Region ',rep(c(1:2), each = nrow(Index.2024)/2)),
         year = rep(sort(unique(new.2024.dat$year)),times=2)) %>%
  rename(SE=Std..Error)

## Index 2024
ggplot(Index.2024  %>%
         mutate(lwr = Estimate-SE, upr = Estimate + SE),
       aes(year, Estimate)) + geom_line(col='lightgray') +
  geom_ribbon(aes(ymin=lwr,ymax=upr),col="white",alpha=0.3) +
  geom_point(size=1) +
  ylab('Relative abundance (1000 mt)') + xlab('Year')+
  theme_bw() + gg.theme + theme(legend.position='none') + 
  facet_grid(Region~., scales='free')

## split new indices into 4 regions to be comparable with 2021
## get estimates by 1x1 degree grid cells
grid.df <- new.2024$spatial_list$latlon_g %>% data.frame() %>%
  mutate(Cell = seq(1:nrow(new.2024$spatial_list$latlon_g)))

new.regions <- 2 ##number of regions in new regional structure
for(i in 1:new.regions){
  df <- new.2024$Report$Index_gctl[,,,(i)] %>% data.frame()
  colnames(df) <- seq(1, ncol(df))
  df %<>% mutate(Cell = seq(1:nrow(new.2024$spatial_list$latlon_g))) %>%
    gather(key = "ts", value = "rel_abundance", -Cell) %>%
    mutate(ts = as.integer(ts))
  df %<>% left_join(grid.df)
  
  if(!exists("newdf")){
    newdf <- df
  } else {
    newdf <- rbind(newdf, df)
  }
  rm(df)
}

Index.1x1 <- newdf %>% group_by(ts, Lat, Lon) %>%
  summarise(rel_abundance = sum(rel_abundance)) %>% data.frame()
rm(newdf)

sp::plot(Index.1x1 %>% select(Lon, Lat) %>% distinct(), col = "blue", pch = 16)
sp::plot(regions.2024, add = T)
graphics.off()

## convert the cell size to 5x5
ji <- function(xy, origin=c(0,0), cellsize=c(5,5)) {
  t(apply(xy, 1, function(z) cellsize/2+origin+
            cellsize*(floor((z - origin)/cellsize))))
}
JI <- ji(cbind(Index.1x1$Lon, Index.1x1$Lat), cellsize = c(5, 5))
Index.1x1$lon_5 <- JI[, 1]
Index.1x1$lat_5 <- JI[, 2]

## need to set year time step
yr.start.new <- 1954 ## start year in model
yr.end.new <- 2022 ## last year in model
time.step.new <- 1 ## annual time step
ts.df.new <- data.frame(year = seq(yr.start.new, yr.end.new, time.step.new))
ts.df.new %<>% mutate(ts =  seq(1, nrow(ts.df.new), 1))
Index.5x5 <- Index.1x1 %>% select(-c(Lat, Lon)) %>%
  group_by(ts, lat_5, lon_5) %>%
  summarise(rel_abundance = sum(rel_abundance)) %>% ungroup() %>%
  left_join(ts.df.new)

## verify it matches with the index calculated in new.2024$Report
## assign regions based on 2024 regional structure
Index.5x5$Region.2024 = NA
for(i in 1:new.regions){
  coords = regions.2024@polygons[[i]]@Polygons[[1]]@coords
  Index.5x5$Region.2024 = ifelse(sp::point.in.polygon(Index.5x5$lat_5,
                                               Index.5x5$lon_5, coords[,2],
                                               coords[,1]) %in%
                            c(1,2),i,Index.5x5$Region.2024)
}
Index.5x5 %<>% mutate(Region.2024 = paste("Region", Region.2024))

## check they agree
ggplot() +
  geom_line(data = Index.5x5 %>% group_by(Region.2024, year) %>%
              summarise(Estimate = sum(rel_abundance)),
            aes(year, Estimate), col='blue') +
  geom_line(data = Index.2024 %>% rename(Region.2024 = Region),
            aes(year, Estimate), col='red') +
  ylab('Relative abundance (1000 mt)') + xlab('Year')+
  theme_bw() + gg.theme + theme(legend.position='none') + 
  facet_grid(Region.2024~., scales='free')

## assign regions based on 2021 regional structure
num.regions <- 4 ## number of regions in 2021
Index.5x5$Region.2021 = NA
for(i in 1:num.regions){
  coords = regions.2021@polygons[[i]]@Polygons[[1]]@coords
  Index.5x5$Region.2021 = ifelse(sp::point.in.polygon(Index.5x5$lat_5,
                                                      Index.5x5$lon_5, coords[,2],
                                                      coords[,1]) %in%
                                   c(1,2),i,Index.5x5$Region.2021)
}
Index.5x5 %<>% mutate(Region.2021 = paste("Region", Region.2021))

## compare with 2021 index

combined.2021 <- Index.5x5 %>% group_by(Region.2021, year) %>%
  summarise(Index = sum(rel_abundance)) %>% ungroup() %>%
  group_by(Region.2021) %>%
  mutate(Index = Index / mean(Index)) %>% mutate(model = "VAST 2024") %>%
  rename(YrQtr = year) %>%
  rbind(Index.2021 %>% filter(Region != "Region All") %>%
          rename(Region.2021 = Region)  %>%
          group_by(Region.2021) %>%
          mutate(Index = Estimate / mean(Estimate)) %>%
          mutate(model = "VAST 2021") %>% select(-c(SE, Estimate)))

cols <- c("salmon", "deepskyblue3")
combined.2021 %>% ggplot() +
  geom_line(aes(YrQtr, Index, col=model), size = 1) +
  geom_point(aes(YrQtr, Index, col=model), size = 1) +
  labs(y = 'Standardized relative abundance (mean-centered)', x = 'Year',
       col = "CPUE model") +
  scale_colour_manual(values = cols) +
  theme_bw() + gg.theme +
  facet_grid(Region.2021~., scales='free')

ggsave('Figures/MC.Std.Indices.2021.updated.new.data.annual.png',
       width=10, height=7, units='in', dpi=300)

## save final 5x5 results for Tom Peatman
#write.csv(Index.5x5 %>% select(-c(Region.2021, Region.2024)),
#          file = "Results/2024.alb.vast.fit.annual.csv", row.names = F)

## rescale indices for EPO and WCPO for Tom Peatman's LF reweigting
#EPO <- read.csv(file = "Results/EPO.quarterly.1954.data2024.index.5x5.csv")
#WCPO <- read.csv(file = "Results/WCPO.quarterly.1954.data2024.index.5x5.csv")

#VAST.2024.scale <- Index.2024 %>% group_by(Region, year) %>%
#  summarise(mean.est = mean(Estimate)) %>% ungroup() %>%
#  pivot_wider(values_from = "mean.est", names_from = "Region") %>%
#  mutate(scale = `Region 2` / `Region 1`) %>% select(year, scale)

#EPO %<>% mutate(year = floor(YrQtr)) %>% left_join(VAST.2024.scale) %>%
#  mutate(rel_abund = rel_abund * scale)

#PO <- WCPO %>% rename(WCPO.rel.abund = rel_abund) %>%
#  full_join(EPO %>% rename(EPO.rel.abund = rel_abund)) %>%
#  mutate(diff = WCPO.rel.abund - EPO.rel.abund)

## take a look to see if it looks right
#PO %>% filter(YrQtr == 2005.5) %>%
#  mutate(EPO.rel.abund = ifelse(!is.na(WCPO.rel.abund), 0, EPO.rel.abund)) %>%
#  mutate(EPO.rel.abund = ifelse(is.na(EPO.rel.abund), 0, EPO.rel.abund)) %>%
#  mutate(WCPO.rel.abund = ifelse(is.na(WCPO.rel.abund), 0, WCPO.rel.abund)) %>%
#  mutate(rel_abund = WCPO.rel.abund + EPO.rel.abund) %>%
#  ggplot() +
#  geom_tile(aes(x = lon_5, y = lat_5, fill = rel_abund)) +
#  gg.theme

#PO %<>%
#  mutate(EPO.rel.abund = ifelse(!is.na(WCPO.rel.abund), 0, EPO.rel.abund)) %>%
#  mutate(EPO.rel.abund = ifelse(is.na(EPO.rel.abund), 0, EPO.rel.abund)) %>%
#  mutate(WCPO.rel.abund = ifelse(is.na(WCPO.rel.abund), 0, WCPO.rel.abund)) %>%
#  mutate(rel_abund = WCPO.rel.abund + EPO.rel.abund) %>% distinct() %>%
#  select(lat_5, lon_5, YrQtr, rel_abund)

#write.csv(PO,
#          file = "Results/ALB.quarterly.1954.data2024.index.5x5.rescaled.csv",
#          row.names = F)

## check that overall trends are the same
## assign regions based on 2024 regional structure
#PO$Region.2024 = NA
#for(i in 1:new.regions){
#  coords = regions.2024@polygons[[i]]@Polygons[[1]]@coords
#  PO$Region.2024 = ifelse(sp::point.in.polygon(PO$lat_5,
#                                                      PO$lon_5, coords[,2],
#                                                      coords[,1]) %in%
#                                   c(1,2),i,PO$Region.2024)
#}
#PO %<>% mutate(Region.2024 = paste("Region", Region.2024))
## compare with VAST model
#comp.vast <- PO %>% group_by(Region.2024, YrQtr) %>%
#  summarise(Index = sum(rel_abund)) %>% ungroup() %>%
#  mutate(model = "sdmTMB") %>%
#  rbind(combined.2021 %>% filter(model == "VAST 2021") %>%
#          mutate(Region.2024 = ifelse(Region.2021 == "Region 4", "Region 2",
#                                      "Region 1")) %>%
#          group_by(Region.2024, YrQtr, model) %>% summarise(Index = sum(Index)) %>%
#          ungroup())

## check they agree
#comp.vast %>% ggplot() +
#  geom_line(aes(YrQtr, Index)) +
#  ylab('Relative abundance (1000 mt)') + xlab('Year')+
#  theme_bw() + gg.theme + theme(legend.position='none') + 
#  facet_grid(model ~ Region.2024, scales='free')

## add year and qtr variables
#Index.2021 %<>% mutate(year = floor(YrQtr)) %>%
#  mutate(qtr = 4 * (YrQtr - year) + 1)

## check correlations between quarters
#rm(cor.df.final)
#regs <- seq(1, num.regions, 1)
#qtrs <- seq(1, 4, 1)
#cor.df <- data.frame(matrix(rep(NA, (length(qtrs) * length(qtrs))),
#             ncol = length(qtrs), nrow = length(qtrs)))
#colnames(cor.df) <- qtrs
#cor.df$Qtr1 <- qtrs

#for(i in 1:length(regs)){
#  cor.df.new <- cor.df
#  for(j in 1:length(qtrs)){
#    for(k in 1:length(qtrs)){
#      if(k > j){
#        ind1 <- Index.2021 %>%
#          filter(qtr == qtrs[j], Region == paste("Region", regs[i])) %>%
#          select(Estimate) %>% unlist() %>% unname()
#        ind2 <- Index.2021 %>%
#          filter(qtr == qtrs[k], Region == paste("Region", regs[i])) %>%
#          select(Estimate) %>% unlist() %>% unname()
#        cor.df.new[j, k] <- round(cor(ind1, ind2), 2)
#      }
#    }
#  }
#  cor.df.new %<>% pivot_longer(cols = !"Qtr1", names_to = "Qtr2",
#                               values_to = "Correlation") %>%
#    mutate(Region = i) %>% filter(!is.na(Correlation))
  
#  if(!exists("cor.df.final")){
#    cor.df.final <- cor.df.new
#  } else {
#    cor.df.final %<>% rbind(cor.df.new)
#  }
#}

#cor.df.final %<>% filter(Correlation < 0.8)

#write.csv(cor.df.final, file = "Results/VAST.2021.correlations.by.qtr.csv",
#          row.names = F)

## plot cpue by quarter
#ind.exp <- Index.2021 %>% select(Estimate, Region, year, qtr) %>%
#  pivot_wider(names_from = "qtr", values_from = "Estimate") %>%
#  filter(Region != "Region All")
#
#cols <- c("1" = "darkred", "2" = "navy",
#          "3" = "darkgreen", "4" = "darkgoldenrod")
#line.size <- 1

#ggplot() +
#  facet_wrap(~Region, ncol = 1, scales = "free", strip.position = "right") +
#  geom_line(data = ind.exp, aes(x = year, y = `1`, col = "1"),
#            size = line.size) +
#  geom_line(data = ind.exp, aes(x = year, y = `2`, col = "2"),
#            size = line.size) +
#  geom_line(data = ind.exp, aes(x = year, y = `3`, col = "3"),
#            size = line.size) +
#  geom_line(data = ind.exp, aes(x = year, y = `4`, col = "4"),
#            size = line.size) +
#  scale_colour_manual(name="Quarter", values = cols) +
#  labs(x = "Year", y = "Relative Abundance") +
#  gg.theme
  
#ggsave('Figures/indices.by.quarter.png',
#       width=10, height=7, units='in', dpi=300)

## check correlations between regions
#rm(cor.df.final)
#cor.df <- data.frame(matrix(rep(NA, (length(regs) * length(regs))),
#                            ncol = length(regs), nrow = length(regs)))
#colnames(cor.df) <- regs
#cor.df$Region1 <- regs

#for(i in 1:length(qtrs)){
#  cor.df.new <- cor.df
#  for(j in 1:length(regs)){
#    for(k in 1:length(regs)){
#      if(k > j){
#        ind1 <- Index.2021 %>%
#          filter(qtr == qtrs[i], Region == paste("Region", regs[j])) %>%
#          select(Estimate) %>% unlist() %>% unname()
#        ind2 <- Index.2021 %>%
#          filter(qtr == qtrs[i], Region == paste("Region", regs[k])) %>%
#          select(Estimate) %>% unlist() %>% unname()
        
#        cor.df.new[j, k] <- round(cor(ind1, ind2), 2)
#      }
#    }
#  }
#  cor.df.new %<>% pivot_longer(cols = !"Region1", names_to = "Region2",
#                               values_to = "Correlation") %>%
#    mutate(Qtr = i) %>% filter(!is.na(Correlation))
  
#  if(!exists("cor.df.final")){
#    cor.df.final <- cor.df.new
#  } else {
#    cor.df.final %<>% rbind(cor.df.new)
#  }
#}

#cor.df.final %<>% filter(Correlation < 0.8)

#write.csv(cor.df.final, file = "Results/VAST.2021.correlations.by.region.csv",
#          row.names = F)

## plot cpue by quarter
#ind.exp <- Index.2021 %>% filter(Region != "Region All") %>%
#  group_by(Region) %>%
#  mutate(Estimate = Estimate / mean(Estimate)) %>% ungroup() %>%
#  select(Estimate, Region, year, qtr) %>%
#  pivot_wider(names_from = "Region", values_from = "Estimate") %>%
#  mutate(qtr = paste("Qtr", qtr))

#ggplot() +
#  facet_wrap(~qtr, ncol = 1, scales = "free", strip.position = "right") +
#  geom_line(data = ind.exp, aes(x = year, y = `Region 1`, col = "1"),
#            size = line.size) +
#  geom_line(data = ind.exp, aes(x = year, y = `Region 2`, col = "2"),
#            size = line.size) +
#  geom_line(data = ind.exp, aes(x = year, y = `Region 3`, col = "3"),
#            size = line.size) +
#  geom_line(data = ind.exp, aes(x = year, y = `Region 4`, col = "4"),
#            size = line.size) +
#  scale_colour_manual(name="Region", values = cols) +
#  labs(x = "Year", y = "Mean-centered Relative Abundance") +
#  gg.theme

#ggsave('Figures/indices.by.region.png',
#       width=10, height=7, units='in', dpi=300)
#graphics.off()

## compare sdmTMB annual indices with and without season
load(file = "Results/index.annual.1954.season.data2024.RData")
index.season <- index.final; rm(index.final)
load(file = "Results/index.annual.1954.data2024.RData")

index.season %<>% mutate(model = "Seasonal.catchability") %>%
  rbind(index.final %>% mutate(model = "Non-seasonal"))
rm(index.final)

load(file = "Results/index.annual.1954.season.as.density.data2024.RData")
index.season %<>%
  rbind(index.final %>% mutate(model = "Seasonal.density"))
rm(index.final)

cols <- c("salmon", "deepskyblue3", "darkgoldenrod")
index.season %>% ggplot() +
  geom_line(aes(year, rel_abund, col=Region), size = 1) +
  geom_point(aes(year, rel_abund, col=Region), size = 1) +
  labs(y = 'Standardized relative abundance', x = 'Year',
       col = "Region") +
  scale_colour_manual(values = cols) +
  theme_bw() + gg.theme +
  facet_wrap(~model, scales='free', nrow = 1) +
  theme(axis.text.y = element_blank())

#ggsave('Figures/Std.Indices.2024.seasonal.new.data.annual.png',
#       width=12, height=7, units='in', dpi=300)

index.season %>%
  group_by(Region, model) %>%
  mutate(Index = rel_abund / mean(rel_abund)) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(year, Index, col=model), size = 1) +
  geom_point(aes(year, Index, col=model), size = 1) +
  labs(y = 'Standardized relative abundance (mean-centered)', x = 'Year',
       col = "CPUE model") +
  scale_colour_manual(values = cols) +
  theme_bw() + gg.theme +
  facet_grid(Region~., scales='free')

#ggsave('Figures/MC.Std.Indices.2024.seasonal.new.data.annual.png',
#       width=10, height=7, units='in', dpi=300)

prop.df <- index.season %>% filter(Region == "1-AB") %>% select(-Region) %>% 
  rename(index1AB = rel_abund) %>%
  left_join(index.season %>% filter(Region == "1-CD") %>% select(-Region) %>%
              rename(index1CD = rel_abund)) %>%
  left_join(index.season %>% filter(Region == "2") %>% select(-Region) %>%
              rename(index2 = rel_abund)) %>%
  mutate("1AB-1CD" = index1AB / index1CD,
         "1AB-2" = index1AB / index2,
         "1CD-2" = index1CD / index2) %>%
  select(-c(index1AB, index1CD, index2)) %>%
  pivot_longer(cols = starts_with("1"), names_to = "regions",
               values_to = "proportion")

prop.df %>% ggplot() +
  geom_line(aes(year, proportion, col=model), size = 1) +
  geom_point(aes(year, proportion, col=model), size = 1) +
  labs(y = 'Proportion of Relative Abundance between Regions',
       x = 'Year', col = "CPUE model") +
  scale_colour_manual(values = cols) +
  facet_wrap(~regions, ncol = 1) +
  theme_bw() + gg.theme

#ggsave('Figures/Regional.ratio.Indices.2024.seasonal.new.data.annual.png',
#       width=10, height=7, units='in', dpi=300)

## compare quarterly indices with sps cluster and hbf
WCPO.sps.clust <- read.csv(file = "Results/index.cv.WCPO.quarterly.1954.data2024.csv")
EPO.sps.clust <- read.csv(file = "Results/index.cv.EPO.quarterly.1954.data2024.csv")
index.sps.clust <- rbind(EPO.sps.clust %>% filter(Region == 2),
                         WCPO.sps.clust %>% filter(Region != 2))
rm(EPO.sps.clust, WCPO.sps.clust)
load(file = "Results/index.HBF.quarterly.1954.data2024.RData")
index.hbf <- index.final
rm(index.final)

load(file = "Results/index.HBF.lat7.quarterly.1954.data2024.RData")
index.hbf.lat7 <- index.final
rm(index.final)

ind.comb <- rbind(index.hbf %>% rename(Index = rel_abund) %>%
                    mutate(model = "hbf"),
                  index.hbf.lat7 %>% rename(Index = rel_abund) %>%
                    mutate(model = "hbf.lat7"),
                  index.sps.clust %>% select(-cv) %>%
                    mutate(model = "sps.clust")) %>%
  group_by(Region, model) %>%
  mutate(Index = Index / mean(Index)) %>% ungroup()

## compare results
ggplot() +
  geom_line(data = ind.comb,
            aes(YrQtr, Index, col = model), size = 0.5) +
  ylab('Std. Mean-centered Relative Abundance') + xlab('YrQtr')+
  scale_color_manual(values = c("deepskyblue3", "brown", "chartreuse3")) +
  theme_bw() + gg.theme +
  facet_grid(Region~., scales='free')

ggsave('Figures/MC.Std.Indices.hbf.sps.clust.png',
       width=10, height=7, units='in', dpi=300)
