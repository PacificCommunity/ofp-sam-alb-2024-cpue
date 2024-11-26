library(data.table)
library(tidyverse)
library(magrittr)
library(randomForest)
library(sp)

gg.theme <- theme_bw() + 
  theme(axis.line = element_line(color="black"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=14),
        axis.title = element_text(size=16),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"))  

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

############################################
setwd('E:/ALB_CPUE/2024/')
############################################

# Load data from: prep.ops.data.R
load('Data/ops.trim.dt.RData')
names(ops.trim.dt)

# Filtering is based off the work done in 2018 
#   Repeating for 2024 analysis

ops.filt = ops.trim.dt[latd>=-50 & latd<=0] # Spatial filter 0 - -50S
ops.filt = ops.filt[lond>=130 ] # Spatial filter > 140  Last assessment: 140 -230
#ops.filt %<>% mutate(drop.ll = ifelse(lond>210 & latd>-5,1,0)) %>% filter(drop.ll==0) %>% select(-drop.ll)

ops.filt = ops.filt[!is.na(logdate)] # No missing logdate
ops.filt = ops.filt[hook!=0] # Hooks >0
ops.filt = ops.filt[year>=1954] # Year >= 1954
ops.filt = ops.filt[alb_cpue>0 | bet_cpue>0 | yft_cpue>0 | swo_cpue>0] #At least 1 fish of one of the 4 main species 


# Check on missing vessel id
#source('./RCode/check.missing.vessel.R')

# Identify strata with insufficient data
# Must have at least 50 sets per year.qtr and 20 sets per 5x5 cell
ops.filt %<>% group_by(yr.qtr) %>% mutate(drop.yr.qtr = ifelse(n()<50,1,0)) %>% group_by(cell.5x5) %>% mutate(drop.cell.5x5 = ifelse(n()<20,1,0))
ops.filt %<>% filter(drop.cell.5x5==0 & drop.yr.qtr==0) # Keep only yr-qtrs with 50+ sets; 5x5 cells with 20+ sets 
ops.filt %<>% select(-c(drop.yr.qtr, drop.cell.5x5))

ops.filt %<>% mutate(YrQtr = year + quarter/4 - 0.25)
ops.filt %<>% mutate(date = as.Date(logdate, format='%m/%d/%Y'))

# Retain only JP in the EPO
load("Background_Data/alb_regions_poly_2024_shp.Rdata")

ops.filt$Region = NA
for(i in 1:2){
  coords = regions.shp@polygons[[i]]@Polygons[[1]]@coords
  ops.filt$Region = ifelse(point.in.polygon(ops.filt$lond, ops.filt$latd, coords[,1], coords[,2]) %in% c(1,2),i,ops.filt$Region)
}

## this was done in 2021 but not done in 2024 because the Japanese data may not
## be very representative of albacore abundance due to targeting issues even
## though it has a long time series
#ops.filt %<>% mutate(drop = ifelse(Region %in% c(4) & flag_id!='JP',1,0)) %>% filter(drop==0) %>% select(-drop)

rm(coords, ops.trim.dt, regions)

##################################################
# Predict missing HBF values 
##################################################
# Turn 0s into NAs and HBF bins into 5 hook bins

ops.filt %<>% mutate(hbf = ifelse(hk_bt_flt==0,NA,hk_bt_flt)) %>% ungroup()
ops.filt$hbf = factor(ceiling(ops.filt$hbf/5)*5, levels= sort(unique(ceiling(ops.filt$hbf/5)*5)) )
ops.filt$hbf.label = as.numeric(ops.filt$hbf)-1

ops.filt %<>% mutate(total.catch = alb_n + bet_n + yft_n + swo_n )
ops.filt %<>% mutate(alb =  alb_n/total.catch,
                            bet = bet_n/total.catch,
                            yft = yft_n/total.catch,
                            swo = swo_n/total.catch )

hbf.good = ops.filt %>% filter(!is.na(hbf)) %>% select(hbf,year,flag_id, month, lond, latd, hhook,total.catch, alb, bet, yft, swo)
hbf.miss =  ops.filt %>% filter(is.na(hbf)) %>% select(hbf,year,flag_id, month, lond, latd, hhook,total.catch, alb, bet, yft, swo)

N =ceiling(nrow(hbf.good)/length(unique(hbf.good$hbf)))*.33
setDT(hbf.good)
hist(as.numeric(as.character(hbf.good$hbf)))

# Repeat with 10 bootstrap samples
pred.hbf = list()

for(i in 1:10){

  train = hbf.good[, .SD[sample(.N, min(N,.N))],by = hbf] %>% data.frame()
  
  train = sample_n(hbf.good, round(1000000)) # about 50% of all samples?
  
  rf = randomForest(hbf~., data=train )
  rf.pred = predict(rf, hbf.miss, predict.all=TRUE)
  
  pred.hbf[[i]] = hbf.miss %>% mutate(pred = rf.pred$aggregate)
  rm(rf, rf.pred)
}

save(pred.hbf, file = "Data/pred.hbf.RData")

##########################################
# Get estimates of HBF and compare
##########################################
hbf.estimates = matrix()

for (i in 1:10){
  hbf.estimates = cbind(hbf.estimates, pred.hbf[[i]][,'pred']  )
}
hbf.estimates = hbf.estimates[,-1] %>% data.frame()

hbf.estimates %<>% mutate(row = 1:nrow(.)) %>% gather(key='iter', value = 'hbf',-row) %>% spread(key=iter, value=hbf)
hbf.estimates$mode = apply(hbf.estimates[,-1], 1,mode)

rm(hbf.good, hbf.miss)

ggplot(hbf.estimates %>% select(-row) %>% gather(key='iter', value='hbf', -mode) , 
       aes(factor(hbf, levels=seq(5,50,5)), fill=factor(mode, levels=seq(5,50,5)) )) +geom_density() +  
  labs(fill='Mode') + gg.theme + xlab('HBF') + ylab('Density') + scale_fill_brewer(palette='Spectral')

ggsave('Figures/hbf_fill_density.png', width=12, height=8, units='in')

##########################################
# Fill in data set
##########################################
tmp = ops.filt 
tmp %<>% mutate(hbf = ifelse(hk_bt_flt==0,NA,hk_bt_flt)) %>% ungroup()
tmp$hbf = factor(ceiling(tmp$hbf/5)*5, levels= sort(unique(ceiling(tmp$hbf/5)*5)) )

fill.hbf = tmp %>% filter(is.na(hbf)) %>% mutate(hbf = hbf.estimates$mode, hbf.fill=1)
good.hbf =  tmp %>% filter(!is.na(hbf)) %>% mutate(hbf.fill=0)

ops.filt= rbind(fill.hbf, good.hbf) %>% mutate(hbf = as.numeric(hbf))

dir.create(path ='Data/', recursive = T)

save(ops.filt, file='Data/ops.filt.Rdata')
