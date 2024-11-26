library(tidyverse)
library(magrittr)
library(data.table)

setwd('E:/ALB_CPUE/2024/')

# Read in filtered data - from filter_data.r
load('Data/ops.filt.Rdata')

##############################################################################
setDT(ops.filt)

cluster.dt = ops.filt[,.(alb=sum(alb_n),bet=sum(bet_n),yft=sum(yft_n),swo=sum(swo_n),mls=sum(mls_n), 
                         hbf = mean(hbf)),by=trip_id]
# Min max conversion to make HBF between 0-1
cluster.dt$hbf = (cluster.dt$hbf-min(cluster.dt$hbf))/(max(cluster.dt$hbf)-min(cluster.dt$hbf))
cluster.dt$total.catch = cluster.dt$alb + cluster.dt$bet + cluster.dt$yft + cluster.dt$swo # + cluster.dt$mls
cluster.dt$alb =  cluster.dt$alb/cluster.dt$total.catch
cluster.dt$bet =  cluster.dt$bet/cluster.dt$total.catch
cluster.dt$yft =  cluster.dt$yft/cluster.dt$total.catch
cluster.dt$swo =  cluster.dt$swo/cluster.dt$total.catch
# cluster.dt$mls =  cluster.dt$mls/cluster.dt$total.catch

# deal with trips that had zero catch of any species
cluster.dt$alb =  ifelse(is.na(cluster.dt$alb),0,cluster.dt$alb)
cluster.dt$bet =  ifelse(is.na(cluster.dt$bet),0,cluster.dt$bet)
cluster.dt$yft =  ifelse(is.na(cluster.dt$yft),0,cluster.dt$yft)
cluster.dt$swo =  ifelse(is.na(cluster.dt$swo),0,cluster.dt$swo)

k2 = kmeans(as.matrix(cluster.dt[,.(alb,bet,yft,swo,hbf)]), centers=2, iter.max = 30, nstart = 30)
k3 = kmeans(as.matrix(cluster.dt[,.(alb,bet,yft,swo,hbf)]), centers=3, iter.max = 30, nstart = 30)
k4 = kmeans(as.matrix(cluster.dt[,.(alb,bet,yft,swo,hbf)]), centers=4, iter.max = 30, nstart = 30)
k5 = kmeans(as.matrix(cluster.dt[,.(alb,bet,yft,swo,hbf)]), centers=5, iter.max = 30, nstart = 30)

tmp  = as.data.frame(cbind(k2.h=k2$cluster, k3.h=k3$cluster, k4.h=k4$cluster, k5.h=k5$cluster, trip_id=cluster.dt$trip_id))
ops.filt %<>% left_join(tmp )

save(ops.filt, file='Data/ops.model.Rdata')
