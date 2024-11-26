
# Nicholas Ducharme-Barth
# 12/11/2020
# Data prep for operational longline data
# In support of 2021 albacore assessment following data cleaning used in 2020 BET & YFT

# Revised 7/05/2021 to include vessel data and oceanography
# T. Vidal

term.yr <- 2022 # terminal year of assessment

# load packages
	library(data.table)

# read in full operational data
	ops.full.dt = fread("E:/tuna_dbs/LOG_DBS/DBF/l_opn_2024-02-01.csv")# this dataset was created on 20/01/2021, will be mostly complete for 2019 data but should get updated

# Make missing vessel proper NA
	ops.full.dt$vessel_id = ifelse(ops.full.dt$vessel_id=='',NA,ops.full.dt$vessel_id)
	ops.rec = nrow(ops.full.dt)

# South Pacific only
  ops.full.dt =  ops.full.dt[latd<=0] 	
# add unique cell id
	ops.full.dt$cell.1x1 = paste0(floor(ops.full.dt$lond),"x",floor(ops.full.dt$latd))
	ops.full.dt$cell.5x5 = paste0(floor(ops.full.dt$lond/5)*5,"x",floor(ops.full.dt$latd/5)*5)
	ops.full.dt$cell.10x10 = paste0(floor(ops.full.dt$lond/10)*10,"x",floor(ops.full.dt$latd/10)*10)
	ops.full.dt$cell.15x15 = paste0(floor(ops.full.dt$lond/15)*15,"x",floor(ops.full.dt$latd/15)*15)

# add other columns
	ops.full.dt$flag.fleet = paste0(ops.full.dt$flag_id,".",ops.full.dt$fleet_id)
	ops.full.dt$year = as.numeric(substr(ops.full.dt$logdate,7,10))
	ops.full.dt$month = substr(ops.full.dt$logdate,4,5)
	ops.full.dt$quarter = c(1,1,1,2,2,2,3,3,3,4,4,4)[as.numeric(ops.full.dt$month)]
	ops.full.dt$yr.qtr = paste0(ops.full.dt$year,"x",ops.full.dt$quarter)
	ops.full.dt$yr.month = paste0(ops.full.dt$year,"_",ops.full.dt$month)
	ops.full.dt$yr.month.cell.1x1 = paste0(ops.full.dt$yr.month,"_",ops.full.dt$cell.1x1)
	ops.full.dt$yr.month.cell.5x5 = paste0(ops.full.dt$yr.month,"_",ops.full.dt$cell.5x5)
	ops.full.dt$yr.5 = floor(ops.full.dt$year/5)*5
	ops.full.dt$yr.10 = floor(ops.full.dt$year/10)*10 
	ops.full.dt = ops.full.dt[!is.na(year)]
	ops.rec = c(ops.rec,nrow(ops.full.dt))
	ops.full.dt = ops.full.dt[year<=term.yr]
	ops.rec = c(ops.rec,nrow(ops.full.dt))

# CPUE related columns
	ops.full.dt$hhook = ops.full.dt$hook/100
	ops.full.dt$alb_cpue = ops.full.dt$alb_n/ops.full.dt$hhook
	ops.full.dt$bet_cpue = ops.full.dt$bet_n/ops.full.dt$hhook
	ops.full.dt$yft_cpue = ops.full.dt$yft_n/ops.full.dt$hhook
	ops.full.dt$swo_cpue = ops.full.dt$swo_n/ops.full.dt$hhook
	ops.full.dt$mls_cpue = ops.full.dt$mls_n/ops.full.dt$hhook

# add column for year-quarter aka ts
	ts.df = expand.grid(quarter=1:4,year=1952:(term.yr + 1))
	ts.df$ts = 1:nrow(ts.df)
	ts.df$yr.qtr = paste0(ts.df$year,"x",ts.df$quarter) 
	ops.full.dt$ts = ts.df$ts[match(ops.full.dt$yr.qtr,ts.df$yr.qtr)]

# remove outliers & limit geographical scope...
	# geographical scope
		load("E:/ALB_CPUE/2024/Background_Data/po.all.RData")
		load("E:/ALB_CPUE/2021/Background_Data/coast.RData")
		po.coast = rgeos::gIntersection(coast,po.all)

		tmp = data.frame(cell.1x1 =  unique(ops.full.dt$cell.1x1),stringsAsFactors = FALSE)
		tmp$lond = sapply(tmp$cell.1x1,function(x)as.numeric(strsplit(x,"x")[[1]][1])+0.5)
		tmp$latd = sapply(tmp$cell.1x1,function(x)as.numeric(strsplit(x,"x")[[1]][2])+0.5)
		tmp = as.data.table(tmp)
		coords = as.matrix(tmp[,.(lond,latd)])
		sp = sp::SpatialPoints(coords,proj4string=sp::CRS(sp::proj4string(po.all)))
		tmp$po.over.index = sp::over(sp,po.all)
		tmp$po.over.index = ifelse(is.na(tmp$po.over.index),0,tmp$po.over.index)
		tmp$coast.over.index = sp::over(sp,po.coast)
		tmp$coast.over.index = ifelse(is.na(tmp$coast.over.index),1,0)
		
		tmp$coast.over.index[which(tmp$lond>=150 & tmp$latd>=-23 & tmp$coast.over.index == 0 & tmp$latd<35 & tmp$lond<230)] = 1
		tmp$keep.index = ifelse(tmp$coast.over.index+tmp$po.over.index == 2,1,0)
		tmp = tmp[keep.index == 1]

		ops.full.dt  = ops.full.dt[cell.1x1 %in% tmp$cell.1x1]
		ops.rec = c(ops.rec,nrow(ops.full.dt))

	# hbf
	# convert bad values for hk_bt_flt
		ops.full.dt$hk_bt_flt[which(ops.full.dt$hk_bt_flt %in% c(-1,"**"))] = NA
		ops.full.dt$hk_bt_flt = as.numeric(ops.full.dt$hk_bt_flt)
		ops.full.dt  = ops.full.dt[hk_bt_flt <= 50]
		ops.rec = c(ops.rec,nrow(ops.full.dt))
		
	# hooks fished
		ops.full.dt  = ops.full.dt[hook <= 5000 & hook >150]
		ops.rec = c(ops.rec,nrow(ops.full.dt))

	# more fish caught than number of hooks fished...
		ops.full.dt = ops.full.dt[-which(rowSums(ops.full.dt[,.(alb_n,bet_n,yft_n,swo_n,mls_n)])-ops.full.dt$hook>0)]
		ops.rec = c(ops.rec,nrow(ops.full.dt))

	# look at vessel_id & anonymize vessel_id
		ops.full.dt$vessel_id = as.character(as.numeric(as.factor(ops.full.dt$vessel_id)))
		ops.full.dt$vessel_id[which(ops.full.dt$vessel_id == "1")] = NA
		vessel.dt = ops.full.dt[,.(first.ts=min(ts),last.ts=max(ts),unique.ts=uniqueN(ts),n.sets.fished=.N),by=vessel_id][order(first.ts,last.ts,unique.ts,n.sets.fished)]
		# remove vessels (with an id) that have less than 10 unique quarters fished or less than 30 sets
		ops.full.dt = ops.full.dt[vessel_id %in% vessel.dt[(n.sets.fished>=30 & unique.ts>=10)|(n.sets.fished>=30 & first.ts>=259)]$vessel_id]
		ops.rec = c(ops.rec,nrow(ops.full.dt))

	# remove data from ES, PA, BZ
		ops.full.dt  = ops.full.dt[!(flag_id %in% c("ES",'PA','BZ'))]
		ops.rec = c(ops.rec,nrow(ops.full.dt))

# assign trip designations
	# add column for vessel_yr.month
	# add column for flag.fleet_yr.month_cell.5x5
	ops.full.dt$trip_id = ifelse(is.na(ops.full.dt$vessel_id) ,paste0(ops.full.dt$flag.fleet,"_",ops.full.dt$yr.month.cell.5x5),paste0(ops.full.dt$vessel_id,"_",ops.full.dt$yr.month))

# clustering analysis (by vessel_yr.month or flag.fleet_yr.month_cell.5x5)
	cluster.dt = ops.full.dt[,.(alb=sum(alb_n),bet=sum(bet_n),yft=sum(yft_n),swo=sum(swo_n),mls=sum(mls_n)),by=trip_id]
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
		# cluster.dt$mls =  ifelse(is.na(cluster.dt$mls),0,cluster.dt$mls)

		set.seed(123)

		k2 = kmeans(as.matrix(cluster.dt[,.(alb,bet,yft,swo)]), centers=2, iter.max = 30, nstart = 30)
		k3 = kmeans(as.matrix(cluster.dt[,.(alb,bet,yft,swo)]), centers=3, iter.max = 30, nstart = 30)
		k4 = kmeans(as.matrix(cluster.dt[,.(alb,bet,yft,swo)]), centers=4, iter.max = 30, nstart = 30)
		k5 = kmeans(as.matrix(cluster.dt[,.(alb,bet,yft,swo)]), centers=5, iter.max = 30, nstart = 30)
		k6 = kmeans(as.matrix(cluster.dt[,.(alb,bet,yft,swo)]), centers=6, iter.max = 30, nstart = 30)
		k7 = kmeans(as.matrix(cluster.dt[,.(alb,bet,yft,swo)]), centers=7, iter.max = 30, nstart = 30)

		 elbow.vec = c(k2$tot.withinss,k3$tot.withinss,k4$tot.withinss,k5$tot.withinss,k6$tot.withinss,k7$tot.withinss)
		 #par(mar=c(4,5,1,1))
		 #png(filename = "Figures/kmeans.elbow.png",width = 6, height = 6, units = "in", res = 300)
		 plot(2:7,elbow.vec/10000,xlab="Number of clusters (k)",ylab="Within cluster sum-of-squares (x 10,000)",cex=1.5,cex.lab=1.5,cex.axis=1.5,las=1,type="b",pch=16)
		 #dev.off()

	# add back in
		ops.full.dt$k2 = k2$cluster[match(ops.full.dt$trip_id,cluster.dt$trip_id)]
		ops.full.dt$k3 = k3$cluster[match(ops.full.dt$trip_id,cluster.dt$trip_id)]
		ops.full.dt$k4 = k4$cluster[match(ops.full.dt$trip_id,cluster.dt$trip_id)]
		#ops.full.dt$k5 = k5$cluster[match(ops.full.dt$trip_id,cluster.dt$trip_id)]
		#ops.full.dt$k6 = k6$cluster[match(ops.full.dt$trip_id,cluster.dt$trip_id)]
		#ops.full.dt$k7 = k7$cluster[match(ops.full.dt$trip_id,cluster.dt$trip_id)]

# rename
	ops.trim.dt = ops.full.dt

# add record_id
	ops.trim.dt$record_id = 1:nrow(ops.trim.dt)

# save "cleaned" data
	save(ops.trim.dt,file = "E:/ALB_CPUE/2024/Data/ops.trim.dt.RData")
	