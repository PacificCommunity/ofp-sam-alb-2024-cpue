
# Nicholas Ducharme-Barth
# 12/11/2020
# a) identify number of observations by year & 5x5 cell in the South Pacific
# b) observations by flag by year & 5x5 cell in the South Pacific

#_____________________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(ggplot2)
	library(ggthemes)
	library(maps)
	library(sf)
	library(magrittr)

#_____________________________________________________________________________________________________________________________________
# read in data
	load("E:/ALB_CPUE/2024/Data/ops.trim.dt.RData")
	load("E:/ALB_CPUE/2024/Background_Data/coast.RData")

#_____________________________________________________________________________________________________________________________________
# define plotting quantities
 	ltb.cpue.col = c("royalblue3","deepskyblue1","gold","orange1","indianred1","firebrick2","#AC2020")
 	all.codes = sort(unique(ops.trim.dt$flag_id))
	dwfn.codes = c("CN","JP","KR","TW","US")
	dwfn.names = c("China", "Japan", "Korea", "Chinese Taipei", "United States")
	pict.codes = c("CK","FJ","FM","GU","KI","MH","NC","NU","PF","PG","PW","SB","TO","TV","VU","WS")

#_____________________________________________________________________________________________________________________________________
# make directory
	save.path = "E:/ALB_CPUE/2024/Figures/Spatial_summaries/"
	dir.create(save.path, recursive = TRUE)

#_____________________________________________________________________________________________________________________________________
# spatiotemporal summaries
	coast = st_as_sf(coast)
	t.block = 1

	# number of records per time block
	tmp.dt = copy(ops.trim.dt)
	g = tmp.dt %>% .[latd<0] %>% .[,latd:=floor(latd/5)*5] %>% .[,lond:=floor(lond/5)*5] %>%
				   .[,year.month := year + (as.numeric(month)-1)/12] %>% 
				   .[,time.block := floor(year.month/t.block)*t.block] %>%	
				   .[,.N,by=.(time.block,lond,latd)] %>%
	ggplot() + theme_few() +
	xlab("Longitude") +
	ylab("Latitude") +
	ggtitle("Effort (Total sets) - All fleets") +
	geom_tile(aes(x=lond,y=latd,fill=N)) + 
	facet_wrap(~time.block,drop=FALSE) + 
	scale_fill_gradientn("Total sets",colors=ltb.cpue.col,trans="log10") +
	# scale_fill_viridis_c("Total effort (Millions of hooks)",option="C", label=label.fn) +
	geom_sf(data=coast) + coord_sf(xlim = c(95, 295), ylim = c(-50, 0), expand = FALSE) 
	ggsave("spo.total.sets.1yr.png",plot=g, device = "png", path = save.path,scale = 1, width = 16, height = 9, units = c("in")); rm(list=c("tmp.dt","g")); gc(verbose=FALSE)

	# number of records per time block by flag group
	tmp.dt = copy(ops.trim.dt)
	g = tmp.dt %>% .[latd<0] %>% .[,latd:=floor(latd/5)*5] %>% .[,lond:=floor(lond/5)*5] %>%
				   .[,flag.grp := flag_id] %>% 
		           .[flag_id %in% pict.codes,flag.grp := "PICT"] %>%
		           .[!(flag_id %in% dwfn.codes | flag_id %in% pict.codes),flag.grp := "OTH"] %>% .[,flag.grp := factor(flag.grp,levels=c("CN","JP","KR","TW","US","PICT","OTH"))] %>%
				   .[,year.month := year + (as.numeric(month)-1)/12] %>% 
				   .[,time.block := floor(year.month/10)*10] %>%	
				   .[,.N,by=.(time.block,lond,latd,flag.grp)] %>%
	ggplot() + theme_few() +
	xlab("Longitude") +
	ylab("Latitude") +
	ggtitle("Decadal Effort (Total sets) by flag group") +
	geom_tile(aes(x=lond,y=latd,fill=N)) + 
	facet_grid(flag.grp~time.block,drop=FALSE) + 
	scale_fill_gradientn("Total sets",colors=ltb.cpue.col,trans="log10") +
	# scale_fill_viridis_c("Total effort (Millions of hooks)",option="C", label=label.fn) +
	geom_sf(data=coast) + coord_sf(xlim = c(95, 295), ylim = c(-50, 0), expand = FALSE) 
	ggsave("spo.total.sets.10yr.flag.png",plot=g, device = "png", path = save.path,scale = 1, width = 16, height = 9, units = c("in")); rm(list=c("tmp.dt","g")); gc(verbose=FALSE)

	# number of hooks fished per time block
	tmp.dt = copy(ops.trim.dt)
	g = tmp.dt %>% .[latd<0] %>% .[,latd:=floor(latd/5)*5] %>% .[,lond:=floor(lond/5)*5] %>%
				   .[,year.month := year + (as.numeric(month)-1)/12] %>% 
				   .[,time.block := floor(year.month/t.block)*t.block] %>%	
				   .[,.(hhook = sum(hook)/1000000),by=.(time.block,lond,latd)] %>%
	ggplot() + theme_few() +
	xlab("Longitude") +
	ylab("Latitude") +
	ggtitle("Effort (Millions of hooks fished) - All fleets") +
	geom_tile(aes(x=lond,y=latd,fill=hhook)) + 
	facet_wrap(~time.block,drop=FALSE) + 
	scale_fill_gradientn("Total hooks (Millions)",colors=ltb.cpue.col,trans="log10") +
	# scale_fill_viridis_c("Total effort (Millions of hooks)",option="C", label=label.fn) +
	geom_sf(data=coast) + coord_sf(xlim = c(95, 295), ylim = c(-50, 0), expand = FALSE) 
	ggsave("spo.total.hooks.1yr.png",plot=g, device = "png", path = save.path,scale = 1, width = 16, height = 9, units = c("in")); rm(list=c("tmp.dt","g")); gc(verbose=FALSE)

	# number of hooks fished per time block by flag group
	tmp.dt = copy(ops.trim.dt)
	g = tmp.dt %>% .[latd<0] %>% .[,latd:=floor(latd/5)*5] %>% .[,lond:=floor(lond/5)*5] %>%
				   .[,flag.grp := flag_id] %>% 
		           .[flag_id %in% pict.codes,flag.grp := "PICT"] %>%
		           .[!(flag_id %in% dwfn.codes | flag_id %in% pict.codes),flag.grp := "OTH"] %>% .[,flag.grp := factor(flag.grp,levels=c("CN","JP","KR","TW","US","PICT","OTH"))] %>%
				   .[,year.month := year + (as.numeric(month)-1)/12] %>% 
				   .[,time.block := floor(year.month/10)*10] %>%	
				   .[,.(hhook = sum(hook)/1000000),by=.(time.block,lond,latd,flag.grp)] %>%
	ggplot() + theme_few() +
	xlab("Longitude") +
	ylab("Latitude") +
	ggtitle("Decadal Effort (Millions of hooks fished) by flag group") +
	geom_tile(aes(x=lond,y=latd,fill=hhook)) + 
	facet_grid(flag.grp~time.block,drop=FALSE) + 
	scale_fill_gradientn("Total hooks (Millions)",colors=ltb.cpue.col,trans="log10") +
	# scale_fill_viridis_c("Total effort (Millions of hooks)",option="C", label=label.fn) +
	geom_sf(data=coast) + coord_sf(xlim = c(95, 295), ylim = c(-50, 0), expand = FALSE) 
	ggsave("spo.total.hooks.10yr.flag.png",plot=g, device = "png", path = save.path,scale = 1, width = 16, height = 9, units = c("in")); rm(list=c("tmp.dt","g")); gc(verbose=FALSE)

# alb cpue per time block
	tmp.dt = copy(ops.trim.dt)
	g = tmp.dt %>% .[latd<0] %>% .[,latd:=floor(latd/5)*5] %>% .[,lond:=floor(lond/5)*5] %>%
				   .[,year.month := year + (as.numeric(month)-1)/12] %>% 
				   .[,time.block := floor(year.month/t.block)*t.block] %>%	
				   .[,.(cpue = mean(alb_cpue)),by=.(time.block,lond,latd)] %>%
	ggplot() + theme_few() +
	xlab("Longitude") +
	ylab("Latitude") +
	ggtitle("ALB CPUE - All fleets") +
	geom_tile(aes(x=lond,y=latd,fill=cpue)) + 
	facet_wrap(~time.block,drop=FALSE) + 
	scale_fill_gradientn("Average CPUE (#s per 100 hooks)",colors=ltb.cpue.col,trans="sqrt") +
	# scale_fill_viridis_c("Total effort (Millions of hooks)",option="C", label=label.fn) +
	geom_sf(data=coast) + coord_sf(xlim = c(95, 295), ylim = c(-50, 0), expand = FALSE) 
	ggsave("spo.alb.cpue.1yr.png",plot=g, device = "png", path = save.path,scale = 1, width = 16, height = 9, units = c("in")); rm(list=c("tmp.dt","g")); gc(verbose=FALSE)

	# alb cpue fished per time block by flag group
	tmp.dt = copy(ops.trim.dt)
	g = tmp.dt %>% .[latd<0] %>% .[,latd:=floor(latd/5)*5] %>% .[,lond:=floor(lond/5)*5] %>%
				   .[,flag.grp := flag_id] %>% 
		           .[flag_id %in% pict.codes,flag.grp := "PICT"] %>%
		           .[!(flag_id %in% dwfn.codes | flag_id %in% pict.codes),flag.grp := "OTH"] %>% .[,flag.grp := factor(flag.grp,levels=c("CN","JP","KR","TW","US","PICT","OTH"))] %>%
				   .[,year.month := year + (as.numeric(month)-1)/12] %>% 
				   .[,time.block := floor(year.month/10)*10] %>%	
				   .[,.(cpue = mean(alb_cpue)),by=.(time.block,lond,latd,flag.grp)] %>%
	ggplot() + theme_few() +
	xlab("Longitude") +
	ylab("Latitude") +
	ggtitle("ALB decadal CPUE by flag group") +
	geom_tile(aes(x=lond,y=latd,fill=cpue)) + 
	facet_grid(flag.grp~time.block,drop=FALSE) + 
	scale_fill_gradientn("Average CPUE (#s per 100 hooks)",colors=ltb.cpue.col,trans="sqrt") +
	# scale_fill_viridis_c("Total effort (Millions of hooks)",option="C", label=label.fn) +
	geom_sf(data=coast) + coord_sf(xlim = c(95, 295), ylim = c(-50, 0), expand = FALSE) 
	ggsave("spo.alb.cpue.10yr.flag.png",plot=g, device = "png", path = save.path,scale = 1, width = 16, height = 9, units = c("in")); rm(list=c("tmp.dt","g")); gc(verbose=FALSE)

	## plot sets per year-quarter by management area
	## identify regions of interest
	library(tidyverse)
	setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
	load("../Background_Data/alb_regions_poly_2024_shp.RData")
	regions.shp_sf = st_as_sf(regions.shp)
	num.regions <- 2
	ops.trim.dt$Region = NA
	for(i in 1:num.regions){
	  coords = regions.shp@polygons[[i]]@Polygons[[1]]@coords
	  ops.trim.dt$Region = ifelse(sp::point.in.polygon(ops.trim.dt$latd,
	                                                ops.trim.dt$lond, coords[,2],
	                                                coords[,1]) %in%
	                             c(1,2),i,ops.trim.dt$Region)
	}

	ops.trim.dt$decade <- 10 * floor(ops.trim.dt$year / 10)
	effort <- ops.trim.dt %>% group_by(latd, lond, decade) %>%
	  summarise(effort = sum(hhook)) %>% ungroup()

eff <- 	effort %>%
	  ggplot() + theme_few() +
	  labs(x = "Longitude", y = "Latitude",
	       title = "Total longline effort (log hundred hooks)") +
	  geom_tile(aes(x=lond,y=latd,fill=log(effort))) + 
    geom_sf(data=regions.shp_sf, fill = NA) +
	  facet_wrap(~decade,drop=FALSE) +
	  scale_fill_gradientn("Effort",colors=ltb.cpue.col) +
    theme(text = element_text(size = 16)) +
	  geom_sf(data=coast) + coord_sf(xlim = c(110, 290), ylim = c(-55, 5), expand = FALSE) 
eff
ggsave("spo.alb.decadal.effort.png",plot=eff, device = "png", path = save.path,scale = 1, width = 16, height = 9, units = c("in")); rm(list=c("tmp.dt","g")); gc(verbose=FALSE)
	
	catch.df <- ops.trim.dt %>% group_by(latd, lond, decade) %>%
	  summarise(catch = sum(alb_n)) %>% ungroup()
	
	cat <- 	
	  catch.df %>%
	  ggplot() + theme_few() +
	  xlab("Longitude") +
	  ylab("Latitude") +
	  ggtitle("Total longline catch (log scale)") +
	  geom_tile(aes(x=lond,y=latd,fill=log(catch))) + 
	  geom_sf(data=regions.shp_sf, fill = NA) +
	  facet_wrap(~decade,drop=FALSE) + 
	  scale_fill_gradientn("# Alb",colors=ltb.cpue.col) +
	  theme(text = element_text(size = 16)) +
	  geom_sf(data=coast) + coord_sf(xlim = c(110, 290), ylim = c(-55, 5), expand = FALSE) 
	cat
	ggsave("spo.alb.decadal.catch.png",plot=cat, device = "png", path = save.path,scale = 1, width = 16, height = 9, units = c("in")); rm(list=c("tmp.dt","g")); gc(verbose=FALSE)
	
	cpue.df <- ops.trim.dt %>% group_by(latd, lond, decade) %>%
	  summarise(cpue = mean(alb_cpue)) %>% ungroup()
	
	cpue.plot <- 	
	  cpue.df %>%
	  ggplot() + theme_few() +
	  xlab("Longitude") +
	  ylab("Latitude") +
	  ggtitle("Average CPUE (# alb per 100 hooks)") +
	  geom_tile(aes(x=lond,y=latd,fill=cpue)) + 
	  geom_sf(data=regions.shp_sf, fill = NA) +
	  facet_wrap(~decade,drop=FALSE) + 
	  scale_fill_gradientn("CPUE",colors=ltb.cpue.col,trans="sqrt") +
	  theme(text = element_text(size = 16)) +
	  geom_sf(data=coast) + coord_sf(xlim = c(110, 290), ylim = c(-55, 5), expand = FALSE) 
	cpue.plot
	ggsave("spo.alb.decadal.cpue.png",plot=cpue.plot, device = "png", path = save.path,scale = 1, width = 16, height = 9, units = c("in")); rm(list=c("tmp.dt","g")); gc(verbose=FALSE)
	
	## get how many NA for HBF over time
	ops.full.dt = fread("E:/tuna_dbs/LOG_DBS/DBF/l_opn_2024-02-01.csv")# this dataset was created on 20/01/2021, will be mostly complete for 2019 data but should get updated
	ops.full.dt$year = as.numeric(substr(ops.full.dt$logdate,7,10))
	ops.full.dt %<>% mutate(decade = 10 * floor(year / 10))

		HBF.df <- ops.full.dt %>%
		  mutate(fleet = ifelse(flag_id %in% pict.codes, "PICT",
		                        ifelse(flag_id %in% c("JP", "TW", "KR", "CN", "US"),
		                               flag_id, "OTH"))) %>%
	  mutate(hbf.pres = ifelse(hk_bt_flt %in% c(0, -1, "**"), 0, 1),
	         total = 1) %>%
	  group_by(fleet, decade) %>%
	  summarise(freq = sum(hbf.pres), total = sum(total)) %>% ungroup() %>%
	  mutate(prop = round(freq / total, 2))
	
		hbf.plot <- 	
		  HBF.df %>%
		  ggplot() + theme_few() +
		  labs(x = "Decade", y = "Proportion Missing HBF", shape = "Flag",
		       col = "Flag") +
		  ggtitle("Missing HBF - by flag") +
		  geom_point(aes(x=decade,y=prop,shape=fleet, col = fleet), size = 6) +
		  scale_shape_manual(values=c(15:19, 15, 17))
		hbf.plot
		ggsave("spo.alb.decadal.hbf.png",plot=hbf.plot, device = "png", path = save.path,scale = 1, width = 16, height = 9, units = c("in")); rm(list=c("tmp.dt","g")); gc(verbose=FALSE)
		