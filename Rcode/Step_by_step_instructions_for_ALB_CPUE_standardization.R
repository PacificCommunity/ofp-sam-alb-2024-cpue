## Thom Teears
## 01/02/2024

##################### Steps for ALB CPUE standardization using VAST
## Run prep.ops.data.r to prepare the raw data
## Run plot.spatial.summaries.r to plot summaries of data
## run alb_filter_data.r to filter the data for albacore
## run target_cluster.r to estimate targeting cluster
## Run plot.spatial.summaries.r to plot summaries of data
## run alb_cpue_run.r to perform cpue analysis in VAST
## **Note: there is some data preparation (sub-sampling) that happens
##        in alb_cpue_run.r
## get results using get_indices.r

##################### Steps for ALB CPUE standardization using sdmTMB
## Run fit.sdmTMB.spawn to fit and predict the spawn indices
## Run fit.sdmTMB.global to fit and predict the global indices
## Run calc.prop.knots.sampled.r to get spatial coverage plot
## Run plot_mesh