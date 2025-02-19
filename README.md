# South Pacific ALB 2024 CPUE Analysis

## Two analytical steps

1. Fit model
2. Calculate predictions

The predictions are calculated with the `get_index()` function, weighted by
area. Some 5x5 rectangles are larger than others, depending on the latitude. The
prediction uses a regular grid called `new.data`.

## Spatial resolution

The observed data and model fit are based on 1x1 degree rectangles.

Model prediction is based on 5x5 degrees, used for calculating the weighting of
length data and for calculating the relative biomass indices.

## Overview of scripts

### VAST analysis

For the ALB 2024 CPUE analysis, Thom first reran the 2021 VAST analysis that
Tiffany had done. As documented in
`Step_by_step_instructions_for_ALB_CPUE_standardization.R`, these are the
scripts involved in the VAST analysis:

- Run `prep.ops.data.r` to prepare the raw data
- Run `plot.spatial.summaries.r` to plot summaries of data
- Run `alb_filter_data.r` to filter the data for albacore
- Run `target_cluster.r` to estimate targeting cluster (target species)
- Run `alb_cpue_run.r` to perform cpue analysis in VAST

*Note:* there is some data preparation (sub-sampling) that happens in
`alb_cpue_run.r`

- Get results using `get_indices.r`

### sdmTMB analysis

Next, Thom developed scripts to conduct `sdmTMB` analysis of ALB 2024 CPUE,
loosely based on YFT 2023 CPUE analysis.

One difference between `VAST` and `sdmTMB` is that `sdmTMB` uses a mesh that
goes well beyond the geographic extent of the data.

These are the scripts to run the `sdmTMB` analysis:

- Run `fit.sdmTMB.spawn.R `to fit and predict the spawning region indices
- Run `fit.sdmTMB.global.R` to fit and predict the global indices
- Run `calc.prop.knots.sampled.r` to get spatial coverage plot
- Run `plot_mesh.r`

Encounter probability is used when a delta model is fitted and the data includes
zero-catch observations.

Hooks between floats: check with Thom.

### Other scripts

Script                      | Purpose
--------------------------- | -------------------------------------------------
`alb_add_regions.R`         | VAST code needed somewhere
`alb_cpue_prep.R`           | VAST code needed somewhere
`compare.catch.quarterly.R` | Plot catch and effort time series by region
`fit.sdmTMB.vessel.R`       | sdmTMB using vessel ID as fixed effect, not used
`fit.sdmTMB.vessel.RE.R`    | sdmTMB using vessel ID as random effect, not used
`select_core_vessels.R`     | Experiment with removing some 5x5, not used
`get_ENSO.R`                | Get ENSO data from web (might want to keep month)
`get_indices.R`             | Compare VAST and sdmTMB results
