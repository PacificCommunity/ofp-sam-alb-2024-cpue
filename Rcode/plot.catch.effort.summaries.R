library(tidyverse)
library(magrittr)

# Remove all previous objects
rm(list=ls())

setwd('E:/ALB_CPUE/2024')

# Read in filtered data - from filter_data.r
load('Data/ops.model.Rdata')

dir.create("Figures/Catch_effort_summaries")

###############*
## Formatting
###############*
gg.theme <- theme_bw() +
  theme(axis.line = element_line(color="black"),
        panel.grid.minor = element_blank(),
        text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.title = element_text(size=14),
        strip.background =element_rect(fill='white'),
        strip.text = element_text(size=14),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"))

library(ggplot2)

ops.filt %>% mutate(Region = ifelse(Region == 1, "WCPFC-CA",
                                    ifelse(Region == 2, "EPO", NA))) %>%
  mutate(Region = factor(Region, levels = c("WCPFC-CA", "EPO"))) %>%
  filter(!is.na(Region),
         flag_id %in% c("CN", "FJ", "JP", "KR", "PF", "TW", "US","VU",
                        "WS")) %>%
  group_by(flag_id, year, Region) %>%
  summarise(alb = sum(alb)) %>% ungroup() %>%
  ggplot() +
  facet_grid(~Region) +
  geom_col(aes(x = year, y = alb, fill = flag_id)) +
  scale_x_continuous(breaks = seq(1950, 2030, 10)) +
  scale_fill_manual(values = c("indianred1", "deepskyblue1", "palegreen3",
                               "orange1", "royalblue3", "wheat2", "#AC2020",
                               "darkgrey", "slateblue")) +
  labs(x = "Year",y = "Catch (in Numbers)", fill = "Flag") +
  gg.theme

ggsave(paste0('Figures/Catch_effort_summaries/catch_by_year_and_flag.png'),
       width=7, height=5, units='in', dpi=300)

effort.df <- ops.filt %>% mutate(effort = 1) %>%
  mutate(Region = ifelse(Region == 1, "WCPFC-CA",
                         ifelse(Region == 2, "EPO", NA))) %>%
  mutate(Region = factor(Region, levels = c("WCPFC-CA", "EPO"))) %>%
  group_by(flag_id, Region, year) %>% summarise(effort = sum(effort)) %>%
  ungroup()

library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
effort.df %>% filter(!is.na(Region)) %>%
ggplot() +
  facet_wrap(~Region, ncol = 1) +
  geom_tile(aes(x = year, y = flag_id, fill = effort)) +
  scale_x_continuous(breaks = seq(1950, 2030, 10)) +
  scale_fill_gradientn(colours = myPalette(100), limits=c(1, 50000)) +
  labs(x = "Year",y = "Flag", fill = "Effort (# sets)") +
  gg.theme

ggsave(paste0('Figures/Catch_effort_summaries/effort_by_year_and_flag.png'),
       width=7, height=8, units='in', dpi=300)

prop.catch <- ops.filt %>% group_by(year, Region) %>%
  summarise(alb = sum(alb_n), bet = sum(bet_n), yft = sum(yft_n),
            swo = sum(swo_n), total.catch = sum(total.catch)) %>%
  ungroup() %>% filter(!is.na(Region)) %>%
  pivot_longer(cols = c("alb", "bet", "yft", "swo"), names_to = "sps",
               values_to = "catch") %>% mutate(prop = catch / total.catch) %>%
  mutate(Region = ifelse(Region == 1, "WCPFC-CA", "EPO")) %>%
  mutate(Region = factor(Region, levels = c("WCPFC-CA", "EPO")))

ggplot(prop.catch, aes(x = year, y = prop, fill = as.factor(sps))) +
  facet_wrap(~Region, ncol = 1) +
  geom_bar(stat = "identity", position = "fill") +
  scale_x_continuous(breaks = seq(1950, 2030, 5)) +
  scale_fill_manual(values = c("palegreen3", "#AC2020",
                               "slateblue", "wheat2")) +
  theme_void() + gg.theme +
  labs(x = "Year", y = "Proportion of Catch", fill = "Species")

ggsave(paste0('Figures//Catch_effort_summaries/proportion_catch_by_year_and_species.png'),
       width=7, height=5, units='in', dpi=300)

prop.hbf <- ops.filt %>% filter(!is.na(Region)) %>% mutate(freq = 1) %>%
  group_by(year, hbf, Region) %>% summarise(freq = sum(freq)) %>%
  ungroup()  %>% left_join(ops.filt %>% filter(!is.na(Region)) %>%
                             mutate(freq = 1) %>% group_by(year, Region) %>%
                             summarise(total = sum(freq)) %>% ungroup()) %>%
  mutate(Region = ifelse(Region == 1, "WCPFC-CA", "EPO")) %>%
  mutate(Region = factor(Region, levels = c("WCPFC-CA", "EPO")),
         prop = freq / total)

ggplot(prop.hbf, aes(x = year, y = prop, fill = as.factor(hbf))) +
  facet_wrap(~Region, ncol = 1) +
  geom_bar(stat = "identity", position = "fill") +
  scale_x_continuous(breaks = seq(1950, 2030, 5)) +
  scale_fill_manual(values = c("indianred1", "deepskyblue1", "palegreen3",
                               "orange1", "royalblue3", "wheat2", "#AC2020",
                               "darkgrey", "slateblue", "gold")) +
  theme_void() + gg.theme +
  labs(x = "Year", y = "Proportion of Catch", fill = "Species")

ggsave(paste0('Figures//Catch_effort_summaries/proportion_hbf_by_year_and_species.png'),
       width=7, height=5, units='in', dpi=300)


catch.df <- ops.filt %>%
  mutate(Region = ifelse(Region == 1, "WCPFC-CA",
                         ifelse(Region == 2, "EPO", NA))) %>%
  mutate(Region = factor(Region, levels = c("WCPFC-CA", "EPO"))) %>%
  group_by(Region, year) %>% summarise(catch = sum(alb)) %>%
  ungroup()

catch.df %>% filter(!is.na(Region)) %>%
  ggplot() +
  facet_wrap(~Region, ncol = 1) +
  geom_col(aes(x = year, y = catch), fill = "darkgreen") +
  scale_x_continuous(breaks = seq(1950, 2030, 10)) +
  labs(x = "Year",y = "Catch (Numbers)") +
  gg.theme

ggsave(paste0('Figures/Catch_effort_summaries/catch_by_year.png'),
       width=7, height=8, units='in', dpi=300)

## plot nominal cpue by flag
## add area to data
library(raster)
r <- raster(ymn = -51, ymx = 1, res = 1)
a <- area(r)
lat <- yFromRow(r, 1:nrow(r))
area <- a[,1]
plot(lat, area, type = "l"); points(lat, area, pch = 16)
area.df <- data.frame(latd = lat,
                      Area_km2 = area)

ops.filt %<>% left_join(area.df)
detach('package:raster')

nom.cpue <- ops.filt %>% filter(!is.na(Region)) %>%
  group_by(year, Region, latd, lond, flag_id) %>%
  summarise(CPUE = mean(alb_cpue),
            area.weight = median(Area_km2)) %>% ungroup() %>%
  mutate(CPUE.weighted = CPUE * area.weight) %>%
  group_by(year, Region, flag_id) %>%
  summarise(Index = sum(CPUE.weighted)) %>% ungroup()

## check correlations
flags <- c("PF", "FJ", "WS", "KR", "TW", "CN", "JP", "US", "VU")
regs <- unique(nom.cpue$Region)

nom.cpue %<>% mutate(time.frame = ifelse(year < 2000, "early", "late"))
time.frames <- unique(nom.cpue$time.frame)

rm(results.df)
for(l in 1:length(time.frames)){
  for(k in 1:length(regs)){
    cors.df <- data.frame(matrix(NA, ncol = length(flags), nrow = length(flags)))
    colnames(cors.df) <- flags
    rownames(cors.df) <- flags
    
    for(i in 1:length(flags)){
      for(j in 1:length(flags)){
        ind1 <- nom.cpue %>% filter(flag_id == flags[i], Region == regs[k],
                                    time.frame == time.frames[l])
        ind2 <- nom.cpue %>% filter(flag_id == flags[j], Region == regs[k],
                                    time.frame == time.frames[l])
        ind1 %<>% filter(year %in% unique(ind2$year))
        ind2 %<>% filter(year %in% unique(ind1$year))
        if(nrow(ind1) < 5){next}
        cors.df[i,j] <- round(cor(ind1$Index, ind2$Index), 2)
        rm(ind1, ind2)
      }}
    cors.df %<>% mutate(Flag1 = rownames(cors.df))
    cors.df$Region <- regs[k]
    cors.df$time.frame <- time.frames[l]
    if(!exists("results.df")){
      results.df <- cors.df
    } else {
      results.df %<>% rbind(cors.df)
    }
    rm(cors.df)
  }
}

write.csv(results.df, file = "Results/nominal.correlations.csv", row.names = F)

flag.cols <- c("indianred1", "deepskyblue1", "palegreen3",
               "orange1", "royalblue3", "#AC2020",
               "darkgrey", "slateblue", "wheat3")
## group non-continuous series
groupFUN <- function(vec){
rr <- rle(vec - seq_along(vec))
rr$values <- seq_along(rr$values)
s <- split(vec, inverse.rle(rr))

for(i in 1:length(s)){
  df <- data.frame(year = s[i], grp = i)
  colnames(df) <- c("year", "grp")
  if(!exists("final.df")){
    final.df <- df
  } else {
    final.df <- rbind(final.df, df)
  }
  rm(df)
}
return(final.df)
}

for(i in 1:length(flags)){
  vec <- nom.cpue %>% filter(flag_id == flags[i]) %>% select(year) %>%
    distinct() %>% unlist() %>% unname()
  res <- groupFUN(vec)
  rm(vec)
  res$flag_id <- flags[i]
  res %<>% mutate(time.frame = ifelse(year < 2000, "early", "late"))
  if(!exists("grp.df")){
    grp.df <- res
  } else {
    res %<>% mutate(grp = grp + max(grp.df$grp))
    grp.df %<>% rbind(res)
  }
  rm(res)
}

## remove some indices to clean up plot
nom.cpue %<>% mutate(id = paste0(Region, flag_id)) %>%
filter(id != "2FJ", id != "2US", id != "2PF") %>% select(-id)
nom.cpue %<>% mutate(id = paste0(time.frame, flag_id)) %>%
  filter(id %in% c("earlyJP", "earlyTW", "earlyKR") | time.frame == "late") %>%
  select(-id)

nom.cpue %<>% left_join(grp.df)

## plot nominal cpue by flag
nom.cpue %>% filter(flag_id %in% flags, time.frame == "early") %>%
  mutate(Region = ifelse(Region == 1, "WCPFC-CA", "EPO")) %>%
  mutate(Region = factor(Region, levels = c("WCPFC-CA", "EPO"))) %>%
  group_by(Region, flag_id) %>%
  mutate(Index = Index / mean(Index)) %>% ungroup() %>%
  ggplot() +
  geom_smooth(aes(x = year, y = Index, col = as.factor(flag_id), group = grp),
            size = 1, se = F) +
  scale_color_manual(values = flag.cols) +
  scale_x_continuous(breaks = seq(1950, 2030, 10)) +
  facet_wrap(~Region, ncol = 1, scales = "free") +
  labs(x = "Year", y = "Mean-centered Nominal CPUE", col = "Flag") +
  gg.theme

ggsave(paste0('Figures/Catch_effort_summaries/nominal.cpue.by.flag.early.png'),
       width=10, height=8, units='in', dpi=300)

nom.cpue %>% filter(flag_id %in% c("CN", "JP", "VU", "KR", "TW", "FJ", "WS",
                                   "US", "PF"), time.frame == "late") %>%
  mutate(Region = ifelse(Region == 1, "WCPFC-CA", "EPO")) %>%
  mutate(Region = factor(Region, levels = c("WCPFC-CA", "EPO"))) %>%
  group_by(Region, flag_id) %>%
  mutate(Index = Index / mean(Index)) %>% ungroup() %>%
  ggplot() +
  geom_smooth(aes(x = year, y = Index, col = as.factor(flag_id), group = grp),
            size = 1, se = F) +
  scale_color_manual(values = flag.cols) +
  scale_x_continuous(breaks = seq(1950, 2030, 10)) +
  facet_wrap(~Region, ncol = 1, scales = "free") +
  labs(x = "Year", y = "Mean-centered Nominal CPUE", col = "Flag") +
  gg.theme

ggsave(paste0('Figures/Catch_effort_summaries/nominal.cpue.by.flag.late.png'),
       width=10, height=8, units='in', dpi=300)
