## extract ENSO data from NOAA

library(dplyr)
library(magrittr)
library(tidyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

fileUrl <- "https://psl.noaa.gov/data/correlation/nina4.data"
download.file(fileUrl, "../Data/ENSO.txt")

ENSO <- read.delim2("../Data/ENSO.txt", sep = "", skip = 1)
colnames(ENSO) <- c("year", paste0("month", seq(1, 12, 1)))
ENSO %<>% mutate(year = as.numeric(year)) %>%
  filter(year >= 1950, year <= as.numeric(format(Sys.time(), "%Y"))) %>%
  pivot_longer(cols = starts_with("month"), names_to = "month",
                       values_to = "ENSO") %>%
  mutate(month = as.numeric(substr(month, 6 , 7)), ENSO = as.numeric(ENSO)) %>%
  filter(ENSO > 0) %>%
  mutate(Qtr = ifelse(month %in% c(1:3), 1,
                      ifelse(month %in% c(4:6), 2,
                             ifelse(month %in% c(7:9), 3,
                                    ifelse(month %in% c(10:12), 4, NA))))) %>%
  mutate(YrQtr = year + (Qtr / 4) - 0.25) %>% select(-Qtr)

write.csv(ENSO, file = "../Data/ENSO.csv", row.names = F)

