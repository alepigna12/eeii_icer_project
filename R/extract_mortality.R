library(here)
library(tidyverse)
male = read.csv(here("data", "USA_mltper_1x1.csv"))
female = read.csv(here("data", "USA_fltper_1x1.csv"))
male = male %>%
  filter(Year == 2018) %>%
  mutate(Age = as.numeric(Age)) %>%
  select(Age, qx) %>%
  rename("d_prob_m" = qx) 
female = female %>%
  filter(Year == 2018) %>%
  mutate(Age = as.numeric(Age)) %>%
  select(Age, qx) %>%
  rename("d_prob_f" = qx)
data = merge(male, female, "Age") %>%
  select(Age, d_prob_m, d_prob_f)
write.csv(data, here("data", "all_cause_death_probabilities.csv"))
