
## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(brms)
library(cowplot)

## LOAD DATASET ----
gall_selection.df <- read_csv("gall_selection_data.csv") %>%
  mutate(Plant_Position = as.factor(Plant_Position),
            Gall_Number = as.factor(Gall_Number),
            Treatment.focus = as.factor(Treatment.focus),
            Treatment = as.factor(Treatment),
            Location = as.factor(Location),
            phenology = as.factor(phenology),
            Genotype = as.factor(Genotype),
            Gall_Letter = as.factor(Gall_Letter)) %>%
  unite(Gall_ID, Gall_Number, Gall_Letter, remove = FALSE) %>%
  
  # subset data for analysis
  filter(phenology == "early", Location == "tree")

## EXPLORATORY PLOTS ----
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

# Gall Density
gall_selection.df %>%
  group_by(Genotype, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots), mean) %>%
  ggplot(., aes(x = Genotype, y = Density_per_100_shoots)) + 
  geom_jitter(alpha = 0.5, width = 0.05) + 
  stat_summary(fun.data = "mean_cl_boot", color = "red", size = 1.5) 

# Number of Individuals per gall
gall_selection.df %>%
  group_by(Genotype, Plant_Position, Gall_Number) %>%
  summarise_at(vars(gall_individuals), mean) %>%
  group_by(Genotype, Plant_Position) %>%
  summarise_at(vars(gall_individuals), funs(mean, n())) %>%
  ggplot(., aes(x = Genotype, y = mean)) + 
  geom_jitter(alpha = 0.5, width = 0.05, aes(size = n)) + 
  stat_summary(fun.data = "mean_cl_boot", color = "red", size = 1.5)

# Gall Size
gall_selection.df %>%
  group_by(Genotype, Plant_Position, Gall_Number) %>%
  summarise_at(vars(Gall_Height_mm), mean) %>%
  group_by(Genotype, Plant_Position) %>%
  summarise_at(vars(Gall_Height_mm), funs(mean, n())) %>%
  ggplot(., aes(x = Genotype, y = mean)) + 
  geom_jitter(alpha = 0.5, width = 0.05, aes(size = n)) + 
  stat_summary(fun.data = "mean_cl_boot", color = "red", size = 1.5)

# Ectoparasitism rates and estimates of how effective the exclusion treatment was
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(ectos, total), sum) %>%
  ggplot(., aes(x = Genotype, y = ectos/total, group = Treatment.focus, color = Treatment.focus)) + 
  geom_jitter(alpha = 0.5, width = 0.05, aes(size = total)) + 
  stat_summary(fun.data = "mean_cl_boot", size = 1.5)

# Platygaster parasitism rates and estimates of intraguild predation 
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(platy, total), sum) %>%
  ggplot(., aes(x = Genotype, y = platy/total, group = Treatment.focus, color = Treatment.focus)) + 
  geom_jitter(alpha = 0.5, width = 0.05, aes(size = total)) + 
  stat_summary(fun.data = "mean_cl_boot", size = 1.5)
