## LOAD REQUIRED LIBRARIES & SET OPTIONS ----
library(tidyverse)
library(cowplot) # pretty default ggplots
library(broom) # for tidying multiple linear models
library(mgcv) # generalized additive models
library(gsg) # calculating selection gradients
library(visreg)
library(viridis)

## LOAD & MANAGE DATASET ----

gall_selection.df <- read_csv("gall_selection_data.csv") %>%
  # convert appropriate variables to characters instead of integers
  mutate(Plant_Position = as.character(Plant_Position),
         Gall_Number = as.character(Gall_Number)) %>%
  unite(Gall_ID, Gall_Number, Gall_Letter, remove = FALSE) %>%
  
  # subset data for analysis
  filter(phenology == "early", Location == "tree",
         platy > 0 | ectos > 0 | pupa > 0) %>%                  # eliminate unknown sources of mortality
  mutate(gall_survival = as.numeric(ifelse(pupa > 0, 1, 0)),
         egg_parasitoid = as.numeric(ifelse(platy > 0, 1, 0))) 

## CREATE SEPARATE DATASETS FOR SELECTION ANALYSES ----
control_df <- gall_selection.df %>%
  filter(Treatment.focus == "Control") %>%
  mutate(c.Gall_Height_mm = Gall_Height_mm - mean(Gall_Height_mm),
         c.gall_individuals = gall_individuals - mean(gall_individuals),
         c.Density_per_100_shoots = Density_per_100_shoots - mean(Density_per_100_shoots))

treatment_df <- gall_selection.df %>%
  filter(Treatment.focus == "Ectoparasitoid exclusion") %>%
  mutate(c.Gall_Height_mm = Gall_Height_mm - mean(Gall_Height_mm),
         c.gall_individuals = gall_individuals - mean(gall_individuals),
         c.Density_per_100_shoots = Density_per_100_shoots - mean(Density_per_100_shoots))


## SUMMARIZE DATA AT TREE-LEVEL

## CONTROL PLANTS
control_tree <- control_df %>%
  group_by(Genotype, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_survival), funs(mean, n())) %>%
  ungroup() %>%
  select(-Density_per_100_shoots_n, -Gall_Height_mm_n)

control_tree %>%
  ggplot(., aes(x = Density_per_100_shoots_mean, y=log(gall_survival_mean), size = gall_survival_n)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x) +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ x, color="red")

control_tree %>%
  filter(Density_per_100_shoots_mean < 75) %>% # excluding the outlying density point reveals no relationship between density and mean fitness
  ggplot(., aes(x = Density_per_100_shoots_mean, y=log(gall_survival_mean), size = gall_survival_n)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x) +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ x, color="red")

control_tree %>%
  ggplot(., aes(x = Gall_Height_mm_mean, y=log(gall_survival_mean), size = gall_survival_n)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x) +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ x, color="red") + 
  geom_smooth(method="lm", formula = y ~ poly(x,2), color="orange") +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ poly(x,2), color="purple")

control_tree %>%
  ggplot(., aes(x = Density_per_100_shoots_mean, y=Gall_Height_mm_mean, fill=log(gall_survival_mean), size = gall_survival_n)) +
  geom_point(shape = 21) 

mcpeek_control.lm <- lm(log(gall_survival_mean+0.1) ~ Gall_Height_mm_mean*Density_per_100_shoots_mean, data = control_tree)
summary(mcpeek_control.lm)
anova(mcpeek_control.lm)

## TREATMENT PLANTS
treatment_tree <- treatment_df %>%
  group_by(Genotype, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_survival), funs(mean, n())) %>%
  ungroup() %>%
  select(-Density_per_100_shoots_n, -Gall_Height_mm_n)

treatment_tree %>%
  ggplot(., aes(x = Density_per_100_shoots_mean, y=log(gall_survival_mean), size = gall_survival_n)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x) +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ x, color="red")

treatment_tree %>%
  ggplot(., aes(x = Gall_Height_mm_mean, y=log(gall_survival_mean), size = gall_survival_n)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x) +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ x, color="red") + 
  geom_smooth(method="lm", formula = y ~ poly(x,2), color="orange") +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ poly(x,2), color="purple")

treatment_tree %>%
  ggplot(., aes(x = Density_per_100_shoots_mean, y=Gall_Height_mm_mean, fill=log(gall_survival_mean), size = gall_survival_n)) +
  geom_point(shape = 21) 

mcpeek_treatment.lm <- lm(log(gall_survival_mean+0.1) ~ Gall_Height_mm_mean*Density_per_100_shoots_mean, data = treatment_tree)
summary(mcpeek_treatment.lm)
anova(mcpeek_treatment.lm)


## FULL STATISTICAL ANALYSIS
all_tree <- gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_survival), funs(mean, n())) %>%
  ungroup() %>%
  select(-Density_per_100_shoots_n, -Gall_Height_mm_n) %>%
  mutate(c.Gall_Height_mm_mean = Gall_Height_mm_mean - mean(Gall_Height_mm_mean),
         c.Density_per_100_shoots_mean = Density_per_100_shoots_mean - mean(Density_per_100_shoots_mean))

all_tree %>%
  ggplot(., aes(x = Density_per_100_shoots_mean, y=log(gall_survival_mean+0.1), size = gall_survival_n)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x) +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ x, color="red") +
  facet_wrap(~Treatment.focus)

all_tree %>%
  ggplot(., aes(x = Gall_Height_mm_mean, y=log(gall_survival_mean+0.1), size = gall_survival_n)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x) +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ x, color="red") + 
  facet_wrap(~Treatment.focus) +
  geom_smooth(method="lm", formula = y ~ poly(x,2), color="orange") +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ poly(x,2), color="purple")

mcpeek_all.lm <- lm(log(gall_survival_mean+0.1) ~ c.Gall_Height_mm_mean*c.Density_per_100_shoots_mean*Treatment.focus, 
                    data = all_tree)
summary(mcpeek_all.lm)
anova(mcpeek_all.lm)
car::Anova(mcpeek_all.lm, type=2)
car::Anova(mcpeek_all.lm, type=3)

hist(all_tree$c.Density_per_100_shoots_mean)
hist(all_tree$gall_survival_mean)
table(all_tree$gall_survival_mean)
summary(all_tree$gall_survival_mean)
hist(log(all_tree$gall_survival_mean))
hist(log(all_tree$gall_survival_mean+0.05))
mcpeek_all.wts.lm <- lm(log(gall_survival_mean+0.05) ~ (c.Gall_Height_mm_mean+c.Density_per_100_shoots_mean)*Treatment.focus + c.Gall_Height_mm_mean:c.Density_per_100_shoots_mean, 
                    data = all_tree, weights = gall_survival_n)
car::vif(mcpeek_all.wts.lm)
summary(mcpeek_all.wts.lm)
anova(mcpeek_all.wts.lm)
car::Anova(mcpeek_all.wts.lm, type=2)
car::Anova(mcpeek_all.wts.lm, type=3)

mcpeek_all.wts.geno.lm <- lm(log(gall_survival_mean+0.05) ~ c.Gall_Height_mm_mean*c.Density_per_100_shoots_mean*Treatment.focus*Genotype, 
                        data = all_tree, weights = gall_survival_n)
summary(mcpeek_all.wts.geno.lm)
anova(mcpeek_all.wts.geno.lm)
car::Anova(mcpeek_all.wts.geno.lm, type=2)
car::Anova(mcpeek_all.wts.geno.lm, type=3)
