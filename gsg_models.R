## LOAD REQUIRED LIBRARIES & SET OPTIONS ----
library(tidyverse)
library(gsg)

## LOAD & MANAGE DATASET ----
gall_selection.df <- read_csv("gall_selection_data.csv") %>%
  mutate(Plant_Position = as.factor(Plant_Position),
         Gall_Number = as.factor(Gall_Number),
         Treatment.focus = as.factor(Treatment.focus),
         Treatment = as.factor(Treatment),
         Genotype = as.factor(Genotype),
         Gall_Letter = as.factor(Gall_Letter)) %>%
  unite(Gall_ID, Gall_Number, Gall_Letter, remove = FALSE) %>%
  
  # subset data for analysis
  filter(phenology == "early", Location == "tree",
         platy > 0 | ectos > 0 | pupa > 0) %>%          # eliminate unknown sources of mortality
  mutate(sc.sqrt.size = scale(sqrt(Gall_Height_mm)),
         sc.log.indiv = scale(log(gall_individuals)),
         sc.sqrt.density = scale(sqrt(Density_per_100_shoots)))


## SUMMARISE AT PLANT LEVEL ----
mean.narm <- function(x) mean(x, na.rm = TRUE)

plant_level.info <- gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position, Gall_Number) %>%
  summarise_at(vars(Density_per_100_shoots, gall_individuals, Gall_Height_mm), mean.narm) %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots, gall_individuals, Gall_Height_mm), mean.narm) %>%
  ungroup()

plant_level.parasitism <- gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(pupa, platy, ectos, platy.ectos, total), sum) %>%
  ungroup()

plant_level.df <- left_join(plant_level.info, plant_level.parasitism) %>%
  mutate(gall_fitness = pupa/total,
         sc.size = (sqrt(Gall_Height_mm) - mean(sqrt(Gall_Height_mm)))/sd(sqrt(Gall_Height_mm)),
         sc.indiv = (log(gall_individuals) - mean(log(gall_individuals)))/sd(log(gall_individuals)),
         sc.density = (sqrt(Density_per_100_shoots) - mean(sqrt(Density_per_100_shoots)))/sd(sqrt(Density_per_100_shoots)))

## GAMM ----

control.df <- filter(gall_selection.df, Treatment.focus == "Control")
library(gamm4)
gall.gam <- gamm4(pupa/total ~ s(Gall_Height_mm, Density_per_100_shoots, gall_individuals),
                  data = control.df,
                  family = binomial,
                  weights = control.df$total,
                  random = ~ (1/Genotype/Plant_Position/Gall_Number))
summary(gall.gam)
gall.gam$mer
gall.gam$gam
plot(gall.gam$gam)
vis.gam(gall.gam$gam, view = c("Gall_Height_mm","Density_per_100_shoots"), type = "response", plot.type = "contour")
vis.gam(gall.gam$gam, view = c("Gall_Height_mm","gall_individuals"), type = "response", plot.type = "contour")
vis.gam(gall.gam$gam, view = c("Density_per_100_shoots","gall_individuals"), type = "response", plot.type = "contour")
