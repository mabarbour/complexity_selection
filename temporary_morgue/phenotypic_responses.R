
## LOAD REQUIRED LIBRARIES & SET OPTIONS ----

library(tidyverse)
library(mgcv)
library(gamm4) # for generalized additive mixed models #
library(gsg)
library(cowplot) # pretty default ggplots

#################   BE CAREFUL IN HOW DATA IS SUBSETTED !!! ###################

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
  mutate(gall_survival = as.numeric(ifelse(pupa > 0, 1, 0)),
         egg_parasitoid = as.numeric(ifelse(platy > 0, 1, 0)),
         sc.Gall_Height_mm = as.numeric(scale(sqrt(Gall_Height_mm))),
         sc.gall_individuals = as.numeric(scale(log(gall_individuals))),
         sc.Density_per_100_shoots = as.numeric(scale(sqrt(Density_per_100_shoots))))

filter(gall_selection.df, pupa == 1) %>%
  group_by(Treatment.focus) %>%
  summarise_at(vars(sc_Gall_Height_mm, sc_gall_individuals, sc_Density_per_100_shoots), funs(mean,sd)) %>%
  data.frame()


filter(gall_selection.df, pupa == 1) %>%
  ggplot(., aes(x=Gall_Height_mm, fill=Treatment.focus)) + geom_density(alpha=0.5)


## EXPLORE WITH PHENOTYPES AS THE RESPONSE ----

#test.control <- filter(gall_selection.df, Treatment.focus == "Control") %>%
#  mutate(selection = "before") %>%
#  bind_rows(., mutate(filter(gall_selection.df, Treatment.focus == "Control", gall_survival == 1), selection = "after")) 

#library(brms)
#gall_traits <- brm(cbind(sc_Gall_Height_mm, sc_gall_individuals, sc_Density_per_100_shoots) ~ selection + (1|Genotype/Plant_Position/Gall_Number), data = test.control)

## Gall size and individuals

filter(gall_selection.df, Treatment.focus == "Control") %>%
  ggplot(., aes(x=sc.Gall_Height_mm, y=sc.gall_individuals)) + stat_density_2d(aes(fill = ..level..), geom="polygon")

filter(gall_selection.df, Treatment.focus == "Control", gall_survival == 1) %>%
  ggplot(., aes(x=sc.Gall_Height_mm, y=sc.gall_individuals)) + stat_density_2d(aes(fill = ..level..), geom="polygon")

filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion") %>%
  ggplot(., aes(x=sc.Gall_Height_mm, y=sc.gall_individuals)) + stat_density_2d(aes(fill = ..level..), geom="polygon")

filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion", gall_survival == 1) %>%
  ggplot(., aes(x=sc.Gall_Height_mm, y=sc.gall_individuals)) + stat_density_2d(aes(fill = ..level..), geom="polygon")

## Gall size and density
filter(gall_selection.df, Treatment.focus == "Control") %>%
  ggplot(., aes(x=sc.Gall_Height_mm, y=sc.Density_per_100_shoots)) + stat_density_2d(aes(fill = ..level..), geom="polygon")

filter(gall_selection.df, Treatment.focus == "Control", gall_survival == 1) %>%
  ggplot(., aes(x=sc.Gall_Height_mm, y=sc.Density_per_100_shoots)) + stat_density_2d(aes(fill = ..level..), geom="polygon")

filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion") %>%
  ggplot(., aes(x=sc.Gall_Height_mm, y=sc.Density_per_100_shoots)) + stat_density_2d(aes(fill = ..level..), geom="polygon")

filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion", gall_survival == 1) %>%
  ggplot(., aes(x=sc.Gall_Height_mm, y=sc.Density_per_100_shoots)) + stat_density_2d(aes(fill = ..level..), geom="polygon")

## Gall individuals and density
filter(gall_selection.df, Treatment.focus == "Control") %>%
  ggplot(., aes(x=sc.gall_individuals, y=sc.Density_per_100_shoots)) + stat_density_2d(aes(fill = ..level..), geom="polygon") 

filter(gall_selection.df, Treatment.focus == "Control", gall_survival == 1) %>%
  ggplot(., aes(x=sc.gall_individuals, y=sc.Density_per_100_shoots)) + stat_density_2d(aes(fill = ..level..), geom="polygon")

filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion") %>%
  ggplot(., aes(x=sc.gall_individuals, y=sc.Density_per_100_shoots)) + stat_density_2d(aes(fill = ..level..), geom="polygon")

filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion", gall_survival == 1) %>%
  ggplot(., aes(x=sc.gall_individuals, y=sc.Density_per_100_shoots)) + stat_density_2d(aes(fill = ..level..), geom="polygon")
