
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
  filter(phenology == "early", Location == "tree") %>%
  mutate(gall_survival = ifelse(pupa > 0, 1, 0))

## EXPLORE DISTRIBUTIONS ----
ggplot(gall_selection.df, aes(x = Gall_Height_mm)) + geom_histogram() + scale_x_sqrt()
summary(gall_selection.df$Gall_Height_mm)
summary(sqrt(gall_selection.df$Gall_Height_mm)) # most symmetric distribution around median and mean.
summary(log(gall_selection.df$Gall_Height_mm))

gall_selection.df %>%
  distinct(Genotype, Plant_Position, Gall_Number, gall_individuals) %>% # avoids pseudo replication
  ggplot(., aes(x = gall_individuals)) + geom_histogram() + scale_x_log10()
summary(gall_selection.df$gall_individuals)
summary(sqrt(gall_selection.df$gall_individuals))
summary(log(gall_selection.df$gall_individuals)) # most symmetric distribution around median and mean.

gall_selection.df %>%
  distinct(Genotype, Plant_Position, Density_per_100_shoots) %>% # avoids pseudo replication
  ggplot(., aes(x = sqrt(Density_per_100_shoots))) + geom_histogram() # more normal, except for low values
summary(gall_selection.df$Density_per_100_shoots)
summary(sqrt(gall_selection.df$Density_per_100_shoots)) # most symmetric distribution around median and mean
summary(log(gall_selection.df$Density_per_100_shoots+1))

gall_selection.df <- mutate(gall_selection.df,
                            sc.sqrt.size = scale(sqrt(Gall_Height_mm)),
                            sc.log.indiv = scale(log(gall_individuals)),
                            sc.sqrt.density = scale(sqrt(Density_per_100_shoots)))

# correlations among predictors are virtually absent, except for between gall size and gall density
car::scatterplotMatrix(select(gall_selection.df, sc.sqrt.size, sc.log.indiv, sc.sqrt.density)) # closer correspondence between loess and linear fits
car::scatterplotMatrix(select(gall_selection.df, Gall_Height_mm, gall_individuals, Density_per_100_shoots)) # similar pattern, but more deviance between loess and linear fits

## SPECIFYING PRIORS ----

## BAYESIAN MODEL ---- 
options(mc.cores = parallel::detectCores()) # use all available cores

brm.gall_fitness <- brm(gall_survival ~ sc.sqrt.size*sc.log.indiv*sc.sqrt.density*Treatment.focus + 
                          (1 | Genotype/Plant_Position/Gall_Number/Gall_ID),
                        data = gall_selection.df,
                        family = bernoulli(link = "logit"),
                        algorithm = "meanfield")
summary(brm.gall_fitness)

