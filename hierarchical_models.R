
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

# PCA
gall_pca <- princomp(select(gall_selection.df, sc.sqrt.size, sc.log.indiv, sc.sqrt.density), cor = TRUE)
summary(gall_pca)
biplot(gall_pca, choices = c(1,2))
biplot(gall_pca, choices = c(1,3))
biplot(gall_pca, choices = c(2,3))
gall_pca$scores

gall_selection.df <- cbind.data.frame(gall_selection.df, gall_pca$scores)

## SPECIFYING PRIORS ----

glmer.gall_fitness <- glmer(gall_survival ~ Comp.1*Comp.2*Comp.3*Treatment.focus + 
                          (1 | Genotype/Plant_Position/Gall_Number/Gall_ID),
                        data = gall_selection.df,
                        family = binomial(link = "logit"))

## BAYESIAN MODEL ---- 
options(mc.cores = parallel::detectCores()) # use all available cores

brm.gall_formula <- brmsformula(gall_survival ~ sc.sqrt.size*sc.log.indiv*sc.sqrt.density*Treatment.focus + 
                                  (1 | Genotype/Plant_Position/Gall_Number/Gall_ID))
get_prior(brm.gall_formula, data = gall_selection.df)

uninf_priors <- c(set_prior("normal(0,2)", class = "b"), set_prior("normal(0,2)", class = "sd"))

brm.gall_fitness <- brm(gall_survival ~ Comp.1*Comp.2*Comp.3*Treatment.focus + 
                          (1 | Genotype/Plant_Position/Gall_Number/Gall_ID),
                        data = gall_selection.df,
                        prior = uninf_priors,
                        family = bernoulli(link = "logit"),
                        algorithm = "sampling",
                        control = list(adapt_delta = 0.99),
                        chains = 4)
summary(brm.gall_fitness)

brm.gall_fitness <- brm(gall_survival ~ sc.sqrt.size*sc.log.indiv*sc.sqrt.density*Treatment.focus + 
                          (1 | Genotype/Plant_Position/Gall_Number/Gall_ID),
                        data = gall_selection.df,
                        prior = gall_priors,
                        family = bernoulli(link = "logit"),
                        algorithm = "sampling",
                        control = list(adapt_delta = 0.95),
                        chains = 4)
summary(brm.gall_fitness)

