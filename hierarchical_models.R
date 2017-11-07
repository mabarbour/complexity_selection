
## LOAD REQUIRED LIBRARIES & SET OPTIONS ----
library(tidyverse)
library(brms)
library(cowplot)
library(bayesplot)

options(mc.cores = parallel::detectCores()) # use all available cores

## LOAD DATASET ----
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
         platy > 0 | ectos > 0 | pupa > 0) %>%     # eliminate unknown sources of mortality
  mutate(gall_survival = ifelse(pupa > 0, 1, 0))

## EXPLORE DISTRIBUTIONS OF PHENOTYPIC TRAITS.
# important implications for selection gradient analysis

trans.ind <- log(gall_selection.df$gall_individuals) # closes transformation to normality
hist(trans.ind)
qqnorm(trans.ind)
qqline(trans.ind, col = "red")
shapiro.test(trans.ind)

trans.size <- gall_selection.df$Gall_Height_mm # no transformation needed
hist(trans.size)
qqnorm(trans.size)
qqline(trans.size, col = "red") 
shapiro.test(trans.size) # doesn't pass the test though

trans.density <- sqrt(gall_selection.df$Density_per_100_shoots) # sqrt is the closest
hist(trans.density)
qqnorm(trans.density)
qqline(trans.density, col = "red") 
shapiro.test(trans.density)

gall_selection.df <- mutate(gall_selection.df,
                            sc.size = scale(Gall_Height_mm),
                            sc.log.ind = scale(log(gall_individuals)),
                            sc.sqrt.density = scale(sqrt(Density_per_100_shoots)))

## SPECIFYING PRIORS ----

uninf_priors <- c(set_prior("normal(0,1)", class = "b"), set_prior("normal(0,2)", class = "sd")) # generic, weakly informative priors

## BAYESIAN MODEL ---- 

# Simplest model predicts that gall survival can actually be worse when ectoparasitoids were excluded.
# This primarily happens at high gall densities and many individuals. While I feel the mortality could be high
# in this scenario, it shouldn't drop below the control, which is what worries me.
# My solution is to try non-linear responses.

# So far, transformed variables are very important for generating more realistic model predictions, although it is still a bit off for individuals

control.df <- filter(gall_selection.df, Treatment.focus == "Control")
treatment.df <- filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion")

gall_fitness3 <- brm(gall_survival ~ sc.size*sc.log.ind*sc.sqrt.density + 
                      (1 | Genotype/Plant_Position/Gall_Number),
                    data = treatment.df, prior = uninf_priors, family = bernoulli(link = "logit"),
                    algorithm = "sampling", chains = 1, control = list(adapt_delta = 0.99))
summary(gall_fitness3)

LOO(gall_fitness, gall_fitness2, gall_fitness3)
# POSTERIOR PREDICTIVE CHECKING
y <- treatment.df$gall_survival
yrep <- posterior_predict(gall_fitness)

ppc_bars(y, yrep)
ppc_bars_grouped(y, yrep, group = treatment.df$Genotype)

# MARGINAL EFFECTS
plot(marginal_effects(gall_fitness), points = TRUE)

plot(marginal_effects(gall_fitness3, effect = "sc.sqrt.density:sc.log.ind", surface = T), 
     stype = "raster") 

plot(marginal_effects(gall_fitness3, effect = "sc.sqrt.density:sc.size", surface = T), 
     stype = "raster") 

plot(marginal_effects(gall_fitness3, effect = "sc.log.ind:sc.size", surface = T), 
     stype = "raster") 


