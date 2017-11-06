
## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
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
  mutate(gall_survival = ifelse(pupa > 0, 1, 0),
         platy_bin = ifelse(platy > 0, 1,0),
         ectos_bin = ifelse(ectos > 0, 1,0))

## EXPLORATORY PLOTS ----

# Plot helpers
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

# explored gam smoothers but they were way too curvy to be informative.
#gam_smooth <- function(...) {
#  geom_smooth(method = "gam", method.args = list(family = "binomial"), formula = y ~ s(x))
#}

viridis6 <- function() {
  # colours taken from the viridis package
  c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725")
}

dodge <- position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.5)


## EXPLORING EFFECT OF GENOTYPE ON GALL PHENOTYPES, SURVIVAL, AND PARASITOID RESPONSES ----

# Gall Density
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots), mean) %>%
  ggplot(., aes(x = Genotype, y = Density_per_100_shoots, group = Treatment.focus, color = Treatment.focus)) + 
  geom_jitter(alpha = 0.5, position = dodge) + 
  stat_summary(fun.data = "mean_cl_boot", size = 1.5, position = dodge) 


# Number of Individuals per gall
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position, Gall_Number) %>%
  summarise_at(vars(gall_individuals), mean) %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(gall_individuals), funs(mean, n())) %>%
  ggplot(., aes(x = Genotype, y = mean, group = Treatment.focus, color = Treatment.focus)) + 
  geom_jitter(alpha = 0.5, position = dodge, aes(size = n)) + 
  stat_summary(fun.data = "mean_cl_boot", size = 1.5, position = dodge)

# Gall Size
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position, Gall_Number) %>%
  summarise_at(vars(Gall_Height_mm), mean) %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Gall_Height_mm), funs(mean, n())) %>%
  ggplot(., aes(x = Genotype, y = mean, group = Treatment.focus, color = Treatment.focus)) + 
  geom_jitter(alpha = 0.5, position = dodge, aes(size = n)) + 
  stat_summary(fun.data = "mean_cl_boot", size = 1.5, position = dodge)


# Pupation rates
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(pupa, total), sum) %>%
  ggplot(., aes(x = Genotype, y = pupa/total, group = Treatment.focus, color = Treatment.focus)) + 
  geom_jitter(alpha = 0.5, position = dodge, aes(size = total)) + 
  stat_summary(fun.data = "mean_cl_boot", size = 1.5, position = dodge)

# Ectoparasitism rates and estimates of how effective the exclusion treatment was
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(ectos, total), sum) %>%
  ggplot(., aes(x = Genotype, y = ectos/total, group = Treatment.focus, color = Treatment.focus)) + 
  geom_jitter(alpha = 0.5, position = dodge, aes(size = total)) + 
  stat_summary(fun.data = "mean_cl_boot", size = 1.5, position = dodge)

# Platygaster parasitism rates and estimates of intraguild predation 
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(platy, total), sum) %>%
  ggplot(., aes(x = Genotype, y = platy/total, group = Treatment.focus, color = Treatment.focus)) + 
  geom_jitter(alpha = 0.5, position = dodge, aes(size = total)) + 
  stat_summary(fun.data = "mean_cl_boot", size = 1.5, position = dodge)

## HOW DOES THE LOSS OF FOOD-WEB COMPLEXITY AFFECT GALL SURVIVAL? ----

# create discrete predictors from continuous variables to examine complex interactions
size.cut <- quantile(gall_selection.df$Gall_Height_mm, probs = c(0,0.33,0.66,1))
density.cut <- quantile(gall_selection.df$Density_per_100_shoots, probs = c(0,0.33,0.66,1))
indiv.cut <- quantile(gall_selection.df$gall_individuals, probs = c(0,0.33,0.66,1))

get_mean.size <- gall_selection.df %>%
  group_by(Plant_Position, Gall_Number) %>%
  summarise_at(vars(Gall_Height_mm), mean) %>%
  ungroup() %>%
  rename(mean.gall.size = Gall_Height_mm)

gall_selection.df <- mutate(gall_selection.df, 
                            size.cut = cut(Gall_Height_mm, size.cut, include.lowest = T),
                            density.cut = cut(Density_per_100_shoots, density.cut, include.lowest = T),
                            indiv.cut = cut(gall_individuals, indiv.cut, include.lowest = T))
gall_selection.df <- left_join(gall_selection.df, get_mean.size) %>%
  mutate(size.enhance = mean.gall.size - Gall_Height_mm)

summary(gall_selection.df$size.enhance)
gall_selection.df$size.enhance.cut <- cut(gall_selection.df$size.enhance, breaks = c(min(gall_selection.df$size.enhance), 0, max(gall_selection.df$size.enhance)), include.lowest = T)

#### GALL SURVIVAL ----
parasitism.df <- gall_selection.df %>%
  filter(platy > 0 | ectos > 0 | pupa > 0) # remove unknown sources of mortality

## GALL SIZE 
base_size <- ggplot(parasitism.df, aes(x = Gall_Height_mm, y = gall_survival, color = Treatment.focus)) + 
  geom_point(alpha = 0.5) + scale_x_sqrt()

# larger galls enhance survival.
# effect of ectoparasitoid exclusion is enhanced in smaller galls
base_size + binomial_smooth(formula = y ~ poly(x, 2)) # no evidence of peak
base_size + binomial_smooth() 

# effect of ectoparasitoid exclusion goes away at high gall densities
base_size + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~density.cut)
base_size + binomial_smooth() + facet_wrap(~density.cut)

# effect of ectoparasitoid exclusion happens primarily in galls with few neighbors
base_size + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~indiv.cut)
base_size + binomial_smooth() + facet_wrap(~indiv.cut)

# weird data at gall densities with many individuals.
# doesn't make anysense that gall survival in higher in the control treatment
base_size + binomial_smooth(formula = y ~ poly(x, 2)) + facet_grid(density.cut~indiv.cut)
base_size + binomial_smooth() + facet_grid(density.cut~indiv.cut)


## GALL NEIGHBORS

base_indiv <- ggplot(parasitism.df, aes(x = gall_individuals, y = gall_survival, color = Treatment.focus)) + 
  geom_point(alpha = 0.5) + scale_x_log10()

# contrasting relationships depending on treatment
# but why is survival higher in the complex treatment with many individuals?
base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) 
base_indiv + binomial_smooth() 

# why is gall survival so low at high gall densities when ectoparasitoids are excluded?
# Just because parasitoids have a positive response doesn't mean they should have 
# lower survival in the control treatment...
base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~density.cut)
base_indiv + binomial_smooth() + facet_wrap(~density.cut)

# gall size doesn't appear to modify the relationship between number of neighbors and treatment.
base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~size.cut)
base_indiv + binomial_smooth() + facet_wrap(~size.cut)

# I find it weird that survival can be lower in the exclusion treatment at high densities...
base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) + facet_grid(density.cut~size.cut)
base_indiv + binomial_smooth() + facet_grid(density.cut~size.cut)

## GALL DENSITY 
base_density <- ggplot(parasitism.df, aes(x = Density_per_100_shoots, y = gall_survival, color = Treatment.focus)) + 
  geom_point(alpha = 0.5) + scale_x_sqrt()

# gall density has a clear negative effect in the simple food web.
# in the complex web, there is either no relationship or a divergent one (better to be at high or low densities)
base_density + binomial_smooth(formula = y ~ poly(x, 2))
base_density + binomial_smooth() 

# weird, super strong negative relationship with many individuals...
base_density + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~indiv.cut)
base_density + binomial_smooth() + facet_wrap(~indiv.cut)

# gall size doesn't appear to alter the relationship between density and survival
base_density + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~size.cut)
base_density + binomial_smooth() + facet_wrap(~size.cut)

# no clear 4-way interaction
base_density + binomial_smooth(formula = y ~ poly(x, 2)) + facet_grid(indiv.cut~size.cut)
base_density + binomial_smooth() + facet_grid(indiv.cut~size.cut)

#### PLATYGASTER PARASITISM ----

## GALL SIZE 
platy.base_size <- ggplot(parasitism.df, aes(x = Gall_Height_mm, y = platy_bin, color = Treatment.focus)) + 
  geom_point(alpha = 0.5)

# intraguild predation occurs primarily in smaller galls
platy.base_size + binomial_smooth(formula = y ~ poly(x, 2)) 
platy.base_size + binomial_smooth() 

# no clear dependence on gall density
platy.base_size + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~density.cut)
platy.base_size + binomial_smooth() + facet_wrap(~density.cut)

# intraguild predation may be more apparent in galls with many individuals
platy.base_size + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~indiv.cut)
platy.base_size + binomial_smooth() + facet_wrap(~indiv.cut)

# error moargins are too large to explore 4-way interactions
platy.base_size + binomial_smooth(formula = y ~ poly(x, 2)) + facet_grid(density.cut~indiv.cut)
platy.base_size + binomial_smooth() + facet_grid(density.cut~indiv.cut)


## GALL NEIGHBORS
platy.base_indiv <- ggplot(parasitism.df, aes(x = gall_individuals, y = platy_bin, color = Treatment.focus)) + 
  geom_point(alpha = 0.5)

# Interesting!
# Same pattern as gall density below. 
platy.base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) 
platy.base_indiv + binomial_smooth() 

# Platygaster response to gall neighbors is exacerbated at high gall densities
# don't understand what is happening at low densities. How could Platygaster parasitism be higher in the control treatment?
platy.base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~density.cut) 
platy.base_indiv + binomial_smooth() + facet_wrap(~density.cut)

# Platygaster response to gall neighbors goes away in large galls
# Perhaps this is because large galls generally provide a refuge from parasitoids...
platy.base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~size.cut)
platy.base_indiv + binomial_smooth() + facet_wrap(~size.cut)

# no clear evidence for a 4-way interaction
# The weird low density pattern is only happening for the largest galls...
platy.base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) + facet_grid(density.cut~size.cut)
platy.base_indiv + binomial_smooth() + facet_grid(density.cut~size.cut)

## GALL DENSITY
platy.base_density <- ggplot(parasitism.df, aes(x = Density_per_100_shoots, y = platy_bin, color = Treatment.focus)) + 
  geom_point(alpha = 0.5) + scale_x_sqrt()

# Interesting!
# In the simple food-web, Platygaster parasitism increases at higher gall densities. 
# However, Platygaster has an optimal intermediate density in the more complex food web. 
# I suspect that this is due to increases in intraguild predation at higher densities.
platy.base_density + binomial_smooth(formula = y ~ poly(x, 2))
platy.base_density + binomial_smooth() 

# above pattern does not depend on gall size, although parasitism is less in larger galls
platy.base_density + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~size.cut)
platy.base_density + binomial_smooth() + facet_wrap(~size.cut)

# above pattern may be intensified in galls with many individuals.
platy.base_density + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~indiv.cut)
platy.base_density + binomial_smooth() + facet_wrap(~indiv.cut)

# margins of error are too wide to suggest much. 
# i.e. not worth exploring a 4-way interaction
platy.base_density + binomial_smooth(formula = y ~ poly(x, 2)) + facet_grid(indiv.cut~size.cut)
platy.base_density + binomial_smooth() + facet_grid(indiv.cut~size.cut)


#### ECTO- PARASITISM ----

ectos.df <- filter(parasitism.df, Treatment.focus == "Control")

## GALL SIZE 
ecto.base_size <- ggplot(ectos.df, aes(x = Gall_Height_mm, y = ectos_bin)) + 
  geom_point(alpha = 0.5)

# ectoparasitoid parasitism decreases linearly with gall size
ecto.base_size + binomial_smooth(formula = y ~ poly(x, 2)) 
ecto.base_size + binomial_smooth() 

# relationship doesn't appear to depend on gall density
ecto.base_size + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~density.cut)
ecto.base_size + binomial_smooth() + facet_wrap(~density.cut)

# pretty wide intervals, tough to judge this relationship
ecto.base_size + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~indiv.cut)
ecto.base_size + binomial_smooth() + facet_wrap(~indiv.cut)

# possible positive relationship with gall size at intermediate densities and individuals...
ecto.base_size + binomial_smooth(formula = y ~ poly(x, 2)) + facet_grid(density.cut~indiv.cut)
ecto.base_size + binomial_smooth() + facet_grid(density.cut~indiv.cut)


## GALL NEIGHBORS
ecto.base_indiv <- ggplot(ectos.df, aes(x = gall_individuals, y = ectos_bin)) + 
  geom_point(alpha = 0.5)

# ectoparasitoid parasitism decreases linearly with number of gall individuals
ecto.base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) 
ecto.base_indiv + binomial_smooth() 

# no clear dependence on gall density
ecto.base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~density.cut)
ecto.base_indiv + binomial_smooth() + facet_wrap(~density.cut)

# no clear dependence on gall size
ecto.base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~size.cut)
ecto.base_indiv + binomial_smooth() + facet_wrap(~size.cut)

# unclear
ecto.base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) + facet_grid(density.cut~size.cut)
ecto.base_indiv + binomial_smooth() + facet_grid(density.cut~size.cut)

## GALL DENSITY
ecto.base_density <- ggplot(ectos.df, aes(x = Density_per_100_shoots, y = ectos_bin)) + 
  geom_point(alpha = 0.5) + scale_x_sqrt()

# ectoparasitism increases at higher densities
ecto.base_density + binomial_smooth(formula = y ~ poly(x, 2))
ecto.base_density + binomial_smooth() 

# relationship is highest at intermediate gall neighbors
ecto.base_density + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~indiv.cut)
ecto.base_density + binomial_smooth() + facet_wrap(~indiv.cut)

# relationship doesn't appear to depend on gall size
ecto.base_density + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~size.cut)
ecto.base_density + binomial_smooth() + facet_wrap(~size.cut)

# difficult to interpret...
ecto.base_density + binomial_smooth(formula = y ~ poly(x, 2)) + facet_grid(indiv.cut~size.cut)
ecto.base_density + binomial_smooth() + facet_grid(indiv.cut~size.cut)


#### PLANT-LEVEL GALL SURVIVAL ----
mean.narm <- function(x) mean(x, na.rm = TRUE)
sum.narm <- function(x) sum(x, na.rm = TRUE)

plant_level.info <- gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position, Gall_Number) %>%
  summarise_at(vars(Density_per_100_shoots, gall_individuals, Gall_Height_mm), mean.narm) %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots, gall_individuals, Gall_Height_mm), mean.narm) %>%
  ungroup()

plant_level.parasitism <- parasitism.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(pupa, platy, ectos, platy.ectos, total), sum.narm) %>%
  ungroup()

plant_level.df <- left_join(plant_level.info, plant_level.parasitism)

plant_level.df$density.cut <- cut(plant_level.df$Density_per_100_shoots, 
                                  quantile(plant_level.df$Density_per_100_shoots, c(0,0.5,1)),
                                  include.lowest = T)
plant_level.df$size.cut <- cut(plant_level.df$Gall_Height_mm, 
                                  quantile(plant_level.df$Gall_Height_mm, c(0,0.5,1)),
                                  include.lowest = T)
plant_level.df$indiv.cut <- cut(plant_level.df$gall_individuals, 
                                  quantile(plant_level.df$gall_individuals, c(0,0.5,1)),
                                  include.lowest = T)

## GALL SIZE 
base_size <- ggplot(plant_level.df, aes(x = log(Gall_Height_mm), y = pupa/total, color = Treatment.focus, size = total, weight = total)) + 
  geom_point(alpha = 0.5) 

# interesting, effect of exclusion increases with gall size...opposite of individual level pattern
base_size + binomial_smooth(formula = y ~ poly(x, 2)) 
base_size + binomial_smooth() 

# treatment effect goes away at high gall densities
base_size + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~density.cut)
base_size + binomial_smooth() + facet_wrap(~density.cut)

# interactive effect with many individuals
base_size + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~indiv.cut)
base_size + binomial_smooth() + facet_wrap(~indiv.cut)

# no clear 4-way interaction
base_size + binomial_smooth(formula = y ~ poly(x, 2)) + facet_grid(density.cut~indiv.cut)
base_size + binomial_smooth() + facet_grid(density.cut~indiv.cut)


## GALL NEIGHBORS

base_indiv <- ggplot(plant_level.df, aes(x = log(gall_individuals), y = pupa/total, color = Treatment.focus, size = total, weight = total)) + 
  geom_point(alpha = 0.5)

# similar relationship at individual level
base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) 
base_indiv + binomial_smooth() 

# effect of treatment reverses at high densities
base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~density.cut)
base_indiv + binomial_smooth() + facet_wrap(~density.cut)

# effect of treatment reverses at large gall sizes
base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~size.cut)
base_indiv + binomial_smooth() + facet_wrap(~size.cut)

# hmmm. Effect of neighbors is qualitatively changes when we remove ectoparasitoids, but primarily in large galls at low densities.
base_indiv + binomial_smooth(formula = y ~ poly(x, 2)) + facet_grid(density.cut~size.cut)
base_indiv + binomial_smooth() + facet_grid(density.cut~size.cut)

## GALL DENSITY 
base_density <- ggplot(plant_level.df, aes(x = log(Density_per_100_shoots+1), y = pupa/total, color = Treatment.focus, size = total, weight = total)) + 
  geom_point(alpha = 0.5) 

# same relationship at individual level
base_density + binomial_smooth(formula = y ~ poly(x, 2))
base_density + binomial_smooth() 

# ?
base_density + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~indiv.cut)
base_density + binomial_smooth() + facet_wrap(~indiv.cut)

# ?
base_density + binomial_smooth(formula = y ~ poly(x, 2)) + facet_wrap(~size.cut)
base_density + binomial_smooth() + facet_wrap(~size.cut)

# ?
base_density + binomial_smooth(formula = y ~ poly(x, 2)) + facet_grid(indiv.cut~size.cut)
base_density + binomial_smooth() + facet_grid(indiv.cut~size.cut)
