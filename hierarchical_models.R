
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

viridis6 <- function() {
  # colours taken from the viridis package
  c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725")
}

# Gall Density. Also checking whether densities vary among treatments.
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots), mean) %>%
  ggplot(., aes(x = Genotype, y = Density_per_100_shoots, group = Treatment.focus, color = Treatment.focus)) + 
  geom_jitter(alpha = 0.5, width = 0.05) + 
  stat_summary(fun.data = "mean_cl_boot", size = 1.5) 

# Number of Individuals per gall
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position, Gall_Number) %>%
  summarise_at(vars(gall_individuals), mean) %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(gall_individuals), funs(mean, n())) %>%
  ggplot(., aes(x = Genotype, y = mean, group = Treatment.focus, color = Treatment.focus)) + 
  geom_jitter(alpha = 0.5, width = 0.05, aes(size = n)) + 
  stat_summary(fun.data = "mean_cl_boot", size = 1.5)

# Gall Size
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position, Gall_Number) %>%
  summarise_at(vars(Gall_Height_mm), mean) %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Gall_Height_mm), funs(mean, n())) %>%
  ggplot(., aes(x = Genotype, y = mean, group = Treatment.focus, color = Treatment.focus)) + 
  geom_jitter(alpha = 0.5, width = 0.05, aes(size = n)) + 
  stat_summary(fun.data = "mean_cl_boot", size = 1.5)


# Pupation rates
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(pupa, total), sum) %>%
  ggplot(., aes(x = Genotype, y = pupa/total, group = Treatment.focus, color = Treatment.focus)) + 
  geom_jitter(alpha = 0.5, width = 0.05, aes(size = total)) + 
  stat_summary(fun.data = "mean_cl_boot", size = 1.5)

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

## KEY QUESTION: HOW DOES THE LOSS OF FOOD-WEB COMPLEXITY AFFECT SELECTION ON INSECT GALLS?
# Excluding ectoparasitoids should increase gall survival. What is more interesting is whether 
# excluding ectoparasitoids alters the selection gradient.

# gall size
gall_selection.df %>%
  mutate(pupa_survival = ifelse(pupa > 0, 1, 0)) %>%
  ggplot(., aes(x = Gall_Height_mm, y = pupa_survival, color = Treatment.focus)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = "gam", method.args = list(family = "binomial"))# + # binomial_smooth() + 
  #facet_wrap(~Genotype)

# number of gall individuals
gall_selection.df %>%
  mutate(pupa_survival = ifelse(pupa > 0, 1, 0)) %>%
  ggplot(., aes(x = gall_individuals, y = pupa_survival, color = Treatment.focus)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = "gam", method.args = list(family = "binomial")) #+ # binomial_smooth() + 
  #facet_wrap(~Genotype)

# gall density
gall_selection.df %>%
  mutate(pupa_survival = ifelse(pupa > 0, 1, 0)) %>%
  ggplot(., aes(x = Density_per_100_shoots, y = pupa_survival, color = Treatment.focus)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = "gam", method.args = list(family = "binomial")) #+ # binomial_smooth() + 
  #facet_wrap(~Genotype)


## POTENTIAL 3-WAY STATISTICAL INTERACTIONS

# when ectoparasitoids are excluded, large galls on trees at low density have the highest survival
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Gall_Height_mm, Density_per_100_shoots, pupa, total), funs(sum, mean)) %>%
  ggplot(., aes(x = Gall_Height_mm_mean, y = log(Density_per_100_shoots_mean+1), fill = pupa_sum/total_sum)) + 
  geom_point(shape = 21, aes(size = total_sum)) + scale_fill_gradientn(colors = viridis6()) + facet_wrap(~Treatment.focus)

# but, I don't understand why this indicates the opposite...right?
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Gall_Height_mm, Density_per_100_shoots, pupa, total), funs(sum, mean, n())) %>%
  mutate(sc.Gall_Height_mm_mean = scale(Gall_Height_mm_mean),
         sc.Density_per_100_shoots_mean = scale(Density_per_100_shoots_mean),
         sc.gall_Height.x.Dens = sc.Gall_Height_mm_mean*sc.Density_per_100_shoots_mean) %>%
  ggplot(., aes(x = sc.gall_Height.x.Dens, y = pupa_sum/total_sum, weight = total_sum, color = Treatment.focus)) + 
  geom_point(alpha = 0.5, aes(size = total_sum)) + 
  geom_smooth(method = "gam", method.args = list(family = "binomial")) #+ # binomial_smooth() + 
#facet_wrap(~Genotype)

# when ectoparasitoids are excluded, large galls with many individuals have the highest survival...
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Gall_Height_mm, gall_individuals, pupa, total), funs(sum, mean)) %>%
  ggplot(., aes(x = Gall_Height_mm_mean, y = gall_individuals_mean, fill = pupa_sum/total_sum)) + 
  geom_point(shape = 21, aes(size = total_sum)) + scale_fill_gradientn(colors = viridis6()) + facet_wrap(~Treatment.focus)

# need to think about this...
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Gall_Height_mm, gall_individuals, pupa, total), funs(sum, mean)) %>%
  mutate(sc.Gall_Height_mm_mean = scale(Gall_Height_mm_mean),
         sc.gall_individuals_mean = scale(gall_individuals_mean),
         sc.gall_Height.x.Ind = sc.Gall_Height_mm_mean*sc.gall_individuals_mean) %>%
  ggplot(., aes(x = sc.gall_Height.x.Ind, y = pupa_sum/total_sum, weight = gall_individuals_sum, color = Treatment.focus)) + 
  geom_point(alpha = 0.5, aes(size = gall_individuals_sum)) + 
  geom_smooth(method = "gam", method.args = list(family = "binomial")) #+ # binomial_smooth() + 
  #facet_wrap(~Genotype) 

# the interpretation from the continuous scaling isn't all that clear here...
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots, gall_individuals, pupa, total), funs(sum, mean)) %>%
  ggplot(., aes(x = Density_per_100_shoots_mean, y = gall_individuals_mean, fill = pupa_sum/total_sum)) + 
  geom_point(shape = 21, aes(size = total_sum)) + scale_fill_gradientn(colors = viridis6()) + facet_wrap(~Treatment.focus)

# interesting, safety in numbers works in a complex environments, but not simple ones.
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots, gall_individuals, pupa, total), funs(sum, mean)) %>%
  mutate(sc.Density_per_100_shoots_mean = scale(Density_per_100_shoots_mean),
         sc.gall_individuals_mean = scale(gall_individuals_mean),
         sc.gall_Density.x.Ind = sc.Density_per_100_shoots_mean*sc.gall_individuals_mean) %>%
  ggplot(., aes(x = sc.gall_Density.x.Ind, y = pupa_sum/total_sum, weight = total_sum, color = Treatment.focus)) + 
  geom_point(alpha = 0.5, aes(size = total_sum)) + 
  geom_smooth(method = "gam", method.args = list(family = "binomial")) #+ # binomial_smooth() + 
#facet_wrap(~Genotype) 


## PCA
pca.df <- gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position, Gall_Number) %>%
  summarise_at(vars(Gall_Height_mm, gall_individuals, Density_per_100_shoots), mean) %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Gall_Height_mm, gall_individuals, Density_per_100_shoots), funs(mean, sd, sum, n())) %>%
  ungroup()

ggplot(pca.df, aes(x = Genotype, y = Gall_Height_mm_sd/Gall_Height_mm_mean)) + geom_boxplot() # gall size in Genotype J is much more variable
ggplot(pca.df, aes(x = Genotype, y = gall_individuals_sd/gall_individuals_mean)) + geom_boxplot() # not much difference in variability

ggplot(pca.df, aes(x = Gall_Height_mm_mean)) + geom_density()
ggplot(pca.df, aes(x = gall_individuals_mean)) + geom_density() + scale_x_log10() # skewed, log transformation best
ggplot(pca.df, aes(x = Density_per_100_shoots_mean)) + geom_density() + scale_x_sqrt() # skewed, sqrt transformation best

pca.df <- pca.df %>%
  mutate(log.gall_individuals_mean = log(gall_individuals_mean),
         sqrt.Density_per_100_shoots_mean = sqrt(Density_per_100_shoots_mean)) %>%
  mutate(sc.log.gall_individuals_mean = (log.gall_individuals_mean - mean(log.gall_individuals_mean))/sd(log.gall_individuals_mean),
         sc.sqrt.Density_per_100_shoots_mean = (sqrt.Density_per_100_shoots_mean - mean(sqrt.Density_per_100_shoots_mean))/sd(sqrt.Density_per_100_shoots_mean),
         sc.Gall_Height_mm_mean = (Gall_Height_mm_mean - mean(Gall_Height_mm_mean))/sd(Gall_Height_mm_mean))

library(GGally)
ggcorr(select(pca.df, Gall_Height_mm_mean, gall_individuals_mean, Density_per_100_shoots_mean))
ggcorr(select(pca.df, Gall_Height_mm_mean, log.gall_individuals_mean, sqrt.Density_per_100_shoots_mean))

ggplot(pca.df, aes(x = Gall_Height_mm_mean, y = gall_individuals_mean)) + geom_point() + geom_smooth(method = "lm") + scale_y_log10()
ggplot(pca.df, aes(x = Gall_Height_mm_mean, y = Density_per_100_shoots_mean)) + geom_point() + geom_smooth(method = "lm") + scale_y_sqrt()
ggplot(pca.df, aes(x = gall_individuals_mean, y = Density_per_100_shoots_mean)) + geom_point() + geom_smooth(method = "lm") + scale_y_sqrt() + scale_x_log10()


gall_pca <- princomp(~Gall_Height_mm_mean + log.gall_individuals_mean + sqrt.Density_per_100_shoots_mean, data = pca.df, cor = TRUE)
summary(gall_pca)
biplot(gall_pca, choices = c(1,2))
biplot(gall_pca, choices = c(1,3))
loadings(gall_pca)

pca.df <- cbind.data.frame(pca.df, gall_pca$scores)

pupa.df <- gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(pupa, total), sum) %>%
  ungroup() %>%
  left_join(., pca.df)

library(lme4)
library(piecewiseSEM)
pupa.surv <- glmer(pupa/total ~ Comp.1*Comp.2*Comp.3*Treatment.focus + (1 | Genotype), data = pupa.df, family = "binomial", weights = total)
summary(pupa.surv)
sem.model.fits(pupa.surv)
library(visreg)

visreg2d(pupa.surv, xvar = "Comp.1", yvar = "Comp.2", scale = "response", cond = list(Treatment.focus = "Control"))
visreg2d(pupa.surv, xvar = "Comp.1", yvar = "Comp.2", scale = "response", cond = list(Treatment.focus = "Ectoparasitoid exclusion"), plot.type = "persp")

pupa.surv2 <- glmer(pupa/total ~ sc.Gall_Height_mm_mean*sc.log.gall_individuals_mean*sc.sqrt.Density_per_100_shoots_mean*Treatment.focus + (1 | Genotype), 
                    data = pupa.df, family = "binomial", weights = total, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e4)))
summary(pupa.surv2)
sem.model.fits(pupa.surv2)
vcov(pupa.surv2)
library(visreg)

visreg2d(pupa.surv2, yvar = "sc.log.gall_individuals_mean", xvar = "sc.sqrt.Density_per_100_shoots_mean", scale = "response", cond = list(Treatment.focus = "Control"))
visreg2d(pupa.surv2, yvar = "sc.log.gall_individuals_mean", xvar = "sc.sqrt.Density_per_100_shoots_mean", scale = "response", cond = list(Treatment.focus = "Ectoparasitoid exclusion"))

visreg2d(pupa.surv2, xvar = "sc.log.gall_individuals_mean", yvar = "sc.Gall_Height_mm_mean", scale = "response", cond = list(Treatment.focus = "Control"))
visreg2d(pupa.surv2, xvar = "sc.log.gall_individuals_mean", yvar = "sc.Gall_Height_mm_mean", scale = "response", cond = list(Treatment.focus = "Ectoparasitoid exclusion"))

################### MAYBE USEFUL BELOW ###############################
# Does the effect of gall size on pupation rates vary among genotypes?
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Gall_Height_mm, pupa, total), funs(sum, mean)) %>%
  ggplot(., aes(x = Gall_Height_mm_mean, y = pupa_sum/total_sum, group = Genotype, color = Genotype, weight = total_sum)) + 
  geom_point(alpha = 0.5, aes(size = total_sum)) + binomial_smooth(se = FALSE) + facet_wrap(~Treatment.focus) # interesting that Genotype J is qualitatively different in both scenarios.

# Does effect of number of gall individuals on pupation rates vary among genotypes?
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(gall_individuals, pupa, total), funs(sum, mean)) %>%
  ggplot(., aes(x = gall_individuals_mean, y = pupa_sum/total_sum, group = Genotype, color = Genotype, weight = total_sum)) + 
  geom_point(alpha = 0.5, aes(size = total_sum)) + binomial_smooth(se = FALSE) + facet_wrap(~Treatment.focus) # interesting that Genotype J is qualitatively different in both scenarios.

# Does effect of gall density on pupation rates vary among genotypes?
gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots, pupa, total), funs(sum, mean)) %>%
  ggplot(., aes(x = log(Density_per_100_shoots_mean+1), y = pupa_sum/total_sum, group = Genotype, color = Genotype, weight = total_sum)) + 
  geom_point(alpha = 0.5, aes(size = total_sum)) + binomial_smooth(se = FALSE) + facet_wrap(~Treatment.focus)




height_all <- gall_selection.df %>%
  filter(Treatment.focus == "Control") %>%
  group_by(Genotype, Plant_Position, Gall_Number) %>%
  summarise_at(vars(Gall_Height_mm), mean) %>%
  group_by(Genotype, Plant_Position) %>%
  summarise_at(vars(Gall_Height_mm), funs(mean, n())) %>%
  mutate(indicator = "all")

height_survive <- gall_selection.df %>%
  filter(Treatment.focus == "Control") %>%
  mutate(pupa_survive = ifelse(pupa > 0, 1, 0)) %>%
  group_by(Genotype, Plant_Position, Gall_Number, pupa_survive) %>%
  summarise_at(vars(Gall_Height_mm), mean) %>%
  group_by(Genotype, Plant_Position, pupa_survive) %>%
  summarise_at(vars(Gall_Height_mm), funs(mean, n())) %>%
  filter(pupa_survive == 1) %>%
  mutate(indicator = "survive")

mean_heights <- bind_rows(height_all, height_survive) %>%
  group_by(Genotype, indicator) %>%
  summarise(mean = mean(mean))
bind_rows(height_all, height_survive) %>%
  ggplot(., aes(x = mean, group = indicator, fill = indicator)) + 
  geom_density(alpha = 0.5) + facet_wrap(~Genotype) + geom_vline(data = mean_heights, aes(xintercept = mean, linetype = indicator))
