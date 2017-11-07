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

## Effect of genotypes on gall landscape
library(vegan)
geno.landscape <- rda(select(pupa.df, sc.log.gall_individuals_mean, sc.sqrt.Density_per_100_shoots_mean, sc.Gall_Height_mm_mean) ~ Genotype, pupa.df)
summary(geno.landscape)
plot(geno.landscape, display = c("sp","cn"))

ggplot(pupa.df, aes(x = Comp.1, y = Comp.2, color = Genotype)) + geom_point(size = 2)

landscape.summary <- pupa.df %>% 
  group_by(Genotype) %>%
  summarise_at(vars(Comp.1, Comp.2, Comp.3, log.gall_individuals_mean, sqrt.Density_per_100_shoots_mean, Gall_Height_mm_mean), funs(mean, sd)) 

ggplot(landscape.summary, aes(x = Comp.1_mean, y = Comp.2_mean, color = Genotype)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin = Comp.1_mean - Comp.1_sd, xmax = Comp.1_mean + Comp.1_sd)) + 
  geom_errorbar(aes(ymin = Comp.2_mean - Comp.2_sd, ymax = Comp.2_mean + Comp.2_sd))

ggplot(landscape.summary, aes(x = Comp.1_mean, y = Comp.3_mean, color = Genotype)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin = Comp.1_mean - Comp.1_sd, xmax = Comp.1_mean + Comp.1_sd)) + 
  geom_errorbar(aes(ymin = Comp.3_mean - Comp.3_sd, ymax = Comp.3_mean + Comp.3_sd))

ggplot(landscape.summary, aes(x = sqrt.Density_per_100_shoots_mean_mean, y = Gall_Height_mm_mean_mean, color = Genotype)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin = sqrt.Density_per_100_shoots_mean_mean - sqrt.Density_per_100_shoots_mean_sd, xmax = sqrt.Density_per_100_shoots_mean_mean + sqrt.Density_per_100_shoots_mean_sd)) + 
  geom_errorbar(aes(ymin = Gall_Height_mm_mean_mean - Gall_Height_mm_mean_sd, ymax = Gall_Height_mm_mean_mean + Gall_Height_mm_mean_sd))

ggplot(landscape.summary, aes(x = sqrt.Density_per_100_shoots_mean_mean, y = log.gall_individuals_mean_mean, color = Genotype)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin = sqrt.Density_per_100_shoots_mean_mean - sqrt.Density_per_100_shoots_mean_sd, xmax = sqrt.Density_per_100_shoots_mean_mean + sqrt.Density_per_100_shoots_mean_sd)) + 
  geom_errorbar(aes(ymin = log.gall_individuals_mean_mean - log.gall_individuals_mean_sd, ymax = log.gall_individuals_mean_mean + log.gall_individuals_mean_sd))

ggplot(landscape.summary, aes(x = Gall_Height_mm_mean_mean, y = log.gall_individuals_mean_mean, color = Genotype)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin = Gall_Height_mm_mean_mean - Gall_Height_mm_mean_sd, xmax = Gall_Height_mm_mean_mean + Gall_Height_mm_mean_sd)) + 
  geom_errorbar(aes(ymin = log.gall_individuals_mean_mean - log.gall_individuals_mean_sd, ymax = log.gall_individuals_mean_mean + log.gall_individuals_mean_sd))

## possible regression models
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


## Bayesian Model

brm.gall_fitness <- brm(gall_survival ~ sc.size*sc.indiv*sc.density*Treatment.focus + 
                          (1 | Genotype/Plant_Position/Gall_Number/Gall_ID),
                        data = gall_selection.df,
                        family = bernoulli(link = "logit"),
                        algorithm = "meanfield")
summary(brm.gall_fitness)



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
