

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

#### OLD BUT MAYBE USEFUL ----
summary(gall_selection.df$sc.gall_individuals)
ggplot(gall_selection.df, aes(x = sc.Gall_Height_mm, y = gall_survival, color = Treatment.focus)) +
  geom_point() + binomial_smooth() +
  facet_wrap(~cut(sc.gall_individuals, breaks = c(-5,0,5)))

summary(gall_selection.df$sc.Density_per_100_shoots)
ggplot(gall_selection.df, aes(x = sc.Gall_Height_mm, y = gall_survival, color = Treatment.focus)) +
  geom_point() + binomial_smooth() +
  facet_wrap(~cut(sc.Density_per_100_shoots, breaks = c(-5,0,5)))

summary(gall_selection.df$sc.Gall_Height_mm)
ggplot(gall_selection.df, aes(x = sc.gall_individuals, y = gall_survival, color = Treatment.focus)) +
  geom_point() + binomial_smooth() +
  facet_wrap(~cut(sc.Gall_Height_mm, breaks = c(-5,0,5)))

summary(gall_selection.df$sc.Density_per_100_shoots)
ggplot(gall_selection.df, aes(x = sc.gall_individuals, y = gall_survival, color = Treatment.focus)) +
  geom_point() + binomial_smooth() +
  facet_wrap(~cut(sc.Density_per_100_shoots, breaks = c(-5,0,5)))

summary(gall_selection.df$sc.gall_individuals)
ggplot(gall_selection.df, aes(x = sc.Density_per_100_shoots, y = gall_survival, color = Treatment.focus)) +
  geom_point() + binomial_smooth() +
  facet_wrap(~cut(sc.gall_individuals, breaks = c(-5,0,5)))

summary(gall_selection.df$sc.Gall_Height_mm)
ggplot(gall_selection.df, aes(x = sc.Density_per_100_shoots, y = gall_survival, color = Treatment.focus)) +
  geom_point() + binomial_smooth() +
  facet_wrap(~cut(sc.Gall_Height_mm, breaks = c(-5,0,5)))

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

viridis6 <- function() {
  # colours taken from the viridis package
  c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725")
}

control.gamm <- gamm4(gall_survival ~ s(sc.Gall_Height_mm) + s(sc.gall_individuals,k=9) + s(sc.Density_per_100_shoots),
                      random = ~(1|Genotype/Plant_Position/Gall_Number/Gall_ID),
                      data = filter(gall_selection.df, Treatment.focus == "Control"),
                      family = "binomial")#,
#method = "GCV.Cp")
summary(control.gamm$gam)
summary(control.gamm$mer)
gam.check(control.gamm$gam)

plot(control.gamm$gam, seWithMean = T, residuals = T)

plot(control.gamm$gam, seWithMean = T, shift = mean(predict(control.gamm$gam)), trans = function(x) {exp(x)/(1+exp(x))})

gam.gradients(control.gamm$gam, phenotype = c("sc.Gall_Height_mm","sc.gall_individuals"), covariates = "sc.Density_per_100_shoots", se.method = 'n')
gam.gradients(control.gamm$gam, phenotype = c("sc.Density_per_100_shoots","sc.gall_individuals"), covariates = "sc.Gall_Height_mm", se.method = 'n')
gam.gradients(control.gamm$gam, phenotype = c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals", se.method = 'n')


control.fl <- fitness.landscape(control.gamm$gam, phenotype = c("sc.Gall_Height_mm","sc.gall_individuals"), covariates = "sc.Density_per_100_shoots", PI.method = 'n')
control.fl_df <- data.frame(control.fl$points, Wbar = control.fl$Wbar)
ggplot(control.fl_df, aes(x = Var1, y = Var2, fill = Wbar)) + geom_raster() + scale_fill_gradientn(colors = viridis6())

control.fl.2 <- fitness.landscape(control.gamm$gam, phenotype = c("sc.Density_per_100_shoots","sc.gall_individuals"), covariates = "sc.Gall_Height_mm", PI.method = 'n')
control.fl.2_df <- data.frame(control.fl.2$points, Wbar = control.fl.2$Wbar)
ggplot(control.fl.2_df, aes(x = Var1, y = Var2, fill = Wbar)) + geom_raster() + scale_fill_gradientn(colors = viridis6())

control.fl.3 <- fitness.landscape(control.gamm$gam, phenotype = c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals", PI.method = 'n')
control.fl.3_df <- data.frame(control.fl.3$points, Wbar = control.fl.3$Wbar)
ggplot(control.fl.3_df, aes(x = Var1, y = Var2, fill = Wbar)) + geom_raster() + scale_fill_gradientn(colors = viridis6())


treatment.gamm <- gamm(gall_survival ~ s(sc.Gall_Height_mm) + s(sc.gall_individuals, k = 9) + s(sc.Density_per_100_shoots),
                       random = list(Genotype=~1, Plant_Position=~1, Gall_Number=~1),
                       data = filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion"),
                       family = "binomial",
                       method = "GCV.Cp")
summary(treatment.gamm$gam)
summary(treatment.gamm$lme)
#plot(treatment.gamm$gam, shift = mean(predict(treatment.gamm$gam)), trans = function(x) {exp(x)/(1+exp(x))})

gam.gradients(treatment.gamm$gam, phenotype = c("sc.Gall_Height_mm","sc.gall_individuals"), covariates = "sc.Density_per_100_shoots", se.method = 'n')
gam.gradients(treatment.gamm$gam, phenotype = c("sc.Density_per_100_shoots","sc.gall_individuals"), covariates = "sc.Gall_Height_mm", se.method = 'n')
gam.gradients(treatment.gamm$gam, phenotype = c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals", se.method = 'n')

treatment.fl <- fitness.landscape(treatment.gamm$gam, phenotype = c("sc.Gall_Height_mm","sc.gall_individuals"), covariates = "sc.Density_per_100_shoots", PI.method = 'n')
treatment.fl_df <- data.frame(treatment.fl$points, Wbar = treatment.fl$Wbar)
ggplot(treatment.fl_df, aes(x = Var1, y = Var2, fill = Wbar)) + geom_raster() + scale_fill_gradientn(colors = viridis6())

treatment.fl.2 <- fitness.landscape(treatment.gamm$gam, phenotype = c("sc.Density_per_100_shoots","sc.gall_individuals"), covariates = "sc.Gall_Height_mm", PI.method = 'n')
treatment.fl.2_df <- data.frame(treatment.fl.2$points, Wbar = treatment.fl.2$Wbar)
ggplot(treatment.fl.2_df, aes(x = Var1, y = Var2, fill = Wbar)) + geom_raster() + scale_fill_gradientn(colors = viridis6())

treatment.fl.3 <- fitness.landscape(treatment.gamm$gam, phenotype = c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals", PI.method = 'n')
treatment.fl.3_df <- data.frame(treatment.fl.3$points, Wbar = treatment.fl.3$Wbar)
ggplot(treatment.fl.3_df, aes(x = Var1, y = Var2, fill = Wbar)) + geom_raster() + scale_fill_gradientn(colors = viridis6())

## STEEPNESS OF FITNESS LANDSCAPE
sd(control.fl_df$Wbar)/mean(control.fl_df$Wbar) # 0.26
sd(control.fl.2_df$Wbar)/mean(control.fl.2_df$Wbar) # 0.01
sd(control.fl.3_df$Wbar)/mean(control.fl.3_df$Wbar) # 0.26
mean(c(0.26,0.01,0.26))

sd(treatment.fl_df$Wbar)/mean(treatment.fl_df$Wbar) # 0.18
sd(treatment.fl.2_df$Wbar)/mean(treatment.fl.2_df$Wbar) # 0.14
sd(treatment.fl.3_df$Wbar)/mean(treatment.fl.3_df$Wbar) # 0.21
mean(c(0.18,0.14,0.21))

## OLD BUT MAYBE USEFUL ----

## SUMMARISE DATA AT GALL & PLANT LEVEL ----
mean.narm <- function(x) mean(x, na.rm = TRUE)

gall_level.info <- gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position, Gall_Number) %>%
  summarise_at(vars(Density_per_100_shoots, gall_individuals, Gall_Height_mm), mean.narm) %>% 
  ungroup()

plant_level.info <- gall_level.info %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots, gall_individuals, Gall_Height_mm), mean.narm) %>%
  ungroup()

gall_level.parasitism <- gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position, Gall_Number) %>%
  summarise_at(vars(pupa, platy, ectos, platy.ectos, total), sum) %>%
  ungroup()

plant_level.parasitism <- gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(pupa, platy, ectos, platy.ectos, total), sum) %>%
  ungroup()

gall_level.df <- left_join(gall_level.info, gall_level.parasitism) %>%
  mutate(gall_survival = pupa/total,
         log.indiv = log(gall_individuals),
         sqrt.size = sqrt(Gall_Height_mm),
         log1.density = log(Density_per_100_shoots+1))

plant_level.df <- left_join(plant_level.info, plant_level.parasitism) %>%
  mutate(gall_survival = pupa/total,
         sc.Gall_Height_mm = scale(Gall_Height_mm),
         sc.gall_individuals = scale(gall_individuals),
         sc.Density_per_100_shoots = scale(Density_per_100_shoots))

plant.cont.glmer <- glm(gall_survival ~ Treatment.focus*(sc.Gall_Height_mm + sc.gall_individuals + sc.Density_per_100_shoots)^2,# + (1|Genotype/Plant_Position/Gall_Number),
                        data = gall_selection.df,
                        family = "binomial")#, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))#,
#weights = total)
summary(plant.cont.glmer)

plant.cont.gam <- gam(gall_survival ~ s(sc.Gall_Height_mm), by = sc.gall_individuals)+ s(Density_per_100_shoots),
#random = ~(1|Genotype),
data = filter(plant_level.df, Treatment.focus == "Ectoparasitoid exclusion"),
family = "binomial",
weights = filter(plant_level.df, Treatment.focus == "Ectoparasitoid exclusion")$total)
summary(plant.cont.gam$gam)
plot(plant.cont.gam$gam, seWithMean = T, shift = mean(predict(plant.cont.gam$gam)), trans = function(x) {exp(x)/(1+exp(x))})
vis.gam(plant.cont.gam$gam, view = c("sc.Gall_Height_mm","sc.gall_individuals"), type = "response", plot.type = "contour")

# wow, still a ton of variability in gall size at the polygall level
gall.size.df <- left_join(select(gall_selection.df, Treatment.focus, Genotype, Plant_Position, gall.size = Gall_Height_mm), 
                          select(gall_level.df, Treatment.focus, Genotype, Plant_Position, polygall.size = Gall_Height_mm)) 
ggplot(gall.size.df, aes(x = polygall.size, y = gall.size)) + 
  geom_point() + geom_smooth(method = "lm")   
summary(lm(gall.size ~ polygall.size, data = gall.size.df)) # only explains about 10% of the variance

## FULL MODEL
control.df <- as.data.frame(filter(gall_selection.df, Treatment.focus == "Control")) %>%
  mutate(Gall_Height_mm = scale(Gall_Height_mm),
         gall_individuals = scale(gall_individuals),
         Density_per_100_shoots = scale(Density_per_100_shoots))
gppr.test <- gppr(y = "gall_survival", xterms = c("Gall_Height_mm","gall_individuals","Density_per_100_shoots"), 
                  data = control.df, nterms = 2)
gppr.test$ppr$alpha
gppr.test$ppr$beta
plot(gppr.test$ppr, ask = T)

control.df$term1 <- control.df$Gall_Height_mm*gppr.test$ppr$alpha[1] + control.df$gall_individuals*gppr.test$ppr$alpha[2] + control.df$Density_per_100_shoots*gppr.test$ppr$alpha[3]
control.df$term2 <- control.df$Gall_Height_mm*gppr.test$ppr$alpha[4] + control.df$gall_individuals*gppr.test$ppr$alpha[5] + control.df$Density_per_100_shoots*gppr.test$ppr$alpha[6]


test.gam <- gamm(gall_survival ~ s(term2), 
                 random = list(Genotype=~1, Plant_Position=~1, Gall_Number=~1),
                 data = control.df, 
                 family = "binomial")
summary(test.gam$gam)
plot(test.gam$gam, shift = mean(predict(test.gam$gam)),
     trans = function(x) {exp(x)/(1+exp(x))})

term1.control.gradient <- gam.gradients(test.gam$gam, phenotype = "term1")
term1.control.gradient$ests

library(lme4)
cont.glmer <- glmer(gall_survival ~ (Gall_Height_mm + gall_individuals + Density_per_100_shoots)^2 + (1|Genotype/Plant_Position/Gall_Number),
                    data = control.df,
                    family = "binomial",
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))
summary(cont.glmer)

treat.glmer <- glmer(gall_survival ~ Gall_Height_mm*gall_individuals*Density_per_100_shoots + (1|Genotype/Plant_Position/Gall_Number),
                     data = treatment.df,
                     family = "binomial",
                     control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))
summary(treat.glmer)

all.glmer <- glmer(gall_survival ~ Treatment.focus*(scale(Gall_Height_mm) + scale(log(gall_individuals)) + scale(sqrt(Density_per_100_shoots)))^2 + (1|Genotype/Plant_Position/Gall_Number),
                   data = gall_selection.df,
                   family = "binomial",
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))
summary(all.glmer)

# why does the plant level data differ so much???
plant.all.glmer <- glm(gall_survival ~ Treatment.focus*(scale(Gall_Height_mm) + scale(gall_individuals) + scale(Density_per_100_shoots))^2,# + (1|Genotype),
                       data = plant_level.df,
                       family = "binomial",
                       weights = total)#,
# control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))
summary(plant.all.glmer)

plant.gam <-  gam(gall_survival ~ s(scale(Gall_Height_mm)) + s(scale(gall_individuals)) + s(scale(Density_per_100_shoots)),# + (1|Genotype),
                  data = filter(plant_level.df, Treatment.focus == "Control"),
                  family = "binomial",
                  weights = total,
                  method = "GCV.Cp")
summary(plant.gam)
plot(plant.gam, shift = mean(predict(plant.gam)),
     trans = function(x) {exp(x)/(1+exp(x))})

cont.gam <- gam(gall_survival ~ s(Gall_Height_mm) + s(Gall_Height_mm, by = gall_individuals) + s(gall_individuals, k = 9) + s(Density_per_100_shoots), 
                #random = list(Genotype=~1, Plant_Position=~1, Gall_Number=~1),
                data = control.df, 
                family = "binomial")
summary(cont.gam)
plot(cont.gam, shift = mean(predict(cont.gam)),
     trans = function(x) {exp(x)/(1+exp(x))})

cont.gam.grad <- gam.gradients(cont.gam$gam, phenotype = "Gall_Height_mm", covariates = c("gall_individuals","Density_per_100_shoots"))
cont.gam.grad$ests

# why does this work when standardized = F but not T? Maybe I should standardize prior to going into the model.
gradients.test <- gppr.gradients(gppr.test, phenotype = c("Gall_Height_mm", "gall_individuals"), covariates = "Density_per_100_shoots", se.method = 'n', standardized = F)
gradients.test

treatment.df <- as.data.frame(filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion")) %>%
  mutate(Gall_Height_mm = scale(Gall_Height_mm),
         gall_individuals = scale(gall_individuals),
         Density_per_100_shoots = scale(Density_per_100_shoots))
gppr.treat.test <- gppr(y = "gall_survival", xterms = c("Gall_Height_mm","gall_individuals","Density_per_100_shoots"), 
                        data = treatment.df, nterms = 3)
gppr.treat.test$ppr$alpha
gppr.treat.test$ppr$beta
plot(gppr.treat.test$ppr, ask = T)

treatment.df$term1 <- treatment.df$Gall_Height_mm*gppr.treat.test$ppr$alpha[1] + treatment.df$gall_individuals*gppr.treat.test$ppr$alpha[2] + treatment.df$Density_per_100_shoots*gppr.treat.test$ppr$alpha[3]

treat.test <- gamm(gall_survival ~ s(term1), 
                   random = list(Genotype=~1, Plant_Position=~1, Gall_Number=~1),
                   data = treatment.df, 
                   family = "binomial")
summary(treat.test$gam)
plot(treat.test$gam, shift = mean(predict(treat.test$gam)),
     trans = function(x) {exp(x)/(1+exp(x))})
term1.gradient <- gam.gradients(treat.test$gam, phenotype = "term1")
term1.gradient$ests

treat.gam <- gamm(gall_survival ~ s(Gall_Height_mm) + s(gall_individuals, k = 9) + s(Density_per_100_shoots), 
                  random = list(Genotype=~1, Plant_Position=~1, Gall_Number=~1),
                  data = treatment.df, 
                  family = "binomial")
summary(treat.gam)#$gam)
plot(treat.gam, shift = mean(predict(treat.gam)),
     trans = function(x) {exp(x)/(1+exp(x))})

treat.gam.grad <- gam.gradients(treat.gam, phenotype = c("Gall_Height_mm","Density_per_100_shoots"), covariates = "gall_individuals")
treat.gam.grad$ests

vis.gam(treat.gam$gam, view = c("Gall_Height_mm","Density_per_100_shoots"), type = "response", plot.type = "contour")

data("SoayLambs")
soay.gppr <- gppr(y = "W", xterms = c("WEIGHT","HINDLEG","HORNLEN","lnKeds"), 
                  data = SoayLambs, nterms = 2)
soay.gppr$ppr$alpha
soay.gppr$ppr$beta

SoayLambs$term1 <- SoayLambs$WEIGHT*soay.gppr$ppr$alpha[1] + SoayLambs$HINDLEG*soay.gppr$ppr$alpha[2] + SoayLambs$HORNLEN*soay.gppr$ppr$alpha[3] + SoayLambs$lnKeds*soay.gppr$ppr$alpha[4] 
major.gam <- gam(W ~ s(term1), SoayLambs, family = "binomial")
summary(major.gam)
gam.gradients(major.gam, phenotype = "term1", se.method = 'n')
plot(major.gam, se = T, seWithMean = T, rug = F, shift = mean(predict(major.gam)),
     trans = function(x) {exp(x)/(1+exp(x))})

test <- gppr.gradients(soay.gppr, phenotype = c("HINDLEG","WEIGHT"), covariates = c("HORNLEN","lnKeds"), se.method = 'n')
test
# CONTROL
z.gam <- gamm(gall_survival ~ (Gall_Height_mm + gall_individuals + Density_per_100_shoots)^2,
              random = list(Genotype=~1, Plant_Position=~1, Gall_Number=~1),
              data = filter(gall_selection.df, Treatment.focus == "Control"), 
              family = binomial(link = "logit"), 
              method = "GCV.Cp")
summary(z.gam$gam)
plot(z.gam$gam, se = T, seWithMean = T, rug = F, shift = mean(predict(z.gam$gam)),
     trans = function(x) {exp(x)/(1+exp(x))})
# not informative # vis.gam(z.gam$gam, view = c("gall_individuals", "Gall_Height_mm"), type = "response", plot.type = "contour")
# not informative # vis.gam(z.gam$gam, view = c("Density_per_100_shoots", "Gall_Height_mm"), type = "response", plot.type = "contour")
# not informative # vis.gam(z.gam$gam, view = c("Density_per_100_shoots","gall_individuals"), type = "response", plot.type = "contour")

gradient.size_indiv <- gam.gradients(z.gam$gam, phenotype = c("Gall_Height_mm", "gall_individuals"), covariates = "Density_per_100_shoots", standardized = T)
gradient.size_indiv$ests

gradient.size_density <- gam.gradients(z.gam$gam, phenotype = c("Gall_Height_mm", "Density_per_100_shoots"), covariates = "gall_individuals", standardized = T)
gradient.size_density$ests

gradient.indiv_density <- gam.gradients(z.gam$gam, phenotype = c("gall_individuals", "Density_per_100_shoots"), covariates = "Gall_Height_mm", standardized = T)
gradient.indiv_density$ests

# TREATMENT
z.gam.treat <- gamm(gall_survival ~ Gall_Height_mm + gall_individuals + Density_per_100_shoots,
                    random = list(Genotype=~1, Plant_Position=~1, Gall_Number=~1),
                    data = filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion"), 
                    family = binomial(link = "logit"), 
                    method = "GCV.Cp")
summary(z.gam.treat$gam)
plot(z.gam.treat$gam, se = T, seWithMean = T, rug = F, shift = mean(predict(z.gam.treat$gam)),
     trans = function(x) {exp(x)/(1+exp(x))})
vis.gam(z.gam.treat$gam, view = c("gall_individuals", "Gall_Height_mm"), type = "response", plot.type = "contour")
vis.gam(z.gam.treat$gam, view = c("Density_per_100_shoots", "Gall_Height_mm"), type = "response", plot.type = "contour")
vis.gam(z.gam.treat$gam, view = c("Density_per_100_shoots","gall_individuals"), type = "response", plot.type = "contour")

gradient.size_indiv <- gam.gradients(z.gam.treat$gam, 
                                     phenotype = c("Gall_Height_mm", "gall_individuals"), 
                                     covariates = "Density_per_100_shoots", 
                                     standardized = T)
gradient.size_indiv$ests

gradient.size_density <- gam.gradients(z.gam.treat$gam, phenotype = c("Gall_Height_mm", "Density_per_100_shoots"), covariates = "gall_individuals", standardized = T)
gradient.size_density$ests

gradient.indiv_density <- gam.gradients(z.gam.treat$gam, phenotype = c("gall_individuals", "Density_per_100_shoots"), covariates = "Gall_Height_mm", standardized = T)
gradient.indiv_density$ests

## SELECTION ON GALL SIZE
# treating each individual gall chamber as an independent data point

z.size <- gam(gall_survival ~ s(Gall_Height_mm), 
              data = filter(gall_selection.df, Treatment.focus == "Control"), 
              family = binomial(link = "logit"), 
              method = "GCV.Cp")
summary(z.size)
plot(z.size, se = T, seWithMean = T, rug = F, shift = mean(predict(z.size)),
     trans = function(x) {exp(x)/(1+exp(x))})
z.size.gradient <- gam.gradients(z.size, phenotype = "Gall_Height_mm", standardized = T)
z.size.gradient$ests


z.size.treat <- gam(gall_survival ~ s(Gall_Height_mm), 
                    data = filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion"),
                    family = binomial(link = "logit"), 
                    method = "GCV.Cp")
summary(z.size.treat)
plot(z.size.treat, se = T, seWithMean = T, rug = F, shift = mean(predict(z.size.treat)),
     trans = function(x) {exp(x)/(1+exp(x))})
z.size.treat.gradient <- gam.gradients(z.size.treat, phenotype = "Gall_Height_mm", standardized = T)
z.size.treat.gradient$ests

# treating each poly-gall as an independent data point
dim(filter(gall_level.df, Treatment.focus == "Control"))[1] # 336 poly galls
z.size.poly <- gam(gall_survival ~ s(Gall_Height_mm), 
                   data = filter(gall_level.df, Treatment.focus == "Control"), 
                   family = binomial(link = "logit"), 
                   method = "GCV.Cp",
                   weight = total)
summary(z.size.poly)
plot(z.size.poly, se = T, seWithMean = T, rug = F, shift = mean(predict(z.size.poly)),
     trans = function(x) {exp(x)/(1+exp(x))})
z.size.poly.gradient <- gam.gradients(z.size.poly, phenotype = "Gall_Height_mm", standardized = T)
z.size.poly.gradient$ests

dim(filter(gall_level.df, Treatment.focus == "Ectoparasitoid exclusion"))[1] # 278 poly galls
z.size.poly.treat <- gam(gall_survival ~ s(Gall_Height_mm), 
                         data = filter(gall_level.df, Treatment.focus == "Ectoparasitoid exclusion"),
                         family = binomial(link = "logit"), 
                         method = "GCV.Cp",
                         weight = total)
summary(z.size.poly.treat)
plot(z.size.poly.treat, se = T, seWithMean = T, rug = F, shift = mean(predict(z.size.poly.treat)),
     trans = function(x) {exp(x)/(1+exp(x))})
z.size.poly.treat.gradient <- gam.gradients(z.size.poly.treat, phenotype = "Gall_Height_mm", standardized = T)
z.size.poly.treat.gradient$ests

## SELECTION ON GALL INDIVIDUALS
ggplot(gall_level.df, aes(x = gall_individuals, y = Gall_Height_mm)) + 
  geom_point() + geom_smooth(method = "lm", formula = y ~ x) + 
  facet_wrap(~Treatment.focus) + scale_x_log10() + scale_y_sqrt()

cor.test(log(gall_level.df$gall_individuals), gall_level.df$Gall_Height_mm)

test <- gamm(gall_survival ~ s(Gall_Height_mm) + s(log.indiv, k = 9),#+ s(log1.density), 
             random = list(Genotype=~1, Plant_Position=~1),
             data = filter(gall_level.df, Treatment.focus == "Control"), 
             family = binomial(link = "logit"), 
             method = "GCV.Cp",
             weights = total)
plot(test$gam, se = T, seWithMean = T, rug = F, shift = mean(predict(test$gam)),
     trans = function(x) {exp(x)/(1+exp(x))})
gam.gradients(test$gam, phenotype = c("Gall_Height_mm", "log.indiv"), standardized = T)

test2 <- gamm(gall_survival ~ s(Gall_Height_mm) + s(log.indiv, k = 9) + s(log1.density), 
              random = list(Genotype=~1, Plant_Position=~1),
              data = filter(gall_level.df, Treatment.focus == "Ectoparasitoid exclusion"), 
              family = binomial(link = "logit"), 
              method = "GCV.Cp",
              weights = total)
plot(test2$gam, se = T, seWithMean = T, rug = F, shift = mean(predict(test2$gam)),
     trans = function(x) {exp(x)/(1+exp(x))})
summary(test2$gam)
gam.gradients(test2$gam, phenotype = c("Gall_Height_mm", "log.indiv"), standardized = T)

# treating each poly-gall as an independent data point

z.indiv <- gam(gall_survival ~ s(sqrt.size) + s(log.indiv, k = 9), # set to 9, highest possible value while still running
               data = filter(gall_level.df, Treatment.focus == "Control"), 
               family = binomial(link = "logit"), 
               method = "GCV.Cp",
               weights = total)
summary(z.indiv)
plot(z.indiv, se = T, seWithMean = T, rug = F, shift = mean(predict(z.indiv)),
     trans = function(x) {exp(x)/(1+exp(x))})
z.indiv.gradient <- gam.gradients(z.indiv, phenotype = c("sqrt.size", "log.indiv"), standardized = T)
z.indiv.gradient$ests

vis.gam(z.indiv, plot.type = 'contour', type = 'response')

z.indiv.landscape <- fitness.landscape(z.indiv, phenotype = c("sqrt.size", "log.indiv"), parallel = 'multicore', ncpus = 32)

ggplot(data.frame(z.indiv.landscape$points, Wbar = z.indiv.landscape$Wbar), aes(x = Var1, y = Var2, fill = Wbar)) +
  geom_raster() #geom_point(shape = 21, size = 5)


# trouble fitting gam.gradient with k = 9
z.indiv.treat <- gam(gall_survival ~ s(sqrt.size) + s(log.indiv, k = 9), 
                     data = filter(gall_level.df, Treatment.focus == "Ectoparasitoid exclusion"), 
                     family = binomial(link = "logit"), 
                     method = "GCV.Cp",
                     weights = total)
summary(z.indiv.treat)
plot(z.indiv.treat, se = T, seWithMean = T, rug = F, shift = mean(predict(z.indiv.treat)),
     trans = function(x) {exp(x)/(1+exp(x))})
z.indiv.treat.gradient <- gam.gradients(z.indiv.treat, phenotype = c("Gall_Height_mm", "log.indiv"), standardized = T)
z.indiv.treat.gradient$ests

vis.gam(z.indiv.treat, plot.type = 'contour', type = 'response')

z.indiv.treat.landscape <- fitness.landscape(z.indiv.treat, phenotype = c("Gall_Height_mm", "log.indiv"), parallel = 'multicore', ncpus = 32)

ggplot(data.frame(z.indiv.treat.landscape$points, Wbar = z.indiv.treat.landscape$Wbar), aes(x = Var1, y = Var2, fill = Wbar)) +
  geom_raster()#geom_point(shape = 21, size = 5)

## SELECTION ON GALL DENSITY
# treating each tree as an independent data point

dim(filter(plant_level.df, Treatment.focus == "Control"))[1] # 56 plants
z.density <- gam(gall_survival ~ s(log1.density), 
                 data = filter(plant_level.df, Treatment.focus == "Control"), 
                 family = binomial(link = "logit"), 
                 method = "GCV.Cp", 
                 weights = total)
summary(z.density)
plot(z.density, se = T, seWithMean = T, rug = F, shift = mean(predict(z.density)),
     trans = function(x) {exp(x)/(1+exp(x))})
z.density.gradient <- gam.gradients(z.density, phenotype = "log1.density", standardized = T)
z.density.gradient$ests

dim(filter(plant_level.df, Treatment.focus == "Ectoparasitoid exclusion"))[1] # 56 plants
z.density.treat <- gam(gall_survival ~ s(log1.density), 
                       data = filter(plant_level.df, Treatment.focus == "Ectoparasitoid exclusion"), 
                       family = binomial(link = "logit"), 
                       method = "GCV.Cp",
                       weights = total)
summary(z.density.treat)
plot(z.density.treat, se = T, seWithMean = T, rug = F, shift = mean(predict(z.density.treat)),
     trans = function(x) {exp(x)/(1+exp(x))})
z.density.treat.gradient <- gam.gradients(z.density.treat, phenotype = "log1.density", standardized = T)
z.density.treat.gradient$ests
