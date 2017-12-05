
## LOAD REQUIRED LIBRARIES & SET OPTIONS ----

library(tidyverse)
library(mgcv)
library(gamm4) # for generalized additive mixed models #
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
  mutate(gall_survival = as.numeric(ifelse(pupa > 0, 1, 0)),
         egg_parasitoid = as.numeric(ifelse(platy > 0, 1, 0)),
         sc.Gall_Height_mm = as.numeric(scale(sqrt(Gall_Height_mm))),
         sc.gall_individuals = as.numeric(scale(log(gall_individuals))),
         sc.Density_per_100_shoots = as.numeric(scale(sqrt(Density_per_100_shoots))))

## GENERALIZED ADDITIVE MIXED MODELS + GENERALIZED PROJECTION PURSUIT REGRESSION ----

gamm.model <- function(formula, data) {
  gamm4(formula, data, 
        random = ~(1|Genotype/Plant_Position/Gall_Number),
        family = binomial(link = logit))
}

gamm.plot <- function(gamm.model) {
  plot(gamm.model$gam, seWithMean = T, shift = mean(predict(gamm.model$gam)), trans = function(x) {exp(x)/(1+exp(x))})
}

# experimenting with boot.case gave bootstrapped distributions that weren't centered near the model mean estimate.
# also, using 'refit.smooth = T' also generated some very skewed distributions
# the recipe below provided the most sensible estimates for my data
gradient.calc <- function(mod, phenotype, covariates = NULL) {
  gam.gradients(mod, phenotype, covariates, parallel = 'multicore', ncpus = 32, 
                refit.smooth = F, 
                se.method = 'boot.para',
                standardized = T)
} 

insect_FL <- function(mod, phenotype, covariates = NULL) {
  fitness.landscape(mod, phenotype, covariates, parallel = 'multicore', ncpus = 32,
                    PI.method = 'boot.para',
                    refit.smooth = F,
                    PI.interval = c(0.025,0.975),
                    plt.density = 100)
}

get_insect_FL.df <- function(FL){
  data.frame(mean_fitness = FL$Wbar, FL$points, t(FL$WbarPI)) %>%
    mutate(lower_2.5 = X2.5., upper_97.5 = X97.5.)
}

ggplot_uniFL <- function(FL.df) {
  ggplot(FL.df, aes(x = Var1, y = mean_fitness)) +
    geom_line() + 
    geom_ribbon(aes(ymin = lower_2.5, ymax = upper_97.5), alpha = 0.5)
}

## GALLS ON CONTROL TREES ----

# subset data
control.df <- as.data.frame(filter(gall_selection.df, Treatment.focus == "Control"))

# get major axis of selection
control.gppr <- gppr(y = "gall_survival", 
                     xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), 
                     data = control.df, nterms = 1)
control.gppr$ppr$alpha  
control.gppr$ppr$beta
control.df$term1 <- control.gppr$ppr$alpha[1]*control.df$sc.Gall_Height_mm + control.gppr$ppr$alpha[2]*control.df$sc.gall_individuals + control.gppr$ppr$alpha[3]*control.df$sc.Density_per_100_shoots
# control.df$term2 <- control.gppr$ppr$alpha[4]*control.df$sc.Gall_Height_mm + control.gppr$ppr$alpha[5]*control.df$sc.gall_individuals + control.gppr$ppr$alpha[6]*control.df$sc.Density_per_100_shoots
# plot(term2 ~ term1, control.df)
# cor.test(control.df$term1, control.df$term2) 

# GAMM
control.major.gam <- gamm.model(gall_survival ~ s(term1), data = control.df)
summary(control.major.gam$gam)
summary(control.major.gam$mer)
gamm.plot(control.major.gam)

# selection gradient for major axis of selection
control.major.gradients <- gradient.calc(mod = control.major.gam$gam, phenotype = c("term1"))
control.major.gradients$ests

# check distributions of bootstrapped estimates
hist(control.major.gradients$boot[,1])
abline(v = control.major.gradients$ests[1,1], col = "red", lty = 2)

hist(control.major.gradients$boot[,2])
abline(v = control.major.gradients$ests[2,1], col = "red", lty = 2)

# GAMM components
control.component.gam <- gamm.model(gall_survival ~ s(sc.Gall_Height_mm) + s(sc.gall_individuals, k = 9) + s(sc.Density_per_100_shoots), data = control.df)
summary(control.component.gam$gam)
summary(control.component.gam$mer)
concurvity(control.component.gam$gam)
gamm.plot(control.component.gam)

# selection gradients
control.height_indiv.gam <- gradient.calc(mod = control.component.gam$gam, phenotype = c("sc.Gall_Height_mm","sc.gall_individuals"), covariates = "sc.Density_per_100_shoots")
control.height_indiv.gam$ests

control.height_density.gam <- gradient.calc(mod = control.component.gam$gam, phenotype = c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals")
control.height_density.gam$ests

control.density_indiv.gam <- gradient.calc(mod = control.component.gam$gam, phenotype = c("sc.gall_individuals","sc.Density_per_100_shoots"), covariates = "sc.Gall_Height_mm")
control.density_indiv.gam$ests

# fitness landscapes

# plots to examine phenotypic distributions. Important for fitness landscapes to only plot where I have most of my data.
ggplot(control.df, aes(x = sc.Gall_Height_mm, y = sc.Density_per_100_shoots)) +
  geom_point(shape = 1)
ggplot(control.df, aes(x = sc.Gall_Height_mm, y = sc.gall_individuals)) +
  geom_jitter(shape = 1, height = 0.05)
ggplot(control.df, aes(x = sc.Density_per_100_shoots, y = sc.gall_individuals)) +
  geom_jitter(shape = 1, height = 0.05)


## GALLS ON ECTOPARSITOID EXCLUSION TREES ----

# subset data
treatment.df <- as.data.frame(filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion"))

# get major axis of selection
treatment.gppr <- gppr(y = "gall_survival", xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), 
                     data = treatment.df, nterms = 1)
treatment.gppr$ppr$alpha # strong loadings of all three traits, with gall height acting in opposing directions to gall individuals and gall density
treatment.gppr$ppr$beta
treatment.df$term1 <- treatment.gppr$ppr$alpha[1]*treatment.df$sc.Gall_Height_mm + treatment.gppr$ppr$alpha[2]*treatment.df$sc.gall_individuals + treatment.gppr$ppr$alpha[3]*treatment.df$sc.Density_per_100_shoots
# treatment.df$term2 <- treatment.gppr$ppr$alpha[4]*treatment.df$sc.Gall_Height_mm + treatment.gppr$ppr$alpha[5]*treatment.df$sc.gall_individuals + treatment.gppr$ppr$alpha[6]*treatment.df$sc.Density_per_100_shoots
# plot(term2 ~ term1, treatment.df)
# cor.test(treatment.df$term1, treatment.df$term2) 

# GAMM
treatment.major.gam <- gamm.model(gall_survival ~ s(term1), data = treatment.df)
summary(treatment.major.gam$gam)
summary(treatment.major.gam$mer)
gamm.plot(treatment.major.gam)

# selection gradient for major axis of selection
treatment.major.gradients <- gradient.calc(mod = treatment.major.gam$gam, phenotype = "term1")
treatment.major.gradients$ests

# check distributions of bootstrapped estimates
hist(treatment.major.gradients$boot[,1])
abline(v = treatment.major.gradients$ests[1,1], col = "red", lty = 2)

hist(treatment.major.gradients$boot[,2])
abline(v = treatment.major.gradients$ests[2,1], col = "red", lty = 2)

# GAMM components
treatment.component.gam <- gamm.model(gall_survival ~ s(sc.Gall_Height_mm) + s(sc.gall_individuals) + s(sc.Density_per_100_shoots), data = treatment.df)
summary(treatment.component.gam$gam)
summary(treatment.component.gam$mer)
concurvity(treatment.component.gam$gam)
gamm.plot(treatment.component.gam)
vis.gam(treatment.component.gam$gam, view = c("sc.Gall_Height_mm","sc.gall_individuals"), type = "response", plot.type = "contour")
vis.gam(treatment.component.gam$gam, view = c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), type = "response", plot.type = "contour")
vis.gam(treatment.component.gam$gam, view = c("sc.gall_individuals","sc.Density_per_100_shoots"), type = "response", plot.type = "contour")


# selection gradients
treatment.height_indiv.gradients <- gradient.calc(mod=treatment.component.gam$gam, phenotype=c("sc.Gall_Height_mm","sc.gall_individuals"), covariates = "sc.Density_per_100_shoots")
treatment.height_indiv.gradients$ests

treatment.height_density.gradients <- gradient.calc(mod=treatment.component.gam$gam, phenotype=c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals")
treatment.height_density.gradients$ests

treatment.density_indiv.gradients <- gradient.calc(mod=treatment.component.gam$gam, phenotype=c("sc.gall_individuals","sc.Density_per_100_shoots"), covariates = "sc.Gall_Height_mm")
treatment.density_indiv.gradients$ests

treatment.height.FL <- insect_FL(mod=treatment.component.gam$gam, phenotype=c("sc.Gall_Height_mm"), covariates = c("sc.gall_individuals","sc.Density_per_100_shoots"))
treatment.height.FL.df <- get_insect_FL.df(treatment.height.FL)
ggplot_uniFL(treatment.height.FL.df)

treatment.density.FL <- insect_FL(mod=treatment.component.gam$gam, phenotype=c("sc.Density_per_100_shoots"), covariates = c("sc.gall_individuals","sc.Gall_Height_mm"))
treatment.density.FL.df <- get_insect_FL.df(treatment.density.FL)
ggplot_uniFL(treatment.density.FL.df) 

treatment.indiv.FL <- insect_FL(mod=treatment.component.gam$gam, phenotype="sc.gall_individuals", covariates = c("sc.Gall_Height_mm","sc.Density_per_100_shoots"))
treatment.indiv.FL.df <- get_insect_FL.df(treatment.indiv.FL)
ggplot_uniFL(treatment.indiv.FL.df)

treatment.all.FL <- fitness.landscape(mod=treatment.all.gam,phenotype=c("sc.Gall_Height_mm","sc.gall_individuals"), covariates = "sc.Density_per_100_shoots", parallel = 'multicore', ncpus = 32)
treatment.all.FL.2 <- fitness.landscape(mod=treatment.all.gam,phenotype=c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals", parallel = 'multicore', ncpus = 32)
treatment.all.FL.3 <- fitness.landscape(mod=treatment.all.gam,phenotype=c("sc.gall_individuals","sc.Density_per_100_shoots"), covariates = "sc.Gall_Height_mm", parallel = 'multicore', ncpus = 32)

viridis6 <- function() {
  # colours taken from the viridis package
  c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725")
}

ggplot(treatment.df, aes(x = sc.Gall_Height_mm, y = sc.Density_per_100_shoots)) +
  geom_point(shape = 1)
ggplot(treatment.df, aes(x = sc.Gall_Height_mm, y = sc.gall_individuals)) +
  geom_jitter(shape = 1, height = 0.05)
ggplot(treatment.df, aes(x = sc.Density_per_100_shoots, y = sc.gall_individuals)) +
  geom_jitter(shape = 1, height = 0.05)

treatment.all.FL_df <- data.frame(treatment.all.FL$points, Wbar = treatment.all.FL$Wbar) %>%
  rename(sc.Gall_Height_mm = Var1, sc.gall_individuals = Var2)
ggplot(treatment.all.FL_df, aes(x = sc.Gall_Height_mm, y = sc.gall_individuals, fill = Wbar)) + geom_raster() + 
  scale_fill_gradientn(colors = viridis6()) +
  geom_jitter(data = filter(treatment.df, 
                           sc.Gall_Height_mm < 1 & sc.Gall_Height_mm > -1 & sc.gall_individuals < 1 & sc.gall_individuals > -1), 
             aes(x = sc.Gall_Height_mm, y = sc.gall_individuals), inherit.aes = F, height = 0.05, shape = 1)

treatment.all.FL.2_df <- data.frame(treatment.all.FL.2$points, Wbar = treatment.all.FL.2$Wbar) %>%
  rename(sc.Gall_Height_mm = Var1, sc.Density_per_100_shoots = Var2)
ggplot(treatment.all.FL.2_df, aes(x = sc.Gall_Height_mm, y = sc.Density_per_100_shoots, fill = Wbar)) + geom_raster() + 
  scale_fill_gradientn(colors = viridis6()) +
  geom_point(data = filter(treatment.df, 
                            sc.Gall_Height_mm < 1 & sc.Gall_Height_mm > -1 & sc.Density_per_100_shoots < 1 & sc.Density_per_100_shoots > -1), 
              aes(x = sc.Gall_Height_mm, y = sc.Density_per_100_shoots), inherit.aes = F, shape = 1)

treatment.all.FL.3_df <- data.frame(treatment.all.FL.3$points, Wbar = treatment.all.FL.3$Wbar) %>%
  rename(sc.gall_individuals = Var1, sc.Density_per_100_shoots = Var2)
ggplot(treatment.all.FL.3_df, aes(x = sc.gall_individuals, y = sc.Density_per_100_shoots, fill = Wbar)) + geom_raster() + 
  scale_fill_gradientn(colors = viridis6()) + 
  geom_jitter(data = filter(treatment.df, 
                           sc.gall_individuals < 1 & sc.gall_individuals > -1 & sc.Density_per_100_shoots < 1 & sc.Density_per_100_shoots > -1), 
             aes(x = sc.gall_individuals, y = sc.Density_per_100_shoots), width = 0.05, inherit.aes = F, shape = 1)

## ECTOPARASITOIDS ON CONTROL TREES
control.ecto.gppr <- gppr(y = "ectos", xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), data = control.df, nterms = 1)
control.ecto.gppr$ppr$alpha
control.ecto.gppr$ppr$beta
control.df$term1_ectos <- control.ecto.gppr$ppr$alpha[1]*control.df$sc.Gall_Height_mm + control.ecto.gppr$ppr$alpha[2]*control.df$sc.gall_individuals + control.ecto.gppr$ppr$alpha[3]*control.df$sc.Density_per_100_shoots

control.ecto.major.gam <- gamm.model(ectos ~ s(term1_ectos), data = control.df)
summary(control.ecto.major.gam$gam)
summary(control.ecto.major.gam$mer)
gamm.plot(control.ecto.major.gam)

control.ecto.component.gam <- gamm.model(ectos ~ s(sc.Gall_Height_mm) + s(sc.gall_individuals, k=9) + s(sc.Density_per_100_shoots), data = control.df)
summary(control.ecto.component.gam$gam)
summary(control.ecto.component.gam$mer)
concurvity(control.ecto.component.gam$gam)
gamm.plot(control.ecto.component.gam)


## PLATYGASTER ON CONTROL TREES
control.platy.gppr <- gppr(y = "egg_parasitoid", xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), data = control.df, nterms = 1) 
control.platy.gppr$ppr$alpha
control.platy.gppr$ppr$beta
control.df$term1_platys <- control.platy.gppr$ppr$alpha[1]*control.df$sc.Gall_Height_mm + control.platy.gppr$ppr$alpha[2]*control.df$sc.gall_individuals + control.platy.gppr$ppr$alpha[3]*control.df$sc.Density_per_100_shoots
# control.df$term2_platys <- control.platy.gppr$ppr$alpha[4]*control.df$sc.Gall_Height_mm + control.platy.gppr$ppr$alpha[5]*control.df$sc.gall_individuals + control.platy.gppr$ppr$alpha[6]*control.df$sc.Density_per_100_shoots
# plot(term2_platys ~ term1_platys, control.df)
# cor.test(control.df$term1_platys, control.df$term2_platys) # cor = -0.69

control.platy.major.gam <- gamm.model(egg_parasitoid ~ s(term1_platys), data = control.df)
summary(control.platy.major.gam$gam)
summary(control.platy.major.gam$mer)
gamm.plot(control.platy.major.gam)

control.platy.component.gam <- gamm.model(egg_parasitoid ~ s(sc.Gall_Height_mm) + s(sc.gall_individuals, k=9) + s(sc.Density_per_100_shoots), data = control.df)
summary(control.platy.component.gam$gam)
summary(control.platy.component.gam$mer)
gamm.plot(control.platy.component.gam)


## PLATYGASTER ON ECTOPARASITOID EXCLUSION TREES
treatment.platy.gppr <- gppr(y = "egg_parasitoid", xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), data = treatment.df, nterms = 1)
treatment.platy.gppr$ppr$alpha
treatment.platy.gppr$ppr$beta
treatment.df$term1_platys <- treatment.platy.gppr$ppr$alpha[1]*treatment.df$sc.Gall_Height_mm + treatment.platy.gppr$ppr$alpha[2]*treatment.df$sc.gall_individuals + treatment.platy.gppr$ppr$alpha[3]*treatment.df$sc.Density_per_100_shoots

treatment.platy.major.gam <- gamm.model(egg_parasitoid ~ s(term1_platys), data = treatment.df)
summary(treatment.platy.major.gam$gam)
summary(treatment.platy.major.gam$mer)
gamm.plot(treatment.platy.major.gam)

treatment.platy.component.gam <- gamm.model(egg_parasitoid ~ s(sc.Gall_Height_mm) + s(sc.gall_individuals) + s(sc.Density_per_100_shoots), data = treatment.df)
summary(treatment.platy.component.gam$gam)
summary(treatment.platy.component.gam$mer)
gamm.plot(treatment.platy.component.gam)
