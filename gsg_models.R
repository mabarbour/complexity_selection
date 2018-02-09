
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

## GENERALIZED ADDITIVE MIXED MODELS + GENERALIZED PROJECTION PURSUIT REGRESSION ----

gamm.model <- function(formula, data) {
  gamm4(formula, data, 
        random = ~(1|Genotype/Plant_Position/Gall_Number),
        family = binomial(link = logit))
}

gamm.plot <- function(gamm.model, ...) {
  plot(gamm.model$gam, seWithMean = T, shift = mean(predict(gamm.model$gam)), trans = function(x) {exp(x)/(1+exp(x))}, ...)
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

## RAW DATA SUMMARIES ----
sum(filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion")$pupa) # 372
sum(filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion")$platy) # 230
sum(filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion")$ectos) # 17
sum(filter(gall_selection.df, Treatment.focus == "Control")$ectos) # 168
sum(filter(gall_selection.df, Treatment.focus == "Control")$platy) # 217
sum(filter(gall_selection.df, Treatment.focus == "Control")$pupa) # 323

## FULL MODEL ----
#full.component.gam <- gamm.model(gall_survival ~ 
#                                   Treatment.focus +
#                                   s(sc.Gall_Height_mm, by=Treatment.focus) + 
#                                   s(sc.gall_individuals, by=Treatment.focus, k=9) + 
#                                   s(sc.Density_per_100_shoots, by=Treatment.focus),
#                                 data = gall_selection.df)


#summary(full.component.gam$gam)
#summary(full.component.gam$mer)
#concurvity(full.component.gam$gam)
#gamm.plot(full.component.gam, pages=1)

## I don't know if this will work for interaction variables...maybe need to subset somehow...
# selection gradients
#full.height_indiv.gam <- gradient.calc(mod = full.component.gam$gam, phenotype = c("sc.Gall_Height_mm"), covariates = c("sc.gall_individuals","sc.Density_per_100_shoots","Treatment.focus"))
#full.height_indiv.gam$ests

#full.height_density.gam <- gradient.calc(mod = full.component.gam$gam, phenotype = c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals")
#full.height_density.gam$ests

#full.density_indiv.gam <- gradient.calc(mod = full.component.gam$gam, phenotype = c("sc.gall_individuals","sc.Density_per_100_shoots"), covariates = "sc.Gall_Height_mm")
#full.density_indiv.gam$ests

## DATA DISTRIBUTIONS
ggplot(distinct(gall_selection.df, Plant_Position, Density_per_100_shoots, Treatment.focus),
       aes(x = Density_per_100_shoots, fill = Treatment.focus)) +
  geom_density(alpha=0.5)
ggplot(distinct(gall_selection.df, Gall_Number, gall_individuals, Treatment.focus),
       aes(x = gall_individuals, fill = Treatment.focus)) +
  geom_density(alpha=0.5)
ggplot(gall_selection.df,
       aes(x = Gall_Height_mm, fill = Treatment.focus)) +
  geom_density(alpha=0.5)

pc_galls <- princomp(select(gall_selection.df, sc.Gall_Height_mm, sc.gall_individuals, sc.Density_per_100_shoots))
biplot(pc_galls)
pc_galls$scores

gall_selection.df <- data.frame(gall_selection.df, pc_galls$scores)

## GALLS ON CONTROL TREES ----

# subset data
control.df <- as.data.frame(filter(gall_selection.df, Treatment.focus == "Control"))

# get major axis of selection
#control.gppr <- gppr(y = "gall_survival", 
#                     xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), 
#                     data = control.df, nterms = 1)
#control.gppr$ppr$alpha  
#control.gppr$ppr$beta
#control.df$term1 <- control.gppr$ppr$alpha[1]*control.df$sc.Gall_Height_mm + control.gppr$ppr$alpha[2]*control.df$sc.gall_individuals + control.gppr$ppr$alpha[3]*control.df$sc.Density_per_100_shoots
# control.df$term2 <- control.gppr$ppr$alpha[4]*control.df$sc.Gall_Height_mm + control.gppr$ppr$alpha[5]*control.df$sc.gall_individuals + control.gppr$ppr$alpha[6]*control.df$sc.Density_per_100_shoots
# plot(term2 ~ term1, control.df)
# cor.test(control.df$term1, control.df$term2) 

control.pc.gam <- gamm.model(gall_survival ~ s(Comp.1) + s(Comp.2) + s(Comp.3), data = control.df)
summary(control.pc.gam$gam)
summary(control.pc.gam$mer)
concurvity(control.pc.gam$gam)
gamm.plot(control.pc.gam, pages=1)

vis.gam(control.pc.gam$gam, view = c("Comp.1","Comp.2"), type = "response", color="topo", theta=45, ticktype="detailed")
vis.gam(control.pc.gam$gam, view = c("Comp.1","Comp.3"), type = "response", color="topo", theta=45, ticktype="detailed")


treatment.pc.gam <- gamm.model(gall_survival ~ s(Comp.1) + s(Comp.2) + s(Comp.3), data = treatment.df)
summary(treatment.pc.gam$gam)
summary(treatment.pc.gam$mer)
concurvity(treatment.pc.gam$gam)
gamm.plot(treatment.pc.gam, pages=1)

vis.gam(treatment.pc.gam$gam, view = c("Comp.1","Comp.2"), type = "response", color="topo", theta=45, ticktype="detailed")
vis.gam(treatment.pc.gam$gam, view = c("Comp.1","Comp.3"), type = "response", color="topo", theta=45, ticktype="detailed")
vis.gam(treatment.pc.gam$gam, view = c("Comp.1","Comp.2"), type = "response", plot.type="contour")
vis.gam(treatment.pc.gam$gam, view = c("Comp.1","Comp.3"), type = "response", plot.type="contour", color="topo")


library(visreg)
devtools::install_github("pbreheny/visreg")
visreg(treatment.pc.gam$gam, xvar="Comp.1",  scale="response")

# GAMM
#control.major.gam <- gamm.model(gall_survival ~ s(term1), data = control.df)
#summary(control.major.gam$gam)
#summary(control.major.gam$mer)
#gamm.plot(control.major.gam)

# selection gradient for major axis of selection
#control.major.gradients <- gradient.calc(mod = control.major.gam$gam, phenotype = c("term1"))
#control.major.gradients$ests

# check distributions of bootstrapped estimates
#hist(control.major.gradients$boot[,1])
#abline(v = control.major.gradients$ests[1,1], col = "red", lty = 2)

#hist(control.major.gradients$boot[,2])
#abline(v = control.major.gradients$ests[2,1], col = "red", lty = 2)

# GAMM components
control.component.gam <- gamm.model(gall_survival ~ s(Gall_Height_mm) + s(gall_individuals, k = 9) + s(Density_per_100_shoots), data = control.df)
summary(control.component.gam$gam)
summary(control.component.gam$mer)
concurvity(control.component.gam$gam)
gamm.plot(control.component.gam, pages=1)

# selection gradients
control.height_indiv.gam <- gradient.calc(mod = control.component.gam$gam, phenotype = c("Gall_Height_mm","gall_individuals"), covariates = "Density_per_100_shoots")
control.height_indiv.gam$ests

control.height_density.gam <- gradient.calc(mod = control.component.gam$gam, phenotype = c("Gall_Height_mm","Density_per_100_shoots"), covariates = "gall_individuals")
control.height_density.gam$ests

control.density_indiv.gam <- gradient.calc(mod = control.component.gam$gam, phenotype = c("gall_individuals","Density_per_100_shoots"), covariates = "Gall_Height_mm")
control.density_indiv.gam$ests


# fitness landscapes

# plots to examine phenotypic distributions. Important for fitness landscapes to only plot where I have most of my data.
#ggplot(control.df, aes(x = sc.Gall_Height_mm, y = sc.Density_per_100_shoots)) +
#  geom_point(shape = 1)
#ggplot(control.df, aes(x = sc.Gall_Height_mm, y = sc.gall_individuals)) +
#  geom_jitter(shape = 1, height = 0.05)
#ggplot(control.df, aes(x = sc.Density_per_100_shoots, y = sc.gall_individuals)) +
#  geom_jitter(shape = 1, height = 0.05)
#pcs <- princomp(x = select(gall_selection.df, Gall_Height_mm, gall_individuals, Density_per_100_shoots), cor=TRUE)
#summary(pcs)
#loadings(pcs)
#pcs$scores
#pc.all <- data.frame(gall_selection.df, pcs$scores)
#pc.all$Treatment.focus <- C(pc.all$Treatment.focus, "contr.sum")

#library(lme4)
#library(lmerTest)
#comp.1 <- lmer(Comp.1 ~ Treatment.focus + (1|Genotype/Plant_Position/Gall_Number), pc.all)
#summary(comp.1)
#comp.2 <- lmer(Comp.2 ~ Treatment.focus + (1|Genotype/Plant_Position/Gall_Number), pc.all)
#summary(comp.2)

#pc.all.component.gam <- gamm.model(gall_survival ~ Treatment.focus + s(Comp.1, by=Treatment.focus) + s(Comp.2, by=Treatment.focus) + s(Comp.3, by=Treatment.focus), data = pc.all)
#summary(pc.all.component.gam$gam)
#summary(pc.all.component.gam$mer)
#concurvity(pc.all.component.gam$gam)
#gamm.plot(pc.all.component.gam, pages = 1)
#vis.gam(pc.all.component.gam$gam, view = c("Comp.1","Comp.2"),cond = list(Treatment.focus = "Control"), type = "response", plot.type = "contour")
#vis.gam(pc.all.component.gam$gam, view = c("Comp.1","Comp.2"),cond = list(Treatment.focus = "Ectoparasitoid exclusion"), type = "response", plot.type = "contour")
#vis.gam(pc.all.component.gam$gam, view = c("Comp.2","Comp.3"),cond = list(Treatment.focus = "Control"), type = "response", plot.type = "contour")
#vis.gam(pc.all.component.gam$gam, view = c("Comp.2","Comp.3"),cond = list(Treatment.focus = "Ectoparasitoid exclusion"), type = "response", plot.type = "contour")
#loadings(pcs)


#pcs <- princomp(x = select(control.df, Gall_Height_mm, gall_individuals, Density_per_100_shoots), cor=TRUE)
#summary(pcs)
#loadings(pcs)
#pcs$scores
#pc.control <- data.frame(control.df, pcs$scores)

#pc.control.component.gam <- gamm.model(gall_survival ~ s(Comp.1) + s(Comp.2), data = pc.control)
#summary(pc.control.component.gam$gam)
#summary(pc.control.component.gam$mer)
#concurvity(pc.control.component.gam$gam)
#gamm.plot(pc.control.component.gam, pages = 1)
#vis.gam(pc.control.component.gam$gam, view = c("Comp.1","Comp.2"), type = "response", plot.type = "contour")
#loadings(pcs)


#control.density <- control.df %>%
#  group_by(Plant_Position) %>%
#  summarise(total_galls = n(), gall_survival = sum(gall_survival)/total_galls,  Density_per_100_shoots = mean(Density_per_100_shoots)) %>%
#  ungroup()
#control.density.gam <- gam(gall_survival ~ s(Density_per_100_shoots), data = control.density, weights = control.density$total_galls, family = binomial(link = logit))
#summary(control.density.gam)
#plot(control.density.gam, seWithMean=T, shift = mean(predict(control.density.gam)), trans = function(x) {exp(x)/(1+exp(x))})
#gradient.calc(mod=control.density.gam, phenotype="Density_per_100_shoots")

#control.indiv <- control.df %>%
#  group_by(Plant_Position, Gall_Number) %>%
#  summarise(total_galls = n(), gall_survival = sum(gall_survival)/total_galls,  gall_individuals = mean(gall_individuals)) %>%
#  ungroup()
#control.indiv.gam <- gamm4(gall_survival ~ s(gall_individuals, k=9), random=~(1|Plant_Position), data = control.indiv, weights = control.indiv$total_galls, family = binomial(link = logit))
#summary(control.indiv.gam$gam)
#plot(control.indiv.gam$gam, seWithMean=T, shift = mean(predict(control.indiv.gam$gam)), trans = function(x) {exp(x)/(1+exp(x))})
#gradient.calc(mod=control.indiv.gam$gam, phenotype="gall_individuals")

#control.size.gam <- gamm4(gall_survival ~ s(Gall_Height_mm), random=~(1|Plant_Position/Gall_Number), data = control.df, family = binomial(link = logit))
#summary(control.size.gam$gam)
#plot(control.size.gam$gam, seWithMean=T, shift = mean(predict(control.size.gam$gam)), trans = function(x) {exp(x)/(1+exp(x))})
#gradient.calc(mod=control.size.gam$gam, phenotype="Gall_Height_mm")

## GALLS ON ECTOPARSITOID EXCLUSION TREES ----

# subset data
treatment.df <- as.data.frame(filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion"))

# get major axis of selection
#treatment.gppr <- gppr(y = "gall_survival", xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), 
#                     data = treatment.df, nterms = 1)
#treatment.gppr$ppr$alpha # strong loadings of all three traits, with gall height acting in opposing directions to gall individuals and gall density
#treatment.gppr$ppr$beta
#treatment.df$term1 <- treatment.gppr$ppr$alpha[1]*treatment.df$sc.Gall_Height_mm + treatment.gppr$ppr$alpha[2]*treatment.df$sc.gall_individuals + treatment.gppr$ppr$alpha[3]*treatment.df$sc.Density_per_100_shoots
# treatment.df$term2 <- treatment.gppr$ppr$alpha[4]*treatment.df$sc.Gall_Height_mm + treatment.gppr$ppr$alpha[5]*treatment.df$sc.gall_individuals + treatment.gppr$ppr$alpha[6]*treatment.df$sc.Density_per_100_shoots
# plot(term2 ~ term1, treatment.df)
# cor.test(treatment.df$term1, treatment.df$term2) 

# GAMM
#treatment.major.gam <- gamm.model(gall_survival ~ s(term1), data = treatment.df)
#summary(treatment.major.gam$gam)
#summary(treatment.major.gam$mer)
#gamm.plot(treatment.major.gam)

# selection gradient for major axis of selection
#treatment.major.gradients <- gradient.calc(mod = treatment.major.gam$gam, phenotype = "term1")
#treatment.major.gradients$ests

# check distributions of bootstrapped estimates
#hist(treatment.major.gradients$boot[,1])
#abline(v = treatment.major.gradients$ests[1,1], col = "red", lty = 2)

#hist(treatment.major.gradients$boot[,2])
#abline(v = treatment.major.gradients$ests[2,1], col = "red", lty = 2)

#treatment.density <- treatment.df %>%
#  group_by(Genotype, Plant_Position) %>%
#  summarise(total_galls = n(), gall_survival = sum(gall_survival)/total_galls,  Density_per_100_shoots = mean(Density_per_100_shoots)) %>%
#  ungroup()
#treatment.density.gam <- gam(gall_survival ~ s(Density_per_100_shoots), random = ~(1|Genotype), data = treatment.density, weights = treatment.density$total_galls, family = binomial(link = logit))
#summary(treatment.density.gam)
#plot(treatment.density.gam, seWithMean=T, shift = mean(predict(treatment.density.gam)), trans = function(x) {exp(x)/(1+exp(x))})
#gradient.calc(mod=treatment.density.gam, phenotype="Density_per_100_shoots")$ests

#treatment.indiv <- treatment.df %>%
#  group_by(Genotype, Plant_Position, Gall_Number) %>%
#  summarise(total_galls = n(), gall_survival = sum(gall_survival)/total_galls,  gall_individuals = mean(gall_individuals)) %>%
#  ungroup()
#treatment.indiv.gam <- gamm4(gall_survival ~ s(gall_individuals), random=~(1|Genotype/Plant_Position), data = treatment.indiv, weights = treatment.indiv$total_galls, family = binomial(link = logit))
#summary(treatment.indiv.gam$gam)
#plot(treatment.indiv.gam$gam, seWithMean=T, shift = mean(predict(treatment.indiv.gam$gam)), trans = function(x) {exp(x)/(1+exp(x))})
#gradient.calc(mod=treatment.indiv.gam$gam, phenotype="gall_individuals")$ests

#treatment.size.gam <- gamm4(gall_survival ~ s(Gall_Height_mm), random=~(1|Genotype/Plant_Position/Gall_Number), data = treatment.df, family = binomial(link = logit))
#summary(treatment.size.gam$gam)
#plot(treatment.size.gam$gam, seWithMean=T, shift = mean(predict(treatment.size.gam$gam)), trans = function(x) {exp(x)/(1+exp(x))})
#gradient.calc(mod=treatment.size.gam$gam, phenotype="Gall_Height_mm")$ests

# GAMM components
#pcs <- princomp(x = select(treatment.df, Gall_Height_mm, gall_individuals, Density_per_100_shoots), cor=TRUE)
#summary(pcs)
#loadings(pcs)
#pcs$scores
#pc.treatment <- data.frame(treatment.df, pcs$scores)

#pc.treatment.component.gam <- gamm.model(gall_survival ~ s(Comp.1) + s(Comp.2), data = pc.treatment)
#summary(pc.treatment.component.gam$gam)
#summary(pc.treatment.component.gam$mer)
#concurvity(pc.treatment.component.gam$gam)
#gamm.plot(pc.treatment.component.gam, pages = 1)
#vis.gam(pc.treatment.component.gam$gam, view = c("Comp.1","Comp.2"), type = "response", plot.type = "contour")
#loadings(pcs)

treatment.component.gam <- gamm.model(gall_survival ~ s(Gall_Height_mm) + s(gall_individuals) + s(Density_per_100_shoots), data = treatment.df)
summary(treatment.component.gam$gam)
summary(treatment.component.gam$mer)
concurvity(treatment.component.gam$gam)
gamm.plot(treatment.component.gam, pages = 1)

#vis.gam(treatment.component.gam$gam, view = c("Gall_Height_mm","gall_individuals"), type = "response", plot.type = "contour")
#vis.gam(treatment.component.gam$gam, view = c("Gall_Height_mm","Density_per_100_shoots"), type = "response", plot.type = "contour")
#vis.gam(treatment.component.gam$gam, view = c("gall_individuals","Density_per_100_shoots"), type = "response", plot.type = "contour")


# selection gradients
treatment.height_indiv.gradients <- gradient.calc(mod=treatment.component.gam$gam, phenotype=c("Gall_Height_mm","gall_individuals"), covariates = "Density_per_100_shoots")
treatment.height_indiv.gradients$ests

treatment.height_density.gradients <- gradient.calc(mod=treatment.component.gam$gam, phenotype=c("Gall_Height_mm","Density_per_100_shoots"), covariates = "gall_individuals")
treatment.height_density.gradients$ests

treatment.density_indiv.gradients <- gradient.calc(mod=treatment.component.gam$gam, phenotype=c("gall_individuals","Density_per_100_shoots"), covariates = "Gall_Height_mm")
treatment.density_indiv.gradients$ests

treatment.height.FL <- insect_FL(mod=treatment.component.gam$gam, phenotype=c("Gall_Height_mm"), covariates = c("gall_individuals","Density_per_100_shoots"))
treatment.height.FL.df <- get_insect_FL.df(treatment.height.FL)
ggplot_uniFL(treatment.height.FL.df)

treatment.density.FL <- insect_FL(mod=treatment.component.gam$gam, phenotype=c("Density_per_100_shoots"), covariates = c("gall_individuals","Gall_Height_mm"))
treatment.density.FL.df <- get_insect_FL.df(treatment.density.FL)
ggplot_uniFL(treatment.density.FL.df) 

treatment.indiv.FL <- insect_FL(mod=treatment.component.gam$gam, phenotype="gall_individuals", covariates = c("Gall_Height_mm","Density_per_100_shoots"))
treatment.indiv.FL.df <- get_insect_FL.df(treatment.indiv.FL)
ggplot_uniFL(treatment.indiv.FL.df)

#treatment.all.FL <- fitness.landscape(mod=treatment.all.gam,phenotype=c("sc.Gall_Height_mm","sc.gall_individuals"), covariates = "sc.Density_per_100_shoots", parallel = 'multicore', ncpus = 32)
#treatment.all.FL.2 <- fitness.landscape(mod=treatment.all.gam,phenotype=c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals", parallel = 'multicore', ncpus = 32)
#treatment.all.FL.3 <- fitness.landscape(mod=treatment.all.gam,phenotype=c("sc.gall_individuals","sc.Density_per_100_shoots"), covariates = "sc.Gall_Height_mm", parallel = 'multicore', ncpus = 32)

viridis6 <- function() {
  # colours taken from the viridis package
  c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725")
}

#ggplot(treatment.df, aes(x = sc.Gall_Height_mm, y = sc.Density_per_100_shoots)) +
#  geom_point(shape = 1)
#ggplot(treatment.df, aes(x = sc.Gall_Height_mm, y = sc.gall_individuals)) +
#  geom_jitter(shape = 1, height = 0.05)
#ggplot(treatment.df, aes(x = sc.Density_per_100_shoots, y = sc.gall_individuals)) +
#  geom_jitter(shape = 1, height = 0.05)

#treatment.all.FL_df <- data.frame(treatment.all.FL$points, Wbar = treatment.all.FL$Wbar) %>%
#  rename(sc.Gall_Height_mm = Var1, sc.gall_individuals = Var2)
#ggplot(treatment.all.FL_df, aes(x = sc.Gall_Height_mm, y = sc.gall_individuals, fill = Wbar)) + geom_raster() + 
#  scale_fill_gradientn(colors = viridis6()) +
#  geom_jitter(data = filter(treatment.df, 
#                           sc.Gall_Height_mm < 1 & sc.Gall_Height_mm > -1 & sc.gall_individuals < 1 & sc.gall_individuals > -1), 
#             aes(x = sc.Gall_Height_mm, y = sc.gall_individuals), inherit.aes = F, height = 0.05, shape = 1)

#treatment.all.FL.2_df <- data.frame(treatment.all.FL.2$points, Wbar = treatment.all.FL.2$Wbar) %>%
#  rename(sc.Gall_Height_mm = Var1, sc.Density_per_100_shoots = Var2)
#ggplot(treatment.all.FL.2_df, aes(x = sc.Gall_Height_mm, y = sc.Density_per_100_shoots, fill = Wbar)) + geom_raster() + 
#  scale_fill_gradientn(colors = viridis6()) +
#  geom_point(data = filter(treatment.df, 
#                            sc.Gall_Height_mm < 1 & sc.Gall_Height_mm > -1 & sc.Density_per_100_shoots < 1 & sc.Density_per_100_shoots > -1), 
#              aes(x = sc.Gall_Height_mm, y = sc.Density_per_100_shoots), inherit.aes = F, shape = 1)

#treatment.all.FL.3_df <- data.frame(treatment.all.FL.3$points, Wbar = treatment.all.FL.3$Wbar) %>%
#  rename(sc.gall_individuals = Var1, sc.Density_per_100_shoots = Var2)
#ggplot(treatment.all.FL.3_df, aes(x = sc.gall_individuals, y = sc.Density_per_100_shoots, fill = Wbar)) + geom_raster() + 
#  scale_fill_gradientn(colors = viridis6()) + 
#  geom_jitter(data = filter(treatment.df, 
#                           sc.gall_individuals < 1 & sc.gall_individuals > -1 & sc.Density_per_100_shoots < 1 & sc.Density_per_100_shoots > -1), 
#             aes(x = sc.gall_individuals, y = sc.Density_per_100_shoots), width = 0.05, inherit.aes = F, shape = 1)

## PLOTS COMPARING CONTROL AND EXCLUSION TREES ----
mean_Gall_Height_mm <- mean(gall_selection.df$Gall_Height_mm)
sd_Gall_Height_mm <- sd(gall_selection.df$Gall_Height_mm)
range_Gall_Height_mm <- range(gall_selection.df$Gall_Height_mm)

mean_gall_individuals <- mean(gall_selection.df$gall_individuals)
sd_gall_individuals <- sd(gall_selection.df$gall_individuals)
range_gall_individuals <- range(gall_selection.df$gall_individuals)

mean_Density_per_100_shoots <- mean(gall_selection.df$Density_per_100_shoots)
sd_Density_per_100_shoots <- sd(gall_selection.df$Density_per_100_shoots)
range_Density_per_100_shoots <- range(gall_selection.df$Density_per_100_shoots)

# control predictions
newdata.height.control <- data.frame(Gall_Height_mm = seq(range_Gall_Height_mm[1], range_Gall_Height_mm[2], length.out=100),
                                     gall_individuals = mean_gall_individuals,
                                     Density_per_100_shoots = mean_Density_per_100_shoots)
predict.height.control <- predict(control.component.gam$gam, newdata = newdata.height.control, type = "response", se.fit = TRUE)
predict.height.control.df <- cbind.data.frame(predict.height.control, newdata.height.control)

newdata.indiv.control <- data.frame(Gall_Height_mm = mean_Gall_Height_mm,
                                     gall_individuals = seq(range_gall_individuals[1], range_gall_individuals[2], length.out=100),
                                     Density_per_100_shoots = mean_Density_per_100_shoots)
predict.indiv.control <- predict(control.component.gam$gam, newdata = newdata.indiv.control, type = "response", se.fit = TRUE)
predict.indiv.control.df <- cbind.data.frame(predict.indiv.control, newdata.indiv.control)

newdata.density.control <- data.frame(Gall_Height_mm = mean_Gall_Height_mm,
                                     gall_individuals = mean_gall_individuals,
                                     Density_per_100_shoots = seq(range_Density_per_100_shoots[1], range_Density_per_100_shoots[2], length.out=100))
predict.density.control <- predict(control.component.gam$gam, newdata = newdata.density.control, type = "response", se.fit = TRUE)
predict.density.control.df <- cbind.data.frame(predict.density.control, newdata.density.control)

# exclusion predictions
newdata.height.treatment <- data.frame(Gall_Height_mm = seq(range_Gall_Height_mm[1], range_Gall_Height_mm[2], length.out=100),
                                     gall_individuals = mean_gall_individuals,
                                     Density_per_100_shoots = mean_Density_per_100_shoots)
predict.height.treatment <- predict(treatment.component.gam$gam, newdata = newdata.height.treatment, type = "response", se.fit = TRUE)
predict.height.treatment.df <- cbind.data.frame(predict.height.treatment, newdata.height.treatment)

newdata.indiv.treatment <- data.frame(Gall_Height_mm = mean_Gall_Height_mm,
                                    gall_individuals = seq(range_gall_individuals[1], range_gall_individuals[2], length.out=100),
                                    Density_per_100_shoots = mean_Density_per_100_shoots)
predict.indiv.treatment <- predict(treatment.component.gam$gam, newdata = newdata.indiv.treatment, type = "response", se.fit = TRUE)
predict.indiv.treatment.df <- cbind.data.frame(predict.indiv.treatment, newdata.indiv.treatment)

newdata.density.treatment <- data.frame(Gall_Height_mm = mean_Gall_Height_mm,
                                      gall_individuals = mean_gall_individuals,
                                      Density_per_100_shoots = seq(range_Density_per_100_shoots[1], range_Density_per_100_shoots[2], length.out=100))
predict.density.treatment <- predict(treatment.component.gam$gam, newdata = newdata.density.treatment, type = "response", se.fit = TRUE)
predict.density.treatment.df <- cbind.data.frame(predict.density.treatment, newdata.density.treatment)

# combine data
predict.height.df <-  bind_rows(mutate(predict.height.control.df, treatment = "Complex"), 
                                mutate(predict.height.treatment.df, treatment = "Simple")) 
predict.indiv.df <-  bind_rows(mutate(predict.indiv.control.df, treatment = "Complex"), 
                                mutate(predict.indiv.treatment.df, treatment = "Simple")) 
predict.density.df <-  bind_rows(mutate(predict.density.control.df, treatment = "Complex"), 
                                mutate(predict.density.treatment.df, treatment = "Simple")) 
predict.all.df <- bind_rows(predict.height.df, predict.indiv.df, predict.density.df)

# plots
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

fitness.height <- ggplot(predict.height.df, aes(x = Gall_Height_mm, y = fit, color = treatment)) + 
  geom_line(size = 1.5) + 
  geom_line(aes(y = fit + se.fit), linetype = "dotted") +
  geom_line(aes(y = fit - se.fit), linetype = "dotted") +
  ylab("") + xlab("Gall diameter (mm)") + scale_color_manual(values = cbPalette[c(2,3)], guide = "none") +
  scale_y_continuous(limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) + 
  scale_x_continuous(limits=c(mean_Gall_Height_mm-sd_Gall_Height_mm, mean_Gall_Height_mm+sd_Gall_Height_mm))

fitness.indiv <- ggplot(predict.indiv.df, aes(x = gall_individuals, y = fit, color = treatment)) + 
  geom_line(size = 1.5) + 
  geom_line(aes(y = fit + se.fit), linetype = "dotted") +
  geom_line(aes(y = fit - se.fit), linetype = "dotted") +
  ylab("") + xlab("Larva per gall") + scale_color_manual(values = cbPalette[c(2,3)], guide = "none") +
  scale_y_continuous(limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) + 
  scale_x_continuous(limits=c(mean_gall_individuals-sd_gall_individuals, mean_gall_individuals+sd_gall_individuals))

fitness.density <- ggplot(predict.density.df, aes(x = Density_per_100_shoots, y = fit, color = treatment)) + 
  geom_line(size = 1.5) + 
  geom_line(aes(y = fit + se.fit), linetype = "dotted") +
  geom_line(aes(y = fit - se.fit), linetype = "dotted") +
  ylab("") + xlab("Larva per 100 shoots") + scale_color_manual(values = cbPalette[c(2,3)], guide = "none") +
  scale_y_continuous(limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) + 
  scale_x_continuous(limits=c(mean_Density_per_100_shoots-sd_Density_per_100_shoots, mean_Density_per_100_shoots+sd_Density_per_100_shoots))

fitness.plot_subGalls <- plot_grid(fitness.height, fitness.indiv, fitness.density, nrow=1)
save_plot("selection_gradients_subGalls.pdf", fitness.plot_subGalls, base_width = 11, base_height = 8.5)

# plot based on principal components
gall_pcs <- princomp(select(gall_selection.df, Gall_Height_mm, gall_individuals, Density_per_100_shoots), cor=TRUE)
biplot(gall_pcs)

gall_loadings <- as.data.frame(gall_pcs$loadings[ ,1:2])
scores_Comp.1 <- gall_loadings$Comp.1[1]*scale(gall_selection.df$Gall_Height_mm) + gall_loadings$Comp.1[2]*scale(gall_selection.df$gall_individuals) + gall_loadings$Comp.1[3]*scale(gall_selection.df$Density_per_100_shoots)
scores_Comp.2 <- gall_loadings$Comp.2[1]*scale(gall_selection.df$Gall_Height_mm) + gall_loadings$Comp.2[2]*scale(gall_selection.df$gall_individuals) + gall_loadings$Comp.2[3]*scale(gall_selection.df$Density_per_100_shoots)

sc.gall_size <- scale(gall_selection.df$Gall_Height_mm)
sc.gall_size * attr(sc.gall_size,"scaled:scale") + attr(sc.gall_size,"scaled:center") - gall_selection.df$Gall_Height_mm # confirms that this is the correct way to transform between scaled and original form

sc.indiv <- scale(gall_selection.df$gall_individuals)
sc.density <- scale(gall_selection.df$Density_per_100_shoots)


gall_trait_combos <- expand.grid(seq((mean_Gall_Height_mm-sd_Gall_Height_mm), (mean_Gall_Height_mm+sd_Gall_Height_mm), length.out = 50), 
            seq((mean_gall_individuals-sd_gall_individuals), (mean_gall_individuals+sd_gall_individuals), length.out = 50),
            seq((mean_Density_per_100_shoots-sd_Density_per_100_shoots), (mean_Density_per_100_shoots+sd_Density_per_100_shoots), length.out = 50))
colnames(gall_trait_combos) <- c("Gall_Height_mm","gall_individuals","Density_per_100_shoots")

predict.control <- predict(control.component.gam$gam, newdata = gall_trait_combos, type = "response", se.fit = TRUE)
predict.control.df <- cbind.data.frame(predict.control, gall_trait_combos)

predict.treatment <- predict(treatment.component.gam$gam, newdata = gall_trait_combos, type = "response", se.fit = TRUE)
predict.treatment.df <- cbind.data.frame(predict.treatment, gall_trait_combos)

predict.combo.df <- bind_rows(mutate(predict.control.df, treatment = "Complex", standard.fit = fit/mean(fit)-1),
                              mutate(predict.treatment.df, treatment = "Simple", standard.fit = fit/mean(fit)-1))

predict.combo.df$Comp.1 <- gall_loadings$Comp.1[1]*(predict.combo.df$Gall_Height_mm - attr(sc.gall_size,"scaled:center"))/attr(sc.gall_size,"scaled:scale") + 
  gall_loadings$Comp.1[2]*(predict.combo.df$gall_individuals - attr(sc.indiv,"scaled:center"))/attr(sc.indiv,"scaled:scale") +
  gall_loadings$Comp.1[3]*(predict.combo.df$Density_per_100_shoots - attr(sc.density,"scaled:center"))/attr(sc.density,"scaled:scale") 

predict.combo.df$Comp.2 <- gall_loadings$Comp.2[1]*(predict.combo.df$Gall_Height_mm - attr(sc.gall_size,"scaled:center"))/attr(sc.gall_size,"scaled:scale") + 
  gall_loadings$Comp.2[2]*(predict.combo.df$gall_individuals - attr(sc.indiv,"scaled:center"))/attr(sc.indiv,"scaled:scale") +
  gall_loadings$Comp.2[3]*(predict.combo.df$Density_per_100_shoots - attr(sc.density,"scaled:center"))/attr(sc.density,"scaled:scale") 

ggplot(predict.combo.df, aes(x=Comp.1, y=Comp.2, z=standard.fit, color=standard.fit)) + 
  geom_point(size=10) +
  scale_color_gradientn(colors = viridis6()) +
  facet_wrap(~treatment)

ggplot(predict.combo.df, aes(x=Gall_Height_mm, y=gall_individuals, z=standard.fit, color=standard.fit)) + 
  geom_point(size=10) +
  scale_color_gradientn(colors = viridis6()) +
  facet_wrap(~treatment)

ggplot(predict.combo.df, aes(x=gall_individuals, y=Density_per_100_shoots, z=standard.fit, color=standard.fit)) + 
  geom_point(size=10) +
  scale_color_gradientn(colors = viridis6()) +
  facet_wrap(~treatment)

ggplot(predict.combo.df, aes(x=Gall_Height_mm, y=Density_per_100_shoots, z=standard.fit, color=standard.fit)) + 
  geom_point(size=10) +
  scale_color_gradientn(colors = viridis6()) +
  facet_wrap(~treatment)


## ECTOPARASITOIDS ON CONTROL TREES ----
#control.ecto.gppr <- gppr(y = "ectos", xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), data = control.df, nterms = 1)
#control.ecto.gppr$ppr$alpha
#control.ecto.gppr$ppr$beta
#control.df$term1_ectos <- control.ecto.gppr$ppr$alpha[1]*control.df$sc.Gall_Height_mm + control.ecto.gppr$ppr$alpha[2]*control.df$sc.gall_individuals + control.ecto.gppr$ppr$alpha[3]*control.df$sc.Density_per_100_shoots

#control.ecto.major.gam <- gamm.model(ectos ~ s(term1_ectos), data = control.df)
#summary(control.ecto.major.gam$gam)
#summary(control.ecto.major.gam$mer)
#gamm.plot(control.ecto.major.gam)

control.ecto.component.gam <- gamm.model(ectos ~ s(Gall_Height_mm) + s(gall_individuals, k=9) + s(Density_per_100_shoots), data = control.df)
summary(control.ecto.component.gam$gam)
summary(control.ecto.component.gam$mer)
concurvity(control.ecto.component.gam$gam)
gamm.plot(control.ecto.component.gam, pages=1)

control.ecto.noplaty.component.gam <- gamm.model(ectos ~ s(Gall_Height_mm) + s(gall_individuals, k=9) + s(Density_per_100_shoots), data = filter(control.df, ectos > 0 | pupa > 0))
summary(control.ecto.noplaty.component.gam$gam)
summary(control.ecto.noplaty.component.gam$mer)
concurvity(control.ecto.noplaty.component.gam$gam)
gamm.plot(control.ecto.noplaty.component.gam, pages=1)

## PLATYGASTER ON CONTROL TREES ----
#control.platy.gppr <- gppr(y = "egg_parasitoid", xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), data = control.df, nterms = 1) 
#control.platy.gppr$ppr$alpha
#control.platy.gppr$ppr$beta
#control.df$term1_platys <- control.platy.gppr$ppr$alpha[1]*control.df$sc.Gall_Height_mm + control.platy.gppr$ppr$alpha[2]*control.df$sc.gall_individuals + control.platy.gppr$ppr$alpha[3]*control.df$sc.Density_per_100_shoots
# control.df$term2_platys <- control.platy.gppr$ppr$alpha[4]*control.df$sc.Gall_Height_mm + control.platy.gppr$ppr$alpha[5]*control.df$sc.gall_individuals + control.platy.gppr$ppr$alpha[6]*control.df$sc.Density_per_100_shoots
# plot(term2_platys ~ term1_platys, control.df)
# cor.test(control.df$term1_platys, control.df$term2_platys) # cor = -0.69

#control.platy.major.gam <- gamm.model(egg_parasitoid ~ s(term1_platys), data = control.df)
#summary(control.platy.major.gam$gam)
#summary(control.platy.major.gam$mer)
#gamm.plot(control.platy.major.gam)

control.platy.component.gam <- gamm.model(egg_parasitoid ~ s(Gall_Height_mm) + s(gall_individuals, k=9) + s(Density_per_100_shoots), data = control.df)
summary(control.platy.component.gam$gam)
summary(control.platy.component.gam$mer)
gamm.plot(control.platy.component.gam, pages=1)

control.platy.noectos.component.gam <- gamm.model(egg_parasitoid ~ s(Gall_Height_mm) + s(gall_individuals, k=9) + s(Density_per_100_shoots), data = filter(control.df, platy > 0 | pupa > 0))
summary(control.platy.noectos.component.gam$gam)
summary(control.platy.noectos.component.gam$mer)
gamm.plot(control.platy.noectos.component.gam, pages=1)


## PLATYGASTER ON ECTOPARASITOID EXCLUSION TREES ----
#treatment.platy.gppr <- gppr(y = "egg_parasitoid", xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), data = treatment.df, nterms = 1)
#treatment.platy.gppr$ppr$alpha
#treatment.platy.gppr$ppr$beta
#treatment.df$term1_platys <- treatment.platy.gppr$ppr$alpha[1]*treatment.df$sc.Gall_Height_mm + treatment.platy.gppr$ppr$alpha[2]*treatment.df$sc.gall_individuals + treatment.platy.gppr$ppr$alpha[3]*treatment.df$sc.Density_per_100_shoots

#treatment.platy.major.gam <- gamm.model(egg_parasitoid ~ s(term1_platys), data = treatment.df)
#summary(treatment.platy.major.gam$gam)
#summary(treatment.platy.major.gam$mer)
#gamm.plot(treatment.platy.major.gam)

treatment.platy.component.gam <- gamm.model(egg_parasitoid ~ s(Gall_Height_mm) + s(gall_individuals) + s(Density_per_100_shoots), data = treatment.df)
summary(treatment.platy.component.gam$gam)
summary(treatment.platy.component.gam$mer)
gamm.plot(treatment.platy.component.gam, pages=1)


test.platy.component.gam <- gamm.model(egg_parasitoid ~ s(Gall_Height_mm, by=Treatment.focus) + s(gall_individuals, by=Treatment.focus) + s(Density_per_100_shoots, by=Treatment.focus), data = gall_selection.df)
summary(test.platy.component.gam$gam)
summary(test.platy.component.gam$mer)
gamm.plot(test.platy.component.gam, pages=1)

## PLOTS COMPARING EGG AND LARVAL PARASITOIDS ON CONTROL TREES 

# egg parasitoid predictions CONTROL
predict.EggPtoid.height.control <- predict(control.platy.component.gam$gam, newdata = newdata.height.control, type = "response", se.fit = TRUE)
predict.EggPtoid.height.control.df <- cbind.data.frame(predict.EggPtoid.height.control, newdata.height.control)

predict.EggPtoid.indiv.control <- predict(control.platy.component.gam$gam, newdata = newdata.indiv.control, type = "response", se.fit = TRUE)
predict.EggPtoid.indiv.control.df <- cbind.data.frame(predict.EggPtoid.indiv.control, newdata.indiv.control)

predict.EggPtoid.density.control <- predict(control.platy.component.gam$gam, newdata = newdata.density.control, type = "response", se.fit = TRUE)
predict.EggPtoid.density.control.df <- cbind.data.frame(predict.EggPtoid.density.control, newdata.density.control)

# larva parasitoid predictions CONTROL
predict.LarvaPtoid.height.control <- predict(control.ecto.component.gam$gam, newdata = newdata.height.control, type = "response", se.fit = TRUE)
predict.LarvaPtoid.height.control.df <- cbind.data.frame(predict.LarvaPtoid.height.control, newdata.height.control)

predict.LarvaPtoid.indiv.control <- predict(control.ecto.component.gam$gam, newdata = newdata.indiv.control, type = "response", se.fit = TRUE)
predict.LarvaPtoid.indiv.control.df <- cbind.data.frame(predict.LarvaPtoid.indiv.control, newdata.indiv.control)

predict.LarvaPtoid.density.control <- predict(control.ecto.component.gam$gam, newdata = newdata.density.control, type = "response", se.fit = TRUE)
predict.LarvaPtoid.density.control.df <- cbind.data.frame(predict.LarvaPtoid.density.control, newdata.density.control)

# egg parasitoid predictions TREATMENT
predict.EggPtoid.height.treatment <- predict(treatment.platy.component.gam$gam, newdata = newdata.height.treatment, type = "response", se.fit = TRUE)
predict.EggPtoid.height.treatment.df <- cbind.data.frame(predict.EggPtoid.height.treatment, newdata.height.treatment)

predict.EggPtoid.indiv.treatment <- predict(treatment.platy.component.gam$gam, newdata = newdata.indiv.treatment, type = "response", se.fit = TRUE)
predict.EggPtoid.indiv.treatment.df <- cbind.data.frame(predict.EggPtoid.indiv.treatment, newdata.indiv.treatment)

predict.EggPtoid.density.treatment <- predict(treatment.platy.component.gam$gam, newdata = newdata.density.treatment, type = "response", se.fit = TRUE)
predict.EggPtoid.density.treatment.df <- cbind.data.frame(predict.EggPtoid.density.treatment, newdata.density.treatment)

# combine data
predict.ptoid.height.df <-  bind_rows(mutate(predict.EggPtoid.height.control.df, treatment = "Complex Platy"),
                                      mutate(predict.LarvaPtoid.height.control.df, treatment = "Complex Larva"),
                                      mutate(predict.EggPtoid.height.treatment.df, treatment = "Simple Platy")) 
predict.ptoid.indiv.df <-  bind_rows(mutate(predict.EggPtoid.indiv.control.df, treatment = "Complex Platy"),
                                     mutate(predict.LarvaPtoid.indiv.control.df, treatment = "Complex Larva"),
                                     mutate(predict.EggPtoid.indiv.treatment.df, treatment = "Simple Platy")) 
predict.ptoid.density.df <-  bind_rows(mutate(predict.EggPtoid.density.control.df, treatment = "Complex Platy"),
                                       mutate(predict.LarvaPtoid.density.control.df, treatment = "Complex Larva"),
                                       mutate(predict.EggPtoid.density.treatment.df, treatment = "Simple Platy")) 
predict.ptoid.all.df <- bind_rows(predict.ptoid.height.df, predict.ptoid.indiv.df, predict.ptoid.density.df)

# plots
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

fitness.ptoid.height <- ggplot(predict.ptoid.height.df, aes(x = Gall_Height_mm, y = fit, color = treatment)) + 
  geom_line(size = 1.5) + 
  geom_line(aes(y = fit + se.fit), linetype = "dotted") +
  geom_line(aes(y = fit - se.fit), linetype = "dotted") +
  ylab("") + xlab("Gall diameter (mm)") + scale_color_manual(values = cbPalette[c(2,7,3)], guide = 'none') +
  scale_y_continuous(limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) + 
  scale_x_continuous(limits=c(mean_Gall_Height_mm-sd_Gall_Height_mm, mean_Gall_Height_mm+sd_Gall_Height_mm))

fitness.ptoid.indiv <- ggplot(predict.ptoid.indiv.df, aes(x = gall_individuals, y = fit, color = treatment)) + 
  geom_line(size = 1.5) + 
  geom_line(aes(y = fit + se.fit), linetype = "dotted") +
  geom_line(aes(y = fit - se.fit), linetype = "dotted") +
  ylab("") + xlab("Larva per gall") + scale_color_manual(values = cbPalette[c(2,7,3)], guide = 'none') +
  scale_y_continuous(limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) + 
  scale_x_continuous(limits=c(mean_gall_individuals-sd_gall_individuals, mean_gall_individuals+sd_gall_individuals))

fitness.ptoid.density <- ggplot(predict.ptoid.density.df, aes(x = Density_per_100_shoots, y = fit, color = treatment)) + 
  geom_line(size = 1.5) + 
  geom_line(aes(y = fit + se.fit), linetype = "dotted") +
  geom_line(aes(y = fit - se.fit), linetype = "dotted") +
  ylab("") + xlab("Larva per 100 shoots") + scale_color_manual(values = cbPalette[c(2,7,3)], guide = 'none') +
  scale_y_continuous(limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) + 
  scale_x_continuous(limits=c(mean_Density_per_100_shoots-sd_Density_per_100_shoots, mean_Density_per_100_shoots+sd_Density_per_100_shoots))

fitness.ptoid.plot_allGalls <- plot_grid(fitness.ptoid.height, fitness.ptoid.indiv, fitness.ptoid.density, nrow=1)
save_plot("ptoid_selection_gradients_allGalls.pdf", fitness.ptoid.plot_allGalls, base_width = 11, base_height = 8.5)
