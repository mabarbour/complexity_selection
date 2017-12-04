
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

## GENERALIZED PROJECTION PURSUIT REGRESSION ANALYSIS ----

## GALLS ON CONTROL TREES
control.df <- as.data.frame(filter(gall_selection.df, Treatment.focus == "Control"))
control.gppr <- gppr(y = "gall_survival", 
                     xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), 
                     data = control.df, nterms = 1)
control.gppr$ppr$alpha  
control.gppr$ppr$beta

# se.method = 'boot.para' takes a while so be prepared. Note that it detected non-linear selection on gall height, which is worth following up on...
# strong evidence for directional selection for larger galls, and some evidence for non-linear selection on gall size (diversifying???)
gppr.gradients(mod=control.gppr,phenotype=c("sc.Gall_Height_mm","sc.gall_individuals"), covariates = "sc.Density_per_100_shoots", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests
gppr.gradients(mod=control.gppr,phenotype=c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests
gppr.gradients(mod=control.gppr,phenotype=c("sc.gall_individuals","sc.Density_per_100_shoots"), covariates = "sc.Gall_Height_mm", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests

control.df$term1 <- control.gppr$ppr$alpha[1]*control.df$sc.Gall_Height_mm + control.gppr$ppr$alpha[2]*control.df$sc.gall_individuals + control.gppr$ppr$alpha[3]*control.df$sc.Density_per_100_shoots
# control.df$term2 <- control.gppr$ppr$alpha[4]*control.df$sc.Gall_Height_mm + control.gppr$ppr$alpha[5]*control.df$sc.gall_individuals + control.gppr$ppr$alpha[6]*control.df$sc.Density_per_100_shoots
# plot(term2 ~ term1, control.df)
# cor.test(control.df$term1, control.df$term2) 

control.major.gam <- gamm4(gall_survival ~ s(term1), random = ~(1|Genotype/Plant_Position/Gall_Number),
                           data = control.df, family = "binomial", method = "GCV.Cp")
summary(control.major.gam)#$gam)
# summary(control.major.gam$mer)
# plot(control.major.gam$gam, seWithMean = T, shift = mean(predict(control.major.gam$gam)), trans = function(x) {exp(x)/(1+exp(x))})
plot(control.major.gam$gam, residuals = T, pch = 19, cex = 0.3, seWithMean = T, shift = mean(predict(control.major.gam$gam)), trans = function(x) {exp(x)/(1+exp(x))})

rsd <- residuals(control.major.gam$gam)
plot(rsd ~ cut(control.df$sc.Gall_Height_mm, breaks = c(-5,-2,-1,1,2,5)))
plot(rsd ~ cut(control.df$sc.gall_individuals, breaks = c(-5,-2,-1,1,2,5)))
plot(rsd ~ cut(control.df$sc.Density_per_100_shoots, breaks = c(-5,-2,-1,1,2,5)))

control.major.gradients <- gam.gradients(mod = control.major.gam$gam, phenotype = c("term1"), refit.smooth = T, se.method = 'boot.para', standardized = T, parallel = 'multicore', ncpus = 32)
control.major.gradients$ests
hist(control.major.gradients$boot[,1])
abline(v = control.major.gradients$ests[1,1])

hist(control.major.gradients$boot[,2])
abline(v = control.major.gradients$ests[2,1])

library(brms)
control.major.brm <- brm(gall_survival ~ s(term1), #random = ~(1|Genotype/Plant_Position/Gall_Number),
                         data = control.df, 
                         family = "bernoulli",
                         control = list(adapt_delta = 0.999))
summary(control.major.brm)
plot(marginal_effects(control.major.brm, spaghetti = T, nsamples = 100), line_args = list(color = "black"))

obs <- control.df$gall_survival
pp_check(control.major.brm, nsamples = 100)

pp_control.major <- posterior_predict(control.major.brm, nsamples = 100)

prob_attack <- function(y) sum(y)/length(y)
ppc_stat_grouped(y = control.df$gall_survival, 
                 yrep = pp_control.major, 
                 group = control.df$Genotype,
                 stat = "prob_attack")
ppc_stat_grouped(y = control.df$gall_survival, 
                 yrep = pp_control.major, 
                 group = cut(control.df$sc.Gall_Height_mm, breaks = c(-5,-1,1,5)),
                 stat = "prob_attack")
ppc_stat_grouped(y = control.df$gall_survival, 
                 yrep = pp_control.major, 
                 group = cut(control.df$sc.gall_individuals, breaks = c(-5,-1,1,5)),
                 stat = "prob_attack")
ppc_stat_grouped(y = control.df$gall_survival, 
                 yrep = pp_control.major, 
                 group = cut(control.df$sc.Density_per_100_shoots, breaks = c(-5,-1,1,5)),
                 stat = "prob_attack")


obs_mean <- mean(obs)
pp_means <- apply(pp_control.major, MARGIN = 1, FUN = mean)

pp_data <- data.frame(t(pp_control.major)) %>% 
  mutate(sc.Gall_Height_mm = control.df$sc.Gall_Height_mm,
         sc.gall_individuals = control.df$gall_individuals,
         sc.Density_per_100_shoots = control.df$sc.Density_per_100_shoots,
         Genotype = control.df$Genotype) %>%
  gather(key = rep, value = gall_survival, X1:X100)# %>%
  #mutate(color = c(rep("obs",703), rep("sims",70300)))

binomial_smooth.PP <- function(...) {
  geom_line(stat = "smooth", method = "glm", method.args = list(family = "binomial"), ...)
}
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}
ggplot(pp_data, 
       aes(x = sc.Gall_Height_mm, y = gall_survival)) +
  binomial_smooth.PP(aes(group = rep), se = F, alpha = 0.1, color = "blue") +
  binomial_smooth(data = control.df, color = "black", se = T) + 
  facet_wrap(~Genotype)

# the model can reproduce the overall average
ggplot(data.frame(pp_means = pp_means), aes(y = pp_means, x = 1)) + geom_boxplot() + geom_hline(yintercept = obs_mean, color = "blue")
obs_mean
boxplot(pp_means)

library(bayesplot)
ppc_rootogram(y = control.df$gall_survival, yrep = pp_control.major)
ppc_violin_grouped(y = control.df$gall_survival, yrep = pp_control.major, group = control.df$Genotype)
ppc_scatter_avg_grouped(y = control.df$gall_survival, 
                        yrep = pp_control.major,
                        group = cut(control.df$sc.Gall_Height_mm, breaks = c(-5,-1,1,5)))
# in combination with the strong loading of Gall_Height_mm, there appears to be a strong directional selection gradient acting primarily on gall size in the complex food web.
gam.gradients(mod = control.major.gam, phenotype = c("term1"), se.method = 'boot.para', standardized = T, parallel = 'multicore', ncpus = 32)$ests

control.Gall_Height.gam <- gam(gall_survival ~ s(Gall_Height_mm), #random = ~(1|Genotype/Plant_Position/Gall_Number),
                         data = control.df, family = "binomial")
plot(control.Gall_Height.gam, seWithMean = T, shift = mean(predict(control.Gall_Height.gam)), trans = function(x) {exp(x)/(1+exp(x))})

# still don't understand the non-linear form of selection going on here...
gam.gradients(mod = control.Gall_Height.gam, phenotype = c("Gall_Height_mm"), se.method = 'boot.para', standardized = T, parallel = 'multicore', ncpus = 32)$ests

# plots to examine phenotypic distributions. Important for fitness landscapes to only plot where I have most of my data.
ggplot(control.df, aes(x = sc.Gall_Height_mm, y = sc.Density_per_100_shoots)) +
  geom_point(shape = 1)
ggplot(control.df, aes(x = sc.Gall_Height_mm, y = sc.gall_individuals)) +
  geom_jitter(shape = 1, height = 0.05)
ggplot(control.df, aes(x = sc.Density_per_100_shoots, y = sc.gall_individuals)) +
  geom_jitter(shape = 1, height = 0.05)


## GALLS ON ECTOPARSITOID EXCLUSION TREES
treatment.df <- as.data.frame(filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion"))
treatment.gppr <- gppr(y = "gall_survival", xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), 
                     data = treatment.df, nterms = 1)
treatment.gppr$ppr$alpha # strong loadings of all three traits, with gall height acting in opposing directions to gall individuals and gall density
treatment.gppr$ppr$beta

# INTERESTING, CALCULATING GRADIENTS WITH 1 OR 2 TERMS MAKES A BIG DIFFERENCE IN TERMS OF SIGNIFICANCE. BUT WITH 2, THE SE SEEM UNECESSARILY LARGE, POSSIBLE DU TO CORRELATIONS???
gppr.size_indiv <- gppr.gradients(mod=treatment.gppr,phenotype=c("sc.Gall_Height_mm","sc.gall_individuals"), covariates = "sc.Density_per_100_shoots", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)#$ests

gppr.size_indiv$ests$estimates[4]
median(gppr.size_indiv$boot[ ,4])
hist(gppr.size_indiv$boot[ ,4])
abline(v=gppr.size_indiv$ests$estimates[4])
sd(gppr.size_indiv$boot[ ,5])
length(which(gppr.size_indiv$boot[ ,5] > gppr.size_indiv$ests$estimates[5]))/1000 # 97% of bootsrapped values are greater than the estimate
quantile(gppr.size_indiv$boot[ ,5], probs = c(0.975,0.025)) # hmmm, estimates falls within 95% interval...
mean(gppr.size_indiv$boot[ ,5])

# current method for calculating p-values. This is very inaccurate if the bootstraps represent actual gradient estimates...
2*min(sum((gppr.size_indiv$boot[ ,5]>0)+0), sum((gppr.size_indiv$boot[ ,5]<0)+0))/length(gppr.size_indiv$boot[ ,5])

2*min(sum(gppr.size_indiv$boot[ ,5]>gppr.size_indiv$ests$estimates[5]), sum(gppr.size_indiv$boot[ ,5]<gppr.size_indiv$ests$estimates[5]))/length(gppr.size_indiv$boot[ ,5])


gppr.gradients(mod=treatment.gppr,phenotype=c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests
gppr.gradients(mod=treatment.gppr,phenotype=c("sc.gall_individuals","sc.Density_per_100_shoots"), covariates = "sc.Gall_Height_mm", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests

treatment.df$term1 <- treatment.gppr$ppr$alpha[1]*treatment.df$sc.Gall_Height_mm + treatment.gppr$ppr$alpha[2]*treatment.df$sc.gall_individuals + treatment.gppr$ppr$alpha[3]*treatment.df$sc.Density_per_100_shoots
# treatment.df$term2 <- treatment.gppr$ppr$alpha[4]*treatment.df$sc.Gall_Height_mm + treatment.gppr$ppr$alpha[5]*treatment.df$sc.gall_individuals + treatment.gppr$ppr$alpha[6]*treatment.df$sc.Density_per_100_shoots
# plot(term2 ~ term1, treatment.df)
# cor.test(treatment.df$term1, treatment.df$term2) 


treatment.major.gam <- gamm4(gall_survival ~ s(term1), random = ~(1|Genotype/Plant_Position/Gall_Number),
                           data = treatment.df, family = "binomial")
summary(treatment.major.gam$gam)
summary(treatment.major.gam$mer)
plot(treatment.major.gam$gam, seWithMean = T, shift = mean(predict(treatment.major.gam$gam)), trans = function(x) {exp(x)/(1+exp(x))})
#plot(treatment.major.gam, seWithMean = T, shift = mean(predict(treatment.major.gam)), trans = function(x) {exp(x)/(1+exp(x))})

#gam.gradients(mod = treatment.major.gam$gam, phenotype = c("term1"), se.method = 'n', standardized = T)
treatment.major.gradients <- gam.gradients(mod = treatment.major.gam$gam, phenotype = c("term1"), refit.smooth = T, se.method = 'boot.para', standardized = T, parallel = 'multicore', ncpus = 32)
treatment.major.gradients$ests
mean(treatment.major.gradients$boot[,1])
mean(treatment.major.gradients$boot[,2])
hist(treatment.major.gradients$boot[,2])

treatment.all.gam <- gamm4(gall_survival ~ s(sc.Gall_Height_mm) + s(sc.gall_individuals) + s(sc.Density_per_100_shoots), 
                           random = ~(1|Genotype/Plant_Position/Gall_Number),
                           data = treatment.df, family = "binomial")
plot(treatment.all.gam$gam, seWithMean = T, shift = mean(predict(treatment.all.gam$gam)), trans = function(x) {exp(x)/(1+exp(x))})

treatment.size_indiv.gradients <- gam.gradients(mod=treatment.all.gam$gam,phenotype=c("sc.Gall_Height_mm","sc.gall_individuals"), refit.smooth = T, covariates = "sc.Density_per_100_shoots", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)
treatment.size_indiv.gradients$ests
hist(treatment.size_indiv.gradients$boot[ ,4])

gam.gradients(mod=treatment.all.gam$gam,phenotype=c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests
gam.gradients(mod=treatment.all.gam,phenotype=c("sc.gall_individuals","sc.Density_per_100_shoots"), covariates = "sc.Gall_Height_mm", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests

vis.gam(treatment.all.gam, view = c("sc.Gall_Height_mm","sc.gall_individuals"), type = "response", plot.type = "contour")
vis.gam(treatment.all.gam, view = c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), type = "response", plot.type = "contour")
vis.gam(treatment.all.gam, view = c("sc.gall_individuals","sc.Density_per_100_shoots"), type = "response", plot.type = "contour")

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
control.ecto.gppr <- gppr(y = "ectos", 
                          xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), 
                          data = control.df, nterms = 1)
control.ecto.gppr$ppr$alpha
control.ecto.gppr$ppr$beta

gppr.gradients(mod=control.ecto.gppr,phenotype=c("sc.Gall_Height_mm","sc.gall_individuals"), covariates = "sc.Density_per_100_shoots", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests
gppr.gradients(mod=control.ecto.gppr,phenotype=c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests
gppr.gradients(mod=control.ecto.gppr,phenotype=c("sc.gall_individuals","sc.Density_per_100_shoots"), covariates = "sc.Gall_Height_mm", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests

control.df$term1_ectos <- control.ecto.gppr$ppr$alpha[1]*control.df$sc.Gall_Height_mm + control.ecto.gppr$ppr$alpha[2]*control.df$sc.gall_individuals + control.ecto.gppr$ppr$alpha[3]*control.df$sc.Density_per_100_shoots

control.ecto.major.gam <- gam(ectos ~ s(term1_ectos), random = ~(1|Genotype/Plant_Position/Gall_Number),
                         data = control.df, family = "binomial")
summary(control.ecto.major.gam$gam)
gam.check(control.ecto.major.gam$gam)

rsd <- residuals(control.ecto.major.gam$gam)
plot(rsd ~ cut(control.df$sc.Gall_Height_mm, breaks = c(-5,-1,1,5)))
plot(rsd ~ cut(control.df$sc.gall_individuals, breaks = c(-5,-1,1,5)))
plot(rsd ~ cut(control.df$sc.Density_per_100_shoots, breaks = c(-5,-1,1,5)))
plot(rsd ~ control.df$sc.gall_individuals)
plot(rsd ~ control.df$sc.Density_per_100_shoots)
plot(rsd ~ control.df$Genotype)
plot(rsd ~ control.df$Plant_Position)
plot(rsd ~ control.df$Gall_Number)
gam(rsd~s(term1_ectos,k=20),data = control.df, gamma = 1.4)

plot(control.ecto.major.gam$gam, residuals = T, pch = 19, cex = 0.3, seWithMean = T, shift = mean(predict(control.ecto.major.gam$gam)), trans = function(x) {exp(x)/(1+exp(x))})
#gam.gradients(mod = control.ecto.major.gam$gam, phenotype = c("term1_ectos"), se.method = 'boot.para', standardized = T, parallel = 'multicore', ncpus = 32)$ests

control.ecto.major.gam$sp
lambda <- exp(seq(-4,0, by=.05))        # fit a range of lambdas >0
gcvscore <- sapply(lambda, function(lambda, control.df){
  gam(ectos ~ s(term1_ectos), data = control.df, family = binomial, 
      sp = lambda, method="GCV.Cp")$gcv.ubre},control.df)

plot(lambda, gcvscore, type = "l") # or
plot(log(lambda), gcvscore, type = "l")

## PLATYGASTER ON CONTROL TREES
control.platy.gppr <- gppr(y = "egg_parasitoid", 
                           xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), 
                           data = control.df, nterms = 1) 
control.platy.gppr$ppr$alpha
control.platy.gppr$ppr$beta

gppr.gradients(mod=control.platy.gppr,phenotype=c("sc.Gall_Height_mm","sc.gall_individuals"), covariates = "sc.Density_per_100_shoots", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests
gppr.gradients(mod=control.platy.gppr,phenotype=c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests
gppr.gradients(mod=control.platy.gppr,phenotype=c("sc.gall_individuals","sc.Density_per_100_shoots"), covariates = "sc.Gall_Height_mm", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests

control.df$term1_platys <- control.platy.gppr$ppr$alpha[1]*control.df$sc.Gall_Height_mm + control.platy.gppr$ppr$alpha[2]*control.df$sc.gall_individuals + control.platy.gppr$ppr$alpha[3]*control.df$sc.Density_per_100_shoots
#control.df$term2_platys <- control.platy.gppr$ppr$alpha[4]*control.df$sc.Gall_Height_mm + control.platy.gppr$ppr$alpha[5]*control.df$sc.gall_individuals + control.platy.gppr$ppr$alpha[6]*control.df$sc.Density_per_100_shoots
#plot(term2_platys ~ term1_platys, control.df)
#cor.test(control.df$term1_platys, control.df$term2_platys) # cor = -0.69

control.platy.major.gam <- gamm4(egg_parasitoid ~ s(term1_platys), random = ~(1|Genotype/Plant_Position/Gall_Number),
                                data = control.df, family = "binomial")
summary(control.platy.major.gam)
gam.check(control.platy.major.gam$gam)
qq.gam(control.platy.major.gam$gam, pch = 19, cex = 0.3)
plot(control.platy.major.gam$gam, residuals = T, pch = 19, cex = 0.3, seWithMean = T, shift = mean(predict(control.platy.major.gam$gam)), trans = function(x) {exp(x)/(1+exp(x))})
gam.gradients(mod = control.platy.major.gam, phenotype = c("term1_platys"), se.method = 'boot.para', standardized = T, parallel = 'multicore', ncpus = 32)$ests


## PLATYGASTER ON ECTOPARASITOID EXCLUSION TREES
treatment.platy.gppr <- gppr(y = "egg_parasitoid", 
                             xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), 
                             data = treatment.df, nterms = 1)
treatment.platy.gppr$ppr$alpha
treatment.platy.gppr$ppr$beta

gppr.gradients(mod=treatment.platy.gppr,phenotype=c("sc.Gall_Height_mm","sc.gall_individuals"), covariates = "sc.Density_per_100_shoots", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests
gppr.gradients(mod=treatment.platy.gppr,phenotype=c("sc.Gall_Height_mm","sc.Density_per_100_shoots"), covariates = "sc.gall_individuals", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests
gppr.gradients(mod=treatment.platy.gppr,phenotype=c("sc.gall_individuals","sc.Density_per_100_shoots"), covariates = "sc.Gall_Height_mm", se.method='boot.para',standardize=FALSE, parallel = 'multicore', ncpus = 32)$ests

treatment.df$term1_platys <- treatment.platy.gppr$ppr$alpha[1]*treatment.df$sc.Gall_Height_mm + treatment.platy.gppr$ppr$alpha[2]*treatment.df$sc.gall_individuals + treatment.platy.gppr$ppr$alpha[3]*treatment.df$sc.Density_per_100_shoots

treatment.platy.major.gam <- gam(egg_parasitoid ~ s(term1_platys), #random = ~(1|Genotype/Plant_Position/Gall_Number),
                                 data = treatment.df, family = "binomial")
summary(treatment.platy.major.gam)
plot(treatment.platy.major.gam, seWithMean = T, shift = mean(predict(treatment.platy.major.gam)), trans = function(x) {exp(x)/(1+exp(x))})
gam.gradients(mod = treatment.platy.major.gam, phenotype = c("term1_platys"), se.method = 'boot.para', standardized = T, parallel = 'multicore', ncpus = 32)$ests

treatment.platy.all.gam <- gam(egg_parasitoid ~ s(sc.Gall_Height_mm) + s(sc.gall_individuals) + s(sc.Density_per_100_shoots), 
                                 #random = ~(1|Genotype/Plant_Position/Gall_Number),
                                 data = treatment.df, family = "binomial")
summary(treatment.platy.all.gam$gam)
plot(treatment.platy.all.gam, seWithMean = T, shift = mean(predict(treatment.platy.all.gam)), trans = function(x) {exp(x)/(1+exp(x))})
#gam.gradients(mod = treatment.platy.all.gam$gam, phenotype = c("term1_platys"), se.method = 'boot.para', standardized = T, parallel = 'multicore', ncpus = 32)$ests


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
