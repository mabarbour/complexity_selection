
## LOAD REQUIRED LIBRARIES & SET OPTIONS ----
library(tidyverse)
library(cowplot) # pretty default ggplots
#library(broom) # for tidying multiple linear models
#library(mgcv) # generalized additive models
#library(gsg) # calculating selection gradients
library(lme4)
library(visreg)
library(viridis)

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

## LOAD & MANAGE DATASET ----

gall_selection.df <- read_csv("gall_selection_data.csv") %>%
  # convert appropriate variables to characters instead of integers
  mutate(Plant_Position = as.character(Plant_Position),
         Gall_Number = as.character(Gall_Number)) %>%
  unite(Gall_ID, Gall_Number, Gall_Letter, remove = FALSE) %>%
  
  # subset data for analysis
  filter(phenology == "early", Location == "tree",
         platy > 0 | ectos > 0 | pupa > 0) %>%                  # eliminate unknown sources of mortality
  mutate(gall_survival = as.numeric(ifelse(pupa > 0, 1, 0)),
         egg_parasitoid = as.numeric(ifelse(platy > 0, 1, 0))) 

# summarize at gall-level for gall individuals
all_gall <- gall_selection.df %>%
  group_by(Genotype, Treatment.focus, Plant_Position, Gall_Number) %>%
  summarise(gall_individuals = mean(gall_individuals)) %>%
  group_by(Genotype, Treatment.focus, Plant_Position) %>%
  summarise(gall_sample_size = n(),
            gall_survival_mean.gall.level = mean(gall_survival), 
            gall_individuals_mean = mean(gall_individuals),
            gall_individuals_var = var(gall_individuals),
            gall_individuals_mean.log = mean(log(gall_individuals)),
            gall_individuals_var.log = var(log(gall_individuals))) %>%
  ungroup() 


# summarize each data set at tree-level
all_tree <- gall_selection.df %>%
  group_by(Genotype, Treatment.focus, Plant_Position) %>%
  summarise(gall_survival_mean = mean(gall_survival), 
            individual_sample_size = n(),
            Density_per_100_shoots = mean(Density_per_100_shoots), # note there is no variance at plant-level for Density_per_100_shoots
            Gall_Height_mm_mean = mean(Gall_Height_mm),
            Gall_Height_mm_var = var(Gall_Height_mm)) %>%               # for testing assumption of no mean-variance relationship
  ungroup() %>%
  left_join(., all_gall)

## TESTING ASSUMPTIONS ----

## 1. Traits are approximately multivariate normal. 
## A necessary, but not sufficient, condition for this is that the phenotypes are normally distributed
## within a population. Here, the population represents the galling insects found on a tree.

hist(all_tree$individual_sample_size)
summary(all_tree$individual_sample_size)

hist(all_tree$gall_sample_size)
summary(all_tree$gall_sample_size)

#### TEST THE ASSUMPTION THAT MEAN IS NOT RELATED TO VARIANCE
all_tree %>%
  ggplot(., aes(x=Gall_Height_mm_mean, y=Gall_Height_mm_var, size=individual_sample_size)) +
  geom_point() + geom_smooth(method="lm")

summary(lm(Gall_Height_mm_var ~ Gall_Height_mm_mean, data=all_tree, weights=individual_sample_size))

# very strong positive correlation...
all_tree %>%
  ggplot(., aes(x=gall_individuals_mean, y=gall_individuals_var)) +
  geom_point() + geom_smooth(method="lm")

all_tree %>%
  ggplot(., aes(x=gall_individuals_mean.log, y=gall_individuals_var.log, size=gall_sample_size, weight=gall_sample_size)) +
  geom_point() + geom_smooth(method="lm")

summary(lm(gall_individuals_var.log ~ gall_individuals_mean.log, data=all_tree))
summary(lm(gall_individuals_var.log ~ gall_individuals_mean.log, data=all_tree, weights = gall_sample_size))

# note that there is no variance in Density_per_100_shoots, since it was measured at the tree-level

## ANALYSES ----
table(all_tree$gall_survival_mean)

all_tree <- all_tree %>%
  mutate(c.Gall_Height_mm_mean = Gall_Height_mm_mean - mean(Gall_Height_mm_mean),
         c.gall_individuals_mean.log = gall_individuals_mean.log - mean(gall_individuals_mean.log),
         c.Density_per_100_shoots = Density_per_100_shoots - mean(Density_per_100_shoots))
all_tree %>%
  group_by(Treatment.focus) %>%
  summarise(mean.fitness = mean(gall_survival_mean))
all_tree$relative.fit <- with(all_tree, ifelse(Treatment.focus=="Control", gall_survival_mean/0.415, gall_survival_mean/0.546))

landscape.lm <- glm(relative.fit ~ 
                      c.Gall_Height_mm_mean + c.gall_individuals_mean.log + c.Density_per_100_shoots +
                      c.Gall_Height_mm_mean:Treatment.focus + c.gall_individuals_mean.log:Treatment.focus + c.Density_per_100_shoots:Treatment.focus,
                   data = all_tree)
summary(landscape.lm)
car::vif(landscape.lm)
plot(landscape.lm)
summary(landscape.lm)
anova(landscape.lm)
car::Anova(landscape.lm, type=2)
car::Anova(landscape.lm, type=3)


## CONTROL TREES 
control_tree <- all_tree %>%
  filter(Treatment.focus=="Control") %>%
  mutate(c.Gall_Height_mm_mean = Gall_Height_mm_mean - mean(Gall_Height_mm_mean),
         c.gall_individuals_mean.log = gall_individuals_mean.log - mean(gall_individuals_mean.log),
         c.Density_per_100_shoots = Density_per_100_shoots - mean(Density_per_100_shoots),
         c.Density_per_100_shoots.sqrt = sqrt(Density_per_100_shoots) - mean(sqrt(Density_per_100_shoots)))

landscape.control.lm <- lm(log(gall_survival_mean+1)~(c.Gall_Height_mm_mean+c.gall_individuals_mean.log+c.Density_per_100_shoots.sqrt)^2,
                   data = control_tree, weights = individual_sample_size)
car::vif(landscape.control.lm)
plot(landscape.control.lm)
summary(landscape.control.lm)
anova(landscape.control.lm)
car::Anova(landscape.control.lm, type=2)
car::Anova(landscape.control.lm, type=3)
visreg(landscape.control.lm)

## TREATMENT TREES
treatment_tree <- all_tree %>%
  filter(Treatment.focus=="Ectoparasitoid exclusion") %>%
  mutate(c.Gall_Height_mm_mean = Gall_Height_mm_mean - mean(Gall_Height_mm_mean),
         c.gall_individuals_mean.log = gall_individuals_mean.log - mean(gall_individuals_mean.log),
         c.Density_per_100_shoots = Density_per_100_shoots - mean(Density_per_100_shoots),
         c.Density_per_100_shoots.sqrt = sqrt(Density_per_100_shoots) - mean(sqrt(Density_per_100_shoots)))

landscape.treatment.lm <- lm(gall_survival_mean~(c.Gall_Height_mm_mean+c.gall_individuals_mean.log+c.Density_per_100_shoots.sqrt),
                           data = treatment_tree)
coef(landscape.treatment.lm)
landscape.treatment.lmer <- lmer(gall_survival_mean~(c.Gall_Height_mm_mean+c.gall_individuals_mean.log+c.Density_per_100_shoots.sqrt)+(1|Genotype),
                                data = treatment_tree)
fixef(landscape.treatment.lmer)
car::vif(landscape.treatment.lm)
overdisp_fun(landscape.treatment.lm)
#plot(landscape.treatment.lm)
summary(landscape.treatment.lm)
anova(landscape.treatment.lm)
car::Anova(landscape.treatment.lm, type=2)
car::Anova(landscape.treatment.lm, type=3)
visreg(landscape.treatment.lm, scale="response")



library(gsg)
library(mgcv)

landscape.gam <- gam(gall_survival_mean ~ s(c.Gall_Height_mm_mean) + s(c.gall_individuals_mean.log) + s(c.Density_per_100_shoots.sqrt), 
                     data = treatment_tree)
landscape.gam <- gam(gall_survival_mean ~ s(c.Gall_Height_mm_mean) + s(c.gall_individuals_mean.log) + s(c.Density_per_100_shoots.sqrt), 
                     data = control_tree)
summary(landscape.gam)
plot(landscape.gam, seWithMean=T, shift = mean(predict(landscape.gam)))
#plot(landscape.gam, seWithMean = T, shift = mean(predict(landscape.gam)), trans = function(x) {exp(x)/(1+exp(x))})

gall_selection.df %>%
  group_by(Treatment.focus) %>%
  summarise(mean.fitness = mean(gall_survival))
gall_selection.df$relative.fit <- with(gall_selection.df, ifelse(Treatment.focus=="Control", gall_survival/0.455, gall_survival/0.599))
chamber <- lm(relative.fit ~ Density_per_100_shoots + gall_individuals + I(Gall_Height_mm-mean(Gall_Height_mm)) + I(Gall_Height_mm-mean(Gall_Height_mm)):Treatment.focus, data=gall_selection.df)
summary(chamber)
anova(chamber)
visreg(chamber, xvar="Gall_Height_mm", by="Treatment.focus")

gall.level <- gall_selection.df %>%
  group_by(Genotype, Treatment.focus, Plant_Position, Gall_Number) %>%
  summarise(gall_individuals = mean(gall_individuals), gall_survival_mean = mean(gall_survival), sample_size = n(), Gall_Height_mm_mean = mean(Gall_Height_mm), Density_per_100_shoots = mean(Density_per_100_shoots)) 
gall.level %>%
  group_by(Treatment.focus) %>%
  summarise(mean.fitness = mean(gall_survival_mean))
gall.level$relative.fit <- with(gall.level, ifelse(Treatment.focus=="Control", gall_survival_mean/0.409, gall_survival_mean/0.611))
gall <- lm(relative.fit ~ 
             I(gall_individuals-mean(gall_individuals)) + I(gall_individuals-mean(gall_individuals)):Treatment.focus +
             I(Gall_Height_mm_mean - mean(Gall_Height_mm_mean)) + I(Density_per_100_shoots - mean(Density_per_100_shoots)), data=gall.level)
summary(gall)
visreg(gall, xvar="gall_individuals", by="Treatment.focus")

 
all_tree %>%
  group_by(Treatment.focus) %>%
  summarise(mean.fitness = mean(gall_survival_mean))
tree.level$relative.fit <- with(tree.level, ifelse(Treatment.focus=="Control", gall_survival_mean/0.415, gall_survival_mean/0.546))
tree <- lm(relative.fit ~ 
             c.Density_per_100_shoots + c.Density_per_100_shoots:Treatment.focus +
             c.gall_individuals_mean.log + c.Gall_Height_mm_mean, data=all_tree)
summary(tree)
visreg(tree, xvar="c.Density_per_100_shoots", by="Treatment.focus", scale="response")




#########################################################################################


ggplot(all_tree, aes(x=c.Gall_Height_mm_mean, y=gall_survival_mean, color=Treatment.focus)) +
  geom_point(aes(size=individual_sample_size), shape=1) +
  geom_smooth(method="lm")

ggplot(all_tree, aes(x=c.gall_individuals_mean.log, y=gall_survival_mean, color=Treatment.focus, weight=individual_sample_size)) +
  geom_point(aes(size=individual_sample_size), shape=1) +
  geom_smooth(method="lm")

ggplot(all_tree, aes(x=c.Density_per_100_shoots, y=gall_survival_mean, color=Treatment.focus, weight=individual_sample_size)) +
  geom_point(aes(size=individual_sample_size), shape=1) +
  geom_smooth(method="glm", method.args=list(family="binomial"))

test2 <- glm(
  gall_survival_mean ~ (c.Gall_Height_mm_mean*c.gall_individuals_mean.log + c.Density_per_100_shoots)*Treatment.focus,
  data=all_tree, weights = individual_sample_size, contrasts=list(Treatment.focus="contr.sum")
)
coef(test2)
W.exp <- predict(test2, type="response")
mean(W.exp*(1-W.exp))*coef(test2)       # is this how you calculate it?
plot(test2)
r.squaredGLMM(test2)
summary(test2)

test <- lmer(
  log(gall_survival_mean+0.05) ~ (c.Gall_Height_mm_mean*c.gall_individuals_mean.log + c.Density_per_100_shoots)*Treatment.focus + 
    (Treatment.focus|Genotype),
  data=all_tree,  contrasts=list(Treatment.focus="contr.sum")
) 

test2 <- lmer(
  log(gall_survival_mean+0.05) ~ (c.Gall_Height_mm_mean*c.gall_individuals_mean.log + c.Density_per_100_shoots) + 
    (1|Genotype),
  data=all_tree,  contrasts=list(Treatment.focus="contr.sum")
) 

model.sel(list(test,test2))

plot(test)
summary(test)
overdisp_fun(test)
library(MuMIn)
r.squaredGLMM(test)
car::vif(test)
car::Anova(test, type=3)

visreg(test, xvar="c.Gall_Height_mm_mean", overlay=TRUE, by="Treatment.focus", scale="response")
visreg(test, xvar="c.gall_individuals_mean.log", overlay=TRUE, by="Treatment.focus", scale="response")
visreg(test, xvar="c.Density_per_100_shoots", overlay=TRUE, by="Treatment.focus", scale="response")

visreg2d(test, xvar="c.Gall_Height_mm_mean", yvar="c.gall_individuals_mean.log", 
         scale="response", cond=list(Treatment.focus="Control"))
visreg2d(test, xvar="c.Gall_Height_mm_mean", yvar="c.gall_individuals_mean.log", 
         scale="response", cond=list(Treatment.focus="Ectoparasitoid exclusion"))

visreg(test, xvar="c.Gall_Height_mm_mean", by="c.gall_individuals_mean.log", 
         scale="response", cond=list(Treatment.focus="Control"), overlay=TRUE, gg=TRUE)
visreg(test, xvar="c.Gall_Height_mm_mean", by="c.gall_individuals_mean.log", 
         scale="response", cond=list(Treatment.focus="Ectoparasitoid exclusion"), overlay=TRUE, gg=TRUE)
