## LOAD REQUIRED LIBRARIES & SET OPTIONS ----
library(tidyverse)
library(cowplot) # pretty default ggplots
#library(broom) # for tidying multiple linear models
#library(mgcv) # generalized additive models
#library(gsg) # calculating selection gradients
library(visreg)
library(viridis)

mean_log <- function(x) mean(log(x))
var_log <- function(x) var(log(x))

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

## CREATE SEPARATE DATASETS FOR SELECTION ANALYSES ----
control_df <- gall_selection.df %>%
  filter(Treatment.focus == "Control")

treatment_df <- gall_selection.df %>%
  filter(Treatment.focus == "Ectoparasitoid exclusion") 

# summarize each data set at tree-level
all_tree <- gall_selection.df %>%
  # summarise first at gall level
  group_by(Genotype, Treatment.focus, Plant_Position, Gall_Number) %>%
  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_individuals, gall_survival), funs(mean)) %>%
  
  # then summarise at the plant level
  group_by(Genotype, Treatment.focus, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_individuals, gall_survival), funs(mean, var, mean_log, var_log, n())) %>%
  ungroup() %>%
  select(-Density_per_100_shoots_n, -Gall_Height_mm_n) %>%
  mutate(c.Gall_Height_mm_mean = Gall_Height_mm_mean - mean(Gall_Height_mm_mean),
         c.gall_individuals_mean = gall_individuals_mean - mean(gall_individuals_mean),
         c.Density_per_100_shoots_mean = Density_per_100_shoots_mean - mean(Density_per_100_shoots_mean))

control_tree <- control_df %>%
  # summarise first at gall level
  group_by(Genotype, Treatment.focus, Plant_Position, Gall_Number) %>%
  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_individuals, gall_survival), funs(mean)) %>%
  
  # then summarise at the plant level
  group_by(Genotype, Treatment.focus, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_individuals, gall_survival), funs(mean, var, n())) %>%
  ungroup() %>%
  select(-Density_per_100_shoots_n, -Gall_Height_mm_n) %>%
  mutate(c.Gall_Height_mm_mean = Gall_Height_mm_mean - mean(Gall_Height_mm_mean),
         c.gall_individuals_mean = gall_individuals_mean - mean(gall_individuals_mean),
         c.Density_per_100_shoots_mean = Density_per_100_shoots_mean - mean(Density_per_100_shoots_mean))

treatment_tree <- treatment_df %>%
  # summarise first at gall level
  group_by(Genotype, Treatment.focus, Plant_Position, Gall_Number) %>%
  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_individuals, gall_survival), funs(mean)) %>%
  
  # then summarise at the plant level
  group_by(Genotype, Treatment.focus, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_individuals, gall_survival), funs(mean, var, n())) %>%
  ungroup() %>%
  select(-Density_per_100_shoots_n, -Gall_Height_mm_n) %>%
  mutate(c.Gall_Height_mm_mean = Gall_Height_mm_mean - mean(Gall_Height_mm_mean),
         c.gall_individuals_mean = gall_individuals_mean - mean(gall_individuals_mean),
         c.Density_per_100_shoots_mean = Density_per_100_shoots_mean - mean(Density_per_100_shoots_mean))


#### TEST THE ASSUMPTION THAT MEAN IS NOT RELATED TO VARIANCE
all_tree %>%
  ggplot(., aes(x=Gall_Height_mm_mean, y=Gall_Height_mm_var)) +
  geom_point() + geom_smooth(method="lm")

# very strong positive correlation...
all_tree %>%
  ggplot(., aes(x=gall_individuals_mean, y=gall_individuals_var)) +
  geom_point() + geom_smooth(method="lm")

all_tree %>%
  ggplot(., aes(x=gall_individuals_mean_log, y=gall_individuals_var_log, size=gall_survival_n, weight=gall_survival_n)) +
  geom_point() + geom_smooth(method="lm")

summary(lm(gall_individuals_var_log ~ gall_individuals_mean_log, data=all_tree))
summary(lm(gall_individuals_var_log ~ gall_individuals_mean_log, data=all_tree, weights = gall_survival_n))

summary(all_tree$Density_per_100_shoots_var)


#### MODEL RATIONALE

## Gall Size * Treatment
# Consistently has a monotonic increase on survival in past studies with a complex community. 
# The strength of this interaction could vary with the predator community (e.g. if we exclude predators causing selection on this trait), and we have evidence of this from the PNAS paper.
# I don't expect other plant traits to interact with gall size to affect survival though.

## Gall Size * Gall Individuals * Treatment
# Currently unkown effects.
# This does represent local density, which we expect to influence insect fitness.
# Also, local density could interact with gall size, by providing more of a refuge to larva that would normally be attacked in smaller galls.
# Likewise, we expect these interactions to potentially depend on the parasitoid community.

## Gall Density * Treatment
# PNAS paper suggested that density could shift the relative dominance of certain parasitoids in the community.
# It's not clear to me whether gall density would interact with gall size or gall individuals. I would expect the effects to be additive.

## Other 3-way or 4-way interactions?
# it is hard to imagine...
library(lme4)
library(piecewiseSEM)
gen.glm <- function(formula, mixed) {
  if(mixed==TRUE){
    glmer(formula, data = all_tree, weights = gall_survival_n, family=binomial(link="logit"), na.action = "na.fail",
          contrasts=list(Genotype="contr.sum",Treatment.focus="contr.sum"), glmerControl(optimizer="bobyqa"))
  } else {
    glm(formula, data = all_tree, weights = gall_survival_n, family=binomial(link="logit"), na.action = "na.fail",
        contrasts=list(Genotype="contr.sum",Treatment.focus="contr.sum"))
  }
  
}
## A pure evolutionary model
## Gall Size
evo <- gen.glm(gall_survival_mean ~ c.Gall_Height_mm_mean+(1|Genotype), mixed=TRUE)
summary(evo)
sem.model.fits(evo)
car::vif(evo)

## A pure ecological model without community context
## Gall Density + Gall Individuals
eco <- gen.glm(gall_survival_mean ~ c.gall_individuals_mean + c.Density_per_100_shoots_mean, mixed=FALSE)
car::vif(eco)
sem.model.fits(eco)

## Community context alone
com <- gen.glm(gall_survival_mean ~ Treatment.focus, mixed=FALSE)
car::vif(com)
sem.model.fits(com)

## A pure ecological model with community context
## (Gall Density + Gall Individuals)*Treatment
eco.com <- gen.glm(gall_survival_mean ~ (c.gall_individuals_mean + c.Density_per_100_shoots_mean)*Treatment.focus, mixed=FALSE)
car::vif(eco.com)
sem.model.fits(eco.com)

## A pure evolutionary model with community context
## Gall size * Treatment
evo.com <- gen.glm(gall_survival_mean ~ c.Gall_Height_mm_mean*Treatment.focus + (1+Treatment.focus|Genotype), mixed=TRUE)
summary(evo.com)
sem.model.fits(evo.com)
r.squaredGLMM(evo.com)

## An eco-evo model without community context
eco.evo <- gen.glm(gall_survival_mean ~ c.Gall_Height_mm_mean*c.gall_individuals_mean + c.Density_per_100_shoots_mean + (1|Genotype), mixed=TRUE)
summary(eco.evo)
r.squaredGLMM(eco.evo)

## An eco-evo model with community context
## (Gall size*Gall Individuals + Gall Density)*Treatment
eco.evo.com <- gen.glm(gall_survival_mean ~ (c.Gall_Height_mm_mean*c.gall_individuals_mean + c.Density_per_100_shoots_mean)*Treatment.focus + (1+Treatment.focus|Genotype), mixed=TRUE)
car::vif(eco.evo.com) # some pretty large values for Genotype and Treatment focus
summary(eco.evo.com)
sem.model.fits(eco.evo.com)
r.squaredGLMM(eco.evo.com)

noRS.eco.evo.com <- gen.glm(gall_survival_mean ~ (c.Gall_Height_mm_mean*c.gall_individuals_mean + c.Density_per_100_shoots_mean)*Treatment.focus + (1|Genotype), mixed=TRUE)

## Compare these models
model.sel(list(evo, eco, com, eco.com, evo.com, eco.evo, eco.evo.com, noRS.eco.evo.com))

## sem.model.fits not working with random slope

sem.model.fits(list(evo, eco, com, eco.com, evo.com))
# there is incredibly strong support for an interaction between ecology, evolution, and community context.

ALT.eco.evo.com <- gen.glm(gall_survival_mean ~ (c.Gall_Height_mm_mean*c.gall_individuals_mean + c.Density_per_100_shoots_mean)*Treatment.focus+Genotype)
AICc(eco.evo.com, ALT.eco.evo.com)
summary(ALT.eco.evo.com)
car::Anova(ALT.eco.evo.com, type=3)
car::vif(ALT.eco.evo.com)

model.sel(list(evo, eco, com, eco.com, evo.com, eco.evo.com, ALT.eco.evo.com))

#### FULL - ECO-EVO FITNESS LANDSCAPE ANALYSES ----

## Summary: As suggested by the piecewise models above, food-web structure alters the abundance (eco) gradient,
## but not the trait (evo) gradient. 

# analysis
mcpeek_all.glm <- glm(gall_survival_mean ~ (c.Gall_Height_mm_mean+c.gall_individuals_mean+c.Density_per_100_shoots_mean+Genotype)*Treatment.focus, 
                      data = all_tree, weights = gall_survival_n, family=binomial(link="logit"),
                      contrasts = list(Treatment.focus="contr.sum",Treatment.focus="contr.sum"), na.action = "na.fail")

# aic model comparison
library(MuMIn)
dredge.test <- dredge(mcpeek_all.glm) %>%
  data.frame() %>%
  filter(delta < 4) 
dredge.test
sum(dredge.test$weight)
par(mar=c(3,5,10,5))
plot(dredge(mcpeek_all.glm))

test <- dredge(mcpeek_all.glm)
test.avg <- model.avg(test)
coefTable(test.avg)

car::vif(mcpeek_all.glm) 
summary(mcpeek_all.glm) # note rather high residual deviance...
car::Anova(mcpeek_all.glm, type=3) # note that there no longer appears to be a density*treatment interaction, but looking at the density*treatment interaction plot makes much more sense now (survival is always higher in ectoparasitoid treatment)

# visualize models
visreg(mcpeek_all.glm, xvar="c.Gall_Height_mm_mean", by="Treatment.focus", gg=TRUE, overlay=TRUE, scale="response")
visreg(mcpeek_all.glm, xvar="c.Density_per_100_shoots_mean", by="Treatment.focus", gg=TRUE, overlay=TRUE, scale="response")
visreg(mcpeek_all.glm, xvar="c.Gall_Height_mm_mean", by="c.gall_individuals_mean", 
       gg=TRUE, overlay=TRUE, scale="response", cond=list(Treatment.focus="Control"))
visreg(mcpeek_all.glm, xvar="c.Gall_Height_mm_mean", by="c.gall_individuals_mean", 
       gg=TRUE, overlay=TRUE, scale="response", cond=list(Treatment.focus="Ectoparasitoid exclusion"))

# standardized effect sizes
coef(mcpeek_all.glm)["c.Gall_Height_mm_mean"] * sd(all_tree$Gall_Height_mm_mean)
coef(mcpeek_all.glm)["c.Density_per_100_shoots_mean"] * sd(all_tree$Density_per_100_shoots_mean)
coef(mcpeek_all.glm)["Treatment.focus1"] * sd(as.numeric(as.factor(all_tree$Treatment.focus))) * 2 # is it correct to multiply this by 2? I'm doing this because the standard deviation for a binary variable is 0.5...need to reread Schielzeth 2010
coef(mcpeek_all.glm)["c.Density_per_100_shoots_mean:Treatment.focus1"] * sd(treatment_tree$Density_per_100_shoots_mean * as.numeric(as.factor(all_tree$Treatment.focus))) # need to double-check whether this makes sense...


best_all.glm <- glm(gall_survival_mean ~ c.Gall_Height_mm_mean*c.gall_individuals_mean+c.Density_per_100_shoots_mean + Treatment.focus, 
                    data = all_tree, weights = gall_survival_n, family=binomial(link="logit"),
                    contrasts = list(Treatment.focus="contr.sum",Treatment.focus="contr.sum"), na.action = "na.fail")
summary(best_all.glm)
visreg(best_all.glm, xvar="c.Gall_Height_mm_mean", by="c.gall_individuals_mean", 
       gg=TRUE, overlay=TRUE, scale="response", cond=list(Treatment.focus="Control"))
visreg(best_all.glm, xvar="c.Gall_Height_mm_mean", by="c.gall_individuals_mean", 
       gg=TRUE, overlay=TRUE, scale="response", cond=list(Treatment.focus="Ectoparasitoid exclusion"))


# the way I initially did it which I don't think is quite accurate
#control_tree <- control_df %>%
#  group_by(Genotype, Plant_Position) %>%
#  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_survival), funs(mean, n())) %>%
#  ungroup() %>%
#  select(-Density_per_100_shoots_n, -Gall_Height_mm_n) %>%
#  mutate(c.Gall_Height_mm_mean = Gall_Height_mm_mean - mean(Gall_Height_mm_mean),
#         c.Density_per_100_shoots_mean = Density_per_100_shoots_mean - mean(Density_per_100_shoots_mean))

#### CAN I TEST FOR GENOTYPE-BY-TRAIT INTERACTION?
geno_mean_traits_control <- control_tree %>%
  group_by(Genotype) %>%
  summarise_at(vars(Density_per_100_shoots_mean, Gall_Height_mm_mean, gall_individuals_mean), funs(mean)) %>%
  rename(geno_mean_density = Density_per_100_shoots_mean, 
         geno_mean_height = Gall_Height_mm_mean, 
         geno_mean_indiv = gall_individuals_mean)

GxTrait_control <- control_tree %>%
  left_join(., geno_mean_traits_control) %>%
  mutate(geno_mean_height_dev = Gall_Height_mm_mean - geno_mean_height,
         geno_mean_ind_dev = gall_individuals_mean - geno_mean_indiv)

library(lmerTest)
GxTrait.lm <- lmer(gall_survival_mean ~ geno_mean_indiv + geno_mean_ind_dev + (geno_mean_indiv | Genotype), 
                  data = GxTrait_control)
#car::vif(GxTrait.lm)
#alias(GxTrait.lm)
summary(GxTrait.lm)
anova(GxTrait.lm)
car::Anova(GxTrait.lm)
car::Anova(GxTrait.lm, type=3)

visreg(GxTrait.lm)

GxTrait.lm2 <- lmer(gall_survival_mean ~ geno_mean_indiv + geno_mean_ind_dev + (1 | Genotype),
                   data = GxTrait_control)
summary(GxTrait.lm2)
0.007393/(0.007393+0.048812)
AIC(GxTrait.lm, GxTrait.lm2)

visreg(GxTrait.lm2)

#### CONTROL - ECO-EVO FITNESS LANDSCAPE ANALYSES ----

## Summary: In the complex food web, trait (evo) gradient is much more important than the 
## abundance (eco) gradient in shaping the fitness landscape.

# analysis
mcpeek_control.glm <- glm(gall_survival_mean ~ c.Gall_Height_mm_mean*c.Density_per_100_shoots_mean*c.gall_individuals_mean + Genotype, 
                         data = control_tree, weights = gall_survival_n, family=binomial(link="logit"), contrasts = list(Genotype="contr.sum"))
car::vif(mcpeek_control.glm)
summary(mcpeek_control.glm)
car::Anova(mcpeek_control.glm, type=3)

# visualize models
visreg(mcpeek_control.glm, xvar="c.Gall_Height_mm_mean", by="c.Density_per_100_shoots_mean", gg=TRUE, overlay=TRUE, scale="response")
visreg(mcpeek_control.glm, xvar="c.Density_per_100_shoots_mean", gg=TRUE, scale="response")
visreg(mcpeek_control.glm, xvar="c.gall_individuals_mean", gg=TRUE, scale="response")
visreg2d(mcpeek_control.glm, xvar="c.Gall_Height_mm_mean", yvar="c.Density_per_100_shoots_mean", 
         scale="response", plot.type = "persp")

# standardized effect sizes
coef(mcpeek_control.glm)["c.Gall_Height_mm_mean"] * sd(control_tree$Gall_Height_mm_mean)
coef(mcpeek_control.glm)["c.Density_per_100_shoots_mean"] * sd(control_tree$Density_per_100_shoots_mean)

#### TREATMENT - ECO-EVO FITNESS LANDSCAPE ANALYSES ----

test_treatment_tree <- mutate(treatment_tree,
                              c.log.gall_individuals_mean = log(gall_individuals_mean)-mean(log(gall_individuals_mean)),
                              c.log1.Density_per_100_shoots_mean = log(Density_per_100_shoots_mean+1) - mean(log(Density_per_100_shoots_mean+1)))
hist(treatment_tree$Gall_Height_mm_mean)
hist(log(treatment_tree$gall_individuals_mean))
hist(log(treatment_tree$Density_per_100_shoots_mean+0.8))
hist(sqrt(treatment_tree$Density_per_100_shoots_mean))

hist(treatment_tree$Density_per_100_shoots_mean)
hist(log(control_tree$Density_per_100_shoots_mean)) 

table(treatment_tree$Density_per_100_shoots_mean)

log(treatment_tree$Density_per_100_shoots_mean+1)
log(treatment_tree$gall_individuals_mean)

## Summary: In the simple food web, trait (evo) gradient and abundance (eco) gradient are similarly
## important in shaping the fitness landscape. 

# analysis
mcpeek_treatment.glm <- glm(gall_survival_mean ~ c.Gall_Height_mm_mean*c.log1.Density_per_100_shoots_mean*c.log.gall_individuals_mean, 
                           data = test_treatment_tree, contrasts=list(Genotype="contr.sum"), weights = gall_survival_n, family=binomial(link="logit"))
car::vif(mcpeek_treatment.glm)
summary(mcpeek_treatment.glm)
car::Anova(mcpeek_treatment.glm, type=3) # potential 3-way interaction...

visreg(mcpeek_treatment.glm, xvar="c.Gall_Height_mm_mean")
visreg(mcpeek_treatment.glm, xvar="c.Gall_Height_mm_mean", by="c.log.gall_individuals_mean", gg=TRUE)
visreg(mcpeek_treatment.glm, xvar="c.Gall_Height_mm_mean", gg=TRUE)
visreg(mcpeek_treatment.glm, xvar="c.gall_individuals_mean", by="c.Gall_Height_mm_mean", gg=TRUE)
visreg(mcpeek_treatment.glm, xvar="c.gall_individuals_mean", by="c.Density_per_100_shoots_mean", gg=TRUE)
visreg2d(mcpeek_treatment.glm, 
         xvar="c.gall_individuals_mean", yvar="c.Density_per_100_shoots_mean",
         cond=list(c.Gall_Height_mm_mean=0),
         scale="response")

# below average gall size, there is a total reversal in the effects of gall indiv and density
visreg2d(mcpeek_treatment.glm, 
         xvar="c.gall_individuals_mean", yvar="c.Density_per_100_shoots_mean",
         cond=list(c.Gall_Height_mm_mean=-2),
         scale="response")

AIC(mcpeek_treatment.glm, update(mcpeek_treatment.glm, .~. -c.Gall_Height_mm_mean:c.Density_per_100_shoots_mean:c.gall_individuals_mean))

library(vegan)
plot(rda(select(treatment_tree, c.Gall_Height_mm_mean,c.Density_per_100_shoots_mean,c.gall_individuals_mean)))#, data=treatment_tree))

# visualize models
visreg(mcpeek_treatment.glm, xvar="c.Gall_Height_mm_mean", by="c.Density_per_100_shoots_mean", gg=TRUE, overlay=TRUE, scale="response")
visreg(mcpeek_treatment.glm, xvar="c.gall_individuals_mean", by="c.Density_per_100_shoots_mean", gg=TRUE, overlay=TRUE, scale="response")
visreg2d(mcpeek_treatment.glm, xvar="c.Gall_Height_mm_mean", yvar="c.Density_per_100_shoots_mean", 
         scale="response", plot.type = "persp")

# standardized effect sizes
coef(mcpeek_treatment.glm)["c.Gall_Height_mm_mean"] * sd(treatment_tree$Gall_Height_mm_mean)
coef(mcpeek_treatment.glm)["c.Density_per_100_shoots_mean"] * sd(treatment_tree$Density_per_100_shoots_mean)


#### FULL - ECO-EVO FITNESS LANDSCAPE ANALYSES ----

## Summary: As suggested by the piecewise models above, food-web structure alters the abundance (eco) gradient,
## but not the trait (evo) gradient. 

# analysis
mcpeek_all.glm <- glm(gall_survival_mean ~ c.Gall_Height_mm_mean*c.Density_per_100_shoots_mean*c.gall_individuals_mean*Treatment.focus, 
                            data = all_tree, weights = gall_survival_n, family=binomial(link="logit"),
                            contrasts = list(Treatment.focus="contr.sum"))
car::vif(mcpeek_all.glm) 
summary(mcpeek_all.glm) # note rather high residual deviance...
car::Anova(mcpeek_all.glm, type=3) # note that there no longer appears to be a density*treatment interaction, but looking at the density*treatment interaction plot makes much more sense now (survival is always higher in ectoparasitoid treatment)

# visualize models
visreg(mcpeek_all.glm, xvar="c.Gall_Height_mm_mean", by="Treatment.focus", overlay=TRUE, scale="response")
visreg(mcpeek_all.glm, xvar="c.Density_per_100_shoots_mean", by="Treatment.focus", overlay=TRUE, scale="response")

# standardized effect sizes
coef(mcpeek_all.glm)["c.Gall_Height_mm_mean"] * sd(all_tree$Gall_Height_mm_mean)
coef(mcpeek_all.glm)["c.Density_per_100_shoots_mean"] * sd(all_tree$Density_per_100_shoots_mean)
coef(mcpeek_all.glm)["Treatment.focus1"] * sd(as.numeric(as.factor(all_tree$Treatment.focus))) * 2 # is it correct to multiply this by 2? I'm doing this because the standard deviation for a binary variable is 0.5...need to reread Schielzeth 2010
coef(mcpeek_all.glm)["c.Density_per_100_shoots_mean:Treatment.focus1"] * sd(treatment_tree$Density_per_100_shoots_mean * as.numeric(as.factor(all_tree$Treatment.focus))) # need to double-check whether this makes sense...


#### WITH GENOTYPE ....
all_tree$test_density <- all_tree$c.Density_per_100_shoots_mean * all_tree$c.gall_individuals_mean

mcpeek_all.geno.glm <- glm(gall_survival_mean ~ c.Gall_Height_mm_mean*c.gall_individuals_mean*c.Density_per_100_shoots_mean*Treatment.focus, 
                           data = all_tree, weights = gall_survival_n, binomial(link="logit"),
                           contrasts = list(Treatment.focus="contr.sum", Genotype="contr.sum"))
summary(mcpeek_all.geno.glm)
car::vif(mcpeek_all.geno.glm) 
car::Anova(mcpeek_all.geno.glm, type=3)

visreg(mcpeek_all.geno.glm, xvar="c.Gall_Height_mm_mean", by="test_density", gg=TRUE, overlay=TRUE, scale="response", cond=list(Treatment.focus="Control"))
visreg(mcpeek_all.geno.glm, xvar="c.Gall_Height_mm_mean", by="test_density", gg=TRUE, overlay=TRUE, scale="response", cond=list(Treatment.focus="Ectoparasitoid exclusion"))

visreg(mcpeek_all.geno.glm, xvar="c.Density_per_100_shoots_mean", by="Treatment.focus", gg=TRUE, overlay=TRUE, scale="response")
visreg(mcpeek_all.geno.glm, xvar="c.gall_individuals_mean", by="Treatment.focus", gg=TRUE, overlay=TRUE, scale="response")

visreg2d(mcpeek_all.geno.glm, xvar="c.gall_individuals_mean", yvar="c.Density_per_100_shoots_mean", cond=list(Treatment.focus="Control"), scale="response")
visreg2d(mcpeek_all.geno.glm, xvar="c.gall_individuals_mean", yvar="c.Density_per_100_shoots_mean", cond=list(Treatment.focus="Ectoparasitoid exclusion"), scale="response")


# note that when I remove GxE, we see that the weird signal where ectoparasitoid survival is less than in the control goes away, which biologiclaly makes the most sense.
visreg(mcpeek_all.geno.glm, xvar="c.Gall_Height_mm_mean", by="Treatment.focus", gg=TRUE, overlay=TRUE, scale="response")
visreg(mcpeek_all.geno.glm, xvar="c.Density_per_100_shoots_mean", by="Treatment.focus", gg=TRUE, overlay=TRUE, scale="response") # hmmm, the fact that ectoparasitoid survival is lower is suspicious to me...Is this a weird genotype effect?
#visreg(mcpeek_all.geno.glm, xvar="Treatment.focus", by="Genotype", gg=TRUE, scale="response")

# exploring basic GxE effects.
# doesn't suggest a GxE effect on gall survival, but strong and independent G and E effects.
mcpeek_GxE.glm <- glm(gall_survival_mean ~ Genotype*Treatment.focus, 
                           data = all_tree, weights = gall_survival_n, family=quasibinomial(), #binomial(link="logit"),
                           contrasts = list(Treatment.focus="contr.sum", Genotype="contr.sum"))
summary(mcpeek_GxE.glm)
car::vif(mcpeek_GxE.glm) 
car::Anova(mcpeek_GxE.glm, type=3)

## exploring
mcpeek_explore.glm <- glm(gall_survival_mean ~ c.Gall_Height_mm_mean*c.Density_per_100_shoots_mean*Treatment.focus*Genotype, 
                           data = all_tree, weights = gall_survival_n, family=binomial(link="logit"),
                           contrasts = list(Treatment.focus="contr.sum", Genotype="contr.sum"))
summary(mcpeek_explore.glm)
car::vif(mcpeek_explore.glm) # worried about incredibly high VIF score
car::Anova(mcpeek_explore.glm, type=3)

test_explore.glm <- update(mcpeek_explore.glm, .~. -c.Gall_Height_mm_mean:c.Density_per_100_shoots_mean:Treatment.focus:Genotype -
                             c.Density_per_100_shoots_mean:Treatment.focus:Genotype - c.Gall_Height_mm_mean:Treatment.focus:Genotype -
                             c.Gall_Height_mm_mean:c.Density_per_100_shoots_mean:Genotype -c.Gall_Height_mm_mean:c.Density_per_100_shoots_mean:Treatment.focus)
car::vif(test_explore.glm)
car::Anova(test_explore.glm, type=3)
summary(test_explore.glm)

## OLD VISUALIZATIONS ----
control_tree %>%
  ggplot(., aes(x = Density_per_100_shoots_mean, y=log(gall_survival_mean), size = gall_survival_n)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x) +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ x, color="red")

control_tree %>%
  filter(Density_per_100_shoots_mean < 75) %>% # excluding the outlying density point reveals no relationship between density and mean fitness
  ggplot(., aes(x = Density_per_100_shoots_mean, y=log(gall_survival_mean), size = gall_survival_n)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x) +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ x, color="red")

control_tree %>%
  ggplot(., aes(x = Gall_Height_mm_mean, y=log(gall_survival_mean), size = gall_survival_n)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x) +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ x, color="red") + 
  geom_smooth(method="lm", formula = y ~ poly(x,2), color="orange") +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ poly(x,2), color="purple")

control_tree %>%
  ggplot(., aes(x = Density_per_100_shoots_mean, y=Gall_Height_mm_mean, fill=log(gall_survival_mean), size = gall_survival_n)) +
  geom_point(shape = 21) 

mcpeek_control.lm <- lm(log(gall_survival_mean+0.1) ~ Gall_Height_mm_mean*Density_per_100_shoots_mean, data = control_tree)
summary(mcpeek_control.lm)
anova(mcpeek_control.lm)

## TREATMENT PLANTS


treatment_tree %>%
  ggplot(., aes(x = Density_per_100_shoots_mean, y=log(gall_survival_mean), size = gall_survival_n)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x) +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ x, color="red")

treatment_tree %>%
  ggplot(., aes(x = Gall_Height_mm_mean, y=log(gall_survival_mean), size = gall_survival_n)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x) +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ x, color="red") + 
  geom_smooth(method="lm", formula = y ~ poly(x,2), color="orange") +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ poly(x,2), color="purple")

treatment_tree %>%
  ggplot(., aes(x = Density_per_100_shoots_mean, y=Gall_Height_mm_mean, fill=log(gall_survival_mean), size = gall_survival_n)) +
  geom_point(shape = 21) 




## FULL STATISTICAL ANALYSIS
all_tree <- gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_survival), funs(mean, n())) %>%
  ungroup() %>%
  select(-Density_per_100_shoots_n, -Gall_Height_mm_n) %>%
  mutate(c.Gall_Height_mm_mean = Gall_Height_mm_mean - mean(Gall_Height_mm_mean),
         c.Density_per_100_shoots_mean = Density_per_100_shoots_mean - mean(Density_per_100_shoots_mean))

all_tree %>%
  ggplot(., aes(x = Density_per_100_shoots_mean, y=log(gall_survival_mean+0.1), size = gall_survival_n)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x) +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ x, color="red") +
  facet_wrap(~Treatment.focus)

all_tree %>%
  ggplot(., aes(x = Gall_Height_mm_mean, y=log(gall_survival_mean+0.1), size = gall_survival_n)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x) +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ x, color="red") + 
  facet_wrap(~Treatment.focus) +
  geom_smooth(method="lm", formula = y ~ poly(x,2), color="orange") +
  geom_smooth(aes(weight = gall_survival_n), method="lm", formula = y ~ poly(x,2), color="purple")

mcpeek_all.lm <- lm(log(gall_survival_mean+0.1) ~ c.Gall_Height_mm_mean*c.Density_per_100_shoots_mean*Treatment.focus, 
                    data = all_tree)
summary(mcpeek_all.lm)
anova(mcpeek_all.lm)
car::Anova(mcpeek_all.lm, type=2)
car::Anova(mcpeek_all.lm, type=3)

hist(all_tree$c.Density_per_100_shoots_mean)
hist(all_tree$gall_survival_mean)
table(all_tree$gall_survival_mean)
summary(all_tree$gall_survival_mean)
hist(log(all_tree$gall_survival_mean))
hist(log(all_tree$gall_survival_mean+0.05))
mcpeek_all.wts.lm <- lm(log(gall_survival_mean+0.05) ~ (c.Gall_Height_mm_mean+c.Density_per_100_shoots_mean)*Treatment.focus + c.Gall_Height_mm_mean:c.Density_per_100_shoots_mean, 
                    data = all_tree, weights = gall_survival_n)
car::vif(mcpeek_all.wts.lm)
summary(mcpeek_all.wts.lm)
anova(mcpeek_all.wts.lm)
car::Anova(mcpeek_all.wts.lm, type=2)
car::Anova(mcpeek_all.wts.lm, type=3)

mcpeek_all.wts.geno.lm <- lm(log(gall_survival_mean+0.05) ~ c.Gall_Height_mm_mean*c.Density_per_100_shoots_mean*Treatment.focus*Genotype, 
                        data = all_tree, weights = gall_survival_n)
summary(mcpeek_all.wts.geno.lm)
anova(mcpeek_all.wts.geno.lm)
car::Anova(mcpeek_all.wts.geno.lm, type=2)
car::Anova(mcpeek_all.wts.geno.lm, type=3)
