## LOAD REQUIRED LIBRARIES & SET OPTIONS ----
library(tidyverse)
library(cowplot) # pretty default ggplots
#library(broom) # for tidying multiple linear models
#library(mgcv) # generalized additive models
#library(gsg) # calculating selection gradients
library(visreg)
library(viridis)

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
  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_individuals, gall_survival), funs(mean, n())) %>%
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
  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_individuals, gall_survival), funs(mean, n())) %>%
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
  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_individuals, gall_survival), funs(mean, n())) %>%
  ungroup() %>%
  select(-Density_per_100_shoots_n, -Gall_Height_mm_n) %>%
  mutate(c.Gall_Height_mm_mean = Gall_Height_mm_mean - mean(Gall_Height_mm_mean),
         c.gall_individuals_mean = gall_individuals_mean - mean(gall_individuals_mean),
         c.Density_per_100_shoots_mean = Density_per_100_shoots_mean - mean(Density_per_100_shoots_mean))

# the way I initially did it which I don't think is quite accurate
#control_tree <- control_df %>%
#  group_by(Genotype, Plant_Position) %>%
#  summarise_at(vars(Density_per_100_shoots, Gall_Height_mm, gall_survival), funs(mean, n())) %>%
#  ungroup() %>%
#  select(-Density_per_100_shoots_n, -Gall_Height_mm_n) %>%
#  mutate(c.Gall_Height_mm_mean = Gall_Height_mm_mean - mean(Gall_Height_mm_mean),
#         c.Density_per_100_shoots_mean = Density_per_100_shoots_mean - mean(Density_per_100_shoots_mean))

#### CONTROL - ECO-EVO FITNESS LANDSCAPE ANALYSES ----

## Summary: In the complex food web, trait (evo) gradient is much more important than the 
## abundance (eco) gradient in shaping the fitness landscape.

# analysis
mcpeek_control.glm <- glm(gall_survival_mean ~ c.Gall_Height_mm_mean + c.Density_per_100_shoots_mean + c.gall_individuals_mean + Genotype, 
                         data = control_tree, weights = gall_survival_n, family=binomial(link="logit"))
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

## Summary: In the simple food web, trait (evo) gradient and abundance (eco) gradient are similarly
## important in shaping the fitness landscape. 

# analysis
mcpeek_treatment.glm <- glm(gall_survival_mean ~ c.Gall_Height_mm_mean*c.Density_per_100_shoots_mean*c.gall_individuals_mean+Genotype, 
                           data = treatment_tree, weights = gall_survival_n, family=binomial(link="logit"))
car::vif(mcpeek_treatment.glm)
summary(mcpeek_treatment.glm)
car::Anova(mcpeek_treatment.glm, type=3) # potential 3-way interaction...

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

mcpeek_all.geno.glm <- glm(gall_survival_mean ~ c.Gall_Height_mm_mean + c.Density_per_100_shoots_mean*Treatment.focus + 
                             Genotype + Treatment.focus, 
                           data = all_tree, weights = gall_survival_n, family=quasibinomial(), #binomial(link="logit"),
                           contrasts = list(Treatment.focus="contr.sum", Genotype="contr.sum"))
summary(mcpeek_all.geno.glm)
car::vif(mcpeek_all.geno.glm) 
car::Anova(mcpeek_all.geno.glm, type=3)

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
