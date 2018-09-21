
## LOAD REQUIRED LIBRARIES & SET OPTIONS ----
library(mgcv)
library(gsg)
library(tidyverse)
library(cowplot) # pretty default ggplots

## LOAD & MANAGE DATASET ----

gall_selection.df <- read_csv("gall_selection_data.csv") %>%
  mutate(Plant_Position = as.character(Plant_Position),
         Gall_Number = as.character(Gall_Number),
         Treatment.focus = as.character(Treatment.focus),
         Treatment = as.character(Treatment),
         Genotype = as.character(Genotype),
         Gall_Letter = as.character(Gall_Letter)) %>%
  unite(Gall_ID, Gall_Number, Gall_Letter, remove = FALSE) %>%
  
  # subset data for analysis
  filter(phenology == "early", Location == "tree",
         platy > 0 | ectos > 0 | pupa > 0) %>%                  # eliminate unknown sources of mortality
  mutate(gall_survival = as.numeric(ifelse(pupa > 0, 1, 0)),
         egg_parasitoid = as.numeric(ifelse(platy > 0, 1, 0)),
         sc.Gall_Height_mm = as.numeric(scale(Gall_Height_mm)), 
         sc.gall_individuals = as.numeric(scale(gall_individuals)), 
         sc.Density_per_100_shoots = as.numeric(scale(Density_per_100_shoots))) 

# aggregate data at gall-level

gall_agg.df <- gall_selection.df %>% 
  group_by(Treatment.focus, Gall_Number) %>%
  summarise(prop_survival = mean(gall_survival), # effectively a proportion
            gall_survival = sum(gall_survival),
            gall_trials = n(),
            max.Gall_Height_mm = max(Gall_Height_mm), # taking max to be consistent with prior work
            gall_individuals = mean(gall_individuals), # all the same number, chose mean for convenience
            Density_per_100_shoots = mean(Density_per_100_shoots)) # all the same, chose mean for convenience

control_agg.df <- filter(gall_agg.df, Treatment.focus == "Control") %>%
  mutate(sc.max.Gall_Height_mm = as.numeric(scale(max.Gall_Height_mm)),
         sc.gall_individuals = as.numeric(scale(gall_individuals)),
         sc.Density_per_100_shoots = as.numeric(scale(Density_per_100_shoots)))
mean.control_survival <- mean(control_agg.df$prop_survival)

treatment_agg.df <- filter(gall_agg.df, Treatment.focus == "Ectoparasitoid exclusion") %>%
  mutate(sc.max.Gall_Height_mm = as.numeric(scale(max.Gall_Height_mm)),
         sc.gall_individuals = as.numeric(scale(gall_individuals)),
         sc.Density_per_100_shoots = as.numeric(scale(Density_per_100_shoots)))
mean.treatment_survival <- mean(treatment_agg.df$prop_survival)

gall_agg.df <- gall_agg.df %>%
  mutate(rel_fitness = ifelse(Treatment.focus == "Control", 
                              prop_survival - mean.control_survival, 
                              prop_survival - mean.treatment_survival),
         sc.max.Gall_Height_mm = as.numeric(scale(max.Gall_Height_mm)),
         sc.gall_individuals = as.numeric(scale(gall_individuals)),
         sc.Density_per_100_shoots = as.numeric(scale(Density_per_100_shoots)))

control.gam <- gam(prop_survival ~ sc.max.Gall_Height_mm + I(sc.max.Gall_Height_mm^2) + sc.gall_individuals + I(sc.gall_individuals^2) + sc.Density_per_100_shoots + I(sc.Density_per_100_shoots^2) + 
                     sc.max.Gall_Height_mm:sc.gall_individuals + sc.max.Gall_Height_mm:sc.Density_per_100_shoots + sc.gall_individuals:sc.Density_per_100_shoots,
                   data = control_agg.df,
                   family = binomial(link = "logit"), weights = gall_trials)
summary(control.gam)
plot(control.gam, se = 1, seWithMean = TRUE, rug = FALSE, shift = mean(predict(control.gam)), trans = function(x){exp(x)/(1+exp(x))})  
gradient.calc(mod = control.gam, 
              phenotype = "max.Gall_Height_mm", 
              covariates = c("gall_individuals", "Density_per_100_shoots"))$ests
gradient.calc(mod = control.gam, 
              phenotype = "gall_individuals", 
              covariates = c("max.Gall_Height_mm", "Density_per_100_shoots"))$ests
gradient.calc(mod = control.gam, 
              phenotype = "Density_per_100_shoots", 
              covariates = c("max.Gall_Height_mm", "gall_individuals"))$ests

treatment.gam <- gam(prop_survival ~ s(max.Gall_Height_mm) + s(gall_individuals, k=9) + s(Density_per_100_shoots) + 
                     max.Gall_Height_mm:gall_individuals + max.Gall_Height_mm:Density_per_100_shoots + gall_individuals:Density_per_100_shoots,
                   data = treatment_agg.df,
                   family = binomial(link = "logit"), weights = gall_trials)
summary(treatment.gam)
plot(treatment.gam, se = 1, seWithMean = TRUE, rug = FALSE, shift = mean(predict(treatment.gam)), trans = function(x){exp(x)/(1+exp(x))})  
gradient.calc(mod = treatment.gam, 
              phenotype = "max.Gall_Height_mm", 
              covariates = c("gall_individuals", "Density_per_100_shoots"))$ests
gradient.calc(mod = treatment.gam, 
              phenotype = "gall_individuals", 
              covariates = c("max.Gall_Height_mm", "Density_per_100_shoots"))$ests
gradient.calc(mod = treatment.gam, 
              phenotype = "Density_per_100_shoots", 
              covariates = c("max.Gall_Height_mm", "gall_individuals"))$ests

library(visreg)
control.lm <- lm(prop_survival ~ sc.max.Gall_Height_mm + I(sc.max.Gall_Height_mm^2) + 
                     sc.gall_individuals + I(sc.gall_individuals^2) + 
                     sc.Density_per_100_shoots + I(sc.Density_per_100_shoots^2) + 
                     sc.max.Gall_Height_mm:sc.gall_individuals + 
                   sc.max.Gall_Height_mm:sc.Density_per_100_shoots + 
                   sc.gall_individuals:sc.Density_per_100_shoots,
                   data = control_agg.df)
summary(control.lm)
plot(control.lm)
visreg(control.lm)

treatment.lm <- lm(prop_survival ~ sc.max.Gall_Height_mm + I(sc.max.Gall_Height_mm^2) + 
                     sc.gall_individuals + I(sc.gall_individuals^2) + 
                     sc.Density_per_100_shoots + I(sc.Density_per_100_shoots^2) +
                   sc.max.Gall_Height_mm:sc.gall_individuals + 
                     sc.max.Gall_Height_mm:sc.Density_per_100_shoots + 
                     sc.gall_individuals:sc.Density_per_100_shoots,
                 data = treatment_agg.df)
summary(treatment.lm)
plot(treatment.lm)
visreg(treatment.lm)

all.lm <- lm(rel_fitness ~ sc.max.Gall_Height_mm + sc.max.Gall_Height_mm:Treatment.focus +
               I(sc.max.Gall_Height_mm^2) + #I(sc.max.Gall_Height_mm^2):Treatment.focus +
               sc.gall_individuals + sc.gall_individuals:Treatment.focus +
               I(sc.gall_individuals^2) + #I(sc.gall_individuals^2):Treatment.focus + 
               sc.Density_per_100_shoots + sc.Density_per_100_shoots:Treatment.focus +
               I(sc.Density_per_100_shoots^2) + #I(sc.Density_per_100_shoots^2):Treatment.focus +
               sc.max.Gall_Height_mm:sc.gall_individuals + #sc.max.Gall_Height_mm:sc.gall_individuals:Treatment.focus +
               sc.max.Gall_Height_mm:sc.Density_per_100_shoots + #sc.max.Gall_Height_mm:sc.Density_per_100_shoots:Treatment.focus +
               sc.gall_individuals:sc.Density_per_100_shoots, #+ sc.gall_individuals:sc.Density_per_100_shoots:Treatment.focus,
              data = gall_agg.df)
summary(all.lm)
plot(all.lm)
visreg(all.lm)
car::vif(all.lm)

## at tree-level
tree_agg.df <- gall_selection.df %>% 
  group_by(Genotype, Treatment.focus, Plant_Position) %>%
  summarise(prop_survival = mean(gall_survival), # effectively a proportion
            gall_survival = sum(gall_survival),
            gall_trials = n(),
            mean.Gall_Height_mm = mean(Gall_Height_mm), 
            mean.gall_individuals = mean(gall_individuals),
            Density_per_100_shoots = mean(Density_per_100_shoots)) %>% # all the same, chose mean for convenience
  mutate(sc.Density_per_100_shoots = as.numeric(scale(Density_per_100_shoots)))

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

ggplot(tree_agg.df, aes(x=Density_per_100_shoots, y=prop_survival, size=gall_trials, color=Treatment.focus)) +
  geom_point() +
  binomial_smooth() 

mean.control_tree <- mean(filter(tree_agg.df, Treatment.focus=="Control")$prop_survival)
mean.treatment_tree <- mean(filter(tree_agg.df, Treatment.focus=="Ectoparasitoid exclusion")$prop_survival)
tree_agg.df <- mutate(tree_agg.df, rel_fitness = ifelse(Treatment.focus=="Control", prop_survival, prop_survival - mean.treatment_tree))

tree.glm <- glm(prop_survival ~ Density_per_100_shoots*Treatment.focus, 
                data = tree_agg.df, 
                family = quasibinomial(), #binomial(link="logit"),
                weights = gall_trials)
summary(tree.glm)

tree.lm <- lm(rel_fitness ~ #Treatment.focus + #mean.Gall_Height_mm + mean.gall_individuals + 
                #mean.Gall_Height_mm:Treatment.focus +
                #mean.gall_individuals:Treatment.focus +
                Density_per_100_shoots + Density_per_100_shoots:Treatment.focus, 
                data = tree_agg.df)
summary(tree.lm)
plot(tree.lm)
visreg(tree.lm, xvar="Density_per_100_shoots", by="Treatment.focus")

library(lme4)
tree.lmer <- lmer(prop_survival ~ Treatment.focus*sc.Density_per_100_shoots + (1|Genotype), 
              data = tree_agg.df) #filter(tree_agg.df, Treatment.focus=="Ectoparasitoid exclusion"))
summary(tree.lmer)

tree.lm1 <- lm(prop_survival ~ Treatment.focus*sc.Density_per_100_shoots, 
                  data = tree_agg.df) #filter(tree_agg.df, Treatment.focus=="Ectoparasitoid exclusion"))
summary(tree.lm1)

control.lmer <- lmer(gall_survival ~ sc.Gall_Height_mm + I(sc.Gall_Height_mm^2) +
                         sc.gall_individuals + I(sc.gall_individuals^2) +
                         sc.Density_per_100_shoots + I(sc.Density_per_100_shoots^2) +
                         sc.Gall_Height_mm:sc.gall_individuals + sc.Gall_Height_mm:sc.Density_per_100_shoots + sc.gall_individuals:sc.Density_per_100_shoots +
                         (1|Genotype/Plant_Position/Gall_Number),
                       data = filter(gall_selection.df, Treatment.focus == "Control"))
summary(control.lmer)
round(fixef(control.lmer),2)

treatment.lmer <- lmer(gall_survival ~ sc.Gall_Height_mm + I(sc.Gall_Height_mm^2) +
                    sc.gall_individuals + I(sc.gall_individuals^2) +
                    sc.Density_per_100_shoots + I(sc.Density_per_100_shoots^2) +
                    sc.Gall_Height_mm:sc.gall_individuals + sc.Gall_Height_mm:sc.Density_per_100_shoots + sc.gall_individuals:sc.Density_per_100_shoots +
                    (1|Genotype/Plant_Position/Gall_Number),
                  data = filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion"))
summary(treatment.lmer)
round(fixef(treatment.lmer),2)

cbind(control = round(fixef(control.lmer),2), treatment = round(fixef(treatment.lmer),2))

mean.control_surv <- mean(filter(gall_selection.df, Treatment.focus=="Control")$gall_survival)
mean.treatment_surv <- mean(filter(gall_selection.df, Treatment.focus=="Ectoparasitoid exclusion")$gall_survival)

gall_selection.df <- mutate(gall_selection.df, 
                            rel_fitness = ifelse(Treatment.focus=="Control", gall_survival/mean.control_surv, gall_survival/mean.treatment_surv))

all.lmer <- lme4::lmer(rel_fitness ~ (Treatment.focus + sc.Gall_Height_mm + sc.gall_individuals + sc.Density_per_100_shoots)^2 +
                         I(sc.Gall_Height_mm^2) + I(sc.gall_individuals^2) + I(sc.Density_per_100_shoots^2) +
                   (1|Genotype/Plant_Position/Gall_Number),
                 data = gall_selection.df)
library(lmerTest)
anova(all.lmer)
summary(all.lmer)
plot(all.lmer)
visreg(all.lmer, xvar="sc.Gall_Height_mm", by="Treatment.focus")
visreg(all.lmer, xvar="sc.gall_individuals", by="Treatment.focus")
visreg(all.lmer, xvar="sc.Density_per_100_shoots", by="Treatment.focus")


all.lm <- lm(rel_fitness ~ (Treatment.focus + sc.Gall_Height_mm + sc.gall_individuals + sc.Density_per_100_shoots)^2 +
               I(sc.Gall_Height_mm^2) + I(sc.gall_individuals^2) + I(sc.Density_per_100_shoots^2),
                 data = gall_selection.df)
summary(all.lm)
plot(all.lm)
visreg(all.lm)
visreg(all.lm, xvar="sc.Gall_Height_mm", by="Treatment.focus")
visreg(all.lm, xvar="sc.gall_individuals", by="Treatment.focus")
visreg(all.lm, xvar="sc.Density_per_100_shoots", by="Treatment.focus")
visreg2d(all.lm, xvar="sc.Density_per_100_shoots", yvar="sc.gall_individuals", cond = list(Treatment.focus="Control"))
visreg2d(all.lm, xvar="sc.Density_per_100_shoots", yvar="sc.gall_individuals", cond = list(Treatment.focus="Ectoparasitoid exclusion"))

cbind(lmer = round(fixef(all.lmer),2), lm = round(coef(all.lm),2))

hist(sqrt(gall_selection.df$Gall_Height_mm))
hist(sqrt(gall_selection.df$Density_per_100_shoots))
hist(log(gall_selection.df$gall_individuals))

test.lmer <- lme4::glmer(gall_survival ~ (Treatment.focus + sc.Gall_Height_mm + sc.gall_individuals + sc.Density_per_100_shoots)^2 +
                          (1|Genotype/Plant_Position/Gall_Number),
                        data = gall_selection.df,
                        family = binomial(link=logit), control = glmerControl(optimizer="bobyqa"))
summary(test.lmer)
visreg(test.lmer, xvar="sc.Gall_Height_mm", by="Treatment.focus", scale="response")
visreg(test.lmer, xvar="sc.gall_individuals", by="Treatment.focus", scale="response")
visreg(test.lmer, xvar="sc.Density_per_100_shoots", by="Treatment.focus", scale="response")
