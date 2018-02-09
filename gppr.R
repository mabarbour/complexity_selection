
## LOAD REQUIRED LIBRARIES & SET OPTIONS ----

library(tidyverse)
#library(mgcv)
#library(gamm4) # for generalized additive mixed models #
library(gsg)
library(cowplot) # pretty default ggplots
library(lme4)

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
         sc.Gall_Height_mm = as.numeric(scale((Gall_Height_mm))),
         sc.gall_individuals = as.numeric(scale((gall_individuals))),
         sc.Density_per_100_shoots = as.numeric(scale((Density_per_100_shoots))))


landscape <- glmer(gall_survival ~ sc.Gall_Height_mm*sc.gall_individuals*sc.Density_per_100_shoots*Treatment.focus +
                     (1|Genotype/Plant_Position/Gall_Number), data = gall_selection.df, family=binomial(link=logit), 
                   control = glmerControl(optimizer = "bobyqa"), contrasts = list(Treatment.focus="contr.sum"))
summary(landscape)
library(visreg)
visreg(landscape)
visreg(landscape, xvar = "sc.gall_individuals", by = "Treatment.focus")
visreg(landscape, xvar = "sc.Density_per_100_shoots", by = "Treatment.focus")

visreg2d(landscape, x="sc.Density_per_100_shoots", y="sc.gall_individuals", print.cond=TRUE, plot.type="image")

## gppr
control.df <- filter(gall_selection.df, Treatment.focus == "Control") %>% as.data.frame()
control.y <- control.df$gall_survival
control.gppr <- gppr(y = "gall_survival", xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"),
                     data = control.df, family="binomial")
control.gppr$ppr
control.gppr$ppr$alpha
control.gppr$ppr$beta

control.df$major_1 <- control.gppr$ppr$alpha[1]*control.df$sc.Gall_Height_mm + control.gppr$ppr$alpha[2]*control.df$sc.gall_individuals + control.gppr$ppr$alpha[3]*control.df$sc.Density_per_100_shoots

control.major <- gamm4::gamm4(gall_survival ~ s(major_1), random = ~(1|Genotype/Plant_Position/Gall_Number), control.df, family=binomial)
summary(control.major$gam)
gamm.plot(control.major, pages=1)
gradient.calc(mod = control.major$gam, phenotype="major_1")$ests

treatment.df <- filter(gall_selection.df, Treatment.focus == "Ectoparasitoid exclusion") %>% as.data.frame()
treatment.gppr <- gppr(y = "gall_survival", xterms = c("sc.Gall_Height_mm","sc.gall_individuals","sc.Density_per_100_shoots"), data = treatment.df, family="poisson")
treatment.gppr$ppr
treatment.gppr$ppr$alpha
treatment.gppr$ppr$beta

treatment.df$major_1 <- treatment.gppr$ppr$alpha[1]*treatment.df$sc.Gall_Height_mm + treatment.gppr$ppr$alpha[2]*treatment.df$sc.gall_individuals + treatment.gppr$ppr$alpha[3]*treatment.df$sc.Density_per_100_shoots

treatment.major <- gamm4::gamm4(gall_survival ~ s(major_1), random = ~(1|Genotype/Plant_Position/Gall_Number), treatment.df, family=binomial)
summary(treatment.major$gam)
gamm.plot(treatment.major, pages=1)
gradient.calc(mod = treatment.major$gam, phenotype="major_1")$ests


pcs <- princomp(x = select(gall_selection.df, sc.Gall_Height_mm, sc.gall_individuals, sc.Density_per_100_shoots), cor=TRUE)
summary(pcs)
loadings(pcs)
pcs$scores
pc.all <- data.frame(gall_selection.df, -1*pcs$scores)

plot(filter(pc.all, Treatment.focus == "Control")$Comp.1, control.df$major_1)
cor.test(filter(pc.all, Treatment.focus == "Control")$Comp.1, control.df$major_1)

plot(filter(pc.all, Treatment.focus == "Control")$Comp.2, control.df$major_1)
cor.test(filter(pc.all, Treatment.focus == "Control")$Comp.2, control.df$major_1)

plot(filter(pc.all, Treatment.focus == "Ectoparasitoid exclusion")$Comp.1, treatment.df$major_1)
cor.test(filter(pc.all, Treatment.focus == "Ectoparasitoid exclusion")$Comp.1, treatment.df$major_1)

plot(filter(pc.all, Treatment.focus == "Ectoparasitoid exclusion")$Comp.2, treatment.df$major_1)
cor.test(filter(pc.all, Treatment.focus == "Ectoparasitoid exclusion")$Comp.2, treatment.df$major_1)

# not working...

cor.comp1 <- lmer(major_1 ~ Comp.1 + (1|Genotype/Plant_Position/Gall_Number), left_join(filter(pc.all, Treatment.focus == "Control"), control.df))


