## LOAD REQUIRED LIBRARIES & SET OPTIONS ----
library(tidyverse)
library(mgcv)
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
  mutate(gall_survival = ifelse(pupa > 0, 1, 0),
         log.size = log(Gall_Height_mm))

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
         log1.density = log(Density_per_100_shoots+1))

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

cont.gam <- gamm(gall_survival ~ s(Gall_Height_mm) + s(gall_individuals, k = 9) + s(Density_per_100_shoots), 
                 random = list(Genotype=~1, Plant_Position=~1, Gall_Number=~1),
                 data = control.df, 
                 family = "binomial")
summary(cont.gam$gam)
plot(cont.gam$gam, shift = mean(predict(cont.gam$gam)),
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
