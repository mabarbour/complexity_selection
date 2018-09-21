
library(RCurl)
library(tidyverse)

gall_network_2012 <-read.csv(text=getURL("https://raw.githubusercontent.com/mabarbour/Genotype_Networks/master/raw_data_management/data/gall_network_data.csv"), header = T)

vLG_network <- gall_network_2012 %>%
  filter(gall.sp == "vLG") %>%
  distinct(Genotype, plant.position, gall.id.nest, g.height, shootEst.all, shootEst.no18, gall_contents, value) %>%
  spread(gall_contents, value, fill = 0)

vLG_network$total <- rowSums(select(vLG_network, Eulo.fem:vLG.pupa))
vLG_network$ectos <- rowSums(select(vLG_network, Eulo.fem, Eulo.mal, Mesopol, Ptero.2, Tory.fem, Tory.mal))

vLG_counts <- vLG_network %>%
  group_by(Genotype, plant.position) %>%
  summarise_at(vars(total), sum) %>%
  rename(vLG_count = total)

vLG_network <- left_join(vLG_network, vLG_counts) %>%
  mutate(Density_per_100_shoots = vLG_count/shootEst.all*100)

## EXPLORATORY PLOTS ----
density.cut.2012 <- quantile(vLG_network$Density_per_100_shoots, probs = c(0,0.33,0.66,1), include.lowest = TRUE)
vLG_network %>%
  mutate(density.cut = cut(Density_per_100_shoots, density.cut.2012)) %>%
  ggplot(., aes(x = g.height, y = vLG.pupa/total, weight = total)) +
  geom_point(alpha = 0.5, aes(size = total)) +
  binomial_smooth(formula = y ~ x) + facet_wrap(~density.cut)
  

## CONTINUED ----

hist(vLG_network$g.height)
hist(sqrt(vLG_network$g.height))

hist(vLG_network$Density_per_100_shoots)
hist(sqrt(vLG_network$Density_per_100_shoots))

vLG_network <- mutate(vLG_network,
                      sc.sqrt.size = (sqrt(g.height) - mean(sqrt(g.height), na.rm = T))/sd(sqrt(g.height), na.rm = T),
                      sc.size = (g.height - mean(g.height, na.rm = T))/sd(g.height, na.rm = T),
                      sc.sqrt.density = (sqrt(Density_per_100_shoots) - mean(sqrt(Density_per_100_shoots)))/sd(sqrt(Density_per_100_shoots)),
                      sc.density = (Density_per_100_shoots - mean(Density_per_100_shoots))/sd(Density_per_100_shoots),
                      Plant_Position = as.factor(plant.position)) 
mean.narm <- function(x) mean(x, na.rm = TRUE)

mean_g.height <- vLG_network %>% 
  group_by(Genotype, Plant_Position) %>%
  summarise_at(vars(g.height), mean.narm) %>%
  rename(mean_g.height.Plant = g.height) 

mean_g.height.Genotype <- mean_g.height %>% 
  group_by(Genotype) %>%
  summarise_at(vars(mean_g.height.Plant), mean.narm) %>%
  rename(mean_g.height.Geno = mean_g.height.Plant)

vLG_network <- left_join(vLG_network, mean_g.height) %>%
  left_join(., mean_g.height.Genotype) %>%
  mutate(c.mean_g.height.Geno = mean_g.height.Geno - mean.narm(mean_g.height.Geno),
         c.mean_g.height.Plant = mean_g.height.Plant - mean_g.height.Geno,
         c.g.height = g.height - mean_g.height.Plant)

test <- glmer(vLG.pupa/total ~ c.mean_g.height.Geno + c.mean_g.height.Plant + c.g.height + (1 | Genotype/Plant_Position), 
                            data = vLG_network, family = binomial(link = "logit"), weights = total)
summary(test)
sem.model.fits(test)
library(lme4)

ectos.2012 <- glmer(ectos/total ~ 1 + (1 | Genotype/Plant_Position), 
                            data = vLG_network, family = binomial(link = "logit"), weights = total)
summary(ectos.2012)
fixef(ectos.2012) # -1.713
# confint.merMod(ectos.2012, method = "boot") # -2.057 to -1.370

gall_survival.2012 <- glmer(vLG.pupa/total ~ sc.sqrt.size*sc.sqrt.density + (1 | Genotype/Plant_Position), 
                            data = vLG_network, family = binomial(link = "logit"), weights = total)
summary(gall_survival.2012)
sem.model.fits(gall_survival.2012)
# gall_survival.95CI <- confint.merMod(gall_survival.2012, method = "boot") # takes a long time

gall_survival.ALT.2012 <- glmer(vLG.pupa/total ~ sc.size*sc.density + (1 | Genotype/Plant_Position), 
                            data = vLG_network, family = binomial(link = "logit"), weights = total)
summary(gall_survival.ALT.2012)
sem.model.fits(gall_survival.ALT.2012)

# pca with these two factors
gall_pca.2012 <- princomp(na.omit(select(vLG_network, sc.sqrt.size, sc.sqrt.density)), cor = TRUE)
summary(gall_pca.2012)
loadings(gall_pca.2012)
biplot(gall_pca.2012, choices = c(1,2))
gall_pca.2012.scores <- data.frame(na.omit(select(vLG_network, Genotype, Plant_Position, gall.id.nest, sc.sqrt.size, sc.sqrt.density)), gall_pca.2012$scores)
#gall_pca.2012.scores$rownames <- rownames(gall_pca.2012.scores)

#vLG_network$rownames <- rownames(vLG_network) 
vLG_network <- left_join(vLG_network, gall_pca.2012.scores)

gall_survival_pca.2012 <- glmer(vLG.pupa/total ~ Comp.1*Comp.2 + (1 | Genotype/Plant_Position), 
                            data = vLG_network, family = binomial(link = "logit"), weights = total)
summary(gall_survival_pca.2012)
library(piecewiseSEM)
sem.model.fits(gall_survival_pca.2012)
# gall_survival_pca.95CI <- confint.merMod(gall_survival_pca.2012, method = "boot") # takes a long time

library(visreg)
visreg2d(gall_survival.2012, xvar = "sc.sqrt.density", yvar = "sc.sqrt.size", scale = "response")

fixef(gall_survival.2012)["sc.sqrt.size"]
fixef(gall_survival.2012)["sc.sqrt.density"]
fixef(gall_survival.2012)["sc.sqrt.size:sc.sqrt.density"]

attr(VarCorr(gall_survival.2012)[["Plant_Position:Genotype"]], "stddev")
attr(VarCorr(gall_survival.2012)[["Genotype"]], "stddev")

# adjusted SD until I adequately covered the 95% confidence intervals of the paramater estimates
quantile(rnorm(n = 1000, fixef(gall_survival.2012)["sc.sqrt.size"], sd = 0.2), probs = c(0.025, 0.975))
quantile(rnorm(n = 1000, fixef(gall_survival.2012)["sc.sqrt.density"], sd = 0.2), probs = c(0.025, 0.975))
quantile(rnorm(n = 1000, fixef(gall_survival.2012)["sc.sqrt.size:sc.sqrt.density"], sd = 0.2), probs = c(0.025, 0.975))

hist(rnorm(n = 1000, mean = 1.713, sd = 0.5)) # this prior seems appropriate for effect of Ectoparasitoid exclusion. I except it to have a strong negative effect on gall survival. 
# I chose a value of -1 because this was close to the effect size of gall size (opposite direction). 
# I chose a relatively wide distribution SD = 1, because I was unsure about the actual magnitude of the effect.
# I imposed an upper bound at zero, becuase it doesn't make sense biologically that excluding parasitoids would increase gall survival.
# I actually tightened up the distribution SD = 0.5, because apparently I can't set a lower bound on a named coefficeint...

# this prior seems appropriate for the unknown fixed-effects. It's possible they could have larger effects than gall size, but I doubt it. 
# I also want to put the mass of this distribution on zero, because I don't know whether these effects will be positive or negative.
hist(rnorm(n = 1000, mean = 0, sd = 1)) 

# for other random effects, I'm centering on zero, but allowing for potentially larger effect sizes.
hist(rnorm(n = 1000, mean = 0, sd = 2)) 

gall_priors <- c(set_prior("normal(1.125, 0.2)", class = "b", coef = "sc.size"),
                 #set_prior("normal(0.087, 0.2)", class = "b", coef = "sc.sqrt.density"),
                 #set_prior("normal(0.332, 0.2)", class = "b", coef = "sc.sqrt.size:sc.sqrt.density"),
                 #set_prior("normal(1.713, 0.5)", class = "b", coef = "Treatment.focusEctoparasitoidexclusion"),
                 set_prior("normal(0, 1)", class = "b"),
                 #set_prior("normal(0.491, 2)", class = "sd", group = "Genotype"),
                 #set_prior("normal(0.489, 2)", class = "sd", group = "Genotype:Plant_Position"),
                 set_prior("normal(0, 2)", class = "sd"))
                 
                 

