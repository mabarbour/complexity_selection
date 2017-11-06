
## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(brms)
library(cowplot)

## LOAD DATASET ----
gall_selection.df <- read_csv("gall_selection_data.csv") %>%
  mutate(Plant_Position = as.factor(Plant_Position),
         Gall_Number = as.factor(Gall_Number),
         Treatment.focus = as.factor(Treatment.focus),
         Treatment = as.factor(Treatment),
         Location = as.factor(Location),
         phenology = as.factor(phenology),
         Genotype = as.factor(Genotype),
         Gall_Letter = as.factor(Gall_Letter)) %>%
  unite(Gall_ID, Gall_Number, Gall_Letter, remove = FALSE) %>%
  
  # subset data for analysis and make useful transformations
  filter(phenology == "early", Location == "tree") %>%
  mutate(gall_survival = ifelse(pupa > 0, 1, 0),
         
         # playing around with effectiveness of transformations in fitting complex models.
         Gall_Height_mm = sqrt(Gall_Height_mm),
         gall_individuals = log(gall_individuals),
         Density_per_100_shoots = sqrt(Density_per_100_shoots)) # interesting that there are occassions with more than one individual within a gall chamber

# get info at different hierarchical levels
mean.NA.rm <- function(x) mean(x, na.rm = TRUE)

size.gall_number <- gall_selection.df %>%
  group_by(Genotype, Plant_Position, Gall_Number) %>%
  summarise_at(vars(Gall_Height_mm, gall_individuals, Density_per_100_shoots), mean.NA.rm) %>%
  ungroup()

size.number_plant <- size.gall_number %>%
  group_by(Genotype, Plant_Position) %>%
  summarise_at(vars(Gall_Height_mm, gall_individuals, Density_per_100_shoots), mean.NA.rm) %>%
  ungroup()

survival_plant <- gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  summarise_at(vars(pupa, total), sum) %>%
  ungroup()

survival_gall <- gall_selection.df %>%
  group_by(Treatment.focus, Genotype, Plant_Position, Gall_Number) %>%
  summarise_at(vars(pupa, total), sum) %>%
  ungroup()

size.number.density_genotype <- size.number_plant %>%
  group_by(Genotype) %>%
  summarise_at(vars(Gall_Height_mm, gall_individuals, Density_per_100_shoots), mean.NA.rm) %>%
  ungroup()

# merge everything back together
merge.Genotype <- left_join(gall_selection.df, 
                            select(size.number.density_genotype, Genotype, size.Geno = Gall_Height_mm, indiv.Geno = gall_individuals, density.Geno = Density_per_100_shoots))
merge.Plant <- left_join(merge.Genotype,
                         select(size.number_plant, Plant_Position, size.Plant = Gall_Height_mm, indiv.Plant = gall_individuals, density.Plant = Density_per_100_shoots))
merge.PolyGall <- left_join(merge.Plant,
                            select(size.gall_number, Gall_Number, size.PolyGall = Gall_Height_mm, indiv.PolyGall = gall_individuals))

# finalize dataset with information at each hierarchical level (gall, polythamolous gall, plant, genotype)
selection.df <- merge.PolyGall %>%
  select(Treatment = Treatment.focus, 
         Genotype, Plant_Position, Gall_Number, Gall_ID, 
         platy, ectos, pupa, gall_survival, total, 
         size.Geno, indiv.Geno, density.Geno,
         size.Plant, indiv.Plant, density.Plant,
         size.PolyGall, indiv.PolyGall,
         size.Gall = Gall_Height_mm) %>%
  mutate(
         # center at Genotype level
         c.size.Geno = size.Geno - mean.NA.rm(size.Geno),
         c.indiv.Geno = indiv.Geno - mean.NA.rm(indiv.Geno),
         c.density.Geno = density.Geno - mean.NA.rm(density.Geno),
         
         # center at Plant level
         c.size.Plant = size.Plant - size.Geno,
         c.indiv.Plant = indiv.Plant - indiv.Geno,
         c.density.Plant = density.Plant - density.Geno,
         
         # center at Polythamolous gall level
         c.size.PolyGall = size.PolyGall - size.Plant,
         c.indiv.PolyGall = indiv.PolyGall - indiv.Plant,
         
         # center at individual level
         c.size.Gall = size.Gall - size.PolyGall)
                       
## PRELIMINARY HIERARCHICAL MODEL

hierarchy_gall.control <-  brm(gall_survival ~ 
                                 (c.size.Geno + c.size.Plant + c.size.PolyGall + c.size.Gall +
                                 c.indiv.Geno + c.indiv.Plant + c.indiv.PolyGall +
                                 c.density.Geno + c.density.Plant)*Treatment +
                                 (1 | Genotype/Plant_Position/Gall_Number/Gall_ID),
                               data = selection.df,
                               family = "bernoulli",
                               algorithm = "sampling",
                               prior = c(set_prior("normal(0,1)", class = "b"),
                                         set_prior("normal(0,2)", class = "sd")),
                               control = list(adapt_delta = 0.99))
summary(hierarchy_gall.control)
plot(marginal_effects(hierarchy_gall.control, effects = c("c.size.PolyGall:Treatment")))
plot(marginal_effects(hierarchy_gall.control, effects = c("c.indiv.PolyGall:Treatment")))
plot(marginal_effects(hierarchy_gall.control, effects = c("c.density.Plant:Treatment")))

hierarchy_gall <-  brm(gall_survival ~ 
                                 c.size.Geno + c.size.Plant + c.size.Gall + c.indiv.Geno + c.indiv.Plant + c.density.Geno + 
                                     c.indiv.PolyGall*c.density.Plant*c.size.PolyGall*Treatment +
                                 (1 | Genotype/Plant_Position/Gall_Number),
                               data = selection.df,
                               family = "bernoulli",
                               algorithm = "sampling",
                               prior = c(set_prior("normal(0,2)", class = "b"),
                                         set_prior("normal(0,2)", class = "sd")),
                               control = list(adapt_delta = 0.99))
summary(hierarchy_gall)
plot(marginal_effects(hierarchy_gall))

plot(marginal_effects(hierarchy_gall, effects = c("c.indiv.PolyGall:c.size.PolyGall"), 
                      conditions = data.frame(Treatment = c("Control", "Ectoparasitoid exclusion")), 
                      surface = TRUE), stype = "raster")

# plant-level
plant_level.df <- left_join(survival_plant, size.number_plant) %>%
  mutate(sc.size = scale(Gall_Height_mm),
         sc.indiv = scale(gall_individuals),
         sc.density = scale(Density_per_100_shoots))

plant_pca <- princomp(select(plant_level.df, sc.size, sc.indiv, sc.density))
summary(plant_pca)
data.frame(loadings(plant_pca)[1])
plant_pca$scale
biplot(plant_pca, choices = c(1,2))
biplot(plant_pca, choices = c(2,3))
biplot(plant_pca, choices = c(1,3))

pca.loadings <- data.frame(Comp.1 = c(-0.622,-0.536,-0.570), Comp.2 = c(0.105, -0.779, 0.618), Comp.3 = c(0.776, -0.325, -0.541))
pca.loadings$label <- c("sc.size","sc.indiv","sc.density")

ggbiplot(data.frame(plant_pca$scores, select(plant_level.df, pupa, total, Treatment.focus)), 
         loadings.data = pca.loadings, 
         loadings = T, size = 0.0001, loadings.label = TRUE, 
         loadings.label.label = pca.loadings$label) + 
  geom_point(aes(size = total, fill = pupa/total), shape = 21) + 
  scale_fill_gradientn(colors = viridis6()) +
  facet_wrap(~Treatment.focus)

autoplot(princomp(select(plant_level.df, sc.size, sc.indiv, sc.density)), 
         loadings = TRUE, label = TRUE,
         loadings.label = TRUE, x = 1, y = 2) +
  geom_point(data = plant_level.df, aes(x = Comp.1, y = Comp.2, fill = pupa/total, size = total), shape = 21) + 
  facet_wrap(~Treatment.focus) + 
  scale_fill_gradientn(colors = viridis6())

autoplot(princomp(select(plant_level.df, sc.size, sc.indiv, sc.density)), 
         loadings = TRUE,
         loadings.label = TRUE, scale = 0.1, x = 2, y = 3) +
  geom_point(data = plant_level.df, aes(x = Comp.2, y = Comp.3, fill = pupa/total, size = total), shape = 21) + 
  facet_wrap(~Treatment.focus) + 
  scale_fill_gradientn(colors = viridis6())

autoplot(princomp(select(plant_level.df, sc.size, sc.indiv, sc.density)), 
         loadings = TRUE,
         loadings.label = TRUE, scale = 0.1, x = 1, y = 3) +
  geom_point(data = plant_level.df, aes(x = Comp.1, y = Comp.2, fill = pupa/total, size = total), shape = 21) + 
  facet_wrap(~Treatment.focus) + 
  scale_fill_gradientn(colors = viridis6())

library(ggfortify)

plant_level.df <- cbind.data.frame(plant_level.df, plant_pca$scores)

car::scatterplotMatrix(select(plant_level.df, sc.size, sc.indiv, sc.density, Treatment.focus))

viridis6 <- function() {
  # colours taken from the viridis package
  c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725")
}

plant_level.df %>%
  ggplot(., aes(x = sc.size, y = sc.indiv, fill = pupa/total, size = total)) +
  geom_point(shape = 21) + scale_fill_gradientn(colors = viridis6()) +
  facet_wrap(~Treatment.focus)

plant_level.df %>%
  ggplot(., aes(x = sc.density, y = sc.indiv, fill = pupa/total, size = total)) +
  geom_point(shape = 21) + scale_fill_gradientn(colors = viridis6()) +
  facet_wrap(~Treatment.focus)

plant_level.df %>%
  ggplot(., aes(x = sc.density, y = sc.size, fill = pupa/total, size = total)) +
  geom_point(shape = 21) + scale_fill_gradientn(colors = viridis6()) +
  facet_wrap(~Treatment.focus)

plant_level.df %>%
  ggplot(., aes(x = Comp.1, y = Comp.2, fill = pupa/total, size = total)) +
  geom_point(shape = 21) + scale_fill_gradientn(colors = viridis6()) +
  facet_wrap(~Treatment.focus)

plant_level.df %>%
  ggplot(., aes(x = Comp.2, y = Comp.3, fill = pupa/total, size = total)) +
  geom_point(shape = 21) + scale_fill_gradientn(colors = viridis6()) +
  facet_wrap(~Treatment.focus)

plant_level.df %>%
  ggplot(., aes(x = Comp.1, y = Comp.3, fill = pupa/total, size = total)) +
  geom_point(shape = 21) + scale_fill_gradientn(colors = viridis6()) +
  facet_wrap(~Treatment.focus)

plant_brm <- brm(pupa | trials(total) ~ sc.size*sc.indiv*sc.density*Treatment.focus + (1|Genotype),
                 data = plant_level.df,
                 family = binomial(link = "logit"),
                 algorithm = "sampling",
                 prior = c(set_prior("normal(0,1)", class = "b"), set_prior("normal(0,2)", class = "sd")))
summary(plant_brm)
summary(plant_level.df$total)

pp_check(plant_brm)
bayes_R2(plant_brm)

plot(marginal_effects(plant_brm))

treat.df <- data.frame(Treatment.focus = c("Control", "Ectoparasitoid exclusion"))
rownames(treat.df) <- treat.df$Treatment.focus

plot(marginal_effects(plant_brm, effects = c("sc.size:sc.indiv"), conditions = treat.df, surface = T), stype = "raster")
plot(marginal_effects(plant_brm, effects = c("sc.indiv:sc.density"), conditions = treat.df, surface = T), stype = "raster")
#plot(marginal_effects(plant_brm, effects = c("Comp.2:Comp.3"), conditions = treat.df, surface = T), stype = "raster")


# gall 
gall_level.df <- left_join(survival_gall, size.gall_number) %>%
  mutate(sc.size = scale(Gall_Height_mm),
         sc.indiv = scale(gall_individuals),
         sc.density = scale(Density_per_100_shoots))

gall_pca <- princomp(select(gall_level.df, sc.size, sc.indiv, sc.density))
summary(gall_pca)
loadings(gall_pca)
biplot(gall_pca, choices = c(1,2))
biplot(gall_pca, choices = c(2,3))
biplot(gall_pca, choices = c(1,3))

gall_level.df <- cbind.data.frame(gall_level.df, gall_pca$scores)

car::scatterplotMatrix(gall_pca$scores)

gall_level.df %>%
  ggplot(., aes(x = sc.size, y = sc.indiv, fill = pupa/total, size = total)) +
  geom_point(shape = 21) + scale_fill_gradientn(colors = viridis6()) +
  facet_wrap(~Treatment.focus)

plant_level.df %>%
  ggplot(., aes(x = sc.density, y = sc.indiv, fill = pupa/total, size = total)) +
  geom_point(shape = 21) + scale_fill_gradientn(colors = viridis6()) +
  facet_wrap(~Treatment.focus)

gall_brm <- brm(pupa | trials(total) ~ Comp.1*Comp.2*Comp.3*Treatment.focus + (1|Genotype/Plant_Position),
                 data = gall_level.df,
                 family = binomial(link = "logit"),
                 algorithm = "sampling",
                control = list(adapt_delta = 0.99),
                 prior = c(set_prior("normal(0,1)", class = "b"), set_prior("normal(0,2)", class = "sd")))
summary(gall_brm)

pp_check(gall_brm)
bayes_R2(gall_brm)

plot(marginal_effects(gall_brm))

treat.df <- data.frame(Treatment.focus = c("Control", "Ectoparasitoid exclusion"))
rownames(treat.df) <- treat.df$Treatment.focus

plot(marginal_effects(gall_brm, effects = c("Comp.1:Comp.2"), conditions = treat.df, surface = T), stype = "raster")
plot(marginal_effects(gall_brm, effects = c("Comp.1:Comp.3"), conditions = treat.df, surface = T), stype = "raster")
plot(marginal_effects(gall_brm, effects = c("Comp.2:Comp.3"), conditions = treat.df, surface = T), stype = "raster")
