## Notes ----
# Restrict analysis to ground, tree, bag, mark, nomark, nobag galls only.
# Focus: Want to determine whether the proportion of galls parasitized by Platygaster varies among willow genotypes and bagging treatment. If Platygaster parasitism is higher in bagged galls, then this is evidence of intraguild predation. This may also depend on willow genotype.
# I don't expect there to be a difference between mark/nomark (paired on same willows), as this would indicate that marking the galls may have affected the foraging behavior of parasitoids. However, an effect could also mean that differences in phenology affect parasitism.
# I do expect that Platygaster parasitism will be higher in bag vs. no-bag galls (paired on same willows) as well as higher in bag vs. mark/nomark galls. This will indicate evidence of intraguild predation.
# I expect there to be a difference between no-bag and nomark galls, because no-bag galls were on trees with highly reduced gall densities. If previous work is correct, there may be higher ectoparasitism rates on these no-bag galls.
# I don't expect a different in Platygaster parasitism between ground vs. tree galls. Note that I only bag and nomark galls have both ground and tree galls, therefore I should restrict my comparison to these ones.

## source functions and other general details ----

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # color-blind palette with grey

## load libraries ----
library(tidyverse)
#library(lme4)
library(cowplot) # for clean multi-panel plots

## read in data frames ----

# parasitism data
vLG.df <- read.csv("MBGalls/VLGCollect2013.csv") %>% tbl_df

# plant information
plant.info <- read.csv("Willow Garden Positions.csv") %>%
  tbl_df %>%
  select(Row = Row.., Plant_Position = Plant.Position, Genotype, Gender) %>%
  mutate(Plant_Position = as.character(Plant_Position))

# upload gall density survey
density <- read.csv("Density_Surveys_Gall_Parasitoid_Bagging_Experiment.csv", skip=1) %>%
  tbl_df %>%
  mutate(Plant_Position = as.character(Plant.Position),
         Total_Shoot_Estimate = Shoots_on_last_sampled_branch*5,
         Density_per_100_shoots = Iteomyia_count/Total_Shoot_Estimate*100,
         Treatment.check = Treatment) %>%
  select(Treatment.check, 
         Plant_Position, Iteomyia_count, Density_per_100_shoots, Total_Shoot_Estimate) %>%
  filter(!is.na(Iteomyia_count))
rowSums(with(density, table(Treatment.check,Genotype)))

density.posit <- names(table(density$Plant_Position)) # extract plant positions
length(density.posit)

## join data frames together ----
vLG.df <- left_join(vLG.df, plant.info, by = "Plant_Position") %>%
  left_join(density, by = "Plant_Position") # %>% # add density data frame

vLG.df$Treatment.focus <- plyr::revalue(vLG.df$Treatment,
                                       c("mark" = "Control", 
                                         "nomark" = "Control", 
                                         "bag" = "Ectoparasitoid exclusion", 
                                         "nobag" = "Reduced density"))
         
vLG.df$phenology <- plyr::revalue(vLG.df$Treatment,
                                  c("mark" = "early",
                                    "nomark" = "unk",
                                    "bag" = "early",
                                    "nobag" = "unk"))
vLG.df.posit <- names(table(vLG.df$Plant_Position)) # extract plant positions in the dataset
length(vLG.df.posit)

# identify potentially erroneous plant positions, that weren't part of original Iteomyia density survey
true.posits <- which(vLG.df.posit %in% density.posit == TRUE)
true.names <- vLG.df.posit[true.posits]

# limit dataset to plant positions we can verify
vLG.df <- filter(vLG.df, Plant_Position %in% true.names)

## visualize data structure ----
table(vLG.df$Plant_Position, vLG.df$Treatment) # see which treatment? I can salvage
table(vLG.df$Treatment, vLG.df$Location)

# view all gall contents
table(vLG.df$Contents, exclude = NULL) # why 4 NA?

## tidy data at individual level ----
IGP.df.ind <- vLG.df %>%
  filter(Location %in% c("ground", "tree"), 
         Treatment %in% c("bag","mark","nobag","nomark"),
         !is.na(Contents)) %>% 
  # clean up content names for future use of select() 
  mutate(Contents = gsub(" ", ".", Contents, fixed = TRUE)) %>% 
  mutate(Contents = gsub("-", "_", Contents, fixed = TRUE)) %>% 
  mutate(Contents = gsub("?", "_", Contents, fixed = TRUE)) %>%
  mutate(Plant_Position = as.factor(Plant_Position)) %>%
  group_by(Genotype, Plant_Position, Location, Treatment, Treatment.check, 
           Treatment.focus, phenology,
           Gall_Number, Gall_Letter, Contents) %>%
  summarise(count = n(),
            Gall_Height_mm = mean(Point_Height_mm),
            Iteomyia_count = mean(Iteomyia_count),
            Density_per_100_shoots = mean(Density_per_100_shoots)) %>% # note that each Gall_Number:Gall_Letter combination is unique, so taking the mean here just prevents me from double counting these galls.
  spread(Contents, count, fill = 0)

IGP.df.ind$total <- IGP.df.ind %>%
  ungroup() %>%
  select(brown.mush:vlg.larva) %>% # all possible contents
  rowSums()

IGP.df.ind$ectos <- IGP.df.ind %>%
  ungroup() %>%
  select(eulo.female:eulo.us, # exclude eury_. because they don't attack vLG
         exit__none_ecto, # remove mymarid because it is an egg parasitoid
         exit__none_no.ecto, # remove Exit_none because it is possible this was caused by an egg parasitoid
         exit_none_ecto:hairy.larva,
         meso.female:meso.us, # remove mymarid because it is an egg parasitoid
         paraunk,
         ptoid.partial,
         tory.female:tory.us) %>%
  rowSums()

IGP.df.ind$platy.ectos <- IGP.df.ind$platy + IGP.df.ind$ectos

gall.individuals.df <- IGP.df.ind %>%
  ungroup() %>%
  group_by(Genotype, Treatment.focus, Treatment.check, Treatment, 
           Location, phenology, Plant_Position, Gall_Number) %>%
  summarise(gall_individuals = n()) %>%
  ungroup() %>%
  select(Gall_Number, gall_individuals)

IGP.df.ind.galln <- left_join(IGP.df.ind, gall.individuals.df, by = "Gall_Number")

IGP.df.ind.galln %>%
  ungroup() %>%
  select(brown.mush:ectos) %>%
  colSums() %>%
  sort() # I think some of the totals and ectos are so high because it is aggregating the unknown Gall_Numbers as the same...This doesn't affect the datasets below though because I'm removing these unknown galls.

## COMPLETE DATASET
ind.tree.data <- IGP.df.ind.galln %>%
  ungroup() %>%
  filter(Gall_Number != "unk") %>%
  select(Treatment.focus, Treatment, Location, phenology, Genotype, Plant_Position, Gall_Number:Density_per_100_shoots, gall_individuals,
         platy, pupa, ectos, platy.ectos, none.none, total) %>%
  droplevels.data.frame() 
write_csv(ind.tree.data, "gall_selection_data.csv")

table(ind.tree.data$Treatment.focus, ind.tree.data$Treatment)
table(ind.tree.data$Treatment, ind.tree.data$Location) # note that there aren't any situations with "marked" ground galls, which makes this comparison difficult.
table(ind.tree.data$Treatment, ind.tree.data$phenology)

######################## EVERYTHING BELOW IS OLD, BUT MAYBE USEFUL ##########################################

## create focused datasets for analysis
ind.tree.IGP <- IGP.df.ind.galln %>%
  ungroup() %>%
  filter(Treatment.focus %in% c("Control", "Ectoparasitoid exclusion"),
         Location == "tree", 
         phenology == "early",
         Gall_Number != "unk") %>% # couldn't link parasitoid to a particular gall
  select(Treatment.focus, Treatment, Location, phenology, Genotype, Plant_Position, Gall_Number:Density_per_100_shoots, gall_individuals,
         platy, pupa, ectos, platy.ectos, none.none, total)

ind.tree.red <- IGP.df.ind.galln %>%
  ungroup() %>%
  filter(Treatment.focus %in% c("Control", "Reduced density"),
         Location == "tree",
         phenology == "unk", # this wasn't originally here, but I think it is necessary to control for potential phenology mismatches...
         Gall_Number != "unk") %>% # couldn't link parasitoid to a particular gall
  select(Treatment.focus, Treatment, phenology, Genotype, Plant_Position, Gall_Number:Density_per_100_shoots, gall_individuals,
         platy, pupa, ectos, platy.ectos, none.none, total)

ind.tree.cont <- IGP.df.ind.galln %>%
  ungroup() %>%
  filter(Treatment.focus %in% c("Control"),
         Location == "tree",
         Gall_Number != "unk") %>% # couldn't link parasitoid to a particular gall
  select(Treatment.focus, Treatment, phenology, Genotype, Plant_Position, Gall_Number:Density_per_100_shoots, gall_individuals,
         platy, pupa, ectos, platy.ectos, none.none, total)

## tidy data at plant level for plotting ----
gall.tree.IGP <- ind.tree.IGP %>%
  group_by(Treatment.focus, Genotype, Plant_Position, Gall_Number) %>%
  # summarise at gall level first
  summarise(Gall_Height_mm = mean(Gall_Height_mm, na.rm = TRUE),
            Iteomyia_count = mean(Iteomyia_count),
            Density_per_100_shoots = mean(Density_per_100_shoots),
            gall_individuals = mean(gall_individuals),
            platy = sum(platy),
            pupa = sum(pupa),
            ectos = sum(ectos),
            none.none = sum(none.none),
            total = sum(total)) %>%
  mutate(prop.platy = platy/total,
         prop.pupa = pupa/total,
         prop.ectos = ectos/total,
         prop.none.none = none.none/total) 

plant.tree.IGP <- gall.tree.IGP %>%
  group_by(Treatment.focus, Genotype, Plant_Position) %>%
  # summarise at plant level now
  summarise(Gall_Height_mm = mean(Gall_Height_mm, na.rm = TRUE),
            Iteomyia_count = mean(Iteomyia_count),
            Density_per_100_shoots = mean(Density_per_100_shoots),
            gall_individuals = mean(gall_individuals),
            prop.platy = weighted.mean(prop.platy, total, na.rm = TRUE),
            prop.pupa = weighted.mean(prop.pupa, total, na.rm = TRUE),
            prop.ectos = weighted.mean(prop.ectos, total, na.rm = TRUE),
            prop.none.none = weighted.mean(prop.none.none, total, na.rm = TRUE),
            total = sum(total)) %>%
  ungroup()

genotype.tree.control.IGP <- plant.tree.IGP %>%
  filter(Treatment.focus == "Control") %>%
  group_by(Genotype) %>%
  summarise(Gall_Height_mm = mean(Gall_Height_mm, na.rm = TRUE),
            Iteomyia_count = mean(Iteomyia_count),
            Density_per_100_shoots = mean(Density_per_100_shoots))
vLG.2011.2012 <- read.csv('leaf gall densities and sizes 2011 and 2012.csv')
vLG.2011.12.13 <- left_join(vLG.2011.2012, genotype.tree.control.IGP, by = "Genotype")

ggplot(filter(vLG.2011.12.13),#, Genotype != "X"), 
       aes(y = Gall_Height_mm,
           x = vLG.height.mean)) +
  geom_text(aes(label = Genotype)) +
  geom_smooth(method = "lm")
with(vLG.2011.12.13, cor.test(Gall_Height_mm, vLG.height.mean), method = c("spearman","pearson")) # highly sig. correlation

ggplot(filter(vLG.2011.12.13),#, Genotype != "X"), 
       aes(y = sqrt(vLG.2011.density),
           x = sqrt(vLG.2012.density))) +
  geom_text(aes(label = Genotype)) +
  geom_smooth(method = "lm")
with(filter(vLG.2011.12.13, Genotype != "T", Genotype != "X"), 
     cor.test(vLG.2012.density, Iteomyia_count), 
     method = c("spearman"))

library(psych)
corr.test(select(plant.tree.IGP, Gall_Height_mm, Iteomyia_count, gall_individuals)) # weak, non-significant correlations among predictor variables


## Hypotheses & Predictions ----

# H1: More complex food webs will reduce overall gall survival due to removal of functionally distinct parasitoids (at least in short-term).
# P1 - gall survival is lower when exposed to more complex food webs (i.e. "control" treatment).
# P2 - ectoparasitism rates are much higher in more complex food web ("control").
# P3 - apparent parasitism from egg parasitoid is reduced in more complex food web (i.e. evidence of intraguild predation due to "control")

# H2: Gall fitness differences among willow genotypes will be reduced in more complex communities due to complimentarity in trophic interactions with parasitoids.
# P4 - Variation in gall survival is lower in more complex communities.
# P5 - Endoparasitoids and ectoparasitoids vary in their response to variation in gall attributes (size, individuals, and density)

## Plots ----

## In support of P1, gall survival tends to be lower when exposed to more complex food webs.
## In support of P4, variation in gall survival appears to be lower in more complex communities.
# order for plotting aesthetics
plant.tree.IGP.order <- plant.tree.IGP
plant.tree.IGP.order$Treatment.focus <- factor(as.character(plant.tree.IGP.order$Treatment.focus),
                                               levels = c("Control","Ectoparasitoid exclusion"), ordered = TRUE)

# create dataframe of weighted means for proportion of galls surviving on each tree.
wt.means.pupa <- plant.tree.IGP.order %>%
  group_by(Treatment.focus, Genotype) %>%
  summarise(wt.means = weighted.mean(prop.pupa, total))
# 4-fold variation in simple food web
# 3.2-fold variation in complex food web
# variation in larval survival among willow genotypes in complex food web is ~27% less variable compared to simple web

# plot
survival_treat_geno_plot <- ggplot(plant.tree.IGP.order, 
                                   aes(x = Treatment.focus, y = prop.pupa,
                                       fill = Genotype)) +
  geom_point(aes(size = total, color = Genotype), shape = 2, 
             position = position_jitterdodge(dodge.width = 0.25)) +
  geom_line(data = wt.means.pupa, aes(x = Treatment.focus, y = wt.means, 
                                      group = Genotype, linetype = Genotype),
            position = position_jitterdodge(dodge.width = 0.25)) +
  geom_point(data = wt.means.pupa, 
             aes(x = Treatment.focus, y = wt.means, fill = Genotype), 
             shape = 21, size = 15,
             position = position_jitterdodge(dodge.width = 0.25)) +
  scale_fill_manual(values = cbPalette, guide = FALSE) +
  scale_linetype_manual(values = c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "12345678"), guide = FALSE) +
  scale_color_manual(values = cbPalette, guide = FALSE) +
  ylab("Probability of gall survival") +
  xlab("") + 
  theme(legend.position = "none",
        axis.text = element_text(size = 25))
#ggsave("survival_treatment_genotype.pdf", survival_treat_geno_plot, height = 8.5, width = 11, units = "in")

## In support of P2, ectoparasitism rates are much higher in more complex food web ("control").
# create dataframe of weighted means for proportion of ectoparasitoids surviving on each tree.
wt.means.ectos <- plant.tree.IGP %>%
  group_by(Treatment.focus, Genotype) %>%
  summarise(wt.means = weighted.mean(prop.ectos, total))

# plot
ggplot(plant.tree.IGP, aes(x = Treatment.focus, y = prop.ectos, fill = Genotype)) +
  geom_point(aes(size = total, color = Genotype), shape = 2, 
             position = position_jitterdodge(dodge.width = 0.1)) +
  geom_line(data = wt.means.ectos, aes(x = Treatment.focus, y = wt.means, 
                                       group = Genotype, linetype = Genotype),
            position = position_jitterdodge(dodge.width = 0.1)) +
  geom_point(data = wt.means.ectos, 
             aes(x = Treatment.focus, y = wt.means, fill = Genotype), 
             shape = 21, size = 8,
             position = position_jitterdodge(dodge.width = 0.1)) +
  ylab("Prop. of galls attacked by ectoparasitoids") +
  xlab("Treatment") +
  theme_bw()

## In support of P3, apparent parasitism from egg parasitoid is reduced in more complex food web (i.e. evidence of intraguild predation due to "control")
# create dataframe of weighted means for proportion of ectoparasitoids surviving on each tree.
wt.means.platy <- plant.tree.IGP %>%
  group_by(Treatment.focus, Genotype) %>%
  summarise(wt.means = weighted.mean(prop.platy, total))

# plot
ggplot(plant.tree.IGP, aes(x = Treatment.focus, y = prop.platy, fill = Genotype)) +
  geom_point(aes(size = total, color = Genotype), shape = 2, 
             position = position_jitterdodge(dodge.width = 0.1)) +
  geom_line(data = wt.means.platy, aes(x = Treatment.focus, y = wt.means, 
                                       group = Genotype, linetype = Genotype),
            position = position_jitterdodge(dodge.width = 0.1)) +
  geom_point(data = wt.means.platy, 
             aes(x = Treatment.focus, y = wt.means, fill = Genotype), 
             shape = 21, size = 8,
             position = position_jitterdodge(dodge.width = 0.1)) +
  ylab("Prop. of galls attacked by endoparasitoids") +
  xlab("Treatment") +
  theme_bw()

## Analyses ----

## gall attributes. Expected no treatment effect, only Genotype effects. ----

# strong effect of Genotype on number of individuals per gall.
individ.per.gall <- glmer(gall_individuals ~ Treatment.focus + Genotype + (1|Plant_Position), data = gall.tree.IGP, family = "poisson")
summary(individ.per.gall)
overdisp_fun_GLMM(individ.per.gall)
drop1(individ.per.gall, test = "Chi")


# refuge? strong effect of genotype
ind.height <- lmer(gall_individuals*Gall_Height_mm ~ Treatment.focus+Genotype + (1|Plant_Position), data = gall.tree.IGP)
summary(ind.height)
drop1(ind.height, test = "Chisq")

plot(gall_individuals*Gall_Height_mm ~ Gall_Height_mm, gall.tree.IGP)

# other? strong effect of genotype
plot((Gall_Height_mm*Iteomyia_count) ~ Genotype, plant.tree.IGP)
height.density <- lm(sqrt(Gall_Height_mm*Iteomyia_count) ~ Treatment.focus*Genotype, data = plant.tree.IGP)
summary(height.density)
anova(height.density)
#plot(height.density)

# other2? strong effect of genotype
ind.density <- lm(sqrt(gall_individuals*Iteomyia_count) ~ Treatment.focus*Genotype, data = plant.tree.IGP)
summary(ind.density)
anova(ind.density)
#plot(ind.density)

# strong effect of genotype on gall height as well as a small, but significant treatment effect (control galls were smaller, suggesting ectoparasitoid may modify gall size)

gall.height <- lmer(Gall_Height_mm ~ Treatment.focus + Genotype + (1|Plant_Position/Gall_Number), data = ind.tree.IGP)
summary(gall.height)
drop1(gall.height, test = "Chisq") 
visreg(gall.height, xvar = "Treatment.focus", by = "Genotype")

# strong effect of genotype on gall density


gall.density <- glmer(Iteomyia_count ~ Treatment.focus+Genotype + (1|Plant_Position),
                      data = plant.tree.IGP, 
                      family = "poisson",
                      control = glmerControl(optimizer="bobyqa", 
                                             optCtrl = list(maxfun = 100000)))
summary(gall.density)
overdisp_fun_GLMM(gall.density)
drop1(gall.density, test = "Chisq")

## gall attribute plots combined ----
dimater.plot <- ggplot(plant.tree.IGP, aes(x = Genotype, y = Gall_Height_mm, fill = Genotype)) +
  geom_boxplot() +
  scale_fill_manual(values = cbPalette, guide = FALSE) +
  ylab("Gall diameter (mm)") +
  xlab("Willow genotype")  + theme(axis.text = element_text(size = 20))

gall.ind_plot <- ggplot(plant.tree.IGP, aes(x = Genotype, y = gall_individuals, fill = Genotype)) +
  geom_boxplot() +
  scale_fill_manual(values = cbPalette, guide = FALSE) +
  ylab("No. of larva per gall") +
  xlab("Willow genotype")  + theme(axis.text = element_text(size = 20))

density.plot <- ggplot(plant.tree.IGP, aes(x = Genotype, y = Iteomyia_count/5, fill = Genotype)) +
  geom_boxplot() +
  scale_fill_manual(values = cbPalette, guide = FALSE) +
  xlab("Willow genotype") +
  ylab("No. of larva per branch") + theme(axis.text = element_text(size = 20))

gall_attr_plot <- plot_grid(dimater.plot, gall.ind_plot, density.plot, ncol = 3, scale = 0.8)
ggsave("gall_attributes.pdf", gall_attr_plot, width = 11, height = 7, units = "in")

## gall survival ----
pupa <- glmer(pupa > 0 ~ Treatment.focus*Genotype + (1|Plant_Position/Gall_Number),
              data = ind.tree.IGP, 
              family = "binomial",
              control = glmerControl(optimizer="bobyqa", 
                                     optCtrl = list(maxfun = 100000)))
summary(pupa)
overdisp_fun_GLMM(pupa)
drop1(pupa, test = "Chisq") # significant GxE
visreg(pupa, by = "Genotype", xvar = "Treatment.focus", scale = "response")

ind.tree.IGP$Treatment.focus <- relevel(ind.tree.IGP$Treatment.focus, ref = "Control")
test <- glmer(pupa > 0 ~ (Treatment.focus + scale(Gall_Height_mm) + scale(Iteomyia_count) + scale(gall_individuals))^2 + (1|Genotype/Plant_Position/Gall_Number),
              data = ind.tree.IGP, 
              family = "binomial",
              control = glmerControl(optimizer="bobyqa", 
                                     optCtrl = list(maxfun = 100000)))
summary(test)

# survival in simple community ----
simple.dredge.df <- ind.tree.IGP %>%
  filter(Treatment.focus == "Ectoparasitoid exclusion")

# models often had large eigenvalue ratios so it was advisable to scale the variables.
simple.dredge.df$Gall_Height_mm <- scale(simple.dredge.df$Gall_Height_mm)
simple.dredge.df$Iteomyia_count <- scale(simple.dredge.df$Iteomyia_count)
simple.dredge.df$gall_individuals <- scale(simple.dredge.df$gall_individuals)

pupa.simp <- glmer(pupa > 0 ~ 
                     (Gall_Height_mm + Iteomyia_count + gall_individuals)^2 +
                     (1|Plant_Position/Gall_Number),
                   data = simple.dredge.df, 
                   family = "binomial", 
                   na.action = na.fail, 
                   control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(pupa.simp)

pupa.simp.dredge <- dredge(pupa.simp, trace = FALSE, rank = "AICc")
models.pupa.simp <- get.models(pupa.simp.dredge, subset = cumsum(weight) <= 0.95)
avg.pupa.simp <- model.avg(models.pupa.simp)
summary(avg.pupa.simp)

## simple model prediction data frame ----
min.height.scale.simp <- min(simple.dredge.df$Gall_Height_mm)
max.height.scale.simp <- max(simple.dredge.df$Gall_Height_mm)
min.count.scale.simp <- min(simple.dredge.df$Iteomyia_count)
max.count.scale.simp <- max(simple.dredge.df$Iteomyia_count)
min.ind.scale.simp <- min(simple.dredge.df$gall_individuals)
max.ind.scale.simp <- max(simple.dredge.df$gall_individuals)

pred.frame.simp <- data.frame(Iteomyia_count = 0,
                              Gall_Height_mm = seq(min.height.scale.simp,
                                                   max.height.scale.simp,
                                                   0.1),
                              gall_individuals = 0)
pred.frame.simp.count <- data.frame(Iteomyia_count = seq(min.count.scale.simp,
                                                         max.count.scale.simp,
                                                         0.1),
                                    Gall_Height_mm = 0,
                                    gall_individuals = 0)
pred.frame.simp.ind <- data.frame(Iteomyia_count = 0,
                                  Gall_Height_mm = 0,
                                  gall_individuals = seq(min.ind.scale.simp,
                                                         max.ind.scale.simp,
                                                         0.1))

pred.frame.simp$prop.pupa <- predict(avg.pupa.simp, newdata = pred.frame.simp, re.form = NA, type = "response")
pred.frame.simp$Gall_Height_mm_backscale <- pred.frame.simp$Gall_Height_mm*attr(simple.dredge.df$Gall_Height_mm, 'scaled:scale') + attr(simple.dredge.df$Gall_Height_mm, 'scaled:center')

pred.frame.simp.count$prop.pupa <- predict(avg.pupa.simp, newdata = pred.frame.simp.count, re.form = NA, type = "response")
pred.frame.simp.count$Iteomyia_count_backscale <- pred.frame.simp.count$Iteomyia_count*attr(simple.dredge.df$Iteomyia_count, 'scaled:scale') + attr(simple.dredge.df$Iteomyia_count, 'scaled:center')

pred.frame.simp.ind$prop.pupa <- predict(avg.pupa.simp, newdata = pred.frame.simp.ind, re.form = NA, type = "response")
pred.frame.simp.ind$gall_individuals_backscale <- pred.frame.simp.ind$gall_individuals*attr(simple.dredge.df$gall_individuals, 'scaled:scale') + attr(simple.dredge.df$gall_individuals, 'scaled:center')

## simple model plot predictions ----
# gall height
summary(filter(plant.tree.IGP, Treatment.focus == "Ectoparasitoid exclusion")$Gall_Height_mm)
simp.height.plot <- ggplot(filter(plant.tree.IGP, Treatment.focus == "Ectoparasitoid exclusion"), 
                           aes(x = Gall_Height_mm, y = prop.pupa)) +
  geom_point(aes(size = total), shape = 1) +
  geom_line(data = pred.frame.simp, aes(x = Gall_Height_mm_backscale, y = prop.pupa), size = 3) +
  scale_x_continuous(limits = c(4.9, 12), breaks = c(5,7,9,11)) +
  #scale_shape_manual(values = c(1,2)) +
  ylab("Probability of larva survival") +
  xlab("Gall diameter (mm)") +
  theme(legend.position = "none",
        axis.text = element_text(size = 20))

# Iteomyia count
summary(filter(plant.tree.IGP, Treatment.focus == "Ectoparasitoid exclusion")$Iteomyia_count)
85/5
simp.count.plot <- ggplot(filter(plant.tree.IGP, Treatment.focus == "Ectoparasitoid exclusion"), aes(x = Iteomyia_count/5, 
                                                                                                     y = prop.pupa)) +
  geom_point(aes(size = total), shape = 1) +
  geom_line(data = pred.frame.simp.count, aes(x = Iteomyia_count_backscale/5, y = prop.pupa), size = 3) +
  scale_x_continuous(limits = c(0,17), breaks = c(0,5,10,15)) +
  #scale_shape_manual(values = c(1,2)) +
  ylab("") +
  xlab("No. of larva per branch") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.position = "none")

# gall individuals
summary(filter(plant.tree.IGP, Treatment.focus == "Ectoparasitoid exclusion")$gall_individuals)
simp.ind.plot <- ggplot(filter(plant.tree.IGP, Treatment.focus == "Ectoparasitoid exclusion"), aes(x = gall_individuals, 
                                                                                                   y = prop.pupa)) +
  geom_point(aes(size = total), shape = 1) +
  geom_line(data = pred.frame.simp.ind, aes(x = gall_individuals_backscale, y = prop.pupa), size = 3) +
  scale_x_continuous(limits = c(1,6), breaks = seq(1,6,1)) +
  ylab("") +
  xlab("No. of larva per gall") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.position = "none")
#theme_bw() + theme(axis.text.y = element_blank())

simp.plots <- plot_grid(simp.height.plot, simp.ind.plot, simp.count.plot, 
                        ncol = 3, align = "vh")
ggsave("gall_survival_simple.pdf", simp.plots, height = 8.5, width = 11, units = "in")

# survival in complex community ----
complex.dredge.df <- ind.tree.IGP %>%
  filter(Treatment.focus == "Control")

# models often had large eigenvalue ratios so it was advisable to scale the variables.
complex.dredge.df$Gall_Height_mm <- scale(complex.dredge.df$Gall_Height_mm)
complex.dredge.df$Iteomyia_count <- scale(complex.dredge.df$Iteomyia_count)
complex.dredge.df$gall_individuals <- scale(complex.dredge.df$gall_individuals)

pupa.comp <- glmer(pupa > 0 ~ 
                     (Gall_Height_mm + Iteomyia_count + gall_individuals)^2 +
                     (1|Plant_Position/Gall_Number),
                   data = complex.dredge.df, 
                   family = "binomial", 
                   na.action = na.fail, 
                   control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(pupa.comp)

pupa.comp.dredge <- dredge(pupa.comp, trace = FALSE, rank = "AICc")
models.pupa.comp <- get.models(pupa.comp.dredge, subset = cumsum(weight) <= 0.95)
avg.pupa.comp <- model.avg(models.pupa.comp )
summary(avg.pupa.comp)

predict(avg.pupa.comp, type = "response")

## data frame for predictions for complex model ----
min.height.scale.comp <- min(complex.dredge.df$Gall_Height_mm)
max.height.scale.comp <- max(complex.dredge.df$Gall_Height_mm)
min.count.scale.comp <- min(complex.dredge.df$Iteomyia_count)
max.count.scale.comp <- max(complex.dredge.df$Iteomyia_count)
min.ind.scale.comp <- min(complex.dredge.df$gall_individuals)
max.ind.scale.comp <- max(complex.dredge.df$gall_individuals)

pred.frame.comp <- data.frame(Iteomyia_count = 0,
                              Gall_Height_mm = seq(min.height.scale.comp,
                                                   max.height.scale.comp,
                                                   0.1),
                              gall_individuals = 0)
pred.frame.comp.count <- data.frame(Iteomyia_count = seq(min.count.scale.comp,
                                                         max.count.scale.comp,
                                                         0.1),
                                    Gall_Height_mm = 0,
                                    gall_individuals = 0)
pred.frame.comp.ind <- data.frame(Iteomyia_count = 0,
                                  Gall_Height_mm = 0,
                                  gall_individuals = seq(min.ind.scale.comp,
                                                         max.ind.scale.comp,
                                                         0.1))

pred.frame.comp$prop.pupa <- predict(avg.pupa.comp, newdata = pred.frame.comp, re.form = NA, type = "response")
pred.frame.comp$Gall_Height_mm_backscale <- pred.frame.comp$Gall_Height_mm*attr(complex.dredge.df$Gall_Height_mm, 'scaled:scale') + attr(complex.dredge.df$Gall_Height_mm, 'scaled:center')

pred.frame.comp.count$prop.pupa <- predict(avg.pupa.comp, newdata = pred.frame.comp.count, re.form = NA, type = "response")
pred.frame.comp.count$Iteomyia_count_backscale <- pred.frame.comp.count$Iteomyia_count*attr(complex.dredge.df$Iteomyia_count, 'scaled:scale') + attr(complex.dredge.df$Iteomyia_count, 'scaled:center')

pred.frame.comp.ind$prop.pupa <- predict(avg.pupa.comp, newdata = pred.frame.comp.ind, re.form = NA, type = "response")
pred.frame.comp.ind$gall_individuals_backscale <- pred.frame.comp.ind$gall_individuals*attr(complex.dredge.df$gall_individuals, 'scaled:scale') + attr(complex.dredge.df$gall_individuals, 'scaled:center')

## complex model plot predictions ----

# gall height
summary(filter(plant.tree.IGP, Treatment.focus == "Control")$Gall_Height_mm)
comp.height.plot <- ggplot(filter(plant.tree.IGP, Treatment.focus == "Control"), 
                           aes(x = Gall_Height_mm, y = prop.pupa)) +
  geom_point(aes(size = total), shape = 1) +
  geom_line(data = pred.frame.comp, aes(x = Gall_Height_mm_backscale, y = prop.pupa), size = 3) +
  scale_x_continuous(limits = c(4.7, 11), breaks = c(5,7,9,11)) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25, 0.5, 0.75, 1.0)) +
  ylab("Probability of larva survival") +
  xlab("Gall diameter (mm)") +
  theme(legend.position = "none",
        axis.text = element_text(size = 20))

# Iteomyia count
summary(filter(plant.tree.IGP, Treatment.focus == "Control")$Iteomyia_count)
90/5
comp.count.plot <- ggplot(filter(plant.tree.IGP, Treatment.focus == "Control"), 
                          aes(x = Iteomyia_count/5, 
                              y = prop.pupa)) +
  geom_point(aes(size = total), shape = 1) +
  geom_line(data = pred.frame.comp.count, aes(x = Iteomyia_count_backscale/5, y = prop.pupa), size = 3, color = "gray75") +
  scale_x_continuous(limits = c(0,18), breaks = c(0,5,10,15)) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25, 0.5, 0.75, 1.0)) +
  #scale_shape_manual(values = c(1,2)) +
  ylab("") +
  xlab("No. of larva per branch") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.position = "none")

# gall individuals
summary(filter(plant.tree.IGP, Treatment.focus == "Control")$gall_individuals)
comp.ind.plot <- ggplot(filter(plant.tree.IGP, Treatment.focus == "Control"), aes(x = gall_individuals, 
                                                                                  y = prop.pupa)) +
  geom_point(aes(size = total), shape = 1) +
  geom_line(data = pred.frame.comp.ind, aes(x = gall_individuals_backscale, y = prop.pupa), size = 3, color = "gray75") +
  scale_x_continuous(limits = c(1,6), breaks = seq(1,6,1)) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25, 0.5, 0.75, 1.0)) +
  #scale_shape_manual(values = c(1,2)) +
  ylab("") +
  xlab("No. of larva per gall") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.position = "none")

comp.plots <- plot_grid(comp.height.plot, comp.ind.plot, comp.count.plot, 
                        ncol = 3, align = "vh")
ggsave("gall_survival_complex.pdf", comp.plots, height = 8.5, width = 11, units = "in")


## IGP strength ----
# Does strength of IGP vary among willow genotypes?
platy <- glmer(platy > 0 ~ Treatment.focus + Genotype + (1|Plant_Position/Gall_Number),
               data = ind.tree.IGP, family = "binomial",
               control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(platy)
overdisp_fun_GLMM(platy)
drop1(platy, test = "Chisq") # likelihood ratio tests
#visreg(igp)
#confint(profile(igp))

# focusing solely on "simple" treatment to identify mechanisms determining gall survival in this simple community

platy.mech <- glmer(platy > 0 ~ (Gall_Height_mm + Iteomyia_count + gall_individuals)^2 + (1|Plant_Position/Gall_Number),
                    data = simple.dredge.df, family = "binomial", 
                    na.action = na.fail,
                    control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(platy.mech)
#visreg(platy.mech, xvar = "Iteomyia_count", by = "gall_individuals", scale = "response")
dredge.platy <- dredge(platy.mech, trace = FALSE, rank = "AICc")
get.platy <- get.models(dredge.platy, subset = cumsum(weight) <= 0.95)
avg.platy <- model.avg(get.platy)
summary(avg.platy)

#summary(platy.mech)
overdisp_fun_GLMM(igp.mech)
#drop1(igp.mech, test = "Chisq")
visreg(igp.mech, xvar = "Iteomyia_count", by = "Treatment.focus", scale = "response") # per-capita strength of IGP is higher at greater productivity. However, abundance of platy (IGPrey) likely is higher at greater productivity.
visreg(igp.mech, xvar = "Gall_Height_mm", by = "Treatment.focus", scale = "response") # per-capita strength of IGP is higher on smaller galls. So larger galls provide a bit of a refuge from IGP?
visreg(igp.mech, xvar = "gall_individuals", by = "Treatment.focus", scale = "response") # per-capita strength of platygaster parasitism increases with the number of individuals in a gall (local productivity)


platy.ectos.mech <- glmer(platy.ectos > 0 ~ (Iteomyia_count + gall_individuals + Gall_Height_mm)^2 + (1|Plant_Position/Gall_Number),
                          data = complex.dredge.df, 
                          family = "binomial",
                          na.action = na.fail,
                          control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(platy.ectos.mech)

dredge.platy.ectos <- dredge(platy.ectos.mech, trace = FALSE, rank = "AICc")
get.platy.ectos <- get.models(dredge.platy.ectos, subset = cumsum(weight) <= 0.95)
avg.platy.ectos <- model.avg(get.platy.ectos)
summary(avg.platy.ectos)


## Ectoparasitism ----
# Predict that treatment will dramatically reduced rates of ectoparasitism
ecto <- glmer(ectos > 0 ~ Treatment.focus + Genotype + (1|Plant_Position/Gall_Number),
              data = ind.tree.IGP, 
              family = "binomial",
              control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(ecto)
overdisp_fun_GLMM(ecto)
drop1(ecto, test = "Chisq") # why is model is nearly unidentifiable when I run the Chi sq test?
visreg(ecto, xvar = "Treatment.focus", by = "Genotype", scale = "response")

#ecto.dredge <- ind.tree.IGP %>%
# filter(Treatment.focus == "Control")

ecto.mech <- glmer(ectos > 0 ~ (Iteomyia_count + gall_individuals + Gall_Height_mm)^2 + (1|Genotype/Plant_Position/Gall_Number),
                   data = complex.dredge.df, family = "binomial",
                   na.action = na.fail, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(ecto.mech)
#visreg(ecto.mech, xvar = "gall_individuals", by = "Gall_Height_mm", scale = "response")
dredge.ecto <- dredge(ecto.mech, trace = FALSE, rank = "AICc")
get.ecto <- get.models(dredge.ecto, subset = cumsum(weight) <= 0.95)
avg.ecto <- model.avg(get.ecto)
summary(avg.ecto)

anova(glmer(ectos > 0 ~ gall_individuals*Gall_Height_mm + (1|Genotype/Plant_Position/Gall_Number),
            data = ecto.dredge, family = "binomial",
            na.action = na.fail),
      glmer(ectos > 0 ~ 1 + (1|Genotype/Plant_Position/Gall_Number),
            data = ecto.dredge, family = "binomial",
            na.action = na.fail))



ecto.mech <- glmer(ectos > 0 ~ gall_individuals*Gall_Height_mm + Iteomyia_count + (1|Plant_Position/Gall_Number),
                   data = filter(ind.tree.cont, gall_individuals > 0, Gall_Height_mm > 0),
                   family = "binomial",
                   control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)), na.action = na.fail)
#dredge(ecto.mech)
summary(ecto.mech)
overdisp_fun_GLMM(ecto.mech)
#drop1(ecto.mech, test = "Chisq") # likelihood ratio tests
visreg(ecto.mech, by = "gall_individuals", xvar = "Gall_Height_mm", scale = "response") # ectoparasitism decreases with increasing number of gall individuals. Perhaps this provides the refuge from predation for Platygaster? 
visreg(ecto.mech, xvar = "Gall_Height_mm", by = "Treatment.focus", scale = "response") # ectoparasitism decreases on larger galls (same as platy)
visreg(ecto.mech, xvar = "Iteomyia_count", by = "Treatment.focus", scale = "response")
#confint(profile(ecto.mech))


## unexpected result that pupa survival was lower at reduced densities (according to previous model), may just ignore this result...
pupa.red <- glmer(pupa > 0 ~ Treatment.focus*Iteomyia_count + gall_individuals + Gall_Height_mm + (1|Plant_Position/Gall_Number),
                  data = ind.tree.red, family = "binomial",
                  control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(pupa.red)
overdisp_fun_GLMM(pupa.red)
drop1(pupa.red, test = "Chisq") # non-significant GxE
visreg(pupa.red, by = "Treatment.focus", xvar = "gall_individuals", scale = "response")
#confint(profile(pupa.red))