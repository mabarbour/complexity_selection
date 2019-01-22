---
title: "Reproduce analyses reported in: *Phenotypic evolution is more constrained in simpler food webs*"
author: "Matthew A. Barbour"
date: "2019-01-22"
output:  
  prettydoc::html_pretty:
    theme: architect
    highlight: github
fontsize: 11pt
---

## Setup

Suppress messages, warnings, and errors for document aesthetics.


```r
knitr::opts_chunk$set(message=FALSE, warning=FALSE, error=FALSE)
```

Load required libraries.


```r
library(car)        # for homogeneity of variance test
library(tidyverse)  # for managing data
library(cowplot)    # pretty default ggplots
library(broom)      # for tidying multiple linear models
library(viridis)    # for color palette
library(lme4)       # for generalized linear mixed models
library(latex2exp)  # using Latex notation for axes labels
library(stringr)
library(MVN)        # assessing multivariate normality assumption
```

Write custom functions to avoid repetitive code.


```r
# Transform logit to probability
inverse_logit <- function(x) exp(x)/(1+exp(x))

# Get confidence intervals from bootstrapped samples
conf.high <- function(x) quantile(x, probs = 0.975)
conf.low <- function(x) quantile(x, probs = 0.025)

# Get data for plotting bootstrapped estimates in Figures 2, 3, and 4 of main text
bootstrap_fitness <- function(logistic_model, newdata, bootstraps, intervals=c(0.025,0.975)){
  
  # Absolute fitness
  get_absolute_fitness <- function(.) predict(., newdata=newdata, type="response", re.form=~0) 
  
  if(is.null(bootstraps) == FALSE){

    get_bootstraps_absolute_fitness <-bootMer(logistic_model, FUN = get_absolute_fitness, nsim=bootstraps, parallel="multicore", ncpus=32, seed=34)
  
   get_intervals_absolute_fitness <- apply(get_bootstraps_absolute_fitness$t, 2, function(x) x[order(x)][c(round(bootstraps*intervals[1],0),round(bootstraps*intervals[2],0))])
   
  absolute_fitness_df <- data.frame(newdata, 
                    average = get_absolute_fitness(logistic_model), 
                    lower = get_intervals_absolute_fitness[1, ],
                    upper = get_intervals_absolute_fitness[2, ],
                    t(get_bootstraps_absolute_fitness$t))
  } else {
    absolute_fitness_df <- data.frame(newdata, 
                    average = get_absolute_fitness(logistic_model))
  }
  
  
  # Relative fitness
  get_relative_fitness <- function(.) get_absolute_fitness(.)/mean(get_absolute_fitness(.))
  
  if(is.null(bootstraps) == FALSE){

    get_bootstraps_relative_fitness <-bootMer(logistic_model, FUN = get_relative_fitness, nsim=bootstraps, parallel="multicore", ncpus=32, seed=34)
  
  get_intervals_relative_fitness <- apply(get_bootstraps_relative_fitness$t, 2, function(x) x[order(x)][c(round(bootstraps*intervals[1],0),round(bootstraps*intervals[2],0))])
    
  relative_fitness_df <- data.frame(newdata, 
                    average = get_relative_fitness(logistic_model), 
                    lower = get_intervals_relative_fitness[1, ],
                    upper = get_intervals_relative_fitness[2, ],
                    t(get_bootstraps_relative_fitness$t))
  } else {
      relative_fitness_df <- data.frame(newdata, 
                    average = get_relative_fitness(logistic_model))
  }
  
  # Organize data
  return(list(absolute_fitness = absolute_fitness_df, relative_fitness = relative_fitness_df))
}
```

Set parameters for analyses and plots.


```r
# Number of bootstrap replicates
n_boots_analysis <- 1000
n_boots_plots <- 100

# Color-blind friendly palette: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/; 
# order: orange, light blue, green, yellow, dark blue, red, pink
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#999999")
complex_color <- cbPalette[6] 
simple_color <- cbPalette[5] 
treatment_colors <- c(complex_color, simple_color)
```

Read and tidy data for analyses.


```r
gall_selection.df <- read_csv("../data/gall_selection_manuscript_dataset.csv") %>%
  # convert appropriate variables to characters instead of integers
  mutate(Plant_Position = as.character(Plant_Position),
         Gall_Number = as.character(Gall_Number)) %>%
  unite(Gall_ID, Gall_Number, Gall_Letter, remove = FALSE) %>%
  
  # eliminate unknown sources of mortality
  filter(egg.ptoid > 0 | larval.ptoid > 0 | pupa > 0) %>% 
  # excluding any larva that were parasitized by an ectoparasitoid in the exclusion experiment.
  filter(Treatment == "Ectoparasitoid exclusion" & larval.ptoid < 1 | Treatment == "Control") %>% 
  mutate(Foodweb = ifelse(Treatment=="Control","Complex","Simple"), # rename Treatment to match presentation in manuscript
         gall_survival = as.numeric(ifelse(pupa > 0, 1, 0)),
         egg.ptoid = as.numeric(ifelse(egg.ptoid > 0, 1, 0)), # there were a couple occassions where two egg parasitoids were found in a single chamber. I set these to 1 to permit an analysis with a binomial model.
         
         # transform and scale traits to approximate multivariate normal assumption
         sc.Diam = as.numeric(scale(chamber.diameter_mm)), 
         sc.log.Clutch = as.numeric(scale(log(clutch.size_individuals))), 
         sc.sqrt.Pref = as.numeric(scale(sqrt(oviposition.preference_individuals.per.100.shoots)))) 
```

\  

## Evaluating assumption of multivariate normality

We used graphical checks to evaluate whether our transformations of trait values resulted in a multivariate normal distribution. The histograms below show that our transformations resulted in approximately normal distributions for each phenotypic trait. Note also that in the multivariate quantile-quantile (Q-Q) plot, most points fall along the expected line, suggesting that our transformations provide a reasonable approximation of a multivariate normal distribution. 


```r
mvn_univariate_hist <- mvn(select(gall_selection.df, sc.Diam, sc.log.Clutch, sc.sqrt.Pref), univariatePlot = "histogram")
```

![\label{fig:Univariate_histograms}Histograms of each phenotypic trait after transformation. The red line illustrates a normal distribution.](manuscript_files/figure-latex/Univariate_histograms-1.pdf) 


```r
mvn_multivariate_qq <- mvn(select(gall_selection.df, sc.Diam, sc.log.Clutch, sc.sqrt.Pref), multivariatePlot = "qq")
```

![Multivariate quantile-quantile (Q-Q) plot to assess deviations from multivariate normality (black line).](manuscript_files/figure-latex/Multivariate_QQ-1.pdf) 

\ 

## Gall trait-fitness relationships: effect of food-web treatment

We write the model in a way that independently estimates the effect of food-web treatment, each trait, and all two-way and three-way statistical interactions, on larval survival. The resulting estimates and confidence intervals are useful for determining whether trait-fitness relationships differ from zero, but not whether they differ between food-web treatments. For the later, we calculate the differences between each food-web treatment from the bootstrapped samples.


```r
foodweb_model <- glmer(
  gall_survival ~ 
    -1 + Foodweb + 
    Foodweb:(sc.Diam + sc.log.Clutch + sc.sqrt.Pref) +
    Foodweb:(I(sc.Diam^2) + I(sc.log.Clutch^2) + I(sc.sqrt.Pref^2)) +
    Foodweb:(sc.Diam:sc.log.Clutch + sc.Diam:sc.sqrt.Pref + sc.log.Clutch:sc.sqrt.Pref) +
    (1|Genotype/Plant_Position/Gall_Number),
  data = gall_selection.df,
  family = binomial(link = logit), control=glmerControl(optimizer = "bobyqa"))
```

We used parametric bootstrapping to calculate 95% confidence intervals of trait-fitness relationships.


```r
boot_foodweb_model <- bootMer(foodweb_model, FUN = fixef, nsim=n_boots_analysis, parallel="multicore", ncpus=32, seed=34)
```

In order to reliably quantify linear trait-fitness relationships, we remove all higher-order terms from the model [Stinchcombe et al. 2008](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1558-5646.2008.00449.x).


```r
linear_foodweb_model <- glmer(
  gall_survival ~ 
    -1 + Foodweb + 
    Foodweb:(sc.Diam + sc.log.Clutch + sc.sqrt.Pref) +
    (1|Genotype/Plant_Position/Gall_Number),
  data = gall_selection.df,
  family = binomial(link = logit), control=glmerControl(optimizer = "bobyqa"))
```

And refit using parametric bootstrapping.


```r
boot_linear_foodweb_model <- bootMer(linear_foodweb_model, FUN = fixef, nsim=n_boots_analysis, parallel="multicore", ncpus=32, seed=34)
```

Chamber diameter, but not the other phenotypes, could be influenced by parasitism itself rather than being under natural selection. To estimate this potential bias, we subset our data to only include multi-chambered galls where there was variability in larval survival. We then fit a reduced model to estimate the bias in the logistic regression coefficient of chamber diameter in each food web.


```r
# Subset data
biased_foodweb_df <- gall_selection.df %>%
  group_by(Foodweb, Gall_Number) %>%
  mutate(mean_survival = mean(gall_survival)) %>%
  filter(mean_survival > 0, mean_survival < 1) %>%
  ungroup()

# Fit linear model with only chamber diameter 
biased_foodweb_model <- glmer(
  gall_survival ~ -1 + Foodweb + 
    Foodweb:sc.Diam + 
    (1|Genotype/Plant_Position/Gall_Number),
  data = biased_foodweb_df,
  family = binomial(link = logit), control=glmerControl(optimizer = "bobyqa"))

# Accuracy of confidence intervals is not a priority here, so we just use the asymptotic estimates rather than parametric bootstrapping.
biased_foodweb_confint <- tidy(biased_foodweb_model, conf.int=TRUE) %>% filter(group=="fixed")
```

\begin{table}[t]

\caption{\label{tab:biased foodweb table}Estimates of bias in trait-fitness relationship of chamber diameter in each food-web treatment.}
\centering
\begin{tabular}{l|r|r|r}
\hline
term & estimate & conf.low & conf.high\\
\hline
FoodwebComplex:sc.Diam & 0.360 & 0.054 & 0.667\\
\hline
FoodwebSimple:sc.Diam & 0.416 & 0.011 & 0.821\\
\hline
\end{tabular}
\end{table}

\ 

We gather regression coefficients from the full and reduced model to display trait-fitness relationships. For chamber diameter, we subtracted the coefficient due to bias to better approximate the trait-fitness relationship. We also used the bootstrapped samples to calculate the differences in estimates between treatments (labelled **Diff**). 


```r
foodweb_bind <- bind_cols(
  as_tibble(boot_linear_foodweb_model$t),
  select(as_tibble(boot_foodweb_model$t), `FoodwebComplex:I(sc.Diam^2)`:`FoodwebSimple:sc.log.Clutch:sc.sqrt.Pref`) # only retain nonlinear coefficients from full model
)

# Adjust coefficients with chamber diameter
foodweb_adj_Diam <- foodweb_bind %>%
  mutate(`FoodwebComplex:sc.Diam` = `FoodwebComplex:sc.Diam` - fixef(biased_foodweb_model)["FoodwebComplex:sc.Diam"],
         `FoodwebSimple:sc.Diam` = `FoodwebSimple:sc.Diam` - fixef(biased_foodweb_model)["FoodwebSimple:sc.Diam"])

# Calculate differences between food-web treatments, then estimate means and confidence intervals
foodweb_alphas <- foodweb_adj_Diam %>%
  mutate(`FoodwebSimple:(Intercept)`=inverse_logit(FoodwebSimple), `FoodwebComplex:(Intercept)`=inverse_logit(FoodwebComplex)) %>% # logit transform to make intercept interpretable as a probability
  select(-FoodwebSimple, -FoodwebComplex) %>% # translated into (Intercept) terms
  mutate(`FoodwebDiff:(Intercept)` = `FoodwebSimple:(Intercept)` - `FoodwebComplex:(Intercept)`,
         `FoodwebDiff:sc.Diam` = `FoodwebSimple:sc.Diam` - `FoodwebComplex:sc.Diam`,
         `FoodwebDiff:sc.log.Clutch` = `FoodwebSimple:sc.log.Clutch` - `FoodwebComplex:sc.log.Clutch`,
         `FoodwebDiff:sc.sqrt.Pref` = `FoodwebSimple:sc.sqrt.Pref` - `FoodwebComplex:sc.sqrt.Pref`,
         `FoodwebDiff:I(sc.Diam^2)` = `FoodwebSimple:I(sc.Diam^2)` - `FoodwebComplex:I(sc.Diam^2)`,
         `FoodwebDiff:I(sc.log.Clutch^2)` = `FoodwebSimple:I(sc.log.Clutch^2)` - `FoodwebComplex:I(sc.log.Clutch^2)`,
         `FoodwebDiff:I(sc.sqrt.Pref^2)` = `FoodwebSimple:I(sc.sqrt.Pref^2)` - `FoodwebComplex:I(sc.sqrt.Pref^2)`,
         `FoodwebDiff:sc.Diam:sc.log.Clutch` = `FoodwebSimple:sc.Diam:sc.log.Clutch` - `FoodwebComplex:sc.Diam:sc.log.Clutch`,
         `FoodwebDiff:sc.Diam:sc.sqrt.Pref` = `FoodwebSimple:sc.Diam:sc.sqrt.Pref` - `FoodwebComplex:sc.Diam:sc.sqrt.Pref`,
         `FoodwebDiff:sc.log.Clutch:sc.sqrt.Pref` = `FoodwebSimple:sc.log.Clutch:sc.sqrt.Pref` - `FoodwebComplex:sc.log.Clutch:sc.sqrt.Pref`) %>%
  summarise_all(funs(mean, conf.high, conf.low))

# Tidy up estimates
tidy_foodweb_alphas <- foodweb_alphas %>%
  t() %>%
  as.data.frame() %>%
  transmute(rows = rownames(.), value=V1) %>%
  separate(col = rows, into = c("term","estimate"), sep="_") %>%
  separate(col = term, into = c("type","term_1","term_2"), sep=":", extra="drop") %>%
  mutate(term = ifelse(is.na(term_2), term_1, paste(term_1,term_2,sep=":"))) %>%
  separate(col = type, into = c("extra","type"), sep = 7) %>%
  select(type, term, estimate, value) %>%
  spread(estimate, value) %>%
  mutate(P_cutoff = ifelse(conf.high*conf.low < 0, "","*"),
         coefficient_type = ifelse(term=="(Intercept)","Mean fitness",
                                   ifelse(grepl("\\^2", .$term)==TRUE, "Quadratic",
                                          ifelse(grepl(":", .$term)==TRUE, "Correlational","Linear"))))
```

\  

\begin{table}[t]

\caption{\label{tab:Table of food web coefficients}Table of mean fitness and trait-fitness relationships in each food-web treatment.}
\centering
\begin{tabular}{l|l|l|r|r|r|l}
\hline
term & coefficient\_type & type & mean & conf.low & conf.high & P\_cutoff\\
\hline
(Intercept) & Mean fitness & Complex & 0.416 & 0.276 & 0.560 & *\\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & Complex & 0.223 & -0.108 & 0.553 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & Complex & -0.090 & -0.447 & 0.297 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & Complex & 0.563 & 0.115 & 1.050 & *\\
\hline
sc.Diam & Linear & Complex & 1.149 & 0.750 & 1.612 & *\\
\hline
sc.Diam:sc.log.Clutch & Correlational & Complex & -0.142 & -0.536 & 0.260 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & Complex & -0.442 & -0.971 & 0.055 & \\
\hline
sc.log.Clutch & Linear & Complex & 0.202 & -0.166 & 0.575 & \\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & Complex & 0.094 & -0.344 & 0.588 & \\
\hline
sc.sqrt.Pref & Linear & Complex & -0.424 & -0.975 & 0.151 & \\
\hline
(Intercept) & Mean fitness & Simple & 0.683 & 0.539 & 0.808 & *\\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & Simple & 0.269 & -0.049 & 0.616 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & Simple & -0.304 & -0.741 & 0.083 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & Simple & 0.044 & -0.378 & 0.468 & \\
\hline
sc.Diam & Linear & Simple & 1.104 & 0.639 & 1.650 & *\\
\hline
sc.Diam:sc.log.Clutch & Correlational & Simple & -0.346 & -0.805 & 0.088 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & Simple & -0.085 & -0.538 & 0.362 & \\
\hline
sc.log.Clutch & Linear & Simple & -0.469 & -0.925 & -0.061 & *\\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & Simple & 0.016 & -0.378 & 0.377 & \\
\hline
sc.sqrt.Pref & Linear & Simple & -0.845 & -1.366 & -0.339 & *\\
\hline
\end{tabular}
\end{table}

\  

\begin{table}[t]

\caption{\label{tab:Table of differences in food web coefficients}Table of differences in mean fitness and trait-fitness relationships between food-web treatments.}
\centering
\begin{tabular}{l|l|r|r|r|l}
\hline
term & coefficient\_type & mean & conf.low & conf.high & P\_cutoff\\
\hline
(Intercept) & Mean fitness & 0.267 & 0.115 & 0.425 & *\\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & 0.047 & -0.424 & 0.529 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & -0.214 & -0.729 & 0.333 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & -0.519 & -1.149 & 0.092 & \\
\hline
sc.Diam & Linear & -0.045 & -0.549 & 0.503 & \\
\hline
sc.Diam:sc.log.Clutch & Correlational & -0.204 & -0.788 & 0.419 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & 0.357 & -0.303 & 1.050 & \\
\hline
sc.log.Clutch & Linear & -0.670 & -1.278 & -0.166 & *\\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & -0.078 & -0.661 & 0.461 & \\
\hline
sc.sqrt.Pref & Linear & -0.421 & -1.152 & 0.274 & \\
\hline
\end{tabular}
\end{table}

\  

The table is good for details, but it is easier to see the results as a figure.


```r
# Plot helpers
dodge_width <- position_dodge(width=0.5)
y_limits_foodweb_coefs <- c(min(filter(tidy_foodweb_alphas, type!="Diff")$conf.low), 
                            max(filter(tidy_foodweb_alphas, type!="Diff")$conf.high))
P_asterisk_foodweb_coefs <- 1.5
P_size <- 8
legend_foodweb_labels <- c("Complex","Simple")
point.size <- 3

# Plot mean fitness
plot_foodweb_mean_fitness <- tidy_foodweb_alphas %>%
  filter(type!="Diff", term=="(Intercept)") %>%
  ggplot(., aes(x=type)) +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_foodweb_alphas, type=="Diff", term=="(Intercept)"), aes(y=0.8, x=1.5, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Larval survival", limits=c(0,1), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
  xlab("Food web") +
  facet_wrap(~coefficient_type) +
  scale_fill_manual(values=cbPalette[c(6,5)], labels=legend_foodweb_labels) +
  scale_color_manual(values=cbPalette[c(6,5)], labels=legend_foodweb_labels) 

# Plot linear trait-fitness relationaships
plot_foodweb_linear_coefs <- tidy_foodweb_alphas %>%
  filter(type!="Diff", term %in% c("sc.Diam","sc.log.Clutch","sc.sqrt.Pref")) %>%
  ggplot(., aes(x=term)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_foodweb_alphas, type=="Diff", term %in% c("sc.Diam","sc.log.Clutch","sc.sqrt.Pref")), 
            aes(y=P_asterisk_foodweb_coefs, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Coefficient", limits=y_limits_foodweb_coefs) +
  scale_x_discrete(name="", labels=parse(text=c("alpha[Diam]","alpha[Clutch]","alpha[Pref]"))) +
  facet_wrap(~coefficient_type) +
  scale_fill_manual(values=cbPalette[c(6,5)], labels=legend_foodweb_labels) +
  scale_color_manual(values=cbPalette[c(6,5)], labels=legend_foodweb_labels)

# Plot quadratic trait-fitness relationships
plot_foodweb_quadratic_coefs <- tidy_foodweb_alphas %>%
  filter(type!="Diff", term %in% c("I(sc.Diam^2)","I(sc.log.Clutch^2)","I(sc.sqrt.Pref^2)")) %>%
  ggplot(., aes(x=term, color=type)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_foodweb_alphas, type=="Diff", term %in% c("I(sc.Diam^2)","I(sc.log.Clutch^2)","I(sc.sqrt.Pref^2)")), 
            aes(y=P_asterisk_foodweb_coefs, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Coefficient", limits=y_limits_foodweb_coefs) +
  scale_x_discrete(name="", labels=parse(text=c("alpha[Diam:Diam]","alpha[Clutch:Clutch]","alpha[Pref:Pref]"))) +
  facet_wrap(~coefficient_type) +
  scale_fill_manual(values=cbPalette[c(6,5)], labels=legend_foodweb_labels) +
  scale_color_manual(values=cbPalette[c(6,5)], labels=legend_foodweb_labels)

# Plot correlational trait-fitness relationaships
plot_foodweb_correlational_coefs <- tidy_foodweb_alphas %>%
  filter(type!="Diff", term %in% c("sc.Diam:sc.log.Clutch","sc.Diam:sc.sqrt.Pref","sc.log.Clutch:sc.sqrt.Pref")) %>%
  ggplot(., aes(x=term, color=type)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_foodweb_alphas, type=="Diff", term %in% c("sc.Diam:sc.log.Clutch","sc.Diam:sc.sqrt.Pref","sc.log.Clutch:sc.sqrt.Pref")), 
            aes(y=P_asterisk_foodweb_coefs, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Coefficient", limits=y_limits_foodweb_coefs) +
  scale_x_discrete(name="", labels=parse(text=c("alpha[Diam:Clutch]","alpha[Diam:Pref]","alpha[Clutch:Pref]"))) +
  facet_wrap(~coefficient_type) +
  scale_fill_manual(values=cbPalette[c(6,5)], labels=legend_foodweb_labels) +
  scale_color_manual(values=cbPalette[c(6,5)], labels=legend_foodweb_labels) 

# Format plots for manuscript
plot_foodweb_coefs <- plot_grid(
  plot_foodweb_mean_fitness + theme(legend.position = "none"),
  plot_foodweb_linear_coefs + theme(legend.position = "none"),
  plot_foodweb_quadratic_coefs + theme(legend.position = "none"),
  plot_foodweb_correlational_coefs + theme(legend.position = "none"), 
  nrow=2, align="hv", labels = "")

figure_foodweb_coefs <- plot_foodweb_coefs 
figure_foodweb_coefs
```

![Plot of mean fitness differences between treatments and trait-fitness relationships. Asterisks denote significant differences (*P* < 0.05) between treatments.](manuscript_files/figure-latex/Plot food-web coefficients-1.pdf) 

\  

## Gall trait selection gradients

Using the methods proposed by [Janzen and Stern, 1998](https://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.1998.tb02237.x), we estimated selection gradients from our trait-fitness relationships.


```r
# Estimate mean fitness and mean "brackets" for each food-web treatment (see Janzen and Stern 1998 equation 4 for details about "brackets")
foodweb_complex_predict <- predict(foodweb_model, newdata=filter(gall_selection.df, Foodweb=="Complex"), type="response")
foodweb_complex_mean_brackets <- mean(foodweb_complex_predict * (1 - foodweb_complex_predict))
foodweb_complex_mean_fitness <- mean(foodweb_complex_predict)

foodweb_simple_predict <- predict(foodweb_model, newdata=filter(gall_selection.df, Foodweb=="Simple"), type="response")
foodweb_simple_mean_brackets <- mean(foodweb_simple_predict * (1 - foodweb_simple_predict))
foodweb_simple_mean_fitness <- mean(foodweb_simple_predict)

foodweb_complex_gradient <- function(x) foodweb_complex_mean_brackets * x / foodweb_complex_mean_fitness
foodweb_simple_gradient <- function(x) foodweb_simple_mean_brackets * x / foodweb_simple_mean_fitness

foodweb_complex_raw_gradients <- select(foodweb_adj_Diam, starts_with("FoodwebComplex:")) %>%
  transmute_all(funs(foodweb_complex_gradient))

foodweb_simple_raw_gradients <- select(foodweb_adj_Diam, starts_with("FoodwebSimple:")) %>%
  transmute_all(funs(foodweb_simple_gradient))

# Calculate differences between food-web treatments, then estimate means and confidence intervals
foodweb_grads <- bind_cols(foodweb_complex_raw_gradients, foodweb_simple_raw_gradients) %>%
  mutate(`FoodwebDiff:sc.Diam` = `FoodwebSimple:sc.Diam` - `FoodwebComplex:sc.Diam`,
         `FoodwebDiff:sc.log.Clutch` = `FoodwebSimple:sc.log.Clutch` - `FoodwebComplex:sc.log.Clutch`,
         `FoodwebDiff:sc.sqrt.Pref` = `FoodwebSimple:sc.sqrt.Pref` - `FoodwebComplex:sc.sqrt.Pref`,
         `FoodwebDiff:I(sc.Diam^2)` = `FoodwebSimple:I(sc.Diam^2)` - `FoodwebComplex:I(sc.Diam^2)`,
         `FoodwebDiff:I(sc.log.Clutch^2)` = `FoodwebSimple:I(sc.log.Clutch^2)` - `FoodwebComplex:I(sc.log.Clutch^2)`,
         `FoodwebDiff:I(sc.sqrt.Pref^2)` = `FoodwebSimple:I(sc.sqrt.Pref^2)` - `FoodwebComplex:I(sc.sqrt.Pref^2)`,
         `FoodwebDiff:sc.Diam:sc.log.Clutch` = `FoodwebSimple:sc.Diam:sc.log.Clutch` - `FoodwebComplex:sc.Diam:sc.log.Clutch`,
         `FoodwebDiff:sc.Diam:sc.sqrt.Pref` = `FoodwebSimple:sc.Diam:sc.sqrt.Pref` - `FoodwebComplex:sc.Diam:sc.sqrt.Pref`,
         `FoodwebDiff:sc.log.Clutch:sc.sqrt.Pref` = `FoodwebSimple:sc.log.Clutch:sc.sqrt.Pref` - `FoodwebComplex:sc.log.Clutch:sc.sqrt.Pref`) %>%
  summarise_all(funs(mean, conf.high, conf.low))


# Tidy up estimates
tidy_foodweb_grads <- foodweb_grads %>%
  t() %>%
  as.data.frame() %>%
  transmute(rows = rownames(.), value=V1) %>%
  separate(col = rows, into = c("term","estimate"), sep="_") %>%
  separate(col = term, into = c("type","term_1","term_2"), sep=":", extra="drop") %>%
  mutate(term = ifelse(is.na(term_2), term_1, paste(term_1,term_2,sep=":"))) %>%
  separate(col = type, into = c("extra","type"), sep = 7) %>%
  select(type, term, estimate, value) %>%
  spread(estimate, value) %>%
  mutate(P_cutoff = ifelse(conf.high*conf.low < 0, "","*"),
         gradient_type = ifelse(grepl("\\^2", .$term)==TRUE, "Quadratic",
                                ifelse(grepl(":", .$term)==TRUE, "Correlational","Directional")),
         mean = ifelse(gradient_type=="Quadratic", 2*mean, mean),
         conf.low = ifelse(gradient_type=="Quadratic", 2*conf.low, conf.low),
         conf.high = ifelse(gradient_type=="Quadratic", 2*conf.high, conf.high))
```

\  

Reproduce numbers reported in **Table 1** of the manuscript.

\begin{table}[t]

\caption{\label{tab:Table of food web selection gradients}Table of selection gradients acting on gall traits in each food-web treatment.}
\centering
\begin{tabular}{l|l|l|r|r|r|l}
\hline
term & gradient\_type & type & mean & conf.low & conf.high & P\_cutoff\\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & Complex & 0.133 & -0.065 & 0.332 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & Complex & -0.054 & -0.268 & 0.178 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & Complex & 0.338 & 0.069 & 0.630 & *\\
\hline
sc.Diam & Directional & Complex & 0.344 & 0.225 & 0.483 & *\\
\hline
sc.Diam:sc.log.Clutch & Correlational & Complex & -0.043 & -0.161 & 0.078 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & Complex & -0.133 & -0.291 & 0.017 & \\
\hline
sc.log.Clutch & Directional & Complex & 0.060 & -0.050 & 0.172 & \\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & Complex & 0.028 & -0.103 & 0.176 & \\
\hline
sc.sqrt.Pref & Directional & Complex & -0.127 & -0.292 & 0.045 & \\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & Simple & 0.101 & -0.019 & 0.232 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & Simple & -0.114 & -0.279 & 0.031 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & Simple & 0.017 & -0.142 & 0.176 & \\
\hline
sc.Diam & Directional & Simple & 0.208 & 0.120 & 0.310 & *\\
\hline
sc.Diam:sc.log.Clutch & Correlational & Simple & -0.065 & -0.151 & 0.017 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & Simple & -0.016 & -0.101 & 0.068 & \\
\hline
sc.log.Clutch & Directional & Simple & -0.088 & -0.174 & -0.011 & *\\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & Simple & 0.003 & -0.071 & 0.071 & \\
\hline
sc.sqrt.Pref & Directional & Simple & -0.159 & -0.257 & -0.064 & *\\
\hline
\end{tabular}
\end{table}

\  

\begin{table}[t]

\caption{\label{tab:Table of differences in selection gradients between food-web treatments}Table of differences in selection gradients acting on gall traits between food-web treatment.}
\centering
\begin{tabular}{l|l|r|r|r|l}
\hline
term & gradient\_type & mean & conf.low & conf.high & P\_cutoff\\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & -0.032 & -0.268 & 0.201 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & -0.060 & -0.317 & 0.203 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & -0.321 & -0.640 & -0.009 & *\\
\hline
sc.Diam & Directional & -0.137 & -0.266 & 0.000 & *\\
\hline
sc.Diam:sc.log.Clutch & Correlational & -0.022 & -0.168 & 0.125 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & 0.117 & -0.050 & 0.305 & \\
\hline
sc.log.Clutch & Directional & -0.149 & -0.293 & -0.031 & *\\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & -0.025 & -0.183 & 0.117 & \\
\hline
sc.sqrt.Pref & Directional & -0.032 & -0.211 & 0.149 & \\
\hline
\end{tabular}
\end{table}

\  

Now visualize the table output as a figure.


```r
# Plot helpers
dodge_width <- position_dodge(width=0.5)
y_limits_foodweb_grads <- c(min(filter(tidy_foodweb_grads, type!="Diff")$conf.low), 
                            max(filter(tidy_foodweb_grads, type!="Diff")$conf.high))
P_asterisk_foodweb_grads <- 0.5
P_size <- 8
legend_foodweb_labels <- c("Complex","Simple")
legend_foodweb_title <- "Food web"
point.size <- 3

# Plot directional selection gradients
plot_foodweb_directional_grads <- tidy_foodweb_grads %>%
  filter(type!="Diff", term %in% c("sc.Diam","sc.log.Clutch","sc.sqrt.Pref")) %>%
  ggplot(., aes(x=term)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_foodweb_grads, type=="Diff", term %in% c("sc.Diam","sc.log.Clutch","sc.sqrt.Pref")), 
            aes(y=P_asterisk_foodweb_grads, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Gradient", limits=y_limits_foodweb_grads) +
  scale_x_discrete(name="", labels=parse(text=c("beta[Diam]","beta[Clutch]","beta[Pref]"))) +
  facet_wrap(~gradient_type) +
  scale_fill_manual(values=cbPalette[c(6,5)], labels=legend_foodweb_labels, name=legend_foodweb_title) +
  scale_color_manual(values=cbPalette[c(6,5)], labels=legend_foodweb_labels, name=legend_foodweb_title)

# Plot quadratic selection gradients
plot_foodweb_quadratic_grads <- tidy_foodweb_grads %>%
  filter(type!="Diff", term %in% c("I(sc.Diam^2)","I(sc.log.Clutch^2)","I(sc.sqrt.Pref^2)")) %>%
  ggplot(., aes(x=term, color=type)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_foodweb_grads, type=="Diff", term %in% c("I(sc.Diam^2)","I(sc.log.Clutch^2)","I(sc.sqrt.Pref^2)")), 
            aes(y=P_asterisk_foodweb_grads, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Gradient", limits=y_limits_foodweb_grads) +
  scale_x_discrete(name="", labels=parse(text=c("gamma[Diam:Diam]","gamma[Clutch:Clutch]","gamma[Pref:Pref]"))) +
  facet_wrap(~gradient_type) +
  scale_fill_manual(values=cbPalette[c(6,5)], labels=legend_foodweb_labels, name=legend_foodweb_title) +
  scale_color_manual(values=cbPalette[c(6,5)], labels=legend_foodweb_labels, name=legend_foodweb_title)

# Plot correlational selection gradients
plot_foodweb_correlational_grads <- tidy_foodweb_grads %>%
  filter(type!="Diff", term %in% c("sc.Diam:sc.log.Clutch","sc.Diam:sc.sqrt.Pref","sc.log.Clutch:sc.sqrt.Pref")) %>%
  ggplot(., aes(x=term, color=type)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_foodweb_grads, type=="Diff", term %in% c("sc.Diam:sc.log.Clutch","sc.Diam:sc.sqrt.Pref","sc.log.Clutch:sc.sqrt.Pref")), 
            aes(y=P_asterisk_foodweb_grads, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Gradient", limits=y_limits_foodweb_grads) +
  scale_x_discrete(name="", labels=parse(text=c("gamma[Diam:Clutch]","gamma[Diam:Pref]","gamma[Clutch:Pref]"))) +
  facet_wrap(~gradient_type) +
  scale_fill_manual(values=cbPalette[c(6,5)], labels=legend_foodweb_labels, name=legend_foodweb_title) +
  scale_color_manual(values=cbPalette[c(6,5)], labels=legend_foodweb_labels, name=legend_foodweb_title) 

# Get legend
legend_foodweb_grads <- get_legend(plot_foodweb_directional_grads)

# Format plots for manuscript
plot_foodweb_grads <- plot_grid(
  plot_foodweb_directional_grads + theme(legend.position = "none", axis.title = element_blank()),
  plot_foodweb_quadratic_grads + theme(legend.position = "none", axis.title = element_blank()),
  plot_foodweb_correlational_grads + theme(legend.position = "none", axis.title = element_blank()), 
  nrow=2, align="hv", labels = "")
plot_legend_foodweb_grads <- ggdraw(plot_foodweb_grads) + draw_grob(legend_foodweb_grads, x = 0.7, y=-0.2)
y_title_foodweb_grads <- ggdraw() + draw_label("Selection gradient (SDs)", angle = 90)

figure_foodweb_grads <- plot_grid(y_title_foodweb_grads, plot_legend_foodweb_grads, ncol=2, rel_widths = c(0.05,1)) 
figure_foodweb_grads
```

![Plot of selection gradients in each food-web treatment. Asterisks denote significant differences (*P* < 0.05) between food-web treatments.](manuscript_files/figure-latex/Plot food-web gradients-1.pdf) 

\  

## Gall trait-fitness relationships: contribution of egg and larval parasitoids

Our simple food-web treatment allows us to estimate the unique contribution of egg parasitoids to selection on *Iteomyia* traits. To estimate the unique contribution of larval parasitoids, we subset our data so that our complex food-web treatment only contained attack by larval parasitoids (and gall survival). We then fit the same models as previously, including one to estimate bias.


```r
# excludes cases of egg-parasitism from Complex food web
egglarval_df <- filter(gall_selection.df, Foodweb == "Simple" | Foodweb == "Complex" & egg.ptoid < 1) 

# fit model with same structure as foodweb_model
egglarval_model <- update(foodweb_model, data=egglarval_df)

# subset data to estimate bias
biased_egglarval_df <- egglarval_df %>%
  group_by(Foodweb, Gall_Number) %>%
  mutate(mean_survival = mean(gall_survival)) %>%
  filter(mean_survival > 0, mean_survival < 1) %>%
  ungroup()

# fit model with same structure as biased_foodweb_model
biased_egglarval_model <- update(biased_foodweb_model, data=biased_egglarval_df)
```

Parametric bootstrapping to estimate confidence intervals.


```r
boot_egglarval_model <- bootMer(egglarval_model, FUN = fixef, nsim=n_boots_analysis, parallel="multicore", ncpus=32, seed=34)
```

Fit a reduced model to reliably estimate linear trait-fitness relationships.


```r
linear_egglarval_model <- update(linear_foodweb_model, data=egglarval_df)
```

Parametric bootstrapping to estimate confidence intervals.


```r
boot_linear_egglarval_model <- bootMer(linear_egglarval_model, FUN = fixef, nsim=n_boots_analysis, parallel="multicore", ncpus=32, seed=34)
```

Tidy trait-fitness relationships for egg and larval parasitoids from above models.


```r
egglarval_bind <- bind_cols(
  as_tibble(boot_linear_egglarval_model$t),
  select(as_tibble(boot_egglarval_model$t), `FoodwebComplex:I(sc.Diam^2)`:`FoodwebSimple:sc.log.Clutch:sc.sqrt.Pref`) # only retain nonlinear coefficients from full model
)

# Adjust coefficients with chamber diameter
egglarval_adj_Diam <- as_tibble(boot_egglarval_model$t) %>%
  mutate(`FoodwebComplex:sc.Diam` = `FoodwebComplex:sc.Diam` - fixef(biased_egglarval_model)["FoodwebComplex:sc.Diam"],
         `FoodwebSimple:sc.Diam` = `FoodwebSimple:sc.Diam` - fixef(biased_egglarval_model)["FoodwebSimple:sc.Diam"])

# Calculate differences between food-web treatments, then estimate means and confidence intervals
egglarval_alphas <- egglarval_adj_Diam %>%
  mutate(`FoodwebSimple:(Intercept)`=inverse_logit(FoodwebSimple), `FoodwebComplex:(Intercept)`=inverse_logit(FoodwebComplex)) %>% # logit transform to make intercept interpretable as a probability
  select(-FoodwebSimple, -FoodwebComplex) %>% # translated into (Intercept) terms
  mutate(`FoodwebDiff:(Intercept)` = `FoodwebSimple:(Intercept)` - `FoodwebComplex:(Intercept)`,
         `FoodwebDiff:sc.Diam` = `FoodwebSimple:sc.Diam` - `FoodwebComplex:sc.Diam`,
         `FoodwebDiff:sc.log.Clutch` = `FoodwebSimple:sc.log.Clutch` - `FoodwebComplex:sc.log.Clutch`,
         `FoodwebDiff:sc.sqrt.Pref` = `FoodwebSimple:sc.sqrt.Pref` - `FoodwebComplex:sc.sqrt.Pref`,
         `FoodwebDiff:I(sc.Diam^2)` = `FoodwebSimple:I(sc.Diam^2)` - `FoodwebComplex:I(sc.Diam^2)`,
         `FoodwebDiff:I(sc.log.Clutch^2)` = `FoodwebSimple:I(sc.log.Clutch^2)` - `FoodwebComplex:I(sc.log.Clutch^2)`,
         `FoodwebDiff:I(sc.sqrt.Pref^2)` = `FoodwebSimple:I(sc.sqrt.Pref^2)` - `FoodwebComplex:I(sc.sqrt.Pref^2)`,
         `FoodwebDiff:sc.Diam:sc.log.Clutch` = `FoodwebSimple:sc.Diam:sc.log.Clutch` - `FoodwebComplex:sc.Diam:sc.log.Clutch`,
         `FoodwebDiff:sc.Diam:sc.sqrt.Pref` = `FoodwebSimple:sc.Diam:sc.sqrt.Pref` - `FoodwebComplex:sc.Diam:sc.sqrt.Pref`,
         `FoodwebDiff:sc.log.Clutch:sc.sqrt.Pref` = `FoodwebSimple:sc.log.Clutch:sc.sqrt.Pref` - `FoodwebComplex:sc.log.Clutch:sc.sqrt.Pref`) %>%
  summarise_all(funs(mean, conf.high, conf.low))

# Tidy up estimates
tidy_egglarval_alphas <- egglarval_alphas %>%
  t() %>%
  as.data.frame() %>%
  transmute(rows = rownames(.), value=V1) %>%
  separate(col = rows, into = c("term","estimate"), sep="_") %>%
  separate(col = term, into = c("type","term_1","term_2"), sep=":", extra="drop") %>%
  mutate(term = ifelse(is.na(term_2), term_1, paste(term_1,term_2,sep=":"))) %>%
  separate(col = type, into = c("extra","type"), sep = 7) %>%
  select(type, term, estimate, value) %>%
  spread(estimate, value) %>%
  mutate(P_cutoff = ifelse(conf.high*conf.low < 0, "","*"),
         coefficient_type = ifelse(term=="(Intercept)","Mean fitness",
                                   ifelse(grepl("\\^2", .$term)==TRUE, "Quadratic",
                                          ifelse(grepl(":", .$term)==TRUE, "Correlational","Linear"))))
```

\  

\begin{table}[t]

\caption{\label{tab:Table of egg vs. larval coefficients}Table of gall trait-fitness relationships imposed by egg and larval parasitoids. Note that type = Simple = Egg parasitoid and Complex = Larval parasitoid.}
\centering
\begin{tabular}{l|l|l|r|r|r|l}
\hline
term & coefficient\_type & type & mean & conf.low & conf.high & P\_cutoff\\
\hline
(Intercept) & Mean fitness & Complex & 0.678 & 0.508 & 0.832 & *\\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & Complex & 0.163 & -0.170 & 0.523 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & Complex & 0.014 & -0.364 & 0.382 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & Complex & 0.467 & -0.042 & 1.073 & \\
\hline
sc.Diam & Linear & Complex & 1.169 & 0.672 & 1.854 & *\\
\hline
sc.Diam:sc.log.Clutch & Correlational & Complex & 0.096 & -0.342 & 0.550 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & Complex & -0.277 & -0.923 & 0.327 & \\
\hline
sc.log.Clutch & Linear & Complex & 0.669 & 0.201 & 1.262 & *\\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & Complex & -0.200 & -0.710 & 0.330 & \\
\hline
sc.sqrt.Pref & Linear & Complex & -0.897 & -1.732 & -0.169 & *\\
\hline
(Intercept) & Mean fitness & Simple & 0.689 & 0.503 & 0.854 & *\\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & Simple & 0.252 & -0.036 & 0.582 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & Simple & -0.301 & -0.669 & 0.056 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & Simple & -0.060 & -0.501 & 0.382 & \\
\hline
sc.Diam & Linear & Simple & 1.072 & 0.634 & 1.648 & *\\
\hline
sc.Diam:sc.log.Clutch & Correlational & Simple & -0.359 & -0.771 & 0.035 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & Simple & -0.112 & -0.481 & 0.269 & \\
\hline
sc.log.Clutch & Linear & Simple & -0.806 & -1.377 & -0.308 & *\\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & Simple & 0.040 & -0.307 & 0.414 & \\
\hline
sc.sqrt.Pref & Linear & Simple & -1.000 & -1.658 & -0.471 & *\\
\hline
\end{tabular}
\end{table}

\ 

\begin{table}[t]

\caption{\label{tab:Table of differences in egg vs. larval coefficients}Table of differences in gall trait-fitness relationships between egg and larval parasitoids.}
\centering
\begin{tabular}{l|l|r|r|r|l}
\hline
term & coefficient\_type & mean & conf.low & conf.high & P\_cutoff\\
\hline
(Intercept) & Mean fitness & 0.010 & -0.192 & 0.206 & \\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & 0.090 & -0.372 & 0.563 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & -0.315 & -0.838 & 0.216 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & -0.527 & -1.269 & 0.155 & \\
\hline
sc.Diam & Linear & -0.097 & -0.793 & 0.542 & \\
\hline
sc.Diam:sc.log.Clutch & Correlational & -0.455 & -1.061 & 0.131 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & 0.165 & -0.543 & 0.917 & \\
\hline
sc.log.Clutch & Linear & -1.474 & -2.276 & -0.790 & *\\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & 0.239 & -0.397 & 0.879 & \\
\hline
sc.sqrt.Pref & Linear & -0.103 & -0.961 & 0.851 & \\
\hline
\end{tabular}
\end{table}

\ 

See the table results in figure form.


```r
# Plot helpers
dodge_width <- position_dodge(width=0.5)
y_limits_egglarval_coefs <- c(min(filter(tidy_egglarval_alphas, type!="Diff")$conf.low), 
                            max(filter(tidy_egglarval_alphas, type!="Diff")$conf.high))
P_asterisk_egglarval_coefs <- 1.3
P_size <- 8
legend_egglarval_labels <- c("Larval","Egg")
point.size <- 3

# Plot mean fitness
plot_egglarval_mean_fitness <- tidy_egglarval_alphas %>%
  filter(type!="Diff", term=="(Intercept)") %>%
  ggplot(., aes(x=type)) +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_egglarval_alphas, type=="Diff", term=="(Intercept)"), aes(y=0.9, x=1.5, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Larval survival", limits=c(0,1), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
  scale_x_discrete(name="Parasitoid", breaks=c("Complex","Simple"), labels=c("Larval","Egg")) +
  facet_wrap(~coefficient_type) +
  scale_fill_manual(values=cbPalette[c(4,5)], labels=legend_egglarval_labels) +
  scale_color_manual(values=cbPalette[c(4,5)], labels=legend_egglarval_labels) 

# Plot linear trait-fitness relationaships
plot_egglarval_linear_coefs <- tidy_egglarval_alphas %>%
  filter(type!="Diff", term %in% c("sc.Diam","sc.log.Clutch","sc.sqrt.Pref")) %>%
  ggplot(., aes(x=term)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_egglarval_alphas, type=="Diff", term %in% c("sc.Diam","sc.log.Clutch","sc.sqrt.Pref")), 
            aes(y=P_asterisk_egglarval_coefs, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Coefficient", limits=y_limits_egglarval_coefs) +
  scale_x_discrete(name="", labels=parse(text=c("alpha[Diam]","alpha[Clutch]","alpha[Pref]"))) +
  facet_wrap(~coefficient_type) +
  scale_fill_manual(values=cbPalette[c(4,5)], labels=legend_egglarval_labels) +
  scale_color_manual(values=cbPalette[c(4,5)], labels=legend_egglarval_labels)

# Plot quadratic trait-fitness relationships
plot_egglarval_quadratic_coefs <- tidy_egglarval_alphas %>%
  filter(type!="Diff", term %in% c("I(sc.Diam^2)","I(sc.log.Clutch^2)","I(sc.sqrt.Pref^2)")) %>%
  ggplot(., aes(x=term, color=type)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_egglarval_alphas, type=="Diff", term %in% c("I(sc.Diam^2)","I(sc.log.Clutch^2)","I(sc.sqrt.Pref^2)")), 
            aes(y=P_asterisk_egglarval_coefs, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Coefficient", limits=y_limits_egglarval_coefs) +
  scale_x_discrete(name="", labels=parse(text=c("alpha[Diam:Diam]","alpha[Clutch:Clutch]","alpha[Pref:Pref]"))) +
  facet_wrap(~coefficient_type) +
  scale_fill_manual(values=cbPalette[c(4,5)], labels=legend_egglarval_labels) +
  scale_color_manual(values=cbPalette[c(4,5)], labels=legend_egglarval_labels)

# Plot correlational trait-fitness relationaships
plot_egglarval_correlational_coefs <- tidy_egglarval_alphas %>%
  filter(type!="Diff", term %in% c("sc.Diam:sc.log.Clutch","sc.Diam:sc.sqrt.Pref","sc.log.Clutch:sc.sqrt.Pref")) %>%
  ggplot(., aes(x=term, color=type)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_egglarval_alphas, type=="Diff", term %in% c("sc.Diam:sc.log.Clutch","sc.Diam:sc.sqrt.Pref","sc.log.Clutch:sc.sqrt.Pref")), 
            aes(y=P_asterisk_egglarval_coefs, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Coefficient", limits=y_limits_egglarval_coefs) +
  scale_x_discrete(name="", labels=parse(text=c("alpha[Diam:Clutch]","alpha[Diam:Pref]","alpha[Clutch:Pref]"))) +
  facet_wrap(~coefficient_type) +
  scale_fill_manual(values=cbPalette[c(4,5)], labels=legend_egglarval_labels) +
  scale_color_manual(values=cbPalette[c(4,5)], labels=legend_egglarval_labels) 

# Format plots for manuscript
plot_egglarval_coefs <- plot_grid(
  plot_egglarval_mean_fitness + theme(legend.position = "none"),
  plot_egglarval_linear_coefs + theme(legend.position = "none"),
  plot_egglarval_quadratic_coefs + theme(legend.position = "none"),
  plot_egglarval_correlational_coefs + theme(legend.position = "none"), 
  nrow=2, align="hv", labels = "")

figure_egglarval_coefs <- plot_egglarval_coefs 
figure_egglarval_coefs
```

![Trait-fitness relationships for egg and larval parasitoids. Asterisks denote significant differences (*P* < 0.05) between parasitoid guilds.](manuscript_files/figure-latex/Plot egglarval coefficients-1.pdf) 

\ 

## Selection gradients acting on gall traits due to egg and larval parasitoids


```r
# Estimate mean fitness and mean "brackets" for each food-web treatment (see Janzen and Stern 1998 equation 4 for details about "brackets")
egglarval_complex_predict <- predict(egglarval_model, newdata=filter(egglarval_df, Foodweb=="Complex"), type="response")
egglarval_complex_mean_brackets <- mean(egglarval_complex_predict * (1 - egglarval_complex_predict))
egglarval_complex_mean_fitness <- mean(egglarval_complex_predict)

egglarval_simple_predict <- predict(egglarval_model, newdata=filter(egglarval_df, Foodweb=="Simple"), type="response")
egglarval_simple_mean_brackets <- mean(egglarval_simple_predict * (1 - egglarval_simple_predict))
egglarval_simple_mean_fitness <- mean(egglarval_simple_predict)

egglarval_complex_gradient <- function(x) egglarval_complex_mean_brackets * x / egglarval_complex_mean_fitness
egglarval_simple_gradient <- function(x) egglarval_simple_mean_brackets * x / egglarval_simple_mean_fitness

egglarval_complex_raw_gradients <- select(egglarval_adj_Diam, starts_with("FoodwebComplex:")) %>%
  transmute_all(funs(egglarval_complex_gradient))

egglarval_simple_raw_gradients <- select(egglarval_adj_Diam, starts_with("FoodwebSimple:")) %>%
  transmute_all(funs(egglarval_simple_gradient))

# Calculate differences between food-web treatments, then estimate means and confidence intervals
egglarval_grads <- bind_cols(egglarval_complex_raw_gradients, egglarval_simple_raw_gradients) %>%
  mutate(`FoodwebDiff:sc.Diam` = `FoodwebSimple:sc.Diam` - `FoodwebComplex:sc.Diam`,
         `FoodwebDiff:sc.log.Clutch` = `FoodwebSimple:sc.log.Clutch` - `FoodwebComplex:sc.log.Clutch`,
         `FoodwebDiff:sc.sqrt.Pref` = `FoodwebSimple:sc.sqrt.Pref` - `FoodwebComplex:sc.sqrt.Pref`,
         `FoodwebDiff:I(sc.Diam^2)` = `FoodwebSimple:I(sc.Diam^2)` - `FoodwebComplex:I(sc.Diam^2)`,
         `FoodwebDiff:I(sc.log.Clutch^2)` = `FoodwebSimple:I(sc.log.Clutch^2)` - `FoodwebComplex:I(sc.log.Clutch^2)`,
         `FoodwebDiff:I(sc.sqrt.Pref^2)` = `FoodwebSimple:I(sc.sqrt.Pref^2)` - `FoodwebComplex:I(sc.sqrt.Pref^2)`,
         `FoodwebDiff:sc.Diam:sc.log.Clutch` = `FoodwebSimple:sc.Diam:sc.log.Clutch` - `FoodwebComplex:sc.Diam:sc.log.Clutch`,
         `FoodwebDiff:sc.Diam:sc.sqrt.Pref` = `FoodwebSimple:sc.Diam:sc.sqrt.Pref` - `FoodwebComplex:sc.Diam:sc.sqrt.Pref`,
         `FoodwebDiff:sc.log.Clutch:sc.sqrt.Pref` = `FoodwebSimple:sc.log.Clutch:sc.sqrt.Pref` - `FoodwebComplex:sc.log.Clutch:sc.sqrt.Pref`) %>%
  summarise_all(funs(mean, conf.high, conf.low))

# Tidy up estimates
tidy_egglarval_grads <- egglarval_grads %>%
  t() %>%
  as.data.frame() %>%
  transmute(rows = rownames(.), value=V1) %>%
  separate(col = rows, into = c("term","estimate"), sep="_") %>%
  separate(col = term, into = c("type","term_1","term_2"), sep=":", extra="drop") %>%
  mutate(term = ifelse(is.na(term_2), term_1, paste(term_1,term_2,sep=":"))) %>%
  separate(col = type, into = c("extra","type"), sep = 7) %>%
  select(type, term, estimate, value) %>%
  spread(estimate, value) %>%
  mutate(P_cutoff = ifelse(conf.high*conf.low < 0, "","*"),
         gradient_type = ifelse(grepl("\\^2", .$term)==TRUE, "Quadratic",
                                ifelse(grepl(":", .$term)==TRUE, "Correlational","Directional")),
         mean = ifelse(gradient_type=="Quadratic", 2*mean, mean),
         conf.low = ifelse(gradient_type=="Quadratic", 2*conf.low, conf.low),
         conf.high = ifelse(gradient_type=="Quadratic", 2*conf.high, conf.high))
```

\  

\begin{table}[t]

\caption{\label{tab:Egg vs. larval selection gradients}Table of selection gradients imposed by egg and larval parasitoids.}
\centering
\begin{tabular}{l|l|l|r|r|r|l}
\hline
term & gradient\_type & type & mean & conf.low & conf.high & P\_cutoff\\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & Complex & 0.063 & -0.066 & 0.202 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & Complex & 0.005 & -0.141 & 0.148 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & Complex & 0.181 & -0.016 & 0.416 & \\
\hline
sc.Diam & Directional & Complex & 0.226 & 0.130 & 0.359 & *\\
\hline
sc.Diam:sc.log.Clutch & Correlational & Complex & 0.019 & -0.066 & 0.106 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & Complex & -0.054 & -0.179 & 0.063 & \\
\hline
sc.log.Clutch & Directional & Complex & 0.129 & 0.039 & 0.244 & *\\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & Complex & -0.039 & -0.138 & 0.064 & \\
\hline
sc.sqrt.Pref & Directional & Complex & -0.174 & -0.336 & -0.033 & *\\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & Simple & 0.100 & -0.014 & 0.230 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & Simple & -0.119 & -0.265 & 0.022 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & Simple & -0.024 & -0.198 & 0.151 & \\
\hline
sc.Diam & Directional & Simple & 0.212 & 0.125 & 0.326 & *\\
\hline
sc.Diam:sc.log.Clutch & Correlational & Simple & -0.071 & -0.153 & 0.007 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & Simple & -0.022 & -0.095 & 0.053 & \\
\hline
sc.log.Clutch & Directional & Simple & -0.159 & -0.272 & -0.061 & *\\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & Simple & 0.008 & -0.061 & 0.082 & \\
\hline
sc.sqrt.Pref & Directional & Simple & -0.198 & -0.328 & -0.093 & *\\
\hline
\end{tabular}
\end{table}

\  

\begin{table}[t]

\caption{\label{tab:Table of differences in selection gradients between egg and larval parasitoids}Table of differences in selection gradients imposed by egg and larval parasitoids.}
\centering
\begin{tabular}{l|l|r|r|r|l}
\hline
term & gradient\_type & mean & conf.low & conf.high & P\_cutoff\\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & 0.037 & -0.143 & 0.222 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & -0.125 & -0.329 & 0.082 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & -0.205 & -0.493 & 0.063 & \\
\hline
sc.Diam & Directional & -0.014 & -0.148 & 0.110 & \\
\hline
sc.Diam:sc.log.Clutch & Correlational & -0.090 & -0.208 & 0.026 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & 0.031 & -0.107 & 0.178 & \\
\hline
sc.log.Clutch & Directional & -0.289 & -0.445 & -0.155 & *\\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & 0.046 & -0.077 & 0.171 & \\
\hline
sc.sqrt.Pref & Directional & -0.024 & -0.191 & 0.163 & \\
\hline
\end{tabular}
\end{table}

\  

See table results in figure form.


```r
# Plot helpers
dodge_width <- position_dodge(width=0.5)
y_limits_egglarval_grads <- c(min(filter(tidy_egglarval_grads, type!="Diff")$conf.low), 
                            max(filter(tidy_egglarval_grads, type!="Diff")$conf.high))
P_asterisk_egglarval_grads <- 0.25
P_size <- 8
legend_egglarval_labels <- c("Larval","Egg")
legend_egglarval_title <- "Parasitoid"
point.size <- 3

# Plot directional selection gradients
plot_egglarval_directional_grads <- tidy_egglarval_grads %>%
  filter(type!="Diff", term %in% c("sc.Diam","sc.log.Clutch","sc.sqrt.Pref")) %>%
  ggplot(., aes(x=term)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_egglarval_grads, type=="Diff", term %in% c("sc.Diam","sc.log.Clutch","sc.sqrt.Pref")), 
            aes(y=P_asterisk_egglarval_grads, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Gradient", limits=y_limits_egglarval_grads) +
  scale_x_discrete(name="", labels=parse(text=c("beta[Diam]","beta[Clutch]","beta[Pref]"))) +
  facet_wrap(~gradient_type) +
  scale_fill_manual(values=cbPalette[c(4,5)], labels=legend_egglarval_labels, name=legend_egglarval_title) +
  scale_color_manual(values=cbPalette[c(4,5)], labels=legend_egglarval_labels, name=legend_egglarval_title)

# Plot quadratic selection gradients
plot_egglarval_quadratic_grads <- tidy_egglarval_grads %>%
  filter(type!="Diff", term %in% c("I(sc.Diam^2)","I(sc.log.Clutch^2)","I(sc.sqrt.Pref^2)")) %>%
  ggplot(., aes(x=term, color=type)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_egglarval_grads, type=="Diff", term %in% c("I(sc.Diam^2)","I(sc.log.Clutch^2)","I(sc.sqrt.Pref^2)")), 
            aes(y=P_asterisk_egglarval_grads, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Gradient", limits=y_limits_egglarval_grads) +
  scale_x_discrete(name="", labels=parse(text=c("gamma[Diam:Diam]","gamma[Clutch:Clutch]","gamma[Pref:Pref]"))) +
  facet_wrap(~gradient_type) +
  scale_fill_manual(values=cbPalette[c(4,5)], labels=legend_egglarval_labels, name=legend_egglarval_title) +
  scale_color_manual(values=cbPalette[c(4,5)], labels=legend_egglarval_labels, name=legend_egglarval_title)

# Plot correlational selection gradients
plot_egglarval_correlational_grads <- tidy_egglarval_grads %>%
  filter(type!="Diff", term %in% c("sc.Diam:sc.log.Clutch","sc.Diam:sc.sqrt.Pref","sc.log.Clutch:sc.sqrt.Pref")) %>%
  ggplot(., aes(x=term, color=type)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_egglarval_grads, type=="Diff", term %in% c("sc.Diam:sc.log.Clutch","sc.Diam:sc.sqrt.Pref","sc.log.Clutch:sc.sqrt.Pref")), 
            aes(y=P_asterisk_egglarval_grads, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Gradient", limits=y_limits_egglarval_grads) +
  scale_x_discrete(name="", labels=parse(text=c("gamma[Diam:Clutch]","gamma[Diam:Pref]","gamma[Clutch:Pref]"))) +
  facet_wrap(~gradient_type) +
  scale_fill_manual(values=cbPalette[c(4,5)], labels=legend_egglarval_labels, name=legend_egglarval_title) +
  scale_color_manual(values=cbPalette[c(4,5)], labels=legend_egglarval_labels, name=legend_egglarval_title) 

# Get legend
legend_egglarval_grads <- get_legend(plot_egglarval_directional_grads)

# Format plots for manuscript
plot_egglarval_grads <- plot_grid(
  plot_egglarval_directional_grads + theme(legend.position = "none", axis.title = element_blank()),
  plot_egglarval_quadratic_grads + theme(legend.position = "none", axis.title = element_blank()),
  plot_egglarval_correlational_grads + theme(legend.position = "none", axis.title = element_blank()), 
  nrow=2, align="hv", labels = "")
plot_legend_egglarval_grads <- ggdraw(plot_egglarval_grads) + draw_grob(legend_egglarval_grads, x = 0.7, y=-0.2)
y_title_egglarval_grads <- ggdraw() + draw_label("Selection gradient (SDs)", angle = 90)

figure_egglarval_grads <- plot_grid(y_title_egglarval_grads, plot_legend_egglarval_grads, ncol=2, rel_widths = c(0.05,1)) 
figure_egglarval_grads
```

![Selection gradients imposed by egg and larval parasitoids. Asterisks denote significant differences (*P* < 0.05) between parasitoid guilds.](manuscript_files/figure-latex/Plot egg vs. larval ptoid gradients-1.pdf) 

\  

## Reproduce Figure 2 in manuscript.

Use the fitted model to generate predicted estimates of mean fitness for changes in the mean trait value of a population. I restrict these predictions to +/- 1 SD because we can only reliably estimate the shape of the adaptive landscape near the mean phenotype of the population [Arnold et al., 2001](https://link.springer.com/article/10.1023/A:1013373907708). Other trait values are held constant at the mean phenotype (i.e. trait = 0).


```r
newdata_Diam <- bind_rows(
  expand.grid(Foodweb = "Complex", sc.Diam = seq(-1,1,length.out=1000), sc.log.Clutch = 0, sc.sqrt.Pref = 0),
  expand.grid(Foodweb = "Simple", sc.Diam = seq(-1,1,length.out=1000), sc.log.Clutch = 0, sc.sqrt.Pref = 0))

RF_Diam <- bootstrap_fitness(
  logistic_model = foodweb_model, 
  newdata = newdata_Diam,
  bootstraps=n_boots_plots)
```


```r
newdata_Clutch <- bind_rows(
  expand.grid(Foodweb = "Complex", sc.Diam = 0, sc.log.Clutch = seq(-1,1,length.out=1000), sc.sqrt.Pref = 0),
  expand.grid(Foodweb = "Simple", sc.Diam = 0, sc.log.Clutch = seq(-1,1,length.out=1000), sc.sqrt.Pref = 0))

RF_Clutch <- bootstrap_fitness(
  logistic_model = foodweb_model,
  newdata = newdata_Clutch,
  bootstraps=n_boots_plots)
```


```r
newdata_Pref <- bind_rows(
  expand.grid(Foodweb = "Complex", sc.Diam = 0, sc.log.Clutch = 0, sc.sqrt.Pref = seq(-1,1,length.out=1000)),
  expand.grid(Foodweb = "Simple", sc.Diam = 0, sc.log.Clutch = 0, sc.sqrt.Pref = seq(-1,1,length.out=1000)))

RF_Pref <- bootstrap_fitness(
  logistic_model = foodweb_model, 
  newdata = newdata_Pref,
  bootstraps=n_boots_plots)
```

Create plots for each panel of Figure 2 in main text.


```r
uni_breaks <- c(0.1,0.5,1.0)

# Gall size 
AF_uni_Diam <- RF_Diam$absolute_fitness %>% 
  select(-sc.log.Clutch, -sc.sqrt.Pref) %>%
  gather(ID, absolute_fitness, -sc.Diam, -Foodweb) %>%
  mutate(ID_group = ifelse(ID == "average", "average", "replicate")) %>% 
  ggplot(., aes(x=sc.Diam, y=absolute_fitness, group=interaction(ID,Foodweb), color=Foodweb, alpha=ID_group, size=ID_group)) +
  geom_line() +
  scale_alpha_manual(values=c(1,0.2), guide=FALSE) +
  scale_size_manual(values=c(1.5,0.25), guide=FALSE) +
  xlab("Chamber diameter (SDs)") +
  ylab("Mean larval survival") +
  scale_color_manual(values = treatment_colors, name = "Food web") +
  scale_x_continuous(breaks=c(-1,-0.5,0,0.5,1), labels = c(-1,-0.5,0,0.5,1)) +
  scale_y_log10(breaks=uni_breaks, labels=uni_breaks) + annotation_logticks(sides = "l")+
  coord_cartesian(ylim=c(0.07,1)) 

# Oviposition preference
AF_uni_Pref <- RF_Pref$absolute_fitness %>% 
  select(-sc.Diam, -sc.log.Clutch) %>%
  gather(ID, absolute_fitness, -sc.sqrt.Pref, -Foodweb) %>%
  mutate(ID_group = ifelse(ID == "average", "average", "replicate")) %>% 
  ggplot(., aes(x=sc.sqrt.Pref, y=absolute_fitness, group=interaction(ID,Foodweb), color=Foodweb, alpha=ID_group, size=ID_group)) +
  geom_line() +
  scale_alpha_manual(values=c(1,0.2), guide=FALSE) +
  scale_size_manual(values=c(1.5,0.25), guide=FALSE) +
  xlab("Oviposition preference (SDs)") +
  ylab("Mean larval survival") +
  scale_color_manual(values = treatment_colors, name = "Food web") +
  scale_x_continuous(breaks=c(-1,-0.5,0,0.5,1), labels = c(-1,-0.5,0,0.5,1)) +
  scale_y_log10(breaks=uni_breaks, labels=uni_breaks) + annotation_logticks(sides = "l") +
  coord_cartesian(ylim=c(0.07,1))

# Clutch size 
AF_uni_Clutch <- RF_Clutch$absolute_fitness %>% 
  select(-sc.Diam, -sc.sqrt.Pref) %>%
  gather(ID, absolute_fitness, -sc.log.Clutch, -Foodweb) %>%
  mutate(ID_group = ifelse(ID == "average", "average", "replicate")) %>% 
  ggplot(., aes(x=sc.log.Clutch, y=absolute_fitness, group=interaction(ID,Foodweb), color=Foodweb, alpha=ID_group, size=ID_group)) +
  geom_line() +
  scale_alpha_manual(values=c(1,0.2), guide=FALSE) +
  scale_size_manual(values=c(1.5,0.25), guide=FALSE) +
  xlab("Clutch size (SDs)") +
  ylab("Mean larval survival") +
  scale_color_manual(values = treatment_colors, name = "Food web") +
  scale_x_continuous(breaks=c(-1,-0.5,0,0.5,1), labels = c(-1,-0.5,0,0.5,1)) +
  scale_y_log10(breaks=uni_breaks, labels=uni_breaks) + annotation_logticks(sides = "l") +
  coord_cartesian(ylim=c(0.07,1))
```

Fit plots together to reproduce Figure 2.


```r
# Get legend
AF_uni_legend <- get_legend(AF_uni_Diam + theme(legend.title.align = 0.5)) 

# Make plots
AF_uni_plots <- plot_grid(AF_uni_Diam + theme(legend.position = "none"), 
                       AF_uni_Pref + theme(legend.position = "none", 
                                           axis.title.y = element_blank()) + xlab("Preference (SDs)"),
                      AF_uni_Clutch + theme(legend.position = "none",
                                             axis.title.y = element_blank()), 
                       labels = "AUTO", nrow = 1, align='hv')

AF_gradients <- plot_grid(AF_uni_plots, AF_uni_legend, ncol = 2, rel_widths = c(0.8,0.2)) 
AF_gradients
```

![Adaptive landscapes of gall phenotypes in complex and simple food webs.](manuscript_files/figure-latex/Univariate-Fitness-Landscapes-1.pdf) 

```r
save_plot(filename = "UV_landscapes.pdf", plot = AF_gradients, base_height = 5, base_width = 8.5)
```

\  

## Replot adaptive landscapes of gall traits in terms of relative fitness


```r
RF_plot_foodweb_Clutch_df <- RF_Clutch$absolute_fitness %>% 
  select(-sc.Diam, -sc.sqrt.Pref) %>%
  gather(ID, absolute_fitness, -sc.log.Clutch, -Foodweb) %>%
  mutate(ID_group = ifelse(ID == "average", "average", "replicate")) %>%
  group_by(Foodweb) %>%
  mutate(relative_fitness = absolute_fitness/mean(absolute_fitness)) %>%
  ungroup()

RF_uni_Clutch <- RF_plot_foodweb_Clutch_df %>%
  ggplot(., aes(x=sc.log.Clutch, y=relative_fitness, group=interaction(ID,Foodweb), color=Foodweb, alpha=ID_group, size=ID_group)) +
  geom_line() +
  scale_alpha_manual(values=c(1,0.2), guide=FALSE) +
  scale_size_manual(values=c(1.5,0.25), guide=FALSE) +
  xlab("Clutch size (SDs)") +
  ylab("Relative fitness") +
  scale_color_manual(values = cbPalette[c(6,5)], name = "Food web") + 
  scale_x_continuous(breaks=c(-1,-0.5,0,0.5,1), labels = c(-1,-0.5,0,0.5,1)) +
  coord_cartesian(ylim=c(0.2,2)) +
  geom_hline(yintercept=1, linetype="dotted")

# Oviposition preference

RF_plot_foodweb_Pref_df <- RF_Pref$absolute_fitness %>% 
  select(-sc.Diam, -sc.log.Clutch) %>%
  gather(ID, absolute_fitness, -sc.sqrt.Pref, -Foodweb) %>%
  mutate(ID_group = ifelse(ID == "average", "average", "replicate")) %>%
  group_by(Foodweb) %>%
  mutate(relative_fitness = absolute_fitness/mean(absolute_fitness))%>%
  ungroup()

RF_uni_Pref <- RF_plot_foodweb_Pref_df %>%
  ggplot(., aes(x=sc.sqrt.Pref, y=relative_fitness, group=interaction(ID,Foodweb), color=Foodweb, alpha=ID_group, size=ID_group)) +
  geom_line() +
  scale_alpha_manual(values=c(1,0.2), guide=FALSE) +
  scale_size_manual(values=c(1.5,0.25), guide=FALSE) +
  xlab("Oviposition preference (SDs)") +
  ylab("Relative fitness") +
  scale_color_manual(values = cbPalette[c(6,5)], name = "Food web") + 
  scale_x_continuous(breaks=c(-1,-0.5,0,0.5,1), labels = c(-1,-0.5,0,0.5,1)) +
  coord_cartesian(ylim=c(0.2,2)) +
  theme(legend.title.align = 0.5) +
  geom_hline(yintercept=1, linetype="dotted")

# Chamber diameter
RF_plot_foodweb_Diam_df <- RF_Diam$absolute_fitness %>% 
  select(-sc.sqrt.Pref, -sc.log.Clutch) %>%
  gather(ID, absolute_fitness, -sc.Diam, -Foodweb) %>%
  mutate(ID_group = ifelse(ID == "average", "average", "replicate")) %>%
  group_by(Foodweb) %>%
  mutate(relative_fitness = absolute_fitness/mean(absolute_fitness))%>%
  ungroup()

RF_uni_Diam <- RF_plot_foodweb_Diam_df %>%
  ggplot(., aes(x=sc.Diam, y=relative_fitness, group=interaction(ID,Foodweb), color=Foodweb, alpha=ID_group, size=ID_group)) +
  geom_line() +
  scale_alpha_manual(values=c(1,0.2), guide=FALSE) +
  scale_size_manual(values=c(1.5,0.25), guide=FALSE) +
  xlab("Chamber diameter (SDs)") +
  ylab("Relative fitness") +
  scale_color_manual(values = cbPalette[c(6,5)], name = "Food web") + 
  scale_x_continuous(breaks=c(-1,-0.5,0,0.5,1), labels = c(-1,-0.5,0,0.5,1)) +
  coord_cartesian(ylim=c(0.2,2)) +
  theme(legend.title.align = 0.5) +
  geom_hline(yintercept=1, linetype="dotted")

# Get legend
RF_uni_legend <- get_legend(RF_uni_Diam)  

# Make plots
RF_uni_plots <- plot_grid(RF_uni_Diam + theme(legend.position = "none", 
                                             axis.text.x = element_text(size=10),
                                             axis.title.x = element_text(size=11)),
                     RF_uni_Pref + theme(legend.position = "none",                          
                                             axis.title.y = element_blank(), 
                                             axis.text.x = element_text(size=10),
                                             axis.title.x = element_text(size=11)),
                     RF_uni_Clutch + theme(legend.position = "none",                          
                                             axis.title.y = element_blank(), 
                                             axis.text.x = element_text(size=10),
                                             axis.title.x = element_text(size=11)),
                     ncol = 3, align='hv')

RF_gradients <- plot_grid(RF_uni_plots, RF_uni_legend, ncol = 2, rel_widths = c(0.8,0.2)) 
RF_gradients
```

![Adaptive landscapes of gall phenotypes in terms of relative fitness.](manuscript_files/figure-latex/Plot Univariate RF Landscapes-1.pdf) 

\ 

## Reproduce Figure 3 in manuscript

Use the fitted model to generate predicted estimates of mean fitness for different phenotypic combinations. I restrict these predictions to +/- 1 SD because we can only reliably estimate the shape of the adaptive landscape near the mean phenotype of the population [Arnold et al., 2001](https://link.springer.com/article/10.1023/A:1013373907708).


```r
newdata_Clutch.Pref <- bind_rows(
  expand.grid(Foodweb = "Complex", sc.Diam = 0, sc.log.Clutch = seq(-1,1,length.out=1000), sc.sqrt.Pref = seq(-1,1,length.out=1000)),
  expand.grid(Foodweb = "Simple", sc.Diam = 0, sc.log.Clutch = seq(-1,1,length.out=1000), sc.sqrt.Pref = seq(-1,1,length.out=1000)))

RF_Clutch.Pref <- bootstrap_fitness(
  logistic_model = foodweb_model, 
  newdata = newdata_Clutch.Pref,
  bootstraps=NULL)
```


```r
newdata_Diam.Clutch <- bind_rows(
  expand.grid(Foodweb = "Complex", sc.Diam = seq(-1,1,length.out=1000), sc.log.Clutch = seq(-1,1,length.out=1000), sc.sqrt.Pref = 0),
  expand.grid(Foodweb = "Simple", sc.Diam = seq(-1,1,length.out=1000), sc.log.Clutch = seq(-1,1,length.out=1000), sc.sqrt.Pref=0))

RF_Diam.Clutch <- bootstrap_fitness(
  logistic_model = foodweb_model, 
  newdata = newdata_Diam.Clutch,
  bootstraps=NULL)
```


```r
newdata_Diam.Pref <- bind_rows(
  expand.grid(Foodweb = "Complex", sc.Diam = seq(-1,1,length.out=1000), sc.log.Clutch = 0, sc.sqrt.Pref = seq(-1,1,length.out=1000)),
  expand.grid(Foodweb = "Simple", sc.Diam = seq(-1,1,length.out=1000), sc.log.Clutch = 0, sc.sqrt.Pref = seq(-1,1,length.out=1000)))

RF_Diam.Pref <- bootstrap_fitness(
  logistic_model = foodweb_model, 
  newdata = newdata_Diam.Pref,
  bootstraps=NULL)
```

Create plots for each panel.


```r
# Plot helpers
AF_range <- bind_rows(RF_Clutch.Pref$absolute_fitness, RF_Diam.Clutch$absolute_fitness, RF_Diam.Pref$absolute_fitness) %>%
  summarise(min = min(average), max = max(average))

my_breaks <- c(0.2,0.4,0.6,0.8)

# Set contour breaks
breaks <- predict(foodweb_model, newdata=data.frame(sc.Diam = c(0,0), sc.log.Clutch = c(0,0), sc.sqrt.Pref = c(0,0), Foodweb = c("Complex","Simple")), type="response", re.form=~0)
names(breaks) <- c("Complex","Simple")

# Size x Clutch
AF_plot_Diam.Clutch <- RF_Diam.Clutch$absolute_fitness %>% 
  filter(Foodweb == "Complex") %>%
  ggplot(., aes(x=sc.Diam, y=sc.log.Clutch)) +
  geom_raster(aes(fill=average)) +
  stat_contour(aes(z = average), breaks = breaks["Complex"], linetype = "dotted", color = "black") +
  geom_raster(data = filter(RF_Diam.Clutch$absolute_fitness, Foodweb == "Simple"), aes(fill=average)) +
  stat_contour(data = filter(RF_Diam.Clutch$absolute_fitness, Foodweb == "Simple"), aes(z = average), breaks = breaks["Simple"], linetype = "dotted", color = "black") +
  facet_wrap(~Foodweb) +
  xlab("Gall diameter (SDs)") +
  ylab("Clutch size (SDs)") +
  scale_x_continuous(labels = c(-1,-0.5,0,0.5,1)) +
  scale_y_continuous(labels = c(-1,-0.5,0,0.5,1)) +
  scale_fill_viridis(name = "Mean larval\nsurvival", trans="log", breaks = my_breaks, labels = my_breaks, limits = c(AF_range$min, AF_range$max)) +
  theme(strip.background = element_blank(),
        strip.text.x=element_text(margin=margin(b=4))) 

# Size x Preference
AF_plot_Diam.Pref <- RF_Diam.Pref$absolute_fitness %>% 
  filter(Foodweb == "Complex") %>%
  ggplot(., aes(x=sc.Diam, y=sc.sqrt.Pref)) +
  geom_raster(aes(fill=average)) +
  stat_contour(aes(z = average), breaks = breaks["Complex"], linetype = "dotted", color = "black") +
  geom_raster(data = filter(RF_Diam.Pref$absolute_fitness, Foodweb == "Simple"), aes(fill=average)) +
  stat_contour(data = filter(RF_Diam.Pref$absolute_fitness, Foodweb == "Simple"), aes(z = average), breaks = breaks["Simple"], linetype = "dotted", color = "black") +
  facet_wrap(~Foodweb) +
  xlab("Chamber diameter (SDs)") +
  ylab("Oviposition preference (SDs)") +
  scale_x_continuous(labels = c(-1,-0.5,0,0.5,1)) +
  scale_y_continuous(labels = c(-1,-0.5,0,0.5,1)) +
  scale_fill_viridis(name = "Mean larval\nsurvival", trans="log", breaks = my_breaks, labels = my_breaks, limits = c(AF_range$min, AF_range$max)) +
  theme(strip.background = element_blank(),
        strip.text.x=element_text(margin=margin(b=4))) 

# Clutch x Preference
AF_plot_Clutch.Pref <- RF_Clutch.Pref$absolute_fitness %>%
  filter(Foodweb == "Complex") %>%
  ggplot(., aes(x=sc.sqrt.Pref, y=sc.log.Clutch)) +
  geom_raster(aes(fill=average)) +
  stat_contour(aes(z = average), breaks = breaks["Complex"], linetype = "dotted", color = "black") +
  geom_raster(data = filter(RF_Clutch.Pref$absolute_fitness, Foodweb == "Simple"), aes(fill=average)) +
  stat_contour(data = filter(RF_Clutch.Pref$absolute_fitness, Foodweb == "Simple"), aes(z = average), breaks = breaks["Simple"], linetype = "dotted", color = "black") +
  facet_wrap(~Foodweb) +
  xlab("Oviposition preference (SDs)") +
  ylab("Clutch size (SDs)") +
  scale_x_continuous(labels = c(-1,-0.5,0,0.5,1)) +
  scale_y_continuous(labels = c(-1,-0.5,0,0.5,1)) +
  scale_fill_viridis(name = "Mean larval\nsurvival", trans="log", breaks = my_breaks, labels = my_breaks, limits = c(AF_range$min, AF_range$max)) +
  theme(strip.background = element_blank(),
        strip.text.x=element_text(margin=margin(b=4)))
```

Fit plots together to reproduce Figure 3.


```r
# Get legend
AF_landscape_legend <- get_legend(AF_plot_Diam.Pref + theme(legend.text = element_text(size=10)))

# Arrow data
gradients_Diam.Clutch <- data.frame(Foodweb = c("Complex","Simple"),
                                    sc.Diam = c(0,0), 
                                    sc.log.Clutch = c(0,0), 
                                    xend = c(filter(tidy_foodweb_grads, type=="Complex", term=="sc.Diam")$mean,
                                             filter(tidy_foodweb_grads, type=="Simple", term=="sc.Diam")$mean), 
                                    yend = c(filter(tidy_foodweb_grads, type=="Complex", term=="sc.log.Clutch")$mean, 
                                             filter(tidy_foodweb_grads, type=="Simple", term=="sc.log.Clutch")$mean))

gradients_Clutch.Pref <- data.frame(Foodweb = c("Complex","Simple"),
                                    sc.sqrt.Pref = c(0,0), 
                                    sc.log.Clutch = c(0,0), 
                                    xend = c(filter(tidy_foodweb_grads, type=="Complex", term=="sc.sqrt.Pref")$mean, 
                                             filter(tidy_foodweb_grads, type=="Simple", term=="sc.sqrt.Pref")$mean), 
                                    yend = c(filter(tidy_foodweb_grads, type=="Complex", term=="sc.log.Clutch")$mean, 
                                             filter(tidy_foodweb_grads, type=="Simple", term=="sc.log.Clutch")$mean))

gradients_Diam.Pref <- data.frame(Foodweb = c("Complex","Simple"),
                                    sc.Diam = c(0,0), 
                                    sc.sqrt.Pref = c(0,0), 
                                    xend = c(filter(tidy_foodweb_grads, type=="Complex", term=="sc.Diam")$mean,
                                             filter(tidy_foodweb_grads, type=="Simple", term=="sc.Diam")$mean), 
                                    yend = c(filter(tidy_foodweb_grads, type=="Complex", term=="sc.sqrt.Pref")$mean, 
                                             filter(tidy_foodweb_grads, type=="Simple", term=="sc.sqrt.Pref")$mean))
                                    

# Format plots
AF_landscape_plots <- plot_grid(
  AF_plot_Diam.Clutch + theme(legend.position = "none", 
                              axis.title.x = element_blank(),
                             axis.text.x = element_blank()) +
    geom_segment(data = gradients_Diam.Clutch, 
                 aes(x = sc.Diam, y = sc.log.Clutch, xend = xend, yend = yend), 
                 arrow = arrow(length = unit(0.03, "npc"))),
  AF_plot_Clutch.Pref + theme(legend.position = "none", 
                              axis.title.y = element_blank(),
                              axis.text.y = element_blank()) + xlab("Preference (SDs)") +
    geom_segment(data = gradients_Clutch.Pref, 
                 aes(x = sc.sqrt.Pref, y = sc.log.Clutch, xend = xend, yend = yend), 
                 arrow = arrow(length = unit(0.03, "npc"))), 
  AF_plot_Diam.Pref + theme(legend.position = "none") + ylab("Preference (SDs)") +
    geom_segment(data = gradients_Diam.Pref, 
                 aes(x = sc.Diam, y = sc.sqrt.Pref, xend = xend, yend = yend), 
                 arrow = arrow(length = unit(0.03, "npc"))), 
  nrow=2, align="hv", labels = "AUTO")

AF_landscape_2d <- ggdraw(AF_landscape_plots) + draw_grob(AF_landscape_legend, x = 0.7, y=-0.2)
AF_landscape_2d
```

![Multivariate adaptive landscapes of gall phenotypes in complex and simple food webs. Dotted lines denote fitness for the mean phenotype. Arrows represent the vector of directional selection gradients acting on each phenotype.](manuscript_files/figure-latex/Figure-Multivariate-Landscapes-1.pdf) 

```r
save_plot(filename = "MV_landscapes.pdf", plot = AF_landscape_2d, base_width=8, base_height = 6)
```

\ 

## Slope and curvature of the adaptive landscape of gall traits

Below is the code I used to get estimates for the slope and curvature of the adaptive landscape in each food-web treatment. If 95% confidence intervals of selection gradients overlap zero, then they are set to zero; otherwise, I retain the mean estimate. 


```r
# Create beta-matrix for complex food-web treatment
complex_betas_df <- tidy_foodweb_grads %>%
  filter(type=="Complex", gradient_type=="Directional") %>%
  select(term, gradient=mean, low=conf.low, high=conf.high)
complex_betas_df$term <- factor(complex_betas_df$term, levels=c("sc.Diam","sc.log.Clutch","sc.sqrt.Pref"), labels=c("Diam","Clutch","Pref"), order=TRUE)
rownames(complex_betas_df) <- complex_betas_df$term

complex_betas <- matrix(nrow = 3, ncol = 1, dimnames = list(c("Diam","Clutch","Pref"), c("")))
complex_betas["Diam",] <- ifelse(complex_betas_df["Diam","low"]*complex_betas_df["Diam","high"] > 0, 
                                        complex_betas_df["Diam","gradient"], 0)
complex_betas["Clutch",] <- ifelse(complex_betas_df["Clutch","low"]*complex_betas_df["Clutch","high"] > 0,
                                            complex_betas_df["Clutch","gradient"], 0)
complex_betas["Pref",] <- ifelse(complex_betas_df["Pref","low"]*complex_betas_df["Pref","high"] > 0,
                                        complex_betas_df["Pref","gradient"], 0)
complex_betas <- round(complex_betas, 2)

# Create gamma-matrix for complex food-web treatment
complex_gammas_df <- tidy_foodweb_grads %>%
  filter(type=="Complex", gradient_type %in% c("Quadratic","Correlational")) %>%
  select(term, gradient=mean, low=conf.low, high=conf.high)
complex_gammas_df$term <- factor(complex_gammas_df$term, levels=c("I(sc.Diam^2)","I(sc.log.Clutch^2)","I(sc.sqrt.Pref^2)","sc.Diam:sc.log.Clutch","sc.Diam:sc.sqrt.Pref","sc.log.Clutch:sc.sqrt.Pref"), labels=c("Diam^2","Clutch^2","Pref^2","Diam:Clutch","Diam:Pref","Clutch:Pref"), order=TRUE)
rownames(complex_gammas_df) <- complex_gammas_df$term
  
complex_gammas <- matrix(nrow = 3, ncol = 3, dimnames = list(c("Diam","Clutch","Pref"), c("Diam","Clutch","Pref")))
complex_gammas["Diam","Diam"] <- ifelse(complex_gammas_df["Diam^2","low"]*complex_gammas_df["Diam^2","high"] > 0, 
                                        complex_gammas_df["Diam^2","gradient"], 0)
complex_gammas["Clutch","Clutch"] <- ifelse(complex_gammas_df["Clutch^2","low"]*complex_gammas_df["Clutch^2","high"] > 0,
                                            complex_gammas_df["Clutch^2","gradient"], 0)
complex_gammas["Pref","Pref"] <- ifelse(complex_gammas_df["Pref^2","low"]*complex_gammas_df["Pref^2","high"] > 0,
                                        complex_gammas_df["Pref^2","gradient"], 0)
complex_gammas["Diam","Clutch"] <- ifelse(complex_gammas_df["Diam:Clutch","low"]*complex_gammas_df["Diam:Clutch","high"] > 0,
                                          complex_gammas_df["Diam:Clutch","gradient"], 0)
complex_gammas["Clutch","Diam"] <- complex_gammas["Diam","Clutch"]
complex_gammas["Diam","Pref"] <- ifelse(complex_gammas_df["Diam:Pref","low"]*complex_gammas_df["Diam:Pref","high"] > 0,
                                            complex_gammas_df["Diam:Pref","gradient"], 0)
complex_gammas["Pref","Diam"] <- complex_gammas["Diam","Pref"]
complex_gammas["Clutch","Pref"] <- ifelse(complex_gammas_df["Clutch:Pref","low"]*complex_gammas_df["Clutch:Pref","high"] > 0, 
                                        complex_gammas_df["Clutch:Pref","gradient"], 0)
complex_gammas["Pref","Clutch"] <- complex_gammas["Clutch","Pref"]
complex_gammas <- round(complex_gammas, 2)

# Create matrix of Betas in simple food-web treatment
simple_betas_df <- tidy_foodweb_grads %>%
  filter(type=="Simple", gradient_type=="Directional") %>%
  select(term, gradient=mean, low=conf.low, high=conf.high)
simple_betas_df$term <- factor(simple_betas_df$term, levels=c("sc.Diam","sc.log.Clutch","sc.sqrt.Pref"), labels=c("Diam","Clutch","Pref"), order=TRUE)
rownames(simple_betas_df) <- simple_betas_df$term

simple_betas <- matrix(nrow = 3, ncol = 1, dimnames = list(c("Diam","Clutch","Pref"), c("")))
simple_betas["Diam",] <- ifelse(simple_betas_df["Diam","low"]*simple_betas_df["Diam","high"] > 0, 
                                        simple_betas_df["Diam","gradient"], 0)
simple_betas["Clutch",] <- ifelse(simple_betas_df["Clutch","low"]*simple_betas_df["Clutch","high"] > 0,
                                            simple_betas_df["Clutch","gradient"], 0)
simple_betas["Pref",] <- ifelse(simple_betas_df["Pref","low"]*simple_betas_df["Pref","high"] > 0,
                                        simple_betas_df["Pref","gradient"], 0)
simple_betas <- round(simple_betas, 2)

# Create gamma-matrix for simple food-web treatment
simple_gammas_df <- tidy_foodweb_grads %>%
  filter(type=="Simple", gradient_type %in% c("Quadratic","Correlational")) %>%
  select(term, gradient=mean, low=conf.low, high=conf.high)
simple_gammas_df$term <- factor(simple_gammas_df$term, levels=c("I(sc.Diam^2)","I(sc.log.Clutch^2)","I(sc.sqrt.Pref^2)","sc.Diam:sc.log.Clutch","sc.Diam:sc.sqrt.Pref","sc.log.Clutch:sc.sqrt.Pref"), labels=c("Diam^2","Clutch^2","Pref^2","Diam:Clutch","Diam:Pref","Clutch:Pref"), order=TRUE)
rownames(simple_gammas_df) <- simple_gammas_df$term

simple_gammas <- matrix(nrow = 3, ncol = 3, dimnames = list(c("Diam","Clutch","Pref"), c("Diam","Clutch","Pref")))
simple_gammas["Diam","Diam"] <- ifelse(simple_gammas_df["Diam^2","low"]*simple_gammas_df["Diam^2","high"] > 0, 
                                        simple_gammas_df["Diam^2","gradient"], 0)
simple_gammas["Clutch","Clutch"] <- ifelse(simple_gammas_df["Clutch^2","low"]*simple_gammas_df["Clutch^2","high"] > 0,
                                            simple_gammas_df["Clutch^2","gradient"], 0)
simple_gammas["Pref","Pref"] <- ifelse(simple_gammas_df["Pref^2","low"]*simple_gammas_df["Pref^2","high"] > 0,
                                        simple_gammas_df["Pref^2","gradient"], 0)
simple_gammas["Diam","Clutch"] <- ifelse(simple_gammas_df["Diam:Clutch","low"]*simple_gammas_df["Diam:Clutch","high"] > 0,
                                          simple_gammas_df["Diam:Clutch","gradient"], 0)
simple_gammas["Clutch","Diam"] <- simple_gammas["Diam","Clutch"]
simple_gammas["Diam","Pref"] <- ifelse(simple_gammas_df["Diam:Pref","low"]*simple_gammas_df["Diam:Pref","high"] > 0,
                                            simple_gammas_df["Diam:Pref","gradient"], 0)
simple_gammas["Pref","Diam"] <- simple_gammas["Diam","Pref"]
simple_gammas["Clutch","Pref"] <- ifelse(simple_gammas_df["Clutch:Pref","low"]*simple_gammas_df["Clutch:Pref","high"] > 0, 
                                        simple_gammas_df["Clutch:Pref","gradient"], 0)
simple_gammas["Pref","Clutch"] <- simple_gammas["Clutch","Pref"]
simple_gammas <- round(simple_gammas, 2)

# Calculate Curvature of the Fitness landscape
complex_curvature <- complex_gammas - complex_betas %*% t(complex_betas)
simple_curvature <- simple_gammas - simple_betas %*% t(simple_betas)
```

I use these estimates to populate the curvature matrices presented in the **Results** of the main text. I reproduce these estimates here for clarity:

$$\textbf{C} = \begin{pmatrix} C_{\text{Diam:Diam}}&& \\ C_{\text{Clutch:Diam}} & C_{\text{Clutch:Clutch}} & \\ C_{\text{Pref:Diam}} & C_{\text{Pref:Clutch}} & C_{\text{Pref:Pref}} \end{pmatrix}$$

$$\textbf{C}_{\text{Complex}} = \begin{pmatrix} 
-0.12 &  &  \\  
0 & 0 &  \\  
0 & 0 & 0.34 \end{pmatrix}$$

$$\textbf{C}_{\text{Simple}} = \begin{pmatrix} 
-0.04 &  &  \\  
0.02 & -0.01 &  \\  
0.03 & -0.01 & -0.03 \end{pmatrix}$$

\  

## Trait-fitness relationships of the extended phenotype of the egg parasitoid *Platygaster* sp.


```r
# convert "gall_survival" to egg parasitoid survival. Note that both Iteomyia pupa and larva parasitoids result in 0. We do not change the name of the response variable to make clear that we are using the same statistical models as for Iteomyia.
eggegg_df <- mutate(gall_selection.df, gall_survival = ifelse(egg.ptoid==1,1,0))

# fit trait-fitness relationship using the same model structure
eggegg_model <- update(foodweb_model, data=eggegg_df)

# subset data for quantifying biased selection on chamber diameter
biased_eggegg_df <- eggegg_df %>%
  group_by(Foodweb, Gall_Number) %>%
  mutate(mean_survival = mean(gall_survival)) %>%
  filter(mean_survival > 0, mean_survival < 1) %>%
  ungroup()

# quantify bias
biased_eggegg_model <- update(biased_foodweb_model, data=biased_eggegg_df)
```

Use parametric bootstrapping to calculating 95% confidence intervals of trait-fitness relationships.



Remove all higher-order terms to quantify linear trait-fitness relationships.


```r
linear_eggegg_model <- update(linear_foodweb_model, data=eggegg_df)
```

And refit using parametric bootstrapping.


```r
boot_linear_eggegg_model <- bootMer(linear_eggegg_model, FUN = fixef, nsim=n_boots_analysis, parallel="multicore", ncpus=32, seed=34)
```

Gather regression coefficients to display trait-fitness relationships.


```r
eggegg_bind <- bind_cols(
  as_tibble(boot_linear_eggegg_model$t),
  select(as_tibble(boot_eggegg_model$t), `FoodwebComplex:I(sc.Diam^2)`:`FoodwebSimple:sc.log.Clutch:sc.sqrt.Pref`) # only retain nonlinear coefficients from full model
)

# Adjust coefficients with chamber diameter
eggegg_adj_Diam <- as_tibble(boot_eggegg_model$t) %>%
  mutate(`FoodwebComplex:sc.Diam` = `FoodwebComplex:sc.Diam` - fixef(biased_eggegg_model)["FoodwebComplex:sc.Diam"],
         `FoodwebSimple:sc.Diam` = `FoodwebSimple:sc.Diam` - fixef(biased_eggegg_model)["FoodwebSimple:sc.Diam"])

# Calculate differences between food-web treatments, then estimate means and confidence intervals
eggegg_alphas <- eggegg_adj_Diam %>%
  mutate(`FoodwebSimple:(Intercept)`=inverse_logit(FoodwebSimple), `FoodwebComplex:(Intercept)`=inverse_logit(FoodwebComplex)) %>% # logit transform to make intercept interpretable as a probability
  select(-FoodwebSimple, -FoodwebComplex) %>% # translated into (Intercept) terms
  mutate(`FoodwebDiff:(Intercept)` = `FoodwebSimple:(Intercept)` - `FoodwebComplex:(Intercept)`,
         `FoodwebDiff:sc.Diam` = `FoodwebSimple:sc.Diam` - `FoodwebComplex:sc.Diam`,
         `FoodwebDiff:sc.log.Clutch` = `FoodwebSimple:sc.log.Clutch` - `FoodwebComplex:sc.log.Clutch`,
         `FoodwebDiff:sc.sqrt.Pref` = `FoodwebSimple:sc.sqrt.Pref` - `FoodwebComplex:sc.sqrt.Pref`,
         `FoodwebDiff:I(sc.Diam^2)` = `FoodwebSimple:I(sc.Diam^2)` - `FoodwebComplex:I(sc.Diam^2)`,
         `FoodwebDiff:I(sc.log.Clutch^2)` = `FoodwebSimple:I(sc.log.Clutch^2)` - `FoodwebComplex:I(sc.log.Clutch^2)`,
         `FoodwebDiff:I(sc.sqrt.Pref^2)` = `FoodwebSimple:I(sc.sqrt.Pref^2)` - `FoodwebComplex:I(sc.sqrt.Pref^2)`,
         `FoodwebDiff:sc.Diam:sc.log.Clutch` = `FoodwebSimple:sc.Diam:sc.log.Clutch` - `FoodwebComplex:sc.Diam:sc.log.Clutch`,
         `FoodwebDiff:sc.Diam:sc.sqrt.Pref` = `FoodwebSimple:sc.Diam:sc.sqrt.Pref` - `FoodwebComplex:sc.Diam:sc.sqrt.Pref`,
         `FoodwebDiff:sc.log.Clutch:sc.sqrt.Pref` = `FoodwebSimple:sc.log.Clutch:sc.sqrt.Pref` - `FoodwebComplex:sc.log.Clutch:sc.sqrt.Pref`) %>%
  summarise_all(funs(mean, conf.high, conf.low))

# Tidy up estimates
tidy_eggegg_alphas <- eggegg_alphas %>%
  t() %>%
  as.data.frame() %>%
  transmute(rows = rownames(.), value=V1) %>%
  separate(col = rows, into = c("term","estimate"), sep="_") %>%
  separate(col = term, into = c("type","term_1","term_2"), sep=":", extra="drop") %>%
  mutate(term = ifelse(is.na(term_2), term_1, paste(term_1,term_2,sep=":"))) %>%
  separate(col = type, into = c("extra","type"), sep = 7) %>%
  select(type, term, estimate, value) %>%
  spread(estimate, value) %>%
  mutate(P_cutoff = ifelse(conf.high*conf.low < 0, "","*"),
         coefficient_type = ifelse(term=="(Intercept)","Mean fitness",
                                   ifelse(grepl("\\^2", .$term)==TRUE, "Quadratic",
                                          ifelse(grepl(":", .$term)==TRUE, "Correlational","Linear"))))
```

\  

\begin{table}[t]

\caption{\label{tab:Egg ptoid coefficients in each food-web treatment}Table of trait-fitness relationships of extended phenotype of egg parasitoids in complex and simple food webs.}
\centering
\begin{tabular}{l|l|l|r|r|r|l}
\hline
term & coefficient\_type & type & mean & conf.low & conf.high & P\_cutoff\\
\hline
(Intercept) & Mean fitness & Complex & 0.162 & 0.022 & 0.313 & *\\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & Complex & -0.260 & -0.690 & 0.126 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & Complex & -0.025 & -0.503 & 0.440 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & Complex & -0.632 & -1.346 & -0.170 & *\\
\hline
sc.Diam & Linear & Complex & -1.189 & -1.895 & -0.692 & *\\
\hline
sc.Diam:sc.log.Clutch & Correlational & Complex & 0.275 & -0.248 & 0.878 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & Complex & 0.640 & 0.025 & 1.385 & *\\
\hline
sc.log.Clutch & Linear & Complex & 0.823 & 0.173 & 1.554 & *\\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & Complex & -0.428 & -1.213 & 0.248 & \\
\hline
sc.sqrt.Pref & Linear & Complex & 0.287 & -0.414 & 1.027 & \\
\hline
(Intercept) & Mean fitness & Simple & 0.258 & 0.075 & 0.471 & *\\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & Simple & -0.274 & -0.724 & 0.150 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & Simple & 0.348 & -0.170 & 0.928 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & Simple & -0.059 & -0.520 & 0.496 & \\
\hline
sc.Diam & Linear & Simple & -1.573 & -2.747 & -0.865 & *\\
\hline
sc.Diam:sc.log.Clutch & Correlational & Simple & 0.384 & -0.183 & 1.108 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & Simple & -0.010 & -0.645 & 0.524 & \\
\hline
sc.log.Clutch & Linear & Simple & 0.915 & 0.233 & 1.846 & *\\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & Simple & 0.060 & -0.443 & 0.658 & \\
\hline
sc.sqrt.Pref & Linear & Simple & 1.314 & 0.666 & 2.361 & *\\
\hline
\end{tabular}
\end{table}

\  

\begin{table}[t]

\caption{\label{tab:Table of differences in egg parasitoid coefficients between food-web treatments}Table of differences in trait-fitness relationships of egg parasitoids between food-web treatment.}
\centering
\begin{tabular}{l|l|r|r|r|l}
\hline
term & coefficient\_type & mean & conf.low & conf.high & P\_cutoff\\
\hline
(Intercept) & Mean fitness & 0.096 & -0.086 & 0.317 & \\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & -0.014 & -0.601 & 0.539 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & 0.373 & -0.362 & 1.106 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & 0.572 & -0.080 & 1.424 & \\
\hline
sc.Diam & Linear & -0.384 & -1.312 & 0.416 & \\
\hline
sc.Diam:sc.log.Clutch & Correlational & 0.110 & -0.705 & 1.003 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & -0.650 & -1.676 & 0.204 & \\
\hline
sc.log.Clutch & Linear & 0.092 & -0.825 & 1.111 & \\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & 0.487 & -0.386 & 1.452 & \\
\hline
sc.sqrt.Pref & Linear & 1.027 & 0.041 & 2.310 & *\\
\hline
\end{tabular}
\end{table}

\ 

Visualize tables in figure form.


```r
# Plot helpers
dodge_width <- position_dodge(width=0.5)
y_limits_eggegg_coefs <- c(min(filter(tidy_eggegg_alphas, type!="Diff")$conf.low), 
                            max(filter(tidy_eggegg_alphas, type!="Diff")$conf.high))
P_asterisk_eggegg_coefs <- 1.7
P_size <- 8
legend_eggegg_labels <- c("Complex","Simple")
point.size <- 3

# Plot mean fitness
plot_eggegg_mean_fitness <- tidy_eggegg_alphas %>%
  filter(type!="Diff", term=="(Intercept)") %>%
  ggplot(., aes(x=type)) +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_eggegg_alphas, type=="Diff", term=="(Intercept)"), aes(y=0.9, x=1.5, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Probability of observing egg parasitoid", limits=c(0,1), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
  scale_x_discrete(name="") +
  facet_wrap(~coefficient_type) +
  scale_fill_manual(values=cbPalette[c(2,5)], labels=legend_eggegg_labels) +
  scale_color_manual(values=cbPalette[c(2,5)], labels=legend_eggegg_labels) 

# Plot linear trait-fitness relationaships
plot_eggegg_linear_coefs <- tidy_eggegg_alphas %>%
  filter(type!="Diff", term %in% c("sc.Diam","sc.log.Clutch","sc.sqrt.Pref")) %>%
  ggplot(., aes(x=term)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_eggegg_alphas, type=="Diff", term %in% c("sc.Diam","sc.log.Clutch","sc.sqrt.Pref")), 
            aes(y=P_asterisk_eggegg_coefs, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Coefficient", limits=y_limits_eggegg_coefs) +
  scale_x_discrete(name="", labels=parse(text=c("alpha[Diam]","alpha[Clutch]","alpha[Pref]"))) +
  facet_wrap(~coefficient_type) +
  scale_fill_manual(values=cbPalette[c(2,5)], labels=legend_eggegg_labels) +
  scale_color_manual(values=cbPalette[c(2,5)], labels=legend_eggegg_labels)

# Plot quadratic trait-fitness relationships
plot_eggegg_quadratic_coefs <- tidy_eggegg_alphas %>%
  filter(type!="Diff", term %in% c("I(sc.Diam^2)","I(sc.log.Clutch^2)","I(sc.sqrt.Pref^2)")) %>%
  ggplot(., aes(x=term, color=type)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_eggegg_alphas, type=="Diff", term %in% c("I(sc.Diam^2)","I(sc.log.Clutch^2)","I(sc.sqrt.Pref^2)")), 
            aes(y=P_asterisk_eggegg_coefs, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Coefficient", limits=y_limits_eggegg_coefs) +
  scale_x_discrete(name="", labels=parse(text=c("alpha[Diam:Diam]","alpha[Clutch:Clutch]","alpha[Pref:Pref]"))) +
  facet_wrap(~coefficient_type) +
  scale_fill_manual(values=cbPalette[c(2,5)], labels=legend_eggegg_labels) +
  scale_color_manual(values=cbPalette[c(2,5)], labels=legend_eggegg_labels)

# Plot correlational trait-fitness relationaships
plot_eggegg_correlational_coefs <- tidy_eggegg_alphas %>%
  filter(type!="Diff", term %in% c("sc.Diam:sc.log.Clutch","sc.Diam:sc.sqrt.Pref","sc.log.Clutch:sc.sqrt.Pref")) %>%
  ggplot(., aes(x=term, color=type)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_eggegg_alphas, type=="Diff", term %in% c("sc.Diam:sc.log.Clutch","sc.Diam:sc.sqrt.Pref","sc.log.Clutch:sc.sqrt.Pref")), 
            aes(y=P_asterisk_eggegg_coefs, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Coefficient", limits=y_limits_eggegg_coefs) +
  scale_x_discrete(name="", labels=parse(text=c("alpha[Diam:Clutch]","alpha[Diam:Pref]","alpha[Clutch:Pref]"))) +
  facet_wrap(~coefficient_type) +
  scale_fill_manual(values=cbPalette[c(2,5)], labels=legend_eggegg_labels) +
  scale_color_manual(values=cbPalette[c(2,5)], labels=legend_eggegg_labels) 

# Format plots for manuscript
plot_eggegg_coefs <- plot_grid(
  plot_eggegg_mean_fitness + theme(legend.position = "none"),
  plot_eggegg_linear_coefs + theme(legend.position = "none"),
  plot_eggegg_quadratic_coefs + theme(legend.position = "none"),
  plot_eggegg_correlational_coefs + theme(legend.position = "none"), 
  nrow=2, align="hv", labels = "")

figure_eggegg_coefs <- plot_eggegg_coefs 
figure_eggegg_coefs
```

![Trait-fitness relationships of egg parasitoid's extended phenotype in complex and simple food webs. Asterisks denote significant differences (*P* < 0.05) between food-web treatments.](manuscript_files/figure-latex/Plot egg vs. egg ptoid coefficients-1.pdf) 

\  

## Selection gradients acting on extended phenotype of the egg parasitoid


```r
# Estimate mean fitness and mean "brackets" for each food-web treatment (see Janzen and Stern 1998 equation 4 for details about "brackets")
eggegg_complex_predict <- predict(eggegg_model, newdata=filter(eggegg_df, Foodweb=="Complex"), type="response")
eggegg_complex_mean_brackets <- mean(eggegg_complex_predict * (1 - eggegg_complex_predict))
eggegg_complex_mean_fitness <- mean(eggegg_complex_predict)

eggegg_simple_predict <- predict(eggegg_model, newdata=filter(eggegg_df, Foodweb=="Simple"), type="response")
eggegg_simple_mean_brackets <- mean(eggegg_simple_predict * (1 - eggegg_simple_predict))
eggegg_simple_mean_fitness <- mean(eggegg_simple_predict)

eggegg_complex_gradient <- function(x) eggegg_complex_mean_brackets * x / eggegg_complex_mean_fitness
eggegg_simple_gradient <- function(x) eggegg_simple_mean_brackets * x / eggegg_simple_mean_fitness

eggegg_complex_raw_gradients <- select(eggegg_adj_Diam, starts_with("FoodwebComplex:")) %>%
  transmute_all(funs(eggegg_complex_gradient))

eggegg_simple_raw_gradients <- select(eggegg_adj_Diam, starts_with("FoodwebSimple:")) %>%
  transmute_all(funs(eggegg_simple_gradient))

# Calculate differences between food-web treatments, then estimate means and confidence intervals
eggegg_grads <- bind_cols(eggegg_complex_raw_gradients, eggegg_simple_raw_gradients) %>%
  mutate(`FoodwebDiff:sc.Diam` = `FoodwebSimple:sc.Diam` - `FoodwebComplex:sc.Diam`,
         `FoodwebDiff:sc.log.Clutch` = `FoodwebSimple:sc.log.Clutch` - `FoodwebComplex:sc.log.Clutch`,
         `FoodwebDiff:sc.sqrt.Pref` = `FoodwebSimple:sc.sqrt.Pref` - `FoodwebComplex:sc.sqrt.Pref`,
         `FoodwebDiff:I(sc.Diam^2)` = `FoodwebSimple:I(sc.Diam^2)` - `FoodwebComplex:I(sc.Diam^2)`,
         `FoodwebDiff:I(sc.log.Clutch^2)` = `FoodwebSimple:I(sc.log.Clutch^2)` - `FoodwebComplex:I(sc.log.Clutch^2)`,
         `FoodwebDiff:I(sc.sqrt.Pref^2)` = `FoodwebSimple:I(sc.sqrt.Pref^2)` - `FoodwebComplex:I(sc.sqrt.Pref^2)`,
         `FoodwebDiff:sc.Diam:sc.log.Clutch` = `FoodwebSimple:sc.Diam:sc.log.Clutch` - `FoodwebComplex:sc.Diam:sc.log.Clutch`,
         `FoodwebDiff:sc.Diam:sc.sqrt.Pref` = `FoodwebSimple:sc.Diam:sc.sqrt.Pref` - `FoodwebComplex:sc.Diam:sc.sqrt.Pref`,
         `FoodwebDiff:sc.log.Clutch:sc.sqrt.Pref` = `FoodwebSimple:sc.log.Clutch:sc.sqrt.Pref` - `FoodwebComplex:sc.log.Clutch:sc.sqrt.Pref`) %>%
  summarise_all(funs(mean, conf.high, conf.low))

# Tidy up estimates
tidy_eggegg_grads <- eggegg_grads %>%
  t() %>%
  as.data.frame() %>%
  transmute(rows = rownames(.), value=V1) %>%
  separate(col = rows, into = c("term","estimate"), sep="_") %>%
  separate(col = term, into = c("type","term_1","term_2"), sep=":", extra="drop") %>%
  mutate(term = ifelse(is.na(term_2), term_1, paste(term_1,term_2,sep=":"))) %>%
  separate(col = type, into = c("extra","type"), sep = 7) %>%
  select(type, term, estimate, value) %>%
  spread(estimate, value) %>%
  mutate(P_cutoff = ifelse(conf.high*conf.low < 0, "","*"),
         gradient_type = ifelse(grepl("\\^2", .$term)==TRUE, "Quadratic",
                                ifelse(grepl(":", .$term)==TRUE, "Correlational","Directional")),
         mean = ifelse(gradient_type=="Quadratic", 2*mean, mean),
         conf.low = ifelse(gradient_type=="Quadratic", 2*conf.low, conf.low),
         conf.high = ifelse(gradient_type=="Quadratic", 2*conf.high, conf.high))
```

\  

\begin{table}[t]

\caption{\label{tab:Egg vs. ptoid selection gradients}Table of selection gradients acting on extended phenotype of egg parasitoids in complex and simple food webs.}
\centering
\begin{tabular}{l|l|l|r|r|r|l}
\hline
term & gradient\_type & type & mean & conf.low & conf.high & P\_cutoff\\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & Complex & -0.201 & -0.534 & 0.098 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & Complex & -0.019 & -0.389 & 0.340 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & Complex & -0.489 & -1.041 & -0.132 & *\\
\hline
sc.Diam & Directional & Complex & -0.460 & -0.733 & -0.268 & *\\
\hline
sc.Diam:sc.log.Clutch & Correlational & Complex & 0.106 & -0.096 & 0.340 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & Complex & 0.247 & 0.010 & 0.536 & *\\
\hline
sc.log.Clutch & Directional & Complex & 0.319 & 0.067 & 0.601 & *\\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & Complex & -0.165 & -0.469 & 0.096 & \\
\hline
sc.sqrt.Pref & Directional & Complex & 0.111 & -0.160 & 0.398 & \\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & Simple & -0.150 & -0.396 & 0.082 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & Simple & 0.190 & -0.093 & 0.508 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & Simple & -0.033 & -0.284 & 0.272 & \\
\hline
sc.Diam & Directional & Simple & -0.430 & -0.752 & -0.237 & *\\
\hline
sc.Diam:sc.log.Clutch & Correlational & Simple & 0.105 & -0.050 & 0.303 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & Simple & -0.003 & -0.176 & 0.143 & \\
\hline
sc.log.Clutch & Directional & Simple & 0.250 & 0.064 & 0.505 & *\\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & Simple & 0.016 & -0.121 & 0.180 & \\
\hline
sc.sqrt.Pref & Directional & Simple & 0.359 & 0.182 & 0.646 & *\\
\hline
\end{tabular}
\end{table}

\  

\begin{table}[t]

\caption{\label{tab:Table of differences in selection gradients acting on egg parasitoids between food-web treatment}Table of differences in selection gradients acting on egg parasitoids between food-web treatment.}
\centering
\begin{tabular}{l|l|r|r|r|l}
\hline
term & gradient\_type & mean & conf.low & conf.high & P\_cutoff\\
\hline
I(sc.Diam\textasciicircum{}2) & Quadratic & 0.051 & -0.335 & 0.433 & \\
\hline
I(sc.log.Clutch\textasciicircum{}2) & Quadratic & 0.210 & -0.260 & 0.678 & \\
\hline
I(sc.sqrt.Pref\textasciicircum{}2) & Quadratic & 0.456 & 0.015 & 1.065 & *\\
\hline
sc.Diam & Directional & 0.030 & -0.247 & 0.304 & \\
\hline
sc.Diam:sc.log.Clutch & Correlational & -0.001 & -0.274 & 0.289 & \\
\hline
sc.Diam:sc.sqrt.Pref & Correlational & -0.250 & -0.595 & 0.045 & \\
\hline
sc.log.Clutch & Directional & -0.068 & -0.382 & 0.257 & \\
\hline
sc.log.Clutch:sc.sqrt.Pref & Correlational & 0.182 & -0.123 & 0.524 & \\
\hline
sc.sqrt.Pref & Directional & 0.248 & -0.086 & 0.641 & \\
\hline
\end{tabular}
\end{table}

\  

Now let's see the table results as a figure.


```r
# Plot helpers
dodge_width <- position_dodge(width=0.5)
y_limits_eggegg_grads <- c(min(filter(tidy_eggegg_grads, type!="Diff")$conf.low), 
                            max(filter(tidy_eggegg_grads, type!="Diff")$conf.high))
P_asterisk_eggegg_grads <- 0.45
P_size <- 8
legend_eggegg_labels <- c("Complex","Simple")
legend_eggegg_title <- "Food web"
point.size <- 3

# Plot directional selection gradients
plot_eggegg_directional_grads <- tidy_eggegg_grads %>%
  filter(type!="Diff", term %in% c("sc.Diam","sc.log.Clutch","sc.sqrt.Pref")) %>%
  ggplot(., aes(x=term)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_eggegg_grads, type=="Diff", term %in% c("sc.Diam","sc.log.Clutch","sc.sqrt.Pref")), 
            aes(y=P_asterisk_eggegg_grads, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Gradient", limits=y_limits_eggegg_grads) +
  scale_x_discrete(name="", labels=parse(text=c("beta[Diam]","beta[Clutch]","beta[Pref]"))) +
  facet_wrap(~gradient_type) +
  scale_fill_manual(values=cbPalette[c(2,5)], labels=legend_eggegg_labels, name=legend_eggegg_title) +
  scale_color_manual(values=cbPalette[c(2,5)], labels=legend_eggegg_labels, name=legend_eggegg_title)

# Plot quadratic selection gradients
plot_eggegg_quadratic_grads <- tidy_eggegg_grads %>%
  filter(type!="Diff", term %in% c("I(sc.Diam^2)","I(sc.log.Clutch^2)","I(sc.sqrt.Pref^2)")) %>%
  ggplot(., aes(x=term, color=type)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_eggegg_grads, type=="Diff", term %in% c("I(sc.Diam^2)","I(sc.log.Clutch^2)","I(sc.sqrt.Pref^2)")), 
            aes(y=P_asterisk_eggegg_grads, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Gradient", limits=y_limits_eggegg_grads) +
  scale_x_discrete(name="", labels=parse(text=c("gamma[Diam:Diam]","gamma[Clutch:Clutch]","gamma[Pref:Pref]"))) +
  facet_wrap(~gradient_type) +
  scale_fill_manual(values=cbPalette[c(2,5)], labels=legend_eggegg_labels, name=legend_eggegg_title) +
  scale_color_manual(values=cbPalette[c(2,5)], labels=legend_eggegg_labels, name=legend_eggegg_title)

# Plot correlational selection gradients
plot_eggegg_correlational_grads <- tidy_eggegg_grads %>%
  filter(type!="Diff", term %in% c("sc.Diam:sc.log.Clutch","sc.Diam:sc.sqrt.Pref","sc.log.Clutch:sc.sqrt.Pref")) %>%
  ggplot(., aes(x=term, color=type)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_linerange(aes(ymax=conf.high, ymin=conf.low, color=type), position=dodge_width) +
  geom_point(aes(y=mean, fill=type), position=dodge_width, shape=21, color="black", size=point.size) +
  geom_text(data=filter(tidy_eggegg_grads, type=="Diff", term %in% c("sc.Diam:sc.log.Clutch","sc.Diam:sc.sqrt.Pref","sc.log.Clutch:sc.sqrt.Pref")), 
            aes(y=P_asterisk_eggegg_grads, label=P_cutoff), size=P_size, color="black") +
  scale_y_continuous(name="Gradient", limits=y_limits_eggegg_grads) +
  scale_x_discrete(name="", labels=parse(text=c("gamma[Diam:Clutch]","gamma[Diam:Pref]","gamma[Clutch:Pref]"))) +
  facet_wrap(~gradient_type) +
  scale_fill_manual(values=cbPalette[c(2,5)], labels=legend_eggegg_labels, name=legend_eggegg_title) +
  scale_color_manual(values=cbPalette[c(2,5)], labels=legend_eggegg_labels, name=legend_eggegg_title) 

# Get legend
legend_eggegg_grads <- get_legend(plot_eggegg_directional_grads)

# Format plots for manuscript
plot_eggegg_grads <- plot_grid(
  plot_eggegg_directional_grads + theme(legend.position = "none", axis.title = element_blank()),
  plot_eggegg_quadratic_grads + theme(legend.position = "none", axis.title = element_blank()),
  plot_eggegg_correlational_grads + theme(legend.position = "none", axis.title = element_blank()), 
  nrow=2, align="hv", labels = "")
plot_legend_eggegg_grads <- ggdraw(plot_eggegg_grads) + draw_grob(legend_eggegg_grads, x = 0.7, y=-0.2)
y_title_eggegg_grads <- ggdraw() + draw_label("Selection gradient (SDs)", angle = 90)

figure_eggegg_grads <- plot_grid(y_title_eggegg_grads, plot_legend_eggegg_grads, ncol=2, rel_widths = c(0.05,1)) 
figure_eggegg_grads
```

![Selection gradients acting on egg parasitoid's extended phenotype in complex and simple food webs. Asterisks denote significant differences (*P* < 0.05) between food-web treatments.](manuscript_files/figure-latex/Plot egg vs. egg ptoid gradients-1.pdf) 

\ 

## Reproduce Figure 4 in main text


```r
# this is used for plotting selection on egg parasitoids extended phenotype. We only did this for preference since this was the only trait under selection.
eggegg_RF_Pref <- bootstrap_fitness(
  logistic_model = eggegg_model, 
  newdata = newdata_Pref,
  bootstraps=n_boots_plots) 

eggptoid_plot_df <- eggegg_RF_Pref$absolute_fitness %>% 
  select(-sc.Diam, -sc.log.Clutch) %>%
  gather(ID, absolute_fitness, -sc.sqrt.Pref, -Foodweb) %>%
  mutate(ID_group = ifelse(ID == "average", "average", "replicate")) %>%
  spread(Foodweb, absolute_fitness) %>%
  mutate(fitness_difference = Complex-Simple)

eggptoid_selection_plot <- eggptoid_plot_df %>%
  ggplot(., aes(x=sc.sqrt.Pref, y=fitness_difference, group=ID, alpha=ID_group, size=ID_group)) +
  geom_line(color=cbPalette[5]) +
  scale_alpha_manual(values=c(1,0.2), guide=FALSE) +
  scale_size_manual(values=c(1.5,0.25), guide=FALSE) +
  xlab("Gall midge preference (SDs)") + 
  ylab(expression(paste(Delta," egg parasitoid survival (complex - simple)"))) +
  scale_x_continuous(breaks=c(-1,-0.5,0,0.5,1), labels = c(-1,-0.5,0,0.5,1)) +
  geom_hline(yintercept=0, linetype="dotted") +
  theme_cowplot(font_size = 12)

eggptoid_selection_plot
```

![Selection on extended phenotype (oviposition preference) of egg parasitoids](manuscript_files/figure-latex/Egg-Parasitoid-Selection-1.pdf) 

```r
save_plot(filename = "selection_on_Platygaster.pdf", plot = eggptoid_selection_plot, base_height = 4) # base_height = 6 is to large for PDF of manuscript
```
