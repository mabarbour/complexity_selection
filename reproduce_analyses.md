---
title: "Supplementary Material"
author: "Matthew A. Barbour"
date: "2018-11-13"
output:
  prettydoc::html_pretty:
    theme: architect
    highlight: github
#output:
#  pdf_document:
#    df_print: paged
#    fig_caption: true
#    citation_package: natbib
#    includes:  
#      in_header: preamble-latex.tex
#bibliography: references_ECD_supp_mat.bib
fontsize: 12pt
---



# Evaluating assumption of multivariate normality

We used graphical checks to evaluate whether our transformations of trait values resulted in a multivariate normal distribution. Figure S\ref{fig:Univariate_histograms} shows that our transformations resulted in approximately normal distributions for each phenotypic trait. Note also that in the multivariate quantile-quantile (Q-Q) plot, most points fall along the expected line (fig. S\ref{fig:Multivariate_QQ}), suggesting that our transformations provide a reasonable approximation of a multivariate normal distribution. 

![\label{fig:Univariate_histograms}Histograms of each phenotypic trait after transformation. The red line illustrates a normal distribution.](manuscript_files/figure-latex/Univariate_histograms-1.pdf) 

![\label{fig:Multivariate_QQ}Multivariate quantile-quantile (Q-Q) plot to assess deviations from multivariate normality (black line).](manuscript_files/figure-latex/Multivariate_QQ-1.pdf) 

# Effect of food-web treatment on trait-fitness relationships and selection gradients

We write the model in a way the independently estimates the effect of food-web treatment, each trait, and all two-way and three-way statistical interactions, on larval survival.


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

Note that the resulting estimates and confidence intervals are useful for determining whether trait-fitness relationships differ from zero, but not whether they differ between food-web treatments. For the later, we calculate the differences between each food-web treatment from the bootstrapped samples.




```r
linear_foodweb_model <- glmer(
  gall_survival ~ 
    -1 + Foodweb + 
    Foodweb:(sc.Diam + sc.log.Clutch + sc.sqrt.Pref) +
    (1|Genotype/Plant_Position/Gall_Number),
  data = gall_selection.df,
  family = binomial(link = logit), control=glmerControl(optimizer = "bobyqa"))
```




To estimate biased selection on chamber diameter, we subset our data to only include multi-chambered galls where there was variability in larval survival. We then fit a reduced model to estimate the bias in the logistic regression coefficient of chamber diameter in each food web.

<!-- We quantified biased selection on chamber diameter in the following way. First, we focus on multi-chambered galls where there is evidence of both survival and parasitism. The idea is that galls from the same clutch should be similar in size; therefore, any selection on diameter would be apparent. If we assume that this is all due to the effect of parasitism on larval development, truncating gall diameter, then we can consider this new selection differential to be the result of confounding factors. Note that this method likely overestimates the confounding effect of parasitism, since we are assuming that any heterogeneity in chamber diameter (for this subset) is due to parasitism. -->


```r
biased_foodweb_df <- gall_selection.df %>%
  group_by(Foodweb, Gall_Number) %>%
  mutate(mean_survival = mean(gall_survival)) %>%
  filter(mean_survival > 0, mean_survival < 1) %>%
  ungroup()

# using linear model, but I'm only interested in the chamber diameter coefficient
biased_foodweb_model <- update(linear_foodweb_model, data=biased_foodweb_df)
  #glmer(
  #gall_survival ~ -1 + Foodweb + 
  #  Foodweb:sc.Diam + 
  #  (1|Genotype/Plant_Position/Gall_Number),
  #data = biased_foodweb_df,
  #family = binomial(link = logit), control=glmerControl(optimizer = "bobyqa"))

biased_foodweb_confint <- tidy(biased_foodweb_model, conf.int=TRUE) %>% filter(group=="fixed")
```

```
## Warning in bind_rows_(x, .id): binding factor and character vector,
## coercing into character vector
```

```
## Warning in bind_rows_(x, .id): binding character and factor vector,
## coercing into character vector
```









# Partitioning the contribution of egg and larval parasitoids to selection gradients

Our simple food-web treatment allows us to estimate the unique contribution of egg parasitoids to selection on *Iteomyia* traits. To estimate the unique contribution of larval parasitoids, we subset our data so that our complex food-web treatment only contained attack by larval parasitoids (and gall survival). We then fit the same models as previously, including one to estimate bias.


```r
# excludes cases of egg-parasitism from Complex food web
egglarval_df <- filter(gall_selection.df, Foodweb == "Simple" | Foodweb == "Complex" & platy < 1) 

egglarval_model <- update(foodweb_model, data=egglarval_df)

biased_egglarval_df <- egglarval_df %>%
  group_by(Foodweb, Gall_Number) %>%
  mutate(mean_survival = mean(gall_survival)) %>%
  filter(mean_survival > 0, mean_survival < 1) %>%
  ungroup()

biased_egglarval_model <- update(biased_foodweb_model, data=biased_egglarval_df)
```




```r
linear_egglarval_model <- update(linear_foodweb_model, data=egglarval_df)
```











# Reproduce Figure 2 of main manuscript

We combine our estimates of selection gradients for each food-web treatment as well as the contribution of larval parasitoids to selection in the complex food web.

![**Selection gradients**. Estimates of standardized selection gradients in complex (orange) and simple (blue) food webs. The contribution of larval parasitoids (yellow) was estimated with a subset of the complex food-web data that only contained attacks from larval parasitoids (and gall survival). Points and lines correspond to estimates of the mean and 95% confidence intervals, respectively. Overlapping confidence intervals with zero (dotted line) indicate no strong evidence of selection.](manuscript_files/figure-latex/Figure-Selection-Gradients-1.pdf) 

# Partitioning the components of selection gradients

Selection gradients are influenced by both trait-fitness relationships and population mean fitness. Here, we partition selection gradients into these underlying components.

![**Partitioning Selection Gradients.**](manuscript_files/figure-latex/Plot foodwebegglarval coefficients-1.pdf) 

# Estimating selection on the egg parasitoid *Platygaster*


```r
# excludes cases of larval-parasitism from Complex food web
#eggegg_df <- #filter(gall_selection.df, Foodweb == "Simple" | Foodweb == "Complex" & ectos < 1) %>%
  #mutate(gall_survival = ifelse(gall_survival==1,0,1)) # attempt to look at survival of egg parasitoid

# convert "gall_survival" to egg parasitoid survival. Note that both Iteomyia pupa and larva parasitoids result in 0. We do not change the name of the response variable to make clear that we are using the same statistical models as for Iteomyia.
eggegg_df <- mutate(gall_selection.df, gall_survival = ifelse(egg_parasitoid==1,1,0))

eggegg_model <- update(foodweb_model, data=eggegg_df)

biased_eggegg_df <- eggegg_df %>%
  group_by(Foodweb, Gall_Number) %>%
  mutate(mean_survival = mean(gall_survival)) %>%
  filter(mean_survival > 0, mean_survival < 1) %>%
  ungroup()

biased_eggegg_model <- update(biased_foodweb_model, data=biased_eggegg_df)
```




```r
linear_eggegg_model <- update(linear_foodweb_model, data=eggegg_df)
```





![](manuscript_files/figure-latex/Plot egg vs. egg ptoid coefficients-1.pdf)<!-- --> 



![](manuscript_files/figure-latex/Plot egg vs. egg ptoid gradients-1.pdf)<!-- --> 




# Reproduce Figure 4 - Univariate adaptive landscape














![](manuscript_files/figure-latex/Figure-Univariate-Landscapes-1.pdf)<!-- --> 



![](manuscript_files/figure-latex/Egg-Parasitoid-Selection-1.pdf)<!-- --> 







# Multivariate fitness landscapes









Reproduce Figure 5 in main manuscript.

![Fitness landscapes of gall traits in complex vs. simple food webs. Each panel corresponds to a different combination of traits: clutch size and gall diameter (A); clutch size and Oviposition preference (B); Oviposition preference and gall diameter (C). Note that traits for all plots range 1 SD below and above the mean (=0).](manuscript_files/figure-latex/Figure-Multivariate-Landscapes-1.pdf) 

<!--
# Effect of Food-Web Treatment on Larval Survival

First we fit the full generalized linear mixed model to our data. This will enable us to evaluate whether there is any evidence that food-web treatment alters nonlinear selection on gall phenotypic variation. Note that we are not interested in directional selection yet, because we need to remove the nonlinear terms to evaluate its effects ([Lande and Arnold 1983](https://doi.org/10.1111/j.1558-5646.1983.tb00236.x)).



We use parametric bootstrapping (n=1000) to evaluate the statistical significance of the effect of food-web treatment on nonlinear terms (i.e. Foodweb:trait~i~^2^ or Foodweb:trait~i~:trait~j~). Food-web treatment does influence the quadratic effect of of Oviposition preference on larval survival. None of the other terms statistically differ from zero.





Therefore, we are going to drop them from the model. Note that we are still fitting the full curvature of the fitness landscape, but that we aren't allowing some aspects of it to vary between food-web treatments (since there was no evidence for it)





This represents our nonlinear regression coefficients.



Note that the non-overlapping 95% CI for *FoodwebSimple:I(sc.sqrt.Pref^2)* indicates that food-web treatment alters quadractic selection on Oviposition preference; however, this coefficient is not useful for actually calculating the selection gradient. To do this, we refit a model so rather than estimating the change in slope (currently), it will actually estimate the slope. See [Schielzeth 2010](https://doi.org/10.1111/j.2041-210X.2010.00012.x) for a detailed explanation of this approach.





Importantly we see now that although quadratic selection varies between treatments, it is not different from zero in the simple food web.



Now we fit a reduced model to assess the effects of whether food-web treatments alters components of directional selection



We then use parametric bootstrapping to test whether food-web treatment alters the linear relationship between each trait and larval survival.







There is clear evidence that food-web treatment alters the linear effect of clutch size on larval survival, but it has no effect on chamber diameter or Oviposition preference. Therefore, we fit a reduced model that removes these non-significant effects: 



This reduced model returns the same qualitative result as before. Specifically, that the linear effect of clutch size depends on food-web treatment. Note also that larval survival is significantly higher in the simple food web compared to the complex food web (dotted line). These are the regression coefficients we will use when calculating directional selection gradients ($\beta$).





As we did with model 2, we refit the model by removing the slope term so that we can fit different slopes for each interaction term. 










So, for plotting the results, we are going to use the focused nonlinear model. This provides the best representation of the fitness landscape (although we can't use it to assess to quantify the effects of directional selection.

Note that while there is no evidence to let the slopes vary among treatments, there is a clear effect of food web treatment on fitness. This could influence estimates of selection by magnifying them in the complex vs. simple food web.

So to calculate the selection gradients, we used Janzen and Stern's method to transform the regression coefficients.

Basically, I need to just multiply all of the coefficients by the simple and complex fitness. Except, however, for I(sc.sqrt.Pref) and sc.log.Clutch. Here, we need the different coefficients.

# Selection Gradients

We used the method of [Janzen and Stern (1998)](https://doi.org/10.1111/j.1558-5646.1998.tb02237.x) to calculate selection gradients. Specifically, we first multiply each regression coefficient (and confidence interval) by the average of $W(z)[1-W(z)]$, where $W(z)$ corresponds to the predicted absolute fitness of an individual (in this case 1 if the larva survived, 0 if not) given its multivariate phenotype ($z$). We then divided this value by the mean absolute fitness ($\bar W$) so our results are in terms of relative fitness ($w$). Note that we doubled all quadratic selection gradients ($\gamma_{z_i,z_i}$) to put them on the same scale as the directional and correlational selection gradients (detailed explanation in [Stinchcombe et al. 2008](https://doi.org/10.1111/j.1558-5646.2008.00449.x)).











Contributions of larval parasitoids to selection gradients

To identify the mechanisms by which food-web treatment altered selection gradients, we refit our models so that the complex food-web treatment only contained cases of parasitisms from larval parasitoids. Therefore, these analyses identify the unique contribution of each parasitoid guild to the selection gradients (egg parasitoid from simple food web; larval parasitoid from complex food web). Note that cases of larval parasitism may have been initially parasitized by an egg parasitoid.

We restrict our model analyses to those where we observed different effects of food-web treatment on larval survival from the prior analyses (i.e. model 2 and model 4).





This represents our nonlinear regression coefficients.


 To put all of the coefficients on the same scale, we refit it as we did in the food-web treatment




This represents our nonlinear regression coefficients.



Now we assess their contribution to the effect of food-web treatment on the linear terms.



This reduced model returns the same qualitative result as before. Specifically, that the linear effect of clutch size depends on food-web treatment. Note also that larval survival is significantly higher in the simple food web compared to the complex food web (dotted line). These are the regression coefficients we will use when calculating directional selection gradients ($\beta$).



Here, we observe clear evidence of conflicting selection pressures on clutch size. This results in net zero selection gradient acting on clutch size in the complex food web.




This reduced model returns the same qualitative result as before. Specifically, that the linear effect of clutch size depends on food-web treatment. Note also that larval survival is significantly higher in the simple food web compared to the complex food web (dotted line). These are the regression coefficients we will use when calculating directional selection gradients ($\beta$).



Here, we observe clear evidence of conflicting selection pressures on clutch size. This results in net zero selection gradient acting on clutch size in the complex food web.


# CALCULATE SELECTION GRADIENTS FOR EGG VS. LARVAL PARASITOIDS











# Quantifying biased selection on chamber diameter

We quantified biased selection on chamber diameter in the following way. First, we focus on multi-chambered galls where there is evidence of both survival and parasitism. The idea is that galls from the same clutch should be similar in size; therefore, any selection on diameter would be apparent. If we assume that this is all due to the effect of parasitism on larval development, truncating gall diameter, then we can consider this new selection differential to be the result of confounding factors. Note that this method likely overestimates the confounding effect of parasitism, since we are assuming that any heterogeneity in chamber diameter (for this subset) is due to parasitism. 

We calculated the bias in selection differential on chamber diameter by calculating the difference in mean chamber diameter before (i.e. all larva within the same gall) and after selection (only surviving larva). We then used a one-sample t-test to determine whether the selection differential was significantly different from zero. We conducted separate t-tests for each food-web treatment.

We can then subtract this "biased selection differential" from the observed directional selection gradient on chamber diameter to adjust for the confounding effects of parasitism. 

















## Figure for manuscript

Only use selection gradients for larval parasitoids, use both complex and simple food webs, and only bias adjusted for chamber diameter.



Standardized selection gradients estimate selection on a trait in terms of the effects on relative fitness in units of (phenotypic) standard deviations of the trait, allowing direct comparisons among traits, fitness components, and study systems (Kingsolver et al. 2001, Am. Nat.).




## Selection on the egg parasitoid *Platygaster*

We cannot calculate selection on the egg parasitoids *Platygaster* with the same method as we did for the herbivore *Iteomyia*. This is because the unique biology of host-parasitoid interactions enables us to determine whether the host (herbivore) has survived or not. With the egg parasitoid, we only have a record of the gall it was found in. But, we can take advantage of our treatment which excludes the guild of larval parasitoids. Specifically, we can compare the mean gall phenotype that egg parasitoids are reared from in the absence of larval parasitoids ("before" selection) to the mean gall phenotype in the presence of larval parasitoids ("after" selection). Therefore, we can use this cross-sectional data to quantify directional and quadratic selection differentials acting on the gall phenotypes. In essence, the gall phenotype becomes the extended phenotype of the egg parasitoid, since its attributes can influence the egg parasitoid's survival to its intraguild predator (larval parasitoid).

To quantify the selection differential acting on the egg parasitoid, we used separate linear mixed models with the gall phenotype as the response variable, food-web treatment as the fixed effect, and then the same suite of random effects as for our primary analyses of *Iteomyia* survival (i.e. gall ID, nested within plant ID, nested within plant genotype).

We found no evidence of directional or quadratic selection on either chamber diameter or clutch size; however, we observed positive directional selection on *Iteomyia* Oviposition preference as well stabilizing selection on preference.

Appears to be selection for Platygaster to attack larva at high densities in more complex food webs (and a decrease in the variance). This may be a result of the poorer searching ability of parasitoids at very high gall densities, potentially due to saturation. This suggests that the fitness landscape may be dynamic because these selective effects (assuming there is heritable variation in the egg parasitoid) will change the following year. This could be an important discussion point to take into account for future work. This is different from simply altering the strength of selection, which will occur as species move along the fitness landscape, but this result suggests that the nature of the fitness landscape may actually be changing. For example, larval parasitoids could be driving selection on Platygaster to be very efficient at foraging, which will dampen once this pressure is removed, and then could also dampen the selection acting on the galling herbivore.













Note that equation 3 in Phillips and Arnold 1998 is good justification for why quadratic selection cofficients should be multiplied by two (but not correlational ones?)

The G-matrix is essentially a scalar. For my data, since there is no correlational selection, it doesn't matter if there is genetic covariance between the traits, as this will not effect the change in genetic covariances within a generation (because there is no selection). Also, the G_matrix does not qualitatively alter the conclusions that there will be a decrease in additive genetic variance in gall diameter due to the strong directional selection; however, there will be an increase in the additive genetic variance in Oviposition preference. Since this is an experiment, we can assume that the G-matrix is the same between treatments

Remember that the diagonals refer to additive genetic CO-variances. Thus, positive or negative values for delta_G give insight to whether selection will act to integrate traits (positive covariance) or create trade-offs (negative covariance).

Assuming that there is positive additive genetic variance and co-variance in these traits (i.e. all values of G-matrix are positive), then the curvature of the fitness landscape can give insight to qualitative changes in the G matrix. This is because the G-matrix acts as a scalar of changes in the fitness landscape.

Create a simple rule. If the confidence interval of the selection gradient overlaps zero, then I set it to zero. Otherwise, I set it to our best estimate of the gradient-->


<!--
## Slope of adaptive landscape

$$\text{Slope} = \beta = \begin{pmatrix} \beta_{diam} \\ \beta_{clutch} \\ \beta_{pref} \end{pmatrix}$$

$$\text{Slope}_{complex} = \begin{pmatrix}  \\  \\  \end{pmatrix}$$
$$\text{Slope}_{simple} =  \begin{pmatrix}  \\  \\  \end{pmatrix}$$

## Curvature of adaptive landscape

$$\text{Curvature} = \gamma - \beta \beta^\text{T}$$

$$\gamma = \begin{pmatrix} \gamma_{diam} & \gamma_{diam,clutch} & \gamma_{diam,pref} \\ \gamma_{clutch,diam} & \gamma_{clutch} & \gamma_{clutch,pref} \\ \gamma_{pref,diam} & \gamma_{pref,clutch} &\gamma_{pref} \end{pmatrix}$$
$$\gamma_{complex} = \begin{pmatrix}  &  &  \\   &  &  \\   &  &  \end{pmatrix}$$
$$\gamma_{simple} = \begin{pmatrix}  &  &  \\   &  &  \\   &  &  \end{pmatrix}$$


$$\text{Curvature}_{complex} = \begin{pmatrix}  &  &  \\   &  &  \\   &  &  \end{pmatrix}$$

$$\text{Curvature}_{simple} = \begin{pmatrix}  &  &  \\   &  &  \\   &  &  \end{pmatrix}$$

## Observed phenotypic variance-covariance matrix



$$\text{P} = \begin{pmatrix}  &  &  \\   &  &  \\   &  &  \end{pmatrix}$$
Note that the variance of each phenotypic trait is standardized to 1. If we assume that each phenotypic trait has a narrow-sense heritability ($h^2$) of , then the additive genetic variance of each trait is equal to its narrow-sense heritability. The resulting G-matrix is:

$$\text{G} = \begin{pmatrix}  &  &  \\   &  &  \\   &  &  \end{pmatrix}$$



# Change in complex G-matrix

$$\text{G}_{complex} = \begin{pmatrix}  &  &  \\   &  &  \\   &  &  \end{pmatrix}$$

# Change in simple G-matrix

$$\text{G}_{simple} = \begin{pmatrix}  &  &  \\   &  &  \\   &  &  \end{pmatrix}$$






## Walk across the adaptive landscape




## Compare evolvabilities of hypothetical G-matrix over time

Since my covariance is standardized, it should be okay to calculate the mean squared correlation as an overall measure of trait integration.




## More walks across the adaptive landscape







## Multivariate trait differences




## How does the adaptive landscape influence the rate of evolution?

Here we are going to see how fast it took for the *Iteomyia* to hit the edge of the adaptive landscape across many different G-matrices.


In general, *Iteomyia* hit the edge of the adaptive landscape faster in the complex vs. simple food web.














$$logit(E(W))=\mu+\beta_{diameter}+\beta_{clutch}+\beta_{preference}+$$
$$logit(E(W))=\mu+treatment+s(diam):treat+s(clutch)+s(pref)+diam:cluch+diam:pref+clutch:pref+(1|Genotype/Plant/Gall)$$

We used generalized linear mixed models (GLMMs, @Bolker2009) to test the effects of food-web complexity on the shape of the fitness landscape. 
larval survival (0 or 1) was our response variable and measure of fitness. 

Need to clean up the language part of this section without it being too complicated and confusing. I should also include tests of non-linear selection gradients.

We specified our food-web treatment, each gall trait, two-way interactions between gall traits, as well interactions between gall trait and food-web treatment, as fixed effects to fully explore the effects of food-web complexity on the fitness landscape. 
This analysis implicitly assumes that selection is linear, which we felt was a necessary trade-off for exploring the shape of the fitness landscape. 
To account for the correlated structure of clutch size (gall level) and Oviposition preference (plant level) as well as any other independent effects of willow genotype on parasitism rates, we specified gall ID nested within plant ID nested within plant genotype as random intercepts in our statistical models. 






## Contributions of egg parasitoids to selection gradients

To identify the mechanisms by which food-web treatment altered selection gradients, we refit our models so that the complex food-web treatment only contained cases of parasitisms from egg parasitoids. Therefore, these analyses identify the change in apparent selection imposed by the egg parasitoid (egg parasitoid from simple food web; egg parasitoid from complex food web).









There is no strong evidence that egg parasitism differs for nonlinear terms between complex and simple food webs. Although close perhaps for sc.Diam.



Now we assess their contribution to the effect of food-web treatment on the linear terms.






Here, we observe clear evidence of conflicting selection pressures on clutch size. This results in net zero selection gradient acting on clutch size in the complex food web.




This reduced model returns the same qualitative result as before. Specifically, that the linear effect of clutch size depends on food-web treatment. Note also that egg survival is significantly higher in the simple food web compared to the complex food web (dotted line). These are the regression coefficients we will use when calculating directional selection gradients ($\beta$).



Here, we observe clear evidence of conflicting selection pressures on clutch size. This results in net zero selection gradient acting on clutch size in the complex food web.


## CALCULATE SELECTION GRADIENTS FOR EGG VS. egg PARASITOIDS










-->
