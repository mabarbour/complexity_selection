---
title: "Reproduce analyses reported in: *Phenotypic evolution is more constrained in simpler food webs*"
author: "Matthew A. Barbour"
date: "2019-01-11"
output: github_document
# Trying github_document for easy viewing on github
#  prettydoc::html_pretty:
#    theme: architect
#    highlight: github
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


```r
biased_foodweb_df <- gall_selection.df %>%
  group_by(Foodweb, Gall_Number) %>%
  mutate(mean_survival = mean(gall_survival)) %>%
  filter(mean_survival > 0, mean_survival < 1) %>%
  ungroup()

# using linear model, but I'm only interested in the chamber diameter coefficient
biased_foodweb_model <- glmer(
  gall_survival ~ -1 + Foodweb + 
    Foodweb:sc.Diam + 
    (1|Genotype/Plant_Position/Gall_Number),
  data = biased_foodweb_df,
  family = binomial(link = logit), control=glmerControl(optimizer = "bobyqa"))
```

```
## singular fit
```

```r
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

```r
knitr::kable(biased_foodweb_confint)
```


\begin{tabular}{l|r|r|r|r|r|r|l}
\hline
term & estimate & std.error & statistic & p.value & conf.low & conf.high & group\\
\hline
FoodwebComplex & 0.0261288 & 0.1596126 & 0.1637011 & 0.8699664 & -0.2867063 & 0.3389638 & fixed\\
\hline
FoodwebSimple & -0.1181035 & 0.2044004 & -0.5778045 & 0.5633961 & -0.5187209 & 0.2825139 & fixed\\
\hline
FoodwebComplex:sc.Diam & 0.3603278 & 0.1564728 & 2.3028133 & 0.0212893 & 0.0536466 & 0.6670089 & fixed\\
\hline
FoodwebSimple:sc.Diam & 0.4158534 & 0.2066780 & 2.0120841 & 0.0442111 & 0.0107721 & 0.8209348 & fixed\\
\hline
\end{tabular}





\begin{tabular}{l|l|r|r|r|l|l}
\hline
type & term & conf.high & conf.low & mean & P\_cutoff & coefficient\_type\\
\hline
Complex & (Intercept) & 0.5561001 & 0.2728743 & 0.4188447 & * & Mean fitness\\
\hline
Complex & I(sc.Diam\textasciicircum{}2) & 0.5458940 & -0.0565636 & 0.2340236 &  & Quadratic\\
\hline
Complex & I(sc.log.Clutch\textasciicircum{}2) & 0.2673715 & -0.4442637 & -0.0934645 &  & Quadratic\\
\hline
Complex & I(sc.sqrt.Pref\textasciicircum{}2) & 1.1078055 & 0.1575109 & 0.5813342 & * & Quadratic\\
\hline
Complex & sc.Diam & 1.6397205 & 0.7475711 & 1.1511781 & * & Linear\\
\hline
Complex & sc.Diam:sc.log.Clutch & 0.2400422 & -0.5655092 & -0.1426436 &  & Correlational\\
\hline
Complex & sc.Diam:sc.sqrt.Pref & 0.0357144 & -0.9376507 & -0.4497678 &  & Correlational\\
\hline
Complex & sc.log.Clutch & 0.5742409 & -0.1782814 & 0.2015568 &  & Linear\\
\hline
Complex & sc.log.Clutch:sc.sqrt.Pref & 0.5409796 & -0.3072356 & 0.1084018 &  & Correlational\\
\hline
Complex & sc.sqrt.Pref & 0.1599006 & -0.9617830 & -0.4056342 &  & Linear\\
\hline
Diff & (Intercept) & 0.4244776 & 0.1154952 & 0.2670227 & * & Mean fitness\\
\hline
Diff & I(sc.Diam\textasciicircum{}2) & 0.4842451 & -0.4243171 & 0.0280245 &  & Quadratic\\
\hline
Diff & I(sc.log.Clutch\textasciicircum{}2) & 0.3296182 & -0.7464913 & -0.2141728 &  & Quadratic\\
\hline
Diff & I(sc.sqrt.Pref\textasciicircum{}2) & 0.0254895 & -1.1655746 & -0.5434675 &  & Quadratic\\
\hline
Diff & sc.Diam & 0.5314249 & -0.6157623 & -0.0390675 &  & Linear\\
\hline
Diff & sc.Diam:sc.log.Clutch & 0.3596399 & -0.7890798 & -0.2172148 &  & Correlational\\
\hline
Diff & sc.Diam:sc.sqrt.Pref & 1.0033748 & -0.2509312 & 0.3671110 &  & Correlational\\
\hline
Diff & sc.log.Clutch & -0.1350917 & -1.2435804 & -0.6762277 & * & Linear\\
\hline
Diff & sc.log.Clutch:sc.sqrt.Pref & 0.4674799 & -0.7392477 & -0.1070403 &  & Correlational\\
\hline
Diff & sc.sqrt.Pref & 0.2620370 & -1.0910036 & -0.4293349 &  & Linear\\
\hline
Simple & (Intercept) & 0.8166066 & 0.5532049 & 0.6858674 & * & Mean fitness\\
\hline
Simple & I(sc.Diam\textasciicircum{}2) & 0.6252214 & -0.0523585 & 0.2620482 &  & Quadratic\\
\hline
Simple & I(sc.log.Clutch\textasciicircum{}2) & 0.1041150 & -0.7351263 & -0.3076373 &  & Quadratic\\
\hline
Simple & I(sc.sqrt.Pref\textasciicircum{}2) & 0.4378315 & -0.3910659 & 0.0378668 &  & Quadratic\\
\hline
Simple & sc.Diam & 1.6299262 & 0.6939471 & 1.1121106 & * & Linear\\
\hline
Simple & sc.Diam:sc.log.Clutch & 0.0772856 & -0.8152498 & -0.3598584 &  & Correlational\\
\hline
Simple & sc.Diam:sc.sqrt.Pref & 0.3179567 & -0.5082598 & -0.0826568 &  & Correlational\\
\hline
Simple & sc.log.Clutch & -0.0787164 & -0.8991793 & -0.4746709 & * & Linear\\
\hline
Simple & sc.log.Clutch:sc.sqrt.Pref & 0.4122045 & -0.4139481 & 0.0013615 &  & Correlational\\
\hline
Simple & sc.sqrt.Pref & -0.3639175 & -1.3831022 & -0.8349691 & * & Linear\\
\hline
\end{tabular}


![](manuscript_files/figure-latex/Plot food-web coefficients-1.pdf)<!-- --> 




\begin{tabular}{l|l|r|r|r|l|l}
\hline
type & term & conf.high & conf.low & mean & P\_cutoff & gradient\_type\\
\hline
Complex & I(sc.Diam\textasciicircum{}2) & 0.3271709 & -0.0339003 & 0.1402575 &  & Quadratic\\
\hline
Complex & I(sc.log.Clutch\textasciicircum{}2) & 0.1602439 & -0.2662608 & -0.0560161 &  & Quadratic\\
\hline
Complex & I(sc.sqrt.Pref\textasciicircum{}2) & 0.6639415 & 0.0944011 & 0.3484113 & * & Quadratic\\
\hline
Complex & sc.Diam & 0.4913672 & 0.2240211 & 0.3449680 & * & Directional\\
\hline
Complex & sc.Diam:sc.log.Clutch & 0.0719323 & -0.1694634 & -0.0427453 &  & Correlational\\
\hline
Complex & sc.Diam:sc.sqrt.Pref & 0.0107023 & -0.2809813 & -0.1347798 &  & Correlational\\
\hline
Complex & sc.log.Clutch & 0.1720800 & -0.0534247 & 0.0603996 &  & Directional\\
\hline
Complex & sc.log.Clutch:sc.sqrt.Pref & 0.1621128 & -0.0920678 & 0.0324842 &  & Correlational\\
\hline
Complex & sc.sqrt.Pref & 0.0479166 & -0.2882129 & -0.1215545 &  & Directional\\
\hline
Diff & I(sc.Diam\textasciicircum{}2) & 0.1766276 & -0.2701726 & -0.0416837 &  & Quadratic\\
\hline
Diff & I(sc.log.Clutch\textasciicircum{}2) & 0.1908814 & -0.3128531 & -0.0597068 &  & Quadratic\\
\hline
Diff & I(sc.sqrt.Pref\textasciicircum{}2) & -0.0401891 & -0.6792307 & -0.3341671 & * & Quadratic\\
\hline
Diff & sc.Diam & -0.0009627 & -0.2885423 & -0.1357985 & * & Directional\\
\hline
Diff & sc.Diam:sc.log.Clutch & 0.1172844 & -0.1668003 & -0.0249380 &  & Correlational\\
\hline
Diff & sc.Diam:sc.sqrt.Pref & 0.2829095 & -0.0403247 & 0.1192334 &  & Correlational\\
\hline
Diff & sc.log.Clutch & -0.0182048 & -0.2927574 & -0.1496772 & * & Directional\\
\hline
Diff & sc.log.Clutch:sc.sqrt.Pref & 0.1127559 & -0.1892195 & -0.0322282 &  & Correlational\\
\hline
Diff & sc.sqrt.Pref & 0.1413158 & -0.2070588 & -0.0354893 &  & Directional\\
\hline
Simple & I(sc.Diam\textasciicircum{}2) & 0.2351875 & -0.0196955 & 0.0985738 &  & Quadratic\\
\hline
Simple & I(sc.log.Clutch\textasciicircum{}2) & 0.0391646 & -0.2765300 & -0.1157229 &  & Quadratic\\
\hline
Simple & I(sc.sqrt.Pref\textasciicircum{}2) & 0.1646976 & -0.1471060 & 0.0142442 &  & Quadratic\\
\hline
Simple & sc.Diam & 0.3065620 & 0.1305199 & 0.2091695 & * & Directional\\
\hline
Simple & sc.Diam:sc.log.Clutch & 0.0145361 & -0.1533349 & -0.0676834 &  & Correlational\\
\hline
Simple & sc.Diam:sc.sqrt.Pref & 0.0598024 & -0.0955952 & -0.0155464 &  & Correlational\\
\hline
Simple & sc.log.Clutch & -0.0148052 & -0.1691207 & -0.0892777 & * & Directional\\
\hline
Simple & sc.log.Clutch:sc.sqrt.Pref & 0.0775288 & -0.0778567 & 0.0002561 &  & Correlational\\
\hline
Simple & sc.sqrt.Pref & -0.0684468 & -0.2601385 & -0.1570438 & * & Directional\\
\hline
\end{tabular}


![](manuscript_files/figure-latex/Plot food-web gradients-1.pdf)<!-- --> 

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

```
## singular fit
```




```r
linear_egglarval_model <- update(linear_foodweb_model, data=egglarval_df)
```






\begin{tabular}{l|l|r|r|r|l|l}
\hline
type & term & conf.high & conf.low & mean & P\_cutoff & coefficient\_type\\
\hline
Complex & (Intercept) & 0.8322366 & 0.4720168 & 0.6700690 & * & Mean fitness\\
\hline
Complex & I(sc.Diam\textasciicircum{}2) & 0.5096068 & -0.1785902 & 0.1643608 &  & Quadratic\\
\hline
Complex & I(sc.log.Clutch\textasciicircum{}2) & 0.3854401 & -0.3636687 & 0.0088227 &  & Quadratic\\
\hline
Complex & I(sc.sqrt.Pref\textasciicircum{}2) & 1.0774509 & -0.0100810 & 0.4726738 &  & Quadratic\\
\hline
Complex & sc.Diam & 1.7594724 & 0.6785354 & 1.1651462 & * & Linear\\
\hline
Complex & sc.Diam:sc.log.Clutch & 0.5539638 & -0.3258649 & 0.0997995 &  & Correlational\\
\hline
Complex & sc.Diam:sc.sqrt.Pref & 0.3476120 & -0.8633409 & -0.2786942 &  & Correlational\\
\hline
Complex & sc.log.Clutch & 1.2008051 & 0.1994744 & 0.6578867 & * & Linear\\
\hline
Complex & sc.log.Clutch:sc.sqrt.Pref & 0.3109048 & -0.7446019 & -0.2056841 &  & Correlational\\
\hline
Complex & sc.sqrt.Pref & -0.1584967 & -1.7124265 & -0.9080523 & * & Linear\\
\hline
Diff & (Intercept) & 0.2252990 & -0.1995580 & 0.0172334 &  & Mean fitness\\
\hline
Diff & I(sc.Diam\textasciicircum{}2) & 0.5719119 & -0.3865508 & 0.0846807 &  & Quadratic\\
\hline
Diff & I(sc.log.Clutch\textasciicircum{}2) & 0.1931778 & -0.8612449 & -0.3012207 &  & Quadratic\\
\hline
Diff & I(sc.sqrt.Pref\textasciicircum{}2) & 0.0930201 & -1.2916097 & -0.5463166 &  & Quadratic\\
\hline
Diff & sc.Diam & 0.5392765 & -0.6909172 & -0.0884426 &  & Linear\\
\hline
Diff & sc.Diam:sc.log.Clutch & 0.1465344 & -1.0612704 & -0.4586126 &  & Correlational\\
\hline
Diff & sc.Diam:sc.sqrt.Pref & 0.8368888 & -0.5233818 & 0.1873353 &  & Correlational\\
\hline
Diff & sc.log.Clutch & -0.7885031 & -2.2567213 & -1.4550976 & * & Linear\\
\hline
Diff & sc.log.Clutch:sc.sqrt.Pref & 0.8736969 & -0.4036615 & 0.2419612 &  & Correlational\\
\hline
Diff & sc.sqrt.Pref & 0.8028004 & -0.9331559 & -0.0939651 &  & Linear\\
\hline
Simple & (Intercept) & 0.8576103 & 0.4955104 & 0.6873024 & * & Mean fitness\\
\hline
Simple & I(sc.Diam\textasciicircum{}2) & 0.6067336 & -0.0419911 & 0.2490416 &  & Quadratic\\
\hline
Simple & I(sc.log.Clutch\textasciicircum{}2) & 0.0652989 & -0.6840806 & -0.2923980 &  & Quadratic\\
\hline
Simple & I(sc.sqrt.Pref\textasciicircum{}2) & 0.3694016 & -0.5284587 & -0.0736428 &  & Quadratic\\
\hline
Simple & sc.Diam & 1.6071684 & 0.6125782 & 1.0767036 & * & Linear\\
\hline
Simple & sc.Diam:sc.log.Clutch & 0.0258376 & -0.8224415 & -0.3588130 &  & Correlational\\
\hline
Simple & sc.Diam:sc.sqrt.Pref & 0.2857568 & -0.4686261 & -0.0913589 &  & Correlational\\
\hline
Simple & sc.log.Clutch & -0.2933414 & -1.3339773 & -0.7972110 & * & Linear\\
\hline
Simple & sc.log.Clutch:sc.sqrt.Pref & 0.3944040 & -0.3566626 & 0.0362771 &  & Correlational\\
\hline
Simple & sc.sqrt.Pref & -0.4224492 & -1.6164371 & -1.0020174 & * & Linear\\
\hline
\end{tabular}


![](manuscript_files/figure-latex/Plot egglarval coefficients-1.pdf)<!-- --> 




\begin{tabular}{l|l|r|r|r|l|l}
\hline
type & term & conf.high & conf.low & mean & P\_cutoff & gradient\_type\\
\hline
Complex & I(sc.Diam\textasciicircum{}2) & 0.1973813 & -0.0691717 & 0.0636604 &  & Quadratic\\
\hline
Complex & I(sc.log.Clutch\textasciicircum{}2) & 0.1492889 & -0.1408564 & 0.0034172 &  & Quadratic\\
\hline
Complex & I(sc.sqrt.Pref\textasciicircum{}2) & 0.4173191 & -0.0039046 & 0.1830763 &  & Quadratic\\
\hline
Complex & sc.Diam & 0.3407401 & 0.1314054 & 0.2256426 & * & Directional\\
\hline
Complex & sc.Diam:sc.log.Clutch & 0.1072808 & -0.0631071 & 0.0193272 &  & Correlational\\
\hline
Complex & sc.Diam:sc.sqrt.Pref & 0.0673187 & -0.1671949 & -0.0539720 &  & Correlational\\
\hline
Complex & sc.log.Clutch & 0.2325483 & 0.0386303 & 0.1274066 & * & Directional\\
\hline
Complex & sc.log.Clutch:sc.sqrt.Pref & 0.0602099 & -0.1441999 & -0.0398329 &  & Correlational\\
\hline
Complex & sc.sqrt.Pref & -0.0306945 & -0.3316291 & -0.1758537 & * & Directional\\
\hline
Diff & I(sc.Diam\textasciicircum{}2) & 0.2270190 & -0.1494159 & 0.0348792 &  & Quadratic\\
\hline
Diff & I(sc.log.Clutch\textasciicircum{}2) & 0.0740182 & -0.3383712 & -0.1191118 &  & Quadratic\\
\hline
Diff & I(sc.sqrt.Pref\textasciicircum{}2) & 0.0380564 & -0.5048781 & -0.2122150 &  & Quadratic\\
\hline
Diff & sc.Diam & 0.1108178 & -0.1306718 & -0.0126303 &  & Directional\\
\hline
Diff & sc.Diam:sc.log.Clutch & 0.0269317 & -0.2087418 & -0.0903139 &  & Correlational\\
\hline
Diff & sc.Diam:sc.sqrt.Pref & 0.1625125 & -0.1019349 & 0.0358978 &  & Correlational\\
\hline
Diff & sc.log.Clutch & -0.1546165 & -0.4409842 & -0.2851248 & * & Directional\\
\hline
Diff & sc.log.Clutch:sc.sqrt.Pref & 0.1704039 & -0.0787670 & 0.0470098 &  & Correlational\\
\hline
Diff & sc.sqrt.Pref & 0.1530003 & -0.1865995 & -0.0223829 &  & Directional\\
\hline
Simple & I(sc.Diam\textasciicircum{}2) & 0.2400693 & -0.0166148 & 0.0985395 &  & Quadratic\\
\hline
Simple & I(sc.log.Clutch\textasciicircum{}2) & 0.0258371 & -0.2706736 & -0.1156946 &  & Quadratic\\
\hline
Simple & I(sc.sqrt.Pref\textasciicircum{}2) & 0.1461630 & -0.2090979 & -0.0291386 &  & Quadratic\\
\hline
Simple & sc.Diam & 0.3179581 & 0.1211909 & 0.2130123 & * & Directional\\
\hline
Simple & sc.Diam:sc.log.Clutch & 0.0051116 & -0.1627098 & -0.0709867 &  & Correlational\\
\hline
Simple & sc.Diam:sc.sqrt.Pref & 0.0565334 & -0.0927118 & -0.0180742 &  & Correlational\\
\hline
Simple & sc.log.Clutch & -0.0580339 & -0.2639107 & -0.1577182 & * & Directional\\
\hline
Simple & sc.log.Clutch:sc.sqrt.Pref & 0.0780279 & -0.0705612 & 0.0071770 &  & Correlational\\
\hline
Simple & sc.sqrt.Pref & -0.0835763 & -0.3197919 & -0.1982366 & * & Directional\\
\hline
\end{tabular}

![](manuscript_files/figure-latex/Plot egg vs. larval ptoid gradients-1.pdf)<!-- --> 


We combine our estimates of selection gradients for each food-web treatment as well as the contribution of larval parasitoids to selection in the complex food web.

![**Selection gradients**. Estimates of standardized selection gradients in complex (orange) and simple (blue) food webs. The contribution of larval parasitoids (yellow) was estimated with a subset of the complex food-web data that only contained attacks from larval parasitoids (and gall survival). Points and lines correspond to estimates of the mean and 95% confidence intervals, respectively. Overlapping confidence intervals with zero (dotted line) indicate no strong evidence of selection.](manuscript_files/figure-latex/Figure-Selection-Gradients-1.pdf) 

# Partitioning the components of selection gradients

Selection gradients are influenced by both trait-fitness relationships and population mean fitness. Here, we partition selection gradients into these underlying components.

![**Partitioning Selection Gradients.**](manuscript_files/figure-latex/Plot foodwebegglarval coefficients-1.pdf) 

# Estimating selection on the egg parasitoid *Platygaster*


```r
# convert "gall_survival" to egg parasitoid survival. Note that both Iteomyia pupa and larva parasitoids result in 0. We do not change the name of the response variable to make clear that we are using the same statistical models as for Iteomyia.
eggegg_df <- mutate(gall_selection.df, gall_survival = ifelse(egg_parasitoid==1,1,0))

eggegg_model <- update(foodweb_model, data=eggegg_df)
```

```
## singular fit
```

```r
biased_eggegg_df <- eggegg_df %>%
  group_by(Foodweb, Gall_Number) %>%
  mutate(mean_survival = mean(gall_survival)) %>%
  filter(mean_survival > 0, mean_survival < 1) %>%
  ungroup()

biased_eggegg_model <- update(biased_foodweb_model, data=biased_eggegg_df)
```

```
## singular fit
```




```r
linear_eggegg_model <- update(linear_foodweb_model, data=eggegg_df)
```

```
## singular fit
```






\begin{tabular}{l|l|r|r|r|l|l}
\hline
type & term & conf.high & conf.low & mean & P\_cutoff & coefficient\_type\\
\hline
Complex & (Intercept) & 0.3114757 & 0.0237400 & 0.1598366 & * & Mean fitness\\
\hline
Complex & I(sc.Diam\textasciicircum{}2) & 0.1365455 & -0.6841817 & -0.2482329 &  & Quadratic\\
\hline
Complex & I(sc.log.Clutch\textasciicircum{}2) & 0.4782877 & -0.5311894 & -0.0138701 &  & Quadratic\\
\hline
Complex & I(sc.sqrt.Pref\textasciicircum{}2) & -0.1237055 & -1.2265495 & -0.6133299 & * & Quadratic\\
\hline
Complex & sc.Diam & -0.6153271 & -1.8173164 & -1.1716018 & * & Linear\\
\hline
Complex & sc.Diam:sc.log.Clutch & 0.8201621 & -0.2571678 & 0.2716512 &  & Correlational\\
\hline
Complex & sc.Diam:sc.sqrt.Pref & 1.3696159 & -0.0217621 & 0.6250497 &  & Correlational\\
\hline
Complex & sc.log.Clutch & 1.4816537 & 0.2205095 & 0.8282898 & * & Linear\\
\hline
Complex & sc.log.Clutch:sc.sqrt.Pref & 0.2966796 & -1.1104209 & -0.4015418 &  & Correlational\\
\hline
Complex & sc.sqrt.Pref & 0.9538929 & -0.4305641 & 0.2832528 &  & Linear\\
\hline
Diff & (Intercept) & 0.3056092 & -0.0886631 & 0.1010864 &  & Mean fitness\\
\hline
Diff & I(sc.Diam\textasciicircum{}2) & 0.6329564 & -0.6818572 & -0.0352761 &  & Quadratic\\
\hline
Diff & I(sc.log.Clutch\textasciicircum{}2) & 1.1041264 & -0.3030515 & 0.3726176 &  & Quadratic\\
\hline
Diff & I(sc.sqrt.Pref\textasciicircum{}2) & 1.3549567 & -0.1352078 & 0.5565992 &  & Quadratic\\
\hline
Diff & sc.Diam & 0.3849687 & -1.5093489 & -0.3962667 &  & Linear\\
\hline
Diff & sc.Diam:sc.log.Clutch & 0.9223406 & -0.7344833 & 0.1057713 &  & Correlational\\
\hline
Diff & sc.Diam:sc.sqrt.Pref & 0.1710301 & -1.6409986 & -0.6422392 &  & Correlational\\
\hline
Diff & sc.log.Clutch & 1.0331306 & -0.8006272 & 0.0938289 &  & Linear\\
\hline
Diff & sc.log.Clutch:sc.sqrt.Pref & 1.3910187 & -0.3889056 & 0.4687916 &  & Correlational\\
\hline
Diff & sc.sqrt.Pref & 2.2007100 & 0.1223735 & 1.0266191 & * & Linear\\
\hline
Simple & (Intercept) & 0.4542118 & 0.0732664 & 0.2609230 & * & Mean fitness\\
\hline
Simple & I(sc.Diam\textasciicircum{}2) & 0.1446221 & -0.7850159 & -0.2835090 &  & Quadratic\\
\hline
Simple & I(sc.log.Clutch\textasciicircum{}2) & 0.9199272 & -0.1590539 & 0.3587474 &  & Quadratic\\
\hline
Simple & I(sc.sqrt.Pref\textasciicircum{}2) & 0.4445187 & -0.5242999 & -0.0567306 &  & Quadratic\\
\hline
Simple & sc.Diam & -0.8860579 & -2.6684745 & -1.5678685 & * & Linear\\
\hline
Simple & sc.Diam:sc.log.Clutch & 0.9657093 & -0.1874869 & 0.3774225 &  & Correlational\\
\hline
Simple & sc.Diam:sc.sqrt.Pref & 0.5161961 & -0.5680212 & -0.0171895 &  & Correlational\\
\hline
Simple & sc.log.Clutch & 1.7542730 & 0.2580699 & 0.9221187 & * & Linear\\
\hline
Simple & sc.log.Clutch:sc.sqrt.Pref & 0.6445268 & -0.4035012 & 0.0672498 &  & Correlational\\
\hline
Simple & sc.sqrt.Pref & 2.2424971 & 0.6427081 & 1.3098720 & * & Linear\\
\hline
\end{tabular}


![](manuscript_files/figure-latex/Plot egg vs. egg ptoid coefficients-1.pdf)<!-- --> 




\begin{tabular}{l|l|r|r|r|l|l}
\hline
type & term & conf.high & conf.low & mean & P\_cutoff & gradient\_type\\
\hline
Complex & I(sc.Diam\textasciicircum{}2) & 0.1056635 & -0.5294428 & -0.1920909 &  & Quadratic\\
\hline
Complex & I(sc.log.Clutch\textasciicircum{}2) & 0.3701151 & -0.4110522 & -0.0107332 &  & Quadratic\\
\hline
Complex & I(sc.sqrt.Pref\textasciicircum{}2) & -0.0957274 & -0.9491452 & -0.4746153 & * & Quadratic\\
\hline
Complex & sc.Diam & -0.2380804 & -0.7031502 & -0.4533124 & * & Directional\\
\hline
Complex & sc.Diam:sc.log.Clutch & 0.3173345 & -0.0995025 & 0.1051064 &  & Correlational\\
\hline
Complex & sc.Diam:sc.sqrt.Pref & 0.5299274 & -0.0084201 & 0.2418422 &  & Correlational\\
\hline
Complex & sc.log.Clutch & 0.5732767 & 0.0853188 & 0.3204792 & * & Directional\\
\hline
Complex & sc.log.Clutch:sc.sqrt.Pref & 0.1147903 & -0.4296405 & -0.1553632 &  & Correlational\\
\hline
Complex & sc.sqrt.Pref & 0.3690772 & -0.1665925 & 0.1095953 &  & Directional\\
\hline
Diff & I(sc.Diam\textasciicircum{}2) & 0.4846209 & -0.3785137 & 0.0369696 &  & Quadratic\\
\hline
Diff & I(sc.log.Clutch\textasciicircum{}2) & 0.6965027 & -0.2558297 & 0.2070211 &  & Quadratic\\
\hline
Diff & I(sc.sqrt.Pref\textasciicircum{}2) & 1.0062894 & -0.0140187 & 0.4435752 &  & Quadratic\\
\hline
Diff & sc.Diam & 0.2771725 & -0.2895892 & 0.0243844 &  & Directional\\
\hline
Diff & sc.Diam:sc.log.Clutch & 0.2604027 & -0.2779800 & -0.0018534 &  & Correlational\\
\hline
Diff & sc.Diam:sc.sqrt.Pref & 0.0402846 & -0.5918523 & -0.2465448 &  & Correlational\\
\hline
Diff & sc.log.Clutch & 0.2310404 & -0.3561994 & -0.0682116 &  & Directional\\
\hline
Diff & sc.log.Clutch:sc.sqrt.Pref & 0.4954332 & -0.1263209 & 0.1737610 &  & Correlational\\
\hline
Diff & sc.sqrt.Pref & 0.6382296 & -0.0594115 & 0.2487516 &  & Directional\\
\hline
Simple & I(sc.Diam\textasciicircum{}2) & 0.0791297 & -0.4295199 & -0.1551213 &  & Quadratic\\
\hline
Simple & I(sc.log.Clutch\textasciicircum{}2) & 0.5033363 & -0.0870260 & 0.1962879 &  & Quadratic\\
\hline
Simple & I(sc.sqrt.Pref\textasciicircum{}2) & 0.2432175 & -0.2868696 & -0.0310401 &  & Quadratic\\
\hline
Simple & sc.Diam & -0.2424024 & -0.7300252 & -0.4289280 & * & Directional\\
\hline
Simple & sc.Diam:sc.log.Clutch & 0.2641929 & -0.0512915 & 0.1032530 &  & Correlational\\
\hline
Simple & sc.Diam:sc.sqrt.Pref & 0.1412178 & -0.1553958 & -0.0047026 &  & Correlational\\
\hline
Simple & sc.log.Clutch & 0.4799234 & 0.0706012 & 0.2522677 & * & Directional\\
\hline
Simple & sc.log.Clutch:sc.sqrt.Pref & 0.1763258 & -0.1103874 & 0.0183978 &  & Correlational\\
\hline
Simple & sc.sqrt.Pref & 0.6134889 & 0.1758282 & 0.3583469 & * & Directional\\
\hline
\end{tabular}


![](manuscript_files/figure-latex/Plot egg vs. egg ptoid gradients-1.pdf)<!-- --> 




# Reproduce Figure 4 - Univariate adaptive landscape










![Selection gradients acting on gall traits in complex vs. simple food webs. Each panel corresponds to a different gall trait: gall diameter (A); clutch size (B); and Oviposition preference (C). Solid lines represent the estimated gradients in complex (orange) and simple (blue) food webs. Transparent lines represent bootstrapped replicates (n=100) to show the uncertainty in estimated gradients. Note that only 100 bootstraps are displayed here, but that inferences are based on 1,000 bootstrapped samples.](manuscript_files/figure-latex/Univariate-Fitness-Landscapes-1.pdf) 





![](manuscript_files/figure-latex/Egg-Parasitoid-Selection-1.pdf)<!-- --> 

![](manuscript_files/figure-latex/Egg-Parasitoid-Selection-linear-1.pdf)<!-- --> 





# Multivariate fitness landscapes









Reproduce Figure 3 in main manuscript.

![Fitness landscapes of gall traits in complex vs. simple food webs. Each panel corresponds to a different combination of traits: clutch size and gall diameter (A); clutch size and Oviposition preference (B); Oviposition preference and gall diameter (C). Note that traits for all plots range 1 SD below and above the mean (=0).](manuscript_files/figure-latex/Figure-Multivariate-Landscapes-1.pdf) 



