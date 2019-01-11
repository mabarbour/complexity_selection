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
Complex & (Intercept) & 0.5577970 & 0.2746507 & 0.4132922 & * & Mean fitness\\
\hline
Complex & I(sc.Diam\textasciicircum{}2) & 0.5383929 & -0.0950983 & 0.2271685 &  & Quadratic\\
\hline
Complex & I(sc.log.Clutch\textasciicircum{}2) & 0.2604850 & -0.4216912 & -0.0884475 &  & Quadratic\\
\hline
Complex & I(sc.sqrt.Pref\textasciicircum{}2) & 1.0744532 & 0.0999931 & 0.5566625 & * & Quadratic\\
\hline
Complex & sc.Diam & 1.5693700 & 0.7298912 & 1.1352672 & * & Linear\\
\hline
Complex & sc.Diam:sc.log.Clutch & 0.2525351 & -0.5715558 & -0.1577118 &  & Correlational\\
\hline
Complex & sc.Diam:sc.sqrt.Pref & 0.0059841 & -0.9494477 & -0.4424437 &  & Correlational\\
\hline
Complex & sc.log.Clutch & 0.5748866 & -0.1309147 & 0.2102483 &  & Linear\\
\hline
Complex & sc.log.Clutch:sc.sqrt.Pref & 0.5866006 & -0.3682643 & 0.0970880 &  & Correlational\\
\hline
Complex & sc.sqrt.Pref & 0.1548968 & -0.9673799 & -0.4054556 &  & Linear\\
\hline
Diff & (Intercept) & 0.4174345 & 0.1187386 & 0.2624413 & * & Mean fitness\\
\hline
Diff & I(sc.Diam\textasciicircum{}2) & 0.5212381 & -0.4324914 & 0.0397940 &  & Quadratic\\
\hline
Diff & I(sc.log.Clutch\textasciicircum{}2) & 0.3120400 & -0.7667277 & -0.2179471 &  & Quadratic\\
\hline
Diff & I(sc.sqrt.Pref\textasciicircum{}2) & 0.0941926 & -1.1762737 & -0.5220837 &  & Quadratic\\
\hline
Diff & sc.Diam & 0.5147400 & -0.5651894 & -0.0397714 &  & Linear\\
\hline
Diff & sc.Diam:sc.log.Clutch & 0.4400598 & -0.7849733 & -0.1938585 &  & Correlational\\
\hline
Diff & sc.Diam:sc.sqrt.Pref & 0.9724702 & -0.2653500 & 0.3573678 &  & Correlational\\
\hline
Diff & sc.log.Clutch & -0.1261059 & -1.2744060 & -0.6819489 & * & Linear\\
\hline
Diff & sc.log.Clutch:sc.sqrt.Pref & 0.4876757 & -0.6934578 & -0.0965558 &  & Correlational\\
\hline
Diff & sc.sqrt.Pref & 0.2549453 & -1.1550758 & -0.4321609 &  & Linear\\
\hline
Simple & (Intercept) & 0.8082259 & 0.5313324 & 0.6757334 & * & Mean fitness\\
\hline
Simple & I(sc.Diam\textasciicircum{}2) & 0.6535475 & -0.0680735 & 0.2669625 &  & Quadratic\\
\hline
Simple & I(sc.log.Clutch\textasciicircum{}2) & 0.0801733 & -0.7093998 & -0.3063946 &  & Quadratic\\
\hline
Simple & I(sc.sqrt.Pref\textasciicircum{}2) & 0.4628712 & -0.3810782 & 0.0345788 &  & Quadratic\\
\hline
Simple & sc.Diam & 1.5833488 & 0.6759297 & 1.0954958 & * & Linear\\
\hline
Simple & sc.Diam:sc.log.Clutch & 0.0807693 & -0.7782059 & -0.3515703 &  & Correlational\\
\hline
Simple & sc.Diam:sc.sqrt.Pref & 0.2986535 & -0.4957363 & -0.0850759 &  & Correlational\\
\hline
Simple & sc.log.Clutch & -0.0892229 & -0.9376601 & -0.4717006 & * & Linear\\
\hline
Simple & sc.log.Clutch:sc.sqrt.Pref & 0.3804199 & -0.3657039 & 0.0005322 &  & Correlational\\
\hline
Simple & sc.sqrt.Pref & -0.3622454 & -1.3513666 & -0.8376166 & * & Linear\\
\hline
\end{tabular}


![](manuscript_files/figure-latex/Plot food-web coefficients-1.pdf)<!-- --> 




\begin{tabular}{l|l|r|r|r|l|l}
\hline
type & term & conf.high & conf.low & mean & P\_cutoff & gradient\_type\\
\hline
Complex & I(sc.Diam\textasciicircum{}2) & 0.3226752 & -0.0569953 & 0.1361490 &  & Quadratic\\
\hline
Complex & I(sc.log.Clutch\textasciicircum{}2) & 0.1561166 & -0.2527324 & -0.0530093 &  & Quadratic\\
\hline
Complex & I(sc.sqrt.Pref\textasciicircum{}2) & 0.6439525 & 0.0599289 & 0.3336248 & * & Quadratic\\
\hline
Complex & sc.Diam & 0.4702856 & 0.2187230 & 0.3402001 & * & Directional\\
\hline
Complex & sc.Diam:sc.log.Clutch & 0.0756760 & -0.1712754 & -0.0472607 &  & Correlational\\
\hline
Complex & sc.Diam:sc.sqrt.Pref & 0.0017932 & -0.2845164 & -0.1325850 &  & Correlational\\
\hline
Complex & sc.log.Clutch & 0.1722735 & -0.0392306 & 0.0630041 &  & Directional\\
\hline
Complex & sc.log.Clutch:sc.sqrt.Pref & 0.1757838 & -0.1103560 & 0.0290939 &  & Correlational\\
\hline
Complex & sc.sqrt.Pref & 0.0464172 & -0.2898901 & -0.1215010 &  & Directional\\
\hline
Diff & I(sc.Diam\textasciicircum{}2) & 0.1971101 & -0.2716844 & -0.0357266 &  & Quadratic\\
\hline
Diff & I(sc.log.Clutch\textasciicircum{}2) & 0.1884530 & -0.3311794 & -0.0622461 &  & Quadratic\\
\hline
Diff & I(sc.sqrt.Pref\textasciicircum{}2) & -0.0035932 & -0.6628308 & -0.3206174 & * & Quadratic\\
\hline
Diff & sc.Diam & 0.0064375 & -0.2627085 & -0.1341556 &  & Directional\\
\hline
Diff & sc.Diam:sc.log.Clutch & 0.1305590 & -0.1583110 & -0.0188638 &  & Correlational\\
\hline
Diff & sc.Diam:sc.sqrt.Pref & 0.2850260 & -0.0457236 & 0.1165836 &  & Correlational\\
\hline
Diff & sc.log.Clutch & -0.0159345 & -0.2928373 & -0.1517231 & * & Directional\\
\hline
Diff & sc.log.Clutch:sc.sqrt.Pref & 0.1228931 & -0.1865889 & -0.0289938 &  & Correlational\\
\hline
Diff & sc.sqrt.Pref & 0.1387417 & -0.2201026 & -0.0360408 &  & Directional\\
\hline
Simple & I(sc.Diam\textasciicircum{}2) & 0.2458428 & -0.0256070 & 0.1004224 &  & Quadratic\\
\hline
Simple & I(sc.log.Clutch\textasciicircum{}2) & 0.0301585 & -0.2668526 & -0.1152554 &  & Quadratic\\
\hline
Simple & I(sc.sqrt.Pref\textasciicircum{}2) & 0.1741167 & -0.1433489 & 0.0130074 &  & Quadratic\\
\hline
Simple & sc.Diam & 0.2978015 & 0.1271311 & 0.2060445 & * & Directional\\
\hline
Simple & sc.Diam:sc.log.Clutch & 0.0151914 & -0.1463676 & -0.0661245 &  & Correlational\\
\hline
Simple & sc.Diam:sc.sqrt.Pref & 0.0561717 & -0.0932397 & -0.0160014 &  & Correlational\\
\hline
Simple & sc.log.Clutch & -0.0167814 & -0.1763583 & -0.0887190 & * & Directional\\
\hline
Simple & sc.log.Clutch:sc.sqrt.Pref & 0.0715506 & -0.0687828 & 0.0001001 &  & Correlational\\
\hline
Simple & sc.sqrt.Pref & -0.0681323 & -0.2541695 & -0.1575417 & * & Directional\\
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




```r
linear_egglarval_model <- update(linear_foodweb_model, data=egglarval_df)
```






\begin{tabular}{l|l|r|r|r|l|l}
\hline
type & term & conf.high & conf.low & mean & P\_cutoff & coefficient\_type\\
\hline
Complex & (Intercept) & 0.8291572 & 0.5018696 & 0.6707670 & * & Mean fitness\\
\hline
Complex & I(sc.Diam\textasciicircum{}2) & 0.5252164 & -0.2124211 & 0.1546565 &  & Quadratic\\
\hline
Complex & I(sc.log.Clutch\textasciicircum{}2) & 0.3763112 & -0.3412520 & 0.0192410 &  & Quadratic\\
\hline
Complex & I(sc.sqrt.Pref\textasciicircum{}2) & 1.0809739 & -0.0286449 & 0.4635128 &  & Quadratic\\
\hline
Complex & sc.Diam & 1.7505608 & 0.6409519 & 1.1493273 & * & Linear\\
\hline
Complex & sc.Diam:sc.log.Clutch & 0.5704469 & -0.3539789 & 0.1020909 &  & Correlational\\
\hline
Complex & sc.Diam:sc.sqrt.Pref & 0.3178538 & -0.8451925 & -0.2724925 &  & Correlational\\
\hline
Complex & sc.log.Clutch & 1.1814515 & 0.1596120 & 0.6742001 & * & Linear\\
\hline
Complex & sc.log.Clutch:sc.sqrt.Pref & 0.2750507 & -0.7387619 & -0.2048680 &  & Correlational\\
\hline
Complex & sc.sqrt.Pref & -0.1391448 & -1.7356777 & -0.8941844 & * & Linear\\
\hline
Diff & (Intercept) & 0.2196657 & -0.1799276 & 0.0174850 &  & Mean fitness\\
\hline
Diff & I(sc.Diam\textasciicircum{}2) & 0.5574665 & -0.3584458 & 0.1022658 &  & Quadratic\\
\hline
Diff & I(sc.log.Clutch\textasciicircum{}2) & 0.1860787 & -0.8691844 & -0.3225297 &  & Quadratic\\
\hline
Diff & I(sc.sqrt.Pref\textasciicircum{}2) & 0.1316942 & -1.3477486 & -0.5351730 &  & Quadratic\\
\hline
Diff & sc.Diam & 0.5371689 & -0.7334456 & -0.0698901 &  & Linear\\
\hline
Diff & sc.Diam:sc.log.Clutch & 0.1098853 & -1.1099933 & -0.4700208 &  & Correlational\\
\hline
Diff & sc.Diam:sc.sqrt.Pref & 0.8535288 & -0.4724468 & 0.1677916 &  & Correlational\\
\hline
Diff & sc.log.Clutch & -0.7742105 & -2.3326708 & -1.4870077 & * & Linear\\
\hline
Diff & sc.log.Clutch:sc.sqrt.Pref & 0.8519946 & -0.3230277 & 0.2409038 &  & Correlational\\
\hline
Diff & sc.sqrt.Pref & 0.8048892 & -0.9803708 & -0.1097636 &  & Linear\\
\hline
Simple & (Intercept) & 0.8534582 & 0.5064552 & 0.6882520 & * & Mean fitness\\
\hline
Simple & I(sc.Diam\textasciicircum{}2) & 0.6009394 & -0.0599399 & 0.2569223 &  & Quadratic\\
\hline
Simple & I(sc.log.Clutch\textasciicircum{}2) & 0.0426746 & -0.6945156 & -0.3032887 &  & Quadratic\\
\hline
Simple & I(sc.sqrt.Pref\textasciicircum{}2) & 0.3458282 & -0.5386174 & -0.0716601 &  & Quadratic\\
\hline
Simple & sc.Diam & 1.6309519 & 0.6443025 & 1.0794372 & * & Linear\\
\hline
Simple & sc.Diam:sc.log.Clutch & 0.0521287 & -0.8171435 & -0.3679299 &  & Correlational\\
\hline
Simple & sc.Diam:sc.sqrt.Pref & 0.2953305 & -0.4636802 & -0.1047010 &  & Correlational\\
\hline
Simple & sc.log.Clutch & -0.3237684 & -1.3682319 & -0.8128076 & * & Linear\\
\hline
Simple & sc.log.Clutch:sc.sqrt.Pref & 0.4021321 & -0.3141034 & 0.0360358 &  & Correlational\\
\hline
Simple & sc.sqrt.Pref & -0.4374714 & -1.6020535 & -1.0039480 & * & Linear\\
\hline
\end{tabular}


![](manuscript_files/figure-latex/Plot egglarval coefficients-1.pdf)<!-- --> 




\begin{tabular}{l|l|r|r|r|l|l}
\hline
type & term & conf.high & conf.low & mean & P\_cutoff & gradient\_type\\
\hline
Complex & I(sc.Diam\textasciicircum{}2) & 0.2034272 & -0.0822751 & 0.0599017 &  & Quadratic\\
\hline
Complex & I(sc.log.Clutch\textasciicircum{}2) & 0.1457531 & -0.1321740 & 0.0074525 &  & Quadratic\\
\hline
Complex & I(sc.sqrt.Pref\textasciicircum{}2) & 0.4186836 & -0.0110947 & 0.1795281 &  & Quadratic\\
\hline
Complex & sc.Diam & 0.3390142 & 0.1241270 & 0.2225791 & * & Directional\\
\hline
Complex & sc.Diam:sc.log.Clutch & 0.1104729 & -0.0685517 & 0.0197710 &  & Correlational\\
\hline
Complex & sc.Diam:sc.sqrt.Pref & 0.0615557 & -0.1636803 & -0.0527710 &  & Correlational\\
\hline
Complex & sc.log.Clutch & 0.2288003 & 0.0309105 & 0.1305658 & * & Directional\\
\hline
Complex & sc.log.Clutch:sc.sqrt.Pref & 0.0532664 & -0.1430689 & -0.0396748 &  & Correlational\\
\hline
Complex & sc.sqrt.Pref & -0.0269468 & -0.3361320 & -0.1731681 & * & Directional\\
\hline
Diff & I(sc.Diam\textasciicircum{}2) & 0.2187304 & -0.1383949 & 0.0417561 &  & Quadratic\\
\hline
Diff & I(sc.log.Clutch\textasciicircum{}2) & 0.0720939 & -0.3400311 & -0.1274562 &  & Quadratic\\
\hline
Diff & I(sc.sqrt.Pref\textasciicircum{}2) & 0.0532460 & -0.5243575 & -0.2078822 &  & Quadratic\\
\hline
Diff & sc.Diam & 0.1107413 & -0.1398044 & -0.0090260 &  & Directional\\
\hline
Diff & sc.Diam:sc.log.Clutch & 0.0204044 & -0.2177643 & -0.0925613 &  & Correlational\\
\hline
Diff & sc.Diam:sc.sqrt.Pref & 0.1660019 & -0.0930250 & 0.0320572 &  & Correlational\\
\hline
Diff & sc.log.Clutch & -0.1512996 & -0.4584534 & -0.2913696 & * & Directional\\
\hline
Diff & sc.log.Clutch:sc.sqrt.Pref & 0.1656268 & -0.0629836 & 0.0468040 &  & Correlational\\
\hline
Diff & sc.sqrt.Pref & 0.1524155 & -0.1934690 & -0.0254505 &  & Directional\\
\hline
Simple & I(sc.Diam\textasciicircum{}2) & 0.2377767 & -0.0237167 & 0.1016577 &  & Quadratic\\
\hline
Simple & I(sc.log.Clutch\textasciicircum{}2) & 0.0168853 & -0.2748024 & -0.1200037 &  & Quadratic\\
\hline
Simple & I(sc.sqrt.Pref\textasciicircum{}2) & 0.1368356 & -0.2131174 & -0.0283541 &  & Quadratic\\
\hline
Simple & sc.Diam & 0.3226634 & 0.1274672 & 0.2135531 & * & Directional\\
\hline
Simple & sc.Diam:sc.log.Clutch & 0.0103130 & -0.1616616 & -0.0727903 &  & Correlational\\
\hline
Simple & sc.Diam:sc.sqrt.Pref & 0.0584274 & -0.0917333 & -0.0207138 &  & Correlational\\
\hline
Simple & sc.log.Clutch & -0.0640535 & -0.2706876 & -0.1608038 & * & Directional\\
\hline
Simple & sc.log.Clutch:sc.sqrt.Pref & 0.0795568 & -0.0621414 & 0.0071292 &  & Correlational\\
\hline
Simple & sc.sqrt.Pref & -0.0865482 & -0.3169462 & -0.1986185 & * & Directional\\
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






\begin{tabular}{l|l|r|r|r|l|l}
\hline
type & term & conf.high & conf.low & mean & P\_cutoff & coefficient\_type\\
\hline
Complex & (Intercept) & 0.3093470 & 0.0342100 & 0.1637632 & * & Mean fitness\\
\hline
Complex & I(sc.Diam\textasciicircum{}2) & 0.1134587 & -0.6750140 & -0.2620791 &  & Quadratic\\
\hline
Complex & I(sc.log.Clutch\textasciicircum{}2) & 0.4333846 & -0.5120566 & -0.0250286 &  & Quadratic\\
\hline
Complex & I(sc.sqrt.Pref\textasciicircum{}2) & -0.1265881 & -1.3942082 & -0.6392281 & * & Quadratic\\
\hline
Complex & sc.Diam & 0.2197196 & -0.9119897 & -0.3078570 &  & Linear\\
\hline
Complex & sc.Diam:sc.log.Clutch & 0.8818965 & -0.2424917 & 0.2790402 &  & Correlational\\
\hline
Complex & sc.Diam:sc.sqrt.Pref & 1.2896701 & 0.0594661 & 0.6355660 & * & Correlational\\
\hline
Complex & sc.log.Clutch & 1.5957475 & 0.2098793 & 0.8263345 & * & Linear\\
\hline
Complex & sc.log.Clutch:sc.sqrt.Pref & 0.2519279 & -1.1311451 & -0.4264012 &  & Correlational\\
\hline
Complex & sc.sqrt.Pref & 0.9876934 & -0.3733873 & 0.2876438 &  & Linear\\
\hline
Diff & (Intercept) & 0.3077723 & -0.0794831 & 0.0980581 &  & Mean fitness\\
\hline
Diff & I(sc.Diam\textasciicircum{}2) & 0.5272038 & -0.6777435 & -0.0207512 &  & Quadratic\\
\hline
Diff & I(sc.log.Clutch\textasciicircum{}2) & 1.1007484 & -0.3052390 & 0.3783894 &  & Quadratic\\
\hline
Diff & I(sc.sqrt.Pref\textasciicircum{}2) & 1.4343500 & -0.1590111 & 0.5703844 &  & Quadratic\\
\hline
Diff & sc.Diam & 0.9350398 & -1.0200247 & 0.0747079 &  & Linear\\
\hline
Diff & sc.Diam:sc.log.Clutch & 0.9744818 & -0.7101728 & 0.1068076 &  & Correlational\\
\hline
Diff & sc.Diam:sc.sqrt.Pref & 0.2253183 & -1.5866545 & -0.6469685 &  & Correlational\\
\hline
Diff & sc.log.Clutch & 1.1184565 & -0.9005062 & 0.0917723 &  & Linear\\
\hline
Diff & sc.log.Clutch:sc.sqrt.Pref & 1.4033743 & -0.3594328 & 0.4985444 &  & Correlational\\
\hline
Diff & sc.sqrt.Pref & 2.2216253 & 0.0253542 & 1.0040638 & * & Linear\\
\hline
Simple & (Intercept) & 0.4615029 & 0.0835218 & 0.2618213 & * & Mean fitness\\
\hline
Simple & I(sc.Diam\textasciicircum{}2) & 0.1457836 & -0.7895109 & -0.2828304 &  & Quadratic\\
\hline
Simple & I(sc.log.Clutch\textasciicircum{}2) & 1.0206323 & -0.1132328 & 0.3533609 &  & Quadratic\\
\hline
Simple & I(sc.sqrt.Pref\textasciicircum{}2) & 0.3531123 & -0.5615823 & -0.0688438 &  & Quadratic\\
\hline
Simple & sc.Diam & 0.4583775 & -1.4036169 & -0.2331491 &  & Linear\\
\hline
Simple & sc.Diam:sc.log.Clutch & 0.9969600 & -0.1532823 & 0.3858478 &  & Correlational\\
\hline
Simple & sc.Diam:sc.sqrt.Pref & 0.4805191 & -0.5843721 & -0.0114025 &  & Correlational\\
\hline
Simple & sc.log.Clutch & 1.8539833 & 0.2214306 & 0.9181068 & * & Linear\\
\hline
Simple & sc.log.Clutch:sc.sqrt.Pref & 0.6062398 & -0.4403961 & 0.0721432 &  & Correlational\\
\hline
Simple & sc.sqrt.Pref & 2.1299312 & 0.6847067 & 1.2917076 & * & Linear\\
\hline
\end{tabular}


![](manuscript_files/figure-latex/Plot egg vs. egg ptoid coefficients-1.pdf)<!-- --> 




\begin{tabular}{l|l|r|r|r|l|l}
\hline
type & term & conf.high & conf.low & mean & P\_cutoff & gradient\_type\\
\hline
Complex & I(sc.Diam\textasciicircum{}2) & 0.0877982 & -0.5223485 & -0.2028056 &  & Quadratic\\
\hline
Complex & I(sc.log.Clutch\textasciicircum{}2) & 0.3353676 & -0.3962466 & -0.0193679 &  & Quadratic\\
\hline
Complex & I(sc.sqrt.Pref\textasciicircum{}2) & -0.0979581 & -1.0788852 & -0.4946562 & * & Quadratic\\
\hline
Complex & sc.Diam & 0.0850132 & -0.3528641 & -0.1191151 &  & Directional\\
\hline
Complex & sc.Diam:sc.log.Clutch & 0.3412206 & -0.0938241 & 0.1079653 &  & Correlational\\
\hline
Complex & sc.Diam:sc.sqrt.Pref & 0.4989950 & 0.0230084 & 0.2459112 & * & Correlational\\
\hline
Complex & sc.log.Clutch & 0.6174215 & 0.0812058 & 0.3197227 & * & Directional\\
\hline
Complex & sc.log.Clutch:sc.sqrt.Pref & 0.0974751 & -0.4376590 & -0.1649818 &  & Correlational\\
\hline
Complex & sc.sqrt.Pref & 0.3821552 & -0.1444698 & 0.1112942 &  & Directional\\
\hline
Diff & I(sc.Diam\textasciicircum{}2) & 0.4181186 & -0.3701468 & 0.0480556 &  & Quadratic\\
\hline
Diff & I(sc.log.Clutch\textasciicircum{}2) & 0.6783308 & -0.2392577 & 0.2127086 &  & Quadratic\\
\hline
Diff & I(sc.sqrt.Pref\textasciicircum{}2) & 1.0956040 & -0.0339619 & 0.4569885 &  & Quadratic\\
\hline
Diff & sc.Diam & 0.3372742 & -0.2578536 & 0.0553315 &  & Directional\\
\hline
Diff & sc.Diam:sc.log.Clutch & 0.2773973 & -0.2705193 & -0.0024074 &  & Correlational\\
\hline
Diff & sc.Diam:sc.sqrt.Pref & 0.0438126 & -0.5659037 & -0.2490306 &  & Correlational\\
\hline
Diff & sc.log.Clutch & 0.2643048 & -0.4093538 & -0.0685526 &  & Directional\\
\hline
Diff & sc.log.Clutch:sc.sqrt.Pref & 0.5003004 & -0.1039114 & 0.1847183 &  & Correlational\\
\hline
Diff & sc.sqrt.Pref & 0.6167639 & -0.0917153 & 0.2420834 &  & Directional\\
\hline
Simple & I(sc.Diam\textasciicircum{}2) & 0.0797652 & -0.4319793 & -0.1547501 &  & Quadratic\\
\hline
Simple & I(sc.log.Clutch\textasciicircum{}2) & 0.5584369 & -0.0619551 & 0.1933407 &  & Quadratic\\
\hline
Simple & I(sc.sqrt.Pref\textasciicircum{}2) & 0.1932047 & -0.3072686 & -0.0376677 &  & Quadratic\\
\hline
Simple & sc.Diam & 0.1254002 & -0.3839931 & -0.0637835 &  & Directional\\
\hline
Simple & sc.Diam:sc.log.Clutch & 0.2727423 & -0.0419341 & 0.1055579 &  & Correlational\\
\hline
Simple & sc.Diam:sc.sqrt.Pref & 0.1314575 & -0.1598690 & -0.0031194 &  & Correlational\\
\hline
Simple & sc.log.Clutch & 0.5072016 & 0.0605777 & 0.2511701 & * & Directional\\
\hline
Simple & sc.log.Clutch:sc.sqrt.Pref & 0.1658514 & -0.1204809 & 0.0197365 &  & Correlational\\
\hline
Simple & sc.sqrt.Pref & 0.5826938 & 0.1873179 & 0.3533776 & * & Directional\\
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



