---
title: "Supplementary Material"
author: "Matthew A. Barbour"
date: "2018-11-14"
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

visreg::visreg(linear_foodweb_model, xvar="sc.sqrt.Pref", by="Foodweb", scale="response")
```

![](manuscript_files/figure-latex/Linear food-web model-1.pdf)<!-- --> 




To estimate biased selection on chamber diameter, we subset our data to only include multi-chambered galls where there was variability in larval survival. We then fit a reduced model to estimate the bias in the logistic regression coefficient of chamber diameter in each food web.

<!-- We quantified biased selection on chamber diameter in the following way. First, we focus on multi-chambered galls where there is evidence of both survival and parasitism. The idea is that galls from the same clutch should be similar in size; therefore, any selection on diameter would be apparent. If we assume that this is all due to the effect of parasitism on larval development, truncating gall diameter, then we can consider this new selection differential to be the result of confounding factors. Note that this method likely overestimates the confounding effect of parasitism, since we are assuming that any heterogeneity in chamber diameter (for this subset) is due to parasitism. -->


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

visreg::visreg(eggegg_model, xvar="sc.sqrt.Pref", by="sc.Diam", scale="response", cond = list(Foodweb="Simple"))
```

![](manuscript_files/figure-latex/Egg vs. egg ptoid model-1.pdf)<!-- --> 

```r
visreg::visreg(eggegg_model, xvar="sc.sqrt.Pref", by="sc.Diam", scale="response", cond = list(Foodweb="Complex"))
```

![](manuscript_files/figure-latex/Egg vs. egg ptoid model-2.pdf)<!-- --> 

```r
visreg::visreg(eggegg_model, xvar="sc.Diam", by="sc.sqrt.Pref", scale="response", cond = list(Foodweb="Simple"))
```

![](manuscript_files/figure-latex/Egg vs. egg ptoid model-3.pdf)<!-- --> 

```r
visreg::visreg(eggegg_model, xvar="sc.Diam", by="sc.sqrt.Pref", scale="response", cond = list(Foodweb="Complex"))
```

![](manuscript_files/figure-latex/Egg vs. egg ptoid model-4.pdf)<!-- --> 




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
Complex & sc.Diam & -0.6522964 & -1.7840057 & -1.1798731 & * & Linear\\
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
Diff & sc.Diam & 0.4654369 & -1.4896276 & -0.3948950 &  & Linear\\
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
Simple & sc.Diam & -0.8832414 & -2.7452358 & -1.5747680 & * & Linear\\
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
Complex & sc.Diam & -0.2523844 & -0.6902618 & -0.4565127 & * & Directional\\
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
Diff & sc.Diam & 0.3076399 & -0.2874880 & 0.0256971 &  & Directional\\
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
Simple & sc.Diam & -0.2416319 & -0.7510251 & -0.4308156 & * & Directional\\
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












![](manuscript_files/figure-latex/Figure-Univariate-Landscapes-1.pdf)<!-- --> 



![](manuscript_files/figure-latex/Egg-Parasitoid-Selection-1.pdf)<!-- --> 

![](manuscript_files/figure-latex/Egg-Parasitoid-Selection-linear-1.pdf)<!-- --> 





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
