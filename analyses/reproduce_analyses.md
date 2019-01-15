Reproduce analyses reported in: *Phenotypic evolution is more constrained in simpler food webs*
================
Matthew A. Barbour
2019-01-15

Evaluating assumption of multivariate normality
===============================================

We used graphical checks to evaluate whether our transformations of trait values resulted in a multivariate normal distribution. Figure S shows that our transformations resulted in approximately normal distributions for each phenotypic trait. Note also that in the multivariate quantile-quantile (Q-Q) plot, most points fall along the expected line (fig. S), suggesting that our transformations provide a reasonable approximation of a multivariate normal distribution.

``` r
mvn_univariate_hist <- mvn(select(gall_selection.df, sc.Diam, sc.log.Clutch, sc.sqrt.Pref), univariatePlot = "histogram")
```

![Histograms of each phenotypic trait after transformation. The red line illustrates a normal distribution.](reproduce_analyses_files/figure-markdown_github/Univariate_histograms-1.pdf)

``` r
mvn_multivariate_qq <- mvn(select(gall_selection.df, sc.Diam, sc.log.Clutch, sc.sqrt.Pref), multivariatePlot = "qq")
```

![Multivariate quantile-quantile (Q-Q) plot to assess deviations from multivariate normality (black line).](reproduce_analyses_files/figure-markdown_github/Multivariate_QQ-1.pdf)

Effect of food-web treatment on trait-fitness relationships and selection gradients
===================================================================================

We write the model in a way that independently estimates the effect of food-web treatment, each trait, and all two-way and three-way statistical interactions, on larval survival. The resulting estimates and confidence intervals are useful for determining whether trait-fitness relationships differ from zero, but not whether they differ between food-web treatments. For the later, we calculate the differences between each food-web treatment from the bootstrapped samples.

``` r
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

``` r
boot_foodweb_model <- bootMer(foodweb_model, FUN = fixef, nsim=n_boots_analysis, parallel="multicore", ncpus=32, seed=34)
```

In order to reliably quantify linear trait-fitness relationships, we remove all higher-order terms from the model [Stinchcombe et al. 2008](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1558-5646.2008.00449.x).

``` r
linear_foodweb_model <- glmer(
  gall_survival ~ 
    -1 + Foodweb + 
    Foodweb:(sc.Diam + sc.log.Clutch + sc.sqrt.Pref) +
    (1|Genotype/Plant_Position/Gall_Number),
  data = gall_selection.df,
  family = binomial(link = logit), control=glmerControl(optimizer = "bobyqa"))
```

And refit using parametric bootstrapping.

``` r
boot_linear_foodweb_model <- bootMer(linear_foodweb_model, FUN = fixef, nsim=n_boots_analysis, parallel="multicore", ncpus=32, seed=34)
```

Chamber diameter, but not the other phenotypes, could be influenced by parasitism itself rather than being under natural selection. To estimate this potential bias, we subset our data to only include multi-chambered galls where there was variability in larval survival. We then fit a reduced model to estimate the bias in the logistic regression coefficient of chamber diameter in each food web.

``` r
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
```

    ## singular fit

``` r
# Accuracy of confidence intervals is not a priority here, so we just use the asymptotic estimates rather than parametric bootstrapping.
biased_foodweb_confint <- tidy(biased_foodweb_model, conf.int=TRUE) %>% filter(group=="fixed")
```

    ## Warning in bind_rows_(x, .id): binding factor and character vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

``` r
# Print table
knitr::kable(biased_foodweb_confint)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
term
</th>
<th style="text-align:right;">
estimate
</th>
<th style="text-align:right;">
std.error
</th>
<th style="text-align:right;">
statistic
</th>
<th style="text-align:right;">
p.value
</th>
<th style="text-align:right;">
conf.low
</th>
<th style="text-align:right;">
conf.high
</th>
<th style="text-align:left;">
group
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
FoodwebComplex
</td>
<td style="text-align:right;">
0.0261288
</td>
<td style="text-align:right;">
0.1596126
</td>
<td style="text-align:right;">
0.1637011
</td>
<td style="text-align:right;">
0.8699664
</td>
<td style="text-align:right;">
-0.2867063
</td>
<td style="text-align:right;">
0.3389638
</td>
<td style="text-align:left;">
fixed
</td>
</tr>
<tr>
<td style="text-align:left;">
FoodwebSimple
</td>
<td style="text-align:right;">
-0.1181035
</td>
<td style="text-align:right;">
0.2044004
</td>
<td style="text-align:right;">
-0.5778045
</td>
<td style="text-align:right;">
0.5633961
</td>
<td style="text-align:right;">
-0.5187209
</td>
<td style="text-align:right;">
0.2825139
</td>
<td style="text-align:left;">
fixed
</td>
</tr>
<tr>
<td style="text-align:left;">
FoodwebComplex:sc.Diam
</td>
<td style="text-align:right;">
0.3603278
</td>
<td style="text-align:right;">
0.1564728
</td>
<td style="text-align:right;">
2.3028133
</td>
<td style="text-align:right;">
0.0212893
</td>
<td style="text-align:right;">
0.0536466
</td>
<td style="text-align:right;">
0.6670089
</td>
<td style="text-align:left;">
fixed
</td>
</tr>
<tr>
<td style="text-align:left;">
FoodwebSimple:sc.Diam
</td>
<td style="text-align:right;">
0.4158534
</td>
<td style="text-align:right;">
0.2066780
</td>
<td style="text-align:right;">
2.0120841
</td>
<td style="text-align:right;">
0.0442111
</td>
<td style="text-align:right;">
0.0107721
</td>
<td style="text-align:right;">
0.8209348
</td>
<td style="text-align:left;">
fixed
</td>
</tr>
</tbody>
</table>
We gather regression coefficients from the full and reduced model to display trait-fitness relationships. For chamber diameter, we subtracted the coefficient due to bias to better approximate the trait-fitness relationship. We also used the bootstrapped samples to calculate the differences in estimates between treatments (labelled **Diff**). The code for tidying this data is not displayed, but is available in the 'reproduce\_analyses.Rmd'.

``` r
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

    ## Warning: Expected 3 pieces. Missing pieces filled with `NA` in 63 rows [1,
    ## 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 19, 20, 21, 22, 23, 24, 25, 26, ...].

Table of trait-fitness relationships:

``` r
knitr::kable(tidy_foodweb_alphas)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
type
</th>
<th style="text-align:left;">
term
</th>
<th style="text-align:right;">
conf.high
</th>
<th style="text-align:right;">
conf.low
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:left;">
P\_cutoff
</th>
<th style="text-align:left;">
coefficient\_type
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Complex
</td>
<td style="text-align:left;">
(Intercept)
</td>
<td style="text-align:right;">
0.5600735
</td>
<td style="text-align:right;">
0.2764114
</td>
<td style="text-align:right;">
0.4157445
</td>
<td style="text-align:left;">
-   </td>
    <td style="text-align:left;">
    Mean fitness
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    I(sc.Diam^2)
    </td>
    <td style="text-align:right;">
    0.5534963
    </td>
    <td style="text-align:right;">
    -0.1078761
    </td>
    <td style="text-align:right;">
    0.2225242
    </td>
    <td style="text-align:left;">
    </td>
    <td style="text-align:left;">
    Quadratic
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    I(sc.log.Clutch^2)
    </td>
    <td style="text-align:right;">
    0.2974010
    </td>
    <td style="text-align:right;">
    -0.4471858
    </td>
    <td style="text-align:right;">
    -0.0902569
    </td>
    <td style="text-align:left;">
    </td>
    <td style="text-align:left;">
    Quadratic
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    I(sc.sqrt.Pref^2)
    </td>
    <td style="text-align:right;">
    1.0504218
    </td>
    <td style="text-align:right;">
    0.1151067
    </td>
    <td style="text-align:right;">
    0.5631589
    </td>
    <td style="text-align:left;">
    -   </td>
        <td style="text-align:left;">
        Quadratic
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Complex
        </td>
        <td style="text-align:left;">
        sc.Diam
        </td>
        <td style="text-align:right;">
        1.6119663
        </td>
        <td style="text-align:right;">
        0.7495632
        </td>
        <td style="text-align:right;">
        1.1492365
        </td>
        <td style="text-align:left;">
        -   </td>
            <td style="text-align:left;">
            Linear
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Complex
            </td>
            <td style="text-align:left;">
            sc.Diam:sc.log.Clutch
            </td>
            <td style="text-align:right;">
            0.2603304
            </td>
            <td style="text-align:right;">
            -0.5360900
            </td>
            <td style="text-align:right;">
            -0.1422864
            </td>
            <td style="text-align:left;">
            </td>
            <td style="text-align:left;">
            Correlational
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Complex
            </td>
            <td style="text-align:left;">
            sc.Diam:sc.sqrt.Pref
            </td>
            <td style="text-align:right;">
            0.0553416
            </td>
            <td style="text-align:right;">
            -0.9712248
            </td>
            <td style="text-align:right;">
            -0.4424390
            </td>
            <td style="text-align:left;">
            </td>
            <td style="text-align:left;">
            Correlational
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Complex
            </td>
            <td style="text-align:left;">
            sc.log.Clutch
            </td>
            <td style="text-align:right;">
            0.5748189
            </td>
            <td style="text-align:right;">
            -0.1657583
            </td>
            <td style="text-align:right;">
            0.2018117
            </td>
            <td style="text-align:left;">
            </td>
            <td style="text-align:left;">
            Linear
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Complex
            </td>
            <td style="text-align:left;">
            sc.log.Clutch:sc.sqrt.Pref
            </td>
            <td style="text-align:right;">
            0.5880679
            </td>
            <td style="text-align:right;">
            -0.3441780
            </td>
            <td style="text-align:right;">
            0.0937689
            </td>
            <td style="text-align:left;">
            </td>
            <td style="text-align:left;">
            Correlational
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Complex
            </td>
            <td style="text-align:left;">
            sc.sqrt.Pref
            </td>
            <td style="text-align:right;">
            0.1508519
            </td>
            <td style="text-align:right;">
            -0.9751303
            </td>
            <td style="text-align:right;">
            -0.4242847
            </td>
            <td style="text-align:left;">
            </td>
            <td style="text-align:left;">
            Linear
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Diff
            </td>
            <td style="text-align:left;">
            (Intercept)
            </td>
            <td style="text-align:right;">
            0.4247075
            </td>
            <td style="text-align:right;">
            0.1147164
            </td>
            <td style="text-align:right;">
            0.2670696
            </td>
            <td style="text-align:left;">
            -   </td>
                <td style="text-align:left;">
                Mean fitness
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                I(sc.Diam^2)
                </td>
                <td style="text-align:right;">
                0.5294304
                </td>
                <td style="text-align:right;">
                -0.4235305
                </td>
                <td style="text-align:right;">
                0.0468140
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Quadratic
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                I(sc.log.Clutch^2)
                </td>
                <td style="text-align:right;">
                0.3326034
                </td>
                <td style="text-align:right;">
                -0.7287654
                </td>
                <td style="text-align:right;">
                -0.2135061
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Quadratic
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                I(sc.sqrt.Pref^2)
                </td>
                <td style="text-align:right;">
                0.0924156
                </td>
                <td style="text-align:right;">
                -1.1487505
                </td>
                <td style="text-align:right;">
                -0.5187575
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Quadratic
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                sc.Diam
                </td>
                <td style="text-align:right;">
                0.5028383
                </td>
                <td style="text-align:right;">
                -0.5493074
                </td>
                <td style="text-align:right;">
                -0.0448235
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Linear
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                sc.Diam:sc.log.Clutch
                </td>
                <td style="text-align:right;">
                0.4185971
                </td>
                <td style="text-align:right;">
                -0.7875250
                </td>
                <td style="text-align:right;">
                -0.2035389
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Correlational
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                sc.Diam:sc.sqrt.Pref
                </td>
                <td style="text-align:right;">
                1.0497497
                </td>
                <td style="text-align:right;">
                -0.3029471
                </td>
                <td style="text-align:right;">
                0.3571216
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Correlational
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                sc.log.Clutch
                </td>
                <td style="text-align:right;">
                -0.1656684
                </td>
                <td style="text-align:right;">
                -1.2777033
                </td>
                <td style="text-align:right;">
                -0.6704726
                </td>
                <td style="text-align:left;">
                -   </td>
                    <td style="text-align:left;">
                    Linear
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.log.Clutch:sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    0.4608604
                    </td>
                    <td style="text-align:right;">
                    -0.6612368
                    </td>
                    <td style="text-align:right;">
                    -0.0776872
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Correlational
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    0.2738247
                    </td>
                    <td style="text-align:right;">
                    -1.1523763
                    </td>
                    <td style="text-align:right;">
                    -0.4210804
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Linear
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Simple
                    </td>
                    <td style="text-align:left;">
                    (Intercept)
                    </td>
                    <td style="text-align:right;">
                    0.8081110
                    </td>
                    <td style="text-align:right;">
                    0.5388474
                    </td>
                    <td style="text-align:right;">
                    0.6828141
                    </td>
                    <td style="text-align:left;">
                    -   </td>
                        <td style="text-align:left;">
                        Mean fitness
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        I(sc.Diam^2)
                        </td>
                        <td style="text-align:right;">
                        0.6164037
                        </td>
                        <td style="text-align:right;">
                        -0.0492869
                        </td>
                        <td style="text-align:right;">
                        0.2693382
                        </td>
                        <td style="text-align:left;">
                        </td>
                        <td style="text-align:left;">
                        Quadratic
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        I(sc.log.Clutch^2)
                        </td>
                        <td style="text-align:right;">
                        0.0825847
                        </td>
                        <td style="text-align:right;">
                        -0.7411326
                        </td>
                        <td style="text-align:right;">
                        -0.3037629
                        </td>
                        <td style="text-align:left;">
                        </td>
                        <td style="text-align:left;">
                        Quadratic
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        I(sc.sqrt.Pref^2)
                        </td>
                        <td style="text-align:right;">
                        0.4682471
                        </td>
                        <td style="text-align:right;">
                        -0.3784503
                        </td>
                        <td style="text-align:right;">
                        0.0444014
                        </td>
                        <td style="text-align:left;">
                        </td>
                        <td style="text-align:left;">
                        Quadratic
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        sc.Diam
                        </td>
                        <td style="text-align:right;">
                        1.6499838
                        </td>
                        <td style="text-align:right;">
                        0.6390364
                        </td>
                        <td style="text-align:right;">
                        1.1044130
                        </td>
                        <td style="text-align:left;">
                        -   </td>
                            <td style="text-align:left;">
                            Linear
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Simple
                            </td>
                            <td style="text-align:left;">
                            sc.Diam:sc.log.Clutch
                            </td>
                            <td style="text-align:right;">
                            0.0881996
                            </td>
                            <td style="text-align:right;">
                            -0.8050264
                            </td>
                            <td style="text-align:right;">
                            -0.3458253
                            </td>
                            <td style="text-align:left;">
                            </td>
                            <td style="text-align:left;">
                            Correlational
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Simple
                            </td>
                            <td style="text-align:left;">
                            sc.Diam:sc.sqrt.Pref
                            </td>
                            <td style="text-align:right;">
                            0.3617315
                            </td>
                            <td style="text-align:right;">
                            -0.5381768
                            </td>
                            <td style="text-align:right;">
                            -0.0853174
                            </td>
                            <td style="text-align:left;">
                            </td>
                            <td style="text-align:left;">
                            Correlational
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Simple
                            </td>
                            <td style="text-align:left;">
                            sc.log.Clutch
                            </td>
                            <td style="text-align:right;">
                            -0.0607767
                            </td>
                            <td style="text-align:right;">
                            -0.9247423
                            </td>
                            <td style="text-align:right;">
                            -0.4686609
                            </td>
                            <td style="text-align:left;">
                            -   </td>
                                <td style="text-align:left;">
                                Linear
                                </td>
                                </tr>
                                <tr>
                                <td style="text-align:left;">
                                Simple
                                </td>
                                <td style="text-align:left;">
                                sc.log.Clutch:sc.sqrt.Pref
                                </td>
                                <td style="text-align:right;">
                                0.3773772
                                </td>
                                <td style="text-align:right;">
                                -0.3782373
                                </td>
                                <td style="text-align:right;">
                                0.0160817
                                </td>
                                <td style="text-align:left;">
                                </td>
                                <td style="text-align:left;">
                                Correlational
                                </td>
                                </tr>
                                <tr>
                                <td style="text-align:left;">
                                Simple
                                </td>
                                <td style="text-align:left;">
                                sc.sqrt.Pref
                                </td>
                                <td style="text-align:right;">
                                -0.3393616
                                </td>
                                <td style="text-align:right;">
                                -1.3661734
                                </td>
                                <td style="text-align:right;">
                                -0.8453650
                                </td>
                                <td style="text-align:left;">
                                -   </td>
                                    <td style="text-align:left;">
                                    Linear
                                    </td>
                                    </tr>
                                    </tbody>
                                    </table>

Plot these trait-fitness relationships:

``` r
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
  scale_y_continuous(name="larval survival", limits=c(0,1), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
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

![](reproduce_analyses_files/figure-markdown_github/Plot%20food-web%20coefficients-1.pdf)

Using the methods proposed by [Janzen and Stern, 1998](https://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.1998.tb02237.x), we estimated selection gradients from our trait-fitness relationships.

``` r
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

    ## Warning: Expected 3 pieces. Missing pieces filled with `NA` in 54 rows
    ## [1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 19, 20, 21, 22, 23, 24, 28,
    ## 29, ...].

Table of selection gradients.

``` r
knitr::kable(tidy_foodweb_grads)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
type
</th>
<th style="text-align:left;">
term
</th>
<th style="text-align:right;">
conf.high
</th>
<th style="text-align:right;">
conf.low
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:left;">
P\_cutoff
</th>
<th style="text-align:left;">
gradient\_type
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Complex
</td>
<td style="text-align:left;">
I(sc.Diam^2)
</td>
<td style="text-align:right;">
0.3317272
</td>
<td style="text-align:right;">
-0.0646534
</td>
<td style="text-align:right;">
0.1333655
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
Quadratic
</td>
</tr>
<tr>
<td style="text-align:left;">
Complex
</td>
<td style="text-align:left;">
I(sc.log.Clutch^2)
</td>
<td style="text-align:right;">
0.1782414
</td>
<td style="text-align:right;">
-0.2680121
</td>
<td style="text-align:right;">
-0.0540937
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
Quadratic
</td>
</tr>
<tr>
<td style="text-align:left;">
Complex
</td>
<td style="text-align:left;">
I(sc.sqrt.Pref^2)
</td>
<td style="text-align:right;">
0.6295498
</td>
<td style="text-align:right;">
0.0689869
</td>
<td style="text-align:right;">
0.3375183
</td>
<td style="text-align:left;">
-   </td>
    <td style="text-align:left;">
    Quadratic
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    sc.Diam
    </td>
    <td style="text-align:right;">
    0.4830502
    </td>
    <td style="text-align:right;">
    0.2246180
    </td>
    <td style="text-align:right;">
    0.3443862
    </td>
    <td style="text-align:left;">
    -   </td>
        <td style="text-align:left;">
        Directional
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Complex
        </td>
        <td style="text-align:left;">
        sc.Diam:sc.log.Clutch
        </td>
        <td style="text-align:right;">
        0.0780120
        </td>
        <td style="text-align:right;">
        -0.1606475
        </td>
        <td style="text-align:right;">
        -0.0426383
        </td>
        <td style="text-align:left;">
        </td>
        <td style="text-align:left;">
        Correlational
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Complex
        </td>
        <td style="text-align:left;">
        sc.Diam:sc.sqrt.Pref
        </td>
        <td style="text-align:right;">
        0.0165840
        </td>
        <td style="text-align:right;">
        -0.2910423
        </td>
        <td style="text-align:right;">
        -0.1325836
        </td>
        <td style="text-align:left;">
        </td>
        <td style="text-align:left;">
        Correlational
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Complex
        </td>
        <td style="text-align:left;">
        sc.log.Clutch
        </td>
        <td style="text-align:right;">
        0.1722532
        </td>
        <td style="text-align:right;">
        -0.0496720
        </td>
        <td style="text-align:right;">
        0.0604759
        </td>
        <td style="text-align:left;">
        </td>
        <td style="text-align:left;">
        Directional
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Complex
        </td>
        <td style="text-align:left;">
        sc.log.Clutch:sc.sqrt.Pref
        </td>
        <td style="text-align:right;">
        0.1762235
        </td>
        <td style="text-align:right;">
        -0.1031382
        </td>
        <td style="text-align:right;">
        0.0280993
        </td>
        <td style="text-align:left;">
        </td>
        <td style="text-align:left;">
        Correlational
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Complex
        </td>
        <td style="text-align:left;">
        sc.sqrt.Pref
        </td>
        <td style="text-align:right;">
        0.0452051
        </td>
        <td style="text-align:right;">
        -0.2922126
        </td>
        <td style="text-align:right;">
        -0.1271434
        </td>
        <td style="text-align:left;">
        </td>
        <td style="text-align:left;">
        Directional
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Diff
        </td>
        <td style="text-align:left;">
        I(sc.Diam^2)
        </td>
        <td style="text-align:right;">
        0.2012363
        </td>
        <td style="text-align:right;">
        -0.2680247
        </td>
        <td style="text-align:right;">
        -0.0320495
        </td>
        <td style="text-align:left;">
        </td>
        <td style="text-align:left;">
        Quadratic
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Diff
        </td>
        <td style="text-align:left;">
        I(sc.log.Clutch^2)
        </td>
        <td style="text-align:right;">
        0.2032530
        </td>
        <td style="text-align:right;">
        -0.3170832
        </td>
        <td style="text-align:right;">
        -0.0601718
        </td>
        <td style="text-align:left;">
        </td>
        <td style="text-align:left;">
        Quadratic
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Diff
        </td>
        <td style="text-align:left;">
        I(sc.sqrt.Pref^2)
        </td>
        <td style="text-align:right;">
        -0.0094026
        </td>
        <td style="text-align:right;">
        -0.6398557
        </td>
        <td style="text-align:right;">
        -0.3208160
        </td>
        <td style="text-align:left;">
        -   </td>
            <td style="text-align:left;">
            Quadratic
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Diff
            </td>
            <td style="text-align:left;">
            sc.Diam
            </td>
            <td style="text-align:right;">
            -0.0003254
            </td>
            <td style="text-align:right;">
            -0.2657059
            </td>
            <td style="text-align:right;">
            -0.1366645
            </td>
            <td style="text-align:left;">
            -   </td>
                <td style="text-align:left;">
                Directional
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                sc.Diam:sc.log.Clutch
                </td>
                <td style="text-align:right;">
                0.1245315
                </td>
                <td style="text-align:right;">
                -0.1684755
                </td>
                <td style="text-align:right;">
                -0.0224057
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Correlational
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                sc.Diam:sc.sqrt.Pref
                </td>
                <td style="text-align:right;">
                0.3046798
                </td>
                <td style="text-align:right;">
                -0.0502934
                </td>
                <td style="text-align:right;">
                0.1165368
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Correlational
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                sc.log.Clutch
                </td>
                <td style="text-align:right;">
                -0.0305383
                </td>
                <td style="text-align:right;">
                -0.2927925
                </td>
                <td style="text-align:right;">
                -0.1486232
                </td>
                <td style="text-align:left;">
                -   </td>
                    <td style="text-align:left;">
                    Directional
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.log.Clutch:sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    0.1166350
                    </td>
                    <td style="text-align:right;">
                    -0.1834444
                    </td>
                    <td style="text-align:right;">
                    -0.0250746
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Correlational
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    0.1486780
                    </td>
                    <td style="text-align:right;">
                    -0.2106648
                    </td>
                    <td style="text-align:right;">
                    -0.0318557
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Directional
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Simple
                    </td>
                    <td style="text-align:left;">
                    I(sc.Diam^2)
                    </td>
                    <td style="text-align:right;">
                    0.2318705
                    </td>
                    <td style="text-align:right;">
                    -0.0185401
                    </td>
                    <td style="text-align:right;">
                    0.1013161
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Quadratic
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Simple
                    </td>
                    <td style="text-align:left;">
                    I(sc.log.Clutch^2)
                    </td>
                    <td style="text-align:right;">
                    0.0310656
                    </td>
                    <td style="text-align:right;">
                    -0.2787894
                    </td>
                    <td style="text-align:right;">
                    -0.1142655
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Quadratic
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Simple
                    </td>
                    <td style="text-align:left;">
                    I(sc.sqrt.Pref^2)
                    </td>
                    <td style="text-align:right;">
                    0.1761390
                    </td>
                    <td style="text-align:right;">
                    -0.1423604
                    </td>
                    <td style="text-align:right;">
                    0.0167023
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Quadratic
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Simple
                    </td>
                    <td style="text-align:left;">
                    sc.Diam
                    </td>
                    <td style="text-align:right;">
                    0.3103345
                    </td>
                    <td style="text-align:right;">
                    0.1201921
                    </td>
                    <td style="text-align:right;">
                    0.2077217
                    </td>
                    <td style="text-align:left;">
                    -   </td>
                        <td style="text-align:left;">
                        Directional
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        sc.Diam:sc.log.Clutch
                        </td>
                        <td style="text-align:right;">
                        0.0165889
                        </td>
                        <td style="text-align:right;">
                        -0.1514121
                        </td>
                        <td style="text-align:right;">
                        -0.0650440
                        </td>
                        <td style="text-align:left;">
                        </td>
                        <td style="text-align:left;">
                        Correlational
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        sc.Diam:sc.sqrt.Pref
                        </td>
                        <td style="text-align:right;">
                        0.0680357
                        </td>
                        <td style="text-align:right;">
                        -0.1012221
                        </td>
                        <td style="text-align:right;">
                        -0.0160468
                        </td>
                        <td style="text-align:left;">
                        </td>
                        <td style="text-align:left;">
                        Correlational
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        sc.log.Clutch
                        </td>
                        <td style="text-align:right;">
                        -0.0114311
                        </td>
                        <td style="text-align:right;">
                        -0.1739286
                        </td>
                        <td style="text-align:right;">
                        -0.0881473
                        </td>
                        <td style="text-align:left;">
                        -   </td>
                            <td style="text-align:left;">
                            Directional
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Simple
                            </td>
                            <td style="text-align:left;">
                            sc.log.Clutch:sc.sqrt.Pref
                            </td>
                            <td style="text-align:right;">
                            0.0709784
                            </td>
                            <td style="text-align:right;">
                            -0.0711401
                            </td>
                            <td style="text-align:right;">
                            0.0030247
                            </td>
                            <td style="text-align:left;">
                            </td>
                            <td style="text-align:left;">
                            Correlational
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Simple
                            </td>
                            <td style="text-align:left;">
                            sc.sqrt.Pref
                            </td>
                            <td style="text-align:right;">
                            -0.0638283
                            </td>
                            <td style="text-align:right;">
                            -0.2569545
                            </td>
                            <td style="text-align:right;">
                            -0.1589991
                            </td>
                            <td style="text-align:left;">
                            -   </td>
                                <td style="text-align:left;">
                                Directional
                                </td>
                                </tr>
                                </tbody>
                                </table>

Plot of selection gradients.

``` r
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

![](reproduce_analyses_files/figure-markdown_github/Plot%20food-web%20gradients-1.pdf)

Partitioning the contribution of egg and larval parasitoids to selection gradients
==================================================================================

Our simple food-web treatment allows us to estimate the unique contribution of egg parasitoids to selection on *Iteomyia* traits. To estimate the unique contribution of larval parasitoids, we subset our data so that our complex food-web treatment only contained attack by larval parasitoids (and gall survival). We then fit the same models as previously, including one to estimate bias.

``` r
# excludes cases of egg-parasitism from Complex food web
egglarval_df <- filter(gall_selection.df, Foodweb == "Simple" | Foodweb == "Complex" & platy < 1) 

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

    ## singular fit

Parametric bootstrapping to estimate confidence intervals

``` r
boot_egglarval_model <- bootMer(egglarval_model, FUN = fixef, nsim=n_boots_analysis, parallel="multicore", ncpus=32, seed=34)
```

Fit a reduced model to reliably estimate linear trait-fitness relationships.

``` r
linear_egglarval_model <- update(linear_foodweb_model, data=egglarval_df)
```

Parametric bootstrapping to estimate confidence intervals.

``` r
boot_linear_egglarval_model <- bootMer(linear_egglarval_model, FUN = fixef, nsim=n_boots_analysis, parallel="multicore", ncpus=32, seed=34)
```

Tidy trait-fitness relationships for egg and larval parasitoids from above models.

``` r
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

    ## Warning: Expected 3 pieces. Missing pieces filled with `NA` in 63 rows [1,
    ## 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 19, 20, 21, 22, 23, 24, 25, 26, ...].

Table of trait-fitness relationships for egg and larval parasitoids.

``` r
knitr::kable(tidy_egglarval_alphas)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
type
</th>
<th style="text-align:left;">
term
</th>
<th style="text-align:right;">
conf.high
</th>
<th style="text-align:right;">
conf.low
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:left;">
P\_cutoff
</th>
<th style="text-align:left;">
coefficient\_type
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Complex
</td>
<td style="text-align:left;">
(Intercept)
</td>
<td style="text-align:right;">
0.8321860
</td>
<td style="text-align:right;">
0.5083217
</td>
<td style="text-align:right;">
0.6784569
</td>
<td style="text-align:left;">
-   </td>
    <td style="text-align:left;">
    Mean fitness
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    I(sc.Diam^2)
    </td>
    <td style="text-align:right;">
    0.5227084
    </td>
    <td style="text-align:right;">
    -0.1704994
    </td>
    <td style="text-align:right;">
    0.1626027
    </td>
    <td style="text-align:left;">
    </td>
    <td style="text-align:left;">
    Quadratic
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    I(sc.log.Clutch^2)
    </td>
    <td style="text-align:right;">
    0.3819561
    </td>
    <td style="text-align:right;">
    -0.3640244
    </td>
    <td style="text-align:right;">
    0.0140399
    </td>
    <td style="text-align:left;">
    </td>
    <td style="text-align:left;">
    Quadratic
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    I(sc.sqrt.Pref^2)
    </td>
    <td style="text-align:right;">
    1.0728129
    </td>
    <td style="text-align:right;">
    -0.0423653
    </td>
    <td style="text-align:right;">
    0.4671612
    </td>
    <td style="text-align:left;">
    </td>
    <td style="text-align:left;">
    Quadratic
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    sc.Diam
    </td>
    <td style="text-align:right;">
    1.8540837
    </td>
    <td style="text-align:right;">
    0.6716393
    </td>
    <td style="text-align:right;">
    1.1691793
    </td>
    <td style="text-align:left;">
    -   </td>
        <td style="text-align:left;">
        Linear
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Complex
        </td>
        <td style="text-align:left;">
        sc.Diam:sc.log.Clutch
        </td>
        <td style="text-align:right;">
        0.5498381
        </td>
        <td style="text-align:right;">
        -0.3421589
        </td>
        <td style="text-align:right;">
        0.0960534
        </td>
        <td style="text-align:left;">
        </td>
        <td style="text-align:left;">
        Correlational
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Complex
        </td>
        <td style="text-align:left;">
        sc.Diam:sc.sqrt.Pref
        </td>
        <td style="text-align:right;">
        0.3268931
        </td>
        <td style="text-align:right;">
        -0.9231589
        </td>
        <td style="text-align:right;">
        -0.2768599
        </td>
        <td style="text-align:left;">
        </td>
        <td style="text-align:left;">
        Correlational
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Complex
        </td>
        <td style="text-align:left;">
        sc.log.Clutch
        </td>
        <td style="text-align:right;">
        1.2616446
        </td>
        <td style="text-align:right;">
        0.2007682
        </td>
        <td style="text-align:right;">
        0.6685121
        </td>
        <td style="text-align:left;">
        -   </td>
            <td style="text-align:left;">
            Linear
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Complex
            </td>
            <td style="text-align:left;">
            sc.log.Clutch:sc.sqrt.Pref
            </td>
            <td style="text-align:right;">
            0.3300190
            </td>
            <td style="text-align:right;">
            -0.7104094
            </td>
            <td style="text-align:right;">
            -0.1995145
            </td>
            <td style="text-align:left;">
            </td>
            <td style="text-align:left;">
            Correlational
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Complex
            </td>
            <td style="text-align:left;">
            sc.sqrt.Pref
            </td>
            <td style="text-align:right;">
            -0.1690631
            </td>
            <td style="text-align:right;">
            -1.7324745
            </td>
            <td style="text-align:right;">
            -0.8965151
            </td>
            <td style="text-align:left;">
            -   </td>
                <td style="text-align:left;">
                Linear
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                (Intercept)
                </td>
                <td style="text-align:right;">
                0.2064694
                </td>
                <td style="text-align:right;">
                -0.1922930
                </td>
                <td style="text-align:right;">
                0.0102961
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Mean fitness
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                I(sc.Diam^2)
                </td>
                <td style="text-align:right;">
                0.5627544
                </td>
                <td style="text-align:right;">
                -0.3723047
                </td>
                <td style="text-align:right;">
                0.0896419
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Quadratic
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                I(sc.log.Clutch^2)
                </td>
                <td style="text-align:right;">
                0.2160491
                </td>
                <td style="text-align:right;">
                -0.8383143
                </td>
                <td style="text-align:right;">
                -0.3153519
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Quadratic
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                I(sc.sqrt.Pref^2)
                </td>
                <td style="text-align:right;">
                0.1552177
                </td>
                <td style="text-align:right;">
                -1.2688888
                </td>
                <td style="text-align:right;">
                -0.5274186
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Quadratic
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                sc.Diam
                </td>
                <td style="text-align:right;">
                0.5416764
                </td>
                <td style="text-align:right;">
                -0.7931427
                </td>
                <td style="text-align:right;">
                -0.0971937
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Linear
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                sc.Diam:sc.log.Clutch
                </td>
                <td style="text-align:right;">
                0.1310230
                </td>
                <td style="text-align:right;">
                -1.0606776
                </td>
                <td style="text-align:right;">
                -0.4551249
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Correlational
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                sc.Diam:sc.sqrt.Pref
                </td>
                <td style="text-align:right;">
                0.9174299
                </td>
                <td style="text-align:right;">
                -0.5432553
                </td>
                <td style="text-align:right;">
                0.1646462
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Correlational
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                sc.log.Clutch
                </td>
                <td style="text-align:right;">
                -0.7901097
                </td>
                <td style="text-align:right;">
                -2.2758162
                </td>
                <td style="text-align:right;">
                -1.4742391
                </td>
                <td style="text-align:left;">
                -   </td>
                    <td style="text-align:left;">
                    Linear
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.log.Clutch:sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    0.8791916
                    </td>
                    <td style="text-align:right;">
                    -0.3970962
                    </td>
                    <td style="text-align:right;">
                    0.2391500
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Correlational
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    0.8507252
                    </td>
                    <td style="text-align:right;">
                    -0.9613207
                    </td>
                    <td style="text-align:right;">
                    -0.1030111
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Linear
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Simple
                    </td>
                    <td style="text-align:left;">
                    (Intercept)
                    </td>
                    <td style="text-align:right;">
                    0.8536756
                    </td>
                    <td style="text-align:right;">
                    0.5030799
                    </td>
                    <td style="text-align:right;">
                    0.6887530
                    </td>
                    <td style="text-align:left;">
                    -   </td>
                        <td style="text-align:left;">
                        Mean fitness
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        I(sc.Diam^2)
                        </td>
                        <td style="text-align:right;">
                        0.5818657
                        </td>
                        <td style="text-align:right;">
                        -0.0359987
                        </td>
                        <td style="text-align:right;">
                        0.2522446
                        </td>
                        <td style="text-align:left;">
                        </td>
                        <td style="text-align:left;">
                        Quadratic
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        I(sc.log.Clutch^2)
                        </td>
                        <td style="text-align:right;">
                        0.0556405
                        </td>
                        <td style="text-align:right;">
                        -0.6690635
                        </td>
                        <td style="text-align:right;">
                        -0.3013121
                        </td>
                        <td style="text-align:left;">
                        </td>
                        <td style="text-align:left;">
                        Quadratic
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        I(sc.sqrt.Pref^2)
                        </td>
                        <td style="text-align:right;">
                        0.3815757
                        </td>
                        <td style="text-align:right;">
                        -0.5008657
                        </td>
                        <td style="text-align:right;">
                        -0.0602575
                        </td>
                        <td style="text-align:left;">
                        </td>
                        <td style="text-align:left;">
                        Quadratic
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        sc.Diam
                        </td>
                        <td style="text-align:right;">
                        1.6477512
                        </td>
                        <td style="text-align:right;">
                        0.6342528
                        </td>
                        <td style="text-align:right;">
                        1.0719856
                        </td>
                        <td style="text-align:left;">
                        -   </td>
                            <td style="text-align:left;">
                            Linear
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Simple
                            </td>
                            <td style="text-align:left;">
                            sc.Diam:sc.log.Clutch
                            </td>
                            <td style="text-align:right;">
                            0.0353198
                            </td>
                            <td style="text-align:right;">
                            -0.7708686
                            </td>
                            <td style="text-align:right;">
                            -0.3590715
                            </td>
                            <td style="text-align:left;">
                            </td>
                            <td style="text-align:left;">
                            Correlational
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Simple
                            </td>
                            <td style="text-align:left;">
                            sc.Diam:sc.sqrt.Pref
                            </td>
                            <td style="text-align:right;">
                            0.2690800
                            </td>
                            <td style="text-align:right;">
                            -0.4808723
                            </td>
                            <td style="text-align:right;">
                            -0.1122137
                            </td>
                            <td style="text-align:left;">
                            </td>
                            <td style="text-align:left;">
                            Correlational
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Simple
                            </td>
                            <td style="text-align:left;">
                            sc.log.Clutch
                            </td>
                            <td style="text-align:right;">
                            -0.3080803
                            </td>
                            <td style="text-align:right;">
                            -1.3766468
                            </td>
                            <td style="text-align:right;">
                            -0.8057270
                            </td>
                            <td style="text-align:left;">
                            -   </td>
                                <td style="text-align:left;">
                                Linear
                                </td>
                                </tr>
                                <tr>
                                <td style="text-align:left;">
                                Simple
                                </td>
                                <td style="text-align:left;">
                                sc.log.Clutch:sc.sqrt.Pref
                                </td>
                                <td style="text-align:right;">
                                0.4139023
                                </td>
                                <td style="text-align:right;">
                                -0.3069665
                                </td>
                                <td style="text-align:right;">
                                0.0396355
                                </td>
                                <td style="text-align:left;">
                                </td>
                                <td style="text-align:left;">
                                Correlational
                                </td>
                                </tr>
                                <tr>
                                <td style="text-align:left;">
                                Simple
                                </td>
                                <td style="text-align:left;">
                                sc.sqrt.Pref
                                </td>
                                <td style="text-align:right;">
                                -0.4711270
                                </td>
                                <td style="text-align:right;">
                                -1.6576630
                                </td>
                                <td style="text-align:right;">
                                -0.9995261
                                </td>
                                <td style="text-align:left;">
                                -   </td>
                                    <td style="text-align:left;">
                                    Linear
                                    </td>
                                    </tr>
                                    </tbody>
                                    </table>

Plot of trait-fitness relationships for egg and larval parasitoids.

``` r
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
  scale_y_continuous(name="larval survival", limits=c(0,1), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
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

![](reproduce_analyses_files/figure-markdown_github/Plot%20egglarval%20coefficients-1.pdf)

Calculate selection gradients from trait-fitness relationships.

``` r
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

    ## Warning: Expected 3 pieces. Missing pieces filled with `NA` in 54 rows
    ## [1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 19, 20, 21, 22, 23, 24, 28,
    ## 29, ...].

Table of selection gradients imposed by egg and larval parasitoids.

``` r
knitr::kable(tidy_egglarval_grads)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
type
</th>
<th style="text-align:left;">
term
</th>
<th style="text-align:right;">
conf.high
</th>
<th style="text-align:right;">
conf.low
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:left;">
P\_cutoff
</th>
<th style="text-align:left;">
gradient\_type
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Complex
</td>
<td style="text-align:left;">
I(sc.Diam^2)
</td>
<td style="text-align:right;">
0.2024558
</td>
<td style="text-align:right;">
-0.0660380
</td>
<td style="text-align:right;">
0.0629794
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
Quadratic
</td>
</tr>
<tr>
<td style="text-align:left;">
Complex
</td>
<td style="text-align:left;">
I(sc.log.Clutch^2)
</td>
<td style="text-align:right;">
0.1479395
</td>
<td style="text-align:right;">
-0.1409942
</td>
<td style="text-align:right;">
0.0054379
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
Quadratic
</td>
</tr>
<tr>
<td style="text-align:left;">
Complex
</td>
<td style="text-align:left;">
I(sc.sqrt.Pref^2)
</td>
<td style="text-align:right;">
0.4155227
</td>
<td style="text-align:right;">
-0.0164090
</td>
<td style="text-align:right;">
0.1809412
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
Quadratic
</td>
</tr>
<tr>
<td style="text-align:left;">
Complex
</td>
<td style="text-align:left;">
sc.Diam
</td>
<td style="text-align:right;">
0.3590625
</td>
<td style="text-align:right;">
0.1300699
</td>
<td style="text-align:right;">
0.2264237
</td>
<td style="text-align:left;">
-   </td>
    <td style="text-align:left;">
    Directional
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    sc.Diam:sc.log.Clutch
    </td>
    <td style="text-align:right;">
    0.1064818
    </td>
    <td style="text-align:right;">
    -0.0662626
    </td>
    <td style="text-align:right;">
    0.0186017
    </td>
    <td style="text-align:left;">
    </td>
    <td style="text-align:left;">
    Correlational
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    sc.Diam:sc.sqrt.Pref
    </td>
    <td style="text-align:right;">
    0.0633062
    </td>
    <td style="text-align:right;">
    -0.1787793
    </td>
    <td style="text-align:right;">
    -0.0536168
    </td>
    <td style="text-align:left;">
    </td>
    <td style="text-align:left;">
    Correlational
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    sc.log.Clutch
    </td>
    <td style="text-align:right;">
    0.2443305
    </td>
    <td style="text-align:right;">
    0.0388808
    </td>
    <td style="text-align:right;">
    0.1294643
    </td>
    <td style="text-align:left;">
    -   </td>
        <td style="text-align:left;">
        Directional
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Complex
        </td>
        <td style="text-align:left;">
        sc.log.Clutch:sc.sqrt.Pref
        </td>
        <td style="text-align:right;">
        0.0639116
        </td>
        <td style="text-align:right;">
        -0.1375781
        </td>
        <td style="text-align:right;">
        -0.0386380
        </td>
        <td style="text-align:left;">
        </td>
        <td style="text-align:left;">
        Correlational
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Complex
        </td>
        <td style="text-align:left;">
        sc.sqrt.Pref
        </td>
        <td style="text-align:right;">
        -0.0327408
        </td>
        <td style="text-align:right;">
        -0.3355116
        </td>
        <td style="text-align:right;">
        -0.1736194
        </td>
        <td style="text-align:left;">
        -   </td>
            <td style="text-align:left;">
            Directional
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Diff
            </td>
            <td style="text-align:left;">
            I(sc.Diam^2)
            </td>
            <td style="text-align:right;">
            0.2220957
            </td>
            <td style="text-align:right;">
            -0.1433090
            </td>
            <td style="text-align:right;">
            0.0368275
            </td>
            <td style="text-align:left;">
            </td>
            <td style="text-align:left;">
            Quadratic
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Diff
            </td>
            <td style="text-align:left;">
            I(sc.log.Clutch^2)
            </td>
            <td style="text-align:right;">
            0.0823091
            </td>
            <td style="text-align:right;">
            -0.3290462
            </td>
            <td style="text-align:right;">
            -0.1246596
            </td>
            <td style="text-align:left;">
            </td>
            <td style="text-align:left;">
            Quadratic
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Diff
            </td>
            <td style="text-align:left;">
            I(sc.sqrt.Pref^2)
            </td>
            <td style="text-align:right;">
            0.0626047
            </td>
            <td style="text-align:right;">
            -0.4925172
            </td>
            <td style="text-align:right;">
            -0.2047836
            </td>
            <td style="text-align:left;">
            </td>
            <td style="text-align:left;">
            Quadratic
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Diff
            </td>
            <td style="text-align:left;">
            sc.Diam
            </td>
            <td style="text-align:right;">
            0.1101885
            </td>
            <td style="text-align:right;">
            -0.1480847
            </td>
            <td style="text-align:right;">
            -0.0143447
            </td>
            <td style="text-align:left;">
            </td>
            <td style="text-align:left;">
            Directional
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Diff
            </td>
            <td style="text-align:left;">
            sc.Diam:sc.log.Clutch
            </td>
            <td style="text-align:right;">
            0.0257850
            </td>
            <td style="text-align:right;">
            -0.2075402
            </td>
            <td style="text-align:right;">
            -0.0896395
            </td>
            <td style="text-align:left;">
            </td>
            <td style="text-align:left;">
            Correlational
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Diff
            </td>
            <td style="text-align:left;">
            sc.Diam:sc.sqrt.Pref
            </td>
            <td style="text-align:right;">
            0.1780003
            </td>
            <td style="text-align:right;">
            -0.1074612
            </td>
            <td style="text-align:right;">
            0.0314167
            </td>
            <td style="text-align:left;">
            </td>
            <td style="text-align:left;">
            Correlational
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Diff
            </td>
            <td style="text-align:left;">
            sc.log.Clutch
            </td>
            <td style="text-align:right;">
            -0.1545573
            </td>
            <td style="text-align:right;">
            -0.4446473
            </td>
            <td style="text-align:right;">
            -0.2888673
            </td>
            <td style="text-align:left;">
            -   </td>
                <td style="text-align:left;">
                Directional
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                sc.log.Clutch:sc.sqrt.Pref
                </td>
                <td style="text-align:right;">
                0.1705295
                </td>
                <td style="text-align:right;">
                -0.0766891
                </td>
                <td style="text-align:right;">
                0.0464794
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Correlational
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                sc.sqrt.Pref
                </td>
                <td style="text-align:right;">
                0.1630412
                </td>
                <td style="text-align:right;">
                -0.1914796
                </td>
                <td style="text-align:right;">
                -0.0241243
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Directional
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Simple
                </td>
                <td style="text-align:left;">
                I(sc.Diam^2)
                </td>
                <td style="text-align:right;">
                0.2302297
                </td>
                <td style="text-align:right;">
                -0.0142438
                </td>
                <td style="text-align:right;">
                0.0998069
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Quadratic
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Simple
                </td>
                <td style="text-align:left;">
                I(sc.log.Clutch^2)
                </td>
                <td style="text-align:right;">
                0.0220155
                </td>
                <td style="text-align:right;">
                -0.2647317
                </td>
                <td style="text-align:right;">
                -0.1192216
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Quadratic
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Simple
                </td>
                <td style="text-align:left;">
                I(sc.sqrt.Pref^2)
                </td>
                <td style="text-align:right;">
                0.1509799
                </td>
                <td style="text-align:right;">
                -0.1981800
                </td>
                <td style="text-align:right;">
                -0.0238424
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Quadratic
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Simple
                </td>
                <td style="text-align:left;">
                sc.Diam
                </td>
                <td style="text-align:right;">
                0.3259870
                </td>
                <td style="text-align:right;">
                0.1254790
                </td>
                <td style="text-align:right;">
                0.2120789
                </td>
                <td style="text-align:left;">
                -   </td>
                    <td style="text-align:left;">
                    Directional
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Simple
                    </td>
                    <td style="text-align:left;">
                    sc.Diam:sc.log.Clutch
                    </td>
                    <td style="text-align:right;">
                    0.0069876
                    </td>
                    <td style="text-align:right;">
                    -0.1525067
                    </td>
                    <td style="text-align:right;">
                    -0.0710378
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Correlational
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Simple
                    </td>
                    <td style="text-align:left;">
                    sc.Diam:sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    0.0532341
                    </td>
                    <td style="text-align:right;">
                    -0.0951346
                    </td>
                    <td style="text-align:right;">
                    -0.0222001
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Correlational
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Simple
                    </td>
                    <td style="text-align:left;">
                    sc.log.Clutch
                    </td>
                    <td style="text-align:right;">
                    -0.0609498
                    </td>
                    <td style="text-align:right;">
                    -0.2723523
                    </td>
                    <td style="text-align:right;">
                    -0.1594030
                    </td>
                    <td style="text-align:left;">
                    -   </td>
                        <td style="text-align:left;">
                        Directional
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        sc.log.Clutch:sc.sqrt.Pref
                        </td>
                        <td style="text-align:right;">
                        0.0818854
                        </td>
                        <td style="text-align:right;">
                        -0.0607295
                        </td>
                        <td style="text-align:right;">
                        0.0078414
                        </td>
                        <td style="text-align:left;">
                        </td>
                        <td style="text-align:left;">
                        Correlational
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        sc.sqrt.Pref
                        </td>
                        <td style="text-align:right;">
                        -0.0932066
                        </td>
                        <td style="text-align:right;">
                        -0.3279479
                        </td>
                        <td style="text-align:right;">
                        -0.1977437
                        </td>
                        <td style="text-align:left;">
                        -   </td>
                            <td style="text-align:left;">
                            Directional
                            </td>
                            </tr>
                            </tbody>
                            </table>

Plot of selection gradients imposed by egg and larval parasitoids.

``` r
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

![](reproduce_analyses_files/figure-markdown_github/Plot%20egg%20vs.%20larval%20ptoid%20gradients-1.pdf)

We combine our estimates of selection gradients for each food-web treatment as well as the contribution of larval parasitoids to selection in the complex food web.

Partitioning the components of selection gradients
==================================================

Selection gradients are influenced by both trait-fitness relationships and population mean fitness. Here, we partition selection gradients into these underlying components.

Reproduce Figure 2 in main text of manuscript.
==============================================

Use the fitted model to generate predicted estimates of mean fitness for changes in the mean trait value of a population. I restrict these predictions to +/- 1 SD because we can only reliably estimate the shape of the adaptive landscape near the mean phenotype of the population [Arnold et al., 2001](https://link.springer.com/article/10.1023/A:1013373907708). Other trait values are held constant at the mean phenotype (i.e. trait = 0).

``` r
newdata_Diam <- bind_rows(
  expand.grid(Foodweb = "Complex", sc.Diam = seq(-1,1,length.out=1000), sc.log.Clutch = 0, sc.sqrt.Pref = 0),
  expand.grid(Foodweb = "Simple", sc.Diam = seq(-1,1,length.out=1000), sc.log.Clutch = 0, sc.sqrt.Pref = 0))
```

    ## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

``` r
RF_Diam <- bootstrap_fitness(
  logistic_model = foodweb_model, 
  newdata = newdata_Diam,
  bootstraps=n_boots_plots)
```

``` r
newdata_Clutch <- bind_rows(
  expand.grid(Foodweb = "Complex", sc.Diam = 0, sc.log.Clutch = seq(-1,1,length.out=1000), sc.sqrt.Pref = 0),
  expand.grid(Foodweb = "Simple", sc.Diam = 0, sc.log.Clutch = seq(-1,1,length.out=1000), sc.sqrt.Pref = 0))
```

    ## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

``` r
RF_Clutch <- bootstrap_fitness(
  logistic_model = foodweb_model,
  newdata = newdata_Clutch,
  bootstraps=n_boots_plots)
```

``` r
newdata_Pref <- bind_rows(
  expand.grid(Foodweb = "Complex", sc.Diam = 0, sc.log.Clutch = 0, sc.sqrt.Pref = seq(-1,1,length.out=1000)),
  expand.grid(Foodweb = "Simple", sc.Diam = 0, sc.log.Clutch = 0, sc.sqrt.Pref = seq(-1,1,length.out=1000)))
```

    ## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

``` r
RF_Pref <- bootstrap_fitness(
  logistic_model = foodweb_model, 
  newdata = newdata_Pref,
  bootstraps=n_boots_plots)
```

Create plots for each panel of Figure 2 in main text.

``` r
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

``` r
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

![Selection gradients acting on gall traits in complex vs. simple food webs. Each panel corresponds to a different gall trait: gall diameter (A); clutch size (B); and Oviposition preference (C). Solid lines represent the estimated gradients in complex (orange) and simple (blue) food webs. Transparent lines represent bootstrapped replicates (n=100) to show the uncertainty in estimated gradients. Note that only 100 bootstraps are displayed here, but that inferences are based on 1,000 bootstrapped samples.](reproduce_analyses_files/figure-markdown_github/Univariate-Fitness-Landscapes-1.pdf)

``` r
save_plot(filename = "UV_landscapes.pdf", plot = AF_gradients, base_height = 5, base_width = 8.5)
```

Replot univariate adaptive landscapes in terms of relative fitness
==================================================================

``` r
RF_plot_foodweb_Clutch_df <- RF_Clutch$absolute_fitness %>% 
  select(-sc.Diam, -sc.sqrt.Pref) %>%
  gather(ID, absolute_fitness, -sc.log.Clutch, -Foodweb) %>%
  mutate(ID_group = ifelse(ID == "average", "average", "replicate")) %>%
  #filter(ID_group == "average") %>%
  group_by(Foodweb) %>%
  mutate(relative_fitness = absolute_fitness/mean(absolute_fitness)) %>%
  ungroup()

#RF_plot_egglarval_Clutch_df <- egglarval_RF_Clutch$absolute_fitness %>% 
#  select(-sc.Diam, -sc.sqrt.Pref) %>%
#  gather(ID, absolute_fitness, -sc.log.Clutch, -Foodweb) %>%
#  mutate(ID_group = ifelse(ID == "average", "average", "replicate")) %>%
#  mutate(Foodweb = ifelse(Foodweb=="Complex","Larval","Simple")) %>%
#  filter(Foodweb=="Larval", ID_group=="average") %>%
#  group_by(Foodweb) %>%
#  mutate(relative_fitness = absolute_fitness/mean(absolute_fitness))%>%
#  ungroup()

RF_uni_Clutch <- RF_plot_foodweb_Clutch_df %>%
#  mutate(Foodweb=factor(Foodweb, levels=c("Complex","Simple","Larval"), labels=c("Egg + Larval (Complex)","Egg (Simple)","Larval"))) %>%
#  filter(ID_group=="average") %>%
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

## Female pref

RF_plot_foodweb_Pref_df <- RF_Pref$absolute_fitness %>% 
  select(-sc.Diam, -sc.log.Clutch) %>%
  gather(ID, absolute_fitness, -sc.sqrt.Pref, -Foodweb) %>%
  mutate(ID_group = ifelse(ID == "average", "average", "replicate")) %>%
  #filter(ID_group == "average") %>%
  group_by(Foodweb) %>%
  mutate(relative_fitness = absolute_fitness/mean(absolute_fitness))%>%
  ungroup()

#RF_plot_egglarval_Pref_df <- egglarval_RF_Pref$absolute_fitness %>%
#  select(-sc.Diam, -sc.log.Clutch) %>%
#  gather(ID, absolute_fitness, -sc.sqrt.Pref, -Foodweb) %>%
#  mutate(ID_group = ifelse(ID == "average", "average", "replicate")) %>%
#  mutate(Foodweb = ifelse(Foodweb=="Complex","Larval","Simple")) %>%
#  filter(Foodweb=="Larval", ID_group=="average") %>%
#  group_by(Foodweb) %>%
#  mutate(relative_fitness = absolute_fitness/mean(absolute_fitness))%>%
#  ungroup()

RF_uni_Pref <- RF_plot_foodweb_Pref_df %>%
#  mutate(Foodweb=factor(Foodweb, levels=c("Complex","Simple","Larval"), labels=c("Egg + Larval (Complex)","Egg (Simple)","Larval"))) %>%
#  filter(ID_group=="average") %>%
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
  #filter(ID_group == "average") %>%
  group_by(Foodweb) %>%
  mutate(relative_fitness = absolute_fitness/mean(absolute_fitness))%>%
  ungroup()

#RF_plot_egglarval_Diam_df <- egglarval_RF_Diam$absolute_fitness %>% 
#  select(-sc.sqrt.Pref, -sc.log.Clutch) %>%
#  gather(ID, absolute_fitness, -sc.Diam, -Foodweb) %>%
#  mutate(ID_group = ifelse(ID == "average", "average", "replicate")) %>%
#  mutate(Foodweb = ifelse(Foodweb=="Complex","Larval","Simple")) %>%
#  filter(Foodweb=="Larval", ID_group=="average") %>%
#  group_by(Foodweb) %>%
#  mutate(relative_fitness = absolute_fitness/mean(absolute_fitness))%>%
#  ungroup()

RF_uni_Diam <- RF_plot_foodweb_Diam_df %>%
#  mutate(Foodweb=factor(Foodweb, levels=c("Complex","Simple","Larval"), labels=c("Egg + Larval (Complex)","Egg (Simple)","Larval"))) %>%
#  filter(ID_group=="average") %>%
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

![](reproduce_analyses_files/figure-markdown_github/Plot%20Univariate%20RF%20Landscapes-1.pdf)

Reproduce Figure 3 presented in main text.
==========================================

Use the fitted model to generate predicted estimates of mean fitness for different phenotypic combinations. I restrict these predictions to +/- 1 SD because we can only reliably estimate the shape of the adaptive landscape near the mean phenotype of the population [Arnold et al., 2001](https://link.springer.com/article/10.1023/A:1013373907708).

``` r
newdata_Clutch.Pref <- bind_rows(
  expand.grid(Foodweb = "Complex", sc.Diam = 0, sc.log.Clutch = seq(-1,1,length.out=1000), sc.sqrt.Pref = seq(-1,1,length.out=1000)),
  expand.grid(Foodweb = "Simple", sc.Diam = 0, sc.log.Clutch = seq(-1,1,length.out=1000), sc.sqrt.Pref = seq(-1,1,length.out=1000)))
```

    ## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

``` r
RF_Clutch.Pref <- bootstrap_fitness(
  logistic_model = foodweb_model, 
  newdata = newdata_Clutch.Pref,
  bootstraps=NULL)
```

``` r
newdata_Diam.Clutch <- bind_rows(
  expand.grid(Foodweb = "Complex", sc.Diam = seq(-1,1,length.out=1000), sc.log.Clutch = seq(-1,1,length.out=1000), sc.sqrt.Pref = 0),
  expand.grid(Foodweb = "Simple", sc.Diam = seq(-1,1,length.out=1000), sc.log.Clutch = seq(-1,1,length.out=1000), sc.sqrt.Pref=0))
```

    ## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

``` r
RF_Diam.Clutch <- bootstrap_fitness(
  logistic_model = foodweb_model, 
  newdata = newdata_Diam.Clutch,
  bootstraps=NULL)
```

``` r
newdata_Diam.Pref <- bind_rows(
  expand.grid(Foodweb = "Complex", sc.Diam = seq(-1,1,length.out=1000), sc.log.Clutch = 0, sc.sqrt.Pref = seq(-1,1,length.out=1000)),
  expand.grid(Foodweb = "Simple", sc.Diam = seq(-1,1,length.out=1000), sc.log.Clutch = 0, sc.sqrt.Pref = seq(-1,1,length.out=1000)))
```

    ## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

``` r
RF_Diam.Pref <- bootstrap_fitness(
  logistic_model = foodweb_model, 
  newdata = newdata_Diam.Pref,
  bootstraps=NULL)
```

Create plots for each panel.

``` r
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

``` r
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

![](reproduce_analyses_files/figure-markdown_github/Figure-Multivariate-Landscapes-1.pdf)

``` r
save_plot(filename = "MV_landscapes.pdf", plot = AF_landscape_2d, base_width=8, base_height = 6)
```

Calculating Slope and Curvature of the Adaptive Landscape
=========================================================

Below is the code I used to get estimates for the slope and curvature of the adaptive landscape in each food-web treatment. If 95% confidence intervals of selection gradients overlap zero, then they are set to zero; otherwise, I retain the mean estimate.

``` r
## Create beta-matrix for complex food-web treatment
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

## Create gamma-matrix for complex food-web treatment
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

## Create matrix of Betas in simple food-web treatment
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

## Create gamma-matrix for simple food-web treatment
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

## Calculate Curvature of the Fitness landscape
complex_curvature <- complex_gammas - complex_betas %*% t(complex_betas)
simple_curvature <- simple_gammas - simple_betas %*% t(simple_betas)
```

I use these estimates to populate the curvature matrices presented in the **Results** of the main text. I reproduce these estimates here for clarity:

$$\\textbf{C} = \\begin{pmatrix} C\_{\\text{Diam:Diam}}&& \\\\ C\_{\\text{Clutch:Diam}} & C\_{\\text{Clutch:Clutch}} & \\\\ C\_{\\text{Pref:Diam}} & C\_{\\text{Pref:Clutch}} & C\_{\\text{Pref:Pref}} \\end{pmatrix}$$

$$\\textbf{C}\_{\\text{Complex}} = \\begin{pmatrix} 
-0.12 &  &  \\\\  
0 & 0 &  \\\\  
0 & 0 & 0.34 \\end{pmatrix}$$

$$\\textbf{C}\_{\\text{Simple}} = \\begin{pmatrix} 
-0.04 &  &  \\\\  
0.02 & -0.01 &  \\\\  
0.03 & -0.01 & -0.03 \\end{pmatrix}$$

Estimating selection on the egg parasitoid *Platygaster*
========================================================

``` r
# convert "gall_survival" to egg parasitoid survival. Note that both Iteomyia pupa and larva parasitoids result in 0. We do not change the name of the response variable to make clear that we are using the same statistical models as for Iteomyia.
eggegg_df <- mutate(gall_selection.df, gall_survival = ifelse(egg_parasitoid==1,1,0))

eggegg_model <- update(foodweb_model, data=eggegg_df)
```

    ## singular fit

``` r
biased_eggegg_df <- eggegg_df %>%
  group_by(Foodweb, Gall_Number) %>%
  mutate(mean_survival = mean(gall_survival)) %>%
  filter(mean_survival > 0, mean_survival < 1) %>%
  ungroup()

biased_eggegg_model <- update(biased_foodweb_model, data=biased_eggegg_df)
```

    ## singular fit

``` r
linear_eggegg_model <- update(linear_foodweb_model, data=eggegg_df)
```

    ## singular fit

``` r
knitr::kable(tidy_eggegg_alphas)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
type
</th>
<th style="text-align:left;">
term
</th>
<th style="text-align:right;">
conf.high
</th>
<th style="text-align:right;">
conf.low
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:left;">
P\_cutoff
</th>
<th style="text-align:left;">
coefficient\_type
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Complex
</td>
<td style="text-align:left;">
(Intercept)
</td>
<td style="text-align:right;">
0.3125329
</td>
<td style="text-align:right;">
0.0215894
</td>
<td style="text-align:right;">
0.1624476
</td>
<td style="text-align:left;">
-   </td>
    <td style="text-align:left;">
    Mean fitness
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    I(sc.Diam^2)
    </td>
    <td style="text-align:right;">
    0.1262707
    </td>
    <td style="text-align:right;">
    -0.6896657
    </td>
    <td style="text-align:right;">
    -0.2601131
    </td>
    <td style="text-align:left;">
    </td>
    <td style="text-align:left;">
    Quadratic
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    I(sc.log.Clutch^2)
    </td>
    <td style="text-align:right;">
    0.4398954
    </td>
    <td style="text-align:right;">
    -0.5026509
    </td>
    <td style="text-align:right;">
    -0.0246095
    </td>
    <td style="text-align:left;">
    </td>
    <td style="text-align:left;">
    Quadratic
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    I(sc.sqrt.Pref^2)
    </td>
    <td style="text-align:right;">
    -0.1702883
    </td>
    <td style="text-align:right;">
    -1.3456039
    </td>
    <td style="text-align:right;">
    -0.6319574
    </td>
    <td style="text-align:left;">
    -   </td>
        <td style="text-align:left;">
        Quadratic
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Complex
        </td>
        <td style="text-align:left;">
        sc.Diam
        </td>
        <td style="text-align:right;">
        -0.6917669
        </td>
        <td style="text-align:right;">
        -1.8948818
        </td>
        <td style="text-align:right;">
        -1.1889025
        </td>
        <td style="text-align:left;">
        -   </td>
            <td style="text-align:left;">
            Linear
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Complex
            </td>
            <td style="text-align:left;">
            sc.Diam:sc.log.Clutch
            </td>
            <td style="text-align:right;">
            0.8778495
            </td>
            <td style="text-align:right;">
            -0.2477370
            </td>
            <td style="text-align:right;">
            0.2746043
            </td>
            <td style="text-align:left;">
            </td>
            <td style="text-align:left;">
            Correlational
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Complex
            </td>
            <td style="text-align:left;">
            sc.Diam:sc.sqrt.Pref
            </td>
            <td style="text-align:right;">
            1.3849867
            </td>
            <td style="text-align:right;">
            0.0251265
            </td>
            <td style="text-align:right;">
            0.6396710
            </td>
            <td style="text-align:left;">
            -   </td>
                <td style="text-align:left;">
                Correlational
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Complex
                </td>
                <td style="text-align:left;">
                sc.log.Clutch
                </td>
                <td style="text-align:right;">
                1.5536961
                </td>
                <td style="text-align:right;">
                0.1726355
                </td>
                <td style="text-align:right;">
                0.8232199
                </td>
                <td style="text-align:left;">
                -   </td>
                    <td style="text-align:left;">
                    Linear
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Complex
                    </td>
                    <td style="text-align:left;">
                    sc.log.Clutch:sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    0.2479499
                    </td>
                    <td style="text-align:right;">
                    -1.2127422
                    </td>
                    <td style="text-align:right;">
                    -0.4277226
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Correlational
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Complex
                    </td>
                    <td style="text-align:left;">
                    sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    1.0274683
                    </td>
                    <td style="text-align:right;">
                    -0.4143469
                    </td>
                    <td style="text-align:right;">
                    0.2869481
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Linear
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    (Intercept)
                    </td>
                    <td style="text-align:right;">
                    0.3170565
                    </td>
                    <td style="text-align:right;">
                    -0.0864192
                    </td>
                    <td style="text-align:right;">
                    0.0955266
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Mean fitness
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    I(sc.Diam^2)
                    </td>
                    <td style="text-align:right;">
                    0.5386612
                    </td>
                    <td style="text-align:right;">
                    -0.6014360
                    </td>
                    <td style="text-align:right;">
                    -0.0141777
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Quadratic
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    I(sc.log.Clutch^2)
                    </td>
                    <td style="text-align:right;">
                    1.1058283
                    </td>
                    <td style="text-align:right;">
                    -0.3624756
                    </td>
                    <td style="text-align:right;">
                    0.3727100
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Quadratic
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    I(sc.sqrt.Pref^2)
                    </td>
                    <td style="text-align:right;">
                    1.4241463
                    </td>
                    <td style="text-align:right;">
                    -0.0796082
                    </td>
                    <td style="text-align:right;">
                    0.5724933
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Quadratic
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.Diam
                    </td>
                    <td style="text-align:right;">
                    0.4157411
                    </td>
                    <td style="text-align:right;">
                    -1.3115623
                    </td>
                    <td style="text-align:right;">
                    -0.3844165
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Linear
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.Diam:sc.log.Clutch
                    </td>
                    <td style="text-align:right;">
                    1.0029878
                    </td>
                    <td style="text-align:right;">
                    -0.7053471
                    </td>
                    <td style="text-align:right;">
                    0.1095376
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Correlational
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.Diam:sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    0.2043927
                    </td>
                    <td style="text-align:right;">
                    -1.6756672
                    </td>
                    <td style="text-align:right;">
                    -0.6495964
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Correlational
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.log.Clutch
                    </td>
                    <td style="text-align:right;">
                    1.1110362
                    </td>
                    <td style="text-align:right;">
                    -0.8246039
                    </td>
                    <td style="text-align:right;">
                    0.0918610
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Linear
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.log.Clutch:sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    1.4515568
                    </td>
                    <td style="text-align:right;">
                    -0.3859602
                    </td>
                    <td style="text-align:right;">
                    0.4874866
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Correlational
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    2.3103307
                    </td>
                    <td style="text-align:right;">
                    0.0406191
                    </td>
                    <td style="text-align:right;">
                    1.0270983
                    </td>
                    <td style="text-align:left;">
                    -   </td>
                        <td style="text-align:left;">
                        Linear
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        (Intercept)
                        </td>
                        <td style="text-align:right;">
                        0.4712831
                        </td>
                        <td style="text-align:right;">
                        0.0750668
                        </td>
                        <td style="text-align:right;">
                        0.2579742
                        </td>
                        <td style="text-align:left;">
                        -   </td>
                            <td style="text-align:left;">
                            Mean fitness
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Simple
                            </td>
                            <td style="text-align:left;">
                            I(sc.Diam^2)
                            </td>
                            <td style="text-align:right;">
                            0.1495383
                            </td>
                            <td style="text-align:right;">
                            -0.7239415
                            </td>
                            <td style="text-align:right;">
                            -0.2742908
                            </td>
                            <td style="text-align:left;">
                            </td>
                            <td style="text-align:left;">
                            Quadratic
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Simple
                            </td>
                            <td style="text-align:left;">
                            I(sc.log.Clutch^2)
                            </td>
                            <td style="text-align:right;">
                            0.9275679
                            </td>
                            <td style="text-align:right;">
                            -0.1696044
                            </td>
                            <td style="text-align:right;">
                            0.3481005
                            </td>
                            <td style="text-align:left;">
                            </td>
                            <td style="text-align:left;">
                            Quadratic
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Simple
                            </td>
                            <td style="text-align:left;">
                            I(sc.sqrt.Pref^2)
                            </td>
                            <td style="text-align:right;">
                            0.4963820
                            </td>
                            <td style="text-align:right;">
                            -0.5198053
                            </td>
                            <td style="text-align:right;">
                            -0.0594641
                            </td>
                            <td style="text-align:left;">
                            </td>
                            <td style="text-align:left;">
                            Quadratic
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Simple
                            </td>
                            <td style="text-align:left;">
                            sc.Diam
                            </td>
                            <td style="text-align:right;">
                            -0.8652251
                            </td>
                            <td style="text-align:right;">
                            -2.7473872
                            </td>
                            <td style="text-align:right;">
                            -1.5733190
                            </td>
                            <td style="text-align:left;">
                            -   </td>
                                <td style="text-align:left;">
                                Linear
                                </td>
                                </tr>
                                <tr>
                                <td style="text-align:left;">
                                Simple
                                </td>
                                <td style="text-align:left;">
                                sc.Diam:sc.log.Clutch
                                </td>
                                <td style="text-align:right;">
                                1.1077034
                                </td>
                                <td style="text-align:right;">
                                -0.1826248
                                </td>
                                <td style="text-align:right;">
                                0.3841420
                                </td>
                                <td style="text-align:left;">
                                </td>
                                <td style="text-align:left;">
                                Correlational
                                </td>
                                </tr>
                                <tr>
                                <td style="text-align:left;">
                                Simple
                                </td>
                                <td style="text-align:left;">
                                sc.Diam:sc.sqrt.Pref
                                </td>
                                <td style="text-align:right;">
                                0.5239441
                                </td>
                                <td style="text-align:right;">
                                -0.6447614
                                </td>
                                <td style="text-align:right;">
                                -0.0099254
                                </td>
                                <td style="text-align:left;">
                                </td>
                                <td style="text-align:left;">
                                Correlational
                                </td>
                                </tr>
                                <tr>
                                <td style="text-align:left;">
                                Simple
                                </td>
                                <td style="text-align:left;">
                                sc.log.Clutch
                                </td>
                                <td style="text-align:right;">
                                1.8463918
                                </td>
                                <td style="text-align:right;">
                                0.2331042
                                </td>
                                <td style="text-align:right;">
                                0.9150808
                                </td>
                                <td style="text-align:left;">
                                -   </td>
                                    <td style="text-align:left;">
                                    Linear
                                    </td>
                                    </tr>
                                    <tr>
                                    <td style="text-align:left;">
                                    Simple
                                    </td>
                                    <td style="text-align:left;">
                                    sc.log.Clutch:sc.sqrt.Pref
                                    </td>
                                    <td style="text-align:right;">
                                    0.6578362
                                    </td>
                                    <td style="text-align:right;">
                                    -0.4432248
                                    </td>
                                    <td style="text-align:right;">
                                    0.0597640
                                    </td>
                                    <td style="text-align:left;">
                                    </td>
                                    <td style="text-align:left;">
                                    Correlational
                                    </td>
                                    </tr>
                                    <tr>
                                    <td style="text-align:left;">
                                    Simple
                                    </td>
                                    <td style="text-align:left;">
                                    sc.sqrt.Pref
                                    </td>
                                    <td style="text-align:right;">
                                    2.3609067
                                    </td>
                                    <td style="text-align:right;">
                                    0.6663466
                                    </td>
                                    <td style="text-align:right;">
                                    1.3140464
                                    </td>
                                    <td style="text-align:left;">
                                    -   </td>
                                        <td style="text-align:left;">
                                        Linear
                                        </td>
                                        </tr>
                                        </tbody>
                                        </table>

``` r
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

![](reproduce_analyses_files/figure-markdown_github/Plot%20egg%20vs.%20egg%20ptoid%20coefficients-1.pdf)

``` r
knitr::kable(tidy_eggegg_grads)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
type
</th>
<th style="text-align:left;">
term
</th>
<th style="text-align:right;">
conf.high
</th>
<th style="text-align:right;">
conf.low
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:left;">
P\_cutoff
</th>
<th style="text-align:left;">
gradient\_type
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Complex
</td>
<td style="text-align:left;">
I(sc.Diam^2)
</td>
<td style="text-align:right;">
0.0977125
</td>
<td style="text-align:right;">
-0.5336865
</td>
<td style="text-align:right;">
-0.2012842
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
Quadratic
</td>
</tr>
<tr>
<td style="text-align:left;">
Complex
</td>
<td style="text-align:left;">
I(sc.log.Clutch^2)
</td>
<td style="text-align:right;">
0.3404058
</td>
<td style="text-align:right;">
-0.3889682
</td>
<td style="text-align:right;">
-0.0190437
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
Quadratic
</td>
</tr>
<tr>
<td style="text-align:left;">
Complex
</td>
<td style="text-align:left;">
I(sc.sqrt.Pref^2)
</td>
<td style="text-align:right;">
-0.1317748
</td>
<td style="text-align:right;">
-1.0412735
</td>
<td style="text-align:right;">
-0.4890299
</td>
<td style="text-align:left;">
-   </td>
    <td style="text-align:left;">
    Quadratic
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Complex
    </td>
    <td style="text-align:left;">
    sc.Diam
    </td>
    <td style="text-align:right;">
    -0.2676562
    </td>
    <td style="text-align:right;">
    -0.7331616
    </td>
    <td style="text-align:right;">
    -0.4600064
    </td>
    <td style="text-align:left;">
    -   </td>
        <td style="text-align:left;">
        Directional
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Complex
        </td>
        <td style="text-align:left;">
        sc.Diam:sc.log.Clutch
        </td>
        <td style="text-align:right;">
        0.3396547
        </td>
        <td style="text-align:right;">
        -0.0958536
        </td>
        <td style="text-align:right;">
        0.1062490
        </td>
        <td style="text-align:left;">
        </td>
        <td style="text-align:left;">
        Correlational
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Complex
        </td>
        <td style="text-align:left;">
        sc.Diam:sc.sqrt.Pref
        </td>
        <td style="text-align:right;">
        0.5358746
        </td>
        <td style="text-align:right;">
        0.0097218
        </td>
        <td style="text-align:right;">
        0.2474995
        </td>
        <td style="text-align:left;">
        -   </td>
            <td style="text-align:left;">
            Correlational
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Complex
            </td>
            <td style="text-align:left;">
            sc.log.Clutch
            </td>
            <td style="text-align:right;">
            0.6011511
            </td>
            <td style="text-align:right;">
            0.0667956
            </td>
            <td style="text-align:right;">
            0.3185176
            </td>
            <td style="text-align:left;">
            -   </td>
                <td style="text-align:left;">
                Directional
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Complex
                </td>
                <td style="text-align:left;">
                sc.log.Clutch:sc.sqrt.Pref
                </td>
                <td style="text-align:right;">
                0.0959360
                </td>
                <td style="text-align:right;">
                -0.4692303
                </td>
                <td style="text-align:right;">
                -0.1654931
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Correlational
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Complex
                </td>
                <td style="text-align:left;">
                sc.sqrt.Pref
                </td>
                <td style="text-align:right;">
                0.3975447
                </td>
                <td style="text-align:right;">
                -0.1603178
                </td>
                <td style="text-align:right;">
                0.1110250
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Directional
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                I(sc.Diam^2)
                </td>
                <td style="text-align:right;">
                0.4325082
                </td>
                <td style="text-align:right;">
                -0.3345143
                </td>
                <td style="text-align:right;">
                0.0512066
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Quadratic
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                I(sc.log.Clutch^2)
                </td>
                <td style="text-align:right;">
                0.6782940
                </td>
                <td style="text-align:right;">
                -0.2599217
                </td>
                <td style="text-align:right;">
                0.2095062
                </td>
                <td style="text-align:left;">
                </td>
                <td style="text-align:left;">
                Quadratic
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                Diff
                </td>
                <td style="text-align:left;">
                I(sc.sqrt.Pref^2)
                </td>
                <td style="text-align:right;">
                1.0651147
                </td>
                <td style="text-align:right;">
                0.0152821
                </td>
                <td style="text-align:right;">
                0.4564942
                </td>
                <td style="text-align:left;">
                -   </td>
                    <td style="text-align:left;">
                    Quadratic
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.Diam
                    </td>
                    <td style="text-align:right;">
                    0.3042380
                    </td>
                    <td style="text-align:right;">
                    -0.2472834
                    </td>
                    <td style="text-align:right;">
                    0.0295872
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Directional
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.Diam:sc.log.Clutch
                    </td>
                    <td style="text-align:right;">
                    0.2894447
                    </td>
                    <td style="text-align:right;">
                    -0.2743468
                    </td>
                    <td style="text-align:right;">
                    -0.0011578
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Correlational
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.Diam:sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    0.0451729
                    </td>
                    <td style="text-align:right;">
                    -0.5952869
                    </td>
                    <td style="text-align:right;">
                    -0.2502148
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Correlational
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.log.Clutch
                    </td>
                    <td style="text-align:right;">
                    0.2569284
                    </td>
                    <td style="text-align:right;">
                    -0.3819310
                    </td>
                    <td style="text-align:right;">
                    -0.0681753
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Directional
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.log.Clutch:sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    0.5239325
                    </td>
                    <td style="text-align:right;">
                    -0.1228667
                    </td>
                    <td style="text-align:right;">
                    0.1818429
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Correlational
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Diff
                    </td>
                    <td style="text-align:left;">
                    sc.sqrt.Pref
                    </td>
                    <td style="text-align:right;">
                    0.6414959
                    </td>
                    <td style="text-align:right;">
                    -0.0860389
                    </td>
                    <td style="text-align:right;">
                    0.2484639
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Directional
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Simple
                    </td>
                    <td style="text-align:left;">
                    I(sc.Diam^2)
                    </td>
                    <td style="text-align:right;">
                    0.0818196
                    </td>
                    <td style="text-align:right;">
                    -0.3961031
                    </td>
                    <td style="text-align:right;">
                    -0.1500777
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Quadratic
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Simple
                    </td>
                    <td style="text-align:left;">
                    I(sc.log.Clutch^2)
                    </td>
                    <td style="text-align:right;">
                    0.5075169
                    </td>
                    <td style="text-align:right;">
                    -0.0927987
                    </td>
                    <td style="text-align:right;">
                    0.1904625
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Quadratic
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Simple
                    </td>
                    <td style="text-align:left;">
                    I(sc.sqrt.Pref^2)
                    </td>
                    <td style="text-align:right;">
                    0.2715944
                    </td>
                    <td style="text-align:right;">
                    -0.2844104
                    </td>
                    <td style="text-align:right;">
                    -0.0325357
                    </td>
                    <td style="text-align:left;">
                    </td>
                    <td style="text-align:left;">
                    Quadratic
                    </td>
                    </tr>
                    <tr>
                    <td style="text-align:left;">
                    Simple
                    </td>
                    <td style="text-align:left;">
                    sc.Diam
                    </td>
                    <td style="text-align:right;">
                    -0.2367031
                    </td>
                    <td style="text-align:right;">
                    -0.7516137
                    </td>
                    <td style="text-align:right;">
                    -0.4304192
                    </td>
                    <td style="text-align:left;">
                    -   </td>
                        <td style="text-align:left;">
                        Directional
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        sc.Diam:sc.log.Clutch
                        </td>
                        <td style="text-align:right;">
                        0.3030388
                        </td>
                        <td style="text-align:right;">
                        -0.0499614
                        </td>
                        <td style="text-align:right;">
                        0.1050912
                        </td>
                        <td style="text-align:left;">
                        </td>
                        <td style="text-align:left;">
                        Correlational
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        sc.Diam:sc.sqrt.Pref
                        </td>
                        <td style="text-align:right;">
                        0.1433375
                        </td>
                        <td style="text-align:right;">
                        -0.1763899
                        </td>
                        <td style="text-align:right;">
                        -0.0027153
                        </td>
                        <td style="text-align:left;">
                        </td>
                        <td style="text-align:left;">
                        Correlational
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        Simple
                        </td>
                        <td style="text-align:left;">
                        sc.log.Clutch
                        </td>
                        <td style="text-align:right;">
                        0.5051248
                        </td>
                        <td style="text-align:right;">
                        0.0637712
                        </td>
                        <td style="text-align:right;">
                        0.2503423
                        </td>
                        <td style="text-align:left;">
                        -   </td>
                            <td style="text-align:left;">
                            Directional
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Simple
                            </td>
                            <td style="text-align:left;">
                            sc.log.Clutch:sc.sqrt.Pref
                            </td>
                            <td style="text-align:right;">
                            0.1799669
                            </td>
                            <td style="text-align:right;">
                            -0.1212548
                            </td>
                            <td style="text-align:right;">
                            0.0163499
                            </td>
                            <td style="text-align:left;">
                            </td>
                            <td style="text-align:left;">
                            Correlational
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Simple
                            </td>
                            <td style="text-align:left;">
                            sc.sqrt.Pref
                            </td>
                            <td style="text-align:right;">
                            0.6458827
                            </td>
                            <td style="text-align:right;">
                            0.1822951
                            </td>
                            <td style="text-align:right;">
                            0.3594889
                            </td>
                            <td style="text-align:left;">
                            -   </td>
                                <td style="text-align:left;">
                                Directional
                                </td>
                                </tr>
                                </tbody>
                                </table>

``` r
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

![](reproduce_analyses_files/figure-markdown_github/Plot%20egg%20vs.%20egg%20ptoid%20gradients-1.pdf)

Reproduce Figure 4 in main text
===============================

``` r
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

![](reproduce_analyses_files/figure-markdown_github/Egg-Parasitoid-Selection-1.pdf)

``` r
save_plot(filename = "selection_on_Platygaster.pdf", plot = eggptoid_selection_plot, base_height = 4) # base_height = 6 is to large for PDF of manuscript
```
