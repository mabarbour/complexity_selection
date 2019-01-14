Reproduce analyses reported in: *Phenotypic evolution is more constrained in simpler food webs*
================
Matthew A. Barbour
2019-01-14

Evaluating assumption of multivariate normality
===============================================

We used graphical checks to evaluate whether our transformations of trait values resulted in a multivariate normal distribution. Figure S shows that our transformations resulted in approximately normal distributions for each phenotypic trait. Note also that in the multivariate quantile-quantile (Q-Q) plot, most points fall along the expected line (fig. S), suggesting that our transformations provide a reasonable approximation of a multivariate normal distribution.

![Histograms of each phenotypic trait after transformation. The red line illustrates a normal distribution.](reproduce_analyses_files/figure-markdown_github/Univariate_histograms-1.pdf)

![Multivariate quantile-quantile (Q-Q) plot to assess deviations from multivariate normality (black line).](reproduce_analyses_files/figure-markdown_github/Multivariate_QQ-1.pdf)

Effect of food-web treatment on trait-fitness relationships and selection gradients
===================================================================================

We write the model in a way the independently estimates the effect of food-web treatment, each trait, and all two-way and three-way statistical interactions, on larval survival.

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

Note that the resulting estimates and confidence intervals are useful for determining whether trait-fitness relationships differ from zero, but not whether they differ between food-web treatments. For the later, we calculate the differences between each food-web treatment from the bootstrapped samples.

``` r
linear_foodweb_model <- glmer(
  gall_survival ~ 
    -1 + Foodweb + 
    Foodweb:(sc.Diam + sc.log.Clutch + sc.sqrt.Pref) +
    (1|Genotype/Plant_Position/Gall_Number),
  data = gall_selection.df,
  family = binomial(link = logit), control=glmerControl(optimizer = "bobyqa"))
```

To estimate biased selection on chamber diameter, we subset our data to only include multi-chambered galls where there was variability in larval survival. We then fit a reduced model to estimate the bias in the logistic regression coefficient of chamber diameter in each food web.

``` r
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

    ## singular fit

``` r
biased_foodweb_confint <- tidy(biased_foodweb_model, conf.int=TRUE) %>% filter(group=="fixed")
```

    ## Warning in bind_rows_(x, .id): binding factor and character vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

``` r
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
0.5564534
</td>
<td style="text-align:right;">
0.2767683
</td>
<td style="text-align:right;">
0.4149125
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
    0.5296517
    </td>
    <td style="text-align:right;">
    -0.0999299
    </td>
    <td style="text-align:right;">
    0.2205316
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
    0.2656116
    </td>
    <td style="text-align:right;">
    -0.4349161
    </td>
    <td style="text-align:right;">
    -0.0879843
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
    1.0476130
    </td>
    <td style="text-align:right;">
    0.1586764
    </td>
    <td style="text-align:right;">
    0.5641922
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
        1.6114227
        </td>
        <td style="text-align:right;">
        0.7362339
        </td>
        <td style="text-align:right;">
        1.1267669
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
            0.2466902
            </td>
            <td style="text-align:right;">
            -0.5503721
            </td>
            <td style="text-align:right;">
            -0.1421779
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
            0.0253427
            </td>
            <td style="text-align:right;">
            -0.9362864
            </td>
            <td style="text-align:right;">
            -0.4365104
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
            0.5698300
            </td>
            <td style="text-align:right;">
            -0.1167954
            </td>
            <td style="text-align:right;">
            0.2104417
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
            0.5716235
            </td>
            <td style="text-align:right;">
            -0.3921574
            </td>
            <td style="text-align:right;">
            0.0848484
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
            0.1570257
            </td>
            <td style="text-align:right;">
            -0.9788150
            </td>
            <td style="text-align:right;">
            -0.3981025
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
            0.4246751
            </td>
            <td style="text-align:right;">
            0.1234339
            </td>
            <td style="text-align:right;">
            0.2701913
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
                0.5490640
                </td>
                <td style="text-align:right;">
                -0.4129338
                </td>
                <td style="text-align:right;">
                0.0448577
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
                0.3492286
                </td>
                <td style="text-align:right;">
                -0.7405645
                </td>
                <td style="text-align:right;">
                -0.2022151
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
                0.0646544
                </td>
                <td style="text-align:right;">
                -1.2067323
                </td>
                <td style="text-align:right;">
                -0.5371566
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
                0.5311993
                </td>
                <td style="text-align:right;">
                -0.5404050
                </td>
                <td style="text-align:right;">
                -0.0227512
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
                0.3724560
                </td>
                <td style="text-align:right;">
                -0.7865475
                </td>
                <td style="text-align:right;">
                -0.2128666
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
                1.0334873
                </td>
                <td style="text-align:right;">
                -0.3039749
                </td>
                <td style="text-align:right;">
                0.3545167
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
                -0.1534212
                </td>
                <td style="text-align:right;">
                -1.2766208
                </td>
                <td style="text-align:right;">
                -0.6771629
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
                    0.5440235
                    </td>
                    <td style="text-align:right;">
                    -0.7081250
                    </td>
                    <td style="text-align:right;">
                    -0.0849055
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
                    0.3133863
                    </td>
                    <td style="text-align:right;">
                    -1.1184941
                    </td>
                    <td style="text-align:right;">
                    -0.4342433
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
                    0.8179927
                    </td>
                    <td style="text-align:right;">
                    0.5449427
                    </td>
                    <td style="text-align:right;">
                    0.6851037
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
                        0.6433138
                        </td>
                        <td style="text-align:right;">
                        -0.0646372
                        </td>
                        <td style="text-align:right;">
                        0.2653893
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
                        0.0850185
                        </td>
                        <td style="text-align:right;">
                        -0.7103899
                        </td>
                        <td style="text-align:right;">
                        -0.2901994
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
                        0.4452962
                        </td>
                        <td style="text-align:right;">
                        -0.3867337
                        </td>
                        <td style="text-align:right;">
                        0.0270357
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
                        1.6423291
                        </td>
                        <td style="text-align:right;">
                        0.6467462
                        </td>
                        <td style="text-align:right;">
                        1.1040157
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
                            0.1192674
                            </td>
                            <td style="text-align:right;">
                            -0.8036069
                            </td>
                            <td style="text-align:right;">
                            -0.3550445
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
                            0.3483629
                            </td>
                            <td style="text-align:right;">
                            -0.4983001
                            </td>
                            <td style="text-align:right;">
                            -0.0819937
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
                            -0.0869617
                            </td>
                            <td style="text-align:right;">
                            -0.9060484
                            </td>
                            <td style="text-align:right;">
                            -0.4667212
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
                                0.4108583
                                </td>
                                <td style="text-align:right;">
                                -0.3924871
                                </td>
                                <td style="text-align:right;">
                                -0.0000571
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
                                -0.3649635
                                </td>
                                <td style="text-align:right;">
                                -1.3398495
                                </td>
                                <td style="text-align:right;">
                                -0.8323458
                                </td>
                                <td style="text-align:left;">
                                -   </td>
                                    <td style="text-align:left;">
                                    Linear
                                    </td>
                                    </tr>
                                    </tbody>
                                    </table>

![](reproduce_analyses_files/figure-markdown_github/Plot%20food-web%20coefficients-1.pdf)

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
0.3174364
</td>
<td style="text-align:right;">
-0.0598910
</td>
<td style="text-align:right;">
0.1321713
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
0.1591891
</td>
<td style="text-align:right;">
-0.2606585
</td>
<td style="text-align:right;">
-0.0527317
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
0.6278664
</td>
<td style="text-align:right;">
0.0950996
</td>
<td style="text-align:right;">
0.3381376
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
    0.4828873
    </td>
    <td style="text-align:right;">
    0.2206237
    </td>
    <td style="text-align:right;">
    0.3376528
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
        0.0739245
        </td>
        <td style="text-align:right;">
        -0.1649274
        </td>
        <td style="text-align:right;">
        -0.0426058
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
        0.0075943
        </td>
        <td style="text-align:right;">
        -0.2805725
        </td>
        <td style="text-align:right;">
        -0.1308070
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
        0.1707582
        </td>
        <td style="text-align:right;">
        -0.0349995
        </td>
        <td style="text-align:right;">
        0.0630621
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
        0.1712957
        </td>
        <td style="text-align:right;">
        -0.1175159
        </td>
        <td style="text-align:right;">
        0.0254261
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
        0.0470551
        </td>
        <td style="text-align:right;">
        -0.2933168
        </td>
        <td style="text-align:right;">
        -0.1192975
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
        0.2056389
        </td>
        <td style="text-align:right;">
        -0.2544854
        </td>
        <td style="text-align:right;">
        -0.0323407
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
        0.2006797
        </td>
        <td style="text-align:right;">
        -0.3107148
        </td>
        <td style="text-align:right;">
        -0.0564317
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
        -0.0383193
        </td>
        <td style="text-align:right;">
        -0.6792632
        </td>
        <td style="text-align:right;">
        -0.3279676
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
            -0.0047796
            </td>
            <td style="text-align:right;">
            -0.2694282
            </td>
            <td style="text-align:right;">
            -0.1300059
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
                0.1198339
                </td>
                <td style="text-align:right;">
                -0.1669181
                </td>
                <td style="text-align:right;">
                -0.0241722
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
                0.2867498
                </td>
                <td style="text-align:right;">
                -0.0545360
                </td>
                <td style="text-align:right;">
                0.1153853
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
                -0.0259685
                </td>
                <td style="text-align:right;">
                -0.2919689
                </td>
                <td style="text-align:right;">
                -0.1508445
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
                    0.1356850
                    </td>
                    <td style="text-align:right;">
                    -0.1870362
                    </td>
                    <td style="text-align:right;">
                    -0.0254368
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
                    0.1559601
                    </td>
                    <td style="text-align:right;">
                    -0.2220200
                    </td>
                    <td style="text-align:right;">
                    -0.0372529
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
                    0.2419932
                    </td>
                    <td style="text-align:right;">
                    -0.0243144
                    </td>
                    <td style="text-align:right;">
                    0.0998306
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
                    0.0319811
                    </td>
                    <td style="text-align:right;">
                    -0.2672250
                    </td>
                    <td style="text-align:right;">
                    -0.1091633
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
                    0.1675056
                    </td>
                    <td style="text-align:right;">
                    -0.1454763
                    </td>
                    <td style="text-align:right;">
                    0.0101699
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
                    0.3088947
                    </td>
                    <td style="text-align:right;">
                    0.1216422
                    </td>
                    <td style="text-align:right;">
                    0.2076470
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
                        0.0224322
                        </td>
                        <td style="text-align:right;">
                        -0.1511451
                        </td>
                        <td style="text-align:right;">
                        -0.0667779
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
                        0.0655213
                        </td>
                        <td style="text-align:right;">
                        -0.0937219
                        </td>
                        <td style="text-align:right;">
                        -0.0154216
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
                        -0.0163560
                        </td>
                        <td style="text-align:right;">
                        -0.1704126
                        </td>
                        <td style="text-align:right;">
                        -0.0877825
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
                            0.0772756
                            </td>
                            <td style="text-align:right;">
                            -0.0738203
                            </td>
                            <td style="text-align:right;">
                            -0.0000107
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
                            -0.0686436
                            </td>
                            <td style="text-align:right;">
                            -0.2520034
                            </td>
                            <td style="text-align:right;">
                            -0.1565504
                            </td>
                            <td style="text-align:left;">
                            -   </td>
                                <td style="text-align:left;">
                                Directional
                                </td>
                                </tr>
                                </tbody>
                                </table>

![](reproduce_analyses_files/figure-markdown_github/Plot%20food-web%20gradients-1.pdf)

Partitioning the contribution of egg and larval parasitoids to selection gradients
==================================================================================

Our simple food-web treatment allows us to estimate the unique contribution of egg parasitoids to selection on *Iteomyia* traits. To estimate the unique contribution of larval parasitoids, we subset our data so that our complex food-web treatment only contained attack by larval parasitoids (and gall survival). We then fit the same models as previously, including one to estimate bias.

``` r
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

    ## singular fit

``` r
linear_egglarval_model <- update(linear_foodweb_model, data=egglarval_df)
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
0.8433768
</td>
<td style="text-align:right;">
0.4886051
</td>
<td style="text-align:right;">
0.6721235
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
    0.5358599
    </td>
    <td style="text-align:right;">
    -0.2086207
    </td>
    <td style="text-align:right;">
    0.1527953
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
    0.3691519
    </td>
    <td style="text-align:right;">
    -0.3576305
    </td>
    <td style="text-align:right;">
    0.0088975
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
    1.0759429
    </td>
    <td style="text-align:right;">
    -0.0623505
    </td>
    <td style="text-align:right;">
    0.4684399
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
    1.7572812
    </td>
    <td style="text-align:right;">
    0.6927731
    </td>
    <td style="text-align:right;">
    1.1443778
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
        0.5414879
        </td>
        <td style="text-align:right;">
        -0.3549199
        </td>
        <td style="text-align:right;">
        0.0989391
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
        0.4181587
        </td>
        <td style="text-align:right;">
        -0.8969578
        </td>
        <td style="text-align:right;">
        -0.2562654
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
        1.2218465
        </td>
        <td style="text-align:right;">
        0.1710512
        </td>
        <td style="text-align:right;">
        0.6714377
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
            0.2597239
            </td>
            <td style="text-align:right;">
            -0.7121276
            </td>
            <td style="text-align:right;">
            -0.2089646
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
            -0.2119552
            </td>
            <td style="text-align:right;">
            -1.7528262
            </td>
            <td style="text-align:right;">
            -0.9110117
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
                0.2227946
                </td>
                <td style="text-align:right;">
                -0.1888661
                </td>
                <td style="text-align:right;">
                0.0151660
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
                0.5724512
                </td>
                <td style="text-align:right;">
                -0.3637220
                </td>
                <td style="text-align:right;">
                0.1002200
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
                0.2006751
                </td>
                <td style="text-align:right;">
                -0.8814329
                </td>
                <td style="text-align:right;">
                -0.3037826
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
                0.1535886
                </td>
                <td style="text-align:right;">
                -1.2993352
                </td>
                <td style="text-align:right;">
                -0.5390193
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
                0.5437150
                </td>
                <td style="text-align:right;">
                -0.6700476
                </td>
                <td style="text-align:right;">
                -0.0607140
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
                0.1513264
                </td>
                <td style="text-align:right;">
                -1.0683525
                </td>
                <td style="text-align:right;">
                -0.4588715
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
                0.8764318
                </td>
                <td style="text-align:right;">
                -0.5751819
                </td>
                <td style="text-align:right;">
                0.1579216
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
                -0.7604633
                </td>
                <td style="text-align:right;">
                -2.3185112
                </td>
                <td style="text-align:right;">
                -1.4820541
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
                    0.8788350
                    </td>
                    <td style="text-align:right;">
                    -0.3537164
                    </td>
                    <td style="text-align:right;">
                    0.2379177
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
                    0.7877674
                    </td>
                    <td style="text-align:right;">
                    -1.0024104
                    </td>
                    <td style="text-align:right;">
                    -0.0911791
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
                    0.8468822
                    </td>
                    <td style="text-align:right;">
                    0.4887100
                    </td>
                    <td style="text-align:right;">
                    0.6872896
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
                        0.5728484
                        </td>
                        <td style="text-align:right;">
                        -0.0454721
                        </td>
                        <td style="text-align:right;">
                        0.2530153
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
                        0.0530138
                        </td>
                        <td style="text-align:right;">
                        -0.6451803
                        </td>
                        <td style="text-align:right;">
                        -0.2948851
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
                        0.3515456
                        </td>
                        <td style="text-align:right;">
                        -0.4988641
                        </td>
                        <td style="text-align:right;">
                        -0.0705794
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
                        1.5959742
                        </td>
                        <td style="text-align:right;">
                        0.6500136
                        </td>
                        <td style="text-align:right;">
                        1.0836639
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
                            0.0651246
                            </td>
                            <td style="text-align:right;">
                            -0.7926376
                            </td>
                            <td style="text-align:right;">
                            -0.3599323
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
                            0.2928635
                            </td>
                            <td style="text-align:right;">
                            -0.4950078
                            </td>
                            <td style="text-align:right;">
                            -0.0983438
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
                            -0.3523461
                            </td>
                            <td style="text-align:right;">
                            -1.3631444
                            </td>
                            <td style="text-align:right;">
                            -0.8106165
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
                                0.3685492
                                </td>
                                <td style="text-align:right;">
                                -0.3284094
                                </td>
                                <td style="text-align:right;">
                                0.0289531
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
                                -0.4843481
                                </td>
                                <td style="text-align:right;">
                                -1.6486545
                                </td>
                                <td style="text-align:right;">
                                -1.0021908
                                </td>
                                <td style="text-align:left;">
                                -   </td>
                                    <td style="text-align:left;">
                                    Linear
                                    </td>
                                    </tr>
                                    </tbody>
                                    </table>

![](reproduce_analyses_files/figure-markdown_github/Plot%20egglarval%20coefficients-1.pdf)

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
0.2075496
</td>
<td style="text-align:right;">
-0.0808031
</td>
<td style="text-align:right;">
0.0591808
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
0.1429802
</td>
<td style="text-align:right;">
-0.1385177
</td>
<td style="text-align:right;">
0.0034462
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
0.4167350
</td>
<td style="text-align:right;">
-0.0241496
</td>
<td style="text-align:right;">
0.1814365
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
0.3403157
</td>
<td style="text-align:right;">
0.1341627
</td>
<td style="text-align:right;">
0.2216206
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
    0.1048647
    </td>
    <td style="text-align:right;">
    -0.0687339
    </td>
    <td style="text-align:right;">
    0.0191606
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
    0.0809808
    </td>
    <td style="text-align:right;">
    -0.1737052
    </td>
    <td style="text-align:right;">
    -0.0496284
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
    0.2366232
    </td>
    <td style="text-align:right;">
    0.0331258
    </td>
    <td style="text-align:right;">
    0.1300309
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
        0.0502982
        </td>
        <td style="text-align:right;">
        -0.1379109
        </td>
        <td style="text-align:right;">
        -0.0404682
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
        -0.0410473
        </td>
        <td style="text-align:right;">
        -0.3394529
        </td>
        <td style="text-align:right;">
        -0.1764268
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
            0.2262719
            </td>
            <td style="text-align:right;">
            -0.1401120
            </td>
            <td style="text-align:right;">
            0.0409310
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
            0.0745750
            </td>
            <td style="text-align:right;">
            -0.3468058
            </td>
            <td style="text-align:right;">
            -0.1201248
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
            0.0647610
            </td>
            <td style="text-align:right;">
            -0.5053067
            </td>
            <td style="text-align:right;">
            -0.2093630
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
            0.1113543
            </td>
            <td style="text-align:right;">
            -0.1259153
            </td>
            <td style="text-align:right;">
            -0.0072313
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
            0.0282850
            </td>
            <td style="text-align:right;">
            -0.2090584
            </td>
            <td style="text-align:right;">
            -0.0903687
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
            0.1713471
            </td>
            <td style="text-align:right;">
            -0.1134273
            </td>
            <td style="text-align:right;">
            0.0301724
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
            -0.1487636
            </td>
            <td style="text-align:right;">
            -0.4545508
            </td>
            <td style="text-align:right;">
            -0.2904012
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
                0.1713788
                </td>
                <td style="text-align:right;">
                -0.0691109
                </td>
                <td style="text-align:right;">
                0.0461962
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
                0.1489210
                </td>
                <td style="text-align:right;">
                -0.2007294
                </td>
                <td style="text-align:right;">
                -0.0218441
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
                0.2266618
                </td>
                <td style="text-align:right;">
                -0.0179922
                </td>
                <td style="text-align:right;">
                0.1001118
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
                0.0209762
                </td>
                <td style="text-align:right;">
                -0.2552817
                </td>
                <td style="text-align:right;">
                -0.1166786
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
                0.1390978
                </td>
                <td style="text-align:right;">
                -0.1973880
                </td>
                <td style="text-align:right;">
                -0.0279265
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
                0.3157435
                </td>
                <td style="text-align:right;">
                0.1285971
                </td>
                <td style="text-align:right;">
                0.2143893
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
                    0.0128841
                    </td>
                    <td style="text-align:right;">
                    -0.1568134
                    </td>
                    <td style="text-align:right;">
                    -0.0712081
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
                    0.0579394
                    </td>
                    <td style="text-align:right;">
                    -0.0979311
                    </td>
                    <td style="text-align:right;">
                    -0.0194561
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
                    -0.0697073
                    </td>
                    <td style="text-align:right;">
                    -0.2696810
                    </td>
                    <td style="text-align:right;">
                    -0.1603703
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
                        0.0729128
                        </td>
                        <td style="text-align:right;">
                        -0.0649717
                        </td>
                        <td style="text-align:right;">
                        0.0057280
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
                        -0.0958222
                        </td>
                        <td style="text-align:right;">
                        -0.3261657
                        </td>
                        <td style="text-align:right;">
                        -0.1982709
                        </td>
                        <td style="text-align:left;">
                        -   </td>
                            <td style="text-align:left;">
                            Directional
                            </td>
                            </tr>
                            </tbody>
                            </table>

![](reproduce_analyses_files/figure-markdown_github/Plot%20egg%20vs.%20larval%20ptoid%20gradients-1.pdf)

We combine our estimates of selection gradients for each food-web treatment as well as the contribution of larval parasitoids to selection in the complex food web.

![**Selection gradients**. Estimates of standardized selection gradients in complex (orange) and simple (blue) food webs. The contribution of larval parasitoids (yellow) was estimated with a subset of the complex food-web data that only contained attacks from larval parasitoids (and gall survival). Points and lines correspond to estimates of the mean and 95% confidence intervals, respectively. Overlapping confidence intervals with zero (dotted line) indicate no strong evidence of selection.](reproduce_analyses_files/figure-markdown_github/Figure-Selection-Gradients-1.pdf)

Partitioning the components of selection gradients
==================================================

Selection gradients are influenced by both trait-fitness relationships and population mean fitness. Here, we partition selection gradients into these underlying components.

![**Partitioning Selection Gradients.**](reproduce_analyses_files/figure-markdown_github/Plot%20foodwebegglarval%20coefficients-1.pdf)

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
0.3041841
</td>
<td style="text-align:right;">
0.0211660
</td>
<td style="text-align:right;">
0.1629208
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
    0.1035432
    </td>
    <td style="text-align:right;">
    -0.6904251
    </td>
    <td style="text-align:right;">
    -0.2625753
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
    0.4504659
    </td>
    <td style="text-align:right;">
    -0.5156402
    </td>
    <td style="text-align:right;">
    -0.0317396
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
    -0.1504103
    </td>
    <td style="text-align:right;">
    -1.3547025
    </td>
    <td style="text-align:right;">
    -0.6267839
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
        -0.6436999
        </td>
        <td style="text-align:right;">
        -1.8222668
        </td>
        <td style="text-align:right;">
        -1.1711583
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
            0.8595340
            </td>
            <td style="text-align:right;">
            -0.2258934
            </td>
            <td style="text-align:right;">
            0.2810014
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
            1.3104580
            </td>
            <td style="text-align:right;">
            0.0092513
            </td>
            <td style="text-align:right;">
            0.6082345
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
                1.6272424
                </td>
                <td style="text-align:right;">
                0.1753339
                </td>
                <td style="text-align:right;">
                0.8242237
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
                    0.2552999
                    </td>
                    <td style="text-align:right;">
                    -1.1274528
                    </td>
                    <td style="text-align:right;">
                    -0.3965201
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
                    0.9651409
                    </td>
                    <td style="text-align:right;">
                    -0.4409892
                    </td>
                    <td style="text-align:right;">
                    0.2572587
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
                    0.3367869
                    </td>
                    <td style="text-align:right;">
                    -0.0977451
                    </td>
                    <td style="text-align:right;">
                    0.0968295
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
                    0.6116945
                    </td>
                    <td style="text-align:right;">
                    -0.5722614
                    </td>
                    <td style="text-align:right;">
                    -0.0059298
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
                    1.1284363
                    </td>
                    <td style="text-align:right;">
                    -0.3402574
                    </td>
                    <td style="text-align:right;">
                    0.3760892
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
                    1.4994696
                    </td>
                    <td style="text-align:right;">
                    -0.1439261
                    </td>
                    <td style="text-align:right;">
                    0.5549636
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
                    0.4083808
                    </td>
                    <td style="text-align:right;">
                    -1.4859239
                    </td>
                    <td style="text-align:right;">
                    -0.4015842
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
                    0.9293746
                    </td>
                    <td style="text-align:right;">
                    -0.6783191
                    </td>
                    <td style="text-align:right;">
                    0.1010129
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
                    0.1958785
                    </td>
                    <td style="text-align:right;">
                    -1.5759066
                    </td>
                    <td style="text-align:right;">
                    -0.6270852
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
                    1.1913401
                    </td>
                    <td style="text-align:right;">
                    -0.8539212
                    </td>
                    <td style="text-align:right;">
                    0.0832914
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
                    1.4070758
                    </td>
                    <td style="text-align:right;">
                    -0.3479323
                    </td>
                    <td style="text-align:right;">
                    0.4697833
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
                    2.3242933
                    </td>
                    <td style="text-align:right;">
                    0.1428940
                    </td>
                    <td style="text-align:right;">
                    1.0495315
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
                        0.4624934
                        </td>
                        <td style="text-align:right;">
                        0.0647210
                        </td>
                        <td style="text-align:right;">
                        0.2597503
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
                            0.1565899
                            </td>
                            <td style="text-align:right;">
                            -0.6975499
                            </td>
                            <td style="text-align:right;">
                            -0.2685050
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
                            0.9194973
                            </td>
                            <td style="text-align:right;">
                            -0.2034549
                            </td>
                            <td style="text-align:right;">
                            0.3443496
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
                            0.4138438
                            </td>
                            <td style="text-align:right;">
                            -0.5510831
                            </td>
                            <td style="text-align:right;">
                            -0.0718203
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
                            -0.8820073
                            </td>
                            <td style="text-align:right;">
                            -2.7385212
                            </td>
                            <td style="text-align:right;">
                            -1.5727425
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
                                1.0363396
                                </td>
                                <td style="text-align:right;">
                                -0.1786224
                                </td>
                                <td style="text-align:right;">
                                0.3820143
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
                                0.5219210
                                </td>
                                <td style="text-align:right;">
                                -0.6277397
                                </td>
                                <td style="text-align:right;">
                                -0.0188507
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
                                1.8220204
                                </td>
                                <td style="text-align:right;">
                                0.1972724
                                </td>
                                <td style="text-align:right;">
                                0.9075150
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
                                    0.6519129
                                    </td>
                                    <td style="text-align:right;">
                                    -0.4240386
                                    </td>
                                    <td style="text-align:right;">
                                    0.0732631
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
                                    2.2940357
                                    </td>
                                    <td style="text-align:right;">
                                    0.6586756
                                    </td>
                                    <td style="text-align:right;">
                                    1.3067902
                                    </td>
                                    <td style="text-align:left;">
                                    -   </td>
                                        <td style="text-align:left;">
                                        Linear
                                        </td>
                                        </tr>
                                        </tbody>
                                        </table>

![](reproduce_analyses_files/figure-markdown_github/Plot%20egg%20vs.%20egg%20ptoid%20coefficients-1.pdf)

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
0.0801252
</td>
<td style="text-align:right;">
-0.5342741
</td>
<td style="text-align:right;">
-0.2031896
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
0.3485856
</td>
<td style="text-align:right;">
-0.3990197
</td>
<td style="text-align:right;">
-0.0245611
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
-0.1163925
</td>
<td style="text-align:right;">
-1.0483143
</td>
<td style="text-align:right;">
-0.4850264
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
    -0.2490583
    </td>
    <td style="text-align:right;">
    -0.7050656
    </td>
    <td style="text-align:right;">
    -0.4531408
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
        0.3325681
        </td>
        <td style="text-align:right;">
        -0.0874019
        </td>
        <td style="text-align:right;">
        0.1087242
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
        0.5070382
        </td>
        <td style="text-align:right;">
        0.0035795
        </td>
        <td style="text-align:right;">
        0.2353361
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
            0.6296074
            </td>
            <td style="text-align:right;">
            0.0678396
            </td>
            <td style="text-align:right;">
            0.3189060
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
                0.0987798
                </td>
                <td style="text-align:right;">
                -0.4362304
                </td>
                <td style="text-align:right;">
                -0.1534203
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
                0.3734292
                </td>
                <td style="text-align:right;">
                -0.1706261
                </td>
                <td style="text-align:right;">
                0.0995377
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
                0.4797535
                </td>
                <td style="text-align:right;">
                -0.3226508
                </td>
                <td style="text-align:right;">
                0.0562776
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
                0.7166843
                </td>
                <td style="text-align:right;">
                -0.2485213
                </td>
                <td style="text-align:right;">
                0.2129713
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
                1.0920899
                </td>
                <td style="text-align:right;">
                -0.0098946
                </td>
                <td style="text-align:right;">
                0.4457301
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
                0.2797324
                </td>
                <td style="text-align:right;">
                -0.3045078
                </td>
                <td style="text-align:right;">
                0.0228794
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
                0.2672369
                </td>
                <td style="text-align:right;">
                -0.2661632
                </td>
                <td style="text-align:right;">
                -0.0042150
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
                0.0396831
                </td>
                <td style="text-align:right;">
                -0.5725511
                </td>
                <td style="text-align:right;">
                -0.2404932
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
                0.2658295
                </td>
                <td style="text-align:right;">
                -0.3802372
                </td>
                <td style="text-align:right;">
                -0.0706335
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
                0.4961143
                </td>
                <td style="text-align:right;">
                -0.1067624
                </td>
                <td style="text-align:right;">
                0.1734632
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
                0.6761936
                </td>
                <td style="text-align:right;">
                -0.0585113
                </td>
                <td style="text-align:right;">
                0.2579661
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
                0.0856778
                </td>
                <td style="text-align:right;">
                -0.3816630
                </td>
                <td style="text-align:right;">
                -0.1469120
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
                0.5031011
                </td>
                <td style="text-align:right;">
                -0.1113199
                </td>
                <td style="text-align:right;">
                0.1884102
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
                0.2264338
                </td>
                <td style="text-align:right;">
                -0.3015240
                </td>
                <td style="text-align:right;">
                -0.0392963
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
                -0.2412942
                </td>
                <td style="text-align:right;">
                -0.7491882
                </td>
                <td style="text-align:right;">
                -0.4302614
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
                    0.2835156
                    </td>
                    <td style="text-align:right;">
                    -0.0488665
                    </td>
                    <td style="text-align:right;">
                    0.1045092
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
                    0.1427840
                    </td>
                    <td style="text-align:right;">
                    -0.1717332
                    </td>
                    <td style="text-align:right;">
                    -0.0051571
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
                    0.4984574
                    </td>
                    <td style="text-align:right;">
                    0.0539686
                    </td>
                    <td style="text-align:right;">
                    0.2482725
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
                        0.1783464
                        </td>
                        <td style="text-align:right;">
                        -0.1160059
                        </td>
                        <td style="text-align:right;">
                        0.0200429
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
                        0.6275885
                        </td>
                        <td style="text-align:right;">
                        0.1801965
                        </td>
                        <td style="text-align:right;">
                        0.3575038
                        </td>
                        <td style="text-align:left;">
                        -   </td>
                            <td style="text-align:left;">
                            Directional
                            </td>
                            </tr>
                            </tbody>
                            </table>

![](reproduce_analyses_files/figure-markdown_github/Plot%20egg%20vs.%20egg%20ptoid%20gradients-1.pdf)

Reproduce Figure 4 - Univariate adaptive landscape
==================================================

![Selection gradients acting on gall traits in complex vs. simple food webs. Each panel corresponds to a different gall trait: gall diameter (A); clutch size (B); and Oviposition preference (C). Solid lines represent the estimated gradients in complex (orange) and simple (blue) food webs. Transparent lines represent bootstrapped replicates (n=100) to show the uncertainty in estimated gradients. Note that only 100 bootstraps are displayed here, but that inferences are based on 1,000 bootstrapped samples.](reproduce_analyses_files/figure-markdown_github/Univariate-Fitness-Landscapes-1.pdf)

![](reproduce_analyses_files/figure-markdown_github/Egg-Parasitoid-Selection-1.pdf)

![](reproduce_analyses_files/figure-markdown_github/Egg-Parasitoid-Selection-linear-1.pdf)

Multivariate fitness landscapes
===============================

Reproduce Figure 3 in main manuscript.

![Fitness landscapes of gall traits in complex vs. simple food webs. Each panel corresponds to a different combination of traits: clutch size and gall diameter (A); clutch size and Oviposition preference (B); Oviposition preference and gall diameter (C). Note that traits for all plots range 1 SD below and above the mean (=0).](reproduce_analyses_files/figure-markdown_github/Figure-Multivariate-Landscapes-1.pdf)
