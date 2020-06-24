Travaglio et al., medRxiv 2020 (v04)
================
Yizhou Yu and Marco Travaglio, MRC Toxicology Unit, University of
Cambridge, (<yzy21@mrc-tox.cam.ac.uk>)
26/05/2020

# General information

This is a data analysis pipeline linked to a manuscript deposited in
medRxiv titled **Links between air pollution and COVID-19 in England**
(link: <https://www.medrxiv.org/content/10.1101/2020.04.16.20067405v2>).
This pipeline relates to data included in **Figure 3** to **Figure 5**
of the manuscript as well as Supplementary tables included in the
supplementary materials.

The input files for this analysis pipeline are on the master branch of
this GitHub page (link: <https://github.com/M1gus/AirPollutionCOVID19>)

The main aims of the workflow presented here were as follows: <br> 1.
Determine a relationship between air pollutants and COVID-19-associated
deaths/cases in England<br> 2. Investigate whether any relationship
between air pollution and COVID-19 remains significant in the presence
of counfounding factors at the regional and subregional level<br> 3.
Determine the effect of air pollutants on infectivity at individual
levels. <br> 4. Determine the main contributors of air pollution at the
subregional level. <br>

Only the UK Biobank data is not available in the repository as they
require separate application. Please visit ukbiobank.ac.uk for more
information. A detailed list of the variables used in the UK Biobank
analysis is available here:

#### Supplementary Table 1. Variables from the UK Biobank.

``` r
library(stargazer)
stargazer(read.csv("data_v4/suppl_table_1.csv"), summary=FALSE, rownames=FALSE, type = "html")
```

<table style="text-align:center">

<tr>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Variable.name

</td>

<td>

UK.Biobank.ID

</td>

<td>

Description.of.the.variable

</td>

</tr>

<tr>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

x\_coord

</td>

<td>

20074

</td>

<td>

Home location at assessment - east co-ordinate

</td>

</tr>

<tr>

<td style="text-align:left">

y\_coord

</td>

<td>

20075

</td>

<td>

Home location at assessment - north co-ordinate

</td>

</tr>

<tr>

<td style="text-align:left">

townsend

</td>

<td>

189

</td>

<td>

Townsend deprivation index at recruitment

</td>

</tr>

<tr>

<td style="text-align:left">

Age

</td>

<td>

34

</td>

<td>

2020 subtracted by year of birth

</td>

</tr>

<tr>

<td style="text-align:left">

whr

</td>

<td>

48, 49

</td>

<td>

Waist circumference divided by hip circumference

</td>

</tr>

<tr>

<td style="text-align:left">

highBP

</td>

<td>

4079, 4080

</td>

<td>

High Blood pressure\*

</td>

</tr>

<tr>

<td style="text-align:left">

WP\_dusty

</td>

<td>

22609

</td>

<td>

Workplace very dusty\*

</td>

</tr>

<tr>

<td style="text-align:left">

WP\_chemicals

</td>

<td>

22610

</td>

<td>

Workplace full of chemical or other fumes\*

</td>

</tr>

<tr>

<td style="text-align:left">

WP\_cig

</td>

<td>

22611

</td>

<td>

Workplace had a lot of cigarette smoke from other people smoking\*

</td>

</tr>

<tr>

<td style="text-align:left">

WP\_diesel

</td>

<td>

22615

</td>

<td>

Workplace had a lot of diesel exhaust\*

</td>

</tr>

<tr>

<td style="text-align:left">

breathing

</td>

<td>

22616

</td>

<td>

Breathing problems during period of job\*

</td>

</tr>

<tr>

<td style="text-align:left">

whistling

</td>

<td>

2316

</td>

<td>

Wheeze or whistling in the chest in last year

</td>

</tr>

<tr>

<td style="text-align:left">

sex

</td>

<td>

31

</td>

<td>

Sex

</td>

</tr>

<tr>

<td style="text-align:left">

fev

</td>

<td>

20153

</td>

<td>

Forced expiratory volume in 1-second (FEV1), predicted\*

</td>

</tr>

<tr>

<td style="text-align:left">

copd

</td>

<td>

22130

</td>

<td>

Doctor diagnosed COPD (chronic obstructive pulmonary disease)\*

</td>

</tr>

<tr>

<td style="text-align:left">

diabetes

</td>

<td>

2443

</td>

<td>

Diabetes diagnosed by doctor

</td>

</tr>

<tr>

<td style="text-align:left">

n\_cancers

</td>

<td>

134

</td>

<td>

Number of self-reported cancers

</td>

</tr>

<tr>

<td style="text-align:left">

smoking

</td>

<td>

20116

</td>

<td>

Smoking status\*

</td>

</tr>

<tr>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

</table>

<br>\*These variables contained too few data and were excluded from
analysis.<br>The detailed information for every variable can be
retrieved by searching for the UK Biobank ID at ukbiobank.ac.uk.

## Aim 1: Preliminary analysis

Deaths/cases \~ air pollutants + population <br> Then, Model comparisons
<br>

### Load packages

``` r
library(MASS)
library(stargazer)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(purrr)
library(AER)
library(reshape)
library(rgdal) # for spTransform, requires sudo apt install libgdal-dev
require(stringr)
library(sp)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
```

Please note: <br> the following packages may need to be installed:
liblapack-dev<br> liblapack3<br> libopenblas-base<br>
libopenblas-dev<br> texlive-fonts-recommended<br> lmodern<br> can be
done with apt-get on ubuntu

### View the data

Read and verify the data:

``` r
preL_dt = read.csv("data_v4/26-4-2020_yyAIR_COVID_PRE_LD_dt.csv")[,-1]
head(preL_dt)
```

    ##                     Region cases_preL deaths_preL Date_cases
    ## 1          East Of England       5356         746         NA
    ## 2                   London      16913        2120         NA
    ## 3                 Midlands      10501        1491         NA
    ## 4 North East And Yorkshire       8004         893         NA
    ## 5               North West       9394         847         NA
    ## 6               South East       8547         783         NA
    ##   Population_size_2018 Average_Pop_density_personkm2 Cases Deaths NO.levels
    ## 1              6201214                         324.0  6499   1448  9.502135
    ## 2              8908081                        5666.0 19511   3522 25.193133
    ## 3             10704906                         380.5 14844   2684 14.529437
    ## 4              8137524                         333.0 10633   1641 16.501209
    ## 5              7292093                         517.0 12093   1801  8.581661
    ## 6              9133625                         479.0 10919   1430 10.646127
    ##   NO2.levels O3.levels
    ## 1   19.55372  54.36748
    ## 2   38.51252  36.91391
    ## 3   24.11985  47.69889
    ## 4   25.15139  45.29534
    ## 5   20.06274  48.95081
    ## 6   20.47013  52.30382

Quickly visualise the distribution of each variable:

``` r
require (reshape)
#normalise per column & reshape the data for plotting
preL_dt_scale = as.data.frame(sapply(preL_dt[,-c(1,4)], scale))
preL_dt_scale$Region = preL_dt$Region
preL_dt_long = melt(preL_dt_scale, in.vars = "Region")

ggplot(preL_dt_long, aes (value)) +
    geom_density() +
  geom_histogram() +
    facet_wrap(~variable) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

### Models

#### Analyse the death data:

``` r
summary(full_lm <- lm(data = preL_dt, deaths_preL ~ Average_Pop_density_personkm2 + NO.levels + NO2.levels + O3.levels))
```

    ## 
    ## Call:
    ## lm(formula = deaths_preL ~ Average_Pop_density_personkm2 + NO.levels + 
    ##     NO2.levels + O3.levels, data = preL_dt)
    ## 
    ## Residuals:
    ##        1        2        3        4        5        6        7 
    ## -101.701   -1.428  192.949 -187.616  -13.232   44.475   66.554 
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                   -1.235e+04  4.844e+03  -2.549   0.1256  
    ## Average_Pop_density_personkm2 -6.310e-01  2.544e-01  -2.480   0.1313  
    ## NO.levels                     -3.326e+02  1.124e+02  -2.958   0.0978 .
    ## NO2.levels                     6.015e+02  1.782e+02   3.375   0.0777 .
    ## O3.levels                      8.824e+01  5.409e+01   1.631   0.2444  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 211.4 on 2 degrees of freedom
    ## Multiple R-squared:  0.956,  Adjusted R-squared:  0.8681 
    ## F-statistic: 10.87 on 4 and 2 DF,  p-value: 0.08599

The model is overdispersed - use negative binomial or poisson
regressions. <br> The poisson and Negative Binomial models work better
with count data. Changing these will not affect the OLS results.

``` r
preL_dt$NO.levels = round(preL_dt$NO.levels * 1000)
preL_dt$NO2.levels = round(preL_dt$NO2.levels * 1000)
preL_dt$O3.levels = round(preL_dt$O3.levels * 1000)
preL_dt$Average_Pop_density_personkm2 = round(preL_dt$Average_Pop_density_personkm2)
```

It is important to ensure that they are integers:

``` r
preL_dt_int = as.data.frame(sapply(preL_dt[,2:ncol(preL_dt)], as.integer))
preL_dt_int$Region = preL_dt$Region
preL_dt_int
```

    ##   cases_preL deaths_preL Date_cases Population_size_2018
    ## 1       5356         746         NA              6201214
    ## 2      16913        2120         NA              8908081
    ## 3      10501        1491         NA             10704906
    ## 4       8004         893         NA              8137524
    ## 5       9394         847         NA              7292093
    ## 6       8547         783         NA              9133625
    ## 7       2898         368         NA              5599735
    ##   Average_Pop_density_personkm2 Cases Deaths NO.levels NO2.levels O3.levels
    ## 1                           324  6499   1448      9502      19554     54367
    ## 2                          5666 19511   3522     25193      38513     36914
    ## 3                           380 14844   2684     14529      24120     47699
    ## 4                           333 10633   1641     16501      25151     45295
    ## 5                           517 12093   1801      8582      20063     48951
    ## 6                           479 10919   1430     10646      20470     52304
    ## 7                           235  3913    608     14047      21991     48054
    ##                     Region
    ## 1          East Of England
    ## 2                   London
    ## 3                 Midlands
    ## 4 North East And Yorkshire
    ## 5               North West
    ## 6               South East
    ## 7               South West

``` r
### Number of death
summary(full_lm.nb <- glm.nb(data = preL_dt_int, deaths_preL ~ Average_Pop_density_personkm2 + NO.levels + NO2.levels + O3.levels))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths_preL ~ Average_Pop_density_personkm2 + 
    ##     NO.levels + NO2.levels + O3.levels, data = preL_dt_int, init.theta = 196.4212515, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##        1         2         3         4         5         6         7  
    ## -1.50927  -0.04204   1.13015  -1.12365  -0.18404   1.44286   0.13514  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                   -1.090e+01  1.890e+00  -5.767 8.07e-09 ***
    ## Average_Pop_density_personkm2 -9.331e-04  9.829e-05  -9.494  < 2e-16 ***
    ## NO.levels                     -4.728e-04  4.487e-05 -10.539  < 2e-16 ***
    ## NO2.levels                     8.135e-04  7.086e-05  11.482  < 2e-16 ***
    ## O3.levels                      1.200e-04  2.068e-05   5.802 6.56e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(196.4213) family taken to be 1)
    ## 
    ##     Null deviance: 294.5029  on 6  degrees of freedom
    ## Residual deviance:   6.9535  on 2  degrees of freedom
    ## AIC: 91.745
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  196 
    ##           Std. Err.:  129 
    ## 
    ##  2 x log-likelihood:  -79.745

``` r
summary(full_lm.p <- glm(data = preL_dt_int, deaths_preL ~ Average_Pop_density_personkm2 + NO.levels + NO2.levels + O3.levels, family = "poisson"))
```

    ## 
    ## Call:
    ## glm(formula = deaths_preL ~ Average_Pop_density_personkm2 + NO.levels + 
    ##     NO2.levels + O3.levels, family = "poisson", data = preL_dt_int)
    ## 
    ## Deviance Residuals: 
    ##       1        2        3        4        5        6        7  
    ## -3.4193  -0.0597   2.5243  -2.8870  -0.4216   3.3526   0.8474  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                   -1.159e+01  9.091e-01  -12.75   <2e-16 ***
    ## Average_Pop_density_personkm2 -9.771e-04  4.570e-05  -21.38   <2e-16 ***
    ## NO.levels                     -4.916e-04  2.301e-05  -21.36   <2e-16 ***
    ## NO2.levels                     8.457e-04  3.571e-05   23.68   <2e-16 ***
    ## O3.levels                      1.247e-04  9.462e-06   13.18   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 1833.017  on 6  degrees of freedom
    ## Residual deviance:   38.539  on 2  degrees of freedom
    ## AIC: 109.09
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
#pchisq(2 * (logLik(full_lm.nb) - logLik(full_lm.p)), df = 1, lower.tail = FALSE)
#(est <- cbind(Estimate = coef(full_lm.nb), confint(full_lm.nb)))
```

#### Analyse the Cases data:

``` r
summary(cases_lm <- lm(data = preL_dt, cases_preL ~ Average_Pop_density_personkm2 + NO.levels + NO2.levels + O3.levels))
```

    ## 
    ## Call:
    ## lm(formula = cases_preL ~ Average_Pop_density_personkm2 + NO.levels + 
    ##     NO2.levels + O3.levels, data = preL_dt)
    ## 
    ## Residuals:
    ##        1        2        3        4        5        6        7 
    ## -1610.04   -62.57   495.92  -375.69  -200.13  2067.35  -314.85 
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                   -5.887e+04  4.405e+04  -1.337    0.313
    ## Average_Pop_density_personkm2 -3.805e+00  2.313e+00  -1.645    0.242
    ## NO.levels                     -2.778e+00  1.022e+00  -2.717    0.113
    ## NO2.levels                     4.118e+00  1.620e+00   2.541    0.126
    ## O3.levels                      2.380e-01  4.921e-01   0.484    0.676
    ## 
    ## Residual standard error: 1923 on 2 degrees of freedom
    ## Multiple R-squared:  0.9365, Adjusted R-squared:  0.8095 
    ## F-statistic: 7.373 on 4 and 2 DF,  p-value: 0.123

``` r
summary(cases_lm.nb <- glm.nb(data = preL_dt_int, cases_preL ~ Average_Pop_density_personkm2 + NO.levels + NO2.levels + O3.levels))
```

    ## 
    ## Call:
    ## glm.nb(formula = cases_preL ~ Average_Pop_density_personkm2 + 
    ##     NO.levels + NO2.levels + O3.levels, data = preL_dt_int, init.theta = 33.53196858, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##        1         2         3         4         5         6         7  
    ## -1.62262  -0.06666  -0.12959   0.24412  -0.18307   1.97477  -0.62489  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                   -4.451e+00  3.968e+00  -1.122 0.261994    
    ## Average_Pop_density_personkm2 -7.942e-04  2.084e-04  -3.811 0.000138 ***
    ## NO.levels                     -4.916e-04  9.215e-05  -5.334 9.60e-08 ***
    ## NO2.levels                     7.404e-04  1.461e-04   5.069 4.01e-07 ***
    ## O3.levels                      6.958e-05  4.430e-05   1.570 0.116300    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(33.532) family taken to be 1)
    ## 
    ##     Null deviance: 55.8888  on 6  degrees of freedom
    ## Residual deviance:  7.0375  on 2  degrees of freedom
    ## AIC: 132.87
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  33.5 
    ##           Std. Err.:  17.9 
    ## 
    ##  2 x log-likelihood:  -120.868

``` r
summary(cases_lm.p <- glm(data = preL_dt_int, cases_preL ~ Average_Pop_density_personkm2 + NO.levels + NO2.levels + O3.levels, family = "poisson"))
```

    ## 
    ## Call:
    ## glm(formula = cases_preL ~ Average_Pop_density_personkm2 + NO.levels + 
    ##     NO2.levels + O3.levels, family = "poisson", data = preL_dt_int)
    ## 
    ## Deviance Residuals: 
    ##        1         2         3         4         5         6         7  
    ## -20.3835   -0.5511   -0.2285    1.7246   -2.0609   28.5914  -10.4853  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                   -2.934e+00  3.006e-01  -9.761   <2e-16 ***
    ## Average_Pop_density_personkm2 -7.301e-04  1.543e-05 -47.324   <2e-16 ***
    ## NO.levels                     -4.524e-04  7.507e-06 -60.273   <2e-16 ***
    ## NO2.levels                     6.787e-04  1.186e-05  57.231   <2e-16 ***
    ## O3.levels                      5.611e-05  3.104e-06  18.080   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 13239.2  on 6  degrees of freedom
    ## Residual deviance:  1350.5  on 2  degrees of freedom
    ## AIC: 1436.1
    ## 
    ## Number of Fisher Scoring iterations: 4

##### Supplementary Table 2. Effect of air pollutants on COVID-19 cases in England at the regional level.

Cases data: table

``` r
stargazer(cases_lm.p, cases_lm.nb,type="html",
          dep.var.labels="Number of cases until 8 April 2020",
          ci=TRUE, ci.level=0.95, single.row=TRUE)
```

<table style="text-align:center">

<tr>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="2">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="2" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="2">

Number of cases until 8 April 2020

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

<em>Poisson</em>

</td>

<td>

<em>negative</em>

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

<em></em>

</td>

<td>

<em>binomial</em>

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

(1)

</td>

<td>

(2)

</td>

</tr>

<tr>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Average\_Pop\_density\_personkm2

</td>

<td>

\-0.001<sup>\*\*\*</sup> (-0.001, -0.001)

</td>

<td>

\-0.001<sup>\*\*\*</sup> (-0.001, -0.0004)

</td>

</tr>

<tr>

<td style="text-align:left">

NO.levels

</td>

<td>

\-0.0005<sup>\*\*\*</sup> (-0.0005, -0.0004)

</td>

<td>

\-0.0005<sup>\*\*\*</sup> (-0.001, -0.0003)

</td>

</tr>

<tr>

<td style="text-align:left">

NO2.levels

</td>

<td>

0.001<sup>\*\*\*</sup> (0.001, 0.001)

</td>

<td>

0.001<sup>\*\*\*</sup> (0.0005, 0.001)

</td>

</tr>

<tr>

<td style="text-align:left">

O3.levels

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.0001, 0.0001)

</td>

<td>

0.0001 (-0.00002, 0.0002)

</td>

</tr>

<tr>

<td style="text-align:left">

Constant

</td>

<td>

\-2.934<sup>\*\*\*</sup> (-3.523, -2.345)

</td>

<td>

\-4.451 (-12.229, 3.327)

</td>

</tr>

<tr>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Observations

</td>

<td>

7

</td>

<td>

7

</td>

</tr>

<tr>

<td style="text-align:left">

Log Likelihood

</td>

<td>

\-713.041

</td>

<td>

\-61.434

</td>

</tr>

<tr>

<td style="text-align:left">

theta

</td>

<td>

</td>

<td>

33.532<sup>\*</sup> (17.929)

</td>

</tr>

<tr>

<td style="text-align:left">

Akaike Inf. Crit.

</td>

<td>

1,436.081

</td>

<td>

132.868

</td>

</tr>

<tr>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="2" style="text-align:right">

<sup>*</sup>p\<0.1; <sup>**</sup>p\<0.05; <sup>***</sup>p\<0.01

</td>

</tr>

</table>

<br> Each column of the table corresponds to a different type of
regression model as indicated at the top. The raw estimate values of
each model are listed with their 95% confidence intervals in
parentheses. The p-values are indicated using the number of asterisks
beside the estimates. OLS, ordinary least square;
Average\_Pop\_densitykm2, average population density per square
kilometer; NO.levels, nitrogen oxide levels; NO2.levels, nitrogen
dioxide levels; O3.levels, ozone levels; Akaike Inf. Crit., Akaike’s
Information Criteria; Residual Std. Error, Residual standard error.

##### Supplementary Table 3. Effect of air pollutants on COVID-19 deaths in England at the regional level.

``` r
stargazer(full_lm.p, full_lm.nb,type="html",
          dep.var.labels="Number of deaths until 8 April 2020",
          ci=TRUE, ci.level=0.95, single.row=TRUE)
```

<table style="text-align:center">

<tr>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="2">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="2" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="2">

Number of deaths until 8 April 2020

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

<em>Poisson</em>

</td>

<td>

<em>negative</em>

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

<em></em>

</td>

<td>

<em>binomial</em>

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

(1)

</td>

<td>

(2)

</td>

</tr>

<tr>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Average\_Pop\_density\_personkm2

</td>

<td>

\-0.001<sup>\*\*\*</sup> (-0.001, -0.001)

</td>

<td>

\-0.001<sup>\*\*\*</sup> (-0.001, -0.001)

</td>

</tr>

<tr>

<td style="text-align:left">

NO.levels

</td>

<td>

\-0.0005<sup>\*\*\*</sup> (-0.001, -0.0004)

</td>

<td>

\-0.0005<sup>\*\*\*</sup> (-0.001, -0.0004)

</td>

</tr>

<tr>

<td style="text-align:left">

NO2.levels

</td>

<td>

0.001<sup>\*\*\*</sup> (0.001, 0.001)

</td>

<td>

0.001<sup>\*\*\*</sup> (0.001, 0.001)

</td>

</tr>

<tr>

<td style="text-align:left">

O3.levels

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.0001, 0.0001)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.0001, 0.0002)

</td>

</tr>

<tr>

<td style="text-align:left">

Constant

</td>

<td>

\-11.590<sup>\*\*\*</sup> (-13.372, -9.809)

</td>

<td>

\-10.900<sup>\*\*\*</sup> (-14.604, -7.195)

</td>

</tr>

<tr>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Observations

</td>

<td>

7

</td>

<td>

7

</td>

</tr>

<tr>

<td style="text-align:left">

Log Likelihood

</td>

<td>

\-49.547

</td>

<td>

\-40.873

</td>

</tr>

<tr>

<td style="text-align:left">

theta

</td>

<td>

</td>

<td>

196.421 (128.791)

</td>

</tr>

<tr>

<td style="text-align:left">

Akaike Inf. Crit.

</td>

<td>

109.094

</td>

<td>

91.745

</td>

</tr>

<tr>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="2" style="text-align:right">

<sup>*</sup>p\<0.1; <sup>**</sup>p\<0.05; <sup>***</sup>p\<0.01

</td>

</tr>

</table>

<br>Each column of the table corresponds to a different type of
regression model as indicated at the top. The raw estimate values of
each model are listed with their 95% confidence intervals in
parentheses. The p-values are indicated using the number of asterisks
beside the estimates. OLS, ordinary least square;
Average\_Pop\_densitykm2, average population density per square
kilometer; NO.levels, nitrogen oxide levels; NO2.levels, nitrogen
dioxide levels; O3.levels, ozone levels; Akaike Inf. Crit., Akaike’s
Information Criteria; Residual Std. Error, Residual standard error.

## Aim 2: subregional level analysis

Aim: Analyse the effect of air pollutants on COVID cases and deaths at
the regional level

### Deaths data

#### Data curation

``` r
pop_dens = read.csv("data_v4/2018_official_popDensity.csv")[c("Code","X2018.people.per.sq..km")]
earnings =read.csv("data_v4/ann_earning_2018_perLA.csv")[c("Code","Mean_ann_earnings")]
age = read.csv("data_v4/processed_median_age_of_population_perLA.csv")[c("Code","median_age_2018","Name")]
covid_deaths = read.csv("data_v4/covid_deaths_until10April_byAreaCode.csv")

covid_deaths$total_deaths = covid_deaths$Home + covid_deaths$Hospital + covid_deaths$Care.home + covid_deaths$Hospice + covid_deaths$Other.communal.establishment + covid_deaths$Elsewhere

covid_deaths = covid_deaths[c("Area.code","total_deaths")]
colnames(covid_deaths) <- c("Code","deaths")

nrow(pop_dens)
```

    ## [1] 435

``` r
nrow(earnings)
```

    ## [1] 422

``` r
nrow(age)
```

    ## [1] 433

``` r
nrow(covid_deaths)
```

    ## [1] 346

``` r
#merge data
merged_covid_dt_LA = merge(covid_deaths, pop_dens, by = "Code")
merged_covid_dt_LA = merge(merged_covid_dt_LA, earnings, by = "Code")
merged_covid_dt_LA = merge(merged_covid_dt_LA, age, by = "Code")
nrow(merged_covid_dt_LA)
```

    ## [1] 334

**Now, we will use a python script (in the form of a Jupyter Notebook)
to match air pollution data to local authorities.** <br> \#\#\#\#
Analysis <br> \#\#\#\#\# Data visualisation for model selection <br>
Load data and set correct formats

``` r
# load data 
colnames(merged_covid_dt_LA)
```

    ## [1] "Code"                    "deaths"                 
    ## [3] "X2018.people.per.sq..km" "Mean_ann_earnings"      
    ## [5] "median_age_2018"         "Name"

``` r
merged_covid_dt_LA$X2018.people.per.sq..km = as.numeric(gsub(",","",merged_covid_dt_LA$X2018.people.per.sq..km))
merged_covid_dt_LA$Mean_ann_earnings = as.numeric(gsub(",","",merged_covid_dt_LA$Mean_ann_earnings))

write.csv(merged_covid_dt_LA, "data_output_v4/merged_covid_cov_dt_LA.csv")
```

re-read here

``` r
covid_air_dt = read.csv("data_output_v4/merged_covidAir_cov_dt_LA.csv", na.strings = "x")
covid_air_dt$X2018.people.per.sq..km = as.numeric(gsub(",","",covid_air_dt$X2018.people.per.sq..km))
covid_air_dt$Mean_ann_earnings = as.numeric(gsub(",","",covid_air_dt$Mean_ann_earnings))
covid_air_dt$X2018.people.per.sq..km = as.numeric(covid_air_dt$X2018.people.per.sq..km)
covid_air_dt$Mean_ann_earnings = as.numeric(covid_air_dt$Mean_ann_earnings)
```

Make sure that only England data is included

``` r
covid_air_dt = covid_air_dt[startsWith(as.character(covid_air_dt$Code), 'E'),]
nrow(covid_air_dt)
```

    ## [1] 312

Visualise all numeric data

``` r
covid_air_dt_vis = subset(covid_air_dt, select=c(deaths, X2018.people.per.sq..km, Mean_ann_earnings,median_age_2018,pm25_val,no2_val,o3_val,
                                                 pm10_val,so2_val, nox_val))
covid_air_dt_vis = as.data.frame(sapply(covid_air_dt_vis, function(x) as.numeric(x) ) )
#normalise per column & reshape the data for plotting
covid_air_dt_vis = as.data.frame(sapply(covid_air_dt_vis, scale))

covid_air_dt_vis$Code = covid_air_dt$Code
covid_air_dt_vis_long = melt(covid_air_dt_vis, in.vars = "Code")
```

    ## Using Code as id variables

``` r
ggplot(covid_air_dt_vis_long, aes (value)) +
  geom_histogram() +
  geom_density() +
  facet_wrap(~variable) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

##### Models

Negative binomial regression model

``` r
summary(pm25_deaths.nb <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                   Mean_ann_earnings + median_age_2018 + pm25_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ X2018.people.per.sq..km + Mean_ann_earnings + 
    ##     median_age_2018 + pm25_val, data = covid_air_dt, init.theta = 2.099197434, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.4419  -0.8381  -0.3209   0.2889   3.2137  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              6.988e+00  6.351e-01  11.002  < 2e-16 ***
    ## X2018.people.per.sq..km  1.040e-04  2.232e-05   4.661 3.15e-06 ***
    ## Mean_ann_earnings        4.670e-06  6.535e-06   0.715    0.475    
    ## median_age_2018         -7.862e-02  1.188e-02  -6.618 3.64e-11 ***
    ## pm25_val                -3.892e-02  3.146e-02  -1.237    0.216    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.0992) family taken to be 1)
    ## 
    ##     Null deviance: 528.24  on 305  degrees of freedom
    ## Residual deviance: 328.09  on 301  degrees of freedom
    ##   (6 observations deleted due to missingness)
    ## AIC: 2801.1
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.099 
    ##           Std. Err.:  0.171 
    ## 
    ##  2 x log-likelihood:  -2789.130

``` r
car::vif(pm25_deaths.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.080989                1.441089                2.135625 
    ##                pm25_val 
    ##                1.950479

``` r
summary(pm10_deaths.nb <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                   Mean_ann_earnings + median_age_2018 + pm10_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ X2018.people.per.sq..km + Mean_ann_earnings + 
    ##     median_age_2018 + pm10_val, data = covid_air_dt, init.theta = 2.108915137, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.4457  -0.8477  -0.3321   0.2782   3.2580  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              7.103e+00  5.970e-01  11.898  < 2e-16 ***
    ## X2018.people.per.sq..km  1.073e-04  2.237e-05   4.798 1.61e-06 ***
    ## Mean_ann_earnings        5.482e-06  6.392e-06   0.858   0.3911    
    ## median_age_2018         -7.866e-02  1.133e-02  -6.943 3.85e-12 ***
    ## pm10_val                -3.566e-02  2.001e-02  -1.782   0.0748 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.1089) family taken to be 1)
    ## 
    ##     Null deviance: 530.51  on 305  degrees of freedom
    ## Residual deviance: 327.92  on 301  degrees of freedom
    ##   (6 observations deleted due to missingness)
    ## AIC: 2799.6
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.109 
    ##           Std. Err.:  0.171 
    ## 
    ##  2 x log-likelihood:  -2787.577

``` r
car::vif(pm10_deaths.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.099725                1.385040                1.950317 
    ##                pm10_val 
    ##                1.667899

``` r
summary(nox_deaths.nb <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                  Mean_ann_earnings + median_age_2018 + nox_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ X2018.people.per.sq..km + Mean_ann_earnings + 
    ##     median_age_2018 + nox_val, data = covid_air_dt, init.theta = 2.142638423, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.3827  -0.8640  -0.3110   0.3171   3.6497  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              5.541e+00  5.422e-01  10.220  < 2e-16 ***
    ## X2018.people.per.sq..km  7.321e-05  2.437e-05   3.004  0.00266 ** 
    ## Mean_ann_earnings       -1.162e-06  6.088e-06  -0.191  0.84861    
    ## median_age_2018         -5.455e-02  1.174e-02  -4.645  3.4e-06 ***
    ## nox_val                  1.602e-02  5.517e-03   2.904  0.00369 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.1426) family taken to be 1)
    ## 
    ##     Null deviance: 538.40  on 305  degrees of freedom
    ## Residual deviance: 327.98  on 301  degrees of freedom
    ##   (6 observations deleted due to missingness)
    ## AIC: 2794.9
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.143 
    ##           Std. Err.:  0.175 
    ## 
    ##  2 x log-likelihood:  -2782.865

``` r
car::vif(nox_deaths.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.529754                1.277034                2.128073 
    ##                 nox_val 
    ##                2.650417

``` r
summary(no2_deaths.nb <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                  Mean_ann_earnings + median_age_2018 + no2_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ X2018.people.per.sq..km + Mean_ann_earnings + 
    ##     median_age_2018 + no2_val, data = covid_air_dt, init.theta = 2.156887438, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.3845  -0.8756  -0.3018   0.3214   3.7677  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              5.268e+00  5.802e-01   9.080  < 2e-16 ***
    ## X2018.people.per.sq..km  7.548e-05  2.375e-05   3.178  0.00148 ** 
    ## Mean_ann_earnings       -1.397e-06  6.037e-06  -0.231  0.81702    
    ## median_age_2018         -5.025e-02  1.221e-02  -4.116 3.86e-05 ***
    ## no2_val                  2.961e-02  9.706e-03   3.051  0.00228 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.1569) family taken to be 1)
    ## 
    ##     Null deviance: 541.72  on 305  degrees of freedom
    ## Residual deviance: 327.96  on 301  degrees of freedom
    ##   (6 observations deleted due to missingness)
    ## AIC: 2792.9
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.157 
    ##           Std. Err.:  0.176 
    ## 
    ##  2 x log-likelihood:  -2780.853

``` r
car::vif(no2_deaths.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.417227                1.263175                2.314129 
    ##                 no2_val 
    ##                2.736928

``` r
summary(o3_deaths.nb <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                 Mean_ann_earnings + median_age_2018 + o3_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ X2018.people.per.sq..km + Mean_ann_earnings + 
    ##     median_age_2018 + o3_val, data = covid_air_dt, init.theta = 2.14216491, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.3916  -0.8383  -0.2967   0.3281   3.1013  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              6.323e+00  4.890e-01  12.931  < 2e-16 ***
    ## X2018.people.per.sq..km  9.218e-05  2.217e-05   4.158 3.21e-05 ***
    ## Mean_ann_earnings        8.137e-06  6.092e-06   1.336   0.1816    
    ## median_age_2018         -6.655e-02  1.089e-02  -6.110 9.99e-10 ***
    ## o3_val                  -3.690e-02  1.252e-02  -2.948   0.0032 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.1422) family taken to be 1)
    ## 
    ##     Null deviance: 538.29  on 305  degrees of freedom
    ## Residual deviance: 327.69  on 301  degrees of freedom
    ##   (6 observations deleted due to missingness)
    ## AIC: 2794.6
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.142 
    ##           Std. Err.:  0.175 
    ## 
    ##  2 x log-likelihood:  -2782.640

``` r
car::vif(o3_deaths.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.094533                1.278499                1.827214 
    ##                  o3_val 
    ##                1.148228

``` r
summary(so2_deaths.nb <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                  Mean_ann_earnings + median_age_2018 + so2_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ X2018.people.per.sq..km + Mean_ann_earnings + 
    ##     median_age_2018 + so2_val, data = covid_air_dt, init.theta = 2.12911154, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.3778  -0.8502  -0.3190   0.2988   3.5655  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              5.498e+00  5.956e-01   9.230  < 2e-16 ***
    ## X2018.people.per.sq..km  9.388e-05  2.222e-05   4.225 2.39e-05 ***
    ## Mean_ann_earnings        5.355e-06  6.001e-06   0.892   0.3722    
    ## median_age_2018         -5.850e-02  1.176e-02  -4.973 6.60e-07 ***
    ## so2_val                  2.016e-01  8.341e-02   2.417   0.0157 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.1291) family taken to be 1)
    ## 
    ##     Null deviance: 535.24  on 305  degrees of freedom
    ## Residual deviance: 328.07  on 301  degrees of freedom
    ##   (6 observations deleted due to missingness)
    ## AIC: 2796.9
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.129 
    ##           Std. Err.:  0.174 
    ## 
    ##  2 x log-likelihood:  -2784.854

``` r
car::vif(so2_deaths.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.090598                1.233099                2.120397 
    ##                 so2_val 
    ##                1.497102

We note that annual earnings in our models is not significant. We
therefore proceed to remove this variable.

``` r
# glm -earnings, since not significant 
summary(pm25_deaths.nb_red <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                    median_age_2018 + pm25_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ X2018.people.per.sq..km + median_age_2018 + 
    ##     pm25_val, data = covid_air_dt, init.theta = 2.126869187, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.4717  -0.8441  -0.3194   0.2899   3.1964  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              6.9995914  0.6180285  11.326  < 2e-16 ***
    ## X2018.people.per.sq..km  0.0001027  0.0000205   5.009 5.46e-07 ***
    ## median_age_2018         -0.0772557  0.0111724  -6.915 4.68e-12 ***
    ## pm25_val                -0.0306711  0.0282245  -1.087    0.277    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.1269) family taken to be 1)
    ## 
    ##     Null deviance: 538.14  on 311  degrees of freedom
    ## Residual deviance: 334.35  on 308  degrees of freedom
    ## AIC: 2850.5
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.127 
    ##           Std. Err.:  0.171 
    ## 
    ##  2 x log-likelihood:  -2840.463

``` r
car::vif(pm25_deaths.nb_red)
```

    ## X2018.people.per.sq..km         median_age_2018                pm25_val 
    ##                1.894731                1.945670                1.639378

``` r
summary(pm10_deaths.nb_red <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                    median_age_2018 + pm10_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ X2018.people.per.sq..km + median_age_2018 + 
    ##     pm10_val, data = covid_air_dt, init.theta = 2.134252985, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.4763  -0.8503  -0.2993   0.3085   3.2222  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              7.109e+00  5.820e-01  12.214  < 2e-16 ***
    ## X2018.people.per.sq..km  1.063e-04  2.058e-05   5.165 2.40e-07 ***
    ## median_age_2018         -7.717e-02  1.079e-02  -7.153 8.47e-13 ***
    ## pm10_val                -2.849e-02  1.829e-02  -1.558    0.119    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.1343) family taken to be 1)
    ## 
    ##     Null deviance: 539.87  on 311  degrees of freedom
    ## Residual deviance: 334.24  on 308  degrees of freedom
    ## AIC: 2849.3
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.134 
    ##           Std. Err.:  0.172 
    ## 
    ##  2 x log-likelihood:  -2839.284

``` r
car::vif(pm10_deaths.nb_red)
```

    ## X2018.people.per.sq..km         median_age_2018                pm10_val 
    ##                1.915727                1.819813                1.458236

``` r
summary(nox_deaths.nb_red <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                   median_age_2018 + nox_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ X2018.people.per.sq..km + median_age_2018 + 
    ##     nox_val, data = covid_air_dt, init.theta = 2.168572732, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.3994  -0.8610  -0.3049   0.3305   3.6400  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              5.614e+00  5.258e-01  10.676  < 2e-16 ***
    ## X2018.people.per.sq..km  6.762e-05  2.311e-05   2.927  0.00343 ** 
    ## median_age_2018         -5.650e-02  1.120e-02  -5.047 4.49e-07 ***
    ## nox_val                  1.519e-02  5.285e-03   2.874  0.00405 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.1686) family taken to be 1)
    ## 
    ##     Null deviance: 547.93  on 311  degrees of freedom
    ## Residual deviance: 334.17  on 308  degrees of freedom
    ## AIC: 2844.3
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.169 
    ##           Std. Err.:  0.175 
    ## 
    ##  2 x log-likelihood:  -2834.334

``` r
car::vif(nox_deaths.nb_red)
```

    ## X2018.people.per.sq..km         median_age_2018                 nox_val 
    ##                2.452140                1.988632                2.557030

``` r
summary(no2_deaths.nb_red <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                   median_age_2018 + no2_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ X2018.people.per.sq..km + median_age_2018 + 
    ##     no2_val, data = covid_air_dt, init.theta = 2.182039301, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.3997  -0.8565  -0.2959   0.3142   3.7547  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              5.344e+00  5.644e-01   9.468  < 2e-16 ***
    ## X2018.people.per.sq..km  6.950e-05  2.234e-05   3.112  0.00186 ** 
    ## median_age_2018         -5.235e-02  1.163e-02  -4.502 6.72e-06 ***
    ## no2_val                  2.817e-02  9.339e-03   3.017  0.00256 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.182) family taken to be 1)
    ## 
    ##     Null deviance: 551.09  on 311  degrees of freedom
    ## Residual deviance: 334.14  on 308  degrees of freedom
    ## AIC: 2842.4
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.182 
    ##           Std. Err.:  0.177 
    ## 
    ##  2 x log-likelihood:  -2832.416

``` r
car::vif(no2_deaths.nb_red)
```

    ## X2018.people.per.sq..km         median_age_2018                 no2_val 
    ##                2.304673                2.157462                2.660320

``` r
summary(o3_deaths.nb_red <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                  median_age_2018 + o3_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ X2018.people.per.sq..km + median_age_2018 + 
    ##     o3_val, data = covid_air_dt, init.theta = 2.163965779, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.4405  -0.8565  -0.3085   0.3545   3.0899  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              6.571e+00  4.641e-01  14.158  < 2e-16 ***
    ## X2018.people.per.sq..km  9.432e-05  1.946e-05   4.847 1.26e-06 ***
    ## median_age_2018         -6.772e-02  1.055e-02  -6.417 1.39e-10 ***
    ## o3_val                  -3.107e-02  1.195e-02  -2.600  0.00933 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.164) family taken to be 1)
    ## 
    ##     Null deviance: 546.85  on 311  degrees of freedom
    ## Residual deviance: 334.05  on 308  degrees of freedom
    ## AIC: 2844.9
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.164 
    ##           Std. Err.:  0.175 
    ## 
    ##  2 x log-likelihood:  -2834.868

``` r
car::vif(o3_deaths.nb_red)
```

    ## X2018.people.per.sq..km         median_age_2018                  o3_val 
    ##                1.736122                1.762935                1.072928

``` r
summary(so2_deaths.nb_red <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                   median_age_2018 + so2_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ X2018.people.per.sq..km + median_age_2018 + 
    ##     so2_val, data = covid_air_dt, init.theta = 2.15583651, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.4164  -0.8464  -0.3192   0.3201   3.5167  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              5.728e+00  5.583e-01  10.260  < 2e-16 ***
    ## X2018.people.per.sq..km  9.475e-05  1.952e-05   4.855 1.21e-06 ***
    ## median_age_2018         -5.959e-02  1.148e-02  -5.189 2.11e-07 ***
    ## so2_val                  1.850e-01  8.147e-02   2.271   0.0231 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.1558) family taken to be 1)
    ## 
    ##     Null deviance: 544.94  on 311  degrees of freedom
    ## Residual deviance: 334.32  on 308  degrees of freedom
    ## AIC: 2846.3
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.156 
    ##           Std. Err.:  0.174 
    ## 
    ##  2 x log-likelihood:  -2836.287

``` r
car::vif(so2_deaths.nb_red)
```

    ## X2018.people.per.sq..km         median_age_2018                 so2_val 
    ##                1.739256                2.080587                1.458499

<br>

##### Calculate odds ratios

``` r
pm25_deaths.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(pm25_deaths.nb), confint(pm25_deaths.nb))), p_value = summary(pm25_deaths.nb)$coefficients[,4]))
pm10_deaths.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(pm10_deaths.nb), confint(pm10_deaths.nb))), p_value = summary(pm10_deaths.nb)$coefficients[,4]))
nox_deaths.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(nox_deaths.nb), confint(nox_deaths.nb))), p_value = summary(nox_deaths.nb)$coefficients[,4]))
no2_deaths.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(no2_deaths.nb), confint(no2_deaths.nb))), p_value = summary(no2_deaths.nb)$coefficients[,4]))
o3_deaths.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(o3_deaths.nb), confint(o3_deaths.nb))), p_value = summary(o3_deaths.nb)$coefficients[,4]))
so2_deaths.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(so2_deaths.nb), confint(so2_deaths.nb))), p_value = summary(so2_deaths.nb)$coefficients[,4]))
#make a df just with the pollutants 

LA_covid_onlyPoll = data.frame(rbind(pm25_deaths.nb_mrr[nrow(pm25_deaths.nb_mrr),],
                                      pm10_deaths.nb_mrr[nrow(pm10_deaths.nb_mrr),],
                                      nox_deaths.nb_mrr[nrow(nox_deaths.nb_mrr),],
                                      no2_deaths.nb_mrr[nrow(no2_deaths.nb_mrr),],
                                      o3_deaths.nb_mrr[nrow(o3_deaths.nb_mrr),],
                                      so2_deaths.nb_mrr[nrow(so2_deaths.nb_mrr),]))
```

##### Plot odds ratios

sort data

``` r
LA_covid_onlyPoll$names = row.names(LA_covid_onlyPoll)
LA_covid_onlyPoll$significance = "p-value > 0.05"
LA_covid_onlyPoll$significance[LA_covid_onlyPoll$p_value < 0.05] <- "p-value < 0.05"
```

plot

``` r
ggplot(LA_covid_onlyPoll, aes(x=reorder(names, OR), y=OR, color=significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_bw() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("Mortality rate ratios") + 
  xlab("Pollutants")+
  ylim(0.75, 1.6)
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
ggsave('fig_out_v4/LA_deaths_MRR.pdf')
```

    ## Saving 7 x 5 in image

### Cases in local authorities

#### Data curation

``` r
cases_LA_raw = read.csv("data_v4/coronavirus-cases_latest-18_5_2020.csv")
cases_LA_dt = subset(cases_LA_raw, Specimen.date == "2020-04-10", select = c(Area.code, Cumulative.lab.confirmed.cases))
cases_LA_dt.agg <-aggregate(cases_LA_dt, by=list(cases_LA_dt$Area.code), 
  FUN=mean, na.rm=TRUE)
nrow(cases_LA_dt.agg)
```

    ## [1] 342

Merge cases data to main data

``` r
deaths_LA_dt = read.csv("data_output_v4/merged_covidAir_cov_dt_LA.csv", na.strings = "x")
deaths_LA_dt$X2018.people.per.sq..km = gsub(",", "", deaths_LA_dt$X2018.people.per.sq..km)
deaths_LA_dt$X2018.people.per.sq..km = as.numeric(deaths_LA_dt$X2018.people.per.sq..km)
deaths_LA_dt$Mean_ann_earnings = gsub(",", "", deaths_LA_dt$Mean_ann_earnings)
deaths_LA_dt$Mean_ann_earnings = as.numeric(deaths_LA_dt$Mean_ann_earnings)

cases_deaths_LA_dt = merge(deaths_LA_dt, cases_LA_dt.agg, by.x = "Code",by.y = "Group.1")
nrow(cases_deaths_LA_dt)
```

    ## [1] 301

There are only 301 observations matched.

#### Analysis of the cases data

Visualisation

``` r
ggplot(data = cases_deaths_LA_dt, aes(x = Cumulative.lab.confirmed.cases))+
  geom_histogram(stat="count")
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
cases_deaths_LA_dt$X2018.people.per.sq..km = as.numeric(cases_deaths_LA_dt$X2018.people.per.sq..km)
cases_deaths_LA_dt$Mean_ann_earnings = as.numeric(cases_deaths_LA_dt$Mean_ann_earnings)
```

Model

``` r
summary(pm25_cases.nb <- glm.nb(data = cases_deaths_LA_dt, Cumulative.lab.confirmed.cases ~ X2018.people.per.sq..km + 
                                   Mean_ann_earnings + median_age_2018 + pm25_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = Cumulative.lab.confirmed.cases ~ X2018.people.per.sq..km + 
    ##     Mean_ann_earnings + median_age_2018 + pm25_val, data = cases_deaths_LA_dt, 
    ##     init.theta = 2.724975873, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.8555  -0.8682  -0.3414   0.2772   2.8464  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              9.593e+00  5.567e-01  17.232  < 2e-16 ***
    ## X2018.people.per.sq..km  9.626e-05  1.953e-05   4.930 8.21e-07 ***
    ## Mean_ann_earnings        6.492e-06  5.788e-06   1.122 0.261996    
    ## median_age_2018         -9.122e-02  1.037e-02  -8.800  < 2e-16 ***
    ## pm25_val                -9.629e-02  2.770e-02  -3.476 0.000509 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.725) family taken to be 1)
    ## 
    ##     Null deviance: 562.55  on 294  degrees of freedom
    ## Residual deviance: 311.82  on 290  degrees of freedom
    ##   (6 observations deleted due to missingness)
    ## AIC: 3583.7
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.725 
    ##           Std. Err.:  0.216 
    ## 
    ##  2 x log-likelihood:  -3571.745

``` r
car::vif(pm25_cases.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.085624                1.427671                2.125130 
    ##                pm25_val 
    ##                1.933232

``` r
summary(pm10_cases.nb <- glm.nb(data = cases_deaths_LA_dt, Cumulative.lab.confirmed.cases ~ X2018.people.per.sq..km + 
                                   Mean_ann_earnings + median_age_2018 + pm10_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = Cumulative.lab.confirmed.cases ~ X2018.people.per.sq..km + 
    ##     Mean_ann_earnings + median_age_2018 + pm10_val, data = cases_deaths_LA_dt, 
    ##     init.theta = 2.745559885, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.8501  -0.8823  -0.3287   0.2934   2.9348  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              9.513e+00  5.216e-01  18.238  < 2e-16 ***
    ## X2018.people.per.sq..km  1.002e-04  1.956e-05   5.122 3.03e-07 ***
    ## Mean_ann_earnings        5.806e-06  5.650e-06   1.028 0.304083    
    ## median_age_2018         -8.743e-02  9.861e-03  -8.866  < 2e-16 ***
    ## pm10_val                -6.748e-02  1.762e-02  -3.831 0.000128 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.7456) family taken to be 1)
    ## 
    ##     Null deviance: 566.73  on 294  degrees of freedom
    ## Residual deviance: 311.67  on 290  degrees of freedom
    ##   (6 observations deleted due to missingness)
    ## AIC: 3581.3
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.746 
    ##           Std. Err.:  0.218 
    ## 
    ##  2 x log-likelihood:  -3569.297

``` r
car::vif(pm10_cases.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.109349                1.370273                1.937009 
    ##                pm10_val 
    ##                1.650101

``` r
summary(nox_cases.nb <- glm.nb(data = cases_deaths_LA_dt, Cumulative.lab.confirmed.cases ~ X2018.people.per.sq..km + 
                                  Mean_ann_earnings + median_age_2018 + nox_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = Cumulative.lab.confirmed.cases ~ X2018.people.per.sq..km + 
    ##     Mean_ann_earnings + median_age_2018 + nox_val, data = cases_deaths_LA_dt, 
    ##     init.theta = 2.721556683, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.7665  -0.8356  -0.3858   0.3034   3.4862  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              7.409e+00  4.954e-01  14.956  < 2e-16 ***
    ## X2018.people.per.sq..km  5.246e-05  2.220e-05   2.363  0.01813 *  
    ## Mean_ann_earnings       -4.812e-06  5.385e-06  -0.894  0.37145    
    ## median_age_2018         -5.803e-02  1.051e-02  -5.522 3.35e-08 ***
    ## nox_val                  1.740e-02  5.344e-03   3.255  0.00113 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.7216) family taken to be 1)
    ## 
    ##     Null deviance: 561.85  on 294  degrees of freedom
    ## Residual deviance: 312.05  on 290  degrees of freedom
    ##   (6 observations deleted due to missingness)
    ## AIC: 3584.4
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.722 
    ##           Std. Err.:  0.216 
    ## 
    ##  2 x log-likelihood:  -3572.364

``` r
car::vif(nox_cases.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.692080                1.233985                2.182823 
    ##                 nox_val 
    ##                2.905970

``` r
summary(no2_cases.nb <- glm.nb(data = cases_deaths_LA_dt, Cumulative.lab.confirmed.cases ~ X2018.people.per.sq..km + 
                                  Mean_ann_earnings + median_age_2018 + no2_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = Cumulative.lab.confirmed.cases ~ X2018.people.per.sq..km + 
    ##     Mean_ann_earnings + median_age_2018 + no2_val, data = cases_deaths_LA_dt, 
    ##     init.theta = 2.729294735, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.7594  -0.8476  -0.3940   0.3071   3.5818  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              7.190e+00  5.306e-01  13.550  < 2e-16 ***
    ## X2018.people.per.sq..km  5.698e-05  2.142e-05   2.661  0.00780 ** 
    ## Mean_ann_earnings       -4.853e-06  5.365e-06  -0.904  0.36574    
    ## median_age_2018         -5.476e-02  1.094e-02  -5.007 5.52e-07 ***
    ## no2_val                  3.005e-02  9.141e-03   3.287  0.00101 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.7293) family taken to be 1)
    ## 
    ##     Null deviance: 563.42  on 294  degrees of freedom
    ## Residual deviance: 312.03  on 290  degrees of freedom
    ##   (6 observations deleted due to missingness)
    ## AIC: 3583.5
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.729 
    ##           Std. Err.:  0.216 
    ## 
    ##  2 x log-likelihood:  -3571.468

``` r
car::vif(no2_cases.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.512648                1.228516                2.370271 
    ##                 no2_val 
    ##                2.910421

``` r
summary(o3_cases.nb <- glm.nb(data = cases_deaths_LA_dt, Cumulative.lab.confirmed.cases ~ X2018.people.per.sq..km + 
                                 Mean_ann_earnings + median_age_2018 + o3_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = Cumulative.lab.confirmed.cases ~ X2018.people.per.sq..km + 
    ##     Mean_ann_earnings + median_age_2018 + o3_val, data = cases_deaths_LA_dt, 
    ##     init.theta = 2.808347379, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.9822  -0.8344  -0.3153   0.3226   2.8354  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              8.206e+00  4.265e-01  19.241  < 2e-16 ***
    ## X2018.people.per.sq..km  6.994e-05  1.934e-05   3.616 0.000299 ***
    ## Mean_ann_earnings        7.689e-06  5.480e-06   1.403 0.160572    
    ## median_age_2018         -6.928e-02  9.456e-03  -7.327 2.36e-13 ***
    ## o3_val                  -5.259e-02  1.111e-02  -4.733 2.21e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.8083) family taken to be 1)
    ## 
    ##     Null deviance: 579.47  on 294  degrees of freedom
    ## Residual deviance: 311.40  on 290  degrees of freedom
    ##   (6 observations deleted due to missingness)
    ## AIC: 3574.1
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.808 
    ##           Std. Err.:  0.223 
    ## 
    ##  2 x log-likelihood:  -3562.114

``` r
car::vif(o3_cases.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.108263                1.318626                1.820923 
    ##                  o3_val 
    ##                1.177405

``` r
summary(so2_cases.nb <- glm.nb(data = cases_deaths_LA_dt, Cumulative.lab.confirmed.cases ~ X2018.people.per.sq..km + 
                                  Mean_ann_earnings + median_age_2018 + so2_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = Cumulative.lab.confirmed.cases ~ X2018.people.per.sq..km + 
    ##     Mean_ann_earnings + median_age_2018 + so2_val, data = cases_deaths_LA_dt, 
    ##     init.theta = 2.717269061, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.7230  -0.8505  -0.3570   0.3183   3.4351  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              7.239e+00  5.311e-01  13.629  < 2e-16 ***
    ## X2018.people.per.sq..km  7.709e-05  1.964e-05   3.925 8.66e-05 ***
    ## Mean_ann_earnings        2.288e-06  5.451e-06   0.420  0.67462    
    ## median_age_2018         -6.018e-02  1.036e-02  -5.811 6.20e-09 ***
    ## so2_val                  2.370e-01  7.467e-02   3.174  0.00151 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.7173) family taken to be 1)
    ## 
    ##     Null deviance: 560.98  on 294  degrees of freedom
    ## Residual deviance: 312.07  on 290  degrees of freedom
    ##   (6 observations deleted due to missingness)
    ## AIC: 3584.9
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.717 
    ##           Std. Err.:  0.215 
    ## 
    ##  2 x log-likelihood:  -3572.864

``` r
car::vif(so2_cases.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.104075                1.262575                2.115436 
    ##                 so2_val 
    ##                1.525447

#### Generate the table

``` r
stargazer(pm25_cases.nb, pm10_cases.nb, nox_cases.nb, no2_cases.nb,
          o3_cases.nb, so2_cases.nb,type="html",out = "fig_out_v4/LA_covid_allPoll_CASES_nb.html",
          dep.var.labels="Number of COVID-19-related cases",
          single.row=TRUE)
```

<table style="text-align:center">

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="6">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="6" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="6">

Number of COVID-19-related cases

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

(1)

</td>

<td>

(2)

</td>

<td>

(3)

</td>

<td>

(4)

</td>

<td>

(5)

</td>

<td>

(6)

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

X2018.people.per.sq..km

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

</tr>

<tr>

<td style="text-align:left">

Mean\_ann\_earnings

</td>

<td>

0.00001 (0.00001)

</td>

<td>

0.00001 (0.00001)

</td>

<td>

\-0.00000 (0.00001)

</td>

<td>

\-0.00000 (0.00001)

</td>

<td>

0.00001 (0.00001)

</td>

<td>

0.00000 (0.00001)

</td>

</tr>

<tr>

<td style="text-align:left">

median\_age\_2018

</td>

<td>

\-0.091<sup>\*\*\*</sup> (0.010)

</td>

<td>

\-0.087<sup>\*\*\*</sup> (0.010)

</td>

<td>

\-0.058<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.055<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.069<sup>\*\*\*</sup> (0.009)

</td>

<td>

\-0.060<sup>\*\*\*</sup> (0.010)

</td>

</tr>

<tr>

<td style="text-align:left">

pm25\_val

</td>

<td>

\-0.096<sup>\*\*\*</sup> (0.028)

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

pm10\_val

</td>

<td>

</td>

<td>

\-0.067<sup>\*\*\*</sup> (0.018)

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

nox\_val

</td>

<td>

</td>

<td>

</td>

<td>

0.017<sup>\*\*\*</sup> (0.005)

</td>

<td>

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

no2\_val

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

0.030<sup>\*\*\*</sup> (0.009)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

o3\_val

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

\-0.053<sup>\*\*\*</sup> (0.011)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

so2\_val

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

0.237<sup>\*\*\*</sup> (0.075)

</td>

</tr>

<tr>

<td style="text-align:left">

Constant

</td>

<td>

9.593<sup>\*\*\*</sup> (0.557)

</td>

<td>

9.513<sup>\*\*\*</sup> (0.522)

</td>

<td>

7.409<sup>\*\*\*</sup> (0.495)

</td>

<td>

7.190<sup>\*\*\*</sup> (0.531)

</td>

<td>

8.206<sup>\*\*\*</sup> (0.426)

</td>

<td>

7.239<sup>\*\*\*</sup> (0.531)

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Observations

</td>

<td>

295

</td>

<td>

295

</td>

<td>

295

</td>

<td>

295

</td>

<td>

295

</td>

<td>

295

</td>

</tr>

<tr>

<td style="text-align:left">

Log Likelihood

</td>

<td>

\-1,786.872

</td>

<td>

\-1,785.649

</td>

<td>

\-1,787.182

</td>

<td>

\-1,786.734

</td>

<td>

\-1,782.057

</td>

<td>

\-1,787.432

</td>

</tr>

<tr>

<td style="text-align:left">

theta

</td>

<td>

2.725<sup>\*\*\*</sup> (0.216)

</td>

<td>

2.746<sup>\*\*\*</sup> (0.218)

</td>

<td>

2.722<sup>\*\*\*</sup> (0.216)

</td>

<td>

2.729<sup>\*\*\*</sup> (0.216)

</td>

<td>

2.808<sup>\*\*\*</sup> (0.223)

</td>

<td>

2.717<sup>\*\*\*</sup> (0.215)

</td>

</tr>

<tr>

<td style="text-align:left">

Akaike Inf. Crit.

</td>

<td>

3,583.745

</td>

<td>

3,581.297

</td>

<td>

3,584.364

</td>

<td>

3,583.468

</td>

<td>

3,574.114

</td>

<td>

3,584.864

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="6" style="text-align:right">

<sup>*</sup>p\<0.1; <sup>**</sup>p\<0.05; <sup>***</sup>p\<0.01

</td>

</tr>

</table>

##### Calculate odds ratios

``` r
pm25_cases.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(pm25_cases.nb), confint(pm25_cases.nb))), p_value = summary(pm25_cases.nb)$coefficients[,4]))
pm10_cases.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(pm10_cases.nb), confint(pm10_cases.nb))), p_value = summary(pm10_cases.nb)$coefficients[,4]))
nox_cases.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(nox_cases.nb), confint(nox_cases.nb))), p_value = summary(nox_cases.nb)$coefficients[,4]))
no2_cases.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(no2_cases.nb), confint(no2_cases.nb))), p_value = summary(no2_cases.nb)$coefficients[,4]))
o3_cases.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(o3_cases.nb), confint(o3_cases.nb))), p_value = summary(o3_cases.nb)$coefficients[,4]))
so2_cases.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(so2_cases.nb), confint(so2_cases.nb))), p_value = summary(so2_cases.nb)$coefficients[,4]))
#make a df just with the pollutants 

LA_covid_CASES_onlyPoll = data.frame(rbind(pm25_cases.nb_mrr[nrow(pm25_cases.nb_mrr),],
                                      pm10_cases.nb_mrr[nrow(pm10_cases.nb_mrr),],
                                      nox_cases.nb_mrr[nrow(nox_cases.nb_mrr),],
                                      no2_cases.nb_mrr[nrow(no2_cases.nb_mrr),],
                                      o3_cases.nb_mrr[nrow(o3_cases.nb_mrr),],
                                      so2_cases.nb_mrr[nrow(so2_cases.nb_mrr),]))
```

##### Plot odds ratios

sort data

``` r
LA_covid_CASES_onlyPoll$names = row.names(LA_covid_CASES_onlyPoll)
LA_covid_CASES_onlyPoll$significance = "p-value > 0.05"
LA_covid_CASES_onlyPoll$significance[LA_covid_CASES_onlyPoll$p_value < 0.05] <- "p-value < 0.05"
```

plot

``` r
ggplot(LA_covid_CASES_onlyPoll, aes(x=reorder(names, OR), y=OR, color=significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_bw() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("Infectivity rate ratios") + 
  xlab("Pollutants")+
  ylim(0.75, 1.6)
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
ggsave('fig_out_v4/LA_CASES_odds_IRR.pdf')
```

Save the numbers

``` r
library(stargazer)
stargazer(LA_covid_CASES_onlyPoll, summary=FALSE ,type="html",out 
          ="fig_out_v4/LA_covid_infectivity_IRR.html")
```

<table style="text-align:center">

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

OR

</td>

<td>

X2.5..

</td>

<td>

X97.5..

</td>

<td>

p\_value

</td>

<td>

names

</td>

<td>

significance

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

pm25\_val

</td>

<td>

0.908

</td>

<td>

0.860

</td>

<td>

0.958

</td>

<td>

0.001

</td>

<td>

pm25\_val

</td>

<td>

p-value \< 0.05

</td>

</tr>

<tr>

<td style="text-align:left">

pm10\_val

</td>

<td>

0.935

</td>

<td>

0.903

</td>

<td>

0.967

</td>

<td>

0.0001

</td>

<td>

pm10\_val

</td>

<td>

p-value \< 0.05

</td>

</tr>

<tr>

<td style="text-align:left">

nox\_val

</td>

<td>

1.018

</td>

<td>

1.008

</td>

<td>

1.028

</td>

<td>

0.001

</td>

<td>

nox\_val

</td>

<td>

p-value \< 0.05

</td>

</tr>

<tr>

<td style="text-align:left">

no2\_val

</td>

<td>

1.031

</td>

<td>

1.014

</td>

<td>

1.047

</td>

<td>

0.001

</td>

<td>

no2\_val

</td>

<td>

p-value \< 0.05

</td>

</tr>

<tr>

<td style="text-align:left">

o3\_val

</td>

<td>

0.949

</td>

<td>

0.928

</td>

<td>

0.970

</td>

<td>

0.00000

</td>

<td>

o3\_val

</td>

<td>

p-value \< 0.05

</td>

</tr>

<tr>

<td style="text-align:left">

so2\_val

</td>

<td>

1.267

</td>

<td>

1.103

</td>

<td>

1.459

</td>

<td>

0.002

</td>

<td>

so2\_val

</td>

<td>

p-value \< 0.05

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

</table>

### Summary tables of models in Aim 2

Add <br> \#\#\#\# Supplementary Table 4. Effect of air pollutants on the
number of COVID-related cases at the subregional level Add <br>

``` r
stargazer(pm25_cases.nb, pm10_cases.nb, nox_cases.nb, no2_cases.nb,
          o3_cases.nb, so2_cases.nb,type="html",out = "fig_out_v4/LA_covid_allPoll_CASES_nb.html",
          dep.var.labels="Number of COVID-19-related cases",
          single.row=TRUE)
```

<table style="text-align:center">

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="6">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="6" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="6">

Number of COVID-19-related cases

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

(1)

</td>

<td>

(2)

</td>

<td>

(3)

</td>

<td>

(4)

</td>

<td>

(5)

</td>

<td>

(6)

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

X2018.people.per.sq..km

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

</tr>

<tr>

<td style="text-align:left">

Mean\_ann\_earnings

</td>

<td>

0.00001 (0.00001)

</td>

<td>

0.00001 (0.00001)

</td>

<td>

\-0.00000 (0.00001)

</td>

<td>

\-0.00000 (0.00001)

</td>

<td>

0.00001 (0.00001)

</td>

<td>

0.00000 (0.00001)

</td>

</tr>

<tr>

<td style="text-align:left">

median\_age\_2018

</td>

<td>

\-0.091<sup>\*\*\*</sup> (0.010)

</td>

<td>

\-0.087<sup>\*\*\*</sup> (0.010)

</td>

<td>

\-0.058<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.055<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.069<sup>\*\*\*</sup> (0.009)

</td>

<td>

\-0.060<sup>\*\*\*</sup> (0.010)

</td>

</tr>

<tr>

<td style="text-align:left">

pm25\_val

</td>

<td>

\-0.096<sup>\*\*\*</sup> (0.028)

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

pm10\_val

</td>

<td>

</td>

<td>

\-0.067<sup>\*\*\*</sup> (0.018)

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

nox\_val

</td>

<td>

</td>

<td>

</td>

<td>

0.017<sup>\*\*\*</sup> (0.005)

</td>

<td>

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

no2\_val

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

0.030<sup>\*\*\*</sup> (0.009)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

o3\_val

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

\-0.053<sup>\*\*\*</sup> (0.011)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

so2\_val

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

0.237<sup>\*\*\*</sup> (0.075)

</td>

</tr>

<tr>

<td style="text-align:left">

Constant

</td>

<td>

9.593<sup>\*\*\*</sup> (0.557)

</td>

<td>

9.513<sup>\*\*\*</sup> (0.522)

</td>

<td>

7.409<sup>\*\*\*</sup> (0.495)

</td>

<td>

7.190<sup>\*\*\*</sup> (0.531)

</td>

<td>

8.206<sup>\*\*\*</sup> (0.426)

</td>

<td>

7.239<sup>\*\*\*</sup> (0.531)

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Observations

</td>

<td>

295

</td>

<td>

295

</td>

<td>

295

</td>

<td>

295

</td>

<td>

295

</td>

<td>

295

</td>

</tr>

<tr>

<td style="text-align:left">

Log Likelihood

</td>

<td>

\-1,786.872

</td>

<td>

\-1,785.649

</td>

<td>

\-1,787.182

</td>

<td>

\-1,786.734

</td>

<td>

\-1,782.057

</td>

<td>

\-1,787.432

</td>

</tr>

<tr>

<td style="text-align:left">

theta

</td>

<td>

2.725<sup>\*\*\*</sup> (0.216)

</td>

<td>

2.746<sup>\*\*\*</sup> (0.218)

</td>

<td>

2.722<sup>\*\*\*</sup> (0.216)

</td>

<td>

2.729<sup>\*\*\*</sup> (0.216)

</td>

<td>

2.808<sup>\*\*\*</sup> (0.223)

</td>

<td>

2.717<sup>\*\*\*</sup> (0.215)

</td>

</tr>

<tr>

<td style="text-align:left">

Akaike Inf. Crit.

</td>

<td>

3,583.745

</td>

<td>

3,581.297

</td>

<td>

3,584.364

</td>

<td>

3,583.468

</td>

<td>

3,574.114

</td>

<td>

3,584.864

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="6" style="text-align:right">

<sup>*</sup>p\<0.1; <sup>**</sup>p\<0.05; <sup>***</sup>p\<0.01

</td>

</tr>

</table>

<br> The value in parentheses represent the standard error. <br> Summary
of the effects of individual air pollutant species on lab-confirmed
COVID-related cases up to April 10th 2020. The three asterisks near the
theta indicate that the models are significantly better than null.

#### Supplementary Table 5. Effect of air pollutants on the number of COVID-related deaths at the subregional level

Add <br>

``` r
stargazer(pm25_deaths.nb, pm10_deaths.nb, nox_deaths.nb, no2_deaths.nb,
          o3_deaths.nb, so2_deaths.nb,type="html",out = "fig_out_v4/LA_covid_allPoll_nb.html",
          dep.var.labels="Number of COVID-19-related deaths",
          single.row=TRUE)
```

<table style="text-align:center">

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="6">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="6" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="6">

Number of COVID-19-related deaths

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

(1)

</td>

<td>

(2)

</td>

<td>

(3)

</td>

<td>

(4)

</td>

<td>

(5)

</td>

<td>

(6)

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

X2018.people.per.sq..km

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

</tr>

<tr>

<td style="text-align:left">

Mean\_ann\_earnings

</td>

<td>

0.00000 (0.00001)

</td>

<td>

0.00001 (0.00001)

</td>

<td>

\-0.00000 (0.00001)

</td>

<td>

\-0.00000 (0.00001)

</td>

<td>

0.00001 (0.00001)

</td>

<td>

0.00001 (0.00001)

</td>

</tr>

<tr>

<td style="text-align:left">

median\_age\_2018

</td>

<td>

\-0.079<sup>\*\*\*</sup> (0.012)

</td>

<td>

\-0.079<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.055<sup>\*\*\*</sup> (0.012)

</td>

<td>

\-0.050<sup>\*\*\*</sup> (0.012)

</td>

<td>

\-0.067<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.058<sup>\*\*\*</sup> (0.012)

</td>

</tr>

<tr>

<td style="text-align:left">

pm25\_val

</td>

<td>

\-0.039 (0.031)

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

pm10\_val

</td>

<td>

</td>

<td>

\-0.036<sup>\*</sup> (0.020)

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

nox\_val

</td>

<td>

</td>

<td>

</td>

<td>

0.016<sup>\*\*\*</sup> (0.006)

</td>

<td>

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

no2\_val

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

0.030<sup>\*\*\*</sup> (0.010)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

o3\_val

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

\-0.037<sup>\*\*\*</sup> (0.013)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

so2\_val

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

0.202<sup>\*\*</sup> (0.083)

</td>

</tr>

<tr>

<td style="text-align:left">

Constant

</td>

<td>

6.988<sup>\*\*\*</sup> (0.635)

</td>

<td>

7.103<sup>\*\*\*</sup> (0.597)

</td>

<td>

5.541<sup>\*\*\*</sup> (0.542)

</td>

<td>

5.268<sup>\*\*\*</sup> (0.580)

</td>

<td>

6.323<sup>\*\*\*</sup> (0.489)

</td>

<td>

5.498<sup>\*\*\*</sup> (0.596)

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Observations

</td>

<td>

306

</td>

<td>

306

</td>

<td>

306

</td>

<td>

306

</td>

<td>

306

</td>

<td>

306

</td>

</tr>

<tr>

<td style="text-align:left">

Log Likelihood

</td>

<td>

\-1,395.565

</td>

<td>

\-1,394.788

</td>

<td>

\-1,392.433

</td>

<td>

\-1,391.427

</td>

<td>

\-1,392.320

</td>

<td>

\-1,393.427

</td>

</tr>

<tr>

<td style="text-align:left">

theta

</td>

<td>

2.099<sup>\*\*\*</sup> (0.171)

</td>

<td>

2.109<sup>\*\*\*</sup> (0.171)

</td>

<td>

2.143<sup>\*\*\*</sup> (0.175)

</td>

<td>

2.157<sup>\*\*\*</sup> (0.176)

</td>

<td>

2.142<sup>\*\*\*</sup> (0.175)

</td>

<td>

2.129<sup>\*\*\*</sup> (0.174)

</td>

</tr>

<tr>

<td style="text-align:left">

Akaike Inf. Crit.

</td>

<td>

2,801.130

</td>

<td>

2,799.577

</td>

<td>

2,794.865

</td>

<td>

2,792.853

</td>

<td>

2,794.640

</td>

<td>

2,796.854

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="6" style="text-align:right">

<sup>*</sup>p\<0.1; <sup>**</sup>p\<0.05; <sup>***</sup>p\<0.01

</td>

</tr>

</table>

<br> Summary of the effects of individual air pollutant species on
lab-confirmed COVID-related deaths up to April 10th 2020. The three
asterisks near the theta indicate that the models are significantly
better than null. The MRR (mortality rate ratios) of O3 is interesting:
more O3 = \~3% less deaths. NOx, NO2 and SO2 are all significant
predictors; they will increase the number of deaths by 1.3, 2.4 and
17.2% respectively, per 1 µg/m3.

## Aim 3: Individual level analysis

### Data curation

#### covid data curation

``` r
covid_df = read.csv("data_v4/covid19_result.txt", sep = '\t')
covid_id = covid_df[,c(1,5,6)]
# I didnt realise that there were repeated values!!!!!
colnames(covid_id)
```

    ## [1] "eid"    "origin" "result"

``` r
covid_id_unique = aggregate(covid_id, by=list(covid_id$eid), 
                            FUN=max)

ggplot(data = covid_id_unique, aes(x = result))+
  geom_histogram(stat="count")
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

``` r
ggplot(data = covid_id_unique, aes(x = origin))+
  geom_histogram(stat="count")
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-31-2.png)<!-- -->
Look at time

``` r
covid_df_time = covid_df
covid_df_time$specdate = as.Date(covid_df_time$specdate, format="%d/%m/%Y")
covid_df_time<-covid_df_time %>% 
group_by(specdate) %>%
mutate(culm_cases=cumsum(result))

ggplot(covid_df_time, aes(x=specdate, y=culm_cases)) +
  geom_bar(stat="identity") + 
  scale_x_date(date_labels = "%d/%m/%Y")
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

#### load non-covid UKB phenotype data

Data curation

``` r
### load non-covid UKB phenotype data ----
ukb_covid = read.csv("data_v4/29_4_2020_ukb41646_covid19_subset.csv")[,-1]
#check distribution

#here, I will aggregate the columns by selecting the last results, for each 
# reshape
ukb_covid_t = melt(ukb_covid, id='eid')
ukb_covid_t <- ukb_covid_t[order(ukb_covid_t$variable),] 
#delete everything after period 
ukb_covid_t$variable = gsub("\\..*","",ukb_covid_t$variable)

#aggregate by last 
ukb_covid_t_na = na.omit(ukb_covid_t)
ukb_covid_t_na = aggregate(ukb_covid_t_na, by=list(ukb_covid_t_na$eid,ukb_covid_t_na$variable), FUN=last)

# put it back straight 
ukb_covid_t_na$variable = as.factor(ukb_covid_t_na$variable)
ukb_covid_t_na = ukb_covid_t_na[,-c(1,2)]
ukb_covid_cur = cast(ukb_covid_t_na, eid~variable)

### label UKB data ----
#get levels and labels 
lvl.0493 <- c(-141,-131,-121,0)
lbl.0493 <- c("Often","Sometimes","Do not know","Rarely/never")

lvl.0007 <- c(0,1)
lbl.0007 <- c("No","Yes")

lvl.100349 <- c(-3,-1,0,1)
lbl.100349 <- c("Prefer not to answer","Do not know","No","Yes")

lvl.0090 <- c(-3,0,1,2)
lbl.0090 <- c("Prefer not to answer","Never","Previous","Current")

lvl.0009 <- c(0,1)
lbl.0009 <- c("Female","Male")

### eid replacement, saved here ----
yy_replace <- function(string, patterns, replacements) {
  for (i in seq_along(patterns))
    string <- gsub(patterns[i], replacements[i], string, perl=TRUE)
  string
}
#label = replacement, lvl = to replace

ukb_covid_cur$X22609 <- yy_replace(ukb_covid_cur$X22609, lvl.0493, lbl.0493)
ukb_covid_cur$X22610 <- yy_replace(ukb_covid_cur$X22610, lvl.0493, lbl.0493)
ukb_covid_cur$X22611 <- yy_replace(ukb_covid_cur$X22611, lvl.0493, lbl.0493)
ukb_covid_cur$X22615 <- yy_replace(ukb_covid_cur$X22615, lvl.0493, lbl.0493)
ukb_covid_cur$X22616 <- yy_replace(ukb_covid_cur$X22616, lvl.0007, lbl.0007)
ukb_covid_cur$X2316 <- yy_replace(ukb_covid_cur$X2316, lvl.100349, lbl.100349)
ukb_covid_cur$X31 <- yy_replace(ukb_covid_cur$X31, lvl.0009, lbl.0009)
ukb_covid_cur$X22130 <- yy_replace(ukb_covid_cur$X22130, lvl.0007, lbl.0007)
ukb_covid_cur$X2443 <- yy_replace(ukb_covid_cur$X2443, lvl.100349, lbl.100349)
ukb_covid_cur$X20116 <- yy_replace(ukb_covid_cur$X20116, lvl.0090, lbl.0090)


### put custom column names----
#get the ID that I designed
UID_names = read.csv("data_v4/UKB_list_columns_AIRcovid.csv")[c("FieldID","my_colname")]
UID_names$FieldID = paste("X",UID_names$FieldID, sep = "")
UID_names$my_colname = lapply(UID_names$my_colname,toString)
UID_names[nrow(UID_names) + 1,] = c("eid","eid")

# sort it according to the dataset

UID_names_match = colnames(ukb_covid_cur)

UID_sort_df = data.frame(UID_names_match)
colnames(UID_sort_df) <- "FieldID"

#merge to sort
UID_sort_df = merge(UID_sort_df, UID_names, by = "FieldID")
# this is now sorted according to the correct df

#replace column names
ukb_covid_df = ukb_covid_cur

colnames(ukb_covid_df) <- UID_sort_df$my_colname

### merge both ukb datasets together ----
ukb_covid_merge = merge(covid_id_unique[,-1], ukb_covid_df, by = "eid")
nrow(ukb_covid_merge)
```

    ## [1] 1474

further processing…

``` r
# do 2020 since COVID occurred in this year
ukb_covid_merge$age = 2020 - ukb_covid_merge$birthYear
ukb_covid_merge$whr = ukb_covid_merge$waist / ukb_covid_merge$hip

#definition from the NHS: Diastolic >= 90 OR Systolic >=140
ukb_covid_merge$highBP[ukb_covid_merge$diaBP>=90 | ukb_covid_merge$sysBP>=140] <-1
ukb_covid_merge$highBP[ukb_covid_merge$diaBP<90 | ukb_covid_merge$sysBP<140] <-0

# add a column: COVID + and inpatient
ukb_covid_merge$inpatient_covid  <-0
ukb_covid_merge$inpatient_covid [ukb_covid_merge$origin == 1 & ukb_covid_merge$origin == 1] <-1

write.csv(ukb_covid_merge, "data_output_v4/processed_29_4_2020_ukb41646_covid19.csv")
```

#### Preliminary analysis with PM2.5

I will load the air pollutants & match the covid patients to their
nearest pollution climate mapping location.

``` r
### load pm2.5 data ----
#note the pop weighted data have the area code, but not the x y coords

pm25 = read.csv("data_v4/processed_PM25_uk-air_annual_mean_mappm252018g.csv", na.strings = "MISSING")[c("x","y","pm252018g")]
nrow(pm25)
```

    ## [1] 281802

``` r
pm25_na = na.omit(pm25)
nrow(pm25_na)
```

    ## [1] 254877

``` r
### convert to lon lat - pm2.5 ----

# transform: easting/northing -> long/lat
### shortcuts (from https://stephendavidgregory.github.io/useful/UKgrid_to_LatLon)
ukgrid <- "+init=epsg:27700"
latlong <- "+init=epsg:4326"

### Create coordinates variable
pm25_coords <- cbind(Easting = as.numeric(as.character(pm25_na$x)),
                Northing = as.numeric(as.character(pm25_na$y)))

### Create the SpatialPointsDataFrame
pm25_SP <- SpatialPointsDataFrame(pm25_coords,
                                 data = pm25_na,
                                 proj4string = CRS("+init=epsg:27700"))

### Convert
pm25_LL <- spTransform(pm25_SP, CRS(latlong))

pm25_ll_df = data.frame('pm25_lon' = coordinates(pm25_LL)[, 1], 'pm25_lat' = coordinates(pm25_LL)[, 2], 'pm25_val' = pm25_LL$pm252018g)


write.csv(pm25_ll_df,"data_output_v4/processed_pm25_lonlat.csv")
```

#### convert UKB data to long lat

``` r
# convert ukb data eastings northing to lon lat 

ukb_covid_coord_na = na.omit(ukb_covid_merge[c("eid","x_coord","y_coord")])
nrow(ukb_covid_coord_na)
```

    ## [1] 1464

``` r
ukb_covid_coords <- cbind(Easting = as.numeric(as.character(ukb_covid_coord_na$x_coord)),
                     Northing = as.numeric(as.character(ukb_covid_coord_na$y_coord)))
### Create the SpatialPointsDataFrame
ukb_covid_SP <- SpatialPointsDataFrame(ukb_covid_coords,
                                  data = ukb_covid_coord_na,
                                  proj4string = CRS("+init=epsg:27700"))

### Convert
ukb_covid_ll <- spTransform(ukb_covid_SP, CRS(latlong))

ukb_covid_ll_df = data.frame('ukb_lon' = coordinates(ukb_covid_ll)[, 1], 
                             'ukb_lat' = coordinates(ukb_covid_ll)[, 2], 
                             'eid' = ukb_covid_ll$eid)
write.csv(ukb_covid_ll_df, ("data_output_v4/ukb_covid_lonlat_df.csv"))
# plot UKB locations
world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf() +
  geom_point(data = ukb_covid_ll_df, aes(x = ukb_lon, y = ukb_lat), size = 1, 
             shape = 23, fill = "darkred") +
  coord_sf(ylim = c(min(ukb_covid_ll_df$ukb_lat)-2, max(ukb_covid_ll_df$ukb_lat)+4), 
           xlim = c(min(ukb_covid_ll_df$ukb_lon)-4, max(ukb_covid_ll_df$ukb_lon)+3), expand = FALSE)
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

``` r
ggsave('fig_out_v4/UKB_COVID_participants_locations.pdf')
```

**I wrote a python script: match\_air\_to\_UKB\_covid**

### Merge and QC data for PM2.5

``` r
ukb_eid_pm25_merged = read.csv('data_output_v4/merged_ukb_pm25.csv')[c("eid","distance","pm25_val")]

ukb_covid_pm25_df = merge(ukb_covid_merge, ukb_eid_pm25_merged, by = 'eid')
colnames(ukb_covid_pm25_df)
```

    ##  [1] "eid"             "origin"          "result"          "n_cancers"      
    ##  [5] "townsend"        "x_coord"         "y_coord"         "smoking"        
    ##  [9] "fev"             "copd"            "WP_dusty"        "WP_chemicals"   
    ## [13] "WP_cig"          "WP_diesel"       "breathing"       "whistling"      
    ## [17] "diabetes"        "sex"             "birthYear"       "diaBP"          
    ## [21] "sysBP"           "waist"           "hip"             "height"         
    ## [25] "age"             "whr"             "highBP"          "inpatient_covid"
    ## [29] "distance"        "pm25_val"

``` r
#omit those without coordinates in the ukb data 

### merge pop density----
ukb_cities_eid = read.csv("data_output_v4/2_5_2020_ukb_eid_match_cities.csv")[c("eid","spec_area","gen_area")]
nrow(ukb_cities_eid)
```

    ## [1] 1464

``` r
popDens_2018 = read.csv("data_v4/2018_official_popDensity.csv")[c("Name","X2018.people.per.sq..km")]
popDens_2018$X2018.people.per.sq..km = gsub(",", "", popDens_2018$X2018.people.per.sq..km)
popDens_2018$X2018.people.per.sq..km = as.numeric(popDens_2018$X2018.people.per.sq..km)

#Make everything lower case
ukb_cities_eid$spec_area = tolower(ukb_cities_eid$spec_area)
ukb_cities_eid$gen_area = tolower(ukb_cities_eid$gen_area)
popDens_2018$Name = tolower(popDens_2018$Name)

#delete repeated strings
ukb_cities_eid$spec_area = gsub("london borough of ", "", ukb_cities_eid$spec_area)
ukb_cities_eid$spec_area = gsub("royal borough of ", "", ukb_cities_eid$spec_area)
ukb_cities_eid$spec_area = gsub("city of ", "", ukb_cities_eid$spec_area)
ukb_cities_eid$spec_area = gsub("st helens", "helens", ukb_cities_eid$spec_area)
ukb_cities_eid$gen_area = gsub(" combined authority", "", ukb_cities_eid$gen_area)
popDens_2018$Name = gsub(", city of", "", popDens_2018$Name)
popDens_2018$Name = gsub("st. ", "", popDens_2018$Name)
#delete west midlands because its redundant (all cities are well defined) and there are several regions called west midlands so it's confusing
popDens_2018$Name = gsub("west midlands", "deleted_west_midlands_because_redundant", popDens_2018$Name)
popDens_2018$Name = gsub("\\s*\\([^\\)]+\\)","",popDens_2018$Name)

ukb_cities_eid_merged = merge(ukb_cities_eid, popDens_2018, by.x="spec_area",by.y="Name", all.x=TRUE)
colnames(ukb_cities_eid_merged)[4] <- "spec_Pop_dens_perkm2"
nrow(ukb_cities_eid_merged)
```

    ## [1] 1464

``` r
#check for duplicated rows in popDens
length(unique(popDens_2018$Name))
```

    ## [1] 434

``` r
length(popDens_2018$Name) # there is a dupicate
```

    ## [1] 435

``` r
popDens_2018$Name[duplicated(popDens_2018$Name)]
```

    ## [1] "deleted_west_midlands_because_redundant"

``` r
ukb_cities_eid_merged2 = merge(ukb_cities_eid_merged, popDens_2018, by.x="gen_area",by.y="Name", all.x=TRUE,all.y=FALSE)
colnames(ukb_cities_eid_merged2)[5] <- "gen_Pop_dens_perkm2"
nrow(ukb_cities_eid_merged2)
```

    ## [1] 1464

``` r
#QC
sum(is.na(ukb_cities_eid_merged2$spec_Pop_dens_perkm2))
```

    ## [1] 214

``` r
sum(is.na(ukb_cities_eid_merged2$gen_Pop_dens_perkm2))
```

    ## [1] 516

``` r
ukb_cities_eid_merged2$merged_popDens_km2 <- ifelse(is.na(ukb_cities_eid_merged2$spec_Pop_dens_perkm2), 
                                                    ukb_cities_eid_merged2$gen_Pop_dens_perkm2, ukb_cities_eid_merged2$spec_Pop_dens_perkm2)
ukb_cities_eid_out = ukb_cities_eid_merged2[c("eid","merged_popDens_km2")]
ukb_covid_pm25_popDens_df = merge(ukb_covid_pm25_df, ukb_cities_eid_out, by = 'eid')
nrow(ukb_covid_pm25_df)
```

    ## [1] 1464

``` r
ukb_covid_pm25_popDens_df = read.csv("data_output_v4/2_5_2020_full_cov_CV_UKB_air_dataset.csv")
ukb_covid_pm25_popDens_df$smoker = ifelse(ukb_covid_pm25_popDens_df$smoking =="Current", "smoking", "not_smoking")

library(MASS)
summary(ukb_covid_pm25.nb <-glm.nb(data = ukb_covid_pm25_popDens_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                               whr + sex + age + pm25_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + pm25_val, 
    ##     data = ukb_covid_pm25_popDens_df, init.theta = 15229.88441, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.2059  -0.9479  -0.8028   0.6589   1.1708  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                   -2.065e+00  8.543e-01  -2.418  0.01562 * 
    ## merged_popDens_km2            -9.956e-06  1.884e-05  -0.529  0.59711   
    ## n_cancers                     -1.258e-01  1.272e-01  -0.989  0.32253   
    ## townsend                       1.354e-02  1.310e-02   1.034  0.30121   
    ## smokingNever                   3.665e-01  1.379e-01   2.658  0.00785 **
    ## smokingPrefer not to answer    7.608e-01  4.478e-01   1.699  0.08929 . 
    ## smokingPrevious                4.056e-01  1.377e-01   2.946  0.00322 **
    ## whistlingNo                    2.671e-01  2.698e-01   0.990  0.32220   
    ## whistlingPrefer not to answer  4.568e-02  1.276e+00   0.036  0.97145   
    ## whistlingYes                   2.533e-01  2.742e-01   0.924  0.35559   
    ## diabetesNo                    -1.619e-01  6.106e-01  -0.265  0.79085   
    ## diabetesPrefer not to answer  -6.877e-01  1.433e+00  -0.480  0.63140   
    ## diabetesYes                   -1.956e-01  6.195e-01  -0.316  0.75215   
    ## whr                            7.856e-01  5.839e-01   1.346  0.17846   
    ## sexMale                        3.340e-02  1.000e-01   0.334  0.73838   
    ## age                           -5.681e-03  4.748e-03  -1.196  0.23153   
    ## pm25_val                       6.035e-02  2.866e-02   2.106  0.03520 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(15229.88) family taken to be 1)
    ## 
    ##     Null deviance: 1040.2  on 1449  degrees of freedom
    ## Residual deviance: 1015.6  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 2365.6
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  15230 
    ##           Std. Err.:  58208 
    ## Warning while fitting theta: iteration limit reached 
    ## 
    ##  2 x log-likelihood:  -2329.646

``` r
summary(ukb_covid_pm25.p <-glm(data = ukb_covid_pm25_popDens_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                     whr + sex + age + pm25_val, family = 'poisson'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + pm25_val, 
    ##     family = "poisson", data = ukb_covid_pm25_popDens_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.2060  -0.9479  -0.8028   0.6589   1.1708  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                   -2.065e+00  8.543e-01  -2.418  0.01562 * 
    ## merged_popDens_km2            -9.956e-06  1.883e-05  -0.529  0.59711   
    ## n_cancers                     -1.258e-01  1.272e-01  -0.989  0.32252   
    ## townsend                       1.354e-02  1.310e-02   1.034  0.30120   
    ## smokingNever                   3.665e-01  1.379e-01   2.658  0.00785 **
    ## smokingPrefer not to answer    7.609e-01  4.478e-01   1.699  0.08928 . 
    ## smokingPrevious                4.056e-01  1.377e-01   2.946  0.00321 **
    ## whistlingNo                    2.671e-01  2.698e-01   0.990  0.32219   
    ## whistlingPrefer not to answer  4.567e-02  1.276e+00   0.036  0.97145   
    ## whistlingYes                   2.533e-01  2.742e-01   0.924  0.35558   
    ## diabetesNo                    -1.619e-01  6.105e-01  -0.265  0.79084   
    ## diabetesPrefer not to answer  -6.877e-01  1.433e+00  -0.480  0.63138   
    ## diabetesYes                   -1.956e-01  6.194e-01  -0.316  0.75214   
    ## whr                            7.856e-01  5.839e-01   1.346  0.17845   
    ## sexMale                        3.340e-02  1.000e-01   0.334  0.73837   
    ## age                           -5.681e-03  4.748e-03  -1.196  0.23151   
    ## pm25_val                       6.035e-02  2.866e-02   2.106  0.03520 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 1040.2  on 1449  degrees of freedom
    ## Residual deviance: 1015.6  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 2363.6
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
#I dont' think that pop density should be an offset here
#check for multicolinearity
car::vif(ukb_covid_pm25.nb)
```

    ##                        GVIF Df GVIF^(1/(2*Df))
    ## merged_popDens_km2 2.443689  1        1.563230
    ## n_cancers          1.019337  1        1.009622
    ## townsend           1.442148  1        1.200895
    ## smoking            1.259363  3        1.039182
    ## whistling          1.875826  3        1.110535
    ## diabetes           2.049216  3        1.127019
    ## whr                1.850854  1        1.360461
    ## sex                1.619879  1        1.272744
    ## age                1.138711  1        1.067104
    ## pm25_val           2.011763  1        1.418366

``` r
# according to http://www.jstor.org/stable/2290467, GVIF^(1/(2*Df)) < 5 is ok 
# model reduction: 

summary(ukb_covid_pm25.nb_red <-glm.nb(data = ukb_covid_pm25_popDens_df, result ~ merged_popDens_km2 + smoking +pm25_val))
```

    ## 
    ## Call:
    ## glm.nb(formula = result ~ merged_popDens_km2 + smoking + pm25_val, 
    ##     data = ukb_covid_pm25_popDens_df, init.theta = 15585.93096, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.1217  -0.9499  -0.8075   0.6755   1.0620  
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                 -1.639e+00  2.629e-01  -6.235 4.52e-10 ***
    ## merged_popDens_km2          -9.261e-07  1.680e-05  -0.055  0.95604    
    ## smokingNever                 3.369e-01  1.344e-01   2.507  0.01216 *  
    ## smokingPrefer not to answer  6.113e-01  4.261e-01   1.435  0.15139    
    ## smokingPrevious              3.899e-01  1.353e-01   2.882  0.00395 ** 
    ## pm25_val                     5.801e-02  2.833e-02   2.048  0.04055 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(15585.93) family taken to be 1)
    ## 
    ##     Null deviance: 1049.0  on 1461  degrees of freedom
    ## Residual deviance: 1032.7  on 1456  degrees of freedom
    ##   (2 observations deleted due to missingness)
    ## AIC: 2370.7
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  15586 
    ##           Std. Err.:  60471 
    ## Warning while fitting theta: iteration limit reached 
    ## 
    ##  2 x log-likelihood:  -2356.739

``` r
summary(ukb_covid_pm25.p_red <-glm(data = ukb_covid_pm25_popDens_df, result ~ merged_popDens_km2 + smoking + pm25_val, family = 'poisson'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + smoking + pm25_val, 
    ##     family = "poisson", data = ukb_covid_pm25_popDens_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.1218  -0.9499  -0.8075   0.6756   1.0620  
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                 -1.639e+00  2.629e-01  -6.235 4.52e-10 ***
    ## merged_popDens_km2          -9.261e-07  1.680e-05  -0.055  0.95604    
    ## smokingNever                 3.369e-01  1.344e-01   2.508  0.01216 *  
    ## smokingPrefer not to answer  6.113e-01  4.261e-01   1.435  0.15138    
    ## smokingPrevious              3.899e-01  1.353e-01   2.882  0.00395 ** 
    ## pm25_val                     5.801e-02  2.832e-02   2.048  0.04055 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 1049.0  on 1461  degrees of freedom
    ## Residual deviance: 1032.7  on 1456  degrees of freedom
    ##   (2 observations deleted due to missingness)
    ## AIC: 2368.7
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
#check 
anova(ukb_covid_pm25.nb, ukb_covid_pm25.nb_red) # p < 0.05 -> not significantly different 
```

    ## Likelihood ratio tests of Negative Binomial Models
    ## 
    ## Response: result
    ##                                                                                                     Model
    ## 1                                                                 merged_popDens_km2 + smoking + pm25_val
    ## 2 merged_popDens_km2 + n_cancers + townsend + smoking + whistling + diabetes + whr + sex + age + pm25_val
    ##      theta Resid. df    2 x log-lik.   Test    df LR stat.   Pr(Chi)
    ## 1 15585.93      1456       -2356.739                                
    ## 2 15229.88      1433       -2329.646 1 vs 2    23  27.0924 0.2520556

Conclusion: <br> 1. PM2.5 is a significant predictor of COVID-19: more
cases with PM2.5 levels are high <br> 2. Smoking interpretation <br> The
probability of infection of people who never smoked is the expected
difference in log count between reference group -\> 1.442676 \* rate of
current smoker<br> That of previous smokers = 1.5001871 \* rate of
current smoker <br> 3. In the actual analysis, I can explore using
binary models <br>

### Convert the rest of the pollutant variables from X Y to long lat

``` r
### convert NO2 x y to lat lon ----
### shortcuts (from https://stephendavidgregory.github.io/useful/UKgrid_to_LatLon)
no2_raw_dt = na.omit(read.csv("data_v4/processed_30_4_2020_NO2_map2018g.csv", na.strings = "MISSING")[c("x","y","no22018")])

ukgrid <- "+init=epsg:27700"
latlong <- "+init=epsg:4326"

### Create coordinates variable
no2_coords <- cbind(Easting = as.numeric(as.character(no2_raw_dt$x)),
                     Northing = as.numeric(as.character(no2_raw_dt$y)))

no2_LL <- spTransform(SpatialPointsDataFrame(no2_coords,
                                  data = no2_raw_dt,
                                  proj4string = CRS("+init=epsg:27700")), CRS(latlong))

no2_LL_df = data.frame('no2_lon' = coordinates(no2_LL)[, 1], 'no2_lat' = coordinates(no2_LL)[, 2], 'no2_val' = no2_LL$no22018)

write.csv(no2_LL_df,"data_output_v4/processed_no2_lonlat.csv")

### convert SO2 x y to lat lon ----
so2_raw_dt = na.omit(read.csv("data_v4/processed_30_4_2020_SO2_map2018g.csv", na.strings = "MISSING")[c("x","y","so22018")])

so2_coords <- cbind(Easting = as.numeric(as.character(so2_raw_dt$x)),
                    Northing = as.numeric(as.character(so2_raw_dt$y)))

so2_LL <- spTransform(SpatialPointsDataFrame(so2_coords,
                                             data = so2_raw_dt,
                                             proj4string = CRS("+init=epsg:27700")), CRS(latlong))

so2_LL_df = data.frame('so2_lon' = coordinates(so2_LL)[, 1], 'so2_lat' = coordinates(so2_LL)[, 2], 'so2_val' = so2_LL$so22018)

write.csv(so2_LL_df,"data_output_v4/processed_so2_lonlat.csv")

### convert 03 x y to lat lon ----

o3_raw_dt = na.omit(read.csv("data_v4/processed_30_4_2020_O3_map2018g.csv", na.strings = "MISSING")[c("x","y","dgt12018")])

o3_coords <- cbind(Easting = as.numeric(as.character(o3_raw_dt$x)),
                    Northing = as.numeric(as.character(o3_raw_dt$y)))

o3_LL <- spTransform(SpatialPointsDataFrame(o3_coords,
                                             data = o3_raw_dt,
                                             proj4string = CRS("+init=epsg:27700")), CRS(latlong))

o3_LL_df = data.frame('o3_lon' = coordinates(o3_LL)[, 1], 'o3_lat' = coordinates(o3_LL)[, 2], 'o3_val' = o3_LL$dgt12018)

write.csv(o3_LL_df,"data_output_v4/processed_o3_lonlat.csv")

### convert PM10 x y to lat lon ----

pm10_raw_dt = na.omit(read.csv("data_v4/processed_30_4_2020_PM10_mappm102018g.csv", na.strings = "MISSING")[c("x","y","pm102018g")])

pm10_coords <- cbind(Easting = as.numeric(as.character(pm10_raw_dt$x)),
                   Northing = as.numeric(as.character(pm10_raw_dt$y)))

pm10_LL <- spTransform(SpatialPointsDataFrame(pm10_coords,
                                            data = pm10_raw_dt,
                                            proj4string = CRS("+init=epsg:27700")), CRS(latlong))

pm10_LL_df = data.frame('pm10_lon' = coordinates(pm10_LL)[, 1], 'pm10_lat' = coordinates(pm10_LL)[, 2], 'pm10_val' = pm10_LL$pm102018g)

write.csv(pm10_LL_df,"data_output_v4/processed_pm10_lonlat.csv")


### convert NOx x y to lat lon ----

nox_raw_dt = na.omit(read.csv("data_v4/processed_30_4_2020_NOX_map2018g.csv", na.strings = "MISSING")[c("x","y","nox2018")])

nox_coords <- cbind(Easting = as.numeric(as.character(nox_raw_dt$x)),
                     Northing = as.numeric(as.character(nox_raw_dt$y)))

nox_LL <- spTransform(SpatialPointsDataFrame(nox_coords,
                                              data = nox_raw_dt,
                                              proj4string = CRS("+init=epsg:27700")), CRS(latlong))

nox_LL_df = data.frame('nox_lon' = coordinates(nox_LL)[, 1], 'nox_lat' = coordinates(nox_LL)[, 2], 'nox_val' = nox_LL$nox2018)

write.csv(nox_LL_df,"data_output_v4/processed_nox_lonlat.csv")
```

### Analysis

load and merge all data

``` r
no2_ukb = read.csv("data_output_v4/merged_ukb_no2.csv")[c('eid','no2_val')]
so2_ukb = read.csv("data_output_v4/merged_ukb_so2.csv")[c('eid','so2_val')]
o3_ukb = read.csv("data_output_v4/merged_ukb_o3.csv")[c('eid','o3_val')]
nox_ukb = read.csv("data_output_v4/merged_ukb_nox.csv")[c('eid','nox_val')]
pm10_ukb = read.csv("data_output_v4/merged_ukb_pm10.csv")[c('eid','pm10_val')]

ukb_additional_pol = merge(no2_ukb, so2_ukb, by = 'eid')
ukb_additional_pol = merge(ukb_additional_pol, o3_ukb, by = 'eid')
ukb_additional_pol = merge(ukb_additional_pol, nox_ukb, by = 'eid')
ukb_additional_pol = merge(ukb_additional_pol, pm10_ukb, by = 'eid')

ukb_covid_allPol_df = merge(ukb_covid_pm25_popDens_df, ukb_additional_pol, by = 'eid')
```

#### Descriptive statistics of the UK Biobank data

``` r
library(arsenal)
ukb_descript_stats <- tableby(result ~ ., data = ukb_covid_allPol_df[,3:ncol(ukb_covid_allPol_df)])
summary(ukb_descript_stats, title = "Descriptive statistics of the UK Biobank data")
```

    ## 
    ## Table: Descriptive statistics of the UK Biobank data
    ## 
    ## |                                       |        0 (N=800)        |        1 (N=664)        |     Total (N=1464)      | p value|
    ## |:--------------------------------------|:-----------------------:|:-----------------------:|:-----------------------:|-------:|
    ## |**origin**                             |                         |                         |                         | < 0.001|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |      0.672 (0.470)      |      0.858 (0.349)      |      0.757 (0.429)      |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |      0.000 - 1.000      |      0.000 - 1.000      |      0.000 - 1.000      |        |
    ## |**n_cancers**                          |                         |                         |                         |   0.173|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |            0            |            1            |            1            |        |
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |      0.119 (0.350)      |      0.095 (0.309)      |      0.108 (0.332)      |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |      0.000 - 2.000      |      0.000 - 2.000      |      0.000 - 2.000      |        |
    ## |**townsend**                           |                         |                         |                         |   0.042|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |     -0.338 (3.506)      |      0.038 (3.551)      |     -0.167 (3.530)      |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |     -6.156 - 8.873      |     -6.079 - 8.940      |     -6.156 - 8.940      |        |
    ## |**x_coord**                            |                         |                         |                         |   0.010|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            | 428271.250 (56003.689)  | 436143.072 (60211.610)  | 431841.530 (58062.493)  |        |
    ## |&nbsp;&nbsp;&nbsp;Range                | 321000.000 - 542000.000 | 327000.000 - 540000.000 | 321000.000 - 542000.000 |        |
    ## |**y_coord**                            |                         |                         |                         |   0.108|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            | 342893.750 (133805.914) | 331689.759 (131765.710) | 337812.158 (132956.294) |        |
    ## |&nbsp;&nbsp;&nbsp;Range                | 149000.000 - 580000.000 | 152000.000 - 580000.000 | 149000.000 - 580000.000 |        |
    ## |**smoking**                            |                         |                         |                         |   0.002|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |            0            |            2            |            2            |        |
    ## |&nbsp;&nbsp;&nbsp;Current              |       135 (16.9%)       |       68 (10.3%)        |       203 (13.9%)       |        |
    ## |&nbsp;&nbsp;&nbsp;Never                |       353 (44.1%)       |       302 (45.6%)       |       655 (44.8%)       |        |
    ## |&nbsp;&nbsp;&nbsp;Prefer not to answer |        4 (0.5%)         |        6 (0.9%)         |        10 (0.7%)        |        |
    ## |&nbsp;&nbsp;&nbsp;Previous             |       308 (38.5%)       |       286 (43.2%)       |       594 (40.6%)       |        |
    ## |**fev**                                |                         |                         |                         |   0.084|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |           559           |           484           |          1043           |        |
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |      2.899 (0.531)      |      2.992 (0.567)      |      2.939 (0.548)      |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |      1.658 - 4.522      |      1.954 - 4.098      |      1.658 - 4.522      |        |
    ## |**copd**                               |                         |                         |                         |   0.055|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |           644           |           564           |          1208           |        |
    ## |&nbsp;&nbsp;&nbsp;No                   |       147 (94.2%)       |       99 (99.0%)        |       246 (96.1%)       |        |
    ## |&nbsp;&nbsp;&nbsp;Yes                  |        9 (5.8%)         |        1 (1.0%)         |        10 (3.9%)        |        |
    ## |**WP_dusty**                           |                         |                         |                         |   0.189|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |           645           |           566           |          1211           |        |
    ## |&nbsp;&nbsp;&nbsp;Often                |        7 (4.5%)         |       10 (10.2%)        |        17 (6.7%)        |        |
    ## |&nbsp;&nbsp;&nbsp;Rarely/never         |       119 (76.8%)       |       73 (74.5%)        |       192 (75.9%)       |        |
    ## |&nbsp;&nbsp;&nbsp;Sometimes            |       29 (18.7%)        |       15 (15.3%)        |       44 (17.4%)        |        |
    ## |**WP_chemicals**                       |                         |                         |                         |   0.182|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |           646           |           566           |          1212           |        |
    ## |&nbsp;&nbsp;&nbsp;Do not know          |        1 (0.6%)         |        2 (2.0%)         |        3 (1.2%)         |        |
    ## |&nbsp;&nbsp;&nbsp;Often                |        2 (1.3%)         |        5 (5.1%)         |        7 (2.8%)         |        |
    ## |&nbsp;&nbsp;&nbsp;Rarely/never         |       136 (88.3%)       |       79 (80.6%)        |       215 (85.3%)       |        |
    ## |&nbsp;&nbsp;&nbsp;Sometimes            |        15 (9.7%)        |       12 (12.2%)        |       27 (10.7%)        |        |
    ## |**WP_cig**                             |                         |                         |                         |   0.843|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |           645           |           566           |          1211           |        |
    ## |&nbsp;&nbsp;&nbsp;Do not know          |        1 (0.6%)         |        0 (0.0%)         |        1 (0.4%)         |        |
    ## |&nbsp;&nbsp;&nbsp;Often                |        8 (5.2%)         |        6 (6.1%)         |        14 (5.5%)        |        |
    ## |&nbsp;&nbsp;&nbsp;Rarely/never         |       115 (74.2%)       |       74 (75.5%)        |       189 (74.7%)       |        |
    ## |&nbsp;&nbsp;&nbsp;Sometimes            |       31 (20.0%)        |       18 (18.4%)        |       49 (19.4%)        |        |
    ## |**WP_diesel**                          |                         |                         |                         |   0.383|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |           645           |           566           |          1211           |        |
    ## |&nbsp;&nbsp;&nbsp;Do not know          |        2 (1.3%)         |        3 (3.1%)         |        5 (2.0%)         |        |
    ## |&nbsp;&nbsp;&nbsp;Often                |        2 (1.3%)         |        4 (4.1%)         |        6 (2.4%)         |        |
    ## |&nbsp;&nbsp;&nbsp;Rarely/never         |       141 (91.0%)       |       85 (86.7%)        |       226 (89.3%)       |        |
    ## |&nbsp;&nbsp;&nbsp;Sometimes            |        10 (6.5%)        |        6 (6.1%)         |        16 (6.3%)        |        |
    ## |**breathing**                          |                         |                         |                         |   0.256|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |           645           |           566           |          1211           |        |
    ## |&nbsp;&nbsp;&nbsp;No                   |       144 (92.9%)       |       87 (88.8%)        |       231 (91.3%)       |        |
    ## |&nbsp;&nbsp;&nbsp;Yes                  |        11 (7.1%)        |       11 (11.2%)        |        22 (8.7%)        |        |
    ## |**whistling**                          |                         |                         |                         |   0.834|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |            0            |            2            |            2            |        |
    ## |&nbsp;&nbsp;&nbsp;Do not know          |        25 (3.1%)        |        16 (2.4%)        |        41 (2.8%)        |        |
    ## |&nbsp;&nbsp;&nbsp;No                   |       527 (65.9%)       |       442 (66.8%)       |       969 (66.3%)       |        |
    ## |&nbsp;&nbsp;&nbsp;Prefer not to answer |        2 (0.2%)         |        1 (0.2%)         |        3 (0.2%)         |        |
    ## |&nbsp;&nbsp;&nbsp;Yes                  |       246 (30.8%)       |       203 (30.7%)       |       449 (30.7%)       |        |
    ## |**diabetes**                           |                         |                         |                         |   0.816|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |            0            |            2            |            2            |        |
    ## |&nbsp;&nbsp;&nbsp;Do not know          |        3 (0.4%)         |        3 (0.5%)         |        6 (0.4%)         |        |
    ## |&nbsp;&nbsp;&nbsp;No                   |       721 (90.1%)       |       588 (88.8%)       |      1309 (89.5%)       |        |
    ## |&nbsp;&nbsp;&nbsp;Prefer not to answer |        2 (0.2%)         |        1 (0.2%)         |        3 (0.2%)         |        |
    ## |&nbsp;&nbsp;&nbsp;Yes                  |        74 (9.2%)        |       70 (10.6%)        |       144 (9.8%)        |        |
    ## |**sex**                                |                         |                         |                         |   0.043|
    ## |&nbsp;&nbsp;&nbsp;Female               |       393 (49.1%)       |       291 (43.8%)       |       684 (46.7%)       |        |
    ## |&nbsp;&nbsp;&nbsp;Male                 |       407 (50.9%)       |       373 (56.2%)       |       780 (53.3%)       |        |
    ## |**birthYear**                          |                         |                         |                         |   0.396|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |    1950.341 (8.665)     |    1950.729 (8.737)     |    1950.517 (8.697)     |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |   1937.000 - 1970.000   |   1937.000 - 1969.000   |   1937.000 - 1970.000   |        |
    ## |**diaBP**                              |                         |                         |                         |   0.045|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |           16            |           21            |           37            |        |
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |     81.922 (10.923)     |     83.075 (10.603)     |     82.441 (10.791)     |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |    46.000 - 125.000     |    57.000 - 120.000     |    46.000 - 125.000     |        |
    ## |**sysBP**                              |                         |                         |                         |   0.361|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |           17            |           21            |           38            |        |
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |    136.498 (19.270)     |    137.440 (19.499)     |    136.923 (19.372)     |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |    92.000 - 219.000     |    91.000 - 229.000     |    91.000 - 229.000     |        |
    ## |**waist**                              |                         |                         |                         |   0.010|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |            7            |            5            |           12            |        |
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |     94.368 (15.035)     |     96.393 (14.517)     |     95.287 (14.831)     |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |    61.000 - 151.000     |    62.000 - 166.000     |    61.000 - 166.000     |        |
    ## |**hip**                                |                         |                         |                         |   0.142|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |            7            |            5            |           12            |        |
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |    104.880 (11.008)     |    105.732 (11.006)     |    105.267 (11.012)     |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |    83.000 - 158.000     |    82.000 - 172.000     |    82.000 - 172.000     |        |
    ## |**height**                             |                         |                         |                         |   0.995|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |            8            |            7            |           15            |        |
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |     168.848 (9.503)     |     168.851 (9.247)     |     168.849 (9.385)     |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |    143.000 - 197.000    |    141.000 - 196.000    |    141.000 - 197.000    |        |
    ## |**age**                                |                         |                         |                         |   0.396|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |     69.659 (8.665)      |     69.271 (8.737)      |     69.483 (8.697)      |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |     50.000 - 83.000     |     51.000 - 83.000     |     50.000 - 83.000     |        |
    ## |**whr**                                |                         |                         |                         |   0.010|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |            7            |            5            |           12            |        |
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |      0.898 (0.092)      |      0.910 (0.090)      |      0.904 (0.091)      |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |      0.642 - 1.211      |      0.660 - 1.194      |      0.642 - 1.211      |        |
    ## |**highBP**                             |                         |                         |                         |   0.199|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |           16            |           21            |           37            |        |
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |      0.189 (0.392)      |      0.216 (0.412)      |      0.201 (0.401)      |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |      0.000 - 1.000      |      0.000 - 1.000      |      0.000 - 1.000      |        |
    ## |**inpatient_covid**                    |                         |                         |                         | < 0.001|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |      0.672 (0.470)      |      0.858 (0.349)      |      0.757 (0.429)      |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |      0.000 - 1.000      |      0.000 - 1.000      |      0.000 - 1.000      |        |
    ## |**distance**                           |                         |                         |                         |   0.070|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |    963.123 (21.328)     |    961.081 (21.525)     |    962.197 (21.434)     |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |   930.406 - 1001.337    |   930.476 - 1001.203    |   930.406 - 1001.337    |        |
    ## |**pm25_val**                           |                         |                         |                         | < 0.001|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |      8.859 (1.853)      |      9.218 (1.924)      |      9.022 (1.893)      |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |     5.630 - 13.650      |     5.769 - 13.667      |     5.630 - 13.667      |        |
    ## |**merged_popDens_km2**                 |                         |                         |                         |   0.008|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |   2770.719 (2909.404)   |   3197.601 (3228.642)   |   2964.332 (3064.644)   |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |   64.000 - 16097.000    |   64.000 - 16097.000    |   64.000 - 16097.000    |        |
    ## |**smoker**                             |                         |                         |                         | < 0.001|
    ## |&nbsp;&nbsp;&nbsp;N-Miss               |            0            |            2            |            2            |        |
    ## |&nbsp;&nbsp;&nbsp;not_smoking          |       665 (83.1%)       |       594 (89.7%)       |      1259 (86.1%)       |        |
    ## |&nbsp;&nbsp;&nbsp;smoking              |       135 (16.9%)       |       68 (10.3%)        |       203 (13.9%)       |        |
    ## |**no2_val**                            |                         |                         |                         | < 0.001|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |     15.485 (6.074)      |     16.789 (6.469)      |     16.076 (6.288)      |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |     4.987 - 48.672      |     5.632 - 48.672      |     4.987 - 48.672      |        |
    ## |**so2_val**                            |                         |                         |                         |   0.027|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |      1.668 (0.546)      |      1.731 (0.544)      |      1.697 (0.545)      |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |      0.661 - 5.108      |      0.700 - 4.578      |      0.661 - 5.108      |        |
    ## |**o3_val**                             |                         |                         |                         |   0.753|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |      5.621 (2.653)      |      5.576 (2.668)      |      5.601 (2.659)      |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |     0.000 - 13.188      |     0.000 - 12.731      |     0.000 - 13.188      |        |
    ## |**nox_val**                            |                         |                         |                         | < 0.001|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |     21.955 (10.770)     |     24.163 (11.841)     |     22.956 (11.318)     |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |     6.254 - 95.880      |     7.116 - 95.880      |     6.254 - 95.880      |        |
    ## |**pm10_val**                           |                         |                         |                         | < 0.001|
    ## |&nbsp;&nbsp;&nbsp;Mean (SD)            |     13.430 (2.683)      |     13.925 (2.816)      |     13.655 (2.755)      |        |
    ## |&nbsp;&nbsp;&nbsp;Range                |     8.244 - 21.178      |     8.735 - 21.207      |     8.244 - 21.207      |        |

``` r
write2pdf(ukb_descript_stats, "ukb_descript_stats_4_2020.pdf")
```

    ##   |                                                                              |                                                                      |   0%  |                                                                              |......................................................................| 100%
    ##   ordinary text without R code
    ## 
    ## 
    ## /usr/local/bin/pandoc +RTS -K512m -RTS ukb_descript_stats_4_2020.pdf.utf8.md --to latex --from markdown+autolink_bare_uris+tex_math_single_backslash --output ukb_descript_stats_4_2020.tex --self-contained --highlight-style tango --pdf-engine pdflatex --variable graphics --lua-filter /Users/yizhouyu/Library/R/4.0/library/rmarkdown/rmd/lua/pagebreak.lua --lua-filter /Users/yizhouyu/Library/R/4.0/library/rmarkdown/rmd/lua/latex-div.lua --variable 'geometry:margin=1in'

### model with all pollutants: glm

``` r
# PM2.5
summary(ukb_covid_pm25.nb <-glm.nb(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                     whr + sex + age + pm25_val))
```

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached
    
    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## 
    ## Call:
    ## glm.nb(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + pm25_val, 
    ##     data = ukb_covid_allPol_df, init.theta = 15229.88441, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.2059  -0.9479  -0.8028   0.6589   1.1708  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                   -2.065e+00  8.543e-01  -2.418  0.01562 * 
    ## merged_popDens_km2            -9.956e-06  1.884e-05  -0.529  0.59711   
    ## n_cancers                     -1.258e-01  1.272e-01  -0.989  0.32253   
    ## townsend                       1.354e-02  1.310e-02   1.034  0.30121   
    ## smokingNever                   3.665e-01  1.379e-01   2.658  0.00785 **
    ## smokingPrefer not to answer    7.608e-01  4.478e-01   1.699  0.08929 . 
    ## smokingPrevious                4.056e-01  1.377e-01   2.946  0.00322 **
    ## whistlingNo                    2.671e-01  2.698e-01   0.990  0.32220   
    ## whistlingPrefer not to answer  4.568e-02  1.276e+00   0.036  0.97145   
    ## whistlingYes                   2.533e-01  2.742e-01   0.924  0.35559   
    ## diabetesNo                    -1.619e-01  6.106e-01  -0.265  0.79085   
    ## diabetesPrefer not to answer  -6.877e-01  1.433e+00  -0.480  0.63140   
    ## diabetesYes                   -1.956e-01  6.195e-01  -0.316  0.75215   
    ## whr                            7.856e-01  5.839e-01   1.346  0.17846   
    ## sexMale                        3.340e-02  1.000e-01   0.334  0.73838   
    ## age                           -5.681e-03  4.748e-03  -1.196  0.23153   
    ## pm25_val                       6.035e-02  2.866e-02   2.106  0.03520 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(15229.88) family taken to be 1)
    ## 
    ##     Null deviance: 1040.2  on 1449  degrees of freedom
    ## Residual deviance: 1015.6  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 2365.6
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  15230 
    ##           Std. Err.:  58208 
    ## Warning while fitting theta: iteration limit reached 
    ## 
    ##  2 x log-likelihood:  -2329.646

``` r
summary(ukb_covid_pm25.p <-glm(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                 whr + sex + age + pm25_val, family = 'poisson'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + pm25_val, 
    ##     family = "poisson", data = ukb_covid_allPol_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.2060  -0.9479  -0.8028   0.6589   1.1708  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                   -2.065e+00  8.543e-01  -2.418  0.01562 * 
    ## merged_popDens_km2            -9.956e-06  1.883e-05  -0.529  0.59711   
    ## n_cancers                     -1.258e-01  1.272e-01  -0.989  0.32252   
    ## townsend                       1.354e-02  1.310e-02   1.034  0.30120   
    ## smokingNever                   3.665e-01  1.379e-01   2.658  0.00785 **
    ## smokingPrefer not to answer    7.609e-01  4.478e-01   1.699  0.08928 . 
    ## smokingPrevious                4.056e-01  1.377e-01   2.946  0.00321 **
    ## whistlingNo                    2.671e-01  2.698e-01   0.990  0.32219   
    ## whistlingPrefer not to answer  4.567e-02  1.276e+00   0.036  0.97145   
    ## whistlingYes                   2.533e-01  2.742e-01   0.924  0.35558   
    ## diabetesNo                    -1.619e-01  6.105e-01  -0.265  0.79084   
    ## diabetesPrefer not to answer  -6.877e-01  1.433e+00  -0.480  0.63138   
    ## diabetesYes                   -1.956e-01  6.194e-01  -0.316  0.75214   
    ## whr                            7.856e-01  5.839e-01   1.346  0.17845   
    ## sexMale                        3.340e-02  1.000e-01   0.334  0.73837   
    ## age                           -5.681e-03  4.748e-03  -1.196  0.23151   
    ## pm25_val                       6.035e-02  2.866e-02   2.106  0.03520 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 1040.2  on 1449  degrees of freedom
    ## Residual deviance: 1015.6  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 2363.6
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
car::vif(ukb_covid_pm25.nb)
```

    ##                        GVIF Df GVIF^(1/(2*Df))
    ## merged_popDens_km2 2.443689  1        1.563230
    ## n_cancers          1.019337  1        1.009622
    ## townsend           1.442148  1        1.200895
    ## smoking            1.259363  3        1.039182
    ## whistling          1.875826  3        1.110535
    ## diabetes           2.049216  3        1.127019
    ## whr                1.850854  1        1.360461
    ## sex                1.619879  1        1.272744
    ## age                1.138711  1        1.067104
    ## pm25_val           2.011763  1        1.418366

``` r
# according to http://www.jstor.org/stable/2290467, GVIF^(1/(2*Df)) < 5 is ok 
# model reduction: 

summary(ukb_covid_pm25.nb_red <-glm.nb(data = ukb_covid_pm25_popDens_df, result ~ merged_popDens_km2 + smoking +pm25_val))
```

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached
    
    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## 
    ## Call:
    ## glm.nb(formula = result ~ merged_popDens_km2 + smoking + pm25_val, 
    ##     data = ukb_covid_pm25_popDens_df, init.theta = 15585.93096, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.1217  -0.9499  -0.8075   0.6755   1.0620  
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                 -1.639e+00  2.629e-01  -6.235 4.52e-10 ***
    ## merged_popDens_km2          -9.261e-07  1.680e-05  -0.055  0.95604    
    ## smokingNever                 3.369e-01  1.344e-01   2.507  0.01216 *  
    ## smokingPrefer not to answer  6.113e-01  4.261e-01   1.435  0.15139    
    ## smokingPrevious              3.899e-01  1.353e-01   2.882  0.00395 ** 
    ## pm25_val                     5.801e-02  2.833e-02   2.048  0.04055 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(15585.93) family taken to be 1)
    ## 
    ##     Null deviance: 1049.0  on 1461  degrees of freedom
    ## Residual deviance: 1032.7  on 1456  degrees of freedom
    ##   (2 observations deleted due to missingness)
    ## AIC: 2370.7
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  15586 
    ##           Std. Err.:  60471 
    ## Warning while fitting theta: iteration limit reached 
    ## 
    ##  2 x log-likelihood:  -2356.739

``` r
summary(ukb_covid_pm25.p_red <-glm(data = ukb_covid_pm25_popDens_df, result ~ merged_popDens_km2 + smoking + pm25_val, family = 'poisson'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + smoking + pm25_val, 
    ##     family = "poisson", data = ukb_covid_pm25_popDens_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.1218  -0.9499  -0.8075   0.6756   1.0620  
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                 -1.639e+00  2.629e-01  -6.235 4.52e-10 ***
    ## merged_popDens_km2          -9.261e-07  1.680e-05  -0.055  0.95604    
    ## smokingNever                 3.369e-01  1.344e-01   2.508  0.01216 *  
    ## smokingPrefer not to answer  6.113e-01  4.261e-01   1.435  0.15138    
    ## smokingPrevious              3.899e-01  1.353e-01   2.882  0.00395 ** 
    ## pm25_val                     5.801e-02  2.832e-02   2.048  0.04055 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 1049.0  on 1461  degrees of freedom
    ## Residual deviance: 1032.7  on 1456  degrees of freedom
    ##   (2 observations deleted due to missingness)
    ## AIC: 2368.7
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
#check 
anova(ukb_covid_pm25.nb, ukb_covid_pm25.nb_red) # p < 0.05 -> not significantly different 
```

    ## Likelihood ratio tests of Negative Binomial Models
    ## 
    ## Response: result
    ##                                                                                                     Model
    ## 1                                                                 merged_popDens_km2 + smoking + pm25_val
    ## 2 merged_popDens_km2 + n_cancers + townsend + smoking + whistling + diabetes + whr + sex + age + pm25_val
    ##      theta Resid. df    2 x log-lik.   Test    df LR stat.   Pr(Chi)
    ## 1 15585.93      1456       -2356.739                                
    ## 2 15229.88      1433       -2329.646 1 vs 2    23  27.0924 0.2520556

``` r
#PM10 
summary(ukb_covid_pm10.nb <-glm.nb(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                     whr + sex + age + pm10_val))
```

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached
    
    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## 
    ## Call:
    ## glm.nb(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + pm10_val, 
    ##     data = ukb_covid_allPol_df, init.theta = 15278.17332, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.2226  -0.9484  -0.8030   0.6608   1.1713  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                   -2.035e+00  8.581e-01  -2.371  0.01774 * 
    ## merged_popDens_km2            -9.198e-06  1.962e-05  -0.469  0.63913   
    ## n_cancers                     -1.255e-01  1.272e-01  -0.987  0.32375   
    ## townsend                       1.338e-02  1.310e-02   1.021  0.30706   
    ## smokingNever                   3.635e-01  1.378e-01   2.637  0.00836 **
    ## smokingPrefer not to answer    7.550e-01  4.479e-01   1.685  0.09189 . 
    ## smokingPrevious                4.008e-01  1.376e-01   2.913  0.00358 **
    ## whistlingNo                    2.608e-01  2.697e-01   0.967  0.33356   
    ## whistlingPrefer not to answer  4.399e-02  1.272e+00   0.035  0.97242   
    ## whistlingYes                   2.481e-01  2.741e-01   0.905  0.36539   
    ## diabetesNo                    -1.706e-01  6.104e-01  -0.279  0.77992   
    ## diabetesPrefer not to answer  -6.936e-01  1.430e+00  -0.485  0.62768   
    ## diabetesYes                   -2.053e-01  6.193e-01  -0.332  0.74021   
    ## whr                            7.986e-01  5.841e-01   1.367  0.17153   
    ## sexMale                        3.228e-02  1.000e-01   0.323  0.74690   
    ## age                           -5.740e-03  4.748e-03  -1.209  0.22678   
    ## pm10_val                       3.829e-02  2.056e-02   1.862  0.06254 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(15278.17) family taken to be 1)
    ## 
    ##     Null deviance: 1040.2  on 1449  degrees of freedom
    ## Residual deviance: 1016.5  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 2366.6
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  15278 
    ##           Std. Err.:  58535 
    ## Warning while fitting theta: iteration limit reached 
    ## 
    ##  2 x log-likelihood:  -2330.591

``` r
summary(ukb_covid_pm10.p <-glm(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                 whr + sex + age + pm10_val, family = 'poisson'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + pm10_val, 
    ##     family = "poisson", data = ukb_covid_allPol_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.2226  -0.9484  -0.8030   0.6608   1.1713  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                   -2.035e+00  8.581e-01  -2.371  0.01774 * 
    ## merged_popDens_km2            -9.198e-06  1.962e-05  -0.469  0.63912   
    ## n_cancers                     -1.255e-01  1.272e-01  -0.987  0.32374   
    ## townsend                       1.338e-02  1.310e-02   1.021  0.30705   
    ## smokingNever                   3.635e-01  1.378e-01   2.637  0.00836 **
    ## smokingPrefer not to answer    7.550e-01  4.479e-01   1.686  0.09188 . 
    ## smokingPrevious                4.008e-01  1.376e-01   2.913  0.00358 **
    ## whistlingNo                    2.608e-01  2.697e-01   0.967  0.33355   
    ## whistlingPrefer not to answer  4.399e-02  1.272e+00   0.035  0.97242   
    ## whistlingYes                   2.481e-01  2.741e-01   0.905  0.36538   
    ## diabetesNo                    -1.706e-01  6.104e-01  -0.279  0.77992   
    ## diabetesPrefer not to answer  -6.936e-01  1.430e+00  -0.485  0.62766   
    ## diabetesYes                   -2.053e-01  6.193e-01  -0.332  0.74020   
    ## whr                            7.986e-01  5.841e-01   1.367  0.17152   
    ## sexMale                        3.228e-02  1.000e-01   0.323  0.74690   
    ## age                           -5.740e-03  4.748e-03  -1.209  0.22676   
    ## pm10_val                       3.829e-02  2.056e-02   1.863  0.06253 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 1040.2  on 1449  degrees of freedom
    ## Residual deviance: 1016.6  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 2364.6
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
car::vif(ukb_covid_pm10.nb)
```

    ##                        GVIF Df GVIF^(1/(2*Df))
    ## merged_popDens_km2 2.690058  1        1.640140
    ## n_cancers          1.019347  1        1.009627
    ## townsend           1.448761  1        1.203645
    ## smoking            1.258582  3        1.039075
    ## whistling          1.863927  3        1.109357
    ## diabetes           2.037411  3        1.125934
    ## whr                1.852457  1        1.361050
    ## sex                1.620823  1        1.273115
    ## age                1.139335  1        1.067396
    ## pm10_val           2.234526  1        1.494833

``` r
#nox
summary(ukb_covid_nox.nb <-glm.nb(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                     whr + sex + age + nox_val))
```

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached
    
    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## 
    ## Call:
    ## glm.nb(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + nox_val, 
    ##     data = ukb_covid_allPol_df, init.theta = 15323.76936, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.3845  -0.9475  -0.8039   0.6603   1.1833  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                   -1.764e+00  8.293e-01  -2.127  0.03339 * 
    ## merged_popDens_km2            -1.289e-05  2.130e-05  -0.605  0.54515   
    ## n_cancers                     -1.220e-01  1.271e-01  -0.960  0.33727   
    ## townsend                       6.366e-03  1.337e-02   0.476  0.63384   
    ## smokingNever                   3.541e-01  1.377e-01   2.571  0.01015 * 
    ## smokingPrefer not to answer    7.226e-01  4.491e-01   1.609  0.10759   
    ## smokingPrevious                3.917e-01  1.374e-01   2.851  0.00436 **
    ## whistlingNo                    2.591e-01  2.697e-01   0.961  0.33670   
    ## whistlingPrefer not to answer  2.364e-03  1.245e+00   0.002  0.99849   
    ## whistlingYes                   2.438e-01  2.741e-01   0.890  0.37362   
    ## diabetesNo                    -1.971e-01  6.103e-01  -0.323  0.74668   
    ## diabetesPrefer not to answer  -7.274e-01  1.408e+00  -0.516  0.60552   
    ## diabetesYes                   -2.279e-01  6.191e-01  -0.368  0.71277   
    ## whr                            8.545e-01  5.847e-01   1.462  0.14387   
    ## sexMale                        1.812e-02  1.004e-01   0.180  0.85680   
    ## age                           -5.503e-03  4.747e-03  -1.159  0.24641   
    ## nox_val                        1.043e-02  5.565e-03   1.874  0.06097 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(15323.77) family taken to be 1)
    ## 
    ##     Null deviance: 1040.2  on 1449  degrees of freedom
    ## Residual deviance: 1016.6  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 2366.6
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  15324 
    ##           Std. Err.:  58831 
    ## Warning while fitting theta: iteration limit reached 
    ## 
    ##  2 x log-likelihood:  -2330.647

``` r
summary(ukb_covid_nox.p <-glm(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                 whr + sex + age + nox_val, family = 'poisson'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + nox_val, 
    ##     family = "poisson", data = ukb_covid_allPol_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.3845  -0.9476  -0.8039   0.6603   1.1833  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                   -1.764e+00  8.293e-01  -2.127  0.03339 * 
    ## merged_popDens_km2            -1.288e-05  2.130e-05  -0.605  0.54515   
    ## n_cancers                     -1.220e-01  1.271e-01  -0.960  0.33726   
    ## townsend                       6.366e-03  1.337e-02   0.476  0.63383   
    ## smokingNever                   3.541e-01  1.377e-01   2.571  0.01015 * 
    ## smokingPrefer not to answer    7.226e-01  4.491e-01   1.609  0.10758   
    ## smokingPrevious                3.917e-01  1.374e-01   2.851  0.00436 **
    ## whistlingNo                    2.591e-01  2.697e-01   0.961  0.33668   
    ## whistlingPrefer not to answer  2.361e-03  1.245e+00   0.002  0.99849   
    ## whistlingYes                   2.438e-01  2.741e-01   0.890  0.37360   
    ## diabetesNo                    -1.971e-01  6.103e-01  -0.323  0.74667   
    ## diabetesPrefer not to answer  -7.274e-01  1.408e+00  -0.517  0.60549   
    ## diabetesYes                   -2.279e-01  6.190e-01  -0.368  0.71277   
    ## whr                            8.545e-01  5.847e-01   1.462  0.14387   
    ## sexMale                        1.812e-02  1.004e-01   0.180  0.85679   
    ## age                           -5.503e-03  4.747e-03  -1.159  0.24640   
    ## nox_val                        1.043e-02  5.565e-03   1.874  0.06097 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 1040.2  on 1449  degrees of freedom
    ## Residual deviance: 1016.6  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 2364.6
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
car::vif(ukb_covid_nox.nb)
```

    ##                        GVIF Df GVIF^(1/(2*Df))
    ## merged_popDens_km2 3.190492  1        1.786195
    ## n_cancers          1.020066  1        1.009983
    ## townsend           1.522514  1        1.233902
    ## smoking            1.265027  3        1.039960
    ## whistling          1.781014  3        1.100976
    ## diabetes           1.961134  3        1.118797
    ## whr                1.848865  1        1.359730
    ## sex                1.633899  1        1.278241
    ## age                1.141394  1        1.068361
    ## nox_val            3.273151  1        1.809185

``` r
#no2
summary(ukb_covid_no2.nb <-glm.nb(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                     whr + sex + age + no2_val))
```

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached
    
    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## 
    ## Call:
    ## glm.nb(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + no2_val, 
    ##     data = ukb_covid_allPol_df, init.theta = 15256.03065, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.3433  -0.9455  -0.8015   0.6594   1.1994  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                   -1.867e+00  8.335e-01  -2.241  0.02506 * 
    ## merged_popDens_km2            -1.812e-05  2.128e-05  -0.852  0.39434   
    ## n_cancers                     -1.211e-01  1.271e-01  -0.953  0.34075   
    ## townsend                       5.084e-03  1.340e-02   0.379  0.70433   
    ## smokingNever                   3.554e-01  1.377e-01   2.580  0.00987 **
    ## smokingPrefer not to answer    7.209e-01  4.493e-01   1.605  0.10860   
    ## smokingPrevious                3.954e-01  1.375e-01   2.876  0.00402 **
    ## whistlingNo                    2.585e-01  2.700e-01   0.957  0.33843   
    ## whistlingPrefer not to answer -1.309e-02  1.238e+00  -0.011  0.99156   
    ## whistlingYes                   2.449e-01  2.744e-01   0.892  0.37217   
    ## diabetesNo                    -1.911e-01  6.113e-01  -0.313  0.75451   
    ## diabetesPrefer not to answer  -7.362e-01  1.403e+00  -0.525  0.59967   
    ## diabetesYes                   -2.218e-01  6.200e-01  -0.358  0.72050   
    ## whr                            8.312e-01  5.848e-01   1.421  0.15525   
    ## sexMale                        1.920e-02  1.003e-01   0.191  0.84815   
    ## age                           -5.471e-03  4.748e-03  -1.152  0.24919   
    ## no2_val                        2.283e-02  1.039e-02   2.197  0.02804 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(15256.03) family taken to be 1)
    ## 
    ##     Null deviance: 1040.2  on 1449  degrees of freedom
    ## Residual deviance: 1015.3  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 2365.3
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  15256 
    ##           Std. Err.:  58376 
    ## Warning while fitting theta: iteration limit reached 
    ## 
    ##  2 x log-likelihood:  -2329.3

``` r
summary(ukb_covid_no2.p <-glm(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                 whr + sex + age + no2_val, family = 'poisson'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + no2_val, 
    ##     family = "poisson", data = ukb_covid_allPol_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.3434  -0.9455  -0.8015   0.6594   1.1994  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                   -1.867e+00  8.335e-01  -2.241  0.02505 * 
    ## merged_popDens_km2            -1.812e-05  2.127e-05  -0.852  0.39433   
    ## n_cancers                     -1.211e-01  1.271e-01  -0.953  0.34074   
    ## townsend                       5.084e-03  1.340e-02   0.379  0.70432   
    ## smokingNever                   3.554e-01  1.377e-01   2.580  0.00987 **
    ## smokingPrefer not to answer    7.209e-01  4.493e-01   1.605  0.10859   
    ## smokingPrevious                3.954e-01  1.375e-01   2.877  0.00402 **
    ## whistlingNo                    2.585e-01  2.700e-01   0.957  0.33841   
    ## whistlingPrefer not to answer -1.310e-02  1.238e+00  -0.011  0.99156   
    ## whistlingYes                   2.449e-01  2.744e-01   0.892  0.37215   
    ## diabetesNo                    -1.911e-01  6.113e-01  -0.313  0.75450   
    ## diabetesPrefer not to answer  -7.362e-01  1.402e+00  -0.525  0.59963   
    ## diabetesYes                   -2.218e-01  6.200e-01  -0.358  0.72049   
    ## whr                            8.312e-01  5.848e-01   1.421  0.15524   
    ## sexMale                        1.920e-02  1.003e-01   0.191  0.84815   
    ## age                           -5.471e-03  4.747e-03  -1.152  0.24918   
    ## no2_val                        2.283e-02  1.039e-02   2.197  0.02804 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 1040.2  on 1449  degrees of freedom
    ## Residual deviance: 1015.3  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 2363.3
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
car::vif(ukb_covid_no2.nb)
```

    ##                        GVIF Df GVIF^(1/(2*Df))
    ## merged_popDens_km2 3.144160  1        1.773178
    ## n_cancers          1.020039  1        1.009970
    ## townsend           1.523933  1        1.234477
    ## smoking            1.267137  3        1.040249
    ## whistling          1.763147  3        1.099128
    ## diabetes           1.944539  3        1.117213
    ## whr                1.847941  1        1.359390
    ## sex                1.628714  1        1.276211
    ## age                1.141009  1        1.068180
    ## no2_val            3.255183  1        1.804212

``` r
#o3
summary(ukb_covid_o3.nb <-glm.nb(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                    whr + sex + age + o3_val))
```

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached
    
    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## 
    ## Call:
    ## glm.nb(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + o3_val, 
    ##     data = ukb_covid_allPol_df, init.theta = 15368.98351, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.2076  -0.9503  -0.8062   0.6636   1.1391  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                   -1.630e+00  8.293e-01  -1.966  0.04931 * 
    ## merged_popDens_km2             1.614e-05  1.399e-05   1.154  0.24867   
    ## n_cancers                     -1.289e-01  1.274e-01  -1.012  0.31151   
    ## townsend                       1.333e-02  1.346e-02   0.990  0.32206   
    ## smokingNever                   3.585e-01  1.381e-01   2.596  0.00942 **
    ## smokingPrefer not to answer    7.494e-01  4.493e-01   1.668  0.09528 . 
    ## smokingPrevious                3.905e-01  1.377e-01   2.837  0.00455 **
    ## whistlingNo                    2.435e-01  2.687e-01   0.906  0.36477   
    ## whistlingPrefer not to answer  2.093e-02  1.278e+00   0.016  0.98693   
    ## whistlingYes                   2.276e-01  2.730e-01   0.834  0.40454   
    ## diabetesNo                    -1.858e-01  6.083e-01  -0.305  0.76009   
    ## diabetesPrefer not to answer  -6.932e-01  1.435e+00  -0.483  0.62909   
    ## diabetesYes                   -2.278e-01  6.172e-01  -0.369  0.71205   
    ## whr                            8.310e-01  5.842e-01   1.422  0.15489   
    ## sexMale                        3.618e-02  1.002e-01   0.361  0.71815   
    ## age                           -5.563e-03  4.749e-03  -1.171  0.24143   
    ## o3_val                         7.573e-03  1.544e-02   0.491  0.62369   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(15368.98) family taken to be 1)
    ## 
    ##     Null deviance: 1040.2  on 1449  degrees of freedom
    ## Residual deviance: 1019.8  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 2369.8
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  15369 
    ##           Std. Err.:  59189 
    ## Warning while fitting theta: iteration limit reached 
    ## 
    ##  2 x log-likelihood:  -2333.797

``` r
summary(ukb_covid_o3.p <-glm(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                whr + sex + age + o3_val, family = 'poisson'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + o3_val, 
    ##     family = "poisson", data = ukb_covid_allPol_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.2076  -0.9503  -0.8062   0.6636   1.1391  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                   -1.630e+00  8.292e-01  -1.966  0.04931 * 
    ## merged_popDens_km2             1.614e-05  1.399e-05   1.154  0.24866   
    ## n_cancers                     -1.289e-01  1.274e-01  -1.012  0.31150   
    ## townsend                       1.333e-02  1.346e-02   0.990  0.32205   
    ## smokingNever                   3.585e-01  1.381e-01   2.596  0.00942 **
    ## smokingPrefer not to answer    7.494e-01  4.492e-01   1.668  0.09527 . 
    ## smokingPrevious                3.905e-01  1.377e-01   2.837  0.00455 **
    ## whistlingNo                    2.435e-01  2.687e-01   0.906  0.36476   
    ## whistlingPrefer not to answer  2.092e-02  1.278e+00   0.016  0.98693   
    ## whistlingYes                   2.276e-01  2.730e-01   0.834  0.40453   
    ## diabetesNo                    -1.858e-01  6.083e-01  -0.305  0.76008   
    ## diabetesPrefer not to answer  -6.932e-01  1.435e+00  -0.483  0.62906   
    ## diabetesYes                   -2.278e-01  6.172e-01  -0.369  0.71205   
    ## whr                            8.310e-01  5.842e-01   1.423  0.15488   
    ## sexMale                        3.618e-02  1.002e-01   0.361  0.71814   
    ## age                           -5.563e-03  4.749e-03  -1.171  0.24142   
    ## o3_val                         7.573e-03  1.544e-02   0.491  0.62368   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 1040.2  on 1449  degrees of freedom
    ## Residual deviance: 1019.8  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 2367.8
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
car::vif(ukb_covid_o3.nb)
```

    ##                        GVIF Df GVIF^(1/(2*Df))
    ## merged_popDens_km2 1.370137  1        1.170528
    ## n_cancers          1.019831  1        1.009867
    ## townsend           1.535336  1        1.239087
    ## smoking            1.267302  3        1.040271
    ## whistling          1.867480  3        1.109710
    ## diabetes           2.047511  3        1.126863
    ## whr                1.856088  1        1.362383
    ## sex                1.627746  1        1.275831
    ## age                1.140751  1        1.068060
    ## o3_val             1.098307  1        1.048001

``` r
#so2
summary(ukb_covid_so2.nb <-glm.nb(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                   whr + sex + age + so2_val))
```

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached
    
    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## 
    ## Call:
    ## glm.nb(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + so2_val, 
    ##     data = ukb_covid_allPol_df, init.theta = 15326.12162, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.2325  -0.9485  -0.8076   0.6604   1.2053  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                   -1.721e+00  8.342e-01  -2.063  0.03915 * 
    ## merged_popDens_km2             1.462e-05  1.409e-05   1.037  0.29962   
    ## n_cancers                     -1.272e-01  1.274e-01  -0.998  0.31824   
    ## townsend                       8.033e-03  1.345e-02   0.597  0.55043   
    ## smokingNever                   3.566e-01  1.378e-01   2.587  0.00967 **
    ## smokingPrefer not to answer    7.159e-01  4.487e-01   1.596  0.11060   
    ## smokingPrevious                3.901e-01  1.375e-01   2.837  0.00456 **
    ## whistlingNo                    2.472e-01  2.689e-01   0.920  0.35778   
    ## whistlingPrefer not to answer -7.566e-03  1.260e+00  -0.006  0.99521   
    ## whistlingYes                   2.310e-01  2.732e-01   0.846  0.39781   
    ## diabetesNo                    -1.926e-01  6.089e-01  -0.316  0.75176   
    ## diabetesPrefer not to answer  -6.932e-01  1.420e+00  -0.488  0.62534   
    ## diabetesYes                   -2.308e-01  6.178e-01  -0.374  0.70864   
    ## whr                            8.300e-01  5.845e-01   1.420  0.15560   
    ## sexMale                        2.495e-02  1.004e-01   0.248  0.80378   
    ## age                           -5.453e-03  4.748e-03  -1.148  0.25078   
    ## so2_val                        8.195e-02  7.699e-02   1.065  0.28709   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(15326.12) family taken to be 1)
    ## 
    ##     Null deviance: 1040.2  on 1449  degrees of freedom
    ## Residual deviance: 1018.9  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 2368.9
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  15326 
    ##           Std. Err.:  58891 
    ## Warning while fitting theta: iteration limit reached 
    ## 
    ##  2 x log-likelihood:  -2332.922

``` r
summary(ukb_covid_so2.p <-glm(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                               whr + sex + age + so2_val, family = 'poisson'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + so2_val, 
    ##     family = "poisson", data = ukb_covid_allPol_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.2325  -0.9485  -0.8076   0.6605   1.2054  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                   -1.721e+00  8.342e-01  -2.063  0.03915 * 
    ## merged_popDens_km2             1.462e-05  1.409e-05   1.037  0.29961   
    ## n_cancers                     -1.272e-01  1.274e-01  -0.998  0.31823   
    ## townsend                       8.033e-03  1.345e-02   0.597  0.55042   
    ## smokingNever                   3.566e-01  1.378e-01   2.587  0.00967 **
    ## smokingPrefer not to answer    7.159e-01  4.487e-01   1.596  0.11059   
    ## smokingPrevious                3.901e-01  1.375e-01   2.837  0.00456 **
    ## whistlingNo                    2.472e-01  2.689e-01   0.920  0.35776   
    ## whistlingPrefer not to answer -7.570e-03  1.260e+00  -0.006  0.99521   
    ## whistlingYes                   2.310e-01  2.732e-01   0.846  0.39779   
    ## diabetesNo                    -1.926e-01  6.089e-01  -0.316  0.75175   
    ## diabetesPrefer not to answer  -6.933e-01  1.420e+00  -0.488  0.62532   
    ## diabetesYes                   -2.309e-01  6.178e-01  -0.374  0.70863   
    ## whr                            8.300e-01  5.845e-01   1.420  0.15559   
    ## sexMale                        2.495e-02  1.004e-01   0.248  0.80377   
    ## age                           -5.453e-03  4.748e-03  -1.148  0.25077   
    ## so2_val                        8.195e-02  7.698e-02   1.065  0.28708   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 1040.2  on 1449  degrees of freedom
    ## Residual deviance: 1018.9  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 2366.9
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
car::vif(ukb_covid_so2.nb)
```

    ##                        GVIF Df GVIF^(1/(2*Df))
    ## merged_popDens_km2 1.386609  1        1.177543
    ## n_cancers          1.019429  1        1.009668
    ## townsend           1.534094  1        1.238585
    ## smoking            1.264270  3        1.039856
    ## whistling          1.816373  3        1.104589
    ## diabetes           1.991278  3        1.121645
    ## whr                1.852603  1        1.361103
    ## sex                1.633541  1        1.278100
    ## age                1.141044  1        1.068196
    ## so2_val            1.185954  1        1.089015

### model with all pollutants: binary

binary model to account for the fact that the response variable is 1/0

``` r
# PM2.5
summary(ukb_covid_pm25.b <-glm(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                 whr + sex + age + pm25_val, family = 'binomial'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + pm25_val, 
    ##     family = "binomial", data = ukb_covid_allPol_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.5015  -1.0980  -0.8723   1.2076   1.7802  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                   -2.495e+00  1.201e+00  -2.078 0.037729 *  
    ## merged_popDens_km2            -1.665e-05  2.673e-05  -0.623 0.533410    
    ## n_cancers                     -2.279e-01  1.656e-01  -1.376 0.168882    
    ## townsend                       2.553e-02  1.803e-02   1.416 0.156877    
    ## smokingNever                   6.317e-01  1.760e-01   3.589 0.000332 ***
    ## smokingPrefer not to answer    1.517e+00  7.460e-01   2.033 0.042046 *  
    ## smokingPrevious                7.049e-01  1.764e-01   3.996 6.45e-05 ***
    ## whistlingNo                    4.731e-01  3.493e-01   1.355 0.175573    
    ## whistlingPrefer not to answer  1.233e-01  1.636e+00   0.075 0.939912    
    ## whistlingYes                   4.484e-01  3.556e-01   1.261 0.207277    
    ## diabetesNo                    -2.793e-01  8.588e-01  -0.325 0.744986    
    ## diabetesPrefer not to answer  -1.350e+00  1.910e+00  -0.707 0.479679    
    ## diabetesYes                   -3.466e-01  8.712e-01  -0.398 0.690727    
    ## whr                            1.465e+00  8.070e-01   1.816 0.069422 .  
    ## sexMale                        6.331e-02  1.375e-01   0.461 0.645148    
    ## age                           -1.069e-02  6.540e-03  -1.634 0.102234    
    ## pm25_val                       1.131e-01  3.964e-02   2.853 0.004331 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1997.4  on 1449  degrees of freedom
    ## Residual deviance: 1952.5  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 1986.5
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
car::vif(ukb_covid_pm25.b)
```

    ##                        GVIF Df GVIF^(1/(2*Df))
    ## merged_popDens_km2 2.324259  1        1.524552
    ## n_cancers          1.020306  1        1.010102
    ## townsend           1.406078  1        1.185782
    ## smoking            1.361982  3        1.052839
    ## whistling          1.946471  3        1.117398
    ## diabetes           2.222505  3        1.142371
    ## whr                1.864962  1        1.365636
    ## sex                1.637204  1        1.279533
    ## age                1.132009  1        1.063959
    ## pm25_val           1.946202  1        1.395064

``` r
#PM10 
summary(ukb_covid_pm10.b <-glm(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                 whr + sex + age + pm10_val, family = 'binomial'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + pm10_val, 
    ##     family = "binomial", data = ukb_covid_allPol_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.5265  -1.1001  -0.8752   1.2106   1.7815  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                   -2.436e+00  1.204e+00  -2.023 0.043114 *  
    ## merged_popDens_km2            -1.452e-05  2.764e-05  -0.525 0.599243    
    ## n_cancers                     -2.270e-01  1.656e-01  -1.371 0.170288    
    ## townsend                       2.519e-02  1.802e-02   1.398 0.162234    
    ## smokingNever                   6.253e-01  1.758e-01   3.556 0.000376 ***
    ## smokingPrefer not to answer    1.503e+00  7.454e-01   2.017 0.043730 *  
    ## smokingPrevious                6.950e-01  1.761e-01   3.946 7.95e-05 ***
    ## whistlingNo                    4.602e-01  3.491e-01   1.318 0.187409    
    ## whistlingPrefer not to answer  1.219e-01  1.635e+00   0.075 0.940573    
    ## whistlingYes                   4.382e-01  3.554e-01   1.233 0.217639    
    ## diabetesNo                    -2.959e-01  8.585e-01  -0.345 0.730345    
    ## diabetesPrefer not to answer  -1.361e+00  1.909e+00  -0.713 0.475962    
    ## diabetesYes                   -3.652e-01  8.709e-01  -0.419 0.674974    
    ## whr                            1.483e+00  8.066e-01   1.839 0.065893 .  
    ## sexMale                        6.211e-02  1.374e-01   0.452 0.651265    
    ## age                           -1.075e-02  6.536e-03  -1.645 0.099922 .  
    ## pm10_val                       7.170e-02  2.826e-02   2.537 0.011168 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1997.4  on 1449  degrees of freedom
    ## Residual deviance: 1954.2  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 1988.2
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
car::vif(ukb_covid_pm10.b)
```

    ##                        GVIF Df GVIF^(1/(2*Df))
    ## merged_popDens_km2 2.479156  1        1.574534
    ## n_cancers          1.020245  1        1.010072
    ## townsend           1.405143  1        1.185387
    ## smoking            1.359251  3        1.052487
    ## whistling          1.938463  3        1.116631
    ## diabetes           2.215023  3        1.141729
    ## whr                1.865545  1        1.365850
    ## sex                1.637883  1        1.279798
    ## age                1.132501  1        1.064190
    ## pm10_val           2.089322  1        1.445449

``` r
#nox
summary(ukb_covid_nox.b <-glm(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                whr + sex + age + nox_val, family = 'binomial'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + nox_val, 
    ##     family = "binomial", data = ukb_covid_allPol_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.8158  -1.0981  -0.8681   1.2109   1.7967  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                   -1.971e+00  1.167e+00  -1.688 0.091373 .  
    ## merged_popDens_km2            -2.821e-05  3.058e-05  -0.922 0.356330    
    ## n_cancers                     -2.220e-01  1.656e-01  -1.341 0.179989    
    ## townsend                       1.120e-02  1.841e-02   0.608 0.542937    
    ## smokingNever                   6.121e-01  1.758e-01   3.483 0.000497 ***
    ## smokingPrefer not to answer    1.418e+00  7.419e-01   1.911 0.056052 .  
    ## smokingPrevious                6.853e-01  1.761e-01   3.893 9.92e-05 ***
    ## whistlingNo                    4.534e-01  3.511e-01   1.291 0.196604    
    ## whistlingPrefer not to answer  3.230e-02  1.630e+00   0.020 0.984188    
    ## whistlingYes                   4.300e-01  3.574e-01   1.203 0.228927    
    ## diabetesNo                    -3.427e-01  8.598e-01  -0.399 0.690220    
    ## diabetesPrefer not to answer  -1.409e+00  1.907e+00  -0.739 0.460161    
    ## diabetesYes                   -4.073e-01  8.721e-01  -0.467 0.640507    
    ## whr                            1.558e+00  8.069e-01   1.930 0.053549 .  
    ## sexMale                        3.870e-02  1.378e-01   0.281 0.778826    
    ## age                           -1.023e-02  6.537e-03  -1.564 0.117746    
    ## nox_val                        2.286e-02  8.580e-03   2.665 0.007708 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1997.4  on 1449  degrees of freedom
    ## Residual deviance: 1953.4  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 1987.4
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
#no2
summary(ukb_covid_no2.b <-glm(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                whr + sex + age + no2_val, family = 'binomial'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + no2_val, 
    ##     family = "binomial", data = ukb_covid_allPol_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.7152  -1.0980  -0.8676   1.2074   1.8112  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                   -2.146e+00  1.175e+00  -1.826 0.067787 .  
    ## merged_popDens_km2            -3.474e-05  3.012e-05  -1.153 0.248722    
    ## n_cancers                     -2.206e-01  1.657e-01  -1.332 0.182872    
    ## townsend                       9.422e-03  1.844e-02   0.511 0.609470    
    ## smokingNever                   6.136e-01  1.759e-01   3.489 0.000485 ***
    ## smokingPrefer not to answer    1.414e+00  7.420e-01   1.906 0.056635 .  
    ## smokingPrevious                6.904e-01  1.762e-01   3.918 8.94e-05 ***
    ## whistlingNo                    4.519e-01  3.508e-01   1.288 0.197732    
    ## whistlingPrefer not to answer  5.918e-03  1.629e+00   0.004 0.997102    
    ## whistlingYes                   4.301e-01  3.571e-01   1.204 0.228481    
    ## diabetesNo                    -3.365e-01  8.622e-01  -0.390 0.696330    
    ## diabetesPrefer not to answer  -1.424e+00  1.908e+00  -0.746 0.455449    
    ## diabetesYes                   -4.023e-01  8.745e-01  -0.460 0.645508    
    ## whr                            1.520e+00  8.073e-01   1.883 0.059643 .  
    ## sexMale                        4.126e-02  1.378e-01   0.300 0.764553    
    ## age                           -1.014e-02  6.542e-03  -1.550 0.121023    
    ## no2_val                        4.586e-02  1.507e-02   3.042 0.002350 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1997.4  on 1449  degrees of freedom
    ## Residual deviance: 1951.3  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 1985.3
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
car::vif(ukb_covid_no2.b)
```

    ##                        GVIF Df GVIF^(1/(2*Df))
    ## merged_popDens_km2 2.908510  1        1.705436
    ## n_cancers          1.020530  1        1.010213
    ## townsend           1.465272  1        1.210484
    ## smoking            1.361503  3        1.052777
    ## whistling          1.873400  3        1.110295
    ## diabetes           2.153442  3        1.136376
    ## whr                1.867270  1        1.366481
    ## sex                1.643660  1        1.282053
    ## age                1.131877  1        1.063897
    ## no2_val            3.004997  1        1.733493

``` r
#o3
summary(ukb_covid_o3.b <-glm(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                               whr + sex + age + o3_val, family = 'binomial'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + o3_val, 
    ##     family = "binomial", data = ukb_covid_allPol_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.5128  -1.1014  -0.8815   1.2147   1.7384  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                   -1.691e+00  1.160e+00  -1.458 0.144974    
    ## merged_popDens_km2             3.266e-05  2.026e-05   1.613 0.106848    
    ## n_cancers                     -2.303e-01  1.657e-01  -1.390 0.164576    
    ## townsend                       2.509e-02  1.857e-02   1.351 0.176671    
    ## smokingNever                   6.132e-01  1.757e-01   3.491 0.000482 ***
    ## smokingPrefer not to answer    1.475e+00  7.434e-01   1.985 0.047180 *  
    ## smokingPrevious                6.728e-01  1.758e-01   3.828 0.000129 ***
    ## whistlingNo                    4.292e-01  3.476e-01   1.235 0.216935    
    ## whistlingPrefer not to answer  7.931e-02  1.637e+00   0.048 0.961358    
    ## whistlingYes                   4.005e-01  3.538e-01   1.132 0.257619    
    ## diabetesNo                    -3.154e-01  8.525e-01  -0.370 0.711418    
    ## diabetesPrefer not to answer  -1.337e+00  1.909e+00  -0.700 0.483791    
    ## diabetesYes                   -3.952e-01  8.650e-01  -0.457 0.647746    
    ## whr                            1.542e+00  8.047e-01   1.917 0.055253 .  
    ## sexMale                        6.833e-02  1.373e-01   0.498 0.618821    
    ## age                           -1.036e-02  6.521e-03  -1.589 0.112170    
    ## o3_val                         1.394e-02  2.102e-02   0.663 0.507271    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1997.4  on 1449  degrees of freedom
    ## Residual deviance: 1960.2  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 1994.2
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
car::vif(ukb_covid_o3.b)
```

    ##                        GVIF Df GVIF^(1/(2*Df))
    ## merged_popDens_km2 1.334395  1        1.155160
    ## n_cancers          1.020762  1        1.010327
    ## townsend           1.497646  1        1.223783
    ## smoking            1.363683  3        1.053058
    ## whistling          1.938926  3        1.116675
    ## diabetes           2.218764  3        1.142050
    ## whr                1.865199  1        1.365723
    ## sex                1.642995  1        1.281794
    ## age                1.132566  1        1.064221
    ## o3_val             1.092175  1        1.045072

``` r
#so2
summary(ukb_covid_so2.b <-glm(data = ukb_covid_allPol_df, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+
                                whr + sex + age + so2_val, family = 'binomial'))
```

    ## 
    ## Call:
    ## glm(formula = result ~ merged_popDens_km2 + n_cancers + townsend + 
    ##     smoking + whistling + diabetes + whr + sex + age + so2_val, 
    ##     family = "binomial", data = ukb_covid_allPol_df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.5419  -1.1008  -0.8828   1.2125   1.8154  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                   -1.865e+00  1.168e+00  -1.597 0.110294    
    ## merged_popDens_km2             2.955e-05  2.042e-05   1.447 0.147827    
    ## n_cancers                     -2.274e-01  1.657e-01  -1.372 0.170114    
    ## townsend                       1.531e-02  1.853e-02   0.826 0.408529    
    ## smokingNever                   6.105e-01  1.755e-01   3.478 0.000505 ***
    ## smokingPrefer not to answer    1.415e+00  7.429e-01   1.905 0.056727 .  
    ## smokingPrevious                6.719e-01  1.756e-01   3.827 0.000130 ***
    ## whistlingNo                    4.384e-01  3.479e-01   1.260 0.207518    
    ## whistlingPrefer not to answer  2.860e-02  1.632e+00   0.018 0.986015    
    ## whistlingYes                   4.079e-01  3.540e-01   1.152 0.249305    
    ## diabetesNo                    -3.297e-01  8.535e-01  -0.386 0.699264    
    ## diabetesPrefer not to answer  -1.338e+00  1.905e+00  -0.703 0.482347    
    ## diabetesYes                   -4.055e-01  8.659e-01  -0.468 0.639609    
    ## whr                            1.542e+00  8.050e-01   1.916 0.055368 .  
    ## sexMale                        4.830e-02  1.377e-01   0.351 0.725700    
    ## age                           -1.015e-02  6.524e-03  -1.556 0.119623    
    ## so2_val                        1.553e-01  1.078e-01   1.441 0.149709    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1997.4  on 1449  degrees of freedom
    ## Residual deviance: 1958.6  on 1433  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 1992.6
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
car::vif(ukb_covid_so2.b)
```

    ##                        GVIF Df GVIF^(1/(2*Df))
    ## merged_popDens_km2 1.352441  1        1.162945
    ## n_cancers          1.020029  1        1.009965
    ## townsend           1.489026  1        1.220257
    ## smoking            1.355795  3        1.052040
    ## whistling          1.916015  3        1.114465
    ## diabetes           2.192588  3        1.139793
    ## whr                1.865722  1        1.365914
    ## sex                1.649086  1        1.284167
    ## age                1.131975  1        1.063943
    ## so2_val            1.177317  1        1.085042

### Supplementary Table 6. Effect of air pollutants on the probability of being infected, using the UK Biobank data

Add <br> plot binary modelsAdd <br>

``` r
stargazer(ukb_covid_pm25.b, ukb_covid_pm10.b, ukb_covid_nox.b, ukb_covid_no2.b,
          ukb_covid_o3.b, ukb_covid_so2.b,type="html",out = "fig_out_v4/ukb_covid_allPoll_binary.html",
          dep.var.labels="COVID positive or not",
          single.row=TRUE)
```

<table style="text-align:center">

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="6">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="6" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="6">

COVID positive or not

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

(1)

</td>

<td>

(2)

</td>

<td>

(3)

</td>

<td>

(4)

</td>

<td>

(5)

</td>

<td>

(6)

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

merged\_popDens\_km2

</td>

<td>

\-0.00002 (0.00003)

</td>

<td>

\-0.00001 (0.00003)

</td>

<td>

\-0.00003 (0.00003)

</td>

<td>

\-0.00003 (0.00003)

</td>

<td>

0.00003 (0.00002)

</td>

<td>

0.00003 (0.00002)

</td>

</tr>

<tr>

<td style="text-align:left">

n\_cancers

</td>

<td>

\-0.228 (0.166)

</td>

<td>

\-0.227 (0.166)

</td>

<td>

\-0.222 (0.166)

</td>

<td>

\-0.221 (0.166)

</td>

<td>

\-0.230 (0.166)

</td>

<td>

\-0.227 (0.166)

</td>

</tr>

<tr>

<td style="text-align:left">

townsend

</td>

<td>

0.026 (0.018)

</td>

<td>

0.025 (0.018)

</td>

<td>

0.011 (0.018)

</td>

<td>

0.009 (0.018)

</td>

<td>

0.025 (0.019)

</td>

<td>

0.015 (0.019)

</td>

</tr>

<tr>

<td style="text-align:left">

smokingNever

</td>

<td>

0.632<sup>\*\*\*</sup> (0.176)

</td>

<td>

0.625<sup>\*\*\*</sup> (0.176)

</td>

<td>

0.612<sup>\*\*\*</sup> (0.176)

</td>

<td>

0.614<sup>\*\*\*</sup> (0.176)

</td>

<td>

0.613<sup>\*\*\*</sup> (0.176)

</td>

<td>

0.610<sup>\*\*\*</sup> (0.176)

</td>

</tr>

<tr>

<td style="text-align:left">

smokingPrefer not to answer

</td>

<td>

1.517<sup>\*\*</sup> (0.746)

</td>

<td>

1.503<sup>\*\*</sup> (0.745)

</td>

<td>

1.418<sup>\*</sup> (0.742)

</td>

<td>

1.414<sup>\*</sup> (0.742)

</td>

<td>

1.475<sup>\*\*</sup> (0.743)

</td>

<td>

1.415<sup>\*</sup> (0.743)

</td>

</tr>

<tr>

<td style="text-align:left">

smokingPrevious

</td>

<td>

0.705<sup>\*\*\*</sup> (0.176)

</td>

<td>

0.695<sup>\*\*\*</sup> (0.176)

</td>

<td>

0.685<sup>\*\*\*</sup> (0.176)

</td>

<td>

0.690<sup>\*\*\*</sup> (0.176)

</td>

<td>

0.673<sup>\*\*\*</sup> (0.176)

</td>

<td>

0.672<sup>\*\*\*</sup> (0.176)

</td>

</tr>

<tr>

<td style="text-align:left">

whistlingNo

</td>

<td>

0.473 (0.349)

</td>

<td>

0.460 (0.349)

</td>

<td>

0.453 (0.351)

</td>

<td>

0.452 (0.351)

</td>

<td>

0.429 (0.348)

</td>

<td>

0.438 (0.348)

</td>

</tr>

<tr>

<td style="text-align:left">

whistlingPrefer not to answer

</td>

<td>

0.123 (1.636)

</td>

<td>

0.122 (1.635)

</td>

<td>

0.032 (1.630)

</td>

<td>

0.006 (1.629)

</td>

<td>

0.079 (1.637)

</td>

<td>

0.029 (1.632)

</td>

</tr>

<tr>

<td style="text-align:left">

whistlingYes

</td>

<td>

0.448 (0.356)

</td>

<td>

0.438 (0.355)

</td>

<td>

0.430 (0.357)

</td>

<td>

0.430 (0.357)

</td>

<td>

0.401 (0.354)

</td>

<td>

0.408 (0.354)

</td>

</tr>

<tr>

<td style="text-align:left">

diabetesNo

</td>

<td>

\-0.279 (0.859)

</td>

<td>

\-0.296 (0.858)

</td>

<td>

\-0.343 (0.860)

</td>

<td>

\-0.337 (0.862)

</td>

<td>

\-0.315 (0.853)

</td>

<td>

\-0.330 (0.853)

</td>

</tr>

<tr>

<td style="text-align:left">

diabetesPrefer not to answer

</td>

<td>

\-1.350 (1.910)

</td>

<td>

\-1.361 (1.909)

</td>

<td>

\-1.409 (1.907)

</td>

<td>

\-1.424 (1.908)

</td>

<td>

\-1.337 (1.909)

</td>

<td>

\-1.338 (1.905)

</td>

</tr>

<tr>

<td style="text-align:left">

diabetesYes

</td>

<td>

\-0.347 (0.871)

</td>

<td>

\-0.365 (0.871)

</td>

<td>

\-0.407 (0.872)

</td>

<td>

\-0.402 (0.875)

</td>

<td>

\-0.395 (0.865)

</td>

<td>

\-0.405 (0.866)

</td>

</tr>

<tr>

<td style="text-align:left">

whr

</td>

<td>

1.465<sup>\*</sup> (0.807)

</td>

<td>

1.483<sup>\*</sup> (0.807)

</td>

<td>

1.558<sup>\*</sup> (0.807)

</td>

<td>

1.520<sup>\*</sup> (0.807)

</td>

<td>

1.542<sup>\*</sup> (0.805)

</td>

<td>

1.542<sup>\*</sup> (0.805)

</td>

</tr>

<tr>

<td style="text-align:left">

sexMale

</td>

<td>

0.063 (0.137)

</td>

<td>

0.062 (0.137)

</td>

<td>

0.039 (0.138)

</td>

<td>

0.041 (0.138)

</td>

<td>

0.068 (0.137)

</td>

<td>

0.048 (0.138)

</td>

</tr>

<tr>

<td style="text-align:left">

age

</td>

<td>

\-0.011 (0.007)

</td>

<td>

\-0.011<sup>\*</sup> (0.007)

</td>

<td>

\-0.010 (0.007)

</td>

<td>

\-0.010 (0.007)

</td>

<td>

\-0.010 (0.007)

</td>

<td>

\-0.010 (0.007)

</td>

</tr>

<tr>

<td style="text-align:left">

pm25\_val

</td>

<td>

0.113<sup>\*\*\*</sup> (0.040)

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

pm10\_val

</td>

<td>

</td>

<td>

0.072<sup>\*\*</sup> (0.028)

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

nox\_val

</td>

<td>

</td>

<td>

</td>

<td>

0.023<sup>\*\*\*</sup> (0.009)

</td>

<td>

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

no2\_val

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

0.046<sup>\*\*\*</sup> (0.015)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

o3\_val

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

0.014 (0.021)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

so2\_val

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

0.155 (0.108)

</td>

</tr>

<tr>

<td style="text-align:left">

Constant

</td>

<td>

\-2.495<sup>\*\*</sup> (1.201)

</td>

<td>

\-2.436<sup>\*\*</sup> (1.204)

</td>

<td>

\-1.971<sup>\*</sup> (1.167)

</td>

<td>

\-2.146<sup>\*</sup> (1.175)

</td>

<td>

\-1.691 (1.160)

</td>

<td>

\-1.865 (1.168)

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Observations

</td>

<td>

1,450

</td>

<td>

1,450

</td>

<td>

1,450

</td>

<td>

1,450

</td>

<td>

1,450

</td>

<td>

1,450

</td>

</tr>

<tr>

<td style="text-align:left">

Log Likelihood

</td>

<td>

\-976.235

</td>

<td>

\-977.093

</td>

<td>

\-976.681

</td>

<td>

\-975.636

</td>

<td>

\-980.105

</td>

<td>

\-979.286

</td>

</tr>

<tr>

<td style="text-align:left">

Akaike Inf. Crit.

</td>

<td>

1,986.469

</td>

<td>

1,988.185

</td>

<td>

1,987.362

</td>

<td>

1,985.271

</td>

<td>

1,994.211

</td>

<td>

1,992.572

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="6" style="text-align:right">

<sup>*</sup>p\<0.1; <sup>**</sup>p\<0.05; <sup>***</sup>p\<0.01

</td>

</tr>

</table>

<br> Each column of the table corresponds to a binary model used to
predict whether one would get infected based on the concentration of an
air pollutant, as well as covariates. We used different models to
characterise each pollutant because the pollutants are highly
correlated. The raw estimate values of each model are listed with their
standard error in parentheses. The p-values are indicated using the
number of asterisks beside the estimates. Average\_Pop\_densitykm2,
average population density per square kilometer; NO.levels, nitrogen
oxides levels; no2\_val, nitrogen dioxide levels; O3.levels, ozone
levels; Akaike Inf. Crit., Akaike’s Information Criteria

### Calculate infectivity odds ratios

### Format odds ratios

``` r
ukb_covid_pm25.b_odds = data.frame(cbind(exp(cbind(OR = coef(ukb_covid_pm25.b), confint(ukb_covid_pm25.b))), p_value = summary(ukb_covid_pm25.b)$coefficients[,4]))
ukb_covid_pm10.b_odds = data.frame(cbind(exp(cbind(OR = coef(ukb_covid_pm10.b), confint(ukb_covid_pm10.b))), p_value = summary(ukb_covid_pm10.b)$coefficients[,4]))
ukb_covid_nox.b_odds = data.frame(cbind(exp(cbind(OR = coef(ukb_covid_nox.b), confint(ukb_covid_nox.b))), p_value = summary(ukb_covid_nox.b)$coefficients[,4]))
ukb_covid_no2.b_odds = data.frame(cbind(exp(cbind(OR = coef(ukb_covid_no2.b), confint(ukb_covid_no2.b))), p_value = summary(ukb_covid_no2.b)$coefficients[,4]))
ukb_covid_o3.b_odds = data.frame(cbind(exp(cbind(OR = coef(ukb_covid_o3.b), confint(ukb_covid_o3.b))), p_value = summary(ukb_covid_o3.b)$coefficients[,4]))
ukb_covid_so2.b_odds = data.frame(cbind(exp(cbind(OR = coef(ukb_covid_so2.b), confint(ukb_covid_so2.b))), p_value = summary(ukb_covid_so2.b)$coefficients[,4]))
#make a df just with the pollutants 

ukb_covid_onlyPoll = data.frame(rbind(ukb_covid_pm25.b_odds[nrow(ukb_covid_pm25.b_odds),],
                           ukb_covid_pm10.b_odds[nrow(ukb_covid_pm10.b_odds),],
                           ukb_covid_nox.b_odds[nrow(ukb_covid_nox.b_odds),],
                           ukb_covid_no2.b_odds[nrow(ukb_covid_no2.b_odds),],
                           ukb_covid_o3.b_odds[nrow(ukb_covid_o3.b_odds),],
                           ukb_covid_so2.b_odds[nrow(ukb_covid_so2.b_odds),]))
```

#### Plot infectivity odds ratios

sort data

``` r
ukb_covid_onlyPoll$names = row.names(ukb_covid_onlyPoll)
ukb_covid_onlyPoll$significance = "p-value > 0.05"
ukb_covid_onlyPoll$significance[ukb_covid_onlyPoll$p_value < 0.05] <- "p-value < 0.05"
```

plot

``` r
ggplot(ukb_covid_onlyPoll, aes(x=reorder(names, OR), y=OR, color=significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_bw() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("Infectivity odds ratios") + 
  xlab("Pollutants")+
  ylim(0.75, 1.6)
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->

``` r
ggsave('fig_out_v4/UKB_infectivityOR_bar.pdf')
```

``` r
stargazer(ukb_covid_onlyPoll, type ="html", single.row=TRUE, summary = FALSE, out = "fig_out_v4/ukb_covid_onlyPoll.html")
```

<table style="text-align:center">

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

OR

</td>

<td>

X2.5..

</td>

<td>

X97.5..

</td>

<td>

p\_value

</td>

<td>

names

</td>

<td>

significance

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

pm25\_val

</td>

<td>

1.120

</td>

<td>

1.036

</td>

<td>

1.211

</td>

<td>

0.004

</td>

<td>

pm25\_val

</td>

<td>

p-value \< 0.05

</td>

</tr>

<tr>

<td style="text-align:left">

pm10\_val

</td>

<td>

1.074

</td>

<td>

1.017

</td>

<td>

1.136

</td>

<td>

0.011

</td>

<td>

pm10\_val

</td>

<td>

p-value \< 0.05

</td>

</tr>

<tr>

<td style="text-align:left">

nox\_val

</td>

<td>

1.023

</td>

<td>

1.006

</td>

<td>

1.041

</td>

<td>

0.008

</td>

<td>

nox\_val

</td>

<td>

p-value \< 0.05

</td>

</tr>

<tr>

<td style="text-align:left">

no2\_val

</td>

<td>

1.047

</td>

<td>

1.017

</td>

<td>

1.079

</td>

<td>

0.002

</td>

<td>

no2\_val

</td>

<td>

p-value \< 0.05

</td>

</tr>

<tr>

<td style="text-align:left">

o3\_val

</td>

<td>

1.014

</td>

<td>

0.973

</td>

<td>

1.057

</td>

<td>

0.507

</td>

<td>

o3\_val

</td>

<td>

p-value \> 0.05

</td>

</tr>

<tr>

<td style="text-align:left">

so2\_val

</td>

<td>

1.168

</td>

<td>

0.946

</td>

<td>

1.444

</td>

<td>

0.150

</td>

<td>

so2\_val

</td>

<td>

p-value \> 0.05

</td>

</tr>

<tr>

<td colspan="7" style="border-bottom: 1px solid black">

</td>

</tr>

</table>

## Aim 4: Identifying the main contributors of air pollution

We downloaded the original files in data\_v4. <br> We then deleted the
irrelevant or repeated columns by hand, and saved in data output.<br>

#### Analysis rationale:

First, we set out to explore how the main contributors of NOx and SO2
emissions relate with each other. <br>Then we wanted to assess which
emission source is best at predicting air pollution.<br>

To deal with the multicollinearity and high dimensionality of these
confounding variables we used two different approaches, based on
previous air pollution studies
(<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6717171/> and
<https://www.sciencedirect.com/science/article/pii/S0048969714016647?via%3Dihub>).
First, we used a principal component analysis (PCA) approach to
transform the original variables into a set of uncorrelated, mutually
orthogonal components called principal components (PCs). PCs represent
linear combinations of the original variables that provide information
on the most meaningful parameter with minimum loss of original
information. <br> PCs were not estimated for gas consumption as only two
variables were included in this analysis. <br>

While effective in showing the overall contribution of each polluting
category, major limitations of PCA include the difficulty to interpret
the results as the components are not in the same unit as the original
variables and the failure to account for the relationship between
exposure and response variables. This is why we also use an interative
approach. <br> \#\#\#\# Supplementary Table 7. Main nitrogen oxides
(NOx) emission sources Add <br>

``` r
NOx_emissions <- read.csv("data_v4/Main_NOx_emission_sector.csv")
library(ggplot2)
NOx.plot <- ggplot(data = NOx_emissions, mapping = aes(y = Emissions, x = Year, fill = Sector)) +
  geom_bar(stat = "identity", position = "fill", width = 8) +
  theme_gray() + scale_fill_brewer() +
  ylab("Emissions of NOx (thousand of tonnes)") +
  theme(
    axis.title.x = element_text(color="black", vjust=-0.35, face = "bold"),
    axis.title.y = element_text(color="black" , vjust=3, face = "bold"),
    axis.text.x = element_text(angle=65, vjust=0.6),
    plot.title = element_text(size=12, face="bold", 
                              margin = margin(10, 10, 0, 0), hjust = 0.5,
                              vjust = 5),
    legend.title = element_blank())
# x axis
NOx.plot + scale_x_continuous(breaks = c(1990, 2005, 2018),labels=c(1990, 2005, 2018))
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->

``` r
ggsave("fig_out_v4/Main_NOx_emission_sector_V2.pdf")
```

The graph shows the contribution of nitrogen oxides emissions for each
sector in the UK over the last three decades. “Road transport” refers to
cumulative emissions recorded for various motor vehicles (buses,
motorcycles, cars, LGV and HGV); “Other category” captures nitrogen
oxides emissions from non-agriculture livestock such as domestic pets as
well as non-agriculture fertilisers; “Non-road transport” represents
emissions from residual transport sources (aviation, rail, shipping);
“Manufacturing industries and construction” refers to cement and lime
production as well as other industrial combustion; and “Energy
industries” included emissions produced by power stations, miscellaneous
commercial and industrial combustion, gas and coke production. <br>

#### Supplementary Table 8. Main SO2 emission sources

Add <br>

``` r
So2_emissions <- read.csv("data_v4/Figure_SO2_sector.csv")
library(ggplot2)
SO2.plot <- ggplot(data = So2_emissions, mapping = aes(y = Emissions, x = Year, fill = Sector)) +
  geom_bar(stat = "identity", position = "fill", width = 8) +
  theme_gray() + scale_fill_brewer() +
  ylab("Emissions of sulphur dioxide (thousand of tonnes) in %") +
  theme(
    axis.title.x = element_text(color="black", vjust=-0.35, face = "bold"),
    axis.title.y = element_text(color="black" , vjust=3, face = "bold"),
    axis.text.x = element_text(angle=65, vjust=0.6),
    plot.title = element_text(size=12, face="bold", 
                              margin = margin(10, 10, 0, 0), hjust = 0.5,
                              vjust = 5),
    legend.title = element_blank())
# x axis
SO2.plot + scale_x_continuous(breaks = c(1990, 2005, 2018),labels=c(1990, 2005, 2018))
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

``` r
ggsave("fig_out_v4/Main_SO2_emission_sector_V2.pdf")
```

This graph shows the contribution of sulphur dioxide emissions for each
sector in the UK over the last three decades. For the explanation of
emission sources, refer to Supplementary Table 7 <br>

### Load and process the data

``` r
LA_merged_dt = read.csv("data_output_v4/merged_covidAir_cov_dt_LA.csv")
LA_merged_code_dt = LA_merged_dt["Code"]
colnames(LA_merged_code_dt) <- "ons_code"
road_fuel_raw = read.csv("data_output_v4/clean_road_fuel_consumption2017.csv")
residuals_raw = read.csv("data_output_v4/clean_residual_consumption2017.csv", na.strings = "NA")
gas_raw = read.csv("data_output_v4/clean_gas_consumption2018.csv", na.strings = "NA")

print(paste("road_fuel_raw: ", colnames(road_fuel_raw)))
```

    ##  [1] "road_fuel_raw:  ons_code"                        
    ##  [2] "road_fuel_raw:  area"                            
    ##  [3] "road_fuel_raw:  personal.buses.Motorways"        
    ##  [4] "road_fuel_raw:  personal.buses.A.roads"          
    ##  [5] "road_fuel_raw:  personal.buses.Minor.roads"      
    ##  [6] "road_fuel_raw:  personal.diesel.cars.Motorways"  
    ##  [7] "road_fuel_raw:  personal.diesel.cars.A.roads"    
    ##  [8] "road_fuel_raw:  personal.diesel.cars.Minor.roads"
    ##  [9] "road_fuel_raw:  personal.petrol.cars.Motorway"   
    ## [10] "road_fuel_raw:  personal.petrol.cars.A.roads"    
    ## [11] "road_fuel_raw:  personal.petrol.cars.Minor.roads"
    ## [12] "road_fuel_raw:  personal.motorcycles.Motorways"  
    ## [13] "road_fuel_raw:  personal.motorcycles.A.roads"    
    ## [14] "road_fuel_raw:  personal.motorcycles.Minor.roads"
    ## [15] "road_fuel_raw:  freight.HGV.Motorways"           
    ## [16] "road_fuel_raw:  freight.HGV.A.roads"             
    ## [17] "road_fuel_raw:  freight.HGV.Minor.roads"         
    ## [18] "road_fuel_raw:  freight.diesel.LGV.Motorways"    
    ## [19] "road_fuel_raw:  freight.diesel.LGV.A.roads"      
    ## [20] "road_fuel_raw:  freight.diesel.LGV.Minor.roads"  
    ## [21] "road_fuel_raw:  freight.petrol.LGV.Motorways"    
    ## [22] "road_fuel_raw:  freight.petrol.LGV.A.roads"      
    ## [23] "road_fuel_raw:  freight.petrol.LGV.Minor.roads"

``` r
print(paste("residuals_raw: ", colnames(residuals_raw)))
```

    ##  [1] "residuals_raw:  ons_code"                           
    ##  [2] "residuals_raw:  area"                               
    ##  [3] "residuals_raw:  petroleum.Industrial"               
    ##  [4] "residuals_raw:  petroleum.Domestic"                 
    ##  [5] "residuals_raw:  petroleum.Rail"                     
    ##  [6] "residuals_raw:  petroleum.Public.Administration"    
    ##  [7] "residuals_raw:  petroleum.Commercial"               
    ##  [8] "residuals_raw:  petroleum.Agriculture"              
    ##  [9] "residuals_raw:  coal.Industrial.Commercial"         
    ## [10] "residuals_raw:  coal.Domestic"                      
    ## [11] "residuals_raw:  coal.Rail"                          
    ## [12] "residuals_raw:  Manufactured.Solid.Fuels.Industrial"
    ## [13] "residuals_raw:  Manufactured.Solid.Fuels.Domestic"  
    ## [14] "residuals_raw:  bioenergy"

Only keep the data with location in England and delete columns that
contain more than 10% NA. <br>

``` r
road_fuel_dt = merge(LA_merged_code_dt, road_fuel_raw, by = "ons_code")
residuals_dt = merge(LA_merged_code_dt, residuals_raw, by = "ons_code")
gas_dt = merge(LA_merged_code_dt, gas_raw, by = "ons_code")

#keep track of deleted columns
road_fuel_rm = road_fuel_dt[,!sapply(road_fuel_dt, function(x) mean(is.na(x)))<0.1]
residuals_rm = residuals_dt[,!sapply(residuals_dt, function(x) mean(is.na(x)))<0.1]
gas_rm = gas_dt[,!sapply(gas_dt, function(x) mean(is.na(x)))<0.1]

#delete columns that contain more than 10% NA
road_fuel_dt = road_fuel_dt[,!sapply(road_fuel_dt, function(x) mean(is.na(x)))>0.1]
residuals_dt = residuals_dt[,!sapply(residuals_dt, function(x) mean(is.na(x)))>0.1]
gas_dt = gas_dt[,!sapply(gas_dt, function(x) mean(is.na(x)))>0.1]

# save 
write.csv(road_fuel_dt, "data_output_v4/processed_road_fuel_consumption2017.csv")
write.csv(residuals_dt, "data_output_v4/processed_residual_consumption2017.csv")
write.csv(gas_dt, "data_output_v4/processed_road_gas_consumption2018.csv")
```

### PCA for variable grouping and multicollinearity issues

Read in and view the data (England only)

``` r
fossilfuel_2017 = read.csv("data_output_v4/AQ_vs_fossilfuels_LA_England.csv")
head(fossilfuel_2017)
```

    ##   Number  ons_code     Region X.people.per.sq.km      Local.Authority   no2_val
    ## 1      1 E06000001 North East                997           Hartlepool 13.105737
    ## 2      2 E06000002 North East               2608        Middlesbrough 20.178191
    ## 3      3 E06000003 North East                558 Redcar and Cleveland  7.749542
    ## 4      4 E06000004 North East                962     Stockton-on-Tees 15.562038
    ## 5      5 E06000005 North East                540           Darlington 11.360341
    ## 6      6 E06000006 North West               1624               Halton 15.611062
    ##     o3_val   nox_val pm25_val pm10_val   so2_val Domestic.mean.gas.consumption
    ## 1 5.303858 17.730895 7.362821 11.50678 1.5755607                         12439
    ## 2 4.251802 29.296034 8.210541 12.56193 3.1332865                         12573
    ## 3 6.851500  9.983241 7.293137 12.85583 0.9988638                         12519
    ## 4 5.050632 21.440737 7.900856 12.25113 1.9343295                         12992
    ## 5 5.198811 15.098272 7.153221 11.12964 1.1290335                         13417
    ## 6 4.393569 21.642167 9.257343 13.26975 3.0365011                         11520
    ##   Non.domestic.mean.gas.consumption personal.buses.Motorways
    ## 1                            698805                        0
    ## 2                            531688                        0
    ## 3                           2904401                        0
    ## 4                            926081                        0
    ## 5                           1039749                       89
    ## 6                           1075639                      120
    ##   personal.buses.A.roads personal.buses.Minor.roads
    ## 1                    545                       1627
    ## 2                    909                       3776
    ## 3                    938                       2473
    ## 4                   1518                       3389
    ## 5                    695                       1808
    ## 6                    481                        443
    ##   personal.diesel.cars.Motorways personal.diesel.cars.A.roads
    ## 1                              0                         7470
    ## 2                              0                        13295
    ## 3                              0                        10349
    ## 4                              0                        20375
    ## 5                           3232                         6702
    ## 6                           5026                         9109
    ##   personal.diesel.cars.Minor.roads personal.petrol.cars.Motorway
    ## 1                             5526                             0
    ## 2                            12563                             0
    ## 3                             8236                             0
    ## 4                            11383                             0
    ## 5                             6224                          2268
    ## 6                             6392                          3552
    ##   personal.petrol.cars.A.roads personal.petrol.cars.Minor.roads
    ## 1                         8479                             6810
    ## 2                        16573                            15706
    ## 3                        11578                            10048
    ## 4                        24160                            13982
    ## 5                         7488                             7395
    ## 6                        11779                             7945
    ##   personal.motorcycles.Motorways personal.motorcycles.A.roads
    ## 1                              0                           50
    ## 2                              0                           84
    ## 3                              0                          114
    ## 4                              0                          141
    ## 5                             22                           68
    ## 6                             23                           97
    ##   personal.motorcycles.Minor.roads freight.HGV.Motorways freight.HGV.A.roads
    ## 1                               78                     0                5470
    ## 2                              176                     0                8967
    ## 3                              119                     0                5558
    ## 4                              162                     0               14478
    ## 5                               90                  5160                3752
    ## 6                              102                  6443                9285
    ##   freight.HGV.Minor.roads freight.diesel.LGV.Motorways
    ## 1                     446                            0
    ## 2                    1018                            0
    ## 3                     678                            0
    ## 4                     927                            0
    ## 5                     506                         1955
    ## 6                     538                         2837
    ##   freight.diesel.LGV.A.roads freight.diesel.LGV.Minor.roads
    ## 1                       4405                           2493
    ## 2                       7272                           5693
    ## 3                       5289                           3802
    ## 4                      11586                           5197
    ## 5                       3676                           2845
    ## 6                       4772                           2666
    ##   freight.petrol.LGV.Motorways freight.petrol.LGV.A.roads
    ## 1                            0                        137
    ## 2                            0                        248
    ## 3                            0                        166
    ## 4                            0                        371
    ## 5                           56                        116
    ## 6                           81                        174
    ##   freight.petrol.LGV.Minor.roads petroleum.Industrial petroleum.Domestic
    ## 1                             99                  5.2                0.4
    ## 2                            227                  3.6                0.3
    ## 3                            149                609.0                1.3
    ## 4                            205                 72.0                1.3
    ## 5                            110                  7.5                1.1
    ## 6                            106                 14.4                0.5
    ##   petroleum.Rail petroleum.Public.Administration petroleum.Commercial
    ## 1            0.2                               0                  0.0
    ## 2            0.6                               0                  0.0
    ## 3            0.9                               0                  0.1
    ## 4            2.1                               0                  0.1
    ## 5            1.9                               0                  0.1
    ## 6            6.5                               0                  0.3
    ##   petroleum.Agriculture coal.Domestic Manufactured.Solid.Fuels.Domestic
    ## 1                   0.4           0.2                               0.1
    ## 2                   0.2           0.2                               0.2
    ## 3                   1.1           0.5                               0.4
    ## 4                   0.8           0.3                               0.2
    ## 5                   1.2           0.2                               0.2
    ## 6                   0.2           0.4                               0.4
    ##   bioenergy
    ## 1       1.3
    ## 2       2.1
    ## 3       2.2
    ## 4       3.1
    ## 5       1.7
    ## 6      48.4

First thing to do here is to analyse the data distribution of the target
or response variable (NO2, NOx, O3 and SO2 concentration). To do this,
we will use the following functions:

``` r
hist(fossilfuel_2017$no2_val, main="Histogram of Yield", xlab="Yield (quintals/ha)") 
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->

``` r
#qqnorm(fossilfuel_2017$no2_val, main="QQplot of NO2")
library(e1071)
skewness(fossilfuel_2017$no2_val) 
```

    ## [1] NA

``` r
hist(fossilfuel_2017$nox_val, main="Histogram of Yield", xlab="Yield (quintals/ha)") 
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-54-2.png)<!-- -->

``` r
#qqnorm(fossilfuel_2017$nox_val, main="QQplot of NOx")
skewness(fossilfuel_2017$nox_val)
```

    ## [1] NA

``` r
hist(fossilfuel_2017$o3_val, main="Histogram of Yield", xlab="Yield (quintals/ha)") 
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-54-3.png)<!-- -->

``` r
#qqnorm(fossilfuel_2017$o3_val, main="QQplot of O3")
skewness(fossilfuel_2017$o3_val)
```

    ## [1] NA

``` r
hist(fossilfuel_2017$so2_val, main="Histogram of Yield", xlab="Yield (quintals/ha)") 
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-54-4.png)<!-- -->

``` r
#qqnorm(fossilfuel_2017$so2_val, main="QQplot of SO2")
skewness(fossilfuel_2017$so2_val) 
```

    ## [1] NA

For NO2, NOx, O3 and SO2, the histogram suggests that data is not
normally distributed.

### Data formatting for PCA

One of the fundamental assumptions of linear modelling is the
independence of predictors. Our independent variables show
multicollinearity. We will adjust for multicollinearity by using a PCA
approach. PCA is ran on the independent variables only, so we must first
subset the data to only include predictors:

``` r
# Before fititng the model, exlude missing data
fossilfuel_2017.omit <- na.omit(fossilfuel_2017)
# Create a dataset that only includes road transport independent variables
data_coll.road.clean <- fossilfuel_2017.omit[,14:34]
data.road.pca <- prcomp(data_coll.road.clean, center = TRUE, scale. = TRUE)
summary(data.road.pca)
```

    ## Importance of components:
    ##                           PC1    PC2     PC3     PC4     PC5     PC6     PC7
    ## Standard deviation     3.2200 2.4035 1.42388 1.19921 0.60912 0.53559 0.45432
    ## Proportion of Variance 0.4937 0.2751 0.09654 0.06848 0.01767 0.01366 0.00983
    ## Cumulative Proportion  0.4937 0.7688 0.86535 0.93383 0.95150 0.96516 0.97499
    ##                           PC8     PC9    PC10    PC11    PC12   PC13    PC14
    ## Standard deviation     0.3693 0.35149 0.27052 0.25962 0.22582 0.1586 0.13393
    ## Proportion of Variance 0.0065 0.00588 0.00348 0.00321 0.00243 0.0012 0.00085
    ## Cumulative Proportion  0.9815 0.98737 0.99085 0.99406 0.99649 0.9977 0.99854
    ##                           PC15    PC16    PC17    PC18    PC19     PC20
    ## Standard deviation     0.12209 0.11030 0.05132 0.02130 0.01884 0.007981
    ## Proportion of Variance 0.00071 0.00058 0.00013 0.00002 0.00002 0.000000
    ## Cumulative Proportion  0.99925 0.99983 0.99996 0.99998 1.00000 1.000000
    ##                            PC21
    ## Standard deviation     0.005863
    ## Proportion of Variance 0.000000
    ## Cumulative Proportion  1.000000

``` r
# gas
data_coll.gas.clean <- fossilfuel_2017.omit[,12:13]
data.gas.pca <- prcomp(data_coll.gas.clean, center = TRUE, scale. = TRUE)
summary(data.gas.pca)
```

    ## Importance of components:
    ##                           PC1    PC2
    ## Standard deviation     1.0833 0.9091
    ## Proportion of Variance 0.5867 0.4133
    ## Cumulative Proportion  0.5867 1.0000

``` r
# Residual fuels
data_coll.rsd.clean <- fossilfuel_2017.omit[,35:43]
data.rsd.pca <- prcomp(data_coll.rsd.clean, center = TRUE, scale. = TRUE)
summary(data.rsd.pca)
```

    ## Importance of components:
    ##                           PC1    PC2    PC3     PC4    PC5    PC6     PC7
    ## Standard deviation     2.1275 1.0342 0.9684 0.94606 0.8044 0.6782 0.56912
    ## Proportion of Variance 0.5029 0.1188 0.1042 0.09945 0.0719 0.0511 0.03599
    ## Cumulative Proportion  0.5029 0.6217 0.7259 0.82539 0.8973 0.9484 0.98438
    ##                            PC8     PC9
    ## Standard deviation     0.32078 0.19400
    ## Proportion of Variance 0.01143 0.00418
    ## Cumulative Proportion  0.99582 1.00000

### Models with PCs for fossil fuel consumption data

Now we can fit glms using PCs 1 and 2 for road transport and fossil fuel
data. Gas consumption data only includes two variables so these will be
added directly into the model. The glms employed here are consistent
with the ones used with the original variables are fitted according to
assumptions of linearity and skewness in the response variable. Add <br>

``` r
## Run model with PC1 and PC2 from each data set
# NO2 glm
summary(mod.pca.no2 <- glm( no2_val~X.people.per.sq.km + data.road.pca$x[,1:2] + data.rsd.pca$x[,1:2] + Domestic.mean.gas.consumption + Non.domestic.mean.gas.consumption,
                       family  = Gamma(link = "log"), data= fossilfuel_2017.omit)) #road and rsd is significant
```

    ## 
    ## Call:
    ## glm(formula = no2_val ~ X.people.per.sq.km + data.road.pca$x[, 
    ##     1:2] + data.rsd.pca$x[, 1:2] + Domestic.mean.gas.consumption + 
    ##     Non.domestic.mean.gas.consumption, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.43119  -0.20277   0.01158   0.15501   1.50927  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        2.437e+00  1.844e-01  13.211   <2e-16 ***
    ## X.people.per.sq.km                 9.517e-05  7.982e-06  11.923   <2e-16 ***
    ## data.road.pca$x[, 1:2]PC1          6.372e-02  7.205e-03   8.844   <2e-16 ***
    ## data.road.pca$x[, 1:2]PC2         -8.225e-03  8.167e-03  -1.007    0.315    
    ## data.rsd.pca$x[, 1:2]PC1          -1.216e-01  1.215e-02 -10.011   <2e-16 ***
    ## data.rsd.pca$x[, 1:2]PC2          -1.152e-02  1.719e-02  -0.670    0.503    
    ## Domestic.mean.gas.consumption     -5.125e-06  1.294e-05  -0.396    0.692    
    ## Non.domestic.mean.gas.consumption -4.402e-09  3.219e-08  -0.137    0.891    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.09243599)
    ## 
    ##     Null deviance: 71.050  on 299  degrees of freedom
    ## Residual deviance: 27.227  on 292  degrees of freedom
    ## AIC: 1650.3
    ## 
    ## Number of Fisher Scoring iterations: 11

``` r
vif(mod.pca.no2) # VIF <3
```

    ##                                       GVIF Df GVIF^(1/(2*Df))
    ## X.people.per.sq.km                1.486531  1        1.219234
    ## data.road.pca$x[, 1:2]            2.118656  2        1.206466
    ## data.rsd.pca$x[, 1:2]             2.202936  2        1.218289
    ## Domestic.mean.gas.consumption     1.261858  1        1.123324
    ## Non.domestic.mean.gas.consumption 1.111511  1        1.054282

``` r
# NOx glm
summary(mod.pca.nox <- glm( nox_val~X.people.per.sq.km + data.road.pca$x[,1:2] + data.rsd.pca$x[,1:2] + Domestic.mean.gas.consumption + Non.domestic.mean.gas.consumption,
                       family  = Gamma(link = "log"), data= fossilfuel_2017.omit)) #road and rsd is significant
```

    ## 
    ## Call:
    ## glm(formula = nox_val ~ X.people.per.sq.km + data.road.pca$x[, 
    ##     1:2] + data.rsd.pca$x[, 1:2] + Domestic.mean.gas.consumption + 
    ##     Non.domestic.mean.gas.consumption, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.55879  -0.23926  -0.00143   0.16741   1.96634  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        2.784e+00  2.196e-01  12.678  < 2e-16 ***
    ## X.people.per.sq.km                 1.147e-04  9.504e-06  12.066  < 2e-16 ***
    ## data.road.pca$x[, 1:2]PC1          6.961e-02  8.579e-03   8.114 1.37e-14 ***
    ## data.road.pca$x[, 1:2]PC2         -8.912e-03  9.724e-03  -0.916    0.360    
    ## data.rsd.pca$x[, 1:2]PC1          -1.317e-01  1.446e-02  -9.105  < 2e-16 ***
    ## data.rsd.pca$x[, 1:2]PC2          -1.306e-02  2.047e-02  -0.638    0.524    
    ## Domestic.mean.gas.consumption     -9.659e-06  1.541e-05  -0.627    0.531    
    ## Non.domestic.mean.gas.consumption -5.029e-09  3.833e-08  -0.131    0.896    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1310457)
    ## 
    ##     Null deviance: 94.264  on 299  degrees of freedom
    ## Residual deviance: 34.589  on 292  degrees of freedom
    ## AIC: 1909.3
    ## 
    ## Number of Fisher Scoring iterations: 12

``` r
vif(mod.pca.nox) # VIF <3
```

    ##                                       GVIF Df GVIF^(1/(2*Df))
    ## X.people.per.sq.km                1.486531  1        1.219234
    ## data.road.pca$x[, 1:2]            2.118656  2        1.206466
    ## data.rsd.pca$x[, 1:2]             2.202936  2        1.218289
    ## Domestic.mean.gas.consumption     1.261858  1        1.123324
    ## Non.domestic.mean.gas.consumption 1.111511  1        1.054282

``` r
# O3
summary(mod.pca.o3 <- glm( o3_val~X.people.per.sq.km + data.road.pca$x[,1:2] + data.rsd.pca$x[,1:2] + Domestic.mean.gas.consumption + Non.domestic.mean.gas.consumption,
                       family  = "gaussian", data= fossilfuel_2017.omit)) # road, rsd and non-domestic gas
```

    ## 
    ## Call:
    ## glm(formula = o3_val ~ X.people.per.sq.km + data.road.pca$x[, 
    ##     1:2] + data.rsd.pca$x[, 1:2] + Domestic.mean.gas.consumption + 
    ##     Non.domestic.mean.gas.consumption, family = "gaussian", data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -7.6768  -2.7207   0.1503   2.7315   9.3601  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        8.925e+00  1.951e+00   4.575 7.06e-06 ***
    ## X.people.per.sq.km                -2.922e-04  8.443e-05  -3.461 0.000619 ***
    ## data.road.pca$x[, 1:2]PC1         -2.078e-01  7.622e-02  -2.727 0.006779 ** 
    ## data.road.pca$x[, 1:2]PC2          4.271e-02  8.639e-02   0.494 0.621387    
    ## data.rsd.pca$x[, 1:2]PC1           3.273e-01  1.285e-01   2.548 0.011355 *  
    ## data.rsd.pca$x[, 1:2]PC2          -4.492e-01  1.819e-01  -2.470 0.014089 *  
    ## Domestic.mean.gas.consumption      2.954e-05  1.369e-04   0.216 0.829293    
    ## Non.domestic.mean.gas.consumption -1.455e-06  3.405e-07  -4.272 2.63e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 10.34288)
    ## 
    ##     Null deviance: 3556.1  on 299  degrees of freedom
    ## Residual deviance: 3020.1  on 292  degrees of freedom
    ## AIC: 1562.1
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
vif(mod.pca.o3) # VIF <3
```

    ##                                       GVIF Df GVIF^(1/(2*Df))
    ## X.people.per.sq.km                1.486531  1        1.219234
    ## data.road.pca$x[, 1:2]            2.118656  2        1.206466
    ## data.rsd.pca$x[, 1:2]             2.202936  2        1.218289
    ## Domestic.mean.gas.consumption     1.261858  1        1.123324
    ## Non.domestic.mean.gas.consumption 1.111511  1        1.054282

``` r
# SO2
summary(mod.pca.so2 <- glm( so2_val~X.people.per.sq.km + data.road.pca$x[,1:2] + data.rsd.pca$x[,1:2] + Domestic.mean.gas.consumption + Non.domestic.mean.gas.consumption,
                       family  = Gamma(link = "log"), data= fossilfuel_2017.omit)) # all significant
```

    ## 
    ## Call:
    ## glm(formula = so2_val ~ X.people.per.sq.km + data.road.pca$x[, 
    ##     1:2] + data.rsd.pca$x[, 1:2] + Domestic.mean.gas.consumption + 
    ##     Non.domestic.mean.gas.consumption, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.05261  -0.24711  -0.04019   0.15673   0.83290  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        7.901e-01  1.951e-01   4.050 6.58e-05 ***
    ## X.people.per.sq.km                 4.349e-05  8.444e-06   5.150 4.79e-07 ***
    ## data.road.pca$x[, 1:2]PC1          5.381e-02  7.622e-03   7.059 1.22e-11 ***
    ## data.road.pca$x[, 1:2]PC2         -7.448e-03  8.639e-03  -0.862 0.389349    
    ## data.rsd.pca$x[, 1:2]PC1          -1.034e-01  1.285e-02  -8.050 2.11e-14 ***
    ## data.rsd.pca$x[, 1:2]PC2           4.960e-02  1.819e-02   2.727 0.006785 ** 
    ## Domestic.mean.gas.consumption     -4.704e-05  1.369e-05  -3.436 0.000675 ***
    ## Non.domestic.mean.gas.consumption  1.070e-07  3.406e-08   3.141 0.001857 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1034392)
    ## 
    ##     Null deviance: 49.192  on 299  degrees of freedom
    ## Residual deviance: 30.357  on 292  degrees of freedom
    ## AIC: 346.66
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
vif(mod.pca.so2) # VIF <3
```

    ##                                       GVIF Df GVIF^(1/(2*Df))
    ## X.people.per.sq.km                1.486531  1        1.219234
    ## data.road.pca$x[, 1:2]            2.118656  2        1.206466
    ## data.rsd.pca$x[, 1:2]             2.202936  2        1.218289
    ## Domestic.mean.gas.consumption     1.261858  1        1.123324
    ## Non.domestic.mean.gas.consumption 1.111511  1        1.054282

### Calculate OR for principal component models:

OR ratio graphs for PCA models

``` r
# Calculate odds ratio (OR) for plotting 
# Format OR
NO2_pca_odds = data.frame(cbind(exp(cbind(OR = coef(mod.pca.no2), confint(mod.pca.no2))), p_value = summary(mod.pca.no2)$coefficients[,4]))
NOx_pca_odds = data.frame(cbind(exp(cbind(OR = coef(mod.pca.nox), confint(mod.pca.nox))), p_value = summary(mod.pca.nox)$coefficients[,4]))
O3_pca_odds = data.frame(cbind(exp(cbind(OR = coef(mod.pca.o3), confint(mod.pca.o3))), p_value = summary(mod.pca.o3)$coefficients[,4]))
SO2_pca_odds = data.frame(cbind(exp(cbind(OR = coef(mod.pca.so2), confint(mod.pca.so2))), p_value = summary(mod.pca.no2)$coefficients[,4]))

#now delete the rows that you dont want - Intercepts
NO2_pca_odds_clean = NO2_pca_odds[-c(1,2),]
NOx_pca_odds_clean = NOx_pca_odds[-c(1,2),]
O3_pca_odds_clean = O3_pca_odds[-c(1,2),]
SO2_pca_odds_clean = SO2_pca_odds[-c(1,2),]

#make a df just with the pollutants 
NO2_pca_onlyPoll = data.frame(NO2_pca_odds_clean)
NOx_pca_onlyPoll = data.frame(NOx_pca_odds_clean)
O3_pca_onlyPoll = data.frame(O3_pca_odds_clean)
SO2_pca_onlyPoll = data.frame(SO2_pca_odds_clean)

# Sort data for NO2
NO2_pca_onlyPoll$names = row.names(NO2_pca_onlyPoll)
NO2_pca_onlyPoll$significance = "p-value > 0.05"
NO2_pca_onlyPoll$significance[NO2_pca_onlyPoll$p_value < 0.05] <- "p-value < 0.05"
# Sort data for NOx
NOx_pca_onlyPoll$names = row.names(NOx_pca_onlyPoll)
NOx_pca_onlyPoll$significance = "p-value > 0.05"
NOx_pca_onlyPoll$significance[NOx_pca_onlyPoll$p_value < 0.05] <- "p-value < 0.05"
# Sort data for O3
O3_pca_onlyPoll$names = row.names(O3_pca_onlyPoll)
O3_pca_onlyPoll$significance = "p-value > 0.05"
O3_pca_onlyPoll$significance[O3_pca_onlyPoll$p_value < 0.05] <- "p-value < 0.05"
# Sort data for SO2
SO2_pca_onlyPoll$names = row.names(SO2_pca_onlyPoll)
SO2_pca_onlyPoll$significance = "p-value > 0.05"
SO2_pca_onlyPoll$significance[SO2_pca_onlyPoll$p_value < 0.05] <- "p-value < 0.05"

# plot for NO2 PCA
ggplot(NO2_pca_onlyPoll, aes(x=reorder(names, OR), y=OR, color= significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_bw() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("AQ odds ratios") + 
  xlab("NO2 PCA")
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-57-1.png)<!-- -->

``` r
# plot for NOx PCA
ggplot(NOx_pca_onlyPoll, aes(x=reorder(names, OR), y=OR, color= significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_bw() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("AQ odds ratios") + 
  xlab("NOx PCA")
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-57-2.png)<!-- -->

``` r
# plot for O3 PCA
ggplot(O3_pca_onlyPoll, aes(x=reorder(names, OR), y=OR, color= significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_bw() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("AQ odds ratios") + 
  xlab("O3 PCA")
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-57-3.png)<!-- -->

``` r
# plot for SO2 PCA
ggplot(SO2_pca_onlyPoll, aes(x=reorder(names, OR), y=OR, color= significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_bw() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("AQ odds ratios") + 
  xlab("SO2 PCA")
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-57-4.png)<!-- -->

#### Supplementary Table 9. Summary of the generalised linear models based on principal components of the effect of fossil fuel consumption on air pollutant levels

<br> <br> ![PC analysis on the effect of fossil fuel consumption on air
pollutant levels](fig_out_v4/GLM_PC_AirPol_sources_fig.png)

PCA-adjusted odds ratios for estimated increase of individual air
pollutant concentrations according to levels of fuel consumption from
road transport, consumption of residual fuels and gas combustion. PCs
were not estimated for gas consumption as only two variables were
included in the analysis. We fit generalized linear models of the gamma
family (log link) to NO2, NOx and SO2 data after testing assumptions of
linearity, multivariate normality and homoscedasticity. A generalized
linear model of the gaussian family was fitted to O3 data. Results were
adjusted for population density in England. Abbreviations: PC, principal
component; Var explained, percentage of variance explained; OR, odds
ratio; CI, confidence interval; NO2, nitrogen dioxide; NOx, nitrogen
oxides; O3, ozone; SO2, sulphur dioxide; NA, Not available, PC,
principal components.

<br> \#\#\# Iterative stepwise regression

#### Aim: Determine contribution of individual fuel-burning sources after variable selection

We will start by looking at NO2. Fitting the gamma identity glm (link =
“identity”) to NO2 data throws a lot of warnings and terminates with
an overflow. This could be because the identity function is not
guaranteed to map the linear predictor to a positive real number,
yielding a valid mean parameter. For this reason, we will need to use
the log link function in the Gamma family as this forces the
predictions/output to “respect the domain” by ensuring all predicted
values are positive.

``` r
# Road transport - log link is used
summary(model_NO2_gamma <- glm(data = fossilfuel_2017.omit, no2_val ~ X.people.per.sq.km + personal.buses.A.roads + personal.buses.Minor.roads + personal.diesel.cars.Motorways + personal.buses.Motorways + personal.diesel.cars.A.roads + personal.diesel.cars.Minor.roads + personal.petrol.cars.Motorway + personal.petrol.cars.A.roads + personal.petrol.cars.Minor.roads + personal.motorcycles.Motorways + personal.motorcycles.A.roads + personal.motorcycles.Minor.roads + freight.HGV.Motorways + freight.HGV.A.roads + freight.HGV.Minor.roads + freight.diesel.LGV.Motorways + freight.diesel.LGV.A.roads + freight.diesel.LGV.Minor.roads + freight.petrol.LGV.Motorways + freight.petrol.LGV.A.roads + freight.petrol.LGV.Minor.roads, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = no2_val ~ X.people.per.sq.km + personal.buses.A.roads + 
    ##     personal.buses.Minor.roads + personal.diesel.cars.Motorways + 
    ##     personal.buses.Motorways + personal.diesel.cars.A.roads + 
    ##     personal.diesel.cars.Minor.roads + personal.petrol.cars.Motorway + 
    ##     personal.petrol.cars.A.roads + personal.petrol.cars.Minor.roads + 
    ##     personal.motorcycles.Motorways + personal.motorcycles.A.roads + 
    ##     personal.motorcycles.Minor.roads + freight.HGV.Motorways + 
    ##     freight.HGV.A.roads + freight.HGV.Minor.roads + freight.diesel.LGV.Motorways + 
    ##     freight.diesel.LGV.A.roads + freight.diesel.LGV.Minor.roads + 
    ##     freight.petrol.LGV.Motorways + freight.petrol.LGV.A.roads + 
    ##     freight.petrol.LGV.Minor.roads, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.36915  -0.17130   0.02453   0.14542   1.29314  
    ## 
    ## Coefficients:
    ##                                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                       2.357e+00  4.010e-02  58.795  < 2e-16 ***
    ## X.people.per.sq.km                6.921e-05  1.286e-05   5.381 1.58e-07 ***
    ## personal.buses.A.roads           -3.024e-05  3.477e-05  -0.870  0.38528    
    ## personal.buses.Minor.roads       -1.096e-05  3.074e-05  -0.357  0.72170    
    ## personal.diesel.cars.Motorways   -5.696e-05  6.922e-05  -0.823  0.41127    
    ## personal.buses.Motorways         -2.168e-04  2.064e-04  -1.051  0.29434    
    ## personal.diesel.cars.A.roads     -7.197e-05  5.405e-05  -1.331  0.18412    
    ## personal.diesel.cars.Minor.roads  1.261e-04  1.198e-04   1.053  0.29330    
    ## personal.petrol.cars.Motorway     7.678e-05  9.567e-05   0.803  0.42294    
    ## personal.petrol.cars.A.roads      7.122e-05  4.771e-05   1.493  0.13667    
    ## personal.petrol.cars.Minor.roads -1.161e-04  9.906e-05  -1.172  0.24211    
    ## personal.motorcycles.Motorways    1.838e-03  6.785e-04   2.709  0.00717 ** 
    ## personal.motorcycles.A.roads      2.609e-04  2.925e-04   0.892  0.37316    
    ## personal.motorcycles.Minor.roads  1.371e-05  3.703e-04   0.037  0.97049    
    ## freight.HGV.Motorways             2.573e-06  4.993e-06   0.515  0.60675    
    ## freight.HGV.A.roads               1.277e-06  4.487e-06   0.285  0.77622    
    ## freight.HGV.Minor.roads          -4.999e-05  9.174e-05  -0.545  0.58628    
    ## freight.diesel.LGV.Motorways     -9.037e-06  3.581e-04  -0.025  0.97989    
    ## freight.diesel.LGV.A.roads        8.356e-05  8.065e-05   1.036  0.30105    
    ## freight.diesel.LGV.Minor.roads   -1.024e-03  3.557e-04  -2.879  0.00430 ** 
    ## freight.petrol.LGV.Motorways      1.264e-04  1.265e-02   0.010  0.99203    
    ## freight.petrol.LGV.A.roads       -2.971e-03  2.689e-03  -1.105  0.27020    
    ## freight.petrol.LGV.Minor.roads    2.718e-02  8.926e-03   3.044  0.00255 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.08175415)
    ## 
    ##     Null deviance: 71.050  on 299  degrees of freedom
    ## Residual deviance: 24.877  on 277  degrees of freedom
    ## AIC: 1652.9
    ## 
    ## Number of Fisher Scoring iterations: 8

In the dev residuals section, we see that 1Q and 3Q are not perfectly
symmetrical as we are not log transforming our data so there is some
sort of skewness. We also get an AIC of 1652. The high AIC value could
be a consequence of a large number of explanatory variables included in
the model: variable selection is required to minimise noise.

To test interdependence of factors we will assess the VIF of each road
transport predictor included in the original model. Based on that, we
will sequentially drop the predictors with the highest VIF until each
VIF is less than 5.

``` r
library(car)
vif(model_NO2_gamma)
```

    ##               X.people.per.sq.km           personal.buses.A.roads 
    ##                         4.364386                         8.658227 
    ##       personal.buses.Minor.roads   personal.diesel.cars.Motorways 
    ##                         5.278376                      1813.546278 
    ##         personal.buses.Motorways     personal.diesel.cars.A.roads 
    ##                        13.674013                       665.083797 
    ## personal.diesel.cars.Minor.roads    personal.petrol.cars.Motorway 
    ##                      2478.069220                      1739.764191 
    ##     personal.petrol.cars.A.roads personal.petrol.cars.Minor.roads 
    ##                       657.078903                      2523.720217 
    ##   personal.motorcycles.Motorways     personal.motorcycles.A.roads 
    ##                        10.635740                        12.440160 
    ## personal.motorcycles.Minor.roads            freight.HGV.Motorways 
    ##                        10.337500                        17.214787 
    ##              freight.HGV.A.roads          freight.HGV.Minor.roads 
    ##                         4.172695                        17.478741 
    ##     freight.diesel.LGV.Motorways       freight.diesel.LGV.A.roads 
    ##                     14314.830468                       551.977298 
    ##   freight.diesel.LGV.Minor.roads     freight.petrol.LGV.Motorways 
    ##                      5702.510145                     14486.032195 
    ##       freight.petrol.LGV.A.roads   freight.petrol.LGV.Minor.roads 
    ##                       589.010803                      5553.482992

``` r
summary(model_NO2_gamma.vif2 <- glm(data = fossilfuel_2017.omit, no2_val ~ X.people.per.sq.km + personal.buses.A.roads + personal.motorcycles.Motorways + personal.motorcycles.A.roads + personal.motorcycles.Minor.roads + freight.HGV.Motorways + freight.HGV.A.roads, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = no2_val ~ X.people.per.sq.km + personal.buses.A.roads + 
    ##     personal.motorcycles.Motorways + personal.motorcycles.A.roads + 
    ##     personal.motorcycles.Minor.roads + freight.HGV.Motorways + 
    ##     freight.HGV.A.roads, family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.46769  -0.23269  -0.00132   0.17465   1.55628  
    ## 
    ## Coefficients:
    ##                                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                       2.318e+00  4.271e-02  54.261  < 2e-16 ***
    ## X.people.per.sq.km                1.189e-04  1.183e-05  10.051  < 2e-16 ***
    ## personal.buses.A.roads            7.821e-05  3.017e-05   2.593  0.01001 *  
    ## personal.motorcycles.Motorways    9.132e-04  3.907e-04   2.338  0.02008 *  
    ## personal.motorcycles.A.roads     -6.207e-04  1.854e-04  -3.348  0.00092 ***
    ## personal.motorcycles.Minor.roads  1.899e-04  2.036e-04   0.933  0.35158    
    ## freight.HGV.Motorways             1.345e-06  2.288e-06   0.588  0.55709    
    ## freight.HGV.A.roads              -6.589e-06  3.367e-06  -1.957  0.05128 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1160752)
    ## 
    ##     Null deviance: 71.050  on 299  degrees of freedom
    ## Residual deviance: 34.457  on 292  degrees of freedom
    ## AIC: 1722.2
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
vif(model_NO2_gamma.vif2) # VIF <5
```

    ##               X.people.per.sq.km           personal.buses.A.roads 
    ##                         2.601107                         4.589445 
    ##   personal.motorcycles.Motorways     personal.motorcycles.A.roads 
    ##                         2.483397                         3.520010 
    ## personal.motorcycles.Minor.roads            freight.HGV.Motorways 
    ##                         2.200378                         2.546692 
    ##              freight.HGV.A.roads 
    ##                         1.654252

The model above only includes predictors displaying modest collinearity.
However, it must be noted that the AIC penalises models with multiple
variables. Therefore, we will apply a stepwise regression method to only
include those variables that have an effect on NO2 concentration in our
final model. To check that multicollinearity is accounted for in the
absence of VIF testing, we will run stepwise regression on the
VIF-corrected and original model and compare the results.

``` r
### Define full and null models and do step procedure
model_0_glm.gamma <- glm(no2_val ~ 1, data = fossilfuel_2017.omit, family = Gamma(link = "log"))
step(model_0_glm.gamma, direction = "both", scope = formula(model_NO2_gamma, test = "Chisq"))
```

    ## Start:  AIC=1931.38
    ## no2_val ~ 1
    ## 
    ##                                    Df Deviance    AIC
    ## + X.people.per.sq.km                1   40.079 1812.6
    ## + personal.buses.A.roads            1   51.658 1857.7
    ## + personal.motorcycles.A.roads      1   61.171 1894.8
    ## + personal.motorcycles.Minor.roads  1   63.913 1905.5
    ## + freight.HGV.A.roads               1   65.739 1912.7
    ## + freight.petrol.LGV.Minor.roads    1   68.848 1924.8
    ## + freight.diesel.LGV.A.roads        1   69.364 1926.8
    ## + personal.buses.Minor.roads        1   69.675 1928.0
    ## + personal.petrol.cars.Minor.roads  1   69.747 1928.3
    ## + freight.diesel.LGV.Minor.roads    1   70.108 1929.7
    ## + personal.diesel.cars.A.roads      1   70.500 1931.2
    ## <none>                                  71.050 1931.4
    ## + personal.diesel.cars.Minor.roads  1   70.716 1932.1
    ## + freight.HGV.Motorways             1   70.936 1932.9
    ## + personal.buses.Motorways          1   70.979 1933.1
    ## + freight.petrol.LGV.Motorways      1   71.007 1933.2
    ## + personal.diesel.cars.Motorways    1   71.009 1933.2
    ## + freight.diesel.LGV.Motorways      1   71.015 1933.2
    ## + personal.petrol.cars.Motorway     1   71.022 1933.3
    ## + freight.HGV.Minor.roads           1   71.027 1933.3
    ## + personal.motorcycles.Motorways    1   71.033 1933.3
    ## + personal.petrol.cars.A.roads      1   71.038 1933.3
    ## + freight.petrol.LGV.A.roads        1   71.048 1933.4
    ## 
    ## Step:  AIC=1756.48
    ## no2_val ~ X.people.per.sq.km
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.motorcycles.Motorways    1   37.618 1738.2
    ## + freight.diesel.LGV.Motorways      1   37.697 1738.9
    ## + personal.petrol.cars.Motorway     1   37.714 1739.0
    ## + freight.petrol.LGV.Motorways      1   37.727 1739.1
    ## + personal.diesel.cars.Motorways    1   37.822 1739.9
    ## + freight.HGV.Motorways             1   38.166 1742.7
    ## + personal.petrol.cars.Minor.roads  1   38.228 1743.2
    ## + personal.buses.Motorways          1   38.426 1744.9
    ## + personal.buses.Minor.roads        1   38.614 1746.4
    ## + freight.petrol.LGV.Minor.roads    1   38.822 1748.1
    ## + personal.diesel.cars.Minor.roads  1   39.103 1750.4
    ## + freight.HGV.A.roads               1   39.160 1750.9
    ## + freight.diesel.LGV.Minor.roads    1   39.335 1752.3
    ## + personal.motorcycles.A.roads      1   39.498 1753.7
    ## + personal.motorcycles.Minor.roads  1   39.553 1754.2
    ## + personal.buses.A.roads            1   39.658 1755.0
    ## + freight.HGV.Minor.roads           1   39.701 1755.4
    ## + freight.diesel.LGV.A.roads        1   39.766 1755.9
    ## <none>                                  40.079 1756.5
    ## + personal.petrol.cars.A.roads      1   39.941 1757.3
    ## + freight.petrol.LGV.A.roads        1   40.054 1758.3
    ## + personal.diesel.cars.A.roads      1   40.060 1758.3
    ## - X.people.per.sq.km                1   71.050 2009.5
    ## 
    ## Step:  AIC=1739.05
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.buses.Minor.roads        1   36.358 1730.1
    ## + personal.petrol.cars.Minor.roads  1   36.505 1731.4
    ## + personal.motorcycles.A.roads      1   36.679 1732.9
    ## + freight.HGV.A.roads               1   36.782 1733.8
    ## + freight.petrol.LGV.Minor.roads    1   36.920 1735.0
    ## + freight.diesel.LGV.A.roads        1   37.164 1737.1
    ## + personal.diesel.cars.Minor.roads  1   37.164 1737.1
    ## + freight.diesel.LGV.Minor.roads    1   37.284 1738.2
    ## <none>                                  37.618 1739.0
    ## + personal.buses.A.roads            1   37.442 1739.5
    ## + personal.buses.Motorways          1   37.454 1739.6
    ## + freight.HGV.Minor.roads           1   37.455 1739.6
    ## + personal.motorcycles.Minor.roads  1   37.460 1739.7
    ## + freight.diesel.LGV.Motorways      1   37.481 1739.9
    ## + freight.petrol.LGV.Motorways      1   37.496 1740.0
    ## + freight.HGV.Motorways             1   37.499 1740.0
    ## + personal.diesel.cars.A.roads      1   37.505 1740.1
    ## + personal.petrol.cars.Motorway     1   37.509 1740.1
    ## + freight.petrol.LGV.A.roads        1   37.513 1740.2
    ## + personal.diesel.cars.Motorways    1   37.556 1740.5
    ## + personal.petrol.cars.A.roads      1   37.596 1740.9
    ## - personal.motorcycles.Motorways    1   40.079 1758.4
    ## - X.people.per.sq.km                1   71.033 2026.3
    ## 
    ## Step:  AIC=1730.62
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Minor.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + freight.diesel.LGV.A.roads        1   34.519 1716.3
    ## + freight.HGV.A.roads               1   34.846 1719.2
    ## + personal.diesel.cars.A.roads      1   34.992 1720.5
    ## + freight.petrol.LGV.A.roads        1   35.147 1721.8
    ## + personal.motorcycles.A.roads      1   35.245 1722.7
    ## + freight.HGV.Minor.roads           1   35.487 1724.9
    ## + personal.petrol.cars.A.roads      1   35.787 1727.5
    ## + personal.buses.Motorways          1   35.865 1728.2
    ## + freight.diesel.LGV.Minor.roads    1   35.899 1728.5
    ## + personal.diesel.cars.Minor.roads  1   36.090 1730.2
    ## <none>                                  36.358 1730.6
    ## + personal.motorcycles.Minor.roads  1   36.195 1731.2
    ## + freight.petrol.LGV.Minor.roads    1   36.303 1732.1
    ## + personal.buses.A.roads            1   36.311 1732.2
    ## + personal.petrol.cars.Minor.roads  1   36.327 1732.3
    ## + personal.diesel.cars.Motorways    1   36.353 1732.6
    ## + freight.diesel.LGV.Motorways      1   36.354 1732.6
    ## + freight.HGV.Motorways             1   36.355 1732.6
    ## + freight.petrol.LGV.Motorways      1   36.356 1732.6
    ## + personal.petrol.cars.Motorway     1   36.358 1732.6
    ## - personal.buses.Minor.roads        1   37.618 1739.8
    ## - personal.motorcycles.Motorways    1   38.614 1748.7
    ## - X.people.per.sq.km                1   69.672 2025.0
    ## 
    ## Step:  AIC=1716.74
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Minor.roads + freight.diesel.LGV.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + freight.petrol.LGV.A.roads        1   32.099 1696.9
    ## + personal.petrol.cars.A.roads      1   32.555 1701.0
    ## + personal.buses.Motorways          1   34.063 1714.6
    ## + personal.petrol.cars.Minor.roads  1   34.160 1715.5
    ## + personal.diesel.cars.A.roads      1   34.186 1715.7
    ## <none>                                  34.519 1716.7
    ## + personal.buses.A.roads            1   34.298 1716.8
    ## + freight.HGV.Minor.roads           1   34.315 1716.9
    ## + freight.petrol.LGV.Minor.roads    1   34.383 1717.5
    ## + personal.motorcycles.Minor.roads  1   34.428 1717.9
    ## + freight.HGV.A.roads               1   34.453 1718.2
    ## + freight.HGV.Motorways             1   34.498 1718.5
    ## + personal.motorcycles.A.roads      1   34.506 1718.6
    ## + personal.diesel.cars.Motorways    1   34.511 1718.7
    ## + freight.diesel.LGV.Motorways      1   34.515 1718.7
    ## + freight.petrol.LGV.Motorways      1   34.516 1718.7
    ## + personal.diesel.cars.Minor.roads  1   34.517 1718.7
    ## + freight.diesel.LGV.Minor.roads    1   34.518 1718.7
    ## + personal.petrol.cars.Motorway     1   34.519 1718.7
    ## - freight.diesel.LGV.A.roads        1   36.358 1731.4
    ## - personal.motorcycles.Motorways    1   36.917 1736.4
    ## - personal.buses.Minor.roads        1   37.164 1738.7
    ## - X.people.per.sq.km                1   65.233 1992.6
    ## 
    ## Step:  AIC=1696.54
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Minor.roads + freight.diesel.LGV.A.roads + 
    ##     freight.petrol.LGV.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.buses.Motorways          1   31.352 1691.3
    ## + freight.HGV.Minor.roads           1   31.558 1693.3
    ## <none>                                  32.099 1696.5
    ## + personal.petrol.cars.A.roads      1   31.898 1696.6
    ## + freight.diesel.LGV.Minor.roads    1   31.945 1697.0
    ## - personal.buses.Minor.roads        1   32.367 1697.1
    ## + personal.diesel.cars.Minor.roads  1   31.970 1697.3
    ## + personal.diesel.cars.Motorways    1   31.972 1697.3
    ## + personal.buses.A.roads            1   31.996 1697.5
    ## + personal.petrol.cars.Motorway     1   32.028 1697.8
    ## + personal.motorcycles.A.roads      1   32.047 1698.0
    ## + freight.diesel.LGV.Motorways      1   32.048 1698.0
    ## + freight.petrol.LGV.Motorways      1   32.049 1698.1
    ## + personal.motorcycles.Minor.roads  1   32.062 1698.2
    ## + freight.petrol.LGV.Minor.roads    1   32.088 1698.4
    ## + freight.HGV.A.roads               1   32.095 1698.5
    ## + personal.petrol.cars.Minor.roads  1   32.095 1698.5
    ## + freight.HGV.Motorways             1   32.099 1698.5
    ## + personal.diesel.cars.A.roads      1   32.099 1698.5
    ## - personal.motorcycles.Motorways    1   33.577 1708.9
    ## - freight.petrol.LGV.A.roads        1   34.519 1718.0
    ## - freight.diesel.LGV.A.roads        1   35.147 1724.1
    ## - X.people.per.sq.km                1   39.277 1764.1
    ## 
    ## Step:  AIC=1691.35
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Minor.roads + freight.diesel.LGV.A.roads + 
    ##     freight.petrol.LGV.A.roads + personal.buses.Motorways
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.petrol.cars.A.roads      1   30.919 1689.0
    ## + freight.HGV.Minor.roads           1   30.944 1689.2
    ## + freight.HGV.Motorways             1   31.105 1690.9
    ## <none>                                  31.352 1691.3
    ## + personal.petrol.cars.Motorway     1   31.197 1691.8
    ## + freight.diesel.LGV.Minor.roads    1   31.233 1692.2
    ## + freight.petrol.LGV.Motorways      1   31.256 1692.4
    ## + freight.diesel.LGV.Motorways      1   31.258 1692.4
    ## + personal.diesel.cars.Motorways    1   31.263 1692.5
    ## + personal.diesel.cars.A.roads      1   31.280 1692.6
    ## + personal.diesel.cars.Minor.roads  1   31.286 1692.7
    ## + personal.motorcycles.A.roads      1   31.294 1692.8
    ## + personal.buses.A.roads            1   31.298 1692.8
    ## - personal.buses.Minor.roads        1   31.720 1693.0
    ## + personal.motorcycles.Minor.roads  1   31.325 1693.1
    ## + personal.petrol.cars.Minor.roads  1   31.332 1693.2
    ## + freight.HGV.A.roads               1   31.339 1693.2
    ## + freight.petrol.LGV.Minor.roads    1   31.347 1693.3
    ## - personal.buses.Motorways          1   32.099 1696.9
    ## - personal.motorcycles.Motorways    1   32.982 1705.7
    ## - freight.petrol.LGV.A.roads        1   34.063 1716.6
    ## - freight.diesel.LGV.A.roads        1   34.713 1723.1
    ## - X.people.per.sq.km                1   37.756 1753.7
    ## 
    ## Step:  AIC=1689.11
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Minor.roads + freight.diesel.LGV.A.roads + 
    ##     freight.petrol.LGV.A.roads + personal.buses.Motorways + personal.petrol.cars.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.diesel.cars.A.roads      1   30.227 1684.2
    ## + freight.HGV.Minor.roads           1   30.510 1687.0
    ## - personal.buses.Minor.roads        1   31.011 1688.0
    ## + freight.HGV.Motorways             1   30.619 1688.1
    ## <none>                                  30.919 1689.1
    ## + personal.diesel.cars.Minor.roads  1   30.754 1689.5
    ## + freight.diesel.LGV.Minor.roads    1   30.788 1689.8
    ## + freight.petrol.LGV.Motorways      1   30.808 1690.0
    ## + freight.diesel.LGV.Motorways      1   30.811 1690.0
    ## + personal.petrol.cars.Motorway     1   30.812 1690.0
    ## + personal.diesel.cars.Motorways    1   30.869 1690.6
    ## + freight.HGV.A.roads               1   30.895 1690.9
    ## + personal.buses.A.roads            1   30.897 1690.9
    ## + personal.motorcycles.Minor.roads  1   30.905 1691.0
    ## + freight.petrol.LGV.Minor.roads    1   30.909 1691.0
    ## + personal.motorcycles.A.roads      1   30.914 1691.1
    ## + personal.petrol.cars.Minor.roads  1   30.918 1691.1
    ## - personal.petrol.cars.A.roads      1   31.352 1691.4
    ## - freight.petrol.LGV.A.roads        1   31.464 1692.5
    ## - personal.buses.Motorways          1   31.898 1696.8
    ## - personal.motorcycles.Motorways    1   32.769 1705.5
    ## - freight.diesel.LGV.A.roads        1   33.094 1708.7
    ## - X.people.per.sq.km                1   37.631 1753.9
    ## 
    ## Step:  AIC=1684.2
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Minor.roads + freight.diesel.LGV.A.roads + 
    ##     freight.petrol.LGV.A.roads + personal.buses.Motorways + personal.petrol.cars.A.roads + 
    ##     personal.diesel.cars.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.motorcycles.A.roads      1   29.588 1680.0
    ## + freight.HGV.Minor.roads           1   29.620 1680.3
    ## - freight.petrol.LGV.A.roads        1   30.228 1682.2
    ## - personal.buses.Minor.roads        1   30.236 1682.3
    ## - freight.diesel.LGV.A.roads        1   30.259 1682.5
    ## + personal.diesel.cars.Minor.roads  1   29.925 1683.3
    ## + freight.HGV.Motorways             1   29.973 1683.7
    ## + freight.diesel.LGV.Minor.roads    1   30.013 1684.1
    ## <none>                                  30.227 1684.2
    ## + personal.buses.A.roads            1   30.057 1684.5
    ## + personal.petrol.cars.Motorway     1   30.138 1685.3
    ## + freight.petrol.LGV.Motorways      1   30.140 1685.4
    ## + freight.diesel.LGV.Motorways      1   30.143 1685.4
    ## + personal.petrol.cars.Minor.roads  1   30.160 1685.5
    ## + freight.petrol.LGV.Minor.roads    1   30.173 1685.7
    ## + personal.diesel.cars.Motorways    1   30.183 1685.8
    ## + personal.motorcycles.Minor.roads  1   30.212 1686.0
    ## + freight.HGV.A.roads               1   30.212 1686.0
    ## - personal.buses.Motorways          1   30.886 1688.6
    ## - personal.diesel.cars.A.roads      1   30.919 1688.9
    ## - personal.petrol.cars.A.roads      1   31.280 1692.4
    ## - personal.motorcycles.Motorways    1   31.730 1696.8
    ## - X.people.per.sq.km                1   37.523 1753.1
    ## 
    ## Step:  AIC=1679.69
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Minor.roads + freight.diesel.LGV.A.roads + 
    ##     freight.petrol.LGV.A.roads + personal.buses.Motorways + personal.petrol.cars.A.roads + 
    ##     personal.diesel.cars.A.roads + personal.motorcycles.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + freight.HGV.Minor.roads           1   28.886 1674.5
    ## + personal.diesel.cars.Minor.roads  1   29.028 1675.9
    ## + freight.diesel.LGV.Minor.roads    1   29.168 1677.4
    ## - personal.buses.Minor.roads        1   29.605 1677.9
    ## + freight.HGV.Motorways             1   29.326 1679.0
    ## + personal.petrol.cars.Minor.roads  1   29.373 1679.5
    ## - freight.diesel.LGV.A.roads        1   29.772 1679.6
    ## <none>                                  29.588 1679.7
    ## + freight.petrol.LGV.Minor.roads    1   29.396 1679.7
    ## + personal.motorcycles.Minor.roads  1   29.443 1680.2
    ## + freight.HGV.A.roads               1   29.460 1680.4
    ## + freight.petrol.LGV.Motorways      1   29.473 1680.5
    ## + freight.diesel.LGV.Motorways      1   29.475 1680.5
    ## + personal.petrol.cars.Motorway     1   29.485 1680.6
    ## - freight.petrol.LGV.A.roads        1   29.899 1680.9
    ## + personal.diesel.cars.Motorways    1   29.528 1681.1
    ## + personal.buses.A.roads            1   29.588 1681.7
    ## - personal.buses.Motorways          1   30.067 1682.6
    ## - personal.motorcycles.A.roads      1   30.227 1684.2
    ## - personal.motorcycles.Motorways    1   30.721 1689.3
    ## - personal.diesel.cars.A.roads      1   30.914 1691.3
    ## - personal.petrol.cars.A.roads      1   31.219 1694.4
    ## - X.people.per.sq.km                1   33.659 1719.4
    ## 
    ## Step:  AIC=1674.37
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Minor.roads + freight.diesel.LGV.A.roads + 
    ##     freight.petrol.LGV.A.roads + personal.buses.Motorways + personal.petrol.cars.A.roads + 
    ##     personal.diesel.cars.A.roads + personal.motorcycles.A.roads + 
    ##     freight.HGV.Minor.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + freight.HGV.Motorways             1   28.598 1673.3
    ## <none>                                  28.886 1674.4
    ## + freight.petrol.LGV.Minor.roads    1   28.710 1674.5
    ## + freight.petrol.LGV.Motorways      1   28.714 1674.5
    ## + freight.diesel.LGV.Motorways      1   28.715 1674.5
    ## + personal.petrol.cars.Motorway     1   28.740 1674.8
    ## + freight.HGV.A.roads               1   28.747 1674.9
    ## - freight.diesel.LGV.A.roads        1   29.150 1675.2
    ## + personal.petrol.cars.Minor.roads  1   28.777 1675.2
    ## + personal.diesel.cars.Motorways    1   28.788 1675.3
    ## - personal.buses.Motorways          1   29.194 1675.6
    ## + personal.buses.A.roads            1   28.834 1675.8
    ## - personal.buses.Minor.roads        1   29.233 1676.0
    ## - freight.petrol.LGV.A.roads        1   29.250 1676.2
    ## + personal.diesel.cars.Minor.roads  1   28.872 1676.2
    ## + freight.diesel.LGV.Minor.roads    1   28.881 1676.3
    ## + personal.motorcycles.Minor.roads  1   28.885 1676.3
    ## - freight.HGV.Minor.roads           1   29.588 1679.8
    ## - personal.motorcycles.A.roads      1   29.620 1680.1
    ## - personal.motorcycles.Motorways    1   29.829 1682.3
    ## - personal.diesel.cars.A.roads      1   30.495 1689.3
    ## - personal.petrol.cars.A.roads      1   30.804 1692.6
    ## - X.people.per.sq.km                1   32.413 1709.5
    ## 
    ## Step:  AIC=1673.31
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Minor.roads + freight.diesel.LGV.A.roads + 
    ##     freight.petrol.LGV.A.roads + personal.buses.Motorways + personal.petrol.cars.A.roads + 
    ##     personal.diesel.cars.A.roads + personal.motorcycles.A.roads + 
    ##     freight.HGV.Minor.roads + freight.HGV.Motorways
    ## 
    ##                                    Df Deviance    AIC
    ## + freight.petrol.LGV.Minor.roads    1   28.335 1672.5
    ## <none>                                  28.598 1673.3
    ## + personal.petrol.cars.Minor.roads  1   28.428 1673.5
    ## - freight.diesel.LGV.A.roads        1   28.871 1674.2
    ## - freight.HGV.Motorways             1   28.886 1674.4
    ## + freight.HGV.A.roads               1   28.530 1674.6
    ## - personal.buses.Minor.roads        1   28.915 1674.7
    ## + personal.buses.A.roads            1   28.550 1674.8
    ## + freight.diesel.LGV.Minor.roads    1   28.564 1674.9
    ## + personal.motorcycles.Minor.roads  1   28.578 1675.1
    ## + personal.diesel.cars.Motorways    1   28.583 1675.2
    ## + freight.petrol.LGV.Motorways      1   28.585 1675.2
    ## + freight.diesel.LGV.Motorways      1   28.586 1675.2
    ## + personal.petrol.cars.Motorway     1   28.598 1675.3
    ## + personal.diesel.cars.Minor.roads  1   28.598 1675.3
    ## - freight.petrol.LGV.A.roads        1   28.988 1675.5
    ## - personal.buses.Motorways          1   29.161 1677.3
    ## - freight.HGV.Minor.roads           1   29.326 1679.1
    ## - personal.motorcycles.A.roads      1   29.347 1679.3
    ## - personal.motorcycles.Motorways    1   29.613 1682.1
    ## - personal.diesel.cars.A.roads      1   30.179 1688.2
    ## - personal.petrol.cars.A.roads      1   30.522 1691.8
    ## - X.people.per.sq.km                1   32.267 1710.4
    ## 
    ## Step:  AIC=1672.49
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Minor.roads + freight.diesel.LGV.A.roads + 
    ##     freight.petrol.LGV.A.roads + personal.buses.Motorways + personal.petrol.cars.A.roads + 
    ##     personal.diesel.cars.A.roads + personal.motorcycles.A.roads + 
    ##     freight.HGV.Minor.roads + freight.HGV.Motorways + freight.petrol.LGV.Minor.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + freight.diesel.LGV.Minor.roads    1   25.168 1640.9
    ## + personal.diesel.cars.Minor.roads  1   27.916 1670.0
    ## - personal.buses.Minor.roads        1   28.447 1671.7
    ## <none>                                  28.335 1672.5
    ## - freight.diesel.LGV.A.roads        1   28.524 1672.5
    ## - freight.petrol.LGV.Minor.roads    1   28.598 1673.3
    ## + freight.HGV.A.roads               1   28.224 1673.3
    ## + personal.motorcycles.Minor.roads  1   28.224 1673.3
    ## - freight.petrol.LGV.A.roads        1   28.621 1673.5
    ## + personal.buses.A.roads            1   28.293 1674.0
    ## + personal.diesel.cars.Motorways    1   28.322 1674.4
    ## + freight.petrol.LGV.Motorways      1   28.329 1674.4
    ## + freight.diesel.LGV.Motorways      1   28.329 1674.4
    ## - freight.HGV.Motorways             1   28.710 1674.5
    ## + personal.petrol.cars.Minor.roads  1   28.334 1674.5
    ## + personal.petrol.cars.Motorway     1   28.335 1674.5
    ## - personal.motorcycles.A.roads      1   28.778 1675.2
    ## - personal.buses.Motorways          1   28.923 1676.7
    ## - freight.HGV.Minor.roads           1   29.162 1679.3
    ## - personal.motorcycles.Motorways    1   29.216 1679.8
    ## - personal.diesel.cars.A.roads      1   29.498 1682.8
    ## - personal.petrol.cars.A.roads      1   29.744 1685.4
    ## - X.people.per.sq.km                1   31.834 1707.6
    ## 
    ## Step:  AIC=1638.41
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Minor.roads + freight.diesel.LGV.A.roads + 
    ##     freight.petrol.LGV.A.roads + personal.buses.Motorways + personal.petrol.cars.A.roads + 
    ##     personal.diesel.cars.A.roads + personal.motorcycles.A.roads + 
    ##     freight.HGV.Minor.roads + freight.HGV.Motorways + freight.petrol.LGV.Minor.roads + 
    ##     freight.diesel.LGV.Minor.roads
    ## 
    ##                                    Df Deviance    AIC
    ## - personal.buses.Minor.roads        1   25.203 1636.8
    ## - freight.HGV.Motorways             1   25.209 1636.9
    ## - personal.motorcycles.A.roads      1   25.272 1637.7
    ## - freight.HGV.Minor.roads           1   25.289 1637.9
    ## - freight.petrol.LGV.A.roads        1   25.290 1637.9
    ## - personal.diesel.cars.A.roads      1   25.304 1638.1
    ## - freight.diesel.LGV.A.roads        1   25.304 1638.1
    ## - personal.petrol.cars.A.roads      1   25.313 1638.2
    ## <none>                                  25.169 1638.4
    ## - personal.buses.Motorways          1   25.387 1639.1
    ## + personal.petrol.cars.Minor.roads  1   25.106 1639.6
    ## + personal.diesel.cars.Minor.roads  1   25.122 1639.8
    ## + freight.diesel.LGV.Motorways      1   25.127 1639.9
    ## + freight.petrol.LGV.Motorways      1   25.131 1640.0
    ## + personal.diesel.cars.Motorways    1   25.131 1640.0
    ## + personal.petrol.cars.Motorway     1   25.142 1640.1
    ## + personal.buses.A.roads            1   25.148 1640.2
    ## + freight.HGV.A.roads               1   25.151 1640.2
    ## + personal.motorcycles.Minor.roads  1   25.168 1640.4
    ## - personal.motorcycles.Motorways    1   25.827 1644.5
    ## - X.people.per.sq.km                1   27.697 1667.5
    ## - freight.diesel.LGV.Minor.roads    1   28.335 1675.3
    ## - freight.petrol.LGV.Minor.roads    1   28.564 1678.1
    ## 
    ## Step:  AIC=1636.84
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     freight.diesel.LGV.A.roads + freight.petrol.LGV.A.roads + 
    ##     personal.buses.Motorways + personal.petrol.cars.A.roads + 
    ##     personal.diesel.cars.A.roads + personal.motorcycles.A.roads + 
    ##     freight.HGV.Minor.roads + freight.HGV.Motorways + freight.petrol.LGV.Minor.roads + 
    ##     freight.diesel.LGV.Minor.roads
    ## 
    ##                                    Df Deviance    AIC
    ## - freight.HGV.Motorways             1   25.242 1635.3
    ## - freight.petrol.LGV.A.roads        1   25.336 1636.5
    ## - freight.HGV.Minor.roads           1   25.336 1636.5
    ## - personal.motorcycles.A.roads      1   25.341 1636.5
    ## - freight.diesel.LGV.A.roads        1   25.354 1636.7
    ## <none>                                  25.203 1636.8
    ## - personal.diesel.cars.A.roads      1   25.367 1636.9
    ## - personal.petrol.cars.A.roads      1   25.369 1636.9
    ## - personal.buses.Motorways          1   25.419 1637.5
    ## + personal.petrol.cars.Minor.roads  1   25.141 1638.1
    ## + personal.diesel.cars.Minor.roads  1   25.159 1638.3
    ## + freight.diesel.LGV.Motorways      1   25.165 1638.4
    ## + freight.petrol.LGV.Motorways      1   25.168 1638.4
    ## + personal.buses.Minor.roads        1   25.169 1638.4
    ## + personal.diesel.cars.Motorways    1   25.171 1638.4
    ## + personal.buses.A.roads            1   25.179 1638.5
    ## + personal.petrol.cars.Motorway     1   25.181 1638.6
    ## + freight.HGV.A.roads               1   25.182 1638.6
    ## + personal.motorcycles.Minor.roads  1   25.196 1638.8
    ## - personal.motorcycles.Motorways    1   25.887 1643.3
    ## - X.people.per.sq.km                1   27.717 1665.8
    ## - freight.diesel.LGV.Minor.roads    1   28.447 1674.8
    ## - freight.petrol.LGV.Minor.roads    1   28.791 1679.1
    ## 
    ## Step:  AIC=1635.3
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     freight.diesel.LGV.A.roads + freight.petrol.LGV.A.roads + 
    ##     personal.buses.Motorways + personal.petrol.cars.A.roads + 
    ##     personal.diesel.cars.A.roads + personal.motorcycles.A.roads + 
    ##     freight.HGV.Minor.roads + freight.petrol.LGV.Minor.roads + 
    ##     freight.diesel.LGV.Minor.roads
    ## 
    ##                                    Df Deviance    AIC
    ## - freight.HGV.Minor.roads           1   25.352 1634.7
    ## - freight.petrol.LGV.A.roads        1   25.369 1634.9
    ## - personal.motorcycles.A.roads      1   25.376 1635.0
    ## - freight.diesel.LGV.A.roads        1   25.392 1635.2
    ## - personal.petrol.cars.A.roads      1   25.398 1635.2
    ## - personal.diesel.cars.A.roads      1   25.400 1635.2
    ## <none>                                  25.242 1635.3
    ## - personal.buses.Motorways          1   25.422 1635.5
    ## + personal.petrol.cars.Minor.roads  1   25.174 1636.5
    ## + personal.diesel.cars.Minor.roads  1   25.194 1636.7
    ## + freight.HGV.Motorways             1   25.203 1636.8
    ## + personal.buses.Minor.roads        1   25.209 1636.9
    ## + freight.HGV.A.roads               1   25.210 1636.9
    ## + personal.buses.A.roads            1   25.217 1637.0
    ## + freight.petrol.LGV.Motorways      1   25.236 1637.2
    ## + personal.motorcycles.Minor.roads  1   25.236 1637.2
    ## + freight.diesel.LGV.Motorways      1   25.237 1637.2
    ## + personal.petrol.cars.Motorway     1   25.240 1637.3
    ## + personal.diesel.cars.Motorways    1   25.241 1637.3
    ## - personal.motorcycles.Motorways    1   25.908 1641.5
    ## - X.people.per.sq.km                1   27.719 1663.9
    ## - freight.diesel.LGV.Minor.roads    1   28.869 1678.2
    ## - freight.petrol.LGV.Minor.roads    1   29.171 1681.9
    ## 
    ## Step:  AIC=1634.62
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     freight.diesel.LGV.A.roads + freight.petrol.LGV.A.roads + 
    ##     personal.buses.Motorways + personal.petrol.cars.A.roads + 
    ##     personal.diesel.cars.A.roads + personal.motorcycles.A.roads + 
    ##     freight.petrol.LGV.Minor.roads + freight.diesel.LGV.Minor.roads
    ## 
    ##                                    Df Deviance    AIC
    ## - freight.petrol.LGV.A.roads        1   25.492 1634.4
    ## - freight.diesel.LGV.A.roads        1   25.513 1634.6
    ## <none>                                  25.352 1634.6
    ## - personal.diesel.cars.A.roads      1   25.513 1634.6
    ## - personal.petrol.cars.A.roads      1   25.514 1634.6
    ## - personal.motorcycles.A.roads      1   25.536 1634.9
    ## + personal.petrol.cars.Minor.roads  1   25.241 1635.2
    ## + freight.HGV.Minor.roads           1   25.242 1635.3
    ## - personal.buses.Motorways          1   25.570 1635.3
    ## + personal.diesel.cars.Minor.roads  1   25.272 1635.6
    ## + personal.buses.Minor.roads        1   25.306 1636.1
    ## + personal.buses.A.roads            1   25.307 1636.1
    ## + personal.motorcycles.Minor.roads  1   25.313 1636.1
    ## + freight.HGV.Motorways             1   25.336 1636.4
    ## + freight.HGV.A.roads               1   25.337 1636.4
    ## + personal.diesel.cars.Motorways    1   25.350 1636.6
    ## + personal.petrol.cars.Motorway     1   25.351 1636.6
    ## + freight.diesel.LGV.Motorways      1   25.352 1636.6
    ## + freight.petrol.LGV.Motorways      1   25.352 1636.6
    ## - personal.motorcycles.Motorways    1   26.120 1642.2
    ## - X.people.per.sq.km                1   27.956 1664.9
    ## - freight.petrol.LGV.Minor.roads    1   29.422 1683.1
    ## - freight.diesel.LGV.Minor.roads    1   29.541 1684.6
    ## 
    ## Step:  AIC=1634.3
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     freight.diesel.LGV.A.roads + personal.buses.Motorways + personal.petrol.cars.A.roads + 
    ##     personal.diesel.cars.A.roads + personal.motorcycles.A.roads + 
    ##     freight.petrol.LGV.Minor.roads + freight.diesel.LGV.Minor.roads
    ## 
    ##                                    Df Deviance    AIC
    ## - freight.diesel.LGV.A.roads        1   25.515 1632.6
    ## - personal.petrol.cars.A.roads      1   25.517 1632.6
    ## - personal.diesel.cars.A.roads      1   25.520 1632.7
    ## - personal.motorcycles.A.roads      1   25.545 1633.0
    ## <none>                                  25.492 1634.3
    ## + freight.petrol.LGV.A.roads        1   25.352 1634.6
    ## + freight.HGV.Minor.roads           1   25.369 1634.8
    ## + personal.buses.A.roads            1   25.383 1635.0
    ## + personal.buses.Minor.roads        1   25.433 1635.6
    ## - personal.buses.Motorways          1   25.763 1635.7
    ## + personal.petrol.cars.Minor.roads  1   25.448 1635.8
    ## + personal.motorcycles.Minor.roads  1   25.456 1635.8
    ## + personal.diesel.cars.Minor.roads  1   25.464 1636.0
    ## + freight.HGV.Motorways             1   25.481 1636.2
    ## + freight.HGV.A.roads               1   25.486 1636.2
    ## + personal.diesel.cars.Motorways    1   25.489 1636.3
    ## + freight.diesel.LGV.Motorways      1   25.490 1636.3
    ## + freight.petrol.LGV.Motorways      1   25.491 1636.3
    ## + personal.petrol.cars.Motorway     1   25.491 1636.3
    ## - personal.motorcycles.Motorways    1   26.361 1643.1
    ## - X.people.per.sq.km                1   27.967 1663.1
    ## - freight.petrol.LGV.Minor.roads    1   29.796 1685.8
    ## - freight.diesel.LGV.Minor.roads    1   29.890 1687.0
    ## 
    ## Step:  AIC=1632.57
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Motorways + personal.petrol.cars.A.roads + 
    ##     personal.diesel.cars.A.roads + personal.motorcycles.A.roads + 
    ##     freight.petrol.LGV.Minor.roads + freight.diesel.LGV.Minor.roads
    ## 
    ##                                    Df Deviance    AIC
    ## - personal.diesel.cars.A.roads      1   25.526 1630.7
    ## - personal.petrol.cars.A.roads      1   25.534 1630.8
    ## - personal.motorcycles.A.roads      1   25.573 1631.3
    ## <none>                                  25.515 1632.6
    ## + personal.buses.A.roads            1   25.392 1633.0
    ## + freight.HGV.Minor.roads           1   25.395 1633.1
    ## + personal.buses.Minor.roads        1   25.450 1633.8
    ## + personal.petrol.cars.Minor.roads  1   25.454 1633.8
    ## + personal.diesel.cars.Minor.roads  1   25.471 1634.0
    ## + personal.motorcycles.Minor.roads  1   25.473 1634.0
    ## + freight.diesel.LGV.A.roads        1   25.492 1634.3
    ## + freight.HGV.A.roads               1   25.498 1634.4
    ## + freight.HGV.Motorways             1   25.499 1634.4
    ## + personal.diesel.cars.Motorways    1   25.510 1634.5
    ## + freight.petrol.LGV.A.roads        1   25.513 1634.5
    ## + personal.petrol.cars.Motorway     1   25.513 1634.5
    ## + freight.diesel.LGV.Motorways      1   25.513 1634.6
    ## + freight.petrol.LGV.Motorways      1   25.514 1634.6
    ## - personal.buses.Motorways          1   25.870 1635.0
    ## - personal.motorcycles.Motorways    1   26.491 1642.8
    ## - X.people.per.sq.km                1   27.988 1661.4
    ## - freight.petrol.LGV.Minor.roads    1   29.960 1686.0
    ## - freight.diesel.LGV.Minor.roads    1   30.071 1687.4
    ## 
    ## Step:  AIC=1630.7
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Motorways + personal.petrol.cars.A.roads + 
    ##     personal.motorcycles.A.roads + freight.petrol.LGV.Minor.roads + 
    ##     freight.diesel.LGV.Minor.roads
    ## 
    ##                                    Df Deviance    AIC
    ## - personal.petrol.cars.A.roads      1   25.570 1629.2
    ## - personal.motorcycles.A.roads      1   25.574 1629.3
    ## <none>                                  25.526 1630.7
    ## + freight.HGV.Minor.roads           1   25.410 1631.2
    ## + personal.buses.A.roads            1   25.424 1631.4
    ## + personal.buses.Minor.roads        1   25.454 1631.8
    ## + personal.petrol.cars.Minor.roads  1   25.471 1632.0
    ## + personal.motorcycles.Minor.roads  1   25.477 1632.1
    ## + personal.diesel.cars.Minor.roads  1   25.486 1632.2
    ## + freight.HGV.Motorways             1   25.512 1632.5
    ## + personal.diesel.cars.A.roads      1   25.515 1632.6
    ## + freight.HGV.A.roads               1   25.518 1632.6
    ## + freight.diesel.LGV.A.roads        1   25.520 1632.6
    ## + personal.diesel.cars.Motorways    1   25.520 1632.6
    ## + personal.petrol.cars.Motorway     1   25.524 1632.7
    ## + freight.diesel.LGV.Motorways      1   25.524 1632.7
    ## + freight.petrol.LGV.A.roads        1   25.525 1632.7
    ## + freight.petrol.LGV.Motorways      1   25.525 1632.7
    ## - personal.buses.Motorways          1   25.879 1633.1
    ## - personal.motorcycles.Motorways    1   26.502 1640.9
    ## - X.people.per.sq.km                1   28.360 1664.2
    ## - freight.diesel.LGV.Minor.roads    1   35.052 1747.8
    ## - freight.petrol.LGV.Minor.roads    1   35.499 1753.5
    ## 
    ## Step:  AIC=1629.22
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Motorways + personal.motorcycles.A.roads + 
    ##     freight.petrol.LGV.Minor.roads + freight.diesel.LGV.Minor.roads
    ## 
    ##                                    Df Deviance    AIC
    ## - personal.motorcycles.A.roads      1   25.709 1629.0
    ## <none>                                  25.570 1629.2
    ## + freight.HGV.Minor.roads           1   25.461 1629.9
    ## + personal.buses.A.roads            1   25.472 1630.0
    ## + freight.diesel.LGV.A.roads        1   25.520 1630.6
    ## + personal.petrol.cars.A.roads      1   25.526 1630.7
    ## + freight.petrol.LGV.A.roads        1   25.526 1630.7
    ## + personal.motorcycles.Minor.roads  1   25.527 1630.7
    ## + personal.buses.Minor.roads        1   25.527 1630.7
    ## + freight.HGV.A.roads               1   25.532 1630.8
    ## + personal.diesel.cars.A.roads      1   25.534 1630.8
    ## + personal.petrol.cars.Minor.roads  1   25.541 1630.9
    ## + personal.diesel.cars.Minor.roads  1   25.547 1630.9
    ## + freight.HGV.Motorways             1   25.552 1631.0
    ## + personal.diesel.cars.Motorways    1   25.567 1631.2
    ## + freight.diesel.LGV.Motorways      1   25.569 1631.2
    ## + personal.petrol.cars.Motorway     1   25.569 1631.2
    ## + freight.petrol.LGV.Motorways      1   25.569 1631.2
    ## - personal.buses.Motorways          1   25.900 1631.4
    ## - personal.motorcycles.Motorways    1   26.508 1639.1
    ## - X.people.per.sq.km                1   28.528 1664.5
    ## - freight.diesel.LGV.Minor.roads    1   35.052 1746.8
    ## - freight.petrol.LGV.Minor.roads    1   35.567 1753.3
    ## 
    ## Step:  AIC=1628.88
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Motorways + freight.petrol.LGV.Minor.roads + 
    ##     freight.diesel.LGV.Minor.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + freight.HGV.Minor.roads           1   25.540 1628.8
    ## <none>                                  25.709 1628.9
    ## + freight.diesel.LGV.A.roads        1   25.552 1628.9
    ## + freight.petrol.LGV.A.roads        1   25.559 1629.0
    ## + personal.motorcycles.A.roads      1   25.570 1629.2
    ## + personal.petrol.cars.A.roads      1   25.574 1629.2
    ## + personal.diesel.cars.A.roads      1   25.575 1629.2
    ## + personal.motorcycles.Minor.roads  1   25.593 1629.5
    ## + freight.HGV.A.roads               1   25.622 1629.8
    ## + personal.petrol.cars.Minor.roads  1   25.639 1630.0
    ## + personal.buses.Minor.roads        1   25.650 1630.2
    ## + personal.diesel.cars.Minor.roads  1   25.663 1630.3
    ## + freight.HGV.Motorways             1   25.684 1630.6
    ## + personal.buses.A.roads            1   25.705 1630.8
    ## + personal.diesel.cars.Motorways    1   25.706 1630.8
    ## + personal.petrol.cars.Motorway     1   25.708 1630.9
    ## + freight.diesel.LGV.Motorways      1   25.708 1630.9
    ## + freight.petrol.LGV.Motorways      1   25.708 1630.9
    ## - personal.buses.Motorways          1   26.056 1631.2
    ## - personal.motorcycles.Motorways    1   26.709 1639.2
    ## - X.people.per.sq.km                1   34.372 1733.7
    ## - freight.diesel.LGV.Minor.roads    1   36.516 1760.1
    ## - freight.petrol.LGV.Minor.roads    1   36.959 1765.6
    ## 
    ## Step:  AIC=1628.87
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Motorways + freight.petrol.LGV.Minor.roads + 
    ##     freight.diesel.LGV.Minor.roads + freight.HGV.Minor.roads
    ## 
    ##                                    Df Deviance    AIC
    ## <none>                                  25.540 1628.9
    ## - freight.HGV.Minor.roads           1   25.709 1629.0
    ## + freight.diesel.LGV.A.roads        1   25.412 1629.3
    ## + freight.petrol.LGV.A.roads        1   25.417 1629.3
    ## + personal.petrol.cars.A.roads      1   25.427 1629.5
    ## + personal.diesel.cars.A.roads      1   25.434 1629.6
    ## + freight.HGV.A.roads               1   25.437 1629.6
    ## + personal.motorcycles.A.roads      1   25.461 1629.9
    ## + freight.HGV.Motorways             1   25.484 1630.2
    ## - personal.buses.Motorways          1   25.822 1630.4
    ## + personal.buses.Minor.roads        1   25.502 1630.4
    ## + personal.motorcycles.Minor.roads  1   25.507 1630.5
    ## + personal.petrol.cars.Minor.roads  1   25.513 1630.5
    ## + personal.diesel.cars.Minor.roads  1   25.523 1630.7
    ## + freight.petrol.LGV.Motorways      1   25.533 1630.8
    ## + freight.diesel.LGV.Motorways      1   25.534 1630.8
    ## + personal.petrol.cars.Motorway     1   25.536 1630.8
    ## + personal.buses.A.roads            1   25.539 1630.9
    ## + personal.diesel.cars.Motorways    1   25.539 1630.9
    ## - personal.motorcycles.Motorways    1   26.373 1637.2
    ## - X.people.per.sq.km                1   31.977 1706.6
    ## - freight.diesel.LGV.Minor.roads    1   35.062 1744.8
    ## - freight.petrol.LGV.Minor.roads    1   36.780 1766.1

    ## 
    ## Call:  glm(formula = no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.Motorways + freight.petrol.LGV.Minor.roads + 
    ##     freight.diesel.LGV.Minor.roads + freight.HGV.Minor.roads, 
    ##     family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Coefficients:
    ##                    (Intercept)              X.people.per.sq.km  
    ##                      2.375e+00                       7.034e-05  
    ## personal.motorcycles.Motorways        personal.buses.Motorways  
    ##                      1.750e-03                      -2.753e-04  
    ## freight.petrol.LGV.Minor.roads  freight.diesel.LGV.Minor.roads  
    ##                      1.722e-02                      -6.420e-04  
    ##        freight.HGV.Minor.roads  
    ##                     -1.058e-04  
    ## 
    ## Degrees of Freedom: 299 Total (i.e. Null);  293 Residual
    ## Null Deviance:       71.05 
    ## Residual Deviance: 25.54     AIC: 1629

``` r
step(model_0_glm.gamma, direction = "both", scope = formula(model_NO2_gamma.vif2, test = "Chisq"))
```

    ## Start:  AIC=1931.38
    ## no2_val ~ 1
    ## 
    ##                                    Df Deviance    AIC
    ## + X.people.per.sq.km                1   40.079 1812.6
    ## + personal.buses.A.roads            1   51.658 1857.7
    ## + personal.motorcycles.A.roads      1   61.171 1894.8
    ## + personal.motorcycles.Minor.roads  1   63.913 1905.5
    ## + freight.HGV.A.roads               1   65.739 1912.7
    ## <none>                                  71.050 1931.4
    ## + freight.HGV.Motorways             1   70.936 1932.9
    ## + personal.motorcycles.Motorways    1   71.033 1933.3
    ## 
    ## Step:  AIC=1756.48
    ## no2_val ~ X.people.per.sq.km
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.motorcycles.Motorways    1   37.618 1738.2
    ## + freight.HGV.Motorways             1   38.166 1742.7
    ## + freight.HGV.A.roads               1   39.160 1750.9
    ## + personal.motorcycles.A.roads      1   39.498 1753.7
    ## + personal.motorcycles.Minor.roads  1   39.553 1754.2
    ## + personal.buses.A.roads            1   39.658 1755.0
    ## <none>                                  40.079 1756.5
    ## - X.people.per.sq.km                1   71.050 2009.5
    ## 
    ## Step:  AIC=1739.05
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.motorcycles.A.roads      1   36.679 1732.9
    ## + freight.HGV.A.roads               1   36.782 1733.8
    ## <none>                                  37.618 1739.0
    ## + personal.buses.A.roads            1   37.442 1739.5
    ## + personal.motorcycles.Minor.roads  1   37.460 1739.7
    ## + freight.HGV.Motorways             1   37.499 1740.0
    ## - personal.motorcycles.Motorways    1   40.079 1758.4
    ## - X.people.per.sq.km                1   71.033 2026.3
    ## 
    ## Step:  AIC=1733.32
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.motorcycles.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.buses.A.roads            1   34.938 1720.6
    ## + personal.motorcycles.Minor.roads  1   35.899 1728.7
    ## + freight.HGV.A.roads               1   36.429 1733.2
    ## <none>                                  36.679 1733.3
    ## + freight.HGV.Motorways             1   36.517 1733.9
    ## - personal.motorcycles.A.roads      1   37.618 1739.3
    ## - personal.motorcycles.Motorways    1   39.498 1755.2
    ## - X.people.per.sq.km                1   60.990 1937.5
    ## 
    ## Step:  AIC=1720.44
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.motorcycles.A.roads + personal.buses.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + freight.HGV.A.roads               1   34.594 1719.5
    ## <none>                                  34.938 1720.4
    ## + personal.motorcycles.Minor.roads  1   34.908 1722.2
    ## + freight.HGV.Motorways             1   34.934 1722.4
    ## - personal.buses.A.roads            1   36.679 1733.2
    ## - personal.motorcycles.Motorways    1   37.379 1739.2
    ## - personal.motorcycles.A.roads      1   37.442 1739.7
    ## - X.people.per.sq.km                1   50.602 1851.4
    ## 
    ## Step:  AIC=1719.41
    ## no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.motorcycles.A.roads + personal.buses.A.roads + freight.HGV.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## <none>                                  34.594 1719.4
    ## - freight.HGV.A.roads               1   34.938 1720.4
    ## + personal.motorcycles.Minor.roads  1   34.498 1720.6
    ## + freight.HGV.Motorways             1   34.558 1721.1
    ## - personal.motorcycles.A.roads      1   36.140 1730.8
    ## - personal.buses.A.roads            1   36.429 1733.3
    ## - personal.motorcycles.Motorways    1   36.809 1736.6
    ## - X.people.per.sq.km                1   44.673 1804.6

    ## 
    ## Call:  glm(formula = no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.motorcycles.A.roads + personal.buses.A.roads + freight.HGV.A.roads, 
    ##     family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Coefficients:
    ##                    (Intercept)              X.people.per.sq.km  
    ##                      2.337e+00                       1.164e-04  
    ## personal.motorcycles.Motorways    personal.motorcycles.A.roads  
    ##                      1.124e-03                      -6.497e-04  
    ##         personal.buses.A.roads             freight.HGV.A.roads  
    ##                      9.757e-05                      -5.540e-06  
    ## 
    ## Degrees of Freedom: 299 Total (i.e. Null);  294 Residual
    ## Null Deviance:       71.05 
    ## Residual Deviance: 34.59     AIC: 1719

Stepwise variable selection led to similar results for the VIF-corrected
and original model. The results from both analysis suggest to include a
smaller number of variables in our final model.

``` r
summary(model_NO2_gamma.final <- glm(formula = no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    personal.motorcycles.A.roads + personal.buses.A.roads + freight.HGV.A.roads, 
    family = Gamma(link = "log"), data = fossilfuel_2017.omit))
```

    ## 
    ## Call:
    ## glm(formula = no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.motorcycles.A.roads + personal.buses.A.roads + freight.HGV.A.roads, 
    ##     family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.46058  -0.22360  -0.00582   0.17644   1.52810  
    ## 
    ## Coefficients:
    ##                                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                     2.337e+00  3.832e-02  60.995  < 2e-16 ***
    ## X.people.per.sq.km              1.164e-04  1.165e-05   9.998  < 2e-16 ***
    ## personal.motorcycles.Motorways  1.124e-03  2.557e-04   4.395 1.55e-05 ***
    ## personal.motorcycles.A.roads   -6.497e-04  1.818e-04  -3.574 0.000410 ***
    ## personal.buses.A.roads          9.757e-05  2.510e-05   3.888 0.000125 ***
    ## freight.HGV.A.roads            -5.540e-06  3.196e-06  -1.734 0.084038 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1155419)
    ## 
    ##     Null deviance: 71.050  on 299  degrees of freedom
    ## Residual deviance: 34.594  on 294  degrees of freedom
    ## AIC: 1719.4
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
# freight.HGV.A.roads is not significant so it can be dropped from the model
summary(model_NO2_gamma.final.v2 <- glm(formula = no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    personal.motorcycles.A.roads + personal.buses.A.roads, 
    family = Gamma(link = "log"), data = fossilfuel_2017.omit))
```

    ## 
    ## Call:
    ## glm(formula = no2_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.motorcycles.A.roads + personal.buses.A.roads, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.47768  -0.22557  -0.01073   0.18751   1.61676  
    ## 
    ## Coefficients:
    ##                                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                     2.298e+00  3.238e-02  70.977  < 2e-16 ***
    ## X.people.per.sq.km              1.273e-04  1.019e-05  12.485  < 2e-16 ***
    ## personal.motorcycles.Motorways  1.168e-03  2.569e-04   4.546 8.00e-06 ***
    ## personal.motorcycles.A.roads   -7.643e-04  1.663e-04  -4.596 6.39e-06 ***
    ## personal.buses.A.roads          9.500e-05  2.534e-05   3.749 0.000213 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1178213)
    ## 
    ##     Null deviance: 71.050  on 299  degrees of freedom
    ## Residual deviance: 34.938  on 295  degrees of freedom
    ## AIC: 1720.4
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
vif(model_NO2_gamma.final.v2)
```

    ##             X.people.per.sq.km personal.motorcycles.Motorways 
    ##                       1.901923                       1.057647 
    ##   personal.motorcycles.A.roads         personal.buses.A.roads 
    ##                       2.790539                       3.189693

``` r
summary(model_NO2_gamma.rsd <- glm(data = fossilfuel_2017.omit, no2_val ~ X.people.per.sq.km + petroleum.Industrial + petroleum.Domestic + petroleum.Rail + petroleum.Public.Administration + petroleum.Commercial + petroleum.Agriculture + coal.Domestic + Manufactured.Solid.Fuels.Domestic + bioenergy, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = no2_val ~ X.people.per.sq.km + petroleum.Industrial + 
    ##     petroleum.Domestic + petroleum.Rail + petroleum.Public.Administration + 
    ##     petroleum.Commercial + petroleum.Agriculture + coal.Domestic + 
    ##     Manufactured.Solid.Fuels.Domestic + bioenergy, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.39314  -0.20144  -0.00765   0.14356   1.13692  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        2.367e+00  3.167e-02  74.730  < 2e-16 ***
    ## X.people.per.sq.km                 8.051e-05  7.054e-06  11.414  < 2e-16 ***
    ## petroleum.Industrial              -1.836e-04  2.173e-04  -0.845 0.398903    
    ## petroleum.Domestic                -1.357e-02  6.693e-03  -2.027 0.043542 *  
    ## petroleum.Rail                     2.200e-02  9.333e-03   2.358 0.019052 *  
    ## petroleum.Public.Administration   -1.075e-01  6.342e-02  -1.695 0.091152 .  
    ## petroleum.Commercial               2.462e-01  8.670e-02   2.839 0.004842 ** 
    ## petroleum.Agriculture             -6.000e-02  9.493e-03  -6.321 9.83e-10 ***
    ## coal.Domestic                     -2.744e-01  7.624e-02  -3.599 0.000376 ***
    ## Manufactured.Solid.Fuels.Domestic  5.847e-01  1.129e-01   5.179 4.18e-07 ***
    ## bioenergy                         -1.662e-03  1.034e-03  -1.608 0.108933    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.07949522)
    ## 
    ##     Null deviance: 71.050  on 299  degrees of freedom
    ## Residual deviance: 24.228  on 289  degrees of freedom
    ## AIC: 1620.8
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
vif(model_NO2_gamma.rsd)
```

    ##                X.people.per.sq.km              petroleum.Industrial 
    ##                          1.349871                          1.042394 
    ##                petroleum.Domestic                    petroleum.Rail 
    ##                          6.433641                          1.588282 
    ##   petroleum.Public.Administration              petroleum.Commercial 
    ##                          1.172803                          2.362592 
    ##             petroleum.Agriculture                     coal.Domestic 
    ##                          5.424378                         15.846945 
    ## Manufactured.Solid.Fuels.Domestic                         bioenergy 
    ##                         13.300006                          1.134588

``` r
summary(model_NO2_gamma.rsd.vif2 <- glm(data = fossilfuel_2017.omit, no2_val ~ X.people.per.sq.km + petroleum.Industrial + petroleum.Rail + petroleum.Public.Administration + petroleum.Commercial + petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + bioenergy, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = no2_val ~ X.people.per.sq.km + petroleum.Industrial + 
    ##     petroleum.Rail + petroleum.Public.Administration + petroleum.Commercial + 
    ##     petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     bioenergy, family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.41673  -0.22125  -0.02595   0.15814   1.13330  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        2.376e+00  3.220e-02  73.792  < 2e-16 ***
    ## X.people.per.sq.km                 9.043e-05  6.839e-06  13.224  < 2e-16 ***
    ## petroleum.Industrial              -2.506e-04  2.223e-04  -1.127   0.2607    
    ## petroleum.Rail                     2.267e-02  9.533e-03   2.378   0.0180 *  
    ## petroleum.Public.Administration   -1.261e-01  6.303e-02  -2.000   0.0464 *  
    ## petroleum.Commercial               2.010e-01  8.092e-02   2.484   0.0136 *  
    ## petroleum.Agriculture             -8.371e-02  6.968e-03 -12.013  < 2e-16 ***
    ## Manufactured.Solid.Fuels.Domestic  2.041e-01  5.161e-02   3.955 9.61e-05 ***
    ## bioenergy                         -1.165e-03  1.041e-03  -1.119   0.2640    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.0833698)
    ## 
    ##     Null deviance: 71.050  on 299  degrees of freedom
    ## Residual deviance: 25.555  on 291  degrees of freedom
    ## AIC: 1633
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
vif(model_NO2_gamma.rsd.vif2)
```

    ##                X.people.per.sq.km              petroleum.Industrial 
    ##                          1.209843                          1.040351 
    ##                    petroleum.Rail   petroleum.Public.Administration 
    ##                          1.580011                          1.104873 
    ##              petroleum.Commercial             petroleum.Agriculture 
    ##                          1.962510                          2.787182 
    ## Manufactured.Solid.Fuels.Domestic                         bioenergy 
    ##                          2.650664                          1.096700

Same as before, we tested all predictors of the original model for
collinearity and dropped predictors with VIF\>5.We will apply a stepwise
regression method to only include those variables that have an effect on
NO2 concentration in our final model.

``` r
### Define full and null models and do step procedure
model_0_glm.gamma.2 <- glm(no2_val ~ 1, data = fossilfuel_2017.omit, family = Gamma(link = "log"))
step(model_0_glm.gamma, direction = "both", scope = formula(model_NO2_gamma.rsd.vif2, test = "Chisq"))
```

    ## Start:  AIC=1931.38
    ## no2_val ~ 1
    ## 
    ##                                     Df Deviance    AIC
    ## + X.people.per.sq.km                 1   40.079 1812.6
    ## + petroleum.Agriculture              1   47.651 1842.1
    ## + Manufactured.Solid.Fuels.Domestic  1   64.818 1909.1
    ## + petroleum.Rail                     1   68.638 1924.0
    ## + bioenergy                          1   68.869 1924.9
    ## + petroleum.Industrial               1   70.084 1929.6
    ## <none>                                   71.050 1931.4
    ## + petroleum.Commercial               1   70.549 1931.4
    ## + petroleum.Public.Administration    1   70.638 1931.8
    ## 
    ## Step:  AIC=1756.48
    ## no2_val ~ X.people.per.sq.km
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Agriculture              1   29.661 1672.7
    ## + petroleum.Commercial               1   38.672 1746.9
    ## + Manufactured.Solid.Fuels.Domestic  1   38.977 1749.4
    ## + bioenergy                          1   39.492 1753.6
    ## + petroleum.Public.Administration    1   39.676 1755.2
    ## <none>                                   40.079 1756.5
    ## + petroleum.Industrial               1   39.910 1757.1
    ## + petroleum.Rail                     1   39.995 1757.8
    ## - X.people.per.sq.km                 1   71.050 2009.5
    ## 
    ## Step:  AIC=1666.43
    ## no2_val ~ X.people.per.sq.km + petroleum.Agriculture
    ## 
    ##                                     Df Deviance    AIC
    ## + Manufactured.Solid.Fuels.Domestic  1   26.973 1640.4
    ## + petroleum.Rail                     1   28.187 1653.1
    ## + petroleum.Commercial               1   28.337 1654.6
    ## <none>                                   29.661 1666.4
    ## + petroleum.Public.Administration    1   29.641 1668.2
    ## + petroleum.Industrial               1   29.650 1668.3
    ## + bioenergy                          1   29.656 1668.4
    ## - petroleum.Agriculture              1   40.079 1773.1
    ## - X.people.per.sq.km                 1   47.651 1852.1
    ## 
    ## Step:  AIC=1639.49
    ## no2_val ~ X.people.per.sq.km + petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Commercial               1   26.429 1635.5
    ## + petroleum.Rail                     1   26.536 1636.6
    ## <none>                                   26.973 1639.5
    ## + bioenergy                          1   26.822 1639.8
    ## + petroleum.Public.Administration    1   26.864 1640.3
    ## + petroleum.Industrial               1   26.914 1640.8
    ## - Manufactured.Solid.Fuels.Domestic  1   29.661 1667.4
    ## - petroleum.Agriculture              1   38.977 1770.9
    ## - X.people.per.sq.km                 1   45.283 1841.0
    ## 
    ## Step:  AIC=1635.29
    ## no2_val ~ X.people.per.sq.km + petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     petroleum.Commercial
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Rail                     1   26.084 1633.3
    ## + petroleum.Public.Administration    1   26.233 1635.0
    ## <none>                                   26.429 1635.3
    ## + bioenergy                          1   26.310 1635.9
    ## + petroleum.Industrial               1   26.353 1636.4
    ## - petroleum.Commercial               1   26.973 1639.6
    ## - Manufactured.Solid.Fuels.Domestic  1   28.337 1655.5
    ## - petroleum.Agriculture              1   38.379 1772.5
    ## - X.people.per.sq.km                 1   42.031 1815.0
    ## 
    ## Step:  AIC=1633.29
    ## no2_val ~ X.people.per.sq.km + petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     petroleum.Commercial + petroleum.Rail
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Public.Administration    1   25.774 1631.7
    ## <none>                                   26.084 1633.3
    ## + bioenergy                          1   25.971 1634.0
    ## + petroleum.Industrial               1   25.979 1634.0
    ## - petroleum.Rail                     1   26.429 1635.3
    ## - petroleum.Commercial               1   26.536 1636.6
    ## - Manufactured.Solid.Fuels.Domestic  1   27.296 1645.5
    ## - petroleum.Agriculture              1   38.223 1773.6
    ## - X.people.per.sq.km                 1   41.945 1817.2
    ## 
    ## Step:  AIC=1631.65
    ## no2_val ~ X.people.per.sq.km + petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     petroleum.Commercial + petroleum.Rail + petroleum.Public.Administration
    ## 
    ##                                     Df Deviance    AIC
    ## <none>                                   25.774 1631.7
    ## + petroleum.Industrial               1   25.651 1632.2
    ## + bioenergy                          1   25.660 1632.3
    ## - petroleum.Public.Administration    1   26.084 1633.3
    ## - petroleum.Rail                     1   26.233 1635.1
    ## - petroleum.Commercial               1   26.317 1636.1
    ## - Manufactured.Solid.Fuels.Domestic  1   26.999 1644.2
    ## - petroleum.Agriculture              1   38.030 1775.7
    ## - X.people.per.sq.km                 1   41.729 1819.7

    ## 
    ## Call:  glm(formula = no2_val ~ X.people.per.sq.km + petroleum.Agriculture + 
    ##     Manufactured.Solid.Fuels.Domestic + petroleum.Commercial + 
    ##     petroleum.Rail + petroleum.Public.Administration, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Coefficients:
    ##                       (Intercept)                 X.people.per.sq.km  
    ##                         2.369e+00                          9.114e-05  
    ##             petroleum.Agriculture  Manufactured.Solid.Fuels.Domestic  
    ##                        -8.388e-02                          1.924e-01  
    ##              petroleum.Commercial                     petroleum.Rail  
    ##                         2.026e-01                          2.190e-02  
    ##   petroleum.Public.Administration  
    ##                        -1.228e-01  
    ## 
    ## Degrees of Freedom: 299 Total (i.e. Null);  293 Residual
    ## Null Deviance:       71.05 
    ## Residual Deviance: 25.77     AIC: 1632

Based on the results of stepwise regression from the VIF-corrected
model, we will now call the following final model:

``` r
### Define final model based on stepwise regression results
summary(model_NO2_gamma.final.rsd <- glm(formula = no2_val ~ X.people.per.sq.km + petroleum.Agriculture + petroleum.Commercial + Manufactured.Solid.Fuels.Domestic + petroleum.Rail, family = Gamma(link = "log"), 
    data = fossilfuel_2017.omit))
```

    ## 
    ## Call:
    ## glm(formula = no2_val ~ X.people.per.sq.km + petroleum.Agriculture + 
    ##     petroleum.Commercial + Manufactured.Solid.Fuels.Domestic + 
    ##     petroleum.Rail, family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.41166  -0.22105  -0.02232   0.16275   1.15941  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        2.370e+00  3.226e-02  73.470  < 2e-16 ***
    ## X.people.per.sq.km                 9.088e-05  6.905e-06  13.161  < 2e-16 ***
    ## petroleum.Agriculture             -8.316e-02  7.008e-03 -11.866  < 2e-16 ***
    ## petroleum.Commercial               1.827e-01  8.076e-02   2.262 0.024420 *  
    ## Manufactured.Solid.Fuels.Domestic  1.907e-01  5.161e-02   3.696 0.000261 ***
    ## petroleum.Rail                     1.871e-02  9.452e-03   1.979 0.048722 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.08530947)
    ## 
    ##     Null deviance: 71.050  on 299  degrees of freedom
    ## Residual deviance: 26.084  on 294  degrees of freedom
    ## AIC: 1633.3
    ## 
    ## Number of Fisher Scoring iterations: 8

``` r
vif(model_NO2_gamma.final.rsd) # VIF is <5
```

    ##                X.people.per.sq.km             petroleum.Agriculture 
    ##                          1.205406                          2.754995 
    ##              petroleum.Commercial Manufactured.Solid.Fuels.Domestic 
    ##                          1.909997                          2.590281 
    ##                    petroleum.Rail 
    ##                          1.518114

``` r
#plot(model_NO2_gamma.final.rsd) # no significant outlier
```

Our results suggest that an increase in petroleum consumption Commercial
and Agricultural activities is associated with increased levels of NO2.
Manufactured solid fuel consumption for domestic purposes also appear as
a strong predictor of NO2 levels.

Last, we will look at the effect of domestic and non-domestic gas
consumption:

``` r
summary(model_NO2_gamma.gas <- glm(data = fossilfuel_2017.omit, no2_val ~  X.people.per.sq.km + Domestic.mean.gas.consumption + Non.domestic.mean.gas.consumption, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = no2_val ~ X.people.per.sq.km + Domestic.mean.gas.consumption + 
    ##     Non.domestic.mean.gas.consumption, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.52285  -0.25481  -0.00762   0.20148   1.38140  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        1.919e+00  2.024e-01   9.479   <2e-16 ***
    ## X.people.per.sq.km                 1.246e-04  8.014e-06  15.542   <2e-16 ***
    ## Domestic.mean.gas.consumption      3.268e-05  1.412e-05   2.314   0.0214 *  
    ## Non.domestic.mean.gas.consumption -3.483e-08  3.658e-08  -0.952   0.3418    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1223536)
    ## 
    ##     Null deviance: 71.050  on 299  degrees of freedom
    ## Residual deviance: 39.149  on 296  degrees of freedom
    ## AIC: 1753.3
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
# Non domestic gas consumption is not significant so will be dropped
summary(model_NO2_gamma.gas.v2 <- glm(data = fossilfuel_2017.omit, no2_val ~  X.people.per.sq.km + Domestic.mean.gas.consumption, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = no2_val ~ X.people.per.sq.km + Domestic.mean.gas.consumption, 
    ##     family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.52758  -0.25816  -0.01262   0.20260   1.37855  
    ## 
    ## Coefficients:
    ##                                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   1.847e+00  1.892e-01   9.763   <2e-16 ***
    ## X.people.per.sq.km            1.264e-04  7.816e-06  16.174   <2e-16 ***
    ## Domestic.mean.gas.consumption 3.600e-05  1.375e-05   2.618   0.0093 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1223722)
    ## 
    ##     Null deviance: 71.050  on 299  degrees of freedom
    ## Residual deviance: 39.255  on 297  degrees of freedom
    ## AIC: 1752.1
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
vif(model_NO2_gamma.gas.v2) # no collinearity 
```

    ##            X.people.per.sq.km Domestic.mean.gas.consumption 
    ##                       1.07663                       1.07663

``` r
#plot(model_NO2_gamma.gas.v2) # no outlier
```

No stepwise regression is needed with only three explanatory variables.

#### NOx

Same as above: analyse the data distribution of the target or response
variable (NO concentration). To do this, we will use the following
functions:

``` r
hist(fossilfuel_2017.omit$nox_val, main="Histogram of Yield", xlab="Yield (quintals/ha)") 
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->

``` r
#qqnorm(fossilfuel_2017.omit$nox_val, main="QQplot of NOx")
```

Data is non-parametric and gamma distributed.

``` r
library(e1071)
skewness(fossilfuel_2017.omit$nox_val) 
```

    ## [1] 2.102793

Based on the results of VIF correction we call the following model:

``` r
summary(model_NOx_gamma.vif2 <- glm(data = fossilfuel_2017.omit, nox_val ~ X.people.per.sq.km + personal.buses.A.roads + personal.motorcycles.Motorways + personal.motorcycles.A.roads + personal.motorcycles.Minor.roads + freight.HGV.Motorways + freight.HGV.A.roads, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = nox_val ~ X.people.per.sq.km + personal.buses.A.roads + 
    ##     personal.motorcycles.Motorways + personal.motorcycles.A.roads + 
    ##     personal.motorcycles.Minor.roads + freight.HGV.Motorways + 
    ##     freight.HGV.A.roads, family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.59280  -0.26674  -0.01698   0.18419   1.98180  
    ## 
    ## Coefficients:
    ##                                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                       2.602e+00  5.001e-02  52.027  < 2e-16 ***
    ## X.people.per.sq.km                1.348e-04  1.385e-05   9.728  < 2e-16 ***
    ## personal.buses.A.roads            9.069e-05  3.532e-05   2.567  0.01074 *  
    ## personal.motorcycles.Motorways    9.346e-04  4.574e-04   2.043  0.04194 *  
    ## personal.motorcycles.A.roads     -5.913e-04  2.170e-04  -2.724  0.00684 ** 
    ## personal.motorcycles.Minor.roads  1.694e-04  2.384e-04   0.711  0.47788    
    ## freight.HGV.Motorways             1.741e-06  2.679e-06   0.650  0.51623    
    ## freight.HGV.A.roads              -7.974e-06  3.942e-06  -2.023  0.04400 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1591281)
    ## 
    ##     Null deviance: 94.264  on 299  degrees of freedom
    ## Residual deviance: 43.137  on 292  degrees of freedom
    ## AIC: 1977
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
vif(model_NOx_gamma.vif2)
```

    ##               X.people.per.sq.km           personal.buses.A.roads 
    ##                         2.601107                         4.589445 
    ##   personal.motorcycles.Motorways     personal.motorcycles.A.roads 
    ##                         2.483397                         3.520010 
    ## personal.motorcycles.Minor.roads            freight.HGV.Motorways 
    ##                         2.200378                         2.546692 
    ##              freight.HGV.A.roads 
    ##                         1.654252

Stepwise regression:

``` r
### Define full and null models and do step procedure
model_NOx.0_glm.gamma <- glm(nox_val ~ 1, data = fossilfuel_2017.omit, family = Gamma(link = "log"))
step(model_NOx.0_glm.gamma, direction = "both", scope = formula(model_NOx_gamma.vif2, test = "Chisq"))
```

    ## Start:  AIC=2205.96
    ## nox_val ~ 1
    ## 
    ##                                    Df Deviance    AIC
    ## + X.people.per.sq.km                1   49.659 2092.1
    ## + personal.buses.A.roads            1   64.744 2131.3
    ## + personal.motorcycles.A.roads      1   77.444 2164.3
    ## + personal.motorcycles.Minor.roads  1   83.395 2179.7
    ## + freight.HGV.A.roads               1   87.282 2189.8
    ## <none>                                  94.264 2206.0
    ## + freight.HGV.Motorways             1   93.943 2207.1
    ## + personal.motorcycles.Motorways    1   94.257 2207.9
    ## 
    ## Step:  AIC=2008.29
    ## nox_val ~ X.people.per.sq.km
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.motorcycles.Motorways    1   46.721 1992.1
    ## + freight.HGV.Motorways             1   47.248 1995.3
    ## + freight.HGV.A.roads               1   48.659 2004.1
    ## + personal.buses.A.roads            1   48.876 2005.4
    ## + personal.motorcycles.Minor.roads  1   48.898 2005.6
    ## + personal.motorcycles.A.roads      1   49.267 2007.8
    ## <none>                                  49.659 2008.3
    ## - X.people.per.sq.km                1   94.264 2283.0
    ## 
    ## Step:  AIC=1991.5
    ## nox_val ~ X.people.per.sq.km + personal.motorcycles.Motorways
    ## 
    ##                                    Df Deviance    AIC
    ## + freight.HGV.A.roads               1   45.813 1987.7
    ## + personal.motorcycles.A.roads      1   46.006 1988.9
    ## + personal.buses.A.roads            1   46.322 1990.9
    ## <none>                                  46.721 1991.5
    ## + personal.motorcycles.Minor.roads  1   46.461 1991.8
    ## + freight.HGV.Motorways             1   46.531 1992.3
    ## - personal.motorcycles.Motorways    1   49.659 2008.4
    ## - X.people.per.sq.km                1   94.257 2295.9
    ## 
    ## Step:  AIC=1987.47
    ## nox_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     freight.HGV.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.buses.A.roads            1   44.729 1982.4
    ## + personal.motorcycles.Minor.roads  1   44.933 1983.7
    ## + freight.HGV.Motorways             1   45.389 1986.7
    ## <none>                                  45.813 1987.5
    ## + personal.motorcycles.A.roads      1   45.622 1988.2
    ## - freight.HGV.A.roads               1   46.721 1991.4
    ## - personal.motorcycles.Motorways    1   48.659 2004.0
    ## - X.people.per.sq.km                1   87.276 2255.3
    ## 
    ## Step:  AIC=1982.1
    ## nox_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     freight.HGV.A.roads + personal.buses.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.motorcycles.A.roads      1   43.280 1974.3
    ## <none>                                  44.729 1982.1
    ## + freight.HGV.Motorways             1   44.514 1982.7
    ## + personal.motorcycles.Minor.roads  1   44.607 1983.3
    ## - personal.buses.A.roads            1   45.813 1987.4
    ## - freight.HGV.A.roads               1   46.322 1990.8
    ## - personal.motorcycles.Motorways    1   46.909 1994.8
    ## - X.people.per.sq.km                1   56.042 2056.2
    ## 
    ## Step:  AIC=1973.98
    ## nox_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     freight.HGV.A.roads + personal.buses.A.roads + personal.motorcycles.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## <none>                                  43.280 1974.0
    ## - freight.HGV.A.roads               1   43.814 1975.3
    ## + personal.motorcycles.Minor.roads  1   43.205 1975.5
    ## + freight.HGV.Motorways             1   43.217 1975.6
    ## - personal.motorcycles.A.roads      1   44.729 1981.1
    ## - personal.buses.A.roads            1   45.622 1986.8
    ## - personal.motorcycles.Motorways    1   45.761 1987.7
    ## - X.people.per.sq.km                1   56.041 2052.6

    ## 
    ## Call:  glm(formula = nox_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     freight.HGV.A.roads + personal.buses.A.roads + personal.motorcycles.A.roads, 
    ##     family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Coefficients:
    ##                    (Intercept)              X.people.per.sq.km  
    ##                      2.620e+00                       1.321e-04  
    ## personal.motorcycles.Motorways             freight.HGV.A.roads  
    ##                      1.191e-03                      -6.888e-06  
    ##         personal.buses.A.roads    personal.motorcycles.A.roads  
    ##                      1.098e-04                      -6.268e-04  
    ## 
    ## Degrees of Freedom: 299 Total (i.e. Null);  294 Residual
    ## Null Deviance:       94.26 
    ## Residual Deviance: 43.28     AIC: 1974

Based on the results of the stepwise regression, we proceed to build our
final model:

``` r
### Define final model based on stepwise regression results
summary(model_NOx_gamma_final <- glm(formula = nox_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    freight.HGV.A.roads + personal.buses.A.roads + personal.motorcycles.A.roads, 
    family = Gamma(link = "log"), data = fossilfuel_2017.omit))
```

    ## 
    ## Call:
    ## glm(formula = nox_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     freight.HGV.A.roads + personal.buses.A.roads + personal.motorcycles.A.roads, 
    ##     family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.58647  -0.25701  -0.02948   0.18730   1.95719  
    ## 
    ## Coefficients:
    ##                                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                     2.620e+00  4.486e-02  58.409  < 2e-16 ***
    ## X.people.per.sq.km              1.321e-04  1.363e-05   9.687  < 2e-16 ***
    ## personal.motorcycles.Motorways  1.191e-03  2.993e-04   3.979 8.71e-05 ***
    ## freight.HGV.A.roads            -6.888e-06  3.741e-06  -1.841 0.066593 .  
    ## personal.buses.A.roads          1.098e-04  2.938e-05   3.737 0.000224 ***
    ## personal.motorcycles.A.roads   -6.268e-04  2.128e-04  -2.946 0.003480 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1583252)
    ## 
    ##     Null deviance: 94.264  on 299  degrees of freedom
    ## Residual deviance: 43.280  on 294  degrees of freedom
    ## AIC: 1974
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
# Drop freight.HGV.A.roads because not significant
summary(model_NOx_gamma_final <-glm(formula = nox_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + personal.buses.A.roads + personal.motorcycles.A.roads, 
    family = Gamma(link = "log"), data = fossilfuel_2017.omit))
```

    ## 
    ## Call:
    ## glm(formula = nox_val ~ X.people.per.sq.km + personal.motorcycles.Motorways + 
    ##     personal.buses.A.roads + personal.motorcycles.A.roads, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.60711  -0.26532  -0.02905   0.20043   2.07780  
    ## 
    ## Coefficients:
    ##                                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                     2.571e+00  3.820e-02  67.309  < 2e-16 ***
    ## X.people.per.sq.km              1.456e-04  1.203e-05  12.109  < 2e-16 ***
    ## personal.motorcycles.Motorways  1.244e-03  3.030e-04   4.106 5.22e-05 ***
    ## personal.buses.A.roads          1.062e-04  2.989e-05   3.552 0.000445 ***
    ## personal.motorcycles.A.roads   -7.651e-04  1.962e-04  -3.900 0.000119 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1639815)
    ## 
    ##     Null deviance: 94.264  on 299  degrees of freedom
    ## Residual deviance: 43.814  on 295  degrees of freedom
    ## AIC: 1975.7
    ## 
    ## Number of Fisher Scoring iterations: 7

Our final NOx model indicates significant effects of same on-road
vehicles reported for NO2. Let’s have a look at fossil fuel data from
manufacturing and construction site using VIF-corrected predictors:

``` r
summary(model_NOx_gamma.rsd.vif <- glm(data = fossilfuel_2017.omit, nox_val ~ X.people.per.sq.km + petroleum.Industrial + petroleum.Rail + petroleum.Public.Administration + petroleum.Commercial + petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + bioenergy, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = nox_val ~ X.people.per.sq.km + petroleum.Industrial + 
    ##     petroleum.Rail + petroleum.Public.Administration + petroleum.Commercial + 
    ##     petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     bioenergy, family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.54504  -0.24882  -0.04069   0.16589   1.46284  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        2.659e+00  3.725e-02  71.366  < 2e-16 ***
    ## X.people.per.sq.km                 1.091e-04  7.912e-06  13.788  < 2e-16 ***
    ## petroleum.Industrial              -2.831e-04  2.573e-04  -1.100 0.272067    
    ## petroleum.Rail                     2.389e-02  1.103e-02   2.166 0.031148 *  
    ## petroleum.Public.Administration   -1.453e-01  7.293e-02  -1.992 0.047269 *  
    ## petroleum.Commercial               2.735e-01  9.362e-02   2.921 0.003757 ** 
    ## petroleum.Agriculture             -9.208e-02  8.062e-03 -11.422  < 2e-16 ***
    ## Manufactured.Solid.Fuels.Domestic  2.135e-01  5.971e-02   3.576 0.000408 ***
    ## bioenergy                         -1.185e-03  1.204e-03  -0.984 0.326053    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1115973)
    ## 
    ##     Null deviance: 94.264  on 299  degrees of freedom
    ## Residual deviance: 32.422  on 291  degrees of freedom
    ## AIC: 1891.5
    ## 
    ## Number of Fisher Scoring iterations: 8

``` r
vif(model_NOx_gamma.rsd.vif) # no collinearity
```

    ##                X.people.per.sq.km              petroleum.Industrial 
    ##                          1.209843                          1.040351 
    ##                    petroleum.Rail   petroleum.Public.Administration 
    ##                          1.580011                          1.104873 
    ##              petroleum.Commercial             petroleum.Agriculture 
    ##                          1.962510                          2.787182 
    ## Manufactured.Solid.Fuels.Domestic                         bioenergy 
    ##                          2.650664                          1.096700

Stepwise regression:

``` r
### Define full and null models and do step procedure
model_NOx.0_glm.gamma <- glm(nox_val ~ 1, data = fossilfuel_2017.omit, family = Gamma(link = "log"))
step(model_NOx.0_glm.gamma, direction = "both", scope = formula(model_NOx_gamma.rsd.vif, test = "Chisq"))
```

    ## Start:  AIC=2205.96
    ## nox_val ~ 1

    ## Warning: glm.fit: algorithm did not converge

    ##                                     Df Deviance    AIC
    ## + X.people.per.sq.km                 1   49.659 2092.1
    ## + petroleum.Agriculture              1   64.933 2131.8
    ## + Manufactured.Solid.Fuels.Domestic  1   85.988 2186.5
    ## + petroleum.Rail                     1   90.872 2199.1
    ## + bioenergy                          1   91.400 2200.5
    ## + petroleum.Industrial               1   92.957 2204.6
    ## <none>                                   94.264 2206.0
    ## + petroleum.Public.Administration    1   93.824 2206.8
    ## + petroleum.Commercial               1   94.154 2207.7
    ## 
    ## Step:  AIC=2008.29
    ## nox_val ~ X.people.per.sq.km
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Agriculture              1   37.616 1935.6
    ## + Manufactured.Solid.Fuels.Domestic  1   48.458 2002.8
    ## + petroleum.Commercial               1   48.478 2003.0
    ## + bioenergy                          1   49.012 2006.3
    ## + petroleum.Public.Administration    1   49.169 2007.2
    ## <none>                                   49.659 2008.3
    ## + petroleum.Industrial               1   49.460 2009.0
    ## + petroleum.Rail                     1   49.572 2009.8
    ## - X.people.per.sq.km                 1   94.264 2283.0
    ## 
    ## Step:  AIC=1924.96
    ## nox_val ~ X.people.per.sq.km + petroleum.Agriculture

    ## Warning: glm.fit: algorithm did not converge

    ##                                     Df Deviance    AIC
    ## + Manufactured.Solid.Fuels.Domestic  1   34.462 1902.7
    ## + petroleum.Commercial               1   35.569 1911.2
    ## + petroleum.Rail                     1   35.895 1913.7
    ## <none>                                   37.616 1925.0
    ## + petroleum.Public.Administration    1   37.588 1926.8
    ## + petroleum.Industrial               1   37.603 1926.9
    ## + bioenergy                          1   37.613 1926.9
    ## - petroleum.Agriculture              1   49.659 2015.7
    ## - X.people.per.sq.km                 1   64.933 2133.4
    ## 
    ## Step:  AIC=1900.16
    ## nox_val ~ X.people.per.sq.km + petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Commercial               1   33.458 1894.1
    ## + petroleum.Rail                     1   33.964 1898.2
    ## <none>                                   34.462 1900.2
    ## + bioenergy                          1   34.290 1900.8
    ## + petroleum.Public.Administration    1   34.321 1901.0
    ## + petroleum.Industrial               1   34.389 1901.6
    ## - Manufactured.Solid.Fuels.Domestic  1   37.616 1923.4
    ## - petroleum.Agriculture              1   48.458 2010.1
    ## - X.people.per.sq.km                 1   62.350 2121.2
    ## 
    ## Step:  AIC=1893.12
    ## nox_val ~ X.people.per.sq.km + petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     petroleum.Commercial
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Rail                     1   33.089 1891.9
    ## + petroleum.Public.Administration    1   33.184 1892.7
    ## <none>                                   33.458 1893.1
    ## + bioenergy                          1   33.331 1894.0
    ## + petroleum.Industrial               1   33.359 1894.3
    ## - petroleum.Commercial               1   34.462 1899.9
    ## - Manufactured.Solid.Fuels.Domestic  1   35.569 1909.6
    ## - petroleum.Agriculture              1   48.031 2018.5
    ## - X.people.per.sq.km                 1   56.781 2094.9
    ## 
    ## Step:  AIC=1891.74
    ## nox_val ~ X.people.per.sq.km + petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     petroleum.Commercial + petroleum.Rail
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Public.Administration    1   32.677 1890.1
    ## <none>                                   33.089 1891.7
    ## + petroleum.Industrial               1   32.956 1892.6
    ## + bioenergy                          1   32.969 1892.7
    ## - petroleum.Rail                     1   33.458 1893.0
    ## - petroleum.Commercial               1   33.964 1897.4
    ## - Manufactured.Solid.Fuels.Domestic  1   34.416 1901.3
    ## - petroleum.Agriculture              1   47.869 2019.0
    ## - X.people.per.sq.km                 1   56.716 2096.3
    ## 
    ## Step:  AIC=1889.92
    ## nox_val ~ X.people.per.sq.km + petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     petroleum.Commercial + petroleum.Rail + petroleum.Public.Administration
    ## 
    ##                                     Df Deviance    AIC
    ## <none>                                   32.677 1889.9
    ## + petroleum.Industrial               1   32.522 1890.5
    ## + bioenergy                          1   32.556 1890.8
    ## - petroleum.Public.Administration    1   33.089 1891.6
    ## - petroleum.Rail                     1   33.184 1892.4
    ## - petroleum.Commercial               1   33.692 1897.0
    ## - Manufactured.Solid.Fuels.Domestic  1   34.022 1899.9
    ## - petroleum.Agriculture              1   47.597 2020.8
    ## - X.people.per.sq.km                 1   56.461 2099.8

    ## 
    ## Call:  glm(formula = nox_val ~ X.people.per.sq.km + petroleum.Agriculture + 
    ##     Manufactured.Solid.Fuels.Domestic + petroleum.Commercial + 
    ##     petroleum.Rail + petroleum.Public.Administration, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Coefficients:
    ##                       (Intercept)                 X.people.per.sq.km  
    ##                         2.6513472                          0.0001099  
    ##             petroleum.Agriculture  Manufactured.Solid.Fuels.Domestic  
    ##                        -0.0922126                          0.2012087  
    ##              petroleum.Commercial                     petroleum.Rail  
    ##                         0.2748911                          0.0229789  
    ##   petroleum.Public.Administration  
    ##                        -0.1416225  
    ## 
    ## Degrees of Freedom: 299 Total (i.e. Null);  293 Residual
    ## Null Deviance:       94.26 
    ## Residual Deviance: 32.68     AIC: 1890

Based on the results of the stepwise regression on the original
VIF-corrected model, we build our final model:

``` r
### Define final model based on stepwise regression results
summary(model_NOx_gamma.final.rsd <- glm(formula = nox_val ~ X.people.per.sq.km + petroleum.Agriculture + 
    Manufactured.Solid.Fuels.Domestic + petroleum.Commercial + 
    petroleum.Rail + petroleum.Public.Administration, family = Gamma(link = "log"), 
    data = fossilfuel_2017.omit))
```

    ## 
    ## Call:
    ## glm(formula = nox_val ~ X.people.per.sq.km + petroleum.Agriculture + 
    ##     Manufactured.Solid.Fuels.Domestic + petroleum.Commercial + 
    ##     petroleum.Rail + petroleum.Public.Administration, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.5454  -0.2522  -0.0370   0.1679   1.4670  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        2.651e+00  3.700e-02  71.649  < 2e-16 ***
    ## X.people.per.sq.km                 1.099e-04  7.921e-06  13.872  < 2e-16 ***
    ## petroleum.Agriculture             -9.221e-02  8.045e-03 -11.462  < 2e-16 ***
    ## Manufactured.Solid.Fuels.Domestic  2.012e-01  5.920e-02   3.399 0.000771 ***
    ## petroleum.Commercial               2.749e-01  9.350e-02   2.940 0.003542 ** 
    ## petroleum.Rail                     2.298e-02  1.103e-02   2.084 0.038029 *  
    ## petroleum.Public.Administration   -1.416e-01  7.305e-02  -1.939 0.053496 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1122551)
    ## 
    ##     Null deviance: 94.264  on 299  degrees of freedom
    ## Residual deviance: 32.677  on 293  degrees of freedom
    ## AIC: 1889.9
    ## 
    ## Number of Fisher Scoring iterations: 8

``` r
# Drop Public Administration because not significant
summary(model_NOx_gamma.final.rsd <- glm(formula = nox_val ~ X.people.per.sq.km + petroleum.Agriculture + 
    Manufactured.Solid.Fuels.Domestic + petroleum.Commercial + 
    petroleum.Rail, family = Gamma(link = "log"), 
    data = fossilfuel_2017.omit))
```

    ## 
    ## Call:
    ## glm(formula = nox_val ~ X.people.per.sq.km + petroleum.Agriculture + 
    ##     Manufactured.Solid.Fuels.Domestic + petroleum.Commercial + 
    ##     petroleum.Rail, family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.53944  -0.25474  -0.03266   0.17177   1.49509  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        2.652e+00  3.735e-02  71.012  < 2e-16 ***
    ## X.people.per.sq.km                 1.095e-04  7.995e-06  13.700  < 2e-16 ***
    ## petroleum.Agriculture             -9.140e-02  8.115e-03 -11.263  < 2e-16 ***
    ## Manufactured.Solid.Fuels.Domestic  1.993e-01  5.976e-02   3.335 0.000962 ***
    ## petroleum.Commercial               2.522e-01  9.351e-02   2.697 0.007390 ** 
    ## petroleum.Rail                     1.932e-02  1.094e-02   1.766 0.078511 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1143779)
    ## 
    ##     Null deviance: 94.264  on 299  degrees of freedom
    ## Residual deviance: 33.089  on 294  degrees of freedom
    ## AIC: 1891.7
    ## 
    ## Number of Fisher Scoring iterations: 8

``` r
# Drop Petroleum Rail because not significant
summary(model_NOx_gamma.final.rsd3 <- glm(formula = nox_val ~ X.people.per.sq.km + petroleum.Agriculture + 
    Manufactured.Solid.Fuels.Domestic + petroleum.Commercial, family = Gamma(link = "log"), 
    data = fossilfuel_2017.omit))
```

    ## 
    ## Call:
    ## glm(formula = nox_val ~ X.people.per.sq.km + petroleum.Agriculture + 
    ##     Manufactured.Solid.Fuels.Domestic + petroleum.Commercial, 
    ##     family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.52518  -0.26620  -0.03299   0.17772   1.46828  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        2.658e+00  3.711e-02  71.628  < 2e-16 ***
    ## X.people.per.sq.km                 1.086e-04  7.970e-06  13.623  < 2e-16 ***
    ## petroleum.Agriculture             -9.033e-02  8.114e-03 -11.133  < 2e-16 ***
    ## Manufactured.Solid.Fuels.Domestic  2.361e-01  5.621e-02   4.201 3.52e-05 ***
    ## petroleum.Commercial               2.701e-01  9.292e-02   2.907  0.00393 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1144435)
    ## 
    ##     Null deviance: 94.264  on 299  degrees of freedom
    ## Residual deviance: 33.458  on 295  degrees of freedom
    ## AIC: 1893.1
    ## 
    ## Number of Fisher Scoring iterations: 8

``` r
#Check for collinearity
vif(model_NOx_gamma.final.rsd3) # VIF <5
```

    ##                X.people.per.sq.km             petroleum.Agriculture 
    ##                          1.197093                          2.752645 
    ## Manufactured.Solid.Fuels.Domestic              petroleum.Commercial 
    ##                          2.290569                          1.884869

Last, we will look at the effect of domestic and non-domestic gas
consumption:

``` r
summary(model_NOx_gamma.gas <- glm(data = fossilfuel_2017.omit, nox_val ~  X.people.per.sq.km + Domestic.mean.gas.consumption + Non.domestic.mean.gas.consumption, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = nox_val ~ X.people.per.sq.km + Domestic.mean.gas.consumption + 
    ##     Non.domestic.mean.gas.consumption, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.65609  -0.30099  -0.02535   0.21248   1.80623  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        2.216e+00  2.346e-01   9.447   <2e-16 ***
    ## X.people.per.sq.km                 1.468e-04  9.287e-06  15.806   <2e-16 ***
    ## Domestic.mean.gas.consumption      3.194e-05  1.637e-05   1.951    0.052 .  
    ## Non.domestic.mean.gas.consumption -3.762e-08  4.239e-08  -0.888    0.376    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1643036)
    ## 
    ##     Null deviance: 94.264  on 299  degrees of freedom
    ## Residual deviance: 48.733  on 296  degrees of freedom
    ## AIC: 2006.5
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
# Non domestic gas consumption is not significant so will be dropped
summary(model_NOx_gamma.gas.v2 <- glm(data = fossilfuel_2017.omit, nox_val ~  X.people.per.sq.km + Domestic.mean.gas.consumption, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = nox_val ~ X.people.per.sq.km + Domestic.mean.gas.consumption, 
    ##     family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.66099  -0.31106  -0.02799   0.21552   1.80287  
    ## 
    ## Coefficients:
    ##                                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   2.138e+00  2.192e-01   9.751   <2e-16 ***
    ## X.people.per.sq.km            1.488e-04  9.057e-06  16.428   <2e-16 ***
    ## Domestic.mean.gas.consumption 3.556e-05  1.594e-05   2.231   0.0264 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1643312)
    ## 
    ##     Null deviance: 94.264  on 299  degrees of freedom
    ## Residual deviance: 48.855  on 297  degrees of freedom
    ## AIC: 2005.3
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
vif(model_NOx_gamma.gas.v2) # no collinearity 
```

    ##            X.people.per.sq.km Domestic.mean.gas.consumption 
    ##                       1.07663                       1.07663

``` r
#plot(model_NOx_gamma.gas.v2) # no significant outlier
```

Based on the results above, it can be concluded that domestic gas
consumption exerts a significant effect on NOx levels in England. <br>
\#\#\#\# O3

Same as above: analyse the data distribution of the target or response
variable (O3 concentration). To do this, we will use the following
functions:

``` r
hist(fossilfuel_2017.omit$o3_val, main="Histogram of Yield", xlab="Yield (quintals/ha)") 
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-75-1.png)<!-- -->

``` r
#qqnorm(fossilfuel_2017.omit$o3_val, main="QQplot of O3")
# Visualise refernce line
library("car")
#qqPlot(fossilfuel_2017.omit$o3_val,  main="QQplot of NO2")
skewness(fossilfuel_2017.omit$o3_val) 
```

    ## [1] 0.01756001

Intuitively, the skewness is a measure of symmetry. For normally
distributed data, all the data points should fall along the line. This
is clearly not the case but the deviation is not substantial. According
to Webster and Oliver (2007) if the skewness is below 0.5, we can
consider the deviation from normality not big enough to transform the
data.

Unlike data for NO2 and NOx, because the dv values of O3 values is
overdispersed and slightly skewed, we have two options: gamma glm model
or simply gaussian glm for O3 data. Our data does not appear
significantly skewed so applying a gaussian model could be more
appropriate. We will compare relative AIC values to choose the best
model.

``` r
# Gamma
summary(model_O3_gamma.vif2 <- glm(data = fossilfuel_2017.omit, o3_val ~ X.people.per.sq.km + personal.buses.A.roads + personal.motorcycles.Motorways + personal.motorcycles.A.roads + personal.motorcycles.Minor.roads + freight.HGV.Motorways + freight.HGV.A.roads, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = o3_val ~ X.people.per.sq.km + personal.buses.A.roads + 
    ##     personal.motorcycles.Motorways + personal.motorcycles.A.roads + 
    ##     personal.motorcycles.Minor.roads + freight.HGV.Motorways + 
    ##     freight.HGV.A.roads, family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -2.45610  -0.33991  -0.01237   0.23605   1.05964  
    ## 
    ## Coefficients:
    ##                                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                       1.978e+00  4.669e-02  42.356  < 2e-16 ***
    ## X.people.per.sq.km                1.210e-06  1.293e-05   0.094    0.926    
    ## personal.buses.A.roads           -3.870e-04  3.297e-05 -11.737  < 2e-16 ***
    ## personal.motorcycles.Motorways    2.260e-03  4.270e-04   5.292 2.38e-07 ***
    ## personal.motorcycles.A.roads      1.121e-03  2.026e-04   5.533 6.99e-08 ***
    ## personal.motorcycles.Minor.roads  1.392e-03  2.225e-04   6.255 1.41e-09 ***
    ## freight.HGV.Motorways            -1.612e-05  2.501e-06  -6.445 4.77e-10 ***
    ## freight.HGV.A.roads              -7.568e-07  3.680e-06  -0.206    0.837    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1386831)
    ## 
    ##     Null deviance: 79.082  on 299  degrees of freedom
    ## Residual deviance: 47.822  on 292  degrees of freedom
    ## AIC: 1476.8
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
# Gaussian glm
summary(model_O3_gaussian.vif <- glm(data = fossilfuel_2017.omit, o3_val ~ X.people.per.sq.km + personal.buses.A.roads + personal.motorcycles.Motorways + personal.motorcycles.A.roads + personal.motorcycles.Minor.roads + freight.HGV.Motorways + freight.HGV.A.roads, family = "gaussian"))
```

    ## 
    ## Call:
    ## glm(formula = o3_val ~ X.people.per.sq.km + personal.buses.A.roads + 
    ##     personal.motorcycles.Motorways + personal.motorcycles.A.roads + 
    ##     personal.motorcycles.Minor.roads + freight.HGV.Motorways + 
    ##     freight.HGV.A.roads, family = "gaussian", data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -8.2330  -2.2562  -0.1874   2.0642   7.9001  
    ## 
    ## Coefficients:
    ##                                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                       7.451e+00  3.363e-01  22.157  < 2e-16 ***
    ## X.people.per.sq.km               -1.439e-04  9.315e-05  -1.545    0.123    
    ## personal.buses.A.roads           -2.251e-03  2.375e-04  -9.478  < 2e-16 ***
    ## personal.motorcycles.Motorways    1.697e-02  3.076e-03   5.518 7.58e-08 ***
    ## personal.motorcycles.A.roads      9.844e-03  1.459e-03   6.745 8.20e-11 ***
    ## personal.motorcycles.Minor.roads  8.073e-03  1.603e-03   5.037 8.31e-07 ***
    ## freight.HGV.Motorways            -1.168e-04  1.802e-05  -6.482 3.85e-10 ***
    ## freight.HGV.A.roads              -3.214e-05  2.651e-05  -1.213    0.226    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 7.194641)
    ## 
    ##     Null deviance: 3556.1  on 299  degrees of freedom
    ## Residual deviance: 2100.8  on 292  degrees of freedom
    ## AIC: 1453.3
    ## 
    ## Number of Fisher Scoring iterations: 2

Above analysis suggests that a gaussian model (AIC 1453) performs better
than a gamma with log link (AIC 1476) based on goodness of fit.
Nonetheless we need to check that our model meets all assumptions for
linear regression:

``` r
# check residuals
#plot(model_O3_gaussian.vif)
```

We will proceed with stepwise regression on VIF-corrected variables:

``` r
### Define full and null models and do step procedure
model_O3.0_gamma <- glm(o3_val ~ 1, data = fossilfuel_2017.omit, family = "gaussian")
step(model_O3.0_gamma, direction = "both", scope = formula(model_O3_gaussian.vif, test = "Chisq"))
```

    ## Start:  AIC=1597.15
    ## o3_val ~ 1
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.buses.A.roads            1   3175.6 1565.2
    ## + X.people.per.sq.km                1   3378.9 1583.8
    ## + personal.motorcycles.Motorways    1   3490.8 1593.6
    ## + freight.HGV.Motorways             1   3491.2 1593.6
    ## <none>                                  3556.1 1597.2
    ## + freight.HGV.A.roads               1   3551.3 1598.8
    ## + personal.motorcycles.A.roads      1   3551.6 1598.8
    ## + personal.motorcycles.Minor.roads  1   3554.2 1599.0
    ## 
    ## Step:  AIC=1565.2
    ## o3_val ~ personal.buses.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.motorcycles.A.roads      1   2688.4 1517.2
    ## + personal.motorcycles.Minor.roads  1   2913.3 1541.3
    ## + freight.HGV.Motorways             1   3109.8 1560.9
    ## + personal.motorcycles.Motorways    1   3121.0 1562.0
    ## <none>                                  3175.6 1565.2
    ## + freight.HGV.A.roads               1   3155.5 1565.3
    ## + X.people.per.sq.km                1   3175.2 1567.2
    ## - personal.buses.A.roads            1   3556.1 1597.2
    ## 
    ## Step:  AIC=1517.24
    ## o3_val ~ personal.buses.A.roads + personal.motorcycles.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.motorcycles.Minor.roads  1   2474.5 1494.4
    ## + personal.motorcycles.Motorways    1   2642.6 1514.1
    ## + freight.HGV.Motorways             1   2643.7 1514.2
    ## <none>                                  2688.4 1517.2
    ## + X.people.per.sq.km                1   2672.3 1517.4
    ## + freight.HGV.A.roads               1   2681.8 1518.5
    ## - personal.motorcycles.A.roads      1   3175.6 1565.2
    ## - personal.buses.A.roads            1   3551.6 1598.8
    ## 
    ## Step:  AIC=1494.37
    ## o3_val ~ personal.buses.A.roads + personal.motorcycles.A.roads + 
    ##     personal.motorcycles.Minor.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + freight.HGV.Motorways             1   2380.2 1484.7
    ## + freight.HGV.A.roads               1   2421.0 1489.8
    ## <none>                                  2474.5 1494.4
    ## + personal.motorcycles.Motorways    1   2460.8 1494.7
    ## + X.people.per.sq.km                1   2474.4 1496.3
    ## - personal.motorcycles.Minor.roads  1   2688.4 1517.2
    ## - personal.motorcycles.A.roads      1   2913.3 1541.3
    ## - personal.buses.A.roads            1   3551.6 1600.8
    ## 
    ## Step:  AIC=1484.71
    ## o3_val ~ personal.buses.A.roads + personal.motorcycles.A.roads + 
    ##     personal.motorcycles.Minor.roads + freight.HGV.Motorways
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.motorcycles.Motorways    1   2120.8 1452.1
    ## + freight.HGV.A.roads               1   2347.5 1482.6
    ## <none>                                  2380.2 1484.7
    ## + X.people.per.sq.km                1   2374.4 1486.0
    ## - freight.HGV.Motorways             1   2474.5 1494.4
    ## - personal.motorcycles.Minor.roads  1   2643.7 1514.2
    ## - personal.motorcycles.A.roads      1   2784.1 1529.7
    ## - personal.buses.A.roads            1   3483.1 1596.9
    ## 
    ## Step:  AIC=1452.09
    ## o3_val ~ personal.buses.A.roads + personal.motorcycles.A.roads + 
    ##     personal.motorcycles.Minor.roads + freight.HGV.Motorways + 
    ##     personal.motorcycles.Motorways
    ## 
    ##                                    Df Deviance    AIC
    ## <none>                                  2120.8 1452.1
    ## + X.people.per.sq.km                1   2111.4 1452.8
    ## + freight.HGV.A.roads               1   2118.0 1453.7
    ## - personal.motorcycles.Minor.roads  1   2327.8 1478.0
    ## - personal.motorcycles.Motorways    1   2380.2 1484.7
    ## - personal.motorcycles.A.roads      1   2454.7 1494.0
    ## - freight.HGV.Motorways             1   2460.8 1494.7
    ## - personal.buses.A.roads            1   3036.3 1557.7

    ## 
    ## Call:  glm(formula = o3_val ~ personal.buses.A.roads + personal.motorcycles.A.roads + 
    ##     personal.motorcycles.Minor.roads + freight.HGV.Motorways + 
    ##     personal.motorcycles.Motorways, family = "gaussian", data = fossilfuel_2017.omit)
    ## 
    ## Coefficients:
    ##                      (Intercept)            personal.buses.A.roads  
    ##                        7.2051799                        -0.0023398  
    ##     personal.motorcycles.A.roads  personal.motorcycles.Minor.roads  
    ##                        0.0088348                         0.0081175  
    ##            freight.HGV.Motorways    personal.motorcycles.Motorways  
    ##                       -0.0001175                         0.0178090  
    ## 
    ## Degrees of Freedom: 299 Total (i.e. Null);  294 Residual
    ## Null Deviance:       3556 
    ## Residual Deviance: 2121  AIC: 1452

Based on the results of the stepwise regression and VIF test, we proceed
to build our final model:

``` r
### Define final model based on stepwise regression results
summary(model_O3_gamma_final <-glm(formula = o3_val ~ personal.buses.A.roads + personal.motorcycles.A.roads + 
    personal.motorcycles.Minor.roads + freight.HGV.Motorways + 
    personal.motorcycles.Motorways, family = "gaussian", data = fossilfuel_2017.omit))
```

    ## 
    ## Call:
    ## glm(formula = o3_val ~ personal.buses.A.roads + personal.motorcycles.A.roads + 
    ##     personal.motorcycles.Minor.roads + freight.HGV.Motorways + 
    ##     personal.motorcycles.Motorways, family = "gaussian", data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -7.8405  -2.2501  -0.1434   2.2126   7.4716  
    ## 
    ## Coefficients:
    ##                                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                       7.205e+00  3.017e-01  23.886  < 2e-16 ***
    ## personal.buses.A.roads           -2.340e-03  2.077e-04 -11.265  < 2e-16 ***
    ## personal.motorcycles.A.roads      8.835e-03  1.299e-03   6.803 5.73e-11 ***
    ## personal.motorcycles.Minor.roads  8.117e-03  1.515e-03   5.358 1.71e-07 ***
    ## freight.HGV.Motorways            -1.175e-04  1.712e-05  -6.866 3.93e-11 ***
    ## personal.motorcycles.Motorways    1.781e-02  2.970e-03   5.996 5.91e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 7.213572)
    ## 
    ##     Null deviance: 3556.1  on 299  degrees of freedom
    ## Residual deviance: 2120.8  on 294  degrees of freedom
    ## AIC: 1452.1
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
vif(model_O3_gamma_final) # VIF <5
```

    ##           personal.buses.A.roads     personal.motorcycles.A.roads 
    ##                         3.500653                         2.779741 
    ## personal.motorcycles.Minor.roads            freight.HGV.Motorways 
    ##                         1.961206                         2.293140 
    ##   personal.motorcycles.Motorways 
    ##                         2.309594

``` r
#plot(model_O3_gamma_final)
```

Construction and industrial sites:

``` r
summary(model_O3_gamma.rsd.vif <- glm(data = fossilfuel_2017.omit, o3_val ~ X.people.per.sq.km + petroleum.Industrial + petroleum.Rail + petroleum.Public.Administration + petroleum.Commercial + petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + bioenergy, family = "gaussian"))
```

    ## 
    ## Call:
    ## glm(formula = o3_val ~ X.people.per.sq.km + petroleum.Industrial + 
    ##     petroleum.Rail + petroleum.Public.Administration + petroleum.Commercial + 
    ##     petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     bioenergy, family = "gaussian", data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -8.7224  -2.8783   0.0898   2.7742   8.3195  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        8.468e+00  3.679e-01  23.018  < 2e-16 ***
    ## X.people.per.sq.km                -3.145e-04  7.813e-05  -4.025 7.27e-05 ***
    ## petroleum.Industrial              -1.241e-03  2.540e-03  -0.489   0.6255    
    ## petroleum.Rail                    -3.088e-01  1.089e-01  -2.835   0.0049 ** 
    ## petroleum.Public.Administration    1.174e+00  7.202e-01   1.630   0.1041    
    ## petroleum.Commercial               1.664e+00  9.246e-01   1.799   0.0730 .  
    ## petroleum.Agriculture              7.358e-02  7.962e-02   0.924   0.3562    
    ## Manufactured.Solid.Fuels.Domestic -1.465e-01  5.896e-01  -0.248   0.8040    
    ## bioenergy                         -2.294e-02  1.189e-02  -1.929   0.0547 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 10.88305)
    ## 
    ##     Null deviance: 3556.1  on 299  degrees of freedom
    ## Residual deviance: 3167.0  on 291  degrees of freedom
    ## AIC: 1578.4
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
vif(model_O3_gamma.rsd.vif) # no collinearity
```

    ##                X.people.per.sq.km              petroleum.Industrial 
    ##                          1.209843                          1.040351 
    ##                    petroleum.Rail   petroleum.Public.Administration 
    ##                          1.580011                          1.104873 
    ##              petroleum.Commercial             petroleum.Agriculture 
    ##                          1.962510                          2.787182 
    ## Manufactured.Solid.Fuels.Domestic                         bioenergy 
    ##                          2.650664                          1.096700

We will proceed to run a stepwise regression to drop unused variables:

``` r
model_O3.0_gamma <- glm(o3_val ~ 1, data = fossilfuel_2017.omit, family = "gaussian")
step(model_O3.0_gamma, direction = "both", scope = formula(model_O3_gamma.rsd.vif, test = "Chisq"))
```

    ## Start:  AIC=1597.15
    ## o3_val ~ 1
    ## 
    ##                                     Df Deviance    AIC
    ## + X.people.per.sq.km                 1   3378.9 1583.8
    ## + petroleum.Agriculture              1   3492.3 1593.7
    ## + petroleum.Commercial               1   3505.1 1594.8
    ## + petroleum.Public.Administration    1   3528.6 1596.8
    ## <none>                                   3556.1 1597.2
    ## + bioenergy                          1   3534.0 1597.3
    ## + Manufactured.Solid.Fuels.Domestic  1   3543.5 1598.1
    ## + petroleum.Rail                     1   3546.4 1598.3
    ## + petroleum.Industrial               1   3554.6 1599.0
    ## 
    ## Step:  AIC=1583.82
    ## o3_val ~ X.people.per.sq.km
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Commercial               1   3339.6 1582.3
    ## + bioenergy                          1   3340.2 1582.4
    ## + petroleum.Rail                     1   3341.8 1582.5
    ## + petroleum.Public.Administration    1   3355.1 1583.7
    ## <none>                                   3378.9 1583.8
    ## + petroleum.Agriculture              1   3366.0 1584.7
    ## + petroleum.Industrial               1   3372.7 1585.3
    ## + Manufactured.Solid.Fuels.Domestic  1   3378.9 1585.8
    ## - X.people.per.sq.km                 1   3556.1 1597.2
    ## 
    ## Step:  AIC=1582.31
    ## o3_val ~ X.people.per.sq.km + petroleum.Commercial
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Rail                     1   3247.7 1575.9
    ## + bioenergy                          1   3290.9 1579.9
    ## <none>                                   3339.6 1582.3
    ## + Manufactured.Solid.Fuels.Domestic  1   3318.2 1582.4
    ## + petroleum.Public.Administration    1   3327.4 1583.2
    ## + petroleum.Industrial               1   3329.2 1583.4
    ## - petroleum.Commercial               1   3378.9 1583.8
    ## + petroleum.Agriculture              1   3339.1 1584.3
    ## - X.people.per.sq.km                 1   3505.1 1594.8
    ## 
    ## Step:  AIC=1575.94
    ## o3_val ~ X.people.per.sq.km + petroleum.Commercial + petroleum.Rail
    ## 
    ##                                     Df Deviance    AIC
    ## + bioenergy                          1   3208.1 1574.2
    ## + petroleum.Public.Administration    1   3218.3 1575.2
    ## <none>                                   3247.7 1575.9
    ## + petroleum.Industrial               1   3242.1 1577.4
    ## + petroleum.Agriculture              1   3245.2 1577.7
    ## + Manufactured.Solid.Fuels.Domestic  1   3247.1 1577.9
    ## - petroleum.Rail                     1   3339.6 1582.3
    ## - petroleum.Commercial               1   3341.8 1582.5
    ## - X.people.per.sq.km                 1   3459.7 1592.9
    ## 
    ## Step:  AIC=1574.25
    ## o3_val ~ X.people.per.sq.km + petroleum.Commercial + petroleum.Rail + 
    ##     bioenergy
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Public.Administration    1   3179.5 1573.6
    ## <none>                                   3208.1 1574.2
    ## + petroleum.Agriculture              1   3200.1 1575.5
    ## + petroleum.Industrial               1   3204.3 1575.9
    ## - bioenergy                          1   3247.7 1575.9
    ## + Manufactured.Solid.Fuels.Domestic  1   3207.8 1576.2
    ## - petroleum.Rail                     1   3290.9 1579.9
    ## - petroleum.Commercial               1   3310.8 1581.7
    ## - X.people.per.sq.km                 1   3434.7 1592.7
    ## 
    ## Step:  AIC=1573.57
    ## o3_val ~ X.people.per.sq.km + petroleum.Commercial + petroleum.Rail + 
    ##     bioenergy + petroleum.Public.Administration
    ## 
    ##                                     Df Deviance    AIC
    ## <none>                                   3179.5 1573.6
    ## - petroleum.Public.Administration    1   3208.1 1574.2
    ## + petroleum.Agriculture              1   3170.3 1574.7
    ## - bioenergy                          1   3218.3 1575.2
    ## + petroleum.Industrial               1   3176.7 1575.3
    ## + Manufactured.Solid.Fuels.Domestic  1   3179.1 1575.5
    ## - petroleum.Commercial               1   3265.4 1579.6
    ## - petroleum.Rail                     1   3278.7 1580.8
    ## - X.people.per.sq.km                 1   3410.3 1592.6

    ## 
    ## Call:  glm(formula = o3_val ~ X.people.per.sq.km + petroleum.Commercial + 
    ##     petroleum.Rail + bioenergy + petroleum.Public.Administration, 
    ##     family = "gaussian", data = fossilfuel_2017.omit)
    ## 
    ## Coefficients:
    ##                     (Intercept)               X.people.per.sq.km  
    ##                       8.4374741                       -0.0003364  
    ##            petroleum.Commercial                   petroleum.Rail  
    ##                       2.0726020                       -0.3011018  
    ##                       bioenergy  petroleum.Public.Administration  
    ##                      -0.0217660                        1.1645811  
    ## 
    ## Degrees of Freedom: 299 Total (i.e. Null);  294 Residual
    ## Null Deviance:       3556 
    ## Residual Deviance: 3179  AIC: 1574

Final model for residual fuels:

``` r
### Define final model based on stepwise regression results
summary(model_O3_gamma_final.rsd <-glm(formula = o3_val ~ X.people.per.sq.km + petroleum.Commercial + 
    petroleum.Rail, 
    family = "gaussian", data = fossilfuel_2017.omit))
```

    ## 
    ## Call:
    ## glm(formula = o3_val ~ X.people.per.sq.km + petroleum.Commercial + 
    ##     petroleum.Rail, family = "gaussian", data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -9.3487  -2.9115  -0.0239   2.8478   8.3086  
    ## 
    ## Coefficients:
    ##                        Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)           8.285e+00  3.195e-01  25.928  < 2e-16 ***
    ## X.people.per.sq.km   -3.210e-04  7.304e-05  -4.395 1.54e-05 ***
    ## petroleum.Commercial  2.140e+00  7.310e-01   2.928  0.00367 ** 
    ## petroleum.Rail       -2.836e-01  9.802e-02  -2.894  0.00409 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 10.97208)
    ## 
    ##     Null deviance: 3556.1  on 299  degrees of freedom
    ## Residual deviance: 3247.7  on 296  degrees of freedom
    ## AIC: 1575.9
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
## Dropped bioenergy and Petroleum Admin - ns
vif(model_O3_gamma_final.rsd)
```

    ##   X.people.per.sq.km petroleum.Commercial       petroleum.Rail 
    ##             1.048588             1.216676             1.269279

``` r
#plot(model_O3_gamma_final.rsd)
```

Our results suggest that petroleum consumed from rail transport and
commercial activities is significantly associated with O3 AQ levels in
England.

Last, we will look at the effect of domestic and non-domestic gas
consumption:

``` r
summary(model_O3_gamma.gas <- glm(data = fossilfuel_2017.omit, o3_val ~  X.people.per.sq.km + Domestic.mean.gas.consumption + Non.domestic.mean.gas.consumption, family  = "gaussian"))
```

    ## 
    ## Call:
    ## glm(formula = o3_val ~ X.people.per.sq.km + Domestic.mean.gas.consumption + 
    ##     Non.domestic.mean.gas.consumption, family = "gaussian", data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -7.1558  -3.0658   0.5701   2.7090   9.4516  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        1.066e+01  1.899e+00   5.615 4.54e-08 ***
    ## X.people.per.sq.km                -3.514e-04  7.519e-05  -4.673 4.51e-06 ***
    ## Domestic.mean.gas.consumption     -9.308e-05  1.325e-04  -0.702    0.483    
    ## Non.domestic.mean.gas.consumption -1.443e-06  3.432e-07  -4.204 3.47e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 10.76952)
    ## 
    ##     Null deviance: 3556.1  on 299  degrees of freedom
    ## Residual deviance: 3187.8  on 296  degrees of freedom
    ## AIC: 1570.4
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
# Domestic gas consumption is not significant so will be dropped
summary(model_O3_gamma.gas.v2 <- glm(data = fossilfuel_2017.omit, o3_val ~  X.people.per.sq.km + Non.domestic.mean.gas.consumption, family  = "gaussian"))
```

    ## 
    ## Call:
    ## glm(formula = o3_val ~ X.people.per.sq.km + Non.domestic.mean.gas.consumption, 
    ##     family = "gaussian", data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##    Min      1Q  Median      3Q     Max  
    ## -7.007  -3.008   0.586   2.644   9.612  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        9.352e+00  3.439e-01  27.192  < 2e-16 ***
    ## X.people.per.sq.km                -3.353e-04  7.157e-05  -4.685 4.26e-06 ***
    ## Non.domestic.mean.gas.consumption -1.388e-06  3.338e-07  -4.157 4.22e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 10.75115)
    ## 
    ##     Null deviance: 3556.1  on 299  degrees of freedom
    ## Residual deviance: 3193.1  on 297  degrees of freedom
    ## AIC: 1568.9
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
vif(model_O3_gamma.gas.v2) # no collinearity 
```

    ##                X.people.per.sq.km Non.domestic.mean.gas.consumption 
    ##                          1.027561                          1.027561

``` r
#plot(model_O3_gamma.gas.v2) # no significant outlier
```

Results suggest that gas consumption has no effect on O3 AQ levels in
England.

#### SO2

Same as above: analyse the data distribution of the target or response
variable (SO2 concentration). To do this, we will use the following
functions:

``` r
hist(fossilfuel_2017.omit$so2_val, main="Histogram of Yield", xlab="Yield (quintals/ha)") 
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-84-1.png)<!-- -->

``` r
#qqnorm(fossilfuel_2017.omit$so2_val, main="QQplot of SO2")
skewness(fossilfuel_2017.omit$so2_val) 
```

    ## [1] 1.035942

  
Data is non-parametric and gamma distributed. Based on the results of
VIF correction we call the following model:

``` r
summary(model_SO2_gamma.vif <- glm(data = fossilfuel_2017.omit, so2_val ~ X.people.per.sq.km + personal.buses.A.roads + personal.motorcycles.Motorways + personal.motorcycles.A.roads + personal.motorcycles.Minor.roads + freight.HGV.Motorways + freight.HGV.A.roads, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = so2_val ~ X.people.per.sq.km + personal.buses.A.roads + 
    ##     personal.motorcycles.Motorways + personal.motorcycles.A.roads + 
    ##     personal.motorcycles.Minor.roads + freight.HGV.Motorways + 
    ##     freight.HGV.A.roads, family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.19260  -0.24930  -0.05693   0.14666   1.03844  
    ## 
    ## Coefficients:
    ##                                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                       2.220e-01  4.358e-02   5.093 6.32e-07 ***
    ## X.people.per.sq.km                7.621e-05  1.207e-05   6.314 1.01e-09 ***
    ## personal.buses.A.roads            1.173e-04  3.078e-05   3.811 0.000169 ***
    ## personal.motorcycles.Motorways   -3.598e-04  3.986e-04  -0.903 0.367414    
    ## personal.motorcycles.A.roads     -1.079e-03  1.891e-04  -5.706 2.84e-08 ***
    ## personal.motorcycles.Minor.roads  4.953e-05  2.077e-04   0.238 0.811664    
    ## freight.HGV.Motorways             6.387e-06  2.335e-06   2.736 0.006607 ** 
    ## freight.HGV.A.roads              -1.045e-06  3.435e-06  -0.304 0.761245    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1208095)
    ## 
    ##     Null deviance: 49.192  on 299  degrees of freedom
    ## Residual deviance: 33.464  on 292  degrees of freedom
    ## AIC: 376.41
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
vif(model_SO2_gamma.vif) # VIF <5
```

    ##               X.people.per.sq.km           personal.buses.A.roads 
    ##                         2.601107                         4.589445 
    ##   personal.motorcycles.Motorways     personal.motorcycles.A.roads 
    ##                         2.483397                         3.520010 
    ## personal.motorcycles.Minor.roads            freight.HGV.Motorways 
    ##                         2.200378                         2.546692 
    ##              freight.HGV.A.roads 
    ##                         1.654252

Stepwise regression:

``` r
### Define full and null models and do step procedure
model_SO2.0_glm.gamma <- glm(so2_val ~ 1, data = fossilfuel_2017.omit, family = Gamma(link = "log"))
step(model_SO2.0_glm.gamma, direction = "both", scope = formula(model_SO2_gamma.vif, test = "Chisq"))
```

    ## Start:  AIC=480.61
    ## so2_val ~ 1
    ## 
    ##                                    Df Deviance    AIC
    ## + X.people.per.sq.km                1   40.711 433.96
    ## + personal.buses.A.roads            1   43.295 448.78
    ## + personal.motorcycles.Minor.roads  1   47.199 471.18
    ## + freight.HGV.A.roads               1   47.565 473.27
    ## + personal.motorcycles.A.roads      1   48.659 479.55
    ## <none>                                  49.192 480.61
    ## + freight.HGV.Motorways             1   48.991 481.45
    ## + personal.motorcycles.Motorways    1   49.187 482.58
    ## 
    ## Step:  AIC=424.43
    ## so2_val ~ X.people.per.sq.km
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.motorcycles.A.roads      1   38.694 412.67
    ## + freight.HGV.Motorways             1   39.069 415.23
    ## + personal.buses.A.roads            1   40.307 423.67
    ## + freight.HGV.A.roads               1   40.375 424.14
    ## + personal.motorcycles.Motorways    1   40.416 424.42
    ## <none>                                  40.711 424.43
    ## + personal.motorcycles.Minor.roads  1   40.482 424.86
    ## - X.people.per.sq.km                1   49.192 480.24
    ## 
    ## Step:  AIC=410.84
    ## so2_val ~ X.people.per.sq.km + personal.motorcycles.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + personal.buses.A.roads            1   34.729 384.43
    ## + freight.HGV.Motorways             1   36.460 396.84
    ## + personal.motorcycles.Minor.roads  1   37.169 401.92
    ## + personal.motorcycles.Motorways    1   38.162 409.03
    ## <none>                                  38.694 410.84
    ## + freight.HGV.A.roads               1   38.672 412.69
    ## - personal.motorcycles.A.roads      1   40.711 423.31
    ## - X.people.per.sq.km                1   48.659 480.26
    ## 
    ## Step:  AIC=379.76
    ## so2_val ~ X.people.per.sq.km + personal.motorcycles.A.roads + 
    ##     personal.buses.A.roads
    ## 
    ##                                    Df Deviance    AIC
    ## + freight.HGV.Motorways             1   33.562 372.45
    ## + personal.motorcycles.Motorways    1   34.394 379.08
    ## <none>                                  34.729 379.76
    ## + personal.motorcycles.Minor.roads  1   34.704 381.56
    ## + freight.HGV.A.roads               1   34.724 381.71
    ## - personal.buses.A.roads            1   38.694 409.36
    ## - X.people.per.sq.km                1   39.285 414.08
    ## - personal.motorcycles.A.roads      1   40.307 422.23
    ## 
    ## Step:  AIC=371.31
    ## so2_val ~ X.people.per.sq.km + personal.motorcycles.A.roads + 
    ##     personal.buses.A.roads + freight.HGV.Motorways
    ## 
    ##                                    Df Deviance    AIC
    ## <none>                                  33.562 371.31
    ## + personal.motorcycles.Motorways    1   33.478 372.61
    ## + personal.motorcycles.Minor.roads  1   33.561 373.30
    ## + freight.HGV.A.roads               1   33.562 373.30
    ## - freight.HGV.Motorways             1   34.729 379.01
    ## - personal.buses.A.roads            1   36.460 393.39
    ## - personal.motorcycles.A.roads      1   38.978 414.31
    ## - X.people.per.sq.km                1   39.193 416.10

    ## 
    ## Call:  glm(formula = so2_val ~ X.people.per.sq.km + personal.motorcycles.A.roads + 
    ##     personal.buses.A.roads + freight.HGV.Motorways, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Coefficients:
    ##                  (Intercept)            X.people.per.sq.km  
    ##                    2.167e-01                     7.687e-05  
    ## personal.motorcycles.A.roads        personal.buses.A.roads  
    ##                   -1.116e-03                     1.242e-04  
    ##        freight.HGV.Motorways  
    ##                    4.858e-06  
    ## 
    ## Degrees of Freedom: 299 Total (i.e. Null);  295 Residual
    ## Null Deviance:       49.19 
    ## Residual Deviance: 33.56     AIC: 371.3

Based on the results of the stepwise regression, we proceed to build our
final model:

``` r
### Define final model based on stepwise regression results
summary(model_SO2_gamma_final <-glm(formula = so2_val ~ X.people.per.sq.km + personal.motorcycles.A.roads + 
    personal.buses.A.roads + freight.HGV.Motorways, family = Gamma(link = "log"), 
    data = fossilfuel_2017.omit))
```

    ## 
    ## Call:
    ## glm(formula = so2_val ~ X.people.per.sq.km + personal.motorcycles.A.roads + 
    ##     personal.buses.A.roads + freight.HGV.Motorways, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.18985  -0.26115  -0.05662   0.15096   1.07876  
    ## 
    ## Coefficients:
    ##                                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   2.167e-01  3.298e-02   6.570 2.28e-10 ***
    ## X.people.per.sq.km            7.687e-05  1.057e-05   7.274 3.16e-12 ***
    ## personal.motorcycles.A.roads -1.116e-03  1.677e-04  -6.651 1.41e-10 ***
    ## personal.buses.A.roads        1.242e-04  2.597e-05   4.781 2.76e-06 ***
    ## freight.HGV.Motorways         4.858e-06  1.542e-06   3.150   0.0018 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1203377)
    ## 
    ##     Null deviance: 49.192  on 299  degrees of freedom
    ## Residual deviance: 33.562  on 295  degrees of freedom
    ## AIC: 371.31
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
vif(model_SO2_gamma_final) # VIF <5
```

    ##           X.people.per.sq.km personal.motorcycles.A.roads 
    ##                     2.001236                     2.779994 
    ##       personal.buses.A.roads        freight.HGV.Motorways 
    ##                     3.280882                     1.115310

Our final model with SO2 data does indicate significant effects of
several on-road vehicles including A-road buses motorcycles and HGVs
from motorways. Let’s have a look at fossil fuel data from manufacturing
and construction site using VIF-corrected predictors:

``` r
summary(model_SO2_gamma.rsd.vif <- glm(data = fossilfuel_2017.omit, so2_val ~ X.people.per.sq.km + petroleum.Industrial + petroleum.Rail + petroleum.Public.Administration + petroleum.Commercial + petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + bioenergy, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = so2_val ~ X.people.per.sq.km + petroleum.Industrial + 
    ##     petroleum.Rail + petroleum.Public.Administration + petroleum.Commercial + 
    ##     petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     bioenergy, family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.11941  -0.21048  -0.06192   0.12808   0.89757  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        1.716e-01  3.460e-02   4.960 1.20e-06 ***
    ## X.people.per.sq.km                 4.771e-05  7.348e-06   6.492 3.64e-10 ***
    ## petroleum.Industrial               5.385e-04  2.389e-04   2.254  0.02495 *  
    ## petroleum.Rail                     3.009e-02  1.024e-02   2.937  0.00357 ** 
    ## petroleum.Public.Administration   -1.472e-01  6.773e-02  -2.173  0.03059 *  
    ## petroleum.Commercial              -2.507e-02  8.695e-02  -0.288  0.77328    
    ## petroleum.Agriculture             -7.587e-02  7.488e-03 -10.132  < 2e-16 ***
    ## Manufactured.Solid.Fuels.Domestic  3.068e-01  5.545e-02   5.532 7.07e-08 ***
    ## bioenergy                         -5.978e-04  1.118e-03  -0.535  0.59337    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.09625709)
    ## 
    ##     Null deviance: 49.192  on 299  degrees of freedom
    ## Residual deviance: 26.529  on 291  degrees of freedom
    ## AIC: 307.59
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
vif(model_SO2_gamma.rsd.vif) # no collinearity
```

    ##                X.people.per.sq.km              petroleum.Industrial 
    ##                          1.209843                          1.040351 
    ##                    petroleum.Rail   petroleum.Public.Administration 
    ##                          1.580011                          1.104873 
    ##              petroleum.Commercial             petroleum.Agriculture 
    ##                          1.962510                          2.787182 
    ## Manufactured.Solid.Fuels.Domestic                         bioenergy 
    ##                          2.650664                          1.096700

Stepwise regression:

``` r
### Define full and null models and do step procedure
model_SO2.0_glm.gamma <- glm(so2_val ~ 1, data = fossilfuel_2017.omit, family = Gamma(link = "log"))
step(model_SO2.0_glm.gamma, direction = "both", scope = formula(model_SO2_gamma.rsd.vif, test = "Chisq"))
```

    ## Start:  AIC=480.61
    ## so2_val ~ 1
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Agriculture              1   37.004 412.69
    ## + X.people.per.sq.km                 1   40.711 433.96
    ## + petroleum.Commercial               1   47.810 474.68
    ## + Manufactured.Solid.Fuels.Domestic  1   48.456 478.38
    ## + petroleum.Public.Administration    1   48.692 479.74
    ## <none>                                   49.192 480.61
    ## + petroleum.Industrial               1   48.952 481.23
    ## + bioenergy                          1   48.959 481.27
    ## + petroleum.Rail                     1   49.167 482.46
    ## 
    ## Step:  AIC=395.17
    ## so2_val ~ petroleum.Agriculture
    ## 
    ##                                     Df Deviance    AIC
    ## + Manufactured.Solid.Fuels.Domestic  1   32.110 361.50
    ## + X.people.per.sq.km                 1   33.474 371.45
    ## + petroleum.Rail                     1   34.546 379.26
    ## + petroleum.Commercial               1   36.121 390.73
    ## + petroleum.Industrial               1   36.142 390.88
    ## <none>                                   37.004 395.17
    ## + bioenergy                          1   36.904 396.44
    ## + petroleum.Public.Administration    1   36.970 396.92
    ## - petroleum.Agriculture              1   49.192 481.98
    ## 
    ## Step:  AIC=353.79
    ## so2_val ~ petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic
    ## 
    ##                                     Df Deviance    AIC
    ## + X.people.per.sq.km                 1   28.431 325.12
    ## + petroleum.Rail                     1   31.517 350.85
    ## + petroleum.Industrial               1   31.526 350.93
    ## <none>                                   32.110 353.79
    ## + petroleum.Commercial               1   31.885 353.92
    ## + petroleum.Public.Administration    1   31.913 354.16
    ## + bioenergy                          1   32.095 355.67
    ## - Manufactured.Solid.Fuels.Domestic  1   37.004 392.60
    ## - petroleum.Agriculture              1   48.456 488.09
    ## 
    ## Step:  AIC=318.68
    ## so2_val ~ petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     X.people.per.sq.km
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Rail                     1   27.591 312.71
    ## + petroleum.Industrial               1   27.626 313.05
    ## + petroleum.Public.Administration    1   28.137 317.89
    ## <none>                                   28.431 318.68
    ## + petroleum.Commercial               1   28.430 320.67
    ## + bioenergy                          1   28.430 320.67
    ## - X.people.per.sq.km                 1   32.110 351.57
    ## - Manufactured.Solid.Fuels.Domestic  1   33.474 364.52
    ## - petroleum.Agriculture              1   40.697 433.03
    ## 
    ## Step:  AIC=311.54
    ## so2_val ~ petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     X.people.per.sq.km + petroleum.Rail
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Industrial               1   26.982 307.53
    ## + petroleum.Public.Administration    1   27.115 308.84
    ## <none>                                   27.591 311.54
    ## + petroleum.Commercial               1   27.572 313.35
    ## + bioenergy                          1   27.590 313.53
    ## - petroleum.Rail                     1   28.431 317.83
    ## - Manufactured.Solid.Fuels.Domestic  1   30.543 338.69
    ## - X.people.per.sq.km                 1   31.517 348.31
    ## - petroleum.Agriculture              1   40.260 434.64
    ## 
    ## Step:  AIC=306.74
    ## so2_val ~ petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     X.people.per.sq.km + petroleum.Rail + petroleum.Industrial
    ## 
    ##                                     Df Deviance    AIC
    ## + petroleum.Public.Administration    1   26.562 304.43
    ## <none>                                   26.982 306.74
    ## + petroleum.Commercial               1   26.954 308.46
    ## + bioenergy                          1   26.963 308.54
    ## - petroleum.Industrial               1   27.591 311.00
    ## - petroleum.Rail                     1   27.626 311.36
    ## - Manufactured.Solid.Fuels.Domestic  1   29.876 334.45
    ## - X.people.per.sq.km                 1   31.075 346.76
    ## - petroleum.Agriculture              1   39.708 435.37
    ## 
    ## Step:  AIC=303.96
    ## so2_val ~ petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     X.people.per.sq.km + petroleum.Rail + petroleum.Industrial + 
    ##     petroleum.Public.Administration
    ## 
    ##                                     Df Deviance    AIC
    ## <none>                                   26.562 303.96
    ## + bioenergy                          1   26.538 305.72
    ## + petroleum.Commercial               1   26.554 305.89
    ## - petroleum.Public.Administration    1   26.982 306.35
    ## - petroleum.Industrial               1   27.115 307.74
    ## - petroleum.Rail                     1   27.364 310.33
    ## - Manufactured.Solid.Fuels.Domestic  1   29.569 333.35
    ## - X.people.per.sq.km                 1   30.803 346.22
    ## - petroleum.Agriculture              1   39.100 432.80

    ## 
    ## Call:  glm(formula = so2_val ~ petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     X.people.per.sq.km + petroleum.Rail + petroleum.Industrial + 
    ##     petroleum.Public.Administration, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Coefficients:
    ##                       (Intercept)              petroleum.Agriculture  
    ##                         1.683e-01                         -7.684e-02  
    ## Manufactured.Solid.Fuels.Domestic                 X.people.per.sq.km  
    ##                         3.005e-01                          4.724e-05  
    ##                    petroleum.Rail               petroleum.Industrial  
    ##                         2.986e-02                          5.232e-04  
    ##   petroleum.Public.Administration  
    ##                        -1.493e-01  
    ## 
    ## Degrees of Freedom: 299 Total (i.e. Null);  293 Residual
    ## Null Deviance:       49.19 
    ## Residual Deviance: 26.56     AIC: 304

Based on the results of the stepwise regression on the original
VIF-corrected model, we build our final model:

``` r
### Define final model based on stepwise regression results
summary(model_SO2_gamma.final.rsd <- glm(formula = so2_val ~ petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    X.people.per.sq.km + petroleum.Rail + petroleum.Industrial + 
    petroleum.Public.Administration, family = Gamma(link = "log"), 
    data = fossilfuel_2017.omit))
```

    ## 
    ## Call:
    ## glm(formula = so2_val ~ petroleum.Agriculture + Manufactured.Solid.Fuels.Domestic + 
    ##     X.people.per.sq.km + petroleum.Rail + petroleum.Industrial + 
    ##     petroleum.Public.Administration, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.11650  -0.21117  -0.06106   0.12561   0.90408  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        1.683e-01  3.409e-02   4.937 1.34e-06 ***
    ## petroleum.Agriculture             -7.684e-02  6.775e-03 -11.342  < 2e-16 ***
    ## Manufactured.Solid.Fuels.Domestic  3.005e-01  5.381e-02   5.585 5.35e-08 ***
    ## X.people.per.sq.km                 4.724e-05  7.134e-06   6.622 1.69e-10 ***
    ## petroleum.Rail                     2.986e-02  1.018e-02   2.933  0.00362 ** 
    ## petroleum.Industrial               5.232e-04  2.376e-04   2.202  0.02843 *  
    ## petroleum.Public.Administration   -1.493e-01  6.694e-02  -2.231  0.02647 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.09582971)
    ## 
    ##     Null deviance: 49.192  on 299  degrees of freedom
    ## Residual deviance: 26.562  on 293  degrees of freedom
    ## AIC: 303.96
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
#Check for collinearity
vif(model_SO2_gamma.final.rsd) # VIF <5
```

    ##             petroleum.Agriculture Manufactured.Solid.Fuels.Domestic 
    ##                          2.291896                          2.507430 
    ##                X.people.per.sq.km                    petroleum.Rail 
    ##                          1.145312                          1.567892 
    ##              petroleum.Industrial   petroleum.Public.Administration 
    ##                          1.033574                          1.083978

Data suggests important contributions of petroleum consumption from rail
network and industrial activities as well as manufactured solid fuels
used for domestic purposes. Last, we will look at the effect of domestic
and non-domestic gas consumption:

``` r
summary(model_SO2_gamma.gas <- glm(data = fossilfuel_2017.omit, so2_val ~  X.people.per.sq.km + Domestic.mean.gas.consumption + Non.domestic.mean.gas.consumption, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = so2_val ~ X.people.per.sq.km + Domestic.mean.gas.consumption + 
    ##     Non.domestic.mean.gas.consumption, family = Gamma(link = "log"), 
    ##     data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.13774  -0.27879  -0.07583   0.13735   1.04062  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        3.531e-01  2.158e-01   1.636   0.1028    
    ## X.people.per.sq.km                 6.716e-05  8.544e-06   7.861 7.14e-14 ***
    ## Domestic.mean.gas.consumption     -1.570e-05  1.506e-05  -1.042   0.2981    
    ## Non.domestic.mean.gas.consumption  9.314e-08  3.900e-08   2.388   0.0176 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1390659)
    ## 
    ##     Null deviance: 49.192  on 299  degrees of freedom
    ## Residual deviance: 39.643  on 296  degrees of freedom
    ## AIC: 420.28
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
# Domestic gas consumption is not significant so will be dropped
summary(model_SO2_gamma.gas.v2 <- glm(data = fossilfuel_2017.omit, so2_val ~  X.people.per.sq.km + Non.domestic.mean.gas.consumption, family  = Gamma(link = "log")))
```

    ## 
    ## Call:
    ## glm(formula = so2_val ~ X.people.per.sq.km + Non.domestic.mean.gas.consumption, 
    ##     family = Gamma(link = "log"), data = fossilfuel_2017.omit)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.16453  -0.28096  -0.08361   0.13715   1.03577  
    ## 
    ## Coefficients:
    ##                                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                       1.304e-01  3.936e-02   3.312  0.00104 ** 
    ## X.people.per.sq.km                7.023e-05  8.191e-06   8.574  5.6e-16 ***
    ## Non.domestic.mean.gas.consumption 1.040e-07  3.821e-08   2.722  0.00688 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1408283)
    ## 
    ##     Null deviance: 49.192  on 299  degrees of freedom
    ## Residual deviance: 39.782  on 297  degrees of freedom
    ## AIC: 419.35
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
vif(model_SO2_gamma.gas.v2) # no collinearity 
```

    ##                X.people.per.sq.km Non.domestic.mean.gas.consumption 
    ##                          1.027561                          1.027561

Based on the results above, it can be concluded that non-domestic gas
consumption exerts a significant effect on SO2 levels in England.

### Results output

Now we can generate tables using the following formula:

#### Supplementary Table 10. Effect of fossil fuel consumption sources annual mean nitrogen dioxide AQ levels at local authority level

<br>

``` r
library(stargazer)
stargazer(model_NO2_gamma.final.v2, model_NO2_gamma.final.rsd, model_NO2_gamma.gas.v2, type ="html", single.row=TRUE, out = "fig_out_v4/Table_ff_NO2.html")
```

<table style="text-align:center">

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="3">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="3">

no2\_val

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

(1)

</td>

<td>

(2)

</td>

<td>

(3)

</td>

</tr>

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

X.people.per.sq.km

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00001)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00001)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00001)

</td>

</tr>

<tr>

<td style="text-align:left">

personal.motorcycles.Motorways

</td>

<td>

0.001<sup>\*\*\*</sup> (0.0003)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

personal.motorcycles.A.roads

</td>

<td>

\-0.001<sup>\*\*\*</sup> (0.0002)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

personal.buses.A.roads

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00003)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

petroleum.Agriculture

</td>

<td>

</td>

<td>

\-0.083<sup>\*\*\*</sup> (0.007)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

petroleum.Commercial

</td>

<td>

</td>

<td>

0.183<sup>\*\*</sup> (0.081)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

Manufactured.Solid.Fuels.Domestic

</td>

<td>

</td>

<td>

0.191<sup>\*\*\*</sup> (0.052)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

petroleum.Rail

</td>

<td>

</td>

<td>

0.019<sup>\*\*</sup> (0.009)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

Domestic.mean.gas.consumption

</td>

<td>

</td>

<td>

</td>

<td>

0.00004<sup>\*\*\*</sup> (0.00001)

</td>

</tr>

<tr>

<td style="text-align:left">

Constant

</td>

<td>

2.298<sup>\*\*\*</sup> (0.032)

</td>

<td>

2.370<sup>\*\*\*</sup> (0.032)

</td>

<td>

1.847<sup>\*\*\*</sup> (0.189)

</td>

</tr>

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Observations

</td>

<td>

300

</td>

<td>

300

</td>

<td>

300

</td>

</tr>

<tr>

<td style="text-align:left">

Log Likelihood

</td>

<td>

\-855.220

</td>

<td>

\-810.644

</td>

<td>

\-873.052

</td>

</tr>

<tr>

<td style="text-align:left">

Akaike Inf. Crit.

</td>

<td>

1,720.440

</td>

<td>

1,633.288

</td>

<td>

1,752.104

</td>

</tr>

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="3" style="text-align:right">

<sup>*</sup>p\<0.1; <sup>**</sup>p\<0.05; <sup>***</sup>p\<0.01

</td>

</tr>

</table>

<br> Summary of reported contributions of fossil fuel consumption
sources on annual mean nitrogen dioxide concentration. Each column
represents an individual single-pollutant generalized linear model
(gamma) used to assess the effects of road transport fuels (1), residual
fuels (2) and gas consumption (3) on annual average concentration of
nitrogen dioxide. Each value corresponds to the estimate of the model,
with the standard error in parentheses. The p-values are indicated using
the number of asterisks beside the estimates. glm, generalized linear
model; no2\_val, annual mean nitrogen dioxide values; Akaike Inf. Crit.,
Akaike’s Information Criteria.

``` r
# Calculate odds ratio (OR) for plotting 
# Format OR
NO2_rd_odds = data.frame(cbind(exp(cbind(OR = coef(model_NO2_gamma.final.v2), confint(model_NO2_gamma.final.v2))), p_value = summary(model_NO2_gamma.final.v2)$coefficients[,4]))
NO2_rsd_odds = data.frame(cbind(exp(cbind(OR = coef(model_NO2_gamma.final.rsd), confint(model_NO2_gamma.final.rsd))), p_value = summary(model_NO2_gamma.final.rsd)$coefficients[,4]))
NO2_gas_odds = data.frame(cbind(exp(cbind(OR = coef(model_NO2_gamma.gas.v2), confint(model_NO2_gamma.gas.v2))), p_value = summary(model_NO2_gamma.gas.v2)$coefficients[,4]))
#now delete the rows that you dont want  - Intercepts
NO2_rd_odds_clean = NO2_rd_odds[-c(1,2),]
NO2_gas_odds_clean = NO2_gas_odds[-c(1,2),]
NO2_rsd_odds_clean = NO2_rsd_odds[-c(1,2),]

#make a df just with the pollutants 
NO2_ff_onlyPoll = data.frame(rbind(NO2_rd_odds_clean, NO2_gas_odds_clean, NO2_rsd_odds_clean))

# Sort data
NO2_ff_onlyPoll$names = row.names(NO2_ff_onlyPoll)
NO2_ff_onlyPoll$significance = "p-value > 0.05"
NO2_ff_onlyPoll$significance[NO2_ff_onlyPoll$p_value < 0.05] <- "p-value < 0.05"

#subset
NO2_ff_onlyPoll <- subset(NO2_ff_onlyPoll, OR >= 1.01 | OR <= 0.99)
# plot
ggplot(NO2_ff_onlyPoll, aes(x=reorder(names, OR), y=OR, color= significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_bw() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("AQ odds ratios") + 
  xlab("NO2")+
  ylim(0.75, 1.6)
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-93-1.png)<!-- -->

``` r
ggsave("fig_out_v4/NO2_sources.pdf")
```

#### Supplementary Table 11. Effect of fossil fuel consumption sources on annual mean nitrogen oxides AQ levels at local authority level

Add <br>

``` r
stargazer(model_NOx_gamma_final, model_NOx_gamma.final.rsd3, model_NOx_gamma.gas.v2, type="html",ci=TRUE, ci.level=0.95, single.row=TRUE, out = "fig_out_v4/Table_ff_NOx.html")
```

<table style="text-align:center">

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="3">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="3">

nox\_val

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

(1)

</td>

<td>

(2)

</td>

<td>

(3)

</td>

</tr>

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

X.people.per.sq.km

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.0001, 0.0002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.0001, 0.0001)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.0001, 0.0002)

</td>

</tr>

<tr>

<td style="text-align:left">

personal.motorcycles.Motorways

</td>

<td>

0.001<sup>\*\*\*</sup> (0.001, 0.002)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

personal.buses.A.roads

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00005, 0.0002)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

personal.motorcycles.A.roads

</td>

<td>

\-0.001<sup>\*\*\*</sup> (-0.001, -0.0004)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

petroleum.Agriculture

</td>

<td>

</td>

<td>

\-0.090<sup>\*\*\*</sup> (-0.106, -0.074)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

Manufactured.Solid.Fuels.Domestic

</td>

<td>

</td>

<td>

0.236<sup>\*\*\*</sup> (0.126, 0.346)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

petroleum.Commercial

</td>

<td>

</td>

<td>

0.270<sup>\*\*\*</sup> (0.088, 0.452)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

Domestic.mean.gas.consumption

</td>

<td>

</td>

<td>

</td>

<td>

0.00004<sup>\*\*</sup> (0.00000, 0.0001)

</td>

</tr>

<tr>

<td style="text-align:left">

Constant

</td>

<td>

2.571<sup>\*\*\*</sup> (2.496, 2.646)

</td>

<td>

2.658<sup>\*\*\*</sup> (2.586, 2.731)

</td>

<td>

2.138<sup>\*\*\*</sup> (1.708, 2.567)

</td>

</tr>

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Observations

</td>

<td>

300

</td>

<td>

300

</td>

<td>

300

</td>

</tr>

<tr>

<td style="text-align:left">

Log Likelihood

</td>

<td>

\-982.874

</td>

<td>

\-941.562

</td>

<td>

\-999.629

</td>

</tr>

<tr>

<td style="text-align:left">

Akaike Inf. Crit.

</td>

<td>

1,975.747

</td>

<td>

1,893.123

</td>

<td>

2,005.257

</td>

</tr>

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="3" style="text-align:right">

<sup>*</sup>p\<0.1; <sup>**</sup>p\<0.05; <sup>***</sup>p\<0.01

</td>

</tr>

</table>

<br> Summary of reported contributions of fossil fuel consumption
sources on annual mean nitrogen oxides concentration. Each column
represents an individual single-pollutant generalized linear model
(gamma) used to assess the effects of road transport fuels (1), residual
fuels (2) and gas consumption (3) on annual average concentration of
nitrogen oxide. Each value in the columns corresponds to the estimate of
the model, with the standard error in parentheses. The p-values are
indicated using the number of asterisks beside the estimates. glm,
generalized linear model; nox\_val, annual mean nitrogen oxide values;
Akaike Inf. Crit., Akaike’s Information Criteria.

OR graphs with individual variables:

``` r
# Calculate odds ratio (OR) for plotting 
# Format OR
NOx_rd_odds = data.frame(cbind(exp(cbind(OR = coef(model_NOx_gamma_final), confint(model_NOx_gamma_final))), p_value = summary(model_NOx_gamma_final)$coefficients[,4]))
NOx_rsd_odds = data.frame(cbind(exp(cbind(OR = coef(model_NOx_gamma.final.rsd3), confint(model_NOx_gamma.final.rsd3))), p_value = summary(model_NOx_gamma.final.rsd3)$coefficients[,4]))
NOx_gas_odds = data.frame(cbind(exp(cbind(OR = coef(model_NOx_gamma.gas.v2), confint(model_NOx_gamma.gas.v2))), p_value = summary(model_NOx_gamma.gas.v2)$coefficients[,4]))

#now delete the rows that you dont want - Intercepts
NOx_rd_odds_clean = NOx_rd_odds[-c(1,2),]
NOx_gas_odds_clean = NOx_gas_odds[-c(1,2),]
NOx_rsd_odds_clean = NOx_rsd_odds[-c(1,2),]

#make a df just with the pollutants 
NOx_ff_onlyPoll = data.frame(rbind(NOx_rd_odds_clean, NOx_gas_odds_clean, NOx_rsd_odds_clean))

# Sort data
NOx_ff_onlyPoll$names = row.names(NOx_ff_onlyPoll)
NOx_ff_onlyPoll$significance = "p-value > 0.05"
NOx_ff_onlyPoll$significance[NOx_ff_onlyPoll$p_value < 0.05] <- "p-value < 0.05"

#subset
NOx_ff_onlyPoll <- subset(NOx_ff_onlyPoll, OR >= 1.01 | OR < 0.99)
                  
# plot
ggplot(NOx_ff_onlyPoll, aes(x=reorder(names, OR), y=OR, color= significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_bw() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("AQ odds ratios") + 
  xlab("NOx")+
  ylim(0.75, 1.6)
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-95-1.png)<!-- -->

``` r
ggsave("fig_out_v4/NOx_sources.pdf")
```

#### Supplementary Table 12. Effect of fossil fuel consumption sources on annual mean ozone AQ levels at local authority level

Add <br>

``` r
stargazer(model_O3_gamma_final, model_O3_gamma_final.rsd, model_O3_gamma.gas.v2, type="html",ci=TRUE, ci.level=0.95, single.row=TRUE, out = "fig_out_v4/Table_ff_O3.html")
```

<table style="text-align:center">

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="3">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="3">

o3\_val

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

(1)

</td>

<td>

(2)

</td>

<td>

(3)

</td>

</tr>

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

personal.buses.A.roads

</td>

<td>

\-0.002<sup>\*\*\*</sup> (-0.003, -0.002)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

personal.motorcycles.A.roads

</td>

<td>

0.009<sup>\*\*\*</sup> (0.006, 0.011)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

personal.motorcycles.Minor.roads

</td>

<td>

0.008<sup>\*\*\*</sup> (0.005, 0.011)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

freight.HGV.Motorways

</td>

<td>

\-0.0001<sup>\*\*\*</sup> (-0.0002, -0.0001)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

personal.motorcycles.Motorways

</td>

<td>

0.018<sup>\*\*\*</sup> (0.012, 0.024)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

X.people.per.sq.km

</td>

<td>

</td>

<td>

\-0.0003<sup>\*\*\*</sup> (-0.0005, -0.0002)

</td>

<td>

\-0.0003<sup>\*\*\*</sup> (-0.0005, -0.0002)

</td>

</tr>

<tr>

<td style="text-align:left">

petroleum.Commercial

</td>

<td>

</td>

<td>

2.140<sup>\*\*\*</sup> (0.708, 3.573)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

petroleum.Rail

</td>

<td>

</td>

<td>

\-0.284<sup>\*\*\*</sup> (-0.476, -0.092)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

Non.domestic.mean.gas.consumption

</td>

<td>

</td>

<td>

</td>

<td>

\-0.00000<sup>\*\*\*</sup> (-0.00000, -0.00000)

</td>

</tr>

<tr>

<td style="text-align:left">

Constant

</td>

<td>

7.205<sup>\*\*\*</sup> (6.614, 7.796)

</td>

<td>

8.285<sup>\*\*\*</sup> (7.658, 8.911)

</td>

<td>

9.352<sup>\*\*\*</sup> (8.678, 10.026)

</td>

</tr>

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Observations

</td>

<td>

300

</td>

<td>

300

</td>

<td>

300

</td>

</tr>

<tr>

<td style="text-align:left">

Log Likelihood

</td>

<td>

\-720.046

</td>

<td>

\-783.971

</td>

<td>

\-781.426

</td>

</tr>

<tr>

<td style="text-align:left">

Akaike Inf. Crit.

</td>

<td>

1,452.092

</td>

<td>

1,575.942

</td>

<td>

1,568.852

</td>

</tr>

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="3" style="text-align:right">

<sup>*</sup>p\<0.1; <sup>**</sup>p\<0.05; <sup>***</sup>p\<0.01

</td>

</tr>

</table>

<br> Summary of reported contributions of fossil fuel consumption
sources on annual mean ozone concentration. Each column represents an
individual single-pollutant generalized linear model (gaussian) used to
assess the effects of road transport fuels (1), residual fuels (2) and
gas consumption (3) on annual average concentration of ozone. Each value
corresponds to the estimate of the model, with the standard error in
parentheses. The p-values are indicated using the number of asterisks
beside the estimates. glm, generalized linear model; o3\_val, annual
mean ozone values; Akaike Inf. Crit., Akaike’s Information Criteria.

``` r
# Calculate odds ratio (OR) for plotting 
# Format OR
O3_rd_odds = data.frame(cbind(exp(cbind(OR = coef(model_O3_gamma_final), confint(model_O3_gamma_final))), p_value = summary(model_O3_gamma_final)$coefficients[,4]))
O3_rsd_odds = data.frame(cbind(exp(cbind(OR = coef(model_O3_gamma_final.rsd), confint(model_O3_gamma_final.rsd))), p_value = summary(model_O3_gamma_final.rsd)$coefficients[,4]))
O3_gas_odds = data.frame(cbind(exp(cbind(OR = coef(model_O3_gamma.gas.v2), confint(model_O3_gamma.gas.v2))), p_value = summary(model_O3_gamma.gas.v2)$coefficients[,4]))

#now delete the rows that you dont want - Intercepts
O3_rd_odds_clean = O3_rd_odds[-c(1,2),]
O3_gas_odds_clean = O3_gas_odds[-c(1,2),]
O3_rsd_odds_clean = O3_rsd_odds[-c(1,2),]

#make a df just with the pollutants 
O3_ff_onlyPoll = data.frame(rbind(O3_rd_odds_clean, O3_gas_odds_clean, O3_rsd_odds_clean))

# Sort data
O3_ff_onlyPoll$names = row.names(O3_ff_onlyPoll)
O3_ff_onlyPoll$significance = "p-value > 0.05"
O3_ff_onlyPoll$significance[O3_ff_onlyPoll$p_value < 0.05] <- "p-value < 0.05"

#subset
O3_ff_onlyPoll <- subset(O3_ff_onlyPoll, OR >= 1.01 | OR <= 0.99)
# plot
ggplot(O3_ff_onlyPoll, aes(x=reorder(names, OR), y=OR, color= significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_bw() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("AQ odds ratios") + 
  xlab("O3")
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-97-1.png)<!-- -->

``` r
ggsave("fig_out_v4/O3_sources.pdf")
```

#### Supplementary Table 13. Effect of fossil fuel consumption sources on annual mean SO2 AQ levels at local authority level

Add <br>

``` r
library(stargazer)
stargazer(model_SO2_gamma_final, model_SO2_gamma.final.rsd, model_SO2_gamma.gas.v2, type ="html", ci=TRUE, ci.level=0.95, single.row=TRUE, out = "fig_out_v4/Table_ff_SO2.html")
```

<table style="text-align:center">

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="3">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="3" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="3">

so2\_val

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td>

(1)

</td>

<td>

(2)

</td>

<td>

(3)

</td>

</tr>

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

petroleum.Agriculture

</td>

<td>

</td>

<td>

\-0.077<sup>\*\*\*</sup> (-0.090, -0.064)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

Manufactured.Solid.Fuels.Domestic

</td>

<td>

</td>

<td>

0.301<sup>\*\*\*</sup> (0.195, 0.406)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

X.people.per.sq.km

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.0001, 0.0001)

</td>

<td>

0.00005<sup>\*\*\*</sup> (0.00003, 0.0001)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.0001, 0.0001)

</td>

</tr>

<tr>

<td style="text-align:left">

personal.motorcycles.A.roads

</td>

<td>

\-0.001<sup>\*\*\*</sup> (-0.001, -0.001)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

personal.buses.A.roads

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.0001, 0.0002)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

freight.HGV.Motorways

</td>

<td>

0.00000<sup>\*\*\*</sup> (0.00000, 0.00001)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

petroleum.Rail

</td>

<td>

</td>

<td>

0.030<sup>\*\*\*</sup> (0.010, 0.050)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

petroleum.Industrial

</td>

<td>

</td>

<td>

0.001<sup>\*\*</sup> (0.0001, 0.001)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

petroleum.Public.Administration

</td>

<td>

</td>

<td>

\-0.149<sup>\*\*</sup> (-0.281, -0.018)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

Non.domestic.mean.gas.consumption

</td>

<td>

</td>

<td>

</td>

<td>

0.00000<sup>\*\*\*</sup> (0.00000, 0.00000)

</td>

</tr>

<tr>

<td style="text-align:left">

Constant

</td>

<td>

0.217<sup>\*\*\*</sup> (0.152, 0.281)

</td>

<td>

0.168<sup>\*\*\*</sup> (0.101, 0.235)

</td>

<td>

0.130<sup>\*\*\*</sup> (0.053, 0.208)

</td>

</tr>

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Observations

</td>

<td>

300

</td>

<td>

300

</td>

<td>

300

</td>

</tr>

<tr>

<td style="text-align:left">

Log Likelihood

</td>

<td>

\-180.654

</td>

<td>

\-144.982

</td>

<td>

\-206.674

</td>

</tr>

<tr>

<td style="text-align:left">

Akaike Inf. Crit.

</td>

<td>

371.308

</td>

<td>

303.964

</td>

<td>

419.348

</td>

</tr>

<tr>

<td colspan="4" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="3" style="text-align:right">

<sup>*</sup>p\<0.1; <sup>**</sup>p\<0.05; <sup>***</sup>p\<0.01

</td>

</tr>

</table>

<br> Summary of reported contributions of fossil fuel consumption
sources on annual mean SO2 concentration. Each column represents an
individual single-pollutant generalized linear model (gamma) used to
assess the effects of road transport fuels (1), residual fuels (2) and
gas consumption (3) on annual average concentration of SO2. The p-values
are indicated using the number of asterisks beside the estimates. glm,
generalized linear model; so2\_val, annual mean sulphur dioxide values;
Akaike Inf. Crit., Akaike’s Information Criteria.

### Calculate the odds ratios

``` r
# Format OR
SO2_rd_odds = data.frame(cbind(exp(cbind(OR = coef(model_SO2_gamma_final), confint(model_SO2_gamma_final))), p_value = summary(model_SO2_gamma_final)$coefficients[,4]))
SO2_rsd_odds = data.frame(cbind(exp(cbind(OR = coef(model_SO2_gamma.final.rsd), confint(model_SO2_gamma.final.rsd))), p_value = summary(model_SO2_gamma.final.rsd)$coefficients[,4]))
SO2_gas_odds = data.frame(cbind(exp(cbind(OR = coef(model_SO2_gamma.gas.v2), confint(model_SO2_gamma.gas.v2))), p_value = summary(model_SO2_gamma.gas.v2)$coefficients[,4]))
#now delete the rows that you dont want  - Intercepts
SO2_rd_odds_clean = SO2_rd_odds[-c(1,2),]
SO2_gas_odds_clean = SO2_gas_odds[-c(1,2),]
SO2_rsd_odds_clean = SO2_rsd_odds[-c(1,2),]

#make a df just with the pollutants 
SO2_ff_onlyPoll = data.frame(rbind(SO2_rd_odds_clean, SO2_gas_odds_clean, SO2_rsd_odds_clean))

# Sort data
SO2_ff_onlyPoll$names = row.names(SO2_ff_onlyPoll)
SO2_ff_onlyPoll$significance = "p-value > 0.05"
SO2_ff_onlyPoll$significance[NO2_ff_onlyPoll$p_value < 0.05] <- "p-value < 0.05"

#take out incorrect variables 
SO2_ff_onlyPoll <- subset(SO2_ff_onlyPoll, OR >= 1.01 | OR <= 0.99)

# plot
ggplot(SO2_ff_onlyPoll, aes(x=reorder(names, OR), y=OR, color= significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_bw() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("AQ odds ratios") + 
  xlab("SO2") +
  ylim(0.75, 1.6)
```

![](Analysis_Workflow_complete_COVID_air_notebook_v04_files/figure-gfm/unnamed-chunk-99-1.png)<!-- -->

``` r
ggsave("fig_out_v4/SO2_sources.pdf")
```
