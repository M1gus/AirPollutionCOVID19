Travaglio et al., 2020
================
Yizhou Yu, MRC Toxicology Unit, University of
Cambridge, (<yzy21@mrc-tox.cam.ac.uk>)
September 2020

# General information

This pipeline relates to data included in **Figures 3 to 5**
of our manuscript titled **Links between air pollution and COVID-19 in
England**, as well as Supplementary tables included in the supplementary
materials.

The input files for this analysis pipeline are on the master branch of
this GitHub page (link: <https://github.com/M1gus/AirPollutionCOVID19>)

The main aims of the workflow presented here were as follows: <br> 1.
Determine a relationship between air pollutants and COVID-19-associated
deaths/cases in England<br> 2. Investigate whether any relationship
between air pollution and COVID-19 remains significant in the presence
of confounding factors at the regional and subregional level<br> 3.
Determine the effect of air pollutants on infectivity at individual
levels. <br>

Only the UK Biobank data is not available in the repository as they
require separate application. Please visit ukbiobank.ac.uk for more
information. A detailed list of the variables used in the UK Biobank
analysis is available here:

#### Supplementary Table 1. Variables from the UK Biobank.

``` r
library(stargazer)
stargazer(read.csv("data/suppl_table_1.csv"), summary=FALSE, rownames=FALSE, type = "html")
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
libopenblas-dev<br> texlive-fonts-recommended<br> lmodern<br> Can be
done with apt-get on ubuntu

### View the data

Read and verify the data:

``` r
preL_dt = read.csv("data/26-4-2020_yyAIR_COVID_PRE_LD_dt.csv")[,-1]
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

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

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
pop_dens = read.csv("data/2018_official_popDensity.csv")[c("Code","X2018.people.per.sq..km")]
earnings =read.csv("data/ann_earning_2018_perLA.csv")[c("Code","Mean_ann_earnings")]
age = read.csv("data/processed_median_age_of_population_perLA.csv")[c("Code","median_age_2018","Name")]
covid_deaths = read.csv("data/covid_deaths_until10April_byAreaCode.csv")

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

write.csv(merged_covid_dt_LA, "data_output/merged_covid_cov_dt_LA.csv")
```

re-read here

``` r
covid_air_dt = read.csv("data_output/merged_covidAir_cov_dt_LA.csv", na.strings = "x")
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

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

##### Models

Negative binomial regression model

``` r
pm25_deaths.nb <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                   Mean_ann_earnings + median_age_2018 + pm25_val)
car::vif(pm25_deaths.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.080989                1.441089                2.135625 
    ##                pm25_val 
    ##                1.950479

``` r
pm10_deaths.nb <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                   Mean_ann_earnings + median_age_2018 + pm10_val)
car::vif(pm10_deaths.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.099725                1.385040                1.950317 
    ##                pm10_val 
    ##                1.667899

``` r
nox_deaths.nb <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                  Mean_ann_earnings + median_age_2018 + nox_val)
car::vif(nox_deaths.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.529754                1.277034                2.128073 
    ##                 nox_val 
    ##                2.650417

``` r
no2_deaths.nb <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                  Mean_ann_earnings + median_age_2018 + no2_val)
car::vif(no2_deaths.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.417227                1.263175                2.314129 
    ##                 no2_val 
    ##                2.736928

``` r
o3_deaths.nb <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                 Mean_ann_earnings + median_age_2018 + o3_val)
car::vif(o3_deaths.nb)
```

    ## X2018.people.per.sq..km       Mean_ann_earnings         median_age_2018 
    ##                2.094533                1.278499                1.827214 
    ##                  o3_val 
    ##                1.148228

``` r
#summary(so2_deaths.nb <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + #
#                                  Mean_ann_earnings + median_age_2018 + so2_val))
#car::vif(so2_deaths.nb)
```

We note that annual earnings in our models is not significant. We
therefore proceed to remove this variable.

``` r
# glm -earnings, since not significant 
summary(pm25_deaths.nb_red <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                    median_age_2018 + pm25_val))
car::vif(pm25_deaths.nb_red)
summary(pm10_deaths.nb_red <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                    median_age_2018 + pm10_val))
car::vif(pm10_deaths.nb_red)
summary(nox_deaths.nb_red <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                   median_age_2018 + nox_val))
car::vif(nox_deaths.nb_red)
summary(no2_deaths.nb_red <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                   median_age_2018 + no2_val))
car::vif(no2_deaths.nb_red)
summary(o3_deaths.nb_red <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                  median_age_2018 + o3_val))
car::vif(o3_deaths.nb_red)
summary(so2_deaths.nb_red <- glm.nb(data = covid_air_dt, deaths ~ X2018.people.per.sq..km + 
                                   median_age_2018 + so2_val))
car::vif(so2_deaths.nb_red)
```

<br>

##### Calculate odds ratios

``` r
pm25_deaths.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(pm25_deaths.nb), confint(pm25_deaths.nb))), p_value = summary(pm25_deaths.nb)$coefficients[,4]))
pm10_deaths.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(pm10_deaths.nb), confint(pm10_deaths.nb))), p_value = summary(pm10_deaths.nb)$coefficients[,4]))
nox_deaths.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(nox_deaths.nb), confint(nox_deaths.nb))), p_value = summary(nox_deaths.nb)$coefficients[,4]))
no2_deaths.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(no2_deaths.nb), confint(no2_deaths.nb))), p_value = summary(no2_deaths.nb)$coefficients[,4]))
o3_deaths.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(o3_deaths.nb), confint(o3_deaths.nb))), p_value = summary(o3_deaths.nb)$coefficients[,4]))
#so2_deaths.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(so2_deaths.nb), confint(so2_deaths.nb))), p_value = summary(so2_deaths.nb)$coefficients[,4]))
#make a df just with the pollutants 

LA_covid_onlyPoll = data.frame(rbind(pm25_deaths.nb_mrr[nrow(pm25_deaths.nb_mrr),],
                                      pm10_deaths.nb_mrr[nrow(pm10_deaths.nb_mrr),],
                                      nox_deaths.nb_mrr[nrow(nox_deaths.nb_mrr),],
                                      no2_deaths.nb_mrr[nrow(no2_deaths.nb_mrr),],
                                      o3_deaths.nb_mrr[nrow(o3_deaths.nb_mrr),]))
```

##### Plot odds ratios

sort data

``` r
LA_covid_onlyPoll$names = row.names(LA_covid_onlyPoll)
LA_covid_onlyPoll$significance = "p-value > 0.05"
LA_covid_onlyPoll$significance[LA_covid_onlyPoll$p_value < 0.05] <- "p-value < 0.05"
write.csv(LA_covid_onlyPoll, "data_output/LA_covid_onlyPoll_preL_AP2018.csv", row.names = F)
```

plot

``` r
ggplot(LA_covid_onlyPoll, aes(x=reorder(names, OR), y=OR, color=significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_classic() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("Mortality rate ratios") + 
  xlab("Pollutants")+
  ylim(0.75, 1.6)
```

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
ggsave('fig_out/LA_deaths_MRR.pdf')
```

### Cases in local authorities

#### Data curation

``` r
cases_LA_raw = read.csv("data/coronavirus-cases_latest-18_5_2020.csv")
cases_LA_dt = subset(cases_LA_raw, Specimen.date == "2020-04-10", select = c(Area.code, Cumulative.lab.confirmed.cases))
cases_LA_dt.agg <-aggregate(cases_LA_dt, by=list(cases_LA_dt$Area.code), 
  FUN=mean, na.rm=TRUE)
nrow(cases_LA_dt.agg)
```

    ## [1] 342

Merge cases data to main data

``` r
deaths_LA_dt = read.csv("data_output/merged_covidAir_cov_dt_LA.csv", na.strings = "x")
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

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

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
          o3_cases.nb, so2_cases.nb,type="html",out = "fig_out/LA_covid_allPoll_CASES_nb.html",
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
pm25_cases.nb_irr = data.frame(cbind(exp(cbind(OR = coef(pm25_cases.nb), confint(pm25_cases.nb))), p_value = summary(pm25_cases.nb)$coefficients[,4]))
pm10_cases.nb_irr = data.frame(cbind(exp(cbind(OR = coef(pm10_cases.nb), confint(pm10_cases.nb))), p_value = summary(pm10_cases.nb)$coefficients[,4]))
nox_cases.nb_irr = data.frame(cbind(exp(cbind(OR = coef(nox_cases.nb), confint(nox_cases.nb))), p_value = summary(nox_cases.nb)$coefficients[,4]))
no2_cases.nb_irr = data.frame(cbind(exp(cbind(OR = coef(no2_cases.nb), confint(no2_cases.nb))), p_value = summary(no2_cases.nb)$coefficients[,4]))
o3_cases.nb_irr = data.frame(cbind(exp(cbind(OR = coef(o3_cases.nb), confint(o3_cases.nb))), p_value = summary(o3_cases.nb)$coefficients[,4]))
so2_cases.nb_irr = data.frame(cbind(exp(cbind(OR = coef(so2_cases.nb), confint(so2_cases.nb))), p_value = summary(so2_cases.nb)$coefficients[,4]))
#make a df just with the pollutants 

LA_covid_CASES_onlyPoll = data.frame(rbind(pm25_cases.nb_irr[nrow(pm25_cases.nb_irr),],
                                      pm10_cases.nb_irr[nrow(pm10_cases.nb_irr),],
                                      nox_cases.nb_irr[nrow(nox_cases.nb_irr),],
                                      no2_cases.nb_irr[nrow(no2_cases.nb_irr),],
                                      o3_cases.nb_irr[nrow(o3_cases.nb_irr),],
                                      so2_cases.nb_irr[nrow(so2_cases.nb_irr),]))
```

##### Plot odds ratios

sort data

``` r
LA_covid_CASES_onlyPoll$names = row.names(LA_covid_CASES_onlyPoll)
LA_covid_CASES_onlyPoll$significance = "p-value > 0.05"
LA_covid_CASES_onlyPoll$significance[LA_covid_CASES_onlyPoll$p_value < 0.05] <- "p-value < 0.05"

write.csv(LA_covid_CASES_onlyPoll, "data_output/LA_covid_CASES_onlyPoll_preL_2018AP.csv", row.names = F)
```

plot

``` r
ggplot(LA_covid_CASES_onlyPoll, aes(x=reorder(names, OR), y=OR, color=significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_classic() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("Infectivity rate ratios") + 
  xlab("Pollutants")+
  ylim(0.75, 1.6)
```

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
ggsave('fig_out/LA_CASES_odds_IRR.pdf')
```

Save the numbers

``` r
library(stargazer)
stargazer(LA_covid_CASES_onlyPoll, summary=FALSE ,type="html",out 
          ="fig_out/LA_covid_infectivity_IRR.html")
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
          o3_cases.nb, so2_cases.nb,type="html",out = "fig_out/LA_covid_allPoll_CASES_nb.html",
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
          o3_deaths.nb,type="html",out = "fig_out/LA_covid_allPoll_nb.html",
          dep.var.labels="Number of COVID-19-related deaths",
          single.row=TRUE)
```

<table style="text-align:center">

<tr>

<td colspan="6" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="5">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="5" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="5">

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

</tr>

<tr>

<td colspan="6" style="border-bottom: 1px solid black">

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

</tr>

<tr>

<td colspan="6" style="border-bottom: 1px solid black">

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

</tr>

<tr>

<td colspan="6" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="5" style="text-align:right">

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

### Analysis using the 5-year average of air pollution values

#### Average PCM data for all pollutants

NO2

``` r
no2_2014 = read.csv("data/raw_mapno22014.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'no22014')]
no2_2015 = read.csv("data/raw_mapno22015.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'no22015')]
no2_2016 = read.csv("data/raw_mapno22016.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'no22016')]
no2_2017 = read.csv("data/raw_mapno22017.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'no22017')]
no2_2018 = read.csv("data/raw_mapno22018.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'no22018')]
no2_5y = Reduce(merge, list(no2_2014,no2_2015,no2_2016,no2_2017,no2_2018))
no2_5y$no2_5yAvg = (no2_5y$no22014+no2_5y$no22015+no2_5y$no22016+no2_5y$no22017+no2_5y$no22018)/5
no2_5y = no2_5y[ ,c('ukgridcode', 'no22018','no2_5yAvg')]
```

NOX

``` r
nox_2014 = read.csv("data/raw_mapnox2014.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'nox2014')]
nox_2015 = read.csv("data/raw_mapnox2015.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'nox2015')]
nox_2016 = read.csv("data/raw_mapnox2016.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'nox2016')]
nox_2017 = read.csv("data/raw_mapnox2017.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'nox2017')]
nox_2018 = read.csv("data/raw_mapnox2018.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'nox2018')]
nox_5y = Reduce(merge, list(nox_2014,nox_2015,nox_2016,nox_2017,nox_2018))
nox_5y$nox_5yAvg = (nox_5y$nox2014+nox_5y$nox2015+nox_5y$nox2016+nox_5y$nox2017+nox_5y$nox2018)/5
nox_5y = nox_5y[ ,c('ukgridcode', 'nox2018','nox_5yAvg')]
```

PM25

``` r
pm25_2014 = read.csv("data/raw_mappm252014.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'pm252014g')]
pm25_2015 = read.csv("data/raw_mappm252015.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'pm252015g')]
pm25_2016 = read.csv("data/raw_mappm252016.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'pm252016g')]
pm25_2017 = read.csv("data/raw_mappm252017.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'pm252017g')]
pm25_2018 = read.csv("data/raw_mappm252018.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'pm252018g')]
pm25_5y = Reduce(merge, list(pm25_2014,pm25_2015,pm25_2016,pm25_2017,pm25_2018))
pm25_5y$pm25_5yAvg = (pm25_5y$pm252014g+pm25_5y$pm252015g+pm25_5y$pm252016g+pm25_5y$pm252017g+pm25_5y$pm252018g)/5
pm25_5y = pm25_5y[ ,c('ukgridcode', 'pm252018g','pm25_5yAvg')]
```

PM10

``` r
pm10_2014 = read.csv("data/raw_mappm102014.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'pm102014g')]
pm10_2015 = read.csv("data/raw_mappm102015.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'pm102015g')]
pm10_2016 = read.csv("data/raw_mappm102016.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'pm102016g')]
pm10_2017 = read.csv("data/raw_mappm102017.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'pm102017g')]
pm10_2018 = read.csv("data/raw_mappm102018.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'pm102018g')]
pm10_5y = Reduce(merge, list(pm10_2014,pm10_2015,pm10_2016,pm10_2017,pm10_2018))
pm10_5y$pm10_5yAvg = (pm10_5y$pm102014g+pm10_5y$pm102015g+pm10_5y$pm102016g+pm10_5y$pm102017g+pm10_5y$pm102018g)/5
pm10_5y = pm10_5y[ ,c('ukgridcode', 'pm102018g','pm10_5yAvg')]
```

O3

``` r
o3_2014 = read.csv("data/raw_mapo312014.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'dgt12014')]
o3_2015 = read.csv("data/raw_mapo312015.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'dgt12015')]
o3_2016 = read.csv("data/raw_mapo312016.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'dgt12016')]
o3_2017 = read.csv("data/raw_mapo32017.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'dgt12017')]
o3_2018 = read.csv("data/raw_mapo312018.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'dgt12018')]
o3_5y = Reduce(merge, list(o3_2014,o3_2015,o3_2016,o3_2017,o3_2018))
o3_5y$o3_5yAvg = (o3_5y$dgt12014+o3_5y$dgt12015+o3_5y$dgt12016+o3_5y$dgt12017+o3_5y$dgt12018)/5
o3_5y = o3_5y[ ,c('ukgridcode', 'dgt12018','o3_5yAvg')]
```

##### Create merged df with 2018 and 5yAvg

``` r
ap_dt = Reduce(merge, list(no2_5y,nox_5y,pm25_5y,pm10_5y,o3_5y))
write.csv(ap_dt, "data_output/5Yaverage_AP_PCMdata.csv", row.names = F)
write.csv(na.omit(ap_dt), "data_output/5Yaverage_AP_PCMdata_na.csv", row.names = F)
```

Convert the XY coordinates to lon-lat

``` r
ap_dt = read.csv("data_output/5Yaverage_AP_PCMdata_na.csv")

# Add X Y coordinates from any of the data 

o3_2014 = read.csv("data/raw_mapo312014.csv", skip = 5, header = T,na.strings="MISSING")[ ,c('ukgridcode', 'x','y')]
ap_dt = merge(ap_dt,o3_2014, by = "ukgridcode")
ap_dt_out = ap_dt[ , -which(names(ap_dt) %in% c("x","y"))]
### convert to lon lat ----

# transform: easting/northing -> long/lat
### shortcuts (from https://stephendavidgregory.github.io/useful/UKgrid_to_LatLon)
ukgrid <- "+init=epsg:27700"
latlong <- "+init=epsg:4326"

### Create coordinates variable
ap_coords <- cbind(Easting = as.numeric(as.character(ap_dt$x)),
                Northing = as.numeric(as.character(ap_dt$y)))

### Create the SpatialPointsDataFrame
ap_SP <- SpatialPointsDataFrame(ap_coords,
                                 data = ap_dt,
                                 proj4string = CRS("+init=epsg:27700"))

### Convert
ap_LL <- spTransform(ap_SP, CRS(latlong))

ap_ll_df = data.frame('lon' = coordinates(ap_LL)[, 1], 'lat' = coordinates(ap_LL)[, 2], 'ukgridcode' = ap_LL$ukgridcode)

ap_dt_ll = merge(ap_dt_out,ap_ll_df, by = "ukgridcode")

write.csv(ap_dt_ll,"data_output/processed_allAP_lonlat.csv", row.names = F)
```

#### Model of Mortality based on 5YA, pre-lockdown

\#\#\#\#\# Data curation

``` r
avgAP_dt = read.csv("data_output/merged_covidAir_cov_LA_5YA.csv", na.strings = "x")

avgAP_dt$X2018.people.per.sq..km = as.numeric(gsub(",","",avgAP_dt$X2018.people.per.sq..km))
avgAP_dt$Mean_ann_earnings = as.numeric(gsub(",","",avgAP_dt$Mean_ann_earnings))
avgAP_dt$X2018.people.per.sq..km = as.numeric(avgAP_dt$X2018.people.per.sq..km)
avgAP_dt$Mean_ann_earnings = as.numeric(avgAP_dt$Mean_ann_earnings)
```

Additional variables for 5YA

``` r
avgEarnings_2014 = na.omit(read.csv("data_output/avgEarnings_2014.csv", na.strings = "x")[,c("Code","avgEarn_2014")])
avgEarnings_2015 = na.omit(read.csv("data_output/avgEarnings_2015.csv", na.strings = "x")[,c("Code","avgEarn_2015")])
avgEarnings_2016 = na.omit(read.csv("data_output/avgEarnings_2016.csv", na.strings = "x")[,c("Code","avgEarn_2016")])
avgEarnings_2017 = na.omit(read.csv("data_output/avgEarnings_2017.csv", na.strings = "x")[,c("Code","avgEarn_2017")])
avgEarnings_2018 = na.omit(read.csv("data_output/avgEarnings_2018.csv", na.strings = "x")[,c("Code","avgEarn_2018")])

avgEarnings_5ya = Reduce(merge, list(avgEarnings_2014,avgEarnings_2015,avgEarnings_2016,avgEarnings_2017,avgEarnings_2018))

avgEarnings_5ya$avgEarn_2014 = as.numeric(gsub(",","",avgEarnings_5ya$avgEarn_2014))
avgEarnings_5ya$avgEarn_2015 = as.numeric(gsub(",","",avgEarnings_5ya$avgEarn_2015))
avgEarnings_5ya$avgEarn_2016 = as.numeric(gsub(",","",avgEarnings_5ya$avgEarn_2016))
avgEarnings_5ya$avgEarn_2017 = as.numeric(gsub(",","",avgEarnings_5ya$avgEarn_2017))
avgEarnings_5ya$avgEarn_2018 = as.numeric(gsub(",","",avgEarnings_5ya$avgEarn_2018))

avgEarnings_5ya$earnings_5ya = (avgEarnings_5ya$avgEarn_2014 + avgEarnings_5ya$avgEarn_2015 + avgEarnings_5ya$avgEarn_2016 + avgEarnings_5ya$avgEarn_2017 + avgEarnings_5ya$avgEarn_2018)/5

avgEarnings_5ya = avgEarnings_5ya[,c("Code","earnings_5ya")]

Pop_density_5ya = read.csv("data/Pop_density_full_dt.csv")[,c("Code","pop_dens_5ya")]

age_5ya = read.csv("data_output/age_5ya.csv")[,c("Code","age_5ya")]

covariates_5ya = Reduce(merge, list(avgEarnings_5ya,Pop_density_5ya,age_5ya))
#avgAP_dt = merge(avgAP_dt, covariates_5ya, by = "Code")
avgAP_dt = merge(avgAP_dt, covariates_5ya[,-3], by = "Code")
avgAP_dt$pop_dens_5ya = as.numeric(gsub(",","",avgAP_dt$pop_dens_5ya))
write.csv(avgAP_dt,"data_output/merged_covidAir_cov_LA_5YA.csv", row.names = F)
```

``` r
summary(pm25_deaths_5YA.nb <- glm.nb(data = avgAP_dt, deaths ~ pop_dens_5ya +
                                   earnings_5ya + age_5ya + pm25_5yAvg))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ pop_dens_5ya + earnings_5ya + age_5ya + 
    ##     pm25_5yAvg, data = avgAP_dt, init.theta = 2.228067951, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.5294  -0.8876  -0.3636   0.2695   3.3767  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   6.390e+00  5.943e-01  10.752  < 2e-16 ***
    ## pop_dens_5ya  9.994e-05  2.273e-05   4.397 1.10e-05 ***
    ## earnings_5ya -2.086e-07  6.848e-06  -0.030    0.976    
    ## age_5ya      -7.057e-02  1.204e-02  -5.861 4.59e-09 ***
    ## pm25_5yAvg    7.720e-03  1.938e-02   0.398    0.690    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.2281) family taken to be 1)
    ## 
    ##     Null deviance: 487.47  on 282  degrees of freedom
    ## Residual deviance: 300.11  on 278  degrees of freedom
    ## AIC: 2612.1
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.228 
    ##           Std. Err.:  0.188 
    ## 
    ##  2 x log-likelihood:  -2600.111

``` r
car::vif(pm25_deaths_5YA.nb)
```

    ## pop_dens_5ya earnings_5ya      age_5ya   pm25_5yAvg 
    ##     2.142889     1.189121     1.982406     1.134312

``` r
summary(pm10_deaths_5YA.nb <- glm.nb(data = avgAP_dt, deaths ~ pop_dens_5ya +
                                   earnings_5ya + age_5ya + pm10_5yAvg))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ pop_dens_5ya + earnings_5ya + age_5ya + 
    ##     pm10_5yAvg, data = avgAP_dt, init.theta = 2.227020703, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.5260  -0.8855  -0.3567   0.2658   3.3505  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   6.456e+00  5.880e-01  10.980  < 2e-16 ***
    ## pop_dens_5ya  9.992e-05  2.276e-05   4.390 1.13e-05 ***
    ## earnings_5ya  2.476e-08  6.845e-06   0.004    0.997    
    ## age_5ya      -7.137e-02  1.193e-02  -5.983 2.20e-09 ***
    ## pm10_5yAvg    2.242e-03  1.366e-02   0.164    0.870    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.227) family taken to be 1)
    ## 
    ##     Null deviance: 487.25  on 282  degrees of freedom
    ## Residual deviance: 300.11  on 278  degrees of freedom
    ## AIC: 2612.2
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.227 
    ##           Std. Err.:  0.187 
    ## 
    ##  2 x log-likelihood:  -2600.239

``` r
car::vif(pm10_deaths_5YA.nb)
```

    ## pop_dens_5ya earnings_5ya      age_5ya   pm10_5yAvg 
    ##     2.148181     1.187764     1.945466     1.101626

``` r
summary(nox_deaths_5YA.nb <- glm.nb(data = avgAP_dt, deaths ~ pop_dens_5ya +
                                   earnings_5ya + age_5ya + nox_5yAvg))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ pop_dens_5ya + earnings_5ya + age_5ya + 
    ##     nox_5yAvg, data = avgAP_dt, init.theta = 2.330563559, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.4959  -0.8776  -0.3217   0.2737   3.8186  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   5.264e+00  5.935e-01   8.870  < 2e-16 ***
    ## pop_dens_5ya  8.440e-05  2.331e-05   3.621 0.000293 ***
    ## earnings_5ya  2.352e-06  6.711e-06   0.350 0.726007    
    ## age_5ya      -5.045e-02  1.244e-02  -4.055    5e-05 ***
    ## nox_5yAvg     1.469e-02  4.179e-03   3.515 0.000440 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.3306) family taken to be 1)
    ## 
    ##     Null deviance: 508.44  on 282  degrees of freedom
    ## Residual deviance: 299.44  on 278  degrees of freedom
    ## AIC: 2598.9
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.331 
    ##           Std. Err.:  0.198 
    ## 
    ##  2 x log-likelihood:  -2586.947

``` r
car::vif(nox_deaths_5YA.nb)
```

    ## pop_dens_5ya earnings_5ya      age_5ya    nox_5yAvg 
    ##     2.352705     1.191368     2.207728     1.827694

``` r
summary(no2_deaths_5YA.nb <- glm.nb(data = avgAP_dt, deaths ~ pop_dens_5ya +
                                   earnings_5ya + age_5ya + no2_5yAvg))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ pop_dens_5ya + earnings_5ya + age_5ya + 
    ##     no2_5yAvg, data = avgAP_dt, init.theta = 2.327684382, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.5071  -0.8829  -0.3292   0.2712   3.8454  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   5.232e+00  6.097e-01   8.582  < 2e-16 ***
    ## pop_dens_5ya  8.883e-05  2.288e-05   3.882 0.000103 ***
    ## earnings_5ya  2.010e-06  6.712e-06   0.299 0.764600    
    ## age_5ya      -5.038e-02  1.261e-02  -3.995 6.47e-05 ***
    ## no2_5yAvg     2.384e-02  7.044e-03   3.385 0.000713 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.3277) family taken to be 1)
    ## 
    ##     Null deviance: 507.85  on 282  degrees of freedom
    ## Residual deviance: 299.49  on 278  degrees of freedom
    ## AIC: 2599.3
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.328 
    ##           Std. Err.:  0.197 
    ## 
    ##  2 x log-likelihood:  -2587.340

``` r
car::vif(no2_deaths_5YA.nb)
```

    ## pop_dens_5ya earnings_5ya      age_5ya    no2_5yAvg 
    ##     2.263950     1.190001     2.265576     1.747459

``` r
summary(o3_deaths_5YA.nb <- glm.nb(data = avgAP_dt, deaths ~ pop_dens_5ya +
                                   earnings_5ya + age_5ya + o3_5yAvg))
```

    ## 
    ## Call:
    ## glm.nb(formula = deaths ~ pop_dens_5ya + earnings_5ya + age_5ya + 
    ##     o3_5yAvg, data = avgAP_dt, init.theta = 2.353153533, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.3590  -0.8889  -0.3458   0.3287   3.3772  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   6.250e+00  5.195e-01  12.031  < 2e-16 ***
    ## pop_dens_5ya  7.587e-05  2.285e-05   3.320 0.000899 ***
    ## earnings_5ya  9.901e-06  6.856e-06   1.444 0.148705    
    ## age_5ya      -5.912e-02  1.171e-02  -5.047 4.49e-07 ***
    ## o3_5yAvg     -1.793e-01  4.201e-02  -4.268 1.97e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.3532) family taken to be 1)
    ## 
    ##     Null deviance: 513.04  on 282  degrees of freedom
    ## Residual deviance: 299.22  on 278  degrees of freedom
    ## AIC: 2596.1
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.353 
    ##           Std. Err.:  0.200 
    ## 
    ##  2 x log-likelihood:  -2584.060

``` r
car::vif(o3_deaths_5YA.nb)
```

    ## pop_dens_5ya earnings_5ya      age_5ya     o3_5yAvg 
    ##     2.282899     1.257509     1.969455     1.311695

``` r
stargazer(pm25_deaths_5YA.nb, pm10_deaths_5YA.nb, nox_deaths_5YA.nb, no2_deaths_5YA.nb,o3_deaths_5YA.nb,type="html",out = "fig_out/LA_covid_allPoll_DEATHS_5YA_nb.html",
          dep.var.labels="Number of COVID-19-related deaths, 5-year average",
          single.row=TRUE)
```

<table style="text-align:center">

<tr>

<td colspan="6" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="5">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="5" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="5">

Number of COVID-19-related deaths, 5-year average

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

</tr>

<tr>

<td colspan="6" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

pop\_dens\_5ya

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

earnings\_5ya

</td>

<td>

\-0.00000 (0.00001)

</td>

<td>

0.00000 (0.00001)

</td>

<td>

0.00000 (0.00001)

</td>

<td>

0.00000 (0.00001)

</td>

<td>

0.00001 (0.00001)

</td>

</tr>

<tr>

<td style="text-align:left">

age\_5ya

</td>

<td>

\-0.071<sup>\*\*\*</sup> (0.012)

</td>

<td>

\-0.071<sup>\*\*\*</sup> (0.012)

</td>

<td>

\-0.050<sup>\*\*\*</sup> (0.012)

</td>

<td>

\-0.050<sup>\*\*\*</sup> (0.013)

</td>

<td>

\-0.059<sup>\*\*\*</sup> (0.012)

</td>

</tr>

<tr>

<td style="text-align:left">

pm25\_5yAvg

</td>

<td>

0.008 (0.019)

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

pm10\_5yAvg

</td>

<td>

</td>

<td>

0.002 (0.014)

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

nox\_5yAvg

</td>

<td>

</td>

<td>

</td>

<td>

0.015<sup>\*\*\*</sup> (0.004)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

no2\_5yAvg

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

0.024<sup>\*\*\*</sup> (0.007)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

o3\_5yAvg

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

\-0.179<sup>\*\*\*</sup> (0.042)

</td>

</tr>

<tr>

<td style="text-align:left">

Constant

</td>

<td>

6.390<sup>\*\*\*</sup> (0.594)

</td>

<td>

6.456<sup>\*\*\*</sup> (0.588)

</td>

<td>

5.264<sup>\*\*\*</sup> (0.594)

</td>

<td>

5.232<sup>\*\*\*</sup> (0.610)

</td>

<td>

6.250<sup>\*\*\*</sup> (0.519)

</td>

</tr>

<tr>

<td colspan="6" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Observations

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

</tr>

<tr>

<td style="text-align:left">

Log Likelihood

</td>

<td>

\-1,301.056

</td>

<td>

\-1,301.119

</td>

<td>

\-1,294.474

</td>

<td>

\-1,294.670

</td>

<td>

\-1,293.030

</td>

</tr>

<tr>

<td style="text-align:left">

theta

</td>

<td>

2.228<sup>\*\*\*</sup> (0.188)

</td>

<td>

2.227<sup>\*\*\*</sup> (0.187)

</td>

<td>

2.331<sup>\*\*\*</sup> (0.198)

</td>

<td>

2.328<sup>\*\*\*</sup> (0.197)

</td>

<td>

2.353<sup>\*\*\*</sup> (0.200)

</td>

</tr>

<tr>

<td style="text-align:left">

Akaike Inf. Crit.

</td>

<td>

2,612.111

</td>

<td>

2,612.239

</td>

<td>

2,598.947

</td>

<td>

2,599.340

</td>

<td>

2,596.060

</td>

</tr>

<tr>

<td colspan="6" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="5" style="text-align:right">

<sup>*</sup>p\<0.1; <sup>**</sup>p\<0.05; <sup>***</sup>p\<0.01

</td>

</tr>

</table>

##### Plot odds ratios

``` r
pm25_deaths_5YA.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(pm25_deaths_5YA.nb), confint(pm25_deaths_5YA.nb))), p_value = summary(pm25_deaths_5YA.nb)$coefficients[,4]))
pm10_deaths_5YA.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(pm10_deaths_5YA.nb), confint(pm10_deaths_5YA.nb))), p_value = summary(pm10_deaths_5YA.nb)$coefficients[,4]))
nox_deaths_5YA.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(nox_deaths_5YA.nb), confint(nox_deaths_5YA.nb))), p_value = summary(nox_deaths_5YA.nb)$coefficients[,4]))
no2_deaths_5YA.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(no2_deaths_5YA.nb), confint(no2_deaths_5YA.nb))), p_value = summary(no2_deaths_5YA.nb)$coefficients[,4]))
o3_deaths_5YA.nb_mrr = data.frame(cbind(exp(cbind(OR = coef(o3_deaths_5YA.nb), confint(o3_deaths_5YA.nb))), p_value = summary(o3_deaths_5YA.nb)$coefficients[,4]))
#make a df just with the pollutants 

LA_covid_DEATHS_onlyPoll_5YA = data.frame(rbind(pm25_deaths_5YA.nb_mrr[nrow(pm25_deaths_5YA.nb_mrr),],
                pm10_deaths_5YA.nb_mrr[nrow(pm10_deaths_5YA.nb_mrr),],
                nox_deaths_5YA.nb_mrr[nrow(nox_deaths_5YA.nb_mrr),],
                no2_deaths_5YA.nb_mrr[nrow(no2_deaths_5YA.nb_mrr),],
                o3_deaths_5YA.nb_mrr[nrow(o3_deaths_5YA.nb_mrr),]))
```

sort data

``` r
LA_covid_DEATHS_onlyPoll_5YA$names = row.names(LA_covid_DEATHS_onlyPoll_5YA)
LA_covid_DEATHS_onlyPoll_5YA$significance = "p-value > 0.05"
LA_covid_DEATHS_onlyPoll_5YA$significance[LA_covid_DEATHS_onlyPoll_5YA$p_value < 0.05] <- "p-value < 0.05"
write.csv(LA_covid_DEATHS_onlyPoll_5YA,"data_output/LA_covid_DEATHS_onlyPoll_5YA.csv", row.names = F)
```

plot

``` r
ggplot(LA_covid_DEATHS_onlyPoll_5YA, aes(x=reorder(names, OR), y=OR, color=significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("Mortality rate ratios, pre-lockdown") + 
  xlab("Pollutants, 5 year average (2014-2018)") + theme_classic()
```

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

``` r
ggsave('fig_out/LA_CASES_odds_MRR_5YA.pdf')
```

#### Infectivity analysis with 5YA air pollution data, pre lockdown

``` r
cases_LA_raw = read.csv("data/coronavirus-cases_latest-18_5_2020.csv")
cases_LA_dt = subset(cases_LA_raw, Specimen.date == "2020-04-10", select = c(Area.code, Cumulative.lab.confirmed.cases))
library(dplyr)
cases_LA_dt = cases_LA_dt %>% 
group_by(Area.code) %>% 
summarise(cases = mean(Cumulative.lab.confirmed.cases))
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
nrow(cases_LA_dt)
```

    ## [1] 342

##### Model

``` r
avgAP_dt = read.csv("data_output/merged_covidAir_cov_LA_5YA.csv", na.strings = "x")
avgAP_dt$pop_dens_5ya = as.numeric(gsub(",","",avgAP_dt$pop_dens_5ya))
avgAP_dt$X2018.people.per.sq..km = as.numeric(gsub(",","",avgAP_dt$X2018.people.per.sq..km))
avgAP_dt$Mean_ann_earnings = as.numeric(gsub(",","",avgAP_dt$Mean_ann_earnings))
avgAP_dt$X2018.people.per.sq..km = as.numeric(avgAP_dt$X2018.people.per.sq..km)
avgAP_dt$Mean_ann_earnings = as.numeric(avgAP_dt$Mean_ann_earnings)
#avgAP_dt = merge(avgAP_dt,cases_LA_dt, by.x = "Code",by.y = "Area.code")
#write.csv(avgAP_dt,"data_output/merged_covidAir_cov_LA_5YA.csv", row.names = F)
```

``` r
library(MASS)
pm25_cases_5YA.nb <- glm.nb(data = avgAP_dt, cases ~ pop_dens_5ya +
                                   earnings_5ya + age_5ya + pm25_5yAvg)
car::vif(pm25_cases_5YA.nb)
```

    ## pop_dens_5ya earnings_5ya      age_5ya   pm25_5yAvg 
    ##     2.127469     1.183202     1.976028     1.133826

``` r
pm10_cases_5YA.nb <- glm.nb(data = avgAP_dt, cases ~ pop_dens_5ya +
                                   earnings_5ya + age_5ya + pm10_5yAvg)
car::vif(pm10_cases_5YA.nb)
```

    ## pop_dens_5ya earnings_5ya      age_5ya   pm10_5yAvg 
    ##     2.132595     1.181741     1.939209     1.100561

``` r
nox_cases_5YA.nb <- glm.nb(data = avgAP_dt, cases ~ pop_dens_5ya +
                                   earnings_5ya + age_5ya + nox_5yAvg)
car::vif(nox_cases_5YA.nb)
```

    ## pop_dens_5ya earnings_5ya      age_5ya    nox_5yAvg 
    ##     2.334119     1.184576     2.204639     1.828731

``` r
no2_cases_5YA.nb <- glm.nb(data = avgAP_dt, cases ~ pop_dens_5ya +
                                   earnings_5ya + age_5ya + no2_5yAvg)
car::vif(no2_cases_5YA.nb)
```

    ## pop_dens_5ya earnings_5ya      age_5ya    no2_5yAvg 
    ##     2.246629     1.183251     2.261649     1.748043

``` r
o3_cases_5YA.nb <- glm.nb(data = avgAP_dt, cases ~ pop_dens_5ya +
                                   earnings_5ya + age_5ya + o3_5yAvg)
car::vif(o3_cases_5YA.nb)
```

    ## pop_dens_5ya earnings_5ya      age_5ya     o3_5yAvg 
    ##     2.261263     1.250229     1.966950     1.314072

``` r
stargazer(pm25_cases_5YA.nb, pm10_cases_5YA.nb, nox_cases_5YA.nb, no2_cases_5YA.nb,o3_cases_5YA.nb,type="html",out = "fig_out/LA_covid_allPoll_CASES_5YA_nb.html",
          dep.var.labels="Number of COVID-19-related cases, pre-lockdown, 5-year averages",
          single.row=TRUE)
```

<table style="text-align:center">

<tr>

<td colspan="6" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="5">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="5" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="5">

Number of COVID-19-related cases, pre-lockdown, 5-year averages

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

</tr>

<tr>

<td colspan="6" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

pop\_dens\_5ya

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

0.00005<sup>\*\*</sup> (0.00002)

</td>

</tr>

<tr>

<td style="text-align:left">

earnings\_5ya

</td>

<td>

\-0.00000 (0.00001)

</td>

<td>

\-0.00000 (0.00001)

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

</tr>

<tr>

<td style="text-align:left">

age\_5ya

</td>

<td>

\-0.082<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.082<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.062<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.062<sup>\*\*\*</sup> (0.012)

</td>

<td>

\-0.063<sup>\*\*\*</sup> (0.010)

</td>

</tr>

<tr>

<td style="text-align:left">

pm25\_5yAvg

</td>

<td>

\-0.028 (0.018)

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

pm10\_5yAvg

</td>

<td>

</td>

<td>

\-0.023<sup>\*</sup> (0.012)

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

nox\_5yAvg

</td>

<td>

</td>

<td>

</td>

<td>

0.012<sup>\*\*\*</sup> (0.004)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

no2\_5yAvg

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

0.019<sup>\*\*\*</sup> (0.006)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

o3\_5yAvg

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

\-0.228<sup>\*\*\*</sup> (0.037)

</td>

</tr>

<tr>

<td style="text-align:left">

Constant

</td>

<td>

8.855<sup>\*\*\*</sup> (0.540)

</td>

<td>

8.900<sup>\*\*\*</sup> (0.533)

</td>

<td>

7.538<sup>\*\*\*</sup> (0.544)

</td>

<td>

7.527<sup>\*\*\*</sup> (0.559)

</td>

<td>

8.194<sup>\*\*\*</sup> (0.460)

</td>

</tr>

<tr>

<td colspan="6" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Observations

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

</tr>

<tr>

<td style="text-align:left">

Log Likelihood

</td>

<td>

\-1,729.358

</td>

<td>

\-1,728.939

</td>

<td>

\-1,725.318

</td>

<td>

\-1,725.733

</td>

<td>

\-1,713.342

</td>

</tr>

<tr>

<td style="text-align:left">

theta

</td>

<td>

2.598<sup>\*\*\*</sup> (0.209)

</td>

<td>

2.605<sup>\*\*\*</sup> (0.210)

</td>

<td>

2.667<sup>\*\*\*</sup> (0.215)

</td>

<td>

2.660<sup>\*\*\*</sup> (0.215)

</td>

<td>

2.882<sup>\*\*\*</sup> (0.234)

</td>

</tr>

<tr>

<td style="text-align:left">

Akaike Inf. Crit.

</td>

<td>

3,468.715

</td>

<td>

3,467.878

</td>

<td>

3,460.637

</td>

<td>

3,461.467

</td>

<td>

3,436.684

</td>

</tr>

<tr>

<td colspan="6" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="5" style="text-align:right">

<sup>*</sup>p\<0.1; <sup>**</sup>p\<0.05; <sup>***</sup>p\<0.01

</td>

</tr>

</table>

##### Plot odds ratios

``` r
pm25_cases_5YA.nb_irr = data.frame(cbind(exp(cbind(OR = coef(pm25_cases_5YA.nb), confint(pm25_cases_5YA.nb))), p_value = summary(pm25_cases_5YA.nb)$coefficients[,4]))
pm10_cases_5YA.nb_irr = data.frame(cbind(exp(cbind(OR = coef(pm10_cases_5YA.nb), confint(pm10_cases_5YA.nb))), p_value = summary(pm10_cases_5YA.nb)$coefficients[,4]))
nox_cases_5YA.nb_irr = data.frame(cbind(exp(cbind(OR = coef(nox_cases_5YA.nb), confint(nox_cases_5YA.nb))), p_value = summary(nox_cases_5YA.nb)$coefficients[,4]))
no2_cases_5YA.nb_irr = data.frame(cbind(exp(cbind(OR = coef(no2_cases_5YA.nb), confint(no2_cases_5YA.nb))), p_value = summary(no2_cases_5YA.nb)$coefficients[,4]))
o3_cases_5YA.nb_irr = data.frame(cbind(exp(cbind(OR = coef(o3_cases_5YA.nb), confint(o3_cases_5YA.nb))), p_value = summary(o3_cases_5YA.nb)$coefficients[,4]))
#make a df just with the pollutants 

LA_covid_cases_onlyPoll_5YA = data.frame(rbind(pm25_cases_5YA.nb_irr[nrow(pm25_cases_5YA.nb_irr),],
                pm10_cases_5YA.nb_irr[nrow(pm10_cases_5YA.nb_irr),],
                nox_cases_5YA.nb_irr[nrow(nox_cases_5YA.nb_irr),],
                no2_cases_5YA.nb_irr[nrow(no2_cases_5YA.nb_irr),],
                o3_cases_5YA.nb_irr[nrow(o3_cases_5YA.nb_irr),]))
```

sort data

``` r
LA_covid_cases_onlyPoll_5YA$names = row.names(LA_covid_cases_onlyPoll_5YA)
LA_covid_cases_onlyPoll_5YA$significance = "p-value > 0.05"
LA_covid_cases_onlyPoll_5YA$significance[LA_covid_cases_onlyPoll_5YA$p_value < 0.05] <- "p-value < 0.05"
write.csv(LA_covid_cases_onlyPoll_5YA,"data_output/LA_covid_cases_onlyPoll_5YA.csv", row.names = F)
```

plot

``` r
ggplot(LA_covid_cases_onlyPoll_5YA, aes(x=reorder(names, OR), y=OR, color=significance)) + 
    geom_point(fill="white", shape=21, size = 2) +
    geom_errorbar(aes(ymin=X2.5.., ymax=X97.5..),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("Infectivity rate ratios, pre-lockdown") + 
  xlab("Pollutants, 5 year average (2014-2018)") + theme_classic()
```

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

``` r
ggsave('fig_out/LA_CASES_odds_IRR_5YA.pdf')
```

## Re-analysis with the cases and deaths numbers until the end of April

### Prepare data: cases and death until the end of April

``` r
april_deaths = read.csv("data_output/LA_deaths_April2020.csv")[,c("Code","april_deaths")]

april_cases = read.csv("data_output/COVID_LA_cases.csv")[,c("lad19_cd","april_26_cases")]
colnames(april_cases) <- c("Code","april_cases")

april_cases= april_cases %>% 
  group_by(Code) %>% 
  summarise(april_cases = sum(april_cases))
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
resub_dt = merge(april_cases,april_deaths, by = "Code")
resub_dt = merge(avgAP_dt,resub_dt, by = "Code", all.x=T)

#ADD 2018 data
resub_dt = merge(resub_dt,
                 read.csv("data_output/merged_covidAir_cov_dt_LA.csv")[,c("Code",
                                                                         "pm25_val",
                                                                         "no2_val",
                                                                         "o3_val",
                                                                         "pm10_val",
                                                                         "nox_val")], by = "Code", all.x=T)

write.csv(resub_dt, "data_output/LA_dt_resub_full.csv", row.names = F)
```

### April cases models

``` r
LA_covid_pm25 <-glm.nb(data = resub_dt, april_cases ~ X2018.people.per.sq..km + Mean_ann_earnings + median_age_2018 +
                          pm25_val)

LA_covid_no2 <-glm.nb(data = resub_dt, april_cases ~ X2018.people.per.sq..km + Mean_ann_earnings + median_age_2018 +
                          no2_val)

LA_covid_nox <-glm.nb(data = resub_dt, april_cases ~ X2018.people.per.sq..km + Mean_ann_earnings + median_age_2018 +
                          nox_val)

LA_covid_pm10 <-glm.nb(data = resub_dt, april_cases ~ X2018.people.per.sq..km + Mean_ann_earnings + median_age_2018 +
                          pm10_val)

LA_covid_o3 <-glm.nb(data = resub_dt, april_cases ~ X2018.people.per.sq..km + Mean_ann_earnings + median_age_2018 +
                          o3_val)


LA_covid_pm25_5yAvg <-glm.nb(data = resub_dt, april_cases ~ pop_dens_5ya + earnings_5ya + age_5ya +
                          pm25_5yAvg)

LA_covid_no2_5yAvg <-glm.nb(data = resub_dt, april_cases ~ pop_dens_5ya + earnings_5ya + age_5ya +
                          no2_5yAvg)

LA_covid_nox_5yAvg <-glm.nb(data = resub_dt, april_cases ~ pop_dens_5ya + earnings_5ya + age_5ya +
                          nox_5yAvg)

LA_covid_pm10_5yAvg <-glm.nb(data = resub_dt, april_cases ~ pop_dens_5ya + earnings_5ya + age_5ya +
                          pm10_5yAvg)

LA_covid_o3_5yAvg <-glm.nb(data = resub_dt, april_cases ~ pop_dens_5ya + earnings_5ya + age_5ya +
                          o3_5yAvg)

stargazer::stargazer(LA_covid_pm25,LA_covid_pm10,LA_covid_no2,LA_covid_nox,LA_covid_o3,
                     LA_covid_pm25_5yAvg,LA_covid_pm10_5yAvg,LA_covid_no2_5yAvg,
                         LA_covid_nox_5yAvg,LA_covid_o3_5yAvg,
                     type="html",out = "fig_out/LA_covid_NB_model_resub_cases.html",
          dep.var.labels="COVID positive or not",
          single.row=TRUE)
```

<table style="text-align:center">

<tr>

<td colspan="11" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="10">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="10" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="10">

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

<td>

(7)

</td>

<td>

(8)

</td>

<td>

(9)

</td>

<td>

(10)

</td>

</tr>

<tr>

<td colspan="11" style="border-bottom: 1px solid black">

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

0.00004 (0.00002)

</td>

<td>

0.00003 (0.00003)

</td>

<td>

0.00004<sup>\*\*</sup> (0.00002)

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

Mean\_ann\_earnings

</td>

<td>

\-0.00000 (0.00001)

</td>

<td>

\-0.00000 (0.00001)

</td>

<td>

\-0.00002<sup>\*\*</sup> (0.00001)

</td>

<td>

\-0.00002<sup>\*\*</sup> (0.00001)

</td>

<td>

0.00000 (0.00001)

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

median\_age\_2018

</td>

<td>

\-0.099<sup>\*\*\*</sup> (0.012)

</td>

<td>

\-0.095<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.055<sup>\*\*\*</sup> (0.013)

</td>

<td>

\-0.059<sup>\*\*\*</sup> (0.012)

</td>

<td>

\-0.072<sup>\*\*\*</sup> (0.011)

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

pm25\_val

</td>

<td>

\-0.120<sup>\*\*\*</sup> (0.031)

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

\-0.087<sup>\*\*\*</sup> (0.020)

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

0.033<sup>\*\*\*</sup> (0.010)

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

</td>

<td>

0.019<sup>\*\*\*</sup> (0.006)

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

\-0.065<sup>\*\*\*</sup> (0.012)

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

pop\_dens\_5ya

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

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.00005<sup>\*\*</sup> (0.00002)

</td>

<td>

0.00005<sup>\*</sup> (0.00002)

</td>

<td>

0.00002 (0.00002)

</td>

</tr>

<tr>

<td style="text-align:left">

earnings\_5ya

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

\-0.00001<sup>\*\*</sup> (0.00001)

</td>

<td>

\-0.00001<sup>\*\*</sup> (0.00001)

</td>

<td>

\-0.00001<sup>\*\*</sup> (0.00001)

</td>

<td>

\-0.00001<sup>\*</sup> (0.00001)

</td>

<td>

\-0.00000 (0.00001)

</td>

</tr>

<tr>

<td style="text-align:left">

age\_5ya

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

\-0.089<sup>\*\*\*</sup> (0.012)

</td>

<td>

\-0.089<sup>\*\*\*</sup> (0.012)

</td>

<td>

\-0.065<sup>\*\*\*</sup> (0.013)

</td>

<td>

\-0.066<sup>\*\*\*</sup> (0.013)

</td>

<td>

\-0.067<sup>\*\*\*</sup> (0.011)

</td>

</tr>

<tr>

<td style="text-align:left">

pm25\_5yAvg

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

\-0.038<sup>\*\*</sup> (0.019)

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

pm10\_5yAvg

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

</td>

<td>

\-0.032<sup>\*\*</sup> (0.014)

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

no2\_5yAvg

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

</td>

<td>

</td>

<td>

0.020<sup>\*\*\*</sup> (0.007)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

nox\_5yAvg

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

</td>

<td>

</td>

<td>

</td>

<td>

0.012<sup>\*\*\*</sup> (0.004)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

o3\_5yAvg

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

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

\-0.256<sup>\*\*\*</sup> (0.041)

</td>

</tr>

<tr>

<td style="text-align:left">

Constant

</td>

<td>

10.872<sup>\*\*\*</sup> (0.629)

</td>

<td>

10.801<sup>\*\*\*</sup> (0.586)

</td>

<td>

7.974<sup>\*\*\*</sup> (0.605)

</td>

<td>

8.233<sup>\*\*\*</sup> (0.564)

</td>

<td>

9.129<sup>\*\*\*</sup> (0.480)

</td>

<td>

10.074<sup>\*\*\*</sup> (0.594)

</td>

<td>

10.149<sup>\*\*\*</sup> (0.586)

</td>

<td>

8.463<sup>\*\*\*</sup> (0.617)

</td>

<td>

8.493<sup>\*\*\*</sup> (0.601)

</td>

<td>

9.243<sup>\*\*\*</sup> (0.505)

</td>

</tr>

<tr>

<td colspan="11" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Observations

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

</tr>

<tr>

<td style="text-align:left">

Log Likelihood

</td>

<td>

\-1,867.267

</td>

<td>

\-1,865.012

</td>

<td>

\-1,869.058

</td>

<td>

\-1,869.567

</td>

<td>

\-1,861.414

</td>

<td>

\-1,873.401

</td>

<td>

\-1,872.673

</td>

<td>

\-1,870.542

</td>

<td>

\-1,870.358

</td>

<td>

\-1,856.444

</td>

</tr>

<tr>

<td style="text-align:left">

theta

</td>

<td>

2.217<sup>\*\*\*</sup> (0.176)

</td>

<td>

2.248<sup>\*\*\*</sup> (0.179)

</td>

<td>

2.193<sup>\*\*\*</sup> (0.174)

</td>

<td>

2.186<sup>\*\*\*</sup> (0.174)

</td>

<td>

2.300<sup>\*\*\*</sup> (0.184)

</td>

<td>

2.134<sup>\*\*\*</sup> (0.169)

</td>

<td>

2.143<sup>\*\*\*</sup> (0.170)

</td>

<td>

2.172<sup>\*\*\*</sup> (0.173)

</td>

<td>

2.175<sup>\*\*\*</sup> (0.173)

</td>

<td>

2.375<sup>\*\*\*</sup> (0.190)

</td>

</tr>

<tr>

<td style="text-align:left">

Akaike Inf. Crit.

</td>

<td>

3,744.535

</td>

<td>

3,740.024

</td>

<td>

3,748.115

</td>

<td>

3,749.134

</td>

<td>

3,732.829

</td>

<td>

3,756.802

</td>

<td>

3,755.346

</td>

<td>

3,751.084

</td>

<td>

3,750.716

</td>

<td>

3,722.888

</td>

</tr>

<tr>

<td colspan="11" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="10" style="text-align:right">

<sup>*</sup>p\<0.1; <sup>**</sup>p\<0.05; <sup>***</sup>p\<0.01

</td>

</tr>

</table>

##### Plot odds ratios

``` r
list_df_full = do.call("rbind", list(
                summary(LA_covid_pm25)$coefficients[nrow(summary(LA_covid_pm25)$coefficients),1:4],
                 summary(LA_covid_pm10)$coefficients[nrow(summary(LA_covid_pm10)$coefficients),1:4],
                 summary(LA_covid_no2)$coefficients[nrow(summary(LA_covid_no2)$coefficients),1:4],
                 summary(LA_covid_nox)$coefficients[nrow(summary(LA_covid_nox)$coefficients),1:4],
                 summary(LA_covid_o3)$coefficients[nrow(summary(LA_covid_o3)$coefficients),1:4],
                 summary(LA_covid_pm25_5yAvg)$coefficients[nrow(summary(LA_covid_pm25_5yAvg)$coefficients),1:4],
                 summary(LA_covid_pm10_5yAvg)$coefficients[nrow(summary(LA_covid_pm10_5yAvg)$coefficients),1:4],
                 summary(LA_covid_no2_5yAvg)$coefficients[nrow(summary(LA_covid_no2_5yAvg)$coefficients),1:4],
                 summary(LA_covid_nox_5yAvg)$coefficients[nrow(summary(LA_covid_nox_5yAvg)$coefficients),1:4],
                 summary(LA_covid_o3_5yAvg)$coefficients[nrow(summary(LA_covid_o3_5yAvg)$coefficients),1:4]))
#add labels 

pollutants = c("pm25","pm10","no2","nox","o3",
         "pm25_5yAvg","pm10_5yAvg","no2_5yAvg","nox_5yAvg","o3_5yAvg")
pollutants = data.frame(pollutants)

list_df_full = cbind(pollutants, list_df_full)

colnames(list_df_full) <- c("pollutants", "Estimate","Std", "z","p")

list_df_full$OR = exp(list_df_full$Estimate)
list_df_full$upper = exp(list_df_full$Estimate + list_df_full$Std)
list_df_full$lower = exp(list_df_full$Estimate - list_df_full$Std)

list_df_full$significance = "p-value<0.05"
list_df_full$significance[list_df_full$p > 0.05] <- "p-value>0.05"

write.csv(list_df_full, "data_output/LA_covidInfection_RRs_fullData_resub.csv", row.names=F)

## odds ratios and 95% CI

library(ggplot2)
ggplot(list_df_full, 
       aes(x=reorder(pollutants, OR), y=OR, color = significance)) + 
    geom_point(shape=21, size = 2) +
    geom_errorbar(aes(ymin=lower, ymax=upper),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_classic() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("Infectivity rate ratios, resub") + 
  xlab("pollutants")
```

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->

``` r
ggsave('fig_out/LA_covidInfection_RRs_fullData_resub.pdf', width = 6, height = 4)
```

### April deaths model

``` r
LA_covid_pm25 <-glm.nb(data = resub_dt, april_deaths ~ X2018.people.per.sq..km + Mean_ann_earnings + median_age_2018 +
                          pm25_val)

LA_covid_no2 <-glm.nb(data = resub_dt, april_deaths ~ X2018.people.per.sq..km + Mean_ann_earnings + median_age_2018 +
                          no2_val)

LA_covid_nox <-glm.nb(data = resub_dt, april_deaths ~ X2018.people.per.sq..km + Mean_ann_earnings + median_age_2018 +
                          nox_val)

LA_covid_pm10 <-glm.nb(data = resub_dt, april_deaths ~ X2018.people.per.sq..km + Mean_ann_earnings + median_age_2018 +
                          pm10_val)

LA_covid_o3 <-glm.nb(data = resub_dt, april_deaths ~ X2018.people.per.sq..km + Mean_ann_earnings + median_age_2018 +
                          o3_val)


LA_covid_pm25_5yAvg <-glm.nb(data = resub_dt, april_deaths ~ pop_dens_5ya + earnings_5ya + age_5ya +
                          pm25_5yAvg)

LA_covid_no2_5yAvg <-glm.nb(data = resub_dt, april_deaths ~ pop_dens_5ya + earnings_5ya + age_5ya +
                          no2_5yAvg)

LA_covid_nox_5yAvg <-glm.nb(data = resub_dt, april_deaths ~ pop_dens_5ya + earnings_5ya + age_5ya +
                          nox_5yAvg)

LA_covid_pm10_5yAvg <-glm.nb(data = resub_dt, april_deaths ~ pop_dens_5ya + earnings_5ya + age_5ya +
                          pm10_5yAvg)

LA_covid_o3_5yAvg <-glm.nb(data = resub_dt, april_deaths ~ pop_dens_5ya + earnings_5ya + age_5ya +
                          o3_5yAvg)

stargazer::stargazer(LA_covid_pm25,LA_covid_pm10,LA_covid_no2,LA_covid_nox,LA_covid_o3,
                     LA_covid_pm25_5yAvg,LA_covid_pm10_5yAvg,LA_covid_no2_5yAvg,
                         LA_covid_nox_5yAvg,LA_covid_o3_5yAvg,
                     type="html",out = "fig_out/LA_covid_NB_model_resub_deaths.html",
          dep.var.labels="COVID-19 death",
          single.row=TRUE)
```

<table style="text-align:center">

<tr>

<td colspan="11" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="10">

<em>Dependent variable:</em>

</td>

</tr>

<tr>

<td>

</td>

<td colspan="10" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

</td>

<td colspan="10">

COVID-19 death

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

<td>

(7)

</td>

<td>

(8)

</td>

<td>

(9)

</td>

<td>

(10)

</td>

</tr>

<tr>

<td colspan="11" style="border-bottom: 1px solid black">

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

0.00003 (0.00002)

</td>

<td>

0.00003 (0.00002)

</td>

<td>

0.00005<sup>\*\*</sup> (0.00002)

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

Mean\_ann\_earnings

</td>

<td>

0.00000 (0.00001)

</td>

<td>

0.00000 (0.00001)

</td>

<td>

\-0.00001 (0.00001)

</td>

<td>

\-0.00001 (0.00001)

</td>

<td>

0.00000 (0.00001)

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

median\_age\_2018

</td>

<td>

\-0.076<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.074<sup>\*\*\*</sup> (0.010)

</td>

<td>

\-0.043<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.046<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.059<sup>\*\*\*</sup> (0.010)

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

pm25\_val

</td>

<td>

\-0.060<sup>\*\*</sup> (0.029)

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

\-0.047<sup>\*\*</sup> (0.018)

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

0.031<sup>\*\*\*</sup> (0.009)

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

</td>

<td>

0.018<sup>\*\*\*</sup> (0.005)

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

\-0.042<sup>\*\*\*</sup> (0.011)

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

pop\_dens\_5ya

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

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.0001<sup>\*\*\*</sup> (0.00002)

</td>

<td>

0.00004<sup>\*\*</sup> (0.00002)

</td>

<td>

0.00004<sup>\*</sup> (0.00002)

</td>

<td>

0.00003 (0.00002)

</td>

</tr>

<tr>

<td style="text-align:left">

earnings\_5ya

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

\-0.00001 (0.00001)

</td>

<td>

\-0.00000 (0.00001)

</td>

<td>

\-0.00000 (0.00001)

</td>

<td>

\-0.00000 (0.00001)

</td>

<td>

0.00000 (0.00001)

</td>

</tr>

<tr>

<td style="text-align:left">

age\_5ya

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

\-0.069<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.069<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.046<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.046<sup>\*\*\*</sup> (0.011)

</td>

<td>

\-0.056<sup>\*\*\*</sup> (0.010)

</td>

</tr>

<tr>

<td style="text-align:left">

pm25\_5yAvg

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

0.005 (0.018)

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

pm10\_5yAvg

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

</td>

<td>

0.0004 (0.012)

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

no2\_5yAvg

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

</td>

<td>

</td>

<td>

0.025<sup>\*\*\*</sup> (0.006)

</td>

<td>

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

nox\_5yAvg

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

</td>

<td>

</td>

<td>

</td>

<td>

0.016<sup>\*\*\*</sup> (0.004)

</td>

<td>

</td>

</tr>

<tr>

<td style="text-align:left">

o3\_5yAvg

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

</td>

<td>

</td>

<td>

</td>

<td>

</td>

<td>

\-0.183<sup>\*\*\*</sup> (0.038)

</td>

</tr>

<tr>

<td style="text-align:left">

Constant

</td>

<td>

8.145<sup>\*\*\*</sup> (0.577)

</td>

<td>

8.170<sup>\*\*\*</sup> (0.540)

</td>

<td>

6.087<sup>\*\*\*</sup> (0.544)

</td>

<td>

6.310<sup>\*\*\*</sup> (0.508)

</td>

<td>

7.170<sup>\*\*\*</sup> (0.442)

</td>

<td>

7.466<sup>\*\*\*</sup> (0.538)

</td>

<td>

7.536<sup>\*\*\*</sup> (0.533)

</td>

<td>

6.150<sup>\*\*\*</sup> (0.548)

</td>

<td>

6.196<sup>\*\*\*</sup> (0.534)

</td>

<td>

7.254<sup>\*\*\*</sup> (0.467)

</td>

</tr>

<tr>

<td colspan="11" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

Observations

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

<td>

283

</td>

</tr>

<tr>

<td style="text-align:left">

Log Likelihood

</td>

<td>

\-1,535.037

</td>

<td>

\-1,534.045

</td>

<td>

\-1,530.920

</td>

<td>

\-1,531.456

</td>

<td>

\-1,530.728

</td>

<td>

\-1,537.505

</td>

<td>

\-1,537.548

</td>

<td>

\-1,528.701

</td>

<td>

\-1,528.580

</td>

<td>

\-1,526.697

</td>

</tr>

<tr>

<td style="text-align:left">

theta

</td>

<td>

2.683<sup>\*\*\*</sup> (0.220)

</td>

<td>

2.700<sup>\*\*\*</sup> (0.221)

</td>

<td>

2.759<sup>\*\*\*</sup> (0.227)

</td>

<td>

2.749<sup>\*\*\*</sup> (0.226)

</td>

<td>

2.761<sup>\*\*\*</sup> (0.227)

</td>

<td>

2.641<sup>\*\*\*</sup> (0.216)

</td>

<td>

2.641<sup>\*\*\*</sup> (0.216)

</td>

<td>

2.799<sup>\*\*\*</sup> (0.230)

</td>

<td>

2.801<sup>\*\*\*</sup> (0.230)

</td>

<td>

2.837<sup>\*\*\*</sup> (0.234)

</td>

</tr>

<tr>

<td style="text-align:left">

Akaike Inf. Crit.

</td>

<td>

3,080.074

</td>

<td>

3,078.090

</td>

<td>

3,071.840

</td>

<td>

3,072.912

</td>

<td>

3,071.456

</td>

<td>

3,085.011

</td>

<td>

3,085.096

</td>

<td>

3,067.402

</td>

<td>

3,067.160

</td>

<td>

3,063.394

</td>

</tr>

<tr>

<td colspan="11" style="border-bottom: 1px solid black">

</td>

</tr>

<tr>

<td style="text-align:left">

<em>Note:</em>

</td>

<td colspan="10" style="text-align:right">

<sup>*</sup>p\<0.1; <sup>**</sup>p\<0.05; <sup>***</sup>p\<0.01

</td>

</tr>

</table>

##### Plot odds ratios

``` r
list_df_full = do.call("rbind", list(
                summary(LA_covid_pm25)$coefficients[nrow(summary(LA_covid_pm25)$coefficients),1:4],
                 summary(LA_covid_pm10)$coefficients[nrow(summary(LA_covid_pm10)$coefficients),1:4],
                 summary(LA_covid_no2)$coefficients[nrow(summary(LA_covid_no2)$coefficients),1:4],
                 summary(LA_covid_nox)$coefficients[nrow(summary(LA_covid_nox)$coefficients),1:4],
                 summary(LA_covid_o3)$coefficients[nrow(summary(LA_covid_o3)$coefficients),1:4],
                 summary(LA_covid_pm25_5yAvg)$coefficients[nrow(summary(LA_covid_pm25_5yAvg)$coefficients),1:4],
                 summary(LA_covid_pm10_5yAvg)$coefficients[nrow(summary(LA_covid_pm10_5yAvg)$coefficients),1:4],
                 summary(LA_covid_no2_5yAvg)$coefficients[nrow(summary(LA_covid_no2_5yAvg)$coefficients),1:4],
                 summary(LA_covid_nox_5yAvg)$coefficients[nrow(summary(LA_covid_nox_5yAvg)$coefficients),1:4],
                 summary(LA_covid_o3_5yAvg)$coefficients[nrow(summary(LA_covid_o3_5yAvg)$coefficients),1:4]))
#add labels 

pollutants = c("pm25","pm10","no2","nox","o3",
         "pm25_5yAvg","pm10_5yAvg","no2_5yAvg","nox_5yAvg","o3_5yAvg")
pollutants = data.frame(pollutants)

list_df_full = cbind(pollutants, list_df_full)

colnames(list_df_full) <- c("pollutants", "Estimate","Std", "z","p")

list_df_full$OR = exp(list_df_full$Estimate)
list_df_full$upper = exp(list_df_full$Estimate + list_df_full$Std)
list_df_full$lower = exp(list_df_full$Estimate - list_df_full$Std)

list_df_full$significance = "p-value<0.05"
list_df_full$significance[list_df_full$p > 0.05] <- "p-value>0.05"

write.csv(list_df_full, "data_output/LA_covidDEATH_RRs_fullData_resub.csv", row.names=F)

## odds ratios and 95% CI

library(ggplot2)
ggplot(list_df_full, 
       aes(x=reorder(pollutants, OR), y=OR, color = significance)) + 
    geom_point(shape=21, size = 2) +
    geom_errorbar(aes(ymin=lower, ymax=upper),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_classic() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("Mortality rate ratios, resub") + 
  xlab("pollutants")
```

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-56-1.png)<!-- -->

``` r
ggsave('fig_out/LA_covidDEATH_RRs_fullData_resub.pdf', width = 6, height = 4)
```

## Aim 3: Individual level analysis

### Data curation

#### covid data curation

``` r
covid_df = read.csv("data/covid19_result.txt", sep = '\t')
covid_id = covid_df[,c(1,5,6)]
colnames(covid_id)
```

    ## [1] "eid"    "origin" "result"

``` r
covid_id_unique = aggregate(covid_id, by=list(covid_id$eid), 
                            FUN=max)

ggplot(data = covid_id_unique, aes(x = result))+
  geom_histogram(stat="count")
```

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-57-1.png)<!-- -->

``` r
ggplot(data = covid_id_unique, aes(x = origin))+
  geom_histogram(stat="count")
```

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-57-2.png)<!-- -->
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

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->

#### load non-covid UKB phenotype data

Data curation

``` r
### load non-covid UKB phenotype data ----
ukb_covid = read.csv("data/29_4_2020_ukb41646_covid19_subset.csv")[,-1]
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
UID_names = read.csv("data/UKB_list_columns_AIRcovid.csv")[c("FieldID","my_colname")]
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

write.csv(ukb_covid_merge, "data_output/processed_29_4_2020_ukb41646_covid19.csv")
```

#### Preliminary analysis with PM2.5

I will load the air pollutants & match the covid patients to their
nearest pollution climate mapping location.

``` r
### load pm2.5 data ----
#note the pop weighted data have the area code, but not the x y coords

pm25 = read.csv("data/processed_PM25_uk-air_annual_mean_mappm252018g.csv", na.strings = "MISSING")[c("x","y","pm252018g")]
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


write.csv(pm25_ll_df,"data_output/processed_pm25_lonlat.csv")
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
write.csv(ukb_covid_ll_df, ("data_output/ukb_covid_lonlat_df.csv"))
# plot UKB locations
world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf() +
  geom_point(data = ukb_covid_ll_df, aes(x = ukb_lon, y = ukb_lat), size = 1, 
             shape = 23, fill = "darkred") +
  coord_sf(ylim = c(min(ukb_covid_ll_df$ukb_lat)-2, max(ukb_covid_ll_df$ukb_lat)+4), 
           xlim = c(min(ukb_covid_ll_df$ukb_lon)-4, max(ukb_covid_ll_df$ukb_lon)+3), expand = FALSE)
```

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-62-1.png)<!-- -->

``` r
ggsave('fig_out/UKB_COVID_participants_locations.pdf')
```

**I wrote a python script: match\_air\_to\_UKB\_covid**

### Merge and QC data for PM2.5

``` r
ukb_eid_pm25_merged = read.csv('data_output/merged_ukb_pm25.csv')[c("eid","distance","pm25_val")]

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
ukb_cities_eid = read.csv("data_output/2_5_2020_ukb_eid_match_cities.csv")[c("eid","spec_area","gen_area")]
nrow(ukb_cities_eid)
```

    ## [1] 1464

``` r
popDens_2018 = read.csv("data/2018_official_popDensity.csv")[c("Name","X2018.people.per.sq..km")]
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
ukb_covid_pm25_popDens_df = read.csv("data_output/2_5_2020_full_cov_CV_UKB_air_dataset.csv")
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
no2_raw_dt = na.omit(read.csv("data/processed_30_4_2020_NO2_map2018g.csv", na.strings = "MISSING")[c("x","y","no22018")])

ukgrid <- "+init=epsg:27700"
latlong <- "+init=epsg:4326"

### Create coordinates variable
no2_coords <- cbind(Easting = as.numeric(as.character(no2_raw_dt$x)),
                     Northing = as.numeric(as.character(no2_raw_dt$y)))

no2_LL <- spTransform(SpatialPointsDataFrame(no2_coords,
                                  data = no2_raw_dt,
                                  proj4string = CRS("+init=epsg:27700")), CRS(latlong))

no2_LL_df = data.frame('no2_lon' = coordinates(no2_LL)[, 1], 'no2_lat' = coordinates(no2_LL)[, 2], 'no2_val' = no2_LL$no22018)

write.csv(no2_LL_df,"data_output/processed_no2_lonlat.csv")

### convert SO2 x y to lat lon ----
so2_raw_dt = na.omit(read.csv("data/processed_30_4_2020_SO2_map2018g.csv", na.strings = "MISSING")[c("x","y","so22018")])

so2_coords <- cbind(Easting = as.numeric(as.character(so2_raw_dt$x)),
                    Northing = as.numeric(as.character(so2_raw_dt$y)))

so2_LL <- spTransform(SpatialPointsDataFrame(so2_coords,
                                             data = so2_raw_dt,
                                             proj4string = CRS("+init=epsg:27700")), CRS(latlong))

so2_LL_df = data.frame('so2_lon' = coordinates(so2_LL)[, 1], 'so2_lat' = coordinates(so2_LL)[, 2], 'so2_val' = so2_LL$so22018)

write.csv(so2_LL_df,"data_output/processed_so2_lonlat.csv")

### convert 03 x y to lat lon ----

o3_raw_dt = na.omit(read.csv("data/processed_30_4_2020_O3_map2018g.csv", na.strings = "MISSING")[c("x","y","dgt12018")])

o3_coords <- cbind(Easting = as.numeric(as.character(o3_raw_dt$x)),
                    Northing = as.numeric(as.character(o3_raw_dt$y)))

o3_LL <- spTransform(SpatialPointsDataFrame(o3_coords,
                                             data = o3_raw_dt,
                                             proj4string = CRS("+init=epsg:27700")), CRS(latlong))

o3_LL_df = data.frame('o3_lon' = coordinates(o3_LL)[, 1], 'o3_lat' = coordinates(o3_LL)[, 2], 'o3_val' = o3_LL$dgt12018)

write.csv(o3_LL_df,"data_output/processed_o3_lonlat.csv")

### convert PM10 x y to lat lon ----

pm10_raw_dt = na.omit(read.csv("data/processed_30_4_2020_PM10_mappm102018g.csv", na.strings = "MISSING")[c("x","y","pm102018g")])

pm10_coords <- cbind(Easting = as.numeric(as.character(pm10_raw_dt$x)),
                   Northing = as.numeric(as.character(pm10_raw_dt$y)))

pm10_LL <- spTransform(SpatialPointsDataFrame(pm10_coords,
                                            data = pm10_raw_dt,
                                            proj4string = CRS("+init=epsg:27700")), CRS(latlong))

pm10_LL_df = data.frame('pm10_lon' = coordinates(pm10_LL)[, 1], 'pm10_lat' = coordinates(pm10_LL)[, 2], 'pm10_val' = pm10_LL$pm102018g)

write.csv(pm10_LL_df,"data_output/processed_pm10_lonlat.csv")


### convert NOx x y to lat lon ----

nox_raw_dt = na.omit(read.csv("data/processed_30_4_2020_NOX_map2018g.csv", na.strings = "MISSING")[c("x","y","nox2018")])

nox_coords <- cbind(Easting = as.numeric(as.character(nox_raw_dt$x)),
                     Northing = as.numeric(as.character(nox_raw_dt$y)))

nox_LL <- spTransform(SpatialPointsDataFrame(nox_coords,
                                              data = nox_raw_dt,
                                              proj4string = CRS("+init=epsg:27700")), CRS(latlong))

nox_LL_df = data.frame('nox_lon' = coordinates(nox_LL)[, 1], 'nox_lat' = coordinates(nox_LL)[, 2], 'nox_val' = nox_LL$nox2018)

write.csv(nox_LL_df,"data_output/processed_nox_lonlat.csv")
```

### Analysis

load and merge all data

``` r
no2_ukb = read.csv("data_output/merged_ukb_no2.csv")[c('eid','no2_val')]
so2_ukb = read.csv("data_output/merged_ukb_so2.csv")[c('eid','so2_val')]
o3_ukb = read.csv("data_output/merged_ukb_o3.csv")[c('eid','o3_val')]
nox_ukb = read.csv("data_output/merged_ukb_nox.csv")[c('eid','nox_val')]
pm10_ukb = read.csv("data_output/merged_ukb_pm10.csv")[c('eid','pm10_val')]

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
          ukb_covid_o3.b, ukb_covid_so2.b,type="html",out = "fig_out/ukb_covid_allPoll_binary.html",
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

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-73-1.png)<!-- -->

``` r
ggsave('fig_out/UKB_infectivityOR_bar.pdf')
```

``` r
stargazer(ukb_covid_onlyPoll, type ="html", single.row=TRUE, summary = FALSE, out = "fig_out/ukb_covid_onlyPoll.html")
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

\#\#\# UKB analysis with 5YA

\#\#\#\# Adding the 5YA air pollution data

``` r
ukb_dt = read.csv("data_output/2_5_2020_full_cov_CV_UKB_air_dataset.csv")
ukb_dt= ukb_dt[ , -which(names(ukb_dt) %in% c("pm25_val"))]
ukb_5yAP = read.csv("data_output/ukb_covid_dt_out.csv")[,c("eid","no2_5yAvg","nox_5yAvg","pm10_5yAvg","o3_5yAvg","pm25_5yAvg","no2_val",    "nox_val",  "pm10_val", "o3_val",   "pm25_val")]
ukb_dt = merge(ukb_dt,ukb_5yAP, by = "eid", all.x = T)
nrow(ukb_dt)
```

    ## [1] 1464

#### models - cases only

``` r
ukb_covid.b_pm25 <-glm(data = ukb_dt, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+whr + sex + age+
                          pm25_val, family = 'binomial')

ukb_covid.b_no2 <-glm(data = ukb_dt, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+whr + sex + age+
                          no2_val, family = 'binomial')

ukb_covid.b_nox <-glm(data = ukb_dt, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+whr + sex + age+
                          nox_val, family = 'binomial')

ukb_covid.b_pm10 <-glm(data = ukb_dt, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+whr + sex + age+
                          pm10_val, family = 'binomial')

ukb_covid.b_o3 <-glm(data = ukb_dt, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+whr + sex + age+
                          o3_val, family = 'binomial')


ukb_covid.b_pm25_5yAvg <-glm(data = ukb_dt, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+whr + sex + age+
                          pm25_5yAvg, family = 'binomial')

ukb_covid.b_no2_5yAvg <-glm(data = ukb_dt, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+whr + sex + age+
                          no2_5yAvg, family = 'binomial')

ukb_covid.b_nox_5yAvg <-glm(data = ukb_dt, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+whr + sex + age+
                          nox_5yAvg, family = 'binomial')

ukb_covid.b_pm10_5yAvg <-glm(data = ukb_dt, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+whr + sex + age+
                          pm10_5yAvg, family = 'binomial')

ukb_covid.b_o3_5yAvg <-glm(data = ukb_dt, result ~ merged_popDens_km2 + n_cancers +townsend + smoking + whistling + diabetes+whr + sex + age+
                          o3_5yAvg, family = 'binomial')

stargazer::stargazer(ukb_covid.b_pm25,ukb_covid.b_pm10,ukb_covid.b_no2,ukb_covid.b_nox,ukb_covid.b_o3,
                     ukb_covid.b_pm25_5yAvg,ukb_covid.b_pm10_5yAvg,ukb_covid.b_no2_5yAvg,
                         ukb_covid.b_nox_5yAvg,ukb_covid.b_o3_5yAvg,
                     type="html",out = "fig_out/ukb_covidInfections_binomial_resub.html",
          dep.var.labels="COVID-19 positive or not",
          single.row=TRUE)
```

    ## 
    ## <table style="text-align:center"><tr><td colspan="11" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="10"><em>Dependent variable:</em></td></tr>
    ## <tr><td></td><td colspan="10" style="border-bottom: 1px solid black"></td></tr>
    ## <tr><td style="text-align:left"></td><td colspan="10">COVID-19 positive or not</td></tr>
    ## <tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td><td>(3)</td><td>(4)</td><td>(5)</td><td>(6)</td><td>(7)</td><td>(8)</td><td>(9)</td><td>(10)</td></tr>
    ## <tr><td colspan="11" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">merged_popDens_km2</td><td>-0.00002 (0.00003)</td><td>-0.00002 (0.00003)</td><td>-0.00004 (0.00003)</td><td>-0.00003 (0.00003)</td><td>0.00003 (0.00002)</td><td>-0.00001 (0.00003)</td><td>-0.00001 (0.00003)</td><td>-0.00003 (0.00003)</td><td>-0.00003 (0.00003)</td><td>0.00003 (0.00002)</td></tr>
    ## <tr><td style="text-align:left">n_cancers</td><td>-0.229 (0.166)</td><td>-0.227 (0.166)</td><td>-0.214 (0.166)</td><td>-0.216 (0.166)</td><td>-0.228 (0.166)</td><td>-0.227 (0.166)</td><td>-0.226 (0.165)</td><td>-0.216 (0.166)</td><td>-0.219 (0.166)</td><td>-0.224 (0.166)</td></tr>
    ## <tr><td style="text-align:left">townsend</td><td>0.025 (0.018)</td><td>0.024 (0.018)</td><td>0.008 (0.018)</td><td>0.010 (0.018)</td><td>0.024 (0.019)</td><td>0.023 (0.018)</td><td>0.023 (0.018)</td><td>0.008 (0.018)</td><td>0.010 (0.018)</td><td>0.020 (0.019)</td></tr>
    ## <tr><td style="text-align:left">smokingNever</td><td>0.635<sup>***</sup> (0.176)</td><td>0.627<sup>***</sup> (0.176)</td><td>0.614<sup>***</sup> (0.176)</td><td>0.612<sup>***</sup> (0.176)</td><td>0.611<sup>***</sup> (0.176)</td><td>0.624<sup>***</sup> (0.176)</td><td>0.619<sup>***</sup> (0.176)</td><td>0.609<sup>***</sup> (0.176)</td><td>0.608<sup>***</sup> (0.176)</td><td>0.603<sup>***</sup> (0.175)</td></tr>
    ## <tr><td style="text-align:left">smokingPrefer not to answer</td><td>1.520<sup>**</sup> (0.746)</td><td>1.505<sup>**</sup> (0.746)</td><td>1.411<sup>*</sup> (0.742)</td><td>1.413<sup>*</sup> (0.742)</td><td>1.470<sup>**</sup> (0.743)</td><td>1.500<sup>**</sup> (0.746)</td><td>1.489<sup>**</sup> (0.745)</td><td>1.408<sup>*</sup> (0.743)</td><td>1.411<sup>*</sup> (0.742)</td><td>1.435<sup>*</sup> (0.742)</td></tr>
    ## <tr><td style="text-align:left">smokingPrevious</td><td>0.703<sup>***</sup> (0.177)</td><td>0.690<sup>***</sup> (0.176)</td><td>0.682<sup>***</sup> (0.176)</td><td>0.677<sup>***</sup> (0.176)</td><td>0.668<sup>***</sup> (0.176)</td><td>0.690<sup>***</sup> (0.176)</td><td>0.682<sup>***</sup> (0.176)</td><td>0.677<sup>***</sup> (0.176)</td><td>0.674<sup>***</sup> (0.176)</td><td>0.658<sup>***</sup> (0.176)</td></tr>
    ## <tr><td style="text-align:left">whistlingNo</td><td>0.475 (0.349)</td><td>0.461 (0.349)</td><td>0.462 (0.351)</td><td>0.461 (0.351)</td><td>0.427 (0.348)</td><td>0.466 (0.349)</td><td>0.458 (0.349)</td><td>0.462 (0.351)</td><td>0.462 (0.351)</td><td>0.422 (0.348)</td></tr>
    ## <tr><td style="text-align:left">whistlingPrefer not to answer</td><td>0.123 (1.637)</td><td>0.127 (1.636)</td><td>0.031 (1.629)</td><td>0.053 (1.629)</td><td>0.082 (1.637)</td><td>0.118 (1.633)</td><td>0.126 (1.632)</td><td>0.014 (1.629)</td><td>0.038 (1.629)</td><td>0.070 (1.635)</td></tr>
    ## <tr><td style="text-align:left">whistlingYes</td><td>0.458 (0.356)</td><td>0.446 (0.356)</td><td>0.450 (0.357)</td><td>0.447 (0.358)</td><td>0.400 (0.354)</td><td>0.448 (0.356)</td><td>0.442 (0.356)</td><td>0.448 (0.357)</td><td>0.445 (0.358)</td><td>0.396 (0.354)</td></tr>
    ## <tr><td style="text-align:left">diabetesNo</td><td>-0.275 (0.859)</td><td>-0.296 (0.858)</td><td>-0.343 (0.861)</td><td>-0.351 (0.859)</td><td>-0.317 (0.852)</td><td>-0.293 (0.858)</td><td>-0.311 (0.857)</td><td>-0.341 (0.860)</td><td>-0.350 (0.858)</td><td>-0.323 (0.852)</td></tr>
    ## <tr><td style="text-align:left">diabetesPrefer not to answer</td><td>-1.368 (1.912)</td><td>-1.377 (1.910)</td><td>-1.422 (1.907)</td><td>-1.407 (1.906)</td><td>-1.338 (1.909)</td><td>-1.466 (1.909)</td><td>-1.462 (1.908)</td><td>-1.435 (1.906)</td><td>-1.428 (1.906)</td><td>-1.344 (1.907)</td></tr>
    ## <tr><td style="text-align:left">diabetesYes</td><td>-0.324 (0.871)</td><td>-0.345 (0.871)</td><td>-0.389 (0.874)</td><td>-0.396 (0.872)</td><td>-0.379 (0.865)</td><td>-0.337 (0.870)</td><td>-0.356 (0.870)</td><td>-0.380 (0.872)</td><td>-0.389 (0.870)</td><td>-0.381 (0.865)</td></tr>
    ## <tr><td style="text-align:left">whr</td><td>1.487<sup>*</sup> (0.808)</td><td>1.511<sup>*</sup> (0.808)</td><td>1.555<sup>*</sup> (0.808)</td><td>1.595<sup>**</sup> (0.808)</td><td>1.580<sup>**</sup> (0.806)</td><td>1.492<sup>*</sup> (0.808)</td><td>1.508<sup>*</sup> (0.808)</td><td>1.547<sup>*</sup> (0.808)</td><td>1.582<sup>*</sup> (0.808)</td><td>1.588<sup>**</sup> (0.806)</td></tr>
    ## <tr><td style="text-align:left">sexMale</td><td>0.059 (0.138)</td><td>0.057 (0.138)</td><td>0.034 (0.138)</td><td>0.030 (0.138)</td><td>0.061 (0.138)</td><td>0.050 (0.138)</td><td>0.050 (0.138)</td><td>0.031 (0.138)</td><td>0.028 (0.138)</td><td>0.056 (0.138)</td></tr>
    ## <tr><td style="text-align:left">age</td><td>-0.011 (0.007)</td><td>-0.011<sup>*</sup> (0.007)</td><td>-0.010 (0.007)</td><td>-0.011 (0.007)</td><td>-0.010 (0.007)</td><td>-0.011<sup>*</sup> (0.007)</td><td>-0.011<sup>*</sup> (0.007)</td><td>-0.011 (0.007)</td><td>-0.011 (0.007)</td><td>-0.010 (0.007)</td></tr>
    ## <tr><td style="text-align:left">pm25_val</td><td>0.120<sup>***</sup> (0.040)</td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr>
    ## <tr><td style="text-align:left">pm10_val</td><td></td><td>0.076<sup>***</sup> (0.028)</td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr>
    ## <tr><td style="text-align:left">no2_val</td><td></td><td></td><td>0.049<sup>***</sup> (0.015)</td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr>
    ## <tr><td style="text-align:left">nox_val</td><td></td><td></td><td></td><td>0.024<sup>***</sup> (0.009)</td><td></td><td></td><td></td><td></td><td></td><td></td></tr>
    ## <tr><td style="text-align:left">o3_val</td><td></td><td></td><td></td><td></td><td>0.012 (0.021)</td><td></td><td></td><td></td><td></td><td></td></tr>
    ## <tr><td style="text-align:left">pm25_5yAvg</td><td></td><td></td><td></td><td></td><td></td><td>0.120<sup>***</sup> (0.044)</td><td></td><td></td><td></td><td></td></tr>
    ## <tr><td style="text-align:left">pm10_5yAvg</td><td></td><td></td><td></td><td></td><td></td><td></td><td>0.077<sup>**</sup> (0.030)</td><td></td><td></td><td></td></tr>
    ## <tr><td style="text-align:left">no2_5yAvg</td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td>0.044<sup>***</sup> (0.014)</td><td></td><td></td></tr>
    ## <tr><td style="text-align:left">nox_5yAvg</td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td>0.021<sup>***</sup> (0.008)</td><td></td></tr>
    ## <tr><td style="text-align:left">o3_5yAvg</td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td>-0.019 (0.079)</td></tr>
    ## <tr><td style="text-align:left">Constant</td><td>-2.577<sup>**</sup> (1.203)</td><td>-2.502<sup>**</sup> (1.205)</td><td>-2.197<sup>*</sup> (1.175)</td><td>-2.006<sup>*</sup> (1.167)</td><td>-1.701 (1.160)</td><td>-2.604<sup>**</sup> (1.213)</td><td>-2.532<sup>**</sup> (1.212)</td><td>-2.203<sup>*</sup> (1.175)</td><td>-2.009<sup>*</sup> (1.167)</td><td>-1.581 (1.171)</td></tr>
    ## <tr><td colspan="11" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>1,447</td><td>1,447</td><td>1,447</td><td>1,447</td><td>1,447</td><td>1,447</td><td>1,447</td><td>1,447</td><td>1,447</td><td>1,447</td></tr>
    ## <tr><td style="text-align:left">Log Likelihood</td><td>-973.678</td><td>-974.683</td><td>-972.953</td><td>-974.207</td><td>-978.070</td><td>-974.439</td><td>-975.008</td><td>-973.235</td><td>-974.345</td><td>-978.204</td></tr>
    ## <tr><td style="text-align:left">Akaike Inf. Crit.</td><td>1,981.357</td><td>1,983.366</td><td>1,979.906</td><td>1,982.413</td><td>1,990.139</td><td>1,982.878</td><td>1,984.015</td><td>1,980.469</td><td>1,982.690</td><td>1,990.408</td></tr>
    ## <tr><td colspan="11" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="10" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
    ## </table>

``` r
list_df_full = do.call("rbind", list(
                summary(ukb_covid.b_pm25)$coefficients[nrow(summary(ukb_covid.b_pm25)$coefficients),1:4],
                 summary(ukb_covid.b_pm10)$coefficients[nrow(summary(ukb_covid.b_pm10)$coefficients),1:4],
                 summary(ukb_covid.b_no2)$coefficients[nrow(summary(ukb_covid.b_no2)$coefficients),1:4],
                 summary(ukb_covid.b_nox)$coefficients[nrow(summary(ukb_covid.b_nox)$coefficients),1:4],
                 summary(ukb_covid.b_o3)$coefficients[nrow(summary(ukb_covid.b_o3)$coefficients),1:4],
                 summary(ukb_covid.b_pm25_5yAvg)$coefficients[nrow(summary(ukb_covid.b_pm25_5yAvg)$coefficients),1:4],
                 summary(ukb_covid.b_pm10_5yAvg)$coefficients[nrow(summary(ukb_covid.b_pm10_5yAvg)$coefficients),1:4],
                 summary(ukb_covid.b_no2_5yAvg)$coefficients[nrow(summary(ukb_covid.b_no2_5yAvg)$coefficients),1:4],
                 summary(ukb_covid.b_nox_5yAvg)$coefficients[nrow(summary(ukb_covid.b_nox_5yAvg)$coefficients),1:4],
                 summary(ukb_covid.b_o3_5yAvg)$coefficients[nrow(summary(ukb_covid.b_o3_5yAvg)$coefficients),1:4]))
#add labels 

pollutants = c("pm25","pm10","no2","nox","o3",
         "pm25_5yAvg","pm10_5yAvg","no2_5yAvg","nox_5yAvg","o3_5yAvg")
pollutants = data.frame(pollutants)

list_df_full = cbind(pollutants, list_df_full)

colnames(list_df_full) <- c("pollutants", "Estimate","Std", "z","p")

list_df_full$OR = exp(list_df_full$Estimate)
list_df_full$upper = exp(list_df_full$Estimate + list_df_full$Std)
list_df_full$lower = exp(list_df_full$Estimate - list_df_full$Std)

list_df_full$significance = "p-value<0.05"
list_df_full$significance[list_df_full$p > 0.05] <- "p-value>0.05"

write.csv(list_df_full, "data_output/UKB_covidInfection_IORs_resub.csv", row.names=F)

## odds ratios and 95% CI

library(ggplot2)
ggplot(list_df_full, 
       aes(x=reorder(pollutants, OR), y=OR, color = significance)) + 
    geom_point(shape=21, size = 2) +
    geom_errorbar(aes(ymin=lower, ymax=upper),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
  theme_classic() + 
  geom_hline(yintercept = 1, linetype="dotted") +
  coord_flip()+ylab("Infectivity odds ratios, resub") + 
  xlab("pollutants")
```

![](Analysis_Workflow_COVID_air_notebook_files/figure-gfm/unnamed-chunk-77-1.png)<!-- -->

``` r
ggsave('fig_out/UKB_covidInfection_IORs_resub.pdf', width = 6, height = 4)
```
