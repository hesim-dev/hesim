---
title: Personalized cost-effectiveness analysis (pCEA) with `hesim`
author: "Devin Incerti"
date: "January 30, 2017"
output: 
  html_document:
bibliography: cea-references.bib
---

* TOC
{:toc}


# Overview
Personalized cost-effectiveness analysis (pCEA) evaluates the cost-effectiveness of treatments at the individual (or subgroup) level. This has two major implications:

* Optimal treatments will vary across patients.
* Treatments will be more cost-effective for some patients than others.

The `hesim` package help facilitate pCEA by providing a number of functions for analyzing subugroup level health and cost outcomes from simulation models that quantify parameter uncertainty using probabilistic sensitivity analysis (PSA). These functions take simulation output and generate measures commonly used for technology assessment including:

* net benefits.
* incremental cost-effectiveness ratios (ICERs).
* cost-effectivenss acceptability curves (CEACs).
* the expected value of perfect information (EVPI).

The rest of this document provides an overview of pCEA and how it can be conducted using `hesim`. The perspective is Bayesian in nature, in that it is concerned with estimating the entire distribution of outcomes rather than just expected values [@baio2012bayesian; @baio2015probabilistic]. It also stresses that both optimal treatments and the cost-effectiveness of those treatments vary accross individuals [@basu2007value; @espinoza2014value].

# Net Benefits
Decision analysis provides a formal framework for making treatment decisions based on the utility that a therapy provides to a patient population. Traditionally, the optimal treatment arm is the one that maximizes expected net benefits. The expected net benefit is calculated by averaging over the patient population and uncertain parameters $\theta$. For a given subgroup $g$ and parameter set $\theta$, these net benefits are computed as the difference between the monetized health gains from an intervention less costs, or,

$$
\begin{aligned}
NB_g(j,\theta) = e_{gj}\cdot k- c_{gj},
\end{aligned}
$$

where $e_{gj}$ and $c_{gj}$ are measures of clinical effectiveness (e.g. QALYs) and costs in subgroup $g$ using treatment $j$ respectively, and $k$ is a decision makers willingness to pay per unit of clinical effectiveness. The optimal treatment for a given subgroup is the one that maximizes expected net benefits,

$$
\begin{aligned}
j^{*}_g = \text{argmax}_j E_{\theta} \left[NB_g(j,\theta)\right].
\end{aligned}
$$

In practice, new interventions are usually compared to a standard treatment often referred to as the comparator. In these cases, a new treatment in a given subgroup is preferred to the comparator if the expected incremental net benefit of the new treatment is positive; that is, treatment 1 is preferred to treatment 0 in subgroup $g$ if $EINB_g > 0$ where the incremental net benefit (INB) is given by

$$
\begin{aligned}
INB_g = NB_g(j = 1, \theta)] - NB_g(j = 0, \theta),
\end{aligned}
$$
and $EINB_g =E_\theta \left[INB_g\right]$. Equivalently, treatment $1$ is preferred to treatment $0$ in subgroup $g$ if the incremental cost-effectiveness ratio (ICER) is greater than the willingness to pay $k$,

$$
\begin{aligned}
k > \frac{c_{g1} - c_{s0}}{e_{g1} - e_{g0}} = ICER_g.
\end{aligned}
$$


# Probabilistic sensitivity analysis
Expected net benefits are based entirely on expected values and ignore parameter uncertainty. This implies that net benefits are uncertain and that optimal treatment arms may be selected incorrectly. This uncertainty can be quantified using PSA, which uses Bayesian and quasi-Bayesian techniques to estimate the distribution of net benefits given the distribution of the parameters for each treatment arm.

Since the joint distribution of the model parameters cannot be derived analytically (except in the simplest cases), the distribution of $\theta$ is approximated by simulating the parameters from their joint posterior distribution and calculating relevant quantities of interest as a function of the simulated parameters. For each treatment arm and subgroup, PSA therefore produces $n$ random draws from the posterior distribution of clinical effectiveness and costs,

$$
\begin{aligned}
e_{gj} &= [e_{gj}^1, e_{gj}^2, \dots, e_{gj}^n] \\
c_{gj} &= [c_{gj}^1, c_{gj}^2, \dots, c_{gj}^n].
\end{aligned}
$$

Below we simulate costs and QALYs for three treatment arms and two subgroups (in a real world analysis, this output would be derived from a detailed health-economic simulation model). Arm 1 is the current standard of care; it is the cheapest therapy, but also the least efficacious. Arms 2 and 3 are equally costly, but Arm 2 is more effective in subgroup 1 while Arm 3 is more effective in subgroup 2. 


```r
nsims <- 1000

# cost
c <- vector(mode = "list", length = 6)
names(c) <- c("Arm 1, Grp 1", "Arm 1, Grp 2", "Arm 2, Grp 1",
              "Arm 2, Grp 2", "Arm 3, Grp 1", "Arm 3, Grp 2")
c[[1]] <- rlnorm(nsims, 2, .1)
c[[2]] <- rlnorm(nsims, 2, .1)
c[[3]] <- rlnorm(nsims, 11, .15)
c[[4]] <- rlnorm(nsims, 11, .15)
c[[5]] <- rlnorm(nsims, 11, .15)
c[[6]] <- rlnorm(nsims, 11, .15)

# effectiveness
e <- c
e[[1]] <- rnorm(nsims, 8, .2)
e[[2]] <- rnorm(nsims, 8, .2)
e[[3]] <- rnorm(nsims, 10, .8)
e[[4]] <- rnorm(nsims, 10.5, .8)
e[[5]] <- rnorm(nsims, 8.5, .6)
e[[6]] <- rnorm(nsims, 11, .6)

# cost and effectiveness by arm and simulation
library("data.table")
ce <- data.table(sim = rep(seq(nsims), length(e)),
                             arm = rep(paste0("Arm ", seq(1, 3)), 
                                       each = nsims * 2),
                             grp = rep(rep(c("Group 1", "Group 2"),
                                           each = nsims), 3),
                             cost = do.call("c", c), qalys = do.call("c", e))
head(ce)
```

```
##    sim   arm     grp     cost    qalys
## 1:   1 Arm 1 Group 1 6.735913 8.449514
## 2:   2 Arm 1 Group 1 7.543181 8.266314
## 3:   3 Arm 1 Group 1 7.739197 8.009014
## 4:   4 Arm 1 Group 1 6.569844 7.756109
## 5:   5 Arm 1 Group 1 7.126620 7.830618
## 6:   6 Arm 1 Group 1 6.950060 8.000831
```

For any given willingness to pay $k$, expected net benefits can be calculated by arm, subgroup, and simulation number. For example, with $k=150,000$, a reasonable estimate of the value of a life-year in the United States, Arm 2 provides the most expected net benefits in subgroup 2 while Arm 3 provides the most expected net benefits in subgroup 2.


```r
ce[, nb := 150000 * qalys - cost]
enb <- ce[, .(enb = mean(nb)), by = c("arm", "grp")]
enb <- dcast(enb, arm ~ grp, value.var = "enb")
print(enb)
```

```
##      arm Group 1 Group 2
## 1: Arm 1 1199431 1198292
## 2: Arm 2 1440149 1511732
## 3: Arm 3 1214023 1591438
```

A number of measures have been proposed in the health economics literature to summarize the uncertainty estimated using PSA. Below we describe the most common measures, which are estimated using the functions `psa` and `psa_pw`. The `psa` function summarizes results by taking into account each treatment arm in the analysis, while the function `psa_pw` summarizes "pairwise" results in which each treatment is compared to a comparator. 

We can use the  `psa` function to summarize results from our `data.table` object of simulated output for a range of willingness to pay values,


```r
library("hesim")
ktop <- 200000
psa.dt <-  psa(ce, k = seq(0, ktop, 500), sim = "sim", arm = "arm",
              grp = "grp", e = "qalys", c = "cost")
```

The most important input in `psa` is the `data.table` object (`x`) containing columns for simulation number (`sim`), treatment arm (`arm`), subgroup (`grp`), clinical effectiveness (`e`), and costs (`c`). Users specify the names of the relevant columns in their output table as strings. The other relevant parameter is $k$, which is a range of willingness to pay values to use for estimating net benefits. 

Likewise, we can use `psa_pw` to summarize the PSA when directly comparing the two treatment arms (Arm 2 and Arm 3) to the comparator (Arm 1). 


```r
psa.pw.dt <-  psa_pw(ce,  k = seq(0, ktop, 500), control = "Arm 1",
                     sim = "sim", arm = "arm", e = "qalys", c = "cost")
```

The same inputs are used as in `psa` except users must specify the name of the control arm. 

## Probability most cost-effective
A useful summary measure for quantifying uncertainty is the probability that each treatment arm is the most cost effective. For a particular subgroup, this is estimated from simulation output as the proportion of simulation draws that each arm has the highest net benefits. For example, consider a random sample of 10 draws from the PSA simulation output and suppose $k$ is again equal to $150,000$. 


```r
library("knitr")
ce.nb <- dcast(ce[sim %in% sample(1:nsims, 10) & grp == "Group 2"], 
               sim ~ arm, value.var = "nb")
setnames(ce.nb, colnames(ce.nb), c("sim", "nb1", "nb2", "nb3"))
ce.nb[, maxj := apply(ce.nb[, .(nb1, nb2, nb3)], 1, which.max)]
ce.nb[, maxj := factor(maxj, levels = c(1, 2, 3))]
```


| sim|     nb1|     nb2|     nb3|maxj |
|---:|-------:|-------:|-------:|:----|
| 125| 1185274| 1659903| 1576115|2    |
| 207| 1154559| 1393306| 1635672|3    |
| 237| 1154464| 1298546| 1432624|3    |
| 293| 1216423| 1432172| 1612988|3    |
| 324| 1173978| 1494560| 1507685|3    |
| 375| 1194971| 1495585| 1692310|3    |
| 516| 1176348| 1467132| 1521095|3    |
| 527| 1247409| 1731416| 1550940|2    |
| 843| 1191053| 1481843| 1710277|3    |
| 907| 1175678| 1544974| 1681721|3    |

```r
mce <- prop.table(table(ce.nb$maxj))
print(mce)
```

```
## 
##   1   2   3 
## 0.0 0.2 0.8
```

In this example, treatments 1, 2, and 3 have the highest net benefits a fraction 0, 0.2, and 0.8 of the time respectively. The `psa` function performs this same calculations for a range of values of $k$ and all `nsims` random draws of the simulation output. The output is a tidy `data.table` which facilitates plotting with `ggplot`.

MCE plot

```r
library("ggplot2")
library("scales")
theme_set(theme_bw())
ggplot(psa.dt$mce, aes(x = k, y = prob, col = factor(arm))) +
  geom_line() + facet_wrap(~grp) + xlab("Willingess to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Arm")
```

![Probability each arm is most cost-effective by subgroup](figs/mce_plot-1.png)

In group 1, Arm 2 provides the greatest net benefits with high probability for almost all reasonable values of k. In group 2, the results are less certain, although Arm 3 provides the greatest net benefits with a higher probability than Arm 2. 

## Value of perfect information
One draw back of the previous measure is that it ignores the magnitude of cost or QALY gains. A measure which combines the probability of being most effective with the magnitude of the expected net benefit is the expected value of perfect information (EVPI). Intuitively, the EVPI provides an estimate of the amount that a decision maker would be willing to pay to collect additional data and completely eliminate uncertainty. Mathematically, the EVPI is defined as the difference between the maximum expected net benefit given perfect information and the maximum expected net benefit given current information. In other words, we calculate the net benefit for the optimal treatment arm for each random draw of the parameters and compare that to the net benefit for the treatment arm that is optimal when averaging across all parameters, 

$$
\begin{aligned}
EVPI &= \int_\theta NB_s(j_s^{*}, \theta^s) - NB_s(j_s^{*}, \theta). \\
\end{aligned}
$$

To illustrate consider the same random sample of 10 draws from our simulation output used above.


```r
armmax.g2 <- which.max(enb[[3]])
ce.nb[, nbpi := apply(ce.nb[, .(nb1, nb2, nb3)], 1, max)]
ce.nb[, nbci := ce.nb[[armmax.g2 + 1]]]
kable(ce.nb, digits = 0)
```



| sim|     nb1|     nb2|     nb3|maxj |    nbpi|    nbci|
|---:|-------:|-------:|-------:|:----|-------:|-------:|
| 125| 1185274| 1659903| 1576115|2    | 1659903| 1576115|
| 207| 1154559| 1393306| 1635672|3    | 1635672| 1635672|
| 237| 1154464| 1298546| 1432624|3    | 1432624| 1432624|
| 293| 1216423| 1432172| 1612988|3    | 1612988| 1612988|
| 324| 1173978| 1494560| 1507685|3    | 1507685| 1507685|
| 375| 1194971| 1495585| 1692310|3    | 1692310| 1692310|
| 516| 1176348| 1467132| 1521095|3    | 1521095| 1521095|
| 527| 1247409| 1731416| 1550940|2    | 1731416| 1550940|
| 843| 1191053| 1481843| 1710277|3    | 1710277| 1710277|
| 907| 1175678| 1544974| 1681721|3    | 1681721| 1681721|

To calculate EVPI, we average net benefits given current information and net benefits given perfect information accross simulation draws. 


```r
enbpi <- mean(ce.nb$nbpi)
enbci <- mean(ce.nb$nbci)
print(enbpi)
```

```
## [1] 1618569
```

```r
print(enbci)
```

```
## [1] 1592143
```

```r
print(enbpi - enbci)
```

```
## [1] 26426.37
```

The `psa` function peforms this same calculation accross all simulation draws from the PSA and for a number of values of willingess to pay values $k$. A plot by group of the the EVPI for different values of $k$ is shown below. The kinks in the plot represent values of $k$ where the optimal arm changes.


```r
ggplot(psa.dt$evpi, aes(x = k, y = evpi)) +
  geom_line() + facet_wrap(~grp) + xlab("Willingess to pay") +
  ylab("Expected value of perfect information") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  scale_y_continuous(label = scales::dollar) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Arm")
```

![Expected value of perfect information by subgroup](figs/evpi_plot-1.png)

## Distribution of health and cost outcomes

```r
ce[, lys := qalys * 1.5]
cea.fun <- function(x) list(mean = mean(x), quant = quantile(x, c(.025, .975)))
psa.custom.dt <- psa(ce, k = seq(0, ktop, 500), sim = "sim", arm = "arm",
               grp = "grp", e = "qalys", c = "cost",
               custom_vars = c("cost", "lys", "qalys"), 
               custom_fun = cea.fun)
```


```r
psa.custom.dt$summary
```

```
##      arm     grp qalys_mean qalys_lower qalys_upper    cost_mean
## 1: Arm 1 Group 1   7.996257    7.588626    8.430087     7.395474
## 2: Arm 2 Group 1  10.003105    8.378064   11.560632 60316.489409
## 3: Arm 3 Group 1   8.496872    7.283156    9.678301 60508.068527
## 4: Arm 1 Group 2   7.988661    7.582103    8.396929     7.452108
## 5: Arm 2 Group 2  10.481121    8.845461   12.004806 60436.512864
## 6: Arm 3 Group 2  11.010114    9.846793   12.163617 60078.680655
##      cost_lower   cost_upper
## 1:     6.076042     8.908801
## 2: 45348.413567 80314.637656
## 3: 45401.087646 80398.671574
## 4:     6.138874     9.117188
## 5: 45257.998841 81490.980326
## 6: 44489.642752 78849.935696
```


```r
psa.custom.dt$custom.table
```

```
##      arm     grp    cost.mean cost.quant.2.5% cost.quant.97.5% lys.mean
## 1: Arm 1 Group 1     7.395474        6.076042         8.908801 11.99439
## 2: Arm 2 Group 1 60316.489409    45348.413567     80314.637656 15.00466
## 3: Arm 3 Group 1 60508.068527    45401.087646     80398.671574 12.74531
## 4: Arm 1 Group 2     7.452108        6.138874         9.117188 11.98299
## 5: Arm 2 Group 2 60436.512864    45257.998841     81490.980326 15.72168
## 6: Arm 3 Group 2 60078.680655    44489.642752     78849.935696 16.51517
##    lys.quant.2.5% lys.quant.97.5% qalys.mean qalys.quant.2.5%
## 1:       11.38294        12.64513   7.996257         7.588626
## 2:       12.56710        17.34095  10.003105         8.378064
## 3:       10.92473        14.51745   8.496872         7.283156
## 4:       11.37315        12.59539   7.988661         7.582103
## 5:       13.26819        18.00721  10.481121         8.845461
## 6:       14.77019        18.24542  11.010114         9.846793
##    qalys.quant.97.5%
## 1:          8.430087
## 2:         11.560632
## 3:          9.678301
## 4:          8.396929
## 5:         12.004806
## 6:         12.163617
```

## Cost-effectiveness plane
The cost-effectivenss plane plots the incremental effectivenss of a treatment arm (relative to a comparator) against the incremental cost of the treatment arm. The plot is useful because it demonstrates both the uncertainty and the magnitutde of the estimates. Each point on the plot is from a particular random draw from the PSA.  

Data for plotting a cost-effectiveness plane comes from the `delta` output generated from the `psa_pw` function, which, for each sampled parameter set, estimtes incremental differences in costs, clinical effectivenss, and any other variables specifies in `custom_vars` relative to the comparator. The dotted line in the plot is the willingess to pay line, with slope equal to the value of $k$. For a given $k$, points below the line are cost-effective while those above are not. 


```r
head(psa.pw.dt$delta)
```

```
##    sim   arm     grp   iqalys    icost
## 1:   1 Arm 2 Group 1 2.134922 42705.08
## 2:   2 Arm 2 Group 1 2.022850 80448.94
## 3:   3 Arm 2 Group 1 3.012861 61797.37
## 4:   4 Arm 2 Group 1 1.996887 66735.06
## 5:   5 Arm 2 Group 1 2.907354 52049.11
## 6:   6 Arm 2 Group 1 3.031366 66954.19
```

```r
ylim <- max(psa.pw.dt$delta[, icost]) * 1.1
xlim <- ceiling(max(psa.pw.dt$delta[, iqalys]) * 1.1)
ggplot(psa.pw.dt$delta, aes(x = iqalys, y = icost, col = factor(arm))) + 
  geom_jitter(size = .5) + facet_wrap(~grp) + 
  xlab("Incremental QALYs") + ylab("Incremental cost") +
  scale_y_continuous(label = dollar, limits = c(-ylim, ylim)) +
  scale_x_continuous(limits = c(-xlim, xlim), breaks = seq(-6, 6, 2)) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Arm") +
  geom_abline(slope = 150000, linetype = "dashed") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
```

![Cost-effectiveness plane by subgroup](figs/ceplane_plot-1.png)

## Cost-effectiveness acceptability curve (CEAC)
The cost-effectivenss acceptability curve (CEAC) is similar to the MCE plot. The difference is that the CEAC compares each arm to a single comparator wheras the MCE plot considers all arms simultaneously. Output to produce a CEAC is generated from the `psa_pw` function.

The plot shows that, in subgroup 1, Arm 2 has large net benefits than Arm 1 with very high probability for reasonable values of $k$. Arm 3 also has higher net benefits than Arm 1 with probability over 1/2 for values of $k$ larger than 119,500. In group 2, both Arm 2 and Arm 3 have higher net benefits than Arm 1 for almost all values of $k$, although this probability is larger form Arm 2 than Arm 3 when $k$ is smaller.


```r
ggplot(psa.pw.dt$ceac, aes(x = k, y = prob, col = factor(arm))) +
  geom_line() + facet_wrap(~grp) + xlab("Willingess to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Arm")
```

![Cost-effectiveness acceptability curve by subgroup](figs/ceac_plot-1.png)

## Credible intervals for incremental costs and effectiveness
`psa_pw` contains summary output with expected costs and expected efficacy as well as their respective 95% credible intervals. The table also contains the ICER, which is equal to expected costs divided by the measure of expected efficacy.  


```r
print(psa.pw.dt$summary)
```

```
##      arm     grp iqalys_mean iqalys_lower iqalys_upper icost_mean
## 1: Arm 2 Group 1   2.0068483    0.2374925     3.632925   60309.09
## 2: Arm 3 Group 1   0.5006157   -0.8198998     1.780500   60500.67
## 3: Arm 2 Group 2   2.4924598    0.8422443     4.021185   60429.06
## 4: Arm 3 Group 2   3.0214535    1.7988992     4.237330   60071.23
##    icost_lower icost_upper      icer
## 1:    45341.18    80307.81  30051.65
## 2:    45392.98    80391.27 120852.53
## 3:    45250.43    81484.04  24244.75
## 4:    44481.62    78843.60  19881.57
```

If the user would like to examine the distribution of outcomes other than those specified this summary table, then they can also generate a custom table of summary output. The custom table can contain any quantities of interest (QOIs) as long as they are specified in addition to the `sim`, `arm`, `e`, and `c` columns in the posterior distribution data table. The default is to estimate means, the 2.5% quantile, and the 97.5% quantile for each variable, but any custom function can used. Below, we create a hypothetical variable for life-years and create table summarizing our estimates of costs, QALYs and life-years. A custom function, identical to the default option, is entered into the function for illustrative purposes.


```r
ce[, lys := qalys * 1.5]
cea.fun <- function(x) list(mean = mean(x), quant = quantile(x, c(.025, .975)))
psa.custom.dt <- psa(ce, k = seq(0, ktop, 500), sim = "sim", arm = "arm",
               grp = "grp", e = "qalys", c = "cost",
               custom_vars = c("cost", "lys", "qalys"), 
               custom_fun = cea.fun)
```


```r
psa.custom.dt$custom.table
```

```
##      arm     grp    cost.mean cost.quant.2.5% cost.quant.97.5% lys.mean
## 1: Arm 1 Group 1     7.395474        6.076042         8.908801 11.99439
## 2: Arm 2 Group 1 60316.489409    45348.413567     80314.637656 15.00466
## 3: Arm 3 Group 1 60508.068527    45401.087646     80398.671574 12.74531
## 4: Arm 1 Group 2     7.452108        6.138874         9.117188 11.98299
## 5: Arm 2 Group 2 60436.512864    45257.998841     81490.980326 15.72168
## 6: Arm 3 Group 2 60078.680655    44489.642752     78849.935696 16.51517
##    lys.quant.2.5% lys.quant.97.5% qalys.mean qalys.quant.2.5%
## 1:       11.38294        12.64513   7.996257         7.588626
## 2:       12.56710        17.34095  10.003105         8.378064
## 3:       10.92473        14.51745   8.496872         7.283156
## 4:       11.37315        12.59539   7.988661         7.582103
## 5:       13.26819        18.00721  10.481121         8.845461
## 6:       14.77019        18.24542  11.010114         9.846793
##    qalys.quant.97.5%
## 1:          8.430087
## 2:         11.560632
## 3:          9.678301
## 4:          8.396929
## 5:         12.004806
## 6:         12.163617
```

# Value of personalized care
The previous analyses allow net benefits and optimal treatment decisions to vary by subgroup. In contrast, most CEAs estimate the treatment, $j^{*}$, that is optimal when averaging net benefits over the entire population. In particular, if the population is broken up into $G$ distinct subgroups, the optimal treatment is given by,

$$
\begin{aligned}
j^{*} = \text{argmax}_j \sum_{g=1}^{G} w_g E_{\theta}\left[NB_g(j,\theta)\right],
\end{aligned}
$$

@basu2007value have shown that selecting subgroup specific treatments increases expected net benefits relative to this one-size fits all approach. They refer to additional net benefit as the expected value of individualized care (EPIC), which can be computed in using the subgroup approach illustrated here as,

$$
\begin{aligned}
\sum_{g=1}^G w_g E_{\theta}\left[NB_g(j^{*}_s,\theta)\right] - \sum_{g=1}^G w_g  E_{\theta}\left[NB_g(j^{*},\theta)\right],
\end{aligned}
$$

where $w_g \in (0, 1)$ is a weight denoting that proportion of the population represented by subgroup $g$ and $\sum_{g=1}^{G} w_g = 1$. We can use `hesim` to demonstrate the value of individualized care.


```r
w.dt <- data.table(grp = paste0("Group ", seq(1, 2)), w = c(0.25, .75))
evpi <- psa.dt$evpi
evpi <- merge(evpi, w.dt, by = "grp")
totevpi <- evpi[,lapply(.SD, weighted.mean, w = w),
                by = "k", .SDcols = c("evpi")]
ggplot(totevpi, aes(x = k, y = evpi)) +
  geom_line() + xlab("Willingess to pay") +
  ylab("Total EVPI") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  scale_y_continuous(label = scales::dollar) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Arm")
```

![Total expected value of perfect information](figs/totevpi-1.png)


```r
ce <- merge(ce, w.dt, by = "grp")
totenb <- ce[, .(totenb = weighted.mean(nb, w = w)), by = c("arm")]
```


```r
ptenb.grp.max <- apply(as.matrix(enb[, -1]), 2, max)
ptenb.max <- sum(ptenb.grp.max * w.dt$w)
tenb.max <- max(totenb$totenb)
tnb <- c(ptenb.max, tenb.max)
names(tnb) <- c("Personalized total TENB", "One-size fits all TENB")
```


```r
evic <- tnb[1] - tnb[2]
names(evic) <- "EVIC"
print(evic)
```

```
##     EVIC 
## 56531.62
```

```r
print(evic/150000)
```

```
##      EVIC 
## 0.3768774
```

Our estimate of the EVIC is \$56,531, or in terms of net health benefits, 0.377 QALYs. 

# References

