---
title: "Estimating the overall ATT in a staggered DID design using a bayesian model"
output: 
  html_document:
    code_folding: hide
    keep_md: yes
bibliography: bayes_did_agg_bib.bib
date: "2024-01-23"
---

## Disclaimer

I am currently a Ph.D student in sociology and a very "applied" kind of Ph.D student at that. As such, any claims I make here about statistics, econometrics, and bayesian statistics in particular, are subject to scrutiny, to say the least. In fact, the whole reason for making this short document here is that I am unsure whether or not it makes sense at all to do the kind of aggregation I propose here. So I am simply making this to, hopefully, show it to people who are smarter than me so that we might find out if doing this kind of aggregation is even the slightest bit sensible. In other words, this document is not so much me trying to say, "Look, I did this thing, I think it's kinda smart" as it is me saying, "Look, I did this thing, I think I might be wrong, can you tell me how to do it right?".

# Intro

In this short document, I am going to try and see what happens when we aggregate the group-time ATTs, $ATT(g,t)$, that result from the Callaway & Sant'Anna (CSA) estimator[@callaway2021difference] for Differences-in-differences(DID), when treatment timing is staggered. I will compare this aggregation to the aggregation that can be done via the `did` package[@didpack], which implements the CSA estimator.

In this comparison, I will take a page from the book on bayesian meta-analysis. While meta-analysis is concerned with aggregating estimates of treatment effects from different independent studies, the procedure, at its absolute most basic, is simply a weighted mean of effect size estimates. The weights are the inverse of the estimated variance of the estimated effects, e.g., $1/v_i$. In that sense, the procedure is quite similar to the aggregation procedure used in the `did` package, given that the overall treatment effect that results from `did::aggte(mod, type = "group")` is also a weighted mean. Here, the overall ATT is a weighted mean of the ATTs of each treatment cohort, where the weights are the proportions of observations that belong to each treatment cohort, i.e., the size of the cohorts. Given that sample size and variance are closely related, I figured that using inverse-variance weights might yield results that are quite comparable.

## Why bayesian?

However, settling for "quite comparable" seems like a bad argument for not using methods that have now become standard repertoire in doing DID with staggered treatment timing. This is where the bayesian part comes in. The thing is that I like the kind of inference that can be drawn from bayesian models(e.g., being able to say, "There is a 95% chance that the true treatment effect is within this credible interval" rather than trying to figure out what is both a valid and useful interpretation of a frequentist confidence interval). The problem is that I have only begun dipping my toes into bayesian inference and modeling, so I cannot figure out how one would "port" the CSA estimator to a fully bayesian approach. Instead, I thought I might opt for a "shortcut" by simply aggregating the ATTs of each treatment cohort, which I get from the CSA estimator, using a bayesian approach. This way, I could get, e.g., a credible interval around the overall ATT and get an estimate of the overall ATT that is, hopefully, "quite comparable" to the estimated overall ATT from the CSA estimator. 
Another possible advantage of using a bayesian approach is that we can specify a prior using either theoretical knowledge or results from previous studies to inform our belief, and subsequently our model, of treatment effects. Secondly, using an informative prior may increase the precision of our results[@banner2020use]. Here, I would argue that we are likely to have more informative priors on the aggregated effect(or perhaps the immediate effect) of a treatment in a staggered DID design rather than having good and informative priors on the full distributions of $ATT(g,t)$.

<style type="text/css">
.table caption {
    font-weight: bold;
}
</style>





# Example using simulated data

In the first example i'll use simulated data to see if the two methods of aggregation are somewhat comparable

The procedure for simulating data is snatched directly from the DID package vignette, https://cran.r-project.org/web/packages/did/vignettes/TWFE.html


``` r
make_data3 <- function(nobs = 1000, 
                       nstates = 40,
                       effect_1=2,
                       effect_2=1,
                       effect_3=3) {
  
  # unit fixed effects (unobservd heterogeneity)
  unit <- tibble(
    unit = 1:nobs,
    # generate state
    state = sample(1:nstates, nobs, replace = TRUE),
    unit_fe = rnorm(nobs, state/5, 1),
    # generate instantaneous treatment effect
    #mu = rnorm(nobs, true_mu, 0.2)
    mu = true_mu
  )
  
  # year fixed effects (first part)
  year <- tibble(
    year = 1980:2010,
    year_fe = rnorm(length(year), 0, 1)
  )
  
  # Put the states into treatment groups
  treat_taus <- tibble(
    # sample the states randomly
    state = sample(1:nstates, nstates, replace = FALSE),
    # place the randomly sampled states into 1\{t \ge g \}G_g
    cohort_year = sort(rep(c(1986, 1992, 1998, 2004), 10))
  )
  
  # make main dataset
  # full interaction of unit X year 
  expand_grid(unit = 1:nobs, year = 1980:2010) %>% 
    left_join(., unit) %>% 
    left_join(., year) %>% 
    left_join(., treat_taus) %>% 
    # make error term and get treatment indicators and treatment effects
    # Also get cohort specific trends (modify time FE)
    mutate(error = rnorm(nobs*31, 0, 1),
           treat = ifelse((year >= cohort_year)* (cohort_year != 2004), 1, 0),
           mu = ifelse(cohort_year==1992, effect_1, ifelse(cohort_year==1998, effect_2, effect_3)),
           tau = ifelse(treat == 1, mu, 0),
           year_fe = year_fe + 0.1*(year - cohort_year)
    ) %>% 
    # calculate cumulative treatment effects
    group_by(unit) %>% 
    mutate(tau_cum = cumsum(tau)) %>% 
    ungroup() %>% 
    # calculate the dep variable
    mutate(dep_var = (2010 - cohort_year) + unit_fe + year_fe + tau_cum + error) %>%
    # Relabel 2004 cohort as never-treated
    mutate(cohort_year = ifelse(cohort_year == 2004, Inf, cohort_year))
  
}
```

## CSA on simulated data

Next we'll generate the data and run the CSA estimator to get an estimated overall ATT


``` r
#generate data
iseed  = 20201221
nrep <- 20  
true_mu <- 1
set.seed(iseed)


data <- make_data3(effect_1 = 2, effect_2 = 1, effect_3 = 3)

#data <- make_data3(effect_1 = .02, effect_2 = .001, effect_3 = .03, nobs=100)


#recode never-treated to 0
data$cohort_year[data$cohort_year==Inf] <- 0


#run CSA estimator
mod <- did::att_gt(yname = "dep_var", 
                   tname = "year",
                   idname = "unit",
                   gname = "cohort_year",
                   control_group= "notyettreated",
                   bstrap = TRUE,
                   data = data,
                   print_details = FALSE)

#Do aggregation at the grou/cohort level
agg_CSA <- did::aggte(mod, type = "group")
```


## Aggregated effect via bayesian model on simulated data


Much of the code below is snatched directly from https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/bayesian-ma.html (which BTW is an excellent resource for doing meta-analysis in R in general, 10/10 would recommend!)


``` r
#we'll go with the default prior from brms for now
#priors <- c(prior(normal(20,10), class = Intercept))

#make the data for the bayesian aggregation, which is basically the output from the CSA aggregation
df_bayes <- data.frame(TE=agg_CSA$att.egt, seTE=agg_CSA$se.egt, unit=c("1986","1992","1998"))

#run the model via the brms package
m.brm <- brm(TE|se(seTE) ~ 1,
             data = df_bayes,
             iter = 4000, silent=2, refresh=0)



# m.brm <- brm(TE|se(seTE) ~ 1,
#              data = df_bayes,
#              iter = 4000, silent=2, refresh=0,
#              prior = priors)


sum_brm <- summary(m.brm)
```


Lets compare the results from CSA and the bayesian aggregation


``` r
#we make a table with results from both CSA and the bayesian aggregation.
knitr::kable(
  data.frame(method=c("CSA", "bayesian"), 
             ATT=c(agg_CSA$overall.att, sum_brm$fixed$Estimate), 
             SE=c(agg_CSA$overall.se, sum_brm$fixed$Est.Error),
             CI=c(paste(round(agg_CSA$overall.att-1.96*agg_CSA$overall.se,2), 
                        round(agg_CSA$overall.att+1.96*agg_CSA$overall.se,2), collapse=";"),
                  paste(round(sum_brm$fixed[3],2), round(sum_brm$fixed[4],2), collapse=";"))),
  caption="<span style='font-size:200%'>Aggregated ATT of simulated data</span>"
  
  
  )
```



Table: <span style='font-size:200%'>Aggregated ATT of simulated data</span>

|method   |      ATT|        SE|CI          |
|:--------|--------:|---------:|:-----------|
|CSA      | 22.46637| 0.0576607|22.35 22.58 |
|bayesian | 22.57214| 0.0507905|22.47 22.67 |


# Using real data

In the next example i'll use a real dataset which is the "castle" dataset from @cheng2013does. The dataset contains information on homicide rates in US states before and after the implementation of the castle-doctrine, the so-called "stand you ground" laws(https://rdrr.io/cran/did2s/man/castle.html). The data and code for fitting the CSA estimator is downloaded from a Github repository maintained by Scott Cunningham, so much thanks to Scott! Alternately the data is also included in the `did2s` package[@did2s].


``` r
#taken straight from https://github.com/scunning1975/causal-inference-class/blob/master/Master_do_files/castle%20_cs.R

library(readstata13)


castle <- data.frame(read.dta13('https://github.com/scunning1975/mixtape/raw/master/castle.dta'))
castle$effyear[is.na(castle$effyear)] <- 0 # untreated units have effective year of 0

# Estimating the effect on log(homicide)
atts <- att_gt(yname = "l_homicide", # LHS variable
               tname = "year", # time variable
               idname = "sid", # id variable
               gname = "effyear", # first treatment period variable
               data = castle, # data
               xformla = NULL, # no covariates
               #xformla = ~ l_police, # with covariates
               est_method = "dr", # "dr" is doubly robust. "ipw" is inverse probability weighting. "reg" is regression
               control_group = "nevertreated", # set the comparison group which is either "nevertreated" or "notyettreated" 
               bstrap = TRUE, # if TRUE compute bootstrapped SE
               biters = 1000, # number of bootstrap iterations
               print_details = FALSE, # if TRUE, print detailed results
               clustervars = "sid", # cluster level
               panel = TRUE) # whether the data is panel or repeated cross-sectional

# Aggregate ATT
agg_effects <- aggte(atts, type = "group")
```

Next we do the bayesian aggregation, same procedure as for the simulated data


``` r
#make the data for bayesian aggregation
df_bayes <- data.frame(TE=agg_effects$att.egt, seTE=agg_effects$se.egt, unit=agg_effects$DIDparams$glist)

#we still roll with the default prior
#priors <- c(prior(normal(10,10), class = Intercept))



m.brm_castle <- brm(TE|se(seTE) ~ 1,
             data = df_bayes,
             iter = 4000, silent=2, refresh=0)



sum_brm_castle <- summary(m.brm_castle)
```

And we compare results


``` r
knitr::kable(
  data.frame(method=c("CSA", "bayesian"), 
             ATT=c(agg_effects$overall.att, sum_brm_castle$fixed$Estimate), 
             SE=c(agg_effects$overall.se, sum_brm_castle$fixed$Est.Error),
             CI=c(paste(round(agg_effects$overall.att-1.96*agg_effects$overall.se,2), 
                        round(agg_effects$overall.att+1.96*agg_effects$overall.se,2), collapse=";"),
                  paste(round(sum_brm_castle$fixed[3],2), round(sum_brm_castle$fixed[4],2), collapse=";"))),
  caption="<span style='font-size:200%'>Aggregated ATT of Castle-doctrine on homicide rates</span>"
  )
```



Table: <span style='font-size:200%'>Aggregated ATT of Castle-doctrine on homicide rates</span>

|method   |       ATT|        SE|CI        |
|:--------|---------:|---------:|:---------|
|CSA      | 0.1084475| 0.0374177|0.04 0.18 |
|bayesian | 0.0795568| 0.0199254|0.04 0.12 |


# more real data: Effect of minimum wage on youth employment

In this section we'll do the whole thing one last time on a dataset that is included in the `did` package concerning the effect of minimum wage implementation on youth employment


``` r
#this piece of code is taken straigt from https://cran.r-project.org/web/packages/did/vignettes/did-basics.html

#load data

data(mpdta)


#fit model

mw.attgt <- att_gt(yname = "lemp",
                   gname = "first.treat",
                   idname = "countyreal",
                   tname = "year",
                   xformla = ~1,
                   data = mpdta,
                   )



agg_mw <- aggte(mw.attgt, type="group")
```

and do the bayesian aggregation one last time


``` r
#make the data for bayesian aggregation
df_bayes <- data.frame(TE=agg_mw$att.egt, seTE=agg_mw$se.egt, unit=agg_mw$DIDparams$glist)

#we still roll with the default prior
#priors <- c(prior(normal(10,10), class = Intercept))



m.brm_mw <- brm(TE|se(seTE) ~ 1,
             data = df_bayes,
             iter = 4000, silent=2, refresh=0)



sum_brm_mw <- summary(m.brm_mw)
```


and do one last comparison


``` r
knitr::kable(
  data.frame(method=c("CSA", "bayesian"), 
             ATT=c(agg_mw$overall.att, sum_brm_mw$fixed$Estimate), 
             SE=c(agg_mw$overall.se, sum_brm_mw$fixed$Est.Error),
             CI=c(paste(round(agg_mw$overall.att-1.96*agg_mw$overall.se,2), 
                        round(agg_mw$overall.att+1.96*agg_mw$overall.se,2), collapse=";"),
                  paste(round(sum_brm_mw $fixed[3],2), round(sum_brm_mw$fixed[4],2), collapse=";"))),
  caption="<span style='font-size:200%'>Aggregated ATT of minimum wage on youth unemployment</span>"
  )
```



Table: <span style='font-size:200%'>Aggregated ATT of minimum wage on youth unemployment</span>

|method   |        ATT|        SE|CI          |
|:--------|----------:|---------:|:-----------|
|CSA      | -0.0310183| 0.0122481|-0.06 -0.01 |
|bayesian | -0.0322365| 0.0111533|-0.05 -0.01 |


# Why bayesian? Probing the posterior distribution

So far, we haven't gotten much more out of the bayesian aggregation, compared to the CSA aggregation, other than a credible interval instead of a confidence interval(which is already something, in my opinion). However, now that we have aggregated the ATTs using a bayesian model, we can now draw from the posterior distribution of the aggregated ATT and get a full probability distribution of the ATT. This allows us to use other kinds of inference that I think wouldn't be so easy to achieve using frequentist inference. 

For instance, we can take a look at the posterior distribution of the effects of implementing the Castle-doctrine on homicide rates. 

In the figure below, we see the distribution of the estimated ATTs from the bayesian model of the impact of the Castle-doctrine, which are sampled from the posterior distribution of the model. In other words, the plot below shows a _probability distribution_ of the estimated ATT. Additionally, the black error-bar shows the credible interval from the bayesian model, while the red error bar represents the confidence interval from the CSA model. The vertical dashed lines represent the point estimates of the ATTs, with the black and red dashed lines representing the bayesian and CSA model estimates, respectively.




``` r
post.samples <- posterior_samples(m.brm_castle)

library(ggplot2)

p <- ggplot(post.samples, aes(b_Intercept))
p <- p+geom_density()
p <- p+geom_errorbarh(aes(xmin=quantile(post.samples$b_Intercept, probs = .025), 
                      xmax=quantile(post.samples$b_Intercept, probs = .975), y=.5))
p <- p+geom_errorbarh(aes(xmin=agg_effects$overall.att-1.96*agg_effects$overall.se, 
                      xmax=agg_effects$overall.att+1.96*agg_effects$overall.se, y=1.5, color="csa"))
p <- p+geom_vline(xintercept=sum_brm_castle$fixed$Estimate, linetype="dashed")
p <- p+geom_vline(aes(xintercept=agg_effects$overall.att, color="csa"), linetype="dashed")
p <- p+scale_x_continuous(breaks=seq(0,.2,.01))
p <- p+scale_color_manual(values=c("csa"="red"))
p <- p+theme_minimal()
p <- p+labs(x="Aggregated ATT of Castle-doctrine on homicide rates", y="Density")



#Empirical cummulative distribution function
f <- ecdf(post.samples$b_Intercept)
```


``` r
p
```

![]([plot](https://github.com/rasmusklokker/staggered-did-bayesian-aggregate/blob/master/staggered-DID-bayesian-aggregate-document_files/figure-html/posterior%20plot-1.png?raw=true "title"))




However, what if we had some kind of minimum effect of interest in mind, that we wanted to assess? For instance, what if we considered an increase in homicides of 5% a practically important effect size? Well, we can easily assess the probability of the ATT being 0.05 or higher, as well. In this case, the probability that the Castle-doctrine increased homicide rates of 5% or higher is 93%. What if we thought a 10% increase was the limit of practical importance? The probability that the policy resulted in an increase of homicide rates of 10%, or higher, is then 15%.

As such, I think that using a bayesian approach makes it easier to assess the uncertainty associated with an estimated parameter and to evaluate how the uncertainty relates to research questions that are more informed by either practical or scholarly information than the question "How likely is it that the effect is different from zero?". While I am a big fan of choosing confidence intervals over p-values to quantify the uncertainty associated with a parameter estimate, I also think being able to "query" the posterior distribution has some advantages over confidence intervals. For instance, if one had obtained a confidence interval on the estimated ATT of the Castle-doctrine that ranged from -0.02 to 0.25, we wouldn't be able to reject that Castle-doctrine could have reduced homicide rates by 2% or that it could have increased homicide rates by as much as 25%. One would be tempted to say that "most" of the confidence interval concerns an increase in homicide rates, e.g., that most of the evidence points towards an increase in homicide rates. However, it is my impression that confidence intervals do not afford that kind of inference. In such situations, I have usually commented that results were "inconclusive due to high uncertainty," which is an anti-climactic way to write up your results. In contrast, had we obtained a credible interval with the same bounds, we could have assessed the probability that, e.g., the change in homicide rates was positive, which would probably have been very likely. In such a situation, which I have found myself in quite a few times, I think that a bayesian approach can provide a more informative quantification of uncertainty, which means that the conclusion of all our hard work and spent taxpayer money doesn't necessarily have to be 🤷 


## Some caveats

I know that some of the use cases I present here are also fully available within a frequentist framework. For instance, a one-sided hypothesis test could have been used to test the hypothesis that the Castle-doctrine increased the homicide rate by 5% or more. However, this still doesn't get around the issues with interpreting p-values, one-sided hypothesis or not. 

Lastly, I'd like to mention that I have not addressed the issue of dependent effect sizes in this application. In this application, the issue arises since the cohorts are compared against the same control group, the never-treated group. This creates a dependence in the errors of the effect sizes, which means that we can't treat the ATTs from each treatment cohort as independent, which is what I have done here[@cheung2019guide]. However, I think that addressing this is relatively easy. My proposed solution would be to take a page from the literature on robust variance estimation in meta-analysis, where one constructs a variance-covariance matrix based on the assumed correlation between effect sizes[@pustejovsky2022meta]. This, I think, can be done relatively easy via the `fcor` function in `brms`[@brms]. However, there is likely also a more sophisticated and more bayesian, solution where the structure of errors is modeled directly, priors and all.
However, if we wanted a bayesian approach that could fully emulate the CSA approach, things might get a bit more hairy. The CSA estimator allows you also to use the "not-yet-treated" groups as well as the "never-treated" group. This potentially buys you much more power to estimate the ATTs but also complicates the question of dependence on errors. For instance, if we have four treated cohorts in 1996, 1998, and 2001, and the last in 2005, we could compute the ATT for the 1996 cohort by comparing it to the never-treated group and the units treated in 1998, 2001, and 2005. For the group treated in 1998, the control group similarly consists of the never treated group and those treated in 2001 and 2005. There is clearly an overlap between the control group used for computing the ATTs for the 1996 and 1998 cohorts. However, the control groups are neither identical nor independent, making the case more complex. One might note that the overlap in control groups follows a pattern. For instance, the 1996 cohort has more overlap in control groups with the 1998 cohort, less with the 2001 cohort, and for the 2005 cohort, only the never-treated group is in common. Such a pattern can be modeled, maybe even credibly so, and so we are, assumingly, back to this being a skill-issue(i.e., the lack of my skills in bayesian modeling) rather than some insurmountable "3-body-problem"-esque challenge.

# Bibliography
