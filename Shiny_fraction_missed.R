
# - - - - - - - - - - - - - - - - - - - - - - -
# Simulate the fraction of the population caught or missed in a growing epidemic.
# - - - - - - - - - - - - - - - - - - - - - - -
rm(list = ls())
source("Code_Screening_model.R")
source("Code_distribution_functions.R")
library(ggplot2)
library(grid)
library(gridExtra)
library(tidyverse)
pathogen=c("nCoV")
pathtablab=c("2019-nCoV")
par(mfrow = c(length(pathtablab), 1))

ftime=24          #Travel time (hours)
flatA=0           #Growing (0) or flat (1) epidemic
relT=0            # Overall (0) or arrival (1) screening
sd=0.7            #fever screening efficacy at departure
sa=sd             #fever screening efficacy at arrival
del.d = ftime/24  #convert fligth time from hours to days
nboot = 50        # n bootstrap samples
popn = 100        # population size of infected travelers
meanIncubate = 5.5
shapeIncubate = NULL
scaleIncubate = 2.73


## Get incubation period distribution
if(length(meanIncubate)>0){         ## If mean incubation period specified, calculate shape
  shapeIncubate = meanIncubate/scaleIncubate
}else if(length(shapeIncubate)==0){ ## Else, use explicit shape and scale inputs
  stop('Error - Must specify meanIncubate or shape and scale!')
}
# Set the incubation function
incubation.d<-(function(x){dgamma(x, shape = shapeIncubate, scale = scaleIncubate)})


## \\\\\\\\\\\\\\\\\\ SIMULATE ///////////////////
one_bootstrap = function(){
## For each individual in the population, draw times since exposre from the appropriate distribution
infAge=sapply(c(1:popn),function(x){exposure.distn(runif(1, 0, 1),pathogen,flat=flatA)})

## Given times since exposure, screen using fever, risk or both
screenWrapper = function(x, ff, gg){
  screen.passengers(x, del.d, ff, gg, sd, sa, rd=.25 , ra=.25, phi.d, incubation.d,pathogen,0,split1 = 2)}

## Fever and risk screening
outcomesBoth=sapply(infAge,function(x){
  f0=rbinom(1,fg.distn("f",pathogen)[2],fg.distn("f",pathogen)[1])/fg.distn("f",pathogen)[2]  ## Draw binomial probability of fever
  g0=rbinom(1,fg.distn("g",pathogen)[2],fg.distn("g",pathogen)[1])/fg.distn("g",pathogen)[2]  ## Draw binomial probability of known risk  
  screenWrapper(x, f0, g0)})                                                                  ## Output individual probability missed
pCaught_both = colSums(outcomesBoth[1:4,])
caughtBoth = sapply(pCaught_both,function(x){ifelse(x<runif(1, 0, 1),0,1)})                   ## Draw whether individual was missed


## Fever only
outcomesFever=sapply(infAge,function(x){
  f0=rbinom(1,fg.distn("f",pathogen)[2],fg.distn("f",pathogen)[1])/fg.distn("f",pathogen)[2]  ## Draw binomial probability of fever
  g0=0  ## Draw binomial probability of known risk  
  screenWrapper(x, f0, g0)})     
pCaught_fever = colSums(outcomesFever[1:4,])                                                    ## Output individual probability missed
caughtFever = sapply(pCaught_fever,function(x){ifelse(x<runif(1, 0, 1),0,1)})                   ## Draw whether individual was missed


## Risk only
outcomesRisk=sapply(infAge,function(x){
  f0=0  ## Draw binomial probability of fever
  g0=rbinom(1,fg.distn("g",pathogen)[2],fg.distn("g",pathogen)[1])/fg.distn("g",pathogen)[2]  ## Draw binomial probability of known risk  
  screenWrapper(x, f0, g0)})                                                                  ## Output individual probability missed
pCaught_risk = colSums(outcomesRisk[1:4,])
caughtRisk = sapply(pCaught_risk,function(x){ifelse(x<runif(1, 0, 1),0,1)})                   ## Draw whether individual was missed

## Return frac.missed.both, frac.missed.fever, frac.missed.risk
# return(c(frac.missed.both = 1- sum(caughtBoth)/popn,
#          frac.missed.fever = 1- sum(caughtFever)/popn,
#          frac.missed.risk = 1- sum(caughtRisk)/popn))

return(list(outcomesBoth = outcomesBoth,
            outcomesFever = outcomesFever,
            outcomesRisk = outcomesRisk,
            caught = c(frac.missed.both = 1- sum(caughtBoth)/popn,
                       frac.missed.fever = 1- sum(caughtFever)/popn,
                       frac.missed.risk = 1- sum(caughtRisk)/popn)))
}
# ## Test
# one_bootstrap()

replicate(nboot, expr = one_bootstrap()) %>% t() %>%
  as.data.frame() -> bootList

## Plot fraction missed, with CIs
sapply(bootList$caught, function(xx) xx) %>% t() %>% as.data.frame() %>%
  pivot_longer(1:3, names_to = c('screen_type'), names_pattern = 'frac.missed.(\\w+)', values_to = 'frac.missed') %>%
  group_by(screen_type) %>%
  summarise(med = median(frac.missed),
            lower = quantile(frac.missed, probs = .025),
            upper = quantile(frac.missed, probs = .975)) %>%
  ggplot()+
  geom_point(aes(x = screen_type, y = med, color = screen_type), size = 3)+
  geom_segment(aes(x = screen_type, xend = screen_type, y = lower, yend = upper, color = screen_type))+
  theme_bw()+
  ylim(c(0,1))+
  xlab('Screening type')+
  ylab('Fraction missed')+
  theme(legend.position = "none") -> fracMissed
fracMissed
ggsave('2020_nCov/Fraction_missed.pdf', width = 4, height = 4, units = 'in')  


## Plot stacked bar showing the probability of being caught, or missed in different ways
# extract_median_probs = function(bootListElement){
#   sapply(bootListElement, function(xx){apply(xx, MARGIN = 1, median)}) -> run_medians
#   apply(run_medians, MARGIN = 1, median)
# }

extract_mean_probs = function(bootListElement){
  sapply(bootListElement, function(xx){xx %>% t() %>% colMeans}) %>%
    t() %>% colMeans()
}

bind_rows(
  bothMedians = extract_mean_probs(bootList$outcomesBoth),
  feverMedians = extract_mean_probs(bootList$outcomesFever),
  riskMedians = extract_mean_probs(bootList$outcomesRisk),
) %>%
  mutate(strategy = c('both', 'fever', 'risk')) -> meanOutcomes
names(meanOutcomes) = c('d.fever', 'd.risk', 'a.fever', 'a.risk', 'm.b', 'm.f', 'm.r', 'nd', 'strategy')



meanOutcomes %>%
  pivot_longer(cols = 1:8) %>%
  mutate(outcome = factor(name, levels = (c('d.fever', 'd.risk', 'a.fever', 'a.risk', 'm.b', 'm.f', 'm.r', 'nd')),
                          labels =(c('stopped: departure fever screen', "stopped: departure risk screen", "stoppd: arrival fever screen",
                                     "stopped: arrival risk screen", 'cleared: missed both', 'cleared: missed fever', 'cleared: missed risk', 'cleared: not detectable')))) %>%
  ggplot()+
  geom_bar(aes(x = strategy, y = value, fill = outcome), stat = 'identity')+
  scale_fill_manual(values = c('darkseagreen2', 'deepskyblue', 'seagreen4', 'royalblue3', 'bisque', 'darksalmon', 'brown3', 'firebrick1')) +
  theme_bw() -> stackedBars
stackedBars
ggsave('2020_nCov/outcome_stacked_bars.pdf', width = 6, height = 4, units = 'in')  
