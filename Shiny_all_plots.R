# Code for 3 x 3 grid
rm(list = ls())
# Source external functions and load packages
source("Code_Screening_model.R")
source("Code_distribution_functions.R")
library(ggplot2)
library(grid)
library(gridExtra)
library(tidyverse)
cols = c('darkseagreen2', 'deepskyblue', 'seagreen4', 'royalblue3', 'bisque', 'darksalmon', 'brown3', 'firebrick1')

# - - - - - - - - - - - - - - - - - - - - - - -
# Plot proportion of travellers identified vs time from exposure to departure (Figure 3)
# - - - - - - - - - - - - - - - - - - - - - - -
#Specifiticty of pathogens of interest. Here we include pathogen-specific distributions 
#   to run the script for H7N9
pathogen=c("nCoV")
pathtablab=c("2019-nCoV")
par(mfrow = c(length(pathtablab), 1))

#Specify efficacy parameters. These will be fed in to the external functions.
rd=0.25 #efficacy of departure questionnaire (proportion of travelers that report honestly)
ra=0.25 #efficacy of arrival questionnaire
sd=0.7 #efficacy of departure fever screening (based on fever detection sensitivity)
sa=0.7 #efficacy of arrival fever screening

#Growing (0) or flat (1) epidemic
flatA=0
# Overall (0) or arrival (1) screening
relT=0

if(flatA==0&relT==0){kk=1}
if(flatA==0&relT==1){kk=2}
if(flatA==1&relT==0){kk=3}
if(flatA==1&relT==1){kk=4}


## Set a grid of fever and incubation inputs
input_grid = expand.grid(ffs = c(.6, .75, .95),
                         meanIncs = c(4, 5.5, 7))



# -------------------------------
# \\\\\\\\\\\\\\\  WRITE FUNCTIONS AND WRAPPERS ///////////////////
# -------------------------------

# -------------------------------
## Write a function to input f, g, R0, mean onset->admit, mean incubation period, days in flight (fraction)
get_frac_caught = function(tSinceExposed, ff, gg, R0, meanToAdmit, meanIncubate = NULL, daysInFlight, shapeIncubate = NULL, scaleIncubate = 2.73){
  ## INPUTS
  ## f1 - probability of fever at onset
  ## g1 - probability aware of risk
  ## meanToAdmit - mean days from onset to admission (this is the detectable window) 
  ## meanIncubate - mean incubation period (exposure ot onset). If meanIncubate specified, assume a gamma distribution with scale = 2.73 and calculate appropriate shape. Alternatively, can specify shape and scale explicitly.
  ## OUTPUTS - vector. Fraction caught by: departure risk, departure fever, arrival risk, arrival fever, cleared.
  
  if(length(meanIncubate)>0){         ## If mean incubation period specified, calculate shape
    shapeIncubate = meanIncubate/scaleIncubate
  }else if(length(shapeIncubate)==0){ ## Else, use explicit shape and scale inputs
    stop('Error - Must specify meanIncubate or shape and scale!')
  }
  # Get probability that symptoms have developed.
  incubation.d<-(function(d)pgamma(d, shape = shapeIncubate, scale = scaleIncubate))
  ## Outputs
  screen.passengers(tSinceExposed, del.d = 1, ff, gg, sd, sa, rd, ra, phi.d, incubation.d, pathogen, relT, split1 = 2)
}
# ## Test function
# get_frac_caught(tSinceExposed = 2, ff = .7, gg = .2, R0 = 2, meanToAdmit = 5.5, meanIncubate = 4.5, daysInFlight = 1)
# -------------------------------

# -------------------------------
## Write a function that repeats get_frac_caught over a grid of times since exposure
get_frac_caught_over_time = function(ff, gg, R0, meanToAdmit, meanIncubate = NULL, daysInFlight = NULL, shapeIncubate = NULL, scaleIncubate = 2.73){
  arrive.times=seq(0,15,0.1)
  sapply(arrive.times, function(tt){get_frac_caught(tSinceExposed = tt, ff, gg, R0, meanToAdmit, meanIncubate = 4.5, daysInFlight = 1)}) %>% t() %>% as.data.frame() %>% mutate(
    days.since.exposed = arrive.times)
}
# Test
# get_frac_caught_over_time(ff = .7, gg = .2, R0 = 2, meanToAdmit = 5.5)
# -------------------------------




# -------------------------------
##  \\\\\\\\\ RUN ON A GRID OF FEVER PROBS AND INCUBATION TIMES /////////
# -------------------------------
apply(X = input_grid, MARGIN = 1, FUN = function(ii){get_frac_caught_over_time(ff = ii[1], gg = .25, R0 = 2, meanToAdmit = ii[2])}) %>%                                                                        ## This would bind outputs into a long data frame instead of a list
  bind_rows() %>% as.tbl() %>%
  mutate(fever = rep(input_grid$ffs, each = nrow(.)/nrow(input_grid)),
         meanIncubate = rep(input_grid$meanIncs, each = nrow(.)/nrow(input_grid))) -> gridOutputs

## Reformat for plotting
gridOutputs %>%
  ## Specify minim and maximum y values in each band
  mutate(dFeverMin = 0, dFeverMax = dFeverMin + caught.dpt.fever,
         dRiskMin = dFeverMax, dRiskMax = dRiskMin + caught.dpt.risk,
         aFeverMin = dRiskMax, aFeverMax = aFeverMin + caught.arv.fever,
         aRiskMin = aFeverMax, aRiskMax = aRiskMin +caught.arv.risk,
         mbMin = aRiskMax, mbMax = mbMin+missed.both,
         mfMin = mbMax, mfMax = mfMin+missed.fever.only,
         mrMin = mfMax, mrMax = mrMin+missed.risk.only,
         ndMin = mrMax, ndMax = 1) %>%
  select(days.since.exposed, contains('Min'), contains('Max'), fever, meanIncubate) %>%
  ## Pivot to long data frame
  pivot_longer(cols = dFeverMin:ndMax, names_to = c('outcome', 'minOrMax'), names_pattern = '(\\w+)(M\\w\\w)', values_to = 'yy') -> temp
## Use full join to create columns for time, ymin, ymax, and band type
full_join(filter(temp, minOrMax == 'Min'), filter(temp, minOrMax == 'Max'), by = c('days.since.exposed', 'outcome', 'fever', 'meanIncubate'), suffix = c('min', 'max'))%>%
  select(-starts_with('min')) %>% 
  ## Clean up categorical variables so that plot labels are publication quality 
  mutate(fever = factor(fever, levels = unique(fever), labels = paste0('P(fever)=', unique(fever))),  ## Rename levels for nice plotting
         meanIncubate = factor(meanIncubate, levels = unique(meanIncubate), labels = paste0('Mean inc. pd.= ', unique(meanIncubate), 'd')),
         outcome = factor(outcome, levels =c('dFever', 'dRisk', 'aFever', 'aRisk', 'mb', 'mf', 'mr', 'nd'), 
                          labels =(c('stopped: departure fever screen', "stopped: departure risk screen", "stoppd: arrival fever screen",
                                     "stopped: arrival risk screen", 'cleared: missed both', 'cleared: missed fever', 'cleared: missed risk', 'cleared: not detectable')))) %>%
  ## Plot
  ggplot()+
  geom_ribbon(aes(x = days.since.exposed, ymin = yymin, ymax = yymax, fill = outcome))+
  facet_grid(fever~meanIncubate) +
  scale_fill_manual(values = cols)+
  theme_minimal() +
  ylab('Prob. exposed individual is detained or cleared')+
  xlab('Days since exposure')-> ribbon
ribbon
ggsave('2020_nCov/Fig1_grid_of_ribbon_plots.pdf', width = 7, height = 4.5, units = 'in')


# - - - - - - - - - - - - - - - - - - - - - - -
# Simulate the fraction of the population caught or missed in a growing epidemic.
# - - - - - - - - - - - - - - - - - - - - - - -
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
incubation.d<-(function(x){pgamma(x, shape = shapeIncubate, scale = scaleIncubate)})


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
  riskMedians = extract_mean_probs(bootList$outcomesRisk)) %>%
  mutate(strategy = c('both', 'fever', 'risk')) -> meanOutcomes
names(meanOutcomes) = c('d.fever', 'd.risk', 'a.fever', 'a.risk', 'm.b', 'm.f', 'm.r', 'nd', 'strategy')



meanOutcomes %>%
  pivot_longer(cols = 1:8) %>%
  mutate(outcome = factor(name, levels = (c('d.fever', 'd.risk', 'a.fever', 'a.risk', 'm.b', 'm.f', 'm.r', 'nd')),
                          labels =(c('stopped: departure fever screen', "stopped: departure risk screen", "stoppd: arrival fever screen",
                                     "stopped: arrival risk screen", 'cleared: missed both', 'cleared: missed fever', 'cleared: missed risk', 'cleared: not detectable')))) %>%
  ggplot()+
  geom_bar(aes(x = strategy, y = value, fill = outcome), stat = 'identity')+
  scale_fill_manual(values = cols) +
  xlab('Screening type')+
  ylab('Fraction screened') +
  theme_bw() -> stackedBars
stackedBars
ggsave('2020_nCov/outcome_stacked_bars.pdf', width = 6, height = 4, units = 'in')  

blank = ggplot()



pdf('2020_nCov/arranged.pdf', width = 12, height = 12)
grid.arrange(ribbon, blank, blank, fracMissed, stackedBars, layout_matrix = matrix(c(1,2,3,4,5,5,6,7,7,8,9,9), byrow = T, nrow = 4, ncol = 3))
dev.off()
