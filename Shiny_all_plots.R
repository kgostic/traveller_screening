## Generate a ribbon plot, a plot of the fraction missed, and stacked bars for any of the following scenarios:
##    Arrival screening only
##    Departure screening only
##    Both

# Source external functions and load packages
source("Code_Screening_model.R")
source("Code_distribution_functions.R")
library(ggplot2)
#library(grid)
library(gridExtra)
library(tidyverse)
#library(doParallel)
cols = c('darkseagreen2', 'deepskyblue', 'seagreen4', 'royalblue3', 'bisque', 'darksalmon', 'brown3', 'firebrick1')





## Parameters for input
arrival_screen = TRUE
departure_screen = TRUE
ff = .8
gg = .1
R0 = 2
meanToAdmit = 5.5
meanIncubate = 5.5
hrs.in.flight = 24
del.d = hrs.in.flight/24 ## days in flight - can be a fraction
shapeIncubate = NULL
scaleIncubate = 2.73


## Get incubation period distribution
if(length(meanIncubate)>0){         ## If mean incubation period specified, calculate shape
  shapeIncubate = meanIncubate/scaleIncubate
}else if(length(shapeIncubate)==0){ ## Else, use explicit shape and scale inputs
  stop('Error - Must specify meanIncubate or shape and scale!')
}
# Set the incubation cdf
incubation.d<-(function(x){pgamma(x, shape = shapeIncubate, scale = scaleIncubate)})


# -------------------------------
# \\\\\\\\\\\\\\\  WRITE FUNCTIONS AND WRAPPERS ///////////////////
# -------------------------------

# -------------------------------
#get_frac_caught = function(tSinceExposed, ff, gg, R0, meanToAdmit, meanIncubate = NULL, dscreen, ascreen, shapeIncubate = NULL, scaleIncubate = 2.73){
get_frac_caught = function(tSinceExposed, ff, gg, R0, meanToAdmit, dscreen, ascreen, incubation.d){
## Calculate the fraction/probability of each screening outcome given a fixed time since exposure
  ## INPUTS
  ## f1 - probability of fever at onset
  ## g1 - probability aware of risk
  ## meanToAdmit - mean days from onset to admission (this is the detectable window) 
  ## meanIncubate - mean incubation period (exposure ot onset). If meanIncubate specified, assume a gamma distribution with scale = 2.73 and calculate appropriate shape. Alternatively, can specify shape and scale explicitly.
  ## OUTPUTS - vector. Fraction caught by: departure risk, departure fever, arrival risk, arrival fever, cleared.
  
  # if(length(meanIncubate)>0){         ## If mean incubation period specified, calculate shape
  #   shapeIncubate = meanIncubate/scaleIncubate
  # }else if(length(shapeIncubate)==0){ ## Else, use explicit shape and scale inputs
  #   stop('Error - Must specify meanIncubate or shape and scale!')
  # }
  # # Get probability that symptoms have developed.
  # incubation.d<-(function(d)pgamma(d, shape = shapeIncubate, scale = scaleIncubate))
  ## Outputs
  screen.passengers(tSinceExposed, del.d, ff, gg, sd, sa, rd, ra, phi.d, incubation.d, pathogen, relative = 0, split1 = 2, arrival_screen = ascreen, departure_screen = dscreen)
}
# ## Test function
# get_frac_caught(tSinceExposed = 2, ff = .7, gg = .2, R0 = 2, meanToAdmit = 5.5, dscreen = TRUE, ascreen = FALSE, incubation.d)
# -------------------------------

# -------------------------------
## Write a function that repeats get_frac_caught over a grid of times since exposure
get_frac_caught_over_time = function(ff, gg, R0, meanToAdmit, ascreen, dscreen, incubation.d){
  arrive.times=seq(0,15,0.1)
  sapply(arrive.times, function(tt){get_frac_caught(tSinceExposed = tt, ff, gg, R0, meanToAdmit, meanIncubate, dscreen, ascreen)}) %>% t() %>% as.data.frame() %>% mutate(
    days.since.exposed = arrive.times)
}
# Test
# get_frac_caught_over_time(ff = .7, gg = .2, R0 = 2, meanToAdmit = 5.5, ascreen = TRUE, dscreen = FALSE, incubation.d)
# -------------------------------

# -------------------------------
# Simulate the fraction of the population caught or missed in a growing epidemic.
one_bootstrap = function(){
## For each individual in the population, draw times since exposre from the appropriate distribution
infAge=sapply(c(1:popn),function(x){exposure.distn(runif(1, 0, 1),pathogen,flat=0)})

## Given times since exposure, screen using fever, risk or both
screenWrapper = function(x, ff, gg){
  screen.passengers(x, del.d, ff, gg, sd, sa, rd, ra, phi.d, incubation.d, pathogen, relative = 0, split1 = 2, arrival_screen, departure_screen)}

fever.sample.size = fg.distn("f",pathogen)[2]
fever.prob = fg.distn("f",pathogen)[1]
risk.sample.size = fg.distn("g",pathogen)[2]
risk.prob = fg.distn("g",pathogen)[1]

f0 = rbinom(popn,fever.sample.size,fever.prob)/fever.sample.size                          ## Draw invidual prob fever
g0 = rbinom(popn, risk.sample.size, risk.prob)/risk.sample.size                           ## Draw individual prob known risk 

## Fever and risk screening
# outcomesBoth=sapply(infAge,function(x){
#   f0=rbinom(1,fever.sample.size,fever.prob)/fever.sample.size                                 ## Draw binomial probability of fever
#   g0=rbinom(1,risk.sample.size, risk.prob)/risk.sample.size                                   ## Draw binomial probability of known risk  
#   screenWrapper(x, f0, g0)})           
outcomesBoth=mapply(FUN = screenWrapper, x = infAge, ff = f0, gg = g0)
## Output individual probability missed
pCaught_both = colSums(outcomesBoth[1:4,])
caughtBoth = sapply(pCaught_both,function(x){ifelse(x<runif(1, 0, 1),0,1)})                   ## Draw whether individual was missed


## Fever only
# outcomesFever=sapply(infAge,function(x){
#   f0=rbinom(1,fever.sample.size,fever.prob)/fever.sample.size  ## Draw binomial probability of fever
#   g0=0  ## Draw binomial probability of known risk  
#   screenWrapper(x, f0, g0)})     
outcomesFever=mapply(FUN = screenWrapper, x = infAge, ff = f0, gg = 0)
pCaught_fever = colSums(outcomesFever[1:4,])                                                    ## Output individual probability missed
caughtFever = sapply(pCaught_fever,function(x){ifelse(x<runif(1, 0, 1),0,1)})                   ## Draw whether individual was missed


## Risk only
# outcomesRisk=sapply(infAge,function(x){
#   f0=0  ## Draw binomial probability of fever
#   g0=rbinom(1,risk.sample.size, risk.prob)/risk.sample.size                                   ## Draw binomial probability of known risk  
#   screenWrapper(x, f0, g0)})                                                                  ## Output individual probability missed
outcomesRisk=mapply(FUN = screenWrapper, x = infAge, ff = 0, gg = g0)
pCaught_risk = colSums(outcomesRisk[1:4,])
caughtRisk = sapply(pCaught_risk,function(x){ifelse(x<runif(1, 0, 1),0,1)})                   ## Draw whether individual was missed

## Return frac.missed.both, frac.missed.fever, frac.missed.risk
# return(c(frac.missed.both = 1- sum(caughtBoth)/popn,
#          frac.missed.fever = 1- sum(caughtFever)/popn,
#          frac.missed.risk = 1- sum(caughtRisk)/popn))

# if(outtype == 'caught'){
#   return(c(frac.missed.both = 1- sum(caughtBoth)/popn,
#                     frac.missed.fever = 1- sum(caughtFever)/popn,
#                     frac.missed.risk = 1- sum(caughtRisk)/popn))
# }else if(outtype == 'meanOutcomes'){
return(list(outcomesBoth = rowMeans(outcomesBoth),
            outcomesFever = rowMeans(outcomesFever),
            outcomesRisk = rowMeans(outcomesRisk),
            caught = c(frac.missed.both = 1- sum(caughtBoth)/popn,
                                            frac.missed.fever = 1- sum(caughtFever)/popn,
                                            frac.missed.risk = 1- sum(caughtRisk)/popn)))
}
# ## Test
# one_bootstrap()
# -------------------------------
replicate(nboot, one_bootstrap()) -> bootList  ## Simulate population outcomes











# -------------------------------
#  \\\\\\\\\\\\\ Make Plots /////////////////
# -------------------------------
# aa = Sys.time()

# bb = Sys.time()
# bb-aa

# cl <- makeCluster(detectCores())
# registerDoParallel(cl)
# foreach(i=1:nboot,
#         .combine = 'cbind',
#         .packages = 'dplyr') %dopar% {one_bootstrap()} -> xx
# xx
# cc = Sys.time()


# -------------------------------
### Fraction missed
# -------------------------------
sapply(xx[4,], function(yy){yy})  %>% t() %>% as.data.frame() %>%
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
#fracMissed
#ggsave('2020_nCov/Shiny_fraction_missed.pdf', width = 4, height = 4, units = 'in')  

# -------------------------------
## Stacked barplot
# -------------------------------
bind_rows(
  bothMedians = sapply(bootList[1,], function(yy){yy}) %>% rowMeans,
  feverMedians = sapply(bootList[2,], function(yy){yy}) %>% rowMeans,
  riskMedians = sapply(bootList[3,], function(yy){yy}) %>% rowMeans
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
  scale_fill_manual(values = cols) +
  xlab('Screening type')+
  ylab('Fraction screened') +
  theme_bw() -> stackedBars
#stackedBars
#ggsave('2020_nCov/Shiny_stacked_bars.pdf', width = 6, height = 4, units = 'in')  




# -------------------------------
# Ribbon Plot
# -------------------------------
get_frac_caught_over_time(ff, gg, R0, meanToAdmit, arrival_screen, departure_screen, incubation.d) %>%
  ## Specify minim and maximum y values in each band
  mutate(dFeverMin = 0, dFeverMax = dFeverMin + caught.dpt.fever,
         dRiskMin = dFeverMax, dRiskMax = dRiskMin + caught.dpt.risk,
         aFeverMin = dRiskMax, aFeverMax = aFeverMin + caught.arv.fever,
         aRiskMin = aFeverMax, aRiskMax = aRiskMin +caught.arv.risk,
         mbMin = aRiskMax, mbMax = mbMin+missed.both,
         mfMin = mbMax, mfMax = mfMin+missed.fever.only,
         mrMin = mfMax, mrMax = mrMin+missed.risk.only,
         ndMin = mrMax, ndMax = 1) %>%
  select(days.since.exposed, contains('Min'), contains('Max')) %>%
  ## Pivot to long data frame
  pivot_longer(cols = dFeverMin:ndMax, names_to = c('outcome', 'minOrMax'), names_pattern = '(\\w+)(M\\w\\w)', values_to = 'yy') -> temp
## Use full join to create columns for time, ymin, ymax, and band type
full_join(filter(temp, minOrMax == 'Min'), filter(temp, minOrMax == 'Max'), by = c('days.since.exposed', 'outcome'), suffix = c('min', 'max'))%>%
  select(-starts_with('min')) %>% 
  ## Clean up categorical variables so that plot labels are publication quality 
  mutate(outcome = factor(outcome, levels =c('dFever', 'dRisk', 'aFever', 'aRisk', 'mb', 'mf', 'mr', 'nd'), 
                          labels =(c('stopped: departure fever screen', "stopped: departure risk screen", "stoppd: arrival fever screen",
                                     "stopped: arrival risk screen", 'cleared: missed both', 'cleared: missed fever', 'cleared: missed risk', 'cleared: not detectable')))) %>%
  ## Plot
  ggplot()+
  geom_ribbon(aes(x = days.since.exposed, ymin = yymin, ymax = yymax, fill = outcome))+
  scale_fill_manual(values = cols)+
  theme_minimal() +
  ylab('Prob. exposed individual is detained or cleared')+
  xlab('Days since exposure') +
  theme(legend.position = 'none')-> ribbon
# ribbon
# ggsave('2020_nCov/Shiny_ribbon_plot.pdf', width = 7, height = 4.5, units = 'in')





# 
# 
# 
# 
# 
# ## Set a grid of fever and incubation inputs
# input_grid = expand.grid(ffs = c(.6, .75, .95),
#                          meanIncs = c(4, 5.5, 7))
# 
# 
# apply(X = input_grid, MARGIN = 1, FUN = function(ii){get_frac_caught_over_time(ff = ii[1], gg = .25, R0 = 2, meanToAdmit = ii[2])}) %>%                                                                        ## This would bind outputs into a long data frame instead of a list
#   bind_rows() %>% as.tbl() %>%
#   mutate(fever = rep(input_grid$ffs, each = nrow(.)/nrow(input_grid)),
#          meanIncubate = rep(input_grid$meanIncs, each = nrow(.)/nrow(input_grid))) -> gridOutputs
# 
# ## Reformat for plotting
# gridOutputs %>%
#   ## Specify minim and maximum y values in each band
#   mutate(dFeverMin = 0, dFeverMax = dFeverMin + caught.dpt.fever,
#          dRiskMin = dFeverMax, dRiskMax = dRiskMin + caught.dpt.risk,
#          aFeverMin = dRiskMax, aFeverMax = aFeverMin + caught.arv.fever,
#          aRiskMin = aFeverMax, aRiskMax = aRiskMin +caught.arv.risk,
#          mbMin = aRiskMax, mbMax = mbMin+missed.both,
#          mfMin = mbMax, mfMax = mfMin+missed.fever.only,
#          mrMin = mfMax, mrMax = mrMin+missed.risk.only,
#          ndMin = mrMax, ndMax = 1) %>%
#   select(days.since.exposed, contains('Min'), contains('Max'), fever, meanIncubate) %>%
#   ## Pivot to long data frame
#   pivot_longer(cols = dFeverMin:ndMax, names_to = c('outcome', 'minOrMax'), names_pattern = '(\\w+)(M\\w\\w)', values_to = 'yy') -> temp
# ## Use full join to create columns for time, ymin, ymax, and band type
# full_join(filter(temp, minOrMax == 'Min'), filter(temp, minOrMax == 'Max'), by = c('days.since.exposed', 'outcome', 'fever', 'meanIncubate'), suffix = c('min', 'max'))%>%
#   select(-starts_with('min')) %>% 
#   ## Clean up categorical variables so that plot labels are publication quality 
#   mutate(fever = factor(fever, levels = unique(fever), labels = paste0('P(fever)=', unique(fever))),  ## Rename levels for nice plotting
#          meanIncubate = factor(meanIncubate, levels = unique(meanIncubate), labels = paste0('Mean inc. pd.= ', unique(meanIncubate), 'd')),
#          outcome = factor(outcome, levels =c('dFever', 'dRisk', 'aFever', 'aRisk', 'mb', 'mf', 'mr', 'nd'), 
#                           labels =(c('stopped: departure fever screen', "stopped: departure risk screen", "stoppd: arrival fever screen",
#                                      "stopped: arrival risk screen", 'cleared: missed both', 'cleared: missed fever', 'cleared: missed risk', 'cleared: not detectable')))) %>%
#   ## Plot
#   ggplot()+
#   geom_ribbon(aes(x = days.since.exposed, ymin = yymin, ymax = yymax, fill = outcome))+
#   facet_grid(fever~meanIncubate) +
#   scale_fill_manual(values = cols)+
#   theme_minimal() +
#   ylab('Prob. exposed individual is detained or cleared')+
#   xlab('Days since exposure')-> ribbon
# ribbon
# ggsave('2020_nCov/Fig1_grid_of_ribbon_plots.pdf', width = 7, height = 4.5, units = 'in')
# 
# 
# 
# 
