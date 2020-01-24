# Code for 3 x 3 grid

# Source external functions and load packages
source("Code_Screening_model.R")
source("Code_distribution_functions.R")
library(ggplot2)
library(grid)
library(gridExtra)
library(tidyverse)
coltab=c("red","blue","orange","green","purple","black")

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
  incubation.d<-dgamma(x = tSinceExposed, shape = shapeIncubate, scale = scaleIncubate) 
  ## Outputs
  screen.passengers(tSinceExposed, del.d = 1, ff, gg, sd, sa, rd, ra, phi.d, incubation.d, pathogen, relT, split1 = 1)
}
# ## Test function
# get_frac_caught(tSinceExposed = 2, ff = fps[1], gg = .2, R0 = 2, meanToAdmit = 5.5, meanIncubate = 4.5, daysInFlight = 1)
# -------------------------------

# -------------------------------
## Write a function that repeats get_frac_caught over a grid of times since exposure
get_frac_caught_over_time = function(ff, gg, R0, meanToAdmit, meanIncubate = NULL, daysInFlight = NULL, shapeIncubate = NULL, scaleIncubate = 2.73){
  arrive.times=seq(0,15,0.1)
  sapply(arrive.times, function(tt){get_frac_caught(tSinceExposed = tt, ff, gg, R0, meanToAdmit, meanIncubate = 4.5, daysInFlight = 1)}) %>% t() %>% as.data.frame() %>% mutate(
    days.since.exposed = arrive.times)
}
# # Test
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
         clearedMin = aRiskMax, clearedMax = 1) %>%
  select(days.since.exposed, contains('Min'), contains('Max'), fever, meanIncubate) %>%
  ## Pivot to long data frame
  pivot_longer(cols = dFeverMin:clearedMax, names_to = c('outcome', 'minOrMax'), names_pattern = '(\\w+)(M\\w\\w)', values_to = 'yy') -> temp
## Use full join to create columns for time, ymin, ymax, and band type
full_join(filter(temp, minOrMax == 'Min'), filter(temp, minOrMax == 'Max'), by = c('days.since.exposed', 'outcome', 'fever', 'meanIncubate'), suffix = c('min', 'max'))%>%
  select(-starts_with('min')) %>% 
  ## Clean up categorical variables so that plot labels are publication quality 
  mutate(fever = factor(fever, levels = unique(fever), labels = paste0('P(fever)=', unique(fever))),  ## Rename levels for nice plotting
         meanIncubate = factor(meanIncubate, levels = unique(meanIncubate), labels = paste0('Mean inc. pd.= ', unique(meanIncubate), 'd')),
         outcome = factor(outcome, levels =c('dFever', 'dRisk', 'aFever', 'aRisk', 'cleared'), 
                          labels = c('stopped: departure fever screen', "stopped: departure risk screen", "stoppd: arrival fever screen",
                                     "stopped: arrival risk screen", 'cleared'))) %>%
  ## Plot
  ggplot()+
  geom_ribbon(aes(x = days.since.exposed, ymin = yymin, ymax = yymax, fill = outcome))+
  facet_grid(fever~meanIncubate) +
  scale_fill_manual(values = c('darkseagreen2', 'deepskyblue', 'seagreen4', 'royalblue3', 'darksalmon'))+
  theme_minimal() +
  ylab('Prob. exposed individual is detained or cleared')+
  xlab('Days since exposure') -> ribbon
ggsave('2020_nCov/Fig1_grid_of_ribbon_plots.pdf', width = 7, height = 4.5, units = 'in')

