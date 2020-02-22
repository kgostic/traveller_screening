library(ggplot2)
library(grid)
library(gridExtra)
library(tidyverse)
library(viridis)
set.seed(38)
reset = FALSE ## If FALSE, load saved .RData files to generate Figs. 3-4. 
             ## If TRUE, rerun simulations before generating Figs. 3-4. Simulations take a few minutes to run.


## ------------------------------------------------------------------------------------------------------------------------
## 
## \\\\\\\\\\\\\\  SECTION 1: SETUP AND FUNCTION DEFINITIONS
## 
## ------------------------------------------------------------------------------------------------------------------------

## ------------------------------------------------------------
## \\\\\\\\\\\\\\ DEFINE GLOBAL VARIABLES ////////////////// ##
## ------------------------------------------------------------
flight.hrs = 24
flight.time = flight.hrs/24  ## Set assumed duration of flight (all figures)


###  SET FIG. 2 INPUTS
NOTsubclinical = c(.5, .75, .95) ## Set the fraction of cases with detectable fever or cough (1- assumed fraction subclincal cases) (Fig. 2)
incubationMeans = c(4.5, 5.5, 6.5) ## Set mean incubation periods (in days) (Fig. 2)


## SET FIG. 3 INPUTS
low.vals = c(gg = .05, f.sens = .6, g.sens = .05, mInc = 4.5, R0 = 1.5, meanToAdmit = 3)  ## Define low and high ranges for LHS (as in Table 1)
high.vals = c(gg = .40, f.sens = .90, g.sens = .25, mInc = 6.5, R0 = 4, meanToAdmit = 7)
## gg = prob. aware of exposure risk
## f.sens = sensitivity of thermal scanners for fever [0, 1]
## g.sens = "sensitivity" (i.e. prob of thruthful self-reporting) on risk questionnaire [0, 1]
## mInc = mean incubation period (days)
## R0 - R0
## meanToAdmit - mean days from onset of symptoms to isolation (affects maximum time from exposure to iniiation of travel)
nsim = 1000    ## n simulation trials, each corresponding to a unique LHS paramter set
popn = 30       ## population size of infected travelers
incVariance = 2.25^2 ## Fixed variance of the incubation period distribution
## Labels for detection categories
cat.labels = c('detected: departure fever screen', "detected: departure risk screen", "detected: arrival fever screen",
               "detected: arrival risk screen", 'missed: had both', 'missed: had fever', 'missed: had risk awareness', 'missed: undetectable')
## Set colors
# cols = c('darkseagreen2', 'deepskyblue', 'seagreen4', 'royalblue3', 'bisque', 'brown4', 'salmon2', 'firebrick1')
cols = c('skyblue', 'deepskyblue', 'royalblue1', 'mediumblue', 'khaki1', 'goldenrod1', 'darkorange', 'firebrick1') # Move to a more colorblind friendly palette





## -----------------------------------------------------
## \\\\\\\\\\\\\\ DEFINE FUNCTIONS ////////////////// ##
## -----------------------------------------------------

# ---------
get.inc.dist <- function(meanInc){
  ## Input mean of gamma-distributed incubation times
  ## Output a function, incCDF, that inputs a time since exposure and outputs the cumulative density
    scale.in = incVariance/meanInc
    shape.in = meanInc/scale.in
  
  return(incCDF <- function(d){pgamma(d, shape = shape.in, scale = scale.in)})
}

# ---------
# Calculate infection age distributions for a growing epidemic
pdf.cdf.travel<-function(x,r0,gen.time,type){
  # Input: x - time since exposure, r0 - pathogen R_0 value, gen.time - pathogen generation time,
  # type - [pdf = 1, cdf = 2]
  # Output: if type = 1, output the density of infection age x among a population of 
  # exposed individuals attempting travel. If type = 2, output the cumulative density.
  alpha=r0/gen.time
  f01<-function(tt){exp(-alpha*tt)}                  #dI/dt
  f11<-function(tt){(1/alpha)*(1-exp(-alpha*tt))}    #I(t)
  if(type==1){
    sapply(x,function(a){if(a<gen.time){f01(a)/f11(gen.time)}else{0}}) #pdf
  }else{
    sapply(x,function(a){if(a<gen.time){f11(a)/f11(gen.time)}else{1}}) #cdf
  }
}
## --------------

## -------------------------
exposure.distn = function(x,r0, meanToAdmit, meanIncubate){
  ## Draw from the infection age cdf using an input uniform random variable, x
  (data.frame(xx = seq(0, 25, by = 0.1)) %>%             # Evaluate cdf at a range of values
     mutate(c.dens = pdf.cdf.travel(xx, r0, meanToAdmit+meanIncubate, type = 2)) %>%
     filter(c.dens>x) %>%                                # Extract the first value at which the cum.dens > x
     pull(xx))[1]
}
## -------------------------



## --------------
screen.passengers = function(d, del.d, f, g, sd = 1, sa =1, rd = 1, ra = 1, 
                             incubation.d, relative=0, split1=0, arrival_screen, departure_screen, frac_evade=0){
  ## Arrival Screening Decision Tree Model
  #   INPUTS:
  #   d = days since onset
  #   del.d = days spent in flight (can be a fraction)
  #   f = probability that patients present with fever at onset
  #   g = probability that patients are aware of exposure risk factors
  #   sd = symptom screen effectiveness on departure
  #   sa = symptom screen effectiveness on arrival
  #   rd = risk factor screen effectiveness on departure
  #   ra = risk factor screen effectiveness on arrival
  #   incubation.d = call to function describing pdf of time from exposure to onset
  #   relative: this takes a logical value. The default is 0 or FALSE, which tells the function to return
  #                the absolute proprtion of travellers detected at arrival. If relative is set to 1 or TRUE, 
  #                the script will return the proportion of infected travellers detected at arrival given 
  #                that they were missed during departure screening
  #                (i.e. [# detected at arrival]/[# missed at departure])
  #   split1: If 1, the function outputs the following vector:
  #                [stopped.at.departure.fever,
  #                 stopped.at.departure.risk,
  #                 stopped.at.arrival.fever,
  #                 stopped.at.arrival.risk,
  #                 cleared.at.arrival]
  #           If 0 (default), the function outputs only: 
  #                 1-[cleared.at.arriva]
  #           If 2, outputs  [stopped.at.departure.fever,
  #                 stopped.at.departure.risk,
  #                 stopped.at.arrival.fever,
  #                 stopped.at.arrival.risk,
  #                 cleared.at.arrival,
  #                 missed.both,
  #                 missed.fever,
  #                 missed.risk,
  #                 not.detectable]
  #   arrival_screen: logical, should the model screen on arrival
  #   departure_screen: logical as above
  #   frac_evade: fration of infected travellers to intentially evade screening. Default is 0.
  #
  #           ALL OUTPUTS ARE GIVEN AS THE PROPORTION OF INFECTED TRAVELLERS IN EACH OUTCOME CLASS
  
  ##Define an internal function to pass travellers in any detection class
  # through the model
  screen.cases = function(case){
    
    #Case 1 = has fever and aware of risk factors
    #Case 2 = has fever and NOT aware of risk factors
    #Case 3 = does NOT have fever and aware of risk factors
    #Case 4 = does NOT have fever and NOT aware of risk factors
    #The proportion of travellers that fall into each case (detection class)
    #is determined by values of f and g
    
    #Split in to groups with and without symptoms upon departure
    Sd = incubation.d(d) #Incubation.d(d) is the CDF of exposure -> onset
    NSd = (1-incubation.d(d))
    
    #SYMPTOM SCREEN
    #First screen symptomatic patients
    #Sd.... denotes symptom onset at departure
    if(!departure_screen | case %in% c(3,4)){          #If no departure screen, or no fever present at onset, skip symptom screen
      Sd.sspass = Sd #If no fever at onset, all pass           #Move on
      Sd.ssfail = 0                                            #Detained
    }else{                                             #If fever at onset, perform symptom screen
      Sd.sspass = Sd*(1-sd) #(1-sd) pass                       #Move on 
      Sd.ssfail = Sd*sd     # (sd) fail                        #Detained
    }
    
    #No symptom screen for asymptomatic patients (all pass)
    #NSd.... denotes NO symptom onset at departure
    NSd.sspass = NSd                                           #Move on
    NSd.ssfail = 0                                             #Detained
    
    
    ## RISK SCREEN - only those in sspass categories move on
    if(!departure_screen | case %in% c(2,4)){        ## Don't screen
      Sd.sspass.rspass = Sd.sspass                        # Move on
      Sd.sspass.rsfail = 0                                # Detained
      NSd.sspass.rspass = NSd.sspass                      # Move on
      NSd.sspass.rsfail = 0                               # Detained
    }else{                                            ## Do screen
      Sd.sspass.rspass = Sd.sspass*(1-rd)                  # Move on
      Sd.sspass.rsfail = Sd.sspass*rd                      # Detained
      NSd.sspass.rspass = NSd.sspass*(1-rd)                # Move on
      NSd.sspass.rsfail = NSd.sspass*rd                    # Detained
    }
    
    
    ### \\\\\\\\\\\\\\\\\\\\\\\\\\\\ FLIGHT TAKES OFF /////////////////////////////// ###
    #TOTAL FLYING AND DETAINED
    #Passengers can be stoped if (1) Symptomatic and fail symptom screen, (2) Asymptomatic and fail SS (this should be 0),
    #  (3) Symptomatic, pass SS but fail RS, (4) Asymptomatic and pass SS but fail RS
    stopped.at.departure.fever = Sd.ssfail + NSd.ssfail  
    stopped.at.departure.risk = Sd.sspass.rsfail + NSd.sspass.rsfail
    stopped.at.departure = stopped.at.departure.fever+stopped.at.departure.risk
    #Passengers are only cleared if they pass both screens
    cleared.at.departure = Sd.sspass.rspass + NSd.sspass.rspass
    
    #With symptoms at departure
    Sd.arrive = Sd.sspass.rspass
    #Without symptoms at departure
    NSd.arrive = NSd.sspass.rspass
    
    ### \\\\\\\\\\\\\\\\\\\\\\\\\\\\ FLIGHT LANDS /////////////////////////////// ###
    #SYMPTOM DEVELOPMENT IN FLIGHT
    #calculate the conditional probability of developing symptoms on flight given no
    #symptoms at departure
    p.new.symptoms = ifelse(arrival_screen, ((incubation.d(d+del.d)-incubation.d(d))/(1-incubation.d(d))), 0)
    Sa.arrive = NSd.arrive*p.new.symptoms #Sa.... denotes symptoms on arrival, but not at departure
    NS.arrive = NSd.arrive*(1-p.new.symptoms) #NS... denotes no symptoms on arrival
    #No modulation of Sd.arrive because this group was already symptomatic
    #Now there are three branches: Sa (symptoms on arrival), Sd (symptoms on departure), NS (no symptoms)
    #Remember we are still in a function that applies this decision tree to each of 4 f,g cases
    
    #SYMPTOM SCREEN
    #First screen symptomatic patients
    if(!arrival_screen | case %in% c(3,4)){          #If no arrival screen, or no fever present at onset, skip symptom screen
      Sd.arrive.sspass = Sd.arrive                                             #Move on 
      Sd.arrive.ssfail = 0                                                     #Detained
      Sa.arrive.sspass = Sa.arrive                                             #Move on
      Sa.arrive.ssfail = 0                                                     #Detained
    }else{                                            # Do screen
      Sd.arrive.sspass = Sd.arrive*(1-sa)    #(1-sa) pass                      #Move on 
      Sd.arrive.ssfail = Sd.arrive*sa        # (sa) faill                      #Detained
      Sa.arrive.sspass = Sa.arrive*(1-sa)                                      #Move on
      Sa.arrive.ssfail = Sa.arrive*sa                                          #Detained
    }
    
    
    #No symptom screen for asymptomatic patients (all pass)
    NS.arrive.sspass = NS.arrive                                                   #Move on
    NS.arrive.ssfail = 0                                                           #Detained
    
    
    #RISK SCREEN
    if(!arrival_screen | case %in% c(2,4)){          #If no arrival screen, or no risk awareness, don't screen
      Sd.arrive.sspass.rspass = Sd.arrive.sspass                           # Move on
      Sd.arrive.sspass.rsfail = 0                                          # Detained
      Sa.arrive.sspass.rspass = Sa.arrive.sspass                           # Move on
      Sa.arrive.sspass.rsfail = 0                                          # Detained
      NS.arrive.sspass.rspass = NS.arrive.sspass                           # Move on
      NS.arrive.sspass.rsfail = 0                                          # Detained
    }else{                                           #If passengers are aware of risk factors, perform screen
      Sd.arrive.sspass.rspass = Sd.arrive.sspass*(1-ra)                    # Move on
      Sd.arrive.sspass.rsfail = Sd.arrive.sspass*ra                        # Detained
      Sa.arrive.sspass.rspass = Sa.arrive.sspass*(1-ra)                    # Move on
      Sa.arrive.sspass.rsfail = Sa.arrive.sspass*ra                        # Detained
      NS.arrive.sspass.rspass = NS.arrive.sspass*(1-ra)                    # Move on
      NS.arrive.sspass.rsfail = NS.arrive.sspass*ra                        # Detained
    }
    
    
    #TOTAL DETAINED AND FREE
    #Passengers in each of three classes (Sd, Sa, NS) are detained if they fail SS or pass SS and fail RS
    stopped.at.arrival.fever=(Sd.arrive.ssfail + Sa.arrive.ssfail + NS.arrive.ssfail)
    stopped.at.arrival.risk=(Sd.arrive.sspass.rsfail + Sa.arrive.sspass.rsfail + NS.arrive.sspass.rsfail)
    stopped.at.arrival = stopped.at.arrival.risk+stopped.at.arrival.fever
    
    #Passengers in each of three classes are cleared only if they pass both screens
    cleared.at.arrival = (Sd.arrive.sspass.rspass + Sa.arrive.sspass.rspass + NS.arrive.sspass.rspass)
    
    
    #specify whether relative or absolute
    # Give proportion caught
    if(relative==1){
      outputs=1-cleared.at.arrival/cleared.at.departure #
      names(outputs) = 'excess.frac.caught.at.arrival'
    }else{
      outputs=1-cleared.at.arrival #overall fraction missed
      names(outputs) = 'p.missed'
    }
    
    if(split1%in%c(1,2)){
      #den1=1#cleared.at.departure
      if(case %in% c(1,2)){ ## If symptoms could have developed
        outputs=c(stopped.at.departure.fever, stopped.at.departure.risk, stopped.at.arrival.fever, stopped.at.arrival.risk, 
                  Sd.arrive.sspass.rspass+Sa.arrive.sspass.rspass, NS.arrive.sspass.rspass)
      }else{
        outputs=c(stopped.at.departure.fever, stopped.at.departure.risk, stopped.at.arrival.fever, stopped.at.arrival.risk, 0, cleared.at.arrival)
      }
      names(outputs) = c('caught.dpt.fever', 'caught.dpt.risk', 'caught.arv.fever', 'caught.arv.risk', 'missed.sd', 'missed.nsd')
    }
    
    return(outputs)
    
  }
  
  
  ff1=f
  gg1=g
  
  #Define the proportion of travellers that fall into each detection class (case)
  cases1 = c((1-frac_evade)*c(ff1*gg1, ff1*(1-gg1), (1-ff1)*gg1, (1-ff1)*(1-gg1)))
  
  if(split1 == 2){
    c1 = cases1[1]*screen.cases(1) %>% as.numeric() # Had both
    c2 = cases1[2]*screen.cases(2) %>% as.numeric() # Had fever only
    c3 = cases1[3]*screen.cases(3) %>% as.numeric() # Had risk only
    c4 = cases1[4]*screen.cases(4) %>% as.numeric() # Had neither (not detectable)
    
    return(c('caught.dpt.fever' = c1[1]+c2[1]+c3[1]+c4[1], 
             'caught.dpt.risk'= c1[2]+c2[2]+c3[2]+c4[2], 
             'caught.arv.fever' = c1[3]+c2[3]+c3[3]+c4[3], 
             'caught.arv.risk' = c1[4]+c2[4]+c3[4]+c4[4], 
             'missed.both' = c1[5],
             'missed.fever.only' = c2[5],
             'missed.risk.only' = c3[6]+c1[6],
             'not.detectable' = c2[6]+c4[6],
             'evaded.screening' = frac_evade))
  }
  
  #Run the screen.cases function for the appropriate case and weight by the 
  #  appropriate proportion of travellers
  cases1[1]*screen.cases(1)+cases1[2]*screen.cases(2)+cases1[3]*screen.cases(3)+cases1[4]*screen.cases(4)
  
}
## --------------


# -------------------------------
# Wrapper for screen.passengers: Calculate the fraction/probability of each screening outcome given a fixed time since exposure
# Used in Fig. 2
get_frac_caught = function(tSinceExposed, ff, gg, dscreen, ascreen, incubation.d, frac_evaded){
  ## Calculate the fraction/probability of each screening outcome given a fixed time since exposure
  ## INPUTS
  ## f1 - probability of fever at onset
  ## g1 - probability aware of risk
  ## meanToAdmit - mean days from onset to admission (this is the detectable window) 
  ## meanIncubate - mean incubation period (exposure ot onset). 
  ## OUTPUTS - vector. Fraction caught by: departure risk, departure fever, arrival risk, arrival fever, cleared.
  screen.passengers(tSinceExposed, del.d=flight.time, ff, gg, sd=.7, sa=.7, rd=.25, ra=.25, incubation.d, relative = 0, split1 = 2, arrival_screen = ascreen, departure_screen = dscreen, frac_evaded)
}
# ## Test function
# get_frac_caught(tSinceExposed = 2, ff = .7, gg = .2, dscreen = TRUE, ascreen = TRUE, incubation.d = get.inc.dist(meanInc = 4.5), frac_evaded = 0)
# -------------------------------

# -------------------------------
## Wrapper for get_frac_caught
## Repeats get_frac_caught over a grid of times since exposure
## Used in Fig. 2 analyses
get_frac_caught_over_time = function(ff, gg, ascreen, dscreen, incubation.d, frac_evaded = 0){
  arrive.times=seq(0,15,0.1)
  sapply(arrive.times, function(tt){get_frac_caught(tSinceExposed = tt, ff, gg, dscreen, ascreen,incubation.d, frac_evaded)}) %>% t() %>% as.data.frame() %>% mutate(
    days.since.exposed = arrive.times)
}
# Test
# get_frac_caught_over_time(ff = .7, gg = .2, ascreen = TRUE, dscreen = TRUE, incubation.d = get.inc.dist(meanInc = 4.5), frac_evaded = 0)
# -------------------------------

# -------------------------------
# Wrapper that simulates a population of individuals, each with different times since exposure
# Simulate the fraction of the population caught or missed in a growing epidemic.
# Below, call this function repeatedly with different parameter inputs. (Fig. 3)
one_sim = function(meanInc, R0, f0, g0, f.sens, g.sens, gg, del.d, as, ds, meanToAdmit, growing = TRUE){
  
  ## For each individual in the population, draw times since exposure from the appropriate distribution
  if(growing){
  infAge=sapply(c(1:popn),function(x){exposure.distn(runif(1, 0, 1),r0 = R0, meanToAdmit = meanToAdmit, meanIncubate = meanInc)})
  }else if(!growing){
  infAge = runif(n = popn, min = 0, max = meanToAdmit+meanInc)
  }
  
  this.inc.cdf = get.inc.dist(meanInc)           ## Set incubation period distribution
  
  ## Get probability of a given screening outcome for each traveller in the population
  outcomes=sapply(infAge, FUN = function(x){screen.passengers(x, del.d, f0, g0, f.sens, f.sens, g.sens, g.sens, this.inc.cdf, relative = 0, split1 = 2, arrival_screen=as, departure_screen=ds)})
  ## Output individual probability caught
  pCaught = colSums(outcomes[1:4,])
  binaryOutcome = sapply(pCaught,function(x){ifelse(x<runif(1, 0, 1),0,1)})                   ## Draw whether the individual in question was detected
  ## Output a logical vector: were the first n individauls detained?
  firstNCaught = sapply(1:popn, function(nn) sum(1-binaryOutcome[1:nn]) == 0)
  
  return(list(outcomes = rowMeans(outcomes),
              caught = c(frac.missed.both = 1- sum(binaryOutcome)/popn),
              indivCaught = binaryOutcome,
              firstNCaught = firstNCaught,
              expFracCaught = mean(pCaught)))
}
## Test
# one_sim(meanInc = 5.5, R0 = 3, f0 = .7, g0 = .2, f.sens = .7, g.sens = .25, gg = .1, del.d = 1, as = TRUE, ds = FALSE, meanToAdmit = 6)
# -------------------------------








## ------------------------------------------------------------------------------------------------------------------------
## 
## \\\\\\\\\\\\\\  SECTION 2: RUN ANALYSES AND PLOT
## 
## ------------------------------------------------------------------------------------------------------------------------

## ----------------------------------------------------
## \\\\\\\\\\\\\\       Fig. 2  & associated supplementary figs   ////////////////// ##
## ----------------------------------------------------
## Fig. 2. Plot the individual probability of different screening outcomes vs. times since exposure
##   Run across a grid of probabilities of detectable symptoms and mean incubation periods
##   Assume no one evades screening
## Set grid of values to test for ff (fraction with symptoms, 1-frac subclinical), and mean incubation (days)
input_grid = expand.grid(ffs = NOTsubclinical, meanIncs = incubationMeans)
input_grid
## Wrapper to repeat fig. 2 analysis across input grid parameter values
gridWrapper = function(ff.in, mInc.in){
  incFun = get.inc.dist(meanInc = mInc.in)
  get_frac_caught_over_time(ff = ff.in, gg = 0.2, ascreen = TRUE, dscreen = TRUE, incubation.d = incFun, frac_evaded = 0)
}


apply(X = input_grid, MARGIN = 1, FUN = function(ii){gridWrapper(ii[1], ii[2])}) %>%                                                                
  bind_rows() %>% as.tbl() %>%
  mutate(fever = rep(input_grid$ffs, each = nrow(.)/nrow(input_grid)),
         meanIncubate = rep(input_grid$meanIncs, each = nrow(.)/nrow(input_grid))) -> gridOutputs

## Reformat for plotting
gridOutputs %>%
  mutate(dFeverMin = 0, dFeverMax = dFeverMin + caught.dpt.fever,
         dRiskMin = dFeverMax, dRiskMax = dRiskMin + caught.dpt.risk,
         aFeverMin = dRiskMax, aFeverMax = aFeverMin + caught.arv.fever,
         aRiskMin = aFeverMax, aRiskMax = aRiskMin +caught.arv.risk,
         mbMin = aRiskMax, mbMax = mbMin+missed.both,
         mfMin = mbMax, mfMax = mfMin+missed.fever.only,
         mrMin = mfMax, mrMax = mrMin+missed.risk.only,
         ndMin = mrMax, ndMax = ndMin+not.detectable,
         ndeMin = ndMax, ndeMax = ndeMin+evaded.screening) %>%
  select(days.since.exposed, fever, meanIncubate, contains('Min'), contains('Max')) %>%
  ## Pivot to long data frame
  pivot_longer(cols = dFeverMin:ndeMax, names_to = c('outcome', 'minOrMax'), names_pattern = '(\\w+)(M\\w\\w)', values_to = 'yy') -> temp

## Use full join to create columns for time, ymin, ymax, and band type
full_join(filter(temp, minOrMax == 'Min'), filter(temp, minOrMax == 'Max'), by = c('days.since.exposed', 'fever', 'meanIncubate', 'outcome'), suffix = c('min', 'max'))%>%
  select(-starts_with('min')) %>% 
  filter(outcome !='nde') %>%
  ## Clean up categorical variables so that plot labels are publication quality 
  mutate(fever = factor(fever, levels = rev(unique(fever)), labels = rev(paste0((1-unique(fever))*100,"% subclinical"))),  ## Rename levels for nice plotting
         meanIncubate = factor(meanIncubate, levels = unique(meanIncubate), labels = paste0('Mean incubation ', unique(meanIncubate), 'd')),
         outcome = factor(outcome, levels = rev(c('dFever', 'dRisk', 'aFever', 'aRisk', 'mb', 'mf', 'mr', 'nd')), 
                          labels =(rev(cat.labels)))) -> rib
blackline <- filter(rib,outcome=="detected: arrival risk screen") # Extract height of dotted line


## ---------------------------------Plot Fig. 2 ---------------------------------
ggplot(rib)+
  geom_ribbon(aes(x = days.since.exposed, ymin = yymin, ymax = yymax, fill = outcome))+
  facet_grid(fever~meanIncubate) +
  scale_fill_manual(values = cols[8:1])+
  theme_bw() +
  geom_line(data=blackline,aes(x=days.since.exposed,y=yymax),lty=2)+
  ylab('Percentage of exposed individuals detected or missed')+
  scale_y_continuous(breaks = seq(0,1,.25),labels=paste(seq(0,100,25),"%",sep=""))+
  xlab('Days since exposure')  
ggsave('2020_nCov/Fig2_grid_of_ribbon_plots.png', width = 8, height = 4.5, units = 'in')
write.csv(rib, file = '2020_nCov/Fig2_sourceData.csv', row.names = FALSE)




## --------------------------------- Fig 2. Supplemtary figure 1. Departure screening only. ---------------------------------
gridWrapper = function(ff.in, mInc.in){
  incFun = get.inc.dist(meanInc = mInc.in)
  get_frac_caught_over_time(ff = ff.in, gg = 0.2, ascreen = FALSE, dscreen = TRUE, incubation.d = incFun, frac_evaded = 0)
}
apply(X = input_grid, MARGIN = 1, FUN = function(ii){gridWrapper(ii[1], ii[2])}) %>%                                                                
  bind_rows() %>% as.tbl() %>%
  mutate(fever = rep(input_grid$ffs, each = nrow(.)/nrow(input_grid)),
         meanIncubate = rep(input_grid$meanIncs, each = nrow(.)/nrow(input_grid))) -> gridOutputs
## Reformat for plotting
gridOutputs %>%
  mutate(dFeverMin = 0, dFeverMax = dFeverMin + caught.dpt.fever,
         dRiskMin = dFeverMax, dRiskMax = dRiskMin + caught.dpt.risk,
         aFeverMin = dRiskMax, aFeverMax = aFeverMin + caught.arv.fever,
         aRiskMin = aFeverMax, aRiskMax = aRiskMin +caught.arv.risk,
         mbMin = aRiskMax, mbMax = mbMin+missed.both,
         mfMin = mbMax, mfMax = mfMin+missed.fever.only,
         mrMin = mfMax, mrMax = mrMin+missed.risk.only,
         ndMin = mrMax, ndMax = ndMin+not.detectable,
         ndeMin = ndMax, ndeMax = ndeMin+evaded.screening) %>%
  select(days.since.exposed, fever, meanIncubate, contains('Min'), contains('Max')) %>%
  ## Pivot to long data frame
  pivot_longer(cols = dFeverMin:ndeMax, names_to = c('outcome', 'minOrMax'), names_pattern = '(\\w+)(M\\w\\w)', values_to = 'yy') -> temp
## Use full join to create columns for time, ymin, ymax, and band type
full_join(filter(temp, minOrMax == 'Min'), filter(temp, minOrMax == 'Max'), by = c('days.since.exposed', 'fever', 'meanIncubate', 'outcome'), suffix = c('min', 'max'))%>%
  select(-starts_with('min')) %>% 
  filter(outcome !='nde') %>%
  ## Clean up categorical variables so that plot labels are publication quality 
  mutate(fever = factor(fever, levels = rev(unique(fever)), labels = rev(paste0((1-unique(fever))*100,"% subclinical"))),  ## Rename levels for nice plotting
         meanIncubate = factor(meanIncubate, levels = unique(meanIncubate), labels = paste0('Mean incubation ', unique(meanIncubate), 'd')),
         outcome = factor(outcome, levels = rev(c('dFever', 'dRisk', 'aFever', 'aRisk', 'mb', 'mf', 'mr','nd')), 
                          labels =(rev(cat.labels))))-> rib
blackline <- filter(rib,outcome=="detected: arrival risk screen")
  ## Plot
  ggplot(rib)+
  geom_ribbon(aes(x = days.since.exposed, ymin = yymin, ymax = yymax, fill = outcome))+
  geom_line(data=blackline,aes(x=days.since.exposed,y=yymax),lty=2)+
  facet_grid(fever~meanIncubate) +
  scale_fill_manual(values = cols[8:1])+
  theme_bw() +
  ylab('Percentage of exposed individuals detected or missed')+
  scale_y_continuous(breaks = seq(0,1,.25),labels=paste(seq(0,100,25),"%",sep=""))+
  xlab('Days since exposure')  
ggsave('2020_nCov/Fig2S1_grid_of_ribbon_plots_departure_only.png', width = 8, height = 4.5, units = 'in')
write.csv(rib, file = '2020_nCov/Fig2S1_sourceData.csv', row.names = FALSE)



## --------------------------------- Fig 2. Supplemtary figure 2. Arrival screening only. ---------------------------------
gridWrapper = function(ff.in, mInc.in){
  incFun = get.inc.dist(meanInc = mInc.in)
  get_frac_caught_over_time(ff = ff.in, gg = 0.2, ascreen = TRUE, dscreen = FALSE, incubation.d = incFun, frac_evaded = 0)
}
apply(X = input_grid, MARGIN = 1, FUN = function(ii){gridWrapper(ii[1], ii[2])}) %>%                                                                
  bind_rows() %>% as.tbl() %>%
  mutate(fever = rep(input_grid$ffs, each = nrow(.)/nrow(input_grid)),
         meanIncubate = rep(input_grid$meanIncs, each = nrow(.)/nrow(input_grid))) -> gridOutputs
## Reformat for plotting
gridOutputs %>%
  mutate(dFeverMin = 0, dFeverMax = dFeverMin + caught.dpt.fever,
         dRiskMin = dFeverMax, dRiskMax = dRiskMin + caught.dpt.risk,
         aFeverMin = dRiskMax, aFeverMax = aFeverMin + caught.arv.fever,
         aRiskMin = aFeverMax, aRiskMax = aRiskMin +caught.arv.risk,
         mbMin = aRiskMax, mbMax = mbMin+missed.both,
         mfMin = mbMax, mfMax = mfMin+missed.fever.only,
         mrMin = mfMax, mrMax = mrMin+missed.risk.only,
         ndMin = mrMax, ndMax = ndMin+not.detectable,
         ndeMin = ndMax, ndeMax = ndeMin+evaded.screening) %>%
  select(days.since.exposed, fever, meanIncubate, contains('Min'), contains('Max')) %>%
  ## Pivot to long data frame
  pivot_longer(cols = dFeverMin:ndeMax, names_to = c('outcome', 'minOrMax'), names_pattern = '(\\w+)(M\\w\\w)', values_to = 'yy') -> temp
## Use full join to create columns for time, ymin, ymax, and band type
full_join(filter(temp, minOrMax == 'Min'), filter(temp, minOrMax == 'Max'), by = c('days.since.exposed', 'fever', 'meanIncubate', 'outcome'), suffix = c('min', 'max'))%>%
  select(-starts_with('min')) %>% 
  filter(outcome !='nde') %>%
  ## Clean up categorical variables so that plot labels are publication quality 
  mutate(fever = factor(fever, levels = rev(unique(fever)), labels = rev(paste0((1-unique(fever))*100,"% symptomatic"))),  ## Rename levels for nice plotting
         meanIncubate = factor(meanIncubate, levels = unique(meanIncubate), labels = paste0('Mean incubation ', unique(meanIncubate), 'd')),
         outcome = factor(outcome, levels = rev(c('dFever', 'dRisk', 'aFever', 'aRisk', 'mb', 'mf', 'mr','nd')), 
                          labels =(rev(cat.labels)))) -> rib
blackline <- filter(rib,outcome=="detected: arrival risk screen")
  ## Plot
  ggplot(rib)+
  geom_ribbon(aes(x = days.since.exposed, ymin = yymin, ymax = yymax, fill = outcome))+
  geom_line(data=blackline,aes(x=days.since.exposed,y=yymax),lty=2)+
  facet_grid(fever~meanIncubate) +
  scale_fill_manual(values = cols[8:1])+
  theme_bw() +
  ylab('Percentage of exposed individuals detected or missed')+
  scale_y_continuous(breaks = seq(0,1,.25),labels=paste(seq(0,100,25),"%",sep=""))+
  xlab('Days since exposure')  
ggsave('2020_nCov/Fig2S2_grid_of_ribbon_plots_arrival_only.png', width = 8, height = 4.5, units = 'in')
write.csv(rib, file = '2020_nCov/Fig2S2_sourceData.csv', row.names = FALSE)






## -------------------------------------------------------------------------
## \\\\\\\\\\\\\\       Fig. 3 & supplementary figs    ////////////////// ##
##                        Growing epidemic
## -------------------------------------------------------------------------
## Generate a range of par combos to test
## Use Latin Hypercube Sampling to span plausible parameter ranges
parsets = pomp::sobolDesign(lower = low.vals, upper = high.vals, nseq = nsim)  # sobolDesign from package pomp draws LHS samples
## Replicate the list of parsets across each subclinical case fraction tested
parsets = bind_rows(parsets, parsets, parsets) %>% mutate(ff = rep(NOTsubclinical, each = nsim))


## --------------------------------- Fig. 3 Analyses  ---------------------------------

## Simulate and save population outcomes
## This takes about 10 mins to run. Could be parallelized easliy.
## Get outcomes for both arrival and departure
cl = makeCluster(detectCores()-1)

if(!file.exists('bootList_ad.RData')|reset){
  bootWrapper = function(f.in, g.in, f.sens, g.sens, mInc, r0, mToAdmit){ one_sim(meanInc = mInc, R0 = r0, f0 = f.in, g0 = g.in, f.sens, g.sens, del.d=flight.time, as=TRUE, ds=TRUE, meanToAdmit = mToAdmit)}
  ## Simulate one population for each plausible paramter set
  mapply(FUN = bootWrapper,
         f.in = parsets$ff,
         g.in = parsets$gg,
         f.sens = parsets$f.sens,
         g.sens = parsets$g.sens,
         mInc = parsets$mInc,
         r0 = parsets$R0, 
         mToAdmit = parsets$meanToAdmit) -> bootList_ad
  save(bootList_ad, file = 'bootList_ad.RData')
}else{
  load('bootList_ad.RData')
}
## Get outcomes for departure only
if(!file.exists('bootList_d.RData')|reset){
  bootWrapper = function(f.in, g.in, f.sens, g.sens, mInc, r0, mToAdmit){ one_sim(meanInc = mInc, R0 = r0, f0 = f.in, g0 = g.in, f.sens, g.sens, del.d=flight.time, as=FALSE, ds=TRUE, meanToAdmit = mToAdmit)}
  ## Simulate one population for each plausible paramter set
  mapply(FUN = bootWrapper,
         f.in = parsets$ff,
         g.in = parsets$gg,
         f.sens = parsets$f.sens,
         g.sens = parsets$g.sens,
         mInc = parsets$mInc, 
         r0 = parsets$R0,
         mToAdmit = parsets$meanToAdmit) -> bootList_d
  save(bootList_d, file = 'bootList_d.RData')
}else{
  load('bootList_d.RData')
}
## Get outcomes for arrival only
if(!file.exists('bootList_a.RData')|reset){
  bootWrapper = function(f.in, g.in, f.sens, g.sens, mInc, r0, mToAdmit){ one_sim(meanInc = mInc, R0 = r0, f0 = f.in, g0 = g.in, f.sens, g.sens, del.d=flight.time, as=TRUE, ds=FALSE, meanToAdmit = mToAdmit)}
  ## Simulate one population for each plausible paramter set
  mapply(FUN = bootWrapper,
         f.in = parsets$ff,
         g.in = parsets$gg,
         f.sens = parsets$f.sens,
         g.sens = parsets$g.sens,
         mInc = parsets$mInc,
         r0 = parsets$R0,
         mToAdmit = parsets$meanToAdmit) -> bootList_a
  save(bootList_a, file = 'bootList_a.RData')
}else{
  load('bootList_a.RData')
}




## --------------------------------- Fig. 3 Plots  ---------------------------------
# -------------------------------
### Fraction Caught (Fig. 3A)
# -------------------------------
data.frame(departure.only = sapply(bootList_d[2,], function(yy){yy}),
           arrival.only = sapply(bootList_a[2,], function(yy){yy}),
           both = sapply(bootList_ad[2,], function(yy){yy})) %>%
  mutate(ff = parsets$ff) %>%
  pivot_longer(1:3, names_to = c('screen_type'), values_to = 'frac.missed') %>%
  group_by(screen_type, ff)  %>% 
  mutate(sc_type = factor(screen_type, levels = c('departure.only', 'arrival.only', 'both'),
                          labels = c('departure', 'arrival', 'both')),
         scenario = factor(ff, levels = rev(NOTsubclinical), 
                           labels = c('5% subclinical', '25% subclinical', '50% subclinical'))) -> long_frac_caught

long_frac_caught %>% summarise(med = median(1-frac.missed),
                   lower = quantile(1-frac.missed, probs = .025),
                   upper = quantile(1-frac.missed, probs = .975)) %>%
  ungroup() %>%
  mutate(sc_type = factor(screen_type, levels = c('departure.only', 'arrival.only', 'both'),
                          labels = c('departure', 'arrival', 'both')),
         scenario = factor(ff, levels = rev(NOTsubclinical), 
                           labels = c('5% subclinical', '25% subclinical',
                                      '50% subclinical'))) -> frac  


ggplot(long_frac_caught)+
  geom_violin(aes(x=sc_type,y=1-frac.missed))+
  geom_point(data=frac,aes(x = sc_type, y = med), size = 3)+
  geom_segment(data=frac,aes(x = sc_type, xend = sc_type, y = lower, yend = upper))+
  geom_text(data = frac, aes(x = sc_type, y = upper+.4, label = sprintf('%1.2f', med))) +
  geom_text(data = frac, aes(x = sc_type, y = upper+.3, label = sprintf('(%1.2f-%1.2f)', lower, upper)), size = 3.5) +
  theme_bw()+
  xlab('Screening type')+
  ylab('Fraction detected')+
  ylim(c(0,1))+
  #geom_dotplot(binaxis='y', binwidth = .005,stackdir='center',aes(x=sc_type,y=1-frac.missed),dotsize=.5,alpha=.5)+
  facet_wrap(~scenario)   -> fracCaught
fracCaught
write.csv(long_frac_caught, file = '2020_nCov/Fig3A_sourceData.csv', row.names = FALSE)


# -------------------------------
## Stacked barplot (Fig. 3B)
# -------------------------------
bind_rows(
  departureMeans = (sapply(bootList_d[1,], function(yy){yy}) %>% t() %>% as.data.frame()),
  arrivalMeans = (sapply(bootList_a[1,], function(yy){yy}) %>% t() %>% as.data.frame()),
  bothMeans = (sapply(bootList_ad[1,], function(yy){yy}) %>% t() %>% as.data.frame()) 
) %>%
  mutate(scenario = factor(rep(c('50% subclinical', '25% subclinical', '5% subclinical'), each = nsim) %>% rep(times = 3), 
                           levels = rev(c('50% subclinical', '25% subclinical', '5% subclinical')))) %>%
  mutate(strategy = rep(c('departure', 'arrival', 'both'), each = nsim*3)) %>%
  group_by(scenario, strategy) %>%
  summarise_all(mean) %>% ungroup() -> meanOutcomes
names(meanOutcomes) = c('scenario', 'strategy', 'd.fever', 'd.risk', 'a.fever', 'a.risk', 'm.b', 'm.f', 'm.r', 'nd', 'nde')

meanOutcomes %>%
  select(-nde) %>%
  pivot_longer(cols = 3:10) %>%
  mutate(outcome = factor(name, levels =rev(c('d.fever', 'd.risk', 'a.fever', 'a.risk', 'm.b', 'm.f', 'm.r', 'nd')), 
                          labels =rev(cat.labels)),
         strategy = factor(strategy, levels = c('departure', 'arrival', 'both'), labels = c('departure', 'arrival', 'both'))) -> stackedB

dashedLine = group_by(stackedB, strategy, scenario) %>%
  filter(name %in% c('d.fever', 'd.risk', 'a.fever', 'a.risk')) %>%
  summarise(yy = sum(value))

  ggplot(stackedB)+
  geom_bar(aes(x = strategy, y = value, fill = outcome), stat = 'identity')+
  geom_segment(data = dashedLine, aes(x = as.numeric(strategy)-.5, xend = as.numeric(strategy)+.5, y = yy, yend = yy), linetype = 2)+
  scale_fill_manual(values = rev(cols)) +
  xlab('Screening type')+
  ylab('Fraction travellers') +
  facet_grid(.~scenario)+
  theme_bw() +
  guides(fill=guide_legend(nrow = 4))+
  theme(legend.position = 'bottom') -> stackedBars
stackedBars
write.csv(stackedB, file = '2020_nCov/Fig3B_sourceData.csv', row.names = FALSE)

# -------------------------------
## n detected before one missed (Fig. 3C)
# -------------------------------
bind_rows(
  departureMeans = (sapply(bootList_d[4,], function(yy){yy}) %>% t() %>% as.data.frame()),
  arrivalMeans = (sapply(bootList_a[4,], function(yy){yy}) %>% t() %>% as.data.frame()),
  bothMeans = (sapply(bootList_ad[4,], function(yy){yy}) %>% t() %>% as.data.frame()) 
) %>%
  mutate(scenario = factor(rep(c('50% subclinical', '25% subclinical', '5% subclinical'), each = nsim) %>% rep(times = 3), 
                           levels = rev(c('50% subclinical', '25% subclinical', '5% subclinical')))) %>%
  mutate(strategy = rep(c('departure', 'arrival', 'both'), each = nsim*3)) %>%
  pivot_longer(V1:V30, names_to = 'nthTraveller', values_to = 'outcome') %>%
  group_by(strategy, scenario, nthTraveller) %>%
  summarise(frac = sum(outcome)/n()) %>%
  extract(col = nthTraveller, into = 'nthTraveller', regex = 'V(\\d+)') %>%
  mutate(nthTraveller = as.numeric(nthTraveller)) %>%
  ungroup() %>%
  mutate(strategy = factor(strategy, levels = c('both', 'arrival', 'departure'))) -> nBeforeMiss

xmax = nBeforeMiss %>% filter(frac>=0.05) %>% pull(nthTraveller) %>% max() + 2
lbls = as.character(1:xmax); lbls[xmax] = paste0(xmax, '+', collapse = '')
nBeforeMiss %>%
  ungroup() %>% rowwise() %>%
  ## Group all nthTraveller rows above the 95th percentile into the same bin
  mutate(nthTraveller = ifelse(nthTraveller > xmax, xmax, nthTraveller)) %>%
  ungroup() %>% group_by(strategy, scenario, nthTraveller) %>%
  summarise(frac = sum(frac)) %>%
  ungroup() %>% 
  mutate(nthTraveller = factor(nthTraveller, labels = lbls)) -> truncated
# Get text labels
txtlbs = sprintf('%1.2f', truncated$frac); txtlbs = sapply(txtlbs, function(tt){if(tt == '0.00'){'<0.01'} else{ tt} })
# Plot
ggplot(truncated)+
  geom_bar(aes(x = nthTraveller, y = frac, color = strategy, fill = strategy), stat = 'identity', alpha = .5, position = 'dodge')+
  geom_text(aes(x = as.numeric(nthTraveller), y = (-.038)+(-(as.numeric(strategy)-2)*.025), label = txtlbs, color = strategy), hjust = .5, size = 3, show.legend = FALSE) +
  facet_grid(.~scenario) +
  scale_fill_viridis_d(end = .8)+
  scale_color_viridis_d(end = .8)+
  theme_bw() +
  xlab('At least n detected before the first case importation') +
  ylab('Fraction of simulations')+
  ylim(-.07, .35)+
  theme(legend.position = 'bottom') -> nMissFig
nMissFig
write.csv(truncated, file = '2020_nCov/Fig3C_sourceData.csv', row.names = FALSE)


## Summarise n detected before first miss
tibble(
  departure = sapply(bootList_d[4,], function(yy){yy %>% which.min()-1}),
  arrival = sapply(bootList_a[4,], function(yy){yy %>% which.min()-1}),
  both = sapply(bootList_ad[4,], function(yy){yy %>% which.min()-1})
) %>%
  mutate(scenario = factor(rep(c('50% subclinical', '25% subclinical', '5% subclinical'), each = nsim),
                           levels = rev(c('50% subclinical', '25% subclinical', '5% subclinical')))) %>%
  pivot_longer(departure:both, names_to = 'strategy', values_to = 'nBeforeMiss') %>%
  group_by(scenario, strategy) %>%
  summarise(median = median(nBeforeMiss),
            upper = quantile(x = nBeforeMiss, .95))



## Layout fig. 3 and save
png('2020_nCov/Fig3_populationOutcomes.png', width = 7, height = 10, units = 'in', res = 480)
grid.arrange(fracCaught, stackedBars, nMissFig, nrow = 3, heights = c(2,3,2.5))
dev.off()





## --------------------------------- Plot Fig. 3S2 ---------------------------------
## Plot plausible incubation period distributions 
#shape=k scale= theta
#scale = variance/mean
#shape = mean/scale
xx = seq(0, 25, by = .01)
dat <- NA
for (m in incubationMeans){
  theta = incVariance/m
  k = m/theta       
  dat <- c(dat,dgamma(xx,shape = k,scale = theta))
}
dat_all <- data.frame(xval=rep(xx,length(means)),vals=dat[-1],
                      meangamma=as.factor(rep(means,each=length(xx))))
png(filename = '2020_nCov/Fig3S2_parRanges.png', width = 6, height = 4, units = 'in', res = 480)
ggplot(dat_all)+
  geom_line(aes(x=xval,y=vals,color=meangamma),lwd=1)+
  scale_color_viridis(discrete=T,name="mean (days)")+
  scale_x_continuous(limits = c(0,25))+
  xlab("incubation period (days)")+
  theme_classic(base_size=14)+
  ylab("density")
dev.off()
write.csv(dat_all, file = '2020_nCov/Fig3S2_sourceData.csv', row.names = FALSE)






# -------------------------------
## PRCC analysis (Fig. 4)
# -------------------------------
## Write a function to perform PRCC on a given set of inputs and outputs
get_prcc = function(input_df, output_df){
  ## input_df is a data frame whose columns contain input parameter values
  ## output_df is a data frame with one column, and the same number of rows as input_df
  ## PRCC method is described in https://reader.elsevier.com/reader/sd/pii/S0022519308001896?token=63EDCF653C192107235289FFA257089BE0BD8CFEF8501D34689D22AEAA56F88094A7B85F746CDD7BE96BC0828FBC926D
  
  rank_inputs <- mutate_all(input_df, rank)
  rank_output <- mutate_all(output_df, rank)

get_one_PRCC = function(ii){    
    ## Regression 1
    out1 = names(rank_inputs)[ii]
    ins = names(rank_inputs)[-ii]
    form1 = sprintf("%s ~ %s", out1, (paste0(ins, collapse = '+')))
    reg1 = lm(form1, data = rank_inputs)
    resid1 = reg1$residuals
    
    
    ## Regression 2
    out2 = names(rank_output)
    dd = cbind(rank_inputs[,-ii], rank_output)
    form2 = sprintf("%s ~ %s", out2, (paste0(ins, collapse = '+')))
    reg2 = lm(form2, data = dd)
    resid2 = reg2$residuals
    
    PRCC = cor.test(resid1, resid2, method = 'pearson')
    list(par = out1, coef = PRCC$estimate, pval = PRCC$p.value)
}
sapply(1:ncol(input_df), get_one_PRCC)
}


## Perform PRCC for each of the following combinations of conditions:
### a. Fraction subclinical 
### b. Screening type
### c. Output = fraction detected or p(first 3 stopped)

## Get a data set of all relevant inputs and outputs
full_join(  ## Join outcomes with underlying parameter values
  tibble(departure.only = sapply(bootList_d[2,], function(yy){yy}),
         arrival.only = sapply(bootList_a[2,], function(yy){yy}),
         both = sapply(bootList_ad[2,], function(yy){yy})) %>%
    mutate(parset = 1:nrow(parsets)),
  parsets %>% mutate(parset = 1:nrow(.)),
  by = 'parset') %>%
  pivot_longer(cols = departure.only:both, names_to = 'screening_type', values_to = 'frac_missed') %>%
  mutate(frac_detected = 1-frac_missed,
         fracSubclinical = factor(1-ff,
                                     labels = c('5%', '25%', '50%')),
        # parameter = factor(parameter, 
        #                    labels = c('Fever sensitivity', 'Risk sensitivity', 'Fraction aware of risk', 'Mean onset to admit (d)', 'Mean incubation (d)', 'R0')),
         screening_type = factor(screening_type,
                                 labels = c('arrival', 'departure', 'both'))) -> fDetainedPRCC_wide

## Apply PRCC function to all combinations
mapply(FUN = function(st, fsc){
  inDf = filter(fDetainedPRCC_wide, screening_type == st & fracSubclinical == fsc) %>%
    select(gg:meanToAdmit)
  outDf = filter(fDetainedPRCC_wide, screening_type == st & fracSubclinical == fsc) %>%
    select(frac_detected)
  list(get_prcc(input_df = inDf, output_df = outDf) %>% t() %>% as.data.frame() %>% mutate(NOTsubclinical = fsc, screeningType = st))
}, 
st = rep(c('arrival', 'departure', 'both'), each = 3),
fsc = rep(c('50%', '25%', '5%'), times = 3),
SIMPLIFY = 'list') %>%
  bind_rows() %>%
  mutate(screeningType = factor(screeningType, levels = c('both', 'arrival', 'departure')),
         par = unlist(par) %>% as.factor(),
         pval = unlist(pval),
         coef = unlist(coef)) -> allPRCC

allPRCC %>%
  filter(NOTsubclinical == '25%') %>%
  mutate(is.signif.wBonferroni = pval < 0.05/nrow(allPRCC),
         ctxt = sprintf('%1.2f%s', coef, ifelse(is.signif.wBonferroni, '*', '')),
         par = factor(par, levels = c('R0', 'mInc', 'gg', 'f.sens', 'g.sens', 'meanToAdmit'),
                      labels = c('R0', 'mean incubation\n(days)', 'frac aware\nof risk', 'thermal scanner\nsensitivity', 'questionnaire\nsensitivity', 'onset to\nisolation (days)')),
         NOTsubclinical = factor(NOTsubclinical, levels = c('50%', '25%', '5%'), labels = c('50% subclinical', '25% subclinical', '5% subclinical'))) -> mainPRCC

mainPRCC %>%
  ggplot() +
  geom_bar(aes(x = par, y = coef, fill = screeningType), stat = 'identity', position = 'dodge', color = 'black', alpha = .5) +
  scale_fill_viridis_d(end = .8, name = 'Screening type')+
  scale_color_viridis_d(end = .8)+
  geom_text(aes(x = (as.numeric(par)+((as.numeric(screeningType)-2)*.35)), y = coef*1.4, label = ctxt), size = 2.5)+
  facet_grid(NOTsubclinical~.)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=0, hjust = .5)) +
  xlab('Parameter')+
  ylab('Partial rank correlation coefficient') +
  theme(legend.position = 'top')
ggsave('2020_nCov/Fig4_PRCC.png', width = 7, height = 3.5, units = 'in')


allPRCC %>%
  mutate(is.signif.wBonferroni = pval < 0.05/nrow(allPRCC),
         ctxt = sprintf('%1.2f%s', coef, ifelse(is.signif.wBonferroni, '*', '')),
         par = factor(par, levels = c('R0', 'mInc', 'gg', 'f.sens', 'g.sens', 'meanToAdmit'),
                      labels = c('R0', 'mean incubation\n(days)', 'frac aware\nof risk', 'thermal scanner\nsensitivity', 'questionnaire\nsensitivity', 'onset to\nisolation (days)')),
         NOTsubclinical = factor(NOTsubclinical, levels = c('50%', '25%', '5%'), labels = c('50% subclinical', '25% subclinical', '5% subclinical'))) -> suppPRCC

suppPRCC %>%
  ggplot() +
  geom_bar(aes(x = par, y = coef, fill = screeningType), stat = 'identity', position = 'dodge', color = 'black', alpha = .5) +
  scale_fill_viridis_d(end = .8, name = 'Screening type')+
  scale_color_viridis_d(end = .8)+
  geom_text(aes(x = (as.numeric(par)+((as.numeric(screeningType)-2)*.35)), y = coef*1.4, label = ctxt), size = 2.5)+
  facet_grid(NOTsubclinical~.)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=0, hjust = .5)) +
  xlab('Parameter')+
  ylab('Partial rank correlation coefficient') +
  theme(legend.position = 'top')
ggsave('2020_nCov/Fig4S1_PRCC.png', width = 7, height = 9, units = 'in')
write.csv(suppPRCC, file = '2020_nCov/Fig4_4S1_sourceData.csv', row.names = FALSE)
















## -------------------------------------------------------------------------
## \\\\\\\\\\\\\\       Fig. 3 supplement    ////////////////// ##
##                       Flat epidemic
## -------------------------------------------------------------------------


## --------------------------------- Fig. 3S (flat epidemic) Analyses  ---------------------------------

## Simulate and save population outcomes
## This takes about 10 mins to run. Could be parallelized easliy.
## Get outcomes for both arrival and departure
cl = makeCluster(detectCores()-1)

if(!file.exists('bootList_ad_flat.RData')|reset){
  bootWrapper = function(f.in, g.in, f.sens, g.sens, mInc, r0, mToAdmit){ one_sim(meanInc = mInc, R0 = r0, f0 = f.in, g0 = g.in, f.sens, g.sens, del.d=flight.time, as=TRUE, ds=TRUE, meanToAdmit = mToAdmit, growing = FALSE)}
  ## Simulate one population for each plausible paramter set
  mapply(FUN = bootWrapper,
         f.in = parsets$ff,
         g.in = parsets$gg,
         f.sens = parsets$f.sens,
         g.sens = parsets$g.sens,
         mInc = parsets$mInc,
         r0 = parsets$R0, 
         mToAdmit = parsets$meanToAdmit) -> bootList_ad_flat
  save(bootList_ad_flat, file = 'bootList_ad_flat.RData')
}else{
  load('bootList_ad_flat.RData')
}
## Get outcomes for departure only
if(!file.exists('bootList_d_flat.RData')|reset){
  bootWrapper = function(f.in, g.in, f.sens, g.sens, mInc, r0, mToAdmit){ one_sim(meanInc = mInc, R0 = r0, f0 = f.in, g0 = g.in, f.sens, g.sens, del.d=flight.time, as=FALSE, ds=TRUE, meanToAdmit = mToAdmit, growing = FALSE)}
  ## Simulate one population for each plausible paramter set
  mapply(FUN = bootWrapper,
         f.in = parsets$ff,
         g.in = parsets$gg,
         f.sens = parsets$f.sens,
         g.sens = parsets$g.sens,
         mInc = parsets$mInc, 
         r0 = parsets$R0,
         mToAdmit = parsets$meanToAdmit) -> bootList_d_flat
  save(bootList_d_flat, file = 'bootList_d_flat.RData')
}else{
  load('bootList_d_flat.RData')
}
## Get outcomes for arrival only
if(!file.exists('bootList_a_flat.RData')|reset){
  bootWrapper = function(f.in, g.in, f.sens, g.sens, mInc, r0, mToAdmit){ one_sim(meanInc = mInc, R0 = r0, f0 = f.in, g0 = g.in, f.sens, g.sens, del.d=flight.time, as=TRUE, ds=FALSE, meanToAdmit = mToAdmit, growing = FALSE)}
  ## Simulate one population for each plausible paramter set
  mapply(FUN = bootWrapper,
         f.in = parsets$ff,
         g.in = parsets$gg,
         f.sens = parsets$f.sens,
         g.sens = parsets$g.sens,
         mInc = parsets$mInc,
         r0 = parsets$R0,
         mToAdmit = parsets$meanToAdmit) -> bootList_a_flat
  save(bootList_a_flat, file = 'bootList_a_flat.RData')
}else{
  load('bootList_a_flat.RData')
}



## --------------------------------- Fig. 3S (flat epidemic) Plots  ---------------------------------
# -------------------------------
### Fraction Caught (Fig. 3A)
# -------------------------------
data.frame(departure.only = sapply(bootList_d_flat[2,], function(yy){yy}),
           arrival.only = sapply(bootList_a_flat[2,], function(yy){yy}),
           both = sapply(bootList_ad_flat[2,], function(yy){yy})) %>%
  mutate(ff = parsets$ff) %>%
  pivot_longer(1:3, names_to = c('screen_type'), values_to = 'frac.missed') %>%
  group_by(screen_type, ff)  %>% 
  mutate(sc_type = factor(screen_type, levels = c('departure.only', 'arrival.only', 'both'),
                          labels = c('departure', 'arrival', 'both')),
         scenario = factor(ff, levels = rev(NOTsubclinical), 
                           labels = c('5% subclinical', '25% subclinical', '50% subclinical'))) -> long_frac_caught

long_frac_caught %>% summarise(med = median(1-frac.missed),
                               lower = quantile(1-frac.missed, probs = .025),
                               upper = quantile(1-frac.missed, probs = .975)) %>%
  ungroup() %>%
  mutate(sc_type = factor(screen_type, levels = c('departure.only', 'arrival.only', 'both'),
                          labels = c('departure', 'arrival', 'both')),
         scenario = factor(ff, levels = rev(NOTsubclinical), 
                           labels = c('5% subclinical', '25% subclinical',
                                      '50% subclinical'))) -> frac  


ggplot(long_frac_caught)+
  geom_violin(aes(x=sc_type,y=1-frac.missed))+
  geom_point(data=frac,aes(x = sc_type, y = med), size = 3)+
  geom_segment(data=frac,aes(x = sc_type, xend = sc_type, y = lower, yend = upper))+
  geom_text(data = frac, aes(x = sc_type, y = upper+.25, label = sprintf('%1.2f', med))) +
  geom_text(data = frac, aes(x = sc_type, y = upper+.15, label = sprintf('(%1.2f-%1.2f)', lower, upper)), size = 3.5) +
  theme_bw()+
  xlab('Screening type')+
  ylab('Fraction detected')+
  ylim(c(0,1))+
  #geom_dotplot(binaxis='y', binwidth = .005,stackdir='center',aes(x=sc_type,y=1-frac.missed),dotsize=.5,alpha=.5)+
  facet_wrap(~scenario)   -> fracCaught
fracCaught
write.csv(long_frac_caught, file = '2020_nCov/Fig3S1A_sourceData.csv', row.names = FALSE)

# -------------------------------
## Stacked barplot (Fig. 3B)
# -------------------------------
bind_rows(
  departureMeans = (sapply(bootList_d_flat[1,], function(yy){yy}) %>% t() %>% as.data.frame()),
  arrivalMeans = (sapply(bootList_a_flat[1,], function(yy){yy}) %>% t() %>% as.data.frame()),
  bothMeans = (sapply(bootList_ad_flat[1,], function(yy){yy}) %>% t() %>% as.data.frame()) 
) %>%
  mutate(scenario = factor(rep(c('50% subclinical', '25% subclinical', '5% subclinical'), each = nsim) %>% rep(times = 3), 
                           levels = rev(c('50% subclinical', '25% subclinical', '5% subclinical')))) %>%
  mutate(strategy = rep(c('departure', 'arrival', 'both'), each = nsim*3)) %>%
  group_by(scenario, strategy) %>%
  summarise_all(mean) %>% ungroup() -> meanOutcomes
names(meanOutcomes) = c('scenario', 'strategy', 'd.fever', 'd.risk', 'a.fever', 'a.risk', 'm.b', 'm.f', 'm.r', 'nd', 'nde')

meanOutcomes %>%
  select(-nde) %>%
  pivot_longer(cols = 3:10) %>%
  mutate(outcome = factor(name, levels =rev(c('d.fever', 'd.risk', 'a.fever', 'a.risk', 'm.b', 'm.f', 'm.r', 'nd')), 
                          labels =rev(cat.labels)),
         strategy = factor(strategy, levels = c('departure', 'arrival', 'both'), labels = c('departure', 'arrival', 'both'))) -> stackedB

dashedLine = group_by(stackedB, strategy, scenario) %>%
  filter(name %in% c('d.fever', 'd.risk', 'a.fever', 'a.risk')) %>%
  summarise(yy = sum(value))

ggplot(stackedB)+
  geom_bar(aes(x = strategy, y = value, fill = outcome), stat = 'identity')+
  geom_segment(data = dashedLine, aes(x = as.numeric(strategy)-.5, xend = as.numeric(strategy)+.5, y = yy, yend = yy), linetype = 2)+
  scale_fill_manual(values = rev(cols)) +
  xlab('Screening type')+
  ylab('Fraction travellers') +
  facet_grid(.~scenario)+
  theme_bw() +
  guides(fill=guide_legend(nrow = 4))+
  theme(legend.position = 'bottom') -> stackedBars
stackedBars
write.csv(stackedB, file = '2020_nCov/Fig3S1B_sourceData.csv', row.names = FALSE)

# -------------------------------
## n detected before one missed (Fig. 3C)
# -------------------------------
bind_rows(
  departureMeans = (sapply(bootList_d_flat[4,], function(yy){yy}) %>% t() %>% as.data.frame()),
  arrivalMeans = (sapply(bootList_a_flat[4,], function(yy){yy}) %>% t() %>% as.data.frame()),
  bothMeans = (sapply(bootList_ad_flat[4,], function(yy){yy}) %>% t() %>% as.data.frame()) 
) %>%
  mutate(scenario = factor(rep(c('50% subclinical', '25% subclinical', '5% subclinical'), each = nsim) %>% rep(times = 3), 
                           levels = rev(c('50% subclinical', '25% subclinical', '5% subclinical')))) %>%
  mutate(strategy = rep(c('departure', 'arrival', 'both'), each = nsim*3)) %>%
  pivot_longer(V1:V30, names_to = 'nthTraveller', values_to = 'outcome') %>%
  group_by(strategy, scenario, nthTraveller) %>%
  summarise(frac = sum(outcome)/n()) %>%
  extract(col = nthTraveller, into = 'nthTraveller', regex = 'V(\\d+)') %>%
  mutate(nthTraveller = as.numeric(nthTraveller)) %>%
  ungroup() %>%
  mutate(strategy = factor(strategy, levels = c('both', 'arrival', 'departure'))) -> nBeforeMiss_flat

xmax = nBeforeMiss_flat %>% filter(frac>=0.05) %>% pull(nthTraveller) %>% max() + 2
lbls = as.character(1:xmax); lbls[xmax] = paste0(xmax, '+', collapse = '')

nBeforeMiss_flat %>%
  ungroup() %>% rowwise() %>%
  ## Group all nthTraveller rows above the 95th percentile into the same bin
  mutate(nthTraveller = ifelse(nthTraveller > xmax, xmax, nthTraveller)) %>%
  ungroup() %>% group_by(strategy, scenario, nthTraveller) %>%
  summarise(frac = sum(frac)) %>%
  ungroup() %>% 
  mutate(nthTraveller = factor(nthTraveller, labels = lbls)) -> truncated
# Get text labels
txtlbs = sprintf('%1.2f', truncated$frac); txtlbs = sapply(txtlbs, function(tt){if(tt == '0.00'){'<0.01'} else{ tt} })
# Plot
ggplot(truncated)+
  geom_bar(aes(x = nthTraveller, y = frac, color = strategy, fill = strategy), stat = 'identity', alpha = .5, position = 'dodge')+
  geom_text(aes(x = as.numeric(nthTraveller), y = (-.06)+(-(as.numeric(strategy)-2)*.035), label = txtlbs, color = strategy), hjust = .5, size = 3, show.legend = FALSE) +
  facet_grid(.~scenario) +
  scale_fill_viridis_d(end = .8)+
  scale_color_viridis_d(end = .8)+
  theme_bw() +
  xlab('At least n detected before the first case importation') +
  ylab('Fraction of simulations')+
  ylim(-.1, .55)+
  theme(legend.position = 'bottom') -> nMissFig
nMissFig
write.csv(truncated, file = '2020_nCov/Fig3S1C_sourceData.csv', row.names = FALSE)


## Layout fig. 3 and save
png('2020_nCov/Fig3_populationOutcomes_flat.png', width = 7, height = 10, units = 'in', res = 480)
grid.arrange(fracCaught, stackedBars, nMissFig, nrow = 3, heights = c(2,3,2.5))
dev.off()







## -------------------------------------------------------------------------
## \\\\\\\\\\\\\\       Fig. 4 supplement    ////////////////// ##
##                       Flat epidemic
## -------------------------------------------------------------------------

## Perform PRCC for each of the following combinations of conditions:
### a. Fraction subclinical 
### b. Screening type
### c. Output = fraction detected or p(first 3 stopped)

## Get a data set of all relevant inputs and outputs
full_join(  ## Join outcomes with underlying parameter values
  tibble(departure.only = sapply(bootList_d_flat[2,], function(yy){yy}),
         arrival.only = sapply(bootList_a_flat[2,], function(yy){yy}),
         both = sapply(bootList_ad_flat[2,], function(yy){yy})) %>%
    mutate(parset = 1:nrow(parsets)),
  parsets %>% mutate(parset = 1:nrow(.)),
  by = 'parset') %>%
  pivot_longer(cols = departure.only:both, names_to = 'screening_type', values_to = 'frac_missed') %>%
  mutate(frac_detected = 1-frac_missed,
         fracSubclinical = factor(1-ff,
                                  labels = c('5%', '25%', '50%')),
         # parameter = factor(parameter, 
         #                    labels = c('Fever sensitivity', 'Risk sensitivity', 'Fraction aware of risk', 'Mean onset to admit (d)', 'Mean incubation (d)', 'R0')),
         screening_type = factor(screening_type,
                                 labels = c('arrival', 'departure', 'both'))) -> fDetainedPRCC_wide_flatEpidemic

## Apply PRCC function to all combinations
mapply(FUN = function(st, fsc){
  inDf = filter(fDetainedPRCC_wide_flatEpidemic, screening_type == st & fracSubclinical == fsc) %>%
    select(gg:meanToAdmit)
  outDf = filter(fDetainedPRCC_wide_flatEpidemic, screening_type == st & fracSubclinical == fsc) %>%
    select(frac_detected)
  list(get_prcc(input_df = inDf, output_df = outDf) %>% t() %>% as.data.frame() %>% mutate(NOTsubclinical = fsc, screeningType = st))
}, 
st = rep(c('arrival', 'departure', 'both'), each = 3),
fsc = rep(c('50%', '25%', '5%'), times = 3),
SIMPLIFY = 'list') %>%
  bind_rows() %>%
  mutate(screeningType = factor(screeningType, levels = c('both', 'arrival', 'departure')),
         par = unlist(par) %>% as.factor(),
         pval = unlist(pval),
         coef = unlist(coef)) -> allPRCC_flat

allPRCC_flat %>%
  mutate(is.signif.wBonferroni = pval < 0.05/nrow(allPRCC),
         ctxt = sprintf('%1.2f%s', coef, ifelse(is.signif.wBonferroni, '*', '')),
         par = factor(par, levels = c('R0', 'mInc', 'gg', 'f.sens', 'g.sens', 'meanToAdmit'),
                      labels = c('R0', 'mean incubation\n(days)', 'frac aware\nof risk', 'thermal scanner\nsensitivity', 'questionnaire\nsensitivity', 'onset to\nisolation (days)')),
         NOTsubclinical = factor(NOTsubclinical, levels = c('50%', '25%', '5%'), labels = c('50% subclinical', '25% subclinical', '5% subclinical'))) -> PRCC_flat

PRCC_flat %>%
  ggplot() +
  geom_bar(aes(x = par, y = coef, fill = screeningType), stat = 'identity', position = 'dodge', color = 'black', alpha = .5) +
  scale_fill_viridis_d(end = .8, name = 'Screening type')+
  scale_color_viridis_d(end = .8)+
  geom_text(aes(x = (as.numeric(par)+((as.numeric(screeningType)-2)*.35)), y = coef*1.4, label = ctxt), size = 2.5)+
  facet_grid(NOTsubclinical~.)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=0, hjust = .5)) +
  xlab('Parameter')+
  ylab('Partial rank correlation coefficient') +
  theme(legend.position = 'top')
ggsave('2020_nCov/Fig4S2_PRCC_flat.png', width = 7, height = 9, units = 'in')
write.csv(PRCC_flat, file = '2020_nCov/Fig4S2_sourceData.csv', row.names = FALSE)
