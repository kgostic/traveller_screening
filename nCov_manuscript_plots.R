library(ggplot2)
library(grid)
library(gridExtra)
library(tidyverse)
library(pomp)



## ------------------------------------------------------------
## \\\\\\\\\\\\\\ DEFINE GLOBAL VARIABLES ////////////////// ##
## ------------------------------------------------------------
## Set Global Vars
pathogen=c("nCoV")
pathtablab=c("2019-nCoV")
par(mfrow = c(length(pathtablab), 1))

#Specify efficacy parameters. These will be fed in to the external functions.
rd=0.25 #efficacy of departure questionnaire (proportion of travelers that report honestly)
ra=0.25 #efficacy of arrival questionnaire
sd=0.7 #efficacy of departure fever screening (based on fever detection sensitivity)
sa=0.7 #efficacy of arrival fever screening
nboot = 1000        # n sim samples
popn = 100        # population size of infected travelers
#Growing (0) or flat (1) epidemic
flatA=0
scale.in = 1.8
mToAdmit = 4.5 ## days from onset to hospitalization. (Assume people don't travel after admit)


cat.labels = c('detected: departure fever screen', "detected: departure risk screen", "detected: arrival fever screen",
               "detected: arrival risk screen", 'missed: had both', 'missed: had fever', 'missed: had risk awareness', 'missed: undetectable')




## -----------------------------------------------------
## \\\\\\\\\\\\\\ DEFINE FUNCTIONS ////////////////// ##
## -----------------------------------------------------

# --------------
# Distributions for fever and known exposure
# --------------
fg.distn<-function(type,pathogen){
  # Fever study sample sizes
  d1fn=c(74)
  # Proportions that display fever
  d1fv=c(0.77)
  
  ######## THIS NEEDS AN UPDATE! 
  # Exposure study sample sizes
  d1gn=c(1)
  # Proportions with known exposure
  d1gv=c(.1)
  
  if(type=="f"){
    c(sum(d1fv*d1fn)/sum(d1fn),sum(d1fn))}else{
      c(sum(d1gv*d1gn)/sum(d1gn),sum(d1gn))}
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


## --------------
exposure.outcome<-function(x,pathogen,r0 = 2,type=1,flat=0){
  #pathogen-specific median time to outcome (median exposure -> onset + onset -> outcome)
  if(pathogen == "H7N9"){f1<-3.1+4.2}
  if(pathogen == "MERS"){f1<-(5.2*7+5.5*23)/(7+23)+5}  ## (5.2*7+5.5*23)/(7+23) is a weighted average based on sample sizes of individual studies
  if(pathogen == "SARS"){f1<-6.4+4.85}
  if(pathogen == "nCoV"){f1<-4.5+5.5}
  median=f1
  # Stable epidemic
  if(flat==1){
    if(type==1){
      ifelse(x<=f1,1/(f1),0) #pdf
    }else{
      ifelse(x<=f1,x/(f1),1) #cdf
    }
  }else{
    # Growing epidemic
    if(type==1){
      pdf.cdf.travel(x,r0,median,type)
    }else{
      pdf.cdf.travel(x,r0,median,type)
    }
  }
}
## --------------


## --------------
exposure.distn<-function(x,pathogen,type=2,flat=0){
  # Use uniform random variable x to make random draw from specified infection age distribution
  a=0
  while(x>exposure.outcome(a,pathogen,type,flat)){a=a+0.1}
  a
}
## --------------



## --------------
screen.passengers = function(d, del.d, f, g, sd = 1, sa =1, rd = 1, ra = 1, 
                             phi.d, incubation.d, pathogen, relative=0, split1=0, arrival_screen, departure_screen, frac_evade=0){
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
  #   phi.d = call to function describing pdf of time from exposure to outcome
  #   incubation.d = call to function describing pdf of time from exposure to onset
  #   pathogen = name of the pathogen (e.g. "H7N9")
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
    p.new.symptoms = ((incubation.d(d+del.d)-incubation.d(d))/(1-incubation.d(d)))
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
#get_frac_caught = function(tSinceExposed, ff, gg, R0, meanToAdmit, meanIncubate = NULL, dscreen, ascreen, shapeIncubate = NULL, scaleIncubate = 2.73){
get_frac_caught = function(tSinceExposed, ff, gg, R0, meanToAdmit, dscreen, ascreen, incubation.d, frac_evaded){
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
  screen.passengers(tSinceExposed, del.d=1, ff, gg, sd, sa, rd, ra, phi.d, incubation.d, pathogen, relative = 0, split1 = 2, arrival_screen = ascreen, departure_screen = dscreen, frac_evaded)
}
# ## Test function
# get_frac_caught(tSinceExposed = 2, ff = .7, gg = .2, R0 = 2, meanToAdmit = 5.5, dscreen = TRUE, ascreen = FALSE, incubation.d)
# -------------------------------

# -------------------------------
## Write a function that repeats get_frac_caught over a grid of times since exposure
get_frac_caught_over_time = function(ff, gg, R0, meanToAdmit, ascreen, dscreen, incubation.d, frac_evaded = 0){
  arrive.times=seq(0,15,0.1)
  sapply(arrive.times, function(tt){get_frac_caught(tSinceExposed = tt, ff, gg, R0, meanToAdmit, dscreen, ascreen,incubation.d, frac_evaded)}) %>% t() %>% as.data.frame() %>% mutate(
    days.since.exposed = arrive.times)
}
# Test
# get_frac_caught_over_time(ff = .7, gg = .2, R0 = 2, meanToAdmit = 5.5, ascreen = TRUE, dscreen = FALSE, incubation.d, frac_evaded = .1)
# -------------------------------



# -------------------------------
# Simulate the fraction of the population caught or missed in a growing epidemic.
one_sim = function(meanInc, R0, f0, g0, f.sens, g.sens, gg, del.d, as, ds){
  ## For each individual in the population, draw times since exposure from the appropriate distribution
  infAge=sapply(c(1:popn),function(x){exposure.distn(runif(1, 0, 1),pathogen,flat=0)})
  
  this.inc.cdf = (function(x){pgamma(q = x, shape = meanInc/scale.in, scale = scale.in)})           ## Set incubation period distribution
  
  ## Fever and risk screening
  # outcomesBoth=sapply(infAge,function(x){
  #   f0=rbinom(1,fever.sample.size,fever.prob)/fever.sample.size                                 ## Draw binomial probability of fever
  #   g0=rbinom(1,risk.sample.size, risk.prob)/risk.sample.size                                   ## Draw binomial probability of known risk  
  #   screenWrapper(x, f0, g0)})           
  outcomesBoth=sapply(infAge, FUN = function(x){screen.passengers(x, del.d, f0, g0, f.sens, f.sens, g.sens, g.sens, 0, this.inc.cdf, pathogen, relative = 0, split1 = 2, arrival_screen=as, departure_screen=ds)})
  ## Output individual probability missed
  pCaught_both = colSums(outcomesBoth[1:4,])
  caughtBoth = sapply(pCaught_both,function(x){ifelse(x<runif(1, 0, 1),0,1)})                   ## Draw whether individual was missed
  
  # 
  # ## Fever only
  # # outcomesFever=sapply(infAge,function(x){
  # #   f0=rbinom(1,fever.sample.size,fever.prob)/fever.sample.size  ## Draw binomial probability of fever
  # #   g0=0  ## Draw binomial probability of known risk  
  # #   screenWrapper(x, f0, g0)})     
  # outcomesFever=sapply(infAge, FUN = function(x){screen.passengers(x, del.d, f0, 0, sd, sa, rd, ra, phi.d, this.inc.cdf, pathogen, relative = 0, split1 = 2, arrival_screen, departure_screen)})
  # pCaught_fever = colSums(outcomesFever[1:4,])                                                    ## Output individual probability missed
  # caughtFever = sapply(pCaught_fever,function(x){ifelse(x<runif(1, 0, 1),0,1)})                   ## Draw whether individual was missed
  # 
  # 
  # ## Risk only
  # # outcomesRisk=sapply(infAge,function(x){
  # #   f0=0  ## Draw binomial probability of fever
  # #   g0=rbinom(1,risk.sample.size, risk.prob)/risk.sample.size                                   ## Draw binomial probability of known risk  
  # #   screenWrapper(x, f0, g0)})                                                                  ## Output individual probability missed
  # outcomesRisk=sapply(infAge, FUN = function(x){screen.passengers(x, del.d, 0, g0, sd, sa, rd, ra, phi.d, this.inc.cdf, pathogen, relative = 0, split1 = 2, arrival_screen, departure_screen)})
  # pCaught_risk = colSums(outcomesRisk[1:4,])
  # caughtRisk = sapply(pCaught_risk,function(x){ifelse(x<runif(1, 0, 1),0,1)})                   ## Draw whether individual was missed
  # 
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
              # outcomesFever = rowMeans(outcomesFever),
              # outcomesRisk = rowMeans(outcomesRisk),
              caught = c(frac.missed.both = 1- sum(caughtBoth)/popn)))
}
# ## Test
# one_sim(meanInc = 5.5, R0 = 2,f0 = .7, g0 = .1, f.sens = .7, g.sens = .2, del.d = 1, as = TRUE, ds = FALSE)
# -------------------------------








## ------------------------------------------------------------
## \\\\\\\\\\\\\\       MAKE PLOTS     ////////////////// ##
## ------------------------------------------------------------
## Fig. 2. Plot the individual probability of different screening outcomes vs. times since exposure
##   Run across a grid of probabilities of detectable symptoms and mean incubation periods
##   Assume no one evades screening

input_grid = expand.grid(ffs = c(.5, .75, .98),                  ## Set grid of values to test
                         meanIncs = c(3, 5, 7))
gridWrapper = function(ff.in, mInc.in){
  incFun = function(x){pgamma(x, shape = mInc.in/scale.in, scale = scale.in)}
  get_frac_caught_over_time(ff = ff.in, gg = 0.2, R0 = 2, meanToAdmit = mToAdmit, ascreen = TRUE, dscreen = TRUE, incubation.d = incFun, frac_evaded = 0)
}


apply(X = input_grid, MARGIN = 1, FUN = function(ii){gridWrapper(ii[1], ii[2])}) %>%                                                                
  bind_rows() %>% as.tbl() %>%
  mutate(fever = rep(input_grid$ffs, each = nrow(.)/nrow(input_grid)),
         meanIncubate = rep(input_grid$meanIncs, each = nrow(.)/nrow(input_grid))) -> gridOutputs

## Reformat for plotting
cols = c('darkseagreen2', 'deepskyblue', 'seagreen4', 'royalblue3', 'bisque', 'brown4', 'salmon2', 'firebrick1')
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
  mutate(fever = factor(fever, levels = unique(fever), labels = paste0((1-unique(fever))*100,"% symptomatic")),  ## Rename levels for nice plotting
         meanIncubate = factor(meanIncubate, levels = unique(meanIncubate), labels = paste0('Mean incubation ', unique(meanIncubate), 'd')),
         outcome = factor(outcome, levels = rev(c('dFever', 'dRisk', 'aFever', 'aRisk', 'mb', 'mf', 'mr', 'nd')), 
                          labels =(rev(cat.labels)))) -> rib
  ## Plot
  ## Plot
  blackline <- filter(rib,outcome=="detected: arrival risk screen")
  ggplot(rib)+
  geom_ribbon(aes(x = days.since.exposed, ymin = yymin, ymax = yymax, fill = outcome))+
  facet_grid(fever~meanIncubate) +
  scale_fill_manual(values = cols[8:1])+
  theme_bw() +
  geom_line(data=blackline,aes(x=days.since.exposed,y=yymax),lty=2)+
  ylab('Percentage of exposed individuals detained or cleared')+
  scale_y_continuous(breaks = seq(0,1,.25),labels=paste(seq(0,100,25),"%",sep=""))+
  xlab('Days since exposure')  
ggsave('2020_nCov/Fig2_grid_of_ribbon_plots.png', width = 8, height = 4.5, units = 'in')




## Fig 2. Supplemtary figure 1. Departure screening only.
gridWrapper = function(ff.in, mInc.in){
  incFun = function(x){pgamma(x, shape = mInc.in/scale.in, scale = scale.in)}
  get_frac_caught_over_time(ff = ff.in, gg = 0.2, R0 = 2, meanToAdmit = mToAdmit, ascreen = FALSE, dscreen = TRUE, incubation.d = incFun, frac_evaded = 0)
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
  mutate(fever = factor(fever, levels = unique(fever), labels = paste0((1-unique(fever))*100,"% symptomatic")),  ## Rename levels for nice plotting
         meanIncubate = factor(meanIncubate, levels = unique(meanIncubate), labels = paste0('Mean incubation ', unique(meanIncubate), 'd')),
         outcome = factor(outcome, levels = rev(c('dFever', 'dRisk', 'aFever', 'aRisk', 'mb', 'mf', 'mr','nd')), 
                          labels =(rev(cat.labels)))) %>%
  ## Plot
  ## Plot
  ggplot()+
  geom_ribbon(aes(x = days.since.exposed, ymin = yymin, ymax = yymax, fill = outcome))+
  facet_grid(fever~meanIncubate) +
  scale_fill_manual(values = cols[8:1])+
  theme_bw() +
  ylab('Percentage of exposed individuals detained or cleared')+
  scale_y_continuous(breaks = seq(0,1,.25),labels=paste(seq(0,100,25),"%",sep=""))+

  xlab('Days since exposure')  
ggsave('2020_nCov/Fig2S1_grid_of_ribbon_plots_departure_only.png', width = 8, height = 4.5, units = 'in')




## Fig 2. Supplemtary figure 2. Arrival screening only.
gridWrapper = function(ff.in, mInc.in){
  incFun = function(x){pgamma(x, shape = mInc.in/scale.in, scale = scale.in)}
  get_frac_caught_over_time(ff = ff.in, gg = 0.2, R0 = 2, meanToAdmit = mToAdmit, ascreen = TRUE, dscreen = FALSE, incubation.d = incFun, frac_evaded = 0)
}
apply(X = input_grid, MARGIN = 1, FUN = function(ii){gridWrapper(ii[1], ii[2])}) %>%                                                                
  bind_rows() %>% as.tbl() %>%
  mutate(fever = rep(input_grid$ffs, each = nrow(.)/nrow(input_grid)),
         meanIncubate = rep(input_grid$meanIncs, each = nrow(.)/nrow(input_grid))) -> gridOutputs
## Reformat for plotting
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
  mutate(fever = factor(fever, levels = unique(fever), labels = paste0((1-unique(fever))*100,"% symptomatic")),  ## Rename levels for nice plotting
         meanIncubate = factor(meanIncubate, levels = unique(meanIncubate), labels = paste0('Mean incubation ', unique(meanIncubate), 'd')),
         outcome = factor(outcome, levels = rev(c('dFever', 'dRisk', 'aFever', 'aRisk', 'mb', 'mf', 'mr','nd')), 
                          labels =(rev(cat.labels)))) %>%
  ## Plot
  ## Plot
  ggplot()+
  geom_ribbon(aes(x = days.since.exposed, ymin = yymin, ymax = yymax, fill = outcome))+
  facet_grid(fever~meanIncubate) +
  scale_fill_manual(values = cols[8:1])+
  theme_bw() +
  ylab('Percentage of exposed individuals detained or cleared')+
  scale_y_continuous(breaks = seq(0,1,.25),labels=paste(seq(0,100,25),"%",sep=""))+
  xlab('Days since exposure')  
ggsave('2020_nCov/Fig2S2_grid_of_ribbon_plots_arrival_only.png', width = 8, height = 4.5, units = 'in')





## Generate a range of par combos to test
## Use Latin Hypercube Sampling to span plausible parameter ranges
low.vals = c(gg = .05, f.sens = .6, g.sens = .05, mInc = 3, R0 = 1.5)
high.vals = c(gg = .5, f.sens = .95, g.sens = .4, mInc = 7, R0 = 4)
parsets = sobolDesign(lower = low.vals, upper = high.vals, nseq = nboot)
parsets = bind_rows(parsets, parsets, parsets) %>% mutate(ff = rep(c(.5, .75, .95), each = nboot))


## Plot tested parameter ranges
data.frame(par = factor(names(low.vals)[-4], levels = names(low.vals[-4]), labels = c('Frac. aware of exposure', 'Sensitivy of symptom screen', 'Sensitivity of risk screen', 'R0')),
           lows = low.vals[-4],
           highs = high.vals[-4]) %>%
  mutate(labs = sprintf('%2.2f-%2.2f', lows, highs)) %>%
  ggplot()+
  geom_segment(aes(x = lows, xend = highs, y = par, yend = par), size = 2, color = 'darkblue') +
  geom_text(aes(x = (lows+highs)/2, y = as.numeric(par)+.2, label = labs, hjust = 0.5), size = 2.5)+
  theme_bw() +
  xlab('Assumed range')+
  ylab('Parameter') -> parRanges

## Plot plausible incubation period distributions
xx = seq(0, 25, by = .01)
sapply(X = seq(low.vals['mInc'], high.vals['mInc'], by = .5), FUN = function(mm){dgamma(xx, shape = mm/scale.in, scale = scale.in)}) -> gammaFits
colnames(gammaFits) = seq(low.vals['mInc'], high.vals['mInc'], by = .5)
as.data.frame(gammaFits) %>%
  mutate(x = xx) %>%
  melt(id.vars = 'x') %>%
  ggplot()+
  geom_line(aes(x = x, y = value, color = variable), size = .6, alpha = .5, show.legend = TRUE)+
  xlab('days')+
  ylab('density')+
  scale_color_discrete(name='mean days \nincubation')+
  theme_classic() -> incPeriods
  




## Simulate
reset = FALSE
## Get outcomes for both arrival and departure
if(!file.exists('bootList_ad.RData')|reset){
bootWrapper = function(f.in, g.in, f.sens, g.sens, mInc, r0){ one_sim(meanInc = mInc, R0 = r0, f0 = f.in, g0 = g.in, f.sens, g.sens, del.d=1, as=TRUE, ds=TRUE)}
  ## Simulate one population for each plausible paramter set
  mapply(FUN = bootWrapper,
         f.in = parsets$ff,
         g.in = parsets$gg,
         f.sens = parsets$f.sens,
         g.sens = parsets$g.sens,
         mInc = parsets$mInc,
         r0 = parsets$R0) -> bootList_ad
save(bootList_ad, file = 'bootList_ad.RData')
}else{
  load('bootList_ad.RData')
}
## Get outcomes for departure only
if(!file.exists('bootList_d.RData')|reset){
  bootWrapper = function(f.in, g.in, f.sens, g.sens, mInc, r0){ one_sim(meanInc = mInc, R0 = r0, f0 = f.in, g0 = g.in, f.sens, g.sens, del.d=1, as=FALSE, ds=TRUE)}
  ## Simulate one population for each plausible paramter set
  mapply(FUN = bootWrapper,
         f.in = parsets$ff,
         g.in = parsets$gg,
         f.sens = parsets$f.sens,
         g.sens = parsets$g.sens,
         mInc = parsets$mInc, 
         r0 = parsets$R0) -> bootList_d
  save(bootList_d, file = 'bootList_d.RData')
}else{
  load('bootList_d.RData')
}
## Get outcomes for arrival only
if(!file.exists('bootList_a.RData')|reset){
  bootWrapper = function(f.in, g.in, f.sens, g.sens, mInc, r0){ one_sim(meanInc = mInc, R0 = r0, f0 = f.in, g0 = g.in, f.sens, g.sens, del.d=1, as=TRUE, ds=FALSE)}
  ## Simulate one population for each plausible paramter set
  mapply(FUN = bootWrapper,
         f.in = parsets$ff,
         g.in = parsets$gg,
         f.sens = parsets$f.sens,
         g.sens = parsets$g.sens,
         mInc = parsets$mInc,
         r0 = parsets$R0) -> bootList_a
  save(bootList_a, file = 'bootList_a.RData')
}else{
  load('bootList_a.RData')
}

# -------------------------------
### Fraction missed
# -------------------------------
data.frame(departure.only = sapply(bootList_d[2,], function(yy){yy}),
           arrival.only = sapply(bootList_a[2,], function(yy){yy}),
           both = sapply(bootList_ad[2,], function(yy){yy})) %>%
  mutate(ff = parsets$ff) %>%
    pivot_longer(1:3, names_to = c('screen_type'), values_to = 'frac.missed') %>%
    group_by(screen_type, ff) %>%
    summarise(med = median(1-frac.missed),
              lower = quantile(1-frac.missed, probs = .025),
              upper = quantile(1-frac.missed, probs = .975)) %>%
  ungroup() %>%
  mutate(screen_type = factor(screen_type, levels = c('departure.only', 'arrival.only', 'both'),
                              labels = c('departure', 'arrival', 'both')),
         scenario = factor(ff, levels = c(.5, .75, .95), labels = c('50% subclinical', '25% subclinical', '5% subclinical'))) -> frac
    
    
  ggplot(frac)+
    geom_point(aes(x = screen_type, y = med), size = 3)+
    geom_segment(aes(x = screen_type, xend = screen_type, y = lower, yend = upper))+
    theme_bw()+
  facet_grid(.~scenario)+
    ylim(c(0,1))+
    xlab('Screening type')+
    ylab('Fraction caught')+
    theme(legend.position = "none") -> fracCaught
  fracCaught
  #ggsave('2020_nCov/Shiny_fraction_missed.pdf', width = 4, height = 4, units = 'in')  
  
  # -------------------------------
  ## Stacked barplot
  # -------------------------------
  cols = c('darkseagreen2', 'deepskyblue', 'seagreen4', 'royalblue3', 'bisque', 'brown4', 'salmon2', 'firebrick1')
  bind_rows(
    departureMeans = (sapply(bootList_d[1,], function(yy){yy}) %>% t() %>% as.data.frame()),
    arrivalMeans = (sapply(bootList_a[1,], function(yy){yy}) %>% t() %>% as.data.frame()),
    bothMeans = (sapply(bootList_ad[1,], function(yy){yy}) %>% t() %>% as.data.frame()) 
  ) %>%
    mutate(scenario = rep(c('50% subclinical', '25% subclinical', '5% subclinical'), each = nboot) %>% rep(times = 3)) %>%
    mutate(strategy = rep(c('departure', 'arrival', 'both'), each = nboot*3)) %>%
    group_by(scenario, strategy) %>%
    summarise_all(mean) %>% ungroup() -> meanOutcomes
  names(meanOutcomes) = c('scenario', 'strategy', 'd.fever', 'd.risk', 'a.fever', 'a.risk', 'm.b', 'm.f', 'm.r', 'nd', 'nde')
  
  meanOutcomes %>%
    select(-nde) %>%
    pivot_longer(cols = 3:10) %>%
    mutate(outcome = factor(name, levels =rev(c('d.fever', 'd.risk', 'a.fever', 'a.risk', 'm.b', 'm.f', 'm.r', 'nd')), 
                            labels =rev(cat.labels)),
           strategy = factor(strategy, levels = c('departure', 'arrival', 'both'), labels = c('departure', 'arrival', 'both'))) %>%
    ggplot()+
    geom_bar(aes(x = strategy, y = value, fill = outcome), stat = 'identity')+
    scale_fill_manual(values = rev(cols)) +
    xlab('Screening type')+
    ylab('Fraction') +
    facet_grid(.~scenario)+
    theme_bw() +
    guides(fill=guide_legend(nrow = 4))+
    theme(legend.position = 'bottom') -> stackedBars
  
  png(filename = '2020_nCov/Fig3S1_parRanges.png', width = 5.5, height = 7, units = 'in', res = 480)
  grid.arrange(parRanges, incPeriods, nrow = 2, ncol = 1)
  dev.off()
  
  png('2020_nCov/Fig3_populationOutcomes.png', width = 7, height = 7, units = 'in', res = 480)
  grid.arrange(fracCaught, stackedBars, nrow = 2, heights = c(2,3))
  dev.off()


#make_plots(meanIncubate = 5.5, meanToAdmit = 4.5, R0 = 2, ff = .7, gg = .2, flight.hrs = 24, screenType = 'both', nboot = 100, popn = 100)

  