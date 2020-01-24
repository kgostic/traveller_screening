
# - - - - - - - - - - - - - - - - - - - - - - -
# Simulate the fraction of the population caught or missed in a growing epidemic.
# - - - - - - - - - - - - - - - - - - - - - - -

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
# Get probability that symptoms have developed.
incubation.d<-dgamma(x = tSinceExposed, shape = shapeIncubate, scale = scaleIncubate) 


## \\\\\\\\\\\\\\\\\\ SIMULATE ///////////////////
one_bootstrap = function(){
## For each individual in the population, draw times since exposre from the appropriate distribution
infAge=sapply(c(1:popn),function(x){exposure.distn(runif(1, 0, 1),pathogen,flat=flatA)})

## Given times since exposure, screen using fever, risk or both
screenWrapper = function(x, ff, gg){screen.passengers(x, del.d, ff, gg, sd, sa, rd=.25 , ra=.25, phi.d, incubation.d,pathogen,0,0)}

## Fever and risk screening
pMissed_both=sapply(infAge,function(x){
  f0=rbinom(1,fg.distn("f",pathogen)[2],fg.distn("f",pathogen)[1])/fg.distn("f",pathogen)[2]  ## Draw binomial probability of fever
  g0=rbinom(1,fg.distn("g",pathogen)[2],fg.distn("g",pathogen)[1])/fg.distn("g",pathogen)[2]  ## Draw binomial probability of known risk  
  screenWrapper(x, f0, g0)})                                                                  ## Output individual probability missed
caughtBoth = sapply(pMissed_both,function(x){ifelse(x<runif(1, 0, 1),0,1)})                   ## Draw whether individual was missed

## Fever only
pMissed_fever=sapply(infAge,function(x){
  f0=rbinom(1,fg.distn("f",pathogen)[2],fg.distn("f",pathogen)[1])/fg.distn("f",pathogen)[2]  ## Draw binomial probability of fever
  g0=0  ## Draw binomial probability of known risk  
  screenWrapper(x, f0, g0)})                                                                  ## Output individual probability missed
caughtFever = sapply(pMissed_fever,function(x){ifelse(x<runif(1, 0, 1),0,1)})                   ## Draw whether individual was missed

## Risk only
pMissed_risk=sapply(infAge,function(x){
  f0=0  ## Draw binomial probability of fever
  g0=rbinom(1,fg.distn("g",pathogen)[2],fg.distn("g",pathogen)[1])/fg.distn("g",pathogen)[2]  ## Draw binomial probability of known risk  
  screenWrapper(x, f0, g0)})                                                                  ## Output individual probability missed
caughtRisk = sapply(pMissed_risk,function(x){ifelse(x<runif(1, 0, 1),0,1)})                   ## Draw whether individual was missed

## Return frac.missed.both, frac.missed.fever, frac.missed.risk
return(c(frac.missed.both = 1- sum(caughtBoth)/popn,
         frac.missed.fever = 1- sum(caughtFever)/popn,
         frac.missed.risk = 1- sum(caughtRisk)/popn))
}
## Test
# one_bootstrap()

replicate(nboot, expr = one_bootstrap()) %>% t() %>%
  as.data.frame() -> bootOutputs

bootOutputs %>%
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
  theme(legend.position = "none")
ggsave('2020_nCov/Fraction_missed.pdf', width = 4, height = 4, units = 'in')  


## Split into missed due to ...


