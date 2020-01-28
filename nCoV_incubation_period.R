## Estimate the incubation period distribution
#source('nCoV_Parmaeters.R')

## \\\\\\\\\\\\ INCUBATION PERIOD //////////////
outsideWuhan = read.csv('2020_nCov/outside_wuhan.csv', stringsAsFactors = FALSE) 
## Get full linelist and extract cases that meet our criteria:
##   Known time of onset
##   Known time of arrival and departure from source city
##   Not resident of source city
##   Not employed in source city
outsideWuhan %>% 
  select(age, sex, date_onset_symptoms, date_admission_hospital, date_arrive_Wuhan, date_depart_Wuhan, city, province, Home, source) %>%
  filter(date_arrive_Wuhan != "" & date_depart_Wuhan != "") %>%
  mutate_at(vars(starts_with('date')), to_date_outsideWuhan) %>%
  mutate(max_incubation = (date_onset_symptoms-date_arrive_Wuhan) %>% as.numeric(),
         min_incubation =   (date_onset_symptoms-date_depart_Wuhan) %>% as.numeric()) %>%
  rowwise() %>%
  mutate(exact = max_incubation == min_incubation) -> incubationPeriods
fixid = which(incubationPeriods$min_incubation < 0) # Clean up one case with symptom onset before departure 
incubationPeriods$min_incubation[fixid] = 0 
incubationPeriods$date_depart_Wuhan[fixid] = incubationPeriods$date_onset_symptoms[fixid]
incubationPeriods %>% ungroup() %>%
  mutate(mid_incubation = (min_incubation + max_incubation)/2 ,
         id = 1:nrow(.)) -> incubationPeriods
## Save valid entries
write.csv(x = incubationPeriods, file = '2020_nCov/known_exposures.csv', row.names = FALSE)


## Write functions to calculate MLEs
## Fit incubation period to distribution
nll_gamma = function(pars, dd, prof.par = 'none', prof.val = NULL){
  ## INPUTS - dd is the data (a vector of times elapsed),
  ##  pars is a named vector containing the shape (shp) and scale (scl) parameter of the gamma distribution
  ##  prof.val is an optional numeric input, which gives the value of the fixed parameter when profiling
  ## OUTPUTS - the negative log likelihood of the data given the parameters
  shp = ifelse(prof.par == 'shp', prof.val, pars['shp'])
  scl = ifelse(prof.par == 'scl', prof.val, pars['scl'])
  - (dgamma(dd, shape = shp, scale = scl, log = TRUE) %>% sum())
}

nll_weibull = function(pars, dd){
  ## INPUTS - dd is the data (a vector of times elapsed),
  ##  pars is a named vector containing the shape (shp) and scale (scl) parameter of the gamma distribution
  ##  prof.val is an optional numeric input, which gives the value of the fixed parameter when profiling
  ## OUTPUTS - the negative log likelihood of the data given the parameters
  shp = pars['shp']
  scl = pars['scl']
  - (dweibull(dd, shape = shp, scale = scl, log = TRUE) %>% sum())
}


nll_lognormal = function(pars, dd){
  mu = pars[1]
  sd = pars[2]
  ## INPUTS - dd is the data (a vector of times elapsed),
  ##  pars is a named vector containing the shape (shp) and scale (scl) parameter of the gamma distribution
  ##  prof.val is an optional numeric input, which gives the value of the fixed parameter when profiling
  ## OUTPUTS - the negative log likelihood of the data given the parameters
  - (dlnorm(dd, meanlog = log(mu), sdlog = log(sd), log = TRUE) %>% sum())
}

## Bootstrap to represent possible times of incubation
nBoot = 1000  ## Repeat 1000 times
nPerBoot = nrow(incubationPeriods)
set.seed(12)

## One_boot does the following:
##     1. Bootstrap individuals from the above line list
##     2. For each individual, randomly sample an exposure date from the range of dates spent in Wuhan or some other source city
##     3. Calculate individual incubation period as date of exposure (imputed) - date of onset (known)
##     4. Output MLE gamma and other parameters from steps 1-3
## Below, use function one_boot to repeat steps 1-4 1000 times. 
one_boot = function(){
incubationPeriods[sample(incubationPeriods$id, size = nPerBoot, replace = TRUE),] %>%  ## Sample individual rows. Times of onset should be known.
                    rowwise() %>%
                    mutate(date_exposed = if(exact == TRUE){
                      date_depart_Wuhan
                      }else{ ## If date of exposure is uncertain, sample from all known possible dates
                        date_arrive_Wuhan + sample(0:as.numeric(date_depart_Wuhan-date_arrive_Wuhan), size = 1)
                      },
                      boot_incubation = max(((date_onset_symptoms-date_exposed) %>% as.numeric()), .5)) -> iBoot
  boot_fits_raw = list(
    gamma = optim(par = c(shp = 1, scl = 2), fn = nll_gamma, dd = iBoot$boot_incubation),
    weibull = optim(par = c(shp = 1, scl = 2), fn = nll_weibull, dd = iBoot$boot_incubation),
    lognormal = optim(par = c(mu = 1, sd = 2), fn = nll_lognormal, dd = iBoot$boot_incubation)
  )

bootMLEs = data.frame(type = names(boot_fits_raw),
                      par1 = sapply(boot_fits_raw, function(xx) xx$par[1]),
                      par2 = sapply(boot_fits_raw, function(xx) xx$par[2]),
                      nll = sapply(boot_fits_raw, function(xx) xx$value)) %>%
  mutate(AIC = 2*2+2*nll) %>%
  arrange(AIC) %>%
  mutate(delAIC = AIC - min(AIC))
bootMLEs
}

## Repeat resampling procedure 1000 times
lapply(X = 1:nBoot, FUN = function(xx){one_boot()}) -> bootResults

## Extract gamma results
gamma_results = sapply(bootResults, function(ll) ll %>% as.data.frame() %>% filter(type == 'gamma'))
## Summarise gamma results to get the mean, and central 95% of shape and scale estimates
grs = data.frame(par1 = gamma_results['par1',] %>% unlist(), par2 = gamma_results['par2',] %>% unlist()) %>%
  mutate(mean = par1*par2) %>%
  summarise(shape.med = median(par1),
            shape.low = quantile(par1, .025),
            shape.high = quantile(par1, .975),
            scale.med = median(par2),
            scale.low = quantile(par2, .025),
            scale.high = quantile(par2, .975),
            meanInc.med = median(mean),
            meanInc.low = quantile(mean, .025),
            meanInc.high = quantile(mean, .975))
grs
sprintf('Gamma: The mean incubation period is %1.2f (%1.2f-%1.2f). MLE shape is %1.2f (%1.2f, %1.2f). MLE scale is %1.2f (%1.f, %1.2f)',
        grs[7], grs[8], grs[9], grs[1], grs[2], grs[3], grs[4], grs[5], grs[6])


## Extract lognormal results
lognormal_results = sapply(bootResults, function(ll) ll %>% as.data.frame() %>% filter(type == 'lognormal'))
## Summarise lognormal results to get the mean, and central 95% of shape and scale estimates
rs = data.frame(par1 = lognormal_results['par1',] %>% unlist(), par2 = lognormal_results['par2',] %>% unlist()) %>%
  mutate(mean = par1*par2) %>%
  summarise(mean.med = median(par1),
            mean.low = quantile(par1, .025),
            mean.high = quantile(par1, .975),
            sd.med = median(par2),
            sd.low = quantile(par2, .025),
            sd.high = quantile(par2, .975))
rs
sprintf('lognormal: The mean incubation period is %1.2f (%1.2f-%1.2f). Standard dev. is %1.2f (%1.2f, %1.2f).',
        rs[1], rs[2], rs[3], rs[4], rs[5], rs[6])


## Extract weibull results
weibull_results = sapply(bootResults, function(ll) ll %>% as.data.frame() %>% filter(type == 'weibull'))
## Summarise weibull results to get the mean, and central 95% of shape and scale estimates
rs = data.frame(par1 = weibull_results['par1',] %>% unlist(), par2 = weibull_results['par2',] %>% unlist()) %>%
  summarise(shape.med = median(par1),
            shape.low = quantile(par1, .025),
            shape.high = quantile(par1, .975),
            scale.med = median(par2),
            scale.low = quantile(par2, .025),
            scale.high = quantile(par2, .975))
rs
sprintf('weibull: MLE shape is %1.2f (%1.2f, %1.2f). MLE scale is %1.2f (%1.f, %1.2f)',
        rs[1], rs[2], rs[3], rs[4], rs[5], rs[6])



## Compare AIC for each distribution fit across runs
lapply(bootResults, function(ll) ll[, c('type', 'delAIC', 'AIC')]) %>% bind_rows() %>%
  group_by(type) %>%
  summarize(nbest = sum(delAIC == 0),
            fracBest = nbest/n(),
            meanDelAIC = median(delAIC)) %>%
  arrange(meanDelAIC)
sprintf('All models perform well, with median delta AIC < 2, although gamma and Weibull slightly outperform lognormal.')
sprintf('We proceed using a gamma distribution for computational convenience (its parameters can easily be related to the mean).')




