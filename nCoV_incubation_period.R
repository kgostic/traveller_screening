## Estimate the incubation period distribution
source('nCoV_Parmaeters.R')

## \\\\\\\\\\\\ INCUBATION PERIOD //////////////
outsideWuhan = read.csv('2020_nCov/outside_wuhan.csv', stringsAsFactors = FALSE) 
outsideWuhan %>% 
  select(age, sex, date_onset_symptoms, date_admission_hospital, date_arrive_Wuhan, date_depart_Wuhan) %>%
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


## Bootstrap to represent possible times of incubation
nBoot = 100
incubationPeriods[sample(incubationPeriods$id, size = nBoot, replace = TRUE),] %>%  ## Sample individual rows. Times of onset should be known.
                    rowwise() %>%
                    mutate(date_exposed = if(exact == TRUE){
                      date_depart_Wuhan
                      }else{ ## If date of exposure is uncertain, sample from all known possible dates
                        date_arrive_Wuhan + sample(0:as.numeric(date_depart_Wuhan-date_arrive_Wuhan), size = 1)
                      },
                      boot_incubation = max(((date_onset_symptoms-date_exposed) %>% as.numeric()), .5)) -> iBoot


incubationPeriods %>%
  pivot_longer(cols = c('max_incubation', 'min_incubation', 'mid_incubation'), values_to = 'period', names_to = 'minOrMax') %>%
  select(exact, minOrMax, period) %>%
  ggplot()+
  geom_density(aes(x = period, color = minOrMax, fill = minOrMax), alpha = .5, breaks = seq(-.5, 25.5, by = 1))

incubationPeriods %>%
  ggplot()+
  geom_histogram(aes(x = mid_incubation), alpha = .5, breaks = seq(-.5, 15.5, by = 1))


## Fit incubation period to distribution
nll_gamma = function(pars, dd){
  ## INPUTS - dd is the data (a vector of times elapsed),
  ##  pars is a named vector containing the shape (shp) and scale (scl) parameter of the gamma distribution
  ##  prof.val is an optional numeric input, which gives the value of the fixed parameter when profiling
  ## OUTPUTS - the negative log likelihood of the data given the parameters
  shp = pars['shp']
  scl = pars['scl']
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


## Fit to midpoints of raw data
incubation_fits_raw = list(
  gamma = optim(par = c(shp = 1, scl = 2), fn = nll_gamma, dd = incubationPeriods$mid_incubation),
  weibull = optim(par = c(shp = 1, scl = 2), fn = nll_weibull, dd = incubationPeriods$mid_incubation),
  lognormal = optim(par = c(mu = 1, sd = 2), fn = nll_lognormal, dd = incubationPeriods$mid_incubation)
)

incubationMLEs = data.frame(type = names(incubation_fits_raw),
                            par1 = sapply(incubation_fits_raw, function(xx) xx$par[1]),
                            par2 = sapply(incubation_fits_raw, function(xx) xx$par[2]),
                            nll = sapply(incubation_fits_raw, function(xx) xx$value)) %>%
  mutate(AIC = 2*2+2*nll) %>%
  arrange(AIC)
incubationMLEs
## Predictions for plotting
incubationFits = data.frame(xx = seq(0, 20, by = .01)) %>%
  mutate(gamma = dgamma(xx, shape = incubationMLEs %>% filter(type == 'gamma') %>% pull(par1), 
                        scale = incubationMLEs %>% filter(type == 'gamma') %>% pull(par2)),
         weibull = dweibull(xx, shape = incubationMLEs %>% filter(type == 'weibull') %>% pull(par1),
                            scale = incubationMLEs %>% filter(type == 'weibull') %>% pull(par2)),
         lognormal = dlnorm(xx, meanlog = log(incubationMLEs %>% filter(type == 'lognormal') %>% pull(par1)),
                            sdlog = log(incubationMLEs %>% filter(type == 'lognormal') %>%  pull(par2)))) %>%
  pivot_longer(gamma:lognormal, values_to = 'yy', names_to = 'fit')

## Fit to bootstrapped data
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
## Predictions for plotting
bootFits = data.frame(xx = seq(0, 20, by = .01)) %>%
  mutate(gamma = dgamma(xx, shape = bootMLEs %>% filter(type == 'gamma') %>% pull(par1), 
                        scale = bootMLEs %>% filter(type == 'gamma') %>% pull(par2)),
         weibull = dweibull(xx, shape = bootMLEs %>% filter(type == 'weibull') %>% pull(par1),
                            scale = bootMLEs %>% filter(type == 'weibull') %>% pull(par2)),
         lognormal = dlnorm(xx, meanlog = log(bootMLEs %>% filter(type == 'lognormal') %>% pull(par1)),
                            sdlog = log(bootMLEs %>% filter(type == 'lognormal') %>%  pull(par2)))) %>%
  pivot_longer(gamma:lognormal, values_to = 'yy', names_to = 'fit')


incubationPeriods %>%
  pivot_longer(cols = c(max_incubation, min_incubation, mid_incubation), names_to = 'kind', values_to = 'period') %>%
  group_by(kind, period) %>%
  summarise(n = n()) %>%
  ungroup() %>% group_by(kind) %>%
  mutate(freq = n/n()) %>%
  filter(kind == 'mid_incubation') %>%
  ggplot() +
  geom_bar(aes(x = period, y = freq), stat = "identity", fill = 'purple3', color = 'purple3')+
  theme_bw()+
  geom_line(data = incubationFits, aes(x = xx, y = yy, color = fit)) +
  ggtitle('Incubation period') +
  ylab('density')


iBoot %>% rowwise() %>%
  mutate(boot_incubation = ifelse(boot_incubation == .5, 0, boot_incubation)) %>%ungroup() %>%
  group_by(boot_incubation) %>%
  summarise(n = n(),
            freq = n()/nBoot) %>% 
  ungroup() %>%
  ggplot()+
  geom_bar(aes(x=boot_incubation, y = freq), stat = "identity", fill = 'purple3', color = 'purple3')+
  theme_bw()+
  geom_line(data = incubationFits, aes(x = xx, y = yy, color = fit)) +
  ggtitle('Bootstrapped incubation period') +
  ylab('density')
ggsave('2020_nCov/incubation_pd_bootstrap.png', width = 5, height = 5, units = 'in', dpi = 320)


sprintf('The mean incubation period is %2.1f days, with sd %2.2f. The upper bound is %2.1f days and the lower bound is %2.1f days',
        incubationPeriods %>% pull(mid_incubation) %>% mean(),
        incubationPeriods %>% pull(mid_incubation) %>% sd(),
        incubationPeriods %>% pull(max_incubation) %>% max(),
        incubationPeriods %>% pull(min_incubation) %>% min())


sprintf('The bootstrapped incubation period is %2.1f days, with sd %2.2f. The upper bound is %2.1f days and the lower bound is %2.1f days',
        iBoot %>% pull(boot_incubation) %>% mean(),
        iBoot %>% pull(boot_incubation) %>% sd(),
        iBoot %>% pull(boot_incubation) %>% max(),
        iBoot %>% pull(boot_incubation) %>% min())
