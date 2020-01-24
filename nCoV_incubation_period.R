## Estimate the incubation period distribution


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
incubationPeriods$min_incubation[which(incubationPeriods$min_incubation < 0)] = 0 # Clean up one case with symptom onset before departure 
incubationPeriods %>%
  mutate(mid_incubation = (min_incubation + max_incubation)/2 ) -> incubationPeriods

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

incubationFits = data.frame(xx = seq(0, 20, by = .01)) %>%
  mutate(gamma = dgamma(xx, shape = incubationMLEs %>% filter(type == 'gamma') %>% pull(par1), 
                        scale = incubationMLEs %>% filter(type == 'gamma') %>% pull(par2)),
         weibull = dweibull(xx, shape = incubationMLEs %>% filter(type == 'weibull') %>% pull(par1),
                            scale = incubationMLEs %>% filter(type == 'weibull') %>% pull(par2)),
         lognormal = dlnorm(xx, meanlog = log(incubationMLEs %>% filter(type == 'lognormal') %>% pull(par1)),
                            sdlog = log(incubationMLEs %>% filter(type == 'lognormal') %>%  pull(par2)))) %>%
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

