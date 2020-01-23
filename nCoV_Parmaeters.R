## USING LINELIST SO FAR, PLOT DISTRIBUTION OF TIME FROM ONSET TO HOSPITAL
library(dplyr)

## Define a funciton to convert from Excel date format to R date format
to_date = function(xx){
  as.character(xx) %>% as.Date(xx, format = '%d-%b-%y')
}

## Define a function to get days elapsed since 2020-01-01
days_from_origin = function(xx){
  (as.Date(xx)-as.Date('2020-01-01')) %>% as.numeric() -> elapsed
  2020+elapsed/365
}

## Import line list data
dat = read.csv(file = '2020_nCov/linelist.csv', stringsAsFactors = FALSE) %>% as.tbl() %>%
  select(-starts_with("X")) %>% ## Drop trailing columns
  mutate_at(vars(starts_with('date')), to_date) ## Convert to R date format
dat

## View histogram of times from onset to hospital admit
dat %>%
  filter(!is.na(dateOnset) & !is.na(dateAdmit)) %>%
  mutate(decimalOnset = days_from_origin(dateOnset),
         decimalAdmit = days_from_origin(dateAdmit),
         elapsed = dateAdmit-dateOnset) -> elapsedVec
elapsedVec %>%
  group_by(elapsed) %>%
  summarise(n = n()) -> elapsedSummary

## Plot
elapsedSummary %>%
  ggplot() +
  geom_bar(aes(x = elapsed, y = n), stat = "identity", color = 'purple3', fill = 'purple3')+
  theme_bw()+
  ggtitle('Days from onset to admission')
ggsave(filename = 'Onset_to_admit.pdf', width = 5, height = 5, units = 'in')

## Calculate the mean time from onset to hospitalization
elapsedSummary %>% pull(n) %>% summary()
sprintf(fmt = "The mean time from onset to admit is %2.1f days. The standard deviation is %2.1f days. n = %2.0f",
        elapsedSummary %>% pull(n) %>% mean(),
        elapsedSummary %>% pull(n) %>% sd(),
        elapsedSummary %>% pull(n) %>% sum())


## function to evaluate the neg log likelihood
nll = function(pars, dd, prof.val = NULL, prof.par = NULL){
  ## INPUTS - dd is the data (a vector of times elapsed),
  ##  pars is a named vector containing the shape (shp) and scale (scl) parameter of the gamma distribution
  ##  prof.val is an optional numeric input, which gives the value of the fixed parameter when profiling
  ## OUTPUTS - the negative log likelihood of the data given the parameters
  shp = pars['shp']
  scl = pars['scl']
  - (dgamma(dd, shape = shp, scale = scl, log = TRUE) %>% sum())
}


for(ii in 1:length(scaleVals)){
  optim(par = c(MLEs$par[1]), fn = nll, dd = fitDat, prof.val = scaleVals[ii], prof.par = 'scl', method = 'Brent', lower = .001, upper = 100)
  nll(pars =  c(MLEs$par[1]), dd = fitDat, prof.val = scaleVals[ii], prof.par = 'scl')
}

## Clean elapsed vec. If 0 days reported between onset and admit, round up to 0.5 (assume at least 12 hours elapsed)
elapsedVec %>%
  rowwise() %>% transmute(dd = ifelse(elapsed < 1, .5, elapsed)) %>%
  pull(dd) -> fitDat
## Use optim() to maximize the likelihood
optim(par = c(shp = 1.1, scl = 2), fn = nll, dd = fitDat) -> MLEs
MLEs
# 
# ## Run profiles on shape and scale pars
# scaleVals = seq(MLEs$par['scl']*.1, MLEs$par['scl']*.5, by = 0.001)
# scaleProf = sapply(scaleVals, FUN = function(xx){
#   optim(par = c('shp' = 1), fn = nll, fitDat, xx, 'scl', method = 'Brent', lower = 0, upper = 1e10)$value
# })
# 
# shapeProf = seq(MLEs$par['shp']*.1, MLEs$par['shp']*.5, by = 0.001)
# shapeProf = sapply(scaleVals, FUN = function(xx){
#   optim(par = c('scl' = 1), fn = nll, fitDat, xx, 'shp', method = 'Brent', lower = 0, upper = 1e10)$value
# })
# 
# par(mfrow = c(1,2))
# plot(scaleVals, scaleProf, 
#      main = "Profile on scale",
#      xlab = "scale value",
#      ylab = "profile neg log lk",
#      type = 'l')
# abline(h = MLEs$value, col = 'blue')
# abline(v = MLEs$par, col = 'red')
# 
# plot(shapeVals, shapeProf, 
#      main = "Profile on shape",
#      xlab = "shape value",
#      ylab = "profile neg log lk",
#      ylim = c(MLEs$value+c(-5, 50)), 
#      xlim = range(shapeVals[scaleProf < MLEs$value+50], na.rm = TRUE),
#      type = 'l')
# abline(h = MLEs$value, col = 'blue')
# abline(v = MLEs$par, col = 'red')


## Calculate the fraction of cases with fever
dat %>% 
  mutate(Fever = ifelse(is.na(Fever) | Fever != "Y", FALSE, TRUE)) %>%
  group_by(Fever) %>%
  summarise(n = n()) -> feverSummary
feverSummary

sprintf('%2.0f of %2.0f subjects (%2.2f percent) had a fever', feverSummary[2,2], sum(feverSummary[,2]), (feverSummary[2,2]/sum(feverSummary[,2]))*100)
