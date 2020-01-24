## USING LINELIST SO FAR, PLOT DISTRIBUTION OF TIME FROM ONSET TO HOSPITAL
library(dplyr)

## Define a funciton to convert from Excel date format to R date format
to_date = function(xx){
  as.character(xx) %>% as.Date(xx, format = '%d-%b-%y')
}

to_date_outsideWuhan = function(xx){
  as.character(xx) %>% as.Date(xx, format = '%d.%m.%Y')
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



## Calculate the mean time from onset to hospitalization
elapsedSummary %>% pull(n) %>% summary()
sprintf(fmt = "The mean time from onset to admit is %2.1f days. The standard deviation is %2.1f days. n = %2.0f",
        elapsedSummary %>% pull(n) %>% mean(),
        elapsedSummary %>% pull(n) %>% sd(),
        elapsedSummary %>% pull(n) %>% sum())


## function to evaluate the neg log likelihood
nll_gamma = function(pars, dd, prof.val = NULL, prof.par = NULL){
  ## INPUTS - dd is the data (a vector of times elapsed),
  ##  pars is a named vector containing the shape (shp) and scale (scl) parameter of the gamma distribution
  ##  prof.val is an optional numeric input, which gives the value of the fixed parameter when profiling
  ## OUTPUTS - the negative log likelihood of the data given the parameters
  shp = pars['shp']
  scl = pars['scl']
  - (dgamma(dd, shape = shp, scale = scl, log = TRUE) %>% sum())
}


## Clean elapsed vec. If 0 days reported between onset and admit, round up to 0.5 (assume at least 12 hours elapsed)
elapsedVec %>%
  rowwise() %>% transmute(dd = ifelse(elapsed < 1, .5, elapsed)) %>%
  pull(dd) -> fitDat
## Use optim() to maximize the likelihood
optim(par = c(shp = 1.1, scl = 2), fn = nll_gamma, dd = fitDat) -> MLEs
MLEs

## Run profiles on shape and scale pars
scaleVals = seq(MLEs$par['scl']*.1, MLEs$par['scl']*5, by = 0.001) ## Set a grid of scale values on which to profile
## Write function to calculat profile value
scaleProfFun = function(pars, scale.val){
  - (dgamma(fitDat, shape = pars[1], scale = scale.val, log = TRUE) %>% sum()) 
}
## Calculate profile
scaleProf = sapply(scaleVals, FUN = function(xx){
  optim(par = c('shp' = 1), fn = scaleProfFun, scale.val = xx, method = 'Brent', lower = 0, upper = 1000)
}) %>% t() %>% as.data.frame()


## Set a grid of shape values on which to profile
shapeVals = seq(MLEs$par['shp']*.1, MLEs$par['shp']*5, by = 0.001) ## Set a grid of scale values on which to profile
## Write a function to calculate profile
shapeProfFun = function(pars, shape.val){
  - (dgamma(fitDat, shape = shape.val, scale = pars[1], log = TRUE) %>% sum()) 
}
## Calculate profile
shapeProf = sapply(shapeVals, FUN = function(xx){
  optim(par = c('scl' = 1), fn = shapeProfFun, shape.val = xx, method = 'Brent', lower = 0, upper = 1000)
}) %>% t() %>% as.data.frame()


## Plot
pdf('Onset_to_admit_profiles.pdf')
par(mfrow = c(1,2))
plot(scaleVals, scaleProf$value,
     main = "Profile on scale",
     xlab = "scale value",
     ylab = "profile neg log lk",
     type = 'l')
abline(h = MLEs$value, col = 'blue')
abline(v = MLEs$par['scl'], col = 'red')

plot(shapeVals, shapeProf$value,
     main = "Profile on shape",
     xlab = "shape value",
     ylab = "profile neg log lk",
     type = 'l')
abline(h = MLEs$value, col = 'blue')
abline(v = MLEs$par['shp'], col = 'red')
dev.off()

sprintf('Time from onset to admission follows a Gamma distribution with shape %2.2f, and scale %2.2f', MLEs$par[1], MLEs$par[2])

gammaFits = data.frame(xx = seq(0, 20, by = .01)) %>%
  mutate(pdf = dgamma(xx, shape = MLEs$par[1], scale = MLEs$par[2]),
         cdf = pgamma(xx, shape = MLEs$par[1], scale = MLEs$par[2]))

## Plot
elapsedSummary %>%
  mutate(n = n/sum(n)) %>%
  ggplot() +
  geom_bar(aes(x = elapsed, y = n), stat = "identity", color = 'purple3', fill = 'purple3')+
  theme_bw()+
  geom_line(data = gammaFits, aes(x = xx, y = pdf)) +
  ggtitle('Days from onset to admission') +
  ylab('density')
ggsave(filename = 'Onset_to_admit.pdf', width = 5, height = 5, units = 'in')


## Calculate the fraction of cases with fever
dat %>% 
  mutate(Fever = ifelse(is.na(Fever) | Fever != "Y", FALSE, TRUE)) %>%
  group_by(Fever) %>%
  summarise(n = n()) -> feverSummary
feverSummary

sprintf('%2.0f of %2.0f subjects (%2.2f percent) had a fever', feverSummary[2,2], sum(feverSummary[,2]), (feverSummary[2,2]/sum(feverSummary[,2]))*100)

