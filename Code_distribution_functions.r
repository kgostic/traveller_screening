### This script defines functions to be called by the master script, "Plot_results.R"

# --------------
# Distribution for incubation period (H7N9 given)
# ------------
incubation.function<-function(x,pathogen,type=1){
  #Input x, the time since exposure
  #Output the incubation cdf (type=1) or pdf (type=2)
  
  # \\\\\\\\\\  H7N9 inputs //////
  if(pathogen == "H7N9"){  
    H7N9.c.nsamp=32
    param.H7N9.fn=c(4.93, 0.63)
    H7N9.fn<-function(x){pgamma(x, param.H7N9.fn[1] ,scale = param.H7N9.fn[2], lower.tail = TRUE,log.p = FALSE)}
    H7N9.fnD<-function(x){dgamma(x, param.H7N9.fn[1] ,scale = param.H7N9.fn[2])}
    
    if(type==1){
      functionpath<-H7N9.fn
      #functionpath(x)
    }else{
      functionpath<-H7N9.fnD
      #functionpath(x)
    }
  }
  
 # \\\\\\  MERS inputs //////
  if(pathogen == "MERS"){  
    MERS.c.nsamp=7+23
    param.MERS.fn=c(log(5.5), log(2.5))
    MERS.fn<-function(x){plnorm(x, meanlog = param.MERS.fn[1], sdlog = param.MERS.fn[2])}
    MERS.fnD<-function(x){dlnorm(x, meanlog = param.MERS.fn[1], sdlog =  param.MERS.fn[2])}
    
    if(type==1){
      functionpath<-MERS.fn
      #functionpath(x)
    }else{
      functionpath<-MERS.fnD
      #functionpath(x)
    }
  }
  
  
  # \\\\\\  SARS inputs //////
  if(pathogen == "SARS"){  
    SARS.c.nsamp=57
    param.SARS.fn=c(2.62, 2.43)
    SARS.fn<-function(x){pgamma(x, param.SARS.fn[1], scale = param.SARS.fn[2])}
    SARS.fnD<-function(x){dgamma(x, param.SARS.fn[1], scale = param.SARS.fn[2])}
    
    if(type==1){
      functionpath<-SARS.fn
      #functionpath(x)
    }else{
      functionpath<-SARS.fnD
      #functionpath(x)
    }
  }
  functionpath(x)
}


# --------------
# Distributions for fever and known exposure
# --------------
fg.distn<-function(type,pathogen){
  #input pathogen name; parameter type (1 = fever, 2 = known risk)
  #output c(value, total sample size)
  #"value" is calculated as the mean proportion reported across studies reviewed, weighted
  # by the study-specific sample size
  #"total sample size" is the sum of sample sizes from all studies considered
  
  # Input pathogen-specific data:
  if(pathogen == "H7N9"){ 
    # Fever study sample sizes
    d1fn=c(85,46,111,46)
    # Proportions that display fever
    d1fv=c(0.7,1,1,1)
    
    # Exposure study sample sizes
    d1gn=c(123,111,46)
    # Proportions with known exposure
    d1gv=c(0.75,0.56,0.78)
  }
  
  
  # \\\\\\ MERS ///////
  if(pathogen == "MERS"){ 
    # Fever study sample sizes
    d1fn=c(23)
    # Proportions that display fever
    d1fv=c(.87)
    
    ######## THIS NEEDS AN UPDATE! 
    # Exposure study sample sizes
    d1gn=c(1000)
    # Proportions with known exposure
    d1gv=c(0)
  }
  
  # \\\\\\ SARS ///////
  if(pathogen == "SARS"){ 
    # Fever study sample sizes
    d1fn=c(1452)
    # Proportions that display fever
    d1fv=c(0.94)
    
    ######## THIS NEEDS AN UPDATE! 
    # Exposure study sample sizes
    d1gn=c(1192)
    # Proportions with known exposure
    d1gv=c(0.29)
  }
  
  
  if(type=="f"){
    c(sum(d1fv*d1fn)/sum(d1fn),sum(d1fn))}else{
      c(sum(d1gv*d1gn)/sum(d1gn),sum(d1gn))}
}


# --------------
# Binomial density function
# --------------
fg.binomial<-function(propn,type,pathogen){
  #Input propn = a vector of possible outcomes when sampling individuals from the population
  # type (1 = fever, 2 = risk awareness (g)), pathogen = "H7N9"
  #Output the binomial density corresponding to each value in "propn"
  #I.e., this function estimates the proability of observing some value of f_hat or
  # g_hat, given that f and g are the true population parameters
  
  aa=fg.distn(type,pathogen) #call the appropriate f, g values from the function above
  dbinom(round(propn*aa[2]),round(aa[2]),aa[1])
}


# ---------
# Calculate infection age distributions for a growing epidemic
# ---------
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

# ---------
# Output infection age ditributions for either a growing or stable epidemic
# ---------
exposure.outcome<-function(x,pathogen,type=1,flat=0){
  #pathogen-specific median time to outcome (median exposure -> onset + onset -> outcome)
  if(pathogen == "H7N9"){f1<-3.1+4.2}
  if(pathogen == "MERS"){f1<-(5.2*7+5.5*23)/(7+23)+5}  ## (5.2*7+5.5*23)/(7+23) is a weighted average based on sample sizes of individual studies
  if(pathogen == "SARS"){f1<-6.4+4.85}
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

# ---------
# Use uniform random variable x to make random draw from specified infection age distribution
# ---------
exposure.distn<-function(x,pathogen,type=2,flat=0){
  a=0
  while(x>exposure.outcome(a,pathogen,type,flat)){a=a+0.1}
  a
}

