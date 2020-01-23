# This is code for Figure 3 and 4
# This script runs the external functions "Code_Screening_model.R" and "Code_distribution_functions.R"

# Source external functions and load packages
source("Code_Screening_model.R")
source("Code_distribution_functions.R")
library(ggplot2)
library(grid)
library(gridExtra)
coltab=c("red","blue","orange","green","purple","black")

# - - - - - - - - - - - - - - - - - - - - - - -
# Plot proportion of travellers identified vs time from exposure to departure (Figure 3)
# - - - - - - - - - - - - - - - - - - - - - - -
#Specifiticty of pathogens of interest. Here we include pathogen-specific distributions 
#   to run the script for H7N9
pathtab=c("SARS", "MERS")
pathtablab=c("SARS CoV", "MERS CoV")

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

#This loop repeats for each pathogen listed in 'pathtab' above
for(i in 1:length(pathtab)){
  
pathogen=pathtab[i]

#Define a series of fixed times since exposure
arrive.time=seq(0,15,0.1)
#detect.break stores the following outputs for fixed time since exposure: "stop.d.fever","stop.d.risk","stop.a.fever","stop.a.risk","c.a"
detect.break=matrix(NA,nrow=length(arrive.time),ncol=5)

for(jj in 1:length(arrive.time)){

  incubation.d<-function(d){incubation.function(d,pathogen)} #'incubation.function' is defined in "Code_distribution_functions.R"

    del.d=1
    f1=fg.distn("f",pathogen)[1]
    g1=fg.distn("g",pathogen)[1]
    pp1split<-function(x,dCI,pathogen,relative){screen.passengers(x, del.d, f1, g1, sd, sa, rd, ra, phi.d, incubation.d,pathogen,relT,1)}
    # Calculate median based on incubation and fever
    detect.break[jj,]=pp1split(arrive.time[jj],1,pathogen,relT)

}

colnames(detect.break)=c("stop.d.fever","stop.d.risk","stop.a.fever","stop.a.risk","not.caught")

#Figure 3: plot the fraction of travellers caught during each screening stage for each 
#  fixed time since exposure
plot0=0*detect.break[,1]
plot1=detect.break[,1]
plot2=plot1+detect.break[,2]
plot3=plot2+detect.break[,3]
plot4=plot3+detect.break[,4]
plot5=plot4+detect.break[,5]
plotT=0*detect.break[,1]

colors = c(rgb(0,1,0),rgb(0,0,1),rgb(0,0.5,0),rgb(0,0,0.5),rgb(0.7,0.2,0.2)) 

plot(arrive.time,plot1,type="l",main=pathtablab[i],ylim=c(0,1),xlab="exposure to departure (days)",ylab="proportion of travellers")

polygon(c(arrive.time, rev(arrive.time)), c(plot1, rev(plot0)),col = colors[1], border = FALSE)
polygon(c(arrive.time, rev(arrive.time)), c(plot2, rev(plot1)),col = colors[2], border = FALSE)
polygon(c(arrive.time, rev(arrive.time)), c(plot3, rev(plot2)),col = colors[3], border = FALSE)
polygon(c(arrive.time, rev(arrive.time)), c(plot4, rev(plot3)),col = colors[4], border = FALSE)
polygon(c(arrive.time, rev(arrive.time)), c(plot5, rev(plot4)),col = colors[5], border = FALSE)
lines(arrive.time,plot2,col="white",lwd=3)
outcome1=exposure.distn(1,pathogen,flat=1)
lines(c(outcome1,outcome1),c(0,1),col="black",lwd=2,lty=2)

}

dev.copy(pdf,paste0("Proportion_identified_", pathogen, ".pdf"),width=10,height=7)
dev.off()


# - - - - - - - - - - - - - - - - - - - - - - -
# Plot distribution of time from exposure to departure in growing/stable epidemic (Figure 4, S1)
# - - - - - - - - - - - - - - - - - - - - - - -

r0=2 # Reproduction number

pathtab=c("H7N9")
xx=seq(0,15,0.1)
typ=1
par(mfrow=c(1,2))

plot(xx,exposure.outcome(xx,pathogen,typ,flat=0),lty=1,xlab="time from exposure to departure",ylab="probability density",ylim=c(0,0.3),col="white")
for(ii in 1:length(pathtab)){
  lines(xx,exposure.outcome(xx,pathtab[ii],typ,flat=0),col=coltab[ii])
}

plot(xx,exposure.outcome(xx,pathogen,typ,flat=1),lty=1,xlab="time from exposure to departure",ylab="probability density",ylim=c(0,0.3),col="white")
for(ii in 1:length(pathtab)){
  lines(xx,exposure.outcome(xx,pathtab[ii],typ,flat=1),col=coltab[ii])
}

dev.copy(pdf,"Time_from_exposure_to_departure.pdf",width=10,height=7)
dev.off()


# - - - - - - - - - - - - - - - - - - - - - - -
# Plot boostrap confidence intervals for effectiveness of different screening measures
# - - - - - - - - - - - - - - - - - - - - - - -

ftime=24   # Travel time

#Growing (0) or flat (1) epidemic
flatA=0
# Overall (0) or arrival (1) screening
relT=0

# Store median and CI for fever only screening outcomes
detect.mat=matrix(NA,nrow=length(pathtab),ncol=length(ftime))
detect.matA=matrix(NA,nrow=length(pathtab),ncol=length(ftime))
detect.matB=matrix(NA,nrow=length(pathtab),ncol=length(ftime))
# Store median and CI for fever & exposure screening outcomes
detect.Smat=matrix(NA,nrow=length(pathtab),ncol=length(ftime))
detect.SmatA=matrix(NA,nrow=length(pathtab),ncol=length(ftime))
detect.SmatB=matrix(NA,nrow=length(pathtab),ncol=length(ftime))
# Store median and CI for exposure only screening outcomes
detect.Cmat=matrix(NA,nrow=length(pathtab),ncol=length(ftime))
detect.CmatA=matrix(NA,nrow=length(pathtab),ncol=length(ftime))
detect.CmatB=matrix(NA,nrow=length(pathtab),ncol=length(ftime))



for(i in 1:length('pathtab')){
  pathogen=pathtab[i]
  incubation.d<-function(d){incubation.function(d,pathogen)}
  
  for(j in 1:length(ftime)){
    del.d=ftime[j]/24 #Convert time in flight from hours to days
    sd=0.7            #fever screening efficacy at departure
    sa=sd             #fever screening efficacy at arrival
    
    blen=20 # Bootstrap samples
    popn=50 # Number of infected travellers
    bootstrap=c(1:blen)  #Output for fever screening only
    bootstrapS=c(1:blen) #Output for risk AND fever screening
    bootstrapC=c(1:blen) #Output for risk screening only
    
    for(ii in 1:blen){
      pp1<-function(x,dCI,pathogen,relative,f,g){screen.passengers(x, del.d, f, g, sd, sa, 0 , 0, phi.d, incubation.d,pathogen,relative)}
      pp1S<-function(x,dCI,pathogen,relative,f,g){screen.passengers(x, del.d, f, g, sd, sa, 0.25 , 0.25, phi.d, incubation.d,pathogen,relative)}
      pp1C<-function(x,dCI,pathogen,relative,f,g){screen.passengers(x, del.d, f, g, 0, 0, 0.25 , 0.25, phi.d, incubation.d,pathogen,relative)}
      arrive.n=sapply(c(1:popn),function(x){exposure.distn(runif(1, 0, 1),pathogen,flat=flatA)})
      bb01=sapply(arrive.n,function(x){
        f0=rbinom(1,fg.distn("f",pathogen)[2],fg.distn("f",pathogen)[1])/fg.distn("f",pathogen)[2]
        g0=rbinom(1,fg.distn("g",pathogen)[2],fg.distn("g",pathogen)[1])/fg.distn("g",pathogen)[2]
        pp1(x,1,pathogen,relT,f0,g0)})
      bb02=sapply(bb01,function(x){ifelse(x<runif(1, 0, 1),0,1)})
      bootstrap[ii]=sum(bb02)/popn
      bb01=sapply(arrive.n,function(x){
        f0=rbinom(1,fg.distn("f",pathogen)[2],fg.distn("f",pathogen)[1])/fg.distn("f",pathogen)[2]
        g0=rbinom(1,fg.distn("g",pathogen)[2],fg.distn("g",pathogen)[1])/fg.distn("g",pathogen)[2]
        pp1S(x,1,pathogen,relT,f0,g0)})
      bb02=sapply(bb01,function(x){ifelse(x<runif(1, 0, 1),0,1)})
      bootstrapS[ii]=sum(bb02)/popn
      bb01=sapply(arrive.n,function(x){
        f0=rbinom(1,fg.distn("f",pathogen)[2],fg.distn("f",pathogen)[1])/fg.distn("f",pathogen)[2]
        g0=rbinom(1,fg.distn("g",pathogen)[2],fg.distn("g",pathogen)[1])/fg.distn("g",pathogen)[2]
        pp1C(x,1,pathogen,relT,f0,g0)})
      bb02=sapply(bb01,function(x){ifelse(x<runif(1, 0, 1),0,1)})
      bootstrapC[ii]=sum(bb02)/popn
    }
    
    #Store results for fever screening only
    s1=sort(bootstrap)
    detect.mat[i,j]=1-s1[ceiling(0.5*blen)]
    detect.matA[i,j]=1-s1[ceiling(0.025*blen)]
    detect.matB[i,j]=1-s1[ceiling(0.975*blen)]
    
    #Store results for fever AND risk screening
    s1S=sort(bootstrapS)
    detect.Smat[i,j]=1-s1S[ceiling(0.5*blen)]
    detect.SmatA[i,j]=1-s1S[ceiling(0.025*blen)]
    detect.SmatB[i,j]=1-s1S[ceiling(0.975*blen)]
    
    #Store results for risk screening only
    s1C=sort(bootstrapC)
    detect.Cmat[i,j]=1-s1C[ceiling(0.5*blen)]
    detect.CmatA[i,j]=1-s1C[ceiling(0.025*blen)]
    detect.CmatB[i,j]=1-s1C[ceiling(0.975*blen)]
  }
  
}


rownames(detect.mat)=pathtab
coltab=c("red","blue","orange","green","purple","black")


# Plot confidence intervals for specific time (uses above data)
picktime=24 #This should be a vector of all flight times tested. 
            #As a default we use a single flight time: 24 hours
nnF=c(1:length(ftime))[ftime==picktime]

dataFlight=matrix(NA,nrow=3*length(pathtab),5)
dataFlight=as.data.frame(dataFlight)

detect.Smat[3,]=detect.mat[3,]
detect.SmatA[3,]=detect.matA[3,]
detect.SmatB[3,]=detect.matB[3,]

dataFlight[,3]=c(detect.mat[,nnF],detect.Smat[,nnF],detect.Cmat[,nnF])
dataFlight[,4]=c(detect.matA[,nnF],detect.SmatA[,nnF],detect.CmatA[,nnF])
dataFlight[,5]=c(detect.matB[,nnF],detect.SmatB[,nnF],detect.CmatB[,nnF])
dataFlight[1:length(pathtab),2]="Fever"
dataFlight[(length(pathtab)+1):(2*length(pathtab)),2]="Fever & risk"
dataFlight[(2*length(pathtab)+1):(3*length(pathtab)),2]="Risk"
dataFlight[,1]=c(pathtab,pathtab,pathtab)

# Save data
write.csv(dataFlight,paste("Dataflight_",flatA,"_rel",relT,"_pop",popn,".csv",sep=""))

colnames(dataFlight)=c("pathogen","Screening","mid","ci1","ci2")
dataFlight$pathogen <- as.character(dataFlight$pathogen)
dataFlight$pathogen<- factor(dataFlight$pathogen, levels=unique(dataFlight$pathogen))
titlen=c("Total screening (growing epidemic)","Arrival screening (growing epidemic)","Total screening (stable epidemic)","Arrival screening (stable epidemic)")



# - - - - - - - - - - - - - - - - - - - - - - -
# Plot figure 4
# - - - - - - - - - - - - - - - - - - - - - - -
pd <- position_dodge(.4)
plot6<-ggplot(dataFlight, aes(x=pathogen, y=mid,colour=Screening)) + 
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(size = 1,rgb(0.95,0.95,0.95)),panel.background=element_blank(),panel.border=element_rect(colour = "black", fill=FALSE))+
  geom_errorbar(aes(ymin=ci1, ymax=ci2),  width=.1,position=pd) +
  geom_point(aes(ymin=mid, ymax=mid),size=3,position=pd) +
  xlab("") +
  ggtitle(titlen[kk]) +
  ylim(0,1) +
  ylab("proportion missed") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(filename="Proportion_of_cases_missed.pdf",plot=plot6,width=6,height=4)





