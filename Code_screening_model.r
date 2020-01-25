

screen.passengers = function(d, del.d, f, g, sd = 1, sa =1, rd = 1, ra = 1, 
                             phi.d, incubation.d, pathogen, relative=0, split1=0, arrival_screen, departure_screen){
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
  screen.cases = function(case, arrival_screen, departure_screen){

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
    outputs=c(stopped.at.departure.fever, stopped.at.departure.risk, stopped.at.arrival.fever, stopped.at.arrival.risk,(cleared.at.arrival))
    names(outputs) = c('caught.dpt.fever', 'caught.dpt.risk', 'caught.arv.fever', 'caught.arv.risk', 'missed')
  }
  
  return(outputs)

  }


ff1=f
gg1=g

#Define the proportion of travellers that fall into each detection class (case)
cases1 = c(ff1*gg1, ff1*(1-gg1), (1-ff1)*gg1, (1-ff1)*(1-gg1))

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
         'missed.risk.only' = c3[5],
         'not.detectable' = c4[5]))
}

#Run the screen.cases function for the appropriate case and weight by the 
#  appropriate proportion of travellers
cases1[1]*screen.cases(1)+cases1[2]*screen.cases(2)+cases1[3]*screen.cases(3)+cases1[4]*screen.cases(4)

}


#Case 1 = has fever and aware of risk factors
#Case 2 = has fever and NOT aware of risk factors
#Case 3 = does NOT have fever and aware of risk factors
#Case 4 = does NOT have fever and NOT aware of risk factors
#The proportion of travellers that fall into each case (detection class)
#is determined by values of f and g
   