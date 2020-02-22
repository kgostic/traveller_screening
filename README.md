## These scripts run the model and generate output from: 

Estimated effectiveness of symptom and risk screening to prevent the spread of COVID-19
Authors: Gostic, Gomez, Mumah, Kucharski, Lloyd-Smith.
https://www.medrxiv.org/content/10.1101/2020.01.28.20019224v2

This code by K. Gostic, adapted from a previous version implemented by A. Kucharski.


Model and methods adapted from  Gostic KM, Kucharski AJ, Lloyd-Smith JO. (2015) eLife 4:e05564.


## Repository contents:
All code runs in R.
* The nCov_manuscript_plots.R script runs all analyses and generates all plots and data outputs.
    * .RData files contain raw, simulated data.
        * "a" indicates outputs generated based on arrival screening only. "d" indicates departure screening only. "ad" indicates both. 
        * "flat" indicates outputs generated assuming a stable epidemic. Outputs without the "flat" flag were generated assuming a growing source epidemic.
  * 2002_nCoV/ contains all plot and source data outputs.
      * Source data .csv files contain processed data frames used to generate the corresponding plot.
      * Plot names corresond to their ID in the manuscript (e.g. Fig. 3, Fig 3S1).
      
## See https://lloyd-smithlab.shinyapps.io/travelScreeningModel/ for an interactive Shiny app.
Code for the shiny app is hosted elsewhere. These scripts here generate figures shown in the manuscript.

------
contact: kgostic [at] uchicago.edu
