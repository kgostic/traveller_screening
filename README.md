## Estimated effectiveness of symptom and risk screening to prevent the spread of COVID-19
**Authors: Gostic, Gomez, Mumah, Kucharski, Lloyd-Smith.**
Last update: Feb. 21, 2020
**https://www.medrxiv.org/content/10.1101/2020.01.28.20019224v2**

### :arrow_right: [[Read the study]](https://elifesciences.org/articles/55570).

### :arrow_right: [[Use our model to generate your own outputs]](https://lloyd-smithlab.shinyapps.io/travelScreeningModel/).

### The source code we used to perform these analyses runs in R, and is all in this repo. If you want to use or modify our code, keep reading.
(The exact release we used to generate the figures in the manuscript is [[here]](https://github.com/kgostic/traveller_screening/releases/tag/v2.1).

Feel free to use or modify this code. We just ask that you cite us:
 * Gostic K, Kucharski AJ, Lloyd-Smith JO. (2015) eLife 4:e05564.
 * Gostic K; Gomez, ACR; Mummah, RO; Kucharski, AJ; Lloyd-Smith JO. (2020). medRxiv: https://doi.org/10.1101/2020.01.28.20019224.
 
 -----
 
 contact kgostic [at] uchicago.edu with questions
 
 ----




## These scripts run the model and generate output from: 



This code by K. Gostic, adapted from a previous version implemented by A. Kucharski.





## Repository contents:
All code runs in R.
* The nCov_manuscript_plots.R script runs all analyses and generates all plots and data outputs.
    * .RData files contain raw, simulated data.
        * "a" indicates outputs generated based on arrival screening only. "d" indicates departure screening only. "ad" indicates both. 
        * "flat" indicates outputs generated assuming a stable epidemic. Outputs without the "flat" flag were generated assuming a growing source epidemic.
  * 2002_nCoV/ contains all plot and source data outputs.
      * Source data .csv files contain processed data frames used to generate the corresponding plot.
      * Plot names corresond to their ID in the manuscript (e.g. Fig. 3, Fig 3S1).
      
Code for the shiny app is hosted elsewhere. These scripts here generate figures shown in the manuscript.

------

