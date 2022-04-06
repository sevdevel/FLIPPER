# FLIPPER
R-files associated with the 'FLexible Interpretation of Porewater Profiles and Estimation of Rates'

# Required packages
FLIPPER requires the following packages to be installed (the packages in brackets are dependencies and should be installed automatically 
together with the other packages)
- marelac (shape)
- signal
- fractaldim (abind)
- tcltk
- ReacTran (rootSolve,deSolve)
- FME (coda)
- wavelets

# User manual
There is a static version of the user manual, which is available in html and pdf version ('Manual_Static').
This version is convenient if you do not have R installed, and want to have a quick overview of how the FLIPPER package works.

If you are already accustomed to R and Rstudio and want to spice up your life, you can also run the interactive version of the Manual in R markdown 
('Manual.rmd'). Make sure you have the 'knitr' package installed from CRAN. Next, simply open the "Manual.Rmd" file in Rstudio, and press the 
'Run document' button at the top of the page. 

A 'userInterface example.R' script is provided, which contains the calling of FLIPPER as it is outlined in the manual. 