# FLIPPER
R-files associated with the 'FLexible Interpretation of Porewater Profiles and Estimation of Rates'

# Required packages
FLIPPER requires the following packages to be installed (the indentet packages are dependencies and should be installed together with the other packages)
- marelac
  -> shape
- signal
- fractaldim
  -> abind
- tcltk
- ReacTran
  -> rootSolve
  -> deSolve
- FME
  -> coda
- wavelets

If you want to run through the manual, you will need to install the 'knitr' package from CRAN.
Simply open the "Manual.Rmd" file in Rstudio, set the working directory to the file location, and press the 'Run document' button at the top of the page. 