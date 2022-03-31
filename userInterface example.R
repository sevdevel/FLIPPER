################################################################################
# Author: Sebastiaan van de Velde (sebastiaan.van.de.velde@ulb.be)
#         Laurine Burdorf         
#
# Example how to use the user interface function of fitprofile to analyse porewater profiles 
#
################################################################################

source("user interface function.R")
source("FLIPPER_plotfunction.R")

#===============================================================================
# Preparing a dataset
#===============================================================================

# load datasets 

load("datasets to test/O2_testprofiles.RData")

# show dataset contents
head(profile.A)

# create dataframe that can be analyzed by FLIPPER
input.profile <- as.data.frame(cbind(profile.A$C*1e3,profile.A$x*1e-2))
colnames(input.profile) <- c('C','x')

# show dataset contents
head(input.profile)

# prepare other required input
env.parms <- list(TC=10.0,S=35.0)
por.cte <- 0.8

#===============================================================================
# The gradient method
#===============================================================================

# Run FLIPPER
test <- FLIPPER.func(input=input.profile,species=c("O2"),por.cte=por.cte,tort.dep=1, 
                     env.parms=env.parms,
                     method="gradient")
# plot result
plot.FLIPPER(test)

# print output
print(test)

#===============================================================================
# The discrete method
#===============================================================================

# prepare discrete parameter list
discrete.parms <- list()

# Run FLIPPER
test <- FLIPPER.func(input=input.profile,species=c("O2"),por.cte=por.cte, 
                     env.parms=env.parms,
                     discrete.parms=discrete.parms,
                     method="discrete")
# plot result
plot.FLIPPER(test)

# print output
print(test)


# change boundary conditions
discrete.parms <- list(UBC="flux.up",LBC="conc.down",flux.up=-8.0)
test <- FLIPPER.func(input=input.profile,species=c("O2"),por.cte=por.cte, 
                     env.parms=env.parms,
                     discrete.parms=discrete.parms,
                     method="discrete")
plot.FLIPPER(test)

# change depth of analysis
discrete.parms <- list(L.down=0.003)
test <- FLIPPER.func(input=input.profile,species=c("O2"),por.cte=por.cte, 
                     env.parms=env.parms,
                     discrete.parms=discrete.parms,
                     method="discrete")
plot.FLIPPER(test)

# change number of zones to test
discrete.parms <- list(i.end=3)
test <- FLIPPER.func(input=input.profile,species=c("O2"),por.cte=por.cte, 
                     env.parms=env.parms,
                     discrete.parms=discrete.parms,
                     method="discrete")
plot.FLIPPER(test)

# immediately start lumping from a given number of zones
discrete.parms <- list(initial.zones=5)
test <- FLIPPER.func(input=input.profile,species=c("O2"),por.cte=por.cte, 
                     env.parms=env.parms,
                     discrete.parms=discrete.parms,
                     method="discrete")
plot.FLIPPER(test)

#===============================================================================
# The continuous method
#===============================================================================

# uniform filter window of 21 (n=15)
continuous.parms <- list(n.C=15,n.J=15,n.R=15)
test <- FLIPPER.func(input=input.profile,species=c("O2"),por.cte=por.cte, 
                     env.parms=env.parms,
                     continuous.parms=continuous.parms,
                     method="continuous")
plot.FLIPPER(test)

# Run FLIPPER continuous and interactive, without uniform filter window
continuous.parms <- list(optimal.window.size="interactive",n.uniform=FALSE)
test <- FLIPPER.func(input=input.profile,species=c("O2"),por.cte=por.cte, 
                     env.parms=env.parms,
                     continuous.parms=continuous.parms,
                     method="continuous")
plot.FLIPPER(test)

# Run FLIPPER continuous and interactive, with uniform filter window
continuous.parms <- list(optimal.window.size="interactive",n.uniform=TRUE)
test <- FLIPPER.func(input=input.profile,species=c("O2"),por.cte=por.cte, 
                     env.parms=env.parms,
                     continuous.parms=continuous.parms,
                     method="continuous")
plot.FLIPPER(test)

#===============================================================================
# Run all methods
#===============================================================================

gradient.parms <- c()
discrete.parms <- list()
continuous.parms <- list(optimal.window.size="interactive",n.uniform=TRUE)
test <- FLIPPER.func(input=input.profile,species=c("O2"),por.cte=por.cte, 
                     env.parms=env.parms,
                     gradient.parms=gradient.parms,
                     discrete.parms=discrete.parms,
                     continuous.parms=continuous.parms,
                     method="all")
plot.FLIPPER(test)

#===============================================================================
# Test dataset of classical porewater data with an electrical field
#===============================================================================

# load dataset
load("datasets to test/test.data.1.Rdata")

# gradient method
test <- FLIPPER.func(input=as.data.frame(test.data$input),species=test.data$species,E.cte=c(-0.1,0,0.03),
                     env.parms=test.data$env.parms,
                     method="gradient")

x11(height=20, width=20)
plot.FLIPPER(test)

# discrete method
discrete.parms <- list(LBC="conc.down",C.down=test.data$input$C[nrow(test.data$input)])
test <- FLIPPER.func(input=as.data.frame(test.data$input),species=test.data$species,E.cte=c(-0.1,0,0.03),
                     env.parms=test.data$env.parms,
                     discrete.parms=discrete.parms,
                     method="discrete")

x11(height=20, width=20)
plot.FLIPPER(test)

# continuous method
continuous.parms <- list()
test <- FLIPPER.func(input=as.data.frame(test.data$input),species=test.data$species,E.cte=c(-0.1,0,0.03),
                     env.parms=test.data$env.parms,
                     continuous.parms=continuous.parms,
                     method="continuous")

x11(height=20, width=20)
plot.FLIPPER(test)

