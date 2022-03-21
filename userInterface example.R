################################################################################
# Author: Sebastiaan van de Velde (sebastiaan.van.de.velde@ulb.be)
#         Laurine Burdorf         
#
# Example how to use the user interface function of fitprofile to analyse an O2 depth profile. 
#
################################################################################

source("user interface function.R")
source("FLIPPER_plotfunction.R")

#===============================================================================
# Test dataset of classical porewater data
#===============================================================================

load("datasets to test/test.data.1.Rdata")

test <- FLIPPER.func(input=as.data.frame(test.data$input),species=test.data$species,E.cte=c(-0.1,0,0.03),
                     env.parms=test.data$env.parms,discrete.parms=list(LBC="conc.down",C.down=test.data$input$C[nrow(test.data$input)]),
                     continuous.parms=list(optimal.window.size="interactive"),
                     method="all")

x11(height=20, width=20)
plot.FLIPPER(test)

test <- FLIPPER.func(input=as.data.frame(test.data$input),species=test.data$species,E.cte=c(-0.1,0,0.03),
                     env.parms=test.data$env.parms,discrete.parms=list(LBC="conc.down",C.down=test.data$input$C[nrow(test.data$input)]),
                     continuous.parms=list(optimal.window.size="interactive"),
                     method="gradient")

x11(height=20, width=20)
plot.FLIPPER(test)

test <- FLIPPER.func(input=as.data.frame(test.data$input),species=test.data$species,E.cte=c(-0.1,0,0.03),
                     env.parms=test.data$env.parms,discrete.parms=list(LBC="conc.down",C.down=test.data$input$C[nrow(test.data$input)]),
                     continuous.parms=list(optimal.window.size="interactive"),
                     method="continuous")

x11(height=20, width=20)
plot.FLIPPER(test)


test <- FLIPPER.func(input=as.data.frame(test.data$input),species=test.data$species,E.cte=c(-0.1,0,0.03),
                     env.parms=test.data$env.parms,discrete.parms=list(LBC="conc.down",C.down=test.data$input$C[nrow(test.data$input)]),
                     continuous.parms=list(optimal.window.size="interactive"),
                     method="discrete")

x11(height=20, width=20)
plot.FLIPPER(test)

#===============================================================================
# Test dataset of microsensor O2 profiling
#===============================================================================

load("datasets to test/O2_testprofiles.RData")

test <- FLIPPER.func(input=as.data.frame(profile.A[-(1:26),]),species=c("O2"),
             discrete.parms=list(LBC= "no.flux"),
             continuous.parms=list(optimal.window.size="interactive"),
             method="all")

test <- FLIPPER.func(input=as.data.frame(profile.A),species=c("O2"),
                     discrete.parms=list(LBC= "no.flux"),
                     continuous.parms=list(optimal.window.size="interactive"),
                     method="all")

x11(height=20, width=60)
plot.FLIPPER(test)


