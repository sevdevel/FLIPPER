################################################################################
# Author: Sebastiaan van de Velde (sevdevel@vub.ac.be)
#         Laurine Burdorf         (Laurine.Burdorf@nioz.nl)
#
# Example how to use the user interface function of fitprofile to analyse an O2 depth profile. 
#
################################################################################

source("user interface function.R")
source("FLIPPER_plotfunction.R")



load("datasets to test/test.data.5.Rdata")

# as.data.frame(summ$input) 

test <- FLIPPER.func(input=as.data.frame(summ$input),species=summ$species,E.cte=c(-0.1,0,0.03),
             env.parms=summ$env.parms[1:5],discrete.parms=list(LBC="conc.down",C.down=summ$input$C[length(summ$input$C)]),
             continuous.parms=list(optimal.window.size="interactive"),
             method="all")

x11(height=20, width=60)
plot.FLIPPER(test)


test <- FLIPPER.func(input=as.data.frame(summ$input),species=summ$species,#E.cte=c(-0.1,0,0.03),
              env.parms=summ$env.parms[1:5],discrete.parms=list(LBC="conc.down",C.down=summ$input$C[length(summ$input$C)]),
              continuous.parms=list(optimal.window.size="interactive"),
              method="gradient")

x11(height=20, width=20)
plot.FLIPPER(test)

test <- FLIPPER.func(input=as.data.frame(summ$input),species=summ$species,#E.cte=c(-0.1,0,0.03),
                     env.parms=summ$env.parms[1:5],discrete.parms=list(LBC="conc.down",C.down=summ$input$C[length(summ$input$C)]),
                     continuous.parms=list(optimal.window.size="interactive"),
                     method="continuous")

x11(height=20, width=20)
plot.FLIPPER(test)



test <- FLIPPER.func(input=as.data.frame(summ$input),species=summ$species,#E.cte=c(-0.1,0,0.03),
                     env.parms=summ$env.parms[1:5],discrete.parms=list(LBC="conc.down",C.down=summ$input$C[length(summ$input$C)]),
                     continuous.parms=list(optimal.window.size="interactive"),
                     method="discrete")

x11(height=20, width=20)
plot.FLIPPER(test)
