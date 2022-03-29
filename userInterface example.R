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

# uniform filter window of 5 (n=2)
continuous.parms <- list(n.C=2,n.J=2,n.R=2)
test <- FLIPPER.func(input=input.profile,species=c("O2"),por.cte=por.cte, 
                     env.parms=env.parms,
                     continuous.parms=continuous.parms,
                     method="continuous")
plot.FLIPPER(test)

# uniform filter window of 21 (n=15)
continuous.parms <- list(n.C=15,n.J=15,n.R=15)
test <- FLIPPER.func(input=input.profile,species=c("O2"),por.cte=por.cte, 
                     env.parms=env.parms,
                     continuous.parms=continuous.parms,
                     method="continuous")
plot.FLIPPER(test)

# uniform filter window of 81 (n=40)
continuous.parms <- NULL#list(n.C=40,n.J=40,n.R=40)
test <- FLIPPER.func(input=input.profile,species=c("O2"),por.cte=por.cte, 
                     env.parms=env.parms,
                     continuous.parms=continuous.parms,
                     method="continuous")
plot.FLIPPER(test)

# Run FLIPPER continuous and interactive
continuous.parms <- list()
test <- FLIPPER.func(input=input.profile,species=c("O2"),por.cte=por.cte, 
                     env.parms=env.parms,
                     continuous.parms=continuous.parms,
                     method="continuous")
plot.FLIPPER(test)




output.example <- list()
#output.example$gradient <- test
#output.example$discrete <- test


save(output.example,file="datasets to test/output.example.Rdata")

input.profile <- as.data.frame(cbind(profile.A$C*1e3,profile.A$x*1e-2,profile.A$por))
colnames(input.profile) <- c('C','x','por')
input.profile <- as.data.frame(cbind(profile.A$C*1e3,profile.A$x*1e-2,profile.A$por))
colnames(input.profile) <- c('C','x','por')

test <- FLIPPER.func(input=input.profile,species=c("O2"),
             discrete.parms=list(LBC= "no.flux"),
             continuous.parms=list(optimal.window.size="interactive"),
             method="all")

x11(height=20, width=20)
plot.FLIPPER(test)

input.profile <- as.data.frame(cbind(profile.B$C*1e3,profile.B$x*1e-2,profile.B$por))
colnames(input.profile) <- c('C','x','por')

test <- FLIPPER.func(input=input.profile,species=c("O2"),
                     discrete.parms=list(LBC= "no.flux"),
                     continuous.parms=list(optimal.window.size="interactive"),
                     method="all")

x11(height=20, width=60)
plot.FLIPPER(test)


input.profile <- as.data.frame(cbind(profile.C$C*1e3,profile.C$x*1e-2,profile.C$por))
colnames(input.profile) <- c('C','x','por')

test <- FLIPPER.func(input=input.profile,species=c("O2"),
                     discrete.parms=list(LBC= "no.flux"),
                     continuous.parms=list(optimal.window.size="interactive",n.uniform=TRUE),
                     method="all")

x11(height=20, width=60)
plot.FLIPPER(test)

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