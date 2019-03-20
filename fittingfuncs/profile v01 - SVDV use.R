###############################################################################
# Rscript to interpret vertical profiles
# possible to include: ionic drift, advection, irrigation

# Based on script by Lieke Mulder and Karline Soetaert  
# Author: Sebastiaan van de Velde (sevdevel@vub.ac.be)
# Affiliation: VUB, Pleinlaan 2, 1050, Brussel 
###############################################################################

require(tcltk)
require(ReacTran)
require(marelac)
require(FME)
source("C:/Users/SVDV/OneDrive - Vrije Universiteit Brussel/PhD/Modelling/FLIPPER/fittingfuncs/user_input.R")

# =============================================================================
# 1. Basic diagenetic model => calculates C or Prod based on 1 set of C or P
# =============================================================================

Basic.model  <- function (t=0, C, parms, Prod) {
  with (as.list(parms),{
    
    # Selection of which boundary conditions to use for tran.1D
#     print(str(Ds))
#     print(str(por.grid))
    
    if (UBC == "flux.up" & LBC == "no.flux"){
      tran.summ <- tran.1D(C = C, flux.up = flux.up, C.up = C.up, D = Ds,
                           v=v, VF = por.grid, dx = grid, full.output = T)
    }
    if (UBC == "conc.up" & LBC == "no.flux"){
      tran.summ <- tran.1D(C = C, C.up = C.up, D = Ds,v = v, 
                           VF = por.grid, dx = grid, full.output = T)
    }
    if (UBC == "flux.up" & LBC == "conc.down"){
      tran.summ <- tran.1D(C = C, flux.up = flux.up, C.up = C.up, C.down = C.down, 
                           D = Ds, v=v, VF = por.grid, dx = grid,
                           full.output = T)
    }
    if (UBC == "conc.up" & LBC == "conc.down"){
      tran.summ <- tran.1D(C = C, C.up = C.up, C.down = C.down, 
                           D = Ds, v = v, VF = por.grid, dx = grid,
                           full.output = T)
    }
    if (UBC == "flux.up" & LBC == "flux.down"){
      tran.summ <- tran.1D(C = C, flux.up = flux.up, flux.down = flux.down,
                           D = Ds, v = v, VF = por.grid, 
                           dx = grid, full.output = T)
    }
    if (UBC == "conc.up" & LBC == "flux.down"){
      tran.summ <- tran.1D(C = C, C.up = C.up, flux.down = flux.down,
                           D = Ds, v = v, VF = por.grid, dx = grid,
                           full.output = T)
    }
    
    # Definition of transport parameters: diffusion + irrigation
    
    tran <- tran.summ$dC # M L-3 T-1 POREWATer
    
    # Return differential equation and other
    
    return(list(dCdt       = tran + Prod/por.grid$mid, # M L-3 T-1 POREWATER
                production = Prod,                     # M L-3 T-1 SEDIMENT
                dif.flux   = tran.summ$dif.flux,
                flux.up    = tran.summ$flux.up,
                flux.down  = tran.summ$flux.down,
                adv.flux   = tran.summ$adv.flux
    ))
  })
}

# ========================================================================================
# 2. Steady state solution for Production => Solves basic.model using 1 Production for C 
# ========================================================================================

SteadyState.production  <- function(Production, parms) {
  
  #needs production in every layer
  fprod <- approxfun(Production, rule = 2, method = "constant") 
  
  # needs upper boundary for every layer=> see Zup 
  Prod                       <- fprod(parms$grid$x.mid) 
  Prod[parms$grid$x.mid < 0] <- 0
  
  steady.1D(y     = rep(0, parms$grid$N), 
            func  = Basic.model, 
            nspec = 1, 
            parms = parms,     
            names = "C", 
            Prod  = Prod, 
            atol  = 1e-8)
}

# ========================================================================================
# 3. Costfunction => Calculates SSR en SST between model (2) and data
# ========================================================================================

Cost.func <- function(estimate, depth.up, DATA = DATA, parms, turn=1) {
  
  # Run model with current parameters
  
  Production <- cbind(depth.up,estimate)
  model      <- cbind(x = parms$grid$x.mid, 
                      C = SteadyState.production(Production,parms)$y )
  
  # Calculate cost only from model fit
  
  Cost  <- modCost(model, DATA, x = "x") 
  if (turn==1)  return(Cost)
  else return(model)
}

# ========================================================================================
# 4. Fitting function => Calculates which P gives lowest Cost
# ========================================================================================

#bnd.up <- depth.up(n=n,z=data$x)

Fitting.func <-  function(bnd.up, DATA, guess = NULL , parms, ...){ 
  
  pp <- rep(1, length(bnd.up))
  BB <- modFit(Cost.func, pp, depth.up = bnd.up, DATA = DATA, parms = parms)
  
  return(BB)
}

# ========================================================
# auxiliary functions
# ========================================================

# ============================================================================
# depth of upper boundary of every layer; only works with equidistant layers
# ============================================================================

depth.up <- function(n,z){  
  diff <- (max(z) - min(z)) / n
  c(1:n) * diff - diff
}


# ========================================================
# F-test for PROFILE
# Laurine
# ========================================================

Ftest.profile <- function(nzones1,nzones2, SSE1, SSE2, N){
  
  testF <- ((SSE1 - SSE2)/(nzones2-nzones1)) / (SSE2 / (N - nzones2))
  p     <- pf(testF, nzones2 - nzones1, N - nzones2, lower.tail=F)
  
  return(p)
}

# ========================================================
# Interactive question; adapted from http://www.r-bloggers.com/user-input-using-tcltk/
# Laurine
# ========================================================
PROFILE.interactive <- function(maxzone=10, guess=4){
  
  varEntryDialog(vars=c('Number of Zones'), 
                 labels=paste("Choose number of zones between 1 and",maxzone,"\n PROFILE proposes:",guess), fun=c(
    function(x) {
      x <- as.integer(x)
      if(x >= 0 & x <= maxzone) {
        return(x)
      } else {
        stop(paste("\nPlease enter a number between 1 and",maxzone))
      }
    } ))
}

#========================================================
# Plot production, concentration and fit line
#========================================================

plot.production <- function(depth, conc, modelfit, zones, prod, 
                            conc.limits=NULL, prod.limits=NULL, y.limits=NULL, stat=T, 
                            p.value=NA, R.square=NA){
  
  zones <- c(zones, max(depth))
  if(is.null(y.limits))    y.limits <- c(max(depth, na.rm=T)*1.25,min(depth, na.rm = T)*0.75)
  if(is.null(prod.limits)) prod.limits <- c(range(c(prod*1.25,prod*0.75)))
  if(prod.limits[1]>0) prod.limits[1] <- 0
  if(prod.limits[2]<0) prod.limits[2] <- 0
  if(is.null(conc.limits)) conc.limits <- c(range(c(conc*1.25,conc*0.75)))
  
  
  par(new=F)
  plot(x=prod, y=zones[1:(length(zones)-1)], type="n", ylim=y.limits, axes=F, ylab="", xlab="",
       xlim = prod.limits)
  rect(xleft=0, ybottom=zones[2:length(zones)], ytop=zones[1:(length(zones)-1)], xright=prod, 
       col=gray(level=0.95))
  
  #abline(v=0)
  
  axis(3, cex.axis=1.2, lwd=1.5)
  # abline(h=par()$yaxp[2])
  mtext(text="Production", side=3, line=2.5, cex=1, lwd=1.5)
  
  
  par(new=T)
  
  plot(x=conc, y=depth, ylim=y.limits, pch=21, cex=1, bg=gray(level=0.2), xlab="",ylab="", axes=F,
       xlim = conc.limits)
  lines(x=modelfit[,2], y=modelfit[,1], lwd=2, lty=1, col="red")
  
  axis(1, cex.axis=1.2, lwd=1.5, pos=par()$yaxp[1])
  #abline(h=par()$yaxp[1])
  
  mtext(side=1, line=1.5, "Concentration", cex=1)
  axis(2, cex.axis=1.2, lwd=1.5, pos=par()$xaxp[1])
  mtext(side=2, line=1.5, "Depth", cex=1)
  
  abline(h=0, lty=2)  
  
  if(stat==T){
    plotstat <- plot.statistics(depth, conc, modelfit, zones, prod)
    text(x=mean(conc.limits), y=mean(y.limits), 
         paste(" N Zones = ", plotstat$N.zones,
               "\n int.R = ", round(plotstat$int.R,2),
               "\n R-square = ", round(R.square,4),
               "\n p value = ", round(p.value,4)),
         adj=c(0,1))
    
    
  }  
  
}


#========================================================
# Plot output of fit.profile
#========================================================
plot.fitprofile <- function(x){
  plot.production(depth = x$data$x, conc= x$data$C, prod = x$prod$Prod, modelfit = x$modelfit$C,
                  zones = x$prod$depth)
  
}



integrate.profile <- function(depth,zones,prod){
  
  
  int.R <- sum( diff(zones)* prod)
  
  return(int.R)
}

plot.statistics <- function(depth, conc, modelfit, zones, prod){
  
  
  
  return(
    data.frame(
      N.zones  = length(prod),
      int.R    = integrate.profile(depth,zones,prod)
    ))
  
  
}


#========================================================
# Main function 
# Laurine Burdorf
# Based on script Karline Soetaert, Lieke Mulder,
# Sebastiaan van de Velde
#========================================================

fit.profile <- function(input, 
                        parms, 
                        func          = NULL, 
                        i.end         = 12, 
                        initial.zones = NULL,
                        full.output,
                        ...){
  
  # Read depth and concentration from input file 
  obs <- input[,c("x","C")]
  
  # Load function to estimate best-fit 
  if (is.null(func)) func <- Cost.func
  
  # Start first step 
  print("Finding initial configuration" )
  
  #=============================================================================================
  # Step 1: Calculate best fit production and concentration for 1 to user defined maximum zones
  #=============================================================================================
  
  # 1.1 Prepare production, concentration and fit result dataframes
  Conc.fit    <- data.frame(x=parms$grid$x.mid, V2 = NA)
  
  SSR         <- data.frame(zones = seq(1,i.end),
                            ssr   = NA )
  
  Prod.zone   <- list()
  
  # 1.2 Calculates best production fit for 1:i.end zones
  for (i in 1 : i.end){
    
    Zone.up       <- depth.up(n = i, z=c(0, parms$L.down))
    
    Zone.fit      <- Fitting.func(bnd.up = Zone.up, 
                                  DATA   = obs,
                                  parms  = parms)
    
    Prod.zone[[i]] <- data.frame(depth = Zone.up, 
                                 Prod  = Zone.fit$par)
    
    
    Conc.fit[, i+1] <- Cost.func(estimate=Zone.fit$par, depth.up=Zone.up, DATA=obs, parms, turn=2)[,2]
    
    SSR[i,2]        <- Zone.fit$ssr
    
  }
  
  colnames(Conc.fit)[2:(i.end+1)] <-  paste("zone",seq(1,i.end), sep="")
  
  # 1.5 Add R-square to fit dataframe
  SST         <- sum((obs$C-mean(obs$C))^2)
  SSR$Rsquare <- 1-SSR$ssr/SST
  # print(SST)
  # print(SSR)
  #===========================================================================================
  # Step 2: Calculate p-values between zones (cfr. Table 1, Berg et al 1998)
  #===========================================================================================  
  
  if(is.null(initial.zones)){
    
    
    # 2.2 Creation of empty matrix to store p-values 
    #     p-values calculates between all zone; upper triangle matrix
    Zone.table <- matrix(data=NA,nrow=i.end,ncol=i.end)
    
    for(i in 1:(i.end-1)){
      
      for(j in (i+1):i.end){
        Zone.table[i,j] <- Ftest.profile(nzones1 = i, nzones2 = j, 
                                         SSE1 = SSR[i,2], SSE2 = SSR[j,2], N=nrow(obs))
        
      }
    }
    
    
    
    #===========================================================================================
    # Step 3: Guess best number of zones to start lumping; 
    # Highest amount of zones that is a significant better fit compared to zone x - 1
    #===========================================================================================  
    
    # Get diagonal out of p-value table; get the highest amount that is better than x-1
    GUESS.zone <- max(seq(2:i.end)[diag(Zone.table[,2:i.end]) < parms$prob]) + 1
  }else{
    GUESS.zone <- initial.zones
    Zone.table <- NULL}
  
  print(paste("best initial layers", GUESS.zone))
  
  #===========================================================================================
  # Step 4; lumping ; starts with number of zones as defined in GUESS.zone
  # combines adjecent zones with similar production rates to 1 zone
  #=========================================================================================== 
  
  # 4.1 Create empty data.frames to store all information
  
  Conc2.fit         <- data.frame(x=parms$grid$x.mid, V2 = Conc.fit[,GUESS.zone+1])
  
  Prod2.zone        <- list()
  Prod2.zone[[GUESS.zone]]   <- data.frame(depth = Prod.zone[[GUESS.zone]]$depth, 
                                           Prod  = Prod.zone[[GUESS.zone]]$Prod)
  
  SSR2              <- data.frame(zones = seq((GUESS.zone),1),
                                  ssr   = NA )
  
  SSR2[1,]          <- SSR[GUESS.zone,(1:2)]
  
  # Start values for loop
  zone.production   <- Prod.zone[[GUESS.zone]]$Prod
  zone.up           <- Prod.zone[[GUESS.zone]]$depth
  
  # 4.2 Loop that lump together the two most similar adjecent production zones, untill one 
  # zone is left
  # To be skipped if best initial layers = 1
  if(GUESS.zone > 1){
    print("now lumping")
    
    for (i in 1:(GUESS.zone-1)){
      
      # Find two adjecent zones with similar production;
      # Group two adjecent zones
      Prod.diff     <- abs(diff(zone.production))
      zone.up       <- zone.up[order(Prod.diff, decreasing = T)][1:(length(Prod.diff))]
      zone.up       <- sort(zone.up)
      Zone.fit      <- Fitting.func(bnd.up = zone.up, 
                                    DATA   = obs,
                                    parms  = parms)
      
      zone.production <- Zone.fit$par
      
      SSR2[i+1,2]     <- Zone.fit$ssr
      
      Prod2.zone[[length(zone.up)]] <- data.frame(depth = zone.up,
                                                  Prod  = Zone.fit$par) 
      
      Conc2.fit[, i+2] <- Cost.func(Zone.fit$par, zone.up, DATA=obs, parms, turn=2)[,2]
      
    }
    
    # Add Rsquare values // Rename columns concentration dataframe
    SSR2$Rsquare <- 1 - SSR2$ssr/SST
    
    colnames(Conc2.fit)[2:ncol(Conc2.fit)] <-  paste("zone",seq(GUESS.zone,1), sep="")
    
    
    #===========================================================================================
    # Step 5; calculate p-values of F-test comparing x zones and x-1 zones
    #===========================================================================================  
    
    Zone.lump    <- data.frame(Zones   = (GUESS.zone):2,
                               p.value = NA )
    
    
    # Ouput gives the F-test between zones (x) and zones (x-1)
    for (i in 1:(GUESS.zone-1)){
      Zone.lump[i,2] <-  Ftest.profile(nzones1 = SSR2[i+1,1], nzones2 = SSR2[i,1], 
                                       SSE1 = SSR2[i+1,2], SSE2 = SSR2[i,2], N=nrow(obs)) 
    }
    
    
    #===========================================================================================
    # Step 6; decide on best lumping amount
    # Highest amount of zones that is a better fit than x+1?
    # Plot +/- 3 profiles from initial guess
    #===========================================================================================
    
    Significant.zones <- Zone.lump$Zones[Zone.lump$p.value < 0.01]
    
    final.zones       <- max(Significant.zones)
    
    x11(height = 20,width=40)
    par(mfrow=c(2,4))
    for(i in max(c(1,(final.zones-3))): min(c((final.zones+3),GUESS.zone))){
      
      plot.production(depth    = obs[,1], conc = obs[,2], 
                      modelfit = Conc2.fit[,c("x",paste("zone",i,sep=""))],
                      prod     = Prod2.zone[[i]]$Prod, 
                      zones    = Prod2.zone[[i]]$depth, 
                      p.value  = Zone.lump[Zone.lump$Zones==(i+1),2],
                      R.square = SSR2$Rsquare[SSR2$zones==i])
      
      
    }
    
    final.zones <- PROFILE.interactive(maxzone=GUESS.zone, guess=final.zones)[[1]]
    dev.off()
  }
  else{
    final.zones            <- GUESS.zone
    SSR2$Rsquare           <- 1 - SSR2$ssr/SST
    colnames(Conc2.fit)[2] <-  paste("zone",GUESS.zone, sep="")
    Zone.table             <- NULL
  }
  
  modeloutput             <- SteadyState.production(Production = data.frame(depth.up=Prod2.zone[[final.zones]][,1], estimate=Prod2.zone[[final.zones]][,2]), parms=parms)

  profile.output           <- list()
  profile.output$J.dif.up   <- modeloutput$dif.flux[1]
  profile.output$J.adv.up   <- modeloutput$adv.flux[1]
  profile.output$J.dif.down <- modeloutput$dif.flux[length(modeloutput$dif.flux)]
  profile.output$J.adv.down <- modeloutput$adv.flux[length(modeloutput$dif.flux)]

  profile.output$R.vol      <- Prod2.zone[[final.zones]]
  
  profile.output$R.int      <- integrate.profile(zones = c(profile.output$R.vol[,"depth"], max(obs[,1])),
                                                 prod  = c(profile.output$R.vol[,"Prod"]))
    
  profile.output$fit        <- data.frame(x = Conc2.fit[,"x"],
                                        C = Conc2.fit[,paste("zone",final.zones, sep="")])
  
  
  if (full.output == FALSE){ return(profile.output) } else{
      profile.output$J.all      <- data.frame(x.int        = parms$grid$x.int,
                                              J.adv = modeloutput$adv.flux,
                                              J.dif = modeloutput$dif.flux,
                                              J.tot = modeloutput$dif.flux + modeloutput$adv.flux)
      profile.output$SSR      <- SSR2[SSR2$zones==final.zones,"ssr"]
      profile.output$R.square <- SSR2[SSR2$zones==final.zones,"Rsquare"]
      profile.output$Zone.Table <- Zone.table
      
      return(profile.output)
    }
  
}

