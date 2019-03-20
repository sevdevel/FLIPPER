################################################################################
# function name: gradient.fit
#
# Author: Filip meysman (filip.meysmna@nioz.nl)
#         Laurine Burdorf (laurine.burdorf@nioz.nl)
#         Sebastiaan van de Velde (sevdevel@vub.ac.be)
#
# Package dedicated to geochemical flux analysis, based on the gradient near the SWI
################################################################################

##############################################################################
# Regression analysis (public)                                                      
##############################################################################

gradient.fit <- function(profile, x.limits = NA, env.parms, full.output)
{
# This function plots a concentration depth profile, determines the mean
# concentration gradient over a given depth interval by linear regression,
# and calculates the diffusive flux.
#
  #=============================================================================
  # Arguments 
  #=============================================================================
  #   profile   a data.frame object that contains the concentration depth profile. 
  #             This data.frame should contain the columns 
  #                x    = depth  [L]
  #                C    = concentration [M L-3]
  #                por  = porosity [dimensionless]
  #                tort = tortuosity [dimensionless]
  #                E    = electrical field [V m-1] (optional)
  #   x.limits  depth limits of the profile
  #   env.parms a list containing temperature TC [degrees Celcius]
  #                salinity S [PSU] 
  #                pressure P [bar]
  #                diffusion coefficient Dmol [cm2 d-1]
  #                the charge of the ion z [] (e.g. +2 for Ca2+)
  #
  #
  #=============================================================================
  # Value 
  #=============================================================================
  #  A list containing:
  #
  #    diff.flux    diffusive flux value calculated as flux = -Dmol*(por/tort)*fit$slope  
  #    adv.flux     advective flux due to ionic drift (electric field), calculated as flux = Dmol*por*z*F*E.SWI/R/T*C.SWI
  #    tot.flux     total flux value, calculated as diff.flux + adv.flux
  #    fit     list containing fitting parameters   
  #=============================================================================
  
  # Set limits of profile to be plotted
  if(is.na(x.limits[1])) x.limits <- sort(range(profile$x),decreasing=TRUE)
  selection <- (profile$x >= min(x.limits))&(profile$x <= max(x.limits))
  Depth <- profile$x[selection]
  Value <- profile$C[selection]
  
  # Satisfaction loop
  satisfied <- FALSE
  while (!satisfied){
    
    x11(width=20, height=20)
    par(mfrow=c(1,1))
    
    # Plot depth profile
    plot(Value, Depth,
      ylab="Depth", xlab="", ylim = x.limits, axes=FALSE, 
      frame.plot=FALSE, type ="p", col="red", pch=15,
      main = "")
    
    axis(pos=par()$xaxp[1], side=2,lwd=1)
    abline(h=0, lwd=1, lty=1)
    axis(pos=par()$yaxp[2], side=3,lwd=1)
    
    # Retrieve marker coordinates
    xdot <- (max(Value) + min(Value))/2
    ydot <- 0.9*max(Depth)
    text(x=xdot,y=ydot,pos=1,labels="Click on start point and end point")

    start.mark <- identifyMark(x=Value, y=Depth, labels = c("point 1"),n = 1,col="blue",pch=19)
    end.mark <- identifyMark(x=Value, y=Depth, labels = c("point 2"),n = 1,col="blue",pch=19)
    
    # Calculate linear lit

    fit    <- grad.fit(x=Depth[start.mark:end.mark],C=Value[start.mark:end.mark])
    fit$n  <- length(Value[start.mark:end.mark])
    lines(fit$C,fit$x,lwd="2",col="blue")
    
    satisfied <- plot.and.ask.for.OK(Value,Depth, xdot=c((max(Value) + min(Value))/2,(max(Value) + min(Value))/2), 
                                     ydot=c(0.5*max(Depth),0.6*max(Depth)))
  }
  dev.off()
  
  if (!is.null(profile$E)){
  Faraday <- as.numeric(Constants$F[1])      # Faraday constant [C mol-1]
  R       <- as.numeric(Constants$gasCt2[1]) # Molar gas constant [J mol-1 K-1]
  selection <- min(which(profile$x >= 0))
  C.SWI    <- profile$C[selection]
  E.SWI    <- profile$E[selection]
  por.SWI  <- profile$por[selection]
  tort.SWI <- 1-2*log(por.SWI) 
  v          <- env.parms$Dmol/tort.SWI*env.parms$z*Faraday*E.SWI/(R*(env.parms$TC+273.15)) # m d-1 
  adv.flux   <- por.SWI*v*C.SWI                                                             # mmol m-2 d-1
  } else {adv.flux <- 0}
  
  por     <- mean(profile$por[start.mark:end.mark])
  tort    <- mean(profile$tort[start.mark:end.mark])
  
  diff.flux  <- - por*(env.parms$Dmol/tort)*fit$slope # mmol m-2 d-1
  tot.flux   <- diff.flux + adv.flux                  # mmol m-2 d-1
    
  if (full.output == FALSE){
    return(list(J.dif.up = diff.flux, J.adv.up = adv.flux, R.int = -diff.flux, fit = fit))
  } else {
  return(list(J.dif.up = diff.flux, J.adv.up = adv.flux, R.int = -diff.flux, fit = fit))
  }
}

##############################################################################
# Gradient fitting procedure (private)
##############################################################################

grad.fit <- function (x,C,x.new=NULL)
{
  a <- lsfit(x = x, y = C)$coef["X"][[1]]
  b <- lsfit(x = x, y = C)$coef["Intercept"][[1]]
  if(!is.null(x.new)) {y = a*x.new + b; x <- x.new}
  if(is.null(x.new)) {y = a*x + b}
  return(list(x=x,C=y,slope=a,intercept=b))
}  

##############################################################################
# Identification of points on a graph (private functions)                                       #
##############################################################################

# Indentify break point 

click.break.point <- function(x,y,xlab="",ylab="",ylim=c(0.5,2), message=NULL)
{
  # Initialisation  
  
  satisfied <- FALSE
  while (!satisfied){
    
    plot(x, y,
      ylab=ylab, xlab=xlab, ylim = ylim,  
      type ="p", col="black", pch=15, main=message)
    
    mark <- identifyMark(x, y, labels = c("mark"),n = 1,col="blue",pch=19)
    
    satisfied <- plot.and.ask.for.OK(x,y,xdot=0.75*c(mean(x),mean(x)),ydot=c(0.8*ylim[2],0.9*ylim[2]))
  }
  
  return(mark)
}

identifyMark <- function(x, y=NULL, n=length(x), col="red", pch=19, ...)
{
  xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
  sel <- rep(FALSE, length(x)); res <- integer(0)
  while(sum(sel) < n) {
    ans <- identify(x[!sel], y[!sel], n=1, plot=TRUE, col=col,...)
    if(!length(ans)) break
    ans <- which(!sel)[ans]
    points(x[ans], y[ans], col = col, pch=pch)
    sel[ans] <- TRUE
    res <- c(res, ans)
  }
  res
}

# Check whether satisfied
plot.and.ask.for.OK <- function(x,y,xdot=NULL,ydot=NULL)
{
  satisfied <- FALSE
  
  if(is.null(xdot))
  {
    xdot <- vector(length=2)
    xdot[1] <- min(x)+0.8*(max(x)-min(x))
    xdot[2] <- xdot[1]
  }
    
    if(is.null(ydot))
  {
    ydot <- vector(length=2)
    ydot[1] <- min(y)+ 0.8*(max(y)-min(y))
    ydot[2] <- min(y)+ 0.9*(max(y)-min(y))
  }
      
  points(x=xdot[1], y=ydot[1],pch=21,col="blue")
  points(x=xdot[2], y=ydot[2],pch=21,col="red")
  text(x=xdot[1],y=ydot[1],labels=c("OK"),col="blue", pos=4)
  text(x=xdot[2],y=ydot[2],labels=c("REDO"),col="red", pos=4)
  mark <- identify(x=xdot, y=ydot, n = 1, plot = FALSE)
  if (mark == 1) points(x=xdot[mark], y=ydot[mark], pch=19, col="blue")
  if (mark == 2) points(x=xdot[mark], y=ydot[mark], pch=19, col="red")
  if (mark == 1) satisfied <- TRUE 
  
  return(satisfied)
}


# Click on a point
click.point <- function(profile, lim = NA, method.label = NA)
{
  # Initialisation  

  if(is.na(lim[1])) lim <- sort(range(profile$x),T)
  
  Depth <- profile$x[profile$x <= max(lim)]
  Value <- profile$C[profile$x <= max(lim)]
  
  satisfied <- FALSE
  while (!satisfied){
    
    plot(Value, Depth,
         ylab=(""), xlab="", ylim = lim, axes=FALSE, 
         frame.plot=FALSE, type ="p", col="red", pch=15,
         main = method.label)
    
    axis(pos=par()$xaxp[1], side=2,lwd=1)
    abline(h=0, lwd=1, lty=1)
    axis(pos=par()$yaxp[2], side=3,lwd=1)
    
    
    xdot <- ( max(Value) + min(Value) ) / 2
    ydot <- max(Depth)
    
    
    text(xdot,0.8*ydot,labels="Indicate mark")
    #text(xdot,ydot,labels="")
    mark <- identifyMark(x=Value, y=Depth, labels = c("mark"),n = 1,col="blue",pch=19)
    
    satisfied <- plot.and.ask.for.OK(Value,Depth)
  }
  
  return(mark)
}

