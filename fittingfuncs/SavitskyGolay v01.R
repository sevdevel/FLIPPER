################################################################################
# Package name: fitProfile
#
# Author: Filip meysman (filip.meysmna@nioz.nl)
#
# Package dedicated to geochemical flux and rate analysis, i.e, the robust
# numerical extraction of fluxes and rates
# from concentration depth profiles in aquatic sediments. 
################################################################################

require(signal)
require(fractaldim)

##############################################################################
# Savitzky Golay filter with selection of the optimal filter window based on
# fractal dimension analysis (public)
##############################################################################

sgolayfilt.optim <- function(
  x,
  p = 3, 
  n = NULL,
  m = 0,
  ts = 1,
  optimal.window.detection = c("interactive","automated"), 
  min.n = p - p%%2,
  max.n = length(x)%/%2-1,
  show.graphics = FALSE,
  save.graphics = FALSE,
  graphics.title = NULL,
  full.output.extra = FALSE
)
  
{

  # Issue warning when filter length is not null
  if (!is.null(n)){

    x.filtered <- sgolayfilt(x=x, p=p, n=n, m=m, ts=ts) 
    warning("sgolayfilt.optim: n is not NULL, classical sgolay filtering is applied")
    return (x.filtered)
  }

  # Preparation of data structures 
  
  SG <- matrix(nrow=length(x),ncol=max.n) # matrix containing filtered SG signal for each window n  
  FD <- rep(0,max.n)
  
  # Optimal filter window is determined based on fractal dimension analysis, see
  # "Estimators of Fractal Dimension: Assessing the Roughness of Time Series
  # and Spatial Data",  Gneiting T. et al (2011)
  #
  # sgolayfilt from package "signal"
  # fd.get and fd.estimate from package "fractaldim"
  
  methods <- c("madogram", "variogram", "hallwood", "boxcount","wavelet")
  for (i in min.n:max.n) {
    SG[,i] <- sgolayfilt(x=x, p=p, n=(2*i+1), m=m, ts=ts) 
    FD[i] <- fd.get(fd.estimate(SG[,i], methods = methods),methods[1])$fd
  }  
  
  if (optimal.window.detection[1] == "interactive")
  {
    # Interactive selection
    x11()
    n <- click.break.point(x=min.n:max.n, y=FD[min.n:max.n],
      xlab="filter window",ylab="Fractal dimension",message=paste("Indicate break point","graphics.title"))
    if (save.graphics) savePlot(filename = graphics.title,type ="jpeg")   
    if (!show.graphics) dev.off(which = dev.cur())   
    
  } else {
    # Automated selection
    if (max(FD) < 1.05){
      n <- min.n
      # plot output
      x11()
      plot(min.n:max.n,FD[min.n:max.n],xlab="filter window",ylab="Fractal dimension",main=graphics.title)  
      points(x=n,y=FD[n],pch=20,lwd=4,col="blue")
      if (save.graphics) savePlot(filename = graphics.title,type ="jpeg")   
      if (!show.graphics) dev.off(which = dev.cur())   
    }else{
      test <- nls(y ~ f.fdts(x,a,b),
        data = data.frame(x=min.n:max.n,y=FD[min.n:max.n]),
        start = list(a = 10, b = 1))
      ind <- which(predict(test) == min(predict(test)))[1]
      n <- (min.n:max.n)[ind]
      # plot output
      x11()
      plot(min.n:max.n,FD[min.n:max.n],xlab="filter window",ylab="Fractal dimension",main=graphics.title)  
      lines(min.n:max.n,predict(test),col="red",lwd=2)
      points(x=n,y=FD[n],pch=20,lwd=4,col="blue")
      if (save.graphics) savePlot(filename = graphics.title,type ="jpeg")   
      if (!show.graphics) dev.off(which = dev.cur())   
    }
  }
  
  # Return statement
  if (full.output.extra)  output <- list(y=SG[,n],n=n,SG=SG,FD=FD)
  else output<- SG[,n]
  
  return(output)
}


##############################################################################
# function f.fdts (private)
##############################################################################

f.fdts <- function(x,a=10,b=1,n=11) 
{
  res <- pmax(b*(1-x/a),0)
  # res <- 1 + res
  res <- 1+smooth.spline(x=res,df=n)$y
  return(res)
} 

# plot(0:35,f.fdts(0:35))

##############################################################################
# Savitzky Golay filter extended with boundary conditions (public)
##############################################################################

sgolayfilt.ext <- function(
  x,
  p = 3, 
  n,
  m = 0,
  ts = 1,
  bnd.upper = 1,
  bnd.lower = 1)

#=============================================================================
# Arguments 
#=============================================================================

#  p : order of the polynomial that is fitted (typically 2 or 3) 
#  bnd.upper : integer determining filter behaviour at the upstream (upper, left) boundary. Value = 1: no constraint on flux (default value). Value = 2: constant flux imposed. Value = 3: zero flux imposed.   
#  bnd.lower : integer determining filter behaviour at the downstream (lower or right) boundary. Value = 1: no constraint on flux (default value). Value = 2: constant flux imposed. Value = 3: zero flux imposed.   

{
  # Preparation of data series
  
  if (n%%2 == 1) stop("filter length n must be an odd integer")
  n.fw <- n%/%2 # filter window

  aux.x <- x
  selection <- 1:length(x)
  
  # Upstream boundary adaptation
  
  if (bnd.upper==2) # d2C_dx2 = 0; constant flux boundary 
  {
    new.x <- grad.fit(x=ts*1:n,C=x[1:n.fw],x.new=(-ts*(n.fw-1):0))$C
    aux.x <- c(new.x,aux.x)
    selection <- 1:length(x)+n.fw 
  } 
  
  if (bnd.upper==3) # dC_dx = 0 and d2C_dx2 = 0; zero flux boundary
  {
    aux.x <- c(rep(aux.x[1],n.fw),aux.x)
    selection <- 1:length(x)+n.fw 
  } 
  
  if (bnd.upper==4) # d2C_dx2 = 0;  constant flux, mirrored profile
  {
    x0 <- aux.x[1] + 0.5*(aux.x[1]-aux.x[2])
    aux.x <- c(2*x0 - aux.x[n.fw:1],aux.x)
    selection <- 1:length(x)+n.fw 
  } 
  
  # Lower Upper boundary adaptation
  
  if (bnd.lower==2) # d2C_dx2 = 0; constant flux boundary 
  {
    new.x <- grad.fit(x=ts*(n.data-n.fw+1):n.data,C=x[(n.data-n.fw+1):n.data],x.new=(ts*((n.data+1):(n.data+n.fw))))$C
    aux.x <- c(aux.x,new.x)
  } 
  
  if (bnd.lower==3) # dC_dx = 0 and d2C_dx2 = 0; zero flux boundary
  {
    aux.x <- c(aux.x,rep(aux.x[length(aux.x)],n.fw))
  } 
  
  if (bnd.lower==4) # d2C_dx2 = 0;  constant flux, mirrored profile
  {
    xn <- aux.x[length(aux.x)] + 0.5*(aux.x[length(aux.x)]-aux.x[length(aux.x)-1])
    aux.x <- c(2*xn - aux.x[length(aux.x):(length(aux.x)-n.fw+1)],aux.x)
  } 
  
  # Lower Upper boundary adaptation

  x.filtered <- sgolayfilt(x=aux.x, p=p, n=n, m=m, ts=ts) 

  if (full.output)  output <- list(y=x.filtered[selection],y.aux=aux.x)
  else output <- x.filtered[selection]
  
  return(output)
  
}

##############################################################################
# Savitzky Golay analysis (public)
#
# Function that analyses the data series representing the one-dimensional 
# concentration depth profile C(x) of a geochemical variable in an aquatic sediment. 
# Based on Savitzky Golay filtering, the function calculates at each point 
#   the smoothed depth concentration profile
#   the first order derivative (dC_dx),
#   the second order derivative (d2C_dx2), 
#   the diffusive flux (J = -por*Ds*dC_dx),
#   the production rate (R = -por*Ds*d2C_dx2- (por*dDs_dx+dDs*dpor_dx)*dC_dx),
# where por denotes the porosity and DS the effective diffusion coefficient. 
#
##############################################################################

SavGolay.analysis <- function(profile,
  p=3,
  n.C=NULL,n.J=NULL,n.R=NULL,
  bnd.upper = 1,
  bnd.lower = 1,
  n.uniform = FALSE,
  optimal.window.size = c("interactive","automated"), 
  min.n = p - p%%2,
  max.n = nrow(profile)%/%2-1,
  keep.graphics = FALSE,
  full.output)
{
  
  #=============================================================================
  # Arguments 
  #=============================================================================
  #
  #  profile : a data.frame object that contains the columns 
  #              x = depth  [m]
  #              C  = concentration [mmol m-3]
  #              por = porosity [dimensionless]
  #              V = advective velocity [m d-1]               
  #              Ds = effective diffusion coefficient [m2 d-1]               
  #              dpor_dx = derivative of porosity [m-1]
  #              dv_dx = derivative of advective velocity [d-1]
  #              dDs_dx = derivative of effective diffusion coefficient [m d-1]
  #
  #  p : filter order, i.e,  order of the polynomial that is fitted (typically 2 or 3) 
  #  n.C : window size used in filtering the concentration C, representing the number of data points to the
  #        left and right of the data midpoint (hence the filter length  = 2*n.C+1)
  #  n.J : window size used in filtering the first order derivative dC_dx, representing the number of data points to the
  #        left and right of the data midpoint (hence the filter length  = 2*n.J+1)
  #  n.R : window size used in filtering the second order derivative d2C_dx2, representing the number of data points to the
  #        left and right of the data midpoint (hence the filter length  = 2*n.R+1)
  #  bnd.upper : integer determining filter behaviour at the upstream (upper, left) boundary. Value = 1: no constraint on flux (default value). Value = 2: constant flux imposed. Value = 3: zero flux imposed.   
  #  bnd.lower : integer determining filter behaviour at the downstream (lower or right) boundary. Value = 1: no constraint on flux (default value). Value = 2: constant flux imposed. Value = 3: zero flux imposed.   
  #  n.uniform : logical, if TRUE then n.C, n.J, and n.R are all set uniform to max(n.C,n.J,n.R)
  #  min.n : minimum value of the window size, used when scanning the optimal window size 
  #  max.n : maximum value of the window size, used when scanning the optimal window size
  #  full.output : logical, TRUE gives all output, FALSE gives selected output (without info, see below)
  #
  #=============================================================================
  # Value 
  #=============================================================================
  # A list containing following components:
  #
  # out : data.frame with all calculated output, which 
  #         x coordinates
  #         C smoothed concentration
  #         J diffusive flux
  #         R production rate
  #         dC_dx first order derivative
  #         d2C_dx2 second order derivative
  #
  # filter.windows : a list containing
  #         n.C filter window for C
  #         n.J filter window for dC_dx 
  #         n.R filter window for d2C_dx2 
  #
  # h : filter distance calculated from tha data series
  # J.up : diffusive flux at the upstream boundary
  # J.down : diffusive flux at the downstream boundary
  # R.int : integrated production rate over the whole domain
  #
  # info: auxiliary matrices and vectors generated during the selection of optimal filter windows, a list containing: 
  #          SG.C : matrix containing the smoothed concentration for all scanned window sizes  
  #          SG.J : matrix containing the first order derivative for all scanned window sizes  
  #          SG.R : matrix containing the second order derivative for all scanned window sizes  
  #          C.fdts : vector containing the fractal dimensions of SG.C for all scanned window sizes 
  #          J.fdts : vector containing the fractal dimensions of SG.J for all scanned window sizes 
  #          R.fdts : vector containing the fractal dimensions of SG.R for all scanned window sizes 
  #         
  #=============================================================================
  # Filter distance
  #=============================================================================
  
  h <- unique(signif(diff(profile$x),digits=2))
  if (length(h)>1) stop("data points are not at equidistant depth intervals")
  
  #=============================================================================
  # Determine the optimal filter windows when they're not specified 
  #=============================================================================

  # Optimal filter window is determined based on fractal dimension analysis, see
  # "Estimators of Fractal Dimension: Assessing the Roughness of Time Series
  # and Spatial Data",  Gneiting T. et al (2011)
  
  methods <- c("madogram", "variogram", "hallwood", "boxcount","wavelet")
  N.data <- length(profile$C)

  start <- min.n
  end <- max.n

  info <- list()
  
  if (is.null(n.C))
  {
    C.mat <- matrix(nrow=N.data,ncol=end)
    C.fdts <- rep(0,end)
    for (i in start:end) {
      C.mat[,i] <- sgolayfilt(x=profile$C, p=p, n=(2*i+1), m=0, ts=h) 
      C.fdts[i] <- fd.get(fd.estimate(C.mat[,i], methods = methods),methods[1])$fd
    }  
    
    if (optimal.window.size[1] == "interactive")
    {
      x11()
      # Interactive selection
      n.C <- click.break.point(x=start:end, y=C.fdts[start:end],message="C plot: Indicate break point")

    } else {
      
      x11()
      # Automated selection
      if (max(C.fdts) < 1.05){
        n.C <- start
        plot(start:end,C.fdts[start:end])  
      
      }else{
        test <- nls(y ~ f.fdts(x,a,b),
          data = data.frame(x=start:end,y=C.fdts[start:end]),
          start = list(a = 10, b = 1))
        ind <- which(predict(test) == min(predict(test)))[1]
        n.C <- (start:end)[ind]
#        n.C <- round(coef(test)["a"][[1]],digits=0)
        plot(start:end,C.fdts[start:end])  
        lines(start:end,predict(test),col="red")
        points(x=n.C,y=1,pch=20,lwd=4)
        
      }
        
    }
    
    info$SG.C <- C.mat
    info$C.fdts <- C.fdts
    if (!keep.graphics) dev.off(which = dev.cur())   
  }
  
  if (is.null(n.J))
  {  
    J.mat <- matrix(nrow=N.data,ncol=end)
    J.fdts <- rep(0,end)
    for (i in start:end)
    {
      J.mat[,i] <- sgolayfilt(x=profile$C, p=p, n=(2*i+1), m=1, ts=h) 
      J.fdts[i] <- fd.get(fd.estimate(J.mat[,i], methods = methods),methods[1])$fd
    }
    
    if (optimal.window.size[1] == "interactive")
    {
      x11()
      # Interactive selection
      n.J <- click.break.point(x=start:end, y=J.fdts[start:end],message="J plot: Indicate break point")
      
    } else {
      
      x11()
      # Automated selection
      if (max(J.fdts) < 1.05){
        n.J <- start
        plot(start:end,J.fdts[start:end])  
        
      }else{
        test <- nls(y ~ f.fdts(x,a,b,n=15),
          data = data.frame(x=start:end,y=J.fdts[start:end]),
          start = list(a = 10, b = 1))
        ind <- which(predict(test) == min(predict(test)))[1]
        n.J <- (start:end)[ind]
        #      n.J <- round(coef(test)["a"][[1]],digits=0)
        plot(start:end,J.fdts[start:end])  
        lines(start:end,predict(test),col="red")
        points(x=n.J,y=1,pch=20,lwd=4)
      }
    }

    info$SG.J <- J.mat
    info$J.fdts <- J.fdts
    if (!keep.graphics) dev.off(which = dev.cur())   
  }
  
  if (is.null(n.R))
  {
    
    R.mat <- matrix(nrow=N.data,ncol=end)
    R.fdts <- rep(0,end)
    for (i in start:end)
    {
      R.mat[,i] <- sgolayfilt(x=profile$C, p=p, n=(2*i+1), m=2, ts=h) 
      R.fdts[i] <- fd.get(fd.estimate(R.mat[,i], methods = methods),methods[1])$fd
    }
    
    if (optimal.window.size[1] == "interactive")
    {
      x11()
      # Interactive selection
    n.R <- click.break.point(x=start:end, y=R.fdts[start:end],message="R plot: Indicate break point")
      
    } else {
      
      # Automated selection
      x11()
      if (max(R.fdts) < 1.05){
        n.R <- start
        plot(start:end,R.fdts[start:end])  
      } else {
        test <- nls(y ~ f.fdts(x,a,b,n),
          data = data.frame(x=start:end,y=R.fdts[start:end]),
          start = list(a = 10, b = 1, n = 11))
        ind <- which(predict(test) == min(predict(test)))[1]
        n.R <- (start:end)[ind]
        #      n.R <- round(coef(test)["a"][[1]],digits=0)
        plot(start:end,R.fdts[start:end])  
        lines(start:end,predict(test),col="red")
        points(x=n.R,y=1,pch=20,lwd=4)
      }
    }
    
    info$SG.R <- R.mat
    info$R.fdts <- R.fdts
    if (!keep.graphics) dev.off(which = dev.cur())   
  }
  
  if (n.uniform) n.C <- n.J <- n.R <- max(c(n.C,n.J,n.R))

  
  #=============================================================================
  # Adapt boundaries if needed
  #=============================================================================

  # Preparation of data series

  aux.profile <- list(x=profile$x,C=profile$C)
  N  <- length(profile$C)
  selection <- 1:N
  
  # Upper boundary adaptation

  if (bnd.upper==2) # d2C_dx2 = 0; constant flux boundary 
  {
    n.fit <- n.C
    new.x <- (aux.profile$x[1]-h*1:N)[N:1]
    new.C <- grad.fit(x=aux.profile$x[1:n.fit],C=aux.profile$C[1:n.fit],x.new=new.x)$C
    aux.profile$x <- c(new.x,aux.profile$x)
    aux.profile$C <- c(new.C,aux.profile$C)
    selection <- (N+1):(2*N) 
  } 

  if (bnd.upper==3) # dC_dx = 0 and d2C_dx2 = 0; zero flux boundary
  {
    aux.profile$x <- c(aux.profile$x[1]-h*1:N,aux.profile$x)
    aux.profile$C <- c(rep(aux.profile$C[1],N),aux.profile$C)
    selection <- (N+1):(2*N) 
  } 
  
  if (bnd.upper==4) # d2C_dx2 = 0;  constant flux, mirrored profile
  {
    C0 <- aux.profile$C[1] + 0.5*(aux.profile$C[1]-aux.profile$C[2])
    aux.profile$x <- c((aux.profile$x[1]-h*1:N)[N:1],aux.profile$x)
    aux.profile$C <- c(2*C0 - aux.profile$C[N:1],aux.profile$C)
    selection <- (N+1):(2*N) 
  } 
  
  # Lower boundary adaptation

  N.data  <- length(aux.profile$C)
  
  if (bnd.lower==2) # d2C_dx2 = 0; constant flux boundary 
  {
    n.fit <- n.C
    new.x <- aux.profile$x[N.data]+h*1:N
    new.C <- grad.fit(x=aux.profile$x[(N.data-n.fit+1):(N.data)],
      C=aux.profile$C[(N.data-n.fit+1):(N.data)],x.new=new.x)$C
    aux.profile$x <- c(aux.profile$x,new.x)
    aux.profile$C <- c(aux.profile$C,new.C)
  } 

  if (bnd.lower==3)# dC_dx = 0 and d2C_dx2 = 0; zero flux boundary
  {
    aux.profile$x <- c(aux.profile$x,aux.profile$x[N.data]+h*1:N)
    aux.profile$C <- c(aux.profile$C,rep(aux.profile$C[N.data],N))
  } 
  
  if (bnd.lower==4) # d2C_dx2 = 0;  constant flux, mirrored profile
  {
    CN <- aux.profile$C[N.data] + 0.5*(aux.profile$C[N.data]-aux.profile$C[N.data-1])
    aux.profile$x <- c(aux.profile$x,aux.profile$x[N.data]+h*1:N)
    aux.profile$C <- c(aux.profile$C,2*CN - aux.profile$C[N.data:(N.data-N+1)])
  } 

  #=============================================================================
  # Actual SG analysis
  #=============================================================================

  C <- sgolayfilt(x=aux.profile$C, p=p, n=2*n.C+1, m=0, ts=h) 
  dC_dx <- sgolayfilt(x=aux.profile$C, p=p, n=2*n.J+1, m=1, ts=h) 
  d2C_dx2 <- sgolayfilt(x=aux.profile$C, p=p, n=2*n.R+1, m=2, ts=h) 

  #=============================================================================
  # Return statement
  #=============================================================================

if (full.output == TRUE){
  out <- data.frame(x=profile$x,
                    C      = C[selection],
                    J.tot  = -profile$por*profile$Ds*dC_dx[selection] + profile$por*profile$v*C[selection],
                    J.dif  = -profile$por*profile$Ds*dC_dx[selection],
                    J.adv  = profile$por*profile$v*C[selection],
                    R.vol = -profile$por*profile$Ds*d2C_dx2[selection] +
                      (profile$por*profile$v - profile$por*profile$dDs_dx - profile$dDs*profile$dpor_dx)*dC_dx[selection] + 
                      (profile$por*profile$dv_dx + profile$dpor_dx*profile$v)*C[selection],
                    dC_dx = dC_dx[selection],
                    d2C_dx2 = d2C_dx2[selection])
  
  output <- list()
  
  output$J.dif.up    <- out$J.dif[1]
  output$J.adv.up    <- out$J.adv[1]
  output$J.dif.down  <- out$J.dif[length(out$J.dif)]
  output$J.adv.down  <- out$J.adv[length(out$J.adv)]
  
  N.out <- length(out$R.vol)
  output$R.int       <- sum(diff(out$x)*0.5*(out$R.vol[2:N.out]+out$R.vol[1:(N.out-1)]))
  
  output$overview <- out
  
  output$filter.windows <- list(n.C=n.C,n.J=n.J,n.R=n.R)
  output$info <- info
  
  return(output)
} else{

  out <- data.frame(x=profile$x,
                    C      = C[selection],
                    J.tot  = -profile$por*profile$Ds*dC_dx[selection] + profile$por*profile$v*C[selection],
                    J.dif  = -profile$por*profile$Ds*dC_dx[selection],
                    J.adv  = profile$por*profile$v*C[selection],
                    R.vol = -profile$por*profile$Ds*d2C_dx2[selection] +
                      (profile$por*profile$v - profile$por*profile$dDs_dx - profile$dDs*profile$dpor_dx)*dC_dx[selection] + 
                      (profile$por*profile$dv_dx + profile$dpor_dx*profile$v)*C[selection])
  
  output <- list()
    
  output$J.dif.up    <- out$J.dif[1]
  output$J.adv.up    <- out$J.adv[1]
  output$J.dif.down  <- out$J.dif[length(out$J.dif)]
  output$J.adv.down  <- out$J.adv[length(out$J.adv)]
  
  N.out <- length(out$R.vol)
  output$R.int       <- sum(diff(out$x)*0.5*(out$R.vol[2:N.out]+out$R.vol[1:(N.out-1)]))
  
  output$overview <- out
  return(output)
}
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

