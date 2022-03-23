#=========================================================
# easy to use interface function for package FLIPPER
# written by Laurine Burdorf (Laurine.Burdorf@nioz.nl)
#            Sebastiaan van de Velde (sevdevel@vub.ac.be)
# version 02 - data: 19/07/2016
#=========================================================

#=========================================================
# adaptions made by Seb since last version
# adapted R#324 -> round digits increased to 8, otherwise gives error for small depths (e.g. O2 microsensor profiles)
# adapted R#196 -> the porosity extracted from the constant value now gives 1 for x-values lower than 0 (i.e. in the overlying water)
# adapted R#331 -> made R skip NA values to calculate max and min value
# adapted R#192 -> added line to remove NA concentration values from the input frame
# 22/03/2022 SVDV: -> moved sourcing of other functions to top of file (to avoid having to stay in right working dir)
#                  -> changed default window selection parm for Savitzky-Golay to "interactive"
#                  -> changed default i.end parameter so that it is at most 2 less than the number of unique 
#                     depths supplied
# 23/03/2022 SVDV: made sure the values analyzed with the discrete function are only positive depths (in the sediment)
#=========================================================

# =========================================================
# load required packages
# =========================================================

require(marelac)
#require(signal) -> is included in fitprofile package
#require(fractaldim) -> is included in fitprofile package
source("fittingfuncs/gradientfit v01.R")
source("fittingfuncs/SavitskyGolay v01.R")
source("fittingfuncs/profile v01.R")
source("fittingfuncs/user_input.R")

# =========================================================
# auxiliary function to generate default parameter sets
# private
# =========================================================

generate.default.parms <- function(input,species,set){
  if (set == "continuous"){
    p                      <- 3
    bnd.upper <- bnd.lower <- 1
    n.C <- n.J <- n.R      <- NULL
    n.uniform              <- FALSE
    optimal.window.size    <- "interactive" 
    min.n                  <- p - p%%2
    max.n                  <- nrow(input)%/%2-1
    keep.graphics          <- FALSE
    interpolation          <- "interpolate"
    continuous.parms.default <- list(p=p,bnd.upper=bnd.upper,bnd.lower=bnd.lower,n.C=n.C,n.J=n.J,n.R=n.R,n.uniform=n.uniform,
                             optimal.window.size=optimal.window.size,min.n=min.n,max.n=max.n,keep.graphics=keep.graphics,
                             interpolation=interpolation)
  return(continuous.parms.default)  
  }
 if (set == "gradient"){
   x.limits   <- NA  
   #E.SWI      <- 0
   gradient.parms.default <- list(x.limits = x.limits)#, E.SWI = E.SWI)
   
  return(gradient.parms.default)
 }
 if(set =="discrete"){
   L.down <- max(input[(input$x>=0.0),"x"]) 
   x.up   <- min(input[(input$x>=0.0),"x"])
   irr    <- 0 
   irr.att<- 0.03
   N      <- 200
   i.end  <- min(12,length(unique(input$x[!is.na(input$C)]))-2)
   initial.zones <- NULL
   prob      <- 0.01
   UBC       <- "conc.up"
   LBC       <- "no.flux"
   flux.up   <- NULL
   flux.down <- NULL
   discrete.parms.default <- list(L.down = L.down, x.up = x.up, irr = irr, irr.att = irr.att,
                                 prob = prob, UBC = UBC, LBC = LBC, flux.up = flux.up, flux.down = flux.down,
                                 N = N, i.end = i.end, initial.zones = initial.zones)
   
  return(discrete.parms.default) 
 }
 if (set == "environmental"){
   TC   <- 10    # temp [deg C]
   S    <- 30    # salinity [PSU]
   P    <- 1.013 # pressure [bar]
   z    <- 0     # charge of the ion at hand []
   #Dmol <- diffcoeff(S=S,t=TC,P=P,species=species)[[species]]*1E4*3600*24 # [cm2 d-1]
   env.parms.default <- list(TC = TC, S = S, P = P, z = z)#, Dmol = Dmol)
   
  return(env.parms.default) 
 }
}

# =========================================================
# auxiliary function to fill in user supplied parameter values in default parameter lists
# private
# =========================================================

fill.in.supplied.values <- function(parms.input,parms.default){
  parms.new <- parms.default
  parms.names <- names(parms.default)
  for (i in 1:length(parms.names)){
    if (!is.null(parms.input[[parms.names[[i]]]])) {parms.new[[parms.names[[i]]]] <- parms.input[[parms.names[[i]]]]}
  }
  return(parms.new)
}

# =========================================================
# MAIN FUNCTION
# =========================================================

FLIPPER.func <- function(input,por.cte=NA,E.cte=NULL,tort.dep=1,species,method=NULL,env.parms=NULL,full.output=FALSE,
                         # optional input for discrete
                         discrete.parms=NULL,
                         # optional input for continuous
                         continuous.parms=NULL,
                         # optional input for gradient
                         gradient.parms=NULL
                         ){

# =============================================================================
# Arguments 
# =============================================================================
  
# input: data for inputting: a dataframe containing 
#                           depth x         [m]
#                           concentration C [mmol m-3]
#                optionally porosity (if depth dependent) [dimensionless]
#                optionally electric field E (if depth dependent) [V m-1]
#        optionally for continuous: 
#               advective velociy (v),
#               first derivatives of porosity (dpor_dx),
#               advective velociy (dv_dx), effective diffusion coefficient (dDs_dx)
# por.cte: a constant value for porosity -> if not defined, porosity should be included 
#          in the input dataframe
# E.cte: a constant value for the electric field (optional),      [V m-1]
#     supplied as vector: c("value E","start depth","end depth")
#                                start depth = SWI (generally 0) [m]
#                                end depth   = deepest value     [m]
# tort.dep: the estimation method for tortuosity:1 (default) = 1-2ln(por), 
#                                                2           = por^-1,
#                                                3           = por^-2,
#                                                4           = 1+3(1-por)
# species: species of interest (input as in marelac; e.g. c("O2"), c("Fe") ...)  
# method: prefered method, can be "discrete", "continuous", "gradient"
#         if NULL or "all" -> all methods are run (and a suggestion of the fit-for-purpose method is given)
# full.output: logical, TRUE gives all possible output, FALSE gives a cleaned output (see below)
# env.parms: (optional) a list containing temperature TC         [degrees Celcius]
#                                     salinity S                 [PSU] 
#                                     pressure P                 [bar]
#                                     diffusion coefficient Dmol [m2 d-1]
#                                     the charge of the ion z    [dimensionless]
#
# parameters for discrete (optional), supplied as list (discrete.parms, only the the parms that have to be changed from default value have to be supplied) 
#        i.end: maximum number of zones to be used in discrete
#        initial.zones: start number of zones for lumping 
#
# parameters for continuous (optional), supplied as list (continuous.parms, only the the parms that have to be changed from default value have to be supplied) 
#        p: order of the polynomial that is fitted (typically 2 or 3) 
#        bnd.upper: integer determining filter behaviour at the upstream (upper, left) boundary. Value = 1: no constraint on flux (default value). Value = 2: constant flux imposed. Value = 3: zero flux imposed.   
#        bnd.lower: integer determining filter behaviour at the downstream (lower or right) boundary. Value = 1: no constraint on flux (default value). Value = 2: constant flux imposed. Value = 3: zero flux imposed.   
#        n.C : window size used in filtering the concentration C, representing the number of data points to the
#           left and right of the data midpoint (hence the filter length  = 2*n.C+1)
#        n.J : window size used in filtering the first order derivative dC_dx, representing the number of data points to the
#           left and right of the data midpoint (hence the filter length  = 2*n.J+1)
#        n.R : window size used in filtering the second order derivative d2C_dx2, representing the number of data points to the
#           left and right of the data midpoint (hence the filter length  = 2*n.R+1)
#        n.uniform : logical, if TRUE then n.C, n.J, and n.R are all set uniform to max(n.C,n.J,n.R)
#        min.n : minimum value of the window size, used when scanning the optimal window size 
#        max.n : maximum value of the window size, used when scanning the optimal window size
#        optimal.window.size : either "automated" or "interactive" (default) 
#                       interactive allows user to select ideal filter size
#                       automated lets function select ideal filter size
#        keep.grapics: either TRUE or FALSE (default) - keeps windows created by automated function
#        full.output : logical, TRUE gives all output, FALSE gives selected output (without info, see fitprofile_package)
#        interpolation: which mode of interpolation of the data is used (if necessary): average (take maximum stepsize) or interpolate (take minimum stepsize)
#
# parameters for gradient fit function (optional) supplied as list (gradient.parms, only the the parms that have to be changed from default value have to be supplied) 
#        x.limits: depth limits of the profile (a vector of 2 -> the upper depth and the lower depth)
#=============================================================================
# OUTPUT
#=============================================================================
# a list with (if applicable):
# input:             a list with user.input: the user supplied dataframe (x, C, por, tort ...)
#                                continuous.input:  the interpolated dataframe for use in the continuous function
# env.parms:         a list with the inputted environmental parameters (default with user supplied)
# gradient.method:   the output of the gradient function
# gradient.parms:    the parameters used in the gradient function (default with user supplied)
# continuous.method:         the output of the continuous function
# continuous.parms:          the parameters used in the continuous function
# discrete.method:    the output of the discrete function
# discrete.parms:     the parameters supplied in the discrete function
#=============================================================================
# ---------------------
# error 
# ---------------------

if (is.null(method)){
  stop("error: method should be one of the following possibilities: \n \"gradient\", \"discrete\", \"continuous\", \"all\"")
}

if (method != "gradient" & method != "discrete" & method != "continuous" & method != "all") {
  stop("error: method should be one of the following possibilities: \n \"gradient\", \"discrete\", \"continuous\", \"all\"")
}

# ---------------------
# prepare input dataframe for general use
# ---------------------

# remove rows with NA concentrations

if (any(is.na(input$C))){
input <- input[-c(which(is.na(input$C))),]  
}

# remove values in the watercolumn
  
  #if (any(input$x<0.0)){
  #  input <- input[(input$x>=0.0),]  
#}
  
# check input of porosity + create por column in input dataframe

if (!is.null(input$por)){
  if(!is.na(por.cte)){
  print("Warning: you have supplied porosity in the input dataframe AND as a constant. FLIPPER will only take into account the porosity supplied in the input dataframe")
}}
  
if (is.null(input$por)) {
  if (is.na(por.cte)){
  stop("porosity should be defined as a constant (por.cte) or included in the input dataframe (input$por)")
}
  if (!is.na(por.cte)){
  input$por <- rep(por.cte,length(input$x))
  input$por[which(input$x < 0)] <- 1
}}

# create tor column in input dataframe

if (tort.dep == 1){ input$tort <- 1-2*log(as.numeric(input$por))   } 
if (tort.dep == 2){ input$tort <- (input$por)^(-1)     }
if (tort.dep == 3){ input$tort <- (input$por)^(-2)     }
if (tort.dep == 4){ input$tort <- 1+3*(1-(input$por))  }

# create E column in input dataframe

if (!is.null(input$E)){
  if(!is.null(E.cte)){
  print("Warning: you have supplied the electrical field in the input dataframe AND as a constant. FLIPPER will only take into account the electrical field supplied in the input dataframe")
}}
  
if (is.null(input$E)) {
  if (!is.null(E.cte)){
    start <- min(which(input$x >= E.cte[2]))
    end   <- max(which(input$x <= E.cte[3]))
    input$E <- rep(0,length(input$x))
    input$E[start:end] <- E.cte[1]
  }
}
  
# ---------------------
# Create environmental parameter list (defaults with user supplied values)
# ---------------------

env.parms.default <- generate.default.parms(input,species,set="environmental")

if (!is.null(env.parms)){
  env.parms <- fill.in.supplied.values(env.parms,env.parms.default); rm(env.parms.default)
} else {env.parms <- env.parms.default; rm(env.parms.default)}

if (is.null(env.parms$Dmol)){
env.parms$Dmol <- diffcoeff(S=env.parms$S,t=env.parms$TC,P=env.parms$P,species=species)[[species]]*3600*24 # [m2 d-1]
}

# ---------------------
# Create gradient parameter list (defaults with user supplied values)
# ---------------------

if (method == "gradient" | method == "all"){
    
gradient.parms.default <- generate.default.parms(input,species,set="gradient")

if (!is.null(gradient.parms)){
  gradient.parms <- fill.in.supplied.values(gradient.parms,gradient.parms.default); rm(gradient.parms.default)
} else {gradient.parms <- gradient.parms.default; rm(gradient.parms.default)}
}

# ---------------------
# Create PROFILE parameter list (defaults with user supplied values)
# ---------------------

if (method == "discrete" | method == "all"){
  
discrete.parms.default <- generate.default.parms(input,species,set="discrete")
  
  if (!is.null(discrete.parms)){
    discrete.parms <- fill.in.supplied.values(discrete.parms,discrete.parms.default); rm(discrete.parms.default)
  } else {discrete.parms <- discrete.parms.default; rm(discrete.parms.default)}

if (discrete.parms$LBC == "conc.down" & is.null(discrete.parms$C.down)){
discrete.parms$C.down     <- mean(input$C[which(input$x == discrete.parms$L.down)])
}

discrete.parms$C.up     <- mean(input$C[which(input$x == discrete.parms$x.up)])
discrete.parms$grid     <- setup.grid.1D(x.up = discrete.parms$x.up, x.down = discrete.parms$L.down, N = discrete.parms$N)
discrete.parms$por.grid <- setup.prop.1D(xy=cbind(input$x,input$por),interpolate="linear", grid = discrete.parms$grid)
discrete.parms$irr.grid <- setup.prop.1D(func=p.exp,grid=discrete.parms$grid,y.0=discrete.parms$irr,y.inf=0,x.att=discrete.parms$irr.att) 
discrete.parms$Ds       <- setup.prop.1D(xy=cbind(input$x,(env.parms$Dmol/input$tort)),interpolate="linear", grid = discrete.parms$grid)


if (is.null(input$v)) {
  if (!is.null(input$E)){
    Faraday <- as.numeric(Constants$F[1])      # Faraday constant [C mol-1]
    R       <- as.numeric(Constants$gasCt2[1]) # Molar gas constant [J mol-1 K-1]
    v.E <- env.parms$Dmol/input$tort*env.parms$z*Faraday*input$E/(R*(env.parms$TC+273.15)) # m d-1
    input$v <- v.E }                                       # advective velocity due to electrical field [m d-1]
  else{v <- 0; input$v <- rep.int(v,length(input$x))}}   # advective velocity [m d-1]

discrete.parms$v <- setup.prop.1D(xy=cbind(input$x,input$v),interpolate="linear", grid = discrete.parms$grid)
}

# ---------------------
# Prepare data.frame for continuous list (defaults with user supplied values)
# ---------------------

if (method == "continuous" | method == "all"){

# Generate parameter list
  
continuous.parms.default <- generate.default.parms(input,species,set="continuous")
  
if (!is.null(continuous.parms)){
continuous.parms <- fill.in.supplied.values(continuous.parms,continuous.parms.default); rm(continuous.parms.default)
} else {continuous.parms <- continuous.parms.default; rm(continuous.parms.default)}
  
# prepare input dataframe for specific use 
# generate extra necessary columns

if (is.null(input$Ds)) {input$Ds  <- env.parms$Dmol/input$tort} # effective diffusion coefficient [L2 T-1]

if (is.null(input$v)) {
  if (!is.null(input$E)){
  Faraday <- as.numeric(Constants$F[1])      # Faraday constant [C mol-1]
  R       <- as.numeric(Constants$gasCt2[1]) # Molar gas constant [J mol-1 K-1]
  v.E <- env.parms$Dmol/input$tort*env.parms$z*Faraday*input$E/(R*(env.parms$TC+273.15)) # m d-1
  input$v <- v.E }                                       # advective velocity due to electrical field [m d-1]
  else{v <- 0; input$v <- rep.int(v,length(input$x))}}   # advective velocity [m d-1]
if (is.null(input$dpor_dx)){dpor_dx <- 0; input$dpor_dx <- rep.int(dpor_dx,length(input$x))} # derivative of porosity [m-1]              
if (is.null(input$dv_dx))  {dv_dx <- 0; input$dv_dx   <- rep.int(dv_dx,length(input$x))}     # derivative of advective velocity [d-1]
if (is.null(input$dDs_dx)) {dDs_dx <- 0; input$dDs_dx  <- rep.int(dDs_dx,length(input$x))}   # derivative of effective diffusion coefficient [m d-1]

# All points must be equidistant for continuous
# Make all points equidistant by linear interpolation
# interpolate = estimate values between data points
# average     = average out values

intervals <- unique(round(diff(input$x),digits=8))
if (length(intervals) > 1){
h.max <- max(intervals)
h.min <- min(intervals)

if (continuous.parms$interpolation == "interpolate"){
  new.x <- seq(min(input$x,na.rm=T),max(input$x,na.rm=T),h.min)
  rm(h.max)}
if (continuous.parms$interpolation == "average"){
  new.x <- seq(min(input$x,na.rm=T),max(input$x,na.rm=T),h.max)
  rm(h.min)}

input.continuous <- data.frame(
  x    = new.x,
  C    = approx(input$x,input$C,new.x,rule=2)$y,
  por  = approx(input$x,input$por,new.x,rule=2)$y,
  tort = approx(input$x,input$tort,new.x,rule=2)$y,
  Ds   = approx(input$x,input$Ds,new.x,rule=2)$y,
  v    = approx(input$x,input$v,new.x,rule=2)$y,
  dpor_dx = approx(input$x,input$dpor_dx,new.x,rule=2)$y,
  dv_dx   = approx(input$x,input$dv_dx,new.x,rule=2)$y,
  dDs_dx  = approx(input$x,input$dDs_dx,new.x,rule=2)$y
)} else {input.continuous <- input}

if (continuous.parms$max.n == nrow(input)%/%2-1){
  continuous.parms$max.n <- nrow(input.continuous)%/%2-1
}

}

# ---------------------
# run requested method
# ---------------------

  if (method == "discrete"){
    output.pr <- fit.profile(input=input,parms=discrete.parms,i.end=discrete.parms$i.end,initial.zones=discrete.parms$initial.zones, full.output = full.output)
    return(list("method" = method,"input" = input, "output" = output.pr, "parms" = list("env.parms" = env.parms, "discrete.parms" = discrete.parms)))
    
  }
  if (method == "continuous"){
    
    # Run continuous method  
    output.continuous <- SavGolay.analysis(profile = input.continuous, p=continuous.parms$p, n.C=continuous.parms$n.C,
                                   n.J=continuous.parms$n.J,n.R=continuous.parms$n.R, bnd.upper = continuous.parms$bnd.upper, 
                                   bnd.lower = continuous.parms$bnd.lower, n.uniform = continuous.parms$n.uniform,
                                   optimal.window.size = continuous.parms$optimal.window.size, 
                                   min.n = continuous.parms$min.n, max.n = continuous.parms$max.n, 
                                   keep.graphics = continuous.parms$keep.graphics, full.output = full.output)
    return(list("method" = method, 
                "input" = list("user.input" = input,"continuous.input"= input.continuous), 
                "output" = output.continuous,
                "parms" = list("env.parms" = env.parms, "continuous.parms" = continuous.parms) 
                 
                ))
  }

  if (method == "gradient"){
    # Run and save gradient method  
    output.gf <- gradient.fit(profile = input, x.limits = gradient.parms$x.limits, env.parms = env.parms, full.output = full.output)
    return(list("method" = method, "input" = input, "output" = output.gf, "parms" = list("env.parms" = env.parms, "gradient.parms" = gradient.parms)))
  }

  if (method == "all"){
    # Run and save gradient method  
    print("Run gradient method")
    output.gf <- try(gradient.fit(profile=input,x.limits=gradient.parms$x.limits,env.parms=env.parms, full.output = full.output))
    
    print("Run discrete method")
    # Run and save discrete method  
    output.ds <- try(fit.profile(input=input,parms=discrete.parms,i.end=discrete.parms$i.end,initial.zones=discrete.parms$initial.zones, full.output = full.output))
    
    print("Run continuous method")
    # Run and save continuous method  
    output.cs <- try(SavGolay.analysis(profile = input.continuous, p=continuous.parms$p, n.C=continuous.parms$n.C,
                                   n.J=continuous.parms$n.J,n.R=continuous.parms$n.R, bnd.upper = continuous.parms$bnd.upper, 
                                   bnd.lower = continuous.parms$bnd.lower, n.uniform = continuous.parms$n.uniform,
                                   optimal.window.size = continuous.parms$optimal.window.size, 
                                   min.n = continuous.parms$min.n, max.n = continuous.parms$max.n, 
                                   keep.graphics = continuous.parms$keep.graphics, full.output = full.output))
    
    if (typeof(output.gf) != "list"){output.gf <- c("FAILED")}
    if (typeof(output.ds) != "list"){output.ds <- c("FAILED")}
    if (typeof(output.cs) != "list"){output.cs <- c("FAILED")}
    return(list("method" = method,"input" = list("user.input" = input,"continuous.input"= input.continuous), 
                "output"= list("gradient.output" = output.gf, "discrete" = output.ds,"continuous" = output.cs), 
                "parms"= list("env.parms" = env.parms, "gradient.parms" = gradient.parms, "discrete.parms" = discrete.parms, "continuous.parms" = continuous.parms)))
    
    }
}
