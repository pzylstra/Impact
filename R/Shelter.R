#' Finds the LAI for a horizontal slice of a plant
#'
#' @param base.params A parameter file
#' @param Sp The number of the species in the parameter file
#' @param yu The upper height of the horizontal slice (m)
#' @param yl The lower height of the horizontal slice (m)
#' @return Numeric value
#' @export

LAIp <- function(base.params, sp = 1, yu = 100, yl = 0)
{
  
  # Subset species
  spPar <- subset(base.params, species == sp)
  
  # Collect parameters
  Cd <- as.numeric(spPar$value[spPar$param == "clumpDiameter"])
  Cs <- as.numeric(spPar$value[spPar$param == "clumpSeparation"])
  ram <- as.numeric(spPar$value[spPar$param == "stemOrder"])
  hc <- as.numeric(spPar$value[spPar$param == "hc"])
  he <- as.numeric(spPar$value[spPar$param == "he"])
  ht <- as.numeric(spPar$value[spPar$param == "ht"])
  hp <- as.numeric(spPar$value[spPar$param == "hp"])
  W <- as.numeric(spPar$value[spPar$param == "w"])
  ll <- as.numeric(spPar$value[spPar$param == "leafLength"])
  lw <- as.numeric(spPar$value[spPar$param == "leafWidth"])
  ls <- as.numeric(spPar$value[spPar$param == "leafSeparation"])
  
  if (yu <= min(hc,he)) {
    l <- 0
  } else if (yl >= max(ht,hp)) {
    l <- 0
  } else {
    
    # Find w at cutoff cones
    topRise <- abs(hp-ht)
    baseFall <- abs(he-hc)
    # Top of slice, plant top
    wa <- ifelse(topRise==0,
                 W,
                 (W/topRise)*max(0,(hp-max(ht,yu))))
    # Base of slice, plant top
    wb <- ifelse(topRise==0,
                 W,
                 (W/topRise)*abs(hp-max(yl,ht)))
    # Top of slice, plant base
    wc <- ifelse(baseFall==0,
                 W,
                 (W/baseFall)*abs(min(he,yu)-hc))
    # Base of slice, plant base
    wd <- ifelse(baseFall==0,
                 W,
                 (W/baseFall)*max(0,(min(he,yl)-hc)))
    
    # Leaf area per clump
    clumpLeaves <- 0.88*(Cd*ram/ls)^1.18
    laClump <- (ll*lw/2) * clumpLeaves
    
    # Plant volume in slice
    upperAboveSlice <- max(0,(hp-max(ht,yu))*(pi*(wa/2)^2))/3
    upperFull <- ((hp-max(yl,ht))*(pi*(wb/2)^2))/3
    upper <- max(0,upperFull-upperAboveSlice)
    centre <- max(0,(min(yu,ht)-max(yl,he)))*(pi*(W/2)^2)
    lowerFull <- ((min(yu,he)-hc)*(pi*(wc/2)^2))/3
    lowerBelowSlice <- (max(0,(min(he,yl)-hc))*(pi*(wd/2)^2))/3
    lower <- max(0,lowerFull-lowerBelowSlice)
    vol <- upper+centre+lower
    
    # Clumps in slice
    volC <- ((4/3)*pi*((Cd+Cs)/2)^3)
    nClumps <- vol/volC
    
    #LAI per plant in slice
    LA <- laClump * nClumps
    l<- ifelse(max(wa,wb,wc,wd)==0,
               0,
               LA/(pi*(max(wa,wb,wc,wd)/2)^2))
  }
  return(l)
}


#####################################################################

#' Finds the LAI for a horizontal slice of a community
#'
#' @param base.params A parameter file
#' @param yu The upper height of the horizontal slice (m)
#' @param yl The lower height of the horizontal slice (m)
#' @return Numeric value
#' @export

LAI <- function(base.params, yu = 100, yl = 0)
{
  # Collect plant figures
  l <- data.frame()
  c <- data.frame()
  s <- data.frame()
  w <- data.frame()
  N <- count(species(base.params))
  str <- strata(base.params)
  n <- 1
  while(n <= N[1]) {
    spPar <- subset(base.params, species == n)
    laiN <- LAIp(base.params, sp = n, yu = yu, yl = yl)
    l <- rbind(l,laiN)
    c <- rbind(c,as.numeric(spPar$value[spPar$param == "composition"]))
    s <- rbind(s, as.numeric(str$separation[as.numeric(spPar$stratum[1])]))
    w <- rbind(w,as.numeric(spPar$value[spPar$param == "w"]))
    n <- n + 1
  }
  
  # Construct table
  colnames(l) <- c("LAIp")
  l$ID <- seq.int(nrow(l))
  colnames(c) <- c("Weight")
  c$ID <- seq.int(nrow(c))
  colnames(s) <- c("Separation")
  s$ID <- seq.int(nrow(s))
  colnames(w) <- c("Width")
  w$ID <- seq.int(nrow(w))
  LAIplant <- left_join(l,c)%>%
    left_join(s)%>%
    mutate(Weight = ifelse(LAIp>0,
                           Weight,
                           0))%>%
    left_join(w)
  
  # Calculate LAI
  all <- sum(LAIplant$Weight)
  LAIplant <- LAIplant %>% 
    mutate(Weight = Weight/all,
           Cover = (Width^2/Separation^2)*Weight,
           LAIw = LAIp*Cover)
  LAI <- sum(LAIplant$LAIw)
  LAI[is.nan(LAI)] <- 0  
  return(LAI)
}

#####################################################################

#' Finds the deterministic vertical wind profile for a community
#'
#' @param base.params A parameter file
#' @param slices The number of vertical slices to measure
#' @return table
#' @export

profileDet <- function(base.params, slices = 10)
{

# Collect slice details
top <- max(species(base.params)$hp)
slice <- top/slices
yu <- top

# Loop through slices
l <- data.frame("l"=0)
gam <- data.frame("gam"=0)
W <- data.frame("w"=1)
w <- 1
n <- 1

while(n <= slices) {
  yl <- yu-slice
  LAIslice <- LAI(base.params = base.params, yl=yl, yu = yu)
  g <- 1.785*LAIslice^0.372
  w <- w*exp(g*((yl/yu)-1))
  l <- rbind(l,LAIslice)
  gam <- rbind(gam,g)
  W <- rbind(W,w)
  yu <- yl
  n = n+1
}

# Construct table
l$Slice <- seq.int(nrow(l))
gam$Slice <- seq.int(nrow(gam))
W$Slice <- seq.int(nrow(W))
wind <- left_join(l,gam)%>%
  left_join(W)%>%
  mutate(z = 1-((Slice-1)*(1/slices)),
         hm = z*top)
return(wind)
}

#####################################################################

#' Finds the deterministic vertical wind profile for a community, above a defined screen level
#'
#' @param base.params A parameter file
#' @param slices The number of vertical slices to measure above the measurement point
#' @param height The height (m) of the lower wind measurement point
#' @return table
#' @export

shelter <- function(base.params, slices = 10, SL = 1.2)
{
  
  # Collect slice details
  top <- max(species(base.params)$hp)
  slice <- (top-SL)/slices
  yu <- top
  
  # Loop through slices
  l <- data.frame("l"=0)
  gam <- data.frame("gam"=0)
  W <- data.frame("w"=1)
  w <- 1
  n <- 1
  
  while(n <= slices) {
    yl <- yu-slice
    LAIslice <- LAI(base.params = base.params, yl=yl, yu = yu)
    g <- 1.785*LAIslice^0.372
    w <- w*exp(g*((yl/yu)-1))
    l <- rbind(l,LAIslice)
    gam <- rbind(gam,g)
    W <- rbind(W,w)
    yu <- yl
    n = n+1
  }
  
  # Construct table
  l$Slice <- seq.int(nrow(l))
  gam$Slice <- seq.int(nrow(gam))
  W$Slice <- seq.int(nrow(W))
  wind <- left_join(l,gam)%>%
    left_join(W)%>%
    mutate(z = 1-((Slice-1)*(1/slices)),
           hm = z*top)
  return(wind)
}

##########################################################################

#' Calculates a non-deterministic wind profile
#'
#' @param base.params Parameter input table
#' @param reps Number of repetitions
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Variation A database of plant variability in traits, with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' Hs - Standard deviation of plant height variations
#' Hr - Truncates plant height variability by +/- Hr * height
#' @return dataframe
#' @export

windProfile <- function(base.params, Variation, reps = 10, slices = 10, l = 0.1,
                        Ms = 0.01, Pm = 1, Mr = 1.001)
{
  Strata <- strata(base.params)
  Species <- species(base.params)
  profW <- data.frame()
  
  if (reps > 0) {
    for (j in 1:reps) {
      base.params <- plantVarS(base.params, Strata, Species, Variation, l = l,
                               Ms = Ms, Pm = Pm, Mr = Mr)
      prof <- profileDet(base.params, slices = slices)%>%
        mutate(run = j)
      profW <- rbind(profW, prof)
    }
  }
  profW  <- profW[complete.cases(profW), ]
  return(profW)
}

#####################################################################

#' Calculates a non-deterministic wind profile on parallel cores
#' Useful when running a large number of reps
#'
#' @param base.params Parameter input table
#' @param reps Number of repetitions
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Variation A database of plant variability in traits, with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' Hs - Standard deviation of plant height variations
#' Hr - Truncates plant height variability by +/- Hr * height
#' @return dataframe
#' @export

windProfileP <- function(base.params, Variation, reps = 1000, slices = 100, l = 0.1,
                         Ms = 0.01, Pm = 1, Mr = 1.001)
{
  Strata <- impact::strata(base.params)
  Species <- impact::species(base.params)
  profW <- data.frame()
  f <- function(i) {
    base.params <- plantVarS(base.params, Strata, Species, Variation, l = l,
                             Ms = Ms, Pm = Pm, Mr = Mr)
    prof <- profileDet(base.params, slices = slices)
    profW <- rbind(profW, prof)
    return(profW)
  }
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl,
               {library(frame)
                 library(dplyr)
                 library(extraDistr)
                 library(impact)})
  # Setting "envir=environment()" allows clusterExport to access the 3 data frames created at the start of the function
  clusterExport(cl, c("slices", "l", "Ms", "Pm", "Mr", "base.params", "Variation", "Strata", "Species", "profW"), envir=environment())
  results <- parLapply(cl, 1:reps, f)
  stopCluster(cl)
  return(results)
}

#####################################################################
