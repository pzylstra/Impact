#' Finds the LAI for a plant above a height threshold
#'
#' @param base.params A parameter file
#' @param Sp The number of the species in the parameter file
#' @param y The height above which the LAI will be measured (m)
#' @return Numeric value
#' @export


LAIp <- function(base.params, sp = 1, y = 0)
{
  # Subset species
  spPar <- subset(base.params, species == sp)
  
  # Collect parameters
  Cd <- as.numeric(spPar$value[spPar$param == "clumpDiameter"])
  Cs <- as.numeric(spPar$value[spPar$param == "clumpSeparation"])
  ram <- as.numeric(spPar$value[spPar$param == "stemOrder"])
  comp <- as.numeric(spPar$value[spPar$param == "composition"])
  hc <- max(y, as.numeric(spPar$value[spPar$param == "hc"]))
  he <- max(y, as.numeric(spPar$value[spPar$param == "he"]))
  ht <- max(y, as.numeric(spPar$value[spPar$param == "ht"]))
  hp <- max(y, as.numeric(spPar$value[spPar$param == "hp"]))
  W <- as.numeric(spPar$value[spPar$param == "w"])
  ll <- as.numeric(spPar$value[spPar$param == "leafLength"])
  lw <- as.numeric(spPar$value[spPar$param == "leafWidth"])
  ls <- as.numeric(spPar$value[spPar$param == "leafSeparation"])
  
  # Find w at cutoff cones
  topRise <- abs(as.numeric(spPar$value[spPar$param == "hp"])-as.numeric(spPar$value[spPar$param == "ht"]))
  w <- ifelse(topRise==0,
              W,
              (W/topRise)*abs(hp-ht))
  
  baseFall <- abs(as.numeric(spPar$value[spPar$param == "he"])-as.numeric(spPar$value[spPar$param == "hc"]))
  wl <- ifelse(baseFall==0,
               W,
               (W/baseFall)*abs(he-hc))
  
  # Leaf area per clump
  clumpLeaves <- 0.88*(Cd*ram/ls)^1.18
  laClump <- (ll*lw/2) * clumpLeaves
  
  # Plant volume above height y
  upper <- ((hp-ht)*(pi*(w/2)^2))/3
  centre <- (ht-he)*(pi*(w/2)^2)
  lower <- ((he-hc)*(pi*(wl/2)^2))/3 
  vol <- upper+centre+lower
  
  # Clumps above height y
  volC <- ((4/3)*pi*((Cd+Cs)/2)^3)
  nClumps <- vol/volC
  
  #LAI per plant above height y
  LA <- laClump * nClumps
  return(LA/(pi*(w/2)^2))
}

#####################################################################