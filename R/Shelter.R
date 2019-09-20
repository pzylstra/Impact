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
  #  hc <- min(yu,max(yl, as.numeric(spPar$value[spPar$param == "hc"])))
  hc <- as.numeric(spPar$value[spPar$param == "hc"])
  #  he <- min(yu,max(yl, as.numeric(spPar$value[spPar$param == "he"])))
  he <- as.numeric(spPar$value[spPar$param == "he"])
  #  ht <- min(yu,max(yl, as.numeric(spPar$value[spPar$param == "ht"])))
  ht <- as.numeric(spPar$value[spPar$param == "ht"])
  #  hp <- min(yu,max(yl, as.numeric(spPar$value[spPar$param == "hp"])))
  hp <- as.numeric(spPar$value[spPar$param == "hp"])
  W <- as.numeric(spPar$value[spPar$param == "w"])
  ll <- as.numeric(spPar$value[spPar$param == "leafLength"])
  lw <- as.numeric(spPar$value[spPar$param == "leafWidth"])
  ls <- as.numeric(spPar$value[spPar$param == "leafSeparation"])
  
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
  upperAboveSlice <- max(0,(hp-yu)*(pi*(wa/2)^2))/3
  upperFull <- ((hp-max(yl,ht))*(pi*(wb/2)^2))/3
  upper <- upperFull-upperAboveSlice
  centre <- (min(yu,ht)-max(yl,he))*(pi*(W/2)^2)
  lowerFull <- ((min(yu,he)-hc)*(pi*(wc/2)^2))/3
  lowerBelowSlice <- (max(0,(yl-hc))*(pi*(wd/2)^2))/3 
  lower <- lowerFull-lowerBelowSlice
  vol <- upper+centre+lower
  
  # Clumps in slice
  volC <- ((4/3)*pi*((Cd+Cs)/2)^3)
  nClumps <- vol/volC
  
  #LAI per plant in slice
  LA <- laClump * nClumps
  l<- ifelse(max(wa,wb,wc,wd)==0,
             0,
             LA/(pi*(max(wa,wb,wc,wd)/2)^2))
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
  return(LAI)
}

#####################################################################