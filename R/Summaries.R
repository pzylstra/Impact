#' Summary table of surface results
#'
#' Summarises FRaME generated surface fire behaviour by RepId
#'
#' @param surface The dataframe res$SurfaceResults
#' @return dataframe

surf <- function(surface)
{
  out <- surface %>%
    group_by(repId) %>%
    mutate(lengthSurface = flameLength,
           heightSurface = flameHeight,
           angleSurface = flameAngle)%>%
    select(repId, lengthSurface, heightSurface, angleSurface) %>%
    summarize_all(max)
  
  return(out)
}

#####################################################################
#' Summary table of stratum results
#'
#' Summarises FRaME generated fire behaviour by stratum and RepId
#'
#' @param flames The dataframe res$FlameSummaries
#' @param sites The dataframe res$Sites
#' @param ros The dataframe res$ROS
#' @param surf The dataframe produced by the extension surf
#' @return stratum dataframe

stratum <- function(flames, sites, ros, surf)
{
  y <- ros%>%
    select(repId, level, ros)
  z <- flames %>%
    select(repId, level, flameLength, flameAngle, flameHeight)
  
  a <- y %>%
    left_join(z) %>%
    left_join(sites) %>%
    
    # Strata without ros will end up with NA values
    # after doing the join above. Convert these missing values to zero.
    mutate(ros = ifelse(is.na(ros), 0.0, ros),
           flameHeight = ifelse(is.na(flameHeight), 0.0, flameHeight),
           flameLength = ifelse(is.na(flameLength), 0.0, flameLength),
           flameAngle = ifelse(is.na(flameAngle), 0.0, flameAngle),
           #PATCH TO COVER MOISTURE EXTINCTION UNTIL FIXED IN SCALA
           #Duplicate DFMC, then create binary for spread/no spread
           extinct = deadFuelMoistureProp,
           extinct = ifelse(extinct == 0.199, 0.0, 1.0)) %>%
    select(repId, level, fuelLoad, flameHeight, flameLength, flameAngle, ros, windSpeed,
           deadFuelMoistureProp, temperature, slope, extinct) %>%
    mutate(oHorizon = fuelLoad * 10,
           slope_degrees = slope * 180 / pi,
           flameA_degrees = flameAngle * 180 / pi,
           ros_kph = extinct * ros * 3.6,
           heightPlant = flameHeight * extinct,
           lengthPlant = flameLength * extinct,
           wind_kph = windSpeed * 3.6,
           spread = ifelse(ros > 0, 1, 0),
           has.flame = spread + (extinct * flameHeight) > 0)
  
  # Add in surface flame descriptors
  rep <- max(a$repId)
  st <- as.numeric(count(a))/rep
  i <- 1
  for(loop in 1:rep) {
    a$flameHeight[i]=surf$heightSurface[loop]
    a$flameLength[i]=surf$lengthSurface[loop]
    a$flameAngle[i]=surf$angleSurface[loop]
    i <- i + st
  }
  return(a)
}

#####################################################################
#' Summary table of fire behaviour
#'
#' Summarises FRaME generated fire behaviour by RepId
#'
#' @param surf The dataframe surf
#' @param stratum The dataframe stratum
#' @return summary dataframe

summary <- function(stratum, surf)
{
  return(stratum %>%
           select(repId, slope_degrees, wind_kph, deadFuelMoistureProp, temperature,
                  heightPlant, lengthPlant, flameAngle, ros_kph, extinct) %>%
           group_by(repId) %>%
           summarize_all(max) %>%
           left_join(surf) %>%
           mutate(fh = pmax(heightSurface, heightPlant) * extinct,
                  fl = pmax(lengthSurface, lengthPlant) * extinct,
                  zeta = 2.5*ros_kph,
                  epsilon = 1-exp(-0.72*zeta)))
}

#####################################################################
#' Representative flame descriptors
#'
#' Summarises FRaME generated flame segments into a combined,
#' representative plant flame for each repId where plants ignited
#'
#' @param IP The dataframe res$IgnitionPaths
#' @return dataframe

repFlame <- function(IP)
{
  top <- IP %>%
    mutate(angle = atan((y1 - y0)/(x1 - x0)),
           repHeight = flameLength*sin(angle)+y0)%>%
    group_by(repId) %>%
    summarize_all(max) %>%
    select(repId, repHeight)
  
  repFlame <- IP %>%
    mutate(repAngle = atan((y1 - y0)/(x1 - x0))
    ) %>%
    select(repId, repAngle)%>%
    group_by(repId) %>%
    summarize_all(mean) %>%
    left_join(top) %>%
    mutate(repLength = repHeight/sin(repAngle)) %>%
    select(repId, repHeight, repLength, repAngle)%>%
    right_join(IP)
  
  return(repFlame)
}


#####################################################################
#' Stratum descriptors from a param file
#'
#' For each stratum, finds mean crown width, plant separation, and number of species
#'
#' @param Param A parameter dataframe used for FRaME,
#' such as produced using readLegacyParamFile
#' @return dataframe of stratum descriptors

strata <- function(Param)
{
  #Number of strata
  StL <- count(Param)-13
  StN <- Param$stratum[max(StL$n)]
  
  #Count species per stratum
  Sp <- numeric(StN)
  for(sn in 1:StN){
    strat <- filter(Param, stratum == sn)
    strat <- na.omit(strat)
    Sp[sn] <- (as.numeric(max(strat$species))+1)-as.numeric(min(strat$species))
  }
  
  #COLLECT DIMENSIONS
  width <- Param[Param$param == "w", ]
  comp <- Param[Param$param == "composition", ]
  sep <- Param[Param$param == "plantSeparation", ]
  peak <- Param[Param$param == "hp", ]
  top <- Param[Param$param == "ht", ]
  edge <- Param[Param$param == "he", ]
  base <- Param[Param$param == "hc", ]
  level <- Param[Param$param == "levelName", ]
  name <- Param[Param$param == "levelName", ]
  
  #BUILD TABLE
  n <- as.data.frame(list('stratum'=name$stratum, 'name'=name$value, 'speciesN'=Sp))
  s <- as.data.frame(list('stratum'=width$stratum, 'comp'=comp$value, 'width'=width$value, 'Hp'=peak$value,
                          'Ht'=top$value, 'He'=edge$value, 'Hc'=base$value))%>%
    mutate(Co = as.numeric(as.character(comp)),
           Ww = as.numeric(as.character(width))*as.numeric(as.character(comp)),
           Wp = as.numeric(as.character(Hp))*as.numeric(as.character(comp)),
           Wt = as.numeric(as.character(Ht))*as.numeric(as.character(comp)),
           We = as.numeric(as.character(He))*as.numeric(as.character(comp)),
           Wc = as.numeric(as.character(Hc))*as.numeric(as.character(comp)),
           top = pmax(Wp,Wt),
           base = pmin(We,Wc))%>%
    group_by(stratum) %>%
    summarize_if(is.numeric,sum)%>%
    mutate(width = Ww/Co,
           top = top/Co,
           base = base/Co)%>%
    left_join(sep, by = "stratum")%>%
    mutate(separation = as.numeric(value),
           cover = width^2/separation^2)%>%
    select(stratum, separation, cover, width, base, top)
  
  strata <- as.data.frame(s)%>%
    left_join(n, by="stratum")
  return(strata)
}
#####################################################################
#' Species descriptors from a param file
#'
#' Finds dimensions and moisture of each species
#'
#' @param Param A parameter dataframe used for FRaME,
#' such as produced using readLegacyParamFile
#' @return dataframe of species descriptors

species <- function(Param)
{
  #Collect traits
  sp <- Param[Param$param == "name", ]
  lfmc <- Param[Param$param == "liveLeafMoisture", ]
  Peak <- Param[Param$param == "hp", ]
  Top <- Param[Param$param == "ht", ]
  Edge <- Param[Param$param == "he", ]
  Base <- Param[Param$param == "hc", ]
  Width <- Param[Param$param == "w", ]
  
  species <- as.data.frame(list('name'=sp$value, 'hp'=as.numeric(Peak$value),'ht'=as.numeric(Top$value),
                                'hc'=as.numeric(Base$value), 'he'=as.numeric(Edge$value),
                                'w'=as.numeric(Width$value), 'lfmc'=as.numeric(lfmc$value))) %>%
    mutate(htR = ht/hp,
           hcR = hc/hp,
           heR = he/hp,
           wR = w/hp)
  return(species)
}

