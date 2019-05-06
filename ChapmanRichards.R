#' Models plant height from time since fire
#'
#' Uses either standard Chapman-Richards negative exponential or linear functions to grow hp
#' All other parameters are altered to maintain original proportions to hp
#'
#' @param growth A dataframe with the six fields:
#' Species - Name of the species consistent with other tables
#' max - Maximum plant height (m)
#' rate	- A constant describing the rate of growth for a Chapman Richards function
#' @param stn The number of the stratum
#' @param a The number of the record to be modelled
#' @param Param A parameter file
#' @param sp The name of the species being modelled
#' @param age The number of years since last fire
#' @return dataframe
#'

growPlants <- function(Param, a, sp, stn, growth, age)
{
  # Find species in param
  stid <- filter(Param, stratum == stn)
  spid <- filter(stid, value == sp)
  specName <- filter(stid, species == spid$species)

  # Collect vals
  hOrigr <- filter(specName, param == "hp")
  hOrig <- as.numeric(hOrigr$value)
  htOrigr <- filter(specName, param == "ht")
  htOrig <- as.numeric(htOrigr$value)
  heOrigr <- filter(specName, param == "he")
  heOrig <- as.numeric(heOrigr$value)
  bOrigr <- filter(specName, param == "hc")
  bOrig <- as.numeric(bOrigr$value)
  wOrigr <- filter(specName, param == "w")
  wOrig <- as.numeric(wOrigr$value)
  cOrigr <- filter(specName, param == "clumpDiameter")
  cOrig <- as.numeric(cOrigr$value)
  sOrigr <- filter(specName, param == "clumpSeparation")
  sOrig <- as.numeric(sOrigr$value)
  # Measure ratios
  baseR <- bOrig/hOrig
  htR <- htOrig/hOrig
  heR <- heOrig/hOrig
  widthR <- wOrig/hOrig
  clumping <- cOrig/pmin(wOrig,(hOrig-bOrig))
  openness <- sOrig/cOrig

  # Find species growth traits
  recGrowth <- filter(growth, record == a)
  stGrowth <- filter(recGrowth, stratum == stn)
  spGrowth <- filter(stGrowth, Species == sp)

  # Model growth
  height <- ifelse(!is.na(spGrowth$max),
                   spGrowth$max*(1-exp(-spGrowth$rate*age)),
                   ifelse(!is.na(spGrowth$aLin),
                          spGrowth$aLin*age+spGrowth$bLin,
                          0))
  ht <- height * htR
  base <- height * baseR
  he <- height * heR
  width <- height * widthR
  clumpDiameter <- clumping * pmin(width,(height-base))
  clumpSeparation <- openness*clumpDiameter

  # Build table
  dim <- as.data.frame(list('hp' = height, 'ht' = ht, 'he' = he, 'hc' = base,
                            'w'= width, 'clumpDiameter' = clumpDiameter, 'clumpSeparation' = clumpSeparation))
  return(dim)
}


#########################################

#' Updates a parameter file with the growth parameters modelled for a species
#'
#' Note: w is currently not changed due to an error:
#' "Error in .match_param(param, section, no.match.error = TRUE, single = TRUE) :
#' w matches more than one parameter"
#'
#' @param Param A parameter file
#' @param sp The name of the species being modelled
#' @param current A dataframe with one row and the 5 fields
#' hp, ht, he, hc, w, clumpDiameter, clumpSeparation; values in m
#' @return dataframe
#'

applyGrowth <- function(Param, sp, current)
{
  spid <- filter(Param, value == sp)
  stratumN <- spid$stratum
  SpeciesN <- spid$species

  Param <- ffm_set_species_param(Param, stratumN, SpeciesN, "hp", current$hp)
  Param <- ffm_set_species_param(Param, stratumN, SpeciesN, "ht", current$ht)
  Param <- ffm_set_species_param(Param, stratumN, SpeciesN, "he", current$he)
  Param <- ffm_set_species_param(Param, stratumN, SpeciesN, "hc", current$hc)
#  Param <- ffm_set_species_param(Param, stratumN, SpeciesN, "w", current$w)
  Param <- ffm_set_species_param(Param, stratumN, SpeciesN, "clumpDiameter", current$clumpDiameter)
  Param <- ffm_set_species_param(Param, stratumN, SpeciesN, "clumpSeparation", current$clumpSeparation)

  return(Param)
}

###################################################################################
#' Models the weighted mean of plant separation from time since fire for a named stratum
#'
#' Checks for alternate growth models
#' @param st The name of the stratum
#' @param a The number of the record to be modelled
#' @param cover A dataframe with the fields:
#' Species - Name of the species consistent with other tables
#' constant - Mean plant separation (m) that does not change with age
#' exp_a	- The first constant in an exponential function describing plant separation with tsf
#' exp_b	- The second constant in an exponential function describing plant separation with tsf
#' @param age The number of years since last fire
#' @return dataframe
#'

coverChange <- function(st, a, cover, Flora, age)
{
  coverA <- filter(cover, record == a)
  strat <- filter(coverA, sName == st)
  rec <- filter(Flora, record == a)
  stratF <- filter(rec, stratum == strat$stratum[1])
  spChange <- left_join(strat,stratF)
  nSp <- as.numeric(count(spChange))

  for(n in 1:nSp){
    spChange$wsep[n] <- ifelse(!is.na(spChange$exp_a[n]),
                               pmax(0,spChange$exp_a[n] * exp(spChange$exp_b[n] * age) * spChange$comp[n]),
                               ifelse(!is.na(spChange$aQ[n]),
                                      pmax(0,(spChange$aQ[n] * age^2 + spChange$bQ[n] * age + spChange$cQ[n])* spChange$comp[n]),
                                      ifelse(!is.na(spChange$aLin[n]),
                                             pmax(0,(spChange$aLin[n] * age + spChange$bLin[n])* spChange$comp[n]),
                                             1000)))
  }

  sep <- sum(spChange$wsep) / sum(spChange$comp)

  return(sep)
}

##############################################################################################
#' Models the weight of the o_horizon from time since fire and
#' updates param file
#'
#' Olsen negative exponential function
#'
#' @param base.params Parameter input table
#' @param growth A dataframe with the six fields:
#' Species - Name of the species consistent with other tables
#' max - Maximum plant height (m)
#' rate	- A constant describing the rate of growth for a Chapman Richards function
#' @param age The number of years since last fire
#' @return dataframe
#'

olsen <- function(base.params, growth, age)
{
  # Find growth curve
  ols <- filter(growth, Species == "O_horizon")

  # Model accumulation
  oHor <- max(4, ols$max*(1-exp(-ols$rate*age)))
  base.params <- ffm_set_site_param(base.params, "fuelLoad", oHor)

  return(base.params)
}

#################################################################################################
#' Models flammability dynamics from a weather set
#'
#' Grows the plant community,
#' modelling fire behaviour against a weather dataset for every age, then summarising all
#' behaviour into risk and impact statistics. Allows trees an older age than understorey to
#' represent lower severity fire histories
#'
#' @param base.params Parameter input table
#' @param suspended Name of the fuel in the suspended layer
#' @param density density of the wood in the suspended layer (kg/m3)
#' @param weather A dataframe with the four fields:
#' tm - Sequential numbering of the records
#' T - Air temperature (deg C)
#' W - Wind velocity (km/h)
#' DFMC - Dead fuel moisture content (proportion ODW)
#' @param growth A dataframe with the six fields:
#' Species - Name of the species consistent with other tables
#' max - Maximum plant height (m)
#' @param cover A dataframe with the fields:
#' Species - Name of the species consistent with other tables
#' constant - Mean plant separation (m) that does not change with age
#' exp_a	- The first constant in an exponential function describing plant separation with tsf
#' exp_b	- The second constant in an exponential function describing plant separation with tsf
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'DefaultSpeciesParams'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' base, he, ht, top & w - canopy dimensions for that species (m). he and ht are optional
#' clump - mean ratio of clump diameter to crown diameter
#' openness - ratio of gap to clump size
#' @param ageStep The number of years between age classes to be examined
#' @param firstAge The youngest age class to be examined
#' @param steps The number of age classes to be examined after firstAge
#' @param tAge The age of trees at the commencement of the run
#' @param test An isotherm for which the height will be reported
#' @param hKill The percentage of tree canopy above which a risk boundary is crossed (0 - 100)
#' @param jitters Number of repetitions for each row in the weather table
#' @param a The number of the record to be modelled
#' @return dataframe

fireDynamics <- function(base.params, weather, growth, cover, Flora, jitters = 5, ageStep = 5, firstAge = 1, steps = 10, tAge = 50, l = 0.1,
                         DefaultSpeciesParams = DefaultSpeciesParams, Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41, a, suspNS,
                         density = 300, test = 80, hKill = 90)
{
  age <- firstAge

  # AGE THE STAND
  nTable <- subset(base.params, param=="name")
  nSp <- as.numeric(count(nTable))
  strat <- filter(base.params, param == "levelName")
  nSt <- as.numeric(count(strat))
  nCanopy <- subset(nTable, stratum==nSt)
  nCsp <- as.numeric(count(nCanopy))
  nLow <- nSp-nCsp

  # Weight of the O-horizon
  base.params <- olsen(base.params, growth, age)

  # Structure of suspended dead material
  if(suspNS != ""){
    base.params <- susp(base.params, a, suspNS, Flora, growth,
                        DefaultSpeciesParams, age, density = density)
  }

  for (stNum in 1:nSt) {
    st <- strat$value[stNum]
    sep <- coverChange(st, a, cover, Flora, age)
    base.params <- ffm_set_stratum_param(base.params, stNum, "plantSeparation", sep)
    spList <- filter(nTable, stratum == stNum)
    n_a <- as.integer(spList$species[1])
    n_b <- as.integer(max(spList$species))
    nSusp <- as.integer((subset(base.params, value ==suspNS))$species[1])

    if (stNum < nSt) {
      for (spName in n_a:n_b)
        if(spList$value != suspNS){
          current <- growPlants(base.params, a, sp = nTable$value[spName], stn = stNum, growth, age)
          base.params <- applyGrowth(base.params, nTable$value[spName], current)
        }
    } else {
      for (spName in n_a:n_b)
        current <- growPlants(base.params, a, nTable$value[spName], stn = stNum, growth, tAge)
      base.params <- applyGrowth(base.params, nTable$value[spName], current)}
  }

  # RUN THE MODEL
  weatherSet(base.params, weather, jitters = jitters, l = l,
             Ms = Ms, Pm = Pm, Mr = Mr, Hs = Hs, Hr = Hr)
  res<-ffm_db_load("out_mc.db")

  # SUMMARISE RESULTS
  #Build tables
  surf <- surf(res$SurfaceResults)
  x <- stratum(res$FlameSummaries, res$Sites, res$ROS, surf)
  runs <- summary(x, surf)%>%
    mutate(time = ceiling(repId/jitters))
  IP <- repFlame(res$IgnitionPaths)

  # Find spread a: likelihood of fire spread
  S1 <- x%>%
    mutate(time = ceiling(repId/jitters),
           spread = ifelse(level == "Surface",
                           ifelse(oHorizon >=4,spread,0),
                           spread))%>%
    group_by(repId)%>%
    select(repId, time, spread)%>%
    summarise_all(max)#%>%

  # Find spread b: likelihood of crown loss
  runsB <- runs%>%
    right_join(S1)%>%
    mutate(lengthSurface = lengthSurface * spread,
           heightSurface = heightSurface * spread,
           angleSurface = angleSurface * spread,
           fh = pmax(heightSurface, heightPlant),
           fl = pmax(lengthSurface, lengthPlant))
  S2<-flora(runsB, IP, base.params, Test = test)%>%
    mutate(time = ceiling(repId/jitters),
           age = age,
           b1 = ifelse(b1>=hKill, 1, 0),
           b2 = ifelse(b2>=hKill, 1, 0),
           b3 = ifelse(b3>=hKill, 1, 0),
           b4 = ifelse(b4>=hKill, 1, 0),
           sc1 = ifelse(sc1>=hKill, 1, 0),
           sc2 = ifelse(sc2>=hKill, 1, 0),
           sc3 = ifelse(sc3>=hKill, 1, 0),
           sc4 = ifelse(sc4>=hKill, 1, 0))%>%
    group_by(time)%>%
    select(time, age, Height, b1, b2, b3, b4, sc1, sc2, sc3, sc4)%>%
    summarise_all(mean)

  S1 <- S1%>%
    group_by(time)%>%
    select(time, spread)%>%
    summarise_all(mean)

  Likelihood <- left_join(S1,S2)%>%
    mutate(Height = Height * spread,
           b1 = b1 * spread,
           b2 = b2 * spread,
           b3 = b3 * spread,
           b4 = b4 * spread,
           sc1 = sc1 * spread,
           sc2 = sc2 * spread,
           sc3 = sc3 * spread,
           sc4 = sc4 * spread)

  # Times table
  Times <- runs %>%
    select(repId, time, fh, fl, ros_kph,
           wind_kph, deadFuelMoistureProp, temperature, slope_degrees) %>%
    group_by(time) %>%
    summarize_all(mean) %>%
    right_join(Likelihood) %>%
    mutate(fh = fh * spread,
           fl = fl * spread,
           ros_kph = ros_kph * spread)

  cat(1, "time step complete")

  # Loop through remaining ages
  for (step in 1:steps) {
    age <- firstAge + ageStep * step
    tAge <- tAge + ageStep * step

    # Weight of the O-horizon
    base.params <- olsen(base.params, growth, age)

    # Structure of suspended dead material
    if(suspNS != ""){
      base.params <- susp(base.params, a, suspNS, Flora, growth,
                          DefaultSpeciesParams, age, density = density)
    }

    for (stNum in 1:nSt) {
      st <- strat$value[stNum]
      sep <- coverChange(st, a, cover, Flora, age)
      base.params <- ffm_set_stratum_param(base.params, stNum, "plantSeparation", sep)
      spList <- filter(nTable, stratum == stNum)
      n_a <- as.integer(spList$species[1])
      n_b <- as.integer(max(spList$species))
      nSusp <- as.integer((subset(base.params, value ==suspNS))$species[1])

      if (stNum < nSt) {
        for (spName in n_a:n_b)
          if(spList$value != suspNS){
            current <- growPlants(base.params, a, sp = nTable$value[spName], stn = stNum, growth, age)
            base.params <- applyGrowth(base.params, nTable$value[spName], current)
          }
      } else {
        for (spName in n_a:n_b)
          current <- growPlants(base.params, a, nTable$value[spName], stn = stNum, growth, tAge)
        base.params <- applyGrowth(base.params, nTable$value[spName], current)}
    }

    # RUN THE MODEL
    weatherSet(base.params, weather, jitters = jitters, l = l,
               Ms = Ms, Pm = Pm, Mr = Mr, Hs = Hs, Hr = Hr)
    res<-ffm_db_load("out_mc.db")

    # SUMMARISE RESULTS
    #Build tables
    surf <- surf(res$SurfaceResults)
    x <- stratum(res$FlameSummaries, res$Sites, res$ROS, surf)
    runs <- summary(x, surf)%>%
      mutate(time = ceiling(repId/jitters))
    IP <- repFlame(res$IgnitionPaths)

    # Find spread a: likelihood of fire spread
    S1 <- x%>%
      mutate(time = ceiling(repId/jitters),
             spread = ifelse(level == "Surface",
                             ifelse(oHorizon >=4,spread,0),
                             spread))%>%
      group_by(repId)%>%
      select(repId, time, spread)%>%
      summarise_all(max)#%>%

    # Find spread b: likelihood of crown loss
    runsB <- runs%>%
      right_join(S1)%>%
      mutate(lengthSurface = lengthSurface * spread,
             heightSurface = heightSurface * spread,
             angleSurface = angleSurface * spread,
             fh = pmax(heightSurface, heightPlant),
             fl = pmax(lengthSurface, lengthPlant))
    S2<-flora(runsB, IP, base.params, Test = test)%>%
      mutate(time = ceiling(repId/jitters),
             age = age,
             b1 = ifelse(b1>=hKill, 1, 0),
             b2 = ifelse(b2>=hKill, 1, 0),
             b3 = ifelse(b3>=hKill, 1, 0),
             b4 = ifelse(b4>=hKill, 1, 0),
             sc1 = ifelse(sc1>=hKill, 1, 0),
             sc2 = ifelse(sc2>=hKill, 1, 0),
             sc3 = ifelse(sc3>=hKill, 1, 0),
             sc4 = ifelse(sc4>=hKill, 1, 0))%>%
      group_by(time)%>%
      select(time, age, Height, b1, b2, b3, b4, sc1, sc2, sc3, sc4)%>%
      summarise_all(mean)

    S1 <- S1%>%
      group_by(time)%>%
      select(time, spread)%>%
      summarise_all(mean)

    Likelihood <- left_join(S1,S2)%>%
      mutate(Height = Height * spread,
             b1 = b1 * spread,
             b2 = b2 * spread,
             b3 = b3 * spread,
             b4 = b4 * spread,
             sc1 = sc1 * spread,
             sc2 = sc2 * spread,
             sc3 = sc3 * spread,
             sc4 = sc4 * spread)

    # Times table
    Times1 <- runs %>%
      select(repId, time, fh, fl, ros_kph,
             wind_kph, deadFuelMoistureProp, temperature, slope_degrees) %>%
      group_by(time) %>%
      summarize_all(mean) %>%
      right_join(Likelihood) %>%
      mutate(fh = fh * spread,
             fl = fl * spread,
             ros_kph = ros_kph * spread)

    Times <- rbind(Times, Times1)

    cat(step, "time steps complete")
  }
  return(Times)
}

########################################################################################
#' Models the packing of materials in a suspended layer
#'
#' Based on an Olsen negative exponential function;updates base.params
#'
#' @param base.params Parameter input table
#' @param a The number of the record to be modelled
#' @param suspended Name of the fuel in the suspended layer
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'DefaultSpeciesParams'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' base, he, ht, top & w - canopy dimensions for that species (m). he and ht are optional
#' clump - mean ratio of clump diameter to crown diameter
#' openness - ratio of gap to clump size
#' @param growth A dataframe with the six fields:
#' Species - Name of the species consistent with other tables
#' max - Maximum plant height (m)
#' rate	- A constant describing the rate of growth for a Chapman Richards function
#' @param DefaultSpeciesParams Leaf traits database
#' @param age The number of years since last fire
#' @param density Wood density (kg/m3)
#' @return dataframe
#'

susp <- function(base.params, a, suspNS, Flora, growth, DefaultSpeciesParams, age, density = 300)
{
  # Find growth curve
  olsS <- filter(growth, Species == suspNS)

  if(count(olsS) > 0) {

    FloraR <- filter(Flora, record == a)
    olsF <- filter(FloraR, species == suspNS)
    olsT <- filter(DefaultSpeciesParams, name == suspNS)
    spN <- filter(base.params, value == suspNS)

    # Model packing
    suspNS <- olsS$max*(1-exp(-olsS$rate*age))
    lengthS <- (0.6*((0.1 * suspNS) / (olsF$top * density))) / (pi * (olsT$leafThickness/2)^2)
    sepS <- mean(sqrt(sqrt(olsF$top/lengthS)^2*2),sqrt(olsF$top/lengthS))

    #Update base.params
    base.params <- ffm_set_species_param(base.params, olsF$stratum[1], spN$species[1], "leafSeparation", sepS)
  }

  return(base.params)
}
