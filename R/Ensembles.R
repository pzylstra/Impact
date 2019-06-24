#' Updates parameter files with weather from a dataset,
#' then models fire from non-deterministic plant parameters
#'
#' @param base.params Parameter input table
#' @param weather A dataframe with the four fields:
#' tm - Sequential numbering of the records
#' T - Air temperature (deg C)
#' W - Wind velocity (km/h)
#' DFMC - Dead fuel moisture content (proportion ODW)
#' @param db.path Name of the exported database
#' @param jitters Number of repetitions
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height
#' @param updateProgress Progress bar for use in the dashboard
#' @return dataframe

weatherSet <- function(base.params, weather, db.path = "out_mc.db", jitters = 10, l = 0.1,
                       Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41,updateProgress = NULL)
{
  
  # This creates a reference to a Scala Database object
  # which will handle writing of model results to the
  # output file.First, clear any old dataframe.
  
  
  # Run the model, updating the base parameter table
  # with MC values at each iteration
  
  pbar <- txtProgressBar(max = max(weather$tm), style = 3)
  for (i in 1:max(weather$tm)) {
    ## Create database and delete the last part
    db.recreate <- i == 1
    
    # Read weather values from the table
    
    w <- weather$W[[i]]
    t <- weather$T[[i]]
    d <- max(0.01,min(0.199,weather$DFMC[[i]]))
    
    # Update parameter table
    tbl <- base.params %>%
      ffm_set_site_param("windSpeed", w, "km/h") %>%
      ffm_set_site_param("temperature", t, "degc") %>%
      ffm_set_site_param("deadFuelMoistureProp", d)
    
    Strata <- strata(base.params)
    Species <- species(base.params)
    
    if (jitters > 0) {
      for (j in 1:jitters) {
        tbl <- plantVar(tbl, Strata, Species, l = l,
                        Ms = Ms, Pm = Pm, Mr = Mr, Hs = Hs, Hr = Hr)
        # Run the model
        ffm_run(tbl, db.path, db.recreate = db.recreate)
      }
    }
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ",max(weather$tm) - i )
      updateProgress(detail = text)
    }
    setTxtProgressBar(pbar, i)
  }
  
  
  cat("Finished.  Output written to", db.path)
}
#######################################################################
#' Models flammability dynamics from a weather set
#'
#' Grows the plant community,
#' modelling fire behaviour against a weather dataset for every age, then summarising all
#' behaviour into risk and impact statistics. Allows trees an older age than understorey to
#' represent lower severity fire histories
#'
#' @param base.params Parameter input table
#' @param suspNS Name of the fuel in the suspended layer
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
#' species - the name of the species, which will call trait data from 'default.species.params'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' base, he, ht, top & w - canopy dimensions for that species (m). he and ht are optional
#' clump - mean ratio of clump diameter to crown diameter
#' openness - ratio of gap to clump size
#' @param default.species.params Leaf traits database
#' @param ageStep The number of years between age classes to be examined
#' @param firstAge The youngest age class to be examined
#' @param steps The number of age classes to be examined after firstAge
#' @param tAge The age of trees at the commencement of the run
#' @param test An isotherm for which the height will be reported
#' @param hKill The percentage of tree canopy above which a risk boundary is crossed (0 - 100)
#' @param jitters Number of repetitions for each row in the weather table
#' @param a The number of the record to be modelled
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height
#' @param updateProgress Progress bar for use in the dashboard
#' @return dataframe

fireDynamics <- function(base.params, weather, growth, cover, Flora, jitters = 5, ageStep = 5, firstAge = 1, steps = 10, tAge = 50, l = 0.1,
                         default.species.params = default.species.params, Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41, a, suspNS,
                         density = 300, test = 80, hKill = 90,updateProgress = NULL)
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
                        default.species.params, age, density = density)
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
  S1 <- x %>%
    mutate(time = ceiling(repId/jitters),
           spread = ifelse(level == "Surface",
                           ifelse(oHorizon >=4,spread,0),
                           spread))%>%
    group_by(repId)%>%
    select(repId, time, spread)%>%
    summarise_all(max)

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
                          default.species.params, age, density = density)
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
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ", steps - step)
      updateProgress(detail = text)
    }
  }
  return(Times)
}

####################################################################
#' Models fires using sites constructed from imported tables

#' @param site A dataframe with the six fields:
#' record - a unique, consecutively numbered identifier per site
#' site - the name of each record
#' slope - slope in degrees
#' wind - velocity in km/h
#' temp - ambient temperature deg. C
#' dfmc - moisture content of fine dead fuels in whole numbers (eg 0.1 for 10%)
#' oHorizon - weight in t/ha of fine dead organic material forming the O horizon
#' fline - the fireline length in m
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. Acceptable values
#' are t, f, or blank, where the outcome will be decided by the relative stratum heights.
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' hc, he, ht, hp & w - canopy dimensions for that species (m)
#' clump - mean ratio of clump diameter to crown diameter
#' openness - proportion of plant canopy occupied by gaps between clumps
#' @param traits A dataframe of plant traits, with the fields:
#' name - the name of the species
#' propDead - the proportion (0-1) of the foliage that is dead
#' leafForm - allowable values are 'Flat' or 'Round' for terete or needle leaves
#' leafThickness, leafWidth, leafLength, leafSeparation - leaf dimensions (m)
#' stemOrder - stem ramification
#' ignitionTemp - minimum temperature at which the leaf can ignite (deg C)
#' @param default.species.params Leaf traits database
#' @return dataframe

fireSet <- function(site, Structure, Flora, traits = default.species.params)
{
  # Number of sites to model
  sn <- max(site$record)
  
  # Enter trait database
  default.species.params <- traits
  
  # Initial record
  # Build parameter files
  param <- paramBuilder(site, Structure, Flora, 1)
  #Run the model
  ffm_run(param, "Behav.db")
  
  #Load and display the outputs
  res<-ffm_db_load("Behav.db")
  
  #Build tables
  surf <- surf(res$SurfaceResults)%>%
    mutate(fireNo = 1)
  x <- stratum(res$FlameSummaries, res$Sites, res$ROS, surf)%>%
    mutate(fireNo = 1)
  runs <- summary(x, surf)%>%
    mutate(cat = round(wind_kph,0),
           fireNo = 1)
  IP <- repFlame(res$IgnitionPaths)%>%
    mutate(fireNo = 1)
  
  if(sn >1){
    
    # Loop through records
    for (a in 2:sn) {
      
      # Build parameter files
      param <- paramBuilder(site, Structure, Flora, a)
      #Run the model
      ffm_run(param, "Behav.db")
      
      # Build a table with a row for each output. Consider a script for different levels of analysis (eg simple to comprehensive),
      # Then list these in If Then scenarios that can be identified in the initial function.
      # Build wider level functions for weather, growth etc that can call this one.
      
      #Load and display the outputs
      res<-ffm_db_load("Behav.db")
      
      #Build tables
      surfa <- surf(res$SurfaceResults)%>%
        mutate(fireNo = a)
      xa <- stratum(res$FlameSummaries, res$Sites, res$ROS, surf)%>%
        mutate(fireNo = a)
      runsa <- summary(xa, surfa)%>%
        mutate(cat = round(wind_kph,0),
               fireNo = a)
      IPa <- repFlame(res$IgnitionPaths)%>%
        mutate(fireNo = a)
      
      # Append to existing table
      surf <- rbind(surf, surfa)
      x <- rbind(x, xa)
      runs <- rbind(runs, runsa)
      IP <- rbind(IP, IPa)
    }
  }
  # Export csv files
  write.csv(surf,"Surface.csv")
  write.csv(IP,"IP.csv")
  write.csv(x,"All.csv")
  return(runs)
}

#####################################################################

#' Randomly modifies plant traits within defined ranges for non-deterministic predictions
#'
#' @param base.params Parameter input table
#' @param Strata - A dataframe of stratum properties as output by the function 'strata'
#' @param Species - A dataframe of species properties as output by the function 'species'
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height
#' @return dataframe

plantVar <- function (base.params, Strata, Species,
                      l = 0.1, Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41)
{
  
  tbl <- base.params
  
  tbl <- ffm_param_variance(tbl, max.prop = l, method = "uniform")
  SpeciesN <- 1
  SpeciesP <- 1
  
  StN <- as.numeric(count(Strata))
  for (si in 1:StN) {
    if (runif(1) <= Strata$cover[si]) {
      for (t in 1:Strata$speciesN[si]) {
        Mrand <- Pm * rtnorm(n = 1, mean = Species$lfmc[SpeciesN],
                             sd = Ms, a = Species$lfmc[SpeciesN]/Mr, b = Species$lfmc[SpeciesN] *
                               Mr)
        tbl <- tbl %>% ffm_set_species_param(si, SpeciesN,
                                             "liveLeafMoisture", Mrand)
        (SpeciesN = SpeciesN + 1)
      }
    }
    else {
      for (f in 1:Strata$speciesN[si]) {
        tbl <- tbl %>% ffm_set_species_param(si, SpeciesN,
                                             "liveLeafMoisture", 100)
        SpeciesN = SpeciesN + 1
      }
    }
    for (p in 1:Strata$speciesN[si]) {
      peak <- rtnorm(n = 1, mean = Species$hp[SpeciesP],
                     sd = Hs, a = Species$hp[SpeciesP]/Hr, b = Species$hp[SpeciesP] *
                       Hr)
      tbl <- tbl %>%
        ffm_set_species_param(si, SpeciesP, "hp", peak) %>%
        ffm_set_species_param(si, SpeciesP, "ht", peak * Species$htR[SpeciesP]) %>%
        ffm_set_species_param(si, SpeciesP, "he", peak * Species$heR[SpeciesP]) %>%
        ffm_set_species_param(si, SpeciesP, "hc", peak *
                                Species$hcR[SpeciesP])
      SpeciesP = SpeciesP + 1
    }
  }
  
  return(tbl)
}
