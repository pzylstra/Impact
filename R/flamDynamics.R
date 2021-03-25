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
#' species - Name of the species consistent with other tables
#' max - Maximum plant height (m)
#' @param cover A dataframe with the fields:
#' species - Name of the species consistent with other tables
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
#' @export

fireDynamics <- function(base.params, weather, growth, cover, Flora, jitters = 5, ageStep = 5, firstAge = 1, steps = 10, tAge = 50, l = 0.1,
                         default.species.params = default.species.params, Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41, a, suspNS = "suspNS",
                         density = 300, test = 80, hKill = 90, db.path = "out_mc.db", updateProgress = NULL)
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
  susPres <- dplyr::filter(nTable, value == suspNS)
  
  # Weight of the O-horizon
  base.params <- olson(base.params, growth, age)
  
  # Structure of suspended dead material
  if(!is.na(susPres$value[1])){
    base.params <- susp(base.params, a, suspNS = suspNS, Flora, growth,
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
             Ms = Ms, Pm = Pm, Mr = Mr, Hs = Hs, Hr = Hrv, db.path = db.path)
  res<-ffm_db_load(db.path)
  
  # SUMMARISE RESULTS
  #Build tables
  surf <- surf(res$SurfaceResults)
  x <- stratum(res$FlameSummaries, res$Sites, res$ROS, surf)
  runs <- summary(x, surf)%>%
    mutate(time = ceiling(repId/jitters))
  IP <- repFlame(res$IgnitionPaths)
  
  # Limit fire spread to litter >= 4t/ha, find max of all strata
  S1 <- x %>%
    mutate(time = ceiling(repId/jitters),
           spread = ifelse(level == "Surface",
                           ifelse(litter >=4,spread,0),
                           spread))%>%
    group_by(repId)%>%
    select(repId, time, spread)%>%
    summarise_all(max)
  
  # Find likelihood of crown loss
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
    #    group_by(time)%>%
    select(time, age, Height, b1, b2, b3, b4, sc1, sc2, sc3, sc4)#%>%
  #    summarise_all(mean)
  
  S1 <- S1%>%
    #    group_by(time)%>%
    select(time, spread)#%>%
  #    summarise_all(mean)
  
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
    #    group_by(time) %>%
    #    summarize_all(mean) %>%
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
    base.params <- olson(base.params, growth, age)
    
    # Structure of suspended dead material
    if(!is.na(susPres$value[1])){
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
               Ms = Ms, Pm = Pm, Mr = Mr, Hs = Hs, Hr = Hr, db.path = db.path)
    res<-ffm_db_load(db.path)
    
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
                             ifelse(litter >=4,spread,0),
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
      #      group_by(time)%>%
      select(time, age, Height, b1, b2, b3, b4, sc1, sc2, sc3, sc4)#%>%
    #      summarise_all(mean)
    
    S1 <- S1%>%
      #      group_by(time)%>%
      select(time, spread)#%>%
    #      summarise_all(mean)
    
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
      #      group_by(time) %>%
      #      summarize_all(mean) %>%
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
