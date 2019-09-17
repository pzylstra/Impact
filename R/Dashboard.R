#' Models probabilistic fire behaviour
#' @param base.params Input parameter file
#' @param db.path Name of the exported database
#' @param jitters Number of repetitions for each row in the weather table
#' @param slope Mean slope (deg)
#' @param slopeSD Standard deviation of slope
#' @param slopeRange Truncates variability by +/- mean * range
#' @param temp Mean ambient temperature (deg.C)
#' @param tempSD Standard deviation of temp
#' @param tempRange Truncates variability by +/- mean * range
#' @param DFMC Mean DFMC (%ODW)
#' @param DFMCSD Standard deviation of DFMC
#' @param DFMCRange Truncates variability by +/- mean * range
#' @param wind Mean wind velocity (km/h)
#' @param windSD Standard deviation of wind velocity
#' @param windRange Truncates variability by +/- mean * range
#' @param moistureMultiplier Multiplies all LFMC values by this number
#' @param moistureSD Standard deviation of moisture
#' @param moistureRange Truncates variability by +/- mean * range
#' @param heightSD Standard deviation of plant height
#' @param heightRange Truncates variability by +/- mean * range
#' @param leafVar Variation around input leaf dimensions, equivalent to l
#' @param updateProgress Progress bar for use in the dashboard
#' @return dataframe
#' @export

probFire <- function(base.params, db.path = "out_mc.db", jitters,
                   slope, slopeSD, slopeRange, temp, tempSD, tempRange,
                   DFMC, DFMCSD, DFMCRange, wind, windSD, windRange,
                   moistureMultiplier, moistureSD, moistureRange,
                   heightSD, heightRange, leafVar,updateProgress = NULL) {
  
  # Collect original descriptors
  Strata <- strata(base.params)
  Species <- species(base.params)
  
  #Set input limits
  DFMCRange <- pmax(1.0001, DFMCRange)
  
  
  pbar <- txtProgressBar(max = jitters, style = 3)
  if (jitters > 2) {
    for (j in 1:jitters) {
      db.recreate <- j == 1
      #Update parameters
      s <- rtnorm(n = 1, mean = slope, sd = slopeSD,
                  a = slope-(slopeRange/2), b = slope+(slopeRange/2))
      t <- rtnorm(n = 1, mean = temp, sd = tempSD,
                  a = temp-(tempRange/2), b = temp+(tempRange/2)) 
      d <- rtnorm(n = 1, mean = DFMC, sd = DFMCSD,
                  a = pmax(0.02, DFMC-(DFMCRange/2)), b = pmin(0.199,DFMC+(DFMCRange/2))) 
      w <- rtnorm(n = 1, mean = wind, sd = windSD,
                  a = wind-(windRange/2), b = wind+(windRange/2))   
      
      base.params <- base.params %>%
        ffm_set_site_param("slope", s, "deg") %>%
        ffm_set_site_param("temperature", t) %>%
        ffm_set_site_param("deadFuelMoistureProp", d) %>%
        ffm_set_site_param("windSpeed", w)
      
      base.params <- plantVar(base.params, Strata, Species, 
                              l = leafVar, Ms = moistureSD, Pm = moistureMultiplier, Mr = moistureRange, 
                              Hs = heightSD, Hr = heightRange)
      ffm_run(base.params, db.path, db.recreate = db.recreate)
      Sys.sleep(0.25)
      ####UpdateProgress
      if (is.function(updateProgress)) {
        text <- paste0("Number of remaining steps is ", jitters - j)
        updateProgress(detail = text)
      }
    }
    
    setTxtProgressBar(pbar, j)
    
  }
  
  else
  {
    print("Probabilistic analysis requires a minimum of 3 replicates")
  }
  
}


######
## drivers function
#####

#' Models fire behaviour across ranged variables
#' @param base.params Input parameter file
#' @param db.path Name of the exported database
#' @param jitters Number of repetitions for each row in the weather table
#' @param windMin Lowest wind velocity tested (km/h)
#' @param windReps Number of wind speeds tested
#' @param windStep Gap (km/h) between wind steps
#' @param moistureMultiplier Multiplies all LFMC values by this number
#' @param moistureSD Standard deviation of moisture
#' @param moistureRange Truncates variability by +/- mean * range
#' @param heightSD Standard deviation of plant height
#' @param heightRange Truncates variability by +/- mean * range
#' @param leafVar Variation around input leaf dimensions, equivalent to l
#' @param updateProgress Progress bar for use in the dashboard
#' @export

drivers <- function(base.params, db.path = "out_mc.db", jitters, windMin, windReps, windStep,
                    moistureMultiplier, moistureSD, moistureRange, heightSD, heightRange, 
                    leafVar,updateProgress = NULL) {
  
  # Collect original descriptors
  Strata <- strata(base.params)
  Species <- species(base.params)
  
  #Range environmental values
  slopes <- c(-10,0,10,20,30)
  DFMCs <- c(0.05, 0.1, 0.15)
  winds <- seq(windMin, (windReps*windStep+windMin), windStep)
  
  #Dataframe of orthogonal combinations
  dat <- expand.grid(slope = slopes, DFMC = DFMCs, wind = winds)
  Niter <- nrow(dat) * jitters
  
  #Loop through combinations
  pbar <- txtProgressBar(max = Niter, style = 3)
  for (i in 1:Niter) {
    set <- ceiling(i / jitters)
    db.recreate <- i == 1
    s <- dat[set, "slope"]
    d <- dat[set, "DFMC"]
    w <- dat[set, "wind"]
    
    #Update environmental parameters if on a new row
    if (set > ceiling((i-1) / jitters)) {
      base.params <- base.params %>%
        ffm_set_site_param("slope", s, "deg") %>%
        ffm_set_site_param("temperature", 30) %>%
        ffm_set_site_param("deadFuelMoistureProp", d) %>%
        ffm_set_site_param("windSpeed", w)
    }
    
    #  for (j in 1:jitters) {
    #Randomise plant parameters
    base.params <- plantVar(base.params, Strata, Species, l = leafVar, Ms = moistureSD, 
                            Pm = moistureMultiplier, Mr = moistureRange, 
                            Hs = heightSD, Hr = heightRange)
    ffm_run(base.params, db.path, db.recreate = db.recreate)
    
    #  }
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ", Niter - i)
      updateProgress(detail = text)
    }
    setTxtProgressBar(pbar, i)
  }
}

######
## drivers function with species-specific changes
#####

#' Models fire behaviour across ranged variables
#' @param base.params Input parameter file
#' @param db.path Name of the exported database
#' @param jitters Number of repetitions for each row in the weather table
#' @param windMin Lowest wind velocity tested (km/h)
#' @param windReps Number of wind speeds tested
#' @param windStep Gap (km/h) between wind steps
#' @param moistureMultiplier Multiplies all LFMC values by this number
#' @param moistureSD Standard deviation of moisture
#' @param moistureRange Truncates variability by +/- mean * range
#' @param leafVar Variation around input leaf dimensions, equivalent to l
#' @param Variation A database of plant variability in traits, with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' Hs - Standard deviation of plant height variations
#' Hr - Truncates plant height variability by +/- Hr * height
#' @param updateProgress Progress bar for use in the dashboard
#' @export

driversS <- function(base.params, db.path = "out_mc.db", jitters, windMin, windReps, windStep,
                    moistureMultiplier, moistureSD, moistureRange, Variation,
                    leafVar,updateProgress = NULL) {
  
  # Collect original descriptors
  Strata <- strata(base.params)
  Species <- species(base.params)
  
  #Range environmental values
  slopes <- c(-10,0,10,20,30)
  DFMCs <- c(0.05, 0.1, 0.15)
  winds <- seq(windMin, (windReps*windStep+windMin), windStep)
  
  #Dataframe of orthogonal combinations
  dat <- expand.grid(slope = slopes, DFMC = DFMCs, wind = winds)
  Niter <- nrow(dat) * jitters
  
  #Loop through combinations
  pbar <- txtProgressBar(max = Niter, style = 3)
  for (i in 1:Niter) {
    set <- ceiling(i / jitters)
    db.recreate <- i == 1
    s <- dat[set, "slope"]
    d <- dat[set, "DFMC"]
    w <- dat[set, "wind"]
    
    #Update environmental parameters if on a new row
    if (set > ceiling((i-1) / jitters)) {
      base.params <- base.params %>%
        ffm_set_site_param("slope", s, "deg") %>%
        ffm_set_site_param("temperature", 30) %>%
        ffm_set_site_param("deadFuelMoistureProp", d) %>%
        ffm_set_site_param("windSpeed", w)
    }
    
    #  for (j in 1:jitters) {
    #Randomise plant parameters
    base.params <- plantVarS(base.params, Strata, Species, l = leafVar, Ms = moistureSD, 
                            Pm = moistureMultiplier, Mr = moistureRange, Variation)
    ffm_run(base.params, db.path, db.recreate = db.recreate)
    
    #  }
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ", Niter - i)
      updateProgress(detail = text)
    }
    setTxtProgressBar(pbar, i)
  }
}
