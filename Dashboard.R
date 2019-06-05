probIn <- function(base.params, out.db = "out_mc.db", jitters,
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
      
      base.params <- plantVar(base.params, out.db = "out_mc.db", 
                              Strata, Species, l = leafVar, Ms = moistureSD, Pm = moistureMultiplier, Mr = moistureRange, 
                              Hs = heightSD, Hr = heightRange)
      ffm_run(base.params, out.db,db.recreate)
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

drivers <- function(base.params, out.db = "out_mc.db", jitters, windMin, windReps, windStep,
                    moistureMultiplier, moistureSD, 
                    moistureRange, heightSD, heightRange, leafVar,updateProgress = NULL) {
  
  # N.B. This is based on the old FRaME version without Michael's recent changes to database format. 
  # You may need to adjust these as I don't really know how they work yet
  
  
  # Collect original descriptors
  Strata <- strata(base.params)
  Species <- species(base.params)
  
  #Range environmental values
  slopes <- c(-10,0,10,20,30)
  DFMCs <- c(0.05, 0.1, 0.15)
  winds <- seq(windMin, (windReps*windStep+windMin), windStep)
  
  #Dataframe of orthogonal combinations
  dat <- expand.grid(slope = slopes, DFMC = DFMCs, wind = winds)
  Niter <- nrow(dat)
  
  #Loop through combinations
  pbar <- txtProgressBar(max = Niter, style = 3)
  for (i in 1:Niter) {
    db.recreate <- i == 1
    s <- dat[i, "slope"]
    d <- dat[i, "DFMC"]
    w <- dat[i, "wind"]
    
    #Update environmental parameters
    base.params <- base.params %>%
      ffm_set_site_param("slope", s, "deg") %>%
      ffm_set_site_param("temperature", 30) %>%
      ffm_set_site_param("deadFuelMoistureProp", d) %>%
      ffm_set_site_param("windSpeed", w)
    
    for (j in 1:jitters) {
      #Randomise plant parameters
      base.params <- plantVar(base.params, out.db = "out_mc.db", 
                              Strata, Species, l = leafVar, Ms = moistureSD, 
                              Pm = moistureMultiplier, Mr = moistureRange, 
                              Hs = heightSD, Hr = heightRange)
      ffm_run(base.params, out.db, db.recreate)
      
    }
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ", Niter - i)
      updateProgress(detail = text)
    }
    setTxtProgressBar(pbar, i)
  }
}
