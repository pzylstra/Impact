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
#' @export

weatherSet <- function(base.params, weather, db.path = "out_mc.db", jitters = 10, l = 0.1,
                       Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41,updateProgress = NULL)
{
  # Collect initial veg descriptors
  
  Strata <- strata(base.params)
  Species <- species(base.params)
  
  # Run the model, updating the base parameter table
  # with MC values at each iteration
  
  pbar <- txtProgressBar(max = max(weather$tm), style = 3)
  for (i in 1:max(weather$tm)) {
    
    # Read weather values from the table
    
    w <- weather$W[[i]]
    t <- weather$T[[i]]
    d <- max(0.01,min(0.199,weather$DFMC[[i]]))
    
    # Update parameter table
    tbl <- base.params %>%
      ffm_set_site_param("windSpeed", w, "km/h") %>%
      ffm_set_site_param("temperature", t, "degc") %>%
      ffm_set_site_param("deadFuelMoistureProp", d)
    
    if (jitters > 0) {
      for (j in 1:jitters) {
        ## Create database and delete the last part
        db.recreate <- (i * j) == 1
        
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

#' Updates parameter files with weather from a dataset,
#' then models fire from non-deterministic plant parameters
#' using plantVarS to modify individual species with their own measured variability
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
#' @param Variation A database of plant variability in traits, with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' Hs - Standard deviation of plant height variations
#' Hr - Truncates plant height variability by +/- Hr * height
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. Acceptable values
#' are t, f, or blank, where the outcome will be decided by the relative stratum heights.
#' nsR, eR, mR, cR - maximum species richness recorded for each stratum
#' @param updateProgress Progress bar for use in the dashboard
#' @return dataframe
#' @export

weatherSetS <- function(base.params, weather, Variation, Structure, a, db.path = "out_mc.db", jitters = 10, l = 0.1,
                        Ms = 0.01, Pm = 1, Mr = 1.001, updateProgress = NULL)
{
  pbar <- txtProgressBar(max = max(weather$tm), style = 3)
  
  # Read weather values from the table
  for (i in 1:max(weather$tm)) {
    w <- weather$W[[i]]
    t <- weather$T[[i]]
    d <- max(0.01,min(0.199,weather$DFMC[[i]]))
    
    # Select species for a random point in the forest and import the weather parameters
    tbl <- specPoint(base.params, Structure, a) %>%
      ffm_set_site_param("windSpeed", w, "km/h") %>%
      ffm_set_site_param("temperature", t, "degc") %>%
      ffm_set_site_param("deadFuelMoistureProp", d)
    
    Strata <- strata(tbl)
    Species <- species(tbl)
    
    if (jitters > 0) {
      for (j in 1:jitters) {
        # Recreate database on first run, then add following runs to this
        db.recreate <- (i * j) == 1
        
        # Vary plant traits for each species within their range
        tbl <- plantVarS(tbl, Strata, Species, Variation, a, l = l,
                         Ms = Ms, Pm = Pm, Mr = Mr)
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

####################################################################
#' Models fires using sites constructed from imported tables
#' 
#' Private function in development

#' @param site A dataframe with the six fields:
#' record - a unique, consecutively numbered identifier per site
#' site - the name of each record
#' slope - slope in degrees
#' wind - velocity in km/h
#' temp - ambient temperature deg. C
#' dfmc - moisture content of fine dead fuels in whole numbers (eg 0.1 for 10%)
#' litter - weight in t/ha of fine dead organic material forming the O horizon
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
#' @param base.params Parameter input table
#' @param Strata Strata descriptor table output by the function 'strata'
#' @param Species Species descriptor table output by the function 'species'
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height
#' @return dataframe
#' @export

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
        ffm_set_species_param(si, SpeciesP, "hc", peak *Species$hcR[SpeciesP]) %>%
        ffm_set_species_param(si, SpeciesP, "w", peak * Species$wR[SpeciesP])
      SpeciesP = SpeciesP + 1
    }
  }
  
  return(tbl)
}

#####################################################################

#' Randomly modifies plant traits within defined ranges for non-deterministic predictions
#' Differs from plantVar by modifying individual species by their own rules
#' @param base.params Parameter input table
#' @param Strata Strata descriptor table output by the function 'strata'
#' @param Species Species descriptor table output by the function 'species'
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

plantVarS <- function (base.params, Strata, Species, Variation, a, l = 0.1, Ms = 0.01, Pm = 1, Mr = 1.001)
{
  
  tbl <- base.params
  # Vary leaf traits
  tbl <- ffm_param_variance(tbl, max.prop = l, method = "uniform")
  SpeciesN <- 1
  SpeciesP <- 1
  
  #Filter variation table to record
  varRec <- Variation[Variation$record ==a, ]
  
  # Loop through plant strata
  
  StN <- as.numeric(count(Strata))
  
  for (si in 1:StN) {
    
    # Vary leaf moisture to randomly place points in the community and decide which plants will be present.
    # Where the random number is > stratum cover, all species are made too moist to burn and excluded from modelling.
    # Otherwise, species are varied by Ms & Mr, and multiplied by Pm
    
    if (runif(1) <= Strata$cover[si]) {
      for (t in 1:Strata$speciesN[si]) {
        Mrand <- Pm * rtnorm(n = 1, mean = Species$lfmc[SpeciesN],
                             sd = Ms, a = Species$lfmc[SpeciesN]/Mr, b = Species$lfmc[SpeciesN] * Mr)
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
    
    # Modify plant dimansions for each species within the stratum
    
    for (p in 1:Strata$speciesN[si]) {
      Hr <- varRec$Hr[SpeciesP]
      peak <- rtnorm(n = 1, mean = Species$hp[SpeciesP],
                     sd = varRec$Hs[SpeciesP], a = Species$hp[SpeciesP]/Hr, b = Species$hp[SpeciesP] *
                       Hr)
      tbl <- tbl %>%
        ffm_set_species_param(si, SpeciesP, "hp", peak) %>%
        ffm_set_species_param(si, SpeciesP, "ht", peak * Species$htR[SpeciesP]) %>%
        ffm_set_species_param(si, SpeciesP, "he", peak * Species$heR[SpeciesP]) %>%
        ffm_set_species_param(si, SpeciesP, "hc", peak *Species$hcR[SpeciesP]) %>%
        ffm_set_species_param(si, SpeciesP, "w", peak * Species$wR[SpeciesP])
      SpeciesP = SpeciesP + 1
    }
  }
  
  return(tbl)
}

#####################################################################
#' Selects random species from the available list, weighted by their frequency
#' Modifies a parameter table to the shortened list
#'
#' @param base.params A parameter file
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. Acceptable values
#' are t, f, or blank, where the outcome will be decided by the relative stratum heights.
#' nsR, eR, mR, cR - maximum species richness recorded for each stratum
#' @return dataframe
#' @export

specPoint <- function(base.params, Structure, a)
{
  
  Species <- species(base.params)
  
  # Add temporary fields
  Species$wComp = 0
  Species$include = 1
  
  # Count strata
  StN <- as.numeric(max(base.params$stratum[!is.na(base.params$stratum)]))
  
  
  # For each stratum, identify the species being considered, then choose how many of these
  # will be modelled from a range set by the recorded maximum point richness of each stratum
  
  richList <- c(as.numeric(Structure[Structure$record == a, ]$nsR), as.numeric(Structure[Structure$record == a, ]$eR),
                as.numeric(Structure[Structure$record == a, ]$mR), as.numeric(Structure[Structure$record == a, ]$cR))
  richList <- richList[!is.na(richList)]
  
  for (StratNo in 1:StN) {
    
    SpL <- as.numeric(min(Species[Species$st == StratNo, ]$sp))
    SpU <- as.numeric(max(Species[Species$st == StratNo, ]$sp))
    SpN <- SpU-SpL+1
    
    # Species richness for point in the stratum
    R <- richList[StratNo]
    choose <- round(runif(n=1)*(min(R,SpN)-1),0)+1
    
    # Limit stratum species list to a random selection weighted by species occurrence
    for (a in SpL:SpU) {
      Species$wComp[a] = runif(n=1)*Species$comp[a]
    }
    # Identify unneeded records
    low <- Rfast::nth(Species[Species$st == StratNo, ]$wComp, choose, descending = TRUE)
    
    for (sp in SpL:SpU) {
      Species$include[sp] = if (Species$wComp[sp] < low) { 0 }
      else {  1 }
    }
  }
  Species <- Species%>%
    mutate(species = as.character(sp))
  Species[,"new"] <- cumsum(Species$include)
  
  param <- left_join(base.params, Species, by="species")%>%
    subset(include != 0 | is.na(include))%>%
    mutate(species = new)%>%
    select(stratum, species, param, value, units)
  rownames(param) <- seq(length=nrow(param))
  
  return(param)
}


#####################################################################
#' Models fire behaviour in each species of a parameter table,
#' summarising the combustibility from the length of flame divided
#' by the length of segment ignited
#'
#' @param base.params A parameter file
#' @return dataframe
#' @export

spComb <- function(base.params)
{
  specflam <- function (base.params, st, sp) 
  {
    SP <- base.params %>% filter(stratum == st & species == sp) %>%
      add_row(stratum = st, param = "levelName", value = "canopy") %>%
      add_row(stratum = st, param = "plantSeparation", value = 100)
    end <- base.params %>% filter(is.na(stratum))
    tab <- rbind(SP, end)
    
    # Adjust plant heights to a 0.5m base
    base <- min(as.numeric((tab[which(tab$param == "hc"), ])$value[1]),
                as.numeric((tab[which(tab$param == "he"), ])$value[1]))
    dif <- 0.5-base
    hc <- as.numeric((tab[which(tab$param == "hc"), ])$value[1])+dif
    he <- as.numeric((tab[which(tab$param == "he"), ])$value[1])+dif
    ht <- as.numeric((tab[which(tab$param == "ht"), ])$value[1])+dif
    hp <- as.numeric((tab[which(tab$param == "hp"), ])$value[1])+dif
    tab <- tab %>% 
      ffm_set_species_param(stratum.id = st, species.id = sp, "hc", hc) %>% 
      ffm_set_species_param(stratum.id = st, species.id = sp, "he", he) %>% 
      ffm_set_species_param(stratum.id = st, species.id = sp, "ht", ht) %>% 
      ffm_set_species_param(stratum.id = st, species.id = sp, "hp", hp) %>%
      ffm_set_site_param("windSpeed", 2) %>%
      ffm_set_site_param("temperature", 30) %>%
      ffm_set_site_param("deadFuelMoistureProp", 0.05) %>%
      ffm_set_site_param("slope", 0)
    
    # Run for surface fuel loads from 4 to 20 t/ha
    for (i in 4:20) {
      db.recreate <- i == 4
      tab <- ffm_set_site_param(tab, "fuelLoad", i)
      ffm_run(params = tab, db.path = "out_mc.db", db.recreate = db.recreate)
    }
    
    res<-ffm_db_load("out_mc.db")
    IP <- repFlame(res$IgnitionPaths) %>%
      mutate(rat = flameLength/length)
    
    name <- (tab[which(tab$param == "name"), ])$value[1]
    flam <- round(mean(IP$rat),2)
    out <- as.data.frame(list('Stratum' = st, 'Species' = sp, 
                              'name' = name,  'Combustibility'=flam))
    return(out)
  }
  
  #Create list of species & strata
  can <- base.params %>%
    filter(!is.na(stratum)) %>%
    select(stratum, species)
  candidates <- distinct(can) %>%
    filter(!is.na(species))
  
  n <- as.numeric(nrow(candidates))
  combustibility <- specflam(base.params, candidates$stratum[1], candidates$species[1])
  
  for (j in 2:n){
    summ <- specflam(base.params, candidates$stratum[j], candidates$species[j])
    combustibility <- rbind(combustibility, summ)
  }
  return(combustibility)
}