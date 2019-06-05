
#' Builds the dataframe site.meta from input tables
#'
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
#' @param a The site number for which to build the table

siteBuilder <- function(site, Structure, a)
{

  # CREATE site.meta
  paramSite <- c('overlapping','overlapping','overlapping','overlapping','overlapping', 'fuelLoad',
                 'meanFuelDiameter','meanFinenessLeaves', 'fireLineLength', 'slope', 'temperature',
                 'deadFuelMoistureProp', 'windSpeed')


  unitsSite <- c(NA, NA, NA, NA, NA, 't/ha', 'm', 'm', 'm', 'deg', 'degC', NA, 'km/h')
  site.meta <- data.frame(matrix(NA, nrow=13, ncol = 5))

  names(site.meta) <- c("stratum", "species", "param", "value", "units")
  site.meta$param <- paramSite
  site.meta$units <- unitsSite
  # ENTER VARIABLES
  site.meta$value[1] <- ifelse(Structure$ns_e[a]=='t',"near surface, elevated, overlapped",
                               ifelse(Structure$ns_e[a]=='f',"near surface, elevated, not overlapped",
                                      "near surface, elevated, automatic"))
  site.meta$value[2] <- ifelse(Structure$ns_m[a]=='t',"near surface, midstorey, overlapped",
                               ifelse(Structure$ns_m[a]=='f',"near surface, midstorey, not overlapped",
                                      "near surface, midstorey, automatic"))
  site.meta$value[3] <- ifelse(Structure$e_m[a]=='t',"elevated, midstorey, overlapped",
                               ifelse(Structure$e_m[a]=='f',"elevated, midstorey, not overlapped",
                                      "elevated, midstorey, automatic"))
  site.meta$value[4] <- ifelse(Structure$e_c[a]=='t',"elevated, canopy, overlapped",
                               ifelse(Structure$e_c[a]=='f',"elevated, canopy, not overlapped",
                                      "elevated, canopy, automatic"))
  site.meta$value[5] <- ifelse(Structure$m_c[a]=='t',"midstorey, canopy, overlapped",
                               ifelse(Structure$m_c[a]=='f',"midstorey, canopy, not overlapped",
                                      "midstorey, canopy, automatic"))
  site.meta$value[6] <- site$oHorizon[a]
  site.meta$value[7] <- 0.005
  site.meta$value[8] <- 0.00025
  site.meta$value[9] <- site$fLine[a]
  site.meta$value[10] <- site$slope[a]
  site.meta$value[11] <- site$temp[a]
  site.meta$value[12] <- site$dfmc[a]
  site.meta$value[13] <- site$wind[a]

  return(site.meta)
}

#####################################################################

#' Builds the dataframe strata.meta from input tables
#'
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - the name of each record
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. Acceptable values
#' are t, f, or blank, where the outcome will be decided by the relative stratum heights.
#' @param a The record number for which to build the table

strataBuilder <- function(Structure, Flora, a)
{
  # Collect subsets for record
  st <- Structure[Structure$record==a,]
  fl <- Flora[Flora$record==a,]

  # CREATE strata.meta
  strata.meta <- data.frame(matrix(NA, nrow=8, ncol = 5))

  names(strata.meta) <- c("stratum", "species", "param", "value", "units")

  # Fill levelNames
  strata.meta$param[1] <- "levelName"
  strata.meta$value[1] <- "near surface"
  strata.meta$param[2] <- "plantSeparation"
  strata.meta$value[2] <- st$NS[1]
  strata.meta$units[2] <- "m"

  strata.meta$param[3] <- "levelName"
  strata.meta$value[3] <- "elevated"
  strata.meta$param[4] <- "plantSeparation"
  strata.meta$value[4] <- st$El[1]
  strata.meta$units[4] <- "m"

  strata.meta$param[5] <- "levelName"
  strata.meta$value[5] <- "midstorey"
  strata.meta$param[6] <- "plantSeparation"
  strata.meta$value[6] <- st$Mid[1]
  strata.meta$units[6] <- "m"

  strata.meta$param[7] <- "levelName"
  strata.meta$value[7] <- "canopy"
  strata.meta$param[8] <- "plantSeparation"
  strata.meta$value[8] <- st$Can[1]
  strata.meta$units[8] <- "m"

  # Find empty strata and remove them
  deleteVector <- vector()
  for (ro in 2:8) {
    if(is.na(strata.meta$value[ro])){
      deleteVector<- c(deleteVector, ro-1, ro)
    }
  }
  if(length(deleteVector)>0){
    strata.meta <- strata.meta[ -deleteVector, ]}

  # Number strata
  rows <- as.numeric(count(strata.meta))
  num <- 0.5
  for (ro in 1:rows) {
    numb <- ceiling(num)
    strata.meta$stratum[ro] <- numb
    num = num + 0.5
  }

  return(strata.meta)
}
#####################################################################

#' Builds the dataframe species.values from input tables
#'
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'DefaultSpeciesParams'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' base, he, ht, top & w - canopy dimensions for that species (m). he and ht are optional
#' clump - mean ratio of clump diameter to crown diameter
#' openness - ratio of gap to clump size
#' @param site A dataframe with the six fields:
#' record - a unique, consecutively numbered identifier per site
#' site - the name of each record
#' slope - slope in degrees
#' wind - velocity in km/h
#' temp - ambient temperature deg. C
#' dfmc - moisture content of fine dead fuels in whole numbers (eg 0.1 for 10%)
#' oHorizon - weight in t/ha of fine dead organic material forming the O horizon
#' fline - the fireline length in m
#' @param a The site number for which to build the table

speciesBuilder <- function(Flora, site, a)
{
  # Collect subsets for site
  fl <- Flora[Flora$record==a,]
  si <- site[site$record==a,]

  # CREATE species.values
  ro <- as.numeric(count(fl))
  species.values <- data.frame(matrix(NA, nrow=ro, ncol = 13))

  names(species.values) <- c("stratum", "species", "name", "clumpDiameter", "clumpSeparation", "composition",
                             "deadLeafMoisture", "hc", "he", "hp", "ht", "w", "liveLeafMoisture")

  # Enter values
  species.values$stratum <- as.numeric(fl$stratum)
  species.values$name <- as.character(fl$species)
  species.values$species <- as.numeric(c(1:ro))
  species.values$liveLeafMoisture <- fl$moisture
  species.values$hc <- fl$base
  species.values$hp <- fl$top
  species.values$he <- fl$he
  species.values$ht <- fl$ht
  species.values$w <- fl$w
  species.values$clumpDiameter <- pmin((fl$top-fl$base),fl$w)*fl$clump
  species.values$clumpSeparation <- fl$openness*species.values$clumpDiameter
  species.values$composition <- as.numeric(fl$comp)
  species.values$deadLeafMoisture <- si$dfmc

  # Calculate gaps
  for(a in 1:ro) {
    if(is.na(species.values$he[a])){
      species.values$he[a] <- species.values$hc[a]
    }
    if(is.na(species.values$ht[a])){
      species.values$ht[a] <- species.values$hp[a]
    }
    if(is.na(species.values$liveLeafMoisture[a])){
      species.values$liveLeafMoisture[a] <- 1
    }
  }

  return(species.values)
}

#####################################################################

#' Builds the dataframe species.units from input tables
#'
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'DefaultSpeciesParams'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' base, he, ht, top & w - canopy dimensions for that species (m). he and ht are optional
#' clump - mean ratio of clump diameter to crown diameter
#' openness - ratio of gaps between clumps to clump size
#' @param a The record number for which to build the table

unitBuilder <- function(Flora, a)
{
  # Collect subsets for site
  fl <- Flora[Flora$record==a,]

  # CREATE species.units
  ro <- as.numeric(count(fl))
  species.units <- data.frame(matrix(NA, nrow=ro, ncol = 13))

  names(species.units) <- c("stratum", "species", "name", "clumpDiameter", "clumpSeparation", "composition",
                            "deadLeafMoisture", "hc", "he", "hp", "ht", "w", "liveLeafMoisture")

  # Enter values
  species.units$stratum <- fl$stratum
  species.units$species <- c(1:ro)
  species.units$clumpDiameter <- "m"
  species.units$clumpSeparation <- "m"
  species.units$hc <- "m"
  species.units$he <- "m"
  species.units$hp <- "m"
  species.units$ht <- "m"
  species.units$w <- "m"

  return(species.units)
}

#####################################################################

#' Constructs parameter files from imported tables

#' @param site A dataframe with the six fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
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
#' species - the name of the species, which will call trait data from 'DefaultSpeciesParams'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition of that species in the stratum. If absent, all species will be considered equally
#' hc, he, ht, hp & w - canopy dimensions for that species (m)
#' clump - mean ratio of clump diameter to crown diameter
#' openness - proportion of plant canopy occupied by gaps between clumps
#' @param DefaultSpeciesParams Leaf traits database
#' @param a The record number for which to build the table

paramBuilder <- function(site, Structure, Flora, DefaultSpeciesParams=DefaultSpeciesParams, a)
{
  # Construct component tables
  site.meta <- siteBuilder(site, Structure, a)
  strata.meta <- strataBuilder(Structure, Flora, a)
  species.values <- speciesBuilder(Flora, site, a)
  species.units <- unitBuilder(Flora, a)

  # Build parameter file
  components <- list(site.meta, strata.meta, species.values, species.units)
  names(components) <- c("site.meta", "strata.meta", "species.values", "species.units")
  param <- ffm_assemble_table(components)

  # Adjust sep to width
  species.values$weightedW <- species.values$composition * species.values$w
  ww <- species.values%>%
    select(stratum, composition, weightedW)%>%
    group_by(stratum) %>%
    summarize_all(sum) %>%
    mutate(mw = weightedW/composition)

  for (stNum in 1:max(ww$stratum)) {
    sep <- filter(strata.meta, stratum == stNum)
    param <- ffm_set_stratum_param(param, stNum, "plantSeparation", 
                                   pmax(sep$value[2],ww$mw[stNum]))

  # Change species param leafForm to characters
  DefaultSpeciesParams$leafForm <- as.character(DefaultSpeciesParams$leafForm)

  # Fill empty traits
  param <- ffm_complete_params(param,DefaultSpeciesParams)

  return(param)
}

#####################################################################

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
#' species - the name of the species, which will call trait data from 'DefaultSpeciesParams'
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
#' @param DefaultSpeciesParams Leaf traits database
#' @param siteV An optional dataframe with the same fields as 'site', but providing standard deviations of numerical values
#' @param commV An optional dataframe with the same fields as 'community', but providing standard deviations of numerical values
#' @param traitV An optional dataframe with the same fields as 'traits', but providing ranges of numerical values
#' @return dataframe
#' @export 3 csv files
#' @examples fireSet(site, Structure, Flora, traits)

fireSet <- function(site, Structure, Flora, traits = DefaultSpeciesParams)
{
  # Number of sites to model
  sn <- max(site$record)

  # Enter trait database
  DefaultSpeciesParams <- traits

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

#' Randomly modifies plant traits within defined ranges for non-deterministic predictions,
#' then runs the model
#'
#' uses Monte Carlo, for set wind, slope, DFMC & temperature
#' LFMC, heights and leaf traits randomly varied within defined ranges
#'
#' @param base-params Parameter input table
#' @param out.db Name of the exported database
#' @param jitters Number of repetitions
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height
#' @examples run_ffm_ND(Param, "MC", jitters = 100, s = 0, l = 0.1,
#' tempm = 30, Dm = 0.07, windl = 0, windst = 50,
#' Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41)

run_ffm_plantVar <- function (base.params, out.db = "out_mc.db", jitters = 100,
                              l = 0.1, Ms = 0.01,
                              Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41)
{
  strataN <- strata(base.params)
  species <- species(base.params)
  db <- ffm_create_database(out.db, delete.existing = TRUE,
                            use.transactions = TRUE)
  tbl <- base.params
  pbar <- txtProgressBar(max = jitters, style = 3)
  for (i in 1:jitters) {
    tbl <- ffm_param_variance(tbl, max.prop = l, method = "uniform")
    SpeciesN <- 1
    SpeciesP <- 1
    si <- 1
    StN <- as.numeric(count(strataN))
    for (loop in 1:StN) {
      if (runif(1) <= strataN$cover[si]) {
        for (t in 1:strataN$speciesN[si]) {
          Mrand <- Pm * rtnorm(n = 1, mean = species$lfmc[SpeciesN],
                               sd = Ms, a = species$lfmc[SpeciesN]/Mr, b = species$lfmc[SpeciesN] *
                                 Mr)
          tbl <- tbl %>% ffm_set_species_param(si, SpeciesN,
                                               "liveLeafMoisture", Mrand)
          (SpeciesN = SpeciesN + 1)
        }
      }
      else {
        for (f in 1:strataN$speciesN[si]) {
          tbl <- tbl %>% ffm_set_species_param(si, SpeciesN,
                                               "liveLeafMoisture", 100)
          SpeciesN = SpeciesN + 1
        }
      }
      for (p in 1:strataN$speciesN[si]) {
        peak <- rtnorm(n = 1, mean = species$hp[SpeciesP],
                       sd = Hs, a = species$hp[SpeciesP]/Hr, b = species$hp[SpeciesP] *
                         Hr)
        tbl <- tbl %>%
          ffm_set_species_param(si, SpeciesP, "hp", peak) %>%
          ffm_set_species_param(si, SpeciesP, "ht", peak * species$htR[SpeciesP]) %>%
          ffm_set_species_param(si, SpeciesP, "he", peak * species$heR[SpeciesP]) %>%
          ffm_set_species_param(si, SpeciesP, "hc", peak *
                                  species$hcR[SpeciesP])
        SpeciesP = SpeciesP + 1
      }
      si = si + 1
    }
    ffm_run(tbl, db)
    setTxtProgressBar(pbar, i)
  }
  db$close()
  cat("Finished.  Output written to", out.db)
}

#####################################################################

#' Randomly modifies plant traits within defined ranges for non-deterministic predictions
#'
#' @param base-params Parameter input table
#' @param Strata - A dataframe of stratum properties as output by the function 'strata'
#' @param Species - A dataframe of species properties as output by the function 'species'
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height
#' @examples run_ffm_ND(Param, "MC", l = 0.1,
#' Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41)
plantVar <- function (base.params, out.db = "out_mc.db", Strata, Species,
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

############################################################

#' Updates parameter files with weather from a dataset,
#' then models fire from non-deterministic plant parameters
#'
#' @param base.params Parameter input table
#' @param weather A dataframe with the four fields:
#' tm - Sequential numbering of the records
#' T - Air temperature (deg C)
#' W - Wind velocity (km/h)
#' DFMC - Dead fuel moisture content (proportion ODW)
#' @param Strata - A dataframe of stratum properties as output by the function 'strata'
#' @param Species - A dataframe of species properties as output by the function 'species'
#' @param out.db Name of the exported database
#' @param jitters Number of repetitions
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height
#' @return dataframe

weatherSet <- function(base.params, weather, out.db = "out_mc.db", jitters = 10, l = 0.1,
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
        tbl <- plantVar(tbl, out.db = "out_mc.db", Strata, Species, l = l,
                        Ms = Ms, Pm = Pm, Mr = Mr, Hs = Hs, Hr = Hr)
        # Run the model
        ffm_run(tbl, out.db, db.recreate = db.recreate)
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


  cat("Finished.  Output written to", out.db)
}


#' Builds the dataframe site.meta from input tables
#'
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
#' @param a The site number for which to build the table

siteBuilder <- function(site, Structure, a)
{

  # CREATE site.meta
  paramSite <- c('overlapping','overlapping','overlapping','overlapping','overlapping', 'fuelLoad',
                 'meanFuelDiameter','meanFinenessLeaves', 'fireLineLength', 'slope', 'temperature',
                 'deadFuelMoistureProp', 'windSpeed')


  unitsSite <- c(NA, NA, NA, NA, NA, 't/ha', 'm', 'm', 'm', 'deg', 'degC', NA, 'km/h')
  site.meta <- data.frame(matrix(NA, nrow=13, ncol = 5))

  names(site.meta) <- c("stratum", "species", "param", "value", "units")
  site.meta$param <- paramSite
  site.meta$units <- unitsSite
  # ENTER VARIABLES
  site.meta$value[1] <- ifelse(Structure$ns_e[a]=='t',"near surface, elevated, overlapped",
                               ifelse(Structure$ns_e[a]=='f',"near surface, elevated, not overlapped",
                                      "near surface, elevated, automatic"))
  site.meta$value[2] <- ifelse(Structure$ns_m[a]=='t',"near surface, midstorey, overlapped",
                               ifelse(Structure$ns_m[a]=='f',"near surface, midstorey, not overlapped",
                                      "near surface, midstorey, automatic"))
  site.meta$value[3] <- ifelse(Structure$e_m[a]=='t',"elevated, midstorey, overlapped",
                               ifelse(Structure$e_m[a]=='f',"elevated, midstorey, not overlapped",
                                      "elevated, midstorey, automatic"))
  site.meta$value[4] <- ifelse(Structure$e_c[a]=='t',"elevated, canopy, overlapped",
                               ifelse(Structure$e_c[a]=='f',"elevated, canopy, not overlapped",
                                      "elevated, canopy, automatic"))
  site.meta$value[5] <- ifelse(Structure$m_c[a]=='t',"midstorey, canopy, overlapped",
                               ifelse(Structure$m_c[a]=='f',"midstorey, canopy, not overlapped",
                                      "midstorey, canopy, automatic"))
  site.meta$value[6] <- site$oHorizon[a]
  site.meta$value[7] <- 0.005
  site.meta$value[8] <- 0.00025
  site.meta$value[9] <- site$fLine[a]
  site.meta$value[10] <- site$slope[a]
  site.meta$value[11] <- site$temp[a]
  site.meta$value[12] <- site$dfmc[a]
  site.meta$value[13] <- site$wind[a]

  return(site.meta)
}

#####################################################################

#' Builds the dataframe strata.meta from input tables
#'
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - the name of each record
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. Acceptable values
#' are t, f, or blank, where the outcome will be decided by the relative stratum heights.
#' @param a The record number for which to build the table

strataBuilder <- function(Structure, Flora, a)
{
  # Collect subsets for record
  st <- Structure[Structure$record==a,]
  fl <- Flora[Flora$record==a,]

  # CREATE strata.meta
  strata.meta <- data.frame(matrix(NA, nrow=8, ncol = 5))

  names(strata.meta) <- c("stratum", "species", "param", "value", "units")

  # Fill levelNames
  strata.meta$param[1] <- "levelName"
  strata.meta$value[1] <- "near surface"
  strata.meta$param[2] <- "plantSeparation"
  strata.meta$value[2] <- st$NS[1]
  strata.meta$units[2] <- "m"

  strata.meta$param[3] <- "levelName"
  strata.meta$value[3] <- "elevated"
  strata.meta$param[4] <- "plantSeparation"
  strata.meta$value[4] <- st$El[1]
  strata.meta$units[4] <- "m"

  strata.meta$param[5] <- "levelName"
  strata.meta$value[5] <- "midstorey"
  strata.meta$param[6] <- "plantSeparation"
  strata.meta$value[6] <- st$Mid[1]
  strata.meta$units[6] <- "m"

  strata.meta$param[7] <- "levelName"
  strata.meta$value[7] <- "canopy"
  strata.meta$param[8] <- "plantSeparation"
  strata.meta$value[8] <- st$Can[1]
  strata.meta$units[8] <- "m"

  # Find empty strata and remove them
  deleteVector <- vector()
  for (ro in 2:8) {
    if(is.na(strata.meta$value[ro])){
      deleteVector<- c(deleteVector, ro-1, ro)
    }
  }
  if(length(deleteVector)>0){
    strata.meta <- strata.meta[ -deleteVector, ]}

  # Number strata
  rows <- as.numeric(count(strata.meta))
  num <- 0.5
  for (ro in 1:rows) {
    numb <- ceiling(num)
    strata.meta$stratum[ro] <- numb
    num = num + 0.5
  }

  return(strata.meta)
}
#####################################################################

#' Builds the dataframe species.values from input tables
#'
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'DefaultSpeciesParams'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' base, he, ht, top & w - canopy dimensions for that species (m). he and ht are optional
#' clump - mean ratio of clump diameter to crown diameter
#' openness - ratio of gap to clump size
#' @param site A dataframe with the six fields:
#' record - a unique, consecutively numbered identifier per site
#' site - the name of each record
#' slope - slope in degrees
#' wind - velocity in km/h
#' temp - ambient temperature deg. C
#' dfmc - moisture content of fine dead fuels in whole numbers (eg 0.1 for 10%)
#' oHorizon - weight in t/ha of fine dead organic material forming the O horizon
#' fline - the fireline length in m
#' @param a The site number for which to build the table

speciesBuilder <- function(Flora, site, a)
{
  # Collect subsets for site
  fl <- Flora[Flora$record==a,]
  si <- site[site$record==a,]

  # CREATE species.values
  ro <- as.numeric(count(fl))
  species.values <- data.frame(matrix(NA, nrow=ro, ncol = 13))

  names(species.values) <- c("stratum", "species", "name", "clumpDiameter", "clumpSeparation", "composition",
                             "deadLeafMoisture", "hc", "he", "hp", "ht", "w", "liveLeafMoisture")

  # Enter values
  species.values$stratum <- as.numeric(fl$stratum)
  species.values$name <- as.character(fl$species)
  species.values$species <- as.numeric(c(1:ro))
  species.values$liveLeafMoisture <- fl$moisture
  species.values$hc <- fl$base
  species.values$hp <- fl$top
  species.values$he <- fl$he
  species.values$ht <- fl$ht
  species.values$w <- fl$w
  species.values$clumpDiameter <- pmin((fl$top-fl$base),fl$w)*fl$clump
  species.values$clumpSeparation <- fl$openness*species.values$clumpDiameter
  species.values$composition <- as.numeric(fl$comp)
  species.values$deadLeafMoisture <- si$dfmc

  # Calculate gaps
  for(a in 1:ro) {
    if(is.na(species.values$he[a])){
      species.values$he[a] <- species.values$hc[a]
    }
    if(is.na(species.values$ht[a])){
      species.values$ht[a] <- species.values$hp[a]
    }
    if(is.na(species.values$liveLeafMoisture[a])){
      species.values$liveLeafMoisture[a] <- 1
    }
  }

  return(species.values)
}

#####################################################################

#' Builds the dataframe species.units from input tables
#'
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'DefaultSpeciesParams'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' base, he, ht, top & w - canopy dimensions for that species (m). he and ht are optional
#' clump - mean ratio of clump diameter to crown diameter
#' openness - ratio of gaps between clumps to clump size
#' @param a The record number for which to build the table

unitBuilder <- function(Flora, a)
{
  # Collect subsets for site
  fl <- Flora[Flora$record==a,]

  # CREATE species.units
  ro <- as.numeric(count(fl))
  species.units <- data.frame(matrix(NA, nrow=ro, ncol = 13))

  names(species.units) <- c("stratum", "species", "name", "clumpDiameter", "clumpSeparation", "composition",
                            "deadLeafMoisture", "hc", "he", "hp", "ht", "w", "liveLeafMoisture")

  # Enter values
  species.units$stratum <- fl$stratum
  species.units$species <- c(1:ro)
  species.units$clumpDiameter <- "m"
  species.units$clumpSeparation <- "m"
  species.units$hc <- "m"
  species.units$he <- "m"
  species.units$hp <- "m"
  species.units$ht <- "m"
  species.units$w <- "m"

  return(species.units)
}

#####################################################################

#' Constructs parameter files from imported tables

#' @param site A dataframe with the six fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
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
#' species - the name of the species, which will call trait data from 'DefaultSpeciesParams'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition of that species in the stratum. If absent, all species will be considered equally
#' hc, he, ht, hp & w - canopy dimensions for that species (m)
#' clump - mean ratio of clump diameter to crown diameter
#' openness - proportion of plant canopy occupied by gaps between clumps
#' @param DefaultSpeciesParams Leaf traits database
#' @param a The record number for which to build the table

paramBuilder <- function(site, Structure, Flora, DefaultSpeciesParams=DefaultSpeciesParams, a)
{
  #Construct component tables
  site.meta <- siteBuilder(site, Structure, a)
  strata.meta <- strataBuilder(Structure, Flora, a)
  species.values <- speciesBuilder(Flora, site, a)
  species.units <- unitBuilder(Flora, a)
  
  #Build parameter file
  components <- list(site.meta, strata.meta, species.values, species.units)
  names(components) <- c("site.meta", "strata.meta", "species.values", "species.units")
  param <- ffm_assemble_table(components)
  
  #Adjust sep to width
  species.values$weightedW <- species.values$composition * species.values$w
  ww <- species.values %>% 
    select(stratum, composition, weightedW) %>% 
    group_by(stratum) %>% 
    summarize_all(sum) %>% 
    mutate(mw = weightedW/composition)
  
 for (stNum in 1:max(ww$stratum)) {
    sep <- filter(strata.meta, stratum == stNum)
    param <- ffm_set_stratum_param(param, stNum, "plantSeparation", 
                                   pmax(sep$value[2],ww$mw[stNum]))
  }
  DefaultSpeciesParams$leafForm <- as.character(DefaultSpeciesParams$leafForm)
  param <- ffm_complete_params(param)
  return(param)
}

#####################################################################

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
#' species - the name of the species, which will call trait data from 'DefaultSpeciesParams'
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
#' @param DefaultSpeciesParams Leaf traits database
#' @param siteV An optional dataframe with the same fields as 'site', but providing standard deviations of numerical values
#' @param commV An optional dataframe with the same fields as 'community', but providing standard deviations of numerical values
#' @param traitV An optional dataframe with the same fields as 'traits', but providing ranges of numerical values
#' @return dataframe
#' @export 3 csv files
#' @examples fireSet(site, Structure, Flora, traits)

fireSet <- function(site, Structure, Flora, traits = DefaultSpeciesParams)
{
  # Number of sites to model
  sn <- max(site$record)

  # Enter trait database
  DefaultSpeciesParams <- traits

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

#' Randomly modifies plant traits within defined ranges for non-deterministic predictions,
#' then runs the model
#'
#' uses Monte Carlo, for set wind, slope, DFMC & temperature
#' LFMC, heights and leaf traits randomly varied within defined ranges
#'
#' @param base-params Parameter input table
#' @param out.db Name of the exported database
#' @param jitters Number of repetitions
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height
#' @examples run_ffm_ND(Param, "MC", jitters = 100, s = 0, l = 0.1,
#' tempm = 30, Dm = 0.07, windl = 0, windst = 50,
#' Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41)

run_ffm_plantVar <- function (base.params, out.db = "out_mc.db", jitters = 100,
                              l = 0.1, Ms = 0.01,
                              Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41)
{
  strataN <- strata(base.params)
  species <- species(base.params)
  db <- ffm_create_database(out.db, delete.existing = TRUE,
                            use.transactions = TRUE)
  tbl <- base.params
  pbar <- txtProgressBar(max = jitters, style = 3)
  for (i in 1:jitters) {
    tbl <- ffm_param_variance(tbl, max.prop = l, method = "uniform")
    SpeciesN <- 1
    SpeciesP <- 1
    si <- 1
    StN <- as.numeric(count(strataN))
    for (loop in 1:StN) {
      if (runif(1) <= strataN$cover[si]) {
        for (t in 1:strataN$speciesN[si]) {
          Mrand <- Pm * rtnorm(n = 1, mean = species$lfmc[SpeciesN],
                               sd = Ms, a = species$lfmc[SpeciesN]/Mr, b = species$lfmc[SpeciesN] *
                                 Mr)
          tbl <- tbl %>% ffm_set_species_param(si, SpeciesN,
                                               "liveLeafMoisture", Mrand)
          (SpeciesN = SpeciesN + 1)
        }
      }
      else {
        for (f in 1:strataN$speciesN[si]) {
          tbl <- tbl %>% ffm_set_species_param(si, SpeciesN,
                                               "liveLeafMoisture", 100)
          SpeciesN = SpeciesN + 1
        }
      }
      for (p in 1:strataN$speciesN[si]) {
        peak <- rtnorm(n = 1, mean = species$hp[SpeciesP],
                       sd = Hs, a = species$hp[SpeciesP]/Hr, b = species$hp[SpeciesP] *
                         Hr)
        tbl <- tbl %>%
          ffm_set_species_param(si, SpeciesP, "hp", peak) %>%
          ffm_set_species_param(si, SpeciesP, "ht", peak * species$htR[SpeciesP]) %>%
          ffm_set_species_param(si, SpeciesP, "he", peak * species$heR[SpeciesP]) %>%
          ffm_set_species_param(si, SpeciesP, "hc", peak *
                                  species$hcR[SpeciesP])
        SpeciesP = SpeciesP + 1
      }
      si = si + 1
    }
    ffm_run(tbl, db)
    setTxtProgressBar(pbar, i)
  }
  db$close()
  cat("Finished.  Output written to", out.db)
}

#####################################################################

#' Randomly modifies plant traits within defined ranges for non-deterministic predictions
#'
#' @param base-params Parameter input table
#' @param Strata - A dataframe of stratum properties as output by the function 'strata'
#' @param Species - A dataframe of species properties as output by the function 'species'
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height
#' @examples run_ffm_ND(Param, "MC", l = 0.1,
#' Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41)
plantVar <- function (base.params, out.db = "out_mc.db", Strata, Species,
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

############################################################

#' Updates parameter files with weather from a dataset,
#' then models fire from non-deterministic plant parameters
#'
#' @param base.params Parameter input table
#' @param weather A dataframe with the four fields:
#' tm - Sequential numbering of the records
#' T - Air temperature (deg C)
#' W - Wind velocity (km/h)
#' DFMC - Dead fuel moisture content (proportion ODW)
#' @param Strata - A dataframe of stratum properties as output by the function 'strata'
#' @param Species - A dataframe of species properties as output by the function 'species'
#' @param out.db Name of the exported database
#' @param jitters Number of repetitions
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height
#' @return dataframe

weatherSet <- function(base.params, weather, out.db = "out_mc.db", jitters = 10, l = 0.1,
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
        tbl <- plantVar(tbl, out.db = "out_mc.db", Strata, Species, l = l,
                        Ms = Ms, Pm = Pm, Mr = Mr, Hs = Hs, Hr = Hr)
        # Run the model
        ffm_run(tbl, out.db, db.recreate = db.recreate)
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


  cat("Finished.  Output written to", out.db)
}

