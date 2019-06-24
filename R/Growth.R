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
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' base, he, ht, top & w - canopy dimensions for that species (m). he and ht are optional
#' clump - mean ratio of clump diameter to crown diameter
#' openness - ratio of gap to clump size
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

########################################################################################
#' Models the packing of materials in a suspended layer
#'
#' Based on an Olsen negative exponential function;updates base.params
#'
#' @param base.params Parameter input table
#' @param a The number of the record to be modelled
#' @param suspNS Name of the fuel in the suspended layer
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
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
#' @param default.species.params Leaf traits database
#' @param age The number of years since last fire
#' @param density Wood density (kg/m3)
#' @return dataframe
#'

susp <- function(base.params, a, suspNS, Flora, growth, default.species.params, age, density = 300)
{
  # Find growth curve
  olsS <- filter(growth, Species == suspNS)
  
  if(count(olsS) > 0) {
    
    FloraR <- filter(Flora, record == a)
    olsF <- filter(FloraR, species == suspNS)
    olsT <- filter(default.species.params, name == suspNS)
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

#########################################################################################

#' Updates the parameter file to a designated age

#'
#' @param base.params A parameter file
#' @param age The number of years since last fire
#' @param tAge The age of trees at the commencement of the run
#' @param growth A dataframe with the six fields:
#' Species - Name of the species consistent with other tables
#' max - Maximum plant height (m)
#' rate	- A constant describing the rate of growth for a Chapman Richards function
#' @param cover A dataframe with the fields:
#' Species - Name of the species consistent with other tables
#' constant - Mean plant separation (m) that does not change with age
#' exp_a	- The first constant in an exponential function describing plant separation with tsf
#' exp_b	- The second constant in an exponential function describing plant separation with tsf
#' @return dataframe
#'

ageCommunity <- function(base.params, age, tAge, growth, cover)
{
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
  
  # Structure of suspended dead material - still working out what to do with this
  suspNS <- ""
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
  
  return(base.params)
}
#################################################################################################