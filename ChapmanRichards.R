#' Models plant height from time since fire
#'
#' Standard Chapman-Richards negative exponential function used to grow hp
#' All other parameters are altered to maintain original proportions to hp
#'
#' @param growth A dataframe with the six fields:
#' Species - Name of the species consistent with other tables
#' max - Maximum plant height (m)
#' rate	- A constant describing the rate of growth for a Chapman Richards function
#' @param Param A parameter file
#' @param sp The name of the species being modelled
#' @param age The number of years since last fire
#' @return dataframe
#'

chapmanRichards <- function(Param, sp, growth, age)
{
  # Find species in param
  spid <- filter(Param, value == sp)
  specName <- filter(Param, species == spid$species)

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
  spGrowth <- filter(growth, Species == sp)

  # Model growth
  height <- spGrowth$max*(1-exp(-spGrowth$rate*age))
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
