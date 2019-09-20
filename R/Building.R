# Builds the dataframe site.meta from input tables

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
  site.meta$value[6] <- site$litter[a]
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

# Builds the dataframe strata.meta from input tables

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

# Builds the dataframe species.values from input tables


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

# Builds the dataframe species.units from input tables

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
#' comp - % composition of that species in the stratum. If absent, all species will be considered equally
#' hc, he, ht, hp & w - canopy dimensions for that species (m)
#' clump - mean ratio of clump diameter to crown diameter
#' openness - proportion of plant canopy occupied by gaps between clumps
#' @param default.species.params Leaf traits database
#' @param a The record number for which to build the table
#' @export

paramBuilder <- function(site, Structure, Flora, default.species.params, a)
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
  
  # Find weighted mean of crown widths per stratum
  species.values$weightedW <- species.values$composition * species.values$w
  ww <- species.values%>%
    select(stratum, composition, weightedW)%>%
    group_by(stratum) %>%
    summarize_all(sum) %>%
    mutate(mw = weightedW/composition)
  
  # Ensure separation is greater than weighted mean of crown widths
  for (stNum in 1:max(ww$stratum)) {
    sep <- filter(strata.meta, stratum == stNum)
    param <- ffm_set_stratum_param(param, stNum, "plantSeparation", 
                                   pmax(as.numeric(sep$value[2]),ww$mw[stNum]))
  }
  
    # Change species param leafForm to characters
    default.species.params$leafForm <- as.character(default.species.params$leafForm)
    
    # Fill empty traits
    param <- ffm_complete_params(param,default.species.params)
    
    return(param)
}
