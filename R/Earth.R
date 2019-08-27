#' Soil heating
#'
#' Calculates the dynamic heating of a soil at 1cm increments, to 5cm depth
#'
#' Assumes all to be A horizon, with constant, uncompacted density
#'
#' Utilises the output tables from 'threat' and 'radiation', and adds to these
#' the Reynolds Number, heat transfer coefficients, Newton's convective energy transfer coefficient,
#' and the temperature of the object each second.
#'
#' Reynolds Number utilises a standard formulation (e.g. Gordon, N. T., McMahon, T. A. & Finlayson, B. L.
#' Stream hydrology: an introduction for ecologists. (Wiley, 1992))
#'
#' Convective heat transfer coefficients use the widely adopted formulations of
#' Williams, F. A. Urban and wildland fire phenomenology. Prog. Energy Combust. Sci. 8, 317–354 (1982),
#' and Drysdale, D. An introduction to fire dynamics. (John Wiley and Sons, 1985)
#' utilising a Prandtl number of 0.7.
#'
#' Heat is transferred into the earth using Fourier's Law. Spread continues for a period after the
#' passage of the fire front, equal to the duration of the surface flame, as determined using
#' Burrows, N. D. Flame residence times and rates of weight loss of eucalypt forest fuel particles.
#' Int. J. Wildl. Fire 10, 137–143 (2001).
#'
#' Default temperature of the resident flame is the average of the surface maximums in
#' Cawson, J. G., Nyman, P., Smith, H. G., Lane, P. N. J. & Sheridan, G. J.
#' How soil temperatures during prescribed burning affect soil water repellency,
#' infiltration and erosion. Geoderma 278, 12–22 (2016).
#'
#' Heating area is set to 1m2, flat, with a characteristic length of 1m
#'
#' Broad germination and seed death temperatures are based on
#' Auld, T. D. & O’Connel, M. A. Predicting patterns of post‐fire germination in
#' 35 eastern Australian Fabaceae. Aust. J. Ecol. 16, 53–70 (1991).
#'
#'
#' @param Surf The dataframe 'runs' exported from Monte Carlos as 'Summary.csv'
#' @param Plant The dataframe 'IP' exported from Monte Carlos as 'IP.csv'.
#' @param diameter Diameter of the surface fuels burning (mm)
#' @param surface Temperature at the surface of the soil, under burning fuels
#' @param percentile defines which heating statistics are used for each second, from 0 (min) to 1 (max)
#' @param RH The relative humidity (0-1)
#' @param moisture The proportion oven-dry weight of moisture in the bark and wood
#' @param distance The furthest horizontal distance between the flame origin and the point (m)
#' @param trail Number of seconds to continue modelling after the front has passed
#' @param var The angle in degrees that the plume spreads above/below a central vector;defaults to 10
#' @param Pressure Sea level atmospheric pressure (hPa)
#' @param Altitude Height above sea level (m)
#' @param texture Soil texture. Allowable values are: "sand", "loamy sand", "sandy loam", "sandy clay loam",
#' "sand clay", "loam", "clay loam", "silt loam", "clay", "silty clay", "silty clay loam", "silt"
#' @param peat Organic proportion of the soil
#' @param grain Allowable values are "fine" or "coarse"
#' @param unfrozen Proportion of soil unfrozen, between 0 and 1
#' @param soilTemp The starting temperature under the ground (deg C)
#' @param updateProgress Progress bar for use in the dashboard
#' @return dataframe
#' @export


soil <- function(Surf, Plant, diameter = 6, surface = 677, percentile = 0.95, RH = 0.2,
                 moisture = 0.1, distance = 50, trail = 600, var = 10, Pressure = 1013.25,
                 Altitude = 0, texture = "sand", peat = 0, 
                 grain = "fine", unfrozen = 1, soilTemp = 25,updateProgress = NULL)
  
{
  # Collect step distance, time, and total distance
  residence <- 0.871*diameter^1.875
  ROS <- mean(Surf$ros_kph)/3.6
  Ta <- round(distance/ROS+residence)
  TIME <- Ta + trail
  Horiz <- distance
  
  # Description of the protection
  densityD <- denSoil(texture)
  mass <- 0.01 * densityD
  R <- sqrt(1/pi)
  startM <- moisture
  
  #Starting values
  Ca <- threat(Surf, Plant, Horiz, Height=0, var, Pressure, Altitude)%>%
    mutate(t = 1,
           #Convective transfer
           Re = (Plume_velocity*Density)/viscosity,
           h = ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888),
           #Incoming heat from surface
           tempS = ifelse(Horiz <0, pmax(tempAir, surface), tempAir),
           qc = h * (tempS - soilTemp),
           att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
           qr = 0.86*qr*att,
           Qi = pmax(0, qc)+qr,
           
           #1ST CM________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*mass,
           # Energy removed by current water quantity
           drain = ifelse(soilTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiA = max(Qi-drain,0),
           #Thermal values
           cpSoilA = cpSoil((soilTemp+273.15), texture, peat, moisture),
           saturationA = satSoil(texture, moisture),
           kSoilA = kSoil(texture, saturationA),
           # Fourier conduction
           fourierA = ifelse(Horiz>0, pmin(QiA,(kSoilA * pmax(0,tempS - soilTemp)) / 0.01),
                             (kSoilA * pmax(0,tempS - soilTemp)) / 0.01),
           soilTempA = pmax(soilTemp, (fourierA / (mass * cpSoilA) + soilTemp)),
           # Change in proportion water this step
           moistureA = ifelse(moisture>0,max(0,moisture-((fourierA/2256400)/mWater)),
                              moisture),
           
           #2ND CM________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*mass,
           # Energy removed by current water quantity
           drain = ifelse(soilTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiB = max(fourierA-drain,0),
           #Thermal values
           cpSoilB = cpSoil((soilTemp+273.15), texture, peat, moistureA),
           saturationB = satSoil(texture, moistureA),
           kSoilB = kSoil(texture, saturationB),
           # Fourier conduction
           fourierB = pmin(fourierA,(kSoilB * pmax(0,soilTempA - soilTemp)) / 0.01),
           soilTempB = pmax(soilTemp, (fourierB / (mass * cpSoilB) + soilTemp)),
           # Change in proportion water this step
           moistureB = ifelse(moisture>0,max(0,moisture-((fourierB/2256400)/mWater)),
                              moisture),
           # Heat draw-down from above: pmin and soilTempA - ... makes sense, but this inexplicably seems to be doing it
           soilTempA = pmax(soilTempA, soilTempA + (fourierB / (mass * cpSoilA))),
           
           #3RD CM________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*mass,
           # Energy removed by current water quantity
           drain = ifelse(soilTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiC = max(fourierB-drain,0),
           #Thermal values
           cpSoilC = cpSoil((soilTemp+273.15), texture, peat, moistureB),
           saturationC = satSoil(texture, moistureB),
           kSoilC = kSoil(texture, saturationC),
           # Fourier conduction
           fourierC = pmin(fourierB,(kSoilC * pmax(0,soilTempB - soilTemp)) / 0.01),
           soilTempC = pmax(soilTemp, (fourierC / (mass * cpSoilC) + soilTemp)),
           # Change in proportion water this step
           moistureC = ifelse(moisture>0,max(0,moisture-((fourierC/2256400)/mWater)),
                              moisture),
           # Heat draw-down from above
           soilTempB = pmax(soilTempB, soilTempB + (fourierC / (mass * cpSoilB))),
           
           #4TH CM________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*mass,
           # Energy removed by current water quantity
           drain = ifelse(soilTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiD = max(fourierC-drain,0),
           #Thermal values
           cpSoilD = cpSoil((soilTemp+273.15), texture, peat, moistureC),
           saturationD = satSoil(texture, moistureC),
           kSoilD = kSoil(texture, saturationD),
           # Fourier conduction
           fourierD = pmin(fourierC,(kSoilD * pmax(0,soilTempC - soilTemp)) / 0.01),
           soilTempD = pmax(soilTemp, (fourierD / (mass * cpSoilD) + soilTemp)),
           # Change in proportion water this step
           moistureD = ifelse(moisture>0,max(0,moisture-((fourierD/2256400)/mWater)),
                              moisture),
           # Heat draw-down from above
           soilTempC = pmax(soilTempC, soilTempC + (fourierD / (mass * cpSoilC))),
           
           #5TH CM________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*mass,
           # Energy removed by current water quantity
           drain = ifelse(soilTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiE = max(fourierD-drain,0),
           #Thermal values
           cpSoilE = cpSoil((soilTemp+273.15), texture, peat, moistureD),
           saturationE = satSoil(texture, moistureD),
           kSoilE = kSoil(texture, saturationE),
           # Fourier conduction
           fourierE = pmin(fourierD,(kSoilE * pmax(0,soilTempD - soilTemp)) / 0.01),
           soilTempE = pmax(soilTemp, (fourierE / (mass * cpSoilE) + soilTemp)),
           # Change in proportion water this step
           moistureE = ifelse(moisture>0,max(0,moisture-((fourierE/2256400)/mWater)),
                              moisture),
           # Heat draw-down from above
           soilTempD = pmax(soilTempD, soilTempD + (fourierE / (mass * cpSoilD))),
           
           #Outgoing 1st________________________________________________________
           fourierOA = pmin(0,(kSoilA * (tempS - soilTempA)) / 0.01),
           soilTempA = pmin(soilTempA, soilTempA + (fourierOA / (mass * cpSoilA))),
           qrO = pmin(0,0.86*0.0000000000567*((tempS+273.15)^4 - (soilTempA+273.15)^4)),
           qR = qr + qrO,
           Q = qc + qR,
           
           
           #Outgoing 2nd________________________________________________________
           fourierOB = pmin(0,(kSoilB * (soilTempA - soilTempB)) / 0.01),
           soilTempB = pmin(soilTempB, soilTempB + (fourierOB / (mass * cpSoilB))),
           #Outgoing 3rd________________________________________________________
           fourierOC = pmin(0,(kSoilC * (soilTempB - soilTempC)) / 0.01),
           soilTempC = pmin(soilTempC, soilTempC + (fourierOC / (mass * cpSoilC))),
           #Outgoing 4th________________________________________________________
           fourierOD = pmin(0,(kSoilD * (soilTempC - soilTempD)) / 0.01),
           soilTempD = pmin(soilTempD, soilTempD + (fourierOD / (mass * cpSoilD))),
           #Outgoing 5th________________________________________________________
           fourierOE = pmin(0,(kSoilE * (soilTempD - soilTempE)) / 0.01),
           soilTempE = pmin(soilTempE, soilTempE + (fourierOE / (mass * cpSoilE))))
  
  
  soilTempA <- quantile(Ca$soilTempA, percentile)
  moistureA <- quantile(Ca$moistureA, percentile)
  soilTempB <- quantile(Ca$soilTempB, percentile)
  moistureB <- quantile(Ca$moistureB, percentile)
  soilTempC <- quantile(Ca$soilTempC, percentile)
  moistureC <- quantile(Ca$moistureC, percentile)
  soilTempD <- quantile(Ca$soilTempD, percentile)
  moistureD <- quantile(Ca$moistureD, percentile)
  soilTempE <- quantile(Ca$soilTempE, percentile)
  moistureE <- quantile(Ca$moistureE, percentile)
  
  # Advance one second's travel
  Horiz = Horiz - ROS
  pbar <- txtProgressBar(max = TIME, style = 3)
  # Loop through each time step and collect outputs
  for(t in 2:TIME){
    Cb <-threat(Surf, Plant, Horiz, Height=0, var, Pressure, Altitude) %>%
      mutate(t = t,
             #Convective transfer
             Re = (Plume_velocity*Density)/viscosity,
             h = ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888),
             #Incoming heat from surface
             tempS = ifelse(t>Ta, tempAir, ifelse(Horiz <0, pmax(tempAir, surface), tempAir)),
             qc = h * (tempS - soilTempA),
             att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
             qr = 0.86*qr*att,
             Qi = pmax(0, qc)+qr,
             
             #1ST CM________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureA*mass,
             # Energy removed by current water quantity
             drain = ifelse(soilTempA>95,
                            ifelse(moistureA>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiA = max(Qi-drain,0),
             #Thermal values
             cpSoilA = cpSoil((soilTemp+273.15), texture, peat, moistureA),
             saturationA = satSoil(texture, moistureA),
             kSoilA = kSoil(texture, saturationA),
             # Fourier conduction
             fourierA = ifelse(Horiz>0, pmin(QiA,(kSoilA * pmax(0,tempS - soilTempA)) / 0.01),
                               (kSoilA * pmax(0,tempS - soilTempA)) / 0.01),
             soilTempA = pmax(soilTempA, (fourierA / (mass * cpSoilA) + soilTempA)),
             # Change in proportion water this step
             moistureA = ifelse(moistureA>0,max(0,moistureA-((fourierA/2256400)/mWater)),
                                moistureA),
             
             #2ND CM________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureB*mass,
             # Energy removed by current water quantity
             drain = ifelse(soilTempB>95,
                            ifelse(moistureB>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiB = max(fourierA-drain,0),
             #Thermal values
             cpSoilB = cpSoil((soilTemp+273.15), texture, peat, moistureB),
             saturationB = satSoil(texture, moistureB),
             kSoilB = kSoil(texture, saturationB),
             # Fourier conduction
             fourierB = pmin(fourierA,(kSoilB * pmax(0,soilTempA - soilTempB)) / 0.01),
             soilTempB = pmax(soilTempB, (fourierB / (mass * cpSoilB) + soilTempB)),
             # Change in proportion water this step
             moistureB = ifelse(moistureB>0,max(0,moistureB-((fourierB/2256400)/mWater)),
                                moistureB),
             # Heat draw-down from above
             soilTempA = pmax(soilTempA, soilTempA + (fourierB / (mass * cpSoilA))),
             
             #3RD CM________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureC*mass,
             # Energy removed by current water quantity
             drain = ifelse(soilTempC>95,
                            ifelse(moistureC>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiC = max(fourierB-drain,0),
             #Thermal values
             cpSoilC = cpSoil((soilTemp+273.15), texture, peat, moistureC),
             saturationC = satSoil(texture, moistureC),
             kSoilC = kSoil(texture, saturationC),
             # Fourier conduction
             fourierC = pmin(fourierB,(kSoilC * pmax(0,soilTempB - soilTempC)) / 0.01),
             soilTempC = pmax(soilTempC, (fourierC / (mass * cpSoilC) + soilTempC)),
             # Change in proportion water this step
             moistureC = ifelse(moistureC>0,max(0,moistureC-((fourierC/2256400)/mWater)),
                                moistureC),
             # Heat draw-down from above
             soilTempB = pmax(soilTempB, soilTempB + (fourierC / (mass * cpSoilB))),
             
             #4TH CM________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureD*mass,
             # Energy removed by current water quantity
             drain = ifelse(soilTempD>95,
                            ifelse(moistureD>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiD = max(fourierC-drain,0),
             #Thermal values
             cpSoilD = cpSoil((soilTemp+273.15), texture, peat, moistureD),
             saturationD = satSoil(texture, moistureD),
             kSoilD = kSoil(texture, saturationD),
             # Fourier conduction
             fourierD = pmin(fourierC,(kSoilD * pmax(0,soilTempC - soilTempD)) / 0.01),
             soilTempD = pmax(soilTempD, (fourierD / (mass * cpSoilD) + soilTempD)),
             # Change in proportion water this step
             moistureD = ifelse(moistureD>0,max(0,moistureD-((fourierD/2256400)/mWater)),
                                moistureD),
             # Heat draw-down from above
             soilTempC = pmax(soilTempC, soilTempC + (fourierD / (mass * cpSoilC))),
             
             #5TH CM________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureD*mass,
             # Energy removed by current water quantity
             drain = ifelse(soilTempE>95,
                            ifelse(moistureE>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiE = max(fourierD-drain,0),
             #Thermal values
             cpSoilE = cpSoil((soilTemp+273.15), texture, peat, moistureE),
             saturationE = satSoil(texture, moistureE),
             kSoilE = kSoil(texture, saturationE),
             # Fourier conduction
             fourierE = pmin(fourierD,(kSoilE * pmax(0,soilTempD - soilTempE)) / 0.01),
             soilTempE = pmax(soilTempE, (fourierE / (mass * cpSoilE) + soilTempE)),
             # Change in proportion water this step
             moistureE = ifelse(moistureE>0,max(0,moistureE-((fourierE/2256400)/mWater)),
                                moistureE),
             # Heat draw-down from above
             soilTempD = pmax(soilTempD, soilTempD + (fourierE / (mass * cpSoilD))),
             
             #Outgoing 1st________________________________________________________
             fourierOA = pmin(0,(kSoilA * (tempS - soilTempA)) / 0.01),
             soilTempA = pmin(soilTempA, soilTempA + (fourierOA / (mass * cpSoilA))),
             qrO = pmin(0,0.86*0.0000000000567*((tempS+273.15)^4 - (soilTempA+273.15)^4)),
             qR = qr + qrO,
             Q = qc + qR,
             
             
             #Outgoing 2nd________________________________________________________
             fourierOB = pmin(0,(kSoilB * (soilTempA - soilTempB)) / 0.01),
             soilTempB = pmin(soilTempB, soilTempB + (fourierOB / (mass * cpSoilB))),
             #Outgoing 3rd________________________________________________________
             fourierOC = pmin(0,(kSoilC * (soilTempB - soilTempC)) / 0.01),
             soilTempC = pmin(soilTempC, soilTempC + (fourierOC / (mass * cpSoilC))),
             #Outgoing 4th________________________________________________________
             fourierOD = pmin(0,(kSoilD * (soilTempC - soilTempD)) / 0.01),
             soilTempD = pmin(soilTempD, soilTempD + (fourierOD / (mass * cpSoilD))),
             #Outgoing 5th________________________________________________________
             fourierOE = pmin(0,(kSoilE * (soilTempD - soilTempE)) / 0.01),
             soilTempE = pmin(soilTempE, soilTempE + (fourierOE / (mass * cpSoilE))))
    
    Ca <- rbind(Ca, Cb)
    
    soilTempA <- quantile(Cb$soilTempA, percentile)
    moistureA <- quantile(Cb$moistureA, percentile)
    soilTempB <- quantile(Cb$soilTempB, percentile)
    moistureB <- quantile(Cb$moistureB, percentile)
    soilTempC <- quantile(Cb$soilTempC, percentile)
    moistureC <- quantile(Cb$moistureC, percentile)
    soilTempD <- quantile(Cb$soilTempD, percentile)
    moistureD <- quantile(Cb$moistureD, percentile)
    soilTempE <- quantile(Cb$soilTempE, percentile)
    moistureE <- quantile(Cb$moistureE, percentile)
    
    setTxtProgressBar(pbar,t)
    ##  progress bar
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ", TIME - t)
      updateProgress(detail = text)
    }
    t = t + 1
    Horiz = Horiz - ROS
  }
  
  # Create table
  Ca <- Ca %>%
    select(t, repId, tempS, soilTempA, soilTempB, soilTempC, soilTempD, soilTempE,
           moistureA, moistureB, moistureC, moistureD, moistureE, fourierOA,
           fourierA, fourierB, fourierC, fourierD, fourierE)%>%
    mutate(startM = startM,
           seedDa = ifelse(soilTempA>100, 1, 0),
           seedDb = ifelse(soilTempB>100, 1, 0),
           seedDc = ifelse(soilTempC>100, 1, 0),
           seedDd = ifelse(soilTempD>100, 1, 0),
           seedDe = ifelse(soilTempE>100, 1, 0),
           seedga = ifelse(soilTempA>60, 1, 0)-seedDa,
           seedgb = ifelse(soilTempB>60, 1, 0)-seedDb,
           seedgc = ifelse(soilTempC>60, 1, 0)-seedDc,
           seedgd = ifelse(soilTempD>60, 1, 0)-seedDd,
           seedge = ifelse(soilTempE>60, 1, 0)-seedDe,
           orgA = pmax(0,pmin(1,0.0038*soilTempA-0.7692)),
           orgB = pmax(0,pmin(1,0.0038*soilTempB-0.7692)),
           orgC = pmax(0,pmin(1,0.0038*soilTempC-0.7692)),
           orgD = pmax(0,pmin(1,0.0038*soilTempD-0.7692)),
           orgE = pmax(0,pmin(1,0.0038*soilTempE-0.7692)),
           repelA = ifelse(soilTempA>175, 1, 0),
           repelB = ifelse(soilTempB>175, 1, 0),
           repelC = ifelse(soilTempC>175, 1, 0),
           repelD = ifelse(soilTempD>175, 1, 0),
           repelE = ifelse(soilTempE>175, 1, 0)
    )
  
  
  return(Ca)
}

#####################################################################

# Soil density
#
# Finds soil density from texture
# 
# Porosity of soils taken from Rawls, W. J., Brakensiek, D. L. & Saxton, K. E.
# Estimation of soil water properties. Transactions of the ASAE 25, 1316–1320 & 1328 (1982).
# 
# Density equation is from Peters-Lidard, C. D., Blackburn, E., Liang, X. & Wood, E. F.
# The effect of soil thermal conductivity parameterization on surface energy fluxes and temperatures.
# J. Atmos. Sci. 55, 1209–1224 (1998).


denSoil <- function(texture="loam")
{
  #Dry density
  porosity <- ifelse(texture=="sand",0.437,
                     ifelse(texture=="loamy sand",0.437,
                            ifelse(texture=="sandy loam", 0.453,
                                   ifelse(texture=="sandy clay loam",0.398,
                                          ifelse(texture=="sand clay",0.43,
                                                 ifelse(texture=="loam",0.463,
                                                        ifelse(texture=="clay loam",0.464,
                                                               ifelse(texture=="silt loam",0.501,
                                                                      ifelse(texture=="clay",0.475,
                                                                             ifelse(texture=="silty clay", 0.479,
                                                                                    ifelse(texture=="silty clay loam",0.471,0.5)))))))))))
  return((1-porosity)*2700)
}

#####################################################################

# Thermal conductivity of soil
#
# Model drawn from Johansen, O.
# Thermal conductivity of soils. PhD Thesis (University of Trondheim, 1971)
# Modified by Farouki, O. Thermal properties of soils. (Trans Tech, 1986)
# 
# Porosity taken from Rawls, W. J., Brakensiek, D. L. & Saxton, K. E.
# Estimation of soil water properties. Transactions of the ASAE 25, 1316–1320 & 1328 (1982).

kSoil <- function(texture="loam", saturation=0.3, grain="fine", unfrozen=1)
{
  saturation <-max(saturation, 0.1)
  porosity <- ifelse(texture=="sand",0.437,
                     ifelse(texture=="loamy sand",0.437,
                            ifelse(texture=="sandy loam", 0.453,
                                   ifelse(texture=="sandy clay loam",0.398,
                                          ifelse(texture=="sand clay",0.43,
                                                 ifelse(texture=="loam",0.463,
                                                        ifelse(texture=="clay loam",0.464,
                                                               ifelse(texture=="silt loam",0.501,
                                                                      ifelse(texture=="clay",0.475,
                                                                             ifelse(texture=="silty clay", 0.479,
                                                                                    ifelse(texture=="silty clay loam",0.471,0.5)))))))))))
  quartz <- ifelse(texture=="sand",0.92,
                   ifelse(texture=="loamy sand",0.82,
                          ifelse(texture=="sandy loam", 0.6,
                                 ifelse(texture=="sandy clay loam",0.6,
                                        ifelse(texture=="sand clay",0.52,
                                               ifelse(texture=="loam",0.4,
                                                      ifelse(texture=="clay loam",0.35,
                                                             ifelse(texture=="silt loam",0.25,
                                                                    ifelse(texture=="clay",0.25,
                                                                           ifelse(texture=="silty clay", 0.1,
                                                                                  ifelse(texture=="silty clay loam",0.1,
                                                                                         ifelse(texture=="silt",0.1,0))))))))))))
  minerals <- ifelse(texture=="sand",2,
                     ifelse(texture=="loamy sand",2,
                            ifelse(texture=="sandy loam", 2,
                                   ifelse(texture=="sandy clay loam",2,
                                          ifelse(texture=="sand clay",2,
                                                 ifelse(texture=="loam",2,
                                                        ifelse(texture=="clay loam",2,
                                                               ifelse(texture=="silt loam",2,
                                                                      ifelse(texture=="clay",2,3)))))))))
  
  kersten <- ifelse(grain == "fine",log10(saturation)+1,0.7*log10(saturation)+1)
  
  #Dry density (kg/m3)
  densityD <- (1-porosity)*2700
  # Dry thermal conductivity
  kDry <- (0.137*densityD+64.7)/(2700-0.947*densityD)
  # Solids thermal conductivity
  kS <- 7.7^quartz*minerals^(1-quartz)
  # Saturated thermal conductivity
  kSat <- kS^(1-porosity)*2.2^(porosity-unfrozen)*0.57^unfrozen
  
  return(kersten*(kSat-kDry)+kDry)
}


#####################################################################

# Specific heat of soil
#
# Finds volumetric specific heat from the mineral, organic and water components of the soil
#
# Specific heats of soil components estimated from figures 108 & 111 in
# Farouki, O. Thermal properties of soils. (Trans Tech, 1981).
#
# Water specific heat 4185 J/kg.K


cpSoil <- function(temp = 300, texture="loam", peat = 0.2, moisture=0.3)
{
  cpClay <- 3.66*temp-188
  cpSand <- 2.305*temp+67
  cpSilt <- 5.58*temp-854
  cpPeat <- 5.025*temp-151
  
  clay <- ifelse(texture=="sand",0.05,
                 ifelse(texture=="loamy sand",0.05,
                        ifelse(texture=="sandy loam", 0.15,
                               ifelse(texture=="sandy clay loam",0.25,
                                      ifelse(texture=="sand clay",0.4,
                                             ifelse(texture=="loam",0.2,
                                                    ifelse(texture=="clay loam",0.35,
                                                           ifelse(texture=="silt loam",0.15,
                                                                  ifelse(texture=="clay",0.7,
                                                                         ifelse(texture=="silty clay", 0.5,
                                                                                ifelse(texture=="silty clay loam",0.35,
                                                                                       ifelse(texture=="silt",0.05,0))))))))))))
  sand <- ifelse(texture=="sand",0.9,
                 ifelse(texture=="loamy sand",0.85,
                        ifelse(texture=="sandy loam", 0.65,
                               ifelse(texture=="sandy clay loam",0.6,
                                      ifelse(texture=="sand clay",0.5,
                                             ifelse(texture=="loam",0.4,
                                                    ifelse(texture=="clay loam",0.3,
                                                           ifelse(texture=="silt loam",0.2,
                                                                  ifelse(texture=="clay",0.15,
                                                                         ifelse(texture=="silty clay", 0.05,
                                                                                ifelse(texture=="silty clay loam",0.1,
                                                                                       ifelse(texture=="silt",0.05,0))))))))))))
  silt <- ifelse(texture=="sand",0.05,
                 ifelse(texture=="loamy sand",0.1,
                        ifelse(texture=="sandy loam", 0.2,
                               ifelse(texture=="sandy clay loam",0.15,
                                      ifelse(texture=="sand clay",0.1,
                                             ifelse(texture=="loam",0.4,
                                                    ifelse(texture=="clay loam",0.35,
                                                           ifelse(texture=="silt loam",0.65,
                                                                  ifelse(texture=="clay",0.15,
                                                                         ifelse(texture=="silty clay", 0.45,
                                                                                ifelse(texture=="silty clay loam",0.55,
                                                                                       ifelse(texture=="silt",0.9,0))))))))))))
  
  cpMineral <- clay*cpClay + sand*cpSand + silt*cpSilt
  cpDry <- peat*cpPeat+(1-peat)*cpMineral
  
  return(moisture*4185+(1-moisture)*cpDry)
}

#####################################################################

# Soil saturation
#
# Finds saturation from ODW moisture and texture
#
# Field capacity of soils taken from Salter, P. J. & Williams, J. B.
# The influence of texture on the moisture characteristics of soil.
# V. Relationships between particle-size composition and moisturecontents
# at the upper and lower limits of available-water. J. Soil Sci. 20, 126–131 (1969).

satSoil <- function(texture="loam", moisture=0.3)
{
  #Field capacity
  fcWW <- ifelse(texture=="sand",0.14,
                 ifelse(texture=="loamy sand",0.18,
                        ifelse(texture=="sandy loam", 0.26,
                               ifelse(texture=="sandy clay loam",0.26,
                                      ifelse(texture=="sand clay",0.29,
                                             ifelse(texture=="loam",0.30,
                                                    ifelse(texture=="clay loam",0.34,
                                                           ifelse(texture=="silt loam",0.39,
                                                                  ifelse(texture=="clay",0.42,
                                                                         ifelse(texture=="silty clay", 0.47,
                                                                                ifelse(texture=="silty clay loam",0.43,
                                                                                       ifelse(texture=="silt",0.45,0.45))))))))))))
  fc <- 1/((1/fcWW)-1)
  
  return(min(1,moisture/fc))
}


