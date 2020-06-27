#' Scorch height using an isotherm
#'
#' Calculates the height to which vegetation will be consumed,
#' and the height to a designated temperature isotherm reached for one second.
#'
#' Output fields are:
#' ht - scorch height from the surface flame (m)
#' hp - scorch height from burning plant (m)
#' Height - overall scorch height (m)
#' ns, e, m, c - height of consumption in each stratum (m)
#' b1, b2, b3, b4 - percentage of strata 1 (ns) to strata 4 (c) consumed
#' sc1, sc2, sc3, sc4 - percentage of strata 1 (ns) to strata 4 (c) scorched
#' d1, d2, d3, d4 - death (1) or survival (0) of standing foliage per stratum (50% or more scorch)
#'
#' @param Surf The dataframe 'runs' exported from Monte Carlos as 'Summary.csv'
#' @param Plant The dataframe 'IP' exported from Monte Carlos as 'IP.csv'.
#' @param Param A parameter dataframe used for FRaME,
#' such as produced using readLegacyParamFile
#' @param Test The temperature of the flora, default 60 degC
#' @return dataframe
#' @export

flora <- function(Surf, Plant, Param = Param, Test = 70)
{
  S <- strata(Param)
  n <- max(S$stratum)
  
  #Max burn height
  c <- Plant[Plant$level == "Canopy", ]%>%
    group_by(repId) %>%
    summarize_if(is.numeric,max)%>%
    mutate(c = y1)%>%
    select(repId, c)%>%
    right_join(Surf)
  c[is.na(c)] <- 0
  m <- Plant[Plant$level == "MidStorey", ]%>%
    group_by(repId) %>%
    summarize_if(is.numeric,max)%>%
    mutate(m = y1)%>%
    select(repId, m)%>%
    right_join(c)%>%
    mutate(m = ifelse(c>0, S$top[3], m))
  m[is.na(m)] <- 0
  e <- Plant[Plant$level == "Elevated", ]%>%
    group_by(repId) %>%
    summarize_if(is.numeric,max)%>%
    mutate(e = y1)%>%
    select(repId, e)%>%
    right_join(m)%>%
    mutate(e = ifelse(m>0, S$top[2], e))
  e[is.na(e)] <- 0
  ns <- Plant[Plant$level == "NearSurface", ]%>%
    group_by(repId) %>%
    summarize_if(is.numeric,max)%>%
    mutate(ns = y1)%>%
    select(repId, ns)%>%
    right_join(e)%>%
    mutate(ns = ifelse(e>0, S$top[1], ns))
  ns[is.na(ns)] <- 0
  
  #Heating from surface fire
  ground <- ns%>%
    mutate(Alpha = 1/(2*lengthSurface^2),
           C = 950*lengthSurface*exp(-Alpha*lengthSurface^2),
           pAlpha = abs(C/(Test-temperature)),
           R = pAlpha*cos(angleSurface),
           El = R * tan((slope_degrees * pi)/180),
           ht = (sin(angleSurface)*pAlpha - El) * extinct,
           ns = ns * extinct,
           e = e * extinct,
           m = m * extinct,
           c = c * extinct
    )%>%
    select(repId, extinct, temperature, slope_degrees, wind_kph, ht, ns, e, m, c)
  
  #Heating from plants
  plant<-Plant%>%
    left_join(ground)%>%
    mutate(Angle = atan((y1-y0)/(x1-x0)),
           Alpha = 1/(2*flameLength*(flameLength-length)),
           C = 950*flameLength*exp(-Alpha*(flameLength-length)^2),
           pAlpha = abs(C/(Test-temperature)),
           Reach = pAlpha*cos(Angle),
           El = Reach * tan((slope_degrees * pi)/180),
           htP = (sin(Angle)*pAlpha + y0 - El) * extinct)%>%
    select(repId, htP)%>%
    group_by(repId) %>%
    summarize_all(max)%>%
    right_join(ground)
  plant[is.na(plant)] <- 0
  
  #Find the corresponding number of each stratum
  NS <- S[S$name == "near surface", ]
  nNS <-ifelse(count(NS)==0, 0,as.numeric(as.character(NS$stratum)))
  E <- S[S$name == "elevated", ]
  nE <- ifelse(count(E)==0, 0,as.numeric(as.character(E$stratum)))
  M <- S[S$name == "midstorey", ]
  nM <-ifelse(count(M)==0, 0,as.numeric(as.character(M$stratum)))
  C <- S[S$name == "canopy", ]
  nC <- ifelse(count(C)==0, 0,as.numeric(as.character(C$stratum)))
  
  
  #Final heating
  Iso <- plant%>%
    mutate(Height = pmax(ht, htP),
           nsBase = ifelse(count(NS)==0, 0, S$base[which(S$name == "near surface")]),
           eBase = ifelse(count(E)==0, 0, S$base[which(S$name == "elevated")]),
           mBase = ifelse(count(M)==0, 0, S$base[which(S$name == "midstorey")]),
           cBase = ifelse(count(C)==0, 0, S$base[which(S$name == "canopy")]),
           nsTop = ifelse(count(NS)==0, 0, S$top[which(S$name == "near surface")]),
           eTop = ifelse(count(E)==0, 0, S$top[which(S$name == "elevated")]),
           mTop = ifelse(count(M)==0, 0, S$top[which(S$name == "midstorey")]),
           cTop = ifelse(count(C)==0, 0, S$top[which(S$name == "canopy")]),
           
           #Calculate scorch and burn
           b1 = ifelse(nsTop == 0, 0, round(pmin(100,pmax(0,100*(ns-nsBase)/(nsTop-nsBase))),0)),
           b2 = ifelse(eTop == 0, 0, round(pmin(100,pmax(0,100*(e-eBase)/(eTop-eBase))),0)),
           b3 = ifelse(mTop == 0, 0, round(pmin(100,pmax(0,100*(m-mBase)/(mTop-mBase))),0)),
           b4 = ifelse(cTop == 0, 0, round(pmin(100,pmax(0,100*(c-cBase)/(cTop-cBase))),0)),
           sc1 = ifelse(nsTop == 0, 0, round(pmin(100,pmax(0,100*(Height-nsBase)/(nsTop-nsBase))),0)),
           sc2 = ifelse(eTop == 0, 0, round(pmin(100,pmax(0,100*(Height-eBase)/(eTop-eBase))),0)),
           sc3 = ifelse(mTop == 0, 0, round(pmin(100,pmax(0,100*(Height-mBase)/(mTop-mBase))),0)),
           sc4 = ifelse(cTop == 0, 0, round(pmin(100,pmax(0,100*(Height-cBase)/(cTop-cBase))),0)))%>%
    select(repId, wind_kph, Height, ns, e, m, c, b1, b2, b3, b4, sc1, sc2, sc3, sc4)
  
  return(Iso)
}

#####################################################################

#' Finds radial bole necrosis depth
#'
#' Depth to which the cambium of a tree recieves lethal heating
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
#' Heat is transferred into the bark and timber using Fourier's Law
#'
#' Thermal conductivity of bark is modelled as per Martin, R. E.
#' Thermal properties of bark. For. Prod. J. 13, 419–426 (1963)
#'
#' Specific heat of bark is modelled using Kain, G., Barbu, M. C., Hinterreiter, S., Richter, K. & Petutschnigg, A.
#' Using bark as a heat insulation material. BioResources 8, 3718–3731 (2013)
#'
#' Thermal conductivity of wood is modelled using an approach from Kollmann, F. F. P. & Cote, W. A.
#' Principles of wood science and technology I. Solid wood. (Springer-Verlag, 1968)
#'
#' Evaporates water at 100 degrees C
#'
#' Specific heat of wood is derived from an established empirical relationship in Volbehr, B.
#' Swelling of wood fiber. PhD Thesis. (University of Kiel, 1896)
#'
#' Continues heating of the bole in the wake of the front for the duration of the surface fire for a period determined using
#' Burrows, N. D. Flame residence times and rates of weight loss of eucalypt forest fuel particles.
#' Int. J. Wildl. Fire 10, 137–143 (2001). Flame lengths are decreased exponentially over this period
#'
#' Bark is assumed to ignite, and burn for an average of resBark seconds, with the default value of 45s used as a mean for
#' figure 6c in Penman, T. D., Cawson, J. G., Murphy, S. & Duff, T. J.
#' Messmate Stringybark: bark ignitability and burning sustainability in relation to fragment dimensions,
#' hazard and time since fire. Int. J. Wildl. Fire 26, 866–876 (2017).
#'
#' Heats an area of 0.01m2
#'
#'
#'
#' @param Surf The dataframe 'runs' exported from Monte Carlos as 'Summary.csv'
#' @param Plant The dataframe 'IP' exported from Monte Carlos as 'IP.csv'.
#' @param percentile defines which heating statistics are used for each second, from 0 (min) to 1 (max)
#' @param Height The height on the bole directly over ground (m)
#' @param woodDensity The density of wood in the tree or log housing the hollow (kg/m3)
#' @param barkDensity The density of bark in the tree or log housing the hollow (kg/m3)
#' @param comBark Temperature directly under the burning bark (C)
#' @param bark The thickness of bark on the thinnest side of the hollow (m)
#' @param resBark Flame residence in the tree bark (s)
#' @param cambThick Thickness of the tree cambium (m)
#' @param xylemThick Thickness of the active xylem (m)
#' @param RH The relative humidity (0-1)
#' @param moisture The proportion oven-dry weight of moisture in the wood
#' @param bMoisture The proportion oven-dry weight of moisture in the bark
#' @param distance The furthest horizontal distance between the flame origin and the point (m)
#' @param trail The number of seconds to continue modelling after all flames have extinguished
#' @param var The angle in degrees that the plume spreads above/below a central vector
#' @param diameter depth of the litter layer (mm)
#' @param Pressure Sea level atmospheric pressure (hPa)
#' @param Altitude Height above sea level (m)
#' @param startTemp The starting temperature of wood and bark (deg C)
#' @param necT Temperature of necrosis (deg C)
#' @param surfDecl adusts the rate at which surface flame length declines after the front
#' @param updateProgress Progress bar for use in the dashboard
#' @return dataframe
#' @export

cambium <- function(Surf, Plant, percentile = 0.95, Height = 0.1, woodDensity = 700, barkDensity = 500,
                    bark = 0.04, comBark = 700, resBark = 45, cambThick = 0.01, xylemThick = 0.01, RH = 0.2,
                    moisture = 1, bMoisture = 0.2, distance = 50, trail = 100, var = 10, diameter = 20, Pressure = 1013.25,
                    Altitude = 0, startTemp = 25, necT = 60, surfDecl = 10,updateProgress=NULL)
{
  
  # Post-front surface flame heating
  lengthSurface <- mean(Surf$lengthSurface)
  residence <- 0.871*diameter^1.875
  depth <- diameter/1000
  
  # Collect step distance, time, and total distance
  ROS <- mean(Surf$ros_kph)/3.6
  Ta <- round(distance/ROS+residence)
  Tb <- round(distance/ROS)
  TIME <- Ta + trail
  Horiz <- distance
  
  # Description of the protection
  massB <- 0.01 * bark * barkDensity
  massW <- 0.0001 * woodDensity
  R <- sqrt(0.01/pi)
  startM <- moisture
  barkTemp <- startTemp
  woodTemp <- startTemp
  
  
  #Starting values
  Ca <- threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude)%>%
    mutate(t = 1,
           #Convective transfer
           Re = (Plume_velocity*Density)/viscosity,
           h = 0.35 + 0.47*Re^(1/2)*0.837,
           #Incoming heat from surface
           pt = pmax(0, t-Tb),
           comBark = ifelse(pt <= resBark, comBark, 0),
           postS = bole(lengthSurface,residence, depth, h = Height,
                        surfDecl = 10, t = pt),
           tempS = ifelse(Horiz <=0, pmax(tempAir, postS, comBark), tempAir),
           qc = h * (tempS - barkTemp),
           att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
           qr = 0.86*qr*att,
           Qi = pmax(0, qc)+qr,
           
           #BARK________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = bMoisture*massB,
           # Energy removed by current water quantity
           drain = ifelse(barkTemp>95,
                          ifelse(bMoisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiA = max(Qi-drain,0),
           #Bark thermal values
           rhoM = (bMoisture+bMoisture^2)*barkDensity,
           cpBark = ((1105+4.85*barkTemp)*(1-bMoisture)+bMoisture*4185+1276*bMoisture)/1000,
           kBark = (2.104*barkDensity+5.544*rhoM+3.266*barkTemp-166.216)*10^-4,
           # Fourier conduction
           fourierA = ifelse(Horiz>0, pmin(QiA,(kBark * pmax(0,tempS - barkTemp)) / bark),
                             (kBark * pmax(0,tempS - barkTemp)) / bark),
           barkTemp = pmax(barkTemp, (0.001 * fourierA / (massB * cpBark) + barkTemp)),
           # Change in proportion water this step
           moistureA = ifelse(bMoisture>0,max(0,bMoisture-((fourierA/2256400)/mWater)),
                              bMoisture),
           
           #1ST CM________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*massW,
           # Energy removed by current water quantity
           drain = ifelse(woodTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiB = max(fourierA-drain,0),
           #Thermal values
           kAir = 0.00028683*(woodTemp+273.15)^0.7919,
           cpWoodB = 1.08+0.00408*(100*moisture)+0.00253*woodTemp+0.0000628*(100*moisture)*woodTemp,
           kWoodB = kWood(woodTemp, woodDensity, kAir),
           # Fourier conduction
           fourierB = pmin(QiB,(kWoodB * pmax(0,barkTemp - woodTemp)) / 0.01),
           woodTempB = pmax(woodTemp, (0.001 * fourierB / (massW * cpWoodB) + woodTemp)),
           # Change in proportion water this step
           moistureB = ifelse(moisture>0,max(0,moisture-((fourierB/2256400)/mWater)),
                              moisture),
           # Heat draw-down from above
           barkTemp = pmax(barkTemp, barkTemp + (0.001 * fourierB / (massB * cpBark))),
           
           #2ND CM________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*massW,
           # Energy removed by current water quantity
           drain = ifelse(woodTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiC = max(fourierB-drain,0),
           #Thermal values
           kAir = 0.00028683*(woodTemp+273.15)^0.7919,
           cpWoodC = 1.08+0.00408*(100*moisture)+0.00253*woodTemp+0.0000628*(100*moisture)*woodTemp,
           kWoodC = kWood(woodTemp, woodDensity, kAir),
           # Fourier conduction
           fourierC = pmin(QiC,(kWoodC * pmax(0,woodTempB - woodTemp)) / 0.01),
           woodTempC = pmax(woodTemp, (0.001 * fourierC / (massW * cpWoodC) + woodTemp)),
           # Change in proportion water this step
           moistureC = ifelse(moisture>0,max(0,moisture-((fourierC/2256400)/mWater)),
                              moisture),
           # Heat draw-down from above
           woodTempB = pmax(woodTempB, woodTempB + (0.001 * fourierC / (massW * cpWoodB))),
           
           #3RD CM________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*massW,
           # Energy removed by current water quantity
           drain = ifelse(woodTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiD = max(fourierC-drain,0),
           #Thermal values
           kAir = 0.00028683*(woodTemp+273.15)^0.7919,
           cpWoodD = 1.08+0.00408*(100*moisture)+0.00253*woodTemp+0.0000628*(100*moisture)*woodTemp,
           kWoodD = kWood(woodTemp, woodDensity, kAir),
           # Fourier conduction
           fourierD = pmin(QiD,(kWoodD * pmax(0,woodTempC - woodTemp)) / 0.01),
           woodTempD = pmax(woodTemp, (0.001 * fourierD / (massW * cpWoodD) + woodTemp)),
           # Change in proportion water this step
           moistureD = ifelse(moisture>0,max(0,moisture-((fourierD/2256400)/mWater)),
                              moisture),
           # Heat draw-down from above
           woodTempC = pmax(woodTempC, woodTempC + (0.001 * fourierD / (massW * cpWoodC))),
           
           #4TH CM________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*massW,
           # Energy removed by current water quantity
           drain = ifelse(woodTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiE = max(fourierD-drain,0),
           #Thermal values
           kAir = 0.00028683*(woodTemp+273.15)^0.7919,
           cpWoodE = 1.08+0.00408*(100*moisture)+0.00253*woodTemp+0.0000628*(100*moisture)*woodTemp,
           kWoodE = kWood(woodTemp, woodDensity, kAir),
           # Fourier conduction
           fourierE = pmin(QiE,(kWoodE * pmax(0,woodTempD - woodTemp)) / 0.01),
           woodTempE = pmax(woodTemp, (0.001 * fourierE / (massW * cpWoodE) + woodTemp)),
           # Change in proportion water this step
           moistureE = ifelse(moisture>0,max(0,moisture-((fourierE/2256400)/mWater)),
                              moisture),
           # Heat draw-down from above
           woodTempD = pmax(woodTempD, woodTempD + (0.001 * fourierE / (massW * cpWoodD))),
           
           #Outgoing BARK________________________________________________________
           fourierOA = pmin(0,(kBark * (tempS - barkTemp)) / bark),
           barkTemp = pmin(barkTemp, barkTemp + (0.001 * fourierOA / (massB * cpBark))),
           qrO = pmin(0,0.86*0.0000000000567*((tempS+273.15)^4 - (barkTemp+273.15)^4)),
           qR = qr + qrO,
           Q = qc + qR,
           
           
           #Outgoing 1ST________________________________________________________
           fourierOB = pmin(0,(kWoodB * (barkTemp - woodTempB)) / 0.01),
           woodTempB = pmin(woodTempB, woodTempB + (0.001 * fourierOB / (massW * cpWoodB))),
           #Outgoing 2ND________________________________________________________
           fourierOC = pmin(0,(kWoodC * (woodTempB - woodTempC)) / 0.01),
           woodTempC = pmin(woodTempC, woodTempC + (0.001 * fourierOC / (massW * cpWoodC))),
           #Outgoing 3RD________________________________________________________
           fourierOD = pmin(0,(kWoodD * (woodTempC - woodTempD)) / 0.01),
           woodTempD = pmin(woodTempD, woodTempD + (0.001 * fourierOD / (massW * cpWoodD))),
           #Outgoing 4TH________________________________________________________
           fourierOE = pmin(0,(kWoodE * (woodTempD - woodTempE)) / 0.01),
           woodTempE = pmin(woodTempE, woodTempE + (0.001 * fourierOE / (massW * cpWoodE))))
  
  
  barkTemp <- quantile(Ca$barkTemp, percentile)
  moistureA <- quantile(Ca$moistureA, percentile)
  woodTempB <- quantile(Ca$woodTempB, percentile)
  moistureB <- quantile(Ca$moistureB, percentile)
  woodTempC <- quantile(Ca$woodTempC, percentile)
  moistureC <- quantile(Ca$moistureC, percentile)
  woodTempD <- quantile(Ca$woodTempD, percentile)
  moistureD <- quantile(Ca$moistureD, percentile)
  woodTempE <- quantile(Ca$woodTempE, percentile)
  moistureE <- quantile(Ca$moistureE, percentile)
  
  # Advance one second's travel
  Horiz = Horiz - ROS
  pbar <-  txtProgressBar(max = TIME, style = 3)
  # Loop through each time step and collect outputs
  for(t in 2:TIME){
    Cb <-threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude) %>%
      mutate(t = t,
             #Convective transfer
             Re = (Plume_velocity*Density)/viscosity,
             h = 0.35 + 0.47*Re^(1/2)*0.837,
             #Incoming heat from surface
             pt = pmax(0, t-Tb),
             comBark = ifelse(pt <= resBark, comBark, 0),
             postS = bole(lengthSurface,residence, depth, h = Height,
                          surfDecl = 10, t = pt),
             tempS = ifelse(t>Ta, tempAir, ifelse(Horiz <=0, pmax(tempAir, postS, comBark), tempAir)),
             qc = h * (tempS - barkTemp),
             att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
             qr = 0.86*qr*att,
             Qi = pmax(0, qc)+qr,
             
             
             #BARK________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureA*massB,
             # Energy removed by current water quantity
             drain = ifelse(barkTemp>95,
                            ifelse(moistureA>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiA = max(Qi-drain,0),
             #Bark thermal values
             rhoM = (moistureA+moistureA^2)*barkDensity,
             cpBark = ((1105+4.85*barkTemp)*(1-moistureA)+moistureA*4185+1276*moistureA)/1000,
             kBark = (2.104*barkDensity+5.544*rhoM+3.266*barkTemp-166.216)*10^-4,
             # Fourier conduction
             fourierA = ifelse(Horiz>0, pmin(QiA,(kBark * pmax(0,tempS - barkTemp)) / bark),
                               (kBark * pmax(0,tempS - barkTemp)) / bark),
             barkTemp = pmax(barkTemp, (0.001 * fourierA / (massB * cpBark) + barkTemp)),
             # Change in proportion water this step
             moistureA = ifelse(moistureA>0,max(0,moistureA-((fourierA/2256400)/mWater)),
                                moistureA),
             
             
             #1ST CM________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureB*massW,
             # Energy removed by current water quantity
             drain = ifelse(woodTempB>95,
                            ifelse(moistureB>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiB = max(fourierA-drain,0),
             #Thermal values
             kAir = 0.00028683*(woodTempB+273.15)^0.7919,
             cpWoodB = 1.08+0.00408*(100*moistureB)+0.00253*woodTempB+0.0000628*(100*moistureB)*woodTempB,
             kWoodB = kWood(woodTempB, woodDensity, kAir),
             # Fourier conduction
             fourierB = pmin(QiB,(kWoodB * pmax(0,barkTemp - woodTempB)) / 0.01),
             woodTempB = pmax(woodTempB, (0.001 * fourierB / (massW * cpWoodB) + woodTempB)),
             # Change in proportion water this step
             moistureB = ifelse(moistureB>0,max(0,moistureB-((fourierB/2256400)/mWater)),
                                moistureB),
             # Heat draw-down from above
             barkTemp = pmin(barkTemp, barkTemp - (0.001 * fourierB / (massB * cpBark))),
             
             #2ND CM________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureC*massW,
             # Energy removed by current water quantity
             drain = ifelse(woodTempC>95,
                            ifelse(moistureC>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiC = max(fourierB-drain,0),
             #Thermal values
             kAir = 0.00028683*(woodTempC+273.15)^0.7919,
             cpWoodC = 1.08+0.00408*(100*moistureC)+0.00253*woodTempC+0.0000628*(100*moistureC)*woodTempC,
             kWoodC = kWood(woodTempC, woodDensity, kAir),
             # Fourier conduction
             fourierC = pmin(QiC,(kWoodC * pmax(0,woodTempB - woodTempC)) / 0.01),
             woodTempC = pmax(woodTempC, (0.001 * fourierC / (massW * cpWoodC) + woodTempC)),
             # Change in proportion water this step
             moistureC = ifelse(moistureC>0,max(0,moistureC-((fourierC/2256400)/mWater)),
                                moistureC),
             # Heat draw-down from above
             woodTempB = pmin(woodTempB, woodTempB - (0.001 * fourierC / (massW * cpWoodB))),
             
             #3RD CM________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureD*massW,
             # Energy removed by current water quantity
             drain = ifelse(woodTempD>95,
                            ifelse(moistureD>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiD = max(fourierC-drain,0),
             #Thermal values
             kAir = 0.00028683*(woodTempD+273.15)^0.7919,
             cpWoodD = 1.08+0.00408*(100*moistureD)+0.00253*woodTempD+0.0000628*(100*moistureD)*woodTempD,
             kWoodD = kWood(woodTempD, woodDensity, kAir),
             # Fourier conduction
             fourierD = pmin(QiD,(kWoodD * pmax(0,woodTempC - woodTempD)) / 0.01),
             woodTempD = pmax(woodTempD, (0.001 * fourierD / (massW * cpWoodD) + woodTempD)),
             # Change in proportion water this step
             moistureD = ifelse(moistureD>0,max(0,moistureD-((fourierD/2256400)/mWater)),
                                moistureD),
             # Heat draw-down from above
             woodTempC = pmin(woodTempC, woodTempC - (0.001 * fourierD / (massW * cpWoodC))),
             
             #4TH CM________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureE*massW,
             # Energy removed by current water quantity
             drain = ifelse(woodTempE>95,
                            ifelse(moistureE>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiE = max(fourierD-drain,0),
             #Thermal values
             kAir = 0.00028683*(woodTempE+273.15)^0.7919,
             cpWoodE = 1.08+0.00408*(100*moistureE)+0.00253*woodTempE+0.0000628*(100*moistureE)*woodTempE,
             kWoodE = kWood(woodTempE, woodDensity, kAir),
             # Fourier conduction
             fourierE = pmin(QiE,(kWoodE * pmax(0,woodTempD - woodTempE)) / 0.01),
             woodTempE = pmax(woodTempE, (0.001 * fourierE / (massW * cpWoodE) + woodTempE)),
             # Change in proportion water this step
             moistureE = ifelse(moistureE>0,max(0,moistureE-((fourierE/2256400)/mWater)),
                                moistureE),
             # Heat draw-down from above
             woodTempD = pmin(woodTempD, woodTempD - (0.001 * fourierE / (massW * cpWoodD))),
             
             #Outgoing 1st________________________________________________________
             fourierOA = pmin(0,(kBark * (tempS - barkTemp)) / bark),
             barkTemp = pmin(barkTemp, barkTemp + (0.001 * fourierOA / (massB * cpBark))),
             qrO = pmin(0,0.86*0.0000000000567*((tempS+273.15)^4 - (barkTemp+273.15)^4)),
             qR = qr + qrO,
             Q = qc + qR,
             
             
             #Outgoing 1ST________________________________________________________
             fourierOB = pmin(0,(kWoodB * (barkTemp - woodTempB)) / 0.01),
             woodTempB = pmin(woodTempB, woodTempB + (0.001 * fourierOB / (massW * cpWoodB))),
             #Outgoing 2ND________________________________________________________
             fourierOC = pmin(0,(kWoodC * (woodTempB - woodTempC)) / 0.01),
             woodTempC = pmin(woodTempC, woodTempC + (0.001 * fourierOC / (massW * cpWoodC))),
             #Outgoing 3RD________________________________________________________
             fourierOD = pmin(0,(kWoodD * (woodTempC - woodTempD)) / 0.01),
             woodTempD = pmin(woodTempD, woodTempD + (0.001 * fourierOD / (massW * cpWoodD))),
             #Outgoing 4TH________________________________________________________
             fourierOE = pmin(0,(kWoodE * (woodTempC - woodTempE)) / 0.01),
             woodTempE = pmin(woodTempE, woodTempD + (0.001 * fourierOE / (massW * cpWoodE))))
    
    Ca <- rbind(Ca, Cb)
    
    barkTemp <- quantile(Cb$barkTemp, percentile)
    moistureA <- quantile(Cb$moistureA, percentile)
    woodTempB <- quantile(Cb$woodTempB, percentile)
    moistureB <- quantile(Cb$moistureB, percentile)
    woodTempC <- quantile(Cb$woodTempC, percentile)
    moistureC <- quantile(Cb$moistureC, percentile)
    woodTempD <- quantile(Cb$woodTempD, percentile)
    moistureD <- quantile(Cb$moistureD, percentile)
    woodTempE <- quantile(Cb$woodTempE, percentile)
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
    select(t, repId, tempS, barkTemp, woodTempB, woodTempC, woodTempD, woodTempE,
           moistureA, moistureB, moistureC, moistureD, moistureE, fourierOA,
           fourierA, fourierB, fourierC, fourierD, fourierE)%>%
    mutate(necrosis = ifelse(woodTempE>=necT, 4,
                             ifelse(woodTempD>=necT, 3,
                                    ifelse(woodTempC>=necT, 2,
                                           ifelse(woodTempB>=necT, 1, 0)))),
           ringbark = ifelse(necrosis>=cambThick, 1, 0),
           girdle = ifelse(necrosis>=cambThick+xylemThick, 1, 0))
  
  return(Ca)
}


#####################################################################

# Air temperature above ambient at the tree bole behind the flame front
#
# Dynamic air temperature at bole height, declining flame length exponentially
#
# Air temperature is modelled from dynamic flame segments using
# Weber R.O., Gill A.M., Lyons P.R.A., Moore P.H.R., Bradstock R.A., Mercer G.N. (1995)
# Modelling wildland fire temperatures. CALMScience Supplement, 4, 23–26.
#
# pAlphas is set to bole height - depth of surface litter

bole <- function(lengthSurface = 2, residence = 300, depth = 0.05, h = 0.1, surfDecl = 10, t = 1)
{
  surfPost <- lengthSurface*exp(-(surfDecl/residence)*t)
  Alpha <- 1/(2 * surfPost^2)
  C <- 950 * surfPost * exp(-Alpha * surfPost^2)
  pAlphas <- pmax(0,h-depth)
  return(ifelse(pAlphas < surfPost,
                950 + exp(-Alpha * pAlphas^2),
                C/pAlphas))
}


#####################################################################
# Thermal conductivity of dry wood
#
# Model drawn from Kollmann, F. F. P. & Cote, W. A.
# Principles of wood science and technology I. Solid wood. (Springer-Verlag, 1968)


kWood <- function(T=100, rhoW=700, kAir = 0.026)
{
  T <- T+273.15
  bridge <- 0.0019*T+0.0503
  r <- 1-(rhoW/1500)
  kA <- (1-r)*0.766+r*kAir
  kB <- 1/((1-r)/0.43+(r/kAir))
  return(bridge*kA+(1-bridge)*kB)
}

