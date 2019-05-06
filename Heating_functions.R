#' Heat and plume produced by a fire front
#'
#' Finds the temperature, velocity, dynamic viscosity, atmospheric pressure, and density of a plume at a point
#' Gives the mean flame temperature (K) and emissive power of each flame
#'
#' tempAir: Air temperature is modelled from dynamic flame segments using
#' Weber R.O., Gill A.M., Lyons P.R.A., Moore P.H.R., Bradstock R.A., Mercer G.N. (1995)
#' Modelling wildland fire temperatures. CALMScience Supplement, 4, 23–26.
#'
#' cpAir: Specific heat of air is found from an empirical function fit to data from
#' Hilsenrath, J. et al.
#' Circular of the Bureau of Standards no. 564: tables of thermal properties of gases comprising tables of
#' thermodynamic and transport properties of air, argon, carbon dioxide, carbon monoxide hydrogen, nitrogen,
#' oxygen, and steam. (Department of Commerce, 1955), and
#' Kyle, B. G. Chemical and process thermodynamics. (Englewood Cliffs / Prentice Hall, 1984)
#'
#' Density: Air density (kg/m^3) is found using the Ideal Gas Law.
#'
#' presAtm: Air pressure (hPa) is calculated using the Barometric Formula with standard Values.
#'
#' viscosity: Dynamic viscosity (10^-6 Pa.s) is calculated using an empirical relationship fit to a subset of
#' air pressures within an order of magnitude of sea level (R2 = 0.999), using data from
#' Kadoya, K., N, M. & Nagashima, A. Viscosity and thermal conductivity of dry air in the gaseous phase.
#' J. Phys. Chem. Ref. Data 14, 947–970 (1985)
#'
#' Plume_velocity: Plume velocity (m/s) is calculated allowing for crosswinds, using
#' Oka Y., Sugawa O., Imamura T. (2008)
#' Correlation of temperature rise and velocity along an inclined fire plume axis in crosswinds.
#' Fire Safety Journal, 43, 391–400.
#'
#' The Froude number is set to 1.5 as per Ma T.G., Quintiere J.G. (2003)
#' Numerical simulation of axi-symmetric fire plumes: Accuracy and limitations.
#' Fire Safety Journal, 38, 467–492.
#'
#' flameTemp: Average flame temperature is integrated from the modelled temperature over the length of the flame
#'
#' epsilon: The emissivity of a flame is calculated from mean flame thickness using
#' Àgueda, A., Pastor, E., Pérez, Y. & Planas, E.
#' Experimental study of the emissivity of flames resulting from the combustion of forest fuels.
#' Int. J. Therm. Sci. 49, 543–554 (2010)
#'
#' phi: The configuration factor is calculated using the function phi, from
#' Tan, Z., Midgley, S. & Douglas, G.
#' A computerised model for bushfire attack assessment and its applications in bushfire protection planning.
#' Int. Congr. Model. Simul. Adv. Appl. Manag. Decis. Making, MODSIM05 538–545 (2005)
#'
#'
#' @param Surf The dataframe 'runs' exported from Monte Carlos as 'Summary.csv'
#' @param repFlame The dataframe 'repFlame'
#' @param Pressure Sea level atmospheric pressure (hPa)
#' @param Altitude Height above sea level (m)
#' @param Horizontal The horizontal distance in metres from the flame origin to the point
#' @param Height The vertical distance in metres from the flame origin to the point
#' @param var The angle in degrees that the plume spreads above/below a central vector;defaults to 10
#' @return dataframe
#' @export Iso.csv A .csv file listing the flora height in m from each model run
#' @examples threat(Surf, repFlame, Horizontal = 10, Height = 10, var = 10, Pressure = 1013.25, Altitude = 0)

threat <- function (Surf, repFlame, Horizontal = 10, Height = 10, var = 10, Pressure = 1013.25, Altitude = 0)
{
  # Calculate surface flame heating
  s <- Surf %>%
    mutate(hor = Horizontal,
           Alpha = 1/(2 * lengthSurface^2),
           C = 950 * lengthSurface * exp(-Alpha * lengthSurface^2),
           pAlphas = abs(Horizontal/cos(angleSurface)),
           El = Horizontal * tan((slope_degrees * pi)/180),
           Yl = Horizontal * (angleSurface - tan((var * pi)/180)),
           Yu = Horizontal * (angleSurface + tan((var * pi)/180)),
           IntL = ifelse(Height + El < Yl, 0, 1),
           IntU = ifelse(Height + El > Yu, 0, 1),
           Intercept = IntL * IntU,
           temp_pointS = ifelse(pAlphas < lengthSurface,
                 950 + exp(-Alpha * pAlphas^2),
                 C/pAlphas) * Intercept + temperature) %>%
    select(repId, hor, ros_kph, wind_kph, angleSurface, lengthSurface, pAlphas, temperature, slope_degrees, temp_pointS, epsilon)

  # Calculate heating from burning plants
  p <- repFlame %>%
    left_join(s) %>%
    mutate(Angle = atan((y1 - y0)/(x1 - x0)),
           Alpha = 1/(2 * flameLength * (flameLength - length)),
           C = 950 * flameLength * exp(-Alpha * (flameLength - length)^2),
           pAlphap = abs(Horizontal/cos(Angle)),
           El = Horizontal * tan((slope_degrees * pi)/180),
           Yl = Horizontal * (Angle - tan((var * pi)/180)) + y0,
           Yu = Horizontal * (Angle + tan((var * pi)/180) + y0),
           IntL = ifelse(Height + El < Yl, 0, 1),
           IntU = ifelse(Height + El > Yu, 0, 1),
           Intercept = IntL * IntU,
           temp_pointP = ifelse(pAlphap < flameLength, 950 + exp(-Alpha * (pAlphap - length)^2),
                                C/pAlphap) * Intercept + temperature,
           flameTemp = 1045*exp(0.155*length/flameLength)) %>%
    group_by(repId) %>%
    filter(temp_pointP == max(temp_pointP)) %>%
    select(repId, repHeight, repHeight, repLength, repAngle, runIndex, segIndex, x0, y0, x1, y1,
           length, flameLength, hor, ros_kph, wind_kph, angleSurface, lengthSurface, pAlphas,
           temperature, slope_degrees, temp_pointS, epsilon, Angle, Alpha, C, pAlphap, El, Yl,
           Yu, IntL, IntU, Intercept, temp_pointP, flameTemp)%>%
    group_by(repId) %>%
    summarize_all(max)%>%
    right_join(s)
  p[is.na(p)] <- 0

  # Calculate heat transfer inputs
  Con <- p %>%
    mutate(tempAir = pmax(temp_pointS, temp_pointP),
           cpAir = 0.00019*(tempAir+273.13)+0.9405,
           pAlpha=max(0.05, ifelse(tempAir == temp_pointS, pAlphas, pAlphap)),
           viscosity = 0.000000388 * (tempAir + 273.15)^0.6818,
           presAtm = (Pressure/10)*exp(-0.00012*Altitude),
           Density = (presAtm * 1000) / (287.05 * (tempAir + 273.15)),
           Plume_velocity = (((1.24*(pAlpha/(pAlpha+2*hor))) *
                                (sqrt(((tempAir-temperature)/(temperature+273.15))*9.81*pAlpha)))^2+
                               (((wind_kph/3.6)*0.47*((pAlpha+hor)/pAlpha)^0.25)^2))^0.5,
           flameTemp = ifelse(flameTemp==0, 1045, flameTemp),
           E = pmax(0, epsilon*0.0000000567*(flameTemp^4-(tempAir+273.15)^4)),
           repAngle = ifelse(repLength>lengthSurface, repAngle, angleSurface),
           repLength = max(repLength, lengthSurface),
           phi = phi(repLength,repAngle,slope_degrees,Horizontal,Height),
           qr = E * phi) %>%
    select(repId, ros_kph, wind_kph, temperature, pAlpha, tempAir, cpAir,
           viscosity, presAtm, Density, Plume_velocity, flameTemp, epsilon, E, phi, qr)
return(Con)
}


#####################################################################

#' Fire risk for an exposed arboreal animal
#'
#' Calculates the degree of injury or likelihood of mortality to an animal
#' caused by an approaching fire front
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
#' Mammal pelage is given a standardised emissivity of 0.86, based on:
#' McGowan, N. E., Scantlebury, D. M., Maule, A. G. & Marks, N. J.
#' Measuring the emissivity of mammal pelage. Quant. Infrared Thermogr. J. 6733, 1–9 (2018).
#'
#' Skin Cp is set to 3.5, averaged from Duck, F. A.
#' Physical properties of tissues: a comprehensive reference book. (Elsevier Science, 1990).
#'
#' Finds furDensity - the density of fur or feathers on the animal using standard values of 66 fibres per mm2,
#' fibre diameter of 0.01mm (10μm), and α-keratin density of 1300 kg.m-3. Fibre count averaged from
#' Liwanag, H. E. M., Berta, A., Costa, D. P., Abney, M. & Williams, T. M.
#' Morphological and thermal properties of mammalian insulation: the evolution of fur for aquatic living.
#' Biol. J. Linn. Soc. 106, 926–939 (2012)
#'
#' Finds thermal dose from all heat inputs, using the formula from
#' Ciesielski, M., Mochnacki, B. & Szopa, R. Numerical modeling of biological tissue heating.
#' Admissible thermal dose. Sci. Res. Inst. Math. Comput. Sci. 10, 11–20 (2011).
#'
#'
#' @param Surf The dataframe 'runs' exported from Monte Carlos as 'Summary.csv'
#' @param Plant The dataframe 'IP' exported from Monte Carlos as 'IP.csv'.
#' @param percentile defines which heating statistics are used for each second, from 0 (min) to 1 (max)
#' @param Height The height directly ove ground (m) at which the species is expected to shelter from a fire.
#' @param low The closest horizontal distance between the flame origin and the point (m)
#' @param high The furthest horizontal distance between the flame origin and the point (m)
#' @param var The angle in degrees that the plume spreads above/below a central vector
#' @param Pressure Sea level atmospheric pressure (hPa)
#' @param Altitude Height above sea level (m)
#' @param RH The relative humidity (0-1)
#' @param Class Class of animal. Allowable values are "mammalia", "aves", "amphibia" and "reptilia"
#' @param Dimension The "Characteristic length" of the animal (m)
#' @param Area The surface area of the animal (m^2)
#' @param protection The thickness of fur, feather or scales covering the animal (m)
#' @param count The number of fibres per square mm
#' @param fibre The mean fibre diameter of hairs in mm
#' @param Specific_heat The specific heat of the fur, feather or scale material only (kJ/kg/deg C)
#' @param skinCp The specific heat of the animal skin (kJ/kg/C)
#' @param skinK Thermal conductivity of the animal skin (W/m/C)
#' @param objectTemp The body temperature of the animal (deg C)
#' @param Shape The approximate shape of the animal - either "Flat", "Sphere", or "Cylinder"
#' @return dataframe
#' @export arboreal.csv
#' @examples fauna(Surf, Plant, percentile = 0.5, Height = 1, low = 1, high = 50, var = 10, Pressure = 1013.25,
#'           Altitude = 0, RH = 0.51, Class = "mammalia", Dimension = 0.1, Area = 0.2,
#'           protection = 0.02, count = 66, fibre = 0.01, Specific_heat = 2.5, skinCp = 3.5,
#'           skinK = 2, objectTemp = 38, Shape = "Flat")

arboreal <- function(Surf, Plant, percentile = 0.5, Height = 1, low = 1, high = 50, var = 10, Pressure = 1013.25,
                  Altitude = 0, RH = 0.51, Class = "mammalia", Dimension = 0.1, Area = 0.2,
                  protection = 0.02, count = 66, fibre = 0.01, Specific_heat = 2.5, skinCp = 3.5,
                  skinK = 2, objectTemp = 38, Shape = "Flat")
{
  # Collect testing stats
  ROS <- mean(Surf$ros_kph)/3.6
  TIME <- round((high - low)/ROS)
  Horiz <- high
  dens <- (count*pi*(fibre/2)^2)
  Volume <- Area * protection
  R <- sqrt(Area/pi)
  cl <- ifelse(Class == "amphibia", 1, 0)
  tempR <- objectTemp

  #Starting values
  Ca <- threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude)%>%
    mutate(t = 1,
           furDensity = (1300*dens)+(1-dens)*Density,
           furCp = (Specific_heat*dens)+(1-dens)*cpAir,
           specificHeat = (cl*skinCp)+abs(cl-1)*furCp,
           Mass = Volume * furDensity,
           Re = (Plume_velocity*Density*Dimension)/viscosity,
           hFlat = ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888),
           hSphere = 2 + 0.6*Re^(0.5)*0.888,
           hCylinder = 0.35 + 0.47*Re^(1/2)*0.837,
           h = ifelse(Shape == "Flat", hFlat, ifelse(Shape == "Sphere", hSphere, hCylinder)),
           #Incoming
           qc = h * Area *(tempAir - tempR),
           att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
           qr = 0.86*qr*att,
           Qi = pmax(0, qc)+qr,
           # Thermal conductivity
           kAir = 0.00028683*(tempAir+273.15)^0.7919,
           kFur = kAir+0.004853*protection,
           k = (cl*skinK)+abs(cl-1)*kFur,
           # Fourier conduction through pelage
           fourier = pmin(Qi,(Area * k * (tempAir - tempR)) / protection),
           tempFourier = pmax(tempR, (0.001 * fourier / (Mass * specificHeat) + tempR)),
           # Thermal dose
           thermalDose = 0.69 * fourier * R,
           fatalityLikelihood = pmax(0,pmin(100, 75.405 * log(thermalDose) - 518.12)),
           #Outgoing
           fourierO = pmin(0,(Area * k * (tempAir - tempFourier)) / protection),
           tempR = pmax(objectTemp, tempFourier + (0.001 * fourierO / (Mass * specificHeat))),
           qrO = pmin(0,Area*0.86*0.0000000000567*((tempAir+273.15)^4 - (tempR+273.15)^4)),
           qR = qr + qrO,
           Q = qc + qR,
           # Conservation of energy
           deltaT = 0.001*Q / (Mass * specificHeat),
           tempConservation = pmax(objectTemp, (tempR + deltaT)))

  # Average the object temperature. These become starting values for each step of the loop
#  SiteM <- mean(Ca$tempConservation)
#  SiteMCond <- mean(Ca$tempR)
#  TDM <- mean(Ca$thermalDose)
  SiteM <- quantile(Ca$tempConservation, percentile)
  SiteMCond <- quantile(Ca$tempR, percentile)
  TDM <- quantile(Ca$thermalDose, percentile)

  # Advance one second's travel
  Horiz = Horiz - ROS

  # Loop through each time step and collect outputs
  for(t in 2:TIME){
    Cb <-threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude) %>%
      mutate(t = t,
             furDensity = (1300*dens)+(1-dens)*Density,
             furCp = (Specific_heat*dens)+(1-dens)*cpAir,
             specificHeat = (cl*skinCp)+abs(cl-1)*furCp,
             Mass = Volume * furDensity,
             tempConservation = SiteM,
             tempR = SiteMCond,
             thermalDose = TDM,
             Re = (Plume_velocity*Density*Dimension)/viscosity,
             hFlat = ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888),
             hSphere = 2 + 0.6*Re^(0.5)*0.888,
             hCylinder = 0.35 + 0.47*Re^(1/2)*0.837,
             h = ifelse(Shape == "Flat", hFlat, ifelse(Shape == "Sphere", hSphere, hCylinder)),
             #Incoming
             qc = h * Area *(tempAir - tempR),
             att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
             qr = 0.86*qr*att,
             Qi = pmax(0, qc)+qr,
             # Thermal conductivity
             kAir = 0.00028683*(tempAir+273.15)^0.7919,
             kFur = kAir+0.004853*protection,
             k = (cl*skinK)+abs(cl-1)*kFur,
             # Fourier conduction through pelage
             fourier = pmin(Qi,(Area * k * (tempAir - tempR)) / protection),
             tempFourier = pmax(tempR, (0.001 * fourier / (Mass * specificHeat) + tempR)),
             # Thermal dose
             thermalDose = thermalDose + 0.69 * pmax(fourier, 0) * R,
             fatalityLikelihood = pmax(0,pmin(100, 75.405 * log(thermalDose) - 518.12)),
             #Outgoing
             fourierO = pmin(0,(Area * k * (tempAir - tempFourier)) / protection),
             tempR = pmax(objectTemp, tempFourier + (0.001 * fourierO / (Mass * specificHeat))),
             qrO = pmin(0,Area*0.86*0.0000000000567*((tempAir+273.15)^4 - (tempR+273.15)^4)),
             qR = qr + qrO,
             Q = qc + qR,
             # Conservation of energy
             deltaT = 0.001*Q / (Mass * specificHeat),
             tempConservation = pmax(objectTemp, (tempConservation + deltaT)))
    Ca <- rbind(Ca, Cb)

    # Average the object temperature
#    SiteM <- mean(Cb$tempConservation)
#    SiteMCond <- mean(Cb$tempR)
#    TDM <- mean(Cb$thermalDose)
    SiteM <- quantile(Cb$tempConservation, percentile)
    SiteMCond <- quantile(Cb$tempR, percentile)
    TDM <- quantile(Cb$thermalDose, percentile)

    t = t + 1
    Horiz = Horiz - ROS
  }

  # Create and export table
  Ca <- Ca %>%
    mutate(plume_kph = Plume_velocity * 3.6) %>%
    select(t, repId, ros_kph, wind_kph, plume_kph, presAtm, tempR,
           tempAir, viscosity, pAlpha, Re, h,
           Mass, furDensity, specificHeat,
           flameTemp, epsilon, E, phi, att,
           qc, qr, qrO, qR, Q,
           deltaT, tempConservation,
           kAir, kFur, k, fourier, tempFourier, tempR, fourierO,
           thermalDose, fatalityLikelihood)
  write.csv(Ca, "arboreal.csv")
  return(Ca)
}

#####################################################################

#' Fire risk for an animal sheltering in a wooden hollow
#'
#' Calculates the likelihood of mortality to an animal
#' caused by an approaching fire front
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
#' Finds animal mortality within a hollow based on the maximum tolerable temperature for a given
#' vapour pressure deficit, based on data from Lawrence, G. E.
#' Ecology of vertebrate animals in relation to chaparral fire in the Sierra Nevada foothills.
#' Ecology 47, 278–291 (1966)
#'
#' Heat is transferred into the hollow using Fourier's Law
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
#'
#' @param Surf The dataframe 'runs' exported from Monte Carlos as 'Summary.csv'
#' @param Plant The dataframe 'IP' exported from Monte Carlos as 'IP.csv'.
#' @param percentile defines which heating statistics are used for each second, from 0 (min) to 1 (max)
#' @param Height The height directly over ground (m) at which the species is expected to shelter from a fire.
#' @param woodDensity The density of wood in the tree or log housing the hollow (kg/m3)
#' @param barkDensity The density of bark in the tree or log housing the hollow (kg/m3)
#' @param wood The thickness of wood on the thinnest side of the hollow (m)
#' @param bark The thickness of bark on the thinnest side of the hollow (m)
#' @param RH The relative humidity (0-1)
#' @param water The proportion oven-dry weight of moisture in the bark and wood
#' @param low The closest horizontal distance between the flame origin and the point (m)
#' @param high The furthest horizontal distance between the flame origin and the point (m)
#' @param var The angle in degrees that the plume spreads above/below a central vector;defaults to 10
#' @param Pressure Sea level atmospheric pressure (hPa)
#' @param Altitude Height above sea level (m)
#' @param Dimension The "Characteristic length" of the hollow (m)
#' @param Area The surface area of the thinnest side of the hollow (m^2)
#' @param hollowTemp The starting temperature inside the hollow (deg C)
#' @param Shape The approximate shape of the hollow exterior - either "Flat", "Sphere", or "Cylinder"
#' @return dataframe
#' @export hollow.csv
#' @examples hollow(Surf, Plant, percentile = 0.5, Height = 1, woodDensity = 700, barkDensity = 500,
#'           RH = 0.2, water = 0.2, low = 1, high = 50, var = 10, Pressure = 1013.25, Altitude = 0,
#'           Dimension = 0.1, Area = 0.2, wood = 0.1, bark = 0.02, hollowTemp = 25, Shape = "Flat")

hollow <- function(Surf, Plant, percentile = 0.5, Height = 1, woodDensity = 700, barkDensity = 500,
                   wood = 0.1, bark = 0.02, RH = 0.2, water = 0.2, low = 1, high = 50, var = 10, Pressure = 1013.25,
                   Altitude = 0, Dimension = 0.3, Area = 0.03, hollowTemp = 25, Shape = "Flat")
{
  # Collect step distance, time, and total distance
  ROS <- mean(Surf$ros_kph)/3.6
  TIME <- round((high - low)/ROS)
  Horiz <- high
  # Description of the protection
  vWood <- Area * wood
  mWood <- vWood * woodDensity
  vBark <- Area * bark
  mBark <- vBark * barkDensity
  R <- sqrt(Area/pi)
  tempW <- hollowTemp
  tempB <- hollowTemp
  # Moisture effects on bark
  rhoM <- (water+water^2)*barkDensity

  #Starting values
  Ca <- threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude)%>%
    mutate(t = 1,
           #Convective transfer
           Re = (Plume_velocity*Density*Dimension)/viscosity,
           hFlat = ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888),
           hSphere = 2 + 0.6*Re^(0.5)*0.888,
           hCylinder = 0.35 + 0.47*Re^(1/2)*0.837,
           h = ifelse(Shape == "Flat", hFlat, ifelse(Shape == "Sphere", hSphere, hCylinder)),
           #Incoming heat
           qc = h * Area *(tempAir - tempB),
           att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
           qr = 0.86*qr*att,
           Qi = pmax(0, qc)+qr,
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWaterB = water*mBark,
           # Energy removed by current water quantity
           drain = ifelse(tempB>95,
                          ifelse(water>0,mWaterB*2256400,0),0),
           # Change in proportion bark water this step
           waterB = ifelse(tempB>95,
                           ifelse(water>0,max(0,water-((Qi/2256400)/mWaterB)),
                                  water),water),
           # Adjusted incoming energy and moisture effect on bark
           Qi = max(Qi-drain,0),
           rhoM <- (waterB+waterB^2)*barkDensity,
           #Bark thermal values
           cpBark = ((1105+4.85*tempAir)*(1-waterB)+waterB*4185+1276*waterB)/1000,
           kAir = 0.00028683*(tempAir+273.15)^0.7919,
           kBark = (2.104*barkDensity+5.544*rhoM+3.266*tempB-166.216)*10^-4,
           # Fourier conduction through bark
           fourierB = pmin(Qi,(Area * kBark * (tempAir - tempB)) / bark),
           tempBark = pmax(tempB, (0.001 * fourierB / (mBark * cpBark) + tempB)),
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWaterW = water*mWood,
           # Energy removed by current water quantity
           drain = ifelse(tempW>95,
                          ifelse(water>0,mWaterW*2256400,0),0),
           # Change in proportion wood water this step
           waterW = ifelse(tempW>95,
                           ifelse(water>0,max(0,water-((Qi/2256400)/mWaterW)),
                                  water),water),
           # Adjusted incoming energy
           Qi = max(Qi-drain,0),
           #Wood thermal values
           cpWood = 1.08+0.00408*(100*waterW)+0.00253*tempAir+0.0000628*(100*waterW)*tempAir,
           kWood = kWood(tempW, woodDensity, kAir),
           # Fourier conduction through wood
           fourierW = pmin(fourierB,(Area * kWood * (tempBark - tempW)) / wood),
           tempWood = pmax(tempW, (0.001 * fourierW / (mWood * cpWood) + tempW)),
           # Mortality
           mortality = ifelse(tempWood < 67.5-30.17*RH, 0, 1),
           #Outgoing wood
           fourierOW = pmin(0,(Area * kWood * (tempBark - tempWood)) / wood),
           tempW = pmax(tempW, tempWood + (0.001 * fourierOW / (mWood * cpWood))),
           #Outgoing bark
           fourierOB = pmin(0,(Area * kBark * (tempAir - tempBark)) / bark),
           tempB = pmax(tempB, tempBark + (0.001 * fourierOB / (mBark * cpBark))),
           qrO = pmin(0,Area*0.86*0.0000000000567*((tempAir+273.15)^4 - (tempB+273.15)^4)),
           qR = qr + qrO,
           Q = qc + qR)

  woodM <- quantile(Ca$tempW, percentile)
  barkM <- quantile(Ca$tempB, percentile)
  waterW <- quantile(Ca$waterW, percentile)
  waterB <- quantile(Ca$waterB, percentile)

  # Advance one second's travel
  Horiz = Horiz - ROS

  # Loop through each time step and collect outputs
  for(t in 2:TIME){
    Cb <-threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude) %>%
      mutate(t = t,
             #Convective transfer
             tempW = woodM,
             tempB = barkM,
             Re = (Plume_velocity*Density*Dimension)/viscosity,
             hFlat = ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888),
             hSphere = 2 + 0.6*Re^(0.5)*0.888,
             hCylinder = 0.35 + 0.47*Re^(1/2)*0.837,
             h = ifelse(Shape == "Flat", hFlat, ifelse(Shape == "Sphere", hSphere, hCylinder)),
             #Incoming
             qc = h * Area *(tempAir - tempB),
             att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
             qr = 0.86*qr*att,
             Qi = pmax(0, qc)+qr,
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWaterB = waterB*mBark,
             # Energy removed by current water quantity
             drain = ifelse(tempB>95,
                            ifelse(waterB>0,mWaterB*2256400,0),0),
             # Change in proportion water this step
             waterB = ifelse(tempB>95,
                             ifelse(waterB>0,max(0,waterB-((Qi/2256400)/mWaterB)),
                                    waterB),waterB),
             # Adjusted incoming energy and moisture effect on bark
             Qi = max(Qi-drain,0),
             rhoM <- (waterB+waterB^2)*barkDensity,
             #Bark thermal values
             cpBark = ((1105+4.85*tempAir)*(1-waterB)+waterB*4185+1276*waterB)/1000,
             kAir = 0.00028683*(tempAir+273.15)^0.7919,
             kBark = (2.104*barkDensity+5.544*rhoM+3.266*tempB-166.216)*10^-4,
             # Fourier conduction through bark
             fourierB = pmin(Qi,(Area * kBark * (tempAir - tempB)) / bark),
             tempBark = pmax(tempB, (0.001 * fourierB / (mBark * cpBark) + tempB)),
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWaterW = water*mWood,
             # Energy removed by current water quantity
             drain = ifelse(tempW>95,
                            ifelse(waterW>0,mWaterW*2256400,0),0),
             # Change in proportion wood water this step
             waterW = ifelse(tempW>95,
                             ifelse(waterW>0,max(0,waterW-((Qi/2256400)/mWaterW)),
                                    waterW),waterW),
             # Adjusted incoming energy
             Qi = max(Qi-drain,0),
             #Wood thermal values
             cpWood = 1.08+0.00408*(100*waterW)+0.00253*tempAir+0.0000628*(100*waterW)*tempAir,
             kWood = kWood(tempW, woodDensity, kAir),
             # Fourier conduction through wood
             fourierW = pmin(fourierB,(Area * kWood * (tempBark - tempW)) / wood),
             tempWood = pmax(tempW, (0.001 * fourierW / (mWood * cpWood) + tempW)),
             # Mortality
             mortality = ifelse(tempWood < 67.5-0.3017*30.17*RH, 0, 1),
             #Outgoing wood
             fourierOW = pmin(0,(Area * kWood * (tempBark - tempWood)) / wood),
             tempW = pmax(tempW, tempWood + (0.001 * fourierOW / (mWood * cpWood))),
             #Outgoing bark
             fourierOB = pmin(0,(Area * kBark * (tempAir - tempBark)) / bark),
             tempB = pmax(tempB, tempBark + (0.001 * fourierOB / (mBark * cpBark))),
             qrO = pmin(0,Area*0.86*0.0000000000567*((tempAir+273.15)^4 - (tempB+273.15)^4)),
             qR = qr + qrO,
             Q = qc + qR)
    Ca <- rbind(Ca, Cb)

    woodM <- quantile(Cb$tempW, percentile)
    barkM <- quantile(Cb$tempB, percentile)
    waterW <- quantile(Cb$waterW, percentile)
    waterB <- quantile(Cb$waterB, percentile)

    t = t + 1
    Horiz = Horiz - ROS
  }

  # Create and export table
  Ca <- Ca %>%
    select(t, repId, ros_kph, tempAir, tempBark, tempWood, waterB, waterW,
           cpBark, cpWood, att, qc, qr, qrO, qR, Q, kAir, kBark, kWood, mortality)
  write.csv(Ca, "hollow.csv")
  return(Ca)
}

#####################################################################

#' Fire risk for an animal sheltering underground
#'
#' Calculates the likelihood of mortality to an animal
#' caused by an approaching fire front
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
#' Finds animal mortality within a hollow based on the maximum tolerable temperature for a given
#' vapour pressure deficit, based on data from Lawrence, G. E.
#' Ecology of vertebrate animals in relation to chaparral fire in the Sierra Nevada foothills.
#' Ecology 47, 278–291 (1966)
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
#' @param depth The depth at which the animal shelters beneath the soil
#' @param soilTemp The starting temperature under the ground (deg C)
#' @return dataframe
#' @export underground.csv
#' @examples underground(Surf, Plant, diameter = 6, percentile = 0.5, RH = 0.2, moisture = 0.2, distance = 50,
#' Pressure = 1013.25, Altitude = 0, texture = "clay", peat = 0.1, depth = 0.1, soilTemp = 25)
#'

underground <- function(Surf, Plant, diameter = 6, surface = 677, percentile = 0.5, RH = 0.2,
                        moisture = 0.2, distance = 50, trail = 300, var = 10, Pressure = 1013.25,
                   Altitude = 0, texture = "clay", peat = 0.1, grain = "fine", unfrozen = 1, depth = 0.1, soilTemp = 25)

{
  # Collect step distance, time, and total distance
  residence <- 0.871*diameter^1.875
  ROS <- mean(Surf$ros_kph)/3.6
  Ta <- round(distance/ROS+residence)
  TIME <- Ta + trail
  Horiz <- distance
  HA <- abs(Horiz)

  # Description of the protection
  densityD <- denSoil(texture)
  mass <- depth * densityD
  R <- sqrt(1/pi)
  soilTemp <- soilTemp

  #Starting values
  Ca <- threat(Surf, Plant, HA, Height=0, var, Pressure, Altitude)%>%
    mutate(t = 1,
           #Convective transfer
           Re = (Plume_velocity*Density)/viscosity,
           h = ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888),
           #Incoming heat
           tempS = ifelse(Horiz <0, pmax(tempAir, surface), tempAir),
           qc = h * (tempS - soilTemp),
           att = tau(D=HA, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
           qr = qr*att,
           Qi = pmax(0, qc)+qr,
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*mass,
           # Energy removed by current water quantity
           drain = ifelse(soilTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           Qi = max(Qi-drain,0),
           #Thermal values
           cpSoil = cpSoil((soilTemp+273.15), texture, peat, moisture),
           saturation = satSoil(texture, moisture),
           kSoil = kSoil(texture, saturation, grain = "grain", unfrozen = unfrozen),
           # Fourier conduction
           fourier = ifelse(Horiz>0, pmin(Qi,(kSoil * pmax(0,tempS - soilTemp)) / depth),
                            (kSoil * pmax(0,tempS - soilTemp)) / depth),
           tempSoil = pmax(soilTemp, (fourier / (mass * cpSoil) + soilTemp)),
           # Change in proportion wood water this step
           moisture = ifelse(moisture>0,max(0,moisture-((fourier/2256400)/mWater)),
                             moisture),
           # Mortality
           mortality = ifelse(tempSoil < 67.5-30.17*RH, 0, 1),
           #Outgoing
           fourierO = pmin(0,(kSoil * (tempS - tempSoil)) / depth),
           soilTemp = tempSoil + (fourierO / (mass * cpSoil)),
           qrO = pmin(0,0.0000000000567*((tempS+273.15)^4 - (soilTemp+273.15)^4)),
           qR = qr + qrO,
           Q = qc + qR)

  soilTemp <- quantile(Ca$soilTemp, percentile)
  moisture <- quantile(Ca$moisture, percentile)

  # Advance one second's travel
  Horiz = Horiz - ROS
  HA <- abs(Horiz)

  # Loop through each time step and collect outputs
  for(t in 2:TIME){
    Cb <-threat(Surf, Plant, HA, Height=0, var, Pressure, Altitude) %>%
      mutate(t = t,
             #Convective transfer
             Re = (Plume_velocity*Density)/viscosity,
             h = ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888),
             #Incoming
             tempS = ifelse(t>Ta, tempAir, ifelse(Horiz <0, pmax(tempAir, surface), tempAir)),
             qc = h * (tempS - soilTemp),
             att = tau(D=HA, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
             qr = qr*att,
             Qi = pmax(0, qc)+qr,
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moisture*mass,
             # Energy removed by current water quantity
             drain = ifelse(soilTemp>95,
                            ifelse(moisture>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             Qi = max(Qi-drain,0),
             #Thermal values
             cpSoil = cpSoil((soilTemp+273.15), texture, peat, moisture),
             saturation = satSoil(texture, moisture),
             kSoil = kSoil(texture, saturation, grain = "grain", unfrozen = unfrozen),
             # Fourier conduction
             fourier = ifelse(Horiz>0, pmin(Qi,(kSoil * pmax(0,tempS - soilTemp)) / depth),
                              (kSoil * pmax(0,tempS - soilTemp)) / depth),
             tempSoil = pmax(soilTemp, (fourier / (mass * cpSoil) + soilTemp)),
             # Change in proportion wood water this step
             moisture = ifelse(moisture>0,max(0,moisture-((fourier/2256400)/mWater)),
                               moisture),
             # Mortality
             mortality = ifelse(tempSoil < 67.5-0.3017*30.17*RH, 0, 1),
             #Outgoing
             fourierO = pmin(0,(kSoil * (tempS - tempSoil)) / depth),
             soilTemp = tempSoil + (fourierO / (mass * cpSoil)),
             qrO = pmin(0,0.0000000000567*((tempS+273.15)^4 - (soilTemp+273.15)^4)),
             qR = qr + qrO,
             Q = qc + qR)
    Ca <- rbind(Ca, Cb)

    soilTemp <- quantile(Cb$soilTemp, percentile)
    moisture <- quantile(Cb$moisture, percentile)

    t = t + 1
    Horiz = Horiz - ROS
    HA <- abs(Horiz)
  }

  # Create table
  Ca <- Ca %>%
    select(t, repId, ros_kph, tempAir, tempS, tempSoil, moisture,
           cpSoil, kSoil, att, qc, qr, qrO, qR, Q, fourier, fourierO, mortality)
  write.csv(Ca, "underground.csv")
  return(Ca)
}
#####################################################################
#' Configuration or view factor
#'
#' Calculates the configuration factor for radiative heat transfer as used
#' in AS3959 for bushfire risk assessment of built structures
#'
#' @param Fl Flame length (metres)
#' @param Fa Flame angle (radians)
#' @param D Distance to the receiver (metres)
#' @param H Height of the receiver (metres)
#' @param S Slope (degrees)
#' @return Phi - the configuration factor
#' @examples phi(1, 1.3, 20, 5, 2)

phi <- function(Fl, Fa, S, D, H)
{
  # Distances
  Fx <- (Fl/2*cos(Fa))
  Sep <- abs(D - Fx)
  Sl <- (S * pi)/180

  # Configuration factor
  X1 <- (Fl*sin(Fa)-Fx*tan(Sl)-D*tan(Sl)-H)/Sep
  X2 <- (H+Sep*tan(Sl))/Sep
  Y1 <- 50/Sep
  X1a <- sqrt(1+X1^2)
  X2a <- sqrt(1+X2^2)
  Y1a <- sqrt(1+Y1^2)
  return(ifelse(Fl<0.1, 0, ifelse(Sep<=0, 1,
        (1/pi)*((X1/X1a)*atan(Y1/X1a)+(Y1/Y1a)*atan(X1/Y1a)+
        (X2/X2a)*atan(Y1/X2a)+(Y1/Y1a)*atan(X2/Y1a)))))
}

#####################################################################
#' Atmospheric attenuation
#'
#' Calculates the atmospheric transmissivity of radiation as used
#' in AS3959 for bushfire risk assessment of built structures. Equations from
#' Fuss, S P & Hamins, A
#' An estimate of the correction applied to radiant flame measurements due to
#' attenuation by atmospheric CO2and H2O.
#' Fire Saf. J. 37, 181–190 (2002)
#'
#' @param D Distance from flame to receiver (metres)
#' @param flameTemp Integrated temperature across the flame (kelvin)
#' @param temperature Ambient air temperature (kelvin)
#' @param rh Relative humidity (proportion)
#' @return tau - the atmospheric attenuation
#' @examples tau(D = 200, flameTemp = 1300, temperature = 288, rh = 0.51)

tau <- function(D = 200, flameTemp = 1300, temperature = 288, rh = 0.51)
{

  a0 <-1.486-0.002003*temperature+0.0000468*flameTemp-0.06052*rh
  a1 <-0.01225-0.000059*temperature+0.00000166*flameTemp-0.001759*rh
  a2 <--0.0001489+0.0000006893*temperature-0.00000001922*flameTemp+0.00002092*rh
  a3 <-0.0000008381-0.000000003823*temperature+0.00000000010511*flameTemp-0.0000001166*rh
  a4 <--0.000000001685+0.000000000007637*temperature-0.0000000000002085*flameTemp+0.0000000002350*rh
  return(a0+a1*D+a2*D^2+a3*D^3+a4*D^4)
}

#####################################################################
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
#' @export Iso.csv A .csv file listing the impacts from each model run
#' @examples flora(runs, IP, Param, 70)

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
  surf <- ns%>%
    mutate(Alpha = 1/(2*lengthSurface^2),
           C = 950*lengthSurface*exp(-Alpha*lengthSurface^2),
           pAlpha = abs(C/(Test-temperature)),
           R = pAlpha*cos(angleSurface),
           El = R * tan((slope_degrees * pi)/180),
           ht = sin(angleSurface)*pAlpha - El
    )%>%
    select(repId, temperature, slope_degrees, wind_kph, ht, ns, e, m, c)

  #Heating from plants
  plant<-Plant%>%
    left_join(surf)%>%
    mutate(Angle = atan((y1-y0)/(x1-x0)),
           Alpha = 1/(2*flameLength*(flameLength-length)),
           C = 950*flameLength*exp(-Alpha*(flameLength-length)^2),
           pAlpha = abs(C/(Test-temperature)),
           Reach = pAlpha*cos(Angle),
           El = Reach * tan((slope_degrees * pi)/180),
           htP = sin(Angle)*pAlpha + y0 - El)%>%
    select(repId, htP)%>%
    group_by(repId) %>%
    summarize_all(max)%>%
    right_join(surf)
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
#' @return dataframe
#' @export soilHeating.csv
#' @examples soil(Surf, Plant, diameter = 6, percentile = 0.5, RH = 0.2, moisture = 0.2, distance = 50,
#' Pressure = 1013.25, Altitude = 0, texture = "clay", peat = 0.1, soilTemp = 25)
#'

soil <- function(Surf, Plant, diameter = 6, surface = 677, percentile = 0.95, RH = 0.2,
                        moisture = 0.1, distance = 50, trail = 600, var = 10, Pressure = 1013.25,
                        Altitude = 0, texture = "sand", peat = 0, grain = "fine", unfrozen = 1, soilTemp = 25)

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

  write.csv(Ca, "soilHeating.csv")
  return(Ca)
}


#####################################################################
#' Summary table of surface results
#'
#' Summarises FRaME generated surface fire behaviour by RepId
#'
#' @param Surface The dataframe res$SurfaceResults
#' @return surface dataframe
#' @examples surf(res$SurfaceResults)

surf <- function(surface)
{

return(surface %>%
         group_by(repId) %>%
         mutate(lengthSurface = flameLength,
                heightSurface = flameHeight,
                angleSurface = flameAngle)%>%
         select(repId, lengthSurface, heightSurface, angleSurface) %>%
         summarize_all(max))
}

#####################################################################
#' Summary table of stratum results
#'
#' Summarises FRaME generated fire behaviour by stratum and RepId
#'
#' @param flames The dataframe res$FlameSummaries
#' @param sites The dataframe res$Sites
#' @param ros The dataframe res$ROS
#' @param surf The dataframe produced by the extension surf
#' @return stratum dataframe
#' @examples stratum(res$FlameSummaries, res$Sites, res$ROS, res$SurfaceResults)

stratum <- function(flames, sites, ros, surf)
{
  y <- ros%>%
    select(repId, level, ros)
  z <- flames %>%
    select(repId, level, flameLength, flameAngle, flameHeight)

  a <- y %>%
    left_join(z) %>%
    left_join(sites) %>%

    # Strata without ros will end up with NA values
    # after doing the join above. Convert these missing values to zero.
    mutate(ros = ifelse(is.na(ros), 0.0, ros),
           flameHeight = ifelse(is.na(flameHeight), 0.0, flameHeight),
           flameLength = ifelse(is.na(flameLength), 0.0, flameLength),
           flameAngle = ifelse(is.na(flameAngle), 0.0, flameAngle),
           #PATCH TO COVER MOISTURE EXTINCTION UNTIL FIXED IN SCALA
           #Duplicate DFMC, then create binary for spread/no spread
           extinct = deadFuelMoistureProp,
           extinct = ifelse(extinct == 0.199, 0.0, 1.0)) %>%
    select(repId, level, fuelLoad, flameHeight, flameLength, flameAngle, ros, windSpeed,
           deadFuelMoistureProp, temperature, slope, extinct) %>%
    mutate(oHorizon = fuelLoad * 10,
           slope_degrees = slope * 180 / pi,
           flameA_degrees = flameAngle * 180 / pi,
           ros_kph = extinct * ros * 3.6,
           heightPlant = flameHeight * extinct,
           lengthPlant = flameLength * extinct,
           wind_kph = windSpeed * 3.6,
           spread = ifelse(ros > 0, 1, 0),
           has.flame = spread + (extinct * flameHeight) > 0)

  # Add in surface flame descriptors
  rep <- max(a$repId)
  st <- as.numeric(count(a))/rep
  i <- 1
  for(loop in 1:rep) {
    a$flameHeight[i]=surf$heightSurface[loop]
    a$flameLength[i]=surf$lengthSurface[loop]
    a$flameAngle[i]=surf$angleSurface[loop]
    i <- i + st
  }
  return(a)
}

#####################################################################
#' Summary table of fire behaviour
#'
#' Summarises FRaME generated fire behaviour by RepId
#'
#' @param surf The dataframe surf
#' @param stratum The dataframe stratum
#' @return summary dataframe
#' @examples summary(stratum, surf)

summary <- function(stratum, surf)
{
  return(stratum %>%
           select(repId, slope_degrees, wind_kph, deadFuelMoistureProp, temperature,
                  heightPlant, lengthPlant, flameAngle, ros_kph, extinct) %>%
           group_by(repId) %>%
           summarize_all(max) %>%
           left_join(surf) %>%
           mutate(fh = pmax(heightSurface, heightPlant) * extinct,
                  fl = pmax(lengthSurface, lengthPlant) * extinct,
                  zeta = 2.5*ros_kph,
                  epsilon = 1-exp(-0.72*zeta)))
}

#####################################################################
#' Representative flame descriptors
#'
#' Summarises FRaME generated flame segments into a combined,
#' representative plant flame for each repId where plants ignited
#'
#' @param IP The dataframe res$IgnitionPaths
#' @return dataframe
#' @examples RepFlame(res$IgnitionPaths)

repFlame <- function(IP)
{
top <- IP %>%
  mutate(angle = atan((y1 - y0)/(x1 - x0)),
         repHeight = flameLength*sin(angle)+y0)%>%
  group_by(repId) %>%
  summarize_all(max) %>%
  select(repId, repHeight)

repFlame <- IP %>%
  mutate(repAngle = atan((y1 - y0)/(x1 - x0))
  ) %>%
  select(repId, repAngle)%>%
  group_by(repId) %>%
  summarize_all(mean) %>%
  left_join(top) %>%
  mutate(repLength = repHeight/sin(repAngle)) %>%
  select(repId, repHeight, repLength, repAngle)%>%
  right_join(IP)

return(repFlame)
}


#####################################################################
#' Stratum descriptors from a param file
#'
#' For each stratum, finds mean crown width, plant separation, and number of species
#'
#' @param Param A parameter dataframe used for FRaME,
#' such as produced using readLegacyParamFile
#' @return dataframe of stratum descriptors
#' @examples strata(Param)

strata <- function(Param)
{
  #Number of strata
  StL <- count(Param)-13
  StN <- Param$stratum[max(StL$n)]

  #Count species per stratum
  Sp <- numeric(StN)
  for(sn in 1:StN){
    strat <- filter(Param, stratum == sn)
    strat <- na.omit(strat)
    Sp[sn] <- (as.numeric(max(strat$species))+1)-as.numeric(min(strat$species))
  }

  #COLLECT DIMENSIONS
  width <- Param[Param$param == "w", ]
  comp <- Param[Param$param == "composition", ]
  sep <- Param[Param$param == "plantSeparation", ]
  peak <- Param[Param$param == "hp", ]
  top <- Param[Param$param == "ht", ]
  edge <- Param[Param$param == "he", ]
  base <- Param[Param$param == "hc", ]
  level <- Param[Param$param == "levelName", ]
  name <- Param[Param$param == "levelName", ]

  #BUILD TABLE
  n <- as.data.frame(list('stratum'=name$stratum, 'name'=name$value, 'speciesN'=Sp))
  s <- as.data.frame(list('stratum'=width$stratum, 'comp'=comp$value, 'width'=width$value, 'Hp'=peak$value,
                          'Ht'=top$value, 'He'=edge$value, 'Hc'=base$value))%>%
    mutate(Co = as.numeric(as.character(comp)),
           Ww = as.numeric(as.character(width))*as.numeric(as.character(comp)),
           Wp = as.numeric(as.character(Hp))*as.numeric(as.character(comp)),
           Wt = as.numeric(as.character(Ht))*as.numeric(as.character(comp)),
           We = as.numeric(as.character(He))*as.numeric(as.character(comp)),
           Wc = as.numeric(as.character(Hc))*as.numeric(as.character(comp)),
           top = pmax(Wp,Wt),
           base = pmin(We,Wc))%>%
    group_by(stratum) %>%
    summarize_if(is.numeric,sum)%>%
    mutate(width = Ww/Co,
           top = top/Co,
           base = base/Co)%>%
    left_join(sep, by = "stratum")%>%
    mutate(separation = as.numeric(value),
           cover = width^2/separation^2)%>%
    select(stratum, separation, cover, width, base, top)

  strata <- as.data.frame(s)%>%
    left_join(n, by="stratum")
  return(strata)
}
#####################################################################
#' Species descriptors from a param file
#'
#' Finds dimensions and moisture of each species
#'
#' @param Param A parameter dataframe used for FRaME,
#' such as produced using readLegacyParamFile
#' @return dataframe of species descriptors
#' @examples species(Param)

species <- function(Param)
{
  #Collect traits
  sp <- Param[Param$param == "name", ]
  lfmc <- Param[Param$param == "liveLeafMoisture", ]
  Peak <- Param[Param$param == "hp", ]
  Top <- Param[Param$param == "ht", ]
  Edge <- Param[Param$param == "he", ]
  Base <- Param[Param$param == "hc", ]
  Width <- Param[Param$param == "w", ]

  species <- as.data.frame(list('name'=sp$value, 'hp'=as.numeric(Peak$value),'ht'=as.numeric(Top$value),
                                'hc'=as.numeric(Base$value), 'he'=as.numeric(Edge$value),
                                'w'=as.numeric(Width$value), 'lfmc'=as.numeric(lfmc$value))) %>%
    mutate(htR = ht/hp,
           hcR = hc/hp,
           heR = he/hp,
           wR = w/hp)
  return(species)
}

#####################################################################
#' Thermal conductivity of dry wood
#'
#'Model drawn from Kollmann, F. F. P. & Cote, W. A.
#'Principles of wood science and technology I. Solid wood. (Springer-Verlag, 1968)
#'
#' @param T Air temperature (deg C)
#' @param rhoW Timber density (kg/m3)
#' @param kAir Thermal conductivity of the air at temperature T (W/M.K)
#' @return Thermal conductivity of wood (W/mK)
#' @examples kWood(50, 700)

kWood <- function(T=100, rhoW=700, kAir = 0.026)
{
  T <- T+273.15
  bridge <- 0.0019*T+0.0503
  r <- 1-(rhoW/1500)
  kA <- (1-r)*0.766+r*kAir
  kB <- 1/((1-r)/0.43+(r/kAir))
  return(bridge*kA+(1-bridge)*kB)
}

#####################################################################

#' Thermal conductivity of soil
#'
#'Model drawn from Johansen, O.
#'Thermal conductivity of soils. PhD Thesis (University of Trondheim, 1971)
#'Modified by Farouki, O. Thermal properties of soils. (Trans Tech, 1986)
#'
#'Porosity taken from Rawls, W. J., Brakensiek, D. L. & Saxton, K. E.
#'Estimation of soil water properties. Transactions of the ASAE 25, 1316–1320 & 1328 (1982).
#'
#' @param texture Soil texture. Allowable values are: "sand", "loamy sand", "sandy loam", "sandy clay loam",
#' "sand clay", "loam", "clay loam", "silt loam", "clay", "silty clay", "silty clay loam", "silt"
#' @param saturation Water saturation of the soil, value between 0.1 and 1
#' @param grain Allowable values are "fine" or "coarse"
#' @param unfrozen Proportion of soil unfrozen, between 0 and 1
#' @return Thermal conductivity of soil (W/mK)
#' @examples kSoil("loam", 0.3)

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
  kDry <- (0.135*densityD+64.7)/(2700-0.947*densityD)
  # Solids thermal conductivity
  kS <- 7.7^quartz*minerals^(1-quartz)
  # Saturated thermal conductivity
  kSat <- kS^(1-porosity)*2.2^(porosity-unfrozen)*0.57^unfrozen

  return(kersten*(kSat-kDry)+kDry)
}


#####################################################################

#' Specific heat of soil
#'
#' Finds volumetric specific heat from the mineral, organic and water components of the soil
#'
#'Specific heats of soil components estimated from figures 108 & 111 in
#'Farouki, O. Thermal properties of soils. (Trans Tech, 1981).
#'
#'Water specific heat 4185 J/kg.K
#'
#' @param temp soil temperature in K
#' @param texture Soil texture. Allowable values are: "sand", "loamy sand", "sandy loam", "sandy clay loam",
#' "sand clay", "loam", "clay loam", "silt loam", "clay", "silty clay", "silty clay loam", "silt"
#' @param peat Organic proportion of the soil
#' @param moisture Water proportion ODW of the soil
#' @return Specific heat of soil (J/kg.K)
#' @examples cpSoil(300, "loam", 0.3, 0.2)

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

#' Soil saturation
#'
#' Finds saturation from ODW moisture and texture
#'
#'Field capacity of soils taken from Salter, P. J. & Williams, J. B.
#'The influence of texture on the moisture characteristics of soil.
#'V. Relationships between particle-size composition and moisturecontents
#'at the upper and lower limits of available-water. J. Soil Sci. 20, 126–131 (1969).
#'
#'
#' @param texture Soil texture. Allowable values are: "sand", "loamy sand", "sandy loam", "sandy clay loam",
#' "sand clay", "loam", "clay loam", "silt loam", "clay", "silty clay", "silty clay loam", "silt"
#' @param moisture Water proportion ODW of the soil
#' @examples satSoil("loam", 0.1)

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


#####################################################################

#' Soil density
#'
#' Finds soil density from texture
#'
#'Porosity of soils taken from Rawls, W. J., Brakensiek, D. L. & Saxton, K. E.
#'Estimation of soil water properties. Transactions of the ASAE 25, 1316–1320 & 1328 (1982).
#'
#'Density equation is from Peters-Lidard, C. D., Blackburn, E., Liang, X. & Wood, E. F.
#'The effect of soil thermal conductivity parameterization on surface energy fluxes and temperatures.
#'J. Atmos. Sci. 55, 1209–1224 (1998).
#'
#'
#' @param texture Soil texture. Allowable values are: "sand", "loamy sand", "sandy loam", "sandy clay loam",
#' "sand clay", "loam", "clay loam", "silt loam", "clay", "silty clay", "silty clay loam", "silt"
#' @examples denSoil("loam")

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
#' @return dataframe
#' @examples cambium(Surf, Plant, percentile = 0.95, Height = 0.1, woodDensity = 700, barkDensity = 500,
#'  bark = 0.04, comBark = 700, resBark = 45, cambThick = 0.01, xylemThick = 0.01, RH = 0.2,
#'  moisture = 0.2, distance = 50, trail = 100, var = 10, diameter = 20, Pressure = 1013.25,
#'  Altitude = 0, hollowTemp = 25, necT = 60, surfDecl = 10)

cambium <- function(Surf, Plant, percentile = 0.95, Height = 0.1, woodDensity = 700, barkDensity = 500,
                   bark = 0.04, comBark = 700, resBark = 45, cambThick = 0.01, xylemThick = 0.01, RH = 0.2,
                   moisture = 1, bMoisture = 0.2, distance = 50, trail = 100, var = 10, diameter = 20, Pressure = 1013.25,
                   Altitude = 0, startTemp = 25, necT = 60, surfDecl = 10)
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

#' Air temperature above ambient at the tree bole behind the flame front
#'
#' Dynamic air temperature at bole height, declining flame length exponentially
#'
#' Air temperature is modelled from dynamic flame segments using
#' Weber R.O., Gill A.M., Lyons P.R.A., Moore P.H.R., Bradstock R.A., Mercer G.N. (1995)
#' Modelling wildland fire temperatures. CALMScience Supplement, 4, 23–26.
#'
#' pAlphas is set to bole height - depth of surface litter
#'
#' @param lengthSurface Surface flame length in m
#' @param residence Surface flame residence in seconds
#' @param surfDecl Constant for the rate of decline in flame length
#' @param h Height of the point being measured
#' @param depth Depth of surface litter in m
#' @param t Time in seconds behind the flame front
#' @return temperature
#' @examples bole(2, 300, 0.05, 0.1, 10, 200)

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

#' Finds non-deterministic fire behaviour predictions
#'
#' uses Monte Carlo, for set slope, DFMC & temperature
#' wind, LFMC, heights and leaf traits randomly varied within defined ranges
#'
#' @param base-params Parameter input table
#' @param out.db Name of the exported database
#' @param jitters Number of repetitions
#' @param s Slope of the site (degrees)
#' @param l Variation around input leaf dimensions
#' @param tempM Mean air teperature (degC)
#' @param Dm Proportion dead fuel moisture
#' @param windl Lower limit of wind speed (km/h)
#' @param windst Range of wind speeds tested (km/h)
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height
#' @examples run_ffm_ND(Param, "MC", jitters = 100, s = 0, l = 0.1,
#' tempm = 30, Dm = 0.07, windl = 0, windst = 50,
#' Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41)

run_ffm_ND <- function(base.params, out.db = "out_mc.db", jitters = 100, s = 0, l = 0.1,
                       tempm = 30, Dm = 0.07, windl = 0, windst = 60,
                       Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41)
{

  # Collect original parameters
  strata <- strata(base.params)
  species <- species(base.params)

  # This creates a reference to a Scala Database object
  # which will handle writing of model results to the
  # output file.
  #
  db <- ffm_create_database(out.db,
                            delete.existing = TRUE,
                            use.transactions = TRUE)

  # Update parameter table
  tbl <- base.params %>%
    ffm_set_site_param("slope", s, "deg") %>%
    ffm_set_site_param("temperature", tempm, "degc") %>%
    ffm_set_site_param("deadFuelMoistureProp", Dm)

  # Run the model, updating the base parameter table
  # with MC values at each iteration

  pbar <- txtProgressBar(max = jitters, style = 3)
  for (i in 1:jitters) {

    # Choose wind speed within a uniform distribution
    winds <- runif(1)*windst+windl

    # Choose random leaf traits within a uniform distribution
    tbl <- ffm_param_variance(tbl,max.prop = l, method = "uniform")%>%
      ffm_set_site_param("windSpeed", winds, "km/h")

    # MODIFY SPECIES-SPECIFIC TRAITS
    SpeciesN <- 1
    SpeciesP <- 1
    si <- 1
    #Replace old parameters
    StN <- as.numeric(count(strata))

    # Change lfmc values

    # If within the cover range, set lfmc to random value
    for (loop in 1:StN) {
      if (runif(1)<=strata$cover[si]) {
        for (t in 1:strata$speciesN[si]) {
          # Introduce Pm as a multiplier of LFMC
          Mrand <- Pm*rtnorm(n = 1, mean = species$lfmc[SpeciesN], sd = Ms,
                             a = species$lfmc[SpeciesN]/Mr, b = species$lfmc[SpeciesN]*Mr)
          tbl <- tbl%>%
            ffm_set_species_param(si, SpeciesN, "liveLeafMoisture", Mrand)
          (SpeciesN = SpeciesN+1)
        }
        # If outside of cover range, set LFMC to value that will not ignite
      } else {
        for (f in 1:strata$speciesN[si]) {
          tbl <- tbl%>%
            ffm_set_species_param(si, SpeciesN, "liveLeafMoisture", 100)
          SpeciesN = SpeciesN+1
        }
      }
      #Set other species params here
      for (p in 1:strata$speciesN[si]) {
        peak <- rtnorm(n = 1, mean = species$hp[SpeciesP], sd = Hs,
                       a = species$hp[SpeciesP]/Hr, b = species$hp[SpeciesP]*Hr)
        #Give other heights as proportion of peak
        tbl <- tbl%>%
          ffm_set_species_param(si, SpeciesP, "hp", peak) %>%
          ffm_set_species_param(si, SpeciesP, "ht", peak*species$htR[SpeciesP]) %>%
          ffm_set_species_param(si, SpeciesP, "he", peak*species$heR[SpeciesP]) %>%
          ffm_set_species_param(si, SpeciesP, "hc", peak*species$hcR[SpeciesP])
        SpeciesP = SpeciesP+1
      }
      si = si+1
    }

    # Run the model
    ffm_run(tbl, db)
    setTxtProgressBar(pbar, i)
  }

  # Tell the Scala database object to close the output file.
  db$close()

  cat("Finished.  Output written to", out.db)
}

########################################################################

#' Finds non-deterministic fire behaviour predictions for the Silvertop validation
#'
#' uses Monte Carlo, for set slope, DFMC, temperature and shrub percent cover
#' wind, LFMC, heights and leaf traits randomly varied within defined distributions
#'
#' @param base-params Parameter input table
#' @param out.db Name of the exported database
#' @param jitters Number of repetitions
#' @param s Slope of the site (degrees)
#' @param l Variation around input leaf dimensions
#' @param tempM Mean air teperature (degC)
#' @param Dm Proportion dead fuel moisture
#' @param f Surface fuel load (t/ha)
#' @param fLine Fireline length (m)
#' @param windm Mean wind speed (km/h)
#' @param windsd Standard deviation of wind speed (km/h)
#' @param cover Shrub cover (percent)
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height


silvertop <- function(base.params, out.db = "out_mc.db", jitters = 100, s = 0, l = 0.1,
                      tempm = 30, Dm = 0.07, f = 10, fLine = 100, windm = 10, windsd = 1, cover = 30,
                      Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41)
{

  #Collect original parameters
  strata <- strata(base.params)
  species <- species(base.params)

  # Find plant separation and replace value in table
  sepe <- sqrt(as.numeric(strata$width[2])^2/(0.01*cover))
  strata$separation[strata$name=="elevated"] = sepe
  strata$cover[strata$name=="elevated"] = cover/100

  # This creates a reference to a Scala Database object
  # which will handle writing of model results to the
  # output file.
  #
  db <- ffm_create_database(out.db,
                            delete.existing = TRUE,
                            use.transactions = TRUE)


  # Update parameter table
  tbl <- base.params %>%
    ffm_set_site_param("slope", s, "deg") %>%
    ffm_set_site_param("temperature", tempm, "degc") %>%
    ffm_set_site_param("deadFuelMoistureProp", Dm)%>%
    ffm_set_site_param("fuelLoad", f)%>%
    ffm_set_site_param("fireLineLength", fLine, "m")


  # Run the model, updating the base parameter table
  # with MC values at each iteration

  pbar <- txtProgressBar(max = jitters, style = 3)
  for (i in 1:jitters) {

    # Choose wind speed within a normal distribution
    winds <- rnorm(n=1, mean = windm, sd = windsd)

    # Choose random leaf traits within a uniform distribution
    tbl <- ffm_param_variance(tbl,max.prop = l, method = "uniform")%>%
      # Set wind speed
      ffm_set_site_param("windSpeed", winds, "km/h")

    # MODIFY SPECIES-SPECIFIC TRAITS
    # Set LFMC, and block some plants from ignition based on cover

    SpeciesN <- 1
    si <- 1
    StN <- as.numeric(count(strata))

    # Change lfmc values

    # If within the cover range, set lfmc to random value
    for (loop in 1:StN) {
      if (runif(1)<=strata$cover[si]) {
        for (t in 1:strata$speciesN[si]) {
          # Introduce Pm as a multiplier of LFMC
          Mrand <- Pm*rtnorm(n = 1, mean = species$lfmc[SpeciesN], sd = Ms,
                             a = species$lfmc[SpeciesN]/Mr, b = species$lfmc[SpeciesN]*Mr)
          tbl <- tbl%>%
            ffm_set_species_param(si, SpeciesN, "liveLeafMoisture", Mrand)
          (SpeciesN = SpeciesN+1)
        }
        # If outside of cover range, set LFMC to value that will not ignite
      } else {
        for (f in 1:strata$speciesN[si]) {
          tbl <- tbl%>%
            ffm_set_species_param(si, SpeciesN, "liveLeafMoisture", 100)
          SpeciesN = SpeciesN+1
        }
      }
      si = si+1
    }
    # Modify shrub heights
    numE <- strata$speciesN[strata$name=="near surface"]+1
    endE <- strata$speciesN[strata$name=="elevated"]+strata$speciesN[strata$name=="near surface"]
    for (p in numE:endE) {
      peak <- rtnorm(n = 1, mean = species$hp[numE], sd = Hs,
                     a = species$hp[numE]/Hr, b = species$hp[numE]*Hr)
      # Give other heights as proportion of peak
      tbl <- tbl%>%
        ffm_set_species_param(2, numE, "hp", peak) %>%
        ffm_set_species_param(2, numE, "ht", peak*species$htR[numE]) %>%
        ffm_set_species_param(2, numE, "he", peak*species$heR[numE]) %>%
        ffm_set_species_param(2, numE, "hc", peak*species$hcR[numE])
      numE = numE+1
    }

    # Run the model
    ffm_run(tbl, db)
    setTxtProgressBar(pbar, i)
  }

  # Tell the Scala database object to close the output file.
  db$close()

  cat("Finished.  Output written to", out.db)
  res<-ffm_db_load("MC")

  # Build tables
  surf <- surf(res$SurfaceResults)
  x <- stratum(res$FlameSummaries, res$Sites, res$ROS, surf)
  runs <- summary(x, surf)

  SFC <- filter(x, level == "Surface")
  sfc <- round(mean(SFC$flameHeight), 2)
  NS <- filter(x, level == "NearSurface")
  ns <- round(mean(NS$flameHeight), 2)
  E <- filter(x, level == "Elevated")
  e <- round(mean(E$flameHeight), 2)
  M <- filter(x, level == "MidStorey")
  m <- round(mean(M$flameHeight), 2)
  C <- filter(x, level == "Canopy")
  c <- round(mean(C$flameHeight), 2)

  ROSm <- round(mean(runs$ros_kph),4)
  ROSsd <- round(sd(runs$ros_kph),4)
  ROSse <- round(ROSsd/sqrt(jitters),4)
  Hm <- round(mean(runs$fh),2)
  Hsd <- round(sd(runs$fh),2)
  Hse <- round(Hsd/sqrt(jitters),4)

  out <- as.data.frame(list('ROS_mean' = ROSm, 'ROS_se' = ROSse, 'ROS_sd' = ROSsd, 'FH_mean'=Hm, 'FH_se' = Hse, 'FH_sd'=Hsd,
                            'FH_s'=sfc, 'FH_ns'=ns, 'FH_e'=e, 'FH_m'=m, 'FH_c'=c))
  return(out)
}

########################################################################

#' Finds non-deterministic fire behaviour predictions for the Eucumbene validation
#'
#' uses Monte Carlo, for set slope, DFMC, temperature and shrub percent cover
#' wind, LFMC, heights and leaf traits randomly varied within defined distributions
#'
#' @param base-params Parameter input table
#' @param out.db Name of the exported database
#' @param jitters Number of repetitions
#' @param s Slope of the site (degrees)
#' @param l Variation around input leaf dimensions
#' @param tempM Mean air teperature (degC)
#' @param Dm Proportion dead fuel moisture
#' @param f Surface fuel load (t/ha)
#' @param fLine Fireline length (m)
#' @param windm Mean wind speed (km/h)
#' @param windsd Standard deviation of wind speed (km/h)
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Pm * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height


eucumbene <- function(base.params, out.db = "out_mc.db", jitters = 100, s = 0, l = 0.1,
                      tempm = 30, Dm = 0.07, f = 10, fLine = 100, windm = 10, windsd = 1,
                      Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41)
{

  #Collect original parameters
  strata <- strata(base.params)
  species <- species(base.params)

  # This creates a reference to a Scala Database object
  # which will handle writing of model results to the
  # output file.
  #
  db <- ffm_create_database(out.db,
                            delete.existing = TRUE,
                            use.transactions = TRUE)

  # Update parameter table
  tbl <- base.params %>%
    ffm_set_site_param("slope", s, "deg") %>%
    ffm_set_site_param("temperature", tempm, "degc") %>%
    ffm_set_site_param("deadFuelMoistureProp", Dm)%>%
    ffm_set_site_param("fuelLoad", f)%>%
    ffm_set_site_param("fireLineLength", fLine, "m")


  # Run the model, updating the base parameter table
  # with MC values at each iteration

  pbar <- txtProgressBar(max = jitters, style = 3)
  for (i in 1:jitters) {

    # Choose wind speed within a normal distribution
    winds <- rnorm(n=1, mean = windm, sd = windsd)

    # Choose random leaf traits within a uniform distribution
    tbl <- ffm_param_variance(tbl,max.prop = l, method = "uniform")%>%
      # Set wind speed
      ffm_set_site_param("windSpeed", winds, "km/h")

    # MODIFY SPECIES-SPECIFIC TRAITS
    # Set LFMC, and block some plants from ignition based on cover

    SpeciesN <- 1
    si <- 1
    StN <- as.numeric(count(strata))

    # Change lfmc values

    # If within the cover range, set lfmc to random value
    for (loop in 1:StN) {
      if (runif(1)<=strata$cover[si]) {
        for (t in 1:strata$speciesN[si]) {
          # Introduce Pm as a multiplier of LFMC
          Mrand <- Pm*rtnorm(n = 1, mean = species$lfmc[SpeciesN], sd = Ms,
                             a = species$lfmc[SpeciesN]/Mr, b = species$lfmc[SpeciesN]*Mr)
          tbl <- tbl%>%
            ffm_set_species_param(si, SpeciesN, "liveLeafMoisture", Mrand)
          (SpeciesN = SpeciesN+1)
        }
        # If outside of cover range, set LFMC to value that will not ignite
      } else {
        for (f in 1:strata$speciesN[si]) {
          tbl <- tbl%>%
            ffm_set_species_param(si, SpeciesN, "liveLeafMoisture", 100)
          SpeciesN = SpeciesN+1
        }
      }
      si = si+1
    }
    # Modify shrub heights
      peak <- rtnorm(n = 1, mean = species$hp[1], sd = Hs,
                     a = species$hp[1]/Hr, b = species$hp[1]*Hr)
      # Give other heights as proportion of peak
      tbl <- tbl%>%
        ffm_set_species_param(1, 1, "hp", peak) %>%
        ffm_set_species_param(1, 1, "ht", peak*species$htR[1]) %>%
        ffm_set_species_param(1, 1, "he", peak*species$heR[1]) %>%
        ffm_set_species_param(1, 1, "hc", peak*species$hcR[1])

    # Run the model
    ffm_run(tbl, db)
    setTxtProgressBar(pbar, i)
  }

  # Tell the Scala database object to close the output file.
  db$close()

  cat("Finished.  Output written to", out.db)
  res<-ffm_db_load("MC")

  # Build tables
  surf <- surf(res$SurfaceResults)
  x <- stratum(res$FlameSummaries, res$Sites, res$ROS, surf)
  runs <- summary(x, surf)

  SFC <- filter(x, level == "Surface")
  sfc <- round(mean(SFC$flameHeight), 2)
  E <- filter(x, level == "Elevated")
  e <- round(mean(E$flameHeight), 2)
  C <- filter(x, level == "Canopy")
  c <- round(mean(C$flameHeight), 2)

  ROSm <- round(mean(runs$ros_kph),4)
  ROSsd <- round(sd(runs$ros_kph),4)
  ROSse <- round(ROSsd/sqrt(jitters),4)
  Hm <- round(mean(runs$fh),2)
  Hsd <- round(sd(runs$fh),2)
  Hse <- round(Hsd/sqrt(jitters),4)

  out <- as.data.frame(list('ROS_mean' = ROSm, 'ROS_se' = ROSse, 'ROS_sd' = ROSsd, 'FH_mean'=Hm, 'FH_se' = Hse, 'FH_sd'=Hsd,
                            'FH_s'=sfc, 'FH_e'=e, 'FH_c'=c))
  return(out)
}

