#' Heat and plume produced by a fire front
#' 
#' Private function, documented for current development
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