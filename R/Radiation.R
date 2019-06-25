# Configuration or view factor
#
# Calculates the configuration factor for radiative heat transfer as used
# in AS3959 for bushfire risk assessment of built structures


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
# Atmospheric attenuation
#
# Calculates the atmospheric transmissivity of radiation as used
# in AS3959 for bushfire risk assessment of built structures. Equations from
# Fuss, S P & Hamins, A
# An estimate of the correction applied to radiant flame measurements due to
# attenuation by atmospheric CO2and H2O.
# Fire Saf. J. 37, 181–190 (2002)


tau <- function(D = 200, flameTemp = 1300, temperature = 288, rh = 0.51)
{

  a0 <-1.486-0.002003*temperature+0.0000468*flameTemp-0.06052*rh
  a1 <-0.01225-0.000059*temperature+0.00000166*flameTemp-0.001759*rh
  a2 <--0.0001489+0.0000006893*temperature-0.00000001922*flameTemp+0.00002092*rh
  a3 <-0.0000008381-0.000000003823*temperature+0.00000000010511*flameTemp-0.0000001166*rh
  a4 <--0.000000001685+0.000000000007637*temperature-0.0000000000002085*flameTemp+0.0000000002350*rh
  return(a0+a1*D+a2*D^2+a3*D^3+a4*D^4)
}
