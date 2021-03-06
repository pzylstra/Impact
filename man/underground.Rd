% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fauna.R
\name{underground}
\alias{underground}
\title{Fire risk for an animal sheltering underground}
\usage{
underground(
  Surf,
  Plant,
  diameter = 6,
  surface = 677,
  percentile = 0.5,
  RH = 0.2,
  moisture = 0.2,
  distance = 50,
  trail = 300,
  var = 10,
  Pressure = 1013.25,
  Altitude = 0,
  texture = "clay",
  peat = 0.1,
  grain = "fine",
  unfrozen = 1,
  depth = 0.1,
  soilTemp = 25,
  updateProgress = NULL
)
}
\arguments{
\item{Surf}{The dataframe 'runs' exported from Monte Carlos as 'Summary.csv'}

\item{Plant}{The dataframe 'IP' exported from Monte Carlos as 'IP.csv'.}

\item{diameter}{Diameter of the surface fuels burning (mm)}

\item{surface}{Temperature at the surface of the soil, under burning fuels}

\item{percentile}{defines which heating statistics are used for each second, from 0 (min) to 1 (max)}

\item{RH}{The relative humidity (0-1)}

\item{moisture}{The proportion oven-dry weight of moisture in the bark and wood}

\item{distance}{The furthest horizontal distance between the flame origin and the point (m)}

\item{trail}{Number of seconds to continue modelling after the front has passed}

\item{var}{The angle in degrees that the plume spreads above/below a central vector;defaults to 10}

\item{Pressure}{Sea level atmospheric pressure (hPa)}

\item{Altitude}{Height above sea level (m)}

\item{texture}{Soil texture. Allowable values are: "sand", "loamy sand", "sandy loam", "sandy clay loam",
"sand clay", "loam", "clay loam", "silt loam", "clay", "silty clay", "silty clay loam", "silt"}

\item{peat}{Organic proportion of the soil}

\item{grain}{Allowable values are "fine" or "coarse"}

\item{unfrozen}{Proportion of soil unfrozen, between 0 and 1}

\item{depth}{The depth at which the animal shelters beneath the soil}

\item{soilTemp}{The starting temperature under the ground (deg C)}

\item{updateProgress}{Progress bar for use in the dashboard}
}
\value{
dataframe
}
\description{
Calculates the likelihood of mortality to an animal
caused by an approaching fire front
}
\details{
Utilises the output tables from 'threat' and 'radiation', and adds to these
the Reynolds Number, heat transfer coefficients, Newton's convective energy transfer coefficient,
and the temperature of the object each second.

Reynolds Number utilises a standard formulation (e.g. Gordon, N. T., McMahon, T. A. & Finlayson, B. L.
Stream hydrology: an introduction for ecologists. (Wiley, 1992))

Convective heat transfer coefficients use the widely adopted formulations of
Williams, F. A. Urban and wildland fire phenomenology. Prog. Energy Combust. Sci. 8, 317–354 (1982),
and Drysdale, D. An introduction to fire dynamics. (John Wiley and Sons, 1985)
utilising a Prandtl number of 0.7.

Finds animal mortality within a hollow based on the maximum tolerable temperature for a given
vapour pressure deficit, based on data from Lawrence, G. E.
Ecology of vertebrate animals in relation to chaparral fire in the Sierra Nevada foothills.
Ecology 47, 278–291 (1966)

Heat is transferred into the earth using Fourier's Law. Spread continues for a period after the
passage of the fire front, equal to the duration of the surface flame, as determined using
Burrows, N. D. Flame residence times and rates of weight loss of eucalypt forest fuel particles.
Int. J. Wildl. Fire 10, 137–143 (2001).

Default temperature of the resident flame is the average of the surface maximums in
Cawson, J. G., Nyman, P., Smith, H. G., Lane, P. N. J. & Sheridan, G. J.
How soil temperatures during prescribed burning affect soil water repellency,
infiltration and erosion. Geoderma 278, 12–22 (2016).

Heating area is set to 1m2, flat, with a characteristic length of 1m
}
