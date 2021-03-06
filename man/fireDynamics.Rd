% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flamDynamics.R
\name{fireDynamics}
\alias{fireDynamics}
\title{Models flammability dynamics from a weather set}
\usage{
fireDynamics(
  base.params,
  weather,
  growth,
  cover,
  Flora,
  jitters = 5,
  ageStep = 5,
  firstAge = 1,
  steps = 10,
  tAge = 50,
  l = 0.1,
  default.species.params = default.species.params,
  Ms = 0.01,
  Pm = 1,
  Mr = 1.001,
  Hs = 0.2,
  Hr = 1.41,
  a,
  suspNS = "suspNS",
  density = 300,
  test = 80,
  hKill = 90,
  db.path = "out_mc.db",
  updateProgress = NULL
)
}
\arguments{
\item{base.params}{Parameter input table}

\item{weather}{A dataframe with the four fields:
tm - Sequential numbering of the records
T - Air temperature (deg C)
W - Wind velocity (km/h)
DFMC - Dead fuel moisture content (proportion ODW)}

\item{growth}{A dataframe with the six fields:
species - Name of the species consistent with other tables
max - Maximum plant height (m)}

\item{cover}{A dataframe with the fields:
species - Name of the species consistent with other tables
constant - Mean plant separation (m) that does not change with age
exp_a    - The first constant in an exponential function describing plant separation with tsf
exp_b    - The second constant in an exponential function describing plant separation with tsf}

\item{Flora}{A dataframe with the fields:
record - a unique, consecutively numbered identifier per site
species - the name of the species, which will call trait data from 'default.species.params'
moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
stratum - numeric value from 1 to 4, counting from lowest stratum
comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
base, he, ht, top & w - canopy dimensions for that species (m). he and ht are optional
clump - mean ratio of clump diameter to crown diameter
openness - ratio of gap to clump size}

\item{jitters}{Number of repetitions for each row in the weather table}

\item{ageStep}{The number of years between age classes to be examined}

\item{firstAge}{The youngest age class to be examined}

\item{steps}{The number of age classes to be examined after firstAge}

\item{tAge}{The age of trees at the commencement of the run}

\item{l}{Variation around input leaf dimensions}

\item{default.species.params}{Leaf traits database}

\item{Ms}{Standard deviation of LFMC}

\item{Pm}{Multiplier of mean LFMC}

\item{Mr}{Truncates LFMC variability by +/- Pm * LFMC}

\item{Hs}{Standard deviation of plant height variations}

\item{Hr}{Truncates plant height variability by +/- Hr * height}

\item{a}{The number of the record to be modelled}

\item{suspNS}{Name of the fuel in the suspended layer}

\item{density}{density of the wood in the suspended layer (kg/m3)}

\item{test}{An isotherm for which the height will be reported}

\item{hKill}{The percentage of tree canopy above which a risk boundary is crossed (0 - 100)}

\item{updateProgress}{Progress bar for use in the dashboard}
}
\value{
dataframe
}
\description{
Grows the plant community,
modelling fire behaviour against a weather dataset for every age, then summarising all
behaviour into risk and impact statistics. Allows trees an older age than understorey to
represent lower severity fire histories
}
