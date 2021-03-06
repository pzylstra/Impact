% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Growth.R
\name{coverChange}
\alias{coverChange}
\title{Models the weighted mean of plant separation from time since fire for a named stratum}
\usage{
coverChange(st, a, cover, Flora, age)
}
\arguments{
\item{st}{The name of the stratum}

\item{a}{The number of the record to be modelled}

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

\item{age}{The number of years since last fire}
}
\value{
dataframe
}
\description{
Checks for alternate growth models
}
