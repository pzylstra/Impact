
library(frame)
library(dplyr)
library(extraDistr)
source("Fires_functions.R")
source("Heating_functions.R")
source("ChapmanRichards.R")
a <- 1  #Record to be modelled
ageStep <- 1 #Years between age classes
firstAge <- 52 #Age of the youngest age class
steps <-18 #Number of age class steps to take
tAge <- 75 #Age of the trees at the starting point
suspNS <- ""
SLOPE <- 0

site<-read.csv("C:/Users/soroush/Documents/Phil Impact package/site_mon.csv")
Structure<-read.csv("C:/Users/soroush/Documents/Phil Impact package/structure_mon.csv")
Flora<-read.csv("C:/Users/soroush/Documents/Phil Impact package/flora_mon.csv")
DefaultSpeciesParams <- read.csv("C:/Users/soroush/Documents/Phil Impact package/Traits_mon.csv")
growth <- read.csv("C:/Users/soroush/Documents/Phil Impact package/growth_mon.csv")
cover <- read.csv("C:/Users/soroush/Documents/Phil Impact package/cover_mon.csv")
weather <- read.csv("C:/Users/soroush/Documents/Phil Impact package/Cooma_FebWF.csv")
base.params <- paramBuilder(site, Structure, Flora,  DefaultSpeciesParams = DefaultSpeciesParams[,-1], 1)
base.params$value[base.params$param == "leafForm"] <- "flat"
# Set slope
base.params <- base.params %>%
  ffm_set_site_param("slope", SLOPE, "deg")
dyn <- fireDynamics(base.params, weather, growth, cover, Flora, jitters = 10, ageStep = ageStep, 
                    firstAge = firstAge, 
                    steps = steps, tAge = tAge, l = 0.1, DefaultSpeciesParams = DefaultSpeciesParams, Ms = 0.01, Pm = 1,
                    Mr = 1.001, Hs = 1.25, Hr = 1.42, a = a, suspNS = suspNS, density = 300)
