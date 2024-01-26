rm(list = ls())
library(raster)
library(tidyverse)
library(GGally)
library(mgcv)
library(sp)

setwd("D:/Study new/PhD/R code/Spatial_Merror")

##-----------------
## Two predictors
##-----------------
rmin = raster("us_tmin.txt")/100 # To celsius. Might be worth dividing by 1000 to examine changes in tens of a degree
rmax = raster("us_tmax.txt")/100 # To celsius. Might be worth dividing by 1000 to examine changes in tens of a degree
r <- raster::stack(rmin, rmax)
par(mfrow = c(2,1), las = 1)
plot(rmin, main = "Min Temp.")
plot(rmax, main = "Max Temp.")


coords <- read.csv("2010_coord.csv") %>% 
   dplyr::select(Longitude, Latitude)
coordinates(coords) ~ Longitude + Latitude

rasValue <- raster::extract(r, coords) 
combinePointValue <- cbind(coords, rasValue)
combinePointValue %>% 
   head

rm(rmin, rmax, r, rasValue)


##-----------------
## Response
##-----------------
resp <- read.csv("Wren_PA_2010.csv")
combinePointValue$resp <- resp$x


##-----------------
## Explore data
##-----------------
combinePointValue %>% 
   head

ggpairs(combinePointValue)

