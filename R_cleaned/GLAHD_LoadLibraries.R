# do this once
# devtools::install_bitbucket("remkoduursma/plotby")
library(plotBy)

r <- require(HIEv)
if(!r)stop("Install the HIEv R package from bitbucket.org/remkoduursma/hiev")





#---------------------------------------------------------------------
#- function to load a package, and install it if necessary
Library <- function(pkg, ...){
  
  PACK <- .packages(all.available=TRUE)
  pkgc <- deparse(substitute(pkg))
  
  if(pkgc %in% PACK){
    library(pkgc, character.only=TRUE)
  } else {
    install.packages(pkgc, ...)
    library(pkgc, character.only=TRUE)
  }
  
}
#---------------------------------------------------------------------





#---------------------------------------------------------------------
#- load all the libraries (and install them if needed)
Library(doBy)
Library(magicaxis)
Library(RColorBrewer)
Library(propagate)
Library(gplots)
Library(scales)
Library(readxl)
Library(maps)
Library(mapdata)
Library(rgeos)
Library(sp)
Library(raster)
Library(nlstools)
Library(rgl)
Library(nlme)
Library(oz)
Library(RNCEP)
Library(mapplots)
Library(plotrix)
Library(plotBy)
Library(dismo)
Library(rgdal)
Library(stringr)
Library(plantecophys)
Library(dplyr)
Library(car)
Library(effects)
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#- load the custom functions that do most of the heavy lifting

source("R_cleaned/generic_functions.R")
source("R_cleaned/GLAHDAnalysisFunctions.R")
source("R_cleaned/gamplotfunctions.R")
#---------------------------------------------------------------------
