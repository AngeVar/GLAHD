#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- Central analysis script for the GLAHD manuscript.
#  The idea is to keep this script nice and tidy, but reproducibly do all the
#  analysis and make all of the figures for the manuscript.
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- load the packages and custom functions that do all the work
source("R_cleaned/1. GLAHD_LoadLibraries.R")
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- download all of the required data from HIEv.
setToken(tokenfile="HIEv_token.txt") #- set HIEv token. See ?setToken
source("R_cleaned/2. downloadData.R")
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Make datasets.
source("R_cleaned/3. Create_datasets.R")
#-------------------------------------------------------------------------------------
       
#-------------------------------------------------------------------------------------
#- Make figure 1.Figure of experimental design.
source("R/"
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Make figure 2.
source("R_cleaned/Figure 2.R") #biomass over time
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Make figure 3.
source("R_cleaned/Figure 3.R") #biomass and RGR response ratios
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Make figure 4.
source("R_cleaned/Figure 4.R")
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Make figure 5.
source("R/"
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Make figure 6.
source("R_cleaned/Figure 6.R")
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Make Table 1.
source("R/make_table")
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Make Table 2. 
source("R/make_table")
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Make Figure S1. Location of provenances
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Make Figure S2. Growth model comparisons
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Make Table S1. Taxa-specific allometries.
source("R/make_table")
#-------------------------------------------------------------------------------------
