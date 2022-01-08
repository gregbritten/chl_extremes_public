rm(list=ls())

library(lubridate)
library(extRemes)

source('r/functions.R')
load_data(modis=TRUE, oc=TRUE, bathy=TRUE)

#source('process_bathy.r')
#source('load_data.r')

###################################################################################
##--FIT GEVDS--####################################################################
###################################################################################
##-Warnings are thrown due to bad first-try initial conditions
DAT    <- fit(chl=chl,   nlon=nlon,nlat=nlat,date=date,   time=time,   oc=FALSE)
DAT_oc <- fit(chl=chl_oc,nlon=nlon,nlat=nlat,date=date_oc,time=time_oc,oc=TRUE)

save(file='data/DAT.rdata',    DAT) 
save(file='data/DAT_oc.rdata', DAT_oc) 

