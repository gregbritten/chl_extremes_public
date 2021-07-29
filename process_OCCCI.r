
library(ncdf4)
library(lubridate)

source('~/dropbox/CODE/functions/resize_bilinear().R')

#FOLDER CONTAINING DATA DOWNLOADED FROM: https://www.oceancolour.org/
setwd('/Volumes/external/DATA/OC_CCI/v5/8day/chlor_a/chlor_a/')

yrs    <- 1997:2020
Nfiles <- length(yrs)*365/8

chl_oc          <- array(NA,dim=c(360*4,180*4,Nfiles))
date_oc=time_oc <- numeric(Nfiles)

k <- 1
for(i in 1:length(yrs)){
  print(paste0('i = ',i))
  yr <- yrs[i]
  
  setwd(paste('/Volumes/external/DATA/OC_CCI/v5/8day/chlor_a/chlor_a/',yr,'/',sep=''))  
  files <- list.files()
  
  for(j in 1:length(files)){
    print(paste0('k = ',k))
    
    nc          <- nc_open(files[j])
    chl_tmp     <- ncvar_get(nc, 'chlor_a')
    mat         <- resize_bilinear(chl_tmp,xin=8640,yin=4320,xout=360*4,yout=180*4)[,(180*4):1]
    chl_oc[,,k] <- mat
    
    date_tmp    <- ncvar_get(nc,'time')
    date_oc[k]  <- as.character(as_date(date_tmp,origin=lubridate::origin))
    time_oc[k]  <- decimal_date(as_date(date_tmp,origin=lubridate::origin))
    
    nc_close(nc)
    
    k <- k+1
  }
}

chl_oc  <- chl_oc[,,1:k]
date_oc <- date_oc[1:k]
time_oc <- time_oc[1:k]

save(file='~/dropbox/working/chlorophyll/data/chl_oc.rdata', chl_oc)
save(file='~/dropbox/working/chlorophyll/data/date_oc.rdata',date_oc)
save(file='~/dropbox/working/chlorophyll/data/time_oc.rdata',time_oc)


