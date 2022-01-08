library(rgdal)
library(lubridate)

source('r/functions.r')

##FOLDER CONTAINING HDF FILES DOWNLOADED FROM: http://sites.science.oregonstate.edu/ocean.productivity/
chl_files  <- list.files('~/dropbox/data/ocean_productivity/modis_chl',full.names=TRUE)
sst_files  <- list.files('~/dropbox/data/ocean_productivity/sst',full.names=TRUE)

Nfiles <- 851
time=date   <- numeric(Nfiles)

lats <- seq(-90,90,length.out=720)
lons <- seq(-180,180,length.out=1440)

latsi <- which(lats > 35 & lats < 80)
lonsi <- which(lons > -75 & lons < 35)
nlat  <- length(latsi)
nlon  <- length(lonsi)

chl=sst     <- array(NA,dim=c(nlon,nlat,Nfiles))

for(i in 1:Nfiles){
print(i)
  file <- chl_files[i]
  yr   <- substr(file,69,72)
  doy  <- as.numeric(substr(file,73,75))
  date[i] <- as.character(as.Date(doy, origin = paste(yr,"01","01",sep='-')))
  time[i] <- decimal_date(as.Date(doy, origin = paste(yr,"01","01",sep='-')))
    
  tmp_chl    <- as.matrix(readGDAL(chl_files[i]))
  tmp_sst    <- as.matrix(readGDAL(sst_files[i]))

  tmp_chl[tmp_chl == -9999] <- NA
  tmp_sst[tmp_sst == -9999] <- NA

  mat_chl      <- resize_bilinear(z=tmp_chl,xin=2160,xout=360*4,yin=1080,yout=180*4)[,720:1]
  mat_sst      <- resize_bilinear(tmp_sst,xin=2160,xout=360*4,yin=1080,yout=180*4)[,720:1]

  chl[,,i] <- mat_chl[lonsi,latsi]
  sst[,,i] <- mat_sst[lonsi,latsi]
}

save(file='data/time.rdata',time)	
save(file='data/date.rdata',date)	
save(file='data/chl.rdata',chl)	
save(file='data/sst.rdata',sst)	


