library(rgdal)
library(lubridate)

source('~/dropbox/code/functions/resize_bilinear().r')

##FOLDER CONTAINING HDF FILES DOWNLOADED FROM: http://sites.science.oregonstate.edu/ocean.productivity/
chl_files  <- list.files('~/dropbox/data/ocean_productivity/modis_chl')
sst_files  <- list.files('~/dropbox/data/ocean_productivity/sst')

Nfiles <- 851

chl=sst     <- array(NA,dim=c(360*4,180*4,Nfiles))
time=date   <- numeric(Nfiles)

for(i in 1:Nfiles){
print(i)
  file <- chl_files[i],sep='')
  yr   <- substr(file,5,8)
  doy  <- as.numeric(substr(file,9,11))
  date[i] <- as.character(as.Date(doy, origin = paste(yr,"01","01",sep='-')))
  time[i] <- decimal_date(as.Date(doy, origin = paste(yr,"01","01",sep='-')))
    
  # tmp_chl    <- as.matrix(readGDAL(paste('~/dropbox/data/ocean_productivity/modis_chl/',  chl_files[i],sep='')))
  # tmp_sst    <- as.matrix(readGDAL(paste('~/dropbox/data/ocean_productivity/sst/',        sst_files[i],sep='')))
  # 
  # tmp_chl[tmp_chl == -9999] <- NA
  # tmp_sst[tmp_sst == -9999] <- NA
  # 
  # mat_chl      <- resize_bilinear(tmp_chl,xin=2160,xout=360*4,yin=1080,yout=180*4)[,720:1]
  # mat_sst      <- resize_bilinear(tmp_sst,xin=2160,xout=360*4,yin=1080,yout=180*4)[,720:1]
  # 
  # chl[,,i] <- mat_chl
  # sst[,,i] <- mat_sst
}

save(file='~/dropbox/working/chlorophyll/data/time.rdata',time)	
save(file='~/dropbox/working/chlorophyll/data/date.rdata',date)	
save(file='~/dropbox/working/chlorophyll/data/chl.rdata',chl)	
save(file='~/dropbox/working/chlorophyll/data/sst.rdata',sst)	


