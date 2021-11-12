library(maps)
library(ncdf4)
library(fields)
library(emdbook)
library(lubridate)
library(extRemes)
library(climextRemes)
library(viridis)
rm(list=ls())

source('~/dropbox/working/chlorophyll/github_extremes_public/functions.R')

##--LOAD PRE-PROCESSED DATA--###################
##uncomment the following for analysis of MODIS chlorophyll
load('~/dropbox/working/chlorophyll/data/chl.rdata')
load('~/dropbox/working/chlorophyll/data/sst.rdata')
load('~/dropbox/working/chlorophyll/data/time.rdata')
load('~/dropbox/working/chlorophyll/data/date.rdata')

lats <- seq(-90,90,length.out=720)
lons <- seq(-180,180,length.out=1440)

latsi <- which(lats > 35 & lats < 80)
lonsi <- which(lons > -75 & lons < 35)
nlat  <- length(latsi)
nlon  <- length(lonsi)

##--PROCES ATLANTIC DATA--####################################
chl_atl=sst_atl <- array(NA,dim=c(nlon,nlat,dim(chl)[3]))
for(i in 1:dim(sst)[3]){
  print(i)
  chl_atl[,,i] <- chl[lonsi,latsi,i]  
  sst_atl[,,i] <- sst[lonsi,latsi,i]
}


###################################################################################
##--FIND MEANS AND STANDARD DEVIATIONS--###########################################
###################################################################################
means=sds=medians  <- matrix(NA,nrow=nlon,ncol=nlat)

for(i in 1:nlon){ print(i)
  for(j in 1:nlat){
    chl_tmp <- chl_atl[i,j,]
    
    if(sum(!is.na(chl_tmp))>200){
      chl_tmp[is.na(chl_tmp)] <- 0
      
      #compute maxima and month
      proc          <- block_maxima(chl_tmp,date=date)
      maxima        <- proc$maxima
      
      means[i,j]   <- mean(maxima)
      sds[i,j]     <- sd(maxima)
      medians[i,j] <- median(maxima)
      
    }
  }
}

med_sd <- stats.bin(x=c(medians),y=c(sds),breaks=seq(0,10,length.out=100))
  medsdx  <- med_sd$centers
  medsdmu <- med_sd$stats[2,]
  medsdsd <- med_sd$stats[3,]
  medsdN  <- med_sd$stats[1,]
mean_sd <- stats.bin(x=c(means),y=c(sds),breaks=seq(0,10,length.out=100))
  meansdx  <- mean_sd$centers
  meansdmu <- mean_sd$stats[2,]
  meansdsd <- mean_sd$stats[3,]
  meansdN  <- mean_sd$stats[1,]

  pdf('~/dropbox/working/chlorophyll/plots/gevd_parameter_maps_oc_revision.pdf',height=4.5,width=14)
  par(mfrow=c(1,3),mar=c(1,1,2,4),oma=c(4,4,2,2))
  
  # image.plot2(x=lons[lonsi],y=lats[latsi],z=bathy_atl/1000,zlim=c(0,5),col=cols)
  #   axis(side=1,labels=NA); axis(side=2)
  #   #contour(x=lons[lonsi],y=lats[latsi],bathy_atl,add=TRUE,levels=c(50,500,1000,1500,2000,2500,3000))
  #   map(add=TRUE,col='grey',fill=TRUE)
  #   box(lwd=3)
  #   mtext('a) Bathymetry',adj=0,line=0.5)
  #   mtext('[km]',adj=1.1,line=0.25)
  
pdf('~/dropbox/working/chlorophyll/plots/extremes_mean_median_sd_maps.pdf',height=4.5,width=14)
par(mfrow=c(1,3),oma=c(2,2,2,2),mar=c(2,2,2,2))
  image.plot2(x=lons[lonsi],y=lats[latsi],z=means,zlim=c(0,10),col=turbo(11))
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2)
  box(lwd=3)
  mtext(expression('a) Mean Extreme'),adj=0,line=0.5)

  image.plot2(x=lons[lonsi],y=lats[latsi],z=medians,zlim=c(0,10),col=turbo(11))
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2)
  box(lwd=3)
  mtext(expression('b) Median Extreme'),adj=0,line=0.5)
  
  image.plot2(x=lons[lonsi],y=lats[latsi],z=sds,zlim=c(0,10),col=turbo(11))
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2)
  box(lwd=3)
  mtext(expression('c) Standard Deviation of Extremes'),adj=0,line=0.5)
  
  mtext('Latitude',outer=TRUE,side=2,line=0.5,cex=1.2)
  mtext('Longitude',outer=TRUE,side=1,line=0.5,cex=1.2)
dev.off()
  
  

pdf('~/dropbox/working/chlorophyll/plots/extremes_mean_median_sd_scatter.pdf',height=4.5,width=10)
par(mfrow=c(1,2),oma=c(2,2,2,2),mar=c(2,2,2,2))
plot(meansdx,meansdmu,xlim=c(0,10),ylim=c(0,14),pch=19)
  segments(x0=meansdx,x1=meansdx,y0=meansdmu - meansdsd, y1=meansdmu + meansdsd)
  mtext(side=1,line=2.5,'Mean Extreme')
  mtext(side=2,line=2.5,'Standard Deviation of Extremes')
  mtext(adj=0,'a)')
plot(medsdx, medsdmu,pch=19,xlim=c(0,10),ylim=c(0,14))
  segments(x0=medsdx, x1=medsdx,y0=medsdmu - medsdsd, y1=medsdmu + medsdsd)
  mtext(side=1,line=2.5,'Median Extreme')
  mtext(adj=0,'b)')
dev.off()  

