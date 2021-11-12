library(sf)
library(R.matlab)
library(ggplot2)
source('~/dropbox/working/chlorophyll/github_extremes_public/functions.r')

long <- readMat('~/dropbox/data/longhurst/Longhurst_180.mat')$Longhurst
image.plot(long)

longdeg <- resize_bilinear(long,xin=360,xout=360*4, yin=180, yout=180*4)

image.plot(round(longdeg))

lons <- seq(-180,180,length.out=360*4)
lats <- seq(-90,90,length.out=180*4)

latsi <- which(lats > 35 & lats < 80)
lonsi <- which(lons > -75 & lons < 35)

image.plot(x=lons[lonsi], y=lats[latsi], round(longdeg)[lonsi,latsi])

#########################################################
## LONGHURST SHAPEFILES--################################
#########################################################
##--READ SHAPE FILE--#############################
longshp <- st_read("~/dropbox/data/longhurst/longhurst_v4_2010/Longhurst_world_v4_2010.shp")
longshp$NUM <- as.numeric(row.names(longshp))

##--RASTERIZE SHAPE FILE--#######################
r0       <- raster(nrow=180*4,ncol=360*4)
longrast <- fasterize(longshp,r0,field="NUM")
longmat  <- t(as.matrix(longrast))[,(180*4):1]

##--PLOT--############################
#image.plot(lons[lonsi],lats[latsi],longmat[lonsi,latsi],col=turbo(20))

pdf('~/dropbox/working/chlorophyll/plots/gevd_parameter_maps_longhurst.pdf',height=4.5,width=14)
par(mfrow=c(1,3),mar=c(1,1,2,4),oma=c(4,4,2,2))
image.plot2(x=lons[lonsi],y=lats[latsi],DAT$location,zlim=c(0,4),col=turbo(20))
  plot(st_geometry(longshp),add=TRUE,border='white',lwd=2)
  map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2)
  box(lwd=3)
  #mtext('Doubling Sample Size',line=0.5)
  mtext(expression('a) Location'~mu),adj=0,line=0.5)
  mtext(expression('[mgChl/m'^3*']'),adj=1.1,line=0.25)
  
image.plot2(x=lons[lonsi],y=lats[latsi],scale,zlim=c(0,2),col=turbo(20))
  plot(st_geometry(longshp),add=TRUE,border='white',lwd=2)
  map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2,labels=NA)
  box(lwd=3)
  #mtext('Doubling Sample Size',line=0.5)
  mtext(expression('b) Scale'~sigma),adj=0,line=0.5)
  mtext(expression('[mgChl/m'^3*']'),adj=1.1,line=0.25)
  
image.plot2(x=lons[lonsi],y=lats[latsi],shape,zlim=c(0,2.5),col=turbo(20))
  plot(st_geometry(longshp),add=TRUE,border='white',lwd=2)
  map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2,labels=NA)
  box(lwd=3)
  #mtext('Doubling Sample Size',line=0.5)
  mtext(expression('c) Shape'~xi),adj=0,line=0.5)
  mtext(expression('[unitless]'),adj=1.1,line=0.25)
  
  mtext('Latitude',outer=TRUE,side=2,line=1.5,cex=1.2)
  mtext('Longitude',outer=TRUE,side=1,line=1.5,cex=1.2)
  
dev.off()


  
  
  
