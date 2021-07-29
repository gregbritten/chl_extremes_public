library(raster)

# load data saved by atlantic_extremes_analysis.r
load('~/dropbox/working/chlorophyll/data/ATL_DAT.rdata')

####################################################################
##--AREA--##########################################################
dBIC_loc             <- DAT$fitBIC - DAT$fitlocBIC
dBIC_loc[dBIC_loc<0] <- 0
dBIC_scl             <- DAT$fitBIC - DAT$fitsclBIC
dBIC_scl[dBIC_scl<0] <- 0
dBIC_shp             <- DAT$fitBIC - DAT$fitshpBIC
dBIC_shp[dBIC_shp<0] <- 0

r <- raster(nrow=180*4, ncol=360*4)
a <- t(as.matrix(area(r)))

lats <- seq(-90,90,length.out=720)
lons <- seq(-180,180,length.out=1440)

latsi <- which(lats > 35 & lats < 80)
lonsi <- which(lons > -75 & lons < 35)

a_ATL <- a[lonsi,latsi]
a_ATL[is.na(chl_atl[,,1])] <- NA

a_tot <- sum(a_ATL,na.rm=TRUE)

##--PERCENT AREA SUPPORTED BY BIC--#################################
1 - sum(a_ATL[dBIC_loc==0],na.rm=TRUE)/a_tot
1 - sum(a_ATL[dBIC_scl==0],na.rm=TRUE)/a_tot
1 - sum(a_ATL[dBIC_shp==0],na.rm=TRUE)/a_tot

###################################################################
## COEFFICIENT OF VARIATION #######################################
###################################################################
#area-unweighted mean
mean(DAT$location_se/DAT$location, na.rm=TRUE)
mean(DAT$scale_se/DAT$scale, na.rm=TRUE)
mean(DAT$shape_se/DAT$shape, na.rm=TRUE)

#area-weighted
sum(a_ATL*(DAT$location_se/DAT$location),na.rm=TRUE)/a_tot
sum(a_ATL*(DAT$scale_se/DAT$scale),na.rm=TRUE)/a_tot
sum(a_ATL*(DAT$shape_se/DAT$shape),na.rm=TRUE)/a_tot



