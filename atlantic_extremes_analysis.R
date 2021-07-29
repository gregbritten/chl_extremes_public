library(maps)
library(ncdf4)
library(fields)
library(emdbook)
library(lubridate)
library(extRemes)
library(climextRemes)
library(viridis)

rm(list=ls())

##--LOAD PRE-PROCESSED DATA--###################
##uncomment the following for analysis of MODIS chlorophyll
load('~/dropbox/working/chlorophyll/data/chl.rdata')
load('~/dropbox/working/chlorophyll/data/sst.rdata')
load('~/dropbox/working/chlorophyll/data/time.rdata')
load('~/dropbox/working/chlorophyll/data/date.rdata')

##uncomment the following for analysis of OC-CCI chlorophyll
load('~/dropbox/working/chlorophyll/data/chl_oc.rdata')
load('~/dropbox/working/chlorophyll/data/time_oc.rdata')
load('~/dropbox/working/chlorophyll/data/date_oc.rdata')
chl <- chl_oc[,,224:1070] 

lats <- seq(-90,90,length.out=720)
lons <- seq(-180,180,length.out=1440)

latsi <- which(lats > 35 & lats < 80)
lonsi <- which(lons > -75 & lons < 35)
nlat  <- length(latsi)
nlon  <- length(lonsi)

##--READ IN BATHYMETRIC DATA--#################
## DATA DOWNLOADED FROM: https://www.gebco.net/data_and_products/gridded_bathymetry_data/
file  <- '~/dropbox/working/chlorophyll/data/GEBCO_BATHY_2002-01-01_rgb_1440x720.csv' 
bathy <- as.matrix(read.csv(file,header=FALSE))
bathy[bathy==99999] <- NA
bathy     <- -t(bathy)[,720:1]
bathy_atl <- bathy[lonsi,latsi]

##--PROCES ATLANTIC DATA--####################################
chl_atl=sst_atl <- array(NA,dim=c(nlon,nlat,dim(chl)[3]))
for(i in 1:dim(chl)[3]){
print(i)
    chl_atl[,,i] <- chl[lonsi,latsi,i]  
    sst_atl[,,i] <- sst[lonsi,latsi,i]
}

##--FUNCTION TO EXTRACT BLOCK MAXIMA--###################################
block_maxima <- function(x,date){
  years <- year(date)
  year_ind <- unique(years)
  maxima <- numeric(length(year_ind))
  for(i in 1:length(year_ind)){
    maxima[i] <- max(x[years==year_ind[i]],na.rm=TRUE)
  }
  return(maxima)
}

##--COMPUTE MEAN, MEDIANS, SDS--########################################
# chl_avg     <- apply(chl,c(1,2),function(x) mean(x, na.rm=TRUE))
# chl_atl_avg <- apply(chl_atl,c(1,2),function(x) mean(x, na.rm=TRUE))
# chl_atl_sd  <- apply(chl_atl,c(1,2),function(x) sd(x, na.rm=TRUE))
# chl_atl_med <- apply(chl_atl,c(1,2),function(x) median(x, na.rm=TRUE))

###################################################################################
##--FIT GEVDS--####################################################################
###################################################################################
##-Warnings are thrown due to bad first-try initial conditions
location=shape=scale=location_se=shape_se=scale_se=record=record_se=cormap <- matrix(NA,nrow=nlon,ncol=nlat)
location_trnd=scale_trnd=shape_trnd           <- matrix(NA,nrow=nlon,ncol=nlat)
location_trnd_se=scale_trnd_se=shape_trnd_se  <- matrix(NA,nrow=nlon,ncol=nlat)
fitBIC=fitlocBIC=fitsclBIC=fitshpBIC          <- matrix(NA,nrow=nlon,ncol=nlat)

periods              <- seq(2,100)
return               <- array(NA,dim=c(nlon,nlat,length(periods)))
qq                   <- array(NA,dim=c(nlon,nlat,2,19))
years                <- 1:19

for(i in 1:nlon){ print(i)
  for(j in 1:nlat){
    chl_tmp <- chl_atl[i,j,]
    
    if(sum(!is.na(chl_tmp))>200){
      chl_tmp[is.na(chl_tmp)] <- 0
      
      maxima        <- block_maxima(chl_tmp,date=date)
      fit           <- fevd(maxima)
      fitloc        <- fevd(maxima,location.fun=~years)
      fitscl        <- fevd(maxima,scale.fun=   ~years)
      fitshp        <- fevd(maxima,shape.fun=   ~years)
      
      location[i,j] <- fit$results$par[1] 
      shape[i,j]    <- fit$results$par[2]
      scale[i,j]    <- fit$results$par[3]
      
      location_trnd[i,j]    <- fitloc$results$par[2]
      loc_trnd_try          <- try(sqrt(diag(parcov.fevd(fitloc)))[2],silent=TRUE)
      location_trnd_se[i,j] <- ifelse(class(loc_trnd_try)=='try-error',NA,loc_trnd_try)

      scale_trnd[i,j]       <- fitscl$results$par[3]
      scl_trnd_try          <- try(sqrt(diag(parcov.fevd(fitscl)))[3],silent=TRUE)
      scale_trnd_se[i,j]    <- ifelse(class(scl_trnd_try)=='try-error',NA,scl_trnd_try)
      
      shape_trnd[i,j]       <- fitshp$results$par[4]
      shp_trnd_try          <- try(sqrt(diag(parcov.fevd(fitshp)))[4],silent=TRUE)
      shape_trnd_se[i,j]    <- ifelse(class(shp_trnd_try)=='try-error',NA,shp_trnd_try)   
      
      ses <- sqrt(diag(parcov.fevd(fit))) #extract standard errors
      location_se[i,j] <- ses[1]
      shape_se[i,j]    <- ses[2]
      scale_se[i,j]    <- ses[3]
      #return[i,j,] <- as.numeric(return.level(fit,return.period=periods))
      
      fitBIC[i,j]    <- 2*fit$results$value    + log(19)*3
      fitlocBIC[i,j] <- 2*fitloc$results$value + log(19)*4
      fitsclBIC[i,j] <- 2*fitscl$results$value + log(19)*4
      fitshpBIC[i,j] <- 2*fitshp$results$value + log(19)*4
      
      rectmp   <- try(exp(calc_logReturnPeriod_fevd(fit,returnValue=max(maxima))$logReturnPeriod),silent=TRUE)
      recsetmp <- try(exp(calc_logReturnPeriod_fevd(fit,returnValue=max(maxima))$se_logReturnPeriod),silent=TRUE)
      
      record[i,j]    <- ifelse(class(rectmp)=='try-error',NA,rectmp)
      record_se[i,j] <- ifelse(class(recsetmp)=='try-error',NA,recsetmp)

      ##--STORE QUANTILES--##################
      vv <- plot(fit,type='qq')
      
      qq[i,j,1,]  <- vv[,1]  #observed quantiles
      qq[i,j,2,]  <- vv[,2]  #fitted quantiles
      cormap[i,j] <- cor(vv[,1],vv[,2])
    }
  }
}

#############################################################################
##--TREND ANALYSES--#########################################################
#############################################################################
chl_trend=sst_trend <- matrix(NA,nrow=nlon,ncol=nlat)
for(i in 1:nlon){ 
  print(i)
  for(j in 1:nlat){
    chl_tmp <- chl_atl[i,j,]
    sst_tmp <- sst_atl[i,j,]
    
    if(sum(!is.na(chl_tmp))>200 & sum(!is.na(sst_tmp))>200){
      chl_lm  <- lm(chl_tmp ~ time)    
      sst_lm  <- lm(sst_tmp ~ time)    
      chl_trend[i,j] <- summary(chl_lm)$coefficients[2,1]    
      sst_trend[i,j] <- summary(sst_lm)$coefficients[2,1] 
    }
  }
}

#############################################################################
##--SAVE RESULTS--###########################################################
#############################################################################
DAT_oc <- list(location=location,           scale=scale,           shape=shape, #uncomment to organize results from OC-CCI data 
#DAT <- list(location=location,           scale=scale,           shape=shape,   #uncomment to organize results from MODIS data
            location_se=location_se,     scale_se=scale_se,     shape_se=shape_se,
            location_trnd=location_trnd, scale_trnd=scale_trnd, shape_trnd=shape_trnd,
            fitlocBIC=fitlocBIC,         fitsclBIC=fitsclBIC,   fitshpBIC=fitshpBIC,   fitBIC=fitBIC,
            record=record, 
            record_se=record_se,
            cormap=cormap,
            chl_trend=chl_trend, sst_trend=sst_trend,
            qq=qq)

#save(file='~/dropbox/working/chlorophyll/data/ATL_DAT.rdata', DAT) #uncomment to save MODIS results to disk
save(file='~/dropbox/working/chlorophyll/data/ATL_DAT_oc.rdata', DAT_oc) #uncomment to save OC-CCI results to disk

